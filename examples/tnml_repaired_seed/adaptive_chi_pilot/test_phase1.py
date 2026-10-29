#!/usr/bin/env python3
import hashlib
import importlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import core
import manifest
import plan_task


def level(chi, probabilities, log_joint, *, status="complete", valid=True):
    document = dict(chi=chi, status=status, valid=valid)
    if status == "complete":
        document.update(
            probabilities=probabilities,
            log_joint_syndrome=log_joint,
            logical_mask=max(range(4), key=probabilities.__getitem__),
        )
    return document


class StoppingRuleTests(unittest.TestCase):
    def setUp(self):
        self.ladder = [32, 64, 128, 240]
        self.stable = [
            level(32, [0.7, 0.1, 0.1, 0.1], -10.0),
            level(64, [0.7000004, 0.0999998, 0.0999999, 0.0999999], -9.9999996),
            level(128, [0.7000005, 0.09999975, 0.099999875, 0.099999875], -9.9999995),
        ]

    def test_accepts_only_after_two_stable_successive_differences(self):
        decision = core.evaluate_stopping(self.stable, self.ladder)
        self.assertEqual(decision["status"], "accepted")
        self.assertEqual(decision["accepted_chi"], 128)
        self.assertEqual(decision["logical_mask"], 0)
        self.assertLessEqual(decision["posterior_bound"], 1e-6)

    def test_unstable_class_escalates_then_becomes_unknown_at_maximum(self):
        unstable = [
            level(32, [0.6, 0.2, 0.1, 0.1], -9.0),
            level(64, [0.2, 0.6, 0.1, 0.1], -9.0),
            level(128, [0.6, 0.2, 0.1, 0.1], -9.0),
        ]
        decision = core.evaluate_stopping(unstable, self.ladder)
        self.assertEqual(decision, {
            "status": "continue", "reason": "not_converged", "next_chi": 240
        })
        decision = core.evaluate_stopping(
            unstable + [level(240, [0.2, 0.6, 0.1, 0.1], -9.0)],
            self.ladder,
        )
        self.assertEqual(decision["status"], "unknown")
        self.assertEqual(decision["reason"], "max_chi_without_convergence")

    def test_near_tie_does_not_pass_margin_guard(self):
        levels = [
            level(32, [0.400000004, 0.399999996, 0.1, 0.1], -8.0),
            level(64, [0.400000004, 0.399999996, 0.1, 0.1], -8.0),
            level(128, [0.400000004, 0.399999996, 0.1, 0.1], -8.0),
        ]
        self.assertEqual(
            core.evaluate_stopping(levels, self.ladder)["status"], "continue"
        )

    def test_exact_tie_is_valid_but_forces_escalation(self):
        tied = level(32, [0.4, 0.4, 0.1, 0.1], -8.0)
        self.assertEqual(
            core.evaluate_stopping([tied], self.ladder),
            {"status": "continue", "reason": "not_converged", "next_chi": 64},
        )

    def test_forced_anchor_preserves_first_acceptance_until_maximum(self):
        provisional = core.evaluate_stopping(
            self.stable, self.ladder, forced_anchor=True
        )
        self.assertEqual(provisional["status"], "continue_anchor")
        self.assertEqual(provisional["provisional_acceptance"]["accepted_chi"], 128)
        self.assertEqual(provisional["next_chi"], 240)
        final = core.evaluate_stopping(
            self.stable
            + [level(240, [0.70000055, 0.099999725, 0.0999998625, 0.0999998625], -9.99999945)],
            self.ladder,
            forced_anchor=True,
        )
        self.assertEqual(final["status"], "accepted")
        self.assertEqual(final["accepted_chi"], 128)
        self.assertTrue(final["anchor_completed"])
        self.assertEqual(final["anchor_final_chi"], 240)
        self.assertFalse(final["anchor_class_changed"])

        changed = core.evaluate_stopping(
            self.stable + [level(240, [0.1, 0.7, 0.1, 0.1], -9.5)],
            self.ladder,
            forced_anchor=True,
        )
        self.assertEqual(changed["status"], "accepted")
        self.assertEqual(changed["logical_mask"], 0)
        self.assertTrue(changed["anchor_class_changed"])
        self.assertAlmostEqual(changed["anchor_posterior_linf"], 0.6000005)

    def test_bad_or_missing_level_propagates_unknown(self):
        timed_out = self.stable[:1] + [
            dict(chi=64, status="timeout", valid=False)
        ]
        decision = core.evaluate_stopping(timed_out, self.ladder)
        self.assertEqual(decision["status"], "unknown")
        self.assertEqual(decision["reason"], "timeout")
        with self.assertRaisesRegex(ValueError, "prefix"):
            core.evaluate_stopping([self.stable[1]], self.ladder)

    def test_malformed_complete_result_propagates_unknown(self):
        malformed = [
            self.stable[0],
            level(64, [float("nan"), 0.1, 0.1, 0.1], -10.0),
        ]
        decision = core.evaluate_stopping(malformed, self.ladder)
        self.assertEqual(decision, {"status": "unknown", "reason": "invalid"})

    def test_summary_keeps_every_nonaccepted_outcome_unknown(self):
        summary = core.summarize_outcomes([
            dict(status="accepted", failure=True),
            dict(status="accepted", failure=False),
            dict(status="unknown", reason="timeout"),
            dict(status="continue", reason="not_converged"),
        ])
        self.assertEqual(summary["requested"], 4)
        self.assertEqual(summary["accepted"], 2)
        self.assertEqual(summary["known_failures"], 1)
        self.assertEqual(summary["unknown"], 2)
        self.assertEqual(summary["failure_interval"], [0.25, 0.75])
        self.assertEqual(summary["unknown_reasons"], {
            "not_converged": 1, "timeout": 1
        })


class AnchorRankingTests(unittest.TestCase):
    def test_golden_digest_and_probability_stratum_validation(self):
        record_id = f"{2:064x}"
        ranked = core.rank_record_ids(3, 0.006, [record_id], count=1)
        self.assertEqual(
            ranked[0]["digest_hex"],
            "0241abd8c8cdc1df24c998c05494c8ce7d05189b0e9d13803b7c22dc013ea854",
        )
        with self.assertRaisesRegex(ValueError, "six preregistered"):
            core.rank_record_ids(3, 0.0061, [record_id], count=1)

    def test_sha256_big_endian_ranking_is_exact_and_input_order_independent(self):
        record_ids = [f"{value:064x}" for value in (7, 2, 11, 5, 3)]
        ranked = core.rank_record_ids(5, 0.035, record_ids, count=3)
        expected = []
        for record_id in record_ids:
            message = (
                "tnml-adaptive-chi-v1|L=5|p=0.0350|record_id=" + record_id
            ).encode("utf-8")
            digest = hashlib.sha256(message).digest()
            expected.append((int.from_bytes(digest, "big"), record_id, digest.hex()))
        expected.sort()
        self.assertEqual(
            [item["record_id"] for item in ranked],
            [item[1] for item in expected[:3]],
        )
        self.assertEqual(ranked[0]["digest_hex"], expected[0][2])
        self.assertEqual(ranked[0]["digest_integer"], str(expected[0][0]))
        self.assertEqual([item["rank"] for item in ranked], [1, 2, 3])
        self.assertEqual(
            ranked,
            core.rank_record_ids(5, 0.035, reversed(record_ids), count=3),
        )

    def test_ranking_rejects_duplicates_and_bad_record_ids(self):
        good = "ab" * 32
        with self.assertRaisesRegex(ValueError, "duplicate"):
            core.rank_record_ids(5, 0.006, [good, good], count=1)
        with self.assertRaisesRegex(ValueError, "64 lowercase"):
            core.rank_record_ids(5, 0.006, ["ABC"], count=1)
        with self.assertRaisesRegex(ValueError, "64 lowercase"):
            core.rank_record_ids(5, 0.006, [["not", "hashable"]], count=1)


class ManifestAndPlannerTests(unittest.TestCase):
    def test_manifest_enumerates_six_strata_and_is_immutable(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            document = manifest.create_manifest(root, seed_base=202610600000)
            self.assertEqual(document["p_values"], list(manifest.P_VALUES))
            l3 = manifest.tasks_for_stage(document, "l3_exact")
            l5 = manifest.tasks_for_stage(document, "l5_feasibility")
            self.assertEqual(len(l3), 30)
            self.assertEqual(len(l5), 60)
            for p_key in manifest.P_VALUES:
                self.assertEqual(sum(task["p_key"] == p_key for task in l3), 5)
                self.assertEqual(sum(task["p_key"] == p_key for task in l5), 10)
            self.assertTrue(all(task["shots"] == 20 and task["exact"] for task in l3))
            self.assertTrue(all(task["shots"] == 1 and not task["exact"] for task in l5))
            l3_history_seeds = [seed for task in l3 for seed in task["history_seeds"]]
            l5_history_seeds = [seed for task in l5 for seed in task["history_seeds"]]
            self.assertEqual(len(l3_history_seeds), 600)
            self.assertEqual(len(set(l3_history_seeds)), 600)
            self.assertEqual(len(l5_history_seeds), 60)
            self.assertEqual(len(set(l5_history_seeds)), 60)
            self.assertTrue(set(l3_history_seeds).isdisjoint(l5_history_seeds))
            self.assertEqual([task["task_id"] for task in l3], list(range(30)))
            self.assertEqual([task["task_id"] for task in l5], list(range(60)))
            self.assertEqual(document["runtime_policy"]["requested_cpus_per_task"], 4)
            self.assertEqual(document["runtime_policy"]["blas_threads"], 4)
            self.assertEqual(document["runtime_policy"]["julia_application_threads"], 1)
            self.assertEqual(document["runtime_policy"]["max_rss_bytes"], 16 * 2**30)
            self.assertEqual(document["planned_batches"]["l5_feasibility"], [
                {"batch": 0, "first_task_id": 0, "last_task_id": 29,
                 "max_concurrent": 20},
                {"batch": 1, "first_task_id": 30, "last_task_id": 59,
                 "max_concurrent": 20},
            ])
            for first, last in ((0, 29), (30, 59)):
                batch = l5[first:last + 1]
                for p_key in manifest.P_VALUES:
                    self.assertEqual(sum(task["p_key"] == p_key for task in batch), 5)
            self.assertEqual(manifest.load_manifest(root), document)
            self.assertEqual(manifest.create_manifest(root, seed_base=202610600000), document)
            with self.assertRaisesRegex(ValueError, "different seed_base"):
                manifest.create_manifest(root, seed_base=1)
            for invalid_seed in (True, 1.0):
                with self.subTest(invalid_seed=invalid_seed):
                    with self.assertRaisesRegex(ValueError, "nonnegative integer"):
                        manifest.create_manifest(root, seed_base=invalid_seed)

    def test_manifest_identity_detects_tampering(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest.create_manifest(root, seed_base=1234)
            path = root / "manifest.json"
            document = json.loads(path.read_text())
            document["stages"]["l3_exact"]["tasks"][0]["p_key"] = "0.9999"
            path.write_text(json.dumps(document))
            with self.assertRaisesRegex(ValueError, "identity"):
                manifest.load_manifest(root)

    def test_dry_run_planner_is_read_only_and_bounds_checked(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            document = manifest.create_manifest(root, seed_base=99)
            task = manifest.task_for_id(document, "l5_feasibility", 0)
            output = root / task["output"]
            self.assertFalse(output.exists())
            plan = plan_task.plan_task(root, "l5_feasibility", 0)
            self.assertEqual(plan["status"], "dry_run")
            self.assertEqual(plan["execution"], "disabled_phase1")
            self.assertEqual(plan["task"], task)
            self.assertEqual(plan["output_directory"], str(output.resolve()))
            self.assertFalse(output.exists())
            with self.assertRaisesRegex(ValueError, "outside stage"):
                plan_task.plan_task(root, "l5_feasibility", 60)

    def test_dry_run_cli_emits_json_without_creating_task_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            document = manifest.create_manifest(root, seed_base=77)
            task = manifest.task_for_id(document, "l3_exact", 3)
            completed = subprocess.run(
                [sys.executable, str(HERE / "plan_task.py"),
                 "--root", str(root), "--stage", "l3_exact", "--task-id", "3"],
                check=True, capture_output=True, text=True,
            )
            payload = json.loads(completed.stdout)
            self.assertEqual(payload["status"], "dry_run")
            self.assertFalse((root / task["output"]).exists())

    def test_modules_are_importable_as_a_package(self):
        imported_manifest = importlib.import_module(
            "examples.tnml_repaired_seed.adaptive_chi_pilot.manifest"
        )
        imported_planner = importlib.import_module(
            "examples.tnml_repaired_seed.adaptive_chi_pilot.plan_task"
        )
        self.assertEqual(imported_manifest.P_VALUES, manifest.P_VALUES)
        self.assertTrue(callable(imported_planner.plan_task))


class CheckpointTests(unittest.TestCase):
    def test_complete_checkpoint_validates_identity_and_payload(self):
        checkpoint = core.make_level_checkpoint(
            manifest_id="1" * 64,
            task_id=4,
            record_id="2" * 64,
            chi=128,
            cutoff=1e-18,
            status="complete",
            valid=True,
            probabilities=[0.7, 0.1, 0.1, 0.1],
            logical_mask=0,
            log_joint_syndrome=-12.5,
            elapsed_seconds=1.25,
            peak_rss_bytes=123456,
        )
        core.validate_level_checkpoint(
            checkpoint,
            manifest_id="1" * 64,
            task_id=4,
            record_id="2" * 64,
            chi=128,
            cutoff=1e-18,
        )
        tampered = dict(checkpoint)
        tampered["probabilities"] = [0.6, 0.2, 0.1, 0.1]
        with self.assertRaisesRegex(ValueError, "identity"):
            core.validate_level_checkpoint(
                tampered,
                manifest_id="1" * 64,
                task_id=4,
                record_id="2" * 64,
                chi=128,
                cutoff=1e-18,
            )

    def test_unknown_checkpoint_has_no_decoded_class(self):
        checkpoint = core.make_level_checkpoint(
            manifest_id="a" * 64,
            task_id=0,
            record_id="b" * 64,
            chi=64,
            cutoff=1e-18,
            status="timeout",
            valid=False,
            elapsed_seconds=300.0,
            peak_rss_bytes=1024,
        )
        validated = core.validate_level_checkpoint(
            checkpoint,
            manifest_id="a" * 64,
            task_id=0,
            record_id="b" * 64,
            chi=64,
            cutoff=1e-18,
        )
        self.assertEqual(validated["status"], "timeout")
        self.assertNotIn("logical_mask", validated)
        self.assertNotIn("probabilities", validated)

    def test_checkpoint_rejects_unrecognized_status(self):
        with self.assertRaisesRegex(ValueError, "status"):
            core.make_level_checkpoint(
                manifest_id="a" * 64,
                task_id=0,
                record_id="b" * 64,
                chi=64,
                cutoff=1e-18,
                status="tiemout",
                valid=False,
                elapsed_seconds=1.0,
                peak_rss_bytes=1024,
            )

    def test_checkpoint_accepts_a_tied_maximizer(self):
        checkpoint = core.make_level_checkpoint(
            manifest_id="a" * 64,
            task_id=0,
            record_id="b" * 64,
            chi=64,
            cutoff=1e-18,
            status="complete",
            valid=True,
            probabilities=[0.4, 0.4, 0.1, 0.1],
            logical_mask=1,
            log_joint_syndrome=-1.0,
            elapsed_seconds=1.0,
            peak_rss_bytes=1024,
        )
        self.assertEqual(checkpoint["logical_mask"], 1)

    def test_checkpoint_rejects_nonmaximizing_or_noninteger_logical_class(self):
        for probabilities, logical_mask in (
            ([0.4, 0.4, 0.1, 0.1], 2),
            ([0.7, 0.1, 0.1, 0.1], 0.0),
            ([0.1, 0.7, 0.1, 0.1], True),
        ):
            with self.subTest(probabilities=probabilities, logical_mask=logical_mask):
                with self.assertRaisesRegex(ValueError, "logical_mask"):
                    core.make_level_checkpoint(
                        manifest_id="a" * 64,
                        task_id=0,
                        record_id="b" * 64,
                        chi=64,
                        cutoff=1e-18,
                        status="complete",
                        valid=True,
                        probabilities=probabilities,
                        logical_mask=logical_mask,
                        log_joint_syndrome=-1.0,
                        elapsed_seconds=1.0,
                        peak_rss_bytes=1024,
                    )


if __name__ == "__main__":
    unittest.main()
