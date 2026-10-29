#!/usr/bin/env python3
"""Bounded tests for the isolated adaptive HPC control plane."""

from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

if __package__:
    from . import adaptive_runner, core, execution_manifest, manifest, task_runner
    from .hpc import control, report
else:
    HERE = Path(__file__).resolve().parent
    sys.path.insert(0, str(HERE))
    import adaptive_runner
    import core
    import execution_manifest
    import manifest
    import task_runner
    from hpc import control, report


def fake_source_identity():
    files = {
        path: hashlib.sha256(path.encode()).hexdigest()
        for path in sorted(
            {"src/TopoNoise.jl"} | set(execution_manifest.REQUIRED_SOURCE_FILES)
        )
    }
    return {
        "project_sha256": "1" * 64,
        "manifest_sha256": "2" * 64,
        "package_source_tree_sha256": "3" * 64,
        "files": files,
    }


def fake_runtime(blas_threads):
    return {
        "julia_version": "1.12.4",
        "python_version": "3.12.14",
        "blas_vendor": "openblas64",
        "cpu_model": "test-aarch64",
        "julia_threads": 1,
        "blas_threads": blas_threads,
    }


def fake_record(shot_index, salt):
    payload = {
        "shot_index": shot_index,
        "raw_fault_count": 12,
        "raw_fault_indices": [1] if salt % 2 else [],
        "seed_syndrome": [bool(salt % 2), False],
        "seed_recovery": [bool(salt % 2), False],
        "recovery_mask": salt % 4,
        "final_syndrome": [False, bool(salt % 3)],
        "true_mask": (salt + 1) % 4,
        "mwpm_failure": bool(salt % 5 == 0),
        "raw_seed_mask": salt % 4,
        "transfer_mask": (salt + 2) % 4,
    }
    payload["record_id"] = hashlib.sha256(
        f"julia:{shot_index}:{salt}".encode()
    ).hexdigest()
    return payload


class HPCControlTests(unittest.TestCase):
    @staticmethod
    def _l3_rows(sampling):
        rows = []
        salt = 1
        for task in manifest.tasks_for_stage(sampling, "l3_exact"):
            for shot_index, seed in enumerate(task["history_seeds"], 1):
                rows.append({"task_id": task["task_id"], "shot_index": shot_index,
                             "history_seed": seed,
                             "record": fake_record(shot_index, salt)})
                salt += 1
        return rows

    def test_sampling_plan_uses_fresh_per_history_seed_and_safe_stage_roots(self):
        sampling = manifest.build_manifest(202610600000)
        plan = control.plan_corpus_task(
            Path("/study"), sampling, stage="l3_exact", task_id=7
        )
        task = manifest.task_for_id(sampling, "l3_exact", 7)
        self.assertEqual(plan["history_seeds"], task["history_seeds"])
        self.assertEqual(plan["seed_mode"], "fresh_xoshiro_per_history_seed")
        self.assertEqual(plan["record_count"], 20)
        self.assertTrue(plan["shard_path"].endswith("task-007.json"))
        self.assertNotEqual(
            control.stage_root(Path("/study"), "l3_exact", replay=False),
            control.stage_root(Path("/study"), "l3_exact", replay=True),
        )

    def test_shards_and_collection_bind_all_locators_source_and_runtime(self):
        sampling = manifest.build_manifest(101)
        source = fake_source_identity()
        profile = control.build_environment_profile(
            source_identity=source,
            primary_runtime=fake_runtime(4),
            replay_runtime=fake_runtime(1),
            extended_runtime={"platform_machine": "aarch64", "python_executable": "/x/python"},
        )
        shards = []
        salt = 1
        for task in manifest.tasks_for_stage(sampling, "l5_feasibility"):
            rows = []
            for shot_index, seed in enumerate(task["history_seeds"], 1):
                rows.append({"history_seed": seed, "record": fake_record(shot_index, salt)})
                salt += 1
            shards.append(
                control.build_corpus_shard(
                    sampling, "l5_feasibility", task["task_id"], rows,
                    source_identity=source,
                    environment_profile_id=profile["environment_profile_id"],
                )
            )
        corpus, provenance = control.collect_corpus_documents(
            sampling, "l5_feasibility", shards,
            source_identity=source,
            environment_profile_id=profile["environment_profile_id"],
        )
        self.assertEqual(corpus["record_count"], 60)
        self.assertEqual(len({item["record_id"] for item in corpus["records"]}), 60)
        self.assertEqual(provenance["corpus_id"], corpus["corpus_id"])
        self.assertEqual(len(provenance["shards"]), 60)
        damaged = json.loads(json.dumps(shards))
        damaged[0]["records"][0]["history_seed"] += 1
        with self.assertRaisesRegex(ValueError, "identity|history_seed"):
            control.collect_corpus_documents(
                sampling, "l5_feasibility", damaged,
                source_identity=source,
                environment_profile_id=profile["environment_profile_id"],
            )

    def test_primary_and_replay_execution_identities_are_separate(self):
        sampling = manifest.build_manifest(7)
        source = fake_source_identity()
        rows = self._l3_rows(sampling)
        corpus = execution_manifest.build_frozen_corpus(sampling, "l3_exact", rows)
        primary = execution_manifest.build_execution_manifest(
            sampling, corpus, source_identity=source, runtime_identity=fake_runtime(4)
        )
        replay = execution_manifest.build_execution_manifest(
            sampling, corpus, source_identity=source, runtime_identity=fake_runtime(1),
            execution_variant="blas1_replay",
        )
        self.assertNotEqual(primary["execution_manifest_id"], replay["execution_manifest_id"])
        self.assertEqual(sum(len(t["records"]) for t in primary["tasks"]), 600)
        self.assertEqual(sum(len(t["records"]) for t in replay["tasks"]), 30)
        self.assertEqual(replay["runtime_identity"]["blas_threads"], 1)
        self.assertTrue(replay["replay_comparison"]["required"])

    def test_slurm_schedule_preserves_manifest_parallelism(self):
        sampling = manifest.build_manifest(3)
        self.assertEqual(control.slurm_array_spec(sampling, "l3_exact"), "0-29%20")
        self.assertEqual(
            control.slurm_array_spec(sampling, "l5_feasibility", batch=0), "0-29%20"
        )
        self.assertEqual(
            control.slurm_array_spec(sampling, "l5_feasibility", batch=1), "30-59%20"
        )
        self.assertEqual(control.replay_array_spec(), "0-29%20")

    def test_hashed_gate_is_write_once_and_l5_requires_pass(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling = manifest.create_manifest(root, seed_base=9)
            repository_root = Path(__file__).resolve().parents[3]
            source = execution_manifest.capture_source_identity(repository_root)
            profile = control.build_environment_profile(
                source_identity=source,
                primary_runtime=fake_runtime(4),
                replay_runtime=fake_runtime(1),
                extended_runtime={"platform_machine": "test", "python_executable": "/x/python"},
            )
            control._write_once_json(
                root / control.PROFILE_NAME, profile, "environment profile"
            )
            all_rows = self._l3_rows(sampling)
            shards = []
            for task in manifest.tasks_for_stage(sampling, "l3_exact"):
                rows = [
                    {"history_seed": row["history_seed"], "record": row["record"]}
                    for row in all_rows if row["task_id"] == task["task_id"]
                ]
                shards.append(control.build_corpus_shard(
                    sampling, "l3_exact", task["task_id"], rows,
                    source_identity=source,
                    environment_profile_id=profile["environment_profile_id"],
                ))
            corpus, corpus_provenance = control.collect_corpus_documents(
                sampling, "l3_exact", shards, source_identity=source,
                environment_profile_id=profile["environment_profile_id"],
            )
            hpc_provenance = {
                "environment_profile_id": profile["environment_profile_id"],
                "corpus_provenance_id": corpus_provenance["corpus_provenance_id"],
                "l3_gate_id": None,
            }
            primary_root = control.stage_root(root, "l3_exact")
            replay_root = control.stage_root(root, "l3_exact", replay=True)
            for shard in shards:
                control._write_once_json(
                    primary_root / "corpus-shards" / f"task-{shard['task_id']:03d}.json",
                    shard,
                    "corpus shard",
                )
            execution_manifest.write_frozen_corpus(primary_root, corpus)
            execution_manifest.write_frozen_corpus(replay_root, corpus)
            control._write_once_json(
                primary_root / control.PROVENANCE_NAME,
                corpus_provenance,
                "corpus provenance",
            )
            control._write_once_json(
                replay_root / control.PROVENANCE_NAME,
                corpus_provenance,
                "replay corpus provenance",
            )
            primary = execution_manifest.create_execution_manifest(
                primary_root, sampling, source_identity=source,
                runtime_identity=fake_runtime(4), hpc_provenance=hpc_provenance,
            )
            replay = execution_manifest.create_execution_manifest(
                replay_root, sampling, source_identity=source,
                runtime_identity=fake_runtime(1), execution_variant="blas1_replay",
                hpc_provenance=hpc_provenance,
            )
            class ProbeExecutor:
                executor_mode = "persistent"

                def __init__(self):
                    self.identities = []
                    self.closed = False

                def preflight(self, identity):
                    self.identities.append(identity)

                def close(self):
                    self.closed = True

            probe = ProbeExecutor()
            with mock.patch.object(
                control, "verify_live_environment_profile", return_value=profile
            ):
                smoke = control.smoke_plan(
                    root, sampling, target="l3-primary", probe=True,
                    repository_root=repository_root, executor=probe,
                )
            self.assertEqual(
                smoke["runtime_probe"], "authenticated_worker_ready_no_contraction"
            )
            self.assertEqual(len(probe.identities), 1)
            self.assertTrue(probe.closed)
            criteria = {name: True for name in report.gate.REQUIRED_L3_CRITERIA}
            calibration_report = {
                "schema_version": 1,
                "kind": "adaptive_chi_l3_calibration_report",
                "sampling_manifest_id": sampling["manifest_id"],
                "primary_execution_manifest_id": primary["execution_manifest_id"],
                "replay_execution_manifest_id": replay["execution_manifest_id"],
                "source_bundle_id": primary["source_bundle_id"],
                "corpus_id": primary["corpus_id"],
                "corpus_provenance_id": hpc_provenance["corpus_provenance_id"],
                "environment_profile_id": profile["environment_profile_id"],
                "criteria": criteria,
                "terminal_complete": True,
                "passed": True,
                "report_source_sha256": source["files"][
                    "examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/report.py"
                ],
            }
            calibration_report["report_id"] = core.identity_hash(
                calibration_report, "report_id"
            )
            primary_root.mkdir(parents=True, exist_ok=True)
            (primary_root / report.REPORT_NAME).write_text(
                json.dumps(calibration_report), encoding="utf-8"
            )
            gate_document = report.build_l3_gate(
                sampling_manifest_id=sampling["manifest_id"],
                primary_execution_manifest_id=primary["execution_manifest_id"],
                replay_execution_manifest_id=replay["execution_manifest_id"],
                report_id=calibration_report["report_id"],
                passed=True,
                criteria=criteria,
                source_bundle_id=primary["source_bundle_id"],
                corpus_id=primary["corpus_id"],
                corpus_provenance_id=hpc_provenance["corpus_provenance_id"],
                environment_profile_id=profile["environment_profile_id"],
            )
            report.write_l3_gate(root, gate_document)
            self.assertEqual(
                report.require_l3_gate(root)["gate_id"], gate_document["gate_id"]
            )
            provenance_path = primary_root / control.PROVENANCE_NAME
            frozen_provenance = provenance_path.read_bytes()
            provenance_path.write_text("{}", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "corpus provenance"):
                report.require_l3_gate(root)
            provenance_path.write_bytes(frozen_provenance)
            l5_execution = {
                "stage": "l5_feasibility",
                "execution_manifest_id": "d" * 64,
                "hpc_provenance": {"l3_gate_id": gate_document["gate_id"]},
            }
            real_load_execution_manifest = execution_manifest.load_execution_manifest

            def load_execution_manifest(path, current_sampling):
                if Path(path) == control.stage_root(root, "l5_feasibility"):
                    return l5_execution
                return real_load_execution_manifest(path, current_sampling)

            with mock.patch.object(
                execution_manifest,
                "load_execution_manifest",
                side_effect=load_execution_manifest,
            ), mock.patch.object(
                control,
                "load_corpus_provenance",
                return_value={
                    "corpus_provenance_id": corpus_provenance["corpus_provenance_id"]
                },
            ):
                authenticated = control.require_l5_execution(root, sampling)
            self.assertEqual(authenticated["status"], "authenticated")
            self.assertEqual(
                authenticated["execution_manifest_id"], "d" * 64
            )
            with self.assertRaisesRegex(ValueError, "different L3 gate"):
                task_runner.verify_l5_gate_binding(
                    root / "l5-primary",
                    {
                        "stage": "l5_feasibility",
                        "hpc_provenance": {"l3_gate_id": "f" * 64},
                    },
                )
            with self.assertRaisesRegex(
                adaptive_runner.IdentityMismatchError, "missing.*HPC provenance"
            ):
                task_runner.verify_l5_gate_binding(
                    root / "l5-primary",
                    {"stage": "l5_feasibility", "hpc_provenance": None},
                    require_hpc=True,
                )
            with tempfile.TemporaryDirectory() as unrelated:
                unrelated_root = Path(unrelated)
                manifest.create_manifest(unrelated_root, seed_base=9)
                report.write_l3_gate(unrelated_root, gate_document)
                with self.assertRaisesRegex(FileNotFoundError, "profile|report|manifest"):
                    report.require_l3_gate(unrelated_root)
            changed = report.build_l3_gate(
                sampling_manifest_id=sampling["manifest_id"],
                primary_execution_manifest_id=primary["execution_manifest_id"],
                replay_execution_manifest_id=replay["execution_manifest_id"],
                report_id=calibration_report["report_id"],
                passed=False,
                criteria={**criteria, "primary_exact_gate": False},
                source_bundle_id=primary["source_bundle_id"],
                corpus_id=primary["corpus_id"],
                corpus_provenance_id=hpc_provenance["corpus_provenance_id"],
                environment_profile_id=profile["environment_profile_id"],
            )
            with self.assertRaisesRegex(ValueError, "different|identity"):
                report.write_l3_gate(root, changed)
            with self.assertRaisesRegex(ValueError, "criteria"):
                report.build_l3_gate(
                    sampling_manifest_id=sampling["manifest_id"],
                    primary_execution_manifest_id=primary["execution_manifest_id"],
                    replay_execution_manifest_id=replay["execution_manifest_id"],
                    report_id=calibration_report["report_id"], passed=True,
                    criteria={"primary_exact_gate": True},
                    source_bundle_id=primary["source_bundle_id"],
                    corpus_id=primary["corpus_id"],
                    corpus_provenance_id=hpc_provenance["corpus_provenance_id"],
                    environment_profile_id=profile["environment_profile_id"],
                )

    def test_julia_corpus_envelope_binds_history_seed_order_and_rng(self):
        sampling = manifest.build_manifest(17)
        plan = control.plan_corpus_task(
            Path("/study"), sampling, stage="l3_exact", task_id=0
        )
        records = [fake_record(index, index) for index in range(1, 21)]
        envelope = {
            "schema_version": 1,
            "kind": "adaptive_chi_julia_corpus_task",
            "L": plan["L"],
            "p": float(plan["p_key"]),
            "p_key": plan["p_key"],
            "task_id": plan["task_id"],
            "seed_mode": plan["seed_mode"],
            "rng_algorithm": "Random.Xoshiro",
            "runtime_identity": {
                key: value for key, value in fake_runtime(4).items()
                if key != "python_version"
            },
            "history_seeds": plan["history_seeds"],
            "records": records,
        }
        rows = control._rows_from_julia(envelope, plan, fake_runtime(4))
        self.assertEqual([row["history_seed"] for row in rows], plan["history_seeds"])
        damaged = json.loads(json.dumps(envelope))
        damaged["history_seeds"][0], damaged["history_seeds"][1] = (
            damaged["history_seeds"][1], damaged["history_seeds"][0]
        )
        with self.assertRaisesRegex(ValueError, "history seed"):
            control._rows_from_julia(damaged, plan, fake_runtime(4))
        damaged = json.loads(json.dumps(envelope))
        damaged["rng_algorithm"] = "MersenneTwister"
        with self.assertRaisesRegex(ValueError, "envelope|RNG"):
            control._rows_from_julia(damaged, plan, fake_runtime(4))
        damaged = json.loads(json.dumps(envelope))
        damaged["runtime_identity"]["blas_threads"] = 1
        with self.assertRaisesRegex(ValueError, "runtime identity"):
            control._rows_from_julia(damaged, plan, fake_runtime(4))

    def test_live_environment_verifier_rejects_python_lock_drift(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            repository = Path(__file__).resolve().parents[3]
            source = execution_manifest.capture_source_identity(repository)
            packages = {
                "numpy": "1", "scipy": "2", "pymatching": "3", "ldpc": "4"
            }
            freeze_text = "numpy==1\n"
            profile = control.build_environment_profile(
                source_identity=source,
                primary_runtime=fake_runtime(4),
                replay_runtime=fake_runtime(1),
                extended_runtime={
                    "platform_machine": platform.machine(),
                    "platform_system": platform.system(),
                    "python_executable": str(Path(sys.executable).resolve()),
                    "python_implementation": platform.python_implementation(),
                    "pip_freeze_sha256": hashlib.sha256(freeze_text.encode()).hexdigest(),
                    "repository_realpath": str(repository),
                    "julia_executable": "/x/julia",
                    "python_packages": packages,
                },
            )
            control._write_once_json(root / control.PROFILE_NAME, profile, "profile")
            completed = subprocess.CompletedProcess([], 0, freeze_text, "")
            with mock.patch.object(control.subprocess, "run", return_value=completed), \
                 mock.patch.object(control, "_package_versions", return_value=packages), \
                 mock.patch.object(control.platform, "python_version", return_value="3.12.14"):
                self.assertEqual(
                    control.verify_live_environment_profile(
                        root, repository_root=repository, julia="/x/julia"
                    )["environment_profile_id"],
                    profile["environment_profile_id"],
                )
                completed.stdout = "numpy==9\n"
                with self.assertRaisesRegex(ValueError, "live Python"):
                    control.verify_live_environment_profile(
                        root, repository_root=repository, julia="/x/julia"
                    )
            with mock.patch.object(
                control.platform, "python_version", return_value="3.13.0"
            ):
                with self.assertRaisesRegex(ValueError, "Python version"):
                    control.verify_live_environment_profile(
                        root, repository_root=repository, julia="/x/julia"
                    )

    def test_clopper_pearson_boundary_convention(self):
        expected = 1.0 - 0.05 ** (1.0 / 600.0)
        self.assertAlmostEqual(report.clopper_pearson_upper(0, 600), expected, places=12)
        self.assertEqual(report.clopper_pearson_upper(600, 600), 1.0)
        self.assertIsNone(report.clopper_pearson_upper(0, 0))

    def test_checkpoint_summary_charges_worker_startup_to_l5_wall(self):
        with tempfile.TemporaryDirectory() as temporary:
            record_dir = Path(temporary)
            (record_dir / "chi_0032.json").write_text(
                json.dumps(
                    {
                        "elapsed_seconds": 10.0,
                        "worker_startup_seconds": 7.0,
                        "decoder_build_seconds": 2.0,
                        "decode_seconds": 8.0,
                        "decoder_cache_hit": False,
                        "peak_rss_bytes": 2048,
                    }
                ),
                encoding="utf-8",
            )
            summary = report._checkpoint_summary(
                [
                    {
                        "record_directory": str(record_dir),
                        "decision": {"tested_chi": [32]},
                    }
                ]
            )
            self.assertEqual(summary["worker_startup_seconds"]["median"], 7.0)
            self.assertEqual(summary["record_total_elapsed_seconds"], [17.0])

    def test_l5_anchor_completion_rule_rejects_unknown_anchor(self):
        anchors = [
            {
                "final_status": "accepted",
                "adaptive_decision": {
                    "anchor_completed": True,
                    "anchor_final_chi": 480,
                    "anchor_class_changed": False,
                },
            }
            for _ in range(24)
        ]
        self.assertTrue(report._all_forced_anchors_complete(anchors, final_chi=480))
        anchors[0]["final_status"] = "unknown"
        anchors[0]["adaptive_decision"]["anchor_completed"] = False
        self.assertFalse(report._all_forced_anchors_complete(anchors, final_chi=480))
        anchors[0]["final_status"] = "accepted"
        anchors[0]["adaptive_decision"]["anchor_completed"] = True
        anchors[0]["adaptive_decision"]["anchor_class_changed"] = True
        self.assertFalse(report._all_forced_anchors_complete(anchors, final_chi=480))


if __name__ == "__main__":
    unittest.main()
