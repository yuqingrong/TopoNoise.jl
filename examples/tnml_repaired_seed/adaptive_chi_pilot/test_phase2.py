#!/usr/bin/env python3
"""Bounded phase-2 tests. No Julia process or contraction is launched."""

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

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import adaptive_runner
import core
import execution_manifest
import manifest


def source_identity():
    return {
        "project_sha256": "1" * 64,
        "manifest_sha256": "2" * 64,
        "package_source_tree_sha256": "c" * 64,
        "files": {
            "src/TopoNoise.jl": "3" * 64,
            "examples/tnml_repaired_seed/model.jl": "4" * 64,
            "examples/tnml_repaired_seed/run.jl": "5" * 64,
            "examples/tnml_repaired_seed/dense_frontier.jl": "6" * 64,
            "examples/tnml_repaired_seed/adaptive_chi_pilot/core.py": "7" * 64,
            "examples/tnml_repaired_seed/adaptive_chi_pilot/manifest.py": "8" * 64,
            "examples/tnml_repaired_seed/adaptive_chi_pilot/execution_manifest.py": "9" * 64,
            "examples/tnml_repaired_seed/adaptive_chi_pilot/adaptive_runner.py": "a" * 64,
            "examples/tnml_repaired_seed/adaptive_chi_pilot/julia_adapter.jl": "b" * 64,
            "examples/tnml_repaired_seed/adaptive_chi_pilot/task_runner.py": "d" * 64,
            "examples/tnml_repaired_seed/adaptive_chi_pilot/persistent_worker.jl": "e" * 64,
        },
    }


def runtime_identity():
    return {
        "julia_version": "1.12.4",
        "python_version": "3.12.14",
        "blas_vendor": "openblas64",
        "cpu_model": "test-aarch64",
        "julia_threads": 1,
        "blas_threads": 4,
    }


def corpus_rows(sampling, stage):
    rows = []
    counter = 1 if stage == "l3_exact" else 1001
    for task in manifest.tasks_for_stage(sampling, stage):
        for shot_index, history_seed in enumerate(task["history_seeds"], start=1):
            record_id = hashlib.sha256(
                f"{stage}:{task['task_id']}:{shot_index}:{counter}".encode()
            ).hexdigest()
            record = {
                "shot_index": shot_index,
                "raw_fault_count": 12,
                "raw_fault_indices": [1] if counter % 2 else [],
                "seed_syndrome": [bool(counter % 2), False],
                "seed_recovery": [bool(counter % 2), False],
                "recovery_mask": counter % 4,
                "final_syndrome": [False, bool(counter % 3)],
                "true_mask": (counter + 1) % 4,
                "mwpm_failure": bool(counter % 5 == 0),
                "raw_seed_mask": counter % 4,
                "transfer_mask": (counter + 2) % 4,
                "record_id": record_id,
            }
            rows.append(
                {
                    "task_id": task["task_id"],
                    "shot_index": shot_index,
                    "history_seed": history_seed,
                    "record": record,
                }
            )
            counter += 1
    return rows


def make_bundle(root: Path, stage="l3_exact"):
    sampling = manifest.build_manifest(202610600000)
    (root / "sampling-manifest.json").write_text(
        json.dumps(sampling, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    corpus = execution_manifest.build_frozen_corpus(
        sampling, stage, corpus_rows(sampling, stage)
    )
    execution_manifest.write_frozen_corpus(root, corpus)
    execution = execution_manifest.create_execution_manifest(
        root,
        sampling,
        source_identity=source_identity(),
        runtime_identity=runtime_identity(),
    )
    return sampling, corpus, execution


def tn_result(chi, *, probabilities=None, log_joint=None, status="complete"):
    if status != "complete":
        return {
            "status": status,
            "valid": False,
            "elapsed_seconds": 0.25,
            "peak_rss_bytes": 1024,
        }
    probabilities = probabilities or {
        32: [0.7, 0.1, 0.1, 0.1],
        64: [0.7000004, 0.0999998, 0.0999999, 0.0999999],
        128: [0.7000005, 0.09999975, 0.099999875, 0.099999875],
    }.get(chi, [0.70000055, 0.099999725, 0.0999998625, 0.0999998625])
    default_log_joint = {
        32: -10.0,
        64: -9.9999996,
        128: -9.9999995,
        240: -9.99999945,
        480: -9.999999425,
        960: -9.9999994125,
    }[chi]
    return {
        "status": "complete",
        "valid": True,
        "probabilities": probabilities,
        "logical_mask": max(range(4), key=probabilities.__getitem__),
        "log_joint_syndrome": log_joint if log_joint is not None else default_log_joint,
        "elapsed_seconds": 0.25,
        "peak_rss_bytes": 1024,
        "rss_scope": "process_lifetime_high_water",
    }


def exact_result(*, probabilities=None, log_joint=None, status="complete"):
    if status != "complete":
        return {
            "status": status,
            "valid": False,
            "elapsed_seconds": 0.1,
            "peak_rss_bytes": 512,
        }
    probabilities = probabilities or [0.7000005, 0.09999975, 0.099999875, 0.099999875]
    return {
        "status": "complete",
        "valid": True,
        "probabilities": probabilities,
        "logical_mask": max(range(4), key=probabilities.__getitem__),
        "log_joint_syndrome": log_joint if log_joint is not None else -9.9999995,
        "elapsed_seconds": 0.1,
        "peak_rss_bytes": 512,
        "rss_scope": "process_lifetime_high_water",
    }


class FakeExecutor:
    executor_mode = "one_shot"

    def __init__(self, *, tn_factory=tn_result, exact_factory=exact_result):
        self.tn_factory = tn_factory
        self.exact_factory = exact_factory
        self.tn_calls = []
        self.exact_calls = []

    def run_tn(self, record, *, chi, cutoff, execution_identity):
        self.tn_calls.append((record["record_id"], chi, cutoff, execution_identity))
        return self.tn_factory(chi)

    def run_exact(self, record, *, execution_identity):
        self.exact_calls.append((record["record_id"], execution_identity))
        return self.exact_factory()


class InterruptOnceExecutor(FakeExecutor):
    def __init__(self):
        super().__init__()
        self.interrupted = False

    def run_tn(self, record, *, chi, cutoff, execution_identity):
        self.tn_calls.append((record["record_id"], chi, cutoff, execution_identity))
        if chi == 64 and not self.interrupted:
            self.interrupted = True
            raise KeyboardInterrupt()
        return self.tn_factory(chi)


class ExecutionManifestTests(unittest.TestCase):
    def test_duplicate_julia_history_hashes_get_unique_phase2_record_identities(self):
        sampling = manifest.build_manifest(202610600000)
        rows = corpus_rows(sampling, "l3_exact")
        candidates = [
            row
            for row in rows
            if row["shot_index"] == 1
            and manifest.task_for_id(sampling, "l3_exact", row["task_id"])["p_key"]
            == "0.0060"
        ]
        self.assertGreaterEqual(len(candidates), 2)
        candidates[1]["record"] = json.loads(json.dumps(candidates[0]["record"]))
        corpus = execution_manifest.build_frozen_corpus(
            sampling, "l3_exact", rows
        )
        duplicated = [
            entry
            for entry in corpus["records"]
            if entry["julia_record_id"] == candidates[0]["record"]["record_id"]
        ]
        self.assertEqual(len(duplicated), 2)
        self.assertEqual(len({entry["record_id"] for entry in duplicated}), 2)
        execution = execution_manifest.build_execution_manifest(
            sampling,
            corpus,
            source_identity=source_identity(),
            runtime_identity=runtime_identity(),
        )
        ranked = {
            item["record_id"]
            for item in execution["selection"]["0.0060"]
        }
        self.assertTrue(ranked <= {entry["record_id"] for entry in corpus["records"]})

    def test_source_identity_capture_and_executor_verification_detect_drift(self):
        with tempfile.TemporaryDirectory() as temporary:
            repository = Path(temporary)
            paths = [
                "Project.toml",
                "Manifest.toml",
                "src/TopoNoise.jl",
                "examples/tnml_repaired_seed/model.jl",
                "examples/tnml_repaired_seed/run.jl",
                "examples/tnml_repaired_seed/dense_frontier.jl",
                "examples/tnml_repaired_seed/adaptive_chi_pilot/core.py",
                "examples/tnml_repaired_seed/adaptive_chi_pilot/manifest.py",
                "examples/tnml_repaired_seed/adaptive_chi_pilot/execution_manifest.py",
                "examples/tnml_repaired_seed/adaptive_chi_pilot/adaptive_runner.py",
                "examples/tnml_repaired_seed/adaptive_chi_pilot/julia_adapter.jl",
                "examples/tnml_repaired_seed/adaptive_chi_pilot/task_runner.py",
                "examples/tnml_repaired_seed/adaptive_chi_pilot/persistent_worker.jl",
            ]
            for index, relative in enumerate(paths):
                path = repository / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(f"source-{index}\n", encoding="utf-8")
            captured = execution_manifest.capture_source_identity(repository)
            self.assertEqual(set(captured["files"]), set(paths[2:]))
            executor = adaptive_runner.JuliaExecutor(repository_root=repository)
            identity = {
                "source_identity": captured,
                "runtime_identity": {
                    "julia_version": "1.12.4",
                    "python_version": platform.python_version(),
                    "blas_vendor": "lbt",
                    "cpu_model": "test-cpu",
                    "julia_threads": 1,
                    "blas_threads": 4,
                },
            }
            executor.verify_execution_environment(identity)
            (repository / "examples/tnml_repaired_seed/run.jl").write_text(
                "changed\n", encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "source identity"):
                executor.verify_execution_environment(identity)

    def test_julia_runtime_evidence_must_match_frozen_profile(self):
        expected = {
            "julia_version": "1.12.4",
            "python_version": platform.python_version(),
            "blas_vendor": "lbt",
            "cpu_model": "arm-test",
            "julia_threads": 1,
            "blas_threads": 4,
        }
        actual = {key: value for key, value in expected.items() if key != "python_version"}
        adaptive_runner.JuliaExecutor.verify_julia_runtime(actual, expected)
        actual["blas_threads"] = 8
        with self.assertRaisesRegex(ValueError, "runtime identity"):
            adaptive_runner.JuliaExecutor.verify_julia_runtime(actual, expected)

    def test_l3_manifest_binds_corpus_identity_replays_sources_runtime_and_cp_rule(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling, corpus, execution = make_bundle(root, "l3_exact")
            self.assertEqual(execution["sampling_manifest_id"], sampling["manifest_id"])
            self.assertEqual(execution["corpus_id"], corpus["corpus_id"])
            self.assertEqual(execution["source_identity"], source_identity())
            self.assertEqual(execution["runtime_identity"], runtime_identity())
            self.assertEqual(
                execution["executor_policy"],
                {
                    "persistent_protocol": "persistent_toml_mailbox_v1",
                    "protocol_version": 1,
                    "scheduling": "serial_one_worker_per_L_p",
                    "decoder_cache": "one_decoder_per_chi_fixed_L_p_cutoff",
                    "allowed_executors": ["one_shot", "persistent"],
                    "timing_scopes": {
                        "one_shot": {
                            "wall_time_scope": "whole_julia_process",
                            "decoder_cache_scope": "none",
                        },
                        "persistent": {
                            "wall_time_scope": "request_after_worker_ready",
                            "decoder_cache_scope": "worker_lifetime_by_chi",
                        },
                    },
                },
            )
            self.assertEqual(execution["chi_ladder"], [32, 64, 128, 240, 480, 960])
            self.assertEqual(execution["exact_audit"]["ml_tie_tolerance"], 1e-9)
            self.assertEqual(execution["exact_audit"]["posterior_linf_gate"], 1e-5)
            self.assertEqual(execution["exact_audit"]["log_joint_gate"], 1e-5)
            cp = execution["clopper_pearson"]
            self.assertEqual(cp["confidence_level"], 0.95)
            self.assertEqual(cp["sided"], "one_sided_upper")
            self.assertEqual(cp["denominator"], "accepted_with_complete_exact_reference")
            resources = execution["resource_policy"]
            self.assertEqual(resources["max_rss_bytes"], 16 * 2**30)
            self.assertEqual(
                resources["wall_seconds_by_chi"],
                {
                    "32": 120,
                    "64": 300,
                    "128": 900,
                    "240": 1800,
                    "480": 3600,
                    "960": 7200,
                },
            )
            self.assertEqual(resources["exact_wall_seconds"], 300)
            self.assertEqual(resources["preflight_wall_seconds"], 300)
            self.assertEqual(
                resources["enforcement"],
                "wall_timeout_internal_memory_cap_external_with_posthoc_rejection",
            )
            selected = [record for task in execution["tasks"] for record in task["records"]
                        if record["replay_rank"] is not None]
            self.assertEqual(len(selected), 30)
            for p_key in manifest.P_VALUES:
                ranks = sorted(
                    record["replay_rank"] for record in selected if record["p_key"] == p_key
                )
                self.assertEqual(ranks, [1, 2, 3, 4, 5])
            self.assertTrue(all(task["exact_required"] for task in execution["tasks"]))
            self.assertEqual(
                execution_manifest.load_execution_manifest(root, sampling), execution
            )

    def test_l5_manifest_selects_four_forced_anchors_per_stratum(self):
        with tempfile.TemporaryDirectory() as temporary:
            _, _, execution = make_bundle(Path(temporary), "l5_feasibility")
            selected = [record for task in execution["tasks"] for record in task["records"]
                        if record["anchor_rank"] is not None]
            self.assertEqual(len(selected), 24)
            for p_key in manifest.P_VALUES:
                ranks = sorted(
                    record["anchor_rank"] for record in selected if record["p_key"] == p_key
                )
                self.assertEqual(ranks, [1, 2, 3, 4])
            self.assertTrue(all(not task["exact_required"] for task in execution["tasks"]))
            self.assertEqual(
                execution["resource_policy"]["wall_seconds_by_chi"],
                {"32": 120, "64": 300, "128": 900, "240": 1800, "480": 3600},
            )
            identity = adaptive_runner._execution_identity(
                execution, execution["tasks"][0], execution["tasks"][0]["records"][0]
            )
            self.assertEqual(
                adaptive_runner.JuliaExecutor._wall_timeout("tn", 32, identity),
                120.0,
            )
            self.assertEqual(
                adaptive_runner.JuliaExecutor._wall_timeout("exact", 1, identity),
                300.0,
            )

    def test_corpus_and_execution_identity_reject_tampering_and_rewrite(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling, corpus, execution = make_bundle(root)
            self.assertEqual(
                execution_manifest.create_execution_manifest(
                    root,
                    sampling,
                    source_identity=source_identity(),
                    runtime_identity=runtime_identity(),
                ),
                execution,
            )
            tampered = json.loads(json.dumps(execution))
            tampered["tasks"][0]["records"][0]["record_id"] = "f" * 64
            tampered["execution_manifest_id"] = core.identity_hash(
                tampered, "execution_manifest_id"
            )
            with self.assertRaisesRegex(ValueError, "derived|identity"):
                execution_manifest.validate_execution_manifest(tampered, sampling, corpus)
            bad_corpus = json.loads(json.dumps(corpus))
            bad_corpus["records"][0]["record"]["true_mask"] ^= 1
            with self.assertRaisesRegex(ValueError, "identity"):
                execution_manifest.validate_frozen_corpus(bad_corpus, sampling)


class AdaptiveRunnerTests(unittest.TestCase):
    def first_record(self, execution):
        return execution["tasks"][0]["records"][0]

    def test_julia_nonzero_exit_accepts_only_explicit_memory_limit_output(self):
        expected_runtime = runtime_identity()
        actual_runtime = {
            key: expected_runtime[key]
            for key in (
                "julia_version",
                "blas_vendor",
                "cpu_model",
                "julia_threads",
                "blas_threads",
            )
        }
        identity = {
            "L": 3,
            "p": 0.006,
            "runtime_identity": expected_runtime,
            "resource_policy": {
                "max_rss_bytes": 16 * 2**30,
                "wall_seconds_by_chi": {},
                "exact_wall_seconds": None,
            },
        }
        executor = adaptive_runner.JuliaExecutor(
            repository_root=Path(__file__).resolve().parents[3], julia="julia"
        )
        executor.preflight = lambda _identity: None

        def completed(command, status):
            output = Path(command[7])
            lines = [
                f'status = "{status}"',
                f"valid = {'true' if status == 'complete' else 'false'}",
                "elapsed_seconds = 0.1",
                "peak_rss_bytes = 1024",
            ]
            if status == "complete":
                lines.extend(
                    [
                        "probabilities = [0.7, 0.1, 0.1, 0.1]",
                        "logical_mask = 0",
                        "log_joint_syndrome = -10.0",
                    ]
                )
            lines.extend(["", "[runtime_identity]"])
            for key, value in actual_runtime.items():
                encoded = f'"{value}"' if isinstance(value, str) else str(value)
                lines.append(f"{key} = {encoded}")
            output.write_text("\n".join(lines) + "\n", encoding="utf-8")
            return subprocess.CompletedProcess(command, 1, stdout="", stderr="boom")

        with mock.patch.object(
            adaptive_runner.subprocess,
            "run",
            side_effect=lambda command, **_kwargs: completed(command, "memory_limit"),
        ):
            result = executor._invoke(
                "tn",
                {"record_id": "a" * 64},
                chi=32,
                cutoff=1e-18,
                execution_identity=identity,
            )
        self.assertEqual(result["status"], "memory_limit")
        self.assertFalse(result["valid"])

        for label, status in (
            ("reported_complete", "complete"),
            ("reported_exception", "exception"),
            ("no_output", None),
        ):
            with self.subTest(label=label), mock.patch.object(
                adaptive_runner.subprocess,
                "run",
                side_effect=(
                    (lambda command, **_kwargs: completed(command, status))
                    if status is not None
                    else (
                        lambda command, **_kwargs: subprocess.CompletedProcess(
                            command, 1, stdout="", stderr="boom"
                        )
                    )
                ),
            ):
                result = executor._invoke(
                    "tn",
                    {"record_id": "a" * 64},
                    chi=32,
                    cutoff=1e-18,
                    execution_identity=identity,
                )
                self.assertEqual(result["status"], "exception")
                self.assertFalse(result["valid"])

    def test_stable_record_accepts_at_128_audits_exact_and_resume_calls_nothing(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling, _, execution = make_bundle(root)
            record = self.first_record(execution)
            executor = FakeExecutor()
            decision = adaptive_runner.run_record(
                root, sampling, task_id=0, record_id=record["record_id"], executor=executor
            )
            self.assertEqual([call[1] for call in executor.tn_calls], [32, 64, 128])
            self.assertEqual(len(executor.exact_calls), 1)
            self.assertEqual(decision["final_status"], "accepted")
            self.assertEqual(decision["adaptive_decision"]["accepted_chi"], 128)
            self.assertTrue(decision["exact_audit"]["passed"])
            self.assertFalse(decision["exact_audit"]["false_accept"])
            record_dir = Path(decision["record_directory"])
            self.assertTrue((record_dir / "chi_0032.json").is_file())
            self.assertTrue((record_dir / "chi_0064.json").is_file())
            self.assertTrue((record_dir / "chi_0128.json").is_file())
            self.assertTrue((record_dir / "exact.json").is_file())
            self.assertTrue((record_dir / "decision.json").is_file())
            first_level = json.loads(
                (record_dir / "chi_0032.json").read_text(encoding="utf-8")
            )
            exact_checkpoint = json.loads(
                (record_dir / "exact.json").read_text(encoding="utf-8")
            )
            for checkpoint in (first_level, exact_checkpoint):
                self.assertEqual(checkpoint["executor_mode"], "one_shot")
                self.assertEqual(
                    checkpoint["wall_time_scope"], "whole_julia_process"
                )
                self.assertEqual(checkpoint["decoder_cache_scope"], "none")

            resumed = FakeExecutor()
            self.assertEqual(
                adaptive_runner.run_record(
                    root,
                    sampling,
                    task_id=0,
                    record_id=record["record_id"],
                    executor=resumed,
                ),
                decision,
            )
            self.assertEqual(resumed.tn_calls, [])
            self.assertEqual(resumed.exact_calls, [])

            class PersistentModeFake(FakeExecutor):
                executor_mode = "persistent"

            with self.assertRaisesRegex(ValueError, "executor.*mismatch"):
                adaptive_runner.run_record(
                    root,
                    sampling,
                    task_id=0,
                    record_id=record["record_id"],
                    executor=PersistentModeFake(),
                )

    def test_worker_lifetime_rss_scope_is_preserved_in_exact_and_tn_checkpoints(self):
        def scoped_tn(chi):
            result = tn_result(chi)
            result["rss_scope"] = "worker_lifetime_high_water"
            result.update(
                decoder_cache_hit=False,
                decoder_build_seconds=0.05,
                decode_seconds=0.20,
            )
            return result

        def scoped_exact():
            result = exact_result()
            result["rss_scope"] = "worker_lifetime_high_water"
            result.update(
                decoder_cache_hit=False,
                decoder_build_seconds=0.0,
                decode_seconds=0.1,
            )
            return result

        class PersistentScopedFake(FakeExecutor):
            executor_mode = "persistent"

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling, _, execution = make_bundle(root)
            record = self.first_record(execution)
            decision = adaptive_runner.run_record(
                root,
                sampling,
                task_id=0,
                record_id=record["record_id"],
                executor=PersistentScopedFake(
                    tn_factory=scoped_tn, exact_factory=scoped_exact
                ),
            )
            record_dir = Path(decision["record_directory"])
            level = json.loads((record_dir / "chi_0032.json").read_text())
            exact = json.loads((record_dir / "exact.json").read_text())
            self.assertEqual(
                level["rss_scope"], "worker_lifetime_high_water"
            )
            self.assertEqual(
                exact["rss_scope"], "worker_lifetime_high_water"
            )
            self.assertFalse(level["decoder_cache_hit"])
            self.assertEqual(level["decoder_build_seconds"], 0.05)
            self.assertEqual(level["decode_seconds"], 0.20)

    def test_exact_disagreement_is_recorded_as_false_accept_not_hidden_unknown(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling, _, execution = make_bundle(root)
            record = self.first_record(execution)
            executor = FakeExecutor(
                exact_factory=lambda: exact_result(
                    probabilities=[0.1, 0.7, 0.1, 0.1], log_joint=-9.0
                )
            )
            decision = adaptive_runner.run_record(
                root, sampling, task_id=0, record_id=record["record_id"], executor=executor
            )
            self.assertEqual(decision["final_status"], "accepted")
            self.assertFalse(decision["exact_audit"]["passed"])
            self.assertTrue(decision["exact_audit"]["false_accept"])
            self.assertFalse(decision["exact_audit"]["ml_tie_aware_match"])
            self.assertGreater(decision["exact_audit"]["posterior_linf"], 1e-5)
            self.assertGreater(decision["exact_audit"]["log_joint_error"], 1e-5)

    def test_timeout_or_missing_exact_propagates_unknown(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling, _, execution = make_bundle(root)
            record = self.first_record(execution)
            timeout_executor = FakeExecutor(
                tn_factory=lambda chi: tn_result(chi, status="timeout") if chi == 64 else tn_result(chi)
            )
            decision = adaptive_runner.run_record(
                root,
                sampling,
                task_id=0,
                record_id=record["record_id"],
                executor=timeout_executor,
            )
            self.assertEqual(decision["final_status"], "unknown")
            self.assertEqual(decision["unknown_reason"], "timeout")
            self.assertEqual([call[1] for call in timeout_executor.tn_calls], [32, 64])
            self.assertNotIn("failure", decision)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling, _, execution = make_bundle(root)
            record = self.first_record(execution)
            exact_timeout = FakeExecutor(
                exact_factory=lambda: exact_result(status="timeout")
            )
            decision = adaptive_runner.run_record(
                root,
                sampling,
                task_id=0,
                record_id=record["record_id"],
                executor=exact_timeout,
            )
            self.assertEqual(decision["final_status"], "unknown")
            self.assertEqual(decision["unknown_reason"], "exact_timeout")
            self.assertEqual(len(exact_timeout.exact_calls), 1)
            self.assertEqual(exact_timeout.tn_calls, [])

    def test_memory_and_malformed_complete_output_are_checkpointed_unknown(self):
        def over_cap(chi):
            result = tn_result(chi)
            result["peak_rss_bytes"] = 16 * 2**30 + 1
            return result

        for label, factory, expected_reason in (
            (
                "memory",
                lambda chi: tn_result(chi, status="memory_limit"),
                "memory_limit",
            ),
            (
                "malformed",
                lambda chi: {
                    "status": "complete",
                    "valid": True,
                    "probabilities": [float("nan"), 0.1, 0.1, 0.1],
                    "logical_mask": 0,
                    "log_joint_syndrome": -1.0,
                    "elapsed_seconds": 0.1,
                    "peak_rss_bytes": 1,
                },
                "invalid",
            ),
            ("posthoc_rss_cap", over_cap, "memory_limit"),
        ):
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                sampling, _, execution = make_bundle(root)
                record = self.first_record(execution)
                decision = adaptive_runner.run_record(
                    root,
                    sampling,
                    task_id=0,
                    record_id=record["record_id"],
                    executor=FakeExecutor(tn_factory=factory),
                )
                self.assertEqual(decision["final_status"], "unknown")
                self.assertEqual(decision["unknown_reason"], expected_reason)
                checkpoint = json.loads(
                    (Path(decision["record_directory"]) / "chi_0032.json").read_text()
                )
                self.assertEqual(checkpoint["status"], expected_reason)
                self.assertFalse(checkpoint["valid"])

    def test_partial_prefix_resume_does_not_recompute_saved_exact_or_chi(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling, _, execution = make_bundle(root)
            record = self.first_record(execution)
            interrupted = InterruptOnceExecutor()
            with self.assertRaises(KeyboardInterrupt):
                adaptive_runner.run_record(
                    root,
                    sampling,
                    task_id=0,
                    record_id=record["record_id"],
                    executor=interrupted,
                )
            self.assertEqual([call[1] for call in interrupted.tn_calls], [32, 64])
            resumed = FakeExecutor()
            decision = adaptive_runner.run_record(
                root,
                sampling,
                task_id=0,
                record_id=record["record_id"],
                executor=resumed,
            )
            self.assertEqual([call[1] for call in resumed.tn_calls], [64, 128])
            self.assertEqual(resumed.exact_calls, [])
            self.assertEqual(decision["final_status"], "accepted")

    def test_resume_rejects_phase2_identity_tampering_and_checkpoint_gaps(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling, _, execution = make_bundle(root)
            record = self.first_record(execution)
            decision = adaptive_runner.run_record(
                root,
                sampling,
                task_id=0,
                record_id=record["record_id"],
                executor=FakeExecutor(),
            )
            record_dir = Path(decision["record_directory"])
            first_path = record_dir / "chi_0032.json"
            first = json.loads(first_path.read_text())
            first["runtime_profile_id"] = "0" * 64
            first["checkpoint_id"] = core.identity_hash(first, "checkpoint_id")
            first_path.write_text(json.dumps(first, indent=2, sort_keys=True) + "\n")
            with self.assertRaisesRegex(ValueError, "phase-2 identity"):
                adaptive_runner.run_record(
                    root,
                    sampling,
                    task_id=0,
                    record_id=record["record_id"],
                    executor=FakeExecutor(),
                )

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling, _, execution = make_bundle(root)
            record = self.first_record(execution)
            decision = adaptive_runner.run_record(
                root,
                sampling,
                task_id=0,
                record_id=record["record_id"],
                executor=FakeExecutor(),
            )
            record_dir = Path(decision["record_directory"])
            (record_dir / "decision.json").unlink()
            (record_dir / "chi_0064.json").unlink()
            with self.assertRaisesRegex(ValueError, "contiguous chi prefix"):
                adaptive_runner.run_record(
                    root,
                    sampling,
                    task_id=0,
                    record_id=record["record_id"],
                    executor=FakeExecutor(),
                )

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling, _, execution = make_bundle(root)
            record = self.first_record(execution)
            decision = adaptive_runner.run_record(
                root,
                sampling,
                task_id=0,
                record_id=record["record_id"],
                executor=FakeExecutor(),
            )
            decision_path = Path(decision["record_directory"]) / "decision.json"
            tampered = json.loads(decision_path.read_text())
            tampered["final_status"] = "unknown"
            tampered["unknown_reason"] = "forged"
            tampered.pop("failure", None)
            tampered.pop("logical_mask", None)
            tampered["decision_id"] = core.identity_hash(tampered, "decision_id")
            decision_path.write_text(
                json.dumps(tampered, indent=2, sort_keys=True) + "\n"
            )
            with self.assertRaisesRegex(ValueError, "derived identity"):
                adaptive_runner.run_record(
                    root,
                    sampling,
                    task_id=0,
                    record_id=record["record_id"],
                    executor=FakeExecutor(),
                )

    def test_exact_tie_boundary_is_inclusive(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling, _, execution = make_bundle(root)
            record = self.first_record(execution)
            stable = [0.45000002, 0.44999998, 0.05, 0.05]
            executor = FakeExecutor(
                tn_factory=lambda chi: tn_result(
                    chi, probabilities=stable, log_joint=-8.0
                ),
                exact_factory=lambda: exact_result(
                    probabilities=[0.45, 0.450000001, 0.05, 0.049999999],
                    log_joint=-8.0,
                ),
            )
            decision = adaptive_runner.run_record(
                root,
                sampling,
                task_id=0,
                record_id=record["record_id"],
                executor=executor,
            )
            self.assertEqual(decision["final_status"], "accepted")
            self.assertIn(0, decision["exact_audit"]["exact_maximizers"])
            self.assertTrue(decision["exact_audit"]["ml_tie_aware_match"])
            self.assertTrue(decision["exact_audit"]["passed"])

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling, _, execution = make_bundle(root)
            record = self.first_record(execution)
            executor = FakeExecutor(exact_factory=lambda: exact_result(status="exception"))
            decision = adaptive_runner.run_record(
                root, sampling, task_id=0, record_id=record["record_id"], executor=executor
            )
            self.assertEqual(decision["final_status"], "unknown")
            self.assertEqual(decision["unknown_reason"], "invalid_exact_reference")
            self.assertEqual(decision["adaptive_decision"]["status"], "unknown")
            self.assertEqual(executor.tn_calls, [])
            self.assertNotIn("failure", decision)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling, _, execution = make_bundle(root)
            record = self.first_record(execution)
            malformed_exact = FakeExecutor(
                exact_factory=lambda: {
                    "status": "complete",
                    "valid": True,
                    "probabilities": [float("nan"), 0.1, 0.1, 0.1],
                    "logical_mask": 0,
                    "log_joint_syndrome": -8.0,
                    "elapsed_seconds": 0.1,
                    "peak_rss_bytes": 512,
                }
            )
            decision = adaptive_runner.run_record(
                root,
                sampling,
                task_id=0,
                record_id=record["record_id"],
                executor=malformed_exact,
            )
            self.assertEqual(decision["final_status"], "unknown")
            self.assertEqual(decision["unknown_reason"], "invalid_exact_reference")
            exact_checkpoint = json.loads(
                (Path(decision["record_directory"]) / "exact.json").read_text()
            )
            self.assertEqual(exact_checkpoint["status"], "invalid")
            self.assertFalse(exact_checkpoint["valid"])

    def test_tied_posterior_is_eligible_but_never_converges(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling, _, execution = make_bundle(root)
            record = self.first_record(execution)

            def tied(chi):
                return tn_result(
                    chi,
                    probabilities=[0.4, 0.4, 0.1, 0.1],
                    log_joint=-8.0,
                )

            executor = FakeExecutor(tn_factory=tied)
            decision = adaptive_runner.run_record(
                root, sampling, task_id=0, record_id=record["record_id"], executor=executor
            )
            self.assertEqual([call[1] for call in executor.tn_calls], [32, 64, 128, 240, 480, 960])
            self.assertEqual(decision["final_status"], "unknown")
            self.assertEqual(decision["unknown_reason"], "max_chi_without_convergence")

    def test_l5_forced_anchor_runs_full_ladder_without_rewriting_early_acceptance(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling, _, execution = make_bundle(root, "l5_feasibility")
            task = next(
                task
                for task in execution["tasks"]
                if any(record["forced_anchor"] for record in task["records"])
            )
            record = next(record for record in task["records"] if record["forced_anchor"])
            executor = FakeExecutor()
            decision = adaptive_runner.run_record(
                root,
                sampling,
                task_id=task["task_id"],
                record_id=record["record_id"],
                executor=executor,
            )
            self.assertEqual(
                [call[1] for call in executor.tn_calls], [32, 64, 128, 240, 480]
            )
            self.assertEqual(executor.exact_calls, [])
            self.assertEqual(decision["final_status"], "accepted")
            adaptive = decision["adaptive_decision"]
            self.assertEqual(adaptive["accepted_chi"], 128)
            self.assertTrue(adaptive["anchor_completed"])
            self.assertEqual(adaptive["anchor_final_chi"], 480)
            self.assertFalse(adaptive["anchor_class_changed"])

    def test_dry_run_is_read_only_and_julia_command_targets_existing_entry_adapter(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling, _, execution = make_bundle(root)
            record = self.first_record(execution)
            before = sorted(str(path.relative_to(root)) for path in root.rglob("*"))
            executor = FakeExecutor()
            plan = adaptive_runner.run_record(
                root,
                sampling,
                task_id=0,
                record_id=record["record_id"],
                executor=executor,
                dry_run=True,
            )
            after = sorted(str(path.relative_to(root)) for path in root.rglob("*"))
            self.assertEqual(before, after)
            self.assertEqual(plan["status"], "dry_run")
            self.assertEqual(plan["chi_ladder"], [32, 64, 128, 240, 480, 960])
            self.assertEqual(executor.tn_calls, [])
            self.assertEqual(executor.exact_calls, [])
            command = adaptive_runner.JuliaExecutor(
                repository_root=Path(__file__).resolve().parents[3], julia="julia"
            ).command("tn", Path("record.toml"), Path("result.toml"), L=3,
                      p=0.006, chi=32, cutoff=1e-18)
            self.assertIn(str(HERE / "julia_adapter.jl"), command)
            self.assertEqual(command[0:3], ["julia", "--project=" + str(Path(__file__).resolve().parents[3]), "--startup-file=no"])

            completed = subprocess.run(
                [
                    sys.executable,
                    str(HERE / "adaptive_runner.py"),
                    "--study-root",
                    str(root),
                    "--sampling-manifest",
                    str(root / "sampling-manifest.json"),
                    "--task-id",
                    "0",
                    "--record-id",
                    record["record_id"],
                    "--dry-run",
                ],
                text=True,
                capture_output=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertEqual(json.loads(completed.stdout)["status"], "dry_run")

    def test_provenance_mismatch_refuses_before_writing_any_scientific_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling, _, execution = make_bundle(root)
            record = self.first_record(execution)
            repository = Path(__file__).resolve().parents[3]
            with self.assertRaises(adaptive_runner.IdentityMismatchError):
                adaptive_runner.run_record(
                    root,
                    sampling,
                    task_id=0,
                    record_id=record["record_id"],
                    executor=adaptive_runner.JuliaExecutor(
                        repository_root=repository, julia="julia"
                    ),
                )
            self.assertFalse((root / "outputs").exists())

    def test_exact_audit_posterior_and_log_gates_are_independent_and_inclusive(self):
        adaptive = {
            "status": "accepted",
            "probabilities": [0.7, 0.1, 0.1, 0.1],
            "logical_mask": 0,
            "log_joint_syndrome": -10.0,
        }
        policy = {
            "ml_tie_tolerance": 1e-9,
            "posterior_linf_gate": 1e-5,
            "log_joint_gate": 1e-5,
        }
        boundary = {
            "status": "complete",
            "valid": True,
            "probabilities": [0.7 - 1e-5, 0.1 + 1e-5, 0.1, 0.1],
            "logical_mask": 0,
            "log_joint_syndrome": -10.0 + 1e-5,
        }
        audit = adaptive_runner._exact_audit(adaptive, boundary, policy)
        self.assertLessEqual(audit["posterior_linf"], 1e-5)
        self.assertLessEqual(audit["log_joint_error"], 1e-5)
        self.assertTrue(audit["passed"])

        posterior_fail = json.loads(json.dumps(boundary))
        posterior_fail["probabilities"] = [
            0.7 - 1.0001e-5,
            0.1 + 1.0001e-5,
            0.1,
            0.1,
        ]
        posterior_fail["log_joint_syndrome"] = -10.0
        audit = adaptive_runner._exact_audit(adaptive, posterior_fail, policy)
        self.assertFalse(audit["posterior_gate_passed"])
        self.assertTrue(audit["log_joint_gate_passed"])
        self.assertTrue(audit["false_accept"])

        log_fail = json.loads(json.dumps(boundary))
        log_fail["probabilities"] = list(adaptive["probabilities"])
        log_fail["log_joint_syndrome"] = -10.0 + 1.0001e-5
        audit = adaptive_runner._exact_audit(adaptive, log_fail, policy)
        self.assertTrue(audit["posterior_gate_passed"])
        self.assertFalse(audit["log_joint_gate_passed"])
        self.assertTrue(audit["false_accept"])


if __name__ == "__main__":
    unittest.main()
