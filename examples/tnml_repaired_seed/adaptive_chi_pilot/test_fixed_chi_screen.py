#!/usr/bin/env python3
"""Bounded saved-file tests for the paired fixed-chi screening report."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

try:
    from . import adaptive_runner, inspect_fixed_chi_screen as screen
except ImportError:
    import adaptive_runner
    import inspect_fixed_chi_screen as screen

HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parents[2]


def write_toml(path, document):
    flat = {key: value for key, value in document.items() if not isinstance(value, dict)}
    adaptive_runner._write_record_toml(path, flat)
    with path.open("a", encoding="utf-8") as stream:
        for key, table in document.items():
            if isinstance(table, dict):
                stream.write(f"\n[{key}]\n")
                for name, value in sorted(table.items()):
                    stream.write(json.dumps(name) + " = " + adaptive_runner._toml_value(value) + "\n")


def source_inventory():
    names = {str(path.relative_to(REPOSITORY)) for path in (REPOSITORY / "src").rglob("*.jl")}
    names.update({"Project.toml", "examples/tnml_repaired_seed/model.jl",
                  "examples/tnml_repaired_seed/run.jl", "examples/tnml_repaired_seed/dense_frontier.jl"})
    names.update("examples/tnml_repaired_seed/adaptive_chi_pilot/" + name for name in (
        "persistent_worker.jl", "julia_adapter.jl", "fixed_chi_screen.jl",
        "inspect_fixed_chi_screen.py", "core.py", "adaptive_runner.py", "inspect_adaptive_ladder.py",
    ))
    if (REPOSITORY / "Manifest.toml").exists():
        names.add("Manifest.toml")
    return {name: hashlib.sha256((REPOSITORY / name).read_bytes()).hexdigest() for name in names}


def fixture(directory, shots=3):
    run = {
        "schema_version": 1, "kind": "standalone_fixed_chi_screen", "L": 3,
        "p": 0.006, "p_key": "0.0060", "p_index": 0, "shots": shots,
        "chis": [64, 128], "seed_base": 202610700000, "cutoff": 1e-18,
        "startup_seconds": 7.0, "production_gate": False,
        "startup_timing_scope": "source_snapshot_backend_load_context_after_driver_import",
        "timing_policy": "first_history_cold_excluded_from_warm_medians_all_histories_scored",
        "runtime_identity": {"julia_version": "1.12.4", "blas_vendor": "lbt", "cpu_model": "arm",
                             "julia_threads": 1, "blas_threads": 4},
        "source_sha256": source_inventory(),
    }
    write_toml(directory / "run.toml", run)
    write_toml(directory / "completion.toml", {
        "status": "complete", "requested_shots": shots, "completed_shots": shots,
        "elapsed_seconds": 100.0, "source_identity_verified": True, "production_gate": False,
    })
    for index in range(1, shots + 1):
        folder = directory / "shots" / f"shot_{index:04d}"
        folder.mkdir(parents=True)
        record_id = f"{index:064x}"
        write_toml(folder / "record.toml", {
            "record_id": record_id, "shot_index": index, "true_mask": 1 if index == 2 else 0,
        })
        write_toml(folder / "metadata.toml", {
            "shot_index": index, "seed": 202610700000 + index - 1, "record_id": record_id,
            "record_file_sha256": hashlib.sha256((folder / "record.toml").read_bytes()).hexdigest(),
            "record_authenticated": True, "cold_start": index == 1,
            "chi_order": [64, 128] if index % 2 else [128, 64],
        })
        exact = {
            "record_id": record_id, "shot_index": index, "status": "complete", "valid": True,
            "probabilities": [0.7, 0.1, 0.1, 0.1], "logical_mask": 0,
            "log_joint_syndrome": -2.0, "elapsed_seconds": 0.2, "cold_start": index == 1,
        }
        write_toml(folder / "exact.toml", exact)
        for chi in (64, 128):
            decode = {64: [20.0, 2.0, 4.0], 128: [30.0, 3.0, 5.0]}[chi][index - 1]
            level = {**exact, "chi": chi, "decode_seconds": decode,
                     "decoder_build_seconds": 0.5 if index == 1 else 0.0,
                     "elapsed_seconds": decode + 0.5, "decoder_cache_hit": index != 1}
            if chi == 64 and index == 1:
                level["probabilities"] = [0.69, 0.11, 0.1, 0.1]
            if chi == 64 and index == 2:
                level.update(status="invalid", valid=False)
                for key in ("probabilities", "logical_mask", "log_joint_syndrome"):
                    level.pop(key)
            if chi == 64 and index == 3:
                level.update(probabilities=[0.1, 0.7, 0.1, 0.1], logical_mask=1)
            write_toml(folder / f"chi_{chi:04d}.toml", level)


class FixedChiScreenTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        fixture(self.directory)

    def change(self, name, update):
        path = self.directory / name
        document = adaptive_runner._read_toml(path)
        update(document)
        write_toml(path, document)

    def test_invalid_is_unknown_and_class_accuracy_is_separate_from_full_audit(self):
        report = screen.inspect_directory(self.directory)
        self.assertIn("by_chi", report)
        chi64 = report["by_chi"]["64"]
        self.assertEqual((chi64["requested_shots"], chi64["valid_shots"], chi64["unknown_shots"]), (3, 2, 1))
        self.assertEqual(chi64["class_disagreements"], 1)
        self.assertEqual(chi64["strict_class_disagreements"], 1)
        self.assertEqual(chi64["full_audit_failures"], 2)
        self.assertEqual(chi64["logical_failures"], 1)
        self.assertEqual(chi64["failure_interval"], [1 / 3, 2 / 3])
        self.assertEqual(report["exact_reference"]["logical_failures"], 1)
        chi128 = report["by_chi"]["128"]
        self.assertEqual((chi128["class_disagreements"], chi128["full_audit_failures"]), (0, 0))
        self.assertEqual(chi128["logical_failures"], 1)
        self.assertFalse(report["production_gate"])

    def test_pair_counts_and_warm_timing_keep_unknown_and_cold_histories_visible(self):
        report = screen.inspect_directory(self.directory)
        self.assertIn("paired", report)
        self.assertEqual(report["paired"]["both_valid_shots"], 2)
        self.assertEqual(report["paired"]["unknown_shots"], 1)
        self.assertEqual(report["paired"]["decoder_class_disagreements"], 1)
        chi64 = report["by_chi"]["64"]
        self.assertEqual(chi64["warm_median_decode_seconds"], 3.0)
        self.assertEqual(chi64["warm_timed_shots"], 2)
        self.assertEqual(chi64["total_decode_seconds"], 26.0)
        self.assertEqual(chi64["total_elapsed_seconds"], 27.5)
        self.assertEqual(chi64["first_history_timing"]["decode_seconds"], 20.0)
        self.assertEqual((chi64["decoder_cache_hits"], chi64["decoder_cache_misses"]), (2, 1))

    def test_tied_exact_maximizers_are_not_tie_aware_class_disagreements(self):
        for filename, mask in (("exact.toml", 0), ("chi_0064.toml", 1)):
            self.change("shots/shot_0001/" + filename,
                        lambda doc, mask=mask: doc.update(probabilities=[0.4, 0.4, 0.1, 0.1], logical_mask=mask))
        report = screen.inspect_directory(self.directory)
        self.assertIn("by_chi", report)
        chi64 = report["by_chi"]["64"]
        self.assertEqual(chi64["class_disagreements"], 1)
        self.assertEqual(chi64["strict_class_disagreements"], 2)

    def test_missing_completion_and_missing_checkpoint_are_rejected(self):
        (self.directory / "completion.toml").unlink()
        with self.assertRaises((ValueError, FileNotFoundError)):
            screen.inspect_directory(self.directory)
        self.assertFalse((self.directory / "report.json").exists())

    def test_missing_checkpoint_after_completion_is_rejected(self):
        (self.directory / "shots/shot_0003/chi_0128.toml").unlink()
        with self.assertRaises((ValueError, FileNotFoundError)):
            screen.inspect_directory(self.directory)

    def test_changed_record_hash_and_result_identity_are_rejected(self):
        with (self.directory / "shots/shot_0001/record.toml").open("a") as stream:
            stream.write("changed = true\n")
        with self.assertRaisesRegex(ValueError, "record.*(hash|content)"):
            screen.inspect_directory(self.directory)

    def test_wrong_output_identity_is_rejected(self):
        self.change("shots/shot_0002/chi_0064.toml", lambda doc: doc.update(record_id="f" * 64))
        with self.assertRaisesRegex(ValueError, "record identity"):
            screen.inspect_directory(self.directory)

    def test_changed_source_and_missing_source_inventory_are_rejected(self):
        self.change("run.toml", lambda doc: doc["source_sha256"].update({"Project.toml": "0" * 64}))
        with self.assertRaisesRegex(ValueError, "source.*changed"):
            screen.inspect_directory(self.directory)
        self.change("run.toml", lambda doc: doc["source_sha256"].pop("Project.toml"))
        with self.assertRaisesRegex(ValueError, "source.*inventory"):
            screen.inspect_directory(self.directory)

    def test_malformed_scientific_result_cannot_be_reclassified_as_invalid(self):
        self.change("shots/shot_0002/chi_0128.toml", lambda doc: doc.update(probabilities=[0.1, 0.1, 0.1, 0.1]))
        with self.assertRaisesRegex(ValueError, "malformed complete"):
            screen.inspect_directory(self.directory)

    def test_unexpected_worker_exception_and_invalid_with_error_fail(self):
        self.change("shots/shot_0002/chi_0064.toml", lambda doc: doc.update(status="exception", error="backend failure"))
        with self.assertRaisesRegex(ValueError, "exception"):
            screen.inspect_directory(self.directory)
        self.change("shots/shot_0002/chi_0064.toml", lambda doc: doc.update(status="invalid"))
        with self.assertRaisesRegex(ValueError, "invalid.*error"):
            screen.inspect_directory(self.directory)

    def test_bad_numeric_timing_and_seed_metadata_are_rejected(self):
        self.change("shots/shot_0003/chi_0128.toml", lambda doc: doc.update(decode_seconds=-1.0))
        with self.assertRaisesRegex(ValueError, "decode_seconds"):
            screen.inspect_directory(self.directory)
        self.change("shots/shot_0003/chi_0128.toml", lambda doc: doc.update(decode_seconds=5.0))
        self.change("shots/shot_0003/metadata.toml", lambda doc: doc.update(seed=2))
        with self.assertRaisesRegex(ValueError, "seed"):
            screen.inspect_directory(self.directory)

    def test_one_history_has_no_warm_median_but_is_still_scored(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            fixture(directory, shots=1)
            report = screen.inspect_directory(directory)
            for result in report["by_chi"].values():
                self.assertEqual(result["requested_shots"], 1)
                self.assertEqual(result["valid_shots"], 1)
                self.assertEqual(result["warm_timed_shots"], 0)
                self.assertIsNone(result["warm_median_decode_seconds"])
                self.assertGreater(result["total_decode_seconds"], 0.0)

    def test_nonstring_error_is_not_a_numerical_invalidity(self):
        for error in (False, 0, []):
            with self.subTest(error=error):
                self.change("shots/shot_0002/chi_0064.toml", lambda doc: doc.update(error=error))
                with self.assertRaisesRegex(ValueError, "invalid.*error"):
                    screen.inspect_directory(self.directory)

    def test_cli_retains_mismatch_counts_exits_zero_and_never_overwrites(self):
        command = [sys.executable, str(HERE / "inspect_fixed_chi_screen.py"), str(self.directory)]
        first = subprocess.run(command, capture_output=True, text=True, timeout=10)
        self.assertEqual(first.returncode, 0, first.stderr)
        self.assertTrue((self.directory / "report.json").is_file())
        compact = json.loads(first.stdout.splitlines()[0])
        self.assertEqual(compact["by_chi"]["64"]["class_disagreements"], 1)
        self.assertIn("SCREEN p=0.0060 chi=64 requested=3 valid=2 unknown=1", first.stdout)
        saved = (self.directory / "report.json").read_bytes()
        second = subprocess.run(command, capture_output=True, text=True, timeout=10)
        self.assertNotEqual(second.returncode, 0)
        self.assertEqual((self.directory / "report.json").read_bytes(), saved)


if __name__ == "__main__":
    unittest.main()
