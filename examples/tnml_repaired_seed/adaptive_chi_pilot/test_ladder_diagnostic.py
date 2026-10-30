"""Bounded checks of full-ladder reporting against the production policy."""

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

try:
    from . import inspect_adaptive_ladder as diagnostic
except ImportError:
    import inspect_adaptive_ladder as diagnostic


def documents():
    run = {
        "record_id": "a" * 64, "L": 3, "p": 0.006, "p_key": "0.0060",
        "seed": 202610600000, "chi_ladder": [32, 64, 128, 240, 480, 960],
        "cutoff": 1e-18, "startup_seconds": 7.0,
        "record_authenticated": True,
    }
    exact = {
        "record_id": run["record_id"], "status": "complete", "valid": True,
        "probabilities": [0.7, 0.1, 0.1, 0.1], "logical_mask": 0,
        "log_joint_syndrome": -2.0, "elapsed_seconds": 0.2,
    }
    levels = [
        {**exact, "chi": chi, "elapsed_seconds": 2.0,
         "decoder_build_seconds": 0.1, "decode_seconds": 1.9}
        for chi in run["chi_ladder"]
    ]
    return run, exact, levels


def save_fixture(directory):
    run, exact, levels = documents()
    record_path = directory / "record.toml"
    record_path.write_text(f'record_id = "{run["record_id"]}"\n')
    run["record_file_sha256"] = hashlib.sha256(record_path.read_bytes()).hexdigest()
    diagnostic.adaptive_runner._write_record_toml(directory / "run.toml", run)
    here = Path(diagnostic.__file__).parent
    with (directory / "run.toml").open("a") as stream:
        stream.write('\n[runtime_identity]\njulia_threads = 1\nblas_threads = 4\n')
        stream.write('\n[policy_source_sha256]\n')
        for name in ("core.py", "adaptive_runner.py", "inspect_adaptive_ladder.py"):
            digest = hashlib.sha256((here / name).read_bytes()).hexdigest()
            stream.write(f'"{name}" = "{digest}"\n')
    diagnostic.adaptive_runner._write_record_toml(directory / "exact.toml", exact)
    for level in levels:
        diagnostic.adaptive_runner._write_record_toml(directory / f'chi_{level["chi"]:04d}.toml', level)


class LadderDiagnosticTests(unittest.TestCase):
    def test_full_ladder_cost_is_separate_from_policy_prefix_cost(self):
        summary = diagnostic.summarize_documents(*documents())
        self.assertEqual(summary["adaptive_decision"]["accepted_chi"], 128)
        self.assertEqual(summary["full_ladder_tn_seconds"], 12.0)
        self.assertEqual(summary["measured_prefix_tn_seconds"], 6.0)
        self.assertEqual(summary["startup_seconds"], 7.0)
        self.assertTrue(summary["accepted_exact_audit"]["passed"])
        self.assertTrue(summary["highest_chi_exact_audit"]["passed"])

    def test_exact_class_match_does_not_hide_posterior_gate_failure(self):
        run, exact, levels = documents()
        for level in levels:
            level["probabilities"] = [0.7001, 0.0999, 0.1, 0.1]
        summary = diagnostic.summarize_documents(run, exact, levels)
        audit = summary["accepted_exact_audit"]
        self.assertTrue(audit["ml_tie_aware_match"])
        self.assertFalse(audit["posterior_gate_passed"])
        self.assertTrue(audit["false_accept"])

    def test_mismatched_record_and_incomplete_ladder_are_rejected(self):
        run, exact, levels = documents()
        changed = deepcopy(levels)
        changed[1]["record_id"] = "b" * 64
        with self.assertRaisesRegex(ValueError, "record identity"):
            diagnostic.summarize_documents(run, exact, changed)
        with self.assertRaisesRegex(ValueError, "full ladder"):
            diagnostic.summarize_documents(run, exact, levels[:-1])

    def test_invalid_level_is_retained_but_never_reported_as_accepted(self):
        run, exact, levels = documents()
        levels[0] = {
            "record_id": run["record_id"], "chi": 32, "status": "invalid",
            "valid": False, "elapsed_seconds": 2.0,
            "decoder_build_seconds": 0.1, "decode_seconds": 1.9,
        }
        summary = diagnostic.summarize_documents(run, exact, levels)
        self.assertEqual(summary["invalid_chis"], [32])
        self.assertNotEqual(summary["prefix_decisions"][0]["decision"]["status"], "accepted")
        self.assertEqual(summary["full_ladder_tn_seconds"], 12.0)
        self.assertTrue(summary["highest_chi_exact_audit"]["passed"])

    def test_exception_cannot_be_conflated_with_low_chi_invalidity(self):
        run, exact, levels = documents()
        levels[0].update(status="exception", valid=False)
        with self.assertRaisesRegex(ValueError, "exception"):
            diagnostic.summarize_documents(run, exact, levels)

    def test_changed_saved_history_is_rejected(self):
        with tempfile.TemporaryDirectory() as root:
            directory = Path(root)
            save_fixture(directory)
            self.assertEqual(diagnostic.inspect_directory(directory)["record_id"], "a" * 64)
            with (directory / "record.toml").open("a") as stream:
                stream.write('changed = true\n')
            with self.assertRaisesRegex(ValueError, "record content changed"):
                diagnostic.inspect_directory(directory)

    def test_changed_policy_is_rejected(self):
        with tempfile.TemporaryDirectory() as root:
            directory = Path(root)
            save_fixture(directory)
            path = directory / "run.toml"
            document = path.read_text()
            here = Path(diagnostic.__file__).parent
            digest = hashlib.sha256((here / "core.py").read_bytes()).hexdigest()
            path.write_text(document.replace(digest, "0" * 64))
            with self.assertRaisesRegex(ValueError, "policy source changed"):
                diagnostic.inspect_directory(directory)

    def test_cli_saves_summary_and_refuses_to_overwrite(self):
        with tempfile.TemporaryDirectory() as root:
            directory = Path(root)
            save_fixture(directory)
            command = [sys.executable, diagnostic.__file__, str(directory)]
            first = subprocess.run(command, capture_output=True, text=True, timeout=10)
            self.assertEqual(first.returncode, 0, first.stderr)
            summary = json.loads(first.stdout)
            self.assertFalse(summary["production_gate"])
            saved = (directory / "summary.json").read_bytes()
            second = subprocess.run(command, capture_output=True, text=True, timeout=10)
            self.assertNotEqual(second.returncode, 0)
            self.assertEqual((directory / "summary.json").read_bytes(), saved)


if __name__ == "__main__":
    unittest.main()
