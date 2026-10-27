"""Integrity-validator regressions using preserved, real runner artifacts."""
import copy
import csv
import importlib.util
import io
import json
from pathlib import Path
import shutil
import tempfile
import tomllib
import unittest


HERE = Path(__file__).resolve().parent
STUDY = HERE.parent


class ProductionAuditTests(unittest.TestCase):
    def checker(self):
        path = HERE / "check_production.py"
        self.assertTrue(path.is_file(), "independent production checker is missing")
        spec = importlib.util.spec_from_file_location("production_checker", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def fixture(self):
        directory = STUDY / "benchmark-x"
        return tomllib.loads((directory / "checkpoint.toml").read_text())

    def test_missing_shard_and_duplicate_cell_are_rejected(self):
        module = self.checker()
        selection = json.loads((STUDY / "selection.json").read_text())
        jobs = tomllib.loads((STUDY / "jobs.toml").read_text())["jobs"]
        self.assertEqual(module.validate_selection(selection, jobs),
                         {"x_only": 106000, "z_only": 98000})
        with self.assertRaises(module.AuditError):
            module.validate_selection(selection, jobs[:-1])
        broken = copy.deepcopy(selection)
        broken["jobs"][0]["rates"].append(0.04)
        with self.assertRaises(module.AuditError):
            module.validate_selection(broken, broken["jobs"])

    def test_checkpoint_content_tampering_is_rejected(self):
        module = self.checker()
        raw = (STUDY / "benchmark-x/checkpoint.toml").read_bytes()
        checkpoint = module.parse_checkpoint(raw, "fixture")
        self.assertEqual(checkpoint["states"][0]["completed_shots"], 20)
        with self.assertRaises(module.AuditError):
            module.parse_checkpoint(raw.replace(b"completed_shots = 20", b"completed_shots = 19"), "tampered")

    def test_logical_union_and_paired_histogram_are_enforced(self):
        module = self.checker()
        state = self.fixture()["states"][0]
        module.validate_state(state, 20, "fixture")
        broken = copy.deepcopy(state)
        broken["logical_1_counts"][0] = 20
        with self.assertRaises(module.AuditError):
            module.validate_state(broken, 20, "bad logical union")
        broken = copy.deepcopy(state)
        broken["joint_counts"][0] += 1
        with self.assertRaises(module.AuditError):
            module.validate_state(broken, 20, "bad histogram")

    def test_false_completion_and_raw_count_disagreement_are_rejected(self):
        module = self.checker()
        checkpoint = self.fixture()
        rows = list(csv.DictReader(io.StringIO((STUDY / "benchmark-x/raw.csv").read_text())))
        module.validate_raw(rows, checkpoint, "fixture")
        broken = copy.deepcopy(rows)
        broken[0]["full_prefix_failures"] = "5"
        with self.assertRaises(module.AuditError):
            module.validate_raw(broken, checkpoint, "wrong raw count")
        broken = copy.deepcopy(rows)
        broken[0]["complete"] = "false"
        with self.assertRaises(module.AuditError):
            module.validate_raw(broken, checkpoint, "incomplete point")

    def test_seed_reuse_against_excluded_runs_is_rejected(self):
        module = self.checker()
        module.validate_seed_sets(["0000000000000001", "0000000000000002"], {"0000000000000003"})
        for new, old in [(["0000000000000001"] * 2, set()),
                         (["0000000000000001"], {"0000000000000001"})]:
            with self.assertRaises(module.AuditError):
                module.validate_seed_sets(new, old)

    def test_missing_immutable_snapshot_is_rejected(self):
        module = self.checker()
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = STUDY / "benchmark-x"
            shutil.copyfile(source / "checkpoint.toml", directory / "checkpoint.toml")
            shutil.copytree(source / "checkpoints", directory / "checkpoints")
            self.assertEqual(module.validate_archive(directory, self.fixture()), 4)
            next((directory / "checkpoints").glob("checkpoint-000002-*.toml")).unlink()
            with self.assertRaises(module.AuditError):
                module.validate_archive(directory, self.fixture())


if __name__ == "__main__":
    unittest.main()
