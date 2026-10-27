"""HPC task tests use synthetic checkpoint files; no Julia/Slurm jobs run."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch


DIRECTORY = Path(__file__).resolve().parent


def import_file(name):
    path = DIRECTORY / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"hpc_test_{name}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TaskTests(unittest.TestCase):
    def setUp(self):
        for name in ("manifest.py", "run_task.py"):
            self.assertTrue((DIRECTORY / name).is_file(), f"missing HPC implementation: {name}")
        self.hpc = import_file("manifest")
        self.worker = import_file("run_task")
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        self.manifest = self.hpc.create_manifest(self.root)

    def task(self, stage="zero", task_id=0):
        return self.hpc.task_for_id(self.manifest, stage, task_id)

    def test_full_plan_has_1000_shots_per_point_and_separate_diagnostics(self):
        self.assertEqual(len(self.hpc.stage_tasks(self.manifest, "zero")), 5)
        smoke = self.hpc.stage_tasks(self.manifest, "smoke")
        self.assertEqual([t["p"] for t in smoke], [.006, .01, .02])
        self.assertTrue(all(t["shots"] == 2 and t["exact"] for t in smoke))
        production = self.hpc.stage_tasks(self.manifest, "production")
        self.assertEqual(len(production), 440)
        self.assertEqual(sum(t["shots"] for t in production), 22000)
        self.assertEqual({t["L"] for t in production}, {3})
        rates = {t["p"] for t in production}
        self.assertEqual(len(rates), 22)
        self.assertEqual(min(rates), .0025)
        self.assertEqual(max(rates), .0425)
        for p in rates:
            tasks = [t for t in production if t["p"] == p]
            self.assertEqual(len(tasks), 20)
            self.assertEqual(sum(t["shots"] for t in tasks), 1000)
            self.assertEqual({t["shard_id"] for t in tasks}, set(range(20)))
        self.assertTrue(all(t["shots"] == 50 and t["exact"] for t in production))

    def test_diagnostic_seeds_are_paired_only_across_contraction_dimensions(self):
        pairs = {}
        for task in self.hpc.stage_tasks(self.manifest, "feasibility"):
            pairs.setdefault(task["pair_id"], []).append(task)
        self.assertEqual(len(pairs), 3)
        for tasks in pairs.values():
            self.assertEqual({t["chi"] for t in tasks}, {240, 480})
            self.assertEqual(len({t["seed"] for t in tasks}), 1)
            self.assertTrue(all(t["L"] == 5 and t["shots"] == 3 and not t["exact"] for t in tasks))
        seed_sets = {stage: {t["seed"] for t in self.hpc.stage_tasks(self.manifest, stage)}
                     for stage in self.manifest["stages"]}
        for a, seeds in seed_sets.items():
            for b, others in seed_sets.items():
                if a != b:
                    self.assertTrue(seeds.isdisjoint(others))
        self.assertEqual(len(seed_sets["production"]), 440)

    def test_manifest_freeze_is_idempotent_and_conflicting_settings_are_refused(self):
        path = self.root / "manifest.json"
        before = path.read_bytes()
        self.assertEqual(self.hpc.create_manifest(self.root), self.manifest)
        self.assertEqual(path.read_bytes(), before)
        with self.assertRaisesRegex(ValueError, "frozen"):
            self.hpc.create_manifest(self.root, cutoff=1e-18)
        altered = json.loads(before)
        altered["cutoff"] = 1e-18
        path.write_text(json.dumps(altered))
        with self.assertRaisesRegex(ValueError, "identity"):
            self.hpc.load_manifest(self.root)

    def test_duplicate_or_unsafe_task_destinations_are_rejected(self):
        for mode in ("duplicate", "absolute", "parent", "seed"):
            bad = copy.deepcopy(self.manifest)
            tasks = bad["stages"]["production"]["tasks"]
            if mode == "duplicate":
                tasks[1]["output"] = tasks[0]["output"]
            elif mode == "absolute":
                tasks[0]["output"] = "/tmp/untrusted-output"
            elif mode == "parent":
                tasks[0]["output"] = "production/../outside"
            else:
                tasks[1]["seed"] = tasks[0]["seed"]
            bad["manifest_id"] = self.hpc.manifest_identity(bad)
            with self.subTest(mode=mode), self.assertRaises(ValueError):
                self.hpc.validate_manifest(bad)

    def test_frozen_scientific_and_wrapper_sources_are_checked(self):
        for field in ("scientific_source_hashes", "wrapper_hashes"):
            with self.subTest(field=field), patch.object(self.hpc, field, return_value={"changed": "source"}):
                with self.assertRaisesRegex(ValueError, "changed"):
                    self.hpc.load_manifest(self.root)

    def test_python_configuration_and_optional_validated_path_are_frozen(self):
        self.assertIn("CondaPkg.toml", self.manifest["wrapper_sha256"])
        self.assertNotIn("CondaPkg.toml", self.manifest["scientific_source_sha256"])
        repository = self.root / "mock-repository"
        experiment = repository / "examples" / "experiment"
        hpc = experiment / "hpc"
        hpc.mkdir(parents=True)
        for name in ("run_full_grid.py", "supervise.py", "plot.py"):
            (experiment / name).write_text("# synthetic wrapper\n")
        config = repository / "CondaPkg.toml"
        config.write_text('[deps]\npython = ">=3.11,<3.13"\n')
        with patch.object(self.hpc, "REPOSITORY", repository), patch.object(self.hpc, "EXPERIMENT", experiment), patch.object(self.hpc, "HERE", hpc):
            original = self.hpc.wrapper_hashes()
            self.assertNotIn(".hpc-python-path", original)
            config.write_text('[deps]\npython = "3.12"\n')
            changed_config = self.hpc.wrapper_hashes()
            self.assertNotEqual(original["CondaPkg.toml"], changed_config["CondaPkg.toml"])
            marker = repository / ".hpc-python-path"
            marker.write_text(str(repository / "native-python") + "\n")
            with_marker = self.hpc.wrapper_hashes()
            self.assertIn(".hpc-python-path", with_marker)
            marker.write_text(str(repository / "other-python") + "\n")
            self.assertNotEqual(with_marker[".hpc-python-path"], self.hpc.wrapper_hashes()[".hpc-python-path"])

    def test_dry_run_requires_no_slurm_and_creates_no_task_outputs(self):
        result = self.worker.execute_task(self.root, "smoke", 0, julia="julia", dry_run=True,
                                          environ={}, run=lambda *a, **k: self.fail("dry run must not launch"))
        self.assertEqual(result["status"], "dry_run")
        self.assertIn("cutoff=1e-26", result["command"])
        self.assertIn("exact=true", result["command"])
        self.assertEqual(result["task"]["shots"], 2)
        self.assertFalse((self.root / self.task("smoke")["output"]).exists())

    def test_real_dispatch_requires_slurm_and_matching_array_id(self):
        with self.assertRaisesRegex(ValueError, "SLURM_JOB_ID"):
            self.worker.execute_task(self.root, "smoke", 0, environ={})
        with self.assertRaisesRegex(ValueError, "array"):
            self.worker.execute_task(self.root, "smoke", 0,
                                     environ={"SLURM_JOB_ID": "42", "SLURM_ARRAY_TASK_ID": "1"})

    def artifacts(self, task, *, bad_exact=False):
        output = self.root / task["output"]
        output.mkdir(parents=True, exist_ok=True)
        (output / "shots").mkdir(exist_ok=True)
        (output / "exact").mkdir(exist_ok=True)
        base = (f'L = {task["L"]}\np = {task["p"]}\nshots = {task["shots"]}\nchi = {task["chi"]}\n'
                f'cutoff = {task["cutoff"]}\nnegative_tolerance = 1e-10\nrng_seed = {task["seed"]}\n'
                f'exact = {str(task["exact"]).lower()}\nrun_id = "test-{task["stage"]}-{task["task_id"]}"\n')
        hashes = self.manifest["scientific_source_sha256"]
        base += f'project_sha256 = "{hashes["Project.toml"]}"\n'
        if "Manifest.toml" in hashes:
            base += f'manifest_sha256 = "{hashes["Manifest.toml"]}"\n'
        base += "[source_sha256]\n" + "".join(f'"{p}" = "{v}"\n' for p,v in hashes.items() if p not in ("Project.toml", "Manifest.toml"))
        (output / "metadata.toml").write_text(base)
        run_id = f'test-{task["stage"]}-{task["task_id"]}'
        corpus = f'run_id = "{run_id}"\n'
        for i in range(1, task["shots"] + 1):
            record_id = f'record-{task["stage"]}-{task["task_id"]}-{i}'
            corpus += (f'[[records]]\nshot_index = {i}\nrecord_id = "{record_id}"\nraw_fault_indices = []\n')
            text = (f'run_id = "{run_id}"\nrecord_id = "{record_id}"\nshot_index = {i}\n'
                    f'L = {task["L"]}\np = {task["p"]}\nchi = {task["chi"]}\ncutoff = {task["cutoff"]}\n'
                    'status = "complete"\nvalid = true\nlogical_mask = 0\ntrue_mask = 0\nTNfailure = false\n'
                    'probabilities = [1.0,0.0,0.0,0.0]\n')
            if task["exact"]:
                (output / "exact" / f"shot_{i:04d}.toml").write_text(text)
                text += '[exact_comparison]\nml_tie_aware_match = true\n'
            (output / "shots" / f"shot_{i:04d}.toml").write_text(text)
        (output / "input-corpus.toml").write_text(corpus)
        (output / "summary.toml").write_text(
            f'run_id = "{run_id}"\nL = {task["L"]}\np = {task["p"]}\nrequested_shots = {task["shots"]}\n'
            f'valid_shots = {task["shots"]}\ntn_failures = 0\nunknown_shots = 0\n')
        (output / "supervision.json").write_text('{"status":"completed","exit_code":0}')
        if bad_exact:
            (output / "exact" / "shot_0001.toml").unlink()
        return output

    def test_production_requires_every_audited_smoke_task(self):
        with self.assertRaisesRegex(ValueError, "smoke"):
            self.worker.execute_task(self.root, "production", 0, dry_run=True, environ={})
        for task in self.hpc.stage_tasks(self.manifest, "smoke"):
            self.artifacts(task)
        result = self.worker.execute_task(self.root, "production", 0, dry_run=True, environ={})
        self.assertEqual(result["status"], "dry_run")
        task = self.task("smoke", 1)
        (self.root / task["output"] / "exact" / "shot_0001.toml").unlink()
        with self.assertRaisesRegex(ValueError, "smoke"):
            self.worker.execute_task(self.root, "production", 0, dry_run=True, environ={})

    def test_metadata_and_saved_exact_files_are_bound_to_task(self):
        task = self.task("smoke")
        output = self.artifacts(task)
        audit = self.hpc.audit_task(self.root, task, self.manifest)
        self.assertTrue(audit["passed"])
        metadata = output / "metadata.toml"
        metadata.write_text(metadata.read_text().replace(f'rng_seed = {task["seed"]}', 'rng_seed = 1'))
        self.assertFalse(self.hpc.audit_task(self.root, task, self.manifest)["passed"])
        self.artifacts(task, bad_exact=True)
        self.assertFalse(self.hpc.audit_task(self.root, task, self.manifest)["passed"])

    def test_per_task_lock_refuses_duplicate_worker(self):
        output = self.root / self.task()["output"]
        with self.worker.task_lock(output):
            with self.assertRaisesRegex(ValueError, "lock"):
                self.worker.execute_task(self.root, "zero", 0, environ={"SLURM_JOB_ID": "42"},
                                         run=lambda *a, **k: self.fail("locked task must not launch"))

    def test_complete_task_records_audit_without_shared_status(self):
        task = self.task()
        seen = []
        def fake(command, **kwargs):
            seen.append((command, kwargs))
            self.artifacts(task)
            return SimpleNamespace(returncode=0)
        result = self.worker.execute_task(self.root, "zero", 0,
                                          environ={"SLURM_JOB_ID": "42"}, run=fake)
        self.assertEqual(result["status"], "completed")
        self.assertEqual(len(seen), 1)
        self.assertEqual(seen[0][1]["env"]["OPENBLAS_NUM_THREADS"], "1")
        self.assertIn("--threads=1,0", seen[0][0])
        self.assertTrue(seen[0][1]["pass_fds"])
        saved = json.loads((self.root / task["output"] / "task-status.json").read_text())
        self.assertEqual(saved["manifest_id"], self.manifest["manifest_id"])
        self.assertTrue(saved["attempts"][-1]["audit"]["passed"])
        self.assertFalse((self.root / "scan-status.json").exists())

    def test_scientific_halt_is_not_automatically_retried(self):
        task = self.task()
        output = self.root / task["output"]
        output.mkdir(parents=True)
        (output / "halt.toml").write_text('reason = "invalid_or_exception_tn_output"\n')
        with self.assertRaisesRegex(ValueError, "halt"):
            self.worker.execute_task(self.root, "zero", 0, environ={"SLURM_JOB_ID": "42"},
                                     run=lambda *a, **k: self.fail("scientific halt must not be retried"))

    def test_task_ids_and_existing_symlink_outputs_are_rejected(self):
        for stage, task_id in (("unknown", 0), ("zero", -1), ("zero", 5)):
            with self.subTest(stage=stage, task_id=task_id), self.assertRaises(ValueError):
                self.hpc.task_for_id(self.manifest, stage, task_id)
        outside = self.root.parent / (self.root.name + "-outside")
        outside.mkdir()
        self.addCleanup(outside.rmdir)
        (self.root / "zero").symlink_to(outside, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, "outside"):
            self.worker.execute_task(self.root, "zero", 0, dry_run=True, environ={})


if __name__ == "__main__":
    unittest.main()
