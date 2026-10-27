"""Small, deterministic tests for full-grid orchestration; no Julia jobs run."""
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch


DRIVER = Path(__file__).with_name("run_full_grid.py")


class FullGridTests(unittest.TestCase):
    def setUp(self):
        self.assertTrue(DRIVER.is_file(), "the full-grid driver must exist")
        spec = importlib.util.spec_from_file_location("full_grid_driver", DRIVER)
        self.driver = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.driver)
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        self.selection = self.driver.initialize_selection(self.root)

    def gate(self, L=3, **changes):
        evidence = self.root / "validation" / "evidence.json"
        evidence.parent.mkdir(exist_ok=True)
        evidence.write_text('{"note":"synthetic evidence for orchestration unit tests"}')
        gate = dict(status="passed", chi=240, cutoff=1e-18,
                    negative_tolerance=1e-10, evidence=["validation/evidence.json"],
                    exact_reference_for_every_shot=L == 3)
        gate.update(changes)
        path = self.root / "validation" / "gates.json"
        document = json.loads(path.read_text()) if path.exists() else {"schema_version": 1, "gates": {}}
        document["gates"][str(L)] = gate
        path.write_text(json.dumps(document))

    def test_full_grid_keeps_reference_rates_and_zero_is_not_1000_shots(self):
        selection = self.selection
        self.assertEqual(selection["sizes"], [3, 5, 7, 9, 11])
        self.assertEqual(len(selection["p_values"]), 23)
        self.assertEqual(selection["p_values"][0], 0.0)
        self.assertEqual(selection["p_values"][-1], 0.0425)
        self.assertIn(0.0025, selection["p_values"])
        self.assertIn(0.01, selection["p_values"])
        self.assertEqual(selection["shots_per_nonzero_point"], 1000)
        zero = self.driver.selected_points(selection, "zero", [3, 5], None)
        self.assertEqual([point["shots"] for point in zero], [1, 1])
        self.assertEqual([point["seed"] for point in zero], [202610280003, 202610280005])
        self.assertTrue(zero[0]["exact"])
        self.assertFalse(zero[1]["exact"])
        production = self.driver.selected_points(selection, "production", None, None)
        self.assertEqual(len(production), 110)
        self.assertEqual(len({point["seed"] for point in production}), 110)
        self.assertTrue(all(point["shots"] == 1000 and point["p"] > 0 for point in production))
        self.assertTrue(all(point["seed"] not in range(2026102711, 2026102715) for point in production))

    def test_selection_is_frozen_and_rejects_changes(self):
        path = self.root / "scan-selection.json"
        before = path.read_bytes()
        self.driver.initialize_selection(self.root)
        self.assertEqual(path.read_bytes(), before)
        modified = json.loads(before)
        modified["shots_per_nonzero_point"] = 20
        path.write_text(json.dumps(modified))
        with self.assertRaisesRegex(ValueError, "frozen"):
            self.driver.initialize_selection(self.root)

    def test_subsets_cannot_add_rates_or_sizes_or_duplicate_points(self):
        for sizes, rates in [([4], None), ([3], [0]), ([3], [.003]), ([3, 3], None), ([3], [.01, .01])]:
            with self.subTest(sizes=sizes, rates=rates), self.assertRaises(ValueError):
                self.driver.selected_points(self.selection, "production", sizes, rates)

    def test_unvalidated_sizes_and_mismatched_settings_are_refused(self):
        point = self.driver.selected_points(self.selection, "production", [5], [.01])[0]
        with self.assertRaisesRegex(ValueError, "validation"):
            self.driver.require_gate(self.root, self.selection, point)
        self.gate(5, status="blocked")
        with self.assertRaisesRegex(ValueError, "blocked"):
            self.driver.require_gate(self.root, self.selection, point)
        self.gate(5, chi=120)
        with self.assertRaisesRegex(ValueError, "chi"):
            self.driver.require_gate(self.root, self.selection, point)
        self.gate(5, evidence=["absent.json"])
        with self.assertRaisesRegex(ValueError, "evidence"):
            self.driver.require_gate(self.root, self.selection, point)
        self.gate(5)
        self.driver.require_gate(self.root, self.selection, point)

    def test_l3_gate_requires_exact_per_shot(self):
        self.gate(3, exact_reference_for_every_shot=False)
        point = self.driver.selected_points(self.selection, "production", [3], [.01])[0]
        with self.assertRaisesRegex(ValueError, "exact"):
            self.driver.require_gate(self.root, self.selection, point)

    def test_command_preserves_shots_limits_and_exact_policy(self):
        point = self.driver.selected_points(self.selection, "production", [3], [.01])[0]
        command = self.driver.build_command(self.root, self.selection, point, "/test/julia")
        self.assertEqual(command[command.index("--decode-seconds") + 1], "600")
        self.assertEqual(command[command.index("--total-seconds") + 1], "86400")
        self.assertEqual(command[-8:-2], ["3", "0.01", "1000", "240", str(point["seed"]), str(self.root / point["output"])])
        self.assertEqual(command[-2], "exact=true")
        zero = self.driver.selected_points(self.selection, "zero", [5], None)[0]
        command = self.driver.build_command(self.root, self.selection, zero, "/test/julia")
        self.assertEqual(command[command.index("--decode-seconds") + 1], "300")
        self.assertEqual(command[command.index("--total-seconds") + 1], "600")
        self.assertEqual(command[-2], "exact=false")

    def artifacts(self, point, *, bad=None):
        """Write tiny realistic runner outputs for audit tests."""
        output = self.root / point["output"]
        output.mkdir(parents=True, exist_ok=True)
        (output / "shots").mkdir(exist_ok=True)
        (output / "exact").mkdir(exist_ok=True)
        shots = point["shots"]
        (output / "metadata.toml").write_text(
            f'L = {point["L"]}\np = {point["p"]}\nshots = {shots}\nchi = 240\n'
            f'rng_seed = {point["seed"]}\nexact = {str(point["exact"]).lower()}\n'
            f'cutoff = {point["cutoff"]}\nnegative_tolerance = 1e-10\nrun_id = "test-run"\n')
        hashes = self.selection["scientific_source_sha256"]
        with (output / "metadata.toml").open("a") as stream:
            stream.write(f'project_sha256 = "{hashes["Project.toml"]}"\n')
            if "Manifest.toml" in hashes:
                stream.write(f'manifest_sha256 = "{hashes["Manifest.toml"]}"\n')
            stream.write("[source_sha256]\n")
            for name, digest in hashes.items():
                if name not in ("Project.toml", "Manifest.toml"):
                    stream.write(f'"{name}" = "{digest}"\n')
        (output / "summary.toml").write_text(
            f'L = {point["L"]}\np = {point["p"]}\nrequested_shots = {shots}\n'
            f'saved_histories = {shots}\nsaved_outputs = {shots}\nvalid_shots = {shots}\n'
            'unknown_shots = 0\ninvalid_shots = 0\nexception_shots = 0\n'
            f'tn_failures = {1 if bad == "disagreement" else 0}\nrun_id = "test-run"\n')
        for i in range(1, shots + 1):
            text = (f'L = {point["L"]}\np = {point["p"]}\nchi = 240\ncutoff = {point["cutoff"]}\nshot_index = {i}\n'
                    f'record_id = "record-{i}"\nrun_id = "test-run"\n'
                    'status = "complete"\nvalid = true\nlogical_mask = 0\ntrue_mask = 0\n'
                    'TNfailure = false\nprobabilities = [1.0,0.0,0.0,0.0]\n')
            if bad == "invalid" and i == 1:
                text = text.replace('valid = true', 'valid = false')
            if point["exact"]:
                exact = text.replace('TNfailure = false\n', '')
                (output / "exact" / f"shot_{i:04d}.toml").write_text(exact)
                text += '[exact_comparison]\nml_tie_aware_match = true\nml_strict_match = true\nposterior_linf = 0.0\n'
                if bad == "disagreement" and i == 1:
                    text = (text.replace('logical_mask = 0', 'logical_mask = 1')
                            .replace('TNfailure = false', 'TNfailure = true')
                            .replace('probabilities = [1.0,0.0,0.0,0.0]', 'probabilities = [0.0,1.0,0.0,0.0]')
                            .replace('ml_tie_aware_match = true', 'ml_tie_aware_match = false'))
            (output / "shots" / f"shot_{i:04d}.toml").write_text(text)
        (output / "input-corpus.toml").write_text('run_id = "test-run"\n' + ''.join(
            f'[[records]]\nshot_index = {i}\nrecord_id = "record-{i}"\nraw_fault_indices = []\n' for i in range(1, shots + 1)))
        if bad == "missing":
            (output / "shots" / "shot_0001.toml").unlink()
        if bad == "missing_exact":
            (output / "exact" / "shot_0001.toml").unlink()
        return output

    def small_point(self):
        point = self.driver.selected_points(self.selection, "production", [3], [.01])[0]
        point["shots"] = 2
        return point

    def test_complete_valid_shots_have_independent_exact_audit(self):
        point = self.small_point()
        self.artifacts(point)
        audit = self.driver.audit_point(self.root, point)
        self.assertTrue(audit["passed"])
        self.assertEqual(audit["exact_checked"], 2)
        self.assertEqual(audit["valid_shots"], 2)
        self.assertEqual(audit["unknown_shots"], 0)

    def test_audit_does_not_accept_invalid_missing_or_exact_disagreement(self):
        for bad in ("invalid", "missing", "missing_exact", "disagreement"):
            with self.subTest(bad=bad):
                point = self.small_point()
                self.artifacts(point, bad=bad)
                audit = self.driver.audit_point(self.root, point)
                self.assertFalse(audit["passed"])
                if bad == "disagreement":
                    self.assertEqual(audit["exact_disagreements"], 1)

    def test_exact_decision_is_recomputed_and_ties_are_accepted(self):
        point = self.small_point()
        output = self.artifacts(point, bad="disagreement")
        shot = output / "shots" / "shot_0001.toml"
        shot.write_text(shot.read_text().replace('ml_tie_aware_match = false', 'ml_tie_aware_match = true'))
        self.assertEqual(self.driver.audit_point(self.root, point)["exact_disagreements"], 1)
        exact = output / "exact" / "shot_0001.toml"
        exact.write_text(exact.read_text().replace('probabilities = [1.0,0.0,0.0,0.0]', 'probabilities = [0.5,0.5,0.0,0.0]'))
        audit = self.driver.audit_point(self.root, point)
        self.assertTrue(audit["passed"])
        self.assertEqual(audit["strict_matches"], 1)

    def test_zero_check_rejects_nonzero_fault_history(self):
        point = self.driver.selected_points(self.selection, "zero", [5], None)[0]
        output = self.artifacts(point)
        self.assertTrue(self.driver.audit_point(self.root, point)["passed"])
        corpus = output / "input-corpus.toml"
        corpus.write_text(corpus.read_text().replace("raw_fault_indices = []", "raw_fault_indices = [1]"))
        self.assertFalse(self.driver.audit_point(self.root, point)["passed"])

    def fake_worker(self, command, **kwargs):
        L, p, shots, chi, seed, output, exact, cutoff = command[-8:]
        output = Path(output)
        point = dict(L=int(L), p=float(p), shots=int(shots), chi=int(chi), seed=int(seed),
                     exact=exact == "exact=true", cutoff=float(cutoff.split("=",1)[1]), output=str(output.relative_to(self.root)),
                     phase="zero" if float(p) == 0 else "production")
        self.artifacts(point)
        (output / "supervision.json").write_text(json.dumps({"status": "completed", "exit_code": 0}))
        return SimpleNamespace(returncode=0)

    def test_completed_resume_calls_runner_for_identity_check(self):
        self.gate()
        calls = []
        def worker(command, **kwargs):
            calls.append(command)
            return self.fake_worker(command, **kwargs)
        for _ in range(2):
            result = self.driver.execute(self.root, self.selection, "zero", [3], None, "julia", run=worker)
            self.assertEqual(result, 0)
        self.assertEqual(len(calls), 2)
        status = json.loads((self.root / "scan-status.json").read_text())
        point = next(iter(status["points"].values()))
        self.assertEqual(point["status"], "completed")
        self.assertEqual(len(point["attempts"]), 2)

    def test_timeout_blocks_size_without_automatic_retry(self):
        calls = []
        def timeout(command, **kwargs):
            calls.append(command)
            output = Path(command[-3])
            output.mkdir(parents=True, exist_ok=True)
            (output / "supervision.json").write_text(json.dumps({"status": "contraction_time_limit"}))
            return SimpleNamespace(returncode=2)
        result = self.driver.execute(self.root, self.selection, "zero", [3], None, "julia", run=timeout)
        self.assertEqual(result, 2)
        result = self.driver.execute(self.root, self.selection, "zero", [3], None, "julia", run=timeout)
        self.assertEqual(result, 2)
        self.assertEqual(len(calls), 1)
        result = self.driver.execute(self.root, self.selection, "zero", [3], None, "julia", run=self.fake_worker, retry_blocked=True)
        self.assertEqual(result, 0)

    def test_existing_manual_timeout_is_not_automatically_retried(self):
        point = self.driver.selected_points(self.selection, "zero", [5], None)[0]
        output = self.root / point["output"]
        output.mkdir(parents=True)
        (output / "supervision.json").write_text('{"status":"contraction_time_limit"}')
        result = self.driver.execute(self.root, self.selection, "zero", [5], None, "julia",
                                     run=lambda *a, **k: self.fail("timeout must not be retried"))
        self.assertEqual(result, 2)

    def test_stop_file_stops_before_dispatch_and_preserves_pending_points(self):
        (self.root / "STOP").write_text("stop after current point")
        result = self.driver.execute(self.root, self.selection, "zero", [3, 5], None, "julia",
                                     run=lambda *a, **k: self.fail("STOP must prevent dispatch"))
        self.assertEqual(result, 3)
        status = json.loads((self.root / "scan-status.json").read_text())
        self.assertEqual(status["status"], "stopped")
        self.assertTrue(all(p["status"] == "pending" for p in status["points"].values()))

    def test_failed_audit_halts_remaining_points_of_size(self):
        self.gate()
        calls = []
        def worker(command, **kwargs):
            calls.append(command)
            result = self.fake_worker(command, **kwargs)
            first = Path(command[-3]) / "shots" / "shot_0001.toml"
            first.write_text(first.read_text().replace('valid = true', 'valid = false', 1))
            return result
        result = self.driver.execute(self.root, self.selection, "production", [3], [.0025, .01], "julia", run=worker)
        self.assertEqual(result, 2)
        self.assertEqual(len(calls), 1)
        status = json.loads((self.root / "scan-status.json").read_text())
        self.assertEqual([row["status"] for row in status["points"].values()], ["blocked", "pending"])
        self.assertEqual(status["points"][next(iter(status["points"]))]["completed_counts"]["unknown_shots"], 1)

    def test_unvalidated_subset_is_rejected_before_any_worker(self):
        self.gate(3)
        with self.assertRaisesRegex(ValueError, "validation"):
            self.driver.execute(self.root, self.selection, "production", [3, 5], [.01], "julia",
                                run=lambda *a, **k: self.fail("no worker before all requested gates pass"))

    def test_source_changes_and_status_identity_changes_are_refused(self):
        with patch.object(self.driver, "scientific_source_hashes", return_value={"src/changed.jl": "different"}):
            with self.assertRaisesRegex(ValueError, "frozen"):
                self.driver.initialize_selection(self.root)
        (self.root / "scan-status.json").write_text('{"selection_sha256":"different"}')
        with self.assertRaisesRegex(ValueError, "identity"):
            self.driver.execute(self.root, self.selection, "zero", [3], None, "julia",
                                run=lambda *a, **k: self.fail("identity mismatch must prevent dispatch"))

    def test_single_worker_lock_is_enforced(self):
        with self.driver.scan_lock(self.root):
            with self.assertRaisesRegex(ValueError, "lock"):
                self.driver.execute(self.root, self.selection, "zero", [3], None, "julia",
                                    run=lambda *a, **k: self.fail("concurrent worker must be refused"))

    def test_source_change_during_long_scan_stops_before_next_point(self):
        calls = []
        def worker(command, **kwargs):
            calls.append(command)
            return self.fake_worker(command, **kwargs)
        original = self.selection["scientific_source_sha256"]
        with patch.object(self.driver, "scientific_source_hashes", side_effect=[original, {"changed": "source"}]):
            result = self.driver.execute(self.root, self.selection, "zero", [3, 5], None, "julia", run=worker)
        self.assertEqual(result, 2)
        self.assertEqual(len(calls), 1)
        status = json.loads((self.root / "scan-status.json").read_text())
        self.assertEqual(status["status"], "source_changed")

    def test_dispatch_audit_rejects_saved_project_hash_different_from_frozen(self):
        def worker(command, **kwargs):
            result = self.fake_worker(command, **kwargs)
            metadata = Path(command[-3]) / "metadata.toml"
            original = self.selection["scientific_source_sha256"]["Project.toml"]
            metadata.write_text(metadata.read_text().replace(original, "different-project"))
            return result
        result = self.driver.execute(self.root, self.selection, "zero", [3], None, "julia", run=worker)
        self.assertEqual(result, 2)

    def test_metadata_identity_binds_all_scientific_hashes(self):
        point = self.small_point()
        output = self.artifacts(point)
        metadata = self.driver.read_toml(output / "metadata.toml")
        helper = getattr(self.driver, "validate_metadata_identity", None)
        self.assertTrue(callable(helper), "saved metadata needs a frozen-source identity validator")
        helper(metadata, point, self.selection)
        for field in ("source_sha256", "project_sha256", "manifest_sha256"):
            if field == "manifest_sha256" and "Manifest.toml" not in self.selection["scientific_source_sha256"]:
                continue
            altered = dict(metadata)
            altered[field] = {"src/other.jl": "different"} if field == "source_sha256" else "different"
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, "hash"):
                helper(altered, point, self.selection)


    def test_precision_selection_freezes_cutoff_and_fresh_seed_offset(self):
        root = self.root / "precision"
        selection = self.driver.initialize_selection(root, cutoff=1e-26, rng_seed_offset=100000)
        self.assertEqual(selection["cutoff"], 1e-26)
        self.assertEqual(selection["rng_seed_offset"], 100000)
        before = (root / "scan-selection.json").read_bytes()
        self.assertEqual(self.driver.initialize_selection(root), selection)
        self.assertEqual((root / "scan-selection.json").read_bytes(), before)
        for options in ({"cutoff": 1e-18}, {"rng_seed_offset": 0}):
            with self.subTest(options=options), self.assertRaisesRegex(ValueError, "frozen"):
                self.driver.initialize_selection(root, **options)
        old = self.driver.selected_points(self.selection, "production", [3], [.0025])[0]
        new = self.driver.selected_points(selection, "production", [3], [.0025])[0]
        self.assertEqual(new["seed"], old["seed"] + 100000)
        self.assertEqual(new["cutoff"], 1e-26)
        self.assertEqual(self.driver.build_command(root, selection, new, "julia")[-1], "cutoff=1e-26")
        old_zero = self.driver.selected_points(self.selection, "zero", [3])[0]
        new_zero = self.driver.selected_points(selection, "zero", [3])[0]
        self.assertEqual(new_zero["seed"], old_zero["seed"] + 100000)

    def test_legacy_selection_without_seed_offset_is_not_rewritten(self):
        path = self.root / "scan-selection.json"
        legacy = json.loads(path.read_text())
        legacy.pop("rng_seed_offset")
        path.write_text(json.dumps(legacy))
        before = path.read_bytes()
        self.assertEqual(self.driver.initialize_selection(self.root), legacy)
        self.assertEqual(path.read_bytes(), before)

    def test_precision_audit_accepts_selected_cutoff_and_rejects_copied_checkpoints(self):
        self.selection["cutoff"] = 1e-26
        point = self.small_point()
        output = self.artifacts(point)
        self.assertTrue(self.driver.audit_point(self.root, point, self.selection)["passed"])
        for component in ("shots", "exact"):
            path = output / component / "shot_0001.toml"
            original = path.read_text()
            path.write_text(original.replace("cutoff = 1e-26", "cutoff = 1e-18"))
            self.assertFalse(self.driver.audit_point(self.root, point, self.selection)["passed"])
            path.write_text(original)

    def test_invalid_precision_configuration_is_rejected(self):
        for cutoff in (float("nan"), float("inf"), -1e-26, 1.0):
            with self.subTest(cutoff=cutoff), self.assertRaises(ValueError):
                self.driver.selection_template(cutoff=cutoff)
        for offset in (-1, True, 0.5, 2**63):
            with self.subTest(offset=offset), self.assertRaises(ValueError):
                self.driver.selection_template(rng_seed_offset=offset)


if __name__ == "__main__":
    unittest.main()
