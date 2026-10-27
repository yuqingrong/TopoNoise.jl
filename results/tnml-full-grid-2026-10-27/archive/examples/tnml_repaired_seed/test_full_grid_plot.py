"""Full-grid reporting tests; artificial data remain in temporary directories."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


def module(name, filename):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(filename))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


PLOT = module("full_grid_plot", "plot_full_grid.py")
OLD_FIXTURE = module("old_plot_fixtures", "test_plot.py").fixture
DRIVER = module("grid_driver_for_plot_tests", "run_full_grid.py")


def FIXTURE(root, *, complete=False, phase="production", p=.01):
    selection_file = root / "scan-selection.json"
    if not selection_file.exists() or phase == "diagnostics":
        return OLD_FIXTURE(root, complete=complete, phase=phase, p=p)
    selection = json.loads(selection_file.read_text())
    if p not in selection["p_values"]:
        return OLD_FIXTURE(root, complete=complete, phase=phase, p=p)
    point, = DRIVER.selected_points(selection, phase, [3], [p])
    if phase == "production":
        gate_path = root / "validation/gates.json"
        gates = json.loads(gate_path.read_text())
        gates["gates"]["3"]["exact_reference_for_every_shot"] = point["exact"]
        gate_path.write_text(json.dumps(gates))
    N = point["shots"]
    path = root / point["output"]
    (path / "shots").mkdir(parents=True, exist_ok=True)
    (path / "exact").mkdir(exist_ok=True)
    run_id = "synthetic-" + point["key"]
    hashes = selection["scientific_source_sha256"]
    (path / "metadata.toml").write_text(
        f'L=3\np={p}\nshots={N}\nchi=240\nrng_seed={point["seed"]}\nrun_id="{run_id}"\n'
        f'tn_uses_both_records=true\nexact={str(point["exact"]).lower()}\ncutoff=1e-18\nnegative_tolerance=1e-10\n'
        + ''.join(f'{field}="{hashes[filename]}"\n' for field, filename in (("project_sha256", "Project.toml"), ("manifest_sha256", "Manifest.toml")) if filename in hashes)
        + '[source_sha256]\n' + ''.join(f'{json.dumps(k)}="{v}"\n' for k, v in hashes.items() if k not in ("Project.toml", "Manifest.toml")))
    corpus = [f'run_id="{run_id}"\nrequested_shots={N}\n']
    valid_count = failures = output_count = 0
    for i in range(1, N+1):
        truth = i % 2 if p else 0
        mwpm = i == 2 and p > 0
        identity = f'run_id="{run_id}"\nrecord_id="record-{i}"\nshot_index={i}\n'
        corpus.append(f'[[records]]\n{identity}true_mask={truth}\nmwpm_failure={str(mwpm).lower()}\nraw_fault_indices=[]\n')
        if not complete and i > 3:
            continue
        valid = complete or i != 3
        text = (identity + f'L=3\np={p}\nchi=240\ntrue_mask={truth}\nmwpm_failure={str(mwpm).lower()}\n'
                f'status="{"complete" if valid else "invalid"}"\nvalid={str(valid).lower()}\n'
                f'logical_mask=0\nTNfailure={str(bool(truth)).lower()}\nelapsed_seconds=1.0\n'
                'probabilities=[1.0,0.0,0.0,0.0]\n')
        if point["exact"] and valid:
            (path / "exact" / f"shot_{i:04d}.toml").write_text(text)
            text += '[exact_comparison]\nml_tie_aware_match=true\nposterior_linf=0.0\n'
        (path / "shots" / f"shot_{i:04d}.toml").write_text(text)
        output_count += 1
        valid_count += valid
        failures += truth if valid else 0
    (path / "input-corpus.toml").write_text(''.join(corpus))
    (path / "summary.toml").write_text(
        f'run_id="{run_id}"\nL=3\np={p}\nrequested_shots={N}\nvalid_shots={valid_count}\n'
        f'tn_failures={failures}\nunknown_shots={N-valid_count}\nsaved_histories={N}\nsaved_outputs={output_count}\n')
    return path


def zero_fixture(root, probabilities=None):
    path = FIXTURE(root, phase="zero", complete=True, p=0)
    for file in [path / "input-corpus.toml", *sorted((path / "shots").glob("*.toml"))]:
        text = file.read_text().replace("true_mask=1", "true_mask=0")
        text = text.replace("mwpm_failure=true", "mwpm_failure=false").replace("TNfailure=true", "TNfailure=false")
        if file.name != "input-corpus.toml":
            text = text.replace("probabilities=[1.0,0.0,0.0,0.0]", "probabilities=" + json.dumps(probabilities or [1., 0., 0., 0.]))
        file.write_text(text)
    return path


class FullGridIntegrity(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="tnml-full-grid-synthetic-")
        self.root = Path(self.temporary.name)
        self.selection = DRIVER.selection_template()
        self.selection.update(sizes=[3, 5], p_values=[0, .006, .01], shots_per_nonzero_point=1000,
                              exact_reference_sizes=[], zero_checks_per_size=1)
        (self.root / "scan-selection.json").write_text(json.dumps(self.selection))
        (self.root / "validation").mkdir()
        (self.root / "validation/evidence.json").write_text('{"purpose":"synthetic test fixture only"}')
        (self.root / "validation/gates.json").write_text(json.dumps(dict(schema_version=1, gates={
            "3": dict(status="passed", chi=240, cutoff=1e-18, negative_tolerance=1e-10,
                      exact_reference_for_every_shot=False, evidence=["validation/evidence.json"]),
            "5": dict(status="blocked", reason="runtime gate")
        })))

    def tearDown(self):
        self.temporary.cleanup()

    def test_analytical_zero_has_no_statistical_shots_or_interval(self):
        report = PLOT.build_report(self.root)
        zeros = [point for point in report["grid"] if point["p"] == 0]
        self.assertEqual(len(zeros), 2)
        for point in zeros:
            self.assertEqual((point["status"], point["analytical_rate"]), ("analytical", 0.0))
            self.assertEqual(point["requested_shots"], 0)
            self.assertIsNone(point["tn_rate"])
            self.assertIsNone(point["tn_wilson95"])
        self.assertEqual(report["totals"]["planned_production_shots"], 4000)
        self.assertEqual(report["totals"]["valid_production_shots"], 0)

    def test_pending_and_blocked_points_do_not_get_failure_estimates(self):
        report = PLOT.build_report(self.root)
        positives = [point for point in report["grid"] if point["p"] > 0]
        self.assertEqual([p["status"] for p in positives], ["pending", "pending", "blocked", "blocked"])
        for point in positives:
            self.assertEqual((point["unknown_shots"], point["tn_lower_bound"], point["tn_upper_bound"]), (1000, 0, 1))
            self.assertIsNone(point["tn_rate"])

    def test_numerical_zero_checks_are_never_statistical_production_shots(self):
        zero_fixture(self.root)
        report = PLOT.build_report(self.root)
        self.assertTrue(report["zero_checks"][0]["regression_passed"])
        self.assertIsNone(report["zero_checks"][0]["tn_wilson95"])
        self.assertIsNone(report["zero_checks"][0]["mwpm_wilson95"])
        self.assertEqual(report["totals"]["allocated_production_shots"], 0)
        self.assertEqual(report["totals"]["valid_production_shots"], 0)
        self.assertTrue(all(point["requested_shots"] == 0 for point in report["grid"] if point["p"] == 0))

    def test_zero_regression_requires_correct_posterior_not_just_label(self):
        zero_fixture(self.root, [.5, .5, 0., 0.])
        check, = PLOT.build_report(self.root)["zero_checks"]
        self.assertFalse(check["regression_passed"])

    def test_partial_run_counts_all_unfinished_requested_histories(self):
        FIXTURE(self.root, p=.01)
        report = PLOT.build_report(self.root)
        point = next(p for p in report["grid"] if p["L"] == 3 and p["p"] == .01)
        self.assertEqual((point["valid_shots"], point["tn_failures"], point["unknown_shots"]), (2, 1, 998))
        self.assertEqual((point["tn_lower_bound"], point["tn_upper_bound"]), (.001, .999))
        self.assertIsNone(point["tn_rate"])
        self.assertIsNone(point["tn_wilson95"])

    def test_previous_pilot_is_never_pooled_into_new_target(self):
        previous = self.root / "old-study"
        FIXTURE(previous, complete=True, p=.006)
        report = PLOT.build_report(self.root, previous)
        self.assertEqual(report["pilot"]["totals"]["requested_shots"], 4)
        point = next(p for p in report["grid"] if p["L"] == 3 and p["p"] == .006)
        self.assertEqual((point["valid_shots"], point["unknown_shots"]), (0, 1000))
        self.assertEqual(report["totals"]["valid_production_shots"], 0)

    def test_complete_target_has_unconditional_wilson_interval(self):
        self.selection["shots_per_nonzero_point"] = 4
        (self.root / "scan-selection.json").write_text(json.dumps(self.selection))
        FIXTURE(self.root, complete=True, p=.01)
        report = PLOT.build_report(self.root)
        point = next(p for p in report["grid"] if p["L"] == 3 and p["p"] == .01)
        self.assertEqual((point["status"], point["tn_rate"], point["unknown_shots"]), ("complete", .5, 0))
        self.assertEqual(point["tn_wilson95"], PLOT.BASE.wilson(2, 4))

    def test_exact_reference_disagreement_blocks_completed_point(self):
        self.selection.update(shots_per_nonzero_point=4, exact_reference_sizes=[3])
        (self.root / "scan-selection.json").write_text(json.dumps(self.selection))
        path = FIXTURE(self.root, complete=True, p=.01)
        shot = path / "shots/shot_0003.toml"
        shot.write_text(shot.read_text().replace("ml_tie_aware_match=true", "ml_tie_aware_match=false"))
        report = PLOT.build_report(self.root)
        point = next(p for p in report["grid"] if p["L"] == 3 and p["p"] == .01)
        self.assertEqual(point["status"], "blocked")
        self.assertIn("exact", point["reason"])
        self.assertIsNone(point["tn_rate"])

    def test_missing_exact_reference_does_not_get_validated_rate(self):
        self.selection.update(shots_per_nonzero_point=4, exact_reference_sizes=[3])
        (self.root / "scan-selection.json").write_text(json.dumps(self.selection))
        path = FIXTURE(self.root, complete=True, p=.01)
        (path / "exact/shot_0001.toml").unlink()
        point = next(p for p in PLOT.build_report(self.root)["grid"] if p["L"] == 3 and p["p"] == .01)
        self.assertEqual((point["status"], point["unknown_shots"]), ("blocked", 0))
        self.assertIsNone(point["tn_rate"])

    def test_copied_run_with_wrong_seed_is_rejected(self):
        path = FIXTURE(self.root, p=.01)
        metadata = path / "metadata.toml"
        metadata.write_text(metadata.read_text().replace(f'rng_seed={self.selection["production_rng_seeds"]["3"]["0.0100"]}', 'rng_seed=1'))
        with self.assertRaisesRegex(ValueError, "metadata|seed|configuration"):
            PLOT.build_report(self.root)

    def test_copied_run_with_wrong_scientific_source_hash_is_rejected(self):
        path = FIXTURE(self.root, p=.01)
        metadata = path / "metadata.toml"
        source = self.selection["scientific_source_sha256"]["examples/tnml_repaired_seed/run.jl"]
        metadata.write_text(metadata.read_text().replace(source, "wrong-source"))
        with self.assertRaisesRegex(ValueError, "source|hash"):
            PLOT.build_report(self.root)

    def test_unrequested_probability_is_rejected(self):
        FIXTURE(self.root, complete=True, p=.02)
        with self.assertRaisesRegex(ValueError, "outside"):
            PLOT.build_report(self.root)

    def test_driver_block_does_not_discard_saved_outcomes(self):
        FIXTURE(self.root, p=.01)
        (self.root / "scan-status.json").write_text(json.dumps(dict(points={"some-key": dict(
            phase="production", L=3, p=.01, status="blocked", reason="contraction limit")
        })))
        point = next(p for p in PLOT.build_report(self.root)["grid"] if p["L"] == 3 and p["p"] == .01)
        self.assertEqual((point["status"], point["valid_shots"], point["unknown_shots"]), ("blocked", 2, 998))

    def test_size_block_marks_pending_points_but_preserves_earlier_complete(self):
        self.selection["shots_per_nonzero_point"] = 4
        (self.root / "scan-selection.json").write_text(json.dumps(self.selection))
        FIXTURE(self.root, complete=True, p=.006)
        (self.root / "scan-status.json").write_text(json.dumps(dict(size_blocks={
            "production:L3": dict(point="production:L3-p0.0100-chi240", reason="contraction limit")
        })))
        points = [point for point in PLOT.build_report(self.root)["grid"] if point["L"] == 3 and point["p"] > 0]
        self.assertEqual([point["status"] for point in points], ["complete", "blocked"])
        self.assertIn("contraction limit", points[1]["reason"])

    def test_terminal_supervision_overrides_stale_running_status(self):
        path = FIXTURE(self.root, p=.01)
        (path / "supervision.json").write_text(json.dumps(dict(status="contraction_time_limit", exit_code=-15)))
        (self.root / "scan-status.json").write_text(json.dumps(dict(points={"key": dict(
            phase="production", L=3, p=.01, status="running", reason="")
        })))
        point = next(p for p in PLOT.build_report(self.root)["grid"] if p["L"] == 3 and p["p"] == .01)
        self.assertEqual(point["status"], "blocked")
        self.assertIn("contraction_time_limit", point["reason"])

    def test_diagnostics_are_not_production_and_render_is_marked_synthetic(self):
        FIXTURE(self.root, phase="diagnostics", complete=True, p=.006)
        output = self.root / "synthetic-figures"
        report = PLOT.main(["--root", str(self.root), "--output", str(output), "--no-pilot", "--synthetic-preview"])
        self.assertEqual(report["runs"], [])
        self.assertEqual(report["totals"]["valid_production_shots"], 0)
        self.assertTrue(report["synthetic_preview"])
        for filename in ("tnml-full-grid.png", "tnml-full-grid.pdf", "report.json"):
            self.assertGreater((output / filename).stat().st_size, 0)


if __name__ == "__main__":
    unittest.main()
