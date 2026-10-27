"""Reporting integrity checks using explicitly synthetic temporary fixtures."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location("repaired_plot", Path(__file__).with_name("plot.py"))
PLOT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PLOT)


def fixture(root, name="partial", *, complete=False, phase="production", p=.01, seed=9):
    """Artificial counts solely for testing; no decoder is executed."""
    path = root / phase / name
    (path / "shots").mkdir(parents=True)
    run = "synthetic-" + name
    (path / "metadata.toml").write_text(
        f'L=3\np={p}\nshots=4\nchi=240\nrng_seed={seed}\nrun_id="{run}"\ntn_uses_both_records=true\n')
    records = [f'run_id="{run}"\nrequested_shots=4\n']
    for i in range(1, 5):
        truth = i % 2
        records.append(f'[[records]]\nshot_index={i}\nrecord_id="record-{i}"\ntrue_mask={truth}\nmwpm_failure={str(i == 2).lower()}\n')
        if i == 4 and not complete:
            continue
        valid = i != 3 or complete
        # The invalid third output deliberately has a stale true TNfailure.
        text = (f'run_id="{run}"\nrecord_id="record-{i}"\nshot_index={i}\n'
                f'L=3\np={p}\nchi=240\ntrue_mask={truth}\nmwpm_failure={str(i == 2).lower()}\n'
                f'status="{"complete" if valid else "invalid"}"\nvalid={str(valid).lower()}\n'
                f'logical_mask=0\nTNfailure={str(bool(truth)).lower()}\nelapsed_seconds=1.0\n')
        (path / "shots" / f"shot_{i:04d}.toml").write_text(text)
    (path / "input-corpus.toml").write_text("\n".join(records))
    return path


class ReportIntegrity(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="tnml-plot-synthetic-")
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def test_invalid_and_absent_outputs_remain_unknown(self):
        run = PLOT.load_run(fixture(self.root), "production")
        self.assertEqual((run["valid_shots"], run["tn_failures"], run["unknown_shots"]), (2, 1, 2))
        self.assertEqual((run["tn_lower_bound"], run["tn_upper_bound"]), (.25, .75))
        self.assertIsNone(run["tn_rate"])
        self.assertIsNone(run["tn_wilson95"])
        self.assertEqual(run["mwpm_rate"], .25)

    def test_complete_run_gets_unconditional_binomial_rate(self):
        run = PLOT.load_run(fixture(self.root, complete=True), "production")
        self.assertEqual((run["tn_rate"], run["unknown_shots"]), (.5, 0))
        self.assertEqual(run["tn_wilson95"], PLOT.wilson(2, 4))

    def test_identity_mismatch_is_rejected(self):
        path = fixture(self.root)
        shot = path / "shots/shot_0001.toml"
        shot.write_text(shot.read_text().replace('record_id="record-1"', 'record_id="wrong"'))
        with self.assertRaisesRegex(ValueError, "identity mismatch"):
            PLOT.load_run(path, "production")

    def test_repeated_seed_not_counted_as_independent(self):
        runs = [PLOT.load_run(fixture(self.root, name), "production") for name in ("one", "two")]
        with self.assertRaisesRegex(ValueError, "repeated RNG seed"):
            PLOT.aggregate_runs(runs)

    def test_independent_shards_combine_full_denominators(self):
        runs = [PLOT.load_run(fixture(self.root, "one"), "production"),
                PLOT.load_run(fixture(self.root, "two", complete=True, seed=10), "production")]
        point, = PLOT.aggregate_runs(runs)
        self.assertEqual((point["requested_shots"], point["tn_failures"], point["unknown_shots"]), (8, 3, 2))
        self.assertEqual((point["tn_lower_bound"], point["tn_upper_bound"]), (3/8, 5/8))

    def test_renderer_excludes_diagnostics_and_exports_artifacts(self):
        fixture(self.root)
        fixture(self.root, "other", complete=True, phase="diagnostics", p=.02)
        output = self.root / "synthetic-figures"
        report = PLOT.main(["--root", str(self.root), "--output", str(output), "--synthetic-preview"])
        self.assertEqual(len(report["runs"]), 1)
        self.assertEqual(report["totals"], {"requested_shots":4,"valid_shots":2,"unknown_shots":2})
        for name in ("tnml-matched-figure.png", "tnml-matched-figure.pdf", "summary.csv", "report.json"):
            self.assertGreater((output / name).stat().st_size, 0)
        saved = json.loads((output / "report.json").read_text())
        self.assertIsNone(saved["points"][0]["tn_rate"])


if __name__ == "__main__":
    unittest.main()
