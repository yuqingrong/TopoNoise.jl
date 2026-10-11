"""Statistical checks using independent seeded binomial fixtures.

Run: .CondaPkg/.pixi/envs/default/bin/python -m unittest discover \
    -s test/toric3d -p 'test_z_threshold_fit.py'
"""

import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np


SCRIPT = Path(__file__).resolve().parents[2] / "examples" / "fit_toric3d_z_threshold.py"
if SCRIPT.exists():
    SPEC = importlib.util.spec_from_file_location("z_threshold_fit", SCRIPT)
    FIT = importlib.util.module_from_spec(SPEC)
    SPEC.loader.exec_module(FIT)
else:
    FIT = None


def synthetic_rows(seed=941):
    rng = np.random.default_rng(seed)
    rows = []
    for size in (5, 7, 9, 11):
        for probability in np.linspace(0.014, 0.036, 9):
            x = ((probability - 0.025) / 0.022) * (size / 7) ** (1 / 1.2)
            rate = 1 / (1 + np.exp(-(-1 + 2 * x + 0.25 * x * x)))
            failures = int(rng.binomial(8000, rate))
            rows.append(dict(steps=size, Ly=size, Lz=size, p=float(probability),
                             shots=8000, failures=failures,
                             logical_1_failures=failures, logical_2_failures=0,
                             seed=seed, decoder="synthetic", elapsed_seconds=0))
    return rows


class ThresholdFitTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(FIT, "the finite-size threshold fitter is not implemented")

    def test_wilson_retains_uncertainty_at_zero_failures(self):
        low, high = FIT.wilson_interval(0, 1000)
        self.assertAlmostEqual(low, 0.0, places=14)
        self.assertGreater(high, 0.003)
        self.assertLess(high, 0.004)
        low, high = FIT.wilson_interval(1000, 1000)
        self.assertGreater(low, 0.996)
        self.assertAlmostEqual(high, 1.0)

    def test_seeded_binomial_crossing_recovers_known_threshold(self):
        result = FIT.analyze(synthetic_rows(), p_min=0.014, p_max=0.036,
                             min_size=5, bootstrap=0, seed=1234, sensitivity=False)
        self.assertEqual(result["status"], "ok", result)
        self.assertLess(abs(result["fit"]["pc"] - 0.025), 0.0015)
        self.assertLess(abs(result["fit"]["nu"] - 1.2), 0.5)
        self.assertGreaterEqual(len(result["crossings"]), 3)
        self.assertGreater(result["fit"]["degrees_of_freedom"], 0)
        self.assertGreater(result["fit"]["deviance_p_value"], 0.001)

    def test_non_crossing_curves_are_not_given_a_threshold(self):
        rows = synthetic_rows()
        for row in rows:
            row["failures"] = int(row["shots"] * (0.2 + row["p"] - 0.005 * row["steps"]))
        result = FIT.analyze(rows, p_min=0.014, p_max=0.036, min_size=5,
                             bootstrap=3, seed=1234, sensitivity=False)
        self.assertEqual(result["status"], "no_crossing")
        self.assertIsNone(result["fit"])
        self.assertIsNone(result["pc_interval"])

    def test_all_zero_failures_are_not_a_crossing(self):
        rows = synthetic_rows()
        for row in rows:
            row["failures"] = 0
        result = FIT.analyze(rows, p_min=0.014, p_max=0.036, min_size=5,
                             bootstrap=0, seed=1234, sensitivity=False)
        self.assertEqual(result["status"], "no_crossing")
        self.assertIsNone(result["fit"])

    def test_seeded_bootstrap_reproduces_intervals(self):
        arguments = dict(p_min=0.014, p_max=0.036, min_size=5,
                         bootstrap=6, seed=52, sensitivity=False)
        one = FIT.analyze(synthetic_rows(), **arguments)
        two = FIT.analyze(synthetic_rows(), **arguments)
        self.assertEqual(one["bootstrap_samples"], two["bootstrap_samples"])
        self.assertEqual(one["pc_interval"], two["pc_interval"])
        self.assertEqual(one["crossings"], two["crossings"])
        self.assertEqual(len(one["bootstrap_samples"]), 6)

    def test_bound_hitting_fit_is_reported_as_inconclusive(self):
        result = FIT.analyze(synthetic_rows(), p_min=0.014, p_max=0.036,
                             bootstrap=0, nu_bounds=(0.3, 0.4), sensitivity=False)
        self.assertEqual(result["status"], "inconclusive")
        self.assertIn("nu", result["fit"]["bound_hits"])
        self.assertTrue(any("bound" in item for item in result["diagnostic"]))

    def test_invalid_counts_or_non_cubic_geometry_are_rejected(self):
        rows = synthetic_rows()
        rows[0]["failures"] = 8001
        with self.assertRaises(ValueError):
            FIT.analyze(rows, p_min=0.014, p_max=0.036)
        rows = synthetic_rows()
        rows[0]["Ly"] = 4
        with self.assertRaises(ValueError):
            FIT.analyze(rows, p_min=0.014, p_max=0.036)

    def test_cli_saves_data_uncertainties_figures_and_metadata(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            source = folder / "raw.csv"
            rows = synthetic_rows()
            with source.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            (folder / "scan-metadata.toml").write_text('noise = "test fixture"\n')
            output = folder / "analysis"
            code = FIT.main([str(source), "--p-min", "0.014", "--p-max", "0.036",
                             "--bootstrap", "2", "--output-dir", str(output)])
            self.assertEqual(code, 0)
            report = json.loads((output / "threshold-fit.json").read_text())
            self.assertEqual(report["observable"], "any residual transverse logical Z")
            self.assertTrue(report["sensitivity"])
            self.assertTrue((output / "metadata" / "scan-metadata.toml").is_file())
            for extension in ("pdf", "svg", "png"):
                self.assertGreater((output / f"threshold-fit.{extension}").stat().st_size, 1000)
            self.assertTrue((output / "threshold-bootstrap.csv").is_file())
            self.assertTrue((output / "threshold-points.csv").is_file())


if __name__ == "__main__":
    unittest.main()
