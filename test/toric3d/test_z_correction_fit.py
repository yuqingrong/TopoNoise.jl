"""Independent numerical checks for the finite-size correction analysis."""

import importlib.util
from pathlib import Path
import unittest

import numpy as np


SCRIPT = Path(__file__).resolve().parents[2] / "examples" / "fit_toric3d_z_corrections.py"
if SCRIPT.exists():
    SPEC = importlib.util.spec_from_file_location("z_correction_fit", SCRIPT)
    FIT = importlib.util.module_from_spec(SPEC)
    SPEC.loader.exec_module(FIT)
else:
    FIT = None


def synthetic_rows():
    rng = np.random.default_rng(7021)
    sizes = np.array([9, 13, 17, 21, 25])
    reference = np.exp(np.mean(np.log(sizes)))
    rows = []
    for size in sizes:
        for p in np.linspace(0.0045, 0.0065, 7):
            x = (p - 0.0052) / 0.002 * (size / reference) ** (1 / 1.05)
            eta = -5 + 5 * x - 1.2 * x * x + 0.4 * reference / size
            probability = 1 / (1 + np.exp(-eta))
            rows.append(dict(steps=int(size), Ly=int(size), Lz=int(size), p=float(p),
                             shots=500000, failures=int(rng.binomial(500000, probability))))
    return rows


class CorrectionFitTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(FIT, "the finite-size correction fitter is not implemented")

    def test_likelihood_gradient_matches_independent_finite_differences(self):
        rows = synthetic_rows()
        model = FIT.CorrectionModel(rows, p_min=0.0045, p_max=0.0065, min_size=9, omega=1)
        parameters = np.array([0.43, np.log(1.2), -4.8, 4.5, -0.7, 0.6])
        loss, gradient = model.objective(parameters)
        self.assertTrue(np.isfinite(loss))
        for index in range(len(parameters)):
            shift = np.zeros(len(parameters))
            shift[index] = 1e-6
            numerical = (model.objective(parameters + shift)[0] - model.objective(parameters - shift)[0]) / 2e-6
            self.assertAlmostEqual(gradient[index], numerical, places=8)

    def test_seeded_binomial_data_recovers_threshold_and_correction(self):
        model = FIT.CorrectionModel(synthetic_rows(), p_min=0.0045, p_max=0.0065, min_size=9, omega=1)
        fit = model.fit()
        self.assertTrue(fit["accepted"], fit)
        self.assertLess(abs(fit["pc"] - 0.0052), 0.0001)
        self.assertLess(abs(fit["nu"] - 1.05), 0.15)
        self.assertLess(abs(fit["b"] - 0.4), 0.2)
        self.assertGreater(fit["deviance_p_value"], 0.001)

    def test_parameter_bound_is_not_an_accepted_fit(self):
        model = FIT.CorrectionModel(synthetic_rows(), p_min=0.0045, p_max=0.0065,
                                     min_size=9, omega=1, nu_bounds=(0.3, 0.4))
        fit = model.fit()
        self.assertFalse(fit["accepted"])
        self.assertIn("nu", fit["bound_hits"])

    def test_bootstrap_is_seeded_and_intervals_use_only_accepted_fits(self):
        model = FIT.CorrectionModel(synthetic_rows(), p_min=0.0045, p_max=0.0065, min_size=9, omega=1)
        fit = model.fit()
        one = FIT.bootstrap(model, fit, count=4, seed=723)
        two = FIT.bootstrap(model, fit, count=4, seed=723)
        self.assertEqual(one, two)
        self.assertEqual(one["accepted_count"], 4)
        self.assertEqual(len(one["samples"]), 4)


if __name__ == "__main__":
    unittest.main()
