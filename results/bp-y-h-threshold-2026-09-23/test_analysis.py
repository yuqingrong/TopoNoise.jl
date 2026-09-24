"""Independent finite-size-scaling and crossing-selection fixtures."""
import importlib
import json
from pathlib import Path
import unittest

import numpy as np


def fixture_rows(center=0.032, shots=1_000_000, flat=False):
    rows = []
    for d in (9, 11, 13, 15):
        for p in np.arange(0, 0.10001, 0.005):
            x = (p - center) * d ** (1 / 1.45)
            rate = 0.48 if flat else float(np.clip(0.12 + 1.4 * x + 0.8 * x*x, 0, 0.5))
            k = int(round(shots * rate))
            rows.append(dict(channel="y_only", distance=d, p=round(float(p), 6),
                             shots=shots, counts=[shots-k, k, 0, 0]))
    return rows


class AnalysisTests(unittest.TestCase):
    def setUp(self):
        try:
            self.module = importlib.import_module("analysis")
        except ModuleNotFoundError:
            self.module = None
        self.assertIsNotNone(self.module, "Local scaling analysis is not implemented")

    def test_recovers_known_threshold_and_exponent(self):
        p = np.tile(np.arange(0.022, 0.04201, 0.002), 4)
        d = np.repeat([9, 11, 13, 15], 11)
        x = (p - 0.032) * d ** (1 / 1.45)
        y = 0.12 + 1.4*x + 0.8*x*x
        data = dict(p=p, d=d, n=np.full(len(p), 50000), y=y)
        fit = self.module.fit_local(data, 0.022, 0.042)
        self.assertAlmostEqual(fit["p_c"], 0.032, places=8)
        self.assertAlmostEqual(fit["nu"], 1.45, places=6)
        reverse = self.module.fit_local({k: v[::-1] for k, v in data.items()}, 0.022, 0.042)
        self.assertAlmostEqual(fit["p_c"], reverse["p_c"], places=10)

    def test_known_crossing_is_selected_but_saturation_is_not(self):
        selected = self.module.select_crossing(fixture_rows(), "y_only", "x")
        self.assertEqual(selected["status"], "selected")
        self.assertEqual(selected["center"], 0.032)
        flat = self.module.select_crossing(fixture_rows(flat=True), "y_only", "x")
        self.assertEqual(flat["status"], "unresolved")
        zeros = fixture_rows()
        for row in zeros:
            row["counts"] = [row["shots"], 0, 0, 0]
        self.assertEqual(self.module.select_crossing(zeros, "y_only", "x")["status"], "unresolved")

    def test_joint_bootstrap_preserves_perfect_correlation(self):
        rows = [dict(p=0.02, distance=9, shots=1000, counts=[700, 0, 0, 300])]
        counts = self.module.resample_joint(rows, np.random.default_rng(123))
        self.assertEqual(counts[0].sum(), 1000)
        self.assertEqual(counts[0, 1], 0)
        self.assertEqual(counts[0, 2], 0)
        self.assertEqual(counts[0, 1] + counts[0, 3], counts[0, 2] + counts[0, 3])

    def test_dispersed_pilot_crossings_trigger_refinement(self):
        path = Path(__file__).resolve().parent / "fixtures" / "dispersed_h_pilot.json"
        rows = json.loads(path.read_text())
        selected = self.module.select_crossing(rows, "hadamard", "x")
        self.assertEqual(selected["status"], "ambiguous")
        self.assertIn(0.05, selected["refine_rates"])
        phase = self.module.select_crossing(rows, "hadamard", "z")
        self.assertEqual(phase["status"], "unresolved")

    def test_unidentifiable_flat_fit_is_rejected(self):
        p = np.tile(np.arange(0.022, 0.04201, 0.002), 4)
        d = np.repeat([9, 11, 13, 15], 11)
        data = dict(p=p, d=d, n=np.full(len(p), 50000), y=np.full(len(p), 0.25))
        with self.assertRaises(ValueError):
            self.module.fit_local(data, 0.022, 0.042)


if __name__ == "__main__":
    unittest.main(verbosity=2)
