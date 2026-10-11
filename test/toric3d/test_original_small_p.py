import importlib.util
import math
from pathlib import Path
import unittest


SOURCE = Path(__file__).parents[2] / "examples" / "analyze_toric3d_original_small_p.py"


class OriginalSmallPAnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location("original_small_p", SOURCE)
        cls.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.module)

    def test_linear_target_crossing_uses_bracketing_points(self):
        points = [(1e-5, 0.004), (3e-5, 0.016), (1e-4, 0.05)]
        self.assertAlmostEqual(
            self.module.target_crossing(points, 0.01), 2e-5)

    def test_target_crossing_rejects_unbracketed_target(self):
        with self.assertRaises(ValueError):
            self.module.target_crossing([(1e-5, 0.001), (3e-5, 0.003)], 0.01)

    def test_power_law_recovers_inverse_square(self):
        values = [(3, 0.02 / 3**2), (5, 0.02 / 5**2),
                  (7, 0.02 / 7**2), (9, 0.02 / 9**2)]
        fit = self.module.fit_power_law(values)
        self.assertAlmostEqual(fit["exponent"], 2.0, places=12)
        self.assertAlmostEqual(fit["coefficient"], 0.02, places=12)


if __name__ == "__main__":
    unittest.main()
