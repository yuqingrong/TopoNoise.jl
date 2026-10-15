"""Tests for the H+6 final-only Z aspect-ratio analysis.

Run with the managed Python environment:
  .CondaPkg/.pixi/envs/default/bin/python -m unittest discover \
      -s test/toric3d -p 'test_z_aspect_ratio.py'
"""

import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


SCRIPT = (Path(__file__).resolve().parents[2]
          / "examples" / "analyze_toric3d_z_aspect_ratio.py")
if SCRIPT.exists():
    SPEC = importlib.util.spec_from_file_location("z_aspect_ratio", SCRIPT)
    ANALYSIS = importlib.util.module_from_spec(SPEC)
    SPEC.loader.exec_module(ANALYSIS)
else:
    ANALYSIS = None


FIXTURE = """steps,Ly,Lz,p,shots,failures,logical_1_failures,logical_2_failures,seed,decoder,elapsed_seconds
3,3,3,0.0045,100,10,6,5,1001,mwpm,1.0
6,3,3,0.0045,100,15,9,8,1002,mwpm,2.0
3,3,3,0.0055,100,20,12,11,1003,mwpm,1.1
6,3,3,0.0055,100,30,18,16,1004,mwpm,2.1
5,5,5,0.0045,100,8,5,4,1005,mwpm,3.0
10,5,5,0.0045,100,12,7,6,1006,mwpm,6.0
5,5,5,0.0055,100,18,11,9,1007,mwpm,3.1
10,5,5,0.0055,100,25,15,13,1008,mwpm,6.1
"""


class ZAspectRatioAnalysisTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(
            ANALYSIS, "the H+6 final-only aspect-ratio analysis is not implemented")

    @staticmethod
    def _write_csv(folder, text=FIXTURE):
        source = Path(folder) / "toric3d-z-threshold-raw.csv"
        source.write_text(text)
        return source

    def test_wilson_retains_uncertainty_at_boundary_counts(self):
        low, high = ANALYSIS.wilson_interval(0, 1000)
        self.assertAlmostEqual(low, 0.0, places=14)
        self.assertGreater(high, 0.003)
        self.assertLess(high, 0.004)

        low, high = ANALYSIS.wilson_interval(1000, 1000)
        self.assertGreater(low, 0.996)
        self.assertAlmostEqual(high, 1.0, places=14)

    def test_cli_writes_points_differences_figures_and_no_threshold_summary(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            source = self._write_csv(folder)
            output = folder / "analysis"

            self.assertEqual(
                ANALYSIS.main([str(source), "--output-dir", str(output)]), 0)

            for csv_name in ("aspect-ratio-points.csv",
                             "aspect-ratio-differences.csv"):
                self.assertNotIn(b"\r\n", (output / csv_name).read_bytes())

            with (output / "aspect-ratio-points.csv").open(newline="") as stream:
                points = list(csv.DictReader(stream))
            self.assertEqual(len(points), 8)
            point = next(row for row in points
                         if row["L"] == "3" and row["aspect_ratio"] == "2"
                         and row["p"] == "0.0045")
            self.assertAlmostEqual(float(point["failure_rate"]), 0.15)
            self.assertLess(float(point["wilson_low"]), 0.15)
            self.assertGreater(float(point["wilson_high"]), 0.15)

            with (output / "aspect-ratio-differences.csv").open(newline="") as stream:
                differences = list(csv.DictReader(stream))
            self.assertEqual(len(differences), 4)
            difference = next(row for row in differences
                              if row["L"] == "3" and row["p"] == "0.0045")
            self.assertAlmostEqual(float(difference["failure_rate_difference"]), 0.05)
            short_low, short_high = ANALYSIS.wilson_interval(10, 100)
            long_low, long_high = ANALYSIS.wilson_interval(15, 100)
            self.assertAlmostEqual(
                float(difference["conservative_low"]), long_low - short_high)
            self.assertAlmostEqual(
                float(difference["conservative_high"]), long_high - short_low)
            self.assertEqual(difference["difference_direction"], "T=2L minus T=L")

            summary = json.loads(
                (output / "aspect-ratio-summary.json").read_text())
            self.assertEqual(summary["experiment"]["circuit"], "H+6 CNOTs")
            self.assertEqual(summary["experiment"]["recovery"], "final-only")
            self.assertEqual(summary["experiment"]["feedback"], "none (no feedback)")
            self.assertFalse(summary["threshold"]["fitted"])
            self.assertIsNone(summary["threshold"]["estimate"])
            self.assertIn("diagnostic", summary["threshold"]["reason"].lower())
            self.assertEqual(len(summary["input"]["sha256"]), 64)
            self.assertEqual(len(summary["analysis_script"]["sha256"]), 64)
            self.assertEqual(summary["sizes"], [3, 5])
            self.assertEqual(summary["rates"], [0.0045, 0.0055])

            for extension in ("pdf", "svg", "png"):
                figure = output / f"aspect-ratio-comparison.{extension}"
                self.assertGreater(figure.stat().st_size, 1000)
            svg_lines = (output / "aspect-ratio-comparison.svg").read_bytes().splitlines()
            self.assertTrue(all(line == line.rstrip() for line in svg_lines))

    def test_missing_aspect_pair_is_rejected(self):
        rows = FIXTURE.strip().splitlines()
        incomplete = "\n".join(
            row for row in rows if not row.startswith("6,3,3,0.0045,")) + "\n"
        with tempfile.TemporaryDirectory() as folder:
            source = self._write_csv(folder, incomplete)
            with self.assertRaisesRegex(ValueError, "exactly T=L and T=2L"):
                ANALYSIS.load_rows(source)

    def test_duplicate_size_rate_aspect_is_rejected(self):
        duplicate = FIXTURE + FIXTURE.strip().splitlines()[1] + "\n"
        with tempfile.TemporaryDirectory() as folder:
            source = self._write_csv(folder, duplicate)
            with self.assertRaisesRegex(ValueError, "duplicate"):
                ANALYSIS.load_rows(source)

    def test_incomplete_common_rate_grid_is_rejected(self):
        rows = FIXTURE.strip().splitlines()
        incomplete = "\n".join(
            row for row in rows
            if not (row.startswith("3,3,3,0.0055,")
                    or row.startswith("6,3,3,0.0055,"))) + "\n"
        with tempfile.TemporaryDirectory() as folder:
            source = self._write_csv(folder, incomplete)
            with self.assertRaisesRegex(ValueError, "common probability grid"):
                ANALYSIS.load_rows(source)

    def test_non_square_transverse_geometry_is_rejected(self):
        malformed = FIXTURE.replace(
            "3,3,3,0.0045,100,10", "3,3,4,0.0045,100,10", 1)
        with tempfile.TemporaryDirectory() as folder:
            source = self._write_csv(folder, malformed)
            with self.assertRaisesRegex(ValueError, "Ly = Lz"):
                ANALYSIS.load_rows(source)

    def test_invalid_counts_are_rejected(self):
        malformed = FIXTURE.replace(
            "3,3,3,0.0045,100,10", "3,3,3,0.0045,100,101", 1)
        with tempfile.TemporaryDirectory() as folder:
            source = self._write_csv(folder, malformed)
            with self.assertRaisesRegex(ValueError, "0 <= failures <= shots"):
                ANALYSIS.load_rows(source)


if __name__ == "__main__":
    unittest.main()
