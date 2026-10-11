"""Guarded finite-size analysis tests for diagonal X-feedback scans."""

import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np


SCRIPT = Path(__file__).resolve().parents[2] / "examples" / "fit_toric3d_x_feedback_threshold.py"
if SCRIPT.exists():
    SPEC = importlib.util.spec_from_file_location("x_feedback_fit", SCRIPT)
    FIT = importlib.util.module_from_spec(SPEC)
    SPEC.loader.exec_module(FIT)
else:
    FIT = None


def synthetic_rows(seed=781, *, shots=30000, sizes=(5, 7, 9, 11, 13),
                   probabilities=None, decreasing=False, bad_fit=False):
    rng = np.random.default_rng(seed)
    if probabilities is None:
        probabilities = np.linspace(0.012, 0.042, 9)
    rows = []
    specifications = {
        "final_only": (0.0235, 1.25, -1.0),
        "slice_feedback": (0.0300, 1.45, -1.1),
    }
    counts = {protocol: {} for protocol in specifications}
    for protocol, (critical, exponent, intercept) in specifications.items():
        for size in sizes:
            for index, probability in enumerate(probabilities):
                x = ((probability - critical) / 0.030) * (size / 9) ** (1 / exponent)
                slope = -2.2 if decreasing else 2.2
                eta = intercept + slope * x + 0.12 * x * x
                rate = 1 / (1 + np.exp(-eta))
                if bad_fit:
                    rate = np.clip(rate + (0.10 if index % 2 else -0.10), 0.005, 0.95)
                counts[protocol][(size, float(probability))] = int(rng.binomial(shots, rate))
    for size in sizes:
        for probability in probabilities:
            probability = float(probability)
            final = counts["final_only"][(size, probability)]
            feedback = counts["slice_feedback"][(size, probability)]
            rows.append(dict(
                steps=size, Ly=size, Lz=size, qubits=3 * (size + 1) * size * size,
                p=probability, metric="state_failure", shots=shots, seed=seed,
                point_seed=f"{size:02x}{int(probability * 1e7):014x}",
                routing="one_step", final_only_failures=final,
                slice_feedback_failures=feedback,
                both_failures=min(final, feedback) // 2,
                max_iter=50, bp_method="minimum_sum", ms_scaling_factor=1.0,
                schedule="parallel", omp_thread_count=1,
                osd_method="OSD_CS", osd_order=2, ldpc_version="synthetic"))
    return rows


def write_rows(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


class XFeedbackThresholdFitTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(FIT, "the diagonal X-feedback fitter is not implemented")

    def test_read_points_filters_state_failure_and_rejects_bad_grids(self):
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder) / "raw.csv"
            rows = synthetic_rows()
            extra = dict(rows[0], metric="logical_1")
            write_rows(source, rows + [extra])
            points = FIT.read_points(source)
            self.assertEqual(len(points), 5 * 9)
            self.assertEqual({point["routing"] for point in points}, {"one_step"})

            write_rows(source, rows + [dict(rows[0])])
            with self.assertRaisesRegex(ValueError, "duplicate"):
                FIT.read_points(source)
            write_rows(source, rows[:-1])
            with self.assertRaisesRegex(ValueError, "rectangular|missing"):
                FIT.read_points(source)

    def test_protocols_recover_distinct_known_thresholds(self):
        points = FIT.validate_points(synthetic_rows())
        final = FIT.fit_protocol(points, "final_only", p_min=0.012, p_max=0.042)
        feedback = FIT.fit_protocol(points, "slice_feedback", p_min=0.012, p_max=0.042)
        self.assertEqual(final["status"], "ok", final)
        self.assertEqual(feedback["status"], "ok", feedback)
        self.assertLess(abs(final["pc"] - 0.0235), 0.0015)
        self.assertLess(abs(final["nu"] - 1.25), 0.45)
        self.assertLess(abs(feedback["pc"] - 0.0300), 0.0015)
        self.assertLess(abs(feedback["nu"] - 1.45), 0.45)
        self.assertGreater(feedback["pc"], final["pc"])

    def test_parametric_bootstrap_is_seeded_and_protocol_specific(self):
        points = FIT.validate_points(synthetic_rows())
        final_fit = FIT.fit_protocol(points, "final_only", p_min=0.012, p_max=0.042)
        feedback_fit = FIT.fit_protocol(points, "slice_feedback", p_min=0.012, p_max=0.042)
        first = FIT.bootstrap_protocol(points, "final_only", final_fit,
                                       replicates=8, seed=1234)
        second = FIT.bootstrap_protocol(points, "final_only", final_fit,
                                        replicates=8, seed=1234)
        other = FIT.bootstrap_protocol(points, "slice_feedback", feedback_fit,
                                       replicates=8, seed=1234)
        self.assertEqual(first, second)
        self.assertEqual(len(first["samples"]), 8)
        self.assertGreaterEqual(first["accepted"], 7)
        self.assertIsNotNone(first["pc_interval"])
        self.assertNotEqual(first["samples"], other["samples"])

    def test_no_crossing_does_not_emit_a_threshold(self):
        rows = synthetic_rows()
        for row in rows:
            rate = 0.08 + 1.7 * row["p"] - 0.002 * row["steps"]
            row["final_only_failures"] = int(row["shots"] * rate)
        result = FIT.fit_protocol(FIT.validate_points(rows), "final_only",
                                  p_min=0.012, p_max=0.042)
        self.assertEqual(result["status"], "no_crossing")
        self.assertIsNone(result["pc"])
        self.assertIsNone(result["nu"])

    def test_optimizer_bound_does_not_emit_a_threshold(self):
        points = FIT.validate_points(synthetic_rows())
        result = FIT.fit_protocol(points, "final_only", p_min=0.012, p_max=0.042,
                                  nu_bounds=(0.30, 0.40))
        self.assertEqual(result["status"], "inconclusive", result)
        self.assertIsNone(result["pc"])
        self.assertIn("nu", result["candidate"]["bound_hits"])

    def test_decreasing_curve_does_not_emit_a_threshold(self):
        points = FIT.validate_points(synthetic_rows(decreasing=True))
        result = FIT.fit_protocol(points, "final_only", p_min=0.012, p_max=0.042)
        self.assertEqual(result["status"], "inconclusive", result)
        self.assertIsNone(result["pc"])
        self.assertTrue(any("monotone increasing" in issue for issue in result["diagnostic"]))

    def test_too_few_sizes_or_probabilities_are_rejected(self):
        two_sizes = FIT.validate_points(synthetic_rows(sizes=(5, 7)))
        result = FIT.fit_protocol(two_sizes, "final_only", p_min=0.012, p_max=0.042)
        self.assertEqual(result["status"], "insufficient_data")
        self.assertIsNone(result["pc"])
        four_probabilities = FIT.validate_points(synthetic_rows(
            probabilities=np.linspace(0.012, 0.042, 4)))
        result = FIT.fit_protocol(four_probabilities, "final_only",
                                  p_min=0.012, p_max=0.042)
        self.assertEqual(result["status"], "insufficient_data")
        self.assertIsNone(result["pc"])

    def test_bad_deviance_does_not_emit_a_threshold(self):
        points = FIT.validate_points(synthetic_rows(bad_fit=True, shots=100000))
        result = FIT.fit_protocol(points, "final_only", p_min=0.012, p_max=0.042)
        self.assertEqual(result["status"], "inconclusive", result)
        self.assertIsNone(result["pc"])
        self.assertLess(result["candidate"]["deviance_p_value"], 0.01)

    def test_cli_writes_guarded_report_sensitivities_hashes_and_figures(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            source = folder / "raw.csv"
            full_path = folder / "full-path.csv"
            write_rows(source, synthetic_rows())
            full_rows = synthetic_rows(seed=872)
            for row in full_rows:
                row["routing"] = "full_path"
            write_rows(full_path, full_rows)
            output = folder / "analysis"
            code = FIT.main([str(source), "--full-path-csv", str(full_path),
                             "--p-min", "0.012", "--p-max", "0.042",
                             "--min-size", "5", "--bootstrap", "2",
                             "--seed", "1234", "--output-dir", str(output)])
            self.assertEqual(code, 0)
            report = json.loads((output / "threshold-fit.json").read_text())
            self.assertEqual(report["input_sha256"], hashlib.sha256(source.read_bytes()).hexdigest())
            self.assertEqual(report["source_sha256"], hashlib.sha256(SCRIPT.read_bytes()).hexdigest())
            self.assertEqual(set(report["protocols"]), {"final_only", "slice_feedback"})
            self.assertIn("full_path_fits", report)
            self.assertIn("numpy", report["package_versions"])
            self.assertIn("scipy", report["package_versions"])
            self.assertTrue(report["sensitivity"])
            self.assertTrue((output / "threshold-sensitivity.csv").is_file())
            self.assertTrue((output / "threshold-bootstrap.csv").is_file())
            self.assertTrue((output / "threshold-points.csv").is_file())
            for extension in ("pdf", "svg", "png"):
                self.assertGreater((output / f"threshold-fit.{extension}").stat().st_size, 1000)


if __name__ == "__main__":
    unittest.main()
