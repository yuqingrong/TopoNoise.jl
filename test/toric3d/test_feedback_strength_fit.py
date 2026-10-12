"""Guarded analysis tests for the diagonal three-feedback experiment."""

import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import types
import unittest
from unittest import mock

import numpy as np


SCRIPT = (Path(__file__).resolve().parents[2] / "examples" /
          "fit_toric3d_diagonal_feedback.py")
if SCRIPT.exists():
    SPEC = importlib.util.spec_from_file_location("feedback_strength_fit", SCRIPT)
    FIT = importlib.util.module_from_spec(SPEC)
    SPEC.loader.exec_module(FIT)
else:
    FIT = None


JOINT_COLUMNS = (
    "joint_n0_s0_f0", "joint_n1_s0_f0", "joint_n0_s1_f0",
    "joint_n1_s1_f0", "joint_n0_s0_f1", "joint_n1_s0_f1",
    "joint_n0_s1_f1", "joint_n1_s1_f1",
)
PROTOCOL_BITS = {"none": 1, "slice_local": 2, "full_prefix": 4}
TRUTH = {
    ("x_only", "none"): (0.020, 1.15),
    ("x_only", "slice_local"): (0.025, 1.25),
    ("x_only", "full_prefix"): (0.030, 1.35),
    ("z_only", "none"): (0.018, 1.10),
    ("z_only", "slice_local"): (0.023, 1.20),
    ("z_only", "full_prefix"): (0.028, 1.30),
}


def _marginal(joint, bit):
    return sum(value for mask, value in enumerate(joint) if mask & bit)


def synthetic_rows(seed=9182, *, shots=120000,
                   sizes=(5, 7, 9, 11, 13),
                   probabilities=None, channels=("x_only", "z_only"),
                   decreasing=False, bad_fit=False):
    if probabilities is None:
        probabilities = np.linspace(0.010, 0.040, 9)
    rng = np.random.default_rng(seed)
    reference = float(np.exp(np.mean(np.log(sizes))))
    rows = []
    for channel in channels:
        for size in sizes:
            for p_index, probability in enumerate(probabilities):
                probability = float(probability)
                marginal_rates = []
                for protocol in PROTOCOL_BITS:
                    critical, exponent = TRUTH[(channel, protocol)]
                    x = ((probability - critical) / 0.030
                         * (size / reference) ** (1 / exponent))
                    slope = -2.8 if decreasing else 2.8
                    eta = -1.15 + slope * x + 0.08 * x * x
                    rate = 1 / (1 + np.exp(-eta))
                    if bad_fit:
                        rate = np.clip(
                            rate + (0.075 if p_index % 2 else -0.075),
                            0.002, 0.998)
                    marginal_rates.append(float(rate))
                pattern_probabilities = []
                for mask in range(8):
                    value = 1.0
                    for index, bit in enumerate((1, 2, 4)):
                        rate = marginal_rates[index]
                        value *= rate if mask & bit else 1 - rate
                    pattern_probabilities.append(value)
                joint = rng.multinomial(shots, pattern_probabilities).tolist()
                marginals = {
                    protocol: _marginal(joint, bit)
                    for protocol, bit in PROTOCOL_BITS.items()
                }
                row = dict(
                    steps=size, Ly=size, Lz=size,
                    qubits=3 * (size + 1) * size * size,
                    p=probability, channel=channel, metric="state_failure",
                    shots=shots, seed=seed,
                    point_seed=f"{seed:04x}{size:04x}{p_index:08x}",
                    x_routing="one_step",
                    none_failures=marginals["none"],
                    slice_local_failures=marginals["slice_local"],
                    full_prefix_failures=marginals["full_prefix"],
                    max_iter=50, bp_method="minimum_sum",
                    ms_scaling_factor=1.0, schedule="parallel",
                    omp_thread_count=1, osd_method="OSD_CS", osd_order=2,
                    invalid_recoveries=0,
                    pymatching_version="synthetic-pymatching",
                    ldpc_version="synthetic-ldpc")
                row.update(dict(zip(JOINT_COLUMNS, joint)))
                rows.append(row)
    return rows


def write_rows(path, rows):
    fields = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


class FeedbackStrengthFitTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(FIT, "the diagonal-feedback fitter is not implemented")

    def test_parser_reads_disjoint_shards_and_derives_marginals(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            rows = synthetic_rows()
            x_rows = [row for row in rows if row["channel"] == "x_only"]
            z_rows = [row for row in rows if row["channel"] == "z_only"]
            x_rows.append(dict(x_rows[0], metric="logical_1"))
            write_rows(folder / "x.csv", x_rows)
            write_rows(folder / "z.csv", z_rows)
            points = FIT.read_feedback_points([folder / "z.csv", folder / "x.csv"])
            self.assertEqual(len(points), 2 * 5 * 9)
            self.assertEqual(points, sorted(
                points, key=lambda row: (row["channel"], row["steps"], row["p"])))
            for point in points:
                self.assertEqual(sum(point["joint_counts"]), point["shots"])
                for protocol, bit in PROTOCOL_BITS.items():
                    self.assertEqual(
                        point[f"{protocol}_failures"],
                        _marginal(point["joint_counts"], bit))
                self.assertEqual(set(point["paired_differences"]), {
                    "none_vs_slice_local", "none_vs_full_prefix",
                    "slice_local_vs_full_prefix"})

    def test_parser_rejects_corrupt_counts_duplicates_and_missing_grid(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            rows = synthetic_rows(channels=("x_only",))
            source = folder / "raw.csv"
            write_rows(source, rows)
            with self.assertRaisesRegex(ValueError, "duplicate"):
                FIT.read_feedback_points([source, source])

            broken = [dict(row) for row in rows]
            broken[0][JOINT_COLUMNS[0]] += 1
            write_rows(source, broken)
            with self.assertRaisesRegex(ValueError, "sum to shots"):
                FIT.read_feedback_points([source])

            broken = [dict(row) for row in rows]
            broken[0]["none_failures"] += 1
            write_rows(source, broken)
            with self.assertRaisesRegex(ValueError, "marginal"):
                FIT.read_feedback_points([source])

            write_rows(source, rows[:-1])
            with self.assertRaisesRegex(ValueError, "rectangular|missing"):
                FIT.read_feedback_points([source])

            broken = [dict(row) for row in rows]
            broken[0]["invalid_recoveries"] = 1
            write_rows(source, broken)
            with self.assertRaisesRegex(ValueError, "invalid recover"):
                FIT.read_feedback_points([source])

            broken = [dict(row) for row in rows]
            broken[0]["ldpc_version"] = "different-ldpc"
            write_rows(source, broken)
            with self.assertRaisesRegex(ValueError, "decoder versions"):
                FIT.read_feedback_points([source])

    def test_six_seeded_crossings_recover_known_parameters(self):
        points = FIT.validate_feedback_points(synthetic_rows())
        report = FIT.analyze_feedback_scan(
            points, bootstrap_replicates=0, include_sensitivities=False,
            p_windows={"x_only": (0.010, 0.040),
                       "z_only": (0.010, 0.040)})
        self.assertEqual(set(report["fits"]), {
            f"{channel}/{protocol}"
            for channel in ("x_only", "z_only")
            for protocol in PROTOCOL_BITS})
        for key, fit in report["fits"].items():
            channel, protocol = key.split("/")
            critical, exponent = TRUTH[(channel, protocol)]
            self.assertEqual(fit["status"], "ok", fit)
            self.assertEqual(fit["degree"], 2)
            self.assertLess(abs(fit["pc"] - critical), 0.0015)
            self.assertLess(abs(fit["nu"] - exponent), 0.40)

    def test_bootstrap_1000_is_reproducible(self):
        points = FIT.validate_feedback_points(
            synthetic_rows(channels=("x_only",)))
        fit = FIT.fit_feedback_protocol(
            points, "x_only", "none", p_min=0.010, p_max=0.040)
        one = FIT.bootstrap_feedback_protocol(
            points, "x_only", "none", fit, replicates=1000, seed=4567)
        two = FIT.bootstrap_feedback_protocol(
            points, "x_only", "none", fit, replicates=1000, seed=4567)
        self.assertEqual(one, two)
        self.assertEqual(len(one["samples"]), 1000)
        self.assertGreaterEqual(one["accepted"], 800)
        self.assertIsNotNone(one["pc_interval"])
        self.assertLess(one["pc_interval"][0], one["pc_interval"][1])
        other = FIT.bootstrap_feedback_protocol(
            points, "x_only", "none", fit, replicates=3, seed=4568)
        self.assertNotEqual(one["samples"][:3], other["samples"])

    def test_data_guards_withhold_thresholds(self):
        cases = [
            ("few sizes", synthetic_rows(
                sizes=(5, 7, 9), channels=("x_only",)),
             "insufficient_data", "four sizes"),
            ("few probabilities", synthetic_rows(
                probabilities=np.linspace(0.010, 0.040, 4),
                channels=("x_only",)), "insufficient_data", "five probabilities"),
            ("decreasing", synthetic_rows(
                decreasing=True, channels=("x_only",)),
             "inconclusive", "monotone increasing"),
            ("bad deviance", synthetic_rows(
                bad_fit=True, shots=300000, channels=("x_only",)),
             "inconclusive", "lack of fit"),
        ]
        for name, rows, status, diagnostic in cases:
            with self.subTest(name=name):
                result = FIT.fit_feedback_protocol(
                    FIT.validate_feedback_points(rows), "x_only", "none",
                    p_min=0.010, p_max=0.040)
                self.assertEqual(result["status"], status, result)
                self.assertIsNone(result["pc"])
                self.assertIsNone(result["nu"])
                self.assertTrue(any(diagnostic in item
                                    for item in result["diagnostic"]), result)

        no_crossing = synthetic_rows(channels=("x_only",))
        for row in no_crossing:
            failure = int(row["shots"] * (0.10 + 1.5 * row["p"]
                                           - 0.001 * row["steps"]))
            joint = [row["shots"] - failure, failure, 0, 0, 0, 0, 0, 0]
            row.update(dict(zip(JOINT_COLUMNS, joint)))
            row["none_failures"] = failure
            row["slice_local_failures"] = 0
            row["full_prefix_failures"] = 0
        result = FIT.fit_feedback_protocol(
            FIT.validate_feedback_points(no_crossing), "x_only", "none",
            p_min=0.010, p_max=0.040)
        self.assertEqual(result["status"], "no_crossing", result)
        self.assertIsNone(result["pc"])
        self.assertIsNone(result["nu"])
        self.assertIsNotNone(result["candidate"], result)
        self.assertTrue(math.isfinite(result["candidate"]["pc"]), result)

        missing_pair = synthetic_rows(channels=("x_only",))
        size_eleven = {
            row["p"]: row for row in missing_pair if row["steps"] == 11
        }
        for row in missing_pair:
            if row["steps"] != 13:
                continue
            source = size_eleven[row["p"]]
            joint = [source[column] for column in JOINT_COLUMNS]
            row.update(dict(zip(JOINT_COLUMNS, joint)))
            for protocol, bit in PROTOCOL_BITS.items():
                row[f"{protocol}_failures"] = _marginal(joint, bit)
        result = FIT.fit_feedback_protocol(
            FIT.validate_feedback_points(missing_pair), "x_only", "none",
            p_min=0.010, p_max=0.040)
        self.assertEqual(result["status"], "inconclusive", result)
        self.assertIsNone(result["pc"])
        self.assertIsNone(result["nu"])
        self.assertIsNotNone(result["candidate"], result)
        self.assertTrue(math.isfinite(result["candidate"]["pc"]), result)
        self.assertTrue(any("missing adjacent-size crossings" in item
                            for item in result["diagnostic"]), result)

    def test_numerical_guards_retain_rejected_candidate(self):
        points = FIT.validate_feedback_points(
            synthetic_rows(channels=("x_only",)))
        baseline = FIT.fit_feedback_protocol(
            points, "x_only", "none", p_min=0.010, p_max=0.040)
        candidate = dict(baseline["candidate"])
        cases = []
        optimizer = dict(candidate, optimizer_success=False,
                         optimizer_message="forced failure")
        cases.append((optimizer, "optimizer did not converge"))
        outside = dict(candidate, pc=0.050)
        cases.append((outside, "outside the fit window"))
        bound = dict(candidate, bound_hits=["nu"])
        cases.append((bound, "parameter bound reached"))
        for injected, diagnostic in cases:
            injected["issues"] = []
            with self.subTest(diagnostic=diagnostic), mock.patch.object(
                    FIT, "_candidate_fit", return_value=injected):
                result = FIT.fit_feedback_protocol(
                    points, "x_only", "none", p_min=0.010, p_max=0.040)
                self.assertEqual(result["status"], "inconclusive", result)
                self.assertIsNone(result["pc"])
                self.assertIsNone(result["nu"])
                self.assertEqual(result["candidate"]["pc"], injected["pc"])
                self.assertTrue(any(diagnostic in item
                                    for item in result["diagnostic"]), result)

    def test_primary_fit_exception_is_reported_as_inconclusive(self):
        points = FIT.validate_feedback_points(
            synthetic_rows(channels=("x_only",)))
        with mock.patch.object(
                FIT, "_candidate_fit",
                side_effect=RuntimeError("forced optimizer failure")):
            result = FIT.fit_feedback_protocol(
                points, "x_only", "none", p_min=0.010, p_max=0.040)
        self.assertEqual(result["status"], "inconclusive", result)
        self.assertIsNone(result["pc"])
        self.assertIsNone(result["nu"])
        self.assertIsNone(result["candidate"])
        self.assertTrue(any(
            "RuntimeError: forced optimizer failure" in item
            for item in result["diagnostic"]), result)

    def test_bootstrap_fit_exceptions_are_rejected_and_counted(self):
        points = FIT.validate_feedback_points(
            synthetic_rows(channels=("x_only",)))
        fit = FIT.fit_feedback_protocol(
            points, "x_only", "none", p_min=0.010, p_max=0.040)
        candidate = fit["candidate"]
        attempts = 0

        def fail_three_replicates(*_args, **_kwargs):
            nonlocal attempts
            replicate = attempts
            attempts += 1
            if replicate in (1, 4, 7):
                raise FloatingPointError(f"forced numerical failure {replicate}")
            return candidate

        with mock.patch.object(
                FIT, "_candidate_fit", side_effect=fail_three_replicates):
            bootstrap = FIT.bootstrap_feedback_protocol(
                points, "x_only", "none", fit, replicates=10, seed=4567)
        self.assertEqual(bootstrap["replicates"], 10)
        self.assertEqual(bootstrap["accepted"], 7)
        rejected = [sample for sample in bootstrap["samples"]
                    if not sample["accepted"]]
        self.assertEqual([sample["replicate"] for sample in rejected], [1, 4, 7])
        for sample in rejected:
            self.assertEqual(sample["status"], "rejected")
            self.assertIsNone(sample["pc"])
            self.assertIsNone(sample["nu"])
            self.assertIn("FloatingPointError: forced numerical failure",
                          sample["diagnostic"])
        self.assertIsNone(bootstrap["pc_interval"])
        self.assertIsNone(bootstrap["nu_interval"])
        guarded = FIT._finalize_bootstrap_guard(fit, bootstrap)
        self.assertEqual(guarded["status"], "inconclusive", guarded)
        self.assertTrue(any("80%" in item
                            for item in guarded["diagnostic"]), guarded)

    def test_fit_guards_do_not_swallow_process_control_exceptions(self):
        points = FIT.validate_feedback_points(
            synthetic_rows(channels=("x_only",)))
        fit = FIT.fit_feedback_protocol(
            points, "x_only", "none", p_min=0.010, p_max=0.040)
        for exception in (KeyboardInterrupt, SystemExit):
            with self.subTest(boundary="primary", exception=exception.__name__), \
                    mock.patch.object(FIT, "_candidate_fit",
                                      side_effect=exception("forced")):
                with self.assertRaises(exception):
                    FIT.fit_feedback_protocol(
                        points, "x_only", "none",
                        p_min=0.010, p_max=0.040)
            with self.subTest(boundary="bootstrap", exception=exception.__name__), \
                    mock.patch.object(FIT, "_candidate_fit",
                                      side_effect=exception("forced")):
                with self.assertRaises(exception):
                    FIT.bootstrap_feedback_protocol(
                        points, "x_only", "none", fit,
                        replicates=1, seed=4567)

    def test_bootstrap_guard_has_exact_eighty_percent_boundary(self):
        fit = {"status": "ok", "pc": 0.02, "nu": 1.0, "diagnostic": []}
        accepted = FIT._finalize_bootstrap_guard(
            fit, {"replicates": 1000, "accepted": 800, "samples": [],
                  "pc_interval": [0.019, 0.021], "nu_interval": [0.9, 1.1]})
        self.assertEqual(accepted["status"], "ok")
        rejected = FIT._finalize_bootstrap_guard(
            fit, {"replicates": 1000, "accepted": 799, "samples": [],
                  "pc_interval": None, "nu_interval": None})
        self.assertEqual(rejected["status"], "inconclusive")
        self.assertIsNone(rejected["pc"])
        self.assertTrue(any("80%" in item for item in rejected["diagnostic"]))

    def test_sensitivities_cover_required_models_and_keep_spreads_separate(self):
        points = FIT.validate_feedback_points(synthetic_rows())
        report = FIT.analyze_feedback_scan(
            points, bootstrap_replicates=0, include_sensitivities=True,
            p_windows={"x_only": (0.010, 0.040),
                       "z_only": (0.010, 0.040)})
        for key, fit in report["fits"].items():
            self.assertEqual(fit["degree"], 2)
            rows = [row for row in report["sensitivity"]
                    if row["fit_key"] == key]
            names = {row["name"] for row in rows}
            self.assertIn("cubic", names)
            self.assertTrue(any(name.startswith("min_size_") for name in names))
            self.assertIn("trim_lowest_p", names)
            self.assertIn("trim_highest_p", names)
            self.assertIn("fixed_L_inverse_correction", names)
            self.assertIn(key, report["accepted_model_pc_spread"])
            self.assertNotEqual(report["accepted_model_pc_spread"][key],
                                fit.get("bootstrap", {}).get("pc_interval"))

    def test_model_spread_is_withheld_when_primary_bootstrap_is_rejected(self):
        points = FIT.validate_feedback_points(
            synthetic_rows(channels=("x_only",)))
        rejected_bootstrap = {
            "replicates": 1000, "accepted": 799,
            "pc_interval": None, "nu_interval": None, "samples": [],
            "method": "forced bootstrap rejection",
        }
        def primary(_points, channel, protocol, **_kwargs):
            return {"channel": channel, "protocol": protocol,
                    "status": "ok", "pc": 0.03, "nu": 1.0,
                    "diagnostic": []}

        def sensitivity(_points, channel, protocol, _primary):
            return [{"fit_key": f"{channel}/{protocol}", "pc": 0.031}]

        with mock.patch.object(FIT, "fit_feedback_protocol",
                               side_effect=primary), mock.patch.object(
                FIT, "bootstrap_feedback_protocol",
                return_value=rejected_bootstrap), mock.patch.object(
                    FIT, "sensitivity_fits", side_effect=sensitivity):
            report = FIT.analyze_feedback_scan(
                points, bootstrap_replicates=1000,
                p_windows={"x_only": (0.010, 0.040)})
        self.assertTrue(all(fit["status"] == "inconclusive"
                            for fit in report["fits"].values()))
        self.assertTrue(all(spread is None for spread in
                            report["accepted_model_pc_spread"].values()))

    def test_sensitivity_figure_distinguishes_accepted_and_withheld_estimates(self):
        protocols = {
            "none": "No feedback",
            "slice_local": "Slice-local",
            "full_prefix": "Full-prefix",
        }
        fits = {}
        sensitivities = []
        for channel_index, channel in enumerate(("x_only", "z_only")):
            for protocol_index, protocol in enumerate(protocols):
                key = f"{channel}/{protocol}"
                pc = 0.02 + 0.004 * channel_index + 0.002 * protocol_index
                fits[key] = {
                    "channel": channel, "protocol": protocol,
                    "status": "ok", "pc": pc,
                    "candidate": {"pc": pc}, "diagnostic": [],
                    "bootstrap": {"pc_interval": [pc - 0.001, pc + 0.001]},
                }
                sensitivities.extend([
                    {"fit_key": key, "name": "cubic", "status": "ok",
                     "pc": pc + 0.0003, "candidate": {"pc": pc + 0.0003},
                     "diagnostic": []},
                    {"fit_key": key, "name": "trim_lowest_p",
                     "status": "inconclusive", "pc": None,
                     "candidate": {"pc": pc + 0.0015},
                     "diagnostic": ["scaling-model lack of fit"]},
                    {"fit_key": key, "name": "trim_highest_p",
                     "status": "no_crossing", "pc": None,
                     "candidate": None,
                     "diagnostic": ["no adjacent-size crossing"]},
                ])
        report = {"fits": fits, "sensitivity": sensitivities}
        fits["x_only/none"].update(status="inconclusive", pc=None,
                                   diagnostic=["bootstrap guard failed"])

        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder)
            FIT._make_sensitivity_plot(report, output)
            for extension, signature in (("pdf", b"%PDF"),
                                         ("png", b"\x89PNG"),
                                         ("svg", b"<?xml")):
                path = output / f"threshold-sensitivity.{extension}"
                self.assertGreater(path.stat().st_size, 1000)
                self.assertTrue(path.read_bytes().startswith(signature))
            svg = (output / "threshold-sensitivity.svg").read_text()
            for channel in ("X-only", "Z-only"):
                for protocol in protocols.values():
                    self.assertIn(f"{channel} / {protocol}", svg)
            for label in ("Primary accepted",
                          "Sensitivity accepted",
                          "Diagnostic candidate (WITHHELD)",
                          "No numerical candidate",
                          "primary withheld",
                          "not threshold estimates"):
                self.assertIn(label, svg)

    def test_cli_reports_are_byte_stable_and_figures_are_complete(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            rows = synthetic_rows()
            x_source = folder / "x.csv"
            z_source = folder / "z.csv"
            write_rows(x_source, [row for row in rows
                                  if row["channel"] == "x_only"])
            write_rows(z_source, [row for row in rows
                                  if row["channel"] == "z_only"])
            output = folder / "analysis"
            arguments = [str(x_source), str(z_source), "--x-p-min", "0.010",
                         "--x-p-max", "0.040", "--z-p-min", "0.010",
                         "--z-p-max", "0.040", "--bootstrap", "2",
                         "--seed", "1234", "--output-dir", str(output)]
            self.assertEqual(FIT.main(arguments), 0)
            stable_names = ["threshold-fit.json", "threshold-summary.csv",
                            "threshold-sensitivity.csv", "threshold-bootstrap.csv",
                            "threshold-points.csv", "paired-differences.csv"]
            first = {name: (output / name).read_bytes() for name in stable_names}
            self.assertEqual(FIT.main(arguments), 0)
            self.assertEqual(first, {name: (output / name).read_bytes()
                                     for name in stable_names})

            report = json.loads((output / "threshold-fit.json").read_text())
            self.assertEqual(len(report["fits"]), 6)
            self.assertEqual(report["input_sha256"], {
                str(x_source.resolve()): hashlib.sha256(x_source.read_bytes()).hexdigest(),
                str(z_source.resolve()): hashlib.sha256(z_source.read_bytes()).hexdigest(),
            })
            self.assertEqual(report["source_sha256"],
                             hashlib.sha256(SCRIPT.read_bytes()).hexdigest())
            self.assertNotIn("timestamp", report)
            self.assertIn("marginal binomial", report["likelihood"])
            self.assertIn("paired", report["protocol_dependence"])
            self.assertIn("joint counts", report["protocol_dependence"])
            self.assertIn("Bonferroni-Wilson",
                          report["paired_rate_difference_interval"])
            for stem in ("threshold-fit", "threshold-sensitivity"):
                for extension, signature in (("pdf", b"%PDF"),
                                             ("png", b"\x89PNG"),
                                             ("svg", b"<?xml")):
                    path = output / f"{stem}.{extension}"
                    self.assertGreater(path.stat().st_size, 1000)
                    self.assertTrue(path.read_bytes().startswith(signature))
            svg = (output / "threshold-fit.svg").read_text()
            for label in ("X-only", "Z-only", "No feedback",
                          "Slice-local", "Full-prefix", "Paired differences",
                          "Scaling collapse"):
                self.assertIn(label, svg)
            sensitivity_svg = (
                output / "threshold-sensitivity.svg").read_text()
            for label in ("X-only / No feedback", "Z-only / Full-prefix",
                          "Primary accepted", "Sensitivity accepted",
                          "not threshold estimates"):
                self.assertIn(label, sensitivity_svg)


if __name__ == "__main__":
    unittest.main()
