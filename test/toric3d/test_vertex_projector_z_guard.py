import copy
import importlib.util
import json
import math
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "h5_z_guard", ROOT / "examples/analyze_toric3d_vertex_projector.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def audit_entry(size, protocol="none", *, mixed=0, unavoidable=0, control=0):
    counts = [1, 0, 0, 0]
    group = {
        "syndrome_indices": [1, 2],
        "upper_parity_odd": False,
        "class_counts": counts,
        "selected_class": [False, False],
        "baseline_class": [False, False],
        "mixed": False,
        "unavoidable_count": 0,
        "baseline_failure_count": 0,
        "control_failure_count": 0,
    }
    parity = {
        "mechanism_count": 1,
        "group_count": 1,
        "mixed_group_count": mixed,
        "unavoidable_count": unavoidable,
        "baseline_failure_count": 0,
        "control_failure_count": control,
    }
    empty = {key: 0 for key in parity}
    return {
        "dimensions": [size, size, size],
        "protocol": protocol,
        "final_layer_only": True,
        "encoder_operation_count": 6 * size**3,
        "active_output_fault_count": 1,
        "mechanism_count": 1,
        "group_count": 1,
        "mixed_group_count": mixed,
        "unavoidable_count": unavoidable,
        "unavoidable_rate": unavoidable,
        "baseline_failure_count": 0,
        "control_failure_count": control,
        "parity": {"even": parity, "odd": empty},
        "groups": [group],
        "mixed_groups": [],
    }


def certificate(entries):
    return {
        "format": "h5-z-boundary-ambiguity-v1",
        "status": "diagnostic_only",
        "single_fault_prior": "equal single-location Z-fault priors",
        "protocols": list(MODULE.PROTOCOLS),
        "source_hashes": {"src/example.jl": "a" * 64},
        "backend_versions": {"pymatching": "2.4.0", "ldpc": "2.4.1"},
        "audits": entries,
    }


def rows(channel="z_only"):
    return [
        {
            "channel": channel,
            "T": size,
            "Lu": size,
            "Lv": size,
            "p": probability,
            "shots": 1000,
            "none_failures": 100,
        }
        for size in (3, 5, 7, 9)
        for probability in (.01, .02, .03, .04, .05)
    ]


class H5ZBoundaryGuardTests(unittest.TestCase):
    def test_guard_returns_before_crossing_optimizer_and_bootstrap(self):
        with mock.patch.object(MODULE, "crossings") as crossing, \
                mock.patch.object(MODULE.engine, "_candidate_fit") as fit:
            result = MODULE.fit_window(
                rows(), "z_only", "none", [.01, .05], bootstrap=1000,
                z_guard_issues=["Z boundary ambiguity: unavoidable=8"])
        crossing.assert_not_called()
        fit.assert_not_called()
        self.assertEqual(result["status"], "diagnostic_candidate")
        self.assertEqual(result["bootstrap"]["accepted"], 0)
        self.assertEqual(result["bootstrap"]["samples"], [])
        self.assertIn("unavoidable=8", " ".join(result["diagnostic"]))

    def test_loader_and_clean_guard(self):
        data = certificate([
            audit_entry(size, protocol)
            for size in (3, 5, 7, 9)
            for protocol in MODULE.PROTOCOLS
        ])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "audit.json"
            path.write_text(json.dumps(data))
            loaded = MODULE.load_z_boundary_audit(path)
        self.assertEqual(
            MODULE.z_boundary_guard(
                loaded, [3, 5, 7, 9], "none",
                source_hashes={"src/example.jl": "a" * 64}), [])

    def test_guard_reports_missing_entries_counts_and_source_mismatch(self):
        data = certificate([audit_entry(3)])
        issues = MODULE.z_boundary_guard(data, [3, 5], "one_step")
        self.assertTrue(any("missing" in issue.lower() for issue in issues))

        bad = certificate([audit_entry(3, mixed=1, unavoidable=2, control=3)])
        issues = MODULE.z_boundary_guard(bad, [3], "none")
        self.assertTrue(any("mixed" in issue.lower() for issue in issues))
        self.assertTrue(any("unavoidable" in issue.lower() for issue in issues))
        self.assertTrue(any("control" in issue.lower() for issue in issues))
        issues = MODULE.z_boundary_guard(
            data, [3], "none", source_hashes={"src/example.jl": "b" * 64})
        self.assertTrue(any("source" in issue.lower() for issue in issues))

    def test_loader_rejects_malformed_counts_and_nonfinite_values(self):
        malformed = certificate([audit_entry(3)])
        malformed["audits"][0]["groups"][0]["class_counts"] = [1, 0]
        nonfinite = certificate([audit_entry(3)])
        nonfinite["audits"][0]["unavoidable_rate"] = math.nan
        duplicate = certificate([audit_entry(3), audit_entry(3)])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "audit.json"
            for value in (malformed, nonfinite, duplicate):
                path.write_text(json.dumps(value))
                with self.assertRaises(ValueError):
                    MODULE.load_z_boundary_audit(path)

    def test_clean_guard_allows_fit_and_x_only_ignores_z_issues(self):
        candidate = {
            "fitted_probabilities": [.1] * 20,
            "parameters": [0.0],
            "pc": .03,
            "nu": 1.0,
        }
        with mock.patch.object(MODULE.engine, "_candidate_fit", return_value=candidate) as fit, \
                mock.patch.object(MODULE.engine, "_candidate_guard_issues", return_value=[]):
            MODULE.fit_window(rows(), "z_only", "none", [.01, .05], bootstrap=0)
            self.assertTrue(fit.called)
        with mock.patch.object(MODULE.engine, "_candidate_fit", return_value=candidate) as fit, \
                mock.patch.object(MODULE.engine, "_candidate_guard_issues", return_value=[]):
            MODULE.fit_window(
                rows("x_only"), "x_only", "none", [.01, .05], bootstrap=0,
                z_guard_issues=["must be ignored for X"])
            self.assertTrue(fit.called)


if __name__ == "__main__":
    unittest.main()
