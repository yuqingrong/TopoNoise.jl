#!/usr/bin/env python3
"""Regression tests for retrying explicit numerical invalidity at larger chi."""

from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import core
import manifest
import adaptive_runner
from test_phase2 import FakeExecutor, make_bundle, tn_result


def complete(chi, *, logical_mask=0):
    probabilities = [0.1, 0.1, 0.1, 0.1]
    probabilities[logical_mask] = 0.7
    return {
        "chi": chi, "status": "complete", "valid": True,
        "probabilities": probabilities, "logical_mask": logical_mask,
        "log_joint_syndrome": -10.0,
    }


def numerical_invalid(chi):
    return {"chi": chi, "status": "invalid", "valid": False}


class InvalidEscalationTests(unittest.TestCase):
    def setUp(self):
        self.ladder = [32, 64, 128, 240, 480, 960]

    def test_explicit_numerical_invalidity_tries_the_next_chi(self):
        for checkpoint in (numerical_invalid(32),
                           {**numerical_invalid(32), "error": ""}):
            with self.subTest(checkpoint=checkpoint):
                decision = core.evaluate_stopping([checkpoint], self.ladder)
                self.assertEqual(decision["status"], "continue")
                self.assertEqual(decision["next_chi"], 64)

    def test_repeated_invalid_levels_do_not_shift_the_ladder_index(self):
        levels = [numerical_invalid(32), numerical_invalid(64),
                  complete(128), complete(240)]
        decision = core.evaluate_stopping(levels, self.ladder)
        self.assertEqual(decision["status"], "continue")
        self.assertEqual(decision["next_chi"], 480)
        decision = core.evaluate_stopping(levels + [complete(480)], self.ladder)
        self.assertEqual(decision["status"], "accepted")
        self.assertEqual(decision["accepted_chi"], 480)

    def test_invalidity_breaks_the_consecutive_valid_window(self):
        levels = [complete(32), numerical_invalid(64),
                  complete(128), complete(240)]
        decision = core.evaluate_stopping(levels, self.ladder)
        self.assertEqual(decision["status"], "continue")
        self.assertEqual(decision["next_chi"], 480)
        decision = core.evaluate_stopping(levels + [complete(480)], self.ladder)
        self.assertEqual(decision["status"], "accepted")
        self.assertEqual(decision["accepted_chi"], 480)

    def test_malformed_results_and_diagnostic_invalidity_remain_terminal(self):
        malformed_complete = complete(32)
        malformed_complete["probabilities"] = [float("nan"), 0.1, 0.1, 0.1]
        cases = [
            malformed_complete,
            {**numerical_invalid(32), "error": "malformed TN output"},
            {"chi": 32, "status": "invalid", "valid": True},
            {"chi": 32, "status": "invalid"},
        ]
        for checkpoint in cases:
            with self.subTest(checkpoint=checkpoint):
                decision = core.evaluate_stopping([checkpoint], self.ladder)
                self.assertEqual(decision["status"], "unknown")
                self.assertNotIn("next_chi", decision)

    def test_exceptions_and_resource_failures_remain_terminal(self):
        for status in ("exception", "timeout", "memory_limit", "interrupted",
                       "missing_checkpoint"):
            with self.subTest(status=status):
                decision = core.evaluate_stopping(
                    [{"chi": 32, "status": status, "valid": False}], self.ladder
                )
                self.assertEqual(decision, {"status": "unknown", "reason": status})

    def test_malformed_error_field_cannot_be_normalized_into_a_retry(self):
        for error in ({"message": "bad output"}, ["bad output"], 1, 0, None, False):
            with self.subTest(error=error):
                raw = {**numerical_invalid(32), "error": error}
                self.assertEqual(core.evaluate_stopping([raw], self.ladder)["status"],
                                 "unknown")
                normalized = adaptive_runner._normalize_executor_result(raw, label="TN")
                self.assertIn("malformed", normalized["error"])
                decision = core.evaluate_stopping([{**normalized, "chi": 32}], self.ladder)
                self.assertEqual(decision["status"], "unknown")

    def test_invalid_maximum_cannot_resurrect_provisional_anchor_acceptance(self):
        ladder = [32, 64, 128, 240]
        levels = [complete(32), complete(64), complete(128)]
        provisional = core.evaluate_stopping(levels, ladder, forced_anchor=True)
        self.assertEqual(provisional["status"], "continue_anchor")
        final = core.evaluate_stopping(
            levels + [numerical_invalid(240)], ladder, forced_anchor=True
        )
        self.assertEqual(final["status"], "unknown")
        self.assertNotIn("accepted_chi", final)

    def test_anchor_requires_new_three_level_window_after_invalidity(self):
        ladder = [32, 64, 128, 240, 480, 960, 1920]
        levels = [complete(32), complete(64), complete(128), numerical_invalid(240),
                  complete(480, logical_mask=1), complete(960, logical_mask=1)]
        provisional = core.evaluate_stopping(levels, ladder, forced_anchor=True)
        self.assertEqual(provisional["status"], "continue")
        self.assertEqual(provisional["next_chi"], 1920)
        self.assertNotIn("provisional_acceptance", provisional)
        final = core.evaluate_stopping(
            levels + [complete(1920, logical_mask=1)], ladder, forced_anchor=True
        )
        self.assertEqual(final["status"], "accepted")
        self.assertEqual(final["accepted_chi"], 1920)
        self.assertEqual(final["logical_mask"], 1)
        self.assertTrue(final["anchor_completed"])
        self.assertFalse(final["anchor_class_changed"])

    def test_manifest_binds_retry_policy_and_rejects_old_policy_identity(self):
        current = manifest.build_manifest(202610600000)
        policy = current["stopping_policy"]
        self.assertEqual(policy.get("version"), "tnml-adaptive-stopping-v2")
        self.assertEqual(policy.get("invalid_action"),
                         "retry_larger_chi_reset_consecutive_window")
        self.assertEqual(policy["minimum_completed_levels"], 3)
        self.assertEqual(policy["minimum_consecutive_valid_levels"], 3)
        self.assertEqual(policy["posterior_bound"], 1e-6)
        self.assertEqual(policy["log_joint_bound"], 1e-6)
        self.assertEqual(manifest.validate_manifest(current), current)
        old = copy.deepcopy(current)
        old["stopping_policy"].pop("version")
        old["stopping_policy"].pop("invalid_action")
        old["stopping_policy"].pop("minimum_consecutive_valid_levels")
        old["manifest_id"] = core.identity_hash(old, "manifest_id")
        self.assertNotEqual(old["manifest_id"], current["manifest_id"])
        with self.assertRaisesRegex(ValueError, "specification"):
            manifest.validate_manifest(old)

    def test_runner_preserves_invalid_checkpoint_retries_and_resumes_without_work(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling, _, execution = make_bundle(root)
            record = execution["tasks"][0]["records"][0]
            executor = FakeExecutor(tn_factory=lambda chi: (
                tn_result(chi, status="invalid") if chi == 32 else tn_result(chi)
            ))
            decision = adaptive_runner.run_record(
                root, sampling, task_id=0, record_id=record["record_id"],
                executor=executor,
            )
            self.assertEqual([call[1] for call in executor.tn_calls],
                             [32, 64, 128, 240])
            self.assertEqual(len(executor.exact_calls), 1)
            self.assertEqual(decision["tested_chi"], [32, 64, 128, 240])
            self.assertEqual(decision["final_status"], "accepted")
            self.assertEqual(decision["adaptive_decision"]["accepted_chi"], 240)
            self.assertTrue(decision["exact_audit"]["passed"])
            self.assertFalse(decision["exact_audit"]["false_accept"])

            directory = Path(decision["record_directory"])
            invalid = json.loads((directory / "chi_0032.json").read_text())
            self.assertEqual(invalid["status"], "invalid")
            self.assertFalse(invalid["valid"])
            self.assertFalse(invalid.get("error"))
            self.assertNotIn("probabilities", invalid)
            for chi in (64, 128, 240):
                checkpoint = json.loads((directory / f"chi_{chi:04d}.json").read_text())
                self.assertEqual(checkpoint["status"], "complete")
                self.assertTrue(checkpoint["valid"])
            saved = {path.name: path.read_bytes() for path in directory.glob("*.json")}

            resumed = FakeExecutor()
            result = adaptive_runner.run_record(
                root, sampling, task_id=0, record_id=record["record_id"],
                executor=resumed,
            )
            self.assertEqual(result, decision)
            self.assertEqual(resumed.tn_calls, [])
            self.assertEqual(resumed.exact_calls, [])
            self.assertEqual(
                {path.name: path.read_bytes() for path in directory.glob("*.json")}, saved
            )


if __name__ == "__main__":
    unittest.main()
