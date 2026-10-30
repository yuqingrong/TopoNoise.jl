#!/usr/bin/env python3
"""Inspect one saved diagnostic ladder using the actual production stopping rule.

This deliberately runs only after all six levels have been measured. Prefix
timing is reconstructed from that diagnostic, not a measured adaptive runtime.
No production manifest, acceptance gate, or existing study is created/edited.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

try:
    from . import adaptive_runner, core
except ImportError:
    import adaptive_runner
    import core


# The L3 exact-reference audit gates in execution_manifest.py. These are audit
# gates, separate from core.evaluate_stopping's stricter convergence tolerances.
EXACT_AUDIT_POLICY = {
    "ml_tie_tolerance": 1e-9,
    "posterior_linf_gate": 1e-5,
    "log_joint_gate": 1e-5,
}
LADDER = [32, 64, 128, 240, 480, 960]


def summarize_documents(
    run: Mapping[str, Any], exact: Mapping[str, Any],
    levels: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    if run.get("record_authenticated") is not True:
        raise ValueError("Julia record authentication did not complete")
    if run.get("chi_ladder") != LADDER or [item.get("chi") for item in levels] != LADDER:
        raise ValueError("diagnostic must contain the full ladder in order")
    record_id = run["record_id"]
    if any(item.get("record_id") != record_id for item in [exact, *levels]):
        raise ValueError("diagnostic record identity mismatch")
    for item in [exact, *levels]:
        if item.get("status") not in ("complete", "invalid"):
            raise ValueError(f"diagnostic failed with {item.get('status')}: {item.get('error', '')}")
        if item.get("status") == "complete":
            checked = core._eligible_level({"chi": item.get("chi", 1), **item})
            if checked["status"] != "complete":
                raise ValueError("malformed complete diagnostic result")
        elif item.get("valid") is not False:
            raise ValueError("invalid diagnostic level must have valid=false")
        elapsed = item.get("elapsed_seconds")
        if isinstance(elapsed, bool) or not isinstance(elapsed, (float, int)) or not math.isfinite(elapsed) or elapsed < 0:
            raise ValueError("elapsed_seconds must be finite and nonnegative")
    if exact.get("status") != "complete":
        raise ValueError("exact reference is invalid")

    prefix_decisions = []
    for count in range(1, len(levels) + 1):
        decision = core.evaluate_stopping(levels[:count], LADDER)
        prefix_decisions.append({"through_chi": LADDER[count - 1], "decision": decision})
    # An actual adaptive run would stop at its first accepted prefix. Later
    # diagnostic levels must not retroactively change that recorded decision.
    adaptive = next(
        (item["decision"] for item in prefix_decisions if item["decision"]["status"] == "accepted"),
        prefix_decisions[-1]["decision"],
    )
    prefix_seconds = None
    if adaptive["status"] == "accepted":
        stop_index = LADDER.index(adaptive["accepted_chi"]) + 1
        prefix_seconds = sum(item["elapsed_seconds"] for item in levels[:stop_index])
    # Apply the same class/posterior/log-evidence audit as the production runner.
    accepted_audit = adaptive_runner._exact_audit(adaptive, exact, EXACT_AUDIT_POLICY)
    highest = levels[-1]
    high_decision = (
        {**highest, "status": "accepted"}
        if highest["status"] == "complete" and highest["valid"] is True
        else {"status": "unknown", "reason": "highest_chi_invalid"}
    )
    high_audit = adaptive_runner._exact_audit(high_decision, exact, EXACT_AUDIT_POLICY)
    return {
        "kind": "standalone_full_ladder_diagnostic",
        "L": run["L"], "p": run["p"], "p_key": run["p_key"],
        "seed": run["seed"], "record_id": record_id,
        "invalid_chis": [item["chi"] for item in levels if item["status"] == "invalid"],
        "adaptive_decision": adaptive,
        "prefix_decisions": prefix_decisions,
        "accepted_exact_audit": accepted_audit,
        "highest_chi_exact_audit": high_audit,
        "exact_audit_policy": EXACT_AUDIT_POLICY,
        "full_ladder_tn_seconds": sum(item["elapsed_seconds"] for item in levels),
        "measured_prefix_tn_seconds": prefix_seconds,
        "startup_seconds": run["startup_seconds"],
        "startup_timing_scope": run.get("startup_timing_scope", "unspecified"),
        "exact_seconds": exact["elapsed_seconds"],
        "timing_scope": "full_ladder_diagnostic_build_plus_decode; prefix_is_retrospective_not_adaptive_runtime",
        "first_level_includes_tn_compilation": True,
        "production_gate": False,
    }


def inspect_directory(directory: Path) -> dict[str, Any]:
    run = adaptive_runner._read_toml(directory / "run.toml")
    record_path = directory / "record.toml"
    if hashlib.sha256(record_path.read_bytes()).hexdigest() != run["record_file_sha256"]:
        raise ValueError("saved record content changed after Julia authentication")
    record = adaptive_runner._read_toml(record_path)
    if record["record_id"] != run["record_id"]:
        raise ValueError("saved record identity differs from run")
    here = Path(__file__).resolve().parent
    expected_sources = {"core.py", "adaptive_runner.py", "inspect_adaptive_ladder.py"}
    if set(run["policy_source_sha256"]) != expected_sources:
        raise ValueError("policy source inventory is incomplete")
    for name, expected in run["policy_source_sha256"].items():
        if hashlib.sha256((here / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"policy source changed after diagnostic started: {name}")
    exact = adaptive_runner._read_toml(directory / "exact.toml")
    levels = [adaptive_runner._read_toml(directory / f"chi_{chi:04d}.toml") for chi in LADDER]
    summary = summarize_documents(run, exact, levels)
    summary["diagnostic_directory"] = str(directory)
    summary["runtime_identity"] = run["runtime_identity"]
    summary["policy_source_sha256"] = run["policy_source_sha256"]
    summary["level_seconds"] = [
        {key: level[key] for key in ("chi", "status", "elapsed_seconds", "decode_seconds", "decoder_build_seconds")}
        for level in levels
    ]
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    directory = args.directory.resolve()
    summary = inspect_directory(directory)
    serialized = json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n"
    # The Slurm wrapper creates a fresh directory; never replace a saved result.
    with (directory / "summary.json").open("x", encoding="utf-8") as stream:
        stream.write(serialized)
    print(serialized, end="", flush=True)
    audit = summary["accepted_exact_audit"]
    return 2 if audit["false_accept"] or summary["highest_chi_exact_audit"].get("passed") is not True else 0


if __name__ == "__main__":
    raise SystemExit(main())
