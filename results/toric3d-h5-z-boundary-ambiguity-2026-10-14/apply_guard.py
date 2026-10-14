#!/usr/bin/env python3
"""Apply and independently verify the H+5 Z-boundary pre-fit guard."""
import csv
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SPEC = importlib.util.spec_from_file_location(
    "h5_analysis", ROOT / "examples/analyze_toric3d_vertex_projector.py")
ANALYSIS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ANALYSIS)
EXPECTED = {
    3: (99, 5, 34, 8, 27, 72),
    5: (275, 9, 131, 14, 75, 200),
    7: (539, 13, 276, 20, 147, 392),
}


def verify(audit):
    identities = {}
    for entry in audit["audits"]:
        key = (entry["dimensions"][0], entry["protocol"])
        if not entry["final_layer_only"] or key in identities:
            raise ValueError("audit identities are not unique final-layer entries")
        identities[key] = entry
    for size, expected in EXPECTED.items():
        for protocol in ANALYSIS.PROTOCOLS:
            entry = identities[(size, protocol)]
            measured = (
                entry["mechanism_count"], entry["mixed_group_count"],
                entry["baseline_failure_count"], entry["control_failure_count"],
                entry["parity"]["even"]["mechanism_count"],
                entry["parity"]["odd"]["mechanism_count"],
            )
            if measured != expected:
                raise ValueError(
                    f"audit count mismatch for {(size, protocol)}: "
                    f"{measured} != {expected}")

    with (HERE / "z-boundary-groups.csv").open(newline="") as handle:
        groups = list(csv.DictReader(handle))
    if len(groups) != sum(entry["group_count"] for entry in audit["audits"]):
        raise ValueError("group CSV row count differs from JSON")

    with (HERE / "parity-erasure-control.csv").open(newline="") as handle:
        controls = list(csv.DictReader(handle))
    if len(controls) != 3 * 2 * 4:
        raise ValueError("control CSV must contain 24 rows")
    for row in controls:
        shots = int(row["shots"])
        accepted = int(row["accepted_even"])
        rejected = int(row["rejected_odd"])
        if row["label"] != "parity-erasure diagnostic" or shots != 1000:
            raise ValueError("control row label/shot count changed")
        if accepted + rejected != shots:
            raise ValueError("control acceptance counts do not sum to shots")
        for key in ("baseline_failures", "control_failures",
                    "accepted_baseline_failures", "accepted_control_failures"):
            if not 0 <= int(row[key]) <= shots:
                raise ValueError(f"invalid control count {key}")
        if int(row["accepted_baseline_failures"]) > accepted or \
                int(row["accepted_control_failures"]) > accepted:
            raise ValueError("accepted failures exceed accepted shots")

    report = (HERE / "z-boundary-report.md").read_text()
    required = ("diagnostic_only", "not `|00_L>` infidelity",
                "cannot establish a deterministic threshold")
    if not all(text in report for text in required):
        raise ValueError("Markdown report is missing required scientific labels")
    return {
        "audit_entries": len(audit["audits"]),
        "group_rows": len(groups),
        "control_rows": len(controls),
        "counts_match": True,
    }


def main():
    audit_path = HERE / "z-boundary-audit.json"
    audit = ANALYSIS.load_z_boundary_audit(audit_path)
    verification = verify(audit)
    guards = []
    for protocol in ANALYSIS.PROTOCOLS:
        issues = ANALYSIS.z_boundary_guard(
            audit, [3, 5, 7], protocol,
            source_hashes=audit["source_hashes"])
        if not issues:
            raise ValueError(f"ambiguity guard unexpectedly passed for {protocol}")
        guards.append({"protocol": protocol, "issues": issues})
    result = {
        "format": "h5-z-boundary-guard-v1",
        "status": "diagnostic_only",
        "bootstrap_requested": 1000,
        "bootstrap_executed": 0,
        "audit_sha256": ANALYSIS.sha(audit_path),
        "control_csv_sha256": ANALYSIS.sha(HERE / "parity-erasure-control.csv"),
        "guards": guards,
        "verification": verification,
    }
    ANALYSIS.atomic_json(HERE / "guard.json", result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
