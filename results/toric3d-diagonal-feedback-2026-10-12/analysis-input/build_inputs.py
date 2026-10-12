#!/usr/bin/env python3
"""Build the immutable rectangular CSV set used by the guarded fit."""

from __future__ import annotations

import csv
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import shutil


REPOSITORY = Path(__file__).resolve().parents[3]
RESULT_ROOT = REPOSITORY / "results/toric3d-diagonal-feedback-2026-10-12"
DESTINATION = RESULT_ROOT / "analysis-input"
X_RATES = tuple(map(Decimal, (
    "0.025", "0.03", "0.04", "0.045", "0.05", "0.055", "0.06",
    "0.07", "0.08",
)))
Z_RATES = tuple(map(Decimal, (
    "0.025", "0.03", "0.0325", "0.035", "0.0375", "0.04", "0.045",
)))
METRICS = ("logical_1", "logical_2", "state_failure")
JOINT_COLUMNS = (
    "joint_n0_s0_f0", "joint_n1_s0_f0", "joint_n0_s1_f0",
    "joint_n1_s1_f0", "joint_n0_s0_f1", "joint_n1_s0_f1",
    "joint_n0_s1_f1", "joint_n1_s1_f1",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.relative_to(REPOSITORY))


def build_one(label: str, source_dir: Path, channel: str, size: int,
              probabilities, shots: int):
    source = source_dir / "toric3d-diagonal-feedback-raw.csv"
    metadata = source_dir / "toric3d-diagonal-feedback-metadata.toml"
    assert source.is_file() and metadata.is_file()
    with source.open(newline="") as stream:
        reader = csv.DictReader(stream)
        fieldnames = reader.fieldnames
        source_rows = list(reader)
    assert fieldnames is not None

    selected = [row for row in source_rows
                if Decimal(row["p"]) in probabilities]
    expected = {
        (size, probability, metric)
        for probability in probabilities for metric in METRICS
    }
    observed = {
        (int(row["steps"]), Decimal(row["p"]), row["metric"])
        for row in selected
    }
    assert observed == expected
    assert len(selected) == len(probabilities) * len(METRICS)
    assert all(
        int(row["steps"]) == int(row["Ly"]) == int(row["Lz"]) == size
        and row["channel"] == channel
        and int(row["shots"]) == shots
        and int(row["invalid_recoveries"]) == 0
        and sum(int(row[column]) for column in JOINT_COLUMNS) == shots
        for row in selected
    )

    output = DESTINATION / f"{label}.csv"
    if len(selected) == len(source_rows):
        shutil.copyfile(source, output)
    else:
        with output.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fieldnames,
                                    lineterminator="\n")
            writer.writeheader()
            writer.writerows(selected)

    return {
        "label": label,
        "channel": channel,
        "size": size,
        "probabilities": [str(value) for value in probabilities],
        "shots_per_point": shots,
        "source_csv": relative(source),
        "source_csv_sha256": sha256(source),
        "source_metadata": relative(metadata),
        "source_metadata_sha256": sha256(metadata),
        "source_rows": len(source_rows),
        "derived_csv": relative(output),
        "derived_csv_sha256": sha256(output),
        "derived_rows": len(selected),
        "byte_copy": len(selected) == len(source_rows),
    }


def main():
    DESTINATION.mkdir(parents=True, exist_ok=True)
    inputs = []
    for size in (5, 7, 9, 11):
        inputs.append(build_one(
            f"x-L{size}", RESULT_ROOT / f"production/x-L{size}",
            "x_only", size, X_RATES, 1000,
        ))
    for size in (9, 13, 17):
        inputs.append(build_one(
            f"z-L{size}", RESULT_ROOT / f"production/z-L{size}",
            "z_only", size, Z_RATES, 2000,
        ))
    inputs.append(build_one(
        "z-L21", RESULT_ROOT / "production/z-L21-p045",
        "z_only", 21, Z_RATES, 2000,
    ))

    manifest = {
        "format_version": 1,
        "purpose": "immutable rectangular inputs for the guarded production fit",
        "selection_file": relative(RESULT_ROOT / "scan-selection.toml"),
        "selection_sha256": sha256(RESULT_ROOT / "scan-selection.toml"),
        "runtime_scope_decision": (
            "Use complete four-size grids. Retain all nine X probabilities and "
            "the first seven frozen Z probabilities through p=0.045; the latter "
            "preserve an above-transition guard point and both endpoint trims."
        ),
        "x_state_failure_points": 36,
        "x_paired_fault_records": 36000,
        "z_state_failure_points": 28,
        "z_paired_fault_records": 56000,
        "inputs": inputs,
    }
    path = DESTINATION / "manifest.json"
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(relative(path))


if __name__ == "__main__":
    main()
