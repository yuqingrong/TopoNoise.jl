#!/usr/bin/env python3
"""Assemble the staged diagonal X-feedback production scan.

The expensive L=9,11,13 points were split into independent 200-shot scout
and 800-shot production streams.  L=5,7 reuse complete 1000-shot pilot
points where available and fill the interleaved probabilities with a second
1000-shot stream.  This script pools only those declared components, checks
the resulting rectangular grid, and recomputes every reported statistic from
the integer paired counts.
"""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


SIZES = (5, 7, 9, 11, 13)
PROBABILITIES = (0.01, 0.0125, 0.015, 0.0175, 0.02, 0.0225, 0.025)
METRICS = ("logical_1", "logical_2", "state_failure")
SETTING_COLUMNS = (
    "ldpc_version", "max_iter", "bp_method", "ms_scaling_factor",
    "schedule", "omp_thread_count", "osd_method", "osd_order",
)
COUNT_COLUMNS = (
    "final_only_failures", "slice_feedback_failures", "both_failures",
)
CSV_COLUMNS = (
    "steps", "Ly", "Lz", "qubits", "p", "metric", "shots", "seed",
    "point_seed", "routing", "final_only_failures",
    "slice_feedback_failures", "both_failures", "final_only_only",
    "slice_feedback_only", "final_only_rate", "final_only_ci_low",
    "final_only_ci_high", "slice_feedback_rate", "slice_feedback_ci_low",
    "slice_feedback_ci_high", "difference_feedback_minus_final",
    "difference_standard_error", "difference_ci_low", "difference_ci_high",
    "elapsed_seconds", *SETTING_COLUMNS,
)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def wilson(failures, shots, z=1.959963984540054):
    rate = failures / shots
    denominator = 1 + z * z / shots
    center = (rate + z * z / (2 * shots)) / denominator
    half = z * math.sqrt(
        rate * (1 - rate) / shots + z * z / (4 * shots * shots)
    ) / denominator
    low = 0.0 if failures == 0 else max(0.0, center - half)
    high = 1.0 if failures == shots else min(1.0, center + half)
    return low, high


def metric_statistics(shots, final, feedback, both):
    if not (0 <= both <= min(final, feedback) <= max(final, feedback) <= shots):
        raise ValueError("invalid paired failure counts")
    if final + feedback - both > shots:
        raise ValueError("paired failure union exceeds shot count")
    final_only = final - both
    feedback_only = feedback - both
    difference = (feedback_only - final_only) / shots
    standard_error = math.sqrt(max(
        0.0,
        ((feedback_only + final_only) / shots - difference * difference)
        / (shots - 1),
    )) if shots > 1 else math.nan
    # Bonferroni combination of 97.5% Wilson intervals for discordances.
    feedback_interval = wilson(feedback_only, shots, z=2.241402727604947)
    final_interval = wilson(final_only, shots, z=2.241402727604947)
    difference_interval = (
        max(-1.0, feedback_interval[0] - final_interval[1]),
        min(1.0, feedback_interval[1] - final_interval[0]),
    )
    return {
        "final_only_only": final_only,
        "slice_feedback_only": feedback_only,
        "final_only_rate": final / shots,
        "final_only_ci_low": wilson(final, shots)[0],
        "final_only_ci_high": wilson(final, shots)[1],
        "slice_feedback_rate": feedback / shots,
        "slice_feedback_ci_low": wilson(feedback, shots)[0],
        "slice_feedback_ci_high": wilson(feedback, shots)[1],
        "difference_feedback_minus_final": difference,
        "difference_standard_error": standard_error,
        "difference_ci_low": difference_interval[0],
        "difference_ci_high": difference_interval[1],
    }


def source_manifest(base):
    """Return (path, allowed sizes) for every declared independent stream."""
    return (
        (base / "pilot/toric3d-x-feedback-raw.csv", {5, 7}),
        (base / "production-shards/L5-L7/toric3d-x-feedback-raw.csv", {5, 7}),
        (base / "extended-scout/L9/toric3d-x-feedback-raw.csv", {9}),
        (base / "production-shards/L9/toric3d-x-feedback-raw.csv", {9}),
        (base / "extended-scout/L11/toric3d-x-feedback-raw.csv", {11}),
        (base / "production-shards/L11/toric3d-x-feedback-raw.csv", {11}),
        (base / "extended-scout/L13/toric3d-x-feedback-raw.csv", {13}),
        (base / "production-shards/L13a/toric3d-x-feedback-raw.csv", {13}),
        (base / "production-shards/L13b/toric3d-x-feedback-raw.csv", {13}),
        (base / "production-shards/L13c/toric3d-x-feedback-raw.csv", {13}),
    )


def read_components(base):
    grouped = {}
    source_records = []
    settings = set()
    for path, allowed_sizes in source_manifest(base):
        if not path.is_file():
            raise FileNotFoundError(f"missing production component: {path}")
        selected_rows = 0
        selected_shots = 0
        with path.open(newline="") as stream:
            for row in csv.DictReader(stream):
                size = int(row["steps"])
                probability = float(row["p"])
                metric = row["metric"]
                if (size not in allowed_sizes or size not in SIZES
                        or probability not in PROBABILITIES or metric not in METRICS):
                    continue
                if int(row["Ly"]) != size or int(row["Lz"]) != size:
                    raise ValueError(f"noncubic component row in {path}")
                if row["routing"] != "one_step":
                    raise ValueError(f"unexpected routing in {path}: {row['routing']}")
                key = (size, probability, metric)
                setting = tuple(row[column] for column in SETTING_COLUMNS)
                settings.add(setting)
                component = {
                    "shots": int(row["shots"]),
                    "elapsed_seconds": float(row["elapsed_seconds"]),
                    "source": str(path.relative_to(base)),
                    "seed": row["seed"],
                    "point_seed": row["point_seed"],
                    **{column: int(row[column]) for column in COUNT_COLUMNS},
                }
                grouped.setdefault(key, []).append(component)
                selected_rows += 1
                selected_shots += component["shots"]
        source_records.append({
            "path": str(path.relative_to(base)),
            "sha256": sha256(path),
            "selected_metric_rows": selected_rows,
            "selected_metric_shots": selected_shots,
        })
    if len(settings) != 1:
        raise ValueError(f"component decoder settings differ: {settings}")
    return grouped, source_records, dict(zip(SETTING_COLUMNS, next(iter(settings))))


def assemble(base):
    grouped, source_records, settings = read_components(base)
    expected = {
        (size, probability, metric)
        for size in SIZES for probability in PROBABILITIES for metric in METRICS
    }
    if set(grouped) != expected:
        missing = sorted(expected - set(grouped))
        extra = sorted(set(grouped) - expected)
        raise ValueError(f"production grid mismatch; missing={missing[:5]}, extra={extra[:5]}")

    rows = []
    provenance = {}
    for key in sorted(expected):
        size, probability, metric = key
        components = grouped[key]
        shots = sum(item["shots"] for item in components)
        if shots != 1000:
            raise ValueError(f"{key} pools to {shots} shots instead of 1000")
        counts = {column: sum(item[column] for item in components)
                  for column in COUNT_COLUMNS}
        values = metric_statistics(
            shots, counts["final_only_failures"],
            counts["slice_feedback_failures"], counts["both_failures"],
        )
        point_key = f"L{size}:p={probability:.4f}"
        provenance.setdefault(point_key, [
            {
                "source": item["source"], "shots": item["shots"],
                "seed": item["seed"], "point_seed": item["point_seed"],
            }
            for item in components
        ])
        rows.append({
            "steps": size,
            "Ly": size,
            "Lz": size,
            "qubits": 3 * (size + 1) * size * size,
            "p": probability,
            "metric": metric,
            "shots": shots,
            "seed": "pooled",
            "point_seed": "pooled",
            "routing": "one_step",
            **counts,
            **values,
            "elapsed_seconds": sum(item["elapsed_seconds"] for item in components),
            **settings,
        })
    return rows, source_records, provenance, settings


def write_csv(rows, path):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def write_plot(rows, stem):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    figure, axes = plt.subplots(1, 3, figsize=(15, 5.7), constrained_layout=True)
    colors = plt.get_cmap("viridis")(np.linspace(0.08, 0.92, len(SIZES)))
    protocols = (
        ("final_only", "o", "-", "Final-only"),
        ("slice_feedback", "s", "--", "Slice feedback"),
    )
    for size, color in zip(SIZES, colors):
        state = sorted(
            (row for row in rows if row["steps"] == size and row["metric"] == "state_failure"),
            key=lambda row: row["p"],
        )
        for prefix, marker, line, label in protocols:
            x = np.asarray([row["p"] for row in state])
            y = np.asarray([row[f"{prefix}_rate"] for row in state])
            low = np.asarray([row[f"{prefix}_ci_low"] for row in state])
            high = np.asarray([row[f"{prefix}_ci_high"] for row in state])
            axes[0].errorbar(
                x, y, yerr=np.maximum(0, np.vstack((y - low, high - y))),
                color=color, marker=marker, ls=line, ms=4, lw=1.3,
                capsize=2, label=f"L={size}, {label}",
            )

        for metric, marker in (("logical_1", "o"), ("logical_2", "D")):
            logical = sorted(
                (row for row in rows if row["steps"] == size and row["metric"] == metric),
                key=lambda row: row["p"],
            )
            for prefix, _, line, _ in protocols:
                axes[1].plot(
                    [row["p"] for row in logical],
                    [row[f"{prefix}_rate"] for row in logical],
                    color=color, marker=marker, ls=line, ms=3, lw=1,
                )

        difference = np.asarray([row["difference_feedback_minus_final"] for row in state])
        low = np.asarray([row["difference_ci_low"] for row in state])
        high = np.asarray([row["difference_ci_high"] for row in state])
        axes[2].errorbar(
            [row["p"] for row in state], difference,
            yerr=np.maximum(0, np.vstack((difference - low, high - difference))),
            color=color, marker="o", lw=1.3, capsize=2, label=f"L={size}",
        )

    axes[0].set(
        title="Logical |00> state failure",
        xlabel="X-fault probability p per gate output",
        ylabel="Failure probability",
    )
    axes[0].legend(fontsize=6.5, ncol=2)
    axes[1].set(
        title="Individual logical-X parities\n(circle: logical 1, diamond: logical 2)",
        xlabel="X-fault probability p per gate output",
        ylabel="Residual logical-X probability",
    )
    axes[2].axhline(0, color="0.35", ls=":", lw=1)
    axes[2].set(
        title="Paired protocol difference",
        xlabel="X-fault probability p per gate output",
        ylabel="Slice feedback minus final-only",
    )
    axes[2].legend(fontsize=7, ncol=2)
    for axis in axes:
        axis.grid(alpha=0.2)
    figure.suptitle(
        "Encoded diagonal 3D toric-code sweep: 1,000 paired shots per point\n"
        "Ideal slice measurements/corrections; final BP+OSD uses uniform-p priors; 95% intervals",
        fontsize=11,
    )
    for extension in ("pdf", "svg", "png"):
        figure.savefig(stem.with_suffix(f".{extension}"), dpi=180)
    plt.close(figure)


def write_metadata(base, output, rows, sources, provenance, settings, script_path):
    state_rows = [row for row in rows if row["metric"] == "state_failure"]
    metadata = {
        "format_version": 1,
        "status": "complete",
        "experiment": "encoded diagonal 3D toric-code X-feedback threshold",
        "sizes": list(SIZES),
        "error_rates": list(PROBABILITIES),
        "paired_shots_per_point": 1000,
        "unique_points": len(state_rows),
        "unique_paired_shots": sum(row["shots"] for row in state_rows),
        "metric_rows": len(rows),
        "metrics": list(METRICS),
        "routing": "one_step",
        "protocols": ["final_only", "slice_feedback"],
        "logical_state": "encoded logical |00> on two logical qubits",
        "noise": "independent Bernoulli X faults after each participating H/CNOT output",
        "ideal_operations": [
            "initialization", "idle qubits", "syndrome measurements",
            "classical processing", "Pauli corrections",
        ],
        "decoder_prior": "uniform declared physical p on every final data qubit",
        "decoder_limitations": "neglects nonuniform exposure and circuit/feedback correlations",
        "decoder_settings": settings,
        "assembly_script": str(script_path.relative_to(script_path.parents[1])),
        "assembly_script_sha256": sha256(script_path),
        "production_csv_sha256": sha256(output / "toric3d-x-feedback-raw.csv"),
        "sources": sources,
        "point_provenance": provenance,
    }
    json_path = output / "toric3d-x-feedback-metadata.json"
    json_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")

    # A compact human-readable TOML companion; the JSON retains full provenance.
    toml = [
        "format_version = 1",
        'status = "complete"',
        'experiment = "encoded diagonal 3D toric-code X-feedback threshold"',
        f"sizes = {list(SIZES)}",
        f"error_rates = {list(PROBABILITIES)}",
        "paired_shots_per_point = 1000",
        f"unique_points = {len(state_rows)}",
        f"unique_paired_shots = {sum(row['shots'] for row in state_rows)}",
        f"metric_rows = {len(rows)}",
        'routing = "one_step"',
        'logical_state = "encoded logical |00> on two logical qubits"',
        f'assembly_script_sha256 = "{metadata["assembly_script_sha256"]}"',
        f'production_csv_sha256 = "{metadata["production_csv_sha256"]}"',
        "",
        "[decoder_settings]",
    ]
    for key, value in settings.items():
        numeric = value.replace(".", "", 1).isdigit()
        toml.append(f"{key} = {value if numeric else json.dumps(value)}")
    (output / "toric3d-x-feedback-metadata.toml").write_text("\n".join(toml) + "\n")


def parse_arguments(arguments=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "base_dir", type=Path,
        help="experiment directory containing pilot, scout, and production-shards",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        help="destination (default: BASE_DIR/production)",
    )
    return parser.parse_args(arguments)


def main(arguments=None):
    options = parse_arguments(arguments)
    base = options.base_dir.resolve()
    output = (options.output_dir or base / "production").resolve()
    output.mkdir(parents=True, exist_ok=True)
    rows, sources, provenance, settings = assemble(base)
    csv_path = output / "toric3d-x-feedback-raw.csv"
    write_csv(rows, csv_path)
    stem = output / "toric3d-x-feedback"
    write_plot(rows, stem)
    script_path = Path(__file__).resolve()
    write_metadata(base, output, rows, sources, provenance, settings, script_path)
    print(f"assembled {len(rows)} metric rows and {len(rows) // len(METRICS)} points")
    print(f"unique paired shots: {sum(row['shots'] for row in rows if row['metric'] == 'state_failure')}")
    print(f"saved: {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
