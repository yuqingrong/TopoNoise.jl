#!/usr/bin/env python3
"""Analyze the original Cartesian-sweep 3D toric-code scan at small Z-noise p."""

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile

import numpy as np


PROTOCOL_COLUMNS = {
    "final_only": "final_only_failures",
    "layer_feedback": "layer_feedback_failures",
}


def wilson_interval(failures, shots, z=1.959963984540054):
    rate = failures / shots
    denominator = 1 + z * z / shots
    center = (rate + z * z / (2 * shots)) / denominator
    halfwidth = z * math.sqrt(rate * (1 - rate) / shots + z * z / (4 * shots * shots)) / denominator
    return max(0.0, center - halfwidth), min(1.0, center + halfwidth)


def target_crossing(points, target):
    """Piecewise-linear p at which a monotone sampled curve reaches target."""
    if not 0 < target < 0.5:
        raise ValueError("target must lie strictly between zero and one half")
    ordered = sorted((float(p), float(rate)) for p, rate in points)
    for p, rate in ordered:
        if rate == target:
            return p
    for (p0, q0), (p1, q1) in zip(ordered[:-1], ordered[1:]):
        if q0 <= target <= q1 and q1 > q0:
            return p0 + (target - q0) * (p1 - p0) / (q1 - q0)
    raise ValueError("sampled probabilities do not bracket the target failure rate")


def fit_power_law(values):
    """Fit p_target(L)=coefficient*L^(-exponent) in log coordinates."""
    if len(values) < 2 or any(L <= 0 or p <= 0 for L, p in values):
        raise ValueError("power-law fitting needs at least two positive (L,p) values")
    slope, intercept = np.polyfit(
        np.log([L for L, _ in values]), np.log([p for _, p in values]), 1)
    predicted = np.exp(intercept + slope * np.log([L for L, _ in values]))
    observed = np.array([p for _, p in values])
    residual = np.log(observed) - np.log(predicted)
    return {
        "exponent": float(-slope),
        "coefficient": float(np.exp(intercept)),
        "log_residual_rms": float(np.sqrt(np.mean(residual * residual))),
    }


def bound_target_probability(L, target):
    """p at which [1-(1-2p)^(L^2)]/2 equals target."""
    return (1 - (1 - 2 * target) ** (1 / (L * L))) / 2


def load_rows(path):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    selected = []
    for row in rows:
        if row["channel"].strip('"') != "z_only" or row["metric"].strip('"') != "logical_z":
            continue
        selected.append({
            "L": int(row["Lx"]),
            "p": float(row["p"]),
            "shots": int(row["shots"]),
            "final_only_failures": int(row["final_only_failures"]),
            "layer_feedback_failures": int(row["layer_feedback_failures"]),
            "both_failures": int(row["both_failures"]),
            "bound": float(row["logical_z_diagnostic_bound"]),
        })
    if not selected:
        raise ValueError("CSV has no z_only/logical_z rows")
    keys = [(row["L"], row["p"]) for row in selected]
    if len(keys) != len(set(keys)):
        raise ValueError("CSV contains duplicate size/probability rows")
    return sorted(selected, key=lambda row: (row["L"], row["p"]))


def estimates(rows, target):
    sizes = sorted({row["L"] for row in rows})
    output = {}
    for protocol, column in PROTOCOL_COLUMNS.items():
        values = []
        for L in sizes:
            points = [(row["p"], row[column] / row["shots"])
                      for row in rows if row["L"] == L]
            values.append((L, target_crossing(points, target)))
        output[protocol] = {
            "values": [{"L": L, "p_target": p, "p_target_times_L2": p * L * L}
                       for L, p in values],
            "power_law": fit_power_law(values),
        }
    return output


def bootstrap(rows, target, count, seed):
    rng = np.random.default_rng(seed)
    sizes = sorted({row["L"] for row in rows})
    samples = []
    for replicate in range(count):
        record = {"replicate": replicate, "accepted": True}
        try:
            for protocol, column in PROTOCOL_COLUMNS.items():
                values = []
                for L in sizes:
                    points = []
                    for row in rows:
                        if row["L"] != L:
                            continue
                        probability = row[column] / row["shots"]
                        failures = rng.binomial(row["shots"], probability)
                        points.append((row["p"], failures / row["shots"]))
                    values.append((L, target_crossing(points, target)))
                fit = fit_power_law(values)
                record[f"{protocol}_exponent"] = fit["exponent"]
                for L, p in values:
                    record[f"{protocol}_p_L{L}"] = p
        except ValueError:
            record["accepted"] = False
        samples.append(record)
    return samples


def quantile_interval(values):
    return [float(value) for value in np.quantile(values, [0.025, 0.5, 0.975])]


def make_plot(rows, report, output):
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "toric3d-matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sizes = sorted({row["L"] for row in rows})
    colors = plt.get_cmap("viridis")(np.linspace(0.08, 0.9, len(sizes)))
    figure, axes = plt.subplots(1, 3, figsize=(15.2, 4.8), constrained_layout=True)
    for axis, (protocol, column), title in zip(
            axes[:2], PROTOCOL_COLUMNS.items(),
            ("Final-only recovery", "Recovery after every x-layer")):
        for L, color in zip(sizes, colors):
            points = [row for row in rows if row["L"] == L]
            p = np.array([row["p"] for row in points])
            rate = np.array([row[column] / row["shots"] for row in points])
            intervals = np.array([wilson_interval(row[column], row["shots"])
                                  for row in points])
            errors = np.vstack([rate - intervals[:, 0], intervals[:, 1] - rate])
            axis.errorbar(p, rate, yerr=errors, marker="o", ms=4, capsize=2,
                          lw=1.4, color=color, label=f"L={L}")
        axis.axhline(report["target_failure_rate"], color="black", ls=":", lw=1,
                     label=f"{100 * report['target_failure_rate']:g}% failure")
        axis.set_xscale("log")
        axis.set(xlabel="Z fault probability per active gate output, p",
                 ylabel="Residual logical-Z probability", title=title)
        axis.grid(alpha=0.22)
    axes[0].legend(fontsize=8, ncol=2)

    axis = axes[2]
    markers = {"final_only": "o", "layer_feedback": "s"}
    labels = {"final_only": "Final-only", "layer_feedback": "Layer feedback"}
    for protocol, color in (("final_only", "#d1495b"), ("layer_feedback", "#2364aa")):
        summary = report["estimates"][protocol]
        L = np.array([item["L"] for item in summary["values"]], dtype=float)
        p = np.array([item["p_target"] for item in summary["values"]])
        low = np.array([item["bootstrap_interval"][0] for item in summary["values"]])
        high = np.array([item["bootstrap_interval"][2] for item in summary["values"]])
        axis.errorbar(L, p, yerr=np.vstack([p - low, high - p]), marker=markers[protocol],
                      color=color, capsize=3, ls="none", label=labels[protocol])
        fit = summary["power_law"]
        grid = np.linspace(min(L), max(L), 150)
        axis.plot(grid, fit["coefficient"] * grid ** (-fit["exponent"]),
                  color=color, lw=1.5, label=f"fit L^(-{fit['exponent']:.2f})")
    bound = np.array([bound_target_probability(int(L), report["target_failure_rate"]) for L in sizes])
    axis.plot(sizes, bound, color="black", ls="--", lw=1.3,
              label="analytic-bound ceiling")
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set(xlabel="Cubic size L", ylabel="p at 1% logical-Z failure",
             title="Pseudothreshold moves toward zero")
    axis.grid(alpha=0.22, which="both")
    axis.legend(fontsize=7.6)
    figure.suptitle(
        "Original Cartesian x-sweep: small-p Z-noise scan\n"
        "10,000 paired shots per point; logical Z is a frame diagnostic for |0_L>",
        fontsize=11)
    for extension in ("pdf", "svg", "png"):
        figure.savefig(output / f"original-small-p.{extension}", dpi=190)
    plt.close(figure)


def write_csv(path, rows, fields):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main(arguments=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--target", type=float, default=0.01)
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=24680)
    args = parser.parse_args(arguments)
    if args.bootstrap < 0:
        parser.error("--bootstrap must be nonnegative")
    rows = load_rows(args.csv)
    primary = estimates(rows, args.target)
    samples = bootstrap(rows, args.target, args.bootstrap, args.seed)
    accepted = [sample for sample in samples if sample["accepted"]]
    if args.bootstrap and len(accepted) < 0.8 * args.bootstrap:
        raise RuntimeError("fewer than 80% of bootstrap replicas bracket the target")

    for protocol in PROTOCOL_COLUMNS:
        for item in primary[protocol]["values"]:
            key = f"{protocol}_p_L{item['L']}"
            item["bootstrap_interval"] = quantile_interval([sample[key] for sample in accepted])
        primary[protocol]["power_law"]["bootstrap_exponent_interval"] = quantile_interval(
            [sample[f"{protocol}_exponent"] for sample in accepted])

    sensitivity = []
    for target in (0.005, 0.01, 0.02):
        try:
            result = estimates(rows, target)
            sensitivity.append({
                "target": target,
                "final_only_exponent": result["final_only"]["power_law"]["exponent"],
                "layer_feedback_exponent": result["layer_feedback"]["power_law"]["exponent"],
            })
        except ValueError:
            sensitivity.append({"target": target, "status": "not bracketed"})

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    metadata = sorted(args.csv.parent.glob("*metadata.toml"))
    report = {
        "status": "zero_threshold_supported",
        "conclusion": "No positive logical-Z diagnostic threshold: the fixed-failure pseudothreshold decreases with L, and the analytic lower bound tends to one half for every fixed p>0.",
        "scope": "Original Cartesian x-sweep, Z-only circuit faults, final-only and full-prefix layer feedback; logical Z is not |0_L> preparation infidelity.",
        "target_failure_rate": args.target,
        "interpolation": "piecewise linear in p between adjacent simulated rates bracketing the target",
        "estimates": primary,
        "sensitivity": sensitivity,
        "bootstrap": {
            "requested": args.bootstrap,
            "accepted": len(accepted),
            "seed": args.seed,
            "method": "independent binomial resampling at observed point rates; intervals are sampling uncertainty conditional on interpolation and power-law fit",
        },
        "analytic_bound": "P_LZ >= [1-(1-2p)^(L^2)]/2 for cubic L; hence lim_{L->infinity} P_LZ >= 1/2 at every fixed p>0",
        "input_csv": str(args.csv.resolve()),
        "input_sha256": hashlib.sha256(args.csv.read_bytes()).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "metadata": [{"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                     for path in metadata],
    }
    (output / "original-small-p-analysis.json").write_text(json.dumps(report, indent=2) + "\n")
    write_csv(output / "original-small-p-bootstrap.csv", samples,
              sorted({key for sample in samples for key in sample}))
    point_rows = []
    for row in rows:
        item = dict(row)
        for protocol, column in PROTOCOL_COLUMNS.items():
            low, high = wilson_interval(row[column], row["shots"])
            item[f"{protocol}_rate"] = row[column] / row["shots"]
            item[f"{protocol}_wilson_low"] = low
            item[f"{protocol}_wilson_high"] = high
        point_rows.append(item)
    write_csv(output / "original-small-p-points.csv", point_rows, list(point_rows[0]))
    make_plot(rows, report, output)
    print(json.dumps({
        "status": report["status"],
        "target": args.target,
        "final_only": primary["final_only"],
        "layer_feedback": primary["layer_feedback"],
        "accepted_bootstraps": len(accepted),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
