#!/usr/bin/env python3
"""Compare T=L and T=2L for the H+6-CNOT final-only Z experiment.

This is an aspect-ratio diagnostic. It reports fixed-size binomial estimates
and the unpaired difference ``T=2L minus T=L``. It does not fit or report a
threshold.

Example:
  python examples/analyze_toric3d_z_aspect_ratio.py raw.csv \
      --output-dir results/aspect-ratio-analysis
"""

import argparse
import csv
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile


REQUIRED_COLUMNS = (
    "steps",
    "Ly",
    "Lz",
    "p",
    "shots",
    "failures",
    "logical_1_failures",
    "logical_2_failures",
    "seed",
    "decoder",
    "elapsed_seconds",
)


def wilson_interval(failures, shots, z=1.959963984540054):
    """Return a two-sided 95% Wilson score interval for binomial counts."""
    if shots <= 0 or not 0 <= failures <= shots:
        raise ValueError(
            "Wilson counts require 0 <= failures <= shots and shots > 0")
    probability = failures / shots
    denominator = 1 + z * z / shots
    center = (probability + z * z / (2 * shots)) / denominator
    half = z * math.sqrt(
        probability * (1 - probability) / shots
        + z * z / (4 * shots * shots)) / denominator
    return max(0.0, center - half), min(1.0, center + half)


def _integer(row, name, row_number):
    try:
        value = float(row[name])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"row {row_number}: {name} must be an integer") from error
    if not math.isfinite(value) or not value.is_integer():
        raise ValueError(f"row {row_number}: {name} must be an integer")
    return int(value)


def _finite_float(row, name, row_number):
    try:
        value = float(row[name])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"row {row_number}: {name} must be finite") from error
    if not math.isfinite(value):
        raise ValueError(f"row {row_number}: {name} must be finite")
    return value


def load_rows(path):
    """Read and validate the rectangular T=L/T=2L comparison grid."""
    path = Path(path)
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        missing = [name for name in REQUIRED_COLUMNS
                   if name not in (reader.fieldnames or ())]
        if missing:
            raise ValueError("CSV is missing required columns: "
                             + ", ".join(missing))
        source_rows = list(reader)

    if not source_rows:
        raise ValueError("CSV contains no data rows")

    rows = []
    keys = set()
    for row_number, source in enumerate(source_rows, start=2):
        steps = _integer(source, "steps", row_number)
        Ly = _integer(source, "Ly", row_number)
        Lz = _integer(source, "Lz", row_number)
        shots = _integer(source, "shots", row_number)
        failures = _integer(source, "failures", row_number)
        logical_1 = _integer(source, "logical_1_failures", row_number)
        logical_2 = _integer(source, "logical_2_failures", row_number)
        seed = _integer(source, "seed", row_number)
        probability = _finite_float(source, "p", row_number)
        elapsed = _finite_float(source, "elapsed_seconds", row_number)
        decoder = source["decoder"].strip()

        if Ly <= 0 or Ly != Lz:
            raise ValueError(
                f"row {row_number}: comparison requires Ly = Lz = L > 0")
        if steps not in (Ly, 2 * Ly):
            raise ValueError(
                f"row {row_number}: steps must equal T=L or T=2L")
        if not 0 <= probability <= 1:
            raise ValueError(f"row {row_number}: p must lie in [0, 1]")
        if shots <= 0 or not 0 <= failures <= shots:
            raise ValueError(
                f"row {row_number}: counts require 0 <= failures <= shots "
                "and shots > 0")
        for name, count in (("logical_1_failures", logical_1),
                            ("logical_2_failures", logical_2)):
            if not 0 <= count <= shots:
                raise ValueError(
                    f"row {row_number}: {name} must lie between 0 and shots")
        if elapsed < 0:
            raise ValueError(
                f"row {row_number}: elapsed_seconds must be nonnegative")
        if not decoder:
            raise ValueError(f"row {row_number}: decoder must be nonempty")

        aspect_ratio = steps // Ly
        key = (Ly, probability, aspect_ratio)
        if key in keys:
            raise ValueError(
                "CSV contains a duplicate size/probability/aspect row: "
                f"L={Ly}, p={probability:g}, T/L={aspect_ratio}")
        keys.add(key)
        rows.append({
            "steps": steps,
            "Ly": Ly,
            "Lz": Lz,
            "L": Ly,
            "aspect_ratio": aspect_ratio,
            "aspect_family": "T=L" if aspect_ratio == 1 else "T=2L",
            "p": probability,
            "shots": shots,
            "failures": failures,
            "logical_1_failures": logical_1,
            "logical_2_failures": logical_2,
            "seed": seed,
            "decoder": decoder,
            "elapsed_seconds": elapsed,
        })

    groups = {}
    for row in rows:
        groups.setdefault((row["L"], row["p"]), set()).add(
            row["aspect_ratio"])
    invalid = [key for key, aspects in groups.items() if aspects != {1, 2}]
    if invalid:
        L, probability = sorted(invalid)[0]
        raise ValueError(
            "each (L,p) must contain exactly T=L and T=2L; "
            f"incomplete pair at L={L}, p={probability:g}")

    sizes = sorted({row["L"] for row in rows})
    probability_grids = {
        L: {row["p"] for row in rows if row["L"] == L} for L in sizes
    }
    reference_grid = probability_grids[sizes[0]]
    if any(grid != reference_grid for grid in probability_grids.values()):
        raise ValueError("all sizes must use a common probability grid")

    return sorted(rows, key=lambda row: (
        row["L"], row["p"], row["aspect_ratio"]))


def point_estimates(rows):
    """Add the binomial rate and 95% Wilson interval to every input point."""
    output = []
    for row in rows:
        low, high = wilson_interval(row["failures"], row["shots"])
        point = dict(row)
        point.update({
            "failure_rate": row["failures"] / row["shots"],
            "wilson_low": low,
            "wilson_high": high,
            "interval_method": "two-sided 95% Wilson score interval",
        })
        output.append(point)
    return output


def aspect_differences(points):
    """Compute unpaired T=2L minus T=L differences and interval envelopes."""
    grouped = {}
    for point in points:
        grouped.setdefault((point["L"], point["p"]), {})[
            point["aspect_ratio"]] = point

    output = []
    for (L, probability), pair in sorted(grouped.items()):
        short = pair[1]
        long = pair[2]
        output.append({
            "L": L,
            "p": probability,
            "short_steps": short["steps"],
            "long_steps": long["steps"],
            "short_shots": short["shots"],
            "long_shots": long["shots"],
            "short_failures": short["failures"],
            "long_failures": long["failures"],
            "short_seed": short["seed"],
            "long_seed": long["seed"],
            "decoder": short["decoder"],
            "short_failure_rate": short["failure_rate"],
            "long_failure_rate": long["failure_rate"],
            "failure_rate_difference": (
                long["failure_rate"] - short["failure_rate"]),
            "short_wilson_low": short["wilson_low"],
            "short_wilson_high": short["wilson_high"],
            "long_wilson_low": long["wilson_low"],
            "long_wilson_high": long["wilson_high"],
            "conservative_low": (
                long["wilson_low"] - short["wilson_high"]),
            "conservative_high": (
                long["wilson_high"] - short["wilson_low"]),
            "difference_direction": "T=2L minus T=L",
            "interval_method": (
                "conservative envelope from separate 95% Wilson intervals; "
                "unpaired comparison"),
        })
    return output


def _write_csv(path, rows):
    if not rows:
        raise ValueError("cannot write an empty CSV")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def make_plot(points, differences, output_directory):
    """Create a deterministic, headless two-panel comparison figure."""
    os.environ.setdefault(
        "MPLCONFIGDIR", str(Path(tempfile.gettempdir())
                            / "toric3d-z-aspect-matplotlib"))
    import matplotlib

    matplotlib.use("Agg")
    matplotlib.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "svg.fonttype": "none",
        "svg.hashsalt": "toric3d-h6-z-aspect-ratio",
    })
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter

    sizes = sorted({point["L"] for point in points})
    colors = plt.get_cmap("viridis")(
        [0.12 + 0.76 * index / max(1, len(sizes) - 1)
         for index in range(len(sizes))])
    styles = {
        1: {"linestyle": "-", "marker": "o", "label": "T=L"},
        2: {"linestyle": "--", "marker": "s", "label": "T=2L"},
    }

    figure, (upper, lower) = plt.subplots(
        2, 1, figsize=(8.2, 7.2), sharex=True,
        gridspec_kw={"height_ratios": [1.6, 1.0]}, constrained_layout=True)

    for L, color in zip(sizes, colors):
        for aspect_ratio in (1, 2):
            selected = sorted(
                (point for point in points
                 if point["L"] == L
                 and point["aspect_ratio"] == aspect_ratio),
                key=lambda point: point["p"])
            x = [point["p"] for point in selected]
            y = [point["failure_rate"] for point in selected]
            yerr = [
                [point["failure_rate"] - point["wilson_low"]
                 for point in selected],
                [point["wilson_high"] - point["failure_rate"]
                 for point in selected],
            ]
            style = styles[aspect_ratio]
            upper.errorbar(
                x, y, yerr=yerr, color=color,
                linestyle=style["linestyle"], marker=style["marker"],
                markersize=4.5, linewidth=1.4, capsize=2.5,
                label=f"L={L}, {style['label']}")

        selected_differences = sorted(
            (item for item in differences if item["L"] == L),
            key=lambda item: item["p"])
        x = [item["p"] for item in selected_differences]
        y = [item["failure_rate_difference"]
             for item in selected_differences]
        yerr = [
            [item["failure_rate_difference"] - item["conservative_low"]
             for item in selected_differences],
            [item["conservative_high"] - item["failure_rate_difference"]
             for item in selected_differences],
        ]
        lower.errorbar(
            x, y, yerr=yerr, color=color, marker="o", markersize=4.5,
            linewidth=1.4, capsize=2.5, label=f"L={L}")

    upper.set_ylabel("Final logical-Z failure probability")
    upper.yaxis.set_major_formatter(PercentFormatter(xmax=1.0))
    upper.set_title("Fixed-size failure curves with 95% Wilson intervals")
    upper.grid(alpha=0.22)
    upper.legend(fontsize=7.6, ncol=2)

    lower.axhline(0.0, color="black", linewidth=1.0, linestyle=":")
    lower.set_xlabel("Z fault probability per active gate output, p")
    lower.set_ylabel("T=2L minus T=L")
    lower.yaxis.set_major_formatter(PercentFormatter(xmax=1.0))
    lower.set_title(
        "Unpaired difference; conservative envelope of separate intervals")
    lower.grid(alpha=0.22)
    lower.legend(fontsize=8, ncol=max(1, len(sizes)))

    figure.suptitle(
        "H+6-CNOT 3D toric-code Z faults: final-only decode, no feedback\n"
        "Aspect-ratio diagnostic only; no threshold fit",
        fontsize=11)

    fixed_date = datetime.datetime(
        2000, 1, 1, tzinfo=datetime.timezone.utc)
    figure.savefig(
        output_directory / "aspect-ratio-comparison.pdf",
        metadata={
            "Title": "H+6 final-only Z aspect-ratio diagnostic",
            "Author": "TopoNoise.jl",
            "Creator": "TopoNoise.jl",
            "CreationDate": fixed_date,
            "ModDate": fixed_date,
        })
    svg_path = output_directory / "aspect-ratio-comparison.svg"
    figure.savefig(
        svg_path,
        metadata={"Title": "H+6 final-only Z aspect-ratio diagnostic",
                  "Creator": "TopoNoise.jl", "Date": "2000-01-01"})
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_path.read_text().splitlines())
        + "\n")
    figure.savefig(
        output_directory / "aspect-ratio-comparison.png", dpi=200,
        metadata={"Software": "TopoNoise.jl"})
    plt.close(figure)


def main(arguments=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(arguments)

    source = args.csv.resolve()
    rows = load_rows(source)
    points = point_estimates(rows)
    differences = aspect_differences(points)

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    _write_csv(output / "aspect-ratio-points.csv", points)
    _write_csv(output / "aspect-ratio-differences.csv", differences)

    sizes = sorted({row["L"] for row in rows})
    rates = sorted({row["p"] for row in rows})
    summary = {
        "status": "diagnostic_complete",
        "analysis": "H+6-CNOT Z aspect-ratio comparison",
        "experiment": {
            "circuit": "H+6 CNOTs",
            "recovery": "final-only",
            "feedback": "none (no feedback)",
            "observable": "any residual transverse logical Z failure",
        },
        "geometry": {
            "comparison": "T=L versus T=2L at fixed Ly=Lz=L",
            "difference_direction": "T=2L minus T=L",
        },
        "sizes": sizes,
        "rates": rates,
        "point_count": len(points),
        "difference_count": len(differences),
        "uncertainty": {
            "points": "two-sided 95% Wilson score intervals",
            "differences": (
                "conservative interval (long Wilson low - short Wilson high, "
                "long Wilson high - short Wilson low)"),
            "paired_statistics": False,
        },
        "threshold": {
            "fitted": False,
            "estimate": None,
            "reason": (
                "This three-size aspect-ratio study is a fixed-size "
                "diagnostic; no finite-size threshold fit was requested or "
                "performed."),
        },
        "input": {
            "path": str(source),
            "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "columns": list(REQUIRED_COLUMNS),
        },
        "analysis_script": {
            "path": str(Path(__file__).resolve()),
            "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        "decoders": sorted({row["decoder"] for row in rows}),
        "seeds": sorted({row["seed"] for row in rows}),
        "outputs": {
            "points": "aspect-ratio-points.csv",
            "differences": "aspect-ratio-differences.csv",
            "figures": [
                "aspect-ratio-comparison.pdf",
                "aspect-ratio-comparison.svg",
                "aspect-ratio-comparison.png",
            ],
        },
    }
    (output / "aspect-ratio-summary.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    make_plot(points, differences, output)

    print(json.dumps({
        "status": summary["status"],
        "sizes": sizes,
        "rates": rates,
        "threshold_fitted": False,
        "output_directory": str(output),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
