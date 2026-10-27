#!/usr/bin/env python3
"""Plot the matched-preparation TN pilot without discarding unknown outcomes.

Usage: python examples/tnml_repaired_seed/plot.py --root RESULTS --output FIGURES

Only RESULTS/production is used by default. Diagnostics can be included explicitly
with --include-diagnostics; they remain separate series and are labelled as such.
No threshold fit is performed. This script does not run a decoder.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import time
import tomllib


REPO = Path(__file__).resolve().parents[2]
REFERENCE = REPO / "results/toric3d-h6-matched-2026-10-17/data/matched-results.csv"
SIZES = (3, 5, 7, 9, 11)
COLORS = {3: "#7F7F7F", 5: "#0072B2", 7: "#E69F00", 9: "#009E73", 11: "#CC79A7"}


def read_toml(path):
    with Path(path).open("rb") as stream:
        return tomllib.load(stream)


def integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def boolean(value, name):
    if not isinstance(value, bool):
        raise ValueError(f"{name} must be Boolean")
    return value


def wilson(failures, shots):
    """Pointwise 95% Wilson interval; exact counts, no pseudocount outcome."""
    if shots <= 0 or not 0 <= failures <= shots:
        raise ValueError("Wilson interval requires 0 <= failures <= shots and shots > 0")
    z = 1.959963984540054
    rate = failures / shots
    denominator = 1 + z * z / shots
    center = (rate + z * z / (2 * shots)) / denominator
    half = z * math.sqrt(rate * (1 - rate) / shots + z * z / (4 * shots * shots)) / denominator
    return max(0.0, center - half), min(1.0, center + half)


def clean_json(value):
    if isinstance(value, dict):
        return {str(key): clean_json(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [clean_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def load_run(directory, phase):
    directory = Path(directory)
    metadata = read_toml(directory / "metadata.toml")
    L = integer(metadata["L"], "L", 2)
    p = float(metadata["p"])
    shots = integer(metadata["shots"], "shots", 1)
    chi = integer(metadata["chi"], "chi", 1)
    if not math.isfinite(p) or not 0 <= p < .5:
        raise ValueError(f"invalid p in {directory}")
    if metadata.get("tn_uses_both_records") is not True:
        raise ValueError(f"this plot requires the seed-plus-final conditioning model: {directory}")
    run_id = metadata["run_id"]
    records = {}
    corpus_path = directory / "input-corpus.toml"
    if corpus_path.exists():
        corpus = read_toml(corpus_path)
        if corpus.get("run_id") != run_id or corpus.get("requested_shots") != shots:
            raise ValueError(f"corpus identity/count mismatch: {directory}")
        for record in corpus["records"]:
            index = integer(record["shot_index"], "shot_index", 1)
            if index > shots or index in records:
                raise ValueError(f"duplicate or out-of-range corpus shot: {directory}/{index}")
            boolean(record["mwpm_failure"], "mwpm_failure")
            if integer(record["true_mask"], "true_mask") > 3:
                raise ValueError(f"invalid two-logical-bit truth in {directory}")
            records[index] = record
        # The runner atomically saves all histories before contraction starts.
        if len(records) != shots:
            raise ValueError(f"partial input corpus (expected an atomic complete corpus): {directory}")

    known_failures = valid = invalid = exceptions = 0
    output_indices, elapsed, exact_linf = set(), [], []
    exact_compared = exact_disagreements = 0
    for path in sorted((directory / "shots").glob("shot_*.toml")):
        output = read_toml(path)
        index = integer(output["shot_index"], "shot_index", 1)
        record = records.get(index)
        if index in output_indices or record is None:
            raise ValueError(f"duplicate output or output without saved history: {path}")
        output_indices.add(index)
        if output.get("run_id") != run_id or output.get("record_id") != record["record_id"]:
            raise ValueError(f"checkpoint identity mismatch: {path}")
        for key, expected in (("L", L), ("p", p), ("chi", chi),
                              ("true_mask", record["true_mask"]),
                              ("mwpm_failure", record["mwpm_failure"])):
            if output.get(key) != expected:
                raise ValueError(f"checkpoint {key} mismatch: {path}")
        status = output["status"]
        is_valid = boolean(output["valid"], "valid")
        if status not in ("complete", "invalid", "exception") or is_valid != (status == "complete"):
            raise ValueError(f"inconsistent completion status: {path}")
        if is_valid:
            failure = boolean(output["TNfailure"], "TNfailure")
            mask = integer(output["logical_mask"], "logical_mask")
            if mask > 3 or failure != (record["true_mask"] != mask):
                raise ValueError(f"TN outcome does not match the physical logical mask: {path}")
            valid += 1
            known_failures += int(failure)
        elif status == "invalid":
            invalid += 1
        else:
            exceptions += 1
        duration = output.get("elapsed_seconds")
        if isinstance(duration, (int, float)) and math.isfinite(duration) and duration >= 0:
            elapsed.append(duration)
        comparison = output.get("exact_comparison", {})
        if is_valid and "ml_tie_aware_match" in comparison:
            exact_compared += 1
            exact_disagreements += not boolean(comparison["ml_tie_aware_match"], "ml_tie_aware_match")
        if math.isfinite(comparison.get("posterior_linf", math.nan)):
            exact_linf.append(comparison["posterior_linf"])

    supervision_path = directory / "supervision.json"
    supervision = json.loads(supervision_path.read_text()) if supervision_path.exists() else None
    unknown = shots - valid
    mwpm_failures = sum(record["mwpm_failure"] for record in records.values()) if records else None
    return dict(phase=phase, directory=str(directory.resolve()), run_id=run_id,
                L=L, p=p, chi=chi, rng_seed=metadata["rng_seed"], requested_shots=shots,
                saved_histories=len(records), saved_outputs=len(output_indices), valid_shots=valid,
                tn_failures=known_failures, unknown_shots=unknown, invalid_shots=invalid,
                exception_shots=exceptions, missing_outputs=shots-len(output_indices),
                tn_complete=unknown == 0, tn_rate=known_failures/shots if unknown == 0 else None,
                tn_lower_bound=known_failures/shots,
                tn_upper_bound=(known_failures+unknown)/shots,
                tn_wilson95=wilson(known_failures, shots) if unknown == 0 else None,
                mwpm_failures=mwpm_failures, mwpm_rate=mwpm_failures/shots if records else None,
                mwpm_wilson95=wilson(mwpm_failures, shots) if records else None,
                exact_compared=exact_compared, exact_ml_disagreements=exact_disagreements,
                exact_max_posterior_linf=max(exact_linf, default=None),
                completed_decode_seconds=sum(elapsed), median_decode_seconds=statistics.median(elapsed) if elapsed else None,
                metadata=metadata, supervision=supervision)


def load_reference(path):
    rows = []
    with Path(path).open(newline="") as stream:
        for raw in csv.DictReader(stream):
            if raw["channel"] != "z_only":
                continue
            L, shots, failures = int(raw["steps"]), int(raw["shots"]), int(raw["final_only_failures"])
            if L not in SIZES:
                continue
            if int(raw["Ly"]) != L or int(raw["Lz"]) != L or shots <= 0 or not 0 <= failures <= shots:
                raise ValueError("reference must have T = Ly = Lz and valid counts")
            rows.append(dict(L=L, p=float(raw["p"]), shots=shots, failures=failures,
                             rate=failures/shots, wilson95=wilson(failures, shots)))
    if len({(row["L"], row["p"]) for row in rows}) != len(rows):
        raise ValueError("duplicate archival reference points")
    return sorted(rows, key=lambda row: (row["L"], row["p"]))


def aggregate_runs(runs):
    """Combine independent shards at a fixed L,p,chi; never pool chi or phases."""
    groups = {}
    for run in runs:
        key = run["phase"], run["L"], run["p"], run["chi"]
        groups.setdefault(key, []).append(run)
    points = []
    for (phase, L, p, chi), members in sorted(groups.items()):
        if len({member["rng_seed"] for member in members}) != len(members):
            raise ValueError(f"repeated RNG seed at {phase}, L={L}, p={p}, chi={chi}; do not double-count histories")
        N = sum(member["requested_shots"] for member in members)
        F = sum(member["tn_failures"] for member in members)
        U = sum(member["unknown_shots"] for member in members)
        mwpm_complete = all(member["saved_histories"] == member["requested_shots"] for member in members)
        mwpm = sum(member["mwpm_failures"] for member in members) if mwpm_complete else None
        points.append(dict(phase=phase, L=L, p=p, chi=chi, requested_shots=N,
                           valid_shots=N-U, tn_failures=F, unknown_shots=U,
                           tn_rate=F/N if not U else None,
                           tn_wilson95=wilson(F, N) if not U else None,
                           tn_lower_bound=F/N, tn_upper_bound=(F+U)/N,
                           mwpm_failures=mwpm, mwpm_rate=mwpm/N if mwpm_complete else None,
                           mwpm_wilson95=wilson(mwpm, N) if mwpm_complete else None,
                           runs=[member["directory"] for member in members]))
    return points


def draw_errorbar(axis, x, rate, interval, *, color, marker="o", filled=True, **kwargs):
    return axis.errorbar(x, rate, yerr=[[max(0, rate-interval[0])], [max(0, interval[1]-rate)]],
                         color=color, marker=marker, markerfacecolor=color if filled else "white",
                         markersize=5, capsize=3, linewidth=1.25, linestyle="none", **kwargs)


def render(points, reference, output, include_diagnostics, synthetic_preview=False):
    os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/toponoise-repaired-seed-mpl")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    plt.rcParams.update({"font.family": "sans-serif", "font.size": 11, "axes.titlesize": 13,
                         "axes.labelsize": 12, "axes.linewidth": 1.1,
                         "xtick.direction": "out", "ytick.direction": "out",
                         "savefig.dpi": 200, "pdf.fonttype": 42, "axes.grid": False})
    fig, axes = plt.subplots(1, 2, figsize=(13.4, 6.7), sharex=True, sharey=True)
    fig.subplots_adjust(left=.075, right=.98, bottom=.30, top=.83, wspace=.13)
    fig.suptitle("No feedback after preparation", fontsize=21, y=.975)
    fig.text(.5, .918, r"Z errors · diagonal H+6CNOT circuit · $T=L$ · complete initial seed repair", ha="center", fontsize=12)
    axes[0].set_title("Uploaded-figure reference: final-syndrome MWPM", pad=12)
    title = "Fresh finite-bond TN pilot" if not include_diagnostics else "TN pilot + labelled diagnostic histories"
    axes[1].set_title("SYNTHETIC LAYOUT TEST — no TN observations" if synthetic_preview else title, pad=12)

    for L in SIZES:
        rows = [row for row in reference if row["L"] == L]
        axes[0].plot([100*r["p"] for r in rows], [r["rate"] for r in rows], color=COLORS[L], lw=1.2)
        for row in rows:
            draw_errorbar(axes[0], 100*row["p"], row["rate"], row["wilson95"], color=COLORS[L])
    size_handles = [Line2D([], [], color=COLORS[L], marker="o", lw=1.2, ms=4, label=f"L={L}") for L in SIZES]
    axes[0].legend(handles=size_handles, loc="upper left", frameon=False, fontsize=10, ncol=2)

    series = sorted({(point["phase"], point["L"], point["chi"]) for point in points})
    markers = ("o", "s", "D", "^", "v")
    chi_values = sorted({point["chi"] for point in points})
    for phase, L, chi in series:
        rows = sorted([r for r in points if (r["phase"], r["L"], r["chi"]) == (phase, L, chi)], key=lambda r:r["p"])
        color = COLORS.get(L, "#333333")
        marker = markers[chi_values.index(chi) % len(markers)]
        # NaN breaks the line at any requested but unfinished point.
        axes[1].plot([100*r["p"] for r in rows], [r["tn_rate"] if r["tn_rate"] is not None else math.nan for r in rows],
                     color=color, lw=1.3, linestyle=":" if phase == "diagnostics" else "-")
        for row in rows:
            if row["tn_rate"] is not None:
                draw_errorbar(axes[1], 100*row["p"], row["tn_rate"], row["tn_wilson95"], color=color, marker=marker)
            else:
                lo, hi = row["tn_lower_bound"], row["tn_upper_bound"]
                axes[1].vlines(100*row["p"], lo, hi, color=color, alpha=.42, linewidth=5)
                axes[1].plot([100*row["p"]]*2, [lo,hi], "_", color=color, ms=9)
        # This comparison is paired with the full saved corpus for each TN point.
        mwpm_rows = [row for row in rows if row["mwpm_rate"] is not None]
        axes[1].plot([100*r["p"] for r in mwpm_rows], [r["mwpm_rate"] for r in mwpm_rows], color=color, lw=1, ls="--", alpha=.8)
        for row in mwpm_rows:
            draw_errorbar(axes[1], 100*row["p"], row["mwpm_rate"], row["mwpm_wilson95"], color=color, marker=marker, filled=False, alpha=.8)

    series_handles = [Line2D([], [], color=COLORS.get(L, "#333333"), marker=markers[chi_values.index(chi) % len(markers)],
                            lw=1.2, label=rf"L={L}, $\chi$={chi}" + (" (diagnostic)" if phase == "diagnostics" else ""))
                      for phase, L, chi in series]
    if series_handles:
        axes[1].legend(handles=series_handles, frameon=False, loc="upper left", fontsize=10)
    else:
        axes[1].text(.5, .48, "No production TN histories available", ha="center", va="center", transform=axes[1].transAxes, color=".4")
    method_handles = [Line2D([], [], color=".3", marker="o", lw=1.3, label="TN: seed + final syndromes"),
                      Line2D([], [], color=".3", marker="o", mfc="white", ls="--", lw=1.1, label="Paired MWPM: final syndrome")]
    if any(row["unknown_shots"] for row in points):
        method_handles.append(Line2D([], [], color=".5", lw=5, alpha=.42, label="Unfinished TN: outcome bounds"))
    fig.legend(handles=method_handles, loc="lower center", bbox_to_anchor=(.53,.17), frameon=False, ncol=len(method_handles), fontsize=10)

    maximum_p = max([r["p"] for r in reference] + [r["p"] for r in points] + [.0425])
    for axis in axes:
        axis.set_xlim(0, 100*maximum_p*1.035)
        axis.set_ylim(0, 1)
        axis.set_xlabel("Gate-output error probability p (%)", labelpad=8)
        axis.tick_params(labelsize=11)
    axes[0].set_ylabel("Z channel\nEither-logical failure probability", labelpad=9)
    reference_N = sorted({r["shots"] for r in reference})
    Ns = sorted({r["requested_shots"] for r in points})
    ntext = str(Ns[0]) if len(Ns) == 1 else (f"{min(Ns)}–{max(Ns)}" if Ns else "0")
    completed_sizes = sorted({r["L"] for r in points if r["unknown_shots"] == 0})
    missing_sizes = [L for L in SIZES if L not in completed_sizes]
    available = ", ".join(map(str, completed_sizes)) or "none"
    missing = ", ".join(map(str, missing_sizes)) or "none"
    reference_ntext = str(reference_N[0]) if len(reference_N) == 1 else "/".join(map(str,reference_N))
    interval_note = "Bars: pointwise 95% Wilson intervals for complete outcomes."
    if any(row["unknown_shots"] for row in points):
        interval_note += " Unfinished TN bounds include every unknown history; they are not confidence intervals."
    elif points:
        interval_note += f" All {sum(row['requested_shots'] for row in points)} TN outcomes completed; no histories were discarded."
    lines = [f"Reference: {reference_ntext} shots/point. Fresh pilot: N={ntext} requested shots/point; exact counts in summary.csv.",
             "Same physical preparation. TN also uses the measured seed syndrome; this is not an equal-information decoder comparison.",
             interval_note,
             f"Requested sizes: 3, 5, 7, 9, 11. Sizes with complete TN points: {available}; without complete TN points: {missing}. No threshold fit."]
    fig.text(.075, .029, "\n".join(lines), fontsize=8.9, va="bottom", linespacing=1.55)
    paths = [output / "tnml-matched-figure.png", output / "tnml-matched-figure.pdf"]
    for path in paths:
        fig.savefig(path, facecolor="white")
    plt.close(fig)
    return [str(path.resolve()) for path in paths]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--reference", type=Path, default=REFERENCE)
    parser.add_argument("--include-diagnostics", action="store_true")
    parser.add_argument("--synthetic-preview", action="store_true", help="mark temporary layout fixtures explicitly as synthetic")
    args = parser.parse_args(argv)
    root = args.root.resolve()
    output = (args.output or root / "figures").resolve()
    output.mkdir(parents=True, exist_ok=True)
    phases = ["production"] + (["diagnostics"] if args.include_diagnostics else [])
    runs = [load_run(path.parent, phase) for phase in phases
            for path in sorted((root / phase).rglob("metadata.toml"))]
    if len({run["run_id"] for run in runs}) != len(runs):
        raise ValueError("duplicate run identities in selected directories")
    points = aggregate_runs(runs)
    reference = load_reference(args.reference)
    figures = render(points, reference, output, args.include_diagnostics, args.synthetic_preview)
    columns = ["phase", "L", "p", "chi", "requested_shots", "saved_histories", "saved_outputs",
               "valid_shots", "tn_failures", "unknown_shots", "invalid_shots", "exception_shots",
               "missing_outputs", "tn_complete", "tn_rate", "tn_lower_bound", "tn_upper_bound",
               "tn_wilson95_lower", "tn_wilson95_upper", "mwpm_failures", "mwpm_rate",
               "mwpm_wilson95_lower", "mwpm_wilson95_upper", "exact_compared", "exact_ml_disagreements",
               "exact_max_posterior_linf", "completed_decode_seconds", "median_decode_seconds",
               "rng_seed", "run_id", "directory"]
    with (output / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for run in runs:
            flat = dict(run)
            for method in ("tn", "mwpm"):
                interval = run[method + "_wilson95"]
                flat[method + "_wilson95_lower"] = interval[0] if interval else None
                flat[method + "_wilson95_upper"] = interval[1] if interval else None
            writer.writerow(flat)
    report = dict(schema_version=1, generated_at_unix=time.time(), root=str(root),
                  diagnostics_included=args.include_diagnostics, synthetic_preview=args.synthetic_preview, requested_sizes=list(SIZES),
                  threshold_fit=None, threshold_policy="No threshold inference from this pilot.",
                  information_difference="TN conditions on seed and final syndromes; MWPM uses the final syndrome only. Both use the same physical initial seed repair.",
                  invalid_policy="Every requested history remains in N. With F known failures and U unknown outcomes, TN is bounded by [F/N,(F+U)/N]. No rate is estimated by conditioning on valid contractions.",
                  intervals="Pointwise 95% Wilson for complete binary outcomes; unknown-outcome bounds are not confidence intervals.",
                  reference=dict(path=str(args.reference.resolve()), sha256=hashlib.sha256(args.reference.read_bytes()).hexdigest(), points=reference),
                  totals=dict(requested_shots=sum(r["requested_shots"] for r in runs),
                              valid_shots=sum(r["valid_shots"] for r in runs),
                              unknown_shots=sum(r["unknown_shots"] for r in runs)),
                  runs=runs, points=points, figures=figures,
                  plot_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (output / "report.json").write_text(json.dumps(clean_json(report), indent=2, allow_nan=False) + "\n")
    print(json.dumps({"runs": len(runs), "points": len(points), **report["totals"], "figures": figures}))
    return report


if __name__ == "__main__":
    main()
