#!/usr/bin/env python3
"""Report the expanded TN scan without presenting pending runs as observations.

Only complete, validated production targets get empirical rates. The zero-noise
point is an analytical model prediction, independent of numerical regression
checks. The earlier pilot is always a separate dataset and panel.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time


_SPEC = importlib.util.spec_from_file_location("matched_pilot_plot", Path(__file__).with_name("plot.py"))
BASE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(BASE)
_DRIVER_SPEC = importlib.util.spec_from_file_location("full_grid_auditor", Path(__file__).with_name("run_full_grid.py"))
DRIVER = importlib.util.module_from_spec(_DRIVER_SPEC)
_DRIVER_SPEC.loader.exec_module(DRIVER)
DEFAULT_PILOT = BASE.REPO / "results/tnml-matched-figure-2026-10-27"
STATUS_COLORS = {"analytical": "#6a51a3", "complete": "#238b45", "running": "#3182bd",
                 "partial": "#9ecae1", "blocked": "#d95f0e", "pending": "#e5e5e5"}


def read_json(path, default=None):
    return json.loads(Path(path).read_text()) if Path(path).exists() else default


def key(L, p):
    return int(L), round(float(p), 12)


def load_phase(root, phase):
    runs = [BASE.load_run(path.parent, phase) for path in sorted((root / phase).rglob("metadata.toml"))]
    if len({run["run_id"] for run in runs}) != len(runs):
        raise ValueError(f"duplicate run identity in {root / phase}")
    return runs


def load_selection(root):
    selection = read_json(root / "scan-selection.json")
    if not isinstance(selection, dict):
        raise ValueError("scan-selection.json is required; the plot never invents the requested grid")
    sizes = selection["sizes"]
    rates = selection["p_values"]
    if (not sizes or len(set(sizes)) != len(sizes)
            or any(isinstance(L, bool) or not isinstance(L, int) or L < 2 for L in sizes)):
        raise ValueError("sizes must be distinct integers >= 2")
    if (not rates or rates[0] != 0 or sorted(rates) != rates or len(set(rates)) != len(rates)
            or any(not math.isfinite(p) or not 0 <= p < .5 for p in rates)):
        raise ValueError("p_values must be an increasing probability grid starting at zero")
    BASE.integer(selection["shots_per_nonzero_point"], "shots_per_nonzero_point", 1)
    return selection


def zero_regression_errors(run):
    """Check the deterministic distribution, not merely the tie-broken label."""
    errors, deviations = [], []
    directory = Path(run["directory"])
    corpus_path = directory / "input-corpus.toml"
    records = BASE.read_toml(corpus_path).get("records", []) if corpus_path.exists() else []
    if len(records) != run["requested_shots"]:
        errors.append("Zero-noise raw histories are incomplete.")
    for record in records:
        if record.get("raw_fault_indices") != [] or record.get("true_mask") != 0:
            errors.append(f"Shot {record['shot_index']}: the raw history is not the noiseless history.")
    for path in sorted((directory / "shots").glob("shot_*.toml")):
        shot = BASE.read_toml(path)
        probabilities = shot.get("probabilities", [])
        if (len(probabilities) != 4
                or any(not isinstance(v, (int, float)) or not math.isfinite(v) for v in probabilities)):
            errors.append(f"{path.name}: missing or nonfinite four-sector posterior.")
            continue
        deviation = max(abs(value - (i == 0)) for i, value in enumerate(probabilities))
        deviations.append(deviation)
        if deviation > 1e-8:
            errors.append(f"{path.name}: posterior differs from [1, 0, 0, 0].")
    return errors, max(deviations, default=None)


def build_report(root, pilot_root=None):
    root = Path(root).resolve()
    selection = load_selection(root)
    target = selection["shots_per_nonzero_point"]
    gates = read_json(root / "validation/gates.json", {}).get("gates", {})
    status_document = read_json(root / "scan-status.json", {})
    statuses = {}
    for item in status_document.get("points", {}).values():
        if item.get("phase") == "production":
            point_key = key(item["L"], item["p"])
            if point_key in statuses:
                raise ValueError("duplicate production point in scan-status.json")
            statuses[point_key] = item

    runs = load_phase(root, "production")
    expected_points = {key(point["L"], point["p"]): point
                       for point in DRIVER.selected_points(selection, "production")}
    allowed = {key(L, p) for L in selection["sizes"] for p in selection["p_values"] if p > 0}
    if any(key(run["L"], run["p"]) not in allowed for run in runs):
        raise ValueError("a production run is outside the frozen nonzero-p grid")
    for run in runs:
        point = expected_points[key(run["L"], run["p"])]
        if Path(run["directory"]) != (root / point["output"]).resolve():
            raise ValueError("production directory differs from the frozen selected point")
        DRIVER.validate_metadata_identity(run["metadata"], point, selection)
    gate_evidence = {}
    for L in selection["sizes"]:
        if gates.get(str(L), {}).get("status") == "passed":
            point = next(point for point in expected_points.values() if point["L"] == L)
            gate_evidence[str(L)] = DRIVER.require_gate(root, selection, point)
    aggregates = BASE.aggregate_runs(runs)
    points_by_key = {}
    for point in aggregates:
        point_key = key(point["L"], point["p"])
        if point_key in points_by_key:
            raise ValueError("different contraction dimensions cannot be pooled at a production point")
        gate = gates.get(str(point["L"]), {})
        if gate.get("chi") is not None and point["chi"] != gate["chi"]:
            raise ValueError("production chi differs from its recorded validation gate")
        if point["requested_shots"] > target:
            raise ValueError("production histories exceed the frozen target; do not silently pool extra shots")
        points_by_key[point_key] = point

    grid = []
    for L in selection["sizes"]:
        gate = gates.get(str(L), {})
        for p in selection["p_values"]:
            if p == 0:
                grid.append(dict(L=L, p=0.0, status="analytical", analytical_rate=0.0,
                                 requested_shots=0, valid_shots=0, unknown_shots=0,
                                 tn_failures=None, tn_rate=None, tn_wilson95=None,
                                 tn_lower_bound=None, tn_upper_bound=None,
                                 reason="Exact noiseless model: ideal measurements and corrections, all fault probabilities zero. No statistical trials or interval.",
                                 runs=[]))
                continue
            item = points_by_key.get(key(L, p))
            source_runs = [run for run in runs if key(run["L"], run["p"]) == key(L, p)]
            valid = item["valid_shots"] if item else 0
            failures = item["tn_failures"] if item else 0
            allocated = item["requested_shots"] if item else 0
            unknown = target - valid
            state = statuses.get(key(L, p), {})
            size_block = status_document.get("size_blocks", {}).get(f"production:L{L}")
            attempts = state.get("attempts", [])
            latest_attempt = attempts[-1] if attempts else {}
            supervisions = [run["supervision"] for run in source_runs if run["supervision"]]
            if state.get("output"):
                state_supervision = read_json(root / state["output"] / "supervision.json")
                if state_supervision and state_supervision not in supervisions:
                    supervisions.append(state_supervision)
            terminal_failures = [summary.get("status", "unknown") for summary in supervisions
                                 if summary.get("status") != "completed"
                                 and not (state.get("status") == "running"
                                          and not latest_attempt.get("finished_at_utc")
                                          and latest_attempt.get("previous_supervision") == summary)]
            exact_count = sum(run["exact_compared"] for run in source_runs)
            exact_disagreements = sum(run["exact_ml_disagreements"] for run in source_runs)
            requires_exact = L in selection.get("exact_reference_sizes", [])
            audit = DRIVER.audit_point(root, expected_points[key(L, p)], selection) if unknown == 0 else None
            status, reason = "pending", "No production observations yet."
            if gate.get("status") != "passed":
                status, reason = "blocked", gate.get("reason", "Contraction validation gate has not passed.")
            elif state.get("status") == "blocked":
                status, reason = "blocked", state.get("reason", "Driver reported a blocked point.")
            elif exact_disagreements:
                status, reason = "blocked", "TN decisions disagree with the exact reference."
            elif terminal_failures:
                status, reason = "blocked", "Terminal supervisor result: " + ", ".join(terminal_failures)
            elif audit is not None and not audit["passed"]:
                status, reason = "blocked", "Independent saved-file audit: " + "; ".join(audit["errors"])
            elif audit is not None and audit["passed"] and (not requires_exact or exact_count == target):
                status, reason = "complete", "All requested histories produced valid, validated outcomes."
            elif size_block:
                status, reason = "blocked", "Same-size production stopped: " + size_block.get("reason", "validation or runtime failure")
            elif state.get("status") == "running":
                if any(summary.get("status") == "completed" for summary in supervisions):
                    status, reason = "partial", "Supervised worker finished; target or audit remains incomplete."
                else:
                    status, reason = "running", "Driver records production as running; no partial failure-rate estimate."
            elif allocated:
                status, reason = "partial", "Requested target or exact-reference validation is unfinished; no partial failure-rate estimate."
            elif state.get("status") == "interrupted":
                status, reason = "partial", state.get("reason", "Production was interrupted before an observation completed.")
            complete = status == "complete"
            grid.append(dict(L=L, p=p, status=status, reason=reason, analytical_rate=None,
                             chi=item["chi"] if item else gate.get("chi"), requested_shots=target,
                             allocated_shots=allocated, valid_shots=valid, unknown_shots=unknown,
                             tn_failures=failures, tn_rate=failures/target if complete else None,
                             tn_wilson95=BASE.wilson(failures, target) if complete else None,
                             tn_lower_bound=failures/target, tn_upper_bound=(failures+unknown)/target,
                             exact_compared=exact_count, exact_ml_disagreements=exact_disagreements,
                             validation_audit=audit,
                             runs=item["runs"] if item else []))

    # These are numerical regression checks of a deterministic prediction, never
    # additional successes in the statistical nonzero-p experiment.
    zero_runs = load_phase(root, "zero")
    if any(run["p"] != 0 or run["L"] not in selection["sizes"] for run in zero_runs):
        raise ValueError("zero/ contains a run outside the requested zero-noise regression checks")
    zero_checks = []
    expected_zeros = {point["L"]: point for point in DRIVER.selected_points(selection, "zero")}
    for run in zero_runs:
        point = expected_zeros[run["L"]]
        if Path(run["directory"]) != (root / point["output"]).resolve():
            raise ValueError("zero-check directory differs from the frozen selected point")
        DRIVER.validate_metadata_identity(run["metadata"], point, selection)
        check = dict(run)
        check["validation_audit"] = DRIVER.audit_point(root, point, selection)
        check["regression_errors"], check["posterior_error_from_analytical"] = zero_regression_errors(run)
        check["regression_passed"] = (run["tn_complete"] and run["tn_failures"] == 0
                                      and run["mwpm_failures"] == 0 and run["exact_ml_disagreements"] == 0
                                      and not check["regression_errors"] and check["validation_audit"]["passed"])
        for field in ("tn_rate", "tn_wilson95", "mwpm_rate", "mwpm_wilson95", "tn_lower_bound", "tn_upper_bound"):
            check[field] = None
        zero_checks.append(check)

    pilot = dict(root=str(Path(pilot_root).resolve()) if pilot_root else None, runs=[], points=[],
                 totals=dict(requested_shots=0, valid_shots=0, unknown_shots=0))
    if pilot_root:
        pilot["runs"] = load_phase(Path(pilot_root), "production")
        pilot["points"] = BASE.aggregate_runs(pilot["runs"])
        pilot["totals"] = {field: sum(run[field] for run in pilot["runs"])
                           for field in ("requested_shots", "valid_shots", "unknown_shots")}
    positives = [point for point in grid if point["p"] > 0]
    return dict(schema_version=1, generated_at_unix=time.time(), root=str(root),
                selection=selection, selection_sha256=hashlib.sha256((root/"scan-selection.json").read_bytes()).hexdigest(),
                gates=gates, gate_evidence=gate_evidence, driver_status=status_document, grid=grid, runs=runs,
                zero_checks=zero_checks, pilot=pilot,
                audit_source_sha256={name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                     for name in ("run_full_grid.py", "plot.py")},
                frozen_scientific_sources_match_current=DRIVER.scientific_source_hashes() == selection["scientific_source_sha256"],
                totals=dict(planned_production_shots=target*len(positives),
                            allocated_production_shots=sum(run["requested_shots"] for run in runs),
                            valid_production_shots=sum(point["valid_shots"] for point in positives),
                            unknown_production_shots=sum(point["unknown_shots"] for point in positives),
                            complete_production_points=sum(point["status"] == "complete" for point in positives),
                            planned_production_points=len(positives), zero_regression_checks=len(zero_checks)),
                status_counts=dict(Counter(point["status"] for point in grid)),
                diagnostics_included=False, threshold_fit=None,
                zero_policy="p=0 is the exact analytical model value, not 1000 observed successes. Numerical zero checks are separate regressions with no binomial interval.",
                incomplete_policy="No rate until the entire requested target is complete and validated. F known failures and U unknown histories give [F/N,(F+U)/N] over the full target. These are outcome bounds, not confidence intervals.",
                pilot_policy="Earlier pilot observations are shown separately and never pooled into the newly requested target.",
                intervals="Pointwise 95% Wilson intervals only for complete validated nonzero-p targets.",
                model="Z-only gate-output faults, diagonal H+6CNOT preparation, initial seed MWPM repair, no later preparation feedback; TN uses both seed and final syndromes.")


def render(report, output, synthetic_preview=False):
    os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/toponoise-full-grid-mpl")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch, Rectangle

    plt.rcParams.update({"font.size": 11, "font.family": "sans-serif", "pdf.fonttype": 42,
                         "axes.titlesize": 13, "axes.labelsize": 12, "savefig.dpi": 180})
    selection = report["selection"]
    target = selection["shots_per_nonzero_point"]
    sizes, rates = selection["sizes"], selection["p_values"]
    fig = plt.figure(figsize=(15.5, 10.0))
    grid_spec = fig.add_gridspec(2, 2, left=.07, right=.975, top=.84, bottom=.23,
                               height_ratios=[1.65, 1], hspace=.54, wspace=.19)
    prior_ax = fig.add_subplot(grid_spec[0, 0])
    new_ax = fig.add_subplot(grid_spec[0, 1], sharex=prior_ax, sharey=prior_ax)
    coverage = fig.add_subplot(grid_spec[1, :])
    fig.suptitle("No feedback after preparation", fontsize=22, y=.97)
    fig.text(.5, .923, "Z-only errors · diagonal H+6CNOT circuit · initial seed repair · TN uses seed and final syndromes",
             ha="center", fontsize=11.5)
    fig.text(.5, .884, "SYNTHETIC REPORT TEST — no scientific observations" if synthetic_preview else
             f"Expanded scan: L = {', '.join(map(str, sizes))}; {target:,} new shots per nonzero-p point",
             ha="center", fontsize=12, color="#a33" if synthetic_preview else ".2")
    prior_ax.set_title("Earlier pilot — separate data", pad=10)
    new_ax.set_title(f"New scan — complete {target:,}-shot points only", pad=10)

    pilot_points = report["pilot"]["points"]
    for L in sorted({point["L"] for point in pilot_points}):
        rows = sorted([point for point in pilot_points if point["L"] == L], key=lambda p: p["p"])
        complete = [point for point in rows if point["tn_rate"] is not None]
        if not complete:
            continue
        counts = sorted({point["requested_shots"] for point in complete})
        counts_text = "/".join(map(str, counts))
        chi_text = "/".join(map(str, sorted({point["chi"] for point in complete})))
        color = BASE.COLORS.get(L, ".3")
        prior_ax.plot([100*p["p"] for p in rows],
                      [p["tn_rate"] if p["tn_rate"] is not None else math.nan for p in rows],
                      color=color, lw=1.1, label=rf"L={L}, $\chi$={chi_text}, {counts_text} shots/point")
        for point in complete:
            BASE.draw_errorbar(prior_ax, 100*point["p"], point["tn_rate"], point["tn_wilson95"], color=color)
    if prior_ax.get_legend_handles_labels()[0]:
        prior_ax.legend(frameon=False, loc="upper left", fontsize=10)
    else:
        prior_ax.text(.5, .5, "No earlier pilot selected", transform=prior_ax.transAxes, ha="center", color=".4")

    # A shared zero star denotes all selected sizes at the same exact coordinate;
    # visually offsetting it would falsely suggest nonzero probabilities.
    new_ax.plot(0, 0, marker="*", markersize=13, color="black", linestyle="none", clip_on=False,
                label="p=0: exact model, all sizes (no trials / CI)")
    for L in sizes:
        rows = [p for p in report["grid"] if p["L"] == L and p["p"] > 0]
        complete = [p for p in rows if p["tn_rate"] is not None]
        if not complete:
            continue
        color = BASE.COLORS.get(L, ".3")
        new_ax.plot([100*p["p"] for p in rows],
                    [p["tn_rate"] if p["tn_rate"] is not None else math.nan for p in rows],
                    color=color, lw=1.2, label=f"L={L}")
        for point in complete:
            BASE.draw_errorbar(new_ax, 100*point["p"], point["tn_rate"], point["tn_wilson95"], color=color)
    if not report["totals"]["complete_production_points"]:
        new_ax.text(.5, .49, f"No complete {target:,}-shot points yet\nSee coverage below for current progress",
                    ha="center", va="center", transform=new_ax.transAxes, color=".4", linespacing=1.6)
    new_ax.legend(frameon=False, loc="upper left", fontsize=9.5)
    for axis in (prior_ax, new_ax):
        axis.set_xlim(0, 100*max(rates)*1.035)
        axis.set_ylim(0, 1)
        axis.set_xlabel("Gate-output error probability p (%)", labelpad=7)
        axis.spines[["right", "top"]].set_visible(False)
    prior_ax.set_ylabel("Either-logical failure probability", labelpad=8)

    lookup = {key(point["L"], point["p"]): point for point in report["grid"]}
    for i, L in enumerate(sizes):
        for j, p in enumerate(rates):
            point = lookup[key(L, p)]
            color = STATUS_COLORS[point["status"]]
            coverage.add_patch(Rectangle((j-.5, i-.5), 1, 1, facecolor=color, edgecolor="white", lw=2))
            if p == 0:
                label = "0*"
            elif point["valid_shots"]:
                label = f"{point['valid_shots']}"
            else:
                label = "—"
            coverage.text(j, i, label, ha="center", va="center", fontsize=9,
                          color="white" if point["status"] in ("analytical", "complete", "running", "blocked") else ".3")
    coverage.set_xlim(-.5, len(rates)-.5)
    coverage.set_ylim(len(sizes)-.5, -.5)
    coverage.set_xticks(range(len(rates)), [f"{100*p:g}" for p in rates], rotation=55, ha="right", fontsize=9)
    coverage.set_yticks(range(len(sizes)), [f"L={L}" for L in sizes])
    coverage.set_xlabel("Requested gate-output error probability p (%)", labelpad=6)
    coverage.set_title(f"Coverage: numbers = completed valid TN outcomes / {target:,}; 0* = analytical zero, not shot count", fontsize=12, pad=10)
    coverage.tick_params(length=0)
    for spine in coverage.spines.values():
        spine.set_visible(False)
    legends = [Patch(facecolor=STATUS_COLORS[state], label=label) for state, label in (
        ("analytical", "Analytical p=0"), ("complete", "Complete + validated"),
        ("running", "Running"), ("partial", "Partial / audit pending"),
        ("blocked", "Blocked by validation / runtime"), ("pending", "Pending"))]
    fig.legend(handles=legends, loc="lower center", bbox_to_anchor=(.52, .085), frameon=False, ncol=3, fontsize=9.5)
    zero_passed = len({check["L"] for check in report["zero_checks"] if check["regression_passed"]})
    totals = report["totals"]
    fig.text(.07, .021,
             f"New scan: {totals['valid_production_shots']:,} valid outcomes of {totals['planned_production_shots']:,} requested; "
             f"{totals['complete_production_points']}/{totals['planned_production_points']} points complete. "
             f"Numerical p=0 checks passed: {zero_passed}/{len(sizes)} (separate regressions).\n"
             "Bars: pointwise 95% Wilson intervals. Incomplete points have no rate; outcome bounds and blocking reasons are in report.json.\n"
             f"Prior pilot is not pooled into the new {target:,}-shot targets. No threshold fitted; no uncomputed size curves are inferred.",
             fontsize=8.7, va="bottom", linespacing=1.5)
    paths = [output / "tnml-full-grid.png", output / "tnml-full-grid.pdf"]
    for path in paths:
        fig.savefig(path, facecolor="white")
    plt.close(fig)
    return [str(path.resolve()) for path in paths]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--pilot-root", type=Path, default=DEFAULT_PILOT)
    parser.add_argument("--no-pilot", action="store_true")
    parser.add_argument("--synthetic-preview", action="store_true")
    args = parser.parse_args(argv)
    output = (args.output or args.root / "figures").resolve()
    output.mkdir(parents=True, exist_ok=True)
    report = build_report(args.root, None if args.no_pilot else args.pilot_root)
    report["synthetic_preview"] = args.synthetic_preview
    report["figures"] = render(report, output, args.synthetic_preview)
    report["plot_script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (output / "report.json").write_text(json.dumps(BASE.clean_json(report), indent=2, allow_nan=False) + "\n")
    print(json.dumps({**report["totals"], "figures": report["figures"]}))
    return report


if __name__ == "__main__":
    main()
