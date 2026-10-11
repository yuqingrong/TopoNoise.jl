#!/usr/bin/env python3
"""Binomial finite-size analysis of the slice-transfer Z experiment.

The input counts failures of either transverse logical Z after final data-error
MWPM. The reported uncertainty is Monte Carlo uncertainty conditional on the
scaling model; the separate sensitivity fits probe finite-size/model effects.
No fit is attempted unless adjacent sizes cross inside the selected window.

Example:
  python examples/fit_toric3d_z_threshold.py scan.csv --p-min .01 --p-max .03 \
      --min-size 5 --bootstrap 200 --seed 1234 --output-dir results/analysis
"""

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import tempfile

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, xlogy
from scipy.stats import chi2


def wilson_interval(failures, shots, z=1.959963984540054):
    """Two-sided 95% Wilson score interval, including boundary counts."""
    if shots <= 0 or not 0 <= failures <= shots:
        raise ValueError("Wilson counts require 0 <= failures <= shots and shots > 0")
    probability = failures / shots
    denominator = 1 + z * z / shots
    center = (probability + z * z / (2 * shots)) / denominator
    half = z * math.sqrt(probability * (1 - probability) / shots + z * z / (4 * shots * shots)) / denominator
    return max(0.0, center - half), min(1.0, center + half)


def _validated_rows(rows):
    combined = {}
    decoders = set()
    for source in rows:
        row = dict(source)
        for key in ("steps", "Ly", "Lz", "shots", "failures"):
            value = float(row[key])
            if not math.isfinite(value) or not value.is_integer():
                raise ValueError(f"{key} must be an integer")
            row[key] = int(value)
        row["p"] = float(row["p"])
        if not math.isfinite(row["p"]) or not 0 <= row["p"] <= 1:
            raise ValueError("p must be finite and in [0, 1]")
        if row["steps"] < 2 or row["steps"] != row["Ly"] or row["steps"] != row["Lz"]:
            raise ValueError("This scaling analysis requires steps = Ly = Lz >= 2")
        if row["shots"] <= 0 or not 0 <= row["failures"] <= row["shots"]:
            raise ValueError("Counts require 0 <= failures <= shots and shots > 0")
        if "decoder" in row:
            decoders.add(row["decoder"])
        key = (row["steps"], row["p"])
        if key in combined:
            combined[key]["shots"] += row["shots"]
            combined[key]["failures"] += row["failures"]
        else:
            combined[key] = row
    if not combined:
        raise ValueError("The scan contains no rows")
    if len(decoders) > 1:
        raise ValueError("Do not pool different decoder settings in one scaling fit")
    return [combined[key] for key in sorted(combined)]


def _arrays(rows):
    return tuple(np.array([row[key] for row in rows], dtype=float)
                 for key in ("steps", "p", "shots", "failures"))


def adjacent_crossings(rows):
    """Positive-direction, linearly interpolated crossings of adjacent sizes.

    A crossing needs a sign reversal: equal zero-error or saturated curves do
    not provide a threshold. Curves are interpolated only over their overlap.
    Multiple crossings are preserved as a diagnostic, not cherry-picked.
    """
    sizes = sorted({row["steps"] for row in rows})
    output = []
    for small, large in zip(sizes[:-1], sizes[1:]):
        curves = []
        for size in (small, large):
            points = sorted((row["p"], row["failures"] / row["shots"])
                            for row in rows if row["steps"] == size)
            curves.append(tuple(np.array(values) for values in zip(*points)))
        low = max(curves[0][0][0], curves[1][0][0])
        high = min(curves[0][0][-1], curves[1][0][-1])
        grid = np.unique(np.concatenate([curve[0] for curve in curves]))
        grid = grid[(grid >= low) & (grid <= high)]
        if len(grid) < 3:
            continue
        difference = np.interp(grid, *curves[1]) - np.interp(grid, *curves[0])
        candidates = []
        for index in range(len(grid) - 1):
            left, right = difference[index:index + 2]
            if left < 0 < right:
                crossing = grid[index] - left * (grid[index + 1] - grid[index]) / (right - left)
                candidates.append(float(crossing))
        for index in range(1, len(grid) - 1):
            if difference[index] == 0 and difference[index - 1] < 0 < difference[index + 1]:
                candidates.append(float(grid[index]))
        for crossing in sorted(set(candidates)):
            output.append(dict(sizes=[small, large], p=crossing))
    return output


def _fit_model(rows, p_min, p_max, degree, nu_bounds, warm_start=None):
    sizes, probabilities, shots, failures = _arrays(rows)
    scale = p_max - p_min
    reference = float(np.exp(np.mean(np.log(np.unique(sizes)))))
    relative_log = np.log(sizes / reference)
    total_shots = float(np.sum(shots))
    bounds = [(0.0, 1.0), tuple(np.log(nu_bounds)), (-15.0, 15.0), (0.0001, 30.0)]
    bounds += [(-30.0, 30.0)] * (degree - 1)

    def objective(parameters):
        center, log_nu = parameters[:2]
        nu = math.exp(log_nu)
        stretch = np.exp(relative_log / nu)
        x = ((probabilities - p_min) / scale - center) * stretch
        coefficients = parameters[2:]
        eta = np.polynomial.polynomial.polyval(x, coefficients)
        derivative = np.polynomial.polynomial.polyval(x, np.arange(1, len(coefficients)) * coefficients[1:])
        residual = shots * expit(eta) - failures
        loss = np.sum(shots * np.logaddexp(0, eta) - failures * eta) / total_shots
        jacobian = [-stretch * derivative, -x * relative_log / nu * derivative]
        jacobian.extend(x ** power for power in range(degree + 1))
        gradient = np.array([np.dot(residual, value) / total_shots for value in jacobian])
        return float(loss), gradient

    empirical = np.clip(np.sum(failures) / total_shots, 0.001, 0.999)
    intercept = math.log(empirical / (1 - empirical))
    if warm_start is not None:
        starts = [np.array(warm_start, dtype=float)]
    else:
        crossings = adjacent_crossings(rows)
        estimate = np.median([item["p"] for item in crossings]) if crossings else (p_min + p_max) / 2
        centers = sorted(set([0.3, 0.7, float(np.clip((estimate - p_min) / scale, 0.05, 0.95))]))
        starts = [np.array([center, math.log(nu), intercept, 2.0] + [0.0] * (degree - 1))
                  for center in centers for nu in (0.7, 1.5, 3.0)]
    optimizations = [minimize(objective, start, method="L-BFGS-B", jac=True, bounds=bounds,
                              options=dict(maxiter=1500, ftol=1e-13, gtol=1e-9)) for start in starts]
    successful = [fit for fit in optimizations if fit.success and np.isfinite(fit.fun)]
    best = min(successful or optimizations, key=lambda fit: fit.fun)
    pc = p_min + float(best.x[0]) * scale
    nu = math.exp(float(best.x[1]))
    coefficients = best.x[2:]
    x = ((probabilities - pc) / scale) * (sizes / reference) ** (1 / nu)
    fitted = expit(np.polynomial.polynomial.polyval(x, coefficients))
    safe_fitted = np.clip(fitted, 1e-14, 1 - 1e-14)
    observed = failures / shots
    deviance = float(2 * np.sum(xlogy(failures, observed / safe_fitted)
                                + xlogy(shots - failures, (1 - observed) / (1 - safe_fitted))))
    degrees = len(rows) - len(best.x)
    issues = []
    names = ["pc", "nu"] + [f"a{power}" for power in range(degree + 1)]
    bound_hits = [name for name, value, (low, high) in zip(names, best.x, bounds)
                  if min(value - low, high - value) <= 0.001 * (high - low)]
    if not best.success:
        issues.append("optimizer did not converge")
    if bound_hits:
        issues.append("parameter bound reached: " + ", ".join(bound_hits))
    grid = np.linspace(float(np.min(x)), float(np.max(x)), 200)
    derivative = np.polynomial.polynomial.polyval(grid, np.arange(1, len(coefficients)) * coefficients[1:])
    if np.min(derivative) <= 0:
        issues.append("fitted failure probability is not monotone over the fit window")
    goodness = float(chi2.sf(deviance, degrees)) if degrees > 0 else None
    if goodness is not None and goodness < 0.01:
        issues.append("scaling-model lack of fit (deviance p < 0.01)")
    return dict(pc=pc, nu=nu, coefficients=coefficients.tolist(), polynomial_degree=degree,
                p_scale=scale, L_reference=reference, p_min=p_min, p_max=p_max,
                nu_bounds=list(nu_bounds), deviance=deviance, degrees_of_freedom=degrees,
                deviance_p_value=goodness, optimizer_success=bool(best.success),
                optimizer_message=str(best.message), bound_hits=bound_hits, issues=issues,
                status="ok" if not issues else "inconclusive", parameters=best.x.tolist(),
                fitted_probabilities=fitted.tolist())


def _crossing_bootstrap(rows, crossings, count, rng):
    intervals = [[] for _ in crossings]
    for _ in range(count):
        sample = [dict(row, failures=int(rng.binomial(row["shots"], row["failures"] / row["shots"]))) for row in rows]
        found = adjacent_crossings(sample)
        for index, base in enumerate(crossings):
            matching = [item["p"] for item in found if item["sizes"] == base["sizes"]]
            if matching:
                intervals[index].append(min(matching, key=lambda value: abs(value - base["p"])))
    for crossing, values in zip(crossings, intervals):
        crossing["bootstrap_crossing_fraction"] = len(values) / count if count else None
        crossing["interval"] = np.quantile(values, [0.025, 0.975]).tolist() if len(values) >= max(2, 0.8 * count) else None
        crossing["bootstrap_method"] = "independent binomial sampling at observed point rates; nearest same-pair crossing"


def analyze(rows, *, p_min, p_max, min_size=2, bootstrap=200, seed=1234,
            degree=2, nu_bounds=(0.3, 5.0), sensitivity=True):
    """Return serializable estimates and all diagnostics; never force a crossing."""
    if not 0 <= p_min < p_max <= 1:
        raise ValueError("The fit window must satisfy 0 <= p_min < p_max <= 1")
    if degree not in (2, 3) or bootstrap < 0 or not 0 < nu_bounds[0] < nu_bounds[1]:
        raise ValueError("Use degree 2 or 3, nonnegative bootstrap count, and positive ordered nu bounds")
    all_rows = _validated_rows(rows)
    selected = [row for row in all_rows if row["steps"] >= min_size and p_min <= row["p"] <= p_max]
    sizes = sorted({row["steps"] for row in selected})
    result = dict(observable="any residual transverse logical Z", status="no_crossing", fit=None,
                  pc_interval=None, confidence_level=0.95, fit_window=[p_min, p_max], min_size=min_size,
                  selected_sizes=sizes, selected_points=len(selected), bootstrap_count=bootstrap,
                  seed=seed, crossings=[], bootstrap_samples=[], sensitivity=[],
                  diagnostic=[], all_rows=all_rows,
                  uncertainty_scope="Monte Carlo uncertainty conditional on scaling model; not finite-size systematic error",
                  likelihood="independent binomial counts with logistic polynomial finite-size scaling",
                  scaling_variable="x = ((p - pc) / p_scale) * (L / L_reference) ** (1 / nu)")
    if len(sizes) < 2 or len(selected) <= degree + 3:
        result["status"] = "insufficient_data"
        result["diagnostic"].append("Need at least two sizes and more points than fitted parameters")
        return result
    crossings = adjacent_crossings(selected)
    result["crossings"] = crossings
    if not crossings:
        result["diagnostic"].append("No adjacent-size sign-changing crossing inside the fit window; no threshold fitted")
        return result
    rng_fit, rng_crossing = [np.random.default_rng(sequence) for sequence in np.random.SeedSequence(seed).spawn(2)]
    _crossing_bootstrap(selected, crossings, bootstrap, rng_crossing)
    fit = _fit_model(selected, p_min, p_max, degree, nu_bounds)
    result["fit"] = fit
    result["status"] = fit["status"]
    result["diagnostic"].extend(fit["issues"])
    if len({tuple(item["sizes"]) for item in crossings}) < len(sizes) - 1:
        result["diagnostic"].append("Not every adjacent-size pair crosses inside this fit window")
    for iteration in range(bootstrap):
        sample = [dict(row, failures=int(rng_fit.binomial(row["shots"], probability)))
                  for row, probability in zip(selected, fit["fitted_probabilities"])]
        replica = _fit_model(sample, p_min, p_max, degree, nu_bounds, warm_start=fit["parameters"])
        # A goodness-of-fit tail in a bootstrap replica is not a failed optimizer;
        # exclude nonconvergence, nonmonotonic curves, and boundary solutions only.
        valid = replica["optimizer_success"] and not replica["bound_hits"] and not any("not monotone" in issue for issue in replica["issues"])
        result["bootstrap_samples"].append(dict(replicate=iteration, pc=replica["pc"], nu=replica["nu"],
                                                accepted=bool(valid), status=replica["status"]))
    accepted = [item["pc"] for item in result["bootstrap_samples"] if item["accepted"]]
    result["bootstrap_accepted"] = len(accepted)
    if bootstrap and len(accepted) >= max(2, 0.8 * bootstrap):
        result["pc_interval"] = np.quantile(accepted, [0.025, 0.975]).tolist()
    elif bootstrap:
        result["status"] = "inconclusive"
        result["diagnostic"].append("Fewer than 80% usable bootstrap fits; threshold uncertainty is not reliable")
    if sensitivity:
        variants = []
        if len(sizes) >= 3:
            variants.append(("exclude_smallest_size", p_min, p_max, sizes[1], degree))
        levels = sorted({row["p"] for row in selected})
        if len(levels) >= 5:
            variants.append(("remove_lowest_p", levels[1], p_max, min_size, degree))
            variants.append(("remove_highest_p", p_min, levels[-2], min_size, degree))
        variants.append(("polynomial_degree_3" if degree == 2 else "polynomial_degree_2", p_min, p_max, min_size, 5 - degree))
        for name, low, high, minimum, order in variants:
            other = analyze(all_rows, p_min=low, p_max=high, min_size=minimum, bootstrap=0,
                            seed=seed, degree=order, nu_bounds=nu_bounds, sensitivity=False)
            result["sensitivity"].append(dict(name=name, status=other["status"], fit_window=[low, high],
                                               min_size=minimum, fit=other["fit"], diagnostic=other["diagnostic"]))
    return result


def _make_plot(result, output):
    # Defer plotting imports so the statistical API works without a display.
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "toric3d-matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8), constrained_layout=True)
    fit = result["fit"]
    sizes = sorted({row["steps"] for row in result["all_rows"]})
    colors = plt.get_cmap("viridis")(np.linspace(0.1, 0.9, len(sizes)))
    for size, color in zip(sizes, colors):
        rows = [row for row in result["all_rows"] if row["steps"] == size]
        probabilities = np.array([row["p"] for row in rows])
        rates = np.array([row["failures"] / row["shots"] for row in rows])
        intervals = np.array([wilson_interval(row["failures"], row["shots"]) for row in rows]).T
        bars = np.maximum(0, np.vstack([rates - intervals[0], intervals[1] - rates]))
        for axis in axes[:2]:
            axis.errorbar(probabilities, rates, yerr=bars, fmt="o", ms=4, capsize=2, color=color, label=f"L = {size}")
        if fit and size in result["selected_sizes"]:
            grid = np.linspace(*result["fit_window"], 200)
            x_grid = ((grid - fit["pc"]) / fit["p_scale"]) * (size / fit["L_reference"]) ** (1 / fit["nu"])
            curve = expit(np.polynomial.polynomial.polyval(x_grid, fit["coefficients"]))
            for axis in axes[:2]:
                axis.plot(grid, curve, color=color, lw=1.2)
            keep = (probabilities >= result["fit_window"][0]) & (probabilities <= result["fit_window"][1])
            x = ((probabilities[keep] - fit["pc"]) / fit["p_scale"]) * (size / fit["L_reference"]) ** (1 / fit["nu"])
            axes[2].errorbar(x, rates[keep], yerr=bars[:, keep], fmt="o", ms=4, capsize=2, color=color)
            axes[2].plot(x_grid, curve, color="black", alpha=0.2, lw=1)
    axes[0].set(xlabel="Z fault probability per active gate output, p", ylabel="Probability of either residual logical Z")
    axes[0].axvspan(*result["fit_window"], color="gray", alpha=0.08, label="fit window")
    axes[1].set(xlabel="Z fault probability per active gate output, p", ylabel="Probability of either residual logical Z")
    axes[0].set_title("All simulated probabilities", fontsize=11)
    axes[1].set_title("Crossing region", fontsize=11)
    axes[2].set_title("Finite-size scaling collapse", fontsize=11)
    full_ymax = min(1, max(0.05, max(row["failures"] / row["shots"] for row in result["all_rows"]) * 1.15))
    for axis in axes:
        axis.grid(alpha=0.2)
        axis.set_ylim(0, full_ymax)
    if fit:
        interval = result["pc_interval"]
        for axis in axes[:2]:
            if interval:
                axis.axvspan(*interval, color="#e9ad36", alpha=0.30, zorder=0,
                             label="95% bootstrap interval")
            axis.axvline(fit["pc"], color="black", ls="--", lw=1.2, label="candidate p_c")
        critical_rate = float(expit(fit["coefficients"][0]))
        width = result["fit_window"][1] - result["fit_window"][0]
        levels = sorted({row["p"] for row in result["all_rows"]
                         if row["steps"] in result["selected_sizes"]
                         and result["fit_window"][0] <= row["p"] <= result["fit_window"][1]})
        spacing = float(np.median(np.diff(levels))) if len(levels) > 1 else width / 4
        half_width = max(0.3 * width, 1.25 * spacing)
        if interval:
            half_width = max(half_width, 1.5 * max(abs(value - fit["pc"]) for value in interval))
        axes[1].set_xlim(max(result["fit_window"][0], fit["pc"] - half_width),
                         min(result["fit_window"][1], fit["pc"] + half_width))
        axes[1].set_ylim(0, min(1, max(0.005, 3 * critical_rate)))
        axes[1].plot(fit["pc"], critical_rate, "D", color="black", ms=4, zorder=5)
        axes[1].ticklabel_format(axis="x", style="plain", useOffset=False)
        from matplotlib.ticker import MaxNLocator
        axes[1].xaxis.set_major_locator(MaxNLocator(nbins=5))
        # A dedicated key keeps a narrow interval visible and identifiable even
        # when the full-range panel compresses the low-error crossing.
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch
        key = [Line2D([], [], color="black", ls="--", lw=1.2, label="candidate p_c")]
        if interval:
            key.append(Patch(facecolor="#e9ad36", alpha=0.5, label="95% bootstrap interval"))
        axes[1].legend(handles=key, fontsize=8, loc="upper left")
        subtitle = f"candidate p_c = {fit['pc']:.5g}, nu = {fit['nu']:.3g}; fit status: {result['status']}"
        if interval:
            subtitle += f"; 95% bootstrap interval [{interval[0]:.5g}, {interval[1]:.5g}]"
    else:
        axes[1].set_xlim(*result["fit_window"])
        axes[2].text(0.5, 0.5, "No supported crossing fit", transform=axes[2].transAxes, ha="center")
        subtitle = "No supported crossing fit"
    axes[0].legend(fontsize=8, ncol=2)
    axes[2].set(xlabel="Scaled distance from candidate threshold, x", ylabel="Probability of either residual logical Z")
    fig.suptitle("Slice-transfer Z noise: final MWPM and finite-size scaling\n" + subtitle
                 + "\nWilson 95% bars; fit uncertainty excludes finite-size systematic effects", fontsize=10.5)
    for suffix in ("pdf", "svg", "png"):
        fig.savefig(output / f"threshold-fit.{suffix}", dpi=180)
    plt.close(fig)


def main(arguments=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv", type=Path)
    parser.add_argument("--p-min", type=float, required=True)
    parser.add_argument("--p-max", type=float, required=True)
    parser.add_argument("--min-size", type=int, default=2)
    parser.add_argument("--bootstrap", type=int, default=200)
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--degree", type=int, choices=(2, 3), default=2)
    parser.add_argument("--nu-min", type=float, default=0.3)
    parser.add_argument("--nu-max", type=float, default=5.0)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, action="append", default=[])
    parser.add_argument("--no-sensitivity", action="store_true")
    args = parser.parse_args(arguments)
    with args.csv.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    result = analyze(rows, p_min=args.p_min, p_max=args.p_max, min_size=args.min_size,
                     bootstrap=args.bootstrap, seed=args.seed, degree=args.degree,
                     nu_bounds=(args.nu_min, args.nu_max), sensitivity=not args.no_sensitivity)
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    result["input_csv"] = str(args.csv.resolve())
    result["input_sha256"] = hashlib.sha256(args.csv.read_bytes()).hexdigest()
    result["metadata_files"] = []
    metadata = set(args.metadata) | set(args.csv.parent.glob("*.toml")) | set(args.csv.parent.glob("*metadata*.json"))
    for source in sorted(metadata):
        target = output / "metadata" / source.name
        target.parent.mkdir(exist_ok=True)
        if source.resolve() != target.resolve():
            shutil.copy2(source, target)
        result["metadata_files"].append(dict(source=str(source.resolve()), copied=str(target),
                                              sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
    with (output / "threshold-bootstrap.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["replicate", "pc", "nu", "accepted", "status"])
        writer.writeheader()
        writer.writerows(result["bootstrap_samples"])
    with (output / "threshold-points.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["L", "p", "shots", "failures", "rate", "wilson_low", "wilson_high", "used_in_fit"])
        writer.writeheader()
        for row in result["all_rows"]:
            low, high = wilson_interval(row["failures"], row["shots"])
            writer.writerow(dict(L=row["steps"], p=row["p"], shots=row["shots"], failures=row["failures"],
                                 rate=row["failures"] / row["shots"], wilson_low=low, wilson_high=high,
                                 used_in_fit=row["steps"] >= args.min_size and args.p_min <= row["p"] <= args.p_max))
    (output / "threshold-fit.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    _make_plot(result, output)
    compact = {key: result[key] for key in ("status", "pc_interval", "diagnostic")}
    if result["fit"]:
        compact.update(pc=result["fit"]["pc"], nu=result["fit"]["nu"], deviance=result["fit"]["deviance"],
                       degrees_of_freedom=result["fit"]["degrees_of_freedom"])
    print(json.dumps(compact, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
