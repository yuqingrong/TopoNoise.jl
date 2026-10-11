#!/usr/bin/env python3
"""Guarded finite-size analysis for the encoded diagonal X-feedback scan.

The input is the raw CSV written by ``save_toric3d_x_feedback_threshold``.
Only ``metric=state_failure`` rows are fitted. Final-only and slice-feedback
counts are fitted independently with a binomial finite-size scaling model.
Numerical candidates remain available for diagnosis, but ``pc`` and ``nu``
are emitted only when crossing, optimizer, monotonicity, and goodness-of-fit
guards all pass.
"""

import argparse
import copy
import csv
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import sys
import tempfile

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, xlogy
from scipy.stats import chi2


PROTOCOL_COLUMNS = {
    "final_only": "final_only_failures",
    "slice_feedback": "slice_feedback_failures",
}
SETTING_COLUMNS = (
    "max_iter", "bp_method", "ms_scaling_factor", "schedule",
    "omp_thread_count", "osd_method", "osd_order", "ldpc_version",
)


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def wilson_interval(failures, shots, z=1.959963984540054):
    if shots <= 0 or not 0 <= failures <= shots:
        raise ValueError("Wilson counts require 0 <= failures <= shots and shots > 0")
    probability = failures / shots
    denominator = 1 + z * z / shots
    center = (probability + z * z / (2 * shots)) / denominator
    half = z * math.sqrt(
        probability * (1 - probability) / shots + z * z / (4 * shots * shots)
    ) / denominator
    return max(0.0, center - half), min(1.0, center + half)


def _integer(row, key):
    try:
        value = float(row[key])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{key} must be an integer") from error
    if not math.isfinite(value) or not value.is_integer():
        raise ValueError(f"{key} must be an integer")
    return int(value)


def validate_points(rows):
    """Validate and normalize one rectangular cubic state-failure grid."""
    normalized = []
    keys = set()
    routings = set()
    settings = set()
    for source in rows:
        if str(source.get("metric", "state_failure")) != "state_failure":
            continue
        row = dict(source)
        steps = _integer(row, "steps")
        Ly = _integer(row, "Ly")
        Lz = _integer(row, "Lz")
        shots = _integer(row, "shots")
        final = _integer(row, "final_only_failures")
        feedback = _integer(row, "slice_feedback_failures")
        try:
            probability = float(row["p"])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError("p must be a finite probability") from error
        if not math.isfinite(probability) or not 0 <= probability < 0.5:
            raise ValueError("p must be finite and lie in [0, 0.5)")
        if steps < 2 or not steps == Ly == Lz:
            raise ValueError("finite-size scaling requires steps = Ly = Lz >= 2")
        if shots <= 0 or not 0 <= final <= shots or not 0 <= feedback <= shots:
            raise ValueError("failure counts must lie between zero and shots")
        key = (steps, probability)
        if key in keys:
            raise ValueError(f"duplicate state-failure point L={steps}, p={probability}")
        keys.add(key)
        routing = str(row.get("routing", "unspecified"))
        routings.add(routing)
        setting = tuple(str(row.get(name, "unspecified")) for name in SETTING_COLUMNS)
        settings.add(setting)
        row.update(steps=steps, Ly=Ly, Lz=Lz, p=probability, shots=shots,
                   final_only_failures=final,
                   slice_feedback_failures=feedback,
                   metric="state_failure", routing=routing)
        normalized.append(row)
    if not normalized:
        raise ValueError("the CSV contains no state_failure rows")
    if len(routings) != 1:
        raise ValueError("do not pool different feedback routings in one fit")
    if len(settings) != 1:
        raise ValueError("do not pool different decoder settings in one fit")
    sizes = sorted({row["steps"] for row in normalized})
    probabilities = sorted({row["p"] for row in normalized})
    expected = {(size, probability) for size in sizes for probability in probabilities}
    if keys != expected:
        missing = sorted(expected - keys)
        raise ValueError(f"scan grid is not rectangular; missing points: {missing[:5]}")
    return sorted(normalized, key=lambda row: (row["steps"], row["p"]))


def read_points(path):
    """Read Task 5 raw CSV data and retain only combined state-failure rows."""
    source = Path(path)
    with source.open(newline="") as stream:
        rows = [row for row in csv.DictReader(stream)
                if row.get("metric") == "state_failure"]
    return validate_points(rows)


def _selected(points, p_min, p_max, min_size):
    return [row for row in points
            if row["steps"] >= min_size and p_min <= row["p"] <= p_max]


def adjacent_crossings(points, protocol):
    """Return every sign-changing interpolated crossing of adjacent sizes."""
    if protocol not in PROTOCOL_COLUMNS:
        raise ValueError(f"unknown protocol {protocol}")
    column = PROTOCOL_COLUMNS[protocol]
    sizes = sorted({row["steps"] for row in points})
    output = []
    for small, large in zip(sizes[:-1], sizes[1:]):
        curves = []
        for size in (small, large):
            pairs = sorted((row["p"], row[column] / row["shots"])
                           for row in points if row["steps"] == size)
            if not pairs:
                curves = []
                break
            curves.append(tuple(np.asarray(values, dtype=float)
                                for values in zip(*pairs)))
        if len(curves) != 2:
            continue
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
            if left * right < 0:
                crossing = grid[index] - left * (grid[index + 1] - grid[index]) / (right - left)
                candidates.append(float(crossing))
        for index in range(1, len(grid) - 1):
            if (difference[index] == 0 and difference[index - 1] != 0
                    and difference[index + 1] != 0
                    and difference[index - 1] * difference[index + 1] < 0):
                candidates.append(float(grid[index]))
        for crossing in sorted(set(candidates)):
            output.append({"sizes": [small, large], "p": crossing})
    return output


def _candidate_fit(points, protocol, *, p_min, p_max, degree, nu_bounds,
                   goodness_alpha, correction=False, warm_start=None):
    column = PROTOCOL_COLUMNS[protocol]
    sizes = np.asarray([row["steps"] for row in points], dtype=float)
    probabilities = np.asarray([row["p"] for row in points], dtype=float)
    shots = np.asarray([row["shots"] for row in points], dtype=float)
    failures = np.asarray([row[column] for row in points], dtype=float)
    observed = failures / shots
    scale = p_max - p_min
    reference = float(np.exp(np.mean(np.log(np.unique(sizes)))))
    relative_log = np.log(sizes / reference)
    correction_basis = reference / sizes
    total_shots = float(np.sum(shots))

    bounds = [(0.0, 1.0), tuple(np.log(nu_bounds))]
    bounds.extend([(-15.0, 15.0), (-30.0, 30.0)])
    bounds.extend([(-30.0, 30.0)] * (degree - 1))
    if correction:
        bounds.append((-15.0, 15.0))

    def objective(parameters):
        center, log_nu = parameters[:2]
        nu = math.exp(log_nu)
        stretch = np.exp(relative_log / nu)
        x = ((probabilities - p_min) / scale - center) * stretch
        coefficients = parameters[2:2 + degree + 1]
        eta = np.polynomial.polynomial.polyval(x, coefficients)
        if correction:
            eta = eta + parameters[-1] * correction_basis
        derivative = np.polynomial.polynomial.polyval(
            x, np.arange(1, len(coefficients)) * coefficients[1:])
        residual = shots * expit(eta) - failures
        loss = np.sum(shots * np.logaddexp(0, eta) - failures * eta) / total_shots
        jacobian = [-stretch * derivative,
                    -x * relative_log / nu * derivative]
        jacobian.extend(x ** power for power in range(degree + 1))
        if correction:
            jacobian.append(correction_basis)
        gradient = np.asarray([np.dot(residual, value) / total_shots
                               for value in jacobian], dtype=float)
        return float(loss), gradient

    empirical = float(np.clip(np.sum(failures) / total_shots, 0.001, 0.999))
    intercept = math.log(empirical / (1 - empirical))
    if warm_start is not None:
        starts = [np.asarray(warm_start, dtype=float)]
    else:
        crossings = adjacent_crossings(points, protocol)
        estimate = (float(np.median([item["p"] for item in crossings]))
                    if crossings else (p_min + p_max) / 2)
        center = float(np.clip((estimate - p_min) / scale, 0.02, 0.98))
        correlation = np.corrcoef(probabilities, observed)[0, 1]
        slope = -2.0 if np.isfinite(correlation) and correlation < 0 else 2.0
        centers = sorted(set((0.25, 0.5, 0.75, center)))
        starts = []
        for candidate_center in centers:
            for exponent in (0.6, 1.2, 2.5):
                parameters = [candidate_center, math.log(exponent),
                              intercept, slope] + [0.0] * (degree - 1)
                correction and parameters.append(0.0)
                starts.append(np.asarray(parameters, dtype=float))
    optimizations = [minimize(
        objective, start, method="L-BFGS-B", jac=True, bounds=bounds,
        options={"maxiter": 2000, "ftol": 1e-13, "gtol": 1e-9},
    ) for start in starts]
    converged = [result for result in optimizations
                 if result.success and np.isfinite(result.fun)]
    best = min(converged or optimizations, key=lambda result: result.fun)
    center, log_nu = best.x[:2]
    pc = p_min + float(center) * scale
    nu = math.exp(float(log_nu))
    coefficients = np.asarray(best.x[2:2 + degree + 1], dtype=float)
    x = ((probabilities - pc) / scale) * (sizes / reference) ** (1 / nu)
    eta = np.polynomial.polynomial.polyval(x, coefficients)
    correction_coefficient = float(best.x[-1]) if correction else None
    if correction:
        eta = eta + correction_coefficient * correction_basis
    fitted = expit(eta)
    safe = np.clip(fitted, 1e-14, 1 - 1e-14)
    deviance = float(2 * np.sum(
        xlogy(failures, np.divide(observed, safe))
        + xlogy(shots - failures, np.divide(1 - observed, 1 - safe))))
    degrees = len(points) - len(best.x)
    goodness = float(chi2.sf(deviance, degrees)) if degrees > 0 else None
    names = ["pc", "nu"] + [f"a{power}" for power in range(degree + 1)]
    correction and names.append("correction_omega1")
    bound_hits = [name for name, value, (low, high) in zip(names, best.x, bounds)
                  if min(value - low, high - value) <= 0.002 * (high - low)]
    x_grid = np.linspace(float(np.min(x)), float(np.max(x)), 400)
    derivative = np.polynomial.polynomial.polyval(
        x_grid, np.arange(1, len(coefficients)) * coefficients[1:])
    monotone = bool(np.min(derivative) > 1e-7)
    issues = []
    if not best.success:
        issues.append("optimizer did not converge")
    if bound_hits:
        issues.append("parameter bound reached: " + ", ".join(bound_hits))
    if not monotone:
        issues.append("fitted failure probability is not monotone increasing in p")
    if goodness is None:
        issues.append("nonpositive goodness-of-fit degrees of freedom")
    elif goodness < goodness_alpha:
        issues.append(f"scaling-model lack of fit (deviance p < {goodness_alpha:g})")
    return {
        "pc": pc,
        "nu": nu,
        "coefficients": coefficients.tolist(),
        "polynomial_degree": degree,
        "p_scale": scale,
        "L_reference": reference,
        "p_min": p_min,
        "p_max": p_max,
        "nu_bounds": list(nu_bounds),
        "correction_omega": 1.0 if correction else None,
        "correction_coefficient": correction_coefficient,
        "deviance": deviance,
        "degrees_of_freedom": degrees,
        "deviance_p_value": goodness,
        "optimizer_success": bool(best.success),
        "optimizer_message": str(best.message),
        "bound_hits": bound_hits,
        "monotone_increasing": monotone,
        "issues": issues,
        "parameters": best.x.tolist(),
        "fitted_probabilities": fitted.tolist(),
    }


def fit_protocol(points, protocol, *, p_min=None, p_max=None, min_size=None,
                 degree=2, nu_bounds=(0.3, 5.0), goodness_alpha=0.01,
                 correction=False):
    """Fit one protocol and expose a threshold only after every guard passes."""
    if protocol not in PROTOCOL_COLUMNS:
        raise ValueError(f"protocol must be one of {tuple(PROTOCOL_COLUMNS)}")
    points = validate_points(points)
    all_probabilities = sorted({row["p"] for row in points})
    all_sizes = sorted({row["steps"] for row in points})
    p_min = all_probabilities[0] if p_min is None else float(p_min)
    p_max = all_probabilities[-1] if p_max is None else float(p_max)
    min_size = all_sizes[0] if min_size is None else int(min_size)
    if not 0 <= p_min < p_max < 0.5:
        raise ValueError("fit window must satisfy 0 <= p_min < p_max < 0.5")
    if degree not in (2, 3):
        raise ValueError("polynomial degree must be 2 or 3")
    if not 0 < nu_bounds[0] < nu_bounds[1]:
        raise ValueError("nu bounds must be positive and ordered")
    if not 0 < goodness_alpha < 1:
        raise ValueError("goodness alpha must lie in (0, 1)")
    selected = _selected(points, p_min, p_max, min_size)
    sizes = sorted({row["steps"] for row in selected})
    probabilities = sorted({row["p"] for row in selected})
    base = {
        "protocol": protocol,
        "status": "insufficient_data",
        "pc": None,
        "nu": None,
        "candidate": None,
        "diagnostic": [],
        "crossings": [],
        "selected_sizes": sizes,
        "selected_probabilities": probabilities,
        "selected_points": len(selected),
        "fit_window": [p_min, p_max],
        "min_size": min_size,
        "degree": degree,
        "nu_bounds": list(nu_bounds),
        "goodness_alpha": goodness_alpha,
        "correction_omega": 1.0 if correction else None,
    }
    if len(sizes) < 3 or len(probabilities) < 5:
        base["diagnostic"].append(
            "need at least three sizes and five probabilities in the fit window")
        return base
    parameter_count = degree + 3 + int(correction)
    if len(selected) <= parameter_count:
        base["diagnostic"].append("need more points than fitted parameters")
        return base
    crossings = adjacent_crossings(selected, protocol)
    base["crossings"] = crossings
    observed_pairs = {tuple(item["sizes"]) for item in crossings}
    expected_pairs = set(zip(sizes[:-1], sizes[1:]))
    if not crossings:
        base["status"] = "no_crossing"
        base["diagnostic"].append(
            "no adjacent-size sign-changing crossing inside the fit window")
        return base
    if observed_pairs != expected_pairs:
        base["status"] = "inconclusive"
        missing = sorted(expected_pairs - observed_pairs)
        base["diagnostic"].append(
            f"missing adjacent-size crossings for {missing}; threshold withheld")
        return base
    candidate = _candidate_fit(
        selected, protocol, p_min=p_min, p_max=p_max, degree=degree,
        nu_bounds=nu_bounds, goodness_alpha=goodness_alpha,
        correction=correction,
    )
    base["candidate"] = candidate
    base["diagnostic"].extend(candidate["issues"])
    if candidate["issues"]:
        base["status"] = "inconclusive"
        return base
    base.update(status="ok", pc=candidate["pc"], nu=candidate["nu"])
    return base


def bootstrap_protocol(points, protocol, fit, *, replicates=1000, seed=1234):
    """Parametric binomial bootstrap conditional on one accepted fit."""
    if replicates < 0:
        raise ValueError("bootstrap replicates must be nonnegative")
    if fit.get("status") != "ok" or fit.get("candidate") is None:
        return {
            "replicates": replicates, "accepted": 0, "pc_interval": None,
            "nu_interval": None, "samples": [],
            "method": "not run because the primary fit was not accepted",
        }
    points = validate_points(points)
    selected = _selected(points, *fit["fit_window"], fit["min_size"])
    candidate = fit["candidate"]
    probabilities = candidate["fitted_probabilities"]
    rng = np.random.default_rng(seed)
    samples = []
    correction = fit["correction_omega"] is not None
    for replicate in range(replicates):
        sampled = [dict(row) for row in selected]
        column = PROTOCOL_COLUMNS[protocol]
        for row, probability in zip(sampled, probabilities):
            row[column] = int(rng.binomial(row["shots"], probability))
        replica = _candidate_fit(
            sampled, protocol, p_min=fit["fit_window"][0],
            p_max=fit["fit_window"][1], degree=fit["degree"],
            nu_bounds=tuple(fit["nu_bounds"]),
            goodness_alpha=fit["goodness_alpha"], correction=correction,
            warm_start=candidate["parameters"],
        )
        accepted = (replica["optimizer_success"] and not replica["bound_hits"]
                    and replica["monotone_increasing"])
        samples.append({
            "replicate": replicate,
            "accepted": bool(accepted),
            "pc": replica["pc"] if accepted else None,
            "nu": replica["nu"] if accepted else None,
            "status": "accepted" if accepted else "rejected",
            "diagnostic": "; ".join(replica["issues"]),
        })
    accepted_samples = [sample for sample in samples if sample["accepted"]]
    required = max(2, math.ceil(0.8 * replicates)) if replicates else 0
    intervals_available = replicates > 0 and len(accepted_samples) >= required
    pc_interval = (np.quantile([sample["pc"] for sample in accepted_samples],
                               [0.025, 0.975]).tolist()
                   if intervals_available else None)
    nu_interval = (np.quantile([sample["nu"] for sample in accepted_samples],
                               [0.025, 0.975]).tolist()
                   if intervals_available else None)
    return {
        "replicates": replicates,
        "accepted": len(accepted_samples),
        "pc_interval": pc_interval,
        "nu_interval": nu_interval,
        "samples": samples,
        "method": "seeded parametric binomial bootstrap conditional on the accepted scaling model",
    }


def _crossing_drift_is_monotone(crossings):
    grouped = {}
    for crossing in crossings:
        grouped.setdefault(tuple(crossing["sizes"]), []).append(crossing["p"])
    sequence = [float(np.median(grouped[pair])) for pair in sorted(grouped)]
    if len(sequence) < 3:
        return False
    differences = np.diff(sequence)
    tolerance = 1e-12
    return bool((np.all(differences >= -tolerance)
                 or np.all(differences <= tolerance))
                and np.any(np.abs(differences) > tolerance))


def sensitivity_fits(points, protocol, primary):
    """Repeat guarded fits across sizes/windows and optional omega=1 correction."""
    sizes = primary["selected_sizes"]
    probabilities = primary["selected_probabilities"]
    low, high = primary["fit_window"]
    variants = []
    for minimum in sizes[1:-2]:
        variants.append((f"min_size_{minimum}", low, high, minimum, False))
    if len(probabilities) >= 7:
        variants.extend([
            ("trim_lowest_p", probabilities[1], high, primary["min_size"], False),
            ("trim_highest_p", low, probabilities[-2], primary["min_size"], False),
        ])
    if _crossing_drift_is_monotone(primary["crossings"]):
        variants.append(("fixed_omega_1_correction", low, high,
                         primary["min_size"], True))
    output = []
    for name, variant_low, variant_high, minimum, correction in variants:
        fit = fit_protocol(
            points, protocol, p_min=variant_low, p_max=variant_high,
            min_size=minimum, degree=primary["degree"],
            nu_bounds=tuple(primary["nu_bounds"]),
            goodness_alpha=primary["goodness_alpha"], correction=correction,
        )
        output.append({
            "protocol": protocol,
            "name": name,
            "fit_window": [variant_low, variant_high],
            "min_size": minimum,
            "correction_omega": 1.0 if correction else None,
            "status": fit["status"],
            "pc": fit["pc"],
            "nu": fit["nu"],
            "diagnostic": fit["diagnostic"],
            "candidate": fit["candidate"],
        })
    return output


def _finalize_bootstrap_guard(fit, bootstrap):
    output = copy.deepcopy(fit)
    output["bootstrap"] = bootstrap
    if output["status"] == "ok" and bootstrap["replicates"]:
        required = max(2, math.ceil(0.8 * bootstrap["replicates"]))
        if bootstrap["accepted"] < required:
            output["status"] = "inconclusive"
            output["pc"] = None
            output["nu"] = None
            output["diagnostic"].append(
                "fewer than 80% of parametric bootstrap fits were usable")
    return output


def _model_probability(candidate, p, size):
    x = ((np.asarray(p) - candidate["pc"]) / candidate["p_scale"]
         * (size / candidate["L_reference"]) ** (1 / candidate["nu"]))
    eta = np.polynomial.polynomial.polyval(x, candidate["coefficients"])
    if candidate["correction_omega"] is not None:
        eta = eta + candidate["correction_coefficient"] * candidate["L_reference"] / size
    return expit(eta)


def _make_plot(report, points, output, full_path_points=None):
    os.environ.setdefault(
        "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "toric3d-matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.2), constrained_layout=True)
    axes = axes.ravel()
    sizes = sorted({row["steps"] for row in points})
    colors = plt.get_cmap("viridis")(np.linspace(0.08, 0.92, len(sizes)))
    styles = {
        "final_only": ("o", "-", "Final-only"),
        "slice_feedback": ("s", "--", "Slice feedback"),
    }
    for size, color in zip(sizes, colors):
        selected = [row for row in points if row["steps"] == size]
        p = np.asarray([row["p"] for row in selected])
        for protocol, (marker, line, label) in styles.items():
            column = PROTOCOL_COLUMNS[protocol]
            rates = np.asarray([row[column] / row["shots"] for row in selected])
            intervals = np.asarray([wilson_interval(row[column], row["shots"])
                                    for row in selected]).T
            bars = np.maximum(0, np.vstack((rates - intervals[0], intervals[1] - rates)))
            axes[0].errorbar(p, rates, yerr=bars, marker=marker, ls=line,
                             color=color, ms=3.5, lw=1, capsize=1.5,
                             alpha=0.85, label=(f"L={size}, {label}"))
    axes[0].set(title="Raw state-failure curves", xlabel="X-fault probability p",
                ylabel="Logical |00> state-failure probability")
    axes[0].legend(fontsize=6, ncol=2)

    collapse_colors = {"final_only": "#2364AA", "slice_feedback": "#D1495B"}
    for protocol, fit in report["fits"].items():
        candidate = fit.get("candidate")
        if fit["status"] != "ok" or candidate is None:
            continue
        selected = _selected(points, *fit["fit_window"], fit["min_size"])
        column = PROTOCOL_COLUMNS[protocol]
        x = np.asarray([((row["p"] - candidate["pc"]) / candidate["p_scale"])
                        * (row["steps"] / candidate["L_reference"]) ** (1 / candidate["nu"])
                        for row in selected])
        y = np.asarray([row[column] / row["shots"] for row in selected])
        axes[1].scatter(x, y, s=13, alpha=0.6, color=collapse_colors[protocol],
                        label=protocol.replace("_", " "))
        grid = np.linspace(float(np.min(x)), float(np.max(x)), 300)
        curve = expit(np.polynomial.polynomial.polyval(grid, candidate["coefficients"]))
        axes[1].plot(grid, curve, color=collapse_colors[protocol], lw=2)
    axes[1].set(title="Finite-size scaling collapse",
                xlabel="x = ((p-pc)/sp)(L/Lref)^(1/nu)",
                ylabel="State-failure probability")
    axes[1].legend(fontsize=8)

    for protocol, fit in report["fits"].items():
        crossings = fit["crossings"]
        if not crossings:
            continue
        x = [(item["sizes"][0] + item["sizes"][1]) / 2 for item in crossings]
        y = [item["p"] for item in crossings]
        axes[2].plot(x, y, "o-", color=collapse_colors[protocol],
                     label=protocol.replace("_", " "))
        if fit["pc"] is not None:
            axes[2].axhline(fit["pc"], color=collapse_colors[protocol],
                            ls="--", lw=1)
    axes[2].set(title="Adjacent-size crossing drift", xlabel="Mean adjacent size",
                ylabel="Interpolated crossing p")
    axes[2].legend(fontsize=8)

    if full_path_points is not None:
        largest = min(max(row["steps"] for row in points),
                      max(row["steps"] for row in full_path_points))
        for source_points, label, line in (
                (points, "one_step", "-"), (full_path_points, "full_path", "--")):
            selected = [row for row in source_points if row["steps"] == largest]
            p = [row["p"] for row in selected]
            y = [row["slice_feedback_failures"] / row["shots"] for row in selected]
            axes[3].plot(p, y, marker="o", ls=line, label=label)
        axes[3].set(title=f"Routing sensitivity at L={largest}",
                    xlabel="X-fault probability p",
                    ylabel="Slice-feedback state failure")
    else:
        for protocol in PROTOCOL_COLUMNS:
            variants = [row for row in report["sensitivity"]
                        if row["protocol"] == protocol and row["pc"] is not None]
            axes[3].scatter(range(len(variants)), [row["pc"] for row in variants],
                            label=protocol.replace("_", " "))
        axes[3].set(title="Accepted model/window sensitivities",
                    xlabel="Accepted sensitivity index", ylabel="Candidate pc")
    axes[3].legend(fontsize=8)

    for axis in axes:
        axis.grid(alpha=0.2)
    fig.suptitle(
        "Encoded diagonal 3D toric-code X faults: final-only versus slice feedback\n"
        "Wilson bars; threshold values are withheld whenever a fit guard fails",
        fontsize=11,
    )
    for extension in ("pdf", "svg", "png"):
        fig.savefig(output / f"threshold-fit.{extension}", dpi=180)
    plt.close(fig)


def _write_bootstrap(path, fits):
    columns = ["protocol", "replicate", "accepted", "pc", "nu", "status", "diagnostic"]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for protocol, fit in fits.items():
            for sample in fit["bootstrap"]["samples"]:
                writer.writerow(dict(protocol=protocol, **sample))


def _write_sensitivity(path, rows):
    columns = ["protocol", "name", "p_min", "p_max", "min_size",
               "correction_omega", "status", "pc", "nu", "diagnostic"]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "protocol": row["protocol"], "name": row["name"],
                "p_min": row["fit_window"][0], "p_max": row["fit_window"][1],
                "min_size": row["min_size"],
                "correction_omega": row["correction_omega"],
                "status": row["status"], "pc": row["pc"], "nu": row["nu"],
                "diagnostic": "; ".join(row["diagnostic"]),
            })


def _write_points(path, points):
    columns = ["L", "p", "shots", "routing", "protocol", "failures",
               "rate", "wilson_low", "wilson_high"]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for row in points:
            for protocol, column in PROTOCOL_COLUMNS.items():
                failures = row[column]
                low, high = wilson_interval(failures, row["shots"])
                writer.writerow({
                    "L": row["steps"], "p": row["p"], "shots": row["shots"],
                    "routing": row["routing"], "protocol": protocol,
                    "failures": failures, "rate": failures / row["shots"],
                    "wilson_low": low, "wilson_high": high,
                })


def _package_versions():
    versions = {
        "python": platform.python_version(),
        "numpy": np.__version__,
    }
    for package in ("scipy", "matplotlib"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "unavailable"
    return versions


def main(arguments=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv", type=Path)
    parser.add_argument("--full-path-csv", type=Path)
    parser.add_argument("--p-min", type=float, required=True)
    parser.add_argument("--p-max", type=float, required=True)
    parser.add_argument("--min-size", type=int, default=5)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--degree", type=int, choices=(2, 3), default=2)
    parser.add_argument("--nu-min", type=float, default=0.3)
    parser.add_argument("--nu-max", type=float, default=5.0)
    parser.add_argument("--goodness-alpha", type=float, default=0.01)
    parser.add_argument("--no-sensitivity", action="store_true")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(arguments)
    if args.bootstrap < 0:
        parser.error("--bootstrap must be nonnegative")

    points = read_points(args.csv)
    protocol_seeds = np.random.SeedSequence(args.seed).spawn(len(PROTOCOL_COLUMNS))
    fits = {}
    sensitivity = []
    for (protocol, seed_sequence) in zip(PROTOCOL_COLUMNS, protocol_seeds):
        fit = fit_protocol(
            points, protocol, p_min=args.p_min, p_max=args.p_max,
            min_size=args.min_size, degree=args.degree,
            nu_bounds=(args.nu_min, args.nu_max),
            goodness_alpha=args.goodness_alpha,
        )
        bootstrap = bootstrap_protocol(
            points, protocol, fit, replicates=args.bootstrap,
            seed=int(seed_sequence.generate_state(1, dtype=np.uint64)[0]),
        )
        fits[protocol] = _finalize_bootstrap_guard(fit, bootstrap)
        if not args.no_sensitivity:
            sensitivity.extend(sensitivity_fits(points, protocol, fit))

    full_path_points = None
    full_path_fits = None
    if args.full_path_csv is not None:
        full_path_points = read_points(args.full_path_csv)
        full_path_fits = {
            protocol: fit_protocol(
                full_path_points, protocol, p_min=args.p_min, p_max=args.p_max,
                min_size=args.min_size, degree=args.degree,
                nu_bounds=(args.nu_min, args.nu_max),
                goodness_alpha=args.goodness_alpha,
            ) for protocol in PROTOCOL_COLUMNS
        }

    accepted_sensitivity = [row["pc"] for row in sensitivity if row["pc"] is not None]
    report = {
        "format_version": 1,
        "observable": "logical |00> state failure: either residual logical-X membrane",
        "protocols": list(PROTOCOL_COLUMNS),
        "fits": fits,
        "sensitivity": sensitivity,
        "rejected_sensitivity_fits": [
            {"protocol": row["protocol"], "name": row["name"],
             "status": row["status"], "diagnostic": row["diagnostic"]}
            for row in sensitivity if row["status"] != "ok"],
        "accepted_model_pc_spread": ([min(accepted_sensitivity), max(accepted_sensitivity)]
                                     if accepted_sensitivity else None),
        "full_path_fits": full_path_fits,
        "input_csv": str(args.csv.resolve()),
        "input_sha256": _sha256(args.csv),
        "full_path_input": (str(args.full_path_csv.resolve())
                            if args.full_path_csv is not None else None),
        "full_path_sha256": (_sha256(args.full_path_csv)
                             if args.full_path_csv is not None else None),
        "source_file": str(Path(__file__).resolve()),
        "source_sha256": _sha256(__file__),
        "command_arguments": list(arguments if arguments is not None else sys.argv[1:]),
        "package_versions": _package_versions(),
        "seed": args.seed,
        "bootstrap_replicates": args.bootstrap,
        "likelihood": "independent binomial counts for each protocol",
        "scaling_variable": "x=((p-pc)/sp)(L/Lref)^(1/nu)",
        "sampling_uncertainty": "parametric bootstrap conditional on each accepted scaling model",
        "model_uncertainty": "reported separately through size, window, and fixed-omega sensitivity fits",
    }

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    _write_bootstrap(output / "threshold-bootstrap.csv", fits)
    _write_sensitivity(output / "threshold-sensitivity.csv", sensitivity)
    _write_points(output / "threshold-points.csv", points)
    (output / "threshold-fit.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n")
    _make_plot(report, points, output, full_path_points)
    print(json.dumps({
        protocol: {"status": fit["status"], "pc": fit["pc"],
                   "nu": fit["nu"], "diagnostic": fit["diagnostic"]}
        for protocol, fit in fits.items()
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
