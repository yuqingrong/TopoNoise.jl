#!/usr/bin/env python3
"""Guarded threshold analysis for diagonal 3D toric-code feedback strengths.

The input is one or more disjoint raw CSV shards written by
``save_toric3d_diagonal_feedback``. Only ``metric=state_failure`` rows are
fitted. Each channel/protocol marginal uses its own binomial likelihood. The
eight joint outcome counts are retained for paired rate-difference plots; the
six threshold estimates remain correlated because they replay the same faults.

Numerical candidates are retained for diagnosis, but ``pc`` and ``nu`` are
published only after data, crossing, optimizer, parameter-bound,
monotonicity, fit-window, goodness-of-fit, and bootstrap-usability guards pass.
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


PROTOCOL_BITS = {"none": 1, "slice_local": 2, "full_prefix": 4}
JOINT_COLUMNS = (
    "joint_n0_s0_f0", "joint_n1_s0_f0", "joint_n0_s1_f0",
    "joint_n1_s1_f0", "joint_n0_s0_f1", "joint_n1_s0_f1",
    "joint_n0_s1_f1", "joint_n1_s1_f1",
)
PAIR_SPECS = {
    "none_vs_slice_local": ("none", "slice_local"),
    "none_vs_full_prefix": ("none", "full_prefix"),
    "slice_local_vs_full_prefix": ("slice_local", "full_prefix"),
}
SETTING_COLUMNS = (
    "max_iter", "bp_method", "ms_scaling_factor", "schedule",
    "omp_thread_count", "osd_method", "osd_order",
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
        probability * (1 - probability) / shots
        + z * z / (4 * shots * shots)
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


def _pair_estimate(joint, first, second):
    first_bit, second_bit = PROTOCOL_BITS[first], PROTOCOL_BITS[second]
    neither = first_only = second_only = both = 0
    for mask, count in enumerate(joint):
        a, b = bool(mask & first_bit), bool(mask & second_bit)
        if a and b:
            both += count
        elif a:
            first_only += count
        elif b:
            second_only += count
        else:
            neither += count
    shots = sum(joint)
    difference = (second_only - first_only) / shots
    standard_error = (math.sqrt(max(
        0.0, ((second_only + first_only) / shots - difference * difference)
        / (shots - 1))) if shots > 1 else None)
    second_interval = wilson_interval(second_only, shots, z=2.241402727604947)
    first_interval = wilson_interval(first_only, shots, z=2.241402727604947)
    interval = (max(-1.0, second_interval[0] - first_interval[1]),
                min(1.0, second_interval[1] - first_interval[0]))
    return {
        "first": first, "second": second, "neither": neither,
        "first_only": first_only, "second_only": second_only, "both": both,
        "difference": difference, "standard_error": standard_error,
        "ci": list(interval),
    }


def _normalize_feedback_row(source):
    row = dict(source)
    steps = _integer(row, "steps")
    Ly = _integer(row, "Ly")
    Lz = _integer(row, "Lz")
    shots = _integer(row, "shots")
    try:
        probability = float(row["p"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("p must be a finite probability") from error
    if not math.isfinite(probability) or not 0 <= probability < 0.5:
        raise ValueError("p must be finite and lie in [0, 0.5)")
    if steps < 2 or not steps == Ly == Lz:
        raise ValueError("finite-size scaling requires steps = Ly = Lz >= 2")
    if shots <= 0:
        raise ValueError("shots must be positive")
    channel = str(row.get("channel", ""))
    if channel not in ("x_only", "z_only"):
        raise ValueError("channel must be x_only or z_only")
    joint = tuple(_integer(row, column) for column in JOINT_COLUMNS)
    if any(count < 0 for count in joint):
        raise ValueError("joint counts must be nonnegative")
    if sum(joint) != shots:
        raise ValueError("joint counts must sum to shots")
    marginals = {
        protocol: sum(count for mask, count in enumerate(joint) if mask & bit)
        for protocol, bit in PROTOCOL_BITS.items()
    }
    for protocol, failures in marginals.items():
        column = f"{protocol}_failures"
        if column in row and str(row[column]).strip() != "":
            declared = _integer(row, column)
            if declared != failures:
                raise ValueError(
                    f"{protocol} marginal does not match the joint counts")
        row[column] = failures
    paired = {
        name: _pair_estimate(joint, *protocols)
        for name, protocols in PAIR_SPECS.items()
    }
    row.update(
        steps=steps, Ly=Ly, Lz=Lz, p=probability, shots=shots,
        channel=channel, metric="state_failure", joint_counts=joint,
        paired_differences=paired,
        x_routing=str(row.get("x_routing", "unspecified")),
    )
    return row


def validate_feedback_points(rows):
    """Validate state-failure rows and rectangular grids independently by channel."""
    normalized = []
    keys = set()
    routings = set()
    settings = set()
    for source in rows:
        if str(source.get("metric", "state_failure")) != "state_failure":
            continue
        row = _normalize_feedback_row(source)
        key = (row["steps"], row["Ly"], row["Lz"], row["p"], row["channel"],
               row["metric"])
        if key in keys:
            raise ValueError(f"duplicate feedback point {key}")
        keys.add(key)
        routings.add(row["x_routing"])
        settings.add(tuple(str(row.get(name, "unspecified"))
                           for name in SETTING_COLUMNS))
        normalized.append(row)
    if not normalized:
        raise ValueError("the CSV contains no state_failure rows")
    if len(routings) != 1:
        raise ValueError("do not pool different X routing policies in one analysis")
    if len(settings) != 1:
        raise ValueError("do not pool different decoder settings in one analysis")
    for channel in sorted({row["channel"] for row in normalized}):
        selected = [row for row in normalized if row["channel"] == channel]
        dimensions = sorted({(row["steps"], row["Ly"], row["Lz"])
                             for row in selected})
        probabilities = sorted({row["p"] for row in selected})
        observed = {((row["steps"], row["Ly"], row["Lz"]), row["p"])
                    for row in selected}
        expected = {(dimensions_value, probability)
                    for dimensions_value in dimensions
                    for probability in probabilities}
        if observed != expected:
            missing = sorted(expected - observed)
            raise ValueError(
                f"{channel} scan grid is not rectangular; missing points: {missing[:5]}")
    return sorted(normalized,
                  key=lambda row: (row["channel"], row["steps"], row["p"]))


def read_feedback_points(paths):
    """Read disjoint CSV shards, reject duplicate row keys, and derive marginals."""
    if isinstance(paths, (str, os.PathLike, Path)):
        paths = [paths]
    sources = [Path(path) for path in paths]
    if not sources:
        raise ValueError("at least one raw CSV path is required")
    state_rows = []
    all_keys = set()
    for source in sources:
        with source.open(newline="") as stream:
            for row in csv.DictReader(stream):
                metric = str(row.get("metric", "state_failure"))
                try:
                    key = (_integer(row, "steps"), _integer(row, "Ly"),
                           _integer(row, "Lz"), float(row["p"]),
                           str(row.get("channel", "")), metric)
                except (KeyError, TypeError, ValueError) as error:
                    raise ValueError(f"invalid row identity in {source}: {error}") from error
                if key in all_keys:
                    raise ValueError(f"duplicate feedback point {key}")
                all_keys.add(key)
                if metric == "state_failure":
                    state_rows.append(row)
    return validate_feedback_points(state_rows)


def _selected(points, channel, p_min, p_max, min_size):
    return [row for row in points
            if row["channel"] == channel and row["steps"] >= min_size
            and p_min <= row["p"] <= p_max]


def adjacent_crossings(points, channel, protocol):
    """Return sign-changing interpolated crossings for adjacent cubic sizes."""
    if protocol not in PROTOCOL_BITS:
        raise ValueError(f"unknown protocol {protocol}")
    selected = [row for row in points if row["channel"] == channel]
    sizes = sorted({row["steps"] for row in selected})
    column = f"{protocol}_failures"
    output = []
    for small, large in zip(sizes[:-1], sizes[1:]):
        curves = []
        for size in (small, large):
            pairs = sorted((row["p"], row[column] / row["shots"])
                           for row in selected if row["steps"] == size)
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
                crossing = (grid[index]
                            - left * (grid[index + 1] - grid[index])
                            / (right - left))
                candidates.append(float(crossing))
        for index in range(1, len(grid) - 1):
            if (difference[index] == 0 and difference[index - 1] != 0
                    and difference[index + 1] != 0
                    and difference[index - 1] * difference[index + 1] < 0):
                candidates.append(float(grid[index]))
        for crossing in sorted(set(candidates)):
            output.append({"sizes": [small, large], "p": crossing})
    return output


def _minimum_polynomial_derivative(coefficients, x_low, x_high):
    derivative = np.polynomial.polynomial.polyder(coefficients)
    candidates = [float(x_low), float(x_high)]
    second = np.polynomial.polynomial.polyder(derivative)
    if len(second) > 1 and np.any(np.abs(second[1:]) > 0):
        roots = np.roots(np.asarray(second[::-1], dtype=float))
        candidates.extend(float(root.real) for root in roots
                          if abs(root.imag) < 1e-10
                          and x_low < root.real < x_high)
    return float(min(np.polynomial.polynomial.polyval(value, derivative)
                     for value in candidates))


def _candidate_fit(points, channel, protocol, *, p_min, p_max, degree,
                   nu_bounds, goodness_alpha, correction=False,
                   warm_start=None):
    column = f"{protocol}_failures"
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
    bounds.extend([(-15.0, 15.0)] + [(-30.0, 30.0)] * degree)
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
        crossings = adjacent_crossings(points, channel, protocol)
        estimate = (float(np.median([item["p"] for item in crossings]))
                    if crossings else (p_min + p_max) / 2)
        center = float(np.clip((estimate - p_min) / scale, 0.02, 0.98))
        correlation = np.corrcoef(probabilities, observed)[0, 1]
        slope = -2.0 if np.isfinite(correlation) and correlation < 0 else 2.0
        starts = []
        for candidate_center in sorted(set((0.25, 0.5, 0.75, center))):
            for exponent in (0.6, 1.2, 2.5):
                parameters = [candidate_center, math.log(exponent),
                              intercept, slope] + [0.0] * (degree - 1)
                if correction:
                    parameters.append(0.0)
                starts.append(np.asarray(parameters, dtype=float))
    optimizations = [minimize(
        objective, start, method="L-BFGS-B", jac=True, bounds=bounds,
        options={"maxiter": 1500, "ftol": 1e-12, "gtol": 1e-8},
    ) for start in starts]
    finite = [result for result in optimizations if np.isfinite(result.fun)]
    if not finite:
        raise RuntimeError("all likelihood optimizations returned nonfinite values")
    converged = [result for result in finite if result.success]
    best = min(converged or finite, key=lambda result: result.fun)
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
    residuals = ((observed - fitted)
                 / np.sqrt(np.maximum(fitted * (1 - fitted) / shots, 1e-30)))
    names = ["pc", "nu"] + [f"a{power}" for power in range(degree + 1)]
    if correction:
        names.append("correction_omega1")
    bound_hits = [name for name, value, (low, high) in zip(names, best.x, bounds)
                  if min(value - low, high - value) <= 0.002 * (high - low)]
    derivative_min = _minimum_polynomial_derivative(
        coefficients, float(np.min(x)), float(np.max(x)))
    return {
        "pc": pc, "nu": nu, "coefficients": coefficients.tolist(),
        "polynomial_degree": degree, "p_scale": scale,
        "L_reference": reference, "p_min": p_min, "p_max": p_max,
        "nu_bounds": list(nu_bounds),
        "correction_omega": 1.0 if correction else None,
        "correction_coefficient": correction_coefficient,
        "deviance": deviance, "degrees_of_freedom": degrees,
        "deviance_p_value": goodness,
        "optimizer_success": bool(best.success),
        "optimizer_message": str(best.message),
        "bound_hits": bound_hits,
        "monotone_increasing": bool(derivative_min > 1e-7),
        "derivative_min": derivative_min,
        "parameters": best.x.tolist(),
        "fitted_probabilities": fitted.tolist(),
        "standardized_residuals": residuals.tolist(),
        "largest_absolute_residual": float(np.max(np.abs(residuals))),
        "issues": [],
    }


def _candidate_guard_issues(candidate, *, p_min, p_max, goodness_alpha,
                            include_goodness=True):
    issues = list(candidate.get("issues", []))
    if not candidate.get("optimizer_success", False):
        issues.append("optimizer did not converge")
    bound_hits = candidate.get("bound_hits", [])
    if bound_hits:
        issues.append("parameter bound reached: " + ", ".join(bound_hits))
    if not candidate.get("monotone_increasing", False):
        issues.append("fitted failure probability is not monotone increasing in p")
    pc = candidate.get("pc")
    if pc is None or not math.isfinite(pc) or not p_min < pc < p_max:
        issues.append("candidate p_c is outside the fit window")
    if include_goodness:
        goodness = candidate.get("deviance_p_value")
        if goodness is None:
            issues.append("nonpositive goodness-of-fit degrees of freedom")
        elif goodness < goodness_alpha:
            issues.append(
                f"scaling-model lack of fit (deviance p < {goodness_alpha:g})")
    return list(dict.fromkeys(issues))


def fit_feedback_protocol(points, channel, protocol, *, p_min=None, p_max=None,
                          min_size=None, degree=2, nu_bounds=(0.3, 5.0),
                          goodness_alpha=0.01, correction=False):
    """Fit one channel/protocol and expose a threshold only after every guard."""
    if channel not in ("x_only", "z_only"):
        raise ValueError("channel must be x_only or z_only")
    if protocol not in PROTOCOL_BITS:
        raise ValueError(f"protocol must be one of {tuple(PROTOCOL_BITS)}")
    points = validate_feedback_points(points)
    channel_points = [row for row in points if row["channel"] == channel]
    if not channel_points:
        raise ValueError(f"no points for channel {channel}")
    all_probabilities = sorted({row["p"] for row in channel_points})
    all_sizes = sorted({row["steps"] for row in channel_points})
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
    selected = _selected(points, channel, p_min, p_max, min_size)
    sizes = sorted({row["steps"] for row in selected})
    probabilities = sorted({row["p"] for row in selected})
    base = {
        "channel": channel, "protocol": protocol,
        "status": "insufficient_data", "pc": None, "nu": None,
        "candidate": None, "diagnostic": [], "crossings": [],
        "selected_sizes": sizes, "selected_probabilities": probabilities,
        "selected_points": len(selected), "fit_window": [p_min, p_max],
        "min_size": min_size, "degree": degree,
        "nu_bounds": list(nu_bounds), "goodness_alpha": goodness_alpha,
        "correction_omega": 1.0 if correction else None,
    }
    if len(sizes) < 4:
        base["diagnostic"].append("need at least four sizes in the fit window")
        return base
    if len(probabilities) < 5:
        base["diagnostic"].append(
            "need at least five probabilities in the fit window")
        return base
    parameter_count = degree + 3 + int(correction)
    if len(selected) <= parameter_count:
        base["diagnostic"].append("need more points than fitted parameters")
        return base
    crossings = adjacent_crossings(selected, channel, protocol)
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
        selected, channel, protocol, p_min=p_min, p_max=p_max,
        degree=degree, nu_bounds=nu_bounds, goodness_alpha=goodness_alpha,
        correction=correction)
    issues = _candidate_guard_issues(
        candidate, p_min=p_min, p_max=p_max, goodness_alpha=goodness_alpha)
    candidate = dict(candidate, issues=issues)
    base["candidate"] = candidate
    base["diagnostic"].extend(issues)
    if issues:
        base["status"] = "inconclusive"
        return base
    base.update(status="ok", pc=candidate["pc"], nu=candidate["nu"])
    return base


def bootstrap_feedback_protocol(points, channel, protocol, fit, *,
                                replicates=1000, seed=1234):
    """Seeded marginal-binomial bootstrap conditional on one accepted model."""
    if replicates < 0:
        raise ValueError("bootstrap replicates must be nonnegative")
    if fit.get("status") != "ok" or fit.get("candidate") is None:
        return {
            "replicates": replicates, "accepted": 0,
            "pc_interval": None, "nu_interval": None, "samples": [],
            "method": "not run because the primary fit was not accepted",
        }
    points = validate_feedback_points(points)
    selected = _selected(points, channel, *fit["fit_window"], fit["min_size"])
    candidate = fit["candidate"]
    probabilities = candidate["fitted_probabilities"]
    rng = np.random.default_rng(seed)
    samples = []
    correction = fit["correction_omega"] is not None
    for replicate in range(replicates):
        sampled = [dict(row) for row in selected]
        column = f"{protocol}_failures"
        for row, probability in zip(sampled, probabilities):
            row[column] = int(rng.binomial(row["shots"], probability))
        replica = _candidate_fit(
            sampled, channel, protocol, p_min=fit["fit_window"][0],
            p_max=fit["fit_window"][1], degree=fit["degree"],
            nu_bounds=tuple(fit["nu_bounds"]),
            goodness_alpha=fit["goodness_alpha"], correction=correction,
            warm_start=candidate["parameters"])
        numerical_issues = _candidate_guard_issues(
            replica, p_min=fit["fit_window"][0], p_max=fit["fit_window"][1],
            goodness_alpha=fit["goodness_alpha"], include_goodness=False)
        accepted = not numerical_issues
        samples.append({
            "replicate": replicate, "accepted": accepted,
            "pc": replica["pc"] if accepted else None,
            "nu": replica["nu"] if accepted else None,
            "status": "accepted" if accepted else "rejected",
            "diagnostic": "; ".join(numerical_issues),
        })
    accepted_samples = [sample for sample in samples if sample["accepted"]]
    required = max(2, math.ceil(0.8 * replicates)) if replicates else 0
    intervals_available = replicates > 0 and len(accepted_samples) >= required
    pc_interval = (np.quantile(
        [sample["pc"] for sample in accepted_samples], [0.025, 0.975]).tolist()
        if intervals_available else None)
    nu_interval = (np.quantile(
        [sample["nu"] for sample in accepted_samples], [0.025, 0.975]).tolist()
        if intervals_available else None)
    return {
        "replicates": replicates, "accepted": len(accepted_samples),
        "pc_interval": pc_interval, "nu_interval": nu_interval,
        "samples": samples,
        "method": "seeded parametric marginal-binomial bootstrap conditional on the accepted scaling model",
    }


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


def sensitivity_fits(points, channel, protocol, primary):
    """Record every required finite-size/model sensitivity, accepted or not."""
    sizes = primary["selected_sizes"]
    probabilities = primary["selected_probabilities"]
    low, high = primary["fit_window"]
    variants = [("cubic", low, high, primary["min_size"], 3, False)]
    if len(sizes) >= 5:
        variants.append((f"min_size_{sizes[1]}", low, high, sizes[1],
                         primary["degree"], False))
    else:
        variants.append(("min_size_unavailable", low, high,
                         sizes[-1] if sizes else primary["min_size"],
                         primary["degree"], False))
    if len(probabilities) >= 6:
        variants.extend([
            ("trim_lowest_p", probabilities[1], high, primary["min_size"],
             primary["degree"], False),
            ("trim_highest_p", low, probabilities[-2], primary["min_size"],
             primary["degree"], False),
        ])
    else:
        variants.extend([
            ("trim_lowest_p", low, high, primary["min_size"],
             primary["degree"], False),
            ("trim_highest_p", low, high, primary["min_size"],
             primary["degree"], False),
        ])
    variants.append(("fixed_L_inverse_correction", low, high,
                     primary["min_size"], primary["degree"], True))
    output = []
    for name, variant_low, variant_high, minimum, degree, correction in variants:
        fit = fit_feedback_protocol(
            points, channel, protocol, p_min=variant_low, p_max=variant_high,
            min_size=minimum, degree=degree,
            nu_bounds=tuple(primary["nu_bounds"]),
            goodness_alpha=primary["goodness_alpha"], correction=correction)
        output.append({
            "fit_key": f"{channel}/{protocol}", "channel": channel,
            "protocol": protocol, "name": name,
            "fit_window": [variant_low, variant_high], "min_size": minimum,
            "degree": degree, "correction_omega": 1.0 if correction else None,
            "status": fit["status"], "pc": fit["pc"], "nu": fit["nu"],
            "diagnostic": fit["diagnostic"], "candidate": fit["candidate"],
        })
    return output


def analyze_feedback_scan(points, *, bootstrap_replicates=1000, seed=1234,
                          p_windows=None, min_sizes=None, degree=2,
                          nu_bounds=(0.3, 5.0), goodness_alpha=0.01,
                          include_sensitivities=True):
    """Analyze every selected channel/protocol with isolated deterministic seeds."""
    points = validate_feedback_points(points)
    channels = sorted({row["channel"] for row in points})
    p_windows = {} if p_windows is None else dict(p_windows)
    min_sizes = {} if min_sizes is None else dict(min_sizes)
    combinations = [(channel, protocol) for channel in channels
                    for protocol in PROTOCOL_BITS]
    seed_sequences = np.random.SeedSequence(seed).spawn(len(combinations))
    fits = {}
    sensitivity = []
    for (channel, protocol), seed_sequence in zip(combinations, seed_sequences):
        channel_points = [row for row in points if row["channel"] == channel]
        probabilities = sorted({row["p"] for row in channel_points})
        sizes = sorted({row["steps"] for row in channel_points})
        p_min, p_max = p_windows.get(
            channel, (probabilities[0], probabilities[-1]))
        min_size = min_sizes.get(channel, sizes[0])
        primary = fit_feedback_protocol(
            points, channel, protocol, p_min=p_min, p_max=p_max,
            min_size=min_size, degree=degree, nu_bounds=nu_bounds,
            goodness_alpha=goodness_alpha)
        bootstrap = bootstrap_feedback_protocol(
            points, channel, protocol, primary,
            replicates=bootstrap_replicates,
            seed=int(seed_sequence.generate_state(1, dtype=np.uint64)[0]))
        key = f"{channel}/{protocol}"
        fits[key] = _finalize_bootstrap_guard(primary, bootstrap)
        if include_sensitivities:
            sensitivity.extend(
                sensitivity_fits(points, channel, protocol, primary))
    model_spread = {}
    for key, fit in fits.items():
        accepted = ([fit["pc"]] if fit["pc"] is not None else [])
        accepted.extend(row["pc"] for row in sensitivity
                        if row["fit_key"] == key and row["pc"] is not None)
        model_spread[key] = ([min(accepted), max(accepted)]
                             if accepted else None)
    return {
        "fits": fits, "sensitivity": sensitivity,
        "accepted_model_pc_spread": model_spread,
        "bootstrap_replicates": bootstrap_replicates, "seed": seed,
    }


def _model_probability(candidate, p, size):
    x = ((np.asarray(p) - candidate["pc"]) / candidate["p_scale"]
         * (size / candidate["L_reference"]) ** (1 / candidate["nu"]))
    eta = np.polynomial.polynomial.polyval(x, candidate["coefficients"])
    if candidate["correction_omega"] is not None:
        eta = (eta + candidate["correction_coefficient"]
               * candidate["L_reference"] / size)
    return expit(eta)


def _make_plot(report, points, output):
    os.environ.setdefault(
        "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "toponoise-mpl"))
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams["svg.fonttype"] = "none"
    import matplotlib.pyplot as plt

    channels = sorted({row["channel"] for row in points})
    figure, axes = plt.subplots(
        len(channels), 3, figsize=(15.2, 4.8 * len(channels)),
        squeeze=False, constrained_layout=True)
    protocol_styles = {
        "none": ("o", "-", "No feedback (final decode only)", "#2364AA"),
        "slice_local": ("s", "--", "Slice-local feedback", "#D1495B"),
        "full_prefix": ("^", ":", "Full-prefix feedback", "#258477"),
    }
    pair_styles = {
        "none_vs_slice_local": ("o", "-", "Slice-local - no feedback"),
        "none_vs_full_prefix": ("s", "--", "Full-prefix - no feedback"),
        "slice_local_vs_full_prefix": ("^", ":", "Full-prefix - slice-local"),
    }
    for row_index, channel in enumerate(channels):
        raw_axis, collapse_axis, paired_axis = axes[row_index]
        selected_channel = [row for row in points if row["channel"] == channel]
        sizes = sorted({row["steps"] for row in selected_channel})
        size_colors = plt.get_cmap("viridis")(
            np.linspace(0.08, 0.92, len(sizes)))
        for size, size_color in zip(sizes, size_colors):
            selected = [row for row in selected_channel if row["steps"] == size]
            p = np.asarray([row["p"] for row in selected])
            for protocol, (marker, line, label, _) in protocol_styles.items():
                column = f"{protocol}_failures"
                rates = np.asarray([row[column] / row["shots"] for row in selected])
                intervals = np.asarray([
                    wilson_interval(row[column], row["shots"])
                    for row in selected]).T
                bars = np.maximum(
                    0, np.vstack((rates - intervals[0], intervals[1] - rates)))
                raw_axis.errorbar(
                    p, rates, yerr=bars, marker=marker, ls=line,
                    color=size_color, ms=3.2, lw=0.9, capsize=1.5,
                    alpha=0.82,
                    label=f"L={size}, {label}" if size == sizes[0] else f"L={size}")
            for pair_name, (marker, line, label) in pair_styles.items():
                differences = np.asarray([
                    point["paired_differences"][pair_name]["difference"]
                    for point in selected])
                intervals = np.asarray([
                    point["paired_differences"][pair_name]["ci"]
                    for point in selected]).T
                bars = np.maximum(0, np.vstack(
                    (differences - intervals[0], intervals[1] - differences)))
                paired_axis.errorbar(
                    p, differences, yerr=bars, marker=marker, ls=line,
                    color=size_color, ms=3.2, lw=0.9, capsize=1.5,
                    alpha=0.82,
                    label=f"L={size}, {label}" if size == sizes[0] else f"L={size}")

        for protocol, (_, _, label, color) in protocol_styles.items():
            key = f"{channel}/{protocol}"
            fit = report["fits"][key]
            candidate = fit.get("candidate")
            if fit["status"] != "ok" or candidate is None:
                continue
            fitted_points = _selected(
                points, channel, *fit["fit_window"], fit["min_size"])
            x = np.asarray([
                ((point["p"] - candidate["pc"]) / candidate["p_scale"])
                * (point["steps"] / candidate["L_reference"])
                ** (1 / candidate["nu"])
                for point in fitted_points])
            y = np.asarray([
                point[f"{protocol}_failures"] / point["shots"]
                for point in fitted_points])
            collapse_axis.scatter(x, y, s=12, alpha=0.55, color=color,
                                  label=label)
            grid = np.linspace(float(np.min(x)), float(np.max(x)), 300)
            curve = expit(np.polynomial.polynomial.polyval(
                grid, candidate["coefficients"]))
            collapse_axis.plot(grid, curve, color=color, lw=2)

        channel_label = "X-only" if channel == "x_only" else "Z-only"
        state = "|00>" if channel == "x_only" else "|++>"
        raw_axis.set(
            title=f"{channel_label}: raw state-failure curves",
            xlabel="Fault probability p per participating gate output",
            ylabel=f"Logical {state} state-failure probability")
        collapse_axis.set(
            title=f"{channel_label}: Scaling collapse (accepted fits only)",
            xlabel="x = ((p-pc)/sp)(L/Lref)^(1/nu)",
            ylabel="State-failure probability")
        paired_axis.axhline(0, color="black", lw=0.8, ls="--")
        paired_axis.set(
            title=f"{channel_label}: Paired differences",
            xlabel="Fault probability p per participating gate output",
            ylabel="Second protocol - first protocol")
        raw_axis.legend(fontsize=5.5, ncol=2)
        if collapse_axis.collections:
            collapse_axis.legend(fontsize=7)
        else:
            collapse_axis.text(
                0.5, 0.5, "Threshold withheld by fit guards",
                transform=collapse_axis.transAxes, ha="center", va="center")
        paired_axis.legend(fontsize=5.5, ncol=2)
        for axis in (raw_axis, collapse_axis, paired_axis):
            axis.grid(alpha=0.2)
    figure.suptitle(
        "Diagonal 3D toric code: No feedback, Slice-local, and Full-prefix feedback\n"
        "Marginal threshold fits use paired shots separately; joint counts determine Paired differences",
        fontsize=11)
    for extension in ("pdf", "svg", "png"):
        target = output / f"threshold-fit.{extension}"
        temporary = output / f"threshold-fit.tmp.{extension}"
        figure.savefig(temporary, dpi=180)
        temporary.replace(target)
    plt.close(figure)


def _write_summary(path, fits):
    columns = ["fit_key", "channel", "protocol", "status", "pc", "nu",
               "bootstrap_pc_low", "bootstrap_pc_high", "degree", "p_min",
               "p_max", "min_size", "deviance", "deviance_p_value",
               "diagnostic"]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for key in sorted(fits):
            fit = fits[key]
            interval = fit["bootstrap"]["pc_interval"]
            candidate = fit.get("candidate") or {}
            writer.writerow({
                "fit_key": key, "channel": fit["channel"],
                "protocol": fit["protocol"], "status": fit["status"],
                "pc": fit["pc"], "nu": fit["nu"],
                "bootstrap_pc_low": interval[0] if interval else None,
                "bootstrap_pc_high": interval[1] if interval else None,
                "degree": fit["degree"], "p_min": fit["fit_window"][0],
                "p_max": fit["fit_window"][1], "min_size": fit["min_size"],
                "deviance": candidate.get("deviance"),
                "deviance_p_value": candidate.get("deviance_p_value"),
                "diagnostic": "; ".join(fit["diagnostic"]),
            })


def _write_bootstrap(path, fits):
    columns = ["fit_key", "channel", "protocol", "replicate", "accepted",
               "pc", "nu", "status", "diagnostic"]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for key in sorted(fits):
            fit = fits[key]
            for sample in fit["bootstrap"]["samples"]:
                writer.writerow({"fit_key": key, "channel": fit["channel"],
                                 "protocol": fit["protocol"], **sample})


def _write_sensitivity(path, rows):
    columns = ["fit_key", "channel", "protocol", "name", "p_min", "p_max",
               "min_size", "degree", "correction_omega", "status", "pc",
               "nu", "diagnostic"]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "fit_key": row["fit_key"], "channel": row["channel"],
                "protocol": row["protocol"], "name": row["name"],
                "p_min": row["fit_window"][0], "p_max": row["fit_window"][1],
                "min_size": row["min_size"], "degree": row["degree"],
                "correction_omega": row["correction_omega"],
                "status": row["status"], "pc": row["pc"], "nu": row["nu"],
                "diagnostic": "; ".join(row["diagnostic"]),
            })


def _write_points(path, points):
    columns = ["channel", "L", "p", "shots", "protocol", "failures",
               "rate", "wilson_low", "wilson_high"]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for row in points:
            for protocol in PROTOCOL_BITS:
                failures = row[f"{protocol}_failures"]
                low, high = wilson_interval(failures, row["shots"])
                writer.writerow({
                    "channel": row["channel"], "L": row["steps"],
                    "p": row["p"], "shots": row["shots"],
                    "protocol": protocol, "failures": failures,
                    "rate": failures / row["shots"],
                    "wilson_low": low, "wilson_high": high,
                })


def _write_paired_differences(path, points):
    columns = ["channel", "L", "p", "shots", "pair", "first", "second",
               "neither", "first_only", "second_only", "both", "difference",
               "standard_error", "ci_low", "ci_high"]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for row in points:
            for pair_name in PAIR_SPECS:
                pair = row["paired_differences"][pair_name]
                writer.writerow({
                    "channel": row["channel"], "L": row["steps"],
                    "p": row["p"], "shots": row["shots"], "pair": pair_name,
                    "first": pair["first"], "second": pair["second"],
                    "neither": pair["neither"],
                    "first_only": pair["first_only"],
                    "second_only": pair["second_only"], "both": pair["both"],
                    "difference": pair["difference"],
                    "standard_error": pair["standard_error"],
                    "ci_low": pair["ci"][0], "ci_high": pair["ci"][1],
                })


def _package_versions():
    versions = {"python": platform.python_version(), "numpy": np.__version__}
    for package in ("scipy", "matplotlib"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "unavailable"
    return versions


def main(arguments=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv", type=Path, nargs="+")
    parser.add_argument("--x-p-min", type=float)
    parser.add_argument("--x-p-max", type=float)
    parser.add_argument("--z-p-min", type=float)
    parser.add_argument("--z-p-max", type=float)
    parser.add_argument("--x-min-size", type=int)
    parser.add_argument("--z-min-size", type=int)
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
    if args.seed < 0:
        parser.error("--seed must be nonnegative")
    for channel, low, high in (
            ("x_only", args.x_p_min, args.x_p_max),
            ("z_only", args.z_p_min, args.z_p_max)):
        if (low is None) != (high is None):
            parser.error(f"{channel} requires both probability-window bounds")
        if low is not None and not 0 <= low < high < 0.5:
            parser.error(f"{channel} probability window must satisfy 0 <= min < max < 0.5")

    points = read_feedback_points(args.csv)
    channels = {row["channel"] for row in points}
    p_windows = {}
    if "x_only" in channels and args.x_p_min is not None:
        p_windows["x_only"] = (args.x_p_min, args.x_p_max)
    if "z_only" in channels and args.z_p_min is not None:
        p_windows["z_only"] = (args.z_p_min, args.z_p_max)
    min_sizes = {}
    if args.x_min_size is not None:
        min_sizes["x_only"] = args.x_min_size
    if args.z_min_size is not None:
        min_sizes["z_only"] = args.z_min_size
    analysis = analyze_feedback_scan(
        points, bootstrap_replicates=args.bootstrap, seed=args.seed,
        p_windows=p_windows, min_sizes=min_sizes, degree=args.degree,
        nu_bounds=(args.nu_min, args.nu_max),
        goodness_alpha=args.goodness_alpha,
        include_sensitivities=not args.no_sensitivity)
    resolved_inputs = sorted(Path(path).resolve() for path in args.csv)
    report = {
        "format_version": 1,
        **analysis,
        "observable": {
            "x_only": "logical |00> state failure: either residual logical-X membrane",
            "z_only": "logical |++> state failure: either residual logical-Z winding",
        },
        "protocols": list(PROTOCOL_BITS),
        "channels": sorted(channels),
        "input_csv": [str(path) for path in resolved_inputs],
        "input_sha256": {str(path): _sha256(path) for path in resolved_inputs},
        "source_file": str(Path(__file__).resolve()),
        "source_sha256": _sha256(__file__),
        "command_arguments": list(arguments if arguments is not None else sys.argv[1:]),
        "package_versions": _package_versions(),
        "likelihood":
            "separate marginal binomial likelihood for every channel/protocol combination",
        "protocol_dependence":
            "protocol observations are paired through replayed faults; joint counts determine paired rate differences, so threshold estimates are correlated and are not compared as independent fits",
        "scaling_variable": "x=((p-pc)/sp)(L/Lref)^(1/nu)",
        "sampling_uncertainty":
            "95% parametric marginal-bootstrap interval conditional on each accepted model",
        "model_uncertainty":
            "per-channel/protocol accepted sensitivity spread, reported separately from bootstrap intervals",
        "threshold_difference_policy":
            "no threshold-shift confidence interval is reported without a joint paired bootstrap",
    }
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    _write_summary(output / "threshold-summary.csv", report["fits"])
    _write_bootstrap(output / "threshold-bootstrap.csv", report["fits"])
    _write_sensitivity(output / "threshold-sensitivity.csv", report["sensitivity"])
    _write_points(output / "threshold-points.csv", points)
    _write_paired_differences(output / "paired-differences.csv", points)
    (output / "threshold-fit.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    _make_plot(report, points, output)
    print(json.dumps({
        key: {"status": fit["status"], "pc": fit["pc"],
              "nu": fit["nu"], "diagnostic": fit["diagnostic"]}
        for key, fit in sorted(report["fits"].items())
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
