#!/usr/bin/env python3
"""Supplemental finite-size scaling with an explicitly assumed correction.

The binomial model is logit(P) = a0 + a1*x + a2*x^2 + b*(Lref/L)^omega,
where x=((p-pc)/p_scale)*(L/Lref)^(1/nu). The correction exponent is fixed,
not estimated. Bootstrap intervals are conditional on this model; window,
size-cut and fixed-exponent sensitivity are reported separately.
"""

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import tempfile

import numpy as np
import scipy
from scipy.optimize import minimize
from scipy.special import expit, xlogy
from scipy.stats import chi2


_SPEC = importlib.util.spec_from_file_location("toric3d_baseline_fit", Path(__file__).with_name("fit_toric3d_z_threshold.py"))
_BASELINE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_BASELINE)


class CorrectionModel:
    def __init__(self, rows, *, p_min=0.0045, p_max=0.0065, min_size=9,
                 omega=1.0, nu_bounds=(0.3, 5.0)):
        if not 0 <= p_min < p_max <= 1 or omega <= 0 or not 0 < nu_bounds[0] < nu_bounds[1]:
            raise ValueError("Use an ordered probability window and positive omega/nu bounds")
        self.all_rows = _BASELINE._validated_rows(rows)
        self.rows = [row for row in self.all_rows if row["steps"] >= min_size and p_min <= row["p"] <= p_max]
        self.sizes = sorted({row["steps"] for row in self.rows})
        if len(self.sizes) < 3 or len(self.rows) <= 6:
            raise ValueError("Correction fits require at least three sizes and more than six points")
        self.L, self.p, self.n, self.k = _BASELINE._arrays(self.rows)
        self.p_min, self.p_max, self.min_size = p_min, p_max, min_size
        self.p_scale = p_max - p_min
        self.pc_bounds = (max(0.0, p_min - self.p_scale / 2), min(1.0, p_max + self.p_scale / 2))
        self.pc_span = self.pc_bounds[1] - self.pc_bounds[0]
        self.reference = float(np.exp(np.mean(np.log(self.sizes))))
        self.logL = np.log(self.L / self.reference)
        self.omega, self.nu_bounds = omega, nu_bounds
        self.correction = (self.reference / self.L) ** omega
        self.bounds = [(0.0, 1.0), tuple(np.log(nu_bounds)), (-15.0, 15.0),
                       (0.0001, 30.0), (-30.0, 30.0), (-15.0, 15.0)]

    def objective(self, parameters, counts=None):
        """Average negative binomial log likelihood and its analytic Jacobian."""
        counts = self.k if counts is None else np.asarray(counts, dtype=float)
        q, log_nu, a0, a1, a2, b = parameters
        nu = math.exp(log_nu)
        pc = self.pc_bounds[0] + self.pc_span * q
        stretch = np.exp(self.logL / nu)
        x = (self.p - pc) / self.p_scale * stretch
        eta = a0 + a1 * x + a2 * x * x + b * self.correction
        residual = self.n * expit(eta) - counts
        derivative = a1 + 2 * a2 * x
        jacobian = [-self.pc_span / self.p_scale * stretch * derivative,
                    -x * self.logL / nu * derivative, np.ones_like(x), x, x * x, self.correction]
        total = np.sum(self.n)
        loss = np.sum(self.n * np.logaddexp(0, eta) - counts * eta) / total
        gradient = np.array([np.dot(residual, value) / total for value in jacobian])
        return float(loss), gradient

    def fit(self, *, counts=None, warm_start=None, pc_fixed=None):
        counts = self.k if counts is None else np.asarray(counts, dtype=float)
        if counts.shape != self.k.shape or np.any(counts < 0) or np.any(counts > self.n):
            raise ValueError("Bootstrap counts must match the scan and satisfy 0 <= failures <= shots")
        bounds = list(self.bounds)
        if pc_fixed is not None:
            if not self.pc_bounds[0] <= pc_fixed <= self.pc_bounds[1]:
                raise ValueError("Profile pc is outside the declared parameter bounds")
            q = (pc_fixed - self.pc_bounds[0]) / self.pc_span
            bounds[0] = (q, q)
        starts = [np.asarray(warm_start, dtype=float)] if warm_start is not None else [
            np.array([q, math.log(nu), -5.0, 5.0, 0.0, b])
            for q in (0.3, 0.5, 0.7) for nu in (0.7, 1.5, 3.0) for b in (0.0, 0.5)]
        fits = [minimize(lambda pars: self.objective(pars, counts), start, method="L-BFGS-B", jac=True,
                         bounds=bounds, options=dict(maxiter=3000, ftol=1e-14, gtol=1e-10)) for start in starts]
        finite = [fit for fit in fits if np.isfinite(fit.fun)]
        if not finite:
            raise RuntimeError("All likelihood optimizations returned nonfinite objectives")
        converged = [fit for fit in finite if fit.success]
        best = min(converged or finite, key=lambda fit: fit.fun)
        q, log_nu, a0, a1, a2, b = best.x
        pc, nu = self.pc_bounds[0] + self.pc_span * q, math.exp(log_nu)
        x = (self.p - pc) / self.p_scale * (self.L / self.reference) ** (1 / nu)
        fitted = expit(a0 + a1 * x + a2 * x * x + b * self.correction)
        safe = np.clip(fitted, 1e-14, 1 - 1e-14)
        observed = counts / self.n
        deviance = float(2 * np.sum(xlogy(counts, observed / safe)
                                    + xlogy(self.n - counts, (1 - observed) / (1 - safe))))
        dof = len(self.rows) - 6 + int(pc_fixed is not None)
        goodness = float(chi2.sf(deviance, dof)) if dof > 0 else None
        names = ["pc", "nu", "a0", "a1", "a2", "b"]
        hits = [name for name, value, (low, high) in zip(names, best.x, bounds)
                if high != low and min(value - low, high - value) <= 0.001 * (high - low)]
        # The quadratic derivative is linear in x, so endpoints suffice.
        derivative_min = float(min(a1 + 2 * a2 * np.min(x), a1 + 2 * a2 * np.max(x)))
        issues = []
        if not best.success:
            issues.append("optimizer did not converge")
        if hits:
            issues.append("parameter bound reached: " + ", ".join(hits))
        if derivative_min <= 0:
            issues.append("fitted probability is not increasing throughout the fitted data range")
        if dof <= 0:
            issues.append("nonpositive goodness-of-fit degrees of freedom")
        if not np.all(np.isfinite(fitted)) or not math.isfinite(deviance):
            issues.append("nonfinite fitted probabilities or deviance")
        numerical = not issues
        if goodness is not None and goodness < 0.01:
            issues.append("scaling-model lack of fit (deviance p < 0.01)")
        extrapolation = not self.p_min <= pc <= self.p_max
        return dict(pc=float(pc), nu=float(nu), b=float(b), coefficients=[float(a0), float(a1), float(a2)],
                    omega=self.omega, L_reference=self.reference, p_scale=self.p_scale,
                    fit_window=[self.p_min, self.p_max], min_size=self.min_size, sizes=self.sizes,
                    selected_points=len(self.rows), pc_bounds=list(self.pc_bounds), nu_bounds=list(self.nu_bounds),
                    parameters=best.x.tolist(), fitted_probabilities=fitted.tolist(),
                    nll_without_binomial_constants=float(best.fun * np.sum(self.n)), deviance=deviance,
                    degrees_of_freedom=dof, deviance_p_value=goodness, derivative_min=derivative_min,
                    optimizer_success=bool(best.success), optimizer_message=str(best.message), bound_hits=hits,
                    numerically_valid=numerical, accepted=not issues, issues=issues,
                    pc_outside_fit_window=extrapolation,
                    status="inconclusive" if issues else "extrapolation" if extrapolation else "ok")


def bootstrap(model, fit, *, count=1000, seed=1234):
    if count < 0:
        raise ValueError("Bootstrap count must be nonnegative")
    rng = np.random.default_rng(seed)
    samples = []
    if fit["accepted"]:
        for iteration in range(count):
            counts = rng.binomial(model.n.astype(int), fit["fitted_probabilities"])
            sample = model.fit(counts=counts, warm_start=fit["parameters"])
            # A chi-square tail of a random replica is not an optimizer failure.
            # Reject numerical, monotonicity and boundary failures, not that tail.
            samples.append(dict(replicate=iteration, pc=sample["pc"], nu=sample["nu"], b=sample["b"],
                                accepted=sample["numerically_valid"], status=sample["status"],
                                issues="; ".join(sample["issues"])))
    accepted = [sample for sample in samples if sample["accepted"]]
    intervals = None
    if count and len(accepted) >= max(2, 0.8 * count):
        intervals = {key: np.quantile([sample[key] for sample in accepted], [0.025, 0.975]).tolist()
                     for key in ("pc", "nu", "b")}
    return dict(requested_count=count, seed=seed, accepted_count=len(accepted), samples=samples,
                intervals=intervals, confidence_level=0.95,
                method="parametric binomial bootstrap at fitted probabilities; likelihood refit for each replica",
                uncertainty_scope="conditional on fixed correction exponent and scaling model; excludes finite-size systematic effects")


def profile(model, fit, *, points=81):
    samples = []
    for pc in np.linspace(*model.pc_bounds, points):
        sample = model.fit(pc_fixed=float(pc), warm_start=fit["parameters"])
        samples.append(dict(pc=float(pc), delta_deviance=2 * (sample["nll_without_binomial_constants"]
                                                            - fit["nll_without_binomial_constants"]),
                            nu=sample["nu"], b=sample["b"], accepted=sample["numerically_valid"],
                            issues="; ".join(sample["issues"])))
    cutoff = float(chi2.ppf(0.95, 1))
    crossings = []
    for left, right in zip(samples[:-1], samples[1:]):
        a, b = left["delta_deviance"] - cutoff, right["delta_deviance"] - cutoff
        if left["accepted"] and right["accepted"] and a * b < 0:
            crossings.append(left["pc"] - a * (right["pc"] - left["pc"]) / (b - a))
    below = [value for value in crossings if value < fit["pc"]]
    above = [value for value in crossings if value > fit["pc"]]
    interval = [max(below), min(above)] if below and above else None
    return dict(samples=samples, cutoff=cutoff, interval=interval,
                interval_method="linearly interpolated profile-likelihood chi-square cutoff; conditional on fixed model")


def sensitivities(model):
    variants = []
    levels = sorted({row["p"] for row in model.rows})
    windows = [(model.p_min, model.p_max)]
    if len(levels) >= 4:
        windows += [(levels[1], model.p_max), (model.p_min, levels[-2]), (levels[1], levels[-2])]
    for low, high in windows:
        for minimum in model.sizes[:max(1, len(model.sizes) - 2)]:
            if low < high:
                variants.append((low, high, minimum, model.omega))
    variants += [(model.p_min, model.p_max, model.min_size, omega) for omega in (0.5, 2.0) if omega != model.omega]
    output = []
    for low, high, minimum, omega in variants:
        try:
            other = CorrectionModel(model.all_rows, p_min=low, p_max=high, min_size=minimum,
                                    omega=omega, nu_bounds=model.nu_bounds)
            output.append(other.fit())
        except ValueError as error:
            output.append(dict(fit_window=[low, high], min_size=minimum, omega=omega,
                               accepted=False, status="insufficient_data", issues=[str(error)]))
    return output


def make_plot(report, output):
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "toric3d-matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    fit, rows = report["fit"], report["all_rows"]
    intervals = report["bootstrap"]["intervals"]
    interval = intervals["pc"] if intervals else None
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.9), constrained_layout=True)
    sizes = sorted({row["steps"] for row in rows})
    colors = plt.get_cmap("viridis")(np.linspace(0.1, 0.9, len(sizes)))
    critical_rates = []
    for size, color in zip(sizes, colors):
        points = [row for row in rows if row["steps"] == size]
        p = np.array([row["p"] for row in points])
        rate = np.array([row["failures"] / row["shots"] for row in points])
        wilson = np.array([_BASELINE.wilson_interval(row["failures"], row["shots"]) for row in points]).T
        bars = np.maximum(0, np.vstack([rate - wilson[0], wilson[1] - rate]))
        for axis in axes:
            axis.errorbar(p, rate, yerr=bars, fmt="o", ms=4, color=color, capsize=2, label=f"L = {size}")
        if size in fit["sizes"]:
            grid = np.linspace(*fit["fit_window"], 250)
            x = (grid - fit["pc"]) / fit["p_scale"] * (size / fit["L_reference"]) ** (1 / fit["nu"])
            correction = fit["b"] * (fit["L_reference"] / size) ** fit["omega"]
            curve = expit(np.polynomial.polynomial.polyval(x, fit["coefficients"]) + correction)
            critical_rates.append(float(expit(fit["coefficients"][0] + correction)))
            for axis in axes:
                axis.plot(grid, curve, lw=1.2, color=color)
    for axis in axes:
        axis.axvline(fit["pc"], color="black", ls="--", lw=1.2, label="asymptotic p_c estimate")
        if interval:
            axis.axvspan(*interval, color="#e9ad36", alpha=0.3, label="conditional 95% bootstrap interval")
        axis.set(xlabel="Z fault probability per active gate output, p", ylabel="Probability of either residual logical Z")
        axis.grid(alpha=0.2)
        axis.ticklabel_format(axis="x", style="plain", useOffset=False)
        axis.xaxis.set_major_locator(MaxNLocator(nbins=5))
    axes[0].set_title("All simulated probabilities", fontsize=11)
    axes[0].set_ylim(0, min(1, 1.15 * max(row["failures"] / row["shots"] for row in rows)))
    axes[0].legend(fontsize=8, ncol=2, loc="upper left")
    axes[1].set_title("Threshold region: finite-size corrections retained", fontsize=10)
    width = fit["fit_window"][1] - fit["fit_window"][0]
    half_width = max(0.4 * width, 1.6 * max(abs(value - fit["pc"]) for value in interval) if interval else 0)
    axes[1].set_xlim(max(fit["fit_window"][0], fit["pc"] - half_width),
                     min(fit["fit_window"][1], fit["pc"] + half_width))
    axes[1].set_ylim(0, min(1, max(0.005, 3 * max(critical_rates))))
    subtitle = f"fixed correction exponent omega = {fit['omega']:g}; p_c = {fit['pc']:.6f}"
    if interval:
        subtitle += f", conditional 95% interval [{interval[0]:.6f}, {interval[1]:.6f}]"
    fig.suptitle("Slice-transfer Z noise with final MWPM\n" + subtitle
                 + "\nWilson 95% point bars; finite-size curves need not cross at the asymptotic p_c", fontsize=10.5)
    for extension in ("pdf", "svg", "png"):
        fig.savefig(output / f"corrected-fit.{extension}", dpi=180)
    plt.close(fig)


def _csv(path, rows, fields):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main(arguments=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--p-min", type=float, default=0.0045)
    parser.add_argument("--p-max", type=float, default=0.0065)
    parser.add_argument("--min-size", type=int, default=9)
    parser.add_argument("--omega", type=float, default=1.0)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--metadata", type=Path, action="append", default=[])
    args = parser.parse_args(arguments)
    with args.csv.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    model = CorrectionModel(rows, p_min=args.p_min, p_max=args.p_max, min_size=args.min_size, omega=args.omega)
    fit = model.fit()
    boot = bootstrap(model, fit, count=args.bootstrap, seed=args.seed)
    profiles = profile(model, fit) if fit["accepted"] else dict(samples=[], interval=None)
    sensitivity = sensitivities(model)
    issues = list(fit["issues"])
    if args.bootstrap and boot["intervals"] is None:
        issues.append("Fewer than 80% usable bootstrap fits; no reliable interval")
    if fit["accepted"] and profiles["interval"] is None:
        issues.append("Profile interval does not close within the declared parameter bounds")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    report = dict(status="ok" if not issues else "inconclusive", diagnostic=issues,
                  observable="any residual transverse logical Z after unit-data MWPM", fit=fit,
                  model="logit(P) = a0 + a1*x + a2*x^2 + b*(L_reference/L)^omega",
                  scaling_variable="x = ((p-pc)/p_scale)*(L/L_reference)^(1/nu)",
                  assumptions="Fixed omega is an explicit model assumption; finite-size curves need not share a crossing",
                  bootstrap=boot, profile=profiles, sensitivity=sensitivity, all_rows=model.all_rows,
                  uncertainty_scope="Bootstrap/profile intervals are conditional statistical intervals, not total uncertainty; sensitivity ranges are not confidence intervals",
                  input_csv=str(args.csv.resolve()), input_sha256=hashlib.sha256(args.csv.read_bytes()).hexdigest(),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  versions=dict(numpy=np.__version__, scipy=scipy.__version__), metadata_files=[])
    metadata = set(args.metadata) | set(args.csv.parent.glob("*.toml")) | set(args.csv.parent.glob("*metadata*.json"))
    for source in sorted(metadata):
        target = output / "metadata" / source.name
        target.parent.mkdir(exist_ok=True)
        if source.resolve() != target.resolve():
            shutil.copy2(source, target)
        report["metadata_files"].append(dict(source=str(source.resolve()), copied=str(target),
                                              sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
    (output / "corrected-fit.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    _csv(output / "corrected-bootstrap.csv", boot["samples"], ["replicate", "pc", "nu", "b", "accepted", "status", "issues"])
    _csv(output / "corrected-profile.csv", profiles["samples"], ["pc", "delta_deviance", "nu", "b", "accepted", "issues"])
    _csv(output / "corrected-sensitivity.csv", [dict(row, p_min=row["fit_window"][0], p_max=row["fit_window"][1],
                                                    issues="; ".join(row["issues"])) for row in sensitivity],
         ["p_min", "p_max", "min_size", "omega", "pc", "nu", "b", "deviance", "degrees_of_freedom", "deviance_p_value", "accepted", "pc_outside_fit_window", "status", "issues"])
    points = []
    for row in model.all_rows:
        low, high = _BASELINE.wilson_interval(row["failures"], row["shots"])
        points.append(dict(row, rate=row["failures"] / row["shots"], wilson_low=low, wilson_high=high,
                           used_in_fit=row in model.rows))
    _csv(output / "corrected-points.csv", points, ["steps", "p", "shots", "failures", "rate", "wilson_low", "wilson_high", "used_in_fit"])
    make_plot(report, output)
    print(json.dumps(dict(status=report["status"], pc=fit["pc"], nu=fit["nu"], b=fit["b"],
                          interval=boot["intervals"], accepted_bootstraps=boot["accepted_count"],
                          deviance=fit["deviance"], dof=fit["degrees_of_freedom"], goodness=fit["deviance_p_value"],
                          profile_interval=profiles["interval"], diagnostic=issues), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
