"""Pilot crossing selection and joint near-transition finite-size fits.

The local quadratic fitter follows the previous plus-input study in this
repository, with explicit joint counts and failure reporting for every panel.
"""
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import least_squares
from scipy.stats import chi2

from sampler import atomic_json, stable_seed

RUN = Path(__file__).resolve().parent
DISTANCES = (9, 11, 13, 15)
CHANNELS = ("y_only", "hadamard")
P_SCALE = 0.01
D_REF = 11.0
REFERENCE = "https://panqec.readthedocs.io/en/latest/tutorials/Computing%20threshold.html"


def counts_and_rate(row, component):
    counts = row["counts"]
    k = counts[1 if component == "x" else 2] + counts[3]
    n = row["shots"]
    rate = k / n
    return k, rate, np.sqrt(max(rate * (1 - rate), 0.5 / n) / n)


def select_crossing(rows, channel, component):
    """Select an interior cluster with the expected distance ordering on both sides.

    All candidate and rejection information is saved. p=0 and rates outside
    (0.005,0.45) cannot nominate a crossing. Three adjacent-distance crossings
    must lie within 0.010; d=9 vs d=15 must reverse ordering at >=1.5 combined
    standard errors at center +/-0.010. These are pilot selection rules, not
    the statistical significance or uncertainty of the final fit.
    """
    curves = {}
    for d in DISTANCES:
        selected = sorted((r for r in rows if r["channel"] == channel and r["distance"] == d),
                          key=lambda r: r["p"])
        if not selected:
            return dict(status="unresolved", reason="Missing pilot distance", centers=[])
        curves[d] = {r["p"]: counts_and_rate(r, component)[1:] for r in selected}
    grid = sorted(set.intersection(*(set(c) for c in curves.values())))
    candidates = []
    for low, high in zip(DISTANCES[:-1], DISTANCES[1:]):
        crosses = []
        for p0, p1 in zip(grid[:-1], grid[1:]):
            a0, a1 = curves[low][p0][0], curves[low][p1][0]
            b0, b1 = curves[high][p0][0], curves[high][p1][0]
            delta0, delta1 = a0 - b0, a1 - b1
            if p0 == 0 or not (delta0 >= 0 and delta1 <= 0 and delta0 > delta1):
                continue
            fraction = delta0 / (delta0 - delta1)
            value = (1 - fraction) * (a0 + b0) / 2 + fraction * (a1 + b1) / 2
            if 0.005 < value < 0.45:
                crosses.append(p0 + fraction * (p1 - p0))
        candidates.append(sorted(set(crosses)))
    common = dict(channel=channel, component=component,
                  adjacent_pair_candidates=candidates, pilot_rate_cut=[0.005, 0.45])

    def unresolved_or_refine(reason):
        # A noisy adjacent pair can hide a crossing already bracketed by the
        # extreme distances. This nominates MORE PILOT SAMPLING, never a fit.
        centers = []
        for p0, p1 in zip(grid[:-1], grid[1:]):
            a0, a1 = curves[9][p0][0], curves[9][p1][0]
            b0, b1 = curves[15][p0][0], curves[15][p1][0]
            delta0, delta1 = a0-b0, a1-b1
            if p0 == 0 or not (delta0 >= 0 and delta1 <= 0 and delta0 > delta1):
                continue
            fraction = delta0/(delta0-delta1)
            value = (1-fraction)*(a0+b0)/2 + fraction*(a1+b1)/2
            if not 0.005 < value < 0.45:
                continue
            center = round(round((p0+fraction*(p1-p0))/0.002)*0.002, 6)
            lower = max((p for p in grid if 0 < p <= center-0.020+1e-12), default=None)
            upper = min((p for p in grid if p >= center+0.020-1e-12), default=None)
            if lower is None or upper is None:
                continue
            zscores = []
            for p, sign in ((lower, 1), (upper, -1)):
                a, sa = curves[9][p]
                b, sb = curves[15][p]
                zscores.append(sign*(a-b)/np.hypot(sa, sb))
            if min(zscores) >= 1.5:
                centers.append(center)
        centers = sorted(set(centers))
        if centers:
            return dict(**common, status="ambiguous", reason=reason + "; extreme distances bracket a candidate",
                        centers=centers, refine_rates=[p for p in grid if any(
                            abs(p-center) <= 0.025+1e-12 for center in centers)])
        return dict(**common, status="unresolved", reason=reason, centers=[], refine_rates=[])

    if any(not c for c in candidates):
        return unresolved_or_refine("No unsaturated crossing for every adjacent distance pair")
    clusters = []
    for values in itertools.product(*candidates):
        if max(values) - min(values) <= 0.010 + 1e-12:
            center = round(round(float(np.median(values)) / 0.002) * 0.002, 6)
            lower = max((p for p in grid if p <= center - 0.010 + 1e-12), default=None)
            upper = min((p for p in grid if p >= center + 0.010 - 1e-12), default=None)
            supported = False
            contrasts = []
            if lower is not None and upper is not None:
                for p, sign in ((lower, 1), (upper, -1)):
                    a, sa = curves[9][p]
                    b, sb = curves[15][p]
                    contrasts.append(sign * (a - b) / np.hypot(sa, sb))
                supported = min(contrasts) >= 1.5
            clusters.append(dict(center=center, crossings=list(values),
                                 spread=max(values) - min(values),
                                 ordering_z_scores=contrasts, supported=supported))
    if not clusters:
        return unresolved_or_refine("Adjacent crossings do not form a compact cluster")
    clusters.sort(key=lambda c: (c["spread"], c["center"]))
    supported = [c for c in clusters if c["supported"]]
    groups = []
    for cluster in supported:
        if not any(abs(cluster["center"] - c["center"]) <= 0.004 + 1e-12 for c in groups):
            groups.append(cluster)
    if len(groups) == 1:
        chosen = groups[0]
        return dict(**common, status="selected", center=chosen["center"],
                    crossings=chosen["crossings"], ordering_z_scores=chosen["ordering_z_scores"],
                    centers=[chosen["center"]], refine_rates=[])
    centers = sorted(set(c["center"] for c in clusters))
    refine = [p for p in grid if any(abs(p - center) <= 0.015 + 1e-12 for center in centers)]
    return dict(**common, status="ambiguous", reason="Pilot ordering or crossing choice requires more shots",
                centers=centers, refine_rates=refine, clusters=clusters)


def model(theta, p, d):
    u = (p / P_SCALE - theta[0]) * (d / D_REF) ** (1 / theta[1])
    return np.polynomial.polynomial.polyval(u, theta[2:])


def fit_local(data, lower, upper, degree=2, initial=None, multistart=True):
    p, d, n, y = (np.asarray(data[k]) for k in ("p", "d", "n", "y"))
    if len(y) <= degree + 3 or len(set(d)) < 3 or np.ptp(y) < 1e-10:
        raise ValueError("Insufficient variation or data to identify a threshold")
    bounded = np.clip(y, 0.5 / n, 1 - 0.5 / n)
    sigma = np.sqrt(bounded * (1 - bounded) / n)
    lo = [lower / P_SCALE, 0.3, 1e-8, 0.0] + [-np.inf] * (degree - 1)
    hi = [upper / P_SCALE, 5.0, 1 - 1e-8, np.inf] + [np.inf] * (degree - 1)

    def residual(theta):
        return (model(theta, p, d) - y) / sigma

    starts = [] if initial is None else [np.array(initial)]
    if multistart or not starts:
        starts.extend(np.array([(lower + upper) / (2 * P_SCALE), nu, float(np.median(y)), 0.1]
                               + [0.0] * (degree - 1)) for nu in (0.7, 1.3, 2.0, 3.0))
    solutions = [least_squares(residual, s, bounds=(lo, hi), x_scale="jac",
                               ftol=1e-11, xtol=1e-11, gtol=1e-9, max_nfev=2000) for s in starts]
    converged = [s for s in solutions if s.success and np.isfinite(s.cost)]
    if not converged:
        raise ValueError("Optimizer did not converge")
    result = min(converged, key=lambda s: s.cost)
    theta = result.x
    predicted = model(theta, p, d)
    if not (lower + 1e-7 < theta[0] * P_SCALE < upper - 1e-7 and
            0.30001 < theta[1] < 4.99999 and np.all((predicted > 0) & (predicted < 1))):
        raise ValueError("Fit reached a parameter boundary or nonphysical probability")
    u = (p / P_SCALE - theta[0]) * (d / D_REF) ** (1 / theta[1])
    slope = np.polynomial.polynomial.polyval(np.linspace(u.min(), u.max(), 101),
                                            np.polynomial.polynomial.polyder(theta[2:]))
    if not np.all(slope > 0):
        raise ValueError("Scaling function is not increasing throughout the fitted window")
    if np.linalg.matrix_rank(result.jac) != len(theta):
        raise ValueError("Unidentifiable fit parameters")
    covariance = np.linalg.inv(result.jac.T @ result.jac)
    errors = np.sqrt(np.diag(covariance))
    dof = len(y) - len(theta)
    norm = P_SCALE * D_REF ** (1 / theta[1])
    return dict(status="fitted", window=[lower, upper], degree=degree,
                distances=sorted(set(d.astype(int).tolist())), points=len(y),
                parameters=len(theta), theta_internal=theta.tolist(),
                coefficients=[float(c / norm**i) for i, c in enumerate(theta[2:])],
                p_c=float(theta[0]*P_SCALE), nu=float(theta[1]),
                p_c_covariance_se=float(errors[0]*P_SCALE), nu_covariance_se=float(errors[1]),
                chi_square=float(2*result.cost), degrees_of_freedom=dof,
                reduced_chi_square=float(2*result.cost/dof),
                goodness_of_fit_p=float(chi2.sf(2*result.cost, dof)),
                start_objective_range=float(max(s.cost for s in converged)-result.cost),
                jacobian_condition=float(np.linalg.cond(result.jac)))


def resample_joint(rows, rng):
    return np.array([rng.multinomial(r["shots"], np.array(r["counts"]) / r["shots"])
                     for r in rows], dtype=np.int64)


def data_for(rows, component, counts=None):
    counts = np.array([r["counts"] for r in rows]) if counts is None else counts
    n = np.array([r["shots"] for r in rows])
    return dict(p=np.array([r["p"] for r in rows]), d=np.array([r["distance"] for r in rows]),
                n=n, y=(counts[:, 1 if component == "x" else 2] + counts[:, 3]) / n)


def load_points(directory=RUN):
    return sorted([json.loads(path.read_text()) for path in (Path(directory) / "checkpoints").glob("*.json")],
                  key=lambda r: (r["stage"], r["channel"], r["distance"], r["p"]))


def write_csv(path, rows):
    if not rows:
        Path(path).write_text("")
        return
    columns = list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def fit_all(replicates=1000):
    selection = json.loads((RUN / "selection.json").read_text())
    rows = [r for r in load_points() if r["stage"] == "production"]
    primary, sensitivity, samples = {}, [], []
    for channel in CHANNELS:
        channel_rows = [r for r in rows if r["channel"] == channel]
        masks = {}
        for component in ("x", "z"):
            key = f"{channel}/{component}"
            choice = selection["panels"][key]
            if choice["status"] != "selected":
                primary[key] = dict(status="unresolved", reason=choice["reason"])
                continue
            center = choice["center"]
            for degree, halfwidth in itertools.product((2, 3), (0.004, 0.006, 0.008, 0.010)):
                lower, upper = max(0.0, round(center-halfwidth, 6)), round(center+halfwidth, 6)
                selected = [r for r in channel_rows if lower-1e-12 <= r["p"] <= upper+1e-12]
                try:
                    fitted = fit_local(data_for(selected, component), lower, upper, degree=degree)
                except (ValueError, np.linalg.LinAlgError) as error:
                    fitted = dict(status="invalid", reason=str(error), window=[lower, upper], degree=degree)
                sensitivity.append(dict(panel=key, kind="window", halfwidth=halfwidth, **fitted))
                if degree == 2 and halfwidth == 0.006:
                    primary[key] = dict(fitted)
                    masks[component] = np.array([lower-1e-12 <= r["p"] <= upper+1e-12
                                                for r in channel_rows])
            for distances in ((11, 13, 15), (9, 11, 13)):
                lower, upper = max(0.0, round(center-0.006, 6)), round(center+0.006, 6)
                selected = [r for r in channel_rows if r["distance"] in distances and
                            lower-1e-12 <= r["p"] <= upper+1e-12]
                try:
                    fitted = fit_local(data_for(selected, component), lower, upper)
                except (ValueError, np.linalg.LinAlgError) as error:
                    fitted = dict(status="invalid", reason=str(error), distances=list(distances))
                sensitivity.append(dict(panel=key, kind="distance", **fitted))
        # A single multinomial resample supplies BOTH logical marginals for a panel pair.
        seed = stable_seed(20260923, "joint-bootstrap", channel)
        rng = np.random.default_rng(seed)
        active = [c for c in ("x", "z") if primary[f"{channel}/{c}"]["status"] == "fitted"]
        for replicate in range(replicates if active else 0):
            joint = resample_joint(channel_rows, rng)
            for component in active:
                key = f"{channel}/{component}"
                original = primary[key]
                mask = masks[component]
                chosen = [r for r, keep in zip(channel_rows, mask) if keep]
                data = data_for(chosen, component, joint[mask])
                try:
                    try:
                        fitted = fit_local(data, *original["window"], initial=original["theta_internal"],
                                           multistart=False)
                    except (ValueError, np.linalg.LinAlgError):
                        fitted = fit_local(data, *original["window"])
                    samples.append(dict(panel=key, replicate=replicate,
                                        p_c=fitted["p_c"], nu=fitted["nu"]))
                except (ValueError, np.linalg.LinAlgError):
                    continue
            if (replicate+1) % 100 == 0:
                print(f"{channel}: bootstrap {replicate+1}/{replicates}", flush=True)
        for component in active:
            key = f"{channel}/{component}"
            fit = primary[key]
            values = np.array([[s["p_c"], s["nu"]] for s in samples if s["panel"] == key])
            fit.update(bootstrap_seed=seed, bootstrap_replicates=replicates,
                       bootstrap_successes=len(values), bootstrap_failures=replicates-len(values))
            if len(values) >= 2:
                fit.update(p_c_standard_error=float(values[:, 0].std(ddof=1)),
                           nu_standard_error=float(values[:, 1].std(ddof=1)),
                           p_c_interval_95=np.quantile(values[:, 0], [0.025, 0.975]).tolist(),
                           nu_interval_95=np.quantile(values[:, 1], [0.025, 0.975]).tolist())
            fit["status"] = ("success" if len(values) >= max(2, 0.95*replicates)
                             and fit["goodness_of_fit_p"] >= 0.01 else "unreliable")
            if fit["status"] != "success":
                fit["reason"] = "Insufficient bootstrap convergence or primary goodness-of-fit p < 0.01"
            print(key, json.dumps(fit, sort_keys=True), flush=True)
    output = dict(method="Joint weighted local quadratic finite-size scaling", reference=REFERENCE,
                  model="P_L=A+B*x+C*x^2; x=(p-p_c)*d^(1/nu)",
                  bootstrap="Multinomial joint-count resampling, shared across X/Z marginals",
                  uncertainty="Statistical, conditional on frozen pilot windows and scaling model",
                  nu_search_bounds=[0.3, 5.0], goodness_of_fit_minimum_p=0.01,
                  raw_sha256=hashlib.sha256((RUN / "raw_counts.csv").read_bytes()).hexdigest(),
                  selection_sha256=hashlib.sha256((RUN / "selection.json").read_bytes()).hexdigest(),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  numpy_version=np.__version__, scipy_version=scipy.__version__, primary_fits=primary)
    atomic_json(RUN / "local_fits.json", output)
    atomic_json(RUN / "sensitivity.json", sensitivity)
    write_csv(RUN / "bootstrap_samples.csv", samples)
    write_csv(RUN / "fit_summary.csv", [dict(panel=k, **v) for k, v in primary.items()])
    write_csv(RUN / "sensitivity.csv", sensitivity)
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bootstrap", type=int, default=1000)
    arguments = parser.parse_args()
    if arguments.bootstrap < 2:
        parser.error("--bootstrap must be at least 2")
    fit_all(arguments.bootstrap)
