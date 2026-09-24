"""Joint near-transition scaling fits to the existing logical-error counts.

No simulation or circuit changes. The primary fit is a quadratic scaling
function, weighted by binomial standard errors. Window, polynomial degree,
and distance sensitivity are retained alongside count-bootstrap intervals.
"""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import least_squares
from scipy.stats import chi2

RUN = Path(__file__).resolve().parent
RAW = RUN.parent / "plus_input_logical_errors-raw.csv"
P_SCALE = 0.01
D_REF = 11.0
REPLICATES = 1000
PANELS = (("as", "z_only", 0.026, 20260924),
          ("bp", "x_only", 0.026, 20260923))
REFERENCE = "https://panqec.readthedocs.io/en/latest/tutorials/Computing%20threshold.html"


def model(theta, p, d):
    # Internal rescaling improves conditioning; it is exactly equivalent to
    # P_L = A + B*x + C*x^2 (+ D*x^3), x = (p-p_c)*d^(1/nu).
    u = (p / P_SCALE - theta[0]) * (d / D_REF) ** (1.0 / theta[1])
    return np.polynomial.polynomial.polyval(u, theta[2:])


def fit(data, lower, upper, degree=2, initial=None, multistart=True):
    p, d, n, y = (data[k] for k in ("p", "d", "n", "y"))
    # The floor matters only for a resampled zero count, never for the actual
    # selected data. Recompute weights in every bootstrap replicate.
    bounded_y = np.clip(y, 0.5 / n, 1.0 - 0.5 / n)
    sigma = np.sqrt(bounded_y * (1.0 - bounded_y) / n)
    lo = [lower / P_SCALE, 0.3, 1e-8, 0.0] + [-np.inf] * (degree - 1)
    hi = [upper / P_SCALE, 5.0, 1.0 - 1e-8, np.inf] + [np.inf] * (degree - 1)

    def residual(theta):
        return (model(theta, p, d) - y) / sigma

    starts = [] if initial is None else [np.array(initial)]
    if multistart or not starts:
        starts.extend(np.array([(lower + upper) / (2 * P_SCALE), nu, 0.075, 0.1]
                               + [0.0] * (degree - 1)) for nu in (0.7, 1.3, 2.0, 3.0))
    solutions = [least_squares(residual, start, bounds=(lo, hi), x_scale="jac",
                               ftol=1e-11, xtol=1e-11, gtol=1e-9, max_nfev=2000)
                 for start in starts]
    converged = [s for s in solutions if s.success and np.isfinite(s.cost)]
    if not converged:
        raise ValueError("Optimizer did not converge")
    result = min(converged, key=lambda s: s.cost)
    theta = result.x
    predicted = model(theta, p, d)
    if not (lower + 1e-7 < theta[0] * P_SCALE < upper - 1e-7 and
            0.30001 < theta[1] < 4.99999 and np.all((predicted > 0) & (predicted < 1))):
        raise ValueError("Fit at a search boundary or outside probability bounds")
    u = (p / P_SCALE - theta[0]) * (d / D_REF) ** (1.0 / theta[1])
    slope = np.polynomial.polynomial.polyval(
        np.linspace(u.min(), u.max(), 101),
        np.polynomial.polynomial.polyder(theta[2:]))
    if not np.all(slope > 0):
        raise ValueError("Scaling function is not increasing inside the fit window")
    if np.linalg.matrix_rank(result.jac) != len(theta):
        raise ValueError("Unidentifiable fit parameters")
    covariance = np.linalg.inv(result.jac.T @ result.jac)
    errors = np.sqrt(np.diag(covariance))
    dof = len(y) - len(theta)
    assert dof > 0
    norm = P_SCALE * D_REF ** (1.0 / theta[1])
    return {
        "window": [lower, upper], "degree": degree,
        "distances": sorted(set(d.astype(int).tolist())), "points": len(y),
        "parameters": len(theta), "theta_internal": theta.tolist(),
        "coefficients": [float(c / norm**i) for i, c in enumerate(theta[2:])],
        "p_c": float(theta[0] * P_SCALE), "nu": float(theta[1]),
        "p_c_covariance_se": float(errors[0] * P_SCALE),
        "nu_covariance_se": float(errors[1]),
        "chi_square": float(2 * result.cost), "degrees_of_freedom": dof,
        "reduced_chi_square": float(2 * result.cost / dof),
        "goodness_of_fit_p": float(chi2.sf(2 * result.cost, dof)),
        "start_objective_range": float(max(s.cost for s in converged) - result.cost),
        "jacobian_condition": float(np.linalg.cond(result.jac)),
    }


def select(rows, panel, lower, upper, distances=(9, 11, 13, 15)):
    selected = [r for r in rows
                if (r["construction"], r["error_channel"]) == panel
                and int(r["distance"]) in distances
                and lower - 1e-12 <= float(r["p_x"]) + float(r["p_z"]) <= upper + 1e-12]
    selected.sort(key=lambda r: (int(r["distance"]), float(r["p_x"]) + float(r["p_z"])))
    assert selected
    data = {
        "p": np.array([float(r["p_x"]) + float(r["p_z"]) for r in selected]),
        "d": np.array([int(r["distance"]) for r in selected]),
        "n": np.array([int(r["shots"]) for r in selected]),
        "k": np.array([int(r["logical_failures"]) for r in selected]),
    }
    data["y"] = data["k"] / data["n"]
    assert all(r["logical_observable"] == ("logical_z" if panel[1] == "z_only" else "logical_x")
               and r["clock"] == "gate_layer" for r in selected)
    assert np.all(data["n"] == 50000)
    assert np.allclose(data["y"], [float(r["logical_failure_rate"]) for r in selected], atol=1e-14)
    assert sorted(set(data["d"])) == sorted(distances)
    return data


def bootstrap(data, primary, seed):
    rng = np.random.default_rng(seed)
    samples = []
    failures = 0
    for replicate in range(REPLICATES):
        sampled = {**data, "y": rng.binomial(data["n"], data["y"]) / data["n"]}
        try:
            result = fit(sampled, *primary["window"], degree=primary["degree"],
                         initial=primary["theta_internal"], multistart=False)
        except ValueError:
            try:
                result = fit(sampled, *primary["window"], degree=primary["degree"])
            except ValueError:
                failures += 1
                continue
        samples.append((replicate, result["p_c"], result["nu"]))
    if len(samples) < 0.95 * REPLICATES:
        raise ValueError("Too many unsuccessful bootstrap fits")
    values = np.array(samples)[:, 1:]
    primary.update({
        "bootstrap_replicates": REPLICATES, "bootstrap_successes": len(samples),
        "bootstrap_failures": failures, "bootstrap_seed": seed,
        "p_c_standard_error": float(values[:, 0].std(ddof=1)),
        "nu_standard_error": float(values[:, 1].std(ddof=1)),
        "p_c_interval_95": np.quantile(values[:, 0], [0.025, 0.975]).tolist(),
        "nu_interval_95": np.quantile(values[:, 1], [0.025, 0.975]).tolist(),
    })
    return samples


def self_check():
    # Independent physical-coordinate fixture: recover its known threshold
    # and exponent despite the internal p and distance rescaling.
    p = np.tile(np.arange(0.016, 0.0321, 0.002), 4)
    d = np.repeat([9, 11, 13, 15], 9)
    x = (p - 0.023) * d ** (1.0 / 1.45)
    y = 0.08 + 1.2 * x + 3.4 * x**2
    data = {"p": p, "d": d, "n": np.full(len(p), 50000), "y": y}
    result = fit(data, 0.016, 0.032)
    assert abs(result["p_c"] - 0.023) < 1e-9
    assert abs(result["nu"] - 1.45) < 1e-7
    assert np.allclose(result["coefficients"], [0.08, 1.2, 3.4], atol=1e-7)
    assert result["chi_square"] < 1e-12
    reversed_result = fit({k: v[::-1] for k, v in data.items()}, 0.016, 0.032)
    assert abs(reversed_result["p_c"] - result["p_c"]) < 1e-10


def summary_row(panel, result):
    return dict(construction=panel[0], error_channel=panel[1],
                p_min=result["window"][0], p_max=result["window"][1],
                degree=result["degree"], distances=";".join(map(str, result["distances"])),
                **{k: result[k] for k in ("points", "p_c", "p_c_covariance_se", "nu",
                    "nu_covariance_se", "chi_square", "degrees_of_freedom",
                    "reduced_chi_square", "goodness_of_fit_p")})


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    self_check()
    print("Known-parameter and row-order checks passed.", flush=True)
    raw_hash = hashlib.sha256(RAW.read_bytes()).hexdigest()
    with RAW.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    primary_fits = {}
    window_rows, distance_rows, bootstrap_rows = [], [], []
    for construction, channel, center, seed in PANELS:
        panel = (construction, channel)
        candidates = {}
        for degree in (2, 3):
            for halfwidth in (0.004, 0.006, 0.008, 0.010):
                lower, upper = round(center - halfwidth, 6), round(center + halfwidth, 6)
                data = select(rows, panel, lower, upper)
                result = fit(data, lower, upper, degree=degree)
                candidates[degree, halfwidth] = result
                window_rows.append(summary_row(panel, result))
        # Report a representative narrow quadratic fit; the full sensitivity
        # table, including poorly fitting wider quadratic windows, is kept.
        primary = dict(candidates[2, 0.006])
        data = select(rows, panel, *primary["window"])
        samples = bootstrap(data, primary, seed)
        primary_fits[f"{construction}/{channel}"] = primary
        bootstrap_rows.extend(dict(construction=construction, error_channel=channel,
                                   replicate=i, p_c=pc, nu=nu) for i, pc, nu in samples)
        for distances in ((11, 13, 15), (9, 11, 13)):
            result = fit(select(rows, panel, *primary["window"], distances=distances),
                         *primary["window"])
            distance_rows.append(summary_row(panel, result))
        print(f"{construction}/{channel}: p_c={primary['p_c']:.8f} +/- "
              f"{primary['p_c_standard_error']:.8f}; nu={primary['nu']:.4f} +/- "
              f"{primary['nu_standard_error']:.4f}; chi2/dof="
              f"{primary['reduced_chi_square']:.3f}; bootstrap="
              f"{primary['bootstrap_successes']}/{REPLICATES}", flush=True)
    assert raw_hash == hashlib.sha256(RAW.read_bytes()).hexdigest()
    output = {
        "method": "Joint weighted least-squares quadratic finite-size scaling",
        "model": "P_L = A + B*x + C*x^2, x=(p-p_c)*d^(1/nu)",
        "method_reference": REFERENCE,
        "weights": "Inverse observed binomial variances; recomputed in every count-bootstrap replicate",
        "bootstrap": "Independent Binomial(N, observed K/N) resampling at each data point",
        "uncertainty_scope": "Statistical only, conditional on window and scaling model; no finite-size correction term",
        "raw_csv": str(RAW.relative_to(RUN.parent)), "raw_sha256": raw_hash,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "numpy_version": np.__version__, "scipy_version": scipy.__version__,
        "internal_p_scale": P_SCALE, "internal_d_reference": D_REF,
        "primary_fits": primary_fits,
    }
    (RUN / "local_fits.json").write_text(json.dumps(output, indent=2) + "\n")
    write_csv(RUN / "window_sensitivity.csv", window_rows)
    write_csv(RUN / "distance_sensitivity.csv", distance_rows)
    write_csv(RUN / "bootstrap_samples.csv", bootstrap_rows)
    print("Saved primary fits, all 16 window/model fits, distance checks, and bootstrap samples.")


if __name__ == "__main__":
    main()
