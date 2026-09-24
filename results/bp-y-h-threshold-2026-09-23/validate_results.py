"""Verify provenance, every count/error bar, fit statistics, and figure exports."""
import csv
import hashlib
import json
from pathlib import Path
import tomllib

import numpy as np
from PIL import Image

from analysis import data_for, load_points, model
from run_experiment import production_specs, provenance
from sampler import atomic_json

RUN = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def close(a, b):
    return bool(np.allclose(a, b, rtol=1e-9, atol=1e-12))


def verify_fit(fit, channel, component, points):
    lower, upper = fit["window"]
    selected = [p for p in points if p["stage"] == "production" and p["channel"] == channel
                and p["distance"] in fit["distances"] and lower-1e-12 <= p["p"] <= upper+1e-12]
    require(len(selected) == fit["points"], "Fit point selection differs from saved metadata")
    require(all(p["shots"] == 50000 for p in selected), "Fit used a non-production shot count")
    data = data_for(selected, component)
    theta = np.array(fit["theta_internal"])
    predicted = model(theta, data["p"], data["d"])
    bounded = np.clip(data["y"], 0.5/data["n"], 1-0.5/data["n"])
    objective = np.sum((predicted-data["y"])**2/(bounded*(1-bounded)/data["n"]))
    require(close(objective, fit["chi_square"]), "Saved chi-square differs from raw counts")
    x = (data["p"]-fit["p_c"]) * data["d"]**(1/fit["nu"])
    require(close(predicted, np.polynomial.polynomial.polyval(x, fit["coefficients"])),
            "Physical and internally rescaled fit coefficients disagree")
    require(lower < fit["p_c"] < upper and 0.3 < fit["nu"] < 5,
            "Fit parameter at a search boundary")


def main():
    meta = json.loads((RUN / "metadata.json").read_text())
    require(meta["sampling_complete"], "Sampling has not completed")
    require(meta["fingerprint"] == provenance()["fingerprint"], "Sampler/environment fingerprint changed")
    sources = tomllib.loads((RUN / "source_hashes.toml").read_text())["sha256"]
    for name, expected in sources.items():
        require(sha(RUN / "source_snapshot" / name) == expected, f"Archived source changed: {name}")
    points = load_points()
    keys = [(r["stage"], r["channel"], r["distance"], r["p"]) for r in points]
    require(len(keys) == len(set(keys)), "Duplicate simulation point")
    for r in points:
        require(sum(r["counts"]) == r["shots"] and min(r["counts"]) >= 0, "Invalid joint counts")
        require(r["fingerprint"] == meta["fingerprint"], "Mixed sampler fingerprints")
        if r["p"] == 0:
            require(r["counts"] == [r["shots"], 0, 0, 0], "Nonzero noiseless logical failure")
    require(meta["saved_points"] == len(points), "Metadata point count mismatch")
    require(meta["saved_shots"] == sum(r["shots"] for r in points), "Metadata shot count mismatch")
    selection = json.loads((RUN / "selection.json").read_text())
    require(selection["frozen"], "Production windows were not frozen")
    pilot = [{k: r[k] for k in ("channel", "distance", "p", "shots", "counts")}
             for r in points if r["stage"] == "pilot"]
    digest = hashlib.sha256(json.dumps(pilot, sort_keys=True).encode()).hexdigest()
    require(digest == selection["pilot_counts_sha256"], "Pilot counts changed after window selection")
    required = {(s, c, d, p): n for s, c, d, p, n in production_specs(selection)}
    production = {k: r for k, r in zip(keys, points) if r["stage"] == "production"}
    require(set(required) == set(production), "Missing or unexpected production grid point")
    require(all(production[k]["shots"] == n for k, n in required.items()), "Incomplete production sampling")
    initial_path = RUN / "selection-initial.json"
    if initial_path.exists():
        initial = json.loads(initial_path.read_text())
        for key, choice in initial["panels"].items():
            if choice["status"] == "selected":
                require(choice == selection["panels"][key], "Previously frozen selected window changed")
    with (RUN / "raw_counts.csv").open(newline="") as stream:
        csv_rows = list(csv.DictReader(stream))
    require(len(csv_rows) == len(points), "Raw CSV point count mismatch")
    point_map = dict(zip(keys, points))
    for row in csv_rows:
        r = point_map[(row["stage"], row["channel"], int(row["distance"]), float(row["p"]))]
        counts = [int(row[k]) for k in ("n00", "n10", "n01", "n11")]
        require(counts == r["counts"] and int(row["shots"]) == r["shots"], "CSV counts differ from checkpoint")
        for component, index in (("x", 1), ("z", 2)):
            failures = counts[index]+counts[3]
            rate = failures/r["shots"]
            require(int(row[f"logical_{component}_failures"]) == failures, "Incorrect marginal count")
            require(close(float(row[f"logical_{component}_rate"]), rate), "Incorrect marginal rate")
            require(close(float(row[f"logical_{component}_standard_error"]),
                          np.sqrt(rate*(1-rate)/r["shots"])), "Incorrect binomial error bar")
    fits = json.loads((RUN / "local_fits.json").read_text())
    require(fits["raw_sha256"] == sha(RUN / "raw_counts.csv"), "Raw data changed since fitting")
    require(fits["script_sha256"] == sha(RUN / "analysis.py"), "Fitting script changed since fitting")
    require(fits["selection_sha256"] == sha(RUN / "selection.json"), "Fit selection hash mismatch")
    with (RUN / "bootstrap_samples.csv").open(newline="") as stream:
        bootstrap = list(csv.DictReader(stream))
    for key, fit in fits["primary_fits"].items():
        if "theta_internal" not in fit:
            require(fit["status"] in ("unresolved", "invalid"), "Unavailable fit carries an invalid status")
            continue
        channel, component = key.split("/")
        verify_fit(fit, channel, component, points)
        samples = [s for s in bootstrap if s["panel"] == key]
        require(len(samples) == fit["bootstrap_successes"], "Bootstrap count mismatch")
        if len(samples) >= 2:
            for parameter in ("p_c", "nu"):
                values = np.array([float(s[parameter]) for s in samples])
                require(close(values.std(ddof=1), fit[parameter+"_standard_error"]), "Bootstrap SE mismatch")
                require(close(np.quantile(values, [0.025, 0.975]), fit[parameter+"_interval_95"]),
                        "Bootstrap confidence interval mismatch")
        if fit["status"] == "success":
            require(len(samples) >= 950 and fit["bootstrap_replicates"] == 1000 and
                    fit["goodness_of_fit_p"] >= 0.01, "Unsupported successful-fit designation")
    sensitivity = json.loads((RUN / "sensitivity.json").read_text())
    for fit in sensitivity:
        if "theta_internal" in fit:
            verify_fit(fit, *fit["panel"].split("/"), points)
    plot = json.loads((RUN / "plot_data.json").read_text())
    require(plot["raw_sha256"] == sha(RUN / "raw_counts.csv"), "Plot used different raw data")
    require(plot["fit_sha256"] == sha(RUN / "local_fits.json"), "Plot used different fits")
    plotted = 0
    for series in plot["series"]:
        channel, component = series["panel"].split("/")
        selected = sorted([r for r in csv_rows if r["channel"] == channel and
                           r["stage"] == series["stage"] and int(r["distance"]) == series["distance"]],
                          key=lambda r: float(r["p"]))
        require(close(series["p"], [float(r["p"]) for r in selected]), "Plot p values mismatch")
        for name, field in (("rate", "rate"), ("standard_error", "standard_error")):
            require(close(series[name], [float(r[f"logical_{component}_{field}"]) for r in selected]),
                    "Plotted rate/error bar mismatch")
        plotted += len(selected)
    require(plotted == 2*len(points), "A marginal data series is missing from the plot")
    figures = []
    for name in ("logical_errors", "near_critical_scaling"):
        for ext in ("png", "pdf", "svg"):
            path = RUN / f"{name}.{ext}"
            require(path.stat().st_size > 1000, "Figure export is missing or empty")
            if ext == "png":
                with Image.open(path) as img:
                    img.verify()
            elif ext == "pdf":
                content = path.read_bytes()
                require(content.startswith(b"%PDF") and content.rstrip().endswith(b"%%EOF"), "Invalid PDF framing")
            else:
                require("<svg" in path.read_text(), "Invalid SVG framing")
            figures.append(dict(path=path.name, bytes=path.stat().st_size, sha256=sha(path)))
    tests = json.loads((RUN / "test_status.json").read_text())
    require(tests["python_exit_code"] == 0 and tests["julia_exit_code"] == 0, "Regression tests did not pass")
    require(tests["python_log_sha256"] == sha(RUN / "tests.log") and
            tests["julia_log_sha256"] == sha(RUN / "package-tests-retry.log"), "Test log hash changed")
    result = dict(status="passed", points=len(points), shots=sum(r["shots"] for r in points),
                  plotted_marginals=plotted, source_files=len(sources), sensitivity_fits=len(sensitivity),
                  bootstrap_fits=len(bootstrap), tests=tests, figures=figures,
                  visual_inspection="Recorded separately in README after viewing rendered PNGs")
    atomic_json(RUN / "validation.json", result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
