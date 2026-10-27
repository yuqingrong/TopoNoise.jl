"""Independent saved-output arithmetic audit; no analysis/model-helper imports.

Run only after analysis reports completion. Reconstructs likelihoods and saved
profile evaluations, regenerates bootstrap count arrays, and checks summaries.
Does not optimize, simulate circuits, or rerun bootstrap/profile fits.
"""
import argparse
from collections import Counter
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.special import expit, gammaln, xlogy
from scipy.stats import binomtest, chi2

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
ROOT = STUDY.parent.parent
OMEGA = {"corr_w0p5": .5, "corr_w1": 1., "corr_w2": 2., "corr_linear_w1": 1.}
CHECKS = 0


class AuditError(ValueError):
    pass


def require(value, message):
    global CHECKS
    CHECKS += 1
    if not value:
        raise AuditError(message)


def close(actual, expected, message, tolerance=2e-8):
    if expected is None:
        require(actual is None, message)
    else:
        a, b = np.asarray(actual), np.asarray(expected)
        require(a.shape == b.shape and np.allclose(a, b, atol=tolerance, rtol=0), message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def key(row):
    return row["channel"], int(row.get("steps", row.get("size"))), float(row["p"])


def likelihood_stats(counts, shots, eta):
    probability = expit(eta); safe = np.clip(probability, 1e-14, 1-1e-14)
    observed = counts / shots
    nll = float(np.sum(shots * np.logaddexp(0., eta) - counts * eta))
    ll = float(np.sum(gammaln(shots+1)-gammaln(counts+1)-gammaln(shots-counts+1))) - nll
    cell_deviance = 2 * (xlogy(counts, observed/safe) + xlogy(shots-counts, (1-observed)/(1-safe)))
    pearson = (observed-probability) / np.sqrt(np.maximum(probability*(1-probability)/shots, 1e-30))
    return dict(nll=nll, log_likelihood=ll, deviance=float(np.sum(cell_deviance)), probability=probability,
                observed=observed, pearson=pearson,
                deviance_residual=np.sign(observed-probability)*np.sqrt(np.maximum(0., cell_deviance)))


def parameter_spec(model):
    require(model in {"quadratic", "cubic", *OMEGA}, "unknown model")
    names = ["pc", "nu", "a0", "a1", "a2"]
    bounds = [(0., 1.), (math.log(.3), math.log(5.)), (-15., 15.), (-30., 30.), (-30., 30.)]
    if model == "cubic": names += ["a3"]; bounds += [(-30., 30.)]
    elif model == "corr_linear_w1": names += ["c0", "c1"]; bounds += [(-15., 15.)]*2
    elif model in OMEGA: names += ["c"]; bounds += [(-15., 15.)]
    return names, bounds


def evaluate_parameters(theta, sizes, rates, limits, reference, model):
    low, high = limits
    x = ((rates-low)/(high-low)-theta[0]) * (sizes/reference)**(1/math.exp(theta[1]))
    degree = 3 if model == "cubic" else 2
    eta = sum(theta[2+j]*x**j for j in range(degree+1))
    if model in OMEGA:
        correction = theta[5] + (theta[6]*x if model == "corr_linear_w1" else 0.)
        eta += correction * (reference/sizes)**OMEGA[model]
    return eta


def parameter_diagnostics(theta, model, sizes, limits, reference):
    names, bounds = parameter_spec(model)
    require(len(theta) == len(bounds), "parameter count differs")
    require(all(math.isfinite(v) and lo-1e-10 <= v <= hi+1e-10 for v, (lo, hi) in zip(theta, bounds)), "parameter outside bounds")
    hits = [name for name, v, (lo, hi) in zip(names, theta, bounds) if min(v-lo, hi-v) <= .002*(hi-lo)]
    minimums = []
    for size in sorted(set(sizes)):
        stretch = (size/reference)**(1/math.exp(theta[1]))
        left, right = -theta[0]*stretch, (1-theta[0])*stretch
        a, b = theta[3], 2*theta[4]
        c = 3*theta[5] if model == "cubic" else 0.
        if model == "corr_linear_w1": a += theta[6]*reference/size
        candidates = [left, right]
        if c != 0 and left <= -b/(2*c) <= right: candidates.append(-b/(2*c))
        minimums.append(min(a+b*x+c*x*x for x in candidates))
    return dict(bound_hits=hits, monotone_increasing=bool(min(minimums)>1e-7), derivative_min=float(min(minimums)))


def audit_fit(fit, rows, limits):
    label = fit["fit_id"]; model = fit["model"]; theta = fit["parameters"]
    sizes = np.array([r["steps"] for r in rows], float); rates = np.array([r["p"] for r in rows])
    shots = np.array([r["shots"] for r in rows], float); counts = np.array([r["failures"] for r in rows], float)
    reference = float(np.exp(np.mean(np.log(np.unique(sizes)))))
    names, bounds = parameter_spec(model)
    close(fit["L_reference"], reference, label+" reference"); close([fit["p_min"], fit["p_max"]], limits, label+" bounds")
    close(fit["p_scale"], limits[1]-limits[0], label+" scale")
    close(fit["pc"], limits[0]+theta[0]*(limits[1]-limits[0]), label+" physical pc")
    close(fit["nu"], math.exp(theta[1]), label+" physical nu")
    require(fit["parameter_names"] == names, label+" parameter names")
    close(fit["parameter_bounds"], bounds, label+" parameter bounds")
    close(fit["coefficients"], theta[2:6 if model == "cubic" else 5], label+" coefficients")
    require(fit["correction_omega"] == OMEGA.get(model), label+" correction exponent")
    close(fit["correction_coefficients"], theta[5:] if model in OMEGA else [], label+" correction coefficients")
    diagnostics = parameter_diagnostics(theta, model, sizes, limits, reference)
    for field in ("bound_hits", "monotone_increasing"):
        require(fit[field] == diagnostics[field], label+" "+field)
    close(fit["derivative_min"], diagnostics["derivative_min"], label+" monotone derivative")
    stats = likelihood_stats(counts, shots, evaluate_parameters(theta, sizes, rates, limits, reference, model))
    for field in ("nll", "log_likelihood", "deviance"): close(fit[field], stats[field], label+" "+field, 1e-7)
    close(fit["fitted_probabilities"], stats["probability"], label+" probabilities")
    close(fit["standardized_residuals"], stats["pearson"], label+" residuals")
    close(fit["largest_absolute_residual"], float(np.max(np.abs(stats["pearson"]))), label+" largest residual")
    cells, parameters = len(rows), len(theta); df = cells-parameters
    require(fit["n_cells"] == cells and fit["n_parameters"] == parameters and fit["degrees_of_freedom"] == df, label+" degrees of freedom")
    require(fit["diagnostic_three_sizes_only"] == (len(set(sizes)) == 3), label+" three-size diagnostic")
    close(fit["deviance_p_value"], float(chi2.sf(stats["deviance"], df)) if df>0 else None, label+" goodness")
    aic = -2*stats["log_likelihood"]+2*parameters
    close(fit["aic"], aic, label+" AIC", 2e-7)
    close(fit["aicc"], aic+2*parameters*(parameters+1)/(cells-parameters-1) if cells>parameters+1 else None, label+" AICc", 2e-7)
    require(fit["aicc_sample_size"] == cells, label+" AICc sample size")
    require(len(fit["predictions"]) == cells, label+" prediction count")
    for i, (saved, original) in enumerate(zip(fit["predictions"], rows)):
        require((saved["size"], saved["p"], saved["shots"], saved["failures"]) ==
                (original["steps"], original["p"], original["shots"], original["failures"]), label+" prediction source cell")
        for field, array in (("predicted", "probability"), ("observed", "observed"), ("pearson_residual", "pearson"), ("deviance_residual", "deviance_residual")):
            close(saved[field], stats[array][i], label+" prediction "+field)


def source_points(snapshot, plan):
    fresh = snapshot["fresh_points"]; pooled = snapshot["pooled_points"]
    require(len(fresh) == len(pooled) == 77, "77 fresh and pooled cells")
    require({c: sum(r["shots"] for r in fresh if r["channel"] == c) for c in ("x_only", "z_only")} == {"x_only":106000, "z_only":98000}, "fresh channel totals")
    expected_jobs = {job["name"]: job for job in plan["jobs"]}
    require({Path(s["directory"]).name for s in snapshot["fresh_shards"]} == set(expected_jobs)
            and len(snapshot["fresh_shards"]) == 13, "13 bound production shard inputs")
    old_root = ROOT/"results/toric3d-h6-full-prefix-2026-10-15"
    manifest_path = old_root/"ARTIFACTS.sha256.json"
    require(sha(manifest_path) == plan["old_manifest_sha256"] == snapshot["old_production"]["manifest_sha256"], "frozen historical manifest binding")
    old_manifest = json.loads(manifest_path.read_text())
    raw_lookup = {}; old = []
    for provenance, is_old in [(s, False) for s in snapshot["fresh_shards"]] + [(s, True) for s in snapshot["old_production"]["shards"]]:
        directory = Path(provenance["directory"])
        for filename, bound in provenance["files_sha256"].items():
            require(sha(directory/filename) == bound, "source shard file binding")
            if is_old: require(bound == old_manifest[f"{directory.name}/{filename}"], "old input file frozen-manifest binding")
        require(sha(directory/provenance["snapshot"]) == provenance["snapshot_sha256"], "source snapshot binding")
        with (directory/"raw.csv").open() as stream:
            for row in csv.DictReader(stream):
                for field in ("steps", "Ly", "Lz", "shots", "target_shots", "seed", "final_only_failures", "full_prefix_failures",
                              "full_prefix_logical_1_failures", "full_prefix_logical_2_failures", "final_only_logical_1_failures", "final_only_logical_2_failures"):
                    row[field] = int(row[field])
                row["p"] = float(row["p"]); row["joint"] = [int(row[f"joint_{i:02d}"]) for i in range(4)]
                row["source_job"] = directory.name
                if is_old: old.append(row)
                else:
                    require(key(row) not in raw_lookup, "duplicate raw fresh cell")
                    raw_lookup[key(row)] = row
    require(len(raw_lookup) == 77 and len(old) == 70 and sum(r["shots"] for r in old) == 70000, "input grids")
    require(len({key(r) for r in fresh}) == 77 and {key(r) for r in fresh} == set(raw_lookup), "fresh snapshot cell coverage")
    seeds = [r["point_seed"] for r in fresh]
    excluded = {r["point_seed"] for r in old}
    for saved in snapshot["excluded_streams"]:
        require(sha(saved["path"]) == saved["sha256"], "excluded stream source binding")
        excluded.update(saved["point_seeds"])
    require(len(set(seeds)) == 77 and not (set(seeds)&excluded), "new stream reuse exclusion")
    fields = ("shots", "target_shots", "seed", "point_seed", "source_job", "full_prefix_failures", "final_only_failures",
              "full_prefix_logical_1_failures", "full_prefix_logical_2_failures", "final_only_logical_1_failures", "final_only_logical_2_failures", "joint")
    for row in fresh:
        require(key(row) in raw_lookup, "unknown fresh cell")
        raw = raw_lookup[key(row)]
        require(all(row[f] == raw[f] for f in fields), "fresh snapshot differs from committed raw counts")
        require(row["failures"] == row["none_failures"] == row["full_prefix_failures"], "fresh failure aliases")
        job = expected_jobs[row["source_job"]]
        require(row["shots"] == job["shots"] and row["seed"] == job["seed"] and row["steps"] in job["sizes"] and row["p"] in job["rates"], "fresh frozen design")
    old_lookup = {key(row): row for row in old}
    matched = [row for row in fresh if key(row) in old_lookup]
    require(len(matched) == 36 and sum(old_lookup[key(r)]["shots"] for r in matched) == 36000, "old pool 36 cells/36k records")
    require(len({key(r) for r in pooled}) == 77 and sum(r["shots"] for r in pooled) == 240000, "pooled counts")
    fresh_lookup = {key(r):r for r in fresh}
    sums = [f for f in fields if "failures" in f]
    for row in pooled:
        r = fresh_lookup[key(row)]; prior = old_lookup.get(key(row)); nold = prior["shots"] if prior else 0
        require(row["shots"] == r["shots"]+nold and row["fresh_shots"] == r["shots"] and row["old_shots"] == nold, "pooled sample counts")
        for field in sums: require(row[field] == r[field]+(prior[field] if prior else 0), "pooled failure count")
        require(row["joint"] == [a+(b if prior else 0) for a,b in zip(r["joint"], prior["joint"] if prior else [0]*4)], "pooled histogram")
        require(row["failures"] == row["none_failures"] == row["full_prefix_failures"], "pooled failure aliases")
        require([s["point_seed"] for s in row["stream_components"]] == [r["point_seed"]]+([prior["point_seed"]] if prior else []), "pooled stream provenance")
    require(snapshot["validation"]["pooling"]["old_records_added"] == 36000, "pool summary")
    return fresh, pooled


def exact_interval(k, n, alpha):
    ci = binomtest(k, n).proportion_ci(1-alpha, method="exact")
    return ci.low, ci.high


def audit_crossings(x, points, plan):
    lookup = {(r["steps"], r["p"]):r for r in points}
    anchor = x["primary_crossing"]; spec = plan["analysis"]["x"]["primary_crossing"]
    wanted = {(L,p) for L in spec["size_pair"] for p in spec["rates"]}
    require({(r["L"],r["p"]) for r in anchor["cells"]} == wanted and len(anchor["cells"]) == 4, "four anchor cells")
    close(anchor["cell_alpha"], .05/4, "anchor multiplicity")
    bands = {}
    for saved in anchor["cells"]:
        row = lookup[saved["L"],saved["p"]]; interval = exact_interval(row["failures"],row["shots"],.05/4)
        require((saved["shots"],saved["failures"]) == (row["shots"],row["failures"]), "anchor counts")
        close([saved["ci_low"],saved["ci_high"]],interval,"exact anchor interval",5e-12)
        bands[saved["L"],saved["p"]] = interval
    expected_signs = []
    for saved,p in zip(anchor["contrasts"],spec["rates"]):
        small,large = spec["size_pair"]; a,b = bands[small,p],bands[large,p]
        close([saved["ci_low"],saved["ci_high"]],[b[0]-a[1],b[1]-a[0]],"anchor contrast")
        close(saved["difference"],lookup[large,p]["failures"]/lookup[large,p]["shots"]-lookup[small,p]["failures"]/lookup[small,p]["shots"],"anchor difference")
        expected_signs.append((b[0]-a[1],b[1]-a[0]))
    require(anchor["success"] == (expected_signs[0][1]<0 and expected_signs[1][0]>0), "anchor success")
    secondary = x["secondary_crossing"]; bands = {}
    require(len(secondary["cells"]) == secondary["number_cells"] == 28, "28 secondary cells")
    require({(r["L"],r["p"]) for r in secondary["cells"]} == set(lookup), "secondary cell coverage")
    close(secondary["cell_alpha"],.05/28,"secondary multiplicity")
    for saved in secondary["cells"]:
        row = lookup[saved["L"],saved["p"]]; interval = exact_interval(row["failures"],row["shots"],.05/28)
        require((saved["shots"],saved["failures"]) == (row["shots"],row["failures"]), "secondary counts")
        close([saved["ci_low"],saved["ci_high"]],interval,"secondary exact CI",5e-12)
        bands[saved["L"],saved["p"]] = interval
    require(len(secondary["bands"]) == 42, "42 pair/rate difference bands")
    expected_bands = {(a,b,p) for a,b in itertools.combinations(sorted({r["steps"] for r in points}),2)
                      for p in sorted({r["p"] for r in points})}
    require({(r["L_small"],r["L_large"],r["p"]) for r in secondary["bands"]} == expected_bands, "secondary pair/rate coverage")
    for saved in secondary["bands"]:
        p,small,large = saved["p"],saved["L_small"],saved["L_large"]
        a,b = bands[small,p],bands[large,p]
        close([saved["ci_low"],saved["ci_high"]],[b[0]-a[1],b[1]-a[0]],"secondary difference band")
        ra,rb = lookup[small,p],lookup[large,p]; qa,qb = ra["failures"]/ra["shots"],rb["failures"]/rb["shots"]
        close(saved["difference"],qb-qa,"secondary difference")
        se = math.sqrt(qa*(1-qa)/ra["shots"]+qb*(1-qb)/rb["shots"])
        close(saved["standard_error"],se,"secondary standard error")


def audit_bootstrap(x, points, plan):
    if "bootstrap" not in x:
        require(any(f["stage"] == "bootstrap" for f in x["failures"]), "missing unreported bootstrap")
        return {"available":False,"declared_failure":True}
    common, alternative, boot = x["common_model"],x["independent_size_model"],x["bootstrap"]
    sizes = np.array([r["steps"] for r in points],float); rates = np.array([r["p"] for r in points])
    n = np.array([r["shots"] for r in points],int); k = np.array([r["failures"] for r in points],int)
    coord = (rates-common["pc"])/common["p_scale"]*(sizes/common["L_reference"])**(1/common["nu"])
    eta = sum(c*coord**j for j,c in enumerate(common["coefficients"])); stats = likelihood_stats(k,n,eta)
    close(common["fitted_probabilities"],stats["probability"],"bootstrap null predictions")
    close(common["deviance"],stats["deviance"],"bootstrap null deviance",1e-7)
    curves = {r["L"]:r for r in alternative["curves"]}; t = 2*(rates-rates.min())/(rates.max()-rates.min())-1
    alteta = np.array([sum(c*t[i]**j for j,c in enumerate(curves[int(L)]["coefficients"])) for i,L in enumerate(sizes)])
    altstats = likelihood_stats(k,n,alteta)
    close(alternative["fitted_probabilities"],altstats["probability"],"independent size predictions")
    close(alternative["deviance"],altstats["deviance"],"independent size deviance",1e-7)
    observed_lr = max(0.,common["deviance"]-alternative["deviance"])
    close(x["observed_likelihood_ratio"],observed_lr,"observed likelihood ratio")
    repetitions,seed = plan["analysis"]["x"]["bootstrap_replicates"],plan["analysis"]["x"]["bootstrap_seed"]
    records = boot["records"]
    require(len(records) == boot["replicates"] == repetitions == 1000 and boot["seed"] == seed,"bootstrap frozen repetition/seed")
    require([r["replicate"] for r in records] == list(range(repetitions)),"bootstrap index sequence")
    draws = np.random.default_rng(seed).binomial(n,np.array(common["fitted_probabilities"]),size=(repetitions,len(points)))
    for i,row in enumerate(records):
        require(row["simulated_failures"] == draws[i].tolist(),"regenerated bootstrap count set")
        if row["common_valid"] and row["alternative_valid"]:
            close(row["raw_likelihood_ratio"],row["common_deviance"]-row["alternative_deviance"],"bootstrap raw likelihood ratio")
            close(row["likelihood_ratio"],max(0.,row["raw_likelihood_ratio"]),"bootstrap truncated likelihood ratio")
    valid = [r for r in records if r["common_valid"]]; joint = [r for r in valid if r["alternative_valid"]]
    for name,rows,field,threshold in (("goodness_of_fit",valid,"common_deviance",common["deviance"]),
                                     ("likelihood_ratio",joint,"likelihood_ratio",observed_lr),
                                     ("old_raw_crossing_rule",valid,"old_crossing_rejected",True)):
        tail = sum(r[field] >= threshold-1e-10 for r in rows); summary = boot[name]; nv = len(rows)
        require(summary["tail_count"] == tail and summary["valid_replicates"] == nv and summary["attempted_replicates"] == repetitions,"bootstrap tail counts")
        require(summary["numerical_failures"] == repetitions-nv,"bootstrap failed-replicate denominator")
        close(summary["p_plus_one"],(tail+1)/(nv+1) if nv else None,"bootstrap tail probability")
        close(summary["mc_binomial_interval"],exact_interval(tail,nv,.05) if nv else None,"bootstrap MC CI",5e-12)
        close(summary["mc_standard_error"],math.sqrt((tail/nv)*(1-tail/nv)/nv) if nv else None,"bootstrap MC SE")
        close(summary["all_attempt_failure_bounds"],[(tail+1)/(repetitions+1),(tail+repetitions-nv+1)/(repetitions+1)],"bootstrap failure bounds")
    for field in ("pc","nu"):
        close(boot["conditional_percentile_intervals"][field],np.quantile([r[field] for r in valid],[.025,.975]) if len(valid)>=2 else None,"bootstrap percentile "+field)
    require(boot["bound_hit_replicates"] == sum(bool(r["bound_hits"]) for r in valid),"bootstrap boundary flags")
    require(boot["nonmonotone_replicates"] == sum(not r["monotone"] for r in valid),"bootstrap monotonicity flags")
    return {"available":True,"regenerated_count_sets":repetitions,"valid_common":len(valid),"valid_common_and_alternative":len(joint)}


def validate_regions(profile, successful):
    cutoff = float(chi2.ppf(.95,1)); close(profile["lr_cutoff"],cutoff,"profile LR cutoff")
    last = -math.inf
    for region in profile["regions"]:
        require(last <= region["lower"] <= region["upper"],"profile region order"); last = region["upper"]
        for side,bound in (("lower",profile["pc_bounds"][0]),("upper",profile["pc_bounds"][1])):
            if region[side+"_bound_truncated"]:
                close(region[side],bound,"truncated profile boundary")
            elif not region[side+"_numerical_gap"]:
                require(bool(successful),"profile endpoint lacks saved successful evaluation")
                point = min(successful,key=lambda r:abs(r["pc"]-region[side]))
                close(point["pc"],region[side],"saved refined profile endpoint",2e-12)
                close(point["lr"],cutoff,"profile endpoint LR cutoff",2e-5)


def audit_profile(profile,fit):
    if profile.get("critical_helper_exception"):
        require(profile["minimum_agreement"] is False and not profile["regions"] and profile["numerical_failures"],"critical profile failure must emit no interval")
        return {"fit_id":profile["fit_id"],"critical_helper_exception":True,"successful_evaluations":0,"numerically_invalid":True}
    rows = fit["predictions"]; sizes = np.array([r["size"] for r in rows],float); rates = np.array([r["p"] for r in rows])
    n = np.array([r["shots"] for r in rows],float); k = np.array([r["failures"] for r in rows],float)
    limits = [fit["p_min"],fit["p_max"]]
    close(profile["pc_bounds"],limits,"profile original bounds"); close(profile["free_nll"],fit["nll"],"profile free likelihood")
    close(profile["free_pc"],fit["pc"],"profile free pc")
    require(profile["window"] == fit["window"] and profile["model"] == fit["model"],"profile model/window labels")
    successful = [r for r in profile["evaluations"] if r["optimizer_success"]]
    require(successful,"noncritical profile has no successful evaluation")
    for row in successful:
        theta = row["parameters"]
        close(theta[0],(row["pc"]-limits[0])/(limits[1]-limits[0]),"profile fixed parameter")
        diagnostics = parameter_diagnostics(theta,fit["model"],sizes,limits,fit["L_reference"])
        require(row["bound_hits"] == [name for name in diagnostics["bound_hits"] if name != "pc"],"profile nuisance boundary flags")
        stats = likelihood_stats(k,n,evaluate_parameters(theta,sizes,rates,limits,fit["L_reference"],fit["model"]))
        close(row["nll"],stats["nll"],"profile independently reconstructed NLL",1e-7)
        close(row["lr"],2*(row["nll"]-fit["nll"]),"profile total LR",1e-7)
        require(row["lr"] >= -4e-4-1e-8,"profile improves free minimum materially")
    minimum = min(successful,key=lambda r:r["nll"]); free = profile["profile_at_free_pc"]
    close(profile["minimum_profile_nll"],minimum["nll"],"minimum profile NLL")
    close(profile["minimum_profile_pc"],minimum["pc"],"minimum profile pc")
    close(profile["minimum_nll_difference"],minimum["nll"]-fit["nll"],"profile minimum difference")
    agree = abs(free["nll"]-fit["nll"])<=2e-4 and minimum["nll"]>=fit["nll"]-2e-4
    require(profile["minimum_agreement"] == agree,"profile agreement flag")
    require(profile["nuisance_bound_evaluation_count"] == sum(bool(r.get("bound_hits")) for r in profile["evaluations"]),"profile nuisance bound count")
    require(profile["polish_failure_count"] == sum(r.get("polish_success") is False for r in profile["evaluations"]),"profile polish count")
    validate_regions(profile,successful)
    invalid = bool(profile["numerical_failures"] or not agree or len(successful)!=len(profile["evaluations"]) or
                   any(r.get("lower_numerical_gap") or r.get("upper_numerical_gap") or r.get("refinement_failures") for r in profile["regions"]))
    return {"fit_id":profile["fit_id"],"successful_evaluations":len(successful),"numerically_invalid":invalid}


def audit(output):
    report_path,input_path = output/"report.json",output/"input.json"
    report,snapshot = json.loads(report_path.read_text()),json.loads(input_path.read_text())
    require(report["status"] in ("complete","complete_with_numerical_failures") and "completed_at_unix" in report,"analysis is not complete; no outcome audit performed")
    plan = json.loads((STUDY/"selection.json").read_text()); plan_sha = sha(STUDY/"selection.json")
    require(plan_sha == (STUDY/"selection.sha256").read_text().split()[0] == report["provenance"]["selection_sha256"] == snapshot["selection_sha256"],"frozen selection binding")
    require(sha(input_path) == report["provenance"]["input_sha256"],"analysis input binding")
    require(plan["frozen_before_production_at_unix"] < report["created_at_unix"] <= report["completed_at_unix"],"analysis chronology")
    require(report["provenance"]["simulation_sources_sha256"] == plan["simulation_source_sha256"],"simulation source set")
    for mapping in (report["provenance"]["simulation_sources_sha256"],report["provenance"]["numerical_sources_sha256"]):
        for name,bound in mapping.items(): require(sha(ROOT/name) == bound,"source hash "+name)
    fresh,pooled = source_points(snapshot,plan)
    require(snapshot["validation"] == report["validation"],"input/report validation equality")
    fitted_count = 0; exceptions = 0; bad_fits = 0; profile_checks = []
    for dataset,points in (("fresh",fresh),("pooled_secondary",pooled)):
        for channel,label in (("x_only","x"),("z_only","z")):
            group = report[dataset][label]; cp = plan["analysis"][label]
            require(group["points"] == [r for r in points if r["channel"] == channel],"report point dataset binding")
            expected = {f"Lmin{m}-{w}-{model}":(m,w,model) for m,w,model in itertools.product(cp["min_sizes"],cp["windows"],cp["models"])}
            fits = {f["fit_id"]:f for f in group["fits"]}; failed = {f["fit_id"] for f in group["failures"] if f["stage"] == "fit"}
            require(len(fits) == len(group["fits"]) and not (set(fits)&failed) and set(fits)|failed == set(expected),"all planned fit outcomes retained")
            require(group["fit_count"] == len(fits) and group["expected_fit_count"] == len(expected),"fit count summary")
            exceptions += len(group["failures"])
            for identity,fit in fits.items():
                minimum,window,model = expected[identity]; limits = cp["windows"][window]
                require((fit["min_size"],fit["window"],fit["model"],fit["channel"],fit["dataset"]) == (minimum,window,model,channel,dataset),"fit identifiers")
                rows = sorted([r for r in points if r["channel"] == channel and r["steps"]>=minimum and limits[0]<=r["p"]<=limits[1]],key=key)
                audit_fit(fit,rows,limits); fitted_count += 1; bad_fits += fit["optimizer_success"] is not True
                same = [f["aicc"] for f in fits.values() if f["min_size"] == minimum and f["window"] == window and f["aicc"] is not None]
                close(fit["delta_aicc_same_subset"],fit["aicc"]-min(same) if same and fit["aicc"] is not None else None,"within-subset AICc",2e-7)
    x = report["fresh"]["x"]; xp = [r for r in fresh if r["channel"] == "x_only"]
    audit_crossings(x,xp,plan); bootstrap = audit_bootstrap(x,xp,plan)
    z = report["fresh"]["z"]; cp = plan["analysis"]["z"]
    wanted = {f"Lmin{m}-{w}-{model}" for m,w,model in itertools.product(cp["min_sizes"],cp["windows"],cp["profile_models"])}
    profiles = {p["fit_id"]:p for p in z["profiles"]}; failed = {f["fit_id"] for f in z["failures"] if f["stage"] == "profile"}
    require(len(profiles) == len(z["profiles"]) and not (set(profiles)&failed) and set(profiles)|failed == wanted and len(wanted) == 12,"all 12 planned profile outcomes")
    fits = {f["fit_id"]:f for f in z["fits"]}
    for identity,profile in profiles.items(): profile_checks.append(audit_profile(profile,fits[identity]))
    require(z["profile_count"] == len(profiles) and z["expected_profile_count"] == 12,"profile count summary")
    numerical_failures = bool(exceptions or bad_fits or any(p["numerically_invalid"] for p in profile_checks) or not bootstrap["available"] or
                              bootstrap.get("valid_common_and_alternative",0) != 1000)
    require(report["status"] == ("complete_with_numerical_failures" if numerical_failures else "complete"),"independent numerical completion classification")
    require(report["numerical_summary"]["has_numerical_failures"] == numerical_failures,"numerical summary failure flag")
    return {"status":"passed", "checked_at_unix":time.time(), "assertions":CHECKS,
            "scope":"Independent source/input bindings and saved arithmetic only; no helper imports, optimization, circuit simulation or full bootstrap/profile reruns.",
            "input_sha256":sha(input_path),"report_sha256":sha(report_path),"checker_sha256":sha(Path(__file__)),"selection_sha256":plan_sha,
            "fresh_records":204000,"pooled_records":240000,"old_records_pooled":36000,"planned_fits":120,"returned_fits_verified":fitted_count,
            "declared_numerical_exceptions":exceptions,"optimizer_failure_count":bad_fits,"analysis_numerical_failures_retained":numerical_failures,
            "x_primary_anchor_cells":4,"x_secondary_cells":28,"x_secondary_difference_bands":42,
            "bootstrap":bootstrap,"planned_profiles":12,"profiles":profile_checks}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis",type=Path,default=STUDY/"analysis")
    parser.add_argument("--output",type=Path)
    args = parser.parse_args()
    try: result = audit(args.analysis)
    except (AuditError,OSError,KeyError,ValueError,TypeError) as error:
        result = {"status":"failed","checked_at_unix":time.time(),"assertions":CHECKS,"error":str(error)}
    text = json.dumps(result,indent=2,allow_nan=False)+"\n"
    if args.output: args.output.write_text(text)
    print(text,end="")
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__": raise SystemExit(main())
