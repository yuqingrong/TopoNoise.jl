#!/usr/bin/env python3
"""Gather only manifest-planned HPC tasks; never import local-study samples.

JSON collection uses the standard library. --plot additionally uses matplotlib.
Missing planned shards stay in the denominator. Rates and Wilson intervals are
available only after every requested history and independent audit qualifies.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import tempfile
import time
import tomllib


HERE = Path(__file__).resolve().parent


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


HPC = module("hpc_manifest_for_collect", HERE / "manifest.py")
BASE = module("hpc_plot_helpers", HERE.parent / "plot.py")
STAGES = ("zero", "smoke", "feasibility", "production")
RECORD_FIELDS = {"shot_index", "raw_fault_count", "raw_fault_indices", "seed_syndrome", "seed_recovery",
                 "recovery_mask", "final_syndrome", "true_mask", "mwpm_failure", "raw_seed_mask", "transfer_mask"}
RUNTIME_FIELDS = ("julia_version", "ITensors_version", "ITensorMPS_version", "PythonCall_version",
                  "python_version", "python_executable", "pymatching_version", "TopoNoise_version",
                  "TopoNoise_project_path", "active_project_path", "julia_threads", "blas_threads")
MODEL_FIELDS = ("physical_noise", "preparation", "tn_uses_both_records", "tn_conditioning", "mwpm_conditioning",
                "logical_mask_order", "truth_definition", "invalid_policy", "posterior_cache", "rng_algorithm",
                "exact_max_width", "ml_tie_tolerance")


def read_toml(path):
    with Path(path).open("rb") as stream:
        return tomllib.load(stream)


def integer(value, minimum=0):
    return type(value) is int and value >= minimum


def record_identity(record):
    """Runner identity_hash for its flat, sorted bool/int record TOML schema.

    For these fields Julia TOML.print(sorted=true) emits one sorted `key = value`
    line per field; JSON bool/int arrays have the same spelling and spacing.
    Do not extend this to arbitrary TOML documents or floating-point metadata.
    """
    fields = {k:v for k,v in record.items() if k != "record_id"}
    if set(fields) != RECORD_FIELDS:
        raise ValueError("unsupported raw record schema")
    for value in fields.values():
        if type(value) not in (int, bool, list):
            raise ValueError("unsupported raw record identity value")
        if isinstance(value, list) and any(type(v) not in (int, bool) for v in value):
            raise ValueError("unsupported raw record identity array")
    text = ''.join(f"{key} = {json.dumps(value)}\n" for key,value in sorted(fields.items()))
    return hashlib.sha256(text.encode()).hexdigest()


def validate_record(record, index):
    if record_identity(record) != record.get("record_id"):
        raise ValueError(f"raw record {index} identity hash differs")
    if record["shot_index"] != index or not integer(record["shot_index"], 1):
        raise ValueError("raw record indices must be consecutive and unique")
    count, faults = record["raw_fault_count"], record["raw_fault_indices"]
    if (not integer(count) or not isinstance(faults, list) or faults != sorted(set(faults))
            or any(not integer(i,1) or i > count for i in faults)):
        raise ValueError("invalid raw fault indices")
    for name in ("seed_syndrome", "seed_recovery", "final_syndrome"):
        if not isinstance(record[name], list) or any(type(v) is not bool for v in record[name]):
            raise ValueError("syndrome and recovery records must contain Boolean bits")
    for name in ("true_mask", "raw_seed_mask", "recovery_mask", "transfer_mask"):
        if not integer(record[name]) or record[name] > 3:
            raise ValueError("invalid two-bit logical mask in raw record")
    if record["true_mask"] != record["raw_seed_mask"] ^ record["recovery_mask"] ^ record["transfer_mask"]:
        raise ValueError("raw record logical truth does not equal the propagated component XOR")
    if type(record["mwpm_failure"]) is not bool:
        raise ValueError("raw record MWPM outcome must be Boolean")


def validate_plan(manifest):
    outputs, production_seeds = set(), set()
    groups = defaultdict(list)
    if set(manifest["stages"]) != set(STAGES):
        raise ValueError("manifest must contain exactly zero, smoke, feasibility, and production stages")
    for stage in STAGES:
        seen_ids = set()
        for task in HPC.stage_tasks(manifest,stage):
            if task["stage"] != stage or not integer(task["task_id"]) or task["task_id"] in seen_ids:
                raise ValueError("duplicate or invalid planned task identity")
            seen_ids.add(task["task_id"])
            output = Path(task["output"])
            if output.is_absolute() or ".." in output.parts or not output.parts:
                raise ValueError("planned output must be a relative path inside the study")
            if str(output) in outputs:
                raise ValueError("duplicate planned output directory")
            outputs.add(str(output))
            if not integer(task["shots"],1) or not integer(task["seed"]):
                raise ValueError("invalid planned history count or seed")
            if stage == "production":
                if task["seed"] in production_seeds:
                    raise ValueError("duplicate production seed; shards must be independent")
                production_seeds.add(task["seed"])
                if task["p"] <= 0:
                    raise ValueError("zero noise is not statistical production")
                if task["L"] == 3 and task["exact"] is not True:
                    raise ValueError("L3 production requires independent exact references")
                groups[(task["L"],task["p"],task["chi"],task["cutoff"])].append(task)
    for tasks in groups.values():
        if len({task["shard_id"] for task in tasks}) != len(tasks):
            raise ValueError("duplicate planned shard identity")
        if (len(tasks) != manifest["shards_per_point"]
                or any(task["shots"] != manifest["shots_per_shard"] for task in tasks)
                or sum(task["shots"] for task in tasks) != manifest["target_shots_per_point"]):
            raise ValueError("planned production shards do not sum to the frozen point target")


def posterior(document, tolerance):
    probabilities = document.get("probabilities")
    label = document.get("logical_mask")
    if (not isinstance(probabilities,list) or len(probabilities) != 4
            or any(type(v) not in (int,float) or not math.isfinite(v) or v < -tolerance for v in probabilities)
            or abs(sum(probabilities)-1) > 1e-7):
        raise ValueError("invalid normalized four-sector posterior")
    if not integer(label) or label > 3 or max(probabilities)-probabilities[label] > 1e-9:
        raise ValueError("logical decision does not maximize its posterior")
    return probabilities, label


def inspect_task(root, task, manifest):
    root = Path(root).resolve()
    path = (root / task["output"]).resolve()
    if not path.is_relative_to(root):
        raise ValueError("task output escapes its study root")
    N = task["shots"]
    row = dict(task_id=task["task_id"], stage=task["stage"], phase=task["phase"], L=task["L"], p=task["p"],
               chi=task["chi"], cutoff=task["cutoff"], seed=task["seed"], output=task["output"],
               requested_shots=N, identity_valid=False, status="missing", complete=False,
               saved_histories=0, saved_outputs=0, tn_valid_shots=0, validated_shots=0,
               unknown_shots=N, tn_failures=0, exact_checked=0, exact_disagreements=0,
               strict_exact_matches=0, max_posterior_linf=None, record_ids=[],
               outcomes=[], errors=[], shot_errors=[], audit=None, supervision=None)
    if not (path / "metadata.toml").exists():
        return row
    try:
        metadata = read_toml(path / "metadata.toml")
        HPC.validate_task_metadata(metadata,task,manifest)
        if metadata.get("tn_uses_both_records") is not True:
            raise ValueError("metadata does not identify the approved seed-plus-final conditioning model")
        for field in RUNTIME_FIELDS:
            value = metadata.get(field)
            if field in ("julia_threads", "blas_threads"):
                if not integer(value,1) or value != manifest[field]:
                    raise ValueError(f"runtime {field} differs from the frozen thread count")
            elif not isinstance(value,str) or not value:
                raise ValueError(f"missing runtime identity: {field}")
        run_id = metadata.get("run_id")
        if not isinstance(run_id,str) or not run_id:
            raise ValueError("missing run identity")
        corpus = read_toml(path / "input-corpus.toml")
        if (corpus.get("run_id") != run_id or corpus.get("requested_shots") != N
                or corpus.get("rng_seed") != task["seed"] or len(corpus.get("records",[])) != N):
            raise ValueError("raw corpus identity, seed, or count differs from the planned shard")
        records = corpus["records"]
        for i,record in enumerate(records,1):
            validate_record(record,i)
        expected_filenames = {f"shot_{i:04d}.toml" for i in range(1,N+1)}
        if any(p.name not in expected_filenames for p in (path / "shots").glob("shot_*.toml")):
            raise ValueError("duplicate or unexpected shot output filename")
        row.update(identity_valid=True,run_id=run_id,saved_histories=N,
                   record_ids=[r["record_id"] for r in records],
                   runtime_model_identity={field:metadata.get(field) for field in RUNTIME_FIELDS+MODEL_FIELDS})
    except (OSError,ValueError,KeyError,TypeError) as error:
        row.update(status="invalid_identity",errors=[str(error)])
        return row

    deviations = []
    for i,record in enumerate(records,1):
        output = path / "shots" / f"shot_{i:04d}.toml"
        if output.exists():
            row["saved_outputs"] += 1
        try:
            shot = read_toml(output)
            expected = dict(run_id=run_id,record_id=record["record_id"],shot_index=i,
                            L=task["L"],p=task["p"],chi=task["chi"],cutoff=task["cutoff"],
                            true_mask=record["true_mask"],mwpm_failure=record["mwpm_failure"])
            if any(shot.get(k) != value for k,value in expected.items()):
                raise ValueError("shot identity, configuration, or raw-truth binding differs")
            if shot.get("status") != "complete" or shot.get("valid") is not True:
                raise ValueError("TN output is invalid or incomplete")
            probabilities,label = posterior(shot,manifest["negative_tolerance"])
            if type(shot.get("TNfailure")) is not bool or shot["TNfailure"] != (label != record["true_mask"]):
                raise ValueError("TN failure score differs from raw logical truth")
            row["tn_valid_shots"] += 1
            if task["stage"] == "zero" and (record["raw_fault_indices"] or record["true_mask"] != 0
                    or max(abs(value-(j==0)) for j,value in enumerate(probabilities)) > 1e-8):
                raise ValueError("zero-noise history or posterior differs from the analytical prediction")
            if task["exact"]:
                exact = read_toml(path / "exact" / f"shot_{i:04d}.toml")
                if (exact.get("status") != "complete" or exact.get("valid") is not True
                        or any(exact.get(k) != value for k,value in expected.items() if k != "mwpm_failure")):
                    raise ValueError("independent exact reference identity/configuration is invalid")
                reference,exact_label = posterior(exact,manifest["negative_tolerance"])
                row["exact_checked"] += 1
                row["strict_exact_matches"] += label == exact_label
                deviations.append(max(abs(a-b) for a,b in zip(probabilities,reference)))
                if max(reference)-reference[label] > 1e-9:
                    row["exact_disagreements"] += 1
                    raise ValueError("TN decision disagrees with independent exact ML")
            row["validated_shots"] += 1
            row["tn_failures"] += shot["TNfailure"]
            row["outcomes"].append(dict(shot_index=i,record_id=record["record_id"],logical_mask=label,
                                        probabilities=probabilities,elapsed_seconds=shot.get("elapsed_seconds")))
        except (OSError,ValueError,KeyError,TypeError) as error:
            row["shot_errors"].append(dict(shot_index=i,reason=str(error)))
    row["unknown_shots"] = N-row["validated_shots"]
    row["max_posterior_linf"] = max(deviations,default=None)
    try:
        supervision_path = path / "supervision.json"
        if supervision_path.exists():
            row["supervision"] = json.loads(supervision_path.read_text())
        row["audit"] = HPC.audit_task(root,task,manifest)
    except (OSError,ValueError,KeyError,TypeError) as error:
        row["errors"].append(str(error))
    process_ok = (isinstance(row["supervision"],dict) and row["supervision"].get("status") == "completed"
                  and type(row["supervision"].get("exit_code")) is int and row["supervision"]["exit_code"] == 0)
    row["complete"] = (row["unknown_shots"] == 0 and not row["errors"] and process_ok
                       and bool(row["audit"] and row["audit"].get("passed")))
    row["status"] = "complete" if row["complete"] else ("blocked" if row["supervision"] and not process_ok else "incomplete")
    return row


def totals(rows):
    return {field:sum(row[field] for row in rows) for field in (
        "requested_shots","saved_histories","saved_outputs","tn_valid_shots","validated_shots",
        "unknown_shots","tn_failures","exact_checked","exact_disagreements")}


def feasibility_pairs(tasks, rows):
    grouped = defaultdict(list)
    lookup = {row["task_id"]:row for row in rows}
    for task in tasks:
        grouped[task["pair_id"]].append(task)
    pairs = []
    for pair_id,members in sorted(grouped.items()):
        pair = dict(pair_id=pair_id,task_ids=[t["task_id"] for t in members],same_raw_histories=None,
                    posterior_linf=None,logical_decisions_match=None,status="pending")
        if len(members) != 2 or len({(t["L"],t["p"],t["seed"],t["shots"],t["cutoff"]) for t in members}) != 1:
            raise ValueError("feasibility pair must vary only contraction dimension")
        left,right = [lookup[task["task_id"]] for task in members]
        if left["identity_valid"] and right["identity_valid"]:
            pair["same_raw_histories"] = left["record_ids"] == right["record_ids"]
            if not pair["same_raw_histories"]:
                pair["status"] = "mismatched_histories"
            elif left["runtime_model_identity"] != right["runtime_model_identity"]:
                pair["status"] = "mismatched_runtime_or_model"
            elif left["complete"] and right["complete"]:
                pair.update(status="compared",posterior_linf=max(
                    abs(a-b) for l,r in zip(left["outcomes"],right["outcomes"])
                    for a,b in zip(l["probabilities"],r["probabilities"])),
                    logical_decisions_match=all(l["logical_mask"] == r["logical_mask"]
                                                for l,r in zip(left["outcomes"],right["outcomes"])))
            else:
                pair["status"] = "matched_histories_incomplete_decodes"
        pairs.append(pair)
    return pairs


def collect_report(root):
    root = Path(root).resolve()
    manifest = HPC.load_manifest(root,verify_sources=True)
    validate_plan(manifest)
    by_stage = {stage:[inspect_task(root,task,manifest) for task in HPC.stage_tasks(manifest,stage)] for stage in STAGES}
    # Equal record IDs across independent seeds can be ordinary repeated physical
    # outcomes. Only run identities are unique; never deduplicate raw outcomes.
    run_ids = [row["run_id"] for rows in by_stage.values() for row in rows if row["identity_valid"]]
    if len(run_ids) != len(set(run_ids)):
        raise ValueError("duplicate saved run_id across planned tasks")
    identities = [row["runtime_model_identity"] for row in by_stage["production"] if row["identity_valid"]]
    if identities:
        mixed = sorted({field for identity in identities for field in RUNTIME_FIELDS+MODEL_FIELDS
                        if identity[field] != identities[0][field]})
        if mixed:
            raise ValueError("mixed production runtime or model identities: " + ", ".join(mixed))
    groups = defaultdict(list)
    for row in by_stage["production"]:
        groups[(row["L"],row["p"],row["chi"],row["cutoff"])].append(row)
    points = []
    for (L,p,chi,cutoff),rows in sorted(groups.items()):
        counts = totals(rows)
        N,F,U = counts["requested_shots"],counts["tn_failures"],counts["unknown_shots"]
        complete = all(row["complete"] for row in rows)
        points.append(dict(L=L,p=p,chi=chi,cutoff=cutoff,**counts,complete=complete,
                           planned_shards=len(rows),complete_shards=sum(row["complete"] for row in rows),
                           failure_rate=F/N if complete else None,
                           wilson95=list(BASE.wilson(F,N)) if complete else None,
                           outcome_bounds=[F/N,(F+U)/N],task_ids=[row["task_id"] for row in rows]))
    zeros = [dict(L=task["L"],p=0.0,analytical_rate=0.0,requested_shots=0,wilson95=None)
             for task in HPC.stage_tasks(manifest,"zero")]
    return dict(schema_version=1,generated_at_unix=time.time(),root=str(root),manifest_id=manifest["manifest_id"],
                manifest_sha256=hashlib.sha256((root/"manifest.json").read_bytes()).hexdigest(),
                scientific_source_sha256=manifest["scientific_source_sha256"],wrapper_sha256=manifest["wrapper_sha256"],
                production_runtime_model_identity=identities[0] if identities else None,
                production_points=points,production_totals=totals(by_stage["production"]),
                production_tasks=by_stage["production"],analytical_zero=zeros,
                diagnostics={stage:{**totals(by_stage[stage]),"tasks":by_stage[stage]} for stage in STAGES if stage != "production"},
                feasibility_pairs=feasibility_pairs(HPC.stage_tasks(manifest,"feasibility"),by_stage["feasibility"]),
                threshold_fit=None,
                counting_policy="Every manifest-planned production history remains in N, including absent whole shards. Only histories with bound identities, valid TN output/scoring, and required independently checked exact agreement are validated; all others remain unknown.",
                rate_policy="Failure rates and pointwise 95% Wilson intervals require every planned shard to pass its full audit. No rate from the valid subset; outcome bounds are not confidence intervals.",
                zero_policy="p=0 is the exact model prediction with no statistical shots or confidence interval. Numerical zero regressions remain separate diagnostics.",
                pooling_policy="Only this manifest's production tasks contribute to curves. No local-study data or diagnostic histories are pooled.")


def render(report, output):
    os.environ.setdefault("MPLCONFIGDIR",str(Path(tempfile.gettempdir()) / "tnml-hpc-matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size":11,"pdf.fonttype":42,"savefig.dpi":180})
    fig,ax = plt.subplots(figsize=(10,7))
    fig.subplots_adjust(left=.12,right=.96,top=.82,bottom=.25)
    fig.suptitle("No feedback after preparation",fontsize=21,y=.97)
    fig.text(.5,.916,"Z-only errors · diagonal H+6CNOT · initial seed repair · TN uses seed + final syndromes",ha="center",fontsize=10)
    points = report["production_points"]
    settings = sorted({(point["chi"],point["cutoff"]) for point in points})
    fig.text(.5,.873,"Finite-bond TN: " + "; ".join(f"χ={chi}, SVD cutoff={cutoff:g}" for chi,cutoff in settings),
             ha="center",fontsize=10)
    for L in sorted({point["L"] for point in points}):
        rows = sorted([point for point in points if point["L"] == L],key=lambda p:p["p"])
        complete = [row for row in rows if row["complete"]]
        if not complete:
            continue
        color = BASE.COLORS.get(L,"#333333")
        ax.plot([100*r["p"] for r in rows],[r["failure_rate"] if r["complete"] else math.nan for r in rows],
                color=color,lw=1.2,label=f"L={L}, N={rows[0]['requested_shots']:,}/point")
        for row in complete:
            BASE.draw_errorbar(ax,100*row["p"],row["failure_rate"],row["wilson95"],color=color)
    ax.plot(0,0,"k*",ms=13,clip_on=False,label="p=0: exact model, all requested sizes (no CI)")
    complete_count = sum(point["complete"] for point in points)
    if not complete_count:
        ax.text(.5,.5,"No fully audited production points yet",ha="center",transform=ax.transAxes,color=".4")
    ax.set(xlim=(0,100*max(p["p"] for p in points)*1.035),ylim=(0,1),
           xlabel="Gate-output error probability p (%)",ylabel="Either-logical failure probability")
    ax.spines[["top","right"]].set_visible(False)
    ax.legend(frameon=False,loc="upper left",fontsize=10)
    totals_ = report["production_totals"]
    fig.text(.12,.075,
             f"HPC production only: {complete_count}/{len(points)} complete points; "
             f"{totals_['validated_shots']:,}/{totals_['requested_shots']:,} validated histories.\n"
             "Bars: pointwise 95% Wilson intervals for complete targets. Incomplete targets have no plotted rate.\n"
             "Absent shards, invalid outputs, and unverified exact references remain unknown. Diagnostics are separate.\n"
             "No threshold fit or uncomputed size curves; analytical zero anchors are not Monte Carlo successes.",
             fontsize=9,linespacing=1.6)
    figures=[]
    for extension in ("png","pdf"):
        path=output/f"tnml-hpc.{extension}"
        fig.savefig(path,facecolor="white")
        figures.append(str(path.resolve()))
    plt.close(fig)
    return figures


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root",required=True,type=Path)
    parser.add_argument("--output",type=Path)
    parser.add_argument("--plot",action="store_true")
    args=parser.parse_args(argv)
    report=collect_report(args.root)
    output=(args.output or args.root/"collected").resolve()
    output.mkdir(parents=True,exist_ok=True)
    report["figures"]=render(report,output) if args.plot else []
    target=output/"report.json"
    with tempfile.NamedTemporaryFile("w",dir=output,delete=False) as stream:
        temporary=Path(stream.name)
        json.dump(BASE.clean_json(report),stream,indent=2,allow_nan=False)
        stream.write("\n")
    temporary.replace(target)
    print(json.dumps({**report["production_totals"],"report":str(target),"figures":report["figures"]}))
    return report


if __name__ == "__main__":
    main()
