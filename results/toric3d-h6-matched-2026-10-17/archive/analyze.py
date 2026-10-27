#!/usr/bin/env python3
"""Validate matched H+6 checkpoints and publish one CSV, JSON report, and PDF.

The frozen H6 inference engine remains unmodified. Incomplete grids are progress
reports only. Each record contributes three paired observations, never three
independent samples. Production input is accepted only from frozen shards.
"""
import argparse
import csv
import hashlib
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import re
import sys
import tempfile
import time
import tomllib

import numpy as np
from scipy.special import expit
from scipy.stats import norm

ARCHIVE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("_matched_h6_fitter", ARCHIVE / "analysis_engine/fit_toric3d_h6_threshold.py")
H6 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(H6)
ARMS = ("final_only", "one_step", "two_slice")
CHANNELS = ("x_only", "z_only")
LABELS = {"final_only": "No feedback after preparation", "one_step": "One-step", "two_slice": "Two-slice"}
# A given L has one color across sectors, pages, and future progress exports.
SIZE_COLORS = {3: "#7F7F7F", 5: "#0072B2", 7: "#E69F00", 9: "#009E73",
               11: "#CC79A7", 13: "#D55E00", 17: "#56B4E9", 21: "#7B3294", 25: "#7F8C2D"}
BASE_COLUMNS = ["channel", "steps", "Ly", "Lz", "p", "shots", "target_shots", "seed", "point_seed", "complete", "elapsed_seconds"]
COUNT_COLUMNS = [f"{arm}_{metric}" for arm in ARMS for metric in ("logical_1_failures", "logical_2_failures", "failures")]
JOINT_COLUMNS = [f"joint_{mask:02}" for mask in range(8)]
COLUMNS = BASE_COLUMNS + COUNT_COLUMNS + JOINT_COLUMNS + ["source_shard"]


class ActiveSnapshot(RuntimeError):
    """A worker is atomically replacing its checkpoint mirrors; retry later."""


def digest(value):
    return hashlib.sha256(value).hexdigest()


def integer(value, field):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{field} must be an integer, not Boolean")
    return H6.integer(value, field)


def boolean(value, field):
    if value is True or str(value).lower() == "true":
        return True
    if value is False or str(value).lower() == "false":
        return False
    raise ValueError(f"{field} must be Boolean")


def grid_key(row):
    return row["channel"], int(row["steps"]), float(row["p"])


def expected_grid(selection):
    return {(ch, L, p) for ch in CHANNELS for L in selection["channels"][ch]["sizes"]
            for p in selection["channels"][ch]["rates"]}


def validate_selection(selection):
    if selection.get("stage") != "production" or selection.get("frozen_before_sampling") is not True:
        raise ValueError("selection must be frozen before production sampling")
    if selection.get("arms") != list(ARMS) or set(selection.get("channels", {})) != set(CHANNELS):
        raise ValueError("selection arm order or channels differ")
    if integer(selection["shots"], "shots") <= 0 or integer(selection["batch_size"], "batch_size") <= 0:
        raise ValueError("positive target shots and batch size required")
    if selection.get("bootstrap_replicates") != 200:
        raise ValueError("the frozen analysis requires 200 bootstrap replicates")
    if selection.get("common_z_seed_repair") is not True or selection.get("pilot_excluded") is not True:
        raise ValueError("common Z seed preparation and pilot exclusion must be declared")
    for ch in CHANNELS:
        spec = selection["channels"][ch]
        sizes, rates = spec["sizes"], spec["rates"]
        if sizes != sorted(set(sizes)) or any(integer(L, "size") < 2 for L in sizes):
            raise ValueError("channel sizes must be increasing integers >=2")
        if rates != sorted(set(rates)) or any(not math.isfinite(p) or not 0 < p < .5 for p in rates):
            raise ValueError("channel rates must be increasing probabilities in (0,.5)")
        if set(spec["fit_rates"]) != set(ARMS):
            raise ValueError("every arm requires a frozen fit window")
        for arm in ARMS:
            fit = spec["fit_rates"][arm]
            if fit != sorted(set(fit)) or not set(fit) <= set(rates) or len(fit) < 7 or len(sizes) < 5:
                raise ValueError("fit windows need at least seven selected rates and five sizes for unchanged sensitivity guards")
    keys, seeds, point_seeds = set(), set(), set()
    for name, shard in selection["shards"].items():
        if Path(name).name != name or name in (".", ".."):
            raise ValueError("unsafe shard name")
        key = shard["channel"], shard["size"], shard["rate"]
        seed = integer(shard["seed"], "master seed")
        if key in keys or seed in seeds or seed < 0:
            raise ValueError("duplicate shard grid point or master seed")
        point_seed = shard.get("point_seed_hex")
        if not isinstance(point_seed, str) or not re.fullmatch(r"[0-9a-f]{16}", point_seed) or point_seed in point_seeds:
            raise ValueError("frozen point RNG seeds must be unique 16-digit hexadecimal strings")
        point_seeds.add(point_seed)
        keys.add(key); seeds.add(seed)
    if keys != expected_grid(selection):
        raise ValueError("shards do not exactly cover the frozen channel grids")
    return selection


def validate_row(original, selection):
    row = dict(original)
    for field in ("steps", "Ly", "Lz", "shots", "target_shots", "seed", *COUNT_COLUMNS, *JOINT_COLUMNS):
        row[field] = integer(row.get(field), field)
    row["p"] = float(row["p"])
    row["elapsed_seconds"] = float(row["elapsed_seconds"])
    row["complete"] = boolean(row.get("complete"), "complete")
    if row["channel"] not in CHANNELS or not row["steps"] == row["Ly"] == row["Lz"]:
        raise ValueError("matched geometry requires T=Ly=Lz=L")
    key = grid_key(row)
    if key not in expected_grid(selection):
        raise ValueError(f"point outside frozen grid: {key}")
    shard = [s for s in selection["shards"].values() if (s["channel"], s["size"], s["rate"]) == key]
    if len(shard) != 1 or row["seed"] != shard[0]["seed"]:
        raise ValueError("point master seed differs from frozen selection")
    if not re.fullmatch(r"[0-9a-f]{16}", str(row["point_seed"])):
        raise ValueError("point seed must be 16 lowercase hexadecimal digits")
    if "point_seed_hex" in shard[0] and row["point_seed"] != shard[0]["point_seed_hex"]:
        raise ValueError("point RNG seed differs from frozen selection")
    n = row["shots"]
    if row["target_shots"] != selection["shots"] or not 0 <= n <= selection["shots"]:
        raise ValueError("shot count or target differs from frozen selection")
    if row["complete"] != (n == selection["shots"]):
        raise ValueError("complete flag differs from shot count")
    if not math.isfinite(row["elapsed_seconds"]) or row["elapsed_seconds"] < 0:
        raise ValueError("elapsed time must be finite and nonnegative")
    bins = [row[k] for k in JOINT_COLUMNS]
    if any(k < 0 or k > n for k in bins) or sum(bins) != n:
        raise ValueError("joint histogram must sum to the shared shot count")
    for index, arm in enumerate(ARMS):
        a, b, k = [row[f"{arm}_{metric}"] for metric in ("logical_1_failures", "logical_2_failures", "failures")]
        if not 0 <= a <= k <= n or not 0 <= b <= k or k > a+b:
            raise ValueError(f"invalid logical union for {arm}")
        if sum(count for mask, count in enumerate(bins) if mask & (1 << index)) != k:
            raise ValueError(f"joint marginal differs from {arm} failures")
    return row


def validate_grid(rows, selection):
    result, seen = [], set()
    for row in rows:
        checked = validate_row(row, selection)
        key = grid_key(checked)
        if key in seen:
            raise ValueError(f"duplicate matched grid point: {key}")
        seen.add(key); result.append(checked)
    return sorted(result, key=grid_key)


def validate_payload(content):
    marker = re.findall(rb'^payload_sha256 = "([0-9a-f]{64})"\n', content, re.M)
    if len(marker) != 1:
        raise ValueError("checkpoint payload checksum missing or ambiguous")
    payload = re.sub(rb'^payload_sha256 = "[0-9a-f]{64}"\n', b'', content, count=1, flags=re.M)
    if digest(payload) != marker[0].decode():
        raise ValueError("checkpoint payload checksum mismatch")


def verify_frozen_sources(study, selection):
    archive = study / "archive"
    expected = {"runtime/"+key: value for key, value in selection["runtime_source_sha256"].items()}
    expected["run_study.jl"] = selection["runner_sha256"]
    analysis_hashes = selection["analysis_source_sha256"]
    if not isinstance(analysis_hashes, dict):
        raise ValueError("analysis_source_sha256 must map archive-relative files to hashes")
    for relative, claimed in {**expected, **analysis_hashes}.items():
        path = (archive / relative).resolve()
        if not path.is_relative_to(archive.resolve()) or not path.is_file() or digest(path.read_bytes()) != claimed:
            raise ValueError(f"frozen source hash mismatch: {relative}")
    current = {"runtime/"+str(p.relative_to(archive/"runtime")) for p in (archive/"runtime/src").rglob("*") if p.is_file()}
    current |= {"runtime/Project.toml", "runtime/Manifest.toml", "run_study.jl"}
    if set(expected) != current:
        raise ValueError("runtime source inventory differs from frozen selection")
    return expected


def stable_mirrors(folder, retries=8):
    paths = [folder / name for name in ("checkpoint.toml", "metadata.toml", "raw.csv")]
    error = None
    for attempt in range(retries):
        try:
            checkpoint_bytes = paths[0].read_bytes()
            validate_payload(checkpoint_bytes)
            metadata_bytes, raw_bytes = paths[1].read_bytes(), paths[2].read_bytes()
            metadata = tomllib.loads(metadata_bytes.decode())
            if checkpoint_bytes != paths[0].read_bytes():
                raise ActiveSnapshot("checkpoint changed during read")
            if metadata.get("checkpoint_sha256") != digest(checkpoint_bytes) or metadata.get("raw_sha256") != digest(raw_bytes):
                raise ActiveSnapshot("checkpoint mirrors are not yet synchronized")
            checkpoint = tomllib.loads(checkpoint_bytes.decode())
            return checkpoint, metadata, raw_bytes, dict(checkpoint_sha256=digest(checkpoint_bytes),
                metadata_sha256=digest(metadata_bytes), raw_sha256=digest(raw_bytes))
        except (OSError, ValueError, ActiveSnapshot) as exc:
            error = exc
            if attempt+1 < retries:
                time.sleep(.15)
    newest = max((path.stat().st_mtime for path in paths if path.exists()), default=0)
    if isinstance(error, (OSError, ActiveSnapshot)) and time.time()-newest < 10:
        raise ActiveSnapshot(f"{folder.name}: active mirror update; keep previous exports") from error
    raise ValueError(f"{folder.name}: invalid checkpoint/mirrors: {error}") from error


def validate_runtime(runtime):
    if runtime.get("julia_threads") != 1 or runtime.get("blas_threads") != 1:
        raise ValueError("runtime must use one Julia and BLAS thread")
    env = runtime["environment"]
    if any(env.get(k) != "1" for k in ("JULIA_NUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")):
        raise ValueError("runtime thread environment differs")
    if env.get("JULIA_CONDAPKG_BACKEND") != "Null" or Path(env["JULIA_PYTHONCALL_EXE"]).resolve() != Path(runtime["python_executable"]).resolve():
        raise ValueError("runtime Python signature differs")
    if not {"numpy", "scipy", "ldpc", "pymatching"} <= set(runtime["packages"]):
        raise ValueError("runtime dependency versions missing")


def load_inputs(study, selection, expected_sources):
    rows, provenance, common_runtime, point_seeds = [], {}, None, set()
    runs = study/"archive/runs"
    for folder in runs.iterdir() if runs.exists() else []:
        if folder.is_dir() and folder.name not in selection["shards"]:
            raise ValueError(f"unselected shard directory: {folder.name}")
    for name, shard in selection["shards"].items():
        folder = runs/name
        if not (folder/"checkpoint.toml").exists():
            if (folder/"raw.csv").exists() or (folder/"metadata.toml").exists():
                raise ValueError(f"{name}: mirrors exist without authoritative checkpoint")
            continue
        checkpoint, metadata, raw_bytes, hashes = stable_mirrors(folder)
        L, p, channel = shard["size"], shard["rate"], shard["channel"]
        config = dict(channel=channel, sizes=[L], rates=[p], shots=selection["shots"],
                      seed=shard["seed"], batch_size=selection["batch_size"], dimensions=[[L,L,L]])
        if checkpoint.get("config") != config:
            raise ValueError(f"{name}: checkpoint config differs from selection")
        for key in ("config", "runtime", "source_hashes", "protocol", "point_seeds_hex", "generation",
                    "started_at_unix", "updated_at_unix", "decoder_build_seconds_by_size", "decoder_build_count_by_size"):
            if metadata.get(key) != checkpoint.get(key):
                raise ValueError(f"{name}: metadata/checkpoint {key} mismatch")
        if checkpoint.get("source_hashes") != expected_sources or checkpoint.get("protocol") != selection["protocols"][channel]:
            raise ValueError(f"{name}: frozen source or protocol mismatch")
        if checkpoint.get("format_version") != 1 or metadata.get("format_version") != 1:
            raise ValueError("unsupported checkpoint version")
        runtime = checkpoint["runtime"]
        validate_runtime(runtime)
        if runtime != selection["expected_runtime"]:
            raise ValueError("shard runtime differs from the frozen runtime signature")
        if common_runtime is not None and runtime != common_runtime:
            raise ValueError("runtime differs between matched shards")
        common_runtime = runtime
        timestamps = (checkpoint["started_at_unix"], checkpoint["updated_at_unix"])
        if not all(math.isfinite(t) for t in timestamps) or timestamps[0] < selection["frozen_at_unix"] or timestamps[1] < timestamps[0]:
            raise ValueError("checkpoint predates frozen selection or has invalid timestamps")
        if integer(checkpoint["generation"], "checkpoint generation") < 0:
            raise ValueError("checkpoint generation cannot be negative")
        for field in ("decoder_build_seconds_by_size", "decoder_build_count_by_size"):
            if set(checkpoint[field]) != {str(L)}:
                raise ValueError("decoder build metadata size differs from selected point")
        if not math.isfinite(checkpoint["decoder_build_seconds_by_size"][str(L)]) or checkpoint["decoder_build_seconds_by_size"][str(L)] < 0:
            raise ValueError("invalid decoder build timing")
        if integer(checkpoint["decoder_build_count_by_size"][str(L)], "decoder build count") < 0:
            raise ValueError("invalid decoder build count")
        records = list(csv.DictReader(io.StringIO(raw_bytes.decode())))
        if len(records) != 1 or len(checkpoint["states"]) != 1 or len(checkpoint["point_seeds_hex"]) != 1:
            raise ValueError("one frozen point required in each shard")
        row = validate_row(records[0], selection)
        if grid_key(row) != (channel, L, p) or row["point_seed"] != checkpoint["point_seeds_hex"][0]:
            raise ValueError("CSV point differs from checkpoint")
        if row["point_seed"] in point_seeds:
            raise ValueError("duplicate point RNG stream across shards")
        point_seeds.add(row["point_seed"])
        state = checkpoint["states"][0]
        if state["completed_shots"] != row["shots"] or state["elapsed_seconds"] != row["elapsed_seconds"]:
            raise ValueError("CSV sample count/time differs from checkpoint")
        for state_key, suffix in (("logical_1_counts", "logical_1_failures"), ("logical_2_counts", "logical_2_failures"), ("failures", "failures")):
            if state[state_key] != [row[f"{arm}_{suffix}"] for arm in ARMS]:
                raise ValueError(f"CSV {state_key} differ from checkpoint")
        if state["joint_counts"] != [row[key] for key in JOINT_COLUMNS]:
            raise ValueError("CSV joint histogram differs from checkpoint")
        if metadata["expected_points"] != 1 or metadata["completed_points"] != int(row["complete"]) or metadata["complete"] != row["complete"] or metadata["completed_shots"] != row["shots"]:
            raise ValueError("metadata progress differs from checkpoint counts")
        row["source_shard"] = name
        rows.append(row)
        provenance[name] = dict(hashes, generation=checkpoint["generation"], point_seed=row["point_seed"],
                                master_seed=shard["seed"], completed_shots=row["shots"])
    return validate_grid(rows, selection), provenance, common_runtime


def wilson(k, n, z=None):
    if n <= 0 or not 0 <= k <= n:
        raise ValueError("Wilson interval requires 0<=k<=n and n>0")
    z = norm.ppf(.975) if z is None else z
    denominator = 1+z*z/n
    center = (k/n+z*z/(2*n))/denominator
    half = z*math.sqrt(k/n*(1-k/n)/n+z*z/(4*n*n))/denominator
    return (0. if k == 0 else max(0., center-half), 1. if k == n else min(1., center+half))


def paired_interval(positive, negative, shots):
    """Bonferroni-Wilson interval for p(comparison only)-p(reference only).

    Each discordance rate gets a two-sided 97.5% Wilson interval. Subtracting
    endpoints gives an approximate pointwise 95% paired interval; it is not a
    simultaneous confidence band over the grid and is not a threshold CI.
    """
    if shots <= 0 or min(positive, negative) < 0 or positive+negative > shots:
        raise ValueError("invalid paired discordance counts")
    z = norm.ppf(.9875)
    plus, minus = wilson(positive, shots, z), wilson(negative, shots, z)
    return dict(difference=(positive-negative)/shots, ci_low=max(-1.,plus[0]-minus[1]),
                ci_high=min(1.,plus[1]-minus[0]), comparison_only=positive, reference_only=negative,
                interval_method="pointwise approximate 95% Bonferroni-Wilson from paired discordant counts")


def paired_contrasts(row):
    if row["shots"] == 0:
        return []
    result = []
    for a, b in ((0,1),(0,2),(1,2)):
        positive = sum(row[f"joint_{m:02}"] for m in range(8) if m & (1<<b) and not m & (1<<a))
        negative = sum(row[f"joint_{m:02}"] for m in range(8) if m & (1<<a) and not m & (1<<b))
        result.append(dict(channel=row["channel"], L=row["steps"], p=row["p"], shots=row["shots"],
            complete=row["complete"], reference=ARMS[a], comparison=ARMS[b], **paired_interval(positive,negative,row["shots"])))
    return result


def fit_all(rows, selection, *, complete):
    if not complete:
        return {}
    fits = {}
    for arm_index, arm in enumerate(ARMS):
        adapted = []
        for row in rows:
            ch = row["channel"]
            if row["p"] not in selection["channels"][ch]["fit_rates"][arm]:
                continue
            adapted.append(dict(channel=ch, steps=row["steps"], Ly=row["Ly"], Lz=row["Lz"], p=row["p"],
                shots=row["shots"], failures=row[f"{arm}_failures"],
                logical_1_failures=row[f"{arm}_logical_1_failures"], logical_2_failures=row[f"{arm}_logical_2_failures"],
                seed=row["seed"], complete=row["complete"], elapsed_seconds=row["elapsed_seconds"], decoder=H6.DECODERS[ch]))
        # Both real sectors are supplied, each restricted to its own frozen arm
        # window. Do not synthesize a second channel to bypass engine validation.
        for ch_index, ch in enumerate(CHANNELS):
            fit = H6.analyze_channel(adapted, ch, bootstrap=selection["bootstrap_replicates"],
                                     seed=selection["bootstrap_seed"]+2*arm_index+ch_index)
            fit["arm"] = arm
            fit["fit_rates"] = selection["channels"][ch]["fit_rates"][arm]
            fits[f"{ch}:{arm}"] = fit
    return fits


def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k,v in value.items()}
    if isinstance(value, (list,tuple,np.ndarray)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def make_report(rows, selection, provenance, runtime, selection_hash):
    complete = len(rows) == len(expected_grid(selection)) and all(r["complete"] for r in rows)
    observed = {grid_key(r) for r in rows}
    completed = sum(r["shots"] for r in rows)
    report = dict(format_version=1, status="complete" if complete else "incomplete", complete=complete,
        generated_at_unix=time.time(), expected_points=len(expected_grid(selection)),
        available_points=len(rows), completed_points=sum(r["complete"] for r in rows),
        completed_fault_records=completed, target_fault_records=selection["shots"]*len(expected_grid(selection)),
        completed_protocol_outcomes=3*completed, target_shots_per_point=selection["shots"],
        missing_points=[dict(channel=ch,L=L,p=p) for ch,L,p in sorted(expected_grid(selection)-observed)],
        arm_order=list(ARMS), sizes_by_channel={ch:selection["channels"][ch]["sizes"] for ch in CHANNELS},
        fit_policy=H6.GUARD_POLICY, bootstrap_replicates=selection["bootstrap_replicates"],
        fits=fit_all(rows,selection,complete=complete),
        paired_contrasts=[contrast for row in rows for contrast in paired_contrasts(row)],
        provenance=dict(selection_sha256=selection_hash, source_sha256=selection["runtime_source_sha256"],
            runner_sha256=selection["runner_sha256"], analysis_source_sha256=selection["analysis_source_sha256"],
            runtime=runtime, shards=provenance),
        interpretation={"failure":"Either residual logical parity after the common final decoder.",
            "pairing":"Three outcomes replay the same original fault record; the eight-bin histogram preserves all pairwise contrasts.",
            "seed":"X uses an ideal encoded seed. Z applies the same full r=0 seed repair in every arm, retaining seed logical errors.",
            "two_slice_x":"With ideal idles, the logical X outcome is fixed after the first corrected two-slice window; this is preparation damage, not a noisy-memory threshold along T.",
            "intervals":"Marginals use pointwise 95% Wilson intervals. Paired differences use approximate pointwise 95% Bonferroni-Wilson intervals for discordant counts; no simultaneous grid coverage is claimed.",
            "threshold_differences":"The fitted thresholds are correlated. Do not infer threshold differences from independent marginal confidence intervals.",
            "incomplete":"No fit is attempted until every frozen point reaches the selected shot target."})
    return json_safe(report)


def plot_observed(axis, rows, channel, arm, sizes, *, local_rates=None):
    for L in sizes:
        points = sorted([r for r in rows if r["channel"] == channel and r["steps"] == L and r["shots"] > 0
                         and (local_rates is None or r["p"] in local_rates)], key=lambda r:r["p"])
        if not points:
            continue
        rate = np.array([r[f"{arm}_failures"]/r["shots"] for r in points])
        intervals = np.array([wilson(r[f"{arm}_failures"],r["shots"]) for r in points])
        color = SIZE_COLORS[L]
        axis.errorbar([100*r["p"] for r in points], rate,
            yerr=np.maximum(0,np.array([rate-intervals[:,0],intervals[:,1]-rate])),
            color=color, marker="o", markersize=3.7, linewidth=1.15, capsize=2,
            markerfacecolor=color, label=f"L={L}")
        partial = [r for r in points if not r["complete"]]
        if partial:
            axis.plot([100*r["p"] for r in partial], [r[f"{arm}_failures"]/r["shots"] for r in partial],
                      "o", color=color, markerfacecolor="white", markersize=4, linestyle="none")


def footer(figure, report, selection, page):
    count = report["completed_fault_records"]
    target = report["target_fault_records"]
    text = (f"T = L. Target: {selection['shots']:,} shared fault records per point; completed: {count:,}/{target:,}. "
            "Measurements, corrections and idles are ideal.\n"
            "X: ideal encoded seed. Z: identical complete seed repair before transfer 1 in all three arms; seed logical errors are retained.\n"
            "Two-slice X measures first-window preparation damage, not a noisy-memory threshold. Error bars: pointwise 95% Wilson intervals.\n"
            "Paired rate contrasts are retained in the report. Correlated threshold fits must not be compared using independent confidence intervals.")
    figure.text(.06,.025,text,fontsize=8.1,va="bottom",linespacing=1.45)
    figure.text(.965,.025,f"{page}/2",ha="right",fontsize=9)


def make_pdf(rows, selection, report, path, *, style_path=None):
    os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/toponoise-matched-mpl")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.lines import Line2D
    plt.style.use(style_path or ARCHIVE/"paper-no-grid.mplstyle")
    sizes = sorted({L for ch in CHANNELS for L in selection["channels"][ch]["sizes"]})
    if not set(sizes) <= set(SIZE_COLORS):
        raise ValueError("explicit size colors required for the selected sizes")
    handles = [Line2D([],[],color=SIZE_COLORS[L],marker="o",lw=1.3,ms=4,label=f"L={L}") for L in sizes]
    complete = report["complete"]
    status = "Completed matched experiment" if complete else "INCOMPLETE - progress only; no threshold fits"
    if report.get("synthetic_preview"):
        status = "SYNTHETIC LAYOUT PREVIEW - no production observations"
    with PdfPages(path) as pdf:
        fig, axes = plt.subplots(2,3,figsize=(16,10),sharex="row",sharey="row")
        fig.subplots_adjust(left=.065,right=.975,top=.84,bottom=.19,wspace=.19,hspace=.37)
        fig.suptitle("H+6 matched feedback comparison",fontsize=19,y=.97)
        fig.text(.5,.924,status,ha="center",fontsize=12,color="#137555" if complete else "#A23B32")
        fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.5,.91),ncol=len(sizes),frameon=False,fontsize=10)
        for ri,ch in enumerate(CHANNELS):
            channel_rows = [r for r in rows if r["channel"] == ch and r["shots"]]
            rates = selection["channels"][ch]["rates"]
            maximum = max((wilson(r[arm+"_failures"],r["shots"])[1] for r in channel_rows for arm in ARMS),default=.8)
            for ci,arm in enumerate(ARMS):
                axis=axes[ri,ci]
                plot_observed(axis,rows,ch,arm,selection["channels"][ch]["sizes"])
                axis.set_title(LABELS[arm],fontsize=12,pad=10)
                axis.set_xlim(0,100*max(rates)*1.035)
                axis.set_ylim(0,min(1.,max(.12,maximum*1.08)))
                axis.set_xlabel("Gate-output error probability p (%)",fontsize=10.5)
                axis.tick_params(labelsize=9)
                if ci == 0:
                    axis.set_ylabel(f"{ch[0].upper()} channel\nEither-logical failure probability",fontsize=11)
                if not channel_rows:
                    axis.text(.5,.5,"Awaiting sampled records",transform=axis.transAxes,ha="center",fontsize=11,color=".45")
        if not complete:
            fig.text(.5,.13,"Open markers: unfinished points. Missing points are not imputed. All three arms use equal counts at each displayed point.",
                     ha="center",fontsize=9,color="#A23B32")
        footer(fig,report,selection,1)
        pdf.savefig(fig); plt.close(fig)

        fig, axes = plt.subplots(2,3,figsize=(16,10))
        fig.subplots_adjust(left=.065,right=.975,top=.82,bottom=.20,wspace=.26,hspace=.46)
        fig.suptitle("Frozen local windows and scaling fits",fontsize=19,y=.97)
        fig.text(.5,.924,"Same H6 likelihood, 200 bootstrap replicates and unchanged acceptance guards"
                 if complete and not report.get("synthetic_preview") else status,
                 ha="center",fontsize=11,color=".25" if complete else "#A23B32")
        fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.5,.91),ncol=len(sizes),frameon=False,fontsize=10)
        for ri,ch in enumerate(CHANNELS):
            for ci,arm in enumerate(ARMS):
                axis=axes[ri,ci]
                rates=selection["channels"][ch]["fit_rates"][arm]
                selected=[r for r in rows if r["channel"] == ch and r["p"] in rates and r["shots"]]
                plot_observed(axis,rows,ch,arm,selection["channels"][ch]["sizes"],local_rates=rates)
                padding=(rates[-1]-rates[0])*.035
                axis.set_xlim(100*(rates[0]-padding),100*(rates[-1]+padding))
                axis.set_ylim(0,1)
                axis.set_xlabel("Gate-output error probability p (%)",fontsize=10.5)
                axis.tick_params(labelsize=9)
                if ci == 0:
                    axis.set_ylabel(f"{ch[0].upper()} channel\nEither-logical failure probability",fontsize=11)
                if not complete:
                    axis.set_title(LABELS[arm]+"\nFit withheld: incomplete experiment",fontsize=10.5,color="#A23B32")
                    continue
                fit=report["fits"][f"{ch}:{arm}"]
                accepted=fit["status"] == "accepted"
                axis.set_title(LABELS[arm]+("\nAccepted by all fit guards" if accepted else "\nDiagnostic only - threshold withheld"),
                               fontsize=10.5,color="#137555" if accepted else "#A23B32")
                candidate=fit.get("candidate")
                if candidate is None:
                    axis.text(.5,.5,"No numerical fit",ha="center",transform=axis.transAxes,color=".45")
                    continue
                inset=axis.inset_axes([.12,.57,.43,.34])
                inset.set_title("Accepted collapse" if accepted else "Diagnostic collapse",fontsize=7.5,pad=3)
                all_scaled=[]
                for L in selection["channels"][ch]["sizes"]:
                    points=[r for r in selected if r["steps"] == L]
                    x=np.array([(r["p"]-candidate["pc"])/candidate["p_scale"]*(L/candidate["L_reference"])**(1/candidate["nu"]) for r in points])
                    all_scaled.extend(x)
                    inset.plot(x,[r[arm+"_failures"]/r["shots"] for r in points],"o",color=SIZE_COLORS[L],ms=2.4)
                    p_grid=np.linspace(rates[0],rates[-1],180)
                    scaled=(p_grid-candidate["pc"])/candidate["p_scale"]*(L/candidate["L_reference"])**(1/candidate["nu"])
                    axis.plot(100*p_grid,expit(np.polynomial.polynomial.polyval(scaled,candidate["coefficients"])),
                              color=SIZE_COLORS[L],lw=.9,alpha=.8)
                if all_scaled:
                    grid=np.linspace(min(all_scaled),max(all_scaled),180)
                    inset.plot(grid,expit(np.polynomial.polynomial.polyval(grid,candidate["coefficients"])),color=".15",lw=1)
                inset.tick_params(labelsize=6.5,length=2)
                inset.set_xlabel("Scaled p",fontsize=7,labelpad=1)
                inset.set_ylabel("Failure",fontsize=7,labelpad=1)
                if accepted:
                    interval=fit["bootstrap"]["pc_interval"]
                    label=(f"p_c = {100*fit['pc']:.4f}%   nu = {fit['nu']:.3f}\n"
                           f"95% conditional CI: {100*interval[0]:.4f}-{100*interval[1]:.4f}%")
                else:
                    label=f"{len(fit['failed_guards'])} acceptance guard(s) failed; see report."
                axis.text(.98,.035,label,transform=axis.transAxes,ha="right",va="bottom",fontsize=7.7,
                          bbox=dict(facecolor="white",edgecolor="none",alpha=.9,pad=2))
        fig.text(.5,.135,"Fit windows were frozen before sampling. A numerical collapse alone does not establish a threshold.",ha="center",fontsize=9)
        footer(fig,report,selection,2)
        pdf.savefig(fig); plt.close(fig)


def export(study, rows, selection, report):
    data, figures = study/"data", study/"figures"
    data.mkdir(parents=True,exist_ok=True); figures.mkdir(parents=True,exist_ok=True)
    # Stage every output before replacing any prior export. A validation or
    # plotting failure leaves the complete previous export set untouched.
    with tempfile.TemporaryDirectory(prefix=".matched-export-",dir=data) as temporary:
        temporary=Path(temporary)
        csv_path=temporary/"matched-results.csv"
        with csv_path.open("w",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=COLUMNS,extrasaction="ignore")
            writer.writeheader(); writer.writerows(rows)
        pdf_path=temporary/"h6-matched-comparison.pdf"
        make_pdf(rows,selection,report,pdf_path,style_path=study/"archive/paper-no-grid.mplstyle")
        report["output_sha256"]={"data/matched-results.csv":digest(csv_path.read_bytes()),
                                  "figures/h6-matched-comparison.pdf":digest(pdf_path.read_bytes())}
        report_path=temporary/"report.json"
        report_path.write_text(json.dumps(json_safe(report),indent=2,allow_nan=False)+"\n")
        os.replace(csv_path,data/csv_path.name)
        os.replace(pdf_path,figures/pdf_path.name)
        os.replace(report_path,data/report_path.name)


def main(arguments=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study",required=True,type=Path)
    parser.add_argument("--progress",action="store_true",help="permit an incomplete progress export; never fit incomplete data")
    args=parser.parse_args(arguments)
    try:
        study=args.study.resolve()
        selection_bytes=(study/"data/selection.json").read_bytes()
        selection=validate_selection(json.loads(selection_bytes))
        sources=verify_frozen_sources(study,selection)
        rows,provenance,runtime=load_inputs(study,selection,sources)
        complete=len(rows) == len(expected_grid(selection)) and all(r["complete"] for r in rows)
        if not complete and not args.progress:
            raise ValueError("experiment incomplete; use --progress for an explicitly incomplete report")
        report=make_report(rows,selection,provenance,runtime,digest(selection_bytes))
        export(study,rows,selection,report)
        print(json.dumps({key:report[key] for key in ("status","completed_points","expected_points","completed_fault_records","target_fault_records")},sort_keys=True))
        return 0
    except ActiveSnapshot as error:
        print(str(error),file=sys.stderr)
        return 0 if args.progress else 2
    except (ValueError,KeyError,OSError) as error:
        print(f"analysis rejected input: {error}",file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
