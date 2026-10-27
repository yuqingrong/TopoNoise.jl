#!/usr/bin/env python3
"""Validate the physical-x scout and compare it with frozen diagonal records.

Only the three arms within a direction are paired. Directions use independent
records. This report deliberately performs no finite-size threshold fitting.
"""
import argparse
import csv
import hashlib
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

from scipy.stats import norm

ARCHIVE = Path(__file__).resolve().parent
ARMS = ("final_only", "one_step", "two_slice")
CHANNELS = ("x_only", "z_only")
LABELS = {"final_only": "No feedback after preparation", "one_step": "One-step feedback", "two_slice": "Two-slice feedback"}
SIZE_COLORS = {3: "#7F7F7F", 5: "#0072B2", 7: "#E69F00", 9: "#009E73", 11: "#CC79A7"}
BASE_COLUMNS = ["channel", "steps", "Ly", "Lz", "p", "shots", "target_shots", "seed", "point_seed", "complete", "elapsed_seconds"]
COUNT_COLUMNS = [f"{arm}_{metric}" for arm in ARMS for metric in ("logical_1_failures", "logical_2_failures", "failures")]
JOINT_COLUMNS = [f"joint_{mask:02}" for mask in range(8)]
COLUMNS = BASE_COLUMNS + COUNT_COLUMNS + JOINT_COLUMNS + ["source_shard"]
THREAD_KEYS = ("JULIA_NUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
GEOMETRY_CAVEATS = [
    "Equal T=Ly=Lz=L matches L^3 growth cells and 13 L^3 noisy transfer gate outputs, not code distance or seed geometry.",
    "Physical-x growth has a square yz seed with 2 L^2 seed qubits and d_Z=L; diagonal growth has 3 L^2 seed qubits and d_Z=2L. Both have d_X=(L+1)L.",
    "X seeds are ideal. Z has independent seed-H faults at 2 L^2 physical-x versus 3 L^2 diagonal seed outputs, then common seed syndrome repair within each direction; logical seed errors are retained.",
    "Physical-x Z one-step feedback deterministically pushes retiring star defects one +x edge. Diagonal Z one-step uses MWPM endpoint routing; these are distinct orientation-adapted policies.",
    "Two-slice feedback includes the final transfer, whereas Z one-step intermediate feedback ends at T-1. All arms use one shared final decoder within their direction and sector.",
    "Two-slice X failure is a preparation outcome fixed after the first corrected window, not a noisy-memory threshold in transfer count.",
    "This is an orientation scout. No axis threshold fit is attempted; historical diagonal fit windows are not presumed valid and a visual crossing does not establish a threshold.",
]


class ActiveSnapshot(RuntimeError):
    """Checkpoint mirrors are currently being replaced; a progress read may skip."""


def digest(content):
    return hashlib.sha256(content).hexdigest()


def integer(value, field):
    if isinstance(value, bool):
        raise ValueError(f"{field} must be an integer, not Boolean")
    if isinstance(value, int):
        return value
    if isinstance(value, str) and re.fullmatch(r"[+-]?\d+", value):
        return int(value)
    raise ValueError(f"{field} must be an integer")


def finite(value, field):
    if isinstance(value, bool):
        raise ValueError(f"{field} must be finite numeric data")
    try:
        value = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{field} must be finite numeric data") from error
    if not math.isfinite(value):
        raise ValueError(f"{field} must be finite numeric data")
    return value


def boolean(value, field):
    if value is True or isinstance(value, str) and value.lower() == "true":
        return True
    if value is False or isinstance(value, str) and value.lower() == "false":
        return False
    raise ValueError(f"{field} must be Boolean")


def grid_key(row):
    return row["channel"], int(row["steps"]), float(row["p"])


def expected_grid(selection):
    return {(ch, L, p) for ch in CHANNELS for L in selection["channels"][ch]["sizes"]
            for p in selection["channels"][ch]["rates"]}


def validate_selection(selection):
    if selection.get("stage") != "direction_comparison" or selection.get("frozen_before_sampling") is not True:
        raise ValueError("direction comparison selection must be frozen before sampling")
    if selection.get("arms") != list(ARMS) or set(selection.get("channels", {})) != set(CHANNELS):
        raise ValueError("selection arms or channels differ")
    if integer(selection["shots"], "shots") != 1000 or integer(selection["batch_size"], "batch_size") <= 0:
        raise ValueError("1000 target shots and positive batch size required")
    if finite(selection["frozen_at_unix"], "freeze timestamp") < 0:
        raise ValueError("invalid freeze timestamp")
    if selection.get("common_z_seed_repair") is not True or selection.get("pilot_excluded") is not True:
        raise ValueError("common Z preparation and excluded pilots must be declared")
    for ch in CHANNELS:
        spec = selection["channels"][ch]
        sizes, rates = spec["sizes"], spec["rates"]
        if not sizes or sizes != sorted(set(sizes)) or any(integer(L, "size") < 2 for L in sizes):
            raise ValueError("channel sizes must be strictly increasing integers >=2")
        if not rates or rates != sorted(set(rates)) or any(not 0 <= finite(p, "rate") < .5 for p in rates):
            raise ValueError("channel rates must be strictly increasing probabilities in [0,.5)")
    points, master_seeds, point_seeds = set(), set(), set()
    for name, shard in selection["shards"].items():
        if Path(name).name != name or name in (".", ".."):
            raise ValueError("unsafe shard name")
        ch, L, rates = shard["channel"], integer(shard["size"], "size"), shard["rates"]
        if ch not in CHANNELS or L not in selection["channels"][ch]["sizes"] or rates != selection["channels"][ch]["rates"]:
            raise ValueError("each shard must contain one selected size and its entire channel rate grid")
        seed = integer(shard["seed"], "master seed")
        if seed < 0 or seed in master_seeds:
            raise ValueError("master seeds must be nonnegative and unique")
        master_seeds.add(seed)
        seeds = shard["point_seeds_hex"]
        if len(seeds) != len(rates):
            raise ValueError("point seed count differs from rate count")
        for p, point_seed in zip(rates, seeds):
            key = ch, L, p
            if key in points:
                raise ValueError("duplicate shard grid point")
            if not isinstance(point_seed, str) or not re.fullmatch(r"[0-9a-f]{16}", point_seed) or point_seed in point_seeds:
                raise ValueError("point seeds must be unique 16-digit lowercase hexadecimal strings")
            points.add(key); point_seeds.add(point_seed)
    if points != expected_grid(selection):
        raise ValueError("shards do not exactly cover the frozen grid")
    target = len(points)*selection["shots"]
    if integer(selection["total_fault_records"], "total records") != target or integer(selection["total_protocol_outcomes"], "total outcomes") != 3*target:
        raise ValueError("selection sample totals disagree with the grid")
    return selection


def validate_counts(original):
    row = dict(original)
    for field in ("steps", "Ly", "Lz", "shots", "target_shots", "seed", *COUNT_COLUMNS, *JOINT_COLUMNS):
        row[field] = integer(row.get(field), field)
    row["p"] = finite(row["p"], "p")
    row["elapsed_seconds"] = finite(row["elapsed_seconds"], "elapsed_seconds")
    row["complete"] = boolean(row.get("complete"), "complete")
    if row["channel"] not in CHANNELS or not row["steps"] == row["Ly"] == row["Lz"] or row["steps"] < 2:
        raise ValueError("geometry requires T=Ly=Lz=L>=2")
    if not 0 <= row["p"] < .5 or row["elapsed_seconds"] < 0 or row["seed"] < 0:
        raise ValueError("invalid rate, seed, or elapsed time")
    if not re.fullmatch(r"[0-9a-f]{16}", str(row["point_seed"])):
        raise ValueError("invalid point seed encoding")
    n = row["shots"]
    if not 0 <= n <= row["target_shots"] or row["target_shots"] <= 0 or row["complete"] != (n == row["target_shots"]):
        raise ValueError("invalid shot count or completion flag")
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


def validate_row(original, selection):
    row = validate_counts(original)
    key = grid_key(row)
    if key not in expected_grid(selection):
        raise ValueError("point outside frozen grid")
    matches = [s for s in selection["shards"].values() if s["channel"] == row["channel"] and s["size"] == row["steps"] and row["p"] in s["rates"]]
    if len(matches) != 1 or row["seed"] != matches[0]["seed"]:
        raise ValueError("point master seed differs from frozen selection")
    shard = matches[0]
    if row["point_seed"] != shard["point_seeds_hex"][shard["rates"].index(row["p"])]:
        raise ValueError("point RNG seed differs from frozen selection")
    if row["target_shots"] != selection["shots"]:
        raise ValueError("target shots differ from selection")
    return row


def verify_frozen_sources(study, selection):
    archive = (study/"archive").resolve()
    sources, auxiliary = selection["source_hashes"], selection["auxiliary_source_hashes"]
    required = {"runtime/Project.toml", "runtime/Manifest.toml", "run_study.jl", "AxisH6.jl", "axis_geometry.jl", "axis_feedback.jl"}
    required |= {str(p.relative_to(archive)) for p in (archive/"runtime/src").rglob("*") if p.is_file()}
    if set(sources) != required:
        raise ValueError("frozen runtime/axis source inventory differs")
    if not {"analyze.py", "supervise.py", "paper-no-grid.mplstyle"} <= set(auxiliary) or set(sources) & set(auxiliary):
        raise ValueError("auxiliary source inventory differs")
    for relative, claimed in {**sources, **auxiliary}.items():
        path = (archive/relative).resolve()
        if not path.is_relative_to(archive) or not path.is_file() or digest(path.read_bytes()) != claimed:
            raise ValueError(f"frozen source hash mismatch: {relative}")
    return sources


def validate_payload(content):
    marker = re.findall(rb'^payload_sha256 = "([0-9a-f]{64})"\n', content, re.M)
    if len(marker) != 1:
        raise ValueError("checkpoint payload checksum missing or ambiguous")
    payload = re.sub(rb'^payload_sha256 = "[0-9a-f]{64}"\n', b'', content, count=1, flags=re.M)
    if digest(payload) != marker[0].decode():
        raise ValueError("checkpoint payload checksum mismatch")


def stable_mirrors(folder, retries=8):
    paths = [folder/name for name in ("checkpoint.toml", "metadata.toml", "raw.csv")]
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
            return tomllib.loads(checkpoint_bytes.decode()), metadata, raw_bytes, dict(
                checkpoint_sha256=digest(checkpoint_bytes), metadata_sha256=digest(metadata_bytes), raw_sha256=digest(raw_bytes))
        except (OSError, ActiveSnapshot) as exc:
            error = exc
            if attempt+1 < retries:
                time.sleep(.15)
    newest = max((path.stat().st_mtime for path in paths if path.exists()), default=0)
    if time.time()-newest < 10:
        raise ActiveSnapshot(f"{folder.name}: active mirror update") from error
    raise ValueError(f"{folder.name}: invalid checkpoint/mirrors: {error}") from error


def validate_runtime(runtime):
    if integer(runtime["julia_threads"], "Julia threads") != 1 or integer(runtime["blas_threads"], "BLAS threads") != 1:
        raise ValueError("single-threaded runtime required")
    env = runtime["environment"]
    if any(env.get(k) != "1" for k in THREAD_KEYS) or env.get("JULIA_CONDAPKG_BACKEND") != "Null":
        raise ValueError("runtime environment differs")
    if not env.get("JULIA_PYTHONCALL_EXE") or Path(env["JULIA_PYTHONCALL_EXE"]).resolve() != Path(runtime["python_executable"]).resolve():
        raise ValueError("runtime Python signature differs")
    if not {"numpy", "scipy", "ldpc", "pymatching"} <= set(runtime["packages"]):
        raise ValueError("runtime dependency versions missing")


def load_inputs(study, selection, expected_sources, *, progress=False):
    rows, provenance, common_runtime, skipped = [], {}, None, []
    runs = study/"archive/runs"
    if runs.exists():
        for folder in runs.iterdir():
            if folder.is_dir() and folder.name not in selection["shards"]:
                raise ValueError(f"unselected shard directory: {folder.name}")
    for name, shard in selection["shards"].items():
        folder = runs/name
        if not (folder/"checkpoint.toml").exists():
            if (folder/"raw.csv").exists() or (folder/"metadata.toml").exists():
                raise ValueError(f"{name}: mirrors exist without authoritative checkpoint")
            continue
        try:
            checkpoint, metadata, raw_bytes, hashes = stable_mirrors(folder)
        except ActiveSnapshot:
            if not progress:
                raise
            skipped.append(name)
            continue
        L, rates, ch = shard["size"], shard["rates"], shard["channel"]
        expected_config = dict(channel=ch, sizes=[L], rates=rates, shots=selection["shots"],
            seed=shard["seed"], batch_size=selection["batch_size"], dimensions=[[L,L,L]])
        if checkpoint["config"] != expected_config:
            raise ValueError(f"{name}: checkpoint config differs from selection")
        for key in ("config", "runtime", "source_hashes", "protocol", "point_seeds_hex", "generation", "started_at_unix", "updated_at_unix", "decoder_build_seconds_by_size", "decoder_build_count_by_size"):
            if metadata.get(key) != checkpoint.get(key):
                raise ValueError(f"{name}: checkpoint/metadata {key} mismatch")
        if checkpoint["source_hashes"] != expected_sources or checkpoint["protocol"] != selection["protocols"][ch]:
            raise ValueError(f"{name}: source or protocol differs from frozen selection")
        if integer(checkpoint["format_version"], "format version") != 1 or integer(metadata["format_version"], "format version") != 1:
            raise ValueError("unsupported checkpoint version")
        if checkpoint["point_seeds_hex"] != shard["point_seeds_hex"]:
            raise ValueError(f"{name}: checkpoint point seeds differ from frozen selection")
        runtime = checkpoint["runtime"]; validate_runtime(runtime)
        if "expected_runtime" in selection and runtime != selection["expected_runtime"]:
            raise ValueError("runtime differs from frozen signature")
        if common_runtime is not None and common_runtime != runtime:
            raise ValueError("runtime differs between shards")
        common_runtime = runtime
        started, updated = [finite(checkpoint[k], k) for k in ("started_at_unix", "updated_at_unix")]
        if started < selection["frozen_at_unix"] or updated < started or integer(checkpoint["generation"], "generation") < 1:
            raise ValueError("invalid generation or checkpoint timestamps")
        for field in ("decoder_build_seconds_by_size", "decoder_build_count_by_size"):
            if set(checkpoint[field]) != {str(L)}:
                raise ValueError("decoder build size metadata differs")
        if finite(checkpoint["decoder_build_seconds_by_size"][str(L)], "build time") < 0 or integer(checkpoint["decoder_build_count_by_size"][str(L)], "build count") < 0:
            raise ValueError("invalid decoder build statistics")
        records = list(csv.DictReader(io.StringIO(raw_bytes.decode())))
        if len(records) != len(rates) or len(checkpoint["states"]) != len(rates):
            raise ValueError(f"{name}: checkpoint/CSV rate-grid length differs")
        checked, unfinished = [], False
        for p, state, original in zip(rates, checkpoint["states"], records):
            row = validate_row(original, selection)
            if grid_key(row) != (ch,L,p):
                raise ValueError("CSV rate order differs from checkpoint grid")
            n = integer(state["completed_shots"], "completed shots")
            if n != row["shots"] or finite(state["elapsed_seconds"], "elapsed time") != row["elapsed_seconds"]:
                raise ValueError("CSV sample count/time differs from checkpoint")
            if unfinished and n != 0:
                raise ValueError("checkpoint sampling order is not sequential")
            unfinished |= n < selection["shots"]
            for state_key, suffix in (("logical_1_counts", "logical_1_failures"), ("logical_2_counts", "logical_2_failures"), ("failures", "failures")):
                values = [integer(v, state_key) for v in state[state_key]]
                if values != [row[f"{arm}_{suffix}"] for arm in ARMS]:
                    raise ValueError(f"CSV {state_key} differs from checkpoint")
            if [integer(v,"joint count") for v in state["joint_counts"]] != [row[k] for k in JOINT_COLUMNS]:
                raise ValueError("CSV joint histogram differs from checkpoint")
            row["source_shard"] = name; checked.append(row)
        if integer(metadata["expected_points"], "expected points") != len(rates) or integer(metadata["completed_points"], "completed points") != sum(r["complete"] for r in checked) or integer(metadata["completed_shots"], "completed shots") != sum(r["shots"] for r in checked) or boolean(metadata["complete"], "complete") != all(r["complete"] for r in checked):
            raise ValueError("metadata progress differs from checkpoint")
        rows.extend(checked)
        provenance[name] = dict(hashes, generation=checkpoint["generation"], completed_shots=sum(r["shots"] for r in checked), point_seeds_hex=checkpoint["point_seeds_hex"])
    keys = [grid_key(r) for r in rows]
    if len(keys) != len(set(keys)):
        raise ValueError("duplicate axis point")
    return sorted(rows,key=grid_key), provenance, common_runtime, skipped


def load_reference(selection):
    reference = selection["reference"]
    csv_bytes, report_bytes = Path(reference["csv"]).read_bytes(), Path(reference["report"]).read_bytes()
    if digest(csv_bytes) != reference["csv_sha256"] or digest(report_bytes) != reference["report_sha256"]:
        raise ValueError("diagonal reference hash mismatch")
    report = json.loads(report_bytes)
    if report.get("complete") is not True or report.get("output_sha256",{}).get("data/matched-results.csv") != digest(csv_bytes):
        raise ValueError("diagonal reference report is incomplete or does not authenticate its CSV")
    rows = [validate_counts(r) for r in csv.DictReader(io.StringIO(csv_bytes.decode()))]
    keys = [grid_key(r) for r in rows]
    if len(keys) != len(set(keys)) or set(keys) != expected_grid(selection) or any(not r["complete"] or r["shots"] != selection["shots"] for r in rows):
        raise ValueError("diagonal reference must be complete on the identical comparison grid")
    axis_seeds = {s for shard in selection["shards"].values() for s in shard["point_seeds_hex"]}
    axis_masters = {shard["seed"] for shard in selection["shards"].values()}
    reference_seeds = [r["point_seed"] for r in rows]
    if len(reference_seeds) != len(set(reference_seeds)) or axis_seeds.intersection(reference_seeds) or axis_masters.intersection(r["seed"] for r in rows):
        raise ValueError("direction comparison requires independent RNG streams")
    return sorted(rows,key=grid_key)


def require_complete(rows, selection, *, progress):
    complete = {grid_key(r) for r in rows} == expected_grid(selection) and all(r["complete"] for r in rows)
    if not complete and not progress:
        raise ValueError("experiment incomplete; use --progress for an explicitly incomplete report")
    return complete


def wilson(k,n,z=None):
    if n <= 0 or not 0 <= k <= n:
        raise ValueError("Wilson interval requires 0<=k<=n and n>0")
    z = float(norm.ppf(.975)) if z is None else z
    denominator = 1+z*z/n
    center = (k/n+z*z/(2*n))/denominator
    half = z*math.sqrt(k/n*(1-k/n)/n+z*z/(4*n*n))/denominator
    return (0. if k==0 else max(0.,center-half), 1. if k==n else min(1.,center+half))


def paired_contrasts(row):
    if row["shots"] == 0:
        return []
    result = []
    for a,b in ((0,1),(0,2),(1,2)):
        positive = sum(row[f"joint_{m:02}"] for m in range(8) if m & (1<<b) and not m & (1<<a))
        negative = sum(row[f"joint_{m:02}"] for m in range(8) if m & (1<<a) and not m & (1<<b))
        plus, minus = wilson(positive,row["shots"],float(norm.ppf(.9875))), wilson(negative,row["shots"],float(norm.ppf(.9875)))
        result.append(dict(channel=row["channel"],L=row["steps"],p=row["p"],shots=row["shots"],
            reference=ARMS[a],comparison=ARMS[b],comparison_only=positive,reference_only=negative,
            difference=(positive-negative)/row["shots"],ci_low=max(-1.,plus[0]-minus[1]),ci_high=min(1.,plus[1]-minus[0]),
            sampling="paired within physical-x direction",interval_method="pointwise approximate 95% Bonferroni-Wilson from paired discordant counts"))
    return result


def independent_difference(axis_failures,axis_shots,diagonal_failures,diagonal_shots):
    if axis_shots <= 0 or diagonal_shots <= 0 or not 0 <= axis_failures <= axis_shots or not 0 <= diagonal_failures <= diagonal_shots:
        raise ValueError("invalid independent binomial counts")
    a,d = axis_failures/axis_shots,diagonal_failures/diagonal_shots
    difference = a-d
    standard_error = math.sqrt(a*(1-a)/axis_shots+d*(1-d)/diagonal_shots)
    axis_low,axis_high = wilson(axis_failures,axis_shots)
    diagonal_low,diagonal_high = wilson(diagonal_failures,diagonal_shots)
    low = difference-math.hypot(a-axis_low,diagonal_high-d)
    high = difference+math.hypot(axis_high-a,d-diagonal_low)
    return dict(axis_rate=a,diagonal_rate=d,axis_shots=axis_shots,diagonal_shots=diagonal_shots,
        difference=difference,standard_error=standard_error,ci_low=max(-1.,low),ci_high=min(1.,high),
        sampling="independent across directions",interval_method="pointwise 95% Newcombe unpooled Wilson-score difference interval; standard_error is the independent plug-in binomial SE")


def make_report(rows, reference_rows, selection, provenance, runtime, selection_sha256, *, skipped=None):
    complete = require_complete(rows,selection,progress=True)
    reference = {grid_key(r):r for r in reference_rows}
    differences = []
    for row in rows:
        if not row["shots"]:
            continue
        diagonal = reference[grid_key(row)]
        for arm in ARMS:
            differences.append(dict(channel=row["channel"],L=row["steps"],p=row["p"],arm=arm,
                axis_complete=row["complete"],**independent_difference(row[arm+"_failures"],row["shots"],diagonal[arm+"_failures"],diagonal["shots"])))
    return dict(format_version=1,status="complete" if complete else "incomplete_progress",complete=complete,
        generated_at_unix=time.time(),expected_points=len(expected_grid(selection)),available_points=len(rows),
        completed_points=sum(r["complete"] for r in rows),completed_fault_records=sum(r["shots"] for r in rows),
        target_fault_records=selection["total_fault_records"],completed_protocol_outcomes=3*sum(r["shots"] for r in rows),
        target_protocol_outcomes=selection["total_protocol_outcomes"],target_shots_per_point=selection["shots"],
        missing_points=[dict(channel=ch,L=L,p=p) for ch,L,p in sorted(expected_grid(selection)-{grid_key(r) for r in rows})],
        skipped_updating_shards=skipped or [],arm_order=list(ARMS),direction_order=["diagonal","physical_x"],
        sizes_by_channel={ch:selection["channels"][ch]["sizes"] for ch in CHANNELS},
        threshold_status="No physical-x threshold fit attempted or validated; orientation scout only.",fits={},
        paired_within_axis_contrasts=[item for row in rows for item in paired_contrasts(row)],
        independent_axis_minus_diagonal=differences,geometry_caveats=GEOMETRY_CAVEATS,
        provenance=dict(selection_sha256=selection_sha256,source_hashes=selection["source_hashes"],
            auxiliary_source_hashes=selection["auxiliary_source_hashes"],reference=selection["reference"],runtime=runtime,shards=provenance),
        interpretation="Across-direction differences use independent-sample Newcombe Wilson-score intervals and unpooled binomial SE. Within-axis arm differences use shared-record discordances. Non-detection does not prove zero failure probability or equal methods. Intervals are pointwise, not simultaneous; geometry and correction-policy differences preclude a pure direction effect at fixed distance.")


def plot_observed(axis,rows,ch,arm,sizes):
    for L in sizes:
        points = sorted((r for r in rows if r["channel"]==ch and r["steps"]==L and r["shots"]>0),key=lambda r:r["p"])
        if not points:
            continue
        y = [r[arm+"_failures"]/r["shots"] for r in points]
        intervals = [wilson(r[arm+"_failures"],r["shots"]) for r in points]
        axis.errorbar([100*r["p"] for r in points],y,
            yerr=[[max(0.,v-lo) for v,(lo,hi) in zip(y,intervals)],[max(0.,hi-v) for v,(lo,hi) in zip(y,intervals)]],
            fmt="o-",color=SIZE_COLORS[L],markersize=3.3,linewidth=1.1,capsize=1.8)
        partial = [r for r in points if not r["complete"]]
        axis.plot([100*r["p"] for r in partial],[r[arm+"_failures"]/r["shots"] for r in partial],
            "o",color=SIZE_COLORS[L],markerfacecolor="white",markersize=4,linestyle="none")


def make_pdf(rows,reference_rows,selection,report,path,*,style_path=None):
    os.environ.setdefault("MPLCONFIGDIR","/private/tmp/toponoise-axis-mpl")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.lines import Line2D
    plt.style.use(style_path or ARCHIVE/"paper-no-grid.mplstyle")
    sizes = sorted({L for ch in CHANNELS for L in selection["channels"][ch]["sizes"]})
    if not set(sizes) <= set(SIZE_COLORS):
        raise ValueError("explicit consistent colors required for every size")
    handles = [Line2D([],[],color=SIZE_COLORS[L],marker="o",lw=1.2,ms=4,label=f"L={L}") for L in sizes]
    with PdfPages(path) as pdf:
        for page,ch in enumerate(CHANNELS,1):
            fig,axes = plt.subplots(3,2,figsize=(14,13),sharex="row",sharey="row")
            fig.subplots_adjust(left=.09,right=.975,top=.86,bottom=.19,wspace=.14,hspace=.43)
            fig.suptitle(f"H+6 growth direction comparison | {ch[0].upper()}-only faults",fontsize=19,y=.974)
            status = "Complete orientation scout; no axis threshold fit validated" if report["complete"] else "INCOMPLETE physical-x scout; no threshold fits"
            if report.get("synthetic_preview"):
                status = "SYNTHETIC LAYOUT PREVIEW - not production evidence"
            fig.text(.5,.94,status,ha="center",fontsize=11.5,color=".25" if report["complete"] else "#A23B32")
            fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.5,.921),ncol=len(sizes),frameon=False,fontsize=10)
            for ri,arm in enumerate(ARMS):
                rates = selection["channels"][ch]["rates"]
                for ci,(direction,data) in enumerate((("Diagonal growth",reference_rows),("Physical x-axis growth",rows))):
                    axis = axes[ri,ci]
                    plot_observed(axis,data,ch,arm,selection["channels"][ch]["sizes"])
                    axis.set_title(f"{direction} | {LABELS[arm]}",fontsize=11.2,pad=8)
                    axis.set_xlim(0,100*max(rates)*1.035); axis.set_ylim(0,1)
                    axis.set_xlabel("Gate-output fault probability p (%)",fontsize=10)
                    axis.tick_params(labelsize=9)
                    if ci==0:
                        axis.set_ylabel("Either-logical failure probability",fontsize=10)
                    if not any(r["channel"]==ch and r["shots"] for r in data):
                        axis.text(.5,.5,"Awaiting sampled records",ha="center",transform=axis.transAxes,color=".45")
            footer = (
                f"T = Ly = Lz = L. Axis records: {report['completed_fault_records']:,}/{report['target_fault_records']:,}; target {selection['shots']:,} per (L,p). 95% Wilson bars; open markers: unfinished points.\n"
                "Square seed versus diagonal seed: 2 L^2 versus 3 L^2 qubits; d_Z = L versus 2L. Both use 13 L^3 noisy transfer outputs and d_X = (L+1)L.\n"
                "X seeds are ideal. Z seed-H faults and common seed repair occur in both directions; undetectable logical seed errors remain.\n"
                "Axis Z one-step: deterministic +x defect push. Diagonal Z one-step: MWPM endpoint routing. These policies and code distances differ.\n"
                "Arms are paired within a direction; axis minus diagonal is independent-sample. Two-slice X is preparation damage, not a memory threshold.\n"
                "Orientation scout only: no axis threshold fit validated; historical diagonal fit windows are not assumed valid. No missing samples are imputed.")
            fig.text(.06,.035,footer,fontsize=8.1,va="bottom",linespacing=1.55)
            fig.text(.975,.023,f"{page}/2",ha="right",fontsize=9)
            pdf.savefig(fig); plt.close(fig)


def export(study,rows,reference_rows,selection,report):
    data,figures = study/"data",study/"figures"
    data.mkdir(parents=True,exist_ok=True); figures.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".direction-export-",dir=data) as temporary:
        temporary = Path(temporary)
        csv_path = temporary/"matched-results.csv"
        with csv_path.open("w",newline="") as stream:
            writer = csv.DictWriter(stream,fieldnames=COLUMNS,extrasaction="ignore")
            writer.writeheader(); writer.writerows(rows)
        pdf_path = temporary/"h6-direction-comparison.pdf"
        make_pdf(rows,reference_rows,selection,report,pdf_path,style_path=study/"archive/paper-no-grid.mplstyle")
        report["output_sha256"] = {"data/matched-results.csv":digest(csv_path.read_bytes()),"figures/h6-direction-comparison.pdf":digest(pdf_path.read_bytes())}
        report_path = temporary/"report.json"
        report_path.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
        os.replace(csv_path,data/csv_path.name); os.replace(pdf_path,figures/pdf_path.name)
        os.replace(report_path,data/report_path.name)


def main(arguments=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study",required=True,type=Path)
    parser.add_argument("--progress",action="store_true")
    args = parser.parse_args(arguments)
    try:
        study = args.study.resolve()
        selection_bytes = (study/"data/selection.json").read_bytes()
        selection = validate_selection(json.loads(selection_bytes))
        sources = verify_frozen_sources(study,selection)
        reference_rows = load_reference(selection)
        rows,provenance,runtime,skipped = load_inputs(study,selection,sources,progress=args.progress)
        require_complete(rows,selection,progress=args.progress)
        report = make_report(rows,reference_rows,selection,provenance,runtime,digest(selection_bytes),skipped=skipped)
        export(study,rows,reference_rows,selection,report)
        print(json.dumps({key:report[key] for key in ("status","completed_points","expected_points","completed_fault_records","target_fault_records")},sort_keys=True))
        return 0
    except (ValueError,KeyError,TypeError,OSError,ActiveSnapshot) as error:
        print(f"direction analysis rejected input: {error}",file=sys.stderr)
        return 2


if __name__=="__main__":
    raise SystemExit(main())
