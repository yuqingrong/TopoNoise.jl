"""Independent, read-only audit of the frozen follow-up simulation artifacts.

The default audit requires every production point to be complete. --design-only
checks the frozen design, source and preservation bindings, and RNG independence;
it never represents incomplete production as a completed audit. No simulation
module or analysis module is imported. Julia is used only for Random/TOML seed
regeneration, preserving Julia's exact RNG version and scalar draw semantics.
"""
import argparse
from collections import Counter
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import re
import subprocess
import time
import tomllib


HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
ROOT = STUDY.parent.parent
OLD = ROOT / "results/toric3d-h6-full-prefix-2026-10-15"
REANALYSIS = ROOT / "results/toric3d-h6-full-prefix-reanalysis-2026-10-16"
EXPECTED_RECORDS = {"x_only": 106000, "z_only": 98000}
ARMS = ("final_only", "full_prefix")


class AuditError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise AuditError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def file_digest(path):
    require(path.is_file(), f"missing file: {path}")
    return digest(path.read_bytes())


def point_grid(config):
    return [(config["channel"], size, p) for size in config["sizes"] for p in config["rates"]]


def validate_selection(selection, jobs):
    require(jobs == selection["jobs"], "jobs.toml differs from frozen selection")
    require(len(jobs) == 13, "expected exactly 13 production shards")
    require(len({job["name"] for job in jobs}) == 13, "duplicate shard names")
    require(len({job["seed"] for job in jobs}) == 13, "duplicate production base seeds")
    counts = Counter()
    cells = []
    for job in jobs:
        require(job["channel"] in EXPECTED_RECORDS, "unrecognized channel")
        require(job["worker"] in (1, 2), "job assigned outside the two frozen workers")
        for key in ("sizes", "rates"):
            values = job[key]
            require(values and values == sorted(set(values)), f"{job['name']}: invalid {key}")
        require(job["shots"] > 0 and job["batch_size"] > 0 and job["seed"] >= 0, "invalid positive configuration")
        grid = point_grid(job)
        cells.extend(grid)
        counts[job["channel"]] += len(grid) * job["shots"]
    require(len(cells) == len(set(cells)) == 77, "expected 77 unique production cells")
    require(dict(counts) == EXPECTED_RECORDS == selection["new_fault_records"], "frozen record totals differ")
    require(sum(counts.values()) == selection["total_new_fault_records"] == 204000, "total record count differs")
    for channel, key in (("x_only", "x"), ("z_only", "z")):
        planned = selection["analysis"][key]
        expected = {(channel, size, p) for size in planned["sizes"] for p in planned["rates"]}
        require({cell for cell in cells if cell[0] == channel} == expected, f"{key} analysis grid differs from jobs")
    return dict(counts)


def parse_checkpoint(raw, label):
    matches = re.findall(rb'^payload_sha256 = "([0-9a-f]{64})"\n', raw, re.MULTILINE)
    require(len(matches) == 1, f"{label}: noncanonical or missing payload checksum")
    payload = re.sub(rb'^payload_sha256 = "[0-9a-f]{64}"\n', b"", raw, count=1, flags=re.MULTILINE)
    require(digest(payload) == matches[0].decode(), f"{label}: checkpoint payload checksum mismatch")
    checkpoint = tomllib.loads(raw.decode())
    require(checkpoint["format_version"] == 1, f"{label}: unknown checkpoint version")
    return checkpoint


def validate_state(state, target, label):
    n = state["completed_shots"]
    require(type(n) is int and 0 <= n <= target, f"{label}: invalid sample count")
    for key, length in (("logical_1_counts", 2), ("logical_2_counts", 2), ("failures", 2), ("joint_counts", 4)):
        values = state[key]
        require(len(values) == length and all(type(value) is int and 0 <= value <= n for value in values),
                f"{label}: invalid {key}")
    require(sum(state["joint_counts"]) == n, f"{label}: histogram total differs from sample count")
    for arm in range(2):
        first, second, failure = (state[key][arm] for key in ("logical_1_counts", "logical_2_counts", "failures"))
        require(max(first, second) <= failure <= min(n, first + second), f"{label}: invalid logical union")
        marginal = sum(count for mask, count in enumerate(state["joint_counts"]) if mask & (1 << arm))
        require(marginal == failure, f"{label}: histogram marginal differs")
    require(math.isfinite(state["elapsed_seconds"]) and state["elapsed_seconds"] >= 0, f"{label}: invalid timing")


def validate_raw(rows, checkpoint, label):
    config = checkpoint["config"]
    grid = point_grid(config)
    require(len(rows) == len(grid), f"{label}: raw row count differs")
    for index, (row, (channel, size, p), state) in enumerate(zip(rows, grid, checkpoint["states"])):
        prefix = f"{label} row {index + 1}"
        require(row["channel"] == channel and all(int(row[key]) == size for key in ("steps", "Ly", "Lz"))
                and float(row["p"]) == p, f"{prefix}: grid differs")
        require(row["complete"] == "true" and int(row["shots"]) == config["shots"] == state["completed_shots"]
                and int(row["target_shots"]) == config["shots"], f"{prefix}: incomplete point")
        require(int(row["seed"]) == config["seed"] and row["point_seed"] == checkpoint["point_seeds_hex"][index],
                f"{prefix}: seed differs")
        require(float(row["elapsed_seconds"]) == state["elapsed_seconds"], f"{prefix}: timing differs")
        for arm, name in enumerate(ARMS):
            for column, key in (("logical_1_failures", "logical_1_counts"), ("logical_2_failures", "logical_2_counts"), ("failures", "failures")):
                require(int(row[f"{name}_{column}"]) == state[key][arm], f"{prefix}: {name} {column} differs")
        require([int(row[f"joint_{mask:02d}"]) for mask in range(4)] == state["joint_counts"], f"{prefix}: raw histogram differs")


def validate_seed_sets(new_seeds, historical_seeds):
    require(all(re.fullmatch(r"[0-9a-f]{16}", seed) for seed in new_seeds), "invalid point seed encoding")
    require(len(set(new_seeds)) == len(new_seeds), "duplicate new point seed")
    require(not (set(new_seeds) & set(historical_seeds)), "new point seed reuses a historical or benchmark stream")


def regenerate_seeds(julia):
    # A separate implementation of the documented scalar MersenneTwister draws.
    program = r'''
using Random, TOML
for job in TOML.parsefile(ARGS[1])["jobs"]
    rng = MersenneTwister(job["seed"])
    count = length(job["sizes"]) * length(job["rates"])
    seen = Set{UInt64}(); seeds = String[]
    while length(seeds) < count
        value = rand(rng, UInt64)
        value in seen && continue
        push!(seen, value); push!(seeds, string(value; base=16, pad=16))
    end
    println(job["name"], " ", join(seeds, ","))
end
'''
    completed = subprocess.run([julia, "--startup-file=no", "--history-file=no", "-e", program, str(STUDY / "jobs.toml")],
                               check=True, text=True, capture_output=True, timeout=60)
    return {line.split()[0]: line.split()[1].split(",") for line in completed.stdout.splitlines() if line.strip()}


def preservation(selection):
    manifest_path = OLD / "ARTIFACTS.sha256.json"
    require(file_digest(manifest_path) == selection["old_manifest_sha256"], "historical artifact manifest changed")
    manifest = json.loads(manifest_path.read_text())
    for name, expected in manifest.items():
        require(file_digest(OLD / name) == expected, f"historical file changed: {name}")
    for name, expected in selection["prior_analysis_report_sha256"].items():
        require(file_digest(REANALYSIS / name) == expected, f"prior reanalysis changed: {name}")
    prior = json.loads((OLD / "verification/prior-study-integrity.json").read_text())
    require(file_digest(ROOT / "document/note.tex") == prior["note_sha256"], "document/note.tex changed")
    return {"historical_files_verified": len(manifest), "prior_reanalysis_reports_verified": 2, "note_unchanged": True}


def validate_archive(directory, checkpoint):
    snapshots = sorted((directory / "checkpoints").glob("checkpoint-*.toml"))
    require(len(snapshots) == checkpoint["generation"], f"{directory.name}: missing or additional snapshots")
    previous = None
    for generation, path in enumerate(snapshots, 1):
        snap = parse_checkpoint(path.read_bytes(), str(path))
        require(snap["generation"] == generation, f"{path}: archive generation gap")
        expected_name = f"checkpoint-{generation:06d}-{snap['payload_sha256'][:12]}.toml"
        require(path.name == expected_name, f"{path}: snapshot filename differs")
        for key in ("config", "runtime", "source_hashes", "protocol", "point_seeds_hex", "started_at_unix"):
            require(snap[key] == checkpoint[key], f"{path}: immutable {key} changed")
        require(len(snap["states"]) == len(checkpoint["states"]), f"{path}: snapshot state count differs")
        for index, state in enumerate(snap["states"]):
            validate_state(state, checkpoint["config"]["shots"], f"{path} state {index}")
            if previous is not None:
                prior = previous["states"][index]
                require(state["completed_shots"] >= prior["completed_shots"] and state["elapsed_seconds"] >= prior["elapsed_seconds"],
                        f"{path}: sample/timing regression")
                for key in ("logical_1_counts", "logical_2_counts", "failures", "joint_counts"):
                    require(all(a >= b for a, b in zip(state[key], prior[key])), f"{path}: count regression")
        if previous is not None:
            require(snap["updated_at_unix"] >= previous["updated_at_unix"], f"{path}: time regression")
        previous = snap
    require(snapshots[-1].read_bytes() == (directory / "checkpoint.toml").read_bytes(), f"{directory}: final snapshot differs")
    return len(snapshots)


def audit_shard(job, selection, regenerated, historical_runtime, historical_protocol):
    directory = STUDY / job["name"]
    require(directory.is_dir(), f"missing shard: {job['name']}")
    raw_bytes = (directory / "raw.csv").read_bytes()
    checkpoint_bytes = (directory / "checkpoint.toml").read_bytes()
    checkpoint = parse_checkpoint(checkpoint_bytes, job["name"])
    metadata = tomllib.loads((directory / "metadata.toml").read_text())
    config = {key: job[key] for key in ("channel", "sizes", "rates", "shots", "seed", "batch_size")}
    config["dimensions"] = [[size, size, size] for size in job["sizes"]]
    require(checkpoint["config"] == config, f"{job['name']}: configuration differs")
    require(checkpoint["source_hashes"] == selection["simulation_source_sha256"], f"{job['name']}: source differs")
    require(checkpoint["runtime"] == historical_runtime, f"{job['name']}: runtime differs from established engine")
    require(checkpoint["protocol"] == historical_protocol[job["channel"]], f"{job['name']}: protocol differs")
    require(checkpoint["point_seeds_hex"] == regenerated[job["name"]], f"{job['name']}: exact seed regeneration differs")
    require(checkpoint["started_at_unix"] >= selection["frozen_before_production_at_unix"], f"{job['name']}: started before selection was frozen")
    require(len(checkpoint["states"]) == len(point_grid(config)), f"{job['name']}: state count differs")
    for index, state in enumerate(checkpoint["states"]):
        validate_state(state, job["shots"], f"{job['name']} state {index}")
    for key in ("format_version", "config", "runtime", "source_hashes", "protocol", "point_seeds_hex", "generation",
                "decoder_build_seconds_by_size", "decoder_build_count_by_size", "started_at_unix", "updated_at_unix"):
        require(metadata[key] == checkpoint[key], f"{job['name']}: metadata {key} differs")
    npoints = len(checkpoint["states"])
    require(metadata["complete"] is True and metadata["completed_points"] == metadata["expected_points"] == npoints
            and metadata["completed_shots"] == npoints * job["shots"], f"{job['name']}: production is incomplete")
    require(metadata["raw_sha256"] == digest(raw_bytes), f"{job['name']}: raw checksum differs")
    require(metadata["checkpoint_sha256"] == digest(checkpoint_bytes), f"{job['name']}: checkpoint file checksum differs")
    expected_snapshot = f"checkpoints/checkpoint-{checkpoint['generation']:06d}-{checkpoint['payload_sha256'][:12]}.toml"
    require(metadata["current_snapshot"] == expected_snapshot, f"{job['name']}: metadata snapshot differs")
    validate_raw(list(csv.DictReader(io.StringIO(raw_bytes.decode()))), checkpoint, job["name"])
    snapshots = validate_archive(directory, checkpoint)
    return {"name": job["name"], "channel": job["channel"], "points": npoints, "records": npoints * job["shots"],
            "snapshots_verified": snapshots, "raw_sha256": digest(raw_bytes), "checkpoint_sha256": digest(checkpoint_bytes),
            "metadata_sha256": file_digest(directory / "metadata.toml")}


def audit(julia="julia", design_only=False):
    selection_bytes = (STUDY / "selection.json").read_bytes()
    selection = json.loads(selection_bytes)
    claimed = (STUDY / "selection.sha256").read_text().split()[0]
    require(digest(selection_bytes) == claimed, "frozen selection checksum differs")
    jobs = tomllib.loads((STUDY / "jobs.toml").read_text())["jobs"]
    totals = validate_selection(selection, jobs)
    for name, expected in selection["simulation_source_sha256"].items():
        require(file_digest(ROOT / name) == expected, f"simulation source changed: {name}")
    expected_sources = {"Project.toml", "results/toric3d-h6-full-prefix-2026-10-15/run_study.jl"}
    if (ROOT / "Manifest.toml").exists():
        expected_sources.add("Manifest.toml")
    expected_sources.update(str(path.relative_to(ROOT)) for path in (ROOT / "src").rglob("*.jl"))
    require(set(selection["simulation_source_sha256"]) == expected_sources, "source manifest does not cover exact current runner source set")
    preserved = preservation(selection)
    prior_paths = sorted(OLD.glob("*/checkpoint.toml"))
    prior_paths += [STUDY / name / "checkpoint.toml" for name in selection["excluded_timing_runs"]]
    prior_seeds = set()
    for path in prior_paths:
        prior = parse_checkpoint(path.read_bytes(), str(path))
        prior_seeds.update(prior["point_seeds_hex"])
    regenerated = regenerate_seeds(julia)
    require(set(regenerated) == {job["name"] for job in jobs}, "seed regeneration returned wrong shard set")
    new_seeds = [seed for job in jobs for seed in regenerated[job["name"]]]
    require(len(new_seeds) == 77, "seed regeneration returned wrong point count")
    validate_seed_sets(new_seeds, prior_seeds)
    failures = list(STUDY.rglob("failed-record.toml"))
    require(not failures, "fatal simulation record dump exists: " + ", ".join(map(str, failures)))
    report = {"status": "design_verified" if design_only else "passed", "checked_at_unix": time.time(),
              "audit_scope": "frozen design only" if design_only else "complete production and all immutable snapshots",
              "selection_sha256": claimed, "planned_shards": 13, "planned_points": 77, "planned_records": totals,
              "simulation_source_files_verified": len(expected_sources), "preservation": preserved,
              "independent_seed_regeneration": True, "new_point_seeds": len(new_seeds),
              "excluded_checkpoint_streams": len(prior_paths), "excluded_unique_point_seeds": len(prior_seeds),
              "seed_sets_disjoint": True, "failed_record_files": 0}
    if design_only:
        return report
    names = {job["name"] for job in jobs}
    require({path.name for path in STUDY.glob("production-*") if path.is_dir()} <= names, "unexpected production shard directory")
    protocols = {}
    runtime = None
    for channel, name in (("x_only", "production-x-small"), ("z_only", "production-z")):
        historical = tomllib.loads((OLD / name / "metadata.toml").read_text())
        protocols[channel] = historical["protocol"]
        if runtime is None:
            runtime = historical["runtime"]
        require(runtime == historical["runtime"], "historical channel runtimes differ")
    shards, errors = [], []
    for job in jobs:
        try:
            shards.append(audit_shard(job, selection, regenerated, runtime, protocols))
        except (AuditError, OSError, KeyError, ValueError) as error:
            errors.append(f"{job['name']}: {error}")
    report["shards"] = shards
    report["completed_shards_verified"] = len(shards)
    report["completed_points_verified"] = sum(shard["points"] for shard in shards)
    report["completed_records_verified"] = dict(Counter({channel: sum(shard["records"] for shard in shards if shard["channel"] == channel)
                                                        for channel in EXPECTED_RECORDS}))
    report["errors"] = errors
    if errors or len(shards) != 13 or report["completed_points_verified"] != 77 or report["completed_records_verified"] != EXPECTED_RECORDS:
        report["status"] = "failed"
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--julia", default="julia")
    parser.add_argument("--design-only", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        report = audit(args.julia, args.design_only)
    except (AuditError, OSError, KeyError, ValueError, subprocess.SubprocessError) as error:
        report = {"status": "failed", "checked_at_unix": time.time(), "errors": [str(error)]}
    rendered = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end="")
    return 0 if report["status"] in ("passed", "design_verified") else 1


if __name__ == "__main__":
    raise SystemExit(main())
