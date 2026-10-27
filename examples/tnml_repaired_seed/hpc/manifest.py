#!/usr/bin/env python3
"""Freeze a portable, explicitly enumerated HPC study after environment setup.

This module only prepares and audits files. It never launches Julia or submits
Slurm jobs. Python >= 3.11 is required; use the bootstrap's managed Python.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import importlib.util
import json
import math
from pathlib import Path, PurePosixPath
import sys

if sys.version_info < (3, 11):
    raise SystemExit("Use the bootstrapped Python >= 3.11, not the system Python 3.9.")

HERE = Path(__file__).resolve().parent
EXPERIMENT = HERE.parent
REPOSITORY = HERE.parents[2]
_SPEC = importlib.util.spec_from_file_location("hpc_existing_auditor", EXPERIMENT / "run_full_grid.py")
DRIVER = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(DRIVER)

STAGES = ("zero", "smoke", "feasibility", "production")
STAGE_LIMITS = {
    "zero": dict(per_contraction_seconds=1800, per_run_total_seconds=7200, worker_rss_gib=8),
    "smoke": dict(per_contraction_seconds=300, per_run_total_seconds=3600, worker_rss_gib=8),
    "feasibility": dict(per_contraction_seconds=1800, per_run_total_seconds=7200, worker_rss_gib=16),
    "production": dict(per_contraction_seconds=600, per_run_total_seconds=21600, worker_rss_gib=8),
}
DEFAULT_CHI = 240
DEFAULT_CUTOFF = 1e-26
DEFAULT_SEED_BASE = 202610500000


def scientific_source_hashes():
    hashes = DRIVER.scientific_source_hashes()
    if not {"Project.toml", "Manifest.toml"}.issubset(hashes):
        raise ValueError("Instantiate the cluster environment before freezing Project.toml and Manifest.toml")
    return hashes


def wrapper_hashes():
    paths = [EXPERIMENT / name for name in ("run_full_grid.py", "supervise.py", "plot.py")]
    # The numerical runner's source identity does not include Python setup.
    # Freeze that configuration here without changing its existing metadata.
    paths.append(REPOSITORY / "CondaPkg.toml")
    python_marker = REPOSITORY / ".hpc-python-path"
    if python_marker.is_file():
        paths.append(python_marker)
    paths += [p for p in HERE.rglob("*") if p.is_file() and p.suffix in (".py", ".sh", ".sbatch", ".jl")
              and not p.name.startswith("test_") and "__pycache__" not in p.parts]
    return {str(p.relative_to(REPOSITORY)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(set(paths))}


def manifest_identity(document):
    return DRIVER.document_hash({k: v for k, v in document.items() if k != "manifest_id"})


def safe_relative(value):
    if not isinstance(value, str) or not value or "\\" in value or "\0" in value:
        raise ValueError("task output and frozen file paths must be safe relative paths")
    path = PurePosixPath(value)
    if path.is_absolute() or ".." in path.parts or str(path) != value or value == ".":
        raise ValueError("task output and frozen file paths must be safe relative paths")
    return path


def output_directory(root, task):
    root = Path(root).resolve()
    output = (root / safe_relative(task["output"])).resolve()
    if not output.is_relative_to(root) or output == root:
        raise ValueError("task output resolves outside the study root")
    return output


def _settings(chi, cutoff, seed_base):
    if isinstance(chi, bool) or not isinstance(chi, int) or chi <= 0:
        raise ValueError("chi must be a positive integer")
    if isinstance(cutoff, bool) or not isinstance(cutoff, (int, float)) or not math.isfinite(cutoff) or not 0 <= cutoff < 1:
        raise ValueError("cutoff must be finite and in [0,1)")
    if isinstance(seed_base, bool) or not isinstance(seed_base, int) or not 0 <= seed_base < 2**63 - 500000:
        raise ValueError("seed base must be nonnegative and leave room for distinct Int64 stage seeds")


def _stage_plan(chi, cutoff, seed_base):
    _settings(chi, cutoff, seed_base)
    stages = {name: {"resource_limits": dict(STAGE_LIMITS[name]), "tasks": []} for name in STAGES}

    def add(stage, L, p, shots, dimension, seed, **extra):
        tasks = stages[stage]["tasks"]
        point_id = f"L{L}-p{p:.4f}-chi{dimension}"
        output = f"{stage}/{point_id}"
        if "shard_id" in extra:
            output += f"/shard-{extra['shard_id']:03d}"
        tasks.append(dict(task_id=len(tasks), stage=stage, phase=stage, L=L, p=p,
                          chi=dimension, cutoff=cutoff, shots=shots, seed=seed,
                          exact=L == 3, output=output, point_id=point_id, **extra))

    for index, L in enumerate((3, 5, 7, 9, 11), 1):
        add("zero", L, 0.0, 1, chi, seed_base + 100000 + index)
    for index, p in enumerate((.006, .01, .02), 1):
        add("smoke", 3, p, 2, chi, seed_base + 200000 + index)
        for dimension in (chi, 2 * chi):
            add("feasibility", 5, p, 3, dimension, seed_base + 300000 + index,
                pair_id=f"L5-p{p:.4f}")
    for p in DRIVER.POSITIVE_RATES:
        for shard_id in range(20):
            task_id = len(stages["production"]["tasks"])
            add("production", 3, p, 50, chi, seed_base + 400000 + task_id + 1,
                shard_id=shard_id)
    return stages


def validate_manifest(document):
    if document.get("schema_version") != 1 or document.get("study_kind") != "tnml_hpc_shards":
        raise ValueError("unsupported HPC manifest schema")
    if document.get("manifest_id") != manifest_identity(document):
        raise ValueError("frozen manifest identity differs; do not edit a running study")
    chi, cutoff, seed_base = document["chi"], document["cutoff"], document["seed_base"]
    _settings(chi, cutoff, seed_base)
    if document.get("negative_tolerance") != 1e-10:
        raise ValueError("the existing runner fixes negative tolerance at 1e-10")
    if (document.get("target_shots_per_point"), document.get("shards_per_point"), document.get("shots_per_shard")) != (1000, 20, 50):
        raise ValueError("the production plan must have 20 independent 50-shot shards per point")
    stages = document.get("stages", {})
    if set(stages) != set(STAGES):
        raise ValueError("manifest must explicitly contain all four stages")
    outputs, production_seeds = set(), set()
    for name in STAGES:
        for index, task in enumerate(stages[name].get("tasks", [])):
            if isinstance(task.get("task_id"), bool) or task.get("task_id") != index or task.get("stage") != name or task.get("phase") != name:
                raise ValueError("task IDs must be ordered, unique integers within each stage")
            output = str(safe_relative(task.get("output")))
            if output in outputs or PurePosixPath(output).parts[0] != name:
                raise ValueError("duplicate or incorrectly staged task output")
            outputs.add(output)
            if name == "production":
                if task.get("seed") in production_seeds:
                    raise ValueError("production seeds must be unique")
                production_seeds.add(task.get("seed"))
    if stages != _stage_plan(chi, cutoff, seed_base):
        raise ValueError("task rows or resource limits differ from the frozen study design")
    for key in ("scientific_source_sha256", "wrapper_sha256"):
        hashes = document.get(key)
        if not isinstance(hashes, dict) or not hashes:
            raise ValueError(f"missing frozen {key}")
        for path, digest in hashes.items():
            safe_relative(path)
            if not isinstance(digest, str) or len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
                raise ValueError(f"invalid SHA256 in {key}")
    return document


def load_manifest(root, verify_sources=True):
    document = validate_manifest(json.loads((Path(root) / "manifest.json").read_text()))
    if verify_sources:
        if document["scientific_source_sha256"] != scientific_source_hashes():
            raise ValueError("scientific sources or dependency lockfiles changed after manifest freeze")
        if document["wrapper_sha256"] != wrapper_hashes():
            raise ValueError("HPC wrappers changed after manifest freeze")
    return document


def create_manifest(root, chi=None, cutoff=None, seed_base=None):
    """Create once after setup; identical calls only verify and never rewrite."""
    root = Path(root)
    if not root.is_absolute():
        raise ValueError("study root must be an absolute path")
    root = root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    with (root / ".manifest.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        path = root / "manifest.json"
        if path.exists():
            document = load_manifest(root)
            requested = dict(chi=chi, cutoff=cutoff, seed_base=seed_base)
            if any(value is not None and document[key] != value for key, value in requested.items()):
                raise ValueError("settings differ from the frozen manifest; create a new study")
            return document
        chi = DEFAULT_CHI if chi is None else chi
        cutoff = DEFAULT_CUTOFF if cutoff is None else cutoff
        seed_base = DEFAULT_SEED_BASE if seed_base is None else seed_base
        stages = _stage_plan(chi, cutoff, seed_base)
        document = dict(
            schema_version=1, study_kind="tnml_hpc_shards",
            frozen_at_utc=datetime.now(timezone.utc).isoformat(),
            chi=chi, cutoff=cutoff, negative_tolerance=1e-10, seed_base=seed_base,
            target_shots_per_point=1000, shards_per_point=20, shots_per_shard=50,
            p_values=[0.0] + DRIVER.POSITIVE_RATES, requested_sizes=[3, 5, 7, 9, 11], production_sizes=[3],
            cpus_per_task=1, julia_threads=1, blas_threads=1,
            physical_preparation="Noisy seed-H outputs; complete initial seed MWPM repair; diagonal H+6CNOT transfer with independent Z gate-output faults; no later feedback.",
            decoder_information="TN conditions on seed and final syndromes plus the known seed recovery; paired MWPM uses the final syndrome only.",
            seed_policy="Distinct stage seed ranges. Only labeled feasibility pairs share a seed across chi. Production shards use distinct fixed Random.Xoshiro seeds.",
            unknown_policy="Every planned production point retains all 1,000 requested outcomes. Unstarted, missing, invalid, or unvalidated outcomes are never discarded or treated as successes.",
            zero_policy="The p=0 failure rate is analytically zero. Single-history checks are implementation regressions, not statistical trials.",
            production_gate="All three manifested L3 smoke tasks must complete with independent exact-reference audits before production dispatch.",
            scientific_source_sha256=scientific_source_hashes(), wrapper_sha256=wrapper_hashes(), stages=stages)
        document["manifest_id"] = manifest_identity(document)
        validate_manifest(document)
        DRIVER.atomic_json(path, document)
        return document


def stage_tasks(document, stage):
    if stage not in STAGES:
        raise ValueError(f"unknown stage: {stage}")
    return document["stages"][stage]["tasks"]


def task_for_id(document, stage, task_id):
    tasks = stage_tasks(document, stage)
    if isinstance(task_id, bool) or not isinstance(task_id, int) or not 0 <= task_id < len(tasks):
        raise ValueError(f"task ID is outside stage {stage}")
    return tasks[task_id]


def point_for_audit(task):
    return {name: task[name] for name in ("phase", "L", "p", "chi", "cutoff", "shots", "seed", "exact", "output")}


def validation_selection(document):
    return {name: document[name] for name in ("scientific_source_sha256", "cutoff", "negative_tolerance")}


def validate_task_metadata(metadata, task, document):
    if task != task_for_id(document, task["stage"], task["task_id"]):
        raise ValueError("task differs from its immutable manifest row")
    DRIVER.validate_metadata_identity(metadata, point_for_audit(task), validation_selection(document))


def audit_task(root, task, document):
    if task != task_for_id(document, task["stage"], task["task_id"]):
        raise ValueError("task differs from its immutable manifest row")
    output_directory(root, task)
    return DRIVER.audit_point(Path(root), point_for_audit(task), validation_selection(document))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    create = commands.add_parser("create")
    create.add_argument("--root", type=Path, required=True)
    create.add_argument("--chi", type=int)
    create.add_argument("--cutoff", type=float)
    create.add_argument("--seed-base", type=int)
    args = parser.parse_args()
    try:
        document = create_manifest(args.root, args.chi, args.cutoff, args.seed_base)
    except (OSError, ValueError, KeyError) as error:
        parser.error(str(error))
    print(json.dumps({"manifest": str(args.root.resolve() / "manifest.json"), "manifest_id": document["manifest_id"],
                      "tasks": {stage: len(stage_tasks(document, stage)) for stage in STAGES}}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
