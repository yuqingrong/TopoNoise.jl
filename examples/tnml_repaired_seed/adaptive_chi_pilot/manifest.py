#!/usr/bin/env python3
"""Immutable phase-1 manifest for the adaptive-chi pilot."""

from __future__ import annotations

import json
import os
from pathlib import Path, PurePosixPath
from typing import Any, Mapping

try:  # Package import.
    from . import core
except ImportError:  # Direct script/test invocation from this directory.
    import core


SCHEMA_VERSION = 1
MANIFEST_NAME = "manifest.json"
P_VALUES = core.P_STRATA
L3_LADDER = (32, 64, 128, 240, 480, 960)
L5_LADDER = (32, 64, 128, 240, 480)
CUTOFF = 1.0e-18


def _seed_base(value: Any) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError("seed_base must be a nonnegative integer")
    if value + 1_000_059 >= 2**63:
        raise ValueError("seed_base is too large")
    return value


def _task(
    *,
    task_id: int,
    stage: str,
    L: int,
    p_key: str,
    shots: int,
    chi_ladder: tuple[int, ...],
    exact: bool,
    history_seeds: list[int],
    output: str,
    history_offset: int,
) -> dict[str, Any]:
    return {
        "task_id": task_id,
        "stage": stage,
        "L": L,
        "p": float(p_key),
        "p_key": p_key,
        "shots": shots,
        "history_offset": history_offset,
        "chi_ladder": list(chi_ladder),
        "cutoff": CUTOFF,
        "exact": exact,
        "seed": history_seeds[0],
        "history_seeds": history_seeds,
        "output": output,
    }


def _enumerate_l3(seed_base: int) -> list[dict[str, Any]]:
    tasks = []
    task_id = 0
    for shard in range(5):
        for p_index, p_key in enumerate(P_VALUES):
            history_offset = p_index * 100 + shard * 20
            tasks.append(
                _task(
                    task_id=task_id,
                    stage="l3_exact",
                    L=3,
                    p_key=p_key,
                    shots=20,
                    chi_ladder=L3_LADDER,
                    exact=True,
                    history_seeds=list(
                        range(seed_base + history_offset, seed_base + history_offset + 20)
                    ),
                    output=f"l3_exact/L3-p{p_key}/shard-{shard:03d}",
                    history_offset=history_offset,
                )
            )
            task_id += 1
    return tasks


def _enumerate_l5(seed_base: int) -> list[dict[str, Any]]:
    tasks = []
    task_id = 0
    for history in range(10):
        for p_index, p_key in enumerate(P_VALUES):
            history_offset = p_index * 10 + history
            tasks.append(
                _task(
                    task_id=task_id,
                    stage="l5_feasibility",
                    L=5,
                    p_key=p_key,
                    shots=1,
                    chi_ladder=L5_LADDER,
                    exact=False,
                    history_seeds=[seed_base + 1_000_000 + history_offset],
                    output=f"l5_feasibility/L5-p{p_key}/shot-{history:03d}",
                    history_offset=history_offset,
                )
            )
            task_id += 1
    return tasks


def build_manifest(seed_base: int) -> dict[str, Any]:
    """Build the deterministic document without touching the filesystem."""
    seed = _seed_base(seed_base)
    document: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "kind": "adaptive_chi_phase1_manifest",
        "seed_base": seed,
        "p_values": list(P_VALUES),
        "cutoff": CUTOFF,
        "stopping_policy": {
            "minimum_completed_levels": 3,
            "posterior_bound": core.POSTERIOR_TOLERANCE,
            "log_joint_bound": core.LOG_JOINT_TOLERANCE,
            "normalization_tolerance": core.NORMALIZATION_TOLERANCE,
            "negative_probability_tolerance": core.NEGATIVE_PROBABILITY_TOLERANCE,
            "minimum_margin": core.MINIMUM_MARGIN,
            "margin_multiplier": 4.0,
            "drift_floor": core.DRIFT_FLOOR,
            "drift_growth": core.DRIFT_GROWTH,
        },
        "selection": {
            "version": "tnml-adaptive-chi-v1",
            "message": (
                "tnml-adaptive-chi-v1|L=<decimal L>|p=<four-decimal p>|"
                "record_id=<lowercase hex>"
            ),
            "digest": "sha256",
            "digest_integer_order": "unsigned_big_endian",
            "tie_break": "record_id_lexical",
            "l3_replay_records_per_p": 5,
            "l5_forced_anchors_per_p": 4,
        },
        "runtime_policy": {
            "execution": "disabled_phase1",
            "requested_cpus_per_task": 4,
            "julia_application_threads": 1,
            "blas_threads": 4,
            "max_rss_bytes": 16 * 2**30,
            "one_thread_replay_required": True,
            "exact_wall_seconds": 300,
            "l3_wall_seconds_by_chi": {
                "32": 120,
                "64": 300,
                "128": 900,
                "240": 1800,
                "480": 3600,
                "960": 7200,
            },
            "l5_wall_seconds_by_chi": {
                "32": 120,
                "64": 300,
                "128": 900,
                "240": 1800,
                "480": 3600,
            },
        },
        "planned_batches": {
            "l3_exact": [
                {"batch": 0, "first_task_id": 0, "last_task_id": 29,
                 "max_concurrent": 20}
            ],
            "l5_feasibility": [
                {"batch": 0, "first_task_id": 0, "last_task_id": 29,
                 "max_concurrent": 20},
                {"batch": 1, "first_task_id": 30, "last_task_id": 59,
                 "max_concurrent": 20},
            ],
        },
        "stages": {
            "l3_exact": {
                "description": "exact calibration; 100 histories per p in 20-shot shards",
                "L": 3,
                "histories_per_p": 100,
                "tasks": _enumerate_l3(seed),
            },
            "l5_feasibility": {
                "description": "runtime/escalation feasibility smoke; ten one-shot tasks per p",
                "L": 5,
                "histories_per_p": 10,
                "forced_anchors_per_p": 4,
                "classification": "feasibility_smoke",
                "tasks": _enumerate_l5(seed),
            },
        },
    }
    document["manifest_id"] = core.identity_hash(document, "manifest_id")
    return document


def _safe_output(output: Any) -> bool:
    if not isinstance(output, str) or not output:
        return False
    path = PurePosixPath(output)
    return not path.is_absolute() and ".." not in path.parts and "." not in path.parts


def validate_manifest(document: Mapping[str, Any]) -> dict[str, Any]:
    """Validate both the document hash and exact deterministic task enumeration."""
    result = dict(document)
    manifest_id = result.get("manifest_id")
    if not isinstance(manifest_id, str) or manifest_id != core.identity_hash(result, "manifest_id"):
        raise ValueError("manifest identity mismatch")
    if result.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("manifest schema identity mismatch")
    if result.get("kind") != "adaptive_chi_phase1_manifest":
        raise ValueError("manifest kind identity mismatch")
    expected = build_manifest(_seed_base(result.get("seed_base")))
    if result != expected:
        raise ValueError("manifest identity does not match the phase-1 specification")
    for stage in result["stages"].values():
        if any(not _safe_output(task.get("output")) for task in stage["tasks"]):
            raise ValueError("manifest contains an unsafe output path")
    return result


def _manifest_path(root: os.PathLike[str] | str) -> Path:
    return Path(root).expanduser().resolve() / MANIFEST_NAME


def create_manifest(
    root: os.PathLike[str] | str, *, seed_base: int
) -> dict[str, Any]:
    """Create the manifest once; an existing valid manifest is never rewritten."""
    seed = _seed_base(seed_base)
    path = _manifest_path(root)
    if path.exists():
        existing = load_manifest(path.parent)
        if existing["seed_base"] != seed:
            raise ValueError("existing manifest has a different seed_base")
        return existing
    document = build_manifest(seed)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    serialized = json.dumps(document, indent=2, sort_keys=True, allow_nan=False) + "\n"
    try:
        with temporary.open("x", encoding="utf-8") as stream:
            stream.write(serialized)
            stream.flush()
            os.fsync(stream.fileno())
        try:
            os.link(temporary, path)
        except FileExistsError:
            existing = load_manifest(path.parent)
            if existing["seed_base"] != seed:
                raise ValueError("existing manifest has a different seed_base")
            return existing
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
    return document


def load_manifest(root: os.PathLike[str] | str) -> dict[str, Any]:
    path = _manifest_path(root)
    if not path.is_file():
        raise FileNotFoundError(f"manifest not found: {path}")
    with path.open("r", encoding="utf-8") as stream:
        document = json.load(stream)
    if not isinstance(document, dict):
        raise ValueError("manifest identity requires a JSON object")
    return validate_manifest(document)


def tasks_for_stage(document: Mapping[str, Any], stage: str) -> list[dict[str, Any]]:
    validated = validate_manifest(document)
    if stage not in validated["stages"]:
        raise ValueError(f"unknown stage: {stage}")
    return list(validated["stages"][stage]["tasks"])


def task_for_id(
    document: Mapping[str, Any], stage: str, task_id: int
) -> dict[str, Any]:
    tasks = tasks_for_stage(document, stage)
    if isinstance(task_id, bool) or not isinstance(task_id, int) or not 0 <= task_id < len(tasks):
        raise ValueError(f"task_id {task_id!r} is outside stage {stage}")
    task = tasks[task_id]
    if task["task_id"] != task_id:
        raise ValueError("manifest task identity is not contiguous")
    return task
