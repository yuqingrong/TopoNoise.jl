#!/usr/bin/env python3
"""Frozen corpus and phase-2 execution-manifest identities."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
from typing import Any, Iterable, Mapping

try:
    from . import core, manifest
except ImportError:
    import core
    import manifest


SCHEMA_VERSION = 1
CORPUS_NAME = "frozen-corpus.json"
EXECUTION_MANIFEST_NAME = "execution-manifest.json"
HEX64 = re.compile(r"^[0-9a-f]{64}$")
REQUIRED_SOURCE_FILES = frozenset(
    {
        "examples/tnml_repaired_seed/model.jl",
        "examples/tnml_repaired_seed/run.jl",
        "examples/tnml_repaired_seed/dense_frontier.jl",
        "examples/tnml_repaired_seed/adaptive_chi_pilot/core.py",
        "examples/tnml_repaired_seed/adaptive_chi_pilot/manifest.py",
        "examples/tnml_repaired_seed/adaptive_chi_pilot/execution_manifest.py",
        "examples/tnml_repaired_seed/adaptive_chi_pilot/adaptive_runner.py",
        "examples/tnml_repaired_seed/adaptive_chi_pilot/julia_adapter.jl",
    }
)


def _write_once_json(path: Path, document: Mapping[str, Any], identity_name: str) -> None:
    serialized = json.dumps(
        document, indent=2, sort_keys=True, allow_nan=False
    ) + "\n"
    if path.exists():
        existing = json.loads(path.read_text(encoding="utf-8"))
        if existing != document:
            raise ValueError(f"existing {identity_name} has a different identity")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("x", encoding="utf-8") as stream:
            stream.write(serialized)
            stream.flush()
            os.fsync(stream.fileno())
        try:
            os.link(temporary, path)
        except FileExistsError:
            existing = json.loads(path.read_text(encoding="utf-8"))
            if existing != document:
                raise ValueError(f"existing {identity_name} has a different identity")
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def _record_payload(record: Any, *, shot_index: int) -> dict[str, Any]:
    if not isinstance(record, Mapping):
        raise ValueError("each corpus record must be an object")
    result = dict(record)
    required = {
        "shot_index",
        "raw_fault_count",
        "raw_fault_indices",
        "seed_syndrome",
        "seed_recovery",
        "recovery_mask",
        "final_syndrome",
        "true_mask",
        "mwpm_failure",
        "raw_seed_mask",
        "transfer_mask",
        "record_id",
    }
    if set(result) != required:
        raise ValueError("corpus record fields do not match the frozen Julia record schema")
    if result["shot_index"] != shot_index:
        raise ValueError("corpus record shot_index differs from its frozen locator")
    record_id = result["record_id"]
    if not isinstance(record_id, str) or HEX64.fullmatch(record_id) is None:
        raise ValueError("record_id must be 64 lowercase hexadecimal characters")
    count = result["raw_fault_count"]
    if type(count) is not int or count <= 0:
        raise ValueError("raw_fault_count must be positive")
    indices = result["raw_fault_indices"]
    if (
        not isinstance(indices, list)
        or any(type(value) is not int or not 1 <= value <= count for value in indices)
        or indices != sorted(set(indices))
    ):
        raise ValueError("raw_fault_indices must be sorted unique Julia indices")
    for name in ("seed_syndrome", "seed_recovery", "final_syndrome"):
        values = result[name]
        if not isinstance(values, list) or any(type(value) is not bool for value in values):
            raise ValueError(f"{name} must be a Boolean array")
    for name in ("recovery_mask", "true_mask", "raw_seed_mask", "transfer_mask"):
        if type(result[name]) is not int or not 0 <= result[name] < 4:
            raise ValueError(f"{name} must be an integer in 0:3")
    if type(result["mwpm_failure"]) is not bool:
        raise ValueError("mwpm_failure must be Boolean")
    return result


def _expected_locators(sampling: Mapping[str, Any], stage: str) -> list[dict[str, int]]:
    locators = []
    for task in manifest.tasks_for_stage(sampling, stage):
        for shot_index, history_seed in enumerate(task["history_seeds"], start=1):
            locators.append(
                {
                    "task_id": task["task_id"],
                    "shot_index": shot_index,
                    "history_seed": history_seed,
                }
            )
    return locators


def build_frozen_corpus(
    sampling_manifest: Mapping[str, Any], stage: str, rows: Iterable[Mapping[str, Any]]
) -> dict[str, Any]:
    sampling = manifest.validate_manifest(sampling_manifest)
    tasks = manifest.tasks_for_stage(sampling, stage)
    task_by_id = {task["task_id"]: task for task in tasks}
    row_values = list(rows)
    expected = _expected_locators(sampling, stage)
    if len(row_values) != len(expected):
        raise ValueError(
            f"frozen corpus has {len(row_values)} records; expected {len(expected)}"
        )
    indexed: dict[tuple[int, int], Mapping[str, Any]] = {}
    for row in row_values:
        if not isinstance(row, Mapping):
            raise ValueError("corpus rows must be objects")
        task_id, shot_index = row.get("task_id"), row.get("shot_index")
        if type(task_id) is not int or type(shot_index) is not int:
            raise ValueError("corpus locator must contain integer task_id and shot_index")
        locator = (task_id, shot_index)
        if locator in indexed:
            raise ValueError("frozen corpus contains a duplicate task/shot locator")
        indexed[locator] = row

    records = []
    instance_ids: set[str] = set()
    for locator in expected:
        key = (locator["task_id"], locator["shot_index"])
        if key not in indexed:
            raise ValueError(f"frozen corpus is missing task/shot locator {key}")
        row = indexed[key]
        if row.get("history_seed") != locator["history_seed"]:
            raise ValueError("frozen corpus history_seed differs from the sampling plan")
        task = task_by_id[locator["task_id"]]
        record = _record_payload(row.get("record"), shot_index=locator["shot_index"])
        julia_record_id = record["record_id"]
        content_sha256 = core.identity_hash(record, "record_content_sha256")
        instance_identity = {
            "version": "adaptive_chi_record_instance_v1",
            "sampling_manifest_id": sampling["manifest_id"],
            "stage": stage,
            "task_id": locator["task_id"],
            "shot_index": locator["shot_index"],
            "history_seed": locator["history_seed"],
            "julia_record_id": julia_record_id,
            "record_content_sha256": content_sha256,
        }
        record_id = core.identity_hash(instance_identity, "record_id")
        if record_id in instance_ids:
            raise ValueError("frozen corpus contains a duplicate record instance identity")
        instance_ids.add(record_id)
        records.append(
            {
                "task_id": locator["task_id"],
                "shot_index": locator["shot_index"],
                "history_seed": locator["history_seed"],
                "L": task["L"],
                "p_key": task["p_key"],
                "record_id": record_id,
                "julia_record_id": julia_record_id,
                "record_content_sha256": content_sha256,
                "record": record,
            }
        )
    document: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "kind": "adaptive_chi_frozen_corpus",
        "sampling_manifest_id": sampling["manifest_id"],
        "stage": stage,
        "record_count": len(records),
        "records": records,
    }
    document["corpus_id"] = core.identity_hash(document, "corpus_id")
    return document


def validate_frozen_corpus(
    document: Mapping[str, Any], sampling_manifest: Mapping[str, Any]
) -> dict[str, Any]:
    if not isinstance(document, Mapping):
        raise ValueError("corpus identity requires an object")
    result = dict(document)
    corpus_id = result.get("corpus_id")
    if not isinstance(corpus_id, str) or corpus_id != core.identity_hash(result, "corpus_id"):
        raise ValueError("frozen corpus identity mismatch")
    if result.get("schema_version") != SCHEMA_VERSION or result.get("kind") != "adaptive_chi_frozen_corpus":
        raise ValueError("frozen corpus schema identity mismatch")
    sampling = manifest.validate_manifest(sampling_manifest)
    if result.get("sampling_manifest_id") != sampling["manifest_id"]:
        raise ValueError("frozen corpus sampling-manifest identity mismatch")
    rows = [
        {
            "task_id": entry.get("task_id"),
            "shot_index": entry.get("shot_index"),
            "history_seed": entry.get("history_seed"),
            "record": entry.get("record"),
        }
        for entry in result.get("records", [])
    ]
    rebuilt = build_frozen_corpus(sampling, result.get("stage"), rows)
    if result != rebuilt:
        raise ValueError("frozen corpus derived identity mismatch")
    return result


def write_frozen_corpus(root: os.PathLike[str] | str, corpus: Mapping[str, Any]) -> Path:
    path = Path(root).expanduser().resolve() / CORPUS_NAME
    # Full plan validation happens when the execution manifest is built.
    if not isinstance(corpus.get("corpus_id"), str) or corpus["corpus_id"] != core.identity_hash(corpus, "corpus_id"):
        raise ValueError("frozen corpus identity mismatch")
    _write_once_json(path, corpus, "frozen corpus")
    return path


def load_frozen_corpus(
    root: os.PathLike[str] | str, sampling_manifest: Mapping[str, Any]
) -> dict[str, Any]:
    path = Path(root).expanduser().resolve() / CORPUS_NAME
    if not path.is_file():
        raise FileNotFoundError(f"frozen corpus not found: {path}")
    return validate_frozen_corpus(json.loads(path.read_text(encoding="utf-8")), sampling_manifest)


def _safe_relative_path(value: Any) -> bool:
    if not isinstance(value, str) or not value:
        return False
    path = PurePosixPath(value)
    return not path.is_absolute() and "." not in path.parts and ".." not in path.parts


def _validate_source_identity(value: Any) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError("source_identity must be an object")
    result = dict(value)
    if set(result) != {
        "project_sha256",
        "manifest_sha256",
        "package_source_tree_sha256",
        "files",
    }:
        raise ValueError("source_identity fields are incomplete")
    for name in (
        "project_sha256",
        "manifest_sha256",
        "package_source_tree_sha256",
    ):
        if not isinstance(result[name], str) or HEX64.fullmatch(result[name]) is None:
            raise ValueError(f"{name} must be a SHA-256 hex digest")
    files = result["files"]
    if not isinstance(files, Mapping) or not files:
        raise ValueError("source_identity.files must be a nonempty object")
    normalized = {}
    for path, digest in files.items():
        if not _safe_relative_path(path):
            raise ValueError("source identity contains an unsafe file path")
        if not isinstance(digest, str) or HEX64.fullmatch(digest) is None:
            raise ValueError("source file identity must be a SHA-256 hex digest")
        normalized[path] = digest
    if not REQUIRED_SOURCE_FILES <= normalized.keys():
        missing = sorted(REQUIRED_SOURCE_FILES - normalized.keys())
        raise ValueError(f"source_identity.files is missing required sources: {missing}")
    if not any(path.startswith("src/") and path.endswith(".jl") for path in normalized):
        raise ValueError("source_identity.files must bind the package Julia sources")
    result["files"] = dict(sorted(normalized.items()))
    return result


def _sha256_file(path: Path) -> str:
    if not path.is_file():
        raise ValueError(f"source identity file is missing: {path}")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def capture_source_identity(
    repository_root: os.PathLike[str] | str,
) -> dict[str, Any]:
    """Capture the exact source convention that the execution adapter verifies."""
    root = Path(repository_root).expanduser().resolve()
    package_sources = {
        path.relative_to(root).as_posix(): _sha256_file(path)
        for path in sorted((root / "src").rglob("*.jl"))
        if path.is_file()
    }
    source_paths = sorted(set(package_sources) | set(REQUIRED_SOURCE_FILES))
    result = {
        "project_sha256": _sha256_file(root / "Project.toml"),
        "manifest_sha256": _sha256_file(root / "Manifest.toml"),
        "package_source_tree_sha256": core.identity_hash(
            {"files": package_sources}, "package_source_tree_sha256"
        ),
        "files": {relative: _sha256_file(root / relative) for relative in source_paths},
    }
    return _validate_source_identity(result)


def verify_source_identity(
    repository_root: os.PathLike[str] | str, frozen_identity: Mapping[str, Any]
) -> None:
    """Refuse execution when any frozen scientific source has changed."""
    root = Path(repository_root).expanduser().resolve()
    expected = _validate_source_identity(frozen_identity)
    package_sources = {
        path.relative_to(root).as_posix(): _sha256_file(path)
        for path in sorted((root / "src").rglob("*.jl"))
        if path.is_file()
    }
    current = {
        "project_sha256": _sha256_file(root / "Project.toml"),
        "manifest_sha256": _sha256_file(root / "Manifest.toml"),
        "package_source_tree_sha256": core.identity_hash(
            {"files": package_sources}, "package_source_tree_sha256"
        ),
        "files": {
            relative: _sha256_file(root / relative)
            for relative in expected["files"]
        },
    }
    if current != expected:
        raise ValueError("current repository differs from the frozen source identity")


def _validate_runtime_identity(value: Any, sampling: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError("runtime_identity must be an object")
    result = dict(value)
    required = {
        "julia_version",
        "python_version",
        "blas_vendor",
        "cpu_model",
        "julia_threads",
        "blas_threads",
    }
    if set(result) != required:
        raise ValueError("runtime_identity fields are incomplete")
    for name in ("julia_version", "python_version", "blas_vendor", "cpu_model"):
        if not isinstance(result[name], str) or not result[name].strip():
            raise ValueError(f"runtime_identity.{name} must be nonempty")
    policy = sampling["runtime_policy"]
    if result["julia_threads"] != policy["julia_application_threads"]:
        raise ValueError("runtime Julia-thread identity differs from the sampling policy")
    if result["blas_threads"] != policy["blas_threads"]:
        raise ValueError("runtime BLAS-thread identity differs from the sampling policy")
    return result


def _execution_tasks(
    sampling: Mapping[str, Any], corpus: Mapping[str, Any], selected: Mapping[str, list[dict[str, Any]]]
) -> list[dict[str, Any]]:
    stage = corpus["stage"]
    entries_by_task: dict[int, list[dict[str, Any]]] = {}
    for entry in corpus["records"]:
        entries_by_task.setdefault(entry["task_id"], []).append(entry)
    selected_rank = {
        (p_key, item["record_id"]): item["rank"]
        for p_key, items in selected.items()
        for item in items
    }
    tasks = []
    for planned in manifest.tasks_for_stage(sampling, stage):
        records = []
        for entry in entries_by_task[planned["task_id"]]:
            rank = selected_rank.get((entry["p_key"], entry["record_id"]))
            records.append(
                {
                    "task_id": entry["task_id"],
                    "shot_index": entry["shot_index"],
                    "history_seed": entry["history_seed"],
                    "L": entry["L"],
                    "p_key": entry["p_key"],
                    "record_id": entry["record_id"],
                    "julia_record_id": entry["julia_record_id"],
                    "record_content_sha256": entry["record_content_sha256"],
                    "replay_rank": rank if stage == "l3_exact" else None,
                    "anchor_rank": rank if stage == "l5_feasibility" else None,
                    "forced_anchor": stage == "l5_feasibility" and rank is not None,
                }
            )
        tasks.append(
            {
                "task_id": planned["task_id"],
                "L": planned["L"],
                "p_key": planned["p_key"],
                "chi_ladder": list(planned["chi_ladder"]),
                "cutoff": planned["cutoff"],
                "exact_required": planned["exact"],
                "output": planned["output"],
                "records": records,
            }
        )
    return tasks


def build_execution_manifest(
    sampling_manifest: Mapping[str, Any],
    corpus_document: Mapping[str, Any],
    *,
    source_identity: Mapping[str, Any],
    runtime_identity: Mapping[str, Any],
    corpus_path: str = CORPUS_NAME,
) -> dict[str, Any]:
    sampling = manifest.validate_manifest(sampling_manifest)
    corpus = validate_frozen_corpus(corpus_document, sampling)
    sources = _validate_source_identity(source_identity)
    runtime = _validate_runtime_identity(runtime_identity, sampling)
    if corpus_path != CORPUS_NAME:
        raise ValueError(f"corpus_path must be the frozen path {CORPUS_NAME!r}")
    stage = corpus["stage"]
    L = sampling["stages"][stage]["L"]
    selection_count = 5 if stage == "l3_exact" else 4
    selected = {}
    for p_key in manifest.P_VALUES:
        record_ids = [
            entry["record_id"] for entry in corpus["records"] if entry["p_key"] == p_key
        ]
        selected[p_key] = core.rank_record_ids(
            L, float(p_key), record_ids, count=selection_count
        )
    tasks = _execution_tasks(sampling, corpus, selected)
    ladders = {tuple(task["chi_ladder"]) for task in tasks}
    if len(ladders) != 1:
        raise ValueError("execution stage must have one frozen chi ladder")
    chi_ladder = list(next(iter(ladders)))
    if stage == "l3_exact" and chi_ladder != [32, 64, 128, 240, 480, 960]:
        raise ValueError("L3 execution ladder differs from the preregistered policy")
    runtime_policy = sampling["runtime_policy"]
    wall_seconds = (
        dict(runtime_policy["l5_wall_seconds_by_chi"])
        if stage == "l5_feasibility"
        else {}
    )
    document: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "kind": "adaptive_chi_execution_manifest",
        "sampling_manifest_id": sampling["manifest_id"],
        "corpus_id": corpus["corpus_id"],
        "corpus_path": corpus_path,
        "stage": stage,
        "purpose": "exact_calibration" if stage == "l3_exact" else "feasibility_smoke",
        "source_identity": sources,
        "source_bundle_id": core.identity_hash(sources, "source_bundle_id"),
        "runtime_identity": runtime,
        "runtime_profile_id": core.identity_hash(runtime, "runtime_profile_id"),
        "chi_ladder": chi_ladder,
        "cutoff": sampling["cutoff"],
        "resource_policy": {
            "max_rss_bytes": runtime_policy["max_rss_bytes"],
            "wall_seconds_by_chi": wall_seconds,
            "exact_wall_seconds": None,
            "preflight_wall_seconds": 300,
            "enforcement": (
                "wall_timeout_internal_memory_cap_external_with_posthoc_rejection"
            ),
        },
        "selection": selected,
        "exact_audit": {
            "required": stage == "l3_exact",
            "ml_tie_tolerance": 1.0e-9,
            "posterior_linf_gate": 1.0e-5,
            "log_joint_gate": 1.0e-5,
            "gate_comparison": "inclusive",
            "false_accept_definition": "failure_of_any_class_posterior_or_log_joint_gate",
        },
        "clopper_pearson": {
            "method": "exact_beta_quantile",
            "confidence_level": 0.95,
            "sided": "one_sided_upper",
            "formula_when_events_less_than_trials": "BetaQuantile(0.95; events+1, trials-events)",
            "all_events": 1.0,
            "zero_trials": None,
            "denominator": "accepted_with_complete_exact_reference",
            "multiplicity_adjustment": "none",
            "stratification": "overall_and_each_probability",
        },
        "tasks": tasks,
    }
    document["execution_manifest_id"] = core.identity_hash(
        document, "execution_manifest_id"
    )
    return document


def validate_execution_manifest(
    document: Mapping[str, Any],
    sampling_manifest: Mapping[str, Any],
    corpus_document: Mapping[str, Any],
) -> dict[str, Any]:
    if not isinstance(document, Mapping):
        raise ValueError("execution manifest identity requires an object")
    result = dict(document)
    identity = result.get("execution_manifest_id")
    if not isinstance(identity, str) or identity != core.identity_hash(result, "execution_manifest_id"):
        raise ValueError("execution manifest identity mismatch")
    if result.get("schema_version") != SCHEMA_VERSION or result.get("kind") != "adaptive_chi_execution_manifest":
        raise ValueError("execution manifest schema identity mismatch")
    rebuilt = build_execution_manifest(
        sampling_manifest,
        corpus_document,
        source_identity=result.get("source_identity"),
        runtime_identity=result.get("runtime_identity"),
        corpus_path=result.get("corpus_path"),
    )
    if result != rebuilt:
        raise ValueError("execution manifest derived identity mismatch")
    return result


def create_execution_manifest(
    root: os.PathLike[str] | str,
    sampling_manifest: Mapping[str, Any],
    *,
    source_identity: Mapping[str, Any],
    runtime_identity: Mapping[str, Any],
) -> dict[str, Any]:
    root_path = Path(root).expanduser().resolve()
    corpus = load_frozen_corpus(root_path, sampling_manifest)
    document = build_execution_manifest(
        sampling_manifest,
        corpus,
        source_identity=source_identity,
        runtime_identity=runtime_identity,
    )
    _write_once_json(
        root_path / EXECUTION_MANIFEST_NAME, document, "execution manifest"
    )
    return document


def load_execution_manifest(
    root: os.PathLike[str] | str, sampling_manifest: Mapping[str, Any]
) -> dict[str, Any]:
    root_path = Path(root).expanduser().resolve()
    path = root_path / EXECUTION_MANIFEST_NAME
    if not path.is_file():
        raise FileNotFoundError(f"execution manifest not found: {path}")
    corpus = load_frozen_corpus(root_path, sampling_manifest)
    document = json.loads(path.read_text(encoding="utf-8"))
    return validate_execution_manifest(document, sampling_manifest, corpus)


def task_for_id(document: Mapping[str, Any], task_id: int) -> dict[str, Any]:
    tasks = document.get("tasks", [])
    if type(task_id) is not int or not 0 <= task_id < len(tasks):
        raise ValueError("task_id is outside the execution stage")
    task = tasks[task_id]
    if task.get("task_id") != task_id:
        raise ValueError("execution task identity is not contiguous")
    return task
