#!/usr/bin/env python3
"""Prepare immutable adaptive-chi corpora and execution manifests on HPC4."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
from typing import Any, Iterable, Mapping

try:
    import tomllib
except ImportError:  # Dry-run remains available on the HPC login-node Python 3.9.
    tomllib = None

if __package__:
    from .. import adaptive_runner, core, execution_manifest, manifest, task_runner
else:
    HERE = Path(__file__).resolve().parent
    sys.path.insert(0, str(HERE.parent))
    import adaptive_runner, core, execution_manifest, manifest, task_runner


SCHEMA_VERSION = 1
PROFILE_NAME = "environment-profile.json"
PROVENANCE_NAME = "corpus-provenance.json"
STAGE_ROOTS = {
    ("l3_exact", False): "l3-primary",
    ("l3_exact", True): "l3-blas1-replay",
    ("l5_feasibility", False): "l5-primary",
}


def _write_once_json(path: Path, document: Mapping[str, Any], label: str) -> Path:
    serialized = json.dumps(document, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if path.exists():
        if json.loads(path.read_text(encoding="utf-8")) != document:
            raise ValueError(f"existing {label} has a different identity")
        return path
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
            if json.loads(path.read_text(encoding="utf-8")) != document:
                raise ValueError(f"concurrent {label} has a different identity")
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
    return path


def stage_root(study_root: os.PathLike[str] | str, stage: str, *, replay: bool = False) -> Path:
    key = (stage, replay)
    if key not in STAGE_ROOTS:
        raise ValueError("unsupported adaptive HPC stage root")
    return Path(study_root).expanduser().resolve() / STAGE_ROOTS[key]


def plan_corpus_task(
    study_root: os.PathLike[str] | str,
    sampling_manifest: Mapping[str, Any],
    *,
    stage: str,
    task_id: int,
) -> dict[str, Any]:
    sampling = manifest.validate_manifest(sampling_manifest)
    task = manifest.task_for_id(sampling, stage, task_id)
    shard = stage_root(study_root, stage) / "corpus-shards" / f"task-{task_id:03d}.json"
    return {
        "status": "dry_run",
        "sampling_manifest_id": sampling["manifest_id"],
        "stage": stage,
        "task_id": task_id,
        "L": task["L"],
        "p_key": task["p_key"],
        "record_count": task["shots"],
        "history_seeds": list(task["history_seeds"]),
        "seed_mode": "fresh_xoshiro_per_history_seed",
        "shard_path": str(shard),
    }


def _runtime(value: Mapping[str, Any], blas_threads: int) -> dict[str, Any]:
    result = dict(value)
    required = {"julia_version", "python_version", "blas_vendor", "cpu_model",
                "julia_threads", "blas_threads"}
    if set(result) != required or result.get("julia_threads") != 1 or result.get("blas_threads") != blas_threads:
        raise ValueError("runtime evidence does not match the frozen thread profile")
    return result


def build_environment_profile(
    *,
    source_identity: Mapping[str, Any],
    primary_runtime: Mapping[str, Any],
    replay_runtime: Mapping[str, Any],
    extended_runtime: Mapping[str, Any],
) -> dict[str, Any]:
    source = execution_manifest._validate_source_identity(source_identity)
    primary = _runtime(primary_runtime, 4)
    replay = _runtime(replay_runtime, 1)
    if not isinstance(extended_runtime, Mapping) or not extended_runtime:
        raise ValueError("extended runtime evidence must be nonempty")
    document = {
        "schema_version": SCHEMA_VERSION,
        "kind": "adaptive_chi_hpc_environment_profile",
        "source_identity": source,
        "source_bundle_id": core.identity_hash(source, "source_bundle_id"),
        "primary_runtime": primary,
        "primary_runtime_profile_id": core.identity_hash(primary, "runtime_profile_id"),
        "replay_runtime": replay,
        "replay_runtime_profile_id": core.identity_hash(replay, "runtime_profile_id"),
        "extended_runtime": dict(sorted(extended_runtime.items())),
    }
    document["environment_profile_id"] = core.identity_hash(
        document, "environment_profile_id"
    )
    return document


def validate_environment_profile(document: Mapping[str, Any]) -> dict[str, Any]:
    result = dict(document)
    identity = result.get("environment_profile_id")
    if identity != core.identity_hash(result, "environment_profile_id"):
        raise ValueError("environment profile identity mismatch")
    rebuilt = build_environment_profile(
        source_identity=result.get("source_identity"),
        primary_runtime=result.get("primary_runtime"),
        replay_runtime=result.get("replay_runtime"),
        extended_runtime=result.get("extended_runtime"),
    )
    if result != rebuilt:
        raise ValueError("environment profile derived identity mismatch")
    return result


def load_environment_profile(study_root: os.PathLike[str] | str) -> dict[str, Any]:
    path = Path(study_root).expanduser().resolve() / PROFILE_NAME
    if not path.is_file():
        raise FileNotFoundError(f"environment profile not found: {path}")
    return validate_environment_profile(json.loads(path.read_text(encoding="utf-8")))


def verify_live_environment_profile(
    study_root: os.PathLike[str] | str,
    *,
    repository_root: os.PathLike[str] | str,
    julia: str = "julia",
) -> dict[str, Any]:
    """Reauthenticate the Python/source environment frozen during setup."""
    root = Path(study_root).expanduser().resolve()
    repository = Path(repository_root).expanduser().resolve()
    profile = load_environment_profile(root)
    execution_manifest.verify_source_identity(repository, profile["source_identity"])
    live_python_version = platform.python_version()
    if any(
        profile[name].get("python_version") != live_python_version
        for name in ("primary_runtime", "replay_runtime")
    ):
        raise ValueError("live Python version differs from the frozen runtime profiles")
    freeze = subprocess.run(
        [sys.executable, "-m", "pip", "freeze", "--all"],
        check=False, capture_output=True, text=True, timeout=120,
    )
    if freeze.returncode != 0:
        raise RuntimeError(f"pip freeze failed: {freeze.stderr[-2000:]}")
    julia_path = shutil.which(julia) if not Path(julia).is_absolute() else julia
    if julia_path is None:
        raise RuntimeError("frozen Julia executable is unavailable")
    current = {
        "platform_machine": platform.machine(),
        "platform_system": platform.system(),
        "python_executable": str(Path(sys.executable).resolve()),
        "python_implementation": platform.python_implementation(),
        "pip_freeze_sha256": hashlib.sha256(freeze.stdout.encode()).hexdigest(),
        "repository_realpath": str(repository),
        "julia_executable": str(Path(julia_path).expanduser().resolve()),
        "python_packages": _package_versions(("numpy", "scipy", "pymatching", "ldpc")),
    }
    if current != profile["extended_runtime"]:
        raise ValueError("live Python/platform environment differs from the frozen profile")
    return profile


def _probe_compact_runtime(repository_root: Path, julia: str, blas_threads: int) -> dict[str, Any]:
    if tomllib is None:
        raise RuntimeError("runtime probing requires the validated Python 3.12 environment")
    adapter = repository_root / "examples/tnml_repaired_seed/adaptive_chi_pilot/julia_adapter.jl"
    with tempfile.TemporaryDirectory(prefix="tnml-adaptive-profile-") as temporary:
        output = Path(temporary) / "runtime.toml"
        environment = os.environ.copy()
        environment.update(
            TNML_BLAS_THREADS=str(blas_threads),
            OPENBLAS_NUM_THREADS=str(blas_threads),
            MKL_NUM_THREADS=str(blas_threads),
        )
        completed = subprocess.run(
            [julia, f"--project={repository_root}", "--startup-file=no", "--threads=1",
             str(adapter), "probe", str(output)],
            check=False, capture_output=True, text=True, timeout=300, env=environment,
        )
        if completed.returncode != 0 or not output.is_file():
            raise RuntimeError(f"Julia runtime probe failed: {completed.stderr[-2000:]}")
        evidence = tomllib.loads(output.read_text(encoding="utf-8"))["runtime_identity"]
    evidence["python_version"] = platform.python_version()
    return _runtime(evidence, blas_threads)


def _package_versions(names: Iterable[str]) -> dict[str, str]:
    result = {}
    for name in names:
        try:
            result[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError as error:
            raise RuntimeError(f"required Python package is unavailable: {name}") from error
    return result


def capture_environment_profile(
    repository_root: os.PathLike[str] | str,
    *,
    julia: str = "julia",
) -> dict[str, Any]:
    repository = Path(repository_root).expanduser().resolve()
    source = execution_manifest.capture_source_identity(repository)
    freeze = subprocess.run(
        [sys.executable, "-m", "pip", "freeze", "--all"],
        check=False, capture_output=True, text=True, timeout=120,
    )
    if freeze.returncode != 0:
        raise RuntimeError(f"pip freeze failed: {freeze.stderr[-2000:]}")
    extended = {
        "platform_machine": platform.machine(),
        "platform_system": platform.system(),
        "python_executable": str(Path(sys.executable).resolve()),
        "python_implementation": platform.python_implementation(),
        "pip_freeze_sha256": hashlib.sha256(freeze.stdout.encode()).hexdigest(),
        "repository_realpath": str(repository),
        "julia_executable": str(Path(julia).expanduser().resolve()),
        "python_packages": _package_versions(("numpy", "scipy", "pymatching", "ldpc")),
    }
    return build_environment_profile(
        source_identity=source,
        primary_runtime=_probe_compact_runtime(repository, julia, 4),
        replay_runtime=_probe_compact_runtime(repository, julia, 1),
        extended_runtime=extended,
    )


def build_corpus_shard(
    sampling_manifest: Mapping[str, Any],
    stage: str,
    task_id: int,
    rows: Iterable[Mapping[str, Any]],
    *,
    source_identity: Mapping[str, Any],
    environment_profile_id: str,
) -> dict[str, Any]:
    sampling = manifest.validate_manifest(sampling_manifest)
    task = manifest.task_for_id(sampling, stage, task_id)
    source = execution_manifest._validate_source_identity(source_identity)
    if not isinstance(environment_profile_id, str) or len(environment_profile_id) != 64:
        raise ValueError("environment_profile_id is malformed")
    values = list(rows)
    if len(values) != task["shots"]:
        raise ValueError("corpus shard record count differs from its planned task")
    records = []
    for shot_index, (row, history_seed) in enumerate(zip(values, task["history_seeds"]), 1):
        if row.get("history_seed") != history_seed:
            raise ValueError("corpus shard history_seed differs from its plan")
        record = execution_manifest._record_payload(row.get("record"), shot_index=shot_index)
        records.append({"shot_index": shot_index, "history_seed": history_seed, "record": record})
    document = {
        "schema_version": SCHEMA_VERSION,
        "kind": "adaptive_chi_corpus_shard",
        "sampling_manifest_id": sampling["manifest_id"],
        "stage": stage,
        "task_id": task_id,
        "L": task["L"],
        "p_key": task["p_key"],
        "seed_mode": "fresh_xoshiro_per_history_seed",
        "rng_algorithm": "Random.Xoshiro",
        "julia_record_identity": "TNMLRepairedRunner.sorted_toml_sha256_v1",
        "source_bundle_id": core.identity_hash(source, "source_bundle_id"),
        "environment_profile_id": environment_profile_id,
        "records": records,
    }
    document["shard_id"] = core.identity_hash(document, "shard_id")
    return document


def _validate_shard(document: Mapping[str, Any]) -> dict[str, Any]:
    result = dict(document)
    if result.get("shard_id") != core.identity_hash(result, "shard_id"):
        raise ValueError("corpus shard identity mismatch")
    return result


def collect_corpus_documents(
    sampling_manifest: Mapping[str, Any],
    stage: str,
    shards: Iterable[Mapping[str, Any]],
    *,
    source_identity: Mapping[str, Any],
    environment_profile_id: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    sampling = manifest.validate_manifest(sampling_manifest)
    source = execution_manifest._validate_source_identity(source_identity)
    source_bundle_id = core.identity_hash(source, "source_bundle_id")
    indexed = {}
    for raw in shards:
        shard = _validate_shard(raw)
        task_id = shard.get("task_id")
        if task_id in indexed:
            raise ValueError("duplicate corpus shard task identity")
        indexed[task_id] = shard
    expected = manifest.tasks_for_stage(sampling, stage)
    if set(indexed) != {task["task_id"] for task in expected}:
        raise ValueError("corpus shards are missing or contain an unexpected task")
    rows = []
    shard_entries = []
    for task in expected:
        shard = indexed[task["task_id"]]
        expected_fields = {
            "sampling_manifest_id": sampling["manifest_id"], "stage": stage,
            "task_id": task["task_id"], "L": task["L"], "p_key": task["p_key"],
            "source_bundle_id": source_bundle_id,
            "environment_profile_id": environment_profile_id,
            "seed_mode": "fresh_xoshiro_per_history_seed",
            "rng_algorithm": "Random.Xoshiro",
        }
        if any(shard.get(key) != value for key, value in expected_fields.items()):
            raise ValueError("corpus shard provenance differs from the frozen plan")
        canonical = build_corpus_shard(
            sampling, stage, task["task_id"], shard.get("records", []),
            source_identity=source, environment_profile_id=environment_profile_id,
        )
        if canonical != shard:
            raise ValueError("corpus shard derived identity mismatch")
        for item in shard["records"]:
            rows.append({"task_id": task["task_id"], "shot_index": item["shot_index"],
                         "history_seed": item["history_seed"], "record": item["record"]})
        shard_entries.append({"task_id": task["task_id"], "shard_id": shard["shard_id"]})
    corpus = execution_manifest.build_frozen_corpus(sampling, stage, rows)
    provenance = {
        "schema_version": SCHEMA_VERSION,
        "kind": "adaptive_chi_corpus_provenance",
        "sampling_manifest_id": sampling["manifest_id"],
        "stage": stage,
        "corpus_id": corpus["corpus_id"],
        "source_bundle_id": source_bundle_id,
        "environment_profile_id": environment_profile_id,
        "rng_algorithm": "Random.Xoshiro",
        "shards": shard_entries,
    }
    provenance["corpus_provenance_id"] = core.identity_hash(
        provenance, "corpus_provenance_id"
    )
    return corpus, provenance


def load_corpus_provenance(
    stage_directory: os.PathLike[str] | str,
    sampling_manifest: Mapping[str, Any],
    *,
    profile: Mapping[str, Any],
    execution: Mapping[str, Any],
    shard_directory: os.PathLike[str] | str | None = None,
) -> dict[str, Any]:
    """Rebuild and authenticate a frozen corpus from every retained shard."""
    stage_root_path = Path(stage_directory).expanduser().resolve()
    shard_root = Path(shard_directory or stage_root_path).expanduser().resolve()
    sampling = manifest.validate_manifest(sampling_manifest)
    validated_profile = validate_environment_profile(profile)
    corpus = execution_manifest.load_frozen_corpus(stage_root_path, sampling)
    provenance_path = stage_root_path / PROVENANCE_NAME
    if not provenance_path.is_file():
        raise FileNotFoundError(f"corpus provenance not found: {provenance_path}")
    observed = json.loads(provenance_path.read_text(encoding="utf-8"))
    stage = corpus["stage"]
    shards = []
    for task in manifest.tasks_for_stage(sampling, stage):
        path = shard_root / "corpus-shards" / f"task-{task['task_id']:03d}.json"
        if not path.is_file():
            raise FileNotFoundError(f"corpus shard not found: {path}")
        shards.append(json.loads(path.read_text(encoding="utf-8")))
    rebuilt_corpus, rebuilt_provenance = collect_corpus_documents(
        sampling,
        stage,
        shards,
        source_identity=validated_profile["source_identity"],
        environment_profile_id=validated_profile["environment_profile_id"],
    )
    if observed != rebuilt_provenance:
        raise ValueError("corpus provenance differs from its retained shards")
    if corpus != rebuilt_corpus:
        raise ValueError("frozen corpus differs from its retained shards")
    hpc = execution.get("hpc_provenance")
    if not isinstance(hpc, Mapping):
        raise ValueError("execution manifest lacks frozen HPC provenance")
    expected = {
        "corpus_id": corpus["corpus_id"],
        "source_bundle_id": validated_profile["source_bundle_id"],
        "environment_profile_id": validated_profile["environment_profile_id"],
        "corpus_provenance_id": rebuilt_provenance["corpus_provenance_id"],
    }
    replay = execution.get("execution_variant") == "blas1_replay"
    expected_runtime = validated_profile[
        "replay_runtime" if replay else "primary_runtime"
    ]
    expected_runtime_id = validated_profile[
        "replay_runtime_profile_id" if replay else "primary_runtime_profile_id"
    ]
    if (
        execution.get("corpus_id") != expected["corpus_id"]
        or execution.get("source_bundle_id") != expected["source_bundle_id"]
        or hpc.get("environment_profile_id") != expected["environment_profile_id"]
        or hpc.get("corpus_provenance_id") != expected["corpus_provenance_id"]
        or execution.get("runtime_identity") != expected_runtime
        or execution.get("runtime_profile_id") != expected_runtime_id
    ):
        raise ValueError("execution manifest differs from authenticated corpus provenance")
    return rebuilt_provenance


def _rows_from_julia(
    document: Mapping[str, Any],
    plan: Mapping[str, Any],
    expected_runtime: Mapping[str, Any],
) -> list[dict[str, Any]]:
    expected_fields = {
        "schema_version", "kind", "L", "p", "p_key", "task_id",
        "seed_mode", "rng_algorithm", "runtime_identity", "history_seeds", "records",
    }
    if set(document) != expected_fields:
        raise ValueError("Julia corpus envelope fields differ from the frozen schema")
    if document.get("schema_version") != 1 or document.get("kind") != "adaptive_chi_julia_corpus_task":
        raise ValueError("Julia corpus envelope schema/kind is invalid")
    if document.get("task_id") != plan["task_id"] or document.get("L") != plan["L"]:
        raise ValueError("Julia corpus output differs from the requested task")
    if (
        document.get("p_key") != plan["p_key"]
        or document.get("p") != float(plan["p_key"])
        or document.get("seed_mode") != plan["seed_mode"]
    ):
        raise ValueError("Julia corpus output differs from the requested probability/seed mode")
    if document.get("rng_algorithm") != "Random.Xoshiro":
        raise ValueError("Julia corpus output reports an unexpected RNG algorithm")
    try:
        adaptive_runner.JuliaExecutor.verify_julia_runtime(
            document.get("runtime_identity"), expected_runtime
        )
    except adaptive_runner.IdentityMismatchError as error:
        raise ValueError(str(error)) from error
    if document.get("history_seeds") != plan["history_seeds"]:
        raise ValueError("Julia corpus history seed list/order differs from the frozen plan")
    records = document.get("records")
    if not isinstance(records, list) or len(records) != plan["record_count"]:
        raise ValueError("Julia corpus output has an incorrect record count")
    return [
        {"history_seed": seed, "record": record}
        for seed, record in zip(plan["history_seeds"], records)
    ]


def run_corpus_task(
    study_root: os.PathLike[str] | str,
    sampling_manifest: Mapping[str, Any],
    *,
    stage: str,
    task_id: int,
    repository_root: os.PathLike[str] | str,
    julia: str = "julia",
    dry_run: bool = False,
) -> dict[str, Any]:
    plan = plan_corpus_task(study_root, sampling_manifest, stage=stage, task_id=task_id)
    if dry_run:
        return plan
    if stage == "l5_feasibility":
        if __package__:
            from . import gate
        else:
            import gate
        gate.require_l3_gate(study_root)
    repository = Path(repository_root).expanduser().resolve()
    if tomllib is None:
        raise RuntimeError("corpus generation requires the validated Python 3.12 environment")
    profile = verify_live_environment_profile(
        study_root, repository_root=repository, julia=julia
    )
    current_source = execution_manifest.capture_source_identity(repository)
    if current_source != profile["source_identity"]:
        raise ValueError("corpus generator source differs from the frozen environment profile")
    generator = repository / "examples/tnml_repaired_seed/adaptive_chi_pilot/prepare_corpus.jl"
    with tempfile.TemporaryDirectory(prefix="tnml-corpus-task-") as temporary:
        output = Path(temporary) / "corpus.toml"
        command = [
            julia, f"--project={repository}", "--startup-file=no", "--threads=1",
            str(generator), "generate", str(plan["L"]), plan["p_key"],
            str(plan["task_id"]), str(output), *map(str, plan["history_seeds"]),
        ]
        environment = os.environ.copy()
        environment.update(TNML_BLAS_THREADS="4", OPENBLAS_NUM_THREADS="4", MKL_NUM_THREADS="4")
        completed = subprocess.run(
            command, check=False, capture_output=True, text=True, timeout=3600, env=environment
        )
        if completed.returncode != 0 or not output.is_file():
            raise RuntimeError(f"corpus generation failed: {completed.stderr[-2000:]}")
        document = tomllib.loads(output.read_text(encoding="utf-8"))
    rows = _rows_from_julia(document, plan, profile["primary_runtime"])
    shard = build_corpus_shard(
        sampling_manifest, stage, task_id, rows,
        source_identity=current_source,
        environment_profile_id=profile["environment_profile_id"],
    )
    _write_once_json(Path(plan["shard_path"]), shard, "corpus shard")
    return {**plan, "status": "complete", "shard_id": shard["shard_id"]}


def freeze_stage(
    study_root: os.PathLike[str] | str,
    sampling_manifest: Mapping[str, Any],
    *,
    stage: str,
    repository_root: os.PathLike[str] | str,
) -> dict[str, Any]:
    root = Path(study_root).expanduser().resolve()
    repository = Path(repository_root).expanduser().resolve()
    profile = verify_live_environment_profile(
        root, repository_root=repository, julia="julia"
    )
    source = execution_manifest.capture_source_identity(repository)
    if source != profile["source_identity"]:
        raise ValueError("freeze source differs from the environment profile")
    shard_dir = stage_root(root, stage) / "corpus-shards"
    shards = [json.loads((shard_dir / f"task-{task['task_id']:03d}.json").read_text())
              for task in manifest.tasks_for_stage(sampling_manifest, stage)]
    corpus, provenance = collect_corpus_documents(
        sampling_manifest, stage, shards, source_identity=source,
        environment_profile_id=profile["environment_profile_id"],
    )
    primary_root = stage_root(root, stage)
    execution_manifest.write_frozen_corpus(primary_root, corpus)
    _write_once_json(primary_root / PROVENANCE_NAME, provenance, "corpus provenance")
    gate_id = None
    if stage == "l5_feasibility":
        if __package__:
            from . import gate
        else:
            import gate
        gate_id = gate.require_l3_gate(root)["gate_id"]
    hpc_provenance = {
        "environment_profile_id": profile["environment_profile_id"],
        "corpus_provenance_id": provenance["corpus_provenance_id"],
        "l3_gate_id": gate_id,
    }
    primary = execution_manifest.create_execution_manifest(
        primary_root, sampling_manifest, source_identity=source,
        runtime_identity=profile["primary_runtime"], hpc_provenance=hpc_provenance,
    )
    result = {"primary_root": str(primary_root),
              "primary_execution_manifest_id": primary["execution_manifest_id"],
              "corpus_id": corpus["corpus_id"],
              "corpus_provenance_id": provenance["corpus_provenance_id"]}
    if stage == "l3_exact":
        replay_root = stage_root(root, stage, replay=True)
        execution_manifest.write_frozen_corpus(replay_root, corpus)
        _write_once_json(replay_root / PROVENANCE_NAME, provenance, "replay corpus provenance")
        replay = execution_manifest.create_execution_manifest(
            replay_root, sampling_manifest, source_identity=source,
            runtime_identity=profile["replay_runtime"], execution_variant="blas1_replay",
            hpc_provenance=hpc_provenance,
        )
        result.update(replay_root=str(replay_root),
                      replay_execution_manifest_id=replay["execution_manifest_id"])
    return result


def slurm_array_spec(sampling_manifest: Mapping[str, Any], stage: str, *, batch: int | None = None) -> str:
    sampling = manifest.validate_manifest(sampling_manifest)
    count = len(manifest.tasks_for_stage(sampling, stage))
    if stage == "l3_exact" and batch is None and count == 30:
        return "0-29%20"
    if stage == "l5_feasibility" and batch in (0, 1) and count == 60:
        return "0-29%20" if batch == 0 else "30-59%20"
    raise ValueError("unsupported Slurm stage/batch schedule")


def replay_array_spec() -> str:
    return "0-29%20"


def require_hpc_execution(
    study_root: os.PathLike[str] | str,
    sampling: Mapping[str, Any],
    *,
    target: str,
) -> dict[str, Any]:
    """Authenticate retained corpus evidence and any required L3 gate."""
    root = Path(study_root).expanduser().resolve()
    if target == "l3-primary":
        stage_dir = stage_root(root, "l3_exact")
        shard_dir = stage_dir
    elif target == "l3-replay":
        stage_dir = stage_root(root, "l3_exact", replay=True)
        shard_dir = stage_root(root, "l3_exact")
    elif target == "l5-primary":
        stage_dir = stage_root(root, "l5_feasibility")
        shard_dir = stage_dir
    else:
        raise ValueError("unknown HPC execution target")
    execution = execution_manifest.load_execution_manifest(stage_dir, sampling)
    if target == "l5-primary":
        task_runner.verify_l5_gate_binding(stage_dir, execution, require_hpc=True)
    profile = load_environment_profile(root)
    provenance = load_corpus_provenance(
        stage_dir,
        sampling,
        profile=profile,
        execution=execution,
        shard_directory=shard_dir,
    )
    return {
        "status": "authenticated",
        "target": target,
        "execution_manifest_id": execution["execution_manifest_id"],
        "corpus_provenance_id": provenance["corpus_provenance_id"],
        "l3_gate_id": execution["hpc_provenance"].get("l3_gate_id"),
    }


def require_l5_execution(
    study_root: os.PathLike[str] | str,
    sampling: Mapping[str, Any],
) -> dict[str, Any]:
    return require_hpc_execution(study_root, sampling, target="l5-primary")


def smoke_plan(
    study_root: os.PathLike[str] | str,
    sampling: Mapping[str, Any],
    *,
    target: str,
    probe: bool = False,
    repository_root: os.PathLike[str] | str | None = None,
    julia: str = "julia",
    executor: Any | None = None,
) -> dict[str, Any]:
    root = Path(study_root).expanduser().resolve()
    if target == "l3-primary":
        stage_dir, count = stage_root(root, "l3_exact"), 30
        plans = [task_runner.plan_manifest_task(stage_dir, sampling, task_id=index)
                 for index in range(count)]
    elif target == "l3-replay":
        stage_dir, count = stage_root(root, "l3_exact", replay=True), 30
        plans = [task_runner.plan_replay_record(stage_dir, sampling, record_index=index)
                 for index in range(count)]
    elif target == "l5-primary":
        if __package__:
            from . import gate
        else:
            import gate
        gate.require_l3_gate(root)
        stage_dir, count = stage_root(root, "l5_feasibility"), 60
        plans = [task_runner.plan_manifest_task(stage_dir, sampling, task_id=index)
                 for index in range(count)]
    else:
        raise ValueError("unknown smoke target")
    repository = Path(repository_root or Path(__file__).resolve().parents[4]).resolve()
    verify_live_environment_profile(
        root, repository_root=repository, julia=julia
    )
    require_hpc_execution(root, sampling, target=target)
    execution = execution_manifest.load_execution_manifest(stage_dir, sampling)
    result = {
        "status": "dry_run" if not probe else "complete",
        "target": target,
        "task_count": len(plans),
        "record_count": sum(item["requested_records"] for item in plans),
        "source_identity_verified": True,
        "runtime_probe": "not_requested",
    }
    if probe:
        first_task = next(task for task in execution["tasks"] if task["records"])
        first_record = first_task["records"][0]
        owned = executor or task_runner.PersistentJuliaExecutor(
            repository_root=repository, julia=julia
        )
        executor_identity = adaptive_runner._executor_checkpoint_identity(
            owned, execution
        )
        identity = adaptive_runner._execution_identity(
            execution, first_task, first_record, executor_identity=executor_identity
        )
        try:
            owned.preflight(identity)
            result["runtime_probe"] = "authenticated_worker_ready_no_contraction"
        finally:
            close = getattr(owned, "close", None)
            if callable(close):
                close()
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    setup = commands.add_parser("setup")
    setup.add_argument("--study-root", required=True)
    setup.add_argument("--seed-base", required=True, type=int)
    setup.add_argument("--repository-root", required=True)
    setup.add_argument("--julia", default="julia")
    setup.add_argument("--dry-run", action="store_true")
    corpus = commands.add_parser("corpus-task")
    corpus.add_argument("--study-root", required=True)
    corpus.add_argument("--stage", required=True, choices=("l3_exact", "l5_feasibility"))
    corpus.add_argument("--task-id", required=True, type=int)
    corpus.add_argument("--repository-root", required=True)
    corpus.add_argument("--julia", default="julia")
    corpus.add_argument("--dry-run", action="store_true")
    freeze = commands.add_parser("freeze")
    freeze.add_argument("--study-root", required=True)
    freeze.add_argument("--stage", required=True, choices=("l3_exact", "l5_feasibility"))
    freeze.add_argument("--repository-root", required=True)
    smoke = commands.add_parser("smoke")
    smoke.add_argument("--study-root", required=True)
    smoke.add_argument("--target", required=True, choices=("l3-primary", "l3-replay", "l5-primary"))
    smoke.add_argument("--repository-root", default=str(Path(__file__).resolve().parents[4]))
    smoke.add_argument("--julia", default="julia")
    smoke.add_argument("--probe", action="store_true")
    require_gate = commands.add_parser("require-l3-gate")
    require_gate.add_argument("--study-root", required=True)
    require_execution = commands.add_parser("require-execution")
    require_execution.add_argument("--study-root", required=True)
    require_execution.add_argument(
        "--target", required=True, choices=("l3-primary", "l3-replay", "l5-primary")
    )
    require_execution.add_argument(
        "--repository-root", default=str(Path(__file__).resolve().parents[4])
    )
    require_execution.add_argument("--julia", default="julia")
    args = parser.parse_args(argv)
    root = Path(args.study_root).expanduser().resolve()
    if args.command == "setup":
        if args.dry_run:
            result = {"status": "dry_run", "study_root": str(root),
                      "seed_base": args.seed_base, "blas_profiles": [4, 1]}
        else:
            sampling = manifest.create_manifest(root, seed_base=args.seed_base)
            profile = capture_environment_profile(args.repository_root, julia=args.julia)
            _write_once_json(root / PROFILE_NAME, profile, "environment profile")
            result = {"status": "complete", "manifest_id": sampling["manifest_id"],
                      "environment_profile_id": profile["environment_profile_id"]}
    elif args.command == "require-l3-gate":
        if __package__:
            from . import gate
        else:
            import gate
        result = gate.require_l3_gate(root)
    else:
        sampling = manifest.load_manifest(root)
        if args.command == "require-execution":
            verify_live_environment_profile(
                root, repository_root=args.repository_root, julia=args.julia
            )
            result = require_hpc_execution(root, sampling, target=args.target)
        elif args.command == "corpus-task":
            result = run_corpus_task(
                root, sampling, stage=args.stage, task_id=args.task_id,
                repository_root=args.repository_root, julia=args.julia,
                dry_run=args.dry_run,
            )
        elif args.command == "freeze":
            result = freeze_stage(
                root, sampling, stage=args.stage, repository_root=args.repository_root
            )
        else:
            result = smoke_plan(
                root, sampling, target=args.target, probe=args.probe,
                repository_root=args.repository_root, julia=args.julia,
            )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
