#!/usr/bin/env python3
"""Resumable per-record adaptive orchestration with immutable JSON checkpoints."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import tempfile
from typing import Any, Mapping, Protocol

try:
    from . import core, execution_manifest, manifest
except ImportError:
    import core
    import execution_manifest
    import manifest


SCHEMA_VERSION = 2
RSS_SCOPES = frozenset(
    {"process_lifetime_high_water", "worker_lifetime_high_water"}
)


class IdentityMismatchError(ValueError):
    """The current source or runtime cannot own a frozen study checkpoint."""


class Executor(Protocol):
    def run_tn(
        self,
        record: Mapping[str, Any],
        *,
        chi: int,
        cutoff: float,
        execution_identity: Mapping[str, Any],
    ) -> Mapping[str, Any]: ...

    def run_exact(
        self,
        record: Mapping[str, Any],
        *,
        execution_identity: Mapping[str, Any],
    ) -> Mapping[str, Any]: ...


def _read_json(path: Path) -> dict[str, Any]:
    document = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(document, dict):
        raise ValueError(f"checkpoint must contain a JSON object: {path}")
    return document


def _write_once(path: Path, document: Mapping[str, Any], label: str) -> None:
    serialized = json.dumps(
        document, indent=2, sort_keys=True, allow_nan=False
    ) + "\n"
    if path.exists():
        if _read_json(path) != document:
            raise ValueError(f"existing {label} differs; refusing to overwrite {path}")
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
            if _read_json(path) != document:
                raise ValueError(f"concurrent {label} differs; refusing {path}")
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def _find_record(
    execution: Mapping[str, Any],
    corpus: Mapping[str, Any],
    task_id: int,
    record_id: str,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    task = execution_manifest.task_for_id(execution, task_id)
    matches = [item for item in task["records"] if item["record_id"] == record_id]
    if len(matches) != 1:
        raise ValueError("record_id is not unique within the requested task")
    record_meta = matches[0]
    corpus_matches = [
        item
        for item in corpus["records"]
        if item["task_id"] == task_id
        and item["shot_index"] == record_meta["shot_index"]
        and item["record_id"] == record_id
    ]
    if len(corpus_matches) != 1:
        raise ValueError("execution record does not resolve uniquely in the frozen corpus")
    corpus_entry = corpus_matches[0]
    if corpus_entry["record_content_sha256"] != record_meta["record_content_sha256"]:
        raise ValueError("execution record content identity differs from the frozen corpus")
    if (
        corpus_entry["julia_record_id"] != record_meta["julia_record_id"]
        or corpus_entry["record"]["record_id"] != record_meta["julia_record_id"]
    ):
        raise ValueError("execution record differs from its authenticated Julia identity")
    return task, record_meta, corpus_entry


def _relative_record_output(task: Mapping[str, Any], record_meta: Mapping[str, Any]) -> str:
    return (
        f"outputs/{task['output']}/records/"
        f"shot-{record_meta['shot_index']:04d}-{record_meta['record_id'][:16]}"
    )


def plan_record(
    root: os.PathLike[str] | str,
    sampling_manifest: Mapping[str, Any],
    *,
    task_id: int,
    record_id: str,
) -> dict[str, Any]:
    root_path = Path(root).expanduser().resolve()
    execution = execution_manifest.load_execution_manifest(root_path, sampling_manifest)
    corpus = execution_manifest.load_frozen_corpus(root_path, sampling_manifest)
    task, record_meta, _ = _find_record(execution, corpus, task_id, record_id)
    relative = _relative_record_output(task, record_meta)
    output = (root_path / relative).resolve()
    if root_path not in output.parents:
        raise ValueError("record output resolves outside the study root")
    return {
        "status": "dry_run",
        "execution_manifest_id": execution["execution_manifest_id"],
        "corpus_id": execution["corpus_id"],
        "stage": execution["stage"],
        "task_id": task_id,
        "record_id": record_id,
        "shot_index": record_meta["shot_index"],
        "forced_anchor": record_meta["forced_anchor"],
        "exact_required": task["exact_required"],
        "chi_ladder": list(task["chi_ladder"]),
        "cutoff": task["cutoff"],
        "record_output": relative,
        "record_directory": str(output),
    }


def _execution_identity(
    execution: Mapping[str, Any],
    task: Mapping[str, Any],
    record_meta: Mapping[str, Any],
    executor_identity: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    result = {
        "execution_manifest_id": execution["execution_manifest_id"],
        "corpus_id": execution["corpus_id"],
        "source_bundle_id": execution["source_bundle_id"],
        "runtime_profile_id": execution["runtime_profile_id"],
        "source_identity": execution["source_identity"],
        "runtime_identity": execution["runtime_identity"],
        "resource_policy": execution["resource_policy"],
        "executor_policy": execution["executor_policy"],
        "L": task["L"],
        "p": float(task["p_key"]),
        "p_key": task["p_key"],
        "task_id": task["task_id"],
        "shot_index": record_meta["shot_index"],
        "record_id": record_meta["record_id"],
        "julia_record_id": record_meta["julia_record_id"],
        "cutoff": task["cutoff"],
    }
    if executor_identity is not None:
        result.update(executor_identity)
    return result


def _executor_checkpoint_identity(
    executor: Any, execution: Mapping[str, Any]
) -> dict[str, str]:
    mode = getattr(executor, "executor_mode", None)
    policy = execution.get("executor_policy")
    if not isinstance(policy, Mapping):
        raise IdentityMismatchError("execution manifest lacks executor policy")
    allowed = policy.get("allowed_executors")
    if not isinstance(allowed, list) or mode not in allowed:
        raise IdentityMismatchError("executor mode is not frozen in the manifest")
    scopes = policy.get("timing_scopes")
    if not isinstance(scopes, Mapping) or not isinstance(scopes.get(mode), Mapping):
        raise IdentityMismatchError("executor timing scope is not frozen")
    scope = scopes[mode]
    expected_keys = {"wall_time_scope", "decoder_cache_scope"}
    if set(scope) != expected_keys or any(
        not isinstance(scope[key], str) or not scope[key] for key in expected_keys
    ):
        raise IdentityMismatchError("executor timing scope is malformed")
    return {
        "executor_mode": mode,
        "wall_time_scope": scope["wall_time_scope"],
        "decoder_cache_scope": scope["decoder_cache_scope"],
    }


def _validate_checkpoint_executor_metrics(
    document: Mapping[str, Any], executor_identity: Mapping[str, str]
) -> None:
    if document.get("status") != "complete":
        return
    expected_rss_scope = (
        "worker_lifetime_high_water"
        if executor_identity["executor_mode"] == "persistent"
        else "process_lifetime_high_water"
    )
    if document.get("rss_scope") != expected_rss_scope:
        raise ValueError("complete checkpoint RSS scope does not match its executor")
    if executor_identity["executor_mode"] == "persistent":
        if type(document.get("decoder_cache_hit")) is not bool:
            raise ValueError("persistent checkpoint lacks decoder-cache evidence")
        for name in (
            "decoder_build_seconds", "decode_seconds", "worker_startup_seconds"
        ):
            if _safe_elapsed(document.get(name)) is None:
                raise ValueError("persistent checkpoint timing evidence is invalid")


def _exception_result(reason: str, error: BaseException) -> dict[str, Any]:
    return {
        "status": reason,
        "valid": False,
        "elapsed_seconds": 0.0,
        "peak_rss_bytes": 0,
        "error": f"{type(error).__name__}: {error}",
    }


def _call_executor(callable_result) -> dict[str, Any]:
    try:
        result = callable_result()
    except TimeoutError as error:
        return _exception_result("timeout", error)
    except MemoryError as error:
        return _exception_result("memory_limit", error)
    except IdentityMismatchError:
        raise
    except BaseException as error:
        if isinstance(error, (KeyboardInterrupt, SystemExit)):
            raise
        return _exception_result("exception", error)
    if not isinstance(result, Mapping):
        return _exception_result("exception", TypeError("executor returned a non-object"))
    return dict(result)


def _safe_elapsed(value: Any) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    result = float(value)
    return result if math.isfinite(result) and result >= 0.0 else None


def _safe_rss(value: Any) -> int | None:
    return value if type(value) is int and value >= 0 else None


def _invalid_executor_result(raw: Mapping[str, Any], message: str) -> dict[str, Any]:
    """Return a JSON-safe terminal result for malformed executor output."""
    return {
        "status": "invalid",
        "valid": False,
        "elapsed_seconds": _safe_elapsed(raw.get("elapsed_seconds")) or 0.0,
        "peak_rss_bytes": _safe_rss(raw.get("peak_rss_bytes")) or 0,
        "error": message,
    }


def _normalize_executor_result(
    raw: Mapping[str, Any], *, label: str
) -> dict[str, Any]:
    """Validate external output before it enters an immutable checkpoint hash."""
    status = raw.get("status")
    if status not in core.CHECKPOINT_STATUSES:
        return _exception_result(
            "exception", ValueError(f"executor returned an unrecognized {label} status")
        )
    elapsed = _safe_elapsed(raw.get("elapsed_seconds", 0.0))
    rss = _safe_rss(raw.get("peak_rss_bytes", 0))
    if elapsed is None or rss is None:
        return _invalid_executor_result(raw, f"malformed {label} resource metrics")
    if "error" in raw and not isinstance(raw["error"], str):
        return _invalid_executor_result(raw, f"malformed {label} error field")

    normalized: dict[str, Any] = {
        "status": status,
        "valid": raw.get("valid"),
        "elapsed_seconds": elapsed,
        "peak_rss_bytes": rss,
    }
    rss_scope = raw.get("rss_scope")
    if rss_scope is not None:
        if rss_scope not in RSS_SCOPES:
            return _invalid_executor_result(raw, f"malformed {label} RSS scope")
        normalized["rss_scope"] = rss_scope
    if "decoder_cache_hit" in raw:
        if type(raw["decoder_cache_hit"]) is not bool:
            return _invalid_executor_result(raw, f"malformed {label} cache evidence")
        normalized["decoder_cache_hit"] = raw["decoder_cache_hit"]
    for name in (
        "decoder_build_seconds", "decode_seconds", "worker_startup_seconds"
    ):
        if name in raw:
            value = _safe_elapsed(raw[name])
            if value is None:
                return _invalid_executor_result(raw, f"malformed {label} timing evidence")
            normalized[name] = value
    if status == "complete":
        try:
            probabilities = _validate_probabilities(raw.get("probabilities"))
            logical_mask = raw.get("logical_mask")
            if (
                type(logical_mask) is not int
                or not 0 <= logical_mask < 4
                or probabilities[logical_mask] != max(probabilities)
            ):
                raise ValueError("logical_mask is not a posterior maximizer")
            log_joint = raw.get("log_joint_syndrome")
            if (
                isinstance(log_joint, bool)
                or not isinstance(log_joint, (int, float))
                or not math.isfinite(log_joint)
            ):
                raise ValueError("log_joint_syndrome is not finite")
            if raw.get("valid") is not True:
                raise ValueError("complete output is not marked valid")
        except (TypeError, ValueError) as error:
            return _invalid_executor_result(raw, f"malformed {label} output: {error}")
        normalized.update(
            probabilities=probabilities,
            logical_mask=logical_mask,
            log_joint_syndrome=float(log_joint),
        )
    else:
        if raw.get("valid") is not False:
            return _invalid_executor_result(
                raw, f"non-complete {label} output is not marked invalid"
            )
        if {"probabilities", "logical_mask", "log_joint_syndrome"} & raw.keys():
            return _invalid_executor_result(
                raw, f"non-complete {label} output contains decoded values"
            )
        if isinstance(raw.get("error"), str):
            normalized["error"] = raw["error"]
    return normalized


def _apply_frozen_rss_policy(
    raw: Mapping[str, Any],
    execution: Mapping[str, Any],
    *,
    label: str,
) -> dict[str, Any]:
    """Reject a completed result that reports exceeding the frozen RSS cap.

    The phase-2 runner does not install an OS-specific memory limiter.  The
    cluster supervisor must enforce the cap while the adapter runs.  This
    check ensures that a completed result observed above the preregistered cap
    can never be counted as a scientific result.
    """
    result = dict(raw)
    policy = execution.get("resource_policy")
    if not isinstance(policy, Mapping):
        raise IdentityMismatchError("execution manifest is missing resource policy")
    maximum = policy.get("max_rss_bytes")
    if type(maximum) is not int or maximum <= 0:
        raise IdentityMismatchError("frozen memory limit is invalid")
    if result.get("status") == "complete" and result.get("peak_rss_bytes", 0) > maximum:
        rejected = {
            "status": "memory_limit",
            "valid": False,
            "elapsed_seconds": result.get("elapsed_seconds", 0.0),
            "peak_rss_bytes": result["peak_rss_bytes"],
            "error": (
                f"{label} reported peak RSS {result['peak_rss_bytes']} above "
                f"the frozen {maximum}-byte cap"
            ),
        }
        if "rss_scope" in result:
            rejected["rss_scope"] = result["rss_scope"]
        if "worker_startup_seconds" in result:
            rejected["worker_startup_seconds"] = result["worker_startup_seconds"]
        return rejected
    return result


def _phase2_level_checkpoint(
    raw: Mapping[str, Any],
    *,
    execution: Mapping[str, Any],
    task: Mapping[str, Any],
    record_meta: Mapping[str, Any],
    chi: int,
    executor_identity: Mapping[str, str],
) -> dict[str, Any]:
    raw = _normalize_executor_result(raw, label="TN")
    raw = _apply_frozen_rss_policy(raw, execution, label="TN")
    status = raw.get("status")
    valid = raw.get("valid")
    kwargs: dict[str, Any] = {}
    if status == "complete":
        kwargs = {
            "probabilities": raw.get("probabilities"),
            "logical_mask": raw.get("logical_mask"),
            "log_joint_syndrome": raw.get("log_joint_syndrome"),
        }
    checkpoint = core.make_level_checkpoint(
        manifest_id=execution["execution_manifest_id"],
        task_id=task["task_id"],
        record_id=record_meta["record_id"],
        chi=chi,
        cutoff=task["cutoff"],
        status=status,
        valid=valid,
        elapsed_seconds=raw.get("elapsed_seconds", 0.0),
        peak_rss_bytes=raw.get("peak_rss_bytes", 0),
        **kwargs,
    )
    checkpoint.pop("checkpoint_id")
    checkpoint.update(
        {
            "phase2_schema_version": SCHEMA_VERSION,
            "corpus_id": execution["corpus_id"],
            "source_bundle_id": execution["source_bundle_id"],
            "runtime_profile_id": execution["runtime_profile_id"],
            "shot_index": record_meta["shot_index"],
            "record_content_sha256": record_meta["record_content_sha256"],
            **executor_identity,
        }
    )
    if status != "complete" and isinstance(raw.get("error"), str):
        checkpoint["error"] = raw["error"]
    if "rss_scope" in raw:
        checkpoint["rss_scope"] = raw["rss_scope"]
    for name in (
        "decoder_cache_hit", "decoder_build_seconds", "decode_seconds",
        "worker_startup_seconds",
    ):
        if name in raw:
            checkpoint[name] = raw[name]
    checkpoint["checkpoint_id"] = core.identity_hash(checkpoint, "checkpoint_id")
    return _validate_level_checkpoint(
        checkpoint,
        execution=execution,
        task=task,
        record_meta=record_meta,
        chi=chi,
        executor_identity=executor_identity,
    )


def _validate_level_checkpoint(
    checkpoint: Mapping[str, Any],
    *,
    execution: Mapping[str, Any],
    task: Mapping[str, Any],
    record_meta: Mapping[str, Any],
    chi: int,
    executor_identity: Mapping[str, str],
) -> dict[str, Any]:
    result = core.validate_level_checkpoint(
        checkpoint,
        manifest_id=execution["execution_manifest_id"],
        task_id=task["task_id"],
        record_id=record_meta["record_id"],
        chi=chi,
        cutoff=task["cutoff"],
    )
    expected = {
        "phase2_schema_version": SCHEMA_VERSION,
        "corpus_id": execution["corpus_id"],
        "source_bundle_id": execution["source_bundle_id"],
        "runtime_profile_id": execution["runtime_profile_id"],
        "shot_index": record_meta["shot_index"],
        "record_content_sha256": record_meta["record_content_sha256"],
    }
    if any(result.get(key) != value for key, value in expected.items()):
        raise ValueError("level checkpoint phase-2 identity mismatch")
    if any(result.get(key) != value for key, value in executor_identity.items()):
        raise ValueError("level checkpoint executor identity mismatch")
    if "rss_scope" in result and result["rss_scope"] not in RSS_SCOPES:
        raise ValueError("level checkpoint RSS scope is invalid")
    _validate_checkpoint_executor_metrics(result, executor_identity)
    return result


def _validate_probabilities(probabilities: Any) -> list[float]:
    if not isinstance(probabilities, (list, tuple)) or len(probabilities) != 4:
        raise ValueError("exact probabilities must contain four entries")
    values = []
    for value in probabilities:
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError("exact probabilities must be finite")
        values.append(float(value))
    if min(values) < -core.NEGATIVE_PROBABILITY_TOLERANCE:
        raise ValueError("exact probabilities are below the negative tolerance")
    if abs(sum(values) - 1.0) > core.NORMALIZATION_TOLERANCE:
        raise ValueError("exact probabilities are not normalized")
    return values


def _exact_checkpoint(
    raw: Mapping[str, Any],
    *,
    execution: Mapping[str, Any],
    task: Mapping[str, Any],
    record_meta: Mapping[str, Any],
    executor_identity: Mapping[str, str],
) -> dict[str, Any]:
    raw = _normalize_executor_result(raw, label="exact")
    raw = _apply_frozen_rss_policy(raw, execution, label="exact")
    status = raw["status"]
    document: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "kind": "adaptive_chi_exact_reference",
        "execution_manifest_id": execution["execution_manifest_id"],
        "corpus_id": execution["corpus_id"],
        "source_bundle_id": execution["source_bundle_id"],
        "runtime_profile_id": execution["runtime_profile_id"],
        "task_id": task["task_id"],
        "shot_index": record_meta["shot_index"],
        "record_id": record_meta["record_id"],
        "record_content_sha256": record_meta["record_content_sha256"],
        "status": status,
        "valid": raw.get("valid"),
        "elapsed_seconds": raw.get("elapsed_seconds", 0.0),
        "peak_rss_bytes": raw.get("peak_rss_bytes", 0),
        **executor_identity,
    }
    if status == "complete":
        document.update(
            probabilities=list(raw.get("probabilities", [])),
            logical_mask=raw.get("logical_mask"),
            log_joint_syndrome=raw.get("log_joint_syndrome"),
        )
    elif isinstance(raw.get("error"), str):
        document["error"] = raw["error"]
    if "rss_scope" in raw:
        document["rss_scope"] = raw["rss_scope"]
    for name in (
        "decoder_cache_hit", "decoder_build_seconds", "decode_seconds",
        "worker_startup_seconds",
    ):
        if name in raw:
            document[name] = raw[name]
    document["exact_checkpoint_id"] = core.identity_hash(
        document, "exact_checkpoint_id"
    )
    return _validate_exact_checkpoint(
        document,
        execution=execution,
        task=task,
        record_meta=record_meta,
        executor_identity=executor_identity,
    )


def _validate_exact_checkpoint(
    document: Mapping[str, Any],
    *,
    execution: Mapping[str, Any],
    task: Mapping[str, Any],
    record_meta: Mapping[str, Any],
    executor_identity: Mapping[str, str],
) -> dict[str, Any]:
    result = dict(document)
    identity = result.get("exact_checkpoint_id")
    if not isinstance(identity, str) or identity != core.identity_hash(result, "exact_checkpoint_id"):
        raise ValueError("exact checkpoint identity mismatch")
    expected = {
        "schema_version": SCHEMA_VERSION,
        "kind": "adaptive_chi_exact_reference",
        "execution_manifest_id": execution["execution_manifest_id"],
        "corpus_id": execution["corpus_id"],
        "source_bundle_id": execution["source_bundle_id"],
        "runtime_profile_id": execution["runtime_profile_id"],
        "task_id": task["task_id"],
        "shot_index": record_meta["shot_index"],
        "record_id": record_meta["record_id"],
        "record_content_sha256": record_meta["record_content_sha256"],
    }
    if any(result.get(key) != value for key, value in expected.items()):
        raise ValueError("exact checkpoint phase-2 identity mismatch")
    if any(result.get(key) != value for key, value in executor_identity.items()):
        raise ValueError("exact checkpoint executor identity mismatch")
    elapsed = result.get("elapsed_seconds")
    rss = result.get("peak_rss_bytes")
    if isinstance(elapsed, bool) or not isinstance(elapsed, (int, float)) or not math.isfinite(elapsed) or elapsed < 0:
        raise ValueError("exact elapsed_seconds must be finite and nonnegative")
    if isinstance(rss, bool) or not isinstance(rss, int) or rss < 0:
        raise ValueError("exact peak_rss_bytes must be nonnegative")
    if "rss_scope" in result and result["rss_scope"] not in RSS_SCOPES:
        raise ValueError("exact RSS scope is invalid")
    _validate_checkpoint_executor_metrics(result, executor_identity)
    status = result.get("status")
    if status == "complete":
        if result.get("valid") is not True:
            raise ValueError("complete exact reference must be valid")
        probabilities = _validate_probabilities(result.get("probabilities"))
        logical_mask = result.get("logical_mask")
        if type(logical_mask) is not int or not 0 <= logical_mask < 4 or probabilities[logical_mask] != max(probabilities):
            raise ValueError("exact logical_mask must be a posterior maximizer")
        log_joint = result.get("log_joint_syndrome")
        if isinstance(log_joint, bool) or not isinstance(log_joint, (int, float)) or not math.isfinite(log_joint):
            raise ValueError("exact log_joint_syndrome must be finite")
    elif status in core.CHECKPOINT_STATUSES - {"complete"}:
        if result.get("valid") is not False:
            raise ValueError("non-complete exact reference must be invalid")
        if {"probabilities", "logical_mask", "log_joint_syndrome"} & result.keys():
            raise ValueError("non-complete exact reference cannot contain a decoded class")
    else:
        raise ValueError("exact reference has an unrecognized status")
    return result


def _exact_audit(
    adaptive_decision: Mapping[str, Any],
    exact: Mapping[str, Any] | None,
    policy: Mapping[str, Any],
) -> dict[str, Any]:
    if exact is None:
        return {"status": "not_required", "passed": None, "false_accept": False}
    if exact.get("status") != "complete" or exact.get("valid") is not True:
        return {
            "status": "unavailable",
            "passed": False,
            "false_accept": None,
            "reason": "invalid_exact_reference",
        }
    exact_probabilities = list(exact["probabilities"])
    maximum = max(exact_probabilities)
    maximizers = [
        index
        for index, value in enumerate(exact_probabilities)
        if maximum - value <= policy["ml_tie_tolerance"]
    ]
    audit: dict[str, Any] = {
        "status": "complete",
        "exact_maximizers": maximizers,
        "ml_tie_tolerance": policy["ml_tie_tolerance"],
        "evaluated": adaptive_decision.get("status") == "accepted",
    }
    if not audit["evaluated"]:
        audit.update(passed=None, false_accept=False)
        return audit
    posterior_linf = max(
        abs(a - b)
        for a, b in zip(adaptive_decision["probabilities"], exact_probabilities)
    )
    log_joint_error = abs(
        adaptive_decision["log_joint_syndrome"] - exact["log_joint_syndrome"]
    )
    class_match = adaptive_decision["logical_mask"] in maximizers
    posterior_pass = posterior_linf <= policy["posterior_linf_gate"]
    log_joint_pass = log_joint_error <= policy["log_joint_gate"]
    passed = class_match and posterior_pass and log_joint_pass
    audit.update(
        ml_tie_aware_match=class_match,
        posterior_linf=posterior_linf,
        posterior_gate_passed=posterior_pass,
        log_joint_error=log_joint_error,
        log_joint_gate_passed=log_joint_pass,
        passed=passed,
        false_accept=not passed,
    )
    return audit


def _build_decision(
    *,
    execution: Mapping[str, Any],
    task: Mapping[str, Any],
    record_meta: Mapping[str, Any],
    record: Mapping[str, Any],
    levels: list[Mapping[str, Any]],
    adaptive_decision: Mapping[str, Any],
    exact: Mapping[str, Any] | None,
    executor_identity: Mapping[str, str],
) -> dict[str, Any]:
    audit = _exact_audit(adaptive_decision, exact, execution["exact_audit"])
    document: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "kind": "adaptive_chi_record_decision",
        "execution_manifest_id": execution["execution_manifest_id"],
        "corpus_id": execution["corpus_id"],
        "source_bundle_id": execution["source_bundle_id"],
        "runtime_profile_id": execution["runtime_profile_id"],
        "task_id": task["task_id"],
        "shot_index": record_meta["shot_index"],
        "record_id": record_meta["record_id"],
        "record_content_sha256": record_meta["record_content_sha256"],
        "record_output": _relative_record_output(task, record_meta),
        "forced_anchor": record_meta["forced_anchor"],
        "tested_chi": [level["chi"] for level in levels],
        "level_checkpoint_ids": [level["checkpoint_id"] for level in levels],
        "exact_checkpoint_id": exact.get("exact_checkpoint_id") if exact else None,
        "adaptive_decision": dict(adaptive_decision),
        "exact_audit": audit,
        **executor_identity,
    }
    if adaptive_decision.get("status") != "accepted":
        document.update(
            final_status="unknown",
            unknown_reason=adaptive_decision.get("reason", "not_converged"),
        )
    elif task["exact_required"] and audit["status"] != "complete":
        document.update(final_status="unknown", unknown_reason="invalid_exact_reference")
    else:
        document.update(
            final_status="accepted",
            logical_mask=adaptive_decision["logical_mask"],
            failure=record["true_mask"] != adaptive_decision["logical_mask"],
        )
    document["decision_id"] = core.identity_hash(document, "decision_id")
    return document


def _present_decision(document: Mapping[str, Any], root: Path) -> dict[str, Any]:
    result = dict(document)
    result["record_directory"] = str((root / result["record_output"]).resolve())
    return result


def _evaluate_record_state(
    levels: list[Mapping[str, Any]],
    *,
    task: Mapping[str, Any],
    record_meta: Mapping[str, Any],
    exact: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """Stop immediately when a required exact call cannot be audited."""
    if task["exact_required"] and exact is not None and exact.get("status") != "complete":
        status = exact.get("status")
        reason = (
            f"exact_{status}"
            if status in ("timeout", "memory_limit")
            else "invalid_exact_reference"
        )
        return {"status": "unknown", "reason": reason}
    return core.evaluate_stopping(
        levels, task["chi_ladder"], forced_anchor=record_meta["forced_anchor"]
    )


def _load_levels(
    record_dir: Path,
    *,
    execution: Mapping[str, Any],
    task: Mapping[str, Any],
    record_meta: Mapping[str, Any],
    executor_identity: Mapping[str, str],
) -> list[dict[str, Any]]:
    levels = []
    missing_seen = False
    for chi in task["chi_ladder"]:
        path = record_dir / f"chi_{chi:04d}.json"
        if not path.exists():
            missing_seen = True
            continue
        if missing_seen:
            raise ValueError("level checkpoints do not form a contiguous chi prefix")
        levels.append(
            _validate_level_checkpoint(
                _read_json(path),
                execution=execution,
                task=task,
                record_meta=record_meta,
                chi=chi,
                executor_identity=executor_identity,
            )
        )
    return levels


def _validate_decision(
    document: Mapping[str, Any],
    *,
    root: Path,
    execution: Mapping[str, Any],
    task: Mapping[str, Any],
    record_meta: Mapping[str, Any],
    record: Mapping[str, Any],
    executor_identity: Mapping[str, str],
) -> dict[str, Any]:
    result = dict(document)
    identity = result.get("decision_id")
    if not isinstance(identity, str) or identity != core.identity_hash(result, "decision_id"):
        raise ValueError("decision identity mismatch")
    record_dir = root / _relative_record_output(task, record_meta)
    levels = _load_levels(
        record_dir,
        execution=execution,
        task=task,
        record_meta=record_meta,
        executor_identity=executor_identity,
    )
    exact = None
    exact_path = record_dir / "exact.json"
    if exact_path.exists():
        exact = _validate_exact_checkpoint(
            _read_json(exact_path),
            execution=execution,
            task=task,
            record_meta=record_meta,
            executor_identity=executor_identity,
        )
    adaptive = _evaluate_record_state(
        levels, task=task, record_meta=record_meta, exact=exact
    )
    rebuilt = _build_decision(
        execution=execution,
        task=task,
        record_meta=record_meta,
        record=record,
        levels=levels,
        adaptive_decision=adaptive,
        exact=exact,
        executor_identity=executor_identity,
    )
    if result != rebuilt:
        raise ValueError("decision derived identity mismatch")
    return result


def run_record(
    root: os.PathLike[str] | str,
    sampling_manifest: Mapping[str, Any],
    *,
    task_id: int,
    record_id: str,
    executor: Executor,
    dry_run: bool = False,
) -> dict[str, Any]:
    root_path = Path(root).expanduser().resolve()
    plan = plan_record(
        root_path, sampling_manifest, task_id=task_id, record_id=record_id
    )
    if dry_run:
        return plan
    execution = execution_manifest.load_execution_manifest(root_path, sampling_manifest)
    corpus = execution_manifest.load_frozen_corpus(root_path, sampling_manifest)
    task, record_meta, corpus_entry = _find_record(
        execution, corpus, task_id, record_id
    )
    record = corpus_entry["record"]
    executor_identity = _executor_checkpoint_identity(executor, execution)
    identity = _execution_identity(
        execution, task, record_meta, executor_identity=executor_identity
    )
    verify_environment = getattr(executor, "verify_environment_only", None)
    if not callable(verify_environment):
        verify_environment = getattr(executor, "verify_execution_environment", None)
    if callable(verify_environment):
        verify_environment(identity)
    record_dir = root_path / plan["record_output"]
    decision_path = record_dir / "decision.json"
    if decision_path.exists():
        decision = _validate_decision(
            _read_json(decision_path),
            root=root_path,
            execution=execution,
            task=task,
            record_meta=record_meta,
            record=record,
            executor_identity=executor_identity,
        )
        return _present_decision(decision, root_path)

    preflight = getattr(executor, "preflight", None)
    if callable(preflight):
        preflight(identity)

    exact = None
    if task["exact_required"]:
        exact_path = record_dir / "exact.json"
        if exact_path.exists():
            exact = _validate_exact_checkpoint(
                _read_json(exact_path),
                execution=execution,
                task=task,
                record_meta=record_meta,
                executor_identity=executor_identity,
            )
        else:
            raw_exact = _call_executor(
                lambda: executor.run_exact(record, execution_identity=identity)
            )
            exact = _exact_checkpoint(
                raw_exact,
                execution=execution,
                task=task,
                record_meta=record_meta,
                executor_identity=executor_identity,
            )
            _write_once(exact_path, exact, "exact checkpoint")

    levels = _load_levels(
        record_dir,
        execution=execution,
        task=task,
        record_meta=record_meta,
        executor_identity=executor_identity,
    )
    adaptive = _evaluate_record_state(
        levels, task=task, record_meta=record_meta, exact=exact
    )
    while adaptive["status"] in ("continue", "continue_anchor"):
        chi = adaptive["next_chi"]
        raw = _call_executor(
            lambda chi=chi: executor.run_tn(
                record,
                chi=chi,
                cutoff=task["cutoff"],
                execution_identity=identity,
            )
        )
        checkpoint = _phase2_level_checkpoint(
            raw,
            execution=execution,
            task=task,
            record_meta=record_meta,
            chi=chi,
            executor_identity=executor_identity,
        )
        _write_once(
            record_dir / f"chi_{chi:04d}.json", checkpoint, "level checkpoint"
        )
        levels.append(checkpoint)
        adaptive = _evaluate_record_state(
            levels, task=task, record_meta=record_meta, exact=exact
        )

    decision = _build_decision(
        execution=execution,
        task=task,
        record_meta=record_meta,
        record=record,
        levels=levels,
        adaptive_decision=adaptive,
        exact=exact,
        executor_identity=executor_identity,
    )
    _write_once(decision_path, decision, "record decision")
    return _present_decision(decision, root_path)


def _toml_value(value: Any) -> str:
    if type(value) is bool:
        return "true" if value else "false"
    if type(value) is int:
        return str(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("cannot serialize nonfinite TOML value")
        return repr(value)
    if isinstance(value, str):
        return json.dumps(value)
    if isinstance(value, list):
        return "[" + ", ".join(_toml_value(item) for item in value) + "]"
    raise TypeError(f"unsupported TOML value: {type(value).__name__}")


def _write_record_toml(path: Path, record: Mapping[str, Any]) -> None:
    lines = [f"{key} = {_toml_value(value)}" for key, value in sorted(record.items())]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _read_toml(path: Path) -> dict[str, Any]:
    try:
        import tomllib  # type: ignore

        with path.open("rb") as stream:
            return dict(tomllib.load(stream))
    except ImportError:
        try:
            import toml  # type: ignore
        except ImportError as error:
            raise RuntimeError("Python 3.11+ or the toml package is required") from error
        return dict(toml.loads(path.read_text(encoding="utf-8")))


class JuliaExecutor:
    """Thin subprocess adapter to existing Julia TN and dense-exact entry points."""

    executor_mode = "one_shot"

    def __init__(
        self,
        *,
        repository_root: os.PathLike[str] | str,
        julia: str = "julia",
    ):
        self.repository_root = Path(repository_root).expanduser().resolve()
        self.julia = julia
        self.adapter = Path(__file__).resolve().with_name("julia_adapter.jl")
        self._verified_profiles: set[tuple[str, str]] = set()

    def verify_execution_environment(
        self, execution_identity: Mapping[str, Any]
    ) -> None:
        source = execution_identity.get("source_identity")
        runtime = execution_identity.get("runtime_identity")
        if not isinstance(source, Mapping) or not isinstance(runtime, Mapping):
            raise IdentityMismatchError(
                "execution identity is missing frozen source/runtime data"
            )
        try:
            execution_manifest.verify_source_identity(self.repository_root, source)
        except ValueError as error:
            raise IdentityMismatchError(str(error)) from error
        if runtime.get("python_version") != platform.python_version():
            raise IdentityMismatchError(
                "current Python differs from the frozen runtime identity"
            )

    @staticmethod
    def verify_julia_runtime(
        actual: Mapping[str, Any], expected: Mapping[str, Any]
    ) -> None:
        fields = {
            "julia_version",
            "blas_vendor",
            "cpu_model",
            "julia_threads",
            "blas_threads",
        }
        if not isinstance(actual, Mapping) or set(actual) != fields:
            raise IdentityMismatchError("Julia adapter returned incomplete runtime identity")
        if any(actual.get(name) != expected.get(name) for name in fields):
            raise IdentityMismatchError(
                "Julia adapter runtime identity differs from the frozen profile"
            )

    def _probe_command(self, output_path: Path) -> list[str]:
        return [
            self.julia,
            "--project=" + str(self.repository_root),
            "--startup-file=no",
            "--threads=1",
            str(self.adapter),
            "probe",
            str(output_path),
        ]

    @staticmethod
    def _wall_timeout(
        mode: str, chi: int, execution_identity: Mapping[str, Any]
    ) -> float | None:
        policy = execution_identity["resource_policy"]
        value = (
            policy["exact_wall_seconds"]
            if mode == "exact"
            else policy["wall_seconds_by_chi"].get(str(chi))
        )
        return None if value is None else float(value)

    def preflight(self, execution_identity: Mapping[str, Any]) -> None:
        """Verify frozen source and Julia runtime before any scientific work."""
        self.verify_execution_environment(execution_identity)
        key = (
            str(execution_identity.get("source_bundle_id")),
            str(execution_identity.get("runtime_profile_id")),
        )
        if key in self._verified_profiles:
            return
        probe_timeout = float(
            execution_identity["resource_policy"]["preflight_wall_seconds"]
        )
        with tempfile.TemporaryDirectory(prefix="tnml-adaptive-probe-") as temporary:
            output_path = Path(temporary) / "runtime.toml"
            environment = os.environ.copy()
            blas_threads = str(execution_identity["runtime_identity"]["blas_threads"])
            environment["TNML_BLAS_THREADS"] = blas_threads
            environment["OPENBLAS_NUM_THREADS"] = blas_threads
            environment["MKL_NUM_THREADS"] = blas_threads
            try:
                completed = subprocess.run(
                    self._probe_command(output_path),
                    check=False,
                    capture_output=True,
                    text=True,
                    timeout=probe_timeout,
                    env=environment,
                )
            except subprocess.TimeoutExpired as error:
                raise IdentityMismatchError("Julia runtime preflight timed out") from error
            if completed.returncode != 0 or not output_path.is_file():
                raise IdentityMismatchError(
                    f"Julia runtime preflight failed ({completed.returncode}): "
                    f"{completed.stderr[-2000:]}"
                )
            evidence = _read_toml(output_path).get("runtime_identity")
            self.verify_julia_runtime(
                evidence, execution_identity["runtime_identity"]
            )
        self._verified_profiles.add(key)

    def command(
        self,
        mode: str,
        record_path: Path,
        output_path: Path,
        *,
        L: int,
        p: float,
        chi: int,
        cutoff: float,
    ) -> list[str]:
        if mode not in ("tn", "exact"):
            raise ValueError("Julia adapter mode must be tn or exact")
        return [
            self.julia,
            "--project=" + str(self.repository_root),
            "--startup-file=no",
            "--threads=1",
            str(self.adapter),
            mode,
            str(record_path),
            str(output_path),
            str(L),
            repr(p),
            str(chi),
            repr(cutoff),
        ]

    def _invoke(
        self,
        mode: str,
        record: Mapping[str, Any],
        *,
        chi: int,
        cutoff: float,
        execution_identity: Mapping[str, Any],
    ) -> dict[str, Any]:
        self.preflight(execution_identity)
        with tempfile.TemporaryDirectory(prefix="tnml-adaptive-") as temporary:
            directory = Path(temporary)
            record_path, output_path = directory / "record.toml", directory / "result.toml"
            _write_record_toml(record_path, record)
            command = self.command(
                mode,
                record_path,
                output_path,
                L=execution_identity["L"],
                p=execution_identity["p"],
                chi=chi,
                cutoff=cutoff,
            )
            environment = os.environ.copy()
            blas_threads = str(execution_identity["runtime_identity"]["blas_threads"])
            environment["TNML_BLAS_THREADS"] = blas_threads
            environment["OPENBLAS_NUM_THREADS"] = blas_threads
            environment["MKL_NUM_THREADS"] = blas_threads
            timeout = self._wall_timeout(mode, chi, execution_identity)
            try:
                completed = subprocess.run(
                    command,
                    check=False,
                    capture_output=True,
                    text=True,
                    timeout=timeout,
                    env=environment,
                )
            except subprocess.TimeoutExpired as error:
                return _exception_result("timeout", error)
            result = None
            if output_path.is_file():
                result = _read_toml(output_path)
                runtime = result.pop("runtime_identity", None)
                self.verify_julia_runtime(
                    runtime, execution_identity["runtime_identity"]
                )
            if completed.returncode != 0:
                if (
                    result is not None
                    and result.get("status") == "memory_limit"
                    and result.get("valid") is False
                ):
                    return result
                return {
                    "status": "exception",
                    "valid": False,
                    "elapsed_seconds": 0.0,
                    "peak_rss_bytes": 0,
                    "error": (
                        f"Julia adapter exited {completed.returncode}: "
                        f"{completed.stderr[-2000:]}"
                    ),
                }
            if result is None:
                return {
                    "status": "exception",
                    "valid": False,
                    "elapsed_seconds": 0.0,
                    "peak_rss_bytes": 0,
                    "error": "Julia adapter exited successfully without a result",
                }
            return result

    def run_tn(
        self,
        record: Mapping[str, Any],
        *,
        chi: int,
        cutoff: float,
        execution_identity: Mapping[str, Any],
    ) -> Mapping[str, Any]:
        return self._invoke(
            "tn",
            record,
            chi=chi,
            cutoff=cutoff,
            execution_identity=execution_identity,
        )

    def run_exact(
        self,
        record: Mapping[str, Any],
        *,
        execution_identity: Mapping[str, Any],
    ) -> Mapping[str, Any]:
        return self._invoke(
            "exact",
            record,
            chi=1,
            cutoff=0.0,
            execution_identity=execution_identity,
        )


def _sampling_manifest_from_file(path: os.PathLike[str] | str) -> dict[str, Any]:
    value = json.loads(Path(path).expanduser().resolve().read_text(encoding="utf-8"))
    if not isinstance(value, Mapping):
        raise ValueError("sampling manifest must contain a JSON object")
    return manifest.validate_manifest(value)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Plan or run one record from a frozen adaptive-chi study"
    )
    parser.add_argument("--study-root", required=True)
    parser.add_argument("--sampling-manifest", required=True)
    parser.add_argument("--task-id", required=True, type=int)
    parser.add_argument("--record-id", required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true")
    mode.add_argument("--execute", action="store_true")
    parser.add_argument(
        "--repository-root", default=str(Path(__file__).resolve().parents[3])
    )
    parser.add_argument("--julia", default="julia")
    arguments = parser.parse_args(argv)
    sampling = _sampling_manifest_from_file(arguments.sampling_manifest)
    if arguments.dry_run:
        result = plan_record(
            arguments.study_root,
            sampling,
            task_id=arguments.task_id,
            record_id=arguments.record_id,
        )
    else:
        executor = JuliaExecutor(
            repository_root=arguments.repository_root,
            julia=arguments.julia,
        )
        result = run_record(
            arguments.study_root,
            sampling,
            task_id=arguments.task_id,
            record_id=arguments.record_id,
            executor=executor,
        )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
