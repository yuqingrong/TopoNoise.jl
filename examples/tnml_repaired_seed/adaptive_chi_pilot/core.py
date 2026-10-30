#!/usr/bin/env python3
"""Pure policy and identity helpers for the adaptive-chi phase-1 pilot."""

from __future__ import annotations

from collections import Counter
import hashlib
import json
import math
import re
from typing import Any, Iterable, Mapping, Sequence


SCHEMA_VERSION = 1
POSTERIOR_TOLERANCE = 1.0e-6
LOG_JOINT_TOLERANCE = 1.0e-6
NORMALIZATION_TOLERANCE = 1.0e-8
NEGATIVE_PROBABILITY_TOLERANCE = 1.0e-10
MINIMUM_MARGIN = 1.0e-8
DRIFT_FLOOR = 1.0e-12
DRIFT_GROWTH = 1.25
RECORD_ID_PATTERN = re.compile(r"^[0-9a-f]{64}$")
P_STRATA = ("0.0060", "0.0100", "0.0200", "0.0300", "0.0350", "0.0425")
CHECKPOINT_STATUSES = frozenset(
    {
        "complete",
        "exception",
        "interrupted",
        "invalid",
        "memory_limit",
        "missing_checkpoint",
        "timeout",
    }
)


def _canonical_json(document: Mapping[str, Any]) -> bytes:
    return json.dumps(
        document,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")


def identity_hash(document: Mapping[str, Any], identity_key: str) -> str:
    """Return the SHA-256 identity after omitting the identity field itself."""
    payload = {key: value for key, value in document.items() if key != identity_key}
    return hashlib.sha256(_canonical_json(payload)).hexdigest()


def _finite_number(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite number")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite number")
    return result


def _positive_integer(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _validate_probabilities(values: Any) -> list[float]:
    if not isinstance(values, (list, tuple)) or len(values) != 4:
        raise ValueError("probabilities must contain exactly four entries")
    probabilities = [
        _finite_number(value, f"probabilities[{index}]")
        for index, value in enumerate(values)
    ]
    if min(probabilities) < -NEGATIVE_PROBABILITY_TOLERANCE:
        raise ValueError("probabilities contain an entry below the negative tolerance")
    if abs(sum(probabilities) - 1.0) > NORMALIZATION_TOLERANCE:
        raise ValueError("probabilities are not normalized within tolerance")
    return probabilities


def _unique_argmax(probabilities: Sequence[float]) -> int | None:
    maximum = max(probabilities)
    maximizers = [index for index, value in enumerate(probabilities) if value == maximum]
    return maximizers[0] if len(maximizers) == 1 else None


def _validate_logical_mask(value: Any, probabilities: Sequence[float]) -> int:
    if type(value) is not int or not 0 <= value < 4:
        raise ValueError("logical_mask must be an integer in 0:3")
    if probabilities[value] != max(probabilities):
        raise ValueError("logical_mask must be a posterior maximizer")
    return value


def _validate_ladder_and_levels(
    levels: Sequence[Mapping[str, Any]], ladder: Sequence[int]
) -> tuple[list[Mapping[str, Any]], list[int]]:
    ladder_values = [_positive_integer(value, "chi ladder entry") for value in ladder]
    if not ladder_values:
        raise ValueError("chi ladder must not be empty")
    if any(right <= left for left, right in zip(ladder_values, ladder_values[1:])):
        raise ValueError("chi ladder must be strictly increasing")
    level_values = list(levels)
    if len(level_values) > len(ladder_values):
        raise ValueError("completed levels must form a prefix of the chi ladder")
    observed_chis = [level.get("chi") for level in level_values]
    if observed_chis != ladder_values[: len(level_values)]:
        raise ValueError("completed levels must form a prefix of the chi ladder")
    return level_values, ladder_values


def _eligible_level(level: Mapping[str, Any]) -> dict[str, Any]:
    status = level.get("status")
    if status != "complete":
        return {"status": "unknown", "reason": str(status or "missing_checkpoint")}
    if level.get("valid") is not True:
        return {"status": "unknown", "reason": "invalid"}
    try:
        probabilities = _validate_probabilities(level.get("probabilities"))
        logical_mask = _validate_logical_mask(level.get("logical_mask"), probabilities)
        log_joint = _finite_number(level.get("log_joint_syndrome"), "log_joint_syndrome")
    except (TypeError, ValueError):
        return {"status": "unknown", "reason": "invalid"}
    return {
        "status": "complete",
        "chi": level["chi"],
        "probabilities": probabilities,
        "logical_mask": logical_mask,
        "unique_logical_mask": _unique_argmax(probabilities),
        "log_joint_syndrome": log_joint,
    }


def _comparison_metrics(eligible: Sequence[Mapping[str, Any]], end: int) -> dict[str, Any]:
    first, second, newest = eligible[end - 2 : end + 1]
    d_previous = max(
        abs(a - b) for a, b in zip(second["probabilities"], first["probabilities"])
    )
    d_newest = max(
        abs(a - b) for a, b in zip(newest["probabilities"], second["probabilities"])
    )
    g_previous = abs(second["log_joint_syndrome"] - first["log_joint_syndrome"])
    g_newest = abs(newest["log_joint_syndrome"] - second["log_joint_syndrome"])
    posterior_bound = max(d_previous, d_newest)
    log_joint_bound = max(g_previous, g_newest)
    ordered = sorted(newest["probabilities"], reverse=True)
    margin = ordered[0] - ordered[1]
    unique_masks = [
        item["unique_logical_mask"] for item in (first, second, newest)
    ]
    same_class = (
        all(mask is not None for mask in unique_masks)
        and len(set(unique_masks)) == 1
    )
    converged = (
        same_class
        and posterior_bound <= POSTERIOR_TOLERANCE
        and log_joint_bound <= LOG_JOINT_TOLERANCE
        and margin > max(MINIMUM_MARGIN, 4.0 * posterior_bound)
        and d_newest <= max(DRIFT_FLOOR, DRIFT_GROWTH * d_previous)
        and g_newest <= max(DRIFT_FLOOR, DRIFT_GROWTH * g_previous)
    )
    return {
        "converged": converged,
        "d_previous": d_previous,
        "d_newest": d_newest,
        "g_previous": g_previous,
        "g_newest": g_newest,
        "posterior_bound": posterior_bound,
        "log_joint_bound": log_joint_bound,
        "margin": margin,
        "same_unique_class": same_class,
    }


def _accepted_decision(level: Mapping[str, Any], metrics: Mapping[str, Any]) -> dict[str, Any]:
    logical_mask = level["unique_logical_mask"]
    if logical_mask is None:
        raise ValueError("an accepted level must have a unique logical maximizer")
    return {
        "status": "accepted",
        "accepted_chi": level["chi"],
        "logical_mask": logical_mask,
        "probabilities": list(level["probabilities"]),
        "log_joint_syndrome": level["log_joint_syndrome"],
        **{key: value for key, value in metrics.items() if key != "converged"},
    }


def evaluate_stopping(
    levels: Sequence[Mapping[str, Any]],
    ladder: Sequence[int],
    *,
    forced_anchor: bool = False,
) -> dict[str, Any]:
    """Evaluate the preregistered stopping rule without running a decoder."""
    level_values, ladder_values = _validate_ladder_and_levels(levels, ladder)
    eligible: list[dict[str, Any]] = []
    for level in level_values:
        if (level.get("status") == "invalid" and level.get("valid") is False
                and ("error" not in level or level["error"] == "")):
            # Numerical invalidity at small chi can disappear at larger chi.
            # It cannot contribute to a convergence window or preserve an
            # earlier provisional acceptance during a forced anchor.
            eligible.clear()
            continue
        checked = _eligible_level(level)
        if checked["status"] != "complete":
            return checked
        eligible.append(checked)

    first_acceptance = None
    for end in range(2, len(eligible)):
        metrics = _comparison_metrics(eligible, end)
        if metrics["converged"]:
            first_acceptance = _accepted_decision(eligible[end], metrics)
            break

    at_maximum = len(level_values) == len(ladder_values)
    if level_values and not eligible:
        if at_maximum:
            return {"status": "unknown", "reason": "invalid"}
        return {
            "status": "continue",
            "reason": "retry_numerical_invalidity",
            "next_chi": ladder_values[len(level_values)],
        }
    if forced_anchor and first_acceptance is not None and not at_maximum:
        return {
            "status": "continue_anchor",
            "reason": "forced_high_chi_anchor",
            "next_chi": ladder_values[len(level_values)],
            "provisional_acceptance": first_acceptance,
        }
    if first_acceptance is not None:
        if forced_anchor:
            anchor_final_logical_mask = eligible[-1]["unique_logical_mask"]
            anchor_posterior_linf = max(
                abs(a - b)
                for a, b in zip(
                    first_acceptance["probabilities"], eligible[-1]["probabilities"]
                )
            )
            return {
                **first_acceptance,
                "anchor_completed": True,
                "anchor_final_chi": eligible[-1]["chi"],
                "anchor_final_logical_mask": anchor_final_logical_mask,
                "anchor_final_probabilities": list(eligible[-1]["probabilities"]),
                "anchor_class_changed": (
                    anchor_final_logical_mask is None
                    or first_acceptance["logical_mask"] != anchor_final_logical_mask
                ),
                "anchor_posterior_linf": anchor_posterior_linf,
                "anchor_log_joint_difference": abs(
                    first_acceptance["log_joint_syndrome"]
                    - eligible[-1]["log_joint_syndrome"]
                ),
            }
        return first_acceptance
    if at_maximum:
        return {"status": "unknown", "reason": "max_chi_without_convergence"}
    return {
        "status": "continue",
        "reason": "not_converged",
        "next_chi": ladder_values[len(level_values)],
    }


def summarize_outcomes(outcomes: Iterable[Mapping[str, Any]]) -> dict[str, Any]:
    """Summarize known failures while retaining every nonacceptance as unknown."""
    values = list(outcomes)
    accepted = [item for item in values if item.get("status") == "accepted"]
    for item in accepted:
        if not isinstance(item.get("failure"), bool):
            raise ValueError("each accepted outcome must contain a Boolean failure")
    known_failures = sum(item["failure"] for item in accepted)
    unknown = [item for item in values if item.get("status") != "accepted"]
    reasons = Counter(str(item.get("reason", "unspecified")) for item in unknown)
    requested = len(values)
    lower = known_failures / requested if requested else 0.0
    upper = (known_failures + len(unknown)) / requested if requested else 0.0
    return {
        "requested": requested,
        "accepted": len(accepted),
        "known_failures": known_failures,
        "unknown": len(unknown),
        "failure_interval": [lower, upper],
        "unknown_reasons": dict(sorted(reasons.items())),
    }


def _canonical_probability(value: Any) -> str:
    probability = _finite_number(value, "p")
    if probability < 0.0 or probability > 1.0:
        raise ValueError("p must be between zero and one")
    p_key = f"{probability:.4f}"
    if p_key not in P_STRATA or abs(probability - float(p_key)) > 1.0e-12:
        raise ValueError("p must be one of the six preregistered probability strata")
    return p_key


def rank_record_ids(
    L: int, p: Any, record_ids: Iterable[str], *, count: int
) -> list[dict[str, Any]]:
    """Select anchors by the exact p-stratified SHA-256 ranking."""
    size = _positive_integer(L, "L")
    if isinstance(count, bool) or not isinstance(count, int) or count < 0:
        raise ValueError("count must be a nonnegative integer")
    identifiers = list(record_ids)
    if any(not isinstance(value, str) or RECORD_ID_PATTERN.fullmatch(value) is None for value in identifiers):
        raise ValueError("record IDs must be 64 lowercase hexadecimal characters")
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("record IDs contain a duplicate")
    if count > len(identifiers):
        raise ValueError("count exceeds the number of record IDs")
    p_key = _canonical_probability(p)
    ranked = []
    for record_id in identifiers:
        message = (
            f"tnml-adaptive-chi-v1|L={size}|p={p_key}|record_id={record_id}"
        ).encode("utf-8")
        digest = hashlib.sha256(message).digest()
        ranked.append((int.from_bytes(digest, "big"), record_id, digest.hex()))
    ranked.sort(key=lambda item: (item[0], item[1]))
    return [
        {
            "rank": rank,
            "record_id": record_id,
            "digest_hex": digest_hex,
            "digest_integer": str(digest_integer),
        }
        for rank, (digest_integer, record_id, digest_hex) in enumerate(
            ranked[:count], start=1
        )
    ]


def make_level_checkpoint(
    *,
    manifest_id: str,
    task_id: int,
    record_id: str,
    chi: int,
    cutoff: float,
    status: str,
    valid: bool,
    elapsed_seconds: float,
    peak_rss_bytes: int,
    probabilities: Sequence[float] | None = None,
    logical_mask: int | None = None,
    log_joint_syndrome: float | None = None,
) -> dict[str, Any]:
    """Build an immutable, self-identifying level checkpoint."""
    if status not in CHECKPOINT_STATUSES:
        raise ValueError(f"unrecognized checkpoint status: {status!r}")
    document: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "kind": "adaptive_chi_level",
        "manifest_id": manifest_id,
        "task_id": task_id,
        "record_id": record_id,
        "chi": chi,
        "cutoff": cutoff,
        "status": status,
        "valid": valid,
        "elapsed_seconds": elapsed_seconds,
        "peak_rss_bytes": peak_rss_bytes,
    }
    if status == "complete":
        document.update(
            probabilities=list(probabilities) if probabilities is not None else None,
            logical_mask=logical_mask,
            log_joint_syndrome=log_joint_syndrome,
        )
    elif any(value is not None for value in (probabilities, logical_mask, log_joint_syndrome)):
        raise ValueError("non-complete checkpoints cannot contain a decoded class")
    document["checkpoint_id"] = identity_hash(document, "checkpoint_id")
    return validate_level_checkpoint(
        document,
        manifest_id=manifest_id,
        task_id=task_id,
        record_id=record_id,
        chi=chi,
        cutoff=cutoff,
    )


def validate_level_checkpoint(
    document: Mapping[str, Any],
    *,
    manifest_id: str,
    task_id: int,
    record_id: str,
    chi: int,
    cutoff: float,
) -> dict[str, Any]:
    """Validate checkpoint identity, expected ownership, and scientific payload."""
    result = dict(document)
    checkpoint_id = result.get("checkpoint_id")
    if not isinstance(checkpoint_id, str) or checkpoint_id != identity_hash(result, "checkpoint_id"):
        raise ValueError("checkpoint identity mismatch")
    expected = {
        "manifest_id": manifest_id,
        "task_id": task_id,
        "record_id": record_id,
        "chi": chi,
        "cutoff": cutoff,
    }
    if any(result.get(key) != value for key, value in expected.items()):
        raise ValueError("checkpoint identity does not match the expected task")
    if result.get("schema_version") != SCHEMA_VERSION or result.get("kind") != "adaptive_chi_level":
        raise ValueError("checkpoint schema identity mismatch")
    if RECORD_ID_PATTERN.fullmatch(str(manifest_id)) is None:
        raise ValueError("manifest_id must be 64 lowercase hexadecimal characters")
    if RECORD_ID_PATTERN.fullmatch(str(record_id)) is None:
        raise ValueError("record_id must be 64 lowercase hexadecimal characters")
    if isinstance(task_id, bool) or not isinstance(task_id, int) or task_id < 0:
        raise ValueError("task_id must be a nonnegative integer")
    _positive_integer(chi, "chi")
    cutoff_value = _finite_number(cutoff, "cutoff")
    if cutoff_value < 0.0:
        raise ValueError("cutoff must be nonnegative")
    _finite_number(result.get("elapsed_seconds"), "elapsed_seconds")
    if result["elapsed_seconds"] < 0.0:
        raise ValueError("elapsed_seconds must be nonnegative")
    rss = result.get("peak_rss_bytes")
    if isinstance(rss, bool) or not isinstance(rss, int) or rss < 0:
        raise ValueError("peak_rss_bytes must be a nonnegative integer")
    status = result.get("status")
    if status not in CHECKPOINT_STATUSES:
        raise ValueError(f"unrecognized checkpoint status: {status!r}")
    if status == "complete":
        if result.get("valid") is not True:
            raise ValueError("a complete checkpoint must be valid")
        probabilities = _validate_probabilities(result.get("probabilities"))
        _validate_logical_mask(result.get("logical_mask"), probabilities)
        _finite_number(result.get("log_joint_syndrome"), "log_joint_syndrome")
    else:
        if result.get("valid") is not False:
            raise ValueError("a non-complete checkpoint must be invalid")
        forbidden = {"probabilities", "logical_mask", "log_joint_syndrome"} & result.keys()
        if forbidden:
            raise ValueError("an unknown checkpoint cannot contain a decoded class")
    return result
