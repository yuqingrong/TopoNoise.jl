#!/usr/bin/env python3
"""Collect authenticated adaptive checkpoints and enforce the L3-to-L5 gate."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import statistics
from typing import Any, Mapping

if __package__:
    from .. import adaptive_runner, core, execution_manifest, manifest, task_runner
    from . import control, gate
else:
    HERE = Path(__file__).resolve().parent
    sys.path.insert(0, str(HERE.parent))
    import adaptive_runner, core, execution_manifest, manifest, task_runner
    import control, gate


REPORT_NAME = "report.json"
GATE_NAME = gate.GATE_NAME


class _ReportExecutor:
    executor_mode = "persistent"

    def run_tn(self, *_args, **_kwargs):
        raise RuntimeError("collector never executes a missing TN checkpoint")

    def run_exact(self, *_args, **_kwargs):
        raise RuntimeError("collector never executes a missing exact checkpoint")


def _write_once(path: Path, document: Mapping[str, Any], label: str) -> Path:
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


def build_l3_gate(
    *,
    sampling_manifest_id: str,
    primary_execution_manifest_id: str,
    replay_execution_manifest_id: str,
    report_id: str,
    passed: bool,
    criteria: Mapping[str, Any],
    source_bundle_id: str,
    corpus_id: str,
    corpus_provenance_id: str,
    environment_profile_id: str,
) -> dict[str, Any]:
    return gate.build_l3_gate(
        sampling_manifest_id=sampling_manifest_id,
        primary_execution_manifest_id=primary_execution_manifest_id,
        replay_execution_manifest_id=replay_execution_manifest_id,
        report_id=report_id, passed=passed, criteria=criteria,
        source_bundle_id=source_bundle_id, corpus_id=corpus_id,
        corpus_provenance_id=corpus_provenance_id,
        environment_profile_id=environment_profile_id,
    )


def validate_l3_gate(document: Mapping[str, Any]) -> dict[str, Any]:
    return gate.validate_l3_gate(document)


def write_l3_gate(study_root: os.PathLike[str] | str, document: Mapping[str, Any]) -> Path:
    return gate.write_l3_gate(study_root, document)


def require_l3_gate(study_root: os.PathLike[str] | str) -> dict[str, Any]:
    return gate.require_l3_gate(study_root)


def _decision_path(stage_root: Path, sampling: Mapping[str, Any], task, record) -> Path:
    plan = adaptive_runner.plan_record(
        stage_root, sampling, task_id=task["task_id"], record_id=record["record_id"]
    )
    return Path(plan["record_directory"]) / "decision.json"


def _collect_stage(stage_root: Path, sampling: Mapping[str, Any]) -> dict[str, Any]:
    execution = execution_manifest.load_execution_manifest(stage_root, sampling)
    task_runner.verify_l5_gate_binding(stage_root, execution, require_hpc=True)
    profile = control.load_environment_profile(stage_root.parent)
    shard_root = (
        stage_root.parent / control.STAGE_ROOTS[("l3_exact", False)]
        if stage_root.name == control.STAGE_ROOTS[("l3_exact", True)]
        else stage_root
    )
    control.load_corpus_provenance(
        stage_root,
        sampling,
        profile=profile,
        execution=execution,
        shard_directory=shard_root,
    )
    execution_manifest.verify_source_identity(
        Path(__file__).resolve().parents[4], execution["source_identity"]
    )
    rows = []
    for task in execution["tasks"]:
        for record in task["records"]:
            path = _decision_path(stage_root, sampling, task, record)
            if not path.is_file():
                rows.append({"task_id": task["task_id"], "record_id": record["record_id"],
                             "p_key": record["p_key"], "status": "missing"})
                continue
            try:
                decision = adaptive_runner.run_record(
                    stage_root, sampling, task_id=task["task_id"],
                    record_id=record["record_id"], executor=_ReportExecutor(),
                )
            except Exception as error:
                rows.append({"task_id": task["task_id"], "record_id": record["record_id"],
                             "p_key": record["p_key"], "status": "invalid",
                             "error": f"{type(error).__name__}: {error}"})
            else:
                rows.append({"task_id": task["task_id"], "record_id": record["record_id"],
                             "p_key": record["p_key"], "status": "complete",
                             "decision": decision, "record_directory": str(path.parent)})
    return {"execution": execution, "rows": rows}


def _level_projection(directory: Path, chis: list[int]) -> dict[int, dict[str, Any]]:
    result = {}
    for chi in chis:
        path = directory / f"chi_{chi:04d}.json"
        if not path.is_file():
            continue
        level = json.loads(path.read_text(encoding="utf-8"))
        result[chi] = {
            "status": level.get("status"), "valid": level.get("valid"),
            "probabilities": level.get("probabilities"),
            "logical_mask": level.get("logical_mask"),
            "log_joint_syndrome": level.get("log_joint_syndrome"),
        }
    return result


def _compare_replay(primary: Mapping[str, Any], replay: Mapping[str, Any], policy) -> dict[str, Any]:
    p_decision, r_decision = primary["decision"], replay["decision"]
    fields = ("final_status", "unknown_reason", "tested_chi", "logical_mask")
    decision_match = (
        all(p_decision.get(name) == r_decision.get(name) for name in fields)
        and p_decision.get("adaptive_decision", {}).get("status")
        == r_decision.get("adaptive_decision", {}).get("status")
        and p_decision.get("adaptive_decision", {}).get("accepted_chi")
        == r_decision.get("adaptive_decision", {}).get("accepted_chi")
    )
    p_levels = _level_projection(Path(primary["record_directory"]), p_decision["tested_chi"])
    r_levels = _level_projection(Path(replay["record_directory"]), r_decision["tested_chi"])
    level_keys_match = set(p_levels) == set(r_levels)
    maximum_posterior = 0.0
    maximum_log = 0.0
    level_status_match = level_keys_match
    if level_keys_match:
        for chi in p_levels:
            left, right = p_levels[chi], r_levels[chi]
            if (left["status"], left["valid"], left["logical_mask"]) != (
                right["status"], right["valid"], right["logical_mask"]
            ):
                level_status_match = False
            if left["status"] == "complete" and right["status"] == "complete":
                maximum_posterior = max(
                    maximum_posterior,
                    max(abs(a - b) for a, b in zip(left["probabilities"], right["probabilities"])),
                )
                maximum_log = max(maximum_log, abs(
                    left["log_joint_syndrome"] - right["log_joint_syndrome"]
                ))
    passed = (
        decision_match and level_status_match
        and maximum_posterior <= policy["posterior_linf_gate"]
        and maximum_log <= policy["log_joint_gate"]
        and r_decision.get("exact_audit", {}).get("status") == "complete"
        and r_decision.get("exact_audit", {}).get("false_accept") is not True
    )
    return {"record_id": p_decision["record_id"], "p_key": primary["p_key"],
            "decision_match": decision_match, "level_status_match": level_status_match,
            "posterior_linf": maximum_posterior, "log_joint_error": maximum_log,
            "passed": passed}


def _checkpoint_summary(rows: list[Mapping[str, Any]]) -> dict[str, Any]:
    chi_stops = Counter()
    elapsed_by_chi = defaultdict(list)
    decoder_build_by_chi = defaultdict(list)
    decode_by_chi = defaultdict(list)
    cache_hits = Counter()
    startup_by_record = []
    total_elapsed = []
    maximum_rss = 0
    for row in rows:
        decision = row["decision"]
        tested = decision.get("tested_chi", [])
        chi_stops[str(tested[-1]) if tested else "none"] += 1
        record_total = 0.0
        directory = Path(row["record_directory"])
        record_startup = 0.0
        exact_path = directory / "exact.json"
        if exact_path.is_file():
            exact = json.loads(exact_path.read_text(encoding="utf-8"))
            record_startup += float(exact.get("worker_startup_seconds", 0.0))
        for chi in tested:
            level = json.loads((directory / f"chi_{chi:04d}.json").read_text(encoding="utf-8"))
            elapsed = float(level.get("elapsed_seconds", 0.0))
            elapsed_by_chi[str(chi)].append(elapsed)
            decoder_build_by_chi[str(chi)].append(float(level.get("decoder_build_seconds", 0.0)))
            decode_by_chi[str(chi)].append(float(level.get("decode_seconds", 0.0)))
            cache_hits["hit" if level.get("decoder_cache_hit") is True else "miss"] += 1
            maximum_rss = max(maximum_rss, int(level.get("peak_rss_bytes", 0)))
            record_startup += float(level.get("worker_startup_seconds", 0.0))
            record_total += elapsed
        startup_by_record.append(record_startup)
        record_total += record_startup
        total_elapsed.append(record_total)
    def summarize(values):
        return None if not values else {
            "count": len(values), "median": statistics.median(values),
            "maximum": max(values), "total": sum(values),
        }
    return {
        "chi_stop_histogram": dict(sorted(chi_stops.items())),
        "elapsed_by_chi": {key: summarize(values) for key, values in sorted(elapsed_by_chi.items())},
        "decoder_build_by_chi": {key: summarize(values) for key, values in sorted(decoder_build_by_chi.items())},
        "decode_by_chi": {key: summarize(values) for key, values in sorted(decode_by_chi.items())},
        "decoder_cache": dict(sorted(cache_hits.items())),
        "worker_startup_seconds": summarize(startup_by_record),
        "record_total_elapsed_seconds": total_elapsed,
        "record_total_definition": (
            "sum_tn_request_elapsed_plus_worker_process_and_context_startup; "
            "exact_reference_elapsed_excluded"
        ),
        "maximum_worker_lifetime_rss_bytes": maximum_rss,
        "rss_scope": "worker_lifetime_high_water",
    }


def _binomial_cdf(k: int, n: int, p: float) -> float:
    if p <= 0.0:
        return 1.0
    if p >= 1.0:
        return 1.0 if k >= n else 0.0
    logs = [
        math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
        + i * math.log(p) + (n - i) * math.log1p(-p)
        for i in range(k + 1)
    ]
    maximum = max(logs)
    return math.exp(maximum) * sum(math.exp(value - maximum) for value in logs)


def clopper_pearson_upper(events: int, trials: int, confidence: float = 0.95) -> float | None:
    if trials == 0:
        return None
    if not 0 <= events <= trials:
        raise ValueError("Clopper-Pearson counts are invalid")
    if events == trials:
        return 1.0
    alpha = 1.0 - confidence
    low, high = events / trials, 1.0
    for _ in range(80):
        middle = (low + high) / 2.0
        if _binomial_cdf(events, trials, middle) > alpha:
            low = middle
        else:
            high = middle
    return (low + high) / 2.0


def _all_forced_anchors_complete(
    anchors: list[Mapping[str, Any]], *, final_chi: int
) -> bool:
    return len(anchors) == 24 and all(
        _forced_anchor_complete(item, final_chi=final_chi)
        and item.get("adaptive_decision", {}).get("anchor_class_changed") is False
        for item in anchors
    )


def _forced_anchor_complete(item: Mapping[str, Any], *, final_chi: int) -> bool:
    return (
        item.get("final_status") == "accepted"
        and item.get("adaptive_decision", {}).get("anchor_completed") is True
        and item.get("adaptive_decision", {}).get("anchor_final_chi") == final_chi
        and item.get("adaptive_decision", {}).get("anchor_class_changed") in (True, False)
    )


def collect_l3_report(study_root: os.PathLike[str] | str, sampling: Mapping[str, Any]) -> dict[str, Any]:
    root = Path(study_root).expanduser().resolve()
    control.verify_live_environment_profile(
        root, repository_root=Path(__file__).resolve().parents[4], julia="julia"
    )
    primary = _collect_stage(control.stage_root(root, "l3_exact"), sampling)
    replay = _collect_stage(control.stage_root(root, "l3_exact", replay=True), sampling)
    p_complete = [row for row in primary["rows"] if row["status"] == "complete"]
    r_complete = [row for row in replay["rows"] if row["status"] == "complete"]
    primary_by_id = {row["record_id"]: row for row in p_complete}
    replay_comparisons = []
    for row in r_complete:
        other = primary_by_id.get(row["record_id"])
        if other is not None:
            replay_comparisons.append(_compare_replay(
                other, row, replay["execution"]["replay_comparison"]
            ))
    decisions = [row["decision"] for row in p_complete]
    exact_complete = sum(item.get("exact_audit", {}).get("status") == "complete" for item in decisions)
    false_accepts = sum(item.get("exact_audit", {}).get("false_accept") is True for item in decisions)
    accepted = sum(item.get("final_status") == "accepted" for item in decisions)
    unknown_reasons = Counter(
        item.get("unknown_reason", "unspecified")
        for item in decisions if item.get("final_status") != "accepted"
    )
    unknown_reasons["missing_decision"] += sum(
        row["status"] == "missing" for row in primary["rows"]
    )
    unknown_reasons["invalid_decision"] += sum(
        row["status"] == "invalid" for row in primary["rows"]
    )
    unknown_reasons += Counter()
    per_p = {}
    for p_key in manifest.P_VALUES:
        subset = [item for item in decisions if next(
            row for row in p_complete if row["record_id"] == item["record_id"]
        )["p_key"] == p_key]
        evaluated = [item for item in subset if item.get("exact_audit", {}).get("evaluated") is True]
        failures = sum(item["exact_audit"].get("false_accept") is True for item in evaluated)
        per_p[p_key] = {"requested": 100, "decisions": len(subset),
                        "accepted_with_exact": len(evaluated), "false_accepts": failures,
                        "false_accept_cp95_upper": clopper_pearson_upper(failures, len(evaluated))}
    criteria = {
        "all_600_primary_decisions": len(p_complete) == 600,
        "all_600_exact_references_complete": exact_complete == 600,
        "primary_exact_gate": false_accepts == 0,
        "all_30_blas1_replays": len(r_complete) == 30,
        "blas1_replay_gate": len(replay_comparisons) == 30 and all(
            item["passed"] for item in replay_comparisons
        ),
        "source_bundle_match": (
            primary["execution"]["source_bundle_id"] == replay["execution"]["source_bundle_id"]
        ),
    }
    passed = all(criteria.values())
    terminal_complete = len(p_complete) == 600 and len(r_complete) == 30
    document = {
        "schema_version": 1, "kind": "adaptive_chi_l3_calibration_report",
        "sampling_manifest_id": sampling["manifest_id"],
        "primary_execution_manifest_id": primary["execution"]["execution_manifest_id"],
        "replay_execution_manifest_id": replay["execution"]["execution_manifest_id"],
        "source_bundle_id": primary["execution"]["source_bundle_id"],
        "corpus_id": primary["execution"]["corpus_id"],
        "corpus_provenance_id": primary["execution"]["hpc_provenance"]["corpus_provenance_id"],
        "environment_profile_id": primary["execution"]["hpc_provenance"]["environment_profile_id"],
        "requested_primary": 600, "completed_primary_decisions": len(p_complete),
        "accepted_primary": accepted, "unknown_primary": 600 - accepted,
        "unknown_reasons": dict(sorted(unknown_reasons.items())),
        "exact_complete": exact_complete, "false_accepts": false_accepts,
        "false_accept_cp95_upper": clopper_pearson_upper(false_accepts, accepted),
        "per_probability": per_p, "replay_comparisons": replay_comparisons,
        "primary_checkpoint_summary": _checkpoint_summary(p_complete),
        "replay_checkpoint_summary": _checkpoint_summary(r_complete),
        "criteria": criteria, "terminal_complete": terminal_complete, "passed": passed,
        "calibration_scope": (
            "conditional_false_acceptance_only; no minimum L3 acceptance coverage was "
            "preregistered; L5 feasibility separately requires early acceptance"
        ),
        "authorization_if_passed": "l5_feasibility_smoke_only",
        "report_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    document["report_id"] = core.identity_hash(document, "report_id")
    return document


def collect_l5_report(study_root: os.PathLike[str] | str, sampling: Mapping[str, Any]) -> dict[str, Any]:
    root = Path(study_root).expanduser().resolve()
    control.verify_live_environment_profile(
        root, repository_root=Path(__file__).resolve().parents[4], julia="julia"
    )
    require_l3_gate(root)
    collected = _collect_stage(control.stage_root(root, "l5_feasibility"), sampling)
    complete = [row for row in collected["rows"] if row["status"] == "complete"]
    decisions = [row["decision"] for row in complete]
    accepted = [item for item in decisions if item.get("final_status") == "accepted"]
    nonanchors = [item for item in decisions if item.get("forced_anchor") is not True]
    early = [item for item in nonanchors if item.get("adaptive_decision", {}).get("accepted_chi", 10**9) <= 128]
    anchors = [item for item in decisions if item.get("forced_anchor") is True]
    final_chi = max(collected["execution"]["chi_ladder"])
    anchors_complete = _all_forced_anchors_complete(anchors, final_chi=final_chi)
    completed_anchors = [
        item for item in anchors if _forced_anchor_complete(item, final_chi=final_chi)
    ]
    anchor_changes = sum(
        item.get("adaptive_decision", {}).get("anchor_class_changed") is True
        for item in completed_anchors
    )
    per_p = {}
    for p_key in manifest.P_VALUES:
        rows = [row for row in complete if row["p_key"] == p_key]
        subset = [row["decision"] for row in rows]
        subset_anchors = [item for item in subset if item.get("forced_anchor") is True]
        subset_completed_anchors = [
            item for item in subset_anchors
            if _forced_anchor_complete(item, final_chi=final_chi)
        ]
        subset_changes = sum(
            item.get("adaptive_decision", {}).get("anchor_class_changed") is True
            for item in subset_completed_anchors
        )
        subset_nonanchors = [item for item in subset if item.get("forced_anchor") is not True]
        subset_early = sum(
            item.get("final_status") == "accepted"
            and item.get("adaptive_decision", {}).get("accepted_chi", 10**9) <= 128
            for item in subset_nonanchors
        )
        per_p[p_key] = {
            "requested": 10,
            "decisions": len(subset),
            "accepted": sum(item.get("final_status") == "accepted" for item in subset),
            "unknown": 10 - sum(item.get("final_status") == "accepted" for item in subset),
            "nonanchor_early_accepts": subset_early,
            "anchors_requested": 4,
            "anchors_completed": len(subset_completed_anchors),
            "anchor_class_changes": subset_changes,
            "anchor_change_cp95_upper": clopper_pearson_upper(
                subset_changes, len(subset_completed_anchors)
            ),
        }
    summary = _checkpoint_summary(complete)
    elapsed = sorted(summary["record_total_elapsed_seconds"])
    p90 = None if not elapsed else elapsed[min(len(elapsed) - 1, math.ceil(0.9 * len(elapsed)) - 1)]
    criteria = {
        "all_60_decisions": len(complete) == 60,
        "nonanchor_accept_by_128_at_least_80_percent": (
            bool(nonanchors) and len(early) / len(nonanchors) >= 0.8
        ),
        "no_forced_anchor_class_change": anchors_complete,
        "median_wall_below_60_seconds": bool(elapsed) and statistics.median(elapsed) < 60.0,
        "p90_wall_below_300_seconds": p90 is not None and p90 < 300.0,
        "observed_worker_high_water_below_16_gib": (
            summary["maximum_worker_lifetime_rss_bytes"] < 16 * 2**30
        ),
    }
    document = {
        "schema_version": 1, "kind": "adaptive_chi_l5_feasibility_report",
        "classification": "feasibility_smoke_not_threshold_data",
        "sampling_manifest_id": sampling["manifest_id"],
        "execution_manifest_id": collected["execution"]["execution_manifest_id"],
        "requested": 60, "completed_decisions": len(complete),
        "accepted": len(accepted), "unknown": 60 - len(accepted),
        "anchors_completed": len(completed_anchors),
        "anchor_class_changes": anchor_changes,
        "anchor_change_cp95_upper": clopper_pearson_upper(
            anchor_changes, len(completed_anchors)
        ),
        "per_probability": per_p,
        "checkpoint_summary": summary,
        "p90_record_elapsed_seconds": p90, "criteria": criteria,
        "terminal_complete": len(complete) == 60,
        "passed": all(criteria.values()),
        "authorization_if_passed": "larger_preregistered_anchor_validation_only",
        "report_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    document["report_id"] = core.identity_hash(document, "report_id")
    return document


def progress_report(study_root: os.PathLike[str] | str, sampling: Mapping[str, Any]) -> dict[str, Any]:
    root = Path(study_root).expanduser().resolve()
    targets = (("l3-primary", control.stage_root(root, "l3_exact")),
               ("l3-replay", control.stage_root(root, "l3_exact", replay=True)),
               ("l5-primary", control.stage_root(root, "l5_feasibility")))
    result = {}
    for name, path in targets:
        manifest_path = path / execution_manifest.EXECUTION_MANIFEST_NAME
        if not manifest_path.is_file():
            result[name] = {"status": "not_frozen"}
            continue
        execution = execution_manifest.load_execution_manifest(path, sampling)
        expected = sum(len(task["records"]) for task in execution["tasks"])
        decisions = len(list((path / "outputs").rglob("decision.json"))) if (path / "outputs").exists() else 0
        levels = len(list((path / "outputs").rglob("chi_*.json"))) if (path / "outputs").exists() else 0
        result[name] = {"status": "complete" if decisions == expected else "running_or_incomplete",
                        "expected_records": expected, "decision_files": decisions,
                        "level_checkpoint_files": levels}
    return {"kind": "adaptive_chi_progress", "stages": result}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-root", required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--progress", action="store_true")
    mode.add_argument("--collect-l3", action="store_true")
    mode.add_argument("--collect-l5", action="store_true")
    mode.add_argument("--require-l3-gate", action="store_true")
    args = parser.parse_args(argv)
    root = Path(args.study_root).expanduser().resolve()
    sampling = manifest.load_manifest(root)
    if args.require_l3_gate:
        document = require_l3_gate(root)
    elif args.progress:
        document = progress_report(root, sampling)
    elif args.collect_l3:
        document = collect_l3_report(root, sampling)
        if document["terminal_complete"]:
            _write_once(control.stage_root(root, "l3_exact") / REPORT_NAME,
                        document, "L3 report")
            gate_document = build_l3_gate(
                sampling_manifest_id=sampling["manifest_id"],
                primary_execution_manifest_id=document["primary_execution_manifest_id"],
                replay_execution_manifest_id=document["replay_execution_manifest_id"],
                report_id=document["report_id"], passed=document["passed"],
                criteria=document["criteria"],
                source_bundle_id=document["source_bundle_id"],
                corpus_id=document["corpus_id"],
                corpus_provenance_id=document["corpus_provenance_id"],
                environment_profile_id=document["environment_profile_id"],
            )
            write_l3_gate(root, gate_document)
            document = {**document, "gate_id": gate_document["gate_id"]}
        else:
            document = {**document, "collection_status": "incomplete_no_gate_written"}
    else:
        document = collect_l5_report(root, sampling)
        if document["terminal_complete"]:
            _write_once(control.stage_root(root, "l5_feasibility") / REPORT_NAME,
                        document, "L5 report")
        else:
            document = {**document, "collection_status": "incomplete_report_not_frozen"}
    print(json.dumps(document, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
