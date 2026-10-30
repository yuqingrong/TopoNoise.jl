#!/usr/bin/env python3
"""Report a completed, authenticated paired fixed-chi screen without decoding.

Exit zero means the saved benchmark was read successfully. Class disagreement,
full-probability audit failure, and numerical invalidity remain separate counts;
none is silently promoted to a certification or discarded from the denominator.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics
from typing import Any, Mapping

try:
    from . import adaptive_runner, core
    from .inspect_adaptive_ladder import EXACT_AUDIT_POLICY
except ImportError:
    import adaptive_runner
    import core
    from inspect_adaptive_ladder import EXACT_AUDIT_POLICY

HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parents[2]
CHIS = [64, 128]
SEED_BASE = 202610700000


def _integer(value: Any, name: str, minimum: int, maximum: int) -> int:
    if type(value) is not int or not minimum <= value <= maximum:
        raise ValueError(f"{name} must be an integer in {minimum}:{maximum}")
    return value


def _seconds(value: Any, name: str) -> float:
    if (isinstance(value, bool) or not isinstance(value, (int, float))
            or not math.isfinite(value) or value < 0):
        raise ValueError(f"{name} must be finite and nonnegative")
    return float(value)


def _digest(value: Any, name: str) -> str:
    if not isinstance(value, str) or core.RECORD_ID_PATTERN.fullmatch(value) is None:
        raise ValueError(f"{name} must be a lowercase SHA-256 digest")
    return value


def _verify_sources(run: Mapping[str, Any]) -> None:
    names = {str(path.relative_to(REPOSITORY)) for path in (REPOSITORY / "src").rglob("*.jl")}
    names.update({"Project.toml", "examples/tnml_repaired_seed/model.jl",
                  "examples/tnml_repaired_seed/run.jl", "examples/tnml_repaired_seed/dense_frontier.jl"})
    names.update("examples/tnml_repaired_seed/adaptive_chi_pilot/" + name for name in (
        "persistent_worker.jl", "julia_adapter.jl", "fixed_chi_screen.jl",
        "inspect_fixed_chi_screen.py", "core.py", "adaptive_runner.py", "inspect_adaptive_ladder.py",
    ))
    if (REPOSITORY / "Manifest.toml").is_file():
        names.add("Manifest.toml")
    sources = run.get("source_sha256")
    if not isinstance(sources, dict) or set(sources) != names:
        raise ValueError("source hash inventory differs from the fixed-chi screen")
    for name, expected in sources.items():
        _digest(expected, "source hash")
        actual = hashlib.sha256((REPOSITORY / name).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f"source changed after screening started: {name}")


def _validate_run(run: Mapping[str, Any], completion: Mapping[str, Any]) -> int:
    if (type(run.get("schema_version")) is not int or run["schema_version"] != 1
            or run.get("kind") != "standalone_fixed_chi_screen"
            or type(run.get("L")) is not int or run["L"] != 3
            or run.get("production_gate") is not False):
        raise ValueError("run identity is not the standalone L=3 fixed-chi screen")
    index = _integer(run.get("p_index"), "p_index", 0, 5)
    shots = _integer(run.get("shots"), "shots", 1, 10)
    if (run.get("p_key") != core.P_STRATA[index]
            or type(run.get("p")) not in (int, float)
            or run["p"] != float(core.P_STRATA[index])):
        raise ValueError("run probability identity differs from p_index")
    if run.get("chis") != CHIS or any(type(chi) is not int for chi in run["chis"]):
        raise ValueError("screen must contain the fixed chis 64 and 128")
    if type(run.get("seed_base")) is not int or run["seed_base"] != SEED_BASE:
        raise ValueError("run seed base differs from the fresh screening cohort")
    if type(run.get("cutoff")) not in (int, float) or run["cutoff"] != 1e-18:
        raise ValueError("screen cutoff differs from 1e-18")
    _seconds(run.get("startup_seconds"), "startup_seconds")
    runtime = run.get("runtime_identity")
    if not isinstance(runtime, dict):
        raise ValueError("runtime identity is missing")
    for name in ("julia_version", "blas_vendor", "cpu_model"):
        if not isinstance(runtime.get(name), str) or not runtime[name]:
            raise ValueError(f"runtime identity lacks {name}")
    _integer(runtime.get("julia_threads"), "julia_threads", 1, 1)
    _integer(runtime.get("blas_threads"), "blas_threads", 1, 1024)
    if (completion.get("status") != "complete"
            or completion.get("source_identity_verified") is not True
            or completion.get("production_gate") is not False):
        raise ValueError("worker completion is missing or not authenticated")
    for name in ("requested_shots", "completed_shots"):
        if _integer(completion.get(name), name, 1, 10) != shots:
            raise ValueError("worker completion count differs from requested shots")
    _seconds(completion.get("elapsed_seconds"), "completion elapsed_seconds")
    _verify_sources(run)
    return shots


def _validate_result(result: Mapping[str, Any], record_id: str, index: int,
                     *, chi: int | None = None) -> None:
    if result.get("record_id") != record_id:
        raise ValueError(f"shot {index}: result record identity mismatch")
    if _integer(result.get("shot_index"), "result shot_index", 1, 10) != index:
        raise ValueError(f"shot {index}: result shot index mismatch")
    if result.get("cold_start") is not (index == 1):
        raise ValueError(f"shot {index}: result cold-start timing identity mismatch")
    _seconds(result.get("elapsed_seconds"), "elapsed_seconds")
    status = result.get("status")
    if status == "complete":
        if core._eligible_level({**result, "chi": chi or 1})["status"] != "complete":
            raise ValueError(f"shot {index}: malformed complete scientific result")
    elif status == "invalid" and chi is not None:
        if result.get("valid") is not False:
            raise ValueError(f"shot {index}: malformed invalid result")
        if "error" in result and result["error"] != "":
            raise ValueError(f"shot {index}: invalid result contains an error diagnostic")
        if {"probabilities", "logical_mask", "log_joint_syndrome"} & result.keys():
            raise ValueError(f"shot {index}: invalid result contains decoded values")
    else:
        raise ValueError(f"shot {index}: unexpected worker status {status!r}: {result.get('error', '')}")
    if chi is not None:
        if type(result.get("chi")) is not int or result["chi"] != chi:
            raise ValueError(f"shot {index}: chi identity mismatch")
        for name in ("decode_seconds", "decoder_build_seconds"):
            _seconds(result.get(name), name)
        if type(result.get("decoder_cache_hit")) is not bool:
            raise ValueError(f"shot {index}: decoder cache evidence is missing")


def _load_history(directory: Path, run: Mapping[str, Any], index: int) -> dict[str, Any]:
    folder = directory / "shots" / f"shot_{index:04d}"
    metadata = adaptive_runner._read_toml(folder / "metadata.toml")
    if metadata.get("record_authenticated") is not True:
        raise ValueError(f"shot {index}: Julia record authentication did not complete")
    if _integer(metadata.get("shot_index"), "metadata shot_index", 1, 10) != index:
        raise ValueError(f"shot {index}: metadata shot identity mismatch")
    seed = SEED_BASE + 100 * run["p_index"] + index - 1
    if type(metadata.get("seed")) is not int or metadata["seed"] != seed:
        raise ValueError(f"shot {index}: saved seed differs from the fresh cohort")
    if metadata.get("cold_start") is not (index == 1):
        raise ValueError(f"shot {index}: metadata cold-start flag differs")
    if metadata.get("chi_order") != (CHIS if index % 2 else list(reversed(CHIS))):
        raise ValueError(f"shot {index}: paired chi order differs")
    record_path = folder / "record.toml"
    record_hash = _digest(metadata.get("record_file_sha256"), "record content hash")
    if hashlib.sha256(record_path.read_bytes()).hexdigest() != record_hash:
        raise ValueError(f"shot {index}: record content hash changed after authentication")
    record = adaptive_runner._read_toml(record_path)
    record_id = _digest(metadata.get("record_id"), "record identity")
    if record.get("record_id") != record_id:
        raise ValueError(f"shot {index}: saved record identity mismatch")
    if _integer(record.get("shot_index"), "record shot_index", 1, 10) != index:
        raise ValueError(f"shot {index}: saved record shot index mismatch")
    truth = _integer(record.get("true_mask"), "true_mask", 0, 3)
    exact = adaptive_runner._read_toml(folder / "exact.toml")
    _validate_result(exact, record_id, index)
    levels = {}
    for chi in CHIS:
        level = adaptive_runner._read_toml(folder / f"chi_{chi:04d}.toml")
        _validate_result(level, record_id, index, chi=chi)
        levels[str(chi)] = level
    return {"shot_index": index, "seed": seed, "record_id": record_id,
            "true_mask": truth, "exact": exact, "levels": levels}


def inspect_directory(directory: Path) -> dict[str, Any]:
    directory = Path(directory).resolve()
    if (directory / "failure.toml").exists():
        raise ValueError("worker saved failure.toml; screening did not complete successfully")
    run = adaptive_runner._read_toml(directory / "run.toml")
    completion = adaptive_runner._read_toml(directory / "completion.toml")
    shots = _validate_run(run, completion)
    histories = [_load_history(directory, run, index) for index in range(1, shots + 1)]
    by_chi = {}
    for chi in CHIS:
        levels = [history["levels"][str(chi)] for history in histories]
        audits = [adaptive_runner._exact_audit(
            {**level, "status": "accepted"} if level["status"] == "complete"
            else {"status": "unknown"}, history["exact"], EXACT_AUDIT_POLICY,
        ) for history, level in zip(histories, levels)]
        valid = [index for index, level in enumerate(levels) if level["status"] == "complete"]
        unknown = shots - len(valid)
        failures = sum(levels[index]["logical_mask"] != histories[index]["true_mask"] for index in valid)
        warm = [level["decode_seconds"] for level in levels[1:]]
        by_chi[str(chi)] = {
            "requested_shots": shots, "valid_shots": len(valid), "unknown_shots": unknown,
            "class_disagreements": sum(not audits[index]["ml_tie_aware_match"] for index in valid),
            "strict_class_disagreements": sum(levels[index]["logical_mask"] != histories[index]["exact"]["logical_mask"] for index in valid),
            "full_audit_failures": sum(not audits[index]["passed"] for index in valid),
            "posterior_linf_max": max((audits[index]["posterior_linf"] for index in valid), default=None),
            "log_joint_error_max": max((audits[index]["log_joint_error"] for index in valid), default=None),
            "logical_failures": failures,
            "failure_interval": [failures / shots, (failures + unknown) / shots],
            "warm_median_decode_seconds": statistics.median(warm) if warm else None,
            "warm_timed_shots": len(warm),
            "total_decode_seconds": sum(level["decode_seconds"] for level in levels),
            "total_elapsed_seconds": sum(level["elapsed_seconds"] for level in levels),
            "total_decoder_build_seconds": sum(level["decoder_build_seconds"] for level in levels),
            "decoder_cache_hits": sum(level["decoder_cache_hit"] for level in levels),
            "decoder_cache_misses": sum(not level["decoder_cache_hit"] for level in levels),
            "first_history_timing": {key: levels[0][key] for key in (
                "elapsed_seconds", "decode_seconds", "decoder_build_seconds", "decoder_cache_hit")},
        }
    both_valid = [history for history in histories if all(level["status"] == "complete" for level in history["levels"].values())]
    return {
        "schema_version": 1, "kind": "standalone_fixed_chi_screen_report", "status": "complete",
        "L": run["L"], "p": run["p"], "p_key": run["p_key"], "p_index": run["p_index"],
        "requested_shots": shots, "completed_shots": shots, "chis": CHIS,
        "seed_base": SEED_BASE, "cutoff": run["cutoff"], "production_gate": False,
        "diagnostic_directory": str(directory), "by_chi": by_chi,
        "paired": {"requested_shots": shots, "both_valid_shots": len(both_valid),
                   "unknown_shots": shots - len(both_valid),
                   "decoder_class_disagreements": sum(history["levels"]["64"]["logical_mask"] != history["levels"]["128"]["logical_mask"] for history in both_valid)},
        "exact_reference": {"valid_shots": shots,
                            "logical_failures": sum(history["exact"]["logical_mask"] != history["true_mask"] for history in histories),
                            "total_elapsed_seconds": sum(history["exact"]["elapsed_seconds"] for history in histories)},
        "exact_audit_policy": EXACT_AUDIT_POLICY, "runtime_identity": run["runtime_identity"],
        "source_sha256": run["source_sha256"], "startup_seconds": run["startup_seconds"],
        "startup_timing_scope": run.get("startup_timing_scope", "unspecified"),
        "worker_elapsed_seconds": completion["elapsed_seconds"],
        "timing_scope": "warm medians exclude history 1; totals include all histories, including invalid attempts",
        "failure_interval_scope": "unknown-outcome bounds on the requested histories, not statistical confidence intervals",
        "records": [{"shot_index": history["shot_index"], "seed": history["seed"],
                     "record_id": history["record_id"], "true_mask": history["true_mask"]} for history in histories],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    try:
        report = inspect_directory(args.directory)
        path = args.directory.resolve() / "report.json"
        serialized = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
        with path.open("x", encoding="utf-8") as stream:
            stream.write(serialized)
    except (OSError, ValueError) as error:
        parser.exit(2, f"fixed-chi report failed: {error}\n")
    compact = {key: report[key] for key in (
        "status", "L", "p_key", "requested_shots", "completed_shots", "by_chi", "paired", "exact_reference", "production_gate")}
    compact["report_path"] = str(path)
    print(json.dumps(compact, sort_keys=True, allow_nan=False), flush=True)
    for chi, result in report["by_chi"].items():
        print(f"SCREEN p={report['p_key']} chi={chi} requested={result['requested_shots']} "
              f"valid={result['valid_shots']} unknown={result['unknown_shots']} "
              f"class_disagreements={result['class_disagreements']} "
              f"full_audit_failures={result['full_audit_failures']} "
              f"warm_median_decode_seconds={result['warm_median_decode_seconds']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
