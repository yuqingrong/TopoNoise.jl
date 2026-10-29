"""Immutable L3 calibration gate used by every L5 execution wrapper."""

from __future__ import annotations

import json
import os
from pathlib import Path
import re
from typing import Any, Mapping

try:
    from .. import core, execution_manifest, manifest
except ImportError:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    import core, execution_manifest, manifest


GATE_NAME = "l3-gate.json"
REPORT_NAME = "report.json"
PROFILE_NAME = "environment-profile.json"
PRIMARY_ROOT_NAME = "l3-primary"
REPLAY_ROOT_NAME = "l3-blas1-replay"
HEX64 = re.compile(r"^[0-9a-f]{64}$")
REQUIRED_L3_CRITERIA = frozenset(
    {
        "all_600_primary_decisions",
        "all_600_exact_references_complete",
        "primary_exact_gate",
        "all_30_blas1_replays",
        "blas1_replay_gate",
        "source_bundle_match",
    }
)


def _sha256_identity(value: Any, name: str) -> str:
    if not isinstance(value, str) or HEX64.fullmatch(value) is None:
        raise ValueError(f"{name} must be a SHA-256 identity")
    return value


def _criteria(value: Any) -> dict[str, bool]:
    if not isinstance(value, Mapping) or set(value) != REQUIRED_L3_CRITERIA:
        raise ValueError("L3 gate criteria do not match the frozen required set")
    result = dict(value)
    if any(type(item) is not bool for item in result.values()):
        raise ValueError("L3 gate criteria must be Boolean")
    return result


def build_l3_gate(*, sampling_manifest_id: str, primary_execution_manifest_id: str,
                  replay_execution_manifest_id: str, report_id: str, passed: bool,
                  criteria: Mapping[str, Any], source_bundle_id: str,
                  corpus_id: str, corpus_provenance_id: str,
                  environment_profile_id: str) -> dict[str, Any]:
    criterion_values = _criteria(criteria)
    derived_passed = all(criterion_values.values())
    if type(passed) is not bool or passed != derived_passed:
        raise ValueError("L3 gate passed flag differs from its frozen criteria")
    document = {
        "schema_version": 1, "kind": "adaptive_chi_l3_gate",
        "sampling_manifest_id": _sha256_identity(sampling_manifest_id, "sampling_manifest_id"),
        "primary_execution_manifest_id": _sha256_identity(
            primary_execution_manifest_id, "primary_execution_manifest_id"
        ),
        "replay_execution_manifest_id": _sha256_identity(
            replay_execution_manifest_id, "replay_execution_manifest_id"
        ),
        "report_id": _sha256_identity(report_id, "report_id"),
        "source_bundle_id": _sha256_identity(source_bundle_id, "source_bundle_id"),
        "corpus_id": _sha256_identity(corpus_id, "corpus_id"),
        "corpus_provenance_id": _sha256_identity(
            corpus_provenance_id, "corpus_provenance_id"
        ),
        "environment_profile_id": _sha256_identity(
            environment_profile_id, "environment_profile_id"
        ),
        "passed": passed,
        "authorization": "l5_feasibility_smoke_only" if passed else "none",
        "criteria": criterion_values,
    }
    document["gate_id"] = core.identity_hash(document, "gate_id")
    return document


def validate_l3_gate(document: Mapping[str, Any]) -> dict[str, Any]:
    result = dict(document)
    required = {
        "schema_version", "kind", "sampling_manifest_id",
        "primary_execution_manifest_id", "replay_execution_manifest_id",
        "report_id", "source_bundle_id", "corpus_id", "corpus_provenance_id",
        "environment_profile_id", "passed", "authorization", "criteria", "gate_id",
    }
    if set(result) != required:
        raise ValueError("L3 gate fields do not match the frozen schema")
    if result.get("gate_id") != core.identity_hash(result, "gate_id"):
        raise ValueError("L3 gate identity mismatch")
    if result.get("kind") != "adaptive_chi_l3_gate" or result.get("schema_version") != 1:
        raise ValueError("L3 gate schema mismatch")
    for name in (
        "sampling_manifest_id", "primary_execution_manifest_id",
        "replay_execution_manifest_id", "report_id", "source_bundle_id",
        "corpus_id", "corpus_provenance_id", "environment_profile_id", "gate_id",
    ):
        _sha256_identity(result.get(name), name)
    criteria = _criteria(result.get("criteria"))
    derived_passed = all(criteria.values())
    if type(result.get("passed")) is not bool or result["passed"] != derived_passed:
        raise ValueError("L3 gate passed flag differs from its frozen criteria")
    if result.get("authorization") != (
        "l5_feasibility_smoke_only" if derived_passed else "none"
    ):
        raise ValueError("L3 gate authorization is inconsistent")
    return result


def _load_local_evidence(root: Path) -> tuple[
    dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any]
]:
    sampling = manifest.load_manifest(root)
    profile_path = root / PROFILE_NAME
    if not profile_path.is_file():
        raise FileNotFoundError("local HPC environment profile is missing")
    profile = json.loads(profile_path.read_text(encoding="utf-8"))
    if profile.get("environment_profile_id") != core.identity_hash(
        profile, "environment_profile_id"
    ):
        raise ValueError("local HPC environment profile identity mismatch")
    source = execution_manifest._validate_source_identity(profile.get("source_identity"))
    if profile.get("source_bundle_id") != core.identity_hash(source, "source_bundle_id"):
        raise ValueError("local HPC environment profile source identity mismatch")
    repository_root = Path(__file__).resolve().parents[4]
    execution_manifest.verify_source_identity(repository_root, source)
    primary = execution_manifest.load_execution_manifest(root / PRIMARY_ROOT_NAME, sampling)
    replay = execution_manifest.load_execution_manifest(root / REPLAY_ROOT_NAME, sampling)
    try:
        if __package__:
            from . import control
        else:
            import control
        control.load_corpus_provenance(
            root / PRIMARY_ROOT_NAME,
            sampling,
            profile=profile,
            execution=primary,
        )
        control.load_corpus_provenance(
            root / REPLAY_ROOT_NAME,
            sampling,
            profile=profile,
            execution=replay,
            shard_directory=root / PRIMARY_ROOT_NAME,
        )
    except (OSError, ValueError) as error:
        raise ValueError(f"local L3 corpus provenance is invalid: {error}") from error
    report_path = root / PRIMARY_ROOT_NAME / REPORT_NAME
    if not report_path.is_file():
        raise FileNotFoundError("local L3 calibration report is missing")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report.get("report_id") != core.identity_hash(report, "report_id"):
        raise ValueError("local L3 calibration report identity mismatch")
    if (
        report.get("schema_version") != 1
        or report.get("kind") != "adaptive_chi_l3_calibration_report"
        or report.get("terminal_complete") is not True
    ):
        raise ValueError("local L3 calibration report is not terminal")
    criteria = _criteria(report.get("criteria"))
    if report.get("passed") is not all(criteria.values()):
        raise ValueError("local L3 calibration report pass flag is inconsistent")
    report_source_path = (
        "examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/report.py"
    )
    if report.get("report_source_sha256") != source["files"].get(report_source_path):
        raise ValueError("local L3 report was not produced by the frozen collector source")
    return sampling, profile, primary, replay, report


def write_l3_gate(study_root: os.PathLike[str] | str, document: Mapping[str, Any]) -> Path:
    gate = validate_l3_gate(document)
    path = Path(study_root).expanduser().resolve() / GATE_NAME
    serialized = json.dumps(gate, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if path.exists():
        if json.loads(path.read_text(encoding="utf-8")) != gate:
            raise ValueError("existing L3 gate has a different identity")
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("x", encoding="utf-8") as stream:
            stream.write(serialized); stream.flush(); os.fsync(stream.fileno())
        os.link(temporary, path)
    except FileExistsError:
        if json.loads(path.read_text(encoding="utf-8")) != gate:
            raise ValueError("concurrent L3 gate has a different identity")
    finally:
        try: temporary.unlink()
        except FileNotFoundError: pass
    return path


def require_l3_gate(study_root: os.PathLike[str] | str) -> dict[str, Any]:
    root = Path(study_root).expanduser().resolve()
    path = root / GATE_NAME
    if not path.is_file():
        raise FileNotFoundError("passed L3 calibration gate not found; L5 is blocked")
    result = validate_l3_gate(json.loads(path.read_text(encoding="utf-8")))
    if result["passed"] is not True:
        raise ValueError("L3 calibration did not pass; L5 is blocked")
    sampling, profile, primary, replay, report = _load_local_evidence(root)
    primary_provenance = primary.get("hpc_provenance")
    replay_provenance = replay.get("hpc_provenance")
    if not isinstance(primary_provenance, Mapping) or not isinstance(replay_provenance, Mapping):
        raise ValueError("L3 execution manifests lack frozen HPC provenance")
    expected = {
        "sampling_manifest_id": sampling["manifest_id"],
        "primary_execution_manifest_id": primary["execution_manifest_id"],
        "replay_execution_manifest_id": replay["execution_manifest_id"],
        "report_id": report["report_id"],
        "source_bundle_id": primary["source_bundle_id"],
        "corpus_id": primary["corpus_id"],
        "corpus_provenance_id": primary_provenance["corpus_provenance_id"],
        "environment_profile_id": profile["environment_profile_id"],
        "criteria": report["criteria"],
    }
    if any(result.get(name) != value for name, value in expected.items()):
        raise ValueError("L3 gate does not match this study's local calibration evidence")
    if (
        report.get("sampling_manifest_id") != sampling["manifest_id"]
        or report.get("primary_execution_manifest_id") != primary["execution_manifest_id"]
        or report.get("replay_execution_manifest_id") != replay["execution_manifest_id"]
        or report.get("source_bundle_id") != primary["source_bundle_id"]
        or report.get("corpus_id") != primary["corpus_id"]
        or report.get("corpus_provenance_id") != primary_provenance["corpus_provenance_id"]
        or report.get("environment_profile_id") != profile["environment_profile_id"]
    ):
        raise ValueError("local L3 report does not match its frozen execution evidence")
    if (
        primary["corpus_id"] != replay["corpus_id"]
        or primary["source_bundle_id"] != replay["source_bundle_id"]
        or primary_provenance != replay_provenance
        or primary_provenance.get("environment_profile_id") != profile["environment_profile_id"]
        or primary_provenance.get("l3_gate_id") is not None
    ):
        raise ValueError("primary/replay L3 execution provenance is inconsistent")
    return result
