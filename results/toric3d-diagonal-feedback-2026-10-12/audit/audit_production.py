#!/usr/bin/env python3
"""Independent consistency audit for the diagonal-feedback production data."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import tomllib


REPOSITORY = Path(__file__).resolve().parents[3]
RESULT_ROOT = REPOSITORY / "results/toric3d-diagonal-feedback-2026-10-12"
PRODUCTION = RESULT_ROOT / "production"
OUTPUT = RESULT_ROOT / "audit/production-artifact-audit.json"
SELECTION = RESULT_ROOT / "scan-selection.toml"

BATCH_SIZE = 32
X_ROUTING = "one_step"
RNG_TYPE = "Random.Xoshiro"
POINT_SOURCE_SNAPSHOT_SHA256 = (
    "b1f9509e77a243b8989bd5147391e76a1e8e1d72b36921e38513139cc37dace5"
)

X_RATES = (0.025, 0.03, 0.04, 0.045, 0.05, 0.055, 0.06, 0.07, 0.08)
Z_FULL_RATES = (0.025, 0.03, 0.0325, 0.035, 0.0375, 0.04, 0.045,
                0.05, 0.06, 0.07)
Z_ANALYSIS_RATES = Z_FULL_RATES[:7]
SHARDS = (
    ("x-L5", PRODUCTION / "x-L5", "x_only", 5, X_RATES, 1000,
     1239, 9, 0),
    ("x-L7", PRODUCTION / "x-L7", "x_only", 7, X_RATES, 1000,
     1241, 9, 0),
    ("x-L9", PRODUCTION / "x-L9", "x_only", 9, X_RATES, 1000,
     1243, 9, 0),
    ("x-L11", PRODUCTION / "x-L11", "x_only", 11, X_RATES, 1000,
     1245, 9, 0),
    ("z-L9", PRODUCTION / "z-L9", "z_only", 9, Z_FULL_RATES, 2000,
     11243, 0, 1),
    ("z-L13", PRODUCTION / "z-L13", "z_only", 13, Z_FULL_RATES, 2000,
     11247, 0, 1),
    ("z-L17", PRODUCTION / "z-L17", "z_only", 17, Z_FULL_RATES, 2000,
     11251, 0, 1),
    ("z-L21-p045", PRODUCTION / "z-L21-p045", "z_only", 21,
     Z_ANALYSIS_RATES, 2000, 11255, 0, 2),
)
ANALYSIS_INPUTS = (
    ("x-L5", PRODUCTION / "x-L5", "x_only", 5, X_RATES, 1000,
     27, 27, True),
    ("x-L7", PRODUCTION / "x-L7", "x_only", 7, X_RATES, 1000,
     27, 27, True),
    ("x-L9", PRODUCTION / "x-L9", "x_only", 9, X_RATES, 1000,
     27, 27, True),
    ("x-L11", PRODUCTION / "x-L11", "x_only", 11, X_RATES, 1000,
     27, 27, True),
    ("z-L9", PRODUCTION / "z-L9", "z_only", 9,
     Z_ANALYSIS_RATES, 2000, 30, 21, False),
    ("z-L13", PRODUCTION / "z-L13", "z_only", 13,
     Z_ANALYSIS_RATES, 2000, 30, 21, False),
    ("z-L17", PRODUCTION / "z-L17", "z_only", 17,
     Z_ANALYSIS_RATES, 2000, 30, 21, False),
    ("z-L21", PRODUCTION / "z-L21-p045", "z_only", 21,
     Z_ANALYSIS_RATES, 2000, 21, 21, True),
)
METRICS = ("logical_1", "logical_2", "state_failure")
PROTOCOLS = {"none": 1, "slice_local": 2, "full_prefix": 4}
PAIRS = {
    "none_vs_slice_local": ("none", "slice_local"),
    "none_vs_full_prefix": ("none", "full_prefix"),
    "slice_local_vs_full_prefix": ("slice_local", "full_prefix"),
}
JOINT_COLUMNS = (
    "joint_n0_s0_f0", "joint_n1_s0_f0", "joint_n0_s1_f0",
    "joint_n1_s1_f0", "joint_n0_s0_f1", "joint_n1_s0_f1",
    "joint_n0_s1_f1", "joint_n1_s1_f1",
)
SETTINGS = {
    "max_iter": 50,
    "bp_method": "minimum_sum",
    "ms_scaling_factor": 1.0,
    "schedule": "parallel",
    "omp_thread_count": 1,
    "osd_method": "OSD_CS",
    "osd_order": 2,
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.relative_to(REPOSITORY))


def recorded_path(value: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = REPOSITORY / path
    return path.resolve()


def read_csv(path: Path):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        fieldnames = reader.fieldnames
        rows = list(reader)
    assert fieldnames is not None
    return fieldnames, rows


def read_toml(path: Path):
    with path.open("rb") as stream:
        return tomllib.load(stream)


def fnv1a64(path: Path) -> str:
    value = 0xCBF29CE484222325
    for byte in path.read_bytes():
        value ^= byte
        value = value * 0x100000001B3 & 0xFFFFFFFFFFFFFFFF
    return f"{value:016x}"


def close(left: float, right: float, tolerance: float = 2e-12) -> bool:
    return math.isclose(left, right, rel_tol=tolerance, abs_tol=tolerance)


def wilson(failures: int, shots: int, z: float = 1.959963984540054):
    probability = failures / shots
    denominator = 1 + z * z / shots
    center = (probability + z * z / (2 * shots)) / denominator
    half = z * math.sqrt(
        probability * (1 - probability) / shots
        + z * z / (4 * shots * shots)
    ) / denominator
    return max(0.0, center - half), min(1.0, center + half)


def pair_counts(joint, first: str, second: str):
    first_bit, second_bit = PROTOCOLS[first], PROTOCOLS[second]
    neither = first_only = second_only = both = 0
    for mask, count in enumerate(joint):
        a, b = bool(mask & first_bit), bool(mask & second_bit)
        if a and b:
            both += count
        elif a:
            first_only += count
        elif b:
            second_only += count
        else:
            neither += count
    return neither, first_only, second_only, both


def audit_assembled_lineage():
    main_path = PRODUCTION / "z-L21/toric3d-diagonal-feedback-checkpoint.toml"
    point_directory = PRODUCTION / "z-L21-p045-point"
    point_path = point_directory / "toric3d-diagonal-feedback-checkpoint.toml"
    point_metadata_path = (
        point_directory / "toric3d-diagonal-feedback-metadata.toml"
    )
    point_derivation_path = point_directory / "derivation-provenance.toml"
    combined_directory = PRODUCTION / "z-L21-p045"
    combined_path = (
        combined_directory / "toric3d-diagonal-feedback-checkpoint.toml"
    )
    combined_metadata_path = (
        combined_directory / "toric3d-diagonal-feedback-metadata.toml"
    )
    combined_derivation_path = combined_directory / "derivation-provenance.toml"
    for path in (
        main_path, point_path, point_metadata_path, point_derivation_path,
        combined_path, combined_metadata_path, combined_derivation_path,
    ):
        assert path.is_file()

    main = read_toml(main_path)
    point = read_toml(point_path)
    point_metadata = read_toml(point_metadata_path)
    point_derivation = read_toml(point_derivation_path)
    combined = read_toml(combined_path)
    combined_metadata = read_toml(combined_metadata_path)
    derivation = read_toml(combined_derivation_path)

    assert main["format_version"] == point["format_version"] == 1
    assert combined["format_version"] == 1
    assert point_derivation["format_version"] == derivation["format_version"] == 1

    # Every available content hash in the assembly record must resolve to the
    # exact file that is being audited. The point's source hash refers to the
    # frozen zero-shot snapshot, which was intentionally superseded in place by
    # the resumable main scan, so that historical identity is pinned directly.
    assert sha256(main_path) == derivation["main_input"]["checkpoint_sha256"]
    assert sha256(point_path) == derivation["point_input"]["checkpoint_sha256"]
    assert sha256(point_metadata_path) == derivation["point_input"][
        "metadata_sha256"]
    assert sha256(combined_path) == derivation["combined_checkpoint_sha256"]
    assert sha256(combined_metadata_path) == derivation[
        "combined_metadata"]["sha256"]
    assert sha256(point_derivation_path) == derivation[
        "point_derivation_provenance"]["sha256"]
    assert point_derivation["source_checkpoint_sha256"] == (
        POINT_SOURCE_SNAPSHOT_SHA256
    )

    assert recorded_path(derivation["main_input"]["checkpoint"]) == (
        main_path.resolve()
    )
    assert recorded_path(derivation["point_input"]["checkpoint"]) == (
        point_path.resolve()
    )
    assert recorded_path(derivation["point_input"]["metadata"]) == (
        point_metadata_path.resolve()
    )
    assert recorded_path(derivation["combined_metadata"]["path"]) == (
        combined_metadata_path.resolve()
    )
    assert recorded_path(
        derivation["point_derivation_provenance"]["path"]
    ) == (
        point_derivation_path.resolve()
    )
    assert recorded_path(point_derivation["source_checkpoint"]) == (
        main_path.resolve()
    )
    assert recorded_path(combined_metadata["derivation_provenance_file"]) == (
        combined_derivation_path.resolve()
    )

    main_indices = derivation["main_input"]["selected_state_indices_one_based"]
    point_indices = derivation["point_input"]["selected_state_indices_one_based"]
    assert main_indices == [1, 2, 3, 4, 5, 6]
    assert point_indices == [1]
    selected_main_states = [main["states"][index - 1] for index in main_indices]
    selected_point_states = [
        point["states"][index - 1] for index in point_indices
    ]
    selected_states = selected_main_states + selected_point_states
    assert combined["states"] == selected_states
    assert combined["point_seeds"] == main["point_seeds"][:7]
    assert point["point_seeds"] == main["point_seeds"][6:7]

    expected_combined_config = dict(main["config"])
    expected_combined_config["error_rates"] = list(Z_ANALYSIS_RATES)
    assert combined["config"] == expected_combined_config
    expected_point_config = dict(main["config"])
    expected_point_config["error_rates"] = [Z_ANALYSIS_RATES[-1]]
    assert point["config"] == expected_point_config

    def check_recorded_selection(record, states, seeds):
        assert record["point_seeds"] == seeds
        for field in (
            "completed_shots", "sampling_seconds", "replay_seconds",
            "invalid_recoveries",
        ):
            assert record[field] == [state[field] for state in states]

    check_recorded_selection(
        derivation["main_input"], selected_main_states,
        main["point_seeds"][:6],
    )
    check_recorded_selection(
        derivation["point_input"], selected_point_states,
        point["point_seeds"],
    )

    assert tuple(derivation["rates"]) == Z_ANALYSIS_RATES
    assert derivation["shots_per_point"] == 2000
    assert derivation["invalid_recoveries"] == 0
    assert derivation["decoder_family_builds"] == 2
    assert derivation["point_seed_p045"] == main["point_seeds"][6]
    assert point_derivation["selected_error_rate"] == Z_ANALYSIS_RATES[-1]
    assert point_derivation["selected_point_seed"] == main["point_seeds"][6]
    assert point_derivation["source_point_index_one_based"] == 7
    assert point_derivation["selected_state_completed_shots"] == 0
    assert point_derivation["selected_state_invalid_recoveries"] == 0
    assert point_metadata["point_seeds_hex"] == point["point_seeds"]
    assert combined_metadata["decoder_builds_x"] == 0
    assert combined_metadata["decoder_builds_z"] == derivation[
        "decoder_family_builds"]

    return {
        "status": "pass",
        "main_state_indices_one_based": main_indices,
        "point_state_indices_one_based": point_indices,
        "combined_states_checked": len(selected_states),
        "main_checkpoint_sha256": sha256(main_path),
        "point_checkpoint_sha256": sha256(point_path),
        "combined_checkpoint_sha256": sha256(combined_path),
        "point_derivation_sha256": sha256(point_derivation_path),
        "combined_derivation_sha256": sha256(combined_derivation_path),
    }


def check_row(row, joint, shots: int):
    assert int(row["shots"]) == shots
    assert int(row["invalid_recoveries"]) == 0
    assert tuple(int(row[column]) for column in JOINT_COLUMNS) == tuple(joint)
    assert int(row["max_iter"]) == SETTINGS["max_iter"]
    assert row["bp_method"] == SETTINGS["bp_method"]
    assert close(float(row["ms_scaling_factor"]),
                 SETTINGS["ms_scaling_factor"])
    assert row["schedule"] == SETTINGS["schedule"]
    assert int(row["omp_thread_count"]) == SETTINGS["omp_thread_count"]
    assert row["osd_method"] == SETTINGS["osd_method"]
    assert int(row["osd_order"]) == SETTINGS["osd_order"]

    for protocol, bit in PROTOCOLS.items():
        failures = sum(count for mask, count in enumerate(joint) if mask & bit)
        assert int(row[f"{protocol}_failures"]) == failures
        assert close(float(row[f"{protocol}_rate"]), failures / shots)
        interval = wilson(failures, shots)
        assert close(float(row[f"{protocol}_ci_low"]), interval[0])
        assert close(float(row[f"{protocol}_ci_high"]), interval[1])

    for name, (first, second) in PAIRS.items():
        values = pair_counts(joint, first, second)
        assert int(row[f"{name}_neither"]) == values[0]
        assert int(row[f"{name}_{first}_only"]) == values[1]
        assert int(row[f"{name}_{second}_only"]) == values[2]
        assert int(row[f"{name}_both"]) == values[3]
        difference = (values[2] - values[1]) / shots
        assert close(float(row[f"{name}_difference_{second}_minus_{first}"]),
                     difference)
        standard_error = math.sqrt(max(
            0.0,
            ((values[1] + values[2]) / shots - difference * difference)
            / (shots - 1),
        ))
        assert close(float(row[f"{name}_standard_error"]), standard_error)
        second_interval = wilson(values[2], shots, z=2.241402727604947)
        first_interval = wilson(values[1], shots, z=2.241402727604947)
        interval = (
            max(-1.0, second_interval[0] - first_interval[1]),
            min(1.0, second_interval[1] - first_interval[0]),
        )
        assert close(float(row[f"{name}_ci_low"]), interval[0])
        assert close(float(row[f"{name}_ci_high"]), interval[1])


def audit_shard(label, directory, channel, size, rates, shots, frozen_seed,
                decoder_builds_x, decoder_builds_z):
    raw_path = directory / "toric3d-diagonal-feedback-raw.csv"
    metadata_path = directory / "toric3d-diagonal-feedback-metadata.toml"
    checkpoint_path = directory / "toric3d-diagonal-feedback-checkpoint.toml"
    assert raw_path.is_file() and metadata_path.is_file() and checkpoint_path.is_file()

    _, rows = read_csv(raw_path)
    metadata = read_toml(metadata_path)
    checkpoint = read_toml(checkpoint_path)

    assert metadata["format_version"] == checkpoint["format_version"] == 1
    assert metadata["complete"] is True
    assert metadata["invalid_recoveries"] == 0
    assert metadata["channels"] == [channel]
    assert metadata["dimensions_steps_Ly_Lz"] == [[size, size, size]]
    assert tuple(metadata["error_rates"]) == tuple(rates)
    assert metadata["shots_per_point"] == shots
    assert metadata["expected_points"] == len(rates)
    assert metadata["completed_points"] == len(rates)
    assert metadata["decoder_settings"] == SETTINGS
    assert metadata["seed_metadata"] == frozen_seed
    assert metadata["batch_size"] == BATCH_SIZE
    assert metadata["x_routing"] == X_ROUTING
    assert metadata["rng_type"] == RNG_TYPE
    assert metadata["point_rng_type"] == RNG_TYPE
    assert metadata["decoder_builds_x"] == decoder_builds_x
    assert metadata["decoder_builds_z"] == decoder_builds_z

    config = checkpoint["config"]
    assert config["dimensions"] == [[size, size, size]]
    assert config["channels"] == [channel]
    assert tuple(config["error_rates"]) == tuple(rates)
    assert config["shots"] == shots
    assert config["settings"] == SETTINGS
    assert config["seed"] == frozen_seed
    assert config["batch_size"] == BATCH_SIZE
    assert config["x_routing"] == X_ROUTING
    assert config["rng_type"] == RNG_TYPE
    assert config["point_rng_type"] == RNG_TYPE
    assert checkpoint["invalid_recoveries"] == 0
    assert len(checkpoint["point_seeds"]) == len(rates)
    assert len(checkpoint["states"]) == len(rates)
    assert metadata["point_seeds_hex"] == checkpoint["point_seeds"]
    assert close(metadata["started_at_unix"], checkpoint["started_at_unix"])

    expected_keys = {
        (float(probability), metric) for probability in rates for metric in METRICS
    }
    identities = [(float(row["p"]), row["metric"]) for row in rows]
    assert len(rows) == len(rates) * len(METRICS)
    assert len(identities) == len(set(identities))
    by_key = dict(zip(identities, rows, strict=True))
    assert set(by_key) == expected_keys

    pymatching_versions = set()
    ldpc_versions = set()
    for point_index, probability in enumerate(rates):
        state = checkpoint["states"][point_index]
        assert state["completed_shots"] == shots
        assert state["invalid_recoveries"] == 0
        point_seed = checkpoint["point_seeds"][point_index]
        for timing in (
            "decoder_build_seconds", "sampling_seconds", "replay_seconds",
        ):
            assert math.isfinite(float(state[timing]))
            assert float(state[timing]) >= 0.0
        pymatching_versions.add(state["pymatching_version"])
        ldpc_versions.add(state["ldpc_version"])
        for metric_index, metric in enumerate(METRICS):
            row = by_key[(probability, metric)]
            assert int(row["steps"]) == int(row["Ly"]) == int(row["Lz"]) == size
            assert row["channel"] == channel
            assert int(row["seed"]) == frozen_seed
            assert row["point_seed"] == point_seed
            assert row["x_routing"] == X_ROUTING
            for timing in (
                "decoder_build_seconds", "sampling_seconds", "replay_seconds",
            ):
                assert close(float(row[timing]), float(state[timing]))
            assert row["pymatching_version"] == state["pymatching_version"]
            assert row["ldpc_version"] == state["ldpc_version"]
            joint = tuple(int(value) for value in state["joint_counts"][metric_index])
            assert sum(joint) == shots
            check_row(row, joint, shots)

    assert pymatching_versions == {"2.4.0"}
    assert len(ldpc_versions) == 1
    if channel == "x_only":
        assert next(iter(ldpc_versions))
    else:
        assert ldpc_versions == {""}

    source_differences = []
    for relative, recorded in metadata["source_hashes"].items():
        path = REPOSITORY / relative
        assert path.is_file()
        if fnv1a64(path) != recorded:
            source_differences.append(relative)

    hashes = {
        path.name: sha256(path)
        for path in sorted(directory.iterdir()) if path.is_file()
    }
    return {
        "status": "pass",
        "channel": channel,
        "L": size,
        "probabilities": list(rates),
        "shots_per_point": shots,
        "csv_rows": len(rows),
        "joint_tables_checked": len(rows),
        "marginal_records_checked": 3 * len(rows),
        "paired_records_checked": 3 * len(rows),
        "recorded_git_commit": metadata["git_state"]["commit"],
        "current_source_differences_from_recorded_snapshot": source_differences,
        "sha256": hashes,
    }


def audit_analysis_inputs():
    manifest_path = RESULT_ROOT / "analysis-input/manifest.json"
    assert manifest_path.is_file()
    manifest = json.loads(manifest_path.read_text())

    assert manifest["format_version"] == 1
    assert manifest["selection_file"] == relative(SELECTION)
    assert manifest["selection_sha256"] == sha256(SELECTION)
    assert manifest["x_state_failure_points"] == 36
    assert manifest["x_paired_fault_records"] == 36000
    assert manifest["z_state_failure_points"] == 28
    assert manifest["z_paired_fault_records"] == 56000
    entries = manifest["inputs"]
    expected_labels = [spec[0] for spec in ANALYSIS_INPUTS]
    labels = [entry["label"] for entry in entries]
    assert labels == expected_labels
    assert len(labels) == len(set(labels)) == len(ANALYSIS_INPUTS)

    seen = set()
    rows_checked = 0
    state_failure_points = {"x_only": 0, "z_only": 0}
    paired_fault_records = {"x_only": 0, "z_only": 0}
    for entry, spec in zip(entries, ANALYSIS_INPUTS, strict=True):
        (label, source_directory, channel, size, rates, shots,
         source_row_count, derived_row_count, byte_copy) = spec
        source = source_directory / "toric3d-diagonal-feedback-raw.csv"
        source_metadata = (
            source_directory / "toric3d-diagonal-feedback-metadata.toml"
        )
        derived = RESULT_ROOT / f"analysis-input/{label}.csv"
        assert entry["label"] == label
        assert entry["channel"] == channel
        assert entry["size"] == size
        assert [Decimal(value) for value in entry["probabilities"]] == [
            Decimal(str(value)) for value in rates
        ]
        assert entry["shots_per_point"] == shots
        assert entry["source_csv"] == relative(source)
        assert entry["source_metadata"] == relative(source_metadata)
        assert entry["derived_csv"] == relative(derived)
        assert entry["byte_copy"] is byte_copy

        assert source.is_file() and source_metadata.is_file() and derived.is_file()
        assert sha256(derived) == entry["derived_csv_sha256"]
        assert sha256(source) == entry["source_csv_sha256"]
        assert sha256(source_metadata) == entry["source_metadata_sha256"]

        source_fields, source_rows = read_csv(source)
        derived_fields, derived_rows = read_csv(derived)
        selected_rows = [
            row for row in source_rows
            if Decimal(row["p"]) in {Decimal(str(value)) for value in rates}
        ]
        assert derived_fields == source_fields
        assert derived_rows == selected_rows
        assert entry["source_rows"] == len(source_rows) == source_row_count
        assert entry["derived_rows"] == len(derived_rows) == derived_row_count
        assert derived_row_count == len(rates) * len(METRICS)
        assert byte_copy == (len(selected_rows) == len(source_rows))
        if byte_copy:
            assert derived.read_bytes() == source.read_bytes()

        rows_checked += len(derived_rows)
        for row in derived_rows:
            key = (row["channel"], int(row["steps"]), float(row["p"]), row["metric"])
            assert key not in seen
            seen.add(key)
            assert row["channel"] == channel
            assert int(row["steps"]) == int(row["Ly"]) == int(row["Lz"]) == size
            assert int(row["shots"]) == shots
            joint = tuple(int(row[column]) for column in JOINT_COLUMNS)
            check_row(row, joint, shots)
            if row["metric"] == "state_failure":
                state_failure_points[channel] += 1
                paired_fault_records[channel] += shots

    assert state_failure_points == {"x_only": 36, "z_only": 28}
    assert paired_fault_records == {"x_only": 36000, "z_only": 56000}
    return {
        "status": "pass",
        "manifest": str(manifest_path.relative_to(REPOSITORY)),
        "manifest_sha256": sha256(manifest_path),
        "csv_rows_checked": rows_checked,
        "x_state_failure_points": state_failure_points["x_only"],
        "z_state_failure_points": state_failure_points["z_only"],
        "x_paired_fault_records": paired_fault_records["x_only"],
        "z_paired_fault_records": paired_fault_records["z_only"],
    }


def main():
    result = {
        "format_version": 1,
        "audit_kind": "independent final production and analysis-input audit",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "selection_sha256": sha256(SELECTION),
        "assembled_lineage": audit_assembled_lineage(),
        "shards": {},
    }
    for arguments in SHARDS:
        result["shards"][arguments[0]] = audit_shard(*arguments)
    result["analysis_inputs"] = audit_analysis_inputs()
    result["totals"] = {
        "complete_shards": len(SHARDS),
        "joint_tables_checked": sum(
            shard["joint_tables_checked"] for shard in result["shards"].values()),
        "marginal_records_checked": sum(
            shard["marginal_records_checked"] for shard in result["shards"].values()),
        "paired_records_checked": sum(
            shard["paired_records_checked"] for shard in result["shards"].values()),
        "invalid_recoveries": 0,
    }
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(OUTPUT.relative_to(REPOSITORY))


if __name__ == "__main__":
    main()
