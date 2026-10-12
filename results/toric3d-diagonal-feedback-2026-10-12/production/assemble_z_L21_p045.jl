#!/usr/bin/env julia

# Combine the six completed points from the resumable ten-rate L=21 shard with
# the independently sharded seventh point. No scientific shots are sampled.

ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
ENV["JULIA_PYTHONCALL_EXE"] =
    joinpath(pwd(), ".CondaPkg/.pixi/envs/default/bin/python")

using SHA
using TOML
using TopoNoise

const ROOT = "results/toric3d-diagonal-feedback-2026-10-12/production"
const MAIN_PATH = joinpath(ROOT, "z-L21",
    "toric3d-diagonal-feedback-checkpoint.toml")
const POINT_PATH = joinpath(ROOT, "z-L21-p045-point",
    "toric3d-diagonal-feedback-checkpoint.toml")
const POINT_METADATA = joinpath(ROOT, "z-L21-p045-point",
    "toric3d-diagonal-feedback-metadata.toml")
const POINT_DERIVATION = joinpath(ROOT, "z-L21-p045-point",
    "derivation-provenance.toml")
const OUTPUT = joinpath(ROOT, "z-L21-p045")
const CHECKPOINT = joinpath(OUTPUT,
    "toric3d-diagonal-feedback-checkpoint.toml")
const METADATA = joinpath(OUTPUT,
    "toric3d-diagonal-feedback-metadata.toml")
const DERIVATION = joinpath(OUTPUT, "derivation-provenance.toml")
const RATES = Float64[0.025, 0.03, 0.0325, 0.035, 0.0375, 0.04, 0.045]

sha256(path) = bytes2hex(SHA.sha256(read(path)))

main = TOML.parsefile(MAIN_PATH)
point = TOML.parsefile(POINT_PATH)
point_metadata = TOML.parsefile(POINT_METADATA)
point_derivation = TOML.parsefile(POINT_DERIVATION)

@assert main["format_version"] == point["format_version"] == 1
@assert main["config"]["dimensions"] == point["config"]["dimensions"] ==
    [[21, 21, 21]]
@assert main["config"]["channels"] == point["config"]["channels"] == ["z_only"]
@assert main["config"]["shots"] == point["config"]["shots"] == 2000
@assert main["config"]["batch_size"] == point["config"]["batch_size"] == 32
@assert main["config"]["seed"] == point["config"]["seed"] == 11255
@assert main["config"]["x_routing"] == point["config"]["x_routing"] == "one_step"
@assert main["config"]["settings"] == point["config"]["settings"]
@assert main["config"]["error_rates"][1:7] == RATES
@assert point["config"]["error_rates"] == [0.045]
@assert main["point_seeds"][7] == point["point_seeds"][1] ==
    "dedd14bb48a97cc7"
@assert point_derivation["format_version"] == 1
@assert point_derivation["selected_error_rate"] == 0.045
@assert point_derivation["selected_point_seed"] == "dedd14bb48a97cc7"
@assert point_derivation["source_point_index_one_based"] == 7
@assert point_derivation["selected_state_completed_shots"] == 0
@assert point_derivation["selected_state_invalid_recoveries"] == 0
@assert abspath(point_derivation["source_checkpoint"]) == abspath(MAIN_PATH)
@assert point_metadata["complete"] == true
@assert point_metadata["invalid_recoveries"] == 0
@assert point_metadata["point_seeds_hex"] == ["dedd14bb48a97cc7"]
@assert all(Int(state["completed_shots"]) == 2000
    for state in main["states"][1:6])
@assert Int(point["states"][1]["completed_shots"]) == 2000
@assert all(Int(state["invalid_recoveries"]) == 0
    for state in main["states"][1:6])
@assert Int(point["states"][1]["invalid_recoveries"]) == 0
@assert all(sum(Int.(row)) == 2000
    for state in vcat(main["states"][1:6], point["states"][1:1])
    for row in state["joint_counts"])

combined = deepcopy(main)
combined["config"]["error_rates"] = RATES
combined["point_seeds"] = main["point_seeds"][1:7]
combined["states"] = vcat(main["states"][1:6], point["states"][1:1])
combined["invalid_recoveries"] = 0

mkpath(OUTPUT)
temporary = CHECKPOINT * ".tmp"
open(temporary, "w") do io
    TOML.print(io, combined; sorted=true)
end
mv(temporary, CHECKPOINT; force=true)

point_seeds = TopoNoise._toric3d_feedback_parse_seed.(combined["point_seeds"])
states = [TopoNoise._toric3d_feedback_state(state, 2000)
    for state in combined["states"]]
points = [TopoNoise._toric3d_feedback_point(
    (21, 21, 21), RATES[index], :z_only, 2000,
    point_seeds[index], states[index]) for index in eachindex(RATES)]
@assert all(result.shots == 2000 for result in points)

# Reconstruct the completed scan directly from validated checkpoint states.
# This avoids rebuilding a decoder or sampling a fault during assembly. Two Z
# families were built in the two disjoint production processes.
scan = Toric3DDiagonalFeedbackScan(
    [(21, 21, 21)], RATES, [:z_only], 2000, 11255, 32, :one_step,
    Toric3DDecoderSettings(), points, "Random.Xoshiro", point_seeds,
    Float64(combined["started_at_unix"]), true, CHECKPOINT,
    joinpath(OUTPUT, "toric3d-diagonal-feedback-failure.toml"), 0,
    (x=0, z=2),
)
paths = save_toric3d_diagonal_feedback(scan, OUTPUT)

# save_toric3d_diagonal_feedback records the checkout used to assemble and
# render this combined artifact. The shots came from two earlier processes, so
# label that scope explicitly and point to a separate immutable lineage file.
metadata = TOML.parsefile(METADATA)
metadata["source_hash_role"] =
    "assembly and rendering environment only; not simulation execution"
metadata["derivation_provenance_file"] = DERIVATION
metadata["execution_source_snapshot_status"] = Dict{String,Any}(
    "main" =>
        "partial: the checkpoint did not embed a source snapshot; the loaded scan source is identified by archived pre-fix interrupt behavior as FNV-1a-64 18a0d574a9cff285",
    "point" =>
        "recorded in the independently generated point metadata referenced by the derivation provenance",
    "scientific_delta" =>
        "the main and assembly scan sources differ only by rethrowing InterruptException before bookkeeping; completed scientific batches are unchanged",
)
metadata_temporary = METADATA * ".tmp"
open(metadata_temporary, "w") do io
    TOML.print(io, metadata; sorted=true)
end
mv(metadata_temporary, METADATA; force=true)

provenance = Dict{String,Any}(
    "format_version" => 1,
    "derivation" => "states 1:6 from the ten-rate L=21 shard plus the independently replayed frozen state 7 at p=0.045; no shots resampled during assembly",
    "combined_checkpoint_sha256" => sha256(CHECKPOINT),
    "combined_metadata" => Dict{String,Any}(
        "path" => METADATA,
        "sha256" => sha256(METADATA),
    ),
    "point_derivation_provenance" => Dict{String,Any}(
        "path" => POINT_DERIVATION,
        "sha256" => sha256(POINT_DERIVATION),
    ),
    "point_seed_p045" => main["point_seeds"][7],
    "rates" => RATES,
    "shots_per_point" => 2000,
    "invalid_recoveries" => 0,
    "decoder_family_builds" => 2,
    "assembly_source_hash_role" =>
        "assembly and rendering only; execution lineage is recorded by input",
    "main_input" => Dict{String,Any}(
        "role" => "states 1:6, p=0.025 through p=0.04",
        "checkpoint" => MAIN_PATH,
        "checkpoint_sha256" => sha256(MAIN_PATH),
        "selected_state_indices_one_based" => collect(1:6),
        "point_seeds" => main["point_seeds"][1:6],
        "completed_shots" =>
            [state["completed_shots"] for state in main["states"][1:6]],
        "sampling_seconds" =>
            [state["sampling_seconds"] for state in main["states"][1:6]],
        "replay_seconds" =>
            [state["replay_seconds"] for state in main["states"][1:6]],
        "invalid_recoveries" =>
            [state["invalid_recoveries"] for state in main["states"][1:6]],
        "execution_source_snapshot_status" =>
            "full snapshot not embedded in checkpoint; loaded feedback_strength_scan.jl is the pre-1e4dbe9 source with FNV-1a-64 18a0d574a9cff285, established by the archived pre-fix interrupt behavior",
    ),
    "point_input" => Dict{String,Any}(
        "role" => "state 7, p=0.045",
        "checkpoint" => POINT_PATH,
        "checkpoint_sha256" => sha256(POINT_PATH),
        "metadata" => POINT_METADATA,
        "metadata_sha256" => sha256(POINT_METADATA),
        "selected_state_indices_one_based" => [1],
        "point_seeds" => point["point_seeds"],
        "completed_shots" => [point["states"][1]["completed_shots"]],
        "sampling_seconds" => [point["states"][1]["sampling_seconds"]],
        "replay_seconds" => [point["states"][1]["replay_seconds"]],
        "invalid_recoveries" => [point["states"][1]["invalid_recoveries"]],
        "execution_source_snapshot_status" =>
            "recorded in the independently generated point metadata",
    ),
    "saved_paths" => string(paths),
)
open(DERIVATION, "w") do io
    TOML.print(io, provenance; sorted=true)
end

println(paths)
