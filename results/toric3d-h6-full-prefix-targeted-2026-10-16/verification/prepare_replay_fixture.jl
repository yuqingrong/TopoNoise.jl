#= Prepare a verification-only snapshot fixture; no simulator imports or decoding.

Usage: julia --startup-file=no prepare_replay_fixture.jl SNAPSHOT OUTPUT L P
The selected immutable checkpoint is copied byte-for-byte; raw/metadata mirrors
are reconstructed with the runner's recorded schema, without importing it.
=#
using TOML, SHA

length(ARGS) == 4 || error("Expected SNAPSHOT OUTPUT L P")
source, output = abspath(ARGS[1]), abspath(ARGS[2])
selected_L, selected_p = parse(Int, ARGS[3]), parse(Float64, ARGS[4])
ispath(output) && error("Refusing to replace existing fixture: $output")
digest(bytes) = bytes2hex(SHA.sha256(bytes))
function canonical_toml(data)
    io = IOBuffer(); TOML.print(io, data; sorted=true); String(take!(io))
end
bytes = read(source)
checkpoint = TOML.parse(String(copy(bytes)))
payload = copy(checkpoint)
claimed = pop!(payload, "payload_sha256")
digest(canonical_toml(payload)) == claimed || error("Original payload checksum differs")
snapshot_name = "checkpoint-" * lpad(checkpoint["generation"], 6, '0') * "-" * claimed[1:12] * ".toml"
basename(source) == snapshot_name || error("Original snapshot filename differs")
config = checkpoint["config"]
grid = [(L, p) for L in config["sizes"] for p in config["rates"]]
selected = only(findall(==((selected_L, selected_p)), grid))
state = checkpoint["states"][selected]
0 < state["completed_shots"] <= config["shots"] || error("Selected point has no saved shots")

arms = ("final_only", "full_prefix")
columns = ["channel", "steps", "Ly", "Lz", "p", "shots", "target_shots", "seed", "point_seed", "complete", "elapsed_seconds"]
for arm in arms
    append!(columns, [arm * "_logical_1_failures", arm * "_logical_2_failures", arm * "_failures"])
end
append!(columns, ["joint_" * lpad(mask, 2, '0') for mask in 0:3])
io = IOBuffer(); println(io, join(columns, ','))
for (index, (L, p)) in enumerate(grid)
    saved = checkpoint["states"][index]
    row = Any[config["channel"], L, L, L, p, saved["completed_shots"], config["shots"], config["seed"],
        checkpoint["point_seeds_hex"][index], saved["completed_shots"] == config["shots"], saved["elapsed_seconds"]]
    for arm in 1:2
        append!(row, [saved["logical_1_counts"][arm], saved["logical_2_counts"][arm], saved["failures"][arm]])
    end
    append!(row, saved["joint_counts"])
    println(io, join(row, ','))
end
raw = String(take!(io))
keys_to_copy = ("format_version", "config", "runtime", "source_hashes", "protocol", "point_seeds_hex",
    "decoder_build_seconds_by_size", "decoder_build_count_by_size", "started_at_unix", "updated_at_unix", "generation")
metadata = Dict{String,Any}(key => checkpoint[key] for key in keys_to_copy)
states, target = checkpoint["states"], config["shots"]
merge!(metadata, Dict("expected_points" => length(states),
    "completed_points" => count(s -> s["completed_shots"] == target, states),
    "completed_shots" => sum(s["completed_shots"] for s in states),
    "complete" => all(s["completed_shots"] == target for s in states),
    "raw_sha256" => digest(raw), "checkpoint_sha256" => digest(bytes),
    "current_snapshot" => joinpath("checkpoints", snapshot_name)))
metadata_text = canonical_toml(metadata)
provenance = Dict{String,Any}("fixture_kind" => "verification-only immutable-snapshot mirror",
    "verification_only" => true, "new_independent_fault_records" => 0,
    "must_not_pool_as_new_samples" => true, "decoding_performed_during_preparation" => false,
    "original_snapshot" => source, "original_snapshot_sha256" => digest(bytes),
    "original_payload_sha256" => claimed, "original_generation" => checkpoint["generation"],
    "selected_channel" => config["channel"], "selected_dimensions" => [selected_L, selected_L, selected_L],
    "selected_p" => selected_p, "selected_point_index" => selected,
    "selected_first_shots" => state["completed_shots"], "production_target_shots" => target,
    "selected_point_seed_hex" => checkpoint["point_seeds_hex"][selected],
    "snapshot_completed_shots_by_point" => [s["completed_shots"] for s in states],
    "created_at_unix" => time(), "builder_sha256" => digest(read(@__FILE__)),
    "fixture_checkpoint_sha256" => digest(bytes), "fixture_raw_sha256" => digest(raw),
    "fixture_metadata_sha256" => digest(metadata_text),
    "mirror_note" => "All original grid cells are retained to satisfy independent replay integrity checks. Only the selected point is replayed. One immutable snapshot is copied; this fixture is not a production shard or a complete checkpoint archive.")
mkpath(joinpath(output, "checkpoints"))
write(joinpath(output, "checkpoint.toml"), bytes)
write(joinpath(output, "checkpoints", snapshot_name), bytes)
write(joinpath(output, "raw.csv"), raw)
write(joinpath(output, "metadata.toml"), metadata_text)
write(joinpath(output, "fixture-provenance.toml"), canonical_toml(provenance))
read(source) == bytes == read(joinpath(output, "checkpoint.toml")) == read(joinpath(output, "checkpoints", snapshot_name)) ||
    error("Byte-identical snapshot preservation failed")
println("Prepared verification-only fixture: ", output)
println("Selected original records: ", state["completed_shots"], "; new independent records: 0; decoding performed: false")
println("Original snapshot SHA-256: ", digest(bytes))
