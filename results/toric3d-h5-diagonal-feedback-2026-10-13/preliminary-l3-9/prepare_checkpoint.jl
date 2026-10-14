using Random, TOML, SHA

# Reuse the completed X prefix only when it exactly matches a fresh reduced
# run's point-seed stream. The unstarted Z records get that fresh stream.
root = normpath(joinpath(@__DIR__, ".."))
source = joinpath(root, "production", "checkpoints", "checkpoint.toml")
target = joinpath(@__DIR__, "production", "checkpoints", "checkpoint.toml")
isfile(target) && error("Target exists; refusing to overwrite a preliminary run")
saved = TOML.parsefile(source)
original = deepcopy(saved)
config = saved["config"]
@assert config["shots"] == 1000 && config["seed"] == 1234
@assert all(s["invalid_recoveries"] == 0 for s in saved["states"])
@assert bytes2hex(sha256(read(joinpath(root,"selection.json")))) == "5041dd833aa8f28f318b4f5fb495f5025d61dfeb87f1a8ab8d76e5edcf6d7e00"
specs = [(ch,d,p) for ch in ("x_only","z_only") for d in config["channel_grids"][ch]["dimensions"] for p in config["channel_grids"][ch]["rates"]]
keep = findall(s -> s[2][1] <= 9, specs)
small_specs = specs[keep]
states = saved["states"][keep]
seeds = string.(rand(Xoshiro(config["seed"]), UInt64, length(keep)); base=16, pad=16)
reused = Int[]
for (i,(ch,d,p)) in enumerate(small_specs)
    if ch == "x_only"
        @assert states[i]["completed_shots"] == 1000
        @assert seeds[i] == saved["point_seeds"][keep[i]]
        push!(reused,keep[i])
    else
        @assert states[i]["completed_shots"] == 0
        @assert all(iszero,Iterators.flatten(states[i]["counts"]))
    end
end
@assert length(reused) == 80 && length(states) == 260
for grid in values(config["channel_grids"])
    filter!(d -> d[1] <= 9, grid["dimensions"])
end
selection_hash = bytes2hex(sha256(read(joinpath(@__DIR__,"selection.json"))))
config["provenance"]["selection_sha256"] = selection_hash
saved["states"] = states
saved["point_seeds"] = seeds
mkpath(dirname(target))
open(target*".tmp","w") do io
    TOML.print(io,saved;sorted=true)
end
mv(target*".tmp",target)
@assert TOML.parsefile(target) == saved
@assert TOML.parsefile(source) == original
manifest = Dict(
    "source_checkpoint_path"=>source,
    "source_checkpoint_sha256"=>bytes2hex(sha256(read(source))),
    "initial_checkpoint_sha256"=>bytes2hex(sha256(read(target))),
    "selection_sha256"=>selection_hash,
    "sizes"=>[3,5,7,9], "reused_points"=>80, "reused_paired_records"=>80000,
    "reused_original_indices_1based"=>reused,
    "new_z_points"=>180, "new_z_paired_records"=>180000,
    "seed_rule"=>"Fresh Random.Xoshiro(1234) UInt64 vector of length260; every reused X seed matches original exactly; no original Z data existed",
    "counts_and_digests"=>"Complete reused X states copied exactly, including all23 joint-count arrays and record/signature digests",
    "source_preserved"=>true,
    "reason"=>"User requested smaller lattices first on2026-10-14; original large scan deliberately interrupted and deferred")
open(joinpath(@__DIR__,"checkpoint-reuse.toml"),"w") do io
    TOML.print(io,manifest;sorted=true)
end
println("PASS: 80 completed X points reused exactly; 180 fresh Z points; original checkpoint unchanged")
