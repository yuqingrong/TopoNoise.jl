using Random, TOML, SHA

parent = normpath(joinpath(@__DIR__, "..", "preliminary-l3-9"))
source = joinpath(parent, "production", "checkpoints", "checkpoint.toml")
positional = filter(a -> a != "--fresh", ARGS)
output = isempty(positional) ? joinpath(@__DIR__, "production") : abspath(only(positional))
target = joinpath(output, "checkpoints", "checkpoint.toml")
isfile(target) && error("Refusing to overwrite an existing checkpoint")
saved = TOML.parsefile(source)
original = deepcopy(saved)
@assert bytes2hex(sha256(read(source))) == "c3544f0d3f72ff2cd5fe0e79379ec502b0360cdc15e76a3f7b9cff4a18faa80c"
@assert saved["config"]["shots"] == 1000 && saved["config"]["seed"] == 1234
@assert length(saved["states"]) == 260
@assert saved["point_seeds"] == string.(rand(Xoshiro(1234), UInt64, 260); base=16, pad=16)
@assert all(s["invalid_recoveries"] == 0 for s in saved["states"])
config = saved["config"]
specs = [(ch,d,p) for ch in ("x_only","z_only") for d in config["channel_grids"][ch]["dimensions"] for p in config["channel_grids"][ch]["rates"]]
keep = findall(s -> s[2][1] <= 7, specs)
@assert length(keep) == 195
saved["states"] = saved["states"][keep]
saved["point_seeds"] = saved["point_seeds"][keep]
reused_records = sum(s["completed_shots"] for s in saved["states"])
if "--fresh" in ARGS
    for s in saved["states"]
        s["completed_shots"] = 0
        s["counts"] = [zeros(Int,16) for _ in 1:23]
        for k in ("decoder_build_seconds","sampling_seconds","replay_seconds")
            s[k] = 0.0
        end
        s["record_digest"] = s["signature_digest"] = "cbf29ce484222325"
    end
    reused_records = 0
end
for grid in values(config["channel_grids"])
    filter!(d -> d[1] <= 7, grid["dimensions"])
end
config["provenance"]["stage"] = "preliminary_three_size"
selection_hash = bytes2hex(sha256(read(joinpath(@__DIR__,"selection.json"))))
@assert selection_hash == "1597b9941132d3f24c11c68ec1933ea639f5289225fe80ff5070367e81416e72"
config["provenance"]["selection_sha256"] = selection_hash
mkpath(dirname(target))
open(target*".tmp","w") do io
    TOML.print(io,saved;sorted=true)
end
mv(target*".tmp",target)
@assert TOML.parsefile(target) == saved
@assert TOML.parsefile(source) == original
manifest = Dict("source_checkpoint_path"=>source,
    "source_checkpoint_sha256"=>bytes2hex(sha256(read(source))),
    "initial_checkpoint_sha256"=>bytes2hex(sha256(read(target))),
    "selection_sha256"=>selection_hash,"sizes"=>[3,5,7],
    "parent_indices_1based"=>keep,"reused_paired_records"=>reused_records,
    "target_paired_records"=>195000,"fresh"=>("--fresh" in ARGS),
    "seed_rule"=>"Random.Xoshiro(1234) UInt64 vector of length260, subset at parent_indices_1based; retain every existing point seed and record",
    "source_preserved"=>true,"threshold_eligible"=>false)
open(joinpath(output,"checkpoint-reuse.toml"),"w") do io
    TOML.print(io,manifest;sorted=true)
end
println("PASS: ",reused_records," records retained with exact counts, digests and seeds; ",195000-reused_records," records remain")
