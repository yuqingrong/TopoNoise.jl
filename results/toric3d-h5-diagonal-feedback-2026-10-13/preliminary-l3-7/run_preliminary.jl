using Random, TopoNoise, TOML, SHA

# This is explicitly a three-size diagnostic run. The production CLI's
# minimum-four-size selection guard is intentionally left unchanged.
dir = isempty(ARGS) ? joinpath(@__DIR__,"production") : abspath(only(ARGS))
checkpoint_dir = joinpath(dir,"checkpoints")
saved = TOML.parsefile(joinpath(checkpoint_dir,"checkpoint.toml"))
cfg = saved["config"]
@assert cfg["provenance"]["stage"] == "preliminary_three_size"
@assert cfg["provenance"]["selection_sha256"] == bytes2hex(sha256(read(joinpath(@__DIR__,"selection.json")))) == "1597b9941132d3f24c11c68ec1933ea639f5289225fe80ff5070367e81416e72"
@assert cfg["shots"] == 1000 && cfg["seed"] == 1234
@assert cfg["source_hashes"] == TopoNoise._vertex_source_hashes()
@assert all(g["dimensions"] == [[3,3,3],[5,5,5],[7,7,7]] for g in values(cfg["channel_grids"]))
grids = Dict(Symbol(ch)=>Toric3DVertexScanGrid(g["dimensions"],g["rates"]) for (ch,g) in cfg["channel_grids"])
settings = Toric3DDecoderSettings()
@assert TopoNoise._vertex_settings(settings) == cfg["decoder_settings"]
extras = cfg["provenance"]
@assert TopoNoise._vertex_scan_config(grids,TORIC3D_VERTEX_PROTOCOLS,1000,1234,settings,cfg["batch_size"],extras) == cfg
warm = Dict(ch=>Toric3DVertexScanGrid([(3,3,3)],[0.]) for ch in keys(grids))
zero = compare_toric3d_vertex_projector_feedback(MersenneTwister(1234);channel_grids=warm,shots=1,seed=1234,settings)
@assert all(pt->pt.metrics[:pauli_any].joint_counts[1]==1,zero.points)
pilot = compare_toric3d_vertex_projector_feedback(MersenneTwister(1234);
    channel_grids=Dict(ch=>Toric3DVertexScanGrid([(3,3,3)],[.01]) for ch in keys(grids)),shots=2,seed=1234,settings)
TopoNoise._vertex_atomic_toml(joinpath(dir,"warmup-timing.toml"),Dict(
    "zero_noise_passed"=>true,"paired_records_per_channel"=>2,
    "seconds"=>sum(sum(pt.phase_seconds) for pt in pilot.points),"excluded_from_scientific_counts"=>true))
println("Warmup passed; resuming three-size diagnostics, no threshold eligibility");flush(stdout)
hook = scan -> begin
    save_toric3d_vertex_projector_feedback(scan,dir;figures=false,extra_metadata=extras)
    p=last(scan.points)
    println("Completed ",p.channel," ",p.dimensions," p=",p.p," shots=",p.shots," elapsed=",sum(p.phase_seconds));flush(stdout)
end
scan = compare_toric3d_vertex_projector_feedback(MersenneTwister(1234);
    channel_grids=grids,shots=1000,seed=1234,settings,batch_size=cfg["batch_size"],
    checkpoint_dir,provenance=extras,resume=true,_point_hook=hook)
save_toric3d_vertex_projector_feedback(scan,dir;extra_metadata=extras)
