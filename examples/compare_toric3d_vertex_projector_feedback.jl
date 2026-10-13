using Random,TopoNoise,TOML
function parse_vertex_args(args)
    d=Dict("stage"=>"scout","dimensions"=>"3x3x3,5x5x5,7x7x7","rates"=>"0,0.0001,0.0003,0.001,0.003,0.01,0.02,0.03,0.05,0.075,0.10","channels"=>"x_only,z_only","protocols"=>"none,one_step,window_complete,full_prefix","shots"=>"500","seed"=>"1234","batch-size"=>"50","output"=>"results/toric3d-h5-feedback","resume"=>"false","max-iter"=>"50","osd-order"=>"2")
    allowed=union(Set(keys(d)),Set(("selection","x-dimensions","z-dimensions","x-rates","z-rates","checkpoint-dir")))
    i=1
    while i<=length(args)
        key=replace(args[i],r"^--"=>"");key in allowed||throw(ArgumentError("unknown option $(args[i])"))
        if key=="resume";d[key]="true";i+=1;continue;end
        i<length(args)||throw(ArgumentError("missing value for $key"));d[key]=args[i+1];i+=2
    end
    d["stage"] in ("pilot","scout","extend-scout","production")||throw(ArgumentError("invalid stage"))
    d["protocols"]=="none,one_step,window_complete,full_prefix"||throw(ArgumentError("all four protocols in canonical order are required"))
    d["stage"] in ("production","extend-scout")&&!haskey(d,"selection")&&throw(ArgumentError("stage requires --selection"))
    d
end
_vertex_cli_dims(s)=[Tuple(parse.(Int,split(v,'x'))) for v in split(s,',')]
function main_vertex(args=ARGS)
    d=parse_vertex_args(args);stage=d["stage"];seed=parse(Int,d["seed"]);shots=parse(Int,d["shots"])
    settings=Toric3DDecoderSettings(max_iter=parse(Int,d["max-iter"]),osd_order=parse(Int,d["osd-order"]))
    extras=Dict{String,Any}("stage"=>stage,"extension_round"=>0,"parent_selection_sha256"=>"","parent_scout_inputs"=>Any[])
    if stage in ("production","extend-scout")
        sel=load_toric3d_vertex_production_selection(d["selection"];settings)
        seed==sel.data["seed"]||throw(ArgumentError("selection seed differs"))
        extras["selection_sha256"]=sel.sha256
        if stage=="production"
            sel.data["status"]=="ready"||throw(ArgumentError("production needs a ready selection"))
            # Default is 1000 in production; explicit other shot counts are refused.
            if !("--shots" in args);shots=1000;end
            shots==sel.data["production_shots"]||throw(ArgumentError("production shots mismatch"))
            grids=sel.channel_grids
        else
            sel.data["status"]=="needs_extended_scout"||throw(ArgumentError("extension requires requesting selection"))
            grids=Dict(Symbol(ch)=>Toric3DVertexScanGrid(v["dimensions_T_Lu_Lv"],v["rates"]) for (ch,v) in sel.data["extension_grids"])
            extras["extension_round"]=sel.data["extension_round"]+1;extras["parent_selection_sha256"]=sel.sha256;extras["parent_scout_inputs"]=sel.data["scout_inputs"]
            shots==500||throw(ArgumentError("extension shots must match scout 500"))
        end
    else
        channels=Symbol.(split(d["channels"],','));length(unique(channels))==length(channels)&&all(ch->ch in (:x_only,:z_only),channels)||throw(ArgumentError("invalid channels"))
        grids=Dict(ch=>Toric3DVertexScanGrid(_vertex_cli_dims(get(d,ch==:x_only ? "x-dimensions" : "z-dimensions",d["dimensions"])),parse.(Float64,split(get(d,ch==:x_only ? "x-rates" : "z-rates",d["rates"]),','))) for ch in channels)
    end
    dir=abspath(d["output"])
    checkpoint_dir=get(d,"checkpoint-dir",joinpath(dir,"checkpoints"))
    if d["resume"]=="true"
        checkpoint=joinpath(checkpoint_dir,"checkpoint.toml")
        isfile(checkpoint)||throw(ArgumentError("missing resume checkpoint"))
        expected=TopoNoise._vertex_scan_config(grids,TORIC3D_VERTEX_PROTOCOLS,shots,seed,settings,parse(Int,d["batch-size"]),extras)
        TOML.parsefile(checkpoint)["config"]==expected||throw(ArgumentError("resume configuration/provenance mismatch"))
    end
    mkpath(dir)
    # Zero-noise warmup and a separate timed pilot precede every scientific draw.
    warm=Dict(ch=>Toric3DVertexScanGrid([(3,3,3)],[0.]) for ch in keys(grids))
    zero=compare_toric3d_vertex_projector_feedback(MersenneTwister(seed);channel_grids=warm,shots=1,seed,settings)
    all(pt->pt.metrics[:pauli_any].joint_counts[1]==1,zero.points)||error("zero-noise warmup failed")
    pilot=compare_toric3d_vertex_projector_feedback(MersenneTwister(seed);channel_grids=Dict(ch=>Toric3DVertexScanGrid([(3,3,3)],[.01]) for ch in keys(grids)),shots=2,seed,settings)
    pilot_seconds=sum(sum(pt.phase_seconds) for pt in pilot.points)
    pilot_record=Dict("zero_noise_passed"=>true,"paired_records_per_channel"=>2,"seconds"=>pilot_seconds,"excluded_from_scientific_counts"=>true,"linear_point_projection_seconds"=>pilot_seconds/(2length(grids))*shots*sum(length(g.dimensions)*length(g.error_rates) for g in values(grids)),"projection_limitation"=>"pilot L=3 only; scout size scaling is measured before production selection")
    TopoNoise._vertex_atomic_toml(joinpath(dir,"warmup-timing.toml"),pilot_record)
    println("Warmup passed; timed pilot seconds=",pilot_seconds,"; projected seconds (L=3 baseline)=",pilot_record["linear_point_projection_seconds"]);flush(stdout)
    hook=scan->begin
        save_toric3d_vertex_projector_feedback(scan,dir;figures=false,extra_metadata=extras)
        p=last(scan.points);println("Completed ",p.channel," ",p.dimensions," p=",p.p," shots=",p.shots," elapsed=",sum(p.phase_seconds));flush(stdout)
    end
    scan=compare_toric3d_vertex_projector_feedback(MersenneTwister(seed);channel_grids=grids,shots,seed,settings,
        batch_size=parse(Int,d["batch-size"]),checkpoint_dir,provenance=extras,resume=d["resume"]=="true",_point_hook=hook)
    save_toric3d_vertex_projector_feedback(scan,dir;extra_metadata=extras)
end
abspath(PROGRAM_FILE)==(@__FILE__) && main_vertex()
