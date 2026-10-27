const _VERTEX_METRICS=(:logical_x_1,:logical_x_2,:logical_x_union,:logical_z_1,:logical_z_2,:logical_z_union,:pauli_any,(Symbol("pauli_class_",i) for i in 0:15)...)
struct Toric3DVertexScanGrid
    dimensions::Vector{NTuple{3,Int}}
    error_rates::Vector{Float64}
    function Toric3DVertexScanGrid(dimensions,error_rates)
        ds=NTuple{3,Int}[Tuple(d) for d in dimensions];ps=Float64.(collect(error_rates))
        !isempty(ds)&&length(unique(ds))==length(ds) || throw(ArgumentError("dimensions must be nonempty and unique"))
        foreach(d->Toric3DVertexProjectorCode(d...),ds)
        !isempty(ps)&&length(unique(ps))==length(ps)&&all(p->isfinite(p)&&0<=p<=1,ps) || throw(ArgumentError("rates must be unique finite probabilities"))
        new(ds,sort(ps))
    end
end
struct Toric3DVertexProtocolEstimate
    failures::Int
    rate::Float64
    ci::Tuple{Float64,Float64}
end
struct Toric3DVertexPairEstimate
    first::Symbol
    second::Symbol
    neither::Int
    first_only::Int
    second_only::Int
    both::Int
    difference::Float64
    standard_error::Float64
    ci::Tuple{Float64,Float64}
end
struct Toric3DVertexQuadMetricEstimate
    shots::Int
    joint_counts::NTuple{16,Int}
    protocols::NTuple{4,Toric3DVertexProtocolEstimate}
    pairs::NTuple{6,Toric3DVertexPairEstimate}
end
function _vertex_quad_metric(shots,raw)
    shots>0 && length(raw)==16 && all(>=(0),raw) && sum(raw)==shots || throw(ArgumentError("invalid sixteen-way counts"))
    counts=Tuple(Int.(raw));estimates=map(1:4) do i
        k=sum(counts[m+1] for m in 0:15 if !iszero(m&(1<<(i-1))))
        Toric3DVertexProtocolEstimate(k,k/shots,_toric3d_wilson(k,shots))
    end
    pairs=Toric3DVertexPairEstimate[]
    for i in 1:3,j in i+1:4
        cells=zeros(Int,4)
        for m in 0:15
            a=!iszero(m&(1<<(i-1)));b=!iszero(m&(1<<(j-1)))
            cells[1+Int(a)+2Int(b)]+=counts[m+1]
        end
        n,a,b,ab=cells;delta=(b-a)/shots
        se=shots>1 ? sqrt(max(0.,((a+b)/shots-delta^2)/(shots-1))) : NaN
        aci=_toric3d_wilson(a,shots;z=2.241402727604947);bci=_toric3d_wilson(b,shots;z=2.241402727604947)
        push!(pairs,Toric3DVertexPairEstimate(TORIC3D_VERTEX_PROTOCOLS[i],TORIC3D_VERTEX_PROTOCOLS[j],n,a,b,ab,delta,se,(max(-1.,bci[1]-aci[2]),min(1.,bci[2]-aci[1]))))
    end
    Toric3DVertexQuadMetricEstimate(shots,counts,Tuple(estimates),Tuple(pairs))
end
struct Toric3DVertexFeedbackPoint
    dimensions::NTuple{3,Int}
    p::Float64
    channel::Symbol
    shots::Int
    point_seed::UInt64
    metrics::Dict{Symbol,Toric3DVertexQuadMetricEstimate}
    phase_seconds::NamedTuple{(:decoder_build,:sampling,:replay),Tuple{Float64,Float64,Float64}}
    invalid_recoveries::Int
    record_digest::String
    signature_digest::String
end
struct Toric3DVertexFeedbackScan
    channel_grids::Dict{Symbol,Toric3DVertexScanGrid}
    shots::Int
    seed::Union{Nothing,Int}
    settings::Toric3DDecoderSettings
    points::Vector{Toric3DVertexFeedbackPoint}
    config::Dict{String,Any}
    checkpoint_path::Union{Nothing,String}
    resumed::Bool
end
_vertex_settings(s)=Dict(string(k)=>v for (k,v) in pairs(_toric3d_settings_metadata(s)))
function _vertex_sha256(path)
    f=PythonCall.pybuiltins.open(path,"rb")
    try
        PythonCall.pyconvert(String,PythonCall.pyimport("hashlib").sha256(f.read()).hexdigest())
    finally
        f.close()
    end
end
function _vertex_source_hashes()
    root=normpath(joinpath(@__DIR__,"..",".."));result=Dict{String,String}()
    for (dir,_,files) in walkdir(joinpath(root,"src")),file in sort(files)
        endswith(file,".jl") || continue
        path=joinpath(dir,file);result[relpath(path,root)]=_vertex_sha256(path)
    end
    result["CondaPkg.toml"]=_vertex_sha256(joinpath(root,"CondaPkg.toml"));result
end
function _vertex_versions()
    versions=Dict("pymatching"=>_toric3d_module_version(_load_pymatching(),"pymatching"),"ldpc"=>_toric3d_module_version(_toric3d_load_ldpc(),"ldpc"))
    versions==Dict("pymatching"=>"2.4.0","ldpc"=>"2.4.1") || throw(ArgumentError("H+5 backend versions differ from pins"))
    versions
end
function _vertex_scan_config(grids,protocols,shots,seed,settings,batch_size,provenance=Dict{String,Any}())
    Dict{String,Any}("format_version"=>1,"protocol_version"=>"h5-four-arm-v1","protocols"=>string.(collect(protocols)),
        "channel_grids"=>Dict(string(ch)=>Dict("dimensions"=>collect.(g.dimensions),"rates"=>g.error_rates) for (ch,g) in grids),
        "shots"=>shots,"seed"=>something(seed,"unspecified"),"batch_size"=>batch_size,"decoder_settings"=>_vertex_settings(settings),
        "provenance"=>deepcopy(provenance),"backend_versions"=>_vertex_versions(),"source_hashes"=>_vertex_source_hashes(),"point_rng"=>"Random.Xoshiro","julia_version"=>string(VERSION))
end
_vertex_state()=Dict{String,Any}("completed_shots"=>0,"counts"=>[zeros(Int,16) for _ in _VERTEX_METRICS],"decoder_build_seconds"=>0.,"sampling_seconds"=>0.,"replay_seconds"=>0.,"invalid_recoveries"=>0,"record_digest"=>"cbf29ce484222325","signature_digest"=>"cbf29ce484222325")
function _vertex_digest(hex,bytes)
    h=parse(UInt64,hex;base=16)
    for b in bytes;h=(h⊻UInt64(b))*UInt64(0x100000001b3);end
    string(h;base=16,pad=16)
end
_vertex_record_bytes(r)=UInt8[vcat(vec(r.x),vec(r.z),collect(r.upper_branch))...]
function _vertex_signature_bytes(out)
    # Canonical numeric signatures, with explicit checkpoint/channel delimiters.
    io=IOBuffer()
    print(io,out.protocol,':')
    for h in out.history;print(io,h.checkpoint,':',h.channel,':',h.signature,';');end
    print(io,"final:",out.z_decision.signature)
    take!(io)
end
function _vertex_atomic_toml(path,payload)
    mkpath(dirname(path))
    _toric3d_atomic_write(path) do io;TOML.print(io,payload;sorted=true);end
end
function _vertex_validate_states(states,n,shots)
    length(states)==n || throw(ArgumentError("checkpoint point count differs"))
    for s in states
        c=s["completed_shots"];0<=c<=shots || throw(ArgumentError("invalid completed shots"))
        length(s["counts"])==length(_VERTEX_METRICS) && all(row->length(row)==16&&all(>=(0),row)&&sum(row)==c,s["counts"]) || throw(ArgumentError("invalid checkpoint counts"))
        s["invalid_recoveries"]==0 || throw(ArgumentError("checkpoint contains invalid recovery; diagnose before restarting"))
        all(k->isfinite(s[k])&&s[k]>=0,("decoder_build_seconds","sampling_seconds","replay_seconds")) || throw(ArgumentError("invalid checkpoint timing"))
        for k in ("record_digest","signature_digest");parse(UInt64,s[k];base=16);end
    end
end
function _vertex_metrics(out)
    (out.logical_x...,out.logical_x_failure,out.logical_z...,out.logical_z_failure,out.pauli_class!=0,(out.pauli_class==i for i in 0:15)...)
end
function _vertex_boundary_context(e,r,frame)
    b=complete_toric3d_vertex_boundary(toric3d_vertex_boundary_map(e.code),frame.z,r.upper_branch)
    Dict{String,Any}("upper_branch"=>collect(r.upper_branch),"observed_branch"=>b.observed,"physical_correction"=>b.physical_correction,
        "clean_correction"=>b.clean_correction,"frame_x"=>frame.x,"frame_z"=>frame.z,"fault_x"=>[collect(r.x[:,i]) for i in axes(r.x,2)],"fault_z"=>[collect(r.z[:,i]) for i in axes(r.z,2)])
end
"""Four protocols replay one record per shot, with atomic batches and exact point substreams."""
function compare_toric3d_vertex_projector_feedback(rng::Random.AbstractRNG;
        channel_grids,protocols=TORIC3D_VERTEX_PROTOCOLS,shots::Integer,seed::Union{Nothing,Integer}=nothing,
        settings::Toric3DDecoderSettings=Toric3DDecoderSettings(),batch_size::Integer=50,
        checkpoint_dir=nothing,resume::Bool=false,provenance=Dict{String,Any}(),_draw_hook=nothing,_recovery_hook=nothing,_point_hook=nothing)
    Tuple(protocols)==TORIC3D_VERTEX_PROTOCOLS || throw(ArgumentError("protocols must equal canonical four-arm tuple"))
    shots>0&&batch_size>0 || throw(ArgumentError("shots and batch_size must be positive"))
    grids=Dict{Symbol,Toric3DVertexScanGrid}(channel_grids)
    !isempty(grids)&&all(ch->ch in (:x_only,:z_only),keys(grids)) || throw(ArgumentError("scan channels must be x_only or z_only"))
    allowed=Set(("stage","selection_sha256","extension_round","parent_selection_sha256","parent_scout_inputs"))
    all(k->k in allowed,keys(provenance))||throw(ArgumentError("unknown scientific provenance key"))
    config=_vertex_scan_config(grids,protocols,Int(shots),seed,settings,Int(batch_size),provenance)
    specs=[(ch,d,p) for ch in (:x_only,:z_only) if haskey(grids,ch) for d in grids[ch].dimensions for p in grids[ch].error_rates]
    path=isnothing(checkpoint_dir) ? nothing : joinpath(abspath(checkpoint_dir),"checkpoint.toml")
    resume&&isnothing(path) && throw(ArgumentError("resume needs checkpoint_dir"))
    if resume
        isfile(path) || throw(ArgumentError("missing resume checkpoint"))
        saved=TOML.parsefile(path)
        saved["config"]==config || throw(ArgumentError("resume configuration mismatch"))
        states=saved["states"];seeds=_toric3d_feedback_parse_seed.(saved["point_seeds"])
        length(seeds)==length(specs) || throw(ArgumentError("checkpoint seeds differ"))
        _vertex_validate_states(states,length(specs),shots)
    else
        !isnothing(path)&&isfile(path) && throw(ArgumentError("checkpoint exists; use resume or a fresh directory"))
        states=[_vertex_state() for _ in specs]
        source=isnothing(seed) ? rng : Random.Xoshiro(seed);seeds=rand(source,UInt64,length(specs))
    end
    checkpoint()=isnothing(path) ? nothing : _vertex_atomic_toml(path,Dict("config"=>config,"states"=>states,"point_seeds"=>_toric3d_feedback_seed_string.(seeds)))
    checkpoint();points=Toric3DVertexFeedbackPoint[]
    for (idx,(ch,dims,p)) in enumerate(specs)
        s=states[idx];n=s["completed_shots"]
        if n<shots
            c=Toric3DVertexProjectorCode(dims...);e=toric3d_vertex_projector_encoder(c)
            started=time_ns();f=build_toric3d_vertex_decoders(c;p,settings);s["decoder_build_seconds"]+=(time_ns()-started)/1e9
            prng=Random.Xoshiro(seeds[idx]);sample()=sample_toric3d_vertex_projector_record(prng,e;p_x=ch==:x_only ? p : 0.,p_z=ch==:z_only ? p : 0.)
            for _ in 1:n;sample();end # includes every clean upper branch draw
            for shot in n+1:shots
                started=time_ns();r=sample();sampling=(time_ns()-started)/1e9
                isnothing(_draw_hook)||_draw_hook(idx,shot,r)
                active=:none;outcomes=Toric3DVertexReplayOutcome[]
                started=time_ns()
                try
                    isnothing(_recovery_hook)||_recovery_hook(idx,shot,e,r,f)
                    for protocol in protocols
                        active=protocol;push!(outcomes,replay_toric3d_vertex_projector(e,r,f;protocol,channel=ch))
                    end
                catch err
                    if err isa Toric3DVertexInvalidRecovery
                        s["invalid_recoveries"]+=1;checkpoint()
                        context=merge(_vertex_boundary_context(e,r,PauliFrame(falses(data_qubit_count(c)),falses(data_qubit_count(c)))),err.context)
                        if !isnothing(path)
                            _vertex_atomic_toml(joinpath(dirname(path),"failure.toml"),Dict("error"=>sprint(showerror,err),"dimensions"=>collect(dims),"channel"=>string(ch),"p"=>p,"protocol"=>string(active),"shot"=>shot,"point_seed"=>_toric3d_feedback_seed_string(seeds[idx]),"context"=>context,"config"=>config))
                        end
                    end
                    rethrow()
                end
                replay=(time_ns()-started)/1e9;values=_vertex_metrics.(outcomes)
                for k in eachindex(_VERTEX_METRICS)
                    mask=sum(Int(values[j][k])<<(j-1) for j in 1:4);s["counts"][k][mask+1]+=1
                end
                s["record_digest"]=_vertex_digest(s["record_digest"],_vertex_record_bytes(r))
                for out in outcomes;s["signature_digest"]=_vertex_digest(s["signature_digest"],_vertex_signature_bytes(out));end
                s["completed_shots"]=shot;s["sampling_seconds"]+=sampling;s["replay_seconds"]+=replay
                (shot%batch_size==0||shot==shots)&&checkpoint()
            end
        end
        metrics=Dict(k=>_vertex_quad_metric(shots,s["counts"][i]) for (i,k) in enumerate(_VERTEX_METRICS))
        push!(points,Toric3DVertexFeedbackPoint(dims,p,ch,shots,seeds[idx],metrics,
            (decoder_build=Float64(s["decoder_build_seconds"]),sampling=Float64(s["sampling_seconds"]),replay=Float64(s["replay_seconds"])),s["invalid_recoveries"],s["record_digest"],s["signature_digest"]))
        partial=Toric3DVertexFeedbackScan(grids,shots,seed,settings,copy(points),config,path,resume)
        isnothing(_point_hook)||_point_hook(partial)
    end
    Toric3DVertexFeedbackScan(grids,shots,seed,settings,points,config,path,resume)
end
