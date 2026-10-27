const _VERTEX_BASENAME="toric3d-h5-vertex-projector-feedback"
const _VERTEX_PAIRS=[(a,b) for (i,a) in enumerate(TORIC3D_VERTEX_PROTOCOLS) for b in TORIC3D_VERTEX_PROTOCOLS[i+1:end]]
function _vertex_columns()
    ["T","Lu","Lv","qubits","d_X","d_Z","p","p_eff","channel","metric","shots","seed","point_seed",
     ["joint_mask_$i" for i in 0:15]...,
     ["$(p)_$k" for p in TORIC3D_VERTEX_PROTOCOLS for k in ("failures","rate","ci_low","ci_high")]...,
     ["$(a)_vs_$(b)_$k" for (a,b) in _VERTEX_PAIRS for k in ("neither","first_only","second_only","both","difference","standard_error","ci_low","ci_high")]...,
     "decoder_build_seconds","sampling_seconds","replay_seconds","pymatching_version","ldpc_version","invalid_recoveries","complete","record_digest","signature_digest"]
end
function write_toric3d_vertex_projector_feedback_csv(scan::Toric3DVertexFeedbackScan,path::AbstractString)
    mkpath(dirname(abspath(path)))
    _toric3d_atomic_write(path) do io
        println(io,join(_vertex_columns(),','))
        for pt in scan.points,metric in _VERTEX_METRICS
            t,u,v=pt.dimensions;m=pt.metrics[metric]
            vals=Any[t,u,v,3(t+1)*u*v,(t+1)*min(u,v),2min(u,v),pt.p,clamp(pt.p,1e-12,1-1e-12),pt.channel,metric,pt.shots,something(scan.seed,"unspecified"),_toric3d_feedback_seed_string(pt.point_seed),m.joint_counts...]
            for a in m.protocols;append!(vals,(a.failures,a.rate,a.ci...));end
            for a in m.pairs;append!(vals,(a.neither,a.first_only,a.second_only,a.both,a.difference,a.standard_error,a.ci...));end
            append!(vals,(pt.phase_seconds...,scan.config["backend_versions"]["pymatching"],scan.config["backend_versions"]["ldpc"],pt.invalid_recoveries,true,pt.record_digest,pt.signature_digest))
            length(vals)==length(_vertex_columns()) || error("H+5 CSV schema mismatch")
            println(io,join(_toric3d_csv.(vals),','))
        end
    end
    path
end
function _vertex_write_aggregate(scan,path)
    _toric3d_atomic_write(path) do io
        println(io,"T,Lu,Lv,p,channel,metric,protocol,shots,failures,rate,ci_low,ci_high")
        for pt in scan.points,k in _VERTEX_METRICS,(j,p) in enumerate(TORIC3D_VERTEX_PROTOCOLS)
            a=pt.metrics[k].protocols[j]
            println(io,join(_toric3d_csv.((pt.dimensions...,pt.p,pt.channel,k,p,pt.shots,a.failures,a.rate,a.ci...)),','))
        end
    end
end
function _vertex_metadata(scan,raw)
    expected=sum(length(g.dimensions)*length(g.error_rates) for g in values(scan.channel_grids))
    Dict{String,Any}("format_version"=>1,"circuit"=>"direct-projector H+5","legacy_circuit"=>"plaquette-transfer H+6",
        "fault_locations_per_cell"=>11,"boundary_operations"=>"ideal lower independent white-triad projection; paired ideal upper instrument",
        "boundary_tree"=>"physical-edge-ID BFS excluding union of both logical-X supports; root A_(T-1,0,0)",
        "branch_policy"=>"first M-1 independent clean bits, final even parity; gate faults then branch; same immutable record in all four arms",
        "rng_draw_order"=>"operation, output slot, X then Z Bernoulli, then clean upper branch",
        "matching_signatures"=>"lexicographic sorted (a,b,physical-edge path) after minimum unit physical length; rolling digest includes each ordered checkpoint and final signature",
        "limitations"=>"ideal measurements/feedback; unit-length MWPM and uniform physical-p BP+OSD neglect circuit correlations; paper-motivated circuit, not a supplied 3D reproduction",
        "primary_observables"=>Dict("x_only"=>"logical_x_union","z_only"=>"logical_z_union"),
        "logical_x_label"=>"residual logical X on either pair; also |00_L> preparation failure for X-only noise",
        "logical_z_label"=>"state-independent residual logical Z on either pair; not |00_L> infidelity",
        "complete"=>length(scan.points)==expected,"completed_points"=>length(scan.points),"expected_points"=>expected,
        "invalid_recoveries"=>sum(p.invalid_recoveries for p in scan.points),"raw_csv_sha256"=>_vertex_sha256(raw),
        "source_hashes"=>scan.config["source_hashes"],"decoder_settings"=>scan.config["decoder_settings"],"protocols"=>scan.config["protocols"],
        "backend_versions"=>scan.config["backend_versions"],"seed"=>something(scan.seed,"unspecified"),"shots"=>scan.shots,"config"=>scan.config,
        "bp_priors"=>[Dict("dimensions"=>collect(p.dimensions),"channel"=>string(p.channel),"p"=>p.p,"p_eff"=>clamp(p.p,1e-12,1-1e-12),"instances"=>vcat(["window_$t" for t in 1:p.dimensions[1]-1],["prefix_$t" for t in 1:p.dimensions[1]-1],["final"])) for p in scan.points])
end
function plot_toric3d_vertex_projector_feedback(scan::Toric3DVertexFeedbackScan)
    figs=Dict{Symbol,Any}();colors=[:royalblue,:darkorange,:forestgreen,:crimson,:purple,:sienna];markers=[:circle,:rect,:utriangle,:diamond,:cross,:star5]
    for ch in (:x_only,:z_only)
        haskey(scan.channel_grids,ch)||continue
        metric=ch==:x_only ? :logical_x_union : :logical_z_union
        fig=CairoMakie.Figure(size=(1600,850));title=ch==:x_only ? "X-only: either-pair logical X / |00_L> failure" : "Z-only: either-pair logical Z (state-independent diagnostic)"
        CairoMakie.Label(fig[0,1:4],"Direct-projector H+5 — $title",fontsize=21)
        for j in 1:4;CairoMakie.colsize!(fig.layout,j,CairoMakie.Relative(0.25));end
        axes=[CairoMakie.Axis(fig[1,j],title=string(protocol),xlabel="post-gate p",ylabel="failure probability") for (j,protocol) in enumerate(TORIC3D_VERTEX_PROTOCOLS)]
        deltas=[CairoMakie.Axis(fig[2,j],title="$(TORIC3D_VERTEX_PROTOCOLS[j+1]) − none",xlabel="post-gate p",ylabel="paired rate difference") for j in 1:3]
        for ax in axes;CairoMakie.ylims!(ax,-.025,1.025);end
        for ax in deltas;CairoMakie.hlines!(ax,[0.],color=:gray,linestyle=:dash);end
        for (i,dims) in enumerate(scan.channel_grids[ch].dimensions)
            rows=sort(filter(p->p.channel==ch&&p.dimensions==dims,scan.points),by=p->p.p);isempty(rows)&&continue
            x=[p.p for p in rows];color=colors[mod1(i,6)];marker=markers[mod1(i,6)];label="T,Lu,Lv=$(join(dims,','))"
            for j in 1:4
                vals=[p.metrics[metric].protocols[j] for p in rows];y=[a.rate for a in vals]
                CairoMakie.errorbars!(axes[j],x,y,[a.rate-a.ci[1] for a in vals],[a.ci[2]-a.rate for a in vals];color,whiskerwidth=4)
                CairoMakie.scatterlines!(axes[j],x,y;color,marker,label,markersize=8)
            end
            for j in 1:3
                vals=[p.metrics[metric].pairs[j] for p in rows];y=[a.difference for a in vals]
                CairoMakie.errorbars!(deltas[j],x,y,[a.difference-a.ci[1] for a in vals],[a.ci[2]-a.difference for a in vals];color,whiskerwidth=4)
                CairoMakie.scatterlines!(deltas[j],x,y;color,marker,markersize=8)
            end
        end
        CairoMakie.Legend(fig[2,4],axes[1],"Dimensions (T,Lu,Lv)",framevisible=false)
        CairoMakie.Label(fig[3,1:4],"Wilson 95% rate intervals; paired discordance intervals. Ideal boundaries and feedback. Eleven noisy outputs per cell.",fontsize=14)
        figs[ch]=fig
    end
    figs
end
function save_toric3d_vertex_projector_feedback(scan::Toric3DVertexFeedbackScan,dir::AbstractString;figures::Bool=true,extra_metadata=Dict{String,Any}())
    any(p->p.invalid_recoveries!=0,scan.points)&&error("refusing successful artifact with invalid recovery")
    bound=get(scan.config,"provenance",Dict{String,Any}())
    for (k,v) in bound
        haskey(extra_metadata,k)&&extra_metadata[k]!=v&&throw(ArgumentError("output provenance differs from checkpoint identity"))
    end
    mkpath(dir);stem=joinpath(abspath(dir),_VERTEX_BASENAME)
    raw=write_toric3d_vertex_projector_feedback_csv(scan,stem*"-raw.csv");aggregate=stem*"-aggregate.csv";_vertex_write_aggregate(scan,aggregate)
    metadata=stem*"-metadata.toml";meta=_vertex_metadata(scan,raw)
    isempty(intersect(Set(keys(meta)),Set(keys(extra_metadata)))) || throw(ArgumentError("extra metadata cannot overwrite scientific provenance"))
    merge!(meta,extra_metadata,bound);_vertex_atomic_toml(metadata,meta)
    paths=String[]
    if figures
        for (ch,fig) in plot_toric3d_vertex_projector_feedback(scan),ext in ("pdf","svg","png")
            path=stem*"-$(ch).$ext";CairoMakie.save(path,fig);push!(paths,path)
        end
    end
    (raw=raw,aggregate=aggregate,metadata=metadata,figures=paths)
end
function _vertex_json_value(v)
    if v isa AbstractDict
        d=PythonCall.pybuiltins.dict()
        for (k,value) in v;d[string(k)]=_vertex_json_value(value);end
        return d
    elseif v isa AbstractVector || v isa Tuple
        a=PythonCall.pybuiltins.list()
        for value in v;a.append(_vertex_json_value(value));end
        return a
    end
    PythonCall.Py(v)
end
function _vertex_write_json(path,payload)
    mkpath(dirname(abspath(path)))
    text=PythonCall.pyconvert(String,PythonCall.pyimport("json").dumps(_vertex_json_value(payload);sort_keys=true,indent=2,allow_nan=false))
    _toric3d_atomic_write(path) do io;println(io,text);end
end
struct Toric3DVertexProductionSelection
    data::Dict{String,Any}
    channel_grids::Dict{Symbol,Toric3DVertexScanGrid}
    sha256::String
end
function load_toric3d_vertex_production_selection(path::AbstractString;settings=Toric3DDecoderSettings(),validate_source::Bool=true)
    try
        data=PythonCall.pyconvert(Dict{String,Any},PythonCall.pyimport("json").loads(read(path,String)))
        required=("format_version","rule_version","status","extension_round","scout_inputs","seed","decoder_settings","source_hashes","production_shots","extension_grids","channels")
        all(k->haskey(data,k),required)||throw(ArgumentError("incomplete selection schema"))
        data["format_version"]==1&&data["rule_version"]=="h5-scout-selection-v1"||throw(ArgumentError("unknown selection format/rule"))
        data["status"] in ("ready","needs_extended_scout")&&0<=data["extension_round"]<=3||throw(ArgumentError("invalid selection status/round"))
        data["seed"]==1234&&data["production_shots"]==1000||throw(ArgumentError("selection seed/shots differ from scientific design"))
        data["decoder_settings"]==_vertex_settings(settings)||throw(ArgumentError("selection decoder settings mismatch"))
        !isempty(data["scout_inputs"])||throw(ArgumentError("missing scout ancestry"))
        for input in data["scout_inputs"]
            all(k->haskey(input,k),("csv_sha256","metadata_sha256","path_ids","stage","extension_round","parent_selection_sha256"))||throw(ArgumentError("incomplete scout input"))
            all(k->occursin(r"^[0-9a-f]{64}$",input[k]),("csv_sha256","metadata_sha256"))||throw(ArgumentError("invalid input hash"))
        end
        validate_source&&data["source_hashes"]!=_vertex_source_hashes()&&throw(ArgumentError("selection source hashes mismatch"))
        Set(keys(data["channels"]))==Set(("x_only","z_only"))||throw(ArgumentError("selection must declare both channels"))
        grids=Dict{Symbol,Toric3DVertexScanGrid}();states=String[]
        for (ch,v) in data["channels"]
            all(k->haskey(v,k),("dimensions_T_Lu_Lv","rates","fit_windows","candidate_centers","extension_states","flags","saturation_bounds","timing_model","projected_seconds"))||throw(ArgumentError("incomplete channel selection"))
            for k in ("fit_windows","candidate_centers","extension_states")
                Set(keys(v[k]))==Set(string.(TORIC3D_VERTEX_PROTOCOLS))||throw(ArgumentError("selection must declare each protocol"))
            end
            length(v["dimensions_T_Lu_Lv"])>=4||throw(ArgumentError("production needs four sizes"))
            issorted(v["rates"])||throw(ArgumentError("selection rates must be sorted"))
            for protocol in string.(TORIC3D_VERTEX_PROTOCOLS)
                state=v["extension_states"][protocol];push!(states,state)
                state in ("resolved_candidates","requesting","exhausted_diagnostic")||throw(ArgumentError("invalid arm extension state"))
                wins=v["fit_windows"][protocol];centers=v["candidate_centers"][protocol]
                length(wins)==length(centers)||throw(ArgumentError("candidate/window lengths differ"))
                state=="resolved_candidates"&&isempty(wins)&&throw(ArgumentError("resolved arm has no fit window"))
                for (w,c) in zip(wins,centers)
                    length(w)==2&&0<w[1]<=c<=w[2]<=1||throw(ArgumentError("invalid frozen fit window"))
                    all(p->p in v["rates"],w)||throw(ArgumentError("fit window endpoints absent from grid"))
                end
            end
            !isempty(v["rates"])&&(grids[Symbol(ch)]=Toric3DVertexScanGrid(v["dimensions_T_Lu_Lv"],v["rates"]))
        end
        (data["status"]=="ready")==!("requesting" in states)||throw(ArgumentError("selection status contradicts arm states"))
        data["status"]=="needs_extended_scout"&&isempty(data["extension_grids"])&&throw(ArgumentError("missing extension grids"))
        Toric3DVertexProductionSelection(data,grids,_vertex_sha256(path))
    catch err
        err isa InterruptException&&rethrow()
        err isa ArgumentError&&rethrow()
        throw(ArgumentError("invalid selection: "*sprint(showerror,err)))
    end
end
