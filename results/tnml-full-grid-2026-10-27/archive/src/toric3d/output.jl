function _toric3d_atomic_write(writer, path)
    temporary = path * ".tmp"
    open(writer, temporary, "w")
    mv(temporary, path; force=true)
    return path
end

_toric3d_csv(value) = value isa Union{AbstractString,Symbol} ?
    "\"" * replace(string(value), "\"" => "\"\"") * "\"" : string(value)

function _write_toric3d_raw(comparison, path)
    columns = ["Lx","Ly","Lz","qubits","d_Z","d_X","p","channel","metric","shots","seed",
        "final_only_failures","layer_feedback_failures","both_failures",
        "final_only_only","layer_feedback_only","final_only_rate","final_only_ci_low","final_only_ci_high",
        "layer_feedback_rate","layer_feedback_ci_low","layer_feedback_ci_high",
        "difference_feedback_minus_final","difference_standard_error","difference_ci_low","difference_ci_high",
        "logical_z_diagnostic_bound","elapsed_seconds","x_backend","ldpc_version","pymatching_version",
        "max_iter","bp_method","ms_scaling_factor","schedule","omp_thread_count","osd_method","osd_order"]
    _toric3d_atomic_write(path) do io
        println(io, join(columns, ','))
        for point in comparison.points
            Lx,Ly,Lz = point.dimensions
            settings = point.decoder_metadata.settings
            for (name,m) in pairs(point.metrics)
                bound = point.channel === :z_only && name in (:logical_z,:any_logical) ?
                    toric3d_z_lower_bound(point.p,Ly,Lz) : 0.0
                values = (Lx,Ly,Lz,(3Lx+1)*Ly*Lz,Lx+1,Ly*Lz,point.p,point.channel,name,point.shots,
                    something(comparison.seed,"unspecified"),m.final_only_failures,m.layer_feedback_failures,
                    m.both_failures,m.final_only_only,m.layer_feedback_only,m.final_only_rate,m.final_only_ci...,
                    m.layer_feedback_rate,m.layer_feedback_ci...,m.difference,m.difference_standard_error,
                    m.difference_ci...,bound,point.elapsed_seconds,point.decoder_metadata.x_backend,
                    comparison.x_backend === :bposd ? point.decoder_metadata.x_version : "not_used",
                    point.decoder_metadata.z_version,settings.max_iter,settings.bp_method,
                    settings.ms_scaling_factor,settings.schedule,settings.omp_thread_count,
                    settings.osd_method,settings.osd_order)
                println(io, join(_toric3d_csv.(values), ','))
            end
        end
    end
    return path
end

function _toric3d_metadata(comparison)
    data = Dict{String,Any}(
        "format_version" => 1,
        "dimensions" => [collect(d) for d in comparison.dimensions],
        "error_rates" => comparison.error_rates,
        "shots_per_channel_point" => comparison.shots,
        "seed_metadata" => something(comparison.seed,"unspecified"),
        "rng_type" => comparison.rng_type,
        "rng_policy" => "supplied RNG consumed continuously; seed is metadata only",
        "julia_version" => string(VERSION),
        "started_at_unix" => comparison.started_at_unix,
        "completed_channel_points" => length(comparison.points),
        "expected_channel_points" => length(comparison.channels)*length(comparison.dimensions)*length(comparison.error_rates),
        "complete" => length(comparison.points) == length(comparison.channels)*length(comparison.dimensions)*length(comparison.error_rates),
        "protocols" => ["final_only","layer_feedback"],
        "channels" => string.(comparison.channels),
        "geometry" => "Lx retained vertex planes, periodic y/z, rough x ends with triangular plaquettes",
        "gate_order" => "vertices x,y,z increasing; H outgoing x; CNOT incoming x,+y,-y,+z,-z",
        "noise" => "independent Bernoulli p per active output qubit after each H/CNOT",
        "ideal_operations" => ["initialization","idles","syndrome measurement","Pauli recovery"],
        "feedback" => "entire prefix slab, truncated frontier plaquettes, completed stars only",
        "decoder_prior" => "unit MWPM weights; uniform physical-p BP+OSD data priors; circuit correlations neglected",
        "individual_interval" => "95% Wilson score",
        "paired_difference" => "layer_feedback minus final_only, using paired discordant counts",
        "paired_interval" => "Bonferroni combination of two 97.5% Wilson discordance intervals (approximate conservative 95%)",
        "state_failure" => "residual logical X for prepared |0_L>; residual logical Z is only a frame diagnostic",
        "logical_z_bound" => "(1-(1-2p)^(Ly*Lz))/2; Z-only diagnostic, both protocols",
        "decoder_settings" => Dict(string(k)=>v for (k,v) in pairs(_toric3d_settings_metadata(comparison.settings))))
    if !isempty(comparison.points)
        metadata = first(comparison.points).decoder_metadata
        data["decoder_versions"] = Dict("x_backend"=>string(metadata.x_backend),
            "x_version"=>string(metadata.x_version),"z_backend"=>string(metadata.z_backend),
            "z_version"=>string(metadata.z_version))
    end
    return data
end

function _write_toric3d_checkpoint(comparison, directory; basename="toric3d-feedback")
    isempty(strip(directory)) && throw(ArgumentError("output directory must not be empty"))
    !isempty(basename) && Base.basename(basename) == basename ||
        throw(ArgumentError("basename must be a nonempty filename stem"))
    directory = abspath(directory)
    mkpath(directory)
    csv = _write_toric3d_raw(comparison, joinpath(directory, basename * "-raw.csv"))
    metadata = joinpath(directory, basename * "-metadata.toml")
    _toric3d_atomic_write(metadata) do io
        TOML.print(io, _toric3d_metadata(comparison); sorted=true)
    end
    return (csv=csv,metadata=metadata)
end

"""Two channel-matched logical-error panels with paired protocols and Wilson intervals."""
function plot_toric3d_comparison(comparison::Toric3DComparison)
    isempty(comparison.points) && throw(ArgumentError("comparison has no completed points"))
    figure = CairoMakie.Figure(size=(1200,600))
    colors = ["#2364AA","#D1495B","#258477","#8B5FBF","#A66326"]
    handles = Any[]
    labels = String[]
    for (panel,channel,metric,title) in (
            (1,:x_only,:logical_x,"X-only noise · logical X / state failure"),
            (2,:z_only,:logical_z,"Z-only noise · logical Z frame diagnostic"))
        channel_points = filter(p -> p.channel === channel, comparison.points)
        upper = isempty(channel_points) ? 0.6 : min(1.0, max(0.55, 1.08maximum(
            max(getproperty(p.metrics,metric).final_only_ci[2],
                getproperty(p.metrics,metric).layer_feedback_ci[2]) for p in channel_points)))
        axis = CairoMakie.Axis(figure[1,panel]; xlabel="Fault probability p per gate-output qubit",
            ylabel="Residual logical error probability",title=title,
            limits=(nothing,nothing,-.015,upper))
        for (index,dims) in enumerate(comparison.dimensions)
            color = colors[mod1(index,length(colors))]
            points = sort(filter(p -> p.dimensions == dims && p.channel === channel,
                                 comparison.points); by=p->p.p)
            isempty(points) && continue
            x = [p.p for p in points]
            for (protocol,style,marker) in ((:final_only,:solid,:circle),(:layer_feedback,:dash,:rect))
                stats = [getproperty(p.metrics,metric) for p in points]
                y = [getproperty(m,Symbol(protocol,:_rate)) for m in stats]
                ci = [getproperty(m,Symbol(protocol,:_ci)) for m in stats]
                CairoMakie.errorbars!(axis,x,y,[y[i]-ci[i][1] for i in eachindex(y)],
                    [ci[i][2]-y[i] for i in eachindex(y)]; color=(color,.6),whiskerwidth=4)
                CairoMakie.lines!(axis,x,y; color,linestyle=style,linewidth=2)
                CairoMakie.scatter!(axis,x,y; color,marker,markersize=6)
            end
            if channel === :z_only
                grid = range(minimum(x),maximum(x);length=200)
                CairoMakie.lines!(axis,grid,[toric3d_z_lower_bound(p,dims[2],dims[3]) for p in grid];
                    color,linestyle=:dot,linewidth=1.6)
            end
            if panel == 1
                push!(handles,CairoMakie.LineElement(color=color,linewidth=3))
                push!(labels,all(==(dims[1]),dims) ? "L = $(dims[1])" : join(dims," × "))
            end
        end
    end
    append!(handles,[CairoMakie.LineElement(color=:gray25,linestyle=:solid),
        CairoMakie.LineElement(color=:gray25,linestyle=:dash),
        CairoMakie.LineElement(color=:gray25,linestyle=:dot)])
    append!(labels,["Final-only decoding","Feedback after each layer","Analytic Z diagnostic lower bound"])
    CairoMakie.Legend(figure[2,1:2],handles,labels;nbanks=2,framevisible=false,tellheight=true)
    CairoMakie.Label(figure[3,1:2],
        "$(comparison.shots) paired shots / point · ideal measurement and recovery · error bars: 95% Wilson intervals\n" *
        "The Z diagnostic is not |0_L⟩ preparation infidelity. Decoder priors neglect circuit correlations.";
        fontsize=12,tellwidth=false)
    return figure
end

"""Save raw paired statistics, TOML metadata, and PDF/SVG/PNG comparison figures."""
function save_toric3d_comparison(comparison::Toric3DComparison,directory::AbstractString;
        basename::AbstractString="toric3d-feedback")
    paths = _write_toric3d_checkpoint(comparison,directory;basename)
    figure = plot_toric3d_comparison(comparison)
    stem = joinpath(abspath(directory),basename)
    for extension in ("svg","pdf","png")
        CairoMakie.save(stem * "." * extension,figure)
    end
    return merge(paths,(svg=stem*".svg",pdf=stem*".pdf",png=stem*".png"))
end
