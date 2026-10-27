"""Shot-level means and approximate normal intervals for bounded slice diagnostics."""
struct Toric3DSliceMetricEstimate
    unassisted_mean::Float64
    unassisted_standard_error::Float64
    unassisted_ci::Tuple{Float64,Float64}
    slice_feedback_mean::Float64
    slice_feedback_standard_error::Float64
    slice_feedback_ci::Tuple{Float64,Float64}
    difference::Float64
    difference_standard_error::Float64
    difference_ci::Tuple{Float64,Float64}
end

function _toric3d_slice_metric(a::AbstractVector, b::AbstractVector)
    length(a) == length(b) && !isempty(a) ||
        throw(ArgumentError("paired diagnostics require equal nonempty samples"))
    all(x -> isfinite(x) && 0 <= x <= 1, a) &&
        all(x -> isfinite(x) && 0 <= x <= 1, b) ||
        throw(ArgumentError("shot diagnostics must be finite and in [0,1]"))
    n = length(a)
    sem(x) = n > 1 ? std(x) / sqrt(n) : NaN
    interval(m, se, lo, hi) = (max(lo, m-1.959963984540054se), min(hi, m+1.959963984540054se))
    ma, mb, md = mean(a), mean(b), mean(b .- a)
    sa, sb, sd = sem(a), sem(b), sem(b .- a)
    return Toric3DSliceMetricEstimate(ma,sa,interval(ma,sa,0.,1.),
        mb,sb,interval(mb,sb,0.,1.),md,sd,interval(md,sd,-1.,1.))
end

struct Toric3DSliceComparisonPoint
    dimensions::NTuple{3,Int}
    p::Float64
    channel::Symbol
    shots::Int
    metrics::NamedTuple
    elapsed_seconds::Float64
end

struct Toric3DSliceComparison
    dimensions::Vector{NTuple{3,Int}}
    error_rates::Vector{Float64}
    shots::Int
    seed::Union{Nothing,Int}
    routing::Symbol
    points::Vector{Toric3DSliceComparisonPoint}
    rng_type::String
    decoder_version::String
    started_at_unix::Float64
end

function _toric3d_slice_cluster_extent(code::Toric3DSliceLattice,
        syndrome::AbstractVector{Bool})
    T, Ly, Lz = code.steps, code.Ly, code.Lz
    M = Ly*Lz
    length(syndrome) == 3T*M || throw(DimensionMismatch("incorrect plaquette syndrome length"))
    any(syndrome) || return 0.0
    # Two violated faces are adjacent precisely when they bound a common cube.
    parent = collect(1:length(syndrome))
    function root(i)
        while parent[i] != i
            parent[i] = parent[parent[i]]
            i = parent[i]
        end
        return i
    end
    face(t,u,v,pair) = 3*((t-1)*M + mod(u,Ly)*Lz + mod(v,Lz)) + pair
    for t in 2:T, u in 0:Ly-1, v in 0:Lz-1
        faces = (face(t,u,v,1),face(t,u,v,2),face(t,u,v,3),
                 face(t-1,u,v,1),face(t-1,u,v-1,2),face(t-1,u-1,v,3))
        first_root = 0
        for i in faces
            syndrome[i] || continue
            if first_root == 0
                first_root = root(i)
            else
                parent[root(i)] = root(first_root)
            end
        end
    end
    low, high = fill(T+1,length(syndrome)), zeros(Int,length(syndrome))
    for i in eachindex(syndrome)
        syndrome[i] || continue
        component, t = root(i), div(i-1,3M)+1
        low[component] = min(low[component],t)
        high[component] = max(high[component],t)
    end
    return maximum(high .- low .+ 1) / T
end

"""
    toric3d_slice_diagnostics(code, syndromes)

Final physical syndrome densities and maximum longitudinal extent of any
connected violated-plaquette cluster, divided by the number of transfer steps.
Bulk stars exclude the final three-edge frontier stars. Adjacent faces share a
cube; the cluster statistic is zero for an empty syndrome. No final recovery is
performed, and these are not logical-error or threshold estimators.
"""
function toric3d_slice_diagnostics(code::Toric3DSliceLattice, syndromes)
    M, T = code.Ly*code.Lz, code.steps
    length(syndromes.stars) == T*M || throw(DimensionMismatch("incorrect star syndrome length"))
    length(syndromes.plaquettes) == 3T*M || throw(DimensionMismatch("incorrect plaquette syndrome length"))
    length(syndromes.frontier) == M || throw(DimensionMismatch("incorrect frontier syndrome length"))
    bulk_count = (T-1)*M
    return (bulk_star_density=bulk_count == 0 ? 0.0 : count(view(syndromes.stars,1:bulk_count))/bulk_count,
        plaquette_density=count(syndromes.plaquettes)/(3T*M),
        frontier_density=count(syndromes.frontier)/M,
        plaquette_cluster_extent=_toric3d_slice_cluster_extent(code,syndromes.plaquettes))
end

"""
    compare_toric3d_slice_feedback(rng; sizes=[3,5,7], error_rates=..., shots=200,
        seed=nothing, routing=:one_step, output_dir=nothing, progress_io=nothing)

Compare unassisted diagonal-slice transfer with ideal active-slice matching
feed-forward, sharing every sampled H/CNOT output fault between the protocols.
Sizes are cube integers or `(steps,Ly,Lz)` tuples. Separate X-only and Z-only
channels use the existing independent gate-output noise convention. `seed` is
metadata; the supplied RNG determines all samples. Save completed channel
points incrementally when `output_dir` is supplied. No global final decoder is
applied. The boundary conditions prepare a unique state, so the outputs are
physical syndrome/propagation diagnostics, not logical failure probabilities.
"""
function compare_toric3d_slice_feedback(rng::Random.AbstractRNG;
        sizes=[3,5,7],error_rates=[0.,.001,.003,.01,.03,.05],shots::Integer=200,
        seed=nothing,routing::Symbol=:one_step,output_dir=nothing,progress_io=nothing)
    dimensions = _toric3d_sizes(sizes)
    rates = _scan_error_rates(error_rates)
    count_shots = _positive_machine_int(shots,"shots")
    seed_value = _seed_metadata(seed)
    routing in (:one_step,:full_path) || throw(ArgumentError("routing must be :one_step or :full_path"))
    progress_io === nothing || progress_io isa IO || throw(ArgumentError("progress_io must be nothing or an IO"))
    output_dir === nothing || output_dir isa AbstractString || throw(ArgumentError("output_dir must be nothing or a path"))
    version = _toric3d_module_version(_load_pymatching(),"pymatching")
    comparison = Toric3DSliceComparison(dimensions,rates,count_shots,seed_value,routing,
        Toric3DSliceComparisonPoint[],string(typeof(rng)),version,time())
    names = (:bulk_star_density,:plaquette_density,:frontier_density,:plaquette_cluster_extent)
    for dims in dimensions
        code = Toric3DSliceLattice(dims...)
        encoder = toric3d_slice_encoder(code)
        decoder = build_toric3d_slice_decoder(code)
        for p in rates, channel in (:x_only,:z_only)
            started = time()
            a, b = zeros(count_shots,4), zeros(count_shots,4)
            for shot in 1:count_shots
                faults = sample_toric3d_faults(rng,encoder;
                    p_x=channel === :x_only ? p : 0., p_z=channel === :z_only ? p : 0.)
                try
                    unassisted = replay_toric3d_slice(encoder,faults,decoder;feedback=false,routing)
                    assisted = replay_toric3d_slice(encoder,faults,decoder;feedback=true,routing)
                    a[shot,:] .= values(toric3d_slice_diagnostics(code,unassisted.syndromes))
                    b[shot,:] .= values(toric3d_slice_diagnostics(code,assisted.syndromes))
                catch error
                    throw(ErrorException("Slice-transfer point failed at dimensions=$dims, p=$p, " *
                        "channel=$channel, shot=$shot, seed=$seed_value: " * sprint(showerror,error)))
                end
            end
            metrics = NamedTuple{names}(ntuple(i -> _toric3d_slice_metric(view(a,:,i),view(b,:,i)),4))
            point = Toric3DSliceComparisonPoint(dims,p,channel,count_shots,metrics,time()-started)
            push!(comparison.points,point)
            output_dir === nothing || _write_toric3d_slice_checkpoint(comparison,output_dir)
            progress_io === nothing || println(progress_io,"Slice transfer ",dims," ",channel,
                " p=",p,": ",count_shots," paired shots in ",round(point.elapsed_seconds;digits=2)," s")
            progress_io === nothing || flush(progress_io)
        end
    end
    return comparison
end

function _write_toric3d_slice_raw(comparison,path)
    columns = ["steps","Ly","Lz","qubits","p","channel","metric","shots","seed","routing",
        "unassisted_mean","unassisted_standard_error","unassisted_ci_low","unassisted_ci_high",
        "slice_feedback_mean","slice_feedback_standard_error","slice_feedback_ci_low","slice_feedback_ci_high",
        "difference_feedback_minus_unassisted","difference_standard_error","difference_ci_low","difference_ci_high",
        "elapsed_seconds","pymatching_version"]
    _toric3d_atomic_write(path) do io
        println(io,join(columns,','))
        for point in comparison.points
            T,Ly,Lz = point.dimensions
            for (name,m) in pairs(point.metrics)
                entries = (T,Ly,Lz,3Ly*Lz*(T+1),point.p,point.channel,name,point.shots,
                    something(comparison.seed,"unspecified"),comparison.routing,
                    m.unassisted_mean,m.unassisted_standard_error,m.unassisted_ci...,
                    m.slice_feedback_mean,m.slice_feedback_standard_error,m.slice_feedback_ci...,
                    m.difference,m.difference_standard_error,m.difference_ci...,
                    point.elapsed_seconds,comparison.decoder_version)
                println(io,join(_toric3d_csv.(entries),','))
            end
        end
    end
    return path
end

function _toric3d_slice_metadata(comparison)
    expected = 2length(comparison.dimensions)*length(comparison.error_rates)
    return Dict{String,Any}(
        "format_version"=>1,"dimensions_steps_Ly_Lz"=>[collect(d) for d in comparison.dimensions],
        "error_rates"=>comparison.error_rates,"shots_per_channel_point"=>comparison.shots,
        "seed_metadata"=>something(comparison.seed,"unspecified"),"rng_type"=>comparison.rng_type,
        "rng_policy"=>"supplied RNG consumed continuously; identical faults replayed across each protocol pair",
        "julia_version"=>string(VERSION),"started_at_unix"=>comparison.started_at_unix,
        "completed_channel_points"=>length(comparison.points),"expected_channel_points"=>expected,
        "complete"=>length(comparison.points)==expected,
        "protocols"=>["unassisted","slice_feedback"],"channels"=>["x_only","z_only"],
        "routing"=>string(comparison.routing),"decoder"=>"PyMatching unit-weight active-slice compatibility graph",
        "pymatching_version"=>comparison.decoder_version,
        "geometry"=>"diagonal cubic lattice t=x+y+z,u=x,v=y; periodic u,v; slices 0:steps",
        "initial_state"=>"all-zero initial slice and fresh qubits; fixed boundary stabilizers; k=0",
        "gate_order"=>"vertices t,u,v; H(new1), CNOT(new1,new2), CNOT(new1,new3), two old K12 controls into new2, two old K13 controls into new3",
        "gates_per_vertex"=>7,"active_gate_outputs_per_vertex"=>13,
        "noise"=>"independent Bernoulli p per participating qubit after every H/CNOT, X-only and Z-only separately",
        "ideal_operations"=>["initialization","idles","syndrome measurement","Pauli feed-forward"],
        "feedback"=>"measure old-slice Q; validate complete matching recovery; apply routed image F(C) only to fresh targets",
        "final_recovery"=>"none",
        "bulk_star_density"=>"violated six-edge stars / ((steps-1)*Ly*Lz); excludes final frontier stars; zero at steps=1",
        "plaquette_density"=>"violated physical plaquettes / (3*steps*Ly*Lz)",
        "frontier_density"=>"violated final active-slice compatibility checks / (Ly*Lz)",
        "plaquette_cluster_extent"=>"maximum time extent max(t)-min(t)+1 of any violated-face component, divided by steps; shared-cube adjacency; empty=0",
        "individual_interval"=>"approximate 95% normal interval from independent shot-level sample SEM, clipped to [0,1]",
        "paired_interval"=>"approximate 95% normal interval from within-shot feedback-minus-unassisted differences, clipped to [-1,1]",
        "interval_limitations"=>"not Wilson intervals over correlated checks; constant samples give zero SEM; one shot gives undefined SEM; rare-event coverage is not guaranteed",
        "interpretation"=>"propagation and physical-syndrome diagnostics only; no logical error rate or threshold claim; global classical matching is not a constant-depth local channel",
        "paper_relation"=>"derived 3D slice construction inspired by arXiv:2609.13383; retained gate-output noise differs from the paper's slice-noise model")
end

function _write_toric3d_slice_checkpoint(comparison,directory;basename="toric3d-slice-feedback")
    isempty(strip(directory)) && throw(ArgumentError("output directory must not be empty"))
    !isempty(basename) && Base.basename(basename)==basename || throw(ArgumentError("basename must be a nonempty filename stem"))
    directory = abspath(directory)
    mkpath(directory)
    csv = _write_toric3d_slice_raw(comparison,joinpath(directory,basename*"-raw.csv"))
    metadata = joinpath(directory,basename*"-metadata.toml")
    _toric3d_atomic_write(metadata) do io
        TOML.print(io,_toric3d_slice_metadata(comparison);sorted=true)
    end
    return (csv=csv,metadata=metadata)
end

"""Four physical-diagnostic panels; every interval is based on independent shots."""
function plot_toric3d_slice_comparison(comparison::Toric3DSliceComparison)
    isempty(comparison.points) && throw(ArgumentError("comparison has no completed points"))
    figure = CairoMakie.Figure(size=(1200,950))
    colors = ["#2364AA","#D1495B","#258477","#8B5FBF","#A66326"]
    handles,labels = Any[],String[]
    panels = ((1,1,:x_only,:plaquette_density,"X noise · physical plaquette defects","Mean defect fraction"),
        (1,2,:z_only,:bulk_star_density,"Z noise · bulk star defects","Mean defect fraction"),
        (2,1,:x_only,:frontier_density,"X noise · final slice compatibility defects","Mean defect fraction"),
        (2,2,:x_only,:plaquette_cluster_extent,"X noise · plaquette-cluster longitudinal extent","Mean maximum extent / depth"))
    for (row,col,channel,metric,title,ylabel) in panels
        axis = CairoMakie.Axis(figure[row,col];xlabel="Fault probability p per gate-output qubit",ylabel,title)
        upper = .02
        for (index,dims) in enumerate(comparison.dimensions)
            color = colors[mod1(index,length(colors))]
            points = sort(filter(p -> p.dimensions==dims && p.channel===channel,comparison.points);by=p->p.p)
            isempty(points) && continue
            x = [p.p for p in points]
            for (protocol,style,marker) in ((:unassisted,:solid,:circle),(:slice_feedback,:dash,:rect))
                stats = [getproperty(p.metrics,metric) for p in points]
                y = [getproperty(m,Symbol(protocol,:_mean)) for m in stats]
                ci = [getproperty(m,Symbol(protocol,:_ci)) for m in stats]
                finite = findall(i -> all(isfinite,ci[i]),eachindex(y))
                if !isempty(finite)
                    CairoMakie.errorbars!(axis,x[finite],y[finite],[y[i]-ci[i][1] for i in finite],
                        [ci[i][2]-y[i] for i in finite];color=(color,.6),whiskerwidth=4)
                end
                upper = max(upper,maximum(y),isempty(finite) ? 0. : maximum(ci[i][2] for i in finite))
                CairoMakie.lines!(axis,x,y;color,linestyle=style,linewidth=2)
                CairoMakie.scatter!(axis,x,y;color,marker,markersize=6)
            end
            if row==1 && col==1
                push!(handles,CairoMakie.LineElement(color=color,linewidth=3))
                push!(labels,all(==(dims[1]),dims) ? "L = $(dims[1])" : join(dims," × "))
            end
        end
        CairoMakie.ylims!(axis,-.02upper,min(1.02,1.08upper))
    end
    append!(handles,[CairoMakie.LineElement(color=:gray25,linestyle=:solid),
        CairoMakie.LineElement(color=:gray25,linestyle=:dash)])
    append!(labels,["Unassisted transfer","Slice feedback ($(comparison.routing))"])
    CairoMakie.Legend(figure[3,1:2],handles,labels;nbanks=2,framevisible=false,tellheight=true)
    CairoMakie.Label(figure[4,1:2],
        "$(comparison.shots) paired shots / point · error bars: approximate 95% shot-level normal intervals\n" *
        "No final global recovery. These physical diagnostics do not establish a logical-error threshold.";
        fontsize=13,tellwidth=false)
    return figure
end

"""Save incremental statistics and metadata plus the four-panel PDF/SVG/PNG figure."""
function save_toric3d_slice_comparison(comparison::Toric3DSliceComparison,directory::AbstractString;
        basename::AbstractString="toric3d-slice-feedback")
    paths = _write_toric3d_slice_checkpoint(comparison,directory;basename)
    figure = plot_toric3d_slice_comparison(comparison)
    stem = joinpath(abspath(directory),basename)
    for extension in ("svg","pdf","png")
        CairoMakie.save(stem*"."*extension,figure)
    end
    return merge(paths,(svg=stem*".svg",pdf=stem*".pdf",png=stem*".png"))
end
