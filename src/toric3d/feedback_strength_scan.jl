const _TORIC3D_FEEDBACK_METRICS = (:logical_1,:logical_2,:state_failure)
const _TORIC3D_DIAGONAL_FEEDBACK_CHECKPOINT =
    "toric3d-diagonal-feedback-checkpoint.toml"
const _TORIC3D_DIAGONAL_FEEDBACK_FAILURE =
    "toric3d-diagonal-feedback-failure.toml"

struct Toric3DProtocolEstimate
    failures::Int
    rate::Float64
    ci::Tuple{Float64,Float64}
end

struct Toric3DPairEstimate
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

struct Toric3DTripleMetricEstimate
    shots::Int
    joint_counts::NTuple{8,Int}
    none::Toric3DProtocolEstimate
    slice_local::Toric3DProtocolEstimate
    full_prefix::Toric3DProtocolEstimate
    none_vs_slice_local::Toric3DPairEstimate
    none_vs_full_prefix::Toric3DPairEstimate
    slice_local_vs_full_prefix::Toric3DPairEstimate
end

function _toric3d_protocol_estimate(shots::Int,counts::NTuple{8,Int},bit::Int)
    failures = sum(counts[mask + 1] for mask in 0:7 if (mask & bit) != 0)
    return Toric3DProtocolEstimate(
        failures,failures/shots,_toric3d_wilson(failures,shots))
end

function _toric3d_pair_estimate(shots::Int,counts::NTuple{8,Int},
        first::Symbol,second::Symbol,first_bit::Int,second_bit::Int)
    neither = first_only = second_only = both = 0
    for mask in 0:7
        value = counts[mask + 1]
        a,b = (mask & first_bit) != 0,(mask & second_bit) != 0
        if a && b
            both += value
        elseif a
            first_only += value
        elseif b
            second_only += value
        else
            neither += value
        end
    end
    difference = (second_only-first_only)/shots
    standard_error = shots > 1 ? sqrt(max(0.0,
        ((second_only+first_only)/shots-difference^2)/(shots-1))) : NaN
    second_interval = _toric3d_wilson(second_only,shots;z=2.241402727604947)
    first_interval = _toric3d_wilson(first_only,shots;z=2.241402727604947)
    interval = (max(-1.0,second_interval[1]-first_interval[2]),
        min(1.0,second_interval[2]-first_interval[1]))
    return Toric3DPairEstimate(first,second,neither,first_only,second_only,both,
        difference,standard_error,interval)
end

function _toric3d_triple_metric(shots::Integer,raw_counts)
    count_shots = _positive_machine_int(shots,"shots")
    length(raw_counts) == 8 ||
        throw(ArgumentError("three-protocol joint counts must have eight entries"))
    counts = try
        Tuple(Int.(raw_counts))
    catch error
        throw(ArgumentError("invalid three-protocol joint counts: " * sprint(showerror,error)))
    end
    all(>=(0),counts) || throw(ArgumentError("joint counts must be nonnegative"))
    sum(counts) == count_shots ||
        throw(ArgumentError("joint counts must sum to shots"))
    none = _toric3d_protocol_estimate(count_shots,counts,1)
    slice_local = _toric3d_protocol_estimate(count_shots,counts,2)
    full_prefix = _toric3d_protocol_estimate(count_shots,counts,4)
    return Toric3DTripleMetricEstimate(count_shots,counts,none,slice_local,
        full_prefix,
        _toric3d_pair_estimate(count_shots,counts,:none,:slice_local,1,2),
        _toric3d_pair_estimate(count_shots,counts,:none,:full_prefix,1,4),
        _toric3d_pair_estimate(count_shots,counts,:slice_local,:full_prefix,2,4))
end

struct Toric3DDiagonalFeedbackPoint
    dimensions::NTuple{3,Int}
    p::Float64
    channel::Symbol
    shots::Int
    point_seed::UInt64
    metrics::NamedTuple
    decoder_versions::NamedTuple{(:pymatching,:ldpc),Tuple{String,String}}
    phase_seconds::NamedTuple{(:decoder_build,:sampling,:replay),
        Tuple{Float64,Float64,Float64}}
    invalid_recoveries::Int
end

struct Toric3DDiagonalFeedbackScan
    dimensions::Vector{NTuple{3,Int}}
    error_rates::Vector{Float64}
    channels::Vector{Symbol}
    shots::Int
    seed::Union{Nothing,Int}
    batch_size::Int
    x_routing::Symbol
    settings::Toric3DDecoderSettings
    points::Vector{Toric3DDiagonalFeedbackPoint}
    rng_type::String
    point_seeds::Vector{UInt64}
    started_at_unix::Float64
    resumed::Bool
    checkpoint_path::Union{Nothing,String}
    failure_path::Union{Nothing,String}
    invalid_recoveries::Int
    decoder_builds::NamedTuple{(:x,:z),Tuple{Int,Int}}
end

mutable struct _Toric3DDiagonalFeedbackState
    completed_shots::Int
    joint_counts::Matrix{Int}
    decoder_build_seconds::Float64
    sampling_seconds::Float64
    replay_seconds::Float64
    pymatching_version::String
    ldpc_version::String
    invalid_recoveries::Int
end

_Toric3DDiagonalFeedbackState() = _Toric3DDiagonalFeedbackState(
    0,zeros(Int,3,8),0.0,0.0,0.0,"","",0)

_toric3d_feedback_seed_string(seed::UInt64) = string(seed;base=16,pad=16)
_toric3d_feedback_parse_seed(seed) = parse(UInt64,String(seed);base=16)

function _toric3d_feedback_config(dimensions,rates,channels,shots,seed,
        batch_size,routing,settings,rng)
    return Dict{String,Any}(
        "dimensions"=>[collect(dims) for dims in dimensions],
        "error_rates"=>rates,"channels"=>string.(channels),"shots"=>shots,
        "seed"=>something(seed,"unspecified"),"batch_size"=>batch_size,
        "x_routing"=>string(routing),"rng_type"=>string(typeof(rng)),
        "point_rng_type"=>"Random.Xoshiro",
        "settings"=>Dict(string(key)=>value
            for (key,value) in pairs(_toric3d_settings_metadata(settings))))
end

function _toric3d_feedback_state_dict(state::_Toric3DDiagonalFeedbackState)
    return Dict{String,Any}(
        "completed_shots"=>state.completed_shots,
        "joint_counts"=>[collect(@view state.joint_counts[row,:]) for row in 1:3],
        "decoder_build_seconds"=>state.decoder_build_seconds,
        "sampling_seconds"=>state.sampling_seconds,
        "replay_seconds"=>state.replay_seconds,
        "pymatching_version"=>state.pymatching_version,
        "ldpc_version"=>state.ldpc_version,
        "invalid_recoveries"=>state.invalid_recoveries)
end

function _toric3d_feedback_state(value,shots::Int)
    try
        completed = Int(value["completed_shots"])
        rows = value["joint_counts"]
        length(rows) == 3 && all(length(row)==8 for row in rows) ||
            throw(ArgumentError("invalid checkpoint joint-count shape"))
        counts = reduce(vcat,permutedims.([Int.(row) for row in rows]))
        all(counts .>= 0) || throw(ArgumentError("negative checkpoint joint count"))
        all(sum(@view counts[row,:]) == completed for row in 1:3) ||
            throw(ArgumentError("checkpoint joint totals differ from completed shots"))
        0 <= completed <= shots || throw(ArgumentError("invalid completed shot count"))
        build = Float64(value["decoder_build_seconds"])
        sampling = Float64(value["sampling_seconds"])
        replay = Float64(value["replay_seconds"])
        all(isfinite(x) && x >= 0 for x in (build,sampling,replay)) ||
            throw(ArgumentError("invalid checkpoint phase time"))
        invalid = Int(value["invalid_recoveries"])
        invalid >= 0 || throw(ArgumentError("invalid recovery count is negative"))
        return _Toric3DDiagonalFeedbackState(completed,counts,build,sampling,replay,
            String(value["pymatching_version"]),String(value["ldpc_version"]),invalid)
    catch error
        error isa ArgumentError && rethrow()
        throw(ArgumentError("invalid diagonal-feedback checkpoint state: " *
            sprint(showerror,error)))
    end
end

function _write_toric3d_feedback_checkpoint(path,config,point_seeds,states,started)
    invalid = sum(state.invalid_recoveries for state in states)
    payload = Dict{String,Any}(
        "format_version"=>1,"started_at_unix"=>started,"config"=>config,
        "point_seeds"=>_toric3d_feedback_seed_string.(point_seeds),
        "invalid_recoveries"=>invalid,
        "states"=>_toric3d_feedback_state_dict.(states))
    _toric3d_atomic_write(path) do io
        TOML.print(io,payload;sorted=true)
    end
    return path
end

function _load_toric3d_feedback_checkpoint(path,config,point_count,shots)
    isfile(path) || throw(ArgumentError("resume checkpoint does not exist: $path"))
    saved = try
        TOML.parsefile(path)
    catch error
        throw(ArgumentError("could not parse diagonal-feedback checkpoint: " *
            sprint(showerror,error)))
    end
    get(saved,"format_version",nothing) == 1 ||
        throw(ArgumentError("unsupported diagonal-feedback checkpoint format"))
    get(saved,"config",nothing) == config ||
        throw(ArgumentError("resume configuration does not match the checkpoint"))
    seeds = try
        _toric3d_feedback_parse_seed.(saved["point_seeds"])
    catch error
        throw(ArgumentError("invalid checkpoint point seeds: " * sprint(showerror,error)))
    end
    length(seeds) == point_count || throw(ArgumentError("checkpoint point count differs"))
    raw_states = get(saved,"states",Any[])
    length(raw_states) == point_count || throw(ArgumentError("checkpoint state count differs"))
    states = [_toric3d_feedback_state(value,shots) for value in raw_states]
    sum(state.invalid_recoveries for state in states) ==
        Int(get(saved,"invalid_recoveries",-1)) ||
        throw(ArgumentError("checkpoint invalid-recovery total differs"))
    started = try
        Float64(saved["started_at_unix"])
    catch error
        throw(ArgumentError("invalid checkpoint start time: " * sprint(showerror,error)))
    end
    isfinite(started) || throw(ArgumentError("invalid checkpoint start time"))
    return (point_seeds=seeds,states=states,started_at_unix=started)
end

function _toric3d_feedback_directory(output_dir,resume::Bool)
    if output_dir === nothing
        resume && throw(ArgumentError("resume=true requires output_dir"))
        return nothing
    end
    output_dir isa AbstractString || throw(ArgumentError("output_dir must be a path"))
    isempty(strip(output_dir)) && throw(ArgumentError("output_dir must not be empty"))
    directory = abspath(output_dir)
    isfile(directory) && throw(ArgumentError("output_dir is an existing file"))
    mkpath(directory)
    return directory
end

function _toric3d_feedback_point(dims,p,channel,shots,point_seed,state)
    state.completed_shots == shots || error("cannot finalize an incomplete feedback point")
    metrics = NamedTuple{_TORIC3D_FEEDBACK_METRICS}(ntuple(metric ->
        _toric3d_triple_metric(shots,@view state.joint_counts[metric,:]),3))
    return Toric3DDiagonalFeedbackPoint(dims,p,channel,shots,point_seed,metrics,
        (pymatching=state.pymatching_version,ldpc=state.ldpc_version),
        (decoder_build=state.decoder_build_seconds,sampling=state.sampling_seconds,
            replay=state.replay_seconds),state.invalid_recoveries)
end

function _toric3d_record_scalar_triplet!(counts,triplet)
    flags = ((triplet.none.logical_1,triplet.none.logical_2,triplet.none.failure),
        (triplet.slice_local.logical_1,triplet.slice_local.logical_2,
            triplet.slice_local.failure),
        (triplet.full_prefix.logical_1,triplet.full_prefix.logical_2,
            triplet.full_prefix.failure))
    for metric in 1:3
        mask = Int(flags[1][metric]) + 2Int(flags[2][metric]) + 4Int(flags[3][metric])
        counts[metric,mask+1] += 1
    end
    return counts
end

function _toric3d_record_batch_triplet!(counts,triplet)
    shots = length(triplet.none.failure)
    for shot in 1:shots, metric in 1:3
        field = _TORIC3D_FEEDBACK_METRICS[metric] === :state_failure ? :failure :
            _TORIC3D_FEEDBACK_METRICS[metric]
        mask = Int(getproperty(triplet.none,field)[shot]) +
            2Int(getproperty(triplet.slice_local,field)[shot]) +
            4Int(getproperty(triplet.full_prefix,field)[shot])
        counts[metric,mask+1] += 1
    end
    return counts
end

function _toric3d_test_recovery!(hook,channel,decoders)
    hook === nothing && return nothing
    if channel === :x_only
        decoder = decoders.prefixes[end]
        syndrome = BitVector(decoder.matrix[:,1])
        context = (channel=channel,protocol=:full_prefix,checkpoint=:final,
            recovery_length=size(decoder.matrix,2),syndrome=syndrome)
        recovery = hook(context)
        _toric3d_validate_x_recovery(decoder,syndrome,recovery)
    else
        decoder = decoders.prefixes[end]
        syndrome = BitVector(Vector(decoder.matrix[:,1]) .!= 0)
        context = (channel=channel,protocol=:full_prefix,checkpoint=:final,
            recovery_length=size(decoder.matrix,2),syndrome=syndrome)
        recovery = hook(context)
        _toric3d_validate_z_recovery(decoder,syndrome,recovery)
    end
    return nothing
end

function _write_toric3d_feedback_failure(path,dims,channel,p,protocol,checkpoint,
        batch_start,point_seed,error)
    payload = Dict{String,Any}(
        "dimensions"=>collect(dims),"channel"=>string(channel),"p"=>p,
        "protocol"=>string(protocol),"checkpoint"=>string(checkpoint),
        "batch_start"=>batch_start,
        "point_seed"=>_toric3d_feedback_seed_string(point_seed),
        "error"=>sprint(showerror,error))
    _toric3d_atomic_write(path) do io
        TOML.print(io,payload;sorted=true)
    end
    return path
end

function _toric3d_discard_point_faults!(rng,channel,encoder,p,completed)
    completed == 0 && return nothing
    if channel === :x_only
        for _ in 1:completed
            sample_toric3d_faults(rng,encoder;p_x=p,p_z=0.0)
        end
    else
        sample_toric3d_diagonal_z_fault_batch(rng,encoder;p,shots=completed)
    end
    return nothing
end

"""
Compare no feedback, constant-width slice feedback, and full-prefix feedback
with immutable paired fault records and resumable deterministic point streams.
"""
function compare_toric3d_diagonal_feedback(rng::Random.AbstractRNG;
        sizes=[3,5,7,9],error_rates=[0.,.001,.003,.005,.01,.015,.02],
        channels=[:x_only,:z_only],shots::Integer=1000,seed=nothing,
        batch_size::Integer=256,x_routing::Symbol=:one_step,
        settings::Toric3DDecoderSettings=Toric3DDecoderSettings(),
        output_dir=nothing,resume::Bool=false,progress_io=nothing,
        _test_stop_after_batches=nothing,_test_fault_draw_hook=nothing,
        _test_recovery_hook=nothing)
    dimensions = _toric3d_sizes(sizes)
    rates = _scan_error_rates(error_rates)
    selected_channels = _toric3d_channels(channels)
    count_shots = _positive_machine_int(shots,"shots")
    batch = _positive_machine_int(batch_size,"batch_size")
    seed_value = _seed_metadata(seed)
    x_routing in (:one_step,:full_path) ||
        throw(ArgumentError("x_routing must be :one_step or :full_path"))
    progress_io === nothing || progress_io isa IO ||
        throw(ArgumentError("progress_io must be nothing or an IO"))
    if _test_stop_after_batches !== nothing
        _test_stop_after_batches isa Integer && _test_stop_after_batches > 0 ||
            throw(ArgumentError("_test_stop_after_batches must be positive"))
    end
    _test_fault_draw_hook === nothing || _test_fault_draw_hook isa Function ||
        throw(ArgumentError("_test_fault_draw_hook must be a function"))
    _test_recovery_hook === nothing || _test_recovery_hook isa Function ||
        throw(ArgumentError("_test_recovery_hook must be a function"))

    directory = _toric3d_feedback_directory(output_dir,resume)
    checkpoint_path = directory === nothing ? nothing :
        joinpath(directory,_TORIC3D_DIAGONAL_FEEDBACK_CHECKPOINT)
    failure_path = directory === nothing ? nothing :
        joinpath(directory,_TORIC3D_DIAGONAL_FEEDBACK_FAILURE)
    config = _toric3d_feedback_config(dimensions,rates,selected_channels,
        count_shots,seed_value,batch,x_routing,settings,rng)
    point_count = length(dimensions)*length(rates)*length(selected_channels)
    point_seeds,states,started_at = if resume
        saved = _load_toric3d_feedback_checkpoint(
            checkpoint_path,config,point_count,count_shots)
        (saved.point_seeds,saved.states,saved.started_at_unix)
    else
        seeds = rand(rng,UInt64,point_count)
        initial = [_Toric3DDiagonalFeedbackState() for _ in 1:point_count]
        started = time()
        checkpoint_path === nothing || _write_toric3d_feedback_checkpoint(
            checkpoint_path,config,seeds,initial,started)
        (seeds,initial,started)
    end

    points = Toric3DDiagonalFeedbackPoint[]
    completed_batches = 0
    x_builds = z_builds = 0
    rates_count,channels_count = length(rates),length(selected_channels)
    point_index(di,pi,ci) = ((di-1)*rates_count + pi-1)*channels_count + ci

    for (dimension_index,dims) in enumerate(dimensions)
        code = Toric3DEncodedSlice(dims...)
        encoder = toric3d_encoded_encoder(code)
        z_indices = [point_index(dimension_index,pi,ci)
            for pi in eachindex(rates),ci in eachindex(selected_channels)
            if selected_channels[ci] === :z_only]
        z_family = nothing
        z_build_time = 0.0
        if any(states[index].completed_shots < count_shots for index in z_indices)
            started = time()
            z_family = build_toric3d_z_feedback_decoders(code)
            z_build_time = time()-started
            z_builds += 1
        end

        for (rate_index,p) in enumerate(rates)
            x_indices = [point_index(dimension_index,rate_index,ci)
                for ci in eachindex(selected_channels)
                if selected_channels[ci] === :x_only]
            x_family = nothing
            x_build_time = 0.0
            if any(states[index].completed_shots < count_shots for index in x_indices)
                started = time()
                x_family = build_toric3d_x_feedback_decoders(code;p,settings)
                x_build_time = time()-started
                x_builds += 1
            end

            for (channel_index,channel) in enumerate(selected_channels)
                index = point_index(dimension_index,rate_index,channel_index)
                state = states[index]
                family = channel === :x_only ? x_family : z_family
                if state.completed_shots < count_shots
                    family === nothing && error("missing decoder family for unfinished point")
                    if isempty(state.pymatching_version)
                        state.pymatching_version = channel === :x_only ?
                            family.slice.version : family.seed.version
                        state.ldpc_version = channel === :x_only ?
                            family.prefixes[end].version : ""
                        state.decoder_build_seconds = channel === :x_only ?
                            x_build_time : z_build_time
                    end
                    point_rng = Random.Xoshiro(point_seeds[index])
                    _toric3d_discard_point_faults!(
                        point_rng,channel,encoder,p,state.completed_shots)
                    while state.completed_shots < count_shots
                        batch_start = state.completed_shots + 1
                        batch_count = min(batch,count_shots-state.completed_shots)
                        _test_fault_draw_hook === nothing || _test_fault_draw_hook()
                        sampled_at = time()
                        faults = channel === :x_only ? nothing :
                            sample_toric3d_diagonal_z_fault_batch(
                                point_rng,encoder;p,shots=batch_count)
                        state.sampling_seconds += time()-sampled_at
                        context_protocol,context_checkpoint = :triplet,:unknown
                        replay_at = time()
                        try
                            if channel === :x_only
                                context_protocol,context_checkpoint = :full_prefix,:final
                                _toric3d_test_recovery!(
                                    _test_recovery_hook,channel,family)
                                for _ in 1:batch_count
                                    sampled_at = time()
                                    record = sample_toric3d_faults(
                                        point_rng,encoder;p_x=p,p_z=0.0)
                                    state.sampling_seconds += time()-sampled_at
                                    triplet = toric3d_x_strength_triplet(
                                        encoder,record,family;routing=x_routing)
                                    _toric3d_record_scalar_triplet!(
                                        state.joint_counts,triplet)
                                end
                            else
                                context_protocol,context_checkpoint = :full_prefix,:final
                                _toric3d_test_recovery!(
                                    _test_recovery_hook,channel,family)
                                triplet = toric3d_z_strength_batch(
                                    encoder,faults,family)
                                _toric3d_record_batch_triplet!(
                                    state.joint_counts,triplet)
                            end
                        catch error
                            state.replay_seconds += time()-replay_at
                            state.invalid_recoveries += 1
                            checkpoint_path === nothing ||
                                _write_toric3d_feedback_checkpoint(checkpoint_path,
                                    config,point_seeds,states,started_at)
                            failure_path === nothing ||
                                _write_toric3d_feedback_failure(failure_path,dims,
                                    channel,p,context_protocol,context_checkpoint,
                                    batch_start,point_seeds[index],error)
                            throw(ErrorException("diagonal feedback recovery failed at " *
                                "dimensions=$dims,channel=$channel,p=$p,batch_start=$batch_start," *
                                "point_seed=$(_toric3d_feedback_seed_string(point_seeds[index])): " *
                                sprint(showerror,error)))
                        end
                        state.replay_seconds += time()-replay_at
                        state.completed_shots += batch_count
                        checkpoint_path === nothing ||
                            _write_toric3d_feedback_checkpoint(checkpoint_path,
                                config,point_seeds,states,started_at)
                        completed_batches += 1
                        completed_batches == _test_stop_after_batches &&
                            throw(InterruptException())
                    end
                end
                point = _toric3d_feedback_point(dims,p,channel,count_shots,
                    point_seeds[index],state)
                push!(points,point)
                if progress_io !== nothing
                    println(progress_io,"Diagonal feedback ",dims," p=",p,
                        " channel=",channel,": ",count_shots," paired triplets")
                    flush(progress_io)
                end
            end
        end
    end
    invalid = sum(state.invalid_recoveries for state in states)
    return Toric3DDiagonalFeedbackScan(dimensions,rates,selected_channels,
        count_shots,seed_value,batch,x_routing,settings,points,string(typeof(rng)),
        point_seeds,started_at,resume,checkpoint_path,failure_path,invalid,
        (x=x_builds,z=z_builds))
end
