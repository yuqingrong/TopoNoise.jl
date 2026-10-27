const _TORIC3D_X_FEEDBACK_METRICS = (:logical_1, :logical_2, :state_failure)
const _TORIC3D_X_FEEDBACK_CHECKPOINT = "toric3d-x-feedback-checkpoint.toml"

"""One completed size/probability point from the paired X-feedback scan."""
struct Toric3DXFeedbackPoint
    dimensions::NTuple{3,Int}
    p::Float64
    shots::Int
    point_seed::UInt64
    routing::Symbol
    metrics::NamedTuple
    decoder_version::String
    elapsed_seconds::Float64
end

"""Configuration, provenance, and completed points for an encoded X scan."""
struct Toric3DXFeedbackScan
    dimensions::Vector{NTuple{3,Int}}
    error_rates::Vector{Float64}
    shots::Int
    seed::Union{Nothing,Int}
    batch_size::Int
    routing::Symbol
    settings::Toric3DDecoderSettings
    points::Vector{Toric3DXFeedbackPoint}
    rng_type::String
    point_seeds::Vector{UInt64}
    started_at_unix::Float64
    resumed::Bool
    checkpoint_path::Union{Nothing,String}
end

mutable struct _Toric3DXFeedbackState
    completed_shots::Int
    final_counts::Vector{Int}
    feedback_counts::Vector{Int}
    both_counts::Vector{Int}
    elapsed_seconds::Float64
    decoder_version::String
end

_Toric3DXFeedbackState() = _Toric3DXFeedbackState(0, zeros(Int, 3),
    zeros(Int, 3), zeros(Int, 3), 0.0, "")

_toric3d_x_feedback_seed_string(seed::UInt64) = string(seed; base=16, pad=16)
_toric3d_x_feedback_parse_seed(seed) = parse(UInt64, String(seed); base=16)

function _toric3d_x_feedback_config(dimensions, rates, shots, seed, batch_size,
        routing, settings, rng)
    return Dict{String,Any}(
        "dimensions" => [collect(dims) for dims in dimensions],
        "error_rates" => rates,
        "shots" => shots,
        "seed" => something(seed, "unspecified"),
        "batch_size" => batch_size,
        "routing" => string(routing),
        "rng_type" => string(typeof(rng)),
        "point_rng_type" => "Random.Xoshiro",
        "settings" => Dict(string(key) => value
            for (key, value) in pairs(_toric3d_settings_metadata(settings))))
end

function _toric3d_x_feedback_state_dict(state::_Toric3DXFeedbackState)
    return Dict{String,Any}(
        "completed_shots" => state.completed_shots,
        "final_counts" => state.final_counts,
        "feedback_counts" => state.feedback_counts,
        "both_counts" => state.both_counts,
        "elapsed_seconds" => state.elapsed_seconds,
        "decoder_version" => state.decoder_version)
end

function _toric3d_x_feedback_state(value, shots::Int)
    try
        completed = Int(value["completed_shots"])
        final = Int.(value["final_counts"])
        feedback = Int.(value["feedback_counts"])
        both = Int.(value["both_counts"])
        elapsed = Float64(value["elapsed_seconds"])
        version = String(value["decoder_version"])
        0 <= completed <= shots || throw(ArgumentError("invalid completed shot count"))
        length(final) == length(feedback) == length(both) == 3 ||
            throw(ArgumentError("invalid checkpoint metric count"))
        all(0 .<= both .<= min.(final, feedback)) &&
            all(max.(final, feedback) .<= completed) ||
            throw(ArgumentError("invalid checkpoint paired counts"))
        isfinite(elapsed) && elapsed >= 0 ||
            throw(ArgumentError("invalid checkpoint elapsed time"))
        return _Toric3DXFeedbackState(completed, final, feedback, both, elapsed, version)
    catch error
        error isa ArgumentError && rethrow()
        throw(ArgumentError("invalid X-feedback checkpoint state: " * sprint(showerror, error)))
    end
end

function _write_toric3d_x_feedback_checkpoint(path::AbstractString, config,
        point_seeds, states, started_at_unix)
    payload = Dict{String,Any}(
        "format_version" => 1,
        "started_at_unix" => started_at_unix,
        "config" => config,
        "point_seeds" => _toric3d_x_feedback_seed_string.(point_seeds),
        "states" => _toric3d_x_feedback_state_dict.(states))
    _toric3d_atomic_write(path) do io
        TOML.print(io, payload; sorted=true)
    end
    return path
end

function _load_toric3d_x_feedback_checkpoint(path::AbstractString, config,
        point_count::Int, shots::Int)
    isfile(path) || throw(ArgumentError("resume checkpoint does not exist: $path"))
    saved = try
        TOML.parsefile(path)
    catch error
        throw(ArgumentError("could not parse X-feedback checkpoint: " * sprint(showerror, error)))
    end
    get(saved, "format_version", nothing) == 1 ||
        throw(ArgumentError("unsupported X-feedback checkpoint format"))
    get(saved, "config", nothing) == config ||
        throw(ArgumentError("resume configuration does not match the checkpoint"))
    seeds = try
        _toric3d_x_feedback_parse_seed.(saved["point_seeds"])
    catch error
        throw(ArgumentError("invalid checkpoint point seeds: " * sprint(showerror, error)))
    end
    length(seeds) == point_count || throw(ArgumentError("checkpoint point count differs"))
    raw_states = get(saved, "states", Any[])
    length(raw_states) == point_count || throw(ArgumentError("checkpoint state count differs"))
    states = [_toric3d_x_feedback_state(value, shots) for value in raw_states]
    started = try
        Float64(saved["started_at_unix"])
    catch error
        throw(ArgumentError("invalid checkpoint start time: " * sprint(showerror, error)))
    end
    isfinite(started) || throw(ArgumentError("invalid checkpoint start time"))
    return (point_seeds=seeds, states=states, started_at_unix=started)
end

function _toric3d_x_feedback_point(dimensions, p, shots, point_seed, routing,
        state::_Toric3DXFeedbackState)
    state.completed_shots == shots || error("cannot finalize an incomplete X-feedback point")
    metrics = NamedTuple{_TORIC3D_X_FEEDBACK_METRICS}(ntuple(index ->
        _toric3d_metric(shots, state.final_counts[index],
            state.feedback_counts[index], state.both_counts[index]), 3))
    return Toric3DXFeedbackPoint(dimensions, p, shots, point_seed, routing,
        metrics, state.decoder_version, state.elapsed_seconds)
end

function _toric3d_x_feedback_directory(output_dir, resume::Bool)
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

"""
    compare_toric3d_x_feedback_threshold(rng; kwargs...)

Compare final-only recovery with ideal within-slice feedback under paired
post-gate X faults. Optional batch checkpoints preserve exact point substreams.
"""
function compare_toric3d_x_feedback_threshold(rng::Random.AbstractRNG;
        sizes=[3,5,7,9],
        error_rates=[0.,.0001,.0003,.001,.002,.003,.005,.0075,.01,.015,.02],
        shots::Integer=1000, seed=nothing, batch_size::Integer=256,
        routing::Symbol=:one_step,
        settings::Toric3DDecoderSettings=Toric3DDecoderSettings(),
        output_dir=nothing, resume::Bool=false, progress_io=nothing,
        _test_stop_after_batches=nothing)
    dimensions = _toric3d_sizes(sizes)
    rates = _scan_error_rates(error_rates)
    count_shots = _positive_machine_int(shots, "shots")
    batch = _positive_machine_int(batch_size, "batch_size")
    seed_value = _seed_metadata(seed)
    routing in (:one_step, :full_path) ||
        throw(ArgumentError("routing must be :one_step or :full_path"))
    progress_io === nothing || progress_io isa IO ||
        throw(ArgumentError("progress_io must be nothing or an IO"))
    if _test_stop_after_batches !== nothing
        _test_stop_after_batches isa Integer && _test_stop_after_batches > 0 ||
            throw(ArgumentError("_test_stop_after_batches must be a positive integer"))
    end
    directory = _toric3d_x_feedback_directory(output_dir, resume)
    checkpoint_path = directory === nothing ? nothing :
        joinpath(directory, _TORIC3D_X_FEEDBACK_CHECKPOINT)
    config = _toric3d_x_feedback_config(dimensions, rates, count_shots,
        seed_value, batch, routing, settings, rng)
    point_count = length(dimensions) * length(rates)
    point_seeds, states, started_at = if resume
        saved = _load_toric3d_x_feedback_checkpoint(checkpoint_path, config,
            point_count, count_shots)
        (saved.point_seeds, saved.states, saved.started_at_unix)
    else
        seeds = rand(rng, UInt64, point_count)
        initial_states = [_Toric3DXFeedbackState() for _ in 1:point_count]
        started = time()
        checkpoint_path === nothing || _write_toric3d_x_feedback_checkpoint(
            checkpoint_path, config, seeds, initial_states, started)
        (seeds, initial_states, started)
    end

    scan = Toric3DXFeedbackScan(dimensions, rates, count_shots, seed_value,
        batch, routing, settings, Toric3DXFeedbackPoint[], string(typeof(rng)),
        point_seeds, started_at, resume, checkpoint_path)
    point_index = 0
    completed_batches = 0
    for dims in dimensions
        code = Toric3DEncodedSlice(dims...)
        encoder = toric3d_encoded_encoder(code)
        slice_decoder = build_toric3d_slice_decoder(code.lattice)
        for p in rates
            point_index += 1
            state = states[point_index]
            if state.completed_shots < count_shots
                final_decoder = build_toric3d_x_threshold_decoder(code; p, settings)
                state.decoder_version = final_decoder.version
                point_rng = Random.Xoshiro(point_seeds[point_index])
                for _ in 1:state.completed_shots
                    sample_toric3d_faults(point_rng, encoder; p_x=p, p_z=0.0)
                end
                while state.completed_shots < count_shots
                    batch_started = time()
                    batch_count = min(batch, count_shots - state.completed_shots)
                    for offset in 1:batch_count
                        shot = state.completed_shots + offset
                        faults = sample_toric3d_faults(point_rng, encoder;
                            p_x=p, p_z=0.0)
                        pair = try
                            toric3d_x_feedback_pair(encoder, faults, slice_decoder,
                                final_decoder; routing)
                        catch error
                            throw(ErrorException("X-feedback recovery failed at " *
                                "dimensions=$dims,p=$p,routing=$routing,shot=$shot," *
                                "point_seed=$(_toric3d_x_feedback_seed_string(point_seeds[point_index])): " *
                                sprint(showerror, error)))
                        end
                        final_flags = (pair.final_only.logical_1,
                            pair.final_only.logical_2, pair.final_only.failure)
                        feedback_flags = (pair.slice_feedback.logical_1,
                            pair.slice_feedback.logical_2, pair.slice_feedback.failure)
                        for metric in 1:3
                            state.final_counts[metric] += final_flags[metric]
                            state.feedback_counts[metric] += feedback_flags[metric]
                            state.both_counts[metric] +=
                                final_flags[metric] & feedback_flags[metric]
                        end
                    end
                    state.completed_shots += batch_count
                    state.elapsed_seconds += time() - batch_started
                    checkpoint_path === nothing || _write_toric3d_x_feedback_checkpoint(
                        checkpoint_path, config, point_seeds, states, started_at)
                    completed_batches += 1
                    completed_batches == _test_stop_after_batches && throw(InterruptException())
                end
            end
            point = _toric3d_x_feedback_point(dims, p, count_shots,
                point_seeds[point_index], routing, state)
            push!(scan.points, point)
            if progress_io !== nothing
                println(progress_io, "X feedback ", dims, " p=", p, ": ",
                    count_shots, " paired shots in ",
                    round(point.elapsed_seconds; digits=2), " s")
                flush(progress_io)
            end
        end
    end
    return scan
end
