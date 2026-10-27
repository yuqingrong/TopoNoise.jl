"""Counts, individual Wilson intervals, and paired feedback-minus-final statistics."""
struct Toric3DMetricEstimate
    final_only_failures::Int
    layer_feedback_failures::Int
    both_failures::Int
    final_only_only::Int
    layer_feedback_only::Int
    final_only_rate::Float64
    layer_feedback_rate::Float64
    final_only_ci::Tuple{Float64,Float64}
    layer_feedback_ci::Tuple{Float64,Float64}
    difference::Float64
    difference_standard_error::Float64
    difference_ci::Tuple{Float64,Float64}
end

function _toric3d_wilson(failures::Int, shots::Int; z=1.959963984540054)
    rate = failures / shots
    denominator = 1 + z^2 / shots
    center = (rate + z^2 / (2shots)) / denominator
    halfwidth = z * sqrt(rate * (1-rate) / shots + z^2 / (4shots^2)) / denominator
    return (failures == 0 ? 0.0 : max(0.0, center - halfwidth),
            failures == shots ? 1.0 : min(1.0, center + halfwidth))
end

function _toric3d_metric(shots::Int, final::Int, feedback::Int, both::Int)
    0 <= both <= min(final, feedback) && max(final, feedback) <= shots &&
        final + feedback - both <= shots || throw(ArgumentError("invalid paired counts"))
    final_only = final - both
    feedback_only = feedback - both
    difference = (feedback_only - final_only) / shots
    se = shots > 1 ? sqrt(max(0.0,
        ((feedback_only + final_only) / shots - difference^2) / (shots-1))) : NaN
    # Bonferroni combination of two 97.5% Wilson intervals for discordances.
    # Unlike a Wald interval, this retains uncertainty when no pair differs.
    a = _toric3d_wilson(feedback_only, shots; z=2.241402727604947)
    b = _toric3d_wilson(final_only, shots; z=2.241402727604947)
    return Toric3DMetricEstimate(final, feedback, both, final_only, feedback_only,
        final/shots, feedback/shots, _toric3d_wilson(final, shots),
        _toric3d_wilson(feedback, shots), difference, se,
        (max(-1.0, a[1]-b[2]), min(1.0, a[2]-b[1])))
end

"""One noise-channel/size/probability point, containing both paired protocols."""
struct Toric3DComparisonPoint
    dimensions::NTuple{3,Int}
    p::Float64
    channel::Symbol
    shots::Int
    metrics::NamedTuple
    decoder_metadata::NamedTuple
    elapsed_seconds::Float64
end

"""A reproducible scan of preparation faults with ideal layer feedback."""
struct Toric3DComparison
    dimensions::Vector{NTuple{3,Int}}
    error_rates::Vector{Float64}
    channels::Vector{Symbol}
    shots::Int
    seed::Union{Nothing,Int}
    settings::Toric3DDecoderSettings
    x_backend::Symbol
    points::Vector{Toric3DComparisonPoint}
    rng_type::String
    started_at_unix::Float64
end

function _toric3d_channels(values)
    channels = Symbol.(collect(values))
    isempty(channels) && throw(ArgumentError("channels must not be empty"))
    all(channel -> channel in (:x_only, :z_only), channels) ||
        throw(ArgumentError("channels may contain only :x_only and :z_only"))
    length(unique(channels)) == length(channels) ||
        throw(ArgumentError("channels must be unique"))
    return channels
end

"""Lower bound on the logical-Z frame diagnostic, not |0_L> infidelity."""
function toric3d_z_lower_bound(p::Real, Ly::Integer, Lz::Integer)
    _validate_fault_rate(p, "p")
    Ly >= 2 && Lz >= 2 || throw(ArgumentError("Ly and Lz must be at least two"))
    return -expm1(Float64(Ly) * Float64(Lz) * log1p(-2Float64(p))) / 2
end

function _toric3d_sizes(values)
    dimensions = NTuple{3,Int}[]
    for value in values
        dims = if value isa Integer
            (Int(value), Int(value), Int(value))
        elseif value isa Tuple && length(value) == 3 && all(x -> x isa Integer, value)
            Tuple(Int.(value))
        else
            throw(ArgumentError("sizes must contain integers or (Lx,Ly,Lz) tuples"))
        end
        dims[1] >= 1 && dims[2] >= 2 && dims[3] >= 2 ||
            throw(ArgumentError("sizes require Lx >= 1 and Ly,Lz >= 2"))
        push!(dimensions, dims)
    end
    isempty(dimensions) && throw(ArgumentError("sizes must not be empty"))
    length(unique(dimensions)) == length(dimensions) ||
        throw(ArgumentError("sizes must be unique"))
    return dimensions
end

function _toric3d_apply_recovery!(frame, qubits, recovery)
    for (local_qubit, full_qubit) in enumerate(qubits)
        frame.x[full_qubit] ⊻= recovery.x[local_qubit]
        frame.z[full_qubit] ⊻= recovery.z[local_qubit]
    end
    return frame
end

"""
    replay_toric3d_feedback(encoder, faults, decoders; protocol=:layer_feedback,
                           record_history=false)

Replay a sampled fault record using only syndromes for recovery. `decoders[j]`
acts on `encoder.prefixes[j]`; final-only uses the last decoder. History, when
requested, contains pre-correction syndromes and local recovery vectors.
"""
function replay_toric3d_feedback(encoder::Toric3DEncoder, faults::Toric3DFaultRecord,
        decoders::AbstractVector; protocol::Symbol=:layer_feedback,
        record_history::Bool=false)
    protocol in (:final_only, :layer_feedback) ||
        throw(ArgumentError("protocol must be :final_only or :layer_feedback"))
    length(decoders) == length(encoder.checkpoints) ||
        throw(DimensionMismatch("one decoder is required per x-layer"))
    n = data_qubit_count(encoder.code)
    frame = PauliFrame(falses(n), falses(n))
    history = NamedTuple[]
    first_gate = 1
    for j in eachindex(encoder.checkpoints)
        last_gate = encoder.checkpoints[j]
        propagate_toric3d_segment!(frame, encoder, faults, first_gate, last_gate)
        first_gate = last_gate + 1
        if protocol === :layer_feedback || j == length(encoder.checkpoints)
            qubits = encoder.prefix_qubits[j]
            local_frame = PauliFrame(frame.x[qubits], frame.z[qubits])
            syndrome = toric3d_syndrome(encoder.prefixes[j], local_frame)
            recovery = decode_toric3d(decoders[j], syndrome.stars, syndrome.plaquettes)
            _toric3d_apply_recovery!(frame, qubits, recovery)
            record_history && push!(history, (layer=j, syndrome=syndrome, recovery=recovery))
        end
    end
    # Decoder syndrome consistency is checked at every call; this also checks
    # that prefix mappings applied its actual recovery to the intended qubits.
    final_syndrome = toric3d_syndrome(encoder.code, frame)
    !any(final_syndrome.stars) && !any(final_syndrome.plaquettes) ||
        error("3D recovery left a nonzero final syndrome")
    return (frame=frame, logicals=toric3d_logical_parities(encoder.code, frame),
            history=history)
end

_toric3d_failures(logicals) =
    (logicals.logical_x, logicals.logical_z,
     logicals.logical_x | logicals.logical_z, logicals.logical_x)

"""
    compare_toric3d_feedback(rng; sizes=[3,5,7], error_rates=..., shots=1000,
                            seed=nothing, settings=Toric3DDecoderSettings(),
                            x_backend=:bposd, channels=[:x_only,:z_only],
                            output_dir=nothing)

For each size/probability and each of X-only/Z-only noise, replay the same
preparation faults under final-only and layer feedback. `seed` is metadata;
the supplied RNG alone determines draws. Sizes may be cubes or dimension
triples. If `output_dir` is supplied, atomically save all completed points
after each channel. Recoveries use no access to sampled faults.
"""
function compare_toric3d_feedback(rng::Random.AbstractRNG;
        sizes=[3,5,7], error_rates=[0.,.001,.003,.01,.02,.03,.05,.075,.10,.15],
        shots::Integer=1000, seed=nothing,
        settings::Toric3DDecoderSettings=Toric3DDecoderSettings(),
        x_backend::Symbol=:bposd, channels=[:x_only,:z_only],
        output_dir=nothing, progress_io=nothing)
    dimensions = _toric3d_sizes(sizes)
    rates = _scan_error_rates(error_rates)
    selected_channels = _toric3d_channels(channels)
    count_shots = _positive_machine_int(shots, "shots")
    seed_value = _seed_metadata(seed)
    x_backend in (:bposd, :exact) || throw(ArgumentError("x_backend must be :bposd or :exact"))
    x_backend === :exact && any(d -> !(d[1] in (1,2) && d[2] == d[3] == 2), dimensions) &&
        throw(ArgumentError("exact reference supports only (1,2,2) and (2,2,2)"))
    progress_io === nothing || progress_io isa IO ||
        throw(ArgumentError("progress_io must be nothing or an IO"))
    output_dir === nothing || output_dir isa AbstractString ||
        throw(ArgumentError("output_dir must be nothing or a path"))
    comparison = Toric3DComparison(dimensions, rates, selected_channels, count_shots, seed_value,
        settings, x_backend, Toric3DComparisonPoint[], string(typeof(rng)), time())
    for dims in dimensions
        encoder = toric3d_encoder(ToricCode3D(dims...))
        for p in rates
            decoders = [build_toric3d_decoder(prefix; p, x_backend, settings)
                        for prefix in encoder.prefixes]
            metadata = toric3d_decoder_metadata(last(decoders))
            for channel in selected_channels
                started = time()
                final_counts, feedback_counts, both_counts = zeros(Int,4), zeros(Int,4), zeros(Int,4)
                for shot in 1:count_shots
                    faults = sample_toric3d_faults(rng, encoder;
                        p_x=channel === :x_only ? p : 0., p_z=channel === :z_only ? p : 0.)
                    a, b = try
                        (_toric3d_failures(replay_toric3d_feedback(encoder, faults, decoders;
                            protocol=:final_only).logicals),
                         _toric3d_failures(replay_toric3d_feedback(encoder, faults, decoders;
                            protocol=:layer_feedback).logicals))
                    catch error
                        throw(ErrorException("3D point failed at dimensions=$dims, p=$p, " *
                            "channel=$channel, shot=$shot, seed=$seed_value: " * sprint(showerror,error)))
                    end
                    for i in 1:4
                        final_counts[i] += a[i]
                        feedback_counts[i] += b[i]
                        both_counts[i] += a[i] & b[i]
                    end
                end
                metrics = NamedTuple{(:logical_x,:logical_z,:any_logical,:state_failure)}(
                    ntuple(i -> _toric3d_metric(count_shots, final_counts[i], feedback_counts[i], both_counts[i]), 4))
                point = Toric3DComparisonPoint(dims, p, channel, count_shots, metrics, metadata, time()-started)
                push!(comparison.points, point)
                output_dir === nothing || _write_toric3d_checkpoint(comparison, output_dir)
                progress_io === nothing || println(progress_io,
                    "3D ", dims, " ", channel, " p=", p, ": ", count_shots,
                    " paired shots in ", round(point.elapsed_seconds; digits=2), " s")
                progress_io === nothing || flush(progress_io)
            end
        end
    end
    return comparison
end
