const _COMPARISON_SERIES_ORDER = (
    (:as, :x_only), (:as, :z_only), (:bp, :x_only), (:bp, :z_only),
)

struct ChannelFailureScan
    construction::Symbol
    error_channel::Symbol
    logical_observable::Symbol
    logical_state::Symbol
    boundary_orientation::Symbol
    clock::Symbol
    distances::Vector{Int}
    error_rates::Vector{Float64}
    shots::Int
    batch_size::Int
    master_seed::Union{Nothing,Int}
    series_seed::UInt64
    points::Vector{LogicalFailurePoint}
end

struct ConstructionChannelComparison
    logical_state::Symbol
    boundary_orientation::Symbol
    clock::Symbol
    distances::Vector{Int}
    error_rates::Vector{Float64}
    shots::Int
    batch_size::Int
    master_seed::Union{Nothing,Int}
    series::Vector{ChannelFailureScan}
end

function _channel_observable(channel::Symbol)
    channel === :x_only && return :logical_x
    channel === :z_only && return :logical_z
    throw(ArgumentError("error_channel must be :x_only or :z_only"))
end

function _channel_noise(rate::Float64, channel::Symbol, clock::Symbol)
    channel === :x_only && return CircuitPauliNoise(rate; p_x=rate, p_z=0, clock)
    channel === :z_only && return CircuitPauliNoise(rate; p_x=0, p_z=rate, clock)
    _channel_observable(channel)
end

channel_failure_count(scan::ChannelFailureScan, point::LogicalFailurePoint) =
    scan.logical_observable === :logical_x ?
    point.logical_x_failures : point.logical_z_failures

channel_failure_rate(scan::ChannelFailureScan, point::LogicalFailurePoint) =
    scan.logical_observable === :logical_x ?
    point.logical_x_failure_rate : point.logical_z_failure_rate

channel_failure_standard_error(scan::ChannelFailureScan, point::LogicalFailurePoint) =
    scan.logical_observable === :logical_x ?
    point.logical_x_standard_error : point.logical_z_standard_error

"""
    scan_channel_logical_failure(rng; construction, error_channel, kwargs...)

Estimate logical failures for one construction/channel series. The logical
state is fixed to `:zero`; the supplied RNG is consumed continuously across
the distance-major, rate-minor scan.
"""
function scan_channel_logical_failure(
        rng::Random.AbstractRNG;
        distances,
        error_rates,
        construction::Symbol=:bp,
        error_channel::Symbol,
        boundary_orientation::Symbol=:x_ns,
        clock::Symbol=:gate_layer,
        shots::Integer=10_000,
        batch_size::Integer=10_000,
        seed=nothing,
        series_seed::UInt64=UInt64(0),
        progress_io=nothing,
    )::ChannelFailureScan
    distance_values = _scan_distances(distances)
    error_rate_values = _scan_error_rates(error_rates)
    shot_count = _positive_machine_int(shots, "shots")
    chunk_limit = _positive_machine_int(batch_size, "batch_size")
    seed_value = _seed_metadata(seed)
    logical_observable = _channel_observable(error_channel)
    (progress_io === nothing || progress_io isa IO) ||
        throw(ArgumentError("progress_io must be nothing or an IO stream"))

    validation_code = RotatedPlanarCode(
        first(distance_values); boundary_orientation)
    rotated_planar_encoder(
        validation_code; logical_state=:zero, construction)
    _channel_noise(first(error_rate_values), error_channel, clock)

    points = LogicalFailurePoint[]
    sizehint!(points, length(distance_values) * length(error_rate_values))
    for distance_value in distance_values
        code = RotatedPlanarCode(distance_value; boundary_orientation)
        encoder = rotated_planar_encoder(
            code; logical_state=:zero, construction)
        for error_rate in error_rate_values
            noise = _channel_noise(error_rate, error_channel, clock)
            point = estimate_logical_failure(
                rng, encoder, noise;
                shots=shot_count, batch_size=chunk_limit, seed=seed_value)
            push!(points, point)
            selected_rate = logical_observable === :logical_x ?
                point.logical_x_failure_rate : point.logical_z_failure_rate
            progress_io === nothing || println(
                progress_io,
                "completed construction=$(construction) " *
                "channel=$(error_channel) d=$(point.distance) p=$(error_rate) " *
                "$(logical_observable)=$(selected_rate)")
        end
    end

    return ChannelFailureScan(
        construction, error_channel, logical_observable, :zero,
        boundary_orientation, clock, distance_values, error_rate_values,
        shot_count, chunk_limit, seed_value, series_seed, points)
end

function run_construction_channel_comparison(
        rng::Random.AbstractRNG;
        distances,
        error_rates,
        boundary_orientation::Symbol=:x_ns,
        clock::Symbol=:gate_layer,
        shots::Integer=10_000,
        batch_size::Integer=10_000,
        seed=nothing,
        progress_io=nothing,
    )::ConstructionChannelComparison
    distance_values = _scan_distances(distances)
    error_rate_values = _scan_error_rates(error_rates)
    shot_count = _positive_machine_int(shots, "shots")
    chunk_limit = _positive_machine_int(batch_size, "batch_size")
    seed_value = _seed_metadata(seed)
    (progress_io === nothing || progress_io isa IO) ||
        throw(ArgumentError("progress_io must be nothing or an IO stream"))

    validation_code = RotatedPlanarCode(
        first(distance_values); boundary_orientation)
    for (construction, error_channel) in _COMPARISON_SERIES_ORDER
        rotated_planar_encoder(
            validation_code; logical_state=:zero, construction)
        _channel_noise(first(error_rate_values), error_channel, clock)
    end

    series = ChannelFailureScan[]
    sizehint!(series, length(_COMPARISON_SERIES_ORDER))
    for (construction, error_channel) in _COMPARISON_SERIES_ORDER
        series_seed = rand(rng, UInt64)
        series_rng = MersenneTwister(series_seed)
        push!(series, scan_channel_logical_failure(
            series_rng;
            distances=distance_values,
            error_rates=error_rate_values,
            construction,
            error_channel,
            boundary_orientation,
            clock,
            shots=shot_count,
            batch_size=chunk_limit,
            seed=seed_value,
            series_seed,
            progress_io))
    end

    return ConstructionChannelComparison(
        :zero, boundary_orientation, clock, distance_values, error_rate_values,
        shot_count, chunk_limit, seed_value, series)
end

function comparison_series(
        comparison::ConstructionChannelComparison,
        construction::Symbol,
        error_channel::Symbol)
    match = findfirst(
        scan -> scan.construction === construction &&
            scan.error_channel === error_channel,
        comparison.series)
    match === nothing && throw(ArgumentError(
        "comparison does not contain construction=$(construction), " *
        "error_channel=$(error_channel)"))
    return comparison.series[match]
end
