const _COMPARISON_SERIES_ORDER = (
    (:as, :x_only), (:as, :z_only), (:bp, :x_only), (:bp, :z_only),
)

const _AS_SERIES_TAG = UInt64(0x243f6a8885a308d3)
const _BP_SERIES_TAG = UInt64(0x13198a2e03707344)
const _X_ONLY_SERIES_TAG = UInt64(0xa4093822299f31d0)
const _Z_ONLY_SERIES_TAG = UInt64(0x082efa98ec4e6c89)

@inline function _mix_uint64(value::UInt64)
    value ⊻= value >> 30
    value *= UInt64(0xbf58476d1ce4e5b9)
    value ⊻= value >> 27
    value *= UInt64(0x94d049bb133111eb)
    return value ⊻ (value >> 31)
end

function _stable_symbol_tag(symbol::Symbol)
    tag = UInt64(0x9e3779b97f4a7c15)
    for byte in codeunits(String(symbol))
        tag = _mix_uint64(tag ⊻ UInt64(byte))
    end
    return tag
end

function _construction_series_tag(construction::Symbol)
    construction === :as && return _AS_SERIES_TAG
    construction === :bp && return _BP_SERIES_TAG
    return _stable_symbol_tag(construction)
end

function _channel_series_tag(channel::Symbol)
    channel === :x_only && return _X_ONLY_SERIES_TAG
    channel === :z_only && return _Z_ONLY_SERIES_TAG
    return _stable_symbol_tag(channel)
end

"""Derive a stable independent RNG seed for one construction/channel key."""
function _comparison_series_seed(
        root_seed::UInt64, key::Tuple{Symbol,Symbol})::UInt64
    construction, channel = key
    return _mix_uint64(root_seed ⊻ _construction_series_tag(construction) ⊻
        _mix_uint64(_channel_series_tag(channel)))
end

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
    # :native selects As |+_L> and Bp |0_L>; each series records its actual state.
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

function _comparison_observable(error_channel::Symbol, failure_metric::Symbol)
    failure_metric === :channel_logical && return _channel_observable(error_channel)
    failure_metric === :state_failure && return :state_failure
    throw(ArgumentError(
        "failure_metric must be :channel_logical or :state_failure"))
end

function _comparison_failure_rate(observable::Symbol, point::LogicalFailurePoint)
    observable === :logical_x && return point.logical_x_failure_rate
    observable === :logical_z && return point.logical_z_failure_rate
    observable === :state_failure && return point.state_failure_rate
    throw(ArgumentError("unsupported comparison observable $(observable)"))
end

function _channel_noise(rate::Float64, channel::Symbol, clock::Symbol)
    channel === :x_only && return CircuitPauliNoise(rate; p_x=rate, p_z=0, clock)
    channel === :z_only && return CircuitPauliNoise(rate; p_x=0, p_z=rate, clock)
    _channel_observable(channel)
end

function channel_failure_count(scan::ChannelFailureScan, point::LogicalFailurePoint)
    scan.logical_observable === :logical_x && return point.logical_x_failures
    scan.logical_observable === :logical_z && return point.logical_z_failures
    scan.logical_observable === :state_failure && return point.state_failures
    throw(ArgumentError("unsupported comparison observable $(scan.logical_observable)"))
end

function channel_failure_rate(scan::ChannelFailureScan, point::LogicalFailurePoint)
    return _comparison_failure_rate(scan.logical_observable, point)
end

function channel_failure_standard_error(scan::ChannelFailureScan, point::LogicalFailurePoint)
    scan.logical_observable === :logical_x && return point.logical_x_standard_error
    scan.logical_observable === :logical_z && return point.logical_z_standard_error
    scan.logical_observable === :state_failure && return point.state_standard_error
    throw(ArgumentError("unsupported comparison observable $(scan.logical_observable)"))
end

"""
    scan_channel_logical_failure(rng; construction, error_channel, kwargs...)

Estimate one reported failure metric for a construction/channel series.
`logical_state=:native` prepares As `|+_L>` from ideal all-plus input or Bp
`|0_L>` from ideal all-zero input, with paired H/CNOT schedules and no logical
basis conversion. `failure_metric=:channel_logical` reports
the channel-matched logical-X/logical-Z diagnostic; `:state_failure` reports
the prepared-state failure rate. The supplied RNG is consumed continuously
across the distance-major, rate-minor scan.
"""
function scan_channel_logical_failure(
        rng::Random.AbstractRNG;
        distances,
        error_rates,
        construction::Symbol=:bp,
        logical_state::Symbol=:native,
        error_channel::Symbol,
        failure_metric::Symbol=:channel_logical,
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
    logical_observable = _comparison_observable(error_channel, failure_metric)
    (progress_io === nothing || progress_io isa IO) ||
        throw(ArgumentError("progress_io must be nothing or an IO stream"))

    validation_code = RotatedPlanarCode(
        first(distance_values); boundary_orientation)
    logical_state = rotated_planar_encoder(
        validation_code; logical_state, construction).logical_state
    _channel_noise(first(error_rate_values), error_channel, clock)

    points = LogicalFailurePoint[]
    sizehint!(points, length(distance_values) * length(error_rate_values))
    for distance_value in distance_values
        code = RotatedPlanarCode(distance_value; boundary_orientation)
        encoder = rotated_planar_encoder(
            code; logical_state, construction)
        for error_rate in error_rate_values
            noise = _channel_noise(error_rate, error_channel, clock)
            point = estimate_logical_failure(
                rng, encoder, noise;
                shots=shot_count, batch_size=chunk_limit, seed=seed_value)
            push!(points, point)
            selected_rate = _comparison_failure_rate(logical_observable, point)
            progress_io === nothing || println(
                progress_io,
                "completed construction=$(construction) " *
                "channel=$(error_channel) state=$(logical_state) " *
                "d=$(point.distance) p=$(error_rate) " *
                "$(logical_observable)=$(selected_rate)")
        end
    end

    return ChannelFailureScan(
        construction, error_channel, logical_observable, logical_state,
        boundary_orientation, clock, distance_values, error_rate_values,
        shot_count, chunk_limit, seed_value, series_seed, points)
end

"""
    run_construction_channel_comparison(rng; distances, error_rates, kwargs...)

Compare As/Bp circuits under X-only and Z-only noise. The default
`logical_state=:native` uses As `|+_L>` and Bp `|0_L>`; an explicit state
prepares that state in all four panels. Each series and point stores its
resolved state. The default metric is the residual logical X or Z error,
selected by the noise channel, even when that Pauli leaves the state unchanged.
"""
function run_construction_channel_comparison(
        rng::Random.AbstractRNG;
        distances,
        error_rates,
        logical_state::Symbol=:native,
        boundary_orientation::Symbol=:x_ns,
        clock::Symbol=:gate_layer,
        failure_metric::Symbol=:channel_logical,
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
            validation_code; logical_state, construction)
        _channel_noise(first(error_rate_values), error_channel, clock)
    end
    _comparison_observable(:x_only, failure_metric)

    root_seed = rand(rng, UInt64)
    series = ChannelFailureScan[]
    sizehint!(series, length(_COMPARISON_SERIES_ORDER))
    for (construction, error_channel) in _COMPARISON_SERIES_ORDER
        series_seed = _comparison_series_seed(
            root_seed, (construction, error_channel))
        series_rng = MersenneTwister(series_seed)
        push!(series, scan_channel_logical_failure(
            series_rng;
            distances=distance_values,
            error_rates=error_rate_values,
            construction,
            logical_state,
            error_channel,
            failure_metric,
            boundary_orientation,
            clock,
            shots=shot_count,
            batch_size=chunk_limit,
            seed=seed_value,
            series_seed,
            progress_io))
    end

    return ConstructionChannelComparison(
        logical_state, boundary_orientation, clock, distance_values, error_rate_values,
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
