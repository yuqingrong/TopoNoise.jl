struct PairCrossing
    lower_distance::Int
    upper_distance::Int
    p::Float64
    standard_error::Float64
end

struct ThresholdFit
    construction::Symbol
    error_channel::Symbol
    logical_observable::Symbol
    status::Symbol
    diagnostic::String
    crossings::Vector{PairCrossing}
    p_c::Union{Nothing,Float64}
    p_c_standard_error::Union{Nothing,Float64}
    nu::Union{Nothing,Float64}
    nu_standard_error::Union{Nothing,Float64}
    bootstrap_replicates::Int
    bootstrap_successes::Int
    bootstrap_seed::UInt64
    exploratory::Bool
end

struct FittedConstructionChannelComparison
    raw::ConstructionChannelComparison
    fits::Vector{ThresholdFit}
end

function _scan_curves(scan::ChannelFailureScan)
    curves = Vector{NamedTuple{(:distance, :rates, :values, :errors),
        Tuple{Int, Vector{Float64}, Vector{Float64}, Vector{Float64}}}}()
    for distance_value in scan.distances
        points = filter(point -> point.distance == distance_value, scan.points)
        length(points) == length(scan.error_rates) ||
            throw(ArgumentError("scan must contain one point per physical rate and distance"))
        values = channel_failure_rate.(Ref(scan), points)
        errors = channel_failure_standard_error.(Ref(scan), points)
        push!(curves, (
            distance=distance_value,
            rates=copy(scan.error_rates),
            values=values,
            errors=errors,
        ))
    end
    return curves
end

function _crossing_candidates(lower, upper)
    candidates = PairCrossing[]
    for index in eachindex(lower.rates)
        delta = lower.values[index] - upper.values[index]
        if delta == 0.0
            combined_error = hypot(lower.errors[index], upper.errors[index])
            push!(candidates, PairCrossing(
                lower.distance, upper.distance, lower.rates[index],
                max(combined_error, eps(Float64))))
        end
    end
    for index in 1:(length(lower.rates) - 1)
        delta = lower.values[index] - upper.values[index]
        delta_next = lower.values[index + 1] - upper.values[index + 1]
        if delta * delta_next < 0.0
            fraction = delta / (delta - delta_next)
            p = lower.rates[index] +
                fraction * (lower.rates[index + 1] - lower.rates[index])
            error = (1 - fraction) * hypot(lower.errors[index], upper.errors[index]) +
                fraction * hypot(lower.errors[index + 1], upper.errors[index + 1])
            push!(candidates, PairCrossing(
                lower.distance, upper.distance, p, max(error, eps(Float64))))
        end
    end
    return candidates
end

"""Return one unambiguous crossing for every adjacent distance pair."""
function _adjacent_crossings(scan::ChannelFailureScan)
    curves = _scan_curves(scan)
    crossings = PairCrossing[]
    for index in 1:(length(curves) - 1)
        candidates = _crossing_candidates(curves[index], curves[index + 1])
        if isempty(candidates)
            return PairCrossing[],
                "missing crossing between d=$(curves[index].distance) and d=$(curves[index + 1].distance)"
        elseif length(candidates) != 1
            return PairCrossing[],
                "ambiguous crossings between d=$(curves[index].distance) and d=$(curves[index + 1].distance)"
        end
        push!(crossings, only(candidates))
    end
    return crossings, nothing
end

function _aggregate_crossings(crossings::Vector{PairCrossing})
    length(crossings) >= 2 || return nothing
    weights = 1.0 ./ getfield.(crossings, :standard_error).^2
    p_c = sum(weights .* getfield.(crossings, :p)) / sum(weights)
    return p_c, sqrt(inv(sum(weights)))
end

function _linear_interpolate(xs::Vector{Float64}, ys::Vector{Float64}, x::Float64)
    (x < first(xs) || x > last(xs)) &&
        throw(ArgumentError("interpolation outside curve domain"))
    index = searchsortedlast(xs, x)
    index == length(xs) && return last(ys)
    left_x, right_x = xs[index], xs[index + 1]
    fraction = (x - left_x) / (right_x - left_x)
    return (1 - fraction) * ys[index] + fraction * ys[index + 1]
end

function _collapse_objective(curves, p_c::Float64, nu::Float64)
    transformed = [
        (
            x=(curve.rates .- p_c) .* curve.distance^(1 / nu),
            y=curve.values,
            se=curve.errors,
        )
        for curve in curves
    ]
    lower = maximum(first(curve.x) for curve in transformed)
    upper = minimum(last(curve.x) for curve in transformed)
    isfinite(lower) && isfinite(upper) && lower < upper || return nothing
    common_x = collect(range(lower, upper; length=64))
    length(common_x) >= 4 || return nothing
    objective = 0.0
    for x in common_x
        for lower_index in 1:(length(transformed) - 1)
            lower_curve = transformed[lower_index]
            lower_y = _linear_interpolate(lower_curve.x, lower_curve.y, x)
            lower_se = _linear_interpolate(lower_curve.x, lower_curve.se, x)
            for upper_index in (lower_index + 1):length(transformed)
                upper_curve = transformed[upper_index]
                upper_y = _linear_interpolate(upper_curve.x, upper_curve.y, x)
                upper_se = _linear_interpolate(upper_curve.x, upper_curve.se, x)
                variance = max(lower_se^2 + upper_se^2, eps(Float64))
                objective += (lower_y - upper_y)^2 / variance
            end
        end
    end
    return isfinite(objective) ? objective : nothing
end

function _fit_nu(scan::ChannelFailureScan, p_c::Float64)
    curves = _scan_curves(scan)
    best_nu = nothing
    best_objective = Inf
    for nu in 0.50:0.01:3.00
        objective = _collapse_objective(curves, p_c, nu)
        objective === nothing && continue
        if objective < best_objective
            best_nu = Float64(nu)
            best_objective = objective
        end
    end
    return best_nu
end

function _unavailable_fit(
        scan::ChannelFailureScan, diagnostic::AbstractString;
        crossings::Vector{PairCrossing}=PairCrossing[],
        bootstrap_replicates::Int=0,
        bootstrap_seed::UInt64=UInt64(0))
    return ThresholdFit(
        scan.construction, scan.error_channel, scan.logical_observable,
        :unavailable, String(diagnostic), crossings, nothing, nothing, nothing,
        nothing, bootstrap_replicates, 0, bootstrap_seed, true)
end

function _raw_threshold_fit(scan::ChannelFailureScan)
    crossings, diagnostic = _adjacent_crossings(scan)
    diagnostic === nothing || return _unavailable_fit(scan, diagnostic)
    aggregate = _aggregate_crossings(crossings)
    aggregate === nothing && return _unavailable_fit(
        scan, "at least two adjacent-pair crossings are required"; crossings)
    p_c, p_c_standard_error = aggregate
    nu = _fit_nu(scan, p_c)
    nu === nothing && return _unavailable_fit(
        scan, "scaling-collapse objective has no finite common domain"; crossings)
    return ThresholdFit(
        scan.construction, scan.error_channel, scan.logical_observable,
        :success, "exploratory crossing and collapse fit", crossings, p_c,
        p_c_standard_error, nu, nothing, 0, 0, UInt64(0), true)
end

function _with_selected_channel_count(
        point::LogicalFailurePoint, observable::Symbol, count::Int)
    rate = count / point.shots
    standard_error = _binomial_standard_error(rate, point.shots)
    if observable === :logical_x
        return LogicalFailurePoint(
            point.distance, point.boundary_orientation, point.construction,
            point.logical_state, point.clock, point.p_x, point.p_z,
            point.shots, point.seed, count, rate, standard_error,
            point.logical_z_failures, point.logical_z_failure_rate,
            point.logical_z_standard_error, point.any_logical_failures,
            point.any_logical_failure_rate, point.any_logical_standard_error,
            point.state_failures, point.state_failure_rate,
            point.state_standard_error)
    end
    return LogicalFailurePoint(
        point.distance, point.boundary_orientation, point.construction,
        point.logical_state, point.clock, point.p_x, point.p_z,
        point.shots, point.seed, point.logical_x_failures,
        point.logical_x_failure_rate, point.logical_x_standard_error,
        count, rate, standard_error, point.any_logical_failures,
        point.any_logical_failure_rate, point.any_logical_standard_error,
        point.state_failures, point.state_failure_rate,
        point.state_standard_error)
end

function _bootstrap_scan(
        rng::Random.AbstractRNG, scan::ChannelFailureScan)
    points = LogicalFailurePoint[]
    sizehint!(points, length(scan.points))
    for point in scan.points
        count = rand(rng, Binomial(point.shots, channel_failure_rate(scan, point)))
        push!(points, _with_selected_channel_count(
            point, scan.logical_observable, count))
    end
    return ChannelFailureScan(
        scan.construction, scan.error_channel, scan.logical_observable,
        scan.logical_state, scan.boundary_orientation, scan.clock,
        copy(scan.distances), copy(scan.error_rates), scan.shots,
        scan.batch_size, scan.master_seed, scan.series_seed, points)
end

"""
    fit_channel_threshold(scan; bootstrap_replicates=500, bootstrap_seed=0)

Return an exploratory threshold crossing and scaling-collapse estimate for one
construction/error-channel scan. An unavailable result is returned whenever
adjacent crossings are missing or ambiguous, or the collapse has no common
interpolation domain.
"""
function fit_channel_threshold(
        scan::ChannelFailureScan;
        bootstrap_replicates::Integer=500,
        bootstrap_seed::Integer=0)::ThresholdFit
    bootstrap_replicates >= 0 ||
        throw(ArgumentError("bootstrap_replicates must be nonnegative"))
    bootstrap_seed >= 0 || throw(ArgumentError("bootstrap_seed must be nonnegative"))
    replicate_count = Int(bootstrap_replicates)
    seed = UInt64(bootstrap_seed)
    raw = _raw_threshold_fit(scan)
    raw.status === :success || return _unavailable_fit(
        scan, raw.diagnostic; crossings=raw.crossings,
        bootstrap_replicates=replicate_count, bootstrap_seed=seed)

    bootstrap_rng = MersenneTwister(seed)
    bootstrap_p_c = Float64[]
    bootstrap_nu = Float64[]
    for _ in 1:replicate_count
        bootstrap = _raw_threshold_fit(_bootstrap_scan(bootstrap_rng, scan))
        bootstrap.status === :success || continue
        push!(bootstrap_p_c, something(bootstrap.p_c))
        push!(bootstrap_nu, something(bootstrap.nu))
    end
    successes = length(bootstrap_p_c)
    p_c_standard_error = successes >= 2 ? std(bootstrap_p_c) : raw.p_c_standard_error
    nu_standard_error = successes >= 2 ? std(bootstrap_nu) : nothing
    return ThresholdFit(
        raw.construction, raw.error_channel, raw.logical_observable,
        raw.status, raw.diagnostic, raw.crossings, raw.p_c,
        p_c_standard_error, raw.nu, nu_standard_error, replicate_count,
        successes, seed, true)
end

"""Fit all stored construction/channel scans in their canonical order."""
function fit_construction_channel_comparison(
        rng::Random.AbstractRNG,
        raw::ConstructionChannelComparison;
        bootstrap_replicates::Integer=500)::FittedConstructionChannelComparison
    fits = ThresholdFit[]
    sizehint!(fits, length(raw.series))
    for scan in raw.series
        push!(fits, fit_channel_threshold(
            scan;
            bootstrap_replicates,
            bootstrap_seed=rand(rng, UInt64)))
    end
    return FittedConstructionChannelComparison(raw, fits)
end

"""Rescale a physical error rate about an exploratory threshold estimate."""
function scaled_error_rate(
        error_rate::Real, distance_value::Integer, fit::ThresholdFit)
    fit.status === :success && fit.p_c !== nothing && fit.nu !== nothing ||
        throw(ArgumentError("scaled error rates require a successful threshold fit"))
    distance_value > 0 || throw(ArgumentError("distance must be positive"))
    return (Float64(error_rate) - fit.p_c) *
        Float64(distance_value)^(1 / fit.nu)
end
