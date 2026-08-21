struct IsometricPlanarCapacityPoint
    distance::Int
    error_rate::Float64
    shots::Int
    logical_failure_count::Int
    logical_failure_rate::Float64
    logical_failure_se::Float64
    logical_failure_batches::Vector{Float64}
    diagnostic_only::Bool
    seed::Union{Missing,Int}
end

"""Finite-size scan of decoded planar virtual-bond capacity curves."""
struct IsometricPlanarCapacityScan
    distances::Vector{Int}
    error_rates::Vector{Float64}
    points::Vector{IsometricPlanarCapacityPoint}
end

function _validate_isometric_planar_scan_counts(shots::Integer, batches::Integer)
    shots > 0 || throw(ArgumentError("shots must be positive, got $shots"))
    2 <= batches <= shots || throw(ArgumentError(
        "batches must satisfy 2 <= batches <= shots, got batches=$batches and shots=$shots"))
    shots % batches == 0 || throw(ArgumentError(
        "shots must be divisible by batches, got shots=$shots and batches=$batches"))
    return Int(shots), Int(batches)
end

function _isometric_planar_result_seed(seed)
    (seed === nothing || ismissing(seed)) && return missing
    try
        return Int(seed)
    catch
        throw(ArgumentError("seed must be convertible to Int, got $seed"))
    end
end

"""
    estimate_isometric_planar_capacity(rng, model, error_rate;
                                       shots=10_000, batches=100, seed=nothing)

Sample internal virtual-X frames, calculate exact detector parities, decode
the full batch in PyMatching, and compare its logical prediction with the
held-out west-boundary frame parity.
"""
function estimate_isometric_planar_capacity(
        rng::Random.AbstractRNG, model::IsometricPlanarCodeModel, error_rate::Real;
        shots::Integer=10_000, batches::Integer=100, seed=nothing)
    rate = _isometric_planar_error_rate(error_rate)
    shot_count, batch_count = _validate_isometric_planar_scan_counts(shots, batches)
    batch_size = div(shot_count, batch_count)
    batch_rates = Vector{Float64}(undef, batch_count)
    failure_count = 0
    for batch in 1:batch_count
        detectors, actual = _isometric_planar_detector_batch(
            rng, model, rate, batch_size)
        predicted = _isometric_planar_batch_predictions(model, detectors)
        failures = predicted .!= actual
        batch_failures = count(failures)
        failure_count += batch_failures
        batch_rates[batch] = batch_failures / batch_size
    end
    failure_rate = failure_count / shot_count
    failure_se = failure_rate == 0 || failure_rate == 1 ? 0.0 :
        sqrt(failure_rate * (1 - failure_rate) / shot_count)
    return IsometricPlanarCapacityPoint(
        model.distance, rate, shot_count, failure_count, failure_rate, failure_se,
        batch_rates, model.distance == 2, _isometric_planar_result_seed(seed))
end

function _isometric_planar_capacity_point(
        scan::IsometricPlanarCapacityScan, distance::Int, error_rate::Float64)
    matches = filter(point -> point.distance == distance &&
                              point.error_rate == error_rate, scan.points)
    length(matches) == 1 || throw(ArgumentError(
        "scan must contain exactly one point for d=$distance, p=$error_rate"))
    return only(matches)
end

"""Scan directly specified distances on a shared virtual-bond error-rate grid."""
function scan_isometric_planar_capacity(
        rng::Random.AbstractRNG, distances, error_rates;
        shots::Integer=10_000, batches::Integer=100, seed=nothing, progress_io=nothing)
    distance_values = Int[value for value in distances]
    rate_values = [_isometric_planar_error_rate(value) for value in error_rates]
    _validate_strictly_increasing(distance_values, "distances")
    all(distance -> distance >= 2, distance_values) || throw(ArgumentError(
        "distances must be at least 2"))
    _validate_strictly_increasing(rate_values, "error_rates")
    shot_count, _ = _validate_isometric_planar_scan_counts(shots, batches)
    result_seed = _isometric_planar_result_seed(seed)
    models = Dict(distance => IsometricPlanarCodeModel(distance)
                  for distance in distance_values)
    points = IsometricPlanarCapacityPoint[]
    for distance in distance_values, rate in rate_values
        point = estimate_isometric_planar_capacity(
            rng, models[distance], rate; shots=shots, batches=batches,
            seed=result_seed)
        push!(points, point)
        progress_io === nothing || println(
            progress_io, "completed d=$distance p=$rate ($shot_count shots)")
    end
    return IsometricPlanarCapacityScan(distance_values, rate_values, points)
end

function _unique_crossing(
        rates::Vector{Float64}, first_curve::Vector{Float64},
        second_curve::Vector{Float64})
    differences = first_curve .- second_curve
    all(iszero, differences) && return (estimate=missing, status=:unstable)
    candidates = Float64[]

    for index in 1:(length(rates) - 1)
        left, right = differences[index], differences[index + 1]
        (iszero(left) || iszero(right)) && continue
        left * right < 0 || continue
        fraction = -left / (right - left)
        push!(candidates,
            rates[index] + fraction * (rates[index + 1] - rates[index]))
    end

    ambiguous_plateau = false
    index = 1
    while index <= length(rates)
        if !iszero(differences[index])
            index += 1
            continue
        end
        first_zero = index
        while index < length(rates) && iszero(differences[index + 1])
            index += 1
        end
        last_zero = index
        if first_zero > 1 && last_zero < length(rates)
            left = differences[first_zero - 1]
            right = differences[last_zero + 1]
            if left * right < 0
                first_zero == last_zero ?
                    push!(candidates, rates[first_zero]) :
                    (ambiguous_plateau = true)
            end
        end
        index += 1
    end

    ambiguous_plateau && return (estimate=missing, status=:unstable)
    isempty(candidates) && return (estimate=missing, status=:unbracketed)
    length(candidates) == 1 || return (estimate=missing, status=:unstable)
    return (estimate=only(candidates), status=:ok)
end

"""Bootstrap adjacent-distance crossings, rejecting multiple intersections."""
function estimate_isometric_planar_crossings(
        rng::Random.AbstractRNG, scan::IsometricPlanarCapacityScan;
        bootstrap::Integer=2_000, confidence::Real=0.95)
    length(scan.distances) >= 2 || throw(ArgumentError(
        "at least two distances are required for crossings"))
    length(scan.error_rates) >= 2 || throw(ArgumentError(
        "at least two error rates are required for crossings"))
    bootstrap > 0 || throw(ArgumentError("bootstrap must be positive, got $bootstrap"))
    isfinite(confidence) && 0 < confidence < 1 || throw(ArgumentError(
        "confidence must be strictly between 0 and 1, got $confidence"))
    results = CriticalCrossing[]
    for index in 1:(length(scan.distances) - 1)
        small, large = scan.distances[index], scan.distances[index + 1]
        small_points = [_isometric_planar_capacity_point(scan, small, rate)
                        for rate in scan.error_rates]
        large_points = [_isometric_planar_capacity_point(scan, large, rate)
                        for rate in scan.error_rates]
        small_curve = _isotonic_non_decreasing(
            [point.logical_failure_rate for point in small_points])
        large_curve = _isotonic_non_decreasing(
            [point.logical_failure_rate for point in large_points])
        selected = _unique_crossing(scan.error_rates, small_curve, large_curve)
        if selected.status != :ok
            push!(results, CriticalCrossing(
                small, large, missing, missing, missing, Float64(confidence), 0.0,
                selected.status))
            continue
        end
        estimates = Float64[]
        for _ in 1:Int(bootstrap)
            sampled_small = _isotonic_non_decreasing([
                _bootstrap_batch_mean(rng, point.logical_failure_batches)
                for point in small_points])
            sampled_large = _isotonic_non_decreasing([
                _bootstrap_batch_mean(rng, point.logical_failure_batches)
                for point in large_points])
            candidate = _unique_crossing(
                scan.error_rates, sampled_small, sampled_large)
            candidate.status == :ok && push!(estimates, candidate.estimate)
        end
        valid_fraction = length(estimates) / bootstrap
        if valid_fraction < 0.8
            push!(results, CriticalCrossing(
                small, large, selected.estimate, missing, missing,
                Float64(confidence), valid_fraction, :unstable))
            continue
        end
        tail = (1 - confidence) / 2
        push!(results, CriticalCrossing(
            small, large, selected.estimate,
            Statistics.quantile(estimates, tail),
            Statistics.quantile(estimates, 1 - tail), Float64(confidence),
            valid_fraction, :ok))
    end
    return results
end
