export OpenCodeCapacityScanPoint, OpenCodeCapacityScan,
       scan_open_code_capacity, estimate_open_code_crossings

"""Open rectangular patch used for data-edge code-capacity experiments."""
struct OpenCodeCapacityModel
    rows::Int
    cols::Int

    function OpenCodeCapacityModel(rows::Integer, cols::Integer)
        rows >= 2 || throw(ArgumentError("rows must be at least 2"))
        cols >= 2 || throw(ArgumentError("cols must be at least 2"))
        new(Int(rows), Int(cols))
    end
end

"""Independent data-edge errors on an open rectangular patch."""
struct DataEdgeErrors
    horizontal::BitMatrix
    vertical::BitMatrix

    function DataEdgeErrors(horizontal::BitMatrix, vertical::BitMatrix)
        rows, horizontal_cols = size(horizontal)
        vertical_rows, cols = size(vertical)
        rows >= 2 || throw(ArgumentError("data-edge arrays must have at least 2 rows"))
        cols >= 2 || throw(ArgumentError("data-edge arrays must have at least 2 columns"))
        horizontal_cols == cols - 1 ||
            throw(ArgumentError("horizontal errors must have size (R, C - 1)"))
        vertical_rows == rows - 1 ||
            throw(ArgumentError("vertical errors must have size (R - 1, C)"))
        new(horizontal, vertical)
    end
end

"""Sample independent data-edge errors at a common Bernoulli error rate."""
function sample_data_edge_errors(rng, model::OpenCodeCapacityModel; error_rate::Real)
    0 <= error_rate <= 1 || throw(ArgumentError("error_rate must lie in [0, 1]"))
    horizontal = BitMatrix(rand(rng, model.rows, model.cols - 1) .< error_rate)
    vertical = BitMatrix(rand(rng, model.rows - 1, model.cols) .< error_rate)
    return DataEdgeErrors(horizontal, vertical)
end

"""Return the plaquette syndrome induced by independent data-edge errors."""
function code_capacity_syndrome(errors::DataEdgeErrors)
    rows, cols = size(errors.horizontal, 1), size(errors.vertical, 2)
    syndrome = falses(rows - 1, cols - 1)
    for row in 1:(rows - 1), col in 1:(cols - 1)
        syndrome[row, col] = xor(errors.horizontal[row, col], errors.horizontal[row + 1, col],
                                 errors.vertical[row, col], errors.vertical[row, col + 1])
    end
    return syndrome
end

"""Return a central transverse data-edge mask for the requested logical sector."""
function logical_cut(model::OpenCodeCapacityModel; sector::Symbol)
    horizontal = falses(model.rows, model.cols - 1)
    vertical = falses(model.rows - 1, model.cols)
    if sector === :north_south
        horizontal[cld(model.rows, 2), :] .= true
    elseif sector === :east_west
        vertical[:, cld(model.cols, 2)] .= true
    else
        throw(ArgumentError("sector must be :north_south or :east_west"))
    end
    return (horizontal=horizontal, vertical=vertical)
end

"""Return the data-edge residual obtained by XORing errors and correction."""
function _residual_errors(errors::DataEdgeErrors, correction)
    size(errors.horizontal) == size(correction.horizontal) ||
        throw(DimensionMismatch("correction.horizontal must match horizontal errors"))
    size(errors.vertical) == size(correction.vertical) ||
        throw(DimensionMismatch("correction.vertical must match vertical errors"))
    return DataEdgeErrors(
        BitMatrix(xor.(errors.horizontal, correction.horizontal)),
        BitMatrix(xor.(errors.vertical, correction.vertical)))
end

"""Return whether the residual data-edge chain crosses the logical cut oddly."""
function _data_edge_logical_failure(
        errors::DataEdgeErrors, correction;
        sector::Symbol=:north_south)
    model = OpenCodeCapacityModel(
        size(errors.horizontal, 1), size(errors.vertical, 2))
    cut = logical_cut(model; sector=sector)
    residual = _residual_errors(errors, correction)
    crossings = count(residual.horizontal .& cut.horizontal) +
                count(residual.vertical .& cut.vertical)
    return isodd(crossings)
end

"""Streaming logical-failure summary for one open-patch size and error rate."""
struct OpenCodeCapacityScanPoint
    size::Int
    error_rate::Float64
    shots::Int
    logical_failure_ns_mean::Float64
    logical_failure_ns_se::Float64
    logical_failure_ew_mean::Float64
    logical_failure_ew_se::Float64
    logical_failure_ns_batches::Vector{Float64}
    logical_failure_ew_batches::Vector{Float64}
end

"""Finite-size open code-capacity scan on a shared error-rate grid."""
struct OpenCodeCapacityScan
    sizes::Vector{Int}
    error_rates::Vector{Float64}
    points::Vector{OpenCodeCapacityScanPoint}
end

function _open_code_capacity_error_rate(error_rate::Real)
    isfinite(error_rate) && 0 <= error_rate < 0.5 || throw(ArgumentError(
        "error_rate must be finite and lie in [0, 0.5), got $error_rate"))
    return Float64(error_rate)
end

"""
    scan_open_code_capacity(rng, sizes, error_rates; shots=10_000, batches=100,
                            progress_io=nothing)

Stream independently sampled data-edge errors, decode both logical sectors, and
retain equally sized batch means for finite-size crossing bootstraps.
"""
function scan_open_code_capacity(
        rng::Random.AbstractRNG, sizes, error_rates;
        shots::Integer=10_000, batches::Integer=100, progress_io=nothing)
    size_values = Int[value for value in sizes]
    _validate_strictly_increasing(size_values, "sizes")
    all(>=(2), size_values) || throw(ArgumentError(
        "scan sizes must be at least 2"))
    rate_values = Float64[
        _open_code_capacity_error_rate(error_rate) for error_rate in error_rates]
    _validate_strictly_increasing(rate_values, "error_rates")
    shots > 0 || throw(ArgumentError("shots must be positive, got $shots"))
    2 <= batches <= shots || throw(ArgumentError(
        "batches must satisfy 2 <= batches <= shots, got batches=$batches and shots=$shots"))

    shot_count, batch_count = Int(shots), Int(batches)
    points = OpenCodeCapacityScanPoint[]
    for size in size_values
        model = OpenCodeCapacityModel(size, size)
        for error_rate in rate_values
            ns_accumulator = _WelfordAccumulator()
            ew_accumulator = _WelfordAccumulator()
            ns_batch_totals = zeros(batch_count)
            ew_batch_totals = zeros(batch_count)
            batch_counts = zeros(Int, batch_count)
            for shot in 1:shot_count
                errors = sample_data_edge_errors(rng, model; error_rate=error_rate)
                syndrome = code_capacity_syndrome(errors)
                ns = decode_syndrome(
                    model, syndrome; sector=:north_south, error_rate=error_rate)
                ew = decode_syndrome(
                    model, syndrome; sector=:east_west, error_rate=error_rate)
                failure_ns = logical_failure(errors, ns; sector=:north_south)
                failure_ew = logical_failure(errors, ew; sector=:east_west)
                _push!(ns_accumulator, failure_ns)
                _push!(ew_accumulator, failure_ew)
                batch = fld((shot - 1) * batch_count, shot_count) + 1
                ns_batch_totals[batch] += failure_ns
                ew_batch_totals[batch] += failure_ew
                batch_counts[batch] += 1
            end
            push!(points, OpenCodeCapacityScanPoint(
                size, error_rate, shot_count,
                ns_accumulator.mean, _standard_error(ns_accumulator),
                ew_accumulator.mean, _standard_error(ew_accumulator),
                ns_batch_totals ./ batch_counts, ew_batch_totals ./ batch_counts))
            progress_io === nothing || println(
                progress_io,
                "completed L=$size p=$(round(error_rate; digits=4)) ($shot_count shots)")
        end
    end
    return OpenCodeCapacityScan(size_values, rate_values, points)
end

function _open_code_capacity_scan_point(
        scan::OpenCodeCapacityScan, size::Int, error_rate::Float64)
    matches = filter(
        point -> point.size == size && point.error_rate == error_rate,
        scan.points)
    length(matches) == 1 || throw(ArgumentError(
        "scan must contain exactly one point for L=$size, p=$error_rate"))
    return only(matches)
end

"""
    estimate_open_code_crossings(rng, scan; sector=:north_south,
                                 bootstrap=2_000, confidence=0.95)

Estimate adjacent-size logical-failure crossings using monotone curves and a
batch bootstrap.
"""
function estimate_open_code_crossings(
        rng::Random.AbstractRNG, scan::OpenCodeCapacityScan;
        sector::Symbol=:north_south, bootstrap::Integer=2_000,
        confidence::Real=0.95)
    length(scan.sizes) >= 2 || throw(ArgumentError(
        "at least two sizes are required for crossings"))
    length(scan.error_rates) >= 2 || throw(ArgumentError(
        "at least two error rates are required for crossings"))
    bootstrap > 0 || throw(ArgumentError(
        "bootstrap must be positive, got $bootstrap"))
    isfinite(confidence) && 0 < confidence < 1 || throw(ArgumentError(
        "confidence must be strictly between 0 and 1, got $confidence"))
    mean_getter, batch_getter = if sector === :north_south
        (p -> p.logical_failure_ns_mean, p -> p.logical_failure_ns_batches)
    elseif sector === :east_west
        (p -> p.logical_failure_ew_mean, p -> p.logical_failure_ew_batches)
    else
        throw(ArgumentError(
            "sector must be :north_south or :east_west, got $sector"))
    end
    all(length(batch_getter(point)) >= 2 for point in scan.points) ||
        throw(ArgumentError(
            "crossing estimation requires at least two batches per scan point"))

    results = CriticalCrossing[]
    for pair in 1:(length(scan.sizes) - 1)
        small_size, large_size = scan.sizes[pair], scan.sizes[pair + 1]
        small_points = [_open_code_capacity_scan_point(
            scan, small_size, rate) for rate in scan.error_rates]
        large_points = [_open_code_capacity_scan_point(
            scan, large_size, rate) for rate in scan.error_rates]
        small_curve = _isotonic_non_decreasing(
            [mean_getter(point) for point in small_points])
        large_curve = _isotonic_non_decreasing(
            [mean_getter(point) for point in large_points])
        selected = _selected_crossing(
            scan.error_rates, small_curve, large_curve)
        if selected.status != :ok
            push!(results, CriticalCrossing(
                small_size, large_size, missing, missing, missing,
                Float64(confidence), 0.0, selected.status))
            continue
        end
        estimate = selected.estimate

        bootstrap_estimates = Float64[]
        for _ in 1:Int(bootstrap)
            small_sample = _isotonic_non_decreasing([
                _bootstrap_batch_mean(rng, batch_getter(point))
                for point in small_points
            ])
            large_sample = _isotonic_non_decreasing([
                _bootstrap_batch_mean(rng, batch_getter(point))
                for point in large_points
            ])
            sample_selection = _selected_crossing(
                scan.error_rates, small_sample, large_sample)
            sample_selection.status == :ok && push!(
                bootstrap_estimates, sample_selection.estimate)
        end
        valid_fraction = length(bootstrap_estimates) / bootstrap
        if valid_fraction < 0.8
            push!(results, CriticalCrossing(
                small_size, large_size, estimate, missing, missing,
                Float64(confidence), valid_fraction, :unstable))
            continue
        end
        tail = (1 - confidence) / 2
        push!(results, CriticalCrossing(
            small_size, large_size, estimate,
            Statistics.quantile(bootstrap_estimates, tail),
            Statistics.quantile(bootstrap_estimates, 1 - tail),
            Float64(confidence), valid_fraction, :ok))
    end
    return results
end
