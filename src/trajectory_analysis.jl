"""Raw observables derived from one measured toric-code trajectory."""
struct TrajectoryObservables
    sampled_error_density::Float64
    boundary_error_density::Float64
    mismatch_density::Float64
    frustration_density::Float64
    largest_cluster_fraction::Float64
    spans_horizontal::Bool
    spans_vertical::Bool
end

"""Streaming summary for one square size and one virtual-bond error rate."""
struct TrajectoryScanPoint
    size::Int
    error_rate::Float64
    shots::Int
    sampled_error_mean::Float64
    sampled_error_se::Float64
    boundary_error_mean::Float64
    boundary_error_se::Float64
    mismatch_mean::Float64
    mismatch_se::Float64
    frustration_mean::Float64
    frustration_se::Float64
    largest_cluster_mean::Float64
    largest_cluster_se::Float64
    horizontal_spanning_mean::Float64
    horizontal_spanning_se::Float64
    vertical_spanning_mean::Float64
    vertical_spanning_se::Float64
    horizontal_span_batches::Vector{Float64}
end

"""Finite-size trajectory scan on a shared error-rate grid."""
struct TrajectoryScan
    sizes::Vector{Int}
    error_rates::Vector{Float64}
    points::Vector{TrajectoryScanPoint}
end

"""One adjacent-size crossing estimate and its bootstrap interval."""
struct CriticalCrossing
    small_size::Int
    large_size::Int
    estimate::Union{Missing,Float64}
    ci_low::Union{Missing,Float64}
    ci_high::Union{Missing,Float64}
    confidence::Float64
    valid_bootstrap_fraction::Float64
    status::Symbol
end

function _find_root!(parent::Vector{Int}, vertex::Int)
    while parent[vertex] != vertex
        parent[vertex] = parent[parent[vertex]]
        vertex = parent[vertex]
    end
    return vertex
end

function _union_vertices!(
        parent::Vector{Int}, component_size::Vector{Int},
        first_vertex::Int, second_vertex::Int)
    first_root = _find_root!(parent, first_vertex)
    second_root = _find_root!(parent, second_vertex)
    first_root == second_root && return first_root
    if component_size[first_root] < component_size[second_root]
        first_root, second_root = second_root, first_root
    end
    parent[second_root] = first_root
    component_size[first_root] += component_size[second_root]
    return first_root
end

function _cluster_observables(horizontal::BitMatrix, vertical::BitMatrix)
    rows = size(horizontal, 1)
    cols = size(vertical, 2)
    vertices = rows * cols
    parent = collect(1:vertices)
    component_size = ones(Int, vertices)
    vertex(row, col) = (row - 1) * cols + col

    for row in 1:rows, col in axes(horizontal, 2)
        horizontal[row, col] || continue
        _union_vertices!(
            parent, component_size, vertex(row, col), vertex(row, col + 1))
    end
    for row in axes(vertical, 1), col in 1:cols
        vertical[row, col] || continue
        _union_vertices!(
            parent, component_size, vertex(row, col), vertex(row + 1, col))
    end

    touches_west = falses(vertices)
    touches_east = falses(vertices)
    touches_north = falses(vertices)
    touches_south = falses(vertices)
    for row in 1:rows, col in 1:cols
        root = _find_root!(parent, vertex(row, col))
        touches_west[root] |= col == 1
        touches_east[root] |= col == cols
        touches_north[root] |= row == 1
        touches_south[root] |= row == rows
    end
    roots = unique(_find_root!(parent, vertex) for vertex in 1:vertices)
    largest = maximum(component_size[root] for root in roots) / vertices
    spans_horizontal = cols > 1 && any(
        touches_west[root] && touches_east[root] for root in roots)
    spans_vertical = rows > 1 && any(
        touches_north[root] && touches_south[root] for root in roots)
    return largest, spans_horizontal, spans_vertical
end

function _boundary_error_count(errors::VirtualBondErrors)
    return count(errors.west_boundary) + count(errors.east_boundary) +
           count(errors.south_boundary) + count(errors.north_boundary)
end

"""Compute raw error, frustration, cluster, and spanning observables."""
function trajectory_observables(trajectory::ToricCodeTrajectory)
    rows, cols, _ = size(trajectory.measurements)
    mismatches = bond_mismatches(trajectory)
    mismatches.horizontal == trajectory.errors.horizontal_internal ||
        throw(ArgumentError(
            "horizontal physical-record mismatches disagree with stored errors"))
    mismatches.vertical == trajectory.errors.vertical_internal ||
        throw(ArgumentError(
            "vertical physical-record mismatches disagree with stored errors"))

    internal_bonds = length(mismatches.horizontal) + length(mismatches.vertical)
    mismatch_count = count(mismatches.horizontal) + count(mismatches.vertical)
    mismatch_density = internal_bonds == 0 ? 0.0 : mismatch_count / internal_bonds

    boundary_bonds = 2 * rows + 2 * cols
    boundary_error_density =
        _boundary_error_count(trajectory.errors) / boundary_bonds
    sampled_error_density = count(trajectory.errors) / length(trajectory.errors)

    plaquettes = (rows - 1) * (cols - 1)
    frustrated = 0
    for row in 1:(rows - 1), col in 1:(cols - 1)
        frustrated += xor(
            mismatches.horizontal[row, col],
            mismatches.horizontal[row + 1, col],
            mismatches.vertical[row, col],
            mismatches.vertical[row, col + 1])
    end
    frustration_density = plaquettes == 0 ? 0.0 : frustrated / plaquettes
    largest, spans_horizontal, spans_vertical =
        _cluster_observables(mismatches.horizontal, mismatches.vertical)

    return TrajectoryObservables(
        sampled_error_density, boundary_error_density, mismatch_density,
        frustration_density, largest, spans_horizontal, spans_vertical)
end

function _validate_strictly_increasing(values, name::AbstractString)
    isempty(values) && throw(ArgumentError("$name must not be empty"))
    all(values[index] < values[index + 1] for index in 1:(length(values) - 1)) ||
        throw(ArgumentError("$name must be sorted with no duplicates"))
    return values
end

function _observable_values(observables::TrajectoryObservables)
    return (
        observables.sampled_error_density,
        observables.boundary_error_density,
        observables.mismatch_density,
        observables.frustration_density,
        observables.largest_cluster_fraction,
        Float64(observables.spans_horizontal),
        Float64(observables.spans_vertical),
    )
end

mutable struct _WelfordAccumulator
    count::Int
    mean::Float64
    squared_deviations::Float64
end

_WelfordAccumulator() = _WelfordAccumulator(0, 0.0, 0.0)

function _push!(accumulator::_WelfordAccumulator, value::Real)
    accumulator.count += 1
    delta = Float64(value) - accumulator.mean
    accumulator.mean += delta / accumulator.count
    accumulator.squared_deviations +=
        delta * (Float64(value) - accumulator.mean)
    return accumulator
end

function _standard_error(accumulator::_WelfordAccumulator)
    accumulator.count <= 1 && return 0.0
    variance = accumulator.squared_deviations / (accumulator.count - 1)
    return sqrt(max(0.0, variance) / accumulator.count)
end

"""
    scan_trajectories(rng, sizes, error_rates; shots=10_000, batches=100)

Stream raw trajectory summaries for square patches. Horizontal-spanning batch
means are retained for later bootstrap crossings.
"""
function scan_trajectories(
        rng::Random.AbstractRNG, sizes, error_rates;
        shots::Integer=10_000, batches::Integer=100)
    size_values = Int[values for values in sizes]
    _validate_strictly_increasing(size_values, "sizes")
    all(>=(2), size_values) || throw(ArgumentError(
        "scan sizes must be at least 2"))

    rate_values = Float64[_validate_error_rate(rate) for rate in error_rates]
    _validate_strictly_increasing(rate_values, "error_rates")
    shots > 0 || throw(ArgumentError("shots must be positive, got $shots"))
    batches > 0 || throw(ArgumentError(
        "batches must be positive, got $batches"))
    batches <= shots || throw(ArgumentError(
        "batches must not exceed shots, got batches=$batches and shots=$shots"))
    shot_count, batch_count = Int(shots), Int(batches)

    points = TrajectoryScanPoint[]
    for size in size_values
        model = ToricCodeTrajectoryModel(size, size)
        for error_rate in rate_values
            accumulators = [_WelfordAccumulator() for _ in 1:7]
            horizontal_batch_totals = zeros(batch_count)
            batch_counts = zeros(Int, batch_count)
            for shot in 1:shot_count
                trajectory = sample_trajectory(
                    rng, model; error_rate=error_rate)
                values = _observable_values(trajectory_observables(trajectory))
                for index in eachindex(values)
                    _push!(accumulators[index], values[index])
                end
                batch = fld((shot - 1) * batch_count, shot_count) + 1
                horizontal_batch_totals[batch] += values[6]
                batch_counts[batch] += 1
            end
            means = [accumulator.mean for accumulator in accumulators]
            standard_errors = _standard_error.(accumulators)
            horizontal_batches =
                horizontal_batch_totals ./ batch_counts
            push!(points, TrajectoryScanPoint(
                size, error_rate, shot_count,
                means[1], standard_errors[1],
                means[2], standard_errors[2],
                means[3], standard_errors[3],
                means[4], standard_errors[4],
                means[5], standard_errors[5],
                means[6], standard_errors[6],
                means[7], standard_errors[7],
                horizontal_batches))
        end
    end
    return TrajectoryScan(size_values, rate_values, points)
end

function _scan_point(scan::TrajectoryScan, size::Int, error_rate::Float64)
    matches = filter(
        point -> point.size == size && point.error_rate == error_rate,
        scan.points)
    length(matches) == 1 || throw(ArgumentError(
        "scan must contain exactly one point for L=$size, p=$error_rate"))
    return only(matches)
end

function _isotonic_non_decreasing(values::Vector{Float64})
    levels = Float64[]
    widths = Int[]
    for value in values
        push!(levels, value)
        push!(widths, 1)
        while length(levels) >= 2 && levels[end - 1] > levels[end]
            combined_width = widths[end - 1] + widths[end]
            combined_level = (
                levels[end - 1] * widths[end - 1] +
                levels[end] * widths[end]) / combined_width
            pop!(levels)
            pop!(widths)
            levels[end] = combined_level
            widths[end] = combined_width
        end
    end
    result = Float64[]
    for (level, width) in zip(levels, widths)
        append!(result, fill(level, width))
    end
    return result
end

function _selected_crossing(
        rates::Vector{Float64}, first_curve::Vector{Float64},
        second_curve::Vector{Float64})
    candidates = Tuple{Float64,Float64}[]
    differences = first_curve .- second_curve
    all(iszero, differences) && return (estimate=missing, status=:unstable)

    for index in 1:(length(rates) - 1)
        left, right = differences[index], differences[index + 1]
        (iszero(left) || iszero(right)) && continue
        left * right < 0 || continue
        fraction = -left / (right - left)
        crossing = rates[index] + fraction * (rates[index + 1] - rates[index])
        first_value = first_curve[index] +
                      fraction * (first_curve[index + 1] - first_curve[index])
        second_value = second_curve[index] +
                       fraction * (second_curve[index + 1] - second_curve[index])
        push!(candidates, (crossing, (first_value + second_value) / 2))
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
                if first_zero == last_zero
                    push!(candidates, (
                        rates[first_zero],
                        (first_curve[first_zero] +
                         second_curve[first_zero]) / 2))
                else
                    ambiguous_plateau = true
                end
            end
        end
        index += 1
    end
    ambiguous_plateau && return (estimate=missing, status=:unstable)
    isempty(candidates) && return (estimate=missing, status=:unbracketed)
    best = argmin(abs(candidate[2] - 0.5) for candidate in candidates)
    return (estimate=candidates[best][1], status=:ok)
end

function _bootstrap_batch_mean(rng::Random.AbstractRNG, batches::Vector{Float64})
    isempty(batches) && throw(ArgumentError(
        "horizontal spanning batches must not be empty"))
    return sum(batches[rand(rng, eachindex(batches))] for _ in eachindex(batches)) /
           length(batches)
end

"""Estimate adjacent-size horizontal-spanning crossings with batch bootstrap."""
function estimate_crossings(
        rng::Random.AbstractRNG, scan::TrajectoryScan;
        bootstrap::Integer=2_000, confidence::Real=0.95)
    length(scan.sizes) >= 2 || throw(ArgumentError(
        "at least two sizes are required for crossings"))
    length(scan.error_rates) >= 2 || throw(ArgumentError(
        "at least two error rates are required for crossings"))
    bootstrap > 0 || throw(ArgumentError(
        "bootstrap must be positive, got $bootstrap"))
    isfinite(confidence) && 0 < confidence < 1 || throw(ArgumentError(
        "confidence must be strictly between 0 and 1, got $confidence"))
    all(length(point.horizontal_span_batches) >= 2 for point in scan.points) ||
        throw(ArgumentError(
            "crossing estimation requires at least two batches per scan point"))

    results = CriticalCrossing[]
    for pair in 1:(length(scan.sizes) - 1)
        small_size, large_size = scan.sizes[pair], scan.sizes[pair + 1]
        small_points = [
            _scan_point(scan, small_size, rate) for rate in scan.error_rates]
        large_points = [
            _scan_point(scan, large_size, rate) for rate in scan.error_rates]
        small_curve = _isotonic_non_decreasing(
            [point.horizontal_spanning_mean for point in small_points])
        large_curve = _isotonic_non_decreasing(
            [point.horizontal_spanning_mean for point in large_points])
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
                _bootstrap_batch_mean(rng, point.horizontal_span_batches)
                for point in small_points
            ])
            large_sample = _isotonic_non_decreasing([
                _bootstrap_batch_mean(rng, point.horizontal_span_batches)
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
