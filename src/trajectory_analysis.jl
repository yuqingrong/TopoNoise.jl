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

"""Equal-weight spin moments after marginalizing mismatched internal edges."""
struct MarginalSpinObservables
    component_count::Int
    second_moment::Float64
    fourth_moment::Float64
    absolute_magnetization::Float64
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

"""Marginalized spin summary for one square size and one error rate."""
struct MarginalSpinScanPoint
    size::Int
    error_rate::Float64
    shots::Int
    absolute_magnetization_mean::Float64
    absolute_magnetization_se::Float64
    second_moment_mean::Float64
    second_moment_se::Float64
    fourth_moment_mean::Float64
    fourth_moment_se::Float64
    binder_cumulant::Float64
    binder_cumulant_se::Union{Missing,Float64}
    batch_counts::Vector{Int}
    second_moment_batches::Vector{Float64}
    fourth_moment_batches::Vector{Float64}
end

"""Finite-size trajectory scan on a shared error-rate grid."""
struct TrajectoryScan
    sizes::Vector{Int}
    error_rates::Vector{Float64}
    points::Vector{TrajectoryScanPoint}
    marginal_spin_points::Vector{MarginalSpinScanPoint}
end

TrajectoryScan(
    sizes::Vector{Int}, error_rates::Vector{Float64},
    points::Vector{TrajectoryScanPoint}) =
    TrajectoryScan(sizes, error_rates, points, MarginalSpinScanPoint[])

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

"""
    bond_mismatches(trajectory)

Return the doubled-edge XOR records as `(horizontal, vertical)` bit matrices.
Horizontal entries compare a site's east record with its right neighbor's west
record. Vertical entries compare a lower site's north record with its upper
neighbor's south record.
"""
function bond_mismatches(trajectory::ToricCodeTrajectory)
    rows, cols, directions = size(trajectory.measurements)
    directions == 4 || throw(DimensionMismatch(
        "trajectory measurements must have four E/N/W/S entries per site, " *
        "got size $(size(trajectory.measurements))"))
    _validate_virtual_error_shapes(rows, cols, trajectory.errors)

    horizontal = BitMatrix(undef, rows, cols - 1)
    for row in 1:rows, col in 1:(cols - 1)
        horizontal[row, col] = xor(
            trajectory.measurements[row, col, 1],
            trajectory.measurements[row, col + 1, 3])
    end
    vertical = BitMatrix(undef, rows - 1, cols)
    for row in 1:(rows - 1), col in 1:cols
        vertical[row, col] = xor(
            trajectory.measurements[row + 1, col, 2],
            trajectory.measurements[row, col, 4])
    end
    return (; horizontal, vertical)
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

function _matched_component_sizes(trajectory::ToricCodeTrajectory)
    mismatches = bond_mismatches(trajectory)
    rows = size(mismatches.horizontal, 1)
    cols = size(mismatches.vertical, 2)
    vertices = rows * cols
    parent = collect(1:vertices)
    component_size = ones(Int, vertices)
    vertex(row, col) = (row - 1) * cols + col

    for row in 1:rows, col in axes(mismatches.horizontal, 2)
        mismatches.horizontal[row, col] && continue
        _union_vertices!(
            parent, component_size, vertex(row, col), vertex(row, col + 1))
    end
    for row in axes(mismatches.vertical, 1), col in 1:cols
        mismatches.vertical[row, col] && continue
        _union_vertices!(
            parent, component_size, vertex(row, col), vertex(row + 1, col))
    end

    roots = unique(_find_root!(parent, index) for index in 1:vertices)
    return [component_size[root] for root in roots]
end

"""
    marginal_spin_observables(rng, trajectory; spin_samples=1)

Treat every mismatched internal doubled edge as an equal-weight erasure and
every matched edge as an equal-spin constraint. Return exact conditional
second and fourth magnetization moments together with a Monte Carlo estimate
of the conditional absolute magnetization.
"""
function marginal_spin_observables(
        rng::Random.AbstractRNG, trajectory::ToricCodeTrajectory;
        spin_samples::Integer=1)
    spin_samples > 0 || throw(ArgumentError(
        "spin_samples must be positive, got $spin_samples"))
    component_sizes = _matched_component_sizes(trajectory)
    site_count = sum(component_sizes)
    squared_size_sum = sum(Float64(size)^2 for size in component_sizes)
    fourth_size_sum = sum(Float64(size)^4 for size in component_sizes)
    second_moment = squared_size_sum / site_count^2
    fourth_moment = (
        3 * squared_size_sum^2 - 2 * fourth_size_sum) / site_count^4

    absolute_total = 0.0
    for _ in 1:Int(spin_samples)
        signed_size_sum = sum(
            rand(rng, Bool) ? size : -size for size in component_sizes)
        absolute_total += abs(signed_size_sum) / site_count
    end
    return MarginalSpinObservables(
        length(component_sizes), second_moment, fourth_moment,
        absolute_total / spin_samples)
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

function _binder_cumulant(second_moment::Real, fourth_moment::Real)
    second_moment > 0 || throw(ArgumentError(
        "second magnetization moment must be positive, got $second_moment"))
    return 1 - fourth_moment / (3 * second_moment^2)
end

function _binder_jackknife_se(
        batch_counts::Vector{Int}, second_batches::Vector{Float64},
        fourth_batches::Vector{Float64})
    batch_count = length(batch_counts)
    batch_count < 2 && return missing
    length(second_batches) == batch_count || throw(DimensionMismatch(
        "second-moment batch count does not match batch_counts"))
    length(fourth_batches) == batch_count || throw(DimensionMismatch(
        "fourth-moment batch count does not match batch_counts"))

    total_count = sum(batch_counts)
    second_total = sum(batch_counts .* second_batches)
    fourth_total = sum(batch_counts .* fourth_batches)
    estimates = Float64[]
    for batch in 1:batch_count
        retained_count = total_count - batch_counts[batch]
        retained_count > 0 || throw(ArgumentError(
            "jackknife requires at least one retained trajectory"))
        second = (
            second_total - batch_counts[batch] * second_batches[batch]) /
            retained_count
        fourth = (
            fourth_total - batch_counts[batch] * fourth_batches[batch]) /
            retained_count
        push!(estimates, _binder_cumulant(second, fourth))
    end
    center = Statistics.mean(estimates)
    return sqrt((batch_count - 1) / batch_count *
                sum((estimate - center)^2 for estimate in estimates))
end

"""
    scan_trajectories(
        rng, sizes, error_rates; shots=10_000, batches=100, spin_samples=1)

Stream raw trajectory summaries and equal-weight marginalized spin moments for
square patches. Horizontal-spanning and paired spin-moment batches are retained
for later bootstrap crossings.
"""
function scan_trajectories(
        rng::Random.AbstractRNG, sizes, error_rates;
        shots::Integer=10_000, batches::Integer=100,
        spin_samples::Integer=1)
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
    spin_samples > 0 || throw(ArgumentError(
        "spin_samples must be positive, got $spin_samples"))
    shot_count, batch_count = Int(shots), Int(batches)
    spin_sample_count = Int(spin_samples)
    spin_rng = Random.Xoshiro(rand(rng, UInt64))

    points = TrajectoryScanPoint[]
    marginal_spin_points = MarginalSpinScanPoint[]
    for size in size_values
        model = ToricCodeTrajectoryModel(size, size)
        for error_rate in rate_values
            accumulators = [_WelfordAccumulator() for _ in 1:7]
            spin_accumulators = [_WelfordAccumulator() for _ in 1:3]
            horizontal_batch_totals = zeros(batch_count)
            second_batch_totals = zeros(batch_count)
            fourth_batch_totals = zeros(batch_count)
            batch_counts = zeros(Int, batch_count)
            for shot in 1:shot_count
                trajectory = sample_trajectory(
                    rng, model; error_rate=error_rate)
                values = _observable_values(trajectory_observables(trajectory))
                for index in eachindex(values)
                    _push!(accumulators[index], values[index])
                end
                spin = marginal_spin_observables(
                    spin_rng, trajectory; spin_samples=spin_sample_count)
                spin_values = (
                    spin.absolute_magnetization, spin.second_moment,
                    spin.fourth_moment)
                for index in eachindex(spin_values)
                    _push!(spin_accumulators[index], spin_values[index])
                end
                batch = fld((shot - 1) * batch_count, shot_count) + 1
                horizontal_batch_totals[batch] += values[6]
                second_batch_totals[batch] += spin.second_moment
                fourth_batch_totals[batch] += spin.fourth_moment
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

            spin_means = [
                accumulator.mean for accumulator in spin_accumulators]
            spin_standard_errors = _standard_error.(spin_accumulators)
            second_batches = second_batch_totals ./ batch_counts
            fourth_batches = fourth_batch_totals ./ batch_counts
            binder = _binder_cumulant(spin_means[2], spin_means[3])
            binder_se = _binder_jackknife_se(
                batch_counts, second_batches, fourth_batches)
            push!(marginal_spin_points, MarginalSpinScanPoint(
                size, error_rate, shot_count,
                spin_means[1], spin_standard_errors[1],
                spin_means[2], spin_standard_errors[2],
                spin_means[3], spin_standard_errors[3],
                binder, binder_se,
                copy(batch_counts), second_batches, fourth_batches))
        end
    end
    return TrajectoryScan(
        size_values, rate_values, points, marginal_spin_points)
end

function _scan_point(scan::TrajectoryScan, size::Int, error_rate::Float64)
    matches = filter(
        point -> point.size == size && point.error_rate == error_rate,
        scan.points)
    length(matches) == 1 || throw(ArgumentError(
        "scan must contain exactly one point for L=$size, p=$error_rate"))
    return only(matches)
end

function _marginal_spin_point(
        scan::TrajectoryScan, size::Int, error_rate::Float64)
    matches = filter(
        point -> point.size == size && point.error_rate == error_rate,
        scan.marginal_spin_points)
    length(matches) == 1 || throw(ArgumentError(
        "scan must contain exactly one marginal spin point for " *
        "L=$size, p=$error_rate"))
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

_isotonic_non_increasing(values::Vector{Float64}) =
    -_isotonic_non_decreasing(-values)

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

function _bootstrap_binder(
        rng::Random.AbstractRNG, point::MarginalSpinScanPoint)
    batch_count = length(point.batch_counts)
    indices = rand(rng, 1:batch_count, batch_count)
    sampled_count = sum(point.batch_counts[index] for index in indices)
    second = sum(
        point.batch_counts[index] * point.second_moment_batches[index]
        for index in indices) / sampled_count
    fourth = sum(
        point.batch_counts[index] * point.fourth_moment_batches[index]
        for index in indices) / sampled_count
    return _binder_cumulant(second, fourth)
end

"""Estimate adjacent-size marginalized Binder crossings by batch bootstrap."""
function estimate_binder_crossings(
        rng::Random.AbstractRNG, scan::TrajectoryScan;
        bootstrap::Integer=2_000, confidence::Real=0.95)
    length(scan.sizes) >= 2 || throw(ArgumentError(
        "at least two sizes are required for Binder crossings"))
    length(scan.error_rates) >= 2 || throw(ArgumentError(
        "at least two error rates are required for Binder crossings"))
    bootstrap > 0 || throw(ArgumentError(
        "bootstrap must be positive, got $bootstrap"))
    isfinite(confidence) && 0 < confidence < 1 || throw(ArgumentError(
        "confidence must be strictly between 0 and 1, got $confidence"))
    isempty(scan.marginal_spin_points) && throw(ArgumentError(
        "scan does not contain marginalized spin points"))
    for point in scan.marginal_spin_points
        batch_count = length(point.batch_counts)
        batch_count >= 2 || throw(ArgumentError(
            "Binder crossing estimation requires at least two batches per point"))
        length(point.second_moment_batches) == batch_count ||
            throw(DimensionMismatch(
                "second-moment batches do not match batch counts"))
        length(point.fourth_moment_batches) == batch_count ||
            throw(DimensionMismatch(
                "fourth-moment batches do not match batch counts"))
        all(>(0), point.batch_counts) || throw(ArgumentError(
            "Binder batch counts must be positive"))
    end

    results = CriticalCrossing[]
    for pair in 1:(length(scan.sizes) - 1)
        small_size, large_size = scan.sizes[pair], scan.sizes[pair + 1]
        small_points = [
            _marginal_spin_point(scan, small_size, rate)
            for rate in scan.error_rates]
        large_points = [
            _marginal_spin_point(scan, large_size, rate)
            for rate in scan.error_rates]
        small_curve = _isotonic_non_increasing(
            [point.binder_cumulant for point in small_points])
        large_curve = _isotonic_non_increasing(
            [point.binder_cumulant for point in large_points])
        selected = _unique_crossing(
            scan.error_rates, small_curve, large_curve)
        if selected.status != :ok
            push!(results, CriticalCrossing(
                small_size, large_size, missing, missing, missing,
                Float64(confidence), 0.0, selected.status))
            continue
        end

        bootstrap_estimates = Float64[]
        for _ in 1:Int(bootstrap)
            small_sample = _isotonic_non_increasing([
                _bootstrap_binder(rng, point) for point in small_points])
            large_sample = _isotonic_non_increasing([
                _bootstrap_binder(rng, point) for point in large_points])
            sample_selection = _unique_crossing(
                scan.error_rates, small_sample, large_sample)
            sample_selection.status == :ok && push!(
                bootstrap_estimates, sample_selection.estimate)
        end
        valid_fraction = length(bootstrap_estimates) / bootstrap
        if valid_fraction < 0.8
            push!(results, CriticalCrossing(
                small_size, large_size, selected.estimate, missing, missing,
                Float64(confidence), valid_fraction, :unstable))
            continue
        end
        tail = (1 - confidence) / 2
        push!(results, CriticalCrossing(
            small_size, large_size, selected.estimate,
            Statistics.quantile(bootstrap_estimates, tail),
            Statistics.quantile(bootstrap_estimates, 1 - tail),
            Float64(confidence), valid_fraction, :ok))
    end
    return results
end
