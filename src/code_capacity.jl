export OpenCodeCapacityScanPoint, OpenCodeCapacityScan,
       scan_open_code_capacity, estimate_open_code_crossings

const _stim_ref = Ref{Py}()

_stim() = isassigned(_stim_ref) ? _stim_ref[] :
    (_stim_ref[] = pyimport("stim"))

"""Open rotated planar surface-code patch for code-capacity experiments."""
struct RotatedCodeCapacityModel
    distance::Int
    base_circuit::Py
    data_qubits::Vector{Int}
end

function _data_qubits_from_terminal_measurement(circuit::Py)
    instructions = split(string(circuit), '\n')
    measurement_indices = findall(instructions) do instruction
        startswith(strip(instruction), "M ")
    end
    isempty(measurement_indices) && throw(ErrorException(
        "generated rotated patch has no data-basis measurement"))

    terminal_index = last(measurement_indices)
    terminal_index == length(instructions) || begin
        trailing = instructions[(terminal_index + 1):end]
        all(instruction -> startswith(strip(instruction), "DETECTOR(") ||
                           startswith(strip(instruction), "OBSERVABLE_INCLUDE(") ||
                           isempty(strip(instruction)), trailing) || throw(ErrorException(
            "terminal data-basis measurement is not followed only by annotations"))
    end
    terminal_measurements = split(strip(instructions[terminal_index]))
    length(terminal_measurements) > 1 || throw(ErrorException(
        "terminal data-basis measurement has no data qubits"))
    data_qubits = parse.(Int, terminal_measurements[2:end])
    length(unique(data_qubits)) == length(data_qubits) || throw(ErrorException(
        "terminal data-basis measurement contains duplicate data qubits"))
    return data_qubits
end

function RotatedCodeCapacityModel(distance::Integer)
    distance >= 3 || throw(ArgumentError("distance must be at least 3"))
    base = _stim().Circuit.generated(
        "surface_code:rotated_memory_z"; distance=Int(distance), rounds=1)
    data = _data_qubits_from_terminal_measurement(base)
    length(data) == Int(distance)^2 || throw(ErrorException(
        "generated rotated patch did not expose d² data qubits"))
    return RotatedCodeCapacityModel(Int(distance), base, data)
end

function _circuit_after_first_tick(base::Py)
    instructions = split(chomp(string(base)), '\n')
    tick = findfirst(instruction -> strip(instruction) == "TICK", instructions)
    tick === nothing && throw(ErrorException("generated circuit has no preparation TICK"))
    prefix = join(instructions[1:tick], "\n")
    suffix = join(instructions[(tick + 1):end], "\n")
    return prefix, suffix
end

function _data_x_circuit(
        model::RotatedCodeCapacityModel, operation::String, qubits::Vector{Int})
    prefix, suffix = _circuit_after_first_tick(model.base_circuit)
    injected = isempty(qubits) ? "" : "$(operation) $(join(qubits, ' '))"
    text = join((prefix, injected, suffix), "\n")
    return _stim().Circuit(text)
end

function _circuit_with_data_x_noise(model::RotatedCodeCapacityModel, error_rate::Real)
    isfinite(error_rate) && 0 <= error_rate < 0.5 || throw(ArgumentError(
        "error_rate must be finite and lie in [0, 0.5), got $error_rate"))
    error_rate == 0 && return model.base_circuit
    return _data_x_circuit(
        model, "X_ERROR($(Float64(error_rate)))", model.data_qubits)
end

"""Return the noiseless extraction circuit with independent data-qubit X noise."""
rotated_code_capacity_circuit(model::RotatedCodeCapacityModel, error_rate::Real) =
    _circuit_with_data_x_noise(model, error_rate)

function _deterministic_data_x_circuit(
        model::RotatedCodeCapacityModel, support::AbstractVector{<:Integer})
    support_values = Int[support...]
    all(qubit -> qubit in model.data_qubits, support_values) || throw(ArgumentError(
        "deterministic X support must contain only data qubits"))
    length(unique(support_values)) == length(support_values) || throw(ArgumentError(
        "deterministic X support must not repeat data qubits"))
    isempty(support_values) && return model.base_circuit
    return _data_x_circuit(model, "X_ERROR(1)", support_values)
end

function _sample_deterministic_x_support(
        model::RotatedCodeCapacityModel, support::Vector{Int})
    circuit = _deterministic_data_x_circuit(model, support)
    detectors, observables = circuit.compile_detector_sampler().sample(
        shots=1, separate_observables=true)
    return pyconvert(BitMatrix, detectors), pyconvert(BitMatrix, observables)
end

function _minimum_x_logical_support(model::RotatedCodeCapacityModel)
    coordinates = Dict{Int,Tuple{Int,Int}}()
    for line in split(string(model.base_circuit), '\n')
        match_result = match(r"^QUBIT_COORDS\(([-0-9]+), ([-0-9]+)\) ([0-9]+)$", strip(line))
        match_result === nothing && continue
        x, y, qubit = parse.(Int, match_result.captures)
        qubit in model.data_qubits && (coordinates[qubit] = (x, y))
    end
    length(coordinates) == length(model.data_qubits) || throw(ErrorException(
        "generated rotated patch did not expose coordinates for every data qubit"))

    rows = Dict{Int,Vector{Int}}()
    columns = Dict{Int,Vector{Int}}()
    for (qubit, (x, y)) in coordinates
        push!(get!(rows, y, Int[]), qubit)
        push!(get!(columns, x, Int[]), qubit)
    end
    candidate_families = Vector{Vector{Vector{Int}}}()
    for family in (values(rows), values(columns))
        candidates = Vector{Vector{Int}}()
        for members in family
            support = sort(members)
            length(support) == model.distance || continue
            detectors, observables = _sample_deterministic_x_support(model, support)
            !any(detectors) && only(vec(observables)) && push!(candidates, support)
        end
        isempty(candidates) || push!(candidate_families, candidates)
    end
    length(candidate_families) == 1 || throw(ErrorException(
        "expected exactly one orientation family of minimum X logical strings"))
    support = first(sort(only(candidate_families)))
    length(support) == model.distance || throw(ErrorException(
        "minimum X logical support has the wrong length"))
    return support
end

"""One decoded logical-X failure estimate for a rotated planar patch."""
struct RotatedCodeCapacityScanPoint
    distance::Int
    error_rate::Float64
    shots::Int
    logical_x_failure_count::Int
    logical_x_failure_rate::Float64
    logical_x_failure_se::Float64
    logical_x_failure_batches::Vector{Float64}
    seed::Union{Missing,Int}
end

"""Finite-size rotated-planar code-capacity scan on a shared error-rate grid."""
struct RotatedCodeCapacityScan
    distances::Vector{Int}
    error_rates::Vector{Float64}
    points::Vector{RotatedCodeCapacityScanPoint}
end

function _rotated_code_matching(model::RotatedCodeCapacityModel, p::Float64)
    p > 0 || throw(ArgumentError("matching requires p > 0"))
    dem = rotated_code_capacity_circuit(model, p).detector_error_model(
        decompose_errors=true,
        block_decomposition_from_introducing_remnant_edges=true)
    return _pymatching().Matching.from_detector_error_model(dem)
end

_decode_observables(matching::Py, syndrome) = matching.decode_batch(syndrome)

function _rotated_detector_count(model::RotatedCodeCapacityModel, circuit::Py)
    expected = model.distance^2 - 1
    pyconvert(Int, model.base_circuit.num_detectors) == expected || throw(ErrorException(
        "rotated model has $(pyconvert(Int, model.base_circuit.num_detectors)) " *
        "detectors; expected $expected for distance $(model.distance)"))
    pyconvert(Int, circuit.num_detectors) == expected || throw(ErrorException(
        "rotated code-capacity circuit has $(pyconvert(Int, circuit.num_detectors)) " *
        "detectors; expected $expected for distance $(model.distance)"))
    return expected
end

function _rotated_code_capacity_error_rate(error_rate::Real)
    isfinite(error_rate) && 0 <= error_rate < 0.5 || throw(ArgumentError(
        "error_rate must be finite and lie in [0, 0.5), got $error_rate"))
    return Float64(error_rate)
end

function _rotated_result_seed(seed)
    seed === nothing && return missing
    try
        return Int(seed)
    catch error
        throw(ArgumentError("seed must be convertible to Int, got $seed"))
    end
end

function _sample_rotated_detectors(
        rng::Random.AbstractRNG, circuit::Py, shots::Int)
    sampler = circuit.compile_detector_sampler(seed=rand(rng, UInt64))
    syndrome_py, actual_py = sampler.sample(
        shots=shots, separate_observables=true)
    return syndrome_py, actual_py
end

function _validate_rotated_samples(
        syndrome::BitMatrix, actual::BitMatrix, shots::Int, detector_count::Int)
    size(syndrome) == (shots, detector_count) || throw(ErrorException(
        "Stim returned detector samples with size $(size(syndrome)); expected " *
        "($shots, $detector_count)"))
    size(actual) == (shots, 1) || throw(ErrorException(
        "Stim returned logical-observable samples with size $(size(actual)); " *
        "expected ($shots, 1)"))
    return nothing
end

"""
    estimate_rotated_code_capacity(rng, model, error_rate;
                                   shots=10_000, batches=100, seed=nothing)

Sample data-X code-capacity noise and score detector-only MWPM predictions
against Stim's held-out logical-X observable.
"""
function estimate_rotated_code_capacity(
        rng::Random.AbstractRNG, model::RotatedCodeCapacityModel, error_rate::Real;
        shots::Integer=10_000, batches::Integer=100, seed=nothing)
    rate = _rotated_code_capacity_error_rate(error_rate)
    shots > 0 || throw(ArgumentError("shots must be positive, got $shots"))
    2 <= batches <= shots || throw(ArgumentError(
        "batches must satisfy 2 <= batches <= shots, got batches=$batches and shots=$shots"))
    shots % batches == 0 || throw(ArgumentError(
        "shots must be divisible by batches, got shots=$shots and batches=$batches"))

    shot_count, batch_count = Int(shots), Int(batches)
    result_seed = _rotated_result_seed(seed)
    circuit = rotated_code_capacity_circuit(model, rate)
    detector_count = _rotated_detector_count(model, circuit)
    syndrome_py, actual_py = _sample_rotated_detectors(
        rng, circuit, shot_count)
    if rate == 0
        syndrome = pyconvert(BitMatrix, syndrome_py)
        actual = pyconvert(BitMatrix, actual_py)
        _validate_rotated_samples(syndrome, actual, shot_count, detector_count)
        !any(syndrome) && !any(actual) || throw(ErrorException(
            "noiseless rotated code-capacity sampling produced detector or logical flips"))
        return RotatedCodeCapacityScanPoint(
            model.distance, rate, shot_count, 0, 0.0, 0.0,
            zeros(batch_count), result_seed)
    end

    matching = _rotated_code_matching(model, rate)
    predicted_py = _decode_observables(matching, syndrome_py)
    syndrome = pyconvert(BitMatrix, syndrome_py)
    actual = pyconvert(BitMatrix, actual_py)
    _validate_rotated_samples(syndrome, actual, shot_count, detector_count)
    predicted = pyconvert(BitMatrix, predicted_py)
    size(predicted) == size(actual) || throw(ErrorException(
        "PyMatching returned logical-observable predictions with size $(size(predicted)); " *
        "expected $(size(actual))"))
    failures = vec(any(predicted .!= actual; dims=2))
    failure_count = count(failures)
    failure_rate = failure_count / shot_count
    batch_size = div(shot_count, batch_count)
    batch_rates = [
        count(failures[((batch - 1) * batch_size + 1):(batch * batch_size)]) / batch_size
        for batch in 1:batch_count
    ]
    failure_se = failure_rate == 0 || failure_rate == 1 ? 0.0 :
        sqrt(failure_rate * (1 - failure_rate) / shot_count)
    return RotatedCodeCapacityScanPoint(
        model.distance, rate, shot_count, failure_count, failure_rate, failure_se,
        batch_rates, result_seed)
end

"""
    scan_rotated_code_capacity(rng, distances, error_rates;
                               shots=10_000, batches=100, seed=nothing,
                               progress_io=nothing)

Estimate logical-X decoder failure across strictly increasing rotated-patch
distances and data-X error rates. Each point retains equal-size batch rates
for subsequent bootstrap crossing estimates.
"""
function scan_rotated_code_capacity(
        rng::Random.AbstractRNG, distances, error_rates;
        shots::Integer=10_000, batches::Integer=100, seed=nothing,
        progress_io=nothing)
    distance_values = Int[value for value in distances]
    rate_values = [_rotated_code_capacity_error_rate(value) for value in error_rates]
    _validate_strictly_increasing(distance_values, "distances")
    all(distance -> distance >= 3, distance_values) || throw(ArgumentError(
        "distances must be at least 3"))
    _validate_strictly_increasing(rate_values, "error_rates")
    shots > 0 || throw(ArgumentError("shots must be positive, got $shots"))
    2 <= batches <= shots || throw(ArgumentError(
        "batches must satisfy 2 <= batches <= shots, got batches=$batches and shots=$shots"))
    shots % batches == 0 || throw(ArgumentError(
        "shots must be divisible by batches, got shots=$shots and batches=$batches"))
    scan_seed = seed === nothing ? nothing : _rotated_result_seed(seed)

    shot_count = Int(shots)
    models = Dict(distance => RotatedCodeCapacityModel(distance)
                  for distance in distance_values)
    points = RotatedCodeCapacityScanPoint[]
    for distance in distance_values
        model = models[distance]
        for error_rate in rate_values
            point = estimate_rotated_code_capacity(
                rng, model, error_rate;
                shots=shots, batches=batches, seed=scan_seed)
            push!(points, point)
            progress_io === nothing || println(
                progress_io,
                "completed d=$distance p=$error_rate ($shot_count shots)")
        end
    end
    return RotatedCodeCapacityScan(distance_values, rate_values, points)
end

function _rotated_code_capacity_scan_point(
        scan::RotatedCodeCapacityScan, distance::Int, error_rate::Float64)
    matches = filter(
        point -> point.distance == distance && point.error_rate == error_rate,
        scan.points)
    length(matches) == 1 || throw(ArgumentError(
        "scan must contain exactly one point for d=$distance, p=$error_rate"))
    return only(matches)
end

"""
    estimate_rotated_code_crossings(rng, scan; bootstrap=2_000, confidence=0.95)

Estimate crossings of adjacent rotated-patch logical-X failure curves using
monotone fits and a nonparametric bootstrap over equal-length stored batches.
"""
function estimate_rotated_code_crossings(
        rng::Random.AbstractRNG, scan::RotatedCodeCapacityScan;
        bootstrap::Integer=2_000, confidence::Real=0.95)
    length(scan.distances) >= 2 || throw(ArgumentError(
        "at least two distances are required for crossings"))
    length(scan.error_rates) >= 2 || throw(ArgumentError(
        "at least two error rates are required for crossings"))
    _validate_strictly_increasing(scan.distances, "distances")
    all(distance -> distance >= 3, scan.distances) || throw(ArgumentError(
        "distances must be at least 3"))
    foreach(_rotated_code_capacity_error_rate, scan.error_rates)
    _validate_strictly_increasing(scan.error_rates, "error_rates")
    bootstrap > 0 || throw(ArgumentError(
        "bootstrap must be positive, got $bootstrap"))
    isfinite(confidence) && 0 < confidence < 1 || throw(ArgumentError(
        "confidence must be strictly between 0 and 1, got $confidence"))

    expected_points = length(scan.distances) * length(scan.error_rates)
    length(scan.points) == expected_points || throw(ArgumentError(
        "scan must contain one point for every distance and error rate"))
    batch_lengths = length.(getfield.(scan.points, :logical_x_failure_batches))
    all(>=(2), batch_lengths) || throw(ArgumentError(
        "crossing estimation requires at least two batches per scan point"))
    all(==(first(batch_lengths)), batch_lengths) || throw(ArgumentError(
        "crossing estimation requires equal batch lengths at every scan point"))
    all(point -> point.shots > 0 &&
                 point.shots % length(point.logical_x_failure_batches) == 0,
        scan.points) || throw(ArgumentError(
        "crossing estimation requires batches with integral shot counts"))
    batch_sizes = [div(point.shots, length(point.logical_x_failure_batches))
                   for point in scan.points]
    all(==(first(batch_sizes)), batch_sizes) || throw(ArgumentError(
        "crossing estimation requires equal-size stored batches"))

    results = CriticalCrossing[]
    for pair in 1:(length(scan.distances) - 1)
        small_distance, large_distance =
            scan.distances[pair], scan.distances[pair + 1]
        small_points = [_rotated_code_capacity_scan_point(
            scan, small_distance, rate) for rate in scan.error_rates]
        large_points = [_rotated_code_capacity_scan_point(
            scan, large_distance, rate) for rate in scan.error_rates]
        small_curve = _isotonic_non_decreasing(
            [point.logical_x_failure_rate for point in small_points])
        large_curve = _isotonic_non_decreasing(
            [point.logical_x_failure_rate for point in large_points])
        selected = _selected_crossing(
            scan.error_rates, small_curve, large_curve)
        if selected.status != :ok
            push!(results, CriticalCrossing(
                small_distance, large_distance, missing, missing, missing,
                Float64(confidence), 0.0, selected.status))
            continue
        end
        estimate = selected.estimate

        bootstrap_estimates = Float64[]
        for _ in 1:Int(bootstrap)
            small_sample = _isotonic_non_decreasing([
                _bootstrap_batch_mean(rng, point.logical_x_failure_batches)
                for point in small_points
            ])
            large_sample = _isotonic_non_decreasing([
                _bootstrap_batch_mean(rng, point.logical_x_failure_batches)
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
                small_distance, large_distance, estimate, missing, missing,
                Float64(confidence), valid_fraction, :unstable))
            continue
        end
        tail = (1 - confidence) / 2
        push!(results, CriticalCrossing(
            small_distance, large_distance, estimate,
            Statistics.quantile(bootstrap_estimates, tail),
            Statistics.quantile(bootstrap_estimates, 1 - tail),
            Float64(confidence), valid_fraction, :ok))
    end
    return results
end

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
