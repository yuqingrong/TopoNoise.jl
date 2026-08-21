"""Open isometric planar patch parameterized directly by its code distance."""
struct IsometricPlanarCodeModel
    distance::Int
    rows::Int
    cols::Int
    layers::Vector{Vector{Tuple{Int,Int}}}
    wire_count::Int

    function IsometricPlanarCodeModel(distance::Integer)
        distance >= 2 || throw(ArgumentError(
            "distance must be at least 2, got $distance"))
        side = Int(distance)
        return new(side, side, side, _site_schedule(side, side), 2 * side + 4)
    end
end

"""Independent X events on the internal virtual bonds of an open planar patch."""
struct IsometricPlanarVirtualErrors
    horizontal::BitMatrix
    vertical::BitMatrix

    function IsometricPlanarVirtualErrors(horizontal::BitMatrix, vertical::BitMatrix)
        rows, horizontal_cols = size(horizontal)
        vertical_rows, cols = size(vertical)
        rows >= 2 || throw(ArgumentError(
            "virtual-link arrays must have at least 2 rows"))
        cols >= 2 || throw(ArgumentError(
            "virtual-link arrays must have at least 2 columns"))
        horizontal_cols == cols - 1 || throw(DimensionMismatch(
            "horizontal virtual links must have size (R, C - 1)"))
        vertical_rows == rows - 1 || throw(DimensionMismatch(
            "vertical virtual links must have size (R - 1, C)"))
        return new(horizontal, vertical)
    end
end

"""Detector record: plaquettes plus local north and south boundary checks."""
struct IsometricPlanarSyndrome
    plaquettes::BitMatrix
    north::BitVector
    south::BitVector
end

"""Non-destructive Yao encoder schedule with every terminal carrier retained."""
struct IsometricPlanarEncoder
    model::IsometricPlanarCodeModel
    wire_count::Int
    output_wires::Array{Int,3}
    operations::Vector{Any}
end

"""Perfect Z-parity check supports for the planar encoder's physical outputs."""
struct IsometricPlanarCheckLayer
    detector_supports::Vector{Vector{Int}}
    logical_z_support::Vector{Int}
end

Base.length(errors::IsometricPlanarVirtualErrors) =
    length(errors.horizontal) + length(errors.vertical)
Base.count(errors::IsometricPlanarVirtualErrors) =
    count(errors.horizontal) + count(errors.vertical)

function _isometric_planar_error_rate(error_rate::Real)
    isfinite(error_rate) && 0 <= error_rate < 0.5 || throw(ArgumentError(
        "error_rate must be finite and lie in [0, 0.5), got $error_rate"))
    return Float64(error_rate)
end

function _validate_isometric_planar_errors(
        model::IsometricPlanarCodeModel, errors::IsometricPlanarVirtualErrors)
    size(errors.horizontal) == (model.rows, model.cols - 1) || throw(DimensionMismatch(
        "horizontal virtual links must have size ($(model.rows), $(model.cols - 1))"))
    size(errors.vertical) == (model.rows - 1, model.cols) || throw(DimensionMismatch(
        "vertical virtual links must have size ($(model.rows - 1), $(model.cols))"))
    return nothing
end

"""Sample independent internal virtual-bond X errors."""
function sample_isometric_planar_virtual_errors(
        rng::Random.AbstractRNG, model::IsometricPlanarCodeModel;
        error_rate::Real)
    probability = _isometric_planar_error_rate(error_rate)
    return IsometricPlanarVirtualErrors(
        BitMatrix(rand(rng, model.rows, model.cols - 1) .< probability),
        BitMatrix(rand(rng, model.rows - 1, model.cols) .< probability))
end

"""Return plaquette and local north/south boundary detectors for a bond record."""
function isometric_planar_syndrome(record)
    horizontal, vertical = record.horizontal, record.vertical
    plaquettes = plaquette_syndrome((; horizontal, vertical))
    return IsometricPlanarSyndrome(
        plaquettes, BitVector(horizontal[1, :]), BitVector(horizontal[end, :]))
end

function _isometric_planar_contracted_tensor(kind::Symbol, ::Type{T}) where {T<:Number}
    tensor = toric_code_local_tensor(T)
    if kind === :bulk
        return tensor, one(T)
    elseif kind === :east
        return (selectdim(tensor, 5, 1) + selectdim(tensor, 5, 2)) /
               sqrt(convert(T, 2)), sqrt(convert(T, 2))
    elseif kind === :north
        return selectdim(tensor, 6, 1), sqrt(convert(T, 2))
    elseif kind === :corner
        east_plus = (selectdim(tensor, 5, 1) + selectdim(tensor, 5, 2)) /
                    sqrt(convert(T, 2))
        return selectdim(east_plus, 5, 1), convert(T, 2)
    end
    throw(ArgumentError("local gate kind must be :bulk, :east, :north, or :corner"))
end

function _isometric_planar_complete_gate(
        isometry::AbstractMatrix{T}, qubits::Int) where {T}
    completed = unitary_completion(isometry)
    raw = reshape(completed, ntuple(_ -> 2, 2 * qubits))
    output_indices = collect(1:qubits)
    physical_input_indices = collect((qubits + 3):(2 * qubits))
    return permutedims(raw, (
        output_indices..., physical_input_indices..., qubits + 1, qubits + 2))
end

const _isometric_planar_gate_cache = Dict{Tuple{Symbol,DataType},Any}()

function _isometric_planar_local_gate(
        kind::Symbol, ::Type{T}=ComplexF64) where {T<:Number}
    key = (kind, T)
    return get!(_isometric_planar_gate_cache, key) do
        kind === :bulk && return toric_code_local_gate(T)
        contracted, normalization = _isometric_planar_contracted_tensor(kind, T)
        qubits = kind === :corner ? 4 : 5
        isometry = normalization * reshape(contracted, 1 << qubits, 4)
        gate = _isometric_planar_complete_gate(isometry, qubits)
        return kind === :north ?
            permutedims(gate, (1, 2, 3, 5, 4, 6, 7, 8, 9, 10)) : gate
    end
end

"""
    isometric_planar_encoder(model) -> IsometricPlanarEncoder

Construct the non-destructive local-gate schedule for the fixed west/east
`|+⟩` and south/north `|0⟩` boundary convention. At the east and north
boundaries the terminal carriers occupy the fourth physical-output slot, so
the retained output register has exactly `4d²` wires.
"""
function isometric_planar_encoder(model::IsometricPlanarCodeModel)
    distance = model.distance
    fresh_count = 4 * (distance - 1)^2 + 6 * (distance - 1) + 2
    wire_count = 4 * distance^2
    first_horizontal = fresh_count + 1
    first_vertical = first_horizontal + distance
    horizontal_wire(row) = first_horizontal + row - 1
    vertical_wire(col) = first_vertical + col - 1
    output_wires = Array{Int}(undef, distance, distance, 4)
    operations = Any[Yao.repeat(wire_count, H, first_horizontal:(first_vertical - 1))]
    next_fresh = 1
    matrices = Dict(kind => reshape(_isometric_planar_local_gate(kind),
                                    1 << (kind === :bulk ? 6 : kind === :corner ? 4 : 5),
                                    1 << (kind === :bulk ? 6 : kind === :corner ? 4 : 5))
                    for kind in (:bulk, :east, :north, :corner))

    for layer in model.layers, (row, col) in layer
        horizontal, vertical = horizontal_wire(row), vertical_wire(col)
        east_terminal, north_terminal = col == distance, row == 1
        if east_terminal && north_terminal
            fresh = next_fresh:(next_fresh + 1)
            next_fresh += 2
            wires = (fresh..., horizontal, vertical)
            output_wires[row, col, :] .= (fresh..., horizontal, vertical)
            kind = :corner
        elseif east_terminal
            fresh = next_fresh:(next_fresh + 2)
            next_fresh += 3
            wires = (fresh..., horizontal, vertical)
            output_wires[row, col, :] .= (fresh..., horizontal)
            kind = :east
        elseif north_terminal
            fresh = next_fresh:(next_fresh + 2)
            next_fresh += 3
            wires = (fresh..., horizontal, vertical)
            output_wires[row, col, :] .= (fresh..., vertical)
            kind = :north
        else
            fresh = next_fresh:(next_fresh + 3)
            next_fresh += 4
            wires = (fresh..., horizontal, vertical)
            output_wires[row, col, :] .= fresh
            kind = :bulk
        end
        push!(operations, subroutine(
            wire_count, _yao_local_gate(matrices[kind], (row, col)), wires))
    end
    next_fresh == fresh_count + 1 || throw(ErrorException(
        "planar encoder allocated an inconsistent number of fresh physical wires"))
    return IsometricPlanarEncoder(model, wire_count, output_wires, operations)
end

function _isometric_planar_horizontal_support(
        encoder::IsometricPlanarEncoder, row::Int, col::Int)
    return [encoder.output_wires[row, col, 1], encoder.output_wires[row, col + 1, 3]]
end

function _isometric_planar_vertical_support(
        encoder::IsometricPlanarEncoder, row::Int, col::Int)
    return [encoder.output_wires[row + 1, col, 2], encoder.output_wires[row, col, 4]]
end

"""Build the plaquette, local-boundary, and held-out logical-Z parity supports."""
function isometric_planar_check_layer(encoder::IsometricPlanarEncoder)
    distance = encoder.model.distance
    detectors = Vector{Vector{Int}}()
    for row in 1:(distance - 1), col in 1:(distance - 1)
        support = Int[]
        append!(support, _isometric_planar_horizontal_support(encoder, row, col))
        append!(support, _isometric_planar_horizontal_support(encoder, row + 1, col))
        append!(support, _isometric_planar_vertical_support(encoder, row, col))
        append!(support, _isometric_planar_vertical_support(encoder, row, col + 1))
        push!(detectors, support)
    end
    for col in 1:(distance - 1)
        push!(detectors, _isometric_planar_horizontal_support(encoder, 1, col))
    end
    for col in 1:(distance - 1)
        push!(detectors, _isometric_planar_horizontal_support(encoder, distance, col))
    end
    logical_z = Int[]
    for row in 1:(distance - 1)
        append!(logical_z, _isometric_planar_vertical_support(encoder, row, 1))
    end
    return IsometricPlanarCheckLayer(detectors, logical_z)
end

"""Measured/reset Yao realization retained as a small-distance detector oracle."""
struct IsometricPlanarYaoTrajectory
    measurements::Array{Bool,3}
    errors::IsometricPlanarVirtualErrors
end

function _store_isometric_planar_outcome!(
        measurements::Array{Bool,3}, row::Int, col::Int, outcome)
    for direction in 1:4
        measurements[row, col, direction] = Bool(outcome[direction])
    end
    return measurements
end

"""
    sample_isometric_planar_yao_trajectory(rng, model, errors; max_qubits=24)

Execute the measured/reset circuit only as a small-system reference. Its
reused physical ancillas keep the state-vector register at `2d+4` qubits;
production capacity scans intentionally use the exact detector sampler.
"""
function sample_isometric_planar_yao_trajectory(
        rng::Random.AbstractRNG, model::IsometricPlanarCodeModel,
        errors::IsometricPlanarVirtualErrors; max_qubits::Integer=24)
    _validate_isometric_planar_errors(model, errors)
    max_qubits > 0 || throw(ArgumentError(
        "max_qubits must be positive, got $max_qubits"))
    model.wire_count <= max_qubits || throw(ArgumentError(
        "isometric Yao reference requires $(model.wire_count) qubits, " *
        "exceeding max_qubits=$max_qubits"))

    physical = 1:4
    first_horizontal = 5
    first_vertical = first_horizontal + model.rows
    horizontal_wire(row) = first_horizontal + row - 1
    vertical_wire(col) = first_vertical + col - 1
    register = Yao.zero_state(model.wire_count)
    Yao.apply!(register, Yao.repeat(
        model.wire_count, H, first_horizontal:(first_vertical - 1)))
    matrices = Dict(kind => reshape(_isometric_planar_local_gate(kind),
                                    1 << (kind === :bulk ? 6 : kind === :corner ? 4 : 5),
                                    1 << (kind === :bulk ? 6 : kind === :corner ? 4 : 5))
                    for kind in (:bulk, :east, :north, :corner))
    measurements = Array{Bool}(undef, model.rows, model.cols, 4)

    for layer in model.layers, (row, col) in layer
        horizontal, vertical = horizontal_wire(row), vertical_wire(col)
        east_terminal, north_terminal = col == model.cols, row == 1
        if east_terminal && north_terminal
            wires = (1, 2, horizontal, vertical)
            Yao.apply!(register, subroutine(
                model.wire_count, _yao_local_gate(matrices[:corner], (row, col)), wires))
        elseif east_terminal
            wires = (1, 2, 3, horizontal, vertical)
            Yao.apply!(register, subroutine(
                model.wire_count, _yao_local_gate(matrices[:east], (row, col)), wires))
            wires = (1, 2, 3, horizontal)
        elseif north_terminal
            wires = (1, 2, 3, horizontal, vertical)
            Yao.apply!(register, subroutine(
                model.wire_count, _yao_local_gate(matrices[:north], (row, col)), wires))
            wires = (1, 2, 3, vertical)
        else
            wires = (physical..., horizontal, vertical)
            Yao.apply!(register, subroutine(
                model.wire_count, _yao_local_gate(matrices[:bulk], (row, col)), wires))
            wires = physical
        end
        outcome = Yao.measure!(Yao.ResetTo(0), register, wires; rng=rng)
        _store_isometric_planar_outcome!(measurements, row, col, outcome)

        col < model.cols && errors.horizontal[row, col] &&
            Yao.apply!(register, Yao.repeat(model.wire_count, Yao.X, (horizontal,)))
        row > 1 && errors.vertical[row - 1, col] &&
            Yao.apply!(register, Yao.repeat(model.wire_count, Yao.X, (vertical,)))
    end
    return IsometricPlanarYaoTrajectory(measurements, errors)
end

function sample_isometric_planar_yao_trajectory(
        rng::Random.AbstractRNG, model::IsometricPlanarCodeModel;
        error_rate::Real, max_qubits::Integer=24)
    errors = sample_isometric_planar_virtual_errors(
        rng, model; error_rate=error_rate)
    return sample_isometric_planar_yao_trajectory(
        rng, model, errors; max_qubits=max_qubits)
end

"""Return doubled-edge virtual-bond mismatches from a Yao reference record."""
function isometric_planar_bond_mismatches(trajectory::IsometricPlanarYaoTrajectory)
    rows, cols, directions = size(trajectory.measurements)
    directions == 4 || throw(DimensionMismatch(
        "measurements must store four E/N/W/S bits per tensor site"))
    horizontal = BitMatrix(undef, rows, cols - 1)
    vertical = BitMatrix(undef, rows - 1, cols)
    for row in 1:rows, col in 1:(cols - 1)
        horizontal[row, col] = xor(trajectory.measurements[row, col, 1],
                                   trajectory.measurements[row, col + 1, 3])
    end
    for row in 1:(rows - 1), col in 1:cols
        vertical[row, col] = xor(trajectory.measurements[row + 1, col, 2],
                                 trajectory.measurements[row, col, 4])
    end
    return (; horizontal, vertical)
end

function _isometric_planar_logical_frame(record)
    parity = false
    for row in axes(record.vertical, 1)
        parity = xor(parity, record.vertical[row, 1])
    end
    return parity
end

_isometric_planar_yao_logical_frame(trajectory::IsometricPlanarYaoTrajectory) =
    _isometric_planar_logical_frame(isometric_planar_bond_mismatches(trajectory))

const _isometric_planar_matching_cache = Dict{Int,Py}()

_isometric_planar_plaquette_detector(distance::Int, row::Int, col::Int) =
    (row - 1) * (distance - 1) + (col - 1)
_isometric_planar_north_detector(distance::Int, col::Int) =
    (distance - 1)^2 + (col - 1)
_isometric_planar_south_detector(distance::Int, col::Int) =
    (distance - 1)^2 + (distance - 1) + (col - 1)
_isometric_planar_detector_count(distance::Int) = (distance - 1)^2 + 2 * (distance - 1)
_isometric_planar_h_fault(distance::Int, row::Int, col::Int) =
    (row - 1) * (distance - 1) + (col - 1)
_isometric_planar_v_fault(distance::Int, row::Int, col::Int) =
    distance * (distance - 1) + (row - 1) * distance + (col - 1)
_isometric_planar_fault_count(distance::Int) = 2 * distance * (distance - 1)

function _isometric_planar_matching(model::IsometricPlanarCodeModel)
    distance = model.distance
    distance == 2 && throw(ArgumentError(
        "distance-2 planar decoding uses its deterministic diagnostic lookup"))
    return get!(_isometric_planar_matching_cache, distance) do
        matching = _pymatching().Matching()
        plaquette(row, col) =
            _isometric_planar_plaquette_detector(distance, row, col)
        north(col) = _isometric_planar_north_detector(distance, col)
        south(col) = _isometric_planar_south_detector(distance, col)
        logical_id = _isometric_planar_fault_count(distance)

        for row in 1:distance, col in 1:(distance - 1)
            fault = _isometric_planar_h_fault(distance, row, col)
            first, second = if row == 1
                north(col), plaquette(1, col)
            elseif row == distance
                plaquette(distance - 1, col), south(col)
            else
                plaquette(row - 1, col), plaquette(row, col)
            end
            matching.add_edge(pyint(first), pyint(second);
                fault_ids=pyset([fault]), weight=1.0,
                merge_strategy="independent")
        end

        for row in 1:(distance - 1), col in 1:distance
            fault = _isometric_planar_v_fault(distance, row, col)
            fault_ids = col == 1 ? [fault, logical_id] : [fault]
            if col == 1 || col == distance
                adjacent = col == 1 ? plaquette(row, 1) :
                                         plaquette(row, distance - 1)
                matching.add_boundary_edge(pyint(adjacent);
                    fault_ids=pyset(fault_ids), weight=1.0,
                    merge_strategy="independent")
            else
                matching.add_edge(pyint(plaquette(row, col - 1)),
                    pyint(plaquette(row, col)); fault_ids=pyset(fault_ids),
                    weight=1.0, merge_strategy="independent")
            end
        end
        matching
    end
end

function _validate_isometric_planar_syndrome(
        model::IsometricPlanarCodeModel, syndrome::IsometricPlanarSyndrome)
    expected_plaquettes = (model.distance - 1, model.distance - 1)
    size(syndrome.plaquettes) == expected_plaquettes || throw(DimensionMismatch(
        "plaquette syndrome must have size $expected_plaquettes"))
    expected_boundary = model.distance - 1
    length(syndrome.north) == expected_boundary || throw(DimensionMismatch(
        "north syndrome must have length $expected_boundary"))
    length(syndrome.south) == expected_boundary || throw(DimensionMismatch(
        "south syndrome must have length $expected_boundary"))
    return nothing
end

function _isometric_planar_syndrome_vector(
        model::IsometricPlanarCodeModel, syndrome::IsometricPlanarSyndrome)
    _validate_isometric_planar_syndrome(model, syndrome)
    distance = model.distance
    values = zeros(Int8, _isometric_planar_detector_count(distance))
    for row in 1:(distance - 1), col in 1:(distance - 1)
        detector = _isometric_planar_plaquette_detector(distance, row, col)
        values[detector + 1] = syndrome.plaquettes[row, col]
    end
    for col in 1:(distance - 1)
        values[_isometric_planar_north_detector(distance, col) + 1] = syndrome.north[col]
        values[_isometric_planar_south_detector(distance, col) + 1] = syndrome.south[col]
    end
    return values
end

"""
    decode_isometric_planar_syndrome(model, syndrome) -> Correction

Decode plaquette and local north/south detector bits. Distance two uses a
fixed west-boundary representative for its unavoidable one-fault tie and is
provided only as a diagnostic geometry.
"""
function decode_isometric_planar_syndrome(
        model::IsometricPlanarCodeModel, syndrome::IsometricPlanarSyndrome)
    syndrome_vector = _isometric_planar_syndrome_vector(model, syndrome)
    correction = Correction(model.rows, model.cols)
    any(value -> !iszero(value), syndrome_vector) || return correction
    if model.distance == 2
        correction.vertical[1, 1] = true
        return correction
    end

    prediction = _isometric_planar_matching(model).decode(
        pyimport("numpy").asarray(syndrome_vector))
    values = pyconvert(Vector{Int}, prediction)
    fault_count = _isometric_planar_fault_count(model.distance)
    length(values) >= fault_count || throw(ErrorException(
        "PyMatching returned $(length(values)) fault predictions; expected at least $fault_count"))
    for row in 1:model.rows, col in 1:(model.cols - 1)
        correction.horizontal[row, col] =
            values[_isometric_planar_h_fault(model.distance, row, col) + 1] == 1
    end
    for row in 1:(model.rows - 1), col in 1:model.cols
        correction.vertical[row, col] =
            values[_isometric_planar_v_fault(model.distance, row, col) + 1] == 1
    end
    return correction
end

"""Return whether the error/correction residual has odd west-boundary parity."""
function isometric_planar_logical_failure(
        errors::IsometricPlanarVirtualErrors, correction::Correction)
    size(errors.horizontal) == size(correction.horizontal) || throw(DimensionMismatch(
        "correction.horizontal must match horizontal virtual links"))
    size(errors.vertical) == size(correction.vertical) || throw(DimensionMismatch(
        "correction.vertical must match vertical virtual links"))
    parity = false
    for row in axes(errors.vertical, 1)
        parity = xor(parity, errors.vertical[row, 1], correction.vertical[row, 1])
    end
    return parity
end

function _isometric_planar_minimum_odd_logical_weight(model::IsometricPlanarCodeModel)
    distance = model.distance
    detector_count = _isometric_planar_detector_count(distance)
    boundary = detector_count + 1
    adjacency = [Tuple{Int,Bool}[] for _ in 1:boundary]
    add_edge!(first::Int, second::Int, logical::Bool=false) = begin
        push!(adjacency[first], (second, logical))
        push!(adjacency[second], (first, logical))
    end
    plaquette(row, col) =
        _isometric_planar_plaquette_detector(distance, row, col) + 1
    north(col) = _isometric_planar_north_detector(distance, col) + 1
    south(col) = _isometric_planar_south_detector(distance, col) + 1

    for row in 1:distance, col in 1:(distance - 1)
        row == 1 ? add_edge!(north(col), plaquette(1, col)) :
        row == distance ? add_edge!(plaquette(distance - 1, col), south(col)) :
        add_edge!(plaquette(row - 1, col), plaquette(row, col))
    end
    for row in 1:(distance - 1), col in 1:distance
        col == 1 ? add_edge!(plaquette(row, 1), boundary, true) :
        col == distance ? add_edge!(plaquette(row, distance - 1), boundary) :
        add_edge!(plaquette(row, col - 1), plaquette(row, col))
    end

    distances = fill(-1, boundary, 2)
    queue = Tuple{Int,Bool}[(boundary, false)]
    distances[boundary, 1] = 0
    index = 1
    while index <= length(queue)
        node, parity = queue[index]
        index += 1
        node == boundary && parity && return distances[node, 2]
        for (next_node, edge_parity) in adjacency[node]
            next_parity = xor(parity, edge_parity)
            parity_index = next_parity ? 2 : 1
            distances[next_node, parity_index] >= 0 && continue
            distances[next_node, parity_index] = distances[node, parity ? 2 : 1] + 1
            push!(queue, (next_node, next_parity))
        end
    end
    throw(ErrorException("planar detector graph has no odd-logical path"))
end

function _isometric_planar_detector_batch(
        rng::Random.AbstractRNG, model::IsometricPlanarCodeModel,
        error_rate::Float64, shots::Int)
    distance = model.distance
    horizontal = rand(rng, distance, distance - 1, shots) .< error_rate
    vertical = rand(rng, distance - 1, distance, shots) .< error_rate
    detectors = falses(shots, _isometric_planar_detector_count(distance))
    actual = falses(shots)
    for shot in 1:shots
        for row in 1:(distance - 1), col in 1:(distance - 1)
            detector = _isometric_planar_plaquette_detector(distance, row, col)
            detectors[shot, detector + 1] = xor(
                horizontal[row, col, shot], horizontal[row + 1, col, shot],
                vertical[row, col, shot], vertical[row, col + 1, shot])
        end
        for col in 1:(distance - 1)
            detectors[shot, _isometric_planar_north_detector(distance, col) + 1] =
                horizontal[1, col, shot]
            detectors[shot, _isometric_planar_south_detector(distance, col) + 1] =
                horizontal[distance, col, shot]
        end
        parity = false
        for row in 1:(distance - 1)
            parity = xor(parity, vertical[row, 1, shot])
        end
        actual[shot] = parity
    end
    return detectors, actual
end

function _isometric_planar_batch_predictions(
        model::IsometricPlanarCodeModel, detectors::BitMatrix)
    shots = size(detectors, 1)
    size(detectors, 2) == _isometric_planar_detector_count(model.distance) ||
        throw(DimensionMismatch("detector batch has the wrong detector count"))
    if model.distance == 2
        return BitVector(any(detectors[shot, :]) for shot in 1:shots)
    end
    prediction = _isometric_planar_matching(model).decode_batch(
        pyimport("numpy").asarray(detectors))
    values = pyconvert(BitMatrix, prediction)
    expected = _isometric_planar_fault_count(model.distance) + 1
    size(values) == (shots, expected) || throw(ErrorException(
        "PyMatching returned predictions with size $(size(values)); expected ($shots, $expected)"))
    return BitVector(values[:, end])
end

"""One batched logical-failure estimate for the planar virtual-bond model."""
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

"""Finite-size-scaling fit for decoded planar virtual-bond capacity curves."""
struct IsometricPlanarScalingFit
    p_c::Union{Missing,Float64}
    nu::Union{Missing,Float64}
    p_c_se::Union{Missing,Float64}
    nu_se::Union{Missing,Float64}
    valid_bootstrap_fraction::Float64
    status::Symbol
    failure_window::Tuple{Float64,Float64}
    distances::Vector{Int}
    point_count::Int
    coefficients::Vector{Float64}
end

function _isometric_planar_scaling_unavailable(
        status::Symbol, window::Tuple{Float64,Float64}, distances::Vector{Int},
        point_count::Int)
    return IsometricPlanarScalingFit(
        missing, missing, missing, missing, 0.0, status, window, distances,
        point_count, Float64[])
end

function _isometric_planar_scaling_selected_points(
        scan::IsometricPlanarCapacityScan, window::Tuple{Float64,Float64})
    lower, upper = window
    return [point for point in scan.points if !point.diagnostic_only &&
            point.distance > 2 && lower <= point.logical_failure_rate <= upper]
end

struct IsometricPlanarScalingBootstrapSample
    replicate::Int
    status::Symbol
    input_rms_shift::Float64
    initial_p_c::Union{Missing,Float64}
    initial_nu::Union{Missing,Float64}
    p_c::Union{Missing,Float64}
    nu::Union{Missing,Float64}
    loss::Union{Missing,Float64}
    converged::Bool
    evaluations::Int
    boundary_hit::Bool
end

struct IsometricPlanarScalingBatchDiagnostic
    distance::Int
    error_rate::Float64
    batch_count::Int
    batch_size::Int
    batch_mean::Float64
    batch_std::Float64
    expected_batch_std::Float64
    std_ratio::Float64
    unique_batch_rates::Int
end

struct IsometricPlanarScalingSensitivity
    label::String
    status::Symbol
    p_c::Union{Missing,Float64}
    nu::Union{Missing,Float64}
    loss::Union{Missing,Float64}
    failure_window::Tuple{Float64,Float64}
    distances::Vector{Int}
    boundary_hit::Bool
end

"""Detailed evidence for auditing a finite-size-scaling fit."""
struct IsometricPlanarScalingDiagnostics
    fit::IsometricPlanarScalingFit
    bootstrap_samples::Vector{IsometricPlanarScalingBootstrapSample}
    loss_p_c::Vector{Float64}
    loss_nu::Vector{Float64}
    loss_values::Matrix{Float64}
    optimizer_p_c::Vector{Float64}
    optimizer_nu::Vector{Float64}
    optimizer_loss::Vector{Float64}
    nominal_initial_p_c::Union{Missing,Float64}
    nominal_initial_nu::Union{Missing,Float64}
    nominal_converged::Bool
    nominal_evaluations::Int
    nominal_boundary_hit::Bool
    batch_diagnostics::Vector{IsometricPlanarScalingBatchDiagnostic}
    sensitivities::Vector{IsometricPlanarScalingSensitivity}
end

struct _IsometricPlanarScalingSimplexer <: Optim.Simplexer
    minimum_rate::Float64
    maximum_rate::Float64
end

function Optim.simplexer(simplexer::_IsometricPlanarScalingSimplexer, initial)
    simplex = [copy(initial) for _ in 1:3]
    p_step = max((simplexer.maximum_rate - simplexer.minimum_rate) / 50, 1e-6)
    p_direction = initial[1] + p_step <= simplexer.maximum_rate ? 1.0 : -1.0
    simplex[2][1] += p_direction * p_step
    minimum_log_nu, maximum_log_nu = log(0.25), log(6.0)
    log_nu_step = (maximum_log_nu - minimum_log_nu) / 50
    nu_direction = initial[2] + log_nu_step <= maximum_log_nu ? 1.0 : -1.0
    simplex[3][2] += nu_direction * log_nu_step
    return simplex
end

function _isometric_planar_scaling_problem(
        points::Vector{IsometricPlanarCapacityPoint}, values::Vector{Float64})
    length(points) == length(values) || throw(DimensionMismatch(
        "scaling-fit values must match selected scan points"))
    length(points) >= 6 || return nothing
    all(isfinite, values) || return nothing
    rates = Float64[point.error_rate for point in points]
    minimum_rate, maximum_rate = extrema(rates)
    minimum_rate < maximum_rate || return nothing
    weights = Float64[1 / max(point.logical_failure_se, 1 / sqrt(point.shots))^2
                      for point in points]
    square_roots = sqrt.(weights)

    function objective(parameters)
        p_c, log_nu = parameters
        nu = exp(log_nu)
        minimum_rate <= p_c <= maximum_rate && 0.25 <= nu <= 6.0 || return Inf
        x = Float64[(point.error_rate - p_c) * point.distance^(1 / nu)
                    for point in points]
        design = hcat(ones(length(x)), x, x .^ 2, x .^ 3)
        rank(design) == 4 || return Inf
        coefficients = (square_roots .* design) \ (square_roots .* values)
        residuals = values .- design * coefficients
        return sum(weights .* residuals .^ 2)
    end

    return (; objective, minimum_rate, maximum_rate, square_roots)
end

function _isometric_planar_scaling_boundary_hit(
        p_c::Float64, nu::Float64, minimum_rate::Float64, maximum_rate::Float64)
    p_tolerance = max(1e-10, 1e-6 * (maximum_rate - minimum_rate))
    nu_tolerance = 1e-6 * (6.0 - 0.25)
    return abs(p_c - minimum_rate) <= p_tolerance ||
           abs(p_c - maximum_rate) <= p_tolerance ||
           abs(nu - 0.25) <= nu_tolerance || abs(nu - 6.0) <= nu_tolerance
end

function _isometric_planar_scaling_profile(
        points::Vector{IsometricPlanarCapacityPoint}, values::Vector{Float64};
        record_trace::Bool=false)
    problem = _isometric_planar_scaling_problem(points, values)
    problem === nothing && return nothing
    objective = problem.objective

    p_candidates = range(problem.minimum_rate, problem.maximum_rate; length=13)
    nu_candidates = range(0.5, 3.0; length=11)
    starts = [(p_c, log(nu)) for p_c in p_candidates, nu in nu_candidates]
    losses = [objective(start) for start in starts]
    initial = starts[argmin(losses)]
    initial_loss = minimum(losses)
    evaluations = Ref(0)
    best_loss = Ref(initial_loss)
    trace_p_c = Float64[initial[1]]
    trace_nu = Float64[exp(initial[2])]
    trace_loss = Float64[initial_loss]
    function optimized_objective(parameters)
        evaluations[] += 1
        loss = objective(parameters)
        if record_trace && isfinite(loss) && loss < best_loss[]
            best_loss[] = loss
            push!(trace_p_c, parameters[1])
            push!(trace_nu, exp(parameters[2]))
            push!(trace_loss, loss)
        end
        return loss
    end
    result = Optim.optimize(
        optimized_objective, collect(initial), Optim.NelderMead(
            initial_simplex=_IsometricPlanarScalingSimplexer(
                problem.minimum_rate, problem.maximum_rate)),
        Optim.Options(iterations=2_000, show_trace=false, store_trace=false))
    parameters = Optim.minimizer(result)
    loss = objective(parameters)
    isfinite(loss) || return nothing
    p_c, log_nu = parameters
    nu = exp(log_nu)
    x = Float64[(point.error_rate - p_c) * point.distance^(1 / nu)
                for point in points]
    design = hcat(ones(length(x)), x, x .^ 2, x .^ 3)
    coefficients = (problem.square_roots .* design) \
                   (problem.square_roots .* values)
    if record_trace && (trace_p_c[end] != p_c || trace_nu[end] != nu)
        push!(trace_p_c, p_c)
        push!(trace_nu, nu)
        push!(trace_loss, loss)
    end
    boundary_hit = _isometric_planar_scaling_boundary_hit(
        p_c, nu, problem.minimum_rate, problem.maximum_rate)
    return (; p_c, nu, coefficients, loss,
            initial_p_c=initial[1], initial_nu=exp(initial[2]),
            converged=Optim.converged(result), evaluations=evaluations[],
            boundary_hit, trace_p_c, trace_nu, trace_loss)
end

function _isometric_planar_scaling_window(failure_window)
    window = (Float64(failure_window[1]), Float64(failure_window[2]))
    all(isfinite, window) && 0 <= window[1] < window[2] <= 1 || throw(ArgumentError(
        "failure_window must satisfy 0 <= lower < upper <= 1"))
    return window
end

function _isometric_planar_scaling_fit_from_analysis(
        nominal, p_c_samples::Vector{Float64}, nu_samples::Vector{Float64},
        bootstrap::Int, window::Tuple{Float64,Float64}, distances::Vector{Int},
        point_count::Int)
    valid_fraction = length(p_c_samples) / bootstrap
    valid_fraction >= 0.8 || return IsometricPlanarScalingFit(
        nominal.p_c, nominal.nu, missing, missing, valid_fraction, :unstable,
        window, distances, point_count, nominal.coefficients)
    p_c_se = length(p_c_samples) >= 2 ? Statistics.std(p_c_samples) : 0.0
    nu_se = length(nu_samples) >= 2 ? Statistics.std(nu_samples) : 0.0
    return IsometricPlanarScalingFit(
        nominal.p_c, nominal.nu, p_c_se, nu_se, valid_fraction, :ok,
        window, distances, point_count, nominal.coefficients)
end

function _run_isometric_planar_scaling(
        rng::Random.AbstractRNG, scan::IsometricPlanarCapacityScan;
        bootstrap::Integer, failure_window, collect_diagnostics::Bool)
    bootstrap > 0 || throw(ArgumentError("bootstrap must be positive, got $bootstrap"))
    window = _isometric_planar_scaling_window(failure_window)
    points = _isometric_planar_scaling_selected_points(scan, window)
    distances = sort(unique(point.distance for point in points))
    if length(distances) < 3
        fit = _isometric_planar_scaling_unavailable(
            :insufficient_sizes, window, distances, length(points))
        return (; fit, points, nominal=nothing,
                bootstrap_samples=IsometricPlanarScalingBootstrapSample[])
    end
    if !all(length(point.logical_failure_batches) >= 2 for point in points)
        fit = _isometric_planar_scaling_unavailable(
            :insufficient_batches, window, distances, length(points))
        return (; fit, points, nominal=nothing,
                bootstrap_samples=IsometricPlanarScalingBootstrapSample[])
    end

    nominal_values = Float64[point.logical_failure_rate for point in points]
    nominal = _isometric_planar_scaling_profile(
        points, nominal_values; record_trace=collect_diagnostics)
    if nominal === nothing
        fit = _isometric_planar_scaling_unavailable(
            :fit_failed, window, distances, length(points))
        return (; fit, points, nominal=nothing,
                bootstrap_samples=IsometricPlanarScalingBootstrapSample[])
    end

    p_c_samples = Float64[]
    nu_samples = Float64[]
    records = IsometricPlanarScalingBootstrapSample[]
    for replicate in 1:Int(bootstrap)
        values = Float64[_bootstrap_batch_mean(rng, point.logical_failure_batches)
                         for point in points]
        differences = values .- nominal_values
        rms_shift = all(isfinite, differences) ?
            sqrt(sum(abs2, differences) / length(differences)) : NaN
        sample = _isometric_planar_scaling_profile(points, values)
        if sample === nothing
            collect_diagnostics && push!(records,
                IsometricPlanarScalingBootstrapSample(
                    replicate, :fit_failed, rms_shift, missing, missing,
                    missing, missing, missing, false, 0, false))
            continue
        end
        push!(p_c_samples, sample.p_c)
        push!(nu_samples, sample.nu)
        collect_diagnostics && push!(records,
            IsometricPlanarScalingBootstrapSample(
                replicate, :ok, rms_shift, sample.initial_p_c, sample.initial_nu,
                sample.p_c, sample.nu, sample.loss, sample.converged,
                sample.evaluations, sample.boundary_hit))
    end
    fit = _isometric_planar_scaling_fit_from_analysis(
        nominal, p_c_samples, nu_samples, Int(bootstrap), window,
        distances, length(points))
    return (; fit, points, nominal, bootstrap_samples=records)
end

"""
    fit_isometric_planar_scaling(rng, scan; bootstrap=2_000,
                                 failure_window=(0.05, 0.45))

Fit the finite-size-scaling form ``P_fail = F((p-p_c)d^(1/nu))`` to the
non-diagnostic points in the requested logical-failure transition window. The
master curve ``F`` is a cubic weighted least-squares polynomial. One-sigma
errors of ``p_c`` and ``nu`` are batch-bootstrap standard deviations.
"""
function fit_isometric_planar_scaling(
        rng::Random.AbstractRNG, scan::IsometricPlanarCapacityScan;
        bootstrap::Integer=2_000,
        failure_window::Tuple{<:Real,<:Real}=(0.05, 0.45))
    return _run_isometric_planar_scaling(
        rng, scan; bootstrap, failure_window,
        collect_diagnostics=false).fit
end

function _isometric_planar_batch_diagnostics(
        points::Vector{IsometricPlanarCapacityPoint})
    summaries = IsometricPlanarScalingBatchDiagnostic[]
    for point in points
        batches = point.logical_failure_batches
        isempty(batches) && continue
        batch_count = length(batches)
        batch_size = div(point.shots, batch_count)
        batch_mean = Statistics.mean(batches)
        batch_std = batch_count >= 2 ? Statistics.std(batches) : 0.0
        expected = sqrt(max(0.0, point.logical_failure_rate *
            (1 - point.logical_failure_rate) / batch_size))
        ratio = expected == 0 ? (batch_std == 0 ? 0.0 : Inf) : batch_std / expected
        push!(summaries, IsometricPlanarScalingBatchDiagnostic(
            point.distance, point.error_rate, batch_count, batch_size,
            batch_mean, batch_std, expected, ratio, length(unique(batches))))
    end
    return summaries
end

function _isometric_planar_sensitivity(
        label::String, points::Vector{IsometricPlanarCapacityPoint},
        window::Tuple{Float64,Float64})
    distances = sort(unique(point.distance for point in points))
    if length(distances) < 3
        return IsometricPlanarScalingSensitivity(
            label, :insufficient_sizes, missing, missing, missing,
            window, distances, false)
    end
    profile = _isometric_planar_scaling_profile(
        points, Float64[point.logical_failure_rate for point in points])
    profile === nothing && return IsometricPlanarScalingSensitivity(
        label, :fit_failed, missing, missing, missing, window, distances, false)
    return IsometricPlanarScalingSensitivity(
        label, :ok, profile.p_c, profile.nu, profile.loss, window,
        distances, profile.boundary_hit)
end

function _isometric_planar_sensitivities(
        scan::IsometricPlanarCapacityScan,
        baseline_points::Vector{IsometricPlanarCapacityPoint},
        window::Tuple{Float64,Float64})
    results = IsometricPlanarScalingSensitivity[
        _isometric_planar_sensitivity("baseline", baseline_points, window)]
    distances = sort(unique(point.distance for point in baseline_points))
    for distance in distances
        selected = [point for point in baseline_points if point.distance != distance]
        push!(results, _isometric_planar_sensitivity(
            "drop d=$distance", selected, window))
    end
    for alternative in ((0.05, 0.35), (0.10, 0.45))
        points = _isometric_planar_scaling_selected_points(scan, alternative)
        label = alternative == (0.05, 0.35) ?
            "window 0.05-0.35" : "window 0.10-0.45"
        push!(results, _isometric_planar_sensitivity(label, points, alternative))
    end
    return results
end

function _isometric_planar_loss_grid(
        points::Vector{IsometricPlanarCapacityPoint}, values::Vector{Float64},
        grid_size::Tuple{Int,Int})
    problem = _isometric_planar_scaling_problem(points, values)
    problem === nothing && return (Float64[], Float64[], zeros(0, 0))
    p_c_values = collect(range(
        problem.minimum_rate, problem.maximum_rate; length=grid_size[1]))
    nu_values = collect(range(0.25, 6.0; length=grid_size[2]))
    losses = Matrix{Float64}(undef, length(p_c_values), length(nu_values))
    for (p_index, p_c) in enumerate(p_c_values),
            (nu_index, nu) in enumerate(nu_values)
        losses[p_index, nu_index] = problem.objective((p_c, log(nu)))
    end
    return p_c_values, nu_values, losses
end

"""
    diagnose_isometric_planar_scaling(rng, scan; bootstrap=2_000,
        failure_window=(0.05, 0.45), loss_grid_size=(121, 121), sensitivity=true)

Fit the finite-size collapse and retain the batch, bootstrap, optimizer, loss
surface, and sensitivity evidence needed to audit the reported uncertainties.
"""
function diagnose_isometric_planar_scaling(
        rng::Random.AbstractRNG, scan::IsometricPlanarCapacityScan;
        bootstrap::Integer=2_000,
        failure_window::Tuple{<:Real,<:Real}=(0.05, 0.45),
        loss_grid_size::Tuple{<:Integer,<:Integer}=(121, 121),
        sensitivity::Bool=true)
    all(size -> size >= 2, loss_grid_size) || throw(ArgumentError(
        "loss_grid_size entries must both be at least 2"))
    grid_size = (Int(loss_grid_size[1]), Int(loss_grid_size[2]))
    analysis = _run_isometric_planar_scaling(
        rng, scan; bootstrap, failure_window, collect_diagnostics=true)
    batch_diagnostics = _isometric_planar_batch_diagnostics(analysis.points)
    if analysis.nominal === nothing
        return IsometricPlanarScalingDiagnostics(
            analysis.fit, analysis.bootstrap_samples, Float64[], Float64[],
            zeros(0, 0), Float64[], Float64[], Float64[], missing, missing,
            false, 0, false, batch_diagnostics,
            IsometricPlanarScalingSensitivity[])
    end
    nominal_values = Float64[point.logical_failure_rate for point in analysis.points]
    loss_p_c, loss_nu, loss_values = _isometric_planar_loss_grid(
        analysis.points, nominal_values, grid_size)
    sensitivities = sensitivity ? _isometric_planar_sensitivities(
        scan, analysis.points, analysis.fit.failure_window) :
        IsometricPlanarScalingSensitivity[]
    nominal = analysis.nominal
    return IsometricPlanarScalingDiagnostics(
        analysis.fit, analysis.bootstrap_samples, loss_p_c, loss_nu,
        loss_values, nominal.trace_p_c, nominal.trace_nu, nominal.trace_loss,
        nominal.initial_p_c, nominal.initial_nu, nominal.converged,
        nominal.evaluations, nominal.boundary_hit, batch_diagnostics,
        sensitivities)
end
