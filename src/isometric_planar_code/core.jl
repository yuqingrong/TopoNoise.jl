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

