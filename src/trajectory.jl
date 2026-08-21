"""A lightweight diagonal schedule for measured toric-code trajectories."""
struct ToricCodeTrajectoryModel
    rows::Int
    cols::Int
    layers::Vector{Vector{Tuple{Int,Int}}}

    function ToricCodeTrajectoryModel(rows::Integer, cols::Integer)
        rows > 0 || throw(ArgumentError("rows must be positive, got $rows"))
        cols > 0 || throw(ArgumentError("cols must be positive, got $cols"))
        nrows, ncols = Int(rows), Int(cols)
        return new(nrows, ncols, _site_schedule(nrows, ncols))
    end
end

"""Independent bit-flip events on every internal and dangling virtual bond."""
struct VirtualBondErrors
    horizontal_internal::BitMatrix
    vertical_internal::BitMatrix
    west_boundary::BitVector
    east_boundary::BitVector
    south_boundary::BitVector
    north_boundary::BitVector
end

"""One measured trajectory, with physical records ordered `(E, N, W, S)`."""
struct ToricCodeTrajectory
    measurements::Array{Bool,3}
    errors::VirtualBondErrors
end

function Base.length(errors::VirtualBondErrors)
    return length(errors.horizontal_internal) +
           length(errors.vertical_internal) +
           length(errors.west_boundary) +
           length(errors.east_boundary) +
           length(errors.south_boundary) +
           length(errors.north_boundary)
end

function Base.count(errors::VirtualBondErrors)
    return count(errors.horizontal_internal) +
           count(errors.vertical_internal) +
           count(errors.west_boundary) +
           count(errors.east_boundary) +
           count(errors.south_boundary) +
           count(errors.north_boundary)
end

function _validate_error_rate(error_rate::Real)
    isfinite(error_rate) || throw(ArgumentError(
        "error_rate must be finite, got $error_rate"))
    0 <= error_rate <= 1 || throw(ArgumentError(
        "error_rate must be between 0 and 1, got $error_rate"))
    return float(error_rate)
end

_bernoulli_bits(rng::Random.AbstractRNG, dimensions::Tuple, probability::Real) =
    BitArray(rand(rng, dimensions...) .< probability)

"""Sample one independent bit-flip event on every virtual bond."""
function sample_virtual_errors(
        rng::Random.AbstractRNG, model::ToricCodeTrajectoryModel;
        error_rate::Real)
    probability = _validate_error_rate(error_rate)
    rows, cols = model.rows, model.cols
    return VirtualBondErrors(
        _bernoulli_bits(rng, (rows, cols - 1), probability),
        _bernoulli_bits(rng, (rows - 1, cols), probability),
        _bernoulli_bits(rng, (rows,), probability),
        _bernoulli_bits(rng, (rows,), probability),
        _bernoulli_bits(rng, (cols,), probability),
        _bernoulli_bits(rng, (cols,), probability),
    )
end

function _validate_virtual_error_shapes(
        rows::Int, cols::Int, errors::VirtualBondErrors)
    expected = (
        horizontal_internal=(rows, cols - 1),
        vertical_internal=(rows - 1, cols),
        west_boundary=(rows,),
        east_boundary=(rows,),
        south_boundary=(cols,),
        north_boundary=(cols,),
    )
    actual = (
        horizontal_internal=size(errors.horizontal_internal),
        vertical_internal=size(errors.vertical_internal),
        west_boundary=size(errors.west_boundary),
        east_boundary=size(errors.east_boundary),
        south_boundary=size(errors.south_boundary),
        north_boundary=size(errors.north_boundary),
    )
    for name in keys(expected)
        getproperty(actual, name) == getproperty(expected, name) ||
            throw(DimensionMismatch(
                "$name must have size $(getproperty(expected, name)), " *
                "got $(getproperty(actual, name))"))
    end
    return nothing
end

_validate_virtual_errors(
    model::ToricCodeTrajectoryModel, errors::VirtualBondErrors) =
    _validate_virtual_error_shapes(model.rows, model.cols, errors)

"""
    sample_trajectory(rng, model, errors)

Sample the physical `Z`-measurement record of the measured-and-reset circuit.
Boundary error bits are retained as metadata but leave the record unchanged
because the dangling virtual legs use `|+⟩` input and `⟨+|` output effects.
"""
function sample_trajectory(
        rng::Random.AbstractRNG, model::ToricCodeTrajectoryModel,
        errors::VirtualBondErrors)
    _validate_virtual_errors(model, errors)
    rows, cols = model.rows, model.cols
    measurements = Array{Bool}(undef, rows, cols, 4)
    horizontal_carriers = rand(rng, Bool, rows)
    vertical_carriers = rand(rng, Bool, cols)

    for layer in model.layers, (row, col) in layer
        gamma = horizontal_carriers[row]
        delta = vertical_carriers[col]
        alpha = rand(rng, Bool)
        beta = xor(alpha, gamma, delta)
        measurements[row, col, :] .= (alpha, beta, gamma, delta)

        horizontal_carriers[row] =
            col < cols ? xor(alpha, errors.horizontal_internal[row, col]) : alpha
        vertical_carriers[col] =
            row > 1 ? xor(beta, errors.vertical_internal[row - 1, col]) : beta
    end
    return ToricCodeTrajectory(measurements, errors)
end

function sample_trajectory(
        rng::Random.AbstractRNG, model::ToricCodeTrajectoryModel;
        error_rate::Real)
    errors = sample_virtual_errors(rng, model; error_rate=error_rate)
    return sample_trajectory(rng, model, errors)
end

function _apply_x!(register, wire_count::Int, wire::Int)
    Yao.apply!(register, repeat(wire_count, Yao.X, (wire,)))
    return register
end

"""
    sample_yao_trajectory(rng, model, errors; max_qubits=24)

Run the literal measured-and-reset Yao circuit using four reusable physical
ancillas and `rows + cols` carrier qubits. This reference implementation is
intended for validating small systems; `max_qubits` guards state-vector
allocation.
"""
function sample_yao_trajectory(
        rng::Random.AbstractRNG, model::ToricCodeTrajectoryModel,
        errors::VirtualBondErrors; max_qubits::Integer=24)
    _validate_virtual_errors(model, errors)
    max_qubits > 0 || throw(ArgumentError(
        "max_qubits must be positive, got $max_qubits"))

    rows, cols = model.rows, model.cols
    physical_wires = 1:4
    first_horizontal = 5
    first_vertical = first_horizontal + rows
    wire_count = 4 + rows + cols
    wire_count <= max_qubits || throw(ArgumentError(
        "measured Yao trajectory requires $wire_count qubits, " *
        "exceeding max_qubits=$max_qubits"))

    register = Yao.zero_state(wire_count)
    carrier_wires = first_horizontal:wire_count
    Yao.apply!(register, repeat(wire_count, H, carrier_wires))
    for row in 1:rows
        errors.west_boundary[row] &&
            _apply_x!(register, wire_count, first_horizontal + row - 1)
    end
    for col in 1:cols
        errors.south_boundary[col] &&
            _apply_x!(register, wire_count, first_vertical + col - 1)
    end

    measurements = Array{Bool}(undef, rows, cols, 4)
    local_matrix = reshape(toric_code_local_gate(ComplexF64), 64, 64)
    for layer in model.layers, (row, col) in layer
        horizontal_wire = first_horizontal + row - 1
        vertical_wire = first_vertical + col - 1
        local_gate = subroutine(
            wire_count, _yao_local_gate(local_matrix, (row, col)),
            (physical_wires..., horizontal_wire, vertical_wire))
        Yao.apply!(register, local_gate)
        outcome = Yao.measure!(
            Yao.ResetTo(0), register, physical_wires; rng=rng)
        for direction in 1:4
            measurements[row, col, direction] = Bool(outcome[direction])
        end

        col < cols && errors.horizontal_internal[row, col] &&
            _apply_x!(register, wire_count, horizontal_wire)
        row > 1 && errors.vertical_internal[row - 1, col] &&
            _apply_x!(register, wire_count, vertical_wire)
        col == cols && errors.east_boundary[row] &&
            _apply_x!(register, wire_count, horizontal_wire)
        row == 1 && errors.north_boundary[col] &&
            _apply_x!(register, wire_count, vertical_wire)
    end
    return ToricCodeTrajectory(measurements, errors)
end

log2_postselection_probability(model::ToricCodeTrajectoryModel) =
    -(model.rows + model.cols)

postselection_probability(model::ToricCodeTrajectoryModel) =
    exp2(float(log2_postselection_probability(model)))

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
