"""One elementary operation in a rotated-planar encoder schedule."""
struct EncoderOperation
    gate::Symbol
    qubits::Vector{Int}
end

"""A simultaneously addressable elementary layer and its literal support."""
struct EncoderGateLayer
    block_kind::Symbol
    block_order::Int
    operations::Vector{EncoderOperation}
    active_qubits::Vector{Int}
end

"""
A complete source-check block, or the explicitly marked logical/preparation
block, in a rotated-planar encoder.
"""
struct PlaquetteEncoderBlock
    kind::Symbol
    family::Symbol
    source_check::Symbol
    source_check_index::Int
    source_support::Vector{Int}
    reduced_support::Vector{Int}
    source_rows::Vector{Int}
    pivot::Int
    center2::Tuple{Int,Int}
    diagonal::Int
    geometric_key::NTuple{4,Int}
    geometric_order::Int
    direction::Symbol
    layers::Vector{EncoderGateLayer}
end

"""
    PlaquetteEncoder

Deterministic Clifford preparation schedule for a `RotatedPlanarCode`.
`blocks` retain the check geometry and elimination provenance; `layers` are
the exact elementary schedule replayed by [`yao_encoder`](@ref).
"""
struct PlaquetteEncoder
    code::RotatedPlanarCode
    logical_state::Symbol
    construction::Symbol
    blocks::Vector{PlaquetteEncoderBlock}
    layers::Vector{EncoderGateLayer}
end

function _check_center2(code::RotatedPlanarCode, support::Vector{Int})
    coordinates = data_qubit_coordinate.(Ref(code), support)
    rows = first.(coordinates)
    columns = last.(coordinates)
    return (minimum(rows) + maximum(rows), minimum(columns) + maximum(columns))
end

function _ordered_source_checks(code::RotatedPlanarCode, construction::Symbol)
    supports = construction === :bp ? b_p_checks(code) : a_s_checks(code)
    d = distance(code)
    records = [
        begin
            center2 = _check_center2(code, support)
            diagonal = 2d - center2[1] + center2[2]
            (index=index, support=support, center2=center2,
             key=(diagonal, center2[1], center2[2], index))
        end
        for (index, support) in enumerate(supports)
    ]
    sort!(records; by=record -> record.key)
    return records
end

function _bitrow(support::Vector{Int}, qubit_count::Int)
    row = falses(qubit_count)
    row[support] .= true
    return row
end

"""
Reduce independent source rows without changing their spatial order. The
provenance vectors record which original source checks enter each reduced
row. Later pivots are eliminated backwards, so the result is systematic in
the family-pivot columns while each row still spans the original family.
"""
function _ordered_rref(records, qubit_count::Int)
    row_count = length(records)
    rows = [_bitrow(record.support, qubit_count) for record in records]
    provenance = [falses(row_count) for _ in 1:row_count]
    for row in 1:row_count
        provenance[row][row] = true
    end
    pivots = Int[]

    for row in 1:row_count
        for previous in 1:(row - 1)
            if rows[row][pivots[previous]]
                rows[row] .= xor.(rows[row], rows[previous])
                provenance[row] .= xor.(provenance[row], provenance[previous])
            end
        end
        pivot = findfirst(rows[row])
        isnothing(pivot) && error("source checks must be linearly independent")
        push!(pivots, pivot)
        for previous in 1:(row - 1)
            if rows[previous][pivot]
                rows[previous] .= xor.(rows[previous], rows[row])
                provenance[previous] .= xor.(provenance[previous], provenance[row])
            end
        end
    end
    return rows, pivots, provenance
end

function _reduced_logical_row(
        support::Vector{Int}, qubit_count::Int,
        family_rows::Vector{BitVector}, family_pivots::Vector{Int},
        provenance::Vector{BitVector}, records)
    row = _bitrow(support, qubit_count)
    logical_provenance = falses(length(family_rows))
    for (source_row, (family_row, pivot)) in
        enumerate(zip(family_rows, family_pivots))
        if row[pivot]
            row .= xor.(row, family_row)
            logical_provenance .= xor.(logical_provenance, provenance[source_row])
        end
    end
    pivot = findfirst(row)
    isnothing(pivot) && error("logical generator must be independent of source checks")
    source_rows = sort!([
        records[row_index].index for row_index in eachindex(records)
        if logical_provenance[row_index]
    ])
    return row, pivot, source_rows
end

function _manhattan_path(code::RotatedPlanarCode, source::Int, destination::Int)
    source_row, source_column = data_qubit_coordinate(code, source)
    destination_row, destination_column = data_qubit_coordinate(code, destination)
    path = [source]
    row, column = source_row, source_column
    while column != destination_column
        column += destination_column > column ? 1 : -1
        push!(path, data_qubit_index(code, row, column))
    end
    while row != destination_row
        row += destination_row > row ? 1 : -1
        push!(path, data_qubit_index(code, row, column))
    end
    return path
end

_cnot(control_qubit::Int, target_qubit::Int) =
    EncoderOperation(:CNOT, [control_qubit, target_qubit])

function _swap_operations(first_qubit::Int, second_qubit::Int)
    return [
        _cnot(first_qubit, second_qubit),
        _cnot(second_qubit, first_qubit),
        _cnot(first_qubit, second_qubit),
    ]
end

"""Route one remote CNOT exactly by moving its control out and back."""
function _local_cnot_operations(
        code::RotatedPlanarCode, control_qubit::Int, target_qubit::Int)
    path = _manhattan_path(code, control_qubit, target_qubit)
    length(path) == 2 && return [_cnot(control_qubit, target_qubit)]

    operations = EncoderOperation[]
    for position in 1:(length(path) - 2)
        append!(operations, _swap_operations(path[position], path[position + 1]))
    end
    push!(operations, _cnot(path[end - 1], path[end]))
    for position in (length(path) - 2):-1:1
        append!(operations, _swap_operations(path[position], path[position + 1]))
    end
    return operations
end

function _layers(
        operations::Vector{EncoderOperation}, block_kind::Symbol, block_order::Int)
    return [
        EncoderGateLayer(
            block_kind, block_order, [operation], sort!(unique(copy(operation.qubits))))
        for operation in operations
    ]
end

function _row_cnot_operations(
        code::RotatedPlanarCode, row::BitVector, pivot::Int, construction::Symbol)
    operations = EncoderOperation[]
    for qubit in findall(row)
        qubit == pivot && continue
        control_qubit, target_qubit = construction === :bp ?
            (pivot, qubit) : (qubit, pivot)
        append!(operations,
                _local_cnot_operations(code, control_qubit, target_qubit))
    end
    return operations
end

function _logical_preparation_operations(
        logical_state::Symbol, logical_pivot::Int)
    logical_state === :zero && return EncoderOperation[]
    logical_state === :one && return [EncoderOperation(:X, [logical_pivot])]
    logical_state === :plus && return [EncoderOperation(:H, [logical_pivot])]
    return [
        EncoderOperation(:X, [logical_pivot]),
        EncoderOperation(:H, [logical_pivot]),
    ]
end

function _logical_block(
        code::RotatedPlanarCode, logical_state::Symbol, construction::Symbol,
        family_pivots::Vector{Int}, logical_row::BitVector, logical_pivot::Int,
        source_rows::Vector{Int})
    operations = EncoderOperation[]
    if construction === :as
        pivot_set = Set(vcat(family_pivots, logical_pivot))
        for qubit in 1:data_qubit_count(code)
            qubit in pivot_set || push!(operations, EncoderOperation(:H, [qubit]))
        end
    end
    append!(operations, _logical_preparation_operations(logical_state, logical_pivot))
    append!(operations,
            _row_cnot_operations(code, logical_row, logical_pivot, construction))
    direction = construction === :bp ? :outgoing_control : :incoming_control
    layers = _layers(operations, :logical_sector, 0)
    return PlaquetteEncoderBlock(
        :logical_sector, construction, :logical, 0,
        construction === :bp ? logical_x_support(code) : logical_z_support(code),
        findall(logical_row), source_rows, logical_pivot, (0, 0), 0,
        (typemin(Int), 0, 0, 0), 0, direction, layers)
end

function _plaquette_block(
        code::RotatedPlanarCode, construction::Symbol, record, row::BitVector,
        pivot::Int, source_rows::Vector{Int}, geometric_order::Int)
    operations = EncoderOperation[]
    construction === :bp && push!(operations, EncoderOperation(:H, [pivot]))
    append!(operations, _row_cnot_operations(code, row, pivot, construction))
    direction = construction === :bp ? :outgoing_control : :incoming_control
    layers = _layers(operations, :plaquette, geometric_order)
    return PlaquetteEncoderBlock(
        :plaquette, construction, construction === :bp ? :B_p : :A_s,
        record.index, copy(record.support), findall(row), source_rows, pivot,
        record.center2, record.key[1], record.key, geometric_order, direction, layers)
end

"""
    rotated_planar_encoder(code; logical_state=:zero, construction=:bp)

Construct a deterministic, Manhattan-local Clifford encoder. The B_p family
uses systematic X-check rows and outgoing-control CNOTs. The A_s family is
the Clifford dual: non-pivot inputs are prepared in the X basis and every
CNOT role is reversed. Both schedules prepare the same selected logical
state from the all-zero input.
"""
function rotated_planar_encoder(
        code::RotatedPlanarCode;
        logical_state::Symbol=:zero,
        construction::Symbol=:bp)::PlaquetteEncoder
    logical_state in (:zero, :one, :plus, :minus) ||
        throw(ArgumentError("logical_state must be :zero, :one, :plus, or :minus"))
    construction in (:as, :bp) ||
        throw(ArgumentError("construction must be :as or :bp"))

    records = _ordered_source_checks(code, construction)
    qubit_count = data_qubit_count(code)
    family_rows, family_pivots, provenance = _ordered_rref(records, qubit_count)
    logical_support = construction === :bp ?
        logical_x_support(code) : logical_z_support(code)
    logical_row, logical_pivot, logical_source_rows = _reduced_logical_row(
        logical_support, qubit_count, family_rows, family_pivots,
        provenance, records)

    blocks = PlaquetteEncoderBlock[
        _logical_block(
            code, logical_state, construction,
            family_pivots, logical_row, logical_pivot, logical_source_rows),
    ]
    for order in eachindex(records)
        source_rows = sort!([
            records[row].index for row in eachindex(records) if provenance[order][row]
        ])
        push!(blocks, _plaquette_block(
            code, construction, records[order], family_rows[order],
            family_pivots[order], source_rows, order))
    end
    layers = [layer for block in blocks for layer in block.layers]
    return PlaquetteEncoder(code, logical_state, construction, blocks, layers)
end

"""Return the stable preparation/check blocks in exact execution order."""
plaquette_blocks(encoder::PlaquetteEncoder) = encoder.blocks

"""Return the exact elementary layers in execution order."""
gate_layers(encoder::PlaquetteEncoder) = encoder.layers

function _yao_operation(qubit_count::Int, operation::EncoderOperation)
    operation.gate === :H && return put(qubit_count, only(operation.qubits) => H)
    operation.gate === :X && return put(qubit_count, only(operation.qubits) => X)
    operation.gate === :CNOT || error("unsupported encoder gate $(operation.gate)")
    return control(qubit_count, operation.qubits[1], operation.qubits[2] => X)
end

"""
    yao_encoder(encoder)

Return a Yao block on exactly `d^2` wires that replays the stored elementary
operation order without allocating a statevector.
"""
function yao_encoder(encoder::PlaquetteEncoder)
    qubit_count = data_qubit_count(encoder.code)
    blocks = [
        _yao_operation(qubit_count, operation)
        for layer in gate_layers(encoder) for operation in layer.operations
    ]
    return chain(qubit_count, blocks)
end

function _conjugate_pauli!(
        x::BitVector, z::BitVector, negative::Base.RefValue{Bool},
        operation::EncoderOperation)
    if operation.gate === :H
        qubit = only(operation.qubits)
        negative[] = xor(negative[], x[qubit] & z[qubit])
        x[qubit], z[qubit] = z[qubit], x[qubit]
    elseif operation.gate === :X
        qubit = only(operation.qubits)
        negative[] = xor(negative[], z[qubit])
    elseif operation.gate === :CNOT
        control_qubit, target_qubit = operation.qubits
        phase_flip = x[control_qubit] & z[target_qubit] &
                     xor(x[target_qubit], xor(z[control_qubit], true))
        negative[] = xor(negative[], phase_flip)
        x[target_qubit] = xor(x[target_qubit], x[control_qubit])
        z[control_qubit] = xor(z[control_qubit], z[target_qubit])
    else
        error("unsupported encoder gate $(operation.gate)")
    end
    return nothing
end

function _prepared_eigenvalue(
        encoder::PlaquetteEncoder, pauli::Symbol, support::Vector{Int})
    qubit_count = data_qubit_count(encoder.code)
    x = falses(qubit_count)
    z = falses(qubit_count)
    pauli === :X ? (x[support] .= true) : (z[support] .= true)
    negative = Ref(false)
    operations = [
        operation for layer in gate_layers(encoder) for operation in layer.operations
    ]
    for operation in Iterators.reverse(operations)
        _conjugate_pauli!(x, z, negative, operation)
    end
    any(x) && return 0
    return negative[] ? -1 : 1
end

"""
    verify_encoder_tableau(encoder) -> Bool

Conjugate every requested output stabilizer and logical observable backwards
through the exact stored Clifford schedule. The resulting operators must be
signed Z products with the required eigenvalues on the all-zero input. This
is an algebraic tableau check and never allocates a `2^n` statevector.
"""
function verify_encoder_tableau(encoder::PlaquetteEncoder)
    code = encoder.code
    all(_prepared_eigenvalue(encoder, :Z, support) == 1
        for support in a_s_checks(code)) || return false
    all(_prepared_eigenvalue(encoder, :X, support) == 1
        for support in b_p_checks(code)) || return false
    if encoder.logical_state in (:zero, :one)
        expected = encoder.logical_state === :zero ? 1 : -1
        return _prepared_eigenvalue(
            encoder, :Z, logical_z_support(code)) == expected
    end
    expected = encoder.logical_state === :plus ? 1 : -1
    return _prepared_eigenvalue(
        encoder, :X, logical_x_support(code)) == expected
end
