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
    representative::Int
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
`blocks` retain their literal source-check geometry and local representative;
`layers` are the exact elementary schedule replayed by [`yao_encoder`](@ref).
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

function _rotate_xew_qubit_to_xns(code::RotatedPlanarCode, qubit::Int)
    row, column = data_qubit_coordinate(code, qubit)
    return data_qubit_index(code, column, distance(code) + 1 - row)
end

function _direct_d3_operations(code::RotatedPlanarCode, construction::Symbol)
    distance(code) == 3 || throw(ArgumentError("direct templates require d=3"))
    q(x, y) = data_qubit_index(code, distance(code) - y + 1, x)
    blocks = construction === :bp ? [
        [EncoderOperation(:H, [q(2, 1)]),
         _cnot(q(2, 1), q(1, 1)), _cnot(q(2, 1), q(2, 2)),
         _cnot(q(2, 2), q(1, 2))],
        [EncoderOperation(:H, [q(1, 3)]), _cnot(q(1, 3), q(1, 2))],
        [EncoderOperation(:H, [q(3, 1)]), _cnot(q(3, 1), q(3, 2))],
        [EncoderOperation(:H, [q(3, 3)]),
         _cnot(q(3, 3), q(3, 2)), _cnot(q(3, 3), q(2, 3)),
         _cnot(q(2, 3), q(2, 2))],
    ] : [
        [EncoderOperation(:H, [q(2, 1)]), _cnot(q(1, 1), q(2, 1))],
        [EncoderOperation(:H, [q(2, 2)]),
         _cnot(q(1, 2), q(2, 2)), _cnot(q(1, 3), q(2, 2)),
         _cnot(q(2, 3), q(2, 2))],
        [EncoderOperation(:H, [q(3, 1)]),
         _cnot(q(2, 1), q(3, 1)), _cnot(q(2, 2), q(3, 1)),
         _cnot(q(3, 2), q(3, 1))],
        [EncoderOperation(:H, [q(3, 3)]), _cnot(q(2, 3), q(3, 3))],
    ]
    boundary_orientation(code) === :x_ew && return blocks
    return [[EncoderOperation(
        operation.gate,
        _rotate_xew_qubit_to_xns.(Ref(code), operation.qubits))
        for operation in block] for block in blocks]
end

function _direct_d3_core_blocks(code::RotatedPlanarCode, construction::Symbol)
    source_records = _ordered_source_checks(code, construction)
    source_by_support = Dict(
        Tuple(sort(copy(record.support))) => record for record in source_records)
    block_data = [
        begin
            support = sort!(unique!([
                qubit for operation in operations for qubit in operation.qubits]))
            record = get(source_by_support, Tuple(support), nothing)
            isnothing(record) && error("direct d=3 template has no source check for $support")
            representative = only(first(operations).qubits)
            (record=record, representative=representative, operations=operations)
        end
        for operations in _direct_d3_operations(code, construction)
    ]
    direction = construction === :bp ? :outgoing_control : :incoming_control
    return PlaquetteEncoderBlock[
        PlaquetteEncoderBlock(
            :plaquette, construction, construction === :bp ? :B_p : :A_s,
            data.record.index, copy(data.record.support), data.representative,
            data.record.center2,
            data.record.key[1], data.record.key, order, direction,
            _layers(data.operations, :plaquette, order))
        for (order, data) in enumerate(block_data)
    ]
end

function _adjacent_data_qubits(code::RotatedPlanarCode, first_qubit::Int, second_qubit::Int)
    first_row, first_column = data_qubit_coordinate(code, first_qubit)
    second_row, second_column = data_qubit_coordinate(code, second_qubit)
    return abs(first_row - second_row) + abs(first_column - second_column) == 1
end

function _fresh_tree_pair(
        code::RotatedPlanarCode, support::Vector{Int}, fresh::Vector{Int})
    length(support) == 2 && return isempty(fresh) ? nothing : (first(fresh), 0)
    pairs = sort!([
        (first_qubit, second_qubit)
        for first_qubit in fresh for second_qubit in fresh
        if first_qubit != second_qubit &&
           _adjacent_data_qubits(code, first_qubit, second_qubit)
    ])
    return isempty(pairs) ? nothing : first(pairs)
end

function _reverse_shell_order(code::RotatedPlanarCode, supports::Vector{Vector{Int}})
    qubit_count = data_qubit_count(code)
    remaining = trues(length(supports))
    incidence = zeros(Int, qubit_count)
    for support in supports, qubit in support
        incidence[qubit] += 1
    end
    removed = Int[]

    while any(remaining)
        candidates = Tuple{Int,Int}[]
        for index in eachindex(supports)
            remaining[index] || continue
            fresh = [qubit for qubit in supports[index] if incidence[qubit] == 1]
            isnothing(_fresh_tree_pair(code, supports[index], fresh)) ||
                push!(candidates, (length(fresh), index))
        end
        isempty(candidates) && error(
            "local plaquette shelling stalled for d=$(distance(code))")
        sort!(candidates)
        _, index = first(candidates)
        push!(removed, index)
        remaining[index] = false
        for qubit in supports[index]
            incidence[qubit] -= 1
        end
    end
    return reverse(removed)
end

function _local_tree_operations(
        code::RotatedPlanarCode, construction::Symbol, support::Vector{Int},
        occupied::Set{Int})
    fresh = sort!([qubit for qubit in support if !(qubit in occupied)])
    pair = _fresh_tree_pair(code, support, fresh)
    isnothing(pair) && error(
        "no fresh local tree for support $support at d=$(distance(code))")
    representative, continuation = pair
    operations = EncoderOperation[EncoderOperation(:H, [representative])]
    if length(support) == 2
        other = only(filter(qubit -> qubit != representative, support))
        push!(operations, construction === :bp ?
              _cnot(representative, other) : _cnot(other, representative))
        return representative, operations
    end

    representative_neighbor = only(filter(
        qubit -> qubit != continuation &&
                  _adjacent_data_qubits(code, representative, qubit), support))
    continuation_neighbor = only(filter(
        qubit -> qubit != representative &&
                  _adjacent_data_qubits(code, continuation, qubit), support))
    append!(operations, construction === :bp ? [
        _cnot(representative, continuation),
        _cnot(representative, representative_neighbor),
        _cnot(continuation, continuation_neighbor),
    ] : [
        _cnot(continuation, representative),
        _cnot(representative_neighbor, representative),
        _cnot(continuation_neighbor, continuation),
    ])
    return representative, operations
end

function _reverse_shelled_core_blocks(code::RotatedPlanarCode, construction::Symbol)
    supports = construction === :bp ? b_p_checks(code) : a_s_checks(code)
    records_by_index = Dict(
        record.index => record for record in _ordered_source_checks(code, construction))
    order = _reverse_shell_order(code, supports)
    occupied = Set{Int}()
    direction = construction === :bp ? :outgoing_control : :incoming_control
    blocks = PlaquetteEncoderBlock[]
    for (geometric_order, index) in enumerate(order)
        support = supports[index]
        representative, operations = _local_tree_operations(
            code, construction, support, occupied)
        record = records_by_index[index]
        push!(blocks, PlaquetteEncoderBlock(
            :plaquette, construction, construction === :bp ? :B_p : :A_s,
            index, copy(support), representative,
            record.center2, record.key[1], record.key, geometric_order,
            direction, _layers(operations, :plaquette, geometric_order)))
        union!(occupied, support)
    end
    return blocks
end

_local_core_blocks(code::RotatedPlanarCode, construction::Symbol) =
    distance(code) == 3 ? _direct_d3_core_blocks(code, construction) :
                          _reverse_shelled_core_blocks(code, construction)

_block_operations(blocks) = [
    operation for block in blocks for layer in block.layers for operation in layer.operations
]

function _local_auxiliary_block(
        code::RotatedPlanarCode, construction::Symbol, source_check::Symbol,
        source_support::Vector{Int}, representative::Int, direction::Symbol,
        operations::Vector{EncoderOperation})
    return PlaquetteEncoderBlock(
        :logical_sector, construction, source_check, 0,
        copy(source_support), representative, (0, 0), 0,
        source_check === :logical ? (typemin(Int), 0, 0, 0) :
                                   (typemax(Int), 0, 0, 0),
        0, direction, _layers(operations, :logical_sector, 0))
end

function _local_logical_block(
        code::RotatedPlanarCode, construction::Symbol,
        as_operations::Vector{EncoderOperation}, bp_operations::Vector{EncoderOperation})
    operations = construction === :as ?
        EncoderOperation[EncoderOperation(:H, [qubit])
                         for qubit in 1:data_qubit_count(code)] :
        EncoderOperation[]
    representative = construction === :bp ?
        only(first(bp_operations).qubits) : only(first(as_operations).qubits)
    support = construction === :bp ? logical_x_support(code) : logical_z_support(code)
    return _local_auxiliary_block(
        code, construction, :logical, support, representative,
        construction === :bp ? :outgoing_control : :incoming_control,
        operations)
end

_logical_x_operations(code::RotatedPlanarCode) =
    EncoderOperation[EncoderOperation(:X, [qubit]) for qubit in logical_x_support(code)]

_logical_z_operations(code::RotatedPlanarCode) =
    EncoderOperation[EncoderOperation(:Z, [qubit]) for qubit in logical_z_support(code)]

function _local_encoder(
        code::RotatedPlanarCode, logical_state::Symbol, construction::Symbol)
    as_blocks = _local_core_blocks(code, :as)
    bp_blocks = _local_core_blocks(code, :bp)
    as_operations = _block_operations(as_blocks)
    bp_operations = _block_operations(bp_blocks)
    logical_block = _local_logical_block(
        code, construction, as_operations, bp_operations)
    blocks = PlaquetteEncoderBlock[logical_block]

    if construction === :bp
        append!(blocks, bp_blocks)
        if logical_state in (:plus, :minus)
            conversion_operations = vcat(
                reverse(bp_operations),
                [EncoderOperation(:H, [qubit]) for qubit in 1:data_qubit_count(code)],
                as_operations,
            )
            push!(blocks, _local_auxiliary_block(
                code, :bp, :logical_conversion, logical_x_support(code),
                only(first(bp_operations).qubits), :logical_basis_conversion,
                conversion_operations))
        elseif logical_state === :one
            push!(blocks, _local_auxiliary_block(
                code, :bp, :logical_x, logical_x_support(code),
                first(logical_x_support(code)), :logical_x_correction,
                _logical_x_operations(code)))
        end
    else
        append!(blocks, as_blocks)
        if logical_state in (:zero, :one)
            conversion_operations = vcat(
                reverse(as_operations),
                [EncoderOperation(:H, [qubit]) for qubit in 1:data_qubit_count(code)],
                bp_operations,
            )
            push!(blocks, _local_auxiliary_block(
                code, :as, :logical_conversion, logical_z_support(code),
                only(first(as_operations).qubits), :logical_basis_conversion,
                conversion_operations))
        end
        if logical_state === :one
            push!(blocks, _local_auxiliary_block(
                code, :as, :logical_x, logical_x_support(code),
                first(logical_x_support(code)), :logical_x_correction,
                _logical_x_operations(code)))
        end
    end

    if logical_state === :minus
        push!(blocks, _local_auxiliary_block(
            code, construction, :logical_z, logical_z_support(code),
            first(logical_z_support(code)), :logical_z_correction,
            _logical_z_operations(code)))
    end
    layers = [layer for block in blocks for layer in block.layers]
    return PlaquetteEncoder(code, logical_state, construction, blocks, layers)
end

_cnot(control_qubit::Int, target_qubit::Int) =
    EncoderOperation(:CNOT, [control_qubit, target_qubit])

function _layers(
        operations::Vector{EncoderOperation}, block_kind::Symbol, block_order::Int)
    return [
        EncoderGateLayer(
            block_kind, block_order, [operation], sort!(unique(copy(operation.qubits))))
        for operation in operations
    ]
end

"""
    rotated_planar_encoder(code; logical_state=:zero, construction=:bp)

Construct a deterministic local Clifford encoder from literal source-check
blocks. `:bp` and `:as` select different directed plaquette-growth circuits,
while construction-specific local basis conversion makes both prepare the
same selected logical state from the all-zero input.
"""
function rotated_planar_encoder(
        code::RotatedPlanarCode;
        logical_state::Symbol=:zero,
        construction::Symbol=:bp)::PlaquetteEncoder
    logical_state in (:zero, :one, :plus, :minus) ||
        throw(ArgumentError("logical_state must be :zero, :one, :plus, or :minus"))
    construction in (:as, :bp) ||
        throw(ArgumentError("construction must be :as or :bp"))

    return _local_encoder(code, logical_state, construction)
end

"""Return the stable preparation/check blocks in exact execution order."""
plaquette_blocks(encoder::PlaquetteEncoder) = encoder.blocks

"""Return the exact elementary layers in execution order."""
gate_layers(encoder::PlaquetteEncoder) = encoder.layers

function _yao_operation(qubit_count::Int, operation::EncoderOperation)
    operation.gate === :H && return put(qubit_count, only(operation.qubits) => H)
    operation.gate === :X && return put(qubit_count, only(operation.qubits) => X)
    operation.gate === :Z && return put(qubit_count, only(operation.qubits) => Z)
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
    elseif operation.gate === :Z
        qubit = only(operation.qubits)
        negative[] = xor(negative[], x[qubit])
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
