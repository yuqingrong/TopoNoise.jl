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
`layers` are the exact elementary schedule acting on the ideal product
`input_state` (`:zero` or `:plus`). [`yao_encoder`](@ref) can include the
ideal input preparation for application to an all-zero register.
"""
struct PlaquetteEncoder
    code::RotatedPlanarCode
    logical_state::Symbol
    construction::Symbol
    input_state::Symbol
    blocks::Vector{PlaquetteEncoderBlock}
    layers::Vector{EncoderGateLayer}
end

# Preserve the all-zero input convention for explicitly assembled schedules.
PlaquetteEncoder(code, logical_state, construction, blocks, layers) =
    PlaquetteEncoder(code, logical_state, construction, :zero, blocks, layers)

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

function _rotate_encoder_qubit_clockwise(code::RotatedPlanarCode, qubit::Int)
    row, column = data_qubit_coordinate(code, qubit)
    return data_qubit_index(code, column, distance(code) + 1 - row)
end

function _bp_d3_operations(code::RotatedPlanarCode)
    distance(code) == 3 || throw(ArgumentError("direct templates require d=3"))
    q(x, y) = data_qubit_index(code, distance(code) - y + 1, x)
    blocks = [
        [EncoderOperation(:H, [q(2, 1)]),
         _cnot(q(2, 1), q(1, 1)), _cnot(q(2, 1), q(2, 2)),
         _cnot(q(2, 2), q(1, 2))],
        [EncoderOperation(:H, [q(1, 3)]), _cnot(q(1, 3), q(1, 2))],
        [EncoderOperation(:H, [q(3, 1)]), _cnot(q(3, 1), q(3, 2))],
        [EncoderOperation(:H, [q(3, 3)]),
         _cnot(q(3, 3), q(3, 2)), _cnot(q(3, 3), q(2, 3)),
         _cnot(q(2, 3), q(2, 2))],
    ]
    boundary_orientation(code) === :x_ew && return blocks
    return [[EncoderOperation(
        operation.gate,
        _rotate_encoder_qubit_clockwise.(Ref(code), operation.qubits))
        for operation in block] for block in blocks]
end

function _bp_d3_core_blocks(code::RotatedPlanarCode)
    source_records = _ordered_source_checks(code, :bp)
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
        for operations in _bp_d3_operations(code)
    ]
    return PlaquetteEncoderBlock[
        PlaquetteEncoderBlock(
            :plaquette, :bp, :B_p,
            data.record.index, copy(data.record.support), data.representative,
            data.record.center2,
            data.record.key[1], data.record.key, order, :outgoing_control,
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

function _bp_tree_operations(
        code::RotatedPlanarCode, support::Vector{Int},
        occupied::Set{Int})
    fresh = sort!([qubit for qubit in support if !(qubit in occupied)])
    pair = _fresh_tree_pair(code, support, fresh)
    isnothing(pair) && error(
        "no fresh local tree for support $support at d=$(distance(code))")
    representative, continuation = pair
    operations = EncoderOperation[EncoderOperation(:H, [representative])]
    if length(support) == 2
        other = only(filter(qubit -> qubit != representative, support))
        push!(operations, _cnot(representative, other))
        return representative, operations
    end

    representative_neighbor = only(filter(
        qubit -> qubit != continuation &&
                  _adjacent_data_qubits(code, representative, qubit), support))
    continuation_neighbor = only(filter(
        qubit -> qubit != representative &&
                  _adjacent_data_qubits(code, continuation, qubit), support))
    append!(operations, [
        _cnot(representative, continuation),
        _cnot(representative, representative_neighbor),
        _cnot(continuation, continuation_neighbor),
    ])
    return representative, operations
end

function _bp_shelled_core_blocks(code::RotatedPlanarCode)
    supports = b_p_checks(code)
    records_by_index = Dict(
        record.index => record for record in _ordered_source_checks(code, :bp))
    order = _reverse_shell_order(code, supports)
    occupied = Set{Int}()
    blocks = PlaquetteEncoderBlock[]
    for (geometric_order, index) in enumerate(order)
        support = supports[index]
        representative, operations = _bp_tree_operations(code, support, occupied)
        record = records_by_index[index]
        push!(blocks, PlaquetteEncoderBlock(
            :plaquette, :bp, :B_p,
            index, copy(support), representative,
            record.center2, record.key[1], record.key, geometric_order,
            :outgoing_control, _layers(operations, :plaquette, geometric_order)))
        union!(occupied, support)
    end
    return blocks
end

function _local_core_blocks(code::RotatedPlanarCode, construction::Symbol)
    bp_blocks = distance(code) == 3 ? _bp_d3_core_blocks(code) :
                                     _bp_shelled_core_blocks(code)
    construction === :bp && return bp_blocks

    # One canonical schedule fixes the order and trees for both families.
    # R exchanges the check families; H⊗n reverses each CNOT's arrows while
    # preserving time order. Representative H gates are paired too when As
    # starts in all-plus and Bp starts in all-zero.
    rotation = _rotate_encoder_qubit_clockwise.(Ref(code), 1:data_qubit_count(code))
    source_by_support = Dict(
        Tuple(sort(record.support)) => record
        for record in _ordered_source_checks(code, :as))
    as_blocks = PlaquetteEncoderBlock[]
    for block in bp_blocks
        record = source_by_support[Tuple(sort(rotation[block.source_support]))]
        operations = EncoderOperation[
            operation.gate === :CNOT ?
                _cnot(rotation[operation.qubits[2]], rotation[operation.qubits[1]]) :
                EncoderOperation(operation.gate, rotation[operation.qubits])
            for layer in block.layers for operation in layer.operations]
        push!(as_blocks, PlaquetteEncoderBlock(
            :plaquette, :as, :A_s, record.index, copy(record.support),
            rotation[block.representative], record.center2, record.key[1],
            record.key, block.geometric_order, :incoming_control,
            _layers(operations, :plaquette, block.geometric_order)))
    end
    return as_blocks
end

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
        source_check in (:logical, :logical_parity) ? (typemin(Int), 0, 0, 0) :
                                                    (typemax(Int), 0, 0, 0),
        0, direction, _layers(operations, :logical_sector, 0))
end

function _as_cnot_blocks(code::RotatedPlanarCode)
    blocks = PlaquetteEncoderBlock[]
    touched = falses(data_qubit_count(code))
    for block in _local_core_blocks(code, :as)
        operations = _block_operations((block,))
        first_operation = first(operations)
        # A fresh representative's H commutes past every earlier block and
        # cancels its H from all-plus initialization. Leave that input in |0>.
        !touched[block.representative] && first_operation.gate === :H &&
            first_operation.qubits == [block.representative] &&
            all(op -> op.gate === :CNOT, operations[2:end]) ||
            error("As Hadamard cancellation requires a fresh representative and CNOT core")
        push!(blocks, PlaquetteEncoderBlock(
            block.kind, block.family, block.source_check, block.source_check_index,
            block.source_support, block.representative, block.center2, block.diagonal,
            block.geometric_key, block.geometric_order, block.direction,
            _layers(operations[2:end], block.kind, block.geometric_order)))
        touched[block.source_support] .= true
    end
    return blocks
end

function _as_input_parity(
        code::RotatedPlanarCode, blocks::Vector{PlaquetteEncoderBlock}, free::BitVector)
    # Pull logical Z back through C: C† Z_L C is a Z product. Factors on
    # representatives are already +1, so only the free inputs need even parity.
    z = falses(data_qubit_count(code))
    z[logical_z_support(code)] .= true
    for operation in Iterators.reverse(_block_operations(blocks))
        control_qubit, target_qubit = operation.qubits
        z[control_qubit] = xor(z[control_qubit], z[target_qubit])
    end
    support = findall(z .& free)
    isempty(support) && error("As logical Z must act on at least one free input")
    return support
end

function _parity_pivot(code::RotatedPlanarCode, support::Vector{Int})
    # Choose a central input deterministically. This reduces the seed's reach
    # but does not route its CNOTs; larger patches can have nonlocal seed edges.
    coordinates = data_qubit_coordinate.(Ref(code), support)
    return argmin(support) do qubit
        row, column = data_qubit_coordinate(code, qubit)
        lengths = [abs(row - other_row) + abs(column - other_column)
                   for (other_row, other_column) in coordinates]
        (maximum(lengths), sum(lengths), qubit)
    end
end

function _direct_as_blocks(code::RotatedPlanarCode, logical_state::Symbol)
    core_blocks = _as_cnot_blocks(code)
    free = trues(data_qubit_count(code))
    free[getproperty.(core_blocks, :representative)] .= false
    parity_support = logical_state in (:zero, :one) ?
        _as_input_parity(code, core_blocks, free) : Int[]
    pivot = isempty(parity_support) ? 0 : _parity_pivot(code, parity_support)
    initial_operations = EncoderOperation[
        EncoderOperation(:H, [qubit]) for qubit in findall(free) if qubit != pivot]
    blocks = PlaquetteEncoderBlock[_local_auxiliary_block(
        code, :as, :logical, logical_z_support(code),
        first(core_blocks).representative, :incoming_control, initial_operations)]
    if !isempty(parity_support)
        # The pivot starts at zero and accumulates the other parity bits.
        # Store this preparation separately from physical source-check blocks.
        push!(blocks, _local_auxiliary_block(
            code, :as, :logical_parity, parity_support, pivot, :incoming_control,
            EncoderOperation[_cnot(qubit, pivot) for qubit in parity_support if qubit != pivot]))
    end
    append!(blocks, core_blocks)
    return blocks
end

_logical_x_operations(code::RotatedPlanarCode) =
    EncoderOperation[EncoderOperation(:X, [qubit]) for qubit in logical_x_support(code)]

_logical_z_operations(code::RotatedPlanarCode) =
    EncoderOperation[EncoderOperation(:Z, [qubit]) for qubit in logical_z_support(code)]

function _local_encoder(
        code::RotatedPlanarCode, logical_state::Symbol, construction::Symbol)
    input_state = :zero
    if construction === :as
        if logical_state in (:plus, :minus)
            input_state = :plus
            core_blocks = _local_core_blocks(code, :as)
            blocks = PlaquetteEncoderBlock[_local_auxiliary_block(
                code, :as, :logical, logical_z_support(code),
                first(core_blocks).representative, :incoming_control, EncoderOperation[])]
            append!(blocks, core_blocks)
        else
            blocks = _direct_as_blocks(code, logical_state)
        end
    else
        bp_blocks = _local_core_blocks(code, :bp)
        bp_operations = _block_operations(bp_blocks)
        blocks = PlaquetteEncoderBlock[_local_auxiliary_block(
            code, :bp, :logical, logical_x_support(code),
            first(bp_blocks).representative, :outgoing_control, EncoderOperation[])]
        append!(blocks, bp_blocks)
        if logical_state in (:plus, :minus)
            conversion_operations = vcat(
                reverse(bp_operations),
                [EncoderOperation(:H, [qubit]) for qubit in 1:data_qubit_count(code)],
                _block_operations(_local_core_blocks(code, :as)),
            )
            push!(blocks, _local_auxiliary_block(
                code, :bp, :logical_conversion, logical_x_support(code),
                only(first(bp_operations).qubits), :logical_basis_conversion,
                conversion_operations))
        end
    end

    if logical_state === :one
        push!(blocks, _local_auxiliary_block(
            code, construction, :logical_x, logical_x_support(code),
            first(logical_x_support(code)), :logical_x_correction,
            _logical_x_operations(code)))
    elseif logical_state === :minus
        push!(blocks, _local_auxiliary_block(
            code, construction, :logical_z, logical_z_support(code),
            first(logical_z_support(code)), :logical_z_correction,
            _logical_z_operations(code)))
    end
    layers = [layer for block in blocks for layer in block.layers]
    return PlaquetteEncoder(code, logical_state, construction, input_state, blocks, layers)
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
    rotated_planar_encoder(code; logical_state=:native, construction=:bp)

Construct a deterministic Clifford encoder from directed source-check blocks.
The As CNOT core is the clockwise 90-degree rotation of the Bp core with
every CNOT's control and target exchanged, preserving the full execution
order, including boundary checks. The representative H gates are also
paired. This holds for every supported distance and both boundary orientations.
The default `:native` prepares `|+_L>` for `:as` and `|0_L>` for `:bp`,
starting from ideal all-plus and all-zero product inputs respectively.
Native As uses H on each fresh representative followed by its incoming
CNOT tree, with no parity block. These native noisy schedules are exactly
X/Z dual under the lattice rotation when fault rates are exchanged.
Explicit `:zero`/`:one` As requests retain all-zero input-parity preparation;
that `:logical_parity` block can contain nonlocal CNOTs. Explicit `:plus`/
`:minus` Bp requests use the logical-basis conversion schedule.
"""
function rotated_planar_encoder(
        code::RotatedPlanarCode;
        logical_state::Symbol=:native,
        construction::Symbol=:bp)::PlaquetteEncoder
    construction in (:as, :bp) ||
        throw(ArgumentError("construction must be :as or :bp"))
    logical_state === :native &&
        (logical_state = construction === :as ? :plus : :zero)
    logical_state in (:zero, :one, :plus, :minus) ||
        throw(ArgumentError("logical_state must be :native, :zero, :one, :plus, or :minus"))

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
    yao_encoder(encoder; prepare_input=true)

Return a Yao block on exactly `d^2` wires without allocating a statevector.
By default, prepend ideal input preparation so the block can be applied to
an all-zero register. This prefix is not part of `gate_layers` or the noisy
schedule. Set `prepare_input=false` to replay only the stored gate sequence
on a register already in `encoder.input_state`.
"""
function yao_encoder(encoder::PlaquetteEncoder; prepare_input::Bool=true)
    qubit_count = data_qubit_count(encoder.code)
    prefix = prepare_input ? _input_preparation_operations(encoder) : EncoderOperation[]
    operations = [op for layer in gate_layers(encoder) for op in layer.operations]
    blocks = [
        _yao_operation(qubit_count, operation)
        for operation in vcat(prefix, operations)
    ]
    return chain(qubit_count, blocks)
end

_input_preparation_operations(encoder::PlaquetteEncoder) =
    encoder.input_state === :plus ?
    [EncoderOperation(:H, [q]) for q in 1:data_qubit_count(encoder.code)] :
    EncoderOperation[]

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
    if encoder.input_state === :plus
        any(z) && return 0
    else
        any(x) && return 0
    end
    return negative[] ? -1 : 1
end

"""
    verify_encoder_tableau(encoder) -> Bool

Conjugate every requested output stabilizer and logical observable backwards
through the exact stored Clifford schedule. The resulting operators must be
signed Z products on all-zero input, or signed X products on all-plus input,
with the required eigenvalues. This is an algebraic tableau check and never
allocates a `2^n` statevector.
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
