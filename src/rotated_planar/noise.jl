"""
    CircuitPauliNoise(p; p_x=p, p_z=p, clock=:gate_layer)

Independent post-step Bernoulli X and Z faults on active data qubits. The
eligible steps are elementary encoder layers for `:gate_layer`, or completed
source-check blocks for `:plaquette`.
"""
struct CircuitPauliNoise
    p_x::Float64
    p_z::Float64
    clock::Symbol

    function CircuitPauliNoise(p_x::Float64, p_z::Float64, clock::Symbol)
        _validate_fault_rate(p_x, "p_x")
        _validate_fault_rate(p_z, "p_z")
        clock in (:gate_layer, :plaquette) ||
            throw(ArgumentError("clock must be :gate_layer or :plaquette"))
        return new(p_x, p_z, clock)
    end
end

function CircuitPauliNoise(
        p::Real; p_x::Real=p, p_z::Real=p, clock::Symbol=:gate_layer)
    _validate_fault_rate(p, "p")
    _validate_fault_rate(p_x, "p_x")
    _validate_fault_rate(p_z, "p_z")
    clock in (:gate_layer, :plaquette) ||
        throw(ArgumentError("clock must be :gate_layer or :plaquette"))
    return CircuitPauliNoise(Float64(p_x), Float64(p_z), clock)
end

function _validate_fault_rate(rate::Real, name::AbstractString)
    isfinite(rate) && 0 <= rate < 0.5 ||
        throw(ArgumentError("$name must be finite and satisfy 0 <= $name < 0.5"))
    return nothing
end

"""One ordered post-step data-qubit fault pattern."""
struct CircuitFaultStep
    step_index::Int
    block_kind::Symbol
    block_order::Int
    support::Tuple
    x::Tuple
    z::Tuple

    function CircuitFaultStep(
            step_index::Integer, block_kind::Symbol, block_order::Integer,
            support::Union{AbstractVector,Tuple},
            x::Union{AbstractVector,Tuple}, z::Union{AbstractVector,Tuple})
        step_index >= 1 || throw(ArgumentError("fault-step index must be positive"))
        block_kind in (:logical_sector, :plaquette) ||
            throw(ArgumentError("fault-step block kind is unsupported"))
        block_order >= 0 ||
            throw(ArgumentError("fault-step block order must be nonnegative"))
        all(qubit -> qubit isa Integer && qubit >= 1, support) ||
            throw(ArgumentError("fault-step support must contain positive integers"))
        length(unique(support)) == length(support) ||
            throw(ArgumentError("fault-step support must not contain duplicates"))
        length(x) == length(support) && length(z) == length(support) ||
            throw(ArgumentError("fault-step bit patterns must match its support"))
        all(bit -> bit isa Bool, x) && all(bit -> bit isa Bool, z) ||
            throw(ArgumentError("fault-step patterns must contain Bool values"))
        return new(
            Int(step_index), block_kind, Int(block_order),
            Tuple(Int(qubit) for qubit in support), Tuple(x), Tuple(z))
    end
end

"""An authoritative sampled circuit-fault record tied to one encoder."""
struct CircuitFaultRecord
    clock::Symbol
    distance::Int
    boundary_orientation::Symbol
    construction::Symbol
    logical_state::Symbol
    steps::Tuple

    function CircuitFaultRecord(
            clock::Symbol, distance::Integer, boundary_orientation::Symbol,
            construction::Symbol, logical_state::Symbol,
            steps::Union{AbstractVector,Tuple})
        clock in (:gate_layer, :plaquette) ||
            throw(ArgumentError("fault-record clock must be :gate_layer or :plaquette"))
        distance >= 3 && isodd(distance) ||
            throw(ArgumentError("fault-record distance must be an odd integer at least 3"))
        boundary_orientation in (:x_ns, :x_ew) ||
            throw(ArgumentError("fault-record boundary orientation is unsupported"))
        construction in (:as, :bp) ||
            throw(ArgumentError("fault-record construction is unsupported"))
        logical_state in (:zero, :one, :plus, :minus) ||
            throw(ArgumentError("fault-record logical state is unsupported"))
        all(step -> step isa CircuitFaultStep, steps) ||
            throw(ArgumentError("fault-record steps must be CircuitFaultStep values"))
        all(pair -> pair[2].step_index == pair[1], enumerate(steps)) ||
            throw(ArgumentError("fault-record step indices must be consecutive"))
        return new(
            clock, Int(distance), boundary_orientation, construction,
            logical_state, Tuple(steps))
    end
end

Base.:(==)(left::CircuitFaultStep, right::CircuitFaultStep) =
    left.step_index == right.step_index &&
    left.block_kind === right.block_kind &&
    left.block_order == right.block_order &&
    left.support == right.support && left.x == right.x && left.z == right.z

Base.:(==)(left::CircuitFaultRecord, right::CircuitFaultRecord) =
    left.clock === right.clock && left.distance == right.distance &&
    left.boundary_orientation === right.boundary_orientation &&
    left.construction === right.construction &&
    left.logical_state === right.logical_state && left.steps == right.steps

function _eligible_fault_steps(encoder::PlaquetteEncoder, clock::Symbol)
    if clock === :gate_layer
        return [
            (index, layer.block_kind, layer.block_order, copy(layer.active_qubits))
            for (index, layer) in enumerate(gate_layers(encoder))
        ]
    end
    blocks = filter(block -> block.kind === :plaquette, plaquette_blocks(encoder))
    return [
        (index, block.kind, block.geometric_order, copy(block.source_support))
        for (index, block) in enumerate(blocks)
    ]
end

"""
    sample_fault_record(rng, encoder, noise) -> CircuitFaultRecord

Sample all eligible post-step X/Z events once. The returned record is the
single source of truth for both binary-frame and Yao replay.
"""
function sample_fault_record(
        rng::Random.AbstractRNG, encoder::PlaquetteEncoder,
        noise::CircuitPauliNoise)
    steps = CircuitFaultStep[]
    for (index, block_kind, block_order, support) in
        _eligible_fault_steps(encoder, noise.clock)
        x = BitVector(rand(rng, length(support)) .< noise.p_x)
        z = BitVector(rand(rng, length(support)) .< noise.p_z)
        push!(steps, CircuitFaultStep(
            index, block_kind, block_order, support, x, z))
    end
    code = encoder.code
    return CircuitFaultRecord(
        noise.clock, distance(code), boundary_orientation(code),
        encoder.construction, encoder.logical_state, steps)
end

"""
    with_pauli_fault(record, step_index, qubit, pauli) -> CircuitFaultRecord

Return a new immutable record with `:X`, `:Z`, or `:Y` composed into one
eligible data-qubit event. This is the validated construction path for manual
fault fixtures; the input record is never modified.
"""
function with_pauli_fault(
        record::CircuitFaultRecord, step_index::Integer, qubit::Integer,
        pauli::Symbol)
    pauli in (:X, :Z, :Y) ||
        throw(ArgumentError("pauli must be :X, :Z, or :Y"))
    1 <= step_index <= length(record.steps) ||
        throw(ArgumentError("fault step index is outside the record"))
    step = record.steps[Int(step_index)]
    position = findfirst(==(Int(qubit)), step.support)
    isnothing(position) &&
        throw(ArgumentError("qubit is not active at the selected fault step"))

    x = collect(step.x)
    z = collect(step.z)
    pauli in (:X, :Y) && (x[position] = !x[position])
    pauli in (:Z, :Y) && (z[position] = !z[position])
    replacement = CircuitFaultStep(
        step.step_index, step.block_kind, step.block_order,
        step.support, x, z)
    steps = collect(record.steps)
    steps[Int(step_index)] = replacement
    return CircuitFaultRecord(
        record.clock, record.distance, record.boundary_orientation,
        record.construction, record.logical_state, steps)
end

"""Phase-free binary X/Z support on all data qubits."""
struct PauliFrame
    x::BitVector
    z::BitVector

    function PauliFrame(x::AbstractVector{Bool}, z::AbstractVector{Bool})
        length(x) == length(z) ||
            throw(ArgumentError("Pauli-frame X and Z vectors must have equal length"))
        return new(BitVector(x), BitVector(z))
    end
end

Base.:(==)(left::PauliFrame, right::PauliFrame) =
    left.x == right.x && left.z == right.z

"""Ordered A_s/Z-check and B_p/X-check measurement outcomes."""
struct SyndromeRecord
    a_s::BitVector
    b_p::BitVector

    function SyndromeRecord(
            a_s::AbstractVector{Bool}, b_p::AbstractVector{Bool})
        return new(BitVector(a_s), BitVector(b_p))
    end
end

Base.:(==)(left::SyndromeRecord, right::SyndromeRecord) =
    left.a_s == right.a_s && left.b_p == right.b_p

"""Return syndrome bits in deterministic A_s-then-B_p family order."""
syndrome_bits(syndrome::SyndromeRecord) = vcat(copy(syndrome.a_s), copy(syndrome.b_p))

function _validate_fault_record(
        encoder::PlaquetteEncoder, record::CircuitFaultRecord)
    code = encoder.code
    record.clock in (:gate_layer, :plaquette) ||
        throw(ArgumentError("fault-record clock must be :gate_layer or :plaquette"))
    record.distance == distance(code) ||
        throw(ArgumentError("fault record has a different code distance"))
    record.boundary_orientation === boundary_orientation(code) ||
        throw(ArgumentError("fault record has a different boundary orientation"))
    record.construction === encoder.construction ||
        throw(ArgumentError("fault record has a different encoder construction"))
    record.logical_state === encoder.logical_state ||
        throw(ArgumentError("fault record has a different logical state"))

    expected = _eligible_fault_steps(encoder, record.clock)
    length(record.steps) == length(expected) ||
        throw(ArgumentError("fault record has the wrong number of eligible steps"))
    for (step, (index, block_kind, block_order, support)) in
        zip(record.steps, expected)
        step.step_index == index && step.block_kind === block_kind &&
            step.block_order == block_order && step.support == Tuple(support) ||
            throw(ArgumentError("fault-record step identity or support is incompatible"))
        length(step.x) == length(support) && length(step.z) == length(support) ||
            throw(ArgumentError("fault-record bit patterns do not match step support"))
    end
    return nothing
end

function _propagate_frame_operation!(
        x::BitVector, z::BitVector, operation::EncoderOperation)
    if operation.gate === :H
        qubit = only(operation.qubits)
        x[qubit], z[qubit] = z[qubit], x[qubit]
    elseif operation.gate === :X
        nothing
    elseif operation.gate === :Z
        nothing
    elseif operation.gate === :CNOT
        control_qubit, target_qubit = operation.qubits
        x[target_qubit] = xor(x[target_qubit], x[control_qubit])
        z[control_qubit] = xor(z[control_qubit], z[target_qubit])
    else
        error("unsupported encoder gate $(operation.gate)")
    end
    return nothing
end

function _inject_frame_step!(
        x::BitVector, z::BitVector, step::CircuitFaultStep)
    for (position, qubit) in enumerate(step.support)
        x[qubit] = xor(x[qubit], step.x[position])
        z[qubit] = xor(z[qubit], step.z[position])
    end
    return nothing
end

"""
    propagate_pauli_frame(encoder, faults) -> PauliFrame

Replay the exact stored Clifford schedule and inject each supplied fault only
after its recorded eligible step.
"""
function propagate_pauli_frame(
        encoder::PlaquetteEncoder, faults::CircuitFaultRecord)
    _validate_fault_record(encoder, faults)
    qubit_count = data_qubit_count(encoder.code)
    x = falses(qubit_count)
    z = falses(qubit_count)

    if faults.clock === :gate_layer
        for (layer, step) in zip(gate_layers(encoder), faults.steps)
            for operation in layer.operations
                _propagate_frame_operation!(x, z, operation)
            end
            _inject_frame_step!(x, z, step)
        end
    else
        fault_index = 0
        for block in plaquette_blocks(encoder)
            for layer in block.layers, operation in layer.operations
                _propagate_frame_operation!(x, z, operation)
            end
            if block.kind === :plaquette
                fault_index += 1
                _inject_frame_step!(x, z, faults.steps[fault_index])
            end
        end
    end
    return PauliFrame(x, z)
end

"""
    measure_syndrome(code, frame) -> SyndromeRecord

Measure the phase-free CSS syndrome. A_s/Z checks detect frame X support and
B_p/X checks detect frame Z support.
"""
function measure_syndrome(code::RotatedPlanarCode, frame::PauliFrame)
    qubit_count = data_qubit_count(code)
    length(frame.x) == qubit_count ||
        throw(ArgumentError("Pauli frame must have $qubit_count data qubits"))
    a_s = BitVector(isodd(count(frame.x[support])) for support in a_s_checks(code))
    b_p = BitVector(isodd(count(frame.z[support])) for support in b_p_checks(code))
    return SyndromeRecord(a_s, b_p)
end

function _apply_yao_fault_step!(
        register, qubit_count::Int, step::CircuitFaultStep)
    for (position, qubit) in enumerate(step.support)
        step.x[position] && apply!(register, put(qubit_count, qubit => X))
        step.z[position] && apply!(register, put(qubit_count, qubit => Z))
    end
    return nothing
end

function _apply_yao_encoder_and_faults!(
        register, qubit_count::Int, encoder::PlaquetteEncoder,
        faults::CircuitFaultRecord)
    if faults.clock === :gate_layer
        for (layer, step) in zip(gate_layers(encoder), faults.steps)
            for operation in layer.operations
                apply!(register, _yao_operation(qubit_count, operation))
            end
            _apply_yao_fault_step!(register, qubit_count, step)
        end
    else
        fault_index = 0
        for block in plaquette_blocks(encoder)
            for layer in block.layers, operation in layer.operations
                apply!(register, _yao_operation(qubit_count, operation))
            end
            if block.kind === :plaquette
                fault_index += 1
                _apply_yao_fault_step!(
                    register, qubit_count, faults.steps[fault_index])
            end
        end
    end
    return nothing
end

function _measure_yao_check!(
        rng::Random.AbstractRNG, register, qubit_count::Int, ancilla::Int,
        support, pauli::Symbol)
    if pauli === :Z
        for data_qubit in support
            apply!(register, control(qubit_count, data_qubit, ancilla => X))
        end
    else
        apply!(register, put(qubit_count, ancilla => H))
        for data_qubit in support
            apply!(register, control(qubit_count, ancilla, data_qubit => X))
        end
        apply!(register, put(qubit_count, ancilla => H))
    end
    return Int(measure!(ResetTo(0), register, ancilla; rng=rng)) == 1
end

"""
    sample_yao_syndrome(rng, code, encoder, faults) -> SyndromeRecord

Small-distance statevector oracle that replays the supplied fault record and
measures explicit ideal ancilla extraction circuits. One ancilla is reset and
reused after every check.
"""
function sample_yao_syndrome(
        rng::Random.AbstractRNG, code::RotatedPlanarCode,
        encoder::PlaquetteEncoder, faults::CircuitFaultRecord)
    distance(code) <= 3 ||
        throw(ArgumentError("Yao syndrome extraction supports distance at most 3"))
    distance(encoder.code) == distance(code) &&
        boundary_orientation(encoder.code) === boundary_orientation(code) ||
        throw(ArgumentError("code and encoder geometry are incompatible"))
    _validate_fault_record(encoder, faults)

    data_count = data_qubit_count(code)
    ancilla = data_count + 1
    qubit_count = ancilla
    register = zero_state(qubit_count)
    _apply_yao_encoder_and_faults!(register, qubit_count, encoder, faults)

    a_s = BitVector(
        _measure_yao_check!(rng, register, qubit_count, ancilla, support, :Z)
        for support in a_s_checks(code))
    b_p = BitVector(
        _measure_yao_check!(rng, register, qubit_count, ancilla, support, :X)
        for support in b_p_checks(code))
    return SyndromeRecord(a_s, b_p)
end
