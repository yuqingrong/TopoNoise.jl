"""One original encoder operation or an explicit verification operation."""
struct VerifiedAsStep
    operation::EncoderOperation
    source_gate::Int
    check::Int
end

"""A Z-product readout with known ideal eigenvalue +1 at its checkpoint."""
struct AsVerificationCheck
    kind::Symbol
    support::Vector{Int}
    after_gate::Int
    checkpoint::Int
    round::Int
end

"""Postselected As preparation of known logical zero, using one reused ancilla."""
struct VerifiedAsEncoder
    base::PlaquetteEncoder
    ancilla::Int
    rounds::Int
    steps::Vector{VerifiedAsStep}
    checks::Vector{AsVerificationCheck}
    signature::String
end

function _verified_schedule_signature(base, ancilla, rounds, steps, checks)
    io = IOBuffer()
    code = base.code
    print(io, distance(code), ':', boundary_orientation(code), ':', base.construction,
          ':', base.logical_state, ':', ancilla, ':', rounds)
    for supports in (a_s_checks(code), b_p_checks(code))
        print(io, '|')
        for support in supports
            print(io, join(support, ','), ';')
        end
    end
    print(io, "|base|")
    for operation in _block_operations(base.blocks)
        print(io, operation.gate, ':', join(operation.qubits, ','), ';')
    end
    print(io, "|layers|")
    for layer in base.layers, operation in layer.operations
        print(io, layer.block_kind, ':', layer.block_order, ':', operation.gate,
              ':', join(operation.qubits, ','), ';')
    end
    print(io, "|steps|")
    for step in steps
        print(io, step.operation.gate, ':', join(step.operation.qubits, ','),
              ':', step.source_gate, ':', step.check, ';')
    end
    print(io, "|checks|")
    for check in checks
        print(io, check.kind, ':', join(check.support, ','), ':', check.after_gate,
              ':', check.checkpoint, ':', check.round, ';')
    end
    return String(take!(io))
end

_verified_schedule_intact(encoder::VerifiedAsEncoder) =
    encoder.signature == _verified_schedule_signature(encoder.base, encoder.ancilla,
        encoder.rounds, encoder.steps, encoder.checks)

function _assert_verified_schedule(encoder::VerifiedAsEncoder)
    _verified_schedule_intact(encoder) ||
        throw(ArgumentError("verification schedule was modified; rebuild the encoder"))
    return nothing
end

"""
    verified_as_encoder(code; rounds=2)

Keep the original As logical-zero gates and insert explicit Z-check
verification after parity preparation, after each As block, and at the end.
All recorded checks must be +1 for acceptance. Each check uses reset, incoming
data-to-ancilla CNOTs and Z readout. Initial data resets are explicit too.

This prepares a known state with postselection. It is not an encoder for an
unknown logical qubit, and the ancilla interactions are not SWAP-routed.
"""
function verified_as_encoder(code::RotatedPlanarCode; rounds::Integer=2)
    repetitions = _positive_machine_int(rounds, "rounds")
    base = rotated_planar_encoder(code; construction=:as, logical_state=:zero)
    n = data_qubit_count(code)
    ancilla = n + 1
    steps = [VerifiedAsStep(EncoderOperation(:RZ, [q]), 0, 0) for q in 1:n]
    checks = AsVerificationCheck[]
    parity = only(filter(b -> b.source_check === :logical_parity, base.blocks))
    logical_generator = falses(n)
    logical_generator[parity.representative] = true
    source_gate = 0
    checkpoint = 0

    function add_checkpoint(entries)
        checkpoint += 1
        for round in 1:repetitions, (kind, support) in entries
            push!(checks, AsVerificationCheck(
                kind, copy(support), source_gate, checkpoint, round))
            index = length(checks)
            push!(steps, VerifiedAsStep(EncoderOperation(:RZ, [ancilla]), 0, index))
            for q in support
                push!(steps, VerifiedAsStep(_cnot(q, ancilla), 0, index))
            end
            push!(steps, VerifiedAsStep(EncoderOperation(:MZ, [ancilla]), 0, index))
        end
    end

    for block in base.blocks
        for layer in block.layers, operation in layer.operations
            source_gate += 1
            push!(steps, VerifiedAsStep(operation, source_gate, 0))
            if operation.gate === :CNOT
                c, t = operation.qubits
                logical_generator[c] = xor(logical_generator[c], logical_generator[t])
            elseif operation.gate === :H
                logical_generator[only(operation.qubits)] &&
                    error("verification generator must stay Z-type")
            end
        end
        if block.source_check === :logical_parity
            add_checkpoint([(:input_parity, findall(logical_generator))])
        elseif block.kind === :plaquette
            add_checkpoint([(:A_s, block.source_support),
                            (:growing_logical, findall(logical_generator))])
        end
    end
    final_checks = [(:A_s, support) for support in a_s_checks(code)]
    push!(final_checks, (:logical_z, logical_z_support(code)))
    add_checkpoint(final_checks)
    signature = _verified_schedule_signature(base, ancilla, repetitions, steps, checks)
    encoder = VerifiedAsEncoder(base, ancilla, repetitions, steps, checks, signature)
    verify_encoder_tableau(encoder) || error("verification measures a nontrivial ideal outcome")
    return encoder
end

"""Verify every checkpoint on the original ideal prefix, plus the final state."""
function verify_encoder_tableau(encoder::VerifiedAsEncoder)
    _verified_schedule_intact(encoder) || return false
    verify_encoder_tableau(encoder.base) || return false
    n = data_qubit_count(encoder.base.code)
    operations = _block_operations(encoder.base.blocks)
    for check in encoder.checks
        x, z, negative = falses(n), falses(n), Ref(false)
        z[check.support] .= true
        for i in check.after_gate:-1:1
            _conjugate_pauli!(x, z, negative, operations[i])
        end
        (any(x) || negative[]) && return false
    end
    return true
end

"""
    VerifiedAsNoise(p; p_x=p, p_z=0, p_reset=p, p_measure=p)

Independent X/Z faults after each H/CNOT operand, X faults after data and
ancilla resets, and classical readout flips. No idle faults are sampled.
The default is X-only gate noise. Final recovery is ideal in the evaluation.
"""
struct VerifiedAsNoise
    gates::CircuitPauliNoise
    p_reset::Float64
    p_measure::Float64
    function VerifiedAsNoise(p::Real; p_x::Real=p, p_z::Real=0,
                             p_reset::Real=p, p_measure::Real=p)
        gates = CircuitPauliNoise(p; p_x, p_z, clock=:gate_layer)
        _validate_fault_rate(p_reset, "p_reset")
        _validate_fault_rate(p_measure, "p_measure")
        return new(gates, Float64(p_reset), Float64(p_measure))
    end
end

"""Schedule-bound quantum fault steps and separately recorded readout flips."""
struct VerifiedAsFaultRecord
    signature::String
    quantum::CircuitFaultRecord
    readout_flips::Tuple
end

function _verified_step_identity(encoder::VerifiedAsEncoder, step::VerifiedAsStep)
    if step.source_gate > 0
        layer = encoder.base.layers[step.source_gate]
        return layer.block_kind, layer.block_order, copy(step.operation.qubits)
    end
    support = step.operation.gate === :MZ ? Int[] : copy(step.operation.qubits)
    return :logical_sector, step.check, support
end

function sample_fault_record(rng::Random.AbstractRNG, encoder::VerifiedAsEncoder,
                             noise::VerifiedAsNoise)
    _assert_verified_schedule(encoder)
    faults = CircuitFaultStep[]
    flips = Bool[]
    for (i, step) in enumerate(encoder.steps)
        gate = step.operation.gate
        kind, order, support = _verified_step_identity(encoder, step)
        px = gate === :RZ ? noise.p_reset : noise.gates.p_x
        pz = gate === :RZ ? 0.0 : noise.gates.p_z
        x = Tuple(rand(rng) < px for _ in support)
        z = Tuple(rand(rng) < pz for _ in support)
        push!(faults, CircuitFaultStep(i, kind, order, support, x, z))
        push!(flips, gate === :MZ && rand(rng) < noise.p_measure)
    end
    code = encoder.base.code
    quantum = CircuitFaultRecord(:gate_layer, distance(code), boundary_orientation(code),
                                  :as, :zero, faults)
    return VerifiedAsFaultRecord(encoder.signature, quantum, Tuple(flips))
end

function with_pauli_fault(record::VerifiedAsFaultRecord, step::Integer,
                          qubit::Integer, pauli::Symbol)
    return VerifiedAsFaultRecord(record.signature,
        with_pauli_fault(record.quantum, step, qubit, pauli), record.readout_flips)
end

"""Toggle one classical verification readout; `step` is a full schedule index."""
function with_measurement_fault(record::VerifiedAsFaultRecord, step::Integer)
    1 <= step <= length(record.quantum.steps) ||
        throw(ArgumentError("measurement step is outside the record"))
    isempty(record.quantum.steps[step].support) ||
        throw(ArgumentError("the selected step is not a measurement"))
    flips = collect(record.readout_flips)
    flips[step] = !flips[step]
    return VerifiedAsFaultRecord(record.signature, record.quantum, Tuple(flips))
end

function _validate_verified_record(encoder::VerifiedAsEncoder, record::VerifiedAsFaultRecord)
    _assert_verified_schedule(encoder)
    record.signature == encoder.signature ||
        throw(ArgumentError("fault record belongs to a different verification schedule"))
    code = encoder.base.code
    quantum = record.quantum
    quantum.clock === :gate_layer && quantum.distance == distance(code) &&
        quantum.boundary_orientation === boundary_orientation(code) &&
        quantum.construction === :as && quantum.logical_state === :zero ||
        throw(ArgumentError("fault record has incompatible encoder metadata"))
    length(quantum.steps) == length(record.readout_flips) == length(encoder.steps) ||
        throw(ArgumentError("fault record has the wrong number of steps"))
    all(b -> b isa Bool, record.readout_flips) ||
        throw(ArgumentError("readout flips must be Boolean"))
    for (i, (step, fault)) in enumerate(zip(encoder.steps, quantum.steps))
        kind, order, support = _verified_step_identity(encoder, step)
        fault.step_index == i && fault.block_kind === kind && fault.block_order == order &&
            fault.support == Tuple(support) ||
            throw(ArgumentError("fault step does not match the verification schedule"))
        step.operation.gate === :MZ || !record.readout_flips[i] ||
            throw(ArgumentError("readout flip attached to a quantum operation"))
    end
    return nothing
end

"""
    propagate_verified_as(encoder, faults)

Replay explicit gates, reset, measurement and supplied faults. Acceptance
uses only recorded check outcomes. The returned data frame is for subsequent
ideal recovery/evaluation and is not used to decide acceptance.
"""
function propagate_verified_as(encoder::VerifiedAsEncoder, record::VerifiedAsFaultRecord)
    _validate_verified_record(encoder, record)
    x, z = falses(encoder.ancilla), falses(encoder.ancilla)
    outcomes = falses(length(encoder.checks))
    for (i, (step, fault)) in enumerate(zip(encoder.steps, record.quantum.steps))
        op = step.operation
        if op.gate === :RZ
            q = only(op.qubits)
            x[q] = z[q] = false
            _inject_frame_step!(x, z, fault)
        elseif op.gate === :MZ
            q = only(op.qubits)
            outcomes[step.check] = xor(x[q], record.readout_flips[i])
            x[q] = z[q] = false
        else
            _propagate_frame_operation!(x, z, op)
            _inject_frame_step!(x, z, fault)
        end
    end
    n = data_qubit_count(encoder.base.code)
    return (accepted=!any(outcomes), outcomes=outcomes,
            frame=PauliFrame(x[1:n], z[1:n]))
end
