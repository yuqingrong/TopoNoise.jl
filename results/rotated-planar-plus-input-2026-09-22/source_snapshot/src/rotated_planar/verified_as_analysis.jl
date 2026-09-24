"""Independent d=3 statevector replay of verification and final ideal diagnostics."""
function sample_verified_as_yao(rng::Random.AbstractRNG, encoder::VerifiedAsEncoder,
                                record::VerifiedAsFaultRecord)
    distance(encoder.base.code) <= 3 ||
        throw(ArgumentError("the verification statevector oracle supports d=3 only"))
    _validate_verified_record(encoder, record)
    n = encoder.ancilla
    register = zero_state(n)
    outcomes = falses(length(encoder.checks))
    for (i, (step, fault)) in enumerate(zip(encoder.steps, record.quantum.steps))
        op = step.operation
        if op.gate === :RZ
            measure!(ResetTo(0), register, only(op.qubits); rng)
            _apply_yao_fault_step!(register, n, fault)
        elseif op.gate === :MZ
            bit = Int(measure!(ResetTo(0), register, only(op.qubits); rng)) == 1
            outcomes[step.check] = xor(bit, record.readout_flips[i])
        else
            apply!(register, _yao_operation(n, op))
            _apply_yao_fault_step!(register, n, fault)
        end
    end
    code = encoder.base.code
    a = BitVector(_measure_yao_check!(rng, register, n, n, support, :Z)
                  for support in a_s_checks(code))
    b = BitVector(_measure_yao_check!(rng, register, n, n, support, :X)
                  for support in b_p_checks(code))
    logical_x = _measure_yao_check!(rng, register, n, n, logical_z_support(code), :Z)
    return (accepted=!any(outcomes), outcomes=outcomes,
            syndrome=SyndromeRecord(a,b), logical_x=logical_x)
end

function _verified_logical_failures(code, decoders, x::BitMatrix, z::BitMatrix)
    predictions = decode_logical_parities(decoders,
        _support_parities(x, a_s_checks(code)), _support_parities(z, b_p_checks(code)))
    return (logical_x=xor.(_logical_parities(x, logical_z_support(code)), predictions.logical_x),
            logical_z=xor.(_logical_parities(z, logical_x_support(code)), predictions.logical_z))
end

"""
    audit_verified_as_single_faults(encoder; paulis=(:X,))

Enumerate each requested Pauli after every H/CNOT operand, every reset-X
fault and every readout flip. Report accepted failures after ideal final
syndrome recovery. The default covers the complete single-X-fault experiment;
it does not enumerate multi-operand correlated gate faults as single events.
"""
function audit_verified_as_single_faults(encoder::VerifiedAsEncoder; paulis=(:X,))
    !isempty(paulis) && length(unique(paulis)) == length(paulis) &&
        all(p -> p in (:X,:Y,:Z), paulis) ||
        throw(ArgumentError("paulis must contain distinct X, Y or Z symbols"))
    clean = sample_fault_record(Random.MersenneTwister(0), encoder, VerifiedAsNoise(0))
    locations = NamedTuple{(:step,:qubit,:pauli),Tuple{Int,Int,Symbol}}[]
    for (i, step) in enumerate(encoder.steps)
        op = step.operation
        if op.gate === :MZ
            push!(locations, (step=i, qubit=encoder.ancilla, pauli=:readout))
        elseif op.gate === :RZ
            push!(locations, (step=i, qubit=only(op.qubits), pauli=:X))
        else
            for q in op.qubits, axis in paulis
                push!(locations, (step=i, qubit=q, pauli=axis))
            end
        end
    end
    n = data_qubit_count(encoder.base.code)
    x, z = falses(length(locations), n), falses(length(locations), n)
    accepted = falses(length(locations))
    for (j, loc) in enumerate(locations)
        fault = loc.pauli === :readout ? with_measurement_fault(clean, loc.step) :
            with_pauli_fault(clean, loc.step, loc.qubit, loc.pauli)
        result = propagate_verified_as(encoder, fault)
        accepted[j] = result.accepted
        x[j,:], z[j,:] = result.frame.x, result.frame.z
    end
    failures = _verified_logical_failures(encoder.base.code,
        build_matching_decoders(encoder.base.code), x, z)
    cases = [merge(loc, (accepted=accepted[i], state_failure=failures.logical_x[i],
                        logical_z=failures.logical_z[i],
                        data_weight=count(x[i,:] .| z[i,:])))
             for (i, loc) in enumerate(locations)]
    return (distance=distance(encoder.base.code), rounds=encoder.rounds,
            fault_count=length(locations), accepted=count(accepted), rejected=count(.!accepted),
            accepted_state_failures=count(accepted .& failures.logical_x), cases=cases)
end

function _sample_verified_batch(rng, encoder::VerifiedAsEncoder, noise::VerifiedAsNoise,
                                shots::Int)
    x, z = falses(shots, encoder.ancilla), falses(shots, encoder.ancilla)
    accepted = trues(shots)
    for step in encoder.steps
        op = step.operation
        if op.gate === :RZ
            q = only(op.qubits)
            x[:,q] .= rand(rng, shots) .< noise.p_reset
            z[:,q] .= false
        elseif op.gate === :MZ
            q = only(op.qubits)
            flips = rand(rng, shots) .< noise.p_measure
            @views accepted .&= .!xor.(x[:,q], flips)
            x[:,q] .= false
            z[:,q] .= false
        else
            _propagate_frame_batch_operation!(x, z, op)
            _inject_random_frame_batch!(rng, x, z, op.qubits, noise.gates)
        end
    end
    n = data_qubit_count(encoder.base.code)
    return x[:,1:n], z[:,1:n], accepted
end

"""
    estimate_verified_as_failure(rng, encoder, noise; shots=10000, batch_size=1000, seed=nothing)

Estimate postselection yield and decoded-state failure among accepted shots.
Verification gates, resets and readouts are noisy; final syndrome recovery
uses the same ideal diagnostics as the original scans. Conditional statistics
are `nothing` if no attempts are accepted. Rejected runs are never silently
counted as successful preparations.
"""
function estimate_verified_as_failure(rng::Random.AbstractRNG, encoder::VerifiedAsEncoder,
                                      noise::VerifiedAsNoise; shots::Integer=10_000,
                                      batch_size::Integer=1000, seed=nothing)
    _assert_verified_schedule(encoder)
    count_shots = _positive_machine_int(shots, "shots")
    batch = _positive_machine_int(batch_size, "batch_size")
    seed_value = _seed_metadata(seed)
    code = encoder.base.code
    decoders = build_matching_decoders(code)
    accepted, failures, done = 0, 0, 0
    while done < count_shots
        count_batch = min(batch, count_shots-done)
        x, z, keep = _sample_verified_batch(rng, encoder, noise, count_batch)
        decoded = _verified_logical_failures(code, decoders, x, z)
        accepted += count(keep)
        failures += count(keep .& decoded.logical_x)
        done += count_batch
    end
    acceptance = accepted / count_shots
    conditional = accepted == 0 ? nothing : failures / accepted
    conditional_se = accepted == 0 ? nothing : sqrt(conditional*(1-conditional)/accepted)
    return (distance=distance(code), boundary_orientation=boundary_orientation(code),
            rounds=encoder.rounds, p_x=noise.gates.p_x, p_z=noise.gates.p_z,
            p_reset=noise.p_reset, p_measure=noise.p_measure, seed=seed_value,
            shots=count_shots, accepted=accepted, rejected=count_shots-accepted,
            acceptance_rate=acceptance,
            acceptance_standard_error=sqrt(acceptance*(1-acceptance)/count_shots),
            accepted_state_failures=failures, conditional_failure_rate=conditional,
            conditional_standard_error=conditional_se,
            accepted_failure_per_attempt=failures/count_shots)
end
