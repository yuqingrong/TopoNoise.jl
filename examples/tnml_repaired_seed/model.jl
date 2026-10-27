module RepairedSeedTN

using Random, SparseArrays, TopoNoise

export build_context, build_decoder, sample_record, decode_record, combine_posteriors

"""Group raw columns only after selecting the seed or transfer fault family."""
function topology_from_raw(detectors, logicals)
    D, B = sparse(UInt8.(detectors)), sparse(UInt8.(logicals))
    endpoints, masks, provenance = NTuple{2,Int}[], UInt8[], Vector{Int}[]
    groups = Dict{Tuple{NTuple{2,Int},UInt8},Int}()
    for column in axes(D, 2)
        rows = rowvals(D)[nzrange(D, column)]
        length(rows) == 2 || error("raw fault must have exactly two detector endpoints")
        edge = (rows[1], rows[2])
        mask = B[1,column] | (B[2,column] << 1)
        group = get(groups, (edge, mask), 0)
        if iszero(group)
            push!(endpoints, edge); push!(masks, mask); push!(provenance, [column])
            groups[(edge, mask)] = length(provenance)
        else
            push!(provenance[group], column)
        end
    end
    return Toric3DZTNTopology(D, B, endpoints, masks, length.(provenance), provenance)
end

function validate_rate(p::Real)
    isfinite(p) && 0 <= p < 0.5 || throw(ArgumentError("p must be finite and in [0,0.5)"))
    return Float64(p)
end

function binary_vector(values::AbstractVector, n::Int, name)
    length(values) == n || throw(DimensionMismatch("incorrect $name length"))
    all(x -> x == 0 || x == 1, values) || throw(ArgumentError("$name must be binary"))
    return BitVector(values)
end

logical_mask(bits) = Int(bits[1]) | (Int(bits[2]) << 1)
parity(matrix, bits) = BitVector(isodd.(matrix * UInt8.(bits)))

"""
    build_context(L)

Experiment-only context for noisy seed H outputs, unit-weight seed MWPM repair,
then diagonal H+6CNOT transfer faults and final readout. No intermediate feedback.
Both logical bits use the full final-code X checks, in mask order ell1+2ell2.
"""
function build_context(L::Integer)
    2 <= L <= typemax(Int) || throw(ArgumentError("L must be at least two and fit Int"))
    code = Toric3DEncodedSlice(L, L, L)
    model = build_toric3d_z_fault_model(code)
    full_topology = build_toric3d_z_tn_topology(model)
    family = build_toric3d_z_feedback_decoders(code)
    seed_decoder, final_decoder = family.seed, family.prefixes[end]
    seed_columns = findall(location -> location.kind === :seed, model.locations)
    transfer_columns = findall(location -> location.kind === :gate, model.locations)
    seed_qubits = [model.locations[column].qubit for column in seed_columns]
    seed_qubits == seed_decoder.window.qubits || error("seed qubit order differs from MWPM repair order")
    seed_logicals = full_topology.raw_logical_matrix[:,seed_columns]
    seed_logicals == final_decoder.logicals[:,seed_qubits] || error("seed logical labels changed during transfer")
    seed_topology = topology_from_raw(seed_decoder.matrix, seed_logicals)
    transfer_topology = topology_from_raw(
        full_topology.raw_detector_matrix[:,transfer_columns],
        full_topology.raw_logical_matrix[:,transfer_columns])
    return (;code, model, full_topology, seed_topology, transfer_topology,
             seed_decoder, final_decoder, seed_columns, transfer_columns,
             seed_logicals, seed_qubits)
end

"""Compile two reusable public TN decoders with the same declared settings."""
function build_decoder(context, p::Real; settings::Toric3DZTNSettings=Toric3DZTNSettings())
    rate = validate_rate(p)
    seed = build_toric3d_z_tn_decoder(context.seed_topology, rate; settings)
    transfer = build_toric3d_z_tn_decoder(context.transfer_topology, rate; settings)
    return (;context, p=rate, settings, seed, transfer)
end

"""
Draw all independent raw faults in established seed-then-gate-output order and
apply the actual seed recovery physically before computing the final syndrome.
Truth fields serve scoring only; decode_record never consumes them.
"""
function sample_record(context, p::Real, rng::Random.AbstractRNG)
    rate = validate_rate(p)
    raw_faults = BitVector(rand(rng) < rate for _ in axes(context.model.matrix, 2))
    seed_bits = raw_faults[context.seed_columns]
    seed_syndrome = parity(context.seed_decoder.matrix, seed_bits)
    seed_recovery = decode_toric3d_z_feedback(context.seed_decoder, seed_syndrome)
    corrected_seed = seed_bits .⊻ seed_recovery
    any(parity(context.seed_decoder.matrix, corrected_seed)) && error("seed repair left a syndrome")
    any(parity(context.full_topology.raw_detector_matrix[:,context.seed_columns], corrected_seed)) &&
        error("repaired seed contributes a final detector syndrome")
    frame = parity(context.model.matrix, raw_faults)
    frame[context.seed_qubits] .⊻= seed_recovery
    outcome = toric3d_z_threshold_outcome(context.final_decoder, frame)
    final_syndrome = outcome.syndrome
    true_mask = logical_mask(parity(context.final_decoder.logicals, frame))
    recovery_mask = logical_mask(parity(context.seed_logicals, seed_recovery))
    raw_seed_mask = logical_mask(parity(context.seed_logicals, seed_bits))
    transfer_mask = logical_mask(parity(context.transfer_topology.raw_logical_matrix,
                                        raw_faults[context.transfer_columns]))
    return (;raw_faults, seed_syndrome, seed_recovery, recovery_mask,
             final_syndrome, true_mask, mwpm_failure=outcome.failure,
             raw_seed_mask, transfer_mask)
end

"""
XOR-convolve four signed normalized weights and shift by the actual seed
recovery's logical mask. Entries stay in (00,10,01,11) order; no clipping.
"""
function combine_posteriors(seed::AbstractVector, transfer::AbstractVector, recovery_mask::Integer)
    length(seed) == length(transfer) == 4 || throw(DimensionMismatch("four logical weights are required"))
    0 <= recovery_mask <= 3 || throw(ArgumentError("recovery mask must be in 0:3"))
    all(isfinite, seed) && all(isfinite, transfer) || throw(ArgumentError("logical weights must be finite"))
    isapprox(sum(seed), 1; atol=1e-8, rtol=0) &&
        isapprox(sum(transfer), 1; atol=1e-8, rtol=0) || throw(ArgumentError("logical weights must be normalized"))
    combined = zeros(Float64, 4)
    for a in 0:3, b in 0:3
        combined[xor(xor(a, b), recovery_mask)+1] += seed[a+1]*transfer[b+1]
    end
    return combined
end

"""
Condition on observed initial seed and final syndromes, retaining the known
seed recovery. The MWPM baseline uses only the final syndrome after the common
seed repair. Numerical validity is not a convergence certificate. Exceptions
and signed weights from the public TN API are propagated without a fallback.
"""
function decode_record(decoder, record)
    context = decoder.context
    s0 = binary_vector(record.seed_syndrome, size(context.seed_decoder.matrix, 1), "seed syndrome")
    sf = binary_vector(record.final_syndrome, decoder.transfer.detector_count, "final syndrome")
    recovery = binary_vector(record.seed_recovery, length(context.seed_qubits), "seed recovery")
    parity(context.seed_decoder.matrix, recovery) == s0 || throw(ArgumentError("seed recovery does not match measured syndrome"))
    record.recovery_mask isa Integer && 0 <= record.recovery_mask <= 3 ||
        throw(ArgumentError("recovery mask must be an integer in 0:3"))
    record.recovery_mask == logical_mask(parity(context.seed_logicals, recovery)) ||
        throw(ArgumentError("recovery logical mask differs from the actual seed recovery"))
    seed_result = decode_toric3d_z_tn(decoder.seed, s0)
    transfer_result = decode_toric3d_z_tn(decoder.transfer, sf)
    probabilities = combine_posteriors(seed_result.probabilities,
        transfer_result.probabilities, record.recovery_mask)
    total = sum(probabilities)
    isfinite(total) && total > 0 || throw(DomainError(total, "invalid combined logical mass"))
    probabilities ./= total
    valid = seed_result.valid && transfer_result.valid && all(isfinite, probabilities) &&
        minimum(probabilities) >= -decoder.settings.negative_tolerance
    return (;valid, probabilities, logical_mask=valid ? argmax(probabilities)-1 : -1,
             log_joint_syndrome=seed_result.log_syndrome+transfer_result.log_syndrome,
             seed_result, transfer_result)
end

end # module
