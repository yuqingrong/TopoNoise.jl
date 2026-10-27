function _toric3d_validate_h6_full_prefix_family(encoder::Toric3DSliceEncoder,
        family::Union{Toric3DXFeedbackDecoderFamily,Toric3DZFeedbackDecoderFamily})
    lattice = encoder.code
    dims = (lattice.steps, lattice.Ly, lattice.Lz)
    code = family.code.lattice
    dims == (code.steps, code.Ly, code.Lz) ||
        throw(DimensionMismatch("H+6 full-prefix decoder and encoder dimensions differ"))
    length(family.prefixes) == lattice.steps ||
        throw(DimensionMismatch("H+6 full-prefix decoder family has the wrong depth"))
    for (r, decoder) in enumerate(family.prefixes)
        prefix = decoder.code.lattice
        (prefix.steps, prefix.Ly, prefix.Lz) == (r, lattice.Ly, lattice.Lz) ||
            throw(DimensionMismatch("H+6 full-prefix decoder has incompatible checkpoint geometry"))
    end
    return nothing
end

"""
    replay_toric3d_h6_full_prefix_x(encoder, faults, family;
        feedback=true, record_history=false)

Replay immutable X-only gate-output faults from the ideal encoded seed, using
`build_toric3d_x_feedback_decoders`' growing-prefix BP+OSD decoders. With
feedback, repair every completed physical plaquette after every transfer
`r=1:T`, including the final transfer. Recovery may act anywhere on `Q₀:Qᵣ`
and never on future slices. In particular, later corrections can revisit the
seed slice, so its logical parity need not freeze after the first transfer.

The common final full-code outcome is called exactly once. With feedback it
receives an already clear syndrome and still scores undetectable logical
errors. With `feedback=false`, this is the unchanged final-only protocol.
Recorded histories use the existing full-prefix X helper's fields.
"""
function replay_toric3d_h6_full_prefix_x(encoder::Toric3DSliceEncoder,
        faults::Toric3DFaultRecord, family::Toric3DXFeedbackDecoderFamily;
        feedback::Bool=true, record_history::Bool=false)::NamedTuple
    _toric3d_validate_h6_full_prefix_family(encoder, family)
    _validate_toric3d_fault_record(encoder, faults)
    any(faults.z) && throw(ArgumentError("H+6 full-prefix X replay requires X-only faults"))
    code = family.code
    frame = PauliFrame(falses(data_qubit_count(code)), falses(data_qubit_count(code)))
    history = NamedTuple[]
    first_operation = 1
    for r in 1:code.lattice.steps
        propagate_toric3d_segment!(frame, encoder, faults,
            first_operation, encoder.checkpoints[r])
        first_operation = encoder.checkpoints[r] + 1
        if feedback
            prefix = toric3d_encoded_prefix(code, r)
            entry = _toric3d_apply_x_prefix_feedback!(
                frame, prefix, family.prefixes[r]; record_history)
            record_history && push!(history, entry)
        end
    end
    outcome = toric3d_x_threshold_outcome(family.prefixes[end], frame.x)
    return (frame=frame, syndrome=outcome.syndrome, recovery=outcome.recovery,
        residual=outcome.residual, logical_1=outcome.logical_1,
        logical_2=outcome.logical_2, failure=outcome.failure, history=history)
end

"""
    replay_toric3d_h6_full_prefix_z(encoder, faults, family;
        feedback=true, record_history=false)

Replay H+6 Z seed and transfer-output faults with the growing-prefix MWPM
decoders from `build_toric3d_z_feedback_decoders`. With feedback, clear every
completed physical X-star check after each transfer `r=1:T`, including the
final transfer. Recovery may act on all prepared slices `Q₀:Qᵣ`. No separate
`r=0` correction is applied: the noisy seed propagates into transfer 1 and is
repaired together with that first completed prefix.

Both feedback and final-only (`feedback=false`) paths invoke the same final
full-code outcome exactly once. The feedback path has already cleared its
syndrome; undetectable logical errors remain in the final score. Recorded
histories use the existing full-prefix Z helper's fields, with `kind=:prefix`.
"""
function replay_toric3d_h6_full_prefix_z(encoder::Toric3DSliceEncoder,
        faults::Toric3DDiagonalZFaultRecord, family::Toric3DZFeedbackDecoderFamily;
        feedback::Bool=true, record_history::Bool=false)::NamedTuple
    _toric3d_validate_h6_full_prefix_family(encoder, family)
    _validate_toric3d_fault_record(encoder, faults.transfer)
    length(faults.seed) == 3encoder.code.Ly * encoder.code.Lz ||
        throw(DimensionMismatch("H+6 full-prefix Z seed and encoder dimensions differ"))
    any(faults.transfer.x) &&
        throw(ArgumentError("H+6 full-prefix Z replay requires Z-only faults"))
    code = family.code
    frame = PauliFrame(falses(data_qubit_count(code)), falses(data_qubit_count(code)))
    frame.z[toric3d_slice_qubits(code.lattice, 0)] .= faults.seed
    history = NamedTuple[]
    first_operation = 1
    for r in 1:code.lattice.steps
        propagate_toric3d_segment!(frame, encoder, faults.transfer,
            first_operation, encoder.checkpoints[r])
        first_operation = encoder.checkpoints[r] + 1
        if feedback
            entry = _toric3d_apply_z_prefix_feedback!(frame, code, r, family; record_history)
            record_history && push!(history, entry)
        end
    end
    outcome = toric3d_z_threshold_outcome(family.prefixes[end], frame.z)
    return (frame=frame, syndrome=outcome.syndrome, recovery=outcome.recovery,
        residual=outcome.residual, logical_1=outcome.logical_1,
        logical_2=outcome.logical_2, failure=outcome.failure, history=history)
end
