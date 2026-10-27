"""
    replay_toric3d_h6_two_slice_z(encoder, faults, family;
        feedback=true, record_history=false)

Replay the unchanged H+6 circuit with complete Z recovery on data slices
`r-1,r` after every transfer `r=1:T`, including the last transfer. No separate
seed correction is applied: the first window clears noisy seed and first
transfer syndromes together. Each window measures the existing three X-star
sheets. The complete prepared-prefix syndrome is checked after each correction.

All arms then invoke the same final full-code decoder and logical score
exactly once. With `feedback=false`, all window corrections are disabled,
reproducing the existing final-only protocol. With feedback,
the final window has already cleared the syndrome; the final score still
retains any undetectable logical error.
"""
function replay_toric3d_h6_two_slice_z(encoder::Toric3DSliceEncoder,
        faults::Toric3DDiagonalZFaultRecord, family::Toric3DZFeedbackDecoderFamily;
        feedback::Bool=true, record_history::Bool=false)::NamedTuple
    encoder_dims = (encoder.code.steps, encoder.code.Ly, encoder.code.Lz)
    encoder_dims == _toric3d_z_family_dimensions(family) ||
        throw(DimensionMismatch("Z two-slice decoder and encoder dimensions differ"))
    length(family.prefixes) == encoder.code.steps ||
        throw(DimensionMismatch("Z two-slice prefix decoder family has the wrong depth"))
    _validate_toric3d_fault_record(encoder, faults.transfer)
    length(faults.seed) == 3encoder.code.Ly*encoder.code.Lz ||
        throw(DimensionMismatch("Z two-slice seed record and encoder dimensions differ"))

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
            entry = _toric3d_apply_z_local_feedback!(frame, code, r, family; record_history)
            record_history && push!(history, entry)
        end
    end
    outcome = toric3d_z_threshold_outcome(family.prefixes[end], frame.z)
    return (frame=frame, syndrome=outcome.syndrome, recovery=outcome.recovery,
        residual=outcome.residual, logical_1=outcome.logical_1,
        logical_2=outcome.logical_2, failure=outcome.failure, history=history)
end
