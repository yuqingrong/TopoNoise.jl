"""Cached active-slice and growing-prefix X decoders for three-way replay."""
struct Toric3DXFeedbackDecoderFamily
    code::Toric3DEncodedSlice
    slice::Toric3DSliceDecoder
    prefixes::Vector{Toric3DXThresholdDecoder}
end

function build_toric3d_x_feedback_decoders(code::Toric3DEncodedSlice; p::Real,
        settings::Toric3DDecoderSettings=Toric3DDecoderSettings())
    isfinite(p) && 0 <= p <= 1 || throw(ArgumentError("p must lie in [0, 1]"))
    slice = build_toric3d_slice_decoder(code.lattice)
    prefixes = [build_toric3d_x_threshold_decoder(
        toric3d_encoded_prefix(code,t).code; p, settings)
        for t in 1:code.lattice.steps]
    return Toric3DXFeedbackDecoderFamily(code,slice,prefixes)
end

function _toric3d_x_family_dimensions(decoders::Toric3DXFeedbackDecoderFamily)
    code = decoders.code.lattice
    return (code.steps,code.Ly,code.Lz)
end

"""Apply and audit one full-prefix X recovery without inspecting fault bits."""
function _toric3d_apply_x_prefix_feedback!(frame::PauliFrame,
        prefix::Toric3DEncodedPrefix, decoder::Toric3DXThresholdDecoder;
        record_history::Bool=true)
    prefix_dims = (prefix.code.lattice.steps,prefix.code.lattice.Ly,
        prefix.code.lattice.Lz)
    decoder_dims = (decoder.code.lattice.steps,decoder.code.lattice.Ly,
        decoder.code.lattice.Lz)
    prefix_dims == decoder_dims ||
        throw(DimensionMismatch("X prefix map and decoder dimensions differ"))
    maximum(prefix.qubits) <= length(frame.x) ||
        throw(DimensionMismatch("X prefix map exceeds the Pauli frame"))

    local_frame = BitVector(frame.x[prefix.qubits])
    pre_syndrome = _toric3d_matrix_syndrome(decoder.matrix,local_frame)
    recovery = decode_toric3d_x(decoder,pre_syndrome)
    future_qubits = record_history ? setdiff(collect(eachindex(frame.x)),prefix.qubits) : Int[]
    future_before = record_history ? BitVector(frame.x[future_qubits]) : BitVector()
    for (local_index,full) in enumerate(prefix.qubits)
        frame.x[full] ⊻= recovery[local_index]
    end
    corrected = BitVector(frame.x[prefix.qubits])
    post_syndrome = _toric3d_matrix_syndrome(decoder.matrix,corrected)
    any(post_syndrome) &&
        error("full-prefix X recovery left a nonzero plaquette syndrome")
    frontier_checks = toric3d_slice_checks(prefix.code.lattice,
        prefix.code.lattice.steps)
    frontier_syndrome = _toric3d_slice_parity(frontier_checks,corrected)
    any(frontier_syndrome) &&
        error("full-prefix X recovery left a nonzero frontier syndrome")
    record_history || return nothing
    future_after = BitVector(frame.x[future_qubits])
    return (checkpoint=prefix.code.lattice.steps,
        pre_syndrome=pre_syndrome,recovery=recovery,qubits=copy(prefix.qubits),
        post_syndrome=post_syndrome,frontier_syndrome=frontier_syndrome,
        future_qubits=future_qubits,future_before=future_before,
        future_after=future_after)
end

function _replay_toric3d_x_full_prefix(encoder::Toric3DSliceEncoder,
        faults::Toric3DFaultRecord,decoders::Toric3DXFeedbackDecoderFamily;
        routing::Symbol,record_history::Bool)
    code = decoders.code
    # With no intermediate checkpoint, full-prefix and final-only replay are
    # exactly the same operation, including the established history format.
    code.lattice.steps == 1 && return replay_toric3d_x_feedback(
        encoder,faults,decoders.slice,only(decoders.prefixes);
        feedback=false,routing,record_history)

    frame = PauliFrame(falses(data_qubit_count(code)),falses(data_qubit_count(code)))
    history = NamedTuple[]
    first_operation = 1
    for t in 1:code.lattice.steps
        propagate_toric3d_segment!(frame,encoder,faults,first_operation,
            encoder.checkpoints[t])
        first_operation = encoder.checkpoints[t] + 1
        if t < code.lattice.steps
            prefix = toric3d_encoded_prefix(code,t)
            checkpoint = _toric3d_apply_x_prefix_feedback!(
                frame,prefix,decoders.prefixes[t];record_history)
            record_history && push!(history,checkpoint)
        end
    end
    outcome = toric3d_x_threshold_outcome(decoders.prefixes[end],frame.x)
    return (frame=frame,syndrome=outcome.syndrome,recovery=outcome.recovery,
        residual=outcome.residual,logical_1=outcome.logical_1,
        logical_2=outcome.logical_2,failure=outcome.failure,history=history)
end

"""Replay one of `:none`, `:slice_local`, or `:full_prefix` for X faults."""
function replay_toric3d_x_strength(encoder::Toric3DSliceEncoder,
        faults::Toric3DFaultRecord,decoders::Toric3DXFeedbackDecoderFamily;
        protocol::Symbol,routing::Symbol=:one_step,
        record_history::Bool=false)::NamedTuple
    protocol in (:none,:slice_local,:full_prefix) ||
        throw(ArgumentError("protocol must be :none, :slice_local, or :full_prefix"))
    routing in (:one_step,:full_path) ||
        throw(ArgumentError("routing must be :one_step or :full_path"))
    encoder_dims = (encoder.code.steps,encoder.code.Ly,encoder.code.Lz)
    encoder_dims == _toric3d_x_family_dimensions(decoders) ||
        throw(DimensionMismatch("X feedback decoders and encoder dimensions differ"))
    slice_dims = (decoders.slice.code.steps,decoders.slice.code.Ly,
        decoders.slice.code.Lz)
    slice_dims == encoder_dims ||
        throw(DimensionMismatch("active-slice decoder and encoder dimensions differ"))
    length(decoders.prefixes) == encoder.code.steps ||
        throw(DimensionMismatch("X prefix decoder family has the wrong depth"))
    _validate_toric3d_fault_record(encoder,faults)
    any(faults.z) &&
        throw(ArgumentError("X-strength replay requires an X-only fault record"))

    if protocol === :none
        return replay_toric3d_x_feedback(encoder,faults,decoders.slice,
            decoders.prefixes[end];feedback=false,routing,record_history)
    elseif protocol === :slice_local
        return replay_toric3d_x_feedback(encoder,faults,decoders.slice,
            decoders.prefixes[end];feedback=true,routing,record_history)
    end
    return _replay_toric3d_x_full_prefix(encoder,faults,decoders;
        routing,record_history)
end

"""Replay one immutable X fault record through all three protocols."""
function toric3d_x_strength_triplet(encoder::Toric3DSliceEncoder,
        faults::Toric3DFaultRecord,decoders::Toric3DXFeedbackDecoderFamily;
        routing::Symbol=:one_step,record_history::Bool=false)::NamedTuple{
            (:none,:slice_local,:full_prefix)}
    none = replay_toric3d_x_strength(encoder,faults,decoders;
        protocol=:none,routing,record_history)
    slice_local = replay_toric3d_x_strength(encoder,faults,decoders;
        protocol=:slice_local,routing,record_history)
    full_prefix = replay_toric3d_x_strength(encoder,faults,decoders;
        protocol=:full_prefix,routing,record_history)
    return (none=none,slice_local=slice_local,full_prefix=full_prefix)
end
