"""Seed-H and transfer-output Z faults for one encoded diagonal preparation."""
struct Toric3DDiagonalZFaultRecord
    seed::BitVector
    transfer::Toric3DFaultRecord

    function Toric3DDiagonalZFaultRecord(seed::AbstractVector{Bool},
            transfer::Toric3DFaultRecord)
        expected = 3 * transfer.Ly * transfer.Lz
        length(seed) == expected ||
            throw(DimensionMismatch("seed fault record has the wrong length"))
        any(transfer.x) &&
            throw(ArgumentError("diagonal Z fault records cannot contain X faults"))
        return new(BitVector(seed),transfer)
    end
end

Base.:(==)(left::Toric3DDiagonalZFaultRecord,
           right::Toric3DDiagonalZFaultRecord) =
    left.seed == right.seed && left.transfer == right.transfer

"""Sample seed-H faults first, followed by transfer outputs in schedule order."""
function sample_toric3d_diagonal_z_faults(rng::Random.AbstractRNG,
        encoder::Toric3DSliceEncoder;p::Real)
    isfinite(p) && 0 <= p <= 1 ||
        throw(ArgumentError("p must be finite and lie in [0, 1]"))
    code = encoder.code
    seed = BitVector(rand(rng) < p for _ in 1:3code.Ly*code.Lz)
    x = falses(2,length(encoder.operations))
    z = falses(2,length(encoder.operations))
    for (step,operation) in enumerate(encoder.operations), slot in eachindex(operation.qubits)
        z[slot,step] = rand(rng) < p
    end
    transfer = Toric3DFaultRecord(code.steps,code.Ly,code.Lz,x,z)
    return Toric3DDiagonalZFaultRecord(seed,transfer)
end

function _toric3d_z_family_dimensions(decoders::Toric3DZFeedbackDecoderFamily)
    code = decoders.code.lattice
    return (code.steps,code.Ly,code.Lz)
end

function _toric3d_apply_z_bits!(frame::PauliFrame,qubits,recovery)
    length(qubits) == length(recovery) ||
        throw(DimensionMismatch("Z recovery and mapped support differ in length"))
    for (local_index,full) in enumerate(qubits)
        frame.z[full] ⊻= recovery[local_index]
    end
    return frame
end

function _toric3d_z_prefix_syndrome(prefix::Toric3DEncodedPrefix,
        decoder::Toric3DZThresholdDecoder,frame::PauliFrame)
    length(prefix.qubits) == size(decoder.matrix,2) ||
        throw(DimensionMismatch("Z prefix map and decoder dimensions differ"))
    return BitVector(isodd.(decoder.matrix *
        Vector{UInt8}(frame.z[prefix.qubits])))
end

function _toric3d_z_history(kind::Symbol,checkpoint::Int,syndrome,recovery,
        qubits,post_syndrome)
    mapped = qubits[findall(recovery)]
    return (kind=kind,checkpoint=checkpoint,syndrome=BitVector(syndrome),
        recovery=BitVector(recovery),qubits=copy(qubits),support=mapped,
        post_syndrome=BitVector(post_syndrome))
end

function _toric3d_apply_z_seed_feedback!(frame::PauliFrame,
        decoder::Toric3DZFeedbackDecoder;record_history::Bool)
    window = decoder.window
    syndrome = toric3d_z_window_syndrome(window,frame.z)
    recovery = decode_toric3d_z_feedback(decoder,syndrome)
    _toric3d_apply_z_bits!(frame,window.qubits,recovery)
    post = toric3d_z_window_syndrome(window,frame.z)
    any(post) && error("seed Z feedback left a nonzero seed-star syndrome")
    return record_history ? _toric3d_z_history(
        :seed,0,syndrome,recovery,window.qubits,post) : nothing
end

function _toric3d_apply_z_local_feedback!(frame::PauliFrame,
        code::Toric3DEncodedSlice,r::Int,decoders::Toric3DZFeedbackDecoderFamily;
        record_history::Bool)
    window = toric3d_z_feedback_window(code,r)
    window.matrix == decoders.transfer.window.matrix ||
        error("Z transfer-window ordering differs from the cached decoder")
    syndrome = toric3d_z_window_syndrome(window,frame.z)
    recovery = decode_toric3d_z_feedback(decoders.transfer,syndrome)
    _toric3d_apply_z_bits!(frame,window.qubits,recovery)
    window_post = toric3d_z_window_syndrome(window,frame.z)
    any(window_post) &&
        error("local Z feedback left a nonzero window syndrome")
    prefix = toric3d_encoded_prefix(code,r)
    post = _toric3d_z_prefix_syndrome(prefix,decoders.prefixes[r],frame)
    any(post) &&
        error("local Z feedback violated the complete-prefix star invariant")
    return record_history ? _toric3d_z_history(
        :window,r,syndrome,recovery,window.qubits,post) : nothing
end

function _toric3d_apply_z_prefix_feedback!(frame::PauliFrame,
        code::Toric3DEncodedSlice,r::Int,decoders::Toric3DZFeedbackDecoderFamily;
        record_history::Bool)
    prefix = toric3d_encoded_prefix(code,r)
    decoder = decoders.prefixes[r]
    syndrome = _toric3d_z_prefix_syndrome(prefix,decoder,frame)
    recovery = _decode_toric3d_z(decoder,syndrome)
    _toric3d_apply_z_bits!(frame,prefix.qubits,recovery)
    post = _toric3d_z_prefix_syndrome(prefix,decoder,frame)
    any(post) &&
        error("full-prefix Z feedback left a nonzero star syndrome")
    return record_history ? _toric3d_z_history(
        :prefix,r,syndrome,recovery,prefix.qubits,post) : nothing
end

"""Replay one Z fault record through a selected feedback scope."""
function replay_toric3d_z_strength(encoder::Toric3DSliceEncoder,
        faults::Toric3DDiagonalZFaultRecord,
        decoders::Toric3DZFeedbackDecoderFamily;protocol::Symbol,
        record_history::Bool=false)::NamedTuple
    protocol in (:none,:slice_local,:full_prefix) ||
        throw(ArgumentError("protocol must be :none, :slice_local, or :full_prefix"))
    encoder_dims = (encoder.code.steps,encoder.code.Ly,encoder.code.Lz)
    encoder_dims == _toric3d_z_family_dimensions(decoders) ||
        throw(DimensionMismatch("Z feedback decoders and encoder dimensions differ"))
    length(decoders.prefixes) == encoder.code.steps ||
        throw(DimensionMismatch("Z prefix decoder family has the wrong depth"))
    _validate_toric3d_fault_record(encoder,faults.transfer)
    length(faults.seed) == 3encoder.code.Ly*encoder.code.Lz ||
        throw(DimensionMismatch("seed fault record and encoder dimensions differ"))

    code = decoders.code
    frame = PauliFrame(falses(data_qubit_count(code)),falses(data_qubit_count(code)))
    seed_qubits = toric3d_slice_qubits(code.lattice,0)
    frame.z[seed_qubits] .= faults.seed
    history = NamedTuple[]
    if protocol !== :none
        seed_entry = _toric3d_apply_z_seed_feedback!(
            frame,decoders.seed;record_history)
        record_history && push!(history,seed_entry)
    end

    first_operation = 1
    for r in 1:code.lattice.steps
        propagate_toric3d_segment!(frame,encoder,faults.transfer,
            first_operation,encoder.checkpoints[r])
        first_operation = encoder.checkpoints[r] + 1
        r == code.lattice.steps && continue
        if protocol === :slice_local
            entry = _toric3d_apply_z_local_feedback!(
                frame,code,r,decoders;record_history)
            record_history && push!(history,entry)
        elseif protocol === :full_prefix
            entry = _toric3d_apply_z_prefix_feedback!(
                frame,code,r,decoders;record_history)
            record_history && push!(history,entry)
        end
    end

    outcome = toric3d_z_threshold_outcome(decoders.prefixes[end],frame.z)
    return (frame=frame,syndrome=outcome.syndrome,recovery=outcome.recovery,
        residual=outcome.residual,logical_1=outcome.logical_1,
        logical_2=outcome.logical_2,failure=outcome.failure,history=history)
end

"""Replay one immutable Z fault record through all three protocols."""
function toric3d_z_strength_triplet(encoder::Toric3DSliceEncoder,
        faults::Toric3DDiagonalZFaultRecord,
        decoders::Toric3DZFeedbackDecoderFamily;
        record_history::Bool=false)::NamedTuple{(:none,:slice_local,:full_prefix)}
    none = replay_toric3d_z_strength(
        encoder,faults,decoders;protocol=:none,record_history)
    slice_local = replay_toric3d_z_strength(
        encoder,faults,decoders;protocol=:slice_local,record_history)
    full_prefix = replay_toric3d_z_strength(
        encoder,faults,decoders;protocol=:full_prefix,record_history)
    return (none=none,slice_local=slice_local,full_prefix=full_prefix)
end
