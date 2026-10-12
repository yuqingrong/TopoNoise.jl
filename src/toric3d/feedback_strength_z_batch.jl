"""Columns-as-shots seed and transfer Z faults for batched adaptive replay."""
struct Toric3DDiagonalZFaultBatch
    steps::Int
    Ly::Int
    Lz::Int
    seed::BitMatrix
    transfer::BitArray{3}

    function Toric3DDiagonalZFaultBatch(steps::Integer,Ly::Integer,Lz::Integer,
            seed::AbstractMatrix{Bool},transfer::AbstractArray{Bool,3})
        code = Toric3DSliceLattice(steps,Ly,Lz)
        size(seed,1) == 3code.Ly*code.Lz ||
            throw(DimensionMismatch("seed fault batch has the wrong row count"))
        size(transfer,1) == 2 ||
            throw(DimensionMismatch("transfer fault batch must have two output slots"))
        size(seed,2) == size(transfer,3) ||
            throw(DimensionMismatch("seed and transfer fault batches differ in shots"))
        size(seed,2) > 0 || throw(ArgumentError("fault batch must contain at least one shot"))
        return new(code.steps,code.Ly,code.Lz,BitMatrix(seed),BitArray(transfer))
    end
end

Base.:(==)(left::Toric3DDiagonalZFaultBatch,
           right::Toric3DDiagonalZFaultBatch) =
    (left.steps,left.Ly,left.Lz) == (right.steps,right.Ly,right.Lz) &&
    left.seed == right.seed && left.transfer == right.transfer

"""Sample shots in scalar order: seed locations, then operation/output slots."""
function sample_toric3d_diagonal_z_fault_batch(rng::Random.AbstractRNG,
        encoder::Toric3DSliceEncoder;p::Real,shots::Integer)
    isfinite(p) && 0 <= p <= 1 ||
        throw(ArgumentError("p must be finite and lie in [0, 1]"))
    count_shots = _positive_machine_int(shots,"shots")
    code = encoder.code
    seed = falses(3code.Ly*code.Lz,count_shots)
    transfer = falses(2,length(encoder.operations),count_shots)
    for shot in 1:count_shots
        for q in axes(seed,1)
            seed[q,shot] = rand(rng) < p
        end
        for (step,operation) in enumerate(encoder.operations), slot in eachindex(operation.qubits)
            transfer[slot,step,shot] = rand(rng) < p
        end
    end
    return Toric3DDiagonalZFaultBatch(
        code.steps,code.Ly,code.Lz,seed,transfer)
end

function _toric3d_z_batch_syndrome(matrix,frames::BitMatrix,qubits)
    return BitMatrix(isodd.(matrix * Matrix{UInt8}(frames[qubits,:])))
end

function _toric3d_z_batch_window_syndrome(window::Toric3DZFeedbackWindow,
        frames::BitMatrix)
    size(frames,1) == data_qubit_count(window.code) ||
        throw(DimensionMismatch("Z frame batch has the wrong row count"))
    syndromes = falses(length(window.checks),size(frames,2))
    for (row,check) in enumerate(window.checks), q in check
        @views syndromes[row,:] .= xor.(syndromes[row,:],frames[q,:])
    end
    return syndromes
end

function _toric3d_apply_z_batch!(frames::BitMatrix,qubits,recoveries::BitMatrix)
    size(recoveries) == (length(qubits),size(frames,2)) ||
        throw(DimensionMismatch("mapped Z recovery batch has the wrong dimensions"))
    for (local_index,full) in enumerate(qubits)
        @views frames[full,:] .= xor.(frames[full,:],recoveries[local_index,:])
    end
    return frames
end

function _toric3d_propagate_z_batch!(frames::BitMatrix,
        encoder::Toric3DSliceEncoder,faults::Toric3DDiagonalZFaultBatch,
        first_operation::Int,last_operation::Int)
    for step in first_operation:last_operation
        operation = encoder.operations[step]
        if operation.gate === :H
            q = only(operation.qubits)
            any(@view frames[q,:]) &&
                error("pure-Z batch reached a nonfresh Hadamard target")
        elseif operation.gate === :CNOT
            control,target = operation.qubits
            @views frames[control,:] .= xor.(frames[control,:],frames[target,:])
        else
            error("unsupported gate in diagonal Z batch propagation")
        end
        for (slot,q) in enumerate(operation.qubits)
            @views frames[q,:] .= xor.(frames[q,:],faults.transfer[slot,step,:])
        end
    end
    return frames
end

function _toric3d_z_batch_final_outcome(decoder::Toric3DZThresholdDecoder,
        frames::BitMatrix)
    syndrome = BitMatrix(isodd.(decoder.matrix*Matrix{UInt8}(frames)))
    recovery = decode_toric3d_z_batch(decoder,syndrome)
    residual = frames .⊻ recovery
    logical = BitMatrix(isodd.(decoder.logicals*Matrix{UInt8}(residual)))
    logical_1 = BitVector(logical[1,:])
    logical_2 = BitVector(logical[2,:])
    return (logical_1=logical_1,logical_2=logical_2,
        failure=logical_1 .| logical_2)
end

"""Vectorized three-protocol Z replay with validated checkpoint recoveries."""
function toric3d_z_strength_batch(encoder::Toric3DSliceEncoder,
        faults::Toric3DDiagonalZFaultBatch,
        decoders::Toric3DZFeedbackDecoderFamily)::NamedTuple{
            (:none,:slice_local,:full_prefix)}
    encoder_dims = (encoder.code.steps,encoder.code.Ly,encoder.code.Lz)
    encoder_dims == _toric3d_z_family_dimensions(decoders) ||
        throw(DimensionMismatch("Z feedback decoders and encoder dimensions differ"))
    encoder_dims == (faults.steps,faults.Ly,faults.Lz) ||
        throw(DimensionMismatch("Z fault batch and encoder dimensions differ"))
    size(faults.transfer,2) == length(encoder.operations) ||
        throw(DimensionMismatch("Z fault batch has the wrong operation count"))
    length(decoders.prefixes) == encoder.code.steps ||
        throw(DimensionMismatch("Z prefix decoder family has the wrong depth"))

    code = decoders.code
    shots = size(faults.seed,2)
    n = data_qubit_count(code)
    none = falses(n,shots)
    slice_local = falses(n,shots)
    full_prefix = falses(n,shots)
    seed_qubits = toric3d_slice_qubits(code.lattice,0)
    none[seed_qubits,:] .= faults.seed
    slice_local[seed_qubits,:] .= faults.seed
    full_prefix[seed_qubits,:] .= faults.seed

    seed_syndrome = _toric3d_z_batch_window_syndrome(
        decoders.seed.window,slice_local)
    seed_recovery = decode_toric3d_z_feedback_batch(decoders.seed,seed_syndrome)
    _toric3d_apply_z_batch!(slice_local,decoders.seed.window.qubits,seed_recovery)
    _toric3d_apply_z_batch!(full_prefix,decoders.seed.window.qubits,seed_recovery)
    any(_toric3d_z_batch_window_syndrome(
        decoders.seed.window,slice_local)) &&
        error("batched seed Z feedback left a nonzero syndrome")

    first_operation = 1
    for r in 1:code.lattice.steps
        last_operation = encoder.checkpoints[r]
        _toric3d_propagate_z_batch!(none,encoder,faults,first_operation,last_operation)
        _toric3d_propagate_z_batch!(slice_local,encoder,faults,first_operation,last_operation)
        _toric3d_propagate_z_batch!(full_prefix,encoder,faults,first_operation,last_operation)
        first_operation = last_operation + 1
        r == code.lattice.steps && continue

        window = toric3d_z_feedback_window(code,r)
        window.matrix == decoders.transfer.window.matrix ||
            error("Z transfer-window ordering differs from the cached batch decoder")
        local_syndrome = _toric3d_z_batch_window_syndrome(
            window,slice_local)
        local_recovery = decode_toric3d_z_feedback_batch(
            decoders.transfer,local_syndrome)
        _toric3d_apply_z_batch!(slice_local,window.qubits,local_recovery)

        prefix = toric3d_encoded_prefix(code,r)
        prefix_decoder = decoders.prefixes[r]
        local_post = _toric3d_z_batch_syndrome(
            prefix_decoder.matrix,slice_local,prefix.qubits)
        any(local_post) &&
            error("batched local Z feedback violated the prefix invariant")

        prefix_syndrome = _toric3d_z_batch_syndrome(
            prefix_decoder.matrix,full_prefix,prefix.qubits)
        prefix_recovery = decode_toric3d_z_batch(prefix_decoder,prefix_syndrome)
        _toric3d_apply_z_batch!(full_prefix,prefix.qubits,prefix_recovery)
        prefix_post = _toric3d_z_batch_syndrome(
            prefix_decoder.matrix,full_prefix,prefix.qubits)
        any(prefix_post) &&
            error("batched full-prefix Z feedback left a nonzero syndrome")
    end

    final_decoder = decoders.prefixes[end]
    return (none=_toric3d_z_batch_final_outcome(final_decoder,none),
        slice_local=_toric3d_z_batch_final_outcome(final_decoder,slice_local),
        full_prefix=_toric3d_z_batch_final_outcome(final_decoder,full_prefix))
end
