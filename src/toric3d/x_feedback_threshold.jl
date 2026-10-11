"""Final BP+OSD X recovery and logical scoring for the encoded diagonal slab."""
struct Toric3DXThresholdDecoder
    code::Toric3DEncodedSlice
    p::Float64
    prior::Float64
    settings::Toric3DDecoderSettings
    matrix::BitMatrix
    logicals::BitMatrix
    bposd::PythonCall.Py
    version::String
end

"""Build the fixed-baseline full-plaquette BP+OSD decoder."""
function build_toric3d_x_threshold_decoder(code::Toric3DEncodedSlice; p::Real,
        settings::Toric3DDecoderSettings=Toric3DDecoderSettings())
    isfinite(p) && 0 <= p <= 1 || throw(ArgumentError("p must lie in [0, 1]"))
    matrix = BitMatrix(plaquette_check_matrix(code))
    logicals = BitMatrix(_check_matrix(toric3d_encoded_logical_z(code),
        data_qubit_count(code); sparse=false))
    prior = clamp(Float64(p), 1e-12, 1 - 1e-12)
    numpy = PythonCall.pyimport("numpy")
    ldpc = _toric3d_load_ldpc()
    bposd = ldpc.BpOsdDecoder(_toric3d_numpy_matrix(numpy, matrix);
        error_rate=prior, max_iter=settings.max_iter,
        bp_method=settings.bp_method,
        ms_scaling_factor=settings.ms_scaling_factor,
        schedule=settings.schedule,
        omp_thread_count=settings.omp_thread_count,
        osd_method=settings.osd_method,
        osd_order=settings.osd_order)
    return Toric3DXThresholdDecoder(code, Float64(p), prior, settings, matrix,
        logicals, bposd, _toric3d_module_version(ldpc, "ldpc"))
end

function _toric3d_validate_x_recovery(decoder::Toric3DXThresholdDecoder,
        syndrome::AbstractVector{Bool}, recovery::AbstractVector{Bool})
    length(syndrome) == size(decoder.matrix, 1) ||
        throw(DimensionMismatch("plaquette syndrome has the wrong length"))
    length(recovery) == size(decoder.matrix, 2) ||
        throw(DimensionMismatch("X recovery has the wrong length"))
    checked = BitVector(recovery)
    _toric3d_matrix_syndrome(decoder.matrix, checked) == BitVector(syndrome) ||
        throw(ErrorException("invalid X recovery: its plaquette syndrome differs from the measurement"))
    return checked
end

"""Decode one measured full plaquette syndrome into a validated data-X frame."""
function decode_toric3d_x(decoder::Toric3DXThresholdDecoder,
        plaquettes::AbstractVector{Bool})::BitVector
    length(plaquettes) == size(decoder.matrix, 1) ||
        throw(DimensionMismatch("plaquette syndrome has the wrong length"))
    n = size(decoder.matrix, 2)
    !any(plaquettes) && return falses(n)
    numpy = PythonCall.pyimport("numpy")
    recovery = _toric3d_py_bits(numpy,
        decoder.bposd.decode(_toric3d_numpy_vector(numpy, plaquettes)), n)
    return _toric3d_validate_x_recovery(decoder, plaquettes, recovery)
end

"""Apply final recovery and score both residual logical-X membrane parities."""
function toric3d_x_threshold_outcome(decoder::Toric3DXThresholdDecoder,
        frame::AbstractVector{Bool})::NamedTuple{
            (:syndrome,:recovery,:residual,:logical_1,:logical_2,:failure),
            Tuple{BitVector,BitVector,BitVector,Bool,Bool,Bool}}
    length(frame) == size(decoder.matrix, 2) ||
        throw(DimensionMismatch("X frame has the wrong length"))
    data = BitVector(frame)
    syndrome = _toric3d_matrix_syndrome(decoder.matrix, data)
    recovery = decode_toric3d_x(decoder, syndrome)
    residual = data .⊻ recovery
    logical = _toric3d_matrix_syndrome(decoder.logicals, residual)
    logical_1, logical_2 = logical
    return (syndrome=syndrome, recovery=recovery, residual=residual,
        logical_1=logical_1, logical_2=logical_2,
        failure=logical_1 | logical_2)
end

"""Replay one X-only preparation fault record and apply the common final decoder."""
function replay_toric3d_x_feedback(encoder::Toric3DSliceEncoder,
        faults::Toric3DFaultRecord, slice_decoder::Toric3DSliceDecoder,
        final_decoder::Toric3DXThresholdDecoder; feedback::Bool,
        routing::Symbol=:one_step, record_history::Bool=false)::NamedTuple
    encoded_lattice = final_decoder.code.lattice
    lattice = encoder.code
    (lattice.steps, lattice.Ly, lattice.Lz) ==
        (encoded_lattice.steps, encoded_lattice.Ly, encoded_lattice.Lz) ||
        throw(DimensionMismatch("final decoder and encoder dimensions differ"))
    any(faults.z) && throw(ArgumentError("X-feedback replay requires an X-only fault record"))
    prepared = replay_toric3d_slice(encoder, faults, slice_decoder;
        feedback, routing, record_history)
    outcome = toric3d_x_threshold_outcome(final_decoder, prepared.frame.x)
    return (frame=prepared.frame, syndrome=outcome.syndrome,
        recovery=outcome.recovery, residual=outcome.residual,
        logical_1=outcome.logical_1, logical_2=outcome.logical_2,
        failure=outcome.failure, history=prepared.history)
end

"""Replay identical sampled faults under final-only and slice-feedback protocols."""
function toric3d_x_feedback_pair(encoder::Toric3DSliceEncoder,
        faults::Toric3DFaultRecord, slice_decoder::Toric3DSliceDecoder,
        final_decoder::Toric3DXThresholdDecoder; routing::Symbol=:one_step,
        record_history::Bool=false)::NamedTuple{(:final_only,:slice_feedback)}
    final_only = replay_toric3d_x_feedback(encoder, faults, slice_decoder,
        final_decoder; feedback=false, routing, record_history)
    slice_feedback = replay_toric3d_x_feedback(encoder, faults, slice_decoder,
        final_decoder; feedback=true, routing, record_history)
    return (final_only=final_only, slice_feedback=slice_feedback)
end
