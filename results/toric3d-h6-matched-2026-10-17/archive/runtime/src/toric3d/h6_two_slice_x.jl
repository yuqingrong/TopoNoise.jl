"""Measured physical plaquettes and allowed two-slice X-recovery columns."""
struct Toric3DH6TwoSliceXWindow
    checkpoint::Int
    checks::Vector{Vector{Int}}
    qubits::Vector{Int}
    matrix::BitMatrix
    prefix_checks::Vector{Vector{Int}}
    frontier_checks::Vector{Vector{Int}}
end

"""Cached BP+OSD problem for one first/interior window geometry."""
struct Toric3DH6TwoSliceXKernel
    matrix::BitMatrix
    prior::Float64
    bposd::PythonCall.Py
end

"""Complete X feedback on the last two data slices after every H+6 transfer."""
struct Toric3DH6TwoSliceXDecoder
    code::Toric3DEncodedSlice
    final::Toric3DXThresholdDecoder
    first::Toric3DH6TwoSliceXKernel
    interior::Union{Nothing,Toric3DH6TwoSliceXKernel}
    windows::Vector{Toric3DH6TwoSliceXWindow}
end

function _toric3d_h6_two_slice_x_window(code::Toric3DEncodedSlice, r::Int)
    prefix = toric3d_encoded_prefix(code, r)
    complete = _toric3d_map_prefix_checks(prefix, plaquette_checks(prefix.code))
    cells = code.lattice.Ly * code.lattice.Lz
    first_row = r == 1 ? 1 : 3cells * (r - 2) + 1
    checks = complete[first_row:end]
    qubits = vcat(toric3d_slice_qubits(code.lattice, r - 1),
                  toric3d_slice_qubits(code.lattice, r))
    matrix = _toric3d_restricted_check_matrix(checks, qubits)
    frontier = _toric3d_map_prefix_checks(prefix,
        toric3d_slice_checks(prefix.code.lattice, r))
    return Toric3DH6TwoSliceXWindow(r, checks, qubits, matrix, complete, frontier)
end

function _toric3d_h6_two_slice_x_kernel(matrix::BitMatrix, prior::Float64,
        settings::Toric3DDecoderSettings)
    numpy = PythonCall.pyimport("numpy")
    ldpc = _toric3d_load_ldpc()
    decoder = ldpc.BpOsdDecoder(_toric3d_numpy_matrix(numpy, matrix);
        error_rate=prior, max_iter=settings.max_iter,
        bp_method=settings.bp_method, ms_scaling_factor=settings.ms_scaling_factor,
        schedule=settings.schedule, omp_thread_count=settings.omp_thread_count,
        osd_method=settings.osd_method, osd_order=settings.osd_order)
    return Toric3DH6TwoSliceXKernel(matrix, prior, decoder)
end

"""
    build_toric3d_h6_two_slice_x_decoder(code; p, settings=Toric3DDecoderSettings())

Cache one BP+OSD problem for the first transfer and one for all later
transfers, plus the established final full-code decoder. Recovery columns
belong only to data slices `r-1:r`. Measured checks are the full physical
plaquettes incident to those columns, including the lower-boundary layer
`r-1` for `r>1`; their measurement support reaches data slices `r-2:r`.
"""
function build_toric3d_h6_two_slice_x_decoder(code::Toric3DEncodedSlice; p::Real,
        settings::Toric3DDecoderSettings=Toric3DDecoderSettings())
    isfinite(p) && 0 <= p <= 1 || throw(ArgumentError("p must lie in [0, 1]"))
    final = build_toric3d_x_threshold_decoder(code; p, settings)
    windows = [_toric3d_h6_two_slice_x_window(code, r) for r in 1:code.lattice.steps]
    first = _toric3d_h6_two_slice_x_kernel(windows[1].matrix, final.prior, settings)
    interior = code.lattice.steps == 1 ? nothing :
        _toric3d_h6_two_slice_x_kernel(windows[2].matrix, final.prior, settings)
    if interior !== nothing
        all(window -> window.matrix == interior.matrix, windows[2:end]) ||
            error("H+6 two-slice X interior matrices are not translation equivalent")
    end
    return Toric3DH6TwoSliceXDecoder(code, final, first, interior, windows)
end

function _toric3d_decode_h6_two_slice_x(kernel::Toric3DH6TwoSliceXKernel,
        syndrome::AbstractVector{Bool})
    length(syndrome) == size(kernel.matrix, 1) ||
        throw(DimensionMismatch("two-slice X syndrome has the wrong number of checks"))
    n = size(kernel.matrix, 2)
    !any(syndrome) && return falses(n)
    numpy = PythonCall.pyimport("numpy")
    recovery = _toric3d_py_bits(numpy,
        kernel.bposd.decode(_toric3d_numpy_vector(numpy, syndrome)), n)
    _toric3d_matrix_syndrome(kernel.matrix, recovery) == BitVector(syndrome) ||
        error("two-slice X recovery does not reproduce its measured syndrome")
    return recovery
end

function _toric3d_apply_h6_two_slice_x!(frame::PauliFrame,
        decoder::Toric3DH6TwoSliceXDecoder, r::Int; record_history::Bool)
    1 <= r <= length(decoder.windows) || throw(ArgumentError("X feedback checkpoint is out of range"))
    length(frame.x) == data_qubit_count(decoder.code) ||
        throw(DimensionMismatch("X feedback frame has the wrong length"))
    window = decoder.windows[r]
    kernel = r == 1 ? decoder.first : decoder.interior
    # Only the measured local checks enter BP+OSD, never the hidden fault record,
    # full-prefix audit syndrome, or logical parity.
    measured = _toric3d_slice_parity(window.checks, frame.x)
    recovery = _toric3d_decode_h6_two_slice_x(kernel, measured)
    for (column, q) in enumerate(window.qubits)
        frame.x[q] ⊻= recovery[column]
    end
    post = _toric3d_slice_parity(window.checks, frame.x)
    any(post) && error("two-slice X recovery left a nonzero window syndrome")
    prefix = _toric3d_slice_parity(window.prefix_checks, frame.x)
    any(prefix) && error("two-slice X recovery left a nonzero completed-prefix syndrome")
    frontier = _toric3d_slice_parity(window.frontier_checks, frame.x)
    any(frontier) && error("two-slice X recovery left a nonzero frontier syndrome")
    record_history || return nothing
    return (checkpoint=r, checks=deepcopy(window.checks),
        measurement_support=sort!(unique(vcat(window.checks...))),
        qubits=copy(window.qubits), recovery=recovery,
        support=window.qubits[findall(recovery)], pre_syndrome=measured,
        post_syndrome=post, prefix_syndrome=prefix, frontier_syndrome=frontier,
        logical_after=_toric3d_matrix_syndrome(decoder.final.logicals, frame.x))
end

"""
    replay_toric3d_h6_two_slice_x(encoder, faults, decoder;
        feedback=true, record_history=false)

Replay immutable X-only gate-output faults from the ideal encoded seed.
With feedback, completely repair the newest two data slices after EVERY
transfer, including the last. Audit all completed-prefix checks before
retirement. With feedback disabled, propagate the unchanged final-only
baseline. Both paths call the same final full-code outcome exactly once.
"""
function replay_toric3d_h6_two_slice_x(encoder::Toric3DSliceEncoder,
        faults::Toric3DFaultRecord, decoder::Toric3DH6TwoSliceXDecoder;
        feedback::Bool=true, record_history::Bool=false)::NamedTuple
    code = decoder.code
    (encoder.code.steps, encoder.code.Ly, encoder.code.Lz) ==
        (code.lattice.steps, code.lattice.Ly, code.lattice.Lz) ||
        throw(DimensionMismatch("two-slice X decoder and encoder dimensions differ"))
    _validate_toric3d_fault_record(encoder, faults)
    any(faults.z) && throw(ArgumentError("two-slice X replay requires X-only faults"))
    frame = PauliFrame(falses(data_qubit_count(code)), falses(data_qubit_count(code)))
    history = NamedTuple[]
    first_operation = 1
    for r in 1:code.lattice.steps
        propagate_toric3d_segment!(frame, encoder, faults,
            first_operation, encoder.checkpoints[r])
        first_operation = encoder.checkpoints[r] + 1
        if feedback
            entry = _toric3d_apply_h6_two_slice_x!(frame, decoder, r; record_history)
            record_history && push!(history, entry)
        end
    end
    outcome = toric3d_x_threshold_outcome(decoder.final, frame.x)
    return (frame=frame, syndrome=outcome.syndrome, recovery=outcome.recovery,
        residual=outcome.residual, logical_1=outcome.logical_1,
        logical_2=outcome.logical_2, failure=outcome.failure, history=history)
end
