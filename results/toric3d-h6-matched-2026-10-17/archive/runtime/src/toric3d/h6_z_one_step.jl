"""Cached closed three-sheet physical graph for causal H+6 Z feedback."""
struct Toric3DH6ZOneStepDecoder
    family::Toric3DZFeedbackDecoderFamily
    adjacency::Vector{Vector{Tuple{Int,Int}}}
    checks::Vector{Vector{Int}}
    windows::Vector{Toric3DZFeedbackWindow}
    prefixes::Vector{Toric3DEncodedPrefix}
end

"""Reuse the complete Z decoder family and retain every physical edge ID."""
function build_toric3d_h6_z_one_step_decoder(family::Toric3DZFeedbackDecoderFamily)
    code = family.code
    cells = code.lattice.Ly * code.lattice.Lz
    matrix = family.transfer.window.matrix
    size(matrix) == (3cells, 6cells) ||
        throw(DimensionMismatch("H+6 Z feedback requires a three-sheet transfer graph"))
    adjacency = [Tuple{Int,Int}[] for _ in axes(matrix, 1)]
    for q in axes(matrix, 2)
        ends = findall(@view matrix[:, q])
        length(ends) == 2 || error("H+6 Z graph edges must have two physical endpoints")
        a, b = ends
        if q <= 3cells
            a <= cells && cells < b <= 2cells ||
                error("retiring data must connect retiring checks to middle checks")
        else
            cells < a <= 2cells && 2cells < b <= 3cells ||
                error("fresh data must connect middle checks to newest checks")
        end
        push!(adjacency[a], (b, q))
        push!(adjacency[b], (a, q))
    end
    foreach(sort!, adjacency)
    windows = [toric3d_z_feedback_window(code, r) for r in 1:code.lattice.steps]
    all(window -> window.matrix == matrix, windows) ||
        error("H+6 Z transfer windows differ from the cached physical graph")
    prefixes = [toric3d_encoded_prefix(code, r) for r in 1:code.lattice.steps]
    return Toric3DH6ZOneStepDecoder(family, adjacency,
        [findall(row) for row in eachrow(matrix)], windows, prefixes)
end

"""
    decode_toric3d_h6_z_one_step(decoder, syndrome)

Match all measured defects on the closed three-sheet H+6 Z graph, then lift
each matched pair to a deterministic physical shortest path. `recovery`
contains the full paths. `route` applies only the endpoint edge of a matched
defect in the first (retiring) sheet, moving it into the middle sheet.
Middle and newest endpoints are left in place. The returned residual can
therefore be nonzero, but its retiring sheet is always zero.
"""
function decode_toric3d_h6_z_one_step(decoder::Toric3DH6ZOneStepDecoder,
        syndrome::AbstractVector{Bool})
    transfer = decoder.family.transfer
    rows, columns = size(transfer.matrix)
    length(syndrome) == rows ||
        throw(DimensionMismatch("H+6 Z syndrome has the wrong number of checks"))
    iseven(count(syndrome)) ||
        throw(ArgumentError("closed H+6 Z graph syndromes must have even parity"))
    cells = decoder.family.code.lattice.Ly * decoder.family.code.lattice.Lz
    full = falses(columns)
    route = falses(columns)
    pairs = Tuple{Int,Int}[]
    paths = Vector{Int}[]
    endpoints = falses(rows)
    if any(syndrome)
        numpy = PythonCall.pyimport("numpy")
        matched = transfer.matcher.decode_to_matched_dets_array(
            _toric3d_numpy_vector(numpy, syndrome))
        raw = PythonCall.pyconvert(Vector{Vector{Int}}, matched.tolist())
        for pair in raw
            length(pair) == 2 || error("H+6 Z matching returned an invalid pair")
            a, b = pair .+ 1
            1 <= a <= rows && 1 <= b <= rows && a != b ||
                error("H+6 Z matching returned an invalid detector or boundary")
            endpoints[a] ⊻= true
            endpoints[b] ⊻= true
            path = _toric3d_slice_shortest_path(decoder.adjacency, a, b)
            push!(pairs, (a, b))
            push!(paths, path)
            for q in path
                full[q] ⊻= true
            end
            a <= cells && (route[first(path)] ⊻= true)
            b <= cells && (route[last(path)] ⊻= true)
        end
    end
    endpoints == syndrome || error("H+6 Z matching omitted measured defects")
    _toric3d_validate_z_feedback_recovery(transfer, syndrome, full)
    any(@view route[3cells+1:end]) &&
        error("H+6 Z one-step route left the retiring data slice")
    residual = syndrome .⊻ _toric3d_slice_parity(decoder.checks, route)
    any(@view residual[1:cells]) &&
        error("H+6 Z one-step route left a retiring defect")
    return (recovery=full, route=route, matched_pairs=pairs,
        paths=paths, residual_syndrome=residual)
end

function _toric3d_apply_h6_z_one_step!(frame::PauliFrame,
        decoder::Toric3DH6ZOneStepDecoder, r::Int; record_history::Bool)
    family = decoder.family
    cells = family.code.lattice.Ly * family.code.lattice.Lz
    window = decoder.windows[r]
    # Only this measured three-sheet window is supplied to the decoder.
    syndrome = toric3d_z_window_syndrome(window, frame.z)
    decision = decode_toric3d_h6_z_one_step(decoder, syndrome)
    prefix = decoder.prefixes[r]
    prefix_before = record_history ?
        _toric3d_z_prefix_syndrome(prefix, family.prefixes[r], frame) : nothing
    _toric3d_apply_z_bits!(frame, window.qubits, decision.route)
    window_post = toric3d_z_window_syndrome(window, frame.z)
    window_post == decision.residual_syndrome ||
        error("H+6 Z observed window syndrome differs from the partial-route prediction")
    # Complete-prefix checks are an audit, never a decoder input. Its newest
    # two sheets may retain defects; every earlier sheet must be clear.
    post = _toric3d_z_prefix_syndrome(prefix, family.prefixes[r], frame)
    retired_rows = 1:r*cells
    any(@view post[retired_rows]) &&
        error("H+6 Z one-step feedback violated the retired-prefix invariant")
    record_history || return nothing
    return (kind=:window, checkpoint=r, syndrome=BitVector(syndrome),
        recovery=decision.recovery, route=decision.route,
        matched_pairs=decision.matched_pairs, paths=decision.paths,
        residual_syndrome=decision.residual_syndrome,
        qubits=copy(window.qubits), support=window.qubits[findall(decision.route)],
        full_support=window.qubits[findall(decision.recovery)],
        retiring_syndrome=syndrome[1:cells],
        predicted_retiring_syndrome=decision.residual_syndrome[1:cells],
        observed_retiring_syndrome=window_post[1:cells],
        window_post_syndrome=window_post, prefix_syndrome=prefix_before,
        post_syndrome=post, retired_prefix_rows=retired_rows)
end

"""
    replay_toric3d_h6_z_one_step(encoder, faults, decoder;
        intermediate_feedback=true, record_history=false)

Replay the unchanged H+6 circuit with a complete noisy-seed correction,
causal one-step Z corrections after transfers `1:T-1`, and exactly one
complete final MWPM recovery and logical score. With
`intermediate_feedback=false`, the same seed correction is retained and
intermediate corrections are skipped, giving the seed-only control.
"""
function replay_toric3d_h6_z_one_step(encoder::Toric3DSliceEncoder,
        faults::Toric3DDiagonalZFaultRecord, decoder::Toric3DH6ZOneStepDecoder;
        intermediate_feedback::Bool=true, record_history::Bool=false)::NamedTuple
    family = decoder.family
    encoder_dims = (encoder.code.steps, encoder.code.Ly, encoder.code.Lz)
    encoder_dims == _toric3d_z_family_dimensions(family) ||
        throw(DimensionMismatch("H+6 Z decoder and encoder dimensions differ"))
    length(family.prefixes) == encoder.code.steps ||
        throw(DimensionMismatch("H+6 Z prefix decoder family has the wrong depth"))
    _validate_toric3d_fault_record(encoder, faults.transfer)
    length(faults.seed) == 3encoder.code.Ly*encoder.code.Lz ||
        throw(DimensionMismatch("H+6 Z seed record and encoder dimensions differ"))
    code = family.code
    frame = PauliFrame(falses(data_qubit_count(code)), falses(data_qubit_count(code)))
    frame.z[toric3d_slice_qubits(code.lattice, 0)] .= faults.seed
    history = NamedTuple[]
    seed_entry = _toric3d_apply_z_seed_feedback!(frame, family.seed; record_history)
    record_history && push!(history, seed_entry)
    first_operation = 1
    for r in 1:code.lattice.steps
        propagate_toric3d_segment!(frame, encoder, faults.transfer,
            first_operation, encoder.checkpoints[r])
        first_operation = encoder.checkpoints[r] + 1
        if intermediate_feedback && r < code.lattice.steps
            entry = _toric3d_apply_h6_z_one_step!(frame, decoder, r; record_history)
            record_history && push!(history, entry)
        end
    end
    outcome = toric3d_z_threshold_outcome(family.prefixes[end], frame.z)
    return (frame=frame, syndrome=outcome.syndrome, recovery=outcome.recovery,
        residual=outcome.residual, logical_1=outcome.logical_1,
        logical_2=outcome.logical_2, failure=outcome.failure, history=history)
end
