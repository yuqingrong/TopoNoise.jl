"""Unit-weight matching of active-slice six-edge checks, with physical edge IDs."""
struct Toric3DSliceDecoder
    code::Toric3DSliceLattice
    matrix::BitMatrix
    checks::Vector{Vector{Int}}
    adjacency::Vector{Vector{Tuple{Int,Int}}}
    matcher::PythonCall.Py
    version::String
    stars::Vector{Vector{Int}}
    plaquettes::Vector{Vector{Int}}
    boundary::Vector{Vector{Int}}
    slices::Vector{Vector{Int}}
end

"""Build the triangular periodic matching graph of the active slice."""
function build_toric3d_slice_decoder(code::Toric3DSliceLattice)
    Hq = BitMatrix(toric3d_slice_check_matrix(code))
    adjacency = [Tuple{Int,Int}[] for _ in axes(Hq,1)]
    for q in axes(Hq,2)
        ends = findall(@view Hq[:,q])
        length(ends) == 2 || error("slice qubits must touch exactly two checks")
        a,b = ends
        push!(adjacency[a],(b,q)); push!(adjacency[b],(a,q))
    end
    numpy = PythonCall.pyimport("numpy")
    pymatching = _load_pymatching()
    matcher = pymatching.Matching.from_check_matrix(_toric3d_numpy_matrix(numpy,Hq);
        weights=1.0, merge_strategy="smallest-weight",min_num_fault_ids=size(Hq,2))
    return Toric3DSliceDecoder(code,Hq,[findall(row) for row in eachrow(Hq)],
        adjacency,matcher,_toric3d_module_version(pymatching,"pymatching"),
        star_checks(code),plaquette_checks(code),toric3d_slice_boundary_checks(code),
        [toric3d_slice_qubits(code,t) for t in 0:code.steps])
end

_toric3d_slice_parity(checks, bits) =
    BitVector(isodd(count(bits[c])) for c in checks)

function _toric3d_slice_shortest_path(adjacency, source::Int, target::Int)
    predecessor = zeros(Int,length(adjacency))
    edge = zeros(Int,length(adjacency))
    queue = Int[source]
    predecessor[source] = -1
    cursor = 1
    while cursor <= length(queue) && predecessor[target] == 0
        current = queue[cursor]; cursor += 1
        for (neighbor,q) in adjacency[current]
            predecessor[neighbor] == 0 || continue
            predecessor[neighbor] = current
            edge[neighbor] = q
            push!(queue,neighbor)
        end
    end
    predecessor[target] == 0 && error("matching endpoints are disconnected")
    path = Int[]
    current = target
    while current != source
        push!(path,edge[current])
        current = predecessor[current]
    end
    return reverse!(path)
end

"""
    decode_toric3d_slice(decoder, syndrome; routing=:one_step)

Match the measured point defects and construct physical shortest paths.
`recovery` reproduces the full syndrome. `route` contains only each path's
endpoint edges for `:one_step`, or the full recovery for `:full_path`.
The partial route deliberately leaves `residual_syndrome`; it is not scored
as a failed full recovery. No boundary detector is present on this torus.
"""
function decode_toric3d_slice(decoder::Toric3DSliceDecoder,
        syndrome::AbstractVector{Bool};routing::Symbol=:one_step)
    routing in (:one_step,:full_path) ||
        throw(ArgumentError("routing must be :one_step or :full_path"))
    length(syndrome) == size(decoder.matrix,1) ||
        throw(DimensionMismatch("slice syndrome has the wrong number of checks"))
    iseven(count(syndrome)) || throw(ArgumentError("a periodic slice syndrome must have even parity"))
    full = falses(size(decoder.matrix,2))
    route = copy(full)
    pairs = Tuple{Int,Int}[]
    paths = Vector{Int}[]
    if any(syndrome)
        numpy = PythonCall.pyimport("numpy")
        matched = decoder.matcher.decode_to_matched_dets_array(_toric3d_numpy_vector(numpy,syndrome))
        raw = PythonCall.pyconvert(Vector{Vector{Int}},matched.tolist())
        endpoints = falses(length(syndrome))
        for pair in raw
            length(pair) == 2 || error("matching returned an invalid pair")
            a,b = pair .+ 1
            1 <= a <= length(syndrome) && 1 <= b <= length(syndrome) && a != b ||
                error("matching returned an invalid detector or boundary")
            endpoints[a] ⊻= true; endpoints[b] ⊻= true
            path = _toric3d_slice_shortest_path(decoder.adjacency,a,b)
            push!(pairs,(a,b)); push!(paths,path)
            for q in path
                full[q] ⊻= true
            end
            if routing === :one_step
                route[first(path)] ⊻= true
                length(path) > 1 && (route[last(path)] ⊻= true)
            end
        end
        endpoints == syndrome || error("matching did not account for all measured defects")
    end
    _toric3d_slice_parity(decoder.checks,full) == syndrome ||
        error("invalid slice recovery: its syndrome does not reproduce the measurement")
    routing === :full_path && (route .= full)
    residual = syndrome .⊻ _toric3d_slice_parity(decoder.checks,route)
    return (recovery=full,route=route,matched_pairs=pairs,paths=paths,residual_syndrome=residual)
end

"""Effect on fresh targets of a virtual repair to old controls; old data stay untouched."""
function _toric3d_slice_feedforward(code::Toric3DSliceLattice,route::AbstractVector{Bool})
    block = 3code.Ly*code.Lz
    length(route) == block || throw(DimensionMismatch("slice route has the wrong length"))
    delta = falses(block)
    for u in 0:code.Ly-1, v in 0:code.Lz-1, d in 2:3
        a,b = _toric3d_slice_pair(code,1,u,v,1,d)
        q = data_qubit_index(code,d,0,u,v)
        delta[q] = route[a] ⊻ route[b]
    end
    return delta
end

function _validate_toric3d_fault_record(encoder::Toric3DSliceEncoder,record::Toric3DFaultRecord)
    code = encoder.code
    (record.Lx,record.Ly,record.Lz) == (code.steps,code.Ly,code.Lz) ||
        throw(ArgumentError("fault record has incompatible slice dimensions"))
    size(record.x) == size(record.z) == (2,length(encoder.operations)) ||
        throw(ArgumentError("fault record has incompatible operation count"))
    for (step,operation) in enumerate(encoder.operations)
        if operation.gate === :H
            !(record.x[2,step] || record.z[2,step]) ||
                throw(ArgumentError("H gates have no second active output"))
        elseif operation.gate !== :CNOT
            throw(ArgumentError("unsupported slice encoder gate"))
        end
    end
    return nothing
end

"""Sample the same independent post-H/CNOT-output noise convention as the original scan."""
function sample_toric3d_faults(rng::Random.AbstractRNG,encoder::Toric3DSliceEncoder;
        p_x::Real=0.0,p_z::Real=0.0)
    _validate_fault_rate(p_x,"p_x"); _validate_fault_rate(p_z,"p_z")
    x,z = falses(2,length(encoder.operations)),falses(2,length(encoder.operations))
    for (step,operation) in enumerate(encoder.operations), slot in eachindex(operation.qubits)
        x[slot,step] = rand(rng) < p_x
        z[slot,step] = rand(rng) < p_z
    end
    # Lx stores the number of transfer steps for this encoder dispatch.
    return Toric3DFaultRecord(encoder.code.steps,encoder.code.Ly,encoder.code.Lz,x,z)
end

function with_toric3d_fault(encoder::Toric3DSliceEncoder,record::Toric3DFaultRecord,
        step::Integer,qubit::Integer,pauli::Symbol)
    _validate_toric3d_fault_record(encoder,record)
    pauli in (:X,:Y,:Z) || throw(ArgumentError("pauli must be X, Y, or Z"))
    1 <= step <= length(encoder.operations) || throw(ArgumentError("fault step is out of range"))
    slot = findfirst(==(Int(qubit)),encoder.operations[step].qubits)
    slot === nothing && throw(ArgumentError("qubit is not active at this gate"))
    x,z = copy(record.x),copy(record.z)
    pauli in (:X,:Y) && (x[slot,step] ⊻= true)
    pauli in (:Z,:Y) && (z[slot,step] ⊻= true)
    return Toric3DFaultRecord(record.Lx,record.Ly,record.Lz,x,z)
end

function _toric3d_slice_propagate!(frame,encoder,record,first,last)
    for step in first:last
        operation = encoder.operations[step]
        _propagate_frame_operation!(frame.x,frame.z,operation)
        for (slot,q) in enumerate(operation.qubits)
            frame.x[q] ⊻= record.x[slot,step]
            frame.z[q] ⊻= record.z[slot,step]
        end
    end
    return frame
end

function propagate_toric3d_segment!(frame::PauliFrame,encoder::Toric3DSliceEncoder,
        record::Toric3DFaultRecord,first::Integer,last::Integer)
    _validate_toric3d_fault_record(encoder,record)
    length(frame.x) == data_qubit_count(encoder.code) ||
        throw(ArgumentError("Pauli frame has the wrong size"))
    n = length(encoder.operations)
    1 <= first <= n+1 && 0 <= last <= n && first <= last+1 ||
        throw(ArgumentError("operation segment is out of range"))
    return _toric3d_slice_propagate!(frame,encoder,record,first,last)
end

"""
    replay_toric3d_slice(encoder, faults, decoder; feedback=true,
                        routing=:one_step, record_history=false)

Measure Q on the active old slice before each transfer. Feed forward the
matched route only onto fresh target qubits, after the fixed noisy gate
sequence. This is equivalent to using virtually repaired controls in that
transfer. Old physical data are not corrected; new gate faults remain.
No global final recovery is applied, and no logical-error metric is defined
for these fixed k=0 boundaries.
"""
function replay_toric3d_slice(encoder::Toric3DSliceEncoder,record::Toric3DFaultRecord,
        decoder::Toric3DSliceDecoder;feedback::Bool=true,routing::Symbol=:one_step,
        record_history::Bool=false)
    routing in (:one_step,:full_path) || throw(ArgumentError("invalid routing policy"))
    c,d = encoder.code,decoder.code
    (c.steps,c.Ly,c.Lz) == (d.steps,d.Ly,d.Lz) ||
        throw(DimensionMismatch("slice decoder and encoder dimensions differ"))
    _validate_toric3d_fault_record(encoder,record)
    frame = PauliFrame(falses(data_qubit_count(c)),falses(data_qubit_count(c)))
    history = NamedTuple[]
    first = 1
    for t in 1:c.steps
        old,new = decoder.slices[t],decoder.slices[t+1]
        syndrome = _toric3d_slice_parity(decoder.checks,frame.x[old])
        decision = feedback ? decode_toric3d_slice(decoder,syndrome;routing) : nothing
        delta = feedback ? _toric3d_slice_feedforward(c,decision.route) : nothing
        _toric3d_slice_propagate!(frame,encoder,record,first,encoder.checkpoints[t])
        first = encoder.checkpoints[t]+1
        feedback && (frame.x[new] .⊻= delta)
        if record_history
            push!(history,(layer=t,syndrome=syndrome,
                recovery=feedback ? decision.recovery : nothing,
                route=feedback ? decision.route : falses(length(old)),
                residual_syndrome=feedback ? decision.residual_syndrome : copy(syndrome),
                outgoing_syndrome=_toric3d_slice_parity(decoder.checks,frame.x[new]),
                target_correction=feedback ? delta : falses(length(new))))
        end
    end
    syndromes = (stars=_toric3d_slice_parity(decoder.stars,frame.z),
        plaquettes=_toric3d_slice_parity(decoder.plaquettes,frame.x),
        boundary=_toric3d_slice_parity(decoder.boundary,frame.x),
        frontier=_toric3d_slice_parity(decoder.checks,frame.x[last(decoder.slices)]))
    return (frame=frame,syndromes=syndromes,history=history)
end
