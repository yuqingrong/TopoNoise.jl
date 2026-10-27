"""Raw paired-shot faults. Replays treat these owned arrays as immutable."""
struct AxisRecord
    seed::BitVector
    x::BitMatrix
    z::BitMatrix
    function AxisRecord(seed::AbstractVector{Bool}, x::AbstractMatrix{Bool}, z::AbstractMatrix{Bool})
        new(BitVector(copy(seed)), BitMatrix(copy(x)), BitMatrix(copy(z)))
    end
end

"""A decoder whose only shot-dependent input is its measured syndrome."""
struct AxisKernel
    matrix::BitMatrix
    backend::Symbol
    decoder::PythonCall.Py
end

struct AxisWindow
    checks::Vector{Vector{Int}}
    qubits::Vector{Int}
    matrix::BitMatrix
    prefix_checks::Vector{Vector{Int}}
    kernel::AxisKernel
end

_axis_matrix_syndrome(matrix::BitMatrix, bits::AbstractVector{Bool}) =
    TN._toric3d_matrix_syndrome(matrix, BitVector(bits))
_axis_parity(checks, bits) = TN._toric3d_slice_parity(checks,bits)
_axis_matrix(checks,qubits) = TN._toric3d_restricted_check_matrix(checks,qubits)

function _axis_kernel(matrix::BitMatrix, backend::Symbol, p::Float64,
        settings::TN.Toric3DDecoderSettings)
    numpy = PythonCall.pyimport("numpy")
    H = TN._toric3d_numpy_matrix(numpy,matrix)
    if backend == :mwpm
        all(sum(@view matrix[:,q]) <= 2 for q in axes(matrix,2)) ||
            error("MWPM checks must have at most two endpoints per physical column")
        matcher = TN._load_pymatching().Matching.from_check_matrix(H;
            weights=1.0,merge_strategy="smallest-weight",min_num_fault_ids=size(matrix,2))
        return AxisKernel(matrix,backend,matcher)
    elseif backend == :bposd
        decoder = TN._toric3d_load_ldpc().BpOsdDecoder(H;
            error_rate=clamp(p,1e-12,1-1e-12),max_iter=settings.max_iter,
            bp_method=settings.bp_method,ms_scaling_factor=settings.ms_scaling_factor,
            schedule=settings.schedule,omp_thread_count=settings.omp_thread_count,
            osd_method=settings.osd_method,osd_order=settings.osd_order)
        return AxisKernel(matrix,backend,decoder)
    end
    throw(ArgumentError("unknown decoder backend"))
end

function _axis_decode(kernel::AxisKernel, syndrome::AbstractVector{Bool})
    length(syndrome) == size(kernel.matrix,1) || throw(DimensionMismatch("syndrome length"))
    n = size(kernel.matrix,2)
    !any(syndrome) && return falses(n)
    numpy = PythonCall.pyimport("numpy")
    prediction = kernel.decoder.decode(TN._toric3d_numpy_vector(numpy,syndrome))
    recovery = TN._toric3d_py_bits(numpy,prediction,n)
    _axis_matrix_syndrome(kernel.matrix,recovery) == syndrome ||
        error("decoder recovery does not reproduce its measured syndrome")
    return recovery
end

function _axis_adjacency(matrix::BitMatrix)
    adjacency = [Tuple{Int,Int}[] for _ in axes(matrix,1)]
    for q in axes(matrix,2)
        ends = findall(@view matrix[:,q])
        length(ends) == 2 || error("square torus edge must have two plaquette endpoints")
        a,b = ends
        push!(adjacency[a],(b,q)); push!(adjacency[b],(a,q))
    end
    # Physical edge IDs distinguish parallel period-two edges.
    foreach(sort!,adjacency)
    return adjacency
end

"""
Construct one cube's channel-specific final and causal feedback decoders.
Changing p rebuilds the BP+OSD priors; p-independent geometry/MWPM objects
are reused when previous has identical dimensions and channel.
"""
build_context(L::Integer,p::Real,channel::String;previous=nothing) =
    build_context(AxisCode(L,L,L),p,channel;previous)

function build_context(code::AxisCode,p::Real,channel::String;previous=nothing)
    code.T >= 1 || throw(ArgumentError("replay requires at least one transfer"))
    isfinite(p) && 0 <= p <= 1 || throw(ArgumentError("p must lie in [0,1]"))
    channel in ("x_only","z_only") || throw(ArgumentError("channel must be x_only or z_only"))
    prior = Float64(p)
    same = previous !== nothing && previous.channel == channel &&
        (previous.code.T,previous.code.Ly,previous.code.Lz) == (code.T,code.Ly,code.Lz)
    if same && (channel == "z_only" || previous.p == prior)
        return merge(previous,(p=prior,))
    end
    encoder = same ? previous.encoder : axis_encoder(code)
    settings = TN.Toric3DDecoderSettings()
    backend = channel == "x_only" ? :bposd : :mwpm
    complete = channel == "x_only" ? plaquette_checks(code) : star_checks(code)
    final_matrix = same ? previous.final.matrix : _axis_matrix(complete,collect(1:nqubits(code)))
    final = _axis_kernel(final_matrix,backend,prior,settings)
    seedcode = AxisCode(0,code.Ly,code.Lz)
    seed_checks = star_checks(seedcode)
    seed_kernel = channel == "z_only" ? _axis_kernel(
        _axis_matrix(seed_checks,slice_qubits(code,0)),:mwpm,prior,settings) : nothing
    slice_checks = plaquette_checks(seedcode)
    slice_kernel = channel == "x_only" ? (same ? previous.slice_kernel : _axis_kernel(
        _axis_matrix(slice_checks,plane_edges(code,0)),:mwpm,prior,settings)) : nothing
    adjacency = channel == "x_only" ? (same ? previous.adjacency : _axis_adjacency(slice_kernel.matrix)) : nothing
    windows = AxisWindow[]
    cached = AxisKernel[]
    for r in 1:code.T
        if same
            old = previous.windows[r]
            checks,qubits,matrix,prefix_checks = old.checks,old.qubits,old.matrix,old.prefix_checks
        else
            prefix = AxisCode(r,code.Ly,code.Lz)
            prefix_checks = channel == "x_only" ? plaquette_checks(prefix) : star_checks(prefix)
            qubits = vcat(slice_qubits(code,r-1),slice_qubits(code,r))
            allowed = Set(qubits)
            # Keep every incident physical row, including the older endpoint of
            # incoming x links. Omitted rows are used only in the prefix audit.
            checks = filter(check -> any(q -> q in allowed,check),prefix_checks)
            matrix = _axis_matrix(checks,qubits)
        end
        index = findfirst(kernel -> kernel.matrix == matrix,cached)
        kernel = isnothing(index) ? _axis_kernel(matrix,backend,prior,settings) : cached[index]
        isnothing(index) && push!(cached,kernel)
        push!(windows,AxisWindow(checks,qubits,matrix,prefix_checks,kernel))
    end
    logicals = channel == "x_only" ? logical_z(code) : logical_x(code)
    N = code.Ly*code.Lz
    push_windows = channel == "z_only" ? [(
        checks=windows[r].prefix_checks[(r-1)*N+1:r*N],
        retired_checks=windows[r].prefix_checks[1:r*N],
        qubits=[edge(code,:x,r,u,v) for u in 0:code.Ly-1 for v in 0:code.Lz-1])
        for r in 1:code.T-1] : nothing
    return (code=code,encoder=encoder,p=prior,channel=channel,settings=settings,
        final=final,final_checks=complete,logicals=logicals,seed_checks=seed_checks,
        seed_kernel=seed_kernel,slice_checks=slice_checks,slice_kernel=slice_kernel,
        adjacency=adjacency,windows=windows,push_windows=push_windows)
end

"""One Bernoulli draw per active selected-channel gate output; Z seed first."""
function sample_record(rng::Random.AbstractRNG,context)
    count = length(context.encoder.operations)
    seed = falses(2context.code.Ly*context.code.Lz)
    x,z = falses(2,count),falses(2,count)
    if context.channel == "z_only"
        for q in eachindex(seed)
            seed[q] = rand(rng) < context.p
        end
    end
    selected = context.channel == "x_only" ? x : z
    for (index,operation) in enumerate(context.encoder.operations), slot in eachindex(operation.qubits)
        selected[slot,index] = rand(rng) < context.p
    end
    return AxisRecord(seed,x,z)
end

function _axis_validate_record(context,record::AxisRecord)
    length(record.seed) == 2context.code.Ly*context.code.Lz || throw(DimensionMismatch("seed record length"))
    expected = (2,length(context.encoder.operations))
    size(record.x) == expected && size(record.z) == expected || throw(DimensionMismatch("gate record shape"))
    if context.channel == "x_only"
        (!any(record.seed) && !any(record.z)) || throw(ArgumentError("X replay requires ideal seed and X-only outputs"))
    else
        !any(record.x) || throw(ArgumentError("Z replay requires Z-only outputs"))
    end
    for (index,operation) in enumerate(context.encoder.operations)
        if length(operation.qubits) == 1
            (!record.x[2,index] && !record.z[2,index]) || throw(ArgumentError("unused one-qubit gate slot is nonzero"))
        end
    end
    return nothing
end

function _axis_endpoint_route(context,syndrome::AbstractVector{Bool})
    kernel = context.slice_kernel
    length(syndrome) == size(kernel.matrix,1) || throw(DimensionMismatch("square syndrome length"))
    iseven(count(syndrome)) || throw(ArgumentError("closed square plane syndrome must have even parity"))
    full,route = falses(size(kernel.matrix,2)),falses(size(kernel.matrix,2))
    paths,pairs = Vector{Int}[],Tuple{Int,Int}[]
    if any(syndrome)
        numpy = PythonCall.pyimport("numpy")
        matched = kernel.decoder.decode_to_matched_dets_array(TN._toric3d_numpy_vector(numpy,syndrome))
        raw = PythonCall.pyconvert(Vector{Vector{Int}},matched.tolist())
        accounted = falses(length(syndrome))
        for pair in raw
            length(pair) == 2 || error("matching returned malformed pair")
            a,b = pair .+ 1
            1 <= a <= length(syndrome) && 1 <= b <= length(syndrome) && a != b ||
                error("square matching returned invalid endpoint")
            accounted[a] ⊻= true; accounted[b] ⊻= true
            path = TN._toric3d_slice_shortest_path(context.adjacency,a,b)
            push!(pairs,(a,b)); push!(paths,path)
            for q in path
                full[q] ⊻= true
            end
            route[first(path)] ⊻= true
            length(path) > 1 && (route[last(path)] ⊻= true)
        end
        accounted == syndrome || error("matching did not account for all measured defects")
    end
    _axis_matrix_syndrome(kernel.matrix,full) == syndrome || error("physical route has wrong syndrome")
    return (recovery=full,route=route,paths=paths,matched_pairs=pairs)
end

function _axis_propagate!(frame,context,record,first,last)
    for index in first:last
        operation = context.encoder.operations[index]
        TN._propagate_frame_operation!(frame.x,frame.z,operation)
        for (slot,q) in enumerate(operation.qubits)
            frame.x[q] ⊻= record.x[slot,index]
            frame.z[q] ⊻= record.z[slot,index]
        end
    end
    return nothing
end

function _axis_apply_two!(frame,context,r;record_history::Bool)
    window = context.windows[r]
    bits = context.channel == "x_only" ? frame.x : frame.z
    measured = _axis_parity(window.checks,bits)
    recovery = _axis_decode(window.kernel,measured)
    for (column,q) in enumerate(window.qubits)
        bits[q] ⊻= recovery[column]
    end
    post = _axis_parity(window.checks,bits)
    any(post) && error("two-slice correction left a measured defect")
    prefix = _axis_parity(window.prefix_checks,bits)
    any(prefix) && error("two-slice correction left a completed-prefix defect")
    record_history || return nothing
    return (kind=:two_slice,checkpoint=r,checks=deepcopy(window.checks),
        post_checks=deepcopy(window.checks),qubits=copy(window.qubits),
        recovery=recovery,support=window.qubits[findall(recovery)],
        pre_syndrome=measured,post_syndrome=post,prefix_syndrome=prefix,
        logical_after=_axis_parity(context.logicals,bits))
end

function _axis_apply_z_push!(frame,context,r;record_history::Bool)
    window = context.push_windows[r]
    checks = window.checks
    measured = _axis_parity(checks,frame.z)
    qubits = window.qubits
    # Each outgoing x link meets exactly one retiring star; its other endpoint
    # is at the frontier. Odd retiring parity is valid and is pushed forward.
    for (index,q) in enumerate(qubits)
        frame.z[q] ⊻= measured[index]
    end
    retired = _axis_parity(window.retired_checks,frame.z)
    any(retired) && error("x-directed Z push left a retired-plane defect")
    record_history || return nothing
    return (kind=:z_push,checkpoint=r,checks=deepcopy(checks),post_checks=deepcopy(checks),qubits=copy(qubits),
        recovery=copy(measured),support=qubits[findall(measured)],pre_syndrome=measured,
        post_syndrome=_axis_parity(checks,frame.z),prefix_syndrome=retired,
        logical_after=_axis_parity(context.logicals,frame.z))
end

function _axis_replay(context,record::AxisRecord,protocol::Symbol,seed_entry;record_history::Bool)
    code = context.code
    frame = TN.PauliFrame(falses(nqubits(code)),falses(nqubits(code)))
    frame.z[slice_qubits(code,0)] .= record.seed
    history = NamedTuple[]
    record_history && seed_entry !== nothing && push!(history,deepcopy(seed_entry))
    firstop = 1
    for r in 1:code.T
        route = nothing
        measured = nothing
        checks = nothing
        if context.channel == "x_only" && protocol == :one_step
            old = plane_edges(code,r-1)
            checks = [[old[q] for q in check] for check in context.slice_checks]
            measured = _axis_parity(checks,frame.x)
            route = _axis_endpoint_route(context,measured)
        end
        _axis_propagate!(frame,context,record,firstop,context.encoder.checkpoints[r])
        firstop = context.encoder.checkpoints[r]+1
        if context.channel == "x_only" && protocol == :one_step
            fresh = plane_edges(code,r)
            for (column,q) in enumerate(fresh)
                frame.x[q] ⊻= route.route[column]
            end
            if record_history
                post_checks = [[fresh[q] for q in check] for check in context.slice_checks]
                push!(history,(kind=:x_endpoint,checkpoint=r,checks=checks,
                    post_checks=post_checks,qubits=fresh,recovery=route.route,
                    support=fresh[findall(route.route)],pre_syndrome=measured,
                    full_recovery=route.recovery,paths=route.paths,matched_pairs=route.matched_pairs,
                    post_syndrome=_axis_parity(post_checks,frame.x),
                    logical_after=_axis_parity(context.logicals,frame.x)))
            end
        elseif context.channel == "z_only" && protocol == :one_step && r < code.T
            entry = _axis_apply_z_push!(frame,context,r;record_history)
            record_history && push!(history,entry)
        elseif protocol == :two_slice
            entry = _axis_apply_two!(frame,context,r;record_history)
            record_history && push!(history,entry)
        end
    end
    bits = context.channel == "x_only" ? frame.x : frame.z
    syndrome = _axis_parity(context.final_checks,bits)
    recovery = _axis_decode(context.final,syndrome)
    residual = bits .⊻ recovery
    flags = _axis_parity(context.logicals,residual)
    return (frame=frame,syndrome=syndrome,recovery=recovery,residual=residual,
        logical_1=flags[1],logical_2=flags[2],failure=any(flags),history=history)
end

"""
Replay identical raw faults through final-only, orientation-adapted one-step,
and complete two-slice feedback. Z seed repair is computed exactly once,
shared by all arms, and never discarded when it contains a logical error.
Neither the raw record nor its arrays are mutated by any replay.
"""
function three_outcomes(context,raw::AxisRecord;record_history::Bool=false)
    _axis_validate_record(context,raw)
    record = raw
    seed_entry = nothing
    if context.channel == "z_only"
        measured = _axis_matrix_syndrome(context.seed_kernel.matrix,raw.seed)
        recovery = _axis_decode(context.seed_kernel,measured)
        corrected = raw.seed .⊻ recovery
        record = AxisRecord(corrected,raw.x,raw.z)
        if record_history
            qubits = slice_qubits(context.code,0)
            seed_entry = (kind=:seed,checkpoint=0,checks=deepcopy(context.seed_checks),
                post_checks=deepcopy(context.seed_checks),qubits=qubits,recovery=recovery,
                support=qubits[findall(recovery)],pre_syndrome=measured,
                post_syndrome=_axis_matrix_syndrome(context.seed_kernel.matrix,corrected))
        end
    end
    final_only = _axis_replay(context,record,:final_only,seed_entry;record_history)
    one_step = _axis_replay(context,record,:one_step,seed_entry;record_history)
    two_slice = _axis_replay(context,record,:two_slice,seed_entry;record_history)
    return (final_only=final_only,one_step=one_step,two_slice=two_slice)
end
