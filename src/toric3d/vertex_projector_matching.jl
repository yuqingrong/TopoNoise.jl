"""A scientific recovery invariant failed; user interruption is a separate exception."""
struct Toric3DVertexInvalidRecovery <: Exception
    message::String
    context::Dict{String,Any}
end
Base.showerror(io::IO,e::Toric3DVertexInvalidRecovery)=print(io,e.message)
function _vertex_invalid(w,s,recovery;message="invalid H+5 recovery: returned data do not reproduce measured syndrome",route=recovery,paths=Any[],observed_post_syndrome=nothing)
    throw(Toric3DVertexInvalidRecovery(message,Dict{String,Any}("channel"=>string(w.channel),"policy"=>string(w.policy),"checkpoint"=>w.checkpoint,
        "syndrome"=>collect(s),"recovery"=>collect(recovery),"post_syndrome"=>s .⊻ _vertex_syndrome(w.matrix,route),"predicted_post_syndrome"=>s .⊻ _vertex_syndrome(w.matrix,route),"observed_post_syndrome"=>isnothing(observed_post_syndrome) ? Bool[] : collect(observed_post_syndrome),"qubits"=>w.qubits,"checks"=>w.checks,"route"=>collect(route),"paths"=>paths)))
end

"""A matched endpoint pair and its ordered full-code physical data edges."""
struct Toric3DVertexMatchedPath
    first::Int
    second::Int
    edges::Vector{Int}
end
struct Toric3DVertexMatchingDecision
    paths::Vector{Toric3DVertexMatchedPath}
    full_recovery::BitVector
    route::BitVector
    residual_syndrome::BitVector
    signature::Vector{Tuple{Int,Int,Tuple}}
    backend::String
    version::String
end
struct Toric3DVertexMatchingDecoder
    window::Toric3DVertexFeedbackWindow
    mode::Symbol
    adjacency::Vector{Vector{Tuple{Int,Int}}}
    matcher::PythonCall.Py
    numpy::PythonCall.Py
    oracle_nodes::Int
    cost_cache::Dict{String,Int}
    version::String
end

function _vertex_matching_decoder(w::Toric3DVertexFeedbackWindow,mode::Symbol)
    mode in (:closed,:single_future,:open) || throw(ArgumentError("invalid matching mode"))
    pm=_load_pymatching();version=_toric3d_module_version(pm,"pymatching")
    version=="2.4.0" || error("H+5 matching requires pymatching 2.4.0; found $version")
    np=PythonCall.pyimport("numpy");matcher=pm.Matching()
    nr,nq=size(w.matrix);nv=nr+(mode==:closed ? 0 : 1)
    adj=[Tuple{Int,Int}[] for _ in 1:nv]
    for j in 1:nq
        ends=findall(@view w.matrix[:,j])
        length(ends) in (1,2) || error("a matching column needs one or two detector endpoints")
        outward=length(ends)==1
        outward && mode==:closed && error("closed matcher has an open edge")
        a=ends[1];b=outward ? nr+1 : ends[2]
        q=w.qubits[j];push!(adj[a],(q,b));push!(adj[b],(q,a))
        # Unique degree-two subdivision node preserves parallel physical edges.
        # Integral weights avoid floating-point tie perturbations.
        aux=nv+j
        weight=mode==:single_future && outward ? nr+1 : 1
        matcher.add_edge(a-1,aux-1;weight=weight)
        matcher.add_edge(aux-1,b-1;weight=weight)
    end
    foreach(sort!,adj)
    mode==:open && matcher.set_boundary_nodes(PythonCall.pyimport("builtins").set([nr]))
    Toric3DVertexMatchingDecoder(w,mode,adj,matcher,np,nv+nq,Dict{String,Int}(),version)
end
build_toric3d_vertex_x_one_step_decoder(c::Toric3DVertexProjectorCode,t::Integer)=
    _vertex_matching_decoder(toric3d_vertex_x_slice_window(c,t),:closed)
build_toric3d_vertex_z_one_step_decoder(c::Toric3DVertexProjectorCode,t::Integer)=
    _vertex_matching_decoder(toric3d_vertex_z_strip_window(c,t),:single_future)

function _vertex_matching_cost(d,s)
    !any(s) && return 0
    key=String(UInt8.(s))
    haskey(d.cost_cache,key) && return d.cost_cache[key]
    nr=length(s);bits=zeros(UInt8,d.oracle_nodes);bits[1:nr].=s
    parity=isodd(count(s))
    d.mode==:closed && parity && throw(ArgumentError("closed matching needs even syndrome"))
    d.mode==:single_future && (bits[nr+1]=parity)
    decoded=d.matcher.decode(d.numpy.asarray(bits;dtype=d.numpy.uint8);return_weight=true)
    weight=PythonCall.pyconvert(Float64,decoded[1])
    doubled=round(Int,weight)
    weight==doubled && iseven(doubled) || error("matching oracle returned a nonintegral physical weight")
    cost=doubled÷2-(d.mode==:single_future && parity ? nr : 0)
    length(d.cost_cache)>=4096 && empty!(d.cost_cache)
    d.cost_cache[key]=cost
    cost
end

# BFS visits paths in lexicographic physical-edge order. A virtual endpoint
# may terminate a path but may never serve as an interior shortcut.
function _vertex_bfs(d,source)
    nr=size(d.window.matrix,1);nv=length(d.adjacency)
    parent=zeros(Int,nv);edge=zeros(Int,nv);distance=fill(typemax(Int),nv)
    parent[source]=source;distance[source]=0;queue=[source]
    for a in queue
        a>nr && continue
        for (q,b) in d.adjacency[a]
            parent[b]!=0 && continue
            parent[b]=a;edge[b]=q;distance[b]=distance[a]+1;push!(queue,b)
        end
    end
    parent,edge,distance
end
function _vertex_lift(parent,edge,a,b)
    parent[b]!=0 || error("matching endpoints are disconnected")
    path=Int[]
    while b!=a;push!(path,edge[b]);b=parent[b];end
    reverse!(path)
end

function _vertex_validate_recovery(w,s,recovery;paths=Any[])
    length(s)==size(w.matrix,1) || throw(DimensionMismatch("measured syndrome length differs"))
    length(recovery)==length(w.qubits) || throw(DimensionMismatch("data recovery length differs"))
    _vertex_syndrome(w.matrix,recovery)==s || _vertex_invalid(w,s,recovery;paths)
    recovery
end

function _vertex_matching_decision(d::Toric3DVertexMatchingDecoder,s::AbstractVector{Bool};trim::Symbol=:complete)
    w=d.window;nr=size(w.matrix,1)
    length(s)==nr || throw(DimensionMismatch("matching syndrome length differs"))
    d.mode==:closed && isodd(count(s)) && throw(ArgumentError("closed matching needs even syndrome"))
    trim in (:complete,:x_endpoints,:z_retiring) || throw(ArgumentError("invalid path trim"))
    remaining=BitVector(s);cost=_vertex_matching_cost(d,remaining)
    paths=Toric3DVertexMatchedPath[]
    while any(remaining)
        a=findfirst(remaining);parent,edge,distance=_vertex_bfs(d,a)
        candidates=findall(remaining);popfirst!(candidates)
        if d.mode==:open || (d.mode==:single_future && isodd(count(remaining)))
            push!(candidates,nr+1)
        end
        selected=0
        for b in candidates
            distance[b]==typemax(Int) && continue
            distance[b]>cost && continue
            remaining[a]=false;b<=nr && (remaining[b]=false)
            rest=_vertex_matching_cost(d,remaining)
            remaining[a]=true;b<=nr && (remaining[b]=true)
            if distance[b]+rest==cost
                selected=b;break
            end
        end
        selected!=0 || error("canonical matching could not realize primary optimum")
        b=selected;push!(paths,Toric3DVertexMatchedPath(a,b,_vertex_lift(parent,edge,a,b)))
        remaining[a]=false;b<=nr && (remaining[b]=false);cost-=distance[b]
    end
    cost==0 || error("canonical matching weight accounting failed")
    index=Dict(q=>j for (j,q) in enumerate(w.qubits));full=falses(length(w.qubits));route=copy(full)
    retired=Set(w.retired_rows)
    for p in paths
        for q in p.edges;full[index[q]]⊻=true;end
        if trim==:x_endpoints
            route[index[first(p.edges)]]⊻=true
            length(p.edges)>1 && (route[index[last(p.edges)]]⊻=true)
        elseif trim==:z_retiring
            p.first in retired && (route[index[first(p.edges)]]⊻=true)
            p.second in retired && (route[index[last(p.edges)]]⊻=true)
        end
    end
    _vertex_validate_recovery(w,s,full;paths=[p.edges for p in paths])
    trim==:complete && (route.=full)
    residual=s .⊻ _vertex_syndrome(w.matrix,route)
    if trim==:z_retiring
        !any(residual[w.retired_rows]) || error("causal Z transport left a retiring detector")
        all(q->data_qubit_coordinate(w.code,q)[2]==w.checkpoint-1,w.qubits[findall(route)]) ||
            error("causal Z transport escaped retiring data")
    end
    signature=Tuple{Int,Int,Tuple}[(p.first,p.second,Tuple(p.edges)) for p in paths]
    Toric3DVertexMatchingDecision(paths,full,route,BitVector(residual),signature,"pymatching",d.version)
end
function decode_toric3d_vertex_x_one_step(d::Toric3DVertexMatchingDecoder,s::AbstractVector{Bool})
    d.window.channel==:x_only && d.window.policy==:one_step || throw(ArgumentError("not an X one-step decoder"))
    _vertex_matching_decision(d,s;trim=:x_endpoints)
end
function decode_toric3d_vertex_z_one_step(d::Toric3DVertexMatchingDecoder,s::AbstractVector{Bool})
    d.window.channel==:z_only && d.window.policy==:one_step || throw(ArgumentError("not a Z one-step decoder"))
    _vertex_matching_decision(d,s;trim=:z_retiring)
end
