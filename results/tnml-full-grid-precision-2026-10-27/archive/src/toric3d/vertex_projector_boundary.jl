"""Ideal lower preparation and fixed zero-logical-parity upper correction tree."""
struct Toric3DVertexBoundaryMap
    code::Toric3DVertexProjectorCode
    root::Int
    parent::Vector{Int}
    parent_edge::Vector{Int}
    paths::Vector{Vector{Int}}
    upper_checks::Vector{Vector{Int}}
    lower_qubits::Vector{Int}
end

function toric3d_vertex_boundary_map(c::Toric3DVertexProjectorCode)
    checks=star_checks(c);m=c.Lu*c.Lv;n=data_qubit_count(c)
    endpoints=[Int[] for _ in 1:n]
    for (r,s) in enumerate(checks),q in s;push!(endpoints[q],r);end
    excluded=Set(vcat(toric3d_encoded_logical_x(c)...))
    adj=[Tuple{Int,Int}[] for _ in checks]
    for q in 1:n
        length(endpoints[q])==2 || error("boundary tree requires two star endpoints")
        q in excluded && continue
        a,b=endpoints[q];push!(adj[a],(q,b));push!(adj[b],(q,a))
    end
    root=c.T*m+1;parent=zeros(Int,length(checks));edge=zeros(Int,length(checks))
    parent[root]=root;queue=[root]
    for a in queue
        for (q,b) in adj[a] # physical IDs are inserted in increasing order
            parent[b]!=0 && continue
            parent[b]=a;edge[b]=q;push!(queue,b)
        end
    end
    all(!iszero,parent) || error("logical-cut boundary graph is disconnected")
    paths=Vector{Int}[]
    for upper in (c.T+1)*m+1:length(checks)
        path=Int[];v=upper
        while v!=root
            push!(path,edge[v]);v=parent[v]
        end
        push!(paths,path)
    end
    lower=[data_qubit_index(c,3,0,u,v) for u in 0:c.Lu-1 for v in 0:c.Lv-1]
    Toric3DVertexBoundaryMap(c,root,parent,edge,paths,checks[end-m+1:end],lower)
end

"""Only the physically observed upper outcomes enter this controller."""
function toric3d_vertex_boundary_correction(map::Toric3DVertexBoundaryMap,
        a::AbstractVector{Bool})::BitVector
    length(a)==length(map.paths) || throw(DimensionMismatch("upper outcome length differs"))
    correction=falses(data_qubit_count(map.code))
    for i in findall(a),q in map.paths[i];correction[q]⊻=true;end
    correction
end

"""Uniform clean even branch; draw the first M-1 bits then fix parity."""
function sample_toric3d_vertex_upper_branch(rng::Random.AbstractRNG,map::Toric3DVertexBoundaryMap)::BitVector
    b=falses(length(map.paths))
    for i in 1:length(b)-1;b[i]=rand(rng,Bool);end
    b[end]=isodd(count(b))
    b
end

struct Toric3DVertexBoundaryCompletion
    z::BitVector
    clean_branch::BitVector
    observed::BitVector
    clean_correction::BitVector
    physical_correction::BitVector
end

"""Simulator comparison of noisy and clean physical measurement branches."""
function complete_toric3d_vertex_boundary(map::Toric3DVertexBoundaryMap,
        z::AbstractVector{Bool},b::AbstractVector{Bool})
    length(z)==data_qubit_count(map.code) || throw(DimensionMismatch("boundary frame length differs"))
    length(b)==length(map.paths) || throw(DimensionMismatch("clean branch length differs"))
    iseven(count(b)) || throw(ArgumentError("the clean upper branch must have even parity"))
    observed=b .⊻ _vertex_syndrome(map.upper_checks,z)
    physical=toric3d_vertex_boundary_correction(map,observed)
    clean=toric3d_vertex_boundary_correction(map,b)
    Toric3DVertexBoundaryCompletion(z .⊻ physical .⊻ clean,BitVector(b),observed,clean,physical)
end

"""Check independent lower GHZ projections and their zero-parity Z corrections."""
function verify_toric3d_vertex_lower_boundary(map::Toric3DVertexBoundaryMap)::Bool
    c=map.code;m=c.Lu*c.Lv;w=star_checks(c)[1:m]
    sort(vcat(w...))==collect(1:3m) || return false
    length(unique(map.lower_qubits))==m || return false
    all(map.lower_qubits[i] in w[i] for i in 1:m) || return false
    cut=Set(vcat(toric3d_encoded_logical_x(c)...))
    all(q->q ∉ cut,map.lower_qubits) || return false
    ops=EncoderOperation[]
    all(s->_vertex_prepared(c,ops,:X,s),w) || return false
    all(s->_vertex_prepared(c,ops,:Z,s),toric3d_encoded_logical_z(c))
end
