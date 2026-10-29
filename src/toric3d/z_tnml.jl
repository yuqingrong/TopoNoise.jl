# Numerical contraction adapted from the validated RepairedBoundaryMPS module.
# Reference SHA256: 0321bc105503ccc642ecc8c33bbf9627e6d50f816b715a4dc24b89ee30db1f53

"""
    Toric3DZTNSettings(; max_bond_dimension=240, cutoff=1e-18,
                        negative_tolerance=1e-10)

Controls for syndrome-conditioned boundary-MPS contraction. The defaults match
an empirically validated L=3 configuration, not an error bound or a claim of
convergence at larger sizes. `cutoff` controls local SVD truncation;
`negative_tolerance` only controls numerical admissibility of final weights.
"""
struct Toric3DZTNSettings
    max_bond_dimension::Int
    cutoff::Float64
    negative_tolerance::Float64

    function Toric3DZTNSettings(max_bond_dimension::Integer, cutoff::Real,
            negative_tolerance::Real)
        1 <= max_bond_dimension <= typemax(Int) || throw(ArgumentError(
            "max_bond_dimension must be positive and fit Int"))
        isfinite(cutoff) && 0 <= cutoff < 1 || throw(ArgumentError(
            "cutoff must be finite and in [0,1)"))
        isfinite(negative_tolerance) && 0 <= negative_tolerance < 1 ||
            throw(ArgumentError("negative tolerance must be finite and in [0,1)"))
        return new(Int(max_bond_dimension), Float64(cutoff), Float64(negative_tolerance))
    end
end

Toric3DZTNSettings(; max_bond_dimension::Integer=240, cutoff::Real=1e-18,
    negative_tolerance::Real=1e-10) =
    Toric3DZTNSettings(max_bond_dimension, cutoff, negative_tolerance)

struct _Toric3DZTNNeighborGroup
    neighbor::Int
    weights::Vector{Float64}
end

"""Compiled exact fault-star channels and boundary-MPS contraction settings."""
struct Toric3DZTNDecoder
    detector_count::Int
    groups::Vector{Vector{_Toric3DZTNNeighborGroup}}
    settings::Toric3DZTNSettings
    edge_count::Int
end

# Compile exact outgoing-star channels; no SVD or approximation occurs here.
function _build_toric3d_z_tn_edges(n::Int, edges, qs, logical; settings::Toric3DZTNSettings=Toric3DZTNSettings())
    n > 0 || throw(ArgumentError("detector count must be positive"))
    length(edges) == length(qs) == length(logical) || throw(DimensionMismatch("edge arrays differ in length"))
    buckets = [Dict{Int,Vector{Float64}}() for _ in 1:n]
    for (edge,q0,g0) in zip(edges,qs,logical)
        u,v = edge
        1 <= u < v <= n || throw(ArgumentError("edges must have ordered, distinct endpoints"))
        q=Float64(q0); g=Int(g0)
        isfinite(q) && 0 <= q <= 1 || throw(ArgumentError("invalid edge probability"))
        0 <= g <= 3 || throw(ArgumentError("logical masks must be in 0:3"))
        old=get!(buckets[u],v) do
            z=zeros(8);z[1]=1.;z
        end
        delta=1+2g
        buckets[u][v]=[(1-q)*old[h+1]+q*old[xor(h,delta)+1] for h in 0:7]
    end
    groups=[[_Toric3DZTNNeighborGroup(w,buckets[v][w]) for w in sort!(collect(keys(buckets[v])))] for v in 1:n]
    Toric3DZTNDecoder(n,groups,settings,length(edges))
end

"""
    build_toric3d_z_tn_decoder(source, p; settings::Toric3DZTNSettings=Toric3DZTNSettings())

Build an approximate maximum-likelihood logical-class decoder for the diagonal
H+6CNOT preparation with independent Z gate-output faults and no intermediate
feedback. `source` may be a `Toric3DEncodedSlice`, its fault model, or a
`Toric3DZTNTopology`. The physical fault probability must satisfy `0 <= p < 0.5`.

Fault propagation and local star MPOs are exact; only boundary-MPS compression
is approximate. The defaults were validated for L=3 at p=0.006, 0.01, and 0.02.
Other sizes or rates require their own convergence checks. Building a decoder
does not sample faults or alter global BLAS/thread settings.
"""
function build_toric3d_z_tn_decoder(topology::Toric3DZTNTopology, p::Real;
        settings::Toric3DZTNSettings=Toric3DZTNSettings())
    n = size(topology.raw_detector_matrix, 1)
    qs = toric3d_z_tn_edge_probabilities(topology, p)
    return _build_toric3d_z_tn_edges(n, topology.detector_endpoints, qs,
        topology.logical_masks; settings)
end

build_toric3d_z_tn_decoder(model::Toric3DZFaultModel, p::Real; kwargs...) =
    build_toric3d_z_tn_decoder(build_toric3d_z_tn_topology(model), p; kwargs...)

build_toric3d_z_tn_decoder(code::Toric3DEncodedSlice, p::Real; kwargs...) =
    build_toric3d_z_tn_decoder(build_toric3d_z_fault_model(code), p; kwargs...)

# Exact suffix support of the three-bit (detector parity, logical xor) register.
function _toric3d_z_tn_suffix_support(groups)
    suffix=Vector{Vector{Int}}(undef,length(groups)+1)
    suffix[end]=[0]
    for k in length(groups):-1:1
        support=findall(!iszero,groups[k].weights).-1
        suffix[k]=sort!(collect(Set(xor(a,b) for a in support for b in suffix[k+1])))
    end
    suffix
end

function _toric3d_z_tn_fuse_current(psi)
    ii=ITensorMPS.siteinds(psi)
    C=ITensors.combiner(ii[2],ii[1];tags="Carrier")
    # With detector first in the combiner: c = parity + 2*logical.
    first=psi[1]*psi[2]*C
    # The untouched suffix remains right orthogonal. Recording that gauge is
    # what lets the later prefix truncation avoid sweeping the product tail.
    ITensorMPS.MPS(vcat(ITensor[first],[psi[k] for k in 3:length(psi)]);
        ortho_lims=1:1)
end

function _toric3d_z_tn_star_mpo(psi,groups,vertex,syndrome_bit)
    sites=ITensorMPS.siteinds(psi)
    output_logical=Index(4,"Logical")
    if isempty(groups)
        arr=zeros(4,8)
        for ell in 0:3
            arr[ell+1,Int(syndrome_bit)+2ell+1]=1.
        end
        return ITensorMPS.MPO(ITensor[ITensor(arr,output_logical,sites[1])])
    end
    lastsite=groups[end].neighbor-vertex+1
    tensors=Vector{ITensor}(undef,lastsite)
    suffix=_toric3d_z_tn_suffix_support(groups)
    left_states=suffix[1]
    previous=Index(length(left_states),"Star,v=$vertex,k=1")
    arr=zeros(4,8,length(left_states))
    lookup=Dict(h=>i for (i,h) in enumerate(left_states))
    for ellout in 0:3,c in 0:7
        h=xor(c,Int(syndrome_bit)+2ellout)
        hpos=get(lookup,h,0)
        iszero(hpos) || (arr[ellout+1,c+1,hpos]=1.)
    end
    tensors[1]=ITensor(arr,output_logical,sites[1],previous)
    group_index=1
    for k in 2:lastsite
        neighbor=vertex+k-1
        is_neighbor=groups[group_index].neighbor==neighbor
        weights=is_neighbor ? groups[group_index].weights : [1.,0.,0.,0.,0.,0.,0.,0.]
        right_states=is_neighbor ? suffix[group_index+1] : left_states
        if k==lastsite
            arr=zeros(length(left_states),2,2)
            @assert right_states==[0]
            for (a,h) in enumerate(left_states),b in 0:1
                arr[a,xor(b,h&1)+1,b+1]=weights[h+1]
            end
            tensors[k]=ITensor(arr,previous,ITensors.prime(sites[k]),sites[k])
        else
            next=Index(length(right_states),"Star,v=$vertex,k=$k")
            arr=zeros(length(left_states),length(right_states),2,2)
            for (a,h) in enumerate(left_states),(b,hnext) in enumerate(right_states),bit in 0:1
                delta=xor(h,hnext)
                arr[a,b,xor(bit,delta&1)+1,bit+1]=weights[delta+1]
            end
            tensors[k]=ITensor(arr,previous,next,ITensors.prime(sites[k]),sites[k])
            previous=next
        end
        left_states=right_states
        is_neighbor && (group_index+=1)
    end
    ITensorMPS.MPO(tensors)
end

# Apply a finite-support star MPO without materializing identity tensors over
# the unchanged suffix. Combining the MPO and MPS bond indices reproduces the
# exact, untruncated naive MPO-MPS product on this prefix.
function _toric3d_z_tn_apply_prefix(psi,star)
    support=length(star)
    1 <= support <= length(psi) || throw(ArgumentError(
        "star support must be a nonempty MPS prefix"))
    tensors=ITensor[star[k]*psi[k] for k in 1:support]
    for bond in 1:(support-1)
        star_links=ITensors.commoninds(star[bond],star[bond+1])
        psi_links=ITensors.commoninds(psi[bond],psi[bond+1])
        links=Index[star_links...,psi_links...]
        isempty(links) && continue
        C=ITensors.combiner(links)
        tensors[bond]*=C
        tensors[bond+1]*=ITensors.dag(C)
    end
    append!(tensors,(psi[k] for k in (support+1):length(psi)))
    return ITensorMPS.MPS(tensors;ortho_lims=1:support)
end

function _toric3d_z_tn_mass(psi)
    value=ITensor(1.)
    for (T,s) in zip(psi,ITensorMPS.siteinds(psi))
        value*=T*ITensor(ones(ITensors.dim(s)),s)
    end
    Float64(value[])
end

"""
    decode_toric3d_z_tn(decoder, syndrome)

Condition the fault tensor network on a measured binary final-star syndrome
and return an approximate maximum-likelihood logical-class decision. The four
`probabilities` are ordered by mask `ell1 + 2ell2`: `(0,0), (1,0), (0,1), (1,1)`.
`logical_mask` is the argmax mask in `0:3`, or `-1` when the output is invalid.
This is a logical-class decision, not a physical recovery vector.

The result also contains `log_syndrome`, `raw_probabilities`, and `diagnostics`.
Signed normalized weights are preserved without clipping or an exact fallback.
`valid` checks only numerical admissibility against `negative_tolerance`; it
is not a posterior-accuracy or convergence certificate. The diagnostic
`discarded_weight` is a sum of local relative SVD truncation errors, not lost
probability or a bound on posterior error.

Incorrect syndrome length raises `DimensionMismatch`; nonbinary or odd-parity
syndromes raise `ArgumentError`. A zero/nonfinite conditioned mass raises
`DomainError`. No fault sampling, physical feedback, or global thread changes
are performed.
"""
function decode_toric3d_z_tn(decoder::Toric3DZTNDecoder,syndrome::AbstractVector)
    n=decoder.detector_count
    length(syndrome)==n || throw(DimensionMismatch("incorrect syndrome length"))
    all(x->x==0 || x==1,syndrome) || throw(ArgumentError("syndrome must be binary"))
    iseven(count(!iszero,syndrome)) || throw(ArgumentError("odd parity syndrome is impossible"))
    cfg=decoder.settings
    sites=vcat(Index[Index(4,"Logical")],[Index(2,"Detector,v=$v") for v in 1:n])
    psi=ITensorMPS.MPS(sites,fill(1,n+1))
    log_syndrome=0.
    maximum_bond=1
    maximum_discarded=0.
    sum_discarded=0.
    events=0
    frontier=0
    callback=function(;link,truncation_error)
        err=Float64(truncation_error)
        isfinite(err) && err>=0 || throw(ErrorException("invalid SVD truncation diagnostic"))
        maximum_discarded=max(maximum_discarded,err)
        sum_discarded+=err
        events+=1
    end
    for v in 1:n
        psi=_toric3d_z_tn_fuse_current(psi)
        groups=decoder.groups[v]
        support=isempty(groups) ? 1 : groups[end].neighbor-v+1
        frontier=max(frontier,isempty(groups) ? v : groups[end].neighbor)
        active_sites=frontier-v+1
        @assert 1 <= support <= active_sites <= length(psi)
        star=_toric3d_z_tn_star_mpo(psi,groups,v,syndrome[v])
        # The observed syndrome is already contracted into the entire star.
        # Only the finite star support changes exactly here. Compression still
        # reaches the full active frontier because this nonunitary update may
        # change every Schmidt spectrum before that frontier.
        psi=_toric3d_z_tn_apply_prefix(psi,star)
        psi=ITensors.replaceprime(psi,1=>0)
        ITensorMPS.truncate!(psi;site_range=1:active_sites,
            maxdim=cfg.max_bond_dimension,cutoff=cfg.cutoff,callback)
        maximum_bond=max(maximum_bond,ITensorMPS.maxlinkdim(psi))
        z=_toric3d_z_tn_mass(psi)
        isfinite(z) && z>0 || throw(DomainError(z,"nonpositive or nonfinite conditioned boundary mass at vertex $v"))
        psi[1]/=z
        log_syndrome+=log(z)
    end
    @assert length(psi)==1
    weights=vec(Array(psi[1],only(ITensorMPS.siteinds(psi))))
    all(isfinite,weights) || throw(DomainError(weights,"nonfinite logical class weights"))
    total=sum(weights)
    total>0 || throw(DomainError(total,"nonpositive final logical mass"))
    signed=weights./total
    log_syndrome+=log(total)
    valid=minimum(signed)>=-cfg.negative_tolerance
    # Preserve signed normalized weights, including roundoff-scale negatives.
    # Validity is a declared tolerance check; no clipping changes the output.
    probabilities=signed
    diagnostics=(actual_maxbond=maximum_bond,discarded_weight=sum_discarded,events=events,
        max_bond_dimension=maximum_bond,maximum_local_truncation_error=maximum_discarded,
        sum_local_truncation_error=sum_discarded,local_svd_events=events,
        minimum_probability=minimum(signed),roundoff_correction_l1=0.,
        edge_count=decoder.edge_count,detector_count=n,
        algorithm="conditioned_star_active_boundary_mps",contraction_order="plane_major_detector_order")
    (valid=valid,probabilities=probabilities,logical_mask=valid ? argmax(probabilities)-1 : -1,
     log_syndrome=log_syndrome,raw_probabilities=signed,diagnostics=diagnostics)
end
