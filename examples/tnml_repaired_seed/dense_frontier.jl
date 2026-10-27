# Independent reference copied without algorithm changes on 2026-10-27.
# Source: examples/tndecoder3d_trial/dense_frontier.jl in the
# toric3d-tnml-z research worktree (copied locally for portability).
# Original SHA256: b0b89152a8e131fcc57ac4429444427f6d3e034e2933edd6c1e7a82f502fdf00
# New experiments pass newly built seed/transfer topologies; no old posterior data are reused.

"""
Independent untruncated dense detector-frontier variable elimination.
This enumerates every boundary assignment and four joint logical classes.
There is no top-K pruning, SVD, or approximate contraction. Memory is
4*2^(maximum open detector count) Float64s per boundary, rather than the
full joint distribution over all detectors. The max_width guard is strict.
"""
module TNDecoder3DDenseFrontier
using TopoNoise
export dense_exact_probabilities

function dense_exact_probabilities(topology, syndrome::AbstractVector{Bool}, p::Real; max_width::Int=22)
    V=size(topology.raw_detector_matrix,1)
    length(syndrome)==V || throw(DimensionMismatch("incorrect detector count"))
    qs=TopoNoise.toric3d_z_tn_edge_probabilities(topology,p)
    active=Int[]
    state=zeros(4,1); state[1]=1
    log_marginal=0.0
    maximum_width=0
    for v in 1:V
        outgoing=findall(e -> e[1]==v, topology.detector_endpoints)
        nextactive=sort!(unique(vcat(filter(!=(v),active),[topology.detector_endpoints[j][2] for j in outgoing])))
        width=length(nextactive)
        width <= max_width || throw(ArgumentError("dense exact frontier needs width $width > guard $max_width"))
        maximum_width=max(maximum_width,width)
        lookup=Dict(node=>i for (i,node) in enumerate(nextactive))
        # Cache old-label reindexing, using least-significant-set-bit recurrence.
        mapped=zeros(Int,1 << length(active))
        bitmaps=[node==v ? 0 : 1 << (lookup[node]-1) for node in active]
        for label in 1:length(mapped)-1
            bit=trailing_zeros(label)+1
            mapped[label+1]=xor(mapped[xor(label,1 << (bit-1))+1],bitmaps[bit])
        end
        # Exact independent edge subsets leaving v, merged only if they have
        # identical boundary, local detector and joint-logical effects.
        transitions=Dict{Tuple{Bool,Int,Int},Float64}((false,0,0)=>1.0)
        for j in outgoing
            q=qs[j]
            boundary_delta=1 << (lookup[topology.detector_endpoints[j][2]]-1)
            logical_delta=Int(topology.logical_masks[j])
            updated=Dict{Tuple{Bool,Int,Int},Float64}()
            for ((parity,mask,logical),weight) in transitions
                absent=(parity,mask,logical)
                present=(!parity,xor(mask,boundary_delta),xor(logical,logical_delta))
                updated[absent]=get(updated,absent,0.0)+(1-q)*weight
                updated[present]=get(updated,present,0.0)+q*weight
            end
            transitions=updated
        end
        choices=([entry for entry in transitions if !entry.first[1]],
                 [entry for entry in transitions if entry.first[1]])
        nextstate=zeros(4,1 << width)
        vposition=findfirst(==(v),active)
        vbit=vposition===nothing ? 0 : 1 << (vposition-1)
        for label in 0:size(state,2)-1
            @inbounds a,b,c,d=state[1,label+1],state[2,label+1],state[3,label+1],state[4,label+1]
            iszero(a+b+c+d) && continue
            parity=xor(!iszero(label & vbit),syndrome[v])
            for ((_,delta,logical),weight) in choices[Int(parity)+1]
                col=xor(mapped[label+1],delta)+1
                @inbounds begin
                    nextstate[xor(0,logical)+1,col] += weight*a
                    nextstate[xor(1,logical)+1,col] += weight*b
                    nextstate[xor(2,logical)+1,col] += weight*c
                    nextstate[xor(3,logical)+1,col] += weight*d
                end
            end
        end
        normalization=sum(nextstate)
        normalization>0 || throw(ArgumentError("syndrome has zero probability"))
        nextstate ./= normalization
        log_marginal+=log(normalization)
        active=nextactive
        state=nextstate
    end
    @assert isempty(active) && size(state)==(4,1)
    return (probabilities=vec(state),log_marginal_probability=log_marginal,
            maximum_width=maximum_width,exact=true)
end
end
