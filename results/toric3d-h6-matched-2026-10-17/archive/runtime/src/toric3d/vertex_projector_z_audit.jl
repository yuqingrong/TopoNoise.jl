const _TORIC3D_VERTEX_Z_CLASSES=((false,false),(false,true),(true,false),(true,true))

struct Toric3DVertexZFaultSignature
    step::Int
    slot::Int
    qubit::Int
    final_layer::Bool
    syndrome_indices::Tuple{Vararg{Int}}
    upper_parity_odd::Bool
    logical_z::NTuple{2,Bool}
end

struct Toric3DVertexZObservableGroup
    syndrome_indices::Tuple{Vararg{Int}}
    upper_parity_odd::Bool
    class_counts::NTuple{4,Int}
    selected_class::NTuple{2,Bool}
    baseline_class::NTuple{2,Bool}
    unavoidable_count::Int
    baseline_failure_count::Int
end

struct Toric3DVertexZParitySummary
    mechanism_count::Int
    group_count::Int
    mixed_group_count::Int
    unavoidable_count::Int
    baseline_failure_count::Int
end

struct Toric3DVertexZAmbiguityAudit
    dimensions::NTuple{3,Int}
    protocol::Symbol
    final_layer_only::Bool
    signatures::Vector{Toric3DVertexZFaultSignature}
    groups::Vector{Toric3DVertexZObservableGroup}
    mechanism_count::Int
    group_count::Int
    mixed_group_count::Int
    unavoidable_count::Int
    baseline_failure_count::Int
    parity::NamedTuple{(:even,:odd),Tuple{Toric3DVertexZParitySummary,Toric3DVertexZParitySummary}}
end

function _vertex_z_audit_inputs(e::Toric3DVertexProjectorEncoder,
        f::Toric3DVertexDecoderFamily,protocol::Symbol)
    protocol in TORIC3D_VERTEX_PROTOCOLS || throw(ArgumentError("invalid H+5 feedback protocol"))
    _vertex_dims(e.code)==_vertex_dims(f.code) ||
        throw(DimensionMismatch("decoder and encoder dimensions differ"))
    nothing
end

"""Enumerate immutable single-Z active-output signatures before final recovery."""
function toric3d_vertex_z_fault_signatures(e::Toric3DVertexProjectorEncoder,
        f::Toric3DVertexDecoderFamily;protocol::Symbol,final_layer_only::Bool=false)
    _vertex_z_audit_inputs(e,f,protocol)
    c=e.code;m=c.Lu*c.Lv;nops=length(e.operations)
    zero=Toric3DVertexProjectorFaultRecord(_vertex_dims(c)...,
        falses(2,nops),falses(2,nops),falses(m))
    first_final=c.T==1 ? 1 : e.checkpoints[end-1]+1
    signatures=Toric3DVertexZFaultSignature[]
    logical_x=toric3d_encoded_logical_x(c)
    checks=star_checks(c)
    for (step,op) in enumerate(e.operations)
        final_layer=step>=first_final
        final_layer_only && !final_layer && continue
        for (slot,q) in enumerate(op.qubits)
            record=with_toric3d_vertex_projector_fault(zero,e,step,q,:Z)
            pre=replay_toric3d_vertex_projector_pre_final(e,record,f;
                protocol,channel=:z_only,record_history=false)
            syndrome=_vertex_syndrome(checks,pre.frame.z)
            logical_z=Tuple(_vertex_syndrome(logical_x,pre.frame.z))
            push!(signatures,Toric3DVertexZFaultSignature(step,slot,q,final_layer,
                Tuple(findall(syndrome)),pre.upper_parity_odd,logical_z))
        end
    end
    signatures
end

_vertex_z_mixed(counts::NTuple{4,Int})=count(>(0),counts)>1

function _vertex_z_parity_summary(groups::Vector{Toric3DVertexZObservableGroup},odd::Bool)
    selected=filter(g->g.upper_parity_odd==odd,groups)
    Toric3DVertexZParitySummary(
        sum(sum(g.class_counts) for g in selected;init=0),
        length(selected),
        count(g->_vertex_z_mixed(g.class_counts),selected),
        sum(g.unavoidable_count for g in selected;init=0),
        sum(g.baseline_failure_count for g in selected;init=0))
end

"""Group single-Z signatures by exactly the syndrome and observed parity available to a decoder."""
function audit_toric3d_vertex_z_boundary(e::Toric3DVertexProjectorEncoder,
        f::Toric3DVertexDecoderFamily;protocol::Symbol,final_layer_only::Bool=false)
    signatures=toric3d_vertex_z_fault_signatures(e,f;protocol,final_layer_only)
    grouped=Dict{Tuple{Tuple,Bool},Vector{Toric3DVertexZFaultSignature}}()
    for signature in signatures
        push!(get!(grouped,(signature.syndrome_indices,signature.upper_parity_odd),
            Toric3DVertexZFaultSignature[]),signature)
    end
    keys_sorted=sort!(collect(keys(grouped));by=k->(k[2],k[1]))
    groups=Toric3DVertexZObservableGroup[]
    nchecks=length(star_checks(e.code))
    logical_x=toric3d_encoded_logical_x(e.code)
    for key in keys_sorted
        members=grouped[key]
        counts=ntuple(i->count(s->s.logical_z==_TORIC3D_VERTEX_Z_CLASSES[i],members),4)
        selected_index=findmax(counts)[2]
        syndrome=falses(nchecks);syndrome[collect(key[1])].=true
        base=decode_toric3d_vertex_z_complete(f.final_z,syndrome)
        baseline_class=Tuple(_vertex_syndrome(logical_x,base.full_recovery))
        push!(groups,Toric3DVertexZObservableGroup(key[1],key[2],counts,
            _TORIC3D_VERTEX_Z_CLASSES[selected_index],baseline_class,
            length(members)-counts[selected_index],
            count(s->s.logical_z!=baseline_class,members)))
    end
    even=_vertex_z_parity_summary(groups,false)
    odd=_vertex_z_parity_summary(groups,true)
    Toric3DVertexZAmbiguityAudit(_vertex_dims(e.code),protocol,final_layer_only,
        signatures,groups,length(signatures),length(groups),
        count(g->_vertex_z_mixed(g.class_counts),groups),
        sum(g.unavoidable_count for g in groups;init=0),
        sum(g.baseline_failure_count for g in groups;init=0),
        (even=even,odd=odd))
end
