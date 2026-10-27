struct Toric3DVertexCircuitZDecoder
    code::Toric3DVertexProjectorCode
    protocol::Symbol
    base::Toric3DVertexMatchingDecoder
    groups::Dict{Tuple{Tuple,Bool},Toric3DVertexZObservableGroup}
    logical_adjustments::Dict{NTuple{2,Bool},BitVector}
end

struct Toric3DVertexCircuitZDecision
    base::Toric3DVertexMatchingDecision
    full_recovery::BitVector
    selected_class::NTuple{2,Bool}
    base_class::NTuple{2,Bool}
    class_counts::NTuple{4,Int}
    mixed::Bool
    adjusted::Bool
    fallback::Bool
end

struct Toric3DVertexZControlOutcome
    syndrome::BitVector
    upper_parity_odd::Bool
    decision::Toric3DVertexCircuitZDecision
    baseline_residual::BitVector
    control_residual::BitVector
    residual_syndrome::BitVector
    baseline_logical_z::NTuple{2,Bool}
    control_logical_z::NTuple{2,Bool}
    baseline_logical_z_failure::Bool
    control_logical_z_failure::Bool
end

struct Toric3DVertexParityErasureEstimate
    shots::Int
    accepted::Int
    rejected::Int
    failures::Int
    accepted_failures::Int
    acceptance_rate::Float64
    acceptance_ci::Tuple{Float64,Float64}
    conditional_rate::Union{Nothing,Float64}
    conditional_ci::Union{Nothing,Tuple{Float64,Float64}}
end

_vertex_z_class_xor(a::NTuple{2,Bool},b::NTuple{2,Bool})=
    (xor(a[1],b[1]),xor(a[2],b[2]))

function _vertex_z_logical_adjustments(c::Toric3DVertexProjectorCode)
    n=data_qubit_count(c);checks=star_checks(c)
    logical_x=toric3d_encoded_logical_x(c)
    logical_z=toric3d_encoded_logical_z(c)
    adjustments=Dict{NTuple{2,Bool},BitVector}()
    for mask in 0:3
        adjustment=falses(n)
        for i in 1:2
            !iszero(mask&(1<<(i-1))) && (adjustment[logical_z[i]].⊻=true)
        end
        !any(_vertex_syndrome(checks,adjustment)) ||
            error("declared logical-Z adjustment has nonzero star syndrome")
        cls=Tuple(_vertex_syndrome(logical_x,adjustment))
        haskey(adjustments,cls) &&
            error("declared logical-Z representatives do not span four unique classes")
        adjustments[cls]=adjustment
    end
    length(adjustments)==4 ||
        error("declared logical-Z representatives do not span four unique classes")
    adjustments
end

"""Layer a decoder-visible single-fault class table over the existing final MWPM decoder."""
function build_toric3d_vertex_circuit_z_decoder(f::Toric3DVertexDecoderFamily,
        audit::Toric3DVertexZAmbiguityAudit)
    audit.protocol in TORIC3D_VERTEX_PROTOCOLS ||
        throw(ArgumentError("audit protocol is not canonical"))
    audit.dimensions==_vertex_dims(f.code) ||
        throw(DimensionMismatch("audit and decoder dimensions differ"))
    groups=Dict{Tuple{Tuple,Bool},Toric3DVertexZObservableGroup}()
    for group in audit.groups
        key=(group.syndrome_indices,group.upper_parity_odd)
        haskey(groups,key) && throw(ArgumentError("audit contains duplicate observable groups"))
        all(>=(0),group.class_counts) && sum(group.class_counts)>0 ||
            throw(ArgumentError("audit group has invalid class counts"))
        selected=findfirst(==(maximum(group.class_counts)),group.class_counts)
        group.selected_class==_TORIC3D_VERTEX_Z_CLASSES[selected] ||
            throw(ArgumentError("audit group selected class violates deterministic tie break"))
        groups[key]=group
    end
    Toric3DVertexCircuitZDecoder(f.code,audit.protocol,f.final_z,groups,
        _vertex_z_logical_adjustments(f.code))
end

"""Decode a final Z syndrome, adjusting only the logical class for known audit keys."""
function decode_toric3d_vertex_circuit_z(d::Toric3DVertexCircuitZDecoder,
        syndrome::AbstractVector{Bool},upper_parity_odd::Bool)
    base=decode_toric3d_vertex_z_complete(d.base,syndrome)
    logical_x=toric3d_encoded_logical_x(d.code)
    base_class=Tuple(_vertex_syndrome(logical_x,base.full_recovery))
    key=(Tuple(findall(syndrome)),upper_parity_odd)
    if !haskey(d.groups,key)
        return Toric3DVertexCircuitZDecision(base,copy(base.full_recovery),base_class,
            base_class,(0,0,0,0),false,false,true)
    end
    group=d.groups[key]
    difference=_vertex_z_class_xor(base_class,group.selected_class)
    full=base.full_recovery .⊻ d.logical_adjustments[difference]
    _vertex_syndrome(star_checks(d.code),full)==syndrome ||
        error("circuit-aware logical adjustment changed the requested syndrome")
    Toric3DVertexCircuitZDecision(base,BitVector(full),group.selected_class,
        base_class,group.class_counts,_vertex_z_mixed(group.class_counts),
        any(difference),false)
end

"""Apply baseline and circuit-aware final recoveries to one shared pre-final frame."""
function finalize_toric3d_vertex_z_control(pre::Toric3DVertexPreFinalState,
        d::Toric3DVertexCircuitZDecoder)
    pre.protocol==d.protocol || throw(ArgumentError("pre-final state and control protocol differ"))
    length(pre.frame.z)==data_qubit_count(d.code) ||
        throw(DimensionMismatch("pre-final frame and control dimensions differ"))
    syndrome=_vertex_syndrome(star_checks(d.code),pre.frame.z)
    decision=decode_toric3d_vertex_circuit_z(d,syndrome,pre.upper_parity_odd)
    baseline=pre.frame.z .⊻ decision.base.full_recovery
    control=pre.frame.z .⊻ decision.full_recovery
    residual_syndrome=_vertex_syndrome(star_checks(d.code),control)
    !any(_vertex_syndrome(star_checks(d.code),baseline)) ||
        error("baseline final Z recovery left syndrome")
    !any(residual_syndrome) || error("circuit-aware final Z recovery left syndrome")
    logical_x=toric3d_encoded_logical_x(d.code)
    baseline_logical=Tuple(_vertex_syndrome(logical_x,baseline))
    control_logical=Tuple(_vertex_syndrome(logical_x,control))
    Toric3DVertexZControlOutcome(BitVector(syndrome),pre.upper_parity_odd,decision,
        BitVector(baseline),BitVector(control),BitVector(residual_syndrome),
        baseline_logical,control_logical,any(baseline_logical),any(control_logical))
end

"""Report the cost and conditional failure rate of accepting only even upper parity."""
function toric3d_vertex_parity_erasure_estimate(upper_parity_odd::AbstractVector{Bool},
        logical_failures::AbstractVector{Bool})
    length(upper_parity_odd)==length(logical_failures) ||
        throw(DimensionMismatch("parity and failure vectors differ in length"))
    shots=length(upper_parity_odd)
    shots>0 || throw(ArgumentError("parity-erasure estimate needs at least one shot"))
    accepted=count(!,upper_parity_odd);rejected=shots-accepted
    failures=count(identity,logical_failures)
    accepted_failures=count(i->!upper_parity_odd[i]&&logical_failures[i],eachindex(upper_parity_odd))
    conditional_rate=accepted==0 ? nothing : accepted_failures/accepted
    conditional_ci=accepted==0 ? nothing : _toric3d_wilson(accepted_failures,accepted)
    Toric3DVertexParityErasureEstimate(shots,accepted,rejected,failures,
        accepted_failures,accepted/shots,_toric3d_wilson(accepted,shots),
        conditional_rate,conditional_ci)
end
