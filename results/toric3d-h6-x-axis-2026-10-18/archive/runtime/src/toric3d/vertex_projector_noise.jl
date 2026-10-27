# An immutable byte string backs the public Bool arrays. Copies are mutable,
# but neither indexed writes nor obtaining a field can mutate a replay record.
struct Toric3DVertexFrozenBits{N} <: AbstractArray{Bool,N}
    bytes::String
    shape::NTuple{N,Int}
end
Toric3DVertexFrozenBits(a::AbstractArray{Bool,N}) where N =
    Toric3DVertexFrozenBits{N}(String(UInt8.(vec(a))),size(a))
Base.size(a::Toric3DVertexFrozenBits)=a.shape
Base.IndexStyle(::Type{<:Toric3DVertexFrozenBits})=IndexLinear()
function Base.getindex(a::Toric3DVertexFrozenBits,i::Int)
    @boundscheck checkbounds(a,i)
    @inbounds codeunit(a.bytes,i)==0x01
end
Base.copy(a::Toric3DVertexFrozenBits)=BitArray(reshape([x for x in a],size(a)))

"""Immutable H+5 active-output faults plus a paired, hidden clean upper branch."""
struct Toric3DVertexProjectorFaultRecord
    T::Int
    Lu::Int
    Lv::Int
    x::Toric3DVertexFrozenBits{2}
    z::Toric3DVertexFrozenBits{2}
    upper_branch::Toric3DVertexFrozenBits{1}
    function Toric3DVertexProjectorFaultRecord(T::Integer,Lu::Integer,Lv::Integer,
            x::AbstractMatrix{Bool},z::AbstractMatrix{Bool},b::AbstractVector{Bool})
        c=Toric3DVertexProjectorCode(T,Lu,Lv);m=c.Lu*c.Lv
        size(x)==size(z)==(2,6c.T*m) || throw(DimensionMismatch("fault record shape differs from H+5 schedule"))
        length(b)==m || throw(DimensionMismatch("upper branch size differs"))
        iseven(count(b)) || throw(ArgumentError("clean upper branch must be even"))
        any(x[2,1:6:end]) || any(z[2,1:6:end]) ?
            throw(ArgumentError("H has no second output")) : nothing
        new(c.T,c.Lu,c.Lv,Toric3DVertexFrozenBits(x),Toric3DVertexFrozenBits(z),Toric3DVertexFrozenBits(b))
    end
end
Base.:(==)(a::Toric3DVertexProjectorFaultRecord,b::Toric3DVertexProjectorFaultRecord)=
    (a.T,a.Lu,a.Lv)==(b.T,b.Lu,b.Lv) && a.x==b.x && a.z==b.z && a.upper_branch==b.upper_branch

"""Draw X then Z at each active output, in gate/slot order; draw branch last."""
function sample_toric3d_vertex_projector_record(rng::Random.AbstractRNG,
        e::Toric3DVertexProjectorEncoder;p_x::Real=0,p_z::Real=0)
    all(p->isfinite(p) && 0<=p<=1,(p_x,p_z)) ||
        throw(ArgumentError("H+5 fault probabilities must be finite and in [0,1]"))
    x=falses(2,length(e.operations));z=copy(x)
    for (step,op) in enumerate(e.operations),slot in eachindex(op.qubits)
        x[slot,step]=rand(rng)<p_x;z[slot,step]=rand(rng)<p_z
    end
    # Sampling needs only branch size, not a rebuilt correction tree.
    b=falses(e.code.Lu*e.code.Lv)
    for i in 1:length(b)-1;b[i]=rand(rng,Bool);end
    b[end]=isodd(count(b))
    Toric3DVertexProjectorFaultRecord(_vertex_dims(e.code)...,x,z,b)
end

function _vertex_validate_record(e,r)
    (r.T,r.Lu,r.Lv)==_vertex_dims(e.code) || throw(DimensionMismatch("encoder and record dimensions differ"))
    size(r.x,2)==length(e.operations) || throw(DimensionMismatch("encoder and record operation count differ"))
    nothing
end

function with_toric3d_vertex_projector_fault(r::Toric3DVertexProjectorFaultRecord,
        e::Toric3DVertexProjectorEncoder,step::Integer,q::Integer,pauli::Symbol)
    _vertex_validate_record(e,r)
    pauli in (:X,:Z,:Y) || throw(ArgumentError("fault Pauli must be X, Z, or Y"))
    1<=step<=length(e.operations) || throw(ArgumentError("fault step is out of range"))
    slot=findfirst(==(q),e.operations[step].qubits)
    isnothing(slot) && throw(ArgumentError("fault qubit is not an active output"))
    x=copy(r.x);z=copy(r.z)
    pauli in (:X,:Y) && (x[slot,step]⊻=true)
    pauli in (:Z,:Y) && (z[slot,step]⊻=true)
    Toric3DVertexProjectorFaultRecord(r.T,r.Lu,r.Lv,x,z,r.upper_branch)
end

function propagate_toric3d_vertex_projector_segment!(frame::PauliFrame,
        e::Toric3DVertexProjectorEncoder,r::Toric3DVertexProjectorFaultRecord,
        first::Integer,last::Integer)
    _vertex_validate_record(e,r)
    length(frame.x)==data_qubit_count(e.code) || throw(DimensionMismatch("frame size differs"))
    1<=first<=length(e.operations)+1 && 0<=last<=length(e.operations) && first<=last+1 ||
        throw(ArgumentError("operation segment is out of range"))
    for step in first:last
        op=e.operations[step]
        _propagate_frame_operation!(frame.x,frame.z,op)
        for (slot,q) in enumerate(op.qubits)
            frame.x[q]⊻=r.x[slot,step];frame.z[q]⊻=r.z[slot,step]
        end
    end
    frame
end
