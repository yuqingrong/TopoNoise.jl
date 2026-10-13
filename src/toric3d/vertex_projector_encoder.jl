"""Direct vertex projectors, one H and five CNOTs per diagonal cell."""
struct Toric3DVertexProjectorEncoder
    code::Toric3DVertexProjectorCode
    operations::Vector{EncoderOperation}
    checkpoints::Vector{Int}
end

function toric3d_vertex_projector_encoder(code::Toric3DVertexProjectorCode)
    operations=EncoderOperation[]; checkpoints=Int[]
    q(d,t,u,v)=data_qubit_index(code,d,t,mod(u,code.Lu),mod(v,code.Lv))
    for t in 1:code.T
        for u in 0:code.Lu-1,v in 0:code.Lv-1
            pivot=q(1,t,u+1,v)
            push!(operations,EncoderOperation(:H,[pivot]))
            for target in (q(2,t,u,v+1),q(3,t,u,v),q(1,t-1,u,v),q(2,t-1,u,v),q(3,t-1,u,v))
                push!(operations,EncoderOperation(:CNOT,[pivot,target]))
            end
        end
        push!(checkpoints,length(operations))
    end
    Toric3DVertexProjectorEncoder(code,operations,checkpoints)
end

# The lower state is a product of GHZ triples on the disjoint white stars
# and fresh zeros. Check signed membership after pulling an operator back.
function _vertex_prepared(code,operations,pauli,support)
    n=data_qubit_count(code);m=code.Lu*code.Lv
    x=falses(n);z=falses(n);negative=Ref(false)
    (pauli===:X ? x : z)[support].=true
    for op in Iterators.reverse(operations)
        _conjugate_pauli!(x,z,negative,op)
    end
    any(@view x[3m+1:end]) && return false
    for w in star_checks(code)[1:m]
        all(x[q]==x[w[1]] for q in w) || return false
        iseven(count(z[w])) || return false
    end
    overlap=count(x .& z)
    iseven(overlap) && negative[]==isodd(overlap÷2)
end

"""Signed stabilizer-algebra verification of every completed H+5 prefix."""
function verify_toric3d_vertex_projector_encoder(e::Toric3DVertexProjectorEncoder)::Bool
    c=e.code;m=c.Lu*c.Lv
    expected=toric3d_vertex_projector_encoder(c)
    e.checkpoints==expected.checkpoints || return false
    length(e.operations)==length(expected.operations) || return false
    all(a.gate==b.gate && a.qubits==b.qubits for (a,b) in zip(e.operations,expected.operations)) || return false
    for t in 1:c.T
        ops=@view e.operations[1:e.checkpoints[t]]
        for (pauli,checks) in ((:X,star_checks(c)[1:(t+1)*m]),
                               (:Z,plaquette_checks(c)[1:3m*t]),
                               (:Z,toric3d_encoded_logical_z(c)))
            all(s->_vertex_prepared(c,ops,pauli,s),checks) || return false
        end
        for u in 0:c.Lu-1,v in 0:c.Lv-1
            s=[data_qubit_index(c,d,t,u,v) for d in 1:3]
            _vertex_prepared(c,ops,:X,s) && return false
        end
    end
    true
end
