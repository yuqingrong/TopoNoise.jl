# Build an independent unsigned stabilizer basis by forward propagation.
function h5_test_initial(c; honeycomb=false)
    n=data_qubit_count(c); m=c.Lu*c.Lv
    w=star_checks(c)[1:m]
    xs=TopoNoise._check_matrix(w,n;sparse=false)
    zs=falses(n-m,n); row=0
    for support in w, q in support[2:end]
        row+=1; zs[row,support[1]]=zs[row,q]=true
    end
    for q in 3m+1:n
        row+=1; zs[row,q]=true
    end
    if honeycomb
        seed=toric3d_encoded_seed_checks(c.encoded)
        xs=TopoNoise._check_matrix(seed.x,n;sparse=false)
        zs=TopoNoise._check_matrix(vcat(seed.z,toric3d_encoded_logical_z(c),[[q] for q in 3m+1:n]),n;sparse=false)
    end
    return hcat(vcat(xs,falses(size(zs))),vcat(falses(size(xs)),zs))
end

function h5_test_forward(state,ops)
    a=copy(state); n=size(a,2)÷2
    for op in ops
        if op.gate==:H
            q=only(op.qubits); a[:,q],a[:,n+q]=copy(a[:,n+q]),copy(a[:,q])
        else
            c,t=op.qubits; a[:,t].⊻=a[:,c];a[:,n+c].⊻=a[:,n+t]
        end
    end
    a
end
h5_test_rank(a)=length(TopoNoise._toric3d_encoded_row_basis(BitMatrix(a)).pivots)
function h5_test_contains(state,pauli,s)
    n=size(state,2)÷2; v=falses(2n);v[(pauli==:X ? 0 : n).+s].=true
    h5_test_rank(vcat(state,permutedims(v)))==h5_test_rank(state)
end


using Test, Random, TopoNoise, LinearAlgebra

# Truth tables are derived by actual 2x2/4x4 matrix conjugation, independently
# of the binary propagation formulas in the simulator.
function h5_physical_tables()
    paulis=[ComplexF64[1 0;0 1],ComplexF64[0 1;1 0],ComplexF64[1 0;0 -1],ComplexF64[0 -im;im 0]]
    h=ComplexF64[1 1;1 -1]/sqrt(2);cx=ComplexF64[1 0 0 0;0 1 0 0;0 0 0 1;0 0 1 0]
    htable=zeros(Int,4);ctable=Matrix{Tuple{Int,Int}}(undef,4,4)
    for a in 1:4
        htable[a]=only([j-1 for j in 1:4 if abs(tr(paulis[j]'*h*paulis[a]*h'))≈2])
        for b in 1:4
            out=cx*kron(paulis[a],paulis[b])*cx'
            ctable[a,b]=only([(i-1,j-1) for i in 1:4 for j in 1:4 if abs(tr(kron(paulis[i],paulis[j])'*out))≈4])
        end
    end
    htable,ctable
end
function h5_physical_segment!(labels,e,r,first,last,tables)
    ht,ct=tables
    for step in first:last
        op=e.operations[step]
        if op.gate==:H
            q=only(op.qubits);labels[q]=ht[labels[q]+1]
        else
            c,t=op.qubits;labels[c],labels[t]=ct[labels[c]+1,labels[t]+1]
        end
        for (slot,q) in enumerate(op.qubits)
            labels[q]⊻=Int(r.x[slot,step])+2Int(r.z[slot,step])
        end
    end
    labels
end

