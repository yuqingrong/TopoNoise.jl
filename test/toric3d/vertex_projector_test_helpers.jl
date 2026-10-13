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

