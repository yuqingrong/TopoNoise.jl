using Test, TopoNoise, LinearAlgebra

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

@testset "H+5 target geometry and logical algebra" begin
    @test isdefined(TopoNoise,:Toric3DVertexProjectorCode)
    if isdefined(TopoNoise,:Toric3DVertexProjectorCode)
        for dims in ((1,2,2),(2,2,2),(2,2,3),(2,3,2),(3,3,3))
            c=Toric3DVertexProjectorCode(dims...); t,u,v=dims;n=data_qubit_count(c)
            hx=star_check_matrix(c);hz=plaquette_check_matrix(c)
            @test n==3(t+1)*u*v
            @test n-h5_test_rank(hx)-h5_test_rank(hz)==2
            @test all(iseven,Int.(hx)*transpose(Int.(hz)))
            lx=TopoNoise._check_matrix(toric3d_encoded_logical_x(c),n;sparse=false)
            lz=TopoNoise._check_matrix(toric3d_encoded_logical_z(c),n;sparse=false)
            @test Int.(lx)*transpose(Int.(lz)).%2 == [1 0;0 1]
            @test minimum(count.(eachrow(lz)))==2min(u,v)
            @test minimum(count.(eachrow(lx)))==(t+1)*min(u,v)
        end
        @test_throws ArgumentError Toric3DVertexProjectorCode(0,2,2)
        @test_throws ArgumentError Toric3DVertexProjectorCode(2,1,2)
    end
end

if isdefined(TopoNoise,:Toric3DVertexProjectorCode)
@testset "H+5 direct-projector cells use the exact six-operation schedule" begin
    e=toric3d_vertex_projector_encoder(Toric3DVertexProjectorCode(2,2,2))
    @test length(e.operations)==48
    @test sum(length(o.qubits) for o in e.operations)==88
    @test e.checkpoints==[24,48]
    @test [o.qubits for o in e.operations[1:6]]==[[15],[15,18],[15,21],[15,1],[15,5],[15,9]]
    @test [o.gate for o in e.operations]==repeat([:H,fill(:CNOT,5)...],8)
    bad=copy(e.operations);bad[6]=TopoNoise.EncoderOperation(:CNOT,[15,10])
    @test !verify_toric3d_vertex_projector_encoder(Toric3DVertexProjectorEncoder(e.code,bad,e.checkpoints))
end
@testset "H+5 layers are fresh, disjoint, ordered, and causal" begin
    c=Toric3DVertexProjectorCode(3,2,3);e=toric3d_vertex_projector_encoder(c);m=c.Lu*c.Lv
    seen=Set{Int}()
    for t in 1:c.T
        layer=Set{Int}()
        for u in 0:c.Lu-1,v in 0:c.Lv-1
            j=6*((t-1)*m+u*c.Lv+v)
            block=e.operations[j+1:j+6];support=unique(vcat([o.qubits for o in block]...))
            @test isempty(intersect(layer,support));union!(layer,support)
            @test only(block[1].qubits) ∉ seen
            @test all(q->data_qubit_coordinate(c,q)[2]==t,block[2].qubits)
            @test all(q->data_qubit_coordinate(c,q)[2]==t-1,[block[k].qubits[2] for k in 4:6])
            union!(seen,support)
        end
    end
end
@testset "H+5 clean prefixes prepare exactly the completed checks" begin
    for dims in ((1,2,2),(2,2,3),(3,3,2))
        c=Toric3DVertexProjectorCode(dims...);e=toric3d_vertex_projector_encoder(c);m=c.Lu*c.Lv
        @test verify_toric3d_vertex_projector_encoder(e)
        for t in 1:c.T
            state=h5_test_forward(h5_test_initial(c),e.operations[1:e.checkpoints[t]])
            @test h5_test_rank(state)==data_qubit_count(c)
            @test all(s->h5_test_contains(state,:X,s),star_checks(c)[1:(t+1)*m])
            @test all(s->h5_test_contains(state,:Z,s),plaquette_checks(c)[1:3m*t])
            @test all(s->h5_test_contains(state,:Z,s),toric3d_encoded_logical_z(c))
            frontier=[[data_qubit_index(c,d,t,u,v) for d in 1:3] for u in 0:c.Lu-1 for v in 0:c.Lv-1]
            @test all(s->!h5_test_contains(state,:X,s),frontier)
        end
    end
end
@testset "The H+6 honeycomb seed is invalid for H+5" begin
    c=Toric3DVertexProjectorCode(1,2,2);e=toric3d_vertex_projector_encoder(c)
    state=h5_test_forward(h5_test_initial(c;honeycomb=true),e.operations)
    @test !all(s->h5_test_contains(state,:Z,s),plaquette_checks(c))
end
end
