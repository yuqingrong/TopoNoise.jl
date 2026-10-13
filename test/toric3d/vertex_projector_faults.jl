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

@testset "H+5 records enumerate exactly eleven outputs per vertex" begin
    @test isdefined(TopoNoise,:sample_toric3d_vertex_projector_record)
    if isdefined(TopoNoise,:sample_toric3d_vertex_projector_record)
        e=toric3d_vertex_projector_encoder(Toric3DVertexProjectorCode(2,2,2))
        for (px,pz) in ((0.,0.),(1.,0.),(0.,1.),(1.,1.))
            r=sample_toric3d_vertex_projector_record(MersenneTwister(5),e;p_x=px,p_z=pz)
            @test size(r.x)==size(r.z)==(2,48)
            @test count(r.x)==88px
            @test count(r.z)==88pz
            @test all(!r.x[2,i]&&!r.z[2,i] for i in 1:6:48)
            @test iseven(count(r.upper_branch))
        end
        @test_throws ArgumentError sample_toric3d_vertex_projector_record(MersenneTwister(5),e;p_x=NaN)
        @test_throws ArgumentError sample_toric3d_vertex_projector_record(MersenneTwister(5),e;p_z=-.1)
    end
end
if isdefined(TopoNoise,:sample_toric3d_vertex_projector_record)
@testset "H+5 frame propagation equals a physical Pauli reference" begin
    e=toric3d_vertex_projector_encoder(Toric3DVertexProjectorCode(2,2,2));n=data_qubit_count(e.code)
    zero=sample_toric3d_vertex_projector_record(MersenneTwister(2),e;p_x=0,p_z=0);tables=h5_physical_tables()
    for (step,op) in enumerate(e.operations),q in op.qubits,pauli in (:X,:Z)
        r=with_toric3d_vertex_projector_fault(zero,e,step,q,pauli)
        frame=PauliFrame(falses(n),falses(n));labels=zeros(Int,n)
        propagate_toric3d_vertex_projector_segment!(frame,e,r,1,24)
        h5_physical_segment!(labels,e,r,1,24,tables)
        @test Int.(frame.x).+2Int.(frame.z)==labels
        # A physical correction and a frame correction at the same checkpoint.
        frame.x[[2,7,13]].⊻=true;frame.z[[1,10,19]].⊻=true
        labels[[2,7,13]].⊻=1;labels[[1,10,19]].⊻=2
        propagate_toric3d_vertex_projector_segment!(frame,e,r,25,48)
        h5_physical_segment!(labels,e,r,25,48,tables)
        @test Int.(frame.x).+2Int.(frame.z)==labels
        @test count(r.x)+count(r.z)==1
        @test !any(zero.x)&&!any(zero.z)
    end
    @test_throws ArgumentError with_toric3d_vertex_projector_fault(zero,e,1,1,:X)
    @test_throws ArgumentError with_toric3d_vertex_projector_fault(zero,e,1,15,:I)
    @test_throws ArgumentError propagate_toric3d_vertex_projector_segment!(PauliFrame(falses(n),falses(n)),e,zero,0,1)
end
@testset "H+5 replay records are immutable and channel separated" begin
    e=toric3d_vertex_projector_encoder(Toric3DVertexProjectorCode(2,2,2))
    r=sample_toric3d_vertex_projector_record(MersenneTwister(1234),e;p_x=.2,p_z=.3)
    @test r==sample_toric3d_vertex_projector_record(MersenneTwister(1234),e;p_x=.2,p_z=.3)
    @test_throws Exception setindex!(r.x,true,1,1)
    @test_throws Exception setindex!(r.z,true,1,1)
    @test_throws Exception setindex!(r.upper_branch,true,1)
    x=copy(r.x);z=copy(r.z);b=copy(r.upper_branch)
    rr=Toric3DVertexProjectorFaultRecord(2,2,2,x,z,b)
    x[1,1]⊻=true;b[1]⊻=true
    @test rr==r
    @test r.x!=r.z
    f=PauliFrame(falses(data_qubit_count(e.code)),falses(data_qubit_count(e.code)))
    propagate_toric3d_vertex_projector_segment!(f,e,r,1,48)
    @test r==sample_toric3d_vertex_projector_record(MersenneTwister(1234),e;p_x=.2,p_z=.3)
    bad=copy(r.x);bad[2,1]=true
    @test_throws ArgumentError Toric3DVertexProjectorFaultRecord(2,2,2,bad,r.z,r.upper_branch)
end
end
