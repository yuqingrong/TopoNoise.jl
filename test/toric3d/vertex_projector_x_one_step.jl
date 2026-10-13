using Test,Random,TopoNoise
@testset "H+5 one-step X MWPM returns canonical physical paths" begin
    @test isdefined(TopoNoise,:build_toric3d_vertex_x_one_step_decoder)
    if isdefined(TopoNoise,:build_toric3d_vertex_x_one_step_decoder)
        c=Toric3DVertexProjectorCode(2,2,2);d=build_toric3d_vertex_x_one_step_decoder(c,1)
        a=decode_toric3d_vertex_x_one_step(d,trues(4))
        @test a.signature==[(1,2,(13,)),(3,4,(15,))]
        @test d.version=="2.4.0"
        @test a.backend=="pymatching"
    end
end
if isdefined(TopoNoise,:build_toric3d_vertex_x_one_step_decoder)
@testset "H+5 period-two X decoding exhausts every reachable syndrome" begin
    d=build_toric3d_vertex_x_one_step_decoder(Toric3DVertexProjectorCode(2,2,2),1)
    for mask in 0:15
        s=BitVector([isodd(mask>>k) for k in 0:3])
        if isodd(count(s));@test_throws ArgumentError decode_toric3d_vertex_x_one_step(d,s);continue;end
        a=decode_toric3d_vertex_x_one_step(d,s)
        @test TopoNoise._vertex_syndrome(d.window.matrix,a.full_recovery)==s
        @test a.residual_syndrome==s.⊻TopoNoise._vertex_syndrome(d.window.matrix,a.route)
        @test count(a.full_recovery)==count(s)÷2
        @test !any(a.residual_syndrome)
    end
end
@testset "H+5 one-step X routes predict their residual syndrome" begin
    d=build_toric3d_vertex_x_one_step_decoder(Toric3DVertexProjectorCode(2,5,5),1)
    s=falses(25);s[[1,13]].=true;a=decode_toric3d_vertex_x_one_step(d,s)
    @test length(only(a.paths).edges)>2
    @test count(a.route)==2
    @test count(a.residual_syndrome)==2
    @test a.residual_syndrome==s.⊻TopoNoise._vertex_syndrome(d.window.matrix,a.route)
end
@testset "H+5 old-target X faults remain terminal local noise" begin
    e=toric3d_vertex_projector_encoder(Toric3DVertexProjectorCode(3,2,2));c=e.code
    zero=sample_toric3d_vertex_projector_record(MersenneTwister(1),e)
    for t in 1:2
        w=toric3d_vertex_x_slice_window(c,t)
        for q in toric3d_slice_qubits(c,t-1)
            frame=PauliFrame(falses(data_qubit_count(c)),falses(data_qubit_count(c)));frame.x[q]=true
            @test !any(TopoNoise._vertex_syndrome(w.checks,frame.x))
            @test all(q ∉ op.qubits for op in e.operations[e.checkpoints[t]+1:end])
            propagate_toric3d_vertex_projector_segment!(frame,e,zero,e.checkpoints[t]+1,length(e.operations))
            @test findall(frame.x)==[q]
        end
    end
end
end
