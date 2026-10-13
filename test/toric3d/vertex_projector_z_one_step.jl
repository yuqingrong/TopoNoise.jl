using Test,Random,TopoNoise
# Independent exhaustive data-subset oracle, with at most one outgoing edge.
function h5_strip_oracle(h,m)
    masks=[sum(1<<(i-1) for i in findall(h[:,j])) for j in axes(h,2)]
    costs=fill(typemax(Int),1<<size(h,1))
    for bits in 0:(1<<(3m))-1
        syndrome=0
        for j in 1:3m;isodd(bits>>(j-1)) && (syndrome⊻=masks[j]);end
        w=count_ones(bits);costs[syndrome+1]=min(costs[syndrome+1],w)
        for j in 3m+1:6m
            s=syndrome⊻masks[j];costs[s+1]=min(costs[s+1],w+1)
        end
    end
    costs
end
@testset "H+5 one-step Z matching uses exactly one temporal mate for odd parity" begin
    @test isdefined(TopoNoise,:build_toric3d_vertex_z_one_step_decoder)
    if isdefined(TopoNoise,:build_toric3d_vertex_z_one_step_decoder)
        d=build_toric3d_vertex_z_one_step_decoder(Toric3DVertexProjectorCode(2,2,2),1)
        s=falses(8);s[1]=true;a=decode_toric3d_vertex_z_one_step(d,s)
        @test a.signature==[(1,9,(3,13))]
        @test count(a.full_recovery)==2
        @test findall(a.route)==[3]
        @test findall(a.residual_syndrome)==[7]
        s[2]=true;a=decode_toric3d_vertex_z_one_step(d,s)
        @test a.signature==[(1,2,(6,10))]
    end
end
if isdefined(TopoNoise,:build_toric3d_vertex_z_one_step_decoder)
@testset "H+5 period-two Z decoding exhausts every reachable strip syndrome" begin
    d=build_toric3d_vertex_z_one_step_decoder(Toric3DVertexProjectorCode(2,2,2),1)
    costs=h5_strip_oracle(d.window.matrix,4)
    for mask in 0:255
        s=BitVector([isodd(mask>>k) for k in 0:7]);a=decode_toric3d_vertex_z_one_step(d,s)
        @test TopoNoise._vertex_syndrome(d.window.matrix,a.full_recovery)==s
        @test count(a.full_recovery)==costs[mask+1]
        @test sum(length(p.edges) for p in a.paths;init=0)==costs[mask+1]
        @test count(p->p.second==9,a.paths)==Int(isodd(count(s)))
        @test a.residual_syndrome==s.⊻TopoNoise._vertex_syndrome(d.window.matrix,a.route)
        @test !any(a.residual_syndrome[1:4])
        @test !any(a.route[13:end])
    end
end
@testset "H+5 causal Z trimming clears retiring sheets and confines residuals" begin
    rng=MersenneTwister(72)
    for dims in ((3,2,3),(3,3,2),(3,3,3)),t in 1:2
        c=Toric3DVertexProjectorCode(dims...);d=build_toric3d_vertex_z_one_step_decoder(c,t)
        for _ in 1:12
            s=bitrand(rng,2c.Lu*c.Lv);a=decode_toric3d_vertex_z_one_step(d,s)
            @test !any(a.residual_syndrome[d.window.retired_rows])
            @test all(q->data_qubit_coordinate(c,q)[2]==t-1,d.window.qubits[findall(a.route)])
            @test TopoNoise._vertex_syndrome(d.window.matrix,a.full_recovery)==s
        end
    end
end
@testset "Symmetric Z endpoint trimming is causally invalid" begin
    d=build_toric3d_vertex_z_one_step_decoder(Toric3DVertexProjectorCode(3,2,2),1)
    s=falses(8);s[[1,8]].=true;a=decode_toric3d_vertex_z_one_step(d,s)
    path=only(a.paths).edges
    @test length(path)==3
    symmetric=falses(24);symmetric[path[[1,end]]].=true
    post=s.⊻TopoNoise._vertex_syndrome(d.window.matrix,symmetric)
    @test any(post[1:4])
    @test !any(a.residual_syndrome[1:4])
end
end
