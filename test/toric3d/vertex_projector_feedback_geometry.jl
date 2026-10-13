using Test, Random, TopoNoise
isdefined(@__MODULE__,:h5_test_initial) || include("vertex_projector_test_helpers.jl")

@testset "H+5 Q slices are target stabilizers with even syndrome" begin
    @test isdefined(TopoNoise,:toric3d_vertex_x_slice_window)
    if isdefined(TopoNoise,:toric3d_vertex_x_slice_window)
        for dims in ((2,2,2),(2,2,3),(2,3,2),(3,3,3))
            c=Toric3DVertexProjectorCode(dims...);e=toric3d_vertex_projector_encoder(c)
            for t in 1:c.T-1
                w=toric3d_vertex_x_slice_window(c,t)
                @test all(==(2),vec(sum(w.matrix;dims=1)))
                @test w.qubits==toric3d_slice_qubits(c,t)
                state=h5_test_forward(h5_test_initial(c),e.operations[1:e.checkpoints[t]])
                @test all(s->h5_test_contains(state,:Z,s),w.checks)
                for q in axes(w.matrix,2);@test iseven(count(w.matrix[:,q]));end
            end
        end
    end
end
if isdefined(TopoNoise,:toric3d_vertex_x_slice_window)
@testset "H+5 X and Z windows use prepared causal supports" begin
    c=Toric3DVertexProjectorCode(3,2,3);m=6
    for t in 1:2
        x=toric3d_vertex_x_complete_window(c,t);z=toric3d_vertex_z_complete_window(c,t)
        d=vcat(toric3d_slice_qubits(c,t-1),toric3d_slice_qubits(c,t))
        @test x.qubits==z.qubits==d
        @test x.checks==plaquette_checks(c)[max(1,3m*(t-2)+1):3m*t]
        @test z.checks==star_checks(c)[(t-1)*m+1:(t+1)*m]
        @test z.retired_rows==collect(1:m)
        @test z.frontier_rows==collect(m+1:2m)
        @test all(==(2),vec(sum(z.matrix[:,1:3m];dims=1)))
        @test all(==(1),vec(sum(z.matrix[:,3m+1:end];dims=1)))
        @test toric3d_vertex_z_strip_window(c,t).matrix==z.matrix
        xp=toric3d_vertex_x_prefix_window(c,t);zp=toric3d_vertex_z_prefix_window(c,t)
        @test xp.qubits==zp.qubits==collect(1:3m*(t+1))
        @test xp.checks==plaquette_checks(c)[1:3m*t]
        @test zp.checks==star_checks(c)[1:(t+1)*m]
        for w in (x,z,xp,zp)
            @test all(q->data_qubit_coordinate(c,q)[2]<=t,w.qubits)
            @test all(iseven,Int.(TopoNoise._check_matrix(x.checks,data_qubit_count(c);sparse=false))*transpose(Int.(TopoNoise._check_matrix(z.checks,data_qubit_count(c);sparse=false))))
        end
    end
    for ctor in (toric3d_vertex_x_slice_window,toric3d_vertex_x_complete_window,toric3d_vertex_x_prefix_window,toric3d_vertex_z_strip_window,toric3d_vertex_z_complete_window,toric3d_vertex_z_prefix_window)
        @test_throws ArgumentError ctor(c,0)
        @test_throws ArgumentError ctor(c,3)
    end
end
@testset "H+5 recovery matrices have the specified ranks" begin
    for (u,v) in ((2,2),(2,3),(3,2),(3,3))
        c=Toric3DVertexProjectorCode(3,u,v);m=u*v
        for t in 1:2
            x=toric3d_vertex_x_complete_window(c,t);z=toric3d_vertex_z_complete_window(c,t)
            @test size(x.matrix)==(t==1 ? 3m : 6m,6m)
            @test h5_test_rank(x.matrix)==(t==1 ? 3m-1 : 4m)
            @test size(z.matrix)==(2m,6m)
            @test h5_test_rank(z.matrix)==2m
        end
    end
end
@testset "H+5 period-two geometry preserves physical edge identities" begin
    for dims in ((2,2,2),(2,2,3),(2,3,2))
        c=Toric3DVertexProjectorCode(dims...);w=toric3d_vertex_x_slice_window(c,1);m=c.Lu*c.Lv
        @test length(unique(w.qubits))==3m
        @test length(w.qubits)>length(unique([Tuple(findall(w.matrix[:,j])) for j in axes(w.matrix,2)]))
        @test all(q->data_qubit_index(c,data_qubit_coordinate(c,q)...)==q,w.qubits)
    end
end
end
