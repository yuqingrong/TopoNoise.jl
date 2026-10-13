using Test,Random,TopoNoise
isdefined(@__MODULE__,:h5_test_initial) || include("vertex_projector_test_helpers.jl")

# Explicit CSS measurement reference, retaining signs and deterministic outcomes.
function h5_test_project_x(state,signs,s,outcome)
    n=size(state,2)÷2; check=falses(2n);check[s].=true
    anti=findall(r->isodd(count(state[r,n .+ s])),axes(state,1))
    a=copy(state); b=copy(signs)
    if isempty(anti)
        basis=TopoNoise._toric3d_encoded_row_basis(hcat(a,b))
        @test TopoNoise._toric3d_encoded_in_row_space(vcat(check,outcome),basis)
        return a,b,0
    end
    pivot=first(anti)
    for r in anti[2:end]
        a[r,:].⊻=a[pivot,:];b[r]⊻=b[pivot]
    end
    a[pivot,:].=check;b[pivot]=outcome
    a,b,1
end

@testset "H+5 lower projection and upper outcomes have the specified ranks" begin
    @test isdefined(TopoNoise,:toric3d_vertex_boundary_map)
    if isdefined(TopoNoise,:toric3d_vertex_boundary_map)
        for dims in ((1,2,2),(2,2,3),(3,3,2))
            c=Toric3DVertexProjectorCode(dims...);m=c.Lu*c.Lv;map=toric3d_vertex_boundary_map(c)
            @test verify_toric3d_vertex_lower_boundary(map)
            @test map.root==c.T*m+1
            @test h5_test_rank(star_check_matrix(c)[1:m,:])==m
            e=toric3d_vertex_projector_encoder(c)
            state=h5_test_forward(h5_test_initial(c),e.operations)
            upper=star_check_matrix(c)[end-m+1:end,:]
            @test h5_test_rank(isodd.(Int.(upper)*transpose(Int.(state[:,data_qubit_count(c)+1:end]))))==m-1
        end
    end
end
if isdefined(TopoNoise,:toric3d_vertex_boundary_map)
@testset "H+5 upper completion realizes every clean even branch" begin
    c=Toric3DVertexProjectorCode(2,2,2);n=data_qubit_count(c);map=toric3d_vertex_boundary_map(c)
    e=toric3d_vertex_projector_encoder(c);initial=h5_test_forward(h5_test_initial(c),e.operations)
    for mask in 0:7
        b=BitVector([isodd(mask>>k) for k in 0:2]);push!(b,isodd(count(b)))
        state=copy(initial); signs=falses(n);random_outcomes=0
        for (s,o) in zip(map.upper_checks,b)
            state,signs,r=h5_test_project_x(state,signs,s,o);random_outcomes+=r
        end
        @test random_outcomes==3 # each even branch has probability 1/8
        corr=toric3d_vertex_boundary_correction(map,b)
        for row in 1:n
            signs[row]⊻=isodd(count(state[row,1:n].&corr))
        end
        basis=TopoNoise._toric3d_encoded_row_basis(hcat(state,signs))
        for (pauli,checks) in ((:X,star_checks(c)),(:Z,plaquette_checks(c)),(:Z,toric3d_encoded_logical_z(c)))
            for support in checks
                v=falses(2n+1);v[(pauli==:X ? 0 : n).+support].=true
                @test TopoNoise._toric3d_encoded_in_row_space(v,basis)
            end
        end
    end
end
@testset "H+5 paired boundary map equals the explicit two-branch instrument" begin
    c=Toric3DVertexProjectorCode(2,2,2);map=toric3d_vertex_boundary_map(c);n=data_qubit_count(c)
    for clean in 0:7,flip in 0:15
        b=BitVector([isodd(clean>>k) for k in 0:2]);push!(b,isodd(count(b)))
        u=BitVector([isodd(flip>>k) for k in 0:3]);z=falses(n)
        for i in findall(u);z[map.upper_checks[i][3]]=true;end
        a=b.⊻u;out=complete_toric3d_vertex_boundary(map,z,b)
        @test out.observed==a
        @test out.clean_correction==toric3d_vertex_boundary_correction(map,b)
        @test out.physical_correction==toric3d_vertex_boundary_correction(map,a)
        @test out.z==z.⊻out.physical_correction.⊻out.clean_correction
        @test out.z==z.⊻toric3d_vertex_boundary_correction(map,u)
        @test !any(TopoNoise._vertex_syndrome(map.upper_checks,out.z))
    end
    @test_throws DimensionMismatch toric3d_vertex_boundary_correction(map,falses(3))
    @test_throws ArgumentError complete_toric3d_vertex_boundary(map,falses(n),BitVector([1,0,0,0]))
end
@testset "H+5 tree corrections have zero homology and preserve all logical classes" begin
    for dims in ((1,2,2),(2,2,3),(2,3,2),(3,3,3))
        c=Toric3DVertexProjectorCode(dims...);map=toric3d_vertex_boundary_map(c);n=data_qubit_count(c);m=c.Lu*c.Lv
        lx=toric3d_encoded_logical_x(c);lz=toric3d_encoded_logical_z(c);hx=star_check_matrix(c)
        @test count(!iszero,map.parent_edge)==size(hx,1)-1
        for i in 1:m
            a=falses(m);a[i]=true;corr=toric3d_vertex_boundary_correction(map,a)
            s=falses(size(hx,1));s[map.root]=true;s[(c.T+1)*m+i]=true
            @test TopoNoise._vertex_syndrome(hx,corr)==s
            @test all(l->iseven(count(corr[l])),lx)
        end
        for support in lz
            z=falses(n);z[support].=true
            @test complete_toric3d_vertex_boundary(map,z,falses(m)).z==z
        end
        for support in lx
            frame=PauliFrame(falses(n),falses(n));frame.x[support].=true
            out=complete_toric3d_vertex_boundary(map,frame.z,falses(m))
            @test !any(out.z)
            @test frame.x[support]==trues(length(support))
        end
    end
    c=Toric3DVertexProjectorCode(2,2,2);map=toric3d_vertex_boundary_map(c)
    @test map.paths[4]==[32,15,21] # avoids the ordinary BFS parity defect [36,20,13]
end
@testset "H+5 clean upper branch sampling is seeded and even" begin
    map=toric3d_vertex_boundary_map(Toric3DVertexProjectorCode(2,2,2))
    r=MersenneTwister(91);s=MersenneTwister(91)
    branches=[sample_toric3d_vertex_upper_branch(r,map) for _ in 1:1000]
    @test branches==[sample_toric3d_vertex_upper_branch(s,map) for _ in 1:1000]
    @test all(b->length(b)==4 && iseven(count(b)),branches)
    @test length(unique(branches))==8
end
end
