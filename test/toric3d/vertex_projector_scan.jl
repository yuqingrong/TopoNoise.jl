using Test,Random,TopoNoise,TOML
@testset "Sixteen joint patterns determine all six H+5 protocol contrasts" begin
    @test isdefined(TopoNoise,:Toric3DVertexScanGrid)
    if isdefined(TopoNoise,:Toric3DVertexScanGrid)
        m=TopoNoise._vertex_quad_metric(136,collect(1:16))
        @test [a.failures for a in m.protocols]==[72,76,84,100]
        @test length(m.pairs)==6
        @test [(p.neither,p.first_only,p.second_only,p.both) for p in m.pairs]==[(28,32,36,40),(24,28,40,44),(16,20,48,52),(22,30,38,46),(14,22,46,54),(10,26,42,58)]
        @test (m.pairs[1].neither,m.pairs[1].first_only,m.pairs[1].second_only,m.pairs[1].both)==(28,32,36,40)
        @test m.pairs[1].difference≈4/136
        @test all(a->a.ci[1]<a.rate<a.ci[2],m.protocols)
        @test_throws ArgumentError TopoNoise._vertex_quad_metric(1,zeros(Int,16))
    end
end
if isdefined(TopoNoise,:Toric3DVertexScanGrid)
h5_grid=Dict(:x_only=>Toric3DVertexScanGrid([(1,2,2)],[0.,.02]),:z_only=>Toric3DVertexScanGrid([(2,2,3)],[.01]))
h5_run(;kw...)=compare_toric3d_vertex_projector_feedback(MersenneTwister(999);channel_grids=h5_grid,shots=4,seed=1234,batch_size=2,kw...)
h5_counts(s)=[(p.point_seed,[m.joint_counts for m in values(p.metrics)],p.record_digest,p.signature_digest) for p in s.points]
@testset "H+5 scans are four-way paired, channel complete, and seeded" begin
    a=h5_run();b=h5_run()
    @test h5_counts(a)==h5_counts(b)
    @test [(p.channel,p.dimensions,p.p) for p in a.points]==[(:x_only,(1,2,2),0.),(:x_only,(1,2,2),.02),(:z_only,(2,2,3),.01)]
    @test all(m->sum(m.joint_counts)==4,values(a.points[1].metrics))
    @test a.points[1].metrics[:logical_x_union].joint_counts[1]==4
    @test all(p->sum(p.metrics[Symbol("pauli_class_",i)].protocols[1].failures for i in 0:15)==4,a.points)
    @test all(p->p.invalid_recoveries==0,a.points)
    @test a.config["protocols"]==string.(collect(TORIC3D_VERTEX_PROTOCOLS))
    for protocols in ((:none,),reverse(TORIC3D_VERTEX_PROTOCOLS))
        @test_throws ArgumentError h5_run(;protocols,_draw_hook=(_...)->error("draw occurred"))
    end
end
@testset "Interrupted H+5 scans resume exact point substreams" begin
    mktempdir() do dir
        records=Any[];hook=(idx,shot,r)->push!(records,(idx,shot,copy(r.x),copy(r.z),copy(r.upper_branch)))
        full=h5_run(;_draw_hook=hook)
        partial=Any[]
        @test_throws InterruptException h5_run(;checkpoint_dir=dir,_draw_hook=(i,s,r)->begin push!(partial,(i,s,copy(r.x),copy(r.z),copy(r.upper_branch)));s==3 && throw(InterruptException());end)
        saved=TOML.parsefile(joinpath(dir,"checkpoint.toml"))
        @test saved["states"][1]["completed_shots"]==2
        @test saved["states"][1]["invalid_recoveries"]==0
        resumed_records=Any[]
        resumed=compare_toric3d_vertex_projector_feedback(MersenneTwister(1);channel_grids=h5_grid,shots=4,seed=1234,batch_size=2,checkpoint_dir=dir,resume=true,_draw_hook=(i,s,r)->push!(resumed_records,(i,s,copy(r.x),copy(r.z),copy(r.upper_branch))))
        @test h5_counts(full)==h5_counts(resumed)
        @test resumed_records==records[3:end]
    end
end
@testset "H+5 resume rejects every configuration mismatch before sampling" begin
    mktempdir() do dir
        h5_run(;checkpoint_dir=dir)
        original=read(joinpath(dir,"checkpoint.toml"),String)
        for key in ("channel_grids","protocols","protocol_version","shots","batch_size","seed","decoder_settings","backend_versions","source_hashes")
            d=TOML.parse(original);d["config"][key]="changed"
            open(joinpath(dir,"checkpoint.toml"),"w") do io;TOML.print(io,d);end
            @test_throws ArgumentError h5_run(;checkpoint_dir=dir,resume=true,_draw_hook=(_...)->error("draw occurred"))
        end
    end
end
@testset "H+5 invalid recovery context is atomically saved before abort" begin
    mktempdir() do dir
        hook=(i,s,e,r,f)->begin
            w=f.final_x.window;bad=falses(length(w.qubits));bad[1]=true
            TopoNoise._vertex_validate_recovery(w,falses(size(w.matrix,1)),bad)
        end
        @test_throws TopoNoise.Toric3DVertexInvalidRecovery h5_run(;checkpoint_dir=dir,_recovery_hook=hook)
        failure=TOML.parsefile(joinpath(dir,"failure.toml"))
        @test failure["shot"]==1
        @test failure["point_seed"]!=""
        @test all(k->haskey(failure["context"],k),("syndrome","recovery","checks","qubits","checkpoint","route","paths","upper_branch","physical_correction","clean_correction"))
        @test TOML.parsefile(joinpath(dir,"checkpoint.toml"))["states"][1]["invalid_recoveries"]==1
        @test_throws ArgumentError h5_run(;checkpoint_dir=dir,resume=true)
    end
end
@testset "H+5 user interruption is not an invalid recovery" begin
    mktempdir() do dir
        @test_throws InterruptException h5_run(;checkpoint_dir=dir,_recovery_hook=(_...)->throw(InterruptException()))
        @test !isfile(joinpath(dir,"failure.toml"))
        @test TOML.parsefile(joinpath(dir,"checkpoint.toml"))["states"][1]["invalid_recoveries"]==0
    end
end
end
@testset "Frozen selection and stage are bound to exact resume" begin
    mktempdir() do dir
        bound=Dict("stage"=>"production","selection_sha256"=>repeat("a",64),"extension_round"=>0,"parent_selection_sha256"=>"","parent_scout_inputs"=>Any[])
        @test hasmethod(compare_toric3d_vertex_projector_feedback,Tuple{MersenneTwister})
        h5_run(;checkpoint_dir=dir,provenance=bound)
        for key in keys(bound)
            other=deepcopy(bound);other[key]=key=="extension_round" ? 1 : "changed"
            @test_throws ArgumentError h5_run(;checkpoint_dir=dir,resume=true,provenance=other,_draw_hook=(_...)->error("draw occurred"))
        end
        @test_throws ArgumentError save_toric3d_vertex_projector_feedback(h5_run(;provenance=bound),joinpath(dir,"wrong");figures=false,extra_metadata=Dict("selection_sha256"=>repeat("b",64)))
    end
end
@testset "Invalid recovery retains actual residual and canonical physical paths" begin
    c=Toric3DVertexProjectorCode(2,2,2);w=toric3d_vertex_x_slice_window(c,1);s=falses(length(w.checks));r=falses(length(w.qubits));observed=copy(s);observed[1]=true
    err=try TopoNoise._vertex_invalid(w,s,r;observed_post_syndrome=observed,paths=[[w.qubits[1]]]);nothing catch e;e end
    @test err isa TopoNoise.Toric3DVertexInvalidRecovery
    if err isa TopoNoise.Toric3DVertexInvalidRecovery
        @test err.context["observed_post_syndrome"]==observed
        @test err.context["predicted_post_syndrome"]==s
        @test err.context["paths"]==[[w.qubits[1]]]
    end
end
