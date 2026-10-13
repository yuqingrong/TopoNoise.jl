using Test,Random,TopoNoise,TOML
@testset "H+5 artifacts preserve four-way paired provenance" begin
    @test isdefined(TopoNoise,:save_toric3d_vertex_projector_feedback)
    if isdefined(TopoNoise,:save_toric3d_vertex_projector_feedback)
        grids=Dict(ch=>Toric3DVertexScanGrid([(1,2,2)],[0.,.03]) for ch in (:x_only,:z_only))
        scan=compare_toric3d_vertex_projector_feedback(MersenneTwister(1);channel_grids=grids,shots=2,seed=1234)
        mktempdir() do dir
            paths=save_toric3d_vertex_projector_feedback(scan,dir;figures=false,extra_metadata=Dict("stage"=>"scout"))
            rows=readlines(paths.raw);header=split(first(rows),',')
            @test length(rows)==1+4*23
            @test all("joint_mask_$i" in header for i in 0:15)
            @test all("$(p)_failures" in header for p in TORIC3D_VERTEX_PROTOCOLS)
            @test count(h->endswith(h,"_difference"),header)==6
            @test all(k->k in header,("T","Lu","Lv","p_eff","record_digest","signature_digest","invalid_recoveries","complete"))
            meta=TOML.parsefile(paths.metadata)
            @test meta["circuit"]=="direct-projector H+5"
            @test meta["fault_locations_per_cell"]==11
            @test meta["complete"] && meta["invalid_recoveries"]==0
            @test meta["primary_observables"]==Dict("x_only"=>"logical_x_union","z_only"=>"logical_z_union")
            @test meta["raw_csv_sha256"]==TopoNoise._vertex_sha256(paths.raw)
            @test meta["source_hashes"]==scan.config["source_hashes"]
            @test isfile(paths.aggregate)
            @testset "H+5 figures distinguish logical diagnostics from state failure" begin
                figs=plot_toric3d_vertex_projector_feedback(scan)
                @test Set(keys(figs))==Set((:x_only,:z_only))
                saved=save_toric3d_vertex_projector_feedback(scan,dir)
                @test length(saved.figures)==6
                @test all(p->isfile(p)&&filesize(p)>1000,saved.figures)
                @test occursin("state-independent",meta["logical_z_label"])
            end
        end
    end
end
