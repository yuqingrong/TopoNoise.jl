using Test
const RUNNER=joinpath(@__DIR__,"run_study.jl")
@testset "Axis comparison runner is implemented" begin
    @test isfile(RUNNER)
end
using TopoNoise, Random
import TOML
include(RUNNER)
const Study=H6AxisStudy

function options(directory,channel;resume=false,max_batches=typemax(Int))
    args=["--channel",channel,"--output",directory,"--sizes","3","--rates","0,0.12",
        "--shots","5","--seed","20261018","--batch-size","2","--max-batches",string(max_batches)]
    resume && push!(args,"--resume")
    Study.parse_arguments(args)
end

@testset "Axis circuit and sampler provenance are explicit" begin
    for channel in ("x_only","z_only")
        protocol=Study.protocol_metadata(channel)
        @test protocol["growth_direction"] == "physical_x"
        @test occursin("square",protocol["seed_geometry"])
        @test protocol["transfer_outputs_per_cell"] == 13
        @test protocol["seed_outputs_per_cell"] == (channel == "z_only" ? 2 : 0)
        context=Study.build_context(3,0.12,channel)
        raw=Study.sample_faults(Random.Xoshiro(77),context,0.12,channel)
        saved=deepcopy(raw)
        actual=Study.three_outcomes(context,raw,channel)
        expected=collect(values(Study.AxisH6.three_outcomes(context,raw)))
        @test raw.seed == saved.seed && raw.x == saved.x && raw.z == saved.z
        for (a,b) in zip(actual,expected), key in (:logical_1,:logical_2,:failure,:residual)
            @test getproperty(a,key) == getproperty(b,key)
        end
        @test_throws ArgumentError Study.sample_faults(Random.Xoshiro(77),context,0.15,channel)
    end
    hashes=Study.source_hashes()
    @test all(haskey(hashes,p) for p in ("AxisH6.jl","axis_geometry.jl","axis_feedback.jl","run_study.jl"))
end

@testset "Paired histogram preserves each marginal" begin
    config=options("/private/tmp/unused-axis-runner-test","x_only").config
    config["shots"]=8
    checkpoint=Study.fresh_checkpoint(config,Dict(),Dict())
    state=first(checkpoint["states"])
    for mask in 0:7
        outcomes=[(logical_1=(mask & (1 << (arm-1))) != 0,logical_2=false,
            failure=(mask & (1 << (arm-1))) != 0) for arm in 1:3]
        Study.accumulate!(state,outcomes)
    end
    @test state["joint_counts"] == ones(Int,8)
    @test state["failures"] == [4,4,4]
    @test Study.validate_states(checkpoint) === checkpoint
    damaged=deepcopy(checkpoint); damaged["states"][1]["failures"][1]+=1
    @test_throws ArgumentError Study.validate_states(damaged)
    damaged=deepcopy(checkpoint); damaged["point_seeds_hex"][1]="0000000000000000"
    @test_throws ArgumentError Study.validate_states(damaged)
end

@testset "Multiple-rate resume preserves exact original streams" begin
    for channel in ("x_only","z_only"), interruption in (2,4)
        mktempdir() do directory
            full=joinpath(directory,"full"); resumed=joinpath(directory,"resumed")
            Study.run_study(options(full,channel);io=devnull)
            partial=Study.run_study(options(resumed,channel;max_batches=interruption);io=devnull)
            @test partial["completed_shots"] == (interruption == 2 ? 4 : 7)
            @test !partial["complete"]
            write(joinpath(resumed,"raw.csv"),"stale")
            rm(joinpath(resumed,"metadata.toml"))
            done=Study.run_study(options(resumed,channel;resume=true);io=devnull)
            @test done["complete"]
            a=TOML.parsefile(joinpath(full,"checkpoint.toml"))
            b=TOML.parsefile(joinpath(resumed,"checkpoint.toml"))
            @test a["point_seeds_hex"] == b["point_seeds_hex"]
            for (sa,sb) in zip(a["states"],b["states"]), key in
                    ("completed_shots","logical_1_counts","logical_2_counts","failures","joint_counts")
                @test sa[key] == sb[key]
            end
            @test length(readlines(joinpath(resumed,"raw.csv"))) == 3
            @test isfile(joinpath(resumed,"previous-checkpoint.toml"))
            # A final committed checkpoint repairs stale mirrors on resume.
            write(joinpath(resumed,"raw.csv"),"stale again")
            @test Study.run_study(options(resumed,channel;resume=true);io=devnull)["complete"]
            @test length(readlines(joinpath(resumed,"raw.csv"))) == 3
            original=read(joinpath(resumed,"checkpoint.toml"),String)
            altered=deepcopy(b); altered["states"][1]["failures"][1]+=1
            write(joinpath(resumed,"checkpoint.toml"),Study.toml_text(altered))
            @test_throws ArgumentError Study.run_study(options(resumed,channel;resume=true);io=devnull)
            write(joinpath(resumed,"checkpoint.toml"),original)
            @test_throws ArgumentError Study.load_checkpoint(resumed,b["config"],b["runtime"],
                merge(b["source_hashes"],Dict("axis_feedback.jl"=>"changed")))
            write(joinpath(resumed,"previous-checkpoint.toml"),"[corrupt")
            @test_throws ArgumentError Study.run_study(options(resumed,channel;resume=true);io=devnull)
        end
    end
end
