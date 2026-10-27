using Test, TOML

const runner_path = joinpath(@__DIR__, "run.jl")
@testset "checkpointed runner is available" begin
    @test isfile(runner_path)
end
include(runner_path)
const Runner = TNMLRepairedRunner

@testset "configuration rejects unusable runs before sampling" begin
    cfg = Runner.parse_config(["3", "0.02", "7", "240", "20261027", "/tmp/pilot", "exact=true"])
    @test cfg.L == 3 && cfg.shots == 7 && cfg.exact
    @test !Runner.parse_config(["5", "0.02", "1", "240", "7", "/tmp/pilot"]).exact
    for args in (["1", "0.02", "1", "240", "7", "/tmp/pilot"],
                 ["3", "NaN", "1", "240", "7", "/tmp/pilot"],
                 ["3", "0.5", "1", "240", "7", "/tmp/pilot"],
                 ["3", "0.02", "0", "240", "7", "/tmp/pilot"],
                 ["3", "0.02", "1", "0", "7", "/tmp/pilot"],
                 ["5", "0.02", "1", "240", "7", "/tmp/pilot", "exact=true"])
        @test_throws ArgumentError Runner.parse_config(args)
    end
end

@testset "resumption rejects incompatible source or record identity" begin
    mktempdir() do dir
        path = joinpath(dir, "checkpoint.toml")
        original = Dict("run_id"=>"same-run", "record_id"=>"same-record", "shot_index"=>1,
                        "status"=>"complete", "valid"=>true, "TNfailure"=>false,
                        "probabilities"=>[-1e-12, 1+1e-12, 0., 0.])
        Runner.atomic_toml(path, original)
        @test Runner.resume_document(path, "same-run", "same-record", 1) == original
        @test_throws ArgumentError Runner.resume_document(path, "changed-source", "same-record", 1)
        @test_throws ArgumentError Runner.resume_document(path, "same-run", "different-history", 1)
        @test_throws ArgumentError Runner.resume_document(path, "same-run", "same-record", 2)
        @test TOML.parsefile(path) == original
        @test Runner.identity_hash(Dict("a"=>1, "b"=>[false,true])) ==
              Runner.identity_hash(Dict("b"=>[false,true], "a"=>1))
        @test Runner.identity_hash(Dict("a"=>1)) != Runner.identity_hash(Dict("a"=>2))
    end
end

@testset "saved histories round trip and keep truth out of observed data" begin
    record = (raw_faults=BitVector([true,false,true,false]),
              seed_syndrome=BitVector([true,true]), seed_recovery=BitVector([false,true]),
              recovery_mask=1, final_syndrome=BitVector([false,false]), true_mask=3,
              mwpm_failure=true, raw_seed_mask=2, transfer_mask=0)
    doc = Runner.record_document(record, 1)
    @test doc["raw_fault_indices"] == [1,3]
    @test Runner.restore_record(doc) == record
    observed = Runner.observed_record(record)
    @test Set(keys(observed)) == Set((:seed_syndrome,:seed_recovery,:recovery_mask,:final_syndrome))
    @test observed.recovery_mask == 1
    mktempdir() do dir
        nofaults = merge(record, (raw_faults=falses(4),))
        path = joinpath(dir,"nofaults.toml")
        Runner.atomic_toml(path,Runner.record_document(nofaults,1))
        @test Runner.restore_record(TOML.parsefile(path)) == nofaults
    end
end

@testset "exact decision comparisons tolerate tied maxima without changing scores" begin
    exact = Dict("probabilities"=>[0.5,0.5,0.,0.], "logical_mask"=>0)
    result = (valid=true, probabilities=[0.5-1e-12,0.5+1e-12,0.,0.],logical_mask=1)
    cmp = Runner.compare_exact(result, exact)
    @test cmp["ml_tie_aware_match"]
    @test !cmp["ml_strict_match"]
    @test cmp["posterior_l1"] < 3e-12
    invalid = merge(result, (valid=false,logical_mask=-1))
    @test !haskey(Runner.compare_exact(invalid, exact), "ml_tie_aware_match")
end

@testset "unknown outcomes cannot also count as known TN failures" begin
    mktempdir() do dir
        cfg = Runner.parse_config(["2","0.02","1","240","7",dir])
        records = [Dict("mwpm_failure"=>true)]
        failed = [Dict{String,Any}("status"=>"exception","valid"=>false,"TNfailure"=>true)]
        summary = Runner.summarize(cfg,"run",records,failed)
        @test summary["tn_failures"] == 0
        @test summary["unknown_shots"] == 1
        @test summary["tn_failure_identification_interval"] == [0.,1.]
    end
end

@testset "an initial summary records all saved histories as unknown" begin
    mktempdir() do dir
        cfg = Runner.parse_config(["2","0.02","2","240","7",dir])
        records = [Dict("mwpm_failure"=>true),Dict("mwpm_failure"=>false)]
        summary = Runner.summarize(cfg,"run",records,Dict{String,Any}[])
        @test summary["unknown_shots"] == 2
        @test summary["tn_failure_identification_interval"] == [0.,1.]
        @test summary["attempted_tn_shots"] == 0
        @test summary["completed_tn_shots"] == 0
        @test summary["mwpm_failures"] == 1
        attempted = Runner.summarize(cfg,"run",records,Dict{String,Any}[];attempted_tn_shots=1)
        @test attempted["attempted_tn_shots"] == 1
        @test attempted["unknown_shots"] == 2
    end
end


@testset "cutoff is explicit, validated, and keeps the legacy default" begin
    base = ["3","0.0025","2","240","7","/tmp/precision-config-only"]
    @test Runner.parse_config(base).cutoff == 1e-18
    @test Runner.parse_config(vcat(base,["exact=true"])).cutoff == 1e-18
    cfg = Runner.parse_config(vcat(base,["exact=true","cutoff=1e-26"]))
    @test cfg.exact && cfg.cutoff == 1e-26
    @test Runner.parse_config(vcat(base,["cutoff=1e-26","exact=true"])) == cfg
    @test Runner.parse_config(vcat(base,["cutoff=0"])).cutoff == 0
    for options in (["cutoff=NaN"],["cutoff=Inf"],["cutoff=-1e-26"],["cutoff=1"],
                    ["cutoff=1e-26","cutoff=1e-22"],["exact=true","exact=false"],["other=1"])
        @test_throws ArgumentError Runner.parse_config(vcat(base,options))
    end
    doc = Dict("record_id"=>"record","shot_index"=>1,"true_mask"=>0,"mwpm_failure"=>false)
    @test Runner.shot_base(cfg,"run",doc)["cutoff"] == 1e-26
end

@testset "scientific guards stop failures but accept exact ties" begin
    cfg = Runner.parse_config(["3","0.0025","2","240","7","/tmp/precision-config-only","exact=true"])
    exact = Dict{String,Any}("status"=>"complete","valid"=>true,"probabilities"=>[0.5,0.5,0.,0.])
    output = Dict{String,Any}("status"=>"complete","valid"=>true,"logical_mask"=>1,
        "exact_comparison"=>Dict("ml_tie_aware_match"=>true),"TNfailure"=>false)
    @test Runner.exact_halt_reason(cfg,exact) === nothing
    @test Runner.output_halt_reason(cfg,output,exact) === nothing
    @test Runner.exact_halt_reason(cfg,nothing) == "missing_exact_reference"
    @test Runner.exact_halt_reason(cfg,merge(exact,Dict("status"=>"exception"))) == "invalid_exact_reference"
    @test Runner.exact_halt_reason(cfg,merge(exact,Dict("probabilities"=>[NaN,0.,0.,0.]))) == "invalid_exact_posterior"
    @test Runner.output_halt_reason(cfg,merge(output,Dict("valid"=>false)),exact) == "invalid_or_exception_tn_output"
    @test Runner.output_halt_reason(cfg,merge(output,Dict("status"=>"exception")),exact) == "invalid_or_exception_tn_output"
    @test Runner.output_halt_reason(cfg,merge(output,Dict("logical_mask"=>2)),exact) == "exact_ml_disagreement"
    @test Runner.output_halt_reason(cfg,merge(output,Dict("exact_comparison"=>Dict())),exact) == "missing_or_failed_exact_comparison"
    @test output["TNfailure"] === false # guards never modify logical scores
end

@testset "resumed failures stop before complete-run return and preserve the denominator" begin
    for mode in (:invalid, :missing_exact, :disagreement)
        mktempdir() do dir
            cfg = Runner.parse_config(["3","0.0025","2","240","7",dir,"exact=true","cutoff=1e-26"])
            records = [Dict("mwpm_failure"=>false),Dict("mwpm_failure"=>false)]
            output = Dict{String,Any}("status"=>"complete","valid"=>true,"logical_mask"=>0,
                "shot_index"=>1,"record_id"=>"record","TNfailure"=>false,
                "exact_comparison"=>Dict("ml_tie_aware_match"=>true))
            exact = Dict{String,Any}("run_id"=>"run","record_id"=>"record","shot_index"=>1,
                "status"=>"complete","valid"=>true,"probabilities"=>[1.,0.,0.,0.])
            mode == :invalid && (output["valid"] = false; output["status"] = "invalid")
            mode == :disagreement && (output["logical_mask"] = 1; output["TNfailure"] = true)
            mode != :missing_exact && Runner.atomic_toml(joinpath(dir,"exact",Runner.shot_name(1)),exact)
            outputs = [output]
            @test_throws ErrorException Runner.validate_resumed_outputs(cfg,"run",records,outputs)
            summary = TOML.parsefile(joinpath(dir,"summary.toml"))
            halt = TOML.parsefile(joinpath(dir,"halt.toml"))
            progress = TOML.parsefile(joinpath(dir,"progress.toml"))
            @test summary["requested_shots"] == 2
            @test summary["unknown_shots"] == (mode == :invalid ? 2 : 1)
            @test halt["shot_index"] == 1
            @test progress["stage"] == "blocked"
            @test length(outputs) == 1
        end
    end
end
