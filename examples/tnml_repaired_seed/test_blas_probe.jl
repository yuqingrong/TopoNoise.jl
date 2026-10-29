using Test, TOML, SHA

const probe_path = joinpath(@__DIR__, "blas_probe.jl")
@testset "BLAS diagnostic is available" begin
    @test isfile(probe_path)
end
include(probe_path)
const Probe = TNMLBLASProbe

filehash(path) = bytes2hex(sha256(read(path)))

function fixture!(directory; source_digest=nothing)
    source = joinpath(@__DIR__, "model.jl")
    digest = something(source_digest, filehash(source))
    metadata = Dict{String,Any}(
        "schema_version"=>1, "L"=>5, "p"=>0.01, "shots"=>3,
        "chi"=>240, "cutoff"=>1e-26, "negative_tolerance"=>1e-10,
        "exact"=>false, "rng_seed"=>20261029,
        "tn_uses_both_records"=>true,
        "source_sha256"=>Dict("examples/tnml_repaired_seed/model.jl"=>digest),
        "project_sha256"=>filehash(joinpath(dirname(dirname(@__DIR__)), "Project.toml")))
    metadata["run_id"] = Probe.Runner.identity_hash(metadata)
    records = Dict{String,Any}[]
    for shot in 1:3
        record = Dict{String,Any}(
            "shot_index"=>shot, "raw_fault_count"=>4,
            "raw_fault_indices"=>shot == 1 ? [1] : Int[],
            "seed_syndrome"=>[false,true], "seed_recovery"=>[true,false],
            "recovery_mask"=>1, "final_syndrome"=>[false,false,true],
            "true_mask"=>2, "mwpm_failure"=>false,
            "raw_seed_mask"=>1, "transfer_mask"=>3)
        record["record_id"] = Probe.Runner.identity_hash(record)
        push!(records,record)
    end
    Probe.Runner.atomic_toml(joinpath(directory,"metadata.toml"),metadata)
    Probe.Runner.atomic_toml(joinpath(directory,"input-corpus.toml"),Dict(
        "run_id"=>metadata["run_id"], "rng_seed"=>metadata["rng_seed"],
        "requested_shots"=>3, "records"=>records))
    return metadata,records
end

function synthetic_case(threads; probabilities=[0.5,0.5,0.0,0.0],
                        label=threads == 1 ? 0 : 1, logz=-12.0,
                        actual_threads=threads)
    Dict{String,Any}(
        "status"=>"complete", "valid"=>true, "record_id"=>"record-one",
        "shot_index"=>1, "requested_blas_threads"=>threads,
        "actual_blas_threads"=>actual_threads, "probabilities"=>probabilities,
        "logical_mask"=>label, "log_joint_syndrome"=>logz,
        "elapsed_seconds"=>Float64(threads), "peak_rss_bytes"=>threads*1024)
end

@testset "saved L5 shot one is authenticated and observed-only" begin
    mktempdir() do run_dir
        metadata,records = fixture!(run_dir)
        inputs = Probe.read_inputs(run_dir)
        @test inputs.metadata == metadata
        @test inputs.record == records[1]
        @test Set(keys(inputs.observed)) ==
            Set((:seed_syndrome,:seed_recovery,:recovery_mask,:final_syndrome))
        @test inputs.record["shot_index"] == 1
        @test inputs.record["record_id"] == records[1]["record_id"]
        @test inputs.metadata_sha256 == filehash(joinpath(run_dir,"metadata.toml"))
        @test inputs.corpus_sha256 == filehash(joinpath(run_dir,"input-corpus.toml"))

        changed = deepcopy(metadata)
        changed["source_sha256"]["examples/tnml_repaired_seed/model.jl"] = "0"^64
        changed["run_id"] = Probe.Runner.identity_hash(Dict(k=>v for (k,v) in changed if k != "run_id"))
        Probe.Runner.atomic_toml(joinpath(run_dir,"metadata.toml"),changed)
        @test_throws ErrorException Probe.read_inputs(run_dir)
    end
end

@testset "diagnostic output is separate and identity-bound" begin
    mktempdir() do parent
        run_dir = joinpath(parent,"source")
        output_dir = joinpath(parent,"diagnostic")
        mkpath(run_dir); fixture!(run_dir)
        inputs = Probe.read_inputs(run_dir)
        manifest = Probe.initialize_probe(inputs,output_dir)
        @test TOML.parsefile(joinpath(output_dir,"probe-metadata.toml")) == manifest
        @test manifest["blas_threads"] == [1,4,8]
        @test manifest["record_id"] == inputs.record["record_id"]
        @test Probe.initialize_probe(inputs,output_dir) == manifest
        @test_throws ErrorException Probe.initialize_probe(inputs,run_dir)
        @test !isfile(joinpath(run_dir,"probe-metadata.toml"))

        nested_target = joinpath(run_dir,"diagnostic-target")
        symlink_path = joinpath(parent,"linked-output")
        mkpath(nested_target)
        symlink(nested_target,symlink_path)
        @test_throws ErrorException Probe.initialize_probe(inputs,symlink_path)
        @test !isfile(joinpath(nested_target,"probe-metadata.toml"))
    end
end

@testset "cross-thread comparison is strict but tie aware" begin
    baseline = synthetic_case(1)
    candidate = synthetic_case(4; probabilities=[0.5-2e-12,0.5+2e-12,0.0,0.0],label=1,
                               logz=-12.0+5e-11)
    comparison = Probe.compare_case(baseline,candidate)
    @test comparison["passed"]
    @test comparison["posterior_linf"] <= Probe.POSTERIOR_LINF_LIMIT
    @test comparison["abs_log_joint_syndrome_error"] <= Probe.LOGZ_LIMIT
    @test comparison["tie_aware_label_match"]

    @test !Probe.compare_case(baseline,synthetic_case(8;
        probabilities=[0.51,0.49,0.0,0.0],label=0))["passed"]
    @test !Probe.compare_case(baseline,synthetic_case(8;
        probabilities=[0.5,0.5,0.0,0.0],label=1,logz=-11.0))["passed"]
    @test !Probe.compare_case(baseline,synthetic_case(8;
        probabilities=[0.5,0.5,0.0,0.0],label=1,actual_threads=7))["passed"]
    @test Probe.require_passed_report(Dict("status"=>"passed"))["status"] == "passed"
    @test_throws ErrorException Probe.require_passed_report(Dict("status"=>"failed"))
end

@testset "collector requires all three cases and never changes them" begin
    mktempdir() do parent
        run_dir = joinpath(parent,"source"); output_dir = joinpath(parent,"diagnostic")
        mkpath(run_dir); fixture!(run_dir)
        inputs = Probe.read_inputs(run_dir)
        manifest = Probe.initialize_probe(inputs,output_dir)
        originals = Dict{Int,Dict{String,Any}}()
        for threads in Probe.BLAS_THREAD_CASES
            case = synthetic_case(threads)
            case["probe_id"] = manifest["probe_id"]
            case["record_id"] = inputs.record["record_id"]
            originals[threads] = case
            Probe.Runner.atomic_toml(Probe.case_path(output_dir,threads),case)
        end
        report = Probe.collect_report(inputs,output_dir)
        @test report["status"] == "passed"
        @test length(report["cases"]) == 3
        @test length(report["comparisons"]) == 2
        @test all(row["passed"] for row in report["comparisons"])
        for threads in Probe.BLAS_THREAD_CASES
            @test TOML.parsefile(Probe.case_path(output_dir,threads)) == originals[threads]
        end
        @test_throws ErrorException Probe.collect_report(inputs,output_dir)
    end
end

@testset "worker failures retain authenticated provenance" begin
    mktempdir() do parent
        run_dir = joinpath(parent,"source"); output_dir = joinpath(parent,"diagnostic")
        mkpath(run_dir); fixture!(run_dir)
        inputs = Probe.read_inputs(run_dir)
        manifest = Probe.initialize_probe(inputs,output_dir)
        timeout = Probe.checkpoint_failure(inputs,output_dir,4,
            "contraction_time_limit",124,7200)
        @test timeout["probe_id"] == manifest["probe_id"]
        @test timeout["status"] == "contraction_time_limit"
        @test timeout["exit_code"] == 124
        @test timeout["limit_seconds"] == 7200
        @test Probe.checkpoint_failure(inputs,output_dir,4,
            "contraction_time_limit",124,7200) == timeout

        failed = Probe.checkpoint_failure(inputs,output_dir,8,"worker_failed",2,7200)
        @test failed["status"] == "worker_failed"
        @test_throws ErrorException Probe.checkpoint_failure(inputs,output_dir,8,
            "contraction_time_limit",2,7200)

        baseline = synthetic_case(1)
        baseline["probe_id"] = manifest["probe_id"]
        baseline["record_id"] = inputs.record["record_id"]
        Probe.Runner.atomic_toml(Probe.case_path(output_dir,1),baseline)
        report = Probe.collect_report(inputs,output_dir)
        @test [case["status"] for case in report["cases"]] ==
            ["complete","contraction_time_limit","worker_failed"]
        @test report["status"] == "failed"
    end

    mktempdir() do parent
        run_dir = joinpath(parent,"source"); output_dir = joinpath(parent,"diagnostic")
        mkpath(run_dir); fixture!(run_dir)
        inputs = Probe.read_inputs(run_dir)
        manifest = Probe.initialize_probe(inputs,output_dir)
        complete = synthetic_case(4)
        complete["probe_id"] = manifest["probe_id"]
        complete["record_id"] = inputs.record["record_id"]
        Probe.Runner.atomic_toml(Probe.case_path(output_dir,4),complete)
        @test_throws ErrorException Probe.checkpoint_failure(inputs,output_dir,4,
            "worker_failed",2,7200)
        @test !isfile(Probe.failure_path(output_dir,4))
    end
end

@testset "Slurm wrapper reserves eight cores and runs fresh workers sequentially" begin
    script = read(joinpath(@__DIR__,"blas_probe.sbatch"),String)
    @test occursin("#SBATCH --cpus-per-task=8",script)
    @test occursin("#SBATCH --mem=32G",script)
    @test occursin("#SBATCH --time=08:00:00",script)
    @test occursin("for blas_threads in 1 4 8",script)
    @test occursin("--worker",script)
    @test occursin("--failure",script)
    @test occursin("--collect",script)
    @test occursin("export OMP_NUM_THREADS=1",script)
    @test occursin("export MKL_NUM_THREADS=1",script)
    @test occursin("export BLIS_NUM_THREADS=1",script)
    @test !occursin("OMP_NUM_THREADS=\"\$blas_threads\"",script)
    @test !occursin("MKL_NUM_THREADS=\"\$blas_threads\"",script)
    @test !occursin("BLIS_NUM_THREADS=\"\$blas_threads\"",script)
end
