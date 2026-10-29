using Test, TOML, SHA

const chi_probe_path = joinpath(@__DIR__, "chi_probe.jl")
@testset "chi convergence diagnostic is available" begin
    @test isfile(chi_probe_path)
end
include(chi_probe_path)
const ChiProbe = TNMLChiProbe

filehash(path) = bytes2hex(sha256(read(path)))

function chi_fixture!(directory; source_digest=nothing)
    source = joinpath(@__DIR__, "model.jl")
    digest = something(source_digest, filehash(source))
    metadata = Dict{String,Any}(
        "schema_version"=>1, "L"=>5, "p"=>0.006, "shots"=>3,
        "chi"=>240, "cutoff"=>1e-26, "negative_tolerance"=>1e-10,
        "exact"=>false, "rng_seed"=>20261029,
        "tn_uses_both_records"=>true,
        "source_sha256"=>Dict("examples/tnml_repaired_seed/model.jl"=>digest),
        "project_sha256"=>filehash(joinpath(dirname(dirname(@__DIR__)), "Project.toml")))
    metadata["run_id"] = ChiProbe.Runner.identity_hash(metadata)
    records = Dict{String,Any}[]
    for shot in 1:3
        record = Dict{String,Any}(
            "shot_index"=>shot, "raw_fault_count"=>4,
            "raw_fault_indices"=>shot == 1 ? [1] : Int[],
            "seed_syndrome"=>[false,true], "seed_recovery"=>[true,false],
            "recovery_mask"=>1, "final_syndrome"=>[false,false,true],
            "true_mask"=>2, "mwpm_failure"=>false,
            "raw_seed_mask"=>1, "transfer_mask"=>3)
        record["record_id"] = ChiProbe.Runner.identity_hash(record)
        push!(records,record)
    end
    ChiProbe.Runner.atomic_toml(joinpath(directory,"metadata.toml"),metadata)
    ChiProbe.Runner.atomic_toml(joinpath(directory,"input-corpus.toml"),Dict(
        "run_id"=>metadata["run_id"], "rng_seed"=>metadata["rng_seed"],
        "requested_shots"=>3, "records"=>records))
    return metadata,records
end

function synthetic_chi_case(chi; valid=true,
        probabilities=[0.55,0.30,0.10,0.05], logz=-20.0,
        status=valid ? "complete" : "computed_invalid", actual_threads=4)
    Dict{String,Any}(
        "status"=>status, "valid"=>valid, "record_id"=>"record-one",
        "shot_index"=>1, "requested_chi"=>chi, "source_chi"=>240,
        "requested_blas_threads"=>4, "actual_blas_threads"=>actual_threads,
        "julia_threads"=>1, "probabilities"=>probabilities,
        "minimum_probability"=>minimum(probabilities),
        "raw_argmax"=>argmax(probabilities)-1,
        "logical_mask"=>valid ? argmax(probabilities)-1 : -1,
        "log_joint_syndrome"=>logz, "elapsed_seconds"=>Float64(chi),
        "peak_rss_bytes"=>chi*1024, "actual_maxbond"=>chi)
end

@testset "probe settings are fixed and explicit" begin
    @test ChiProbe.CHI_CASES == (240,480)
    @test ChiProbe.BLAS_THREADS == 4
    @test ChiProbe.CASE_LIMIT_SECONDS == 7200
    @test ChiProbe.POSTERIOR_LINF_LIMIT == 1e-6
    @test ChiProbe.LOGZ_LIMIT == 1e-6
end

@testset "numerical invalidity is a completed diagnostic status" begin
    @test ChiProbe.numerical_status(true) == "complete"
    @test ChiProbe.numerical_status(false) == "computed_invalid"
end

@testset "saved L5 shot is authenticated and exposed as observed-only" begin
    mktempdir() do run_dir
        metadata,records = chi_fixture!(run_dir)
        inputs = ChiProbe.read_inputs(run_dir)
        @test inputs.metadata == metadata
        @test inputs.record == records[1]
        @test Set(keys(inputs.observed)) ==
            Set((:seed_syndrome,:seed_recovery,:recovery_mask,:final_syndrome))
        @test !(:true_mask in keys(inputs.observed))
        @test !(:raw_fault_indices in keys(inputs.observed))

        changed = deepcopy(metadata)
        changed["source_sha256"]["examples/tnml_repaired_seed/model.jl"] = "0"^64
        changed["run_id"] = ChiProbe.Runner.identity_hash(
            Dict(k=>v for (k,v) in changed if k != "run_id"))
        ChiProbe.Runner.atomic_toml(joinpath(run_dir,"metadata.toml"),changed)
        @test_throws ErrorException ChiProbe.read_inputs(run_dir)
    end
end

@testset "output is external, immutable, and bound to both chi cases" begin
    mktempdir() do parent
        run_dir = joinpath(parent,"source")
        output_dir = joinpath(parent,"diagnostic")
        mkpath(run_dir); chi_fixture!(run_dir)
        inputs = ChiProbe.read_inputs(run_dir)
        manifest = ChiProbe.initialize_probe(inputs,output_dir)
        @test TOML.parsefile(joinpath(output_dir,"probe-metadata.toml")) == manifest
        @test manifest["chi_cases"] == [240,480]
        @test manifest["blas_threads"] == 4
        @test manifest["case_limit_seconds"] == 7200
        @test manifest["record_id"] == inputs.record["record_id"]
        @test ChiProbe.initialize_probe(inputs,output_dir) == manifest
        @test_throws ErrorException ChiProbe.initialize_probe(inputs,run_dir)
        @test_throws ErrorException ChiProbe.initialize_probe(inputs,parent)

        incompatible = deepcopy(manifest)
        incompatible["blas_threads"] = 8
        ChiProbe.Runner.atomic_toml(joinpath(output_dir,"probe-metadata.toml"),incompatible)
        @test_throws ErrorException ChiProbe.initialize_probe(inputs,output_dir)
    end
end

@testset "adjacent comparison retains raw invalid posteriors" begin
    lower = synthetic_chi_case(240; valid=false,
        probabilities=[0.0005,1.0007,-0.0014,0.0002],logz=-68.194)
    higher = synthetic_chi_case(480; valid=true,
        probabilities=[0.0004,0.9991,0.0001,0.0004],logz=-68.1935)
    comparison = ChiProbe.compare_adjacent(lower,higher)
    @test comparison["computed"]
    @test !comparison["both_valid"]
    @test isapprox(comparison["posterior_linf"],0.0016;atol=1e-15)
    @test isapprox(comparison["posterior_l1"],0.0034;atol=1e-15)
    @test isapprox(comparison["abs_log_joint_syndrome_error"],0.0005;atol=1e-12)
    @test comparison["raw_argmax_match"]
    @test !comparison["stable"]
end

@testset "stability requires two valid close results" begin
    lower = synthetic_chi_case(240)
    higher = synthetic_chi_case(480;
        probabilities=[0.55-2e-8,0.30+2e-8,0.10,0.05],logz=-20+4e-8)
    comparison = ChiProbe.compare_adjacent(lower,higher)
    @test comparison["stable"]
    @test comparison["both_valid"]
    @test comparison["tie_aware_label_match"]

    far = synthetic_chi_case(480;
        probabilities=[0.54,0.31,0.10,0.05],logz=-20.01)
    @test !ChiProbe.compare_adjacent(lower,far)["stable"]

    wrong_threads = synthetic_chi_case(480; actual_threads=8)
    @test !ChiProbe.compare_adjacent(lower,wrong_threads)["stable"]
end

@testset "collector preserves computed-invalid results without calling them failures" begin
    mktempdir() do parent
        run_dir = joinpath(parent,"source"); output_dir = joinpath(parent,"diagnostic")
        mkpath(run_dir); chi_fixture!(run_dir)
        inputs = ChiProbe.read_inputs(run_dir)
        manifest = ChiProbe.initialize_probe(inputs,output_dir)
        cases = [synthetic_chi_case(240;valid=false,
                    probabilities=[0.0005,1.0007,-0.0014,0.0002]),
                 synthetic_chi_case(480)]
        for case in cases
            case["probe_id"] = manifest["probe_id"]
            case["source_run_id"] = inputs.metadata["run_id"]
            case["record_id"] = inputs.record["record_id"]
            ChiProbe.Runner.atomic_toml(
                ChiProbe.case_path(output_dir,case["requested_chi"]),case)
        end
        report = ChiProbe.collect_report(inputs,output_dir)
        @test report["status"] == "needs_higher_chi"
        @test report["all_cases_computed"]
        @test report["computed_invalid_cases"] == 1
        @test length(report["comparisons"]) == 1
        @test !report["comparisons"][1]["stable"]
        @test [case["status"] for case in report["cases"]] ==
            ["computed_invalid","complete"]
        @test_throws ErrorException ChiProbe.collect_report(inputs,output_dir)
    end
end

@testset "scheduler failures remain distinct from computed-invalid results" begin
    mktempdir() do parent
        run_dir = joinpath(parent,"source"); output_dir = joinpath(parent,"diagnostic")
        mkpath(run_dir); chi_fixture!(run_dir)
        inputs = ChiProbe.read_inputs(run_dir)
        manifest = ChiProbe.initialize_probe(inputs,output_dir)

        invalid = synthetic_chi_case(240;valid=false,
            probabilities=[0.0005,1.0007,-0.0014,0.0002])
        invalid["probe_id"] = manifest["probe_id"]
        invalid["source_run_id"] = inputs.metadata["run_id"]
        invalid["record_id"] = inputs.record["record_id"]
        ChiProbe.Runner.atomic_toml(ChiProbe.case_path(output_dir,240),invalid)
        @test_throws ErrorException ChiProbe.checkpoint_failure(
            inputs,output_dir,240,"worker_failed",2,7200)

        timeout = ChiProbe.checkpoint_failure(
            inputs,output_dir,480,"contraction_time_limit",124,7200)
        @test timeout["status"] == "contraction_time_limit"
        @test timeout["requested_chi"] == 480
        @test timeout["exit_code"] == 124
        @test_throws ErrorException ChiProbe.checkpoint_failure(
            inputs,output_dir,480,"worker_failed",124,7200)

        report = ChiProbe.collect_report(inputs,output_dir)
        @test report["status"] == "incomplete"
        @test !report["all_cases_computed"]
        @test report["cases"][2]["status"] == "contraction_time_limit"
    end
end

@testset "Slurm wrapper fixes four BLAS threads and bounded resources" begin
    script = read(joinpath(@__DIR__,"chi_probe.sbatch"),String)
    @test occursin("#SBATCH --cpus-per-task=8",script)
    @test occursin("#SBATCH --mem=32G",script)
    @test occursin("#SBATCH --time=05:00:00",script)
    @test occursin("for chi in 240 480",script)
    @test occursin("timeout --signal=TERM --kill-after=30s 2h",script)
    @test occursin("export OPENBLAS_NUM_THREADS=4",script)
    @test occursin("--threads=1,0",script)
    @test occursin("--worker",script)
    @test occursin("--failure",script)
    @test occursin("--collect",script)
    @test occursin("export OMP_NUM_THREADS=1",script)
    @test occursin("export MKL_NUM_THREADS=1",script)
    @test occursin("export BLIS_NUM_THREADS=1",script)
    @test !occursin("exact",lowercase(script))
end
