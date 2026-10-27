"""
Validate cutoff=1e-26, chi=240 against all 73 saved production outcomes,
including every formerly invalid outcome. This is a diagnostic population,
not a replacement for the independent 1,000-shot production target.

Usage:
  julia --project=. --threads=1 precision_validate.jl RUN_DIR OUTPUT_DIR [--validate-only]

Reuses precision_probe.jl unchanged. The baseline, raw corpus, and independent
exact references are read only. No resampling, clipping, relaxed negativity
tolerance, or exact substitution. Each diagnostic case is checkpointed.
Stop immediately after checkpointing any invalid/exception outcome or a
decision that disagrees with independent exact ML beyond the tie tolerance.
"""
module TNMLPrecisionValidation

using TOML, LinearAlgebra, TopoNoise, ITensors, ITensorMPS
include(joinpath(@__DIR__, "precision_probe.jl"))
const Probe = TNMLPrecisionProbe
const Runner = Probe.Runner
const Model = Probe.RepairedSeedTN
const EXPECTED_RECORDS = 73
const CHI = 240
const CUTOFF = 1e-26

function read_all_inputs(run_dir)
    # Reuse the frozen probe's metadata/source/corpus checks; its five cases
    # are also required members of the complete, unfiltered 73-case set.
    checked = Probe.read_inputs(run_dir)
    metadata = checked.metadata
    corpus = TOML.parsefile(checked.corpus_path)
    metadata["L"] == 3 && metadata["p"] == 0.0025 || error("unexpected diagnostic run")
    selected_shots = collect(1:EXPECTED_RECORDS)
    expected_files = Set(Probe.filename(i) for i in selected_shots)
    for folder in ("shots","exact")
        actual_files = Set(name for name in readdir(joinpath(run_dir,folder))
            if startswith(name,"shot_") && endswith(name,".toml"))
        actual_files == expected_files || error("$folder must contain exactly saved shots 1:$EXPECTED_RECORDS; no omitted or additional cases")
    end
    selected = Dict{Int,Any}()
    for i in selected_shots
        record = corpus["records"][i]
        record["shot_index"] == i || error("corpus shot order differs")
        Runner.identity_hash(Dict(k=>v for (k,v) in record if k != "record_id")) == record["record_id"] ||
            error("saved history $i failed its identity check")
        exact_path = joinpath(run_dir,"exact",Probe.filename(i))
        baseline_path = joinpath(run_dir,"shots",Probe.filename(i))
        exact, baseline = TOML.parsefile(exact_path), TOML.parsefile(baseline_path)
        for (name,doc) in (("exact",exact),("baseline",baseline))
            for (field,value) in (("run_id",metadata["run_id"]),("record_id",record["record_id"]),
                                  ("shot_index",i),("L",metadata["L"]),("p",metadata["p"]),
                                  ("chi",metadata["chi"]),("cutoff",metadata["cutoff"]),
                                  ("true_mask",record["true_mask"]))
                get(doc,field,nothing) == value || error("$name $i has mismatched $field")
            end
        end
        baseline["status"] in ("complete","invalid") || error("baseline $i has no completed numerical outcome")
        baseline["valid"] === (baseline["status"] == "complete") || error("baseline $i validity/status differs")
        length(baseline["probabilities"]) == 4 && all(isfinite,baseline["probabilities"]) ||
            error("baseline $i probabilities are missing or nonfinite")
        exact["status"] == "complete" && exact["valid"] === true || error("exact $i is incomplete")
        exact["algorithm"] == "independent_untruncated_dense_detector_frontier" || error("unexpected exact algorithm")
        probabilities = exact["probabilities"]
        length(probabilities) == 4 && all(isfinite,probabilities) && minimum(probabilities) >= 0 &&
            isapprox(sum(probabilities),1;atol=1e-10,rtol=0) || error("exact $i posterior is invalid")
        exact["logical_mask"] == argmax(probabilities)-1 || error("exact $i label differs from posterior")
        observed = (seed_syndrome=BitVector(record["seed_syndrome"]),
                    final_syndrome=BitVector(record["final_syndrome"]),
                    seed_recovery=BitVector(record["seed_recovery"]),
                    recovery_mask=Int(record["recovery_mask"]))
        selected[i] = (;record,observed,exact,baseline,exact_path,baseline_path)
    end
    return (;metadata,selected,selected_shots,metadata_path=checked.metadata_path,corpus_path=checked.corpus_path)
end

function report_document(probe_id, cases, setup_seconds, baseline_invalid)
    complete = length(cases) == EXPECTED_RECORDS
    valid = [case for case in cases if case["status"] == "complete" && case["valid"]]
    comparisons = [case["exact_comparison"] for case in cases if haskey(case,"exact_comparison")]
    all_passed = length(valid) == length(cases) && length(comparisons) == length(cases) &&
        all(c->c["valid_ml_tie_aware_match"],comparisons)
    return Dict("probe_id"=>probe_id,"expected_cases"=>EXPECTED_RECORDS,"completed_cases"=>length(cases),
        "status"=>!all_passed ? "blocked" : complete ? "complete" : "running",
        "stop_reason"=>all_passed ? "" : "Invalid/exception outcome or independent exact-ML disagreement; no further cases launched.",
        "valid_cases"=>length(valid),"invalid_cases"=>count(case->case["status"]=="invalid",cases),
        "exception_cases"=>count(case->case["status"]=="exception",cases),
        "exact_tie_aware_matches"=>count(c->c["valid_ml_tie_aware_match"],comparisons),
        "exact_strict_matches"=>count(c->c["valid_ml_strict_match"],comparisons),
        "maximum_posterior_linf"=>maximum((c["posterior_linf"] for c in comparisons);init=0.),
        "cap_reached_cases"=>count(case->get(case,"max_bond_cap_reached",false),cases),
        "baseline_invalid_shots"=>baseline_invalid,"setup_seconds"=>setup_seconds,
        "total_case_seconds"=>sum((get(case,"elapsed_seconds",0.) for case in cases);init=0.),
        "population_complete"=>complete,
        "validation_passed"=>complete && all_passed,
        "interpretation"=>"All 73 previously completed histories, including invalid baseline outcomes; diagnostic validation only, no replacement production outcomes or convergence guarantee.",
        "cases"=>cases)
end

function main(args=ARGS)
    2 <= length(args) <= 3 || error("expected RUN_DIR OUTPUT_DIR [--validate-only]")
    length(args) == 2 || args[3] == "--validate-only" || error("unknown option")
    run_dir, output_dir = realpath(args[1]), abspath(args[2])
    (output_dir == run_dir || startswith(output_dir,run_dir*"/")) && error("diagnostics must be outside the production run")
    Threads.nthreads() == 1 || error("use --threads=1")
    BLAS.set_num_threads(1)
    inputs = read_all_inputs(run_dir)
    baseline_invalid = [i for i in inputs.selected_shots if !inputs.selected[i].baseline["valid"]]
    manifest = Dict("schema_version"=>1,"purpose"=>"Validate precision on all 73 saved outcomes before changing the production runner",
        "production_directory"=>run_dir,"production_run_id"=>inputs.metadata["run_id"],
        "production_metadata_sha256"=>Probe.filehash(inputs.metadata_path),
        "corpus_sha256"=>Probe.filehash(inputs.corpus_path),
        "validation_script_sha256"=>Probe.filehash(@__FILE__),
        "probe_script_sha256"=>Probe.filehash(joinpath(@__DIR__,"precision_probe.jl")),
        "selected_shots"=>inputs.selected_shots,"baseline_invalid_shots"=>baseline_invalid,
        "configurations"=>[Dict("chi"=>CHI,"cutoff"=>CUTOFF)],
        "negative_tolerance"=>Probe.NEGATIVE_TOLERANCE,"ml_tie_tolerance"=>Probe.TIE_TOLERANCE,
        "julia_version"=>string(VERSION),"ITensors_version"=>string(Base.pkgversion(ITensors)),
        "ITensorMPS_version"=>string(Base.pkgversion(ITensorMPS)),"blas_threads"=>BLAS.get_num_threads(),
        "truth_policy"=>"Only seed syndrome, final syndrome, known seed recovery and its logical shift enter decode_record.",
        "population_policy"=>"Every saved production outcome in shots1:73 is included, regardless of baseline validity or agreement. No new histories are generated.",
        "references"=>[Dict("shot_index"=>i,"record_id"=>inputs.selected[i].record["record_id"],
            "exact_sha256"=>Probe.filehash(inputs.selected[i].exact_path),
            "baseline_sha256"=>Probe.filehash(inputs.selected[i].baseline_path)) for i in inputs.selected_shots])
    probe_id = Runner.identity_hash(manifest)
    manifest["probe_id"] = probe_id
    manifest_path = joinpath(output_dir,"probe-metadata.toml")
    if isfile(manifest_path)
        TOML.parsefile(manifest_path) == manifest || error("diagnostic metadata differs; use a fresh output directory")
    else
        Runner.atomic_toml(manifest_path,manifest)
    end
    println("Validated all 73 saved histories/exact references, including baseline invalid shots=$baseline_invalid; probe_id=$probe_id")
    flush(stdout)
    "--validate-only" in args && return manifest
    now = time()
    Runner.atomic_toml(joinpath(output_dir,"progress.toml"),Dict("pid"=>getpid(),"stage"=>"setup",
        "shot"=>0,"stage_started_unix"=>now,"updated_at_unix"=>now))
    started = time_ns()
    context = Model.build_context(inputs.metadata["L"])
    settings = Toric3DZTNSettings(max_bond_dimension=CHI,cutoff=CUTOFF,negative_tolerance=Probe.NEGATIVE_TOLERANCE)
    decoder = Model.build_decoder(context,inputs.metadata["p"];settings)
    setup_seconds = (time_ns()-started)/1e9
    cases = Dict{String,Any}[]
    for i in inputs.selected_shots
        push!(cases,Probe.probe_case(decoder,inputs.selected[i],CHI,CUTOFF,probe_id,output_dir))
        report = report_document(probe_id,cases,setup_seconds,baseline_invalid)
        Runner.atomic_toml(joinpath(output_dir,"report.toml"),report)
        if report["status"] == "blocked"
            now = time()
            Runner.atomic_toml(joinpath(output_dir,"progress.toml"),Dict("pid"=>getpid(),"stage"=>"blocked",
                "shot"=>i,"stage_started_unix"=>now,"updated_at_unix"=>now,"reason"=>report["stop_reason"]))
            error("Precision validation stopped after shot $i: $(report["stop_reason"])")
        end
        GC.gc()
    end
    now = time()
    Runner.atomic_toml(joinpath(output_dir,"progress.toml"),Dict("pid"=>getpid(),"stage"=>"complete",
        "shot"=>length(cases),"stage_started_unix"=>now,"updated_at_unix"=>now))
    return cases
end

end # module

abspath(PROGRAM_FILE) == (@__FILE__) && TNMLPrecisionValidation.main()
