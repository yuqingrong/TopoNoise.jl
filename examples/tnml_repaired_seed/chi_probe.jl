"""
Read-only bond-dimension convergence diagnostic for one saved L=5 history.

The authenticated source corpus is never modified or resampled. Fresh Julia
processes replay its first observed history at chi=240 and chi=480 with one
Julia thread and four BLAS threads. A numerically inadmissible posterior is a
completed diagnostic result (`computed_invalid`), not a worker failure.

Usage:
  julia --project=. --threads=1 chi_probe.jl --validate RUN_DIR OUTPUT_DIR
  julia --project=. --threads=1 chi_probe.jl --worker RUN_DIR OUTPUT_DIR CHI
  julia --project=. --threads=1 chi_probe.jl --failure RUN_DIR OUTPUT_DIR CHI STATUS EXIT LIMIT
  julia --project=. --threads=1 chi_probe.jl --collect RUN_DIR OUTPUT_DIR
"""
module TNMLChiProbe

using TOML, LinearAlgebra, TopoNoise, ITensors, ITensorMPS
include(joinpath(@__DIR__,"blas_probe.jl"))
const BaseProbe = TNMLBLASProbe
const Runner = BaseProbe.Runner
const Model = BaseProbe.Model

const CHI_CASES = (240,480)
const BLAS_THREADS = 4
const CASE_LIMIT_SECONDS = 7200
const POSTERIOR_LINF_LIMIT = 1e-6
const LOGZ_LIMIT = 1e-6
const ML_TIE_TOLERANCE = 1e-9
const NEGATIVE_TOLERANCE = BaseProbe.NEGATIVE_TOLERANCE

numerical_status(valid::Bool) = valid ? "complete" : "computed_invalid"

read_inputs(run_directory) = BaseProbe.read_inputs(run_directory)
filehash(path) = BaseProbe.filehash(path)
canonical_destination(path) = BaseProbe.canonical_destination(path)
path_within(path,parent) = BaseProbe.path_within(path,parent)

case_directory(output_dir,chi) =
    joinpath(canonical_destination(output_dir),"chi-"*lpad(string(chi),4,'0'))
case_path(output_dir,chi) = joinpath(case_directory(output_dir,chi),"result.toml")
failure_path(output_dir,chi) = joinpath(case_directory(output_dir,chi),"failure.toml")

function require_chi(chi)
    chi isa Integer && !(chi isa Bool) && chi in CHI_CASES ||
        error("chi must be one of $(collect(CHI_CASES))")
    return Int(chi)
end

function manifest_document(inputs)
    script_path = @__FILE__
    sbatch_path = joinpath(@__DIR__,"chi_probe.sbatch")
    dependency_path = joinpath(@__DIR__,"blas_probe.jl")
    isfile(sbatch_path) || error("chi probe Slurm wrapper is missing")
    isfile(dependency_path) || error("authenticated input helper is missing")
    document = Dict{String,Any}(
        "schema_version"=>1,
        "purpose"=>"Read-only fixed-history bond-dimension convergence diagnostic; excluded from production statistics",
        "source_run_directory"=>inputs.run_dir,
        "source_run_id"=>inputs.metadata["run_id"],
        "source_metadata_sha256"=>inputs.metadata_sha256,
        "source_corpus_sha256"=>inputs.corpus_sha256,
        "record_id"=>inputs.record["record_id"], "shot_index"=>1,
        "L"=>inputs.metadata["L"], "p"=>inputs.metadata["p"],
        "source_chi"=>inputs.metadata["chi"],
        "chi_cases"=>collect(CHI_CASES),
        "cutoff"=>inputs.metadata["cutoff"],
        "negative_tolerance"=>NEGATIVE_TOLERANCE,
        "blas_threads"=>BLAS_THREADS, "julia_threads"=>1,
        "case_limit_seconds"=>CASE_LIMIT_SECONDS,
        "posterior_linf_limit"=>POSTERIOR_LINF_LIMIT,
        "abs_log_joint_syndrome_limit"=>LOGZ_LIMIT,
        "ml_tie_tolerance"=>ML_TIE_TOLERANCE,
        "probe_script_sha256"=>filehash(script_path),
        "sbatch_sha256"=>filehash(sbatch_path),
        "input_helper_sha256"=>filehash(dependency_path),
        "truth_policy"=>"Only seed syndrome, final syndrome, known seed recovery and its logical shift enter decode_record. Raw faults and truth labels are excluded.",
        "interpretation"=>"Adjacent-chi agreement is an approximate convergence diagnostic, not an exact-ML certificate.")
    document["probe_id"] = Runner.identity_hash(document)
    return document
end

"Create or verify immutable diagnostic metadata outside and disjoint from the source run."
function initialize_probe(inputs,output_directory)
    output_dir = canonical_destination(output_directory)
    (path_within(output_dir,inputs.run_dir) || path_within(inputs.run_dir,output_dir)) &&
        error("diagnostic output and saved source run must be disjoint")
    manifest = manifest_document(inputs)
    path = joinpath(output_dir,"probe-metadata.toml")
    if ispath(output_dir) && !isdir(output_dir)
        error("diagnostic output path is not a directory")
    elseif isdir(output_dir) && !isfile(path) && !isempty(readdir(output_dir))
        error("nonempty diagnostic directory has no compatible manifest")
    elseif isfile(path)
        TOML.parsefile(path) == manifest ||
            error("diagnostic manifest differs; use a fresh output directory")
    else
        Runner.atomic_toml(path,manifest)
    end
    return manifest
end

function progress(case_dir,stage,chi)
    now = time()
    document = Dict{String,Any}(
        "pid"=>getpid(), "stage"=>stage, "shot"=>1,
        "chi"=>chi, "blas_threads"=>BLAS_THREADS,
        "stage_started_unix"=>now, "updated_at_unix"=>now)
    Runner.atomic_toml(joinpath(case_dir,"progress.toml"),document)
    println("stage=$stage shot=1/1 chi=$chi BLAS=$BLAS_THREADS")
    flush(stdout)
end

function component_document(result)
    return Dict("valid"=>result.valid,
        "probabilities"=>collect(result.probabilities),
        "logical_mask"=>result.logical_mask,
        "log_syndrome"=>result.log_syndrome,
        "diagnostics"=>Runner.toml_value(result.diagnostics))
end

"Run one fixed observed history. Numerical invalidity is persisted and returns normally."
function run_worker(inputs,output_directory,requested_chi::Integer)
    chi = require_chi(requested_chi)
    Threads.nthreads() == 1 || error("run this diagnostic with one Julia thread")
    manifest = initialize_probe(inputs,output_directory)
    case_dir = case_directory(output_directory,chi)
    result_path = case_path(output_directory,chi)
    ispath(case_dir) && error("chi=$chi output already exists; use a fresh diagnostic directory")
    mkpath(case_dir)

    BLAS.set_num_threads(BLAS_THREADS)
    actual_threads = BLAS.get_num_threads()
    actual_threads == BLAS_THREADS ||
        error("BLAS applied $actual_threads threads, requested $BLAS_THREADS")
    progress(case_dir,"setup",chi)
    setup_started = time_ns()
    context = Model.build_context(inputs.metadata["L"])
    inputs.record["raw_fault_count"] == size(context.model.matrix,2) ||
        error("saved raw fault count differs from the reconstructed model")
    settings = Toric3DZTNSettings(max_bond_dimension=chi,
        cutoff=inputs.metadata["cutoff"],negative_tolerance=NEGATIVE_TOLERANCE)
    decoder = Model.build_decoder(context,inputs.metadata["p"];settings)
    setup_seconds = (time_ns()-setup_started)/1e9
    progress(case_dir,"tn",chi)

    output = Dict{String,Any}(
        "schema_version"=>1, "probe_id"=>manifest["probe_id"],
        "source_run_id"=>inputs.metadata["run_id"],
        "record_id"=>inputs.record["record_id"], "shot_index"=>1,
        "L"=>inputs.metadata["L"], "p"=>inputs.metadata["p"],
        "source_chi"=>inputs.metadata["chi"], "requested_chi"=>chi,
        "cutoff"=>inputs.metadata["cutoff"],
        "requested_blas_threads"=>BLAS_THREADS,
        "actual_blas_threads"=>actual_threads,
        "julia_threads"=>Threads.nthreads(),
        "blas_vendor"=>string(BLAS.vendor()),
        "blas_config"=>sprint(show,BLAS.get_config()),
        "julia_version"=>string(VERSION),
        "ITensors_version"=>string(Base.pkgversion(ITensors)),
        "ITensorMPS_version"=>string(Base.pkgversion(ITensorMPS)),
        "setup_seconds"=>setup_seconds)
    started = time_ns()
    caught = nothing
    try
        result = Model.decode_record(decoder,inputs.observed)
        probabilities = collect(result.probabilities)
        output["valid"] = result.valid
        output["status"] = numerical_status(result.valid)
        output["probabilities"] = probabilities
        output["minimum_probability"] = minimum(probabilities)
        output["raw_argmax"] = argmax(probabilities)-1
        output["logical_mask"] = result.logical_mask
        output["log_joint_syndrome"] = result.log_joint_syndrome
        output["seed_result"] = component_document(result.seed_result)
        output["transfer_result"] = component_document(result.transfer_result)
        output["actual_maxbond"] = max(result.seed_result.diagnostics.actual_maxbond,
            result.transfer_result.diagnostics.actual_maxbond)
        output["max_bond_cap_reached"] = output["actual_maxbond"] >= chi
    catch exception
        exception isa InterruptException && rethrow()
        caught = (exception,catch_backtrace())
        output["status"] = "exception"
        output["valid"] = false
        merge!(output,Runner.exception_details(caught...))
    end
    output["elapsed_seconds"] = (time_ns()-started)/1e9
    output["peak_rss_bytes"] = Int(Sys.maxrss())
    Runner.atomic_toml(result_path,output)
    progress(case_dir,string(output["status"]),chi)
    caught === nothing || error("chi=$chi worker raised an exception; result retained")
    return output
end

function raw_probabilities(case)
    probabilities = get(case,"probabilities",nothing)
    probabilities isa AbstractVector && length(probabilities) == 4 || return nothing
    all(x->x isa Real && isfinite(x),probabilities) || return nothing
    abs(sum(probabilities)-1) <= 1e-7 || return nothing
    return Float64.(probabilities)
end

function is_computed(case)
    get(case,"status",nothing) in ("complete","computed_invalid") || return false
    get(case,"requested_blas_threads",nothing) == BLAS_THREADS || return false
    get(case,"actual_blas_threads",nothing) == BLAS_THREADS || return false
    get(case,"julia_threads",nothing) == 1 || return false
    return raw_probabilities(case) !== nothing
end

function is_valid_case(case)
    is_computed(case) || return false
    get(case,"status",nothing) == "complete" || return false
    get(case,"valid",false) === true || return false
    probabilities = raw_probabilities(case)
    return minimum(probabilities) >= -NEGATIVE_TOLERANCE
end

"Compare raw adjacent-chi outputs, even when one completed output is inadmissible."
function compare_adjacent(lower,higher)
    lower_probabilities = raw_probabilities(lower)
    higher_probabilities = raw_probabilities(higher)
    computed = is_computed(lower) && is_computed(higher) &&
        lower_probabilities !== nothing && higher_probabilities !== nothing
    linf = computed ? maximum(abs.(lower_probabilities.-higher_probabilities)) : Inf
    l1 = computed ? sum(abs.(lower_probabilities.-higher_probabilities)) : Inf
    lower_logz = get(lower,"log_joint_syndrome",NaN)
    higher_logz = get(higher,"log_joint_syndrome",NaN)
    log_error = lower_logz isa Real && higher_logz isa Real &&
        isfinite(lower_logz) && isfinite(higher_logz) ? abs(lower_logz-higher_logz) : Inf
    lower_argmax = computed ? argmax(lower_probabilities)-1 : -1
    higher_argmax = computed ? argmax(higher_probabilities)-1 : -1
    argmax_match = computed && lower_argmax == higher_argmax
    tie_match = computed &&
        maximum(lower_probabilities)-lower_probabilities[higher_argmax+1] <= ML_TIE_TOLERANCE &&
        maximum(higher_probabilities)-higher_probabilities[lower_argmax+1] <= ML_TIE_TOLERANCE
    both_valid = is_valid_case(lower) && is_valid_case(higher)
    stable = computed && both_valid && tie_match &&
        linf <= POSTERIOR_LINF_LIMIT && log_error <= LOGZ_LIMIT
    return Dict{String,Any}(
        "lower_chi"=>get(lower,"requested_chi",-1),
        "higher_chi"=>get(higher,"requested_chi",-1),
        "computed"=>computed, "both_valid"=>both_valid,
        "posterior_linf"=>linf, "posterior_l1"=>l1,
        "abs_log_joint_syndrome_error"=>log_error,
        "lower_raw_argmax"=>lower_argmax,
        "higher_raw_argmax"=>higher_argmax,
        "raw_argmax_match"=>argmax_match,
        "tie_aware_label_match"=>tie_match,
        "stable"=>stable)
end

"Persist timeout/process failure without replacing a completed numerical result."
function checkpoint_failure(inputs,output_directory,requested_chi::Integer,
        status::AbstractString,exit_code::Integer,limit_seconds::Real)
    chi = require_chi(requested_chi)
    status in ("contraction_time_limit","worker_failed") ||
        error("invalid worker failure status")
    exit_code != 0 || error("worker failure exit code must be nonzero")
    (exit_code == 124) == (status == "contraction_time_limit") ||
        error("exit 124 must be contraction_time_limit, and only exit 124 may use it")
    limit_seconds == CASE_LIMIT_SECONDS || error("case time limit differs")
    manifest = initialize_probe(inputs,output_directory)
    result_path = case_path(output_directory,chi)
    if isfile(result_path)
        result = TOML.parsefile(result_path)
        get(result,"status",nothing) in ("complete","computed_invalid") &&
            error("refusing to label a completed computation as a worker failure")
    end
    document = Dict{String,Any}(
        "schema_version"=>1, "checkpoint_kind"=>"worker_failure",
        "probe_id"=>manifest["probe_id"],
        "source_run_id"=>inputs.metadata["run_id"],
        "record_id"=>inputs.record["record_id"], "shot_index"=>1,
        "source_chi"=>inputs.metadata["chi"], "requested_chi"=>chi,
        "requested_blas_threads"=>BLAS_THREADS, "actual_blas_threads"=>0,
        "status"=>status, "valid"=>false, "exit_code"=>Int(exit_code),
        "limit_seconds"=>Float64(limit_seconds))
    progress_path = joinpath(case_directory(output_directory,chi),"progress.toml")
    if isfile(progress_path)
        saved = TOML.parsefile(progress_path)
        document["last_progress_stage"] = string(get(saved,"stage","unknown"))
        document["last_progress_shot"] = get(saved,"shot",0)
    end
    path = failure_path(output_directory,chi)
    if isfile(path)
        TOML.parsefile(path) == document || error("incompatible worker failure checkpoint")
    else
        Runner.atomic_toml(path,document)
    end
    return document
end

function validate_case_identity(case,chi,manifest,inputs)
    for (field,value) in (("probe_id",manifest["probe_id"]),
                          ("source_run_id",inputs.metadata["run_id"]),
                          ("record_id",inputs.record["record_id"]),
                          ("shot_index",1),("requested_chi",chi))
        get(case,field,nothing) == value || error("chi=$chi case has mismatched $field")
    end
    return case
end

function load_case(output_dir,chi,manifest,inputs)
    result_path = case_path(output_dir,chi)
    failed_path = failure_path(output_dir,chi)
    result = isfile(result_path) ? TOML.parsefile(result_path) : nothing
    if result isa AbstractDict && get(result,"status",nothing) in ("complete","computed_invalid")
        isfile(failed_path) && error("chi=$chi has both completed and failure checkpoints")
        return validate_case_identity(result,chi,manifest,inputs)
    end
    if isfile(failed_path)
        failure = validate_case_identity(TOML.parsefile(failed_path),chi,manifest,inputs)
        result === nothing || (failure["worker_result"] = result)
        return failure
    end
    result === nothing || return validate_case_identity(result,chi,manifest,inputs)
    return Dict{String,Any}(
        "status"=>"missing", "valid"=>false,
        "probe_id"=>manifest["probe_id"],
        "source_run_id"=>inputs.metadata["run_id"],
        "record_id"=>inputs.record["record_id"], "shot_index"=>1,
        "source_chi"=>inputs.metadata["chi"], "requested_chi"=>chi,
        "requested_blas_threads"=>BLAS_THREADS, "actual_blas_threads"=>0)
end

"Collect every requested case; invalid/missing outcomes are never discarded."
function collect_report(inputs,output_directory)
    output_dir = canonical_destination(output_directory)
    manifest = initialize_probe(inputs,output_directory)
    report_path = joinpath(output_dir,"report.toml")
    isfile(report_path) && error("diagnostic report already exists; refusing overwrite")
    cases = [load_case(output_dir,chi,manifest,inputs) for chi in CHI_CASES]
    comparisons = [compare_adjacent(cases[index],cases[index+1])
        for index in 1:length(cases)-1]
    all_computed = all(is_computed, cases)
    stable = all_computed && all(row["stable"] for row in comparisons)
    status = !all_computed ? "incomplete" : stable ? "pairwise_stable" : "needs_higher_chi"
    report = Dict{String,Any}(
        "schema_version"=>1, "probe_id"=>manifest["probe_id"],
        "status"=>status, "purpose"=>manifest["purpose"],
        "record_id"=>inputs.record["record_id"], "shot_index"=>1,
        "L"=>inputs.metadata["L"], "p"=>inputs.metadata["p"],
        "source_chi"=>inputs.metadata["chi"], "chi_cases"=>collect(CHI_CASES),
        "cutoff"=>inputs.metadata["cutoff"], "blas_threads"=>BLAS_THREADS,
        "all_cases_computed"=>all_computed,
        "computed_invalid_cases"=>count(case->get(case,"status",nothing)=="computed_invalid",cases),
        "cases"=>cases, "comparisons"=>comparisons,
        "interpretation"=>"This single-history adjacent-chi comparison diagnoses numerical convergence only. It is not an exact-ML check and is excluded from production statistics.")
    Runner.atomic_toml(report_path,report)
    return report
end

function main(args=ARGS)
    length(args) >= 3 ||
        error("expected --validate|--worker|--failure|--collect RUN_DIR OUTPUT_DIR [arguments]")
    mode,run_dir,output_dir = args[1],args[2],args[3]
    inputs = read_inputs(run_dir)
    if mode == "--validate" && length(args) == 3
        manifest = initialize_probe(inputs,output_dir)
        probe_id = manifest["probe_id"]
        println("Validated saved L5 shot 1 for chi=$(collect(CHI_CASES)); probe_id=$probe_id")
        return manifest
    elseif mode == "--worker" && length(args) == 4
        chi = parse(Int,args[4])
        result = run_worker(inputs,output_dir,chi)
        status = result["status"]
        elapsed_seconds = result["elapsed_seconds"]
        println("RESULT chi=$chi status=$status seconds=$elapsed_seconds")
        return result
    elseif mode == "--failure" && length(args) == 7
        chi = parse(Int,args[4])
        status = args[5]
        exit_code = parse(Int,args[6])
        limit_seconds = parse(Float64,args[7])
        checkpoint = checkpoint_failure(inputs,output_dir,chi,status,exit_code,limit_seconds)
        checkpoint_status = checkpoint["status"]
        println("FAILURE chi=$chi status=$checkpoint_status exit_code=$exit_code")
        return checkpoint
    elseif mode == "--collect" && length(args) == 3
        report = collect_report(inputs,output_dir)
        status = report["status"]
        report_path = joinpath(output_dir,"report.toml")
        println("CHI PROBE status=$status report=$report_path")
        return report
    end
    error("invalid arguments")
end

end # module

abspath(PROGRAM_FILE) == (@__FILE__) && TNMLChiProbe.main()
