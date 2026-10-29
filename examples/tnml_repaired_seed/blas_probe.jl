"""
Read-only BLAS scaling diagnostic for one saved L=5 feasibility history.

The source run is authenticated but never modified. Three fresh Julia
processes should invoke `--worker` with BLAS thread counts 1, 4, and 8; the
Slurm wrapper does this sequentially on one node. `--collect` compares the
four- and eight-thread results with the one-thread result.

Usage:
  julia --project=. --threads=1 blas_probe.jl --validate RUN_DIR OUTPUT_DIR
  julia --project=. --threads=1 blas_probe.jl --worker RUN_DIR OUTPUT_DIR N
  julia --project=. --threads=1 blas_probe.jl --failure RUN_DIR OUTPUT_DIR N STATUS EXIT LIMIT
  julia --project=. --threads=1 blas_probe.jl --collect RUN_DIR OUTPUT_DIR
"""
module TNMLBLASProbe

using TOML, SHA, LinearAlgebra, TopoNoise, ITensors, ITensorMPS
include(joinpath(@__DIR__,"run.jl"))
include(joinpath(@__DIR__,"model.jl"))
const Runner = TNMLRepairedRunner
const Model = RepairedSeedTN

const BLAS_THREAD_CASES = (1,4,8)
const POSTERIOR_LINF_LIMIT = 1e-10
const LOGZ_LIMIT = 1e-9
const ML_TIE_TOLERANCE = 1e-9
const NEGATIVE_TOLERANCE = 1e-10

filehash(path) = bytes2hex(sha256(read(path)))

"Resolve existing symlinks, including those in the nearest existing parent."
function canonical_destination(path)
    candidate = normpath(abspath(path))
    suffix = String[]
    existing = candidate
    while !ispath(existing)
        islink(existing) && error("output path contains a broken symbolic link")
        parent = dirname(existing)
        parent == existing && error("cannot resolve output path")
        pushfirst!(suffix,basename(existing))
        existing = parent
    end
    return normpath(joinpath(realpath(existing),suffix...))
end

case_directory(output_dir,threads) = joinpath(canonical_destination(output_dir),"threads-$(threads)")
case_path(output_dir,threads) = joinpath(case_directory(output_dir,threads),"result.toml")
failure_path(output_dir,threads) = joinpath(case_directory(output_dir,threads),"failure.toml")

function repository_root()
    return realpath(dirname(dirname(pathof(TopoNoise))))
end

function path_within(path,parent)
    child = normpath(abspath(path))
    root = normpath(abspath(parent))
    return child == root || startswith(child,root*string(Base.Filesystem.path_separator))
end

function authenticated_file(root,relative,digest,description)
    relative isa AbstractString || error("$description path is invalid")
    path = normpath(joinpath(root,relative))
    path_within(path,root) || error("$description escapes the repository: $relative")
    isfile(path) || error("$description is missing: $relative")
    filehash(path) == digest || error("$description differs: $relative")
    return path
end

"""Authenticate a frozen L5 run and expose only shot one's observed fields."""
function read_inputs(run_directory)
    run_dir = realpath(run_directory)
    metadata_path = joinpath(run_dir,"metadata.toml")
    corpus_path = joinpath(run_dir,"input-corpus.toml")
    isfile(metadata_path) || error("saved metadata is missing: $metadata_path")
    isfile(corpus_path) || error("saved input corpus is missing: $corpus_path")
    metadata = TOML.parsefile(metadata_path)
    corpus = TOML.parsefile(corpus_path)
    expected_run_id = Runner.identity_hash(Dict(k=>v for (k,v) in metadata if k != "run_id"))
    get(metadata,"run_id",nothing) == expected_run_id ||
        error("saved metadata failed its identity check")
    get(metadata,"L",nothing) == 5 || error("this diagnostic requires a saved L=5 run")
    get(metadata,"exact",nothing) === false ||
        error("expected an L5 feasibility run without exact contraction")
    get(metadata,"tn_uses_both_records",nothing) === true ||
        error("saved run uses an unexpected conditioning model")
    get(metadata,"negative_tolerance",nothing) == NEGATIVE_TOLERANCE ||
        error("negative tolerance differs")
    get(metadata,"chi",0) isa Integer && metadata["chi"] > 0 || error("invalid saved chi")
    get(metadata,"cutoff",nothing) isa Real && isfinite(metadata["cutoff"]) &&
        0 <= metadata["cutoff"] < 1 || error("invalid saved cutoff")
    get(metadata,"p",nothing) isa Real && isfinite(metadata["p"]) &&
        0 <= metadata["p"] < 0.5 || error("invalid saved probability")

    root = repository_root()
    sources = get(metadata,"source_sha256",nothing)
    sources isa AbstractDict && !isempty(sources) || error("saved source hashes are missing")
    for (relative,digest) in sources
        authenticated_file(root,relative,digest,"scientific source")
    end
    authenticated_file(root,"Project.toml",metadata["project_sha256"],"project")
    if haskey(metadata,"manifest_sha256")
        authenticated_file(root,"Manifest.toml",metadata["manifest_sha256"],"manifest")
    end

    get(corpus,"run_id",nothing) == metadata["run_id"] || error("corpus run identity differs")
    records = get(corpus,"records",nothing)
    records isa AbstractVector || error("corpus records are missing")
    length(records) == get(corpus,"requested_shots",nothing) == metadata["shots"] ||
        error("corpus count differs")
    isempty(records) && error("corpus has no saved histories")
    for (index,record) in enumerate(records)
        get(record,"shot_index",nothing) == index || error("corpus shot order differs at $index")
        expected = Runner.identity_hash(Dict(k=>v for (k,v) in record if k != "record_id"))
        get(record,"record_id",nothing) == expected ||
            error("saved history $index failed its identity check")
    end
    record = records[1]
    observed = (seed_syndrome=BitVector(record["seed_syndrome"]),
        seed_recovery=BitVector(record["seed_recovery"]),
        recovery_mask=Int(record["recovery_mask"]),
        final_syndrome=BitVector(record["final_syndrome"]))
    return (;run_dir,metadata,record,observed,metadata_path,corpus_path,
        metadata_sha256=filehash(metadata_path),corpus_sha256=filehash(corpus_path))
end

function manifest_document(inputs)
    script_path = @__FILE__
    sbatch_path = joinpath(@__DIR__,"blas_probe.sbatch")
    isfile(sbatch_path) || error("BLAS probe Slurm wrapper is missing")
    document = Dict{String,Any}(
        "schema_version"=>1,
        "purpose"=>"Read-only BLAS thread scaling diagnostic; excluded from production statistics",
        "source_run_directory"=>inputs.run_dir,
        "source_run_id"=>inputs.metadata["run_id"],
        "source_metadata_sha256"=>inputs.metadata_sha256,
        "source_corpus_sha256"=>inputs.corpus_sha256,
        "record_id"=>inputs.record["record_id"], "shot_index"=>1,
        "L"=>inputs.metadata["L"], "p"=>inputs.metadata["p"],
        "chi"=>inputs.metadata["chi"], "cutoff"=>inputs.metadata["cutoff"],
        "negative_tolerance"=>inputs.metadata["negative_tolerance"],
        "blas_threads"=>collect(BLAS_THREAD_CASES),
        "timing_design"=>"One fresh process per thread count in fixed 1,4,8 order; speedups are scouting estimates without repeated-sample uncertainty.",
        "posterior_linf_limit"=>POSTERIOR_LINF_LIMIT,
        "abs_log_joint_syndrome_limit"=>LOGZ_LIMIT,
        "ml_tie_tolerance"=>ML_TIE_TOLERANCE,
        "julia_threads"=>1,
        "probe_script_sha256"=>filehash(script_path),
        "sbatch_sha256"=>filehash(sbatch_path),
        "truth_policy"=>"Only seed syndrome, final syndrome, known seed recovery and its logical shift enter decode_record. Raw faults and truth labels are excluded.")
    document["probe_id"] = Runner.identity_hash(document)
    return document
end

"""Create or verify the immutable diagnostic manifest outside the source run."""
function initialize_probe(inputs,output_directory)
    output_dir = canonical_destination(output_directory)
    path_within(output_dir,inputs.run_dir) &&
        error("diagnostic output must be outside the saved source run")
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

function component_document(result)
    return Dict("valid"=>result.valid,"probabilities"=>collect(result.probabilities),
        "logical_mask"=>result.logical_mask,"log_syndrome"=>result.log_syndrome,
        "diagnostics"=>Runner.toml_value(result.diagnostics))
end

function progress(case_dir,stage,threads; details=Dict{String,Any}())
    now = time()
    document = Dict{String,Any}("pid"=>getpid(),"stage"=>stage,"shot"=>1,
        "blas_threads"=>threads,"stage_started_unix"=>now,"updated_at_unix"=>now)
    merge!(document,details)
    Runner.atomic_toml(joinpath(case_dir,"progress.toml"),document)
    println("stage=$stage shot=1/1 BLAS=$threads"); flush(stdout)
end

"""Run exactly one authenticated history in a fresh process."""
function run_worker(inputs,output_directory,threads::Integer)
    threads in BLAS_THREAD_CASES || error("BLAS threads must be one of $(collect(BLAS_THREAD_CASES))")
    Threads.nthreads() == 1 || error("run this diagnostic with one Julia thread")
    manifest = initialize_probe(inputs,output_directory)
    case_dir = case_directory(output_directory,threads)
    result_path = case_path(output_directory,threads)
    ispath(case_dir) && error("case output already exists; use a fresh diagnostic directory")
    mkpath(case_dir)
    BLAS.set_num_threads(threads)
    actual_threads = BLAS.get_num_threads()
    actual_threads == threads || error("BLAS applied $actual_threads threads, requested $threads")
    progress(case_dir,"setup",threads)
    setup_started = time_ns()
    context = Model.build_context(inputs.metadata["L"])
    inputs.record["raw_fault_count"] == size(context.model.matrix,2) ||
        error("saved raw fault count differs from the reconstructed model")
    settings = Toric3DZTNSettings(max_bond_dimension=inputs.metadata["chi"],
        cutoff=inputs.metadata["cutoff"],negative_tolerance=NEGATIVE_TOLERANCE)
    decoder = Model.build_decoder(context,inputs.metadata["p"];settings)
    setup_seconds = (time_ns()-setup_started)/1e9
    progress(case_dir,"tn",threads)
    output = Dict{String,Any}(
        "schema_version"=>1,"probe_id"=>manifest["probe_id"],
        "source_run_id"=>inputs.metadata["run_id"],
        "record_id"=>inputs.record["record_id"],"shot_index"=>1,
        "L"=>inputs.metadata["L"],"p"=>inputs.metadata["p"],
        "chi"=>inputs.metadata["chi"],"cutoff"=>inputs.metadata["cutoff"],
        "requested_blas_threads"=>threads,"actual_blas_threads"=>actual_threads,
        "julia_threads"=>Threads.nthreads(),"blas_vendor"=>string(BLAS.vendor()),
        "blas_config"=>sprint(show,BLAS.get_config()),
        "julia_version"=>string(VERSION),
        "ITensors_version"=>string(Base.pkgversion(ITensors)),
        "ITensorMPS_version"=>string(Base.pkgversion(ITensorMPS)),
        "setup_seconds"=>setup_seconds)
    started = time_ns()
    try
        result = Model.decode_record(decoder,inputs.observed)
        output["valid"] = result.valid
        output["status"] = result.valid ? "complete" : "invalid"
        output["probabilities"] = collect(result.probabilities)
        output["logical_mask"] = result.logical_mask
        output["log_joint_syndrome"] = result.log_joint_syndrome
        output["seed_result"] = component_document(result.seed_result)
        output["transfer_result"] = component_document(result.transfer_result)
        output["actual_maxbond"] = max(result.seed_result.diagnostics.actual_maxbond,
            result.transfer_result.diagnostics.actual_maxbond)
    catch exception
        exception isa InterruptException && rethrow()
        output["status"] = "exception"; output["valid"] = false
        merge!(output,Runner.exception_details(exception,catch_backtrace()))
    end
    output["elapsed_seconds"] = (time_ns()-started)/1e9
    output["peak_rss_bytes"] = Int(Sys.maxrss())
    Runner.atomic_toml(result_path,output)
    progress(case_dir,output["status"] == "complete" ? "complete" : "blocked",threads)
    output["status"] == "complete" || error("BLAS $threads diagnostic did not complete")
    return output
end

function valid_probabilities(case)
    probabilities = get(case,"probabilities",nothing)
    return probabilities isa AbstractVector && length(probabilities) == 4 &&
        all(x->x isa Real && isfinite(x),probabilities) &&
        minimum(probabilities) >= -NEGATIVE_TOLERANCE &&
        abs(sum(probabilities)-1) <= 1e-7
end

"""Persist scheduler/subprocess failure provenance without replacing a result."""
function checkpoint_failure(inputs,output_directory,threads::Integer,status::AbstractString,
                            exit_code::Integer,limit_seconds::Real)
    threads in BLAS_THREAD_CASES || error("BLAS threads must be one of $(collect(BLAS_THREAD_CASES))")
    status in ("contraction_time_limit","worker_failed") || error("invalid worker failure status")
    exit_code != 0 || error("worker failure exit code must be nonzero")
    (exit_code == 124) == (status == "contraction_time_limit") ||
        error("exit 124 must be recorded as contraction_time_limit, and only exit 124 may use that status")
    isfinite(limit_seconds) && limit_seconds > 0 || error("failure limit must be positive and finite")
    manifest = initialize_probe(inputs,output_directory)
    result_path = case_path(output_directory,threads)
    if isfile(result_path)
        result = TOML.parsefile(result_path)
        get(result,"status",nothing) == "complete" && get(result,"valid",false) === true &&
            error("refusing to attach failure provenance to a valid completed result")
    end
    document = Dict{String,Any}(
        "schema_version"=>1,"checkpoint_kind"=>"worker_failure",
        "probe_id"=>manifest["probe_id"],"source_run_id"=>inputs.metadata["run_id"],
        "record_id"=>inputs.record["record_id"],"shot_index"=>1,
        "requested_blas_threads"=>threads,"actual_blas_threads"=>0,
        "status"=>status,"valid"=>false,"exit_code"=>Int(exit_code),
        "limit_seconds"=>Float64(limit_seconds))
    progress_path = joinpath(case_directory(output_directory,threads),"progress.toml")
    if isfile(progress_path)
        saved_progress = TOML.parsefile(progress_path)
        document["last_progress_stage"] = string(get(saved_progress,"stage","unknown"))
        document["last_progress_shot"] = get(saved_progress,"shot",0)
    end
    path = failure_path(output_directory,threads)
    if isfile(path)
        TOML.parsefile(path) == document || error("incompatible worker failure checkpoint")
    else
        Runner.atomic_toml(path,document)
    end
    return document
end

"""Compare one threaded result with the one-thread numerical reference."""
function compare_case(baseline,candidate)
    base_probabilities = get(baseline,"probabilities",Float64[])
    candidate_probabilities = get(candidate,"probabilities",Float64[])
    numerical_inputs = valid_probabilities(baseline) && valid_probabilities(candidate)
    linf = numerical_inputs ? maximum(abs.(base_probabilities .- candidate_probabilities)) : Inf
    l1 = numerical_inputs ? sum(abs.(base_probabilities .- candidate_probabilities)) : Inf
    base_logz = get(baseline,"log_joint_syndrome",NaN)
    candidate_logz = get(candidate,"log_joint_syndrome",NaN)
    log_error = base_logz isa Real && candidate_logz isa Real && isfinite(base_logz) &&
        isfinite(candidate_logz) ? abs(base_logz-candidate_logz) : Inf
    base_label = get(baseline,"logical_mask",-1)
    candidate_label = get(candidate,"logical_mask",-1)
    labels_valid = base_label isa Integer && candidate_label isa Integer &&
        base_label in 0:3 && candidate_label in 0:3
    tie_match = numerical_inputs && labels_valid &&
        maximum(base_probabilities)-base_probabilities[candidate_label+1] <= ML_TIE_TOLERANCE &&
        maximum(candidate_probabilities)-candidate_probabilities[base_label+1] <= ML_TIE_TOLERANCE
    requested = get(candidate,"requested_blas_threads",-1)
    actual = get(candidate,"actual_blas_threads",-2)
    thread_count_match = requested == actual
    complete = get(baseline,"status",nothing) == "complete" && get(baseline,"valid",false) === true &&
        get(candidate,"status",nothing) == "complete" && get(candidate,"valid",false) === true
    passed = complete && numerical_inputs && thread_count_match && tie_match &&
        linf <= POSTERIOR_LINF_LIMIT && log_error <= LOGZ_LIMIT
    baseline_elapsed = get(baseline,"elapsed_seconds",NaN)
    candidate_elapsed = get(candidate,"elapsed_seconds",NaN)
    speedup = baseline_elapsed isa Real && candidate_elapsed isa Real &&
        isfinite(baseline_elapsed) && isfinite(candidate_elapsed) && candidate_elapsed > 0 ?
        baseline_elapsed/candidate_elapsed : NaN
    return Dict{String,Any}("blas_threads"=>requested,"passed"=>passed,
        "complete"=>complete,"actual_thread_count_matches"=>thread_count_match,
        "posterior_linf"=>linf,"posterior_l1"=>l1,
        "abs_log_joint_syndrome_error"=>log_error,
        "tie_aware_label_match"=>tie_match,"speedup_vs_one_thread"=>speedup,
        "core_seconds"=>(candidate_elapsed isa Real ? requested*candidate_elapsed : NaN))
end

function load_case(output_dir,threads,manifest,inputs)
    path = case_path(output_dir,threads)
    failure_checkpoint_path = failure_path(output_dir,threads)
    case = isfile(path) ? TOML.parsefile(path) : nothing
    if isfile(failure_checkpoint_path) &&
            !(case isa AbstractDict && get(case,"status",nothing) == "complete" &&
              get(case,"valid",false) === true)
        failure = TOML.parsefile(failure_checkpoint_path)
        for (field,value) in (("probe_id",manifest["probe_id"]),
                              ("record_id",inputs.record["record_id"]),("shot_index",1),
                              ("requested_blas_threads",threads))
            get(failure,field,nothing) == value ||
                error("BLAS $threads failure checkpoint has mismatched $field")
        end
        case === nothing || (failure["worker_result"] = case)
        return failure
    end
    case === nothing && return Dict{String,Any}(
        "status"=>"missing","valid"=>false,"requested_blas_threads"=>threads,
        "actual_blas_threads"=>0,"record_id"=>inputs.record["record_id"],"shot_index"=>1)
    for (field,value) in (("probe_id",manifest["probe_id"]),
                          ("record_id",inputs.record["record_id"]),("shot_index",1),
                          ("requested_blas_threads",threads))
        get(case,field,nothing) == value || error("BLAS $threads case has mismatched $field")
    end
    return case
end

"""Collect all configurations without dropping missing or failed cases."""
function collect_report(inputs,output_directory)
    output_dir = canonical_destination(output_directory)
    manifest = initialize_probe(inputs,output_directory)
    report_path = joinpath(output_dir,"report.toml")
    isfile(report_path) && error("diagnostic report already exists; refusing overwrite")
    cases = [load_case(output_dir,threads,manifest,inputs) for threads in BLAS_THREAD_CASES]
    baseline = cases[1]
    comparisons = [compare_case(baseline,candidate) for candidate in cases[2:end]]
    baseline_ok = get(baseline,"status",nothing) == "complete" &&
        get(baseline,"valid",false) === true &&
        get(baseline,"actual_blas_threads",nothing) == 1 && valid_probabilities(baseline)
    passed = baseline_ok && all(row["passed"] for row in comparisons)
    report = Dict{String,Any}(
        "schema_version"=>1,"probe_id"=>manifest["probe_id"],
        "status"=>passed ? "passed" : "failed",
        "purpose"=>manifest["purpose"],"record_id"=>inputs.record["record_id"],
        "shot_index"=>1,"L"=>inputs.metadata["L"],"p"=>inputs.metadata["p"],
        "chi"=>inputs.metadata["chi"],"cutoff"=>inputs.metadata["cutoff"],
        "baseline_complete"=>baseline_ok,"cases"=>cases,"comparisons"=>comparisons,
        "interpretation"=>"This one-shot fixed-order probe gives scouting estimates of BLAS scaling and checks cross-thread numerical consistency. It does not quantify timing variance, TN convergence, or exact ML accuracy.")
    Runner.atomic_toml(report_path,report)
    return report
end

function require_passed_report(report)
    get(report,"status",nothing) == "passed" ||
        error("BLAS diagnostic failed its completion or numerical-consistency checks")
    return report
end

function main(args=ARGS)
    length(args) >= 3 || error("expected --validate|--worker|--collect RUN_DIR OUTPUT_DIR [BLAS_THREADS]")
    mode,run_dir,output_dir = args[1],args[2],args[3]
    inputs = read_inputs(run_dir)
    if mode == "--validate" && length(args) == 3
        manifest = initialize_probe(inputs,output_dir)
        probe_id = manifest["probe_id"]
        println("Validated saved L5 shot 1; probe_id=$probe_id")
        return manifest
    elseif mode == "--worker" && length(args) == 4
        threads = parse(Int,args[4])
        result = run_worker(inputs,output_dir,threads)
        elapsed_seconds = result["elapsed_seconds"]
        peak_rss_bytes = result["peak_rss_bytes"]
        println("RESULT BLAS=$threads seconds=$elapsed_seconds peak_rss_bytes=$peak_rss_bytes")
        return result
    elseif mode == "--failure" && length(args) == 7
        threads = parse(Int,args[4])
        status = args[5]
        exit_code = parse(Int,args[6])
        limit_seconds = parse(Float64,args[7])
        checkpoint = checkpoint_failure(inputs,output_dir,threads,status,exit_code,limit_seconds)
        checkpoint_status = checkpoint["status"]
        println("FAILURE BLAS=$threads status=$checkpoint_status exit_code=$exit_code")
        return checkpoint
    elseif mode == "--collect" && length(args) == 3
        report = collect_report(inputs,output_dir)
        status = report["status"]
        report_path = joinpath(output_dir,"report.toml")
        println("BLAS PROBE status=$status report=$report_path")
        return require_passed_report(report)
    end
    error("invalid arguments")
end

end # module

abspath(PROGRAM_FILE) == (@__FILE__) && TNMLBLASProbe.main()
