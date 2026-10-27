"""
Checkpointed matched-preparation experiment. All histories are persisted before
any exact or TN contraction. Invoke with --threads=1; the runner fixes BLAS to one
thread. A separate process may enforce resource limits using progress.toml.

Usage: julia --project=. --threads=1 examples/tnml_repaired_seed/run.jl \
    L p shots chi seed output_dir [exact=true|false]
Exact comparison defaults to false and is restricted to L <= 3, width <= 22.
"""
module TNMLRepairedRunner

using Random, TOML, SHA, Dates, LinearAlgebra

const SCHEMA_VERSION = 1
const CUTOFF = 1e-18
const NEGATIVE_TOLERANCE = 1e-10
const EXACT_MAX_WIDTH = 22
const ML_TIE_TOLERANCE = 1e-9

function parse_config(args)
    6 <= length(args) <= 7 || throw(ArgumentError(
        "expected L p shots chi seed output_dir [exact=true|false]"))
    L, p, shots, chi, seed = parse(Int,args[1]), parse(Float64,args[2]),
        parse(Int,args[3]), parse(Int,args[4]), parse(Int,args[5])
    exact_arg = length(args) == 7 ? replace(args[7], "exact="=>"") : "false"
    exact_arg in ("true","false") || throw(ArgumentError("exact must be true or false"))
    exact = exact_arg == "true"
    L >= 2 || throw(ArgumentError("L must be at least 2"))
    isfinite(p) && 0 <= p < 0.5 || throw(ArgumentError("p must be finite and in [0,0.5)"))
    shots > 0 || throw(ArgumentError("shots must be positive"))
    chi > 0 || throw(ArgumentError("chi must be positive"))
    seed >= 0 || throw(ArgumentError("seed must be nonnegative"))
    exact && L > 3 && throw(ArgumentError("exact reference is restricted to L <= 3"))
    return (;L,p,shots,chi,seed,output_dir=abspath(args[6]),exact)
end

"Convert diagnostics to portable TOML values while preserving signed floats."
toml_value(x::NamedTuple) = Dict(string(k)=>toml_value(v) for (k,v) in pairs(x))
toml_value(x::AbstractDict) = Dict(string(k)=>toml_value(v) for (k,v) in pairs(x))
toml_value(x::AbstractArray) = [toml_value(v) for v in x]
toml_value(x::Tuple) = [toml_value(v) for v in x]
toml_value(x::Symbol) = string(x)
toml_value(x::Nothing) = "nothing"
toml_value(x) = x

function identity_hash(doc)
    io = IOBuffer()
    TOML.print(io, toml_value(doc); sorted=true)
    return bytes2hex(sha256(take!(io)))
end

"Replace a checkpoint atomically after flushing its complete TOML payload."
function atomic_toml(path, doc)
    mkpath(dirname(path))
    temporary, io = mktemp(dirname(path))
    try
        TOML.print(io, toml_value(doc); sorted=true)
        flush(io)
        close(io)
        mv(temporary, path; force=true)
    finally
        isopen(io) && close(io)
        isfile(temporary) && rm(temporary)
    end
    return path
end

function resume_document(path, run_id, record_id, shot_index)
    doc = TOML.parsefile(path)
    get(doc,"run_id",nothing) == run_id || throw(ArgumentError(
        "refusing incompatible run/source checkpoint: $path"))
    get(doc,"record_id",nothing) == record_id || throw(ArgumentError(
        "refusing incompatible saved history: $path"))
    get(doc,"shot_index",nothing) == shot_index || throw(ArgumentError(
        "refusing mismatched shot index: $path"))
    get(doc,"status",nothing) in ("complete","invalid","exception") ||
        throw(ArgumentError("refusing incomplete or unrecognized checkpoint: $path"))
    return doc
end

function record_document(record, shot_index)
    doc = Dict{String,Any}(
        "shot_index"=>shot_index,
        "raw_fault_count"=>length(record.raw_faults),
        "raw_fault_indices"=>findall(record.raw_faults),
        "seed_syndrome"=>collect(record.seed_syndrome),
        "seed_recovery"=>collect(record.seed_recovery),
        "recovery_mask"=>Int(record.recovery_mask),
        "final_syndrome"=>collect(record.final_syndrome),
        "true_mask"=>Int(record.true_mask),
        "mwpm_failure"=>record.mwpm_failure,
        "raw_seed_mask"=>Int(record.raw_seed_mask),
        "transfer_mask"=>Int(record.transfer_mask))
    doc["record_id"] = identity_hash(doc)
    return doc
end

function restore_record(doc)
    raw_faults = falses(doc["raw_fault_count"])
    raw_faults[doc["raw_fault_indices"]] .= true
    return (;raw_faults, seed_syndrome=BitVector(doc["seed_syndrome"]),
        seed_recovery=BitVector(doc["seed_recovery"]), recovery_mask=Int(doc["recovery_mask"]),
        final_syndrome=BitVector(doc["final_syndrome"]), true_mask=Int(doc["true_mask"]),
        mwpm_failure=Bool(doc["mwpm_failure"]), raw_seed_mask=Int(doc["raw_seed_mask"]),
        transfer_mask=Int(doc["transfer_mask"]))
end

observed_record(record) = (;record.seed_syndrome,record.seed_recovery,
                           record.recovery_mask,record.final_syndrome)

function load_backend!()
    isdefined(@__MODULE__, :RepairedSeedTN) && return
    @eval using TopoNoise, ITensors, ITensorMPS
    include(joinpath(@__DIR__, "model.jl"))
    include(joinpath(@__DIR__, "dense_frontier.jl"))
end

function run_metadata(cfg)
    root = dirname(dirname(pathof(TopoNoise)))
    project = TOML.parsefile(joinpath(root,"Project.toml"))
    sources = Dict{String,String}()
    for (dir, _, files) in walkdir(joinpath(root,"src")), file in sort(files)
        endswith(file,".jl") || continue
        path = joinpath(dir,file)
        sources[relpath(path,root)] = bytes2hex(sha256(read(path)))
    end
    for file in ("model.jl","run.jl","dense_frontier.jl")
        sources["examples/tnml_repaired_seed/"*file] =
            bytes2hex(sha256(read(joinpath(@__DIR__,file))))
    end
    doc = Dict{String,Any}(
        "schema_version"=>SCHEMA_VERSION,
        "L"=>cfg.L,"p"=>cfg.p,"shots"=>cfg.shots,"chi"=>cfg.chi,
        "cutoff"=>CUTOFF,"negative_tolerance"=>NEGATIVE_TOLERANCE,
        "exact"=>cfg.exact,"exact_max_width"=>EXACT_MAX_WIDTH,
        "ml_tie_tolerance"=>ML_TIE_TOLERANCE,
        "rng_algorithm"=>"Random.Xoshiro","rng_seed"=>cfg.seed,
        "physical_noise"=>"Independent Z faults with probability p on every seed-H output and every H/CNOT output qubit in the diagonal H+6CNOT circuit.",
        "preparation"=>"Noisy seed-H outputs; complete initial seed MWPM syndrome repair; diagonal H+6CNOT transfer; no later feedback.",
        "tn_uses_both_records"=>true,
        "tn_conditioning"=>"Measured preparation syndrome and measured final syndrome, with the known seed-recovery logical shift.",
        "mwpm_conditioning"=>"Final syndrome only, after the same physical initial seed repair.",
        "logical_mask_order"=>["(0,0)","(1,0)","(0,1)","(1,1)"],
        "truth_definition"=>"Logical mask of the full raw fault frame after seed repair, before final decoder recovery; TNfailure is true_mask != logical_mask.",
        "invalid_policy"=>"Invalid, exception, and absent shot outputs are unknown outcomes; retain the requested-history denominator and record attempted and completed TN counts separately.",
        "posterior_cache"=>"disabled; each Monte Carlo history has its own output",
        "TopoNoise_version"=>string(project["version"]),
        "TopoNoise_project_path"=>root,"active_project_path"=>something(Base.active_project(),""),
        "julia_version"=>string(VERSION),"ITensors_version"=>string(Base.pkgversion(ITensors)),
        "ITensorMPS_version"=>string(Base.pkgversion(ITensorMPS)),
        "PythonCall_version"=>string(Base.pkgversion(TopoNoise.PythonCall)),
        "python_version"=>TopoNoise.PythonCall.pyconvert(String,TopoNoise.PythonCall.pyimport("sys").version),
        "python_executable"=>TopoNoise.PythonCall.pyconvert(String,TopoNoise.PythonCall.pyimport("sys").executable),
        "pymatching_version"=>TopoNoise._toric3d_module_version(TopoNoise._load_pymatching(),"pymatching"),
        "julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),
        "source_sha256"=>sources)
    manifest = joinpath(root,"Manifest.toml")
    isfile(manifest) && (doc["manifest_sha256"] = bytes2hex(sha256(read(manifest))))
    doc["project_sha256"] = bytes2hex(sha256(read(joinpath(root,"Project.toml"))))
    doc["run_id"] = identity_hash(doc)
    return doc
end

function prepare_metadata(cfg, metadata)
    path = joinpath(cfg.output_dir,"metadata.toml")
    if isfile(path)
        existing = TOML.parsefile(path)
        existing == metadata || throw(ArgumentError(
            "output directory has different configuration, sources, or runtime; choose a fresh output_dir: $(cfg.output_dir)"))
    else
        leftovers = ["input-corpus.toml","progress.toml","summary.toml","shots","exact"]
        any(name->ispath(joinpath(cfg.output_dir,name)),leftovers) && throw(ArgumentError(
            "output directory has prior checkpoints but no matching metadata; refusing to overwrite"))
        atomic_toml(path,metadata)
    end
    return metadata["run_id"]
end

function progress(cfg, run_id, stage, shot; details=Dict{String,Any}())
    now = time()
    doc = Dict{String,Any}("run_id"=>run_id,"stage"=>stage,"shot"=>shot,
        "stage_started_unix"=>now,"updated_at_unix"=>now,"pid"=>getpid(),
        "requested_shots"=>cfg.shots,"L"=>cfg.L,"p"=>cfg.p)
    merge!(doc,details)
    atomic_toml(joinpath(cfg.output_dir,"progress.toml"),doc)
    println("stage=$stage shot=$shot/$(cfg.shots) L=$(cfg.L) p=$(cfg.p)")
    flush(stdout)
    return nothing
end

function input_corpus(cfg, ctx, run_id)
    path = joinpath(cfg.output_dir,"input-corpus.toml")
    if isfile(path)
        corpus = TOML.parsefile(path)
        get(corpus,"run_id",nothing) == run_id || throw(ArgumentError("corpus run identity differs"))
        records = corpus["records"]
        length(records) == cfg.shots || throw(ArgumentError("corpus shot count differs"))
        for (i,doc) in enumerate(records)
            get(doc,"shot_index",0) == i || throw(ArgumentError("corpus shot order differs"))
            identity_hash(Dict(k=>v for (k,v) in doc if k != "record_id")) == doc["record_id"] ||
                throw(ArgumentError("saved history $i failed its identity check"))
            doc["raw_fault_count"] == size(ctx.model.matrix,2) || throw(ArgumentError("raw fault count differs"))
        end
        return records
    end
    any(name->isdir(joinpath(cfg.output_dir,name)) &&
        !isempty(readdir(joinpath(cfg.output_dir,name))), ("shots","exact")) &&
        throw(ArgumentError("saved outputs exist without their input corpus; refusing regeneration"))
    progress(cfg,run_id,"sampling",0)
    rng = Xoshiro(cfg.seed)
    started = time_ns()
    records = [record_document(RepairedSeedTN.sample_record(ctx,cfg.p,rng),i) for i in 1:cfg.shots]
    atomic_toml(path,Dict("run_id"=>run_id,"rng_seed"=>cfg.seed,
        "requested_shots"=>cfg.shots,"sampling_elapsed_seconds"=>(time_ns()-started)/1e9,
        "records"=>records))
    return records
end

shot_name(i) = "shot_" * lpad(string(i),4,'0') * ".toml"

function shot_base(cfg, run_id, doc)
    return Dict{String,Any}("schema_version"=>SCHEMA_VERSION,"run_id"=>run_id,
        "record_id"=>doc["record_id"],"shot_index"=>doc["shot_index"],
        "L"=>cfg.L,"p"=>cfg.p,"chi"=>cfg.chi,"cutoff"=>CUTOFF,
        "true_mask"=>doc["true_mask"],"mwpm_failure"=>doc["mwpm_failure"],
        "input_corpus"=>"input-corpus.toml")
end

function exception_details(err, backtrace)
    Dict("exception_type"=>string(typeof(err)),"exception_message"=>sprint(showerror,err),
         "exception_trace"=>sprint(showerror,err,backtrace))
end

"Independent exact seed/transfer contraction and explicit logical XOR sum."
function exact_reference(cfg, ctx, record, run_id, doc)
    path = joinpath(cfg.output_dir,"exact",shot_name(doc["shot_index"]))
    isfile(path) && return resume_document(path,run_id,doc["record_id"],doc["shot_index"])
    progress(cfg,run_id,"exact",doc["shot_index"])
    result = shot_base(cfg,run_id,doc)
    result["algorithm"] = "independent_untruncated_dense_detector_frontier"
    result["max_width_guard"] = EXACT_MAX_WIDTH
    started = time_ns()
    try
        component_started = time_ns()
        seed = TNDecoder3DDenseFrontier.dense_exact_probabilities(
            ctx.seed_topology,record.seed_syndrome,cfg.p;max_width=EXACT_MAX_WIDTH)
        result["seed_elapsed_seconds"] = (time_ns()-component_started)/1e9
        result["seed"] = toml_value(seed)
        result["seed_joint_weights"] = seed.probabilities .* exp(seed.log_marginal_probability)
        component_started = time_ns()
        transfer = TNDecoder3DDenseFrontier.dense_exact_probabilities(
            ctx.transfer_topology,record.final_syndrome,cfg.p;max_width=EXACT_MAX_WIDTH)
        result["transfer_elapsed_seconds"] = (time_ns()-component_started)/1e9
        result["transfer"] = toml_value(transfer)
        result["transfer_joint_weights"] = transfer.probabilities .* exp(transfer.log_marginal_probability)
        probabilities = zeros(4)
        # Deliberately independent of the wrapper's combine_posteriors helper.
        for raw_seed in 0:3, transferred in 0:3
            corrected = xor(xor(raw_seed,record.recovery_mask),transferred)
            probabilities[corrected+1] += seed.probabilities[raw_seed+1] * transfer.probabilities[transferred+1]
        end
        probabilities ./= sum(probabilities)
        result["status"] = "complete"
        result["valid"] = true
        result["probabilities"] = probabilities
        result["logical_mask"] = argmax(probabilities)-1
        result["log_joint_syndrome"] = seed.log_marginal_probability + transfer.log_marginal_probability
        result["joint_weights"] = probabilities .* exp(result["log_joint_syndrome"])
        result["exact_failure"] = record.true_mask != result["logical_mask"]
    catch err
        err isa InterruptException && rethrow()
        result["status"] = "exception"
        result["valid"] = false
        merge!(result,exception_details(err,catch_backtrace()))
    end
    result["elapsed_seconds"] = (time_ns()-started)/1e9
    atomic_toml(path,result)
    return result
end

function compare_exact(result, exact)
    probabilities = exact["probabilities"]
    comparison = Dict{String,Any}("posterior_l1"=>sum(abs.(result.probabilities .- probabilities)),
        "posterior_linf"=>maximum(abs.(result.probabilities .- probabilities)),
        "ml_tie_tolerance"=>ML_TIE_TOLERANCE,
        "exact_maximizers"=>findall(x->maximum(probabilities)-x <= ML_TIE_TOLERANCE,probabilities).-1)
    if result.valid
        comparison["ml_strict_match"] = result.logical_mask == exact["logical_mask"]
        comparison["ml_tie_aware_match"] = result.logical_mask in comparison["exact_maximizers"]
    end
    return comparison
end

function decode_shot(cfg, decoder, record, run_id, doc, exact)
    progress(cfg,run_id,"tn",doc["shot_index"])
    output = shot_base(cfg,run_id,doc)
    output["exact_requested"] = cfg.exact
    if exact !== nothing
        output["exact_reference_path"] = joinpath("exact",shot_name(doc["shot_index"]))
        output["exact_status"] = exact["status"]
    end
    started = time_ns()
    try
        result = RepairedSeedTN.decode_record(decoder,observed_record(record))
        output["valid"] = result.valid
        output["status"] = result.valid ? "complete" : "invalid"
        output["probabilities"] = collect(result.probabilities)
        output["raw_probabilities"] = collect(result.probabilities)
        output["logical_mask"] = result.logical_mask
        output["log_joint_syndrome"] = result.log_joint_syndrome
        output["seed_result"] = toml_value(result.seed_result)
        output["transfer_result"] = toml_value(result.transfer_result)
        # true_mask is after the physical seed repair and BEFORE final MWPM.
        # Do not XOR in final-MWPM recovery: that would count it twice.
        result.valid && (output["TNfailure"] = record.true_mask != result.logical_mask)
        if exact !== nothing && exact["status"] == "complete"
            output["exact_comparison"] = compare_exact(result,exact)
            output["exact_comparison"]["log_joint_syndrome_error"] =
                result.log_joint_syndrome - exact["log_joint_syndrome"]
        end
    catch err
        err isa InterruptException && rethrow()
        output["status"] = "exception"
        output["valid"] = false
        delete!(output,"TNfailure")
        merge!(output,exception_details(err,catch_backtrace()))
    end
    output["elapsed_seconds"] = (time_ns()-started)/1e9
    atomic_toml(joinpath(cfg.output_dir,"shots",shot_name(doc["shot_index"])),output)
    return output
end

function summarize(cfg, run_id, records, outputs; attempted_tn_shots=length(outputs))
    valid = count(doc->doc["status"] == "complete" && doc["valid"],outputs)
    failures = count(doc->doc["status"] == "complete" && doc["valid"] &&
        get(doc,"TNfailure",false),outputs)
    unknown = cfg.shots-valid
    summary = Dict{String,Any}("run_id"=>run_id,"L"=>cfg.L,"p"=>cfg.p,
        "requested_shots"=>cfg.shots,"saved_histories"=>length(records),
        "saved_outputs"=>length(outputs),"valid_shots"=>valid,"tn_failures"=>failures,
        "attempted_tn_shots"=>attempted_tn_shots,
        "completed_tn_shots"=>count(d->d["status"] in ("complete","invalid"),outputs),
        "unknown_shots"=>unknown,"invalid_shots"=>count(d->d["status"] == "invalid",outputs),
        "exception_shots"=>count(d->d["status"] == "exception",outputs),
        "mwpm_failures"=>count(d->d["mwpm_failure"],records),
        "tn_failure_identification_interval"=>[failures/cfg.shots,(failures+unknown)/cfg.shots],
        "updated_at_utc"=>string(now(UTC)))
    atomic_toml(joinpath(cfg.output_dir,"summary.toml"),summary)
    return summary
end

function run_loaded(cfg)
    Threads.nthreads() == 1 || throw(ArgumentError("run this experiment with --threads=1"))
    BLAS.set_num_threads(1)
    metadata = run_metadata(cfg)
    run_id = prepare_metadata(cfg,metadata)
    progress(cfg,run_id,"setup",0)
    ctx = RepairedSeedTN.build_context(cfg.L)
    records = input_corpus(cfg,ctx,run_id)
    outputs = Dict{String,Any}[]
    resumed = Set{Int}()
    for doc in records
        i = doc["shot_index"]
        path = joinpath(cfg.output_dir,"shots",shot_name(i))
        if isfile(path)
            push!(outputs,resume_document(path,run_id,doc["record_id"],i))
            push!(resumed,i)
        end
    end
    # Persist an all-unknown initial summary even if the first contraction is
    # killed; a resumed run retains all already validated output counts.
    initial_summary = summarize(cfg,run_id,records,outputs)
    if length(outputs) == cfg.shots
        progress(cfg,run_id,"complete",cfg.shots)
        return initial_summary
    end
    # Compile the decoder only after saving the complete raw-history corpus.
    progress(cfg,run_id,"setup",0;details=Dict("detail"=>"building TN decoder"))
    settings = TopoNoise.Toric3DZTNSettings(max_bond_dimension=cfg.chi,
        cutoff=CUTOFF,negative_tolerance=NEGATIVE_TOLERANCE)
    decoder = RepairedSeedTN.build_decoder(ctx,cfg.p;settings)
    for doc in records
        i = doc["shot_index"]
        i in resumed && continue
        record = restore_record(doc)
        exact = cfg.exact ? exact_reference(cfg,ctx,record,run_id,doc) : nothing
        summarize(cfg,run_id,records,outputs;attempted_tn_shots=length(outputs)+1)
        push!(outputs,decode_shot(cfg,decoder,record,run_id,doc,exact))
        summarize(cfg,run_id,records,outputs)
        GC.gc()
    end
    summary = summarize(cfg,run_id,records,outputs)
    progress(cfg,run_id,"complete",cfg.shots)
    return summary
end

function main(args=ARGS)
    if any(arg->arg in ("-h","--help"),args)
        println("Usage: julia --project=. --threads=1 examples/tnml_repaired_seed/run.jl L p shots chi seed output_dir [exact=true|false]")
        println("Exact defaults to false; L <= 3 only. Cutoff=1e-18, negative tolerance=1e-10, BLAS threads=1.")
        println("Resume requires identical configuration, sources, runtime, and saved history identities.")
        return nothing
    end
    cfg = parse_config(args)
    load_backend!()
    return Base.invokelatest(run_loaded,cfg)
end

end # module

abspath(PROGRAM_FILE) == (@__FILE__) && TNMLRepairedRunner.main()
