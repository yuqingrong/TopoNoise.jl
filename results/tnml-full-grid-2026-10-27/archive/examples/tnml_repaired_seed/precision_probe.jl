"""
Fixed-history precision diagnostic; never modifies the production directory.

Usage:
  julia --project=. --threads=1 precision_probe.jl RUN_DIR OUTPUT_DIR [--validate-only] [--with-chi480]

Default cases are saved shots 1,16,19,38,39 at (chi,cutoff)=(240,1e-22)
and (240,1e-26). --with-chi480 additionally tests (480,1e-26).
The negative tolerance stays 1e-10. No clipping, resampling, exact fallback,
or production outcome replacement occurs. Exact files are read, not recomputed.
Each case and the aggregate report are atomically checkpointed in OUTPUT_DIR.
"""
module TNMLPrecisionProbe

using TOML, SHA, Dates, LinearAlgebra, TopoNoise, ITensors, ITensorMPS
include(joinpath(@__DIR__, "run.jl"))
include(joinpath(@__DIR__, "model.jl"))
const Runner = TNMLRepairedRunner
const SHOTS = [1,16,19,38,39]
const NEGATIVE_TOLERANCE = 1e-10
const TIE_TOLERANCE = 1e-9

filehash(path) = bytes2hex(sha256(read(path)))
filename(i) = "shot_" * lpad(string(i),4,'0') * ".toml"

function read_inputs(run_dir)
    metadata_path = joinpath(run_dir,"metadata.toml")
    corpus_path = joinpath(run_dir,"input-corpus.toml")
    metadata, corpus = TOML.parsefile(metadata_path), TOML.parsefile(corpus_path)
    metadata["run_id"] == Runner.identity_hash(Dict(k=>v for (k,v) in metadata if k != "run_id")) ||
        error("production metadata failed its identity check")
    metadata["exact"] === true || error("saved independent exact references are required")
    metadata["tn_uses_both_records"] === true || error("unexpected conditioning model")
    metadata["negative_tolerance"] == NEGATIVE_TOLERANCE || error("negative tolerance differs")
    metadata["L"] <= 3 || error("this fixed diagnostic requires an exact-audited small run")
    corpus["run_id"] == metadata["run_id"] || error("corpus run identity differs")
    length(corpus["records"]) == corpus["requested_shots"] == metadata["shots"] ||
        error("corpus count differs")
    root = dirname(dirname(pathof(TopoNoise)))
    for (relative,digest) in metadata["source_sha256"]
        filehash(joinpath(root,relative)) == digest || error("scientific source differs: $relative")
    end
    filehash(joinpath(root,"Project.toml")) == metadata["project_sha256"] || error("project differs")
    if haskey(metadata,"manifest_sha256")
        filehash(joinpath(root,"Manifest.toml")) == metadata["manifest_sha256"] || error("manifest differs")
    end
    selected = Dict{Int,Any}()
    for i in SHOTS
        record = corpus["records"][i]
        record["shot_index"] == i || error("corpus shot order differs")
        Runner.identity_hash(Dict(k=>v for (k,v) in record if k != "record_id")) == record["record_id"] ||
            error("saved history $i failed its identity check")
        exact_path, baseline_path = joinpath(run_dir,"exact",filename(i)), joinpath(run_dir,"shots",filename(i))
        exact, baseline = TOML.parsefile(exact_path), TOML.parsefile(baseline_path)
        for (name,doc) in (("exact",exact),("baseline",baseline))
            for (field,value) in (("run_id",metadata["run_id"]),("record_id",record["record_id"]),
                                  ("shot_index",i),("L",metadata["L"]),("p",metadata["p"]),
                                  ("chi",metadata["chi"]),("true_mask",record["true_mask"]))
                get(doc,field,nothing) == value || error("$name $i has mismatched $field")
            end
        end
        exact["status"] == "complete" && exact["valid"] === true || error("exact $i is incomplete")
        exact["algorithm"] == "independent_untruncated_dense_detector_frontier" || error("unexpected exact algorithm")
        probabilities = exact["probabilities"]
        length(probabilities) == 4 && all(isfinite,probabilities) && minimum(probabilities) >= 0 &&
            isapprox(sum(probabilities),1;atol=1e-10,rtol=0) || error("exact $i posterior is invalid")
        exact["logical_mask"] == argmax(probabilities)-1 || error("exact $i label differs from posterior")
        # Build an observed-only record directly; no truth or raw-fault fields
        # are ever passed to RepairedSeedTN.decode_record.
        observed = (seed_syndrome=BitVector(record["seed_syndrome"]),
                    final_syndrome=BitVector(record["final_syndrome"]),
                    seed_recovery=BitVector(record["seed_recovery"]),
                    recovery_mask=Int(record["recovery_mask"]))
        selected[i] = (;record,observed,exact,baseline,exact_path,baseline_path)
    end
    return (;metadata,selected,metadata_path,corpus_path)
end

function comparison(probabilities, label, valid, exact)
    reference = exact["probabilities"]
    raw_argmax = argmax(probabilities)-1
    return Dict("posterior_linf"=>maximum(abs.(probabilities .- reference)),
        "posterior_l1"=>sum(abs.(probabilities .- reference)),
        "exact_probabilities"=>reference,"exact_logical_mask"=>exact["logical_mask"],
        "candidate_argmax"=>raw_argmax,
        "candidate_argmax_tie_aware_match"=>maximum(reference)-reference[raw_argmax+1] <= TIE_TOLERANCE,
        "valid_ml_tie_aware_match"=>valid && maximum(reference)-reference[label+1] <= TIE_TOLERANCE,
        "valid_ml_strict_match"=>valid && label == exact["logical_mask"])
end

function component_report(result)
    return Dict("valid"=>result.valid,"probabilities"=>collect(result.probabilities),
        "minimum_probability"=>minimum(result.probabilities),"logical_mask"=>result.logical_mask,
        "log_syndrome"=>result.log_syndrome,"diagnostics"=>Runner.toml_value(result.diagnostics))
end

function probe_case(decoder, input, chi, cutoff, probe_id, output_dir)
    i = input.record["shot_index"]
    tag = "chi$(chi)-cutoff$(cutoff)"
    path = joinpath(output_dir,tag,filename(i))
    identity = Dict("probe_id"=>probe_id,"record_id"=>input.record["record_id"],
        "shot_index"=>i,"chi"=>chi,"cutoff"=>cutoff,"negative_tolerance"=>NEGATIVE_TOLERANCE,
        "exact_file_sha256"=>filehash(input.exact_path))
    if isfile(path)
        previous = TOML.parsefile(path)
        all(get(previous,k,nothing) == v for (k,v) in identity) || error("incompatible diagnostic checkpoint: $path")
        previous["status"] in ("complete","invalid","exception") || error("incomplete diagnostic checkpoint")
        return previous
    end
    now = time()
    Runner.atomic_toml(joinpath(output_dir,"progress.toml"),Dict("pid"=>getpid(),"stage"=>"tn",
        "shot"=>"$tag/$i","stage_started_unix"=>now,"updated_at_unix"=>now))
    println("BEGIN precision $tag shot=$i"); flush(stdout)
    output = copy(identity)
    started = time_ns()
    try
        result = RepairedSeedTN.decode_record(decoder,input.observed)
        output["valid"] = result.valid
        output["status"] = result.valid ? "complete" : "invalid"
        output["probabilities"] = collect(result.probabilities)
        output["minimum_probability"] = minimum(result.probabilities)
        output["logical_mask"] = result.logical_mask
        output["log_joint_syndrome"] = result.log_joint_syndrome
        output["seed"] = component_report(result.seed_result)
        output["transfer"] = component_report(result.transfer_result)
        output["actual_maxbond"] = max(result.seed_result.diagnostics.actual_maxbond,
                                        result.transfer_result.diagnostics.actual_maxbond)
        output["max_bond_cap_reached"] = output["actual_maxbond"] >= chi
        output["exact_comparison"] = comparison(result.probabilities,result.logical_mask,result.valid,input.exact)
        output["exact_comparison"]["log_joint_syndrome_error"] =
            result.log_joint_syndrome-input.exact["log_joint_syndrome"]
        output["baseline_probabilities"] = input.baseline["probabilities"]
        output["baseline_valid"] = input.baseline["valid"]
        output["baseline_minimum_probability"] = minimum(input.baseline["probabilities"])
        output["baseline_cutoff"] = input.baseline["cutoff"]
    catch err
        err isa InterruptException && rethrow()
        output["valid"] = false
        output["status"] = "exception"
        merge!(output,Runner.exception_details(err,catch_backtrace()))
    end
    output["elapsed_seconds"] = (time_ns()-started)/1e9
    Runner.atomic_toml(path,output)
    println("END precision $tag shot=$i status=$(output["status"]) seconds=$(output["elapsed_seconds"])")
    flush(stdout)
    return output
end

function main(args=ARGS)
    length(args) >= 2 || error("expected RUN_DIR OUTPUT_DIR [--validate-only] [--with-chi480]")
    all(arg in ("--validate-only","--with-chi480") for arg in args[3:end]) || error("unknown option")
    run_dir, output_dir = realpath(args[1]), abspath(args[2])
    (output_dir == run_dir || startswith(output_dir,run_dir*"/")) && error("diagnostics must be outside the production run")
    Threads.nthreads() == 1 || error("use --threads=1")
    BLAS.set_num_threads(1)
    inputs = read_inputs(run_dir)
    configurations = [(240,1e-22),(240,1e-26)]
    "--with-chi480" in args && push!(configurations,(480,1e-26))
    manifest = Dict("schema_version"=>1,"purpose"=>"Fixed-history cutoff diagnostic; excluded from production statistics",
        "production_directory"=>run_dir,"production_run_id"=>inputs.metadata["run_id"],
        "production_metadata_sha256"=>filehash(inputs.metadata_path),"corpus_sha256"=>filehash(inputs.corpus_path),
        "probe_script_sha256"=>filehash(@__FILE__),"selected_shots"=>SHOTS,
        "configurations"=>[Dict("chi"=>chi,"cutoff"=>cutoff) for (chi,cutoff) in configurations],
        "negative_tolerance"=>NEGATIVE_TOLERANCE,"ml_tie_tolerance"=>TIE_TOLERANCE,
        "julia_version"=>string(VERSION),"ITensors_version"=>string(Base.pkgversion(ITensors)),
        "ITensorMPS_version"=>string(Base.pkgversion(ITensorMPS)),"blas_threads"=>BLAS.get_num_threads(),
        "truth_policy"=>"Only seed syndrome, final syndrome, known seed recovery and its logical shift enter decode_record.",
        "references"=>[Dict("shot_index"=>i,"record_id"=>inputs.selected[i].record["record_id"],
            "exact_sha256"=>filehash(inputs.selected[i].exact_path),
            "baseline_sha256"=>filehash(inputs.selected[i].baseline_path)) for i in SHOTS])
    probe_id = Runner.identity_hash(manifest)
    manifest["probe_id"] = probe_id
    manifest_path = joinpath(output_dir,"probe-metadata.toml")
    if isfile(manifest_path)
        TOML.parsefile(manifest_path) == manifest || error("diagnostic metadata differs; use a fresh output directory")
    else
        Runner.atomic_toml(manifest_path,manifest)
    end
    println("Validated five saved records and independent exact references; probe_id=$probe_id");flush(stdout)
    "--validate-only" in args && return manifest
    started = time_ns()
    context = RepairedSeedTN.build_context(inputs.metadata["L"])
    setup_seconds = (time_ns()-started)/1e9
    cases = Dict{String,Any}[]
    for (chi,cutoff) in configurations
        settings = Toric3DZTNSettings(max_bond_dimension=chi,cutoff=cutoff,negative_tolerance=NEGATIVE_TOLERANCE)
        decoder = RepairedSeedTN.build_decoder(context,inputs.metadata["p"];settings)
        for i in SHOTS
            push!(cases,probe_case(decoder,inputs.selected[i],chi,cutoff,probe_id,output_dir))
            Runner.atomic_toml(joinpath(output_dir,"report.toml"),Dict("probe_id"=>probe_id,
                "expected_cases"=>length(configurations)*length(SHOTS),"completed_cases"=>length(cases),
                "setup_seconds"=>setup_seconds,"cases"=>cases))
            GC.gc()
        end
    end
    Runner.atomic_toml(joinpath(output_dir,"progress.toml"),Dict("pid"=>getpid(),"stage"=>"complete",
        "shot"=>length(cases),"stage_started_unix"=>time(),"updated_at_unix"=>time()))
    return cases
end

end # module

abspath(PROGRAM_FILE) == (@__FILE__) && TNMLPrecisionProbe.main()
