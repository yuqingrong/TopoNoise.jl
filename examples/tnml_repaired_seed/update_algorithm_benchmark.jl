"""
Experimental comparison of finite-prefix MPO×MPS update algorithms.

This file is deliberately outside `src/`: it does not change the production
decoder.  The `current` path reproduces the production product-then-truncate
update.  The `zipup` path uses ITensorMPS' public MPO×MPS implementation on an
isolated prefix while retaining the live bond to the untouched suffix.
ITensorMPS' density-matrix implementation replaces that external boundary
index, so it is rejected whenever a live suffix exists.
"""
module TNMLUpdateAlgorithmBenchmark

using LinearAlgebra, TOML, ITensors, ITensorMPS, TopoNoise
include(joinpath(@__DIR__, "model.jl"))
include(joinpath(@__DIR__, "blas_probe.jl"))
const Model = RepairedSeedTN
const InputProbe = TNMLBLASProbe

export UnsupportedFinitePrefix, apply_update, compare_small_decoder,
       collect_replay, combine_component_results, decode_component,
       initialize_replay, initial_boundary, run_replay_worker

struct UnsupportedFinitePrefix <: Exception
    algorithm::Symbol
    message::String
end
Base.showerror(io::IO, error::UnsupportedFinitePrefix) =
    print(io, error.algorithm, ": ", error.message)

function algorithm_symbol(algorithm)
    result = algorithm isa Symbol ? algorithm : Symbol(algorithm)
    result in (:current, :zipup, :densitymatrix) ||
        throw(ArgumentError("algorithm must be current, zipup, or densitymatrix"))
    return result
end

function initial_boundary(decoder)
    n = decoder.detector_count
    sites = vcat(Index[Index(4, "Logical")],
        [Index(2, "Detector,v=$vertex") for vertex in 1:n])
    return ITensorMPS.MPS(sites, fill(1, n + 1))
end

function append_suffix(prefix, original, support)
    tensors = ITensor[prefix[k] for k in eachindex(prefix)]
    append!(tensors, (original[k] for k in (support + 1):length(original)))
    return ITensorMPS.MPS(tensors; ortho_lims=1:support)
end

function live_boundary(psi, support)
    support < length(psi) || return nothing
    return ITensorMPS.linkind(psi, support)
end

function update_stats(algorithm, elapsed, svd_events, callback_events,
        discarded_weight, maximum_local_truncation_error)
    return (; algorithm=String(algorithm), elapsed_seconds=elapsed,
        svd_events, callback_events, discarded_weight,
        maximum_local_truncation_error, peak_rss_bytes=Int(Sys.maxrss()))
end

"""
Apply one finite-support star MPO and compress its prefix.

`zipup` follows the ITensorMPS recommendation: the first sweep uses
`2maxdim` and `cutoff/10`, then the final optimal sweep uses the requested
limits.  Its first-sweep truncation errors are not exposed by the package API;
`callback_events` and discarded-weight fields therefore cover the final sweep,
while `svd_events` records both known SVD sweeps.
"""
function apply_update(psi, star, algorithm; maxdim::Integer, cutoff::Real)
    alg = algorithm_symbol(algorithm)
    maxdim > 0 || throw(ArgumentError("maxdim must be positive"))
    isfinite(cutoff) && 0 <= cutoff < 1 ||
        throw(ArgumentError("cutoff must be finite and in [0,1)"))
    support = length(star)
    1 <= support <= length(psi) ||
        throw(DimensionMismatch("star support must be a nonempty MPS prefix"))

    callback_events = Ref(0)
    discarded_weight = Ref(0.0)
    maximum_local = Ref(0.0)
    callback = function (; link, truncation_error)
        value = Float64(truncation_error)
        callback_events[] += 1
        discarded_weight[] += value
        maximum_local[] = max(maximum_local[], value)
        return nothing
    end

    started = time_ns()
    result = if alg === :current
        updated = TopoNoise._toric3d_z_tn_apply_prefix(psi, star)
        updated = ITensors.replaceprime(updated, 1 => 0)
        ITensorMPS.truncate!(updated; site_range=1:support,
            maxdim, cutoff, callback)
        updated
    elseif alg === :zipup
        boundary = live_boundary(psi, support)
        prefix = ITensorMPS.MPS(
            ITensor[psi[k] for k in 1:support]; ortho_lims=1:support)
        first_maxdim = maxdim > typemax(Int) ÷ 2 ? typemax(Int) : 2Int(maxdim)
        updated_prefix = ITensorMPS.apply(star, prefix; alg="zipup",
            maxdim=first_maxdim, cutoff=Float64(cutoff) / 10,
            truncate_kwargs=(; maxdim=Int(maxdim), cutoff=Float64(cutoff), callback))
        if !isnothing(boundary) && !any(==(boundary), inds(updated_prefix[support]))
            throw(UnsupportedFinitePrefix(:zipup,
                "ITensorMPS did not preserve the live suffix bond"))
        end
        append_suffix(updated_prefix, psi, support)
    else
        boundary = live_boundary(psi, support)
        if !isnothing(boundary)
            throw(UnsupportedFinitePrefix(:densitymatrix,
                "ITensorMPS densitymatrix replaces the live suffix bond"))
        end
        prefix = ITensorMPS.MPS(
            ITensor[psi[k] for k in 1:support]; ortho_lims=1:support)
        ITensorMPS.apply(star, prefix; alg="densitymatrix",
            maxdim=Int(maxdim), cutoff=Float64(cutoff))
    end
    elapsed = (time_ns() - started) / 1e9
    svd_events = if alg === :zipup
        2 * max(0, support - 1)
    else
        max(0, support - 1)
    end
    return result, update_stats(alg, elapsed, svd_events,
        callback_events[], discarded_weight[], maximum_local[])
end

function decode_component(decoder, syndrome, algorithm)
    wall_started = time_ns()
    alg = algorithm_symbol(algorithm)
    n = decoder.detector_count
    length(syndrome) == n || throw(DimensionMismatch("incorrect syndrome length"))
    all(bit -> bit == 0 || bit == 1, syndrome) ||
        throw(ArgumentError("syndrome must be binary"))
    iseven(count(!iszero, syndrome)) ||
        throw(ArgumentError("odd parity syndrome is impossible"))
    cfg = decoder.settings
    psi = initial_boundary(decoder)
    log_syndrome = 0.0
    frontier = 0
    total_elapsed = 0.0
    total_svd_events = 0
    total_callback_events = 0
    total_discarded = 0.0
    maximum_local = 0.0
    maximum_bond = 1

    for vertex in 1:n
        psi = TopoNoise._toric3d_z_tn_fuse_current(psi)
        groups = decoder.groups[vertex]
        support = isempty(groups) ? 1 : groups[end].neighbor - vertex + 1
        frontier = max(frontier,
            isempty(groups) ? vertex : groups[end].neighbor)
        active_sites = frontier - vertex + 1
        star = TopoNoise._toric3d_z_tn_star_mpo(
            psi, groups, vertex, syndrome[vertex])
        psi, stats = apply_update(psi, star, alg;
            maxdim=cfg.max_bond_dimension, cutoff=cfg.cutoff)
        total_elapsed += stats.elapsed_seconds
        total_svd_events += stats.svd_events
        total_callback_events += stats.callback_events
        total_discarded += stats.discarded_weight
        maximum_local = max(maximum_local,
            stats.maximum_local_truncation_error)
        maximum_bond = max(maximum_bond,
            TopoNoise._toric3d_z_tn_active_maxlinkdim(psi, active_sites))
        mass = TopoNoise._toric3d_z_tn_mass(psi, active_sites, 1.0)
        isfinite(mass) && mass > 0 ||
            throw(DomainError(mass, "invalid conditioned mass at vertex $vertex"))
        psi[1] /= mass
        log_syndrome += log(mass)
    end

    length(psi) == 1 || error("component contraction did not eliminate all detectors")
    logical = ITensorMPS.siteind(only, psi, 1)
    weights = vec(Array(psi[1], logical))
    total = sum(weights)
    isfinite(total) && total > 0 ||
        throw(DomainError(total, "invalid final logical mass"))
    probabilities = weights ./ total
    log_syndrome += log(total)
    valid = all(isfinite, probabilities) &&
        minimum(probabilities) >= -cfg.negative_tolerance
    wall_seconds = (time_ns() - wall_started) / 1e9
    return (; valid, probabilities,
        logical_mask=valid ? argmax(probabilities) - 1 : -1,
        raw_argmax=argmax(probabilities) - 1,
        log_syndrome, elapsed_seconds=wall_seconds,
        update_seconds=total_elapsed,
        svd_events=total_svd_events,
        callback_events=total_callback_events,
        discarded_weight=total_discarded,
        maximum_local_truncation_error=maximum_local,
        actual_maxbond=maximum_bond, peak_rss_bytes=Int(Sys.maxrss()))
end

function result_document(result)
    return Dict{String,Any}(
        "status" => "complete", "valid" => result.valid,
        "probabilities" => collect(result.probabilities),
        "logical_mask" => result.logical_mask,
        "raw_argmax" => result.raw_argmax,
        "log_syndrome" => result.log_syndrome,
        "elapsed_seconds" => result.elapsed_seconds,
        "update_seconds" => result.update_seconds,
        "svd_events" => result.svd_events,
        "callback_events" => result.callback_events,
        "discarded_weight" => result.discarded_weight,
        "maximum_local_truncation_error" =>
            result.maximum_local_truncation_error,
        "actual_maxbond" => result.actual_maxbond,
        "peak_rss_bytes" => result.peak_rss_bytes)
end

function combine_component_results(seed, transfer, recovery_mask::Integer,
        negative_tolerance::Real)
    probabilities = Model.combine_posteriors(
        seed.probabilities, transfer.probabilities, recovery_mask)
    total = sum(probabilities)
    isfinite(total) && total > 0 ||
        throw(DomainError(total, "invalid combined logical mass"))
    probabilities ./= total
    valid = seed.valid && transfer.valid && all(isfinite, probabilities) &&
        minimum(probabilities) >= -negative_tolerance
    return (; valid, probabilities,
        logical_mask=valid ? argmax(probabilities) - 1 : -1,
        raw_argmax=argmax(probabilities) - 1,
        log_joint_syndrome=seed.log_syndrome + transfer.log_syndrome,
        elapsed_seconds=seed.elapsed_seconds + transfer.elapsed_seconds,
        update_seconds=seed.update_seconds + transfer.update_seconds,
        svd_events=seed.svd_events + transfer.svd_events,
        callback_events=seed.callback_events + transfer.callback_events,
        discarded_weight=seed.discarded_weight + transfer.discarded_weight,
        maximum_local_truncation_error=max(
            seed.maximum_local_truncation_error,
            transfer.maximum_local_truncation_error),
        actual_maxbond=max(seed.actual_maxbond, transfer.actual_maxbond),
        peak_rss_bytes=max(seed.peak_rss_bytes, transfer.peak_rss_bytes),
        seed_result=seed, transfer_result=transfer)
end

function compare_small_decoder(; L::Integer=2, p::Real=0.02,
        maxdim::Integer=512, cutoff::Real=0.0, syndrome=nothing)
    context = Model.build_context(L)
    settings = Toric3DZTNSettings(max_bond_dimension=maxdim, cutoff=cutoff)
    decoder = Model.build_decoder(context, p; settings).transfer
    observed = isnothing(syndrome) ? falses(decoder.detector_count) : syndrome
    report = Dict{String,Any}(
        "L" => Int(L), "p" => Float64(p), "maxdim" => Int(maxdim),
        "cutoff" => Float64(cutoff),
        "densitymatrix" => Dict{String,Any}(
            "status" => "unsupported_finite_prefix",
            "reason" => "ITensorMPS densitymatrix replaces the live suffix bond"))
    for algorithm in (:current, :zipup)
        report[String(algorithm)] = result_document(
            decode_component(decoder, observed, algorithm))
    end
    current = report["current"]
    zipup = report["zipup"]
    report["comparison"] = Dict{String,Any}(
        "posterior_linf" => maximum(abs.(
            current["probabilities"] .- zipup["probabilities"])),
        "abs_log_syndrome_delta" => abs(
            current["log_syndrome"] - zipup["log_syndrome"]),
        "same_logical_mask" =>
            current["raw_argmax"] == zipup["raw_argmax"])
    return report
end

const REPLAY_ALGORITHMS = (:current, :zipup)
const REPLAY_BLAS_THREADS = 4

input_probe_module() = InputProbe

function replay_algorithm(value)
    algorithm = algorithm_symbol(value)
    algorithm in REPLAY_ALGORITHMS ||
        throw(ArgumentError("replay algorithm must be current or zipup"))
    return algorithm
end

replay_case_directory(output_dir, algorithm) = joinpath(
    input_probe_module().canonical_destination(output_dir), String(algorithm))
replay_result_path(output_dir, algorithm) =
    joinpath(replay_case_directory(output_dir, algorithm), "result.toml")
replay_failure_path(output_dir, algorithm) =
    joinpath(replay_case_directory(output_dir, algorithm), "failure.toml")

function replay_manifest(inputs)
    probe = input_probe_module()
    script_path = @__FILE__
    wrapper_path = joinpath(@__DIR__, "update_algorithm_benchmark.sbatch")
    isfile(wrapper_path) || error("update benchmark Slurm wrapper is missing")
    document = Dict{String,Any}(
        "schema_version" => 1,
        "purpose" => "Read-only finite-prefix update algorithm benchmark; excluded from production statistics",
        "source_run_directory" => inputs.run_dir,
        "source_run_id" => inputs.metadata["run_id"],
        "source_metadata_sha256" => inputs.metadata_sha256,
        "source_corpus_sha256" => inputs.corpus_sha256,
        "record_id" => inputs.record["record_id"],
        "shot_index" => 1,
        "L" => inputs.metadata["L"], "p" => inputs.metadata["p"],
        "chi" => inputs.metadata["chi"],
        "cutoff" => inputs.metadata["cutoff"],
        "negative_tolerance" => inputs.metadata["negative_tolerance"],
        "algorithms" => String.(collect(REPLAY_ALGORITHMS)),
        "blas_threads" => REPLAY_BLAS_THREADS,
        "julia_threads" => 1,
        "probe_script_sha256" => probe.filehash(script_path),
        "sbatch_sha256" => probe.filehash(wrapper_path),
        "densitymatrix_status" => "unsupported_finite_prefix",
        "densitymatrix_reason" =>
            "ITensorMPS densitymatrix replaces the live suffix bond",
        "zipup_schedule" =>
            "Initial zipup sweep uses cutoff/10 and 2*chi; final truncate sweep uses cutoff and chi.",
        "truth_policy" =>
            "Only observed seed/final syndromes and known seed recovery enter the benchmark; raw faults and truth labels are excluded.")
    document["probe_id"] = probe.Runner.identity_hash(document)
    return document
end

function initialize_replay(inputs, output_directory)
    probe = input_probe_module()
    output = probe.canonical_destination(output_directory)
    (probe.path_within(output, inputs.run_dir) ||
        probe.path_within(inputs.run_dir, output)) &&
        error("benchmark output and source run must be disjoint")
    manifest = replay_manifest(inputs)
    path = joinpath(output, "probe-metadata.toml")
    if ispath(output) && !isdir(output)
        error("benchmark output path is not a directory")
    elseif isdir(output) && !isfile(path) && !isempty(readdir(output))
        error("nonempty benchmark directory has no compatible manifest")
    elseif isfile(path)
        TOML.parsefile(path) == manifest ||
            error("benchmark manifest differs; use a fresh output directory")
    else
        probe.Runner.atomic_toml(path, manifest)
    end
    return manifest
end

function component_document_for_replay(component)
    return Dict{String,Any}(
        "valid" => component.valid,
        "probabilities" => collect(component.probabilities),
        "logical_mask" => component.logical_mask,
        "raw_argmax" => component.raw_argmax,
        "log_syndrome" => component.log_syndrome,
        "elapsed_seconds" => component.elapsed_seconds,
        "update_seconds" => component.update_seconds,
        "svd_events" => component.svd_events,
        "callback_events" => component.callback_events,
        "discarded_weight" => component.discarded_weight,
        "maximum_local_truncation_error" =>
            component.maximum_local_truncation_error,
        "actual_maxbond" => component.actual_maxbond,
        "peak_rss_bytes" => component.peak_rss_bytes)
end

function replay_progress(case_dir, stage, algorithm)
    probe = input_probe_module()
    now = time()
    probe.Runner.atomic_toml(joinpath(case_dir, "progress.toml"),
        Dict{String,Any}("pid" => getpid(), "stage" => stage,
            "shot" => 1, "algorithm" => String(algorithm),
            "stage_started_unix" => now, "updated_at_unix" => now))
    println("stage=$stage shot=1/1 algorithm=$algorithm")
    flush(stdout)
end

function run_replay_worker(inputs, output_directory, requested_algorithm)
    algorithm = replay_algorithm(requested_algorithm)
    Threads.nthreads() == 1 || error("run benchmark with one Julia thread")
    manifest = initialize_replay(inputs, output_directory)
    case_dir = replay_case_directory(output_directory, algorithm)
    ispath(case_dir) && error("algorithm output already exists: $case_dir")
    mkpath(case_dir)
    BLAS.set_num_threads(REPLAY_BLAS_THREADS)
    BLAS.get_num_threads() == REPLAY_BLAS_THREADS ||
        error("failed to apply four BLAS threads")

    replay_progress(case_dir, "setup", algorithm)
    context = Model.build_context(inputs.metadata["L"])
    inputs.record["raw_fault_count"] == size(context.model.matrix, 2) ||
        error("saved raw fault count differs from reconstructed model")
    settings = Toric3DZTNSettings(
        max_bond_dimension=inputs.metadata["chi"],
        cutoff=inputs.metadata["cutoff"],
        negative_tolerance=inputs.metadata["negative_tolerance"])
    decoder = Model.build_decoder(context, inputs.metadata["p"]; settings)

    replay_progress(case_dir, "seed", algorithm)
    seed = decode_component(
        decoder.seed, inputs.observed.seed_syndrome, algorithm)
    replay_progress(case_dir, "transfer", algorithm)
    transfer = decode_component(
        decoder.transfer, inputs.observed.final_syndrome, algorithm)
    combined = combine_component_results(seed, transfer,
        inputs.observed.recovery_mask, inputs.metadata["negative_tolerance"])
    document = Dict{String,Any}(
        "schema_version" => 1, "probe_id" => manifest["probe_id"],
        "source_run_id" => inputs.metadata["run_id"],
        "record_id" => inputs.record["record_id"], "shot_index" => 1,
        "algorithm" => String(algorithm),
        "L" => inputs.metadata["L"], "p" => inputs.metadata["p"],
        "chi" => inputs.metadata["chi"],
        "cutoff" => inputs.metadata["cutoff"],
        "blas_threads" => BLAS.get_num_threads(),
        "julia_threads" => Threads.nthreads(),
        "status" => combined.valid ? "complete" : "computed_invalid",
        "valid" => combined.valid,
        "probabilities" => collect(combined.probabilities),
        "logical_mask" => combined.logical_mask,
        "raw_argmax" => combined.raw_argmax,
        "log_joint_syndrome" => combined.log_joint_syndrome,
        "elapsed_seconds" => combined.elapsed_seconds,
        "update_seconds" => combined.update_seconds,
        "svd_events" => combined.svd_events,
        "callback_events" => combined.callback_events,
        "discarded_weight" => combined.discarded_weight,
        "maximum_local_truncation_error" =>
            combined.maximum_local_truncation_error,
        "actual_maxbond" => combined.actual_maxbond,
        "peak_rss_bytes" => combined.peak_rss_bytes,
        "seed_result" => component_document_for_replay(seed),
        "transfer_result" => component_document_for_replay(transfer))
    input_probe_module().Runner.atomic_toml(
        replay_result_path(output_directory, algorithm), document)
    replay_progress(case_dir, document["status"], algorithm)
    return document
end

function checkpoint_replay_failure(inputs, output_directory,
        requested_algorithm, status, exit_code, limit_seconds)
    algorithm = replay_algorithm(requested_algorithm)
    status in ("contraction_time_limit", "worker_failed") ||
        error("invalid failure status")
    manifest = initialize_replay(inputs, output_directory)
    case_dir = replay_case_directory(output_directory, algorithm)
    mkpath(case_dir)
    document = Dict{String,Any}(
        "schema_version" => 1, "probe_id" => manifest["probe_id"],
        "source_run_id" => inputs.metadata["run_id"],
        "record_id" => inputs.record["record_id"], "shot_index" => 1,
        "algorithm" => String(algorithm), "status" => String(status),
        "valid" => false, "exit_code" => Int(exit_code),
        "limit_seconds" => Float64(limit_seconds))
    input_probe_module().Runner.atomic_toml(
        replay_failure_path(output_directory, algorithm), document)
    return document
end

function collect_replay(inputs, output_directory)
    manifest = initialize_replay(inputs, output_directory)
    report = Dict{String,Any}(
        "schema_version" => 1, "probe_id" => manifest["probe_id"],
        "source_run_id" => inputs.metadata["run_id"],
        "record_id" => inputs.record["record_id"],
        "densitymatrix" => Dict{String,Any}(
            "status" => "unsupported_finite_prefix",
            "reason" => "ITensorMPS densitymatrix replaces the live suffix bond"))
    cases = Dict{String,Any}()
    for algorithm in REPLAY_ALGORITHMS
        result_path = replay_result_path(output_directory, algorithm)
        failure_path = replay_failure_path(output_directory, algorithm)
        cases[String(algorithm)] = if isfile(result_path)
            TOML.parsefile(result_path)
        elseif isfile(failure_path)
            TOML.parsefile(failure_path)
        else
            Dict{String,Any}("algorithm" => String(algorithm),
                "status" => "missing", "valid" => false)
        end
    end
    report["cases"] = cases
    if all(get(cases[String(algorithm)], "status", nothing) in
            ("complete", "computed_invalid") for algorithm in REPLAY_ALGORITHMS)
        current, zipup = cases["current"], cases["zipup"]
        comparison = Dict{String,Any}(
            "posterior_linf" => maximum(abs.(
                current["probabilities"] .- zipup["probabilities"])),
            "abs_log_joint_syndrome_delta" => abs(
                current["log_joint_syndrome"] -
                zipup["log_joint_syndrome"]),
            "same_raw_argmax" =>
                current["raw_argmax"] == zipup["raw_argmax"],
            "current_over_zipup_speedup" =>
                current["elapsed_seconds"] / zipup["elapsed_seconds"],
            "current_over_zipup_update_speedup" =>
                current["update_seconds"] / zipup["update_seconds"],
            "zipup_over_current_svd_ratio" =>
                zipup["svd_events"] / current["svd_events"])
        report["comparison"] = comparison
        report["status"] = "complete"
    else
        report["status"] = "incomplete"
    end
    path = joinpath(input_probe_module().canonical_destination(output_directory),
        "report.toml")
    input_probe_module().Runner.atomic_toml(path, report)
    TOML.print(stdout, report; sorted=true)
    println()
    return report
end

function main(args=ARGS)
    length(args) >= 3 || error(
        "usage: --validate|--worker|--failure|--collect SOURCE_RUN OUTPUT [arguments]")
    mode, source_run, output = args[1:3]
    probe = input_probe_module()
    inputs = probe.read_inputs(source_run)
    if mode == "--validate"
        length(args) == 3 || error("--validate takes SOURCE_RUN OUTPUT")
        initialize_replay(inputs, output)
    elseif mode == "--worker"
        length(args) == 4 || error("--worker requires ALGORITHM")
        run_replay_worker(inputs, output, Symbol(args[4]))
    elseif mode == "--failure"
        length(args) == 7 ||
            error("--failure requires ALGORITHM STATUS EXIT LIMIT")
        checkpoint_replay_failure(inputs, output, Symbol(args[4]), args[5],
            parse(Int, args[6]), parse(Float64, args[7]))
    elseif mode == "--collect"
        length(args) == 3 || error("--collect takes SOURCE_RUN OUTPUT")
        report = collect_replay(inputs, output)
        get(report, "status", nothing) == "complete" || exit(2)
    else
        error("unknown mode: $mode")
    end
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    TNMLUpdateAlgorithmBenchmark.main()
end
