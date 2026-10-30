#!/usr/bin/env julia
"""One real L=3 history, all six chi levels, independent exact reference.

Standalone diagnostic only: this writes into a fresh empty output directory.
Invalid low-chi results are saved and escalated; exceptions abort the task.
All six contractions are intentional, so this is not an adaptive runtime test.
"""

include(joinpath(@__DIR__, "persistent_worker.jl"))
using Random, SHA

const DIAGNOSTIC_P_KEYS = ("0.0060", "0.0100", "0.0200", "0.0300", "0.0350", "0.0425")
const DIAGNOSTIC_LADDER = (32, 64, 128, 240, 480, 960)

function save_diagnostic(path, output)
    ispath(path) && throw(ArgumentError("refusing to overwrite diagnostic output: $path"))
    atomic_toml(path, output)
end

function ladder_diagnostic(args=ARGS)
    length(args) == 2 || throw(ArgumentError("usage: check_adaptive_ladder.jl PINDEX OUTPUT_DIR"))
    index = parse(Int, args[1])
    0 <= index < length(DIAGNOSTIC_P_KEYS) || throw(ArgumentError("PINDEX must be 0 through 5"))
    output = abspath(args[2])
    isdir(output) || throw(ArgumentError("OUTPUT_DIR must be an existing fresh directory"))
    # worker.log is opened by tee before Julia starts; no scientific outputs may exist.
    all(name == "worker.log" for name in readdir(output)) ||
        throw(ArgumentError("OUTPUT_DIR already contains scientific outputs"))
    Threads.nthreads() == 1 || error("run this diagnostic with --threads=1")
    BLAS.set_num_threads(configured_blas_threads())
    p_key = DIAGNOSTIC_P_KEYS[index + 1]
    p = parse(Float64, p_key)
    seed = 202610600000 + 100 * index
    println("DIAGNOSTIC output=", output, " L=3 p=", p_key,
            " seed=", seed, " histories=1 chi_levels=6 full_ladder=true")
    flush(stdout)
    started = time_ns()
    # Preserve runtime loading to exercise the original world-age regression.
    Runner.load_backend!()
    context = Base.invokelatest(repaired_seed_function(:build_context), 3)
    sampled = Base.invokelatest(repaired_seed_function(:sample_record), context, p, Xoshiro(seed))
    save_diagnostic(joinpath(output, "record.toml"), Runner.record_document(sampled, 1))
    document = authenticate_record(joinpath(output, "record.toml"))
    record = Runner.restore_record(document)
    record_id = document["record_id"]
    run = Dict{String,Any}(
        "kind" => "standalone_full_ladder_diagnostic", "L" => 3,
        "p" => p, "p_key" => p_key, "seed" => seed,
        "chi_ladder" => collect(DIAGNOSTIC_LADDER), "cutoff" => 1e-18,
        "record_id" => record_id, "record_authenticated" => true,
        "record_file_sha256" => bytes2hex(sha256(read(joinpath(output, "record.toml")))),
        "runtime_identity" => runtime_identity(),
        "startup_seconds" => (time_ns() - started) / 1e9,
        "startup_timing_scope" => "backend_load_context_and_sample_after_driver_import",
        "policy_source_sha256" => Dict(name => bytes2hex(sha256(read(joinpath(@__DIR__, name))))
            for name in ("core.py", "adaptive_runner.py", "inspect_adaptive_ladder.py")),
        "production_gate" => false,
    )
    save_diagnostic(joinpath(output, "run.toml"), run)
    println("stage=exact L=3 p=", p_key, " shot=1/1")
    flush(stdout)
    exact_started = time_ns()
    exact = exact_result(context, record, document, 3, p)
    exact["elapsed_seconds"] = (time_ns() - exact_started) / 1e9
    exact["record_id"] = record_id
    save_diagnostic(joinpath(output, "exact.toml"), exact)
    exact["valid"] === true || error("invalid independent exact reference")
    println("EXACT logical_mask=", exact["logical_mask"], " probabilities=", exact["probabilities"])
    flush(stdout)
    config = (p=p, cutoff=1e-18)
    cache = Dict{Int,Any}()
    for chi in DIAGNOSTIC_LADDER
        println("stage=tn L=3 p=", p_key, " chi=", chi, " shot=1/1")
        flush(stdout)
        level_started = time_ns()
        level_path = joinpath(output, "chi_" * lpad(string(chi), 4, '0') * ".toml")
        result = try
            decoder, cache_hit, build_seconds = cached_decoder!(cache, context, config, chi)
            decode_started = time_ns()
            result = tn_result(decoder, record)
            result["decode_seconds"] = (time_ns() - decode_started) / 1e9
            result["decoder_build_seconds"] = build_seconds
            result["decoder_cache_hit"] = cache_hit
            result
        catch error
            result = exception_result(error, catch_backtrace(), (time_ns() - level_started) / 1e9)
            result["chi"] = chi
            result["record_id"] = record_id
            save_diagnostic(level_path, result)
            rethrow()
        end
        result["chi"] = chi
        result["record_id"] = record_id
        result["elapsed_seconds"] = (time_ns() - level_started) / 1e9
        save_diagnostic(level_path, result)
        println("RESULT chi=", chi, " status=", result["status"],
                " logical_mask=", get(result, "logical_mask", "unavailable"),
                " decode_seconds=", result["decode_seconds"],
                " build_seconds=", result["decoder_build_seconds"])
        flush(stdout)
    end
    println("CHECK complete: full ladder saved; Python will evaluate production stopping policy")
    flush(stdout)
    return nothing
end

abspath(PROGRAM_FILE) == (@__FILE__) && ladder_diagnostic()
