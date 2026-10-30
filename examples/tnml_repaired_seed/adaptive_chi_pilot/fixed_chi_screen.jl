#!/usr/bin/env julia
"""Paired fixed-chi screening on fresh L=3 histories, not a production study.

Every scored history is decoded at chi=64 and chi=128, irrespective of the
result of the other decoder. The first history is marked cold for timing;
there are no unscored warm-up histories. Exceptions preserve evidence and fail
the job; numerical invalidity and exact-reference disagreement do not stop it.
"""

include(joinpath(@__DIR__, "persistent_worker.jl"))
using Random, SHA

const SCREEN_P_KEYS = ("0.0060", "0.0100", "0.0200", "0.0300", "0.0350", "0.0425")
const SCREEN_CHIS = (64, 128)
const SCREEN_SEED_BASE = 202610700000
const SCREEN_REPOSITORY = normpath(joinpath(@__DIR__, "..", "..", ".."))

function screen_source_hashes()
    names = String[]
    for (directory, _, files) in walkdir(joinpath(SCREEN_REPOSITORY, "src"))
        for name in files
            endswith(name, ".jl") && push!(names, relpath(joinpath(directory, name), SCREEN_REPOSITORY))
        end
    end
    append!(names, [
        "Project.toml",
        "examples/tnml_repaired_seed/model.jl",
        "examples/tnml_repaired_seed/run.jl",
        "examples/tnml_repaired_seed/dense_frontier.jl",
        ["examples/tnml_repaired_seed/adaptive_chi_pilot/" * name for name in (
            "persistent_worker.jl", "julia_adapter.jl", "fixed_chi_screen.jl",
            "inspect_fixed_chi_screen.py", "core.py", "adaptive_runner.py",
            "inspect_adaptive_ladder.py",
        )]...,
    ])
    isfile(joinpath(SCREEN_REPOSITORY, "Manifest.toml")) && push!(names, "Manifest.toml")
    return Dict(name => bytes2hex(sha256(read(joinpath(SCREEN_REPOSITORY, name))))
                for name in sort!(unique!(names)))
end

function save_screen(path, document)
    ispath(path) && throw(ArgumentError("refusing to overwrite screening output: $path"))
    atomic_toml(path, document)
end

function screen_exception(path, error, backtrace, started; fields...)
    document = exception_result(error, backtrace, (time_ns() - started) / 1e9)
    for (key, value) in fields
        document[string(key)] = value
    end
    save_screen(path, document)
end

function screen_history(context, cache, config, output, shot_index, index)
    cold_start = shot_index == 1
    seed = SCREEN_SEED_BASE + 100 * index + shot_index - 1
    directory = joinpath(output, "shots", "shot_" * lpad(string(shot_index), 4, '0'))
    mkdir(directory)
    println("stage=sample L=3 p=", config.p_key, " shot=", shot_index,
            "/", config.shots, " seed=", seed, " cold_start=", cold_start)
    flush(stdout)
    sampled = Base.invokelatest(repaired_seed_function(:sample_record), context, config.p, Xoshiro(seed))
    record_path = joinpath(directory, "record.toml")
    save_screen(record_path, Runner.record_document(sampled, shot_index))
    document = authenticate_record(record_path)
    record = Runner.restore_record(document)
    record_id = document["record_id"]
    order = isodd(shot_index) ? (64, 128) : (128, 64)
    save_screen(joinpath(directory, "metadata.toml"), Dict(
        "shot_index" => shot_index, "seed" => seed, "record_id" => record_id,
        "record_file_sha256" => bytes2hex(sha256(read(record_path))),
        "record_authenticated" => true, "cold_start" => cold_start,
        "chi_order" => collect(order),
    ))

    println("stage=exact L=3 p=", config.p_key, " shot=", shot_index, "/", config.shots)
    flush(stdout)
    exact_started = time_ns()
    exact = try
        exact_result(context, record, document, 3, config.p)
    catch error
        screen_exception(joinpath(directory, "exact.toml"), error, catch_backtrace(), exact_started;
                         record_id=record_id, shot_index=shot_index)
        rethrow()
    end
    exact["record_id"] = record_id
    exact["shot_index"] = shot_index
    exact["elapsed_seconds"] = (time_ns() - exact_started) / 1e9
    exact["cold_start"] = cold_start
    save_screen(joinpath(directory, "exact.toml"), exact)
    get(exact, "status", nothing) == "complete" && get(exact, "valid", false) === true ||
        error("invalid independent exact reference for shot $shot_index")
    println("EXACT shot=", shot_index, " logical_mask=", exact["logical_mask"],
            " elapsed_seconds=", exact["elapsed_seconds"])
    flush(stdout)

    for chi in order
        println("stage=tn L=3 p=", config.p_key, " shot=", shot_index,
                "/", config.shots, " chi=", chi, " cold_start=", cold_start)
        flush(stdout)
        started = time_ns()
        path = joinpath(directory, "chi_" * lpad(string(chi), 4, '0') * ".toml")
        result = try
            decoder, cache_hit, build_seconds = cached_decoder!(cache, context, config, chi)
            decode_started = time_ns()
            result = tn_result(decoder, record)
            result["decode_seconds"] = (time_ns() - decode_started) / 1e9
            result["decoder_build_seconds"] = build_seconds
            result["decoder_cache_hit"] = cache_hit
            result
        catch error
            screen_exception(path, error, catch_backtrace(), started;
                             record_id=record_id, shot_index=shot_index, chi=chi,
                             cold_start=cold_start)
            rethrow()
        end
        result["record_id"] = record_id
        result["shot_index"] = shot_index
        result["chi"] = chi
        result["cold_start"] = cold_start
        result["elapsed_seconds"] = (time_ns() - started) / 1e9
        save_screen(path, result)
        println("RESULT shot=", shot_index, " chi=", chi, " status=", result["status"],
                " logical_mask=", get(result, "logical_mask", "unavailable"),
                " exact_mask=", exact["logical_mask"],
                " decode_seconds=", result["decode_seconds"],
                " build_seconds=", result["decoder_build_seconds"],
                " cache_hit=", result["decoder_cache_hit"], " cold_start=", cold_start)
        flush(stdout)
    end
end

function fixed_chi_screen(args=ARGS)
    length(args) in (2, 3) || throw(ArgumentError("usage: fixed_chi_screen.jl PINDEX OUTPUT_DIR [SHOTS=10]"))
    index = parse(Int, args[1])
    0 <= index < length(SCREEN_P_KEYS) || throw(ArgumentError("PINDEX must be 0 through 5"))
    shots = length(args) == 3 ? parse(Int, args[3]) : 10
    1 <= shots <= 10 || throw(ArgumentError("SHOTS must be 1 through 10"))
    output = abspath(args[2])
    isdir(output) || throw(ArgumentError("OUTPUT_DIR must be a fresh existing directory"))
    all(name == "worker.log" for name in readdir(output)) ||
        throw(ArgumentError("OUTPUT_DIR already contains scientific outputs"))
    Threads.nthreads() == 1 || error("run this screen with --threads=1")
    BLAS.set_num_threads(configured_blas_threads())
    config = (p=parse(Float64, SCREEN_P_KEYS[index + 1]), p_key=SCREEN_P_KEYS[index + 1],
              cutoff=1e-18, shots=shots)
    println("FIXED_SCREEN output=", output, " L=3 p=", config.p_key,
            " histories=", shots, " chis=64,128 production_gate=false")
    flush(stdout)
    started = time_ns()
    completed = 0
    try
        sources = screen_source_hashes()
        # Loading inside this function retains coverage of the world-age fix.
        Runner.load_backend!()
        context = Base.invokelatest(repaired_seed_function(:build_context), 3)
        cache = Dict{Int,Any}()
        save_screen(joinpath(output, "run.toml"), Dict(
            "schema_version" => 1, "kind" => "standalone_fixed_chi_screen",
            "L" => 3, "p" => config.p, "p_key" => config.p_key, "p_index" => index,
            "shots" => shots, "chis" => collect(SCREEN_CHIS), "seed_base" => SCREEN_SEED_BASE,
            "cutoff" => config.cutoff, "runtime_identity" => runtime_identity(),
            "startup_seconds" => (time_ns() - started) / 1e9,
            "startup_timing_scope" => "source_snapshot_backend_load_context_after_driver_import",
            "timing_policy" => "first_history_cold_excluded_from_warm_medians_all_histories_scored",
            "source_sha256" => sources, "production_gate" => false,
        ))
        mkdir(joinpath(output, "shots"))
        for shot_index in 1:shots
            screen_history(context, cache, config, output, shot_index, index)
            completed += 1
        end
        screen_source_hashes() == sources || error("source files changed during the fixed-chi screen")
        save_screen(joinpath(output, "completion.toml"), Dict(
            "status" => "complete", "requested_shots" => shots,
            "completed_shots" => completed, "elapsed_seconds" => (time_ns() - started) / 1e9,
            "source_identity_verified" => true, "production_gate" => false,
        ))
    catch error
        screen_exception(joinpath(output, "failure.toml"), error, catch_backtrace(), started;
                         requested_shots=shots, completed_shots=completed)
        rethrow()
    end
    println("FIXED_SCREEN complete L=3 p=", config.p_key, " histories=", completed,
            " paired_fixed_chis=64,128")
    flush(stdout)
    return nothing
end

abspath(PROGRAM_FILE) == (@__FILE__) && fixed_chi_screen()
