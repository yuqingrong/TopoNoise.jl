#!/usr/bin/env julia
"""Persistent, serial Julia worker for one frozen adaptive ``(L,p)`` task."""

include(joinpath(@__DIR__, "julia_adapter.jl"))

const PERSISTENT_PROTOCOL_VERSION = 1

function option(args, name)
    index = findfirst(==(name), args)
    index === nothing && throw(ArgumentError("missing worker option $name"))
    index == length(args) && throw(ArgumentError("worker option $name has no value"))
    return args[index + 1]
end

function request_envelope(request)
    names = ("protocol_version", "session_token", "sequence", "request_id",
             "mode", "source_bundle_id", "runtime_profile_id", "L", "p",
             "p_key", "chi", "cutoff", "record_id", "julia_record_id",
             "executor_mode", "wall_time_scope", "decoder_cache_scope")
    return Dict(name => request[name] for name in names)
end

function ready_document(config; status="ready", valid=true, error=nothing,
                        startup_elapsed_seconds=0.0)
    result = Dict{String,Any}(
        "status" => status,
        "valid" => valid,
        "protocol_version" => PERSISTENT_PROTOCOL_VERSION,
        "session_token" => config.session_token,
        "source_bundle_id" => config.source_bundle_id,
        "runtime_profile_id" => config.runtime_profile_id,
        "L" => config.L,
        "p" => config.p,
        "p_key" => config.p_key,
        "cutoff" => config.cutoff,
        "startup_elapsed_seconds" => startup_elapsed_seconds,
        "peak_rss_bytes" => Int(Sys.maxrss()),
        "rss_scope" => "worker_lifetime_high_water",
        "executor_mode" => "persistent",
        "wall_time_scope" => "request_after_worker_ready",
        "decoder_cache_scope" => "worker_lifetime_by_chi",
        "runtime_identity" => runtime_identity(),
    )
    error === nothing || (result["error"] = error)
    return result
end

function worker_config(args)
    mode = "--serve" in args ? "serve" : "--probe" in args ? "probe" :
        throw(ArgumentError("worker requires --serve or --probe"))
    L = parse(Int, option(args, "--L"))
    p = parse(Float64, option(args, "--p"))
    p_key = option(args, "--p-key")
    cutoff = parse(Float64, option(args, "--cutoff"))
    return (
        mode=mode,
        L=L,
        p=p,
        p_key=p_key,
        cutoff=cutoff,
        ready=abspath(option(args, "--ready")),
        session_token=option(args, "--session-token"),
        source_bundle_id=option(args, "--source-bundle-id"),
        runtime_profile_id=option(args, "--runtime-profile-id"),
    )
end

function validate_request(request, config, expected_sequence, request_path)
    get(request, "protocol_version", nothing) == PERSISTENT_PROTOCOL_VERSION ||
        throw(ArgumentError("request protocol version mismatch"))
    get(request, "session_token", nothing) == config.session_token ||
        throw(ArgumentError("request session token mismatch"))
    get(request, "sequence", nothing) == expected_sequence ||
        throw(ArgumentError("request sequence mismatch"))
    get(request, "source_bundle_id", nothing) == config.source_bundle_id ||
        throw(ArgumentError("request source identity mismatch"))
    get(request, "runtime_profile_id", nothing) == config.runtime_profile_id ||
        throw(ArgumentError("request runtime identity mismatch"))
    get(request, "L", nothing) == config.L ||
        throw(ArgumentError("request L mismatch"))
    get(request, "p", nothing) == config.p ||
        throw(ArgumentError("request p mismatch"))
    get(request, "p_key", nothing) == config.p_key ||
        throw(ArgumentError("request probability stratum mismatch"))
    get(request, "executor_mode", nothing) == "persistent" ||
        throw(ArgumentError("request executor mode mismatch"))
    get(request, "wall_time_scope", nothing) == "request_after_worker_ready" ||
        throw(ArgumentError("request wall-time scope mismatch"))
    get(request, "decoder_cache_scope", nothing) == "worker_lifetime_by_chi" ||
        throw(ArgumentError("request decoder-cache scope mismatch"))
    mode = get(request, "mode", nothing)
    mode in ("tn", "exact") || throw(ArgumentError("request mode is invalid"))
    cutoff = get(request, "cutoff", nothing)
    if mode == "tn"
        cutoff == config.cutoff || throw(ArgumentError("request cutoff mismatch"))
    else
        cutoff == 0.0 || throw(ArgumentError("exact request cutoff must be zero"))
    end
    chi = get(request, "chi", nothing)
    chi isa Integer && chi > 0 || throw(ArgumentError("request chi is invalid"))
    session_dir = dirname(config.ready)
    dirname(abspath(request_path)) == session_dir ||
        throw(ArgumentError("request path is outside the private worker session"))
    for name in ("record_path", "output_path")
        dirname(abspath(request[name])) == session_dir ||
            throw(ArgumentError("$name is outside the private worker session"))
    end
    isfile(request["output_path"]) && throw(ArgumentError("request output already exists"))
    return mode
end

function cached_decoder!(cache, context, config, chi)
    haskey(cache, chi) && return cache[chi], true, 0.0
    started = time_ns()
    settings = tn_settings(chi, config.cutoff)
    decoder = Base.invokelatest(
        repaired_seed_function(:build_decoder), context, config.p; settings=settings
    )
    cache[chi] = decoder
    return decoder, false, (time_ns() - started) / 1e9
end

function handle_request(request, config, context, cache)
    started = time_ns()
    authenticated = try
        document = authenticate_record(request["record_path"])
        document["record_id"] == request["julia_record_id"] ||
            throw(ArgumentError("request Julia record identity mismatch"))
        (document=document, record=Runner.restore_record(document))
    catch error
        error isa InterruptException && rethrow()
        result = exception_result(
            error, catch_backtrace(), (time_ns() - started) / 1e9
        )
        result["status"] = "fatal_identity"
        result["peak_rss_bytes"] = Int(Sys.maxrss())
        result["rss_scope"] = "worker_lifetime_high_water"
        result["runtime_identity"] = runtime_identity()
        merge!(result, request_envelope(request))
        atomic_toml(request["output_path"], result)
        return false
    end
    result = try
        if request["mode"] == "tn"
            decoder, cache_hit, build_seconds =
                cached_decoder!(cache, context, config, request["chi"])
            decode_started = time_ns()
            output = tn_result(decoder, authenticated.record)
            output["decoder_cache_hit"] = cache_hit
            output["decoder_build_seconds"] = build_seconds
            output["decode_seconds"] = (time_ns() - decode_started) / 1e9
            output
        else
            output = exact_result(
                context,
                authenticated.record,
                authenticated.document,
                config.L,
                config.p,
            )
            output["decoder_cache_hit"] = false
            output["decoder_build_seconds"] = 0.0
            output["decode_seconds"] = output["elapsed_seconds"]
            output
        end
    catch error
        error isa InterruptException && rethrow()
        exception_result(error, catch_backtrace(), (time_ns() - started) / 1e9)
    end
    result["elapsed_seconds"] = (time_ns() - started) / 1e9
    result["peak_rss_bytes"] = Int(Sys.maxrss())
    result["rss_scope"] = "worker_lifetime_high_water"
    result["runtime_identity"] = runtime_identity()
    merge!(result, request_envelope(request))
    atomic_toml(request["output_path"], result)
    return true
end

function main(args=ARGS)
    config = worker_config(args)
    Threads.nthreads() == 1 || throw(ArgumentError("worker requires one Julia thread"))
    BLAS.set_num_threads(configured_blas_threads())
    started = time_ns()
    context = try
        Runner.load_backend!()
        Base.invokelatest(repaired_seed_function(:build_context), config.L)
    catch error
        output = ready_document(
            config;
            status="exception",
            valid=false,
            error=sprint(showerror, error, catch_backtrace()),
            startup_elapsed_seconds=(time_ns() - started) / 1e9,
        )
        atomic_toml(config.ready, output)
        return 2
    end
    atomic_toml(
        config.ready,
        ready_document(
            config; startup_elapsed_seconds=(time_ns() - started) / 1e9
        ),
    )
    config.mode == "probe" && return 0

    cache = Dict{Int,Any}()
    expected_sequence = 1
    for raw_line in eachline(stdin)
        line = strip(raw_line)
        isempty(line) && continue
        if startswith(line, "CLOSE ")
            split(line; limit=2)[2] == config.session_token ||
                throw(ArgumentError("worker close token mismatch"))
            return 0
        end
        request_path = abspath(line)
        request = TOML.parsefile(request_path)
        validate_request(request, config, expected_sequence, request_path)
        keep_running = handle_request(request, config, context, cache)
        expected_sequence += 1
        keep_running || return 3
    end
    return 0
end

abspath(PROGRAM_FILE) == (@__FILE__) && exit(main())
