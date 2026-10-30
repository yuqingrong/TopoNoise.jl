#!/usr/bin/env julia
"""
One-record adapter for the adaptive pilot.

This file contains orchestration and serialization only. TN decoding calls the
existing `RepairedSeedTN.decode_record`; exact decoding calls the existing
`TNMLRepairedRunner.exact_reference` in an isolated scratch directory.
"""

using TOML, LinearAlgebra

const HERE = @__DIR__
include(joinpath(HERE, "..", "run.jl"))
const Runner = TNMLRepairedRunner

function repaired_seed_function(name::Symbol)
    backend = Base.invokelatest(getproperty, Runner, :RepairedSeedTN)
    return Base.invokelatest(getproperty, backend, name)
end

function tn_settings(chi, cutoff)
    # load_backend! creates these bindings while the worker is already running.
    # Both the lookup and the keyword constructor must see the latest world.
    backend = Base.invokelatest(getproperty, Runner, :TopoNoise)
    constructor = Base.invokelatest(getproperty, backend, :Toric3DZTNSettings)
    return Base.invokelatest(constructor;
        max_bond_dimension=chi,
        cutoff=cutoff,
        negative_tolerance=Runner.NEGATIVE_TOLERANCE,
    )
end

function atomic_toml(path, document)
    mkpath(dirname(path))
    temporary, io = mktemp(dirname(path))
    try
        TOML.print(io, Runner.toml_value(document); sorted=true)
        flush(io)
        close(io)
        mv(temporary, path; force=true)
    finally
        isopen(io) && close(io)
        isfile(temporary) && rm(temporary)
    end
end

function exception_result(error, backtrace, elapsed)
    status = error isa OutOfMemoryError ? "memory_limit" : "exception"
    return Dict{String,Any}(
        "status" => status,
        "valid" => false,
        "elapsed_seconds" => elapsed,
        "peak_rss_bytes" => Int(Sys.maxrss()),
        "rss_scope" => "process_lifetime_high_water",
        "error" => sprint(showerror, error, backtrace),
    )
end

function authenticate_record(path)
    document = TOML.parsefile(path)
    record_id = get(document, "record_id", nothing)
    record_id isa String || throw(ArgumentError("record_id is missing"))
    payload = Dict(key => value for (key, value) in document if key != "record_id")
    Runner.identity_hash(payload) == record_id ||
        throw(ArgumentError("record failed its established sorted-TOML identity"))
    return document
end

function tn_result(decoder, record)
    observed = Runner.observed_record(record)
    decoded = Base.invokelatest(
        repaired_seed_function(:decode_record), decoder, observed
    )
    output = Dict{String,Any}(
        "status" => decoded.valid ? "complete" : "invalid",
        "valid" => decoded.valid,
        "elapsed_seconds" => 0.0,
        "peak_rss_bytes" => Int(Sys.maxrss()),
        "rss_scope" => "process_lifetime_high_water",
    )
    if decoded.valid
        output["probabilities"] = collect(decoded.probabilities)
        output["logical_mask"] = Int(decoded.logical_mask)
        output["log_joint_syndrome"] = Float64(decoded.log_joint_syndrome)
    end
    return output
end

function tn_result(context, record, p, chi, cutoff)
    settings = tn_settings(chi, cutoff)
    decoder = Base.invokelatest(
        repaired_seed_function(:build_decoder), context, p; settings=settings
    )
    return tn_result(decoder, record)
end

function runtime_identity()
    return Dict{String,Any}(
        "julia_version" => string(VERSION),
        "blas_vendor" => string(BLAS.vendor()),
        "cpu_model" => string(Sys.CPU_NAME),
        "julia_threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(),
    )
end

function configured_blas_threads()
    value = parse(Int, get(ENV, "TNML_BLAS_THREADS", "4"))
    value > 0 || throw(ArgumentError("TNML_BLAS_THREADS must be positive"))
    return value
end

function exact_result(context, record, record_document, L, p)
    return mktempdir() do scratch
        cfg = (L=L, p=p, shots=1, chi=1, seed=0, output_dir=scratch,
               exact=true, cutoff=0.0)
        run_id = Runner.identity_hash(Dict(
            "adapter" => "adaptive_chi_exact_v1",
            "record_id" => record_document["record_id"],
            "L" => L,
            "p" => p,
        ))
        Base.invokelatest(
            Runner.exact_reference,
            cfg,
            context,
            record,
            run_id,
            record_document,
        )
    end
end

function main(args=ARGS)
    if length(args) == 2 && args[1] == "probe"
        BLAS.set_num_threads(configured_blas_threads())
        atomic_toml(abspath(args[2]), Dict("runtime_identity" => runtime_identity()))
        return nothing
    end
    length(args) == 7 || throw(ArgumentError(
        "expected probe result.toml or mode record.toml result.toml L p chi cutoff"))
    mode, record_path, output_path = args[1], abspath(args[2]), abspath(args[3])
    mode in ("tn", "exact") || throw(ArgumentError("mode must be tn or exact"))
    L, p, chi, cutoff = parse(Int, args[4]), parse(Float64, args[5]),
                         parse(Int, args[6]), parse(Float64, args[7])
    Threads.nthreads() == 1 || throw(ArgumentError("adapter requires one Julia thread"))
    BLAS.set_num_threads(configured_blas_threads())
    started = time_ns()
    output = try
        document = authenticate_record(record_path)
        Runner.load_backend!()
        context = Base.invokelatest(repaired_seed_function(:build_context), L)
        record = Runner.restore_record(document)
        result = mode == "tn" ?
            tn_result(context, record, p, chi, cutoff) :
            exact_result(context, record, document, L, p)
        result["elapsed_seconds"] = (time_ns() - started) / 1e9
        result["peak_rss_bytes"] = Int(Sys.maxrss())
        result["rss_scope"] = "process_lifetime_high_water"
        result
    catch error
        error isa InterruptException && rethrow()
        exception_result(error, catch_backtrace(), (time_ns() - started) / 1e9)
    end
    output["runtime_identity"] = runtime_identity()
    atomic_toml(output_path, output)
    return nothing
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
