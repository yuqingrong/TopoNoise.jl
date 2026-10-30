#!/usr/bin/env julia
"""One real L=3 history exercising late backend loading and decoder caching.

This is a diagnostic, not a production sample. It writes only stdout and uses
temporary files for its exact reference; existing frozen studies are untouched.
"""

include(joinpath(@__DIR__, "persistent_worker.jl"))
using Random, Test

function check_real_late_loading()
    Threads.nthreads() == 1 || error("run this diagnostic with --threads=1")
    BLAS.set_num_threads(configured_blas_threads())
    # Load inside this function to retain the older caller world that exposed
    # the original failure. Loading at top level would mask that regression.
    Runner.load_backend!()
    context = Base.invokelatest(repaired_seed_function(:build_context), 3)
    p = 0.006
    record = Base.invokelatest(
        repaired_seed_function(:sample_record), context, p, Xoshiro(202610600000)
    )
    document = Runner.record_document(record, 1)
    exact = exact_result(context, record, document, 3, p)
    println("EXACT logical_mask=", exact["logical_mask"],
            " probabilities=", exact["probabilities"])
    flush(stdout)
    config = (p=p, cutoff=1e-18)
    cache = Dict{Int,Any}()

    @testset "real TN after late backend loading" begin
        for chi in (32, 128, 960)
            decoder, hit, build_seconds = cached_decoder!(cache, context, config, chi)
            @test hit === false
            started = time_ns()
            result = Base.invokelatest(
                repaired_seed_function(:decode_record), decoder,
                Runner.observed_record(record),
            )
            decode_seconds = (time_ns() - started) / 1e9
            posterior_error = maximum(abs.(collect(result.probabilities) .- exact["probabilities"]))
            println("RESULT chi=", chi, " valid=", result.valid,
                    " logical_mask=", result.logical_mask,
                    " exact_mask=", exact["logical_mask"],
                    " posterior_linf=", posterior_error,
                    " build_seconds=", build_seconds,
                    " decode_seconds=", decode_seconds)
            flush(stdout)
            reused, hit, repeat_build = cached_decoder!(cache, context, config, chi)
            @test hit && reused === decoder && repeat_build == 0.0
            # A low-chi invalid probability is reported, never clamped or
            # accepted. Only the high-chi reference is required to pass here.
            if chi == 960
                @test result.valid
                @test result.logical_mask == exact["logical_mask"]
                @test posterior_error < 1e-8
            end
        end
    end
    println("CHECK complete: one history; loading/cache and chi960 exact checks passed")
    return nothing
end

abspath(PROGRAM_FILE) == (@__FILE__) && check_real_late_loading()
