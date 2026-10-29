#!/usr/bin/env julia
"""Deterministically generate authenticated Julia records from frozen seeds."""

using Random, TOML, LinearAlgebra

include(joinpath(@__DIR__, "..", "run.jl"))
const Runner = TNMLRepairedRunner

runtime_identity() = Dict{String,Any}(
    "julia_version" => string(VERSION),
    "blas_vendor" => string(BLAS.vendor()),
    "cpu_model" => string(Sys.CPU_NAME),
    "julia_threads" => Threads.nthreads(),
    "blas_threads" => BLAS.get_num_threads(),
)

function atomic_toml(path, document)
    mkpath(dirname(path))
    temporary, io = mktemp(dirname(path))
    try
        TOML.print(io, Runner.toml_value(document); sorted=true)
        flush(io)
        close(io)
        mv(temporary, path; force=false)
    finally
        isopen(io) && close(io)
        isfile(temporary) && rm(temporary)
    end
end

function main(args=ARGS)
    length(args) >= 6 || throw(ArgumentError(
        "expected generate L p task_id output.toml history_seed..."))
    args[1] == "generate" || throw(ArgumentError("mode must be generate"))
    L = parse(Int, args[2])
    p_key = args[3]
    p = parse(Float64, p_key)
    task_id = parse(Int, args[4])
    output = abspath(args[5])
    seeds = parse.(Int, args[6:end])
    isempty(seeds) && throw(ArgumentError("at least one history seed is required"))
    length(unique(seeds)) == length(seeds) || throw(ArgumentError("history seeds repeat"))
    Threads.nthreads() == 1 || throw(ArgumentError("corpus generator requires one Julia thread"))
    blas_threads = parse(Int, get(ENV, "TNML_BLAS_THREADS", "4"))
    blas_threads == 4 || throw(ArgumentError("corpus generation uses the BLAS=4 primary profile"))
    BLAS.set_num_threads(blas_threads)
    Runner.load_backend!()
    backend = Base.invokelatest(getproperty, Runner, :RepairedSeedTN)
    context = Base.invokelatest(getproperty(backend, :build_context), L)
    records = Dict{String,Any}[]
    for (shot_index, seed) in enumerate(seeds)
        rng = Random.Xoshiro(seed)
        record = Base.invokelatest(getproperty(backend, :sample_record), context, p, rng)
        document = Runner.record_document(record, shot_index)
        payload = Dict(key => value for (key, value) in document if key != "record_id")
        Runner.identity_hash(payload) == document["record_id"] ||
            error("generated record failed established sorted-TOML identity")
        push!(records, document)
    end
    atomic_toml(output, Dict{String,Any}(
        "schema_version" => 1,
        "kind" => "adaptive_chi_julia_corpus_task",
        "L" => L,
        "p" => p,
        "p_key" => p_key,
        "task_id" => task_id,
        "seed_mode" => "fresh_xoshiro_per_history_seed",
        "rng_algorithm" => "Random.Xoshiro",
        "runtime_identity" => runtime_identity(),
        "history_seeds" => seeds,
        "records" => records,
    ))
    return nothing
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
