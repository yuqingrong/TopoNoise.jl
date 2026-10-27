# Fixed study orchestration. The simulation and checkpoint implementation is unchanged.
include(joinpath(@__DIR__, "..", "toric3d-h6-full-prefix-2026-10-15", "run_study.jl"))
using .H6FullPrefixStudy
import TOML, SHA

const STUDY = @__DIR__
const JOBS = TOML.parsefile(joinpath(STUDY, "jobs.toml"))["jobs"]

if ARGS == ["check"]
    generated = Dict{String,Any}()
    all_new = String[]
    for job in JOBS
        seeds = H6FullPrefixStudy.hexseed.(H6FullPrefixStudy.point_seeds(
            job["seed"], length(job["sizes"]) * length(job["rates"])))
        generated[job["name"]] = seeds
        append!(all_new, seeds)
    end
    length(all_new) == length(unique(all_new)) || error("new point seeds overlap")
    previous = Set{String}()
    source_metadata = TOML.parsefile(joinpath(STUDY, "benchmark-x", "metadata.toml"))
    source_metadata["source_hashes"] == H6FullPrefixStudy.source_hashes() || error("source changed since benchmark")
    source_metadata["runtime"] == H6FullPrefixStudy.runtime_metadata() || error("runtime changed since benchmark")
    for (directory, _, files) in walkdir(joinpath(H6FullPrefixStudy.ROOT, "results"))
        startswith(directory, STUDY) && !occursin("benchmark-", directory) && continue
        "metadata.toml" in files || continue
        metadata = TOML.parsefile(joinpath(directory, "metadata.toml"))
        for seed in get(metadata, "point_seeds_hex", String[])
            push!(previous, seed)
        end
    end
    isempty(intersect(Set(all_new), previous)) || error("new point seed overlaps historical/timing stream")
    report = Dict("status" => "passed", "generated_point_seeds" => generated,
        "new_points" => length(all_new), "distinct_previous_point_seeds_checked" => length(previous),
        "source_and_runtime_match_benchmark" => true, "checked_at_unix" => time())
    open(joinpath(STUDY, "verification", "preflight.toml"), "w") do io
        TOML.print(io, report; sorted=true)
    end
    println("Preflight passed for ", length(all_new), " fresh point streams; checked ", length(previous), " prior streams.")
    exit(0)
end

length(ARGS) == 1 || error("Specify check or worker 1/2")
worker = parse(Int, only(ARGS))
worker in (1, 2) || error("Worker must be 1 or 2")
isfile(joinpath(STUDY, "verification", "preflight.toml")) || error("Run preflight first")
for job in JOBS
    job["worker"] == worker || continue
    directory = joinpath(STUDY, job["name"])
    arguments = ["--channel", job["channel"], "--output", directory,
        "--sizes", join(job["sizes"], ','), "--rates", join(job["rates"], ','),
        "--shots", string(job["shots"]), "--seed", string(job["seed"]),
        "--batch-size", string(job["batch_size"])]
    isfile(joinpath(directory, "checkpoint.toml")) && push!(arguments, "--resume")
    mkpath(directory)
    println("START worker=", worker, " job=", job["name"], " unix=", time()); flush(stdout)
    status = open(joinpath(directory, "run.log"), "a") do io
        H6FullPrefixStudy.main(arguments; io, error_io=io)
    end
    status == 0 || error("Simulation job failed: " * job["name"])
    println("DONE worker=", worker, " job=", job["name"], " unix=", time()); flush(stdout)
    GC.gc()
end
println("ALL JOBS COMPLETE worker=", worker)
