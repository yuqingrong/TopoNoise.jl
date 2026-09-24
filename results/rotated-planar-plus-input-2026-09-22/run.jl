# Matched native product inputs and all noisy H/CNOT gates.
# Resample both changed As series; reuse only the verified unchanged Bp series.
using Dates, Random, Serialization, SHA, TOML, Test
using TopoNoise
using PythonCall

const RUN_DIR = @__DIR__
const REPO_DIR = normpath(joinpath(RUN_DIR, "..", ".."))
const ARCHIVE_DIR = joinpath(REPO_DIR, "results", "rotated-planar-native-logical-errors-2026-09-22")
const ARCHIVE_PATH = joinpath(ARCHIVE_DIR, "native_logical_errors-raw.jls")
const STEM = "plus_input_logical_errors"
const DISTANCES = [9, 11, 13, 15]
const ERROR_RATES = collect(0.0:0.002:0.1)
const SHOTS = 50_000
const BATCH_SIZE = 10_000
const RUN_SEED = 1235
const BOOTSTRAP_REPLICATES = 1_000
const SOURCE_FILES = vcat(["Project.toml", "CondaPkg.toml", "src/TopoNoise.jl"],
    [relpath(path, REPO_DIR) for path in
     readdir(joinpath(REPO_DIR, "src", "rotated_planar"); join=true)
     if endswith(path, ".jl")])
const SOURCE_HASHES = Dict(relative => bytes2hex(sha256(read(joinpath(REPO_DIR, relative))))
                           for relative in SOURCE_FILES)
const ARCHIVE_METADATA = TOML.parsefile(joinpath(ARCHIVE_DIR, "metadata.toml"))
const ARCHIVE = deserialize(ARCHIVE_PATH)
const REUSED_KEYS = ((:bp, :x_only), (:bp, :z_only))

# Reuse requires identical geometry, sampling and decoding implementations.
# The As schedule changed; compare the reused Bp gate schedules separately.
for file in ("geometry.jl", "decoder.jl")
    relative = "src/rotated_planar/" * file
    @assert SOURCE_HASHES[relative] == ARCHIVE_METADATA["source_sha256"][relative]
end
# Only the Yao oracle and the leading docstring changed in noise.jl. The
# batch sampler lives in decoder.jl, whose entire source is identical above.
# Check the exact executable noise definitions used by it as well.
function noise_sampling_source(path)
    source = read(path, String)
    return first(split(last(split(source, "struct CircuitPauliNoise"; limit=2)),
                       "function _apply_yao_encoder_and_faults!"; limit=2))
end
@assert noise_sampling_source(joinpath(REPO_DIR, "src/rotated_planar/noise.jl")) ==
        noise_sampling_source(joinpath(ARCHIVE_DIR, "source_snapshot/src/rotated_planar/noise.jl"))
module ArchivedEncoder
using Random, SparseArrays
import Yao
const SOURCE = joinpath(@__DIR__, "..", "rotated-planar-native-logical-errors-2026-09-22",
                        "source_snapshot", "src", "rotated_planar")
for file in ("geometry.jl", "encoder.jl")
    include(joinpath(SOURCE, file))
end
end
function gate_signature(encoder)
    return [(layer.block_kind, layer.block_order, copy(layer.active_qubits),
             [(op.gate, copy(op.qubits)) for op in layer.operations])
            for layer in encoder.layers]
end

const SCAN_RNG = MersenneTwister(RUN_SEED)
const ROOT_SEED = rand(SCAN_RNG, UInt64)
const SCHEDULES = Dict{String,Any}[]
@testset "Native inputs, all noisy gates matched, and Bp archive compatibility" begin
    @test ARCHIVE.distances == DISTANCES
    @test ARCHIVE.error_rates == ERROR_RATES
    @test ARCHIVE.shots == SHOTS
    @test ARCHIVE.batch_size == BATCH_SIZE
    @test ARCHIVE.master_seed == RUN_SEED
    @test ARCHIVE.clock === :gate_layer
    @test ARCHIVE.boundary_orientation === :x_ns
    @test pyconvert(String, pyimport("pymatching").__version__) ==
          ARCHIVE_METADATA["pymatching_version"]
    for d in DISTANCES, construction in (:as, :bp)
        state = construction === :as ? :plus : :zero
        encoder = rotated_planar_encoder(RotatedPlanarCode(d); construction)
        archived = ArchivedEncoder.rotated_planar_encoder(
            ArchivedEncoder.RotatedPlanarCode(d); construction, logical_state=state)
        @test encoder.logical_state === state
        @test encoder.input_state === (construction === :as ? :plus : :zero)
        @test verify_encoder_tableau(encoder)
        if construction === :bp
            @test gate_signature(encoder) == gate_signature(archived)
        else
            @test gate_signature(encoder) != gate_signature(archived)
            bp = rotated_planar_encoder(RotatedPlanarCode(d); construction=:bp)
            rotation = [(column - 1) * d + d + 1 - row
                        for row in 1:d for column in 1:d]
            ops(e) = [op for layer in e.layers for op in layer.operations]
            @test [(op.gate, op.qubits) for op in ops(encoder)] ==
                  [(op.gate, op.gate === :CNOT ? reverse(rotation[op.qubits]) : rotation[op.qubits])
                   for op in ops(bp)]
            as_blocks = filter(b -> b.kind === :plaquette, encoder.blocks)
            bp_blocks = filter(b -> b.kind === :plaquette, bp.blocks)
            @test [sort(b.source_support) for b in as_blocks] ==
                  [sort(rotation[b.source_support]) for b in bp_blocks]
        end
        @test all(b.source_check ∉ (:logical_parity, :logical_conversion) for b in encoder.blocks)
        operations = [only(layer.operations) for layer in encoder.layers]
        push!(SCHEDULES, Dict("distance" => d, "construction" => string(construction),
            "logical_state" => string(state), "input_state" => string(encoder.input_state), "gates" => length(operations),
            "hadamards" => count(op.gate === :H for op in operations),
            "cnots" => count(op.gate === :CNOT for op in operations),
            "fault_locations" => sum(length(op.qubits) for op in operations)))
    end
    for key in REUSED_KEYS
        scan = comparison_series(ARCHIVE, key...)
        @test scan.logical_state === (key[1] === :as ? :plus : :zero)
        @test scan.logical_observable === (key[2] === :x_only ? :logical_x : :logical_z)
        @test scan.series_seed == TopoNoise._comparison_series_seed(ROOT_SEED, key)
        @test all(p.logical_state === scan.logical_state && p.shots == SHOTS for p in scan.points)
    end
end
flush(stdout)

for relative in SOURCE_FILES
    destination = joinpath(RUN_DIR, "source_snapshot", relative)
    mkpath(dirname(destination))
    if isfile(destination)
        @assert bytes2hex(sha256(read(destination))) == SOURCE_HASHES[relative]
    else
        cp(joinpath(REPO_DIR, relative), destination)
    end
end
const METADATA_PATH = joinpath(RUN_DIR, "metadata.toml")
const METADATA = isfile(METADATA_PATH) ? TOML.parsefile(METADATA_PATH) : Dict{String,Any}(
    "started_utc" => string(now(UTC)), "julia_version" => string(VERSION),
    "pymatching_version" => pyconvert(String, pyimport("pymatching").__version__),
    "logical_state" => "native", "as_state" => "plus", "bp_state" => "zero",
    "metric" => "channel_logical", "clock" => "gate_layer", "boundary_orientation" => "x_ns",
    "seed" => RUN_SEED, "distances" => DISTANCES, "error_rates" => ERROR_RATES,
    "shots_per_point" => SHOTS, "batch_size" => BATCH_SIZE,
    "bootstrap_replicates" => BOOTSTRAP_REPLICATES,
    "noise" => "Independent selected Pauli after each representative H and each CNOT operand; ideal native product inputs and final syndrome, no idle/readout noise.",
    "decoder" => "Independent uniform-weight CSS PyMatching decoders using final syndrome only.",
    "as_preparation" => "Ideal all-plus input; H on fresh representatives followed by incoming CNOT trees; no logical-parity or basis-conversion block.",
    "cnot_matching" => "As is the clockwise 90-degree rotation of Bp with every control/target exchanged and the complete execution order preserved, including boundary checks.",
    "initialization_noise" => "Ideal all-plus input for As and all-zero input for Bp. Every noisy H and CNOT is paired under rotation and X/Z exchange; the native experiments are fully basis-dual.",
    "new_series" => ["as/x_only/plus", "as/z_only/plus"],
    "reused_series" => ["bp/x_only/zero", "bp/z_only/zero"],
    "archive_path" => relpath(ARCHIVE_PATH, REPO_DIR),
    "archive_sha256" => bytes2hex(sha256(read(ARCHIVE_PATH))),
    "archive_compatibility_verified" => true, "source_sha256" => SOURCE_HASHES,
    "schedules" => SCHEDULES,
)
@assert METADATA["source_sha256"] == SOURCE_HASHES
@assert METADATA["archive_sha256"] == bytes2hex(sha256(read(ARCHIVE_PATH)))
function save_metadata()
    open(METADATA_PATH, "w") do io
        TOML.print(io, METADATA; sorted=true)
    end
end
function checkpoint(path, object)
    serialize(path * ".tmp", object)
    mv(path * ".tmp", path; force=true)
end
save_metadata()

function sample_as(channel)
    observable = channel === :x_only ? :logical_x : :logical_z
    series_seed = TopoNoise._comparison_series_seed(ROOT_SEED, (:as, channel))
    rng = MersenneTwister(series_seed)
    points = LogicalFailurePoint[]
    for d in DISTANCES
        path = joinpath(RUN_DIR, "as-plus-$(channel)-d$(d).jls")
        if isfile(path)
            saved = deserialize(path)
            @assert saved.source_sha256 == SOURCE_HASHES
            scan, rng = saved.scan, saved.rng
            println("Restored As |+_L>, $(channel) noise, d=$d")
        else
            println("Sampling As |+_L>, $(channel) noise, d=$d: 51 points × 50,000 shots")
            flush(stdout)
            scan = open(joinpath(RUN_DIR, "as-plus-$(channel)-d$(d)-progress.log"), "w") do io
                scan_channel_logical_failure(rng;
                    distances=[d], error_rates=ERROR_RATES, construction=:as,
                    error_channel=channel, failure_metric=:channel_logical,
                    clock=:gate_layer, boundary_orientation=:x_ns,
                    shots=SHOTS, batch_size=BATCH_SIZE, seed=RUN_SEED,
                    series_seed, progress_io=io)
            end
            checkpoint(path, (scan=scan, rng=copy(rng), source_sha256=SOURCE_HASHES))
        end
        @assert scan.logical_state === :plus && scan.logical_observable === observable
        if channel === :x_only
            @assert all(p.logical_z_failures == 0 && p.state_failures == 0 for p in scan.points)
        else
            @assert all(p.logical_x_failures == 0 &&
                        p.state_failures == p.logical_z_failures for p in scan.points)
        end
        @assert first(scan.points).any_logical_failures == 0
        append!(points, scan.points)
        println("Saved d=$d")
        flush(stdout)
    end
    return ChannelFailureScan(:as, channel, observable, :plus, :x_ns, :gate_layer,
        DISTANCES, ERROR_RATES, SHOTS, BATCH_SIZE, RUN_SEED, series_seed, points)
end

const SERIES = [sample_as(:x_only), sample_as(:z_only),
                [comparison_series(ARCHIVE, key...) for key in REUSED_KEYS]...]
const RAW = ConstructionChannelComparison(:native, :x_ns, :gate_layer,
    DISTANCES, ERROR_RATES, SHOTS, BATCH_SIZE, RUN_SEED, SERIES)
checkpoint(joinpath(RUN_DIR, "$(STEM)-raw.jls"), RAW)
write_comparison_raw_csv(joinpath(RUN_DIR, "$(STEM)-raw.csv"),
    FittedConstructionChannelComparison(RAW, ThresholdFit[]))
METADATA["simulation_completed_utc"] = string(now(UTC))
METADATA["new_points"] = 408
METADATA["new_shots"] = 20_400_000
METADATA["reused_points"] = 408
METADATA["reused_shots"] = 20_400_000
save_metadata()
println("All 816 points saved. Local joint scaling fits are computed separately in Python.")
flush(stdout)

@assert sum(length(s.points) for s in SERIES) == 816
@assert all(p.shots == SHOTS for s in SERIES for p in s.points)
METADATA["source_unchanged_during_run"] = all(
    bytes2hex(sha256(read(joinpath(REPO_DIR, relative)))) == digest
    for (relative, digest) in SOURCE_HASHES)
@assert METADATA["source_unchanged_during_run"]
METADATA["completed_utc"] = string(now(UTC))
save_metadata()
println("Completed: $(joinpath(RUN_DIR, "$(STEM)-raw.csv"))")
