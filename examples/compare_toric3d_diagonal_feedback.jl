module Toric3DDiagonalFeedbackExample

using Random
using TopoNoise
import TOML

const USAGE = """
Usage: julia --project=. examples/compare_toric3d_diagonal_feedback.jl [options]
  --sizes LIST          cubic L values or stepsxLyxLz triples
  --rates LIST          physical fault probabilities
  --channels LIST       x_only, z_only, or x_only,z_only (default x_only)
  --shots N             paired scientific triplets per point (default 200)
  --batch-size N        atomic checkpoint batch size (default 64)
  --seed N              master seed and recorded provenance (default 1234)
  --x-routing ROUTE     one_step or full_path (default one_step)
  --output PATH         output directory
  --resume              resume an exactly matching checkpoint
  --timing-shots N      independent timing triplets per channel (default 20)
  --pilot-only          write timing-pilot.toml and stop before the scan
  --help

Without an explicit grid, x_only runs the approved 200-shot X scout and
z_only runs the approved 200-shot Z scout. Both channels require explicit
--sizes and --rates. Warmup and timing records are never included in the
scientific counts.
"""

const _X_SIZES = [3, 5, 7, 9]
const _X_RATES = [0.0, 0.0001, 0.0003, 0.001, 0.003, 0.005, 0.0075,
    0.01, 0.015, 0.02, 0.03, 0.05, 0.075, 0.10]
const _Z_SIZES = [5, 9, 13, 17]
const _Z_RATES = [0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006,
    0.0075, 0.01, 0.015, 0.02, 0.03]

const _VALUE_OPTIONS = Set(["--sizes", "--rates", "--channels", "--shots",
    "--batch-size", "--seed", "--x-routing", "--output", "--timing-shots"])

function _parse_sizes(text::AbstractString)
    isempty(strip(text)) && throw(ArgumentError("sizes must not be empty"))
    values = Union{Int,NTuple{3,Int}}[]
    for raw in split(text, ',')
        token = strip(raw)
        isempty(token) && throw(ArgumentError("sizes contain an empty entry"))
        fields = split(lowercase(token), 'x')
        dimensions = try
            if length(fields) == 1
                size = parse(Int, only(fields))
                (size, size, size)
            elseif length(fields) == 3
                Tuple(parse.(Int, fields))
            else
                throw(ArgumentError("expected L or stepsxLyxLz"))
            end
        catch error
            error isa ArgumentError && contains(sprint(showerror, error),
                "expected L or stepsxLyxLz") && rethrow()
            throw(ArgumentError("could not parse sizes: $(sprint(showerror, error))"))
        end
        dimensions[1] >= 1 && dimensions[2] >= 2 && dimensions[3] >= 2 ||
            throw(ArgumentError("sizes require steps >= 1 and Ly,Lz >= 2"))
        push!(values, dimensions)
    end
    length(unique(values)) == length(values) ||
        throw(ArgumentError("sizes must be unique"))
    all(values[index - 1] < values[index] for index in 2:length(values)) ||
        throw(ArgumentError("sizes must be strictly increasing"))
    return NTuple{3,Int}[value for value in values]
end

function _parse_rates(text::AbstractString)
    isempty(strip(text)) && throw(ArgumentError("rates must not be empty"))
    rates = try
        parse.(Float64, strip.(split(text, ',')))
    catch error
        throw(ArgumentError("could not parse rates: $(sprint(showerror, error))"))
    end
    all(rate -> isfinite(rate) && 0 <= rate < 0.5, rates) ||
        throw(ArgumentError("rates must lie in [0, 0.5)"))
    all(rates[index - 1] < rates[index] for index in 2:length(rates)) ||
        throw(ArgumentError("rates must be strictly increasing"))
    return rates
end

function _parse_channels(text::AbstractString)
    isempty(strip(text)) && throw(ArgumentError("channels must not be empty"))
    channels = Symbol.(strip.(split(text, ',')))
    all(channel -> channel in (:x_only, :z_only), channels) ||
        throw(ArgumentError("channels may contain only x_only and z_only"))
    length(unique(channels)) == length(channels) ||
        throw(ArgumentError("channels must be unique"))
    return channels
end

function _positive_integer(text, name)
    value = try
        parse(Int, text)
    catch error
        throw(ArgumentError("could not parse $name: $(sprint(showerror, error))"))
    end
    value > 0 || throw(ArgumentError("$name must be positive"))
    return value
end

function _parse_arguments(arguments)
    values = Dict{String,String}()
    resume = false
    pilot_only = false
    index = 1
    while index <= length(arguments)
        option = arguments[index]
        if option == "--resume"
            resume = true
            index += 1
            continue
        elseif option == "--pilot-only"
            pilot_only = true
            index += 1
            continue
        end
        option in _VALUE_OPTIONS || throw(ArgumentError("unknown option $option"))
        index < length(arguments) || throw(ArgumentError("$option needs a value"))
        value = arguments[index + 1]
        startswith(value, "--") && throw(ArgumentError("$option needs a value"))
        values[option] = value
        index += 2
    end

    channels = _parse_channels(get(values, "--channels", "x_only"))
    sizes_explicit = haskey(values, "--sizes")
    rates_explicit = haskey(values, "--rates")
    if length(channels) == 2 && !(sizes_explicit && rates_explicit)
        throw(ArgumentError("both channels require explicit --sizes and --rates"))
    end
    defaults_sizes = channels == [:z_only] ? _Z_SIZES : _X_SIZES
    defaults_rates = channels == [:z_only] ? _Z_RATES : _X_RATES
    sizes = sizes_explicit ? _parse_sizes(values["--sizes"]) : copy(defaults_sizes)
    error_rates = rates_explicit ? _parse_rates(values["--rates"]) : copy(defaults_rates)
    shots = _positive_integer(get(values, "--shots", "200"), "shots")
    batch_size = _positive_integer(get(values, "--batch-size", "64"),
        "batch size")
    timing_shots = _positive_integer(get(values, "--timing-shots", "20"),
        "timing shots")
    seed = try
        parse(Int, get(values, "--seed", "1234"))
    catch error
        throw(ArgumentError("could not parse seed: $(sprint(showerror, error))"))
    end
    seed >= 0 || throw(ArgumentError("seed must be nonnegative"))
    x_routing = Symbol(get(values, "--x-routing", "one_step"))
    x_routing in (:one_step, :full_path) ||
        throw(ArgumentError("x-routing must be one_step or full_path"))
    output = get(values, "--output",
        "results/toric3d-diagonal-feedback-2026-10-12/scout")
    isempty(strip(output)) && throw(ArgumentError("output must not be empty"))
    return (; sizes, error_rates, channels, shots, batch_size, seed, x_routing,
        output, resume, timing_shots, pilot_only, sizes_explicit,
        rates_explicit)
end

function _write_timing(path, data)
    temporary = path * ".tmp"
    try
        open(temporary, "w") do io
            TOML.print(io, data; sorted=true)
        end
        mv(temporary, path; force=true)
    finally
        isfile(temporary) && rm(temporary; force=true)
    end
    return path
end

function _timing_dimensions(options, channel)
    if options.sizes_explicit
        value = last(options.sizes)
        return value isa Integer ? (value, value, value) : value
    end
    return channel === :x_only ? (9, 9, 9) : (17, 17, 17)
end

function _timing_probability(options, channel)
    target = channel === :x_only ? 0.015 : 0.006
    options.rates_explicit || return target
    return options.error_rates[argmin(abs.(options.error_rates .- target))]
end

function _warmup(options)
    rng = Xoshiro(UInt64(options.seed) ⊻ UInt64(0x7761726d75700001))
    for channel in options.channels
        code = Toric3DEncodedSlice(1, 2, 2)
        encoder = toric3d_encoded_encoder(code)
        if channel === :x_only
            family = build_toric3d_x_feedback_decoders(code; p=0.01)
            faults = sample_toric3d_faults(rng, encoder; p_x=0.01, p_z=0.0)
            for protocol in (:none, :slice_local, :full_prefix)
                replay_toric3d_x_strength(encoder, faults, family;
                    protocol, routing=options.x_routing)
            end
        else
            family = build_toric3d_z_feedback_decoders(code)
            faults = sample_toric3d_diagonal_z_faults(rng, encoder; p=0.01)
            for protocol in (:none, :slice_local, :full_prefix)
                replay_toric3d_z_strength(encoder, faults, family; protocol)
            end
        end
    end
    return nothing
end

function _timing_pilot(options; io::IO)
    println(io, "Warmup: independent records outside the scientific scan...")
    flush(io)
    warmup_seconds = @elapsed _warmup(options)
    protocol_seconds = Dict(protocol => 0.0
        for protocol in ("none", "slice_local", "full_prefix"))
    construction_by_channel = Dict{String,Float64}()
    sampling_seconds = 0.0
    targets = Dict{String,Any}[]
    timing_rng = Xoshiro(UInt64(options.seed) ⊻ UInt64(0x74696d696e670001))

    println(io, "Timing pilot: $(options.timing_shots) independent triplets per channel...")
    flush(io)
    for channel in options.channels
        dimensions = _timing_dimensions(options, channel)
        p = _timing_probability(options, channel)
        code = Toric3DEncodedSlice(dimensions...)
        encoder = toric3d_encoded_encoder(code)
        family = nothing
        build_seconds = @elapsed family = channel === :x_only ?
            build_toric3d_x_feedback_decoders(code; p) :
            build_toric3d_z_feedback_decoders(code)
        construction_by_channel[string(channel)] = build_seconds
        push!(targets, Dict{String,Any}(
            "channel" => string(channel), "dimensions" => collect(dimensions),
            "p" => p, "decoder_construction_seconds" => build_seconds))

        for _ in 1:options.timing_shots
            faults = nothing
            sampling_seconds += @elapsed faults = channel === :x_only ?
                sample_toric3d_faults(timing_rng, encoder; p_x=p, p_z=0.0) :
                sample_toric3d_diagonal_z_faults(timing_rng, encoder; p)
            for protocol in (:none, :slice_local, :full_prefix)
                elapsed = @elapsed if channel === :x_only
                    replay_toric3d_x_strength(encoder, faults, family;
                        protocol, routing=options.x_routing)
                else
                    replay_toric3d_z_strength(encoder, faults, family; protocol)
                end
                protocol_seconds[string(protocol)] += elapsed
            end
        end
    end

    decoder_construction_seconds = sum(values(construction_by_channel))
    timed_triplets = options.timing_shots * length(options.channels)
    per_triplet = (sampling_seconds + sum(values(protocol_seconds))) /
        timed_triplets
    point_count = length(options.sizes) * length(options.error_rates) *
        length(options.channels)
    estimated_scan_seconds = decoder_construction_seconds * point_count +
        per_triplet * options.shots * point_count
    println(io, "Warmup ", round(warmup_seconds; digits=2),
        " s; decoder construction ", round(decoder_construction_seconds; digits=3),
        " s; timed triplets ", timed_triplets, ".")
    println(io, "Rough scan estimate: ",
        round(estimated_scan_seconds / 60; digits=1), " min.")
    flush(io)
    return Dict{String,Any}(
        "format_version" => 1,
        "warmup_seconds" => warmup_seconds,
        "warmup_shots" => length(options.channels),
        "timing_shots" => options.timing_shots,
        "timed_triplets" => timed_triplets,
        "timing_shots_included_in_production" => false,
        "production_shots_per_point" => options.shots,
        "production_points" => point_count,
        "decoder_construction_seconds" => decoder_construction_seconds,
        "decoder_construction_by_channel_seconds" => construction_by_channel,
        "sampling_seconds" => sampling_seconds,
        "protocol_seconds" => protocol_seconds,
        "estimated_scan_seconds" => estimated_scan_seconds,
        "targets" => targets,
        "x_routing" => string(options.x_routing),
        "timing_policy" =>
            "one immutable sampled record replayed independently by each protocol; all records excluded from production")
end

function main(arguments=ARGS; io::IO=stdout, error_io::IO=stderr)::Int
    if "--help" in arguments
        print(io, USAGE)
        return 0
    end
    try
        options = _parse_arguments(arguments)
        mkpath(options.output)
        timing = _timing_pilot(options; io)
        _write_timing(joinpath(options.output, "timing-pilot.toml"), timing)
        if options.pilot_only
            println(io, "Pilot-only run complete; no scientific shots were sampled.")
            flush(io)
            return 0
        end
        options.resume && println(io,
            "Resuming the exactly matching checkpoint in $(options.output)...")
        scan = compare_toric3d_diagonal_feedback(Xoshiro(options.seed);
            sizes=options.sizes, error_rates=options.error_rates,
            channels=options.channels, shots=options.shots, seed=options.seed,
            batch_size=options.batch_size, x_routing=options.x_routing,
            output_dir=options.output, resume=options.resume, progress_io=io)
        paths = save_toric3d_diagonal_feedback(scan, options.output)
        println(io, "Saved: ", paths)
        flush(io)
        return 0
    catch error
        error isa InterruptException && rethrow()
        print(error_io, "error: ")
        showerror(error_io, error)
        println(error_io)
        print(error_io, USAGE)
        return 1
    end
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(Toric3DDiagonalFeedbackExample.main())
end
