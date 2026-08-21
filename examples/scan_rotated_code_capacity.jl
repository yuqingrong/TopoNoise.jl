module ScanRotatedCodeCapacityExample

using TopoNoise
import CairoMakie
import Random

const DEFAULT_OUTPUT_DIR = normpath(joinpath(
    @__DIR__, "..", "results", "rotated-code-capacity"))

const USAGE = """
Usage: julia --project=. examples/scan_rotated_code_capacity.jl \\
  --p-min P --p-max P --p-step DP [options]

Required:
  --p-min P          first data-X error rate
  --p-max P          last data-X error rate
  --p-step DP        positive grid spacing that exactly reaches p-max

Options:
  --sizes LIST       comma-separated rotated planar distances (default: 3,5,7)
  --shots N          samples per (distance,p) point (default: 10000)
  --batches N        bootstrap batches, at least 2 (default: min(100,shots))
  --bootstrap N      bootstrap replicates (default: 2000)
  --confidence X     confidence level (default: 0.95)
  --seed N           random seed (default: 1234)
  --output-dir PATH  output directory (default: results/rotated-code-capacity)
"""

const _KNOWN_OPTIONS = Set((
    "--p-min", "--p-max", "--p-step", "--sizes", "--shots",
    "--batches", "--bootstrap", "--confidence", "--seed", "--output-dir"))

function _option_dictionary(args)
    iseven(length(args)) || throw(ArgumentError(
        "every option must be followed by a value"))
    options = Dict{String,String}()
    for index in 1:2:length(args)
        option, value = args[index], args[index + 1]
        option in _KNOWN_OPTIONS || throw(ArgumentError("unknown option: $option"))
        haskey(options, option) && throw(ArgumentError(
            "option appears more than once: $option"))
        options[option] = value
    end
    return options
end

function _required(options, name)
    haskey(options, name) || throw(ArgumentError(
        "required option is missing: $name"))
    return options[name]
end

function _parse_number(::Type{T}, value, name) where {T<:Real}
    return try
        parse(T, value)
    catch
        throw(ArgumentError("$name must be a valid $(T) value"))
    end
end

function _positive_integer(options, name, default)
    value = _parse_number(Int, get(options, name, string(default)), name)
    value > 0 || throw(ArgumentError("$name must be positive"))
    return value
end

function _parse_sizes(value)
    isempty(value) && throw(ArgumentError("--sizes must not be empty"))
    sizes = try
        parse.(Int, split(value, ','))
    catch
        throw(ArgumentError("--sizes must be comma-separated integers"))
    end
    length(sizes) >= 2 || throw(ArgumentError(
        "--sizes must contain at least two rotated planar distances"))
    all(>=(3), sizes) || throw(ArgumentError(
        "--sizes entries must be at least 3"))
    all(sizes[index] < sizes[index + 1]
        for index in 1:(length(sizes) - 1)) || throw(ArgumentError(
            "--sizes must be strictly increasing"))
    return sizes
end

function _error_rate_grid(options)
    p_min = _parse_number(Float64, _required(options, "--p-min"), "--p-min")
    p_max = _parse_number(Float64, _required(options, "--p-max"), "--p-max")
    p_step = _parse_number(
        Float64, _required(options, "--p-step"), "--p-step")
    all(isfinite, (p_min, p_max, p_step)) || throw(ArgumentError(
        "data-X error-rate arguments must be finite"))
    0 <= p_min < p_max < 0.5 || throw(ArgumentError(
        "require 0 <= p-min < p-max < 0.5"))
    p_step > 0 || throw(ArgumentError("--p-step must be positive"))
    interval_count = (p_max - p_min) / p_step
    rounded_count = round(Int, interval_count)
    isapprox(interval_count, rounded_count; atol=1e-10, rtol=1e-10) ||
        throw(ArgumentError("--p-step must exactly reach --p-max"))
    rounded_count >= 1 || throw(ArgumentError(
        "--p-step must produce at least two data-X error rates"))
    rates = [p_min + index * p_step for index in 0:rounded_count]
    rates[end] = p_max
    return rates
end

function _parse_options(args)
    options = _option_dictionary(args)
    rates = _error_rate_grid(options)
    sizes = _parse_sizes(get(options, "--sizes", "3,5,7"))
    shots = _positive_integer(options, "--shots", 10_000)
    batches = _positive_integer(options, "--batches", min(100, shots))
    batches >= 2 || throw(ArgumentError(
        "--batches must be at least 2 for bootstrap crossings"))
    batches <= shots || throw(ArgumentError("--batches must not exceed --shots"))
    shots % batches == 0 || throw(ArgumentError(
        "--shots must be divisible by --batches"))
    bootstrap = _positive_integer(options, "--bootstrap", 2_000)
    seed = _parse_number(Int, get(options, "--seed", "1234"), "--seed")
    confidence = _parse_number(
        Float64, get(options, "--confidence", "0.95"), "--confidence")
    isfinite(confidence) && 0 < confidence < 1 || throw(ArgumentError(
        "--confidence must be strictly between 0 and 1"))
    output_dir = get(options, "--output-dir", DEFAULT_OUTPUT_DIR)
    isempty(output_dir) && throw(ArgumentError("--output-dir must not be empty"))
    return (; rates, sizes, shots, batches, bootstrap, seed, confidence,
            output_dir)
end

function _csv_value(value)
    ismissing(value) && return ""
    value isa Symbol && return string(value)
    return repr(value)
end

function _write_scan_csv(path, scan)
    fields = (
        :distance, :error_rate, :shots, :logical_x_failure_count,
        :logical_x_failure_rate, :logical_x_failure_se)
    open(path, "w") do io
        println(io, join((
            "L", "p", "shots", "logical_x_failure_count",
            "logical_x_failure_rate", "logical_x_failure_se", "batches", "seed"), ','))
        for point in scan.points
            values = (_csv_value(getfield(point, field)) for field in fields)
            println(io, join((values..., length(point.logical_x_failure_batches),
                              _csv_value(point.seed)), ','))
        end
    end
    return path
end

function _write_crossing_csv(path, crossings)
    fields = (
        :small_size, :large_size, :estimate, :ci_low, :ci_high,
        :confidence, :valid_bootstrap_fraction, :status)
    open(path, "w") do io
        println(io, join((
            "small_size", "large_size", "estimate", "ci_low", "ci_high",
            "confidence", "valid_bootstrap_fraction", "status"), ','))
        for crossing in crossings
            println(io, join(
                (_csv_value(getfield(crossing, field)) for field in fields), ','))
        end
    end
    return path
end

function _crossing_summary(crossing)
    sizes = "L=$(crossing.small_size)/$(crossing.large_size)"
    if crossing.status == :ok && !ismissing(crossing.estimate) &&
       !ismissing(crossing.ci_low) && !ismissing(crossing.ci_high)
        level = round(Int, 100 * crossing.confidence)
        return "$sizes: p_c = $(crossing.estimate) " *
               "[$(crossing.ci_low), $(crossing.ci_high)] ($level%)"
    elseif crossing.status == :unbracketed
        return "$sizes: unbracketed — extend p range"
    elseif crossing.status == :unstable
        return "$sizes: unstable — increase samples"
    end
    return "$sizes: $(crossing.status)"
end

function run(options; io::IO=stdout)
    rng = Random.MersenneTwister(options.seed)
    scan = scan_rotated_code_capacity(
        rng, options.sizes, options.rates;
        shots=options.shots, batches=options.batches, seed=options.seed,
        progress_io=io)
    crossings = estimate_rotated_code_crossings(
        rng, scan; bootstrap=options.bootstrap, confidence=options.confidence)
    figure = plot_rotated_code_capacity(scan, crossings)

    mkpath(options.output_dir)
    scan_path = abspath(joinpath(
        options.output_dir, "rotated_code_capacity_scan.csv"))
    crossing_path = abspath(joinpath(
        options.output_dir, "rotated_code_capacity_crossings.csv"))
    _write_scan_csv(scan_path, scan)
    _write_crossing_csv(crossing_path, crossings)
    outputs = [scan_path, crossing_path]
    for extension in ("svg", "pdf", "png")
        path = abspath(joinpath(
            options.output_dir, "rotated_code_capacity_scan.$extension"))
        CairoMakie.save(path, figure)
        push!(outputs, path)
    end

    for crossing in crossings
        println(io, _crossing_summary(crossing))
    end
    for path in outputs
        println(io, "wrote: $path")
    end
    return (; scan, crossings, figure, outputs)
end

function main(args=ARGS; io::IO=stdout, error_io::IO=stderr)
    if args == ["--help"]
        println(io, USAGE)
        return 0
    end
    options = try
        _parse_options(args)
    catch error
        error isa ArgumentError || rethrow()
        println(error_io, "error: ", error.msg)
        println(error_io, USAGE)
        return 1
    end
    run(options; io=io)
    return 0
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(ScanRotatedCodeCapacityExample.main())
end
