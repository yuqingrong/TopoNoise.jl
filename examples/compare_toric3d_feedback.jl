module Toric3DFeedbackExample

using Random

const DEFAULT_RATES = [0.0, 0.001, 0.003, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15]
const USAGE = """
Usage: julia --project=. examples/compare_toric3d_feedback.jl [options]

Compare paired final-only and layer-feedback recovery on cubic 3D toric slabs.
Each size L means (L,L,L). The pilot runs at p=0.03 on the largest size.

Options:
  --sizes LIST        comma-separated cube sizes (default: 3,5,7)
  --error-rates LIST  comma-separated rates (default: 0,.001,.003,.01,.02,.03,.05,.075,.1,.15)
  --shots INTEGER     shots per full-scan point (default: 1000)
  --seed INTEGER      nonnegative RNG seed (default: 1234)
  --output-dir PATH   artifact directory (default: results/toric3d-feedback-2026-10-11)
  --pilot-shots N     shots for the timing pilot (default: 10)
  --help              show this help
"""

function _list(value, parser, option)
    entries = split(value, ','; keepempty=true)
    (isempty(entries) || any(isempty, entries)) &&
        throw(ArgumentError("$option needs a comma-separated list without empty entries"))
    try
        return parser.(strip.(entries))
    catch
        throw(ArgumentError("$option has an invalid entry: $value"))
    end
end

function _integer(value, option)
    try
        return parse(Int, value)
    catch
        throw(ArgumentError("$option must be an integer"))
    end
end

function _parse_arguments(arguments)
    options = Dict{Symbol,Any}(
        :sizes => [3, 5, 7],
        :error_rates => copy(DEFAULT_RATES),
        :shots => 1000,
        :seed => 1234,
        :output_dir => joinpath("results", "toric3d-feedback-2026-10-11"),
        :pilot_shots => 10,
    )
    named = Dict(
        "--sizes" => :sizes,
        "--error-rates" => :error_rates,
        "--shots" => :shots,
        "--seed" => :seed,
        "--output-dir" => :output_dir,
        "--pilot-shots" => :pilot_shots,
    )
    index = 1
    while index <= length(arguments)
        option = arguments[index]
        option == "--help" && return nothing
        haskey(named, option) || throw(ArgumentError("unknown option: $option"))
        index < length(arguments) && !startswith(arguments[index + 1], "--") ||
            throw(ArgumentError("$option requires a value"))
        value = arguments[index + 1]
        key = named[option]
        options[key] = if key == :sizes
            _list(value, x -> parse(Int, x), option)
        elseif key == :error_rates
            _list(value, x -> parse(Float64, x), option)
        elseif key == :output_dir
            String(value)
        else
            _integer(value, option)
        end
        index += 2
    end
    sizes = options[:sizes]
    !isempty(sizes) && all(size -> size >= 2, sizes) && all(diff(sizes) .> 0) ||
        throw(ArgumentError("--sizes must contain strictly increasing integers at least 2"))
    rates = options[:error_rates]
    !isempty(rates) && all(rate -> isfinite(rate) && 0 <= rate < 0.5, rates) &&
        all(diff(rates) .> 0) ||
        throw(ArgumentError("--error-rates must be strictly increasing finite rates in [0,0.5)"))
    options[:shots] > 0 || throw(ArgumentError("--shots must be positive"))
    options[:pilot_shots] > 0 || throw(ArgumentError("--pilot-shots must be positive"))
    options[:seed] >= 0 || throw(ArgumentError("--seed must be nonnegative"))
    !isempty(strip(options[:output_dir])) || throw(ArgumentError("--output-dir must not be empty"))
    return options
end

function main(args=ARGS; io::IO=stdout, error_io::IO=stderr)
    options = try
        _parse_arguments(args)
    catch error
        error isa ArgumentError || rethrow()
        println(error_io, "error: ", sprint(showerror, error))
        println(error_io, USAGE)
        return 1
    end
    if options === nothing
        print(io, USAGE)
        return 0
    end

    # Keep --help and argument validation available before loading the package.
    try
        @eval import TopoNoise
    catch error
        println(error_io, "error: ", sprint(showerror, error))
        return 1
    end
    return Base.invokelatest(_run, options, io, error_io)
end

# The import above occurs after this function was defined. Invoke in the latest
# world so freshly loaded TopoNoise methods are visible on the first CLI run.
function _run(options, io::IO, error_io::IO)
    sizes = options[:sizes]
    pilot_size = maximum(sizes)
    seed = options[:seed]
    # The warmup and pilot streams cannot consume any full-scan draws.
    warmup_rng = MersenneTwister(UInt(seed) ⊻ UInt(0x6f70656e))
    pilot_rng = MersenneTwister(UInt(seed) ⊻ UInt(0x746f7269))
    try
        println(io, "Warming up at L=$pilot_size, p=0.03 (one paired shot)...")
        compilation_seconds = @elapsed TopoNoise.compare_toric3d_feedback(
            warmup_rng; sizes=[pilot_size], error_rates=[0.03], shots=1,
            output_dir=nothing)
        println(io, "Warmup including compilation: ", round(compilation_seconds; digits=2), " s")
        pilot_shots = options[:pilot_shots]
        println(io, "Timing pilot: $pilot_shots paired shots at L=$pilot_size, p=0.03...")
        pilot_seconds = @elapsed TopoNoise.compare_toric3d_feedback(
            pilot_rng; sizes=[pilot_size], error_rates=[0.03], shots=pilot_shots,
            output_dir=nothing)
        points = length(sizes) * length(options[:error_rates])
        estimate = pilot_seconds * options[:shots] * points / pilot_shots
        println(io, "Pilot runtime: ", round(pilot_seconds; digits=2), " s")
        println(io, "Rough full-scan estimate: ", round(estimate / 60; digits=1),
            " min for $points size/rate combinations ($((2 * points)) paired channel points; ",
            "rates and decoder cost vary).")

        println(io, "Running full scan; completed points are saved in ", options[:output_dir])
        comparison = TopoNoise.compare_toric3d_feedback(
            MersenneTwister(seed);
            sizes=sizes, error_rates=options[:error_rates], shots=options[:shots],
            seed=seed, output_dir=options[:output_dir], progress_io=io)
        paths = TopoNoise.save_toric3d_comparison(
            comparison, options[:output_dir]; basename="toric3d-feedback")
        println(io, "Saved comparison: ", paths)
    catch error
        println(error_io, "error: ", sprint(showerror, error))
        return 1
    end
    return 0
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(Toric3DFeedbackExample.main())
end
