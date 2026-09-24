module RotatedPlanarScanExample

using Random
using TopoNoise

const USAGE = """
Usage: julia --project=. examples/scan_rotated_planar.jl [options]

Options:
  --distances LIST             comma-separated odd distances (default: 3,5,7)
  --error-rates LIST           comma-separated physical error rates
  --p-min VALUE                error-grid minimum (default: 0)
  --p-max VALUE                error-grid maximum (default: 0.1)
  --p-step VALUE               error-grid step (default: 0.01)
  --state NAME                 native, zero, one, plus, or minus (default: native: As plus, Bp zero)
  --construction NAME          bp or as (default: bp)
  --boundary-orientation NAME  x_ns or x_ew (default: x_ns)
  --clock NAME                 gate_layer, plaquette, or post_encoding (default: gate_layer)
  --shots INTEGER              shots per scan point (default: 10000)
  --batch-size INTEGER         estimator batch size (default: 10000)
  --seed INTEGER               RNG seed, or nothing (default: 1234)
  --output-dir PATH            artifact directory (default: results/rotated-planar)
  --basename NAME              artifact basename (default: rotated_planar_logical_failure)
  --help                       show this help
"""

const _VALUE_OPTIONS = (
    "--distances",
    "--error-rates",
    "--p-min",
    "--p-max",
    "--p-step",
    "--state",
    "--construction",
    "--boundary-orientation",
    "--clock",
    "--shots",
    "--batch-size",
    "--seed",
    "--output-dir",
    "--basename",
)
const _MAX_ERROR_GRID_POINTS = 1_000_000

function _comma_separated(value::AbstractString, parser, name::AbstractString)
    entries = split(value, ','; keepempty=true)
    any(isempty, entries) &&
        throw(ArgumentError("$name must be a comma-separated list without empty entries"))
    try
        return parser.(entries)
    catch
        throw(ArgumentError("$name contains an invalid value: $value"))
    end
end

function _value(arguments, index, option)
    index < length(arguments) ||
        throw(ArgumentError("$option requires a value"))
    return arguments[index + 1]
end

function _parse_integer(value, option)
    try
        return parse(Int, value)
    catch
        throw(ArgumentError("$option must be an integer"))
    end
end

function _parse_float(value, option)
    try
        return parse(Float64, value)
    catch
        throw(ArgumentError("$option must be a number"))
    end
end

function _parse_arguments(arguments)
    settings = Dict{Symbol,Any}(
        :distances => [3, 5, 7],
        :error_rates => nothing,
        :p_min => 0.0,
        :p_max => 0.1,
        :p_step => 0.01,
        :logical_state => :native,
        :construction => :bp,
        :boundary_orientation => :x_ns,
        :clock => :gate_layer,
        :shots => 10_000,
        :batch_size => 10_000,
        :seed => 1234,
        :output_dir => joinpath("results", "rotated-planar"),
        :basename => "rotated_planar_logical_failure",
    )
    explicit_error_rates = false
    explicit_grid = false
    index = 1
    while index <= length(arguments)
        option = arguments[index]
        option == "--help" && return nothing
        option in _VALUE_OPTIONS ||
            throw(ArgumentError("unknown option: $option"))
        value = _value(arguments, index, option)
        if option == "--distances"
            settings[:distances] = _comma_separated(value, x -> parse(Int, x), option)
        elseif option == "--error-rates"
            settings[:error_rates] = _comma_separated(value, x -> parse(Float64, x), option)
            explicit_error_rates = true
        elseif option == "--p-min"
            settings[:p_min] = _parse_float(value, option)
            explicit_grid = true
        elseif option == "--p-max"
            settings[:p_max] = _parse_float(value, option)
            explicit_grid = true
        elseif option == "--p-step"
            settings[:p_step] = _parse_float(value, option)
            explicit_grid = true
        elseif option == "--state"
            settings[:logical_state] = Symbol(value)
        elseif option == "--construction"
            settings[:construction] = Symbol(value)
        elseif option == "--boundary-orientation"
            settings[:boundary_orientation] = Symbol(value)
        elseif option == "--clock"
            settings[:clock] = Symbol(value)
        elseif option == "--shots"
            settings[:shots] = _parse_integer(value, option)
        elseif option == "--batch-size"
            settings[:batch_size] = _parse_integer(value, option)
        elseif option == "--seed"
            settings[:seed] = lowercase(value) in ("nothing", "none") ?
                nothing : _parse_integer(value, option)
        elseif option == "--output-dir"
            settings[:output_dir] = value
        elseif option == "--basename"
            settings[:basename] = value
        end
        index += 2
    end

    explicit_error_rates && explicit_grid && throw(ArgumentError(
        "--error-rates cannot be combined with --p-min, --p-max, or --p-step"))
    if !explicit_error_rates
        p_min = settings[:p_min]
        p_max = settings[:p_max]
        p_step = settings[:p_step]
        isfinite(p_min) && isfinite(p_max) && 0 <= p_min <= p_max < 0.5 ||
            throw(ArgumentError(
                "--p-min and --p-max must be finite and satisfy " *
                "0 <= --p-min <= --p-max < 0.5"))
        isfinite(p_step) && p_step > 0 ||
            throw(ArgumentError("--p-step must be finite and positive"))
        grid, grid_size = try
            grid = p_min:p_step:p_max
            (grid, length(grid))
        catch error
            error isa InexactError || error isa OverflowError || rethrow()
            throw(ArgumentError("the requested error grid is too large"))
        end
        grid_size <= _MAX_ERROR_GRID_POINTS || throw(ArgumentError(
            "the requested error grid has $grid_size points; " *
            "use at most $_MAX_ERROR_GRID_POINTS"))
        settings[:error_rates] = collect(grid)
    end
    return settings
end

function main(args=ARGS; io::IO=stdout, error_io::IO=stderr)
    settings = try
        _parse_arguments(args)
    catch error
        error isa ArgumentError || rethrow()
        println(error_io, "error: ", error.msg)
        println(error_io, USAGE)
        return 1
    end
    if settings === nothing
        println(io, USAGE)
        return 0
    end

    try
        rng = settings[:seed] === nothing ? Random.default_rng() :
            MersenneTwister(settings[:seed])
        scan = scan_logical_failure(
            rng;
            distances=settings[:distances],
            error_rates=settings[:error_rates],
            logical_state=settings[:logical_state],
            construction=settings[:construction],
            boundary_orientation=settings[:boundary_orientation],
            clock=settings[:clock],
            shots=settings[:shots],
            batch_size=settings[:batch_size],
            seed=settings[:seed],
            progress_io=io,
        )
        paths = save_logical_failure_scan(
            scan, settings[:output_dir]; basename=settings[:basename])
        for path in paths
            println(io, "wrote: $path")
        end
    catch error
        error isa ArgumentError || error isa ErrorException || rethrow()
        println(error_io, "error: ", sprint(showerror, error))
        return 1
    end
    return 0
end

end


if abspath(PROGRAM_FILE) == @__FILE__
    exit(RotatedPlanarScanExample.main())
end
