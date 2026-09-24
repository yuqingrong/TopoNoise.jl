module RotatedPlanarConstructionComparison

using Random
using TopoNoise

const USAGE = """
Usage: julia --project=. examples/compare_rotated_planar_constructions.jl [options]

Compare As/Bp encoder constructions under X-only and Z-only circuit noise.
By default As prepares |+_L> from ideal all-plus input and Bp prepares |0_L>
from ideal all-zero input. Their noisy H/CNOT schedules are basis-matched.
An explicit logical state applies to all panels.
`--metric channel_logical` reports channel-matched residual logical parity;
`--metric state_failure` reports failure of the selected prepared state.

Options:
  --distances LIST             comma-separated odd distances (default: 3,5,7)
  --error-rates LIST           comma-separated physical error rates
  --p-min VALUE                error-grid minimum (default: 0)
  --p-max VALUE                error-grid maximum (default: 0.16)
  --p-step VALUE               error-grid step (default: 0.01)
  --shots INTEGER              shots per point (default: 10000)
  --batch-size INTEGER         estimator batch size (default: 10000)
  --seed INTEGER               RNG seed (default: 1234)
  --boundary-orientation NAME  x_ns or x_ew (default: x_ns)
  --logical-state NAME         native, zero, one, plus, or minus (default: native)
  --clock NAME                 gate_layer, plaquette, or post_encoding (default: gate_layer)
  --metric NAME                channel_logical or state_failure (default: channel_logical)
  --bootstrap-replicates N     fitting bootstrap replicates (default: 500)
  --no-fit                     save raw data and figures without fitting
  --output-dir PATH            artifact directory (default: results/rotated-planar-comparison)
  --basename NAME              artifact filename stem (default: rotated_planar_comparison)
  --help                       show this help
"""

const _VALUE_OPTIONS = (
    "--distances",
    "--error-rates",
    "--p-min",
    "--p-max",
    "--p-step",
    "--shots",
    "--batch-size",
    "--seed",
    "--boundary-orientation",
    "--logical-state",
    "--clock",
    "--metric",
    "--bootstrap-replicates",
    "--output-dir",
    "--basename",
)

const _MAX_ERROR_GRID_POINTS = 1_000_000

function _comma_separated(value::AbstractString, parser, name::AbstractString)
    entries = split(value, ','; keepempty=true)
    any(isempty, entries) && throw(ArgumentError(
        "$name must be a comma-separated list without empty entries"))
    try
        return parser.(entries)
    catch
        throw(ArgumentError("$name contains an invalid value: $value"))
    end
end

function _value(arguments, index, option)
    index < length(arguments) || throw(ArgumentError("$option requires a value"))
    value = arguments[index + 1]
    startswith(value, "--") && throw(ArgumentError("$option requires a value"))
    return value
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

function _validate_distances(distances)
    isempty(distances) && throw(ArgumentError("--distances must not be empty"))
    all(distance -> distance >= 3 && isodd(distance), distances) || throw(ArgumentError(
        "--distances must contain odd integers at least 3"))
    all(diff(distances) .> 0) || throw(ArgumentError(
        "--distances must be strictly increasing"))
    return distances
end

function _validate_error_rates(rates)
    isempty(rates) && throw(ArgumentError("--error-rates must not be empty"))
    all(rate -> isfinite(rate) && 0 <= rate < 0.5, rates) || throw(ArgumentError(
        "error rates must be finite and satisfy 0 <= p < 0.5"))
    all(diff(rates) .> 0) || length(rates) == 1 || throw(ArgumentError(
        "error rates must be strictly increasing"))
    return rates
end

function _validate_positive(value::Integer, option::AbstractString)
    value > 0 || throw(ArgumentError("$option must be positive"))
    return value
end

function _validate_basename(basename::AbstractString)
    isempty(basename) && throw(ArgumentError("--basename must not be empty"))
    (basename in (".", "..") || occursin('/', basename) || occursin('\\', basename)) &&
        throw(ArgumentError("--basename must be a filename without path components"))
    return String(basename)
end

function _error_grid(p_min::Float64, p_max::Float64, p_step::Float64)
    isfinite(p_min) && isfinite(p_max) && 0 <= p_min <= p_max < 0.5 ||
        throw(ArgumentError(
            "--p-min and --p-max must be finite and satisfy 0 <= --p-min <= --p-max < 0.5"))
    isfinite(p_step) && p_step > 0 ||
        throw(ArgumentError("--p-step must be finite and positive"))
    grid, grid_size = try
        values = p_min:p_step:p_max
        (values, length(values))
    catch error
        error isa InexactError || error isa OverflowError || rethrow()
        throw(ArgumentError("the requested error grid is too large"))
    end
    grid_size <= _MAX_ERROR_GRID_POINTS || throw(ArgumentError(
        "the requested error grid has $grid_size points; use at most $_MAX_ERROR_GRID_POINTS"))
    return _validate_error_rates(collect(grid))
end

function _parse_arguments(arguments)
    settings = Dict{Symbol,Any}(
        :distances => [3, 5, 7],
        :error_rates => nothing,
        :p_min => 0.0,
        :p_max => 0.16,
        :p_step => 0.01,
        :shots => 10_000,
        :batch_size => 10_000,
        :seed => 1234,
        :boundary_orientation => :x_ns,
        :logical_state => :native,
        :clock => :gate_layer,
        :metric => :channel_logical,
        :bootstrap_replicates => 500,
        :no_fit => false,
        :output_dir => joinpath("results", "rotated-planar-comparison"),
        :basename => "rotated_planar_comparison",
    )
    explicit_error_rates = false
    explicit_grid = false
    index = 1
    while index <= length(arguments)
        option = arguments[index]
        option == "--help" && return nothing
        if option == "--no-fit"
            settings[:no_fit] = true
            index += 1
            continue
        end
        option in _VALUE_OPTIONS || throw(ArgumentError("unknown option: $option"))
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
        elseif option == "--shots"
            settings[:shots] = _parse_integer(value, option)
        elseif option == "--batch-size"
            settings[:batch_size] = _parse_integer(value, option)
        elseif option == "--seed"
            settings[:seed] = _parse_integer(value, option)
        elseif option == "--boundary-orientation"
            settings[:boundary_orientation] = Symbol(value)
        elseif option == "--logical-state"
            settings[:logical_state] = Symbol(value)
        elseif option == "--clock"
            settings[:clock] = Symbol(value)
        elseif option == "--metric"
            settings[:metric] = Symbol(value)
        elseif option == "--bootstrap-replicates"
            settings[:bootstrap_replicates] = _parse_integer(value, option)
        elseif option == "--output-dir"
            settings[:output_dir] = value
        elseif option == "--basename"
            settings[:basename] = value
        end
        index += 2
    end

    explicit_error_rates && explicit_grid && throw(ArgumentError(
        "--error-rates cannot be combined with --p-min, --p-max, or --p-step"))
    _validate_distances(settings[:distances])
    settings[:error_rates] = explicit_error_rates ?
        _validate_error_rates(settings[:error_rates]) :
        _error_grid(settings[:p_min], settings[:p_max], settings[:p_step])
    _validate_positive(settings[:shots], "--shots")
    _validate_positive(settings[:batch_size], "--batch-size")
    _validate_positive(settings[:bootstrap_replicates], "--bootstrap-replicates")
    settings[:boundary_orientation] in (:x_ns, :x_ew) || throw(ArgumentError(
        "--boundary-orientation must be x_ns or x_ew"))
    settings[:logical_state] in (:native, :zero, :one, :plus, :minus) || throw(ArgumentError(
        "--logical-state must be native, zero, one, plus, or minus"))
    settings[:clock] in (:gate_layer, :plaquette, :post_encoding) || throw(ArgumentError(
        "--clock must be gate_layer, plaquette, or post_encoding"))
    settings[:metric] in (:channel_logical, :state_failure) || throw(ArgumentError(
        "--metric must be channel_logical or state_failure"))
    settings[:basename] = _validate_basename(settings[:basename])
    return settings
end

function _disabled_fitted_comparison(raw::ConstructionChannelComparison)
    fits = ThresholdFit[]
    for scan in raw.series
        push!(fits, ThresholdFit(
            scan.construction, scan.error_channel, scan.logical_observable,
            :unavailable, "fit disabled by --no-fit", PairCrossing[],
            nothing, nothing, nothing, nothing, nothing, nothing,
            0, 0, UInt64(0), true))
    end
    return FittedConstructionChannelComparison(raw, fits)
end

function _print_paths(io::IO, paths)
    println(io, "wrote: $(paths.raw_csv)")
    println(io, "wrote: $(paths.fits_csv)")
    for path in values(paths.combined)
        println(io, "wrote: $path")
    end
    for panel in values(paths.panels), path in values(panel)
        println(io, "wrote: $path")
    end
end

function _print_fit_summary(io::IO, fitted::FittedConstructionChannelComparison)
    for fit in fitted.fits
        if fit.status === :success
            println(io,
                "construction=$(fit.construction) channel=$(fit.error_channel) " *
                "fit $(fit.status) p_c=$(fit.p_c) nu=$(fit.nu)")
        else
            println(io,
                "construction=$(fit.construction) channel=$(fit.error_channel) " *
                "fit unavailable: $(fit.diagnostic) p_c=unavailable nu=unavailable")
        end
    end
end

function main(args=ARGS; io::IO=stdout, error_io::IO=stderr)::Int
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
        rng = MersenneTwister(settings[:seed])
        raw = run_construction_channel_comparison(
            rng;
            distances=settings[:distances],
            error_rates=settings[:error_rates],
            logical_state=settings[:logical_state],
            boundary_orientation=settings[:boundary_orientation],
            clock=settings[:clock],
            failure_metric=settings[:metric],
            shots=settings[:shots],
            batch_size=settings[:batch_size],
            seed=settings[:seed],
            progress_io=io,
        )
        fitted = settings[:no_fit] ? _disabled_fitted_comparison(raw) :
            fit_construction_channel_comparison(
                rng, raw; bootstrap_replicates=settings[:bootstrap_replicates])
        paths = save_construction_channel_comparison(
            fitted, settings[:output_dir]; basename=settings[:basename])
        _print_paths(io, paths)
        _print_fit_summary(io, fitted)
    catch error
        error isa ArgumentError || error isa ErrorException || rethrow()
        println(error_io, "error: ", sprint(showerror, error))
        return 1
    end
    return 0
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(RotatedPlanarConstructionComparison.main())
end
