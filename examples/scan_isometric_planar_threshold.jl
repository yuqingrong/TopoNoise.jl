module ScanIsometricPlanarThresholdExample

using TopoNoise
import CairoMakie
import Random

const USAGE = """
Usage: julia --project=. examples/scan_isometric_planar_threshold.jl [options]

Options:
  --p-min P          first virtual-bond X error rate (default: 0.0)
  --p-max P          last virtual-bond X error rate (default: 0.16)
  --p-step DP        grid spacing that reaches p-max (default: 0.005)
  --distances LIST   comma-separated direct distances (default: 3,5,7,9)
  --shots N          samples per (d,p) point (default: 100000)
  --batches N        equal bootstrap batches (default: 100)
  --bootstrap N      bootstrap replicates (default: 2000)
  --confidence X     confidence level (default: 0.95)
  --seed N           RNG seed (default: 1234)
  --output-dir PATH  output directory (default: results/isometric-planar-code)
"""

const _KNOWN_OPTIONS = Set((
    "--p-min", "--p-max", "--p-step", "--distances", "--shots", "--batches",
    "--bootstrap", "--confidence", "--seed", "--output-dir"))

function _option_dictionary(args)
    iseven(length(args)) || throw(ArgumentError("every option must be followed by a value"))
    options = Dict{String,String}()
    for index in 1:2:length(args)
        option, value = args[index], args[index + 1]
        option in _KNOWN_OPTIONS || throw(ArgumentError("unknown option: $option"))
        haskey(options, option) && throw(ArgumentError("option appears more than once: $option"))
        options[option] = value
    end
    return options
end

function _parse(::Type{T}, value, name) where {T<:Real}
    return try
        parse(T, value)
    catch
        throw(ArgumentError("$name must be a valid $(T) value"))
    end
end

function _parse_distances(value)
    distances = try
        parse.(Int, split(value, ','))
    catch
        throw(ArgumentError("--distances must be comma-separated integers"))
    end
    length(distances) >= 2 || throw(ArgumentError(
        "--distances must contain at least two distances"))
    all(>=(2), distances) || throw(ArgumentError(
        "--distances entries must be at least 2"))
    all(distances[index] < distances[index + 1]
        for index in 1:(length(distances) - 1)) || throw(ArgumentError(
        "--distances must be strictly increasing"))
    return distances
end

function _rates(options)
    p_min = _parse(Float64, get(options, "--p-min", "0.0"), "--p-min")
    p_max = _parse(Float64, get(options, "--p-max", "0.16"), "--p-max")
    p_step = _parse(Float64, get(options, "--p-step", "0.005"), "--p-step")
    all(isfinite, (p_min, p_max, p_step)) && 0 <= p_min < p_max < 0.5 ||
        throw(ArgumentError("require finite 0 <= p-min < p-max < 0.5"))
    p_step > 0 || throw(ArgumentError("--p-step must be positive"))
    intervals = (p_max - p_min) / p_step
    count = round(Int, intervals)
    isapprox(intervals, count; atol=1e-10, rtol=1e-10) || throw(ArgumentError(
        "--p-step must exactly reach --p-max"))
    count >= 1 || throw(ArgumentError("--p-step must produce at least two rates"))
    rates = collect(p_min:p_step:p_max)
    length(rates) == count + 1 || throw(ArgumentError(
        "--p-step must produce the expected number of rates"))
    rates[end] = p_max
    return rates
end

function _parse_options(args)
    raw = _option_dictionary(args)
    rates = _rates(raw)
    distances = _parse_distances(get(raw, "--distances", "3,5,7,9"))
    shots = _parse(Int, get(raw, "--shots", "100000"), "--shots")
    batches = _parse(Int, get(raw, "--batches", "100"), "--batches")
    bootstrap = _parse(Int, get(raw, "--bootstrap", "2000"), "--bootstrap")
    seed = _parse(Int, get(raw, "--seed", "1234"), "--seed")
    confidence = _parse(Float64, get(raw, "--confidence", "0.95"), "--confidence")
    shots > 0 && 2 <= batches <= shots && shots % batches == 0 || throw(ArgumentError(
        "require positive shots divisible by batches with 2 <= batches <= shots"))
    bootstrap > 0 || throw(ArgumentError("--bootstrap must be positive"))
    isfinite(confidence) && 0 < confidence < 1 || throw(ArgumentError(
        "--confidence must be strictly between 0 and 1"))
    output_dir = get(raw, "--output-dir", joinpath(
        @__DIR__, "..", "results", "isometric-planar-code"))
    isempty(output_dir) && throw(ArgumentError("--output-dir must not be empty"))
    return (; rates, distances, shots, batches, bootstrap, seed, confidence, output_dir)
end

_csv(value) = ismissing(value) ? "" : value isa Symbol ? string(value) : repr(value)

function _write_scan_csv(path, scan)
    open(path, "w") do io
        println(io, "d,p,shots,logical_failure_count,logical_failure_rate,logical_failure_se,batches,diagnostic_only,seed")
        for point in scan.points
            println(io, join((_csv(point.distance), _csv(point.error_rate),
                _csv(point.shots), _csv(point.logical_failure_count),
                _csv(point.logical_failure_rate), _csv(point.logical_failure_se),
                _csv(length(point.logical_failure_batches)), _csv(point.diagnostic_only),
                _csv(point.seed)), ','))
        end
    end
    return path
end

function _write_batches_csv(path, scan)
    open(path, "w") do io
        println(io, "d,p,batch,batch_shots,failure_count,failure_rate")
        for point in scan.points
            batch_count = length(point.logical_failure_batches)
            batch_count == 0 && continue
            batch_shots = div(point.shots, batch_count)
            for (batch, failure_rate) in enumerate(point.logical_failure_batches)
                failure_count = round(Int, failure_rate * batch_shots)
                println(io, join((
                    _csv(point.distance), _csv(point.error_rate), _csv(batch),
                    _csv(batch_shots), _csv(failure_count), _csv(failure_rate)), ','))
            end
        end
    end
    return path
end

function _write_crossings_csv(path, crossings)
    open(path, "w") do io
        println(io, "small_distance,large_distance,estimate,ci_low,ci_high,confidence,valid_bootstrap_fraction,status")
        for crossing in crossings
            println(io, join((_csv(crossing.small_size), _csv(crossing.large_size),
                _csv(crossing.estimate), _csv(crossing.ci_low),
                _csv(crossing.ci_high), _csv(crossing.confidence),
                _csv(crossing.valid_bootstrap_fraction), _csv(crossing.status)), ','))
        end
    end
    return path
end

function _write_scaling_fit_csv(path, fit, bootstrap)
    open(path, "w") do io
        println(io, join((
            "p_c", "nu", "p_c_se", "nu_se", "valid_bootstrap_fraction",
            "status", "failure_window_low", "failure_window_high", "distances",
            "point_count", "bootstrap"), ','))
        println(io, join((
            _csv(fit.p_c), _csv(fit.nu), _csv(fit.p_c_se), _csv(fit.nu_se),
            _csv(fit.valid_bootstrap_fraction), _csv(fit.status),
            _csv(fit.failure_window[1]), _csv(fit.failure_window[2]),
            join(fit.distances, ';'), _csv(fit.point_count), _csv(bootstrap)), ','))
    end
    return path
end

function _write_scaling_bootstrap_csv(path, diagnostics)
    open(path, "w") do io
        println(io, join((
            "replicate", "status", "input_rms_shift", "initial_p_c", "initial_nu",
            "p_c", "nu", "loss", "converged", "evaluations", "boundary_hit"), ','))
        for sample in diagnostics.bootstrap_samples
            println(io, join((
                _csv(sample.replicate), _csv(sample.status),
                _csv(sample.input_rms_shift), _csv(sample.initial_p_c),
                _csv(sample.initial_nu), _csv(sample.p_c), _csv(sample.nu),
                _csv(sample.loss), _csv(sample.converged),
                _csv(sample.evaluations), _csv(sample.boundary_hit)), ','))
        end
    end
    return path
end

function _write_scaling_loss_surface_csv(path, diagnostics)
    open(path, "w") do io
        println(io, "p_c,nu,loss,delta_loss,log10_one_plus_delta_loss")
        isempty(diagnostics.loss_values) && return path
        minimum_loss = minimum(diagnostics.loss_values)
        for (p_index, p_c) in enumerate(diagnostics.loss_p_c),
                (nu_index, nu) in enumerate(diagnostics.loss_nu)
            loss = diagnostics.loss_values[p_index, nu_index]
            delta = max(0.0, loss - minimum_loss)
            println(io, join((
                _csv(p_c), _csv(nu), _csv(loss), _csv(delta),
                _csv(log10(1 + delta))), ','))
        end
    end
    return path
end

function _write_scaling_optimizer_trace_csv(path, diagnostics)
    open(path, "w") do io
        println(io, join((
            "step", "p_c", "nu", "loss", "initial_p_c", "initial_nu",
            "converged", "evaluations", "boundary_hit"), ','))
        for step in eachindex(diagnostics.optimizer_p_c)
            println(io, join((
                _csv(step), _csv(diagnostics.optimizer_p_c[step]),
                _csv(diagnostics.optimizer_nu[step]),
                _csv(diagnostics.optimizer_loss[step]),
                _csv(diagnostics.nominal_initial_p_c),
                _csv(diagnostics.nominal_initial_nu),
                _csv(diagnostics.nominal_converged),
                _csv(diagnostics.nominal_evaluations),
                _csv(diagnostics.nominal_boundary_hit)), ','))
        end
    end
    return path
end

function _write_scaling_sensitivity_csv(path, diagnostics)
    open(path, "w") do io
        println(io, join((
            "label", "status", "p_c", "nu", "loss", "failure_window_low",
            "failure_window_high", "distances", "boundary_hit"), ','))
        for summary in diagnostics.sensitivities
            println(io, join((
                _csv(summary.label), _csv(summary.status), _csv(summary.p_c),
                _csv(summary.nu), _csv(summary.loss),
                _csv(summary.failure_window[1]), _csv(summary.failure_window[2]),
                join(summary.distances, ';'), _csv(summary.boundary_hit)), ','))
        end
    end
    return path
end

function run(options; io::IO=stdout)
    rng = Random.MersenneTwister(options.seed)
    scan = scan_isometric_planar_capacity(
        rng, options.distances, options.rates; shots=options.shots,
        batches=options.batches, seed=options.seed, progress_io=io)
    crossings = estimate_isometric_planar_crossings(
        rng, scan; bootstrap=options.bootstrap, confidence=options.confidence)
    scaling_diagnostics = diagnose_isometric_planar_scaling(
        rng, scan; bootstrap=options.bootstrap)
    scaling_fit = scaling_diagnostics.fit
    figure = plot_isometric_planar_capacity(
        scan, crossings; scaling_fit=scaling_fit)
    diagnostic_figure = plot_isometric_planar_scaling_diagnostics(
        scan, scaling_diagnostics)
    mkpath(options.output_dir)
    scan_image_paths = (
        joinpath(options.output_dir, "isometric_planar_capacity_scan.svg"),
        joinpath(options.output_dir, "isometric_planar_capacity_scan.pdf"),
        joinpath(options.output_dir, "isometric_planar_capacity_scan.png"))
    diagnostic_image_paths = (
        joinpath(options.output_dir,
            "isometric_planar_capacity_scaling_diagnostics.svg"),
        joinpath(options.output_dir,
            "isometric_planar_capacity_scaling_diagnostics.pdf"),
        joinpath(options.output_dir,
            "isometric_planar_capacity_scaling_diagnostics.png"))
    paths = (
        _write_scan_csv(joinpath(options.output_dir, "isometric_planar_capacity_scan.csv"), scan),
        _write_batches_csv(joinpath(
            options.output_dir, "isometric_planar_capacity_batches.csv"), scan),
        _write_crossings_csv(joinpath(options.output_dir, "isometric_planar_capacity_crossings.csv"), crossings),
        _write_scaling_fit_csv(joinpath(
            options.output_dir, "isometric_planar_capacity_scaling_fit.csv"),
            scaling_fit, options.bootstrap),
        _write_scaling_bootstrap_csv(joinpath(
            options.output_dir, "isometric_planar_capacity_scaling_bootstrap.csv"),
            scaling_diagnostics),
        _write_scaling_loss_surface_csv(joinpath(
            options.output_dir, "isometric_planar_capacity_scaling_loss_surface.csv"),
            scaling_diagnostics),
        _write_scaling_optimizer_trace_csv(joinpath(
            options.output_dir, "isometric_planar_capacity_scaling_optimizer_trace.csv"),
            scaling_diagnostics),
        _write_scaling_sensitivity_csv(joinpath(
            options.output_dir, "isometric_planar_capacity_scaling_sensitivity.csv"),
            scaling_diagnostics),
        scan_image_paths..., diagnostic_image_paths...)
    for path in scan_image_paths
        CairoMakie.save(path, figure)
    end
    for path in diagnostic_image_paths
        CairoMakie.save(path, diagnostic_figure)
    end
    for path in paths
        println(io, "wrote: $(abspath(path))")
    end
    return paths
end

function main(args=ARGS; io::IO=stdout, error_io::IO=stderr)
    try
        run(_parse_options(args); io=io)
        return 0
    catch error
        error isa ArgumentError || rethrow()
        println(error_io, error)
        println(error_io, USAGE)
        return 2
    end
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(ScanIsometricPlanarThresholdExample.main())
end
