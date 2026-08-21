"""Finite-size-scaling fit for decoded planar virtual-bond capacity curves."""
struct IsometricPlanarScalingFit
    p_c::Union{Missing,Float64}
    nu::Union{Missing,Float64}
    p_c_se::Union{Missing,Float64}
    nu_se::Union{Missing,Float64}
    valid_bootstrap_fraction::Float64
    status::Symbol
    failure_window::Tuple{Float64,Float64}
    distances::Vector{Int}
    point_count::Int
    coefficients::Vector{Float64}
end

function _isometric_planar_scaling_unavailable(
        status::Symbol, window::Tuple{Float64,Float64}, distances::Vector{Int},
        point_count::Int)
    return IsometricPlanarScalingFit(
        missing, missing, missing, missing, 0.0, status, window, distances,
        point_count, Float64[])
end

function _isometric_planar_scaling_selected_points(
        scan::IsometricPlanarCapacityScan, window::Tuple{Float64,Float64})
    lower, upper = window
    return [point for point in scan.points if !point.diagnostic_only &&
            point.distance > 2 && lower <= point.logical_failure_rate <= upper]
end

struct IsometricPlanarScalingBootstrapSample
    replicate::Int
    status::Symbol
    input_rms_shift::Float64
    initial_p_c::Union{Missing,Float64}
    initial_nu::Union{Missing,Float64}
    p_c::Union{Missing,Float64}
    nu::Union{Missing,Float64}
    loss::Union{Missing,Float64}
    converged::Bool
    evaluations::Int
    boundary_hit::Bool
end

struct IsometricPlanarScalingBatchDiagnostic
    distance::Int
    error_rate::Float64
    batch_count::Int
    batch_size::Int
    batch_mean::Float64
    batch_std::Float64
    expected_batch_std::Float64
    std_ratio::Float64
    unique_batch_rates::Int
end

struct IsometricPlanarScalingSensitivity
    label::String
    status::Symbol
    p_c::Union{Missing,Float64}
    nu::Union{Missing,Float64}
    loss::Union{Missing,Float64}
    failure_window::Tuple{Float64,Float64}
    distances::Vector{Int}
    boundary_hit::Bool
end

"""Detailed evidence for auditing a finite-size-scaling fit."""
struct IsometricPlanarScalingDiagnostics
    fit::IsometricPlanarScalingFit
    bootstrap_samples::Vector{IsometricPlanarScalingBootstrapSample}
    loss_p_c::Vector{Float64}
    loss_nu::Vector{Float64}
    loss_values::Matrix{Float64}
    optimizer_p_c::Vector{Float64}
    optimizer_nu::Vector{Float64}
    optimizer_loss::Vector{Float64}
    nominal_initial_p_c::Union{Missing,Float64}
    nominal_initial_nu::Union{Missing,Float64}
    nominal_converged::Bool
    nominal_evaluations::Int
    nominal_boundary_hit::Bool
    batch_diagnostics::Vector{IsometricPlanarScalingBatchDiagnostic}
    sensitivities::Vector{IsometricPlanarScalingSensitivity}
end

struct _IsometricPlanarScalingSimplexer <: Optim.Simplexer
    minimum_rate::Float64
    maximum_rate::Float64
end

function Optim.simplexer(simplexer::_IsometricPlanarScalingSimplexer, initial)
    simplex = [copy(initial) for _ in 1:3]
    p_step = max((simplexer.maximum_rate - simplexer.minimum_rate) / 50, 1e-6)
    p_direction = initial[1] + p_step <= simplexer.maximum_rate ? 1.0 : -1.0
    simplex[2][1] += p_direction * p_step
    minimum_log_nu, maximum_log_nu = log(0.25), log(6.0)
    log_nu_step = (maximum_log_nu - minimum_log_nu) / 50
    nu_direction = initial[2] + log_nu_step <= maximum_log_nu ? 1.0 : -1.0
    simplex[3][2] += nu_direction * log_nu_step
    return simplex
end

function _isometric_planar_scaling_problem(
        points::Vector{IsometricPlanarCapacityPoint}, values::Vector{Float64})
    length(points) == length(values) || throw(DimensionMismatch(
        "scaling-fit values must match selected scan points"))
    length(points) >= 6 || return nothing
    all(isfinite, values) || return nothing
    rates = Float64[point.error_rate for point in points]
    minimum_rate, maximum_rate = extrema(rates)
    minimum_rate < maximum_rate || return nothing
    weights = Float64[1 / max(point.logical_failure_se, 1 / sqrt(point.shots))^2
                      for point in points]
    square_roots = sqrt.(weights)

    function objective(parameters)
        p_c, log_nu = parameters
        nu = exp(log_nu)
        minimum_rate <= p_c <= maximum_rate && 0.25 <= nu <= 6.0 || return Inf
        x = Float64[(point.error_rate - p_c) * point.distance^(1 / nu)
                    for point in points]
        design = hcat(ones(length(x)), x, x .^ 2, x .^ 3)
        rank(design) == 4 || return Inf
        coefficients = (square_roots .* design) \ (square_roots .* values)
        residuals = values .- design * coefficients
        return sum(weights .* residuals .^ 2)
    end

    return (; objective, minimum_rate, maximum_rate, square_roots)
end

function _isometric_planar_scaling_boundary_hit(
        p_c::Float64, nu::Float64, minimum_rate::Float64, maximum_rate::Float64)
    p_tolerance = max(1e-10, 1e-6 * (maximum_rate - minimum_rate))
    nu_tolerance = 1e-6 * (6.0 - 0.25)
    return abs(p_c - minimum_rate) <= p_tolerance ||
           abs(p_c - maximum_rate) <= p_tolerance ||
           abs(nu - 0.25) <= nu_tolerance || abs(nu - 6.0) <= nu_tolerance
end

function _isometric_planar_scaling_profile(
        points::Vector{IsometricPlanarCapacityPoint}, values::Vector{Float64};
        record_trace::Bool=false)
    problem = _isometric_planar_scaling_problem(points, values)
    problem === nothing && return nothing
    objective = problem.objective

    p_candidates = range(problem.minimum_rate, problem.maximum_rate; length=13)
    nu_candidates = range(0.5, 3.0; length=11)
    starts = [(p_c, log(nu)) for p_c in p_candidates, nu in nu_candidates]
    losses = [objective(start) for start in starts]
    initial = starts[argmin(losses)]
    initial_loss = minimum(losses)
    evaluations = Ref(0)
    best_loss = Ref(initial_loss)
    trace_p_c = Float64[initial[1]]
    trace_nu = Float64[exp(initial[2])]
    trace_loss = Float64[initial_loss]
    function optimized_objective(parameters)
        evaluations[] += 1
        loss = objective(parameters)
        if record_trace && isfinite(loss) && loss < best_loss[]
            best_loss[] = loss
            push!(trace_p_c, parameters[1])
            push!(trace_nu, exp(parameters[2]))
            push!(trace_loss, loss)
        end
        return loss
    end
    result = Optim.optimize(
        optimized_objective, collect(initial), Optim.NelderMead(
            initial_simplex=_IsometricPlanarScalingSimplexer(
                problem.minimum_rate, problem.maximum_rate)),
        Optim.Options(iterations=2_000, show_trace=false, store_trace=false))
    parameters = Optim.minimizer(result)
    loss = objective(parameters)
    isfinite(loss) || return nothing
    p_c, log_nu = parameters
    nu = exp(log_nu)
    x = Float64[(point.error_rate - p_c) * point.distance^(1 / nu)
                for point in points]
    design = hcat(ones(length(x)), x, x .^ 2, x .^ 3)
    coefficients = (problem.square_roots .* design) \
                   (problem.square_roots .* values)
    if record_trace && (trace_p_c[end] != p_c || trace_nu[end] != nu)
        push!(trace_p_c, p_c)
        push!(trace_nu, nu)
        push!(trace_loss, loss)
    end
    boundary_hit = _isometric_planar_scaling_boundary_hit(
        p_c, nu, problem.minimum_rate, problem.maximum_rate)
    return (; p_c, nu, coefficients, loss,
            initial_p_c=initial[1], initial_nu=exp(initial[2]),
            converged=Optim.converged(result), evaluations=evaluations[],
            boundary_hit, trace_p_c, trace_nu, trace_loss)
end

function _isometric_planar_scaling_window(failure_window)
    window = (Float64(failure_window[1]), Float64(failure_window[2]))
    all(isfinite, window) && 0 <= window[1] < window[2] <= 1 || throw(ArgumentError(
        "failure_window must satisfy 0 <= lower < upper <= 1"))
    return window
end

function _isometric_planar_scaling_fit_from_analysis(
        nominal, p_c_samples::Vector{Float64}, nu_samples::Vector{Float64},
        bootstrap::Int, window::Tuple{Float64,Float64}, distances::Vector{Int},
        point_count::Int)
    valid_fraction = length(p_c_samples) / bootstrap
    valid_fraction >= 0.8 || return IsometricPlanarScalingFit(
        nominal.p_c, nominal.nu, missing, missing, valid_fraction, :unstable,
        window, distances, point_count, nominal.coefficients)
    p_c_se = length(p_c_samples) >= 2 ? Statistics.std(p_c_samples) : 0.0
    nu_se = length(nu_samples) >= 2 ? Statistics.std(nu_samples) : 0.0
    return IsometricPlanarScalingFit(
        nominal.p_c, nominal.nu, p_c_se, nu_se, valid_fraction, :ok,
        window, distances, point_count, nominal.coefficients)
end

function _run_isometric_planar_scaling(
        rng::Random.AbstractRNG, scan::IsometricPlanarCapacityScan;
        bootstrap::Integer, failure_window, collect_diagnostics::Bool)
    bootstrap > 0 || throw(ArgumentError("bootstrap must be positive, got $bootstrap"))
    window = _isometric_planar_scaling_window(failure_window)
    points = _isometric_planar_scaling_selected_points(scan, window)
    distances = sort(unique(point.distance for point in points))
    if length(distances) < 3
        fit = _isometric_planar_scaling_unavailable(
            :insufficient_sizes, window, distances, length(points))
        return (; fit, points, nominal=nothing,
                bootstrap_samples=IsometricPlanarScalingBootstrapSample[])
    end
    if !all(length(point.logical_failure_batches) >= 2 for point in points)
        fit = _isometric_planar_scaling_unavailable(
            :insufficient_batches, window, distances, length(points))
        return (; fit, points, nominal=nothing,
                bootstrap_samples=IsometricPlanarScalingBootstrapSample[])
    end

    nominal_values = Float64[point.logical_failure_rate for point in points]
    nominal = _isometric_planar_scaling_profile(
        points, nominal_values; record_trace=collect_diagnostics)
    if nominal === nothing
        fit = _isometric_planar_scaling_unavailable(
            :fit_failed, window, distances, length(points))
        return (; fit, points, nominal=nothing,
                bootstrap_samples=IsometricPlanarScalingBootstrapSample[])
    end

    p_c_samples = Float64[]
    nu_samples = Float64[]
    records = IsometricPlanarScalingBootstrapSample[]
    for replicate in 1:Int(bootstrap)
        values = Float64[_bootstrap_batch_mean(rng, point.logical_failure_batches)
                         for point in points]
        differences = values .- nominal_values
        rms_shift = all(isfinite, differences) ?
            sqrt(sum(abs2, differences) / length(differences)) : NaN
        sample = _isometric_planar_scaling_profile(points, values)
        if sample === nothing
            collect_diagnostics && push!(records,
                IsometricPlanarScalingBootstrapSample(
                    replicate, :fit_failed, rms_shift, missing, missing,
                    missing, missing, missing, false, 0, false))
            continue
        end
        push!(p_c_samples, sample.p_c)
        push!(nu_samples, sample.nu)
        collect_diagnostics && push!(records,
            IsometricPlanarScalingBootstrapSample(
                replicate, :ok, rms_shift, sample.initial_p_c, sample.initial_nu,
                sample.p_c, sample.nu, sample.loss, sample.converged,
                sample.evaluations, sample.boundary_hit))
    end
    fit = _isometric_planar_scaling_fit_from_analysis(
        nominal, p_c_samples, nu_samples, Int(bootstrap), window,
        distances, length(points))
    return (; fit, points, nominal, bootstrap_samples=records)
end

"""
    fit_isometric_planar_scaling(rng, scan; bootstrap=2_000,
                                 failure_window=(0.05, 0.45))

Fit the finite-size-scaling form ``P_fail = F((p-p_c)d^(1/nu))`` to the
non-diagnostic points in the requested logical-failure transition window. The
master curve ``F`` is a cubic weighted least-squares polynomial. One-sigma
errors of ``p_c`` and ``nu`` are batch-bootstrap standard deviations.
"""
function fit_isometric_planar_scaling(
        rng::Random.AbstractRNG, scan::IsometricPlanarCapacityScan;
        bootstrap::Integer=2_000,
        failure_window::Tuple{<:Real,<:Real}=(0.05, 0.45))
    return _run_isometric_planar_scaling(
        rng, scan; bootstrap, failure_window,
        collect_diagnostics=false).fit
end

function _isometric_planar_batch_diagnostics(
        points::Vector{IsometricPlanarCapacityPoint})
    summaries = IsometricPlanarScalingBatchDiagnostic[]
    for point in points
        batches = point.logical_failure_batches
        isempty(batches) && continue
        batch_count = length(batches)
        batch_size = div(point.shots, batch_count)
        batch_mean = Statistics.mean(batches)
        batch_std = batch_count >= 2 ? Statistics.std(batches) : 0.0
        expected = sqrt(max(0.0, point.logical_failure_rate *
            (1 - point.logical_failure_rate) / batch_size))
        ratio = expected == 0 ? (batch_std == 0 ? 0.0 : Inf) : batch_std / expected
        push!(summaries, IsometricPlanarScalingBatchDiagnostic(
            point.distance, point.error_rate, batch_count, batch_size,
            batch_mean, batch_std, expected, ratio, length(unique(batches))))
    end
    return summaries
end

function _isometric_planar_sensitivity(
        label::String, points::Vector{IsometricPlanarCapacityPoint},
        window::Tuple{Float64,Float64})
    distances = sort(unique(point.distance for point in points))
    if length(distances) < 3
        return IsometricPlanarScalingSensitivity(
            label, :insufficient_sizes, missing, missing, missing,
            window, distances, false)
    end
    profile = _isometric_planar_scaling_profile(
        points, Float64[point.logical_failure_rate for point in points])
    profile === nothing && return IsometricPlanarScalingSensitivity(
        label, :fit_failed, missing, missing, missing, window, distances, false)
    return IsometricPlanarScalingSensitivity(
        label, :ok, profile.p_c, profile.nu, profile.loss, window,
        distances, profile.boundary_hit)
end

function _isometric_planar_sensitivities(
        scan::IsometricPlanarCapacityScan,
        baseline_points::Vector{IsometricPlanarCapacityPoint},
        window::Tuple{Float64,Float64})
    results = IsometricPlanarScalingSensitivity[
        _isometric_planar_sensitivity("baseline", baseline_points, window)]
    distances = sort(unique(point.distance for point in baseline_points))
    for distance in distances
        selected = [point for point in baseline_points if point.distance != distance]
        push!(results, _isometric_planar_sensitivity(
            "drop d=$distance", selected, window))
    end
    for alternative in ((0.05, 0.35), (0.10, 0.45))
        points = _isometric_planar_scaling_selected_points(scan, alternative)
        label = alternative == (0.05, 0.35) ?
            "window 0.05-0.35" : "window 0.10-0.45"
        push!(results, _isometric_planar_sensitivity(label, points, alternative))
    end
    return results
end

function _isometric_planar_loss_grid(
        points::Vector{IsometricPlanarCapacityPoint}, values::Vector{Float64},
        grid_size::Tuple{Int,Int})
    problem = _isometric_planar_scaling_problem(points, values)
    problem === nothing && return (Float64[], Float64[], zeros(0, 0))
    p_c_values = collect(range(
        problem.minimum_rate, problem.maximum_rate; length=grid_size[1]))
    nu_values = collect(range(0.25, 6.0; length=grid_size[2]))
    losses = Matrix{Float64}(undef, length(p_c_values), length(nu_values))
    for (p_index, p_c) in enumerate(p_c_values),
            (nu_index, nu) in enumerate(nu_values)
        losses[p_index, nu_index] = problem.objective((p_c, log(nu)))
    end
    return p_c_values, nu_values, losses
end

"""
    diagnose_isometric_planar_scaling(rng, scan; bootstrap=2_000,
        failure_window=(0.05, 0.45), loss_grid_size=(121, 121), sensitivity=true)

Fit the finite-size collapse and retain the batch, bootstrap, optimizer, loss
surface, and sensitivity evidence needed to audit the reported uncertainties.
"""
function diagnose_isometric_planar_scaling(
        rng::Random.AbstractRNG, scan::IsometricPlanarCapacityScan;
        bootstrap::Integer=2_000,
        failure_window::Tuple{<:Real,<:Real}=(0.05, 0.45),
        loss_grid_size::Tuple{<:Integer,<:Integer}=(121, 121),
        sensitivity::Bool=true)
    all(size -> size >= 2, loss_grid_size) || throw(ArgumentError(
        "loss_grid_size entries must both be at least 2"))
    grid_size = (Int(loss_grid_size[1]), Int(loss_grid_size[2]))
    analysis = _run_isometric_planar_scaling(
        rng, scan; bootstrap, failure_window, collect_diagnostics=true)
    batch_diagnostics = _isometric_planar_batch_diagnostics(analysis.points)
    if analysis.nominal === nothing
        return IsometricPlanarScalingDiagnostics(
            analysis.fit, analysis.bootstrap_samples, Float64[], Float64[],
            zeros(0, 0), Float64[], Float64[], Float64[], missing, missing,
            false, 0, false, batch_diagnostics,
            IsometricPlanarScalingSensitivity[])
    end
    nominal_values = Float64[point.logical_failure_rate for point in analysis.points]
    loss_p_c, loss_nu, loss_values = _isometric_planar_loss_grid(
        analysis.points, nominal_values, grid_size)
    sensitivities = sensitivity ? _isometric_planar_sensitivities(
        scan, analysis.points, analysis.fit.failure_window) :
        IsometricPlanarScalingSensitivity[]
    nominal = analysis.nominal
    return IsometricPlanarScalingDiagnostics(
        analysis.fit, analysis.bootstrap_samples, loss_p_c, loss_nu,
        loss_values, nominal.trace_p_c, nominal.trace_nu, nominal.trace_loss,
        nominal.initial_p_c, nominal.initial_nu, nominal.converged,
        nominal.evaluations, nominal.boundary_hit, batch_diagnostics,
        sensitivities)
end
