function _isometric_planar_crossing_summary(crossing::CriticalCrossing)
    sizes = "d=$(crossing.small_size)/$(crossing.large_size)"
    if crossing.status == :ok && !ismissing(crossing.estimate) &&
       !ismissing(crossing.ci_low) && !ismissing(crossing.ci_high)
        level = round(Int, 100 * crossing.confidence)
        return "$sizes: p_c=$(round(crossing.estimate; digits=4)) " *
               "[$(round(crossing.ci_low; digits=4)), " *
               "$(round(crossing.ci_high; digits=4))] ($level%)"
    elseif crossing.status == :unbracketed
        return "$sizes: unbracketed; extend p range"
    elseif crossing.status == :unstable
        return "$sizes: unstable; no unique crossing"
    end
    return "$sizes: $(crossing.status)"
end

function _isometric_planar_fit_value(value::Float64, error::Float64, digits::Int)
    uncertainty = round(Int, error * 10^digits)
    return "$(round(value; digits=digits))($uncertainty)"
end

function _isometric_planar_scaling_summary(fit::IsometricPlanarScalingFit)
    fit.status == :ok || return "scaling fit unavailable ($(fit.status))"
    return "p_c=$(_isometric_planar_fit_value(fit.p_c, fit.p_c_se, 4))\n" *
           "ν=$(_isometric_planar_fit_value(fit.nu, fit.nu_se, 2))"
end

function _plot_isometric_planar_scaling_inset!(
        curve_axis, scan::IsometricPlanarCapacityScan,
        inset_position, fit::IsometricPlanarScalingFit, colors)
    fit.status == :ok || return CairoMakie.text!(
        curve_axis, 0.97, 0.05;
        text=_isometric_planar_scaling_summary(fit), space=:relative,
        align=(:right, :bottom), fontsize=13, color=:gray30)
    inset_axis = CairoMakie.Axis(
        inset_position;
        width=CairoMakie.Relative(0.34), height=CairoMakie.Relative(0.36),
        halign=0.10, valign=0.90,
        xlabel="(p-p_c)d^(1/ν)", ylabel="logical failure",
        xlabelsize=11, ylabelsize=11, xticklabelsize=10, yticklabelsize=10,
        backgroundcolor=:white, xgridvisible=false, ygridvisible=false)
    CairoMakie.translate!(inset_axis.blockscene, 0, 0, 150)
    lower, upper = fit.failure_window
    xs = Float64[]
    for (index, distance) in enumerate(scan.distances)
        points = [_isometric_planar_capacity_point(scan, distance, rate)
                  for rate in scan.error_rates]
        selected = [point for point in points if !point.diagnostic_only &&
                    lower <= point.logical_failure_rate <= upper]
        isempty(selected) && continue
        x = [(point.error_rate - fit.p_c) * distance^(1 / fit.nu)
             for point in selected]
        y = [point.logical_failure_rate for point in selected]
        errors = [point.logical_failure_se for point in selected]
        color = colors[mod1(index, length(colors))]
        CairoMakie.errorbars!(inset_axis, x, y, errors;
            color=(color, 0.75), whiskerwidth=5)
        CairoMakie.scatter!(inset_axis, x, y;
            color=color, markersize=8)
        append!(xs, x)
    end
    if !isempty(xs)
        xmin, xmax = extrema(xs)
        padding = max(0.05, 0.08 * (xmax - xmin))
        CairoMakie.xlims!(inset_axis, xmin - padding, xmax + padding)
        CairoMakie.text!(inset_axis, xmax + padding, 0.035;
            text=_isometric_planar_scaling_summary(fit),
            align=(:right, :bottom), fontsize=12, color=:gray15)
    end
    CairoMakie.ylims!(inset_axis, 0.0, 0.5)
    return inset_axis
end

"""
    plot_isometric_planar_capacity(scan, crossings; scaling_fit=nothing)

Render direct-distance planar virtual-bond capacity curves. Distance two is
shown only as a diagnostic because its one-fault decoding tie is unavoidable.
"""
function plot_isometric_planar_capacity(
        scan::IsometricPlanarCapacityScan,
        crossings::AbstractVector{<:CriticalCrossing}=CriticalCrossing[];
        scaling_fit::Union{Nothing,IsometricPlanarScalingFit}=nothing)
    figure = CairoMakie.Figure(size=(950, 600), backgroundcolor=:white)
    curve_axis = CairoMakie.Axis(
        figure[1, 1];
        title="Isometric planar-code logical failure",
        subtitle="internal virtual-bond X noise; plaquette plus local N/S checks",
        xlabel="virtual-bond error rate p", ylabel="logical failure rate",
        yticks=0.0:0.1:0.5,
        xgridvisible=false, ygridvisible=false)
    colors = ("#0072B2", "#D55E00", "#009E73", "#CC79A7",
              "#E69F00", "#56B4E9", "#000000")

    for (index, distance) in enumerate(scan.distances)
        points = [_isometric_planar_capacity_point(scan, distance, rate)
                  for rate in scan.error_rates]
        color = colors[mod1(index, length(colors))]
        failures = [point.logical_failure_rate for point in points]
        errors = [point.logical_failure_se for point in points]
        label = distance == 2 ? "d=2 (diagnostic)" : "d=$distance"
        CairoMakie.lines!(curve_axis, scan.error_rates, failures;
            color=color, linewidth=2, label=label)
        CairoMakie.scatter!(curve_axis, scan.error_rates, failures; color=color)
        CairoMakie.errorbars!(curve_axis, scan.error_rates, failures, errors;
            color=color, whiskerwidth=7)
    end
    if scaling_fit !== nothing && scaling_fit.status == :ok &&
       !ismissing(scaling_fit.p_c)
        line_y = [-0.01, 0.51]
        CairoMakie.lines!(curve_axis, fill(scaling_fit.p_c, length(line_y)), line_y;
            color=(:gray40, 0.75), linestyle=:dash, linewidth=1.5)
    end
    CairoMakie.ylims!(curve_axis, -0.01, 0.51)
    CairoMakie.axislegend(
        curve_axis;
        position=scaling_fit === nothing ? :lt : :rb,
        framevisible=false)
    scaling_fit === nothing || _plot_isometric_planar_scaling_inset!(
        curve_axis, scan, figure[1, 1], scaling_fit, colors)
    return figure
end

function _isometric_planar_successful_bootstrap(
        diagnostics::IsometricPlanarScalingDiagnostics)
    return [sample for sample in diagnostics.bootstrap_samples
            if sample.status == :ok && !ismissing(sample.p_c) &&
               !ismissing(sample.nu)]
end

function _isometric_planar_diagnostic_spread(values)
    return length(values) >= 2 ? Statistics.std(values) : 0.0
end

function _plot_isometric_planar_loss_diagnostics!(position, diagnostics)
    layout = CairoMakie.GridLayout()
    position[] = layout
    axis = CairoMakie.Axis(
        layout[1, 1]; title="Profiled loss surface",
        xlabel="p_c", ylabel="ν")
    minimum_loss = minimum(diagnostics.loss_values)
    transformed = log10.(1 .+ max.(diagnostics.loss_values .- minimum_loss, 0.0))
    heatmap = CairoMakie.heatmap!(
        axis, diagnostics.loss_p_c, diagnostics.loss_nu, transformed)
    CairoMakie.Colorbar(
        layout[1, 2], heatmap; label="log₁₀(1 + Δχ²)", width=12)
    if maximum(transformed) > minimum(transformed)
        CairoMakie.contour!(
            axis, diagnostics.loss_p_c, diagnostics.loss_nu, transformed;
            levels=5, color=(:white, 0.7), linewidth=1)
    end
    grid_p_c = collect(range(
        first(diagnostics.loss_p_c), last(diagnostics.loss_p_c); length=13))
    grid_nu = collect(range(0.5, 3.0; length=11))
    CairoMakie.scatter!(
        axis, repeat(grid_p_c, inner=length(grid_nu)),
        repeat(grid_nu, outer=length(grid_p_c));
        color=(:gray30, 0.35), markersize=3, label="initial grid")
    CairoMakie.lines!(
        axis, diagnostics.optimizer_p_c, diagnostics.optimizer_nu;
        color=:orange, linewidth=2, label="running best")
    CairoMakie.scatter!(
        axis, diagnostics.optimizer_p_c, diagnostics.optimizer_nu;
        color=:orange, markersize=5)
    CairoMakie.scatter!(
        axis, [diagnostics.fit.p_c], [diagnostics.fit.nu];
        color=:red, marker=:star5, markersize=16, label="final minimum")
    CairoMakie.axislegend(axis; position=:rt, framevisible=false)
    return axis
end

function _plot_isometric_planar_bootstrap_diagnostics!(position, diagnostics)
    layout = CairoMakie.GridLayout()
    position[] = layout
    top_axis = CairoMakie.Axis(layout[1, 1])
    joint_axis = CairoMakie.Axis(
        layout[2, 1]; title="Bootstrap parameter distribution",
        xlabel="p_c", ylabel="ν")
    right_axis = CairoMakie.Axis(layout[2, 2])
    CairoMakie.rowsize!(layout, 1, CairoMakie.Relative(0.23))
    CairoMakie.colsize!(layout, 2, CairoMakie.Relative(0.23))
    successful = _isometric_planar_successful_bootstrap(diagnostics)
    if isempty(successful)
        CairoMakie.text!(joint_axis, 0.5, 0.5;
            text="no successful bootstrap refits", space=:relative,
            align=(:center, :center), color=:gray30)
    else
        p_c = Float64[sample.p_c for sample in successful]
        nu = Float64[sample.nu for sample in successful]
        bins = max(1, min(20, round(Int, sqrt(length(successful)))))
        CairoMakie.hist!(top_axis, p_c; bins, color=(:steelblue, 0.75))
        CairoMakie.hist!(right_axis, nu;
            bins, direction=:x, color=(:steelblue, 0.75))
        CairoMakie.scatter!(joint_axis, p_c, nu;
            color=(:steelblue, 0.45), markersize=6)
        CairoMakie.scatter!(joint_axis, [diagnostics.fit.p_c], [diagnostics.fit.nu];
            color=:red, marker=:star5, markersize=15)
        unique_count = length(unique(zip(p_c, nu)))
        p_c_spread = _isometric_planar_diagnostic_spread(p_c)
        nu_spread = _isometric_planar_diagnostic_spread(nu)
        summary = "valid=$(length(successful))/$(length(diagnostics.bootstrap_samples)); " *
                  "unique=$unique_count\nσ(p_c)=$(round(p_c_spread; sigdigits=3)); " *
                  "σ(ν)=$(round(nu_spread; sigdigits=3))"
        CairoMakie.text!(joint_axis, 0.02, 0.98;
            text=summary, space=:relative, align=(:left, :top),
            fontsize=12, color=:gray20)
        if unique_count == 1 || p_c_spread == 0 || nu_spread == 0
            CairoMakie.text!(joint_axis, 0.98, 0.02;
                text="degenerate bootstrap", space=:relative,
                align=(:right, :bottom), fontsize=12, color=:red)
        end
    end
    CairoMakie.hidedecorations!(top_axis)
    CairoMakie.hidespines!(top_axis)
    CairoMakie.hidedecorations!(right_axis)
    CairoMakie.hidespines!(right_axis)
    return joint_axis
end

function _plot_isometric_planar_batch_diagnostics!(position, scan, diagnostics)
    axis = CairoMakie.Axis(
        position; title="Batch variation check", xlabel="error rate p",
        ylabel="empirical SD / binomial SD")
    colors = ("#0072B2", "#D55E00", "#009E73", "#CC79A7",
              "#E69F00", "#56B4E9", "#000000")
    for (index, distance) in enumerate(scan.distances)
        summaries = [summary for summary in diagnostics.batch_diagnostics
                     if summary.distance == distance && isfinite(summary.std_ratio)]
        isempty(summaries) && continue
        color = colors[mod1(index, length(colors))]
        CairoMakie.lines!(axis,
            [summary.error_rate for summary in summaries],
            [summary.std_ratio for summary in summaries];
            color, linewidth=1.5, label="d=$distance")
        CairoMakie.scatter!(axis,
            [summary.error_rate for summary in summaries],
            [summary.std_ratio for summary in summaries]; color, markersize=7)
    end
    CairoMakie.hlines!(axis, [1.0];
        color=:gray30, linestyle=:dash, linewidth=1.5,
        label="binomial expectation")
    CairoMakie.axislegend(axis; position=:rt, framevisible=false)
    return axis
end

function _plot_isometric_planar_sensitivity_diagnostics!(position, diagnostics)
    layout = CairoMakie.GridLayout()
    position[] = layout
    summaries = diagnostics.sensitivities
    labels = [summary.label for summary in summaries]
    positions = collect(length(summaries):-1:1)
    p_axis = CairoMakie.Axis(
        layout[1, 1]; title="Sensitivity: p_c", xlabel="p_c",
        yticks=(positions, labels))
    nu_axis = CairoMakie.Axis(
        layout[1, 2]; title="Sensitivity: ν", xlabel="ν",
        yticks=(positions, fill("", length(positions))))
    successful = [(index, summary) for (index, summary) in enumerate(summaries)
                  if summary.status == :ok && !ismissing(summary.p_c) &&
                     !ismissing(summary.nu)]
    if isempty(successful)
        for axis in (p_axis, nu_axis)
            CairoMakie.text!(axis, 0.5, 0.5;
                text="no successful sensitivity fits", space=:relative,
                align=(:center, :center), color=:gray30)
        end
    else
        y = [positions[index] for (index, _) in successful]
        colors = [summary.boundary_hit ? :red : :steelblue
                  for (_, summary) in successful]
        CairoMakie.scatter!(p_axis,
            Float64[summary.p_c for (_, summary) in successful], y;
            color=colors, markersize=10)
        CairoMakie.scatter!(nu_axis,
            Float64[summary.nu for (_, summary) in successful], y;
            color=colors, markersize=10)
        baseline = first(summaries)
        if baseline.status == :ok
            CairoMakie.vlines!(p_axis, [baseline.p_c];
                color=(:gray30, 0.7), linestyle=:dash)
            CairoMakie.vlines!(nu_axis, [baseline.nu];
                color=(:gray30, 0.7), linestyle=:dash)
            if !ismissing(diagnostics.fit.p_c_se) && diagnostics.fit.p_c_se > 0
                CairoMakie.errorbars!(p_axis, [baseline.p_c], [positions[1]],
                    [diagnostics.fit.p_c_se]; direction=:x,
                    color=:gray20, whiskerwidth=8)
            end
            if !ismissing(diagnostics.fit.nu_se) && diagnostics.fit.nu_se > 0
                CairoMakie.errorbars!(nu_axis, [baseline.nu], [positions[1]],
                    [diagnostics.fit.nu_se]; direction=:x,
                    color=:gray20, whiskerwidth=8)
            end
        end
    end
    return p_axis, nu_axis
end

"""
    plot_isometric_planar_scaling_diagnostics(scan, diagnostics)

Render the profiled loss, bootstrap distribution, batch-variation check, and
fit-window/distance sensitivity evidence in a separate audit figure.
"""
function plot_isometric_planar_scaling_diagnostics(
        scan::IsometricPlanarCapacityScan,
        diagnostics::IsometricPlanarScalingDiagnostics)
    figure = CairoMakie.Figure(size=(1450, 1000), backgroundcolor=:white)
    if diagnostics.fit.status != :ok || isempty(diagnostics.loss_values)
        axis = CairoMakie.Axis(
            figure[1, 1]; title="Finite-size-scaling diagnostics")
        CairoMakie.hidedecorations!(axis)
        CairoMakie.hidespines!(axis)
        CairoMakie.text!(axis, 0.5, 0.5;
            text="diagnostics unavailable ($(diagnostics.fit.status))",
            space=:relative, align=(:center, :center),
            fontsize=22, color=:gray30)
        return figure
    end
    _plot_isometric_planar_loss_diagnostics!(figure[1, 1], diagnostics)
    _plot_isometric_planar_bootstrap_diagnostics!(figure[1, 2], diagnostics)
    _plot_isometric_planar_batch_diagnostics!(figure[2, 1], scan, diagnostics)
    _plot_isometric_planar_sensitivity_diagnostics!(figure[2, 2], diagnostics)
    return figure
end
