function _rotated_crossing_summary(crossing::CriticalCrossing)
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

"""
    plot_rotated_code_capacity(scan, crossings)

Render decoded logical-X failure rates for an open rotated planar code under
data-X noise and a text summary of every adjacent-distance crossing result.
"""
function plot_rotated_code_capacity(
        scan::RotatedCodeCapacityScan,
        crossings::AbstractVector{<:CriticalCrossing}=CriticalCrossing[])
    figure = CairoMakie.Figure(size=(1150, 600), backgroundcolor=:white)
    curve_axis = CairoMakie.Axis(
        figure[1, 1];
        title="Logical X failure — open rotated planar code",
        subtitle="data-X noise; perfect Z syndrome",
        xlabel="data-X error rate p", ylabel="logical X failure rate")
    summary_axis = CairoMakie.Axis(
        figure[1, 2]; title="Adjacent-distance crossing summary")
    colors = ("#0072B2", "#D55E00", "#009E73", "#CC79A7",
              "#E69F00", "#56B4E9", "#000000")

    for (distance_index, distance) in enumerate(scan.distances)
        points = [_rotated_code_capacity_scan_point(
            scan, distance, rate) for rate in scan.error_rates]
        color = colors[mod1(distance_index, length(colors))]
        failures = [point.logical_x_failure_rate for point in points]
        errors = [point.logical_x_failure_se for point in points]
        CairoMakie.lines!(curve_axis, scan.error_rates, failures;
            color=color, linewidth=2, label="d=$distance")
        CairoMakie.scatter!(curve_axis, scan.error_rates, failures;
            color=color)
        CairoMakie.errorbars!(curve_axis, scan.error_rates, failures, errors;
            color=color, whiskerwidth=7)
    end
    for crossing in crossings
        crossing.status == :ok || continue
        ismissing(crossing.estimate) && continue
        CairoMakie.vlines!(curve_axis, [crossing.estimate];
            color=(:gray25, 0.55), linestyle=:dash, linewidth=1.5)
    end
    CairoMakie.ylims!(curve_axis, -0.03, 1.03)
    CairoMakie.axislegend(curve_axis; position=:lt, framevisible=false)

    summary_count = max(1, length(crossings))
    CairoMakie.xlims!(summary_axis, 0.0, 1.0)
    CairoMakie.ylims!(summary_axis, 0.0, summary_count + 1.0)
    CairoMakie.hidedecorations!(summary_axis)
    CairoMakie.hidespines!(summary_axis)
    if isempty(crossings)
        CairoMakie.text!(summary_axis, 0.04, 0.5;
            text="No adjacent-distance crossings available",
            align=(:left, :center), fontsize=17, color=:gray30)
    else
        for (index, crossing) in enumerate(crossings)
            CairoMakie.text!(summary_axis, 0.04, summary_count - index + 1;
                text=_rotated_crossing_summary(crossing),
                align=(:left, :center), fontsize=17, color=:gray20)
        end
    end
    return figure
end

