function _trajectory_scan_points(scan::TrajectoryScan, size::Int)
    return [_scan_point(scan, size, rate) for rate in scan.error_rates]
end

"""
    plot_trajectory_scan(scan, crossings=[])

Render the injected-error diagnostics, plaquette frustration, largest
occupied internal-bond cluster, and horizontal spanning probability for a
finite-size trajectory scan. Crossing estimates are shown when available.
"""
function plot_trajectory_scan(
        scan::TrajectoryScan,
        crossings::AbstractVector{<:CriticalCrossing}=CriticalCrossing[])
    figure = CairoMakie.Figure(size=(1100, 820), backgroundcolor=:white)
    injected_axis = CairoMakie.Axis(
        figure[1, 1]; title="Injected errors and mismatches",
        xlabel="virtual-bond error rate p", ylabel="density")
    frustration_axis = CairoMakie.Axis(
        figure[1, 2]; title="Plaquette frustration",
        xlabel="virtual-bond error rate p", ylabel="odd-plaquette density")
    cluster_axis = CairoMakie.Axis(
        figure[2, 1]; title="Largest occupied cluster",
        xlabel="virtual-bond error rate p", ylabel="fraction of sites")
    spanning_axis = CairoMakie.Axis(
        figure[2, 2]; title="Horizontal spanning probability",
        xlabel="virtual-bond error rate p", ylabel="probability")

    colors = ("#0072B2", "#D55E00", "#009E73", "#CC79A7",
              "#E69F00", "#56B4E9", "#000000")
    rates = scan.error_rates
    CairoMakie.lines!(injected_axis, rates, rates;
        color=:gray35, linestyle=:dot, linewidth=2, label="y = p")
    frustration_theory =
        0.5 .* (1 .- (1 .- 2 .* rates) .^ 4)
    CairoMakie.lines!(frustration_axis, rates, frustration_theory;
        color=:gray35, linestyle=:dot, linewidth=2,
        label="½[1-(1-2p)⁴]")

    for (size_index, size) in enumerate(scan.sizes)
        points = _trajectory_scan_points(scan, size)
        color = colors[mod1(size_index, length(colors))]
        total = [point.sampled_error_mean for point in points]
        total_se = [point.sampled_error_se for point in points]
        mismatch = [point.mismatch_mean for point in points]
        mismatch_se = [point.mismatch_se for point in points]
        frustration = [point.frustration_mean for point in points]
        frustration_se = [point.frustration_se for point in points]
        largest = [point.largest_cluster_mean for point in points]
        largest_se = [point.largest_cluster_se for point in points]
        spanning = [point.horizontal_spanning_mean for point in points]
        spanning_se = [point.horizontal_spanning_se for point in points]
        CairoMakie.lines!(injected_axis, rates, total;
            color=color, linewidth=2, label="L=$size injected")
        CairoMakie.scatter!(injected_axis, rates, total; color=color)
        CairoMakie.errorbars!(injected_axis, rates, total, total_se;
            color=color, whiskerwidth=7)
        CairoMakie.lines!(injected_axis, rates, mismatch;
            color=color, linewidth=2, linestyle=:dash,
            label="L=$size mismatch")
        CairoMakie.errorbars!(injected_axis, rates, mismatch, mismatch_se;
            color=color, whiskerwidth=7)

        for (axis, values, errors) in (
                (frustration_axis, frustration, frustration_se),
                (cluster_axis, largest, largest_se),
                (spanning_axis, spanning, spanning_se))
            CairoMakie.lines!(axis, rates, values;
                color=color, linewidth=2, label="L=$size")
            CairoMakie.scatter!(axis, rates, values; color=color)
            CairoMakie.errorbars!(axis, rates, values, errors;
                color=color, whiskerwidth=7)
        end
    end

    for crossing in crossings
        crossing.status == :ok || continue
        ismissing(crossing.estimate) && continue
        CairoMakie.vlines!(spanning_axis, [crossing.estimate];
            color=(:gray25, 0.55), linestyle=:dash, linewidth=1.5)
    end
    for axis in (
            injected_axis, frustration_axis, cluster_axis, spanning_axis)
        CairoMakie.ylims!(axis, -0.03, 1.03)
        CairoMakie.axislegend(axis; position=:lt, framevisible=false)
    end
    return figure
end

