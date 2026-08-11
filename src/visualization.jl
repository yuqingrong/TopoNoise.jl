const _PEPS_NODE_COLOR = "#315B7D"
const _PEPS_BOND_COLOR = "#4C566A"
const _PEPS_PHYSICAL_COLOR = "#C44E52"
const _CIRCUIT_CARRIER_COLOR = "#3B6F8F"
const _CIRCUIT_PHYSICAL_COLOR = "#6B7280"
const _CIRCUIT_GATE_COLOR = "#D88C4B"

_peps_position(rows::Int, site::Tuple{Int,Int}) =
    (Float64(site[2]), Float64(rows - site[1] + 1))

function _compact_index_label(index::Index, direction::String)
    hexadecimal_id = string(ITensors.id(index); base=16)
    suffix = last(hexadecimal_id, min(6, length(hexadecimal_id)))
    return "$direction:#$suffix"
end

function _finish_graph_axis!(axis, xmin, xmax, ymin, ymax)
    CairoMakie.xlims!(axis, xmin, xmax)
    CairoMakie.ylims!(axis, ymin, ymax)
    CairoMakie.hidedecorations!(axis)
    CairoMakie.hidespines!(axis)
    return axis
end

"""
    plot_peps_graph(peps::ToricCodePEPS; show_index_labels=false)

Render a publication-ready PEPS graph. Site `(row, col)` is placed at
`(col, R-row+1)`, so matrix row one appears at the top. Virtual bonds join
tensor nodes and the four physical indices are shown as directional stubs.
When `show_index_labels=true`, each stub includes a compact suffix of its
actual ITensor index identifier.
"""
function plot_peps_graph(
        peps::ToricCodePEPS; show_index_labels::Bool=false)
    network = peps_graph(peps)
    rows, cols = size(peps)
    figure = CairoMakie.Figure(
        size=(max(480, 150 * cols), max(420, 150 * rows)),
        backgroundcolor=:white)
    axis = CairoMakie.Axis(figure[1, 1], aspect=CairoMakie.DataAspect())

    for edge in ITensorNetworks.edges(network)
        first_site, second_site = Tuple(edge)
        x1, y1 = _peps_position(rows, first_site)
        x2, y2 = _peps_position(rows, second_site)
        CairoMakie.lines!(axis, [x1, x2], [y1, y2];
            color=_PEPS_BOND_COLOR, linewidth=3)
    end

    stub_geometry = (
        ((0.22, -0.12), (0.63, -0.12), (0.43, -0.24), "E"),
        ((0.12, 0.22), (0.12, 0.63), (0.25, 0.43), "N"),
        ((-0.22, 0.12), (-0.63, 0.12), (-0.43, 0.24), "W"),
        ((-0.12, -0.22), (-0.12, -0.63), (-0.25, -0.43), "S"),
    )
    for row in 1:rows, col in 1:cols
        x, y = _peps_position(rows, (row, col))
        site_physical = physicalinds(peps, row, col)
        for (direction_number, geometry) in enumerate(stub_geometry)
            start_offset, end_offset, label_position, direction = geometry
            start_x, start_y = x + start_offset[1], y + start_offset[2]
            end_x, end_y = x + end_offset[1], y + end_offset[2]
            CairoMakie.lines!(axis, [start_x, end_x], [start_y, end_y];
                color=_PEPS_PHYSICAL_COLOR, linewidth=2.2)
            CairoMakie.scatter!(axis, [end_x], [end_y];
                color=_PEPS_PHYSICAL_COLOR, markersize=7)
            show_index_labels && CairoMakie.text!(
                axis, x + label_position[1], y + label_position[2];
                text=_compact_index_label(
                    site_physical[direction_number], direction),
                color=_PEPS_PHYSICAL_COLOR, fontsize=10,
                align=(:center, :center))
        end
        CairoMakie.scatter!(axis, [x], [y]; marker=:rect, markersize=58,
            color=_PEPS_NODE_COLOR, strokecolor=:white, strokewidth=1.5)
        CairoMakie.text!(axis, x, y; text="T[$row,$col]", color=:white,
            fontsize=13, align=(:center, :center))
    end

    _finish_graph_axis!(axis, 0.15, cols + 0.85, 0.15, rows + 0.85)
    return figure
end

function _circuit_wire_layout(
        circuit::ToricCodeSequentialCircuit, expand_physical_buses::Bool)
    entries = NamedTuple[]
    for row in 1:circuit.rows
        push!(entries, (kind=:horizontal, label="h$row  W→E",
            site=nothing, direction=nothing,
            wires=(circuit.horizontal_wires[row],)))
    end
    for col in 1:circuit.cols
        push!(entries, (kind=:vertical, label="v$col  S→N",
            site=nothing, direction=nothing,
            wires=(circuit.vertical_wires[col],)))
    end
    directions = ("E", "N", "W", "S")
    for row in 1:circuit.rows, col in 1:circuit.cols
        wires = _physical_wires(circuit.cols, row, col)
        if expand_physical_buses
            for direction in 1:4
                push!(entries, (kind=:physical,
                    label="p[$row,$col,$(directions[direction])]",
                    site=(row, col), direction=directions[direction],
                    wires=(wires[direction],)))
            end
        else
            push!(entries, (kind=:physical,
                label="p[$row,$col]  (E,N,W,S) ×4", site=(row, col),
                direction=nothing, wires=wires))
        end
    end

    wire_y = Dict{Int,Float64}()
    count = length(entries)
    for (position, entry) in enumerate(entries)
        y = Float64(count - position + 1)
        for wire in entry.wires
            wire_y[wire] = y
        end
    end
    return entries, wire_y
end

function _layer_geometry(circuit::ToricCodeSequentialCircuit)
    starts = Float64[]
    ends = Float64[]
    gate_x = Dict{Tuple{Int,Int},Float64}()
    cursor = 1.8
    for layer in circuit.layers
        width = max(1, length(layer)) + 0.8
        push!(starts, cursor)
        push!(ends, cursor + width)
        for (subcolumn, gate) in enumerate(layer)
            gate_x[gate.site] = cursor + 0.4 + subcolumn - 0.5
        end
        cursor += width + 0.25
    end
    return starts, ends, gate_x, cursor
end

"""
    plot_sequential_circuit(circuit::ToricCodeSequentialCircuit;
                            show_layer_labels=true,
                            expand_physical_buses=false)

Render the sequential circuit from left to right. Carrier wires appear above
row-major physical buses; every diagonal layer has a shaded band and its
parallel site gates occupy deterministic subcolumns. Set
`expand_physical_buses=true` to draw all four physical qubits separately.
"""
function plot_sequential_circuit(
        circuit::ToricCodeSequentialCircuit;
        show_layer_labels::Bool=true,
        expand_physical_buses::Bool=false)
    entries, wire_y = _circuit_wire_layout(circuit, expand_physical_buses)
    starts, ends, gate_x, circuit_end = _layer_geometry(circuit)
    wire_count = length(entries)
    figure = CairoMakie.Figure(
        size=(max(900, round(Int, 115 * circuit_end)),
              max(440, 42 * wire_count + 120)),
        backgroundcolor=:white)
    axis = CairoMakie.Axis(figure[1, 1])
    ylow, yhigh = 0.35, wire_count + 0.65

    for layer_number in eachindex(starts)
        color = isodd(layer_number) ? (0.85, 0.91, 0.96, 0.42) :
                                      (0.93, 0.88, 0.82, 0.42)
        rectangle = CairoMakie.Rect2f(
            starts[layer_number], ylow,
            ends[layer_number] - starts[layer_number], yhigh - ylow)
        CairoMakie.poly!(axis, rectangle; color=color, strokewidth=0)
        show_layer_labels && CairoMakie.text!(
            axis, (starts[layer_number] + ends[layer_number]) / 2,
            wire_count + 1.02; text="layer $layer_number", fontsize=13,
            color=:gray35, align=(:center, :bottom))
    end

    input_x, output_x = 0.7, circuit_end + 0.15
    wire_start, wire_end = input_x + 0.58, output_x - 0.58
    for (position, entry) in enumerate(entries)
        y = Float64(wire_count - position + 1)
        line_color = entry.kind == :physical ?
                     _CIRCUIT_PHYSICAL_COLOR : _CIRCUIT_CARRIER_COLOR
        CairoMakie.lines!(axis, [wire_start, wire_end], [y, y];
            color=line_color, linewidth=entry.kind == :physical ? 1.8 : 2.2)
        CairoMakie.text!(axis, input_x - 0.18, y; text=entry.label,
            align=(:right, :center), fontsize=13, color=:gray20)
        if entry.kind == :physical
            CairoMakie.text!(axis, input_x + 0.08, y; text="|0⟩",
                align=(:left, :center), fontsize=13, color=:gray20)
        else
            CairoMakie.text!(axis, input_x + 0.08, y; text="|+⟩",
                align=(:left, :center), fontsize=13,
                color=_CIRCUIT_CARRIER_COLOR)
            CairoMakie.text!(axis, output_x - 0.08, y; text="⟨+|",
                align=(:right, :center), fontsize=13,
                color=_CIRCUIT_CARRIER_COLOR)
        end
    end

    for layer in circuit.layers, gate in layer
        x = gate_x[gate.site]
        ys = sort(unique(wire_y[wire] for wire in gate.wires))
        CairoMakie.lines!(axis, fill(x, 2), [first(ys), last(ys)];
            color=_CIRCUIT_GATE_COLOR, linewidth=2.5)
        CairoMakie.scatter!(axis, fill(x, length(ys)), ys;
            marker=:rect, markersize=16, color=_CIRCUIT_GATE_COLOR,
            strokecolor=:white, strokewidth=1)
        CairoMakie.text!(axis, x, (first(ys) + last(ys)) / 2;
            text="U[$(gate.site[1]),$(gate.site[2])]", fontsize=12,
            rotation=pi / 2, color=:black, align=(:center, :bottom))
    end

    _finish_graph_axis!(
        axis, -2.0, output_x + 0.35, 0.05,
        wire_count + (show_layer_labels ? 1.45 : 0.9))
    return figure
end

function _trajectory_scan_points(scan::TrajectoryScan, size::Int)
    return [_scan_point(scan, size, rate) for rate in scan.error_rates]
end

function _marginal_spin_scan_points(scan::TrajectoryScan, size::Int)
    return [
        _marginal_spin_point(scan, size, rate) for rate in scan.error_rates]
end

"""
    plot_trajectory_scan(scan, crossings; binder_crossings=[])

Render the injected-error diagnostics, plaquette frustration, largest
occupied internal-bond cluster, horizontal spanning probability, marginalized
absolute magnetization, and Binder cumulant for a finite-size trajectory scan.
Crossing estimates are shown on their corresponding panels when available.
"""
function plot_trajectory_scan(
        scan::TrajectoryScan,
        crossings::AbstractVector{<:CriticalCrossing}=CriticalCrossing[];
        binder_crossings::AbstractVector{<:CriticalCrossing}=
            CriticalCrossing[])
    has_marginal_spins = !isempty(scan.marginal_spin_points)
    figure = CairoMakie.Figure(
        size=has_marginal_spins ? (1100, 1180) : (1100, 820),
        backgroundcolor=:white)
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
    magnetization_axis = has_marginal_spins ? CairoMakie.Axis(
        figure[3, 1]; title="Marginalized |M|",
        xlabel="virtual-bond error rate p", ylabel="⟨|M|⟩") : nothing
    binder_axis = has_marginal_spins ? CairoMakie.Axis(
        figure[3, 2]; title="Binder cumulant U₄",
        xlabel="virtual-bond error rate p", ylabel="U₄") : nothing

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
        if has_marginal_spins
            spin_points = _marginal_spin_scan_points(scan, size)
            magnetization = [
                point.absolute_magnetization_mean for point in spin_points]
            magnetization_se = [
                point.absolute_magnetization_se for point in spin_points]
            binder = [point.binder_cumulant for point in spin_points]
            CairoMakie.lines!(magnetization_axis, rates, magnetization;
                color=color, linewidth=2, label="L=$size")
            CairoMakie.scatter!(
                magnetization_axis, rates, magnetization; color=color)
            CairoMakie.errorbars!(
                magnetization_axis, rates, magnetization, magnetization_se;
                color=color, whiskerwidth=7)

            CairoMakie.lines!(binder_axis, rates, binder;
                color=color, linewidth=2, label="L=$size")
            CairoMakie.scatter!(binder_axis, rates, binder; color=color)
            binder_error_indices = findall(
                point -> !ismissing(point.binder_cumulant_se), spin_points)
            if !isempty(binder_error_indices)
                CairoMakie.errorbars!(
                    binder_axis, rates[binder_error_indices],
                    binder[binder_error_indices],
                    Float64[spin_points[index].binder_cumulant_se
                            for index in binder_error_indices];
                    color=color, whiskerwidth=7)
            end
        end
    end

    for crossing in crossings
        crossing.status == :ok || continue
        ismissing(crossing.estimate) && continue
        CairoMakie.vlines!(spanning_axis, [crossing.estimate];
            color=(:gray25, 0.55), linestyle=:dash, linewidth=1.5)
    end
    if has_marginal_spins
        for crossing in binder_crossings
            crossing.status == :ok || continue
            ismissing(crossing.estimate) && continue
            CairoMakie.vlines!(binder_axis, [crossing.estimate];
                color=(:gray25, 0.55), linestyle=:dash, linewidth=1.5)
        end
    end
    for axis in (
            injected_axis, frustration_axis, cluster_axis, spanning_axis)
        CairoMakie.ylims!(axis, -0.03, 1.03)
        CairoMakie.axislegend(axis; position=:lt, framevisible=false)
    end
    if has_marginal_spins
        CairoMakie.ylims!(magnetization_axis, -0.03, 1.03)
        CairoMakie.axislegend(
            magnetization_axis; position=:lt, framevisible=false)
        binder_lower = minimum(
            point.binder_cumulant -
            coalesce(point.binder_cumulant_se, 0.0)
            for point in scan.marginal_spin_points)
        binder_upper = maximum(
            point.binder_cumulant +
            coalesce(point.binder_cumulant_se, 0.0)
            for point in scan.marginal_spin_points)
        binder_padding = max(0.03, 0.05 * (binder_upper - binder_lower))
        CairoMakie.ylims!(
            binder_axis, binder_lower - binder_padding,
            max(0.70, binder_upper + binder_padding))
        CairoMakie.axislegend(binder_axis; position=:lt, framevisible=false)
    end
    return figure
end

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

"""
    plot_open_code_capacity(scan, ns_crossings; ew_crossings=[])

Render north-south and east-west logical-failure rates for an independent
open code-capacity scan. Standard-error bars are shown at each scan point,
and successful adjacent-size crossing estimates are marked by vertical lines.
"""
function plot_open_code_capacity(
        scan::OpenCodeCapacityScan,
        ns_crossings::AbstractVector{<:CriticalCrossing}=CriticalCrossing[];
        ew_crossings::AbstractVector{<:CriticalCrossing}=CriticalCrossing[])
    figure = CairoMakie.Figure(size=(1100, 470), backgroundcolor=:white)
    ns_axis = CairoMakie.Axis(
        figure[1, 1]; title="Logical failure (N-S)",
        xlabel="data-edge error rate p", ylabel="p_fail")
    ew_axis = CairoMakie.Axis(
        figure[1, 2]; title="Logical failure (E-W)",
        xlabel="data-edge error rate p", ylabel="p_fail")
    colors = ("#0072B2", "#D55E00", "#009E73", "#CC79A7",
              "#E69F00", "#56B4E9", "#000000")

    for (size_index, size) in enumerate(scan.sizes)
        points = [_open_code_capacity_scan_point(
            scan, size, rate) for rate in scan.error_rates]
        color = colors[mod1(size_index, length(colors))]
        for (axis, values, errors) in (
                (ns_axis,
                 [point.logical_failure_ns_mean for point in points],
                 [point.logical_failure_ns_se for point in points]),
                (ew_axis,
                 [point.logical_failure_ew_mean for point in points],
                 [point.logical_failure_ew_se for point in points]))
            CairoMakie.lines!(axis, scan.error_rates, values;
                color=color, linewidth=2, label="L=$size")
            CairoMakie.scatter!(axis, scan.error_rates, values; color=color)
            CairoMakie.errorbars!(axis, scan.error_rates, values, errors;
                color=color, whiskerwidth=7)
        end
    end

    for (axis, crossings) in (
            (ns_axis, ns_crossings), (ew_axis, ew_crossings))
        for crossing in crossings
            crossing.status == :ok || continue
            ismissing(crossing.estimate) && continue
            CairoMakie.vlines!(axis, [crossing.estimate];
                color=(:gray25, 0.55), linestyle=:dash, linewidth=1.5)
        end
        CairoMakie.ylims!(axis, -0.03, 1.03)
        CairoMakie.axislegend(axis; position=:lt, framevisible=false)
    end
    return figure
end
