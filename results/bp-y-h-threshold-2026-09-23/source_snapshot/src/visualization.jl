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
