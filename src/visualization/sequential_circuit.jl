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

