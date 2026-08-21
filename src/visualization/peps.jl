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

