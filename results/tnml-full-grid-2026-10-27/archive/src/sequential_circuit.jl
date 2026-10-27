struct _SequentialSiteGate
    site::Tuple{Int,Int}
    layer::Int
    wires::NTuple{6,Int}
end

"""
    ToricCodeSequentialCircuit

A diagonal, lower-left-to-upper-right preparation schedule for a finite
`ToricCodePEPS`. The Yao block returned by [`yao_unitary`](@ref) contains the
unitary local-gate core. Carrier preparation in `|+⟩` and final `⟨+|`
postselection are stored as the horizontal and vertical carrier wire ranges
and are shown explicitly by [`plot_sequential_circuit`](@ref).
"""
struct ToricCodeSequentialCircuit{Block}
    rows::Int
    cols::Int
    physical_wire_count::Int
    horizontal_wires::UnitRange{Int}
    vertical_wires::UnitRange{Int}
    carrier_input_state::Symbol
    carrier_output_effect::Symbol
    layers::Vector{Vector{_SequentialSiteGate}}
    unitary::Block
    log2_success_probability::Int
end

function _physical_wires(cols::Int, row::Int, col::Int)
    first_wire = 4 * ((row - 1) * cols + col - 1) + 1
    return ntuple(offset -> first_wire + offset - 1, 4)
end

function _site_schedule(rows::Int, cols::Int)
    return [
        [(row, layer - rows + row)
         for row in rows:-1:1
         if 1 <= layer - rows + row <= cols]
        for layer in 1:(rows + cols - 1)
    ]
end

function _yao_local_gate(matrix::AbstractMatrix, site::Tuple{Int,Int})
    row, col = site
    return matblock(matrix; tag="U[$row,$col]")
end

_yao_local_gate(::Type{T}, site::Tuple{Int,Int}) where {T<:Number} =
    _yao_local_gate(reshape(toric_code_local_gate(T), 64, 64), site)

"""
    sequential_circuit_graph(peps::ToricCodePEPS)

Construct the diagonal sequential-circuit model for `peps`. Physical wires
come first, in row-major site order and `(E, N, W, S)` order within a site.
They are followed by one horizontal carrier per row and one vertical carrier
per column. Gates in layer `R - row + col` are ordered by descending row.
"""
function sequential_circuit_graph(peps::ToricCodePEPS)
    _validate_peps_graph_input(peps)
    rows, cols = size(peps)
    vertices = rows * cols
    physical_wire_count = 4 * vertices
    horizontal_wires = (physical_wire_count + 1):(physical_wire_count + rows)
    vertical_wires =
        (physical_wire_count + rows + 1):(physical_wire_count + rows + cols)
    wire_count = physical_wire_count + rows + cols
    scheduled_sites = _site_schedule(rows, cols)

    layers = Vector{Vector{_SequentialSiteGate}}(undef, length(scheduled_sites))
    scalar_type = eltype(peps[1, 1])
    local_matrix = reshape(toric_code_local_gate(scalar_type), 64, 64)
    for (layer_number, sites) in enumerate(scheduled_sites)
        layer = _SequentialSiteGate[]
        for site in sites
            row, col = site
            wires = (
                _physical_wires(cols, row, col)...,
                horizontal_wires[row],
                vertical_wires[col],
            )
            push!(layer, _SequentialSiteGate(site, layer_number, wires))
        end
        layers[layer_number] = layer
    end

    blocks = [
        subroutine(
            wire_count, _yao_local_gate(local_matrix, gate.site), gate.wires)
        for layer in layers for gate in layer
    ]
    unitary = chain(wire_count, blocks)
    internal_edges = rows * (cols - 1) + (rows - 1) * cols
    log2_success_probability = internal_edges - 2 * vertices
    return ToricCodeSequentialCircuit(
        rows, cols, physical_wire_count, horizontal_wires, vertical_wires,
        :plus, :plus_projection, layers, unitary, log2_success_probability)
end

"""
    circuit_layers(circuit::ToricCodeSequentialCircuit)

Return the scheduled site coordinates, grouped into diagonal circuit layers.
"""
circuit_layers(circuit::ToricCodeSequentialCircuit) =
    [[gate.site for gate in layer] for layer in circuit.layers]

"""Return the Yao block for the circuit's unitary local-gate core."""
yao_unitary(circuit::ToricCodeSequentialCircuit) = circuit.unitary

"""
    yao_circuit(circuit::ToricCodeSequentialCircuit)
    yao_circuit(peps::ToricCodePEPS)

Return an executable Yao circuit that, assuming the carrier wires initially
contain `|0⟩`, prepares every horizontal and vertical carrier in `|+⟩` and
then applies the diagonal local-gate core. The returned block does not include
carrier measurement or postselection.
"""
function yao_circuit(circuit::ToricCodeSequentialCircuit)
    wire_count = circuit.physical_wire_count +
                 length(circuit.horizontal_wires) +
                 length(circuit.vertical_wires)
    carrier_wires = first(circuit.horizontal_wires):last(circuit.vertical_wires)
    carrier_preparation = repeat(wire_count, H, carrier_wires)
    return chain(wire_count, carrier_preparation, yao_unitary(circuit))
end

yao_circuit(peps::ToricCodePEPS) = yao_circuit(sequential_circuit_graph(peps))

"""Return the exactly stored base-two logarithm of carrier success probability."""
log2_postselection_probability(circuit::ToricCodeSequentialCircuit) =
    circuit.log2_success_probability

"""
Return the carrier postselection success probability as a floating-point
number. For very large lattices this value may underflow; use
[`log2_postselection_probability`](@ref) as the authoritative exact value.
"""
postselection_probability(circuit::ToricCodeSequentialCircuit) =
    exp2(float(log2_postselection_probability(circuit)))
