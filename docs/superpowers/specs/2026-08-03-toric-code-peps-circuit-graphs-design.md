# Toric-Code PEPS and Sequential-Circuit Graphs Design

## Goal

Turn a finite rectangular `ToricCodePEPS` into both an attributed tensor-network
graph and an exact postselected sequential-circuit model, then render both as
publication-quality static figures.

## Package roles

- `ITensorNetworks` stores the PEPS topology with `(row, col)` vertex keys.
- `Yao` stores the unitary core as embedded six-qubit local gates.
- `CairoMakie` renders fixed-layout PEPS and grouped-wire circuit figures.

All three packages are direct dependencies. The first supported input is only
`ToricCodePEPS`; arbitrary tensor-network inference and hardware decomposition
are out of scope.

## Network and circuit conventions

The PEPS graph preserves shared virtual indices as internal edges and the four
E/N/W/S physical indices as external legs. Tensor coordinates are
`x = col`, `y = rows - row + 1`.

For an `R × C` PEPS, the circuit has `4RC + R + C` fixed qubit wires. Physical
wires occupy `1:4RC` in row-major site order and E/N/W/S direction order,
followed by one horizontal carrier per row and one vertical carrier per column.
At site `(r, c)`, the local gate acts on
`(pE, pN, pW, pS, h_r, v_c)`.

The causal layer is

```math
\ell(r,c)=R-r+c.
```

Rows are descending within a layer, so the sweep begins at `(R, 1)`. Gates in
one layer have distinct row and column carriers and are therefore disjoint.

All carrier qubits start in `|+⟩` and are postselected with `⟨+|` after the
unitary core. With `V = RC` and
`E = R(C-1) + (R-1)C`, the exact success exponent is
`log2(p_success) = E - 2V`. Conditional normalization therefore reproduces the
finite PEPS scale `2^(V-E/2)`.

## Public API

```julia
ToricCodeSequentialCircuit

peps_graph(peps::ToricCodePEPS)
sequential_circuit_graph(peps::ToricCodePEPS)

circuit_layers(circuit::ToricCodeSequentialCircuit)
yao_unitary(circuit::ToricCodeSequentialCircuit)
log2_postselection_probability(circuit::ToricCodeSequentialCircuit)
postselection_probability(circuit::ToricCodeSequentialCircuit)

plot_peps_graph(peps::ToricCodePEPS; show_index_labels=false)
plot_sequential_circuit(circuit::ToricCodeSequentialCircuit;
                        show_layer_labels=true,
                        expand_physical_buses=false)
```

Plotting returns a `CairoMakie.Figure`; file export uses CairoMakie's `save`.

## Rendering

The circuit figure uses conventional left-to-right time. Horizontal carrier
wires appear first, vertical carriers second, and grouped four-qubit physical
buses last. Each diagonal layer is a shaded band. Parallel gates receive stable
subcolumns inside that band, with labels `U[r,c]`. Carrier labels show `W→E` or
`S→N`, physical inputs show `|0⟩`, and carrier endpoints show `|+⟩` and `⟨+|`.

## Verification

Validate rectangular shape, directional physical indices, expected nearest-
neighbor shared indices, and absence of unexpected graph edges. Test topology
and schedules through `3×3`, verify local Yao basis ordering, compare exact
postselected states for `1×1`, `1×2`, and `2×1`, render SVG/PDF fixtures, and
keep all pre-existing tests passing. Never form the global circuit matrix.
