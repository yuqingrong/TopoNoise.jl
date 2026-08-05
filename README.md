# TopoNoise

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://yuqingrong.github.io/TopoNoise.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://yuqingrong.github.io/TopoNoise.jl/dev/)
[![Build Status](https://github.com/yuqingrong/TopoNoise.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/yuqingrong/TopoNoise.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/yuqingrong/TopoNoise.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/yuqingrong/TopoNoise.jl)

## Finite toric-code PEPS

`toric_code_peps` constructs the doubled-edge PEPS on a rectangular open
patch. Every site has four physical qubits in east/north/west/south order,
matching `(i, j, k, l)` in the local-tensor diagram.
Virtual indices are shared across internal edges, so a bond state `|i⟩` is
copied to the physical state `|i i⟩` at its endpoints. Dangling virtual legs
are terminated by the x-polarized vector `(|0⟩ + |1⟩)/√2`, and the resulting
finite state is normalized.

```julia
using TopoNoise
using ITensors

peps = toric_code_peps(3, 3)

# Inspect the center tensor and its four physical qubit indices.
center = peps[2, 2]
east, north, west, south = physicalinds(peps, 2, 2)

@show size(peps) inds(center)
```

The underlying normalized rank-eight parity tensor is available separately.
Its virtual-leg order is `(α, β, γ, δ) = (east, north, west, south)`:

```julia
local_tensor = toric_code_local_tensor()
@assert size(local_tensor) == ntuple(_ -> 2, 8)
```

## PEPS graph and diagonal sequential circuit

`peps_graph` converts the same object to an `ITensorNetwork` keyed by
`(row, col)`. Shared virtual indices become graph edges, while all four
directional physical indices remain external legs. The sequential model uses
one horizontal carrier per row and one vertical carrier per column; its Yao
unitary applies sites along diagonals starting at the lower-left corner.

```julia
using TopoNoise

peps = toric_code_peps(3, 3)
network = peps_graph(peps)
circuit = sequential_circuit_graph(peps)

# [(3,1)]; [(3,2),(2,1)]; [(3,3),(2,2),(1,1)]; ...
layers = circuit_layers(circuit)
core = yao_unitary(circuit)

@show log2_postselection_probability(circuit)
@show postselection_probability(circuit)
```

`yao_circuit` adds the carrier `|+⟩` preparation to that core and returns an
executable `Yao.ChainBlock`. Apply it to an all-zero register or draw that exact
block with Yao's native circuit renderer:

```julia
using Yao

block = yao_circuit(peps) # yao_circuit(circuit) is equivalent
register = zero_state(nqubits(block))
apply!(register, block)

Yao.plot(block)
Yao.vizcircuit(block; filename="toric_yao_circuit.svg")
Yao.vizcircuit(block; format=:pdf, filename="toric_yao_circuit.pdf")
Yao.vizcircuit(block; format=:png, filename="toric_yao_circuit.png")
```

### Runnable example

From the repository root, instantiate the project once and run the generator:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. examples/generate_toric_circuit.jl 2 2
```

Omit `2 2` to use the default lattice size. The script prints the diagonal
schedule and writes SVG, PDF, and PNG files to `examples/output/`. It builds
and draws the circuit without allocating the exponentially large state vector.

Yao's native view uses one line per qubit and its standard circuit layout. For
a publication-style diagram with grouped physical buses, diagonal layer bands,
and carrier annotations, `plot_sequential_circuit` returns a
`CairoMakie.Figure`. Save the CairoMakie figures using Makie's standard `save`
function:

```julia
using CairoMakie

peps_figure = plot_peps_graph(peps; show_index_labels=true)
circuit_figure = plot_sequential_circuit(circuit)

save("toric_peps.svg", peps_figure)
save("toric_peps.pdf", peps_figure)
save("toric_sequential_circuit.svg", circuit_figure)
save("toric_sequential_circuit.pdf", circuit_figure)

# To show all four physical qubit wires instead of one grouped bus:
expanded = plot_sequential_circuit(circuit; expand_physical_buses=true)
```

`yao_unitary(circuit)` contains only the reusable local-gate core, while
`yao_circuit(circuit)` also prepares the carriers in `|+⟩`. The final
carrier `⟨+|` projections remain explicit model metadata and are not part of
the executable unitary block. The global circuit matrix is never materialized.
