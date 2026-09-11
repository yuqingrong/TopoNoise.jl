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

## Rotated planar logical-failure scans

The rotated planar implementation uses the authoritative CSS convention

```math
A_s = \prod_{q \in s} Z_q, \qquad B_p = \prod_{q \in p} X_q.
```

It supports one logical qubit on odd square patches with distance `d >= 3`.
The `construction` keyword selects how the Clifford encoder is synthesized;
it is independent of `logical_state`. Both `construction=:as` and `:bp`
prepare the same requested state from `:zero`, `:one`, `:plus`, or `:minus`.
The default is `logical_state=:zero`, namely `|0_L>`.

Conceptually, the two default-state constructions are the projector formulas

```math
|0_L\rangle \propto (I + \bar Z)\prod_s(I + A_s)|+\rangle^{\otimes n},
\qquad
|0_L\rangle \propto \prod_p(I + B_p)|0\rangle^{\otimes n}.
```

The stored encoder is a deterministic local Clifford schedule implementing
these states without materializing either projector. At `d=3` it replays the
literal As/Bp plaquette circuits; larger patches use deterministic reverse
shelling so each plaquette obtains a fresh local representative. The Bp core
directly prepares `|0_L>`; the As core is converted locally from its natural
`|+_L>` preparation to the requested logical basis. `clock=:gate_layer`
injects independent X and Z faults after every elementary encoder layer;
`clock=:plaquette` injects once after every completed source-check block.
Both models use ideal preparation and a final perfect stabilizer measurement:
there are no repeated syndrome rounds or measurement errors.

### PyMatching setup

The production estimator sends only final CSS syndrome batches to independent
PyMatching decoders. Python and PyMatching are managed by `CondaPkg.toml`:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'import CondaPkg; CondaPkg.resolve()'
```

`sample_yao_syndrome` is a small state-vector oracle restricted to `d=3` so
that accidental exponential allocations fail early. Logical-failure estimates
and scans use the binary Pauli-frame engine and support every implemented odd
distance.

### Run a seeded scan

From Julia, pass an RNG explicitly. The optional `seed` keyword records the
seed in the scan and CSV metadata; it never reseeds or replaces that RNG.

```julia
using Random
using TopoNoise

scan = scan_logical_failure(
    MersenneTwister(1234);
    distances=[3, 5, 7], error_rates=collect(0:0.01:0.1),
    shots=10_000, batch_size=10_000, seed=1234,
)
paths = save_logical_failure_scan(scan, "results/rotated-planar")
```

The convenience scan uses equal independent physical rates, `p_x == p_z`, at
each grid point. Use `estimate_logical_failure` directly with
`CircuitPauliNoise(p; p_x=..., p_z=...)` for an unequal-rate experiment.

The example defaults to distances `3,5,7`, physical error rates
`0:0.01:0.1`, `|0_L>`, the B_p construction, north/south X boundaries, the
gate-layer noise clock, 10,000 shots per point, and seed 1234:

```bash
julia --project=. examples/scan_rotated_planar.jl
```

An explicit small scan looks like:

```bash
julia --project=. examples/scan_rotated_planar.jl \
  --distances 3,5,7 --error-rates 0,0.01,0.02 \
  --state zero --construction bp --boundary-orientation x_ns \
  --clock gate_layer --shots 10000 --batch-size 10000 --seed 1234 \
  --output-dir results/rotated-planar --basename logical-failure
```

Use `--p-min`, `--p-max`, and `--p-step` instead of `--error-rates` to build
an evenly spaced grid, and run with `--help` for all options. The command
creates the output directory and writes a tidy CSV plus matching SVG, PDF, and
PNG plots. The three plot panels report logical-X, logical-Z, and
either-logical failure rates, with one error-bar curve per distance.

## Construction/channel scaling comparison

Compare the two encoder constructions under separate single-channel noise with
the seeded production command:

```bash
julia --project=. examples/compare_rotated_planar_constructions.jl \
  --distances 3,5,7 --p-min 0 --p-max 0.16 --p-step 0.01 \
  --shots 10000 --batch-size 10000 --seed 1234 \
  --output-dir results/rotated-planar-comparison
```

All panels prepare `|0_L>`. X-only noise means `p_x=p,p_z=0` and reports the
logical-X failure rate; Z-only noise means `p_x=0,p_z=p` and reports the
logical-Z failure rate. `As` and `Bp` choose the circuit construction only;
they do not choose a different logical state.

The command writes one combined 2×2 figure and four standalone figures
(`As`/`Bp` × X-only/Z-only) in SVG, PDF, and PNG, along with raw and fit CSV
files. When a fit is available, each panel includes an inset with its
exploratory `p_c` and `nu` collapse. The `d=3,5,7` fit values are exploratory
diagnostics, not threshold claims. Use `--no-fit` to preserve the raw curves
and write explicit unavailable-fit records without attempting the fit.

Version 1 intentionally excludes repeated syndrome rounds, measurement noise,
periodic layouts, holes, and multi-logical-qubit patches. Its three-distance
scaling fits are exploratory diagnostics, not threshold claims.
