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
the default `logical_state=:native` uses each circuit's direct preparation:
As prepares `|+_L>` and Bp prepares `|0_L>`. An explicit `logical_state` of
`:zero`, `:one`, `:plus`, or `:minus` overrides that choice.

Conceptually, the two default-state constructions are the projector formulas

```math
|+_L\rangle \propto \prod_s(I + A_s)|+\rangle^{\otimes n},
\qquad
|0_L\rangle \propto \prod_p(I + B_p)|0\rangle^{\otimes n}.
```

The stored encoder is a deterministic Clifford schedule implementing these
states without materializing either projector. Bp follows its diagram at
`d=3`; larger Bp patches use deterministic reverse shelling so each plaquette
obtains a fresh representative. As is derived from this same Bp schedule:
rotate every qubit clockwise by 90 degrees and reverse each CNOT's control
and target, keeping the complete execution order. Boundary checks, trees,
and representatives all follow that mapping at every supported distance
and in both boundary orientations. The old independent As `d=3` template
and independent As shelling are no longer used. The Bp core directly
prepares `|0_L>`.

The default As encoder starts directly in the ideal product state
`|+>^⊗n` (`encoder.input_state == :plus`). Each source-check block applies
H to its fresh representative and then its incoming CNOT tree, preparing
`|+_L>`. Bp starts in ideal `|0>^⊗n` and uses the paired H and outgoing
CNOT tree. Every native H and CNOT noise location is now matched under
rotation and X/Z exchange. An explicit As `:minus` adds a final logical Z.

`yao_encoder(encoder; prepare_input=false)` returns just those stored gates
for a register already in the declared product input. The default
`prepare_input=true` includes ideal input preparation so existing calls on
`zero_state(n)` still work. That ideal prefix is outside the noisy schedule;
the representative H gates inside the schedule remain noisy.

Explicit As `:zero`/`:one` requests retain the separate all-zero-input
parity construction. For `:zero`, one free input stays in `|0>`
and incoming CNOTs set even parity on the input support obtained by
propagating logical Z backwards through the core. `:one` also adds a final
logical X. This optional preparation is separate from the native As circuit.

Gate counts with the default `:x_ns` orientation:

| Distance | Native As (`plus`) | Explicit As `zero` | Native Bp (`zero`) |
|---|---:|---:|---:|
| 3 | 12 (4 H + 8 CNOT) | 13 (4 H + 9 CNOT) | 12 |
| 5 | 40 (12 H + 28 CNOT) | 42 (12 H + 30 CNOT) | 40 |
| 7 | 84 (24 H + 60 CNOT) | 87 (24 H + 63 CNOT) | 84 |

The As parity CNOTs are exposed as a separate `source_check=:logical_parity`
preparation block, with their actual input support and target. At larger
distances they can connect nonadjacent data qubits; the counts above do not
include hardware routing. The plaquette CNOTs stay inside their source
checks. See the [direct As preparation derivation](docs/superpowers/specs/2026-09-12-direct-as-preparation.md)
for the input-parity method; its literal gate examples predate the matched
lattice schedule.

`clock=:gate_layer` injects independent X and Z faults after every elementary
encoder layer, including the initial H and parity CNOTs. `clock=:plaquette`
injects once after each completed source-check block and excludes logical
preparation/correction blocks. Both models use ideal initial product inputs and a
final perfect stabilizer measurement, with no repeated syndrome rounds or
measurement errors. For native states, the complete `:gate_layer` noise,
including H faults, and the `:plaquette` noise exchange logical X and Z
under the lattice rotation when `p_x` and `p_z` are exchanged. Thus As/Z and Bp/X have equal distributions,
as do As/X and Bp/Z. Independent Monte Carlo estimates fluctuate;
properly mapped copies of the same faults give identical decoded outcomes.
Earlier As scans started from all-zero with noisy H gates on free inputs;
those archived results describe a different preparation-noise model.

`clock=:post_encoding` is a separate code-capacity model: prepare the full
encoded state perfectly, then sample one independent X/Z event per data
qubit before the perfect final syndrome round. For IID bit flips only, use:

```julia
noise = CircuitPauliNoise(0.1; p_x=0.1, p_z=0, clock=:post_encoding)
```

There is no preparation, gate, idle, or measurement noise in this model.
For the same geometry and logical state, As and Bp therefore have the same
failure distribution. Replaying the same seeded post-encoding errors gives
identical outcomes; such paired samples must not be counted as independent
extra shots when fitting a threshold.

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

By default (`--logical-state native`), both As panels prepare `|+_L>` directly
and both Bp panels prepare `|0_L>`. Use `--logical-state zero`, `one`, `plus`,
or `minus` to prepare the same explicit state in all panels.
X-only noise means `p_x=p,p_z=0` and reports the
logical X error rate; Z-only noise means `p_x=0,p_z=p` and reports the
logical Z error rate after decoding.

This default `channel_logical` metric diagnoses residual logical parity in the
noise channel. To measure survival of the selected prepared state in both
noise channels, add `--metric state_failure`; it reports residual logical-X
failure for `zero`/`one` and residual logical-Z failure for `plus`/`minus`.
Logical-Z failure can be nonzero on `|0_L>` even though that state is unchanged.
Likewise, logical-X failure can be nonzero on As `|+_L>` while its state
failure is zero. The default plots count those residual logical errors.

The [all-plus-input As comparison](results/rotated-planar-plus-input-2026-09-22/README.md)
uses ideal As `|+>^⊗n` and Bp `|0>^⊗n` inputs with every noisy H/CNOT gate
matched. It resamples As at distances 9, 11, 13, 15 with 50,000 shots per
point and retains the unchanged, independently sampled Bp data.

The archived [matched-CNOT logical-error figure](results/rotated-planar-matched-lattice-2026-09-22/README.md)
uses As `|+_L>` and Bp `|0_L>` at distances 9, 11, 13, 15 with 50,000 shots
per point. Both As panels are resampled with the matched CNOT schedule;
the unchanged Bp data are reused after exact schedule and model checks.
Those As data use the earlier all-zero input with noisy H gates on free
inputs, so they do not represent the current all-plus input. A separate archived
[near-critical fitting analysis](results/rotated-planar-matched-lattice-2026-09-22/critical-window-fit/README.md)
adds joint local fits for As/Z and Bp/X, with bootstrap uncertainty and
window/model/distance sensitivity checks. The earlier
[native-state figure](results/rotated-planar-native-logical-errors-2026-09-22/README.md)
is retained as an archive of the independent As schedule.

For an As-only Z-noise scan, the single-series API avoids sampling the other
panels:

```julia
as_plus_z = scan_channel_logical_failure(MersenneTwister(1235);
    construction=:as, logical_state=:plus, error_channel=:z_only,
    distances=[9, 11, 13, 15], error_rates=0:0.002:0.1,
    shots=50_000, batch_size=10_000, seed=1235)
```

The command writes one combined 2×2 figure and four standalone figures
(`As`/`Bp` × X-only/Z-only) in SVG, PDF, and PNG, along with raw and fit CSV
files. When a fit is available, each panel includes an inset with its
exploratory `p_c` and `nu` collapse. The `d=3,5,7` fit values are exploratory
diagnostics, not threshold claims. Use `--no-fit` to preserve the raw curves
and write explicit unavailable-fit records without attempting the fit.

The comparison API excludes repeated syndrome rounds, measurement noise,
periodic layouts, holes, and multi-logical-qubit patches. Its three-distance
scaling fits are exploratory diagnostics, not threshold claims.
