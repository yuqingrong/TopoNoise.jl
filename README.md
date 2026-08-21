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

## Noisy measured trajectories

`ToricCodeTrajectoryModel` implements the measurement/reset version of the
same diagonal circuit. It uses four reusable physical ancillas and one carrier
per row and column. At every site the physical record is stored in
`(E, N, W, S)` order, the four ancillas are reset to zero, and any internal
carrier error is applied before that carrier reaches its next site.

Every internal and dangling virtual bond receives one independent Bernoulli
`X`-error opportunity. For an `R` by `C` patch there are

```text
R(C-1) + (R-1)C + 2R + 2C
```

such opportunities. The west/south boundary errors occur after the carrier
`|+⟩` preparation; east/north errors occur after the final site on the
corresponding carrier and before its analytic `⟨+|` projection. These boundary
events are deliberately retained in `VirtualBondErrors`, even though they do
not affect physical records because `X|+⟩=|+⟩` and `⟨+|X=⟨+|`.

```julia
using TopoNoise, Random

model = ToricCodeTrajectoryModel(8, 8)
trajectory = sample_trajectory(
    MersenneTwister(1234), model; error_rate=0.1)
mismatches = bond_mismatches(trajectory)
observables = trajectory_observables(trajectory)

# Literal small-system Yao reference: 4 reusable ancillas + R+C carriers.
errors = sample_virtual_errors(MersenneTwister(7), model; error_rate=0.1)
reference = sample_yao_trajectory(
    MersenneTwister(8), model, errors; max_qubits=24)
```

The observable internal doubled-edge syndromes are
`E[r,c] ⊻ W[r,c+1]` horizontally and
`N[r+1,c] ⊻ S[r,c]` vertically. They reproduce the internal error maps
exactly. Dangling boundary errors have no doubled neighbor and therefore
cannot be inferred from these mismatch records.

### Finite-size raw-percolation scan

The trajectory scanner reports injected-error density, boundary density,
internal mismatch density, odd-plaquette frustration, largest occupied-bond
cluster, and horizontal/vertical spanning. Run the example with an explicit
error-rate grid:

```bash
julia --project=. examples/scan_toric_trajectories.jl \
  --p-min 0.35 --p-max 0.65 --p-step 0.01
```

Defaults are sizes `4,8,16,32`, 10,000 shots per point, seed `1234`, and
2,000 bootstrap replicates. The script writes `trajectory_scan.csv`,
`trajectory_crossings.csv`, and `trajectory_scan.{svg,pdf,png}`. Adjacent-size
crossings use monotone-smoothed horizontal-spanning curves and batch bootstrap
intervals; unbracketed or unstable results remain in the output with an
explicit status.

This is a raw bond-percolation diagnostic of the sampled internal errors. It
does not reconstruct spins, compute a Binder cumulant, or solve a Nishimori
decoding problem; those require an additional Gibbs model or decoder.

## Isometric planar-code virtual-bond capacity

`IsometricPlanarCodeModel(d)` is an open planar logical-code experiment whose
public size is the direct distance `d` (not a plaquette count). It prepares
`|0_L⟩` with west/east `|+⟩` and south/north `|0⟩` virtual boundaries. The
noise model is independent `X` noise on internal virtual carriers only;
preparation, local gates, and parity checks are perfect.

The decoder receives `(d-1)^2` plaquette checks plus one **local** north and
south check for each horizontal boundary bond. The held-out logical frame is
the parity of residual west-boundary vertical bonds. This makes a west-to-east
virtual string the shortest odd logical path, with weight exactly `d`.

```julia
using TopoNoise, Random

model = IsometricPlanarCodeModel(5)
point = estimate_isometric_planar_capacity(
    MersenneTwister(1234), model, 0.10; shots=10_000, batches=100)
```

`d=2` is supported as a diagnostic geometry with a fixed tie rule, but it is
excluded from threshold fitting. For a finite-size-scaling fit, use at least
three non-diagnostic distances; a useful production choice is
`d=7,9,11,13,15`:

```bash
julia --project=. examples/scan_isometric_planar_threshold.jl \
  --distances 7,9,11,13,15 \
  --output-dir results/isometric-planar-code/d7-15
```

The command writes CSV, SVG, PNG, and PDF artifacts under
`results/isometric-planar-code/`. Its defaults are `p=0:0.005:0.16`,
100,000 shots per point, 100 batches, 2,000 bootstraps, and seed `1234`.
Crossings with no unique interior intersection are reported as unstable.
The publication figure marks the fitted `p_c` with a gray dashed line and
adds a collapse inset. It fits a cubic master curve to points with logical failure rates from
0.05 to 0.45 under
`P_fail = F((p-p_c)d^(1/nu))`; the inset reports bootstrap one-sigma values
for `p_c` and `nu`. The accompanying
`isometric_planar_capacity_scaling_fit.csv` records the fit, fit window, and
bootstrap count. Two-size smoke scans still render and record an
`insufficient_sizes` fit status without an inset.

For a reproducible fit audit, use
`diagnose_isometric_planar_scaling(rng, scan)`. It returns the ordinary fit
together with the bootstrap refits, the profiled `(p_c, nu)` loss surface,
the running-best optimizer path, empirical batch variation, and nominal
leave-one-distance-out and fit-window sensitivity results. The scan example
writes the raw batch rates and each diagnostic table to CSV and produces a
separate `isometric_planar_capacity_scaling_diagnostics.pdf`; the publication
capacity figure remains uncluttered. A zero parenthetical bootstrap error is
treated as a diagnostic condition and should be checked against the batch,
optimizer, and loss-surface panels rather than interpreted as exact physical
precision.

`isometric_planar_encoder` retains terminal carriers as physical outputs and
therefore has `4d²` output wires. `sample_isometric_planar_yao_trajectory` is
kept only as a small-distance detector reference; the production estimator
calculates the same detector parities algebraically and decodes batches with
PyMatching rather than allocating an exponential Yao state vector.
