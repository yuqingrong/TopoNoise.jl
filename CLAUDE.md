# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Despite the directory name `admin-dashboard`, this repository is **TopoNoise.jl**, a Julia
package for simulating noisy measured trajectories of the doubled-edge toric-code PEPS.

## Commands

All commands assume the repo root and use the local project environment.

```bash
# One-time (or after Project.toml/Manifest.toml changes)
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Full test suite
julia --project=. -e 'using Pkg; Pkg.test()'

# Single test file (mirrors what runtests.jl includes)
julia --project=. -e 'using TopoNoise, ITensors, Test; include("test/trajectory.jl")'

# Runnable examples
julia --project=. examples/generate_toric_circuit.jl 2 2
julia --project=. examples/scan_toric_trajectories.jl --p-min 0.35 --p-max 0.65 --p-step 0.01

# Build docs
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

CI (`.github/workflows/CI.yml`) runs `Pkg.test()` and builds docs; match locally before pushing.

## Architecture

The package builds three progressively richer objects on top of the same finite rectangular
lattice, and all of them share the `(east, north, west, south)` physical-index / bond
convention. Understanding that shared ordering is the key to reading any file here.

### Layered model stack (order matches `src/TopoNoise.jl` includes)

1. **`toric_code_peps.jl` — `ToricCodePEPS`.** The rank-8 parity local tensor and the
   finite doubled-edge PEPS on an open rectangular patch. Dangling virtual legs are
   terminated with `|+⟩ = (|0⟩+|1⟩)/√2` and the state is normalized. `physicalinds(peps, r, c)`
   returns the four physical qubit indices in `(E, N, W, S)` order — every downstream file
   assumes this order.
2. **`local_gate.jl` — `toric_code_local_gate`, `unitary_completion`.** The 64×64
   local unitary used by the sequential circuit. `unitary_completion` extends an isometry
   using `nullspace` from `LinearAlgebra`; changes here ripple into circuit correctness tests.
3. **`network_graph.jl` — `peps_graph`.** Converts a `ToricCodePEPS` to an
   `ITensorNetwork` keyed by `(row, col)`; shared virtual indices become graph edges,
   all four physical indices remain external.
4. **`sequential_circuit.jl` — `ToricCodeSequentialCircuit`.** A diagonal
   lower-left → upper-right preparation schedule. Uses one horizontal carrier per row and
   one vertical carrier per column. `yao_unitary(circuit)` is the reusable local-gate core
   (state never materialized); `yao_circuit(circuit)` prepends `|+⟩` carrier prep and returns
   an executable `Yao.ChainBlock`. Final `⟨+|` carrier projections are model metadata, not
   part of the executable block. `log2_postselection_probability` / `postselection_probability`
   report the analytic success weight.
5. **`trajectory.jl` — `ToricCodeTrajectoryModel`, `VirtualBondErrors`, `ToricCodeTrajectory`.**
   Measurement/reset version of the same diagonal circuit using four reusable physical
   ancillas plus row/column carriers. Every internal and dangling virtual bond gets one
   independent Bernoulli `X`-error opportunity — for an R×C patch:
   `R(C-1) + (R-1)C + 2R + 2C` opportunities total. Boundary errors are retained in
   `VirtualBondErrors` even though `X|+⟩=|+⟩` / `⟨+|X=⟨+|` make them invisible to physical
   records; downstream percolation/spin analysis relies on them being present.
   `sample_yao_trajectory` is a small-system Yao reference (guard with `max_qubits`).
6. **`trajectory_analysis.jl` — observables, scans, crossings.** Derives
   `TrajectoryObservables` (mismatch/frustration/cluster/spanning densities) and
   `MarginalSpinObservables` (equal-weight spin moments after marginalizing mismatched
   edges) from a trajectory. `scan_trajectories` produces `TrajectoryScan` sweeps over
   `(L, p)` grids with batch bootstrapping; `estimate_crossings` / `estimate_binder_crossings`
   report adjacent-size crossings with monotone-smoothed horizontal-spanning curves.
   Unbracketed/unstable results are kept in the output with an explicit status field.
7. **`visualization.jl` — CairoMakie plots.** `plot_peps_graph`, `plot_sequential_circuit`
   (with `expand_physical_buses` option), `plot_trajectory_scan`. Yao's native
   `Yao.vizcircuit` is preferred for raw circuit diagrams.

### Cross-cutting conventions

- **Direction/order everywhere is `(E, N, W, S)`** — for physical indices, per-site
  measurement records `(measurements[i, j, k])`, virtual bond orientation, and boundary
  error fields. If you add fields to `VirtualBondErrors` or `TrajectoryObservables`, update
  `length` / `count` / `_validate_virtual_error_shapes` in `trajectory.jl` and the scan
  aggregators in `trajectory_analysis.jl` in lockstep.
- **Observable syndromes** on internal doubled edges are
  `E[r,c] ⊻ W[r,c+1]` (horizontal) and `N[r+1,c] ⊻ S[r,c]` (vertical); these are the
  ground truth used by tests and analysis.
- **The scan pipeline is a raw bond-percolation diagnostic**, not a Nishimori decoder or
  Binder-with-Gibbs analysis. Don't conflate them when extending observables — the
  marginal-spin path is the separate `MarginalSpinObservables` branch.
- **Never allocate the global state vector.** The sequential circuit deliberately keeps only
  the reusable local-gate core in `yao_unitary`. The measurement/reset trajectory model
  uses four reusable ancillas plus carriers. Reference-comparison paths (`sample_yao_trajectory`)
  are opt-in and gated by `max_qubits`.

### Tests

`test/runtests.jl` includes one file per source module plus `example.jl` and
`trajectory_example.jl` (end-to-end smoke tests of the runnable examples). When adding a
source module, add a matching test file and include it in `runtests.jl`.
