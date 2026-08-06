```@meta
CurrentModule = TopoNoise
```

# Noisy toric-circuit trajectories

## Virtual bonds and timing

For an ``R\times C`` open patch, [`VirtualBondErrors`](@ref) stores one
independent Bernoulli-``p`` event on each of

```math
R(C-1)+(R-1)C+2R+2C
```

virtual bonds. Its arrays have the following locations:

- `horizontal_internal[r,c]` joins ``(r,c)`` to ``(r,c+1)`` and has shape
  ``R\times(C-1)``;
- `vertical_internal[r,c]` joins ``(r+1,c)`` to ``(r,c)`` and has shape
  ``(R-1)\times C``;
- `west_boundary[r]` and `east_boundary[r]` are the dangling ends of the
  horizontal carrier in row `r`;
- `south_boundary[c]` and `north_boundary[c]` are the dangling ends of the
  vertical carrier in column `c`.

Horizontal carriers travel west to east, while vertical carriers travel south
to north. West and south errors are applied after carrier preparation in
``|+\rangle``. At an internal bond, the error is applied after the upstream
site's measurement/reset and before the downstream site's gate. East and
north errors are applied immediately after the last site on their carrier and
before the final analytic ``\langle+|`` effect.

The four physical ancillas enter every local gate in ``|0000\rangle``. At site
``(r,c)``, [`sample_yao_trajectory`](@ref) applies the six-qubit gate, measures
the physical wires in the computational basis in `(E,N,W,S)` order, records
``(i,j,k,l)``, and resets those four wires to zero before reusing them.

Boundary channels are explicit sampled metadata. They do not change the
physical-record distribution or the final projection probability because

```math
X|+\rangle=|+\rangle,\qquad \langle+|X=\langle+|.
```

The final effects are handled analytically rather than by rejecting shots.
Their success probability is ``2^{-(R+C)}``, available from
[`postselection_probability`](@ref).

## Exact scalable sampler

[`sample_trajectory`](@ref) generates the same physical Born distribution
without a state vector. If the incoming horizontal and vertical carrier bits
at a site are ``(\gamma,\delta)``, it samples ``\alpha`` uniformly and sets

```math
\beta=\alpha\oplus\gamma\oplus\delta.
```

It then records
``(i,j,k,l)=(\alpha,\beta,\gamma,\delta)``. This is precisely the even-parity
support of the local toric-code isometry. Internal errors flip the outgoing
carrier after this record and are therefore visible at the next endpoint.
Boundary bits remain attached to the returned trajectory but need not alter
the Boolean carrier state.

For the doubled physical records, [`bond_mismatches`](@ref) computes

```math
m^h_{r,c}=i_{r,c}\oplus k_{r,c+1},\qquad
m^v_{r,c}=j_{r+1,c}\oplus l_{r,c}.
```

With the array orientation above, these maps equal `horizontal_internal` and
`vertical_internal` exactly. A dangling error has no doubled endpoint and
cannot be inferred from a mismatch measurement.

## Raw observables

[`trajectory_observables`](@ref) returns the total sampled error density,
boundary-error density, observable internal mismatch density, odd-plaquette
frustration density, largest occupied-internal-bond cluster fraction, and
horizontal/vertical spanning flags. Clusters and spans are computed by
union--find.

## Equal-weight marginalized spins

[`marginal_spin_observables`](@ref) deliberately treats every mismatched
internal doubled edge (`01` or `10`) as an erasure, even though its XOR is
observable. A matched edge (`00` or `11`) imposes an equal-spin constraint.
Union--find partitions the ``N=RC`` sites into matched components of sizes
``n_a``. Open boundaries leave each component with an independent, equally
weighted sign ``\eta_a\in\{-1,+1\}``, so

```math
M=\frac{1}{N}\sum_a n_a\eta_a.
```

The second and fourth moments are marginalized exactly:

```math
\mathbb E[M^2\mid t]=\frac{\sum_a n_a^2}{N^2},\qquad
\mathbb E[M^4\mid t]=
\frac{3(\sum_a n_a^2)^2-2\sum_a n_a^4}{N^4}.
```

The conditional absolute magnetization is estimated from `spin_samples`
independent component-sign draws. Equal weights make an erased edge impose no
spin constraint, so this is a bond-diluted component-spin model. It is not a
decoder and not a finite-temperature random-bond Ising model. Dangling virtual
legs do not join two sites and are excluded from the spin graph.

## Critical scans

[`scan_trajectories`](@ref) streams these quantities with online Welford
accumulators. From the same physical shots it also streams marginalized
``|M|``, ``M^2``, and ``M^4`` while using an independent random stream for
component signs. The scan forms the Binder cumulant only after averaging the
two exact conditional moments,

```math
U_4=1-\frac{\langle M^4\rangle}{3\langle M^2\rangle^2},
```

and estimates its standard error with a weighted delete-one-batch jackknife.

The scan retains horizontal-spanning and paired ``M^2``/``M^4`` batch means.
[`estimate_crossings`](@ref) monotone-smooths each spanning curve and estimates
adjacent-size crossings with batch bootstrap confidence intervals. A crossing
that is not bracketed or has too few valid bootstrap replicates is preserved
with status `:unbracketed` or `:unstable`.

[`estimate_binder_crossings`](@ref) recomputes ``U_4`` inside every paired
batch-bootstrap replicate and smooths Binder curves in the non-increasing
direction. It accepts one unique interior crossing; missing or ambiguous
crossings retain the same public statuses.

```julia
using TopoNoise, Random

rng = MersenneTwister(1234)
scan = scan_trajectories(
    rng, [4, 8, 16, 32], 0.35:0.01:0.65;
    shots=10_000, batches=100, spin_samples=1)
crossings = estimate_crossings(rng, scan; bootstrap=2_000, confidence=0.95)
binder_crossings = estimate_binder_crossings(
    rng, scan; bootstrap=2_000, confidence=0.95)
figure = plot_trajectory_scan(
    scan, crossings; binder_crossings=binder_crossings)
```

The runnable script requires the three error-grid arguments:

```bash
julia --project=. examples/scan_toric_trajectories.jl \
  --p-min 0.35 --p-max 0.65 --p-step 0.01 --spin-samples 1
```

It produces an aggregate CSV, separate spanning- and Binder-crossing CSVs, and
six-panel SVG, PDF, and PNG figures. The original four raw-percolation panels
remain alongside marginalized ``|M|`` and ``U_4``. The plotted frustration
reference is ``\tfrac12[1-(1-2p)^4]``.
