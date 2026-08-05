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

## Observables and critical scan

[`trajectory_observables`](@ref) returns the total sampled error density,
boundary-error density, observable internal mismatch density, odd-plaquette
frustration density, largest occupied-internal-bond cluster fraction, and
horizontal/vertical spanning flags. Clusters and spans are computed by
union--find.

[`scan_trajectories`](@ref) streams these quantities with online Welford
accumulators and retains horizontal-spanning batch means. Then
[`estimate_crossings`](@ref) monotone-smooths each spanning curve and estimates
adjacent-size crossings with batch bootstrap confidence intervals. A crossing
that is not bracketed or has too few valid bootstrap replicates is preserved
with status `:unbracketed` or `:unstable`.

```julia
using TopoNoise, Random

rng = MersenneTwister(1234)
scan = scan_trajectories(
    rng, [4, 8, 16, 32], 0.35:0.01:0.65;
    shots=10_000, batches=100)
crossings = estimate_crossings(rng, scan; bootstrap=2_000, confidence=0.95)
figure = plot_trajectory_scan(scan, crossings)
```

The runnable script requires the three error-grid arguments:

```bash
julia --project=. examples/scan_toric_trajectories.jl \
  --p-min 0.35 --p-max 0.65 --p-step 0.01
```

It produces aggregate and crossing CSV files plus SVG, PDF, and PNG figures.
The plotted frustration reference is
``\tfrac12[1-(1-2p)^4]``.

These observables describe raw percolation of the sampled internal error
bonds. They are not a spin/Binder analysis and not a Nishimori decoding
threshold. Either interpretation would require an additional reconstruction,
Gibbs distribution, or decoder that is not defined here.
