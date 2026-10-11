# Logical-Z threshold study under gate-output noise

The user requests an actual threshold simulation after the exploratory slice
experiment. Retain independent rate-p Z faults after each active H/CNOT output,
ideal measurements, corrections and idles. This extends the authorized study;
no threshold value will be assumed in advance.

## Scientific definition

Use the diagonal transfer circuit with a natural honeycomb toric-code seed.
There are two transverse logical qubits, n=3(T+1)LyLz edges,
rank(Hx)=(T+2)LyLz−1, rank(Hz)=(2T+1)LyLz−1, and
dZ=2min(Ly,Lz). Include lower white triad stars, six-edge bulk stars and upper
triad stars. The existing physical plaquettes are the Z checks. Logical-X
membranes are direction-1 edges at u=0 and direction-2 edges at v=0 through
all slices. Logical-Z representatives are alternating 1/3 or 2/3 seed loops.

Prepare the seed with H on every initial qubit, each followed by Z-p noise,
then ideal Q measurements and X recovery. This prepares logical |++>; the
noise commutes through those measurements, so its exact frame is iid seed Z.
Transfer gates retain the previous schedule and noise. Slice feedback is
inactive under pure Z noise. Final unit-weight data-edge MWPM receives only
the full star syndrome; every returned data recovery is checked. Failure
means either residual winding bit is nonzero. This is a Z-sector threshold
for the specified preparation/recovery procedure, not a two-sector result.

## Tasks and ownership

- [x] Geometry agent: encoded geometry, checks/logicals, signed seed/prefix
  verification, growing-distance argument, length-two and rectangular tests.
- [x] Root: exact sparse map from independent seed/gate faults to final data
  Z frames; batch sampling and syndrome-only MWPM with full recovery validation;
  public scan API and CLI, incremental CSV/TOML and timing pilot.
- [x] Fitting agent: standalone binomial-likelihood finite-size fitter,
  adjacent-size crossing bootstrap, fit diagnostics, Wilson plots and collapse,
  bootstrap intervals and synthetic/rejection tests.
- [x] Root: pilot sizes 5,7,9, probabilities .002,.004,.006,.008,.010,.015,.020,
  1000 shots. Refine using larger sizes and at least 5000 shots per point in
  the crossing neighborhood; save every point and record seeds/settings.
- [x] Root/reviewer: exhaustive fault-map agreement with explicit propagation,
  graphlike fault syndromes, all eligible single faults, noiseless behavior,
  seeded reproducibility, invalid recovery rejection, old regression suite.
- [x] Root: assess fit window and smallest-size sensitivity; report a rounded
  threshold with statistical and finite-size limitations, or an honest bound
  if the data do not support a crossing. Deliver CSV, metadata, fit JSON,
  PDF/SVG/PNG, runnable examples, and scientific documentation.

## Implementation contracts

`Toric3DEncodedSlice(steps,Ly,Lz)` wraps a `Toric3DSliceLattice` as `lattice`.
Standard geometry/check-matrix functions dispatch on the wrapper.
`toric3d_encoded_logical_x/z` return two supports. Seed checks return `(x,z)`.
`toric3d_encoded_encoder` returns the existing transfer encoder; noisy seed H
faults are modeled separately. `verify_toric3d_encoded_encoder` verifies the
ideal |++> state and every prefix.

The exact fault map has one column per seed H output and each transfer gate
output, with at most three rows per column. Construct it using reverse suffix
propagation; verify independently against explicit gate replay. Batch sampling
draws every fault independently and multiplies this map over GF(2). This is an
acceleration of the declared gate-level model, not an iid data-error substitute.

Raw CSV columns: steps,Ly,Lz,p,shots,failures,logical_1_failures,
logical_2_failures,seed,decoder,elapsed_seconds. Metadata records gate/fault
counts, boundary and seed preparation, noise, decoder version, RNG and exact
sampling method. Fit the probability of either logical Z with a local logistic
quadratic finite-size ansatz, binomial likelihood and parametric bootstrap;
compare exclusions of small sizes and narrower fitting windows.

No unrelated changes, commits, or alteration of earlier result sets.

## Completed result

The 21,000-shot pilot was followed by 330,000 production shots at
L=9,13,17,21,25 (33 points, 10,000 shots each). The corrected finite-size
estimate is p_c,Z=.00517145, approximately 0.52% per gate-output qubit.
The conditional 95% interval from 999 accepted fits out of 1,000 bootstrap
replicas is [.00504584,.00529814]. Model/window/size choices give estimates
roughly 0.49%–0.54%; this is a sensitivity range, not a confidence interval.

The full Julia regression suite passed, including the existing 2D tests and
the new 135 encoded-geometry and 376 fault/decoder/scan assertions. Twelve
Python statistical tests passed. Independent reviews checked geometry,
fault modeling, provenance, numerical gradients and the fitted estimate.
Both baseline and corrected fit JSON/CSV and PDF/SVG/PNG artifacts are
saved under results/toric3d-z-threshold-2026-10-11. The note was updated
and compiled successfully; pre-existing unrelated warnings remain.
