# Rotated planar code-capacity threshold design

## Purpose

Replace the invalid open-boundary code-capacity implementation with a genuine
open rotated planar surface-code experiment.  The experiment estimates the
finite-size crossing of the logical failure rate under independent data-qubit
X noise and perfect Z-check syndrome extraction.  It does not use PEPS virtual
bonds as syndrome bits.

The existing virtual-trajectory workflow has a different purpose and remains
unchanged.

## Model

For each requested distance `d`, construct Stim's
`surface_code:rotated_memory_z` generated circuit with one stabilizer
measurement round and no built-in noise.  The generated patch is an open
rotated planar surface-code patch, not a periodic torus.

The circuit is modified in one place only:

1. Determine the data qubits from the final data-basis measurement instruction.
2. Insert independent `X_ERROR(p)` on precisely those data qubits immediately
   after the preparation block and before syndrome-extraction gates.
3. Leave all preparation, Clifford gates, ancillary qubits, and measurements
   noiseless.

The resulting detection events are the perfect syndrome information from the
Z checks affected by X errors.  The generated circuit may retain deterministic
detectors for the complementary checks; these are harmless zero components of
the syndrome and are retained so the detector model and the circuit have the
same indexing.

The model has one logical observable and one meaningful failure channel in
this noise setting.  A scan therefore reports `logical_x_failure_rate`, not
fabricated north-south/east-west channels.  A coordinate-rotated duplicate,
if later desired, is a validation of symmetry rather than a second physical
observable of one patch.

## Sampling and decoding

For every `(d, p)` point:

1. Obtain a graphlike Stim detector-error model (DEM) from the modified
   circuit.  The only nonzero error mechanisms are the inserted independent
   data-X errors.
2. Build `pymatching.Matching` from that DEM.
3. Sample detector events and logical-observable flips from Stim in one batch.
4. Give only detector events to `decode_batch`.
5. Compare decoder predictions with the separately sampled logical observable;
   their XOR is the logical-X decoder failure for each shot.

No sampled error configuration, virtual-bond value, or actual logical bit is
available to the decoder.  The actual logical bit is evaluation-only ground
truth.

## Public surface and migration

The replacement public types and functions use the `RotatedCodeCapacity`
prefix, including a point/result type, the circuit/model constructor, a
single-point estimator, a scan function, and a crossing estimator.  The
example is `examples/scan_rotated_code_capacity.jl` and writes
`rotated_code_capacity_*` artifacts.

The invalid `OpenCodeCapacity*` types, scanner, example, and output naming are
removed rather than retained as an apparently valid alternative.  Any
transition shim must fail explicitly and direct callers to the new API; it
must never produce an old-geometry estimate.

`src/visualization.jl` is also repaired so the legacy
`plot_trajectory_scan` reads only fields that are present in the committed
trajectory result type.  This is an isolation repair: no dirty trajectory
files are staged and the virtual trajectory workflow remains separate from
the code-capacity experiment.

## Inputs, statistics, and output

Distances are positive planar-code distances accepted by the generated circuit;
the example defaults to odd values such as `3,5,7`.  Physical probabilities
must lie in `[0, 1]`.  Scans require a nonempty increasing probability grid,
positive `shots` and `batches`, and `shots % batches == 0`.  Equal-size batches
make each bootstrap replicate an unambiguous resampling of the same estimator.

Each scan row records at least distance, `p`, shots, logical failure count,
logical failure rate, binomial standard error, batch count, and the seed when
one is supplied.  Batch-level failure counts are retained in memory for the
crossing bootstrap rather than exposed as misleading per-shot data.

Crossings are estimated from neighbouring finite-size curves and a batch
bootstrap.  The result has an explicit status:

- `ok`: a bracketed crossing and finite confidence interval;
- `unbracketed`: no sign change in the scanned probability range;
- `unstable`: bootstrap results do not support one finite interval.

The figure contains one panel of logical-X failure curves and one crossing
summary panel or inset.  It states the model as “open rotated planar code,
data-X noise, perfect Z syndrome” and visibly prints the primary crossing
estimate and interval, or the non-estimate status.  CSV, SVG, PDF, and PNG
outputs carry the same unambiguous model name.

## Dependencies

Add the maintained Python `stim` package to `CondaPkg.toml`; keep the existing
PyMatching dependency.  Julia continues to call both through PythonCall.  No
state-vector simulation or Yao trajectory sampling is used for this scan.

## Required validation

The implementation is accepted only when all of the following hold:

1. A noiseless generated circuit has deterministic zero detector events and
   zero logical-observable flips.
2. The data-qubit identification and insertion position are asserted from the
   generated circuit structure, so no ancilla receives code-capacity noise.
3. Every single injected data-X pattern at small distances is decoded without
   logical failure.
4. A representative minimum boundary-to-boundary X string of weight `d` flips
   the logical observable without a detection event, while local
   stabilizer-equivalent patterns do not.  Stim's graphlike shortest-logical
   diagnostic must also return length `d` for the corresponding circuit/DEM.
5. The decoder test proves logical-observable ground truth is withheld from
   the decoding input.
6. Invalid dimensions, probabilities, empty scans, and non-divisible
   shots/batches are rejected clearly.
7. Synthetic crossing data exercise `ok`, `unbracketed`, and `unstable`; a
   deterministic mixed-status figure test verifies a visible status marker.
8. The example smoke test verifies every promised CSV plus SVG, PDF, and PNG,
   including the crossing summary in the rendered figure.
9. A clean committed checkout passes the focused tests, the example smoke
   test, and the package test suite.  The legacy trajectory visualization is
   checked separately to prove this work did not alter it.

## Non-goals

This work does not claim a universal asymptotic threshold from a small scan,
does not introduce circuit-level gate/measurement noise, and does not reuse
virtual PEPS bonds as physical syndrome data.  A numerical critical-point
claim is reported only after a sufficiently resolved scan produces a stable
finite-size crossing and its bootstrap interval.
