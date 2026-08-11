# Open-Boundary Code-Capacity Threshold Design

## Goal

Add a separate, conventional code-capacity experiment that estimates the decoding threshold of an open planar code on the repository's rectangular bond lattice. It retains the open north, south, east, and west boundaries used by the current sequential circuit; it does not replace that geometry with a periodic torus.

The existing trajectory scan remains a measurement-record and virtual-bond percolation diagnostic. It is deliberately unchanged. Because that experiment's physical record reconstructs every internal virtual error exactly, its current logical-failure value is not a code-capacity threshold observable.

## Model

The new experiment represents data-edge Pauli errors independently of the PEPS trajectory.

- `OpenCodeCapacityModel(rows, cols)` uses the same internal horizontal-edge shape `rows × (cols - 1)` and vertical-edge shape `(rows - 1) × cols` as the existing bond lattice.
- `DataEdgeErrors` stores only those two Boolean edge arrays. It has no measurement record and no virtual-boundary error fields.
- For each shot, each data edge receives an independent `X` error at probability `p`. The simulator retains that error only for scoring and never passes it to the decoder.
- Plaquette-parity syndrome is calculated from the data-edge errors using the existing incidence convention. The syndrome, rather than the error array, is the decoder input.

The matching graph is a planar, open-boundary graph derived from that same incidence convention. Its boundary assignment and two logical-cut masks will be tested against the stabilizer/check matrix before a threshold scan is trusted. The primary sector has a north--south logical direction; the east--west sector is its transposed cross-check. Both use the same open patch, not a periodic identification of edges.

## Decoder and logical verdict

Expose a syndrome-only decoder interface:

```julia
decode_syndrome(model, syndrome; sector=:north_south) -> Correction
```

`Correction` has horizontal and vertical Boolean arrays matching `DataEdgeErrors`. PyMatching may be used internally, but the decoder receives only the binary syndrome and known error-rate/edge weights. It cannot inspect the sampled data-edge error configuration.

For scoring only, compute `error XOR correction`. A logical failure is the parity of that residual chain through an explicit transverse logical cut. The north--south and east--west masks are geometry objects rather than implicit index choices. A zero-syndrome chain connecting the relevant boundaries changes exactly one declared logical parity.

Existing trajectory helpers may share low-level plaquette-incidence or matching construction code when their geometry truly agrees. The new public code-capacity API must not accept `ToricCodeTrajectory`, a mismatch record, or `VirtualBondErrors`; that separation prevents giving the decoder information a syndrome measurement would not reveal.

## Scan and estimate

Add a streaming scan API:

```julia
scan_open_code_capacity(rng; sizes, error_rates, shots, batches)
```

Each point records logical failures and rates in both sectors, their standard errors, and batch-level values. The scanner prints progress after each completed `(size, p)` point so that a long command is visibly active.

Adjacent-size crossings of north--south logical-failure curves are the primary finite-size estimate of `p_c`. Reuse the monotone interpolation and batch-bootstrap machinery after making it generic over the new scan point type. A crossing is reported only when the curves bracket a sign change inside the supplied error range; otherwise it is explicitly `unbracketed`, never coerced to an endpoint.

East--west curves and crossings are a geometry-symmetry check. Agreement within finite-size uncertainty supports the boundary and logical-mask implementation. Disagreement is a validation failure to investigate, not a second threshold to average blindly.

## Command-line example and output

Add `examples/scan_open_code_capacity.jl` with the established options:

```bash
julia --project=. examples/scan_open_code_capacity.jl \
  --p-min 0.05 --p-max 0.16 --p-step 0.005 \
  --sizes 8,12,16,20 --shots 10000 --batches 100
```

It writes separate scan-point and north--south/east--west crossing CSV files under `examples/output/`, without overwriting `trajectory_scan.{csv,pdf}`. Its PDF includes the two logical-failure-versus-`p` panels and crossing summaries, labels the primary north--south estimate, and marks unbracketed intervals clearly.

The example prints output paths, incremental progress, and a final compact summary. A small run is only a smoke test. The documented production command uses enough sizes, a grid focused on crossings, and enough shots for a defensible finite-size estimate. It reports an estimate and uncertainty rather than a universal exact critical point inferred from small lattices.

## Validation

Tests establish the decoding problem before interpreting a scan:

1. Zero noise has all-zero syndrome and zero logical-failure rate.
2. Every single internal data-edge error yields the expected syndrome and is corrected without logical failure in both sectors.
3. Each returned correction has the input syndrome, making the residual syndrome zero.
4. Known boundary-to-boundary zero-syndrome chains toggle the intended logical parity, but not the transverse one.
5. Rotation of a square error configuration exchanges north--south and east--west verdicts.
6. A seeded small scan produces valid point counts, finite standard errors, and correctly labels both bracketed and unbracketed crossings.

Run the full Julia test suite, a small deterministic example scan, and `git diff --check`. Only after those checks pass should a larger scan be used to quote a critical-point estimate.

## Scope boundaries

This work does not change the PEPS local tensor, the sequential Yao circuit, the exact trajectory sampler, or existing virtual-bond plots. It also does not simulate an unmeasured many-qubit state: code-capacity noise and syndrome sampling are classical/stabilizer-compatible and scale with the matching problem rather than a state vector of the whole lattice.
