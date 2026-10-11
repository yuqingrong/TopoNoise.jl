# Encoded diagonal-sweep logical-X feedback study

## Status

Design approved in conversation on 2026-10-11. This specification defines a
new logical-X preparation experiment for the existing diagonal-slice 3D toric
code. It preserves the current gate-output noise convention and keeps the
earlier rough-slab, physical slice-diagnostic, and logical-Z studies intact.

## Scientific question

Does ideal stabilizer measurement and classical feedback within every
diagonal slice give the encoded diagonal-sweep circuit a nonzero logical-X
preparation threshold under independent post-gate-output X faults?

The experiment compares two protocols with identical sampled fault records:

1. **Final-only:** run the diagonal preparation without intermediate
   corrections, then decode the complete final plaquette syndrome.
2. **Slice feedback:** before every transfer, measure the active-slice
   compatibility checks, use minimum-weight perfect matching (MWPM) to pair
   their point defects, move each endpoint one graph edge toward its partner,
   propagate that virtual repair to the fresh targets, and finally decode the
   complete plaquette syndrome.

The primary result is a logical-X state-failure probability and, if the data
support it, a finite-size threshold estimate. Physical syndrome and cluster
statistics may be retained as diagnostics, but they cannot replace logical
failure in the threshold analysis.

## Relationship to the paper and earlier experiments

The measurement-feedback prescription follows the qualitative construction
in *Noise-resilient sequential circuits for generating quantum order*: point
defects are measured in a two-dimensional slice, MWPM pairs them, and local
feedback moves each defect toward its partner. The explicit diagonal lattice,
seven-gate transfer, boundary preparation, routing tie breaks, gate-output
noise, and final decoder are the implementation defined in this repository.
The result is therefore a test of this derived circuit under the repository's
noise model, not a reproduction of a numerical threshold reported by the
paper.

The existing `Toric3DSliceLattice` experiment has fixed product boundaries and
encodes no logical qubits. It remains the propagation diagnostic. The existing
`Toric3DEncodedSlice` logical-Z experiment prepares logical `|++>` and uses no
active feedback under pure Z noise. The new experiment reuses its encoded
geometry but prepares logical `|00>` so residual logical X is an actual state
preparation error.

## Encoded geometry and logical state

For transfer depth `T` and periodic transverse dimensions `Ly,Lz`, use the
existing diagonal lattice with three oriented edge qubits per cell and slices
`0:T`. Its parameters are

\[
n=3(T+1)L_yL_z,\qquad k=2,\qquad
d_X=(T+1)\min(L_y,L_z),\qquad d_Z=2\min(L_y,L_z).
\]

The two logical X operators are transverse membranes through all slices. The
dual logical Z operators are the two independent honeycomb cycles on the
initial slice. Their binary symplectic pairing is the identity matrix.

Prepare the seed boundary from product `|0>` qubits. Product `|0>` already has
eigenvalue `+1` for every seed Z compatibility check and both logical Z
cycles. Ideally measure a complete independent set of seed X stars and apply
a Z correction with the measured star syndrome. The correction commutes with
the logical Z cycles, so the resulting code state is logical `|00>`. The
simulation need not sample the random projection outcomes: ideal measurement
followed by a syndrome-valid ideal correction is exactly represented by the
zero Pauli frame on this prepared seed.

This seed preparation introduces no noisy H/CNOT location. That is a direct
consequence of the agreed model in which initialization, stabilizer
measurements, and corrections are ideal. The metadata must state this
asymmetry with the logical-Z study, whose logical `|++>` seed uses noisy H
outputs.

## Transfer circuit and noise

Use the unchanged diagonal seven-gate transfer at each `(t,u,v)` cell:

1. `H(new1)`;
2. `CNOT(new1,new2)`;
3. `CNOT(new1,new3)`;
4. two old-slice parity controls into `new2`;
5. two old-slice parity controls into `new3`.

There are thirteen participating gate outputs per cell. Independently apply an
X fault with probability `p` to each participating qubit immediately after
each H or CNOT. Initialization, idles, compatibility measurements, seed-star
measurements, classical processing, and Pauli corrections remain ideal. Pure
X output faults stay in the X Pauli frame under the remaining CNOT gates.

Fault records are sampled without access to either decoder and replayed under
both protocols. The supplied RNG controls the samples; an optional integer
seed is recorded only as provenance.

## Active-slice measurement feedback

At the beginning of transfer `t`, measure all six-edge Z compatibility checks
`Q` on slice `t-1`. A physical old-slice edge joins exactly two checks, so the
syndrome defines a periodic triangular matching graph. The decoder receives
only the syndrome and the declared unit-weight graph.

MWPM must return a complete old-slice recovery whose syndrome exactly equals
the measurement. Invalid recoveries stop the scan. For the primary
`one_step` policy, choose deterministic shortest paths for the matched pairs
and toggle the first edge and, for paths longer than one edge, the last edge.
This moves both endpoints one graph step toward each other; overlapping moves
combine modulo two. Map this virtual old-slice route to its propagated target
image `F(C)` and apply that X correction only to the fresh target slice after
the fixed noisy gates. This is the existing frame implementation of ideal
measurement feed-forward and leaves future fault locations unchanged.

The final transfer has no later active slice on which its new defects can be
fed back. Those faults, together with any residual earlier faults, are handled
by the common final decoder. A `full_path` policy that applies the entire
matching route is a labeled sensitivity study because it uses stronger,
spatially nonlocal feedback than the paper-like one-step policy.

## Final X recovery and logical scoring

Measure the complete physical Z-plaquette syndrome `H_Z e_X` after
preparation. Build a BP+OSD decoder from the encoded slab's full plaquette
matrix. The initial baseline uses uniform physical-`p` data priors and the
same fixed settings as the rough-slab X decoder:

- minimum-sum BP;
- 50 iterations;
- scaling factor 1;
- parallel update schedule;
- one thread; and
- `OSD_CS` order 2.

The decoder input contains only the syndrome, geometry, declared `p`, and
these settings. It cannot inspect the sampled circuit faults or intermediate
frames. It must return a full data-qubit X recovery. Recompute `H_Z r_X` and
stop with the dimensions, probability, protocol, shot, and seed if it differs
from the input syndrome.

After applying the validated recovery, score the residual frame against both
dual logical Z cycles. Record logical-1 failure, logical-2 failure, and either
logical failure. Because the seed is logical `|00>`, either residual logical X
is exactly a preparation failure rather than only a Pauli-frame diagnostic.

Uniform data priors neglect nonuniform output exposure, CNOT-induced
correlations, and correlations introduced by syndrome-dependent feedback.
The reported threshold is consequently decoder-dependent. A circuit-level
maximum-likelihood decoder is outside this study.

## Public API and artifacts

Add a separate X-feedback threshold subsystem without changing existing 2D or
3D entry points. Its public API will provide:

- a logical-`|00>` encoded-boundary verifier;
- a BP+OSD full-X-recovery decoder;
- paired final-only and feedback shot/count functions;
- a configurable scan accepting independent `(T,Ly,Lz)` dimensions,
  probabilities, shots, batch/checkpoint settings, RNG, seed metadata,
  decoder settings, routing, output directory, and progress stream; and
- save and plot functions.

Use descriptive X-specific types and names, centered on
`compare_toric3d_x_feedback_threshold`, so this experiment cannot be confused
with the Z-only threshold or fixed-boundary physical-diagnostic scans.

Write completed points atomically after every `(dimensions,p)` pair. The raw
CSV contains both protocol counts, paired discordant counts, individual
logical counts, rates, intervals, elapsed time, routing, and decoder version.
TOML metadata records every geometry, preparation, circuit, noise, feedback,
decoder, RNG, statistical, and interpretation convention. Deliver PDF, SVG,
and PNG figures plus a runnable example and a Python analysis program.

## Staged numerical study

First time a compiled warmup and a short pilot before committing to the
production grid. The broad pilot uses cubic `L=T=Ly=Lz` with

- `L = 3,5,7,9`;
- `p = 0, 1e-4, 3e-4, 1e-3, 2e-3, 3e-3, 5e-3, 7.5e-3,
  1e-2, 1.5e-2, 2e-2`;
- 1,000 paired shots per point;
- one-step feedback; and
- seed 1234.

If the pilot brackets a common crossing, refine around it with at least five
sizes, at least seven probabilities in the critical window, and 10,000 paired
shots per point. Include larger sizes until either the crossing drift is
resolved or runtime/memory scaling makes the limitation explicit. If the
pilot does not bracket a crossing, extend the probability range in the
indicated direction before attempting a fit. Do not force a threshold fit to
data without crossings.

Run full-path feedback on the final critical window as a labeled sensitivity
check. A short pilot determines whether 10,000 shots at every full-path point
is practical; otherwise use enough independent paired shots to make its Monte
Carlo uncertainty smaller than the observed difference from one-step
feedback. Never pool timing pilots into production counts.

## Statistical analysis

For each protocol and logical outcome, report 95% Wilson score intervals.
Use paired discordant counts to report feedback-minus-final-only differences
and uncertainty. Plot both protocols together for every size, and plot raw
adjacent-size crossings and the full-path sensitivity separately.

Only estimate a threshold when the simulated window brackets crossings and
the fitted model passes numerical guards. Use a binomial-likelihood
finite-size scaling model with a quadratic scaling function. Report:

- the fitted threshold and conditional parametric-bootstrap interval;
- fit deviance and degrees of freedom;
- sensitivity to the minimum included size and probability window;
- sensitivity to a leading finite-size correction when the raw crossings
  drift systematically; and
- the spread across acceptable models separately from sampling uncertainty.

The fit treats protocol counts separately; paired statistics quantify the
direct feedback benefit at each simulated point. Do not claim that a fitted
critical exponent is universal from this finite-size study.

## Verification

Automated tests must establish:

1. stabilizer commutation, ranks, `k=2`, logical anticommutation, and the
   stated distances, including periodic transverse length two;
2. logical-`|00>` seed preparation from product zero plus ideal X-star
   measurement and syndrome-valid Z correction;
3. noiseless preparation and every checkpoint stabilizer, including both
   logical Z eigenvalues;
4. fresh targets and equivalence of physical old-slice correction with the
   propagated target-frame correction;
5. every eligible single X gate-output fault on a small instance under both
   protocols and both routing policies;
6. active-slice MWPM and final BP+OSD recovery validity;
7. correct scoring of injected stabilizers and both logical X membranes;
8. zero-noise success, paired replay, seeded reproducibility, and incremental
   checkpoint recovery;
9. statistical analysis on synthetic data, including rejection of missing
   crossings and parameter-bound fits; and
10. the complete existing 2D and 3D regression suite.

## Interpretation and limitations

A nonzero crossing would support a decoder-dependent logical-X threshold for
this encoded diagonal circuit under the specified ideal-measurement noise
model. It would not by itself reproduce the paper's mixed-state phase
criterion or prove a local-channel theorem. The active-slice quantum
corrections are local, but global MWPM classical processing is not a
constant-depth spatially local operation. The study excludes measurement
faults, correction faults, idle faults, combined X/Z noise, and a
circuit-informed final decoder.

If logical failure decreases with size throughout the scan but crossings are
not resolved, report a lower bound on the threshold window rather than a point
estimate. If it increases throughout, extend toward smaller `p` and test for a
zero-threshold pseudothreshold before drawing a conclusion.

## Acceptance criteria

The work is complete when the new scan is reproducible from a documented CLI,
all recoveries are validated, the staged production data are saved
incrementally, the statistical conclusion follows the guards above, the note
and experiment documentation distinguish state failure from diagnostics, all
figures and machine-readable artifacts are present, and the full regression
suite passes.
