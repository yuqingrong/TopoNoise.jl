# Diagonal 3D toric-code feedback-strength comparison

## Status and purpose

Design approved in conversation on 2026-10-12. This specification extends the
existing encoded diagonal-sweep simulations to compare three feedback scopes
under pure X noise and pure Z noise:

1. no intermediate feedback;
2. constant-width slice-local feedback; and
3. full-prefix feedback.

The purpose is to separate the effect of sweep geometry from the strength and
spatial support of measurement feedback. For every shot in a channel, the
simulator samples one immutable fault record, replays that same record through
all three protocols, and scores logical state-preparation failure after a
common final decoder.

The existing rough Cartesian slab, physical slice diagnostic, encoded Z-only
threshold scan, and encoded X one-step threshold scan remain available and
unchanged. The new subsystem provides a paired three-protocol comparison rather
than silently changing those published artifacts.

## Scientific questions

For the encoded diagonal circuit and retained post-gate-output noise model:

- Does full-prefix feedback raise the logical-X threshold above the existing
  one-step slice-feedback estimate?
- Can a valid constant-width Z-feedback window improve on final-only recovery?
- How does its threshold compare with full-prefix Z feedback?
- Which apparent crossings survive larger sizes, guarded finite-size fits, and
  decoder/model sensitivity checks?

The study estimates decoder-dependent thresholds. It does not reproduce a
threshold quoted by the reference paper and does not establish stability under
noisy measurement, correction, or idle operations.

## Common encoded geometry and transfer

Use `Toric3DEncodedSlice(T,Ly,Lz)` with slices `0:T`, periodic transverse
directions, and three oriented edge qubits per cell. It has

\[
n=3(T+1)L_yL_z,\qquad k=2,\qquad
d_X=(T+1)\min(L_y,L_z),\qquad d_Z=2\min(L_y,L_z).
\]

Every transfer cell uses the unchanged seven-gate schedule:

1. H on the first fresh edge;
2. two CNOTs from it to the other fresh edges; and
3. four CNOTs from old-slice parity controls to the last two fresh edges.

There are thirteen participating gate outputs per transfer cell. Independently
apply the selected X or Z fault with probability `p` after each participating
H/CNOT output. Initialization, idles, stabilizer measurements, classical
processing, and Pauli corrections are ideal.

An explicit encoded-prefix helper maps `Toric3DEncodedSlice(t,Ly,Lz)` into the
full code. Its map is an ordered local-decoder-column-to-full-frame-qubit vector:
its length equals `data_qubit_count(prefix)`, and local recovery bit `q` is
applied to full-frame qubit `map[q]`. The present time-major ordering makes this
vector contiguous, but the implementation and tests must not rely on that fact
without checking it. A feedback recovery may touch only mapped prefix qubits
and must never change a future slice.

## Sector-specific logical states

The X study prepares logical `|00>` using the existing ideal seed X-star
projection and ideal syndrome-valid Z correction from product zero. There are
no noisy seed H/CNOT locations. Failure is either residual logical-X membrane.

The Z study prepares logical `|++>` using H on every seed edge, with an
independent Z fault after every seed H, followed by the existing ideal seed
compatibility measurement and X correction. Failure is either residual
logical-Z winding.

These choices turn each channel's residual logical Pauli into actual state
failure. They also mean X and Z use different seed circuits and fault-location
sets. Protocols are directly comparable within one channel; the numerical X
and Z thresholds are not an apples-to-apples comparison of otherwise identical
circuits. Metadata and reports must state this asymmetry.

## Protocol 1: no intermediate feedback

Prepare every transfer slice without syndrome-dependent action. At the end:

- X faults: measure the full plaquette syndrome and apply validated BP+OSD X
  recovery.
- Z faults: measure the full star syndrome and apply validated unit-data-weight
  MWPM Z recovery.

This protocol is called `none`, rather than `final_only`, in the new API. The
artifact labels may say “no feedback (final decode only).”

## Protocol 2X: one-step slice-local X feedback

Retain the established paper-like X protocol so the existing 1.48% result is a
direct baseline.

Before transfer `t`, measure all six-edge Z compatibility checks on active
slice `t-1`. Unit-weight MWPM returns a complete syndrome-valid matching path.
The primary `one_step` route toggles the first edge and, when distinct, the last
edge of every path. Map the propagation image of that partial route to the
fresh targets after the noisy transfer. Old physical data are left in their
current frame and the final decoder handles their residual error.

The complete MWPM recovery must reproduce the measured syndrome even though
the deliberately partial route need not clear it. Record the residual active
syndrome. `full_path` target feed-forward remains a named sensitivity and is
not relabeled as full-prefix feedback.

This protocol is local in quantum support but uses a global classical matching
calculation on the periodic slice.

## Why Z feedback needs a wider local window

The X-star checks of a prepared prefix are ordered as:

- lower seed-white triads;
- bulk six-edge sheets `A_0,...,A_(r-1)`; and
- frontier triads `B_r`.

An edge on physical slice `r` joins `A_(r-1)` to `B_r`. An edge on slice
`r-1` joins `A_(r-2)` to `A_(r-1)`, with the lower seed-white sheet replacing
`A_-1` at the first transfer.

A literal one-physical-slice Z decoder is therefore invalid under the retained
noise model. Transfer faults occur on both new targets and old controls. Old-
control Z faults place one defect on `A_(r-2)`, outside the two check sheets
adjacent to the new slice. The restricted syndrome can be odd and cannot be
reproduced by a recovery on the new slice. Exhaustive enumeration for one
`L=3` transfer gives an odd restricted syndrome for 54 of 117 eligible output
locations. All eligible single-fault syndromes have both endpoints inside a
three-check-sheet constant-width band.

This extra thickness follows from post-gate faults on both outputs. It does not
change the transfer circuit or noise convention.

## Protocol 2Z: constant-width local Z feedback

Use a complete recovery in a moving window:

- At seed time `r=0`, measure both colors of seed X stars and unit-weight
  MWPM-decode on all seed edges.
- After transfer `r`, for `1 <= r < T`, measure the three check sheets
  `{A_(r-2), A_(r-1), B_r}`, with the seed-white sheet replacing `A_-1` at
  `r=1`. Decode using only edges on physical slices `r-1:r`.
- Validate that the local recovery exactly reproduces the measured window
  syndrome, then apply its Z correction to those two slices before transfer
  `r+1`.

For `r>=2`, the measured `A_(r-2)` check sheet itself touches qubit slices
`r-2:r-1`, so the complete measurement support reaches physical slices
`r-2:r`. Candidate recovery edges remain restricted to the newest two slices
`r-1:r`. At `r=1`, the lower-white, `A_0`, and `B_1` checks touch only slices
`0:1`.

Every local-decoder column has two detector endpoints. The invariant has an
explicit induction. The seed recovery first leaves zero seed-star syndrome.
Assuming the full `r-1` prefix has zero star syndrome, all faults introduced by
transfer `r` are supported on qubit slices `r-1:r`, and their syndrome lies in
`{A_(r-2),A_(r-1),B_r}`. A validated window recovery therefore restores zero
full-prefix star syndrome. After applying it, recompute the complete prefix
star syndrome and require zero as an audit assertion; this recomputation is not
decoder input. Tests must establish the invariant for all eligible single
faults and randomized multifault histories.

Zero syndrome does not determine the logical class. Seed or window MWPM may
introduce or leave a logical Z. That parity is propagated in the Pauli frame,
may combine modulo two with later events, and enters final logical scoring; it
is never truth-corrected.

Use the complete MWPM route. A one-step Z route is excluded from the primary
study because it could leave a defect on the oldest check sheet just as that
sheet retires. Such a route would require a separately proved causal rule.

The window contains `O(Ly*Lz)` checks and data independent of transfer depth,
so it is the Z-sector analogue of weak/local feedback. Measurement support
touches at most three physical slices, while correction support uses the newest
two.

## Protocol 3X: full-prefix X feedback

After transfer `t`, for every `1 <= t < T`:

1. construct the encoded prefix through slice `t`;
2. measure all completed prefix Z plaquettes;
3. run the existing validated BP+OSD X decoder on the full prefix;
4. XOR the returned recovery into every mapped prefix qubit; and
5. verify that the complete prefix plaquette syndrome is zero.

For each local bit `q` set in the recovery, toggle only
`frame.x[prefix_map[q]]`; leave `frame.z` and every future qubit unchanged. The
corrected frontier is then used as the old-control slice for transfer `t+1`, so
ordinary Pauli-frame propagation carries its effect through the fixed circuit.
The logical class chosen by the decoder is not corrected using hidden fault
information. A prefix logical X is propagated in the Pauli frame, may combine
modulo two with later events, and enters final logical scoring without being
truth-corrected.

The logical-`|00>` seed has no X fault locations, so no seed-time X decode is
needed. Strong-feedback checkpoints are exactly `t=1:T-1`. At `t=T`, invoke the
common full-code BP+OSD decoder exactly once; this is the final recovery, and
logical scoring uses its residual without a second decoder call.

Full-prefix checkpoint history contains the pre-correction syndrome, validated
recovery, ordered prefix map, and post-correction syndrome. The post-correction
plaquette syndrome and the active-frontier compatibility syndrome must both be
zero.

## Protocol 3Z: full-prefix Z feedback

At seed time, use the same seed-star graph and seed correction as local Z
feedback. After every transfer `t` with `1 <= t < T`:

1. measure every X star of the complete encoded prefix;
2. run unit-data-weight MWPM over all prefix edges;
3. validate the full recovery against the measured syndrome;
4. apply it to all mapped prefix qubits; and
5. require zero complete prefix star syndrome.

At `t=T`, use the common final full-code MWPM decoder. Prefix recoveries may
choose a logical Z class. Its parity is propagated in the frame, may combine
modulo two with later events, and enters final scoring without being
truth-corrected.

## Decoder settings and information boundary

Every decoder receives only measured syndromes, geometry, the declared model
parameters, and its fixed settings. It cannot inspect sampled fault bits.

X BP+OSD uses the established baseline:

- minimum-sum BP;
- 50 iterations;
- scaling factor 1;
- parallel schedule;
- one OpenMP thread; and
- `OSD_CS` order 2.

The prior is uniform declared physical `p` on every data qubit, numerically
clamped to `[1e-12,1-1e-12]` as in the existing decoder. It neglects nonuniform
gate exposure, circuit correlations, and correlations conditioned on earlier
feedback. Z MWPM uses unit data-edge weights and deterministic tie handling,
which likewise neglects circuit correlations. These limitations apply to
prefix and final decoders and must be recorded.

Any full recovery whose syndrome differs from its input stops the point with
dimensions, channel, probability, protocol, checkpoint, shot/batch, and seed
context. Partial one-step X routing records its expected residual syndrome but
still requires the underlying full MWPM recovery to be valid.

## Paired replay and simulation engine

For every shot at every `(dimensions,p,channel)` point, sample one immutable
fault record and replay that record through `none`, `slice_local`, and
`full_prefix`.

X records contain the thirteen transfer-output Bernoulli bits per cell. Z
records contain the three seed-H bits per seed cell plus the same transfer
output locations. Adaptive Z feedback requires checkpoint frames, so the
existing final-frame sparse shortcut cannot serve as the replay engine.
Implement incremental propagation, preferably vectorized over shot batches,
and use PyMatching batch decoding at each Z checkpoint.

Cache the existing X slice matcher, one Z seed-star decoder, one post-transfer
Z-window decoder per transverse geometry, and the Z growing-prefix decoders by
dimensions because they do not depend on `p`. The `r=1` boundary-labeled window
is graph-isomorphic to later windows but must be tested explicitly. Cache X
BP+OSD prefix decoders by `(dimensions,p,settings)`; the `t=T` member is also
the final X decoder. No decoder is rebuilt per shot. The supplied RNG owns all
draws; integer seeds and deterministic point substream seeds are provenance
and resume controls.

## Public API and compatibility

Add a separate diagonal feedback-comparison subsystem. Preserve every existing
public entry point and artifact format.

The new public API will provide:

- encoded-prefix construction and mappings;
- local-window Z and full-prefix X/Z decoder builders;
- validated single-protocol and three-protocol replay functions;
- a configurable scan accepting independent `(T,Ly,Lz)` dimensions,
  probabilities, channels, shot count, batch size, RNG, seed metadata, BP+OSD
  settings, X routing, checkpoint/resume paths, and progress output;
- save and plot functions; and
- a runnable command-line example.

Use protocol symbols `:none`, `:slice_local`, and `:full_prefix`. Do not reuse
the two-protocol `Toric3DMetricEstimate` fields. The new result model stores
protocol marginals and explicit joint outcomes. X routing affects only
`:slice_local`; `:none` and `:full_prefix` ignore it. `:full_path` is a
slice-local sensitivity, not a fourth protocol. Routing is nevertheless part
of saved configuration and exact resume matching.

## Statistics and artifacts

For each logical-1, logical-2, and union state-failure metric, store the eight
joint Boolean outcome patterns for the three protocols. Derive from them:

- each protocol's failures and rate;
- 95% Wilson intervals;
- every pairwise discordant count;
- paired rate differences and uncertainty; and
- exact identities needed to audit the marginals.

Save one raw CSV row per `(dimensions,p,channel,metric)` plus complete metadata,
checkpoint state, and timing. Metadata includes geometry, logical state, check
sets, correction support, protocol timing, fault-location count, RNG policy,
decoder versions/settings, limitations, invalid-recovery count, source hashes,
and Git state. Before aborting on an invalid recovery, atomically persist a
failure record with the recovery context and incremented invalid-recovery
count; successful completed artifacts must report zero.

Figures show all three protocols together for each channel, individual logical
rates, and paired differences. Export PDF, SVG, and PNG. Coincident or
unresolved curves must remain visible and labeled rather than being displaced
for presentation.

The threshold fitter treats each channel/protocol separately using binomial
likelihood. It requires at least four sizes and a bracketed finite-size
crossing. It withholds `p_c` for missing crossings, optimizer or parameter-bound
failures, nonmonotone fits, estimates outside the simulated window, deviance
goodness-of-fit p-value below `0.01`, or fewer than 80% usable bootstrap fits.
Use a quadratic scaling function as the primary model and report cubic,
minimum-size, probability-window, and leading `L^-1` correction sensitivities.
Keep conditional bootstrap uncertainty separate from model spread. Never force
a threshold for a protocol without supported crossings.

## Staged numerical study

First compile and time independent 20-shot points at approximately the known
critical regions:

- X: `L=9,p=0.015`;
- Z: `L=17,p=0.006`.

Record setup, decoder construction, propagation, local feedback, prefix
feedback, and final decode time separately. Full-prefix X requires `O(T)`
BP+OSD calls per strong-protocol shot, while adaptive Z removes the old
final-only sparse shortcut. Do not start a large grid before these timings.

The broad scout uses 200 paired triplets per point:

- X: `L=3,5,7,9` and
  `p={0,.0001,.0003,.001,.003,.005,.0075,.01,.015,.02,.03,.05,.075,.10}`;
- Z: `L=5,9,13,17` and
  `p={0,.001,.002,.003,.004,.005,.006,.0075,.01,.015,.02,.03}`.

Extend upward or downward if a protocol's transition is not bracketed. Freeze
production windows only after the scout. Use at least five sizes and seven
probabilities for a threshold fit when runtime permits, starting with 1,000
paired triplets per point for X and 2,000 for Z. Top up critical points when
Wilson precision or fit diagnostics justify it. Save batches incrementally and
exclude warmup/timing shots from scientific counts.

If the requested precision would require impractical runtime, complete the
largest defensible grid, report timing and interval widths, and withhold any
unsupported threshold rather than extrapolating.

## Verification requirements

Tests must establish:

1. encoded-prefix mappings, ranks, stabilizer commutation, logical pairing, and
   periodic transverse length two;
2. prefix plaquettes span the active X compatibility checks;
3. the Z local-window columns have exactly two endpoints and the stated
   measurement and correction support at the seed and every checkpoint;
4. every seed fault lies in the seed graph, and every transfer-`r` Z fault has
   syndrome in its associated three-sheet band; for `r<T` this is the feedback
   checkpoint, while `r=T` is handled by the common final decoder. This includes
   the 54/117 `L=3` counterexamples to a two-sheet active-slice restriction;
5. noiseless preparation and zero syndrome after every local/full checkpoint;
6. every eligible single X and Z gate-output fault on a small instance under
   all three protocols, with `T>=2` so an intermediate full-prefix X recovery
   is propagated;
7. randomized multifault histories preserve the Z sliding-window invariant;
8. recoveries reproduce exactly their measured syndromes and future qubits are
   untouched;
9. applying physical corrections agrees with applying and propagating their
   Pauli-frame representation;
10. prefix logical parities propagate into the final parity calculation rather
    than being truth-corrected;
11. immutable three-way replay, zero-noise success, seeded reproducibility,
    batching, interruption, and exact resume;
12. joint-pattern counts reproduce all protocol marginals and pairwise
    discordances;
13. threshold fitting accepts valid synthetic crossings and rejects every
    guarded failure mode; and
14. the complete existing 2D and 3D regression suite remains green.

The X tests additionally require bit-for-bit equality between the new `:none`
and `:slice_local` replay and the existing `toric3d_x_feedback_pair` for
identical records, including pre-transfer timing, seed-slice no-op behavior,
path tie handling, and both routing policies. For `T=1`, `:full_prefix` must
equal `:none` because there is no intermediate checkpoint.

## Documentation and interpretation

Add a runnable example, an experiment report, and a concise update to the 3D
section of `document/note.tex`. The report distinguishes:

- final decoding from feedback during preparation;
- one-step X feed-forward from complete local recovery;
- constant-width Z feedback from full-prefix feedback;
- logical state failure from Pauli-frame diagnostics;
- within-channel protocol comparisons from cross-channel threshold numbers;
  and
- conditional Monte Carlo intervals from decoder and scaling-model
  uncertainty.

The strong protocol is a useful upper benchmark under ideal operations, but it
repeatedly measures and decodes a growing volume. Its total check participation
scales as `O(T^2 Ly Lz)`, compared with `O(T Ly Lz)` for constant-width local
feedback. That resource distinction is part of the conclusion rather than a
hidden implementation detail.
