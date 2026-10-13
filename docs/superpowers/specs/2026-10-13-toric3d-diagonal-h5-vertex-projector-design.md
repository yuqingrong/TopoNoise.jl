# Diagonal H+5 vertex-projector feedback experiment

**Date:** 2026-10-13  
**Status:** revised after protocol-consistency review; implementation pending

**Purpose:** determine whether ideal slice feedback gives a nonzero logical
threshold for the direct cubic-vertex projector under the existing
post-gate-output noise model.

## Scientific question

Replace the derived diagonal plaquette-transfer gadget in the primary new
experiment by the literal cubic-vertex projector. Compare three protocols
under paired circuit faults:

1. final recovery only (`none`);
2. constant-width active-frontier recovery (`slice_local`); and
3. recovery over the entire completed prefix (`full_prefix`).

Run X-only and Z-only faults separately and report residual logical Pauli
classes, rather than inferring them from the fidelity of one selected logical
state. The existing H+6 implementation and its saved results remain available
under their current names for reproducibility, but all new artifacts must call
it the *plaquette-transfer* circuit. The H+5 implementation is the primary
direct-projector experiment.

This is not a numerical reproduction of a circuit supplied by
*Noise-resilient sequential circuits for generating quantum order*. That
paper motivates diagonal slicing and within-slice feedback but does not give
an explicit 3D toric-code gate decomposition, boundary preparation, or
threshold scan. This specification defines those missing choices.

## Geometry and target code

Use diagonal coordinates

\[
t=x+y+z,\qquad u=x,\qquad v=y,
\]

with open transfer coordinate `t=0:T` and periodic transverse coordinates
`u=0:L_u-1`, `v=0:L_v-1`. Define

\[
s_1=(1,0),\qquad s_2=(0,1),\qquad s_3=(0,0).
\]

The edge `E[d,t,r]`, with `r=(u,v)`, ends at vertex `V(t,r)` and starts at
`V(t-1,r-s_d)`. Each edge slice contains `3L_uL_v` qubits. The final slab is
the existing `Toric3DEncodedSlice` cell complex:

\[
n=3(T+1)L_uL_v,\qquad k=2,
\]
\[
d_Z=2\min(L_u,L_v),\qquad
d_X=(T+1)\min(L_u,L_v).
\]

Its X checks are:

- lower-boundary white triads on edge slice zero;
- six-edge cubic stars at vertex slices `0:T-1`; and
- upper-boundary triads on edge slice `T`.

Its Z checks are the four-edge cubic plaquettes between adjacent edge slices.
Keep the existing two transverse logical pairs and preserve distinct physical
edges when either transverse period equals two.

The legacy fields `Ly,Lz` currently store `L_u,L_v`. New public metadata and
documentation must use `Lu,Lv`; compatibility wrappers may retain the old
field names internally until a separate API migration.

## Ideal boundary preparation

A literal H+5 projector is incompatible with the existing encoded honeycomb
seed: it preserves both seed-star colors and does not prepare the required
upper frontier or physical plaquettes. Therefore the H+5 experiment uses a
different, symmetric boundary completion:

1. initialize every data edge in product zero;
2. ideally measure all `L_uL_v` independent lower white X triads and apply an
   ideal Z correction that fixes their eigenvalues;
3. apply the noisy bulk vertex projectors described below; and
4. ideally measure all upper frontier X triads and pass all outcomes to the
   completion instrument below; the clean outcomes contain one dependency.

The product-zero input fixes all Z plaquettes and both logical-Z cycles. The
boundary projections, measurements, classical processing, and Pauli
corrections are ideal under the agreed noise convention. They introduce no
H/CNOT fault locations. The final frontier projection is necessary because a
bulk projector always completes the vertex behind the frontier; the upper
triads have no fresh outgoing data edge that could serve as a pivot.

The `L_uL_v` lower white triads are independent at initialization. After the
lower and bulk projectors, the `L_uL_v` upper outcomes have one parity
constraint, so only `L_uL_v-1` are independent. They are otherwise
intrinsically random and are not an error syndrome.

Model upper completion as an explicit physical measurement instrument. Build
a deterministic lexicographic spanning tree of the full star graph, rooted at
one declared final bulk star. For each upper star `U`, let `P_U` be its unique
tree path to that root. For an observed upper-outcome vector `a`, define

\[
C(a)=\bigoplus_{U:a_U=1}P_U.
\]

This map has fixed zero logical-Z homology because its support lies in a tree.
It flips every indicated upper outcome. If `a` has odd parity, it transfers
the unmatched endpoint to the declared bulk root, where the final decoder can
observe it. For a clean even-parity outcome it leaves every completed bulk
star unchanged.

For each shot, sample a clean upper outcome `b` uniformly from the even-parity
subspace (sample the first `L_uL_v-1` bits uniformly and set the last to their
parity), store it with the immutable replay record, and reuse it for all three
protocols. If the current Z frame flips the upper outcomes by `u`, the physical noisy outcome is
`a=b xor u`; the boundary controller receives only `a` and applies `C(a)`.
The simulator compares this branch with the matched clean branch, which
applies `C(b)`. The relative post-completion Z frame is therefore

\[
z'=z\oplus C(b\oplus u)\oplus C(b).
\]

No decoder or boundary controller receives `u`, the clean counterfactual, or
the sampled gate faults. Since the declared tree map is linear, an optimized
implementation may reduce the displayed expression to `z xor C(u)`, but it
must retain `b` in provenance and test exact equality with the explicit
two-branch construction. Record both physical boundary corrections in replay
history, then perform final recovery from the relative full-code syndrome.
Residual logical parities include every intermediate correction, the matched
boundary-correction difference, and final recovery.

Small-system stabilizer algebra must verify both sides of this construction:
before upper projection, the lower stars, all bulk stars, every plaquette, and
both logical-Z cycles are stabilized while the upper triads are not; after a
clean measurement branch and `C(b)`, every target-code generator and logical
Z is stabilized. For noisy branches, all target checks are required to be
restored only after final recovery. Injecting any of the four logical Pauli
generators immediately before completion must preserve its relative logical
class after the paired boundary map.

## Direct vertex-projector circuit

For old vertex `V(t-1,r)`, define its incoming old edges and outgoing fresh
edges by

\[
e_d=E[d,t-1,r],\qquad n_d=E[d,t,r+s_d],\qquad d=1,2,3.
\]

Choose `n_1` as the fresh data pivot and apply the fixed gate order

\[
\begin{aligned}
&H(n_1),\\
&\operatorname{CX}(n_1\!\to n_2),\quad
 \operatorname{CX}(n_1\!\to n_3),\\
&\operatorname{CX}(n_1\!\to e_1),\quad
 \operatorname{CX}(n_1\!\to e_2),\quad
 \operatorname{CX}(n_1\!\to e_3).
\end{aligned}
\]

On fresh zero targets this implements

\[
\frac{I+A_{t-1,r}}{\sqrt 2},\qquad
A_{t-1,r}=\prod_{d=1}^3X_{e_d}X_{n_d},
\]

on the ideal input subspace. There is one H and five CNOTs per vertex, hence
11 participating gate outputs. Cells within one diagonal layer are edge
disjoint and may be regarded as parallel; stored operations use deterministic
increasing `t,u,v` order. Apply feedback after a full layer and before any edge
in its frontier becomes an old target in the next layer.

Implement this as a separate `Toric3DVertexProjectorEncoder`. Do not mutate
the operation list or metadata of `Toric3DSliceEncoder`, because the existing
H+6 artifacts must remain exactly reproducible.

## Noise and paired replay

For each selected channel, independently place a rate-`p` Pauli after every
participating gate output:

- one possible fault after the H output; and
- one possible fault on each of the two CNOT outputs.

Thus every transfer cell has 11 Bernoulli locations. Run X-only and Z-only
channels separately. Initialization, boundary projections, idle qubits,
syndrome measurements, decoding computation, and applied Pauli corrections
are ideal.

Sample faults and the clean upper-boundary branch before decoding and replay
the same immutable record through all three protocols. Decoder inputs contain
only physically available measured outcomes or syndromes, geometry, declared
`p`, and decoder settings. Feedback recoveries are applied to the Pauli frame
before later gates. A physical-Pauli reference on the small lattice must agree
with frame propagation and with the explicit clean/noisy boundary branches.

## Feedback protocols

All intermediate feedback occurs after layers `1:T-1`. Every protocol uses
the same final full-code recovery after ideal upper-boundary completion.

The weak protocol uses one common two-edge-slice data window in both Pauli
sectors,

\[
D_t=E_{t-1}\cup E_t.
\]

This is the smallest syndrome-complete window for the H+5 circuit: a layer-`t`
gate acts only on old targets in `E_{t-1}` and fresh edges in `E_t`, and Pauli
propagation through the rest of that layer never leaves those two slices. The
X- and Z-sector measurements below are commuting CSS checks. Their recoveries
may overlap, but their order changes only an ignored Pauli phase, so the same
window also defines a compatible combined-channel protocol.

### No feedback

Propagate the complete fault record through all H+5 layers, apply the common
ideal upper-boundary instrument, measure the full syndrome, and decode once.
Here `none` means no *intermediate bulk* feedback; it is not a fully unitary
or unassisted boundary preparation.

### Weak X feedback

Write `P_j` for all `3L_uL_v` physical Z plaquettes between edge slices
`E_{j-1}` and `E_j`. After layer `t`, measure

- `P_1` at `t=1`; and
- `P_{t-1} union P_t` at `t>=2`.

Restrict recovery columns to `D_t` and decode the resulting binary incidence
matrix with BP+OSD. Require the recovery to reproduce the complete window
syndrome and, after application, require every completed plaquette in
`P_1:P_t` to have zero syndrome.

The previous active-slice checks `Q_t` may be retained as a diagnostic, but
they do not define weak H+5 feedback: they cannot see all physical-plaquette
defects made by faults on old CNOT targets. Do not reuse either their MWPM
decoder or the H+6 target-feedforward map `F(C)`.

### Strong X feedback

After layer `t`, measure every physical Z plaquette supported on edge slices
`E_0:E_t`. Use BP+OSD to return an X recovery restricted to those prepared
edges. Apply it and require the complete measured prefix plaquette syndrome to
vanish. The decoder must not access future edges.

### Weak Z feedback

Write `W` for the lower white boundary triads and `A_j` for bulk-star sheet
`j`, supported on `E_j union E_{j+1}`. On the same recovery window `D_t`, use
the following constant-width open matching graph after each layer:

- at `t=1`, measure `W union A_0`; and
- at `t>=2`, measure `A_{t-2} union A_{t-1}`.

Columns ending at the unprepared upper frontier have one detector endpoint,
so odd syndrome parity is valid and must be matched to the open boundary.
Apply the recovery only inside `D_t`. Require it to reproduce the complete
window syndrome and, after application, require `W` and every completed star
sheet `A_0:A_{t-1}` to have zero syndrome.

For both sectors, the completed-prefix requirement follows inductively. Rows
older than `P_{t-1}` in the X sector or `A_{t-2}` in the Z sector are disjoint
from `D_t` and cannot change during layer `t` or its weak recovery. The actual
net change of the propagated Pauli frame during layer `t`, including injected
faults, is supported on `D_t` and witnesses that the measured syndrome lies in
the restricted matrix image. Thus any syndrome-valid recovery on `D_t` clears
the only completed rows that could have changed, for arbitrary fault patterns
and not merely for single faults.

For `M=L_uL_v`, the expected binary-matrix ranks are `3M-1` for the `t=1`
weak-X matrix of shape `3M by 6M`, `4M` for a steady weak-X matrix of shape
`6M by 6M`, and full row rank `2M` for both weak-Z matrices of shape
`2M by 6M`. These expectations also apply when either transverse period is
two; duplicated coordinates must still denote distinct physical edges.

### Strong Z feedback

After layer `t`, measure the lower white boundary stars and every completed
bulk-star sheet `0:t-1`. Decode on all prepared edge slices `E_0:E_t` using an
open-boundary MWPM graph. Apply the recovery anywhere in that prefix and
require every completed prefix-star syndrome bit to vanish. No upper-frontier
triad is treated as a stabilizer before final boundary completion.

## Final decoding and observables

After the final layer, apply the explicit relative-outcome upper-triad channel
defined above and use the existing full-code decoders:

- unit-weight MWPM for Z recovery from the full star syndrome; and
- BP+OSD for X recovery from the full plaquette syndrome.

For every intermediate and final BP+OSD instance, use minimum-sum BP, 50
iterations, scaling factor one, parallel schedule, one thread, and `OSD_CS`
order two, with a uniform physical-`p` prior. Record that these baseline
decoders neglect circuit-induced correlations. Every complete intermediate
and final recovery must reproduce its input syndrome. Save an atomic failure
record and stop the affected point on any mismatch.

Let `F_X,F_Z` be the final Pauli frame after all gate faults, intermediate
corrections, paired boundary completion, and final recovery. For logical pair
`i in {u,v}`, score

\[
\ell_{X,i}=\langle F_X,\overline Z_i\rangle,\qquad
\ell_{Z,i}=\langle F_Z,\overline X_i\rangle.
\]

Report both components, either-pair union rates, and the full logical Pauli
class. These are state-independent logical diagnostics. Also label logical-X
as `|00>` state-preparation failure for the selected ideal reference state;
do not label logical Z as `|00>` infidelity.

Weak and strong protocols differ in measured check sets and correction
support. Both strengths use BP+OSD in the X sector and MWPM in the Z sector,
but build different decoder instances for their respective matrices. Their
threshold difference is therefore an end-to-end protocol comparison, not an
isolated estimate of correction range alone.

Likewise, do not compare H+5 numbers directly with historical H+6 thresholds
as a controlled circuit-only change: the seed, boundary instrument, and old
channel-specific seed noise differ.

## Verification

Development follows test-first red/green cycles. Required checks are:

1. exact six-operation sequence and support for every cell, `6TL_uL_v`
   operations, and 11 participating outputs per cell;
2. every pivot and other outgoing edge is fresh when first used, cells within
   a layer are disjoint, and no correction touches future data;
3. boundary-check commutation, stabilizer ranks, `k=2`, distances, logical
   commutation, and length-two periodic edges;
4. noiseless preparation of all final checks and logical Z, plus the expected
   completed-check set at every intermediate checkpoint;
5. an explicit regression demonstrating that the old honeycomb seed is not a
   valid drop-in boundary for the H+5 circuit;
6. exhaustive X and Z faults at all 11 output locations on a small lattice;
7. open-MWPM recovery for both even and odd frontier syndromes;
8. syndrome-image and restricted-matrix-rank checks for both weak windows at
   periods two and three, weak X/Z complete-prefix clearing, commuting weak
   check measurements, and equality of the two X/Z correction orders in the
   binary frame;
9. strong X/Z complete-prefix clearing;
10. equality of physical correction and propagated-frame correction;
11. stabilizer-rank verification that clean upper outcomes are uniform on the
    even-parity subspace, all allowed clean upper branches, physical use of
    `b xor u`, equality of explicit two-branch and optimized relative
    completion, zero homology of the tree correction, and preservation of all
    four logical generators;
12. final BP+OSD/MWPM recovery validity, zero-noise success, paired replay,
    seeded reproducibility, incremental resume, and invalid-recovery aborts;
13. existing H+6 and Cartesian 3D regression tests.

The verification must separately distinguish “prepared stabilizer state” from
“decoder returned a syndrome-valid correction”; neither substitutes for the
other.

## Scan and threshold analysis

Add a runnable CLI and Julia API accepting independent `T,Lu,Lv`, probability
grid, shot count, RNG/seed, protocols, channels, BP+OSD settings, checkpoint
directory, and resume policy. Save each completed point atomically.

Run a timed scout before selecting production windows. The initial scout uses
cubic parameters `T=Lu=Lv=L`, `L=3,5,7`, 500 paired shots, and

\[
p\in\{0,10^{-4},3\!\times\!10^{-4},10^{-3},3\!\times\!10^{-3},
0.01,0.02,0.03,0.05,0.075,0.10\}.
\]

Use its crossings, saturation, and runtime to freeze separate X and Z
production grids with at least four sizes and five probabilities around every
candidate crossing. Do not select a fit window using the production outcomes.
Record the scout-to-production decision in a machine-readable selection file.

Use Wilson intervals for individual rates and paired discordance statistics
for protocol differences. A threshold fit must pass the existing guards for
adjacent-size crossings, interior monotone fit, likelihood goodness of fit,
and bootstrap acceptance. Report a diagnostic candidate rather than a
threshold whenever a guard fails. Do not infer a positive threshold merely
from feedback improvement at fixed size.

## Deliverables

- isolated H+5 geometry/encoder, feedback decoders, paired replay, and scan
  API;
- runnable scout/production CLI with incremental checkpoints;
- raw and aggregate CSV, complete metadata, selection record, and failure
  records;
- combined final-only/weak/strong figures in PDF, SVG, and PNG for X and Z;
- guarded threshold and sensitivity reports;
- an H+5 experiment report and corrected `note.tex` discussion;
- retained, explicitly labeled H+6 implementation and historical artifacts.

## Non-goals

- noisy syndrome extraction or noisy classical feedback;
- depolarizing or correlated hardware noise in the first study;
- claiming the paper specifies this exact 3D circuit;
- treating unit-weight MWPM or uniform-prior BP+OSD as circuit-level maximum
  likelihood;
- deleting or rewriting historical H+6 data.
