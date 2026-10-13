# Diagonal H+5 vertex-projector feedback experiment

**Date:** 2026-10-13  
**Status:** approved four-arm design

**Purpose:** determine whether ideal slice feedback gives a nonzero logical
threshold for the direct cubic-vertex projector under the existing
post-gate-output noise model.

## Scientific question

Replace the derived diagonal plaquette-transfer gadget in the primary new
experiment by the literal cubic-vertex projector. Compare four protocols
under paired circuit faults:

1. final recovery only (`none`);
2. decoder-guided one-step defect transport (`one_step`);
3. syndrome-complete two-slice recovery (`window_complete`); and
4. recovery over the entire completed prefix (`full_prefix`).

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
threshold scan. The X-sector `one_step` arm implements its qualitative rule:
measure paired point defects in one slice, use MWPM to pair them, and move
each endpoint one graph step. The Z-sector `one_step` arm is a causal
open-frontier extension required by this H+5 orientation and noise model.
This specification defines all remaining choices and labels the resulting
experiment as paper-motivated rather than a reproduction of a supplied
circuit.

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
a physical-edge-ID-ordered breadth-first spanning tree on the full star
vertex set, excluding edges in the union of the two declared logical-X
supports, rooted at `A_(T-1,0,0)`. The cut graph stays connected: direction-3
edges connect slices and nonwrapping transverse edges connect the rectangular
transverse domain. For each upper star `U`, let `P_U` be its unique
tree path to that root. For an observed upper-outcome vector `a`, define

\[
C(a)=\bigoplus_{U:a_U=1}P_U.
\]

Every path has zero overlap with both declared logical-X supports by
construction, and therefore so does every correction. A tree alone would not
guarantee this for open paths.
It flips every indicated upper outcome. If `a` has odd parity, it transfers
the unmatched endpoint to the declared bulk root, where the final decoder can
observe it. For a clean even-parity outcome it leaves every completed bulk
star unchanged.

For each shot, sample a clean upper outcome `b` uniformly from the even-parity
subspace (sample the first `L_uL_v-1` bits uniformly and set the last to their
parity), store it with the immutable replay record, and reuse it for all four
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
the same immutable record through all four protocols. Decoder inputs contain
only physically available measured outcomes or syndromes, geometry, declared
`p`, and decoder settings. Feedback recoveries are applied to the Pauli frame
before later gates. A physical-Pauli reference on the small lattice must agree
with frame propagation and with the explicit clean/noisy boundary branches.

## Feedback protocols

All intermediate feedback occurs after layers `1:T-1`. Every protocol uses
the same final full-code recovery after ideal upper-boundary completion.

### No feedback

Propagate the complete fault record through all H+5 layers, apply the common
ideal upper-boundary instrument, measure the full syndrome, and decode once.
Here `none` means no *intermediate bulk* feedback; it is not a fully unitary
or unassisted boundary preparation.

### One-step X transport

Write `Q_t` for the `L_uL_v` six-edge Z compatibility checks contained in
`E_t`. After layer `t`, these are noiseless `+1` stabilizers of the prepared
prefix. Every edge in `E_t` touches exactly two `Q_t` checks, including at
transverse period two, so their point-defect syndrome has even parity on a
periodic triangular matching graph.

Use unit-weight MWPM to pair the measured defects. Lift every pair to a
deterministic shortest path of physical `E_t` edge identities, with
lexicographic tie breaking. The XOR of all complete paths is the decoder's
`full_recovery`; require its `Q_t` syndrome to reproduce the measurement.
The physically applied `route` is deliberately only one step:

- for a path of length one, apply its single edge and annihilate the pair;
- for a longer path, apply only its first and last edge, moving both defects
  one graph step toward one another.

The partial route is not required to reproduce the original syndrome. Record
and verify its predicted residual syndrome instead. It implements the paper's
within-slice MWPM and endpoint-motion rule. It does not reproduce the paper's
iterative transport architecture: under the fresh-control-to-old-target H+5
orientation, residual `Q_t` defects do not advance to `Q_{t+1}` and receive
only this one contraction before final decoding. Likewise, the rule does not
see X faults deposited on old targets in `E_{t-1}`, but those qubits have
finished all H+5 gates and such X faults cannot proliferate; they remain local
noise for final recovery. Do not reuse the H+6 target-feedforward map `F(C)`:
the H+5 CNOT orientation and correction action are different.

### Complete-window X recovery

Write `P_j` for all `3L_uL_v` physical Z plaquettes between edge slices
`E_{j-1}` and `E_j`. After layer `t`, measure

- `P_1` at `t=1`; and
- `P_{t-1} union P_t` at `t>=2`.

Restrict recovery columns to `D_t=E_{t-1} union E_t` and decode the resulting
binary incidence matrix with BP+OSD. Require the recovery to reproduce the
complete window syndrome and, after application, require every completed
plaquette in `P_1:P_t` to have zero syndrome.

This is the smallest syndrome-complete X window for the H+5 layer. It is a
stronger operation than one-step transport: it also handles old-target faults
and clears every completed physical plaquette instead of intentionally
leaving transported point defects.

### Full-prefix X recovery

After layer `t`, measure every physical Z plaquette supported on edge slices
`E_0:E_t`. Use BP+OSD to return an X recovery restricted to those prepared
edges. Apply it and require the complete measured prefix plaquette syndrome to
vanish. The decoder must not access future edges.

### One-step Z transport

Write `W` for the lower white boundary triads and `A_j` for bulk-star sheet
`j`, supported on `E_j union E_{j+1}`. At checkpoint `t`, define retiring and
frontier detector sheets

\[
(R_t,F_t)=
\begin{cases}
(W,A_0),&t=1,\\
(A_{t-2},A_{t-1}),&t\geq2.
\end{cases}
\]

The physical strip graph is

\[
R_t\xleftrightarrow{E_{t-1}}F_t
\xleftrightarrow{E_t}B_t,
\]

where `B_t=A_t` is the unprepared future sheet and is not measured. Thus an
`E_t` column has one measured endpoint, and the measured syndrome may be odd.
Treat `W` as a detector sheet, not as a matching boundary. Compute
detector-to-detector distances and paths in the periodic `R_t`--`F_t` graph
using physical `E_{t-1}` edges. Use a parity-constrained temporal-boundary
matching policy: allow no boundary path for even syndrome and add exactly one
dummy future-boundary mate for odd syndrome. Its distance from a detector is
the shortest path to an `F_t` vertex followed by one physical `E_t` edge. This
is MWPM on a metric closure with an optional single dummy, not ordinary
open-boundary MWPM, which can terminate any number of paths. Lift the matches
to deterministic physical-edge paths and require the XOR of their complete
paths to reproduce the measured `R_t union F_t` syndrome.

The applied route is a causal endpoint trim:

- apply a path's endpoint edge only when that endpoint lies in `R_t`;
- apply no endpoint edge for an endpoint in `F_t` or at `B_t`; and
- a length-one `R_t`--`F_t` path contributes its shared edge once.

The route is therefore supported only on `E_{t-1}`, whose H+5 gates are
finished. It must make the residual syndrome zero on `R_t`; residual defects
may remain on `F_t`. At the next intermediate checkpoint that frontier becomes
the new retiring sheet; after checkpoint `T-1`, the common final decoder
handles it instead. In this Z arm, `one_step` therefore means causal
retiring-sheet transport: each `R_t` endpoint advances one edge, while an
`F_t` endpoint waits until it retires. Applying both the first and last edge
of an arbitrary strip path is forbidden because it can move a frontier
endpoint backward into `R_t`, after which the next window cannot reach it.
This causal rule is the Z-sector H+5 extension of the paper's one-step idea,
rather than the paper's two-endpoint periodic loop-intersection rule.

### Complete-window Z recovery

On `D_t=E_{t-1} union E_t`, use the following syndrome-complete open matching
graph after each layer:

- at `t=1`, measure `W union A_0`; and
- at `t>=2`, measure `A_{t-2} union A_{t-1}`.

Columns ending at the unprepared upper frontier have one detector endpoint,
so odd syndrome parity is valid. Here and in full-prefix Z recovery, use
ordinary unit-weight open-boundary MWPM: any number of detector paths may
terminate on the future boundary at their nonzero physical path cost.
Apply the recovery only inside `D_t`. Require it to reproduce the complete
window syndrome and, after application, require `W` and every completed star
sheet `A_0:A_{t-1}` to have zero syndrome.

For both complete-window sectors, the completed-prefix requirement follows
inductively. Rows older than `P_{t-1}` in the X sector or `A_{t-2}` in the Z
sector are disjoint from `D_t` and cannot change during layer `t` or its window
recovery. The actual net change of the propagated Pauli frame during layer
`t`, including injected faults, is supported on `D_t` and witnesses that the
measured syndrome lies in the restricted matrix image. Thus any
syndrome-valid recovery on `D_t` clears the only completed rows that could
have changed, for arbitrary fault patterns and not merely for single faults.

For `M=L_uL_v`, the expected binary-matrix ranks are `3M-1` for the `t=1`
window-X matrix of shape `3M by 6M`, `4M` for a steady window-X matrix of
shape `6M by 6M`, and full row rank `2M` for both window-Z matrices of shape
`2M by 6M`. These expectations also apply when either transverse period is
two; duplicated coordinates must still denote distinct physical edges.

### Full-prefix Z recovery

After layer `t`, measure the lower white boundary stars and every completed
bulk-star sheet `0:t-1`. Decode on all prepared edge slices `E_0:E_t` using an
open-boundary MWPM graph. Apply the recovery anywhere in that prefix and
require every completed prefix-star syndrome bit to vanish. No upper-frontier
triad is treated as a stabilizer before final boundary completion.

Every MWPM instance uses a deterministic two-level objective: first minimize
total unit physical-edge length, then choose the lexicographically smallest
ordered matching-and-path signature among primary optima. Reconstruct a
physical recovery by XORing the selected path edge IDs. Preserve parallel
periodic edges as distinct IDs when `L_u=2` or `L_v=2`; choose the lowest edge
ID only as the secondary tie break, never by collapsing the columns. Pin and
record the matching backend and version as well as the canonical signature.

At a common checkpoint, the selected X- and Z-type measurements commute. An
X recovery does not change an X-star syndrome, and a Z recovery does not
change a Z-check syndrome. The Pauli recovery order therefore changes only a
discarded global phase. The public replay `channel` values are exactly
`:x_only`, `:z_only`, and `:both`. A pure-channel call rejects a record with
any fault in the opposite component before propagation. For `channel=:both`,
capture both syndrome sets from the same pre-correction frame, compute both
decisions, then XOR both routes into the frame. Although the first study scans
X-only and Z-only noise separately, verify `:both` compatibility for both
one-step and complete-window policies so that the public API has an
unambiguous combined-channel meaning.

## Final decoding and observables

After the final layer, apply the explicit relative-outcome upper-triad channel
defined above and use the existing full-code decoders:

- unit-weight MWPM for Z recovery from the full star syndrome; and
- BP+OSD for X recovery from the full plaquette syndrome.

For every intermediate and final BP+OSD instance, use minimum-sum BP, 50
iterations, scaling factor one, parallel schedule, one thread, and `OSD_CS`
order two. Use the declared physical rate `p` as the uniform prior except at
the numerical endpoints, where

\[
p_{\rm eff}=\operatorname{clamp}(p,10^{-12},1-10^{-12}).
\]

Record both `p` and `p_eff` for every BP+OSD instance. These baseline
decoders neglect circuit-induced correlations. Every decoder's complete
intermediate path or recovery and every final recovery must reproduce its
input syndrome: independently recompute `H*c` from the returned physical
recovery bits and require `H*c == s`. The physically applied one-step routes
are intentionally partial. Apply each route to a copy of the complete frame,
remeasure independently, and require

\[
s_{\rm after}=s_{\rm before}\oplus H r.
\]

For X, save the resulting `Q_t` residual. For Z, require zero syndrome on
every retired completed sheet and allow residual measured syndrome only on
`F_t`; make no syndrome claim about unprepared `B_t`. On failure, atomically
save before/after syndromes, complete paths, applied routes, physical supports,
and replay metadata, then stop the affected point.

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

For each channel and protocol, predeclare failure of either transverse logical
pair as the primary threshold observable. Individual logical components and
the full Pauli class are secondary diagnostics. Report all six paired
protocol contrasts at each common `(T,Lu,Lv,p)` point, with the three
feedback-versus-`none` contrasts highlighted as primary feedback comparisons.

One-step, complete-window, and full-prefix protocols differ in measured check
sets, correction support, and whether residual syndrome is intentionally
retained. The X one-step arm uses MWPM while its complete-window and
full-prefix arms use BP+OSD; all Z arms use channel-appropriate MWPM graphs.
Their threshold differences are therefore end-to-end protocol comparisons,
not isolated estimates of correction range alone.

Likewise, do not compare H+5 numbers directly with historical H+6 thresholds
as a controlled circuit-only change: the seed, boundary instrument, and old
channel-specific seed noise differ.

## Verification

Development follows test-first red/green cycles. Required checks are:

1. exact six-operation sequence and support for every cell, `6TL_uL_v`
   operations, and 11 participating outputs per cell;
2. every pivot and other outgoing edge is fresh when first used, cells within
   a layer are disjoint, and no correction touches unprepared edge slices
   `E_{t+1}:E_T`; the prepared frontier `E_t` is allowed where declared;
3. boundary-check commutation, stabilizer ranks, `k=2`, distances, logical
   commutation, and length-two periodic edges;
4. noiseless preparation of all final checks and logical Z, plus the expected
   completed-check set at every intermediate checkpoint;
5. an explicit regression demonstrating that the old honeycomb seed is not a
   valid drop-in boundary for the H+5 circuit;
6. exhaustive X and Z faults at all 11 output locations on a small lattice in
   every applicable intermediate-feedback arm;
7. `Q_t` stabilizer membership, even syndrome parity, physical-edge path
   identities, complete-path validity, one-step endpoint motion, predicted
   residual syndrome, and terminal behavior of missed old-target X faults;
8. even and odd Z-strip matching, exactly one dummy future-boundary mate for
   odd parity, complete-path validity, zero retiring-sheet syndrome after the
   causal trim, residual confinement to the frontier among measured completed
   sheets, and an explicit counterexample showing why symmetric first/last
   trimming is invalid;
9. exhaustive reachable one-step syndromes for `L_u=L_v=2`, including every
   even `Q_t` syndrome and every Z-strip syndrome, plus randomized multifault
   multi-checkpoint histories at transverse periods two and three;
10. syndrome-image and restricted-matrix-rank checks for both complete windows
   at periods two and three, plus X/Z completed-prefix clearing;
11. full-prefix X/Z completed-prefix clearing;
12. commuting intermediate X/Z measurements and equality of both Pauli
    correction orders in the binary frame;
13. equality of physical correction and propagated-frame correction;
14. stabilizer-rank verification that clean upper outcomes are uniform on the
    even-parity subspace, all allowed clean upper branches, physical use of
    `b xor u`, equality of explicit two-branch and optimized relative
    completion, zero homology of the tree correction, and preservation of all
    four logical generators;
15. final BP+OSD/MWPM recovery validity, zero-noise success, endpoint prior
    clamping, four-way paired replay, seeded reproducibility, incremental
    resume, and invalid-recovery aborts;
16. existing H+6 and Cartesian 3D regression tests.

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
- combined final-only/one-step/window-complete/full-prefix figures in PDF,
  SVG, and PNG for X and Z;
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
