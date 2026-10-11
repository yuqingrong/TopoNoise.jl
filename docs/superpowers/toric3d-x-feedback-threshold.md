# Logical-X threshold of diagonal-slice preparation with feedback

This study tests whether ideal stabilizer measurement and classical feedback
inside each diagonal slice can stabilize logical `|00>` preparation against
independent X faults on gate outputs.  It compares final-only decoding with
one-step slice feedback using the same fault record for both protocols.

The completed finite-size study supports a nonzero X-sector threshold near

\[
p_{c,X}=0.01482\simeq 1.48\%
\]

per participating H/CNOT output.  This value comes from the least aggressive
accepted scaling model: a cubic logistic scaling function with sizes
`L=7,9,11,13`.  Its conditional 95% parametric-bootstrap interval is
`[0.01428,0.01543]`; 999 of 1,000 bootstrap refits passed the numerical
checks.  Accepted size and probability-window variants give
`p_c=0.01476`-`0.01511`.

The qualification on the model matters.  The planned all-size quadratic fit
is rejected by its goodness-of-fit test.  An all-size cubic alternative gives
a candidate `p_c=0.01470`, but also misses the fixed guard narrowly.  These
rejections remain in the saved analysis.  The reported 1.48% value is the
accepted `L>=7` sensitivity estimate, not a relabeled all-size fit.

Final-only preparation has no adjacent-size crossing in the production
window `0.01<=p<=0.025`; its failure probability increases rapidly with size.
No positive final-only threshold is inferred from these data.

## Encoded geometry and state

For transfer depth `T` and periodic transverse dimensions `Ly,Lz`, the
diagonal encoded slab has

\[
n=3(T+1)L_yL_z,\qquad k=2,\qquad
d_X=(T+1)\min(L_y,L_z),\qquad d_Z=2\min(L_y,L_z).
\]

The seed slice is a periodic honeycomb toric code.  Product zero already fixes
its Z compatibility checks and both logical-Z cycles.  Ideal measurement of
an independent set of seed X stars followed by an ideal syndrome-valid Z
correction prepares logical `|00>`.  The zero Pauli frame represents this
ideal projection and correction, so the seed introduces no noisy H/CNOT
locations under the agreed noise model.

Each cell then uses the fixed seven-gate transfer:

1. H on the first fresh target;
2. two CNOTs from that target to the other two fresh targets; and
3. four CNOTs from old-slice parity controls into the last two targets.

There are thirteen participating gate outputs per cell.  Each output receives
an independent Bernoulli-`p` X fault.  Initialization, idle qubits, stabilizer
measurements, classical processing, and Pauli corrections are ideal.

## What the feedback does

Before transfer `t`, the feedback protocol measures the six-edge Z
compatibility checks on active slice `t-1`.  An X error anticommutes with these
checks and creates point defects on their periodic triangular matching graph.
Unit-weight MWPM pairs the defects.  For the primary `one_step` policy, each
endpoint moves one graph edge toward its partner.  The simulator propagates
that virtual old-slice repair through the fixed transfer and applies its image
only to the fresh targets.  Every active recovery is checked against the
measured syndrome before the next slice is prepared.

This operation prevents many old-slice X errors from continuing through later
CNOT controls.  It does not remove the last slice's new errors, so both
protocols still use the same final decoder.  The feedback route is local on
the quantum data, while the classical MWPM computation is global and is not a
constant-depth spatially local channel.

The `full_path` sensitivity applies the complete MWPM route within a slice.
It is a stronger, spatially nonlocal correction policy and is kept separate
from the paper-like one-step result.

## Final recovery and observable

After the transfer, measure the complete physical Z-plaquette syndrome and
decode a full data-qubit X recovery with BP+OSD.  The fixed baseline settings
are:

- minimum-sum BP with 50 iterations;
- scaling factor 1 and parallel updates;
- one OpenMP thread; and
- `OSD_CS` order 2.

The prior is uniform physical `p` on every final data qubit.  The decoder sees
only the syndrome, geometry, declared `p`, and these settings.  Every recovery
must reproduce its input syndrome or the scan stops.

After recovery, the residual frame is paired with both dual logical-Z cycles.
Failure means that either logical-X membrane has odd parity.  Since the seed
is logical `|00>`, this union is an actual state-preparation failure.  It is
different from the earlier rough-slab logical-Z diagnostic, for which a
logical Z acts trivially on `|0_L>`.

## Staged scan

The broad pilot used `L=3,5,7,9`, eleven probabilities from zero through
0.02, and 1,000 paired shots at every point.  Its feedback crossings moved
upward with size, while the final-only crossings moved toward smaller `p`.
An extended 200-shot scout at `L=9,11,13` selected the production window.

Production used cubic `T=Ly=Lz=L` with

\[
L\in\{5,7,9,11,13\},\qquad
p\in\{0.0100,0.0125,0.0150,0.0175,0.0200,0.0225,0.0250\}.
\]

Each of the 35 points contains 1,000 independent paired shots, for 35,000
unique fault records.  The three largest sizes pool an independent 200-shot
scout with 800 new shots; smaller-size components and every seed are recorded
in the selection metadata.  Every recovery validation passed.

The literal planned 10,000-shot, 50-point refinement timed at approximately
143.8 single-core hours.  A profile placed 96.7% of runtime in final BP+OSD,
and the installed decoder exposes no batch-decode interface.  The completed
1,000-shot grid retains all five sizes and seven critical probabilities and
reports its wider Monte Carlo intervals directly.

At the central point `p=0.015`, the state-failure rates are:

| L | Final-only, 95% Wilson interval | One-step feedback, 95% Wilson interval | Paired feedback minus final-only |
|---:|---:|---:|---:|
| 5 | 0.163 [0.141, 0.187] | 0.222 [0.197, 0.249] | +0.059 |
| 7 | 0.287 [0.260, 0.316] | 0.248 [0.222, 0.276] | -0.039 |
| 9 | 0.399 [0.369, 0.430] | 0.233 [0.208, 0.260] | -0.166 |
| 11 | 0.464 [0.433, 0.495] | 0.211 [0.187, 0.237] | -0.253 |
| 13 | 0.591 [0.560, 0.621] | 0.226 [0.201, 0.253] | -0.365 |

Feedback adds overhead at `L=5`, but its benefit grows with size.  That change
of size dependence is what creates the feedback crossing; a pointwise
feedback improvement by itself is not a threshold estimate.

## Scaling analysis

The fitter uses independent binomial likelihoods for the two protocols and

\[
P_{\rm fail}(p,L)=\operatorname{logistic}
\left(\sum_{j=0}^{q}a_j x^j\right),\qquad
x=\frac{p-p_c}{s_p}(L/L_{\rm ref})^{1/\nu}.
\]

It withholds `p_c` for missing adjacent-size crossings, optimizer-bound hits,
nonmonotone fitted curves, or a deviance p-value below 0.01.

The main diagnostics are:

| Model | Status | Candidate or accepted `p_c` | Deviance / degrees of freedom | Deviance p-value |
|---|---|---:|---:|---:|
| Final-only, all sizes | no crossing | unavailable | unavailable | unavailable |
| Quadratic feedback, all sizes | rejected | 0.01494 | 66.79 / 30 | 0.000129 |
| Cubic feedback, all sizes | rejected | 0.01470 | 50.30 / 29 | 0.00838 |
| Cubic feedback, `L>=7` | accepted | **0.01482** | 37.67 / 22 | 0.01996 |

For the accepted model, `nu=1.056`; its conditional 95% bootstrap interval is
`[0.914,1.239]`.  This finite-size study does not establish a universal
critical exponent.  The threshold interval is conditional on the selected
scaling form and excludes decoder-model systematics.

## Routing sensitivity

At `L=13,p=0.015`, the 200-shot full-path run gives 40 failures, or 0.200,
while the one-step production point gives 226/1000, or 0.226.  The observed
difference is 0.026.  Even 10,000 full-path shots would have an estimated 95%
binomial half-width near 0.0078, above one quarter of the difference (0.0065).
None of the prespecified 2,000, 5,000, or 10,000 counts satisfies that precision
rule.  The saved 200-shot result is therefore exploratory and does not support
a quantitative claim that full-path routing is better.

## Relationship to the paper

The construction follows the qualitative mechanism in *Noise-resilient
sequential circuits for generating quantum order*: advance through diagonal
two-dimensional slices, measure slice stabilizers, pair point defects, and
feed corrections forward before those defects propagate farther.  The
explicit honeycomb boundary, seven-gate transfer, matching tie breaks,
post-gate-output noise, and final BP+OSD decoder are the repository's concrete
choices.  This is a simulation of that derived circuit under the retained
noise convention.  It is not a numerical reproduction of a threshold quoted
by the paper, and it does not prove the paper's local-channel stability
statement.

## Reproduction and artifacts

The public Julia API is `compare_toric3d_x_feedback_threshold`; the runnable
scan is [`scan_toric3d_x_feedback_threshold.jl`](../../examples/scan_toric3d_x_feedback_threshold.jl),
the pooler is [`assemble_toric3d_x_feedback_production.py`](../../examples/assemble_toric3d_x_feedback_production.py),
and the guarded fitter is
[`fit_toric3d_x_feedback_threshold.py`](../../examples/fit_toric3d_x_feedback_threshold.py).

```sh
JULIA_CONDAPKG_BACKEND=Null \
JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" \
MPLCONFIGDIR=/tmp/toponoise-mpl \
julia --project=. examples/scan_toric3d_x_feedback_threshold.jl --help

PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/toponoise-mpl \
.CondaPkg/.pixi/envs/default/bin/python \
  examples/fit_toric3d_x_feedback_threshold.py \
  results/toric3d-x-feedback-threshold-2026-10-11/production/toric3d-x-feedback-raw.csv \
  --full-path-csv results/toric3d-x-feedback-threshold-2026-10-11/full-path/toric3d-x-feedback-raw.csv \
  --p-min .01 --p-max .025 --min-size 7 --degree 3 \
  --bootstrap 1000 --seed 1234 \
  --output-dir results/toric3d-x-feedback-threshold-2026-10-11/analysis-cubic-L7
```

Production artifacts:

- [raw paired CSV](../../results/toric3d-x-feedback-threshold-2026-10-11/production/toric3d-x-feedback-raw.csv)
- [full provenance metadata](../../results/toric3d-x-feedback-threshold-2026-10-11/production/toric3d-x-feedback-metadata.json)
- [scan-selection and runtime record](../../results/toric3d-x-feedback-threshold-2026-10-11/production/scan-selection.toml)
- production curves as [PDF](../../results/toric3d-x-feedback-threshold-2026-10-11/production/toric3d-x-feedback.pdf),
  [SVG](../../results/toric3d-x-feedback-threshold-2026-10-11/production/toric3d-x-feedback.svg), and
  [PNG](../../results/toric3d-x-feedback-threshold-2026-10-11/production/toric3d-x-feedback.png)

Accepted analysis:

- [fit and bootstrap report](../../results/toric3d-x-feedback-threshold-2026-10-11/analysis-cubic-L7/threshold-fit.json)
- [bootstrap samples](../../results/toric3d-x-feedback-threshold-2026-10-11/analysis-cubic-L7/threshold-bootstrap.csv)
- [size/window sensitivity table](../../results/toric3d-x-feedback-threshold-2026-10-11/analysis-cubic-L7/threshold-sensitivity.csv)
- fit figures as [PDF](../../results/toric3d-x-feedback-threshold-2026-10-11/analysis-cubic-L7/threshold-fit.pdf),
  [SVG](../../results/toric3d-x-feedback-threshold-2026-10-11/analysis-cubic-L7/threshold-fit.svg), and
  [PNG](../../results/toric3d-x-feedback-threshold-2026-10-11/analysis-cubic-L7/threshold-fit.png)

Guarded rejections are preserved in the
[quadratic report](../../results/toric3d-x-feedback-threshold-2026-10-11/analysis-quadratic/threshold-fit.json)
and [all-size cubic report](../../results/toric3d-x-feedback-threshold-2026-10-11/analysis-cubic/threshold-fit.json).
The routing sensitivity has its own
[raw CSV](../../results/toric3d-x-feedback-threshold-2026-10-11/full-path/toric3d-x-feedback-raw.csv)
and [metadata](../../results/toric3d-x-feedback-threshold-2026-10-11/full-path/toric3d-x-feedback-metadata.toml).

## Limits of the conclusion

The estimate applies to pure X faults, ideal measurements and corrections,
the encoded logical-`|00>` boundary, one-step matching feedback, and the
specified final BP+OSD decoder.  Uniform data priors neglect nonuniform gate
exposure, CNOT-induced correlations, and correlations created by feedback.
The study omits measurement faults, correction faults, idle faults, combined
X/Z noise, circuit-level maximum-likelihood decoding, and a proof of locality
for the classical computation.  The 1,000-shot production design also leaves
visible statistical fluctuations in individual adjacent-size crossings.

Within those assumptions, the result answers the original question: the
diagonal encoded circuit with within-slice X-syndrome feedback shows a
nonzero logical-X preparation threshold near 1.5%, while the same circuit
with final-only decoding does not show a positive crossing in the simulated
window.
