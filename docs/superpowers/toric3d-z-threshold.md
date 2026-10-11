# Logical-Z threshold of diagonal-slice preparation

This study measures logical failure after recovery. It retains independent
rate-p Z faults on each participating qubit after every H or CNOT, ideal idle
qubits, measurements, and Pauli corrections. The diagonal transfer gates are
unchanged. Its initial boundary is now an encoded honeycomb toric-code slice,
prepared with noisy H gates and ideal stabilizer measurements. This boundary
choice replaces the earlier fixed-product, k=0 boundary; it is essential to
defining a logical threshold.

The result concerns the Z sector and the specified final MWPM decoder. Under
pure Z noise the active-slice Z compatibility measurements see no defects,
so the two earlier feedback protocols coincide exactly. This study does not
claim a threshold for both noise sectors or numerical reproduction of the
paper's noise model.

The completed study gives a Z-only threshold estimate of approximately
**p_c,Z = 0.0052 (0.52% per participating gate-output qubit)**. This applies
to the encoded diagonal-slice circuit defined here. The earlier rough-slab
circuit still has its logical-Z obstruction: a change of circuit and boundary,
with the same gate-noise convention, is what makes this new result possible.
Pure-Z feedback and final-only recovery coincide for the new circuit because
the slice compatibility measurements detect X defects. Its Z stability comes
from confining fault propagation, not from those measurements correcting Z.

## Checks and logical qubits

Use the same coordinates and edges as the
[slice-transfer derivation](toric3d-slice-feedback.md). The initial slice
contains 3M edges, M=LyLz, forming a periodic honeycomb lattice. Its X stars
are black triads E[0,d,r] and white triads E[0,d,r+s_d]. Its Z checks are the
six-edge compatibility loops Q. The seed code has two logical qubits.

The final X checks include the lower white triads, six-edge stars at times
0 through T−1, and upper triads at T. Z checks are the existing four-edge
physical plaquettes. Every data edge touches exactly two X checks. The ranks
and distances are

\[
n=3(T+1)M,\qquad r_X=(T+2)M-1,\qquad
r_Z=(2T+1)M-1,\qquad k=2,
\]
\[
d_Z=2\min(L_y,L_z),\qquad d_X=(T+1)\min(L_y,L_z).
\]

Logical Z_u is the alternating direction-1/direction-3 cycle along u on a
fixed seed row; Z_v is the analogous direction-2/direction-3 cycle along v.
Logical X_u consists of all direction-1 edges ending at u=0, through all
slices; X_v consists of all direction-2 edges ending at v=0. These membrane
parities distinguish the two residual Z windings after syndrome removal.

A closed Z path must take equal numbers of forward and backward t steps.
Winding once in u or v therefore requires at least 2Ly or 2Lz edges; the
displayed seed cycles attain those bounds. The X distance follows from the
disjoint parallel logical-Z cycles, one per slice and transverse row. Tests
also compute shortest winding cycles in a four-sheet graph lift, including
periodic size two.

## Noisy preparation of the encoded seed

Start with product zeros. Apply H to each of the 3M seed qubits, including an
independent Z-p fault after each H. Ideally measure the Q checks and apply an
X recovery to their measurement outcomes. The all-plus input and these X
recoveries preserve both logical-X eigenvalues, preparing logical |++>.

Sampling the ideal projection outcomes is unnecessary. If P_s projects onto
the measured Q syndrome s, X^c maps s to zero, and Z^e is the seed fault,

\[
X^c P_s Z^e |+\cdots+\rangle
\propto Z^e P_0 |+\cdots+\rangle.
\]

Here Z^e commutes with every Q, and X^c fixes the all-plus input; commuting
the two Paulis changes only a global phase. Thus independent seed Z-p faults
are an exact frame representation of this noisy-H preparation, without
assuming a perfectly prepared entangled seed. The subsequent transfer adds
the same seven gates and thirteen noisy active outputs per vertex as before.

## Exact accelerated sampling and recovery

Each independent fault has a sparse column recording its final Z support.
There are 3M seed fault locations and 13TM transfer fault locations. A column
has at most three qubits in two adjacent slices. For a shot, independently
sample all fault bits f and compute the final data frame e=Ff over GF(2).
The sparse operation is exactly equivalent to explicit gate-by-gate Pauli
propagation and preserves the correlations induced by CNOTs. Exhaustive
small-system tests compare every column to an independent explicit replay.

Compute the measured star syndrome H_X e. The decoder receives only this
syndrome and the declared lattice. It uses unit-weight data-edge PyMatching
MWPM and returns a full data recovery r. The scan stops if H_X r differs
from H_X e. Score the two parities X_u·(e+r), X_v·(e+r). A shot fails if either
is nonzero; both individual counts are also saved. For this Z-only |++>
preparation, this event is exactly residual logical state failure.

Unit data weights deliberately retain the baseline decoder. They neglect
known short-range circuit correlations. A circuit-informed matching graph
could improve the result; the present threshold is decoder dependent.

## Scan and fitting

The configurable Julia entry point accepts independent (T,Ly,Lz) dimensions,
a caller-owned RNG, probabilities, shot count, and batch size:

```julia
using Random, TopoNoise
scan = scan_toric3d_z_threshold(MersenneTwister(1234);
    sizes=[5,7,9], error_rates=[.002,.004,.006,.008,.01,.015,.02],
    shots=1000, seed=1234, output_dir="results/z-pilot", progress_io=stdout)
```

The CLI `examples/scan_toric3d_z_threshold.jl` runs a separate warmup and
timing pilot before production. CSV/TOML files are saved after every completed
point. The Python analysis example accepts an explicit fitting window and
minimum size. It uses binomial likelihood with

\[
P_{\rm fail}(p,L)=\operatorname{logistic}(a_0+a_1x+a_2x^2),
\qquad x=\frac{p-p_c}{s_p}(L/L_{\rm ref})^{1/\nu}.
\]

The normalizations s_p and L_ref improve numerical conditioning and do not
change p_c or nu. A cubic alternative is available. Wilson intervals use
independent shot counts. Parametric binomial bootstrap intervals quantify
Monte Carlo uncertainty conditional on the scaling ansatz; omitting smaller
sizes and narrowing the fitting window probes systematic finite-size drift.
Missing crossings, parameter-bound hits, nonmonotone fits, and a poor
deviance goodness-of-fit are reported rather than hidden.

```sh
julia --project=. examples/scan_toric3d_z_threshold.jl --help
.CondaPkg/.pixi/envs/default/bin/python examples/fit_toric3d_z_threshold.py --help
```

On this machine Julia uses the existing decoder environment via
`JULIA_CONDAPKG_BACKEND=Null` and
`JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python"`.

## Completed numerical study, 2026-10-11

The exploratory scan used L=5,7,9, seven probabilities, and 1,000 shots per
point (21,000 shots). These pilot counts are excluded from threshold fitting.
Production used T=Ly=Lz=L, 10,000 independent shots per point:

| Sizes | Gate-output probabilities | Seed | Shots |
|---|---|---|---:|
| 9,13,17,21 | .004,.0045,.005,.0055,.006,.0065,.007 | 1234 | 280,000 |
| 25 | .0045,.005,.0055,.006,.0065 | 1235 | 50,000 |

All 33 points were saved incrementally. The combined file contains 330,000
production shots. The seed and transfer H/CNOT faults use the same p.
Final unit-weight PyMatching 2.4.0 removes the measured star syndrome; every
full recovery is checked before either logical winding is scored.

The observed change of size dependence is substantial:

| p | L=9 logical failure, 95% Wilson interval | L=25 logical failure, 95% Wilson interval |
|---|---|---|
| .005 | 0.850% [0.688%,1.050%] | 0.320% [0.227%,0.451%] |
| .006 | 3.180% [2.854%,3.542%] | 10.280% [9.700%,10.891%] |

An uncorrected quadratic fit over .0045<=p<=.0065 drifts from p_c=.005374
to .005276 to .005167 as the minimum size increases from 9 to 13 to 17.
To examine this drift, the supplemental fitter adds a leading size correction:

\[
\operatorname{logit}P_{\rm fail}
=a_0+a_1x+a_2x^2+b(L_{\rm ref}/L)^\omega.
\]

The primary corrected fit fixes omega=1 and uses all five sizes in this
window. Its estimate is p_c=.005171, nu=.944, with binomial deviance
17.72 for 19 degrees of freedom (goodness-of-fit p=.541). Excluding the
smallest sizes gives .005155 and .005153. The correction exponent is not
determined by these data: omega=.5 and 2 give similar fit quality and
estimates .005100 and .005236. Trimming the fit window introduces further
drift. The examined corrected and uncorrected fits give estimates roughly
0.49%–0.54%. This is a model-sensitivity range, not a 95% confidence interval
or a guarantee on the infinite-size value. The appropriate rounded numerical
conclusion is about 0.52%, with finite-size uncertainty retained explicitly.
Of 1,000 parametric bootstrap refits, 999 passed the numerical guards; the
conditional 95% interval is p_c=[.00504584,.00529814], or
**0.5046%–0.5298%**. The profile-likelihood interval [.00504153,.00528963]
agrees. These intervals measure sampling uncertainty under the fixed
correction model; they do not include the model-sensitivity range above.
Neither the fitted nu nor the assumed correction exponent is established as
a universal critical exponent by this finite-size study.

The [corrected fit](../../results/toric3d-z-threshold-2026-10-11/corrected-analysis/corrected-fit.json)
records the primary estimate, conditional bootstrap and profile-likelihood
intervals, numerical diagnostics, and input/source hashes. Its
[sensitivity table](../../results/toric3d-z-threshold-2026-10-11/corrected-analysis/corrected-sensitivity.csv)
keeps the model/window comparisons separate. Figures are available as
[PDF](../../results/toric3d-z-threshold-2026-10-11/corrected-analysis/corrected-fit.pdf),
[SVG](../../results/toric3d-z-threshold-2026-10-11/corrected-analysis/corrected-fit.svg), and
[PNG](../../results/toric3d-z-threshold-2026-10-11/corrected-analysis/corrected-fit.png).
The [uncorrected fit](../../results/toric3d-z-threshold-2026-10-11/analysis/threshold-fit.json)
and its
[minimum-size/window audit](../../results/toric3d-z-threshold-2026-10-11/analysis/minimum-size-sensitivity.json)
provide a separate comparison. Raw adjacent-size crossings are noisy and
drift with size, so no individual crossing is treated as the threshold.

The two fit commands, run from the repository root, are:

```sh
.CondaPkg/.pixi/envs/default/bin/python examples/fit_toric3d_z_corrections.py \
  results/toric3d-z-threshold-2026-10-11/combined/toric3d-z-threshold-raw.csv \
  --p-min .0045 --p-max .0065 --min-size 9 --omega 1 \
  --bootstrap 1000 --seed 1234 \
  --output-dir results/toric3d-z-threshold-2026-10-11/corrected-analysis
.CondaPkg/.pixi/envs/default/bin/python examples/fit_toric3d_z_threshold.py \
  results/toric3d-z-threshold-2026-10-11/combined/toric3d-z-threshold-raw.csv \
  --p-min .0045 --p-max .0065 --min-size 13 --bootstrap 1000 --seed 1234 \
  --output-dir results/toric3d-z-threshold-2026-10-11/analysis
```

Reproduce the production scans with the Julia environment variables described
above, then combine the disjoint size sets:

```sh
julia --project=. examples/scan_toric3d_z_threshold.jl \
  --sizes 9,13,17,21 --error-rates .004,.0045,.005,.0055,.006,.0065,.007 \
  --shots 10000 --seed 1234 --output-dir results/z-refined
julia --project=. examples/scan_toric3d_z_threshold.jl \
  --sizes 25 --error-rates .0045,.005,.0055,.006,.0065 \
  --shots 10000 --seed 1235 --output-dir results/z-largest
```

The archived combined CSV, source hashes, and metadata are in
[the combined results directory](../../results/toric3d-z-threshold-2026-10-11/combined/).
The first timing estimates included lazy decoder initialization; their TOML
files record this limitation. The example now times initialization separately.
Recorded per-point runtimes are the actual production measurements.

Verification includes signed binary stabilizer checks at every preparation
prefix, ranks and logical anticommutation including periodic size two,
shortest nontrivial winding cycles, exhaustive agreement of every eligible
small-system fault with explicit propagation, valid single-fault recovery,
noiseless success, seeded reproducibility, and invalid-recovery rejection.
The entire Julia regression suite, including the existing 2D interfaces,
passed. The new encoded-geometry tests contribute 135 assertions and the
fault/decoder/scan tests contribute 376 assertions. The uncorrected fitter's
eight statistical tests and the corrected fitter's four tests passed. The
latter check the analytic Jacobian, synthetic parameter recovery, rejection
at parameter bounds, and bootstrap reproducibility. Independent numerical
reviews verified both likelihood gradients and reproduced the main estimate.
The updated LaTeX note compiles; its pre-existing duplicate-label and
unrelated overfull-box warnings remain.
