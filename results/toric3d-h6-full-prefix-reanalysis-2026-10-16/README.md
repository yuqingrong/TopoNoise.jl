# H+6 full-prefix feedback: uncertainty-aware follow-up

This is an exploratory reanalysis of the existing 70,000 production fault
records, not an additional circuit simulation. Each channel has five sizes
L=3,5,7,9,11, seven physical fault probabilities, T=L, and 1,000 records per
point. Only the full-prefix arm enters these fits. All original simulation,
fit, and data-only figure files are preserved.

The analysis choices are recorded in `analysis-plan.md`, fixed before the new
fits but after seeing the original results. `analysis-addendum.md` identifies
one later extension: examining all ten size pairs using the same simultaneous
confidence bounds. The secondary extension is explicitly exploratory.

## X: the data are compatible with the common scaling model

The original common quadratic-logit fit is reproduced exactly:

- Candidate p_c = 0.05376756 (5.3768%); nu = 1.45184.
- Conditional 95% bootstrap interval for p_c: [0.05259380, 0.05502110].
- Conditional 95% bootstrap interval for nu: [1.26109, 1.66321].
- Parametric-bootstrap goodness-of-fit p-value: 0.9680.
- Bootstrap likelihood-ratio p-value against a separate quadratic-logit curve
  for every size: 0.6973. The observed likelihood-ratio statistic is 7.37727.

All 1,000 bootstrap datasets were refitted successfully in both models; none
hit a parameter bound or produced a nonmonotone common fit. The Monte Carlo
95% intervals for the two tail frequencies are [0.9551,0.9780] and
[0.6675,0.7254], respectively. These intervals quantify bootstrap simulation
error, not uncertainty in the physical threshold.

The old raw-sign crossing rule rejected 918 of these 1,000 datasets (91.8%),
even though each dataset was generated from the fitted common-scaling model.
Its rejection-frequency Monte Carlo 95% interval is [89.92%,93.43%]. This
demonstrates brittleness under this fitted model and sampling design; it is
not a universal false-rejection rate for every possible code/noise model.

The uncertainty-aware crossing check uses exact binomial intervals with 95%
simultaneous coverage across the 35 X probabilities. Adjacent-size contrasts
alone are compatible with zero throughout the whole measured window, so they
do not localize a crossing. Including all ten size pairs gives the common
compatibility region [0.04618863,0.06066628], conditional on linear interpolation
between measured probabilities. This is broad and contains the fitted p_c.
It neither proves a unique crossing nor an asymptotic threshold.

Thus about 5.38% is a statistically compatible **finite-size model estimate**.
Its narrow bootstrap interval conditions on the quadratic scaling model and
does not include uncertainty from the choice of sizes or correction models.
The original saved acceptance status is not rewritten.

## Z: smallest-size corrections matter, with remaining model uncertainty

All 36 planned combinations are retained: minimum size 3,5,7; full and narrow
rate windows; and six leading/correction models. The full window is
[0.025,0.040], and the narrow window is [0.0275,0.0375]. The leading models
are quadratic or cubic in the scaling variable. The correction models add
terms proportional to L^(-1/2), L^(-1), L^(-2), or an L^(-1) offset and slope
in logit space. The correction exponent is fixed in each model.

Representative full-window results:

| Included sizes | Model | Candidate p_c | Goodness-of-fit p-value |
| --- | --- | ---: | ---: |
| L >= 3 | Leading quadratic | 0.03359795 | 0.00493 |
| L >= 3 | Quadratic + L^(-1) correction | 0.03641195 | 0.10464 |
| L >= 5 | Leading quadratic | 0.03433973 | 0.15623 |
| L >= 5 | Quadratic + L^(-1) correction | 0.03451315 | 0.12464 |

The L >= 3 leading fit has poor agreement. Explicit size-dependent corrections
improve it, whereas simply making the common curve cubic does not. Dropping
L=3 substantially reduces the discrepancy: all six full-window L >= 5 models
fit adequately, with candidate p_c values around 3.43%--3.48%. In the narrower
window, those candidates span about 3.39%--3.57%. These ranges describe model
sensitivity and are **not confidence intervals**.

The leading quadratic has the smallest descriptive AICc among the six models
for L >= 5 in both windows. The extra correction is therefore not clearly
required by that subset. This does not prove that corrections vanish. Models
using only L=7,9,11 are retained as three-size diagnostics and show greater
instability; they are not independent threshold determinations.

The four full-window profile intervals show how much the correction can change
parameter precision even when its best-fit value changes little:

| Included sizes | Model | Conditional 95% profile interval for p_c |
| --- | --- | --- |
| L >= 3 | Leading quadratic (poor fit) | [0.03312360, 0.03409923] |
| L >= 3 | Quadratic + L^(-1) correction | [0.03476279, 0.03854570] |
| L >= 5 | Leading quadratic | [0.03367279, 0.03509221] |
| L >= 5 | Quadratic + L^(-1) correction | [0.03198904, 0.03795083] |

The first interval is diagnostic because that model fails goodness of fit.
The L >= 5 corrected interval is much broader despite a similar best-fit p_c:
the data do not tightly separate the correction amplitude from the critical
probability. Around 3.4% is a useful working value, but a precise estimate
independent of correction-model assumptions is not supported. All intervals
are conditional on their individual models and bounds; none is an uncertainty
interval covering all finite-size model choices.

## Verification

All 36 Z fits converged without parameter-bound or monotonicity failures.
All four final profiles have one detected connected interval, no final
numerical gaps, and minima consistent with the unrestricted fits. Failed
internal optimizer starts/polishes are retained in the logs; they are not
silently interpreted as disconnected statistical regions.

Ten X tests and twelve Z tests passed. Checks include analytic gradients,
synthetic parameter recovery, exact binomial intervals, original-fit
reproduction, nesting, and likelihood-profile endpoints. The independent
review regenerated all 1,000 X bootstrap count sets and checked the rejection
counts and intervals. The unchanged circuit simulation was not rerun.

## Files and reproduction

`input.json` contains the validated counts and original fit records. The input
hash and initial plan hash are in `plan-freeze.json`. X numerical results are
in `x-report.json`, `x-cells.csv`, `x-pair-bands.csv`, `x-crossing-regions.csv`
and `x-bootstrap.csv`; the last retains the generated binomial counts for all
bootstrap replicates. Z results are in `z-report.json`, `z-fits.csv`,
`z-predictions.csv` and the four `z-profile-*.csv` files.

`methods-notes.md` explains interval coverage, refitting, normalization,
profile-likelihood assumptions and the descriptive use of AICc. All PDFs in
`figures/` use the no-grid theme. These are diagnostic figures; the prior
data-only PDFs remain in their separate results folder.

For a clean numerical reproduction, restore `reproduction-source.tar.gz`
into a fresh compatible project root and follow `reproduce.txt`, using the
versions in `reproduction-environment.json`. The archive retains the original
relative source paths required by the X likelihood engine and plotting theme.
Its source hashes are saved separately. The follow-up input is self-contained
for fitting, while its provenance binds the prior raw counts and validated
simulation bundle. No original circuit simulation is needed for this reanalysis.

The exported follow-up is located in the normal project at:
`/Users/rongyuqing/jcode/TopoNoise.jl/results/toric3d-h6-full-prefix-reanalysis-2026-10-16/`.
