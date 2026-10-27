# H+6 full-prefix feedback: targeted independent simulation

Status: completed and verified. All 204,000 fresh records finished. The
prespecified X anchor check resolves a size-order reversal; Z still has
finite-size correction and fit-window ambiguity. Production took about
98 minutes, and the complete numerical analysis took 70.17 seconds.

This follow-up addresses the sampling uncertainty in the X crossing and the
finite-size correction ambiguity in Z found in the prior reanalysis. The
user selected a balanced local run. The design and sample counts were fixed
before production in `selection.json`; its checksum and orchestration files
are bound by `setup-freeze.json`.

## Fixed sampling design

| Channel | L, with T=L | Physical gate-output error probabilities | Fresh records per point |
| --- | --- | --- | --- |
| X | 5 | .0475, .0500, .0525, .0550, .0575, .0600, .0625 | 4,000; 5,000 at .0500 and .0600 |
| X | 7, 9 | same X grid | 4,000 |
| X | 11 | same X grid | 2,000; 5,000 at .0500 and .0600 |
| Z | 5, 7, 9, 11, 13, 15, 17 | .0275, .0300, .0325, .0340, .0350, .0375, .0400 | 2,000 |

There are 28 fresh X cells and 49 fresh Z cells, containing 106,000 X and
98,000 Z records: 204,000 independent new fault records. Each record is
replayed through both final-only and full-prefix arms, yielding 408,000
paired outcomes. The full-prefix arm has 1,896,000 required feedback steps.
The two outcomes from one record are paired, not independent samples.

Two 20-record timing benchmarks are excluded from all estimates. Their
counts are not used to choose the rate window or alter sample allocations.
All new point seeds were checked against 315 previously saved distinct
streams before production. Original results and simulation sources remain
unchanged. The original paired runner is reused with its exact protocol,
source, runtime, checksum, and resume guards.

## Physical protocol

After every transfer r=1,...,T, including T, recover the entire available
prefix Q0 through Qr and verify that every required completed syndrome is
zero. Future slices are untouched. Score one final full-code outcome per
arm. X begins in ideal encoded logical |00>; Z has noisy seed-H preparation
of logical |++> and no separate r=0 repair. Measurements, corrections and
idles are ideal. A failure means either of the two residual logical parities
is nonzero. X and Z are separate noise channels. Existing BP+OSD (X) and
unit-weight MWPM (Z) settings are retained.

## Analysis fixed before observing production outcomes

The primary analysis uses only fresh production records. Old-plus-new pooling
is secondary, restricted to the same declared fresh cells and independent
original production streams. Pilots, timing benchmarks, replay checks, and
archived duplicate runs never contribute extra samples.

For X, the primary crossing check is the L=11 minus L=5 contrast at p=.050
and p=.060. Exact binomial intervals with Bonferroni correction over those
four cells provide simultaneous 95% coverage. The desired sign pattern is
negative at .050 and positive at .060. All-cell bands and adjacent-size
contrasts remain visible secondary diagnostics; there is no strict raw-sign
crossing rejection rule. The leading fit receives 1,000 parametric bootstrap
refits. Six scaling/correction models are examined in windows [.0475,.0625]
and [.050,.060], at minimum sizes 5 and 7. The latter has only three sizes
and is diagnostic.

For Z, the same six models are examined in windows [.0275,.0400] and
[.0300,.0375], at minimum sizes 5,7,9. The largest minimum retains five
sizes. Leading quadratic and inverse-L correction fits receive profile
likelihood intervals in every window/subset. All fits, convergence flags,
boundary hits and failed checks are retained. Model comparisons use identical
data subsets; AICc is descriptive. All intervals condition on their model.

The stopping rule was completion of the fixed counts, never a desired
crossing or p-value. A clearer X sign bracket and narrower Z corrected
intervals stable across size/window changes were the planned objectives.
All remaining drift, broad intervals and boundary limitations are reported.
No outcome-driven extension of the size/rate grid or sample counts was made.
No finite collection of successful fits proves an asymptotic threshold.

## X: the prespecified sampling check succeeds

For the full fresh window [.0475,.0625] and L=5,7,9,11, the leading
quadratic-logit estimate is **p_c = 0.05342984 (5.343%)**, with a conditional
95% bootstrap interval **[0.05248217, 0.05443343] (5.248%--5.443%)**.
The exponent estimate is nu=1.41796, with conditional interval
[1.19709,1.72542]. Bootstrap goodness-of-fit p=0.53047 and the bootstrap
likelihood-ratio p=0.68332 against independent quadratic curves for each
size show no detected incompatibility with this leading model. All 1,000
bootstrap datasets were successfully refitted, with no bound hits or
nonmonotone common refits.

The primary contrast is the larger-size failure probability minus the
smaller-size probability. Its intervals have simultaneous 95% coverage
from exact binomial intervals over the four predeclared cells:

| Physical p | L=11 minus L=5 failure probability | Simultaneous difference interval |
| --- | ---: | --- |
| .050 | -0.027800 | [-0.055284, -0.00026549] |
| .060 | +0.062000 | [+0.028178, +0.095709] |

Both required signs are resolved. The .050 interval only narrowly excludes
zero, so the lower-rate evidence should not be described as overwhelming.
This check establishes the planned separated-size ordering reversal, not a
unique common crossing or an asymptotic threshold.

The secondary simultaneous adjacent-size bands remain compatible with zero
throughout the full [.0475,.0625] window. The all-pair common compatibility
region is [.04754951,.05873118], conditional on linear interpolation between
sampled probabilities. The L>=5 model/window candidates span approximately
5.301%--5.544%; that span is sensitivity evidence, not a confidence interval.
Thus approximately 5.34% is a supported leading-model finite-size estimate,
with remaining model uncertainty beyond its bootstrap interval.

## Z: larger sizes expose remaining correction and window limits

For the fresh full window [.0275,.0400], the leading estimate increases as
small sizes are removed. At L>=5, its goodness-of-fit p=0.01041 is weak;
adding a 1/L correction gives p=0.52379 and improves descriptive AICc by
24.13 on the same data. The larger size range therefore does not support
treating the earlier roughly 3.4% working estimate as settled.

The table reports percentages. Every interval is conditional on its listed
model. **UPPER-TRUNCATED** means the profile region reaches the imposed
p_c=4.0000% boundary; 4.0000% is not a resolved upper confidence endpoint.

| Included sizes | Leading p_c | Leading 95% profile interval | +1/L p_c | +1/L 95% profile region |
| --- | ---: | --- | ---: | --- |
| L>=5 | 3.6852% | [3.6529%, 3.7207%] | 3.9436% | [3.8245%, 4.0000%], **UPPER-TRUNCATED** |
| L>=7 | 3.7369% | [3.6912%, 3.7898%] | 3.9404% | [3.7577%, 4.0000%], **UPPER-TRUNCATED** |
| L>=9 | 3.8060% | [3.7342%, 3.8965%] | 3.7842% | [3.5221%, 4.0000%], **UPPER-TRUNCATED** |

The similar leading and corrected point estimates at L>=9 do not resolve
the correction uncertainty: the corrected profile remains broad and reaches
the upper boundary. In the narrow [.0300,.0375] window, the L>=5 and L>=7
corrected estimates reach p_c=.0375. All three narrow corrected profiles
reach that boundary, as does the L>=9 narrow leading profile. In total,
**all six corrected profiles and the L>=9 narrow leading profile are
upper-truncated**.

The expanded data reveal size drift and a rate window that is too limited
to close the corrected upper profiles. Full-window candidates of roughly
3.7%--4.0% are useful descriptive working values, not a confidence interval
or a settled Z threshold. A future independently frozen Z study needs
sampling above p=.04 to examine this upper-window limitation, together with
the same size/correction sensitivity checks. No such extension was added
after observing these outcomes.

## Numerical completion, model flags and prior comparisons

The numerical pipeline completed all 120 planned fits: 24 X and 36 Z fits
for fresh data, and the same schedules for secondary pooling. All 12 fresh
Z profiles and all 1,000 fresh X bootstrap replicates completed without
numerical failures, missing results or profile-minimum disagreement.
The report's `complete` status describes this numerical completion; it is
not acceptance of every physical scaling model or either asymptotic threshold.

There are 28 fits at or near parameter bounds, all in Z: 10 fresh and
18 pooled. Seven fresh Z profile regions are upper-truncated as described
above. Four goodness-of-fit p-values are below .01; all four are the pooled
Z L>=5 leading quadratic/cubic fits in the full and narrow windows. No fit
is nonmonotone. These flags remain in `analysis/report.json` and are kept
separate from numerical convergence.

Secondary pooling adds 36,000 matching original production records to the
204,000 fresh records. Its leading full-window X estimate is 5.369%; the
L>=5 full-window corrected Z estimate is 3.978%. These results support the
same qualitative interpretation, but pooling shares the fresh records and
is not independent confirmation.

The prior saved X primary bootstrap included L=3 and used [.035,.065],
whereas the new primary analysis excludes L=3 and uses [.0475,.0625].
Comparing their interval widths does not isolate a gain from additional
shots. The old Z windows were [.025,.040] and [.0275,.0375], and its largest
size was L=11. New Z sizes, rates and windows differ. In particular, a
shorter new corrected profile span cut off at a boundary must not be
presented as resolved uncertainty or an equivalent precision improvement.

## Verification and figures

The completed production audit verified 204,000 records across 13 shards,
77 cells and 2,088 immutable snapshots, including 1,896,000 required
full-prefix feedback steps. Independent replay checks reproduced 4,000 X
records with 336,201 assertions and 100 Z records with 21,327 assertions.
The independent numerical audit passed 56,144 assertions. Four focused test
suites passed 33, 6, 5 and 8 tests, respectively. Original simulation and
analysis artifacts and `document/note.tex` remain unchanged; preservation
evidence is saved in `verification/prior-artifact-preservation.json`.

Five PDFs were generated and passed rendered visual and numerical QA,
recorded in `verification/pdf-qa.json`:

- `figures/x-data.pdf` and `figures/z-data.pdf`: fresh measured failure
  probabilities, pointwise 95% Wilson intervals and straight connecting
  segments only; no fitted curves, fit annotations or insets.
- `figures/x-crossing.pdf`, `figures/z-size-sensitivity.pdf` and
  `figures/z-profile.pdf`: separate crossing, finite-size and likelihood
  diagnostics, with boundary limitations shown explicitly.

Every figure uses the established no-grid theme. The five PDFs were checked
page by page for clipping, overlaps, legibility and embedded fonts; all 77
plotted data means and Wilson intervals were independently verified.

## Files and reproduction

`jobs.toml` assigns 13 disjoint production shards to two workers.
`worker-1.log` and `worker-2.log` show job completion. Each `production-*`
directory retains `run.log`, raw counts, runtime/protocol/source metadata,
and immutable batch checkpoints. All production rows reached their target
counts. `analysis/input.json` binds the validated fresh and pooled cells;
`analysis/report.json` retains the complete numerical results and provenance.
CSV files retain fit predictions, profile evaluations and bootstrap counts.
`reproduce.txt` documents reproduction and its source/data prerequisites.

Delivery directory:
`/Users/rongyuqing/jcode/TopoNoise.jl/results/toric3d-h6-full-prefix-targeted-2026-10-16/`.
