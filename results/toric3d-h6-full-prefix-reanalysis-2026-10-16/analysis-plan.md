# Exploratory follow-up of H+6 full-prefix results

This follow-up uses the same 70,000 production fault records. It is exploratory:
the original fit diagnostics were already seen. No additional circuit samples,
changes to noise/decoding, or changes to the original acceptance policy are made.
All original artifacts remain unchanged. The validated snapshot is `input.json`.

## X: uncertainty-aware assessment

Use all five sizes and the original seven probabilities. Report exact
Clopper-Pearson intervals with Bonferroni coverage across all 35 binomial
probabilities; subtract endpoints to bound every adjacent-size difference.
Piecewise-linear interpolation of those bounds supplies crossing-compatible
regions. Keep every connected component, including missing/unbounded-within-grid
cases. This is conditional on linear interpolation between measured rates;
compatibility does not demonstrate a unique crossing or an asymptotic threshold.

Retain the original quadratic-logit scaling model. Compare it with a separate
quadratic-logit curve for each size, which contains the common scaling model.
Use 1,000 parametric binomial bootstrap replicates under the fitted common model,
seed 2026101601, to calibrate the likelihood-ratio statistic and saturated-model
deviance. Refit every replicate; count all numerical failures explicitly. Report
conditional percentile intervals for p_c and nu and Monte Carlo uncertainty of
tail probabilities. Record how often the old raw-sign crossing rule rejects
these common-scaling bootstrap datasets. This calibrates the rule under the
fitted model; it is not a new independent physical experiment.

## Z: finite-size sensitivity matrix

Fit the binomial likelihood in logit space, retaining the original parameter
bounds. Use minimum sizes 3, 5 and 7; the last has only three sizes and is
diagnostic. Use the full original rate window [0.025, 0.040] and the previously
defined narrower window [0.0275, 0.0375]. Evaluate all six models in every subset:

1. Quadratic leading scaling curve.
2. Cubic leading scaling curve.
3. Quadratic plus c (L_ref/L)^0.5.
4. Quadratic plus c (L_ref/L)^1.
5. Quadratic plus c (L_ref/L)^2.
6. Quadratic plus (c0+c1*x) (L_ref/L), where x is the normalized scaling variable.

These are a fixed exploratory sensitivity grid, not a model search stopped on
success. Keep all 36 results, failures and boundary hits. Report deviance,
degrees of freedom, asymptotic goodness, and AICc within identical data subsets
only. AICc is a descriptive heuristic using independent binomial cells as the
sample size; it does not certify the preferred model or a threshold.

For full-window L_min=3 and L_min=5, profile p_c for the leading quadratic and
constant omega=1 correction models. Reoptimize nuisance parameters at each p_c;
report every connected 95% likelihood-ratio region (cutoff chi-square_1 0.95),
boundary truncation and numerical failures. These are approximate conditional
intervals, not intervals that include model uncertainty. Do not pool model
ranges into a confidence interval or estimate omega freely from five sizes.

## Validation and delivery

Check analytic gradients, exact interval boundaries, synthetic known cases,
shared-model nesting, original-fit reproduction, and profile minima. Use fresh
independent review of statistics and saved numerical outputs. Do not rerun the
unchanged circuit suite. Save numerical tables, a concise interpretation, and
PDF figures without grids; diagnostic PDFs are appropriate for this analysis.
Preserve the previously requested data-only PDFs. Export the completed follow-up
to the normal project's results folder and verify copied bytes.
