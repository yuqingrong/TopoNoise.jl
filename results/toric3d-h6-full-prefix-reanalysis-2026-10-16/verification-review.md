# Independent review of the exploratory follow-up

The final numerical outputs pass the bounded independent review. One numerical
reporting defect was found during review and corrected before final results:
a failed tight optimizer polish could replace a converged nuisance fit at
essentially identical likelihood, creating artificial gaps in Z profiles. The
final code preserves the converged result in that case and records the polish
failure separately. All four final profiles have zero failed nuisance fits,
no numerical gaps, and agreement with their free-fit minimum. No unresolved
scientific or numerical blocker was found within the reviewed scope.

This remains an exploratory analysis of already inspected data. The original
plan was frozen before these follow-up fits, but after the original findings.
The all-ten-pair X extension was chosen after inspecting the primary adjacent
bands and is correctly identified in `analysis-addendum.md`. Neither this
follow-up nor the model-dependent parameter intervals change the original
acceptance decision or supply an independent physical experiment.

`independent_spotcheck.py` imports neither analysis implementation. Its final
run completed with **8,433 assertions**, saved in `independent-spotcheck.json`.
It checked the frozen input/plan hashes, the original report and artifact
manifest hashes, and all 70 unchanged source points (70,000 fault records).
Both saved reports match the final implementation hashes. Independent checks
reconstructed exact binomial intervals through a separate SciPy interface,
affine interval intersections, predicted probabilities, likelihoods,
deviances, degrees of freedom, AICc and within-subset comparisons. For every
saved Z profile evaluation they reconstructed the likelihood from its stored
nuisance parameters, checked the fixed p_c, total-likelihood LR scale, and
refined cutoff endpoints. All 1,000 X binomial datasets were regenerated from
the fixed seed and matched the saved counts; the old crossing rule was
implemented separately and reproduced for every replicate. Three bounded
convex coefficient fits at saved (p_c, nu) values independently reproduced
bootstrap likelihoods. No circuit simulation, full bootstrap, or profile was
rerun by this review.

The X simultaneous coverage construction is valid for the 35 measured binomial
probabilities. Between-grid claims additionally require the stated linear
interpolation assumption. All derived pairwise differences inherit the same
joint coverage event, so examining all ten pairs needs no further multiplicity
penalty. Intersecting adjacent-pair compatible regions alone does not ensure a
shared ordinate across all five sizes; the all-pair intersection does, because
it is equivalent to the maximum lower band being below the minimum upper band.
The primary result spans the entire measured [0.035, 0.065] range and therefore
does not localize a crossing. The secondary all-pair region is
[0.0461886250, 0.0606662786]. Compatibility establishes neither existence of a
unique crossing in the true curves nor an asymptotic threshold.

The independent quadratic-logit curves contain the common quadratic scaling
model exactly. Their observed deviance is 10.5686906617 versus 17.9459651518
for the common model, giving LR = 7.3772744900. All 1,000 bootstrap refits are
usable in both models, with no boundary or nonmonotonicity flags. The independent
tail recount gives 968/1,000 for common-model deviance and 697/1,000 for LR;
the saved plus-one probabilities are 0.968032 and 0.697303. These are plug-in
parametric-bootstrap calibrations, not exact finite-sample tests. Their Monte
Carlo intervals quantify simulation error conditional on the fitted null.
The old rule rejects 918/1,000 datasets generated under that fitted model;
this is a conditional diagnostic, not a universal false-rejection rate.
The conditional percentile intervals reproduce as p_c in
[0.0525937984, 0.0550210950] and nu in [1.2610878251, 1.6632089484]. They exclude
uncertainty in finite-size modeling, subset selection and correction terms.

All 36 prescribed Z model/subset combinations and all four profiles are
present. The leading full-window L >= 3 model has deviance 53.7248 on 30 df
(asymptotic p = 0.00493). Constant corrections improve the fit, but their
candidate p_c varies from 0.0354155 to 0.0376081 across fixed exponents. For
L >= 5 the leading model already has deviance 29.7663 on 23 df (p = 0.15623);
adding the constant inverse-size correction changes the deviance only to
29.7510. AICc is compared only within identical subsets and retains the
explicit independent-cell sample-size heuristic. Improved fit or lower AICc
does not itself validate a threshold or a particular correction exponent.

The final approximate 95% conditional profile regions are:

| Full-window subset | Model | p_c region |
| --- | --- | --- |
| L >= 3 | Quadratic | [0.03312360, 0.03409923] |
| L >= 3 | Quadratic + constant inverse-size correction | [0.03476279, 0.03854570] |
| L >= 5 | Quadratic | [0.03367279, 0.03509221] |
| L >= 5 | Quadratic + constant inverse-size correction | [0.03198904, 0.03795083] |

These use the asymptotic one-parameter LR cutoff with nuisance reoptimization
and fixed parameter bounds. They omit model uncertainty. The 301-point scan
and refined endpoints identify every detected connected component; they cannot
guarantee the absence of arbitrarily narrow islands between grid points.
Raw correction coefficients also depend on the subset-specific reference size
and should not be compared directly across subsets.

The machine-readable audit records the exact input, report and review-script
hashes, managed Python executable/package versions, and thread settings.
Its scope is numerical and methodological review; final figure inspection
and export-byte verification are handled separately by the parent task.
