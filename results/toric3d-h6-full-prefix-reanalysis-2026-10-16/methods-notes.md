# Statistical interpretation

Each simulated grid point supplies one binomial count for full-prefix state
failure. The final-only arm is not treated as an additional independent sample.
X and Z, and distinct size/probability points, use distinct random streams.

The leading model is

    logit(P) = a0 + a1*x + a2*x^2
    x = (p-p_c)/(p_max-p_min) * (L/L_ref)^(1/nu).

Here L_ref is the geometric mean of the sizes in that fitted subset. Correction
models add a size-dependent term in logit space. The normalization changes
coefficient units but not the meaning of p_c or nu. Different subset fits use
their own normalization, so their raw polynomial coefficients should not be
directly compared.

For X, Clopper-Pearson intervals are constructed at level 1-0.05/35 for each
measured probability. Bonferroni's union bound gives at least 95% simultaneous
coverage of these 35 probabilities under their binomial model. Subtracting
interval endpoints bounds each difference on the same event. Coverage between
measured probabilities additionally assumes piecewise-linear interpolation;
no unmeasured curvature is bounded by this construction. Exact binomial
intervals follow the [SciPy Clopper-Pearson documentation](https://docs.scipy.org/doc/scipy-1.14.1/reference/generated/scipy.stats._result_classes.BinomTestResult.proportion_ci.html).

The X likelihood comparison tests a common quadratic scaling curve against
separate quadratic-logit curves for each size. The more flexible model contains
the common model. Each bootstrap replicate regenerates binomial counts at the
fitted common-model probabilities and refits the unknown parameters in both
models. The same replicates calibrate the saturated-model deviance. Tail
probabilities use (exceedances+1)/(replicates+1). Refitting each replicate follows
the general [parametric Monte Carlo procedure documented by SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.goodness_of_fit.html);
the implementation here handles the study's binomial regression directly.

Bootstrap parameter intervals condition on the fitted model. A non-small
goodness-of-fit or likelihood-ratio p-value means the data do not show the
specified discrepancy; it does not prove the model or an asymptotic threshold.
The rejection rate of the old crossing rule is conditional on the fitted null
model, not a universal false-rejection rate.

For Z, AICc is compared only among models fitted to the same cells. Its small
sample adjustment uses the number of independent binomial cells as a descriptive
heuristic. Profile intervals use the asymptotic one-parameter likelihood-ratio
cutoff, with nuisance parameters reoptimized. Finite-size model uncertainty,
exploratory model comparisons, bound hits, or weak parameter identification can
limit their interpretation. A range across model estimates is not a confidence
interval, and no model is promoted to a confirmed threshold because its fit or
AICc happens to improve.
