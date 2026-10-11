# Small-p test of the original Cartesian 3D toric-code sweep

This scan tests whether the original rough-boundary, Cartesian x-sweep could
have a very small positive logical-Z diagnostic threshold that was missed by
the first probability grid. It retains the original geometry, vertex order,
H/CNOT schedule, independent post-gate-output Z noise, ideal measurements and
corrections, unit-weight MWPM, and paired final-only/layer-feedback fault
records. It does not use the later diagonal-slice circuit.

The answer is no. Layer feedback substantially reduces finite-size logical-Z
failure, but its fixed-failure pseudothreshold moves toward zero approximately
as 1/L^2. This agrees with the decoder-independent lower bound

\[
P_{L,Z}\geq \frac{1-(1-2p)^{L^2}}{2}
\]

for cubic size L. At every fixed p>0, the right-hand side approaches 1/2 as
L grows. A smaller probability grid can delay the onset of failures at finite
L; it cannot produce a nonzero asymptotic threshold.

## Scan

The production scan used

- L=3,5,7,9;
- p=10^{-6},3\times10^{-6},10^{-5},3\times10^{-5},10^{-4},3\times10^{-4};
- 10,000 paired shots per (L,p), or 240,000 shots total;
- seed 24680;
- Z-only gate-output faults; and
- the same sampled fault record for final-only and layer-feedback recovery.

Selected observed logical-Z diagnostic rates are:

| L | p=10^-5 final / feedback | p=3x10^-5 final / feedback | p=10^-4 final / feedback |
|---:|---:|---:|---:|
| 3 | 0.06% / 0.06% | 0.25% / 0.25% | 0.73% / 0.72% |
| 5 | 0.69% / 0.21% | 1.63% / 0.55% | 5.81% / 2.10% |
| 7 | 1.03% / 0.29% | 3.36% / 1.03% | 9.49% / 3.38% |
| 9 | 3.41% / 0.73% | 8.67% / 1.78% | 22.59% / 6.01% |

At fixed p, failure increases with L under both protocols. Feedback lowers
the prefactor because intermediate recovery removes many correctable partial
strings, but it cannot resolve the indistinguishable first-layer fault pairs.

## Moving 1% pseudothreshold

For each size, define p_1% as the physical rate at which the sampled curve
reaches 1% logical-Z failure. Piecewise-linear interpolation between adjacent
simulated probabilities gives:

| L | final-only p_1% | layer-feedback p_1% |
|---:|---:|---:|
| 3 | 1.3803e-4 | 1.3862e-4 |
| 5 | 1.6596e-5 | 5.0323e-5 |
| 7 | 9.7042e-6 | 2.9189e-5 |
| 9 | 2.9726e-6 | 1.5143e-5 |

A log-log fit to p_1%(L)=A L^{-beta} gives

\[
\beta_{\rm final}=3.34\quad[3.13,3.52],\qquad
\beta_{\rm feedback}=1.97\quad[1.77,2.17],
\]

where brackets are conditional 95% intervals from 5,000 independent
binomial resamples of the measured point rates. All bootstrap replicas
bracketed the target. The final-only exponent is an empirical four-size fit,
not a universal exponent. For feedback, p_1% L^2 is nearly constant across
the scan, directly exposing the approximately 1/L^2 drift. Repeating the
interpolation at 0.5% and 2% failure gives the same qualitative conclusion.

The fixed-failure scale is a pseudothreshold, not a threshold estimate. A
positive threshold would require p_1%(L) to approach a nonzero constant.
Instead it decreases throughout the measured sizes, while the analytic bound
proves that it must ultimately vanish.

## Reproduction and artifacts

The scan API now accepts `channels=[:z_only]`, avoiding the unrelated X-noise
run while preserving the default two-channel interface:

```julia
using Random, TopoNoise
scan = compare_toric3d_feedback(MersenneTwister(24680);
    sizes=[3,5,7,9],
    error_rates=[1e-6,3e-6,1e-5,3e-5,1e-4,3e-4],
    shots=10_000, seed=24680, channels=[:z_only],
    output_dir="results/toric3d-original-small-p-2026-10-11",
    progress_io=stdout)
```

Run the statistical analysis with

```sh
.CondaPkg/.pixi/envs/default/bin/python \
  examples/analyze_toric3d_original_small_p.py \
  results/toric3d-original-small-p-2026-10-11/toric3d-original-small-p-raw.csv \
  --output-dir results/toric3d-original-small-p-2026-10-11/analysis \
  --target .01 --bootstrap 5000 --seed 24680
```

The [raw CSV](../../results/toric3d-original-small-p-2026-10-11/toric3d-original-small-p-raw.csv),
[metadata](../../results/toric3d-original-small-p-2026-10-11/toric3d-original-small-p-metadata.toml),
[analysis JSON](../../results/toric3d-original-small-p-2026-10-11/analysis/original-small-p-analysis.json),
and figures in
[PDF](../../results/toric3d-original-small-p-2026-10-11/analysis/original-small-p.pdf),
[SVG](../../results/toric3d-original-small-p-2026-10-11/analysis/original-small-p.svg), and
[PNG](../../results/toric3d-original-small-p-2026-10-11/analysis/original-small-p.png)
record the complete result. Wilson intervals are saved for every simulated
point. Logical Z remains a Pauli-frame diagnostic for the original prepared
|0_L> state; it is not the |0_L> preparation infidelity.
