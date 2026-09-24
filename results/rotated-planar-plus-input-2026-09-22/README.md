# As starts in all-plus: fully matched native noisy circuits

Native As starts directly in the ideal physical product state
`|+>^⊗(d^2)` and prepares `|+_L>`. Bp starts in ideal `|0>^⊗(d^2)` and
prepares `|0_L>`. Each As source-check block applies H to its fresh
representative, followed by the incoming CNOT tree. The entire As H/CNOT
schedule is the clockwise quarter-turn of Bp, with each CNOT's control and
target exchanged and time order preserved.

Both ideal product inputs are outside the noisy circuit. Every H gate in
the stored circuit remains noisy, as does each CNOT operand. This pairs
the complete noise model under rotation and exchanging X/Z. It replaces
the earlier realization of native As from all-zero using noisy H gates on
free inputs; that realization's figures remain in their original archives.

![Logical errors with fully matched native inputs and gates](critical-window-fit/plus_input_logical_errors-local-fit.png)

[Four-panel PDF](critical-window-fit/plus_input_logical_errors-local-fit.pdf) ·
[Zoomed fits and data collapse](critical-window-fit/near_critical_scaling.pdf) ·
[Raw counts](plus_input_logical_errors-raw.csv)

## Sampling and provenance

- Distances 9, 11, 13, 15; `p = 0:0.002:0.1`; boundary orientation `x_ns`.
- 50,000 shots per point, batches of 10,000, master seed 1235.
- Independent selected Pauli faults after each stored H and each CNOT operand.
- Ideal initial product states and final syndrome measurements; no idle or
  measurement noise. Independent uniform-weight CSS PyMatching decoders.
- Both As panels are independently resampled: 408 points, 20.4 million shots.
- Both Bp panels reuse 408 points, 20.4 million shots from the native-logical-error
  archive after exact gate, input, geometry, sampler, decoder, seed, and version
  checks. The executable noise definitions used by the batch sampler match the
  archive; changes to its leading documentation and the Yao oracle are separate.

The raw CSV reports residual logical X under X-only noise and residual logical
Z under Z-only noise. These are logical Pauli-frame diagnostics, including
logical operators that leave the prepared logical state unchanged.

Independent sampling means the As/Z and Bp/X curves fluctuate separately even
though their distributions are equal. The As data are not copied from Bp.
The same applies to As/X versus Bp/Z. The 20.4 million reused Bp shots are
not new samples relative to the earlier Bp figure.

`metadata.toml` records the source hashes, native product inputs, gate counts,
noise model, archive provenance, and completion status. `source_snapshot/`
retains the source used in the run. Distance checkpoints support reproducible
resumption with matching sources. `tests.log` records 6,464 passing assertions
in 34 testsets. The expanded duality checks cover 7,664 single H/CNOT fault
pairs and 640 paired many-fault histories, including both gate and plaquette
clocks, with identical decoded results under X/Z exchange.

## Local fitting

The fitting scripts in `critical-window-fit/` jointly fit all four distances
to `P_L = A + B*x + C*x^2`, where `x=(p-p_c)*d^(1/nu)`. Both As/Z and Bp/X
use the same primary window, `0.020 ≤ p ≤ 0.032`, with 28 data points and
five fit parameters per panel. This is the near-transition expansion
described in the [PanQEC threshold-fitting documentation](https://panqec.readthedocs.io/en/latest/tutorials/Computing%20threshold.html).

Weights are inverse observed binomial variances. Each primary fit uses
1,000 independent count-bootstrap replicates, refitting the model and
recomputing the weights each time. Fit errors are statistical only,
conditional on the window and model; there is no finite-size correction term.
Four window widths, quadratic/cubic models, and omission of the smallest or
largest distance are retained as sensitivity checks.

| Panel | p_c | nu | chi-square / dof |
|---|---|---|---|
| As/Z | 0.02672 ± 0.00022 | 1.56 ± 0.09 | 0.98 |
| Bp/X | 0.02620 ± 0.00022 | 1.54 ± 0.08 | 0.77 |

All 1,000 bootstrap fits per panel succeed. The independent threshold estimates
differ by 0.000518, approximately 1.64 combined statistical standard errors.
That difference is compatible with independent sampling; the paired-fault
tests establish the exact circuit symmetry. These finite-distance fits remain
conditional on the chosen scaling model and window.

For As/Z, the four quadratic windows give `p_c=0.02655–0.02673`; the widest
window has marginal fit quality. Cubic fits give `p_c=0.02655–0.02672`.
Removing d=9 gives `p_c=0.02666`, while removing d=15 gives `p_c=0.02660`.
Full results, statistical intervals, fit coefficients, and all sensitivity
checks are in `critical-window-fit/local_fits.json` and its companion CSVs.

`validate_results.py` verifies the input/source hashes, exact Bp archive reuse,
fit coefficients and uncertainties, bootstrap samples, test log, and all six
PNG/SVG/PDF exports. `validation.json` also compares the independently sampled
dual curves using a two-sample chi-square statistic, omitting points with fewer
than five expected failures or successes per experiment. Both PNGs were
visually inspected after rendering.

## Reproduce

From the repository root, using the installed project environment:

```sh
MPLCONFIGDIR=/tmp/toponoise-matplotlib \
JULIA_PKG_PRECOMPILE_AUTO=0 \
JULIA_CONDAPKG_BACKEND=Null \
JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" \
julia --compiled-modules=existing --project=. \
  -e 'include("results/rotated-planar-plus-input-2026-09-22/run.jl")'

MPLCONFIGDIR=/tmp/toponoise-matplotlib \
.CondaPkg/.pixi/envs/default/bin/python \
  results/rotated-planar-plus-input-2026-09-22/critical-window-fit/fit_local_scaling.py

MPLCONFIGDIR=/tmp/toponoise-matplotlib \
.CondaPkg/.pixi/envs/default/bin/python \
  results/rotated-planar-plus-input-2026-09-22/critical-window-fit/plot_local_scaling.py
```
