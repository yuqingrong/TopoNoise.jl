# Bp circuit under Y-only and exact Hadamard faults

![Logical error curves and local threshold fits](logical_errors.png)

[Figure PDF](logical_errors.pdf) · [Near-critical fits and collapse](near_critical_scaling.pdf) · [Raw counts](raw_counts.csv) · [Fit summary](fit_summary.csv)

## Results

| Noise | Logical component | p_c | ν | χ²/dof | Fit window |
|---|---|---|---|---|---|
| Y-only | X | 0.026214 ± 0.000250 | 1.611 ± 0.090 | 0.75 | 0.022–0.034 |
| Y-only | Z | No supported crossing in the scanned range | — | — | — |
| Exact H=(X+Z)/√2 | X | 0.053675 ± 0.000533 | 1.575 ± 0.184 | 1.09 | 0.050–0.062 |
| Exact H=(X+Z)/√2 | Z | No supported crossing in the scanned range | — | — | — |

The ± values are one statistical standard error. Each logical-X and logical-Z rate is a component marginal: a joint logical-Y outcome contributes to both. The logical-Z curves do not support a positive finite-size crossing in the scanned range. This observation does not by itself establish a zero asymptotic threshold.

- Y-only/X: 95% statistical interval for p_c [0.025717, 0.026680] and for ν [1.455, 1.796]; 1000/1000 bootstrap fits succeeded; goodness-of-fit p = 0.798.
- Exact H=(X+Z)/√2/X: 95% statistical interval for p_c [0.052584, 0.054647] and for ν [1.299, 1.987]; 1000/1000 bootstrap fits succeeded; goodness-of-fit p = 0.351.

These are finite-size encoding-channel estimates conditional on this Bp schedule, the perfect final syndrome round, the independent CSS matching decoder, and the selected scaling window. Statistical intervals exclude window, polynomial, and distance-selection effects.

## Exact physical model and logical readout

Each active data-qubit operand receives an independent post-gate fault with total probability p. The two channels are E_Y(ρ)=(1−p)ρ+pYρY and E_H(ρ)=(1−p)ρ+pHρH, where H=(X+Z)/√2. Each Y event flips both Pauli-frame components together. H is replayed as an actual Hadamard in a full Stim stabilizer trajectory, retaining its quantum action through all subsequent gates and measurements.

The data use the current native Bp gate schedule, boundary orientation x_ns, at distances 9, 11, 13, and 15. Fault locations are the stored H gates and individual CNOT operands. Input preparation, the reference qubit, and final measurements are ideal; there is no idle or measurement noise.

To characterize both logical components with one fixed noisy schedule, the ideal input is (U_Bp†⊗I_R)|Φ_L+⟩. Final checks are followed by the commuting observables Z_L Z_R and X_L X_R. Their signs are XORed with the X- and Z-sector corrections predicted by the existing uniform-weight PyMatching decoder, with correlations disabled. Only final syndromes enter the decoder. The two failure bits f_X,f_Z give joint counts n00,n10,n01,n11, with P_X=(n10+n11)/N and P_Z=(n01+n11)/N. This Bell-reference experiment characterizes the logical channel; the reference=0 branch reproduces the native all-zero Bp input.

## Sampling and fitting

The archive contains 256 points and 5,056,000 shots: 168 pilot points (656,000 shots) and 88 fresh production points (4,400,000 shots). The pilot uses p=0:0.005:0.1 and 2,000 shots per point; 40 ambiguous-region points were extended to 10,000 shots. Production uses 50,000 independent shots per point and Δp=0.002. Both logical marginals are retained from each trial.

Crossing candidates exclude p=0 and pilot logical rates outside (0.005,0.45). The selected adjacent-distance crossings must form a compact cluster, with supported d=9 versus d=15 ordering reversal. A broad but supported extreme-distance reversal triggers further pilot sampling instead of an immediate no-crossing conclusion. The full numerical rules and candidate lists are in analysis.py and selection.json. The initial Y selection was retained unchanged while the H pilot was refined; no H production data existed before its window was frozen. Both selection records are archived.

Production centers are Y-only/X: 0.028, Exact H=(X+Z)/√2/X: 0.056. Each production window extends ±0.010; the primary fit uses ±0.006.

The local model is P_L=A+B*x+C*x², x=(p−p_c)d^(1/ν), with inverse observed-binomial-variance weights and four optimizer starts. ν is searched over (0.3,5). The [PanQEC threshold-fitting tutorial](https://panqec.readthedocs.io/en/latest/tutorials/Computing%20threshold.html) describes this near-transition expansion. All four distances are fitted jointly. Each of 1,000 multinomial bootstrap draws resamples joint outcome counts and supplies both logical marginals, preserving their correlation. Variances are recomputed each draw. Accepted fits require interior, identifiable parameters, physical increasing probabilities, at least 95% bootstrap success, and goodness-of-fit p≥0.01.

## Sensitivity

- Y-only/X, quadratic windows: p_c=0.026051–0.026264; χ²/dof=0.75–1.54.
- Y-only/X, cubic windows: p_c=0.026116–0.026222; χ²/dof=0.78–0.97.
- Y-only/X, distances 11, 13, 15: p_c=0.026144, ν=1.582, χ²/dof=0.72.
- Y-only/X, distances 9, 11, 13: p_c=0.026416, ν=1.650, χ²/dof=0.90.
- Exact H=(X+Z)/√2/X, quadratic windows: p_c=0.052934–0.053675; χ²/dof=1.03–1.29.
- Exact H=(X+Z)/√2/X, cubic windows: p_c=0.052576–0.053693; χ²/dof=1.01–1.28.
- Exact H=(X+Z)/√2/X, distances 11, 13, 15: p_c=0.053260, ν=1.572, χ²/dof=0.89.
- Exact H=(X+Z)/√2/X, distances 9, 11, 13: p_c=0.053769, ν=1.668, χ²/dof=1.36.

All attempted quadratic/cubic fits at halfwidths 0.004,0.006,0.008,0.010 and both distance omissions are retained in sensitivity.csv and sensitivity.json, including rejected or poorly fitting cases. These alternatives share samples and are not independent estimates.

## Verification and provenance

The standalone tests verify noiseless Bell checks at d=3,9,11,13,15; 73 independent Yao statevector histories at d=3; Y histories against the original Julia Pauli propagation and decoder; exact H batch replay; multiple faults and H²=I; probability endpoints; seeded repetition; checkpoint continuation; changed-source rejection; a synthetic known threshold/exponent; row-order invariance; saturation/flat-fit rejection; and correlated bootstrap behavior. The dispersed H pilot is retained as a regression fixture for refinement.

The complete Julia package regression suite passed on retry after granting access to its normal plotting cache; the original sandbox denial and successful retry logs are both retained. Exact test exit codes and log hashes are in test_status.json. validation.json verifies every checkpoint count, CSV marginal and error bar, source hash, production window, fit objective, bootstrap interval, plotted series, and all six figure exports.

Both final PNG figures were visually inspected for readable labels, complete panels, error bars, fit-window shading, scaling insets, and caption spacing. The PDF and SVG exports use the same figure objects. Styling follows the September 22 native logical-error figures: the same distance colors, boxed axes, per-panel legends, upper-left collapse insets, and typography. figure_style_revision.json records the reference and unchanged data/fit hashes for this presentation-only revision.

metadata.json records runtime versions, seeds, circuit/sampler fingerprints, and completion. source_snapshot retains the actual working-tree Julia sources and dependency manifests. Sampling uses deterministic seeds per 1,000-shot chunk and atomic checkpoints, independent of worker scheduling. Reproducibility is conditional on the recorded Stim version and machine architecture. No existing package API or earlier result archive was modified.

## Reproduce

Run from the repository root in the installed project environment. Completed points resume without resampling; a changed sampler or environment fingerprint is rejected. requirements-python.txt records the exact Python package versions used.

```sh
MPLCONFIGDIR=/tmp/toponoise-bp-yh-mpl OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .CondaPkg/.pixi/envs/default/bin/python \
  results/bp-y-h-threshold-2026-09-23/run_experiment.py --stage all --workers 4

MPLCONFIGDIR=/tmp/toponoise-bp-yh-mpl .CondaPkg/.pixi/envs/default/bin/python \
  results/bp-y-h-threshold-2026-09-23/analysis.py
MPLCONFIGDIR=/tmp/toponoise-bp-yh-mpl .CondaPkg/.pixi/envs/default/bin/python \
  results/bp-y-h-threshold-2026-09-23/plot_results.py
MPLCONFIGDIR=/tmp/toponoise-bp-yh-mpl .CondaPkg/.pixi/envs/default/bin/python \
  results/bp-y-h-threshold-2026-09-23/validate_results.py
```

The sampler reads archived circuits.toml directly. export_circuits.jl regenerates schedules and Yao fixtures from matching Julia sources and refuses to overwrite a different export.

For a fresh replay, make a new directory under results/ and copy the Python scripts, circuits.toml, source_hashes.toml, source_snapshot/, oracle_fixtures.toml, fixtures/, and requirements-python.txt into it. Leave checkpoints/, metadata.json, and selection*.json uncopied. Run the new directory's runner, analysis, and plot scripts with the same commands. The unchanged root seed, chunk size, architecture, and recorded package versions reproduce the simulated counts. Final validation additionally requires that directory's completed test logs and test_status.json, as in this archive.
