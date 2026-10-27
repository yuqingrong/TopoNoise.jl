# Full-grid TN study with tighter SVD cutoff — 27 October 2026

**Launched on 27 October 2026 at 10:37:52 Asia/Shanghai (02:37:52 UTC), starting with zero-noise checks.** The launch snapshot records the first L3, p=0 check; positive-noise production had not started in that snapshot. [Live launcher status](background-status.json) identifies the current phase and completed points; [launch.json](launch.json) records launcher PID 98745 and frozen identities. This study preserves the requested **L=3,5,7,9,11**, **p=0–4.25%**, and **1,000 independent histories per nonzero-p point**, totaling **110,000 requested statistical histories**. The numerical settings are χ=240, SVD cutoff `1e-26`, and the unchanged negative-probability tolerance `1e-10`. A fresh seed offset of 100,000 separates this study's histories from the earlier studies. The full statistical experiment remains unfinished.

The physical model is unchanged: independent Z faults at seed-H outputs, complete initial MWPM seed repair, diagonal H+6CNOT transfer with independent gate-output faults, and no later preparation feedback. Measurements and corrections are ideal. TN conditions on the seed and final syndromes and the known seed recovery. Logical truth is used only for scoring. The tighter SVD cutoff changes numerical truncation, not the physical error probability p.

## Why this is a separate study

The [earlier full-grid attempt](../tnml-full-grid-2026-10-27/README.md) stopped at L3, p=0.25%, after 73 saved outputs: 68 valid and five invalid due to negative probabilities. Those inputs and outputs remain preserved under their original cutoff and source identities. They are not relabeled or pooled into the new 1,000-shot targets. The earlier 80-shot pilot also remains a separate dataset.

Initial precision diagnostics on five fixed histories passed at both `1e-22` and `1e-26`. The completed replay of all 73 saved histories at `1e-26`, including all five originally invalid cases, produced **73 valid outputs and 73 strict exact-ML matches**, with maximum absolute posterior error **2.4463×10⁻¹²**. The [independent validation audit](../tnml-full-grid-2026-10-27/validation/precision-validation73-audit.json) preserves the evidence. Twenty-six cases reached the bond-dimension cap; these results support a monitored L3 restart, not a general convergence guarantee or validation of larger sizes. Diagnostic histories do not count as new production samples. The selection and gate evidence were frozen before launch.

## Execution sequence

The running `run_l3_background.py` uses the frozen `scan-selection.json`; it does not create or change the selection. The sequence is:

1. One actual p=0 regression for each of L=3,5,7,9,11, sequentially.
2. Only after all five checks pass, L3 production at each of the 22 positive rates, with 1,000 fresh histories and independent exact references per point.

The positive probabilities, in percent, are `0.25, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00, 4.25`.

There is one numerical worker and one BLAS thread at a time. The launcher refreshes the figure after each zero check or production point, stops on driver or reporting failure, and respects the driver's `STOP` handling. A launcher lock and the driver's study lock prevent concurrent dispatch. No automatic retry of blocked points is enabled. All metadata, settings, seed, source, project, and manifest identities must agree with the frozen selection; complete points also undergo saved-file and independent exact-reference audits.

L5–11 positive-noise production remains outside this launch sequence until runtime and contraction-dimension checks pass. The previous L5, p=1%, χ=240 diagnostic still had no decision at a 600-second TN limit. That old-setting result is not a timing measurement for this tighter cutoff, and new large-size feasibility is not established. Completing the selected L3 sequence would not complete the requested five-size study.

## Zero-noise interpretation and reporting

At p=0 the ideal model has exactly zero logical failure probability for every size. The plotted zero anchor is analytical, with no statistical confidence interval. The old study passed all five zero regressions; this study schedules its own checks under the new settings and identities. Numerical regressions are not counted as 1,000 statistical successes.

The figure displays an empirical rate only after an entire requested point is complete and audited. Invalid, interrupted, and missing TN outcomes remain unknown. No negative-probability clipping, tolerance relaxation, exact substitution, or valid-only denominator is allowed. For F known failures and U unknown outcomes among N requested histories, reporting retains the bounds `[F/N,(F+U)/N]`; these are not confidence intervals. No threshold is fitted from incomplete curves.

## Status and provenance

These files are created as preparation and execution proceed; an absent file is not evidence that a run occurred:

- [scan-selection.json](scan-selection.json): frozen rates, settings, independent seeds, limits, and source hashes.
- [launch.json](launch.json): recorded launch time, PID, and launcher/selection hashes.
- [validation/gates.json](validation/gates.json): approved numerical evidence and eligible production sizes.
- [background-status.json](background-status.json) and [background.log](background.log): launcher status, current phase, and output.
- [scan-status.json](scan-status.json): driver results, audits, and blocks.
- `zero/` and `production/`: saved raw histories, TN checkpoints, and exact references.
- [figures/report.json](figures/report.json), [figure PNG](figures/tnml-full-grid.png), and [figure PDF](figures/tnml-full-grid.pdf): current counts, provenance, outcome bounds, and requested-grid coverage.

The read-only command `python run_l3_background.py --check` inspects the frozen selection without starting another launcher. Precision validation is complete and the sequence has launched. Positive-noise L3 production follows only after all five new-setting zero checks pass; consult the live status files for whether that transition has occurred.
