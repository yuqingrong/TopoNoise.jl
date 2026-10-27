# Matched-preparation TN pilot — 27 October 2026

The requested initial seed repair is included. The new TN curve is a completed **L=3 pilot**, with **20 fresh shots at each of four noise rates** and bond dimension **χ=240**. All 80 outputs are valid, and all 80 TN logical decisions agree with independent exact ML. The full L=3,5,7,9,11 scan has not been completed; no threshold is fitted.

![TN pilot and original reference](figures/tnml-matched-figure.png)

[PDF figure](figures/tnml-matched-figure.pdf) · [Counts and intervals](figures/summary.csv) · [Full machine-readable report](figures/report.json)

## Observed failure counts

Failure means either of the two logical bits is inferred incorrectly. The paired MWPM column uses the same 20 sampled histories as the TN column. Intervals are pointwise 95% Wilson intervals. Zero observed failures does not establish a zero underlying failure probability.

| Gate-output p (%) | TN failures | TN rate [95% interval] | Paired MWPM failures |
|---:|---:|---:|---:|
| 0.60 | 0/20 | 0.000 [0.000, 0.161] | 0/20 |
| 1.00 | 0/20 | 0.000 [0.000, 0.161] | 2/20 |
| 2.00 | 1/20 | 0.050 [0.009, 0.236] | 5/20 |
| 4.25 | 13/20 | 0.650 [0.433, 0.819] | 16/20 |

## What is matched

The physical protocol matches the uploaded figure: independent Z faults on seed-H outputs, complete initial unit-weight MWPM seed repair, then diagonal H+6CNOT preparation with independent Z gate-output faults and no later feedback. Measurements and corrections are ideal. The seed-repair implementation and PyMatching 2.4.0 match the archived reference.

The TN decoder retains the measured seed syndrome and final syndrome, together with the known seed recovery. The final-MWPM baseline uses the final syndrome only. Consequently, the comparison includes both a decoder change and an information change. The decoder never receives sampled fault locations or the true logical mask.

The left panel is archived MWPM data for five sizes, with 1,000 shots per point. The right panel contains only the 80 newly simulated L3 histories. Earlier TN results with an unrepaired seed were not reused. Diagnostic histories are excluded from these failure-rate estimates.

## Numerical checks and limits

- Every pilot shot has an independent, untruncated dense-frontier exact reference. Strict logical-label matches: 80/80; tie-aware exact-ML matches: 80/80. Maximum absolute error in a normalized logical probability: 0.000270836051.
- A separate 12-shot L3 diagnostic set at the same four rates also gave 12/12 exact-ML agreements. The repaired-seed model passed 16,775 assertions, including all 4,096 L2 seed patterns and an independent physical probability calculation. The maximum L2 posterior discrepancy was 3.12e-11.
- Checkpoint/scoring tests: 35 assertions passed. Plot/report integrity: six tests passed. The supervisor passed normal completion, time-limit termination, and stale-progress checks.
- All 80 requested pilot histories completed. There were no invalid, exception, missing, discarded, or replacement shots. TN decisions were never replaced by exact decisions.
- One L5 trial at p=1%, χ=240 did not finish within the 180-second contraction limit. Its saved history remains an unknown diagnostic outcome. Peak observed memory was 3.09 GiB. L7, L9 and L11 were not launched.
- One numerical worker and one BLAS thread were used. The four pilot workers took 37.3 minutes in total, including setup, compilation and exact references. Peak observed worker memory was 2.74 GiB. Median TN times at increasing p were 4.5 s, 10.2 s, 36.9 s, 39.7 s.
- Numerical agreement is established for these sampled L3 records. Finite bond dimension remains an approximation; convergence at larger sizes is unestablished. This small, single-size pilot cannot determine a threshold.

## Reproduce or inspect

- [Experiment instructions and likelihood factorization](../../examples/tnml_repaired_seed/README.md). The integrated contraction in `src/` was unchanged; the matched-seed wrapper and scan scripts are under `examples/tnml_repaired_seed/`.
- `pilot-selection.json` was frozen before production sampling. It records rates, shot counts, independent RNG seeds, settings and resource limits.
- `production/` contains all sampled histories, per-shot TN outputs and exact references, versions/source hashes, progress, summaries and resource-monitor logs.
- `diagnostics/` and `validation/` retain diagnostic results, tests and resource-limit evidence. `reference/` preserves the original counts and selection. `archive/` contains the final source/configuration snapshot and SHA256 manifest.

From the repository root, regenerate the figure with:

```sh
python examples/tnml_repaired_seed/plot.py --root results/tnml-matched-figure-2026-10-27
```

The saved pilot can be resumed with `run_pilot.py` as described in the experiment instructions. Completed outputs are skipped only after source/configuration/history identity checks. New settings or fresh histories require new output directories.
