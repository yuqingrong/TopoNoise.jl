# Matched H+6-CNOT comparison

This study compares **no feedback after preparation**, **one-step feedback**, and **two-slice feedback** on the same sampled faults. Historical studies remain unchanged.

The approved sampling settings are:

| Sector | Sizes, with T = Ly = Lz = L | Shots per (L, p), per method |
|---|---|---:|
| X | 3, 5, 7, 9, 11 | 1,000 |
| Z | 3, 5, 7, 9, 11 | 1,000 |

Every method uses the entire common probability grid within its sector: 21 X probabilities and 20 Z probabilities. There are **205 points, 205,000 independent sampled records, and 615,000 protocol outcomes**. Each record is replayed by all three methods; the three outcomes are paired, not independent samples. Pilot and verification records are excluded.

The initial 10,000-shot proposal was replaced before any production sampling. This study uses the revised 1,000-shot settings above.

## Results and progress

- `data/progress.json`: live status, completed shots, points, and running workers, available after launch.
- `data/matched-results.csv`: one combined table for all three methods and both sectors.
- `data/report.json`: provenance, paired differences, and guarded fit results.
- `figures/h6-matched-comparison.pdf`: one PDF containing the comparison and fit panels. Interim versions are clearly marked incomplete and do not publish fits.
- `data/selection.json`: sampling grid and fit windows, frozen before production.

Only PDF figures are exported. Reproduction code, frozen simulation sources, excluded pilots, logs, and resumable checkpoints are contained in `archive/`. At most the current and previous checkpoint are retained for each point.

## What is held fixed

Within each sector all three methods have the same lattice sizes, transfer count T=L, physical fault probability, fault locations, number of shots, original sampled record, final decoder, failure definition, plot theme, size colors, and statistical fit policy. All samples count toward the result; there is no postselection.

X starts from the same ideal encoded logical |00> seed. Z starts from the same noisy seed-H logical |++> preparation, then receives the same complete seed syndrome repair at r=0 in every arm, as approved by the user. This repair preserves any undetectable seed logical error. Thus the Z baseline means **no feedback after common preparation** and is different from the historical uncorrected-seed baseline.

The X final decoder is full-code BP+OSD with the existing settings and declared p prior. The Z final decoder is unit-data-weight MWPM. The same final decoder object is used by all three methods at each point.

Feedback rules retain their existing meaning: one-step uses its established timing and limited correction; two-slice applies complete local physical correction after every transfer, including the last. In Z, one-step's local corrections end at T−1, whereas two-slice includes T. These terminal rules are part of the compared methods.

Two-slice X logical failure is fixed after the first corrected window because later operations do not touch the seed logical scoring cycles. Its size dependence describes preparation failure under this circuit and readout; it must not be interpreted as a noisy-memory threshold in the number of transfers.

## Fit policy

Seven probabilities per method and sector define the preselected local fit window. Every method is nevertheless simulated on the full union of these windows, permitting pointwise comparisons. Windows use historical transition ranges; they are not tuned to the new production outcomes.

Fits use the archived binomial-logit finite-size scaling analysis, with the same requirements for all six panels: minimum grid resolution, numerical convergence and parameter bounds, monotonicity, goodness of fit, crossing evidence, 200 bootstrap replicates, and sensitivity to dropping the smallest size, narrowing the window, cubic scaling, and an inverse-size correction. Five sizes leave four after dropping the smallest size. This makes the check possible; it does not guarantee a reliable threshold.

Rejected numerical candidates are retained and labeled diagnostic. A threshold value is reported only when all checks pass. Bootstrap intervals are conditional on the fitted model. Paired failure-rate differences use the eight joint outcome counts; independent threshold intervals are not used to claim a significant threshold improvement.

## Reproduction

The frozen runtime includes the actual modified and previously untracked source files used for these protocols. Hashes, runtime versions, exact per-point seeds, and configuration are checked on resume. An interrupted point samples and discards its already committed original records to restore the exact random stream.

Run or resume the complete study with the existing scientific Python environment:

```sh
/Users/rongyuqing/jcode/TopoNoise.jl/.CondaPkg/.pixi/envs/default/bin/python archive/supervise.py --study . --workers 2
```

The supervisor refuses duplicate launches, uses one thread per worker, saves progress, and generates the combined report and PDF automatically. It prevents idle system sleep while running; closing the laptop or shutting it down still interrupts computation. On interruption, rerun the same command to resume committed work.
