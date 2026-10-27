# Full-range matched-preparation TN experiment — 27 October 2026

Requested scope: **L=3,5,7,9,11**, **p=0–4.25%**, and **1,000 independent histories per nonzero-p point**. This is a long-running experiment, not a completed five-size result.

The frozen grid has 22 positive probabilities (percent):

`0.25, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00, 4.25`.

This preserves every original reference point and adds 0.25% and 1%. The requested statistical total is **110,000 histories**. The prior 80-shot pilot is a separate dataset; its observations are not counted toward these new targets.

## The zero-noise point

For p=0, every fault variable is zero. With ideal measurements and corrections, the logical failure probability is **exactly zero for every L**. The plotted zero point is this analytical prediction and has no binomial error bar. One actual TN contraction per size checks the implementation; these deterministic regressions are not 1,000 independent statistical successes.

Numerical zero checks and any retries are saved in `zero/`. The first L11 attempt exposed a timing issue: the supervisor compared the current wall clock to Julia's stage timestamp, while total elapsed time used a monotonic clock. It reported a 300-second stage timeout after only 133.9 monotonic seconds in total. That stopped run remains an unknown numerical outcome; it is not evidence that contraction required 300 seconds. The supervisor was corrected to use a monotonic stage timer with owned-worker PID checking. All eight timer lifecycle tests passed. Repeating the identical saved L11 history succeeded in 267.2 monotonic seconds including setup, and all five zero-noise checks now pass the strict source/history/posterior audit. The first interrupted attempt is preserved as `zero/L11-p0.0000-chi240/supervision-wall-clock-attempt.json`.

## Extended large-size probe

The same saved L5, p=1%, chi=240 history from the first pilot was retried with a ten-minute TN limit using the corrected monotonic timer. It still had no completed decision after 600.3 seconds in the TN stage (634.5 seconds including setup), with peak observed memory 2.91 GiB. The history remains unknown; no outcome was substituted. L5–11 positive-noise production remains blocked pending runtime and contraction-dimension validation. This is a limitation of the current run and settings, not a proof that those sizes are mathematically undecodable.

## Physical preparation and decoding

The preparation is unchanged from the approved matched experiment: independent Z faults on seed-H outputs, complete initial unit-weight MWPM seed repair, diagonal H+6CNOT transfer with independent Z faults on every gate-output qubit, and no later preparation feedback. Measurements and corrections are ideal.

TN conditions on the seed and final syndromes and the known seed recovery. Paired final-MWPM uses only the final syndrome. Raw sampled faults and logical truth are used for scoring, never as TN decoder inputs. Finite bond dimension χ=240 remains approximate; exact comparisons never replace TN decisions.

## Execution and validation

`scan-selection.json` freezes sizes, rates, 1,000-shot targets, fresh RNG seeds, settings, resource limits, and numerical-source hashes. `validation/gates.json` records which sizes have enough evidence to begin a monitored production scan. L3 is eligible based on the earlier 80 exact comparisons and receives an independent exact reference for **every** new history, including at new probability values. Larger positive-noise sizes require runtime and contraction-dimension checks before production.

One numerical worker and one BLAS thread run at a time. Production limits are 600 seconds per TN call, 24 hours per probability point, and 4 GiB resident memory. The scan checks source hashes before each point and saves every requested input history before decoding. Existing outputs are resumed only after strict configuration/source/history identity checks. Timed-out, invalid, or missing outcomes remain unknown.

Each point is audited after its worker finishes. An unsuccessful point, missing exact reference, or disagreement with exact ML beyond tied maxima blocks later points for that size. This is a per-point audit, not an immediate per-shot cancellation. It does not certify convergence for untested larger sizes. No threshold is fitted.

## Inspect and resume

- `background-status.json` and `background.log`: local launch progress and output.
- `scan-status.json`: driver state and per-point audit results (created when the driver runs).
- `production/*/summary.toml`: current saved/valid/unknown counts; updated after every history.
- `production/*/shots/` and `exact/`: independent TN results and exact references.
- `figures/tnml-full-grid.png` and `.pdf`: snapshot of completed curves and requested-grid coverage.
- `figures/report.json`: counts, unknown-outcome bounds, provenance, and blocking reasons.

The new full-grid scripts extend the experiment only; they do not alter the production numerical decoder in `src/`. The existing runner has no posterior cache, so even repeated syndromes are recomputed. Prior timings imply several days for the full L3 curve alone; larger-size runtimes cannot be inferred reliably from L3.

Use the same Julia/Python environment as [the experiment instructions](../../examples/tnml_repaired_seed/README.md). The driver command for the currently eligible size is:

```sh
python examples/tnml_repaired_seed/run_full_grid.py \
  --root results/tnml-full-grid-2026-10-27 --phase production --sizes 3 \
  --julia /Users/rongyuqing/.juliaup/bin/julia
```

The study-specific `run_l3_background.py` runs this driver once per L3 probability point and refreshes the figure after each completed point. Its runtime records distinguish the requested five-size grid from the selected L3 portion.

Create a file named `STOP` in this study directory to stop dispatch **after the current probability point finishes**. The driver will not start another point while it exists. Remove it to resume deliberately. Resource failures are not retried automatically; inspect the reason before using `--retry-blocked`. A process lock prevents two drivers from dispatching workers for this study at once.

Regenerate the snapshot without running contractions:

```sh
python examples/tnml_repaired_seed/plot_full_grid.py \
  --root results/tnml-full-grid-2026-10-27
```

This plot gives empirical rates only for complete, audited 1,000-shot points. It does not divide by just the successfully decoded subset. For N requested histories, F known failures and U unknown outcomes, it records the bounds `[F/N,(F+U)/N]`; these bounds are not confidence intervals.
