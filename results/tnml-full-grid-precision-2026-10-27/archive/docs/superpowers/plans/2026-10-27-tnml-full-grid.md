# Full-range matched-preparation TN scan

**Goal:** Extend the user-authorized matched preparation to L=3,5,7,9,11, p starting at zero, and 1,000 independent histories per nonzero probability point.

**Design:** Preserve the existing scientific model, runner, and contraction backend. Start a new study rather than changing the frozen pilot. Preserve all reference probabilities; add 0%, 0.25%, and 1% to cover the low-noise gap and existing TN validation point. Thus use 22 nonzero rates from 0.25% through 4.25%. The zero-noise failure probability is analytically zero; run one deterministic implementation check per size rather than claiming Monte Carlo evidence there. Use exact references for every L3 production history. Higher sizes require separate contraction-dimension diagnostics before production, and resource limits must leave outcomes unknown rather than discard histories.

**User decisions:** Match the uploaded preparation (prior turn); all five sizes and p starting at zero; 1,000 shots per point (this turn).

**Implementation decisions:** One numerical worker and one BLAS thread. Preserve the old pilot and prior source hashes. Run finite resource probes before a days-long production launch. Make all production runs restartable with saved input histories, atomic per-shot results, and identity checks. Do not claim that a valid finite-bond result proves convergence. No exact fallback and no threshold fit.

## Tasks

- [x] Freeze full-grid configuration and add a tested resumable driver without changing the existing numerical runner.
- [x] Verify p=0 for all sizes using the actual TN backend; distinguish analytical anchors from simulation counts.
- [x] Extend the previous L5 time probe and record why larger-size production can or cannot proceed.
- [x] Check the scan/report selection and launch the validated portion of the requested 1,000-shot scan with durable logs and bounds.
- [x] Produce a current-state figure starting at zero and document exact completed counts, pending sizes, resource settings, and how to resume.
- [x] Replay all 73 saved cases at the tighter cutoff, compare with exact ML, and freeze a separate study for a monitored restart; preserve the original invalid outcomes and run identities.
- [ ] Finish the new study's five p=0 regressions before dispatching its L3 positive-noise production points.
- [ ] Complete the requested statistical scan where runtime and convergence validation permit; the target remains 110,000 histories and no new 1,000-shot point is complete.

## Execution status

The requested five-size statistical experiment is not complete. The original L3 attempt at p=0.25%, χ=240, and cutoff `1e-18` remains stopped with all data preserved: 73 saved outputs, of which 68 are valid and five invalid because posterior entries fell below `-1e-10`. All 1,000 raw histories remain saved. There are zero known TN failures and 932 unknown outcomes; its outcome bounds `[0, 0.932]` are not a confidence interval. The original `STOP` file and blocked status remain in place.

The completed precision replay at cutoff `1e-26` produced 73 valid outputs and 73 strict exact-ML matches, including all five originally invalid cases, with maximum posterior error 2.4463×10⁻¹². Twenty-six cases reached the bond cap; this is validation of those cases, not a general convergence guarantee. A separate study, `results/tnml-full-grid-precision-2026-10-27`, froze χ=240, cutoff `1e-26`, unchanged tolerance, and seed offset 100,000, then launched at 10:37:52 Asia/Shanghai on 27 October 2026 (02:37:52 UTC). Its initial phase is five sequential new-setting zero checks, followed by the 22 L3 production points only if every zero check passes. The launch snapshot records L3 p=0, not positive-noise production; current progress is in that study's `background-status.json` and `scan-status.json`.

All five original-study p=0 regressions passed; the new-study regressions are separate checks. The repeated identical L5 positive-noise diagnostic reached the corrected monotonic 600-second TN limit with no completed output, so L5–11 positive-noise production remains blocked on runtime and convergence validation. Fifty-seven Python tests passed after the configurable-cutoff update; independent review confirmed saved-file audit and frozen-provenance checks. A wall-clock-based supervisor timeout bug was reproduced, corrected, and tested before the first production launch.

The physical model, five requested sizes, 22 positive rates, 1,000-history targets, and requirement for independent exact references for every L3 production history remain unchanged. The earlier 80-shot pilot, stopped original attempt, diagnostic replays, and fresh precision-study production remain separate datasets. No probability clipping, tolerance relaxation, exact fallback, or replacement of invalid outcomes is allowed. The full 110,000-history target is unfinished.

## Validation

Unit tests cover full grid/rate/shot selection, nonzero-only Monte Carlo counts, resume behavior, refusing unvalidated sizes, worker errors/timeouts, and deterministic p=0 provenance. Numerical checks must use the existing contraction. A reviewer verifies that no pending/blocked point is displayed as observed success and that the plot does not invent missing curves.
