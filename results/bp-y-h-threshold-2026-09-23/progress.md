# Execution ledger — Bp Y/H threshold study, approved plan in task

The user approved exact stochastic Y/H faults, an ideal Bell reference, and
targeted 50,000-shot fits at d=9,11,13,15. Follow the current working-tree Bp
schedule and existing uniform-weight final-syndrome decoder.

Ruling: Use a new experiment archive in the current checkout. The required
physics includes existing uncommitted source changes; no existing package
files or previous run archives will be edited. A HEAD-only worktree would
silently change the experiment. All source files used will be snapshotted.

Tasks:
1. Export current circuits and independent Julia/Yao fixtures; implement and
   test exact-H and correlated-Y samplers and checkpointed counts.
2. Run the 2,000-shot pilot, resolve crossing ambiguity with 10,000 shots,
   freeze windows, and run fresh 50,000-shot production points.
3. Fit local finite-size scaling with joint bootstrap and sensitivity checks.
4. Independently review, validate, render figures, and write the report.

Pre-flight: exporter -> sampler: one-based Julia indices become zero-based
Python indices exactly once. Sampler -> fitter: explicit total p and joint
counts (00,10,01,11), never p_x+p_z. Fitter -> plotter: production-only primary
fits, pilot displayed separately, no fabricated threshold for unresolved panels.

Task 1: complete. Five sampler integration tests pass, including 73 Yao
statevector fixtures, fixed Y histories against the existing Julia propagation
and decoder, fast-H versus explicit replay, noiseless d=3,9,11,13,15, H²,
probability endpoints, repeatability, resume equivalence, and stale-source rejection.
Four analysis tests pass: known p_c/nu recovery, row-order invariance, crossing
selection, saturation rejection, flat-fit rejection, and correlated bootstrap.
The initial absent implementations failed the tests as expected. The exporter
had a first-line Julia docstring/import syntax error; replaced it with a comment,
then successfully exported all schedules and fixtures.

Ruling: Make pilot selection reproducible with declared numerical rules:
interior rates 0.005–0.45, three adjacent crossings within 0.010, and d=9/15
ordering reversal at 1.5 combined standard errors on either side. Unresolved
ambiguity after one 10,000-shot refinement stays unresolved. These rules can
miss a weak transition, so the report will describe the scanned range and will
not infer a zero threshold from absence of selection.

Task 2: running. All-stage scan started with four workers, 1,000-shot atomic
checkpoints, and source/version fingerprints. Full Julia regression suite is
running independently; package source remains untouched.

Task 2 correction: the initial H pilot had a clear extreme-distance ordering
reversal, but its noisy adjacent crossings spanned more than 0.010. The selector
incorrectly classified this as unresolved without using the approved refinement
budget. A saved observed-count regression fixture failed, then passed after a
fallback nominated 10,000-shot pilot refinement when d=9/15 bracket a candidate.
The fallback cannot nominate production directly. The initial selection is
preserved as selection-initial.json; no H production data existed. The Y window
and all Y production data remain fixed. The final H selection will still use
only its independent pilot data.

The first full Julia regression run hit a sandbox denial while Makie wrote its
normal scratch-usage cache. Retrying the same suite with the requested filesystem
permission, retaining both logs. This was not a physics assertion failure.

Refinement complete: 40 H pilot points now have 10,000 shots. Frozen production
centers are Y/X=0.028 and H/X=0.056. Both Z selections remain unresolved. Y's
already-frozen selection is unchanged. Fresh H production is running.

Independent whole-code review: no findings. The reviewer additionally checked
random H/Y histories at all five distances against tableau replay. Final output
completeness, numerical fits, figures/report, and Julia regression status were
pending and declined to judge. Ruling: verify those locally with the artifact
validator, completed-test exit codes, and rendered-image inspection; no further
physics implementation review is required absent a new finding.

Tasks 2–4: complete. The archive contains 256 points and 5,056,000 shots:
168 pilot points (40 H points refined to 10,000) and 88 fresh 50,000-shot
production points. Both fitted X panels completed 1,000/1,000 bootstrap fits.
Y/X gives p_c=0.02621390 ± 0.00025032, nu=1.61147 ± 0.08975;
H/X gives p_c=0.05367550 ± 0.00053258, nu=1.57507 ± 0.18387.
Both Z panels remain unresolved in the scanned range; no zero threshold is
inferred. All 20 fit-window, cubic, and distance-omission checks are archived.

Final verification: all 10 standalone scientific integration tests pass, and
all 6,464 assertions in 34 Julia testsets pass. The artifact validator checks
all 256 checkpoints, 512 plotted marginals, 2,000 bootstrap fits, and six
PNG/PDF/SVG exports. Both rendered figures were visually inspected; the main
caption spacing was corrected and the shading convention added. The final
report includes p_c and nu standard errors and 95% bootstrap intervals,
goodness of fit, sensitivity outcomes, provenance, and reproduction steps.

Requested presentation revision: matched the September 22 native logical-error
figure style, also shared by the matched-lattice archive. Restored d=9/11/13/15
blue/orange/green/pink colors, boxed axes, circle/line series, per-panel legends,
upper-left collapse insets, and the prior typography and layout. The H fit label
is placed above its low-lying curves to avoid covering data. Both rendered PNGs
were inspected. All six PNG/PDF/SVG exports were regenerated; no simulation,
threshold fit, bootstrap, window-selection, or sensitivity data were changed.
