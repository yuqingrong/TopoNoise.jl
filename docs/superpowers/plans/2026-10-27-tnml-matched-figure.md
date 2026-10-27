# TN ML figure with matched seed preparation

**Goal:** Produce an honest logical-failure-versus-gate-noise figure using the integrated TN decoder for the physical preparation shown in the user's uploaded panel.

**Authorization and design:** The user explicitly chose to match the uploaded figure's preparation. It has noisy seed-H outputs, a complete initial seed MWPM syndrome repair, then the diagonal H+6CNOT circuit with no later feedback. The previous TN300 data omit this seed repair and cannot replace this panel. The TN decoder will retain the measured seed syndrome as well as the final syndrome; the reference final-MWPM decoder discards the seed record for its final decision. This information difference must be stated.

**Architecture:** Use the existing, unchanged boundary-MPS engine for two independent likelihoods: raw seed logical class conditional on the seed syndrome, and transfer logical class conditional on the final syndrome. The corrected seed has zero final syndrome. Shift the seed logical class by the actual seed recovery and XOR-convolve the four class weights. Split raw fault columns before rebuilding each reduced topology. Preserve signed outputs and component validity flags.

**Scope:** Experiment scripts under `examples/tnml_repaired_seed/`, provenance and figures under `results/tnml-matched-figure-2026-10-27/`. No changes to the integrated production contraction are planned. Existing dirty files and prior study results remain intact. No commit or publication is requested.

## Work and checks

- [x] Identify original figure, frozen settings, noise and preparation protocol.
- [x] Obtain the user's choice of matching seed preparation.
- [x] Independently review the factorization and logical-label mapping.
- [x] Implement a repaired-seed wrapper with test-first, independent L2 probability checks (16,775 assertions; full 4,096-seed enumeration and direct probability oracle).
- [x] Build a checkpointed experiment runner, saving raw histories before contraction, settings, source hashes, durations and all numerical failures (35 pure checks and a real L2 exact/TN/resume smoke).
- [x] Recompute exact L3 references for the new seed and transfer topologies; compare fixed diagnostic histories with TN results (12/12 decisions agree; maximum posterior absolute error 5.12e-5).
- [x] Probe larger sizes sequentially under explicit resource limits before choosing feasible production points (L5 stopped at the 180-second contraction cap; no decision completed; L7–11 not launched).
- [x] Freeze a pilot grid and sample count based on measured feasibility; retain incomplete and invalid histories in denominators (L3, p=0.006/0.01/0.02/0.0425, 20 fresh shots each; all 80 receive independent exact references).
- [x] Render PNG/PDF plots with 95% Wilson intervals for complete binary-outcome points, explicit finite-bond settings, shot counts and unavailable sizes. No threshold fit. All 80 pilot outcomes completed; TN failures 0/0/1/13, paired MWPM 0/2/5/16, each out of 20.
- [x] Independently review scoring, provenance and plot interpretation; verify final artifacts. Raw records, exact references, CSV/report and figure labels reconciled; 80/80 strict exact-ML decisions agree. Largest posterior absolute error 2.71e-4. Simulation source hashes verified and final sources archived.

## Original panel

Source: `results/toric3d-h6-matched-2026-10-17/data/matched-results.csv`, channel `z_only`, arm `final_only`; selection in `data/selection.json`, protocol in `archive/run_study.jl:342`.

Sizes `L=T=Ly=Lz=3,5,7,9,11`; 1,000 samples per point. Noise percentages: `0.45,0.50,0.55,0.60,0.65,0.70,0.75,1.25,1.50,1.75,2.00,2.25,2.50,2.75,3.00,3.25,3.50,3.75,4.00,4.25`. Colors: gray, blue, orange, green, purple in increasing size order.

## Numerical safeguards and resource policy

Use one numerical worker and one BLAS thread. Start from the validated engine settings `chi=240`, `cutoff=1e-18`, negative tolerance `1e-10`, but establish accuracy anew for this changed likelihood. Old cached posterior values are not valid here. Reuse the independent dense-frontier algorithm only with fresh topology-specific inputs.

An initial larger-size decode gets a 180-second contraction limit and a conservative memory cap where process inspection is available. Save the actual sample before attempting a decode. Resource-limited or invalid outcomes are unknown outcomes, never successes or silently discarded trials. For N attempted histories, F known failures and U unknowns, retain the identification interval `[F/N,(F+U)/N]`. Finite-chi instability must be disclosed; local SVD diagnostics do not bound posterior errors.

The full 100-point, 1,000-shot reference grid is not presumed feasible on this computer. The first deliverable is a measured, explicitly labeled pilot plus resumable scripts and feasibility results. Any reduction in size/grid/count must be visible in the plot and report.
