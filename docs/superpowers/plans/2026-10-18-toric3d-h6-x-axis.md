# Physical x-axis H+6 Comparison Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans. Track completion below.

**Goal:** Run and report a verified physical-x growth comparison against the completed diagonal H+6 study.

**Architecture:** Create an isolated experiment bundle under results/toric3d-h6-x-axis-2026-10-18, using a copy of the frozen matched-study runtime and a separate AxisH6 module. Geometry and feedback are independently verified; a root-owned runner reuses the established checkpoint and paired-count design. Existing source, notes, and studies remain untouched.

**Tech Stack:** Julia, existing TopoNoise primitives, PythonCall, ldpc BP+OSD, PyMatching, Python/Matplotlib PDF reporting.

**Spec:** docs/superpowers/specs/2026-10-18-toric3d-h6-x-axis-design.md

## Global Constraints

- Preserve old scientific evidence and source hashes; no modification of historical studies.
- H+6 CNOTs, thirteen noisy transfer outputs per cell; matched 1,000 shots and L=3,5,7,9,11.
- Common seed repair across arms, no postselection, truthful geometry and policy differences.
- All corrections depend only on declared measured syndromes.
- PDF-only figures, archived paper-no-grid theme.

## Review Focus

- L=2 parallel physical edges must remain distinct in recovery IDs.
- Seed logical faults must survive common syndrome repair.
- Corrections on the oldest allowed slice must not create unmeasured retired defects.
- Partial files and resumed runs must preserve stream and outcome counts exactly.
- Timing pilots and verification records must not enter production data.

## Tasks

- [x] Geometry: add archive/axis_geometry.jl and test_axis_geometry.jl. Define AxisCode(T,Ly,Lz), AxisEncoder, axis_encoder, nqubits, edge, plane_edges, slice_qubits, star_checks, plaquette_checks, logical_x, logical_z, verify_encoder. Establish invariant tests before implementation; verify every prefix in both prepared logical bases.
- [x] Feedback: add archive/axis_feedback.jl and test_axis_feedback.jl. Define build_context(L,p,channel), sample_record(rng,context), three_outcomes(context,record;record_history=false). Reuse existing decoder settings and Pauli propagation. Test the declared feedback support and invariant failures, all noiseless arms, explicit seed logicals, and independent fault propagation.
- [x] Runner: create archive/AxisH6.jl wrapper and adapt a copied matched-study runner and supervisor, preserving atomic checksummed checkpoints and paired eight-bin counts. Add exact resume and corrupt-checkpoint tests. Freeze the final implementation before production.
- [x] Verification/pilot: run bounded geometry/replay tests and excluded timing pilots; independently review the circuit and feedback. Resolve failures before launching production.
- [x] Run/report: freeze selection and launch the matched axis scan. Produce one themed PDF comparing physical directions and all three arms, retain uncertainty and geometry labels, and verify final data and PDF before claiming completion.

## Decisions and progress

2026-10-18: User clarified that direction refers to sequential growth, not noise sector. Existing rough-boundary x-axis H+5 code is unsuitable. Independently derived square-seed H+6 transfer preserves k=2 and transfer gate count. Work is confined to a new experiment bundle rather than changing the user's active source or manuscript.

2026-10-18 verification: 376 geometry + 15,740 feedback + 127 runner Julia assertions and 23 Python tests pass. Independent causal review found no issue. All 18 excluded pilot points completed; interpolation suggests about 128 summed worker minutes before startup/report overhead. Two-worker production expected approximately 1.5–3 elapsed hours. Two-page PDF layout rendered and inspected.

2026-10-18 production dispatched: source/protocol/reference/verification hashes frozen in data/selection.json. Supervisor PID 50645, two workers. Task 5 remains in progress until the full scan and final report complete.

Completion: all 205 points finished at 2026-10-18 15:22:57 Asia/Shanghai after 92.4 minutes. Independent full data/source/statistics audit passed; both final PDF pages rendered and visually inspected. Final audit is archived in final-verification.json. No x-axis threshold fit was attempted or validated.
