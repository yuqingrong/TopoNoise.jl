# H+5 Z-Boundary Ambiguity Controls Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a reproducible single-Z-fault ambiguity audit, a circuit-aware finite-size control decoder, parity-erasure diagnostics, and a hard pre-fit guard for the existing three-dimensional H+5 upper-boundary instrument.

**Architecture:** Refactor replay once so production and diagnostics share the exact pre-final state. Build immutable fault signatures from that state, group them by the physically available `(final star syndrome, observed upper parity)` key, and use the groups both for an ambiguity certificate and for a maximum-likelihood logical-frame adjustment layered on the existing MWPM recovery. Keep all historical scans unchanged; serialize new evidence in a separate result directory and pass its certificate to an optional Python Z-fit guard.

**Tech Stack:** Julia 1.12.6, Python 3.12.13, PyMatching 2.4.0, ldpc 2.4.1, NumPy 2.5.1, existing TopoNoise test and atomic-output helpers.

**Spec:** `docs/superpowers/specs/2026-10-14-toric3d-h5-z-boundary-ambiguity-design.md`

## Global Constraints

- Work in the existing `/Users/rongyuqing/jcode/TopoNoise.jl` checkout on branch `3dTC`; use native execution and local commits only.
- Preserve every existing H+5 selection, checkpoint, CSV, figure, report, and verification artifact byte-for-byte.
- Never modify `document/note.tex` or its compiled note artifacts. Its required SHA-256 is `3a40710c2ab6bd43e53f5f7703b7c8716b209dc1e060a5092f17e0edce073f3a`.
- Keep the canonical four-arm replay and unit-length MWPM outputs as the baseline. New decoding behavior is an explicitly named control.
- The control decoder may use only final star syndrome and observed upper parity. It may not use the clean branch, exact fault location, or exact observed upper bit pattern.
- Do not start a new large threshold scan. A failed single-fault guard stops fitting and bootstrap before numerical optimization.
- Use the recorded single-thread environment: `JULIA_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, `VECLIB_MAXIMUM_THREADS=1`, `JULIA_CONDAPKG_BACKEND=Null`, and `JULIA_PYTHONCALL_EXE=$PWD/.CondaPkg/.pixi/envs/default/bin/python`.
- Follow red-green-refactor for every production behavior; each task records the expected failing test before implementation.
- Do not push, merge, publish, submit remote jobs, or implement a 4D/boundary-circuit redesign.

## Review Focus

- A random even clean branch must not leak through the audit key; tests vary clean branches and require identical relative frames and observable parity.
- Equal-count logical classes must use the specified deterministic order `(00,01,10,11)`; tests pin the tie break.
- A multi-fault syndrome absent from the single-fault model must return the original MWPM recovery and mark a fallback.
- A parity-erasure point with zero accepted shots must serialize missing conditional values, not zero or NaN.
- Existing artifacts and `document/note.tex` must remain immutable; final verification checks their existing delivery audit and the exact note hash.

---

### Task 1: Expose the exact pre-final replay state

**Files:**
- Modify: `src/toric3d/vertex_projector_replay.jl:15-115`
- Modify: `src/TopoNoise.jl:163-166`
- Modify: `test/toric3d/vertex_projector_replay.jl`

**Interfaces:**
- Produces: `Toric3DVertexPreFinalState` with `protocol`, `channel`, `before_boundary`, `boundary`, `frame`, `history`, and `upper_parity_odd`.
- Produces: `replay_toric3d_vertex_projector_pre_final(e, r, f; protocol, channel, record_history=true)::Toric3DVertexPreFinalState`.
- Preserves: `replay_toric3d_vertex_projector` return type, logical outcomes, histories, recovery signatures, and checkpoint digests.

- [ ] **Step 1: Write failing replay tests**

Add a testset that calls the wished-for pre-final API and asserts:

```julia
pre = replay_toric3d_vertex_projector_pre_final(e, r, f;
    protocol=:window_complete, channel=:z_only)
out = replay_toric3d_vertex_projector(e, r, f;
    protocol=:window_complete, channel=:z_only)
@test pre.upper_parity_odd == isodd(count(pre.boundary.observed))
@test pre.frame.z == out.boundary.z
@test pre.before_boundary == out.before_boundary
@test [(h.checkpoint, h.channel, h.signature) for h in pre.history] ==
      [(h.checkpoint, h.channel, h.signature) for h in out.history]
@test out.upper_parity_odd == pre.upper_parity_odd
```

Cover clean even branches and an injected last-layer Z fault. Extend `h5_test_outcome` to validate the parity field for every existing replay test.

- [ ] **Step 2: Run the focused test and verify RED**

Run:

```bash
env JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 julia --project=. --startup-file=no -e 'using Test; include("test/toric3d/vertex_projector_test_helpers.jl"); include("test/toric3d/vertex_projector_replay.jl")'
```

Expected: FAIL because `Toric3DVertexPreFinalState`, `replay_toric3d_vertex_projector_pre_final`, and `out.upper_parity_odd` do not exist.

- [ ] **Step 3: Refactor replay without changing scientific behavior**

Implement the two interfaces in `vertex_projector_replay.jl`. Move propagation, intermediate feedback, and boundary completion from `_vertex_replay` into the pre-final function. Construct `frame` as a copied `PauliFrame` after `frame.z .= boundary.z`. Make `_vertex_replay` call the pre-final function and then the unchanged `_vertex_final_recovery`.

Add `upper_parity_odd::Bool` to `Toric3DVertexReplayOutcome`, populate it from the pre-final state, and export the new struct/function from `src/TopoNoise.jl`.

- [ ] **Step 4: Run replay and reference tests and verify GREEN**

Run the focused command above, followed by:

```bash
env JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 julia --project=. --startup-file=no -e 'using Test; include("test/toric3d/vertex_projector_test_helpers.jl"); include("test/toric3d/vertex_projector_reference.jl")'
```

Expected: PASS with unchanged logical/recovery assertions.

- [ ] **Step 5: Commit**

```bash
git add src/TopoNoise.jl src/toric3d/vertex_projector_replay.jl test/toric3d/vertex_projector_replay.jl
git commit -m "refactor: expose H+5 pre-final boundary state"
```

### Task 2: Build the exhaustive Z-fault ambiguity audit

**Files:**
- Create: `src/toric3d/vertex_projector_z_audit.jl`
- Create: `test/toric3d/vertex_projector_z_audit.jl`
- Modify: `src/TopoNoise.jl:163-170`
- Modify: `test/toric3d/runtests.jl`

**Interfaces:**
- Consumes: `replay_toric3d_vertex_projector_pre_final` from Task 1 and `decode_toric3d_vertex_z_complete`.
- Produces: immutable `Toric3DVertexZFaultSignature` fields `step::Int`, `slot::Int`, `qubit::Int`, `final_layer::Bool`, `syndrome_indices::Tuple{Vararg{Int}}`, `upper_parity_odd::Bool`, and `logical_z::NTuple{2,Bool}`.
- Produces: immutable `Toric3DVertexZObservableGroup` fields `syndrome_indices`, `upper_parity_odd`, `class_counts::NTuple{4,Int}`, `selected_class::NTuple{2,Bool}`, `baseline_class::NTuple{2,Bool}`, `unavoidable_count::Int`, and `baseline_failure_count::Int`.
- Produces: immutable `Toric3DVertexZParitySummary` fields `mechanism_count`, `group_count`, `mixed_group_count`, `unavoidable_count`, and `baseline_failure_count`, all `Int`.
- Produces: immutable `Toric3DVertexZAmbiguityAudit` with `dimensions`, `protocol`, `final_layer_only`, `signatures`, `groups`, total counts matching the parity-summary names, and `parity::NamedTuple{(:even,:odd),Tuple{Toric3DVertexZParitySummary,Toric3DVertexZParitySummary}}`.
- Produces: `toric3d_vertex_z_fault_signatures(e, f; protocol, final_layer_only=false)::Vector{Toric3DVertexZFaultSignature}`.
- Produces: `audit_toric3d_vertex_z_boundary(e, f; protocol, final_layer_only=false)::Toric3DVertexZAmbiguityAudit`.
- Uses the immutable lookup key `(syndrome_indices::Tuple{Vararg{Int}}, upper_parity_odd::Bool)` and class order `(false,false)`, `(false,true)`, `(true,false)`, `(true,true)`.

- [ ] **Step 1: Write failing enumeration and grouping tests**

Create standalone tests that first cover `L=3`:

```julia
audit = audit_toric3d_vertex_z_boundary(e, f;
    protocol=:none, final_layer_only=true)
@test length(audit.signatures) == 99
@test length(audit.groups) == 36
@test audit.mixed_group_count == 5
@test audit.unavoidable_count == 8
@test audit.baseline_failure_count == 34
@test audit.parity.even.mixed_group_count == 0
@test audit.parity.even.mechanism_count == 27
@test audit.parity.odd.mechanism_count == 72
```

Also assert total active-output counts `11*T*Lu*Lv`, exact output slots, deterministic class order, and argument rejection for an invalid protocol.

- [ ] **Step 2: Run the audit test and verify RED**

Run:

```bash
env JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 julia --project=. --startup-file=no test/toric3d/vertex_projector_z_audit.jl
```

Expected: FAIL because the audit types and functions are undefined.

- [ ] **Step 3: Implement fault signatures and audit grouping**

In `vertex_projector_z_audit.jl`:

- construct one all-zero immutable Z-only record with an even clean branch;
- enumerate `(step, slot, qubit)` directly from `e.operations` and inject exactly one `:Z` fault through `with_toric3d_vertex_projector_fault`;
- use the Task 1 pre-final replay, never a copied propagation algorithm;
- store final syndrome indices, observed upper parity, pre-recovery logical-Z tuple, and final-layer membership;
- group by the decoder-visible key, count classes in the specified order, choose the largest count with first-index tie break, and calculate unavoidable count;
- decode each group once with existing final MWPM to calculate its recovery logical class and baseline failures; and
- compute total and parity-stratified summaries for the selected enumeration; callers obtain a separately identified final-layer certificate with `final_layer_only=true` rather than hiding a second enumeration inside one audit.

Validate encoder/decoder dimensions and canonical protocol membership. Include the new file after replay and before scan in `src/TopoNoise.jl`, export the four public types and two public functions, and include the test in `test/toric3d/runtests.jl`.

- [ ] **Step 4: Add branch and protocol invariance tests**

Before changing production code further, add tests that:

- replay selected last-layer faults with at least three different even clean branches and require identical relative post-boundary frame, signature key, and logical class;
- require the exact final-layer signature vector to match across all four protocols; and
- check `L=5` and `L=7` final-layer totals, mixed groups, unavoidable counts, baseline failures, and even/odd mechanism totals against the approved table.

Run the test and expect RED on any missing summary or branch-invariance behavior, then make only the minimal audit changes needed.

- [ ] **Step 5: Run audit and existing replay tests and verify GREEN**

Run:

```bash
env JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 julia --project=. --startup-file=no -e 'using Test; include("test/toric3d/vertex_projector_z_audit.jl"); include("test/toric3d/vertex_projector_replay.jl")'
```

Expected: PASS for `L=3,5,7`, including 99/275/539 locations, 5/9/13 mixed groups, 8/14/20 unavoidable faults, and 34/131/276 baseline failures.

- [ ] **Step 6: Commit**

```bash
git add src/TopoNoise.jl src/toric3d/vertex_projector_z_audit.jl test/toric3d/runtests.jl test/toric3d/vertex_projector_z_audit.jl
git commit -m "feat: audit H+5 Z-boundary ambiguity"
```

### Task 3: Add the circuit-aware control decoder and parity estimate

**Files:**
- Create: `src/toric3d/vertex_projector_z_control.jl`
- Create: `test/toric3d/vertex_projector_z_control.jl`
- Modify: `src/TopoNoise.jl`
- Modify: `test/toric3d/runtests.jl`

**Interfaces:**
- Consumes: `Toric3DVertexZAmbiguityAudit` and `Toric3DVertexPreFinalState`.
- Produces: `Toric3DVertexCircuitZDecoder`, `Toric3DVertexCircuitZDecision`, `Toric3DVertexZControlOutcome`, and `Toric3DVertexParityErasureEstimate`.
- Produces: `build_toric3d_vertex_circuit_z_decoder(f, audit)::Toric3DVertexCircuitZDecoder`.
- Produces: `decode_toric3d_vertex_circuit_z(d, syndrome, upper_parity_odd)::Toric3DVertexCircuitZDecision`.
- Produces: `finalize_toric3d_vertex_z_control(pre, d)::Toric3DVertexZControlOutcome`.
- Produces: `toric3d_vertex_parity_erasure_estimate(upper_parity_odd, failures)::Toric3DVertexParityErasureEstimate`.
- `Toric3DVertexZControlOutcome` reports the common final syndrome plus both baseline and control residual logical-Z tuples/failure flags, so diagnostic scans propagate a record only once.
- `Toric3DVertexParityErasureEstimate` reports `shots`, `accepted`, `rejected`, `failures`, `accepted_failures`, `acceptance_rate`, `acceptance_ci`, and nullable `conditional_rate`/`conditional_ci`.

- [ ] **Step 1: Write failing decoder decision tests**

Build a small synthetic audit fixture and assert:

```julia
decision = decode_toric3d_vertex_circuit_z(decoder, syndrome, true)
@test decision.selected_class == (false, true)
@test decision.class_counts == (1, 3, 3, 0) # first maximum wins
@test TopoNoise._vertex_syndrome(star_checks(c), decision.full_recovery) == syndrome

fallback = decode_toric3d_vertex_circuit_z(decoder, unseen_syndrome, false)
@test fallback.fallback
@test fallback.full_recovery == fallback.base.full_recovery
```

Test construction rejection for mismatched dimensions/protocol and verify that all four logical differences have a zero-syndrome physical adjustment.

- [ ] **Step 2: Run the control test and verify RED**

Run:

```bash
env JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 julia --project=. --startup-file=no test/toric3d/vertex_projector_z_control.jl
```

Expected: FAIL because the decoder/control types and functions are undefined.

- [ ] **Step 3: Implement logical-class adjustment over baseline MWPM**

Build a dictionary from audit keys to observable groups. Precompute the four possible logical-Z support combinations and their two-bit pairing with declared logical-X representatives; require all four parities to be unique instead of assuming an index convention.

For known keys, call existing MWPM, compute its recovery class, XOR the precomputed zero-syndrome logical adjustment needed to reach the selected class, and validate the syndrome. For unknown keys, return the unchanged MWPM result with `fallback=true` and zero class counts. `finalize_toric3d_vertex_z_control` XORs the chosen recovery into `pre.frame.z` and reports residual syndrome and logical class without changing canonical replay.

- [ ] **Step 4: Write failing real-audit floor tests**

For `L=3,5,7`, enumerate final-layer signatures, build the control decoder, finalize every signature, and assert:

```julia
@test count(identity, logical_failures) == audit.unavoidable_count
@test audit.unavoidable_count == expected_floor
@test all(!any(out.residual_syndrome) for out in outcomes)
```

These tests must fail if the implementation still follows the baseline 34/131/276 logical failures.

- [ ] **Step 5: Write failing parity-erasure tests**

Pin total, accepted, rejected, unconditional failures, accepted failures, acceptance Wilson interval, and conditional Wilson interval. For all-odd input, require `conditional_rate === nothing` and `conditional_ci === nothing`.

- [ ] **Step 6: Complete the minimal implementation and verify GREEN**

Run the focused control and audit tests together. Expected: the circuit control reaches exactly 8/14/20 final-layer failures, the unknown-key fallback is unchanged, and parity estimates never emit NaN.

- [ ] **Step 7: Commit**

```bash
git add src/TopoNoise.jl src/toric3d/vertex_projector_z_control.jl test/toric3d/runtests.jl test/toric3d/vertex_projector_z_control.jl
git commit -m "feat: add circuit-aware H+5 Z control decoder"
```

### Task 4: Serialize immutable ambiguity artifacts

**Files:**
- Create: `src/toric3d/vertex_projector_z_output.jl`
- Create: `test/toric3d/vertex_projector_z_output.jl`
- Modify: `src/TopoNoise.jl`
- Modify: `test/toric3d/runtests.jl`

**Interfaces:**
- Consumes: audit/control types from Tasks 2-3 and existing `_vertex_write_json`, `_toric3d_atomic_write`, `_vertex_source_hashes`, `_vertex_versions` helpers.
- Produces: `toric3d_vertex_z_audit_payload(audits; provenance=Dict())::Dict{String,Any}`.
- Produces: `write_toric3d_vertex_z_ambiguity_audit(audits, dir; provenance=Dict())`, returning paths for canonical JSON, group CSV, and Markdown summary.
- Canonical format: `h5-z-boundary-ambiguity-v1` with deterministic geometry/protocol/group ordering.

- [ ] **Step 1: Write failing output tests**

Using an `L=3` final-layer audit, require:

- JSON has format version, equal-prior declaration, source hashes, backend versions, audit totals, parity totals, and every mixed group;
- group CSV row order is geometry, canonical protocol order, parity, syndrome indices;
- Markdown states both the avoidable baseline mismatch and unavoidable floor;
- a second write is byte-identical;
- extra provenance cannot overwrite scientific keys; and
- no serialized numeric value is NaN or infinity.

- [ ] **Step 2: Run the output test and verify RED**

Run the new Julia test file. Expected: FAIL because the payload/writer APIs are undefined.

- [ ] **Step 3: Implement deterministic atomic output**

Keep struct-to-dictionary conversion in `vertex_projector_z_output.jl`. Reuse the established atomic JSON and text writers. Serialize syndrome indices as one-based integer arrays and logical classes as two booleans. Put baseline and circuit-control failure counts beside the unavoidable floor. Reject empty audit lists, duplicate `(dimensions, protocol, final_layer_only)` identities, noncanonical protocols, and provenance collisions.

Include this file in `src/TopoNoise.jl` after `vertex_projector_output.jl`, where `_vertex_write_json`, `_vertex_source_hashes`, and backend-version helpers are already defined; export the payload and writer functions there.

The Markdown writer must label all results as `diagnostic_only`, say logical Z is not `|00_L>` infidelity, and say the control decoder cannot establish a deterministic threshold.

- [ ] **Step 4: Run output, audit, and historical output tests and verify GREEN**

Run:

```bash
env JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 julia --project=. --startup-file=no -e 'using Test; include("test/toric3d/vertex_projector_z_output.jl"); include("test/toric3d/vertex_projector_output.jl")'
```

Expected: PASS and no change to historical H+5 output schemas.

- [ ] **Step 5: Commit**

```bash
git add src/TopoNoise.jl src/toric3d/vertex_projector_z_output.jl test/toric3d/runtests.jl test/toric3d/vertex_projector_z_output.jl
git commit -m "feat: serialize H+5 Z ambiguity diagnostics"
```

### Task 5: Add the Python pre-fit ambiguity guard

**Files:**
- Modify: `examples/analyze_toric3d_vertex_projector.py:183-260`
- Create: `test/toric3d/test_vertex_projector_z_guard.py`

**Interfaces:**
- Produces: `load_z_boundary_audit(path)` with strict format and finite-value validation.
- Produces: `z_boundary_guard(audit, selected_sizes, protocol, source_hashes=None)::list[str]`.
- Extends: `fit_window(..., z_guard_issues=())`; nonempty issues are appended before crossing, optimization, or bootstrap work.
- Extends: `analyze(..., z_boundary_audit=None)`; only `z_only` candidates receive the new guard.
- Extends CLI with optional `--z-boundary-audit`; omission preserves historical behavior.

- [ ] **Step 1: Write failing guard tests**

Create a minimal valid audit fixture and test:

```python
with mock.patch.object(MODULE.engine, "_candidate_fit") as fit:
    result = MODULE.fit_window(rows, "z_only", "none", [.01, .05],
                               bootstrap=1000,
                               z_guard_issues=["Z boundary ambiguity: unavoidable=8"])
fit.assert_not_called()
self.assertEqual(result["status"], "diagnostic_candidate")
self.assertEqual(result["bootstrap"]["accepted"], 0)
self.assertEqual(result["bootstrap"]["samples"], [])
```

Also test missing size/protocol entries, nonzero mixed/unavoidable/control-failure counts, source-hash mismatch, malformed class counts, and a clean certificate that allows the existing fitter to run. Confirm X-only fits ignore the Z guard and existing calls without a certificate are unchanged.

- [ ] **Step 2: Run the Python guard test and verify RED**

Run:

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .CondaPkg/.pixi/envs/default/bin/python -m unittest discover -s test/toric3d -p 'test_vertex_projector_z_guard.py' -v
```

Expected: FAIL because the loader, guard, and new arguments do not exist.

- [ ] **Step 3: Implement validation and early return**

Validate `h5-z-boundary-ambiguity-v1`, unique audit identities, canonical protocols, requested sizes, nonnegative integer counts, class-count sums, finite rates, and optional source-hash equality. Return explicit issues with geometry/protocol/count values.

In `fit_window`, add guard issues before the existing four-size/rate checks and return immediately on any issue. In `analyze`, compute protocol-specific issues once for `z_only`. Parse the optional CLI file and include its SHA-256 in final analysis hashes; do not alter output when the option is absent.

- [ ] **Step 4: Run new and historical Python tests and verify GREEN**

Run:

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .CondaPkg/.pixi/envs/default/bin/python -m unittest discover -s test/toric3d -p 'test_vertex_projector_z_guard.py' -v
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .CondaPkg/.pixi/envs/default/bin/python -m unittest discover -s test/toric3d -p 'test_vertex_projector_fit.py' -v
```

Expected: PASS; the new test proves zero optimizer/bootstrap calls on ambiguity, and every historical H+5 fit test remains unchanged.

- [ ] **Step 5: Commit**

```bash
git add examples/analyze_toric3d_vertex_projector.py test/toric3d/test_vertex_projector_z_guard.py
git commit -m "feat: guard H+5 Z fits on boundary ambiguity"
```

### Task 6: Generate the approved diagnostic evidence

**Files:**
- Create: `results/toric3d-h5-z-boundary-ambiguity-2026-10-14/run_audit.jl`
- Create: `results/toric3d-h5-z-boundary-ambiguity-2026-10-14/apply_guard.py`
- Create/generated: `results/toric3d-h5-z-boundary-ambiguity-2026-10-14/z-boundary-audit.json`
- Create/generated: `results/toric3d-h5-z-boundary-ambiguity-2026-10-14/z-boundary-groups.csv`
- Create/generated: `results/toric3d-h5-z-boundary-ambiguity-2026-10-14/z-boundary-report.md`
- Create: `results/toric3d-h5-z-boundary-ambiguity-2026-10-14/run_control.jl`
- Create/generated: `results/toric3d-h5-z-boundary-ambiguity-2026-10-14/parity-erasure-control.csv`
- Create/generated: `results/toric3d-h5-z-boundary-ambiguity-2026-10-14/guard.json`
- Create: `results/toric3d-h5-z-boundary-ambiguity-2026-10-14/runtime-environment.json`
- Create: `results/toric3d-h5-z-boundary-ambiguity-2026-10-14/README.md`

**Interfaces:**
- Consumes: the public audit/output APIs and Python guard from Tasks 2-5.
- Produces: a separate immutable diagnostic result set for final-layer `L=3,5,7`, all four protocols, plus a paired parity-erasure control at `p=10^-4,10^-3`, with no threshold fit.

- [ ] **Step 1: Write the reproducible Julia audit driver**

The driver builds one decoder family per size with declared `p=0.01`, audits only final-layer single Z faults for every canonical protocol, requires cross-protocol signature equality, builds the control decoder, and asserts the approved counts before writing artifacts atomically.

The output provenance must include spec commit `cb2ec97`, the implementation commit IDs available at runtime, single-thread environment, equal-prior assumption, and links to the completed preliminary result directory. It must not open or write `document/note.tex`.

- [ ] **Step 2: Run the driver and inspect all measured values**

Run with the global single-thread environment. Expected for each protocol:

| L | locations | mixed groups | baseline failures | circuit-control/unavoidable failures | even-sector mixed groups |
|---:|---:|---:|---:|---:|---:|
| 3 | 99 | 5 | 34 | 8 | 0 |
| 5 | 275 | 9 | 131 | 14 | 0 |
| 7 | 539 | 13 | 276 | 20 | 0 |

Any mismatch stops artifact creation and returns to diagnosis; do not edit expected values to match unexpected output.

- [ ] **Step 3: Run the paired parity-erasure control**

`run_control.jl` uses `Random.Xoshiro(1234)` to assign deterministic point seeds, 1000 paired Z-only records per `(L,p)` point, cubic `L=3,5,7`, and `p in (1e-4,1e-3)`. For every record it replays each canonical protocol once to the pre-final state, then obtains both baseline and circuit-control outcomes from `finalize_toric3d_vertex_z_control`.

Write one CSV row per size/rate/protocol with total shots, odd/even counts, unconditional baseline/control failures, accepted baseline/control failures, acceptance rate/interval, and conditional rates/intervals. Label every row `parity-erasure diagnostic`; missing conditional estimates use empty CSV fields. This run has no threshold fitting and does not reuse or overwrite any prior checkpoint.

- [ ] **Step 4: Write and run the Python guard driver**

`apply_guard.py` loads the canonical audit through the production loader, requests 1000 bootstrap replicates for each protocol, and writes `guard.json` showing `status=diagnostic_only`, the exact mixed/unavoidable/control-failure reasons, `bootstrap_requested=1000`, and `bootstrap_executed=0`. It must not construct synthetic threshold estimates.

- [ ] **Step 5: Complete the README and inspect artifacts**

Document:

- the distinction between current MWPM failures and the optimal observable floor;
- why exact upper outcome bits reveal only parity under a uniform even clean branch;
- why odd-parity erasure is a finite-size/postselected control;
- why logical Z is not `|00_L>` infidelity; and
- why a deterministic positive threshold requires a later physical-boundary or dimensional redesign.

Check JSON/audit CSV/control CSV/Markdown consistency with an independent short Python verification embedded in `apply_guard.py`. No LaTeX note or PDF work belongs in this task.

- [ ] **Step 6: Commit the result set**

```bash
git add results/toric3d-h5-z-boundary-ambiguity-2026-10-14
git commit -m "data: report H+5 Z-boundary ambiguity"
```

### Task 7: Run complete verification and preserve evidence

**Files:**
- Create: `results/toric3d-h5-z-boundary-ambiguity-2026-10-14/verification.md`
- Modify only if verification reveals a defect: files owned by Tasks 1-6, with a new failing regression test before each fix.

**Interfaces:**
- Consumes: all implementation and artifact tasks.
- Produces: final regression evidence, immutable-note proof, and a clean local branch.

- [ ] **Step 1: Run the complete Julia regression suite**

Run:

```bash
env JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 julia --project=. --startup-file=no test/runtests.jl
```

Expected: exit 0 with no failed/error/broken tests. Record the assertion count and elapsed time.

- [ ] **Step 2: Run the complete Python regression suite**

Run:

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .CondaPkg/.pixi/envs/default/bin/python -m unittest discover -s test/toric3d -p 'test_*.py' -v
```

Expected: exit 0. Record test count and elapsed time.

- [ ] **Step 3: Reverify the completed historical delivery**

Run:

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .CondaPkg/.pixi/envs/default/bin/python results/toric3d-h5-diagonal-feedback-2026-10-13/preliminary-l3-7/verify_delivery.py
```

Expected: `PASS: final measurements, guards, provenance, ... verified` with no historical artifact changes.

- [ ] **Step 4: Verify source, artifact, and note immutability**

Run:

```bash
git diff --check
shasum -a 256 document/note.tex
git status --short
```

Expected note hash: `3a40710c2ab6bd43e53f5f7703b7c8716b209dc1e060a5092f17e0edce073f3a`. Inspect `git diff --name-only cb2ec97..HEAD` and require that neither `document/note.tex` nor any prior H+5 result artifact appears.

- [ ] **Step 5: Record verification and commit**

Write exact commands, exit codes, counts, artifact SHA-256 values, and note hash to `verification.md`. Then commit:

```bash
git add results/toric3d-h5-z-boundary-ambiguity-2026-10-14/verification.md
git commit -m "docs: verify H+5 Z-boundary controls"
```

- [ ] **Step 6: Final branch check**

Require a clean worktree, list local commits since `cb2ec97`, and report that nothing was pushed, merged, or published. Present the audit JSON, CSV, Markdown report, guard output, and verification file to the user, with the conclusion labeled diagnostic rather than a threshold.
