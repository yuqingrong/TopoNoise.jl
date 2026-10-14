# Diagonal H+5 3D Toric-Code Feedback Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build, verify, simulate, and document the direct H+5 diagonal 3D toric-code preparation with paired `none`, `one_step`, `window_complete`, and `full_prefix` feedback for separate X-only and Z-only circuit noise.

**Architecture:** Add an isolated H+5 subsystem that reuses the final `Toric3DEncodedSlice` cell complex but owns its boundary instrument, six-gate schedule, immutable replay record, feedback windows, decoders, four-arm statistics, and artifacts. Keep every H+6 plaquette-transfer type, operation list, public function, CSV schema, and historical result unchanged. Apply all corrections to binary Pauli frames, validate complete recoveries against independently recomputed syndromes, and preserve partial-route residuals explicitly.

**Tech Stack:** Julia 1.10+, `TopoNoise`, `PythonCall`, PyMatching 2.4.0, `ldpc` 2.4.1 BP+OSD, CairoMakie, TOML/CSV, Python 3.10+ with SciPy/Matplotlib, LaTeX.

**Spec:** `docs/superpowers/specs/2026-10-13-toric3d-diagonal-h5-vertex-projector-design.md`

## Global Constraints

- Public H+5 dimensions are `(T, Lu, Lv)` with `T >= 1` and `Lu,Lv >= 2`; legacy `Ly,Lz` fields may remain internal compatibility details.
- Use the existing final code with `n=3(T+1)LuLv`, `k=2`, `dZ=2min(Lu,Lv)`, and `dX=(T+1)min(Lu,Lv)`.
- Each cell has exactly `H(n1)`, `CX(n1,n2)`, `CX(n1,n3)`, `CX(n1,e1)`, `CX(n1,e2)`, `CX(n1,e3)` in that order: six operations and eleven noisy participating outputs.
- Use independent post-output rate-`p` X-only or Z-only faults; initialization, idle qubits, boundary operations, measurements, classical processing, and corrections are ideal.
- Sample one immutable gate-fault record and one clean even upper-boundary branch before decoding; replay both through all four protocols.
- Intermediate feedback runs only after layers `1:T-1`; no decoder may inspect or correct `E_{t+1}:E_T`.
- Use minimum-sum BP, 50 iterations, scaling factor 1, parallel schedule, one thread, and `OSD_CS` order 2; record declared `p` and `p_eff=clamp(p,1e-12,1-1e-12)`.
- Preserve distinct period-two physical edges; matching minimizes unit physical length first and the canonical lexicographic matching/path signature second. Pin PyMatching 2.4.0 and `ldpc` 2.4.1 in `CondaPkg.toml` and record the loaded versions.
- A complete decoder recovery must satisfy `H*c == s`; a partial one-step route must satisfy `s_after == s_before xor H*r` and its protocol-specific residual invariant.
- Score both transverse logical-X and logical-Z pairs after final recovery. Only residual logical X may also be labeled `|00_L>` preparation failure.
- Do not modify or rename `Toric3DSliceEncoder`; label it and its saved outputs as the H+6 *plaquette-transfer* circuit.
- Scientific runs use seed 1234, atomic point checkpoints, and zero tolerated invalid recoveries.

## Review Focus

- Period-two and asymmetric transverse sizes `(2,2,3)` and `(2,3,2)` must retain parallel edge identities and decoder validity; Task 4 owns these tests.
- Matching ties must select the declared canonical physical paths, rather than backend-dependent pairs or merged incidence columns; Task 5 owns these tests.
- The sampled clean upper branch `b` must remain hidden from decoders, paired across protocols, and included in resume draw order; Tasks 2 and 9 own these tests.
- Partial one-step routes must not be subjected to the complete-recovery equation, while every full path and complete recovery must be; Tasks 5–7 own these tests.
- A syndrome-valid final recovery can leave a zero-syndrome logical operator, so state preparation and decoder validity need separate assertions; Tasks 7 and 8 own this regression.

---

### Task 1: Freeze the legacy circuit and add the H+5 code facade and encoder

**Files:**
- Create: `src/toric3d/vertex_projector_geometry.jl`
- Create: `src/toric3d/vertex_projector_encoder.jl`
- Create: `test/toric3d/vertex_projector_legacy.jl`
- Create: `test/toric3d/vertex_projector_encoder.jl`
- Modify: `src/TopoNoise.jl`
- Modify: `test/toric3d/runtests.jl`

**Interfaces:**
- Consumes: `Toric3DEncodedSlice`, `Toric3DSliceLattice`, `EncoderOperation`, existing check/logical-support functions.
- Produces: `Toric3DVertexProjectorCode(T,Lu,Lv)`, `Toric3DVertexProjectorEncoder`, `toric3d_vertex_projector_encoder(code)`, `verify_toric3d_vertex_projector_encoder(encoder)::Bool`, and public `(T,Lu,Lv)` dimension metadata.

**Testsets:** `"Legacy H+6 plaquette-transfer schedule is frozen"`, `"Legacy H+6 replay and artifact formats remain reproducible"`, `"H+5 target geometry and logical algebra"`, `"H+5 direct-projector cells use the exact six-operation schedule"`, `"H+5 layers are fresh, disjoint, ordered, and causal"`, `"H+5 clean prefixes prepare exactly the completed checks"`, and `"The H+6 honeycomb seed is invalid for H+5"`.

- [ ] **Step 1: Add a passing legacy characterization test**

  Pin the complete `(T,Lu,Lv)=(1,2,2)` H+6 operation/checkpoint signature, thirteen active outputs per cell, a seeded replay result, and the current three-arm CSV header in testsets `"Legacy H+6 plaquette-transfer schedule is frozen"` and `"Legacy H+6 replay and artifact formats remain reproducible"`.

- [ ] **Step 2: Run the legacy characterization**

  Run:
  `JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. -e 'using Test, Random, TopoNoise; include("test/toric3d/vertex_projector_legacy.jl")'`

  Expected: PASS before H+5 source files are included.

- [ ] **Step 3: Write failing H+5 geometry and encoder tests**

  Test `n`, `k=2`, both logical pairs, commutation, distances, and exact gate order/support for every cell. For `(2,2,2)`, assert 48 operations, 88 active outputs, checkpoints `[24,48]`, increasing `t,u,v` order, within-layer edge-disjoint cells, and fresh outgoing edges. At every checkpoint `t`, assert the exact prepared set: lower white triads, bulk stars `A_0:A_(t-1)`, plaquettes `P_1:P_t`, both logical-Z cycles, and no upper/frontier triad. Add an explicit test that the old honeycomb seed is not a valid drop-in H+5 boundary.

- [ ] **Step 4: Run the H+5 encoder test and confirm the missing interface**

  Run: `JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. test/toric3d/vertex_projector_encoder.jl`.

  Expected: FAIL because `Toric3DVertexProjectorCode` and its encoder are undefined.

- [ ] **Step 5: Implement the code facade and direct projector schedule**

  `Toric3DVertexProjectorCode` wraps the existing final cell complex while exposing `T,Lu,Lv`. Delegate qubit indexing, final checks, and logical supports without copying the H+6 encoder. Store one operation vector and one end-of-layer checkpoint vector in `Toric3DVertexProjectorEncoder`.

- [ ] **Step 6: Export/include the new interfaces and run both focused suites**

  Run: `JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. -e 'using Test, Random, TopoNoise; include("test/toric3d/vertex_projector_legacy.jl"); include("test/toric3d/vertex_projector_encoder.jl")'`.

  Expected: both suites PASS; the H+6 signature remains byte-for-byte unchanged.

- [ ] **Step 7: Commit**

  ```bash
  git add src/TopoNoise.jl src/toric3d/vertex_projector_geometry.jl src/toric3d/vertex_projector_encoder.jl test/toric3d/runtests.jl test/toric3d/vertex_projector_legacy.jl test/toric3d/vertex_projector_encoder.jl
  git commit -m "feat: add diagonal H+5 vertex projector encoder"
  ```

### Task 2: Implement the ideal lower preparation and paired upper-boundary instrument

**Files:**
- Create: `src/toric3d/vertex_projector_boundary.jl`
- Create: `test/toric3d/vertex_projector_boundary.jl`
- Modify: `src/TopoNoise.jl`
- Modify: `test/toric3d/runtests.jl`

**Interfaces:**
- Consumes: `Toric3DVertexProjectorCode`, `PauliFrame`, final star/logical supports.
- Produces: `Toric3DVertexBoundaryMap`, `toric3d_vertex_boundary_map(code)` rooted at `A_(T-1,0,0)`, `verify_toric3d_vertex_lower_boundary(map)::Bool`, `toric3d_vertex_boundary_correction(map,a)::BitVector`, `sample_toric3d_vertex_upper_branch(rng,map)::BitVector`, `Toric3DVertexBoundaryCompletion`, and the simulator-only `complete_toric3d_vertex_boundary(map,z_frame,b)`.

**Testsets:** `"H+5 lower projection and upper outcomes have the specified ranks"`, `"H+5 upper completion realizes every clean even branch"`, `"H+5 paired boundary map equals the explicit two-branch instrument"`, `"H+5 tree corrections have zero homology and preserve all logical classes"`, and `"H+5 clean upper branch sampling is seeded and even"`.

- [ ] **Step 1: Write failing boundary-algebra tests**

  Assert independent lower white triads, rank exactly `M-1` for the clean upper-outcome distribution, uniform support over all `2^(M-1)` even-parity branches, a physical-edge-ID-ordered BFS spanning tree excluding the union of both logical-X supports, rooted specifically at final bulk star `A_(T-1,0,0)`, and zero logical-Z homology for every tree correction. On `Lu=Lv=2`, exhaust all eight clean even branches and all sixteen upper flip vectors.

- [ ] **Step 2: Write failing paired-instrument tests**

  Require the physical observed branch `a=b xor u`, an observable-only controller call `toric3d_vertex_boundary_correction(map,a)` that receives neither `b` nor the frame flip `u`, storage of `C(a)` and `C(b)` in simulator provenance, and equality of `z xor C(a) xor C(b)` with the optimized relative map. Inject each of the four logical Pauli generators immediately before completion and require its relative class to be preserved.

- [ ] **Step 3: Run the focused test and confirm failure**

  Run: `JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. test/toric3d/vertex_projector_boundary.jl`.

  Expected: FAIL on the missing boundary types and functions.

- [ ] **Step 4: Implement the deterministic tree map and branch sampler**

  Represent the ideal lower-white-triad measurement and its deterministic zero-homology Z correction as the noiseless reference preparation. Use physical-edge-ID-ordered BFS on the star graph with the union of the two logical-X supports excluded, root it at `A_(T-1,0,0)`, and record that convention and root in provenance. The excluded-edge condition, not acyclicity alone, guarantees zero logical overlaps. The upper sampler draws the first `Lu*Lv-1` bits and fixes the final parity. The controller accepts only the observed physical outcome; keep clean `b` solely in replay provenance and the explicit clean/noisy comparison.

- [ ] **Step 5: Implement stabilizer-state verification helpers**

  Verify lower stars, every completed bulk star, all plaquettes, and both logical Z cycles before upper completion; after each clean branch require every final target generator and logical Z.

- [ ] **Step 6: Run the boundary and encoder suites**

  Run: `JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. -e 'using Test, Random, TopoNoise; include("test/toric3d/vertex_projector_encoder.jl"); include("test/toric3d/vertex_projector_boundary.jl")'`.

  Expected: PASS for all allowed branches, logical injections, and seeded sampling.

- [ ] **Step 7: Commit**

  ```bash
  git add src/TopoNoise.jl src/toric3d/vertex_projector_boundary.jl test/toric3d/runtests.jl test/toric3d/vertex_projector_boundary.jl
  git commit -m "feat: add paired H+5 boundary completion"
  ```

### Task 3: Add immutable fault records and frame propagation

**Files:**
- Create: `src/toric3d/vertex_projector_noise.jl`
- Create: `test/toric3d/vertex_projector_faults.jl`
- Modify: `src/TopoNoise.jl`
- Modify: `test/toric3d/runtests.jl`

**Interfaces:**
- Consumes: H+5 encoder operations and upper-branch sampler.
- Produces: immutable `Toric3DVertexProjectorFaultRecord` with fields `(T,Lu,Lv,x,z,upper_branch)`, `sample_toric3d_vertex_projector_record(rng,encoder;p_x,p_z)`, `with_toric3d_vertex_projector_fault(record,encoder,step,qubit,pauli)`, and `propagate_toric3d_vertex_projector_segment!(frame,encoder,record,first,last)`.

**Testsets:** `"H+5 records enumerate exactly eleven outputs per vertex"`, `"H+5 frame propagation equals a physical Pauli reference"`, and `"H+5 replay records are immutable and channel separated"`.

- [ ] **Step 1: Write failing record-shape and sampling tests**

  Assert exactly eleven eligible outputs per cell, inactive second H slots, `p=0`, `p=1`, independent channels, repeat seeds, immutable arrays, and one stored even `upper_branch` per record.

- [ ] **Step 2: Write the independent physical-Pauli reference test**

  On `(2,2,2)`, inject X and Z at each of the 88 active output locations and compare binary-frame propagation with direct Clifford conjugation. Assert no replay mutates or resamples the record.

- [ ] **Step 3: Run the test and confirm missing interfaces**

  Run: `JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. test/toric3d/vertex_projector_faults.jl`.

  Expected: FAIL.

- [ ] **Step 4: Implement validation, sampling, fault insertion, and segment propagation**

  Use a `2 × operation_count` bit matrix per Pauli, validate every active slot from the stored operation, and sample `upper_branch` after gate faults in a fixed documented RNG order.

- [ ] **Step 5: Run fault, boundary, and encoder suites**

  Run: `JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. -e 'using Test, Random, TopoNoise; include("test/toric3d/vertex_projector_encoder.jl"); include("test/toric3d/vertex_projector_boundary.jl"); include("test/toric3d/vertex_projector_faults.jl")'`.

  Expected: PASS, including all 176 single-location/channel reference cases.

- [ ] **Step 6: Commit**

  ```bash
  git add src/TopoNoise.jl src/toric3d/vertex_projector_noise.jl test/toric3d/runtests.jl test/toric3d/vertex_projector_faults.jl
  git commit -m "feat: add paired H+5 circuit faults"
  ```

### Task 4: Build all H+5 feedback windows and restricted matrices

**Files:**
- Create: `src/toric3d/vertex_projector_feedback_geometry.jl`
- Create: `test/toric3d/vertex_projector_feedback_geometry.jl`
- Modify: `src/TopoNoise.jl`
- Modify: `test/toric3d/runtests.jl`

**Interfaces:**
- Consumes: code facade, `E_t` qubit slices, final checks.
- Produces: `Toric3DVertexFeedbackWindow` with `channel`, `policy`, `checkpoint`, `checks`, `qubits`, `matrix`, `retired_rows`, and `frontier_rows`; constructors `toric3d_vertex_x_slice_window`, `toric3d_vertex_x_complete_window`, `toric3d_vertex_x_prefix_window`, `toric3d_vertex_z_strip_window`, `toric3d_vertex_z_complete_window`, and `toric3d_vertex_z_prefix_window`.

**Testsets:** `"H+5 Q slices are target stabilizers with even syndrome"`, `"H+5 X and Z windows use prepared causal supports"`, `"H+5 recovery matrices have the specified ranks"`, and `"H+5 period-two geometry preserves physical edge identities"`.

- [ ] **Step 1: Write failing X-window tests**

  Verify `Q_t` stabilizer membership, every `E_t` edge incident on two `Q_t` checks, even syndromes, `P_1` at `t=1`, `P_{t-1} union P_t` later, recovery support `E_{t-1} union E_t`, and no future columns.

- [ ] **Step 2: Write failing Z-window tests**

  Verify `(W,A_0)` at `t=1`, `(A_{t-2},A_{t-1})` later, the `R-F-B` strip labels, `W` as a detector rather than boundary, open future columns, and prefix checks excluding unprepared upper triads.

- [ ] **Step 3: Pin ranks and period-two edge identities**

  For `M=LuLv`, assert X ranks `3M-1` for the `3M × 6M` first window and `4M` for steady `6M × 6M`; assert full row rank `2M` for both Z `2M × 6M` windows. Cover `(2,2,2)`, `(2,2,3)`, `(2,3,2)`, and `(3,3,3)` and compare physical edge IDs, not coordinates alone.

- [ ] **Step 4: Run the test and confirm failure**

  Run: `JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. test/toric3d/vertex_projector_feedback_geometry.jl`.

  Expected: FAIL on missing constructors.

- [ ] **Step 5: Implement the common restricted-window representation and constructors**

  Mark retiring and frontier row sets explicitly so later validators do not infer them from row positions.

- [ ] **Step 6: Run geometry, encoder, and feedback-geometry suites**

  Run: `JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. -e 'using Test, Random, TopoNoise; include("test/toric3d/vertex_projector_encoder.jl"); include("test/toric3d/vertex_projector_feedback_geometry.jl")'`.

  Expected: PASS for all ranks, supports, and asymmetric period-two cases.

- [ ] **Step 7: Commit**

  ```bash
  git add src/TopoNoise.jl src/toric3d/vertex_projector_feedback_geometry.jl test/toric3d/runtests.jl test/toric3d/vertex_projector_feedback_geometry.jl
  git commit -m "feat: define H+5 feedback windows"
  ```

### Task 5: Implement deterministic MWPM paths and one-step transport

**Files:**
- Create: `src/toric3d/vertex_projector_matching.jl`
- Create: `test/toric3d/vertex_projector_x_one_step.jl`
- Create: `test/toric3d/vertex_projector_z_one_step.jl`
- Modify: `CondaPkg.toml`
- Modify: `src/TopoNoise.jl`
- Modify: `test/toric3d/runtests.jl`

**Interfaces:**
- Consumes: one-step X and Z windows.
- Produces: `Toric3DVertexMatchedPath`, `Toric3DVertexMatchingDecision` with named path/recovery/route/residual/signature/backend/version fields, `build_toric3d_vertex_x_one_step_decoder`, `build_toric3d_vertex_z_one_step_decoder`, `decode_toric3d_vertex_x_one_step`, and `decode_toric3d_vertex_z_one_step`.

**Testsets:** `"H+5 one-step X MWPM returns canonical physical paths"`, `"H+5 one-step X routes predict their residual syndrome"`, `"H+5 old-target X faults remain terminal local noise"`, `"H+5 period-two X decoding exhausts every reachable syndrome"`, `"H+5 one-step Z matching uses exactly one temporal mate for odd parity"`, `"H+5 causal Z trimming clears retiring sheets and confines residuals"`, `"Symmetric Z endpoint trimming is causally invalid"`, and `"H+5 period-two Z decoding exhausts every reachable strip syndrome"`.

- [ ] **Step 1: Write exhaustive period-two X tests**

  Exhaust all eight even `Q_t` syndromes for `Lu=Lv=2`. Pin tied examples to the lexicographically smallest matching/path signature. Require `H*full_recovery=s`, length-one annihilation, longer first/last endpoint motion, and `s_after=s_before xor H*route`.

- [ ] **Step 2: Write terminal old-target X tests**

  Inject old-target faults outside `Q_t`, show that the one-step measurement does not see them, and prove those qubits undergo no later H+5 gate and cannot proliferate.

- [ ] **Step 3: Write exhaustive period-two Z-strip tests**

  Exhaust all 256 syndromes of the full-rank `2M` strip. Even syndromes use no dummy; odd syndromes use exactly one future-boundary mate. Pin a literal odd case whose dummy distance and lifted support equal the shortest periodic `R -> F` physical path plus exactly one `E_t` edge, ruling out a zero-cost ordinary boundary. Pin a separate equal-primary-weight strip case to the exact lexicographically smallest matching/path signature. Require complete path validity, route support only on `E_{t-1}`, zero retiring rows, and residuals confined to measured frontier rows.

- [ ] **Step 4: Add the symmetric-trim counterexample**

  Pin one path where applying both endpoint edges moves a frontier defect backward into the retiring sheet and leaves it unreachable at the next checkpoint.

- [ ] **Step 5: Run both tests and confirm missing decoders**

  Run: `JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. -e 'using Test, Random, TopoNoise; include("test/toric3d/vertex_projector_x_one_step.jl"); include("test/toric3d/vertex_projector_z_one_step.jl")'`.

  Expected: FAIL.

- [ ] **Step 6: Implement canonical matching, physical path lifting, and protocol trims**

  Pin PyMatching 2.4.0 and `ldpc` 2.4.1 in `CondaPkg.toml`. Construct deterministic graph inputs from physical edge IDs, preserve parallel period-two edges, optimize unit length before the ordered signature, and record the backend/version and chosen signature. Validate the full path before deriving the partial route.

- [ ] **Step 7: Run one-step and geometry suites**

  Run: `JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. -e 'using Test, Random, TopoNoise; include("test/toric3d/vertex_projector_feedback_geometry.jl"); include("test/toric3d/vertex_projector_x_one_step.jl"); include("test/toric3d/vertex_projector_z_one_step.jl")'`.

  Expected: PASS for exhaustive X/Z syndrome classes and the causal trim regression.

- [ ] **Step 8: Commit**

  ```bash
  git add CondaPkg.toml src/TopoNoise.jl src/toric3d/vertex_projector_matching.jl test/toric3d/runtests.jl test/toric3d/vertex_projector_x_one_step.jl test/toric3d/vertex_projector_z_one_step.jl
  git commit -m "feat: add causal H+5 one-step feedback"
  ```

### Task 6: Add complete-window, full-prefix, and final decoders

**Files:**
- Create: `src/toric3d/vertex_projector_decoders.jl`
- Create: `test/toric3d/vertex_projector_complete_recovery.jl`
- Modify: `src/TopoNoise.jl`
- Modify: `test/toric3d/runtests.jl`

**Interfaces:**
- Consumes: all restricted windows, canonical matching, `Toric3DDecoderSettings`.
- Produces: `Toric3DVertexBPOSDDecoder`, `Toric3DVertexDecoderFamily`, `build_toric3d_vertex_decoders(code;p,settings)`, `decode_toric3d_vertex_bposd`, complete Z matching decode, and common final X/Z recovery methods.

**Testsets:** `"H+5 window-complete X recovery clears every completed plaquette"`, `"H+5 window-complete Z recovery clears every completed star sheet"`, `"H+5 full-prefix recovery never accesses future data"`, `"H+5 combined-channel decisions share one pre-correction frame and commute"`, and `"H+5 BP+OSD priors and every complete recovery are validated"`.

- [ ] **Step 1: Write failing X complete-recovery tests**

  Test first and steady windows, every single output fault, randomized multifault histories at transverse periods two and three, full-prefix support, and zero completed plaquette syndrome after recovery.

- [ ] **Step 2: Write failing Z complete-recovery tests**

  Test ordinary open-boundary matching with any number of nonzero-cost boundary terminations, window and prefix supports, odd syndromes, and zero completed star syndrome after recovery. Pin literal equal-primary-weight examples for window, prefix, and final Z recovery and require the exact lexicographically smallest physical matching/path signature in every case.

- [ ] **Step 3: Write decoder-configuration and validation tests**

  Pin minimum-sum/50/1/parallel/one-thread/`OSD_CS`-2, record `p` and `p_eff`, and test endpoint clamping at `p=0` and `p=1`. Flip one returned bit in every complete decoder family and require a hard invalid-recovery error.

- [ ] **Step 4: Run the test and confirm failure**

  Run: `JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. test/toric3d/vertex_projector_complete_recovery.jl`.

  Expected: FAIL.

- [ ] **Step 5: Implement reusable restricted BP+OSD and open-MWPM builders**

  Build X decoders per distinct restricted matrix and declared prior. Build Z window/prefix/final matchers with the same physical-path and deterministic-signature rules as Task 5. Record the neglected circuit correlations in decoder metadata.

- [ ] **Step 6: Implement the decoder family cache**

  Cache per-checkpoint one-step/window/prefix decoders and common final decoders. Reject any recovery accessing an unlisted support column.

- [ ] **Step 7: Run complete-recovery, one-step, and existing decoder suites**

  Run: `JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. -e 'using Test, Random, TopoNoise; include("test/toric3d/vertex_projector_x_one_step.jl"); include("test/toric3d/vertex_projector_z_one_step.jl"); include("test/toric3d/vertex_projector_complete_recovery.jl"); include("test/toric3d/decoder.jl")'`.

  Expected: PASS with every full recovery independently reproducing its syndrome.

- [ ] **Step 8: Commit**

  ```bash
  git add src/TopoNoise.jl src/toric3d/vertex_projector_decoders.jl test/toric3d/runtests.jl test/toric3d/vertex_projector_complete_recovery.jl
  git commit -m "feat: add complete H+5 feedback decoders"
  ```

### Task 7: Integrate four-protocol replay, boundary completion, and logical scoring

**Files:**
- Create: `src/toric3d/vertex_projector_replay.jl`
- Create: `test/toric3d/vertex_projector_replay.jl`
- Modify: `src/TopoNoise.jl`
- Modify: `test/toric3d/runtests.jl`

**Interfaces:**
- Consumes: encoder, immutable record, boundary map, decoder family.
- Produces: `TORIC3D_VERTEX_PROTOCOLS = (:none,:one_step,:window_complete,:full_prefix)`, `TORIC3D_VERTEX_REPLAY_CHANNELS = (:x_only,:z_only,:both)`, `Toric3DVertexFeedbackHistory`, `Toric3DVertexReplayOutcome`, and `replay_toric3d_vertex_projector(encoder,record,decoders;protocol,channel)`.

**Testsets:** `"H+5 four-arm replay shares one immutable fault and boundary record"`, `"H+5 exhaustive single-output faults have validated histories"`, `"H+5 zero noise prepares the target state in both channels"`, `"H+5 final recovery reports state-independent logical Pauli classes"`, and `"H+5 prepared-state verification is distinct from syndrome-valid decoding"`.

- [ ] **Step 1: Write failing four-arm replay tests**

  Assert identical record and clean branch use across all arms, feedback only after `1:T-1`, no intermediate feedback for `T=1`, and exact protocol-specific check/correction supports.

- [ ] **Step 2: Write final scoring tests**

  Require two logical-X bits, two logical-Z bits, their union rates, and the full four-bit Pauli class. Label `state_failure` only for X-only `|00_L>` preparation, while retaining logical Z as a state-independent diagnostic.

- [ ] **Step 3: Write combined-channel and validation tests**

  Accept exactly `channel=:x_only`, `:z_only`, or `:both`; reject every other symbol before propagation. A pure-channel call must reject a record with any fault in the opposite Pauli component. For `:both`, capture X and Z measurements from one pre-correction frame, calculate both decisions, and require both correction orders to yield the same binary frame for `one_step` and `window_complete`. Insert a zero-syndrome logical operator to show that recovery validity does not prove the intended logical state.

- [ ] **Step 4: Run the test and confirm failure**

  Run: `JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. test/toric3d/vertex_projector_replay.jl`.

  Expected: FAIL.

- [ ] **Step 5: Implement layer replay and feedback history**

  Propagate each layer once, apply the selected route before the next layer, independently remeasure after every route, store complete paths/recoveries and residuals, then apply the explicit paired upper completion and common final recovery.

- [ ] **Step 6: Run replay and decoder suites**

  Run: `JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. -e 'using Test, Random, TopoNoise; include("test/toric3d/vertex_projector_complete_recovery.jl"); include("test/toric3d/vertex_projector_replay.jl")'`.

  Expected: zero-noise success for both channels and valid histories for all protocols.

- [ ] **Step 7: Commit**

  ```bash
  git add src/TopoNoise.jl src/toric3d/vertex_projector_replay.jl test/toric3d/runtests.jl test/toric3d/vertex_projector_replay.jl
  git commit -m "feat: replay four H+5 feedback protocols"
  ```

### Task 8: Add exhaustive small-system and randomized reference verification

**Files:**
- Create: `test/toric3d/vertex_projector_reference.jl`
- Modify: `test/toric3d/runtests.jl`

**Interfaces:**
- Consumes: every H+5 interface through replay.
- Produces: an integration verification gate; no new public API.

**Testsets:** `"H+5 exhaustive single faults satisfy every checkpoint invariant"` and `"H+5 randomized period-two histories match physical Pauli evolution"`.

- [ ] **Step 1: Add exhaustive single-fault replay**

  On `(T,Lu,Lv)=(2,2,2)`, run all 88 active output locations for both Paulis through every applicable arm. Require every complete intermediate/final recovery equation and every partial-route residual equation.

- [ ] **Step 2: Add randomized multifault multi-checkpoint histories**

  Use fixed seeds for `(3,2,2)`, `(3,2,3)`, and `(3,3,2)`; assert causality, prefix-clearing invariants, paired upper completion, and physical/reference-frame equality.

- [ ] **Step 3: Run the new integration suite**

  Run: `JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. test/toric3d/vertex_projector_reference.jl`.

  Expected: PASS. Any failure must report the exact gate output, checkpoint, physical support, syndrome before/after, route, and branch.

- [ ] **Step 4: Run the full toric3d suite**

  Run:
  `JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. -e 'using Test, Random, TopoNoise; include("test/toric3d/runtests.jl")'`

  Expected: PASS, including all H+6 and Cartesian 3D tests.

- [ ] **Step 5: Commit**

  ```bash
  git add test/toric3d/runtests.jl test/toric3d/vertex_projector_reference.jl
  git commit -m "test: exhaust H+5 feedback fault locations"
  ```

### Task 9: Add four-arm paired statistics, atomic checkpointing, and resume

**Files:**
- Create: `src/toric3d/vertex_projector_scan.jl`
- Create: `test/toric3d/vertex_projector_scan.jl`
- Modify: `src/TopoNoise.jl`
- Modify: `test/toric3d/runtests.jl`

**Interfaces:**
- Consumes: replay API and decoder settings.
- Produces: `Toric3DVertexScanGrid(dimensions,error_rates)`, `Toric3DVertexProtocolEstimate`, `Toric3DVertexPairEstimate`, `Toric3DVertexQuadMetricEstimate`, `Toric3DVertexFeedbackPoint`, `Toric3DVertexFeedbackScan`, and `compare_toric3d_vertex_projector_feedback(rng; channel_grids,protocols=TORIC3D_VERTEX_PROTOCOLS,shots,seed,settings,batch_size,checkpoint_dir,resume)`. `channel_grids` maps each requested `:x_only`/`:z_only` channel to its own grid. Accept `protocols` as an explicit reproducibility input, validate that it equals the exact canonical four-entry tuple, record it in checkpoints, and reject reordered, missing, or extra entries because sixteen-mask statistics require all four arms.

**Testsets:** `"Sixteen joint patterns determine all six H+5 protocol contrasts"`, `"H+5 scans are four-way paired, channel complete, and seeded"`, `"Interrupted H+5 scans resume exact point substreams"`, `"H+5 resume rejects every configuration mismatch before sampling"`, `"H+5 invalid recovery context is atomically saved before abort"`, and `"H+5 user interruption is not an invalid recovery"`.

- [ ] **Step 1: Write the hand-calculated 16-pattern statistics test**

  From one literal 16-count vector, assert all four marginals, all six paired discordance tables, Wilson intervals, and the three highlighted feedback-versus-`none` contrasts.

- [ ] **Step 2: Write paired, seeded, and channel-complete scan tests**

  Require one record per shot replayed through all four arms, independent `(T,Lu,Lv)` triples, repeatable point substreams, full logical-class counts, recorded timings/settings/versions, and zero-noise success. Give X and Z different dimension/rate grids and prove the point set is their disjoint per-channel union rather than a shared Cartesian grid; within each channel, every selected rate still has all four paired arms and sixteen joint masks. Accept the exact canonical `protocols` tuple and reject a subset or reordering before the first RNG draw.

- [ ] **Step 3: Write exact resume tests**

  Interrupt after a partial batch, resume with a different caller RNG, and compare point seeds, sixteen-way counts, boundary branches, and the next generated record with an uninterrupted scan. Discard complete records, including `b`, when reconstructing the point stream.

- [ ] **Step 4: Write pre-draw mismatch and failure-record tests**

  Change dimensions, rates, channels, the fixed protocol tuple/version, shots, batch size, seed, decoder settings, or pinned backend metadata one at a time; each resume must fail before the draw hook fires. An invalid recovery must atomically save syndromes, paths, route, supports, `b`, both physical boundary corrections, protocol/checkpoint/shot, and point seed, then abort. `InterruptException` must never increment invalid recovery counts.

- [ ] **Step 5: Run the test and confirm failure**

  Run: `JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. test/toric3d/vertex_projector_scan.jl`.

  Expected: FAIL on missing scan types.

- [ ] **Step 6: Implement statistics and atomic point state**

  Keep this format separate from the H+6 three-arm/eight-pattern types. Mark a point complete only after its final batch and write through a temporary file plus atomic rename.

- [ ] **Step 7: Run scan, replay, and legacy scan suites**

  Run: `JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. -e 'using Test, Random, TopoNoise; include("test/toric3d/vertex_projector_replay.jl"); include("test/toric3d/vertex_projector_scan.jl"); include("test/toric3d/feedback_strength_scan.jl")'`.

  Expected: PASS, including deterministic resume and invalid-recovery aborts.

- [ ] **Step 8: Commit**

  ```bash
  git add src/TopoNoise.jl src/toric3d/vertex_projector_scan.jl test/toric3d/runtests.jl test/toric3d/vertex_projector_scan.jl
  git commit -m "feat: scan four paired H+5 protocols"
  ```

### Task 10: Add CSV/metadata/figures and the runnable pilot/scout/production CLI

**Files:**
- Create: `src/toric3d/vertex_projector_output.jl`
- Create: `examples/compare_toric3d_vertex_projector_feedback.jl`
- Create: `test/toric3d/vertex_projector_output.jl`
- Create: `test/toric3d/vertex_projector_example.jl`
- Modify: `src/TopoNoise.jl`
- Modify: `test/toric3d/runtests.jl`

**Interfaces:**
- Consumes: `Toric3DVertexFeedbackScan`.
- Produces: `write_toric3d_vertex_projector_feedback_csv`, `plot_toric3d_vertex_projector_feedback`, `save_toric3d_vertex_projector_feedback`, `Toric3DVertexProductionSelection`, `load_toric3d_vertex_production_selection(path)`, and a CLI with `--stage pilot|scout|extend-scout|production`, per-channel dimensions/rates, `--protocols`, shots/seed/settings/checkpoint/resume controls. `--protocols` defaults to and must equal `none,one_step,window_complete,full_prefix`; any other value fails before sampling.

**Testsets:** `"H+5 artifacts preserve four-way paired provenance"`, `"H+5 figures distinguish logical diagnostics from state failure"`, `"H+5 selection schema supports separate channel grids"`, and `"H+5 direct-projector scout and production CLI"`.

- [ ] **Step 1: Write failing artifact-schema tests**

  Pin raw and aggregate CSV columns for sixteen joint patterns, four marginals, six contrasts, logical components, either-pair unions, and full Pauli classes. Require `p`, every BP+OSD `p_eff`, canonical matching signatures, decoder versions/settings, source hashes, branch policy, timings, completion status, and zero invalid recoveries.

- [ ] **Step 2: Write failing plot tests**

  Require four protocol curves per channel, Wilson intervals, paired-difference panels, PDF/SVG/PNG output, public `T,Lu,Lv` labels, and no logical-Z-as-`|00_L>`-infidelity caption.

- [ ] **Step 3: Write failing CLI tests**

  Pin scout defaults `L=3,5,7`, 500 shots, seed 1234, both channels, the exact protocol string `none,one_step,window_complete,full_prefix`, and probabilities `0,1e-4,3e-4,1e-3,3e-3,.01,.02,.03,.05,.075,.10`. Reject any subset, extra, or reordering before sampling. Define a version-1 `selection.json` fixture before Task 11: top-level source hashes and selection-rule version; separate `x_only` and `z_only` objects containing dimension triples, sorted rates, per-protocol frozen fit windows, candidate centers/flags, saturation bounds, measured timing inputs, projected seconds, and the 1000 production shots. Production mode must validate this complete schema and refuse mismatched hashes/settings before sampling. Require a timed pilot before any scientific draw, incremental point saves, and exact resume.

  Pin these exact schema keys: `format_version=1`, `rule_version="h5-scout-selection-v1"`, `status` (`"ready"` or `"needs_extended_scout"`), `extension_round`, `scout_inputs=[{csv_sha256,metadata_sha256,path_ids,stage,extension_round,parent_selection_sha256}]`, `seed`, `decoder_settings`, `production_shots=1000`, `extension_grids`, and `channels.{x_only,z_only}={dimensions_T_Lu_Lv,rates,fit_windows,candidate_centers,extension_states,flags,saturation_bounds,timing_model,projected_seconds}`. Each `fit_windows`, `candidate_centers`, and `extension_states` object contains exactly `none`, `one_step`, `window_complete`, and `full_prefix`; the first two values are ordered parallel arrays because noisy scout data can retain multiple candidates, while each extension state is exactly `"resolved_candidates"`, `"requesting"`, or `"exhausted_diagnostic"`. Production refuses any status except `"ready"`; `extend-scout` accepts only `"needs_extended_scout"` and runs the declared extension grids into a new shard directory.

- [ ] **Step 4: Run output/example tests and confirm failure**

  Run: `JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. -e 'using Test, Random, TopoNoise; include("test/toric3d/vertex_projector_output.jl"); include("test/toric3d/vertex_projector_example.jl")'`.

  Expected: FAIL.

- [ ] **Step 5: Implement artifact writers and figures**

  Metadata must say “direct-projector H+5,” eleven fault locations per cell, ideal boundary operations, circuit-correlation limitations, and that the historical H+6 code is “plaquette-transfer.”

- [ ] **Step 6: Implement and test the staged CLI**

  Implement the versioned selection parser against the Task 10 fixture. Scout and extension runs each own an exact-resume checkpoint and never append incompatible grids to an existing checkpoint. Production mode must require a `status="ready"` selection emitted by Task 11, instantiate separate X/Z `channel_grids`, rerun all four arms at the unioned per-channel rates, and never derive a grid from production outcomes. Warmup/timing records never contribute scientific counts.

- [ ] **Step 7: Run output/example and legacy artifact suites**

  Run: `JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. -e 'using Test, Random, TopoNoise; include("test/toric3d/vertex_projector_output.jl"); include("test/toric3d/vertex_projector_example.jl"); include("test/toric3d/feedback_strength_output.jl"); include("test/toric3d/feedback_strength_example.jl")'`.

  Expected: PASS without changing old schemas or hashes.

- [ ] **Step 8: Commit**

  ```bash
  git add src/TopoNoise.jl src/toric3d/vertex_projector_output.jl examples/compare_toric3d_vertex_projector_feedback.jl test/toric3d/runtests.jl test/toric3d/vertex_projector_output.jl test/toric3d/vertex_projector_example.jl
  git commit -m "feat: deliver H+5 scan artifacts and CLI"
  ```

### Task 11: Add frozen-grid selection and guarded threshold analysis

**Files:**
- Create: `examples/analyze_toric3d_vertex_projector.py`
- Create: `test/toric3d/test_vertex_projector_fit.py`

**Interfaces:**
- Consumes: one or more paired H+5 scout CSV/metadata shards and the version-1 `Toric3DVertexProductionSelection` schema defined and tested in Task 10. `select` accepts repeatable `--input-csv` and `--input-metadata` arguments with equal lengths plus `--previous-selection` on every reselection after round zero.
- Produces: `selection.json`, guarded per-channel/protocol fit JSON, sensitivity CSV, bootstrap CSV, report Markdown, and fit figures.

- [ ] **Step 1: Write failing four-arm parser/statistics tests**

  Add `test_reads_four_arm_joint_counts_and_six_contrasts` and verify independent reconstruction of every marginal and pair table from raw sixteen-way counts.

- [ ] **Step 2: Write failing frozen-selection tests**

  Add `test_scout_selection_is_frozen_before_production`, `test_scout_shards_merge_exact_points`, `test_multiple_crossings_fund_every_candidate_window`, `test_empty_saturation_extends_then_terminates_at_bounds`, `test_mixed_arm_extension_status_uses_active_union`, `test_extension_rounds_require_cumulative_parent_chain`, and `test_timing_rule_selects_expected_channel_sizes`.

  Process X and Z independently. A rate is unsaturated for one adjacent-size pair and protocol only when both primary rates lie in `[0.02,0.98]`. Linearly interpolate every positive sign change. Sort candidates by `p`, greedily form complete-link clusters whose largest/smallest ratio is at most `1.5`, and rank clusters by: more distinct adjacent-size pairs, then more candidates, then smaller log-span, then smaller geometric-mean center. Retain every cluster as a candidate, order them by that rank, and define each center as `pc=float(format(exp(mean(log(p_i))),".8g"))`.

  Around every retained center independently, use `float(format(x,".8g"))` on `pc .* [0.8,0.9,1.0,1.1,1.2]`, sort/deduplicate, and freeze one fit window to that set's first/last value. `candidate_centers[protocol]` and `fit_windows[protocol]` are parallel arrays. Make the channel production rates the union of all five-rate sets across every candidate and all four arms, so every arm is rerun at every channel rate. During analysis, fit every frozen window; claim a unique threshold only if exactly one candidate passes every guard and its cluster contains every adjacent size pair, otherwise report ordered diagnostic candidates.

  If there is no sign change, use the interior unsaturated rate minimizing the absolute difference for the largest adjacent-size pair and mark `no_scout_crossing`. If there is no unsaturated rate at all, compute the median over measured rates of the mean failure of that largest pair: below `0.5` requests an upper extension, above `0.5` requests a lower extension, and exactly `0.5` requests both. A lower extension proposes `max(1e-8,pmin/10)` and `max(1e-8,pmin/3)`; an upper extension proposes `min(0.49,1.5*pmax)` and `min(0.49,2*pmax)`, all rounded with `.8g` and with already measured/duplicate rates removed. First compute every channel/protocol `extension_state`. Union only new rates from arms in `"requesting"` into each channel's `extension_grids`; set global `status="needs_extended_scout"` if and only if this union is nonempty and `extension_round < 3`.

  If round three is reached or clipping produces no new rate for one arm, set only that arm's state to `"exhausted_diagnostic"`, with `extension_exhausted=true` and `threshold_ineligible=true`; other arms may remain `"requesting"`. Freeze that arm's diagnostic window from the five distinct measured positive rates closest in absolute log-distance to the boundary center (or every available rate if fewer than five); the fitter must never turn that window into a threshold claim. After all arms are classified, set global `status="ready"` exactly when no arm remains `"requesting"`. These `1e-8`/`0.49` bounds and the three-round limit make extension terminal. Add a literal mixed-arm fixture where one arm exhausts at `0.49` while another requests a lower extension and require the global status to remain `"needs_extended_scout"`.

  The initial scout has `extension_round=0`. Every `extend-scout` shard has round `previous.extension_round+1`, records the SHA-256 of the exact parent selection, and records the parent selection's full ordered `scout_inputs`. Every reselection supplies `--previous-selection` and all original plus earlier extension CSV/metadata pairs; reject a missing prior pair, a gap, duplicate round, wrong parent hash, or a new shard whose inputs do not equal the previous inputs plus that shard. Increment the output round from the validated new shard, never from a caller counter. This cumulative ancestry is the source of truth for the three-round limit.

  For timing, at each channel/size take the median across complete scout probabilities of `(sampling_seconds+replay_seconds)/shots` and separately the median `decoder_build_seconds`. Fit two ordinary least-squares models in natural-log coordinates, `log(seconds)=a+b*log(L)`, on `L=3,5,7`. For candidate prefixes `[3,5,7,9]`, `[3,5,7,9,11]`, and `[3,5,7,9,11,13]`, project `N_rates * sum(exp(a_build+b_build*log(L)) + 1000*exp(a_shot+b_shot*log(L)))`. Select the longest prefix not exceeding 43,200 seconds for that channel. If the four-size prefix exceeds the budget, retain it and set `runtime_budget_exceeded=true`. A literal synthetic fixture must select `[3,5,7,9,11]` for X and `[3,5,7,9]` for Z.

  Merge repeatable scout shards only when seed, decoder settings, protocol tuple, source hashes, shot count, schemas, and extension ancestry match. Reject conflicting duplicate points; accept identical duplicates once. Write separate X/Z dimensions, unioned rates, ordered retained candidate arrays, one frozen window per retained candidate, per-arm extension states, saturation decisions, extension grids, timing coefficients/projections, rule version, settings, and every scout hash to the Task 10 schema. Prove that adding or changing production rows cannot alter selection.

- [ ] **Step 3: Write failing guard and bootstrap tests**

  Add `test_seeded_bootstrap_is_reproducible_per_channel_and_protocol`, `test_all_threshold_guards_withhold_threshold_claims`, and `test_reports_hashes_sensitivities_and_diagnostic_candidates`. Preserve adjacent-size crossing, interior monotone fit, likelihood goodness-of-fit, and >=80% bootstrap acceptance guards.

- [ ] **Step 4: Run Python tests and confirm failure**

  Run:
  `PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/toponoise-mpl .CondaPkg/.pixi/envs/default/bin/python -m unittest test/toric3d/test_vertex_projector_fit.py -v`

  Expected: FAIL because the analysis script is missing.

- [ ] **Step 5: Implement `select`, `fit`, and `report` subcommands**

  Merge and validate repeatable scout inputs and their cumulative parent-selection chain, then use only those rows for selection. Emit either a ready production selection or deterministic extension grids. Fit the predeclared either-pair primary observable separately for every channel/protocol/candidate inside that candidate's parallel frozen window. Report ordered “diagnostic candidates” unless exactly one candidate passes all claim guards; keep statistical intervals separate from model/window sensitivity.

- [ ] **Step 6: Run new and existing Python analysis suites**

  Run:

  `PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/toponoise-mpl .CondaPkg/.pixi/envs/default/bin/python -m unittest discover -s test/toric3d -p 'test_*.py' -v`

  Expected: PASS.

- [ ] **Step 7: Commit**

  ```bash
  git add examples/analyze_toric3d_vertex_projector.py test/toric3d/test_vertex_projector_fit.py
  git commit -m "feat: analyze guarded H+5 thresholds"
  ```

### Task 12: Update the note, experiment documentation, and full regressions

**Files:**
- Create: `docs/superpowers/toric3d-h5-vertex-projector.md`
- Modify: `document/note.tex`
- Modify: `README.md`

**Interfaces:**
- Consumes: implemented semantics and later generated artifacts.
- Produces: corrected 3D discussion, exact reproduction commands, and an experiment report template populated after Task 13.

- [ ] **Step 1: Rewrite the diagonal-feedback discussion around the final H+5 implementation**

  State the direct vertex projector, diagonal coordinates, lower/upper boundary instrument, four protocols, X/Z one-step differences, eleven post-gate locations, two logical pairs, state-independent diagnostics, and BP+OSD/MWPM limitations. Retain the corrected Cartesian normalization, boundary checks, and logical-Z lower bound with its assumptions.

- [ ] **Step 2: Relabel historical H+6 text and artifacts**

  Call the seven-gate circuit “plaquette-transfer H+6,” preserve its numerical results as historical, and explicitly state that H+5 is paper-motivated rather than a numerical reproduction of a supplied 3D circuit. Warn that H+5/H+6 threshold values are not a controlled circuit-only comparison because the seed preparation, boundary instrument, and historical channel-specific seed noise also differ.

- [ ] **Step 3: Add runnable commands and artifact interpretation to README/report**

  Explain pilot, scout, frozen selection, production resume, guarded fits, and why feedback-vs-none comparisons do not isolate correction range alone.

- [ ] **Step 4: Run all Julia and Python regressions**

  Run:

  ```bash
  JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. test/runtests.jl
  PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/toponoise-mpl .CondaPkg/.pixi/envs/default/bin/python -m unittest discover -s test/toric3d -p 'test_*.py' -v
  git diff --check
  ```

  Expected: all tests PASS and no whitespace errors.

- [ ] **Step 5: Build the note PDF**

  Run `latexmk -pdf -interaction=nonstopmode -halt-on-error note.tex` from `document/` and inspect the 3D section pages for clipped equations, captions, and missing references.

- [ ] **Step 6: Commit**

  ```bash
  git add README.md document/note.tex docs/superpowers/toric3d-h5-vertex-projector.md
  git commit -m "docs: describe diagonal H+5 feedback experiment"
  ```

### Task 13: Run the timed pilot, scout, frozen production scan, and final analysis

**User amendment, 2026-10-14:** After measured large-size runtimes exceeded the small-scout projection, the user requested smaller lattices first. Finish a separate preliminary run at cubic L=3,5,7,9; defer L11,13. Preserve the original selection and interrupted full checkpoint unchanged. Use `results/toric3d-h5-diagonal-feedback-2026-10-13/preliminary-l3-9/selection.json` (SHA-256 `dca614dc243ec8d918cb91ac27c74425c7af0fbba6f2381c88d33d3ebf025412`) and its `production/` and `analysis/` directories for Steps 5–9. Rates, candidate windows, four protocols, 1000 records per point, seed1234, decoder settings and source hashes remain unchanged. The amendment is motivated only by runtime. Reuse the 80 completed X points only after checking exact counts, digests and fresh-run point seeds; the 180 Z points use the fresh reduced-grid seed stream. `prepare_checkpoint.jl` and `checkpoint-reuse.toml` document this derivation. Retain all original fit guards and label the resulting measurements preliminary. The original six-size run is deliberately deferred, not failed or complete, and must not resume automatically. Finish the audits, report, note, final regressions and local commit for this revised scope, then disable the follow-up.

**Files:**
- Create: `results/toric3d-h5-diagonal-feedback-2026-10-13/pilot/*`
- Create: `results/toric3d-h5-diagonal-feedback-2026-10-13/scout/*`
- Create: `results/toric3d-h5-diagonal-feedback-2026-10-13/production/*`
- Modify: `docs/superpowers/toric3d-h5-vertex-projector.md`
- Modify: `document/note.tex`

**Interfaces:**
- Consumes: staged CLI and analysis script.
- Produces: raw/aggregate CSV, metadata, selection, checkpoints, failure records if any, PDF/SVG/PNG curves, guarded threshold reports, and final documented conclusions.

- [ ] **Step 1: Run a short timed pilot**

  Run both channels at `(T,Lu,Lv)=(3,3,3)`, `p=0.01`, 20 paired records, all four protocols, seed 1234, and save timing separately from science counts:

  `JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. examples/compare_toric3d_vertex_projector_feedback.jl --stage pilot --dimensions 3x3x3 --rates 0.01 --channels x_only,z_only --shots 20 --seed 1234 --output results/toric3d-h5-diagonal-feedback-2026-10-13/pilot`

- [ ] **Step 2: Review pilot validity and runtime**

  Require zero invalid recoveries, zero-noise success from the warmup, all four replay arms, and a recorded runtime projection before launching the scout.

- [ ] **Step 3: Run or resume the exact scout**

  Use cubic `L=3,5,7`, 500 paired shots, both channels, all four protocols, seed 1234, and the eleven-point probability grid from Task 10. Save every completed point atomically:

  `JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. examples/compare_toric3d_vertex_projector_feedback.jl --stage scout --dimensions 3x3x3,5x5x5,7x7x7 --rates 0,0.0001,0.0003,0.001,0.003,0.01,0.02,0.03,0.05,0.075,0.10 --channels x_only,z_only --shots 500 --seed 1234 --output results/toric3d-h5-diagonal-feedback-2026-10-13/scout`

  If interrupted, repeat this exact command with `--resume`.

- [ ] **Step 4: Freeze separate X and Z production selections**

  Run the `select` subcommand on scout artifacts only:

  `.CondaPkg/.pixi/envs/default/bin/python examples/analyze_toric3d_vertex_projector.py select --input-csv results/toric3d-h5-diagonal-feedback-2026-10-13/scout/toric3d-h5-vertex-projector-feedback-raw.csv --input-metadata results/toric3d-h5-diagonal-feedback-2026-10-13/scout/toric3d-h5-vertex-projector-feedback-metadata.toml --output results/toric3d-h5-diagonal-feedback-2026-10-13/selection.json`

  Inspect `selection.json`, its input hashes, separate X/Z four-to-six-size grids, frozen per-protocol fit windows, saturation and runtime decisions, five-or-more rates around every candidate crossing, and its status. If it says `needs_extended_scout`, run exactly its declared points into a separate shard:

  `JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. examples/compare_toric3d_vertex_projector_feedback.jl --stage extend-scout --selection results/toric3d-h5-diagonal-feedback-2026-10-13/selection.json --seed 1234 --output results/toric3d-h5-diagonal-feedback-2026-10-13/scout-extension-1`

  Then rerun:

  `.CondaPkg/.pixi/envs/default/bin/python examples/analyze_toric3d_vertex_projector.py select --previous-selection results/toric3d-h5-diagonal-feedback-2026-10-13/selection.json --input-csv results/toric3d-h5-diagonal-feedback-2026-10-13/scout/toric3d-h5-vertex-projector-feedback-raw.csv --input-metadata results/toric3d-h5-diagonal-feedback-2026-10-13/scout/toric3d-h5-vertex-projector-feedback-metadata.toml --input-csv results/toric3d-h5-diagonal-feedback-2026-10-13/scout-extension-1/toric3d-h5-vertex-projector-feedback-raw.csv --input-metadata results/toric3d-h5-diagonal-feedback-2026-10-13/scout-extension-1/toric3d-h5-vertex-projector-feedback-metadata.toml --output results/toric3d-h5-diagonal-feedback-2026-10-13/selection.json`

  Repeat with a new shard number only if the deterministic rule requests another extension. On every reselection, pass the current selection as `--previous-selection` and include the initial scout plus every prior extension CSV/metadata pair in round order; dropping or replacing any earlier shard is an error. Freeze the first `status="ready"` selection hash before production; never resume an old checkpoint with changed grids.

- [ ] **Step 5: Run or resume the production scan**

  Use the frozen selection, 1000 paired records per point, seed 1234, both channels, all four protocols, and incremental checkpoints:

  `JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" julia --project=. examples/compare_toric3d_vertex_projector_feedback.jl --stage production --selection results/toric3d-h5-diagonal-feedback-2026-10-13/selection.json --shots 1000 --seed 1234 --output results/toric3d-h5-diagonal-feedback-2026-10-13/production`

  If interrupted, repeat this exact command with `--resume`. Stop any affected point immediately if recovery validation creates a failure record; diagnose and fix before resuming.

- [ ] **Step 6: Generate final artifacts and guarded fits**

  Run:

  `.CondaPkg/.pixi/envs/default/bin/python examples/analyze_toric3d_vertex_projector.py fit --selection results/toric3d-h5-diagonal-feedback-2026-10-13/selection.json --csv results/toric3d-h5-diagonal-feedback-2026-10-13/production/toric3d-h5-vertex-projector-feedback-raw.csv --metadata results/toric3d-h5-diagonal-feedback-2026-10-13/production/toric3d-h5-vertex-projector-feedback-metadata.toml --bootstrap 1000 --seed 1234 --output results/toric3d-h5-diagonal-feedback-2026-10-13/analysis`

  Save CSV plus PDF/SVG/PNG curves, reconstruct all sixteen-way counts and six contrasts independently, and report threshold values only for arms that pass every guard.

- [ ] **Step 7: Inspect every final figure and audit provenance**

  Check readable legends, coincident-curve visibility, error bars, correct channel/state labels, source/data hashes, backend versions, `p/p_eff`, zero invalid recoveries, complete point counts, and exclusion of timing/warmup records.

- [ ] **Step 8: Populate the report and note with measured results**

  Record thresholds or guarded diagnostic candidates, uncertainty and sensitivity separately, fixed-size paired feedback effects, total shots/runtime, exact artifact paths, and scientific limitations. Never reuse the historical H+6 threshold values as H+5 results.

- [ ] **Step 9: Run final verification and commit the reproducible delivery**

  Re-run full Julia/Python tests, rebuild the note PDF, run `git diff --check`, and verify all reported numbers directly from saved CSV/JSON.

  ```bash
  git add results/toric3d-h5-diagonal-feedback-2026-10-13 docs/superpowers/toric3d-h5-vertex-projector.md document/note.tex
  git commit -m "data: report diagonal H+5 feedback thresholds"
  ```
