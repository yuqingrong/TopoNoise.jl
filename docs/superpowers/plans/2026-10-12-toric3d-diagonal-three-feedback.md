# Diagonal 3D Toric-Code Three-Feedback Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement and run a paired diagonal 3D toric-code comparison of no intermediate feedback, constant-width slice-local feedback, and full-prefix feedback under separate pure X and pure Z gate-output noise.

**Architecture:** Add a separate feedback-strength subsystem around the existing encoded diagonal geometry and final decoders. It defines ordered prefix mappings, a syndrome-complete sliding Z window, three-way X/Z replay on immutable per-shot faults, joint statistics, resumable scans, artifacts, and guarded threshold analysis while preserving every existing interface and result.

**Tech Stack:** Julia 1.12, binary Pauli frames, PythonCall, ldpc BP+OSD 2.4.1, PyMatching 2.4.0, CairoMakie, Python 3.12, NumPy/SciPy/Matplotlib, TOML/CSV, LaTeX.

**Spec:** `docs/superpowers/specs/2026-10-12-toric3d-diagonal-three-feedback-design.md`

## Global Constraints

- Preserve all current 2D, rough-slab 3D, fixed-boundary slice, encoded-Z, and encoded-X public APIs and artifact formats.
- Use `Toric3DEncodedSlice(T,Ly,Lz)` with `n=3(T+1)LyLz`, `k=2`, `dX=(T+1)min(Ly,Lz)`, and `dZ=2min(Ly,Lz)`.
- Keep the seven-gate transfer and independent rate-`p` faults after all thirteen participating H/CNOT outputs per cell.
- Keep initialization, idle qubits, measurements, processing, and Pauli corrections ideal.
- Prepare logical `|00>` with no noisy seed gates for X; prepare logical `|++>` with `3LyLz` noisy seed-H outputs for Z.
- For every shot, replay one immutable fault record through `:none`, `:slice_local`, and `:full_prefix`.
- X slice-local feedback is the existing one-step compatibility-check target feed-forward; `:full_path` remains only a routing sensitivity.
- Z slice-local feedback measures the three-sheet constant-width band and applies a complete recovery on the newest two data slices.
- X full-prefix recovery uses minimum-sum BP, 50 iterations, scaling 1, parallel schedule, one thread, and `OSD_CS` order 2 with uniform declared `p` clamped to `[1e-12,1-1e-12]`.
- Z seed/window/prefix/final recovery uses deterministic unit-data-weight MWPM.
- Reject and persist context for every full recovery that does not reproduce its measured syndrome; never score it silently.
- Decoders receive syndromes and declared model parameters only, never sampled fault bits.
- Save batch progress atomically and require exact configuration equality on resume.
- Withhold `p_c` when the crossing or any numerical/statistical guard fails.
- Use the managed decoder interpreter through `JULIA_CONDAPKG_BACKEND=Null` and `JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python"` on this machine.
- Preserve the pre-existing modified generated files `document/note.fdb_latexmk`, `document/note.fls`, and `document/note.synctex.gz`; never stage them.

## Review Focus

- `T=1`: X `:full_prefix` must equal `:none`; Z `:full_prefix` must equal `:slice_local` because both perform only seed feedback before the common final decode. Task 4 pins both cases.
- `p=0` with a manually supplied nonzero X syndrome must still use the clamped BP+OSD prior, while zero sampled faults must give zero failures in every protocol. Tasks 2, 3, and 6 pin this.
- `Ly=2` or `Lz=2` must preserve distinct periodic edges, two-endpoint Z-window columns, and both logical parities. Tasks 1 and 2 pin these geometries.
- Resume with different dimensions, rates, channels, shots, batch size, routing, seed provenance, or BP+OSD settings must fail before sampling. Task 6 pins every field.
- An invalid intermediate or final recovery must atomically save its protocol/checkpoint/shot context before aborting; a completed artifact must report zero invalid recoveries. Tasks 2, 6, and 7 pin this behavior.

---

## File structure

- `src/toric3d/feedback_strength_geometry.jl`: ordered encoded-prefix maps and Z seed/window measurement/correction geometry.
- `src/toric3d/feedback_strength_decoders.jl`: validated seed/window MWPM and decoder-family construction.
- `src/toric3d/feedback_strength_x.jl`: three-protocol X replay and full-prefix X history.
- `src/toric3d/feedback_strength_z.jl`: Z fault records, scalar three-protocol replay, and checkpoint histories.
- `src/toric3d/feedback_strength_z_batch.jl`: batched Z sampling/replay used by scans.
- `src/toric3d/feedback_strength_scan.jl`: joint statistics, point substreams, incremental checkpointing, and public scan API.
- `src/toric3d/feedback_strength_output.jl`: CSV/TOML/failure records and CairoMakie figures.
- `examples/compare_toric3d_diagonal_feedback.jl`: warmup, timing pilot, scout/production CLI.
- `examples/fit_toric3d_diagonal_feedback.py`: guarded per-channel/per-protocol scaling analysis.
- `test/toric3d/feedback_strength_*.jl` and `test/toric3d/test_feedback_strength_fit.py`: focused red/green coverage.

### Task 1: Add ordered prefix and Z-window geometry

**Files:**
- Create: `src/toric3d/feedback_strength_geometry.jl`
- Modify: `src/TopoNoise.jl:66-77,103-108`
- Create: `test/toric3d/feedback_strength_geometry.jl`
- Modify: `test/toric3d/runtests.jl:16-23`

**Interfaces:**
- Consumes: `Toric3DEncodedSlice`, `star_checks`, `toric3d_encoded_seed_checks`, `toric3d_slice_qubits`, `data_qubit_count`, `data_qubit_coordinate`.
- Produces: `Toric3DEncodedPrefix`, `toric3d_encoded_prefix(code::Toric3DEncodedSlice, t::Integer)::Toric3DEncodedPrefix`, `Toric3DZFeedbackWindow`, `toric3d_z_feedback_window(code::Toric3DEncodedSlice, r::Integer)::Toric3DZFeedbackWindow`, and `toric3d_z_window_syndrome(window::Toric3DZFeedbackWindow, frame::AbstractVector{Bool})::BitVector`.

- [ ] **Step 1: Write the failing prefix-map tests**

  In `feedback_strength_geometry.jl`, assert for `(T,Ly,Lz)=(3,2,3)` and every `t=1:T` that `length(prefix.qubits)==data_qubit_count(prefix.code)`, local coordinate `q` equals the full-code coordinate at `prefix.qubits[q]`, mappings contain no duplicate/future qubit, and invalid `t=0,T+1` throws `ArgumentError`.

- [ ] **Step 2: Write the failing Z-window tests**

  Assert with `M=LyLz`: `r=0` has a `2M x 3M` matrix and seed-only support; `r=1` has `3M x 6M`; every `r>=1` has correction qubits exactly on slices `r-1:r`, every matrix column has weight two, and measured check support reaches `0:1` at `r=1` and `r-2:r` thereafter. Include `(T,Ly,Lz)=(2,2,2),(2,2,3),(3,3,2)` and verify distinct length-two periodic edges.

- [ ] **Step 3: Run the focused test and observe RED**

  Run the managed Julia environment with `julia --project=. -e 'using Test, TopoNoise; include("test/toric3d/feedback_strength_geometry.jl")'`.

  Expected: FAIL because `Toric3DEncodedPrefix` and `Toric3DZFeedbackWindow` are undefined.

- [ ] **Step 4: Implement the geometry interfaces**

  `Toric3DEncodedPrefix` stores the prefix code plus an ordered local-column-to-full-qubit vector. `Toric3DZFeedbackWindow` stores checkpoint `r`, full-code check supports, ordered correction qubits, and the restricted binary recovery matrix. Use lower-white/A0/B1 at `r=1` and `A_(r-2)/A_(r-1)/B_r` thereafter.

- [ ] **Step 5: Export/include and verify GREEN**

  Include the file after `encoded_slice_geometry.jl`, export the five interfaces, include the test after `encoded_slice_zero.jl`, and rerun Task 1's test plus `test/toric3d/encoded_slice_geometry.jl`.

  Expected: all assertions pass.

- [ ] **Step 6: Commit the geometry**

  ```bash
  git add src/toric3d/feedback_strength_geometry.jl src/TopoNoise.jl \
    test/toric3d/feedback_strength_geometry.jl test/toric3d/runtests.jl
  git commit -m "feat: define diagonal feedback windows"
  ```

### Task 2: Add validated Z graph decoders and single-shot final outcome

**Files:**
- Create: `src/toric3d/feedback_strength_decoders.jl`
- Modify: `src/toric3d/z_threshold.jl:76-141`
- Modify: `src/TopoNoise.jl`
- Create: `test/toric3d/feedback_strength_decoders.jl`
- Modify: `test/toric3d/runtests.jl`

**Interfaces:**
- Consumes: Task 1 windows, `_toric3d_scipy_sparse`, `_load_pymatching`, `_toric3d_py_bits`, `Toric3DZThresholdDecoder`.
- Produces: `Toric3DZFeedbackDecoder`, `Toric3DZFeedbackDecoderFamily`, `build_toric3d_z_feedback_decoders(code::Toric3DEncodedSlice)::Toric3DZFeedbackDecoderFamily`, `decode_toric3d_z_feedback(decoder::Toric3DZFeedbackDecoder, syndrome::AbstractVector{Bool})::BitVector`, `decode_toric3d_z_feedback_batch(decoder::Toric3DZFeedbackDecoder, syndromes::AbstractMatrix{Bool})::BitMatrix`, and `toric3d_z_threshold_outcome(decoder::Toric3DZThresholdDecoder, frame::AbstractVector{Bool})::NamedTuple`.

- [ ] **Step 1: Write failing decoder-construction and recovery tests**

  Require one seed decoder, one reusable transfer-window decoder, and prefix decoders for `t=1:T`. Assert the `r=1` and `r=2` restricted matrices are identical under their ordered local columns. For every column of the seed and transfer matrices, decode its syndrome and require exact syndrome reproduction. Assert wrong-length and odd syndromes throw rather than score.

- [ ] **Step 2: Write failing final-outcome tests**

  Inject zero, one stabilizer, and each of the two encoded logical-Z cycles. Require the per-shot outcome to return measured syndrome, validated recovery, residual, both logical flags, and their union. Corrupt a recovery through the private validator and require an `ErrorException`.

- [ ] **Step 3: Run the decoder test and observe RED**

  Run the managed Julia environment with `include("test/toric3d/feedback_strength_decoders.jl")`.

  Expected: FAIL because the feedback decoder family and single-shot outcome are undefined.

- [ ] **Step 4: Implement the decoder family**

  Build unit-weight PyMatching graphs from the Task 1 seed/window matrices, validate every scalar/batch output by multiplying the recovery matrix, and build growing-prefix `Toric3DZThresholdDecoder` values once per `t`. Keep deterministic matrix ordering and tie handling.

- [ ] **Step 5: Implement `toric3d_z_threshold_outcome`**

  Measure the full star syndrome, decode once, apply recovery, and pair the residual with both encoded logical-X membranes. Refactor `toric3d_z_threshold_counts` only enough to share validation/scoring while preserving its output and batching.

- [ ] **Step 6: Export/include and verify GREEN**

  Include the decoder file after `z_threshold.jl`, export its public types/functions and the single-shot outcome, then run the new test and `test/toric3d/z_threshold.jl`.

  Expected: all tests pass, including transverse length two and odd-syndrome rejection.

- [ ] **Step 7: Commit the decoders**

  ```bash
  git add src/toric3d/feedback_strength_decoders.jl src/toric3d/z_threshold.jl \
    src/TopoNoise.jl test/toric3d/feedback_strength_decoders.jl test/toric3d/runtests.jl
  git commit -m "feat: decode diagonal Z feedback windows"
  ```

### Task 3: Add full-prefix X replay and the X protocol triplet

**Files:**
- Create: `src/toric3d/feedback_strength_x.jl`
- Modify: `src/TopoNoise.jl`
- Create: `test/toric3d/feedback_strength_x.jl`
- Modify: `test/toric3d/runtests.jl`

**Interfaces:**
- Consumes: Task 1 prefixes, `Toric3DSliceDecoder`, `Toric3DXThresholdDecoder`, `replay_toric3d_x_feedback`, `propagate_toric3d_segment!`, `decode_toric3d_x`, `toric3d_x_threshold_outcome`.
- Produces: `Toric3DXFeedbackDecoderFamily`, `build_toric3d_x_feedback_decoders(code::Toric3DEncodedSlice; p::Real, settings::Toric3DDecoderSettings=Toric3DDecoderSettings())::Toric3DXFeedbackDecoderFamily`, `replay_toric3d_x_strength(encoder::Toric3DSliceEncoder, faults::Toric3DFaultRecord, decoders::Toric3DXFeedbackDecoderFamily; protocol::Symbol, routing::Symbol=:one_step, record_history::Bool=false)::NamedTuple`, and `toric3d_x_strength_triplet(encoder::Toric3DSliceEncoder, faults::Toric3DFaultRecord, decoders::Toric3DXFeedbackDecoderFamily; routing::Symbol=:one_step, record_history::Bool=false)::NamedTuple{(:none,:slice_local,:full_prefix)}`.

- [ ] **Step 1: Write the failing compatibility tests for `:none` and `:slice_local`**

  For clean, seeded random, and every eligible single X fault on `(T,Ly,Lz)=(2,2,2)`, require bit-for-bit equality with existing `toric3d_x_feedback_pair`, including frame, syndrome, recovery, residual logicals, history timing, and both `:one_step`/`:full_path` tie behavior.

- [ ] **Step 2: Write the failing full-prefix tests**

  Require zero noise to remain clean; every intermediate history entry to contain pre-syndrome, validated recovery, ordered map, and zero post-plaquette/frontier syndrome; future qubits to remain unchanged at each checkpoint; and every eligible single fault to reach a valid final residual. Inject a prefix logical X and require its parity to enter final scoring without truth information.

- [ ] **Step 3: Pin edge cases from Review Focus**

  Assert `T=1` gives `:full_prefix == :none`, `p=0` can decode a manually nonzero syndrome through the clamped prior, invalid protocol/routing/dimensions throw, and the supplied fault record is unchanged by all three replays.

- [ ] **Step 4: Run the focused X test and observe RED**

  Run the managed Julia environment with `include("test/toric3d/feedback_strength_x.jl")`.

  Expected: FAIL because `Toric3DXFeedbackDecoderFamily` is undefined.

- [ ] **Step 5: Implement cached X decoder-family construction**

  Build the existing active-slice matcher once and one BP+OSD decoder for each encoded prefix `t=1:T`; use the `t=T` member as the final decoder.

- [ ] **Step 6: Implement three-protocol X replay**

  Delegate `:none` and `:slice_local` to the existing replay. For `:full_prefix`, propagate checkpoint segments, decode/apply mapped prefix X recoveries only at `t=1:T-1`, and call the final outcome exactly once at `t=T`.

- [ ] **Step 7: Export/include and verify GREEN**

  Run the new tests plus `test/toric3d/x_feedback_replay.jl` and `test/toric3d/x_feedback_decoder.jl`.

  Expected: all tests pass and old/new baseline outcomes are bit-identical.

- [ ] **Step 8: Commit X replay**

  ```bash
  git add src/toric3d/feedback_strength_x.jl src/TopoNoise.jl \
    test/toric3d/feedback_strength_x.jl test/toric3d/runtests.jl
  git commit -m "feat: replay full-prefix diagonal X feedback"
  ```

### Task 4: Add scalar Z fault replay for all three protocols

**Files:**
- Create: `src/toric3d/feedback_strength_z.jl`
- Modify: `src/TopoNoise.jl`
- Create: `test/toric3d/feedback_strength_z.jl`
- Modify: `test/toric3d/runtests.jl`

**Interfaces:**
- Consumes: Tasks 1-2 decoder/window families, `Toric3DFaultRecord`, `sample_toric3d_faults`, `propagate_toric3d_segment!`, `toric3d_z_threshold_outcome`.
- Produces: `Toric3DDiagonalZFaultRecord`, `sample_toric3d_diagonal_z_faults(rng::Random.AbstractRNG, encoder::Toric3DSliceEncoder; p::Real)::Toric3DDiagonalZFaultRecord`, `replay_toric3d_z_strength(encoder::Toric3DSliceEncoder, faults::Toric3DDiagonalZFaultRecord, decoders::Toric3DZFeedbackDecoderFamily; protocol::Symbol, record_history::Bool=false)::NamedTuple`, and `toric3d_z_strength_triplet(encoder::Toric3DSliceEncoder, faults::Toric3DDiagonalZFaultRecord, decoders::Toric3DZFeedbackDecoderFamily; record_history::Bool=false)::NamedTuple{(:none,:slice_local,:full_prefix)}`.

- [ ] **Step 1: Write failing fault-record tests**

  Require exactly `3LyLz` seed bits plus every eligible transfer output, deterministic seeded sampling, no X component, `p=0` all zero, `p=1` all one, immutability under replay, and dimension/rate validation.

- [ ] **Step 2: Write the failing seed/window/prefix replay tests**

  Require local and full protocols to decode the same seed graph, local recovery support only on `r-1:r`, full recovery support only in the prefix, zero complete prefix star syndrome after every feedback checkpoint, and a single common final MWPM call.

- [ ] **Step 3: Exhaust every small-code fault location**

  On `T=2,Ly=Lz=2`, inject every seed and gate-output Z fault and validate all three outcomes. On one `L=3` transfer, assert exactly 54 of 117 faults have odd two-sheet syndrome while all 117 are syndrome-complete in the three-sheet band. Add randomized multifault histories that preserve the full-prefix zero-syndrome invariant after each local/full checkpoint.

- [ ] **Step 4: Pin Z edge cases**

  Assert `T=1` gives `:slice_local == :full_prefix`, zero syndrome does not erase injected logical parity, invalid recovery/protocol/dimensions throw, and final-transfer faults bypass intermediate feedback and reach the common final decoder.

- [ ] **Step 5: Run the focused Z test and observe RED**

  Run the managed Julia environment with `include("test/toric3d/feedback_strength_z.jl")`.

  Expected: FAIL because `Toric3DDiagonalZFaultRecord` is undefined.

- [ ] **Step 6: Implement scalar Z sampling and replay**

  Sample seed faults before transfer faults, maintain one Z Pauli frame, perform seed feedback only for local/full protocols, decode at `r=1:T-1`, and invoke the full final decoder exactly once. Histories store measured syndrome, recovery, mapped support, and audited full-prefix post-syndrome.

- [ ] **Step 7: Export/include and verify GREEN**

  Run the new test with `test/toric3d/z_threshold.jl`, `slice_encoder.jl`, and `encoded_slice_geometry.jl`.

  Expected: all exhaustive and randomized assertions pass.

- [ ] **Step 8: Commit scalar Z replay**

  ```bash
  git add src/toric3d/feedback_strength_z.jl src/TopoNoise.jl \
    test/toric3d/feedback_strength_z.jl test/toric3d/runtests.jl
  git commit -m "feat: replay local and prefix diagonal Z feedback"
  ```

### Task 5: Add a batched Z engine and prove scalar equivalence

**Files:**
- Create: `src/toric3d/feedback_strength_z_batch.jl`
- Modify: `src/TopoNoise.jl`
- Create: `test/toric3d/feedback_strength_z_batch.jl`
- Modify: `test/toric3d/runtests.jl`

**Interfaces:**
- Consumes: Task 4 scalar records/replay and Task 2 batch decoders.
- Produces: `Toric3DDiagonalZFaultBatch`, `sample_toric3d_diagonal_z_fault_batch(rng::Random.AbstractRNG, encoder::Toric3DSliceEncoder; p::Real, shots::Integer)::Toric3DDiagonalZFaultBatch`, and `toric3d_z_strength_batch(encoder::Toric3DSliceEncoder, faults::Toric3DDiagonalZFaultBatch, decoders::Toric3DZFeedbackDecoderFamily)::NamedTuple{(:none,:slice_local,:full_prefix)}` with one logical-1, logical-2, and union `BitVector` per protocol.

- [ ] **Step 1: Write the failing sampling-order equivalence test**

  With identical `Xoshiro` seeds, require batch columns to equal sequential scalar records exactly, including seed-first and operation/slot iteration order. Test batch sizes `1,7,16`, `p=0`, and `p=1`.

- [ ] **Step 2: Write the failing replay equivalence test**

  For seeded random faults on `(T,Ly,Lz)=(2,2,3)`, compare every batch protocol/logical flag with scalar replay column by column. Require batch feedback recoveries to pass the same syndrome validators.

- [ ] **Step 3: Run the batch test and observe RED**

  Expected: FAIL because `Toric3DDiagonalZFaultBatch` is undefined.

- [ ] **Step 4: Implement vectorized checkpoint propagation and decoding**

  Store `n x shots` Z frames for each protocol, propagate the fixed Clifford operations over columns, inject batched faults, and call PyMatching `decode_batch` once per active checkpoint/protocol. Do not vectorize X BP+OSD, whose backend has no batch API.

- [ ] **Step 5: Verify GREEN and benchmark the test fixture**

  Run the new test twice, once normally and once with `JULIA_NUM_THREADS=1`; require identical flags. Print no benchmark output in the test suite.

- [ ] **Step 6: Commit the batched engine**

  ```bash
  git add src/toric3d/feedback_strength_z_batch.jl src/TopoNoise.jl \
    test/toric3d/feedback_strength_z_batch.jl test/toric3d/runtests.jl
  git commit -m "perf: batch diagonal Z feedback replay"
  ```

### Task 6: Add joint statistics and resumable paired-triplet scans

**Files:**
- Create: `src/toric3d/feedback_strength_scan.jl`
- Modify: `src/TopoNoise.jl`
- Create: `test/toric3d/feedback_strength_scan.jl`
- Modify: `test/toric3d/runtests.jl`

**Interfaces:**
- Consumes: Tasks 3 and 5 triplet replay, `_toric3d_wilson`, `_toric3d_sizes`, `_scan_error_rates`, `Toric3DDecoderSettings`.
- Produces: `Toric3DProtocolEstimate`, `Toric3DPairEstimate`, `Toric3DTripleMetricEstimate`, `Toric3DDiagonalFeedbackPoint`, `Toric3DDiagonalFeedbackScan`, and `compare_toric3d_diagonal_feedback(rng::Random.AbstractRNG; sizes, error_rates, channels=[:x_only,:z_only], shots::Integer, seed=nothing, batch_size::Integer=256, x_routing::Symbol=:one_step, settings::Toric3DDecoderSettings=Toric3DDecoderSettings(), output_dir=nothing, resume::Bool=false, progress_io=nothing)`.

- [ ] **Step 1: Write failing eight-pattern statistic tests**

  Encode joint outcome mask bits as `none=1`, `slice_local=2`, `full_prefix=4`. With hand-computed counts, require eight patterns to reproduce every marginal, Wilson interval, three pairwise discordance tables, rate difference, standard error, and conservative interval. Reject negative counts or totals different from `shots`.

- [ ] **Step 2: Write failing scan tests**

  On `(2,2,2)`, rates `[0,.01]`, both channels, 20 triplets, and batch size 7, require six protocol/channel marginals, three logical metrics, exact joint totals, zero-noise success, one X decoder family per `(size,p)`, one Z decoder family per size, progress output, and deterministic repeated seeds.

- [ ] **Step 3: Write failing interruption/resume tests**

  Interrupt after one internal batch, resume, and require byte-identical joint counts to uninterrupted execution. Separately change dimensions, rates, channels, shots, batch size, routing, seed metadata, and one BP+OSD setting; each mismatch must throw before a new fault draw.

- [ ] **Step 4: Write the invalid-recovery persistence test**

  Through a private test hook, return a wrong-syndrome recovery. Require an atomic failure TOML containing dimensions, channel, `p`, protocol, checkpoint, shot/batch, point seed, and error text before the scan aborts. Completed metadata state must retain invalid count zero.

- [ ] **Step 5: Run the scan test and observe RED**

  Expected: FAIL because `Toric3DTripleMetricEstimate` and the scan API are undefined.

- [ ] **Step 6: Implement statistics, point substreams, and checkpoint state**

  Store one `3 x 8` joint-count table per point, deterministic `UInt64` point seeds, completed shots, elapsed phase times, decoder versions, invalid count, and complete configuration. Resume reconstructs the point RNG and discards exactly the recorded fault records before continuing.

- [ ] **Step 7: Implement channel dispatch and incremental saving hooks**

  Use scalar X triplets and batched Z triplets. Build/copy decoder families outside shot loops. Save internal checkpoint state after every batch and finalize a point only when all joint-pattern totals equal `shots`.

- [ ] **Step 8: Export/include and verify GREEN**

  Run the focused scan test twice, including `JULIA_NUM_THREADS=1`, then run existing `x_feedback_scan.jl` and `z_threshold.jl` tests.

- [ ] **Step 9: Commit the scan API**

  ```bash
  git add src/toric3d/feedback_strength_scan.jl src/TopoNoise.jl \
    test/toric3d/feedback_strength_scan.jl test/toric3d/runtests.jl
  git commit -m "feat: scan three diagonal feedback protocols"
  ```

### Task 7: Add artifacts, figures, and the runnable CLI

**Files:**
- Create: `src/toric3d/feedback_strength_output.jl`
- Create: `examples/compare_toric3d_diagonal_feedback.jl`
- Modify: `src/TopoNoise.jl`
- Create: `test/toric3d/feedback_strength_output.jl`
- Create: `test/toric3d/feedback_strength_example.jl`
- Modify: `test/toric3d/runtests.jl`
- Modify: `README.md`

**Interfaces:**
- Consumes: `Toric3DDiagonalFeedbackScan` and Task 6 statistics/checkpoint fields.
- Produces: `save_toric3d_diagonal_feedback(scan::Toric3DDiagonalFeedbackScan, directory::AbstractString; basename::AbstractString="toric3d-diagonal-feedback")::NamedTuple`, `plot_toric3d_diagonal_feedback(scan::Toric3DDiagonalFeedbackScan)::CairoMakie.Figure`, and `Toric3DDiagonalFeedbackExample.main(arguments=ARGS)::Int`.

- [ ] **Step 1: Write failing CSV/TOML/failure-artifact tests**

  Require one CSV row per point/metric with all eight joint patterns, three marginals, Wilson intervals, pairwise discordances/differences, seeds, routing, phase timings, and versions. Require metadata fields from the spec, source/Git provenance, completion counts, zero invalid recoveries, and exact X/Z seed fault-location differences.

- [ ] **Step 2: Write failing figure tests**

  Require nonempty PDF/SVG/PNG exports. Inspect the Makie figure object for channel facets, all three protocol labels, individual logicals, paired differences, and coincident-curve annotation support.

- [ ] **Step 3: Write failing CLI tests**

  Pin `--help`, unknown/missing arguments, `--sizes`, `--rates`, `--channels`, `--shots`, `--batch-size`, `--seed`, `--x-routing`, `--output`, `--resume`, `--timing-shots`, and `--pilot-only`. A zero-noise smoke run must exclude warmup/timing shots from CSV counts and create `timing-pilot.toml`.

- [ ] **Step 4: Run output/CLI tests and observe RED**

  Expected: FAIL because save, plot, and CLI interfaces are undefined.

- [ ] **Step 5: Implement atomic outputs and plots**

  Reuse `_toric3d_atomic_write` and the existing Wilson conventions. Persist a recovery-failure record before rethrowing. Plot marginal state failure, individual logicals, and the three paired differences; label X one-step versus Z full-window local policies explicitly.

- [ ] **Step 6: Implement the CLI and timing pilot**

  Compile with separate warmup records, run independent timing triplets, record decoder construction and each protocol phase, estimate scout runtime, then run/resume requested scientific points. The default invocation runs the approved 200-shot X scout; selecting only `z_only` without explicit sizes/rates selects the approved Z scout. Selecting both channels requires explicit common sizes/rates. `--pilot-only` stops after timing.

- [ ] **Step 7: Document the public API and run GREEN checks**

  Add a concise README example and interpretation warning. Run output/CLI tests, the zero-noise CLI smoke run, and existing X-output/example tests.

- [ ] **Step 8: Commit artifacts and CLI**

  ```bash
  git add src/toric3d/feedback_strength_output.jl src/TopoNoise.jl \
    examples/compare_toric3d_diagonal_feedback.jl README.md \
    test/toric3d/feedback_strength_output.jl \
    test/toric3d/feedback_strength_example.jl test/toric3d/runtests.jl
  git commit -m "feat: save diagonal feedback strength scans"
  ```

### Task 8: Add guarded three-protocol threshold analysis

**Files:**
- Create: `examples/fit_toric3d_diagonal_feedback.py`
- Create: `test/toric3d/test_feedback_strength_fit.py`

**Interfaces:**
- Consumes: Task 7 raw CSV state-failure rows.
- Produces: `read_feedback_points(paths: Sequence[Path])`, `fit_feedback_protocol`, `bootstrap_feedback_protocol`, `analyze_feedback_scan`, CLI JSON/CSV reports, and PDF/SVG/PNG figures. The CLI accepts one or more disjoint shard CSVs and rejects duplicate `(dimensions,p,channel,metric)` keys.

- [ ] **Step 1: Write failing parser and synthetic-fit tests**

  Generate seeded three-protocol data for two channels with known crossings. Require correct marginal extraction from the eight patterns, parameter recovery within Monte Carlo tolerance, reproducible 1,000-replicate bootstrap, and separate results for all six channel/protocol combinations.

- [ ] **Step 2: Write failing guard tests**

  Require explicit rejection for fewer than four sizes, missing adjacent-size crossings, too few probabilities, optimizer failure, parameter-bound hits, nonmonotone fits, `p_c` outside the window, deviance goodness-of-fit p-value below 0.01, and fewer than 80% usable bootstrap fits.

- [ ] **Step 3: Write failing sensitivity/output tests**

  Require quadratic primary plus cubic, minimum-size, probability-window, and fixed `L^-1` sensitivity rows; keep bootstrap intervals separate from accepted-model spread. Require paired rate-difference panels without pretending paired point counts make protocol fits independent.

- [ ] **Step 4: Run the Python test and observe RED**

  Run `PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/toponoise-mpl .CondaPkg/.pixi/envs/default/bin/python -m unittest test/toric3d/test_feedback_strength_fit.py -v`.

  Expected: FAIL because the fitter module does not exist.

- [ ] **Step 5: Implement the guarded fitter and bootstrap**

  Adapt the accepted existing X/Z logistic likelihood code without weakening guards. Report `status` plus diagnostics for every protocol, including `no_crossing` and rejected candidate values without promoting them to thresholds.

- [ ] **Step 6: Implement reports and figures, then verify GREEN**

  Save full input/source hashes, arguments, versions, accepted/rejected fits, bootstrap samples, sensitivity CSV, fit curves, and scaled collapse only for accepted fits. Rerun the Python test twice and require byte-stable JSON apart from declared timestamps, which should be omitted.

- [ ] **Step 7: Commit statistical analysis**

  ```bash
  git add examples/fit_toric3d_diagonal_feedback.py \
    test/toric3d/test_feedback_strength_fit.py
  git commit -m "feat: fit diagonal feedback thresholds"
  ```

### Task 9: Run complete verification before scientific scans

**Files:**
- Modify only if a verified defect is found, following a new red/green test cycle.

**Interfaces:**
- Consumes: Tasks 1-8.
- Produces: a verified implementation ready for timed scientific execution.

- [ ] **Step 1: Run the focused 3D Julia suite**

  Run the managed Julia environment with `julia --project=. -e 'using Test, TopoNoise; include("test/toric3d/runtests.jl")'`.

  Expected: zero failures and no invalid-recovery warnings.

- [ ] **Step 2: Run all Python analysis tests**

  Run the four existing toric3d Python test files plus `test_feedback_strength_fit.py` with `unittest -v`.

  Expected: zero failures.

- [ ] **Step 3: Run the complete repository suite**

  Run the managed Julia environment with `julia --project=. test/runtests.jl`.

  Expected: zero failures across existing 2D and 3D tests.

- [ ] **Step 4: Inspect the branch diff and repair only demonstrated issues**

  Check `git diff --check`, exported names, checkpoint compatibility, and that generated LaTeX files remain unstaged. If any defect appears, write and observe a failing regression test before changing production code, then repeat Steps 1-3.

- [ ] **Step 5: Commit any verification-only fixes**

  Skip this commit when no files changed. Otherwise use a focused message naming the demonstrated defect.

### Task 10: Run timing pilots and broad paired scouts

**Files:**
- Create: `results/toric3d-diagonal-feedback-2026-10-12/timing/*`
- Create: `results/toric3d-diagonal-feedback-2026-10-12/scout-x/*`
- Create: `results/toric3d-diagonal-feedback-2026-10-12/scout-z/*`
- Create: `results/toric3d-diagonal-feedback-2026-10-12/scan-selection.toml`

**Interfaces:**
- Consumes: Task 7 CLI and Task 8 fitter.
- Produces: runtime measurements, crossing brackets, and a frozen production grid.

- [ ] **Step 1: Run excluded 20-triplet timing points**

  Run X at `L=9,p=.015` and Z at `L=17,p=.006` with `--pilot-only --timing-shots 20`. Save setup and protocol phase times; confirm timing shots do not appear in scientific CSV.

- [ ] **Step 2: Record the production-runtime estimate before starting scouts**

  Save per-channel estimates for the approved 200-shot grids. If one point exceeds 60 seconds, run it as an incremental checkpointed command and report progress at least once per minute rather than holding a blocking call.

- [ ] **Step 3: Run the 200-triplet X scout**

  Use `L={3,5,7,9}` and `p={0,.0001,.0003,.001,.003,.005,.0075,.01,.015,.02,.03,.05,.075,.10}` with point substreams and incremental saving.

- [ ] **Step 4: Run the 200-triplet Z scout**

  Use `L={5,9,13,17}` and `p={0,.001,.002,.003,.004,.005,.006,.0075,.01,.015,.02,.03}` with the same checkpoint discipline.

- [ ] **Step 5: Extend only unbracketed transitions**

  If a protocol has no adjacent-size sign change, add logarithmic/lattice probabilities in the indicated direction and rerun all three paired protocols at those points. Do not fit until each reported threshold has at least one crossing bracket.

- [ ] **Step 6: Freeze the production selection**

  Write `scan-selection.toml` with input hashes, timing, scout grids, observed brackets, and the union critical window. Select X sizes `5,7,9,11,13`, Z sizes `9,13,17,21,25`, and at least seven probabilities spanning the union of protocol brackets with one simulated guard point on each side. Record any protocol lacking a supported production window as `no_crossing` rather than inventing one.

- [ ] **Step 7: Commit scout artifacts and selection**

  ```bash
  git add results/toric3d-diagonal-feedback-2026-10-12
  git commit -m "data: scout diagonal feedback strengths"
  ```

### Task 11: Run production, fit thresholds, and deliver the scientific report

**Files:**
- Create: `results/toric3d-diagonal-feedback-2026-10-12/production/*`
- Create: `results/toric3d-diagonal-feedback-2026-10-12/analysis/*`
- Create: `docs/superpowers/toric3d-diagonal-feedback-strength.md`
- Modify: `document/note.tex`

**Interfaces:**
- Consumes: the frozen Task 10 grid and Tasks 7-8 tooling.
- Produces: final paired CSV/metadata/checkpoints, threshold reports, PDF/SVG/PNG figures, corrected note, and reproducible conclusions.

- [ ] **Step 1: Run checkpointed production scans**

  Start with 1,000 triplets per X point and 2,000 per Z point on the frozen grid. Shard only by disjoint `(channel,size)` sets with recorded point seeds; never pool warmup, timing, or scout shots. Save every completed batch.

- [ ] **Step 2: Audit raw production data before fitting**

  Require every point/metric joint table to sum to shots, recompute all marginals and discordances independently, verify zero-noise success and zero invalid recoveries, and hash every raw/metadata input. Stop on the first mismatch.

- [ ] **Step 3: Run guarded fits and precision top-ups**

  Run 1,000 bootstrap replicates per eligible channel/protocol. Top up only critical points whose Wilson intervals or bootstrap acceptance prevent the preregistered conclusion, preserving their point seeds and checkpoint counts. Do not top up a rejected no-crossing protocol merely to force a fit.

- [ ] **Step 4: Generate and visually inspect every figure**

  Export protocol curves, individual logicals, paired differences, accepted fit/collapse panels, and sensitivity plots as PDF/SVG/PNG. Render each PDF to PNG and inspect labels, coincident curves, clipping, legends, and uncertainty bars.

- [ ] **Step 5: Write the report and update the note**

  Report every accepted threshold with conditional interval and model spread, every rejected/no-crossing result with its diagnostic, paired protocol differences, runtime, resource scaling, seed asymmetry, circuit-correlation limitations, and exact reproduction commands. Update the 3D note without overwriting its existing corrected geometry/results.

- [ ] **Step 6: Compile and inspect the note**

  Run `latexmk -pdf note.tex` from `document`, render the changed pages, and inspect them. Stage `document/note.tex` explicitly with `git add -f`; do not stage `.fdb_latexmk`, `.fls`, `.synctex.gz`, or other generated files.

- [ ] **Step 7: Re-run final verification against delivered artifacts**

  Repeat Task 9's complete Julia/Python suites, run the fitter reproduction command from the report, verify saved hashes, and run `git diff --check`.

- [ ] **Step 8: Commit production and documentation**

  ```bash
  git add results/toric3d-diagonal-feedback-2026-10-12 \
    docs/superpowers/toric3d-diagonal-feedback-strength.md
  git add -f document/note.tex
  git commit -m "data: compare diagonal feedback strengths"
  ```

### Task 12: Request final whole-branch review and integrate findings

**Files:**
- Modify only files required by verified review findings.

**Interfaces:**
- Consumes: all implementation, tests, results, and documentation.
- Produces: a reviewed branch ready for the user's chosen integration action.

- [ ] **Step 1: Dispatch one fresh whole-branch reviewer**

  Ask the reviewer to compare every spec requirement with the branch, inspect decoder information boundaries and checkpoint invariants, audit statistics/results, and rank findings by severity. The reviewer must not edit files.

- [ ] **Step 2: Verify each finding independently**

  Reproduce actionable findings. For code defects, add a failing test before a fix; for analysis defects, add a failing synthetic or artifact audit.

- [ ] **Step 3: Apply verified fixes and rerun affected plus full verification**

  Run the narrow red/green check, Task 9 complete suites, artifact audits, and note compilation as applicable.

- [ ] **Step 4: Commit review fixes and prepare the integration handoff**

  Use focused commits per verified finding. Summarize final thresholds/statuses, validation evidence, material limitations, commits, and the remaining choice to merge, push, or keep the branch.
