# Encoded Diagonal-Sweep Logical-X Feedback Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add and run a reproducible logical-X threshold comparison for the encoded diagonal 3D toric-code sweep with ideal within-slice measurement feedback.

**Architecture:** Reuse `Toric3DEncodedSlice`, the seven-gate diagonal transfer, and the existing active-slice MWPM router. Prepare logical `|00>`, replay each X-only gate-fault record through final-only and one-step-feedback protocols, validate a final BP+OSD X recovery, and score the two residual logical-X membrane parities. Keep simulation, output, statistical fitting, and scientific documentation in focused files.

**Tech Stack:** Julia 1.12, TopoNoise.jl binary Pauli frames, PythonCall, ldpc BP+OSD, PyMatching, CairoMakie, Python 3, NumPy/SciPy/Matplotlib, unittest, LaTeX.

**Spec:** `docs/superpowers/specs/2026-10-11-toric3d-diagonal-x-feedback-design.md`

## Global Constraints

- Preserve the existing rough-slab, fixed-boundary slice, and encoded logical-Z public interfaces and results.
- Use the existing diagonal geometry with `n=3(T+1)LyLz`, `k=2`, `dX=(T+1)min(Ly,Lz)`, and `dZ=2min(Ly,Lz)`.
- Prepare logical `|00>` from product zero using ideal X-star measurements and ideal Z correction; add no seed H/CNOT fault locations.
- Apply independent rate-`p` X faults after all thirteen participating H/CNOT outputs per transfer cell; initialization, idles, measurements, classical processing, and corrections remain ideal.
- Replay identical fault records across final-only and feedback protocols.
- Use one-step active-slice MWPM routing as the primary feedback and full-path routing only as a labeled sensitivity.
- Use minimum-sum BP, 50 iterations, scaling factor 1, parallel schedule, one thread, and `OSD_CS` order 2 for final X recovery.
- Reject any active-slice or final recovery that does not reproduce its input syndrome.
- Treat either residual logical-X membrane as logical `|00>` preparation failure.
- Record that uniform data priors neglect nonuniform exposure and circuit/feedback correlations.
- Save completed points incrementally; never pool warmup or timing-pilot shots into production counts.
- Use the existing managed decoder interpreter through `JULIA_CONDAPKG_BACKEND=Null` and `JULIA_PYTHONCALL_EXE=$PWD/.CondaPkg/.pixi/envs/default/bin/python` on this machine.
- The current workspace contains the uncommitted 3D base required by this feature. Execution must preserve it and must not start from a clean `HEAD` worktree unless that base has first been checkpointed.

## Review Focus

- `p=0` with a manually supplied nonzero syndrome must still decode through a clamped finite BP prior; Task 2 pins this behavior.
- Periodic `Ly=2` or `Lz=2` must preserve distinct parallel physical edges and both logical parities; Tasks 1 and 2 cover these geometries.
- Odd, wrong-length, or inconsistent active/final syndromes must raise a descriptive error rather than produce a score; Tasks 2 and 3 cover each case.
- Resuming a checkpoint with different dimensions, probabilities, routing, shots, decoder settings, or seed provenance must fail before sampling; Task 4 covers configuration mismatches.
- Missing crossings, optimizer-bound fits, nonmonotone fitted curves, and poor fit diagnostics must be reported instead of converted into a threshold; Task 6 covers each guard.

---

### Task 1: Verify the encoded logical-`|00>` boundary and transfer

**Files:**
- Modify: `src/toric3d/encoded_slice_geometry.jl`
- Modify: `src/TopoNoise.jl`
- Create: `test/toric3d/encoded_slice_zero.jl`
- Modify: `test/toric3d/runtests.jl`

**Interfaces:**
- Consumes: `Toric3DEncodedSlice`, `toric3d_encoded_seed_checks`, `toric3d_encoded_logical_z`, `toric3d_encoded_encoder`.
- Produces: `verify_toric3d_encoded_zero_encoder(code::Toric3DEncodedSlice, encoder::Toric3DSliceEncoder=toric3d_encoded_encoder(code))::Bool`.

- [x] **Step 1: Write the failing logical-zero verification tests**

Add `@testset "Encoded transfer prepares logical zero at every checkpoint"` covering `(1,2,2)`, `(2,2,3)`, and `(3,3,3)`. Assert that the new verifier accepts the valid encoder, rejects a missing operation, rejects shifted checkpoints, and rejects a changed first gate. In the same file, assert seed X/Z ranks, both logical-Z cycles outside the seed Z-check row space, and identity anticommutation with the logical-X membranes for periodic length two.

- [x] **Step 2: Run the focused test and observe the missing interface**

Run:

```bash
JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" \
  julia --project=. -e 'using Test, TopoNoise; include("test/toric3d/encoded_slice_zero.jl")'
```

Expected: FAIL because `verify_toric3d_encoded_zero_encoder` is undefined.

- [x] **Step 3: Implement logical-zero signed stabilizer verification**

Generalize the existing pullback helper only enough to select the seed stabilizer group. For logical zero, accept seed X operators in the row space of the seed star matrix and seed Z operators in the span of the compatibility checks plus both logical-Z cycles. At every prefix verify final X stars, physical Z plaquettes, both logical Z operators, and the active compatibility checks. Preserve `verify_toric3d_encoded_encoder` and its logical-plus behavior unchanged.

- [x] **Step 4: Export the verifier and include the new test**

Add `verify_toric3d_encoded_zero_encoder` to the encoded-slice export list and include `encoded_slice_zero.jl` immediately after `encoded_slice_geometry.jl` in the 3D test runner.

- [x] **Step 5: Run the encoded geometry tests**

Run the Task 1 command and then:

```bash
JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" \
  julia --project=. -e 'using Test, TopoNoise; include("test/toric3d/encoded_slice_geometry.jl"); include("test/toric3d/encoded_slice_zero.jl")'
```

Expected: all assertions pass, including `Ly=2` and `Lz=2` cases.

- [x] **Step 6: Commit the algebraic boundary support**

```bash
git add src/toric3d/encoded_slice_geometry.jl src/TopoNoise.jl \
  test/toric3d/encoded_slice_zero.jl test/toric3d/runtests.jl
git commit -m "feat: verify diagonal logical zero preparation"
```

### Task 2: Add validated final BP+OSD X recovery

**Files:**
- Create: `src/toric3d/x_feedback_threshold.jl`
- Modify: `src/TopoNoise.jl`
- Create: `test/toric3d/x_feedback_decoder.jl`
- Modify: `test/toric3d/runtests.jl`

**Interfaces:**
- Consumes: `Toric3DEncodedSlice`, `Toric3DDecoderSettings`, `plaquette_check_matrix`, `toric3d_encoded_logical_z`, `_toric3d_load_ldpc`, `_toric3d_py_bits`.
- Produces: `Toric3DXThresholdDecoder`, `build_toric3d_x_threshold_decoder(code::Toric3DEncodedSlice; p::Real, settings::Toric3DDecoderSettings=Toric3DDecoderSettings())`, `decode_toric3d_x(decoder::Toric3DXThresholdDecoder, plaquettes::AbstractVector{Bool})::BitVector`, and `toric3d_x_threshold_outcome(decoder::Toric3DXThresholdDecoder, frame::AbstractVector{Bool})::NamedTuple{(:syndrome,:recovery,:residual,:logical_1,:logical_2,:failure),Tuple{BitVector,BitVector,BitVector,Bool,Bool,Bool}}`.

- [x] **Step 1: Write failing decoder and logical-scoring tests**

For `Toric3DEncodedSlice(2,2,3)`, assert the builder uses the full plaquette matrix, default settings, and a prior clamped away from zero. Decode every single-qubit X syndrome and assert `H_Z*r == syndrome`. Inject each X logical membrane and assert only its matching logical parity fails; inject one Z-check stabilizer and assert no logical failure. Test wrong-length syndrome and `_toric3d_validate_x_recovery` on wrong-length and syndrome-invalid vectors; retain the existing `_toric3d_py_bits` nonbinary/length tests for Python output conversion.

- [x] **Step 2: Run the focused test and verify it fails**

Run:

```bash
JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" \
  julia --project=. -e 'using Test, TopoNoise; include("test/toric3d/x_feedback_decoder.jl")'
```

Expected: FAIL because the X-threshold decoder API is undefined.

- [x] **Step 3: Implement `Toric3DXThresholdDecoder` and its builder**

Store the encoded code, declared `p`, settings, full binary plaquette matrix, two-row logical-Z matrix, Python BP+OSD object, and ldpc version. Construct `ldpc.BpOsdDecoder` with `error_rate=clamp(p,1e-12,1-1e-12)` and the exact settings in Global Constraints.

- [x] **Step 4: Implement validated decode and logical scoring**

`decode_toric3d_x` returns all-zero recovery for an all-zero syndrome; otherwise convert the backend result with `_toric3d_py_bits` and require `_toric3d_matrix_syndrome(H_Z,recovery)==syndrome`. `toric3d_x_threshold_outcome` measures, decodes, applies recovery, pairs the residual X frame with both logical-Z cycles, and returns the measured syndrome and validated recovery together with that residual and the three failure flags.

- [x] **Step 5: Export/include the API and run the decoder tests**

Include `x_feedback_threshold.jl` after `z_threshold.jl`; export the four public names from the Interfaces block; include the test in `test/toric3d/runtests.jl`. Run the Task 2 command.

Expected: all tests pass, including `p=0` with a nonzero manual syndrome.

- [x] **Step 6: Commit the final X decoder**

```bash
git add src/toric3d/x_feedback_threshold.jl src/TopoNoise.jl \
  test/toric3d/x_feedback_decoder.jl test/toric3d/runtests.jl
git commit -m "feat: decode encoded diagonal X errors"
```

### Task 3: Replay paired X faults through preparation feedback

**Files:**
- Modify: `src/toric3d/x_feedback_threshold.jl`
- Create: `test/toric3d/x_feedback_replay.jl`
- Modify: `test/toric3d/runtests.jl`

**Interfaces:**
- Consumes: `toric3d_encoded_encoder`, `build_toric3d_slice_decoder`, `sample_toric3d_faults`, `replay_toric3d_slice`, `Toric3DXThresholdDecoder`.
- Produces: `replay_toric3d_x_feedback(encoder::Toric3DSliceEncoder, faults::Toric3DFaultRecord, slice_decoder::Toric3DSliceDecoder, final_decoder::Toric3DXThresholdDecoder; feedback::Bool, routing::Symbol=:one_step, record_history::Bool=false)::NamedTuple` and `toric3d_x_feedback_pair(encoder::Toric3DSliceEncoder, faults::Toric3DFaultRecord, slice_decoder::Toric3DSliceDecoder, final_decoder::Toric3DXThresholdDecoder; routing::Symbol=:one_step, record_history::Bool=false)::NamedTuple{(:final_only,:slice_feedback)}`.

- [x] **Step 1: Write failing paired-replay tests**

Assert zero faults give zero final syndrome and no logical failure for both protocols and both routes. For every eligible post-H/CNOT output in `(T,Ly,Lz)=(2,2,2)`, inject a single X fault and assert: no Z component appears; each active MWPM recovery and final BP+OSD recovery is syndrome-valid; both protocols consume the same immutable fault record; and repeating the replay is deterministic. Add errors for odd active syndrome, wrong dimensions, and invalid routing.

- [x] **Step 2: Run the focused test and verify it fails**

Run the same Julia environment command with `include("test/toric3d/x_feedback_replay.jl")`.

Expected: FAIL because `replay_toric3d_x_feedback` is undefined.

- [x] **Step 3: Implement the single-protocol replay wrapper**

Require the encoded code and plain slice encoder/decoder dimensions to agree. Delegate gate propagation and target feed-forward to `replay_toric3d_slice`, then pass its final X frame to `toric3d_x_threshold_outcome`. Return the frame, final syndrome, recovery, logical flags, and optional active-slice history needed for validation.

- [x] **Step 4: Implement paired replay**

Call the single-protocol wrapper twice with the same immutable fault record: `feedback=false` and `feedback=true`. Do not sample faults or inspect them inside either decoder.

- [x] **Step 5: Run the replay and existing slice-feedback tests**

Run:

```bash
JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" \
  julia --project=. -e 'using Test, TopoNoise; include("test/toric3d/slice_feedback.jl"); include("test/toric3d/x_feedback_replay.jl")'
```

Expected: all tests pass.

- [x] **Step 6: Commit paired preparation replay**

```bash
git add src/toric3d/x_feedback_threshold.jl \
  test/toric3d/x_feedback_replay.jl test/toric3d/runtests.jl
git commit -m "feat: replay diagonal X feedback protocols"
```

### Task 4: Add the resumable paired scan and statistics

**Files:**
- Create: `src/toric3d/x_feedback_scan.jl`
- Modify: `src/TopoNoise.jl`
- Create: `test/toric3d/x_feedback_scan.jl`
- Modify: `test/toric3d/runtests.jl`

**Interfaces:**
- Consumes: `toric3d_x_feedback_pair`, `_toric3d_metric`, `_toric3d_sizes`, `_scan_error_rates`, `Toric3DDecoderSettings`.
- Produces: `Toric3DXFeedbackPoint`, `Toric3DXFeedbackScan`, and `compare_toric3d_x_feedback_threshold(rng::Random.AbstractRNG; sizes, error_rates, shots::Integer, seed=nothing, batch_size::Integer=256, routing::Symbol=:one_step, settings::Toric3DDecoderSettings=Toric3DDecoderSettings(), output_dir=nothing, resume::Bool=false, progress_io=nothing)`.

- [x] **Step 1: Write failing scan tests**

Use `(2,2,2)`, rates `[0,.01]`, 20 paired shots, and `batch_size=7`. Assert zero-noise success, exact counts and paired discordances, identical results for repeated seeded runs, one point per `(dimensions,p)`, one final decoder per point, and a progress line per point. Interrupt a point after one batch through the internal test hook, resume it, and require the completed counts to equal an uninterrupted run. Assert invalid empty/duplicate rates, zero shots, zero batch size, invalid routing, invalid output path, and mismatched resume metadata throw before sampling.

- [x] **Step 2: Run the focused test and verify it fails**

Run the Julia environment command with `include("test/toric3d/x_feedback_scan.jl")`.

Expected: FAIL because the scan API is undefined.

- [x] **Step 3: Implement scan types and paired accumulation**

Store per metric (`logical_1`, `logical_2`, `state_failure`) the existing `Toric3DMetricEstimate`. Build the encoder and active-slice decoder once per size and the final BP+OSD decoder once per `(size,p)`. Sample X-only `Toric3DFaultRecord` values and accumulate final-only, feedback, both, and discordant counts.

- [x] **Step 4: Implement completed-point checkpoint resume**

Before sampling, write a private TOML checkpoint containing configuration metadata and deterministic per-point `UInt64` substream seeds drawn from the supplied RNG. Process a point in `batch_size` chunks and atomically checkpoint completed shots plus the three paired count accumulators after every chunk. With `resume=true`, require all saved configuration fields to match, restore completed points, recreate an interrupted point's RNG from its stored substream seed, call `sample_toric3d_faults` and discard exactly the recorded number of shot records to reproduce the draw position, and continue from the next shot. Never merge a checkpoint from a different configuration. Task 5 promotes completed scan state to the public CSV and metadata formats; it does not own this resumable internal checkpoint.

- [x] **Step 5: Export/include and run the scan tests twice**

Run the focused test once normally and once after setting `JULIA_NUM_THREADS=1` explicitly.

Expected: both runs pass with identical saved counts.

- [x] **Step 6: Commit the scan API**

```bash
git add src/toric3d/x_feedback_scan.jl src/TopoNoise.jl \
  test/toric3d/x_feedback_scan.jl test/toric3d/runtests.jl
git commit -m "feat: scan paired diagonal X feedback"
```

### Task 5: Add CSV/TOML/figures and a runnable CLI

**Files:**
- Create: `src/toric3d/x_feedback_output.jl`
- Modify: `src/TopoNoise.jl`
- Create: `examples/scan_toric3d_x_feedback_threshold.jl`
- Create: `test/toric3d/x_feedback_output.jl`
- Create: `test/toric3d/x_feedback_example.jl`
- Modify: `test/toric3d/runtests.jl`

**Interfaces:**
- Consumes: `Toric3DXFeedbackScan`, `_toric3d_atomic_write`, `_toric3d_csv`.
- Produces: `save_toric3d_x_feedback_threshold(scan::Toric3DXFeedbackScan, directory::AbstractString; basename::AbstractString="toric3d-x-feedback")::NamedTuple`, `plot_toric3d_x_feedback_threshold(scan::Toric3DXFeedbackScan)::CairoMakie.Figure`, and `Toric3DXFeedbackExample.main(arguments=ARGS)::Int`.

- [x] **Step 1: Write failing artifact tests**

For a two-point scan, assert CSV columns include dimensions, `p`, shots, point seed, routing, all three paired metrics, both logical counts, Wilson intervals, discordances, elapsed time, ldpc version, and every BP+OSD setting. Assert metadata describes logical `|00>`, thirteen transfer outputs, no seed gate faults, ideal operations, decoder limitations, completion count, and resume configuration. Assert PDF/SVG/PNG files exist and are nonempty.

- [x] **Step 2: Write failing CLI tests**

Test `--help`, unknown/missing options, invalid sizes/rates/shots/batch-size/routing, a one-shot zero-noise smoke run, timing-pilot exclusion, and all expected artifacts. Require a standalone boolean `--resume`, integer `--batch-size`, and `--routing one_step|full_path` options.

- [x] **Step 3: Run both focused tests and verify they fail**

Run each new Julia test file with the managed decoder environment.

Expected: FAIL because output and CLI interfaces are undefined.

- [x] **Step 4: Implement atomic output and the comparison figure**

Plot final-only and feedback state-failure curves together for each size with 95% Wilson bars. Put both individual logical rates in a secondary panel and paired feedback-minus-final differences in a third panel. Label routing, shots, logical `|00>` state failure, ideal measurement/correction assumptions, and decoder-correlation limitation.

- [x] **Step 5: Implement the CLI and timing pilot**

Default to the approved broad pilot grid, seed 1234, batch size 256, one-step routing, 1,000 shots, and `results/toric3d-x-feedback-threshold-2026-10-11/pilot`. Compile a two-shot warmup, time setup/first-decode separately from the independent timing shots, save `timing-pilot.toml`, print the rough scan estimate, then run and save production points incrementally. Parse `--resume` without a following value while every other option consumes exactly one value.

- [x] **Step 6: Export/include and run artifact/CLI tests**

Expected: all tests pass; zero-noise CSV rows contain zero failures in both protocols.

- [x] **Step 7: Commit output and CLI support**

```bash
git add src/toric3d/x_feedback_output.jl src/TopoNoise.jl \
  examples/scan_toric3d_x_feedback_threshold.jl \
  test/toric3d/x_feedback_output.jl test/toric3d/x_feedback_example.jl \
  test/toric3d/runtests.jl
git commit -m "feat: save diagonal X feedback scans"
```

### Task 6: Add guarded finite-size threshold analysis

**Files:**
- Create: `examples/fit_toric3d_x_feedback_threshold.py`
- Create: `test/toric3d/test_x_feedback_fit.py`

**Interfaces:**
- Consumes: the Task 5 raw CSV, one row per metric and point.
- Produces: `read_points`, `fit_protocol`, `bootstrap_protocol`, guarded CLI analysis, JSON fit report, sensitivity CSV, and PDF/SVG/PNG figures.

- [x] **Step 1: Write failing statistical tests**

Generate seeded synthetic binomial data with known `pc` and `nu` for both protocols. Assert parameter recovery within the Monte Carlo tolerance, bootstrap reproducibility, and separate protocol results. Add explicit rejection tests for no adjacent-size crossing, optimizer parameters at bounds, decreasing-in-`p` fitted curves, too few sizes/probabilities, and deviance outside the configured goodness-of-fit guard.

- [x] **Step 2: Run the Python test and verify it fails**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 .CondaPkg/.pixi/envs/default/bin/python \
  test/toric3d/test_x_feedback_fit.py
```

Expected: FAIL because the fitter module is missing.

- [x] **Step 3: Implement input validation and the quadratic scaling fit**

Read only `metric=state_failure`, reject duplicate/missing points, and fit each requested protocol independently by binomial likelihood with `x=((p-pc)/sp)(L/Lref)^(1/nu)`. Return numerical status, parameters, deviance, degrees of freedom, and raw crossings. Do not emit a threshold when a guard fails.

- [x] **Step 4: Implement bootstrap and sensitivity analysis**

Use a seeded parametric binomial bootstrap conditional on the accepted model. Repeat fits across minimum sizes and probability windows; add the leading fixed-`omega=1` size correction when raw crossings drift monotonically. Keep sampling intervals separate from the acceptable-model spread.

- [x] **Step 5: Implement figures and provenance**

Plot both raw protocol curves, scaled collapse for accepted fits, raw crossing drift, and full-path sensitivity when supplied. Record input and source SHA-256 hashes, command arguments, accepted bootstrap count, package versions, and all rejected sensitivity fits.

- [x] **Step 6: Run the Python tests**

Expected: all tests pass and no `__pycache__` directory appears.

- [x] **Step 7: Commit the fitter**

```bash
git add examples/fit_toric3d_x_feedback_threshold.py \
  test/toric3d/test_x_feedback_fit.py
git commit -m "feat: fit diagonal X feedback threshold"
```

### Task 7: Run the staged simulation and analyze the threshold

**Files:**
- Create: `results/toric3d-x-feedback-threshold-2026-10-11/pilot/*`
- Create: `results/toric3d-x-feedback-threshold-2026-10-11/production/*`
- Create: `results/toric3d-x-feedback-threshold-2026-10-11/full-path/*`
- Create: `results/toric3d-x-feedback-threshold-2026-10-11/analysis/*`

**Interfaces:**
- Consumes: Task 5 CLI and Task 6 fitter.
- Produces: timed pilot, completed raw scans, metadata, fit/sensitivity reports, and publication-ready figures.

- [x] **Step 1: Run the approved broad pilot**

```bash
JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" \
  julia --project=. examples/scan_toric3d_x_feedback_threshold.jl \
  --sizes 3,5,7,9 \
  --error-rates 0,.0001,.0003,.001,.002,.003,.005,.0075,.01,.015,.02 \
  --shots 1000 --seed 1234 --routing one_step \
  --output-dir results/toric3d-x-feedback-threshold-2026-10-11/pilot
```

Expected: 44 complete paired points, timing pilot stored separately, zero failures at `p=0`, and all recovery validations pass.

- [x] **Step 2: Select and record the production window**

Compute adjacent-size crossings from the pilot. Choose at least seven probabilities spanning the union of the central adjacent-size crossing intervals, expanded by one pilot grid interval on both sides. Use at least five odd sizes beginning with `L=5,7,9,11,13`; add `L=15` or larger if the two largest-size crossings still drift by more than their Monte Carlo interval. If no crossing exists, extend `p` in the indicated direction and rerun the pilot before fitting. Save the selection and numerical rationale in `production/scan-selection.toml`.

- [x] **Step 3: Time and run 10,000-shot one-step production**

Run the CLI with the selected sizes/probabilities, `--shots 10000`, seed 1234, and `--routing one_step`. Use `--resume` after interruption. Inspect each completed point for monotonicity, invalid-recovery absence, and sufficient critical-window coverage before stopping the size sequence.

Execution deviation: the literal 50-point, 10,000-shot design timed at 143.8
single-core hours, with 96.7% of a profile in final BP+OSD and no batch API in
the installed backend.  The completed production grid therefore retains five
sizes and seven probabilities with 1,000 independent paired shots per point,
pooling the recorded 200-shot scout with 800 new shots for the three largest
sizes.  `production/scan-selection.toml` records every component seed and the
35,000-shot total.

- [x] **Step 4: Run full-path sensitivity**

Time 200 shots at the largest production size and central crossing. If the estimate is practical, run 10,000 shots on the production grid; otherwise choose the smallest shot count of 2,000, 5,000, or 10,000 whose predicted 95% Monte Carlo half-width is below one quarter of the observed one-step/full-path rate difference. Record the rule and selected count in metadata.

Execution result: at `L=13,p=.015`, 200 full-path shots gave 0.200 state
failure versus 0.226 for the 1,000-shot one-step point.  Even 10,000 full-path
shots have an estimated 0.00784 half-width, above one quarter of the observed
0.026 difference (0.0065), so none of the enumerated counts satisfies the
rule.  The 200-shot point is retained only as an exploratory sensitivity.

- [x] **Step 5: Fit, bootstrap, and audit the results**

Run the Task 6 fitter with 1,000 bootstrap replicas and seed 1234. Start with all production sizes and the selected critical window, then apply its minimum-size/window/correction sensitivities. If the primary fit fails a guard, report the failure and widen the scan rather than overriding it.

Execution result: the window already contains both tails and every adjacent
feedback pair crosses.  The all-size quadratic fit fails strongly and the
all-size cubic alternative narrowly fails the fixed goodness guard.  These
failures are retained.  The least aggressive planned size sensitivity,
cubic with `L>=7`, passes all guards and 999/1000 bootstrap refits; it is
reported as a conditional estimate rather than relabeling the rejected
all-size model as accepted.

- [x] **Step 6: Inspect every final figure and verify hashes**

Open the PNG and PDF figures, check labels/error bars/legends, recompute input and script SHA-256 values, and require them to match the JSON report. Confirm the summed unique paired shots equal the metadata total; do not multiply by metric rows.

- [x] **Step 7: Commit reproducible numerical artifacts**

```bash
git add results/toric3d-x-feedback-threshold-2026-10-11
git commit -m "data: run diagonal X feedback threshold scan"
```

### Task 8: Document the result and run full verification

**Files:**
- Create: `docs/superpowers/toric3d-x-feedback-threshold.md`
- Modify: `README.md`
- Modify: `document/note.tex`
- Test: `test/runtests.jl`

**Interfaces:**
- Consumes: verified Task 7 artifacts and fit report.
- Produces: scientific interpretation, reproduction commands, updated note, and final regression evidence.

- [x] **Step 1: Write the experiment report from machine-readable results**

State the circuit, logical `|00>` preparation, X fault locations, paired protocols, one-step and full-path feedback, final BP+OSD settings, scan grid, timing, raw rates, fitted threshold or guarded no-fit conclusion, Wilson/bootstrap intervals, model sensitivity, paper relationship, and all limitations. Link every CSV/TOML/JSON/PDF/SVG/PNG artifact.

- [x] **Step 2: Update README and the 3D section of the note**

Add the public API/CLI and a concise result. Explicitly distinguish logical-X `|00>` preparation failure from the prior logical-Z diagnostic, explain why feedback acts only in the X sector, and state that global classical MWPM is not a constant-depth local channel.

- [x] **Step 3: Run focused Julia and Python suites**

```bash
JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" \
  julia --project=. -e 'using Test, TopoNoise; include("test/toric3d/runtests.jl")'
PYTHONDONTWRITEBYTECODE=1 .CondaPkg/.pixi/envs/default/bin/python \
  test/toric3d/test_x_feedback_fit.py
```

Expected: all focused tests pass.

- [x] **Step 4: Run the complete regression suite**

```bash
JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" \
  julia --project=. test/runtests.jl
```

Expected: exit code 0, including all existing 2D and 3D tests.

- [x] **Step 5: Compile and inspect the note**

```bash
cd document
latexmk -pdf -interaction=nonstopmode -halt-on-error \
  -outdir=/tmp/toric3d-x-feedback-note-check note.tex
```

Expected: PDF generated successfully; no new undefined references or overfull boxes in the edited section.

- [x] **Step 6: Run repository and artifact checks**

```bash
git diff --check
find examples test/toric3d -type d -name __pycache__ -print
```

Expected: no whitespace errors and no generated Python cache directories.

- [x] **Step 7: Commit documentation and final verification updates**

```bash
git add docs/superpowers/toric3d-x-feedback-threshold.md README.md
git commit -m "docs: report diagonal X feedback threshold"
```

`document/note.tex` is intentionally ignored by the repository's `*.tex`
rule; leave its verified workspace update in place rather than forcing an
ignored source into this commit.
