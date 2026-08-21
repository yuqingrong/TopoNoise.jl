# Rotated Planar Code-Capacity Threshold Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the invalid open-edge scan with a verified open rotated planar surface-code code-capacity experiment for data-X noise and perfect Z-check syndrome.

**Architecture:** Stim generates the open rotated-memory-Z stabilizer circuit and supplies detector events plus a held-out logical-observable truth bit. A data-X-only injection is inserted after initialization, PyMatching decodes only the detector events, and the existing finite-size crossing machinery bootstraps equal-sized batches. The virtual-trajectory path is restored to its committed interface and remains independent.

**Tech Stack:** Julia 1.10, PythonCall, CondaPkg, Stim, PyMatching, Random, Statistics, CairoMakie, Test.

## Global Constraints

- Geometry is an open rotated planar surface-code patch, never a periodic torus or PEPS virtual-bond array.
- Noise is independent `X_ERROR(p)` on data qubits only, inserted after preparation and before noiseless extraction; gates, ancillas, and measurements are perfect.
- Scan the `surface_code:rotated_memory_z` task and report exactly one logical-X decoder-failure channel.
- The decoder receives detector syndrome only; the sampled logical observable is evaluation-only ground truth.
- Require `shots % batches == 0`, `2 <= batches <= shots`, positive strictly increasing distances/rates, and finite `0 <= p < 0.5`.
- Retain `CriticalCrossing` status semantics: `:ok`, `:unbracketed`, and `:unstable`.
- The existing virtual-trajectory sampler, its dirty files, and its public result type must not be changed or staged.
- Remove the invalid `OpenCodeCapacity*` public API and old output names; do not leave a path that can generate its false threshold.
- Every implementation task uses TDD, runs its focused test, passes `git diff --check`, and commits only its listed files.

---

## File structure

| File | Responsibility |
| --- | --- |
| `CondaPkg.toml` | Declares the Python `stim` runtime dependency beside PyMatching. |
| `src/code_capacity.jl` | Rotated-patch circuit construction, data-X injection, sampling, scan data, and crossing estimation. |
| `src/decoder.jl` | Retains virtual-trajectory MWPM; removes the obsolete manual open-code decoder during migration. |
| `src/TopoNoise.jl` | Exports only the supported rotated code-capacity API and includes files in their dependency order. |
| `src/visualization.jl` | Restores the legacy trajectory plot and adds the isolated rotated-code plot. |
| `examples/scan_rotated_code_capacity.jl` | Validated CLI, CSV writing, figure export, and textual crossing summary. |
| `test/code_capacity.jl` | Circuit, distance, decoding, scan, validation, and crossing-status tests. |
| `test/code_capacity_example.jl` | CLI output and artifact smoke test. |
| `test/visualization.jl` | Rendered rotated-plot and unchanged trajectory-plot regressions. |

### Shared interfaces introduced by this plan

```julia
struct RotatedCodeCapacityModel
    distance::Int
    base_circuit::Py
    data_qubits::Vector{Int}              # zero-based Stim qubit indices
end

struct RotatedCodeCapacityScanPoint
    distance::Int
    error_rate::Float64
    shots::Int
    logical_x_failure_count::Int
    logical_x_failure_rate::Float64
    logical_x_failure_se::Float64
    logical_x_failure_batches::Vector{Float64}
    seed::Union{Missing,Int}
end

struct RotatedCodeCapacityScan
    distances::Vector{Int}
    error_rates::Vector{Float64}
    points::Vector{RotatedCodeCapacityScanPoint}
end

RotatedCodeCapacityModel(distance::Integer)
rotated_code_capacity_circuit(model::RotatedCodeCapacityModel,
                              error_rate::Real)::Py
estimate_rotated_code_capacity(rng::Random.AbstractRNG,
                               model::RotatedCodeCapacityModel,
                               error_rate::Real;
                               shots::Integer, batches::Integer,
                               seed::Union{Nothing,Integer}=nothing)
scan_rotated_code_capacity(rng::Random.AbstractRNG, distances, error_rates;
                           shots::Integer=10_000, batches::Integer=100,
                           seed::Union{Nothing,Integer}=nothing,
                           progress_io=nothing)::RotatedCodeCapacityScan
estimate_rotated_code_crossings(rng::Random.AbstractRNG,
                                scan::RotatedCodeCapacityScan;
                                bootstrap::Integer=2_000,
                                confidence::Real=0.95)::Vector{CriticalCrossing}
plot_rotated_code_capacity(scan::RotatedCodeCapacityScan,
                           crossings::AbstractVector{<:CriticalCrossing}=CriticalCrossing[])
```

The private helpers `_stim()`, `_data_qubits_from_terminal_measurement`,
`_circuit_with_data_x_noise`, `_deterministic_data_x_circuit`,
`_rotated_code_matching`, `_decode_observables`,
`_rotated_code_capacity_scan_point`, and `_minimum_x_logical_support` are
defined in `src/code_capacity.jl`.  `decode_observables` deliberately has no
argument for actual logical-observable bits.

### Task 1: Add the Stim-backed rotated-patch circuit model

**Files:**
- Modify: `CondaPkg.toml`
- Modify: `src/TopoNoise.jl:12-32`
- Modify: `src/code_capacity.jl:1-4` and append rotated-model helpers
- Modify: `test/code_capacity.jl`

**Consumes:** PythonCall `Py`, `pyimport`, and the current CondaPkg-managed Python environment.

**Produces:** `RotatedCodeCapacityModel`, `rotated_code_capacity_circuit`, the private deterministic-injection helper, and a reproducible model-validation surface for later decoder work.

- [ ] **Step 1: Write the failing circuit-model tests**

  Replace the old manual-edge construction tests with a focused rotated-patch section. It must assert distance validation, data-only injection, a noiseless zero result, and the generated logical distance.

  ```julia
  @testset "rotated planar code-capacity circuit" begin
      @test_throws ArgumentError RotatedCodeCapacityModel(0)
      model = RotatedCodeCapacityModel(3)
      @test model.distance == 3
      @test length(model.data_qubits) == 9

      noisy = rotated_code_capacity_circuit(model, 0.1)
      @test occursin("X_ERROR(0.1)", string(noisy))
      @test all(string(qubit) in string(noisy) for qubit in model.data_qubits)

      clean = rotated_code_capacity_circuit(model, 0.0)
      dets, observables = clean.compile_detector_sampler().sample(
          shots=8, separate_observables=true)
      @test !any(TopoNoise.pyconvert(BitMatrix, dets))
      @test !any(TopoNoise.pyconvert(BitMatrix, observables))
      @test length(noisy.shortest_graphlike_error()) == 3
  end
  ```

- [ ] **Step 2: Run the test to verify it fails**

  Run: `julia --project=. test/code_capacity.jl`

  Expected: FAIL because `RotatedCodeCapacityModel` is undefined and Stim is not yet declared.

- [ ] **Step 3: Add the dependency and model implementation**

  Add `stim = ""` under `[pip.deps]` in `CondaPkg.toml`. Add the complete
  replacement public surface to the `TopoNoise` export list without removing
  the old API yet: `RotatedCodeCapacityModel`, `RotatedCodeCapacityScanPoint`,
  `RotatedCodeCapacityScan`, `rotated_code_capacity_circuit`,
  `estimate_rotated_code_capacity`, `scan_rotated_code_capacity`,
  `estimate_rotated_code_crossings`, and `plot_rotated_code_capacity`.
  Julia permits these forward exports while later tasks define the names; this
  lets the Task 4 example import the API normally.

  In `src/code_capacity.jl`, cache `pyimport("stim")` in a `Ref{Py}`. Construct the base circuit with no noise:

  ```julia
  const _stim_ref = Ref{Py}()
  _stim() = isassigned(_stim_ref) ? _stim_ref[] :
      (_stim_ref[] = pyimport("stim"))

  function RotatedCodeCapacityModel(distance::Integer)
      distance >= 3 || throw(ArgumentError("distance must be at least 3"))
      base = _stim().Circuit.generated(
          "surface_code:rotated_memory_z"; distance=Int(distance), rounds=1)
      data = _data_qubits_from_terminal_measurement(base)
      length(data) == Int(distance)^2 || throw(ErrorException(
          "generated rotated patch did not expose d² data qubits"))
      return RotatedCodeCapacityModel(Int(distance), base, data)
  end
  ```

  Parse the terminal data-basis `M` instruction from Stim's circuit text; reject absent, duplicate, or non-`d^2` terminal data measurements. Rebuild a new Stim circuit by splitting immediately after the first preparation `TICK`, injecting `X_ERROR(p)` on exactly `data_qubits`, then appending the untouched remainder. The deterministic-support helper uses `X_ERROR(1)`, not a literal Clifford `X`: the physical Pauli is the same, but Stim must classify it as an error mechanism in order to report its detector and logical-observable flips. Validate `0 <= p < 0.5`; preserve the no-noise `p == 0` circuit without a zero-probability instruction.

  Create `_minimum_x_logical_support(model)` by grouping data qubits by each coordinate row and column, sampling each length-`d` deterministic X string, and selecting the unique candidate with zero detector events and logical observable `true`. Throw if there is not exactly one orientation family or if support length differs from `d`.

- [ ] **Step 4: Resolve CondaPkg and run the circuit-model tests**

  Run: `julia --project=. -e 'using CondaPkg; CondaPkg.resolve(); using TopoNoise'`

  Then run: `julia --project=. test/code_capacity.jl`

  Expected: PASS. The test confirms a d=3 patch has nine data qubits, p=0 produces no detection/logical events, injected noise targets only data qubits, and the only graphlike noise mechanism has shortest logical length three.

- [ ] **Step 5: Commit the circuit model**

  ```bash
  git add CondaPkg.toml src/TopoNoise.jl src/code_capacity.jl test/code_capacity.jl
  git diff --cached --check
  git commit -m "feat: add rotated planar code-capacity circuit"
  ```

### Task 2: Decode held-out logical observables and certify error correction

**Files:**
- Modify: `src/code_capacity.jl`
- Modify: `test/code_capacity.jl`

**Consumes:** `RotatedCodeCapacityModel`, `rotated_code_capacity_circuit`, the current `_pymatching()` loader in `src/decoder.jl`, and Stim detector sampling.

**Produces:** `estimate_rotated_code_capacity`, private detector-only decoder helpers, and physics validation tests for single errors and a weight-`d` logical string.

- [ ] **Step 1: Write the failing decoding tests**

  Add tests that use deterministic X supports so the expected syndrome and logical truth are known without Monte Carlo.

  ```julia
  @testset "perfect Z-syndrome MWPM decoding" begin
      model = RotatedCodeCapacityModel(3)
      matching = TopoNoise._rotated_code_matching(model, 0.1)

      for qubit in model.data_qubits
          circuit = TopoNoise._deterministic_data_x_circuit(model, [qubit])
          syndrome_py, actual_py = circuit.compile_detector_sampler().sample(
              shots=1, separate_observables=true)
          syndrome = TopoNoise.pyconvert(BitMatrix, syndrome_py)
          actual = TopoNoise.pyconvert(BitMatrix, actual_py)
          predicted = TopoNoise.pyconvert(
              BitMatrix, TopoNoise._decode_observables(matching, syndrome_py))
          @test !any(predicted .!= actual)
      end

      support = TopoNoise._minimum_x_logical_support(model)
      circuit = TopoNoise._deterministic_data_x_circuit(model, support)
      syndrome_py, actual_py = circuit.compile_detector_sampler().sample(
          shots=1, separate_observables=true)
      syndrome = TopoNoise.pyconvert(BitMatrix, syndrome_py)
      actual = TopoNoise.pyconvert(BitMatrix, actual_py)
      @test !any(syndrome)
      @test only(vec(actual))
      @test_throws MethodError TopoNoise._decode_observables(
          matching, syndrome_py, actual_py)
  end

  @testset "single-point estimator" begin
      point = estimate_rotated_code_capacity(
          MersenneTwister(17), RotatedCodeCapacityModel(3), 0.0;
          shots=12, batches=3, seed=17)
      @test point.logical_x_failure_count == 0
      @test point.logical_x_failure_rate == 0.0
      @test length(point.logical_x_failure_batches) == 3
  end
  ```

- [ ] **Step 2: Run the test to verify it fails**

  Run: `julia --project=. test/code_capacity.jl`

  Expected: FAIL because matching and single-point estimator helpers are undefined.

- [ ] **Step 3: Implement detector-only matching and estimation**

  Build a matching graph only from `rotated_code_capacity_circuit(model, p)`:

  ```julia
  function _rotated_code_matching(model, p::Float64)
      p > 0 || throw(ArgumentError("matching requires p > 0"))
      dem = rotated_code_capacity_circuit(model, p).detector_error_model(
          decompose_errors=true,
          block_decomposition_from_introducing_remnant_edges=true)
      return _pymatching().Matching.from_detector_error_model(dem)
  end

  _decode_observables(matching::Py, syndrome) = matching.decode_batch(syndrome)
  ```

  Convert detector and observable matrices to Julia `BitMatrix` only after
  `decode_batch`; do not convert or pass actual observables to the decoder.
  For `p == 0`, sample the clean circuit and assert no detector/logical flips,
  then return a zero-count point without constructing a degenerate matching
  graph. For `p > 0`, use a Stim detector sampler seeded from a `UInt64` drawn
  from the supplied Julia RNG, decode the whole batch, score `any(predicted .!=
  actual; dims=2)`, and aggregate exactly equal contiguous batches. Calculate
  `sqrt(rate * (1-rate) / shots)` with zero at the endpoints.

  Validate positive shots, `2 <= batches <= shots`, divisibility, valid p, and
  a supplied `seed` convertible to `Int`; copy the supplied seed into the
  result row. Add a model-level detector-count check before decoding so a
  malformed circuit cannot silently reach PyMatching.

- [ ] **Step 4: Run the focused tests**

  Run: `julia --project=. test/code_capacity.jl`

  Expected: PASS. Every weight-one data-X error is decoded correctly, the
  selected weight-three boundary string has zero syndrome but flips logical
  truth, and no callable decoder overload accepts that truth bit.

- [ ] **Step 5: Commit the decoder path**

  ```bash
  git add src/code_capacity.jl test/code_capacity.jl
  git diff --cached --check
  git commit -m "feat: decode rotated code-capacity syndromes"
  ```

### Task 3: Add equal-batch finite-size scans and crossing estimates

**Files:**
- Modify: `src/code_capacity.jl`
- Modify: `test/code_capacity.jl`

**Consumes:** `estimate_rotated_code_capacity`, `_validate_strictly_increasing`, `_isotonic_non_decreasing`, `_selected_crossing`, `_bootstrap_batch_mean`, and `CriticalCrossing` from `src/trajectory_analysis.jl`.

**Produces:** `RotatedCodeCapacityScanPoint`, `RotatedCodeCapacityScan`, `scan_rotated_code_capacity`, and `estimate_rotated_code_crossings`.

- [ ] **Step 1: Write failing scan and status tests**

  ```julia
  @testset "rotated finite-size scan" begin
      scan = scan_rotated_code_capacity(
          MersenneTwister(21), [3, 5], [0.0, 0.1];
          shots=12, batches=3, seed=21)
      @test length(scan.points) == 4
      @test all(point -> point.seed == 21, scan.points)
      @test all(point -> length(point.logical_x_failure_batches) == 3, scan.points)
      @test_throws ArgumentError scan_rotated_code_capacity(
          MersenneTwister(1), [3, 5], [0.1, 0.2]; shots=10, batches=3)
  end

  @testset "rotated crossing statuses" begin
      @test only(estimate_rotated_code_crossings(
          MersenneTwister(2), synthetic_rotated_scan(
              [0.08, 0.10, 0.12], [0.2, 0.45, 0.7], [0.3, 0.45, 0.6]);
          bootstrap=20)).status == :ok
      @test only(estimate_rotated_code_crossings(
          MersenneTwister(3), synthetic_rotated_scan(
              [0.08, 0.10, 0.12], [0.2, 0.3, 0.4], [0.1, 0.2, 0.3]);
          bootstrap=20)).status == :unbracketed
      @test only(estimate_rotated_code_crossings(
          MersenneTwister(4), synthetic_rotated_scan(
              [0.08, 0.10, 0.12], [0.2, 0.5, 0.8], [0.2, 0.5, 0.8]);
          bootstrap=20)).status == :unstable
  end
  ```

  `synthetic_rotated_scan` must construct the shared scan-point fields with
  four equal batches per point, so test data does not bypass bootstrap shape
  checks.

- [ ] **Step 2: Run the test to verify it fails**

  Run: `julia --project=. test/code_capacity.jl`

  Expected: FAIL because scan types and functions are undefined.

- [ ] **Step 3: Implement scan and crossing functions**

  Define the two scan structs exactly as in the shared interfaces. Validate
  all distances/rates before doing any sampling, cache one
  `RotatedCodeCapacityModel` per distance inside each scan, call the
  single-point estimator for every grid point, and print
  `completed d=<distance> p=<rate> (<shots> shots)` only after each successful
  point.

  Implement `_rotated_code_capacity_scan_point(scan, distance, error_rate)`
  with the same uniqueness guard used by trajectory scans. Implement crossings
  from `logical_x_failure_rate` and `logical_x_failure_batches`: monotone
  curves, `_selected_crossing`, bootstrap resampling, `0.8` minimum valid
  bootstrap fraction, and the exact `CriticalCrossing` field order. Do not
  accept unequal batch lengths or incomplete grids.

- [ ] **Step 4: Run focused scan tests**

  Run: `julia --project=. test/code_capacity.jl`

  Expected: PASS. The tiny scan returns four points, invalid unequal batch
  divisions fail before sampling, and all three crossing states are covered.

- [ ] **Step 5: Commit scan statistics**

  ```bash
  git add src/code_capacity.jl test/code_capacity.jl
  git diff --cached --check
  git commit -m "feat: scan rotated code-capacity crossings"
  ```

### Task 4: Replace the CLI and render an honest crossing summary

**Files:**
- Create: `examples/scan_rotated_code_capacity.jl`
- Modify: `src/visualization.jl:232-500`
- Modify: `test/code_capacity_example.jl`
- Modify: `test/visualization.jl:14-49`
- Delete: `examples/scan_open_code_capacity.jl`

**Consumes:** `RotatedCodeCapacityScan`, `estimate_rotated_code_crossings`, `plot_rotated_code_capacity`, and `CriticalCrossing`.

**Produces:** the isolated rotated-code figure and CLI output files:
`rotated_code_capacity_scan.csv`, `rotated_code_capacity_crossings.csv`, and
`rotated_code_capacity_scan.{svg,pdf,png}`.

- [ ] **Step 1: Write failing CLI and visualization tests**

  Change the example test to include the new module and assert every promised
  artifact and its model label. Add a deterministic mixed-status figure test.

  ```julia
  status = ScanRotatedCodeCapacityExample.main([
      "--p-min", "0.05", "--p-max", "0.10", "--p-step", "0.05",
      "--sizes", "3,5", "--shots", "12", "--batches", "3",
      "--bootstrap", "12", "--seed", "31", "--output-dir", directory]; io=output)
  @test status == 0
  @test all(isfile(joinpath(directory, name)) for name in (
      "rotated_code_capacity_scan.csv",
      "rotated_code_capacity_crossings.csv",
      "rotated_code_capacity_scan.svg",
      "rotated_code_capacity_scan.pdf",
      "rotated_code_capacity_scan.png"))
  @test occursin("open rotated planar code", read(
      joinpath(directory, "rotated_code_capacity_scan.svg"), String))
  @test occursin("unbracketed", String(take!(output)))
  @test startswith(read(joinpath(directory, "rotated_code_capacity_scan.svg"), String), "<?xml")
  @test startswith(read(joinpath(directory, "rotated_code_capacity_scan.pdf"), String), "%PDF")
  @test read(joinpath(directory, "rotated_code_capacity_scan.png"))[1:8] ==
      UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]
  ```

  In `test/visualization.jl`, build a synthetic scan with one `:ok` and one
  `:unbracketed` crossing, render it, and assert `rendered_text(figure)`
  includes both `p_c` and `unbracketed`. Keep the existing trajectory test
  and make it assert that no "Logical failure (N-S, MWPM)" panel exists.

- [ ] **Step 2: Run tests to verify they fail**

  Run: `julia --project=. test/code_capacity_example.jl`

  Run: `julia --project=. test/visualization.jl`

  Expected: FAIL because the new example, new plot, and status annotations do
  not yet exist; the legacy trajectory test may also expose its current
  incompatible logical-field access.

- [ ] **Step 3: Implement the CLI and two-panel plot**

  Copy the validated option parsing pattern from the old example, changing all
  user-facing wording to “data-X error rate”, “rotated planar distance”, and
  “perfect Z syndrome”. Require increasing distances of at least three and
  `shots % batches == 0` in CLI validation. Pass `seed=options.seed` into the
  scan. Write these exact CSV columns:

  ```text
  L,p,shots,logical_x_failure_count,logical_x_failure_rate,logical_x_failure_se,batches,seed
  small_size,large_size,estimate,ci_low,ci_high,confidence,valid_bootstrap_fraction,status
  ```

  `plot_rotated_code_capacity` creates a curve axis titled
  “Logical X failure — open rotated planar code” with subtitle text
  “data-X noise; perfect Z syndrome”, plus a summary axis. The summary prints
  for every adjacent pair either
  `L=a/b: p_c = e [lo, hi] (95%)` or
  `L=a/b: unbracketed — extend p range` / `unstable — increase samples`.
  Thus an invalid crossing cannot appear as an unlabeled blank plot.

  Restore `plot_trajectory_scan` exactly to the signature and body from
  `ea21dc2` (only raw, marginal-spin, and Binder panels). Do not read logical
  fields absent from committed `TrajectoryScanPoint`, and do not stage any
  dirty trajectory file.

- [ ] **Step 4: Run all focused render and CLI tests**

  Run: `julia --project=. test/code_capacity_example.jl`

  Run: `julia --project=. test/visualization.jl`

  Expected: PASS. SVG contains the model and crossing-status text; CSV, SVG,
  PDF, and PNG are nonempty; legacy trajectory visualization renders without
  logical-field errors.

- [ ] **Step 5: Commit the CLI and plotting layer**

  ```bash
  git add examples/scan_rotated_code_capacity.jl src/visualization.jl \
      test/code_capacity_example.jl test/visualization.jl
  git rm examples/scan_open_code_capacity.jl
  git diff --cached --check
  git commit -m "feat: add rotated code-capacity scan example"
  ```

### Task 5: Remove the invalid manual open-code path and publish only the replacement API

**Files:**
- Modify: `src/code_capacity.jl`
- Modify: `src/decoder.jl:1-310`
- Modify: `src/TopoNoise.jl:20-32`
- Modify: `src/visualization.jl` (remove only `plot_open_code_capacity`)
- Modify: `test/code_capacity.jl`
- Modify: `test/code_capacity_example.jl`
- Modify: `test/visualization.jl`

**Consumes:** the completed rotated model, decoder, scan, CLI, and plotting API.

**Produces:** a clean module with no `OpenCodeCapacity`, `DataEdgeErrors`,
`scan_open_code_capacity`, `estimate_open_code_crossings`,
`plot_open_code_capacity`, or `decode_syndrome` code-capacity path.

- [ ] **Step 1: Write the failing migration/isolation tests**

  ```julia
  @testset "invalid open-code API is removed" begin
      @test !isdefined(TopoNoise, :OpenCodeCapacityModel)
      @test !isdefined(TopoNoise, :scan_open_code_capacity)
      @test !isdefined(TopoNoise, :plot_open_code_capacity)
  end
  ```

  Retain the `decode_uf` virtual-trajectory tests in `test/decoder.jl`; add a
  short invocation of `plot_trajectory_scan` to the visualization test as the
  regression proof that no trajectory interface changed.

- [ ] **Step 2: Run the migration test to verify it fails**

  Run: `julia --project=. test/code_capacity.jl`

  Expected: FAIL because old manual public symbols are still exported.

- [ ] **Step 3: Remove the false geometry without touching virtual trajectories**

  Delete `OpenCodeCapacityModel`, `DataEdgeErrors`, manual plaquette syndrome,
  manual logical cuts, old scan types/functions, and their tests from
  `src/code_capacity.jl`. Delete only the `OpenCodeCapacityModel`-typed
  matching/cache/`decode_syndrome` functions from `src/decoder.jl`; retain
  `Correction`, `plaquette_syndrome`, `_matching`, `decode_uf`,
  `logical_failure(::VirtualBondErrors, ...)`, and `decode_trajectory`.

  Remove `plot_open_code_capacity` from `src/visualization.jl`. Remove all obsolete exports while preserving the complete replacement list
  added in Task 1:

  ```julia
  RotatedCodeCapacityModel, RotatedCodeCapacityScanPoint,
  RotatedCodeCapacityScan, rotated_code_capacity_circuit,
  estimate_rotated_code_capacity, scan_rotated_code_capacity,
  estimate_rotated_code_crossings, plot_rotated_code_capacity
  ```

  Use `rg -n "OpenCodeCapacity|DataEdgeErrors|scan_open_code|plot_open_code|decode_syndrome" src test examples` to verify no stale API remains. Do not create a compatibility alias that could run the invalid model.

- [ ] **Step 4: Run migration and regression tests**

  Run: `julia --project=. test/code_capacity.jl`

  Run: `julia --project=. test/decoder.jl`

  Run: `julia --project=. test/trajectory_analysis.jl`

  Run: `julia --project=. test/visualization.jl`

  Expected: PASS. The old API is absent, the new API is exported, and virtual
  decoder/trajectory behavior still works.

- [ ] **Step 5: Commit the API migration**

  ```bash
  git add src/code_capacity.jl src/decoder.jl src/TopoNoise.jl \
      test/code_capacity.jl test/code_capacity_example.jl test/visualization.jl
  git diff --cached --check
  git commit -m "refactor: replace invalid open code-capacity model"
  ```

### Task 6: Verify a clean tree and produce a qualified threshold scan

**Files:**
- Modify: `test/code_capacity.jl` only if a failed full-suite check exposes a test defect in the new implementation
- Generated, untracked: `examples/output/rotated_code_capacity_*`

**Consumes:** all completed implementation tasks.

**Produces:** verified test evidence and a finite-size threshold output that is reported only when its crossing status is `:ok`.

- [ ] **Step 1: Run clean focused verification**

  Run:

  ```bash
  julia --project=. test/code_capacity.jl
  julia --project=. test/code_capacity_example.jl
  julia --project=. test/decoder.jl
  julia --project=. test/visualization.jl
  julia --project=. -e 'using Pkg; Pkg.test()'
  ```

  Expected: every command exits zero. If the package suite fails, fix only the
  failing new-path test or implementation code, rerun the complete list, and
  commit that exact correction before proceeding.

- [ ] **Step 2: Run the reproducible scientific scan**

  Run:

  ```bash
  julia --project=. examples/scan_rotated_code_capacity.jl \
    --p-min 0.075 --p-max 0.125 --p-step 0.005 \
    --sizes 3,5,7,9 --shots 100000 --batches 100 --bootstrap 2000 \
    --seed 20260811 --output-dir examples/output
  ```

  Verify the five named output files are nonempty, that SVG/PDF/PNG are valid,
  and that every CSV row identifies data-X/perfect-Z model parameters.

- [ ] **Step 3: Apply the reporting gate**

  Read `rotated_code_capacity_crossings.csv`. Report a numerical finite-size
  `p_c` only for rows with `status=ok`, finite `estimate`, `ci_low`, and
  `ci_high`. If all rows are `unbracketed`, report that the scanned interval
  does not bracket a crossing; if any selected row is `unstable`, report that
  the bootstrap is unstable. Do not convert either condition into a threshold
  number.

- [ ] **Step 4: Commit only any verification-driven source correction**

  ```bash
  git status --short
  git diff --check
  ```

  Generated scan files remain untracked/ignored. If source or test changes
  were needed in Step 1, stage exactly those files and commit them with
  `fix: verify rotated code-capacity scan`; otherwise make no empty commit.

## Plan self-review

- **Spec coverage:** Tasks 1–2 implement data-X-only perfect extraction and distance/decoder isolation; Task 3 implements equal-batch statistics and all statuses; Task 4 implements artifacts and visible status; Task 5 removes the invalid model and repairs trajectory isolation; Task 6 runs the clean suite and applies the numerical reporting gate.
- **Placeholder scan:** no `TODO`, `TBD`, deferred implementation, or unspecified validation step remains.
- **Type consistency:** all later tasks use the Task 1 model and Task 3 scan types, `CriticalCrossing` retains its existing fields, and the sole public plotting function is `plot_rotated_code_capacity`.
