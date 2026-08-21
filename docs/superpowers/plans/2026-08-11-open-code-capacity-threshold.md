# Open-Boundary Code-Capacity Threshold Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an independent open-boundary code-capacity sampler, syndrome-only MWPM decoder workflow, finite-size crossing scan, and reproducible command-line figure without altering the virtual-bond trajectory scan.

**Architecture:** A new `code_capacity.jl` module owns data-edge errors, syndrome generation, logical-cut geometry, streaming statistics, and crossing estimation. `decoder.jl` gains a syndrome-only entry point which shares the Boolean `Correction` representation but cannot receive sampled errors. A dedicated example and two-panel CairoMakie plot consume only the new scan type; trajectory APIs and their output names stay intact.

**Tech Stack:** Julia, Random, Statistics, PythonCall/PyMatching, CairoMakie, Test.

## Global Constraints

- Retain the existing open rectangular `R × C` bond geometry; do not add periodic edge identification.
- Data-edge errors are independent code-capacity `X` errors and must never be taken from a PEPS trajectory or virtual mismatch record.
- The decoder accepts only syndrome plus model/sector metadata; true errors are available only after decoding for the verdict.
- North--south is the primary crossing; east--west is a symmetry validation series.
- Existing `scan_toric_trajectories.jl`, `plot_trajectory_scan`, and their outputs are behaviorally unchanged.
- Treat an unbracketed curve pair as `:unbracketed`, never as an endpoint estimate.

---

### Task 1: Isolate data-edge sampling and open-patch geometry

**Files:**

- Create: `src/code_capacity.jl`
- Modify: `src/TopoNoise.jl`
- Modify: `test/runtests.jl`
- Create: `test/code_capacity.jl`

**Interfaces:**

- Produces `OpenCodeCapacityModel(rows::Integer, cols::Integer)` with positive dimensions at least 2.
- Produces `DataEdgeErrors(horizontal::BitMatrix, vertical::BitMatrix)` with shapes `(R, C - 1)` and `(R - 1, C)`.
- Produces `sample_data_edge_errors(rng, model; error_rate::Real)::DataEdgeErrors` and `code_capacity_syndrome(errors)::BitMatrix`.
- Produces `logical_cut(model; sector::Symbol)::NamedTuple` with Boolean masks named `horizontal` and `vertical`.

- [ ] **Step 1: Write failing model, sampler, and incidence tests**

Add `test/code_capacity.jl`, include it from `test/runtests.jl`, and start with deterministic checks:

```julia
using TopoNoise
using Random
using Test

@testset "open code-capacity data edges" begin
    model = OpenCodeCapacityModel(4, 5)
    zero = sample_data_edge_errors(MersenneTwister(1), model; error_rate=0.0)
    one = sample_data_edge_errors(MersenneTwister(1), model; error_rate=1.0)
    @test size(zero.horizontal) == (4, 4)
    @test size(zero.vertical) == (3, 5)
    @test !any(zero.horizontal) && !any(zero.vertical)
    @test all(one.horizontal) && all(one.vertical)
    @test !any(code_capacity_syndrome(zero))
end

@testset "one edge has the plaquette-incidence syndrome" begin
    errors = DataEdgeErrors(falses(4, 3), falses(3, 4))
    errors.horizontal[2, 2] = true
    syndrome = code_capacity_syndrome(errors)
    @test syndrome[1, 2] && syndrome[2, 2]
    @test count(syndrome) == 2
end
```

- [ ] **Step 2: Run the focused test to verify it fails**

Run: `julia --project=. test/code_capacity.jl`

Expected: FAIL because `OpenCodeCapacityModel`, `DataEdgeErrors`, and the sampler are undefined.

- [ ] **Step 3: Implement the standalone model and sampler**

Create `src/code_capacity.jl` before `decoder.jl` in `src/TopoNoise.jl`. Implement validated model/error constructors, independent Bernoulli sampling, and an explicit code-capacity syndrome method. Reuse only the four-edge parity formula, not `bond_mismatches` or any `ToricCodeTrajectory` field:

```julia
struct OpenCodeCapacityModel
    rows::Int
    cols::Int
end

struct DataEdgeErrors
    horizontal::BitMatrix
    vertical::BitMatrix
end

function code_capacity_syndrome(errors::DataEdgeErrors)
    R, C = size(errors.horizontal, 1), size(errors.vertical, 2)
    syndrome = falses(R - 1, C - 1)
    for r in 1:(R - 1), c in 1:(C - 1)
        syndrome[r, c] = xor(errors.horizontal[r, c], errors.horizontal[r + 1, c],
                              errors.vertical[r, c], errors.vertical[r, c + 1])
    end
    syndrome
end
```

Add `logical_cut` with explicit central transverse masks: for `:north_south`, set `horizontal[cld(rows, 2), :]` true; for `:east_west`, set `vertical[:, cld(cols, 2)]` true. All other mask entries are false. A north--south zero-syndrome chain is `horizontal[:, c]`, while an east--west one is `vertical[r, :]`; each therefore crosses its corresponding mask exactly once. Validate `sector in (:north_south, :east_west)`, and export the four public data-edge functions/types. Do not add a method that converts a `VirtualBondErrors` or a trajectory into `DataEdgeErrors`.

- [ ] **Step 4: Run the focused tests to verify they pass**

Run: `julia --project=. test/code_capacity.jl`

Expected: PASS for zero/full sampling and the single-edge incidence test.

- [ ] **Step 5: Commit the isolated data-edge model**

```bash
git add src/code_capacity.jl src/TopoNoise.jl test/runtests.jl test/code_capacity.jl
git commit -m "feat: add open code-capacity data model"
```

### Task 2: Decode syndrome only and verify logical sectors

**Files:**

- Modify: `src/decoder.jl`
- Modify: `src/code_capacity.jl`
- Modify: `test/code_capacity.jl`
- Modify: `test/decoder.jl`

**Interfaces:**

- Consumes `OpenCodeCapacityModel`, `DataEdgeErrors`, `code_capacity_syndrome`, and existing `Correction`.
- Produces `decode_syndrome(model::OpenCodeCapacityModel, syndrome::BitMatrix; sector::Symbol=:north_south, error_rate::Real=0.1)::Correction`.
- Produces `residual_errors(errors::DataEdgeErrors, correction::Correction)::DataEdgeErrors` and `logical_failure(errors::DataEdgeErrors, correction::Correction; sector::Symbol=:north_south)::Bool`.

- [ ] **Step 1: Write failing syndrome-only decoder tests**

Extend `test/code_capacity.jl` with a test that calls the decoder with a syndrome, checks that it restores that syndrome, and exhaustively checks every single data edge on a small patch:

```julia
function single_edge_error_configurations(model)
    errors = DataEdgeErrors[]
    for r in 1:model.rows, c in 1:(model.cols - 1)
        horizontal = falses(model.rows, model.cols - 1)
        horizontal[r, c] = true
        push!(errors, DataEdgeErrors(horizontal, falses(model.rows - 1, model.cols)))
    end
    for r in 1:(model.rows - 1), c in 1:model.cols
        vertical = falses(model.rows - 1, model.cols)
        vertical[r, c] = true
        push!(errors, DataEdgeErrors(falses(model.rows, model.cols - 1), vertical))
    end
    errors
end

@testset "syndrome-only decode" begin
    model = OpenCodeCapacityModel(4, 4)
    for sector in (:north_south, :east_west)
        for errors in single_edge_error_configurations(model)
            syndrome = code_capacity_syndrome(errors)
            correction = decode_syndrome(model, syndrome; sector=sector)
            residual = residual_errors(errors, correction)
            @test !any(code_capacity_syndrome(residual))
            @test !logical_failure(errors, correction; sector=sector)
        end
    end
end
```

Add the explicit zero-syndrome boundary-chain fixtures below. They make the physical geometry of each logical sector reviewable:

```julia
@testset "logical-cut geometry" begin
    model = OpenCodeCapacityModel(4, 4)
    none = Correction(4, 4)
    north_south = DataEdgeErrors(falses(4, 3), falses(3, 4))
    north_south.horizontal[:, 2] .= true
    @test !any(code_capacity_syndrome(north_south))
    @test logical_failure(north_south, none; sector=:north_south)
    @test !logical_failure(north_south, none; sector=:east_west)

    east_west = DataEdgeErrors(falses(4, 3), falses(3, 4))
    east_west.vertical[2, :] .= true
    @test !any(code_capacity_syndrome(east_west))
    @test logical_failure(east_west, none; sector=:east_west)
    @test !logical_failure(east_west, none; sector=:north_south)
end
```

- [ ] **Step 2: Run the focused decoder test to verify it fails**

Run: `julia --project=. test/code_capacity.jl`

Expected: FAIL because `decode_syndrome`, `residual_errors`, and the data-edge `logical_failure` method do not exist.

- [ ] **Step 3: Implement a decoder entry point without true-error access**

Refactor `_build_matching` only enough to expose a cached open-patch matching graph keyed by `(rows, cols, sector)`. Add a method which flattens the supplied `BitMatrix` syndrome in row-major plaquette order, calls PyMatching, and reconstructs a `Correction` from its data-edge fault IDs:

```julia
function decode_syndrome(
        model::OpenCodeCapacityModel, syndrome::BitMatrix;
        sector::Symbol=:north_south, error_rate::Real=0.1)::Correction
    _validate_open_syndrome(model, syndrome)
    0 <= error_rate < 0.5 || throw(ArgumentError("error_rate must be in [0, 0.5)"))
    any(syndrome) || return Correction(model.rows, model.cols)
    error_rate > 0 || throw(ArgumentError("nonzero syndrome is impossible at p=0"))
    matching = _open_code_matching(model, sector; error_rate=Float64(error_rate))
    predicted = matching.decode(_syndrome_vector(syndrome))
    return _correction_from_faults(model, predicted)
end
```

Make the graph's boundary edges and fault logical observable derive from the same `logical_cut` masks used for scoring. Preserve the old `decode_uf(ToricCodeTrajectoryModel, mismatches)` entry point unchanged for compatibility. Implement residual XOR and the data-edge overload of `logical_failure` in `code_capacity.jl`; only the scan harness calls it after decoding.

- [ ] **Step 4: Run decoder and trajectory-regression tests**

Run: `julia --project=. test/code_capacity.jl && julia --project=. test/decoder.jl`

Expected: PASS. Every single edge is corrected, each residual syndrome is zero, the two known logical strings give the declared parities, and the prior trajectory decoder tests still pass.

- [ ] **Step 5: Commit the syndrome-only decoder path**

```bash
git add src/decoder.jl src/code_capacity.jl test/code_capacity.jl test/decoder.jl
git commit -m "feat: decode open code-capacity syndromes"
```

### Task 3: Stream logical-failure scans and bootstrap crossings

**Files:**

- Modify: `src/code_capacity.jl`
- Modify: `test/code_capacity.jl`

**Interfaces:**

- Produces `OpenCodeCapacityScanPoint(size, error_rate, shots, logical_failure_ns_mean, logical_failure_ns_se, logical_failure_ew_mean, logical_failure_ew_se, logical_failure_ns_batches, logical_failure_ew_batches)`.
- Produces `OpenCodeCapacityScan(sizes, error_rates, points)`.
- Produces `scan_open_code_capacity(rng, sizes, error_rates; shots=10_000, batches=100, progress_io=nothing)::OpenCodeCapacityScan`.
- Produces `estimate_open_code_crossings(rng, scan; sector=:north_south, bootstrap=2_000, confidence=0.95)::Vector{CriticalCrossing}`.

- [ ] **Step 1: Write failing scan and crossing tests**

Add a seeded smoke test and deterministic synthetic-curve test:

```julia
function synthetic_open_code_scan(rates, small_curve, large_curve)
    points = OpenCodeCapacityScanPoint[]
    for (size, curve) in ((4, small_curve), (8, large_curve))
        for (rate, failure) in zip(rates, curve)
            batches = fill(Float64(failure), 4)
            push!(points, OpenCodeCapacityScanPoint(
                size, Float64(rate), 100, Float64(failure), 0.0,
                Float64(failure), 0.0, batches, batches))
        end
    end
    return OpenCodeCapacityScan([4, 8], Float64[rates...], points)
end

@testset "open code-capacity scan" begin
    scan = scan_open_code_capacity(MersenneTwister(8), [3, 4], [0.0, 0.2];
                                   shots=12, batches=3)
    @test length(scan.points) == 4
    @test all(point -> point.logical_failure_ns_mean == 0.0,
              filter(point -> point.error_rate == 0.0, scan.points))
    @test all(point -> length(point.logical_failure_ns_batches) == 3, scan.points)
end

@testset "crossing statuses" begin
    bracketed = synthetic_open_code_scan([0.08, 0.10, 0.12],
        [0.20, 0.45, 0.70], [0.30, 0.45, 0.60])
    @test only(estimate_open_code_crossings(MersenneTwister(9), bracketed;
        bootstrap=20)).status == :ok

    unbracketed = synthetic_open_code_scan([0.08, 0.10, 0.12],
        [0.20, 0.30, 0.40], [0.10, 0.20, 0.30])
    crossing = only(estimate_open_code_crossings(MersenneTwister(10), unbracketed;
        bootstrap=20))
    @test ismissing(crossing.estimate)
    @test crossing.status == :unbracketed
end
```

- [ ] **Step 2: Run the focused scan test to verify it fails**

Run: `julia --project=. test/code_capacity.jl`

Expected: FAIL because the scan and crossing types/functions are undefined.

- [ ] **Step 3: Implement streaming statistics and generic crossing reuse**

Use the existing `_WelfordAccumulator`, `_isotonic_non_decreasing`, `_selected_crossing`, and batch-bootstrap behavior rather than duplicating their numerical rules. For each `(size, p)`, sample data errors, calculate syndrome, decode each sector separately, score after decoding, and add to sector-specific accumulators and batches:

```julia
errors = sample_data_edge_errors(rng, model; error_rate=error_rate)
syndrome = code_capacity_syndrome(errors)
ns = decode_syndrome(model, syndrome; sector=:north_south, error_rate=error_rate)
ew = decode_syndrome(model, syndrome; sector=:east_west, error_rate=error_rate)
failure_ns = logical_failure(errors, ns; sector=:north_south)
failure_ew = logical_failure(errors, ew; sector=:east_west)
```

Allow `0 <= error_rate < 0.5`; return the all-false correction directly when the syndrome is zero, including at `p = 0`. A nonzero syndrome at `p = 0` is an invalid decoder request. Validate strictly increasing sizes/rates, positive shots, and `2 <= batches <= shots`. Write one progress line such as `completed L=8 p=0.1000 (10000 shots)` when `progress_io` is supplied. Implement crossings by selecting the declared sector's mean/batches and returning `CriticalCrossing` with the existing status rules.

- [ ] **Step 4: Run the focused scan tests to verify they pass**

Run: `julia --project=. test/code_capacity.jl`

Expected: PASS for seeded zero-noise samples, batches, valid crossing, and unbracketed crossing status.

- [ ] **Step 5: Commit scan and threshold estimation support**

```bash
git add src/code_capacity.jl test/code_capacity.jl
git commit -m "feat: scan open code-capacity thresholds"
```

### Task 4: Add independent plot and reproducible command-line runner

**Files:**

- Modify: `src/visualization.jl`
- Modify: `src/TopoNoise.jl`
- Create: `examples/scan_open_code_capacity.jl`
- Create: `test/code_capacity_example.jl`
- Modify: `test/visualization.jl`
- Modify: `test/runtests.jl`

**Interfaces:**

- Produces `plot_open_code_capacity(scan, ns_crossings; ew_crossings=CriticalCrossing[])::CairoMakie.Figure`.
- Produces the example module `ScanOpenCodeCapacityExample` with `main(args=ARGS; io=stdout, error_io=stderr)::Int` and `run(options; io=stdout)`.

- [ ] **Step 1: Write failing plot and CLI tests**

Create `test/code_capacity_example.jl`, include it from `test/runtests.jl`, and include the example module in the test under a unique module name. Verify output separation and visible progress:

```julia
@testset "open code-capacity example" begin
    mktempdir() do directory
        output = IOBuffer()
        status = ScanOpenCodeCapacityExample.main([
            "--p-min", "0.05", "--p-max", "0.10", "--p-step", "0.05",
            "--sizes", "3,4", "--shots", "12", "--batches", "3",
            "--bootstrap", "12", "--output-dir", directory]; io=output)
        @test status == 0
        @test isfile(joinpath(directory, "open_code_capacity_scan.csv"))
        @test isfile(joinpath(directory, "open_code_capacity_scan.pdf"))
        @test occursin("completed L=3", String(take!(output)))
        @test !isfile(joinpath(directory, "trajectory_scan.csv"))
    end
end
```

Extend `test/visualization.jl` to check that the new figure has titles `"Logical failure (N-S)"` and `"Logical failure (E-W)"` and saves nonempty SVG, PDF, and PNG outputs.

- [ ] **Step 2: Run the example and visualization tests to verify they fail**

Run: `julia --project=. test/code_capacity_example.jl && julia --project=. test/visualization.jl`

Expected: FAIL because the independent example and plotting function do not exist.

- [ ] **Step 3: Implement plot and runner without touching trajectory outputs**

Add a two-panel `plot_open_code_capacity` that draws per-size N--S/E--W logical-failure rates with standard-error bars and adds vertical lines only for `:ok` crossings. Copy only the proven option parsing structure from `examples/scan_toric_trajectories.jl`, excluding trajectory-only `--spin-samples`. The runner must:

```julia
scan = scan_open_code_capacity(rng, options.sizes, options.rates;
    shots=options.shots, batches=options.batches, progress_io=io)
ns_crossings = estimate_open_code_crossings(rng, scan; sector=:north_south,
    bootstrap=options.bootstrap, confidence=options.confidence)
ew_crossings = estimate_open_code_crossings(rng, scan; sector=:east_west,
    bootstrap=options.bootstrap, confidence=options.confidence)
```

Write `open_code_capacity_scan.csv`, `open_code_capacity_ns_crossings.csv`, `open_code_capacity_ew_crossings.csv`, and `open_code_capacity_scan.{svg,pdf,png}`. Print absolute output paths and concise primary/cross-check crossing summaries. Keep `trajectory_scan` filenames and source untouched.

- [ ] **Step 4: Run example/figure tests and a small manual command**

Run:

```bash
julia --project=. test/code_capacity_example.jl
julia --project=. test/visualization.jl
julia --project=. examples/scan_open_code_capacity.jl --p-min 0.05 --p-max 0.10 --p-step 0.05 --sizes 3,4 --shots 20 --batches 4 --bootstrap 20 --output-dir /private/tmp/toponoise-open-code-capacity
```

Expected: PASS tests; the command prints progress, emits the five independent output files, and either reports a crossing or labels it unbracketed.

- [ ] **Step 5: Run complete verification and commit the public workflow**

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
git diff --check
```

Expected: full test suite PASS and no whitespace errors.

Commit only feature files:

```bash
git add src/TopoNoise.jl src/visualization.jl examples/scan_open_code_capacity.jl test/code_capacity_example.jl test/visualization.jl test/runtests.jl
git commit -m "feat: add open code-capacity scan example"
```
