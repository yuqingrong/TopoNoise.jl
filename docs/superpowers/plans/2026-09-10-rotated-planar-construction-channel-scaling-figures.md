# Rotated-Planar Construction/Channel Scaling Figures Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add reproducible As/Bp × X-only/Z-only logical-failure scans, exploratory threshold/collapse fits, and one combined plus four standalone publication-style figures.

**Architecture:** Keep the existing balanced `LogicalFailureScan` API untouched. Add a raw comparison subsystem that owns four independent channel scans, a pure scaling-analysis subsystem that fits each raw series, and a CairoMakie output subsystem that renders all five figures from the same fitted comparison object. A dedicated CLI orchestrates data collection, fitting, and artifacts without changing the existing balanced-scan CLI.

**Tech Stack:** Julia 1.10+, TopoNoise.jl, Random, Statistics, LinearAlgebra, Distributions, CairoMakie, PythonCall/PyMatching, Julia Test.

**Spec:** `docs/superpowers/specs/2026-09-10-rotated-planar-construction-channel-scaling-design.md`

## Global Constraints

- Preserve every existing `LogicalFailureScan` and balanced `scan_logical_failure` behavior; it continues to mean `p_x == p_z == p`.
- Use `logical_state=:zero` for every new series.
- `:x_only` means `CircuitPauliNoise(p; p_x=p, p_z=0, clock=clock)` and reports logical-X failure; `:z_only` means `CircuitPauliNoise(p; p_x=0, p_z=p, clock=clock)` and reports logical-Z failure.
- Run four series in this exact order: `(:as, :x_only)`, `(:as, :z_only)`, `(:bp, :x_only)`, `(:bp, :z_only)`.
- Production CLI defaults: distances `[3, 5, 7]`, error grid `0:0.01:0.16`, 10,000 shots and batch size per point, seed `1234`, `:x_ns`, `:gate_layer`, 500 bootstrap replicates.
- Use caller-owned RNG only to derive independent deterministic per-series streams. Store master and derived seeds in results/CSV; reordering series must not change a fixed series’ samples.
- Fit each panel independently. A threshold/exponent result from three distances is always labeled **exploratory**.
- A fit must return structured `:unavailable` status rather than invent `p_c`/`nu` if crossings are missing or numerical collapse is ill-conditioned.
- Do not add balanced-noise panels, repeated syndrome rounds, measurement noise, correlated/circuit-aware decoding, periodic codes, holes, or multi-logical patches.
- Stage only files in the active task; preserve the user’s document build artifacts and untracked Typst/Conda files.

---

## File structure

| File | Responsibility |
|---|---|
| `src/rotated_planar/comparison.jl` | Channel-sweep validation, four-series data model, deterministic simulation runner, and rate/series accessors. |
| `src/rotated_planar/scaling.jl` | Adjacent-distance crossings, aggregation, collapse objective, seeded binomial bootstrap, and fit-result types. |
| `src/rotated_planar/comparison_plot.jl` | Raw/fitted CSV writing, common panel drawing, combined/standalone CairoMakie figures, and artifact saving. |
| `examples/compare_rotated_planar_constructions.jl` | Dependency-free CLI for the new analysis workflow. |
| `test/rotated_planar/comparison.jl` | Channel construction, raw comparison, reproducibility, and existing balanced-scan compatibility tests. |
| `test/rotated_planar/scaling.jl` | Synthetic exact crossing, collapse, bootstrap, and unavailable-fit tests. |
| `test/rotated_planar/comparison_plot.jl` | CSV schema, five-figure artifact, and new CLI smoke tests. |
| `src/TopoNoise.jl`, `Project.toml`, `test/runtests.jl`, `README.md` | Dependency, public API, test-runner, and user-facing integration. |

### Task 1: Channel scans and four-series comparison data

**Files:**
- Create: `src/rotated_planar/comparison.jl`
- Create: `test/rotated_planar/comparison.jl`
- Modify: `src/TopoNoise.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `RotatedPlanarCode`, `rotated_planar_encoder`, `CircuitPauliNoise`, `estimate_logical_failure`, and `LogicalFailurePoint`.
- Produces: `ChannelFailureScan`, `ConstructionChannelComparison`, `scan_channel_logical_failure`, `run_construction_channel_comparison`, `comparison_series`, `channel_failure_count`, `channel_failure_rate`, and `channel_failure_standard_error`.

- [ ] **Step 1: Write failing public-contract and channel-semantics tests**

Create `test/rotated_planar/comparison.jl` with an API guard and a cheap real-PyMatching `p=0`, two-shot scan. The guard makes RED deliberate before the module exists:

```julia
using Random
using Test

@testset "Construction and channel comparison scans" begin
    required = (
        :ChannelFailureScan, :ConstructionChannelComparison,
        :scan_channel_logical_failure, :run_construction_channel_comparison,
        :comparison_series, :channel_failure_count,
        :channel_failure_rate, :channel_failure_standard_error,
    )
    for name in required
        @test isdefined(TopoNoise, name)
    end

    if all(name -> isdefined(TopoNoise, name), required)
        x_scan = scan_channel_logical_failure(
            MersenneTwister(11); construction=:as, error_channel=:x_only,
            distances=[3], error_rates=[0.0], shots=2, batch_size=2, seed=11)
        z_scan = scan_channel_logical_failure(
            MersenneTwister(12); construction=:bp, error_channel=:z_only,
            distances=[3], error_rates=[0.0], shots=2, batch_size=2, seed=12)

        @test x_scan.logical_state === :zero
        @test x_scan.error_channel === :x_only
        @test x_scan.logical_observable === :logical_x
        @test only(x_scan.points).p_x == 0.0
        @test only(x_scan.points).p_z == 0.0
        @test channel_failure_count(x_scan, only(x_scan.points)) == 0
        @test channel_failure_rate(x_scan, only(x_scan.points)) == 0.0
        @test z_scan.logical_observable === :logical_z
        @test channel_failure_count(z_scan, only(z_scan.points)) == 0

        comparison = run_construction_channel_comparison(
            MersenneTwister(1234); distances=[3], error_rates=[0.0],
            shots=2, batch_size=2, seed=1234)
        @test [(scan.construction, scan.error_channel) for scan in comparison.series] ==
            [(:as, :x_only), (:as, :z_only), (:bp, :x_only), (:bp, :z_only)]
        @test comparison_series(comparison, :bp, :z_only) === comparison.series[4]
        @test_throws ArgumentError scan_channel_logical_failure(
            MersenneTwister(1); construction=:as, error_channel=:balanced,
            distances=[3], error_rates=[0.0], shots=2)
    end
end
```

Append this include after the existing rotated-planar scan include in `test/runtests.jl`:

```julia
include("rotated_planar/comparison.jl")
```

- [ ] **Step 2: Run the focused test and verify RED**

Run:

```bash
JULIA_CONDAPKG_BACKEND=Null \
JULIA_PYTHONCALL_EXE=/Users/rongyuqing/jcode/TopoNoise.jl/.CondaPkg/.pixi/envs/default/bin/python3 \
JULIA_DEPOT_PATH=/private/tmp/toponoise-julia-depot:/Users/rongyuqing/.julia \
julia --project=. -e 'using TopoNoise, Random, Test; include("test/rotated_planar/comparison.jl")'
```

Expected: failures only for the eight absent exported symbols.

- [ ] **Step 3: Add the raw comparison model and simulation runner**

Create `src/rotated_planar/comparison.jl`. Use these exact result types and channel helpers:

```julia
const _COMPARISON_SERIES_ORDER = (
    (:as, :x_only), (:as, :z_only), (:bp, :x_only), (:bp, :z_only),
)

struct ChannelFailureScan
    construction::Symbol
    error_channel::Symbol
    logical_observable::Symbol
    logical_state::Symbol
    boundary_orientation::Symbol
    clock::Symbol
    distances::Vector{Int}
    error_rates::Vector{Float64}
    shots::Int
    batch_size::Int
    master_seed::Union{Nothing,Int}
    series_seed::UInt64
    points::Vector{LogicalFailurePoint}
end

struct ConstructionChannelComparison
    logical_state::Symbol
    boundary_orientation::Symbol
    clock::Symbol
    distances::Vector{Int}
    error_rates::Vector{Float64}
    shots::Int
    batch_size::Int
    master_seed::Union{Nothing,Int}
    series::Vector{ChannelFailureScan}
end

_channel_observable(:x_only) = :logical_x
_channel_observable(:z_only) = :logical_z
_channel_observable(channel::Symbol) = throw(ArgumentError(
    "error_channel must be :x_only or :z_only"))

function _channel_noise(rate::Float64, channel::Symbol, clock::Symbol)
    channel === :x_only && return CircuitPauliNoise(rate; p_x=rate, p_z=0, clock)
    channel === :z_only && return CircuitPauliNoise(rate; p_x=0, p_z=rate, clock)
    _channel_observable(channel)
end

channel_failure_count(scan::ChannelFailureScan, point::LogicalFailurePoint) =
    scan.logical_observable === :logical_x ? point.logical_x_failures : point.logical_z_failures
channel_failure_rate(scan::ChannelFailureScan, point::LogicalFailurePoint) =
    scan.logical_observable === :logical_x ? point.logical_x_failure_rate : point.logical_z_failure_rate
channel_failure_standard_error(scan::ChannelFailureScan, point::LogicalFailurePoint) =
    scan.logical_observable === :logical_x ? point.logical_x_standard_error : point.logical_z_standard_error
```

Implement `scan_channel_logical_failure` by reusing `_scan_distances`, `_scan_error_rates`, `_positive_machine_int`, `_seed_metadata`, `RotatedPlanarCode`, `rotated_planar_encoder`, and `estimate_logical_failure`; traverse distance-major/rate-minor. Implement `run_construction_channel_comparison` by drawing one `UInt64` seed from the caller RNG per tuple in `_COMPARISON_SERIES_ORDER`, constructing `MersenneTwister(series_seed)`, and collecting exactly four `ChannelFailureScan`s. `comparison_series` must throw `ArgumentError` for absent pairs.

Add exports and includes in `src/TopoNoise.jl` after the existing scan exports/includes.

- [ ] **Step 4: Extend the test with real nonzero X-only/Z-only metadata and deterministic streams**

Add these assertions below the initial test:

```julia
nonzero_x = scan_channel_logical_failure(
    MersenneTwister(22); construction=:bp, error_channel=:x_only,
    distances=[3], error_rates=[0.1], shots=4, batch_size=4, seed=22)
nonzero_z = scan_channel_logical_failure(
    MersenneTwister(23); construction=:bp, error_channel=:z_only,
    distances=[3], error_rates=[0.1], shots=4, batch_size=4, seed=23)
@test only(nonzero_x.points).p_x == 0.1
@test only(nonzero_x.points).p_z == 0.0
@test only(nonzero_z.points).p_x == 0.0
@test only(nonzero_z.points).p_z == 0.1

a = run_construction_channel_comparison(
    MersenneTwister(99); distances=[3], error_rates=[0.0], shots=2, batch_size=2, seed=99)
b = run_construction_channel_comparison(
    MersenneTwister(99); distances=[3], error_rates=[0.0], shots=2, batch_size=2, seed=99)
@test [scan.series_seed for scan in a.series] == [scan.series_seed for scan in b.series]
@test [(
    scan.construction, scan.error_channel, scan.logical_observable,
    [(point.p_x, point.p_z, point.logical_x_failures, point.logical_z_failures)
     for point in scan.points],
) for scan in a.series] == [(
    scan.construction, scan.error_channel, scan.logical_observable,
    [(point.p_x, point.p_z, point.logical_x_failures, point.logical_z_failures)
     for point in scan.points],
) for scan in b.series]
```

- [ ] **Step 5: Run focused and legacy balanced-scan tests to verify GREEN**

Run:

```bash
JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE=/Users/rongyuqing/jcode/TopoNoise.jl/.CondaPkg/.pixi/envs/default/bin/python3 JULIA_DEPOT_PATH=/private/tmp/toponoise-julia-depot:/Users/rongyuqing/.julia julia --project=. -e 'using TopoNoise, Random, Test; include("test/rotated_planar/comparison.jl"); include("test/rotated_planar/scan.jl")'
```

Expected: the new comparison test and existing balanced scan/artifact/CLI tests pass.

- [ ] **Step 6: Commit the scoped raw-comparison change**

```bash
git add src/rotated_planar/comparison.jl src/TopoNoise.jl test/rotated_planar/comparison.jl test/runtests.jl
git commit -m "feat: add construction channel comparison scans"
```

### Task 2: Exploratory crossings, collapse fit, and bootstrap

**Files:**
- Create: `src/rotated_planar/scaling.jl`
- Create: `test/rotated_planar/scaling.jl`
- Modify: `Project.toml`
- Modify: `src/TopoNoise.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `ChannelFailureScan`, `ConstructionChannelComparison`, and channel-rate accessors from Task 1.
- Produces: `PairCrossing`, `ThresholdFit`, `FittedConstructionChannelComparison`, `fit_channel_threshold`, `fit_construction_channel_comparison`, and `scaled_error_rate`.

- [ ] **Step 1: Write synthetic fitting tests before adding an analysis module**

Create `test/rotated_planar/scaling.jl`. Define a local helper that creates a `ChannelFailureScan` with three distances and synthetic `LogicalFailurePoint`s whose selected rates are `0.25 + slope[d] * (p - 0.1)`. Use 10,000 shots and populate the selected count/rate/standard-error fields consistently. Then test the public behavior:

```julia
@testset "Exploratory scaling fit" begin
    scan = synthetic_channel_scan(; crossing=0.1)
    fit = fit_channel_threshold(scan; bootstrap_replicates=20, bootstrap_seed=7)
    @test fit.status === :success
    @test isapprox(fit.p_c, 0.1; atol=1e-10)
    @test 0.5 <= fit.nu <= 3.0
    @test length(fit.crossings) == 2
    @test fit.bootstrap_replicates == 20
    @test fit.exploratory
    @test scaled_error_rate(0.11, 5, fit) ==
        (0.11 - fit.p_c) * 5^(1 / fit.nu)

    unavailable = synthetic_channel_scan(; crossing=nothing)
    failed = fit_channel_threshold(unavailable; bootstrap_replicates=10, bootstrap_seed=8)
    @test failed.status === :unavailable
    @test failed.p_c === nothing
    @test failed.nu === nothing
    @test !isempty(failed.diagnostic)
end
```

Include it after `comparison.jl` in `test/runtests.jl`.

- [ ] **Step 2: Run the focused fitting test and verify RED**

Run:

```bash
JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE=/Users/rongyuqing/jcode/TopoNoise.jl/.CondaPkg/.pixi/envs/default/bin/python3 JULIA_DEPOT_PATH=/private/tmp/toponoise-julia-depot:/Users/rongyuqing/.julia julia --project=. -e 'using TopoNoise, Random, Test; include("test/rotated_planar/scaling.jl")'
```

Expected: absent type/function failures only.

- [ ] **Step 3: Declare the bootstrap dependency and add fit result types**

Add `Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"` to `[deps]` and `Distributions = "0.25"` to `[compat]` in `Project.toml`; add `using Distributions: Binomial` in `src/TopoNoise.jl`.

Create `src/rotated_planar/scaling.jl` with these exact data types:

```julia
struct PairCrossing
    lower_distance::Int
    upper_distance::Int
    p::Float64
    standard_error::Float64
end

struct ThresholdFit
    construction::Symbol
    error_channel::Symbol
    logical_observable::Symbol
    status::Symbol
    diagnostic::String
    crossings::Vector{PairCrossing}
    p_c::Union{Nothing,Float64}
    p_c_standard_error::Union{Nothing,Float64}
    nu::Union{Nothing,Float64}
    nu_standard_error::Union{Nothing,Float64}
    bootstrap_replicates::Int
    bootstrap_successes::Int
    bootstrap_seed::UInt64
    exploratory::Bool
end

struct FittedConstructionChannelComparison
    raw::ConstructionChannelComparison
    fits::Vector{ThresholdFit}
end
```

- [ ] **Step 4: Implement crossings, collapse objective, and bootstrap without touching simulation code**

Implement `_adjacent_crossings(scan)` by comparing two adjacent-distance curves at matching physical-rate points. For a difference `delta_i = y_low[i] - y_high[i]`, accept a crossing only when `delta_i == 0` or `delta_i * delta_{i+1} < 0`; use

```julia
fraction = delta_i / (delta_i - delta_next)
p_cross = p_i + fraction * (p_next - p_i)
```

and obtain crossing uncertainty by linear interpolation of `sqrt(se_low^2 + se_high^2)`, with a finite positive floor of `eps(Float64)`.

Aggregate at least two crossings with inverse-variance weights:

```julia
weights = 1.0 ./ getfield.(crossings, :standard_error).^2
p_c = sum(weights .* getfield.(crossings, :p)) / sum(weights)
p_c_se = sqrt(inv(sum(weights)))
```

For `nu`, search `0.50:0.01:3.00`. For each candidate, transform every curve’s x coordinates by `(p - p_c) * d^(1/nu)`, select the intersection of all curve x ranges, interpolate every curve at 64 equally spaced common x locations, and minimize the inverse-variance weighted sum of pairwise rate residuals. Require finite objective values and at least four common interpolation locations; otherwise return `:unavailable`.

For each bootstrap replicate, sample each selected failure count from `Binomial(point.shots, selected_rate)`, rebuild only the selected channel-rate fields, rerun the same crossing/collapse path, and retain successful `p_c`/`nu` values. Use a `MersenneTwister(bootstrap_seed)` and set errors to sample standard deviations only when at least two bootstrap fits succeed. `fit_channel_threshold` must return `:unavailable` without bootstrapping if raw crossings/fit are unavailable. `fit_construction_channel_comparison` fits the four raw scans in their stored order and derives one deterministic bootstrap seed per fit from the supplied RNG.

Export the five public names and include `scaling.jl` after `comparison.jl`.

- [ ] **Step 5: Add bootstrap reproducibility and all-panel wrapper tests**

Extend `test/rotated_planar/scaling.jl`:

```julia
first = fit_channel_threshold(scan; bootstrap_replicates=20, bootstrap_seed=42)
second = fit_channel_threshold(scan; bootstrap_replicates=20, bootstrap_seed=42)
@test (
    first.status, first.p_c, first.p_c_standard_error, first.nu,
    first.nu_standard_error, first.bootstrap_successes,
) == (
    second.status, second.p_c, second.p_c_standard_error, second.nu,
    second.nu_standard_error, second.bootstrap_successes,
)

synthetic_series = [
    ChannelFailureScan(
        construction, channel, _channel_observable(channel), :zero, :x_ns,
        :gate_layer, copy(scan.distances), copy(scan.error_rates), scan.shots,
        scan.batch_size, 3, UInt64(index), copy(scan.points),
    )
    for (index, (construction, channel)) in
        enumerate(((:as, :x_only), (:as, :z_only), (:bp, :x_only), (:bp, :z_only)))
]
raw = ConstructionChannelComparison(
    :zero, :x_ns, :gate_layer, copy(scan.distances), copy(scan.error_rates),
    scan.shots, scan.batch_size, 3, synthetic_series,
)
fitted = fit_construction_channel_comparison(
    MersenneTwister(3), raw; bootstrap_replicates=10)
@test length(fitted.fits) == 4
@test [(fit.construction, fit.error_channel) for fit in fitted.fits] ==
    [(:as, :x_only), (:as, :z_only), (:bp, :x_only), (:bp, :z_only)]
```

The helper must construct all fields named in Task 1; do not use an untyped dictionary or `Any` container.

- [ ] **Step 6: Run focused fit tests and commit**

Run:

```bash
JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE=/Users/rongyuqing/jcode/TopoNoise.jl/.CondaPkg/.pixi/envs/default/bin/python3 JULIA_DEPOT_PATH=/private/tmp/toponoise-julia-depot:/Users/rongyuqing/.julia julia --project=. -e 'using TopoNoise, Random, Test; include("test/rotated_planar/scaling.jl")'
```

Expected: known-crossing, unavailable, bootstrap, and wrapper tests pass.

```bash
git add Project.toml src/TopoNoise.jl src/rotated_planar/scaling.jl test/rotated_planar/scaling.jl test/runtests.jl
git commit -m "feat: add exploratory scaling fits"
```

### Task 3: Comparison CSV and five CairoMakie figures

**Files:**
- Create: `src/rotated_planar/comparison_plot.jl`
- Create: `test/rotated_planar/comparison_plot.jl`
- Modify: `src/TopoNoise.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `FittedConstructionChannelComparison`, `ThresholdFit`, and channel accessors.
- Produces: `write_comparison_raw_csv`, `write_comparison_fit_csv`, `plot_construction_channel_comparison`, `plot_construction_channel_panel`, and `save_construction_channel_comparison`.

- [ ] **Step 1: Write failing CSV and artifact contract tests**

Create `test/rotated_planar/comparison_plot.jl` using a tiny fitted synthetic comparison from Task 2. Test exact CSV headers and the complete artifact surface:

```julia
@testset "Construction/channel comparison artifacts" begin
    fitted = synthetic_fitted_comparison()
    mktempdir() do directory
        paths = save_construction_channel_comparison(
            fitted, directory; basename="comparison")
        @test keys(paths) == (:raw_csv, :fits_csv, :combined, :panels)
        @test isfile(paths.raw_csv)
        @test isfile(paths.fits_csv)
        @test readlines(paths.raw_csv)[1] ==
            "construction,error_channel,logical_observable,distance,p_x,p_z,shots,master_seed,series_seed,logical_failures,logical_failure_rate,logical_failure_standard_error,logical_state,boundary_orientation,clock"
        @test readlines(paths.fits_csv)[1] ==
            "construction,error_channel,logical_observable,status,diagnostic,crossing_lower_distance,crossing_upper_distance,crossing_p,crossing_standard_error,p_c,p_c_standard_error,nu,nu_standard_error,bootstrap_replicates,bootstrap_successes,bootstrap_seed,exploratory"
        @test length(paths.panels) == 4
        for group in (paths.combined, values(paths.panels)...), path in values(group)
            @test isfile(path)
            @test filesize(path) > 100
        end
        @test plot_construction_channel_comparison(fitted) isa CairoMakie.Figure
        @test plot_construction_channel_panel(fitted, :as, :x_only) isa CairoMakie.Figure
    end
end
```

Include it after `scaling.jl` in `test/runtests.jl`.

- [ ] **Step 2: Run the focused artifact test and verify RED**

Run:

```bash
JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE=/Users/rongyuqing/jcode/TopoNoise.jl/.CondaPkg/.pixi/envs/default/bin/python3 JULIA_DEPOT_PATH=/private/tmp/toponoise-julia-depot:/Users/rongyuqing/.julia MPLCONFIGDIR=/private/tmp/toponoise-comparison-matplotlib julia --project=. -e 'using TopoNoise, Random, Test; include("test/rotated_planar/comparison_plot.jl")'
```

Expected: absent artifact-function failures only.

- [ ] **Step 3: Implement exact raw/fitted CSV rows**

Create `src/rotated_planar/comparison_plot.jl`. Define column tuples exactly as tested. `write_comparison_raw_csv(path, fitted)` must loop over `fitted.raw.series`, then each point, and use `channel_failure_count/rate/standard_error` rather than `state_failure_rate`. `write_comparison_fit_csv(path, fitted)` must write one row per `PairCrossing`, duplicating the owning fit fields so the CSV contains both pairwise crossing evidence and the panel estimate; write one row with empty crossing fields when a fit has no crossings. Emit empty fields for `nothing`, and preserve diagnostics in a CSV-safe quoted field by replacing `"` with `""` and wrapping the field in quotes.

- [ ] **Step 4: Implement one reusable panel renderer and both figure layouts**

Implement `_comparison_panel!(position, scan, fit; title)` to draw physical p versus the selected channel rate, one color/error-bar curve per distance, and an `Exploratory` title suffix. When `fit.status === :success`, add `vlines!` at `fit.p_c`, annotate `p_c` and `nu` with bootstrap errors, and create an inset Axis using `scaled_error_rate(point.p_x, point.distance, fit)`. When unavailable, add the diagnostic text and do not create a threshold line or inset.

`plot_construction_channel_comparison(fitted)` must create a `Figure(size=(1400, 980))`, use rows `:as`, `:bp` and columns `:x_only`, `:z_only`, and call the same renderer exactly four times. `plot_construction_channel_panel(fitted, construction, error_channel)` must validate the pair and call the renderer once in `Figure(size=(700, 560))`.

- [ ] **Step 5: Implement safe deterministic save paths and strengthen layout assertions**

Implement:

```julia
function save_construction_channel_comparison(
        fitted::FittedConstructionChannelComparison,
        output_dir::AbstractString;
        basename::AbstractString="rotated_planar_comparison")
```

Reject empty, `.`/`..`, slash-containing, and backslash-containing basenames. Make `output_dir`, write `basename-raw.csv` and `basename-fits.csv`, then save SVG/PDF/PNG for `basename-combined` and the four `basename-as-x-only`, `basename-as-z-only`, `basename-bp-x-only`, and `basename-bp-z-only` figures. Return:

```julia
(
    raw_csv=raw_csv_path,
    fits_csv=fits_csv_path,
    combined=(svg=combined_svg_path, pdf=combined_pdf_path, png=combined_png_path),
    panels=Dict(
        (:as, :x_only) => (svg=as_x_svg_path, pdf=as_x_pdf_path, png=as_x_png_path),
        (:as, :z_only) => (svg=as_z_svg_path, pdf=as_z_pdf_path, png=as_z_png_path),
        (:bp, :x_only) => (svg=bp_x_svg_path, pdf=bp_x_pdf_path, png=bp_x_png_path),
        (:bp, :z_only) => (svg=bp_z_svg_path, pdf=bp_z_pdf_path, png=bp_z_png_path),
    ),
)
```

Add a test that `basename="../escape"` throws `ArgumentError` and a structural test that the combined figure has four primary `Axis` objects.

- [ ] **Step 6: Run artifact tests and commit**

Run the Task 3 command from Step 2. Expected: all CSV, figure, unavailable-fit, and output-path tests pass.

```bash
git add src/rotated_planar/comparison_plot.jl src/TopoNoise.jl test/rotated_planar/comparison_plot.jl test/runtests.jl
git commit -m "feat: plot construction channel scaling comparisons"
```

### Task 4: Dedicated CLI, README, and final integration

**Files:**
- Create: `examples/compare_rotated_planar_constructions.jl`
- Modify: `README.md`
- Modify: `test/rotated_planar/comparison_plot.jl`

**Interfaces:**
- Consumes: `run_construction_channel_comparison`, `fit_construction_channel_comparison`, and `save_construction_channel_comparison`.
- Produces: `RotatedPlanarConstructionComparison.main(args=ARGS; io=stdout, error_io=stderr)::Int` and documented production command.

- [ ] **Step 1: Write the failing CLI smoke/validation tests**

Append to `test/rotated_planar/comparison_plot.jl`:

```julia
const _COMPARISON_CLI = joinpath(
    @__DIR__, "..", "..", "examples", "compare_rotated_planar_constructions.jl")

@testset "Construction/channel comparison CLI" begin
    @test isfile(_COMPARISON_CLI)
    if isfile(_COMPARISON_CLI)
        include(_COMPARISON_CLI)
        mktempdir() do directory
            output, errors = IOBuffer(), IOBuffer()
            result = RotatedPlanarConstructionComparison.main([
                "--distances", "3", "--error-rates", "0", "--shots", "2",
                "--batch-size", "2", "--bootstrap-replicates", "2", "--no-fit",
                "--seed", "5", "--output-dir", directory, "--basename", "smoke",
            ]; io=output, error_io=errors)
            @test result == 0
            @test isfile(joinpath(directory, "smoke-combined.svg"))
            @test isfile(joinpath(directory, "smoke-as-x-only.png"))
            @test contains(String(take!(output)), "fit unavailable")
            @test isempty(String(take!(errors)))
        end
        for arguments in (["--bootstrap-replicates", "0"], ["--error-rates", "0.1,0.1"], ["--unknown"])
            errors = IOBuffer()
            @test RotatedPlanarConstructionComparison.main(arguments; error_io=errors) == 1
            @test contains(String(take!(errors)), "error:")
        end
    end
end
```

- [ ] **Step 2: Run the CLI smoke test and verify RED**

Run the Task 3 focused test command. Expected: the `isfile(_COMPARISON_CLI)` assertion fails before the script exists.

- [ ] **Step 3: Implement the standalone CLI module**

Create `examples/compare_rotated_planar_constructions.jl` with module `RotatedPlanarConstructionComparison`. Reuse the proven dependency-free argument-parser style from `examples/scan_rotated_planar.jl`, but implement these options and defaults:

```text
--distances 3,5,7
--error-rates LIST | --p-min 0 --p-max 0.16 --p-step 0.01
--shots 10000 --batch-size 10000 --seed 1234
--boundary-orientation x_ns --clock gate_layer
--bootstrap-replicates 500 --no-fit
--output-dir results/rotated-planar-comparison
--basename rotated_planar_comparison
--help
```

Reject combining `--error-rates` with any grid argument, require finite strictly increasing rates in `[0, 0.5)`, require odd strictly increasing distances at least 3, require positive shots/batch/bootstrap counts, and validate only leaf basenames before running simulations. `--no-fit` constructs one `ThresholdFit` with `status=:unavailable`, diagnostic `"fit disabled by --no-fit"`, and no bootstrap; it still writes all raw figures and CSVs.

With fitting enabled, construct `MersenneTwister(seed)`, call the raw runner, call the fitting wrapper with the same caller-owned RNG, save artifacts, print each path, then print one construction/channel/status/`p_c`/`nu` summary per panel. The program-file guard must call `exit(main())` only when executed directly.

- [ ] **Step 4: Document scientific interpretation, command, outputs, and caution**

Add a `Construction/channel scaling comparison` README section immediately after the existing balanced scan section. Include:

```bash
julia --project=. examples/compare_rotated_planar_constructions.jl \
  --distances 3,5,7 --p-min 0 --p-max 0.16 --p-step 0.01 \
  --shots 10000 --batch-size 10000 --seed 1234 \
  --output-dir results/rotated-planar-comparison
```

State exactly: all panels prepare `|0_L>`; X-only uses `p_x=p,p_z=0` and reports logical-X failure; Z-only uses `p_x=0,p_z=p` and reports logical-Z failure; As/Bp choose circuit construction only. List the combined and four standalone figures, raw/fits CSV files, per-panel `p_c`/`nu` inset, and the fact that d=3/5/7 fit values are exploratory rather than threshold claims.

- [ ] **Step 5: Run focused CLI/artifact tests, the complete suite, and a CLI preview**

Run:

```bash
JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE=/Users/rongyuqing/jcode/TopoNoise.jl/.CondaPkg/.pixi/envs/default/bin/python3 JULIA_DEPOT_PATH=/private/tmp/toponoise-julia-depot:/Users/rongyuqing/.julia MPLCONFIGDIR=/private/tmp/toponoise-comparison-matplotlib julia --project=. -e 'using TopoNoise, Random, Test; include("test/rotated_planar/comparison_plot.jl")'
```

Then run the repository test suite:

```bash
JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE=/Users/rongyuqing/jcode/TopoNoise.jl/.CondaPkg/.pixi/envs/default/bin/python3 JULIA_DEPOT_PATH=/private/tmp/toponoise-julia-depot:/Users/rongyuqing/.julia MPLCONFIGDIR=/private/tmp/toponoise-comparison-matplotlib julia --project=. test/runtests.jl
```

Finally run a retained, cheap real CLI preview:

```bash
JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE=/Users/rongyuqing/jcode/TopoNoise.jl/.CondaPkg/.pixi/envs/default/bin/python3 JULIA_DEPOT_PATH=/private/tmp/toponoise-julia-depot:/Users/rongyuqing/.julia julia --project=. examples/compare_rotated_planar_constructions.jl --distances 3 --error-rates 0 --shots 2 --batch-size 2 --bootstrap-replicates 2 --no-fit --output-dir results/rotated-planar-comparison-smoke --basename smoke
```

Expected: all tests pass and the preview writes one combined and four standalone PNG/SVG/PDF outputs plus two CSV files under `results/rotated-planar-comparison-smoke`.

- [ ] **Step 6: Perform scope checks and commit**

Run:

```bash
git diff --check
git status --short
```

Confirm that only Task 4 files are staged, then commit:

```bash
git add examples/compare_rotated_planar_constructions.jl README.md test/rotated_planar/comparison_plot.jl
git commit -m "feat: add construction channel scaling CLI"
```

## Final verification checklist

- [ ] Existing balanced `scan_logical_failure` and `plot_logical_failure_scan` tests remain green without behavior changes.
- [ ] Every new X-only point has `p_z == 0`; every new Z-only point has `p_x == 0`.
- [ ] Every comparison starts from `logical_state=:zero`; no state-failure metric is substituted for the selected channel logical metric.
- [ ] Combined and standalone panels use the exact same raw points and fit result for each construction/channel pair.
- [ ] Unavailable fits retain raw data/artifacts and explicit diagnostics.
- [ ] New data output and figure output are reproducible from the configured seed.
- [ ] Full test suite passes in the established PyMatching environment and `git diff --check` is clean.
