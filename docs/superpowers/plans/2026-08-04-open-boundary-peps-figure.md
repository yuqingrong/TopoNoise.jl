# Open-Boundary Toric-Code PEPS Figure Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Change only `plot_peps_graph` so every displayed tensor has four directed virtual half-legs, four northeast physical legs, explicit open-boundary stubs, and the approved representative-site labels.

**Architecture:** Add a small pure geometry model inside `src/visualization.jl` that describes every leg independently of CairoMakie. Test that model exhaustively, then make the existing plot function render it with `arrows2d!`; keep `peps_graph(peps)` as the validation gate and leave all PEPS and circuit data unchanged.

**Tech Stack:** Julia 1.10+, ITensors 0.9, ITensorNetworks 0.21, CairoMakie 0.15, Test stdlib.

## Global Constraints

- This is a figure-only change; do not modify `ToricCodePEPS`, `toric_code_peps`, `peps_graph`, Yao circuit construction, boundary contractions, normalization, or wavefunctions.
- Preserve tensor coordinates `x=col`, `y=R-row+1`.
- Preserve E/N/W/S mappings `(α,β,γ,δ)` for virtual legs and `(i,j,k,l)` for physical legs.
- Horizontal virtual arrows point east; vertical virtual arrows point north.
- Every physical arrow slopes northeast.
- Show symbolic leg labels only at the lower-left site `(R,1)`.
- Preserve `plot_peps_graph(peps; show_index_labels=false)::CairoMakie.Figure` and the optional compact physical-index overlay.
- Preserve vector SVG/PDF export through `CairoMakie.save`.

---

### Task 1: Pure eight-leg geometry

**Files:**
- Modify: `src/visualization.jl:1-90`
- Modify: `test/visualization.jl:1-31`

**Interfaces:**
- Consumes: `peps::ToricCodePEPS`, `size(peps)`, and `_peps_position(rows, site)`.
- Produces: `_PEPSLeg` and `_peps_leg_geometry(peps)::Vector{_PEPSLeg}` for the renderer and tests.

- [ ] **Step 1: Write the failing geometry tests**

Add this testset before the existing vector-figure testset in
`test/visualization.jl`:

```julia
@testset "Open-boundary PEPS leg geometry" begin
    peps = toric_code_peps(2, 2)
    legs = TopoNoise._peps_leg_geometry(peps)
    virtual = filter(leg -> leg.kind == :virtual, legs)
    physical = filter(leg -> leg.kind == :physical, legs)

    @test length(legs) == 8 * 4
    @test length(virtual) == 4 * 4
    @test length(physical) == 4 * 4
    @test count(leg -> leg.boundary, virtual) == 2 * 2 + 2 * 2

    virtual_symbols = Dict(
        :east => "α", :north => "β", :west => "γ", :south => "δ")
    physical_symbols = Dict(
        :east => "i", :north => "j", :west => "k", :south => "l")

    for row in 1:2, col in 1:2
        site = (row, col)
        site_virtual = filter(
            leg -> leg.site == site && leg.kind == :virtual, legs)
        site_physical = filter(
            leg -> leg.site == site && leg.kind == :physical, legs)
        @test Set(leg.direction for leg in site_virtual) ==
              Set((:east, :north, :west, :south))
        @test Set(leg.direction for leg in site_physical) ==
              Set((:east, :north, :west, :south))
        @test all(leg.symbol == virtual_symbols[leg.direction]
                  for leg in site_virtual)
        @test all(leg.symbol == physical_symbols[leg.direction]
                  for leg in site_physical)
    end

    horizontal = filter(
        leg -> leg.kind == :virtual && leg.direction in (:east, :west), legs)
    vertical = filter(
        leg -> leg.kind == :virtual && leg.direction in (:north, :south), legs)
    @test all(leg.tip[1] > leg.tail[1] && leg.tip[2] == leg.tail[2]
              for leg in horizontal)
    @test all(leg.tip[2] > leg.tail[2] && leg.tip[1] == leg.tail[1]
              for leg in vertical)
    @test all(leg.tip[1] > leg.tail[1] && leg.tip[2] > leg.tail[2]
              for leg in physical)

    east_half = only(filter(
        leg -> leg.site == (2, 1) && leg.kind == :virtual &&
               leg.direction == :east, legs))
    west_half = only(filter(
        leg -> leg.site == (2, 2) && leg.kind == :virtual &&
               leg.direction == :west, legs))
    @test east_half.tip == west_half.tail

    north_half = only(filter(
        leg -> leg.site == (2, 1) && leg.kind == :virtual &&
               leg.direction == :north, legs))
    south_half = only(filter(
        leg -> leg.site == (1, 1) && leg.kind == :virtual &&
               leg.direction == :south, legs))
    @test north_half.tip == south_half.tail
end
```

- [ ] **Step 2: Run the focused test and confirm the red state**

Run:

```bash
julia --project=. -e 'using TopoNoise, Test, ITensors; include("test/visualization.jl")'
```

Expected: error or failure because `TopoNoise._peps_leg_geometry` is not
defined.

- [ ] **Step 3: Add the pure leg representation**

Add near `_peps_position` in `src/visualization.jl`:

```julia
struct _PEPSLeg
    site::Tuple{Int,Int}
    kind::Symbol
    direction::Symbol
    tail::NTuple{2,Float64}
    tip::NTuple{2,Float64}
    symbol::String
    boundary::Bool
end

const _VIRTUAL_SYMBOLS = (
    east="α", north="β", west="γ", south="δ")
const _PHYSICAL_SYMBOLS = (
    east="i", north="j", west="k", south="l")
```

Implement virtual half-leg geometry with a node clearance of `0.20`, an
internal midpoint at `0.50`, and an open-boundary endpoint at `0.72`:

```julia
function _virtual_peps_leg(
        rows::Int, cols::Int, site::Tuple{Int,Int}, direction::Symbol)
    row, col = site
    x, y = _peps_position(rows, site)
    symbol = getproperty(_VIRTUAL_SYMBOLS, direction)
    boundary = direction == :east  ? col == cols :
               direction == :north ? row == 1 :
               direction == :west  ? col == 1 : row == rows
    reach = boundary ? 0.72 : 0.50

    tail, tip = if direction == :east
        ((x + 0.20, y), (x + reach, y))
    elseif direction == :north
        ((x, y + 0.20), (x, y + reach))
    elseif direction == :west
        ((x - reach, y), (x - 0.20, y))
    else
        ((x, y - reach), (x, y - 0.20))
    end
    return _PEPSLeg(site, :virtual, direction, tail, tip, symbol, boundary)
end
```

Implement physical legs with one anchor on each local E/N/W/S port and a
positive `x` and `y` displacement:

```julia
function _physical_peps_leg(
        rows::Int, site::Tuple{Int,Int}, direction::Symbol)
    x, y = _peps_position(rows, site)
    symbol = getproperty(_PHYSICAL_SYMBOLS, direction)
    tail = direction == :east  ? (x + 0.22, y - 0.08) :
           direction == :north ? (x + 0.08, y + 0.22) :
           direction == :west  ? (x - 0.22, y + 0.08) :
                                 (x - 0.08, y - 0.22)
    tip = (tail[1] + 0.24, tail[2] + 0.30)
    return _PEPSLeg(site, :physical, direction, tail, tip, symbol, false)
end

function _peps_leg_geometry(peps::ToricCodePEPS)
    rows, cols = size(peps)
    directions = (:east, :north, :west, :south)
    legs = _PEPSLeg[]
    for row in 1:rows, col in 1:cols
        site = (row, col)
        for direction in directions
            push!(legs, _virtual_peps_leg(rows, cols, site, direction))
        end
        for direction in directions
            push!(legs, _physical_peps_leg(rows, site, direction))
        end
    end
    return legs
end
```

Keep the exact ordering virtual E/N/W/S followed by physical E/N/W/S per site.

- [ ] **Step 4: Run the focused geometry tests and confirm green**

Run the same focused command from Step 2.

Expected: the geometry testset passes; the pre-existing figure testset also
continues to pass because the renderer has not changed yet.

- [ ] **Step 5: Commit the geometry model and tests**

```bash
git add src/visualization.jl test/visualization.jl
git commit -m "test: define open-boundary PEPS leg geometry"
```

---

### Task 2: Render the approved open-boundary isometry diagram

**Files:**
- Modify: `src/visualization.jl:20-100`
- Modify: `test/visualization.jl:15-70`
- Modify: `README.md` in the PEPS graph example section

**Interfaces:**
- Consumes: `_peps_leg_geometry(peps)::Vector{_PEPSLeg}` from Task 1.
- Produces: the unchanged public function `plot_peps_graph(peps; show_index_labels=false)::CairoMakie.Figure` with the approved eight-leg rendering.

- [ ] **Step 1: Extend the figure tests for representative labels**

In the existing `"Toric-code vector figures"` testset, construct both default
and diagnostic figures and assert that symbolic labels occur exactly once:

```julia
peps_figure = plot_peps_graph(peps)
peps_index_figure = plot_peps_graph(peps; show_index_labels=true)

peps_labels = rendered_text(peps_figure)
@test all("T[$row,$col]" in peps_labels for row in 1:2 for col in 1:2)
@test all(count(==(symbol), peps_labels) == 1
          for symbol in ("α", "β", "γ", "δ", "i", "j", "k", "l"))

diagnostic_labels = rendered_text(peps_index_figure)
@test all(any(startswith("$direction:"), diagnostic_labels)
          for direction in ("E", "N", "W", "S"))
```

Keep SVG and PDF export assertions, saving `peps_figure` for the default
publication form.

- [ ] **Step 2: Run the focused test and confirm the red state**

Run:

```bash
julia --project=. -e 'using TopoNoise, Test, ITensors; include("test/visualization.jl")'
```

Expected: failures because the current renderer does not emit the eight
representative symbolic labels.

- [ ] **Step 3: Replace bond/stub lines with arrow rendering**

Keep `peps_graph(peps)` at the start of `plot_peps_graph` to validate input,
then obtain `legs = _peps_leg_geometry(peps)`. Remove the old iteration over
`ITensorNetworks.edges(network)` and the old four red stub geometry.

Add a focused drawing helper:

```julia
function _draw_peps_arrows!(axis, legs::Vector{_PEPSLeg}, kind::Symbol)
    selected = filter(leg -> leg.kind == kind, legs)
    tails = [CairoMakie.Point2f(leg.tail...) for leg in selected]
    vectors = [CairoMakie.Vec2f(
        leg.tip[1] - leg.tail[1], leg.tip[2] - leg.tail[2])
        for leg in selected]
    color = kind == :virtual ? _PEPS_BOND_COLOR : _PEPS_PHYSICAL_COLOR
    CairoMakie.arrows2d!(axis, tails, vectors;
        color=color, shaftwidth=2.2, tipwidth=10, tiplength=7)
    return axis
end
```

Call `_draw_peps_arrows!` first for virtual legs and then physical legs. Draw
tensor nodes after all arrows so arrow attachments terminate cleanly beneath
the blue node glyph.

- [ ] **Step 4: Add representative symbolic and optional diagnostic labels**

Add deterministic label placement so labels do not cover arrow shafts:

```julia
function _peps_symbol_label_position(leg::_PEPSLeg)
    midpoint = (
        (leg.tail[1] + leg.tip[1]) / 2,
        (leg.tail[2] + leg.tip[2]) / 2,
    )
    offset = if leg.kind == :virtual
        leg.direction == :east  ? (0.0, -0.11) :
        leg.direction == :north ? (-0.11, 0.0) :
        leg.direction == :west  ? (0.0, 0.11) : (-0.11, 0.0)
    else
        leg.direction == :east  ? (0.07, -0.01) :
        leg.direction == :north ? (0.02, 0.07) :
        leg.direction == :west  ? (-0.07, 0.03) : (0.05, -0.06)
    end
    return (midpoint[1] + offset[1], midpoint[2] + offset[2])
end
```

In `plot_peps_graph`, select `leg.site == (rows, 1)` and draw each
`leg.symbol` once:

```julia
for leg in filter(leg -> leg.site == (rows, 1), legs)
    label_x, label_y = _peps_symbol_label_position(leg)
    CairoMakie.text!(axis, label_x, label_y;
        text=leg.symbol, fontsize=14, color=:gray15,
        align=(:center, :center))
end
```

For `show_index_labels=true`, retain the current compact diagnostic labels for
physical legs. Recover the actual index with
`physicalinds(peps, row, col)[direction_number]`, where direction numbers are
E=1, N=2, W=3, S=4. Use
`direction_labels = ("E", "N", "W", "S")` and place
`_compact_index_label(index, direction_labels[direction_number])` near the
physical arrow tip. Do not add fabricated identifiers to conceptual
open-boundary virtual stubs.

Expand PEPS plot bounds to include open stubs and physical arrow tips:

```julia
_finish_graph_axis!(axis, 0.12, cols + 0.88, 0.12, rows + 0.88)
```

Increase the figure size only if the 3×3 visual check shows collisions; keep
the default compact enough for a paper column.

- [ ] **Step 5: Update the README figure description**

Add this paragraph next to the existing `plot_peps_graph` example:

```markdown
The PEPS figure is an open-boundary tensor diagram: every `T` has four virtual
legs and four physical legs. Horizontal virtual arrows run
`γ → T → α`, vertical arrows run `δ → T → β`, and the lower-left tensor
provides the representative `α,β,γ,δ,i,j,k,l` labels. This is a diagrammatic
view; the finite PEPS data continues to contract its boundary virtual indices
with `|+⟩` vectors.
```

- [ ] **Step 6: Run focused rendering/export tests**

Run the focused command from Step 2.

Expected: geometry and figure testsets pass; temporary SVG/PDF files are
nonempty and include the required text plots.

- [ ] **Step 7: Render and visually inspect 2×2 and 3×3 PNG previews**

Run:

```bash
julia --project=. -e 'using TopoNoise, CairoMakie; for n in (2,3); p=toric_code_peps(n,n); save("/tmp/toponoise-open-peps-$(n)x$(n).png", plot_peps_graph(p)); end'
```

Inspect both previews. Confirm all boundary stubs are visible, internal bonds
are continuous, no arrowheads cover tensor labels, all four physical arrows
are distinguishable, and only `(R,1)` has symbolic leg labels. Adjust only
geometry constants and label offsets if needed, then rerun Step 6.

- [ ] **Step 8: Commit the renderer and documentation**

```bash
git add src/visualization.jl test/visualization.jl README.md
git commit -m "feat: render open-boundary PEPS isometry"
```

---

### Task 3: Full regression verification

**Files:**
- Verify only; no production changes expected.

**Interfaces:**
- Consumes: the completed renderer and the unchanged package test entrypoint.
- Produces: fresh evidence that the figure-only change did not affect PEPS or circuit behavior.

- [ ] **Step 1: Run the complete direct test suite**

```bash
julia --project=. test/runtests.jl
```

Expected: every existing tensor, PEPS, graph, schedule, Yao state, probability,
and vector-figure testset passes with zero failures.

- [ ] **Step 2: Run isolated package tests**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: `Testing TopoNoise tests passed` and exit code zero.

- [ ] **Step 3: Check scope and whitespace**

```bash
git diff --check
git status --short
```

Expected: no whitespace errors. Only the planned renderer, tests,
documentation, and pre-existing unrelated dirty-worktree files may differ.

- [ ] **Step 4: Review the requirement checklist**

Confirm from code, tests, and previews:

- `8RC` total leg records;
- `4RC` virtual half-legs and `4RC` physical legs;
- `2R+2C` boundary virtual stubs;
- east/north virtual orientation;
- northeast physical orientation;
- labels only at `(R,1)`;
- unchanged public signature and valid CairoMakie figure return;
- no edits to PEPS/circuit production files;
- no global circuit matrix construction.
