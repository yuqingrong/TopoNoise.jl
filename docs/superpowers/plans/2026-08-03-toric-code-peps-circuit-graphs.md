# Toric-Code PEPS and Sequential-Circuit Graphs Implementation Plan

> **For agentic workers:** Follow test-driven development task by task. Do not
> materialize the global circuit matrix.

**Goal:** Add tensor-network graph conversion, an exact postselected diagonal
circuit schedule, Yao integration, and static PEPS/circuit figures.

**Architecture:** `ITensorNetworks` owns PEPS topology, a focused TopoNoise
model owns schedule and boundary semantics, `Yao` owns the unitary core, and
`CairoMakie` renders both approved fixed layouts.

**Tech Stack:** Julia 1.10+, ITensors 0.9, ITensorNetworks 0.21, Yao 0.9,
CairoMakie 0.15, Test stdlib.

## Global Constraints

- Accept only `ToricCodePEPS` in version one.
- Preserve E/N/W/S physical-index and gate ordering.
- Use layer `rows - row + col`, sorting rows descending within each layer.
- Use `4rows*cols + rows + cols` circuit wires.
- Represent `|+⟩` input and `⟨+|` postselection explicitly.
- Return figures; use CairoMakie's standard `save` for export.

## Tasks

### 1. Dependencies and PEPS graph

- Add direct package compat entries and module registration.
- Write failing topology and malformed-input tests.
- Implement `peps_graph` with coordinate keys and rectangular validation.
- Run focused and existing tests.

### 2. Sequential schedule and Yao core

- Write failing tests for wire allocation, exact layers, success exponent, and
  local gate embedding.
- Implement `ToricCodeSequentialCircuit`, its accessors, and Yao subroutines.
- Test small exact postselected states without constructing a global matrix.

### 3. Static visualization

- Write failing figure/export tests.
- Implement the fixed PEPS layout and approved grouped-wire circuit layout in
  `src/visualization.jl`.
- Save temporary SVG/PDF files and verify labels and nonempty output.

### 4. Documentation and final verification

- Add the `3×3` generation/export example to the README.
- Run focused tests, `test/runtests.jl`, `Pkg.test()`, and `git diff --check`.
- Review the final diff against every public API and design invariant.
