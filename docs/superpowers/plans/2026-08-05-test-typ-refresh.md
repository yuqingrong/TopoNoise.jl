# `document/test.typ` Full Refresh Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the obsolete test summary with an accurate, technically detailed Typst document covering all 837 passing assertions in 23 test groups.

**Architecture:** Keep one self-contained Typst document organized by subsystem. Preserve useful tensor and PEPS derivations, then add circuit, noisy-trajectory, statistical-analysis, CLI, and rendering verification sections whose claims map directly to the final `Pkg.test()` output.

**Tech Stack:** Typst, Julia 1.12.6, Test, Yao 0.9.3, ITensors 0.9.30, ITensorNetworks 0.21.5, CairoMakie 0.15.13, and PDF rendering for visual inspection.

## Global Constraints

- State exactly 837 assertions across 23 named groups.
- Use `(E,N,W,S)` physical-record order throughout.
- Horizontal carriers run west to east; vertical carriers run south to north.
- Boundary `X` events are sampled metadata but record-invariant for `|+>` and `<+|` boundaries.
- Do not describe raw percolation observables as Binder-cumulant or Nishimori-decoding results.
- Preserve useful existing derivations while removing obsolete totals.
- Leave `document/note.tex` unchanged.

---

### Task 1: Refresh and verify the complete Typst test report

**Files:**
- Modify: `document/test.typ`
- Verify only: `document/test.pdf`
- Preserve: `document/note.tex`

**Interfaces:**
- Consumes: the final `Pkg.test()` group names and counts below.
- Produces: a compiling `document/test.typ` and a legible PDF describing the complete verified suite.

- [ ] **Step 1: Replace the obsolete overview and summary table**

Use this authoritative inventory exactly:

```text
Local toric-code tensor                          262
Finite PEPS structure                            25
Complete small-patch wavefunctions                9
3 x 3 patch                                       6
Unitary completion                               10
Toric-code local unitary gate                    19
Toric-code PEPS graph                            49
Toric-code sequential circuit schedule          122
Executable Yao circuit                           26
Yao local-gate convention and small-patch states 11
Exact small-patch trajectory distributions       25
Virtual-bond trajectory model                    29
Literal Yao measured trajectories                 8
Scalable measured trajectories                   10
Trajectory mismatch observables                  29
Welford streaming accumulator                     3
Streaming trajectory scans                       45
Finite-size spanning crossings                   26
Runnable toric-code circuit example              36
Trajectory critical-scan CLI                     36
Trajectory scan figure                           14
Native Yao circuit figures                       13
Toric-code vector figures                        24
TOTAL                                            837
```

Use a compact four-column Typst table so 23 rows remain readable:

```typst
#table(
  columns: (2.5fr, 0.55fr, 2.5fr, 0.55fr),
  table.header([*Test group*], [*Count*], [*Test group*], [*Count*]),
  // Pair the 23 groups across two halves and place TOTAL in the last row.
)
```

- [ ] **Step 2: Retain and update tensor, gate, and PEPS explanations**

Keep the equations for the normalized local tensor, parity matrix, isometric
cuts, unitary completion, finite-patch configuration count, and memory-safe
`3 x 3` norm. Update prose so these areas no longer claim to be the whole
suite. Retain the tensor definition:

```typst
$
T^(i j k l)_(alpha beta gamma delta)
= cases(
  1 / sqrt(2) delta_(i alpha) delta_(j beta)
    delta_(k gamma) delta_(l delta)
    & alpha + beta + gamma + delta = 0 quad "mod" quad 2,
  0 & "otherwise",
)
$
```

- [ ] **Step 3: Add graph and sequential-circuit verification**

Cover these executed checks:

```text
PEPS graph: site keys, internal graph edges, retained physical legs, invalid input.
Schedule: R+C-1 diagonal layers, lower-left to upper-right order, disjoint layer wires.
Yao circuit: carrier H preparation, local subroutines, qubit counts, state agreement.
Convention: local Yao blocks match the tensor orientation on small patches.
Figures: native Yao and CairoMakie SVG/PDF/PNG artifacts are nonempty and valid.
```

State that full circuit matrices are not materialized.

- [ ] **Step 4: Add complete virtual-bond and trajectory verification**

Document all six arrays and the total count:

```typst
$N_("bond") = R(C-1) + (R-1)C + 2R + 2C.$
```

Explain the tested timeline: prepare carriers in `|+>`, apply west/south
boundary `X`, reset four physical ancillas, apply the local gate, measure and
reset `(E,N,W,S)`, apply eligible internal `X` before the downstream site,
then apply east/north `X` after the final site before analytic `<+|`.

Include the scalable local rule:

```typst
$beta = alpha xor gamma xor delta,$
$(i,j,k,l) = (alpha,beta,gamma,delta).$
```

Describe the exact `1 x 1`, `1 x 2`, and `2 x 1` classical-versus-Yao
enumeration, six independent noise-group checks, every shape validator, every
dangling boundary location, the 24-qubit Yao guard, ancilla reset, and
fixed-seed boundary invariance.

- [ ] **Step 5: Add observables, scans, crossings, CLI, and artifact checks**

Use the implementation's vertical-array convention:

```typst
$m^h_(r,c) = E_(r,c) xor W_(r,c+1),$
$m^v_(r,c) = N_(r+1,c) xor S_(r,c).$
```

Describe exact mismatch/error equality, total and boundary densities,
odd-plaquette frustration, union-find largest cluster, both spanning flags,
and the tested `p=0` and `p=1` limits.

Explain Welford means/standard errors, seeded scans, PAVA smoothing, batch
bootstrap intervals, the two-batch minimum, and all crossing statuses. State
that endpoint-only contacts are `:unbracketed`; identical or sign-changing
equal plateaus and low-valid-bootstrap estimates are `:unstable`.

List the CLI artifacts exactly:

```text
trajectory_scan.csv
trajectory_crossings.csv
trajectory_scan.svg
trajectory_scan.pdf
trajectory_scan.png
```

Mention validation of required `--p-min`, `--p-max`, and `--p-step`, seeded
tiny integration runs, CSV headers/rows, and file signatures/nonempty output.

- [ ] **Step 6: Add scope and execution notes**

End with this scope callout:

```typst
#block(fill: rgb("f5f1e8"), inset: 8pt)[
  These trajectory observables test raw percolation of sampled internal error
  bonds. They do not reconstruct spins, compute a Binder cumulant, or test a
  Nishimori decoder.
]
```

Give the verified command:

```text
julia --project=. -e 'using Pkg; Pkg.test()'
```

- [ ] **Step 7: Compile the Typst source**

Run:

```bash
typst compile document/test.typ document/test.pdf
```

Expected: exit code 0 with no Typst parse, layout, or font error.

- [ ] **Step 8: Render and visually inspect the PDF**

Render every page to PNG using the PDF workflow. Inspect table width, row
wrapping, equation clipping, page breaks, headings, and monospace text. If a
page is crowded or clipped, adjust table inset/font size or pagination,
recompile, and inspect again.

- [ ] **Step 9: Verify counts and repository hygiene**

Re-add the 23 counts from Step 1 and require 837. Run:

```bash
git diff --check
git diff --exit-code -- document/note.tex
git status --short
```

Expected: no whitespace errors, no `note.tex` diff, and only the intended
Typst/planning changes.

- [ ] **Step 10: Commit the refreshed report**

```bash
git add document/test.typ
git commit -m "docs: refresh complete test report"
```
