# `document/test.typ` Full Refresh Design

## Goal

Refresh `document/test.typ` so it accurately documents the complete current
TopoNoise.jl test suite, including the noisy measured-trajectory work. The
document will remain a readable technical explanation rather than becoming a
raw inventory of test source lines.

## Authoritative Evidence

The assertion counts and pass status will come from the final successful
`julia --project=. -e 'using Pkg; Pkg.test()'` run on commit `13f5c72` plus
the reviewed crossing fixes in the current working tree. That run reports 837
assertions across 23 named test groups. The final document will not claim tests
that were not executed.

## Organization

The refreshed Typst document will be organized by subsystem:

1. Overview and a 23-row summary table with 837 total assertions.
2. Local tensor, isometry, unitary completion, and local-gate tests.
3. Finite PEPS structure, exact small-patch wavefunctions, and the memory-safe
   `3 x 3` verification.
4. PEPS graph, sequential schedule, executable Yao circuit, and gate-convention
   tests.
5. Complete virtual-bond noise model, including all six array shapes, total
   bond count, independent sampling, and boundary metadata.
6. Literal Yao measurement/reset trajectories and the scalable parity-rule
   sampler.
7. Exact `1 x 1`, `1 x 2`, and `2 x 1` classical-versus-Yao distributions,
   including every dangling boundary `X` location.
8. Mismatch, frustration, cluster, and spanning observables.
9. Welford scans and crossing tests, including `:ok`, `:unbracketed`, and
   `:unstable` cases, endpoint contacts, plateaus, and batch validation.
10. Runnable examples, CLI CSV and figure artifacts, native Yao figures, and
    CairoMakie vector figures.
11. Commands and environment used to run the suite.

## Content Rules

- Preserve useful mathematical derivations from the existing document while
  removing obsolete totals and claims.
- Use `(E,N,W,S)` consistently for physical records and describe horizontal
  west-to-east and vertical south-to-north carrier orientation.
- State the internal mismatch identities with the same row convention used by
  the implementation.
- Explain that boundary errors are sampled and tested but remain physically
  record-invariant for `|+>` input and `<+|` output boundaries.
- Distinguish raw internal-bond percolation tests from Binder-cumulant or
  Nishimori-decoding analyses, which are not implemented.
- Leave `document/note.tex` unchanged.

## Verification

After editing, run:

1. `typst compile document/test.typ` when the Typst CLI is available.
2. Inspect the generated PDF or rendered pages for table overflow, broken
   equations, and pagination problems.
3. Recalculate the table total and require it to equal 837.
4. Run `git diff --check` and confirm `document/note.tex` has no diff.

## Acceptance Criteria

The refresh is complete when `document/test.typ` is internally consistent,
contains all 23 passing groups and their correct counts, explains the new
trajectory tests at the same technical depth as the existing tensor tests,
compiles cleanly, renders legibly, and makes no unsupported scientific claim.
