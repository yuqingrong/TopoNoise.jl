# Toric-Code Isometry Orientation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Express the normalized local tensor as a `64 × 4` isometry and test the requested identity `T†T = I₄`.

**Architecture:** The dense rank-eight tensor remains unchanged. Only its matrix interpretation changes: `(α, β)` is the four-dimensional input, while `(γ, δ, i, j, k, l)` is the 64-dimensional output. Tests and explanatory documentation will use this orientation consistently.

**Tech Stack:** Julia 1.12.6, Test stdlib, ITensors.jl 0.9.30, Typst 0.14.2.

## Global Constraints

- Preserve stored tensor order `(pE, pN, pW, pS, vE, vN, vW, vS) = (i, j, k, l, α, β, γ, δ)`.
- Preserve all normalized local tensor entries and every finite-PEPS amplitude.
- Define the matrix as `T[(γ,δ,i,j,k,l), (α,β)] ∈ ℂ^(64×4)`.
- Test exactly `adjoint(T) * T ≈ I₄`; do not retain `T * adjoint(T) ≈ I₄` as the isometry assertion.
- Keep the local test count at 262 and the complete suite at 302 assertions.
- Introduce no new exported API or dependency.

---

### Task 1: Reorient the local isometry test and documentation

**Files:**
- Modify: `test/toric_code_peps.jl:94-103`
- Modify: `document/test.typ:37,70,93-106`
- Verify: `docs/superpowers/specs/2026-08-03-toric-code-local-isometry-design.md`

**Interfaces:**
- Consumes: `toric_code_local_tensor()::Array{Float64,8}` in stored E–N–W–S order.
- Produces: a test-local `isometry::Matrix{Float64}` with shape `(64, 4)` and the documented identity `T†T = I₄`.

- [x] **Step 1: Replace the co-isometric matrix interpretation in the local test**

```julia
# Stored order: (pE, pN, pW, pS, vE, vN, vW, vS).
# Rows are (γ, δ, i, j, k, l) and columns are (α, β).
isometry = reshape(permutedims(tensor, (7, 8, 1, 2, 3, 4, 5, 6)), 64, 4)
identity_αβ = [
    1.0 0.0 0.0 0.0
    0.0 1.0 0.0 0.0
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 1.0
]
@test adjoint(isometry) * isometry ≈ identity_αβ
```

- [x] **Step 2: Run the focused test entrypoint**

Run:

```sh
julia --project=. test/runtests.jl
```

Expected: `Local toric-code tensor | 262 262`, followed by all other test groups passing. The assertion is expected to pass immediately because it reorients an already normalized mathematical tensor; no production-code change is required.

- [x] **Step 3: Update the Typst report to the isometric orientation**

Replace “co-isometry” with “isometry,” describe the permutation as
`(γ,δ,i,j,k,l,α,β)`, and document

```text
T_((gamma delta i j k l),(alpha beta)) in CC^(64 times 4)
T^dagger T = I_4.
```

State that each `(α,β)` column contains two entries of magnitude `1/sqrt(2)` and distinct columns have disjoint support.

- [x] **Step 4: Compile and visually verify the Typst report**

Run:

```sh
typst compile document/test.typ /tmp/toponoise-test-summary.pdf
pdfinfo /tmp/toponoise-test-summary.pdf
```

Expected: compilation exits successfully and the report remains a readable three-page A4 document without overflow or orphaned content.

- [x] **Step 5: Run the complete package test suite**

Run:

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected test summaries: local tensor `262/262`, finite structure `25/25`, small patches `9/9`, and `3×3` patch `6/6`, totaling 302 assertions.

- [x] **Step 6: Check consistency and whitespace**

Run:

```sh
rg -n 'local co-isometry|reshaped into a \$4 times 64\$|M M\^dagger' README.md src test document docs
git diff --check
```

Expected: no stale active co-isometry description in tests or user documentation, and no whitespace errors. The design specification may mention the rejected `4×64` convention only to explain why it is out of scope.
