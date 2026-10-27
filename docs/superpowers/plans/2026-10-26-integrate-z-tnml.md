# Integrate the validated Z-only TN decoder

> **For agentic workers:** Use the approved integration contract below. Tests, documentation, and implementation have separate file ownership; request one independent final code review.

**Goal:** Make the validated diagonal H+6CNOT, no-feedback, Z-only approximate ML decoder available through `using TopoNoise` in the user's main checkout.

**Architecture:** Extract the exact propagated-fault topology helpers from the research worktree, then transplant the validated conditioned-star MPO / canonical boundary-MPS algorithm into package source with package-specific names. Preserve the numerical operations and signed output weights. Keep the original experimental source and recorded benchmarks intact.

**Tech stack:** Julia, existing ITensors, explicit ITensorMPS dependency.

## Approved scope and implementation decisions

The user approved moving the validated implementation, its dependencies, API, and tests into the main package. Work directly in `/Users/rongyuqing/jcode/TopoNoise.jl`; preserve existing document edits and unrelated files. No branch switch, commit, new worktree, benchmark resampling, or threshold scan is required. This integration request already authorizes the previously discussed design; no additional approval round is needed.

Public API:

- `Toric3DZTNTopology`, `build_toric3d_z_tn_topology(model)`, `toric3d_z_tn_edge_probabilities(topology, p)`.
- `Toric3DZTNSettings(; max_bond_dimension=240, cutoff=1e-18, negative_tolerance=1e-10)`; reject invalid settings.
- `Toric3DZTNDecoder`, `build_toric3d_z_tn_decoder(code_or_model_or_topology, p; settings=...)`.
- `decode_toric3d_z_tn(decoder, syndrome)` returns the validated named-tuple fields: `valid`, `probabilities`, `logical_mask`, `log_syndrome`, `raw_probabilities`, and `diagnostics`.

Logical masks are `ell1 + 2ell2`, so vector entries correspond to `(0,0),(1,0),(0,1),(1,1)`. A logical-class decision is not a physical correction vector. `valid` reports numerical admissibility, not convergence. Local SVD truncation diagnostics do not bound posterior error. Preserve the existing physical-rate domain `0 <= p < 0.5`. Defaults were validated for L3 at p=.006,.01,.02; larger sizes and other noise models require separate validation.

Files:

- `src/toric3d/z_tnml_topology.jl`: exact topology/provenance construction and effective-edge probabilities only.
- `src/toric3d/z_tnml.jl`: validated MPS algorithm, public build/decode functions, settings and documentation.
- `src/TopoNoise.jl`, `Project.toml`, `Manifest.toml`: dependency and package wiring.
- `test/toric3d/z_tnml.jl`, `test/toric3d/fixtures/z_tnml_l3.toml`, test runner include: enumerated small examples, input/error behavior, topology mapping, saved L3 regression.
- `docs/src/z_tnml.md`, docs navigation, README link, `examples/toric3d_z_tnml.jl`: executable usage and scope.

## Tasks

- [x] Record starting status and run baseline/package availability check.
- [x] Write package-facing tests; observe missing API failure before source integration.
- [x] Add exact topology helpers and transplant the validated contraction with explicit settings.
- [x] Add dependency/import/export wiring and resolve the manifest without unrelated upgrades.
- [x] Run focused tests, then compare integrated output with frozen L3 cases.
- [x] Add and run the user-facing example; finish docs.
- [x] Run the project's complete test suite; report any failure by name.
- [x] Independent review, fix material findings, verify final diff and report completion.

## Review focus

Preserve numerical equivalence to research module SHA256 `0321bc105503ccc642ecc8c33bbf9627e6d50f816b715a4dc24b89ee30db1f53`. Do not import earlier Fourier reconstruction, top-K pruning, mdopt, or an exact fallback. Verify complete detector-plus-logical signature merging, class bit order, zero-probability syndrome behavior, signed probability handling, and no hidden dependency on a worktree or temporary path. Ensure `using TopoNoise` suffices. Package entry points must not set global BLAS threads. Exercise the actual public path, not only internal constructors.

## Execution record

Starting main checkout has unrelated modifications to `.DS_Store` and `document/note.*`, plus existing untracked research notes/results. No source or dependency files were modified before this integration.

Baseline `using TopoNoise` loaded this main checkout successfully. Before implementation, the new focused tests failed on 14 missing API/export assertions with no test errors. Offline dependency resolution added ITensorMPS 0.3.45. The repository ignores local `Manifest.toml`; the explicit dependency and compatibility range are recorded in `Project.toml`.

Verification uses Julia 1.12.6 with one Julia and BLAS thread, a writable temporary depot layered over the user's depot, and the existing local Python environment. These are test-process settings; the decoder does not change global thread settings.

- Focused tests: `using Test, TopoNoise; include("test/toric3d/z_tnml.jl")` passed all 629 assertions, including independent small-graph enumeration and both frozen L3 cases. The L3 posterior regression tolerance is `1e-10`, with a separate posterior-L1 check against the exact reference at `1e-3`.
- Public example: `julia --project=. examples/toric3d_z_tnml.jl` completed with an admissible result, logical mask 1, and actual maximum bond dimension 186 for its deterministic L3 syndrome at p=.01.
- Independent read-only review found no critical or important issues. It verified normalized numerical equivalence to the reference, unchanged topology helpers, fixture hashes and values, documentation statistics, and package wiring. One optional future test could explicitly exercise a known negative-weight/invalid result; existing tests check the signed-output validity rule and the unchanged implementation preserves that behavior.
- Local documentation build passed, including Documenter doctests and page checks, with deployment excluded. Run from `docs/`; an initial invocation from the repository root had the wrong page root and was corrected. Its generated duplicate source directory was verified against `src/` and removed. The generated all-API index emitted a nonfatal size warning (182.18 KiB, below the 200 KiB failure limit).
- Full package validation: `using Pkg; Pkg.test(; allow_reresolve=false)` exited successfully. All 47,142 assertions in 163 test sets passed, including the 629 new decoder assertions and both L3 posterior regressions. Cold dependency compilation dominated the initial wait; the suite subsequently completed normally. Matplotlib used its writable temporary cache fallback.
- Final whitespace check passed. The final working-tree status contains the intended integration files and the same unrelated pre-existing notes/results edits. No commits, publication, threshold scan, or new held-out benchmark was performed.

Verification logs for this run are in `/private/tmp/toponoise-tnml-integration-20261026/`: `focused-red.log`, `focused-green.log`, `example.log`, `full-test.log`, and `docs-build.log`. Permanent regression inputs are in `test/toric3d/fixtures/z_tnml_l3.toml`; no temporary or worktree path is needed to run the package decoder or tests.
