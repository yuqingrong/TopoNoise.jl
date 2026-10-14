# H+5 Z-boundary control verification

**Status:** PASS, diagnostic evidence only. No deterministic Z threshold is claimed.

**Verified implementation head before this record:** `b0d886c7dbb375246dbde439069c097a8b68fa4e`

## Complete regressions

The complete Julia regression command was:

```text
env JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 julia --project=. --startup-file=no test/runtests.jl
```

Exit code: `0`. The retained log contains 151 test summaries, 46,390 passing assertions, zero failures, zero errors, and zero broken tests. Timing was 540.09 s wall, 487.23 s user, and 31.16 s system.

The complete Python regression command was:

```text
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .CondaPkg/.pixi/envs/default/bin/python -m unittest discover -s test/toric3d -p 'test_*.py' -v
```

Exit code: `0`. It ran 56 tests in 25.668 s and reported `OK`; wall time was 27.33 s.

## New-artifact audit and fit guard

The command

```text
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .CondaPkg/.pixi/envs/default/bin/python results/toric3d-h5-z-boundary-ambiguity-2026-10-14/apply_guard.py
```

exited `0` after independently checking 12 audit identities, 1,328 observable-group rows, and 24 paired control rows. Every protocol was rejected before fitting for nonzero mixed groups, unavoidable faults, and control failures. The run records 1,000 requested bootstrap replicates and 0 executed replicates.

## Historical-delivery immutability

Running the older delivery verifier directly at the new implementation head exits `1` at its assertion that current source bytes equal the run's frozen source hashes. This is expected after extending the source tree. A blob-level comparison showed:

- all frozen source hashes match commit `4abe313` exactly;
- the current head differs only in `src/TopoNoise.jl` and `src/toric3d/vertex_projector_replay.jl`, both intentionally modified by this plan; and
- `git diff --name-only cb2ec97..HEAD -- document results/toric3d-h5-diagonal-feedback-2026-10-13` is empty.

The unmodified verifier was therefore run in a temporary detached worktree at `4abe313`, using the same pinned Python environment:

```text
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /Users/rongyuqing/jcode/TopoNoise.jl/.CondaPkg/.pixi/envs/default/bin/python results/toric3d-h5-diagonal-feedback-2026-10-13/preliminary-l3-7/verify_delivery.py
```

Exit code: `0`. It reported:

```text
PASS: final measurements, guards, provenance, 45610 Julia assertions and 51 Python tests verified
```

That check revalidated 195 points, 195,000 paired records, 780,000 arm replays, 186,000 exactly reused records, zero invalid recoveries, checkpoint ancestry, quoted measurements, figures, and frozen source provenance. The temporary worktree was removed afterward; no old result file was changed in this branch.

## Immutable note and branch checks

`git diff --check` exited `0`. The required note hash is unchanged:

```text
3a40710c2ab6bd43e53f5f7703b7c8716b209dc1e060a5092f17e0edce073f3a  document/note.tex
```

Neither `document/note.tex`, compiled note artifacts, nor any prior H+5 result artifact appears in the branch diff. The worktree was clean before creating this verification record.

## New artifact hashes

| Artifact | SHA-256 |
|:---|:---|
| `README.md` | `656154ff4e134596697f95781a2fedb53ca2e0ef8bd5189e57888c688e01d101` |
| `apply_guard.py` | `c4df89e380c36b8a155908a956354e393224bd01ec7d65a05282077770b21d89` |
| `guard.json` | `dc9b3ffb413c0b8fb1ace50366ede957dc05b04ebbaff60089dbc03bc2fb082c` |
| `parity-erasure-control.csv` | `98fcbab9714a90fb5be8e486b6db58e4ef5320bdb7ab769c4572abc0be39b468` |
| `run_audit.jl` | `704509e19d963353460418c98a1a8fa35a9bd8e79ff27bf692b4a46ae9a28327` |
| `run_control.jl` | `ca47b7d744fa1976b3a37adfe0b812026e442cb0ce1b9b27d9c6188181af6d04` |
| `runtime-environment.json` | `fca01ad566032192f5ff50fc60cdcc122c5b4735887911f585747dda0de32218` |
| `z-boundary-audit.json` | `b11e3f8e3779c0f6419ee95726f6d5bfcf48a77003eaa3558989ffbaed0e120b` |
| `z-boundary-groups.csv` | `f7c7b322cb40972b4308ef94b26ed6db1df6c7619fa8e2bba50be511dc1328f7` |
| `z-boundary-report.md` | `9023e9908a2e34b3bf7e41d068d81f4a6a972738916f07fe1c284ce35fc4a255` |

No push, merge, publication, remote job, or large threshold scan was performed.
