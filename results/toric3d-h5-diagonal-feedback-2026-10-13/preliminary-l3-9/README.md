# Preliminary H+5 run: cubic L=3,5,7,9

**Deferred:** The user subsequently confirmed L=3,5,7 only. This run was stopped at206000 saved records. Continue `../preliminary-l3-7/` instead; do not resume this four-size run automatically.

On 2026-10-14 the user requested smaller lattices first because the largest
sizes were taking substantially longer than the scout timing projection.
This is an explicit amendment to the six-size selection, motivated only by
runtime. All original rates, candidate fit windows, protocols, decoder
settings, source hashes, 1000 paired records per point and seed1234 remain.
Results must be labeled preliminary and retain every existing fit guard.

The original `../selection.json` and `../production/` checkpoint are preserved.
That run was deliberately interrupted at106350 saved paired records with no
invalid recovery. L11 and L13 are deferred and must not resume automatically.

The amended selection hash is
`dca614dc243ec8d918cb91ac27c74425c7af0fbba6f2381c88d33d3ebf025412`.
There are80 X points and180 Z points, totaling260000 paired records and1040000
protocol replays. All80 X points already existed as a complete prefix.

`prepare_checkpoint.jl` copies their entire states, including all23 joint-count
arrays, timings and record/signature digests. It generates the fresh260-point
seed stream from `Random.Xoshiro(1234)` and verifies that every reused X seed
matches exactly. Unstarted Z points use this fresh reduced-grid seed stream;
their seeds need not equal the deferred six-size run's Z seeds.
`checkpoint-reuse.toml` records source and initial checkpoint hashes, indices,
and reuse assertions. An independent Python comparison also verified these
states and the configuration transformation before launch.

Run from the repository root with the numerical environment in
`../runtime-environment.json`:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null \
JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" \
MPLCONFIGDIR=/tmp/toponoise-mpl \
julia --project=. examples/compare_toric3d_vertex_projector_feedback.jl \
  --stage production \
  --selection results/toric3d-h5-diagonal-feedback-2026-10-13/preliminary-l3-9/selection.json \
  --shots 1000 --seed 1234 \
  --output results/toric3d-h5-diagonal-feedback-2026-10-13/preliminary-l3-9/production \
  --resume
```

For a completely fresh reproduction, choose a new empty output directory and
omit `--resume`; it recomputes all260 points using the same seed assignment.
Do not rerun the preparation script over an existing preliminary checkpoint.
The original full scan is not complete, and this subset must never be reported
as completion of that original frozen grid.

After this run completes, independently audit `production/`, run the existing
analysis CLI with this selection and `--bootstrap 1000 --seed 1234`, and write
outputs to `analysis/`. Audit the reused X states against their original states
again, inspect all figures, then finish the report, note and final regressions.
