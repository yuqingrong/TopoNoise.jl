# Preliminary H+5 diagnostics: L=3,5,7

The user explicitly confirmed these three sizes on2026-10-14, clarifying an
earlier L=4 typo. L9,11,13 are deferred. Both prior runs are deliberately stopped
with their selections/checkpoints preserved. No scientific source code, rate,
fit window, decoder setting or four-size threshold guard has changed.

This run has195 points and195000 paired records (780000 protocol replays).
Its derived checkpoint preserves186000 existing records exactly, leaving9000
Z-only L7 records at launch. The selection SHA-256 is
`1597b9941132d3f24c11c68ec1933ea639f5289225fe80ff5070367e81416e72`.
Its status is explicitly `diagnostic_only`, and metadata stage is
`preliminary_three_size`.

The seed assignment is inherited from the parent260-point run: generate260
UInt64 seeds with `Random.Xoshiro(1234)`, then take the parent indices listed in
`production/checkpoint-reuse.toml`. Existing records are never reseeded.
`prepare_checkpoint.jl` copies every retained point's complete state, including
all23 joint-count arrays, exact point seeds, timings and record/signature
digests. An independent Python comparison verified the projection before
launch. The parent checkpoint SHA-256 is
`c3544f0d3f72ff2cd5fe0e79379ec502b0360cdc15e76a3f7b9cff4a18faa80c`.

Resume from the repository root (the runner always resumes its checkpoint):

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null \
JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" \
MPLCONFIGDIR=/tmp/toponoise-mpl \
julia --project=. results/toric3d-h5-diagonal-feedback-2026-10-13/preliminary-l3-7/run_preliminary.jl
```

For a fresh reproduction, run `prepare_checkpoint.jl /new/empty/output --fresh`
and pass that same output directory to `run_preliminary.jl`. This resets counts
while preserving the documented parent-grid seed assignment. The preserved
parent checkpoint is required. Never prepare over an existing checkpoint.

After completion, run the independent `../audit_artifacts.py` on `production/`
and audit all already-complete retained states against the parent checkpoint.
Then run `analyze_preliminary.py` with the same one-thread numerical environment.
It calls the unchanged fitting guard with1000 requested bootstrap replicates:
all candidates fail the minimum-four-size check before fitting or bootstrap.
Report fixed-size paired effects and preliminary curves, with no threshold,
scaling exponent or fitted confidence interval. Empty bootstrap/sensitivity
outputs and blank scaling panels must be explained as not performed.

The final report, LaTeX note, PDF, Julia/Python regressions and local data commit
remain part of delivery. The original six-size and four-size runs are deferred,
not completed by this subset.
