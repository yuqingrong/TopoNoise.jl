# H+6 final-only Z aspect-ratio comparison

This study compares `T=L` and `T=2L` for the encoded diagonal
plaquette-transfer H+6-CNOT circuit.  It uses Z-only faults after every active
H/CNOT output, the noisy measurement-prepared logical `|++>` seed, no
intermediate feedback, and one validated unit-data-weight MWPM decode at the
end.

The run completed all 30 requested points: `L=3,5,7`, both aspect ratios,
five probabilities from `0.0045` through `0.0065`, and 1,000 independent
samples per point (30,000 fault records total).  Every recovery passed the
core syndrome validation.  The independent artifact audit passed every grid,
count, metadata, input-hash, and output-presence check.

## Result

At every measured `(L,p)` point, `T=2L` had a larger final logical-Z failure
rate than `T=L`.  Representative values are:

| p | L | T=L | T=2L | Difference (`2L-L`) |
|---:|---:|---:|---:|---:|
| 0.0050 | 3 | 1.6% | 5.1% | +3.5 points |
| 0.0050 | 5 | 1.6% | 4.5% | +2.9 points |
| 0.0050 | 7 | 1.0% | 4.2% | +3.2 points |
| 0.0065 | 3 | 3.0% | 7.8% | +4.8 points |
| 0.0065 | 5 | 3.1% | 10.9% | +7.8 points |
| 0.0065 | 7 | 3.9% | 15.2% | +11.3 points |

The conservative envelope formed from the two separate 95% Wilson intervals
is above zero at 14 of the 15 comparisons.  At `L=7,p=0.0045`, it is
`[-0.296, 2.488]` percentage points and therefore includes zero.  These are
unpaired comparisons; the envelope is descriptive and is not presented as a
multiple-comparison-adjusted confidence interval.

This behavior is expected geometrically.  Changing `T=L` to `T=2L` doubles
the number of noisy transfer cells and increases
`dX=(T+1)L`, but leaves the relevant logical-Z distance `dZ=2L` unchanged.
The small-size data therefore provide no evidence that `T=2L` repairs the Z
sector.  They show worse finite-size final-only performance in this sampled
window.

No threshold or scaling exponent was fitted.  Three small transverse sizes
are insufficient for a guarded threshold claim, and the two aspect-ratio
families should not be combined in one scaling fit.  The result is specific
to the implemented circuit-correlated Z noise and unit-weight MWPM decoder.

## Artifacts

- [`production/toric3d-z-threshold-raw.csv`](production/toric3d-z-threshold-raw.csv): raw counts.
- [`production/toric3d-z-threshold-metadata.toml`](production/toric3d-z-threshold-metadata.toml): core simulation metadata.
- [`production/toric3d-z-aspect-ratio-study.toml`](production/toric3d-z-aspect-ratio-study.toml): comparison scope.
- [`analysis/aspect-ratio-points.csv`](analysis/aspect-ratio-points.csv): rates and Wilson intervals.
- [`analysis/aspect-ratio-differences.csv`](analysis/aspect-ratio-differences.csv): fixed-size unpaired contrasts.
- [`analysis/aspect-ratio-summary.json`](analysis/aspect-ratio-summary.json): analysis provenance.
- [`analysis/aspect-ratio-comparison.png`](analysis/aspect-ratio-comparison.png),
  [`PDF`](analysis/aspect-ratio-comparison.pdf), and
  [`SVG`](analysis/aspect-ratio-comparison.svg): comparison figure.
- [`verification/artifact-audit.json`](verification/artifact-audit.json): independent artifact audit.

## Reproduction

The implementation commit is `2835029`.  From the repository root, using the
recorded single-thread environment:

```bash
JULIA_CONDAPKG_BACKEND=Null \
JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" \
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
VECLIB_MAXIMUM_THREADS=1 \
julia --project=. --startup-file=no \
  examples/compare_toric3d_z_aspect_ratio.jl

.CondaPkg/.pixi/envs/default/bin/python \
  examples/analyze_toric3d_z_aspect_ratio.py \
  results/toric3d-h6-z-aspect-ratio-2026-10-15/production/toric3d-z-threshold-raw.csv \
  --output-dir results/toric3d-h6-z-aspect-ratio-2026-10-15/analysis
```

The historical H+6 and H+5 result directories were not changed, and
`document/note.tex` was not modified.
