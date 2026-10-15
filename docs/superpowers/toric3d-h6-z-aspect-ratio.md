# H+6 Z aspect-ratio diagnostic

The final-only H+6 Z experiment now has a dedicated aspect-ratio runner and
analysis:

```bash
julia --project=. examples/compare_toric3d_z_aspect_ratio.jl
python examples/analyze_toric3d_z_aspect_ratio.py \
  results/toric3d-h6-z-aspect-ratio-2026-10-15/production/toric3d-z-threshold-raw.csv \
  --output-dir results/toric3d-h6-z-aspect-ratio-2026-10-15/analysis
```

The runner uses the existing exact circuit-correlated Z fault model and final
unit-data-weight MWPM decoder.  It does not invoke any intermediate feedback
protocol.  The default grid compares `(T,Ly,Lz)=(L,L,L)` and `(2L,L,L)` for
`L=3,5,7`, probabilities `.0045,.005,.0055,.006,.0065`, 1,000 samples per
point, and seed 1234.

All 30,000 fault records completed.  At every fixed `(L,p)`, the measured
failure rate was higher for `T=2L`.  At `p=.005`, the `T=L -> T=2L` changes
were `1.6% -> 5.1%` for `L=3`, `1.6% -> 4.5%` for `L=5`, and
`1.0% -> 4.2%` for `L=7`.  At `p=.0065`, they were `3.0% -> 7.8%`,
`3.1% -> 10.9%`, and `3.9% -> 15.2%`, respectively.

This is consistent with the exact code geometry

\[
n=3(T+1)L^2,\qquad d_X=(T+1)L,\qquad d_Z=2L.
\]

Doubling `T` adds noisy transfer locations without increasing `dZ`.  The
result therefore does not support using `T=2L` as a remedy for logical-Z
failure in the final-only circuit.

The study is a fixed-size diagnostic.  It does not fit a threshold or scaling
exponent: only three small sizes were requested, the samples for the two
depths are unpaired, and decoder/model systematics are not included in the
Wilson intervals.  Full counts, intervals, provenance, figure, and audit are
in
[`results/toric3d-h6-z-aspect-ratio-2026-10-15`](../../results/toric3d-h6-z-aspect-ratio-2026-10-15/README.md).
