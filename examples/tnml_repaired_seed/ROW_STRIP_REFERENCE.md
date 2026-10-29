# Experimental row-strip TN reference

`row_strip_mps.jl` is a correctness and performance experiment. It is not the
production TN decoder and is not used by `run.jl` or the HPC wrappers.

The most useful variant contracts four Walsh/Fourier sectors of the two-bit
logical class. In sector `k`, an occupied effective edge with logical mask
`g` receives the signed weight

```text
p(edge) * (-1)^popcount(k & g).
```

This removes the four-state logical carrier. Fixing each periodic-column seam
separately leaves only the incoming two-bit cross-column mask and one source
parity bit in the row automaton, so its maximum bond is 8. Each sector MPS is
a signed partial partition function. The implementation normalizes it after
each row by a positive Euclidean MPS norm and tracks the logarithmic scale.
The four final signed sector values are transformed back with the inverse
Hadamard transform.

## Correctness evidence

- The dense per-row reference agrees with direct edge enumeration at `L=2`.
- The logical-carrying row MPO agrees with the dense per-row map.
- Fixed-seam MPO branches sum to the same dense per-row map.
- The Fourier decoder agrees with the independent dense-frontier exact ML
  decoder at `L=2`, including syndromes whose nonzero Fourier sectors have
  negative signs.
- The four-thread implementation is deterministic in the `L=2` test and
  agrees with both the serial implementation and exact ML.
- For the `L=3`, `p=0.02`, zero-syndrome timing case at `chi=240` and
  `cutoff=1e-18`, the Fourier result differed from exact ML by
  `6.03e-12` in posterior L-infinity norm and `1.69e-11` in log syndrome
  probability.

Run the concurrency check explicitly with more than one Julia thread:

```bash
JULIA_CONDAPKG_BACKEND=Null \
JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" \
OPENBLAS_NUM_THREADS=1 \
julia --project=. --startup-file=no --threads=4 \
  examples/tnml_repaired_seed/test_row_strip_fourier_threads.jl
```

The long `L=3`, cutoff-zero, `maxdim>=1024` proof was not run after the
generic implementation failed its performance gate. The independent exact
checks above establish the algebra and small-system implementation, but they
do not turn the finite-`chi` method into an exact ML decoder for larger `L`.

## Measured performance

All numbers below use one BLAS thread on the local Apple Silicon machine.
Cold timings include Julia compilation and should not be compared with a
steady-state shot.

These were exploratory rejection measurements made while the implementation
was still uncommitted; the raw timing harness and a machine-specific report
were not retained. They support the decision not to use this generic backend,
but they are not a reproducible performance claim. A future performance
comparison must use a committed benchmark harness and save Julia/package,
hardware, source-hash, and raw per-row timing metadata.

| L=3 variant | Julia threads | elapsed | allocations | peak RSS | max MPS bond |
|---|---:|---:|---:|---:|---:|
| logical + dynamic seam, cold | 1 | 215.50 s | 59.85 GB | about 2.71 GiB | 240 |
| logical + four fixed seams, cold | 1 | 173.93 s | 65.79 GB | 2.54 GB | 240 |
| Fourier + fixed seams, cold | 1 | 131.66 s | 33.81 GB | 2.51 GB | 170 |
| Fourier + fixed seams, warm | 1 | 39.38 s | 27.34 GB | 2.41 GB | 170 |
| Fourier + fixed seams, warm | 4 | 30.75 s | 27.33 GB | 2.53 GB | 170 |

Four-sector Julia parallelism is safe in the tested case, with BLAS held at
one thread, but yielded only a 1.28x warm speedup. The run is dominated by
memory traffic, tensor construction, SVD work, and garbage collection rather
than by work that scales cleanly across four cores. The slowest `L=3`
sector-row update took 2.42 seconds and allocated 1.23 GB by itself.

## Decision

The generic ITensor backend is rejected for production and no `L=5` timing
was attempted. Although the warm threaded `L=3` case is under one minute, its
per-row allocation and weak parallel scaling make a sub-minute `L=5` result
implausible without replacing the local contraction kernel. Keep this code as
an independent semantic reference for the next implementation.

## Next implementation gate

The next prototype should keep the Fourier and fixed-seam mathematics but
replace generic ITensor row application with a sparse, matrix-free zip-up
kernel:

1. Compile each local row tensor into compact transition tables. For every
   automaton/input state, store only its at-most-eight nonzero transitions.
2. Store each MPS tensor as a dense `Array{Float64,3}` with layout
   `(left_bond, physical_bit, right_bond)`.
3. Apply a local transition by indexed accumulation into preallocated work
   matrices. Never materialize a dense rank-4 MPO tensor or a general-purpose
   ITensor contraction.
4. Reuse thread-local workspaces for QR/SVD buffers, seam branches, branch
   addition, and the next-row MPS. Use `mul!`, `permutedims!`, and in-place
   LAPACK operations where possible.
5. Process the four Fourier sectors on four Julia threads with BLAS fixed to
   one thread. Each thread owns all buffers for one sector.
6. Merge the one-to-four fixed-seam branch MPSs once per row, then normalize
   the signed result by its positive Euclidean norm and update its log scale.
7. Cache geometry, transition tables, syndrome-bit variants, site insertion
   maps, and retirement maps for each strip. Only observed syndrome bits and
   boundary tensor values should change between shots.

The first gate is the slow `L=3` row: reduce it from 2.42 seconds and 1.23 GB
allocated to at most 0.25 seconds and 150 MB allocated while matching this
reference. If that 8-10x row speedup is achieved, the expected increase in
row/site work from `L=3` to `L=5` leaves a credible path to a four-thread
`L=5` shot below 60 seconds. A gain below about 5x is unlikely to overcome
the larger `L=5` boundary and SVD count, so the project should then switch to
a different contraction order or a factorized sequence of bond-2/4 layers.
