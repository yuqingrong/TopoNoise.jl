# Sparse row microkernel gate result

This experiment is isolated from the production decoder and from the
committed row-strip correctness reference.  It applies one Fourier sector,
one fixed periodic-seam branch, and one row with a custom dense-array MPS
kernel.  Local MPO tensors are compiled to CSR-like transition tables and
large accumulation/output buffers are reused.

## Correctness

The untruncated `L=2` test compares the complete coefficient vector with the
ITensor zip-up action and passes at `2e-12` absolute and `2e-11` relative
tolerance.

For the finite-bond representative row (`L=3`, `p=0.02`, Fourier sector 2,
strip 9, seam 0, `chi=240`, cutoff `1e-18`), the custom output differed from
the committed ITensor reference by:

- relative coefficient-vector norm: `6.387671174607e-9`;
- absolute log-norm difference: `9.214851104389e-15`.

The small vector difference is expected because the committed implementation
right-canonicalizes the input MPO/MPS, permits a first zip-up bond of 480,
and then performs a second right-to-left truncation to 240.  This prototype
uses a single raw-automaton left-to-right truncation to 240.

## Warm row measurement

The row has 22 sites and 499 nonzero local transitions.  With one BLAS thread,
after one warm call:

| quantity | result | gate |
|---|---:|---:|
| wall time | `0.393346 s` | at most `0.25 s` |
| allocated bytes | `37,191,760` | at most `150,000,000` |
| persistent workspace | `15,004,416` bytes | informational |
| maximum output bond | `240` | `240` |

One immediate profiled repeat took `0.393721 s` and allocated the same
`37,191,760` bytes.  Its time split was:

| phase | seconds | fraction |
|---|---:|---:|
| sparse transition accumulation | `0.025989` | 6.6% |
| SVD | `0.364967` | 92.7% |
| output/carry compression | `0.002730` | 0.7% |
| final open-boundary contraction | about `0.000001` | below 0.1% |

The committed ITensor branch application took `0.390276 s` once in the same
process.  That differs from the custom repeat by less than 1%, and there was
no multi-sample timing study, so the apparent ordering is measurement noise;
these data do not establish that either implementation is faster.  Preparation
of the true pre-row boundary and coefficient-vector comparison are excluded
from all row timings.

If all SVD time is removed, the measured accumulation, output/carry work, and
final contraction sum to about `0.028720 s`.  This is a favorable kernel-only
lower bound, not an achievable decoder time: the row update still needs a
rank-revealing compression with controlled error.  It does show that sparse
transition bookkeeping is no longer the limiting phase.

## Decision

The allocation gate passes, but the time gate fails.  The custom full-SVD
kernel and the ITensor reference are statistically indistinguishable in these
single-run timings, and both are around `0.39 s` rather than the required
`0.25 s`.  The experiment stops here.  Bookkeeping changes cannot deliver the
required gain because SVD accounts for more than 90% of the measured time.

A randomized, iterative, or Gram-matrix truncated factorization might reduce
the SVD phase, but this experiment did not record the local singular spectra,
benchmark such a method, or establish an error bound for it.  The signed
Fourier sectors can cancel and the requested cutoff is `1e-18`, so any such
replacement needs its own exact-small and finite-chi accuracy audit.  The
present measurements therefore support rejecting this full-SVD kernel only;
they do not support a general claim that truncated methods cannot work.  The
next performance experiment should either test a controlled truncated solver
or change the contraction order/factor the row into bond-2/4 layers.

