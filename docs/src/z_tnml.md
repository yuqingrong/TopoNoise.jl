```@meta
CurrentModule = TopoNoise
```

# Z-only tensor-network decoding

The Z-only tensor-network decoder estimates the posterior probability of each
of the four logical classes of `Toric3DEncodedSlice`. It chooses the class with
the largest approximate posterior. Its output is a two-bit logical-class
decision, not a data-qubit correction vector.

## Physical model

The model is the diagonal H+6CNOT transfer circuit with a honeycomb toric-code
seed prepared in logical `|++>`. Each seed-H output and each active transfer
H/CNOT output has an independent Z fault with probability `p`, where
`0 <= p < 0.5`. Faults propagate through the remaining ideal circuit. Seed
compatibility measurements and the final stabilizer syndrome are perfect;
there is no intermediate feedback. Initialization and idle periods are ideal.

The probabilities are per circuit fault location. Propagation can correlate
the final data errors. `build_toric3d_z_fault_model` retains these correlations
exactly in its propagated-fault matrix.

## Decode a syndrome

Only the package import is needed:

```julia
using TopoNoise

code = Toric3DEncodedSlice(3, 3, 3)  # transfer steps, Ly, Lz
p = 0.01
model = build_toric3d_z_fault_model(code)
topology = build_toric3d_z_tn_topology(model)
settings = Toric3DZTNSettings(
    max_bond_dimension=240,
    cutoff=1e-18,
    negative_tolerance=1e-10,
)
decoder = build_toric3d_z_tn_decoder(topology, p; settings)

# A deterministic nonzero syndrome from one raw circuit fault.
syndrome = BitVector(topology.raw_detector_matrix[:, 1] .!= 0)
result = decode_toric3d_z_tn(decoder, syndrome)

@show result.valid result.probabilities result.logical_mask result.log_syndrome
```

The single raw column only supplies the observed syndrome. The decoder sums
over all compatible fault histories, rather than being told which fault
occurred. `syndrome` uses the row order of `topology.raw_detector_matrix`, which
is the order of `star_check_matrix(code)`. It must have one binary entry per
detector and even parity. Wrong lengths, nonbinary entries, and odd parity
are rejected. A zero-probability conditioned syndrome, for example this
nonzero syndrome at `p=0`, cannot be normalized and raises a `DomainError`.

The builder also accepts the code or its existing propagated-fault model:

```julia
decoder_from_code = build_toric3d_z_tn_decoder(code, p; settings)
decoder_from_model = build_toric3d_z_tn_decoder(model, p; settings)
```

A decoder can be reused for multiple syndromes at the same `p`. The topology
is independent of `p` and can be reused when building decoders at other rates.

The complete executable example defaults to `L=3`, `p=0.01` and prints the
four class weights and the chosen logical bits:

```sh
julia --project=. examples/toric3d_z_tnml.jl
```

Optional positional arguments set `L` and `p`. The script performs one
deterministic decode; it does not run a Monte Carlo experiment or a threshold
scan. First use includes Julia compilation, and increasing `L` or the bond
limit can substantially increase time and memory.

## Result and logical-bit order

`decode_toric3d_z_tn` returns a named tuple:

| Field | Meaning |
|:--|:--|
| `valid` | Whether the final signed weights satisfy the negative tolerance. This is a numerical admissibility check, not an accuracy or convergence guarantee. |
| `probabilities` | Four normalized signed estimates of the class posterior. Small negative values are preserved, not clipped. |
| `logical_mask` | `argmax(probabilities) - 1` when valid, or `-1` when invalid. |
| `log_syndrome` | Approximate natural logarithm of the syndrome probability, accumulated during contraction. |
| `raw_probabilities` | The same signed normalized weights, retained for numerical inspection. |
| `diagnostics` | Bond dimensions, local SVD truncation statistics, and contraction metadata. |

The mask is `ell1 + 2ell2`, with this Julia vector order:

| Vector index | Mask | `(ell1, ell2)` |
|:--|:--|:--|
| 1 | 0 | `(0, 0)` |
| 2 | 1 | `(1, 0)` |
| 3 | 2 | `(0, 1)` |
| 4 | 3 | `(1, 1)` |

The bits are final Z-error parities against the two X membranes returned by
`toric3d_encoded_logical_x(code)`, in that order. Use `result.valid` before
interpreting `logical_mask`. Selecting the class does not construct a
physical recovery operator.

## Exact topology and approximate contraction

`Toric3DZTNTopology` stores the detector signatures
`D = H * F (mod 2)` and logical signatures `B = X * F (mod 2)` of all raw
circuit faults. It merges raw faults only when their two detector endpoints
and both logical bits agree. Each effective edge retains its multiplicity
and original raw-column indices in `original_columns`.

`toric3d_z_tn_edge_probabilities(topology, p)` gives the odd-parity probability
of each effective edge. For multiplicity `k`, this is

```math
q_k = \frac{1-(1-2p)^k}{2}.
```

This reduction is exact. The contraction builds an exact star
matrix-product operator (MPO), including the observed syndrome bit, before
compressing the canonical boundary matrix-product state (MPS). Each star MPO
is applied only through its last nonidentity site. Compression uses that same
finite support. Although a nonunitary update can change Schmidt values across
later cuts, it cannot increase their Schmidt ranks; those later bonds were
already bounded by `max_bond_dimension`. The implementation retains their
updated spectra until a later star reaches each cut, instead of reapplying the
cutoff after every earlier star. This changes the finite-cutoff approximation
policy as well as reducing QR/SVD work, so it needs independent convergence
checks. Normalization and bond diagnostics still contract through the active
boundary. Detector sites beyond it are an untouched, normalized,
dimension-one product suffix and contribute mass one. The boundary compression
introduces the approximation. The implementation preserves signed weights and
uses a fixed plane-major detector order, so there is no per-shot
contraction-path search.

For the two frozen `L=3`, `p=0.02`, bond-240 cases, support-only compression
reduced the number of local SVD events from 714 to 605. In the validation run,
mean decode time decreased from about 51.1 seconds to 43.0 seconds. This timing
is an implementation benchmark rather than a scaling guarantee.

The topology fixes the number of local SVD events before any shot is decoded:

| `L` | Support-only events | Earlier active-frontier events | Reduction |
|:--|--:|--:|--:|
| 3 | 605 | 714 | 15.3% |
| 5 | 7,079 | 8,105 | 12.7% |
| 7 | 36,609 | 40,824 | 10.3% |
| 9 | 125,963 | 137,871 | 8.6% |
| 11 | 339,173 | 366,278 | 7.4% |

Event counts alone do not predict wall time because the SVD matrix dimensions
also vary with the syndrome and retained bond dimensions.

`Toric3DZTNSettings()` uses `max_bond_dimension=240`, `cutoff=1e-18`, and
`negative_tolerance=1e-10`. The bond limit must be a positive integer; the
cutoff and negative tolerance must be finite and in `[0,1)`. Increasing the
bond limit or decreasing the cutoff changes the compression accuracy and
cost. Assess convergence for the sizes, rates, and syndromes of interest.

`diagnostics.max_bond_dimension` records the largest bond dimension used.
`maximum_local_truncation_error`, `sum_local_truncation_error`, and
`local_svd_events` summarize local relative SVD truncation errors on the
current star support. They are not directly comparable with diagnostics from
the earlier active-frontier compression schedule. The legacy name
`discarded_weight` is an alias for
`sum_local_truncation_error`; it is neither lost probability mass nor a bound
on posterior error.

## Validation scope

The numerical settings were selected for the earlier active-frontier
compression by comparing bond limits 240 and 256 on a fixed 48-case `L=3`
diagnostic set, followed by 300 held-out samples at `p=0.006, 0.01, 0.02`.
Those results do not validate the current support-only cutoff schedule. The
support-only implementation has instead been checked on 78 exact small-graph
syndromes and two frozen `L=3`, `p=0.02`, bond-240 cases. On the frozen cases,
the largest posterior L1 error against the independent dense-frontier
reference was approximately `0.000407`.

These are empirical errors on those cases, not bounds for arbitrary syndromes.
Before using support-only compression for a production scan, repeat the saved
diagnostic and held-out comparisons. The present checks do not establish
accuracy at larger sizes, under other noise models, or for feedback protocols.
A successful `valid` check and a small local SVD diagnostic alone do not extend
that validation. This integration does not supply a new threshold estimate.
