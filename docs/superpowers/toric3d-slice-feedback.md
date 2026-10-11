# Diagonal-slice transfer with measurement feedback

This experiment changes the preparation circuit and active-slice feedback while retaining the previous simulation's gate noise: each participating qubit independently receives an X fault or, in a separate run, a Z fault with probability `p` after every H/CNOT. Initialization, measurement, Pauli feed-forward, and idle qubits are ideal. It preserves the existing rough-slab experiment and its results.

The construction is a derived 3D realization motivated by [Noise-resilient sequential circuits for generating quantum order](https://arxiv.org/html/2609.13383v1). It is not an author-supplied circuit or a reproduction of a numerical threshold from that paper. The paper motivates transfer from an old control slice to a new target slice and defect pairing within a slice. Our explicit geometry, boundaries, gate decomposition, routing rule, and retained gate-output noise are specified below.

## Geometry and fixed boundaries

Use diagonal coordinates `t=x+y+z`, `u=x`, `v=y`. The transverse coordinates `u,v` are periodic with lengths `Ly,Lz`; the transfer depth is `T=steps`. Each slice `t=0:T` contains three incoming edges per vertex. An edge `E[t,d,u,v]` ends at `(t,u,v)`, with directions and transverse displacements

\[
s_1=(1,0),\qquad s_2=(0,1),\qquad s_3=(0,0).
\]

All physical edge IDs remain distinct, including when a periodic length is two. The number of qubits is `3 Ly Lz (T+1)`. The initial slice starts in the all-zero state, and its individual Z operators are fixed boundary checks. Six-edge X stars occur at bulk vertices; final frontier stars have three incoming edges. Physical Z plaquettes contain two new edges and two old edges.

These boundaries specify a unique state, `k=0`. The full ranks are

\[
\operatorname{rank}H_X=T L_yL_z,\qquad
\operatorname{rank}H_Z=(2T+3)L_yL_z,
\]

where `H_Z` includes the initial boundary checks. The physical plaquette matrix alone has rank `(2T+1)LyLz-1`. The old rough-slab logical operators and distances do not apply to this construction.

## Transfer gate and compatibility measurement

Write `b_d(r)` for the incoming computational-basis bit of direction `d` at transverse vertex `r` in the old slice. Define

\[
K_{ij}(b,r)=b_j(r-s_i)\oplus b_i(r-s_j),\qquad
Q(b,r)=K_{12}(b,r)\oplus K_{13}(b,r)\oplus K_{23}(b,r).
\]

The six-edge Z compatibility operator `Q` is an actual stabilizer of each ideal prepared prefix. Its syndrome can therefore be measured before the next transfer. The local matrix has `LyLz` rows and `3LyLz` columns, column weight two, and rank `LyLz-1`; the syndrome has even parity.

For each new vertex, in increasing `t,u,v` order, initialize the three new edges to zero and execute:

1. H on `new1`.
2. CNOT from `new1` to `new2`.
3. CNOT from `new1` to `new3`.
4. The two old controls defining `K12` into `new2`.
5. The two old controls defining `K13` into `new3`.

The two controls within each last pair use the displayed order in the definition of `Kij`. There are seven gates and thirteen participating gate outputs per vertex. In the absence of faults the output bits have the form

\[
b'_1(r)=a(r),\quad b'_2(r)=a(r)\oplus K_{12}(b,r),\quad
b'_3(r)=a(r)\oplus K_{13}(b,r),
\]

where H creates the coherent bit `a(r)`. The third plaquette condition follows from `Q=0`. Binary stabilizer conjugation verifies the physical checks, initial boundary checks, and compatibility checks at every checkpoint.

Old qubits are controls throughout this transfer. A Z error already on an old slice does not propagate into future slices. A Z error arising within the fresh triad can spread locally through its remaining gates, but exhaustive single-fault tests check that its final support has at most three qubits. This removes the previous circuit's longitudinal Z-propagation mechanism; it does not by itself prove a threshold.

The depth-independent support bound also follows directly from the gate order. CNOT sends `Z_control` to itself and `Z_target` to `Z_control Z_target`. A fault on `new1` occurs after its only H and remains on that control. A fault on `new2` or `new3` can propagate back to at most the two old controls that have yet to act on that target. Once a qubit becomes an old-slice control, further transfers cannot spread its Z component. Thus a single Z gate-output fault occupies at most three edges in two adjacent slices, at any transfer depth. The compatibility measurements detect X faults, so under pure Z noise the assisted and unassisted results coincide: the change of circuit addresses Z propagation, while slice feedback addresses X propagation.

## Matching and routed target feed-forward

Each physical old-slice edge joins the two compatibility checks it anticommutes with. The resulting periodic graph is triangular. MWPM uses unit edge weights and only the measured compatibility syndrome. Its full data-edge recovery must exactly reproduce that syndrome, or the simulation stops.

The default `routing=:one_step` extracts matched defect pairs and a deterministic shortest path for each pair. Each endpoint advances one graph edge toward its partner. A length-one path annihilates the pair. Overlapping corrections combine modulo two. The stronger comparison `routing=:full_path` uses the entire matching recovery in that transfer and is labeled separately. The residual syndrome of a partial route is reported explicitly; a partial route is not represented as a complete recovery.

For a routed virtual correction `C` on old-slice bits, apply only its propagated target image

\[
F(C)=(0,K_{12}(C),K_{13}(C))
\]

on the new triads, after their fixed transfer gates. This is ideal Pauli feed-forward determined from the pre-transfer syndrome, rather than a correction on the old data or a change to the gate fault locations. With no new faults during that transfer,

\[
Q_{\rm new}(r)=\bigl[Q_{\rm old}\oplus Q(C)\bigr](r-s_1).
\]

New gate faults can add a fresh outgoing syndrome, so neither routing mode is assumed to clear the actual outgoing slice. Both protocols replay the same pre-sampled faults through the same H/CNOT schedule. The decoder cannot inspect those faults.

## Diagnostics and uncertainty

No global final recovery is performed. For each individual shot the simulation records:

| Diagnostic | Definition |
|---|---|
| Bulk star density | Violated six-edge stars divided by `(T-1)LyLz`; excludes the final frontier stars; defined as zero for `T=1` |
| Plaquette density | Violated physical plaquettes divided by `3TLyLz` |
| Final compatibility density | Violated final-slice `Q` checks divided by `LyLz` |
| Plaquette-cluster extent | Maximum longitudinal extent `max(t)-min(t)+1` of any connected violated-plaquette component, divided by `T`; zero when empty |

Two violated faces belong to adjacent cluster sites when they share a cube. At `(t,r)`, the cube has upper faces `B12,B13,B23` at `(t,r)` and lower faces `B12(t-1,r-s3)`, `B13(t-1,r-s2)`, `B23(t-1,r-s1)`. Cubes exist for `t=2:T`. This defines the cluster statistic without relying on a drawing or on decoder matches.

Means and sample standard errors are computed over independent shots. The approximate 95% intervals use the normal multiplier `1.959963984540054`, clipped to the diagnostic's allowed range. Paired differences use the sample variance of each shot's feedback-minus-unassisted difference. Checks within one shot are correlated and are never counted as independent Bernoulli trials. Consequently these are not Wilson intervals. Constant samples have zero sample SEM, which should not be interpreted as an upper confidence bound on unseen rare events; a one-shot timing run has undefined SEM.

These physical diagnostics test propagation and feedback behavior. They do not establish a logical-error threshold, preservation of an arbitrary encoded state, or a spatially local recovery of sublinear depth. In particular, global classical matching with parallel target Paulis is not a constant-depth local channel.

## Running and extending the scan

```julia
using TopoNoise, Random

comparison = compare_toric3d_slice_feedback(MersenneTwister(1234);
    sizes=[3,5,7], error_rates=[0,.001,.003,.01,.03,.05], shots=200,
    seed=1234, routing=:one_step,
    output_dir="results/toric3d-slice-feedback-2026-10-11",
    progress_io=stdout)
save_toric3d_slice_comparison(comparison,
    "results/toric3d-slice-feedback-2026-10-11")
```

Each `sizes` entry may instead be a `(steps,Ly,Lz)` tuple. The supplied RNG controls the samples; `seed` is recorded metadata. The standalone example runs an independently seeded warmup and timing pilot before the scan:

```sh
julia --project=. examples/compare_toric3d_slice_feedback.jl
```

CSV and TOML files are atomically updated after every completed channel point. The final save adds PDF, SVG, and PNG figures. Each CSV row contains both protocol means, SEMs and intervals plus their paired difference. Metadata records gate order, active outputs, noise assumptions, routing, package version, dimensions, seed, completion count, cluster convention, and statistical limitations.

When the local CondaPkg resolver is unavailable but the existing managed Python environment is installed, the documented project workaround also applies:

```sh
JULIA_CONDAPKG_BACKEND=Null \
JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" \
julia --project=. examples/compare_toric3d_slice_feedback.jl
```

## Initial paired scan

The completed exploratory scan uses `L=3,5,7`, `p=0,.001,.003,.01,.03,.05`, 200 paired shots per channel point, seed 1234, and one-step routing. The independent five-shot-per-channel timing pilot at `L=7,p=.03` took approximately 0.01 seconds after warmup, excluding figure compilation. All 36 channel points were saved incrementally.

For `L=7`, the X-noise results are:

| p | Final compatibility density, unassisted → feedback | Mean maximum cluster extent / depth, unassisted → feedback |
|---|---:|---:|
| 0.001 | 0.11510 → 0.01796 | 0.79571 → 0.32786 |
| 0.003 | 0.26265 → 0.04694 | 0.95929 → 0.44786 |
| 0.010 | 0.46490 → 0.14714 | 0.99857 → 0.78786 |

All six paired-difference intervals exclude zero. At `p=.001`, the cluster statistic corresponds to a mean maximum span of 5.57 layers without feedback and 2.30 layers with feedback. The Z-noise results are identical across protocols, as required by the compatibility measurement and the absence of X components. Every zero-noise statistic is zero. These are finite-size propagation results, not evidence sufficient to locate a threshold.

The saved [CSV](../../results/toric3d-slice-feedback-2026-10-11/toric3d-slice-feedback-raw.csv), [metadata](../../results/toric3d-slice-feedback-2026-10-11/toric3d-slice-feedback-metadata.toml), and [PDF figure](../../results/toric3d-slice-feedback-2026-10-11/toric3d-slice-feedback.pdf) contain the full grid, uncertainty estimates, and assumptions. SVG and PNG versions are saved beside the PDF.

