# Direct-projector H+5 diagonal feedback experiment

This study implements the approved diagonal H+5 design with the existing independent post-gate-output noise model. The physical circuit is H(n1), CX(n1,n2), CX(n1,n3), CX(n1,e1), CX(n1,e2), CX(n1,e3). Each cell has eleven noisy outputs. Boundaries, measurements, corrections, initialization, and idles are ideal. The experiment is motivated by diagonal feedback ideas; it is not a reproduction of a supplied 3D gate decomposition.

The public dimensions are `(T,Lu,Lv)`, with `n=3(T+1)LuLv`, `k=2`, `dZ=2min(Lu,Lv)` and `dX=(T+1)min(Lu,Lv)`. Physical period-two parallel edges retain distinct identities. All four arms replay one immutable fault record and one sampled clean even upper branch per shot. Feedback runs only at layers `1:T-1`.

| Protocol | X feedback | Z feedback |
|---|---|---|
| none | final recovery | final recovery |
| one_step | periodic Q-slice MWPM, first/last path edges | retiring/frontier strip matching, one future dummy iff odd; retiring endpoint edges only |
| window_complete | BP+OSD on two data slices and completed plaquette window | ordinary open-boundary MWPM on the retiring/frontier strip |
| full_prefix | BP+OSD on all prepared edges | open-boundary MWPM on all prepared edges |

One-step X contracts each visible defect pair once; old-target X faults remain local terminal noise. One-step Z clears retiring detectors and allows defects only at the frontier. Complete feedback clears all prepared checks. Every complete recovery satisfies `Hc=s`; partial routes satisfy the independently remeasured residual equation. Both Pauli sectors use a common pre-correction frame in combined-channel replay.

The upper controller sees only observed `a=b xor u` and applies `C(a)`. The simulator compares it with clean `C(b)`, giving relative frame `z xor C(a) xor C(b)`. The correction tree is physical-edge-ID BFS rooted at `A_(T-1,0,0)` after deleting both logical-X supports. An arbitrary open tree path need not have zero logical overlap: the naive 2x2x2 upper-(1,1) path `[36,20,13]` overlaps logical X1 once; the declared cut-tree fixture is `[32,15,21]`. This correction convention changes noisy logical-Z diagnostics relative to the naive tree.

Final recovery scores two logical-X and two logical-Z bits, either-pair unions, and all 16 Pauli classes. The X-only logical-X union also denotes `|00_L>` preparation failure. Logical Z is a state-independent diagnostic and is never labeled logical-zero infidelity. Each observable has 16 joint protocol masks and six paired contrasts. Rate intervals are Wilson 95%; paired intervals use discordance counts with simultaneous bounds. Threshold intervals are conditional on the scaling model and frozen window; sensitivity is separate.

PyMatching 2.4.0 minimizes unit physical length, followed by lexicographic feasibility fixing of canonical endpoint/path signatures. ldpc 2.4.1 uses minimum-sum BP, 50 iterations, scaling one, parallel schedule, one thread, OSD_CS order two, and `p_eff=clamp(p,1e-12,1-1e-12)`. Neither decoder is circuit-level maximum likelihood. Ideal measurement and boundary assumptions limit hardware interpretation. The arms also change measured checks, residual constraints, and decoding policy, so their differences do not isolate correction range.

## Reproduction

Run from the repository root with installed dependencies:

```sh
export JULIA_NUM_THREADS=1 JULIA_CONDAPKG_BACKEND=Null
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python"
export MPLCONFIGDIR=/tmp/toponoise-mpl
julia --project=. examples/compare_toric3d_vertex_projector_feedback.jl --stage pilot --dimensions 3x3x3 --rates 0.01 --channels x_only,z_only --shots 20 --seed 1234 --output results/toric3d-h5-diagonal-feedback-2026-10-13/pilot
julia --project=. examples/compare_toric3d_vertex_projector_feedback.jl --stage scout --dimensions 3x3x3,5x5x5,7x7x7 --rates 0,0.0001,0.0003,0.001,0.003,0.01,0.02,0.03,0.05,0.075,0.10 --channels x_only,z_only --shots 500 --seed 1234 --output results/toric3d-h5-diagonal-feedback-2026-10-13/scout
.CondaPkg/.pixi/envs/default/bin/python examples/analyze_toric3d_vertex_projector.py select --input-csv results/toric3d-h5-diagonal-feedback-2026-10-13/scout/toric3d-h5-vertex-projector-feedback-raw.csv --input-metadata results/toric3d-h5-diagonal-feedback-2026-10-13/scout/toric3d-h5-vertex-projector-feedback-metadata.toml --output results/toric3d-h5-diagonal-feedback-2026-10-13/selection.json
julia --project=. examples/compare_toric3d_vertex_projector_feedback.jl --stage production --selection results/toric3d-h5-diagonal-feedback-2026-10-13/selection.json --shots 1000 --seed 1234 --output results/toric3d-h5-diagonal-feedback-2026-10-13/production
.CondaPkg/.pixi/envs/default/bin/python examples/analyze_toric3d_vertex_projector.py fit --selection results/toric3d-h5-diagonal-feedback-2026-10-13/selection.json --csv results/toric3d-h5-diagonal-feedback-2026-10-13/production/toric3d-h5-vertex-projector-feedback-raw.csv --metadata results/toric3d-h5-diagonal-feedback-2026-10-13/production/toric3d-h5-vertex-projector-feedback-metadata.toml --bootstrap 1000 --seed 1234 --output results/toric3d-h5-diagonal-feedback-2026-10-13/analysis
```

If selection requests extensions, run `--stage extend-scout --selection <selection.json> --output <new-shard-directory>`. Reselect with `--previous-selection <selection.json>` and repeated `--input-csv`/`--input-metadata` for the initial scout and every extension in order. Parent selection SHA-256 and full cumulative ancestry are validated and bound to checkpoint identity. The 0.49 bound limits scout extensions; the prescribed exact ±20% production neighborhoods may exceed it within the physical noise domain [0,1]. No-crossing fallbacks use strictly interior measured rates. The process terminates after at most three extension rounds; exhausted arms remain threshold-ineligible diagnostics.

The frozen selection independently chooses X and Z grids, every candidate window, and four to six cubic sizes from scout timing. It retains every ranked crossing cluster and unions all five-rate candidate neighborhoods. Production outcomes cannot choose or change windows. A threshold requires a single passing candidate, a consistent cluster covering every adjacent-size pair, an interior monotone likelihood fit, goodness-of-fit p>=0.01, and at least 80% accepted bootstrap replicates. Fixed-size improvement alone is not evidence of a positive threshold.

Append `--resume` to the identical scan command after interruption. Configuration, source hashes, pinned backend versions, point seeds, and complete random records (including upper branches) are checked. Checkpoints save every batch; completed-point CSV and metadata update incrementally. Invalid recovery creates an atomic `failure.toml` and aborts. Such a checkpoint cannot be resumed as successful data until the cause is resolved. Timing and zero-noise warmup shots never enter scientific counts.

## Verification and measured results

Implementation validation includes exhaustive 88-location X/Z injections, independent matrix-conjugation Pauli evolution through every correction, asymmetric multi-checkpoint histories, all period-two matching syndromes, frozen legacy behavior, and exact interrupted replay.

The final-source Julia regression passed 45,610 assertions across 145 test sets, including the existing 2D and H+6 suites. The Python analysis regression passed all 51 tests. The whole-change review findings were fixed and covered by failing-then-passing tests: checkpoint selection/stage binding, usable high-probability candidate windows and strictly interior fallbacks, report input-hash validation, and preservation of independently observed failure syndromes and physical paths.

<!-- H5_MEASURED_RESULTS_BEGIN -->
The final-source timed pilot completed 20 paired records per channel at (3,3,3), p=0.01, seed 1234, with zero invalid recoveries. X failure counts were [0,0,0,0]; Z counts were [8,7,6,6] in protocol order. The recorded build/sampling/replay times were 0.08279 seconds (X) and 0.06401 seconds (Z), excluding startup, warmup, and plotting. These pilot counts do not establish a threshold. Scout and production results remain pending.
<!-- H5_MEASURED_RESULTS_END -->

Historical `Toric3DSliceEncoder` and its saved H+6 data retain their names, schemas, and numerical results. The seven-gate circuit is labeled **plaquette-transfer H+6**. H+5/H+6 threshold comparisons are not controlled circuit-only comparisons: seed preparation, the upper instrument, and historical channel-dependent seed noise also differ.
