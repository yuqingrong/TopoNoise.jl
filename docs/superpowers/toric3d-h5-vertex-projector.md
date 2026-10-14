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

Final recovery scores two logical-X and two logical-Z bits, either-pair unions, and all 16 Pauli classes. The X-only logical-X union also denotes `|00_L>` preparation failure. Logical Z is a state-independent diagnostic and is never labeled logical-zero infidelity. Each observable has 16 joint protocol masks and six paired contrasts. Rate intervals are Wilson 95%; paired intervals use simultaneous bounds on the two discordance probabilities, giving nominal 95% coverage for each contrast. They are not simultaneous across the six contrasts, rates or sizes. Threshold intervals are conditional on the scaling model and frozen window; sensitivity is separate.

PyMatching 2.4.0 minimizes unit physical length, followed by lexicographic feasibility fixing of canonical endpoint/path signatures. ldpc 2.4.1 uses minimum-sum BP, 50 iterations, scaling one, parallel schedule, one thread, OSD_CS order two, and `p_eff=clamp(p,1e-12,1-1e-12)`. Neither decoder is circuit-level maximum likelihood. Ideal measurement and boundary assumptions limit hardware interpretation. The arms also change measured checks, residual constraints, and decoding policy, so their differences do not isolate correction range.

## Reproduction

The completed delivery is the user-confirmed **L=3,5,7 preliminary subset**. Its runner, inherited seed mapping and independent audits are documented in [the preliminary README](../../results/toric3d-h5-diagonal-feedback-2026-10-13/preliminary-l3-7/README.md). Resume that runner only for an ordinary interruption; completed output needs no further simulation. Larger sizes remain deferred.

The commands below describe the original six-size design for provenance. **They are not instructions to resume the deliberately stopped full scan.** Run any explicitly requested fresh reproduction from the repository root with installed dependencies:

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
Completed on 2026-10-14: cubic **L=3,5,7**, 60 X-only points and 135 Z-only points, 1000 paired records at every point. The 195,000 records produce 780,000 protocol replays. All points are complete, with **zero invalid recoveries**. Sampling, decoder construction and replay total 1532.72 seconds for X and 12168.14 seconds for Z (3.806 hours combined). These phase totals include the original cost of reused records and exclude startup, warmup, plotting and work on deferred sizes; they are not elapsed time for the final resumed segment.

**No H+5 threshold, scaling exponent or fitted interval is reported.** The unchanged minimum-four-size guard rejects all 14 frozen candidate windows across the eight channel/protocol arms before optimization. The requested 1000 bootstrap replicates were **not executed**; no sensitivity fits were performed. Empty bootstrap/sensitivity CSVs mean not performed, not successful fits. Candidate windows in the analysis report remain scout diagnostics, not estimates from these measurements.

The X-only observable is either-pair logical X, which also measures logical-zero preparation failure. The Z-only observable is state-independent either-pair logical Z, not logical-zero infidelity. The following illustrative points were selected after the run to show the rate dependence and the low end of the Z grid; they are not prespecified confirmatory comparisons. Entries are observed percentages, each from 1000 paired records.

| Channel | p | L | None | One step | Window complete | Full prefix |
|---|---:|---:|---:|---:|---:|---:|
| x_only | 0.039448164 | 3 | 13.6% | 14.8% | 17.7% | 18.2% |
| x_only | 0.039448164 | 5 | 6.5% | 12.3% | 17.1% | 17.1% |
| x_only | 0.039448164 | 7 | 4.2% | 8.6% | 15.0% | 15.0% |
| x_only | 0.077446809 | 3 | 57.5% | 58.6% | 49.5% | 50.4% |
| x_only | 0.077446809 | 5 | 71.6% | 72.7% | 58.2% | 58.2% |
| x_only | 0.077446809 | 7 | 73.9% | 73.5% | 62.5% | 62.2% |
| z_only | 0.0096470592 | 3 | 51.6% | 50.2% | 41.0% | 41.1% |
| z_only | 0.0096470592 | 5 | 74.5% | 75.5% | 64.7% | 64.7% |
| z_only | 0.0096470592 | 7 | 76.3% | 74.4% | 71.7% | 71.7% |

Feedback is not uniformly beneficial in these finite systems. For X at L=7, p=0.039448164, the paired changes relative to no feedback are +4.4 percentage points for one step (95% interval +1.38 to +7.38), and +10.8 points for either complete protocol (+6.99 to +14.50). At p=0.077446809 the window-complete change is -11.4 points (-17.07 to -5.62), and full-prefix is -11.7 points (-17.34 to -5.94). These opposite signs illustrate dependence on the error rate, not a threshold result.

For Z at the lowest sampled production rate p=0.0096470592, window-complete changes the diagnostic by -10.6 points at L=3 (-15.39 to -5.70) but by -4.6 points at L=7 (-10.10 to +0.94). The latter interval includes zero. Most Z curves are near the 0.75 union probability expected for two uniformly random logical bits, so this grid provides little unsaturated large-size evidence. This observation is descriptive; it does not establish zero or positive asymptotic threshold. All six paired contrasts at every sampled point are retained in the CSVs. Intervals are per contrast and do not correct for searching across this grid.

The final-source pilot used 20 records per channel at L=3, p=0.01: X counts [0,0,0,0] and Z [8,7,6,6], with zero invalid recoveries and recorded phase times 0.08279/0.06401 seconds. The scout completed 66 points, 33,000 paired records and 132,000 arm replays, with zero invalid recoveries and X/Z phase times 394.21/842.82 seconds. Both stages were independently audited and excluded from final counts.

The original scout selection chose L=3,5,7,9,11,13, 20 X rates and 45 Z rates without extensions. Its median-scout model projected 17.2 hours total but underestimated dense large-size costs; actual first-rate X L13 cost was about three times that projection. At the user's runtime-directed request, the six-size run was stopped at 106,350 records. A four-size precursor was subsequently stopped at 206,000 records when the user explicitly confirmed **L=3,5,7** (clarifying an earlier L=4 typo). The reduction was directed by runtime, not observed failure values. Both selections and checkpoints are preserved and are **deferred, not completed**. L9,11,13 require a new explicit resumption request.

The final 195-point checkpoint reused exactly 186 complete points (186,000 records) from that four-size precursor and completed 9000 additional records. Every retained count array, record/signature digest, timing and seed was compared independently. Seeds retain their indices in the parent 260-element `Random.Xoshiro(1234)` UInt64 stream; subsetting never reseeds existing data. The precursor's 80 reused X points also match the original six-size checkpoint exactly. Rates, protocols, decoder settings and final scientific source hashes remain unchanged.

| Provenance object | SHA-256 |
|---|---|
| Original six-size selection | `5041dd833aa8f28f318b4f5fb495f5025d61dfeb87f1a8ab8d76e5edcf6d7e00` |
| Final three-size selection | `1597b9941132d3f24c11c68ec1933ea639f5289225fe80ff5070367e81416e72` |
| Final raw CSV | `ba3765ff260f09e109510bd574e19aa9a71e2da560c0e6916db1a8f44fbd8390` |
| Final metadata | `11da3b08a5b6fdb1caec613554084077558cdc3f6d6ed386c4cab1c334609107` |
| Final checkpoint | `3180191df43597ed1f9768d7ee4b8cc1c251326045f806761f27939b8d5333d7` |

Artifacts live under `results/toric3d-h5-diagonal-feedback-2026-10-13/preliminary-l3-7/`: `production/` contains raw and aggregate CSVs, metadata, the checkpoint, both independent audit JSONs, and X/Z plots in PDF/SVG/PNG; `analysis/` contains the guarded report, all primary paired contrasts, exact illustrative measurements and measured-curve figures; `verification/` contains the full Julia/Python logs and run log. The independent audit reconstructs all 23 metrics, all 16 joint masks and six pair tables, Wilson intervals, paired intervals and standard errors, logical labels from Pauli classes, digests, timings, source/data hashes and warmup exclusion.

The fresh final regressions passed **45,610 Julia assertions in 145 test sets** and **51 Python tests**. The analysis driver passed after its presentation-only update; no scientific source or guard changed. Final figures and the rebuilt note were visually inspected. Existing unrelated LaTeX duplicate-label and overfull-box warnings remain outside the new H+5 content.

<!-- H5_MEASURED_RESULTS_END -->

Historical `Toric3DSliceEncoder` and its saved H+6 data retain their names, schemas, and numerical results. The seven-gate circuit is labeled **plaquette-transfer H+6**. H+5/H+6 threshold comparisons are not controlled circuit-only comparisons: seed preparation, the upper instrument, and historical channel-dependent seed noise also differ.
