# Feedback strength in diagonal 3D toric-code preparation

This experiment compares three feedback policies in the encoded diagonal sweep:

1. **No intermediate feedback (`none`)**: run the complete transfer circuit and
   decode once at the end.
2. **Slice-local feedback (`slice_local`)**: use a constant-width correction
   region in the sweep direction. For X faults this is the established
   one-step target feed-forward rule; for Z faults it is a syndrome-complete
   three-check-sheet measurement window with recovery on the newest two data
   slices.
3. **Full-prefix feedback (`full_prefix`)**: after every nonfinal transfer,
   measure and decode every relevant syndrome check in the prepared prefix and
   correct anywhere in that prefix.

The comparison keeps the same diagonal transfer circuit and post-gate-output
noise convention in all three protocols. For each Monte Carlo shot, one
immutable fault record is sampled and replayed through all protocols. The
result therefore gives paired differences between the complete implemented
protocol policies. For X, the intermediate check set, routing, correction
support, and decoder differ between slice-local and full-prefix feedback. For
Z, the slice-local and full-prefix policies additionally perform the same
seed-time X-star measurement and Z recovery, while `none` does not. The
comparison must therefore not be interpreted as isolating only the spatial
support of an otherwise identical correction.

## Result status

<!-- FINAL_RESULTS_START -->
The runtime-scoped production study is complete: 36,000 X-channel and 56,000
Z-channel paired fault records, each replayed through all three protocols.
Every complete recovery was syndrome-valid. The independent artifact audit
reconstructed 219 saved joint tables, 657 marginal records, 657 paired records,
and all 192 guarded-fit input rows with zero invalid recoveries. Recorded
decoder construction, sampling, and three-protocol replay total 10.734 compute
hours (5.773 X and 4.962 Z); this is not elapsed wall time.

**No primary fit passes every registered threshold guard.** The table therefore
reports diagnostic finite-size candidates, not threshold estimates:

| Channel | Policy | Registered result | Diagnostic evidence |
|---|---|---|---|
| X | `none` | Inconclusive | Missing the `L=5,7` crossing; the quadratic candidate pins to the lower window edge and is nonmonotone. |
| X | `slice_local` | Inconclusive | The quadratic candidate pins to `p=.025`, is nonmonotone, and has deviance goodness `p=3.35e-44`. |
| X | `full_prefix` | **Withheld** | Adjacent crossings are 5.159%, 5.586%, and 5.684%. The candidate is `p_c=5.297%`, `nu=1.391`, but deviance `54.77/31` gives `p=0.00530`, below the 1% guard. |
| Z | `none` | Inconclusive in this window | The candidate `p_c=2.936%` has `nu=0.3` on its lower bound; the selected high-probability window cannot resolve the known low-probability onset. |
| Z | `slice_local` | **Withheld** | Candidate `p_c=3.934%`, `nu=1.423`; the `L=13,17` crossing is missing, the fit is nonmonotone, and deviance `62.80/23` gives `p=1.49e-5`. |
| Z | `full_prefix` | **Withheld** | Candidate `p_c=3.923%`, `nu=1.417`; the `L=17,21` crossing is missing, the fit is nonmonotone, and deviance `63.68/23` gives `p=1.10e-5`. |

Because every primary was withheld before bootstrap, conditional threshold
intervals were not generated. X full-prefix has the clearest transition-like
finite-size crossing: its individually accepted cubic and endpoint-window sensitivities
give 5.293%--5.323%, while the fixed-`L^-1` correction gives 6.154%. These are
model diagnostics because the primary fit failed. All Z-feedback sensitivities
remain withheld; they place the transition-like onset near 3.9% but do not
establish an asymptotic threshold. Missing crossings are not evidence for a
zero threshold.

The paired outcomes do establish protocol differences at fixed size and
probability. At X `L=11,p=.05`, state-failure rates are 0.736 (`none`), 0.765
(`slice_local`), and 0.177 (`full_prefix`). Full-prefix minus none is -0.559
with conservative 95% interval `[-0.606,-0.506]`; full-prefix minus slice-local
is -0.588 `[-0.632,-0.538]`, while slice-local minus none is compatible with
zero. At Z `L=21,p=.035`, the rates are 0.7545, 0.2445, and 0.2415. Full-prefix
minus none is -0.513 `[-0.549,-0.475]`, while full-prefix minus slice-local is
-0.003 `[-0.0146,0.0086]`. Thus active Z feedback is effective, but this study
does not resolve a threshold difference between its local and prefix forms.

These conclusions do not replace the separate targeted estimates already in
the note. The earlier low-probability Z final-only study estimated 0.517%, and
the earlier one-step X slice-local study estimated 1.482%; the present windows
do not refit either number. See the
[primary curves](../../results/toric3d-diagonal-feedback-2026-10-12/analysis/threshold-fit.png),
[sensitivity chart](../../results/toric3d-diagonal-feedback-2026-10-12/analysis/threshold-sensitivity.png),
[fit record](../../results/toric3d-diagonal-feedback-2026-10-12/analysis/threshold-fit.json),
[paired statistics](../../results/toric3d-diagonal-feedback-2026-10-12/analysis/paired-differences.csv),
and [independent audit](../../results/toric3d-diagonal-feedback-2026-10-12/audit/production-artifact-audit.json).
<!-- FINAL_RESULTS_END -->

The threshold guards require at least four sizes, at least five probabilities,
a sign-changing crossing for every adjacent-size pair in the chosen window,
an interior and monotone likelihood fit, deviance goodness-of-fit
`p >= 0.01`, and at least 80% usable bootstrap fits. A numerical candidate is
retained for diagnosis when any guard fails, but it is not reported as a
threshold.

## Geometry, states, and noise

For `T` transfers and periodic transverse dimensions `Ly,Lz`, the encoded
diagonal slab has

\[
n=3(T+1)L_yL_z,\qquad k=2,\qquad
d_X=(T+1)\min(L_y,L_z),\qquad d_Z=2\min(L_y,L_z).
\]

Each transfer cell applies one H and six CNOTs. Faults are independently
sampled after all thirteen participating one- or two-qubit gate outputs. The
X experiment applies only X faults; the Z experiment applies only Z faults.
Initialization, idle qubits, measurements, classical processing, and Pauli
corrections are ideal.

The two channels deliberately prepare different logical states so that the
reported residual Pauli is state failure:

- The X experiment starts from encoded logical `|00>` and has no noisy seed
  gates. Failure is either residual logical-X membrane.
- The Z experiment applies H to all `3*Ly*Lz` seed edges and puts an
  independent Z fault after each H output. Ideal Z-type compatibility
  measurements followed by X correction complete the intended encoded
  logical `|++>` preparation; these checks commute with the sampled Z faults
  and do not remove them. Failure is either residual logical-Z winding.

Consequently, X and Z threshold numbers are not two sectors of one otherwise
identical state-preparation experiment. Protocols can be compared directly
within a channel; cross-channel numerical comparisons also include the seed,
logical-state, and decoder differences.

## What weak and strong feedback do

The slice-local X protocol measures the old active slice's Z compatibility
checks before the next transfer. Unit-weight MWPM constructs complete matching
paths, but the one-step rule retains only the first and, when distinct, last
edge of each path. After the noisy transfer, it applies the propagation image
of that partial route to fresh targets. This is target feed-forward rather
than a direct recovery on the old data. Its correction has bounded quantum
support in the sweep direction, although the matching computation over the
periodic slice is global. The partial route need not clear the measured slice;
the final BP+OSD decoder handles the remaining physical errors.

A literal one-slice Z decoder is not syndrome-complete for the retained noise
model because faults occur on both old controls and new targets. The local Z
protocol first measures both colors of seed X stars at `r=0` and applies a
complete MWPM Z recovery on the seed edges. After each nonfinal transfer it
measures three adjacent X-star check sheets, applies a complete MWPM Z recovery
on the newest two data slices, and verifies that the complete prepared-prefix
star syndrome is zero. The window width stays constant as `T` grows.

Full-prefix feedback repeatedly decodes the entire prepared volume after each
nonfinal transfer. In the X experiment it measures Z plaquettes and applies a
complete BP+OSD X recovery. This clears the complete prefix syndrome, whereas
the one-step slice-local route can leave residual active-slice syndrome before
the affected qubits control another transfer. In the Z experiment,
full-prefix feedback shares the local protocol's seed X-star/Z-recovery step,
then measures all completed prefix X stars and applies a complete MWPM Z
recovery. The three-sheet local Z protocol already restores zero complete-
prefix star syndrome at every checkpoint, so extending its support mainly
changes available recovery paths and logical-class/tie choices; it need not
outperform local Z feedback. Full-prefix feedback costs total check
participation `O(T^2*Ly*Lz)`, compared with `O(T*Ly*Lz)` for constant-width
feedback.

This active Z feedback is different from the Z-compatibility feedback in the
older encoded Z-only scan. There, “compatibility” means measuring Z-type checks
and applying X corrections, so pure Z faults do not trigger it. Here, the
seed, moving-window, and prefix feedback measure X-type stars and apply Z
corrections, so they do respond to pure Z faults.

Every complete intermediate and final decoder recovery is checked against its
input syndrome. The deliberately partial X one-step feed-forward is not itself
required to clear that syndrome, but its underlying complete MWPM route is
validated before truncation. A complete-recovery mismatch atomically saves the
protocol, checkpoint, shot/batch, dimensions, probability, and point seed,
then aborts the point. No invalid recovery is silently included in a logical-
error rate.

## Decoders and information boundary

X prefix and final recovery use BP+OSD with minimum-sum updates, 50 BP
iterations, scaling factor 1, parallel schedule, one thread, and `OSD_CS`
order 2. The decoder uses a uniform declared physical-`p` prior, clamped to
`[1e-12,1-1e-12]`.

Z feedback's seed-star, moving-window, prefix, and final recovery use
deterministic unit-data-weight MWPM. Every decoder receives only the measured
syndrome, geometry, declared `p` where applicable, and fixed settings. Sampled
fault bits are never decoder input.

These baseline models neglect nonuniform gate exposure, correlations created
by CNOT propagation, and correlations conditioned on earlier feedback. The
thresholds are therefore decoder-dependent estimates under ideal measurement
and feedback, rather than circuit-level maximum-likelihood thresholds.

## Registered numerical study

The 200-shot paired scouts were used only to choose production windows. They
showed no visible finite-size crossing for X `none`, X `slice_local`, or Z
`none`, while bracketing an X `full_prefix` transition near 0.052 and Z
`slice_local` and `full_prefix` transitions near 0.039. These adaptive,
low-count scout observations are not threshold estimates and are excluded from
all production fits. The bootstrap-usability guard applies to the registered
production analysis; it is not the reason the scouts are excluded.

The following target grid was frozen before production counts were generated:

| Channel | Cubic sizes | Physical probabilities | Paired shots per point |
|---|---|---|---:|
| X | 5, 7, 9, 11, 13 | .025, .030, .040, .045, .050, .055, .060, .070, .080 | 1,000 |
| Z | 9, 13, 17, 21, 25 | .025, .030, .0325, .035, .0375, .040, .045, .050, .060, .070 | 2,000 |

Measured runtime set the final size scope to four complete sizes per channel.
The X input set is `L=5,7,9,11` at all nine frozen X rates: 36 size-rate points
and 36,000 paired fault records. The Z input set is `L=9,13,17,21` at the first
seven frozen rates through `p=.045`: 28 size-rate points and 56,000 paired
fault records. This keeps a guard point above the transition and leaves six
probabilities in either endpoint-trimmed sensitivity. The higher frozen Z
rates chiefly probe saturation and were stopped on the measured runtime before
the four-size threshold fit was run. Original checkpoints remain resumable;
the completed rectangular analysis inputs and their source hashes are recorded
separately.

The `L=13` X and `L=25` Z shards are resumable partial provenance only and are
excluded from every fit and aggregate. A deliberate interruption exposed a
now-fixed handling bug: `InterruptException` was caught by recovery-failure
bookkeeping, and scalar X replay could have recorded a few triplets before the
batch completion counter advanced. The original artifacts are archived in
`x-L13-interrupt-diagnostic`, `z-L25-interrupt-diagnostic`, and
`z-L21-interrupt-diagnostic`. The live checkpoints were repaired from atomic
counts or deterministic replay and report zero invalid recoveries. The partial
X `L=13`, Z `L=25`, and original-shard Z `L=21,p=.045` records remain excluded;
the archived interruption records are not evidence of decoder-invalid
recoveries.

The master seed is 1234. Each `(channel,size,p)` point owns a deterministic
substream, and each size is a disjoint resumable shard. Timing and warmup
records are excluded from all scientific counts. The frozen selection and
scout hashes are in
[`scan-selection.toml`](../../results/toric3d-diagonal-feedback-2026-10-12/scan-selection.toml).

Each raw row stores the eight joint Boolean outcome counts for the three
protocols. Protocol marginals, 95% Wilson intervals, discordant pairs, and the
paired difference (second protocol minus first) are derived from this joint
table. The analysis independently reconstructs those quantities before any
fit.

## Reproduction

Use the repository's managed Python decoder environment:

```sh
export JULIA_CONDAPKG_BACKEND=Null
export JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python"
```

The runnable simulation is
[`compare_toric3d_diagonal_feedback.jl`](../../examples/compare_toric3d_diagonal_feedback.jl).
For example, one X production shard is:

```sh
julia --project=. examples/compare_toric3d_diagonal_feedback.jl \
  --sizes 5 \
  --rates .025,.03,.04,.045,.05,.055,.06,.07,.08 \
  --channels x_only --shots 1000 --batch-size 32 --seed 1239 \
  --x-routing one_step \
  --output results/toric3d-diagonal-feedback-2026-10-12/production/x-L5
```

The corresponding Z shard changes the channel, grid, shot count, and seed:

```sh
julia --project=. examples/compare_toric3d_diagonal_feedback.jl \
  --sizes 9 \
  --rates .025,.03,.0325,.035,.0375,.04,.045,.05,.06,.07 \
  --channels z_only --shots 2000 --batch-size 32 --seed 11243 \
  --x-routing one_step \
  --output results/toric3d-diagonal-feedback-2026-10-12/production/z-L9
```

Add `--resume` to continue an exactly matching checkpoint. Configuration
mismatches are rejected before sampling.

For the runtime-limited `L=21` Z shard, the first six probabilities remain in
the original ten-rate checkpoint and the seventh frozen point was run in a
disjoint checkpoint with its original point seed. The checked assembly script
combines those completed states without drawing another fault record:

```sh
julia --project=. \
  results/toric3d-diagonal-feedback-2026-10-12/production/assemble_z_L21_p045.jl

.CondaPkg/.pixi/envs/default/bin/python -B \
  results/toric3d-diagonal-feedback-2026-10-12/analysis-input/build_inputs.py
```

The source checkpoints, point seeds, and hashes remain alongside the combined
artifact. The original ten-rate checkpoint is still resumable.

After collecting every raw shard, run the guarded analysis with
[`fit_toric3d_diagonal_feedback.py`](../../examples/fit_toric3d_diagonal_feedback.py):

The command names the eight frozen analysis inputs explicitly so partial
provenance directories cannot enter through a wildcard. The X files are byte
copies of the complete production rows. The Z files retain the seven complete
probabilities through `p=.045`; their manifest records every source and hash.

```sh
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/toponoise-mpl \
.CondaPkg/.pixi/envs/default/bin/python \
  examples/fit_toric3d_diagonal_feedback.py \
  results/toric3d-diagonal-feedback-2026-10-12/analysis-input/x-L5.csv \
  results/toric3d-diagonal-feedback-2026-10-12/analysis-input/x-L7.csv \
  results/toric3d-diagonal-feedback-2026-10-12/analysis-input/x-L9.csv \
  results/toric3d-diagonal-feedback-2026-10-12/analysis-input/x-L11.csv \
  results/toric3d-diagonal-feedback-2026-10-12/analysis-input/z-L9.csv \
  results/toric3d-diagonal-feedback-2026-10-12/analysis-input/z-L13.csv \
  results/toric3d-diagonal-feedback-2026-10-12/analysis-input/z-L17.csv \
  results/toric3d-diagonal-feedback-2026-10-12/analysis-input/z-L21.csv \
  --x-p-min .025 --x-p-max .08 --x-min-size 5 \
  --z-p-min .025 --z-p-max .045 --z-min-size 9 \
  --degree 2 --nu-min .3 --nu-max 5 --goodness-alpha .01 \
  --bootstrap 1000 --seed 1234 \
  --output-dir results/toric3d-diagonal-feedback-2026-10-12/analysis
```

Finally, independently reconstruct every saved marginal and paired statistic
from the joint tables and compare the raw rows with their checkpoint states:

```sh
.CondaPkg/.pixi/envs/default/bin/python -B \
  results/toric3d-diagonal-feedback-2026-10-12/audit/audit_production.py
```

The command above is the primary analysis for the runtime-scoped complete
grids. The larger frozen target remains recorded in `scan-selection.toml`, and
the scope decision was made before the four-size fit. Required sensitivity
fits vary the polynomial degree, minimum size, probability endpoints, and a
fixed leading `L^-1` correction.

## Interpretation limits

The study uses separate pure-Pauli channels, ideal measurements and
corrections, ideal idles, a particular encoded boundary, and baseline
decoders. It does not estimate a depolarizing threshold, include measurement
or correction faults, prove a local classical implementation, or numerically
reproduce a threshold stated in *Noise-resilient sequential circuits for
generating quantum order*. It tests the paper's diagonal-sweep feedback
mechanism in the explicit circuit and retained noise model used by this
repository.
