# H+6 full-prefix feedback for X and Z

This study extends the previous complete two-slice study by correcting the
entire prepared prefix Q0 union ... union Qr after every transfer r=1,...,T,
including the last. All completed physical syndromes must be zero afterward;
future slices are untouched. The decoder sees only measured syndromes.

The timing, physical noise and scoring conventions match the two-slice study:
no separate Z seed correction before transfer 1; one final full-code outcome
per arm; final-only versus full-prefix replay share an immutable fault record
and final decoder. X prepares ideal logical |00> with 13L^3 noisy transfer
outputs. Z prepares logical |++> with 3L^2 noisy seed-H outputs plus 13L^3
transfer outputs. Each specified gate output independently receives an error
with probability p. The X and Z studies are separate noise channels.
Measurements, physical corrections and idles are ideal. State failure means
either of the two residual logical parities is nonzero.

X uses the existing full-plaquette BP+OSD decoder for every growing prefix,
with the same declared physical-p prior and decoder settings as before. Z
uses unit-weight data-edge MWPM for every growing prefix. Full-prefix feedback
can revisit previously prepared slices; all those slices remain available.
Unlike the previous two-slice X protocol, later corrections can touch Q0 and
its logical-Z scoring cycles, so first-window logical-parity freezing is not
assumed here.

The historical built-in full-prefix Z replay repairs the seed separately.
This study uses a separate replay entry point to retain the immediately
preceding study's no-r=0-repair convention and explicit final-transfer audit.
Historical source and result files are preserved. document/note.tex is not
modified.

The production design was frozen in `selection.json` after 420 broad pilot
records and an additional 1,800-record Z refinement. All 2,220 pilot records
are excluded from production estimates and fitting. The refinement showed
size suppression at p=0.025 and size-order reversal by p=0.040, motivating a
local Z grid around p=0.030--0.035.

| Channel | L values, with T=L | Gate-output probabilities p | Records per point |
| --- | --- | --- | --- |
| X-only | 3, 5, 7, 9, 11 | 0.035, 0.040, 0.045, 0.050, 0.055, 0.060, 0.065 | 1,000 |
| Z-only | 3, 5, 7, 9, 11 | 0.0250, 0.0275, 0.0300, 0.0325, 0.0350, 0.0375, 0.0400 | 1,000 |

This is 70 production points, 70,000 independent fault records and 140,000
protocol outcomes. Within each channel and point, the two arms use the same
record. X and Z use independent streams and are simulated separately, not
as simultaneous noise. The six production directories are disjoint shards
of the declared grid. Scalar UInt64 point seeds, raw counts, four-bin paired
histograms, individual logical marginals and immutable batch checkpoints
are retained.

The initial implementation spent much of the X runtime calculating binary
matrix parity through strided rows. The final implementation computes the
same parity by XOR of packed columns selected by the error frame. Exact
comparisons covered 198 prefix cases, two complete L=11 paired histories,
and all 2,220 pilot records. The focused parity and full-prefix checks passed
16,162 assertions; the runner passed all 105 checks after the change.
This changes calculation cost, not the decoder or physical protocol.

The initial run and its completed 95,339-assertion integration suite are
preserved under `pre-optimization/`, with a separate archive manifest.
Final pilots and production use only the optimized source. Their grids,
seeds, sample counts, fit windows, numerical model and bootstrap settings
are identical to the original frozen selection. Repeated streams are not
pooled or counted as additional independent samples. The transition audit
compares scientific counts at matching immutable checkpoint generations.
That audit passed for all 2,220 pilot records and all 55,250 production
records saved by the initial implementation. The archived production
records are verification evidence only and do not add to the sample count.

The final integrated Julia suite passed 96,525 assertions. Independent
replays regenerated two complete 1,000-record production points from their
point seeds: X at L=5, p=0.05 and Z at L=11, p=0.0325. Those checks covered
16,000 prefix corrections and 4,000 final outcomes, passed 230,000
assertions, and reproduced both arms' logical counts and paired histograms.
The independent replay builds the fault propagation and prefix mappings
without importing the study runner. Its reports and the source bindings
are retained in `verification/`.

The existing guarded local finite-size fitting policy, 200 conditional
bootstrap replicates, drop-L=3 sensitivity, inverse-size correction and
narrower-window checks are retained. All seven production rates in each
channel define its frozen fit window. Thresholds are reported only if the
existing acceptance checks pass; a visually plausible collapse alone is
insufficient. Error bars in the reusable no-grid plotting theme are 95%
Wilson intervals. Paired differences use approximate pointwise 95%
Bonferroni--Wilson intervals.

## Completed results

All 70 production points completed their 1,000 records. This gives 70,000
fault records, 140,000 paired final outcomes and 490,000 full-prefix
feedback checkpoints across the production runs. All required physical
syndromes were checked after each feedback correction.

For example, at L=11:

| Channel and physical p | Final-only state failures | Full-prefix state failures |
| --- | --- | --- |
| X-only, p=0.050 | 743/1,000 (74.3%) | 153/1,000 (15.3%) |
| Z-only, p=0.030 | 744/1,000 (74.4%) | 224/1,000 (22.4%) |

The numerical fits converged, but neither channel passed every frozen
threshold acceptance check. Their candidate parameters are therefore
diagnostic estimates, not established thresholds:

| Channel | All-size candidate p_c | Candidate nu | Goodness-of-fit p-value | Excluding L=3: candidate p_c |
| --- | --- | --- | --- | --- |
| X-only | 0.05377 | 1.452 | 0.959 | 0.05440 |
| Z-only | 0.03360 | 1.172 | 0.00493 | 0.03434 |

For X, local fit quality and every tested sensitivity pass. The strict
crossing check fails: the adjacent-size sign-change brackets have only the
single common endpoint p=0.055, which does not contain the fitted p_c.
Thus approximately 5.4% is a useful provisional value, while the stored
accepted-threshold field remains null.

For Z, the primary local goodness-of-fit check fails, as do the crossing,
cubic-fit quality and inverse-L correction checks. Excluding L=3 improves
the local goodness-of-fit p-value to 0.156. However, adding the inverse-L
correction shifts the candidate to 0.03641, beyond the fixed sensitivity
tolerance. Approximately 3.4% remains provisional. The drop-L=3 outputs
are predefined sensitivity fits; they do not repeat the entire acceptance
policy on the reduced size set or establish a separate threshold.

Both primary fits produced 200 usable conditional bootstrap replicates.
Those intervals describe sampling uncertainty conditional on each fitted
model and do not resolve failed crossing or model checks. The complete
parameters, intervals, residuals, checks and sensitivity results are in
`analysis/report.json`. The global Z fit exactly matches the independent
Z-only recomputation in `verification/z-independent-fit.json`.

All six PNG figures were visually inspected. Each also has PDF and SVG
exports. `full-prefix-comparison` shows both protocols, the paired-difference
figure shows their difference, and each channel has all-size and drop-L=3
figures with scaling-collapse insets. The axes have no grids.

## Reproduction and file locations

`verification/commands.txt` records the commands; `selection.json` binds
the grid, seeds, pilot inputs, protocol and numerical fitting sources.
The simulation and analysis source archives and corresponding SHA-256
manifests are in `verification/`. For a clean reproduction, restore those
archives into a compatible checkout before using the commands. The study
validates the exact source manifest, package versions and one-thread
runtime; merely copying result files into a different current checkout does
not satisfy those checks. Existing checkpoints require `--resume`, which
regenerates and discards completed random records before continuing.

The completed artifact bundle is exported to the normal project directory:
`/Users/rongyuqing/jcode/TopoNoise.jl/results/toric3d-h6-full-prefix-2026-10-15/`.
PNG, PDF and SVG figures are in `analysis/`, alongside raw marginal and
paired-difference tables and the machine-readable fit report. Historical
studies remain separate and unchanged.
The reusable plotting helper and no-grid style are also included in
`theme/`; the analysis source archive retains their original project paths.
