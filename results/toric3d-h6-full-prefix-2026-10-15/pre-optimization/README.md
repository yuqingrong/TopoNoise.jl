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

The existing guarded local finite-size fitting policy, 200 conditional
bootstrap replicates, drop-L=3 sensitivity, inverse-size correction and
narrower-window checks are retained. All seven production rates in each
channel define its frozen fit window. Thresholds are reported only if the
existing acceptance checks pass; a visually plausible collapse alone is
insufficient. Error bars in the reusable no-grid plotting theme are 95%
Wilson intervals. Paired differences use approximate pointwise 95%
Bonferroni--Wilson intervals.

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

Artifacts will be copied to the normal project directory:
`/Users/rongyuqing/jcode/TopoNoise.jl/results/toric3d-h6-full-prefix-2026-10-15/`.
PNG, PDF and SVG figures will be in `analysis/`, alongside raw marginal and
paired-difference tables and the machine-readable fit report. Historical
studies remain separate and unchanged.
