# Rotated-planar construction/channel scaling figures

## Goal

Add a reproducible analysis workflow for comparing the two existing rotated-planar encoder constructions (`:as`, `:bp`) under independent X-only and Z-only circuit-noise sweeps. It must produce raw logical-failure curves, exploratory finite-size threshold estimates, scaling-collapse insets, one combined 2×2 figure, and four standalone panels.

The supplied reference image is a visual target only. Its title, numerical threshold values, and plotted data are not reused.

## Fixed scientific convention

- Every circuit prepares the existing default logical state `|0_L>` (`logical_state=:zero`).
- The construction remains an independent circuit choice: `construction=:as` or `construction=:bp`.
- X-only noise means `p_x=p`, `p_z=0`; its reported observable is the existing logical-X failure rate and standard error.
- Z-only noise means `p_x=0`, `p_z=p`; its reported observable is the existing logical-Z failure rate and standard error. This is a code-level phase-parity failure, even though it does not flip a Z-basis readout of `|0_L>`.
- The initial production sweep uses distances `[3, 5, 7]`, error rates `0:0.01:0.16`, 10,000 shots per point, the `:x_ns` boundary orientation, and the `:gate_layer` clock. Distances, rates, shots, batch size, seed, orientation, and clock remain command-line configurable.
- The four panels are independently fit. With only three distances, all threshold and exponent estimates must be labeled **exploratory**, not publication-quality threshold claims.

## Data model and execution

Keep `LogicalFailureScan` and `scan_logical_failure` unchanged: they continue to mean the existing balanced convenience sweep (`p_x=p_z=p`).

Add a dedicated channel-sweep result that records:

- construction, error channel (`:x_only` or `:z_only`), and its reported logical observable;
- the shared lattice/noise/shot metadata;
- per-distance/per-rate `LogicalFailurePoint` values with their actual `p_x` and `p_z` fields;
- the physical-rate vector, a reproducibility seed, and the per-series derived seed.

Add a comparison result that owns exactly four channel sweeps in this deterministic order: As/X-only, As/Z-only, Bp/X-only, Bp/Z-only. It stores the common configuration plus fit results.

The public comparison runner accepts a caller-owned RNG and a master seed. It derives a deterministic independent RNG stream per construction/channel series so adding or reordering a panel does not perturb another panel's raw data. Every point continues to use the production syndrome-only PyMatching estimator; no threshold information is passed to decoding.

## Fitting procedure

For each construction/channel panel:

1. Form pairwise crossings between adjacent-distance raw curves using linear interpolation of the sampled rates.
2. Require at least two usable crossing estimates in the scanned physical-rate range. Their uncertainty-weighted aggregate is the raw threshold estimate `p_c`.
3. Fit `nu` by minimizing a deterministic scaling-collapse scatter objective in the coordinates
   `(p - p_c) * d^(1 / nu)`.
4. Run a seeded parametric binomial bootstrap over the observed counts/shots. Each replicate repeats crossing aggregation and collapse fitting. Report central estimates plus bootstrap uncertainty/intervals.
5. If crossings cannot be formed or the numerical fit is ill-conditioned, return a structured `fit unavailable` status with a diagnostic. Render raw curves, but omit the threshold line and inset rather than manufacturing a result.

The fitting module is analysis-only. It does not change the encoder, noise model, syndrome extraction, matching graph, or logical-failure estimator.

## Figures and files

Produce five data-identical visualizations:

- one 2×2 comparison figure: rows As/Bp and columns X-only/Z-only;
- four standalone panels: As/X-only, As/Z-only, Bp/X-only, Bp/Z-only.

Each panel has the physical error rate on its x-axis, the channel-appropriate logical failure rate on its y-axis, one colored error-bar curve per distance, construction/channel-aware title, and a panel-specific dashed `p_c` line when fitting succeeded. An inset shows the scaling-collapse x coordinate and reproduces the same curves. The default axes use the supplied rate range and automatic y limits; command-line rate-grid settings intentionally control the physical x scale.

For an output basename `comparison`, save:

- `comparison-raw.csv`: every raw Monte Carlo point, including construction, channel, actual `p_x`/`p_z`, selected logical observable, counts, rates, standard errors, shots, seed, orientation, clock, and state;
- `comparison-fits.csv`: pair crossings, `p_c`, `nu`, bootstrap summaries, fit status, and diagnostics per panel;
- `comparison-combined.{svg,pdf,png}`;
- `comparison-as-x-only.{svg,pdf,png}`, `comparison-as-z-only.{svg,pdf,png}`, `comparison-bp-x-only.{svg,pdf,png}`, and `comparison-bp-z-only.{svg,pdf,png}`.

All figures are generated from the saved in-memory raw data, so the combined and standalone plots cannot disagree. Output file names are leaf names only and stay inside the chosen directory.

## Command line

Add a dedicated example CLI rather than overloading the balanced-noise CLI. Its documented defaults are the initial production sweep above. It accepts comma-separated distances or a `p-min`/`p-max`/`p-step` grid, shot/batch counts, master seed, boundary orientation, clock, output directory, basename, bootstrap replicate count, and explicit `--no-fit` for raw-only previews.

On success it prints every saved artifact and a compact per-panel fit summary. It validates unsupported channels, empty/nonmonotonic grids, invalid rates, bad bootstrap counts, and unsafe output names before simulation begins. It records `fit unavailable` as a normal analysis result, not a process failure.

## Testing and verification

Develop test-first.

- Validate channel rate construction exactly: X-only never injects Z and Z-only never injects X.
- Verify the logical-X/logical-Z observable mapping, fixed `|0_L>` metadata, deterministic four-series order, independent per-series RNG reproducibility, and preservation of existing balanced-scan behavior.
- Exercise CSV schemas and exact row counts for all four series.
- Use small synthetic, known-crossing data to test `p_c` interpolation, collapse fitting, bootstrap determinism, and fit-unavailable diagnostics without relying on noisy Monte Carlo.
- Smoke-test one combined and four standalone SVG/PDF/PNG artifacts in a temporary directory, plus the new CLI at `d=3`, `p=0`, and two shots.
- Run all legacy and new tests with the established PyMatching environment, then run `git diff --check`.

## Non-goals

- No claim of a precision threshold from the initial three-distance data.
- No balanced `p_x=p_z` comparison series in this workflow.
- No repeated syndrome rounds, measurement error, circuit-aware/correlated decoding, periodic geometry, holes, multi-logical patches, or changes to existing A_s/Z and B_p/X conventions.
