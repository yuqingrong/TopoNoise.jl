# H+6 no feedback results

Legacy H+6-CNOT plaquette-transfer circuit with T=Ly=Lz=L, no intermediate feedback, and one final decode. X and Z are separate state-preparation experiments with different seeds and decoders.

Production contains 35,000 X trials and 450,000 Z trials. The Z primary fit passes the declared guards: p_c=0.00534660, with a conditional 95% bootstrap interval [0.00530812, 0.00537976]. Fit variants span 0.00525086–0.00534660; this model sensitivity range is not a confidence interval. X does not pass the threshold guards, so no X threshold or exponent is established.

Both channels have saved fitted models. The X diagnostic candidate is p_0=0.00298108 (0.298108%), with nu_0=0.85970. Its goodness-of-fit p-value is 4.50e-7, below the declared 0.01 cutoff; it also lacks a common adjacent-size crossing bracket and fails the four sensitivity guards. Its plotted collapse is therefore labeled diagnostic. This does not establish the absence of a physical X threshold. The accepted Z fit has nu=0.91509.

## Folder guide

- This README summarizes the protocol and results.
- [data](data/) contains final analysis tables, production counts, and metadata.
- [figures](figures/) contains the final PDF plots.
- [archive](archive/) contains one compressed study archive and its integrity manifest.

## Figures

- [h6-threshold-failure-rates (PDF)](figures/h6-threshold-failure-rates.pdf)

The PDF uses the same boxed-axis, no-grid paper theme as the one-step figures. Both panels show the original data with 95% Wilson error bars and their saved scaling fits. The Z dashed line and shaded band show the accepted estimate and its conditional bootstrap interval; the X candidate is explicitly marked as rejected. This is a presentation-only redraw, with no resampling or refitting. Its [renderer](archive/redraw_no_feedback.py) and [theme snapshot](archive/paper-no-grid.mplstyle) are preserved for reproduction. The original PDF remains in the study archive.

## Data

- [bootstrap-samples.csv](data/bootstrap-samples.csv)
- [sensitivity-fits.csv](data/sensitivity-fits.csv)
- [threshold-points.csv](data/threshold-points.csv)
- [threshold-report.json](data/threshold-report.json)
- [production-x/metadata.toml](data/production-x/metadata.toml)
- [production-x/raw.csv](data/production-x/raw.csv)
- [production-z/metadata.toml](data/production-z/metadata.toml)
- [production-z/raw.csv](data/production-z/raw.csv)
- [selection.json](data/selection.json)
- [fit-guard-policy.json](data/fit-guard-policy.json)
- [runtime-environment.json](data/runtime-environment.json)

The data files retain their original bytes. Paths and hashes inside JSON/TOML files describe
the original run layout and remain provenance records. The presentation redraw manifest,
where present, also records the original image formats.

## Archive and reproduction

[study-records.tar.gz](archive/study-records.tar.gz) preserves 2,264 original files:
the full scientific README, raw data, metadata, checkpoints, pilot runs, scripts, validation
records, and all original PDF plots. The 2 non-PDF figure files were omitted.
[manifest.json](archive/manifest.json) records every archived and omitted file's original
SHA-256 hash, the published file mapping, and checksums for this cleaned export.

Archive members retain their original `results/toric3d-h6-final-only-2026-10-15/` paths. Extract into a separate
reproduction checkout to restore that layout before using the original scripts or checkpoints.
The original README and artifact manifests inside the archive are historical records;
their PNG/SVG links and image checks no longer apply to this PDF-only export. Use the new
manifest for export integrity. Recorded simulation source and environments are still required
to rerun or resume the study.

The untouched canonical source checkout remains at
`/Users/rongyuqing/.codex/worktrees/61d4/TopoNoise.jl`.
