# H+6 one step X feedback results

Paired final-only and one-step X feedback on the legacy H+6-CNOT circuit, with T=L and L=3,5,7,9,11. Production contains 17,500 shared fault records, with 500 pairs per point.

Neither protocol passes every threshold guard. At p=0.015 and L=11, final-only fails in 243/500 trials (48.6%), versus 111/500 (22.2%) for one-step feedback. This is a fixed-size improvement; no threshold or exponent is established.

The visible figures are the later paper-style presentation redraws. Their local inset is exploratory (L=5,7,9,11 and p=0.010–0.020), and does not replace the original guarded analysis. The original analysis PDFs and pilot/shard plots are retained in the archive.

## Folder guide

- This README summarizes the protocol and results.
- [data](data/) contains final analysis tables, production counts, and metadata.
- [figures](figures/) contains the final PDF plots.
- [archive](archive/) contains one compressed study archive and its integrity manifest.

## Figures

- **[Combined X and Z one-step fit panels (PDF)](figures/one-step-feedback-x-z.pdf)** — the main side-by-side figure, using the existing themed panels with their original fit-status labels. X shows the exploratory restricted-window fit; Z shows the rejected primary fit. The two panels describe separate channel experiments.
- [x-one-step-failure-rates (PDF)](figures/x-one-step-failure-rates.pdf)
- [x-one-step-local-fit (PDF)](figures/x-one-step-local-fit.pdf)
- [x-one-step-paired-differences (PDF)](figures/x-one-step-paired-differences.pdf)
- [x-one-step-scaling-collapse (PDF)](figures/x-one-step-scaling-collapse.pdf)

The combined figure preserves the original vector plots without resampling or refitting. Its [assembly script](archive/combine_one_step_figures.py) reads the two individual PDFs.

## Data

- [bootstrap-samples.csv](data/bootstrap-samples.csv)
- [paired-differences.csv](data/paired-differences.csv)
- [points.csv](data/points.csv)
- [report.json](data/report.json)
- [sensitivity-fits.csv](data/sensitivity-fits.csv)
- [production-L11/toric3d-x-feedback-metadata.toml](data/production-L11/toric3d-x-feedback-metadata.toml)
- [production-L11/toric3d-x-feedback-raw.csv](data/production-L11/toric3d-x-feedback-raw.csv)
- [production-small/toric3d-x-feedback-metadata.toml](data/production-small/toric3d-x-feedback-metadata.toml)
- [production-small/toric3d-x-feedback-raw.csv](data/production-small/toric3d-x-feedback-raw.csv)
- [selection.json](data/selection.json)
- [provenance.json](data/provenance.json)
- [presentation/near-critical-fits.json](data/presentation/near-critical-fits.json)
- [presentation/redraw-manifest.json](data/presentation/redraw-manifest.json)

The data files retain their original bytes. Paths and hashes inside JSON/TOML files describe
the original run layout and remain provenance records. The presentation redraw manifest,
where present, also records the original image formats.

## Archive and reproduction

[study-records.tar.gz](archive/study-records.tar.gz) preserves 91 original files:
the full scientific README, raw data, metadata, checkpoints, pilot runs, scripts, validation
records, and all original PDF plots. The 18 non-PDF figure files were omitted.
[manifest.json](archive/manifest.json) records every archived and omitted file's original
SHA-256 hash, the published file mapping, and checksums for this cleaned export.

Archive members retain their original `results/toric3d-h6-one-step-2026-10-15/` paths. Extract into a separate
reproduction checkout to restore that layout before using the original scripts or checkpoints.
The original README and artifact manifests inside the archive are historical records;
their PNG/SVG links and image checks no longer apply to this PDF-only export. Use the new
manifest for export integrity. Recorded simulation source and environments are still required
to rerun or resume the study.

The untouched canonical source checkout remains at
`/Users/rongyuqing/.codex/worktrees/61d4/TopoNoise.jl`.
