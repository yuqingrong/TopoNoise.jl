# H+6 two slice recovery results

Complete syndrome-consistent recovery on the newest two data slices after every transfer, including the last, paired with a final-only control on identical fault records. The H+6 circuit is unchanged; X and Z are separate experiments.

All 70 production points completed with 1,000 shared records each: 35,000 records per channel and 140,000 protocol outcomes. At L=11, the X failure rate at p=0.05 changes from 72.3% (final-only) to 17.4% (two-slice); Z at p=0.03 changes from 76.9% to 20.0%.

Both primary threshold fits remain diagnostic only. The drop-L=3 sensitivity plots are included; they do not establish robust thresholds. Under the ideal-idle model, the X outcome is fixed after the first corrected window, so its curve is not a noisy-memory threshold along T.

## Folder guide

- This README summarizes the protocol and results.
- [data](data/) contains final analysis tables, production counts, and metadata.
- [figures](figures/) contains the final PDF plots.
- [archive](archive/) contains one compressed study archive and its integrity manifest.

## Figures

- [two-slice-comparison (PDF)](figures/two-slice-comparison.pdf)
- [two-slice-paired-differences (PDF)](figures/two-slice-paired-differences.pdf)
- [x-two-slice-drop-L3-fit (PDF)](figures/x-two-slice-drop-L3-fit.pdf)
- [x-two-slice-local-fit (PDF)](figures/x-two-slice-local-fit.pdf)
- [z-two-slice-drop-L3-fit (PDF)](figures/z-two-slice-drop-L3-fit.pdf)
- [z-two-slice-local-fit (PDF)](figures/z-two-slice-local-fit.pdf)

## Data

- [marginals.csv](data/marginals.csv)
- [paired-differences.csv](data/paired-differences.csv)
- [report.json](data/report.json)
- [production-x-L11-a/metadata.toml](data/production-x-L11-a/metadata.toml)
- [production-x-L11-a/raw.csv](data/production-x-L11-a/raw.csv)
- [production-x-L11-b/metadata.toml](data/production-x-L11-b/metadata.toml)
- [production-x-L11-b/raw.csv](data/production-x-L11-b/raw.csv)
- [production-x-L11-c/metadata.toml](data/production-x-L11-c/metadata.toml)
- [production-x-L11-c/raw.csv](data/production-x-L11-c/raw.csv)
- [production-x-L9/metadata.toml](data/production-x-L9/metadata.toml)
- [production-x-L9/raw.csv](data/production-x-L9/raw.csv)
- [production-x-small/metadata.toml](data/production-x-small/metadata.toml)
- [production-x-small/raw.csv](data/production-x-small/raw.csv)
- [production-z/metadata.toml](data/production-z/metadata.toml)
- [production-z/raw.csv](data/production-z/raw.csv)
- [selection.json](data/selection.json)

The data files retain their original bytes. Paths and hashes inside JSON/TOML files describe
the original run layout and remain provenance records. The presentation redraw manifest,
where present, also records the original image formats.

## Archive and reproduction

[study-records.tar.gz](archive/study-records.tar.gz) preserves 1,609 original files:
the full scientific README, raw data, metadata, checkpoints, pilot runs, scripts, validation
records, and all original PDF plots. The 20 non-PDF figure files were omitted.
[manifest.json](archive/manifest.json) records every archived and omitted file's original
SHA-256 hash, the published file mapping, and checksums for this cleaned export.

Archive members retain their original `results/toric3d-h6-two-slice-2026-10-15/` paths. Extract into a separate
reproduction checkout to restore that layout before using the original scripts or checkpoints.
The original README and artifact manifests inside the archive are historical records;
their PNG/SVG links and image checks no longer apply to this PDF-only export. Use the new
manifest for export integrity. Recorded simulation source and environments are still required
to rerun or resume the study.

The untouched canonical source checkout remains at
`/Users/rongyuqing/.codex/worktrees/61d4/TopoNoise.jl`.
