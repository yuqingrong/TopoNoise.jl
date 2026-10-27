# H+6 one step Z feedback results

Four paired protocols on the legacy H+6 circuit: final-only, seed-only, causal one-step, and complete-window recovery. Production contains 50,000 shared fault records (200,000 protocol outcomes), with T=Ly=Lz=L and L=3,5,7,9,11.

At L=11 and p=0.02, logical-state failure rates are 74.8%, 75.6%, 42.6%, and 0.3%, respectively. The primary one-step scaling model fails its checks; no threshold or critical exponent is established. The complete-window arm here is distinct from the later two-slice protocol, which also applies its local recovery after the last transfer.

## Folder guide

- This README summarizes the protocol and results.
- [data](data/) contains final analysis tables, production counts, and metadata.
- [figures](figures/) contains the final PDF plots.
- [archive](archive/) contains one compressed study archive and its integrity manifest.

## Figures

- **[Combined X and Z one-step fit panels (PDF)](../toric3d-h6-one-step-2026-10-15/figures/one-step-feedback-x-z.pdf)** — the main side-by-side figure. X shows its exploratory restricted-window fit; Z shows its rejected primary fit. The separate channel data and fit-status labels are unchanged.
- [z-one-step-four-protocols (PDF)](figures/z-one-step-four-protocols.pdf)
- [z-one-step-local-fit (PDF)](figures/z-one-step-local-fit.pdf)
- [z-one-step-paired-differences (PDF)](figures/z-one-step-paired-differences.pdf)

## Data

- [marginals.csv](data/marginals.csv)
- [paired-differences.csv](data/paired-differences.csv)
- [report.json](data/report.json)
- [production/metadata.toml](data/production/metadata.toml)
- [production/raw.csv](data/production/raw.csv)
- [selection.json](data/selection.json)

The data files retain their original bytes. Paths and hashes inside JSON/TOML files describe
the original run layout and remain provenance records. The presentation redraw manifest,
where present, also records the original image formats.

## Archive and reproduction

[study-records.tar.gz](archive/study-records.tar.gz) preserves 1,077 original files:
the full scientific README, raw data, metadata, checkpoints, pilot runs, scripts, validation
records, and all original PDF plots. The 6 non-PDF figure files were omitted.
[manifest.json](archive/manifest.json) records every archived and omitted file's original
SHA-256 hash, the published file mapping, and checksums for this cleaned export.

Archive members retain their original `results/toric3d-h6-z-one-step-2026-10-15/` paths. Extract into a separate
reproduction checkout to restore that layout before using the original scripts or checkpoints.
The original README and artifact manifests inside the archive are historical records;
their PNG/SVG links and image checks no longer apply to this PDF-only export. Use the new
manifest for export integrity. Recorded simulation source and environments are still required
to rerun or resume the study.

The untouched canonical source checkout remains at
`/Users/rongyuqing/.codex/worktrees/61d4/TopoNoise.jl`.

The top-level `analyze.py` is retained unchanged because older full-prefix analysis scripts import it.
