# Guarded-fit input set

Run `build_inputs.py` from the repository root after the complete production
artifacts exist. It writes the eight CSV files named explicitly in the report
and records their source and derived SHA-256 hashes in `manifest.json`.

The four X files are byte copies of the complete `L=5,7,9,11` production
CSVs. The four Z files retain the complete rectangular seven-probability grid
through `p=0.045` for `L=9,13,17,21`. The original ten-rate smaller-size Z
artifacts and resumable checkpoints are preserved under `production/`; this
directory only fixes the immutable inputs used by the reported fit.
