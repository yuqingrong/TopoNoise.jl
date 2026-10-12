# Production artifact audit

`audit_production.py` independently reconstructs every protocol marginal,
Wilson interval, and paired contingency count from the saved eight-outcome
joint tables. It also checks each raw row against its checkpoint state,
requires zero invalid recoveries, verifies frozen seeds, decoder settings,
versions, and timings, follows the assembled `L=21` checkpoint lineage back to
both source shards, and hashes the complete production and guarded-fit inputs.
Duplicate rows and any analysis-input selection or ordering mismatch abort the
audit.

Run it with the repository's managed Python after `analysis-input/manifest.json`
has been generated. The immutable result is written to
`production-artifact-audit.json`, replacing any earlier in-progress audit.
