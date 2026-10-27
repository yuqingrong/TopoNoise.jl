# Independent review of the full-prefix protocol

This bounded review covers the new H+6 full-prefix replay and its independent
point verifier. It does not establish a threshold or certify production results
that have not yet been replayed. The approved intervention is complete physical
recovery on all available slices Q0 through Qr after every transfer r=1,…,T,
including T. There is no separate r=0 seed correction.

## Protocol source review

`src/toric3d/h6_full_prefix.jl` implements separate X and Z public replay
functions using the established growing-prefix decoder families. Both functions
propagate a transfer, apply its prefix recovery when feedback is enabled, and
repeat through T without an early skip. After the loop, each calls the same
full-code outcome function exactly once for both arms. The family geometry and
every prefix depth are validated before propagation.

The X helper measures every physical plaquette of the current prefix, supplies
only that syndrome to BP+OSD, applies the returned syndrome-consistent recovery
through the prefix map, and checks all prefix plaquettes and frontier checks.
The Z helper similarly measures all current-prefix physical stars, uses
unit-weight data-edge MWPM, validates its recovery syndrome, and checks the
complete prefix after applying the correction. Neither decoder uses the hidden
fault record or logical parities as an input.

The X frame starts from the ideal logical |00> seed, with faults at the
13 L^3 participating transfer outputs. The Z frame starts with the recorded
3 L^2 noisy seed-H outputs for logical |++> and then executes transfer 1
without a seed-only repair. Its transfer record adds 13 L^3 output locations.
Measurements, corrections and idle data remain ideal. These differences in
seed preparation, noise locations and decoder must remain explicit when
comparing X and Z results.

The final feedback syndrome is already clear, and both final decoders return
zero recovery for zero syndrome. Thus the common final outcome does not
silently repair a missed last checkpoint. It still scores both residual
logical parities, retaining undetectable logical failures. The final-only arm
has no physical correction history and receives its sole recovery at the end.

The earlier two-slice argument that X parities freeze after the first window
does not apply. Prefix recovery at later transfers can revisit Q0, which
contains the logical-Z scoring cycles. This study must score both logicals
after the final prefix recovery and must not infer its later-time behavior
from the earlier two-slice results. The old `replay_toric3d_*_strength` protocol
also has different checkpoint/seed timing and is not substituted here.

## Independent replay checks

`verification/replay_point.jl` never imports the study runner. It regenerates
the complete size-major, rate-minor point-seed list from repeated scalar
MersenneTwister UInt64 draws, discarding duplicates, then regenerates the
selected point's Xoshiro fault stream from scalar draws. X consumes its
established zero-probability Z draw after each X draw. Z draws seed faults
first and then one Z bit per participating transfer output. The verifier does
not call either production fault sampler.

For each regenerated record, it calls the two public replay functions for the
appropriate channel, checks fault immutability, and reconstructs both frames
using explicit binary H/CNOT propagation. Feedback corrections are taken from
the recorded history and physically XORed onto a separately constructed
coordinate map. The verifier builds a fresh encoded prefix at each r and
forms its complete physical-check incidence matrix independently of production
prefix builders and decoder matrices. It checks:

- Exactly T ordered feedback entries, including T, with no r=0 entry; no
  correction entries in the final-only arm.
- Recovery columns equal the complete Q0:Qr prefix; no future-slice correction.
- Every measured prefix syndrome matches independent propagation, every
  recovery reproduces that syndrome, and every prefix is clear afterward.
- The X frontier and recorded future-data fields agree with independent
  measurements; Z entries have `kind=:prefix` and the correct physical support.
- Reconstructed final frames equal the replay frames in both Pauli channels.
- Final physical syndrome, recovery syndrome, residual, both logical parities
  and their state-failure union agree with separate matrix products.
- All four paired joint bins, each arm's two logical counts, and each failure
  marginal exactly match the selected checkpoint state.

It also binds the result to the checkpoint's canonical payload checksum and
byte-identical immutable snapshot, the complete current simulation source
manifest, the declared decoder settings, Julia/Python/package versions and
one-thread/Null-backend environment. Raw CSV and metadata mirrors are checked
against the checkpoint, and all input files are checked for unchanged bytes
after replay. Verification output must be outside the input directory.

The verifier deliberately shares the public circuit schedule, physical check
definitions, logical-support definitions and decoder implementations with the
simulation. It therefore independently audits the runner, sampling,
propagation, prefix mapping and scoring plumbing, but is not a separate
implementation of BP+OSD or MWPM and cannot by itself prove the underlying
geometry definitions.

## Validation status

Focused primitive validation passed 145/145 assertions with exit code 0.
`verifier-primitives.log` records checks of exact seed values, independent
X/Z fault records and subsequent RNG position, direct gate propagation,
coordinate maps and parity products. Cases include T=1 and rectangular
geometries, with p=0, 0.05 and 0.3. The tiny history/mutation test passed 52/52
assertions with exit code 0 (`verifier-history-mutations.log`). Valid X/Z
histories and both final-only frames passed; deliberately missing final
entries, spurious r=0 entries, future recovery columns, wrong recoveries,
incorrect final residuals/logicals and altered X future-data or Z entry-kind
fields were rejected.

The completed excluded Z pilot point L=3, p=0.03 was independently replayed:
10 paired records, 30 prefix checkpoints, 20 final outcomes and 979 assertions
all passed, with joint counts `[5, 2, 0, 3]` exactly matching its checkpoint.
`replay-pilot-z-L3-p003.toml` binds the counts and checks to the exact input,
source and runtime; the matching `.log` records exit-success output. The
completed excluded X pilot point L=3, p=0.035 also passed: 10 paired records,
30 prefix checkpoints, 20 final outcomes and 1,009 assertions, with joint
counts `[6, 4, 0, 0]`. Its evidence is
`replay-pilot-x-L3-p0035.toml` and the matching `.log`. Both replays checked
byte-identical inputs after finishing. Production replay remains pending; no
production simulation or replay has been started by this reviewer.

## Frozen design and provenance review

The frozen `selection.json` has SHA256
`a27f27669fb7e2b5d2aa534a9b605f79c4b111a086eeaead323f5cade9fdefdc`.
An independent read-only audit passed 604 assertions; its evidence is
`frozen-provenance-review.json` and `.log`. The six shards partition the two
five-size, seven-rate channel grids exactly, with no duplicated or missing
point: 70 points, 70,000 planned fault records and 140,000 paired outcomes.
Both fit-rate sets, 200 bootstrap replicates and distinct bootstrap seeds were
recorded before sampling. All six actual production checkpoints match their
frozen configurations, source/runtime/protocol declarations and regenerated
point seeds, and their recorded start times follow the production freeze.

All three excluded pilot datasets, including the Z refinement, match the
frozen raw/metadata/checkpoint hashes, payload checksums and immutable
snapshots. Their 60 points contain 2,220 completed records; all finish before
the production freeze and use identical simulation sources and runtime.
Every pilot raw joint bin and per-logical count matches its checkpoint.
Independently regenerated scalar MersenneTwister seeds match all pilot and
production declarations: the 60 pilot point seeds and 70 production point
seeds are individually unique and disjoint. This checks exclusion by record
stream even where pilot and production grid coordinates overlap.

The simulation archive contains exactly the 56 declared files, the numerical
archive exactly its three frozen fitting sources, and the analysis archive
exactly its nine manifest entries. Every archived member matches its declared
SHA256 and the current file; no archive was extracted or modified. The
numerical adapter reads measured full-prefix counts and reuses the existing
fitter and guard policy. Input validation and presentation changes are
separately archived, as the frozen selection permits. This review validates
the frozen design and provenance; it does not certify unfinished production
counts or fit acceptance. The requested full production point replays remain
on hold until their entire input shards are complete.
