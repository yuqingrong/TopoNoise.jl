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

The primitive, mutation, pilot replay and frozen-provenance results in this
section and the next section belong to the original source snapshot. Their
byte-preserved evidence is now under `pre-optimization/verification/`; final
source results must be verified separately after the parity optimization.

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
byte-identical inputs after finishing. Production replay was still pending at
that original-source validation stage; final-source replay evidence appears
below.

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
counts or fit acceptance. The requested production-point replays were pending
at that stage; final-source replay evidence appears below.

## Exact parity optimization transition

`verification/check_transition.py` is a read-only transition audit of
the rerun. It compares the archived original selection with the
new selection, preserving channel grids, shard configurations and seeds,
sample targets, protocol/decoder settings, fit rates, bootstrap settings and
the numerical-source manifest. It requires the new selection to bind the
original selection hash and statistical freeze time. The original and final
simulation manifests must contain the same 56 paths, with exactly
`src/toric3d/decoder.jl` changed. It validates both simulation-source archives,
the current source bytes and the complete 1,364-file baseline archive.

For every one of the nine archived pilot/production authoritative
checkpoints, the script requires exactly one new immutable snapshot at the
same generation. It verifies checkpoint checksums and compares the entire
configuration, exact point-seed list, protocol and every state's completed
records, both per-logical counts, both failure counts and all four joint bins.
Timing and runtime fields are ignored for this scientific-count comparison;
source hashes are separately bound to the appropriate original/final manifest.
Missing or ambiguous new generations fail the audit. It writes JSON only to
standard output and never samples, decodes, resumes or modifies an input.

Seven focused tests passed, including rejection of scientific-design changes,
any source delta except the one allowed file, changed scientific counts,
missing/ambiguous generations, checksum changes and input modification. The
actual archived 1,364-file manifest also passed byte verification. The full
transition audit then passed 7,337 assertions, recorded in
`performance-transition-audit.json`. It exactly reproduced 2,220 excluded
pilot records and 55,250 production records (35,000 Z and 20,250 X) across
all nine archived checkpoints, including zero completed records in two early
X shards. This establishes exact aggregate reproduction at those generations,
not identity of every individual decoder recovery. The archived records remain
verification evidence and are never pooled with the rerun. The independent
final review rechecked all nine reported baseline/new snapshot hashes and
scientific-count comparisons, the original/final selection hashes, and both
source-archive hashes. Exactly `src/toric3d/decoder.jl` differs between the
56-file source manifests. The two final-source full-point replays are recorded
below.

## Final-source X production replay

The completed final-source `production-x-small` point L=5, p=0.05 passed the
independent replay for all 1,000 paired records, with exit code 0 and 84,409
assertions. The verifier checked all 5,000 full-prefix checkpoints, including
1,000 final-transfer corrections, and independently scored 2,000 final
outcomes. Four-bin counts were exactly `[293, 511, 61, 135]`; final-only and
full-prefix failure counts were `[646, 196]`, logical-1 counts `[417, 118]`,
and logical-2 counts `[439, 137]`. Every count matches checkpoint generation
442.

`replay-production-x-L5-p005.toml` and its `.log` preserve the result. The
verifier bound the final optimized 56-file source manifest, managed runtime,
payload checksum, immutable snapshot, raw and metadata mirrors, and confirmed
all input files remained unchanged. The archived original evidence was not
modified. This X result certifies the selected point only; it does not
establish fit acceptance or the full transition audit.

## Final-source Z production replay

The completed final-source `production-z` point L=11, p=0.0325 passed the
independent replay for all 1,000 paired records, with exit code 0 and 145,591
assertions. All 11,000 full-prefix checkpoints, including 1,000 final-transfer
corrections, and 2,000 final outcomes were independently checked. Four-bin
counts were exactly `[165, 469, 121, 245]`; final-only and full-prefix failure
counts were `[714, 366]`, logical-1 counts `[456, 263]`, and logical-2 counts
`[485, 230]`. Every count matches checkpoint generation 706.

`replay-production-z-L11-p00325.toml` and its `.log` preserve the result. The
final optimized source manifest, managed one-thread runtime, checkpoint
checksum, immutable snapshot, raw CSV and metadata all passed the same strict
binding checks as X, and all input bytes remained unchanged. Both requested
final-source point replays are complete: together they checked 2,000 paired
records, 16,000 prefix checkpoints, 4,000 final outcomes and 230,000
assertions. No additional replay was started. The saved replay reports and
logs were rechecked against the completed input checkpoints, immutable
snapshots, raw/metadata hashes, unchanged verifier hash and final source/runtime
manifests during the final evidence review.

## Final suite provenance and review scope

The saved final `full-suite.log` contains 189 printed test summaries. An
independent recount found 96,525 passing assertions and zero failures, errors
or broken tests, matching `full-suite-summary.json` and its recorded exit code
0. The summary's source and test manifest hashes match the current manifests;
all 56 simulation files and 71 test files match their archived and current
bytes. This connects the saved suite result to the same optimized source used
by both independent production-point replays.

The full-suite log has SHA256
`279f22bb4294ce76cbc9399bd7a72661b01b23b09be5d7d7293a772541edb895`,
now recorded in the summary. The parent task's retained execution context
identifies process session 44900 and exit code 0, but the verbatim launch
command and launch environment were not independently recovered. The summary
therefore supplies a clearly labeled reproduction command using the saved
managed environment and `test/runtests.jl`; it does not present that command as
recovered launch evidence. No expensive test or simulation was rerun for this
review.

No scientific-count, source-binding or replay discrepancy was found in this
bounded verification review. The transition and two requested point replays
are complete. The overall production grid still has unfinished large-X points
at this review stage, and analysis has not yet been generated. Final grid
completeness, the generated analysis, fit guard outcomes and figure contents
remain outside this review and must be checked after those artifacts exist.
