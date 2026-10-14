# H+5 Z-boundary ambiguity diagnostics and circuit-aware control

**Date:** 2026-10-14  
**Status:** in-chat design approved; written specification awaiting review

**Purpose:** measure and expose the information-theoretic Z-decoding limit of
the current three-dimensional direct-projector H+5 boundary instrument, reduce
avoidable finite-size failures with a circuit-aware control decoder, and
prevent unsupported Z-threshold claims.

## Intent and success criteria

The completed `L=3,5,7` experiment shows increasing state-independent logical-Z
failure. Exhaustive follow-up diagnostics identify a concrete mechanism: a
single last-layer Z fault can be expanded by upper-boundary completion into a
long relative error, after which the unit-length final MWPM decoder may choose
the wrong homology class. More strongly, some equally likely single faults
produce the same final star syndrome and the same physically informative
upper-outcome parity while having different logical classes. No decoder using
only those observables can correct all such faults.

This work must:

1. turn that diagnosis into a deterministic, reproducible repository audit;
2. quantify the optimal first-order logical floor under the current observable
   set;
3. add a separate circuit-aware maximum-likelihood control decoder that removes
   avoidable, but not information-theoretically unavoidable, failures;
4. expose odd upper-boundary parity as an erasure/postselection diagnostic with
   its acceptance cost reported explicitly; and
5. block Z-threshold fitting whenever the single-fault correctness guard fails.

The existing scans, selections, checkpoints, figures, reports, and commits are
immutable historical evidence. New work goes in a new result directory. In
particular, `document/note.tex` must not be modified. Its pre-work SHA-256 is
`3a40710c2ab6bd43e53f5f7703b7c8716b209dc1e060a5092f17e0edce073f3a`, and final
verification must reproduce that hash.

## Established evidence

For all active outputs in the final H+5 layer, exhaustive Z-fault injection
gave:

| L | single-fault locations | current MWPM logical failures | mixed observable classes | optimal unavoidable failures |
|---:|---:|---:|---:|---:|
| 3 | 99 | 34 | 5 | 8 |
| 5 | 275 | 131 | 9 | 14 |
| 7 | 539 | 276 | 13 | 20 |

The final-layer location count is exactly `11L^2`. The unavoidable counts for
the tested sizes are `3L-1`; this observed formula is not promoted to a theorem
outside the audited sizes. All mixed classes occur in the odd upper-parity
sector. The even sector contains 27, 75, and 147 mechanisms respectively and
has no mixed logical class for this final-layer audit.

The clean upper branch `b` is uniform over even bit strings and the physical
outcome is `a = b xor u`, where `u` is the fault-induced upper syndrome. For a
fixed `u`, `a` is uniform over the parity class selected by `parity(u)`.
Consequently, the exact bit pattern of `a` contains no fault-location
information beyond its parity. A decoder must not treat the hidden clean branch
or the exact observed pattern as an error oracle.

## Scope

### Included

- A core exhaustive single-Z-fault signature enumerator for any H+5 dimensions
  and each canonical feedback protocol.
- A machine-readable ambiguity certificate with observable-class counts,
  logical-class multiplicities, parity stratification, and the uniform-prior
  maximum-likelihood floor.
- A circuit-aware final Z control decoder built from those signatures.
- A parity-erasure analysis that reports acceptance and conditional logical
  failure separately from unconditional deterministic performance.
- A Z-specific pre-fit guard that stops before optimization or bootstrap.
- Focused Julia tests, full Julia and Python regressions, provenance checks, and
  a new Markdown report under a new results directory.

### Excluded

- No modification of `document/note.tex` or its compiled note artifacts.
- No mutation or resumption of the original, four-size, or three-size
  production checkpoints.
- No replacement of the canonical unit-length MWPM result in existing data.
- No claim that postselection or the circuit-aware control creates a positive
  deterministic threshold.
- No new large Monte Carlo threshold scan while the single-fault guard fails.
- No redesign of the physical boundary circuit in this deliverable.
- No 4D implementation, remote/HPC submission, push, merge, or publication.

The physical boundary redesign is a separate research project because the
current upper branch carries only a global parity constraint. A scalable
deterministic repair must change the observable information or the preparation
architecture; decoder tuning cannot supply missing information.

## Observable signature

For one protocol and one injected single Z fault, replay all intermediate
feedback and the paired upper-boundary completion, but stop before final Z
recovery. Record:

- the active-output location `(step, output slot, physical qubit)`;
- the final full-code star syndrome;
- `upper_parity_odd = isodd(count(boundary.observed))`;
- the two logical-Z parity bits of the pre-recovery relative Z frame;
- whether the fault lies in the final H+5 layer; and
- the protocol and dimensions.

The decoder-visible key is

```text
(final star syndrome, upper_parity_odd)
```

The exact observed upper bit string, the clean branch, and sampled fault
location are retained only for invariant testing or failure diagnostics; they
are not decoder inputs.

The replay helper must share the production propagation, feedback, and
boundary code. It must not copy those algorithms into the audit. Final-layer
fault signatures must be identical across all four protocols because no
intermediate feedback runs after the final layer.

## Ambiguity certificate

Group signatures by decoder-visible key. Within each group, count the four
logical classes in deterministic order:

```text
(false,false), (false,true), (true,false), (true,true)
```

For equal first-order location probabilities, the optimal class is the largest
multiplicity, with the listed order as a reproducible tie break. A group's
unavoidable count is its total count minus its largest class count. The audit's
unavoidable count is the sum over groups.

The certificate must include:

- dimensions and protocol;
- encoder operation count and active-output fault count;
- total observable groups and mixed-logical groups;
- total unavoidable count and rate among enumerated single faults;
- the same counts split by even and odd upper parity;
- final-layer-only counts as a separately identified subset;
- every mixed group's syndrome indices and logical multiplicities;
- relevant source hashes, backend versions, and deterministic format version;
  and
- the declaration that the floor assumes equal single-location Z-fault priors.

The certificate is a scientific result, not merely debug output. JSON is the
canonical machine-readable form; a compact CSV table and Markdown summary are
derived from it.

## Circuit-aware control decoder

The control decoder is protocol- and geometry-specific. It contains the
signature table and the existing closed unit-length MWPM decoder.

For a final syndrome and observed upper parity:

1. obtain the physical recovery from the existing MWPM decoder;
2. look up the maximum-likelihood pre-recovery logical class for the observable
   key;
3. compute the MWPM recovery's two logical-Z parity bits;
4. XOR the required declared logical-Z representatives into the physical
   recovery so its logical class matches the selected class; and
5. validate that the adjusted recovery has exactly the requested star
   syndrome.

If the key is absent from the single-fault table, as can occur for multi-fault
events, return the existing MWPM recovery unchanged and mark the decision as a
fallback. The decoder must report the selected class, class counts, whether the
key was mixed, whether a logical-frame adjustment was made, and whether it
fell back.

This decoder is a finite-size control. It should reduce the demonstrated
final-layer single-fault failures from 34/131/276 toward the optimal floors
8/14/20, but those nonzero floors remain failures. The decoder therefore does
not satisfy the threshold-eligibility guard for the current instrument.

## Parity-erasure control

The replay outcome exposes the parity of the physically observed upper branch.
Since the clean branch is always even, odd observed parity is a physically
available flag. The diagnostic reports, for each protocol and scan point:

- total shots;
- accepted even-parity shots;
- rejected odd-parity shots;
- acceptance probability with a binomial interval;
- unconditional logical-Z failure; and
- logical-Z failure conditional on acceptance, with its own numerator and
  denominator.

Rejected shots never disappear from denominators silently. Postselected curves
must be labeled `parity-erasure diagnostic`, and no deterministic threshold may
be inferred from them. Version one does not fit a postselected threshold,
because acceptance may decay with lattice size. If a point has zero accepted
shots, its conditional rate and interval are serialized as missing values,
never as zero.

The control may use a small, new diagnostic run if needed. It must not rewrite
or reinterpret the completed production CSV because those artifacts do not
store enough per-shot boundary information.

## Scientific guard

A Z-threshold candidate is eligible for optimization and bootstrap only when
all of the following hold:

1. at least four lattice sizes are present;
2. every included geometry and protocol has zero mixed observable groups;
3. every included geometry and protocol has zero unavoidable single-fault
   count;
4. exhaustive replay with the proposed decoder has zero single-fault logical-Z
   failures; and
5. the ordinary existing data-completeness and fit-window guards pass.

The ambiguity checks precede fitting. On failure, the result is
`diagnostic_only`, the reason lists the failed counts, and requested bootstrap
replicates execute zero times. Existing generic fit behavior remains unchanged;
the new guard wraps Z-boundary analyses rather than weakening historical APIs.

## Testing and verification

Implementation follows test-driven development. Each behavior receives a
failing test before production code.

Focused tests must cover:

- deterministic enumeration of all 11 active outputs per H+5 cell;
- final-layer counts 99/275/539 for `L=3/5/7`;
- the exact final-layer ambiguity results in the evidence table;
- no mixed final-layer class in the even-parity sector for `L=3/5/7`;
- equality of final-layer signatures across all four protocols;
- invariance of the relative completed frame and parity under multiple clean
  even branches;
- maximum-likelihood tie breaking and unknown-key fallback;
- syndrome preservation after logical-frame adjustment;
- circuit-aware final-layer failures equal the certified floor;
- parity-erasure numerator, denominator, and interval handling, including zero
  accepted shots;
- the Z-fit guard stopping before optimization and bootstrap; and
- stable JSON/CSV serialization and provenance.

After focused tests pass, run the complete Julia suite and complete Python
analysis suite under the recorded single-thread numerical environment. Final
verification must show a clean worktree apart from intended new files, verify
all historical result hashes that the existing delivery verifier covers, and
recompute the exact `document/note.tex` hash stated above.

## Deliverables

Implementation may add focused source and test files under `src/toric3d/` and
`test/toric3d/`, plus a new directory:

```text
results/toric3d-h5-z-boundary-ambiguity-2026-10-14/
```

That directory contains the audit JSON, derived CSV, concise Markdown report,
commands/environment provenance, and any small diagnostic-control outputs. It
must link to the immutable completed preliminary results instead of copying or
editing them.

The final report must distinguish three conclusions:

1. the current decoder has avoidable circuit-model mismatch;
2. the current observable set also has an unavoidable first-order ambiguity;
   and
3. a deterministic positive Z threshold requires a later physical boundary or
   dimensional redesign.
