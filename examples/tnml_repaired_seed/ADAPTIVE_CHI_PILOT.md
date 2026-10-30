# Adaptive bond-dimension pilot design

This document proposes an isolated pilot around the existing repaired-seed TN
decoder.  It does not change `src/`, `model.jl`, `run.jl`, the existing HPC
wrappers, or the scoring of any completed study.  The pilot must live in a
fresh output root and call the same `RepairedSeedTN.decode_record` function
with a sequence of `Toric3DZTNSettings` values.

The purpose is to learn whether most shots can be decided reliably at a small
bond dimension while rare ambiguous shots escalate.  A converged adaptive
result remains an approximate finite-bond TN result.  It must not be called
exact ML unless an independent exact calculation is available for that shot.

Policy revision (2026-10-30): new sampling manifests use
`tnml-adaptive-stopping-v2`. An explicitly numerical-invalid TN result retries
the next chi, clearing the convergence window and provisional anchor decision.
Acceptance requires three consecutive valid levels after any invalid level.
This revises the original blanket invalid-result termination rule below;
exceptions, malformed results, resource failures, and invalidity at maximum chi
remain unknown. Existing frozen studies retain their original policy and
results. The six-history ladder diagnostic in `adaptive_chi_pilot/hpc` checks
this change without granting a production calibration gate.

## Approaches considered

1. **Independent fixed-chi studies followed by an offline comparison.** This
   reuses `run.jl` unchanged, but every history is recomputed at every bond
   dimension.  The current exact mode also halts at the first low-chi class
   disagreement, so it is poorly suited to calibrating an adaptive rule.
2. **A new checkpointed per-shot adaptive wrapper (recommended).** Generate a
   frozen history corpus once, build decoders lazily for a fixed chi ladder,
   and run each history until it passes a predeclared convergence rule or
   reaches the ladder/resource limit.  This reuses the scientific decoder
   while avoiding work at larger chi for easy histories.
3. **Change chi during a single contraction.** This could reuse intermediate
   tensors, but it changes the production algorithm and couples the stopping
   decision to internal implementation details.  It is too invasive for the
   present evidence-gathering step.

The pilot should use approach 2.  No adaptive result should enter a threshold
figure until the `L=3` exact audit below passes and the complete pilot report
is frozen.

## Fixed numerical policy

Use one cutoff at every level so that the comparison isolates bond dimension:

```text
cutoff = 1e-18
chi ladder for L=3: 32, 64, 128, 240, 480, 960
chi ladder for L=5: 32, 64, 128, 240, 480
negative-probability tolerance = 1e-10
ML tie tolerance = 1e-9
```

Julia uses one application thread.  The HPC pilot may request four CPU cores
and set BLAS to four threads, because the earlier isolated BLAS probe found
four better than eight.  Thread count, BLAS vendor, package lockfiles, source
hashes, and CPU model are part of the frozen manifest.  A one-thread numerical
replay subset must confirm that BLAS threading does not change accepted
classes or move posteriors beyond the convergence tolerance.

## Per-shot convergence rule

Let `q_k` be the normalized four-class posterior and `z_k` the log joint
syndrome probability at the `k`th chi.  An output at level `k` is eligible
only if it is valid, finite, normalized within `1e-8`, and has no entry below
`-1e-10`.  Define

```text
d_k = max(abs(q_k - q_(k-1)))
g_k = abs(z_k - z_(k-1))
m_k = largest(q_k) - second_largest(q_k)
B_k = max(d_k, d_(k-1))
G_k = max(g_k, g_(k-1))
```

Acceptance requires all of the following:

1. At least three chi levels have completed, so there are two successive
   comparisons.  The earliest possible accepted level is `chi=128`.
2. The unique maximizing class is identical at the last three levels.
3. `B_k <= 1e-6` and `G_k <= 1e-6`.
4. `m_k > max(1e-8, 4 * B_k)`.  This leaves a factor-two safety margin over
   the largest possible two-class margin shift implied by `B_k`.
5. The newest differences do not increase sharply:
   `d_k <= max(1e-12, 1.25*d_(k-1))` and
   `g_k <= max(1e-12, 1.25*g_(k-1))`.  This rejects an apparent plateau
   followed by renewed drift.  If a preceding difference is zero, the newer
   difference must be at most `1e-12`.

The accepted posterior and logical class are those at level `k`, never an
average across levels.  The rule is intentionally conservative near ties.
If any condition fails, evaluate the next chi.  If the maximum chi, time, or
memory limit is reached, the shot is **unknown**.  An invalid tensor result,
exception, missing checkpoint, or interrupted process is also unknown.
Unknown shots remain in the requested denominator and contribute the usual
failure-rate interval

```text
[known failures / requested shots,
 (known failures + unknown shots) / requested shots].
```

Successive-chi agreement is evidence of numerical stability, not a proof of
the exact posterior.  The exact audit and forced high-chi anchors are what
test whether this empirical stopping rule is safe.

## Checkpoint and identity model

Each task owns one output directory and locks it during execution.  Before
any contraction it writes the complete raw-fault histories and their stable
record hashes.  A shot directory then contains one immutable file per
attempted chi and one final decision file:

```text
task/
  metadata.toml
  input-corpus.toml
  shots/shot_0001/
    chi_0032.toml
    chi_0064.toml
    chi_0128.toml
    ...
    exact.toml             # L=3 audit only
    decision.toml
  summary.toml
  progress.toml
```

Every chi checkpoint records the record ID, chi, cutoff, posterior, logical
class, log joint syndrome, validity, elapsed seconds, peak RSS observation,
and hashes of all scientific sources.  `decision.toml` records every tested
level, all `d_k`, `g_k`, margins, the first accepted chi or the unknown reason,
and whether the shot was selected as a forced high-chi anchor.  Resume accepts
only an identical manifest, record hash, runtime identity, and source hashes.

Decoder objects are cached by chi inside a task.  A shard therefore pays the
decoder-construction cost at most once per level, while every saved history is
still scored independently.

## Predeclared p-stratified selection

Replay cases and forced high-chi anchors are selected independently within
each probability stratum.  After the complete corpus is frozen but before any
contraction, form the UTF-8 byte string

```text
tnml-adaptive-chi-v1|L=<decimal L>|p=<four-decimal p>|record_id=<lowercase hex>
```

and compute its SHA-256 digest.  Interpret the 32 digest bytes as one unsigned
256-bit **big-endian** integer (the first digest byte is the most significant).
Within each `(L,p)` stratum, sort by `(digest_integer, record_id)` and take the
first predeclared number of records.  The lexical record-ID tie break makes the
order total even under a digest collision.  The manifest stores the selected
record IDs and their ranks.  Language-specific `hash` functions and selection
after observing runtime or decoder output are forbidden.

## L=3 exact calibration pilot

Use `p = 0.006, 0.010, 0.020, 0.030, 0.035, 0.0425`, 100 histories per
value, and exact dense-frontier references for all 600 histories.  Use
20-shot shards, giving 30 Slurm array elements.  Every history runs the
adaptive ladder; exact calculation is recorded even when the adaptive class
disagrees.  Unlike the current validation runner, the pilot must not stop at
the first disagreement, because disagreement frequency is the quantity being
measured.  Five hash-ranked records per p stratum are also replayed from a
fresh output directory, giving 30 predeclared determinism checks.

The rule passes this calibration gate only if:

- all accepted classes are members of the exact maximizing set under the
  existing `1e-9` tie tolerance;
- accepted posterior L-infinity error is at most `1e-5` and absolute log-joint
  error is at most `1e-5` for every accepted history;
- no invalid/exception result is silently converted into a class;
- results and accept/unknown decisions are identical in an independent replay
  of all 30 predeclared histories; and
- the report gives the exact disagreement count, unknown count, chi stopping
  histogram, runtime by chi, and one-sided 95% exact Clopper-Pearson upper
  bounds on
  the conditional false-accept rate both overall and separately for every p
  stratum.

Zero observed false accepts is not evidence of a zero error rate.  For
example, zero errors among 600 accepted histories gives an approximate 95%
upper bound of 0.5% by the rule of three.  The report must state the actual
accepted denominator and a Clopper-Pearson upper bound.

## L=5 runtime and escalation feasibility smoke

Begin with ten histories at each of the same six probabilities: 60 one-shot
Slurm array elements.  Submit them as two immutable 30-element batches with
array concurrency at most 20, so the cluster's 50-submitted-job QoS limit is
never exceeded.  One-shot tasks prevent a difficult history from hiding the
runtimes of other histories.  Request four CPUs per task, one Julia
application thread, four BLAS threads, and at most 16 GiB RSS.

Predeclare the first four hash-ranked histories in every p stratum as
**forced anchors**, giving 24 anchors (40% of this small smoke).  Anchors
continue through `chi=480` even if they meet the stopping rule earlier.  Their
early accepted class and posterior are compared with the highest completed
chi, but the high-chi result is not used to rewrite the original stopping
decision.  This measures false early convergence without selecting convenient
histories after seeing the answer.

Suggested per-contraction wall limits are:

| chi | limit |
|---:|---:|
| 32 | 120 s |
| 64 | 300 s |
| 128 | 900 s |
| 240 | 1,800 s |
| 480 | 3,600 s |

The pilot is successful only if all of the following are reported without
discarding unknowns:

- at least 80% of non-anchor histories accept by `chi=128`;
- no forced anchor changes its maximizing class after its early acceptance;
- the median total adaptive wall time is below 60 seconds and the 90th
  percentile is below 5 minutes;
- peak RSS remains below 16 GiB; and
- every timeout, invalid result, exception, and unfinished anchor remains an
  explicit unknown.

The report gives anchor class-change counts and one-sided 95% exact
Clopper-Pearson upper bounds both overall and per p stratum. Even zero changes
among 24 anchors leaves a broad upper bound, and four anchors per stratum give
very little p-resolved precision.  This phase is therefore only a runtime and
escalation feasibility smoke.  Passing it authorizes a larger preregistered
anchor-validation pilot, not a production threshold run.  If it fails, no
larger Monte Carlo run should be launched with this adaptive scheme.

## Slurm layout

Create a new pilot directory rather than editing the existing HPC study:

```text
examples/tnml_repaired_seed/adaptive_chi_pilot/
  runner.jl
  manifest.py
  run_task.py
  supervise.py
  array.sbatch
  collect.py
  test_manifest.py
  test_collect.py
```

The manifest enumerates every task, history seed, p value, chi ladder,
resource limit, exact-audit flag, forced-anchor selection rule, and expected
output path.  Array concurrency stays at or below the user's `MaxJobsPU=20`.
Submit the 30-element L3 exact array first.  Freeze and audit its report before
creating the separate two-batch, 60-element L5 feasibility smoke.  Collection
refuses to publish point estimates when requested histories are missing and
always emits known-failure intervals plus the complete unknown-reason table.

## What remains before implementation

This design deliberately does not create or submit numerical jobs.  The next
implementation step is a new, tested pilot wrapper with synthetic unit tests
for the stopping rule, manifest immutability, resume identity, unknown
propagation, forced-anchor selection, and exact-audit collection.  Only after
those tests pass should a three-history L3 smoke run be authorized.

