# Adaptive-chi pilot control plane

This directory contains the tested phase-1 sampling plan and phase-2 execution
control plane described in `../ADAPTIVE_CHI_PILOT.md`. Nothing here submits a
Slurm job. The phase-2 tests use a fake executor and never run a contraction.

`manifest.py` creates an immutable sampling-plan manifest. It freezes the six
probability strata, every history seed, task layout, chi ladders, stopping
policy, resource limits, and the exact hash-ranking rule. `plan_task.py` only
reads that manifest and prints a dry-run task plan. `core.py` implements the
pure stopping, ranking, checkpoint-validation, and unknown-accounting rules.

`execution_manifest.py` freezes a generated history corpus, the selected
record IDs and ranks, source hashes, and runtime identity into a separate
immutable execution manifest. It also freezes the exact-audit thresholds and
one-sided 95% exact Clopper-Pearson convention. The phase-1 sampling manifest
is never mutated.

`adaptive_runner.py` applies the pure stopping rule after each immutable chi
checkpoint and resumes only a validated contiguous prefix. Invalid, timed-out,
memory-limited, or unavailable exact results remain unknown. Manifest-defined
per-chi wall limits are enforced by the runner. The frozen 16 GiB memory cap
must be enforced by the cluster supervisor; a completed adapter result
reporting a larger peak RSS is rejected as memory-limited. `julia_adapter.jl`
is a thin optional adapter: it authenticates the frozen record and calls the
existing `RepairedSeedTN.decode_record` and `TNMLRepairedRunner.exact_reference`
entry points. It does not reimplement decoder mathematics.

Stopping policy v2 retries an explicitly numerical-invalid TN level at larger
chi, clearing the convergence window and any provisional anchor acceptance.
Acceptance requires three consecutive valid levels after the most recent
invalidity. Exceptions, malformed outputs, timeouts, and memory failures remain
terminal unknowns. The new policy is frozen in new manifests; earlier frozen
studies retain their original evidence and must not be reused under this policy.

`task_runner.py` adds an optional persistent executor. One process owns every
selected record for one complete `(L, p)` stratum, including all of that
stratum's manifest shards. It loads Julia, the backend, and the fixed-L context
once, then caches one decoder for each chi at the frozen `(p, cutoff)`. Records
remain serial within the process and still pass through `adaptive_runner.py`,
so checkpoint validation, stopping, exact audits, resume, and unknown
propagation are identical to the one-shot fallback.

Every exact and TN checkpoint binds `executor_mode`, `wall_time_scope`, and
`decoder_cache_scope`. Persistent checkpoints also retain cache-hit, decoder
construction, and decode timing evidence. A resumed record must use the same
executor/timing identity, so one-shot and persistent checkpoints cannot be
mixed silently. In one-shot mode the wall limit covers the whole Julia
process. In persistent mode it covers one request after the worker is ready,
including decoder construction on a cache miss.

The persistent protocol uses a private TOML file mailbox. Julia's stdout and
stderr go only to a worker log, so package chatter cannot be mistaken for a
response. Every request and response carries the session, sequence, source,
runtime, `(L,p)`, record, and cutoff identity. A per-call wall timeout kills
the entire worker process group and makes that record unknown; the next record
may start a fresh worker. The task runner also takes a nonblocking lock for the
stratum and closes the worker cleanly. `Sys.maxrss()` from this worker is
recorded as `rss_scope = "worker_lifetime_high_water"`: it is evidence for the
whole worker lifetime, not a per-call peak. The external supervisor must still
enforce the frozen 16 GiB cap, followed by the runner's post-hoc cap check.
Record authentication, request/response identity mismatches, broken request
pipes, and worker exits without an authenticated response are fatal provenance
errors. They abort without writing a scientific unknown checkpoint. A reported
memory limit retires the worker before any later record is attempted.

The preregistered L3 per-call limits for chi
`32,64,128,240,480,960` are `120,300,900,1800,3600,7200` seconds. Exact L3
calls have a 300-second limit. L5 retains its separately frozen table. Both the
one-shot and persistent executors read these same values from the immutable
execution manifest.

Persistence removes repeated Julia startup, backend/context construction, and
decoder construction. It does not reduce the boundary-MPS SVD contraction
cost, so it should be treated as orchestration overhead reduction rather than
a new decoding algorithm.

Run the bounded phase-1 tests from the repository root with:

```sh
PYTHONPYCACHEPREFIX=/tmp/tnml-pycache \
  python3 -m unittest \
  examples/tnml_repaired_seed/adaptive_chi_pilot/test_phase1.py -v
```

Run the full bounded suite, including the fake persistent-process protocol,
with:

```sh
PYTHONPYCACHEPREFIX=/tmp/tnml-pycache \
  python3 -m unittest \
  examples/tnml_repaired_seed/adaptive_chi_pilot/test_phase1.py \
  examples/tnml_repaired_seed/adaptive_chi_pilot/test_phase2.py \
  examples/tnml_repaired_seed/adaptive_chi_pilot/test_persistent_worker.py -v
```

The one-record fallback requires `adaptive_runner.py --execute`. The persistent
task runner requires `task_runner.py --execute-task`. `--task-id` preserves
the 20-history L3 shard and one-history L5 schedule; `--record-index` addresses
the 30 independent BLAS=1 replay records. The older `--p-key` grouping remains
available for local experiments but is not the HPC schedule because it would
discard most array-level parallelism.

The isolated, non-submitting HPC4 control plane and exact commands are in
`hpc/README.md`. It creates separate L3 primary, L3 BLAS=1 replay, and L5 roots
without modifying the production HPC wrappers. L5 is blocked by a write-once
gate until all 600 exact audits and 30 numerical replays pass.
The L3 gate calibrates the conditional error of histories that the rule
accepts; the preregistration did not set a minimum accepted fraction, so the
report always exposes the accepted denominator and does not claim production
utility. The L5 feasibility gate separately requires at least 80% early
acceptance among non-anchors. Its median and p90 totals include the one-time
Julia process/backend/context startup charged to each one-record L5 task as
well as every TN request.
