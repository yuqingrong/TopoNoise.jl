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

Run the bounded phase-1 tests from the repository root with:

```sh
PYTHONPYCACHEPREFIX=/tmp/tnml-pycache \
  python3 -m unittest \
  examples/tnml_repaired_seed/adaptive_chi_pilot/test_phase1.py -v
```

Run both bounded phases with:

```sh
PYTHONPYCACHEPREFIX=/tmp/tnml-pycache \
  python3 -m unittest \
  examples/tnml_repaired_seed/adaptive_chi_pilot/test_phase1.py \
  examples/tnml_repaired_seed/adaptive_chi_pilot/test_phase2.py -v
```

Use `--dry-run` to plan without execution. An explicit `--execute` flag is
required before the command can start the Julia adapter.
