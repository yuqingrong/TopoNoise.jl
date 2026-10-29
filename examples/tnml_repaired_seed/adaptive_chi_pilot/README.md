# Adaptive-chi pilot: phase 1

This directory currently contains only the tested control plane for the pilot
described in `../ADAPTIVE_CHI_PILOT.md`. It cannot launch Julia, run a tensor
network contraction, or submit a Slurm job.

`manifest.py` creates an immutable sampling-plan manifest. It freezes the six
probability strata, every history seed, task layout, chi ladders, stopping
policy, resource limits, and the exact hash-ranking rule. `plan_task.py` only
reads that manifest and prints a dry-run task plan. `core.py` implements the
pure stopping, ranking, checkpoint-validation, and unknown-accounting rules.

Record IDs do not exist until the raw history corpus has been generated. A
later phase must therefore write a **new execution manifest**, rather than
mutating this sampling manifest. That execution manifest must bind the frozen
corpus and selected record IDs/ranks along with source and package hashes,
BLAS vendor, compute-node identity, and other runtime identity required by the
design document. Numerical execution remains disabled until that second freeze
is implemented and reviewed.

That phase-2 execution manifest must also freeze the independent-exact-ML tie
tolerance at `1e-9`, the posterior and log-joint audit gates at `1e-5`, and the
one-sided 95% exact Clopper-Pearson confidence-bound convention.

Run the bounded phase-1 tests from the repository root with:

```sh
PYTHONPYCACHEPREFIX=/tmp/tnml-pycache \
  python3 -m unittest \
  examples/tnml_repaired_seed/adaptive_chi_pilot/test_phase1.py -v
```
