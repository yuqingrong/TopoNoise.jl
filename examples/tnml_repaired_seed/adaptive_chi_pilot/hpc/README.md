# Adaptive-chi HPC4 workflow

This directory is isolated from `examples/tnml_repaired_seed/hpc`. It uses the
already validated `.hpc-python-path` and Julia depot, but never edits the older
fixed-chi study.

Before restarting a batch after a worker fix, run the one-history diagnostic:

```bash
cd "$HOME/TopoNoise.jl"
sbatch --partition=a128m512u \
  examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/check_worker.sbatch
```

It runs a real L3 history at p=0.006 with chi 32, 128, and 960, compares with
exact ML, and checks decoder cache reuse. It reports invalid low-chi outputs
and requires chi 960 to pass. Unlike the startup smoke, it actually contracts
the TN. Its output is `slurm-tnml-worker-check-JOBID.out`; existing frozen studies
are untouched. Passing this one-history check does not validate the adaptive
stopping policy. A changed worker requires a new study root for subsequent
frozen batches; retain old exception checkpoints rather than deleting them or
reusing their source identity.

After that loading check passes, exercise the adaptive policy with six histories
(one at each of p=0.006, 0.01, 0.02, 0.03, 0.035, and 0.0425):

```bash
cd "$HOME/TopoNoise.jl"
sbatch --partition=a128m512u --array=0-5%6 \
  examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/check_adaptive.sbatch
```

Each task requests four CPUs and 16 GiB and runs one history at all six chi
levels (32, 64, 128, 240, 480, 960), plus exact ML. This is six histories total,
not six histories per chi. The full ladder tests where the adaptive rule would
stop and checks that decision against exact ML. It is a diagnostic, not the
production calibration gate. The first contraction includes compilation, and
the reported prefix time is reconstructed from the full-ladder diagnostic,
not a separately measured adaptive runtime.

Logs are `slurm-tnml-ladder-check-JOBID_TASK.out` and `.err`. Each task prints
its new output directory under `$HOME/tnml-adaptive-ladder-checks`; that directory
contains the sampled record, exact reference, all chi results, and `summary.json`.
Read `adaptive_decision`, `accepted_exact_audit`, and `highest_chi_exact_audit`
in the summary. Exit zero requires the highest-chi exact audit to pass and no
false adaptive acceptance; the adaptive decision can still be `unknown`.
No old study files are changed. Retain the same code until these tasks finish.

Stopping policy v2 retries explicit numerical invalidity at the next chi and
requires a new window of three consecutive valid levels. Exceptions, malformed
output, and resource failures remain terminal unknowns. The policy version is
frozen in new manifests; subsequent production work needs a fresh study root.

The frozen schedule is:

- L3 primary: 30 array elements, 20 histories each, at most 20 running;
- L3 numerical replay: 30 hash-selected records with BLAS=1, at most 20 running;
- L5 feasibility: two separately submitted 30-element arrays, one history per
  element, at most 20 running.

Every science task requests four CPUs and 16 GiB. Scripts do not name a
partition; pass it to `sbatch`. L3 caches decoders only within one manifest
shard, preserving array parallelism. L5 has one record per task.
Each compute task reauthenticates the frozen source, Python package lock,
runtime profile, retained corpus shards, and corpus-provenance document before
starting Julia. A restarted Julia worker repeats its source/runtime check.
Commands marked `--wait` keep that shell attached until the job or array ends;
run the sequence inside `tmux` so disconnecting SSH does not interrupt it.
If a waited science array returns nonzero, stop the dependent sequence. For an
L3 science/replay failure, submit `collect.sbatch "$TNML_STUDY" l3` directly;
it reports missing work and writes no gate. For an L5 science failure, submit
the analogous `l5` collector. Fix/resubmit the failed array before continuing.

Run from the repository root after pulling the branch:

```bash
cd "$HOME/TopoNoise.jl"
export TNML_STUDY="$HOME/tnml-adaptive-20261030"
export TNML_SEED_BASE=202610600000
export TNML_PARTITION="${TNML_PARTITION:-a128m512u}"
tnml_submit() {
  local tnml_submit_output
  tnml_submit_output=$(sbatch --parsable "$@") || return
  printf '%s\n' "${tnml_submit_output%%;*}"
}

python3 examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/control.py setup \
  --study-root "$TNML_STUDY" --seed-base "$TNML_SEED_BASE" \
  --repository-root "$PWD" --dry-run

tnml_setup_job=$(tnml_submit --partition="$TNML_PARTITION" \
  examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/setup.sbatch \
  "$TNML_STUDY" "$TNML_SEED_BASE")
```

After setup completes, generate all L3 corpus shards, freeze separate BLAS=4
primary and BLAS=1 replay manifests, and run both read-only smokes:

```bash
l3_corpus_job=$(tnml_submit --partition="$TNML_PARTITION" \
  --dependency="afterok:$tnml_setup_job" --array=0-29%20 \
  examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/corpus.sbatch \
  "$TNML_STUDY" l3_exact)

l3_freeze_job=$(tnml_submit --partition="$TNML_PARTITION" \
  --dependency="afterok:$l3_corpus_job" \
  examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/freeze.sbatch \
  "$TNML_STUDY" l3_exact)

l3_smoke_primary_job=$(tnml_submit --partition="$TNML_PARTITION" \
  --dependency="afterok:$l3_freeze_job" \
  examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/smoke.sbatch \
  "$TNML_STUDY" l3-primary)
l3_smoke_replay_job=$(tnml_submit --partition="$TNML_PARTITION" \
  --dependency="afterok:$l3_freeze_job" \
  examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/smoke.sbatch \
  "$TNML_STUDY" l3-replay)
```

These four submissions remain below the 50-submitted-job limit. Wait until
both smoke jobs are `COMPLETED` before adding the 30-element science array:

```bash
sacct -j "$l3_smoke_primary_job,$l3_smoke_replay_job" -X \
  --format=JobID,State,Elapsed,ExitCode
```

Submit primary and replay separately to stay below the normal-QOS limit of 50
submitted jobs:

```bash
l3_primary_job=$(tnml_submit --partition="$TNML_PARTITION" \
  --dependency="afterok:$l3_smoke_primary_job:$l3_smoke_replay_job" \
  --wait --array=0-29%20 \
  examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/l3_calibration.sbatch \
  "$TNML_STUDY")

# --wait above ensures every primary array element finished before this submission.
l3_replay_job=$(tnml_submit --partition="$TNML_PARTITION" \
  --dependency="afterok:$l3_primary_job" --array=0-29%20 \
  examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/l3_replay.sbatch \
  "$TNML_STUDY")

# afterany collects failures as missing/unknown evidence instead of racing them.
l3_collect_job=$(tnml_submit --partition="$TNML_PARTITION" \
  --dependency="afterany:$l3_replay_job" \
  examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/collect.sbatch \
  "$TNML_STUDY" l3)
```

The collector writes a hashed `l3-gate.json`. It passes only when all 600
primary decisions and exact references exist, every accepted result passes the
frozen exact gates, all 30 BLAS=1 decisions match, and their per-chi posterior
and log-joint changes are at most `1e-6`. It authorizes only the L5 feasibility
smoke.

L5 commands fail before work without that passed gate. Submit each 30-element
batch separately:

```bash
l5_corpus_a_job=$(tnml_submit --partition="$TNML_PARTITION" --wait --array=0-29%20 \
  examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/corpus.sbatch \
  "$TNML_STUDY" l5_feasibility)
# Wait, then generate the second half.
l5_corpus_b_job=$(tnml_submit --partition="$TNML_PARTITION" \
  --dependency="afterok:$l5_corpus_a_job" --wait --array=30-59%20 \
  examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/corpus.sbatch \
  "$TNML_STUDY" l5_feasibility)

l5_freeze_job=$(tnml_submit --partition="$TNML_PARTITION" \
  --dependency="afterok:$l5_corpus_b_job" --wait \
  examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/freeze.sbatch \
  "$TNML_STUDY" l5_feasibility)
l5_smoke_job=$(tnml_submit --partition="$TNML_PARTITION" \
  --dependency="afterok:$l5_freeze_job" --wait \
  examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/smoke.sbatch \
  "$TNML_STUDY" l5-primary)

# Wait for the smoke to complete before adding another 30 array elements.
l5_science_a_job=$(tnml_submit --partition="$TNML_PARTITION" \
  --dependency="afterok:$l5_smoke_job" --wait --array=0-29%20 \
  examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/l5_feasibility.sbatch \
  "$TNML_STUDY")
# Wait, then submit the second half.
l5_science_b_job=$(tnml_submit --partition="$TNML_PARTITION" \
  --dependency="afterok:$l5_science_a_job" --array=30-59%20 \
  examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/l5_feasibility.sbatch \
  "$TNML_STUDY")

l5_collect_job=$(tnml_submit --partition="$TNML_PARTITION" \
  --dependency="afterany:$l5_science_b_job" \
  examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/collect.sbatch \
  "$TNML_STUDY" l5)
```

Read-only progress can be printed at any time:

```bash
module load julia/1.12.4
export TNML_BLAS_THREADS=4
source examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/env.sh
"$TNML_ADAPTIVE_PYTHON" \
  examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/report.py \
  --study-root "$TNML_STUDY" --progress
```

Scheduler timeout or preemption is resumable by resubmitting the same array
ID. Completed immutable checkpoints are reused. A per-call timeout, reported
memory limit, exception, or malformed result is a terminal unknown for that
manifest. An explicit numerical-invalid TN result retries the next chi; if
the final level is invalid, the record remains unknown. L5 is
a feasibility smoke, not threshold data. Its gate also requires all 24 forced
anchors to finish through chi 480 with an explicit unchanged logical class.
No L7/L9/L11 scripts are provided.
