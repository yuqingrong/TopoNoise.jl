#!/bin/bash -l
# Submit from the repository root; Slurm copies this script to its spool.
#SBATCH --job-name=tnml
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=08:00:00
#SBATCH --output=slurm-tnml-%A_%a.out
#SBATCH --error=slurm-tnml-%A_%a.err

set -euo pipefail
: "${SLURM_JOB_ID:?Submit this file using sbatch, from the repository root}"
: "${SLURM_SUBMIT_DIR:?Missing submission directory}"
: "${SLURM_ARRAY_TASK_ID:?Submit with --array; see hpc/README.md}"
if [[ $# -ne 2 ]]; then
    echo "Usage: sbatch [resources] --array=... examples/tnml_repaired_seed/hpc/array.sbatch STAGE STUDY_DIRECTORY" >&2
    exit 2
fi
case "$1" in zero|smoke|feasibility|production) ;; *) echo "Unknown stage: $1" >&2; exit 2 ;; esac
if [[ "${SLURM_CPUS_PER_TASK:-1}" != 1 ]]; then
    echo "This runner uses one Julia/BLAS thread. Request one CPU per task; use array concurrency for parallelism." >&2
    exit 2
fi
cd -- "$SLURM_SUBMIT_DIR"
[[ -f Project.toml && -f examples/tnml_repaired_seed/run.jl ]] || {
    echo "Submit from the transferred TopoNoise repository root." >&2; exit 2;
}
module load julia/1.12.4
source examples/tnml_repaired_seed/hpc/env.sh
exec srun --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
    "$TNML_PYTHON" "$TNML_REPO/examples/tnml_repaired_seed/hpc/run_task.py" \
    --root "$2" --stage "$1" --task-id "$SLURM_ARRAY_TASK_ID" --julia "$(command -v julia)"
