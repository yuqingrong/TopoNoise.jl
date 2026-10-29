#!/bin/bash -l
#SBATCH --job-name=tnml-adapt-corpus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-tnml-adapt-corpus-%A_%a.out
#SBATCH --error=slurm-tnml-adapt-corpus-%A_%a.err

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?Submit from the TopoNoise repository root}"
: "${SLURM_ARRAY_TASK_ID:?Submit with the documented array range}"
[[ $# -eq 2 ]] || { echo "Usage: sbatch --partition=... --array=... corpus.sbatch STUDY_ROOT STAGE" >&2; exit 2; }
case "$2" in l3_exact|l5_feasibility) ;; *) echo "Unknown corpus stage" >&2; exit 2 ;; esac
cd -- "$SLURM_SUBMIT_DIR"
module load julia/1.12.4
export TNML_BLAS_THREADS=4
source examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/env.sh
exec srun --ntasks=1 --cpus-per-task=4 --cpu-bind=cores \
  "$TNML_ADAPTIVE_PYTHON" examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/control.py corpus-task \
  --study-root "$1" --stage "$2" --task-id "$SLURM_ARRAY_TASK_ID" \
  --repository-root "$TNML_ADAPTIVE_REPO" --julia "$(command -v julia)"
