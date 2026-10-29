#!/bin/bash -l
#SBATCH --job-name=tnml-adapt-replay
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=7-00:00:00
#SBATCH --output=slurm-tnml-adapt-replay-%A_%a.out
#SBATCH --error=slurm-tnml-adapt-replay-%A_%a.err

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?Submit from the TopoNoise repository root}"
: "${SLURM_ARRAY_TASK_ID:?Submit with --array=0-29%20}"
[[ $# -eq 1 ]] || { echo "Usage: sbatch --partition=... --array=0-29%20 l3_replay.sbatch STUDY_ROOT" >&2; exit 2; }
cd -- "$SLURM_SUBMIT_DIR"
module load julia/1.12.4
export TNML_BLAS_THREADS=1
source examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/env.sh
"$TNML_ADAPTIVE_PYTHON" examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/control.py \
  require-execution --study-root "$1" --target l3-replay \
  --repository-root "$TNML_ADAPTIVE_REPO" --julia "$(command -v julia)" >/dev/null
exec srun --ntasks=1 --cpus-per-task=4 --cpu-bind=cores \
  "$TNML_ADAPTIVE_PYTHON" examples/tnml_repaired_seed/adaptive_chi_pilot/task_runner.py \
  --study-root "$1/l3-blas1-replay" --sampling-manifest "$1/manifest.json" \
  --record-index "$SLURM_ARRAY_TASK_ID" --execute-task \
  --repository-root "$TNML_ADAPTIVE_REPO" --julia "$(command -v julia)"
