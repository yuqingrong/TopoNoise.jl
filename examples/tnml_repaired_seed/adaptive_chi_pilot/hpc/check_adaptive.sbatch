#!/bin/bash -l
# Six standalone histories, one per p. Existing frozen studies are untouched.
#SBATCH --job-name=tnml-ladder-check
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=0-5%6
#SBATCH --output=slurm-tnml-ladder-check-%A_%a.out
#SBATCH --error=slurm-tnml-ladder-check-%A_%a.err

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?Submit from the TopoNoise repository root}"
: "${SLURM_ARRAY_TASK_ID:?Submit as an array with task IDs 0 through 5}"
[[ $# -le 1 ]] || { echo "Usage: sbatch --partition=... check_adaptive.sbatch [OUTPUT_ROOT]" >&2; exit 2; }
[[ "$SLURM_ARRAY_TASK_ID" =~ ^[0-5]$ ]] || { echo "Task ID must be 0 through 5" >&2; exit 2; }
cd -- "$SLURM_SUBMIT_DIR"
module load julia/1.12.4
export TNML_BLAS_THREADS=4
source examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/env.sh
tnml_diagnostic_root="${1:-$HOME/tnml-adaptive-ladder-checks}"
mkdir -p -- "$tnml_diagnostic_root"
tnml_diagnostic_output=$(mktemp -d "$tnml_diagnostic_root/job-${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}-XXXXXX")
printf 'Standalone diagnostic output: %s\n' "$tnml_diagnostic_output"
srun --ntasks=1 --cpus-per-task=4 --cpu-bind=cores \
  julia --project="$TNML_ADAPTIVE_REPO" --startup-file=no --threads=1 \
  examples/tnml_repaired_seed/adaptive_chi_pilot/check_adaptive_ladder.jl \
  "$SLURM_ARRAY_TASK_ID" "$tnml_diagnostic_output" \
  2>&1 | tee "$tnml_diagnostic_output/worker.log"
"$TNML_ADAPTIVE_PYTHON" \
  examples/tnml_repaired_seed/adaptive_chi_pilot/inspect_adaptive_ladder.py \
  "$tnml_diagnostic_output"
