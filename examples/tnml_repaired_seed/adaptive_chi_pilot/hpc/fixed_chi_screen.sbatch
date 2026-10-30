#!/bin/bash -l
# Paired fixed-chi screening: 10 fresh histories per p, not a production gate.
#SBATCH --job-name=tnml-fixed
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --array=0-5%6
#SBATCH --output=slurm-tnml-fixed-%A_%a.out
#SBATCH --error=slurm-tnml-fixed-%A_%a.err

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?Submit from the TopoNoise repository root}"
: "${SLURM_ARRAY_TASK_ID:?Submit as an array with task IDs 0 through 5}"
[[ $# -le 2 ]] || { echo "Usage: sbatch fixed_chi_screen.sbatch [OUTPUT_ROOT [SHOTS]]" >&2; exit 2; }
[[ "$SLURM_ARRAY_TASK_ID" =~ ^[0-5]$ ]] || { echo "Task ID must be 0 through 5" >&2; exit 2; }
tnml_fixed_shots="${2:-10}"
[[ "$tnml_fixed_shots" =~ ^([1-9]|10)$ ]] || { echo "SHOTS must be 1 through 10" >&2; exit 2; }
cd -- "$SLURM_SUBMIT_DIR"
module load julia/1.12.4
export TNML_BLAS_THREADS=4
source examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/env.sh
tnml_fixed_root="${1:-$HOME/tnml-fixed-chi-screen}"
mkdir -p -- "$tnml_fixed_root"
tnml_fixed_output=$(mktemp -d "$tnml_fixed_root/job-${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}-XXXXXX")
printf 'Fixed-chi screening output: %s\n' "$tnml_fixed_output"
srun --ntasks=1 --cpus-per-task=4 --cpu-bind=cores \
  julia --project="$TNML_ADAPTIVE_REPO" --startup-file=no --threads=1 \
  examples/tnml_repaired_seed/adaptive_chi_pilot/fixed_chi_screen.jl \
  "$SLURM_ARRAY_TASK_ID" "$tnml_fixed_output" "$tnml_fixed_shots" \
  2>&1 | tee "$tnml_fixed_output/worker.log"
"$TNML_ADAPTIVE_PYTHON" \
  examples/tnml_repaired_seed/adaptive_chi_pilot/inspect_fixed_chi_screen.py \
  "$tnml_fixed_output"
