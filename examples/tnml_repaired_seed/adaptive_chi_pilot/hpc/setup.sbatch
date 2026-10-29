#!/bin/bash -l
#SBATCH --job-name=tnml-adapt-setup
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-tnml-adapt-setup-%j.out
#SBATCH --error=slurm-tnml-adapt-setup-%j.err

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?Submit from the TopoNoise repository root}"
[[ $# -eq 2 ]] || { echo "Usage: sbatch --partition=... setup.sbatch STUDY_ROOT SEED_BASE" >&2; exit 2; }
cd -- "$SLURM_SUBMIT_DIR"
module load julia/1.12.4
export TNML_BLAS_THREADS=4
source examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/env.sh
exec srun --ntasks=1 --cpus-per-task=4 --cpu-bind=cores \
  "$TNML_ADAPTIVE_PYTHON" examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/control.py setup \
  --study-root "$1" --seed-base "$2" --repository-root "$TNML_ADAPTIVE_REPO" \
  --julia "$(command -v julia)"
