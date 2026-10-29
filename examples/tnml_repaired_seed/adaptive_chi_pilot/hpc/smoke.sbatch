#!/bin/bash -l
#SBATCH --job-name=tnml-adapt-smoke
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=slurm-tnml-adapt-smoke-%j.out
#SBATCH --error=slurm-tnml-adapt-smoke-%j.err

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?Submit from the TopoNoise repository root}"
[[ $# -eq 2 ]] || { echo "Usage: sbatch --partition=... smoke.sbatch STUDY_ROOT TARGET" >&2; exit 2; }
cd -- "$SLURM_SUBMIT_DIR"
module load julia/1.12.4
case "$2" in
  l3-replay) export TNML_BLAS_THREADS=1 ;;
  l3-primary|l5-primary) export TNML_BLAS_THREADS=4 ;;
esac
source examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/env.sh
exec srun --ntasks=1 --cpus-per-task=4 --cpu-bind=cores \
  "$TNML_ADAPTIVE_PYTHON" examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/control.py smoke \
  --study-root "$1" --target "$2" --probe \
  --repository-root "$TNML_ADAPTIVE_REPO" --julia "$(command -v julia)"
