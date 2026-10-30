#!/bin/bash -l
# One-history diagnostic; does not read or modify a frozen study.
#SBATCH --job-name=tnml-worker-check
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=slurm-tnml-worker-check-%j.out
#SBATCH --error=slurm-tnml-worker-check-%j.err

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?Submit from the TopoNoise repository root}"
cd -- "$SLURM_SUBMIT_DIR"
module load julia/1.12.4
export TNML_BLAS_THREADS=4
source examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/env.sh
exec srun --ntasks=1 --cpus-per-task=4 --cpu-bind=cores \
  julia --project="$TNML_ADAPTIVE_REPO" --startup-file=no --threads=1 \
  examples/tnml_repaired_seed/adaptive_chi_pilot/test_real_late_loading.jl
