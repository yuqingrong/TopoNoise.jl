#!/bin/bash -l
# Submit from the transferred repository root.
#SBATCH --job-name=tnml-collect
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-tnml-collect-%j.out
#SBATCH --error=slurm-tnml-collect-%j.err

set -euo pipefail
: "${SLURM_JOB_ID:?Submit this file using sbatch}"
: "${SLURM_SUBMIT_DIR:?Missing submission directory}"
[[ $# -eq 1 ]] || { echo "Usage: sbatch examples/tnml_repaired_seed/hpc/collect.sbatch STUDY_DIRECTORY" >&2; exit 2; }
cd -- "$SLURM_SUBMIT_DIR"
module load julia/1.12.4
source examples/tnml_repaired_seed/hpc/env.sh
export MPLCONFIGDIR="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/tnml-matplotlib-${SLURM_JOB_ID}"
export MPLBACKEND=Agg
exec srun --ntasks=1 --cpus-per-task=1 \
    "$TNML_PYTHON" "$TNML_REPO/examples/tnml_repaired_seed/hpc/collect.py" --root "$1" --plot
