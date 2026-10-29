#!/bin/bash -l
# Frozen collection wrapper; its hash is bound by each execution manifest.
#SBATCH --job-name=tnml-adapt-collect
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-tnml-adapt-collect-%j.out
#SBATCH --error=slurm-tnml-adapt-collect-%j.err

set -euo pipefail
: "${SLURM_SUBMIT_DIR:?Submit from the TopoNoise repository root}"
[[ $# -eq 2 ]] || { echo "Usage: sbatch --partition=... collect.sbatch STUDY_ROOT l3|l5" >&2; exit 2; }
case "$2" in l3) tnml_collect_flag=--collect-l3 ;; l5) tnml_collect_flag=--collect-l5 ;; *) echo "Mode must be l3 or l5" >&2; exit 2 ;; esac
cd -- "$SLURM_SUBMIT_DIR"
module load julia/1.12.4
export TNML_BLAS_THREADS=4
source examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/env.sh
exec srun --ntasks=1 --cpus-per-task=1 \
  "$TNML_ADAPTIVE_PYTHON" examples/tnml_repaired_seed/adaptive_chi_pilot/hpc/report.py \
  --study-root "$1" "$tnml_collect_flag"
