#!/bin/bash -l
# Read-only chi=240/480 convergence diagnostic for one saved L5 history.
# Submit from the TopoNoise repository root.
#SBATCH --job-name=tnml-chi
#SBATCH --partition=a128m512u
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=05:00:00
#SBATCH --output=slurm-tnml-chi-%j.out
#SBATCH --error=slurm-tnml-chi-%j.err

set -uo pipefail
: "${SLURM_JOB_ID:?Submit this file with sbatch}"
: "${SLURM_SUBMIT_DIR:?Missing Slurm submission directory}"
if [[ $# -ne 2 ]]; then
    echo "Usage: sbatch examples/tnml_repaired_seed/chi_probe.sbatch FAILED_L5_TASK_DIR NEW_OUTPUT_DIR" >&2
    exit 2
fi

cd -- "$SLURM_SUBMIT_DIR" || exit 2
[[ -f Project.toml && -f examples/tnml_repaired_seed/chi_probe.jl ]] || {
    echo "Submit from the TopoNoise repository root." >&2
    exit 2
}
module load julia/1.12.4 || exit 2
source examples/tnml_repaired_seed/hpc/env.sh || exit 2

source_run="$1"
probe_output="$2"
julia_executable="$(command -v julia)" || exit 2
julia_command=("$julia_executable" --project="$TNML_REPO" --startup-file=no --threads=1,0 --gcthreads=1
    "$TNML_REPO/examples/tnml_repaired_seed/chi_probe.jl")

export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=4
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

srun --ntasks=1 --cpus-per-task=8 --cpu-bind=cores \
    "${julia_command[@]}" --validate "$source_run" "$probe_output" || exit 2

worker_failure=0
for chi in 240 480; do
    echo "Starting fresh chi=$chi BLAS=4 worker at $(date -Is)"
    if timeout --signal=TERM --kill-after=30s 2h \
        srun --exclusive --ntasks=1 --cpus-per-task=8 --cpu-bind=cores \
        "${julia_command[@]}" --worker "$source_run" "$probe_output" "$chi"; then
        echo "Finished chi=$chi BLAS=4 worker at $(date -Is)"
    else
        worker_exit=$?
        worker_status=worker_failed
        [[ "$worker_exit" -eq 124 ]] && worker_status=contraction_time_limit
        echo "chi=$chi status=$worker_status exit_code=$worker_exit" >&2
        srun --ntasks=1 --cpus-per-task=8 --cpu-bind=cores \
            "${julia_command[@]}" --failure "$source_run" "$probe_output" \
            "$chi" "$worker_status" "$worker_exit" 7200 || true
        worker_failure=1
    fi
done

srun --ntasks=1 --cpus-per-task=8 --cpu-bind=cores \
    "${julia_command[@]}" --collect "$source_run" "$probe_output"
collector_status=$?
(( worker_failure == 0 && collector_status == 0 )) || exit 2
