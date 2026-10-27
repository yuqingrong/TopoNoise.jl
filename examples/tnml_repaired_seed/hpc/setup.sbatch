#!/bin/bash -l
# Submit from the repository root: Slurm executes a spooled script copy.
# Resource requests are provisional; adjust to your allocation if necessary.
#SBATCH --job-name=tnml-setup
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=slurm-tnml-setup-%j.out
#SBATCH --error=slurm-tnml-setup-%j.err

set -euo pipefail
: "${SLURM_JOB_ID:?Submit with sbatch from the transferred repository root}"
: "${SLURM_SUBMIT_DIR:?Missing submission directory}"
cd -- "$SLURM_SUBMIT_DIR"
[[ -f Project.toml && -f CondaPkg.toml && -f examples/tnml_repaired_seed/hpc/env.sh ]] || {
    echo "Submit setup.sbatch from the transferred TopoNoise repository root." >&2; exit 2;
}
[[ "$(uname -s)" == Linux && "$(uname -m)" == aarch64 ]] || {
    echo "This bootstrap targets the confirmed Linux aarch64 HPC4 environment." >&2; exit 2;
}
[[ "${SLURM_CPUS_PER_TASK:-1}" == 1 ]] || {
    echo "Setup is configured for one CPU and one precompile task." >&2; exit 2;
}
[[ ! -e .hpc-python-path ]] || {
    echo "This copy already has a validated runtime. Use a fresh copy to rebuild; do not modify an active study." >&2
    exit 2
}
mkdir .hpc-setup.lock || {
    echo "Another setup owns .hpc-setup.lock. Check that job before removing a stale lock." >&2; exit 2;
}
trap 'rmdir -- "$SLURM_SUBMIT_DIR/.hpc-setup.lock"' EXIT
module load julia/1.12.4
export TNML_HPC_SETUP=1
source examples/tnml_repaired_seed/hpc/env.sh

tnml_log_dir="$TNML_REPO/hpc_environment/setup-${SLURM_JOB_ID}"
mkdir -p -- "$tnml_log_dir"
{
    date -u '+%Y-%m-%dT%H:%M:%SZ'
    uname -a
    module list 2>&1
    command -v julia
    julia --startup-file=no --version
    for tnml_tool in c++ g++ gcc cc cmake make; do
        if command -v "$tnml_tool" >/dev/null 2>&1; then
            printf '%s: %s\n' "$tnml_tool" "$(command -v "$tnml_tool")"
        else
            printf '%s: unavailable (may be needed if an ARM wheel is absent)\n' "$tnml_tool"
        fi
    done
} > "$tnml_log_dir/platform.txt"
cat -- "$tnml_log_dir/platform.txt"

# Only this setup job may instantiate/resolve. The tracked dependency files
# are not intentionally changed; create the frozen study after setup succeeds.
julia --startup-file=no --threads=1 --project="$TNML_REPO" -e '
    VERSION == v"1.12.4" || error("Load the confirmed julia/1.12.4 module")
    import Pkg
    Pkg.instantiate(; allow_autoprecomp=false)
    import CondaPkg
    CondaPkg.backend() in (:Pixi, :MicroMamba) || error("Use a managed native CondaPkg backend; remove conflicting local preferences")
    CondaPkg.resolve()
    Pkg.precompile()
' 2>&1 | tee "$tnml_log_dir/bootstrap.log"

julia --startup-file=no --threads=1 --project="$TNML_REPO" \
    "$TNML_REPO/examples/tnml_repaired_seed/hpc/environment_check.jl" \
    "$tnml_log_dir" 2>&1 | tee "$tnml_log_dir/check.log"
echo "Native runtime validated. Create the frozen study only now. Reports: $tnml_log_dir"
