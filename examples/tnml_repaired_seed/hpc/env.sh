#!/bin/bash
# Source this file after loading julia/1.12.4. Do not execute it directly.
# The repository must stay at the same absolute path after setup.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo "Source this file after loading julia/1.12.4." >&2
    exit 2
fi

TNML_REPO="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd -P)" || return 2
export TNML_REPO
[[ -f "$TNML_REPO/Project.toml" && -f "$TNML_REPO/CondaPkg.toml" ]] || {
    echo "Cannot locate the transferred TopoNoise repository." >&2; return 2;
}

# Isolate this Linux installation from personal Julia/Python environments.
# HOME is intentionally left unchanged.
export JULIA_PROJECT="$TNML_REPO"
export JULIA_LOAD_PATH='@:@stdlib'
export JULIA_DEPOT_PATH="$TNML_REPO/.hpc-julia-depot"
export JULIA_CPU_TARGET=generic
export JULIA_NUM_THREADS=1
export JULIA_NUM_PRECOMPILE_TASKS=1
export JULIA_PKG_PRECOMPILE_AUTO=0
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONNOUSERSITE=1
unset PYTHONHOME PYTHONPATH VIRTUAL_ENV CONDA_PREFIX
unset JULIA_PYTHONCALL_EXE JULIA_PYTHONCALL_LIB JULIA_CONDAPKG_ENV JULIA_CONDAPKG_EXE
unset JULIA_CONDAPKG_BACKEND JULIA_CONDAPKG_OFFLINE JULIA_CONDAPKG_PIP_BACKEND
unset __JULIA_PYTHONCALL_EMBEDDED_LIBPTR__
unset TNML_PYTHON

if [[ "${TNML_HPC_SETUP:-0}" == 1 ]]; then
    # A real pip installation also supplies the immutable pip-freeze record.
    export JULIA_CONDAPKG_PIP_BACKEND=pip
    export JULIA_PKG_OFFLINE=false
else
    [[ ! -d "$TNML_REPO/.hpc-setup.lock" ]] || {
        echo "HPC setup is still active; wait for it to finish." >&2; return 2;
    }
    [[ -f "$TNML_REPO/.hpc-python-path" ]] || {
        echo "Native Python is not validated. Run setup.sbatch successfully first." >&2; return 2;
    }
    TNML_PYTHON="$(cat -- "$TNML_REPO/.hpc-python-path")" || return 2
    case "$TNML_PYTHON" in
        "$TNML_REPO"/*) ;;
        *) echo "The saved Python belongs to another repository path; use a fresh HPC copy." >&2; return 2 ;;
    esac
    [[ "$TNML_PYTHON" != *$'\n'* && -x "$TNML_PYTHON" ]] || {
        echo "The saved native Python path is invalid." >&2; return 2;
    }
    export TNML_PYTHON
    export JULIA_PYTHONCALL_EXE="$TNML_PYTHON"
    # Array tasks must not resolve or modify the shared Python environment.
    export JULIA_CONDAPKG_BACKEND=Null
    export JULIA_CONDAPKG_OFFLINE=true
    export JULIA_PKG_OFFLINE=true
fi
