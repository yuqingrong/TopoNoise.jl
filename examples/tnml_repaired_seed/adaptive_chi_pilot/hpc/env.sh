#!/bin/bash
# Source after loading julia/1.12.4. The adaptive pilot is isolated from the
# existing production HPC wrappers.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo "Source this file; do not execute it." >&2
    exit 2
fi

TNML_ADAPTIVE_REPO="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../../.." && pwd -P)" || return 2
export TNML_ADAPTIVE_REPO
[[ -f "$TNML_ADAPTIVE_REPO/Project.toml" && -f "$TNML_ADAPTIVE_REPO/.hpc-python-path" ]] || {
    echo "Validated HPC environment not found in this repository checkout." >&2
    return 2
}
TNML_ADAPTIVE_PYTHON="$(cat -- "$TNML_ADAPTIVE_REPO/.hpc-python-path")" || return 2
[[ -x "$TNML_ADAPTIVE_PYTHON" ]] || {
    echo "Saved HPC Python is missing or not executable." >&2
    return 2
}
export TNML_ADAPTIVE_PYTHON
export JULIA_PROJECT="$TNML_ADAPTIVE_REPO"
export JULIA_LOAD_PATH='@:@stdlib'
export JULIA_DEPOT_PATH="$TNML_ADAPTIVE_REPO/.hpc-julia-depot"
export JULIA_CPU_TARGET=generic
export JULIA_NUM_THREADS=1
export JULIA_NUM_PRECOMPILE_TASKS=1
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_PKG_OFFLINE=true
export JULIA_CONDAPKG_BACKEND=Null
export JULIA_CONDAPKG_OFFLINE=true
export JULIA_PYTHONCALL_EXE="$TNML_ADAPTIVE_PYTHON"
export PYTHONNOUSERSITE=1
unset PYTHONHOME PYTHONPATH VIRTUAL_ENV CONDA_PREFIX

TNML_BLAS_THREADS="${TNML_BLAS_THREADS:-4}"
case "$TNML_BLAS_THREADS" in 1|4) ;; *) echo "TNML_BLAS_THREADS must be 1 or 4." >&2; return 2 ;; esac
export TNML_BLAS_THREADS
export OPENBLAS_NUM_THREADS="$TNML_BLAS_THREADS"
export MKL_NUM_THREADS="$TNML_BLAS_THREADS"
export OMP_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS="$TNML_BLAS_THREADS"
export NUMEXPR_NUM_THREADS=1
