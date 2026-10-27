"""Run all zero checks, then audited L3 points, with one worker at a time.

Requires an already frozen selection and validation evidence. Importing this
file or using --check does not initialize the study or launch numerical work.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parent
JULIA = "/Users/rongyuqing/.juliaup/bin/julia"


def plan(selection):
    if selection["sizes"] != [3, 5, 7, 9, 11] or selection["shots_per_nonzero_point"] != 1000:
        raise ValueError("This launcher requires the approved five-size, 1000-shot selection")
    if selection["cutoff"] != 1e-26 or selection["rng_seed_offset"] != 100000:
        raise ValueError("This launcher requires cutoff 1e-26 and seed offset 100000")
    if any(selection["chi_by_size"][str(L)] != 240 for L in selection["sizes"]):
        raise ValueError("This launcher requires chi=240 for every size")
    if selection["zero_checks_per_size"] != 1 or selection["exact_reference_sizes"] != [3]:
        raise ValueError("This launcher requires one zero check per size and exact L3 references")
    rates = selection["p_values"]
    if len(rates) != 23 or rates[0] != 0 or rates[-1] != .0425 or rates != sorted(set(rates)):
        raise ValueError("The approved 23-point probability grid is required")
    return ([dict(phase="zero", L=L, p=0.0, requested_shots=1) for L in selection["sizes"]]
            + [dict(phase="production", L=3, p=p, requested_shots=1000) for p in rates if p > 0])


def save_state(root, state):
    state["updated_at_utc"] = datetime.now(timezone.utc).isoformat()
    with tempfile.NamedTemporaryFile("w", dir=root, prefix="background-status-", suffix=".tmp", delete=False) as stream:
        temporary = Path(stream.name)
        json.dump(state, stream, indent=2)
        stream.write("\n")
    temporary.replace(root / "background-status.json")


def execute(root=ROOT, julia=JULIA, run=subprocess.run):
    root = Path(root).resolve()
    repository = root.parent.parent
    scripts = repository / "examples/tnml_repaired_seed"
    selection = json.loads((root / "scan-selection.json").read_text())
    jobs = plan(selection)
    environment = os.environ.copy()
    environment.update({
        "JULIA_DEPOT_PATH": "/tmp/toponoise-julia-depot:/Users/rongyuqing/.julia",
        "JULIA_CONDAPKG_BACKEND": "Null",
        "JULIA_PYTHONCALL_EXE": str(repository / ".CondaPkg/.pixi/envs/default/bin/python"),
        "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "PYTHONDONTWRITEBYTECODE": "1",
        "MPLCONFIGDIR": "/private/tmp/tnml-full-grid-precision-mpl",
    })
    with (root / ".background.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        state = dict(pid=os.getpid(), status="running", zero_check_sizes=selection["sizes"],
                     selected_production_sizes=[3], requested_sizes=selection["sizes"],
                     planned_production_histories=110000, cutoff=selection["cutoff"],
                     rng_seed_offset=selection["rng_seed_offset"], points=[],
                     started_at_utc=datetime.now(timezone.utc).isoformat())
        save_state(root, state)
        try:
            for job in jobs:
                state["current_point"] = dict(job)
                save_state(root, state)
                command = [sys.executable, str(scripts / "run_full_grid.py"), "--root", str(root),
                           "--phase", job["phase"], "--sizes", str(job["L"]), "--p-values", str(job["p"]),
                           "--cutoff", str(selection["cutoff"]), "--rng-seed-offset", str(selection["rng_seed_offset"]),
                           "--julia", julia]
                result = run(command, cwd=repository, env=environment)
                plot = run([sys.executable, str(scripts / "plot_full_grid.py"), "--root", str(root)],
                           cwd=repository, env=environment)
                state["points"].append({**job, "driver_returncode": result.returncode,
                                         "plot_returncode": plot.returncode})
                if result.returncode:
                    state["status"] = "stopped" if result.returncode == 3 else "blocked"
                    save_state(root, state)
                    return result.returncode
                if plot.returncode:
                    state["status"] = "report_error"
                    save_state(root, state)
                    return plot.returncode
                save_state(root, state)
        except KeyboardInterrupt:
            state["status"] = "interrupted"
            save_state(root, state)
            raise
        except Exception as error:
            state.update(status="launcher_error", error=str(error))
            save_state(root, state)
            raise
        state["status"] = "completed_zero_checks_and_selected_L3_grid"
        state["unrun_positive_sizes"] = [5, 7, 9, 11]
        save_state(root, state)
        return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Inspect the existing selection without launching or writing state")
    parser.add_argument("--julia", default=JULIA)
    args = parser.parse_args()
    if args.check:
        selection = json.loads((ROOT / "scan-selection.json").read_text())
        print(json.dumps(dict(status="not_launched_by_check", jobs=plan(selection)), indent=2))
        return 0
    return execute(julia=args.julia)


if __name__ == "__main__":
    raise SystemExit(main())
