"""Sequentially execute the frozen pilot selection, with resumable shot files."""
import argparse
import json
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--julia", default="julia")
    args = parser.parse_args()
    directory = Path(__file__).resolve().parent
    repository = directory.parents[1]
    root = args.root.resolve()
    selection = json.loads((root / "pilot-selection.json").read_text())
    limits = selection["resource_limits"]
    if selection["sizes"] != [3]:
        parser.error("this frozen pilot driver is restricted to the validated L=3 scan")
    rates, seeds = selection["p_values"], selection["rng_seeds"]
    if len(rates) != len(seeds):
        parser.error("each noise rate needs a separately declared RNG seed")
    outcomes = []
    for p, seed in zip(rates, seeds):
        output = root / "production" / f"L3-p{p:.4f}-chi{selection['chi']}"
        command = [sys.executable, str(directory / "supervise.py"),
                   "--output", str(output),
                   "--decode-seconds", str(limits["per_contraction_seconds"]),
                   "--total-seconds", str(limits["per_run_total_seconds"]),
                   "--rss-gib", str(limits["worker_rss_gib"]), "--",
                   args.julia, "--startup-file=no", "--threads=1,0", "--gcthreads=1",
                   f"--project={repository}", str(directory / "run.jl"), "3", str(p),
                   str(selection["shots_per_point"]), str(selection["chi"]), str(seed),
                   str(output), "exact=" + str(selection["exact_reference_for_every_shot"]).lower()]
        print(f"BEGIN L=3 p={p} N={selection['shots_per_point']}", flush=True)
        result = subprocess.run(command, cwd=repository)
        status_file = output / "supervision.json"
        status = json.loads(status_file.read_text()) if status_file.exists() else {}
        outcomes.append(dict(p=p, output=str(output), returncode=result.returncode,
                             status=status.get("status", "supervisor_failed")))
        temporary = root / "pilot-status.json.tmp"
        temporary.write_text(json.dumps(outcomes, indent=2) + "\n")
        temporary.replace(root / "pilot-status.json")
        print(f"END p={p} status={outcomes[-1]['status']}", flush=True)
        if outcomes[-1]["status"] in ("worker_failed", "supervisor_failed"):
            return 2
    return 0 if all(row["status"] == "completed" for row in outcomes) else 2


if __name__ == "__main__":
    raise SystemExit(main())
