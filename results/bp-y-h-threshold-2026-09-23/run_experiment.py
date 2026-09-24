"""Run/resume the approved pilot, ambiguity refinement, and fresh production scan."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import platform
from pathlib import Path
import sys
import time

import numpy as np
import pymatching
import stim

from analysis import CHANNELS, DISTANCES, counts_and_rate, load_points, select_crossing, write_csv
from sampler import atomic_json, load_circuits, run_point

RUN = Path(__file__).resolve().parent
ROOT_SEED = 20260923
BATCH_SIZE = 1000


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def provenance():
    versions = dict(python=sys.version, numpy=np.__version__, pymatching=pymatching.__version__,
                    stim=stim.__version__, machine=platform.machine(), platform=platform.platform())
    hashes = {name: file_hash(RUN / name) for name in
              ("circuits.toml", "source_hashes.toml", "sampler.py")}
    payload = dict(versions=versions, source_hashes=hashes, root_seed=ROOT_SEED,
                   batch_size=BATCH_SIZE, count_order=["00", "10", "01", "11"])
    fingerprint = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
    return dict(**payload, fingerprint=fingerprint, construction="bp", boundary_orientation="x_ns",
                clock="gate_layer", initialization="ideal Bell preimage; ideal reference qubit",
                final_measurements="ideal A_s, B_p, Z_L Z_R, X_L X_R",
                decoder="independent uniform-weight CSS PyMatching; enable_correlations=False",
                channels=dict(y_only="(1-p)rho + p Y rho Y; shared X/Z event",
                              hadamard="(1-p)rho + p H rho H; H=(X+Z)/sqrt(2), exact trajectories"),
                distances=list(DISTANCES), pilot_shots=2000, refinement_shots=10000,
                production_shots=50000, sampling_complete=False)


def _run_task(arguments):
    return run_point(**arguments)


def run_tasks(specifications, circuits, metadata, workers):
    tasks = []
    saved = {(r["stage"], r["channel"], r["distance"], int(round(r["p"]*1e6))): r
             for r in load_points()}
    for stage, channel, distance, p, shots in specifications:
        old = saved.get((stage, channel, distance, int(round(p*1e6))))
        target = max(shots, old["shots"] if old else 0)
        tasks.append(dict(circuit=circuits[distance], channel=channel, p=p, stage=stage,
                          shots=target, root_seed=ROOT_SEED, batch_size=BATCH_SIZE,
                          fingerprint=metadata["fingerprint"], directory=RUN / "checkpoints"))
    started = time.perf_counter()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(_run_task, task) for task in tasks]
        for i, future in enumerate(as_completed(futures), 1):
            row = future.result()
            x, z = (counts_and_rate(row, c)[1] for c in ("x", "z"))
            print(f"{i}/{len(tasks)} {row['stage']} {row['channel']} d={row['distance']} "
                  f"p={row['p']:.3f} N={row['shots']} PX={x:.5f} PZ={z:.5f} "
                  f"elapsed={time.perf_counter()-started:.1f}s", flush=True)


def select_panels():
    rows = [r for r in load_points() if r["stage"] == "pilot"]
    return {f"{channel}/{component}": select_crossing(rows, channel, component)
            for channel in CHANNELS for component in ("x", "z")}


def freeze_selection(panels):
    # Refinement has a finite, declared budget. An ambiguous result stays unresolved.
    for choice in panels.values():
        if choice["status"] == "ambiguous":
            choice["status"] = "unresolved"
            choice["reason"] += "; still ambiguous after the 10,000-shot refinement"
    pilot = [{k: r[k] for k in ("channel", "distance", "p", "shots", "counts")}
             for r in load_points() if r["stage"] == "pilot"]
    selection = dict(frozen=True, panels=panels,
                     pilot_counts_sha256=hashlib.sha256(json.dumps(pilot, sort_keys=True).encode()).hexdigest(),
                     center_grid=0.002, production_halfwidth=0.010, primary_halfwidth=0.006,
                     policy="Frozen from independent pilot before production; see analysis.select_crossing")
    path = RUN / "selection.json"
    if path.exists() and json.loads(path.read_text()) != selection:
        raise ValueError("Refusing to change a previously frozen production-window selection")
    atomic_json(path, selection)
    return selection


def production_specs(selection):
    specs = []
    for channel in CHANNELS:
        centers = [v["center"] for k, v in selection["panels"].items()
                   if k.startswith(channel + "/") and v["status"] == "selected"]
        rates = sorted(set(round(c + i * 0.002, 6) for c in centers for i in range(-5, 6)
                           if 0 <= c + i * 0.002 <= 1))
        specs += [("production", channel, d, p, 50000) for d in DISTANCES for p in rates]
    return specs


def export_raw():
    rows = []
    for point in load_points():
        n = point["shots"]
        counts = point["counts"]
        if sum(counts) != n or min(counts) < 0:
            raise ValueError("Invalid joint counts")
        row = {k: point[k] for k in ("stage", "channel", "distance", "p", "shots", "root_seed")}
        row.update(zip(("n00", "n10", "n01", "n11"), counts))
        row.update(construction="bp", boundary_orientation="x_ns", clock="gate_layer")
        for component in ("x", "z"):
            k, rate, _ = counts_and_rate(point, component)
            z = 1.959963984540054
            denominator = 1 + z*z/n
            center = (rate + z*z/(2*n))/denominator
            half = z*np.sqrt(rate*(1-rate)/n + z*z/(4*n*n))/denominator
            row.update({f"logical_{component}_failures": k,
                        f"logical_{component}_rate": rate,
                        f"logical_{component}_standard_error": np.sqrt(rate*(1-rate)/n),
                        f"logical_{component}_wilson_low": max(0.0, center-half),
                        f"logical_{component}_wilson_high": min(1.0, center+half)})
        row["any_logical_failures"] = n-counts[0]
        rows.append(row)
    write_csv(RUN / "raw_counts.csv", rows)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("pilot", "refine", "production", "all", "export"), default="all")
    parser.add_argument("--workers", type=int, default=4)
    arguments = parser.parse_args()
    if arguments.workers <= 0:
        parser.error("--workers must be positive")
    circuits = load_circuits(RUN / "circuits.toml")
    metadata = provenance()
    if (RUN / "metadata.json").exists():
        previous = json.loads((RUN / "metadata.json").read_text())
        if previous["fingerprint"] != metadata["fingerprint"]:
            raise ValueError("The sampler, environment, or archived circuit changed; use a new run directory")
        metadata = previous
    atomic_json(RUN / "metadata.json", metadata)
    if arguments.stage == "export":
        export_raw()
        return
    frozen = (RUN / "selection.json").exists()
    if arguments.stage in ("pilot", "all") and not frozen:
        specs = [("pilot", channel, d, round(i*0.005, 6), 2000)
                 for channel in CHANNELS for d in DISTANCES for i in range(21)]
        run_tasks(specs, circuits, metadata, arguments.workers)
        atomic_json(RUN / "selection-before-refinement.json", dict(panels=select_panels()))
    if arguments.stage in ("refine", "all") and not frozen:
        panels = select_panels()
        specs = set()
        for channel in CHANNELS:
            for component in ("x", "z"):
                choice = panels[f"{channel}/{component}"]
                if choice["status"] == "ambiguous":
                    for d in DISTANCES:
                        for p in choice["refine_rates"]:
                            specs.add(("pilot", channel, d, p, 10000))
        if specs:
            run_tasks(sorted(specs), circuits, metadata, arguments.workers)
        selection = freeze_selection(select_panels())
        print("Frozen windows:", json.dumps(selection["panels"], sort_keys=True), flush=True)
    if arguments.stage in ("production", "all"):
        selection = json.loads((RUN / "selection.json").read_text())
        run_tasks(production_specs(selection), circuits, metadata, arguments.workers)
        metadata["sampling_complete"] = True
    rows = export_raw()
    metadata.update(saved_points=len(rows), saved_shots=sum(r["shots"] for r in rows))
    atomic_json(RUN / "metadata.json", metadata)
    print(f"Saved {len(rows)} points, {metadata['saved_shots']:,} shots.", flush=True)


if __name__ == "__main__":
    main()
