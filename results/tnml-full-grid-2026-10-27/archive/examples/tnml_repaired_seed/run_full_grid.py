"""Resume a frozen matched-preparation TN scan, one supervised point at a time.

The zero phase performs one deterministic implementation check per size. These
checks are not Monte Carlo estimates. Production requests 1,000 histories at
each positive rate, and requires an explicit passed validation gate per size.
The existing Julia runner owns scientific results and strict resume identities.
STOP is observed between points; it never interrupts an active point.
"""
import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import tomllib


DIRECTORY = Path(__file__).resolve().parent
REPOSITORY = DIRECTORY.parents[1]
SIZES = [3, 5, 7, 9, 11]
POSITIVE_RATES = [.0025, .0045, .005, .0055, .006, .0065, .007, .0075,
                  .01, .0125, .015, .0175, .02, .0225, .025, .0275,
                  .03, .0325, .035, .0375, .04, .0425]
CUTOFF = 1e-18
NEGATIVE_TOLERANCE = 1e-10


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def document_hash(document):
    return hashlib.sha256(json.dumps(document, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def atomic_json(path, document):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", dir=path.parent, delete=False) as stream:
        temporary = Path(stream.name)
        try:
            json.dump(document, stream, indent=2, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    temporary.replace(path)


def scientific_source_hashes():
    paths = list((REPOSITORY / "src").rglob("*.jl"))
    paths += [DIRECTORY / name for name in ("run.jl", "model.jl", "dense_frontier.jl")]
    paths += [REPOSITORY / name for name in ("Project.toml", "Manifest.toml") if (REPOSITORY / name).is_file()]
    return {str(path.relative_to(REPOSITORY)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(paths)}


def selection_template():
    return {
        "schema_version": 1,
        "purpose": "User-requested full-range finite-bond TN scan; no threshold fit",
        "sizes": SIZES,
        "p_values": [0.0] + POSITIVE_RATES,
        "shots_per_nonzero_point": 1000,
        "zero_checks_per_size": 1,
        "zero_interpretation": "The analytical failure probability is zero. One deterministic TN implementation check per size is recorded separately, not as 1,000 sampled histories.",
        "chi_by_size": {str(L): 240 for L in SIZES},
        "cutoff": CUTOFF,
        "negative_tolerance": NEGATIVE_TOLERANCE,
        "exact_reference_sizes": [3],
        "zero_rng_seeds": {str(L): 202610280000 + L for L in SIZES},
        "production_rng_seeds": {str(L): {f"{p:.4f}": 202610290000 + 1000 * L + index
                                         for index, p in enumerate(POSITIVE_RATES, 1)} for L in SIZES},
        "seed_policy": "Unique fixed Random.Xoshiro seeds for every size/rate, separate from the old pilot and diagnostics.",
        "resource_limits": {
            "zero": {"per_contraction_seconds": 300, "per_run_total_seconds": 600, "worker_rss_gib": 4},
            "production": {"per_contraction_seconds": 600, "per_run_total_seconds": 86400, "worker_rss_gib": 4},
            "one_numerical_worker": True, "blas_threads": 1,
        },
        "validation_gate": "validation/gates.json must explicitly pass each production size at the frozen contraction settings. L3 also requires an independent exact reference for every production history.",
        "halt_policy": "A timeout, worker failure, unknown or invalid result, missing exact reference, or non-tied exact ML disagreement blocks further points of that size. Resume never retries a blocked point unless --retry-blocked is explicitly supplied.",
        "unknown_policy": "Retain every requested history; missing or invalid outcomes are unknown, never successes or dropped samples. Never substitute exact output for a TN result.",
        "physical_preparation": "Independent noisy seed-H outputs, complete initial seed MWPM repair, diagonal H+6CNOT transfer with independent Z gate-output faults; no later feedback.",
        "decoder_information": "TN: seed and final syndromes plus known seed recovery. Paired final MWPM: final syndrome only.",
        "threshold_fit": False,
        "scientific_source_sha256": scientific_source_hashes(),
    }


def initialize_selection(root):
    """Create once or verify the frozen selection; never rewrite existing data."""
    root = Path(root)
    expected = selection_template()
    path = root / "scan-selection.json"
    if path.exists():
        selection = json.loads(path.read_text())
        if {k: v for k, v in selection.items() if k != "frozen_at_utc"} != expected:
            raise ValueError("frozen selection or scientific sources differ; use a new study, not overwritten checkpoints")
        return selection
    selection = {"frozen_at_utc": timestamp(), **expected}
    atomic_json(path, selection)
    return selection


def selected_points(selection, phase, sizes=None, p_values=None):
    if phase not in ("zero", "production"):
        raise ValueError("phase must be zero or production")
    sizes = list(selection["sizes"] if sizes is None else sizes)
    allowed = [0.0] if phase == "zero" else [p for p in selection["p_values"] if p > 0]
    rates = list(allowed if p_values is None else p_values)
    if not sizes or len(set(sizes)) != len(sizes) or any(L not in selection["sizes"] for L in sizes):
        raise ValueError("sizes must be a nonempty unique subset of the frozen grid")
    if not rates or len(set(rates)) != len(rates) or any(p not in allowed for p in rates):
        raise ValueError("p-values must be a nonempty unique subset of the phase's frozen probability grid")
    points = []
    for L in sizes:
        chi = selection["chi_by_size"][str(L)]
        for p in rates:
            name = f"L{L}-p{p:.4f}-chi{chi}"
            points.append(dict(
                key=f"{phase}:{name}", phase=phase, L=L, p=p, chi=chi,
                shots=selection["zero_checks_per_size"] if phase == "zero" else selection["shots_per_nonzero_point"],
                seed=selection["zero_rng_seeds"][str(L)] if phase == "zero" else selection["production_rng_seeds"][str(L)][f"{p:.4f}"],
                exact=L in selection["exact_reference_sizes"], output=f"{phase}/{name}"))
    return points


def require_gate(root, selection, point):
    if point["phase"] == "zero":
        return None
    path = root / "validation" / "gates.json"
    if not path.is_file():
        raise ValueError(f"L={point['L']} production has no validation gate")
    document = json.loads(path.read_text())
    if document.get("schema_version") != 1:
        raise ValueError("unsupported validation gate schema")
    gate = document.get("gates", {}).get(str(point["L"]), {})
    if gate.get("status") != "passed":
        raise ValueError(f"L={point['L']} validation is {gate.get('status', 'missing')}; production refused")
    for field, expected in (("chi", point["chi"]), ("cutoff", selection["cutoff"]),
                            ("negative_tolerance", selection["negative_tolerance"])):
        if gate.get(field) != expected:
            raise ValueError(f"L={point['L']} validation {field} differs from frozen settings")
    if gate.get("exact_reference_for_every_shot") is not point["exact"]:
        raise ValueError(f"L={point['L']} validation exact-reference policy differs")
    evidence = gate.get("evidence")
    if not isinstance(evidence, list) or not evidence or any(not isinstance(p, str) or not (root / p).is_file() for p in evidence):
        raise ValueError(f"L={point['L']} validation needs existing evidence files")
    return {"document_sha256": document_hash(document), "gate": gate,
            "evidence_sha256": {p: hashlib.sha256((root / p).read_bytes()).hexdigest() for p in evidence}}


def build_command(root, selection, point, julia):
    limits = selection["resource_limits"][point["phase"]]
    output = root / point["output"]
    return [sys.executable, str(DIRECTORY / "supervise.py"), "--output", str(output),
            "--decode-seconds", str(limits["per_contraction_seconds"]),
            "--total-seconds", str(limits["per_run_total_seconds"]),
            "--rss-gib", str(limits["worker_rss_gib"]), "--", julia,
            "--startup-file=no", "--threads=1,0", "--gcthreads=1", f"--project={REPOSITORY}",
            str(DIRECTORY / "run.jl"), str(point["L"]), str(point["p"]), str(point["shots"]),
            str(point["chi"]), str(point["seed"]), str(output), "exact=" + str(point["exact"]).lower()]


def read_toml(path):
    with path.open("rb") as stream:
        return tomllib.load(stream)


def validate_metadata_identity(metadata, point, selection=None):
    """Bind a saved run to its selected configuration and frozen source files.

    A caller without a selection gets configuration checks only. Production
    dispatch and the full-grid report always supply the frozen selection.
    """
    expected = {"L": point["L"], "p": point["p"], "shots": point["shots"],
                "chi": point["chi"], "rng_seed": point["seed"], "exact": point["exact"],
                "cutoff": CUTOFF, "negative_tolerance": NEGATIVE_TOLERANCE}
    if any(metadata.get(key) != value for key, value in expected.items()):
        raise ValueError("metadata configuration differs from selected point")
    if selection is None:
        return
    hashes = selection["scientific_source_sha256"]
    expected_sources = {name: digest for name, digest in hashes.items()
                        if name not in ("Project.toml", "Manifest.toml")}
    if metadata.get("source_sha256") != expected_sources:
        raise ValueError("saved scientific source hashes differ from frozen selection")
    if metadata.get("project_sha256") != hashes.get("Project.toml"):
        raise ValueError("saved project hash differs from frozen selection")
    if metadata.get("manifest_sha256") != hashes.get("Manifest.toml"):
        raise ValueError("saved manifest hash differs from frozen selection")


def audit_point(root, point, selection=None):
    """Audit saved results independently of the Julia worker's exit status."""
    output = root / point["output"]
    audit = dict(passed=False, valid_shots=0, unknown_shots=point["shots"],
                 exact_checked=0, exact_disagreements=0, strict_matches=0,
                 max_posterior_linf=0.0, scientific_identity_checked=False, errors=[])
    try:
        metadata = read_toml(output / "metadata.toml")
        summary = read_toml(output / "summary.toml")
        corpus = read_toml(output / "input-corpus.toml")
        validate_metadata_identity(metadata, point, selection)
        audit["scientific_identity_checked"] = selection is not None
        run_id = metadata["run_id"]
        records = corpus["records"]
        if not run_id or corpus.get("run_id") != run_id or summary.get("run_id") != run_id:
            raise ValueError("run identities differ")
        if len(records) != point["shots"] or summary.get("requested_shots") != point["shots"]:
            raise ValueError("requested history count differs")
        if summary.get("L") != point["L"] or summary.get("p") != point["p"]:
            raise ValueError("summary selected point differs")
        failures = 0
        for i, record in enumerate(records, 1):
            if record.get("shot_index") != i:
                raise ValueError("saved input order differs")
            try:
                shot = read_toml(output / "shots" / f"shot_{i:04d}.toml")
                if any(shot.get(k) != v for k, v in (("run_id", run_id), ("record_id", record["record_id"]),
                                                    ("shot_index", i), ("L", point["L"]), ("p", point["p"]), ("chi", point["chi"]))):
                    raise ValueError("output identity differs")
                if shot.get("status") != "complete" or shot.get("valid") is not True:
                    raise ValueError("TN result is invalid or incomplete")
                probabilities = shot["probabilities"]
                if len(probabilities) != 4 or any(not math.isfinite(v) or v < -NEGATIVE_TOLERANCE for v in probabilities) or abs(sum(probabilities) - 1) > 1e-7:
                    raise ValueError("TN posterior is not a valid normalized four-sector distribution")
                label = shot["logical_mask"]
                if label not in range(4) or shot.get("TNfailure") is not (label != shot["true_mask"]):
                    raise ValueError("TN logical failure scoring differs")
                audit["valid_shots"] += 1
                failures += shot["TNfailure"]
                if point["phase"] == "zero" and (record.get("raw_fault_indices") != [] or label != 0 or
                    shot["TNfailure"] or max(abs(value - (index == 0)) for index, value in enumerate(probabilities)) > 1e-8):
                    raise ValueError("zero-noise implementation check differs from analytical zero")
                if point["exact"]:
                    exact = read_toml(output / "exact" / f"shot_{i:04d}.toml")
                    if exact.get("status") != "complete" or exact.get("valid") is not True or any(
                        exact.get(k) != shot[k] for k in ("run_id", "record_id", "shot_index", "L", "p", "chi")):
                        raise ValueError("independent exact reference is invalid or has a different identity")
                    reference = exact["probabilities"]
                    if len(reference) != 4 or any(not math.isfinite(v) or v < -NEGATIVE_TOLERANCE for v in reference) or abs(sum(reference) - 1) > 1e-7:
                        raise ValueError("exact posterior is invalid")
                    audit["exact_checked"] += 1
                    audit["max_posterior_linf"] = max(audit["max_posterior_linf"], max(abs(a - b) for a, b in zip(probabilities, reference)))
                    audit["strict_matches"] += label == exact["logical_mask"]
                    comparison = shot.get("exact_comparison", {})
                    if max(reference) - reference[label] > 1e-9 or comparison.get("ml_tie_aware_match") is not True:
                        audit["exact_disagreements"] += 1
                        raise ValueError("TN decision disagrees with independent exact ML beyond tie tolerance")
            except (OSError, ValueError, KeyError, TypeError) as error:
                audit["errors"].append(f"shot {i}: {error}")
        audit["unknown_shots"] = point["shots"] - audit["valid_shots"]
        audit["tn_failures"] = failures
        if summary.get("valid_shots") != audit["valid_shots"] or summary.get("tn_failures") != failures or summary.get("unknown_shots") != audit["unknown_shots"]:
            audit["errors"].append("summary counts differ from saved outputs")
        if audit["valid_shots"] != point["shots"]:
            audit["errors"].append("some requested TN outcomes remain unknown")
        if point["exact"] and audit["exact_checked"] != point["shots"]:
            audit["errors"].append("some requested independent exact references are missing")
    except (OSError, ValueError, KeyError, TypeError) as error:
        audit["errors"].append(str(error))
    audit["passed"] = not audit["errors"]
    return audit


@contextmanager
def scan_lock(root):
    """Prevent two full-grid drivers from dispatching workers in the same study."""
    root.mkdir(parents=True, exist_ok=True)
    with (root / ".scan.lock").open("a+") as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise ValueError("another scan driver holds this study's worker lock") from error
        try:
            yield
        finally:
            fcntl.flock(stream, fcntl.LOCK_UN)


def execute(root, selection, phase, sizes, p_values, julia, *, retry_blocked=False, run=subprocess.run):
    root = Path(root).resolve()
    points = selected_points(selection, phase, sizes, p_values)
    # Refuse an unvalidated requested size before launching any other size.
    for point in points:
        require_gate(root, selection, point)
    with scan_lock(root):
        destination = root / "scan-status.json"
        status = json.loads(destination.read_text()) if destination.exists() else dict(
            schema_version=1, selection_sha256=document_hash(selection), points={}, size_blocks={})
        if status.get("selection_sha256") != document_hash(selection):
            raise ValueError("scan-status identity differs from frozen selection")
        for point in points:
            status["points"].setdefault(point["key"], dict(
                phase=phase, L=point["L"], p=point["p"], requested_shots=point["shots"],
                output=point["output"], status="pending", reason="", attempts=[]))

        def checkpoint(state):
            status["status"] = state
            status["updated_at_utc"] = timestamp()
            atomic_json(destination, status)

        checkpoint("running")
        failed = False
        retried = set()
        for point in points:
            if (root / "STOP").exists():
                checkpoint("stopped")
                return 3
            size_key = f"{phase}:L{point['L']}"
            block = status["size_blocks"].get(size_key)
            if block:
                if not retry_blocked or size_key in retried or block["point"] != point["key"]:
                    failed = True
                    continue
                del status["size_blocks"][size_key]
                retried.add(size_key)
            row = status["points"][point["key"]]
            # A scan can run for days while the checkout is edited. Recheck
            # before every new Julia process, not merely at driver startup.
            if scientific_source_hashes() != selection["scientific_source_sha256"]:
                row.update(status="blocked", reason="frozen scientific sources changed during scan; no further worker launched")
                status["size_blocks"][size_key] = {"point": point["key"], "reason": row["reason"]}
                checkpoint("source_changed")
                return 2
            output = root / point["output"]
            supervision_file = output / "supervision.json"
            previous = json.loads(supervision_file.read_text()) if supervision_file.exists() else None
            if previous and previous.get("status") != "completed" and not retry_blocked:
                row.update(status="blocked", reason="previous supervised run: " + previous.get("status", "unknown"))
                status["size_blocks"][size_key] = {"point": point["key"], "reason": row["reason"]}
                failed = True
                checkpoint("running")
                continue
            gate = require_gate(root, selection, point)
            command = build_command(root, selection, point, julia)
            row.update(status="running", reason="")
            attempt = dict(started_at_utc=timestamp(), command=command, validation=gate)
            if previous:
                attempt["previous_supervision"] = previous
            row["attempts"].append(attempt)
            checkpoint("running")
            print(f"BEGIN {phase} L={point['L']} p={point['p']} N={point['shots']}", flush=True)
            try:
                completed = run(command, cwd=REPOSITORY)
                supervision = json.loads(supervision_file.read_text()) if supervision_file.exists() else {}
                attempt.update(returncode=completed.returncode, supervision=supervision)
                if completed.returncode != 0 or supervision.get("status") != "completed":
                    reason = supervision.get("status", "supervisor_failed")
                    if reason == "completed":
                        reason = "supervisor_failed"
                    row.update(status="blocked", reason=reason)
                else:
                    audit = audit_point(root, point, selection)
                    attempt["audit"] = audit
                    row["completed_counts"] = {key: audit[key] for key in (
                        "valid_shots", "unknown_shots", "exact_checked", "exact_disagreements", "strict_matches", "max_posterior_linf")}
                    row.update(status="completed" if audit["passed"] else "blocked",
                               reason="" if audit["passed"] else "; ".join(audit["errors"][:5]))
            except (OSError, ValueError) as error:
                row.update(status="blocked", reason=f"driver or supervisor error: {error}")
            except KeyboardInterrupt:
                row.update(status="interrupted", reason="user interruption; saved histories retained")
                attempt["finished_at_utc"] = timestamp()
                checkpoint("interrupted")
                raise
            attempt["finished_at_utc"] = timestamp()
            if row["status"] != "completed":
                status["size_blocks"][size_key] = {"point": point["key"], "reason": row["reason"]}
                failed = True
            checkpoint("running")
            print(f"END {phase} L={point['L']} p={point['p']} status={row['status']}", flush=True)
        checkpoint("blocked" if failed else "completed_selected_points")
        return 2 if failed else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--julia", default="julia")
    parser.add_argument("--phase", choices=("zero", "production"))
    parser.add_argument("--sizes", type=int, nargs="+")
    parser.add_argument("--p-values", type=float, nargs="+", help="Probabilities, e.g. 0.006 for 0.6%%; frozen-grid subsets only")
    parser.add_argument("--initialize-only", action="store_true")
    parser.add_argument("--retry-blocked", action="store_true", help="Explicitly allow one new attempt of a blocked selected point; no automatic retry loops")
    args = parser.parse_args()
    if not args.initialize_only and args.phase is None:
        parser.error("provide --phase zero|production or --initialize-only")
    try:
        root = args.root.resolve()
        selection = initialize_selection(root)
        if args.initialize_only:
            print(root / "scan-selection.json")
            return 0
        return execute(root, selection, args.phase, args.sizes, args.p_values, args.julia,
                       retry_blocked=args.retry_blocked)
    except (OSError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    raise SystemExit(main())
