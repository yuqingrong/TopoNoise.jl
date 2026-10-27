#!/usr/bin/env python3
"""Run one immutable HPC task inside an existing Slurm allocation.

There is no study-wide mutable status or lock. Each task owns one output
directory. Re-submission resumes the identical saved histories; a scientific
halt requires investigation and a new validated study, never an automatic retry.
"""
import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import importlib.util
import json
import os
from pathlib import Path
import socket
import subprocess
import sys

if sys.version_info < (3, 11):
    raise SystemExit("Use the bootstrapped Python >= 3.11, not the system Python 3.9.")

_SPEC = importlib.util.spec_from_file_location("hpc_task_manifest", Path(__file__).with_name("manifest.py"))
HPC = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(HPC)


def timestamp():
    return datetime.now(timezone.utc).isoformat()


@contextmanager
def task_lock(output):
    output.mkdir(parents=True, exist_ok=True)
    with (output / ".task.lock").open("a+") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise ValueError("another worker holds this task's output lock") from error
        try:
            yield lock.fileno()
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)


def smoke_evidence(root, document):
    evidence = []
    for task in HPC.stage_tasks(document, "smoke"):
        output = HPC.output_directory(root, task)
        supervision_file = output / "supervision.json"
        supervision = json.loads(supervision_file.read_text()) if supervision_file.is_file() else {}
        audit = HPC.audit_task(root, task, document)
        if supervision.get("status") != "completed" or supervision.get("exit_code") != 0 or not audit["passed"]:
            raise ValueError(f"production blocked: smoke task {task['task_id']} has not passed its complete independent exact audit")
        evidence.append(dict(task_id=task["task_id"], output=task["output"], audit=audit,
                             supervision_sha256=HPC.hashlib.sha256(supervision_file.read_bytes()).hexdigest()))
    return evidence


def build_command(root, document, task, julia):
    output = HPC.output_directory(root, task)
    limits = document["stages"][task["stage"]]["resource_limits"]
    return [sys.executable, str(HPC.EXPERIMENT / "supervise.py"), "--output", str(output),
            "--decode-seconds", str(limits["per_contraction_seconds"]),
            "--total-seconds", str(limits["per_run_total_seconds"]),
            "--rss-gib", str(limits["worker_rss_gib"]), "--", julia,
            "--startup-file=no", "--threads=1,0", "--gcthreads=1", f"--project={HPC.REPOSITORY}",
            str(HPC.EXPERIMENT / "run.jl"), str(task["L"]), str(task["p"]), str(task["shots"]),
            str(task["chi"]), str(task["seed"]), str(output), "exact=" + str(task["exact"]).lower(),
            "cutoff=" + str(task["cutoff"])]


def execute_task(root, stage, task_id, julia="julia", *, dry_run=False, environ=None, run=subprocess.run):
    root = Path(root).resolve()
    environment = dict(os.environ if environ is None else environ)
    if not dry_run:
        if not environment.get("SLURM_JOB_ID"):
            raise ValueError("SLURM_JOB_ID is required; use --dry-run for a read-only preview outside Slurm")
        array_id = environment.get("SLURM_ARRAY_TASK_ID")
        if array_id is not None and array_id != str(task_id):
            raise ValueError("selected task ID differs from the Slurm array task ID")
    document = HPC.load_manifest(root)
    task = HPC.task_for_id(document, stage, task_id)
    output = HPC.output_directory(root, task)
    gate = smoke_evidence(root, document) if stage == "production" else []
    if (output / "halt.toml").exists():
        raise ValueError("saved scientific halt prevents automatic retry; preserve this study and investigate the failed history")
    command = build_command(root, document, task, julia)
    if dry_run:
        return dict(status="dry_run", manifest_id=document["manifest_id"], task=task, command=command,
                    resource_limits=document["stages"][stage]["resource_limits"], smoke_evidence=gate)
    environment.update(JULIA_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    environment.setdefault("JULIA_CONDAPKG_BACKEND", "Null")
    environment.setdefault("JULIA_PYTHONCALL_EXE", str(Path(sys.executable).resolve()))
    with task_lock(output) as lock_fd:
        # Recheck after acquiring ownership, before touching a previous attempt.
        HPC.load_manifest(root)
        if (output / "halt.toml").exists():
            raise ValueError("saved scientific halt prevents automatic retry")
        status_file = output / "task-status.json"
        task_hash = HPC.DRIVER.document_hash(task)
        if status_file.exists():
            status = json.loads(status_file.read_text())
            if status.get("manifest_id") != document["manifest_id"] or status.get("task_sha256") != task_hash:
                raise ValueError("saved task status has a different manifest or task identity")
        else:
            status = dict(manifest_id=document["manifest_id"], task_sha256=task_hash, task=task, attempts=[])
        metadata_file = output / "metadata.toml"
        if metadata_file.exists():
            HPC.validate_task_metadata(HPC.DRIVER.read_toml(metadata_file), task, document)
        supervision_file = output / "supervision.json"
        previous = json.loads(supervision_file.read_text()) if supervision_file.exists() else None
        attempt = dict(started_at_utc=timestamp(), slurm_job_id=environment["SLURM_JOB_ID"],
                       slurm_array_task_id=environment.get("SLURM_ARRAY_TASK_ID"),
                       host=socket.gethostname(), command=command, smoke_evidence=gate)
        if previous:
            attempt["previous_supervision"] = previous
        status["attempts"].append(attempt)
        status.update(status="running", reason="", updated_at_utc=timestamp())
        HPC.DRIVER.atomic_json(status_file, status)
        try:
            completed = run(command, cwd=HPC.REPOSITORY, env=environment, pass_fds=(lock_fd,))
            supervision = json.loads(supervision_file.read_text()) if supervision_file.exists() else {}
            attempt.update(returncode=completed.returncode, supervision=supervision)
            audit = HPC.audit_task(root, task, document)
            attempt["audit"] = audit
            if completed.returncode != 0 or supervision.get("status") != "completed" or supervision.get("exit_code") != 0:
                reason = supervision.get("status", "supervisor_failed")
                if reason == "completed":
                    reason = "supervisor_failed"
                status.update(status="blocked", reason=reason)
            elif not audit["passed"]:
                status.update(status="blocked", reason="; ".join(audit["errors"][:5]))
            else:
                status.update(status="completed", reason="")
        except (OSError, ValueError) as error:
            status.update(status="blocked", reason=str(error))
        except KeyboardInterrupt:
            status.update(status="interrupted", reason="interrupted; saved histories retained")
            raise
        finally:
            attempt["finished_at_utc"] = timestamp()
            status["updated_at_utc"] = timestamp()
            HPC.DRIVER.atomic_json(status_file, status)
        return status


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--stage", choices=HPC.STAGES, required=True)
    parser.add_argument("--task-id", type=int, required=True)
    parser.add_argument("--julia", default="julia")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    try:
        result = execute_task(args.root, args.stage, args.task_id, args.julia, dry_run=args.dry_run)
    except (OSError, ValueError, KeyError) as error:
        parser.error(str(error))
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0 if result["status"] in ("completed", "dry_run") else 2


if __name__ == "__main__":
    raise SystemExit(main())
