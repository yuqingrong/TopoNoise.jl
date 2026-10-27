"""Run one owned Julia worker with contraction, total-time and RSS limits."""
import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import time
import tomllib


def read_progress(directory):
    try:
        return tomllib.loads((directory / "progress.toml").read_text())
    except (FileNotFoundError, tomllib.TOMLDecodeError):
        return {}


def progress_marker(progress):
    # stage_started_unix is an opaque generation tag, never a deadline clock.
    return (progress.get("pid"), progress.get("stage"), progress.get("shot"),
            progress.get("stage_started_unix"))


class OwnedStageTimer:
    """Time a freshly observed TN stage of this supervisor's owned worker."""
    def __init__(self, worker_pid, initial_progress):
        self.worker_pid = worker_pid
        self.initial_marker = progress_marker(initial_progress)
        self.marker = None
        self.started = None

    def observe(self, progress, monotonic_now):
        marker = progress_marker(progress)
        # Ignore another process's checkpoint and an unchanged pre-launch file,
        # including the rare case where the OS has reused its worker PID.
        if progress.get("pid") == self.worker_pid and marker != self.initial_marker:
            if marker != self.marker:
                self.marker = marker
                self.started = monotonic_now if progress.get("stage") == "tn" else None
        return None if self.started is None else monotonic_now - self.started


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--decode-seconds", type=float, default=180)
    parser.add_argument("--total-seconds", type=float, default=900)
    parser.add_argument("--rss-gib", type=float, default=4)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command or min(args.decode_seconds, args.total_seconds, args.rss_gib) <= 0:
        parser.error("provide a command and positive limits")
    args.output.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    initial_progress = read_progress(args.output)
    peak_rss = 0
    rss_available = True
    reason = None
    last_stage = None
    progress = {}
    with (args.output / "worker.log").open("a") as log:
        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                   start_new_session=True)
        stage_timer = OwnedStageTimer(process.pid, initial_progress)
        stage_elapsed = None
        try:
            while process.poll() is None:
                progress = read_progress(args.output)
                marker = (progress.get("stage"), progress.get("shot"))
                if marker != last_stage:
                    print(f"stage={marker[0]} shot={marker[1]}", flush=True)
                    last_stage = marker
                monotonic_now = time.monotonic()
                elapsed = monotonic_now - started
                stage_elapsed = stage_timer.observe(progress, monotonic_now)
                if elapsed > args.total_seconds:
                    reason = "total_time_limit"
                if stage_elapsed is not None and stage_elapsed > args.decode_seconds:
                    reason = "contraction_time_limit"
                if rss_available:
                    try:
                        inspected = subprocess.run(["ps", "-o", "rss=", "-p", str(process.pid)],
                                                   capture_output=True, text=True)
                    except OSError:
                        inspected = None
                    if inspected is not None and inspected.returncode == 0 and inspected.stdout.strip():
                        rss = int(inspected.stdout.strip()) * 1024
                        peak_rss = max(peak_rss, rss)
                        if rss > args.rss_gib * 1024**3:
                            reason = "resident_memory_limit"
                    elif process.poll() is None:
                        rss_available = False
                        print("RSS inspection unavailable; time limits remain active", flush=True)
                if reason:
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                    break
                time.sleep(1)
        finally:
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
            exit_code = process.wait()
    result = {
        "status": reason or ("completed" if exit_code == 0 else "worker_failed"),
        "exit_code": exit_code,
        "elapsed_seconds": time.monotonic() - started,
        "timing_basis": "time.monotonic for total duration and first observation of a new owned-worker TN stage; wall-clock progress timestamps are generation tags only",
        "worker_pid": process.pid,
        "tn_stage_elapsed_seconds_at_last_poll": stage_elapsed,
        "peak_rss_bytes": peak_rss,
        "rss_limit_enforced": rss_available,
        "decode_limit_seconds": args.decode_seconds,
        "total_limit_seconds": args.total_seconds,
        "rss_limit_gib": args.rss_gib,
        "progress_at_stop": read_progress(args.output),
        "command": command,
    }
    destination = args.output / "supervision.json"
    temporary = destination.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(result, indent=2) + "\n")
    temporary.replace(destination)
    print(json.dumps(result), flush=True)
    return 0 if result["status"] == "completed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
