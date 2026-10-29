#!/usr/bin/env python3
"""Persistent Julia execution for one complete adaptive ``(L, p)`` stratum."""

from __future__ import annotations

import argparse
import fcntl
import json
import os
from pathlib import Path
import secrets
import signal
import subprocess
import tempfile
import time
from typing import Any, Callable, Mapping

try:
    from . import adaptive_runner, execution_manifest, manifest
except ImportError:
    import adaptive_runner
    import execution_manifest
    import manifest


PROTOCOL_VERSION = 1
WORKER_POLICY_NAME = "persistent_toml_mailbox_v1"


class PersistentJuliaExecutor(adaptive_runner.JuliaExecutor):
    """Serial file-mailbox client for a Julia process bound to one ``(L,p)``."""

    executor_mode = "persistent"

    def __init__(
        self,
        *,
        repository_root: os.PathLike[str] | str,
        julia: str = "julia",
        popen_factory: Callable[..., Any] = subprocess.Popen,
        monotonic: Callable[[], float] = time.monotonic,
        sleep: Callable[[float], None] = time.sleep,
        poll_interval: float = 0.05,
    ):
        super().__init__(repository_root=repository_root, julia=julia)
        self.worker_script = Path(__file__).resolve().with_name("persistent_worker.jl")
        self._popen_factory = popen_factory
        self._monotonic = monotonic
        self._sleep = sleep
        self._poll_interval = float(poll_interval)
        self._process = None
        self._session = None
        self._session_dir: Path | None = None
        self._session_token: str | None = None
        self._binding: tuple[Any, ...] | None = None
        self._sequence = 0
        self._stdout_log = None

    @staticmethod
    def _binding_for(identity: Mapping[str, Any]) -> tuple[Any, ...]:
        policy = identity.get("executor_policy")
        expected_policy = {
            "persistent_protocol": WORKER_POLICY_NAME,
            "protocol_version": PROTOCOL_VERSION,
            "scheduling": "serial_one_worker_per_L_p",
            "decoder_cache": "one_decoder_per_chi_fixed_L_p_cutoff",
            "allowed_executors": ["one_shot", "persistent"],
            "timing_scopes": {
                "one_shot": {
                    "wall_time_scope": "whole_julia_process",
                    "decoder_cache_scope": "none",
                },
                "persistent": {
                    "wall_time_scope": "request_after_worker_ready",
                    "decoder_cache_scope": "worker_lifetime_by_chi",
                },
            },
        }
        if not isinstance(policy, Mapping) or any(
            policy.get(key) != value for key, value in expected_policy.items()
        ):
            raise adaptive_runner.IdentityMismatchError(
                "execution identity does not match the frozen persistent protocol policy"
            )
        return (
            identity.get("source_bundle_id"),
            identity.get("runtime_profile_id"),
            identity.get("L"),
            identity.get("p_key"),
            identity.get("p"),
            identity.get("cutoff"),
        )

    def _worker_command(
        self, identity: Mapping[str, Any], ready_path: Path, token: str
    ) -> list[str]:
        return [
            self.julia,
            "--project=" + str(self.repository_root),
            "--startup-file=no",
            "--threads=1",
            str(self.worker_script),
            "--serve",
            "--L",
            str(identity["L"]),
            "--p",
            repr(identity["p"]),
            "--p-key",
            str(identity["p_key"]),
            "--cutoff",
            repr(identity["cutoff"]),
            "--ready",
            str(ready_path),
            "--session-token",
            token,
            "--source-bundle-id",
            str(identity["source_bundle_id"]),
            "--runtime-profile-id",
            str(identity["runtime_profile_id"]),
        ]

    def _log_tail(self) -> str:
        if self._session_dir is None:
            return ""
        path = self._session_dir / "worker.log"
        try:
            return path.read_text(encoding="utf-8", errors="replace")[-2000:]
        except OSError:
            return ""

    def _wait_for_file(self, path: Path, timeout: float | None) -> str:
        started = self._monotonic()
        deadline = None if timeout is None else started + timeout
        while True:
            now = self._monotonic()
            if deadline is not None and now >= deadline:
                return "timeout"
            if path.is_file():
                return "ready"
            if self._process is None or self._process.poll() is not None:
                return "exited"
            remaining = (
                self._poll_interval
                if timeout is None
                else min(self._poll_interval, max(0.0, deadline - now))
            )
            self._sleep(remaining)

    def _terminate_worker(self) -> None:
        process = self._process
        if process is None:
            return
        if process.poll() is None:
            try:
                pid = process.pid
            except AttributeError:
                pid = None
            try:
                if pid is None:
                    process.terminate()
                else:
                    os.killpg(pid, signal.SIGTERM)
            except (OSError, ProcessLookupError):
                try:
                    process.terminate()
                except OSError:
                    pass
            try:
                process.wait(timeout=2.0)
            except subprocess.TimeoutExpired:
                try:
                    if pid is None:
                        process.kill()
                    else:
                        os.killpg(pid, signal.SIGKILL)
                except (OSError, ProcessLookupError):
                    try:
                        process.kill()
                    except OSError:
                        pass
                try:
                    process.wait(timeout=2.0)
                except subprocess.TimeoutExpired:
                    pass
        self._process = None

    def _discard_session(self) -> None:
        if self._stdout_log is not None:
            try:
                self._stdout_log.close()
            except OSError:
                pass
        self._stdout_log = None
        if self._session is not None:
            self._session.cleanup()
        self._session = None
        self._session_dir = None
        self._session_token = None
        self._binding = None
        self._sequence = 0

    @staticmethod
    def _ready_expected(identity: Mapping[str, Any], token: str) -> dict[str, Any]:
        return {
            "protocol_version": PROTOCOL_VERSION,
            "session_token": token,
            "source_bundle_id": identity["source_bundle_id"],
            "runtime_profile_id": identity["runtime_profile_id"],
            "L": identity["L"],
            "p": identity["p"],
            "p_key": identity["p_key"],
            "cutoff": identity["cutoff"],
            "executor_mode": "persistent",
            "wall_time_scope": "request_after_worker_ready",
            "decoder_cache_scope": "worker_lifetime_by_chi",
        }

    def _start_worker(self, identity: Mapping[str, Any]) -> None:
        self._session = tempfile.TemporaryDirectory(prefix="tnml-persistent-worker-")
        self._session_dir = Path(self._session.name)
        self._session_token = secrets.token_hex(32)
        ready_path = self._session_dir / "ready.toml"
        self._stdout_log = (self._session_dir / "worker.log").open(
            "w", encoding="utf-8"
        )
        environment = os.environ.copy()
        environment["OPENBLAS_NUM_THREADS"] = "4"
        environment["MKL_NUM_THREADS"] = "4"
        command = self._worker_command(identity, ready_path, self._session_token)
        try:
            self._process = self._popen_factory(
                command,
                stdin=subprocess.PIPE,
                stdout=self._stdout_log,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                env=environment,
                start_new_session=True,
            )
            timeout = float(identity["resource_policy"]["preflight_wall_seconds"])
            observed = self._wait_for_file(ready_path, timeout)
            if observed != "ready":
                tail = self._log_tail()
                self._terminate_worker()
                raise adaptive_runner.IdentityMismatchError(
                    f"persistent Julia worker startup {observed}: {tail}"
                )
            ready = adaptive_runner._read_toml(ready_path)
            runtime = ready.pop("runtime_identity", None)
            self.verify_julia_runtime(runtime, identity["runtime_identity"])
            expected = self._ready_expected(identity, self._session_token)
            if any(ready.get(key) != value for key, value in expected.items()):
                raise adaptive_runner.IdentityMismatchError(
                    "persistent Julia worker ready identity mismatch"
                )
            if ready.get("status") != "ready" or ready.get("valid") is not True:
                raise RuntimeError(
                    "persistent Julia worker failed to initialize: "
                    + str(ready.get("error", self._log_tail()))
                )
            self._binding = self._binding_for(identity)
        except BaseException:
            self._terminate_worker()
            self._discard_session()
            raise

    def preflight(self, execution_identity: Mapping[str, Any]) -> None:
        self.verify_execution_environment(execution_identity)
        binding = self._binding_for(execution_identity)
        if self._process is not None and self._process.poll() is None:
            if self._binding != binding:
                raise adaptive_runner.IdentityMismatchError(
                    "live persistent worker is bound to a different (L,p,cutoff) identity"
                )
            return
        self._terminate_worker()
        self._discard_session()
        self._start_worker(execution_identity)

    @staticmethod
    def _request_envelope(
        identity: Mapping[str, Any],
        *,
        session_token: str,
        sequence: int,
        request_id: str,
        mode: str,
        record_path: Path,
        output_path: Path,
        chi: int,
        cutoff: float,
    ) -> dict[str, Any]:
        return {
            "protocol_version": PROTOCOL_VERSION,
            "session_token": session_token,
            "sequence": sequence,
            "request_id": request_id,
            "mode": mode,
            "record_path": str(record_path),
            "output_path": str(output_path),
            "source_bundle_id": identity["source_bundle_id"],
            "runtime_profile_id": identity["runtime_profile_id"],
            "L": identity["L"],
            "p": identity["p"],
            "p_key": identity["p_key"],
            "chi": chi,
            "cutoff": cutoff,
            "record_id": identity["record_id"],
            "julia_record_id": identity["julia_record_id"],
            "executor_mode": identity["executor_mode"],
            "wall_time_scope": identity["wall_time_scope"],
            "decoder_cache_scope": identity["decoder_cache_scope"],
        }

    def _validate_response_envelope(
        self, result: dict[str, Any], request: Mapping[str, Any]
    ) -> None:
        keys = (
            "protocol_version",
            "session_token",
            "sequence",
            "request_id",
            "mode",
            "source_bundle_id",
            "runtime_profile_id",
            "L",
            "p",
            "p_key",
            "chi",
            "cutoff",
            "record_id",
            "julia_record_id",
            "executor_mode",
            "wall_time_scope",
            "decoder_cache_scope",
        )
        if any(result.get(key) != request[key] for key in keys):
            self._terminate_worker()
            raise adaptive_runner.IdentityMismatchError(
                "persistent Julia worker response identity mismatch"
            )
        for key in keys:
            result.pop(key)

    def _request(
        self,
        mode: str,
        record: Mapping[str, Any],
        *,
        chi: int,
        cutoff: float,
        execution_identity: Mapping[str, Any],
    ) -> dict[str, Any]:
        self.preflight(execution_identity)
        assert self._session_dir is not None
        assert self._session_token is not None
        assert self._process is not None
        self._sequence += 1
        sequence = self._sequence
        request_id = secrets.token_hex(32)
        record_path = self._session_dir / f"record-{sequence:08d}.toml"
        request_path = self._session_dir / f"request-{sequence:08d}.toml"
        output_path = self._session_dir / f"result-{sequence:08d}.toml"
        adaptive_runner._write_record_toml(record_path, record)
        request = self._request_envelope(
            execution_identity,
            session_token=self._session_token,
            sequence=sequence,
            request_id=request_id,
            mode=mode,
            record_path=record_path,
            output_path=output_path,
            chi=chi,
            cutoff=cutoff,
        )
        adaptive_runner._write_record_toml(request_path, request)
        request_started = self._monotonic()
        try:
            self._process.stdin.write(str(request_path) + "\n")
            self._process.stdin.flush()
        except (AttributeError, BrokenPipeError, OSError) as error:
            tail = self._log_tail()
            self._terminate_worker()
            raise adaptive_runner.IdentityMismatchError(
                f"persistent worker request pipe failed before authentication: "
                f"{error}; {tail}"
            )
        timeout = self._wall_timeout(mode, chi, execution_identity)
        if timeout is None:
            self._terminate_worker()
            raise adaptive_runner.IdentityMismatchError(
                "persistent worker request has no frozen wall timeout"
            )
        observed = self._wait_for_file(output_path, timeout)
        if observed == "timeout":
            elapsed = max(0.0, self._monotonic() - request_started)
            self._terminate_worker()
            result = adaptive_runner._exception_result(
                "timeout", TimeoutError(f"persistent worker exceeded {timeout} seconds")
            )
            result["elapsed_seconds"] = elapsed
            return result
        if observed != "ready":
            tail = self._log_tail()
            self._terminate_worker()
            raise adaptive_runner.IdentityMismatchError(
                f"persistent worker exited without an authenticated result: {tail}"
            )
        try:
            result = adaptive_runner._read_toml(output_path)
            runtime = result.pop("runtime_identity", None)
        except Exception as error:
            self._terminate_worker()
            raise adaptive_runner.IdentityMismatchError(
                "persistent worker returned an unreadable result"
            ) from error
        try:
            self.verify_julia_runtime(runtime, execution_identity["runtime_identity"])
            self._validate_response_envelope(result, request)
            if result.get("rss_scope") != "worker_lifetime_high_water":
                raise adaptive_runner.IdentityMismatchError(
                    "persistent worker returned an invalid RSS scope"
                )
            if result.get("status") == "complete" and (
                type(result.get("decoder_cache_hit")) is not bool
                or any(
                    adaptive_runner._safe_elapsed(result.get(name)) is None
                    for name in ("decoder_build_seconds", "decode_seconds")
                )
            ):
                raise adaptive_runner.IdentityMismatchError(
                    "persistent worker returned incomplete cache timing evidence"
                )
            if result.get("status") == "fatal_identity":
                raise adaptive_runner.IdentityMismatchError(
                    "persistent worker rejected record provenance: "
                    + str(result.get("error", "identity failure"))
                )
        except BaseException:
            self._terminate_worker()
            raise
        maximum = execution_identity["resource_policy"]["max_rss_bytes"]
        if (
            result.get("status") == "memory_limit"
            or result.get("peak_rss_bytes", 0) > maximum
        ):
            self._terminate_worker()
        return result

    def run_tn(
        self,
        record: Mapping[str, Any],
        *,
        chi: int,
        cutoff: float,
        execution_identity: Mapping[str, Any],
    ) -> Mapping[str, Any]:
        return self._request(
            "tn",
            record,
            chi=chi,
            cutoff=cutoff,
            execution_identity=execution_identity,
        )

    def run_exact(
        self,
        record: Mapping[str, Any],
        *,
        execution_identity: Mapping[str, Any],
    ) -> Mapping[str, Any]:
        return self._request(
            "exact",
            record,
            chi=1,
            cutoff=0.0,
            execution_identity=execution_identity,
        )

    def close(self) -> None:
        failure = None
        if self._process is not None and self._process.poll() is None:
            try:
                self._process.stdin.write(f"CLOSE {self._session_token}\n")
                self._process.stdin.flush()
                returncode = self._process.wait(timeout=5.0)
                if returncode != 0:
                    failure = RuntimeError(
                        f"persistent worker close exited {returncode}"
                    )
            except (AttributeError, BrokenPipeError, OSError, subprocess.TimeoutExpired) as error:
                failure = RuntimeError(f"persistent worker close failed: {error}")
                self._terminate_worker()
        self._terminate_worker()
        self._discard_session()
        if failure is not None:
            raise failure

    def __enter__(self):
        return self

    def __exit__(self, _type, _value, _traceback):
        self.close()
        return False


def _probability_records(
    execution: Mapping[str, Any], p_key: str
) -> list[tuple[dict[str, Any], dict[str, Any]]]:
    if p_key not in manifest.P_VALUES:
        raise ValueError(f"p_key must be one of {manifest.P_VALUES}")
    records = [
        (task, record)
        for task in execution["tasks"]
        if task["p_key"] == p_key
        for record in task["records"]
    ]
    records.sort(key=lambda item: (item[0]["task_id"], item[1]["shot_index"]))
    if not records:
        raise ValueError("execution manifest has no records for the requested p")
    return records


def plan_probability_task(
    root: os.PathLike[str] | str,
    sampling_manifest: Mapping[str, Any],
    *,
    p_key: str,
) -> dict[str, Any]:
    root_path = Path(root).expanduser().resolve()
    execution = execution_manifest.load_execution_manifest(root_path, sampling_manifest)
    records = _probability_records(execution, p_key)
    task_ids = sorted({task["task_id"] for task, _record in records})
    L_values = {task["L"] for task, _record in records}
    if len(L_values) != 1:
        raise ValueError("one probability task must have a single L")
    return {
        "status": "dry_run",
        "execution_manifest_id": execution["execution_manifest_id"],
        "stage": execution["stage"],
        "L": next(iter(L_values)),
        "p_key": p_key,
        "task_ids": task_ids,
        "requested_records": len(records),
        "records": [
            {
                "task_id": task["task_id"],
                "shot_index": record["shot_index"],
                "record_id": record["record_id"],
            }
            for task, record in records
        ],
    }


def _task_lock(root: Path, execution_id: str, L: int, p_key: str):
    lock_dir = root / "locks"
    lock_dir.mkdir(parents=True, exist_ok=True)
    path = lock_dir / f"{execution_id}-L{L}-p{p_key}.lock"
    stream = path.open("a+", encoding="utf-8")
    try:
        fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError as error:
        stream.close()
        raise RuntimeError(f"probability task is already running: {path}") from error
    return stream


def run_probability_task(
    root: os.PathLike[str] | str,
    sampling_manifest: Mapping[str, Any],
    *,
    p_key: str,
    dry_run: bool = False,
    executor: Any | None = None,
    repository_root: os.PathLike[str] | str | None = None,
    julia: str = "julia",
) -> dict[str, Any]:
    root_path = Path(root).expanduser().resolve()
    plan = plan_probability_task(root_path, sampling_manifest, p_key=p_key)
    if dry_run:
        return plan
    execution = execution_manifest.load_execution_manifest(root_path, sampling_manifest)
    records = _probability_records(execution, p_key)
    owned_executor = executor or PersistentJuliaExecutor(
        repository_root=repository_root or Path(__file__).resolve().parents[3],
        julia=julia,
    )
    lock = _task_lock(
        root_path, execution["execution_manifest_id"], plan["L"], p_key
    )
    decisions = []
    try:
        for task, record in records:
            decisions.append(
                adaptive_runner.run_record(
                    root_path,
                    sampling_manifest,
                    task_id=task["task_id"],
                    record_id=record["record_id"],
                    executor=owned_executor,
                )
            )
    finally:
        close = getattr(owned_executor, "close", None)
        try:
            if callable(close):
                close()
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
            lock.close()
    accepted = sum(item["final_status"] == "accepted" for item in decisions)
    unknown = sum(item["final_status"] == "unknown" for item in decisions)
    false_accepts = sum(
        item.get("exact_audit", {}).get("false_accept") is True
        for item in decisions
    )
    return {
        "status": "complete",
        "execution_manifest_id": execution["execution_manifest_id"],
        "stage": execution["stage"],
        "L": plan["L"],
        "p_key": p_key,
        "task_ids": plan["task_ids"],
        "requested_records": len(decisions),
        "accepted_records": accepted,
        "unknown_records": unknown,
        "false_accepts": false_accepts,
        "decision_ids": [item["decision_id"] for item in decisions],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Plan or execute one persistent adaptive (L,p) task"
    )
    parser.add_argument("--study-root", required=True)
    parser.add_argument("--sampling-manifest", required=True)
    parser.add_argument("--p-key", required=True, choices=manifest.P_VALUES)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true")
    mode.add_argument("--execute-task", action="store_true")
    parser.add_argument(
        "--repository-root", default=str(Path(__file__).resolve().parents[3])
    )
    parser.add_argument("--julia", default="julia")
    arguments = parser.parse_args(argv)
    sampling = adaptive_runner._sampling_manifest_from_file(
        arguments.sampling_manifest
    )
    result = run_probability_task(
        arguments.study_root,
        sampling,
        p_key=arguments.p_key,
        dry_run=arguments.dry_run,
        repository_root=arguments.repository_root,
        julia=arguments.julia,
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
