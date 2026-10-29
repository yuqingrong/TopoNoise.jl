#!/usr/bin/env python3
"""Bounded persistent-worker tests. No Julia contraction is launched."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import adaptive_runner
import task_runner
from test_phase2 import FakeExecutor, make_bundle, runtime_identity


def _runtime_toml_lines(runtime):
    lines = ["", "[runtime_identity]"]
    for key in (
        "julia_version",
        "blas_vendor",
        "cpu_model",
        "julia_threads",
        "blas_threads",
    ):
        value = runtime[key]
        encoded = json.dumps(value) if isinstance(value, str) else str(value)
        lines.append(f"{key} = {encoded}")
    return lines


class FakeClock:
    def __init__(self):
        self.value = 0.0

    def monotonic(self):
        return self.value

    def sleep(self, seconds):
        self.value += max(float(seconds), 0.001)


class FakeStdin:
    def __init__(self, process):
        self.process = process
        self.closed = False

    def write(self, value):
        if self.closed:
            raise BrokenPipeError("fake worker stdin is closed")
        if self.process.broken_pipe:
            raise BrokenPipeError("fake worker died before request authentication")
        self.process.handle_line(value.rstrip("\n"))
        return len(value)

    def flush(self):
        return None

    def close(self):
        self.closed = True


class FakeWorkerProcess:
    def __init__(
        self,
        command,
        runtime,
        *,
        hang=False,
        corrupt_response=False,
        crash=False,
        fatal_identity=False,
        result_status="complete",
        rss_scope="worker_lifetime_high_water",
        close_returncode=0,
        broken_pipe=False,
        malformed_result=False,
        stdout=None,
        **_kwargs,
    ):
        self.command = list(command)
        self.runtime = runtime
        self.hang = hang
        self.corrupt_response = corrupt_response
        self.crash = crash
        self.fatal_identity = fatal_identity
        self.result_status = result_status
        self.rss_scope = rss_scope
        self.close_returncode = close_returncode
        self.broken_pipe = broken_pipe
        self.malformed_result = malformed_result
        self.returncode = None
        self.terminated = False
        self.killed = False
        self.stdin = FakeStdin(self)
        self.requests = []
        self.cached_chi = set()
        self.stdout_target = stdout
        self.environment = dict(_kwargs.get("env", {}))
        if stdout is not None:
            stdout.write("backend startup chatter that is not protocol data\n")
            stdout.flush()
        ready = Path(self.command[self.command.index("--ready") + 1])
        L = int(self.command[self.command.index("--L") + 1])
        p = float(self.command[self.command.index("--p") + 1])
        lines = [
            'status = "ready"',
            "valid = true",
            "protocol_version = 1",
            f'session_token = {json.dumps(self.command[self.command.index("--session-token") + 1])}',
            f'source_bundle_id = {json.dumps(self.command[self.command.index("--source-bundle-id") + 1])}',
            f'runtime_profile_id = {json.dumps(self.command[self.command.index("--runtime-profile-id") + 1])}',
            f"L = {L}",
            f"p = {p!r}",
            f'p_key = {json.dumps(self.command[self.command.index("--p-key") + 1])}',
            f'cutoff = {self.command[self.command.index("--cutoff") + 1]}',
            'executor_mode = "persistent"',
            'wall_time_scope = "request_after_worker_ready"',
            'decoder_cache_scope = "worker_lifetime_by_chi"',
            "startup_elapsed_seconds = 0.25",
            "peak_rss_bytes = 1024",
            'rss_scope = "worker_lifetime_high_water"',
        ] + _runtime_toml_lines(self.runtime)
        ready.write_text("\n".join(lines) + "\n", encoding="utf-8")

    def handle_line(self, line):
        if line.startswith("CLOSE "):
            self.returncode = self.close_returncode
            return
        request = adaptive_runner._read_toml(Path(line))
        self.requests.append(request)
        if self.hang:
            return
        if self.crash:
            self.returncode = 3
            return
        mode = request["mode"]
        chi = int(request["chi"])
        cache_hit = mode == "tn" and chi in self.cached_chi
        mode == "tn" and self.cached_chi.add(chi)
        output = Path(request["output_path"])
        if self.malformed_result:
            output.write_text('status = "unterminated\n', encoding="utf-8")
            return
        status = "fatal_identity" if self.fatal_identity else self.result_status
        lines = [
            f'status = {json.dumps(status)}',
            f"valid = {'true' if status == 'complete' else 'false'}",
            "elapsed_seconds = 0.5",
            "peak_rss_bytes = 2048",
            f'rss_scope = {json.dumps(self.rss_scope)}',
            f'request_id = {json.dumps(request["request_id"])}',
            f"decoder_cache_hit = {'true' if cache_hit else 'false'}",
            "decoder_build_seconds = 0.125",
            "decode_seconds = 0.375",
        ]
        if status == "complete":
            lines.extend(
                [
                    "probabilities = [0.7, 0.1, 0.1, 0.1]",
                    "logical_mask = 0",
                    "log_joint_syndrome = -10.0",
                ]
            )
        else:
            lines.append(f'error = {json.dumps(status + " from fake worker")}')
        for key in (
            "protocol_version",
            "session_token",
            "sequence",
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
        ):
            value = request[key]
            if self.corrupt_response and key == "sequence":
                value += 1
            encoded = json.dumps(value) if isinstance(value, str) else str(value).lower()
            lines.insert(5, f"{key} = {encoded}")
        lines.extend(_runtime_toml_lines(self.runtime))
        output.write_text("\n".join(lines) + "\n", encoding="utf-8")

    def poll(self):
        return self.returncode

    def wait(self, timeout=None):
        if self.returncode is None:
            raise subprocess.TimeoutExpired(self.command, timeout)
        return self.returncode

    def terminate(self):
        self.terminated = True
        self.returncode = -15

    def kill(self):
        self.killed = True
        self.returncode = -9


class FakePopenFactory:
    def __init__(
        self,
        runtime,
        *,
        hang_first=False,
        corrupt_response=False,
        crash=False,
        fatal_identity=False,
        result_status="complete",
        rss_scope="worker_lifetime_high_water",
        close_returncode=0,
        broken_pipe=False,
        malformed_result=False,
        startup_delay=0.0,
        clock=None,
    ):
        self.runtime = runtime
        self.hang_first = hang_first
        self.corrupt_response = corrupt_response
        self.crash = crash
        self.fatal_identity = fatal_identity
        self.result_status = result_status
        self.rss_scope = rss_scope
        self.close_returncode = close_returncode
        self.broken_pipe = broken_pipe
        self.malformed_result = malformed_result
        self.startup_delay = float(startup_delay)
        self.clock = clock
        self.processes = []

    def __call__(self, command, **kwargs):
        if self.startup_delay:
            if self.clock is None:
                raise RuntimeError("startup delay test requires a fake clock")
            self.clock.sleep(self.startup_delay)
        process = FakeWorkerProcess(
            command,
            self.runtime,
            hang=self.hang_first and not self.processes,
            corrupt_response=self.corrupt_response,
            crash=self.crash,
            fatal_identity=self.fatal_identity,
            result_status=self.result_status,
            rss_scope=self.rss_scope,
            close_returncode=self.close_returncode,
            broken_pipe=self.broken_pipe,
            malformed_result=self.malformed_result,
            **kwargs,
        )
        self.processes.append(process)
        return process


def execution_identity(*, L=5, p=0.006, blas_threads=4):
    runtime = runtime_identity()
    runtime["blas_threads"] = blas_threads
    return {
        "source_bundle_id": "a" * 64,
        "runtime_profile_id": "b" * 64,
        "source_identity": {},
        "runtime_identity": runtime,
        "executor_policy": {
            "persistent_protocol": "persistent_toml_mailbox_v1",
            "protocol_version": 1,
            "scheduling": "serial_one_worker_per_manifest_task_or_L_p_group",
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
        },
        "resource_policy": {
            "max_rss_bytes": 16 * 2**30,
            "wall_seconds_by_chi": {
                "32": 0.1,
                "64": 0.2,
                "128": 0.3,
                "240": 0.4,
                "480": 0.5,
            },
            "exact_wall_seconds": 0.1,
            "preflight_wall_seconds": 1,
        },
        "L": L,
        "p": p,
        "p_key": f"{p:.4f}",
        "cutoff": 1e-18,
        "executor_mode": "persistent",
        "wall_time_scope": "request_after_worker_ready",
        "decoder_cache_scope": "worker_lifetime_by_chi",
        "record_id": "c" * 64,
        "julia_record_id": "a" * 64,
    }


class PersistentProtocolTests(unittest.TestCase):
    def make_executor(self, *, factory=None, clock=None):
        factory = factory or FakePopenFactory(runtime_identity())
        clock = clock or FakeClock()
        executor = task_runner.PersistentJuliaExecutor(
            repository_root=Path(__file__).resolve().parents[3],
            popen_factory=factory,
            monotonic=clock.monotonic,
            sleep=clock.sleep,
            poll_interval=0.01,
        )
        executor.verify_execution_environment = lambda _identity: None
        return executor, factory

    def test_worker_reuses_process_and_reports_decoder_cache_hit_despite_stdout_noise(self):
        executor, factory = self.make_executor()
        environment_checks = []
        executor.verify_execution_environment = lambda identity: environment_checks.append(
            identity["source_bundle_id"]
        )
        identity = execution_identity()
        record = {"record_id": "a" * 64}
        executor.preflight(identity)
        executor.preflight(identity)
        self.assertEqual(len(factory.processes), 1)
        self.assertEqual(environment_checks, [identity["source_bundle_id"]])
        first = executor.run_tn(
            record, chi=32, cutoff=1e-18, execution_identity=identity
        )
        second = executor.run_tn(
            record, chi=32, cutoff=1e-18, execution_identity=identity
        )
        self.assertFalse(first["decoder_cache_hit"])
        self.assertTrue(second["decoder_cache_hit"])
        self.assertEqual(first["rss_scope"], "worker_lifetime_high_water")
        self.assertEqual(len(factory.processes[0].requests), 2)
        executor.close()
        self.assertEqual(factory.processes[0].returncode, 0)

    def test_worker_startup_wall_is_charged_once_to_first_record(self):
        clock = FakeClock()
        factory = FakePopenFactory(
            runtime_identity(), startup_delay=2.5, clock=clock
        )
        executor, _ = self.make_executor(factory=factory, clock=clock)
        identity = execution_identity()
        record = {"record_id": "a" * 64}
        first = executor.run_tn(
            record, chi=32, cutoff=1e-18, execution_identity=identity
        )
        second = executor.run_tn(
            record, chi=64, cutoff=1e-18, execution_identity=identity
        )
        self.assertEqual(first["worker_startup_seconds"], 2.5)
        self.assertEqual(second["worker_startup_seconds"], 0.0)
        executor.close()

    def test_worker_uses_frozen_blas_thread_count_for_replay(self):
        runtime = runtime_identity()
        runtime["blas_threads"] = 1
        factory = FakePopenFactory(runtime)
        executor, _ = self.make_executor(factory=factory)
        executor.preflight(execution_identity(blas_threads=1))
        self.assertEqual(factory.processes[0].environment["TNML_BLAS_THREADS"], "1")
        self.assertEqual(factory.processes[0].environment["OPENBLAS_NUM_THREADS"], "1")
        executor.close()

    def test_timeout_terminates_worker_and_next_preflight_restarts_it(self):
        clock = FakeClock()
        factory = FakePopenFactory(runtime_identity(), hang_first=True)
        executor, _ = self.make_executor(factory=factory, clock=clock)
        environment_checks = []
        executor.verify_execution_environment = lambda identity: environment_checks.append(
            identity["source_bundle_id"]
        )
        identity = execution_identity()
        result = executor.run_tn(
            {"record_id": "a" * 64},
            chi=32,
            cutoff=1e-18,
            execution_identity=identity,
        )
        self.assertEqual(result["status"], "timeout")
        self.assertGreaterEqual(result["elapsed_seconds"], 0.1)
        self.assertTrue(factory.processes[0].terminated)
        executor.preflight(identity)
        self.assertEqual(len(factory.processes), 2)
        self.assertEqual(environment_checks, [identity["source_bundle_id"]] * 2)
        executor.close()

    def test_deadline_precedes_a_late_result_file_observation(self):
        class DelayedClock:
            def __init__(self):
                self.calls = 0

            def monotonic(self):
                self.calls += 1
                return 0.0 if self.calls == 1 else 2.0

        with tempfile.TemporaryDirectory() as temporary:
            result = Path(temporary) / "late-result.toml"
            result.write_text('status = "complete"\n', encoding="utf-8")
            clock = DelayedClock()
            executor = task_runner.PersistentJuliaExecutor(
                repository_root=Path(__file__).resolve().parents[3],
                monotonic=clock.monotonic,
                sleep=lambda _seconds: None,
            )
            self.assertEqual(executor._wait_for_file(result, 1.0), "timeout")

    def test_worker_refuses_identity_or_probability_change_while_live(self):
        executor, factory = self.make_executor()
        executor.preflight(execution_identity())
        changed = execution_identity(p=0.01)
        with self.assertRaisesRegex(adaptive_runner.IdentityMismatchError, "worker"):
            executor.preflight(changed)
        self.assertEqual(len(factory.processes), 1)
        executor.close()

    def test_worker_requires_the_complete_frozen_protocol_policy(self):
        executor, factory = self.make_executor()
        identity = execution_identity()
        identity["executor_policy"]["protocol_version"] = 2
        with self.assertRaisesRegex(
            adaptive_runner.IdentityMismatchError, "protocol"
        ):
            executor.preflight(identity)
        self.assertEqual(factory.processes, [])

    def test_response_identity_mismatch_terminates_worker(self):
        factory = FakePopenFactory(
            runtime_identity(), corrupt_response=True
        )
        executor, _ = self.make_executor(factory=factory)
        with self.assertRaisesRegex(
            adaptive_runner.IdentityMismatchError, "response identity"
        ):
            executor.run_tn(
                {"record_id": "a" * 64},
                chi=32,
                cutoff=1e-18,
                execution_identity=execution_identity(),
            )
        self.assertTrue(factory.processes[0].terminated)
        executor.close()

    def test_worker_exit_without_authenticated_result_is_fatal(self):
        factory = FakePopenFactory(runtime_identity(), crash=True)
        executor, _ = self.make_executor(factory=factory)
        with self.assertRaisesRegex(
            adaptive_runner.IdentityMismatchError, "authenticated result"
        ):
            executor.run_tn(
                {"record_id": "a" * 64},
                chi=32,
                cutoff=1e-18,
                execution_identity=execution_identity(),
            )
        executor.close()
        executor.close()

    def test_record_identity_failure_is_fatal_not_scientific_unknown(self):
        factory = FakePopenFactory(runtime_identity(), fatal_identity=True)
        executor, _ = self.make_executor(factory=factory)
        with self.assertRaisesRegex(
            adaptive_runner.IdentityMismatchError, "record provenance"
        ):
            executor.run_tn(
                {"record_id": "a" * 64},
                chi=32,
                cutoff=1e-18,
                execution_identity=execution_identity(),
            )
        self.assertTrue(factory.processes[0].terminated)
        executor.close()

    def test_persistent_response_requires_worker_rss_scope(self):
        factory = FakePopenFactory(
            runtime_identity(), rss_scope="process_lifetime_high_water"
        )
        executor, _ = self.make_executor(factory=factory)
        with self.assertRaisesRegex(
            adaptive_runner.IdentityMismatchError, "RSS scope"
        ):
            executor.run_tn(
                {"record_id": "a" * 64},
                chi=32,
                cutoff=1e-18,
                execution_identity=execution_identity(),
            )
        self.assertTrue(factory.processes[0].terminated)
        executor.close()

    def test_memory_limit_result_retires_worker_and_close_checks_exit(self):
        memory_factory = FakePopenFactory(
            runtime_identity(), result_status="memory_limit"
        )
        executor, _ = self.make_executor(factory=memory_factory)
        result = executor.run_tn(
            {"record_id": "a" * 64},
            chi=32,
            cutoff=1e-18,
            execution_identity=execution_identity(),
        )
        self.assertEqual(result["status"], "memory_limit")
        self.assertTrue(memory_factory.processes[0].terminated)
        executor.close()

        close_factory = FakePopenFactory(runtime_identity(), close_returncode=7)
        close_executor, _ = self.make_executor(factory=close_factory)
        close_executor.preflight(execution_identity())
        with self.assertRaisesRegex(RuntimeError, "close exited 7"):
            close_executor.close()


class TaskRunnerTests(unittest.TestCase):
    def test_malformed_worker_result_is_fatal_and_writes_no_checkpoint(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling, _, execution = make_bundle(root, "l3_exact")
            record = execution["tasks"][0]["records"][0]
            plan = adaptive_runner.plan_record(
                root, sampling, task_id=0, record_id=record["record_id"]
            )
            factory = FakePopenFactory(
                runtime_identity(), malformed_result=True
            )
            clock = FakeClock()
            executor = task_runner.PersistentJuliaExecutor(
                repository_root=Path(__file__).resolve().parents[3],
                popen_factory=factory,
                monotonic=clock.monotonic,
                sleep=clock.sleep,
            )
            executor.verify_execution_environment = lambda _identity: None
            with self.assertRaisesRegex(
                adaptive_runner.IdentityMismatchError, "unreadable result"
            ):
                adaptive_runner.run_record(
                    root,
                    sampling,
                    task_id=0,
                    record_id=record["record_id"],
                    executor=executor,
                )
            self.assertFalse(Path(plan["record_directory"]).exists())
            self.assertTrue(factory.processes[0].terminated)
            executor.close()

    def test_broken_request_pipe_is_fatal_and_writes_no_checkpoint(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling, _, execution = make_bundle(root, "l3_exact")
            record = execution["tasks"][0]["records"][0]
            plan = adaptive_runner.plan_record(
                root, sampling, task_id=0, record_id=record["record_id"]
            )
            factory = FakePopenFactory(runtime_identity(), broken_pipe=True)
            executor = task_runner.PersistentJuliaExecutor(
                repository_root=Path(__file__).resolve().parents[3],
                popen_factory=factory,
                monotonic=FakeClock().monotonic,
                sleep=FakeClock().sleep,
            )
            executor.verify_execution_environment = lambda _identity: None
            with self.assertRaisesRegex(
                adaptive_runner.IdentityMismatchError, "request pipe"
            ):
                adaptive_runner.run_record(
                    root,
                    sampling,
                    task_id=0,
                    record_id=record["record_id"],
                    executor=executor,
                )
            self.assertFalse(Path(plan["record_directory"]).exists())
            self.assertTrue(factory.processes[0].terminated)
            executor.close()

    def test_run_record_timeout_is_immutable_unknown_and_next_record_restarts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling, _, execution = make_bundle(root, "l3_exact")
            records = execution["tasks"][0]["records"][:2]
            clock = FakeClock()
            factory = FakePopenFactory(runtime_identity(), hang_first=True)
            executor = task_runner.PersistentJuliaExecutor(
                repository_root=Path(__file__).resolve().parents[3],
                popen_factory=factory,
                monotonic=clock.monotonic,
                sleep=clock.sleep,
                poll_interval=0.01,
            )
            executor.verify_execution_environment = lambda _identity: None
            first = adaptive_runner.run_record(
                root,
                sampling,
                task_id=0,
                record_id=records[0]["record_id"],
                executor=executor,
            )
            self.assertEqual(first["final_status"], "unknown")
            self.assertEqual(first["unknown_reason"], "exact_timeout")
            exact_path = Path(first["record_directory"]) / "exact.json"
            frozen_timeout = exact_path.read_bytes()

            second = adaptive_runner.run_record(
                root,
                sampling,
                task_id=0,
                record_id=records[1]["record_id"],
                executor=executor,
            )
            self.assertEqual(second["final_status"], "accepted")
            self.assertEqual(exact_path.read_bytes(), frozen_timeout)
            self.assertGreaterEqual(len(factory.processes), 2)
            executor.close()

    def test_probability_lock_rejects_concurrent_owner(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first = task_runner._task_lock(root, "a" * 64, 3, "0.0060")
            try:
                with self.assertRaisesRegex(RuntimeError, "already running"):
                    task_runner._task_lock(root, "a" * 64, 3, "0.0060")
            finally:
                task_runner.fcntl.flock(
                    first.fileno(), task_runner.fcntl.LOCK_UN
                )
                first.close()

    def test_dry_run_groups_every_record_for_one_probability_without_writes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling, _, _ = make_bundle(root, "l3_exact")
            before = sorted(str(path.relative_to(root)) for path in root.rglob("*"))
            plan = task_runner.run_probability_task(
                root, sampling, p_key="0.0060", dry_run=True
            )
            after = sorted(str(path.relative_to(root)) for path in root.rglob("*"))
            self.assertEqual(before, after)
            self.assertEqual(plan["status"], "dry_run")
            self.assertEqual(plan["L"], 3)
            self.assertEqual(plan["p_key"], "0.0060")
            self.assertEqual(plan["requested_records"], 100)
            self.assertEqual(len(plan["task_ids"]), 5)

    def test_manifest_task_mode_preserves_twenty_record_l3_shard(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling, _, _ = make_bundle(root, "l3_exact")
            plan = task_runner.run_manifest_task(
                root, sampling, task_id=0, dry_run=True
            )
            self.assertEqual(plan["task_ids"], [0])
            self.assertEqual(plan["requested_records"], 20)
            self.assertEqual(plan["p_key"], "0.0060")

    def test_one_executor_processes_and_resumes_entire_l5_probability_group(self):
        class ClosingFake(FakeExecutor):
            def __init__(self):
                super().__init__()
                self.close_calls = 0

            def close(self):
                self.close_calls += 1

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampling, _, _ = make_bundle(root, "l5_feasibility")
            executor = ClosingFake()
            report = task_runner.run_probability_task(
                root,
                sampling,
                p_key="0.0060",
                executor=executor,
            )
            self.assertEqual(report["requested_records"], 10)
            self.assertEqual(report["accepted_records"], 10)
            self.assertEqual(report["unknown_records"], 0)
            self.assertEqual(executor.close_calls, 1)
            self.assertGreater(len(executor.tn_calls), 10)

            resumed = ClosingFake()
            resumed_report = task_runner.run_probability_task(
                root,
                sampling,
                p_key="0.0060",
                executor=resumed,
            )
            self.assertEqual(resumed_report, report)
            self.assertEqual(resumed.tn_calls, [])
            self.assertEqual(resumed.exact_calls, [])
            self.assertEqual(resumed.close_calls, 1)


if __name__ == "__main__":
    unittest.main()
