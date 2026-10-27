"""Deterministic supervisor lifecycle tests; all child processes are simulated."""
import importlib.util
import json
from pathlib import Path
import signal
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch


SPEC = importlib.util.spec_from_file_location("tnml_supervisor", Path(__file__).with_name("supervise.py"))
SUPERVISOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SUPERVISOR)


class SupervisorTests(unittest.TestCase):
    def run_supervisor(self, *, finish_at=2, decode_limit=3, total_limit=30,
                       wall_jump=0, stale_pid=False, stale_same_pid=False,
                       stages=None, rss_kib=1):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            clock = SimpleNamespace(elapsed=0.0, launched=False)
            process = SimpleNamespace(pid=43210, returncode=None)
            def poll():
                if process.returncode is None and clock.elapsed >= finish_at:
                    process.returncode = 0
                return process.returncode
            def wait(timeout=None):
                self.assertIsNotNone(poll(), "simulated worker must finish or be terminated")
                return process.returncode
            process.poll, process.wait = poll, wait
            def launch(*args, **kwargs):
                clock.launched = True
                return process
            def killpg(pid, sig):
                self.assertEqual(pid, process.pid)
                self.assertIn(sig, (signal.SIGTERM, signal.SIGKILL))
                process.returncode = -sig
            def sleep(seconds):
                clock.elapsed += seconds
            def progress(_directory):
                pid = process.pid - 1 if stale_pid or not clock.launched and not stale_same_pid else process.pid
                stage, shot, marker = (stages(clock.elapsed) if stages else ("tn", 1, 1000.0))
                return dict(pid=pid, stage=stage, shot=shot, stage_started_unix=marker)
            argv = ["supervise.py", "--output", str(root), "--decode-seconds", str(decode_limit),
                    "--total-seconds", str(total_limit), "--rss-gib", "4", "--", "fake-julia"]
            with patch.object(SUPERVISOR.subprocess, "Popen", side_effect=launch), \
                 patch.object(SUPERVISOR.subprocess, "run", return_value=SimpleNamespace(returncode=0, stdout=str(rss_kib))), \
                 patch.object(SUPERVISOR.os, "killpg", side_effect=killpg), \
                 patch.object(SUPERVISOR.time, "monotonic", side_effect=lambda: clock.elapsed), \
                 patch.object(SUPERVISOR.time, "time", side_effect=lambda: 1000 + clock.elapsed + (wall_jump if clock.elapsed >= 1 else 0)), \
                 patch.object(SUPERVISOR.time, "sleep", side_effect=sleep), \
                 patch.object(SUPERVISOR, "read_progress", side_effect=progress), \
                 patch("sys.argv", argv):
                code = SUPERVISOR.main()
            return code, json.loads((root / "supervision.json").read_text())

    def test_forward_wall_clock_jump_does_not_kill_short_contraction(self):
        code, result = self.run_supervisor(wall_jump=10000)
        self.assertEqual(code, 0)
        self.assertEqual(result["status"], "completed")
        self.assertEqual(result["elapsed_seconds"], 2)
        self.assertIn("monotonic", result["timing_basis"])

    def test_actual_monotonic_contraction_over_limit_is_stopped(self):
        code, result = self.run_supervisor(finish_at=10)
        self.assertEqual(code, 2)
        self.assertEqual(result["status"], "contraction_time_limit")
        self.assertGreater(result["elapsed_seconds"], 3)
        self.assertLessEqual(result["elapsed_seconds"], 4)

    def test_backward_wall_clock_jump_does_not_disable_timeout(self):
        code, result = self.run_supervisor(finish_at=10, wall_jump=-10000)
        self.assertEqual(code, 2)
        self.assertEqual(result["status"], "contraction_time_limit")
        self.assertEqual(result["elapsed_seconds"], 4)

    def test_stale_progress_from_another_pid_does_not_start_timer(self):
        code, result = self.run_supervisor(finish_at=5, stale_pid=True)
        self.assertEqual(code, 0)
        self.assertEqual(result["status"], "completed")

    def test_stale_snapshot_with_reused_pid_needs_a_new_marker(self):
        code, result = self.run_supervisor(finish_at=5, stale_same_pid=True)
        self.assertEqual(code, 0)
        self.assertEqual(result["status"], "completed")

    def test_new_owned_tn_stage_gets_a_new_monotonic_timer(self):
        def stages(t):
            return ("tn", 1, 1000) if t < 2 else (("setup", 2, 1002) if t == 2 else ("tn", 2, 1003))
        code, result = self.run_supervisor(finish_at=6, stages=stages)
        self.assertEqual(code, 0)
        self.assertEqual(result["status"], "completed")

    def test_total_time_guard_still_covers_stale_progress(self):
        code, result = self.run_supervisor(finish_at=10, total_limit=3, stale_pid=True)
        self.assertEqual(code, 2)
        self.assertEqual(result["status"], "total_time_limit")
        self.assertEqual(result["elapsed_seconds"], 4)

    def test_memory_guard_still_covers_owned_worker(self):
        code, result = self.run_supervisor(rss_kib=5 * 1024**2)
        self.assertEqual(code, 2)
        self.assertEqual(result["status"], "resident_memory_limit")


if __name__ == "__main__":
    unittest.main()
