#!/usr/bin/env python3
"""Cross-runner fault tests for bounded Linux descendant containment."""

from __future__ import annotations

import gc
import importlib.util
import os
import signal
import subprocess
import sys
import tempfile
import time
import unittest
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable
from unittest import mock


HERE = Path(__file__).resolve().parent


def load_module(name: str, path: Path) -> Any:
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load containment runner: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


MAIN = load_module(
    "leopard2_containment_main_compare",
    HERE / "main_compare" / "run_abba.py")
LOW = load_module(
    "leopard2_containment_low_copy",
    HERE / "low_encode_copy" / "run_abba.py")
BUTTERFLY = load_module(
    "leopard2_containment_butterfly",
    HERE / "backend_butterfly" / "run_abba.py")


@dataclass(frozen=True)
class Runner:
    name: str
    module: Any
    containment: Any
    invoke: Callable[[list[str], float], Any]
    returncode: Callable[[Any], int]


ENVIRONMENT = dict(os.environ)
RUNNERS = (
    Runner(
        "main_compare", MAIN, MAIN,
        lambda command, timeout: MAIN.run_process_bounded(
            command, environment=ENVIRONMENT, timeout=timeout,
            max_stdout=4096, max_stderr=4096),
        lambda completed: completed.returncode),
    Runner(
        "low_encode_copy", LOW, LOW.SUPPORT,
        lambda command, timeout: LOW.run_bounded(
            command, ENVIRONMENT, timeout),
        lambda completed: completed[0]),
    Runner(
        "backend_butterfly", BUTTERFLY, BUTTERFLY,
        lambda command, timeout: BUTTERFLY.run_benchmark_bounded(
            command, ENVIRONMENT, timeout),
        lambda completed: completed.returncode),
)


def escaping_program(marker_delay: float, close_output: bool = False) -> str:
    output = (
        " null=os.open('/dev/null', os.O_WRONLY)\n"
        " os.dup2(null, 1);os.dup2(null, 2);os.close(null)\n"
        if close_output else ""
    )
    leader_tail = (
        "deadline=time.monotonic()+5\n"
        "while not os.path.exists(sys.argv[2]) and time.monotonic()<deadline:\n"
        " time.sleep(.005)\n"
    )
    if not close_output:
        leader_tail += "while True: time.sleep(1)\n"
    return (
        "import os,pathlib,sys,time\n"
        "pid=os.fork()\n"
        "if pid == 0:\n"
        " os.setsid()\n"
        " daemon=os.fork()\n"
        " if daemon != 0: os._exit(0)\n" +
        output +
        " pathlib.Path(sys.argv[2]).write_text(str(os.getpid()))\n"
        f" time.sleep({marker_delay!r})\n"
        " pathlib.Path(sys.argv[1]).write_text('escaped')\n"
        " os._exit(0)\n" +
        leader_tail)


@unittest.skipUnless(sys.platform.startswith("linux"),
                     "Linux subreaper/procfs/pidfd tests")
class ProcessContainmentFaultTests(unittest.TestCase):
    def fault_patch(self, runner: Runner, fault: str, ready: Path):
        module = runner.containment
        if fault == "observe_spawn":
            def fail_observe(_self: Any, _process: Any) -> None:
                deadline = time.monotonic() + 1.5
                while not ready.exists() and time.monotonic() < deadline:
                    time.sleep(0.005)
                raise module.EvidenceError("injected spawn observation failure")
            return mock.patch.object(
                module.LinuxDescendantContainment,
                "observe_spawn", fail_observe)
        if fault == "attach":
            def fail_attach(_self: Any, _pid: int) -> None:
                deadline = time.monotonic() + 1.5
                while not ready.exists() and time.monotonic() < deadline:
                    time.sleep(0.005)
                raise module.EvidenceError("injected attach failure")
            return mock.patch.object(
                module.LinuxDescendantContainment, "attach", fail_attach)
        if fault == "procfs_snapshot":
            original = module._proc_process_snapshot
            calls = 0

            def fail_after_entry() -> Any:
                nonlocal calls
                calls += 1
                if calls > 1:
                    raise module.EvidenceError("injected procfs snapshot failure")
                return original()
            return mock.patch.object(
                module, "_proc_process_snapshot", fail_after_entry)
        if fault == "pidfd_open":
            original = module._linux_pidfd_open

            def fail_child(pid: int) -> int | None:
                if pid == os.getpid():
                    return original(pid)
                raise module.EvidenceError("injected pidfd_open failure")
            return mock.patch.object(module, "_linux_pidfd_open", fail_child)
        if fault == "pidfd_send_signal":
            original = module._linux_pidfd_signal

            def fail_kill(descriptor: int, number: int) -> None:
                if number == 0:
                    original(descriptor, number)
                    return
                raise module.EvidenceError("injected pidfd_send_signal failure")
            return mock.patch.object(module, "_linux_pidfd_signal", fail_kill)
        if fault == "primary_teardown":
            def fail_teardown(_self: Any, _process: Any) -> None:
                raise module.EvidenceError("injected primary teardown failure")
            return mock.patch.object(
                module.LinuxDescendantContainment,
                "terminate_and_reap", fail_teardown)
        raise AssertionError("unknown fault " + fault)

    def test_post_spawn_fault_matrix_cleans_every_descendant(self) -> None:
        faults = (
            "observe_spawn", "attach", "procfs_snapshot", "pidfd_open",
            "pidfd_send_signal", "primary_teardown")
        combinations = 0
        for runner in RUNNERS:
            for fault in faults:
                combinations += 1
                with self.subTest(runner=runner.name, fault=fault), \
                     tempfile.TemporaryDirectory(
                         prefix=f"leo2-{runner.name}-{fault}-") as directory:
                    root = Path(directory)
                    marker = root / "delayed-marker"
                    ready = root / "daemon-pid"
                    command = [
                        sys.executable, "-c", escaping_program(0.8),
                        str(marker), str(ready)]
                    before = runner.containment._emergency_get_child_subreaper()
                    started = time.monotonic()
                    with warnings.catch_warnings(record=True) as caught:
                        warnings.simplefilter("always", ResourceWarning)
                        with self.fault_patch(runner, fault, ready), \
                             self.assertRaises(Exception) as raised:
                            runner.invoke(command, 0.25)
                        gc.collect()
                    self.assertIn("injected", str(raised.exception))
                    elapsed = time.monotonic() - started
                    self.assertLess(elapsed, 7.0)
                    self.assertTrue(ready.is_file())
                    daemon_pid = int(ready.read_text(encoding="utf-8"))
                    self.assertEqual(
                        runner.containment._emergency_get_child_subreaper(), before)
                    self.assertFalse(Path("/proc", str(daemon_pid)).exists())
                    time.sleep(0.85)
                    self.assertFalse(marker.exists())
                    resources = [item for item in caught
                                 if issubclass(item.category, ResourceWarning)]
                    self.assertEqual(resources, [])
        # This is a strict superset of the audit's twelve runner/fault pairs:
        # all six post-spawn fault classes run through all three real runners.
        self.assertEqual(combinations, 18)

    def test_successful_leader_exit_cannot_leave_daemon(self) -> None:
        for runner in RUNNERS:
            with self.subTest(runner=runner.name), tempfile.TemporaryDirectory(
                    prefix=f"leo2-{runner.name}-successful-daemon-") as directory:
                root = Path(directory)
                marker = root / "delayed-marker"
                ready = root / "daemon-pid"
                command = [
                    sys.executable, "-c", escaping_program(0.8, close_output=True),
                    str(marker), str(ready)]
                before = runner.containment._emergency_get_child_subreaper()
                completed = runner.invoke(command, 2.0)
                self.assertEqual(runner.returncode(completed), 0)
                daemon_pid = int(ready.read_text(encoding="utf-8"))
                self.assertEqual(
                    runner.containment._emergency_get_child_subreaper(), before)
                self.assertFalse(Path("/proc", str(daemon_pid)).exists())
                time.sleep(0.85)
                self.assertFalse(marker.exists())

    def test_pid_reuse_identity_mismatch_never_signals_unrelated_process(self) -> None:
        for runner in RUNNERS:
            module = runner.containment
            with self.subTest(runner=runner.name):
                unrelated = subprocess.Popen(
                    [sys.executable, "-c", "import time;time.sleep(30)"],
                    stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL, start_new_session=True)
                try:
                    record = module._emergency_proc_process_record(unrelated.pid)
                    self.assertIsNotNone(record)
                    if record is None:
                        self.fail(
                            "process record is None after non-None assertion")
                    real_descriptor = module._emergency_pidfd_open(unrelated.pid)
                    self.assertIsNotNone(real_descriptor)
                    if real_descriptor is None:
                        self.fail(
                            "pidfd is None after non-None assertion")
                    reused = tuple(record[:3]) + (record[3] + 1,) + tuple(record[4:])
                    with mock.patch.object(
                            module, "_emergency_proc_process_record",
                            side_effect=[record, reused]), \
                         mock.patch.object(
                             module, "_emergency_pidfd_open",
                             side_effect=lambda _pid: os.dup(real_descriptor)), \
                         mock.patch.object(
                             module, "_emergency_pidfd_signal") as send_signal:
                        module._emergency_signal_identity(
                            (unrelated.pid, record[3]))
                    send_signal.assert_not_called()
                    self.assertIsNone(unrelated.poll())
                    os.close(real_descriptor)
                finally:
                    if unrelated.poll() is None:
                        unrelated.send_signal(signal.SIGKILL)
                    unrelated.wait(timeout=5.0)

    def test_exact_prior_subreaper_state_is_restored_for_zero_and_one(self) -> None:
        original = MAIN._emergency_get_child_subreaper()
        try:
            for desired in (0, 1):
                MAIN._emergency_restore_child_subreaper(desired)
                for runner in RUNNERS:
                    with self.subTest(desired=desired, runner=runner.name):
                        completed = runner.invoke(
                            [sys.executable, "-c", "pass"], 2.0)
                        self.assertEqual(runner.returncode(completed), 0)
                        self.assertEqual(
                            runner.containment._emergency_get_child_subreaper(),
                            desired)
        finally:
            MAIN._emergency_restore_child_subreaper(original)

    def test_every_post_spawn_exception_invokes_emergency_cleanup(self) -> None:
        for runner in RUNNERS:
            module = runner.containment
            with self.subTest(runner=runner.name):
                original = module.LinuxDescendantContainment \
                    .emergency_terminate_and_reap
                calls = 0

                def observed(containment: Any) -> None:
                    nonlocal calls
                    calls += 1
                    original(containment)

                before = module._emergency_get_child_subreaper()
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always", ResourceWarning)
                    with mock.patch.object(
                            module.LinuxDescendantContainment,
                            "emergency_terminate_and_reap", observed), \
                         mock.patch.object(
                             module.os, "set_blocking",
                             side_effect=module.EvidenceError(
                                 "injected post-spawn setup failure")), \
                         self.assertRaises(Exception) as raised:
                        runner.invoke(
                            [sys.executable, "-c",
                             "import time;time.sleep(30)"], 2.0)
                    gc.collect()
                self.assertIn("injected post-spawn", str(raised.exception))
                self.assertEqual(calls, 1)
                self.assertEqual(
                    module._emergency_get_child_subreaper(), before)
                resources = [item for item in caught
                             if issubclass(item.category, ResourceWarning)]
                self.assertEqual(resources, [])


if __name__ == "__main__":
    unittest.main(verbosity=2)
