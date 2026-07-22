#!/usr/bin/env python3
"""Fail-closed scheduler-isolation tests for run_small_direct_abba.py."""

from __future__ import annotations

import argparse
import fcntl
import importlib.util
import json
import os
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock


RUNNER_PATH = Path(__file__).with_name("run_small_direct_abba.py")
SPEC = importlib.util.spec_from_file_location(
    "leopard2_test_direct_odd_runner", RUNNER_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot import direct-repair ABBA runner")
RUNNER = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RUNNER
SPEC.loader.exec_module(RUNNER)


def record() -> dict:
    return {
        "pid": 101, "tid": 102, "command": "worker",
        "affinity": [0, 1, 2], "start_time_ticks": 12345,
        "task_inode": 777,
    }


def runner_record(affinity: list[int]) -> dict:
    return {
        "pid": 1001, "tid": 1001, "affinity": affinity,
        "start_time_ticks": 55555, "task_inode": 999,
    }


class IsolationTests(unittest.TestCase):
    def test_small_direct_matrices_are_complete_and_stable(self) -> None:
        core = RUNNER.make_matrix()
        self.assertEqual(core["schema"], RUNNER.MATRIX_SCHEMA)
        self.assertEqual(core["cell_count"], 160)
        self.assertEqual(
            core["matrix_sha256"], RUNNER.make_matrix()["matrix_sha256"])
        self.assertEqual(
            {cell["K"] for cell in core["cells"]}, {5, 8, 9, 12, 16})
        self.assertEqual(
            {cell["bytes"] for cell in core["cells"]},
            {64, 2048, 4096, 65536})
        self.assertEqual(
            [cell["index"] for cell in core["cells"]],
            list(range(core["cell_count"])))
        self.assertEqual(
            len({cell["id"] for cell in core["cells"]}),
            core["cell_count"])
        self.assertTrue(all(
            cell["loss"] <= min(cell["K"], cell["R"])
            for cell in core["cells"]))

        full = RUNNER.make_large_matrix()
        self.assertEqual(full["schema"], RUNNER.LARGE_MATRIX_SCHEMA)
        self.assertEqual(full["cell_count"], 1264)
        self.assertEqual(
            full["matrix_sha256"],
            RUNNER.make_large_matrix()["matrix_sha256"])
        self.assertEqual(
            [cell["index"] for cell in full["cells"]],
            list(range(full["cell_count"])))
        self.assertEqual(
            len({cell["id"] for cell in full["cells"]}),
            full["cell_count"])
        expected = {
            (k, r, loss, byte_count)
            for k in range(5, 17)
            for r in range(5, 9)
            for loss in range(4, min(k, r) + 1)
            for byte_count in
                (64, 65, 256, 1024, 2048, 2049, 4096, 65536)
        }
        self.assertEqual({
            (cell["K"], cell["R"], cell["loss"], cell["bytes"])
            for cell in full["cells"]
        }, expected)

    def test_small_direct_modes_are_explicit_and_directional(self) -> None:
        self.assertEqual(
            RUNNER.comparison_modes("transform", "output"),
            {"baseline": "transform", "candidate": "output"})
        self.assertEqual(
            RUNNER.comparison_modes("transform", "source"),
            {"baseline": "transform", "candidate": "source"})
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "distinct and known"):
            RUNNER.comparison_modes("source", "source")
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "distinct and known"):
            RUNNER.comparison_modes("unknown", "source")

        loss4 = {"loss": 4}
        loss5 = {"loss": 5}
        self.assertEqual(
            RUNNER.expected_direct_executor("transform", loss4),
            "output_major")
        self.assertEqual(
            RUNNER.expected_direct_executor("output", loss5),
            "output_major")
        self.assertEqual(
            RUNNER.expected_direct_executor("source", loss5),
            "source_major")
        self.assertEqual(
            RUNNER.expected_direct_executor("transform", loss5), "none")

        arguments = ["c++", "-O3", RUNNER.MODE_COMPILE_DEFINITIONS["source"]]
        self.assertEqual(
            RUNNER.strip_mode_definition(arguments, "source", "test"),
            ["c++", "-O3"])
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "exact mode definition once"):
            RUNNER.strip_mode_definition(["c++", "-O3"], "source", "test")

    @staticmethod
    def physical_pairs() -> list[tuple[int, int]]:
        allowed = RUNNER.cgroup_effective_cpus() & set(os.sched_getaffinity(0))
        pairs = []
        seen = set()
        for cpu in sorted(allowed):
            path = Path(
                "/sys/devices/system/cpu/cpu%d/topology/"
                "thread_siblings_list" % cpu)
            siblings = RUNNER.parse_cpu_list(path.read_text()) & allowed
            if len(siblings) == 2:
                sibling = next(value for value in siblings if value != cpu)
                pair = tuple(sorted((cpu, sibling)))
                if pair not in seen:
                    seen.add(pair)
                    pairs.append(pair)
        return pairs

    @classmethod
    def physical_pair(cls) -> tuple[int, int]:
        pairs = cls.physical_pairs()
        if pairs:
            return pairs[0]
        raise unittest.SkipTest("no process-visible two-thread SMT pair")

    @classmethod
    def free_physical_pair(cls, source_root: Path) -> tuple[int, int]:
        PairLease, unused = RUNNER.load_pair_lease(source_root)
        for cpu, sibling in cls.physical_pairs():
            try:
                with PairLease(cpu, sibling):
                    pass
                return cpu, sibling
            except Exception:
                continue
        raise unittest.SkipTest("no unleased process-visible SMT pair")

    def test_exact_masks_are_excluded_and_restored(self) -> None:
        state = {102: {0, 1, 2}}

        def get_affinity(tid: int) -> set[int]:
            return set(state[tid])

        def set_affinity(tid: int, cpus: set[int]) -> None:
            state[tid] = set(cpus)

        with mock.patch.object(
                RUNNER, "cgroup_effective_cpus",
                return_value={0, 1, 2, 3}), \
             mock.patch.object(
                RUNNER, "same_user_pair_affinity",
                side_effect=[[record()], [], []]), \
             mock.patch.object(
                RUNNER, "runner_affinity_record",
                side_effect=[runner_record([0, 1, 2, 3]),
                             runner_record([0, 2])]), \
             mock.patch.object(RUNNER, "require_current_task"), \
             mock.patch.object(
                RUNNER.os, "sched_getaffinity", side_effect=get_affinity), \
             mock.patch.object(
                RUNNER.os, "sched_setaffinity", side_effect=set_affinity):
            exclusion = RUNNER.exclude_same_user_from_pair(1, 3)
        self.assertEqual(exclusion["changed"][0]["affinity"], [0, 1, 2])
        self.assertEqual(exclusion["changed"][0]["after"], [0, 2])
        self.assertEqual(state[102], {0, 2})

        with mock.patch.object(RUNNER, "require_current_task"), \
             mock.patch.object(
                RUNNER, "runner_affinity_record",
                return_value=runner_record([0, 1, 2, 3])), \
             mock.patch.object(
                RUNNER.os, "sched_getaffinity", side_effect=get_affinity), \
             mock.patch.object(
                RUNNER.os, "sched_setaffinity", side_effect=set_affinity):
            restoration = RUNNER.restore_same_user_affinity(exclusion)
        self.assertEqual(state[102], {0, 1, 2})
        self.assertEqual(restoration["restored"][0]["affinity"], [0, 1, 2])
        self.assertEqual(restoration["failures"], [])

    def test_vanished_task_identity_is_rejected(self) -> None:
        with mock.patch.object(
                RUNNER, "current_task_identity",
                side_effect=FileNotFoundError("gone")):
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "vanished before mutation"):
                RUNNER.require_current_task(record())

    def test_pid_reused_task_identity_is_rejected(self) -> None:
        with mock.patch.object(
                RUNNER, "current_task_identity", return_value=(99999, 777)):
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "reused or replaced"):
                RUNNER.require_current_task(record())

    def test_same_tick_reused_task_inode_is_rejected(self) -> None:
        with mock.patch.object(
                RUNNER, "current_task_identity", return_value=(12345, 888)):
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "inode was reused or replaced"):
                RUNNER.require_current_task(record())

    def test_failed_affinity_mutation_rolls_back_and_rejects(self) -> None:
        with mock.patch.object(
                RUNNER, "cgroup_effective_cpus",
                return_value={0, 1, 2, 3}), \
             mock.patch.object(
                RUNNER, "same_user_pair_affinity",
                return_value=[record()]), \
             mock.patch.object(RUNNER, "require_current_task"), \
             mock.patch.object(
                RUNNER.os, "sched_getaffinity", return_value={0, 1, 2}), \
             mock.patch.object(
                RUNNER.os, "sched_setaffinity",
                side_effect=PermissionError("denied")):
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "failed to exclude"):
                RUNNER.exclude_same_user_from_pair(1, 3)

    def test_new_pair_eligible_task_is_captured_and_restored(self) -> None:
        state = {102: {0, 1, 2}, 202: {1}}
        foreign = {
            "pid": 201, "tid": 202, "command": "foreign",
            "affinity": [1], "start_time_ticks": 54321,
            "task_inode": 778,
        }

        def get_affinity(tid: int) -> set[int]:
            return set(state[tid])

        def set_affinity(tid: int, cpus: set[int]) -> None:
            state[tid] = set(cpus)

        with mock.patch.object(
                RUNNER, "cgroup_effective_cpus",
                return_value={0, 1, 2, 3}), \
             mock.patch.object(
                RUNNER, "same_user_pair_affinity",
                side_effect=[[record()], [foreign], [], []]), \
             mock.patch.object(
                RUNNER, "runner_affinity_record",
                side_effect=[runner_record([0, 1, 2, 3]),
                             runner_record([0, 2]),
                             runner_record([0, 1, 2, 3])]), \
             mock.patch.object(RUNNER, "require_current_task"), \
             mock.patch.object(
                RUNNER.os, "sched_getaffinity", side_effect=get_affinity), \
             mock.patch.object(
                RUNNER.os, "sched_setaffinity", side_effect=set_affinity):
            exclusion = RUNNER.exclude_same_user_from_pair(1, 3)
            self.assertEqual(len(exclusion["changed"]), 2)
            self.assertEqual(state[202], {0, 2})
            RUNNER.restore_same_user_affinity(exclusion)
        self.assertEqual(state[102], {0, 1, 2})
        self.assertEqual(state[202], {1})

    def test_restore_failure_is_rejected(self) -> None:
        changed = {**record(), "after": [0, 2]}
        with mock.patch.object(RUNNER, "require_current_task"), \
             mock.patch.object(
                RUNNER.os, "sched_getaffinity", return_value={0, 2}), \
             mock.patch.object(
                RUNNER.os, "sched_setaffinity",
                side_effect=PermissionError("denied")):
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "restoration failed"):
                RUNNER.restore_same_user_affinity({"changed": [changed]})

    def test_cgroup_eligibility_is_distinct_from_inherited_mask(self) -> None:
        identity = RUNNER.validate_smt_pair_identity(
            1, 3, {0, 1, 3}, {0}, {1, 3})
        self.assertEqual(identity["cgroup_effective_cpus"], [0, 1, 3])
        self.assertEqual(identity["launch_affinity"], [0])
        self.assertEqual(identity["housekeeping_affinity"], [0])

    def test_transient_target_runtime_is_rejected(self) -> None:
        evidence = RUNNER.target_runtime_evidence(
            1_000_000, 2_000_000,
            1_000_000 - RUNNER.TARGET_RUNTIME_TOLERANCE_NS - 1)
        self.assertFalse(evidence["accepted"])
        self.assertEqual(
            evidence["unexplained_runtime_ns"],
            RUNNER.TARGET_RUNTIME_TOLERANCE_NS + 1)

    def test_runtime_accounting_tolerance_is_fixed(self) -> None:
        evidence = RUNNER.target_runtime_evidence(
            10, 10 + RUNNER.TARGET_RUNTIME_TOLERANCE_NS, 0)
        self.assertTrue(evidence["accepted"])
        self.assertEqual(evidence["tolerance_ns"], 5_000)

    def test_child_runtime_absent_from_target_is_rejected(self) -> None:
        evidence = RUNNER.target_runtime_evidence(
            100, 200, 200 + RUNNER.TARGET_RUNTIME_TOLERANCE_NS)
        self.assertFalse(evidence["accepted"])
        self.assertLess(evidence["signed_difference_ns"], 0)

    def test_target_interrupt_fields_are_rejected(self) -> None:
        before = (0,) * 10
        for index, name in ((5, "irq"), (6, "softirq"), (7, "steal"),
                            (8, "guest"), (9, "guest_nice")):
            with self.subTest(field=name):
                after = list(before)
                after[index] = 1
                evidence = RUNNER.target_interrupt_evidence(
                    before, tuple(after))
                self.assertFalse(evidence["accepted"])
                self.assertEqual(evidence["rejected_fields"][name], 1)

    def test_cleanup_never_signals_after_reap_even_if_pid_reused(self) -> None:
        with mock.patch.object(
                RUNNER, "require_child_identity_gone_after_reap",
                side_effect=RUNNER.EvidenceError("reused")), \
             mock.patch.object(RUNNER.os, "kill") as kill, \
             mock.patch.object(RUNNER.os, "killpg") as killpg:
            failures = RUNNER.cleanup_gated_child(1234, True, True, 0.01)
        self.assertTrue(any("reused" in failure for failure in failures))
        kill.assert_not_called()
        killpg.assert_not_called()

    def test_cleanup_reaps_when_retained_session_scan_fails(self) -> None:
        with mock.patch.object(
                RUNNER, "signal_retained_child_session", return_value=[]), \
             mock.patch.object(
                RUNNER, "waitid_exit_without_reap", return_value=object()), \
             mock.patch.object(
                RUNNER, "process_group_or_session_members",
                side_effect=PermissionError("proc denied")), \
             mock.patch.object(
                RUNNER.os, "wait4", return_value=(1234, 0, object())) as wait4, \
             mock.patch.object(
                RUNNER, "require_child_identity_gone_after_reap"):
            failures = RUNNER.cleanup_gated_child(1234, False, True, 0.01)
        self.assertTrue(any("proc denied" in failure for failure in failures))
        wait4.assert_called_once_with(1234, 0)
        self.assertFalse(any("leader was not reaped" in failure
                             for failure in failures))

    def test_cleanup_kills_descendant_before_releasing_leader(self) -> None:
        leader = {
            "pid": 1234, "state": "Z", "pgrp": 1234, "session": 1234,
            "start_time_ticks": 10,
        }
        descendant = {
            "pid": 1235, "state": "S", "pgrp": 1234, "session": 1234,
            "start_time_ticks": 11,
        }
        with mock.patch.object(
                RUNNER, "signal_retained_child_session",
                return_value=[]) as signal_child, \
             mock.patch.object(
                RUNNER, "waitid_exit_without_reap", return_value=object()), \
             mock.patch.object(
                RUNNER, "process_group_or_session_members",
                side_effect=[[leader, descendant], [leader]]), \
             mock.patch.object(
                RUNNER.os, "wait4", return_value=(1234, 0, object())), \
             mock.patch.object(
                RUNNER, "require_child_identity_gone_after_reap"):
            failures = RUNNER.cleanup_gated_child(1234, False, True, 0.01)
        self.assertEqual(failures, [])
        self.assertEqual(signal_child.call_count, 2)

    def test_waitid_proof_error_reaps_without_resignaling(self) -> None:
        with mock.patch.object(
                RUNNER, "signal_retained_child_session",
                return_value=[]) as signal_child, \
             mock.patch.object(
                RUNNER, "waitid_exit_without_reap",
                side_effect=RUNNER.EvidenceError("bad waitid")), \
             mock.patch.object(
                RUNNER, "wait4_until", return_value=(1234, 0, object())), \
             mock.patch.object(
                RUNNER, "require_child_identity_gone_after_reap"):
            failures = RUNNER.cleanup_gated_child(1234, False, True, 0.01)
        self.assertTrue(any("bad waitid" in failure for failure in failures))
        signal_child.assert_called_once_with(1234, True)

    def test_pending_signal_aborts_after_child_before_next_slot(self) -> None:
        cell = RUNNER.make_matrix()["cells"][0]
        reservation = {"frozen": True}
        with tempfile.TemporaryDirectory(
                prefix="leo2-pending-slot-") as directory, \
             mock.patch.object(
                RUNNER, "reservation_identity", return_value=reservation), \
             mock.patch.object(
                RUNNER, "same_user_pair_affinity", return_value=[]), \
             mock.patch.object(
                RUNNER, "run_gated_benchmark", return_value={}) as gated, \
             mock.patch.object(
                RUNNER.signal, "sigpending", return_value={signal.SIGTERM}):
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "pending control signal"):
                RUNNER.run_invocation(
                    Path("/bin/true"), "baseline", cell, "transform",
                    1, 17, 1.0,
                    Path(directory), 0, 0, 0, Path("/unused"), reservation,
                    "0" * 64, set())
        gated.assert_called_once()

    def test_control_signals_wait_for_exact_child_cleanup(self) -> None:
        source_root = RUNNER_PATH.parents[3]
        cpu, sibling = self.free_physical_pair(source_root)
        PairLease, unused = RUNNER.load_pair_lease(source_root)
        helper = (
            "import importlib.util,signal,sys;"
            "[signal.signal(s,signal.SIG_DFL) for s in "
            "(signal.SIGHUP,signal.SIGINT,signal.SIGTERM)];"
            "p=sys.argv[1];"
            "s=importlib.util.spec_from_file_location('runner_signal_test',p);"
            "m=importlib.util.module_from_spec(s);"
            "sys.modules[s.name]=m;"
            "s.loader.exec_module(m);"
            "m.run_gated_benchmark(['/bin/sleep','10'],int(sys.argv[2]),"
            "int(sys.argv[3]),0.25,m.Path(sys.argv[4]),m.Path(sys.argv[5]))"
        )
        with PairLease(cpu, sibling):
            for control_signal in (
                    signal.SIGTERM, signal.SIGHUP, signal.SIGINT):
                with self.subTest(signal=control_signal), \
                        tempfile.TemporaryDirectory(
                            prefix="leo2-gated-signal-") as directory:
                    root = Path(directory)
                    process = subprocess.Popen([
                        sys.executable, "-c", helper, str(RUNNER_PATH),
                        str(cpu), str(sibling), str(root / "stdout"),
                        str(root / "stderr"),
                    ])
                    child_pid = None
                    deadline = time.monotonic() + 5.0
                    children_path = Path("/proc") / str(process.pid) / \
                        "task" / str(process.pid) / "children"
                    while time.monotonic() < deadline:
                        if process.poll() is not None:
                            break
                        try:
                            children = children_path.read_text().split()
                            if children:
                                candidate = int(children[0])
                                if os.getsid(candidate) == candidate:
                                    child_pid = candidate
                                    break
                        except (FileNotFoundError, ProcessLookupError):
                            pass
                        time.sleep(0.001)
                    self.assertIsNotNone(
                        child_pid, "gated child did not establish its session")
                    os.kill(process.pid, control_signal)
                    return_code = process.wait(timeout=5.0)
                    self.assertEqual(return_code, -control_signal)
                    self.assertFalse(
                        (Path("/proc") / str(child_pid)).exists())
                    self.assertEqual(
                        RUNNER.process_group_or_session_members(child_pid), [])

    def test_campaign_signal_restores_affinity_and_releases_locks(self) -> None:
        source_root = RUNNER_PATH.parents[3]
        cpu, sibling = self.free_physical_pair(source_root)
        helper = """
import importlib.util,json,os,signal,subprocess,sys
from pathlib import Path
[signal.signal(s, signal.SIG_DFL) for s in
 (signal.SIGHUP, signal.SIGINT, signal.SIGTERM)]
runner_path = Path(sys.argv[1])
spec = importlib.util.spec_from_file_location('runner_campaign_signal', runner_path)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
cpu, sibling = int(sys.argv[2]), int(sys.argv[3])
module.DEFAULT_LOCK = Path(sys.argv[4])
status_path = Path(sys.argv[5])
source_root = Path(sys.argv[6])
stdout_path, stderr_path = Path(sys.argv[7]), Path(sys.argv[8])
victim = subprocess.Popen(['/bin/sleep', '10'])
original = sorted(os.sched_getaffinity(victim.pid))
changed = sorted(set(original) - {cpu, sibling})
if not changed or changed == original:
    raise RuntimeError('victim affinity does not exercise restoration')
PairLease, unused = module.load_pair_lease(source_root)
class Options:
    lock_timeout = 5.0
    lock = module.DEFAULT_LOCK
def fake_campaign(options, lock_stream, lock_identity, child_signal_mask):
    guard = PairLease(cpu, sibling)
    guard.__enter__()
    try:
        os.sched_setaffinity(victim.pid, set(changed))
        module.atomic_json(status_path, {
            'stage': 'mask_changed', 'victim_pid': victim.pid,
            'original': original, 'changed': changed})
        return module.run_gated_benchmark(
            ['/bin/sleep', '10'], cpu, sibling, 0.25,
            stdout_path, stderr_path, child_signal_mask)
    finally:
        os.sched_setaffinity(victim.pid, set(original))
        guard.__exit__(None, None, None)
module.run_campaign_locked = fake_campaign
module.run_campaign(Options())
"""
        for control_signal in (signal.SIGTERM, signal.SIGHUP, signal.SIGINT):
            with self.subTest(signal=control_signal), \
                    tempfile.TemporaryDirectory(
                        prefix="leo2-campaign-signal-") as directory:
                root = Path(directory)
                test_lock = root / "campaign.lock"
                status_path = root / "status.json"
                process = subprocess.Popen([
                    sys.executable, "-c", helper, str(RUNNER_PATH),
                    str(cpu), str(sibling), str(test_lock), str(status_path),
                    str(source_root), str(root / "stdout"),
                    str(root / "stderr"),
                ])
                victim_pid = None
                gated_pid = None
                try:
                    deadline = time.monotonic() + 5.0
                    while time.monotonic() < deadline:
                        if process.poll() is not None:
                            break
                        try:
                            status = json.loads(status_path.read_text())
                            if status.get("stage") == "mask_changed":
                                victim_pid = status["victim_pid"]
                                break
                        except (FileNotFoundError, json.JSONDecodeError):
                            pass
                        time.sleep(0.001)
                    self.assertIsNotNone(
                        victim_pid, "campaign did not mutate victim mask")
                    children_path = Path("/proc") / str(process.pid) / \
                        "task" / str(process.pid) / "children"
                    deadline = time.monotonic() + 5.0
                    while time.monotonic() < deadline:
                        if process.poll() is not None:
                            break
                        try:
                            children = children_path.read_text().split()
                        except FileNotFoundError:
                            break
                        for raw_pid in children:
                            candidate = int(raw_pid)
                            if candidate != victim_pid and \
                                    os.getsid(candidate) == candidate:
                                gated_pid = candidate
                                break
                        if gated_pid is not None:
                            break
                        time.sleep(0.001)
                    self.assertIsNotNone(
                        gated_pid, "campaign gated child is unavailable")
                    os.kill(process.pid, control_signal)
                    self.assertEqual(
                        process.wait(timeout=5.0), -control_signal)
                    self.assertEqual(
                        sorted(os.sched_getaffinity(victim_pid)),
                        status["original"])
                    self.assertFalse(
                        (Path("/proc") / str(gated_pid)).exists())
                    self.assertEqual(
                        RUNNER.process_group_or_session_members(gated_pid), [])

                    descriptor = os.open(
                        test_lock,
                        os.O_RDWR | os.O_CLOEXEC | os.O_NOFOLLOW)
                    try:
                        fcntl.flock(
                            descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
                        fcntl.flock(descriptor, fcntl.LOCK_UN)
                    finally:
                        os.close(descriptor)
                    PairLease, unused = RUNNER.load_pair_lease(source_root)
                    with PairLease(cpu, sibling):
                        pass
                finally:
                    if process.poll() is None:
                        process.kill()
                        process.wait(timeout=5.0)
                    for pid in (gated_pid, victim_pid):
                        if pid is None:
                            continue
                        try:
                            os.kill(pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass

    def test_compile_entry_rejects_sole_nonproduction_target(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-nonproduction-entry-") as directory:
            commands = Path(directory) / "compile_commands.json"
            commands.write_text(json.dumps([{
                "directory": directory,
                "command": (
                    "/usr/bin/c++ -o CMakeFiles/"
                    "leopard2_backend_avx2_test_hooks.dir/"
                    "Leopard2BackendAVX2.cpp.o -c /tmp/"
                    "Leopard2BackendAVX2.cpp"),
                "file": "/tmp/Leopard2BackendAVX2.cpp",
                "output": (
                    "CMakeFiles/leopard2_backend_avx2_test_hooks.dir/"
                    "Leopard2BackendAVX2.cpp.o"),
            }]))
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "exactly one production"):
                RUNNER.compile_entry_for(
                    commands, "Leopard2BackendAVX2.cpp")

    def test_compile_profile_rejects_unexpected_isa_enablement(self) -> None:
        common = [
            "/usr/bin/c++", "-Wall", "-Wextra", "-fopenmp", "-O3",
            "-DNDEBUG", "-std=gnu++11", "-falign-functions=64",
            "-mavx2", "-mno-avx512f",
        ]
        RUNNER.validate_compile_profile(
            common, "leopard2_backend_avx2.dir")
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "wrong ISA profile"):
            RUNNER.validate_compile_profile(
                [*common, "-mavx512bw"],
                "leopard2_backend_avx2.dir")
        with self.assertRaisesRegex(
                RUNNER.EvidenceError, "exact Release profile"):
            RUNNER.validate_compile_profile(
                [*common, "-O0"], "leopard2_backend_avx2.dir")

    def test_archive_identity_binds_exact_member_bytes_and_order(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-archive-closure-") as directory:
            root = Path(directory)
            members = []
            for index, name in enumerate(RUNNER.EXPECTED_ARCHIVE_MEMBERS):
                member = root / name
                member.write_bytes(bytes((index + 1, 0xa5, index ^ 0x5a)))
                members.append(member)
            archive = root / "libleopard.a"
            subprocess.run(
                ["/usr/bin/ar", "qc", str(archive),
                 *[str(member) for member in members]], check=True)
            identity = RUNNER.archive_members_identity(archive)
            self.assertEqual(
                identity["members"], list(RUNNER.EXPECTED_ARCHIVE_MEMBERS))
            for index, name in enumerate(RUNNER.EXPECTED_ARCHIVE_MEMBERS):
                payload = bytes((index + 1, 0xa5, index ^ 0x5a))
                self.assertEqual(
                    identity["member_files"][name]["sha256"],
                    RUNNER.sha256_bytes(payload))
            extra = root / "unexpected.o"
            extra.write_bytes(b"unexpected")
            subprocess.run(
                ["/usr/bin/ar", "q", str(archive), str(extra)], check=True)
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "closure or exact order"):
                RUNNER.archive_members_identity(archive)

    def test_historical_tile_spec_schema_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-old-tile-spec-") as directory:
            root = Path(directory)
            manifest = root / "manifest.json"
            manifest.write_text(json.dumps({
                "source_spec": {
                    "schema": "leopard2-direct-tile-lab-spec/v1",
                },
            }))
            options = argparse.Namespace(
                lock_timeout=0.0, progress_seconds=0.0,
                manifest=manifest, lab=root / "missing-lab.py",
                output_dir=root / "output", rerun_failed=False,
                lock=root / "unused.lock")
            with self.assertRaisesRegex(
                    RUNNER.EvidenceError, "frozen tile-screen spec"):
                RUNNER.run_tile_screen_lab(options)


if __name__ == "__main__":
    unittest.main()
