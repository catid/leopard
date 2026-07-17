#!/usr/bin/env python3
"""Portable and mutation tests for the authoritative C7 runner."""

from __future__ import annotations

import argparse
import copy
import fcntl
import gc
import importlib.util
import json
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import unittest
import warnings
from pathlib import Path
from unittest import mock


MODULE_PATH = Path(__file__).with_name("run_authoritative.py")
SPEC = importlib.util.spec_from_file_location("c7_run_authoritative", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
runner = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = runner
SPEC.loader.exec_module(runner)


def load_peer_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


LEOPARD2_EXPERIMENTS = MODULE_PATH.parents[2]
main_runner = load_peer_module(
    "c7_main_compare_peer",
    LEOPARD2_EXPERIMENTS / "main_compare/run_abba.py")
butterfly_runner = load_peer_module(
    "c7_backend_butterfly_peer",
    LEOPARD2_EXPERIMENTS / "backend_butterfly/run_abba.py")


def resign(value: dict) -> dict:
    payload = copy.deepcopy(value)
    payload.pop("digest", None)
    return runner.signed(payload)


class AuthoritativeRunnerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(
            prefix="leopard2-c7-authoritative-test-")
        self.root = Path(self.temporary.name)
        self.manifest, self.raw, self.result = runner.synthetic_bundle(self.root)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def assert_raw_rejected(self, value: dict) -> None:
        with self.assertRaises(runner.EvidenceError):
            runner.validate_raw(resign(value), self.root)

    def test_valid_fixture_and_portable_verify_mode(self) -> None:
        normalized = runner.validate_manifest(self.manifest, self.root)
        self.assertEqual(normalized["cell_count"], len(runner.EXPECTED_CELLS))
        manifest_path = self.root / "manifest.json"
        runner.write_json_exclusive(manifest_path, self.manifest)
        completed = subprocess.run(
            [sys.executable, str(MODULE_PATH), "verify", "--manifest",
             str(manifest_path)], check=False, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1",
                 "PYTHONHASHSEED": "0"})
        self.assertEqual(completed.returncode, 0, completed.stderr)
        report = json.loads(completed.stdout)
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(report["cells"], 12)

    def test_synthetic_raw_and_manifest_are_path_independent(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard2-c7-authoritative-peer-") as directory:
            peer_root = Path(directory)
            peer_manifest, peer_raw, peer_result = runner.synthetic_bundle(peer_root)
            self.assertEqual(runner.canonical_bytes(peer_manifest),
                             runner.canonical_bytes(self.manifest))
            self.assertEqual(runner.canonical_bytes(peer_raw),
                             runner.canonical_bytes(self.raw))
            self.assertEqual(runner.canonical_bytes(peer_result),
                             runner.canonical_bytes(self.result))
            self.assertEqual((peer_root / "snapshot/raw.json").read_bytes(),
                             (self.root / "snapshot/raw.json").read_bytes())

    def test_child_matrix_and_numeric_mutations_rejected(self) -> None:
        mutations = (
            lambda value: value.update(runtime_backend="scalar"),
            lambda value: value.update(affinity=[1]),
            lambda value: value.update(core_git_sha="f" * 40),
            lambda value: value["profile"].update(version=True),
            lambda value: value["correctness"].update(hot_path_allocations=False),
            lambda value: value["benchmarks"].reverse(),
            lambda value: value["benchmarks"][0].update(exact_field=True),
            lambda value: value["benchmarks"][0].update(
                padded_decode_scratch=False),
            lambda value: value["benchmarks"][0][
                "exact_encode_samples_us"].pop(),
            lambda value: value["benchmarks"][0]["exact_encode"].update(
                median_us=123.0),
        )
        inputs = runner.synthetic_inputs()
        for index, mutation in enumerate(mutations):
            changed = copy.deepcopy(self.result)
            mutation(changed)
            with self.subTest(index=index), self.assertRaises(runner.EvidenceError):
                runner.validate_child_result(changed, 0, inputs)

    def test_request_child_and_input_mutations_rejected(self) -> None:
        changed = copy.deepcopy(self.raw)
        changed["request"]["cell_geometry"].reverse()
        self.assert_raw_rejected(changed)

        changed = copy.deepcopy(self.raw)
        changed["request"]["command"][2] = "1"
        changed["child"]["command"][2] = "1"
        self.assert_raw_rejected(changed)

        changed = copy.deepcopy(self.raw)
        changed["child"]["returncode"] = False
        self.assert_raw_rejected(changed)

        changed = copy.deepcopy(self.raw)
        changed["child"]["environment"]["OMP_NUM_THREADS"] = "2"
        self.assert_raw_rejected(changed)

        changed = copy.deepcopy(self.raw)
        changed["inputs_after"]["archive"]["sha256"] = "e" * 64
        changed["inputs_after"]["binding_sha256"] = runner.sha256_bytes(
            runner.canonical_bytes({key: value for key, value in
                                    changed["inputs_after"].items()
                                    if key != "binding_sha256"}))
        self.assert_raw_rejected(changed)

        request = copy.deepcopy(self.raw["request"])
        request["cpu"] = runner.MAX_CPU_ID + 1
        with self.assertRaises(runner.EvidenceError):
            runner.validate_request(request)
        with self.assertRaises(runner.EvidenceError):
            runner.parse_cpu_list(f"0-{runner.MAX_CPU_SET_SIZE}")

        changed = copy.deepcopy(self.raw)
        changed["inputs_after"]["executable"]["sha256"] = "f" * 64
        changed["inputs_after"]["binding_sha256"] = runner.sha256_bytes(
            runner.canonical_bytes({key: value for key, value in
                                    changed["inputs_after"].items()
                                    if key != "binding_sha256"}))
        self.assert_raw_rejected(changed)

    def test_core_source_closure_mutations_rejected(self) -> None:
        for mutation in (
            lambda value: value["core_source_closure"].reverse(),
            lambda value: value["core_source_closure"][0].update(bytes=True),
            lambda value: value["core_source_closure"][0].update(
                path="LeopardFF8.cpp"),
        ):
            inputs = runner.synthetic_inputs()
            mutation(inputs)
            inputs["binding_sha256"] = runner.sha256_bytes(
                runner.canonical_bytes({key: value for key, value in inputs.items()
                                        if key != "binding_sha256"}))
            with self.assertRaises(runner.EvidenceError):
                runner.validate_input_snapshot(inputs)

    def test_v4_build_attestation_rejects_build_and_binary_mutations(self) -> None:
        inputs = runner.synthetic_inputs()
        original = runner.synthetic_build_manifest(inputs)

        def validate(value: dict) -> None:
            runner.derive_build_attestation(
                value, inputs["build_manifest"], inputs["git"]["tooling_commit"],
                inputs["git"]["core_commit"], inputs["source"],
                inputs["build_runner"], inputs["build_validator"],
                inputs["archive"], inputs["executable"],
                inputs["core_source_closure"], inputs["taskset"],
                inputs["python"])

        mutations = (
            lambda value: value.update(schema="leopard2-c7-build-run-manifest/v3"),
            lambda value: value["reproducibility"].update(
                comparison={"status": "not-run"}),
            lambda value: value["builds"][2]["compile_argv"].remove("-O2"),
            lambda value: value["builds"][2].update(sanitizer=True),
            lambda value: value["builds"][2]["executable"].update(
                sha256="f" * 64),
        )
        for index, mutation in enumerate(mutations):
            changed = copy.deepcopy(original)
            mutation(changed)
            with self.subTest(index=index), self.assertRaises(runner.EvidenceError):
                validate(changed)

    def test_isolation_requires_observed_timing_and_duration_binding(self) -> None:
        pair = runner.synthetic_pair_lease(0, 1)
        inactive = runner.isolation_record(
            0, 1, pair, runner.synthetic_cpu_stat(0, user=0, idle=0),
            runner.synthetic_cpu_stat(0, user=0, idle=1),
            runner.synthetic_cpu_stat(1, user=0, idle=0),
            runner.synthetic_cpu_stat(1, user=0, idle=1), 1, 2)
        with self.assertRaises(runner.EvidenceError):
            runner.validate_isolation(inactive, 0, 1)

        changed = copy.deepcopy(self.raw)
        changed["child"]["duration_ns"] += 1
        self.assert_raw_rejected(changed)

    def test_host_isolation_and_reservation_mutations_rejected(self) -> None:
        changed = copy.deepcopy(self.raw)
        changed["host_after"]["timing_cpu"]["frequency_policy"][
            "scaling_governor"] = "powersave"
        self.assert_raw_rejected(changed)

        changed = copy.deepcopy(self.raw)
        isolation = changed["isolation"]
        after_sibling = copy.deepcopy(isolation["after"]["sibling_cpu"])
        after_sibling["fields"]["system"] += 1
        after_sibling["total_jiffies"] += 1
        changed["isolation"] = runner.isolation_record(
            0, 1, isolation["pair_lease"], isolation["before"]["timing_cpu"],
            isolation["after"]["timing_cpu"],
            isolation["before"]["sibling_cpu"], after_sibling, 1, 2)
        self.assert_raw_rejected(changed)

        changed = copy.deepcopy(self.raw)
        changed["reservation"]["payload"]["status"] = "released"
        encoded = runner.canonical_bytes(changed["reservation"]["payload"])
        changed["reservation"]["bytes"] = len(encoded)
        changed["reservation"]["sha256"] = runner.sha256_bytes(encoded)
        self.assert_raw_rejected(changed)

    def test_retained_artifact_path_hash_and_bytes_mutations_rejected(self) -> None:
        result_path = self.root / "snapshot/child/result.json"
        result_path.write_bytes(b"{}\n")
        with self.assertRaises(runner.EvidenceError):
            runner.validate_manifest(self.manifest, self.root)

        result_path.write_bytes(runner.canonical_bytes(self.result) + b"\n")
        for path in ("../result.json", "child//result.json", "./child/result.json"):
            changed = copy.deepcopy(self.raw)
            changed["child"]["result"]["path"] = path
            with self.subTest(path=path):
                self.assert_raw_rejected(changed)

        changed = copy.deepcopy(self.raw)
        changed["child"]["stderr"]["size"] = True
        self.assert_raw_rejected(changed)

        build_path = self.root / "snapshot/provenance/build-run-manifest-v4.json"
        build_original = build_path.read_bytes()
        build_path.write_bytes(b"{}\n")
        with self.assertRaises(runner.EvidenceError):
            runner.validate_manifest(self.manifest, self.root)
        build_path.write_bytes(build_original)

    def test_replay_rejects_symlinks_hardlinks_and_artifact_explosion(self) -> None:
        result_path = self.root / "snapshot/child/result.json"
        saved_result = result_path.read_bytes()
        result_path.unlink()
        result_path.symlink_to("result.json")
        with self.assertRaises(runner.EvidenceError):
            runner.validate_manifest(self.manifest, self.root)
        result_path.unlink()
        result_path.write_bytes(saved_result)

        raw_path = self.root / "snapshot/raw.json"
        saved_raw = raw_path.read_bytes()
        raw_path.unlink()
        raw_path.symlink_to("raw.json")
        with self.assertRaises(runner.EvidenceError):
            runner.validate_manifest(self.manifest, self.root)
        raw_path.unlink()
        raw_path.write_bytes(saved_raw)

        stdout_path = self.root / "snapshot/child/stdout.bin"
        stdout_path.unlink()
        external = self.root.parent / f"{self.root.name}-external-stdout"
        external.write_bytes(b"")
        os.link(external, stdout_path)
        try:
            with self.assertRaises(runner.EvidenceError):
                runner.validate_manifest(self.manifest, self.root)
        finally:
            stdout_path.unlink()
            external.unlink()
            stdout_path.write_bytes(b"")

        crowded = self.root / "crowded"
        crowded.mkdir()
        for index in range(runner.MAX_ARTIFACT_COUNT + 1):
            (crowded / f"{index:04d}").write_bytes(b"")
        with self.assertRaises(runner.EvidenceError):
            runner.artifact_inventory(crowded)

    def test_fifo_inputs_are_rejected_without_blocking(self) -> None:
        fifo = self.root / "manifest-fifo.json"
        os.mkfifo(fifo)
        completed = subprocess.run(
            [sys.executable, str(MODULE_PATH), "verify", "--manifest", str(fifo)],
            check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            timeout=2.0,
            env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1",
                 "PYTHONHASHSEED": "0"})
        self.assertEqual(completed.returncode, 1)
        self.assertIn("not a single-link regular file", completed.stderr)
        self.assertNotIn("Traceback (most recent call last)", completed.stderr)

        result_path = self.root / "snapshot/child/result.json"
        result_path.unlink()
        os.mkfifo(result_path)
        with self.assertRaises(runner.EvidenceError):
            runner.validate_manifest(self.manifest, self.root)

    def test_file_growth_cpu_syntax_numeric_and_cli_bounds(self) -> None:
        for value in ("1,", ",", "1,,2"):
            with self.subTest(cpu_list=value), self.assertRaises(runner.EvidenceError):
                runner.parse_cpu_list(value)

        huge = runner.synthetic_result(0, runner.synthetic_inputs())
        huge["benchmarks"][0]["padded_encode_scratch"] = 10 ** 10000
        with self.assertRaises(runner.EvidenceError):
            runner.validate_child_result(huge, 0, runner.synthetic_inputs())

        growing = self.root / "growing.bin"
        growing.write_bytes(b"a" * 32)
        real_read = os.read
        appended = False

        def racing_read(descriptor: int, count: int) -> bytes:
            nonlocal appended
            data = real_read(descriptor, count)
            if data and not appended:
                appended = True
                with growing.open("ab") as stream:
                    stream.write(b"b")
                    stream.flush()
                    os.fsync(stream.fileno())
            return data

        with mock.patch.object(runner.os, "read", side_effect=racing_read), \
             self.assertRaises(runner.EvidenceError):
            runner.read_bounded(growing)

        completed = subprocess.run(
            [sys.executable, str(MODULE_PATH), "run", "--cpu", "bad"],
            check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(completed.returncode, 1)
        self.assertIn("invalid command line", completed.stderr)

    def test_subprocess_capture_is_file_backed_and_bounded(self) -> None:
        stdout = self.root / "bounded/stdout.bin"
        stderr = self.root / "bounded/stderr.bin"
        completed, timed_out, _started, _ended = runner.run_child_bounded(
            [sys.executable, "-c",
             "import os; os.write(1, b'x' * (2 * 1024 * 1024))"],
            stdout, stderr, 10.0, os.environ)
        self.assertFalse(timed_out)
        self.assertIsInstance(completed.returncode, int)
        self.assertLessEqual(stdout.stat().st_size, runner.MAX_LOG_BYTES)

    def test_timeout_kills_descendants_and_bounds_the_reap(self) -> None:
        marker = self.root / "descendant-marker"
        ready = self.root / "descendant-ready"
        subreaper_before = runner._get_child_subreaper()
        child = (
            "import os,sys,time\n"
            "pid=os.fork()\n"
            "if pid == 0:\n"
            " os.setsid()\n"
            " daemon=os.fork()\n"
            " if daemon != 0: os._exit(0)\n"
            " open(sys.argv[2], 'w').write('setsid-double-fork-ready')\n"
            " time.sleep(1.5)\n"
            " open(sys.argv[1], 'w').write('escaped')\n"
            " os._exit(0)\n"
            "deadline=time.monotonic()+5\n"
            "while not os.path.exists(sys.argv[2]) and time.monotonic()<deadline:\n"
            " time.sleep(.01)\n"
            "time.sleep(10)\n"
        )
        completed, timed_out, _started, _ended = runner.run_child_bounded(
            [sys.executable, "-c", child, str(marker), str(ready)],
            self.root / "group/stdout.bin", self.root / "group/stderr.bin",
            1.0, os.environ)
        self.assertTrue(timed_out)
        self.assertEqual(completed.returncode, 124)
        self.assertEqual(ready.read_text(encoding="utf-8"),
                         "setsid-double-fork-ready")
        self.assertEqual(runner._get_child_subreaper(), subreaper_before)
        time.sleep(1.75)
        self.assertFalse(marker.exists())

        process = mock.Mock(pid=123456)
        process.wait.side_effect = subprocess.TimeoutExpired(["fixture"], 0.1)
        returncode, mock_timed_out = runner._wait_isolated_child(process, 0.1)
        self.assertEqual(returncode, 124)
        self.assertTrue(mock_timed_out)
        self.assertEqual(process.wait.call_args_list[0], mock.call(timeout=0.1))

    def test_post_spawn_attach_failure_reaps_and_restores_both_call_paths(
            self) -> None:
        for index, cwd in enumerate((None, self.root)):
            with self.subTest(cwd=cwd):
                marker = self.root / f"attach-delayed-{index}"
                child = (
                    "import pathlib,sys,time\n"
                    "time.sleep(.5)\n"
                    "pathlib.Path(sys.argv[1]).write_text('escaped')\n"
                )
                subreaper_before = runner._get_child_subreaper()
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always", ResourceWarning)
                    with mock.patch.object(
                            runner.LinuxDescendantContainment, "attach",
                            side_effect=runner.EvidenceError(
                                "synthetic post-spawn attach failure")), \
                         self.assertRaisesRegex(
                             runner.EvidenceError, "post-spawn attach"):
                        if cwd is None:
                            runner.run_child_bounded(
                                [sys.executable, "-c", child, str(marker)],
                                self.root / f"attach-{index}/stdout.bin",
                                self.root / f"attach-{index}/stderr.bin",
                                2.0, os.environ)
                        else:
                            runner.run_capture_bounded(
                                [sys.executable, "-c", child, str(marker)],
                                cwd=cwd, timeout_seconds=2.0,
                                environment=os.environ)
                    gc.collect()
                self.assertFalse(any(
                    item.category is ResourceWarning for item in caught), caught)
                self.assertEqual(
                    runner._get_child_subreaper(), subreaper_before)
                self.assertEqual(
                    runner._raw_direct_child_pids(os.getpid()), set())
                time.sleep(.65)
                self.assertFalse(marker.exists())

    def test_post_spawn_teardown_failures_use_independent_safe_reap(self) -> None:
        real_snapshot = runner._proc_process_snapshot
        real_pidfd_signal = runner._linux_pidfd_signal

        def fail_snapshot_after_preflight():
            fail_snapshot_after_preflight.calls += 1
            if fail_snapshot_after_preflight.calls == 1:
                return real_snapshot()
            raise runner.EvidenceError("synthetic post-spawn procfs failure")

        fail_snapshot_after_preflight.calls = 0

        def fail_pidfd_kill(descriptor: int, signal_number: int) -> None:
            if signal_number == 0:
                real_pidfd_signal(descriptor, signal_number)
                return
            raise runner.EvidenceError("synthetic pidfd signal failure")

        injections = (
            mock.patch.object(
                runner.LinuxDescendantContainment, "_signal_identity",
                side_effect=runner.EvidenceError(
                    "synthetic identity signal failure")),
            mock.patch.object(
                runner, "_linux_pidfd_signal", side_effect=fail_pidfd_kill),
            mock.patch.object(
                runner, "_proc_process_snapshot",
                side_effect=fail_snapshot_after_preflight),
            mock.patch.object(
                runner.LinuxDescendantContainment, "terminate_and_reap",
                side_effect=runner.EvidenceError(
                    "synthetic primary teardown failure")),
        )
        for index, injection in enumerate(injections):
            with self.subTest(injection=index):
                fail_snapshot_after_preflight.calls = 0
                marker = self.root / f"teardown-delayed-{index}"
                ready = self.root / f"teardown-ready-{index}"
                child = (
                    "import os,pathlib,sys,time\n"
                    "pid=os.fork()\n"
                    "if pid == 0:\n"
                    " os.setsid()\n"
                    " daemon=os.fork()\n"
                    " if daemon != 0: os._exit(0)\n"
                    " pathlib.Path(sys.argv[2]).write_text('ready')\n"
                    " time.sleep(.6)\n"
                    " pathlib.Path(sys.argv[1]).write_text('escaped')\n"
                    " os._exit(0)\n"
                    "deadline=time.monotonic()+2\n"
                    "while not os.path.exists(sys.argv[2]) and "
                    "time.monotonic()<deadline: time.sleep(.005)\n"
                    "time.sleep(10)\n"
                )
                subreaper_before = runner._get_child_subreaper()
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always", ResourceWarning)
                    with injection, self.assertRaisesRegex(
                            runner.EvidenceError,
                            "safe emergency cleanup"):
                        runner.run_child_bounded(
                            [sys.executable, "-c", child,
                             str(marker), str(ready)],
                            self.root / f"teardown-{index}/stdout.bin",
                            self.root / f"teardown-{index}/stderr.bin",
                            .3, os.environ)
                    gc.collect()
                self.assertEqual(ready.read_text(encoding="utf-8"), "ready")
                self.assertFalse(any(
                    item.category is ResourceWarning for item in caught), caught)
                self.assertEqual(
                    runner._get_child_subreaper(), subreaper_before)
                self.assertEqual(
                    runner._raw_direct_child_pids(os.getpid()), set())
                time.sleep(.75)
                self.assertFalse(marker.exists())

    def test_emergency_reap_does_not_signal_disappeared_or_unrelated_pid(
            self) -> None:
        with mock.patch.object(
                runner, "_raw_direct_child_pids",
                side_effect=({456789}, set())) as children, \
             mock.patch.object(
                 runner.os, "pidfd_open", side_effect=OSError("gone")), \
             mock.patch.object(runner.os, "kill") as kill:
            runner.LinuxDescendantContainment._kill_unreaped_direct_child(
                456789, os.getpid())
        self.assertEqual(children.call_count, 2)
        kill.assert_not_called()

        with mock.patch.object(
                runner, "_raw_direct_child_pids", return_value={654321}), \
             mock.patch.object(runner.os, "pidfd_open") as pidfd_open, \
             mock.patch.object(runner.os, "kill") as kill:
            runner.LinuxDescendantContainment._kill_unreaped_direct_child(
                456789, os.getpid())
        pidfd_open.assert_not_called()
        kill.assert_not_called()

    def test_post_spawn_failure_leaves_unrelated_same_uid_process_alive(
            self) -> None:
        pid_path = self.root / "unrelated.pid"
        marker = self.root / "unrelated-alive"
        helper = (
            "import os,pathlib,sys,time\n"
            "pid=os.fork()\n"
            "if pid != 0:\n"
            " pathlib.Path(sys.argv[1]).write_text(str(pid))\n"
            " os._exit(0)\n"
            "os.setsid()\n"
            "time.sleep(.7)\n"
            "pathlib.Path(sys.argv[2]).write_text('alive')\n"
            "time.sleep(10)\n"
        )
        launcher = subprocess.Popen(
            [sys.executable, "-c", helper, str(pid_path), str(marker)],
            stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL, close_fds=True)
        self.assertEqual(launcher.wait(timeout=2.0), 0)
        unrelated_pid = int(pid_path.read_text(encoding="utf-8"))
        unrelated_pidfd = os.pidfd_open(unrelated_pid, 0)
        try:
            self.assertEqual(
                runner._raw_direct_child_pids(os.getpid()), set())
            delayed = self.root / "contained-delayed"
            child = (
                "import pathlib,sys,time\n"
                "time.sleep(.5)\n"
                "pathlib.Path(sys.argv[1]).write_text('escaped')\n"
            )
            with mock.patch.object(
                    runner.LinuxDescendantContainment, "attach",
                    side_effect=runner.EvidenceError(
                        "synthetic post-spawn attach failure")), \
                 self.assertRaises(runner.EvidenceError):
                runner.run_child_bounded(
                    [sys.executable, "-c", child, str(delayed)],
                    self.root / "unrelated-check/stdout.bin",
                    self.root / "unrelated-check/stderr.bin",
                    2.0, os.environ)
            time.sleep(.85)
            self.assertEqual(marker.read_text(encoding="utf-8"), "alive")
            self.assertFalse(delayed.exists())
            signal.pidfd_send_signal(unrelated_pidfd, 0)
        finally:
            try:
                signal.pidfd_send_signal(unrelated_pidfd, signal.SIGKILL)
            except ProcessLookupError:
                pass
            os.close(unrelated_pidfd)

    def test_descendant_containment_fails_closed_before_spawn_when_unavailable(
            self) -> None:
        with mock.patch.object(runner.sys, "platform", "not-linux"), \
             mock.patch.object(runner.subprocess, "Popen") as popen, \
             self.assertRaisesRegex(runner.EvidenceError, "requires Linux"):
            runner.run_child_bounded(
                [sys.executable, "-c", "pass"],
                self.root / "unavailable/stdout.bin",
                self.root / "unavailable/stderr.bin", 1.0, os.environ)
        popen.assert_not_called()

        for index, unavailable in enumerate((
                "_validate_linux_pidfd_support", "_get_child_subreaper",
                "_proc_process_snapshot")):
            with self.subTest(unavailable=unavailable), \
                 mock.patch.object(
                     runner, unavailable,
                     side_effect=runner.EvidenceError(
                         f"{unavailable} unavailable")), \
                 mock.patch.object(runner.subprocess, "Popen") as popen, \
                 self.assertRaisesRegex(runner.EvidenceError, "unavailable"):
                runner.run_child_bounded(
                    [sys.executable, "-c", "pass"],
                    self.root / f"unavailable-{index}/stdout.bin",
                    self.root / f"unavailable-{index}/stderr.bin",
                    1.0, os.environ)
            popen.assert_not_called()

    def test_live_build_validator_invocation_is_fail_closed(self) -> None:
        validator = self.root / runner.BUILD_VALIDATOR_RELATIVE
        validator.parent.mkdir(parents=True)
        validator.write_text("fixture\n", encoding="utf-8")
        build_manifest = self.root / "build-manifest-v4.json"
        build_manifest.write_text("{}\n", encoding="utf-8")
        successful = subprocess.CompletedProcess(
            ["validator"], 0, b"C7 evidence validation passed (live)\n", b"")
        with mock.patch.object(
                runner, "run_capture_bounded",
                return_value=(successful,
                              b"C7 evidence validation passed (live)\n",
                              b"", False)) as call:
            runner.run_build_validator(self.root, build_manifest)
        argv = call.call_args.args[0]
        self.assertIn("--live", argv)
        self.assertIn("--require-checkout-head", argv)
        self.assertEqual(call.call_args.kwargs["environment"],
                         runner.VALIDATOR_ENVIRONMENT)

        changed = subprocess.CompletedProcess(["validator"], 0, b"portable only\n", b"")
        with mock.patch.object(
                runner, "run_capture_bounded",
                return_value=(changed, b"portable only\n", b"", False)), \
             self.assertRaises(runner.EvidenceError):
            runner.run_build_validator(self.root, build_manifest)

    def test_verified_build_attestation_binds_exact_local_avx2_files(self) -> None:
        source = self.root / runner.SOURCE_RELATIVE
        build_runner = self.root / runner.BUILD_RUNNER_RELATIVE
        build_validator = self.root / runner.BUILD_VALIDATOR_RELATIVE
        archive = self.root / ".research/fixture/build/core-avx2/liblibleopard.a"
        executable = self.root / ".research/fixture/build/c7-avx2"
        for path, data in ((source, b"source\n"), (build_runner, b"runner\n"),
                           (build_validator, b"validator\n"),
                           (archive, b"!<arch>\n"), (executable, b"ELF fixture\n")):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        executable.chmod(0o755)
        identities = {
            "git": {"tooling_commit": "b" * 40, "core_commit": "e" * 40},
            "source": {**runner.file_identity(source, "source"),
                       "path": runner.SOURCE_RELATIVE.as_posix()},
            "build_runner": {**runner.file_identity(build_runner, "build-runner"),
                             "path": runner.BUILD_RUNNER_RELATIVE.as_posix()},
            "build_validator": {
                **runner.file_identity(build_validator, "build-validator"),
                "path": runner.BUILD_VALIDATOR_RELATIVE.as_posix()},
            "archive": runner.file_identity(archive, "archive"),
            "executable": runner.file_identity(executable, "executable"),
        }
        fixture_inputs = runner.synthetic_inputs()
        manifest_inputs = copy.deepcopy(fixture_inputs)
        manifest_inputs.update({
            "git": identities["git"],
            "source": {"bytes": identities["source"]["size"],
                       "path": runner.SOURCE_RELATIVE.as_posix(),
                       "sha256": identities["source"]["sha256"]},
            "build_runner": {"bytes": identities["build_runner"]["size"],
                             "path": runner.BUILD_RUNNER_RELATIVE.as_posix(),
                             "sha256": identities["build_runner"]["sha256"]},
            "build_validator": {
                "bytes": identities["build_validator"]["size"],
                "path": runner.BUILD_VALIDATOR_RELATIVE.as_posix(),
                "sha256": identities["build_validator"]["sha256"]},
            "archive": {"bytes": identities["archive"]["size"],
                        "path": ".research/fixture/build/core-avx2/liblibleopard.a",
                        "sha256": identities["archive"]["sha256"]},
            "executable": {"bytes": identities["executable"]["size"],
                           "path": ".research/fixture/build/c7-avx2",
                           "sha256": identities["executable"]["sha256"]},
        })
        build_manifest = self.root / ".research/fixture/results/build-run-manifest.json"
        build_manifest.parent.mkdir(parents=True)
        build_manifest.write_bytes(runner.canonical_pretty_bytes(
            runner.synthetic_build_manifest(manifest_inputs)))
        with mock.patch.object(runner, "run_build_validator"):
            identity, attestation = runner.verified_build_attestation(
                self.root, build_manifest, "b" * 40, "e" * 40,
                identities["source"], identities["build_runner"],
                identities["build_validator"], identities["archive"],
                identities["executable"], manifest_inputs["core_source_closure"],
                manifest_inputs["taskset"], manifest_inputs["python"],
                archive, executable)
        self.assertEqual(attestation["manifest"]["sha256"], identity["sha256"])
        self.assertEqual(attestation["avx2"]["optimization"], "-O2")
        self.assertEqual(identities["build_runner"]["mode"] & 0o111, 0)
        self.assertEqual(identities["build_validator"]["mode"] & 0o111, 0)

        with executable.open("ab") as stream:
            stream.write(b"trailing mutation")
        with mock.patch.object(runner, "run_build_validator"), \
             self.assertRaises(runner.EvidenceError):
            runner.verified_build_attestation(
                self.root, build_manifest, "b" * 40, "e" * 40,
                identities["source"], identities["build_runner"],
                identities["build_validator"], identities["archive"],
                runner.file_identity(executable, "executable"),
                manifest_inputs["core_source_closure"],
                manifest_inputs["taskset"], manifest_inputs["python"],
                archive, executable)

    def test_committed_python_tool_modes_are_nonexecutable_and_enforced(self) -> None:
        source_root = MODULE_PATH.parents[4]
        head = subprocess.check_output(
            ["git", "-C", str(source_root), "rev-parse", "HEAD"],
            text=True).strip()
        for relative, kind in (
                (runner.BUILD_RUNNER_RELATIVE, "build-runner"),
                (runner.BUILD_VALIDATOR_RELATIVE, "build-validator")):
            entry = subprocess.check_output(
                ["git", "-C", str(source_root), "ls-tree", head, "--",
                 relative.as_posix()], text=True)
            self.assertTrue(entry.startswith("100644 blob "), entry)
            identity = runner.committed_file_identity(
                source_root, head, relative, kind)
            self.assertEqual(identity["mode"] & 0o111, 0)

        repository = self.root / "mode-repository"
        repository.mkdir()
        subprocess.run(["git", "init", "-q"], cwd=repository, check=True)
        subprocess.run(["git", "config", "user.email", "fixture@example.com"],
                       cwd=repository, check=True)
        subprocess.run(["git", "config", "user.name", "Fixture"],
                       cwd=repository, check=True)
        tool = repository / "tool.py"
        tool.write_text("print('fixture')\n", encoding="utf-8")
        tool.chmod(0o644)
        subprocess.run(["git", "add", "tool.py"], cwd=repository, check=True)
        subprocess.run(["git", "commit", "-qm", "fixture"],
                       cwd=repository, check=True)
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=repository, text=True).strip()
        identity = runner.committed_file_identity(
            repository, commit, Path("tool.py"), "build-runner")
        self.assertEqual(identity["mode"] & 0o111, 0)
        tool.chmod(0o755)
        with self.assertRaises(runner.EvidenceError):
            runner.committed_file_identity(
                repository, commit, Path("tool.py"), "build-runner")

        inputs = runner.synthetic_inputs()
        self.assertNotIn("mode", inputs["build_runner"])
        self.assertNotIn("mode", inputs["build_validator"])
        runner.validate_input_snapshot(inputs)
        inputs["build_runner"]["mode"] = 0o755
        inputs["binding_sha256"] = runner.sha256_bytes(runner.canonical_bytes(
            {key: value for key, value in inputs.items()
             if key != "binding_sha256"}))
        with self.assertRaises(runner.EvidenceError):
            runner.validate_input_snapshot(inputs)

    def test_manifest_and_raw_coordinated_mutations_rejected(self) -> None:
        changed = copy.deepcopy(self.manifest)
        changed["validated_output"]["cell_count"] = 11
        with self.assertRaises(runner.EvidenceError):
            runner.validate_manifest(resign(changed), self.root)

        raw_path = self.root / "raw.json"
        pretty = (json.dumps(self.raw, indent=2, sort_keys=True) + "\n").encode()
        raw_path.write_bytes(pretty)
        changed = copy.deepcopy(self.manifest)
        changed["raw"] = {"path": "raw.json", "size": len(pretty),
                          "sha256": runner.sha256_bytes(pretty)}
        with self.assertRaises(runner.EvidenceError):
            runner.validate_manifest(resign(changed), self.root)

    def test_verify_rejects_noncanonical_manifest_bytes(self) -> None:
        manifest_path = self.root / "manifest.json"
        manifest_path.write_text(json.dumps(self.manifest, indent=2),
                                 encoding="utf-8")
        options = argparse.Namespace(manifest=manifest_path)
        with self.assertRaises(runner.EvidenceError):
            runner.verify_campaign(options)

    def test_reservation_lock_canonicality_and_current_bytes(self) -> None:
        payload = {"benchmark_cpu": 0, "nonce": "fixture-nonce",
                   "owner": "unit-test", "reserved_sibling": 1,
                   "schema": runner.RESERVATION_SCHEMA, "status": "held"}
        path = self.root / "reservation.json"
        path.write_bytes(runner.canonical_bytes(payload))
        path.chmod(0o600)
        guard = runner.Reservation(path, 0, 1)
        with guard:
            with self.assertRaises(runner.EvidenceError):
                with runner.Reservation(path, 0, 1):
                    pass
            path.write_bytes(runner.canonical_bytes({**payload, "nonce": "changed-nonce"}))
            with self.assertRaises(runner.EvidenceError):
                guard.validate_current()

        path.write_bytes(runner.canonical_bytes(payload) + b"\n")
        with self.assertRaises(runner.EvidenceError):
            with runner.Reservation(path, 0, 1):
                pass
        path.write_bytes(runner.canonical_bytes(payload))
        with runner.Reservation(path, 0, 1):
            pass

    def test_reservation_rejects_path_replacement(self) -> None:
        payload = {"benchmark_cpu": 0, "nonce": "fixture-nonce",
                   "owner": "unit-test", "reserved_sibling": 1,
                   "schema": runner.RESERVATION_SCHEMA, "status": "held"}
        encoded = runner.canonical_bytes(payload)
        path = self.root / "replacement-reservation.json"
        path.write_bytes(encoded)
        path.chmod(0o600)
        guard = runner.Reservation(path, 0, 1)
        with guard as retained:
            path.rename(self.root / "old-reservation.json")
            path.write_bytes(encoded)
            path.chmod(0o600)
            with self.assertRaises(runner.EvidenceError):
                guard.validate_current()
            with self.assertRaises(runner.EvidenceError):
                with runner.Reservation(path, 0, 1):
                    pass
            changed = copy.deepcopy(retained)
            changed["payload"]["owner"] = "forged-owner"
            changed_raw = runner.canonical_bytes(changed["payload"])
            changed["bytes"] = len(changed_raw)
            changed["sha256"] = runner.sha256_bytes(changed_raw)
            with self.assertRaises(runner.EvidenceError):
                guard.validate_current(changed)

    def test_pair_lease_serializes_pair_and_rejects_path_replacement(self) -> None:
        runtime_root = self.root / "runtime-root"
        runtime_root.mkdir(mode=0o700)
        runtime_root.chmod(0o700)
        first = runner.PairLease(7, 8, runtime_root=runtime_root)
        with first as identity:
            self.assertEqual(set(identity), {
                "schema", "authority", "device", "directory_device",
                "directory_inode", "inode", "lock", "path", "payload",
                "sha256"})
            with self.assertRaises(runner.EvidenceError):
                with runner.PairLease(7, 8, runtime_root=runtime_root):
                    pass
            old = first.path.with_suffix(".old")
            first.path.rename(old)
            first.path.write_bytes(runner.canonical_bytes(
                runner.pair_lease_payload(7, 8)))
            first.path.chmod(0o600)
            with self.assertRaises(runner.EvidenceError):
                first.validate_current()

    def test_pair_directory_replacement_cannot_split_kernel_authority(self) -> None:
        runtime_root = self.root / "runtime-directory-replacement"
        runtime_root.mkdir(mode=0o700)
        runtime_root.chmod(0o700)
        first = runner.PairLease(9, 10, runtime_root=runtime_root)
        with first as retained:
            original = runtime_root / "leopard2-cpu-leases"
            replaced = runtime_root / "leopard2-cpu-leases.old"
            original.rename(replaced)
            original.mkdir(mode=0o700)
            original.chmod(0o700)
            with self.assertRaises(runner.EvidenceError):
                with runner.PairLease(9, 10, runtime_root=runtime_root):
                    pass
            changed = copy.deepcopy(retained)
            changed["payload"]["uid"] += 1
            changed["path"] = str(runner.pair_lease_directory(
                changed["payload"]["uid"]) / runner.pair_lease_name(
                    9, 10, changed["payload"]["uid"]))
            changed["sha256"] = runner.sha256_bytes(
                runner.canonical_bytes(changed["payload"]))
            with self.assertRaises(runner.EvidenceError):
                first.validate_current(changed)

    def test_shared_anchor_blocks_cross_runner_inode_splits(self) -> None:
        runtime_root = self.root / "shared-anchor-runtime"
        runtime_root.mkdir(mode=0o700)
        runtime_root.chmod(0o700)
        pair = runner.PairLease(21, 22, runtime_root=runtime_root)
        with pair as retained:
            current = runtime_root / "leopard2-cpu-leases"
            original = runtime_root / "leopard2-cpu-leases.original"
            current.rename(original)
            current.mkdir(mode=0o700)
            peer = os.open(
                runtime_root, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) |
                getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0))
            try:
                with self.assertRaises(BlockingIOError):
                    fcntl.flock(peer, fcntl.LOCK_EX | fcntl.LOCK_NB)
            finally:
                os.close(peer)
            current.rmdir()
            original.rename(current)
            pair.validate_current(retained)

        anchor_root = self.root / "shared-reservation-anchor"
        anchor_root.mkdir(mode=0o700)
        anchor_root.chmod(0o700)
        reservation = self.root / "shared-reservation.json"
        payload = {
            "benchmark_cpu": 23, "nonce": "shared-anchor-fixture",
            "owner": "unit-test", "reserved_sibling": 24,
            "schema": runner.RESERVATION_SCHEMA, "status": "held",
        }
        encoded = runner.canonical_bytes(payload)
        reservation.write_bytes(encoded)
        reservation.chmod(0o600)
        anchor = runner.StableLeaseAnchor(anchor_root)
        with anchor:
            guard = runner.Reservation(reservation, 23, 24, anchor=anchor)
            with guard as retained:
                original = self.root / "shared-reservation.original"
                reservation.rename(original)
                reservation.write_bytes(encoded)
                reservation.chmod(0o600)
                peer = os.open(
                    anchor_root, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) |
                    getattr(os, "O_CLOEXEC", 0) |
                    getattr(os, "O_NOFOLLOW", 0))
                try:
                    with self.assertRaises(BlockingIOError):
                        fcntl.flock(peer, fcntl.LOCK_EX | fcntl.LOCK_NB)
                finally:
                    os.close(peer)
                reservation.unlink()
                original.rename(reservation)
                guard.validate_current(retained)

    def test_current_runners_share_anchor_across_reservation_replacement(
            self) -> None:
        runtime_root = self.root / "cross-runner-runtime"
        runtime_root.mkdir(mode=0o700)
        runtime_root.chmod(0o700)
        reservation_parent = self.root / "cross-runner-reservation-directory"
        reservation_parent.mkdir(mode=0o700)
        reservation_parent.chmod(0o700)
        reservation = reservation_parent / "reservation.json"
        original = reservation_parent / "reservation.original"
        original_parent = self.root / "cross-runner-reservation-directory.original"
        payload = {
            "benchmark_cpu": 25,
            "nonce": "cross-runner-fixture-nonce",
            "owner": "unit-test",
            "reserved_sibling": 26,
            "schema": runner.RESERVATION_SCHEMA,
            "status": "held",
        }
        encoded = runner.canonical_bytes(payload)
        reservation.write_bytes(encoded)
        reservation.chmod(0o600)

        def replace_while_held(blocked_error: type[BaseException], acquire) -> None:
            reservation.rename(original)
            reservation.write_bytes(encoded)
            reservation.chmod(0o600)
            try:
                with self.assertRaises(blocked_error):
                    acquire()
            finally:
                reservation.unlink()
                original.rename(reservation)

        def replace_directory_while_held(
                blocked_error: type[BaseException], acquire) -> None:
            reservation_parent.rename(original_parent)
            reservation_parent.mkdir(mode=0o700)
            reservation_parent.chmod(0o700)
            reservation.write_bytes(encoded)
            reservation.chmod(0o600)
            try:
                with self.assertRaises(blocked_error):
                    acquire()
            finally:
                reservation.unlink()
                reservation_parent.rmdir()
                original_parent.rename(reservation_parent)

        def replace_both(blocked_error: type[BaseException], acquire) -> None:
            replace_while_held(blocked_error, acquire)
            replace_directory_while_held(blocked_error, acquire)

        def acquire_c7() -> None:
            with runner.Reservation(
                    reservation, 25, 26, runtime_root=runtime_root):
                pass

        def acquire_main() -> None:
            with main_runner.Reservation(
                    reservation, 25, 26, runtime_root=runtime_root):
                pass

        def acquire_butterfly() -> None:
            handle, _identity = butterfly_runner.reservation_record(
                reservation, 25, 26, runtime_root=runtime_root)
            handle.close()

        with runner.Reservation(
                reservation, 25, 26, runtime_root=runtime_root):
            replace_both(main_runner.EvidenceError, acquire_main)
        with main_runner.Reservation(
                reservation, 25, 26, runtime_root=runtime_root):
            replace_both(runner.EvidenceError, acquire_c7)
        with runner.Reservation(
                reservation, 25, 26, runtime_root=runtime_root):
            replace_both(
                butterfly_runner.EvidenceError, acquire_butterfly)
        handle, _identity = butterfly_runner.reservation_record(
            reservation, 25, 26, runtime_root=runtime_root)
        try:
            replace_both(runner.EvidenceError, acquire_c7)
        finally:
            handle.close()
        with main_runner.Reservation(
                reservation, 25, 26, runtime_root=runtime_root):
            replace_both(
                butterfly_runner.EvidenceError, acquire_butterfly)
        handle, _identity = butterfly_runner.reservation_record(
            reservation, 25, 26, runtime_root=runtime_root)
        try:
            replace_both(main_runner.EvidenceError, acquire_main)
        finally:
            handle.close()

    def test_current_c7_retains_predecessor_kernel_lease_interop(self) -> None:
        runtime_root = self.root / "legacy-kernel-runtime"
        runtime_root.mkdir(mode=0o700)
        runtime_root.chmod(0o700)
        legacy_pair = runner.KernelNamespaceLease(
            "pair", runner.pair_lease_payload(27, 28))
        legacy_pair.acquire()
        try:
            with self.assertRaises(runner.EvidenceError):
                with runner.PairLease(27, 28, runtime_root=runtime_root):
                    pass
        finally:
            legacy_pair.release()
        with runner.PairLease(27, 28, runtime_root=runtime_root):
            legacy_peer = runner.KernelNamespaceLease(
                "pair", runner.pair_lease_payload(27, 28))
            with self.assertRaises(runner.EvidenceError):
                legacy_peer.acquire()

        reservation = self.root / "legacy-kernel-reservation.json"
        payload = {
            "benchmark_cpu": 27,
            "nonce": "legacy-kernel-fixture",
            "owner": "unit-test",
            "reserved_sibling": 28,
            "schema": runner.RESERVATION_SCHEMA,
            "status": "held",
        }
        reservation.write_bytes(runner.canonical_bytes(payload))
        reservation.chmod(0o600)
        authority = runner.reservation_authority_payload(
            reservation, 27, 28)
        legacy_reservation = runner.KernelNamespaceLease(
            "reservation", authority)
        legacy_reservation.acquire()
        try:
            with self.assertRaises(runner.EvidenceError):
                with runner.Reservation(
                        reservation, 27, 28, runtime_root=runtime_root):
                    pass
        finally:
            legacy_reservation.release()
        with runner.Reservation(
                reservation, 27, 28, runtime_root=runtime_root):
            legacy_peer = runner.KernelNamespaceLease(
                "reservation", authority)
            with self.assertRaises(runner.EvidenceError):
                legacy_peer.acquire()

    def test_git_identity_uses_exact_tree_and_rejects_untracked_files(self) -> None:
        repository = self.root / "repository"
        repository.mkdir()
        subprocess.run(["git", "init", "-q"], cwd=repository, check=True)
        subprocess.run(["git", "config", "user.email", "fixture@example.com"],
                       cwd=repository, check=True)
        subprocess.run(["git", "config", "user.name", "Fixture"],
                       cwd=repository, check=True)
        (repository / "tracked.txt").write_text("fixture\n", encoding="utf-8")
        subprocess.run(["git", "add", "tracked.txt"], cwd=repository, check=True)
        subprocess.run(["git", "commit", "-qm", "fixture"],
                       cwd=repository, check=True)
        core_commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=repository, text=True).strip()
        core_tree = subprocess.check_output(
            ["git", "rev-parse", "HEAD^{tree}"], cwd=repository, text=True).strip()
        (repository / "tooling.txt").write_text("tooling\n", encoding="utf-8")
        subprocess.run(["git", "add", "tooling.txt"], cwd=repository, check=True)
        subprocess.run(["git", "commit", "-qm", "tooling"],
                       cwd=repository, check=True)
        tooling_commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=repository, text=True).strip()
        identity = runner.git_identity(repository, tooling_commit, core_commit)
        self.assertEqual(identity["core_tree"], core_tree)
        self.assertTrue(identity["core_is_ancestor"])
        with self.assertRaises(runner.EvidenceError):
            runner.git_identity(repository, tooling_commit, "f" * 40)
        (repository / "untracked.txt").write_text("untracked\n", encoding="utf-8")
        with self.assertRaises(runner.EvidenceError):
            runner.git_identity(repository, tooling_commit, core_commit)

    def test_preflight_failure_is_retained_without_running_child(self) -> None:
        executable = self.root / "executable"
        archive = self.root / "archive.a"
        executable.write_bytes(b"fixture")
        archive.write_bytes(b"!<arch>\n")
        build_manifest = self.root / "build-manifest.json"
        build_manifest.write_bytes(b"{}\n")
        executable.chmod(0o755)
        output = self.root / "failed-run"
        options = argparse.Namespace(
            output=output, executable=executable, archive=archive,
            source_root=self.root, cpu=0, sibling=1, timeout=1.0,
            expected_tooling_commit="a" * 40,
            expected_core_commit="a" * 40,
            reservation_file=self.root / "missing-reservation.json",
            build_manifest=build_manifest)
        with mock.patch.object(
                runner, "validate_topology",
                side_effect=runner.EvidenceError("synthetic topology failure")):
            with self.assertRaises(runner.EvidenceError):
                runner.run_campaign(options)
        failure_path = output / "failure.json"
        self.assertTrue(failure_path.is_file())
        failure = runner.strict_json(failure_path.read_bytes(), "failure")
        runner.validate_failure(failure, output, check_files=True)
        self.assertEqual(failure["schema"], runner.FAILURE_SCHEMA)
        self.assertEqual(failure["status"], "failed")
        self.assertEqual(failure["failure_code"], "host-capture-failed")
        self.assertEqual(failure["error"],
                         runner.FAILURE_DIAGNOSTICS["host-capture-failed"])
        self.assertIsNone(failure["child"])
        self.assertIsNotNone(failure["request"])
        self.assertFalse((output / "manifest.json").exists())

    def _run_child_failure(self, *, timeout: bool) -> tuple[dict, Path]:
        inputs = runner.synthetic_inputs()
        executable = self.root / ("timeout-executable" if timeout else "exit-executable")
        archive = self.root / ("timeout-archive.a" if timeout else "exit-archive.a")
        build_manifest = self.root / (
            "timeout-build-manifest.json" if timeout else "exit-build-manifest.json")
        reservation = self.root / (
            "timeout-reservation.json" if timeout else "exit-reservation.json")
        output = self.root / ("timeout-run" if timeout else "exit-run")
        executable.write_bytes(b"fixture")
        executable.chmod(0o755)
        archive.write_bytes(b"!<arch>\n")
        build_manifest.write_bytes(runner.canonical_pretty_bytes(
            runner.synthetic_build_manifest(inputs)))
        reservation_payload = {
            "benchmark_cpu": 0, "nonce": "failure-fixture", "owner": "unit-test",
            "reserved_sibling": 1, "schema": runner.RESERVATION_SCHEMA,
            "status": "held",
        }
        reservation.write_bytes(runner.canonical_bytes(reservation_payload))
        reservation.chmod(0o600)
        options = argparse.Namespace(
            output=output, executable=executable, archive=archive,
            build_manifest=build_manifest, source_root=self.root, cpu=0, sibling=1,
            timeout=1.0, expected_tooling_commit="b" * 40,
            expected_core_commit="e" * 40, reservation_file=reservation)
        pair_guard = mock.MagicMock()
        pair_identity = runner.synthetic_pair_lease(0, 1)
        pair_guard.__enter__.return_value = pair_identity
        pair_guard.identity = copy.deepcopy(pair_identity)
        pair_guard.__exit__.return_value = None
        snapshots = [
            runner.synthetic_cpu_stat(0, user=0, idle=0),
            runner.synthetic_cpu_stat(1, user=0, idle=0),
            runner.synthetic_cpu_stat(0, user=1, idle=0),
            runner.synthetic_cpu_stat(1, user=0, idle=1),
        ]

        def child_effect(command: object, stdout_path: Path, stderr_path: Path,
                         *unused_args: object, **unused_kwargs: object) -> object:
            runner.write_exclusive(stdout_path, b"partial" if timeout else b"")
            runner.write_exclusive(stderr_path, b"timed out" if timeout else b"failed\n")
            return (subprocess.CompletedProcess(
                ["fixture"], 124 if timeout else 7, b"", b""), timeout, 1, 2)
        with mock.patch.object(runner, "validate_topology", return_value=({0, 1, 2}, {2})), \
             mock.patch.object(runner, "host_identity",
                               return_value=runner.synthetic_host(0, 1)), \
             mock.patch.object(runner, "PairLease", return_value=pair_guard), \
             mock.patch.object(runner, "input_snapshot", return_value=inputs), \
             mock.patch.object(runner, "cpu_stat_snapshot", side_effect=snapshots), \
             mock.patch.object(runner.os, "sched_getaffinity", return_value={0, 1, 2}), \
             mock.patch.object(runner.os, "sched_setaffinity"), \
             mock.patch.object(runner, "run_child_bounded", side_effect=child_effect):
            with self.assertRaises(runner.EvidenceError):
                runner.run_campaign(options)
        failure = runner.strict_json((output / "failure.json").read_bytes(), "failure")
        runner.validate_failure(failure, output, check_files=True)
        self.assertIsNotNone(failure["request"])
        self.assertIsNotNone(failure["inputs_before"])
        self.assertIsNotNone(failure["host_before"])
        self.assertIsNotNone(failure["reservation"])
        self.assertIsNotNone(failure["build_provenance"])
        self.assertIsNotNone(failure["isolation"])
        self.assertFalse((output / "manifest.json").exists())
        stdout = runner.read_artifact(output, failure["child"]["stdout"], "stdout")
        stderr = runner.read_artifact(output, failure["child"]["stderr"], "stderr")
        self.assertEqual(stdout, b"partial" if timeout else b"")
        self.assertEqual(stderr, b"timed out" if timeout else b"failed\n")
        return failure, output

    def _run_success_variant(self, name: str, *, mutate_on_final_lease: bool = False,
                             exit_failure: bool = False,
                             restore_failure: bool = False,
                             malformed_result: bool = False,
                             post_child_guard_failure: bool = False) -> tuple[Path, dict | None]:
        inputs = runner.synthetic_inputs()
        executable = self.root / f"{name}-executable"
        archive = self.root / f"{name}-archive.a"
        build_manifest = self.root / f"{name}-build-manifest.json"
        reservation = self.root / f"{name}-reservation.json"
        output = self.root / f"{name}-run"
        executable.write_bytes(b"fixture")
        executable.chmod(0o755)
        archive.write_bytes(b"!<arch>\n")
        build_manifest.write_bytes(runner.canonical_pretty_bytes(
            runner.synthetic_build_manifest(inputs)))
        reservation.write_bytes(runner.canonical_bytes({
            "benchmark_cpu": 0, "nonce": f"{name}-fixture", "owner": "unit-test",
            "reserved_sibling": 1, "schema": runner.RESERVATION_SCHEMA,
            "status": "held",
        }))
        reservation.chmod(0o600)
        options = argparse.Namespace(
            output=output, executable=executable, archive=archive,
            build_manifest=build_manifest, source_root=self.root, cpu=0, sibling=1,
            timeout=1.0, expected_tooling_commit="b" * 40,
            expected_core_commit="e" * 40, reservation_file=reservation)
        pair_guard = mock.MagicMock()
        pair_identity = runner.synthetic_pair_lease(0, 1)
        pair_guard.__enter__.return_value = pair_identity
        pair_guard.identity = copy.deepcopy(pair_identity)
        pair_guard.__exit__.return_value = None
        if exit_failure:
            pair_guard.__exit__.side_effect = RuntimeError("synthetic exit failure")
        validation_calls = 0

        def validate_pair(*unused: object) -> None:
            nonlocal validation_calls
            validation_calls += 1
            if post_child_guard_failure and validation_calls == 1:
                raise runner.EvidenceError("synthetic post-child guard failure")
            if mutate_on_final_lease and validation_calls == 2:
                (output / "snapshot/child/result.json").write_bytes(b"{}\n")

        pair_guard.validate_current.side_effect = validate_pair
        snapshots = [
            runner.synthetic_cpu_stat(0, user=0, idle=0),
            runner.synthetic_cpu_stat(1, user=0, idle=0),
            runner.synthetic_cpu_stat(0, user=1, idle=0),
            runner.synthetic_cpu_stat(1, user=0, idle=1),
        ]

        def child_effect(command: object, stdout_path: Path, stderr_path: Path,
                         *unused: object, **unused_kw: object) -> object:
            result_path = output / "snapshot/child/result.json"
            if malformed_result:
                runner.write_exclusive(result_path, b"{malformed\n")
            else:
                runner.write_json_exclusive(result_path, runner.synthetic_result(0, inputs))
            runner.write_exclusive(stdout_path, b"")
            runner.write_exclusive(stderr_path, runner.expected_stderr())
            return subprocess.CompletedProcess(["fixture"], 0, b"", b""), False, 1, 2

        set_calls = 0

        def set_affinity(unused_pid: int, unused_mask: object) -> None:
            nonlocal set_calls
            set_calls += 1
            if restore_failure and set_calls == 2:
                raise OSError("synthetic restore failure")

        with mock.patch.object(runner, "validate_topology", return_value=({0, 1, 2}, {2})), \
             mock.patch.object(runner, "host_identity",
                               return_value=runner.synthetic_host(0, 1)), \
             mock.patch.object(runner, "PairLease", return_value=pair_guard), \
             mock.patch.object(runner, "input_snapshot", return_value=inputs), \
             mock.patch.object(runner, "cpu_stat_snapshot", side_effect=snapshots), \
             mock.patch.object(runner.os, "sched_getaffinity",
                               return_value={0, 1, 2}), \
             mock.patch.object(runner.os, "sched_setaffinity",
                               side_effect=set_affinity), \
             mock.patch.object(runner, "run_child_bounded", side_effect=child_effect):
            if (mutate_on_final_lease or exit_failure or restore_failure or
                    post_child_guard_failure):
                with self.assertRaises(runner.EvidenceError):
                    runner.run_campaign(options)
            else:
                self.assertEqual(runner.run_campaign(options), 0)
        failure = None
        if (output / "failure.json").exists():
            failure = runner.strict_json(
                (output / "failure.json").read_bytes(), "failure")
            runner.validate_failure(failure, output, check_files=True)
        return output, failure

    def test_terminal_publication_waits_for_final_snapshot_guard_exit_and_affinity(self) -> None:
        output, failure = self._run_success_variant(
            "final-replace", mutate_on_final_lease=True)
        self.assertIsNotNone(failure)
        self.assertEqual(failure["failure_code"], "final-snapshot-invalid")
        self.assertFalse((output / "manifest.json").exists())

        output, failure = self._run_success_variant(
            "exit-failure", exit_failure=True)
        self.assertIsNotNone(failure)
        self.assertEqual(failure["failure_code"], "guard-exit-failed")
        self.assertFalse((output / "manifest.json").exists())

        output, failure = self._run_success_variant(
            "restore-failure", restore_failure=True)
        self.assertIsNotNone(failure)
        self.assertEqual(failure["failure_code"], "affinity-restore-failed")
        self.assertFalse((output / "manifest.json").exists())

        output, failure = self._run_success_variant(
            "post-child-guard", malformed_result=True,
            post_child_guard_failure=True)
        self.assertIsNotNone(failure)
        self.assertEqual(failure["failure_code"], "lease-post-child-invalid")
        self.assertFalse((output / "manifest.json").exists())

        output, failure = self._run_success_variant("clean-success")
        self.assertIsNone(failure)
        self.assertTrue((output / "manifest.json").is_file())
        self.assertFalse((output / "failure.json").exists())

    def test_terminal_publication_is_root_bound_and_mutually_exclusive(self) -> None:
        root = self.root / "terminal-root"
        root.mkdir(mode=0o700)
        root.chmod(0o700)
        descriptor, _absolute = runner._open_directory_nofollow(
            root, "terminal test root")
        try:
            identity = runner._root_identity(descriptor, root)
        finally:
            os.close(descriptor)
        original_write = runner._write_json_exclusive_at
        displaced = self.root / "terminal-root.original"

        def replace_root(directory: int, name: str, value: dict) -> tuple[int, int]:
            root.rename(displaced)
            root.mkdir(mode=0o700)
            root.chmod(0o700)
            return original_write(directory, name, value)

        with mock.patch.object(
                runner, "_write_json_exclusive_at", side_effect=replace_root), \
             self.assertRaises(runner.EvidenceError):
            runner._publish_terminal(
                root, "manifest.json", {"probe": "root replacement"}, identity)
        self.assertFalse((root / "manifest.json").exists())
        self.assertFalse((displaced / "manifest.json").exists())

        exclusive = self.root / "terminal-exclusive"
        exclusive.mkdir(mode=0o700)
        exclusive.chmod(0o700)

        def inject_other(directory: int, name: str,
                         value: dict) -> tuple[int, int]:
            original_write(directory, "failure.json", {"probe": "other"})
            return original_write(directory, name, value)

        with mock.patch.object(
                runner, "_write_json_exclusive_at", side_effect=inject_other), \
             self.assertRaises(runner.EvidenceError):
            runner._publish_terminal(
                exclusive, "manifest.json", {"probe": "requested"})
        self.assertFalse((exclusive / "manifest.json").exists())
        self.assertTrue((exclusive / "failure.json").exists())

        manifest_path = self.root / "manifest.json"
        failure_path = self.root / "failure.json"
        runner.write_json_exclusive(manifest_path, self.manifest)
        runner.write_json_exclusive(failure_path, {"probe": "opposite terminal"})
        with self.assertRaises(runner.EvidenceError):
            runner.validate_manifest(self.manifest, self.root)

        failure, output = self._run_child_failure(timeout=True)
        runner.write_json_exclusive(output / "manifest.json", {"probe": "opposite"})
        with self.assertRaises(runner.EvidenceError):
            runner.validate_failure(failure, output, check_files=True)

    def test_nonzero_child_retains_failure_evidence(self) -> None:
        failure, _output = self._run_child_failure(timeout=False)
        self.assertEqual(failure["child"]["returncode"], 7)
        self.assertFalse(failure["child"]["timed_out"])
        self.assertEqual(failure["failure_code"], "child-exit")
        self.assertEqual(failure["error"], "C7 child exited 7")

    def test_timeout_child_retains_failure_evidence(self) -> None:
        failure, output = self._run_child_failure(timeout=True)
        self.assertEqual(failure["child"]["returncode"], 124)
        self.assertTrue(failure["child"]["timed_out"])
        self.assertEqual(failure["failure_code"], "child-timeout")
        self.assertEqual(failure["error"], "C7 child timed out")

        completed = subprocess.run(
            [sys.executable, str(MODULE_PATH), "verify-failure", "--failure",
             str(output / "failure.json")], check=False, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1",
                 "PYTHONHASHSEED": "0"})
        self.assertEqual(completed.returncode, 2, completed.stderr)
        report = json.loads(completed.stdout)
        self.assertEqual(report["status"], "VERIFIED_FAILURE")
        self.assertEqual(report["failure_code"], "child-timeout")

        hardlink = self.root / "failure-hardlink.json"
        os.link(output / "failure.json", hardlink)
        completed = subprocess.run(
            [sys.executable, str(MODULE_PATH), "verify-failure", "--failure",
             str(output / "failure.json")], check=False, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(completed.returncode, 1)
        self.assertIn("single-link regular file", completed.stderr)
        hardlink.unlink()

    def test_failure_coordinated_resigned_mutations_are_rejected(self) -> None:
        failure, output = self._run_child_failure(timeout=True)

        def drop_captured_context(value: dict) -> None:
            for key in runner.FAILURE_CONTEXT_FIELDS:
                value[key] = None
            value["completed_stage"] = runner.completed_failure_stage(-1)
            value["failure_code"] = "arguments-invalid"
            value["error"] = "synthetic preflight failure"

        def timeout_to_missing_result(value: dict) -> None:
            value["child"].update(timed_out=False, returncode=0)
            value["failure_code"] = "child-result-missing"
            value["error"] = "C7 child did not write result JSON"

        def rewrite_duration(value: dict) -> None:
            value["isolation"]["after"]["monotonic_ns"] += 1_000_000
            value["isolation"]["duration_ns"] += 1_000_000
            value["child"]["duration_ns"] += 1_000_000

        def fabricate_post_timeout(value: dict) -> None:
            value["inputs_after"] = copy.deepcopy(value["inputs_before"])

        def mutate_closure(value: dict) -> None:
            value["inputs_before"]["core_source_closure"][0]["sha256"] = "0" * 64
            payload = {key: item for key, item in value["inputs_before"].items()
                       if key != "binding_sha256"}
            value["inputs_before"]["binding_sha256"] = runner.sha256_bytes(
                runner.canonical_bytes(payload))

        def add_removed_tree_claim(value: dict) -> None:
            value["inputs_before"]["git"]["tooling_tree"] = "0" * 40
            payload = {key: item for key, item in value["inputs_before"].items()
                       if key != "binding_sha256"}
            value["inputs_before"]["binding_sha256"] = runner.sha256_bytes(
                runner.canonical_bytes(payload))

        mutations = (
            drop_captured_context,
            timeout_to_missing_result,
            rewrite_duration,
            fabricate_post_timeout,
            mutate_closure,
            add_removed_tree_claim,
            lambda value: value["host_before"]["timing_cpu"].update(cpu=False),
            lambda value: value.update(status="pass"),
            lambda value: value.update(failure_code="child-exit"),
            lambda value: value.update(error="C7 child exited 124"),
            lambda value: value["child"].update(returncode=0),
            lambda value: value["child"].update(timed_out=False),
            lambda value: (
                value.update(failure_code="child-result-missing"),
                value["child"].update(returncode=0, timed_out=False)),
            lambda value: value["child"]["stdout"].update(sha256="f" * 64),
            lambda value: value["child"]["stderr"].update(path="../outside"),
            lambda value: value["build_provenance"]["attestation"]["avx2"].update(
                optimization="-O0"),
            lambda value: value["reservation"]["payload"].update(status="released"),
            lambda value: value["pair_lease"].update(inode=99),
        )
        for index, mutation in enumerate(mutations):
            changed = copy.deepcopy(failure)
            mutation(changed)
            with self.subTest(index=index), self.assertRaises(runner.EvidenceError):
                runner.validate_failure(
                    resign(changed), output, check_files=True)

    def test_rebound_failure_rejects_mask_type_and_provenance_forgeries(self) -> None:
        failure, output = self._run_child_failure(timeout=True)

        def rebind(changed: dict, root: Path) -> dict:
            state = runner.signed(runner.failure_state_payload(
                changed["completed_stage"], changed["failure_code"],
                {key: changed[key] for key in runner.FAILURE_CONTEXT_FIELDS}))
            state_path = root / "failure/state-v2.json"
            state_path.write_bytes(runner.canonical_bytes(state) + b"\n")
            changed["state"] = runner.artifact_record(root, state_path)
            changed["artifact_inventory"] = runner.artifact_inventory(root)
            return resign(changed)

        mutations = []
        changed = copy.deepcopy(failure)
        changed["inputs_after"] = copy.deepcopy(changed["inputs_before"])
        mutations.append(changed)
        changed = copy.deepcopy(failure)
        changed["host_before"]["timing_cpu"]["cpu"] = False
        mutations.append(changed)
        changed = copy.deepcopy(failure)
        inputs = changed["inputs_before"]
        inputs["core_source_closure"][0]["sha256"] = "0" * 64
        inputs["build_attestation"]["core_source_closure_sha256"] = \
            runner.sha256_bytes(runner.canonical_bytes(
                inputs["core_source_closure"]))
        inputs["binding_sha256"] = runner.sha256_bytes(runner.canonical_bytes(
            {key: item for key, item in inputs.items()
             if key != "binding_sha256"}))
        changed["build_provenance"]["attestation"] = copy.deepcopy(
            inputs["build_attestation"])
        mutations.append(changed)
        changed = copy.deepcopy(failure)
        changed["inputs_before"]["git"]["core_tree"] = "1" * 40
        changed["inputs_before"]["binding_sha256"] = runner.sha256_bytes(
            runner.canonical_bytes({key: item for key, item in
                                    changed["inputs_before"].items()
                                    if key != "binding_sha256"}))
        mutations.append(changed)

        for index, changed in enumerate(mutations):
            root = self.root / f"rebound-{index}"
            shutil.copytree(output, root)
            rebound = rebind(changed, root)
            with self.subTest(index=index), self.assertRaises(runner.EvidenceError):
                runner.validate_failure(rebound, root, check_files=True)

    def test_failure_artifact_inventory_is_exact_and_state_path_is_fixed(self) -> None:
        failure, output = self._run_child_failure(timeout=True)
        changed = copy.deepcopy(failure)
        changed["artifact_inventory"].pop()
        with self.assertRaises(runner.EvidenceError):
            runner.validate_failure(resign(changed), output, check_files=True)

        changed = copy.deepcopy(failure)
        changed["state"]["path"] = "failure/renamed-state.json"
        with self.assertRaises(runner.EvidenceError):
            runner.validate_failure(resign(changed), output, check_files=True)

        runner.write_exclusive(output / "unexpected.bin", b"unexpected")
        with self.assertRaises(runner.EvidenceError):
            runner.validate_failure(failure, output, check_files=True)

    def test_failure_cli_bounds_numeric_conversion_and_top_level_reads(self) -> None:
        failure, output = self._run_child_failure(timeout=True)
        numeric_root = self.root / "overflow-timeout"
        shutil.copytree(output, numeric_root)
        changed = copy.deepcopy(failure)
        huge = 10 ** 400
        changed["arguments"]["timeout_seconds"] = huge
        changed["request"]["timeout_seconds"] = huge
        state = runner.signed(runner.failure_state_payload(
            changed["completed_stage"], changed["failure_code"],
            {key: changed[key] for key in runner.FAILURE_CONTEXT_FIELDS}))
        state_path = numeric_root / "failure/state-v2.json"
        state_path.write_bytes(runner.canonical_bytes(state) + b"\n")
        changed["state"] = runner.artifact_record(numeric_root, state_path)
        changed["artifact_inventory"] = runner.artifact_inventory(numeric_root)
        (numeric_root / "failure.json").write_bytes(
            runner.canonical_bytes(resign(changed)) + b"\n")

        completed = subprocess.run(
            [sys.executable, str(MODULE_PATH), "verify-failure", "--failure",
             str(numeric_root / "failure.json")], check=False, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1",
                 "PYTHONHASHSEED": "0"})
        self.assertEqual(completed.returncode, 1)
        self.assertIn("C7 authoritative evidence error:", completed.stderr)
        self.assertNotIn("Traceback (most recent call last)", completed.stderr)

        oversized_root = self.root / "oversized"
        oversized_root.mkdir()
        oversized = oversized_root / "failure.json"
        oversized.write_bytes(b" " * (runner.MAX_TOP_LEVEL_JSON_BYTES + 1))
        completed = subprocess.run(
            [sys.executable, str(MODULE_PATH), "verify-failure", "--failure",
             str(oversized)], check=False, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1",
                 "PYTHONHASHSEED": "0"})
        self.assertEqual(completed.returncode, 1)
        self.assertIn("exceeds the evidence bound", completed.stderr)
        self.assertNotIn("Traceback (most recent call last)", completed.stderr)

    def test_malformed_nested_failure_is_caught_by_cli(self) -> None:
        failure, output = self._run_child_failure(timeout=True)
        malformed_values = []
        changed = copy.deepcopy(failure)
        changed["isolation"] = {"timing_cpu": 0, "sibling_cpu": 1}
        malformed_values.append(changed)
        changed = copy.deepcopy(failure)
        changed["isolation"]["before"]["timing_cpu"] = None
        malformed_values.append(changed)
        changed = copy.deepcopy(failure)
        changed["child"] = {"timed_out": True}
        malformed_values.append(changed)
        changed = copy.deepcopy(failure)
        changed["build_provenance"] = {"manifest": None}
        malformed_values.append(changed)
        changed = copy.deepcopy(failure)
        changed["pair_lease"] = {"payload": None}
        malformed_values.append(changed)

        for index, changed in enumerate(malformed_values):
            malformed_root = self.root / f"malformed-{index}"
            shutil.copytree(output, malformed_root)
            path = malformed_root / "failure.json"
            path.write_bytes(runner.canonical_bytes(resign(changed)) + b"\n")
            completed = subprocess.run(
                [sys.executable, str(MODULE_PATH), "verify-failure", "--failure",
                 str(path)], check=False, text=True, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1",
                     "PYTHONHASHSEED": "0"})
            with self.subTest(index=index):
                self.assertEqual(completed.returncode, 1)
                self.assertIn("C7 authoritative evidence error:", completed.stderr)
                self.assertNotIn("Traceback (most recent call last)", completed.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
