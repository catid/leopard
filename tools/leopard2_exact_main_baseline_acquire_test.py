#!/usr/bin/env python3
"""Deterministic tests for exact-main acquisition and lane sealing."""

from __future__ import annotations

import ast
import contextlib
import copy
import errno
import hashlib
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
from types import SimpleNamespace
import unittest
from unittest import mock


TOOLS = Path(__file__).resolve().parent
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))
import leopard2_exact_main_baseline as identity_contract  # noqa: E402
import leopard2_exact_main_baseline_acquire as acquire  # noqa: E402
import leopard2_exact_main_baseline_record as record_contract  # noqa: E402
import leopard2_exact_main_baseline_verifier as verifier  # noqa: E402


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load test dependency {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


verifier_fixtures = load_module(
    "leopard2_exact_main_baseline_acquire_verifier_fixtures",
    TOOLS / "leopard2_exact_main_baseline_verifier_test.py",
)


def sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def resign_for_root(record: dict, root: Path) -> dict:
    value = copy.deepcopy(record)
    value["lane"]["root"] = str(root)
    value["record_sha256"] = identity_contract.canonical_sha256({
        key: copy.deepcopy(field)
        for key, field in value.items()
        if key != "record_sha256"
    })
    if value["schema"] == record_contract.AUTHORITY_SCHEMA:
        return record_contract.validate_baseline_authority_record(value)
    return record_contract.validate_baseline_failure_record(value)


@contextlib.contextmanager
def fixture_anchors():
    with contextlib.ExitStack() as stack:
        for module in (record_contract, verifier):
            stack.enter_context(mock.patch.object(
                module, "BASELINE_COMMIT",
                verifier_fixtures._TEST_BASELINE_COMMIT))
            stack.enter_context(mock.patch.object(
                module, "BASELINE_TREE",
                verifier_fixtures._TEST_BASELINE_TREE))
        yield


def fixture_for_root(root: Path, kind: str) -> tuple[dict[str, bytes], dict]:
    with fixture_anchors():
        if kind == "authority":
            files, record = verifier_fixtures.authority_files_and_record()
            terminal = "baseline-authority.json"
        else:
            files, record = verifier_fixtures.failure_files_and_record(
                verification_failure=kind == "verification_failure")
            terminal = "FAILED.json"
        relocated = resign_for_root(record, root)
    retained = copy.deepcopy(files)
    retained.pop(terminal)
    return retained, relocated


def lane_plan(root: str = "/tmp/leopard-exact-main-lane-a1") -> acquire.LanePlan:
    return acquire.LanePlan(
        lane_root=root,
        attempt=1,
        repository="/home/catid/leopard",
        verifier=(
            "/home/catid/leopard/tools/"
            "leopard2_exact_main_baseline_verifier.py"),
        canonical_adapter_root="/tmp/leopard-acquire-canonical-adapter",
        canonical_baseline_root="/tmp/leopard-acquire-canonical-baseline",
        canonical_build_root="/tmp/leopard-acquire-canonical-build",
        variant_adapter_root="/tmp/leopard-acquire-variant-adapter",
        variant_baseline_root="/tmp/leopard-acquire-variant-baseline",
        variant_build_root="/tmp/leopard-acquire-variant-build",
    )


class ExactMainBaselineAcquireTest(unittest.TestCase):
    def assertAcquireError(self, function, *arguments, **keywords) -> None:  # noqa: N802
        with self.assertRaises(acquire.AcquisitionError):
            function(*arguments, **keywords)

    def test_lane_plan_is_detached_and_roots_are_disjoint(self) -> None:
        plan = lane_plan()
        self.assertEqual(acquire.validate_lane_plan(plan), plan)
        for mutation in (
                {"attempt": 0}, {"attempt": 4},
                {"lane_root": "relative/lane"},
                {"variant_build_root": plan.canonical_build_root + "/child"},
                {"variant_build_root": "/var" + plan.canonical_build_root +
                 "-suffix"},
                {"canonical_adapter_root": plan.canonical_baseline_root},
                {"lane_root": plan.repository + "/lane"},
                {"canonical_build_root": plan.repository + "/build"},
                {"verifier": plan.canonical_build_root + "/verifier.py"},
                {"verifier": plan.repository}):
            with self.subTest(mutation=mutation):
                value = acquire.LanePlan(**{
                    **plan.__dict__, **mutation,
                })
                self.assertAcquireError(acquire.validate_lane_plan, value)
        self.assertAcquireError(acquire.validate_lane_plan, plan.__dict__)

    def test_canonical_ldd_and_build_closure_helpers(self) -> None:
        rows = [
            {"soname": "libc.so.6", "kind": "file",
             "path": "/usr/lib/libc.so.6"},
            {"soname": "linux-vdso.so.1", "kind": "virtual",
             "path": None},
        ]
        content = acquire.canonical_ldd_text(rows)
        self.assertEqual(
            list(record_contract.parse_canonical_ldd_output(content)), rows)
        self.assertAcquireError(acquire.canonical_ldd_text,
                                list(reversed(rows)))
        self.assertAcquireError(acquire.canonical_ldd_text, [rows[0], rows[0]])
        self.assertAcquireError(
            acquire.canonical_ldd_text,
            [rows[0]] * (record_contract.MAX_DEPENDENCIES + 1))
        oversized_soname = copy.deepcopy(rows[0])
        oversized_soname["soname"] = \
            "a" * (record_contract.MAX_TEXT_LENGTH + 1)
        self.assertAcquireError(
            acquire.canonical_ldd_text, [oversized_soname])
        files = [
            {"relative_path": "leopard_main_benchmark", "size": 5,
             "sha256": sha256(b"main\n")},
            {"relative_path": "libleopard_main_exact.a", "size": 8,
             "sha256": sha256(b"!<arch>\n")},
        ]
        files.sort(key=lambda item: item["relative_path"])
        closure = acquire.build_closure_document(
            "canonical_first", "/tmp/exact-main-build", files)
        self.assertEqual(closure, {
            "schema": acquire.BUILD_CLOSURE_SCHEMA,
            "role": "canonical_first",
            "build_root": "/tmp/exact-main-build",
            "files": files,
            "file_count": 2,
        })
        self.assertAcquireError(
            acquire.build_closure_document, "canonical_first",
            "/tmp/exact-main-build", list(reversed(files)))

    def test_raw_ldd_normalization_is_exact_and_fail_closed(self) -> None:
        raw = (
            b"\tlinux-vdso.so.1 (0x00007ffe01234000)\n"
            b"\tlibc.so.6 => /lib/x86_64-linux-gnu/libc.so.6 "
            b"(0x0000700000000000)\n"
            b"\t/lib64/ld-linux-x86-64.so.2 (0x0000700000100000)\n"
        )
        rows = acquire.normalize_ldd_output(raw)
        self.assertEqual(rows, (
            {
                "soname": "ld-linux-x86-64.so.2", "kind": "file",
                "path": "/lib64/ld-linux-x86-64.so.2",
            },
            {
                "soname": "libc.so.6", "kind": "file",
                "path": "/lib/x86_64-linux-gnu/libc.so.6",
            },
            {
                "soname": "linux-vdso.so.1", "kind": "virtual",
                "path": None,
            },
        ))
        self.assertEqual(
            record_contract.parse_canonical_ldd_output(
                acquire.canonical_ldd_text(rows)), rows)
        invalid = (
            b"libmissing.so => not found\n",
            b"statically linked\n",
            raw.replace(b"\n", b"\r\n"),
            raw.replace(b"linux-vdso.so.1", b"libc.so.6"),
            raw + raw.splitlines(keepends=True)[1],
            (
                b"libalias.so.1 => /lib/x86_64-linux-gnu/libc.so.6 "
                b"(0x0000700000200000)\n" + raw
            ),
            b"\xff\n",
            b"\n",
        )
        for changed in invalid:
            with self.subTest(changed=changed[:80]):
                self.assertAcquireError(
                    acquire.normalize_ldd_output, changed)

    def test_host_command_capture_is_bounded_and_exact(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-exact-main-command.") as temporary:
            environment = acquire.HostEnvironment(
                anchor_path=temporary,
                canonical_lock_path=str(Path(temporary) / "canonical.lock"))
            subreaper_before = acquire._linux_subreaper_state()
            command_environment = {
                "LANG": "C", "LC_ALL": "C", "PATH": "/usr/bin:/bin",
                "TOKEN": "frozen",
            }
            result = environment.run(
                [sys.executable, "-c", (
                    "import os,sys;"
                    "sys.stdout.buffer.write(os.environ['TOKEN'].encode());"
                    "sys.stderr.buffer.write(b'err');"
                    "raise SystemExit(7)")],
                cwd=temporary,
                env=command_environment,
                timeout=5,
                maximum_bytes=1024,
            )
            self.assertEqual(result, acquire.CommandResult(
                exit_status=7, stdout=b"frozen", stderr=b"err"))
            with self.assertRaisesRegex(
                    acquire.AcquisitionError, "combined output bound"):
                environment.run(
                    [sys.executable, "-c",
                     "import sys;sys.stdout.buffer.write(b'x'*65536)"],
                    cwd=temporary, env=command_environment,
                    timeout=5, maximum_bytes=31)
            with self.assertRaisesRegex(
                    acquire.AcquisitionError, "timeout"):
                environment.run(
                    [sys.executable, "-c", "import time;time.sleep(30)"],
                    cwd=temporary, env=command_environment,
                    timeout=0.1, maximum_bytes=1024)
            started = time.monotonic()
            with self.assertRaisesRegex(
                    acquire.AcquisitionError, "timeout"):
                environment.run(
                    [sys.executable, "-c",
                     "import os,time;os.close(1);os.close(2);time.sleep(30)"],
                    cwd=temporary, env=command_environment,
                    timeout=0.2, maximum_bytes=1024)
            self.assertLess(time.monotonic() - started, 2.5)
            escaped_pid = Path(temporary) / "escaped.pid"
            child = (
                "import os,time;"
                "os.setsid();"
                f"p={str(escaped_pid)!r};"
                "fd=os.open(p,os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o600);"
                "os.write(fd,str(os.getpid()).encode());os.fsync(fd);os.close(fd);"
                "time.sleep(30)"
            )
            parent = (
                "import os,subprocess,sys,time;"
                f"subprocess.Popen([sys.executable,'-c',{child!r}]);"
                f"p={str(escaped_pid)!r};"
                "deadline=time.monotonic()+10;"
                "\nwhile not os.path.isfile(p):\n"
                "  assert time.monotonic()<deadline\n"
                "  time.sleep(.01)\n"
                "time.sleep(30)"
            )
            with self.assertRaisesRegex(
                    acquire.AcquisitionError, "timeout"):
                environment.run(
                    [sys.executable, "-c", parent],
                    cwd=temporary, env=command_environment,
                    timeout=2, maximum_bytes=1024)
            self.assertTrue(escaped_pid.is_file())
            escaped = int(escaped_pid.read_text(encoding="ascii"))
            self.assertFalse(Path(f"/proc/{escaped}").exists())
            self.assertEqual(
                acquire._linux_subreaper_state(), subreaper_before)
            for invalid in (
                    {"argv": []},
                    {"timeout": True},
                    {"timeout": 0},
                    {"maximum_bytes": -1},
                    {"env": {"BAD=NAME": "value"}},
                    {"env": {1: "value"}},
                    {"cwd": "relative"}):
                arguments = {
                    "argv": [sys.executable, "-c", "pass"],
                    "cwd": temporary,
                    "env": command_environment,
                    "timeout": 5,
                    "maximum_bytes": 1024,
                }
                arguments.update(invalid)
                with self.subTest(invalid=invalid):
                    self.assertAcquireError(environment.run, **arguments)

    def test_host_facts_and_timestamp_fit_the_record_contract(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-exact-main-host.") as temporary:
            environment = acquire.HostEnvironment(
                anchor_path=temporary,
                canonical_lock_path=str(Path(temporary) / "canonical.lock"))
            facts = environment.host_facts()
            self.assertEqual(set(facts), record_contract.HOST_KEYS)
            self.assertEqual(facts["architecture"], "x86_64")
            self.assertEqual(facts["online_cpus"],
                             sorted(set(facts["online_cpus"])))
            expected_online = acquire._parse_linux_cpu_list(
                acquire._read_bounded_system_file(
                    "/sys/devices/system/cpu/online", 65536,
                    "test online CPU list"))
            self.assertEqual(facts["online_cpus"], expected_online)
            self.assertIn(
                b"model name",
                acquire._read_bounded_system_file(
                    "/proc/cpuinfo", acquire.MAX_CPUINFO_BYTES,
                    "test CPU identity"))
            self.assertEqual(record_contract._host(facts), facts)
            self.assertRegex(
                environment.now_utc(),
                r"^[0-9]{4}-[0-9]{2}-[0-9]{2}T"
                r"[0-9]{2}:[0-9]{2}:[0-9]{2}Z$")

            cpuinfo = (
                b"processor : 0\nmodel name : uniform model\n\n"
                b"processor : 1\nmodel name : uniform model\n")
            online = b"0-1\n"
            with mock.patch.object(
                    acquire, "_read_bounded_system_file",
                    side_effect=(cpuinfo, online)), \
                    mock.patch.object(os, "uname", return_value=type(
                        "Uname", (), {
                            "nodename": "fixture", "release": "fixture-kernel",
                            "machine": "x86_64",
                        })()), \
                    mock.patch.object(os, "sysconf", return_value=100):
                synthetic = environment.host_facts()
            self.assertEqual(synthetic["cpu_model"], "uniform model")
            self.assertEqual(synthetic["online_cpus"], [0, 1])
            mixed = cpuinfo.replace(
                b"processor : 1\nmodel name : uniform model",
                b"processor : 1\nmodel name : other model")
            with mock.patch.object(
                    acquire, "_read_bounded_system_file",
                    side_effect=(mixed, online)), \
                    mock.patch.object(os, "uname", return_value=type(
                        "Uname", (), {
                            "nodename": "fixture", "release": "fixture-kernel",
                            "machine": "x86_64",
                        })()), \
                    mock.patch.object(os, "sysconf", return_value=100):
                self.assertAcquireError(environment.host_facts)

    def test_proc_permission_and_pidfd_faults_are_fail_closed(self) -> None:
        denied = PermissionError(errno.EACCES, "denied")
        with mock.patch.object(Path, "read_bytes", side_effect=denied), \
                mock.patch.object(
                    os, "stat",
                    return_value=SimpleNamespace(st_uid=os.getuid() + 1)):
            self.assertIsNone(acquire._linux_process_identity(123456))
        with mock.patch.object(Path, "read_bytes", side_effect=denied), \
                mock.patch.object(
                    os, "stat",
                    return_value=SimpleNamespace(st_uid=os.getuid())):
            self.assertAcquireError(
                acquire._linux_process_identity, 123456)
        descriptor = os.open("/dev/null", os.O_RDONLY | os.O_CLOEXEC)
        containment = acquire._LinuxChildContainment()
        with mock.patch.object(
                acquire, "_linux_pidfd_open", return_value=descriptor), \
                mock.patch.object(
                    acquire, "_linux_process_identity",
                    side_effect=acquire.AcquisitionError("injected stat fault")):
            self.assertAcquireError(containment._retain, 123456, 1)
        with self.assertRaises(OSError):
            os.fstat(descriptor)

    def test_process_snapshot_tolerates_owner_reparenting(self) -> None:
        owner = os.getpid()
        with mock.patch.object(os, "listdir", return_value=[str(owner)]), \
                mock.patch.object(
                    acquire, "_linux_process_identity",
                    side_effect=((100, 123456), (200, 123456))):
            self.assertEqual(
                acquire._linux_process_snapshot(),
                {owner: (100, 123456)})
        with mock.patch.object(os, "listdir", return_value=[str(owner)]), \
                mock.patch.object(
                    acquire, "_linux_process_identity",
                    side_effect=((100, 123456), (200, 123457))):
            self.assertAcquireError(acquire._linux_process_snapshot)

    def test_containment_capacity_is_dynamic_and_resource_faults_stick(self) \
            -> None:
        with mock.patch.object(
                acquire.resource, "getrlimit", return_value=(1024, 1024)):
            with acquire._LinuxChildContainment() as containment:
                self.assertGreaterEqual(
                    containment.process_bound,
                    acquire.MIN_CONTAINED_PROCESSES)
                self.assertLess(containment.process_bound, 1024)
        with mock.patch.object(
                acquire.resource, "getrlimit", return_value=(32, 32)):
            self.assertAcquireError(
                acquire._LinuxChildContainment().__enter__)

        containment = acquire._LinuxChildContainment()
        containment.process_bound = 0
        containment._retain(os.getpid(), 1)
        self.assertTrue(containment.containment_bound_reached)

        containment = acquire._LinuxChildContainment()
        containment.reserve_descriptor = os.open(
            "/dev/null", os.O_RDONLY | os.O_CLOEXEC)
        retained_descriptor = os.open(
            "/dev/null", os.O_RDONLY | os.O_CLOEXEC)
        with mock.patch.object(
                acquire, "_linux_pidfd_open",
                side_effect=(OSError(errno.EMFILE, "full"),
                             retained_descriptor)), \
                mock.patch.object(
                    acquire, "_linux_process_identity",
                    return_value=(os.getpid(), 1)):
            containment._retain(os.getpid(), 1)
        self.assertTrue(containment.pidfd_resource_exhausted)
        self.assertEqual(containment.reserve_descriptor, -1)
        self.assertIn((os.getpid(), 1), containment.handles)
        for descriptor in containment.handles.values():
            os.close(descriptor)
        containment.handles.clear()

        self.assertAcquireError(
            acquire._parse_linux_cpu_list, b"0-2147483647\n")

    def test_acquisition_locks_reject_contention_and_replacement(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-exact-main-lock-parent.") as parent:
            anchor = Path(parent) / "anchor"
            anchor.mkdir(mode=0o700)
            canonical = Path(parent) / "canonical.lock"
            environment = acquire.HostEnvironment(
                anchor_path=str(anchor),
                canonical_lock_path=str(canonical))
            with acquire.AcquisitionLocks(environment, blocking=False) as locks:
                locks.validate_current()
                self.assertAcquireError(
                    acquire.StableLeaseAnchor(str(anchor)).__enter__)
                self.assertAcquireError(
                    acquire.CanonicalFileLock(
                        str(canonical), blocking=False).__enter__)
                moved = Path(parent) / "old-anchor"
                anchor.rename(moved)
                anchor.mkdir(mode=0o700)
                self.assertAcquireError(locks.validate_current)
            canonical.unlink()
            with acquire.CanonicalFileLock(
                    str(canonical), blocking=False) as lock:
                replacement = Path(parent) / "replacement.lock"
                replacement.write_bytes(b"")
                replacement.chmod(0o600)
                os.replace(replacement, canonical)
                self.assertAcquireError(lock.validate_current)

    def test_lock_metadata_and_path_errors_are_normalized(self) -> None:
        for kind in ("anchor_mode", "missing_anchor", "lock_mode",
                     "lock_hardlink", "lock_directory", "lock_symlink"):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory(
                    prefix="leopard-exact-main-lock-metadata.") as parent:
                anchor = Path(parent) / "anchor"
                if kind != "missing_anchor":
                    anchor.mkdir(mode=0o700)
                if kind == "anchor_mode":
                    anchor.chmod(0o777)
                lock_path = Path(parent) / "canonical.lock"
                if kind == "lock_mode":
                    lock_path.write_bytes(b"")
                    lock_path.chmod(0o644)
                elif kind == "lock_hardlink":
                    lock_path.write_bytes(b"")
                    lock_path.chmod(0o600)
                    os.link(lock_path, Path(parent) / "second-name")
                elif kind == "lock_directory":
                    lock_path.mkdir(mode=0o700)
                elif kind == "lock_symlink":
                    target = Path(parent) / "target.lock"
                    target.write_bytes(b"")
                    target.chmod(0o600)
                    lock_path.symlink_to(target)
                environment = acquire.HostEnvironment(
                    anchor_path=str(anchor),
                    canonical_lock_path=str(lock_path))
                if kind in ("anchor_mode", "missing_anchor"):
                    with self.assertRaises(acquire.AcquisitionError):
                        acquire.StableLeaseAnchor(str(anchor)).__enter__()
                else:
                    with self.assertRaises(acquire.AcquisitionError):
                        acquire.CanonicalFileLock(
                            str(lock_path), blocking=False).__enter__()

    def test_authority_and_failure_lanes_verify_offline(self) -> None:
        for kind, expected_status in (
                ("authority", 0),
                ("acquisition_failure", 3),
                ("verification_failure", 3)):
            with self.subTest(kind=kind), fixture_anchors(), \
                    tempfile.TemporaryDirectory(
                        prefix="leopard-exact-main-acquire-test.") as temporary:
                root = Path(temporary) / "lane"
                retained, record = fixture_for_root(root, kind)
                with acquire.LaneWriter(str(root)) as writer:
                    seal = writer.seal_record(record, retained)
                self.assertEqual(seal["terminal"],
                                 "baseline-authority.json" if kind ==
                                 "authority" else "FAILED.json")
                self.assertEqual(oct(os.stat(root).st_mode & 0o7777), "0o500")
                self.assertFalse(any(
                    path.name.startswith(".leopard-stage-")
                    for path in root.rglob("*")))
                for path in root.rglob("*"):
                    mode = os.lstat(path).st_mode & 0o7777
                    self.assertEqual(
                        mode, acquire.LANE_DIRECTORY_MODE if path.is_dir()
                        else acquire.LANE_FILE_MODE)
                completed = subprocess.run(
                    verifier_fixtures.verifier_cli_command(root),
                    check=False, stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE, timeout=30)
                self.assertEqual(
                    completed.returncode, expected_status,
                    completed.stderr.decode("utf-8", "replace"))

    def test_lane_root_is_never_reused_and_publish_is_no_replace(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-exact-main-acquire-collision.") as temporary:
            root = Path(temporary) / "lane"
            with acquire.LaneWriter(str(root)) as writer:
                collision = root / "collision.bin"
                collision.write_bytes(b"preexisting\n")
                self.assertAcquireError(
                    writer.publish_bytes, "collision.bin", b"replacement\n")
                self.assertEqual(collision.read_bytes(), b"preexisting\n")
                self.assertFalse(any(
                    path.name.startswith(".leopard-stage-")
                    for path in root.iterdir()))
            self.assertAcquireError(acquire.LaneWriter, str(root))

    def test_seal_rejects_inventory_and_byte_drift(self) -> None:
        for mutation in ("missing", "extra", "changed"):
            with self.subTest(mutation=mutation), fixture_anchors(), \
                    tempfile.TemporaryDirectory(
                        prefix="leopard-exact-main-acquire-inventory.") \
                    as temporary:
                root = Path(temporary) / "lane"
                retained, record = fixture_for_root(
                    root, "acquisition_failure")
                if mutation == "missing":
                    retained.pop(sorted(retained)[0])
                elif mutation == "extra":
                    retained["diagnostics/extra"] = b"extra\n"
                else:
                    retained[sorted(retained)[0]] += b"changed\n"
                with acquire.LaneWriter(str(root)) as writer:
                    self.assertAcquireError(
                        writer.seal_record, record, retained)

    def test_terminal_root_and_record_shape_are_bound(self) -> None:
        with fixture_anchors(), tempfile.TemporaryDirectory(
                prefix="leopard-exact-main-acquire-record.") as temporary:
            root = Path(temporary) / "lane"
            retained, record = fixture_for_root(
                root, "acquisition_failure")
            wrong_root = copy.deepcopy(record)
            wrong_root["lane"]["root"] = str(root) + "-other"
            wrong_root["record_sha256"] = identity_contract.canonical_sha256({
                key: copy.deepcopy(value)
                for key, value in wrong_root.items()
                if key != "record_sha256"
            })
            with acquire.LaneWriter(str(root)) as writer:
                self.assertAcquireError(
                    writer.seal_record, wrong_root, retained)

    def test_final_seal_rebinds_every_published_byte(self) -> None:
        with fixture_anchors(), tempfile.TemporaryDirectory(
                prefix="leopard-exact-main-acquire-substitution.") \
                as temporary:
            root = Path(temporary) / "lane"
            retained, record = fixture_for_root(
                root, "acquisition_failure")
            with acquire.LaneWriter(str(root)) as writer:
                apply_modes = writer._apply_final_modes

                def substitute_then_lock(files, directories):
                    terminal = root / "FAILED.json"
                    content = terminal.read_bytes()
                    terminal.write_bytes(b"X" + content[1:])
                    apply_modes(files, directories)

                with mock.patch.object(
                        writer, "_apply_final_modes",
                        side_effect=substitute_then_lock):
                    self.assertAcquireError(
                        writer.seal_record, record, retained)

    def test_production_boundary_is_independent_and_optimized_safe(self) -> None:
        source = (TOOLS / "leopard2_exact_main_baseline_acquire.py").read_text(
            encoding="utf-8")
        tree = ast.parse(source)
        imported = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.update(
                    alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                imported.add((node.module or "").split(".")[0])
        self.assertNotIn(
            "leopard2_exact_main_baseline_verifier", imported)
        contract_calls = [
            node for node in ast.walk(tree)
            if isinstance(node, ast.Call) and
            isinstance(node.func, ast.Name) and
            node.func.id == "_load_local_contract"
        ]
        self.assertEqual(len(contract_calls), 2)
        loaded_contracts = set()
        for call in contract_calls:
            self.assertEqual(call.keywords, [])
            self.assertEqual(len(call.args), 2)
            self.assertTrue(all(
                isinstance(argument, ast.Constant) and
                type(argument.value) is str
                for argument in call.args))
            loaded_contracts.add(call.args[1].value)
        self.assertEqual(loaded_contracts, {
            "leopard2_exact_main_baseline.py",
            "leopard2_exact_main_baseline_record.py",
        })
        self.assertFalse(any(isinstance(node, ast.Assert)
                             for node in ast.walk(tree)))


if __name__ == "__main__":
    unittest.main()
