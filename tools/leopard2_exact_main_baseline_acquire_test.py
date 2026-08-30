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
        scratch_root="/tmp/leopard-acquire-scratch",
    )


def temporary_lane_plan(parent: Path) -> acquire.LanePlan:
    return acquire.LanePlan(
        lane_root=str(parent / "sealed-lane"), attempt=1,
        repository=str(TOOLS.parent),
        verifier=str(TOOLS / "leopard2_exact_main_baseline_verifier.py"),
        canonical_adapter_root=str(parent / "canonical-adapter"),
        canonical_baseline_root=str(parent / "canonical-baseline"),
        canonical_build_root=str(parent / "canonical-build"),
        variant_adapter_root=str(parent / "variant-adapter"),
        variant_baseline_root=str(parent / "variant-baseline"),
        variant_build_root=str(parent / "variant-build"),
        scratch_root=str(parent / "scratch"),
    )


def git_command(root: Path, *arguments: str) -> bytes:
    completed = subprocess.run(
        ["/usr/bin/git", "-C", str(root), *arguments], check=False,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        env={
            **os.environ,
            "GIT_CONFIG_GLOBAL": "/dev/null",
            "GIT_CONFIG_NOSYSTEM": "1",
            "GIT_CONFIG_SYSTEM": "/dev/null",
            "GIT_ATTR_NOSYSTEM": "1",
            "GIT_NO_REPLACE_OBJECTS": "1",
            "GIT_OPTIONAL_LOCKS": "0",
            "GIT_AUTHOR_DATE": "2000-01-01T00:00:00Z",
            "GIT_COMMITTER_DATE": "2000-01-01T00:00:00Z",
        }, timeout=30)
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr.decode("utf-8", "replace"))
    return completed.stdout


def tiny_repository(parent: Path) -> tuple[Path, str, str, Path, str]:
    submodule = parent / "submodule"
    submodule.mkdir(mode=0o700)
    git_command(submodule, "init", "--quiet")
    (submodule / "submodule.txt").write_bytes(b"submodule fixture\n")
    git_command(submodule, "add", "submodule.txt")
    git_command(
        submodule, "-c", "user.name=Fixture", "-c",
        "user.email=fixture@example.invalid", "commit", "--quiet", "-m",
        "submodule fixture")
    submodule_commit = git_command(
        submodule, "rev-parse", "HEAD^{commit}").decode("ascii").strip()

    repository = parent / "repository"
    repository.mkdir(mode=0o700)
    git_command(repository, "init", "--quiet")
    (repository / "tracked.txt").write_bytes(b"tracked fixture\n")
    for index, relative in enumerate(record_contract.ADAPTER_PATHS):
        path = repository / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(f"adapter fixture {index}\n".encode("ascii"))
    verifier_path = repository / "tools/verifier.py"
    verifier_path.parent.mkdir(parents=True, exist_ok=True)
    verifier_path.write_bytes(b"#!/usr/bin/env python3\n")
    git_command(repository, "add", "tracked.txt", "tools/verifier.py",
                *record_contract.ADAPTER_PATHS)
    git_command(
        repository, "-c", "protocol.file.allow=always", "submodule", "add",
        "--quiet", str(submodule), "sse2neon")
    git_command(
        repository, "-c", "user.name=Fixture", "-c",
        "user.email=fixture@example.invalid", "commit", "--quiet", "-m",
        "repository fixture")
    commit = git_command(
        repository, "rev-parse", "HEAD^{commit}").decode("ascii").strip()
    tree = git_command(
        repository, "rev-parse", "HEAD^{tree}").decode("ascii").strip()
    return repository, commit, tree, submodule, submodule_commit


def git_object_inodes(root: Path) -> set[tuple[int, int]]:
    common = Path(git_command(
        root, "rev-parse", "--path-format=absolute", "--git-common-dir").
        decode("ascii").strip()).resolve()
    inodes = set()
    for path in (common / "objects").rglob("*"):
        status = path.lstat()
        if path.is_file() and not path.is_symlink():
            inodes.add((status.st_dev, status.st_ino))
    return inodes


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
                {"scratch_root": plan.variant_build_root},
                {"scratch_root": plan.canonical_build_root + "/scratch"},
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

    def test_frozen_child_environment_is_exact(self) -> None:
        expected = {
            "GIT_ATTR_NOSYSTEM": "1",
            "GIT_CONFIG_GLOBAL": "/dev/null",
            "GIT_CONFIG_NOSYSTEM": "1",
            "GIT_CONFIG_SYSTEM": "/dev/null",
            "GIT_NO_REPLACE_OBJECTS": "1",
            "GIT_OPTIONAL_LOCKS": "0",
            "LANG": "C",
            "LC_ALL": "C",
            "PATH": "/usr/bin:/bin",
            "TZ": "UTC",
        }
        self.assertEqual(acquire.frozen_child_environment(), expected)
        self.assertEqual(
            record_contract.exact_main_build_profile()["environment"],
            [{"name": name, "value": value}
             for name, value in expected.items()],
        )

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

    def test_prepared_roots_are_exclusive_canonical_and_cleaned(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-exact-main-roots.") as temporary:
            parent = Path(temporary).resolve()
            plan = temporary_lane_plan(parent)
            with acquire.prepare_acquisition_roots(plan) as prepared:
                prepared.validate_current()
                for field in acquire._PREPARED_ROOT_FIELDS:
                    path = Path(getattr(plan, field))
                    self.assertTrue(path.is_dir())
                    self.assertEqual(path.resolve(), path)
                    self.assertEqual(path.stat().st_mode & 0o7777, 0o700)
                (Path(plan.scratch_root) / "transient").write_bytes(b"x")
            for field in acquire._PREPARED_ROOT_FIELDS:
                self.assertFalse(Path(getattr(plan, field)).exists())
            self.assertFalse(Path(plan.lane_root).exists())

            stale = Path(plan.canonical_adapter_root)
            stale.mkdir(mode=0o700)
            self.assertAcquireError(
                acquire.prepare_acquisition_roots(plan).__enter__)
            self.assertTrue(stale.is_dir())
            for field in acquire._PREPARED_ROOT_FIELDS[1:]:
                self.assertFalse(Path(getattr(plan, field)).exists())

    def test_prepared_root_open_failures_close_and_normalize(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-exact-main-root-open.") as temporary:
            parent = Path(temporary).resolve()
            plan = temporary_lane_plan(parent)
            first = getattr(plan, acquire._PREPARED_ROOT_FIELDS[0])
            original_mkdir = os.mkdir
            original_open = os.open
            retained_descriptors = []

            def unsafe_mkdir(path, mode=0o777, *, dir_fd=None):
                result = original_mkdir(path, mode, dir_fd=dir_fd)
                if path == first:
                    os.chmod(path, 0o755)
                return result

            def recording_open(path, flags, mode=0o777, *, dir_fd=None):
                descriptor = original_open(
                    path, flags, mode, dir_fd=dir_fd)
                if path == first and flags & os.O_DIRECTORY:
                    retained_descriptors.append(descriptor)
                return descriptor

            with mock.patch.object(os, "mkdir", side_effect=unsafe_mkdir), \
                    mock.patch.object(os, "open", side_effect=recording_open):
                self.assertAcquireError(
                    acquire.prepare_acquisition_roots(plan).__enter__)
            self.assertTrue(retained_descriptors)
            for descriptor in retained_descriptors:
                with self.assertRaises(OSError):
                    os.fstat(descriptor)
            self.assertFalse(Path(first).exists())

        with tempfile.TemporaryDirectory(
                prefix="leopard-exact-main-root-denied.") as temporary:
            parent = Path(temporary).resolve()
            plan = temporary_lane_plan(parent)
            first = getattr(plan, acquire._PREPARED_ROOT_FIELDS[0])
            original_open = os.open
            denied = False

            def denied_open(path, flags, mode=0o777, *, dir_fd=None):
                nonlocal denied
                if (path == first and flags & os.O_DIRECTORY and
                        not denied):
                    denied = True
                    raise OSError(errno.EACCES, "injected root denial")
                return original_open(path, flags, mode, dir_fd=dir_fd)

            with mock.patch.object(os, "open", side_effect=denied_open):
                with self.assertRaisesRegex(
                        acquire.AcquisitionError,
                        "cannot retain acquisition root"):
                    acquire.prepare_acquisition_roots(plan).__enter__()
            self.assertFalse(Path(first).exists())

    def test_streamed_command_output_is_bounded_and_exact(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-exact-main-stream.") as temporary:
            root = Path(temporary).resolve()
            environment = acquire.HostEnvironment(anchor_path=str(root))
            destination = root / "output.bin"
            descriptor = os.open(
                destination,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC,
                0o600)
            try:
                result = environment.run_to_path(
                    [sys.executable, "-c",
                     "import sys;sys.stdout.buffer.write(b'x'*131071)"],
                    cwd=str(root), env={}, timeout=10,
                    destination_fd=descriptor, maximum_bytes=131071)
            finally:
                os.close(descriptor)
            self.assertEqual(result.exit_status, 0)
            self.assertEqual(result.stdout_size, 131071)
            self.assertEqual(
                result.stdout_sha256, sha256(destination.read_bytes()))
            self.assertEqual(result.stderr, b"")

            oversized = root / "oversized.bin"
            descriptor = os.open(
                oversized,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC,
                0o600)
            try:
                self.assertAcquireError(
                    environment.run_to_path,
                    [sys.executable, "-c",
                     "import sys;sys.stdout.buffer.write(b'abcdef')"],
                    cwd=str(root), env={}, timeout=10,
                    destination_fd=descriptor, maximum_bytes=5)
            finally:
                os.close(descriptor)

    def test_detached_source_and_canonical_archive_are_replayed(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-exact-main-source.") as temporary:
            parent = Path(temporary).resolve()
            repository, commit, tree, submodule, submodule_commit = \
                tiny_repository(parent)
            destination = parent / "detached"
            destination.mkdir(mode=0o700)
            scratch = parent / "scratch"
            scratch.mkdir(mode=0o700)
            environment = acquire.HostEnvironment(anchor_path=str(parent))
            child_environment = acquire.frozen_child_environment()
            command_log = []
            acquire.stage_detached_source(
                environment, source_repository=str(repository),
                submodule_repository=str(submodule),
                destination=str(destination), commit=commit, tree=tree,
                submodule_commit=submodule_commit,
                child_environment=child_environment, log=command_log)
            self.assertEqual(
                git_command(destination, "rev-parse", "HEAD^{commit}").
                decode("ascii").strip(), commit)
            self.assertTrue(git_object_inodes(repository))
            self.assertTrue(git_object_inodes(destination))
            self.assertTrue(git_object_inodes(submodule))
            self.assertTrue(git_object_inodes(destination / "sse2neon"))
            self.assertTrue(git_object_inodes(repository).isdisjoint(
                git_object_inodes(destination)))
            self.assertTrue(git_object_inodes(submodule).isdisjoint(
                git_object_inodes(destination / "sse2neon")))
            capture = environment.capture_git_identity(
                str(TOOLS.parent), str(destination), commit,
                require_detached=True)
            self.assertEqual(capture["head"], commit)
            self.assertTrue(capture["detached"])
            self.assertEqual(
                capture["submodules"][0]["object_id"], submodule_commit)
            first = acquire.canonical_git_archive(
                environment, source_repository=str(repository), commit=commit,
                prefix="leopard2-adapter-source/", scratch_root=str(scratch),
                destination_name="first.tar",
                child_environment=child_environment, log=command_log)
            git_dir = repository / ".git"
            with (git_dir / "config").open("ab") as stream:
                stream.write(b"\n[tar]\n\tumask = 0077\n")
            (git_dir / "info" / "attributes").write_bytes(b"* export-ignore\n")
            second = acquire.canonical_git_archive(
                environment, source_repository=str(repository), commit=commit,
                prefix="leopard2-adapter-source/", scratch_root=str(scratch),
                destination_name="second.tar",
                child_environment=child_environment, log=command_log)
            self.assertEqual(
                (first["size"], first["sha256"]),
                (second["size"], second["sha256"]))
            self.assertTrue(acquire._owned_files_identical(
                first["path"], second["path"], "fixture archive"))
            self.assertFalse(any(
                path.name.startswith(".leopard-canonical-git-archive.")
                for path in scratch.iterdir()))

    def test_source_and_archive_adversarial_failures_are_closed(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-exact-main-source-negative.") as temporary:
            parent = Path(temporary).resolve()
            repository, commit, tree, submodule, submodule_commit = \
                tiny_repository(parent)
            child_environment = acquire.frozen_child_environment()

            nonempty = parent / "nonempty"
            nonempty.mkdir(mode=0o700)
            (nonempty / "unexpected").write_bytes(b"x")
            self.assertAcquireError(
                acquire.stage_detached_source,
                acquire.HostEnvironment(anchor_path=str(parent)),
                source_repository=str(repository),
                submodule_repository=str(submodule),
                destination=str(nonempty), commit=commit, tree=tree,
                submodule_commit=submodule_commit,
                child_environment=child_environment, log=[])

            wrong_tree = parent / "wrong-tree"
            wrong_tree.mkdir(mode=0o700)
            self.assertAcquireError(
                acquire.stage_detached_source,
                acquire.HostEnvironment(anchor_path=str(parent)),
                source_repository=str(repository),
                submodule_repository=str(submodule),
                destination=str(wrong_tree), commit=commit, tree="0" * 40,
                submodule_commit=submodule_commit,
                child_environment=child_environment, log=[])

            class RewritingEnvironment(acquire.HostEnvironment):
                def __init__(self, mode):
                    super().__init__(anchor_path=str(parent))
                    self.mode = mode

                def run(self, argv, **keywords):
                    result = super().run(argv, **keywords)
                    if (self.mode == "submodule-nonempty" and
                            "checkout" in argv and
                            argv[2] == str(parent / "submodule-nonempty")):
                        injected = Path(argv[2]) / "sse2neon" / "injected"
                        injected.write_bytes(b"injected\n")
                    if (self.mode == "attached" and
                            list(argv[-3:]) == ["symbolic-ref", "-q", "HEAD"] and
                            argv[2] == str(parent / "attached")):
                        return acquire.CommandResult(
                            exit_status=0, stdout=b"refs/heads/main\n",
                            stderr=b"")
                    if (self.mode == "dirty" and "status" in argv and
                            argv[2] == str(parent / "dirty")):
                        return acquire.CommandResult(
                            exit_status=0, stdout=b"?? injected\n", stderr=b"")
                    return result

            for mode in ("attached", "dirty", "submodule-nonempty"):
                destination = parent / mode
                destination.mkdir(mode=0o700)
                with self.subTest(mode=mode):
                    self.assertAcquireError(
                        acquire.stage_detached_source,
                        RewritingEnvironment(mode),
                        source_repository=str(repository),
                        submodule_repository=str(submodule),
                        destination=str(destination), commit=commit, tree=tree,
                        submodule_commit=submodule_commit,
                        child_environment=child_environment, log=[])

            scratch = parent / "scratch"
            scratch.mkdir(mode=0o700)

            class ArchiveEnvironment(acquire.HostEnvironment):
                def __init__(self, mode):
                    super().__init__(anchor_path=str(parent))
                    self.mode = mode

                def run_to_path(self, argv, **keywords):
                    if self.mode == "exit":
                        return acquire.StreamedCommandResult(
                            exit_status=19, stdout_size=0,
                            stdout_sha256=sha256(b""), stderr=b"fatal\n")
                    result = super().run_to_path(argv, **keywords)
                    if self.mode == "stderr":
                        return acquire.StreamedCommandResult(
                            exit_status=result.exit_status,
                            stdout_size=result.stdout_size,
                            stdout_sha256=result.stdout_sha256,
                            stderr=b"unexpected\n")
                    return result

            for mode in ("exit", "stderr"):
                with self.subTest(archive_mode=mode):
                    self.assertAcquireError(
                        acquire.canonical_git_archive,
                        ArchiveEnvironment(mode),
                        source_repository=str(repository), commit=commit,
                        prefix="leopard2-adapter-source/",
                        scratch_root=str(scratch),
                        destination_name=f"{mode}.tar",
                        child_environment=child_environment, log=[])
            self.assertAcquireError(
                acquire.canonical_git_archive,
                acquire.HostEnvironment(anchor_path=str(parent)),
                source_repository=str(repository), commit=commit,
                prefix="wrong-prefix/", scratch_root=str(scratch),
                destination_name="bad-prefix.tar",
                child_environment=child_environment, log=[])
            collision = scratch / "collision.tar"
            collision.write_bytes(b"existing\n")
            collision.chmod(0o600)
            self.assertAcquireError(
                acquire.canonical_git_archive,
                acquire.HostEnvironment(anchor_path=str(parent)),
                source_repository=str(repository), commit=commit,
                prefix="leopard2-adapter-source/",
                scratch_root=str(scratch),
                destination_name="collision.tar",
                child_environment=child_environment, log=[])

    def test_adapter_inventory_binds_capture_blobs(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-exact-main-adapter.") as temporary:
            root = Path(temporary).resolve()
            records = []
            contents = {}
            for index, relative in enumerate(record_contract.ADAPTER_PATHS):
                content = f"adapter {index}\n".encode("ascii")
                path = root / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(content)
                contents[relative] = content
                records.append({
                    "path": relative, "kind": "regular",
                    "git_mode": "100644", "git_type": "blob",
                    "object_id": hashlib.sha1(
                        f"blob {len(content)}\0".encode("ascii") + content,
                        usedforsecurity=False).hexdigest(),
                })
            adapter, controller = acquire.adapter_inventory(
                str(root), {"tracked_files": records})
            self.assertEqual(
                adapter["files_sha256"],
                identity_contract.canonical_sha256(adapter["files"]))
            self.assertEqual(
                controller, contents[record_contract.ADAPTER_PATHS[2]])
            changed = copy.deepcopy(records)
            changed[0]["object_id"] = "0" * 40
            self.assertAcquireError(
                acquire.adapter_inventory, str(root),
                {"tracked_files": changed})

    def test_real_toolchain_capture_has_exact_roles_and_versions(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-exact-main-toolchain.") as temporary:
            root = Path(temporary).resolve()
            plan = temporary_lane_plan(root)
            environment = acquire.HostEnvironment(anchor_path=str(root))
            command_log = []
            toolchain, retained = acquire.resolve_toolchain(
                environment, plan,
                child_environment=acquire.frozen_child_environment(),
                log=command_log)
            self.assertEqual(
                tuple(item["role"] for item in toolchain["tools"]),
                record_contract.TOOL_ROLES)
            self.assertEqual(
                tuple(item["role"] for item in toolchain["subtools"]),
                record_contract.SUBTOOL_ROLES)
            self.assertEqual(
                tuple(item["role"] for item in toolchain["versions"]),
                record_contract.VERSION_ROLES)
            self.assertEqual(len(retained), 2 * len(
                record_contract.VERSION_ROLES))
            self.assertTrue(all(
                version["argv"][0] in {
                    item["resolved_path"] for item in
                    toolchain["tools"] + toolchain["subtools"]}
                for version in toolchain["versions"]))
            git_tool = next(
                item for item in toolchain["tools"]
                if item["role"] == "git")
            capture = {"git_executable": {"source": {
                key: git_tool[key]
                for key in ("path", "size", "mode", "sha256")
            }}}
            capture["git_executable"]["source"]["path"] = \
                git_tool["resolved_path"]
            acquire._capture_git_tool_binding([capture], toolchain)
            changed = copy.deepcopy(capture)
            changed["git_executable"]["source"]["sha256"] = "0" * 64
            self.assertAcquireError(
                acquire._capture_git_tool_binding, [changed], toolchain)

    def test_complete_source_stage_uses_real_offline_git_and_no_timing(self) \
            -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-exact-main-stage.") as temporary:
            parent = Path(temporary).resolve()
            repository, commit, tree, _submodule, submodule_commit = \
                tiny_repository(parent)
            work = parent / "work"
            work.mkdir(mode=0o700)
            plan = acquire.LanePlan(
                lane_root=str(work / "sealed-lane"), attempt=1,
                repository=str(repository),
                verifier=str(repository / "tools/verifier.py"),
                canonical_adapter_root=str(work / "canonical-adapter"),
                canonical_baseline_root=str(work / "canonical-baseline"),
                canonical_build_root=str(work / "canonical-build"),
                variant_adapter_root=str(work / "variant-adapter"),
                variant_baseline_root=str(work / "variant-baseline"),
                variant_build_root=str(work / "variant-build"),
                scratch_root=str(work / "scratch"),
            )

            class CaptureEnvironment(acquire.HostEnvironment):
                def capture_git_identity(
                        self, controller_repository, source_root,
                        requested_commit, *, require_detached):
                    self.assertions.append((
                        controller_repository, source_root,
                        requested_commit, require_detached))
                    observed_commit = git_command(
                        Path(source_root), "rev-parse", "HEAD^{commit}"). \
                        decode("ascii").strip()
                    observed_tree = git_command(
                        Path(source_root), "rev-parse", "HEAD^{tree}"). \
                        decode("ascii").strip()
                    records = []
                    for relative in record_contract.ADAPTER_PATHS:
                        content = (Path(source_root) / relative).read_bytes()
                        records.append({
                            "path": relative, "kind": "regular",
                            "git_mode": "100644", "git_type": "blob",
                            "object_id": hashlib.sha1(
                                f"blob {len(content)}\0".encode("ascii") +
                                content, usedforsecurity=False).hexdigest(),
                        })
                    git_tool = acquire._tool_record("git", "/usr/bin/git")
                    return {
                        "schema": "scripted-source-capture/v1",
                        "head": observed_commit, "tree": observed_tree,
                        "path": source_root, "detached": True,
                        "tracked_status": "clean", "tracked_files": records,
                        "submodules": [{
                            "path": "sse2neon",
                            "object_id": submodule_commit,
                        }],
                        "git_executable": {"source": {
                            "path": git_tool["resolved_path"],
                            "size": git_tool["size"],
                            "mode": git_tool["mode"],
                            "sha256": git_tool["sha256"],
                        }},
                    }

                def __init__(self):
                    super().__init__(anchor_path=str(parent))
                    self.assertions = []

            environment = CaptureEnvironment()
            with mock.patch.object(
                    record_contract, "BASELINE_COMMIT", commit), \
                    mock.patch.object(
                        record_contract, "BASELINE_TREE", tree), \
                    mock.patch.object(
                        record_contract, "BASELINE_SSE2NEON_COMMIT",
                        submodule_commit), \
                    acquire.prepare_acquisition_roots(plan) as prepared:
                result = acquire.acquire_source_stage(
                    environment, plan, prepared)
                self.assertEqual(
                    result.source["baseline"]["commit"], commit)
                self.assertEqual(
                    result.source["adapter_repository"]["commit"], commit)
                self.assertEqual(len(environment.assertions), 2)
                self.assertEqual(
                    set(result.retained_paths), {
                        "source/leopard-main-source.tar",
                        "source/leopard2-adapter-source.tar",
                    })
                self.assertTrue(all(
                    Path(path).is_file()
                    for path in result.retained_paths.values()))
                self.assertIn(
                    "controllers/test_legacy_main_benchmark.py",
                    result.retained_bytes)
                self.assertNotIn(b"timing", result.log.lower())
                self.assertNotIn(b"duration", result.log.lower())

    def test_source_failure_helper_seals_verifiable_lane(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-exact-main-source-failure.") as temporary:
            parent = Path(temporary).resolve()
            plan = temporary_lane_plan(parent)
            environment = acquire.HostEnvironment(anchor_path=str(parent))
            seal = acquire.seal_source_acquisition_failure(
                environment, plan,
                acquire.CommandExecutionError("fixture clone failed", 17),
                diagnostics={"diagnostics/git.stderr": b"fatal fixture\n"})
            self.assertEqual(seal["terminal"], "FAILED.json")
            completed = subprocess.run(
                verifier_fixtures.verifier_cli_command(Path(plan.lane_root)),
                check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                timeout=30)
            self.assertEqual(
                completed.returncode, 3,
                completed.stderr.decode("utf-8", "replace"))
            failure = record_contract.load_baseline_failure_record(
                (Path(plan.lane_root) / "FAILED.json").read_bytes())
            self.assertEqual(failure["stage"], "source_acquisition")
            self.assertEqual(failure["error"]["exit_status"], 17)

    def test_source_stage_failure_seals_its_captured_transcript(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-exact-main-source-transcript.") as temporary:
            parent = Path(temporary).resolve()
            plan = temporary_lane_plan(parent)

            class FailingEnvironment(acquire.HostEnvironment):
                def run(self, argv, **keywords):
                    return acquire.CommandResult(
                        exit_status=23, stdout=b"", stderr=b"injected\n")

                def now_utc(self):
                    return "2026-08-30T00:00:00Z"

            environment = FailingEnvironment(anchor_path=str(parent))
            with acquire.prepare_acquisition_roots(plan) as prepared:
                with self.assertRaises(acquire.SourceStageError) as raised:
                    acquire.acquire_source_stage(environment, plan, prepared)
                failure_error = raised.exception
                transcript = identity_contract.strict_json_loads(
                    failure_error.log, "source-stage failure transcript")
                self.assertEqual(transcript["status"], "failed")
                self.assertEqual(transcript["command_count"], 1)
                self.assertEqual(
                    transcript["commands"][0]["exit_status"], 23)
                seal = acquire.seal_source_acquisition_failure(
                    environment, plan, failure_error)
            self.assertEqual(seal["terminal"], "FAILED.json")
            retained_log = identity_contract.strict_json_loads(
                (Path(plan.lane_root) /
                 "logs/00-source_acquisition.log").read_bytes(),
                "sealed source-stage failure transcript")
            self.assertEqual(retained_log, transcript)
            failure = record_contract.load_baseline_failure_record(
                (Path(plan.lane_root) / "FAILED.json").read_bytes())
            self.assertEqual(failure["error"], {
                "kind": "command_error",
                "message": "adapter repository commit query failed with "
                           "status 23: injected",
                "exit_status": 23,
            })
            completed = subprocess.run(
                verifier_fixtures.verifier_cli_command(Path(plan.lane_root)),
                check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                timeout=30)
            self.assertEqual(completed.returncode, 3)

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

    def test_path_backed_retained_file_is_streamed_into_the_seal(self) -> None:
        with fixture_anchors(), tempfile.TemporaryDirectory(
                prefix="leopard-exact-main-path-backed.") as temporary:
            parent = Path(temporary).resolve()
            root = parent / "lane"
            retained, record = fixture_for_root(
                root, "acquisition_failure")
            relative = sorted(retained)[0]
            source = parent / "path-backed.bin"
            source.write_bytes(retained.pop(relative))
            source.chmod(0o600)
            with acquire.LaneWriter(str(root)) as writer:
                seal = writer.seal_record(
                    record, retained,
                    retained_paths={relative: str(source)})
            self.assertEqual(seal["terminal"], "FAILED.json")
            self.assertEqual(
                (root / relative).read_bytes(), source.read_bytes())
            completed = subprocess.run(
                verifier_fixtures.verifier_cli_command(root),
                check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                timeout=30)
            self.assertEqual(completed.returncode, 3)

    def test_path_backed_publication_rejects_unsafe_sources(self) -> None:
        for kind in ("symlink", "mode", "hardlink", "empty",
                     "symlink-parent", "overlap", "extra"):
            with self.subTest(kind=kind), fixture_anchors(), \
                    tempfile.TemporaryDirectory(
                        prefix="leopard-exact-main-path-negative.") \
                    as temporary:
                parent = Path(temporary).resolve()
                root = parent / "lane"
                retained, record = fixture_for_root(
                    root, "acquisition_failure")
                relative = sorted(retained)[0]
                content = retained[relative]
                source = parent / "source.bin"
                path_key = relative
                if kind != "overlap":
                    retained.pop(relative)
                if kind == "symlink-parent":
                    real_parent = parent / "real-parent"
                    real_parent.mkdir(mode=0o700)
                    real_source = real_parent / "source.bin"
                    real_source.write_bytes(content)
                    real_source.chmod(0o600)
                    linked_parent = parent / "linked-parent"
                    linked_parent.symlink_to(real_parent, target_is_directory=True)
                    source = linked_parent / "source.bin"
                else:
                    source.write_bytes(b"" if kind == "empty" else content)
                    source.chmod(0o644 if kind == "mode" else 0o600)
                if kind == "symlink":
                    target = parent / "target.bin"
                    source.rename(target)
                    source.symlink_to(target)
                elif kind == "hardlink":
                    os.link(source, parent / "second-name.bin")
                elif kind == "extra":
                    path_key = "diagnostics/unexpected.bin"
                with acquire.LaneWriter(str(root)) as writer:
                    self.assertAcquireError(
                        writer.seal_record, record, retained,
                        retained_paths={path_key: str(source)})

    def test_path_backed_publication_rebinds_after_source_swap(self) -> None:
        with fixture_anchors(), tempfile.TemporaryDirectory(
                prefix="leopard-exact-main-path-swap.") as temporary:
            parent = Path(temporary).resolve()
            root = parent / "lane"
            retained, record = fixture_for_root(
                root, "acquisition_failure")
            relative = sorted(retained)[0]
            content = retained.pop(relative)
            source = parent / "source.bin"
            source.write_bytes(content)
            source.chmod(0o600)
            replacement = parent / "replacement.bin"
            replacement.write_bytes(b"X" * len(content))
            replacement.chmod(0o600)
            with acquire.LaneWriter(str(root)) as writer:
                publish_path = writer.publish_path

                def replace_then_publish(path, source_path):
                    os.replace(replacement, source)
                    return publish_path(path, source_path)

                with mock.patch.object(
                        writer, "publish_path",
                        side_effect=replace_then_publish):
                    self.assertAcquireError(
                        writer.seal_record, record, retained,
                        retained_paths={relative: str(source)})

    def test_acquisition_error_failure_kind_verifies_offline(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-exact-main-acquisition-error.") as temporary:
            parent = Path(temporary).resolve()
            plan = temporary_lane_plan(parent)
            environment = acquire.HostEnvironment(anchor_path=str(parent))
            seal = acquire.seal_source_acquisition_failure(
                environment, plan,
                acquire.AcquisitionError("fixture source identity drift"))
            self.assertEqual(seal["terminal"], "FAILED.json")
            failure = record_contract.load_baseline_failure_record(
                (Path(plan.lane_root) / "FAILED.json").read_bytes())
            self.assertEqual(failure["error"], {
                "kind": "acquisition_error",
                "message": "fixture source identity drift",
                "exit_status": 1,
            })
            transcript = identity_contract.strict_json_loads(
                (Path(plan.lane_root) /
                 "logs/00-source_acquisition.log").read_bytes(),
                "acquisition-error transcript")
            self.assertEqual(set(transcript), {
                "schema", "status", "command_count", "commands",
                "adapter_commit", "adapter_tree", "retained_byte_paths",
                "retained_path_sources", "error",
            })
            self.assertEqual(transcript["command_count"], 0)
            completed = subprocess.run(
                verifier_fixtures.verifier_cli_command(Path(plan.lane_root)),
                check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                timeout=30)
            self.assertEqual(completed.returncode, 3)

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
        repository_calls = [
            node for node in ast.walk(tree)
            if isinstance(node, ast.Call) and
            isinstance(node.func, ast.Name) and
            node.func.id == "_load_repository_module"
        ]
        self.assertEqual(len(repository_calls), 1)
        self.assertEqual(len(repository_calls[0].args), 3)
        self.assertEqual(repository_calls[0].keywords, [])
        self.assertEqual(
            [argument.value for argument in repository_calls[0].args[1:]],
            ["_leopard2_exact_main_git_capture",
             "experiments/leopard2/main_compare/git_capture.py"])
        self.assertFalse(any(isinstance(node, ast.Assert)
                             for node in ast.walk(tree)))


if __name__ == "__main__":
    unittest.main()
