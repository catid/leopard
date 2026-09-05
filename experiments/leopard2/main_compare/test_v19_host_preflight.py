#!/usr/bin/python3
"""Synthetic and filesystem-only tests for the dormant v19 pre-lane gate."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("v19_host_preflight_tested", HERE / "v19_host_preflight.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
SCOPE_MEMBER = "/user.slice/user-1000.slice/test.scope"
SCOPE = module.CGROUP_MOUNT + SCOPE_MEMBER
MOUNTINFO = "29 23 0:26 / /sys/fs/cgroup rw,nosuid,nodev,noexec,relatime - cgroup2 cgroup rw,nsdelegate\n"


def fake_identity(path: str) -> dict:
    return {"device": 27, "inode": 1 + int(hashlib.sha256(path.encode()).hexdigest()[:12], 16)}


class FakeReader:
    def __init__(self) -> None:
        self.file_identities = {}
        self.calls = 0
        self.second_process = None
        self.second_files = {}
        self.identity_drift = False
        self.directory_drift = False
        self.directory_reads = 0
        self.process_value = {
            "pid": 123, "uid": 1000, "euid": 1000, "system": "Linux",
            "hostname": "ripper", "machine": "x86_64", "kernel_release": "6.8.0-test",
            "clock_ticks_per_second": 100, "allowed_cpus": list(range(128)),
            "nofile_soft": 65536, "nofile_hard": 1048576,
            "namespaces": {name: fake_identity(name) for name in ("pid", "user", "mnt", "cgroup", "uts")},
        }
        self.files = {
            "/proc/123/cgroup": "0::" + SCOPE_MEMBER + "\n",
            "/proc/123/mountinfo": MOUNTINFO,
            "/proc/swaps": "Filename Type Size Used Priority\n",
            "/proc/sys/kernel/random/boot_id": "01234567-89ab-cdef-0123-456789abcdef\n",
            "/proc/cpuinfo": "\n\n".join(
                f"processor : {cpu}\nmodel name : {module.CPU_MODEL}\ncpu MHz : 1.23"
                for cpu in range(128)) + "\n",
            module.SYS_CPU + "/online": "0-127\n",
            module.SYS_NODE + "/online": "0\n",
            module.SYS_NODE + "/node0/cpulist": "0-127\n",
            SCOPE + "/cgroup.type": "domain\n",
            SCOPE + "/memory.max": str(module.MEMORY_MAX) + "\n",
            SCOPE + "/memory.swap.max": "0\n",
            SCOPE + "/memory.swap.current": "0\n",
            SCOPE + "/memory.events": "low 0\nhigh 0\nmax 0\noom 0\noom_kill 0\noom_group_kill 0\n",
            module.CGROUP_MOUNT + "/user.slice/memory.max": "max\n",
            module.CGROUP_MOUNT + "/user.slice/user-1000.slice/memory.max": "max\n",
        }
        for cpu in range(128):
            root = f"{module.SYS_CPU}/cpu{cpu}"
            primary = cpu % 64
            group = primary // 8
            domain = list(range(group * 8, group * 8 + 8)) + list(range(group * 8 + 64, group * 8 + 72))
            self.files.update({
                root + "/topology/physical_package_id": "0\n",
                root + "/topology/die_id": "0\n",
                root + "/topology/core_id": f"{primary}\n",
                root + "/topology/thread_siblings_list": f"{primary},{primary + 64}\n",
                root + "/cache/index0/level": "1\n",
                root + "/cache/index3/level": "3\n",
                root + "/cache/index3/type": "Unified\n",
                root + "/cache/index3/shared_cpu_list": ",".join(map(str, domain)) + "\n",
            })

    def process(self) -> dict:
        self.calls += 1
        return copy.deepcopy(self.second_process if self.calls == 2 and self.second_process is not None
                             else self.process_value)

    def read_text(self, path: str, limit: int = 4096) -> str:
        files = dict(self.files, **self.second_files) if self.calls == 2 else self.files
        module.require(path in files, "missing fixture kernel file")
        data = files[path]
        module.require(len(data.encode("ascii")) <= limit, "fixture read exceeds bound")
        self.file_identities[path] = fake_identity(path)
        if self.identity_drift and self.calls == 2:
            self.file_identities[path]["inode"] += 1
        return data

    def directory_identity(self, path: str) -> dict:
        value = fake_identity(path)
        if path == SCOPE:
            self.directory_reads += 1
            if self.directory_drift and self.directory_reads > 1:
                value["inode"] += 1
        return value

    def list_cache_indices(self, path: str) -> list[str]:
        return ["index0", "index3"]


class HostPreflightTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        result = subprocess.run([str(HERE / "run_authoritative_v17_gfni_main_compare.sh"),
                                 "--print-conditioned-v19-preregistration"],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        cls.preregistration = result.stdout
        cls.record = module.capture_host(cls.preregistration, FakeReader())

    def test_capture_is_dormant_and_canonical(self) -> None:
        with mock.patch.object(module.os, "mkdir", side_effect=AssertionError("write")), \
             mock.patch.object(module.os, "system", side_effect=AssertionError("command")), \
             mock.patch.object(module.os, "sched_setaffinity", side_effect=AssertionError("affinity")), \
             mock.patch.object(module.resource, "setrlimit", side_effect=AssertionError("limits")):
            record = module.capture_host(self.preregistration, FakeReader())
        self.assertEqual(module.validate_capture(record, self.preregistration), record)
        self.assertEqual(record["observations"][0], record["observations"][1])
        self.assertFalse(record["live_acquisition_armed"])
        self.assertFalse(record["atomic_snapshot"])
        self.assertFalse(record["resource_lifetime_proved"])
        self.assertLess(len(module.canonical_bytes(record)), module.MAX_RECORD_BYTES)

    def test_preregistration_cannot_be_replaced_or_armed(self) -> None:
        for data in (self.preregistration + b"\n", b"{}\n",
                     self.preregistration.replace(b'"live_acquisition_armed":false',
                                                 b'"live_acquisition_armed":true')):
            with self.subTest(data=data[:20]), self.assertRaises(module.PreflightError):
                module.capture_host(data, FakeReader())

    def test_process_contract_rejects_host_affinity_and_limits(self) -> None:
        for key, replacement in (("hostname", "foureyes"), ("machine", "aarch64"),
                                  ("system", "Darwin"), ("clock_ticks_per_second", 100.0),
                                  ("allowed_cpus", list(range(127))), ("nofile_soft", 1024),
                                  ("nofile_hard", 100), ("euid", 0), ("pid", True)):
            reader = FakeReader()
            reader.process_value[key] = replacement
            with self.subTest(key=key), self.assertRaises(module.PreflightError):
                module.capture_host(self.preregistration, reader)

    def test_cgroup_paths_mounts_and_events_are_closed_world(self) -> None:
        for membership, mount in (
                ("0::/\n", MOUNTINFO), ("0::relative\n", MOUNTINFO),
                ("0::/x/../test.scope\n", MOUNTINFO), ("0:://test.scope\n", MOUNTINFO),
                ("0::/test.scope/\n", MOUNTINFO), ("0::/test.scope\n0::/other\n", MOUNTINFO),
                ("7:memory:/test.scope\n", MOUNTINFO),
                ("0::" + SCOPE_MEMBER, MOUNTINFO + MOUNTINFO),
                ("0::" + SCOPE_MEMBER, MOUNTINFO.replace(" / /sys", " /rebased /sys")),
                ("0::" + SCOPE_MEMBER, MOUNTINFO.replace("/sys/fs/cgroup", "/tmp/fake")),
                ("0::" + SCOPE_MEMBER, MOUNTINFO.replace("29 23", "x 23")),
                ("0::" + SCOPE_MEMBER, MOUNTINFO.replace("rw,nsdelegate", "rw,memory_localevents"))):
            with self.subTest(membership=membership, mount=mount), self.assertRaises(module.PreflightError):
                module.cgroup_location(membership, mount)
        for text in ("max 0\nmax 0\n", "high 0\n", "low 0\nhigh 0\nmax 0\noom 1\noom_kill 0\n",
                     "low 0\nhigh 0\nmax -1\noom 0\noom_kill 0\n"):
            with self.subTest(events=text), self.assertRaises(module.PreflightError):
                module.memory_events(text)

    def test_collector_rejects_resources_swap_and_missing_topology(self) -> None:
        for path, replacement in (
                (SCOPE + "/memory.max", "max\n"), (SCOPE + "/memory.max", "1073741824\n"),
                (SCOPE + "/memory.swap.max", "1\n"), (SCOPE + "/memory.swap.current", "1\n"),
                (SCOPE + "/cgroup.type", "threaded\n"),
                (module.CGROUP_MOUNT + "/user.slice/memory.max", "268435456\n"),
                ("/proc/swaps", "Filename Type Size Used Priority\n/swapfile file 100 0 -2\n"),
                (module.SYS_CPU + "/online", "0-126\n"), (module.SYS_NODE + "/online", "0-1\n"),
                (module.SYS_NODE + "/node0/cpulist", "0-63\n"),
                (module.SYS_CPU + "/cpu0/topology/thread_siblings_list", "0,65\n"),
                (module.SYS_CPU + "/cpu0/cache/index3/type", "Data\n"),
                (module.SYS_CPU + "/cpu0/cache/index3/level", "2\n")):
            reader = FakeReader()
            reader.files[path] = replacement
            with self.subTest(path=path, replacement=replacement), self.assertRaises(module.PreflightError):
                module.capture_host(self.preregistration, reader)

    def test_normalized_cpuinfo_ignores_only_dynamic_fields(self) -> None:
        reader = FakeReader()
        reader.second_files["/proc/cpuinfo"] = reader.files["/proc/cpuinfo"].replace("1.23", "4.56")
        module.capture_host(self.preregistration, reader)
        for text in (reader.files["/proc/cpuinfo"].replace(module.CPU_MODEL, "Other CPU", 1),
                     reader.files["/proc/cpuinfo"].replace("processor : 0", "processor : 1", 1)):
            bad = FakeReader()
            bad.files["/proc/cpuinfo"] = text
            with self.assertRaises(module.PreflightError):
                module.capture_host(self.preregistration, bad)

    def test_capture_detects_observed_identity_drift(self) -> None:
        readers = []
        for key in ("kernel_release", "namespaces"):
            reader = FakeReader()
            reader.second_process = copy.deepcopy(reader.process_value)
            if key == "namespaces":
                reader.second_process[key]["cgroup"]["inode"] += 1
            else:
                reader.second_process[key] = "6.8.0-different"
            readers.append(reader)
        reader = FakeReader()
        reader.second_files["/proc/sys/kernel/random/boot_id"] = "01234567-89ab-cdef-0123-456789abcdee\n"
        readers.append(reader)
        for key in ("identity_drift", "directory_drift"):
            reader = FakeReader()
            setattr(reader, key, True)
            readers.append(reader)
        for index, reader in enumerate(readers):
            with self.subTest(case=index), self.assertRaises(module.PreflightError):
                module.capture_host(self.preregistration, reader)

    def test_coherently_resigned_records_cannot_weaken_contract(self) -> None:
        for name in ("affinity boolean", "siblings float", "L3 cross die", "core alias",
                     "ancestor missing", "ancestor boolean", "memory boolean", "events boolean",
                     "different observations", "false authority claim"):
            record = copy.deepcopy(self.record)
            for observation in record["observations"]:
                topology, resources = observation["topology"], observation["resources"]
                if name == "affinity boolean": observation["process"]["allowed_cpus"][0] = False
                if name == "siblings float": topology[0]["siblings"][0] = 0.0
                if name == "L3 cross die": topology[0]["die"] = topology[64]["die"] = 1
                if name == "core alias": topology[0]["core"] = topology[64]["core"] = 1
                if name == "ancestor missing": resources["ancestors"].pop()
                if name == "ancestor boolean": resources["ancestors"][0]["memory_max"] = True
                if name == "memory boolean": resources["memory_swap_max"] = False
                if name == "events boolean": resources["memory_events"]["max"] = False
            if name == "different observations": record["observations"][1]["process"]["namespaces"]["mnt"]["inode"] += 1
            if name == "false authority claim": record["resource_lifetime_proved"] = True
            record.pop("digest")
            record["digest"] = hashlib.sha256(module.canonical_bytes(record)).hexdigest()
            with self.subTest(splice=name), self.assertRaises(module.PreflightError):
                module.validate_capture(record, self.preregistration)

    def test_native_reader_bounds_regular_files_and_rejects_symlinks(self) -> None:
        reader = module.LinuxReader()
        with tempfile.TemporaryDirectory(prefix="leopard-v19-host-reader-") as temporary:
            root = Path(temporary)
            (root / "value").write_bytes(b"1234")
            self.assertEqual(reader.read_text(str(root / "value"), 4), "1234")
            with self.assertRaises(module.PreflightError): reader.read_text(str(root / "value"), 3)
            (root / "link").symlink_to(root / "value")
            (root / "dirlink").symlink_to(root, target_is_directory=True)
            os.mkfifo(root / "fifo")
            for path in (root / "link", root / "dirlink/value", root / "fifo"):
                with self.subTest(path=path), self.assertRaises((OSError, module.PreflightError)):
                    reader.read_text(str(path), 4096)
            for text in ("../x", "/tmp/../x", "/tmp//x", "/tmp/x/", "/tmp/\x00x"):
                with self.subTest(text=text), self.assertRaises(module.PreflightError):
                    module.canonical_path(text)

    def test_native_reader_detects_file_replacement(self) -> None:
        reader = module.LinuxReader()
        with tempfile.TemporaryDirectory(prefix="leopard-v19-host-replace-") as temporary:
            root = Path(temporary)
            value, replacement = root / "value", root / "replacement"
            value.write_bytes(b"one")
            replacement.write_bytes(b"two")
            real_read = os.read
            swapped = False
            def swapping_read(fd: int, count: int) -> bytes:
                nonlocal swapped
                result = real_read(fd, count)
                if not swapped:
                    replacement.replace(value)
                    swapped = True
                return result
            with mock.patch.object(module.os, "read", side_effect=swapping_read), \
                 self.assertRaisesRegex(module.PreflightError, "identity changed"):
                reader.read_text(str(value))

    def test_record_shapes_are_bounded_before_json_encoding(self) -> None:
        for case in ("affinity", "membership", "mountinfo", "events"):
            value = copy.deepcopy(self.record)
            observation = value["observations"][0]
            if case == "affinity": observation["process"]["allowed_cpus"].append(128)
            if case == "membership": observation["resources"]["membership"] = "/" + "x" * 4097
            if case == "mountinfo": observation["resources"]["mountinfo"] = "x" * (module.MAX_TEXT + 1)
            if case == "events": observation["resources"]["memory_events"] = {"x" * (i + 1): 0 for i in range(65)}
            with self.subTest(case=case), \
                 mock.patch.object(module, "canonical_bytes", side_effect=AssertionError("encoded unbounded input")), \
                 self.assertRaises(module.PreflightError):
                module.validate_capture(value, self.preregistration)

    def test_native_reader_detects_parent_replacement(self) -> None:
        reader = module.LinuxReader()
        with tempfile.TemporaryDirectory(prefix="leopard-v19-host-parent-") as temporary:
            root = Path(temporary)
            parent = root / "parent"
            parent.mkdir()
            (parent / "value").write_bytes(b"one")
            replacement = root / "replacement"
            replacement.mkdir()
            (replacement / "value").write_bytes(b"two")
            real_read = os.read
            swapped = False
            def swapping_read(fd: int, count: int) -> bytes:
                nonlocal swapped
                result = real_read(fd, count)
                if not swapped:
                    parent.rename(root / "retained-original")
                    replacement.rename(parent)
                    swapped = True
                return result
            with mock.patch.object(module.os, "read", side_effect=swapping_read), \
                 self.assertRaisesRegex(module.PreflightError, "parent identity changed"):
                reader.read_text(str(parent / "value"))


if __name__ == "__main__":
    unittest.main()
