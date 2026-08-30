#!/usr/bin/env python3
"""Deterministic, host-free tests for pair qualification acquisition."""

from __future__ import annotations

import ast
import copy
import hashlib
import importlib.util
import os
from pathlib import Path
import stat
import sys
import tempfile
import unittest
from unittest import mock


HERE = Path(__file__).resolve().parent


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load test dependency {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


acquire = load_module(
    "leopard2_pair_qualification_acquire_test",
    HERE / "pair_qualification_acquire.py")
contract = acquire.contract


def policy_fixture() -> dict:
    return contract.qualification_policy_record(
        clock_ticks_per_second=100,
        candidate_primary_cpus=[2, 3],
        excluded_pairs=[],
        domain_mode="pair-and-domain")


def topology_files(*, core_3: int = 1) -> dict[str, bytes]:
    values = {
        acquire.ONLINE_CPUS_PATH: b"2-3,6-7\n",
    }
    for cpu, core, sibling in ((2, 0, 6), (6, 0, 2),
                               (3, core_3, 7), (7, core_3, 3)):
        root = f"{acquire.SYSFS_CPU_ROOT}/cpu{cpu}"
        values[root + "/topology/physical_package_id"] = b"0\n"
        values[root + "/topology/die_id"] = b"0\n"
        values[root + "/topology/core_id"] = f"{core}\n".encode("ascii")
        values[root + "/topology/thread_siblings_list"] = \
            f"{min(cpu, sibling)},{max(cpu, sibling)}\n".encode("ascii")
        values[root + "/cache/index3/shared_cpu_list"] = b"2-3,6-7\n"
    return values


def topology_files_with_out_of_affinity_candidate() -> dict[str, bytes]:
    values = topology_files()
    values[acquire.ONLINE_CPUS_PATH] = b"2-3,6-7,34,38\n"
    for cpu, sibling in ((34, 38), (38, 34)):
        root = f"{acquire.SYSFS_CPU_ROOT}/cpu{cpu}"
        values[root + "/topology/physical_package_id"] = b"1\n"
        values[root + "/topology/die_id"] = b"0\n"
        values[root + "/topology/core_id"] = b"0\n"
        values[root + "/topology/thread_siblings_list"] = b"34,38\n"
        values[root + "/cache/index3/shared_cpu_list"] = b"34,38\n"
    return values


def proc_stat(snapshot_index: int, *, cpu_order=(0, 2, 3, 6, 7)) -> bytes:
    rows = [b"cpu  0 0 0 0 0 0 0 0 0 0"]
    for cpu in cpu_order:
        idle = 1000 + cpu + snapshot_index * 25
        rows.append(
            f"cpu{cpu} 0 0 0 {idle} 0 0 0 0 0 0".encode("ascii"))
    rows.extend((b"intr 0", b"ctxt 0", b"btime 0", b"processes 0"))
    return b"\n".join(rows) + b"\n"


class FakeHostReader:
    def __init__(
        self, *, proc_records=None, topology_variants=None,
        allowed_values=None, tick_values=None, sleep_scale: int = 1,
    ) -> None:
        self.proc_records = list(proc_records or [
            proc_stat(0), proc_stat(1), proc_stat(2)])
        self.topology_variants = list(topology_variants or [topology_files()])
        self.allowed_values = list(allowed_values or [[2, 3, 6, 7]])
        self.tick_values = list(tick_values or [100])
        self.sleep_scale = sleep_scale
        self.now = 1_000_000_000
        self.allowed_calls = 0
        self.tick_calls = 0
        self.proc_calls = 0
        self.topology_pass = -1
        self.sleep_calls: list[int] = []
        self.paths: list[str] = []

    @staticmethod
    def _value(values, index):
        return copy.deepcopy(values[min(index, len(values) - 1)])

    def allowed_cpus(self):
        value = self._value(self.allowed_values, self.allowed_calls)
        self.allowed_calls += 1
        return value

    def clock_ticks_per_second(self):
        value = self._value(self.tick_values, self.tick_calls)
        self.tick_calls += 1
        return value

    def monotonic_ns(self):
        value = self.now
        self.now += 1_000_000
        return value

    def sleep_ns(self, duration_ns: int):
        self.sleep_calls.append(duration_ns)
        self.now += duration_ns * self.sleep_scale

    def read_bytes(self, path: str, maximum_bytes: int):
        self.paths.append(path)
        if path == acquire.PROC_STAT_PATH:
            if self.proc_calls >= len(self.proc_records):
                raise acquire.AcquisitionError("fake /proc/stat is exhausted")
            value = self.proc_records[self.proc_calls]
            self.proc_calls += 1
            return value
        if path == acquire.ONLINE_CPUS_PATH:
            self.topology_pass += 1
        variant = self.topology_variants[
            min(max(self.topology_pass, 0), len(self.topology_variants) - 1)]
        if path not in variant:
            raise acquire.AcquisitionError(f"unmapped fake host path {path}")
        value = variant[path]
        if len(value) > maximum_bytes:
            raise acquire.AcquisitionError("fake host read exceeded bound")
        return value


def acquisition_fixture(reader=None) -> tuple[dict, dict, FakeHostReader]:
    policy = policy_fixture()
    host = reader or FakeHostReader()
    record = acquire.acquire_pair_qualification(
        host, policy_value=policy, window_count=2,
        nominal_window_ns=250_000_000,
        frozen_pair_from_prior_attempt=None)
    return record, policy, host


class PairQualificationAcquireTest(unittest.TestCase):
    def assertAcquireFailure(self, function, *args, **kwargs) -> None:  # noqa: N802
        with self.assertRaises(contract.QualificationError):
            function(*args, **kwargs)

    def validate(self, record, policy):
        return acquire.validate_acquisition_record(
            record, expected_policy=policy, expected_frozen_pair=None,
            expected_window_count=2,
            expected_nominal_window_ns=250_000_000)

    def test_synthetic_acquisition_has_a_canonical_golden_record(self) -> None:
        record, policy, reader = acquisition_fixture()
        self.assertEqual(self.validate(record, policy), record)
        self.assertEqual(reader.sleep_calls, [250_000_000, 250_000_000])
        self.assertEqual(reader.allowed_calls, 2)
        self.assertEqual(reader.tick_calls, 2)
        self.assertEqual(reader.proc_calls, 3)
        self.assertEqual(record["scan"]["scan_window_count"], 2)
        self.assertEqual(record["scan"]["candidate_pair_count"], 2)
        self.assertEqual(record["scan"]["eligible_pair_count"], 2)
        self.assertEqual(record["scan"]["selected"], {
            "benchmark_cpu": 2, "reserved_sibling": 6,
        })
        self.assertEqual(
            contract.canonical_sha256(record),
            "1865f274bdb4ffa7963e2be0d70dfbd1c5984e2ffcd3a0e0ed13a10139ce7c67")

    def test_cpu_list_parser_is_closed_and_canonical(self) -> None:
        self.assertEqual(acquire.parse_cpu_list(b"0-3,8,10-11\n"),
                         [0, 1, 2, 3, 8, 10, 11])
        self.assertEqual(acquire.format_cpu_list([0, 1, 2, 3, 8, 10, 11]),
                         "0-3,8,10-11")
        for invalid in (
                b"", b"0,1\n", b"01\n", b"2-1\n", b"1,,2\n",
                b"1, 2\n", b"1-2,2-3\n", b"1\n2\n", b"1\r\n",
                b"1048576\n", b"x\n"):
            with self.subTest(invalid=invalid):
                self.assertAcquireFailure(acquire.parse_cpu_list, invalid)

    def test_proc_stat_parser_rejects_malformed_missing_and_reordered_rows(self) -> None:
        expected = [2, 3, 6, 7]
        parsed = acquire.parse_proc_stat(proc_stat(0), expected)
        self.assertEqual(sorted(parsed), expected)
        self.assertEqual(parsed[2]["idle"], 1002)
        mutations = (
            proc_stat(0).replace(b"cpu3 ", b"cpu03 "),
            proc_stat(0).replace(b"cpu3 ", b"cpu2 "),
            proc_stat(0).replace(b"cpu3 0 0 0", b"cpu3 x 0 0"),
            proc_stat(0).replace(b"cpu3 0 0 0 1003 0 0 0 0", b"cpu3 0 0"),
            proc_stat(0).replace(b"cpu7 0 0 0 1007 0 0 0 0 0 0\n", b""),
            proc_stat(0, cpu_order=(0, 3, 2, 6, 7)),
            proc_stat(0)[:-1],
            proc_stat(0) + b"cpu8 0 0 0 bad 0 0 0 0\n",
        )
        for mutation in mutations:
            with self.subTest(digest=hashlib.sha256(mutation).hexdigest()):
                self.assertAcquireFailure(
                    acquire.parse_proc_stat, mutation, expected)

    def test_acquisition_rejects_topology_affinity_clock_and_duration_drift(self) -> None:
        changed = topology_files(core_3=9)
        readers = (
            FakeHostReader(topology_variants=[topology_files(), changed]),
            FakeHostReader(allowed_values=[[2, 3, 6, 7], [2, 6, 7]]),
            FakeHostReader(tick_values=[100, 250]),
            FakeHostReader(sleep_scale=0),
        )
        for reader in readers:
            with self.subTest(reader=reader.__dict__):
                self.assertAcquireFailure(
                    acquire.acquire_pair_qualification,
                    reader, policy_value=policy_fixture(), window_count=2,
                    nominal_window_ns=250_000_000,
                    frozen_pair_from_prior_attempt=None)

    def test_online_policy_primary_outside_affinity_is_retained_then_skipped(self) -> None:
        policy = contract.qualification_policy_record(
            clock_ticks_per_second=100,
            candidate_primary_cpus=[2, 34],
            excluded_pairs=[], domain_mode="pair-and-domain")
        reader = FakeHostReader(
            topology_variants=[topology_files_with_out_of_affinity_candidate()])
        record = acquire.acquire_pair_qualification(
            reader, policy_value=policy, window_count=2,
            nominal_window_ns=250_000_000,
            frozen_pair_from_prior_attempt=None)
        self.assertEqual(record["scan"]["candidate_pair_count"], 1)
        self.assertEqual(record["scan"]["selected"], {
            "benchmark_cpu": 2, "reserved_sibling": 6,
        })
        retained_cpus = {
            entry["cpu"] for entry in record["scan"]["topology_before"]["cpus"]
        }
        self.assertTrue({34, 38} <= retained_cpus)

    def test_acquisition_rejects_missing_extra_and_regressing_counter_data(self) -> None:
        missing = proc_stat(0).replace(
            b"cpu7 0 0 0 1007 0 0 0 0 0 0\n", b"")
        noncanonical = proc_stat(0).replace(b"cpu7 ", b"cpu07 ")
        regressed = [proc_stat(0), proc_stat(1), proc_stat(0)]
        for records in ([missing, proc_stat(1), proc_stat(2)],
                        [noncanonical, proc_stat(1), proc_stat(2)], regressed):
            with self.subTest(records=len(records)):
                self.assertAcquireFailure(
                    acquire.acquire_pair_qualification,
                    FakeHostReader(proc_records=records),
                    policy_value=policy_fixture(), window_count=2,
                    nominal_window_ns=250_000_000,
                    frozen_pair_from_prior_attempt=None)

    def test_record_validator_recomputes_every_bound_view(self) -> None:
        record, policy, unused_reader = acquisition_fixture()
        mutations = []
        extra = copy.deepcopy(record)
        extra["extra"] = None
        mutations.append(extra)
        for key, value in (
                ("host_mutation_performed", True),
                ("candidate_timing_performed", 0),
                ("policy_sha256", "0" * 64),
                ("scan_sha256", "0" * 64),
                ("topology_before_sha256", "0" * 64),
                ("requested_window_count", True),
                ("nominal_window_ns", 250_000_001),
                ("clock_ticks_per_second_after_scan", 101)):
            mutation = copy.deepcopy(record)
            mutation[key] = value
            mutations.append(mutation)
        affinity = copy.deepcopy(record)
        affinity["allowed_cpu_set_after_scan"] = [2, 3, 6]
        mutations.append(affinity)
        scan = copy.deepcopy(record)
        scan["scan"]["candidate_timing_performed"] = True
        scan["scan_sha256"] = contract.canonical_sha256(scan["scan"])
        mutations.append(scan)
        for mutation in mutations:
            with self.subTest(keys=set(mutation)):
                self.assertAcquireFailure(self.validate, mutation, policy)

    def test_strict_record_loader_rejects_duplicate_multivalue_and_nonfinite(self) -> None:
        record, policy, unused_reader = acquisition_fixture()
        canonical = contract.canonical_json_bytes(record)
        loaded = acquire.load_acquisition_record(
            canonical, expected_policy=policy, expected_frozen_pair=None,
            expected_window_count=2,
            expected_nominal_window_ns=250_000_000)
        self.assertEqual(loaded, record)
        for invalid in (
                b'{"schema":1,"schema":2}\n',
                canonical + canonical,
                b'{"value":NaN}\n',
                b'{"value":1e309}\n'):
            with self.subTest(invalid=invalid[:32]):
                self.assertAcquireFailure(
                    acquire.load_acquisition_record, invalid,
                    expected_policy=policy, expected_frozen_pair=None,
                    expected_window_count=2,
                    expected_nominal_window_ns=250_000_000)

    def test_fake_acquisition_calls_no_host_mutator_or_process_launcher(self) -> None:
        forbidden = [
            (os, "sched_setaffinity"), (os, "kill"),
            (os, "setpriority"), (os, "nice"),
        ]
        with mock.patch.object(os, "system", side_effect=AssertionError), \
                mock.patch.object(os, "fork", side_effect=AssertionError), \
                mock.patch.object(os, "execve", side_effect=AssertionError):
            with mock.patch.object(
                    forbidden[0][0], forbidden[0][1],
                    side_effect=AssertionError, create=True), \
                    mock.patch.object(
                        forbidden[1][0], forbidden[1][1],
                        side_effect=AssertionError, create=True), \
                    mock.patch.object(
                        forbidden[2][0], forbidden[2][1],
                        side_effect=AssertionError, create=True), \
                    mock.patch.object(
                        forbidden[3][0], forbidden[3][1],
                        side_effect=AssertionError, create=True), \
                    mock.patch.object(
                        os, "open", side_effect=AssertionError), \
                    mock.patch.object(
                        os, "write", side_effect=AssertionError), \
                    mock.patch.object(
                        os, "unlink", side_effect=AssertionError), \
                    mock.patch.object(
                        os, "fsync", side_effect=AssertionError):
                record, policy, unused_reader = acquisition_fixture()
                self.assertEqual(self.validate(record, policy), record)

        tree = ast.parse((HERE / "pair_qualification_acquire.py").read_text("utf-8"))
        direct_imports = {
            (alias.name, alias.asname) for node in ast.walk(tree)
            if isinstance(node, ast.Import) for alias in node.names
        }
        from_imports = {
            (node.module, tuple((alias.name, alias.asname) for alias in node.names))
            for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)
        }
        os_attributes = {
            node.attr for node in ast.walk(tree)
            if isinstance(node, ast.Attribute) and
            isinstance(node.value, ast.Name) and node.value.id == "os"
        }
        self.assertEqual(direct_imports, {
            ("copy", None), ("importlib.util", None), ("os", None),
            ("re", None), ("stat", None), ("sys", None), ("time", None),
        })
        self.assertEqual(from_imports, {
            ("__future__", (("annotations", None),)),
            ("pathlib", (("Path", None),)),
            ("typing", (("Any", None), ("NoReturn", None),
                        ("Protocol", None), ("Sequence", None))),
        })
        self.assertEqual(os_attributes, {
            "O_CLOEXEC", "O_CREAT", "O_DIRECTORY", "O_EXCL", "O_NOFOLLOW",
            "O_RDONLY", "O_WRONLY", "PathLike", "close", "fstat", "fsync",
            "geteuid", "open", "read", "sched_getaffinity", "sysconf",
            "unlink", "write",
        })
        system_reader = next(
            node for node in tree.body
            if isinstance(node, ast.ClassDef) and node.name == "SystemHostReader")
        system_reader_os_attributes = {
            node.attr for node in ast.walk(system_reader)
            if isinstance(node, ast.Attribute) and
            isinstance(node.value, ast.Name) and node.value.id == "os"
        }
        self.assertEqual(system_reader_os_attributes, {
            "O_CLOEXEC", "O_NOFOLLOW", "O_RDONLY", "close", "fstat",
            "open", "read", "sched_getaffinity", "sysconf",
        })

    def test_exclusive_writer_publishes_0600_and_cleans_partial_output(self) -> None:
        record, policy, unused_reader = acquisition_fixture()
        keywords = {
            "expected_policy": policy,
            "expected_frozen_pair": None,
            "expected_window_count": 2,
            "expected_nominal_window_ns": 250_000_000,
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            os.chmod(root, 0o700)
            output = root / "qualification.json"
            written = acquire.write_acquisition_record_exclusive(
                output, record, **keywords)
            self.assertEqual(written, record)
            self.assertEqual(output.read_bytes(), contract.canonical_json_bytes(record))
            self.assertEqual(stat.S_IMODE(output.stat().st_mode), 0o600)
            self.assertAcquireFailure(
                acquire.write_acquisition_record_exclusive,
                output, record, **keywords)

            partial = root / "partial.json"
            real_write = os.write
            calls = 0

            def failing_write(descriptor, content):
                nonlocal calls
                calls += 1
                if calls == 1:
                    return real_write(descriptor, content[:7])
                raise OSError("injected write failure")

            with mock.patch.object(acquire.os, "write", side_effect=failing_write):
                self.assertAcquireFailure(
                    acquire.write_acquisition_record_exclusive,
                    partial, record, **keywords)
            self.assertFalse(partial.exists())

    def test_system_reader_only_reads_bounded_regular_temporary_files(self) -> None:
        reader = acquire.SystemHostReader()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            source = root / "source"
            source.write_bytes(b"fixture\n")
            self.assertEqual(reader.read_bytes(str(source), 8), b"fixture\n")
            self.assertTrue(source.is_file())
            self.assertAcquireFailure(reader.read_bytes, str(source), 7)
            link = root / "link"
            link.symlink_to(source)
            self.assertAcquireFailure(reader.read_bytes, str(link), 8)
            for path in ("relative", str(root / "child" / ".." / "source")):
                with self.subTest(path=path):
                    self.assertAcquireFailure(reader.read_bytes, path, 8)

    def test_system_reader_is_not_needed_for_deterministic_tests(self) -> None:
        with mock.patch.object(
                acquire.SystemHostReader, "__init__",
                side_effect=AssertionError("live host reader instantiated"),
                create=True):
            record, policy, unused_reader = acquisition_fixture()
            self.assertEqual(self.validate(record, policy), record)


if __name__ == "__main__":
    unittest.main()
