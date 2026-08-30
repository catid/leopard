#!/usr/bin/env python3
"""Read-only Linux acquisition for the pure CPU-pair qualification contract.

The public acquisition boundary is deliberately injectable.  Production reads are
available through :class:`SystemHostReader`, while deterministic tests supply a
fake reader and therefore never inspect or mutate the host.  This module does not
launch benchmark children, change affinity, signal processes, alter priorities, or
write cgroups.  Its only write operation is exclusive publication of an already
validated JSON record.
"""

from __future__ import annotations

import copy
import importlib.util
import os
from pathlib import Path
import re
import stat
import sys
import time
from typing import Any, NoReturn, Protocol, Sequence


def _load_contract() -> Any:
    module_name = "leopard2_pair_qualification_contract_for_acquisition"
    expected = Path(__file__).resolve().with_name("pair_qualification_contract.py")
    loaded = sys.modules.get(module_name)
    if loaded is not None:
        if Path(getattr(loaded, "__file__", "")).resolve() != expected:
            raise RuntimeError("pair qualification contract came from another path")
        return loaded
    specification = importlib.util.spec_from_file_location(module_name, expected)
    if specification is None or specification.loader is None:
        raise RuntimeError("cannot load the pair qualification contract")
    module = importlib.util.module_from_spec(specification)
    sys.modules[module_name] = module
    try:
        specification.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(module_name, None)
        raise
    if Path(getattr(module, "__file__", "")).resolve() != expected:
        raise RuntimeError("pair qualification contract resolved to another path")
    return module


contract = _load_contract()

ACQUISITION_SCHEMA = "leopard2-pair-qualification-acquisition/v1"
ACQUISITION_METHOD = \
    "read-only-linux-sysfs-proc-stat-shared-snapshot-chain/v1"
MAX_NOMINAL_WINDOW_NS = 24 * 60 * 60 * 1_000_000_000
MAX_HOST_FILE_BYTES = 16 * 1024 * 1024
MAX_PROC_STAT_BYTES = 4 * 1024 * 1024
OUTPUT_FILE_MODE = 0o600

ONLINE_CPUS_PATH = "/sys/devices/system/cpu/online"
PROC_STAT_PATH = "/proc/stat"
SYSFS_CPU_ROOT = "/sys/devices/system/cpu"

_SOURCE_ITEMS = (
    ("allowed_cpu_set_at_launch", "sched_getaffinity(0)"),
    ("clock_ticks_per_second", "sysconf(SC_CLK_TCK)"),
    ("counter_source", PROC_STAT_PATH),
    ("monotonic_clock", "time.monotonic_ns"),
    ("online_cpus", ONLINE_CPUS_PATH),
    ("physical_package_id", SYSFS_CPU_ROOT +
     "/cpu{cpu}/topology/physical_package_id"),
    ("die_id", SYSFS_CPU_ROOT + "/cpu{cpu}/topology/die_id"),
    ("core_id", SYSFS_CPU_ROOT + "/cpu{cpu}/topology/core_id"),
    ("thread_siblings", SYSFS_CPU_ROOT +
     "/cpu{cpu}/topology/thread_siblings_list"),
    ("domain_cpus", SYSFS_CPU_ROOT +
     "/cpu{cpu}/cache/index3/shared_cpu_list"),
)

ACQUISITION_KEYS = frozenset((
    "schema", "acquisition_method", "sources", "policy", "policy_sha256",
    "requested_window_count", "nominal_window_ns",
    "frozen_pair_from_prior_attempt", "allowed_cpu_set_at_launch",
    "allowed_cpu_set_after_scan", "clock_ticks_per_second_at_launch",
    "clock_ticks_per_second_after_scan", "topology_before_sha256",
    "topology_after_sha256", "scan", "scan_sha256",
    "host_mutation_performed", "candidate_timing_performed",
))
SOURCE_KEYS = frozenset(key for key, unused_value in _SOURCE_ITEMS)
PAIR_KEYS = frozenset(("benchmark_cpu", "reserved_sibling"))
_CPU_LINE = re.compile(rb"cpu([0-9]+)[ \t]+([0-9]+(?:[ \t]+[0-9]+){7,})[ \t]*")
_DECIMAL = re.compile(rb"0|[1-9][0-9]*")


class AcquisitionError(contract.QualificationError):
    """Host acquisition or publication violated the closed contract."""


def _fail(message: str) -> NoReturn:
    raise AcquisitionError(message)


def _require(condition: bool, message: str) -> None:
    if not condition:
        _fail(message)


def _bounded_int(value: Any, minimum: int, maximum: int, label: str) -> int:
    _require(type(value) is int, f"{label} is not an exact integer")
    _require(minimum <= value <= maximum, f"{label} is outside its structural bound")
    return value


def _sources_record() -> dict[str, str]:
    return dict(_SOURCE_ITEMS)


def _canonical_cpu_list(
    value: Any, *, label: str, allow_empty: bool = False,
) -> list[int]:
    _require(type(value) in (list, tuple), f"{label} is not a CPU sequence")
    _require((allow_empty or len(value) > 0) and
             len(value) <= contract.MAX_CPU_COUNT,
             f"{label} has an invalid length")
    cpus = [
        _bounded_int(cpu, 0, contract.MAX_CPU_ID, f"{label} entry")
        for cpu in value
    ]
    _require(len(cpus) == len(set(cpus)), f"{label} contains a duplicate")
    return sorted(cpus)


def _pair(value: Any, label: str) -> dict[str, int] | None:
    if value is None:
        return None
    _require(type(value) is dict and set(value) == PAIR_KEYS,
             f"{label} is not an exact pair object")
    benchmark = _bounded_int(
        value["benchmark_cpu"], 0, contract.MAX_CPU_ID,
        f"{label} benchmark CPU")
    sibling = _bounded_int(
        value["reserved_sibling"], 0, contract.MAX_CPU_ID,
        f"{label} reserved sibling")
    _require(benchmark != sibling, f"{label} repeats one logical CPU")
    return {"benchmark_cpu": benchmark, "reserved_sibling": sibling}


def format_cpu_list(cpus_value: Any) -> str:
    cpus = _canonical_cpu_list(cpus_value, label="CPU list")
    ranges: list[str] = []
    start = previous = cpus[0]
    for cpu in cpus[1:]:
        if cpu == previous + 1:
            previous = cpu
            continue
        ranges.append(str(start) if start == previous else f"{start}-{previous}")
        start = previous = cpu
    ranges.append(str(start) if start == previous else f"{start}-{previous}")
    return ",".join(ranges)


def parse_cpu_list(data: bytes, label: str = "Linux CPU list") -> list[int]:
    _require(type(data) is bytes and 0 < len(data) <= MAX_HOST_FILE_BYTES,
             f"{label} has an invalid byte length")
    raw = data[:-1] if data.endswith(b"\n") else data
    _require(raw and b"\n" not in raw and b"\r" not in raw,
             f"{label} is not one line")
    try:
        text = raw.decode("ascii")
    except UnicodeDecodeError as error:
        raise AcquisitionError(f"{label} is not ASCII") from error
    cpus: list[int] = []
    for item in text.split(","):
        _require(item and item.strip() == item, f"{label} has empty or spaced terms")
        parts = item.split("-")
        _require(len(parts) in (1, 2) and all(
            len(part) <= 7 and
            (part == "0" or (part.isdigit() and not part.startswith("0")))
            for part in parts), f"{label} has a malformed range")
        start = _bounded_int(int(parts[0]), 0, contract.MAX_CPU_ID, label)
        finish = start if len(parts) == 1 else _bounded_int(
            int(parts[1]), 0, contract.MAX_CPU_ID, label)
        _require(start <= finish, f"{label} has a reversed range")
        _require(finish - start < contract.MAX_CPU_COUNT,
                 f"{label} range is too large")
        cpus.extend(range(start, finish + 1))
        _require(len(cpus) <= contract.MAX_CPU_COUNT, f"{label} is too large")
    canonical = _canonical_cpu_list(cpus, label=label)
    _require(cpus == canonical and text == format_cpu_list(canonical),
             f"{label} is not canonical")
    return canonical


def _parse_decimal(data: bytes, label: str) -> int:
    _require(type(data) is bytes and 0 < len(data) <= 64,
             f"{label} has an invalid byte length")
    raw = data[:-1] if data.endswith(b"\n") else data
    _require(_DECIMAL.fullmatch(raw) is not None,
             f"{label} is not canonical decimal")
    return _bounded_int(int(raw), 0, (1 << 31) - 1, label)


def parse_proc_stat(data: bytes, expected_cpus_value: Any) -> dict[int, dict[str, int]]:
    """Parse ordered Linux per-CPU counters and retain the fixed first 8 fields."""
    expected = _canonical_cpu_list(
        expected_cpus_value, label="expected /proc/stat CPUs", allow_empty=True)
    _require(type(data) is bytes and 0 < len(data) <= MAX_PROC_STAT_BYTES,
             "/proc/stat has an invalid byte length")
    _require(data.endswith(b"\n") and b"\x00" not in data,
             "/proc/stat is truncated or contains NUL")
    observed_order: list[int] = []
    retained: dict[int, dict[str, int]] = {}
    expected_set = set(expected)
    for line in data.splitlines():
        if line.startswith(b"cpu "):
            continue
        if not line.startswith(b"cpu"):
            continue
        match = _CPU_LINE.fullmatch(line)
        _require(match is not None, "/proc/stat contains a malformed CPU row")
        raw_cpu = match.group(1)
        _require(len(raw_cpu) <= 7 and
                 (raw_cpu == b"0" or not raw_cpu.startswith(b"0")),
                 "/proc/stat contains a noncanonical CPU number")
        cpu = _bounded_int(
            int(raw_cpu), 0, contract.MAX_CPU_ID, "/proc/stat CPU")
        _require(not observed_order or cpu > observed_order[-1],
                 "/proc/stat CPU rows are duplicated or reordered")
        observed_order.append(cpu)
        values = match.group(2).split()
        _require(len(values) >= len(contract.COUNTER_FIELDS),
                 "/proc/stat CPU row has too few counters")
        first = values[:len(contract.COUNTER_FIELDS)]
        _require(all(len(value) <= 20 and
                     _DECIMAL.fullmatch(value) is not None for value in first),
                 "/proc/stat CPU counter is not canonical decimal")
        if cpu in expected_set:
            fields = {
                field: _bounded_int(
                    int(raw), 0, contract.MAX_COUNTER,
                    f"/proc/stat CPU {cpu} {field}")
                for field, raw in zip(contract.COUNTER_FIELDS, first)
            }
            retained[cpu] = fields
    _require(set(retained) == expected_set,
             "/proc/stat differs from the observed CPU universe")
    return retained


class HostReader(Protocol):
    """The only host-facing operations used by qualification acquisition."""

    def allowed_cpus(self) -> Sequence[int]: ...

    def clock_ticks_per_second(self) -> int: ...

    def monotonic_ns(self) -> int: ...

    def sleep_ns(self, duration_ns: int) -> None: ...

    def read_bytes(self, path: str, maximum_bytes: int) -> bytes: ...


class SystemHostReader:
    """Bounded, read-only Linux implementation of :class:`HostReader`."""

    def allowed_cpus(self) -> Sequence[int]:
        return sorted(os.sched_getaffinity(0))

    def clock_ticks_per_second(self) -> int:
        return os.sysconf("SC_CLK_TCK")

    def monotonic_ns(self) -> int:
        return time.monotonic_ns()

    def sleep_ns(self, duration_ns: int) -> None:
        duration = _bounded_int(
            duration_ns, 1, MAX_NOMINAL_WINDOW_NS, "sleep duration")
        deadline = time.monotonic_ns() + duration
        while True:
            remaining = deadline - time.monotonic_ns()
            if remaining <= 0:
                return
            time.sleep(remaining / 1_000_000_000)

    def read_bytes(self, path: str, maximum_bytes: int) -> bytes:
        _require(type(path) is str and path.startswith("/") and
                 "//" not in path and "/../" not in path and
                 not path.endswith("/.."), "host read path is not canonical")
        maximum = _bounded_int(
            maximum_bytes, 1, MAX_HOST_FILE_BYTES, "host read byte limit")
        descriptor = -1
        try:
            descriptor = os.open(
                path, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW)
            status = os.fstat(descriptor)
            _require(stat.S_ISREG(status.st_mode), "host read target is not regular")
            chunks: list[bytes] = []
            total = 0
            while True:
                chunk = os.read(descriptor, min(65536, maximum + 1 - total))
                if not chunk:
                    break
                chunks.append(chunk)
                total += len(chunk)
                _require(total <= maximum, "host read exceeded its byte limit")
            return b"".join(chunks)
        except AcquisitionError:
            raise
        except OSError as error:
            raise AcquisitionError(f"cannot read host source {path!r}: {error}") from error
        finally:
            if descriptor >= 0:
                os.close(descriptor)


def _topology_path(cpu: int, suffix: str) -> str:
    return f"{SYSFS_CPU_ROOT}/cpu{cpu}/{suffix}"


def _acquire_topology(
    reader: HostReader, allowed_cpus: Sequence[int],
    policy_primary_cpus: Sequence[int],
) -> dict[str, Any]:
    online = set(parse_cpu_list(
        reader.read_bytes(ONLINE_CPUS_PATH, MAX_HOST_FILE_BYTES),
        "online CPU list"))
    required = set(_canonical_cpu_list(allowed_cpus, label="allowed CPU set"))
    _require(required <= online, "allowed CPU set includes an offline CPU")
    policy_primaries = set(_canonical_cpu_list(
        policy_primary_cpus, label="policy primary CPU set", allow_empty=True))
    # Retain every online policy primary even when it is outside the launch
    # affinity closure.  The pure contract then applies its intended
    # online-and-allowed filter instead of treating the primary as missing.
    required.update(policy_primaries & online)
    raw_entries: dict[int, dict[str, Any]] = {}
    while required - set(raw_entries):
        cpu = min(required - set(raw_entries))
        siblings = parse_cpu_list(reader.read_bytes(
            _topology_path(cpu, "topology/thread_siblings_list"),
            MAX_HOST_FILE_BYTES), f"CPU {cpu} thread siblings")
        domain = parse_cpu_list(reader.read_bytes(
            _topology_path(cpu, "cache/index3/shared_cpu_list"),
            MAX_HOST_FILE_BYTES), f"CPU {cpu} domain CPUs")
        raw_entries[cpu] = {
            "cpu": cpu,
            "online": cpu in online,
            "physical_package_id": _parse_decimal(reader.read_bytes(
                _topology_path(cpu, "topology/physical_package_id"), 64),
                f"CPU {cpu} package ID"),
            "die_id": _parse_decimal(reader.read_bytes(
                _topology_path(cpu, "topology/die_id"), 64),
                f"CPU {cpu} die ID"),
            "core_id": _parse_decimal(reader.read_bytes(
                _topology_path(cpu, "topology/core_id"), 64),
                f"CPU {cpu} core ID"),
            "thread_siblings": siblings,
            "domain_cpus": domain,
        }
        required.update(siblings)
        required.update(domain)
        _require(len(required) <= contract.MAX_CPU_COUNT,
                 "topology closure is too large")
    _require(required <= online,
             "qualification topology contains an offline observed CPU")
    return contract.topology_record([raw_entries[cpu] for cpu in sorted(raw_entries)])


def _capture_snapshot(
    reader: HostReader, observed_cpus: Sequence[int],
) -> dict[str, Any]:
    started = reader.monotonic_ns()
    data = reader.read_bytes(PROC_STAT_PATH, MAX_PROC_STAT_BYTES)
    finished = reader.monotonic_ns()
    return contract.shared_snapshot_record(
        read_started_monotonic_ns=started,
        read_finished_monotonic_ns=finished,
        counters=parse_proc_stat(data, observed_cpus),
    )


def _acquisition_record(
    *, policy: dict[str, Any], requested_window_count: int,
    nominal_window_ns: int, frozen_pair: dict[str, int] | None,
    allowed_at_launch: list[int], allowed_after_scan: list[int],
    ticks_at_launch: int, ticks_after_scan: int, scan: dict[str, Any],
) -> dict[str, Any]:
    return {
        "schema": ACQUISITION_SCHEMA,
        "acquisition_method": ACQUISITION_METHOD,
        "sources": _sources_record(),
        "policy": copy.deepcopy(policy),
        "policy_sha256": contract.canonical_sha256(policy),
        "requested_window_count": requested_window_count,
        "nominal_window_ns": nominal_window_ns,
        "frozen_pair_from_prior_attempt": copy.deepcopy(frozen_pair),
        "allowed_cpu_set_at_launch": list(allowed_at_launch),
        "allowed_cpu_set_after_scan": list(allowed_after_scan),
        "clock_ticks_per_second_at_launch": ticks_at_launch,
        "clock_ticks_per_second_after_scan": ticks_after_scan,
        "topology_before_sha256": contract.canonical_sha256(scan["topology_before"]),
        "topology_after_sha256": contract.canonical_sha256(scan["topology_after"]),
        "scan": copy.deepcopy(scan),
        "scan_sha256": contract.canonical_sha256(scan),
        "host_mutation_performed": False,
        "candidate_timing_performed": False,
    }


def validate_acquisition_record(
    value: Any, *, expected_policy: Any, expected_frozen_pair: Any,
    expected_window_count: int, expected_nominal_window_ns: int,
) -> dict[str, Any]:
    _require(type(value) is dict and set(value) == ACQUISITION_KEYS,
             "pair qualification acquisition has an unexpected key set")
    _require(value["schema"] == ACQUISITION_SCHEMA,
             "pair qualification acquisition schema differs")
    _require(value["acquisition_method"] == ACQUISITION_METHOD,
             "pair qualification acquisition method differs")
    _require(type(value["sources"]) is dict and
             set(value["sources"]) == SOURCE_KEYS and
             contract.exact_json_equal(value["sources"], _sources_record()),
             "pair qualification acquisition sources differ")
    policy = contract.validate_qualification_policy(expected_policy)
    _require(contract.exact_json_equal(value["policy"], policy),
             "pair qualification acquisition policy differs")
    _require(value["policy_sha256"] == contract.canonical_sha256(policy),
             "pair qualification acquisition policy hash differs")
    windows = _bounded_int(
        expected_window_count, 1, contract.MAX_WINDOWS,
        "expected acquisition window count")
    nominal = _bounded_int(
        expected_nominal_window_ns, 1, MAX_NOMINAL_WINDOW_NS,
        "expected nominal window duration")
    _require(type(value["requested_window_count"]) is int and
             value["requested_window_count"] == windows,
             "pair qualification acquisition window count differs")
    _require(type(value["nominal_window_ns"]) is int and
             value["nominal_window_ns"] == nominal,
             "pair qualification acquisition nominal duration differs")
    frozen = _pair(expected_frozen_pair, "expected frozen pair")
    _require(contract.exact_json_equal(
        value["frozen_pair_from_prior_attempt"], frozen),
        "pair qualification acquisition frozen pair differs")
    allowed_launch = _canonical_cpu_list(
        value["allowed_cpu_set_at_launch"], label="launch allowed CPU set")
    allowed_after = _canonical_cpu_list(
        value["allowed_cpu_set_after_scan"], label="final allowed CPU set")
    _require(type(value["allowed_cpu_set_at_launch"]) is list and
             type(value["allowed_cpu_set_after_scan"]) is list and
             value["allowed_cpu_set_at_launch"] == allowed_launch and
             value["allowed_cpu_set_after_scan"] == allowed_after and
             allowed_launch == allowed_after,
             "allowed CPU set changed across qualification acquisition")
    ticks_launch = _bounded_int(
        value["clock_ticks_per_second_at_launch"], 1,
        contract.MAX_CLOCK_TICKS_PER_SECOND, "launch clock tick rate")
    ticks_after = _bounded_int(
        value["clock_ticks_per_second_after_scan"], 1,
        contract.MAX_CLOCK_TICKS_PER_SECOND, "final clock tick rate")
    _require(ticks_launch == ticks_after == policy["clock_ticks_per_second"],
             "clock tick rate changed or differs from policy")
    scan = contract.validate_pair_qualification_scan(
        value["scan"], expected_policy=policy, expected_frozen_pair=frozen)
    _require(scan["scan_window_count"] == windows,
             "pair qualification scan has the wrong window count")
    _require(scan["allowed_cpu_set_at_launch"] == allowed_launch,
             "pair qualification scan allowed CPU set differs")
    _require(all(window["window_ns"] >= nominal for window in scan["windows"]),
             "pair qualification scan contains a short observation window")
    _require(value["topology_before_sha256"] ==
             contract.canonical_sha256(scan["topology_before"]) and
             value["topology_after_sha256"] ==
             contract.canonical_sha256(scan["topology_after"]),
             "pair qualification acquisition topology hash differs")
    _require(value["scan_sha256"] == contract.canonical_sha256(scan),
             "pair qualification acquisition scan hash differs")
    _require(type(value["host_mutation_performed"]) is bool and
             value["host_mutation_performed"] is False,
             "pair qualification acquisition reports host mutation")
    _require(type(value["candidate_timing_performed"]) is bool and
             value["candidate_timing_performed"] is False and
             scan["candidate_timing_performed"] is False,
             "pair qualification acquisition reports candidate timing")
    expected = _acquisition_record(
        policy=policy, requested_window_count=windows,
        nominal_window_ns=nominal, frozen_pair=frozen,
        allowed_at_launch=allowed_launch, allowed_after_scan=allowed_after,
        ticks_at_launch=ticks_launch, ticks_after_scan=ticks_after, scan=scan)
    _require(contract.exact_json_equal(value, expected),
             "pair qualification acquisition is not canonical")
    return copy.deepcopy(expected)


def acquire_pair_qualification(
    reader: HostReader, *, policy_value: Any, window_count: int,
    nominal_window_ns: int, frozen_pair_from_prior_attempt: Any,
) -> dict[str, Any]:
    """Acquire and self-validate one read-only qualification scan."""
    policy = contract.validate_qualification_policy(policy_value)
    windows = _bounded_int(
        window_count, 1, contract.MAX_WINDOWS, "acquisition window count")
    nominal = _bounded_int(
        nominal_window_ns, 1, MAX_NOMINAL_WINDOW_NS,
        "acquisition nominal window duration")
    frozen = _pair(frozen_pair_from_prior_attempt, "frozen pair")
    allowed_launch = _canonical_cpu_list(
        reader.allowed_cpus(), label="launch allowed CPU set")
    ticks_launch = _bounded_int(
        reader.clock_ticks_per_second(), 1,
        contract.MAX_CLOCK_TICKS_PER_SECOND, "launch clock tick rate")
    _require(ticks_launch == policy["clock_ticks_per_second"],
             "host clock tick rate differs from qualification policy")
    topology_before = _acquire_topology(
        reader, allowed_launch, policy["candidate_primary_cpus"])
    candidates = contract.derive_candidate_pairs(
        policy, topology_before, allowed_launch)
    observed_cpus = sorted({
        cpu for candidate in candidates for cpu in candidate["domain_cpus"]
    })
    snapshots = [_capture_snapshot(reader, observed_cpus)]
    for unused_index in range(windows):
        reader.sleep_ns(nominal)
        snapshots.append(_capture_snapshot(reader, observed_cpus))
    topology_after = _acquire_topology(
        reader, allowed_launch, policy["candidate_primary_cpus"])
    allowed_after = _canonical_cpu_list(
        reader.allowed_cpus(), label="final allowed CPU set")
    ticks_after = _bounded_int(
        reader.clock_ticks_per_second(), 1,
        contract.MAX_CLOCK_TICKS_PER_SECOND, "final clock tick rate")
    _require(allowed_after == allowed_launch,
             "allowed CPU set changed across qualification acquisition")
    _require(ticks_after == ticks_launch,
             "clock tick rate changed across qualification acquisition")
    scan = contract.pair_qualification_scan_record(
        policy, allowed_cpu_set_at_launch=allowed_launch,
        topology_before=topology_before, topology_after=topology_after,
        snapshots=snapshots, frozen_pair_from_prior_attempt=frozen)
    value = _acquisition_record(
        policy=policy, requested_window_count=windows,
        nominal_window_ns=nominal, frozen_pair=frozen,
        allowed_at_launch=allowed_launch, allowed_after_scan=allowed_after,
        ticks_at_launch=ticks_launch, ticks_after_scan=ticks_after, scan=scan)
    return validate_acquisition_record(
        value, expected_policy=policy, expected_frozen_pair=frozen,
        expected_window_count=windows, expected_nominal_window_ns=nominal)


def load_acquisition_record(
    data: bytes, *, expected_policy: Any, expected_frozen_pair: Any,
    expected_window_count: int, expected_nominal_window_ns: int,
) -> dict[str, Any]:
    return validate_acquisition_record(
        contract.strict_json_loads(data, "pair qualification acquisition"),
        expected_policy=expected_policy,
        expected_frozen_pair=expected_frozen_pair,
        expected_window_count=expected_window_count,
        expected_nominal_window_ns=expected_nominal_window_ns)


def write_acquisition_record_exclusive(
    output_path: str | os.PathLike[str], value: Any, *,
    expected_policy: Any, expected_frozen_pair: Any,
    expected_window_count: int, expected_nominal_window_ns: int,
) -> dict[str, Any]:
    """Publish one canonical 0600 record without following or replacing names."""
    validated = validate_acquisition_record(
        value, expected_policy=expected_policy,
        expected_frozen_pair=expected_frozen_pair,
        expected_window_count=expected_window_count,
        expected_nominal_window_ns=expected_nominal_window_ns)
    content = contract.canonical_json_bytes(validated)
    path = Path(output_path)
    _require(path.is_absolute() and path.name not in ("", ".", ".."),
             "acquisition output path is not canonical and absolute")
    parent = path.parent
    try:
        _require(parent.resolve(strict=True) == parent,
                 "acquisition output parent is not canonical")
    except (OSError, RuntimeError) as error:
        raise AcquisitionError("cannot resolve acquisition output parent") from error
    parent_fd = -1
    file_fd = -1
    created = False
    published = False
    try:
        parent_fd = os.open(
            parent, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW)
        parent_status = os.fstat(parent_fd)
        _require(stat.S_ISDIR(parent_status.st_mode) and
                 parent_status.st_uid == os.geteuid() and
                 stat.S_IMODE(parent_status.st_mode) & 0o022 == 0,
                 "acquisition output parent is not owner-safe")
        file_fd = os.open(
            path.name,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC | os.O_NOFOLLOW,
            OUTPUT_FILE_MODE, dir_fd=parent_fd)
        created = True
        status = os.fstat(file_fd)
        _require(stat.S_ISREG(status.st_mode) and status.st_uid == os.geteuid() and
                 status.st_nlink == 1 and
                 stat.S_IMODE(status.st_mode) == OUTPUT_FILE_MODE,
                 "acquisition output file has unsafe identity or mode")
        offset = 0
        while offset < len(content):
            written = os.write(file_fd, content[offset:])
            _require(type(written) is int and written > 0,
                     "acquisition output write made no progress")
            offset += written
        os.fsync(file_fd)
        os.close(file_fd)
        file_fd = -1
        os.fsync(parent_fd)
        published = True
        return copy.deepcopy(validated)
    except AcquisitionError:
        raise
    except OSError as error:
        raise AcquisitionError(f"cannot publish acquisition record: {error}") from error
    finally:
        if file_fd >= 0:
            try:
                os.close(file_fd)
            except OSError:
                pass
        if created and not published:
            try:
                os.unlink(path.name, dir_fd=parent_fd)
                os.fsync(parent_fd)
            except OSError:
                pass
        if parent_fd >= 0:
            os.close(parent_fd)


__all__ = (
    "ACQUISITION_METHOD", "ACQUISITION_SCHEMA", "AcquisitionError",
    "HostReader", "SystemHostReader", "acquire_pair_qualification",
    "format_cpu_list", "load_acquisition_record", "parse_cpu_list",
    "parse_proc_stat", "validate_acquisition_record",
    "write_acquisition_record_exclusive",
)
