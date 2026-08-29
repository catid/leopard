#!/usr/bin/python3

"""Capture and verify a read-only endpoint census for passive v17 evidence.

This controller never changes affinity and never signals another process.  Its
snapshots disclose affinity eligibility at two non-atomic endpoints; they do
not claim interval completeness or CPU exclusivity.
"""

from __future__ import annotations

import argparse
import errno
import hashlib
import json
import os
import re
import stat
import sys
import time
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA = "leopard2-passive-same-uid-endpoint-census/v1"
POLICY_SCHEMA = "leopard2-passive-shared-host-policy/v1"
SEMANTICS = "non-atomic-endpoint-affinity-eligibility"
MUTATION_POLICY = "read-only; no affinity mutation, signal, or procfs write"
MAX_CPU_ID = 1_048_575
MAX_STATUS_BYTES = 1 << 20
MAX_STAT_BYTES = 1 << 20
MAX_TASKS = 1_000_000
CLOCK_TICKS_PER_SECOND = 100
MAX_NONIDLE_EXCESS_OVER_CHILD_WALL_CEILING_JIFFIES = 16
EXPECTED_INVOCATION_COUNT = 12
HEX256 = re.compile(r"^[0-9a-f]{64}$")
BOOT_ID = re.compile(
    r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$")
CPU_FIELDS = (
    "user", "nice", "system", "idle", "iowait", "irq", "softirq",
    "steal", "guest", "guest_nice",
)
NONIDLE_FIELDS = ("user", "nice", "system", "irq", "softirq", "steal")


class CensusError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise CensusError(message)


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False,
    ).encode("ascii")


def digest_payload(value: Mapping[str, Any]) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def bounded_read(path: Path, maximum: int) -> bytes:
    with path.open("rb", buffering=0) as source:
        data = source.read(maximum + 1)
    require(len(data) <= maximum, f"oversized procfs record: {path}")
    return data


def parse_cpu_list(text: str) -> list[int]:
    require(text != "" and text == text.strip(), "invalid empty/padded CPU list")
    values: list[int] = []
    previous = -1
    for token in text.split(","):
        require(token != "", "empty CPU-list token")
        if "-" in token:
            require(token.count("-") == 1, "invalid CPU range")
            first_text, last_text = token.split("-", 1)
        else:
            first_text = last_text = token
        require(first_text.isdecimal() and last_text.isdecimal(),
                "non-decimal CPU-list token")
        require((first_text == "0" or not first_text.startswith("0")) and
                (last_text == "0" or not last_text.startswith("0")),
                "noncanonical CPU-list integer")
        first, last = int(first_text), int(last_text)
        require(0 <= first <= last <= MAX_CPU_ID and first > previous,
                "overlapping, unordered, or out-of-range CPU list")
        values.extend(range(first, last + 1))
        previous = last
    require(values, "CPU list is empty")
    return values


def format_cpu_list(values: Sequence[int]) -> str:
    ordered = sorted(set(values))
    require(ordered and all(type(value) is int and 0 <= value <= MAX_CPU_ID
                            for value in ordered), "invalid CPU values")
    ranges: list[str] = []
    start = previous = ordered[0]
    for value in ordered[1:]:
        if value == previous + 1:
            previous = value
            continue
        ranges.append(str(start) if start == previous else f"{start}-{previous}")
        start = previous = value
    ranges.append(str(start) if start == previous else f"{start}-{previous}")
    return ",".join(ranges)


def parse_status(data: bytes) -> dict[str, Any]:
    text = data.decode("utf-8", "strict")
    require("\x00" not in text, "NUL in proc status")
    fields: dict[str, str] = {}
    for line in text.splitlines():
        if ":" not in line:
            continue
        name, value = line.split(":", 1)
        if name in ("Uid", "Cpus_allowed_list", "Name"):
            require(name not in fields, f"duplicate proc status field: {name}")
            fields[name] = value.strip()
    require(set(fields) == {"Uid", "Cpus_allowed_list", "Name"},
            "proc status lacks required fields")
    uid_tokens = fields["Uid"].split()
    require(len(uid_tokens) == 4 and all(token.isdecimal()
                                        for token in uid_tokens),
            "invalid proc status Uid")
    cpus = parse_cpu_list(fields["Cpus_allowed_list"])
    require(format_cpu_list(cpus) == fields["Cpus_allowed_list"],
            "noncanonical Cpus_allowed_list")
    name_bytes = fields["Name"].encode("utf-8")
    return {
        "uid": int(uid_tokens[0]),
        "cpus_allowed": cpus,
        "comm_size": len(name_bytes),
        "comm_sha256": hashlib.sha256(name_bytes).hexdigest(),
    }


def parse_starttime(data: bytes) -> int:
    require(b"\x00" not in data and data.endswith(b"\n"),
            "invalid proc stat framing")
    close = data.rfind(b")")
    require(close > 1 and data[0:data.find(b" ")].isdigit(),
            "invalid proc stat identity")
    tail = data[close + 2:].split()
    # tail[0] is field 3 (state), so field 22 (starttime) is tail[19].
    require(len(tail) >= 20 and tail[19].isdigit(),
            "proc stat lacks starttime")
    return int(tail[19])


def identity(path: Path, stat_path: Path) -> dict[str, int]:
    node = path.stat(follow_symlinks=False)
    require(stat.S_ISDIR(node.st_mode), f"proc task is not a directory: {path}")
    return {
        "device": node.st_dev,
        "inode": node.st_ino,
        "starttime_ticks": parse_starttime(bounded_read(stat_path, MAX_STAT_BYTES)),
    }


def namespace_identity(path: Path) -> dict[str, int]:
    link = path.stat(follow_symlinks=False)
    require(stat.S_ISLNK(link.st_mode),
            f"namespace path is not a symlink: {path}")
    # The procfs symlink belongs to the short-lived collector process.  Follow
    # it so endpoint snapshots bind the namespace object shared by both
    # collectors, rather than two necessarily different procfs dentries.
    node = path.stat(follow_symlinks=True)
    require(node.st_dev >= 0 and node.st_ino > 0,
            f"namespace target identity is invalid: {path}")
    return {"device": node.st_dev, "inode": node.st_ino}


def read_cpu_stat(cpus: Sequence[int]) -> dict[str, Any]:
    wanted = set(cpus)
    observed: dict[int, dict[str, Any]] = {}
    data = bounded_read(Path("/proc/stat"), 8 << 20)
    for line in data.splitlines():
        parts = line.split()
        if not parts or not parts[0].startswith(b"cpu") or parts[0] == b"cpu":
            continue
        suffix = parts[0][3:]
        if not suffix.isdigit() or int(suffix) not in wanted:
            continue
        cpu = int(suffix)
        require(cpu not in observed and 8 <= len(parts) - 1 <= len(CPU_FIELDS),
                "invalid or duplicate /proc/stat CPU record")
        values = []
        for token in parts[1:]:
            require(token.isdigit(), "non-integer /proc/stat counter")
            values.append(int(token))
        values.extend([0] * (len(CPU_FIELDS) - len(values)))
        fields = dict(zip(CPU_FIELDS, values))
        observed[cpu] = {
            "cpu": cpu,
            "fields": fields,
            "total_jiffies": sum(values[:8]),
            "nonidle_jiffies": sum(fields[name] for name in NONIDLE_FIELDS),
        }
    require(set(observed) == wanted, "reserved CPU missing from /proc/stat")
    return {str(cpu): observed[cpu] for cpu in sorted(observed)}


def vanished_record(pid: int, tid: int | None, stage: str) -> dict[str, Any]:
    return {"pid": pid, "tid": tid, "stage": stage}


def thread_entry(pid: int, tid: int, uid: int,
                 reserved: Sequence[int]) -> dict[str, Any]:
    process = Path(f"/proc/{pid}")
    task = process / "task" / str(tid)
    process_before = identity(process, process / "stat")
    task_before = identity(task, task / "stat")
    parsed = parse_status(bounded_read(task / "status", MAX_STATUS_BYTES))
    require(parsed["uid"] == uid, "same-UID task changed ownership")
    affinity = sorted(os.sched_getaffinity(tid))
    require(affinity == parsed["cpus_allowed"],
            "procfs and sched_getaffinity disagree")
    task_after = identity(task, task / "stat")
    process_after = identity(process, process / "stat")
    require(process_before == process_after and task_before == task_after,
            "proc task identity changed during scan")
    intersection = sorted(set(affinity).intersection(reserved))
    return {
        "pid": pid,
        "tid": tid,
        "uid": uid,
        "process_identity": process_before,
        "thread_identity": task_before,
        "comm_size": parsed["comm_size"],
        "comm_sha256": parsed["comm_sha256"],
        "cpus_allowed": affinity,
        "reserved_pair_intersection": intersection,
    }


def scan_same_uid(uid: int, reserved: Sequence[int]) -> dict[str, Any]:
    entries: list[dict[str, Any]] = []
    vanished: list[dict[str, Any]] = []
    numeric = sorted(
        (int(path.name) for path in Path("/proc").iterdir()
         if path.name.isdecimal()),
    )
    require(len(numeric) <= MAX_TASKS, "proc process census is unbounded")
    for pid in numeric:
        process = Path(f"/proc/{pid}")
        try:
            parsed = parse_status(bounded_read(process / "status", MAX_STATUS_BYTES))
        except FileNotFoundError:
            vanished.append(vanished_record(pid, None, "process-status"))
            continue
        if parsed["uid"] != uid:
            continue
        try:
            tids = sorted(int(path.name) for path in (process / "task").iterdir()
                          if path.name.isdecimal())
        except FileNotFoundError:
            vanished.append(vanished_record(pid, None, "task-list"))
            continue
        require(len(entries) + len(tids) <= MAX_TASKS,
                "proc thread census is unbounded")
        for tid in tids:
            try:
                entries.append(thread_entry(pid, tid, uid, reserved))
            except OSError as error:
                if error.errno in (errno.ENOENT, errno.ESRCH):
                    vanished.append(vanished_record(pid, tid, "thread-read"))
                    continue
                raise
    entries.sort(key=lambda item: (item["pid"], item["tid"]))
    require(len({(item["pid"], item["tid"]) for item in entries}) == len(entries),
            "duplicate census thread identity")
    by_process: dict[int, list[dict[str, Any]]] = {}
    for entry in entries:
        by_process.setdefault(entry["pid"], []).append(entry)
    pair_processes = sum(any(item["reserved_pair_intersection"] for item in group)
                         for group in by_process.values())
    nonuniform = sum(len({tuple(item["cpus_allowed"]) for item in group}) > 1
                     for group in by_process.values())
    confined_threads = sum(
        set(item["cpus_allowed"]).issubset(set(reserved)) for item in entries)
    confined_processes = sum(
        any(set(item["cpus_allowed"]).issubset(set(reserved)) for item in group)
        for group in by_process.values())
    return {
        "entries": entries,
        "vanished_during_scan": vanished,
        "summary": {
            "retained_process_count": len(by_process),
            "retained_thread_count": len(entries),
            "pair_eligible_process_count": pair_processes,
            "pair_eligible_thread_count": sum(
                bool(item["reserved_pair_intersection"]) for item in entries),
            "nonuniform_process_count": nonuniform,
            "confined_to_reserved_subset_process_count": confined_processes,
            "confined_to_reserved_subset_thread_count": confined_threads,
            "vanished_record_count": len(vanished),
        },
    }


def capture(phase: str, reserved: Sequence[int]) -> dict[str, Any]:
    require(phase in ("pre", "post"), "invalid census phase")
    require(len(reserved) == 2 and len(set(reserved)) == 2 and
            all(type(cpu) is int and 0 <= cpu <= MAX_CPU_ID for cpu in reserved),
            "invalid reserved CPU pair")
    reserved = sorted(reserved)
    uid = os.getuid()
    collector_affinity = sorted(os.sched_getaffinity(0))
    require(not set(collector_affinity).intersection(reserved),
            "census collector affinity includes a reserved CPU")
    clock_ticks = os.sysconf("SC_CLK_TCK")
    require(type(clock_ticks) is int and
            clock_ticks == CLOCK_TICKS_PER_SECOND,
            "passive contract requires SC_CLK_TCK=100")
    if phase == "post":
        boundary_ns = time.clock_gettime_ns(time.CLOCK_MONOTONIC)
        cpu_stat = read_cpu_stat(reserved)
    scan_started = time.clock_gettime_ns(time.CLOCK_MONOTONIC)
    census = scan_same_uid(uid, reserved)
    scan_finished = time.clock_gettime_ns(time.CLOCK_MONOTONIC)
    if phase == "pre":
        cpu_stat = read_cpu_stat(reserved)
        boundary_ns = time.clock_gettime_ns(time.CLOCK_MONOTONIC)
    require(scan_finished >= scan_started, "census monotonic interval regressed")
    summary = census["summary"]
    require(summary["retained_process_count"] > 0 and
            summary["retained_thread_count"] > 0,
            "same-UID census retained no tasks")
    require(summary["confined_to_reserved_subset_thread_count"] == 0,
            "same-UID thread is confined to the reserved pair")
    payload = {
        "schema": SCHEMA,
        "phase": phase,
        "semantics": SEMANTICS,
        "mutation_policy": MUTATION_POLICY,
        "uid": uid,
        "boot_id": bounded_read(Path("/proc/sys/kernel/random/boot_id"), 128)
            .decode("ascii", "strict").strip(),
        "namespaces": {
            "pid": namespace_identity(Path("/proc/self/ns/pid")),
            "user": namespace_identity(Path("/proc/self/ns/user")),
        },
        "clock_ticks_per_second": clock_ticks,
        "reserved_cpus": reserved,
        "scan_started_monotonic_ns": scan_started,
        "scan_finished_monotonic_ns": scan_finished,
        "activity_boundary_monotonic_ns": boundary_ns,
        "collector": {
            "pid": os.getpid(),
            "allowed_cpus": collector_affinity,
            "reserved_pair_excluded": True,
        },
        "proc_stat": cpu_stat,
        "same_uid_processes": census,
        "capabilities": {
            "atomic_snapshot": False,
            "interval_complete": False,
            "records_execution_history": False,
            "establishes_cpu_exclusivity": False,
        },
    }
    return {**payload, "digest": digest_payload(payload)}


def exact_keys(value: Any, keys: set[str], label: str) -> Mapping[str, Any]:
    require(type(value) is dict and set(value) == keys, f"{label} keys differ")
    return value


def validate_snapshot(value: Any, phase: str) -> Mapping[str, Any]:
    value = exact_keys(value, {
        "schema", "phase", "semantics", "mutation_policy", "uid", "boot_id",
        "namespaces", "clock_ticks_per_second", "reserved_cpus",
        "scan_started_monotonic_ns", "scan_finished_monotonic_ns",
        "activity_boundary_monotonic_ns", "collector", "proc_stat",
        "same_uid_processes", "capabilities", "digest",
    }, f"{phase} census")
    payload = dict(value)
    digest = payload.pop("digest")
    require(type(digest) is str and HEX256.fullmatch(digest) is not None and
            digest == digest_payload(payload), f"{phase} census digest differs")
    capabilities = exact_keys(value["capabilities"], {
        "atomic_snapshot", "interval_complete", "records_execution_history",
        "establishes_cpu_exclusivity",
    }, f"{phase} capabilities")
    require(value["schema"] == SCHEMA and value["phase"] == phase and
            value["semantics"] == SEMANTICS and
            value["mutation_policy"] == MUTATION_POLICY and
            type(value["clock_ticks_per_second"]) is int and
            value["clock_ticks_per_second"] == CLOCK_TICKS_PER_SECOND and
            all(capabilities[name] is False for name in capabilities),
            f"{phase} census contract differs")
    reserved_cpus = value["reserved_cpus"]
    require(type(value["uid"]) is int and value["uid"] >= 0 and
            type(value["boot_id"]) is str and BOOT_ID.fullmatch(value["boot_id"])
            is not None and type(reserved_cpus) is list and
            all(type(cpu) is int for cpu in reserved_cpus) and
            reserved_cpus == [52, 116],
            f"{phase} census host identity differs")
    namespaces = exact_keys(value["namespaces"], {"pid", "user"},
                            f"{phase} namespaces")
    for name in ("pid", "user"):
        record = exact_keys(namespaces[name], {"device", "inode"},
                            f"{phase} {name} namespace")
        require(all(type(record[field]) is int and record[field] > 0
                    for field in ("device", "inode")),
                f"{phase} namespace identity is invalid")
    started = value["scan_started_monotonic_ns"]
    finished = value["scan_finished_monotonic_ns"]
    boundary = value["activity_boundary_monotonic_ns"]
    require(all(type(item) is int and item >= 0
                for item in (started, finished, boundary)) and
            started <= finished and
            ((phase == "pre" and finished <= boundary) or
             (phase == "post" and boundary <= started)),
            f"{phase} census interval is invalid")
    collector = exact_keys(value["collector"], {
        "pid", "allowed_cpus", "reserved_pair_excluded"},
        f"{phase} collector")
    collector_cpus = collector["allowed_cpus"]
    require(type(collector["pid"]) is int and collector["pid"] > 0 and
            collector["reserved_pair_excluded"] is True and
            type(collector_cpus) is list and collector_cpus ==
                sorted(set(collector_cpus)) and
            all(type(cpu) is int and 0 <= cpu <= MAX_CPU_ID
                for cpu in collector_cpus) and collector_cpus and
            not set(collector_cpus).intersection(value["reserved_cpus"]),
            f"{phase} collector used a reserved CPU")
    proc_stat = exact_keys(value["proc_stat"], {"52", "116"},
                           f"{phase} proc stat")
    for cpu in (52, 116):
        record = exact_keys(proc_stat[str(cpu)], {
            "cpu", "fields", "total_jiffies", "nonidle_jiffies"},
            f"{phase} CPU {cpu}")
        fields = exact_keys(record["fields"], set(CPU_FIELDS),
                            f"{phase} CPU {cpu} fields")
        require(type(record["cpu"]) is int and record["cpu"] == cpu and
                all(type(fields[name]) is int and fields[name] >= 0
                    for name in CPU_FIELDS) and
                type(record["total_jiffies"]) is int and
                type(record["nonidle_jiffies"]) is int and
                record["total_jiffies"] == sum(fields[name]
                                                for name in CPU_FIELDS[:8]) and
                record["nonidle_jiffies"] == sum(fields[name]
                                                  for name in NONIDLE_FIELDS),
                f"{phase} CPU counters differ")
    census = exact_keys(value["same_uid_processes"], {
        "entries", "vanished_during_scan", "summary"},
        f"{phase} same-UID census")
    entries = census["entries"]
    require(type(entries) is list and 0 < len(entries) <= MAX_TASKS,
            f"{phase} retained thread census is invalid")
    identities: list[tuple[int, int]] = []
    by_process: dict[int, list[Mapping[str, Any]]] = {}
    for index, entry_value in enumerate(entries):
        entry = exact_keys(entry_value, {
            "pid", "tid", "uid", "process_identity", "thread_identity",
            "comm_size", "comm_sha256", "cpus_allowed",
            "reserved_pair_intersection"}, f"{phase} entry {index}")
        require(all(type(entry[name]) is int and entry[name] > 0
                    for name in ("pid", "tid")) and
                type(entry["uid"]) is int and
                entry["uid"] == value["uid"] and
                type(entry["comm_size"]) is int and entry["comm_size"] >= 0 and
                type(entry["comm_sha256"]) is str and
                HEX256.fullmatch(entry["comm_sha256"]) is not None,
                f"{phase} entry identity differs")
        for identity_name in ("process_identity", "thread_identity"):
            identity_value = exact_keys(entry[identity_name], {
                "device", "inode", "starttime_ticks"},
                f"{phase} entry {identity_name}")
            require(all(type(identity_value[name]) is int and
                        identity_value[name] >= 0 for name in identity_value) and
                    identity_value["inode"] > 0,
                    f"{phase} proc identity differs")
        cpus = entry["cpus_allowed"]
        intersection = entry["reserved_pair_intersection"]
        require(type(cpus) is list and cpus == sorted(set(cpus)) and cpus and
                all(type(cpu) is int and 0 <= cpu <= MAX_CPU_ID for cpu in cpus)
                and type(intersection) is list and
                all(type(cpu) is int for cpu in intersection) and
                intersection == sorted(set(cpus).intersection((52, 116))),
                f"{phase} entry affinity differs")
        identities.append((entry["pid"], entry["tid"]))
        by_process.setdefault(entry["pid"], []).append(entry)
    require(identities == sorted(identities) and
            len(set(identities)) == len(identities),
            f"{phase} census entries are not sorted and unique")
    for group in by_process.values():
        require(len({canonical_bytes(item["process_identity"])
                     for item in group}) == 1,
                f"{phase} process identity differs across threads")
    vanished = census["vanished_during_scan"]
    require(type(vanished) is list and len(vanished) <= MAX_TASKS,
            f"{phase} vanished census is invalid")
    for index, item_value in enumerate(vanished):
        item = exact_keys(item_value, {"pid", "tid", "stage"},
                          f"{phase} vanished {index}")
        require(type(item["pid"]) is int and item["pid"] > 0 and
                (item["tid"] is None or
                 (type(item["tid"]) is int and item["tid"] > 0)) and
                item["stage"] in ("process-status", "task-list", "thread-read"),
                f"{phase} vanished record differs")
    summary = exact_keys(census["summary"], {
        "retained_process_count", "retained_thread_count",
        "pair_eligible_process_count", "pair_eligible_thread_count",
        "nonuniform_process_count",
        "confined_to_reserved_subset_process_count",
        "confined_to_reserved_subset_thread_count", "vanished_record_count"},
        f"{phase} summary")
    expected_summary = {
        "retained_process_count": len(by_process),
        "retained_thread_count": len(entries),
        "pair_eligible_process_count": sum(
            any(item["reserved_pair_intersection"] for item in group)
            for group in by_process.values()),
        "pair_eligible_thread_count": sum(
            bool(item["reserved_pair_intersection"]) for item in entries),
        "nonuniform_process_count": sum(
            len({tuple(item["cpus_allowed"]) for item in group}) > 1
            for group in by_process.values()),
        "confined_to_reserved_subset_process_count": sum(
            any(set(item["cpus_allowed"]).issubset((52, 116))
                for item in group) for group in by_process.values()),
        "confined_to_reserved_subset_thread_count": sum(
            set(item["cpus_allowed"]).issubset((52, 116))
            for item in entries),
        "vanished_record_count": len(vanished),
    }
    require(all(type(item) is int and item >= 0 for item in summary.values()) and
            summary == expected_summary and
            summary["confined_to_reserved_subset_thread_count"] == 0,
            f"{phase} census summary differs")
    return value


def cpu_delta(before: Mapping[str, Any], after: Mapping[str, Any]) -> dict[str, Any]:
    require(before["cpu"] == after["cpu"], "CPU identity changed")
    fields = {}
    for name in CPU_FIELDS:
        first = before["fields"][name]
        last = after["fields"][name]
        require(type(first) is int and type(last) is int and 0 <= first <= last,
                "CPU counter regressed")
        fields[name] = last - first
    return {
        "cpu": before["cpu"],
        "fields": fields,
        "total_jiffies": sum(fields[name] for name in CPU_FIELDS[:8]),
        "nonidle_jiffies": sum(fields[name] for name in NONIDLE_FIELDS),
    }


def validate_controller(value: Any) -> Mapping[str, Any]:
    value = exact_keys(value, {
        "schema", "acquisition_generation", "wrapper_pid",
        "before_allowed_cpus", "after_allowed_cpus",
        "runner_launch_allowed_cpus", "benchmark_cpu", "reserved_sibling",
        "affinity_mutation_scope", "active_affinity_supervisor_executed",
    }, "passive controller affinity")
    before = value["before_allowed_cpus"]
    after = value["after_allowed_cpus"]
    launch = value["runner_launch_allowed_cpus"]
    require(value["schema"] == "leopard2-v17-passive-controller-affinity/v1" and
            value["acquisition_generation"] == "passive-v1" and
            type(value["wrapper_pid"]) is int and value["wrapper_pid"] > 0 and
            type(before) is list and before == sorted(set(before)) and
            type(after) is list and after == sorted(set(after)) and after and
            type(launch) is list and launch == sorted(set(launch)) and
            all(type(cpu) is int and 0 <= cpu <= MAX_CPU_ID
                for cpu in before + after + launch) and
            52 in before and 116 in before and
            after == [cpu for cpu in before if cpu not in (52, 116)] and
            launch == before and
            type(value["benchmark_cpu"]) is int and
            value["benchmark_cpu"] == 52 and
            type(value["reserved_sibling"]) is int and
            value["reserved_sibling"] == 116 and
            value["affinity_mutation_scope"] ==
                "wrapper-process-and-owned-descendants-only" and
            value["active_affinity_supervisor_executed"] is False,
            "passive controller affinity contract differs")
    return value


def raw_cpu_stat(value: Any, cpu: int, label: str) -> Mapping[str, Any]:
    value = exact_keys(value, {"cpu", "fields", "total_jiffies"}, label)
    fields = exact_keys(value["fields"], set(CPU_FIELDS[:8]),
                        f"{label} fields")
    require(type(value["cpu"]) is int and value["cpu"] == cpu and
            all(type(fields[name]) is int and fields[name] >= 0
                for name in CPU_FIELDS[:8]) and
            type(value["total_jiffies"]) is int and
            value["total_jiffies"] == sum(fields.values()),
            f"{label} counters differ")
    return value


def endpoint_identity_key(entry: Mapping[str, Any]) -> tuple[int, ...]:
    process = entry["process_identity"]
    thread = entry["thread_identity"]
    return (
        entry["pid"], entry["tid"],
        process["device"], process["inode"], process["starttime_ticks"],
        thread["device"], thread["inode"], thread["starttime_ticks"],
    )


def child_wall_reference(raw: Mapping[str, Any]) -> dict[str, int]:
    """Return the descriptive wall-time reference for the CPU52 screen.

    This arithmetic does not attribute /proc/stat ticks to a process.  It only
    defines a conservative rejection reference that is independently
    recomputed from the twelve retained child wall durations.
    """
    invocations = raw.get("invocations")
    require(type(invocations) is list and
            len(invocations) == EXPECTED_INVOCATION_COUNT,
            "passive campaign invocation census differs")
    durations: list[int] = []
    for index, invocation in enumerate(invocations):
        require(type(invocation) is dict,
                f"passive invocation {index} is not an object")
        duration = invocation.get("duration_ns")
        require(type(duration) is int and duration > 0,
                f"passive invocation {index} duration differs")
        durations.append(duration)
    total = sum(durations)
    numerator = total * CLOCK_TICKS_PER_SECOND
    ceiling = (numerator + 1_000_000_000 - 1) // 1_000_000_000
    return {
        "child_wall_duration_total_ns": total,
        "child_wall_time_ceiling_jiffies": ceiling,
    }


def load_json(path: Path) -> Any:
    data = bounded_read(path, 256 << 20)
    return json.loads(data.decode("utf-8", "strict"),
                      object_pairs_hook=lambda pairs: _unique_object(pairs))


def _unique_object(pairs: Sequence[tuple[str, Any]]) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, item in pairs:
        require(key not in value, f"duplicate JSON key: {key}")
        value[key] = item
    return value


def compare(pre: Mapping[str, Any], post: Mapping[str, Any],
            raw: Mapping[str, Any], controller: Mapping[str, Any]) \
        -> dict[str, Any]:
    pre = validate_snapshot(pre, "pre")
    post = validate_snapshot(post, "post")
    controller = validate_controller(controller)
    require(pre["uid"] == post["uid"] and
            pre["boot_id"] == post["boot_id"] and
            pre["namespaces"] == post["namespaces"] and
            pre["reserved_cpus"] == post["reserved_cpus"],
            "passive census identity changed")
    require(type(raw) is dict and "supervision" in raw and
            raw["supervision"] is None,
            "passive campaign supervision is not null")
    campaign = raw.get("campaign")
    raw_launch_cpus = (campaign.get("allowed_cpu_set_at_launch")
                       if type(campaign) is dict else None)
    require(type(campaign) is dict and type(raw_launch_cpus) is list and
            raw_launch_cpus == sorted(set(raw_launch_cpus)) and
            all(type(cpu) is int and 0 <= cpu <= MAX_CPU_ID
                for cpu in raw_launch_cpus) and
            raw_launch_cpus == controller["runner_launch_allowed_cpus"],
            "passive runner launch affinity differs from the controller")
    isolation = raw.get("isolation")
    require(type(isolation) is dict and
            set(isolation) == {
                "accepted", "after", "before", "benchmark_cpu", "delta",
                "pair_lease", "policy", "reserved_sibling", "schema"} and
            isolation.get("schema") == "leopard2-main-compare-isolation/v1" and
            type(isolation.get("benchmark_cpu")) is int and
            isolation.get("benchmark_cpu") == 52 and
            type(isolation.get("reserved_sibling")) is int and
            isolation.get("reserved_sibling") == 116 and
            isolation.get("accepted") is True,
            "passive campaign isolation was not accepted")
    isolation_before = exact_keys(isolation["before"], {
        "benchmark_cpu", "monotonic_ns", "reserved_sibling"},
        "raw isolation before")
    isolation_after = exact_keys(isolation["after"], {
        "benchmark_cpu", "monotonic_ns", "reserved_sibling"},
        "raw isolation after")
    before_ns = isolation_before["monotonic_ns"]
    after_ns = isolation_after["monotonic_ns"]
    require(type(before_ns) is int and type(after_ns) is int and
            pre["scan_finished_monotonic_ns"] <=
                pre["activity_boundary_monotonic_ns"] <= before_ns <= after_ns <=
                post["activity_boundary_monotonic_ns"] <=
                post["scan_started_monotonic_ns"],
            "passive census does not enclose campaign isolation")
    for cpu, raw_name in ((52, "benchmark_cpu"), (116, "reserved_sibling")):
        raw_before = raw_cpu_stat(
            isolation_before[raw_name], cpu, f"raw CPU {cpu} before")
        raw_after = raw_cpu_stat(
            isolation_after[raw_name], cpu, f"raw CPU {cpu} after")
        outer_before = pre["proc_stat"][str(cpu)]["fields"]
        outer_after = post["proc_stat"][str(cpu)]["fields"]
        require(all(
            outer_before[name] <= raw_before["fields"][name] <=
                raw_after["fields"][name] <= outer_after[name]
            for name in CPU_FIELDS[:8]),
            f"raw CPU {cpu} counters escape the passive outer interval")
    cpus = pre["reserved_cpus"]
    outer = {
        str(cpu): cpu_delta(pre["proc_stat"][str(cpu)],
                            post["proc_stat"][str(cpu)]) for cpu in cpus
    }
    benchmark = outer[str(cpus[0])]
    sibling = outer[str(cpus[1])]
    require(sibling["nonidle_jiffies"] == 0,
            "outer passive sibling activity gate failed")
    wall_reference = child_wall_reference(raw)
    require(after_ns > before_ns and
            wall_reference["child_wall_duration_total_ns"] <=
                after_ns - before_ns,
            "passive child wall durations escape the isolation interval")
    observed = benchmark["nonidle_jiffies"]
    excess = max(
        0, observed - wall_reference["child_wall_time_ceiling_jiffies"])
    require(
        excess <= MAX_NONIDLE_EXCESS_OVER_CHILD_WALL_CEILING_JIFFIES,
        "outer benchmark CPU nonidle excess rejection threshold exceeded")
    outer_contamination = {
        "clock_ticks_per_second": CLOCK_TICKS_PER_SECOND,
        **wall_reference,
        "benchmark_cpu_nonidle_jiffies": observed,
        "benchmark_cpu_auxiliary_class_jiffies": {
            name: benchmark["fields"][name]
            for name in ("nice", "iowait", "irq", "softirq", "steal")
        },
        "benchmark_cpu_nonidle_excess_over_child_wall_ceiling_jiffies":
            excess,
        "policy_max_nonidle_excess_over_child_wall_ceiling_jiffies":
            MAX_NONIDLE_EXCESS_OVER_CHILD_WALL_CEILING_JIFFIES,
        "reserved_sibling_nonidle_jiffies": sibling["nonidle_jiffies"],
        "auxiliary_class_policy": (
            "nice, irq, softirq, and steal are included in aggregate nonidle; "
            "iowait is disclosed as diagnostic idle accounting; no auxiliary "
            "class has a separate zero-tolerance gate"
        ),
        "interpretation": (
            "descriptive one-sided rejection screen only; the excess is not "
            "process ownership attribution or an upper bound on interference"
        ),
    }
    pre_entries = {
        endpoint_identity_key(item): item
        for item in pre["same_uid_processes"]["entries"]
    }
    post_entries = {
        endpoint_identity_key(item): item
        for item in post["same_uid_processes"]["entries"]
    }
    persistent = sorted(set(pre_entries).intersection(post_entries))
    require(all(pre_entries[key]["cpus_allowed"] ==
                post_entries[key]["cpus_allowed"] for key in persistent),
            "persistent same-UID task changed affinity between endpoints")
    wrapper_pid = controller["wrapper_pid"]
    wrapper_keys = [
        key for key in persistent if key[0] == wrapper_pid and key[1] == wrapper_pid
    ]
    require(len(wrapper_keys) == 1 and
            pre_entries[wrapper_keys[0]]["cpus_allowed"] ==
                controller["after_allowed_cpus"] and
            post_entries[wrapper_keys[0]]["cpus_allowed"] ==
                controller["after_allowed_cpus"],
            "wrapper endpoint identity/affinity is not persistent")
    wrapper = pre_entries[wrapper_keys[0]]
    return {
        "schema": POLICY_SCHEMA,
        "status": "complete",
        "policy_evaluation_complete": True,
        "cpu_pair_exclusive": False,
        "interval_complete_task_observation": False,
        "foreign_cpu52_work_attributable": False,
        "causal_performance_claim_eligible": False,
        "promotion_eligible": False,
        "promotion_passed": False,
        "census_collector_foreign_process_affinity_mutation_performed": False,
        "census_collector_foreign_process_signal_sent": False,
        "campaign_supervision": None,
        "outer_cpu_activity": outer,
        "outer_contamination": outer_contamination,
        "persistent_same_uid_thread_count": len(persistent),
        "wrapper_endpoint": {
            "pid": wrapper_pid,
            "process_identity": wrapper["process_identity"],
            "thread_identity": wrapper["thread_identity"],
            "cpus_allowed": wrapper["cpus_allowed"],
        },
        "interpretation": (
            "endpoint affinity eligibility and jiffy counters do not establish "
            "benchmark-CPU exclusivity or freedom from transient interference"
        ),
    }


def write_exclusive(path: Path, value: Mapping[str, Any]) -> None:
    require(path.is_absolute() and path.parent.resolve(strict=True) == path.parent and
            not path.exists(), "output path is not a new canonical absolute file")
    data = canonical_bytes(value) + b"\n"
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL |
                         getattr(os, "O_CLOEXEC", 0), 0o600)
    try:
        view = memoryview(data)
        while view:
            count = os.write(descriptor, view)
            require(count > 0, "short census output write")
            view = view[count:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)


def self_test() -> None:
    require(parse_cpu_list("0-1,4,52,116") == [0, 1, 4, 52, 116],
            "CPU-list parser failed")
    require(format_cpu_list([0, 1, 4, 52, 116]) == "0-1,4,52,116",
            "CPU-list formatter failed")
    for bad in ("", "01", "1-0", "1,1", "1,,2", "1-2,2-3", " 1"):
        try:
            parse_cpu_list(bad)
        except CensusError:
            pass
        else:
            raise CensusError(f"self-test accepted CPU list: {bad!r}")
    base = {
        "cpu": 52,
        "fields": {name: 10 for name in CPU_FIELDS},
        "total_jiffies": 80,
        "nonidle_jiffies": 60,
    }
    after = json.loads(json.dumps(base))
    after["fields"]["user"] += 1
    delta = cpu_delta(base, after)
    require(delta["fields"]["user"] == 1 and delta["nonidle_jiffies"] == 1,
            "CPU delta self-test failed")
    print("passive environment census self-test passed")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    capture_parser = subparsers.add_parser("capture")
    capture_parser.add_argument("--phase", choices=("pre", "post"), required=True)
    capture_parser.add_argument("--reserved-cpus", default="52,116")
    capture_parser.add_argument("--output", type=Path, required=True)
    verify_parser = subparsers.add_parser("verify")
    verify_parser.add_argument("--input", type=Path, required=True)
    verify_parser.add_argument("--phase", choices=("pre", "post"), required=True)
    compare_parser = subparsers.add_parser("compare")
    compare_parser.add_argument("--pre", type=Path, required=True)
    compare_parser.add_argument("--post", type=Path, required=True)
    compare_parser.add_argument("--raw", type=Path, required=True)
    compare_parser.add_argument("--controller", type=Path, required=True)
    compare_parser.add_argument("--output", type=Path, required=True)
    subparsers.add_parser("self-test")
    options = parser.parse_args()
    if options.command == "self-test":
        self_test()
        return 0
    if options.command == "capture":
        reserved = parse_cpu_list(options.reserved_cpus)
        value = capture(options.phase, reserved)
        write_exclusive(options.output, value)
        return 0
    if options.command == "verify":
        validate_snapshot(load_json(options.input), options.phase)
        return 0
    result = compare(
        load_json(options.pre), load_json(options.post), load_json(options.raw),
        load_json(options.controller))
    write_exclusive(options.output, result)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (CensusError, UnicodeError, json.JSONDecodeError, OSError) as error:
        print(f"passive environment census failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error
