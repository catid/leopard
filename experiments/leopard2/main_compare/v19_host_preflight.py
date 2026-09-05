#!/usr/bin/python3
"""Read-only, pre-lane host gate for leopard-79h.38.5.4.8.2.2.1.

This module cannot build, qualify, time, create a lane, or arm v19. Two matching
observations bound one capture; they are not an atomic snapshot, a lease, or
continuous resource authority. The future orchestrator must revalidate before
each consequential stage and keep its children in the observed cgroup.
Replay checks internal consistency, not authenticity of a claimed host. The
file-identity digest is a capture fingerprint, not a standalone attestation.

Linux memory-controller semantics:
https://docs.kernel.org/admin-guide/cgroup-v2.html#memory-interface-files
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
from pathlib import PurePosixPath
import re
import resource
import stat
import sys
import time
from typing import Any


PREREGISTRATION_SHA256 = (
    "27c1b7d76a0ecdbe194d6e6b62c01e48b1c7d10fc8ef99ebad4d76238669f0c1")
SCHEMA = "leopard2-v19-pre-lane-host-capture/v1"
MEMORY_MAX = 536_870_912
CPUS = list(range(128))
CPU_MODEL = "AMD Ryzen Threadripper PRO 9985WX 64-Cores"
CGROUP_MOUNT = "/sys/fs/cgroup"
SYS_CPU = "/sys/devices/system/cpu"
SYS_NODE = "/sys/devices/system/node"
MAX_TEXT = 1 << 20
MAX_RECORD_BYTES = 2 << 20


class PreflightError(RuntimeError):
    """The observation does not establish the frozen pre-lane conditions."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise PreflightError(message)


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=True, allow_nan=False).encode("ascii")


def exact_keys(value: Any, keys: set[str], label: str) -> dict:
    require(type(value) is dict and len(value) == len(keys) and set(value) == keys,
            f"{label} has unexpected or missing fields")
    return value


def integer(value: Any, minimum: int, maximum: int, label: str) -> int:
    require(type(value) is int and minimum <= value <= maximum,
            f"{label} is not a bounded integer")
    return value


def decimal(text: str, label: str) -> int:
    value = text.strip()
    require(re.fullmatch(r"0|[1-9][0-9]{0,18}", value) is not None,
            f"{label} is not a canonical unsigned decimal")
    return int(value)


def cpu_list(text: str) -> list[int]:
    values: list[int] = []
    for token in text.strip().split(","):
        match = re.fullmatch(r"(0|[1-9][0-9]{0,2})(?:-(0|[1-9][0-9]{0,2}))?", token)
        require(match is not None, "malformed CPU list")
        first = int(match[1])
        last = int(match[2]) if match[2] is not None else first
        require(0 <= first <= last < 128, "CPU list is outside the frozen host")
        values.extend(range(first, last + 1))
    require(values and values == sorted(set(values)),
            "CPU list is duplicated or out of order")
    return values


def canonical_path(value: str, *, nonroot: bool = False) -> str:
    if value == "/":
        require(not nonroot, "root cgroup is not a bounded scope")
        return value
    require(type(value) is str and len(value) <= 4096 and
            value.startswith("/") and "\x00" not in value and
            all(part not in ("", ".", "..") for part in value.split("/")[1:]),
            "unsafe or noncanonical absolute path")
    require(not nonroot or value != "/", "root cgroup is not a bounded scope")
    return value


def cgroup_location(membership: str, mountinfo: str) -> tuple[str, str]:
    require(type(membership) is str and type(mountinfo) is str and
            len(membership) <= 4096 and len(mountinfo) <= MAX_TEXT,
            "cgroup membership or mountinfo is unbounded")
    rows = membership.splitlines()
    require(len(rows) == 1 and rows[0].startswith("0::"),
            "exact unified cgroup membership is required")
    member = canonical_path(rows[0][3:], nonroot=True)
    mounts = []
    for line in mountinfo.splitlines():
        halves = line.split(" - ")
        if len(halves) != 2:
            raise PreflightError("malformed mountinfo row")
        left, right = halves[0].split(), halves[1].split()
        require(len(left) >= 6 and len(right) >= 3, "short mountinfo row")
        if right[0] == "cgroup2":
            decimal(left[0], "mount ID")
            decimal(left[1], "parent mount ID")
            require(re.fullmatch(r"[0-9]+:[0-9]+", left[2]) is not None and
                    "memory_localevents" not in (left[5] + "," + right[2]).split(","),
                    "cgroup2 mount device or event semantics differ")
            require(left[3] == "/" and left[4] == CGROUP_MOUNT,
                    "cgroup2 mount is rebased or noncanonical")
            mounts.append(line)
    require(len(mounts) == 1, "cgroup2 mount is absent or ambiguous")
    return member, mounts[0]


def memory_limit(text: str) -> int | str:
    return "max" if text.strip() == "max" else decimal(text, "memory limit")


def memory_events(text: str) -> dict[str, int]:
    events: dict[str, int] = {}
    for line in text.splitlines():
        fields = line.split()
        require(len(events) < 64 and len(fields) == 2 and
                re.fullmatch(r"[a-z_]{1,64}", fields[0]) is not None,
                "malformed memory.events row")
        require(fields[0] not in events, "duplicate memory.events key")
        events[fields[0]] = decimal(fields[1], "memory event")
    require({"low", "high", "max", "oom", "oom_kill"}.issubset(events) and
            all(value == 0 for value in events.values()),
            "cgroup memory events are incomplete or not clean")
    return events


def load_preregistration(data: bytes) -> dict:
    require(type(data) is bytes and len(data) <= MAX_TEXT and
            hashlib.sha256(data).hexdigest() == PREREGISTRATION_SHA256,
            "preregistration bytes differ from frozen NOT ARMED authority")
    value = json.loads(data)
    require(value["live_acquisition_armed"] is False and
            value["host_authority"]["hostname"] == "ripper" and
            value["resource_envelope"]["memory_max_bytes"] == MEMORY_MAX,
            "frozen host/resource preregistration differs")
    return value


def identity(value: Any, label: str) -> dict:
    value = exact_keys(value, {"device", "inode"}, label)
    integer(value["device"], 0, (1 << 64) - 1, label + " device")
    integer(value["inode"], 1, (1 << 64) - 1, label + " inode")
    return value


class LinuxReader:
    """Only reads bounded kernel files; never shells out or changes state."""

    def __init__(self) -> None:
        require(sys.platform.startswith("linux"), "v19 requires Linux")
        self.file_identities: dict[str, dict] = {}

    @staticmethod
    def open_directory(path: str) -> int:
        canonical_path(path)
        flags = os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW
        descriptor = os.open("/", flags)
        try:
            for part in (() if path == "/" else path.split("/")[1:]):
                next_descriptor = os.open(part, flags, dir_fd=descriptor)
                os.close(descriptor)
                descriptor = next_descriptor
            result, descriptor = descriptor, -1
            return result
        finally:
            if descriptor >= 0:
                os.close(descriptor)

    @staticmethod
    def stat_identity(value: os.stat_result) -> dict:
        return {"device": value.st_dev, "inode": value.st_ino}

    def read_text(self, path: str, limit: int = 4096) -> str:
        canonical_path(path)
        require(type(limit) is int and 0 < limit <= MAX_TEXT, "invalid read bound")
        parent, name = path.rsplit("/", 1)
        parent = parent or "/"
        directory = self.open_directory(parent)
        descriptor = -1
        try:
            descriptor = os.open(name, os.O_RDONLY | os.O_CLOEXEC |
                                 os.O_NOFOLLOW | os.O_NONBLOCK, dir_fd=directory)
            before = os.fstat(descriptor)
            require(stat.S_ISREG(before.st_mode), "kernel input is not a regular file")
            data = bytearray()
            while True:
                part = os.read(descriptor, limit + 1 - len(data))
                if not part:
                    break
                data.extend(part)
                require(len(data) <= limit, "kernel input exceeds read bound")
            after = os.fstat(descriptor)
            visible_directory = self.open_directory(parent)
            try:
                require(self.stat_identity(os.fstat(directory)) ==
                        self.stat_identity(os.fstat(visible_directory)),
                        "kernel input parent identity changed during read")
                named = os.stat(name, dir_fd=visible_directory, follow_symlinks=False)
            finally:
                os.close(visible_directory)
            require(stat.S_ISREG(named.st_mode) and
                    (before.st_dev, before.st_ino, before.st_mode) ==
                    (after.st_dev, after.st_ino, after.st_mode) ==
                    (named.st_dev, named.st_ino, named.st_mode),
                    "kernel input identity changed during read")
            self.file_identities[path] = self.stat_identity(before)
            return data.decode("ascii", errors="strict")
        finally:
            if descriptor >= 0:
                os.close(descriptor)
            os.close(directory)

    def directory_identity(self, path: str) -> dict:
        descriptor = self.open_directory(path)
        try:
            return self.stat_identity(os.fstat(descriptor))
        finally:
            os.close(descriptor)

    def list_cache_indices(self, path: str) -> list[str]:
        descriptor = self.open_directory(path)
        try:
            names = []
            with os.scandir(descriptor) as entries:
                for entry in entries:
                    require(len(names) < 32, "cache directory exceeds entry bound")
                    names.append(entry.name)
            indices = sorted(name for name in names if name.startswith("index"))
            require(indices and all(re.fullmatch(r"index(?:0|[1-9][0-9]?)", name)
                                    for name in indices), "invalid cache index names")
            return indices
        finally:
            os.close(descriptor)

    def process(self) -> dict:
        pid = os.getpid()
        directory = self.open_directory(f"/proc/{pid}/ns")
        try:
            namespaces = {}
            for name in ("pid", "user", "mnt", "cgroup", "uts"):
                # These particular kernel namespace links must be followed.
                descriptor = os.open(name, os.O_RDONLY | os.O_CLOEXEC,
                                     dir_fd=directory)
                try:
                    current = os.fstat(descriptor)
                    named = os.stat(name, dir_fd=directory)
                    require(self.stat_identity(current) == self.stat_identity(named),
                            "namespace identity changed during read")
                    namespaces[name] = self.stat_identity(current)
                finally:
                    os.close(descriptor)
        finally:
            os.close(directory)
        uname = os.uname()
        soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        return {
            "pid": pid, "uid": os.getuid(), "euid": os.geteuid(),
            "system": uname.sysname, "hostname": uname.nodename,
            "machine": uname.machine, "kernel_release": uname.release,
            "clock_ticks_per_second": os.sysconf("SC_CLK_TCK"),
            "allowed_cpus": sorted(os.sched_getaffinity(0)),
            "nofile_soft": soft, "nofile_hard": hard,
            "namespaces": namespaces,
        }


def validate_process(value: Any) -> None:
    value = exact_keys(value, {
        "pid", "uid", "euid", "system", "hostname", "machine", "kernel_release",
        "clock_ticks_per_second", "allowed_cpus", "nofile_soft", "nofile_hard",
        "namespaces"}, "process")
    integer(value["pid"], 1, (1 << 31) - 1, "process PID")
    integer(value["uid"], 0, (1 << 32) - 1, "process UID")
    require(type(value["euid"]) is int and value["euid"] == value["uid"] and
            value["system"] == "Linux" and value["hostname"] == "ripper" and
            value["machine"] == "x86_64" and
            type(value["kernel_release"]) is str and
            0 < len(value["kernel_release"]) <= 256,
            "process is not the native ripper host")
    require(type(value["clock_ticks_per_second"]) is int and
            value["clock_ticks_per_second"] == 100,
            "host clock tick rate differs")
    require(type(value["allowed_cpus"]) is list and len(value["allowed_cpus"]) == 128 and
            all(type(cpu) is int for cpu in value["allowed_cpus"]) and
            value["allowed_cpus"] == CPUS, "launch affinity is not the full frozen host")
    soft, hard = value["nofile_soft"], value["nofile_hard"]
    require(type(soft) is int and type(hard) is int and
            (soft == -1 or 65536 <= soft < (1 << 63)) and
            (hard == -1 or soft != -1 and soft <= hard < (1 << 63)),
            "RLIMIT_NOFILE is below the frozen build requirement")
    namespaces = exact_keys(value["namespaces"], {"pid", "user", "mnt", "cgroup", "uts"},
                            "process namespaces")
    for name, item in namespaces.items():
        identity(item, name + " namespace")


def cpu_models(text: str) -> list[str]:
    models = {}
    for block in text.strip().split("\n\n"):
        retained = {}
        for line in block.splitlines():
            fields = line.split(":", 1)
            if len(fields) == 2 and fields[0].strip() in ("processor", "model name"):
                key = fields[0].strip()
                require(key not in retained, "duplicate CPU identity field")
                retained[key] = fields[1].strip()
        require(set(retained) == {"processor", "model name"}, "CPU identity is incomplete")
        cpu = decimal(retained["processor"], "processor ID")
        require(cpu not in models and cpu < 128, "duplicate or unexpected processor")
        models[cpu] = retained["model name"]
    require(sorted(models) == CPUS, "CPU model inventory is incomplete")
    return [models[cpu] for cpu in CPUS]


def collect_observation(reader: LinuxReader) -> dict:
    process = reader.process()
    validate_process(process)  # Reject the wrong host before the detailed scan.
    pid = process["pid"]
    read = reader.read_text
    member, mount = cgroup_location(read(f"/proc/{pid}/cgroup"),
                                    read(f"/proc/{pid}/mountinfo", MAX_TEXT))
    scope = CGROUP_MOUNT + member
    scope_identity = reader.directory_identity(scope)
    cgroup_type = read(scope + "/cgroup.type").strip()
    require(cgroup_type == "domain", "cgroup is not a domain")
    ancestors = []
    parent = str(PurePosixPath(scope).parent)
    while parent != CGROUP_MOUNT:
        require(len(ancestors) < 128, "cgroup ancestor depth exceeds bound")
        require(parent.startswith(CGROUP_MOUNT + "/"), "cgroup parent escaped its mount")
        ancestors.append({"path": parent, "identity": reader.directory_identity(parent),
                          "memory_max": memory_limit(read(parent + "/memory.max"))})
        parent = str(PurePosixPath(parent).parent)
    resources = {
        "membership": member, "mountinfo": mount, "scope_identity": scope_identity,
        "cgroup_type": cgroup_type,
        "memory_max": memory_limit(read(scope + "/memory.max")),
        "memory_swap_max": memory_limit(read(scope + "/memory.swap.max")),
        "memory_swap_current": decimal(read(scope + "/memory.swap.current"), "swap current"),
        "memory_events": memory_events(read(scope + "/memory.events")),
        "ancestors": ancestors,
    }
    swaps = read("/proc/swaps").splitlines()
    require(len(swaps) == 1 and swaps[0].split() ==
            ["Filename", "Type", "Size", "Used", "Priority"], "host swap is enabled or unreadable")
    topology = []
    for cpu in CPUS:
        root = f"{SYS_CPU}/cpu{cpu}"
        domains = []
        for index in reader.list_cache_indices(root + "/cache"):
            cache = root + "/cache/" + index
            if decimal(read(cache + "/level"), "cache level") == 3:
                require(read(cache + "/type").strip() == "Unified", "L3 cache is not unified")
                domains.append(cpu_list(read(cache + "/shared_cpu_list")))
        require(len(domains) == 1, "CPU has absent or ambiguous L3 cache")
        topology.append({
            "cpu": cpu,
            "package": decimal(read(root + "/topology/physical_package_id"), "package"),
            "die": decimal(read(root + "/topology/die_id"), "die"),
            "core": decimal(read(root + "/topology/core_id"), "core"),
            "siblings": cpu_list(read(root + "/topology/thread_siblings_list")),
            "l3_cpus": domains[0],
        })
    require(reader.directory_identity(scope) == scope_identity,
            "cgroup directory changed during capture")
    return {
        "process": process,
        "boot_id": read("/proc/sys/kernel/random/boot_id").strip(),
        "cpu_models": cpu_models(read("/proc/cpuinfo", MAX_TEXT)),
        "online_cpus": cpu_list(read(SYS_CPU + "/online")),
        "online_nodes": cpu_list(read(SYS_NODE + "/online")),
        "node_cpus": cpu_list(read(SYS_NODE + "/node0/cpulist")),
        "topology": topology, "resources": resources, "host_swap_enabled": False,
    }


def validate_observation(value: Any) -> None:
    value = exact_keys(value, {"process", "boot_id", "cpu_models", "online_cpus",
                              "online_nodes", "node_cpus", "topology", "resources",
                              "host_swap_enabled"}, "observation")
    validate_process(value["process"])
    require(type(value["boot_id"]) is str and re.fullmatch(
        r"[0-9a-f]{8}(?:-[0-9a-f]{4}){3}-[0-9a-f]{12}", value["boot_id"]) is not None,
        "host boot identity differs")
    for key, expected in (("online_cpus", CPUS), ("online_nodes", [0]), ("node_cpus", CPUS)):
        require(type(value[key]) is list and len(value[key]) == len(expected) and
                all(type(item) is int for item in value[key]) and
                value[key] == expected, f"host {key} differs")
    require(type(value["cpu_models"]) is list and len(value["cpu_models"]) == 128 and
            value["cpu_models"] == [CPU_MODEL] * 128 and value["host_swap_enabled"] is False,
            "host model or swap differs")
    topology = value["topology"]
    require(type(topology) is list and len(topology) == 128, "topology size differs")
    for cpu, item in enumerate(topology):
        exact_keys(item, {"cpu", "package", "die", "core", "siblings", "l3_cpus"}, "CPU topology")
        for name in ("cpu", "package", "die", "core"):
            integer(item[name], 0, (1 << 31) - 1, "CPU " + name)
        require(item["cpu"] == cpu and type(item["siblings"]) is list and len(item["siblings"]) == 2 and
                all(type(x) is int for x in item["siblings"]) and
                item["siblings"] == [cpu % 64, cpu % 64 + 64], "SMT sibling geometry differs")
        domain = item["l3_cpus"]
        require(type(domain) is list and len(domain) == 16 and
                all(type(x) is int and 0 <= x < 128 for x in domain) and
                domain == sorted(set(domain)) and cpu in domain,
                "L3 domain shape differs")
    cores = {(x["package"], x["die"], x["core"]) for x in topology}
    require(len(cores) == 64 and len({x["package"] for x in topology}) == 1 and
            len({tuple(x["l3_cpus"]) for x in topology}) == 8, "host core/package/L3 counts differ")
    for item in topology:
        for cpu in item["siblings"]:
            other = topology[cpu]
            require(all(other[key] == item[key] for key in ("package", "die", "core", "l3_cpus")),
                    "SMT physical identity differs")
        require(all(topology[cpu]["l3_cpus"] == item["l3_cpus"] for cpu in item["l3_cpus"]),
                "L3 domain membership is asymmetric")
        require(all(topology[cpu][key] == item[key] for cpu in item["l3_cpus"]
                    for key in ("package", "die")), "L3 domain crosses a package or die")
    resources = exact_keys(value["resources"], {"membership", "mountinfo", "scope_identity",
        "memory_max", "memory_swap_max", "memory_swap_current", "memory_events", "ancestors",
        "cgroup_type"},
        "resources")
    require(type(resources["membership"]) is str and len(resources["membership"]) <= 4096 and
            resources["cgroup_type"] == "domain",
            "cgroup membership or domain type differs")
    member, mount = cgroup_location("0::" + resources["membership"], resources["mountinfo"])
    require(mount == resources["mountinfo"], "retained cgroup mount row is not canonical")
    identity(resources["scope_identity"], "cgroup scope")
    require(type(resources["memory_max"]) is int and resources["memory_max"] == MEMORY_MAX and
            type(resources["memory_swap_max"]) is int and resources["memory_swap_max"] == 0 and
            type(resources["memory_swap_current"]) is int and resources["memory_swap_current"] == 0,
            "cgroup hard memory/swap limits differ")
    events = resources["memory_events"]
    require(type(events) is dict and len(events) <= 64 and
            all(type(key) is str and len(key) <= 64 and type(v) is int and 0 <= v < (1 << 63)
                for key, v in events.items()), "memory events types or bounds differ")
    memory_events("\n".join(f"{key} {v}" for key, v in events.items()))
    ancestors = resources["ancestors"]
    require(type(ancestors) is list and len(ancestors) <= 128, "ancestor list differs")
    parent = str(PurePosixPath(CGROUP_MOUNT + member).parent)
    for item in ancestors:
        exact_keys(item, {"path", "identity", "memory_max"}, "ancestor")
        require(parent != CGROUP_MOUNT and item["path"] == parent, "ancestor path chain differs")
        identity(item["identity"], "ancestor")
        limit = item["memory_max"]
        require(limit == "max" or type(limit) is int and MEMORY_MAX <= limit < (1 << 63),
                "ancestor tightens the frozen memory limit")
        parent = str(PurePosixPath(parent).parent)
    require(parent == CGROUP_MOUNT, "ancestor chain is incomplete")


def validate_capture(value: Any, preregistration_bytes: bytes) -> dict:
    load_preregistration(preregistration_bytes)
    value = exact_keys(value, {"schema", "preregistration_sha256", "started_monotonic_ns",
        "finished_monotonic_ns", "observations", "file_identity_sha256", "atomic_snapshot",
        "resource_lifetime_proved", "live_acquisition_armed", "digest"}, "capture")
    require(type(value["digest"]) is str and re.fullmatch(r"[0-9a-f]{64}", value["digest"]) is not None,
            "capture digest is malformed")
    require(value["schema"] == SCHEMA and value["preregistration_sha256"] == PREREGISTRATION_SHA256 and
            value["atomic_snapshot"] is False and value["resource_lifetime_proved"] is False and
            value["live_acquisition_armed"] is False, "capture claim ceiling differs")
    started = integer(value["started_monotonic_ns"], 0, (1 << 63) - 1, "capture start")
    integer(value["finished_monotonic_ns"], started, (1 << 63) - 1, "capture finish")
    require(type(value["file_identity_sha256"]) is str and re.fullmatch(
        r"[0-9a-f]{64}", value["file_identity_sha256"]) is not None, "file identity digest differs")
    observations = value["observations"]
    require(type(observations) is list and len(observations) == 2, "capture needs two observations")
    for observation in observations:
        validate_observation(observation)
    # Bound every variable-size field before allocating its JSON encoding.
    payload = {key: item for key, item in value.items() if key != "digest"}
    data = canonical_bytes(payload)
    require(len(data) + 80 <= MAX_RECORD_BYTES, "capture exceeds size bound")
    require(hashlib.sha256(data).hexdigest() == value["digest"], "capture digest differs")
    require(canonical_bytes(observations[0]) == canonical_bytes(observations[1]),
            "host or resource identity changed between observations")
    return copy.deepcopy(value)


def capture_host(preregistration_bytes: bytes, reader: LinuxReader | None = None) -> dict:
    load_preregistration(preregistration_bytes)
    reader = LinuxReader() if reader is None else reader
    started = time.monotonic_ns()
    observations, identities = [], []
    for _ in range(2):
        reader.file_identities = {}
        observation = collect_observation(reader)
        validate_observation(observation)
        observations.append(observation)
        identities.append(copy.deepcopy(reader.file_identities))
    require(identities[0] == identities[1], "kernel file identities changed between observations")
    payload = {
        "schema": SCHEMA, "preregistration_sha256": PREREGISTRATION_SHA256,
        "started_monotonic_ns": started, "finished_monotonic_ns": time.monotonic_ns(),
        "observations": observations,
        "file_identity_sha256": hashlib.sha256(canonical_bytes(identities[0])).hexdigest(),
        "atomic_snapshot": False, "resource_lifetime_proved": False,
        "live_acquisition_armed": False,
    }
    return validate_capture(dict(payload, digest=hashlib.sha256(canonical_bytes(payload)).hexdigest()),
                            preregistration_bytes)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preregistration", required=True,
                        help="absolute path to the exact NOT ARMED wrapper preregistration")
    options = parser.parse_args()
    try:
        reader = LinuxReader()
        data = reader.read_text(options.preregistration, MAX_TEXT).encode("ascii")
        result = capture_host(data, reader)
        sys.stdout.buffer.write(canonical_bytes(result) + b"\n")
        return 0
    except (PreflightError, OSError, ValueError, TypeError, KeyError) as error:
        print(f"v19 host preflight rejected: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
