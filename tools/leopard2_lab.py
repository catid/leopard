#!/usr/bin/env python3
"""Deterministic, resumable experiment runner for Leopard2.

The runner deliberately uses only the Python standard library.  It records the
CPU topology used to construct a manifest, assigns stable seeds and CPU sets,
applies per-process address-space limits, and captures every job independently.

Run ``python3 tools/leopard2_lab.py --help`` for the command-line interface.
"""

from __future__ import print_function

import argparse
import hashlib
import json
import math
import os
import posixpath
import re
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path

_TOOLS_DIRECTORY = str(Path(__file__).resolve().parent)
if _TOOLS_DIRECTORY not in sys.path:
    sys.path.insert(0, _TOOLS_DIRECTORY)

from leopard2_perf_evidence import (  # noqa: E402
    derive_perf_stat as _derive_perf_stat,
    parse_perf_stat as _parse_perf_stat,
    perf_probe_command as _perf_probe_command,
    probe_command_matches_request as _probe_command_matches_request,
    read_perf_stat_evidence as _read_perf_stat_evidence,
)

try:
    import resource
except ImportError:  # pragma: no cover - Unix/Linux is the production target.
    resource = None


MANIFEST_SCHEMA = "leopard2-lab-manifest/v5"
RESULT_SCHEMA = "leopard2-lab-result/v3"
MERGE_SCHEMA = "leopard2-lab-merged/v3"
JOB_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
TOKEN_REPLACEMENTS = ("seed", "job_id", "cpu_set")
PERF_PROVIDER = "linux-perf-stat"
PERF_EVENT_RE = re.compile(r"^[A-Za-z0-9_.:/=-]+$")
THREAD_POLICY_SCHEMA = "leopard2-lab-thread-policy/v1"
THREAD_OBSERVATION_SCHEMA = "leopard2-lab-thread-observation/v2"

# These variables cover OpenMP and the common nested native thread pools that
# can otherwise multiply a one-process-per-CPU fuzz campaign into hundreds of
# runnable threads.  The complete, effective subset is signed into each job;
# inherited values never leak into a child unnoticed.
THREAD_COUNT_ENV = (
    "OMP_NUM_THREADS",
    "OMP_THREAD_LIMIT",
    "OPENBLAS_NUM_THREADS",
    "GOTO_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "RAYON_NUM_THREADS",
)
THREAD_FALSE_ENV = (
    "OMP_DYNAMIC",
    "OMP_NESTED",
    "MKL_DYNAMIC",
)
THREAD_FIXED_ENV = {
    "OMP_MAX_ACTIVE_LEVELS": "1",
}
THREAD_ENV_KEYS = tuple(sorted(
    THREAD_COUNT_ENV + THREAD_FALSE_ENV + tuple(THREAD_FIXED_ENV)))


class LabError(Exception):
    """An actionable input, platform, or execution error."""


def _read_text(path):
    try:
        return Path(path).read_text(encoding="utf-8").strip()
    except (OSError, UnicodeError):
        return None


def _read_int(path):
    value = _read_text(path)
    if value is None:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def parse_cpu_list(text):
    """Parse Linux cpulist syntax such as ``0-3,8,10-12``."""
    cpus = set()
    if not text:
        return []
    for raw_part in text.split(","):
        part = raw_part.strip()
        if not part:
            continue
        if "-" in part:
            endpoints = part.split("-", 1)
            try:
                first, last = int(endpoints[0]), int(endpoints[1])
            except ValueError:
                raise LabError("invalid CPU range: {!r}".format(part))
            if first < 0 or last < first:
                raise LabError("invalid CPU range: {!r}".format(part))
            cpus.update(range(first, last + 1))
        else:
            try:
                cpu = int(part)
            except ValueError:
                raise LabError("invalid CPU number: {!r}".format(part))
            if cpu < 0:
                raise LabError("CPU numbers must be non-negative")
            cpus.add(cpu)
    return sorted(cpus)


def format_cpu_list(cpus):
    """Return a stable, compact Linux-style representation of CPU numbers."""
    values = sorted(set(int(cpu) for cpu in cpus))
    if not values:
        return ""
    ranges = []
    start = previous = values[0]
    for value in values[1:]:
        if value == previous + 1:
            previous = value
            continue
        ranges.append(str(start) if start == previous else "{}-{}".format(start, previous))
        start = previous = value
    ranges.append(str(start) if start == previous else "{}-{}".format(start, previous))
    return ",".join(ranges)


def _validate_canonical_cpu_set(value, label, allowed=None):
    """Require a non-empty, sorted, duplicate-free list of CPU identifiers."""
    if (not isinstance(value, list) or not value or
            any(isinstance(cpu, bool) or not isinstance(cpu, int) or cpu < 0
                for cpu in value)):
        raise LabError(
            "{} must be a non-empty list of non-negative CPU integers".format(
                label))
    if value != sorted(set(value)):
        raise LabError(
            "{} must be sorted and contain no duplicate CPUs".format(label))
    if allowed is not None:
        unavailable = set(value) - set(allowed)
        if unavailable:
            raise LabError(
                "{} contains unavailable CPUs {}".format(
                    label, format_cpu_list(unavailable)))
    return value


def _allowed_cpus():
    if hasattr(os, "sched_getaffinity"):
        try:
            allowed = sorted(os.sched_getaffinity(0))
            if allowed:
                return allowed, "sched_getaffinity"
        except OSError:
            pass
    count = os.cpu_count() or 1
    return list(range(count)), "os.cpu_count"


def _memory_total_bytes(proc_root="/proc"):
    text = _read_text(Path(proc_root) / "meminfo")
    if text:
        for line in text.splitlines():
            fields = line.split()
            if len(fields) >= 2 and fields[0] == "MemTotal:":
                try:
                    return int(fields[1]) * 1024
                except ValueError:
                    break
    return None


def _decode_mount_field(value):
    """Decode the octal escapes used by Linux mountinfo path fields."""
    return re.sub(
        r"\\([0-7]{3})",
        lambda match: chr(int(match.group(1), 8)),
        value)


def _cgroup_memberships(proc_root):
    text = _read_text(Path(proc_root) / "self" / "cgroup")
    if not text:
        return []
    memberships = []
    for line in text.splitlines():
        fields = line.split(":", 2)
        if len(fields) != 3:
            continue
        controllers = set(filter(None, fields[1].split(",")))
        path = posixpath.normpath(fields[2])
        if not path.startswith("/"):
            continue
        memberships.append((controllers, path))
    return memberships


def _cgroup_mounts(proc_root):
    text = _read_text(Path(proc_root) / "self" / "mountinfo")
    if not text:
        return []
    mounts = []
    for line in text.splitlines():
        try:
            before, after = line.split(" - ", 1)
        except ValueError:
            continue
        before_fields = before.split()
        after_fields = after.split()
        if len(before_fields) < 5 or len(after_fields) < 3:
            continue
        filesystem = after_fields[0]
        if filesystem not in ("cgroup", "cgroup2"):
            continue
        mounts.append({
            "filesystem": filesystem,
            "root": posixpath.normpath(_decode_mount_field(before_fields[3])),
            "mount_point": Path(_decode_mount_field(before_fields[4])),
            "controllers": set(after_fields[2].split(",")),
        })
    return mounts


def _mounted_cgroup_path(mount, membership):
    mount_root = mount["root"]
    if mount_root == "/":
        relative = membership.lstrip("/")
    elif membership == mount_root:
        relative = ""
    elif membership.startswith(mount_root.rstrip("/") + "/"):
        relative = membership[len(mount_root.rstrip("/")) + 1:]
    else:
        return None
    candidate = mount["mount_point"] / relative
    mount_point = mount["mount_point"]
    if candidate != mount_point and mount_point not in candidate.parents:
        return None
    return candidate


def _read_cgroup_limits(proc_root):
    memberships = _cgroup_memberships(proc_root)
    mounts = _cgroup_mounts(proc_root)
    limits = []
    visited = set()
    for mount in mounts:
        if mount["filesystem"] == "cgroup2":
            matching = [path for controllers, path in memberships if not controllers]
            limit_name = "memory.max"
        elif "memory" in mount["controllers"]:
            matching = [path for controllers, path in memberships if "memory" in controllers]
            limit_name = "memory.limit_in_bytes"
        else:
            continue
        for membership in matching:
            current = _mounted_cgroup_path(mount, membership)
            if current is None:
                continue
            while True:
                limit_path = current / limit_name
                if limit_path not in visited:
                    visited.add(limit_path)
                    text = _read_text(limit_path)
                    if text and text != "max":
                        try:
                            value = int(text)
                        except ValueError:
                            value = 0
                        # cgroup v1 uses enormous values as unlimited sentinels.
                        if 0 < value < (1 << 60):
                            limits.append(value)
                if current == mount["mount_point"]:
                    break
                current = current.parent
    return limits


def _memory_capacity_bytes(proc_root="/proc"):
    """Return the smaller of physical memory and readable cgroup limits."""
    candidates = []
    total = _memory_total_bytes(proc_root)
    if total:
        candidates.append(total)
    candidates.extend(_read_cgroup_limits(proc_root))
    return min(candidates) if candidates else None


def detect_topology():
    """Detect the process-visible Linux CPU and NUMA topology without root."""
    allowed, affinity_source = _allowed_cpus()
    allowed_set = set(allowed)

    numa_for_cpu = {}
    numa_nodes = []
    node_root = Path("/sys/devices/system/node")
    if node_root.is_dir():
        for node_path in sorted(node_root.glob("node[0-9]*"), key=lambda p: int(p.name[4:])):
            node_id = int(node_path.name[4:])
            try:
                node_cpus = parse_cpu_list(_read_text(node_path / "cpulist") or "")
            except LabError:
                node_cpus = []
            node_cpus = sorted(allowed_set.intersection(node_cpus))
            if not node_cpus:
                continue
            numa_nodes.append({"node": node_id, "cpus": node_cpus})
            for cpu in node_cpus:
                numa_for_cpu[cpu] = node_id

    cpu_entries = []
    core_keys = set()
    sockets = set()
    for cpu in allowed:
        topology_root = Path("/sys/devices/system/cpu/cpu{}".format(cpu)) / "topology"
        socket_id = _read_int(topology_root / "physical_package_id")
        core_id = _read_int(topology_root / "core_id")
        if socket_id is None:
            socket_id = 0
        if core_id is None:
            core_id = cpu
        sockets.add(socket_id)
        core_keys.add((socket_id, core_id))
        entry = {"cpu": cpu, "socket": socket_id, "core": core_id}
        if cpu in numa_for_cpu:
            entry["numa_node"] = numa_for_cpu[cpu]
        cpu_entries.append(entry)

    if not numa_nodes:
        # Unknown NUMA topology is represented explicitly rather than pretending
        # that the system has one node.
        numa_count = None
    else:
        numa_count = len(numa_nodes)

    return {
        "affinity_source": affinity_source,
        "allowed_cpus": allowed,
        "allowed_cpu_list": format_cpu_list(allowed),
        "logical_cpus": len(allowed),
        "physical_cores": len(core_keys),
        "sockets": len(sockets),
        "numa_node_count": numa_count,
        "numa_nodes": numa_nodes,
        "cpus": cpu_entries,
        "memory_total_bytes": _memory_total_bytes(),
        "memory_capacity_bytes": _memory_capacity_bytes(),
    }


def _canonical_json_bytes(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")


def _digest(value):
    return hashlib.sha256(_canonical_json_bytes(value)).hexdigest()


def _new_run_epoch():
    """Return an opaque identity shared by results from one runner call."""
    return hashlib.sha256(os.urandom(32)).hexdigest()


def _file_identity(path):
    """Return a content identity for an executable, independent of mtime."""
    resolved = Path(path).resolve()
    try:
        if not resolved.is_file():
            raise LabError("executable is not a regular file: {}".format(resolved))
        hasher = hashlib.sha256()
        with resolved.open("rb") as source:
            before = os.fstat(source.fileno())
            while True:
                block = source.read(1024 * 1024)
                if not block:
                    break
                hasher.update(block)
            stat_result = os.fstat(source.fileno())
        if (before.st_ino != stat_result.st_ino or
                before.st_size != stat_result.st_size or
                before.st_mtime_ns != stat_result.st_mtime_ns):
            raise LabError("executable changed while it was being hashed: {}".format(resolved))
    except OSError as error:
        raise LabError("cannot hash executable {}: {}".format(resolved, error))
    return {
        "path": str(resolved),
        "sha256": hasher.hexdigest(),
        "size_bytes": stat_result.st_size,
    }


def _content_identity(path):
    path = Path(path)
    try:
        hasher = hashlib.sha256()
        size = 0
        with path.open("rb") as source:
            while True:
                block = source.read(1024 * 1024)
                if not block:
                    break
                hasher.update(block)
                size += len(block)
    except OSError as error:
        raise LabError("cannot hash output {}: {}".format(path, error))
    return {"sha256": hasher.hexdigest(), "size_bytes": size}


def _resolve_executable(command, root_path, cwd, environment, identity_cache=None):
    token = command[0]
    if any("{" + replacement + "}" in token for replacement in TOKEN_REPLACEMENTS):
        raise LabError("the command executable may not contain runtime tokens")
    working_directory = Path(cwd)
    if not working_directory.is_absolute():
        working_directory = root_path / working_directory
    if os.path.sep in token or (os.path.altsep and os.path.altsep in token):
        candidate = Path(token)
        if not candidate.is_absolute():
            candidate = working_directory / candidate
        executable = str(candidate)
    else:
        search_path = environment.get("PATH", os.environ.get("PATH", ""))
        search_path = os.pathsep.join(
            entry if Path(entry or ".").is_absolute()
            else str(working_directory / (entry or "."))
            for entry in search_path.split(os.pathsep))
        executable = shutil.which(token, path=search_path)
        if executable is None:
            raise LabError("cannot resolve executable {!r}".format(token))
    if not os.access(executable, os.X_OK):
        raise LabError("command executable is not executable: {}".format(executable))
    cache_key = str(Path(executable).resolve())
    if identity_cache is not None and cache_key in identity_cache:
        return dict(identity_cache[cache_key])
    identity = _file_identity(executable)
    if identity_cache is not None:
        identity_cache[cache_key] = dict(identity)
    return identity


def _normalize_performance_counters(
        value, root_path, cwd, environment, identity_cache):
    """Validate and content-address an optional Linux perf-stat request."""
    if value is None:
        return None
    if not isinstance(value, dict):
        raise LabError("performance_counters must be an object")
    provider = value.get("provider", PERF_PROVIDER)
    if provider != PERF_PROVIDER:
        raise LabError("unsupported performance counter provider {!r}".format(provider))
    command = value.get("command", "perf")
    if not isinstance(command, str) or not command:
        raise LabError("performance counter command must be a non-empty string")
    events = value.get("events")
    if (not isinstance(events, list) or not events or
            not all(isinstance(event, str) and PERF_EVENT_RE.match(event)
                    for event in events)):
        raise LabError(
            "performance counter events must be a non-empty list of simple perf event names")
    if len(events) != len(set(events)):
        raise LabError("performance counter events must be unique")
    optional = value.get("optional", True)
    if not isinstance(optional, bool):
        raise LabError("performance counter optional must be boolean")
    executable = _resolve_executable(
        [command], root_path, cwd, environment, identity_cache)
    request = {
        "provider": provider,
        "events": list(events),
        "optional": optional,
        "executable": executable,
    }
    request["probe_command"] = _perf_probe_command(
        executable["path"], events, sys.executable)
    return request


def _json_copy(value):
    """Copy JSON-compatible metadata while rejecting non-JSON input early."""
    try:
        return json.loads(_canonical_json_bytes(value).decode("ascii"))
    except (TypeError, ValueError) as error:
        raise LabError("metadata must be JSON-serializable: {}".format(error))


def _atomic_write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=str(path.parent))
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as output:
            json.dump(value, output, indent=2, sort_keys=True, ensure_ascii=True)
            output.write("\n")
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary_name, str(path))
    except BaseException:
        try:
            os.unlink(temporary_name)
        except OSError:
            pass
        raise


def _load_json(path):
    try:
        with Path(path).open("r", encoding="utf-8") as source:
            return json.load(source)
    except OSError as error:
        raise LabError("cannot read {}: {}".format(path, error))
    except ValueError as error:
        raise LabError("invalid JSON in {}: {}".format(path, error))


def _positive_number(value, label, allow_zero=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise LabError("{} must be a number".format(label))
    if not math.isfinite(value) or value < 0 or (value == 0 and not allow_zero):
        comparator = "non-negative" if allow_zero else "positive"
        raise LabError("{} must be {}".format(label, comparator))
    return value


def _positive_int(value, label, allow_zero=False):
    if isinstance(value, bool) or not isinstance(value, int):
        raise LabError("{} must be an integer".format(label))
    _positive_number(value, label, allow_zero=allow_zero)
    return value


def _thread_environment(max_threads):
    """Return the deterministic nested-runtime cap for one workload job."""
    value = str(max_threads)
    environment = {key: value for key in THREAD_COUNT_ENV}
    environment.update({key: "FALSE" for key in THREAD_FALSE_ENV})
    environment.update(THREAD_FIXED_ENV)
    return dict(sorted(environment.items()))


def _normalize_thread_policy(value, cpu_set, label):
    """Validate an explicit internal-team request and bind its environment."""
    if value is None:
        value = {}
    if not isinstance(value, dict):
        raise LabError("{} thread_runtime must be an object".format(label))
    unknown = set(value) - {"max_threads", "allow_internal_team"}
    if unknown:
        raise LabError(
            "{} thread_runtime has unknown fields: {}".format(
                label, ", ".join(sorted(unknown))))
    max_threads = value.get("max_threads", 1)
    _positive_int(max_threads, "{} thread_runtime max_threads".format(label))
    allow_internal_team = value.get("allow_internal_team", False)
    if not isinstance(allow_internal_team, bool):
        raise LabError(
            "{} thread_runtime allow_internal_team must be boolean".format(
                label))
    if max_threads != 1 and not allow_internal_team:
        raise LabError(
            "{} requests {} internal threads without explicit "
            "allow_internal_team opt-in".format(label, max_threads))
    if max_threads > len(cpu_set):
        raise LabError(
            "{} requests {} internal threads but its CPU set contains only {} "
            "CPUs".format(label, max_threads, len(cpu_set)))
    return {
        "schema": THREAD_POLICY_SCHEMA,
        "max_threads": max_threads,
        "allow_internal_team": allow_internal_team,
        "effective_env": _thread_environment(max_threads),
    }


def _validate_thread_policy(policy, cpu_set, label):
    if not isinstance(policy, dict) or set(policy) != {
            "schema", "max_threads", "allow_internal_team", "effective_env"}:
        raise LabError("{} has an invalid thread runtime policy".format(label))
    if policy.get("schema") != THREAD_POLICY_SCHEMA:
        raise LabError("{} has an unsupported thread runtime policy".format(label))
    normalized = _normalize_thread_policy({
        "max_threads": policy.get("max_threads"),
        "allow_internal_team": policy.get("allow_internal_team"),
    }, cpu_set, label)
    if policy != normalized:
        raise LabError("{} thread runtime environment is not canonical".format(label))
    return policy


def _stable_seed(base_seed, job_id):
    material = "{}\0{}".format(base_seed, job_id).encode("utf-8")
    # Stay below 2^63 so seeds are accepted by libraries with signed APIs.
    return int.from_bytes(hashlib.sha256(material).digest()[:8], "big") & ((1 << 63) - 1)


def _cpu_order(topology, policy):
    entries = topology["cpus"]
    if policy == "logical":
        return [entry["cpu"] for entry in entries]
    if policy == "numa":
        seen = set()
        ordered = []
        for node in topology.get("numa_nodes", []):
            for cpu in node["cpus"]:
                if cpu not in seen:
                    ordered.append(cpu)
                    seen.add(cpu)
        ordered.extend(cpu for cpu in topology["allowed_cpus"] if cpu not in seen)
        return ordered
    if policy != "physical-first":
        raise LabError("unknown CPU policy {!r}".format(policy))

    by_core = {}
    core_order = []
    for entry in entries:
        key = (entry["socket"], entry["core"])
        if key not in by_core:
            by_core[key] = []
            core_order.append(key)
        by_core[key].append(entry["cpu"])
    ordered = []
    sibling_index = 0
    while True:
        added = False
        for key in core_order:
            siblings = by_core[key]
            if sibling_index < len(siblings):
                ordered.append(siblings[sibling_index])
                added = True
        if not added:
            break
        sibling_index += 1
    return ordered


def _expand_jobs(spec_jobs):
    expanded = []
    for raw in spec_jobs:
        if not isinstance(raw, dict):
            raise LabError("each job specification must be an object")
        repeat = raw.get("repeat", 1)
        _positive_int(repeat, "job repeat")
        base_id = raw.get("id")
        if not isinstance(base_id, str) or not JOB_ID_RE.match(base_id):
            raise LabError("invalid job id {!r}; use letters, digits, dot, dash, and underscore".format(base_id))
        for repetition in range(repeat):
            job = dict(raw)
            job.pop("repeat", None)
            if repeat > 1:
                job["id"] = "{}.{:04d}".format(base_id, repetition)
            expanded.append(job)
    ids = [job["id"] for job in expanded]
    if len(ids) != len(set(ids)):
        duplicates = sorted(job_id for job_id, count in Counter(ids).items() if count > 1)
        raise LabError("duplicate job ids: {}".format(", ".join(duplicates)))
    return sorted(expanded, key=lambda job: job["id"])


def _default_memory_mb(topology, workers):
    total = topology.get("memory_capacity_bytes")
    if total:
        per_worker = int((total // (1024 * 1024)) * 0.80 / workers)
        return max(256, min(4096, per_worker))
    return 1024


def build_manifest(spec, root=None, workers=None, base_seed=None):
    """Build and validate a deterministic manifest from a concise job spec."""
    if not isinstance(spec, dict):
        raise LabError("job specification must be a JSON object")
    source_spec = dict(spec)
    supplied_spec_digest = source_spec.pop("spec_digest", None)
    computed_spec_digest = _digest(source_spec)
    if supplied_spec_digest is not None and supplied_spec_digest != computed_spec_digest:
        raise LabError("source specification digest does not match its contents")
    spec_metadata = spec.get("metadata", {})
    if not isinstance(spec_metadata, dict):
        raise LabError("specification metadata must be an object")
    topology = detect_topology()
    if not topology["allowed_cpus"]:
        raise LabError("no CPUs are available in the process affinity mask")
    default_workers = min(topology["logical_cpus"], 128)
    selected_workers = workers if workers is not None else spec.get("workers", default_workers)
    _positive_int(selected_workers, "workers")
    if selected_workers > 128:
        raise LabError("workers may not exceed 128")

    selected_seed = base_seed if base_seed is not None else spec.get("base_seed", 1)
    _positive_int(selected_seed, "base_seed", allow_zero=True)
    root_path = Path(root if root is not None else spec.get("root", os.getcwd())).resolve()
    defaults = spec.get("defaults", {})
    if not isinstance(defaults, dict):
        raise LabError("defaults must be an object")
    default_performance_counters = defaults.get("performance_counters")
    default_timeout = defaults.get("timeout_seconds", 300.0)
    _positive_number(default_timeout, "default timeout_seconds")
    default_memory = defaults.get("memory_mb", _default_memory_mb(topology, selected_workers))
    _positive_int(default_memory, "default memory_mb", allow_zero=True)
    default_rss_limit = defaults.get("rss_limit_mb", 0)
    _positive_int(default_rss_limit, "default rss_limit_mb", allow_zero=True)
    default_minimum_memory = defaults.get("minimum_memory_mb", 0)
    _positive_int(
        default_minimum_memory, "default minimum_memory_mb", allow_zero=True)
    default_cpu_count = defaults.get("cpu_count", 1)
    _positive_int(default_cpu_count, "default cpu_count")
    default_thread_runtime = defaults.get("thread_runtime", {})
    if not isinstance(default_thread_runtime, dict):
        raise LabError("default thread_runtime must be an object")
    cpu_policy = defaults.get("cpu_policy", "physical-first")
    cpu_order = _cpu_order(topology, cpu_policy)
    allowed_set = set(topology["allowed_cpus"])

    raw_jobs = spec.get("jobs")
    if not isinstance(raw_jobs, list) or not raw_jobs:
        raise LabError("specification must contain a non-empty jobs array")
    jobs = []
    cursor = 0
    cpu_groups = {}
    executable_identities = {}
    for raw in _expand_jobs(raw_jobs):
        job_id = raw["id"]
        command = raw.get("command")
        if isinstance(command, str):
            command = shlex.split(command)
        if not isinstance(command, list) or not command or not all(isinstance(arg, str) and arg for arg in command):
            raise LabError("job {} command must be a non-empty string array".format(job_id))

        environment = raw.get("env", {})
        if not isinstance(environment, dict) or not all(isinstance(key, str) and isinstance(value, str)
                                                        for key, value in environment.items()):
            raise LabError("job {} env must map strings to strings".format(job_id))
        timeout_seconds = raw.get("timeout_seconds", default_timeout)
        memory_mb = raw.get("memory_mb", default_memory)
        rss_limit_mb = raw.get("rss_limit_mb", default_rss_limit)
        _positive_number(timeout_seconds, "job {} timeout_seconds".format(job_id))
        _positive_int(memory_mb, "job {} memory_mb".format(job_id), allow_zero=True)
        _positive_int(
            rss_limit_mb, "job {} rss_limit_mb".format(job_id),
            allow_zero=True)

        cpu_group = raw.get("cpu_group")
        if cpu_group is not None and (not isinstance(cpu_group, str) or not cpu_group):
            raise LabError("job {} cpu_group must be a non-empty string".format(job_id))
        resume_group = raw.get("resume_group")
        if resume_group is not None and (
                not isinstance(resume_group, str) or not resume_group):
            raise LabError(
                "job {} resume_group must be a non-empty string".format(job_id))
        if "cpu_set" in raw:
            raw_cpu_set = raw["cpu_set"]
            cpu_set = parse_cpu_list(raw_cpu_set) if isinstance(raw_cpu_set, str) else raw_cpu_set
            if not isinstance(cpu_set, list) or not cpu_set or not all(isinstance(cpu, int) for cpu in cpu_set):
                raise LabError("job {} cpu_set must be a non-empty CPU list".format(job_id))
            cpu_set = sorted(set(cpu_set))
            unavailable = set(cpu_set) - allowed_set
            if unavailable:
                raise LabError("job {} requests unavailable CPUs {}".format(job_id, format_cpu_list(unavailable)))
            if cpu_group is not None and cpu_group in cpu_groups:
                if cpu_set != cpu_groups[cpu_group]:
                    raise LabError("job {} cpu_group changes its explicit CPU set".format(job_id))
            elif cpu_group is not None:
                cpu_groups[cpu_group] = list(cpu_set)
        else:
            cpu_count = raw.get("cpu_count", default_cpu_count)
            _positive_int(cpu_count, "job {} cpu_count".format(job_id))
            if cpu_count > len(cpu_order):
                raise LabError("job {} requests {} CPUs, but only {} are allowed".format(
                    job_id, cpu_count, len(cpu_order)))
            if cpu_group is not None and cpu_group in cpu_groups:
                cpu_set = list(cpu_groups[cpu_group])
                if len(cpu_set) != cpu_count:
                    raise LabError(
                        "job {} cpu_group changes CPU count from {} to {}".format(
                            job_id, len(cpu_set), cpu_count))
            else:
                cpu_set = [
                    cpu_order[(cursor + offset) % len(cpu_order)]
                    for offset in range(cpu_count)]
                cpu_set = sorted(set(cpu_set))
                cursor = (cursor + cpu_count) % len(cpu_order)
                if cpu_group is not None:
                    cpu_groups[cpu_group] = list(cpu_set)

        thread_runtime = _normalize_thread_policy(
            raw.get("thread_runtime", default_thread_runtime), cpu_set,
            "job {}".format(job_id))
        environment = dict(environment)
        for key in THREAD_ENV_KEYS:
            if key not in environment:
                continue
            expected = thread_runtime["effective_env"][key]
            if environment[key].upper() != expected.upper():
                raise LabError(
                    "job {} env {}={!r} conflicts with its signed thread "
                    "runtime value {!r}".format(
                        job_id, key, environment[key], expected))
            # Keep one canonical source of truth in the signed job.  The
            # source-spec digest still records that the author supplied it.
            environment.pop(key)

        cwd = raw.get("cwd", ".")
        if not isinstance(cwd, str) or not cwd:
            raise LabError("job {} cwd must be a non-empty string".format(job_id))
        minimum_memory_mb = raw.get(
            "minimum_memory_mb", default_minimum_memory)
        _positive_int(
            minimum_memory_mb, "job {} minimum_memory_mb".format(job_id),
            allow_zero=True)
        if rss_limit_mb and minimum_memory_mb > rss_limit_mb:
            raise LabError(
                "job {} minimum_memory_mb exceeds rss_limit_mb".format(
                    job_id))
        executable = _resolve_executable(
            command, root_path, cwd, environment, executable_identities)
        performance_counters = _normalize_performance_counters(
            raw.get("performance_counters", default_performance_counters),
            root_path, cwd, environment, executable_identities)
        benchmark_cell = raw.get("benchmark_cell")
        if benchmark_cell is not None and not isinstance(benchmark_cell, dict):
            raise LabError("job {} benchmark_cell must be an object".format(job_id))
        seed = _stable_seed(selected_seed, job_id)
        job = {
            "id": job_id,
            "command": command,
            "cwd": cwd,
            "env": dict(sorted(environment.items())),
            "timeout_seconds": timeout_seconds,
            "memory_mb": memory_mb,
            "rss_limit_mb": rss_limit_mb,
            "minimum_memory_mb": minimum_memory_mb,
            "cpu_set": cpu_set,
            "thread_runtime": thread_runtime,
            "seed": seed,
            "executable": executable,
        }
        if raw.get("cpu_group") is not None:
            job["cpu_group"] = raw["cpu_group"]
        if resume_group is not None:
            job["resume_group"] = resume_group
        if benchmark_cell is not None:
            job["benchmark_cell"] = _json_copy(benchmark_cell)
        if performance_counters is not None:
            job["performance_counters"] = performance_counters
        job["job_digest"] = _digest(job)
        jobs.append(job)

    manifest = {
        "schema": MANIFEST_SCHEMA,
        "root": str(root_path),
        "base_seed": selected_seed,
        "workers": selected_workers,
        "cpu_policy": cpu_policy,
        "thread_safety": {
            "schema": THREAD_POLICY_SCHEMA,
            "controlled_environment": list(THREAD_ENV_KEYS),
            "runtime_observation": THREAD_OBSERVATION_SCHEMA,
        },
        "topology": topology,
        "source_spec": {
            "schema": spec.get("schema"),
            "digest": supplied_spec_digest or computed_spec_digest,
            "metadata": _json_copy(spec_metadata),
        },
        "jobs": jobs,
    }
    manifest["manifest_digest"] = _digest(manifest)
    validate_manifest(manifest)
    return manifest


def validate_manifest(manifest):
    if not isinstance(manifest, dict) or manifest.get("schema") != MANIFEST_SCHEMA:
        raise LabError("unsupported or missing manifest schema")
    expected_manifest_digest = manifest.get("manifest_digest")
    unsigned = dict(manifest)
    unsigned.pop("manifest_digest", None)
    if expected_manifest_digest != _digest(unsigned):
        raise LabError("manifest digest does not match its contents")
    _positive_int(manifest.get("workers"), "manifest workers")
    if manifest["workers"] > 128:
        raise LabError("manifest workers may not exceed 128")
    if not isinstance(manifest.get("root"), str):
        raise LabError("manifest root is missing")
    if manifest.get("thread_safety") != {
            "schema": THREAD_POLICY_SCHEMA,
            "controlled_environment": list(THREAD_ENV_KEYS),
            "runtime_observation": THREAD_OBSERVATION_SCHEMA}:
        raise LabError("manifest thread-safety policy is invalid")
    topology = manifest.get("topology")
    if not isinstance(topology, dict):
        raise LabError("manifest topology is missing")
    allowed_cpus = _validate_canonical_cpu_set(
        topology.get("allowed_cpus"), "manifest topology allowed_cpus")
    logical_cpus = topology.get("logical_cpus")
    cpu_entries = topology.get("cpus")
    if (topology.get("allowed_cpu_list") != format_cpu_list(allowed_cpus) or
            isinstance(logical_cpus, bool) or
            not isinstance(logical_cpus, int) or
            logical_cpus != len(allowed_cpus) or
            not isinstance(cpu_entries, list) or
            len(cpu_entries) != len(allowed_cpus)):
        raise LabError("manifest topology CPU identity is not canonical")
    entry_cpus = []
    for entry in cpu_entries:
        if (not isinstance(entry, dict) or
                any(isinstance(entry.get(name), bool) or
                    not isinstance(entry.get(name), int) or
                    entry[name] < 0
                    for name in ("cpu", "socket", "core")) or
                ("numa_node" in entry and (
                    isinstance(entry["numa_node"], bool) or
                    not isinstance(entry["numa_node"], int) or
                    entry["numa_node"] < 0))):
            raise LabError("manifest topology contains an invalid CPU entry")
        entry_cpus.append(entry["cpu"])
    if entry_cpus != allowed_cpus:
        raise LabError("manifest topology CPU entries do not match allowed_cpus")
    source_spec = manifest.get("source_spec")
    if (not isinstance(source_spec, dict) or
            not isinstance(source_spec.get("digest"), str) or
            not re.match(r"^[0-9a-f]{64}$", source_spec["digest"]) or
            not isinstance(source_spec.get("metadata"), dict)):
        raise LabError("manifest source specification identity is invalid")
    jobs = manifest.get("jobs")
    if not isinstance(jobs, list) or not jobs:
        raise LabError("manifest jobs are missing")
    ids = []
    for job in jobs:
        if not isinstance(job, dict) or not JOB_ID_RE.match(job.get("id", "")):
            raise LabError("manifest contains an invalid job")
        ids.append(job["id"])
        signed = dict(job)
        expected_job_digest = signed.pop("job_digest", None)
        if expected_job_digest != _digest(signed):
            raise LabError("job {} digest does not match its contents".format(job["id"]))
        if not isinstance(job.get("command"), list) or not job["command"]:
            raise LabError("job {} has no command".format(job["id"]))
        _validate_canonical_cpu_set(
            job.get("cpu_set"), "job {} CPU set".format(job["id"]),
            allowed_cpus)
        _validate_thread_policy(
            job.get("thread_runtime"), job["cpu_set"],
            "job {}".format(job["id"]))
        environment = job.get("env")
        if (not isinstance(environment, dict) or
                any(key in environment for key in THREAD_ENV_KEYS)):
            raise LabError(
                "job {} has noncanonical thread controls in env".format(
                    job["id"]))
        _positive_number(job.get("timeout_seconds"), "job {} timeout_seconds".format(job["id"]))
        _positive_int(job.get("memory_mb"), "job {} memory_mb".format(job["id"]), allow_zero=True)
        _positive_int(
            job.get("rss_limit_mb"),
            "job {} rss_limit_mb".format(job["id"]), allow_zero=True)
        _positive_int(
            job.get("minimum_memory_mb"),
            "job {} minimum_memory_mb".format(job["id"]), allow_zero=True)
        if (job["rss_limit_mb"] and
                job["minimum_memory_mb"] > job["rss_limit_mb"]):
            raise LabError(
                "job {} memory estimate exceeds its RSS limit".format(
                    job["id"]))
        executable = job.get("executable")
        if (not isinstance(executable, dict) or
                not isinstance(executable.get("path"), str) or
                not re.match(r"^[0-9a-f]{64}$", str(executable.get("sha256", ""))) or
                not isinstance(executable.get("size_bytes"), int) or
                executable["size_bytes"] < 0):
            raise LabError("job {} has an invalid executable identity".format(job["id"]))
        counters = job.get("performance_counters")
        if counters is not None:
            counter_executable = counters.get("executable") if isinstance(
                counters, dict) else None
            if (not isinstance(counters, dict) or
                    set(counters) != {
                        "provider", "events", "optional", "executable",
                        "probe_command"} or
                    counters.get("provider") != PERF_PROVIDER or
                    not isinstance(counters.get("events"), list) or
                    not counters["events"] or
                    not all(isinstance(event, str) and PERF_EVENT_RE.match(event)
                            for event in counters["events"]) or
                    len(counters["events"]) != len(set(counters["events"])) or
                    not isinstance(counters.get("optional"), bool) or
                    not isinstance(counter_executable, dict) or
                    not isinstance(counter_executable.get("path"), str) or
                    not re.match(r"^[0-9a-f]{64}$", str(
                        counter_executable.get("sha256", ""))) or
                    not isinstance(counter_executable.get("size_bytes"), int) or
                    counter_executable["size_bytes"] < 0 or
                    not _probe_command_matches_request(
                        counters, counters.get("probe_command"))):
                raise LabError(
                    "job {} has an invalid performance counter request".format(
                        job["id"]))
        resume_group = job.get("resume_group")
        if resume_group is not None and (
                not isinstance(resume_group, str) or not resume_group):
            raise LabError(
                "job {} has an invalid resume group".format(job["id"]))
    if ids != sorted(ids) or len(ids) != len(set(ids)):
        raise LabError("manifest jobs must have unique ids in sorted order")
    return manifest


def _job_directory(output_dir, job_id):
    return Path(output_dir) / "jobs" / job_id


def _read_completed_result(output_dir, job, rerun_failed):
    result_path = _job_directory(output_dir, job["id"]) / "result.json"
    if not result_path.is_file():
        return None
    try:
        result = _load_json(result_path)
    except LabError:
        return None
    if (result.get("schema") != RESULT_SCHEMA or result.get("state") != "complete" or
            result.get("job_digest") != job["job_digest"]):
        return None
    _validate_terminal_result(result_path, result, job)
    if rerun_failed and result.get("outcome") != "success":
        return None
    return result


def _resume_candidates(manifest, output_dir, rerun_failed):
    """Return resumable terminal results, enforcing atomic resume groups.

    Jobs without ``resume_group`` retain the historical job-granular behavior.
    A grouped job is resumable only when every manifest member of that group
    has a valid terminal result carrying the same run epoch.  Otherwise every
    member is scheduled again, preventing a timing comparison from combining
    observations made by separate runner invocations.
    """
    candidates = {
        job["id"]: _read_completed_result(output_dir, job, rerun_failed)
        for job in manifest["jobs"]
    }
    groups = {}
    for job in manifest["jobs"]:
        resume_group = job.get("resume_group")
        if resume_group is not None:
            groups.setdefault(resume_group, []).append(job)
    for members in groups.values():
        results = [candidates[job["id"]] for job in members]
        epochs = {
            result.get("run_epoch") for result in results
            if result is not None
        }
        if (any(result is None for result in results) or len(epochs) != 1 or
                not all(isinstance(epoch, str) and
                        re.match(r"^[0-9a-f]{64}$", epoch)
                        for epoch in epochs)):
            for job in members:
                candidates[job["id"]] = None
    return candidates


def _validate_thread_runtime_evidence(result, job):
    runtime = result.get("thread_runtime")
    if (not isinstance(runtime, dict) or
            set(runtime) != {"policy", "observation"} or
            runtime.get("policy") != job.get("thread_runtime")):
        raise LabError(
            "terminal thread policy is invalid for job {}".format(job["id"]))
    observation = runtime.get("observation")
    required = {
        "schema", "status", "declared_max_threads", "allocated_cpu_count",
        "sample_count", "affinity_sample_count", "peak_process_count",
        "peak_thread_count",
        "peak_rss_bytes", "rss_limit_bytes", "rss_exceeded",
        "outside_cpu_set", "oversubscribed"}
    if (not isinstance(observation, dict) or
            not required.issubset(observation) or
            not set(observation).issubset(required | {"detail"}) or
            observation.get("schema") != THREAD_OBSERVATION_SCHEMA or
            observation.get("status") not in {
                "pending", "not_run", "observed", "unavailable"} or
            observation.get("declared_max_threads") !=
                job["thread_runtime"]["max_threads"] or
            observation.get("allocated_cpu_count") != len(job["cpu_set"]) or
            any(isinstance(observation.get(key), bool) or
                not isinstance(observation.get(key), int) or
                observation[key] < 0
                for key in ("sample_count", "affinity_sample_count",
                            "peak_process_count",
                            "peak_thread_count", "peak_rss_bytes",
                            "rss_limit_bytes")) or
            not isinstance(observation.get("outside_cpu_set"), list) or
            observation["outside_cpu_set"] != sorted(set(
                observation["outside_cpu_set"])) or
            not all(isinstance(cpu, int) and not isinstance(cpu, bool) and cpu >= 0
                    for cpu in observation["outside_cpu_set"]) or
            not isinstance(observation.get("oversubscribed"), bool) or
            not isinstance(observation.get("rss_exceeded"), bool) or
            ("detail" in observation and
             (not isinstance(observation["detail"], str) or
              not observation["detail"]))):
        raise LabError(
            "terminal thread observation is invalid for job {}".format(
                job["id"]))
    expected_oversubscribed = (
        observation["peak_thread_count"] >
        observation["declared_max_threads"])
    if observation["oversubscribed"] != expected_oversubscribed:
        raise LabError(
            "terminal oversubscription evidence is inconsistent for job {}".format(
                job["id"]))
    if observation["rss_limit_bytes"] != job.get("rss_limit_mb", 0) * 1024 * 1024:
        raise LabError(
            "terminal RSS limit evidence is inconsistent for job {}".format(
                job["id"]))
    expected_rss_exceeded = bool(
        observation["rss_limit_bytes"] and
        observation["peak_rss_bytes"] > observation["rss_limit_bytes"])
    if observation["rss_exceeded"] != expected_rss_exceeded:
        raise LabError(
            "terminal RSS observation is inconsistent for job {}".format(
                job["id"]))
    if (result.get("outcome") == "success" and
            not _thread_observation_supports_success(observation)):
        raise LabError(
            "successful job {} lacks valid thread-allocation evidence".format(
                job["id"]))
    return runtime


def _validate_terminal_result(result_path, result, job):
    unsigned = dict(result)
    expected_result_digest = unsigned.pop("result_digest", None)
    if expected_result_digest != _digest(unsigned):
        raise LabError("terminal result digest is invalid for job {}".format(job["id"]))
    if result.get("job_id") != job["id"] or result.get("cpu_set") != job["cpu_set"]:
        raise LabError("terminal result identity is invalid for job {}".format(job["id"]))
    run_epoch = result.get("run_epoch")
    if run_epoch is not None and (
            not isinstance(run_epoch, str) or
            not re.match(r"^[0-9a-f]{64}$", run_epoch)):
        raise LabError("terminal result run epoch is invalid for job {}".format(job["id"]))
    if job.get("resume_group") is not None and run_epoch is None:
        raise LabError("terminal result run epoch is missing for grouped job {}".format(
            job["id"]))
    _validate_thread_runtime_evidence(result, job)
    outputs = result.get("outputs")
    if not isinstance(outputs, dict):
        raise LabError("terminal result output identities are missing for job {}".format(job["id"]))
    job_dir = Path(result_path).parent
    for name in ("stdout", "stderr"):
        expected = outputs.get(name)
        if not isinstance(expected, dict) or _content_identity(
                job_dir / (name + ".txt")) != expected:
            raise LabError("{} changed after job {} completed".format(name, job["id"]))
    _validate_performance_evidence(
        job, result.get("performance_counters"), outputs, job_dir)


def _replace_tokens(value, job):
    replacements = {
        "seed": str(job["seed"]),
        "job_id": job["id"],
        "cpu_set": format_cpu_list(job["cpu_set"]),
    }
    result = value
    for token in TOKEN_REPLACEMENTS:
        result = result.replace("{" + token + "}", replacements[token])
    return result


def _expanded_command(job):
    command = [_replace_tokens(argument, job) for argument in job["command"]]
    # Execute the exact content-addressed file recorded by the manifest.  This
    # prevents a later PATH change from silently selecting a different binary.
    command[0] = job["executable"]["path"]
    return command


def _performance_command(job, command, output_path):
    """Wrap a command in the exact content-addressed perf executable."""
    counters = job["performance_counters"]
    wrapped = [
        counters["executable"]["path"], "stat", "--no-big-num", "-x", ";",
    ]
    for event in counters["events"]:
        wrapped.extend(("-e", event))
    wrapped.extend(("-o", str(output_path), "--"))
    wrapped.extend(command)
    return wrapped


def _counter_executable_matches(job):
    counters = job.get("performance_counters")
    if counters is None:
        return True
    current = _file_identity(counters["executable"]["path"])
    expected = counters["executable"]
    return (current["sha256"] == expected["sha256"] and
            current["size_bytes"] == expected["size_bytes"])


def _validate_performance_evidence(job, counters, outputs, job_dir):
    """Re-derive counter evidence from retained raw bytes before trusting it."""
    request = job.get("performance_counters")
    has_output_identity = "performance_counters" in outputs
    output_identity = outputs.get("performance_counters")
    if request is None:
        if counters is not None or has_output_identity:
            raise LabError(
                "job {} has unexpected performance counter evidence".format(
                    job["id"]))
        return
    if (not isinstance(counters, dict) or
            not set(counters).issubset({
                "provider", "events", "optional", "executable", "probe",
                "status", "raw_output", "measurements", "detail"}) or
            counters.get("provider") != request["provider"] or
            counters.get("events") != request["events"] or
            counters.get("optional") != request["optional"] or
            counters.get("executable") != request["executable"] or
            counters.get("status") not in ("available", "partial", "unavailable")):
        raise LabError(
            "terminal performance counter evidence is invalid for job {}".format(
                job["id"]))

    probe = counters.get("probe")
    if (not isinstance(probe, dict) or
            not {"status", "command", "cpu_set", "exit_code",
                 "duration_seconds"}.issubset(probe) or
            not set(probe).issubset({
                "status", "command", "cpu_set", "exit_code",
                "duration_seconds", "stderr_tail", "detail"}) or
            probe.get("status") not in ("available", "unavailable") or
            probe.get("cpu_set") != job["cpu_set"] or
            not _probe_command_matches_request(request, probe.get("command")) or
            (probe.get("exit_code") is not None and
             (isinstance(probe.get("exit_code"), bool) or
              not isinstance(probe.get("exit_code"), int))) or
            (probe.get("status") == "available" and
             (probe.get("exit_code") != 0 or "detail" in probe)) or
            (probe.get("status") == "unavailable" and
             (not isinstance(probe.get("detail"), str) or
              not probe.get("detail"))) or
            isinstance(probe.get("duration_seconds"), bool) or
            not isinstance(probe.get("duration_seconds"), (int, float)) or
            not math.isfinite(float(probe["duration_seconds"])) or
            float(probe["duration_seconds"]) < 0.0 or
            ("stderr_tail" in probe and
             not isinstance(probe["stderr_tail"], str))):
        raise LabError(
            "terminal performance counter probe is invalid for job {}".format(
                job["id"]))

    status = counters["status"]
    raw_name = counters.get("raw_output")
    measurements = counters.get("measurements", [])
    if not isinstance(measurements, list):
        raise LabError(
            "terminal performance counter measurements are invalid for job {}".format(
                job["id"]))
    if raw_name is None:
        if has_output_identity or measurements:
            raise LabError(
                "job {} has counter measurements without retained raw output".format(
                    job["id"]))
        expected_detail = probe.get(
            "detail", "performance counters were unavailable during preflight")
        if status != "unavailable" or counters.get("detail") != expected_detail:
            raise LabError(
                "job {} counter status disagrees with its preflight".format(
                    job["id"]))
        expected_keys = {
            "provider", "events", "optional", "executable", "probe",
            "status", "detail"}
        if set(counters) != expected_keys:
            raise LabError(
                "job {} has noncanonical unavailable counter evidence".format(
                    job["id"]))
        return
    elif raw_name != "perf-stat.txt" or not isinstance(output_identity, dict):
        raise LabError(
            "performance counter output changed after job {} completed".format(
                job["id"]))

    try:
        (raw_identity, derived_measurements, derived_status,
         derived_detail) = _read_perf_stat_evidence(
             Path(job_dir) / raw_name, request["events"])
    except OSError as error:
        raise LabError(
            "cannot read retained counters for job {}: {}".format(
                job["id"], error))
    if raw_identity != output_identity:
        raise LabError(
            "performance counter output changed after job {} completed".format(
                job["id"]))

    if probe["status"] != "available":
        raise LabError(
            "job {} has raw counters from an unavailable preflight".format(
                job["id"]))

    if measurements != derived_measurements or status != derived_status:
        raise LabError(
            "job {} counter JSON differs from retained raw output".format(
                job["id"]))
    if ((derived_detail is None and "detail" in counters) or
            (derived_detail is not None and
             counters.get("detail") != derived_detail)):
        raise LabError(
            "job {} counter detail differs from retained raw output".format(
                job["id"]))
    expected_keys = {
        "provider", "events", "optional", "executable", "probe", "status",
        "raw_output", "measurements"}
    if derived_detail is not None:
        expected_keys.add("detail")
    if set(counters) != expected_keys:
        raise LabError(
            "job {} has noncanonical raw counter evidence".format(job["id"]))


def _probe_performance_counters(job):
    """Check PMU access on the exact CPU set before timing real work."""
    started = time.monotonic()
    counters = job["performance_counters"]
    logical_command = list(counters["probe_command"])
    probe_workload = logical_command[-3:]
    probe = {
        "status": "unavailable",
        "command": logical_command,
        "cpu_set": job["cpu_set"],
        "exit_code": None,
    }
    try:
        with tempfile.TemporaryDirectory(prefix="leopard2-perf-probe-") as temporary:
            output_path = Path(temporary) / "perf-stat.txt"
            command = _performance_command(
                job, probe_workload, output_path)
            completed = subprocess.run(
                command,
                env=_expanded_environment(job),
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=10.0,
                check=False,
                preexec_fn=lambda: _child_setup(job["cpu_set"], 0),
            )
            probe["exit_code"] = completed.returncode
            stderr = completed.stderr.decode("utf-8", errors="replace").strip()
            if stderr:
                probe["stderr_tail"] = stderr[-4096:]
            measurements, parse_detail = _parse_perf_stat(
                output_path, counters["events"])
            if completed.returncode == 0 and parse_detail is None:
                probe["status"] = "available"
                probe["duration_seconds"] = round(
                    time.monotonic() - started, 6)
                return probe
            details = []
            if completed.returncode != 0:
                details.append("perf probe exited {}".format(completed.returncode))
            if parse_detail:
                details.append(parse_detail)
            probe["detail"] = "; ".join(details) or "perf probe failed"
    except (OSError, subprocess.SubprocessError, LabError) as error:
        probe["detail"] = "perf probe could not run: {}".format(error)
    probe["duration_seconds"] = round(time.monotonic() - started, 6)
    return probe


def _performance_result(job, probe, output_path=None):
    request = job.get("performance_counters")
    if request is None:
        return None
    if not isinstance(probe, dict):
        probe = {
            "status": "unavailable",
            "command": [],
            "cpu_set": job["cpu_set"],
            "exit_code": None,
            "detail": "performance counter preflight was not run",
        }
    result = {
        "provider": request["provider"],
        "events": request["events"],
        "optional": request["optional"],
        "executable": request["executable"],
        "probe": probe,
        "status": "unavailable",
    }
    if output_path is None:
        result["detail"] = probe.get(
            "detail", "performance counters were unavailable during preflight")
        return result
    result["raw_output"] = "perf-stat.txt"
    measurements, status, detail = _derive_perf_stat(
        output_path, request["events"])
    result["measurements"] = measurements
    result["status"] = status
    if detail:
        result["detail"] = detail
    return result


def _expanded_environment(job, containment_token=None):
    environment = os.environ.copy()
    environment.update({key: _replace_tokens(value, job) for key, value in job["env"].items()})
    environment.update(job["thread_runtime"]["effective_env"])
    environment["LEO2_LAB_SEED"] = str(job["seed"])
    environment["LEO2_LAB_JOB_ID"] = job["id"]
    environment["LEO2_LAB_CPUSET"] = format_cpu_list(job["cpu_set"])
    if containment_token is not None:
        environment["LEO2_LAB_CONTAINMENT_TOKEN"] = containment_token
    return environment


def _new_thread_observation(job, status="pending", detail=None):
    observation = {
        "schema": THREAD_OBSERVATION_SCHEMA,
        "status": status,
        "declared_max_threads": job["thread_runtime"]["max_threads"],
        "allocated_cpu_count": len(job["cpu_set"]),
        "sample_count": 0,
        "affinity_sample_count": 0,
        "peak_process_count": 0,
        "peak_thread_count": 0,
        "peak_rss_bytes": 0,
        "rss_limit_bytes": job.get("rss_limit_mb", 0) * 1024 * 1024,
        "rss_exceeded": False,
        "outside_cpu_set": [],
        "oversubscribed": False,
    }
    if detail:
        observation["detail"] = detail
    return observation


def _thread_observation_supports_success(observation):
    """Use one predicate when classifying and later replaying a success."""
    return (
        observation.get("status") == "observed" and
        observation.get("sample_count", 0) >= 1 and
        observation.get("affinity_sample_count", 0) ==
            observation.get("sample_count", 0) and
        observation.get("peak_process_count", 0) >= 1 and
        observation.get("peak_thread_count", 0) >= 1 and
        not observation.get("oversubscribed", False) and
        not observation.get("rss_exceeded", False) and
        not observation.get("outside_cpu_set", []))


def _parse_linux_proc_stat(text):
    """Return stable Linux process identity fields, or None.

    Splitting after the final right parenthesis handles spaces and parentheses
    in ``comm``.  ``starttime_ticks`` prevents a cleanup scan from signalling a
    different process that reused a departed child's PID.
    """
    try:
        open_paren = text.index("(")
        close_paren = text.rindex(")")
        pid = int(text[:open_paren].strip())
        fields = text[close_paren + 1:].split()
        # Fields after comm begin with state, ppid, pgrp, session ... starttime.
        if len(fields) < 20:
            return None
        return {
            "pid": pid,
            "state": fields[0],
            "ppid": int(fields[1]),
            "pgrp": int(fields[2]),
            "session": int(fields[3]),
            "starttime_ticks": int(fields[19]),
        }
    except (ValueError, IndexError):
        return None


def _status_cpu_set(path):
    text = _read_text(path)
    if not text:
        return []
    for line in text.splitlines():
        if line.startswith("Cpus_allowed_list:"):
            try:
                return parse_cpu_list(line.split(":", 1)[1].strip())
            except LabError:
                return []
    return []


def _process_containment_token(entry):
    try:
        environment = (entry / "environ").read_bytes().split(b"\0")
    except OSError:
        return None
    prefix = b"LEO2_LAB_CONTAINMENT_TOKEN="
    for value in environment:
        if not value.startswith(prefix):
            continue
        try:
            token = value[len(prefix):].decode("ascii")
        except UnicodeDecodeError:
            return None
        return token if re.match(r"^[0-9a-f]{64}$", token) else None
    return None


def _stable_thread_record(task):
    """Read one TID affinity while proving its identity did not change."""
    before = _parse_linux_proc_stat(_read_text(task / "stat") or "")
    if before is None or before["state"] == "Z":
        return None
    cpu_set = _status_cpu_set(task / "status")
    after = _parse_linux_proc_stat(_read_text(task / "stat") or "")
    if (after is None or after["pid"] != before["pid"] or
            after["starttime_ticks"] != before["starttime_ticks"]):
        return None
    return {
        "tid": before["pid"],
        "starttime_ticks": before["starttime_ticks"],
        "cpu_set": cpu_set,
    }


def _linux_session_snapshot(proc_root="/proc", session_tokens=None):
    """Collect process/TID evidence keyed by each job's original session.

    ``session_tokens`` maps an original session ID to the unique environment
    token injected into that job.  Token matching keeps detached ``setsid``
    descendants inside the same runtime observation.
    """
    root = Path(proc_root)
    sessions = {}
    session_tokens = session_tokens or {}
    token_sessions = {
        token: session for session, token in session_tokens.items()}
    try:
        entries = list(root.iterdir())
    except OSError:
        return None
    for entry in entries:
        if not entry.name.isdigit():
            continue
        stat_text = _read_text(entry / "stat")
        parsed = _parse_linux_proc_stat(stat_text or "")
        if parsed is None or parsed["state"] == "Z":
            continue
        target_session = parsed["session"]
        if session_tokens and target_session not in session_tokens:
            target_session = token_sessions.get(
                _process_containment_token(entry))
            if target_session is None:
                continue
        try:
            tasks = [
                task for task in (entry / "task").iterdir()
                if task.name.isdigit()]
        except OSError:
            continue
        threads = [
            record for record in (
                _stable_thread_record(task) for task in tasks)
            if record is not None]
        rss_bytes = 0
        status_text = _read_text(entry / "status")
        if status_text:
            for line in status_text.splitlines():
                if line.startswith("VmRSS:"):
                    fields = line.split()
                    if len(fields) >= 2:
                        try:
                            rss_bytes = int(fields[1]) * 1024
                        except ValueError:
                            rss_bytes = 0
        after = _parse_linux_proc_stat(_read_text(entry / "stat") or "")
        if (after is None or after["pid"] != parsed["pid"] or
                after["starttime_ticks"] != parsed["starttime_ticks"]):
            continue
        sessions.setdefault(target_session, []).append({
            "pid": parsed["pid"],
            "starttime_ticks": parsed["starttime_ticks"],
            # Preserve every enumerated task in demand evidence.  A TID that
            # exits or becomes unreadable during the scan may lack an affinity
            # record, but it must not disappear from the observed team size.
            "thread_count": len(tasks),
            "threads": threads,
            "rss_bytes": rss_bytes,
        })
    return sessions


def _sample_thread_runtime(active_jobs, proc_root="/proc"):
    """Update all active jobs from one global /proc scan."""
    if not active_jobs:
        return
    session_tokens = {
        active["process"].pid: active["containment_token"]
        for active in active_jobs}
    sessions = _linux_session_snapshot(proc_root, session_tokens)
    for active in active_jobs:
        observation = active["thread_observation"]
        if sessions is None:
            observation["status"] = "unavailable"
            observation["detail"] = "Linux /proc runtime observation is unavailable"
            continue
        processes = list(sessions.get(active["process"].pid, []))
        if active.get("counter_active"):
            # perf stat is runner instrumentation, not part of the workload's
            # declared internal team.  Its descendant command remains counted.
            processes = [
                process for process in processes
                if process["pid"] != active["process"].pid]
        if not processes:
            continue
        process_count = len(processes)
        thread_count = sum(process["thread_count"] for process in processes)
        rss_bytes = sum(process["rss_bytes"] for process in processes)
        requested = set(active["job"]["cpu_set"])
        outside = set(observation["outside_cpu_set"])
        for process in processes:
            for thread in process["threads"]:
                outside.update(set(thread["cpu_set"]) - requested)
        observation.update(
            status="observed",
            sample_count=observation["sample_count"] + 1,
            affinity_sample_count=(
                observation["affinity_sample_count"] +
                (1 if all(
                    len(process["threads"]) == process["thread_count"] and
                    all(thread["cpu_set"] for thread in process["threads"])
                    for process in processes) else 0)),
            peak_process_count=max(
                observation["peak_process_count"], process_count),
            peak_thread_count=max(
                observation["peak_thread_count"], thread_count),
            peak_rss_bytes=max(observation["peak_rss_bytes"], rss_bytes),
            outside_cpu_set=sorted(outside),
        )
        if thread_count > observation["declared_max_threads"]:
            observation["oversubscribed"] = True
        if (observation["rss_limit_bytes"] and
                rss_bytes > observation["rss_limit_bytes"]):
            observation["rss_exceeded"] = True


def _child_setup(cpu_set, memory_mb):
    if hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, set(cpu_set))
    elif cpu_set:
        raise RuntimeError("CPU affinity is unavailable on this platform")
    if memory_mb and resource is not None:
        limit = int(memory_mb) * 1024 * 1024
        resource.setrlimit(resource.RLIMIT_AS, (limit, limit))


def _base_result(
        job, command, duration_seconds, outcome, exit_code=None, detail=None,
        run_epoch=None):
    result = {
        "schema": RESULT_SCHEMA,
        "state": "complete",
        "job_id": job["id"],
        "job_digest": job["job_digest"],
        "outcome": outcome,
        "exit_code": exit_code,
        "duration_seconds": round(max(0.0, duration_seconds), 6),
        "seed": job["seed"],
        "cpu_set": job["cpu_set"],
        "memory_mb": job["memory_mb"],
        "rss_limit_mb": job.get("rss_limit_mb", 0),
        "minimum_memory_mb": job.get("minimum_memory_mb", 0),
        "timeout_seconds": job["timeout_seconds"],
        "command": command,
        "thread_runtime": {
            "policy": _json_copy(job["thread_runtime"]),
            "observation": _new_thread_observation(
                job, status="not_run",
                detail="job did not reach runtime observation"),
        },
        "stdout": "stdout.txt",
        "stderr": "stderr.txt",
    }
    if run_epoch is not None:
        result["run_epoch"] = run_epoch
    if detail:
        result["detail"] = detail
    return result


def _write_terminal_result(job_dir, result):
    job_dir = Path(job_dir)
    result["outputs"] = {
        "stdout": _content_identity(job_dir / "stdout.txt"),
        "stderr": _content_identity(job_dir / "stderr.txt"),
    }
    counters = result.get("performance_counters")
    if isinstance(counters, dict) and counters.get("raw_output") is not None:
        raw_name = counters["raw_output"]
        if raw_name != "perf-stat.txt":
            raise LabError("invalid performance counter output name")
        result["outputs"]["performance_counters"] = _content_identity(
            job_dir / raw_name)
    unsigned = dict(result)
    unsigned.pop("result_digest", None)
    result["result_digest"] = _digest(unsigned)
    _atomic_write_json(job_dir / "result.json", result)
    return result


def _rewrite_evidence_invalid(job, output_dir, result, detail):
    """Invalidate an already written result without changing captured streams."""
    rewritten = _json_copy(result)
    rewritten["outcome"] = "evidence_invalid"
    existing_detail = rewritten.get("detail")
    if not existing_detail:
        rewritten["detail"] = detail
    elif detail not in existing_detail:
        rewritten["detail"] = existing_detail + "; " + detail
    return _write_terminal_result(
        _job_directory(output_dir, job["id"]), rewritten)


def _job_executable_matches(job):
    current = _file_identity(job["executable"]["path"])
    expected = job["executable"]
    return (current["sha256"] == expected["sha256"] and
            current["size_bytes"] == expected["size_bytes"])


def _launch_job(
        job, manifest_root, output_dir, run_epoch, performance_probe=None):
    job_dir = _job_directory(output_dir, job["id"])
    job_dir.mkdir(parents=True, exist_ok=True)
    stdout_path = job_dir / "stdout.txt"
    stderr_path = job_dir / "stderr.txt"
    counter_output_path = job_dir / "perf-stat.txt"
    try:
        counter_output_path.unlink()
    except FileNotFoundError:
        pass
    stdout_handle = stdout_path.open("wb")
    stderr_handle = stderr_path.open("wb")
    command = _expanded_command(job)
    process_command = command
    counter_active = (
        job.get("performance_counters") is not None and
        isinstance(performance_probe, dict) and
        performance_probe.get("status") == "available")
    if counter_active:
        process_command = _performance_command(
            job, command, counter_output_path)
    cwd = Path(job["cwd"])
    if not cwd.is_absolute():
        cwd = Path(manifest_root) / cwd
    started = time.monotonic()
    containment_token = hashlib.sha256(
        (run_epoch + "\0" + job["job_digest"]).encode("utf-8")).hexdigest()
    try:
        if not _job_executable_matches(job):
            raise LabError("executable changed before launch: {}".format(
                job["executable"]["path"]))
        if not _counter_executable_matches(job):
            raise LabError("performance counter executable changed before launch: {}".format(
                job["performance_counters"]["executable"]["path"]))
        process = subprocess.Popen(
            process_command,
            cwd=str(cwd),
            env=_expanded_environment(job, containment_token),
            stdin=subprocess.DEVNULL,
            stdout=stdout_handle,
            stderr=stderr_handle,
            start_new_session=True,
            preexec_fn=lambda: _child_setup(job["cpu_set"], job["memory_mb"]),
        )
    except BaseException as error:
        stdout_handle.close()
        stderr_handle.close()
        result = _base_result(
            job, command, time.monotonic() - started, "launch_error",
            detail=str(error), run_epoch=run_epoch)
        performance = _performance_result(
            job, performance_probe,
            counter_output_path if counter_output_path.is_file() else None)
        if performance is not None:
            result["performance_counters"] = performance
        _write_terminal_result(job_dir, result)
        return None, result
    active = {
        "job": job,
        "command": command,
        "process_command": process_command,
        "process": process,
        "started": started,
        "stdout_handle": stdout_handle,
        "stderr_handle": stderr_handle,
        "timed_out": False,
        "resource_limited": False,
        "thread_limited": False,
        "terminate_started": None,
        "run_epoch": run_epoch,
        "performance_probe": performance_probe,
        "counter_output_path": counter_output_path if counter_active else None,
        "counter_active": counter_active,
        "thread_observation": _new_thread_observation(job),
        "containment_token": containment_token,
    }
    # This is a real /proc observation, not a synthetic one.  Very short jobs
    # that depart before either this scan or the run-loop scan fail closed.
    _sample_thread_runtime([active])
    return active, None


def _signal_group(process, sig):
    try:
        os.killpg(process.pid, sig)
    except ProcessLookupError:
        pass
    except OSError:
        try:
            process.send_signal(sig)
        except OSError:
            pass


def _process_has_containment_token(entry, token):
    """Return whether a process retained this run/job's opaque environment tag."""
    return _process_containment_token(entry) == token


def _contained_processes(active, proc_root="/proc"):
    """Find live descendants by original session or inherited opaque token."""
    root = Path(proc_root)
    matches = []
    try:
        entries = list(root.iterdir())
    except OSError:
        return matches
    original_session = active["process"].pid
    token = active["containment_token"]
    for entry in entries:
        if not entry.name.isdigit():
            continue
        parsed = _parse_linux_proc_stat(_read_text(entry / "stat") or "")
        if (parsed is None or parsed["pid"] == os.getpid() or
                parsed["state"] == "Z"):
            continue
        if (parsed["session"] == original_session or
                _process_has_containment_token(entry, token)):
            matches.append(parsed)
    return matches


def _signal_process_identity(record, sig, proc_root="/proc"):
    """Signal only if PID and kernel start time still identify one process."""
    entry = Path(proc_root) / str(record["pid"])
    current = _parse_linux_proc_stat(_read_text(entry / "stat") or "")
    if (current is None or
            current["starttime_ticks"] != record["starttime_ticks"]):
        return
    try:
        os.kill(record["pid"], sig)
    except (ProcessLookupError, PermissionError):
        pass


def _cleanup_active_descendants(active, proc_root="/proc"):
    """Terminate descendants that outlive their signed job leader.

    A normal session scan covers ordinary children.  The opaque inherited token
    also covers children that call ``setsid()``.  Repeated PID/start-time scans
    close the common fork-during-cleanup race without signalling a reused PID.
    This is fail-closed evidence hygiene for trusted benchmark workloads, not a
    security sandbox for a hostile child that deliberately scrubs its
    environment and daemonizes.
    """
    initial = _contained_processes(active, proc_root)
    if not initial:
        return {"detected": [], "remaining": []}

    deadline = time.monotonic() + 0.5
    current = initial
    while current and time.monotonic() < deadline:
        for record in current:
            _signal_process_identity(record, signal.SIGTERM, proc_root)
        time.sleep(0.02)
        current = _contained_processes(active, proc_root)

    deadline = time.monotonic() + 0.75
    while current and time.monotonic() < deadline:
        for record in current:
            _signal_process_identity(record, signal.SIGKILL, proc_root)
        time.sleep(0.02)
        current = _contained_processes(active, proc_root)

    return {
        "detected": sorted(record["pid"] for record in initial),
        "remaining": sorted(record["pid"] for record in current),
    }


def _finish_active(active, output_dir, forced_outcome=None, detail=None):
    process = active["process"]
    try:
        exit_code = process.wait(timeout=1.0)
    except subprocess.TimeoutExpired:
        _signal_group(process, signal.SIGKILL)
        exit_code = process.wait()
    containment = _cleanup_active_descendants(active)
    active["containment_detected"] = list(containment["detected"])
    active["containment_remaining"] = list(containment["remaining"])
    active["stdout_handle"].close()
    active["stderr_handle"].close()
    executable_detail = None
    try:
        if not _job_executable_matches(active["job"]):
            executable_detail = "executable changed during execution"
        elif not _counter_executable_matches(active["job"]):
            executable_detail = "performance counter executable changed during execution"
    except LabError as error:
        executable_detail = str(error)
    if forced_outcome:
        outcome = forced_outcome
    elif active["resource_limited"]:
        outcome = "memory_limit"
    elif active["timed_out"]:
        outcome = "timeout"
    elif exit_code == 0:
        outcome = "success"
    else:
        outcome = "failed"
    if executable_detail:
        outcome = "evidence_invalid"
        detail = executable_detail if detail is None else detail + "; " + executable_detail
    if containment["detected"]:
        containment_detail = (
            "job leader exited with residual descendants {}; cleanup "
            "remaining {}".format(
                ",".join(str(pid) for pid in containment["detected"]),
                ",".join(str(pid) for pid in containment["remaining"]) or
                "none"))
        if forced_outcome is None and outcome in ("success", "failed"):
            outcome = "evidence_invalid"
        detail = (
            containment_detail if detail is None
            else detail + "; " + containment_detail)
    observation = active["thread_observation"]
    if (observation["oversubscribed"] or
            observation["outside_cpu_set"]):
        runtime_detail = (
            "runtime thread demand exceeded its signed allocation"
            if observation["oversubscribed"] else
            "runtime process affinity escaped its signed CPU set")
        if forced_outcome is None:
            outcome = "evidence_invalid"
        detail = runtime_detail if detail is None else detail + "; " + runtime_detail
    if observation["rss_exceeded"]:
        runtime_detail = "runtime RSS exceeded its signed per-job limit"
        if outcome == "success":
            outcome = "memory_limit"
        detail = runtime_detail if detail is None else detail + "; " + runtime_detail
    elif (outcome == "success" and
            not _thread_observation_supports_success(observation)):
        outcome = "evidence_invalid"
        runtime_detail = (
            "successful workload lacked a complete real runtime thread and "
            "affinity observation")
        detail = runtime_detail if detail is None else detail + "; " + runtime_detail
    result = _base_result(
        active["job"], active["command"], time.monotonic() - active["started"],
        outcome, exit_code=exit_code, detail=detail,
        run_epoch=active["run_epoch"])
    result["thread_runtime"]["observation"] = _json_copy(observation)
    counter_output_path = active.get("counter_output_path")
    performance = _performance_result(
        active["job"], active.get("performance_probe"),
        counter_output_path if (
            counter_output_path is not None and counter_output_path.is_file()) else None)
    if performance is not None:
        result["performance_counters"] = performance
    result_path = _job_directory(output_dir, active["job"]["id"]) / "result.json"
    return _write_terminal_result(result_path.parent, result)


def _can_launch(job, active, allow_cpu_overlap, memory_budget_mb=None):
    if memory_budget_mb is not None:
        active_memory_mb = sum(
            current["job"].get("minimum_memory_mb", 0) for current in active)
        if active_memory_mb + job.get("minimum_memory_mb", 0) > memory_budget_mb:
            return False
    active_cpus = set()
    aggregate_threads = job["thread_runtime"]["max_threads"]
    allocated_cpus = set(job["cpu_set"])
    for current in active:
        active_cpus.update(current["job"]["cpu_set"])
        allocated_cpus.update(current["job"]["cpu_set"])
        aggregate_threads += current["job"]["thread_runtime"]["max_threads"]
    # Even the explicit overlap mode may not oversubscribe the union of CPUs
    # assigned to simultaneously active work.  It can share partially
    # overlapping multi-CPU sets only while aggregate demand still fits.
    if aggregate_threads > len(allocated_cpus):
        return False
    return allow_cpu_overlap or not active_cpus.intersection(job["cpu_set"])


def _validate_runtime_cpus(manifest):
    current, source = _allowed_cpus()
    current_set = set(current)
    for job in manifest["jobs"]:
        unavailable = set(job["cpu_set"]) - current_set
        if unavailable:
            raise LabError(
                "job {} requests CPUs {} that are not in the current {} mask {}".format(
                    job["id"], format_cpu_list(unavailable), source, format_cpu_list(current)))


def _validate_runtime_executables(manifest):
    identities = {}
    for job in manifest["jobs"]:
        expected = job["executable"]
        path = expected["path"]
        if path not in identities:
            identities[path] = _file_identity(path)
        current = identities[path]
        if (current["sha256"] != expected["sha256"] or
                current["size_bytes"] != expected["size_bytes"]):
            raise LabError(
                "job {} executable changed after manifest creation: {}".format(
                    job["id"], expected["path"]))
        counters = job.get("performance_counters")
        if counters is None:
            continue
        expected = counters["executable"]
        path = expected["path"]
        if path not in identities:
            identities[path] = _file_identity(path)
        current = identities[path]
        if (current["sha256"] != expected["sha256"] or
                current["size_bytes"] != expected["size_bytes"]):
            raise LabError(
                "job {} performance counter executable changed after manifest "
                "creation: {}".format(job["id"], expected["path"]))


def _memory_unavailable_reason(job, memory_budget_mb):
    required_mb = job.get("minimum_memory_mb", 0)
    if not required_mb or memory_budget_mb is None:
        return None
    if required_mb > memory_budget_mb:
        return (
            "requires an estimated {} MiB but only {} MiB of the recorded "
            "and current host memory budget is available".format(
                required_mb, memory_budget_mb))
    return None


def _manifest_memory_budget_mb(manifest):
    capacities = [
        value for value in (
            manifest.get("topology", {}).get("memory_capacity_bytes"),
            _memory_capacity_bytes())
        if isinstance(value, int) and value > 0]
    if not capacities:
        return None
    capacity = min(capacities)
    return int((capacity // (1024 * 1024)) * 0.80)


def _record_unavailable(
        job, output_dir, reason, run_epoch, performance_probe=None):
    job_dir = _job_directory(output_dir, job["id"])
    job_dir.mkdir(parents=True, exist_ok=True)
    (job_dir / "stdout.txt").write_text("", encoding="utf-8")
    (job_dir / "stderr.txt").write_text(reason + "\n", encoding="utf-8")
    result = _base_result(
        job, _expanded_command(job), 0.0, "unavailable", detail=reason,
        run_epoch=run_epoch)
    if job.get("performance_counters") is not None:
        if performance_probe is None:
            performance_probe = {
                "status": "unavailable",
                "command": [],
                "cpu_set": job["cpu_set"],
                "exit_code": None,
                "detail": "job was unavailable before performance counter preflight",
            }
        result["performance_counters"] = _performance_result(
            job, performance_probe)
    return _write_terminal_result(job_dir, result)


def _dry_run_plan(manifest, output_dir, rerun_failed, workers=None):
    _validate_runtime_cpus(manifest)
    _validate_runtime_executables(manifest)
    worker_count = manifest["workers"] if workers is None else workers
    _positive_int(worker_count, "workers")
    if worker_count > 128:
        raise LabError("workers may not exceed 128")
    memory_budget_mb = _manifest_memory_budget_mb(manifest)
    candidates = _resume_candidates(manifest, output_dir, rerun_failed)
    planned = []
    for job in manifest["jobs"]:
        completed = candidates[job["id"]]
        unavailable_reason = None if completed else _memory_unavailable_reason(
            job, memory_budget_mb)
        planned.append({
            "id": job["id"],
            "action": "resume" if completed else (
                "unavailable" if unavailable_reason else "run"),
            "cpu_set": job["cpu_set"],
            "seed": job["seed"],
            "timeout_seconds": job["timeout_seconds"],
            "memory_mb": job["memory_mb"],
            "rss_limit_mb": job.get("rss_limit_mb", 0),
            "minimum_memory_mb": job.get("minimum_memory_mb", 0),
            "thread_runtime": job["thread_runtime"],
            "command": _expanded_command(job),
            "performance_counters": job.get("performance_counters"),
            "detail": unavailable_reason,
        })
    return {
        "manifest_digest": manifest["manifest_digest"],
        "workers": worker_count,
        "output_dir": str(Path(output_dir).resolve()),
        "jobs": planned,
    }


def run_manifest(manifest, output_dir, workers=None, rerun_failed=False,
                 allow_cpu_overlap=False, progress_seconds=1.0, quiet=False):
    """Execute a manifest and return a summary; result files are always durable."""
    validate_manifest(manifest)
    _validate_runtime_cpus(manifest)
    _validate_runtime_executables(manifest)
    worker_count = manifest["workers"] if workers is None else workers
    _positive_int(worker_count, "workers")
    if worker_count > 128:
        raise LabError("workers may not exceed 128")
    _positive_number(progress_seconds, "progress_seconds", allow_zero=True)
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    _atomic_write_json(output_dir / "manifest.json", manifest)
    run_epoch = _new_run_epoch()
    run_started = time.monotonic()

    results = {}
    pending = []
    resumed = 0
    preflight_unavailable = 0
    memory_budget_mb = _manifest_memory_budget_mb(manifest)
    candidates = _resume_candidates(manifest, output_dir, rerun_failed)
    for job in manifest["jobs"]:
        completed = candidates[job["id"]]
        if completed is not None:
            results[job["id"]] = completed
            resumed += 1
        else:
            unavailable_reason = _memory_unavailable_reason(job, memory_budget_mb)
            if unavailable_reason:
                results[job["id"]] = _record_unavailable(
                    job, output_dir, unavailable_reason, run_epoch)
                preflight_unavailable += 1
            else:
                pending.append(job)

    performance_probes = {}
    probe_cache = {}
    runnable = []
    for job in pending:
        request = job.get("performance_counters")
        if request is None:
            performance_probes[job["id"]] = None
            runnable.append(job)
            continue
        probe_key = _digest({
            "provider": request["provider"],
            "events": request["events"],
            "executable": request["executable"],
            "cpu_set": job["cpu_set"],
        })
        if probe_key not in probe_cache:
            probe_cache[probe_key] = _probe_performance_counters(job)
        probe = _json_copy(probe_cache[probe_key])
        performance_probes[job["id"]] = probe
        if probe["status"] == "unavailable" and not request["optional"]:
            reason = (
                "required performance counters are unavailable: " +
                probe.get("detail", "perf preflight failed"))
            results[job["id"]] = _record_unavailable(
                job, output_dir, reason, run_epoch, probe)
            preflight_unavailable += 1
        else:
            runnable.append(job)
    pending = runnable

    active = []
    launched_job_ids = set()
    jobs_by_id = {job["id"]: job for job in manifest["jobs"]}
    total = len(manifest["jobs"])
    last_progress = 0.0
    interrupted = False
    containment_unavailable = 0
    containment_failure = None
    try:
        while pending or active:
            launched_any = True
            while pending and len(active) < worker_count and launched_any:
                launched_any = False
                for index, job in enumerate(pending):
                    if not _can_launch(
                            job, active, allow_cpu_overlap, memory_budget_mb):
                        continue
                    pending.pop(index)
                    launched, launch_result = _launch_job(
                        job, manifest["root"], output_dir, run_epoch,
                        performance_probes.get(job["id"]))
                    if launched is not None:
                        active.append(launched)
                        launched_job_ids.add(job["id"])
                    else:
                        results[job["id"]] = launch_result
                    launched_any = True
                    break

            _sample_thread_runtime(active)
            now = time.monotonic()
            for current in list(active):
                process = current["process"]
                exit_code = process.poll()
                elapsed = now - current["started"]
                if current["thread_observation"]["outside_cpu_set"]:
                    containment_failure = (
                        "runner quarantine: job {} escaped its signed CPU set "
                        "onto {}".format(
                            current["job"]["id"],
                            format_cpu_list(
                                current["thread_observation"][
                                    "outside_cpu_set"])))
                    break
                if (exit_code is None and
                        current["thread_observation"]["oversubscribed"] and
                        not current["thread_limited"]):
                    current["thread_limited"] = True
                    current["terminate_started"] = now
                    _signal_group(process, signal.SIGTERM)
                elif (exit_code is None and
                        current["thread_observation"]["rss_exceeded"] and
                        not current["resource_limited"]):
                    current["resource_limited"] = True
                    current["terminate_started"] = now
                    _signal_group(process, signal.SIGTERM)
                elif (exit_code is None and not current["timed_out"] and
                      not current["resource_limited"] and
                      not current["thread_limited"] and
                      elapsed >= current["job"]["timeout_seconds"]):
                    current["timed_out"] = True
                    current["terminate_started"] = now
                    _signal_group(process, signal.SIGTERM)
                elif (exit_code is None and
                      (current["timed_out"] or current["resource_limited"] or
                       current["thread_limited"]) and
                      now - current["terminate_started"] >= 0.5):
                    _signal_group(process, signal.SIGKILL)
                if process.poll() is not None:
                    result = _finish_active(current, output_dir)
                    results[current["job"]["id"]] = result
                    active.remove(current)
                    if current.get("containment_detected"):
                        containment_failure = (
                            "runner quarantine: job {} left residual "
                            "descendants {} (remaining {})".format(
                                current["job"]["id"],
                                ",".join(
                                    str(pid) for pid in
                                    current["containment_detected"]),
                                ",".join(
                                    str(pid) for pid in
                                    current["containment_remaining"]) or
                                "none"))
                        break

            if containment_failure is not None:
                # No later job may become accepted evidence while a known
                # contaminant remains.  Stop concurrent jobs, invalidate every
                # result launched by this invocation, and terminally quarantine
                # work that had not started.
                for current in active:
                    _signal_group(current["process"], signal.SIGTERM)
                for current in list(active):
                    result = _finish_active(
                        current, output_dir,
                        forced_outcome="evidence_invalid",
                        detail=containment_failure)
                    results[current["job"]["id"]] = result
                    active.remove(current)
                for job_id in sorted(launched_job_ids):
                    if job_id not in results:
                        continue
                    results[job_id] = _rewrite_evidence_invalid(
                        jobs_by_id[job_id], output_dir, results[job_id],
                        containment_failure)
                for job in pending:
                    results[job["id"]] = _record_unavailable(
                        job, output_dir, containment_failure, run_epoch,
                        performance_probes.get(job["id"]))
                    containment_unavailable += 1
                pending = []
                break

            if not quiet and (time.monotonic() - last_progress >= progress_seconds or not pending and not active):
                counts = Counter(result["outcome"] for result in results.values())
                print(
                    "lab: {}/{} complete ({} resumed, {} active, {} pending; "
                    "{} success, {} failed, {} timeout, {} launch error, "
                    "{} unavailable)".format(
                        len(results), total, resumed, len(active), len(pending),
                        counts["success"], counts["failed"], counts["timeout"],
                        counts["launch_error"], counts["unavailable"]),
                    file=sys.stderr,
                    flush=True,
                )
                last_progress = time.monotonic()
            if active:
                time.sleep(0.05)
    except KeyboardInterrupt:
        interrupted = True
        for current in active:
            _signal_group(current["process"], signal.SIGTERM)
        deadline = time.monotonic() + 0.5
        while active and time.monotonic() < deadline:
            for current in list(active):
                if current["process"].poll() is not None:
                    result = _finish_active(current, output_dir, forced_outcome="interrupted")
                    results[current["job"]["id"]] = result
                    active.remove(current)
            time.sleep(0.02)
        for current in active:
            _signal_group(current["process"], signal.SIGKILL)
            result = _finish_active(current, output_dir, forced_outcome="interrupted")
            results[current["job"]["id"]] = result
        active = []

    merged = merge_results(manifest, output_dir, allow_missing=True)
    summary = {
        "total": total,
        "resumed": resumed,
        "executed": (
            total - resumed - preflight_unavailable -
            containment_unavailable - len(pending)),
        "preflight_unavailable": preflight_unavailable,
        "containment_unavailable": containment_unavailable,
        "pending": len(pending),
        "elapsed_seconds": round(time.monotonic() - run_started, 6),
        "outcomes": merged["summary"],
        "interrupted": interrupted,
    }
    return summary


def merge_results(manifest, output_dir, output_path=None, allow_missing=False):
    """Merge per-job JSON in manifest order into one canonical JSON document."""
    validate_manifest(manifest)
    output_dir = Path(output_dir).resolve()
    records = []
    missing = []
    for job in manifest["jobs"]:
        result_path = _job_directory(output_dir, job["id"]) / "result.json"
        if not result_path.is_file():
            missing.append(job["id"])
            continue
        result = _load_json(result_path)
        if result.get("schema") != RESULT_SCHEMA or result.get("job_digest") != job["job_digest"]:
            raise LabError("stale or invalid result for job {}".format(job["id"]))
        _validate_terminal_result(result_path, result, job)
        records.append(result)
    if missing and not allow_missing:
        raise LabError("missing results for jobs: {}".format(", ".join(missing)))
    counts = Counter(record.get("outcome", "invalid") for record in records)
    summary = {name: counts[name] for name in sorted(counts)}
    summary["missing"] = len(missing)
    merged = {
        "schema": MERGE_SCHEMA,
        "manifest_digest": manifest["manifest_digest"],
        "summary": summary,
        "missing_jobs": missing,
        "results": records,
    }
    destination = Path(output_path) if output_path else output_dir / "merged-results.json"
    _atomic_write_json(destination, merged)
    return merged


def _demo_spec(root):
    python = sys.executable
    return {
        "root": str(root),
        "base_seed": 20260716,
        "defaults": {"timeout_seconds": 5.0, "memory_mb": 256, "cpu_count": 1},
        "jobs": [
            {
                "id": "demo.hello",
                "command": [
                    python, "-c",
                    "import os,time; "
                    "print(os.environ['LEO2_LAB_JOB_ID'],"
                    "os.environ['LEO2_LAB_SEED']); time.sleep(0.08)"],
            },
            {
                "id": "demo.affinity",
                "command": [
                    python, "-c",
                    "import os,time; "
                    "print(sorted(os.sched_getaffinity(0))); "
                    "time.sleep(0.08)"],
            },
        ],
    }


def self_test():
    with tempfile.TemporaryDirectory(prefix="leopard2-lab-self-test-") as temporary:
        root = Path(temporary)
        python = sys.executable
        fake_proc = root / "fake-proc"
        (fake_proc / "self").mkdir(parents=True)
        (fake_proc / "meminfo").write_text(
            "MemTotal:       8388608 kB\n", encoding="utf-8")
        cgroup2_mount = root / "cgroup2"
        cgroup2_scope = cgroup2_mount / "slice" / "scope"
        cgroup2_scope.mkdir(parents=True)
        (cgroup2_mount / "memory.max").write_text("max\n", encoding="utf-8")
        (cgroup2_mount / "slice" / "memory.max").write_text(
            str(3 * 1024 * 1024 * 1024) + "\n", encoding="utf-8")
        (cgroup2_scope / "memory.max").write_text(
            str(4 * 1024 * 1024 * 1024) + "\n", encoding="utf-8")
        (fake_proc / "self" / "cgroup").write_text(
            "0::/slice/scope\n", encoding="utf-8")
        (fake_proc / "self" / "mountinfo").write_text(
            "29 23 0:26 / {} rw - cgroup2 cgroup rw\n".format(
                cgroup2_mount),
            encoding="utf-8")
        if _memory_capacity_bytes(fake_proc) != 3 * 1024 * 1024 * 1024:
            raise LabError(
                "self-test: cgroup v2 limiting ancestor was not discovered")

        cgroup1_mount = root / "cgroup1-memory"
        cgroup1_scope = cgroup1_mount / "job"
        cgroup1_scope.mkdir(parents=True)
        (cgroup1_mount / "memory.limit_in_bytes").write_text(
            str(5 * 1024 * 1024 * 1024) + "\n", encoding="utf-8")
        (cgroup1_scope / "memory.limit_in_bytes").write_text(
            str(7 * 1024 * 1024 * 1024) + "\n", encoding="utf-8")
        (fake_proc / "self" / "cgroup").write_text(
            "5:cpu,memory:/docker/outer/job\n", encoding="utf-8")
        (fake_proc / "self" / "mountinfo").write_text(
            "30 23 0:27 /docker/outer {} rw - cgroup cgroup rw,memory\n".format(
                cgroup1_mount),
            encoding="utf-8")
        if _memory_capacity_bytes(fake_proc) != 5 * 1024 * 1024 * 1024:
            raise LabError(
                "self-test: cgroup v1 current path or limiting ancestor was not discovered")

        spec = {
            "root": str(root),
            "workers": min(2, len(_allowed_cpus()[0])),
            "base_seed": 424242,
            "metadata": {"campaign": "lab-self-test", "revision": 2},
            "defaults": {"timeout_seconds": 2.0, "memory_mb": 0, "cpu_count": 1},
            "jobs": [
                {
                    "id": "success",
                    "command": [
                        python, "-c",
                        "import os,time; "
                        "print(os.environ['LEO2_LAB_SEED']); "
                        "print(sorted(os.sched_getaffinity(0))); "
                        "time.sleep(0.08)"],
                },
                {
                    "id": "memory-limit",
                    "command": [
                        python, "-c",
                        "import time\ntry:\n bytearray(256 * 1024 * 1024)\n"
                        "except MemoryError:\n print('limited'); time.sleep(0.08)\n"
                        "else:\n raise SystemExit(9)",
                    ],
                    "memory_mb": 96,
                },
                {"id": "failure", "command": [python, "-c", "raise SystemExit(7)"]},
                {
                    "id": "timeout",
                    "command": [python, "-c", "import time; time.sleep(2)"],
                    "timeout_seconds": 0.15,
                },
            ],
        }
        first_manifest = build_manifest(spec)
        second_manifest = build_manifest(spec)
        if _canonical_json_bytes(first_manifest) != _canonical_json_bytes(second_manifest):
            raise LabError("self-test: manifest generation is not deterministic")
        if first_manifest["source_spec"]["metadata"] != spec["metadata"]:
            raise LabError("self-test: source specification metadata was not preserved")
        if first_manifest["source_spec"]["digest"] != _digest(spec):
            raise LabError("self-test: source specification digest is incorrect")
        expected_python = _file_identity(python)
        if any(job["executable"]["sha256"] != expected_python["sha256"]
               for job in first_manifest["jobs"]):
            raise LabError("self-test: executable content was not hashed into each job")
        grouped_manifest = build_manifest({
            "root": str(root),
            "defaults": {"memory_mb": 0, "cpu_count": 1},
            "jobs": [
                {"id": "pair.automatic", "command": [python, "-c", "pass"],
                 "cpu_group": "pair"},
                {"id": "pair.forced-generic", "command": [python, "-c", "pass"],
                 "cpu_group": "pair"},
            ],
        })
        if grouped_manifest["jobs"][0]["cpu_set"] != grouped_manifest["jobs"][1]["cpu_set"]:
            raise LabError("self-test: CPU assignment group did not preserve pair affinity")

        def resign_manifest(candidate):
            for candidate_job in candidate["jobs"]:
                unsigned_job = dict(candidate_job)
                unsigned_job.pop("job_digest", None)
                candidate_job["job_digest"] = _digest(unsigned_job)
            unsigned_manifest = dict(candidate)
            unsigned_manifest.pop("manifest_digest", None)
            candidate["manifest_digest"] = _digest(unsigned_manifest)
            return candidate

        def expect_bad_resigned_manifest(label, mutation):
            candidate = _json_copy(first_manifest)
            mutation(candidate)
            resign_manifest(candidate)
            try:
                validate_manifest(candidate)
            except LabError:
                return
            raise LabError(
                "self-test: re-signed manifest accepted {}".format(label))

        first_cpu = first_manifest["topology"]["allowed_cpus"][0]
        expect_bad_resigned_manifest(
            "a duplicate CPU allocation",
            lambda candidate: candidate["jobs"][0].update(
                cpu_set=[first_cpu, first_cpu]))
        expect_bad_resigned_manifest(
            "a boolean CPU identifier",
            lambda candidate: candidate["jobs"][0].update(cpu_set=[True]))
        expect_bad_resigned_manifest(
            "a negative CPU identifier",
            lambda candidate: candidate["jobs"][0].update(cpu_set=[-1]))
        unavailable_cpu = max(
            first_manifest["topology"]["allowed_cpus"]) + 1000000
        expect_bad_resigned_manifest(
            "an unavailable CPU identifier",
            lambda candidate: candidate["jobs"][0].update(
                cpu_set=[unavailable_cpu]))
        if len(first_manifest["topology"]["allowed_cpus"]) >= 2:
            first_two = first_manifest["topology"]["allowed_cpus"][:2]
            expect_bad_resigned_manifest(
                "an unsorted CPU allocation",
                lambda candidate: candidate["jobs"][0].update(
                    cpu_set=list(reversed(first_two))))
        expect_bad_resigned_manifest(
            "a noncanonical topology CPU list",
            lambda candidate: candidate["topology"].update(
                allowed_cpus=[first_cpu, first_cpu]))
        expect_bad_resigned_manifest(
            "a boolean topology CPU count",
            lambda candidate: candidate["topology"].update(logical_cpus=True))
        expect_bad_resigned_manifest(
            "topology CPU entries outside allowed_cpus",
            lambda candidate: candidate["topology"].update(
                cpus=candidate["topology"]["cpus"] + [{
                    "cpu": unavailable_cpu, "socket": 0,
                    "core": unavailable_cpu}])),

        stat_fields = ["S"] + [str(value) for value in range(1, 20)]
        parsed_stat = _parse_linux_proc_stat(
            "123 (command ) with spaces) " + " ".join(stat_fields))
        if parsed_stat != {
                "pid": 123, "state": "S", "ppid": 1, "pgrp": 2,
                "session": 3, "starttime_ticks": 19}:
            raise LabError(
                "self-test: /proc stat identity parser lost a field")

        def expect_bad_thread_spec(fragment, label):
            candidate = {
                "root": str(root),
                "defaults": {"memory_mb": 0, "cpu_count": 1},
                "jobs": [{
                    "id": "bad-thread-spec",
                    "command": [python, "-c", "pass"],
                }],
            }
            candidate["jobs"][0].update(fragment)
            try:
                build_manifest(candidate)
            except LabError:
                return
            raise LabError("self-test: {} was accepted".format(label))

        expect_bad_thread_spec(
            {"env": {"OMP_NUM_THREADS": "999"}},
            "conflicting explicit OpenMP environment")
        expect_bad_thread_spec(
            {"thread_runtime": {"max_threads": 2}},
            "unacknowledged internal thread team")
        expect_bad_thread_spec(
            {"thread_runtime": {
                "max_threads": 2, "allow_internal_team": True}},
            "internal thread demand larger than its CPU set")

        inherited_code = (
            "import json,os,time; "
            "keys={{}}; "
            "keys.update((key,os.environ.get(key)) for key in {}); "
            "print(json.dumps(keys,sort_keys=True)); time.sleep(0.08)".format(
                repr(list(THREAD_ENV_KEYS))))
        inherited_manifest = build_manifest({
            "root": str(root),
            "workers": 1,
            "defaults": {
                "memory_mb": 0, "rss_limit_mb": 64, "cpu_count": 1},
            "jobs": [{
                "id": "inherited-thread-defaults",
                "command": [python, "-c", inherited_code],
            }],
        })
        saved_thread_environment = {
            key: os.environ.get(key) for key in THREAD_ENV_KEYS}
        try:
            for key in THREAD_COUNT_ENV:
                os.environ[key] = "999"
            for key in THREAD_FALSE_ENV:
                os.environ[key] = "TRUE"
            os.environ["OMP_MAX_ACTIVE_LEVELS"] = "999"
            inherited_output = root / "inherited-thread-results"
            inherited_summary = run_manifest(
                inherited_manifest, inherited_output, quiet=True)
        finally:
            for key, value in saved_thread_environment.items():
                if value is None:
                    os.environ.pop(key, None)
                else:
                    os.environ[key] = value
        if inherited_summary["outcomes"] != {"missing": 0, "success": 1}:
            raise LabError(
                "self-test: inherited runtime defaults escaped the cap")
        inherited_job = inherited_manifest["jobs"][0]
        inherited_result = _load_json(
            _job_directory(inherited_output, inherited_job["id"]) /
            "result.json")
        inherited_child_env = json.loads((
            _job_directory(inherited_output, inherited_job["id"]) /
            "stdout.txt").read_text(encoding="utf-8"))
        if (inherited_child_env !=
                inherited_job["thread_runtime"]["effective_env"] or
                inherited_result["thread_runtime"]["observation"][
                    "peak_thread_count"] != 1):
            raise LabError(
                "self-test: effective nested-runtime environment was not recorded")

        evidence_manifest = build_manifest({
            "root": str(root), "workers": 1,
            "defaults": {"memory_mb": 0, "cpu_count": 1},
            "jobs": [{
                "id": "missing-real-observation",
                "command": [
                    python, "-c",
                    "import threading,time; "
                    "threads=[threading.Thread(target=lambda:time.sleep(0.08)) "
                    "for _ in range(3)]; "
                    "[thread.start() for thread in threads]; "
                    "[thread.join() for thread in threads]"],
            }],
        })
        original_snapshot = globals()["_linux_session_snapshot"]
        try:
            globals()["_linux_session_snapshot"] = (
                lambda *_args, **_kwargs: {})
            evidence_output = root / "missing-observation-results"
            evidence_summary = run_manifest(
                evidence_manifest, evidence_output, quiet=True)
        finally:
            globals()["_linux_session_snapshot"] = original_snapshot
        evidence_result = _load_json(
            _job_directory(
                evidence_output, "missing-real-observation") / "result.json")
        if (evidence_summary["outcomes"] != {
                "evidence_invalid": 1, "missing": 0} or
                evidence_result["thread_runtime"]["observation"][
                    "sample_count"] != 0):
            raise LabError(
                "self-test: an unsampled success received fabricated evidence")

        original_status_cpu_set = globals()["_status_cpu_set"]
        try:
            globals()["_status_cpu_set"] = lambda _path: []
            incomplete_output = root / "incomplete-affinity-results"
            incomplete_summary = run_manifest(
                evidence_manifest, incomplete_output, quiet=True)
        finally:
            globals()["_status_cpu_set"] = original_status_cpu_set
        incomplete_result = _load_json(
            _job_directory(
                incomplete_output, "missing-real-observation") /
            "result.json")
        incomplete_observation = incomplete_result[
            "thread_runtime"]["observation"]
        if (incomplete_summary["outcomes"] != {
                "evidence_invalid": 1, "missing": 0} or
                incomplete_observation["sample_count"] < 1 or
                incomplete_observation["affinity_sample_count"] != 0):
            raise LabError(
                "self-test: incomplete real affinity evidence poisoned merge")

        unstable_tid_manifest = build_manifest({
            "root": str(root), "workers": 1,
            "defaults": {"memory_mb": 0, "cpu_count": 1},
            "jobs": [{
                "id": "unstable-tid-observation",
                "command": [
                    python, "-c",
                    "import threading,time; "
                    "threads=[threading.Thread(target=lambda:time.sleep(0.3)) "
                    "for _ in range(3)]; "
                    "[thread.start() for thread in threads]; "
                    "time.sleep(0.2); "
                    "[thread.join() for thread in threads]"],
            }],
        })
        original_stable_thread_record = globals()["_stable_thread_record"]
        try:
            globals()["_stable_thread_record"] = lambda _task: None
            unstable_tid_output = root / "unstable-tid-results"
            unstable_tid_summary = run_manifest(
                unstable_tid_manifest, unstable_tid_output, quiet=True)
        finally:
            globals()["_stable_thread_record"] = original_stable_thread_record
        unstable_tid_result = _load_json(
            _job_directory(
                unstable_tid_output, "unstable-tid-observation") /
            "result.json")
        unstable_tid_observation = unstable_tid_result[
            "thread_runtime"]["observation"]
        if (unstable_tid_summary["outcomes"] != {
                "evidence_invalid": 1, "missing": 0} or
                unstable_tid_observation["peak_thread_count"] < 4 or
                unstable_tid_observation["affinity_sample_count"] != 0):
            raise LabError(
                "self-test: unreadable TIDs disappeared from demand evidence")

        oversubscribe_code = (
            "import threading,time; "
            "gate=threading.Event(); "
            "threads=[threading.Thread(target=gate.wait) "
            "for _ in range(3)]; "
            "[thread.start() for thread in threads]; time.sleep(0.8); "
            "gate.set(); "
            "[thread.join() for thread in threads]")
        oversubscribe_manifest = build_manifest({
            "root": str(root), "workers": 1,
            "defaults": {"memory_mb": 0, "cpu_count": 1},
            "jobs": [{
                "id": "observed-oversubscription",
                "command": [python, "-c", oversubscribe_code],
            }],
        })
        oversubscribe_output = root / "oversubscribe-results"
        oversubscribe_summary = run_manifest(
            oversubscribe_manifest, oversubscribe_output, quiet=True)
        if oversubscribe_summary["outcomes"] != {
                "evidence_invalid": 1, "missing": 0}:
            raise LabError(
                "self-test: observed oversubscription was accepted: {}".format(
                    oversubscribe_summary["outcomes"]))
        oversubscribe_result = _load_json(
            _job_directory(
                oversubscribe_output, "observed-oversubscription") /
            "result.json")
        if (not oversubscribe_result["thread_runtime"]["observation"][
                "oversubscribed"] or
                oversubscribe_result["thread_runtime"]["observation"][
                    "peak_thread_count"] < 4):
            raise LabError(
                "self-test: oversubscription evidence did not record the team")

        def residual_command(marker, detached):
            child = (
                "import pathlib,time; time.sleep(0.45); "
                "pathlib.Path({!r}).write_text('escaped',encoding='utf-8')"
            ).format(str(marker))
            return (
                "import subprocess,sys; "
                "subprocess.Popen([sys.executable,'-c',{!r}],"
                "start_new_session={!r})"
            ).format(child, detached)

        for label, detached in (
                ("same-session", False), ("detached-session", True)):
            marker = root / ("residual-" + label + ".txt")
            residual_manifest = build_manifest({
                "root": str(root), "workers": 1,
                "defaults": {"memory_mb": 0, "cpu_count": 1},
                "jobs": [{
                    "id": "residual-" + label,
                    "command": [
                        python, "-c", residual_command(marker, detached)],
                }],
            })
            residual_output = root / ("residual-" + label + "-results")
            residual_summary = run_manifest(
                residual_manifest, residual_output, quiet=True)
            residual_result = _load_json(
                _job_directory(
                    residual_output, "residual-" + label) / "result.json")
            time.sleep(0.5)
            if (residual_summary["outcomes"] != {
                    "evidence_invalid": 1, "missing": 0} or
                    "residual descendants" not in
                    residual_result.get("detail", "") or marker.exists()):
                raise LabError(
                    "self-test: {} residual process escaped cleanup: {}".format(
                        label, residual_summary["outcomes"]))

        quarantine_pid = root / "quarantine-child.pid"
        quarantine_later = root / "quarantine-later-ran.txt"
        quarantine_child = (
            "import os,pathlib,time; "
            "pathlib.Path({!r}).write_text(str(os.getpid()),"
            "encoding='utf-8'); time.sleep(5)"
        ).format(str(quarantine_pid))
        quarantine_parent = (
            "import subprocess,sys; "
            "subprocess.Popen([sys.executable,'-c',{!r}],"
            "start_new_session=True)"
        ).format(quarantine_child)
        quarantine_manifest = build_manifest({
            "root": str(root), "workers": 1,
            "defaults": {"memory_mb": 0, "cpu_count": 1},
            "jobs": [
                {
                    "id": "quarantine-a-contaminant",
                    "command": [python, "-c", quarantine_parent],
                },
                {
                    "id": "quarantine-b-must-not-run",
                    "command": [
                        python, "-c",
                        "import pathlib; pathlib.Path({!r}).write_text("
                        "'ran',encoding='utf-8')".format(
                            str(quarantine_later))],
                },
            ],
        })
        original_signal_identity = globals()["_signal_process_identity"]
        try:
            globals()["_signal_process_identity"] = (
                lambda _record, _sig, _proc_root="/proc": None)
            quarantine_output = root / "quarantine-results"
            quarantine_summary = run_manifest(
                quarantine_manifest, quarantine_output, quiet=True)
        finally:
            globals()["_signal_process_identity"] = original_signal_identity
            if quarantine_pid.is_file():
                try:
                    escaped_pid = int(
                        quarantine_pid.read_text(encoding="utf-8"))
                    os.kill(escaped_pid, signal.SIGKILL)
                except (OSError, ValueError):
                    pass
        if (quarantine_summary["outcomes"] != {
                "evidence_invalid": 1, "missing": 0, "unavailable": 1} or
                quarantine_summary["containment_unavailable"] != 1 or
                quarantine_summary["executed"] != 1 or
                quarantine_later.exists()):
            raise LabError(
                "self-test: residual containment failure did not quarantine "
                "later work")

        if len(_allowed_cpus()[0]) >= 3:
            available = _allowed_cpus()[0]
            escape_cpu = available[2]
            detached_peer_marker = root / "detached-peer-finished.txt"
            detached_pending_marker = root / "detached-pending-ran.txt"
            detached_child = (
                "import os,time; "
                # Keep the escaped process alive until the runner observes and
                # quarantines it.  A short sleep makes this test depend on the
                # supervisor receiving CPU time within that window, which is
                # not guaranteed during a highly parallel ctest run.
                "os.sched_setaffinity(0,{{{}}}); time.sleep(60)"
            ).format(escape_cpu)
            detached_waiter = (
                "import subprocess,sys; "
                "child=subprocess.Popen([sys.executable,'-c',{!r}],"
                "start_new_session=True); child.wait()"
            ).format(detached_child)
            detached_manifest = build_manifest({
                "root": str(root), "workers": 2,
                # The workloads below are deliberately long-lived but the
                # self-test must still fail in bounded time if observation is
                # broken.
                "defaults": {
                    "memory_mb": 0, "cpu_count": 1,
                    "timeout_seconds": 10,
                },
                "jobs": [
                    {
                        "id": "detached-a-runtime-affinity-escape",
                        "command": [python, "-c", detached_waiter],
                        "cpu_set": available[:2],
                        "thread_runtime": {
                            "max_threads": 2,
                            "allow_internal_team": True,
                        },
                    },
                    {
                        "id": "detached-b-active-peer",
                        "command": [
                            python, "-c",
                            # This peer must remain active while the escaped
                            # process is observable.  It writes only if the
                            # quarantine fails to terminate it.
                            "import pathlib,time; time.sleep(60); "
                            "pathlib.Path({!r}).write_text("
                            "'finished',encoding='utf-8')".format(
                                str(detached_peer_marker))],
                        "cpu_set": [escape_cpu],
                    },
                    {
                        "id": "detached-c-pending",
                        "command": [
                            python, "-c",
                            "import pathlib; pathlib.Path({!r}).write_text("
                            "'ran',encoding='utf-8')".format(
                                str(detached_pending_marker))],
                        "cpu_set": [available[0]],
                    },
                ],
            })
            detached_output = root / "detached-runtime-results"
            detached_summary = run_manifest(
                detached_manifest, detached_output, quiet=True)
            detached_result = _load_json(
                _job_directory(
                    detached_output,
                    "detached-a-runtime-affinity-escape") /
                "result.json")
            if (detached_summary["outcomes"] != {
                    "evidence_invalid": 2, "missing": 0,
                    "unavailable": 1} or
                    detached_summary["containment_unavailable"] != 1 or
                    escape_cpu not in
                    detached_result["thread_runtime"]["observation"][
                        "outside_cpu_set"] or
                    detached_peer_marker.exists() or
                    detached_pending_marker.exists()):
                raise LabError(
                    "self-test: detached affinity escape did not quarantine "
                    "active and pending peers: outcomes={!r}, "
                    "containment_unavailable={!r}, outside_cpu_set={!r}, "
                    "active_peer_finished={!r}, pending_ran={!r}".format(
                        detached_summary["outcomes"],
                        detached_summary["containment_unavailable"],
                        detached_result["thread_runtime"]["observation"][
                            "outside_cpu_set"],
                        detached_peer_marker.exists(),
                        detached_pending_marker.exists()))

        rss_manifest = build_manifest({
            "root": str(root), "workers": 1,
            "defaults": {
                "memory_mb": 0, "rss_limit_mb": 32, "cpu_count": 1},
            "jobs": [{
                "id": "observed-rss-limit",
                "command": [
                    python, "-c",
                    "import time; payload=bytearray(96*1024*1024); "
                    "time.sleep(0.3)"],
            }],
        })
        rss_output = root / "rss-limit-results"
        rss_summary = run_manifest(rss_manifest, rss_output, quiet=True)
        rss_result = _load_json(
            _job_directory(rss_output, "observed-rss-limit") /
            "result.json")
        if (rss_summary["outcomes"] != {"memory_limit": 1, "missing": 0} or
                not rss_result["thread_runtime"]["observation"][
                    "rss_exceeded"] or
                rss_result["thread_runtime"]["observation"][
                    "peak_rss_bytes"] <= 32 * 1024 * 1024):
            raise LabError(
                "self-test: sampled RSS limit did not bound the job")

        if len(_allowed_cpus()[0]) >= 2:
            intentional_code = (
                "import threading,time; "
                "thread=threading.Thread(target=lambda:time.sleep(0.2)); "
                "thread.start(); time.sleep(0.12); thread.join()")
            intentional_manifest = build_manifest({
                "root": str(root), "workers": 1,
                "defaults": {"memory_mb": 0, "cpu_count": 2},
                "jobs": [{
                    "id": "intentional-team",
                    "command": [python, "-c", intentional_code],
                    "thread_runtime": {
                        "max_threads": 2,
                        "allow_internal_team": True,
                    },
                }],
            })
            intentional_output = root / "intentional-team-results"
            intentional_summary = run_manifest(
                intentional_manifest, intentional_output, quiet=True)
            intentional_result = _load_json(
                _job_directory(intentional_output, "intentional-team") /
                "result.json")
            if (intentional_summary["outcomes"] != {
                    "missing": 0, "success": 1} or
                    intentional_result["thread_runtime"]["observation"][
                        "peak_thread_count"] != 2):
                raise LabError(
                    "self-test: explicit internal thread team was not honored")

        if len(_allowed_cpus()[0]) >= 3:
            available = _allowed_cpus()[0]
            escape_cpu = available[2]
            affinity_escape_code = (
                "import os,threading,time\n"
                "ready=threading.Event()\n"
                "def worker():\n"
                " os.sched_setaffinity(threading.get_native_id(),{{{}}})\n"
                " ready.set()\n"
                " time.sleep(0.4)\n"
                "thread=threading.Thread(target=worker)\n"
                "thread.start()\n"
                "ready.wait()\n"
                "time.sleep(0.3)\n"
                "thread.join()\n"
            ).format(escape_cpu)
            affinity_escape_manifest = build_manifest({
                "root": str(root), "workers": 1,
                "defaults": {"memory_mb": 0, "cpu_count": 1},
                "jobs": [{
                    "id": "per-thread-affinity-escape",
                    "command": [python, "-c", affinity_escape_code],
                    "cpu_set": available[:2],
                    "thread_runtime": {
                        "max_threads": 2,
                        "allow_internal_team": True,
                    },
                }],
            })
            affinity_escape_output = root / "affinity-escape-results"
            affinity_escape_summary = run_manifest(
                affinity_escape_manifest, affinity_escape_output, quiet=True)
            affinity_escape_result = _load_json(
                _job_directory(
                    affinity_escape_output, "per-thread-affinity-escape") /
                "result.json")
            if (affinity_escape_summary["outcomes"] != {
                    "evidence_invalid": 1, "missing": 0} or
                    escape_cpu not in
                    affinity_escape_result["thread_runtime"]["observation"][
                        "outside_cpu_set"]):
                raise LabError(
                    "self-test: a non-leader thread affinity escape was "
                    "accepted")

        fake_perf = root / "fake-perf.py"
        fake_perf.write_text(
            "#!/usr/bin/env python3\n"
            "import subprocess, sys\n"
            "args = sys.argv[1:]\n"
            "if not args or args[0] != 'stat' or '--' not in args:\n"
            "    raise SystemExit(64)\n"
            "events = [args[index + 1] for index, value in enumerate(args[:-1]) "
            "if value == '-e']\n"
            "output = args[args.index('-o') + 1]\n"
            "command = args[args.index('--') + 1:]\n"
            "completed = subprocess.run(command, check=False)\n"
            "with open(output, 'w', encoding='utf-8') as destination:\n"
            "    for index, event in enumerate(events):\n"
            "        destination.write('{};;{};1.000;100.00;\\n'.format("
            "1000 + index, event))\n"
            "raise SystemExit(completed.returncode)\n",
            encoding="utf-8")
        fake_perf.chmod(0o700)
        counter_request = {
            "provider": PERF_PROVIDER,
            "command": str(fake_perf),
            "events": ["cycles", "instructions"],
            "optional": True,
        }
        counter_manifest = build_manifest({
            "root": str(root),
            "workers": 1,
            "defaults": {
                "memory_mb": 0,
                "cpu_count": 1,
                "performance_counters": counter_request,
            },
            "jobs": [{
                "id": "counter-success",
                "command": [
                    python, "-c",
                    # The counter wrapper is instrumentation and is excluded
                    # from workload-team accounting.  Keep its child alive
                    # long enough for a loaded supervisor to observe the real
                    # workload instead of making success depend on an 80-ms
                    # scheduling window.
                    "import time; print('counter child'); time.sleep(2)"],
            }],
        })
        counter_job = counter_manifest["jobs"][0]
        if (counter_job["performance_counters"]["executable"] !=
                _file_identity(fake_perf)):
            raise LabError(
                "self-test: performance counter executable was not content-addressed")
        counter_output = root / "counter-results"
        counter_summary = run_manifest(
            counter_manifest, counter_output, quiet=True)
        counter_result_path = (
            _job_directory(counter_output, counter_job["id"]) / "result.json")
        counter_result = _load_json(counter_result_path)
        if counter_summary["outcomes"] != {"missing": 0, "success": 1}:
            raise LabError(
                "self-test: fake performance counters did not run: "
                "outcomes={!r}, outcome={!r}, detail={!r}, observation={!r}, "
                "counters={!r}".format(
                    counter_summary["outcomes"],
                    counter_result.get("outcome"),
                    counter_result.get("detail"),
                    counter_result.get("thread_runtime", {}).get(
                        "observation"),
                    counter_result.get("performance_counters")))
        counter_evidence = counter_result.get("performance_counters", {})
        if (counter_evidence.get("status") != "available" or
                [entry.get("event") for entry in
                 counter_evidence.get("measurements", [])] !=
                counter_request["events"] or
                "performance_counters" not in counter_result.get("outputs", {})):
            raise LabError(
                "self-test: counted events were not retained as signed evidence")
        counter_second = run_manifest(
            counter_manifest, counter_output, quiet=True)
        if counter_second["resumed"] != 1 or counter_second["executed"] != 0:
            raise LabError(
                "self-test: performance counter result was not resumable")

        def expect_invalid_counter_result(label, mutation):
            candidate = _json_copy(counter_result)
            candidate.pop("result_digest", None)
            mutation(candidate)
            candidate["result_digest"] = _digest(candidate)
            original_bytes = counter_result_path.read_bytes()
            _atomic_write_json(counter_result_path, candidate)
            try:
                _dry_run_plan(counter_manifest, counter_output, False)
            except LabError:
                return
            finally:
                counter_result_path.write_bytes(original_bytes)
            raise LabError(
                "self-test: invalid counter evidence was resumed: " + label)

        def remove_counter_raw_identity(candidate):
            candidate["performance_counters"].pop("raw_output")
            candidate["outputs"].pop("performance_counters")

        expect_invalid_counter_result(
            "available measurements without retained raw output",
            remove_counter_raw_identity)
        expect_invalid_counter_result(
            "available result with unavailable probe",
            lambda candidate: candidate["performance_counters"]["probe"].update(
                status="unavailable"))
        expect_invalid_counter_result(
            "reported event differs from requested event",
            lambda candidate: candidate["performance_counters"]["measurements"][
                0].update(reported_event="instructions"))
        expect_invalid_counter_result(
            "JSON counter value differs from retained raw output",
            lambda candidate: candidate["performance_counters"]["measurements"][
                0].update(raw_value="999999999999", value=999999999999.0))
        expect_invalid_counter_result(
            "probe argv and successful exit were forged",
            lambda candidate: candidate["performance_counters"]["probe"].update(
                command=["bogus-probe"], exit_code=99))

        with (_job_directory(counter_output, counter_job["id"]) /
              "perf-stat.txt").open("a", encoding="utf-8") as output:
            output.write("corruption\n")
        try:
            _dry_run_plan(counter_manifest, counter_output, False)
        except LabError:
            pass
        else:
            raise LabError(
                "self-test: changed performance counter evidence was accepted")

        denied_perf = root / "denied-perf.py"
        denied_perf.write_text(
            "#!/usr/bin/env python3\n"
            "import sys\n"
            "sys.stderr.write('permission denied by test policy\\n')\n"
            "raise SystemExit(255)\n",
            encoding="utf-8")
        denied_perf.chmod(0o700)

        def denied_counter_manifest(optional):
            return build_manifest({
                "root": str(root),
                "workers": 1,
                "defaults": {
                    "memory_mb": 0,
                    "cpu_count": 1,
                    "performance_counters": {
                        "command": str(denied_perf),
                        "events": ["cycles"],
                        "optional": optional,
                    },
                },
                "jobs": [{
                    "id": "counter-denied-{}".format(
                        "optional" if optional else "required"),
                    "command": [
                        python, "-c",
                        "import time; print('ran without counters'); "
                        "time.sleep(0.08)"],
                }],
            })

        optional_manifest = denied_counter_manifest(True)
        optional_output = root / "counter-denied-optional-results"
        optional_summary = run_manifest(
            optional_manifest, optional_output, quiet=True)
        optional_result = _load_json(
            _job_directory(optional_output, optional_manifest["jobs"][0]["id"]) /
            "result.json")
        if (optional_summary["outcomes"] != {"missing": 0, "success": 1} or
                optional_result["performance_counters"]["status"] !=
                "unavailable" or
                "performance_counters" in optional_result["outputs"]):
            raise LabError(
                "self-test: optional denied counters did not preserve bare execution")

        required_manifest = denied_counter_manifest(False)
        required_output = root / "counter-denied-required-results"
        required_summary = run_manifest(
            required_manifest, required_output, quiet=True)
        if (required_summary["outcomes"] != {
                "missing": 0, "unavailable": 1} or
                required_summary["executed"] != 0 or
                required_summary["preflight_unavailable"] != 1):
            raise LabError(
                "self-test: required denied counters were not preflight unavailable")

        wrong_event_perf = root / "wrong-event-perf.py"
        wrong_event_perf.write_text(
            "#!/usr/bin/env python3\n"
            "import subprocess, sys\n"
            "args = sys.argv[1:]\n"
            "output = args[args.index('-o') + 1]\n"
            "command = args[args.index('--') + 1:]\n"
            "completed = subprocess.run(command, check=False)\n"
            "with open(output, 'w', encoding='utf-8') as destination:\n"
            "    destination.write('123;;instructions;1.000;100.00;\\n')\n"
            "raise SystemExit(completed.returncode)\n",
            encoding="utf-8")
        wrong_event_perf.chmod(0o700)
        wrong_event_manifest = build_manifest({
            "root": str(root),
            "workers": 1,
            "defaults": {
                "memory_mb": 0,
                "cpu_count": 1,
                "performance_counters": {
                    "command": str(wrong_event_perf),
                    "events": ["cycles"],
                    "optional": False,
                },
            },
            "jobs": [{
                "id": "counter-wrong-reported-event",
                "command": [python, "-c", "print('must not execute')"],
            }],
        })
        wrong_event_output = root / "counter-wrong-event-results"
        wrong_event_summary = run_manifest(
            wrong_event_manifest, wrong_event_output, quiet=True)
        wrong_event_result = _load_json(
            _job_directory(
                wrong_event_output, wrong_event_manifest["jobs"][0]["id"]) /
            "result.json")
        if (wrong_event_summary["outcomes"] != {
                "missing": 0, "unavailable": 1} or
                wrong_event_summary["executed"] != 0 or
                wrong_event_result["performance_counters"]["probe"][
                    "status"] != "unavailable" or
                "performance_counters" in wrong_event_result["outputs"]):
            raise LabError(
                "self-test: positional event relabeling made wrong PMU data available")

        mutable_counter_manifest = build_manifest({
            "root": str(root),
            "defaults": {
                "memory_mb": 0,
                "performance_counters": {
                    "command": str(denied_perf),
                    "events": ["cycles"],
                },
            },
            "jobs": [{"id": "mutable-counter", "command": [python, "-c", "pass"]}],
        })
        denied_perf.write_text(
            "#!/usr/bin/env python3\nraise SystemExit(254)\n",
            encoding="utf-8")
        try:
            _dry_run_plan(
                mutable_counter_manifest, root / "mutable-counter-results", False)
        except LabError:
            pass
        else:
            raise LabError(
                "self-test: changed performance counter executable was accepted")

        atomic_manifest = build_manifest({
            "root": str(root),
            "workers": 1,
            "defaults": {"memory_mb": 0, "cpu_count": 1},
            "jobs": [
                {"id": "atomic-a", "command": [
                    python, "-c",
                    "import time; print('a'); time.sleep(0.08)"],
                 "cpu_group": "atomic", "resume_group": "atomic"},
                {"id": "atomic-b", "command": [
                    python, "-c",
                    "import time; print('b'); time.sleep(0.08)"],
                 "cpu_group": "atomic", "resume_group": "atomic"},
            ],
        })
        if any(job.get("resume_group") != "atomic"
               for job in atomic_manifest["jobs"]):
            raise LabError("self-test: resume group was not signed into the manifest")
        atomic_output = root / "atomic-results"
        atomic_first = run_manifest(atomic_manifest, atomic_output, quiet=True)
        if atomic_first["executed"] != 2 or atomic_first["resumed"] != 0:
            raise LabError("self-test: atomic group did not execute initially")

        def atomic_results():
            return {
                job["id"]: _load_json(
                    _job_directory(atomic_output, job["id"]) / "result.json")
                for job in atomic_manifest["jobs"]
            }

        first_atomic_results = atomic_results()
        first_epochs = {
            result.get("run_epoch") for result in first_atomic_results.values()}
        if (len(first_epochs) != 1 or
                not all(re.match(r"^[0-9a-f]{64}$", epoch or "")
                        for epoch in first_epochs)):
            raise LabError("self-test: atomic group did not share a valid run epoch")
        first_atomic_bytes = {
            job_id: (_job_directory(atomic_output, job_id) / "result.json").read_bytes()
            for job_id in first_atomic_results}
        atomic_second = run_manifest(atomic_manifest, atomic_output, quiet=True)
        if atomic_second["resumed"] != 2 or atomic_second["executed"] != 0:
            raise LabError("self-test: homogeneous atomic group was not resumed")
        if first_atomic_bytes != {
                job_id: (_job_directory(atomic_output, job_id) / "result.json").read_bytes()
                for job_id in first_atomic_results}:
            raise LabError("self-test: atomic resume rewrote homogeneous results")

        missing_atomic_job = atomic_manifest["jobs"][0]["id"]
        shutil.rmtree(str(_job_directory(atomic_output, missing_atomic_job)))
        atomic_plan = _dry_run_plan(
            atomic_manifest, atomic_output, rerun_failed=False)
        if any(job["action"] != "run" for job in atomic_plan["jobs"]):
            raise LabError("self-test: partial atomic group was not wholly rescheduled")
        atomic_third = run_manifest(atomic_manifest, atomic_output, quiet=True)
        if atomic_third["resumed"] != 0 or atomic_third["executed"] != 2:
            raise LabError("self-test: partial atomic group did not rerun every member")
        third_atomic_results = atomic_results()
        third_epochs = {
            result.get("run_epoch") for result in third_atomic_results.values()}
        if len(third_epochs) != 1 or third_epochs == first_epochs:
            raise LabError("self-test: atomic rerun did not establish a new shared epoch")

        mixed_job = atomic_manifest["jobs"][0]
        mixed_path = _job_directory(atomic_output, mixed_job["id"]) / "result.json"
        mixed_result = _load_json(mixed_path)
        mixed_result["run_epoch"] = "f" * 64
        mixed_unsigned = dict(mixed_result)
        mixed_unsigned.pop("result_digest", None)
        mixed_result["result_digest"] = _digest(mixed_unsigned)
        _atomic_write_json(mixed_path, mixed_result)
        mixed_plan = _dry_run_plan(
            atomic_manifest, atomic_output, rerun_failed=False)
        if any(job["action"] != "run" for job in mixed_plan["jobs"]):
            raise LabError("self-test: mixed-epoch atomic group was not wholly rescheduled")
        atomic_fourth = run_manifest(atomic_manifest, atomic_output, quiet=True)
        if atomic_fourth["resumed"] != 0 or atomic_fourth["executed"] != 2:
            raise LabError("self-test: mixed-epoch group did not rerun every member")
        if len({result.get("run_epoch")
                for result in atomic_results().values()}) != 1:
            raise LabError("self-test: mixed-epoch rerun did not restore one epoch")

        granular_manifest = build_manifest({
            "root": str(root),
            "workers": 1,
            "defaults": {"memory_mb": 0, "cpu_count": 1},
            "jobs": [
                {"id": "granular-a", "command": [
                    python, "-c",
                    "import time; print('a'); time.sleep(0.08)"]},
                {"id": "granular-b", "command": [
                    python, "-c",
                    "import time; print('b'); time.sleep(0.08)"]},
            ],
        })
        granular_output = root / "granular-results"
        run_manifest(granular_manifest, granular_output, quiet=True)
        shutil.rmtree(str(_job_directory(granular_output, "granular-a")))
        granular_second = run_manifest(
            granular_manifest, granular_output, quiet=True)
        if granular_second["resumed"] != 1 or granular_second["executed"] != 1:
            raise LabError("self-test: ungrouped resume stopped being job-granular")

        output_dir = root / "results"
        first_summary = run_manifest(first_manifest, output_dir, progress_seconds=0.05, quiet=True)
        expected = {"success": 2, "failed": 1, "timeout": 1, "missing": 0}
        if first_summary["outcomes"] != expected:
            raise LabError("self-test: unexpected outcomes {}".format(first_summary["outcomes"]))
        success_job = next(job for job in first_manifest["jobs"] if job["id"] == "success")
        stdout = (_job_directory(output_dir, "success") / "stdout.txt").read_text(encoding="utf-8")
        if str(success_job["seed"]) not in stdout:
            raise LabError("self-test: child did not receive its stable seed")
        result_files_before = {
            job["id"]: (_job_directory(output_dir, job["id"]) / "result.json").read_bytes()
            for job in first_manifest["jobs"]
        }
        second_summary = run_manifest(first_manifest, output_dir, progress_seconds=0.05, quiet=True)
        if second_summary["resumed"] != 4:
            raise LabError("self-test: completed jobs were not resumed")
        result_files_after = {
            job["id"]: (_job_directory(output_dir, job["id"]) / "result.json").read_bytes()
            for job in first_manifest["jobs"]
        }
        if result_files_before != result_files_after:
            raise LabError("self-test: resume rewrote completed result files")
        first_merge = output_dir / "merge-a.json"
        second_merge = output_dir / "merge-b.json"
        merge_results(first_manifest, output_dir, first_merge)
        merge_results(first_manifest, output_dir, second_merge)
        if first_merge.read_bytes() != second_merge.read_bytes():
            raise LabError("self-test: merge output is not deterministic")
        plan = _dry_run_plan(first_manifest, output_dir, rerun_failed=False)
        if any(job["action"] != "resume" for job in plan["jobs"]):
            raise LabError("self-test: dry-run did not identify resumable jobs")
        corrupted_stdout = _job_directory(output_dir, "success") / "stdout.txt"
        with corrupted_stdout.open("a", encoding="utf-8") as output:
            output.write("corruption\n")
        try:
            _dry_run_plan(first_manifest, output_dir, rerun_failed=False)
        except LabError:
            pass
        else:
            raise LabError("self-test: post-run output corruption was accepted")

        stale_manifest = _json_copy(first_manifest)
        stale_manifest["source_spec"]["metadata"]["revision"] = 3
        try:
            validate_manifest(stale_manifest)
        except LabError:
            pass
        else:
            raise LabError("self-test: stale manifest digest was accepted")

        signed_spec = _json_copy(spec)
        signed_spec["spec_digest"] = _digest(signed_spec)
        signed_spec["metadata"]["revision"] = 3
        try:
            build_manifest(signed_spec)
        except LabError:
            pass
        else:
            raise LabError("self-test: stale source-spec digest was accepted")

        mutable_executable = root / "mutable-executable.py"
        mutable_executable.write_text("#!/usr/bin/env python3\nprint('one')\n", encoding="utf-8")
        mutable_executable.chmod(0o700)
        executable_manifest = build_manifest({
            "root": str(root),
            "defaults": {"memory_mb": 0},
            "jobs": [{"id": "mutable", "command": [str(mutable_executable)]}],
        })
        mutable_executable.write_text("#!/usr/bin/env python3\nprint('two')\n", encoding="utf-8")
        try:
            _dry_run_plan(executable_manifest, root / "mutable-results", False)
        except LabError:
            pass
        else:
            raise LabError("self-test: changed executable content was accepted")

        queued_target = root / "queued-target.py"
        queued_target.write_text("#!/usr/bin/env python3\nprint('original')\n", encoding="utf-8")
        queued_target.chmod(0o700)
        mutate_code = (
            "import time; from pathlib import Path; "
            "Path({!r}).write_text({!r}); time.sleep(0.08)").format(
            str(queued_target), "#!/usr/bin/env python3\nprint('changed')\n")
        queued_manifest = build_manifest({
            "root": str(root), "workers": 1,
            "defaults": {"memory_mb": 0, "cpu_count": 1},
            "jobs": [
                {"id": "queued-a-mutate", "command": [python, "-c", mutate_code]},
                {"id": "queued-b-target", "command": [str(queued_target)]},
            ],
        })
        queued_summary = run_manifest(
            queued_manifest, root / "queued-results", quiet=True)
        if queued_summary["outcomes"] != {
                "launch_error": 1, "missing": 0, "success": 1}:
            raise LabError(
                "self-test: delayed executable replacement was not rejected: {}".format(
                    queued_summary["outcomes"]))

        if len(_allowed_cpus()[0]) >= 2:
            available = _allowed_cpus()[0]
            try:
                build_manifest({
                    "root": str(root),
                    "jobs": [
                        {"id": "group-a", "command": [python, "-c", "pass"],
                         "cpu_group": "explicit", "cpu_set": [available[0]]},
                        {"id": "group-b", "command": [python, "-c", "pass"],
                         "cpu_group": "explicit", "cpu_set": [available[1]]},
                    ],
                })
            except LabError:
                pass
            else:
                raise LabError("self-test: explicit cpu_group mismatch was accepted")

        total_mb = max(1, (_memory_total_bytes() or (1024 * 1024)) // (1024 * 1024))
        unavailable_manifest = build_manifest({
            "root": str(root),
            "defaults": {"memory_mb": 0},
            "jobs": [{
                "id": "unavailable-memory",
                "command": [python, "-c", "raise SystemExit(99)"],
                "minimum_memory_mb": total_mb + 1,
            }],
        })
        unavailable_output = root / "unavailable-results"
        unavailable_summary = run_manifest(
            unavailable_manifest, unavailable_output, quiet=True)
        if unavailable_summary["outcomes"] != {"missing": 0, "unavailable": 1}:
            raise LabError(
                "self-test: memory preflight did not record unavailable: {}".format(
                    unavailable_summary["outcomes"]))
        if (unavailable_summary["executed"] != 0 or
                unavailable_summary["preflight_unavailable"] != 1):
            raise LabError(
                "self-test: unavailable work was counted as executed: {}".format(
                    unavailable_summary))

        capacity_bytes = unavailable_manifest["topology"].get("memory_capacity_bytes")
        if len(_allowed_cpus()[0]) >= 2 and isinstance(capacity_bytes, int):
            capacity_mb = max(4, capacity_bytes // (1024 * 1024))
            budget_mb = int(capacity_mb * 0.80)
            per_job_mb = max(1, (budget_mb * 3) // 5)
            aggregate_manifest = build_manifest({
                "root": str(root),
                "workers": 2,
                "defaults": {"memory_mb": 0, "cpu_count": 1},
                "jobs": [
                    {"id": "aggregate-a", "command": [
                        python, "-c", "import time; time.sleep(0.08)"],
                     "minimum_memory_mb": per_job_mb},
                    {"id": "aggregate-b", "command": [
                        python, "-c", "import time; time.sleep(0.08)"],
                     "minimum_memory_mb": per_job_mb},
                ],
            })
            aggregate_summary = run_manifest(
                aggregate_manifest, root / "aggregate-results", quiet=True)
            if (aggregate_summary["outcomes"] != {"missing": 0, "success": 2} or
                    aggregate_summary["executed"] != 2):
                raise LabError(
                    "self-test: aggregate-memory scheduling failed: {}".format(
                        aggregate_summary))
    print("leopard2_lab self-test: PASS")


def _build_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(dest="command", required=True)

    topology_parser = subparsers.add_parser("topology", help="print process-visible CPU/NUMA topology")
    topology_parser.add_argument("--compact", action="store_true", help="emit one-line JSON")

    manifest_parser = subparsers.add_parser("manifest", help="create a deterministic manifest from a JSON job spec")
    manifest_parser.add_argument("--spec", required=True, help="input job-spec JSON")
    manifest_parser.add_argument("--output", required=True, help="output manifest JSON")
    manifest_parser.add_argument("--root", help="override the manifest working root")
    manifest_parser.add_argument("--workers", type=int, help="maximum concurrent jobs (default min(allowed CPUs, 128))")
    manifest_parser.add_argument("--base-seed", type=int, help="override the stable seed namespace")

    run_parser = subparsers.add_parser("run", help="run or resume a manifest")
    run_parser.add_argument("--manifest", required=True)
    run_parser.add_argument("--output-dir", required=True)
    run_parser.add_argument("--workers", type=int, help="override manifest concurrency, maximum 128")
    run_parser.add_argument("--rerun-failed", action="store_true", help="rerun non-success terminal results")
    run_parser.add_argument("--allow-cpu-overlap", action="store_true", help="allow simultaneous jobs to share CPUs")
    run_parser.add_argument("--progress-seconds", type=float, default=1.0)
    run_parser.add_argument("--dry-run", action="store_true", help="validate and print the schedule without running jobs")
    run_parser.add_argument("--quiet", action="store_true", help="suppress live progress")

    merge_parser = subparsers.add_parser("merge", help="merge per-job result JSON deterministically")
    merge_parser.add_argument("--manifest", required=True)
    merge_parser.add_argument("--output-dir", required=True)
    merge_parser.add_argument("--output", help="output JSON (default OUTPUT_DIR/merged-results.json)")
    merge_parser.add_argument("--allow-missing", action="store_true")

    demo_parser = subparsers.add_parser("demo", help="write a small runnable demonstration manifest")
    demo_parser.add_argument("--output-dir", required=True)
    demo_parser.add_argument("--run", action="store_true", help="execute the demonstration")
    demo_parser.add_argument("--dry-run", action="store_true", help="print the demonstration schedule")

    subparsers.add_parser("self-test", help="exercise manifests, pinning, timeout, resume, and merge")
    return parser


def main(argv=None):
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "topology":
            topology = detect_topology()
            if args.compact:
                print(json.dumps(topology, sort_keys=True, separators=(",", ":")))
            else:
                print(json.dumps(topology, sort_keys=True, indent=2))
            return 0
        if args.command == "manifest":
            spec = _load_json(args.spec)
            manifest = build_manifest(spec, root=args.root, workers=args.workers, base_seed=args.base_seed)
            _atomic_write_json(args.output, manifest)
            print("wrote {} jobs to {} (digest {})".format(
                len(manifest["jobs"]), args.output, manifest["manifest_digest"]))
            return 0
        if args.command == "run":
            manifest = validate_manifest(_load_json(args.manifest))
            if args.dry_run:
                print(json.dumps(
                    _dry_run_plan(manifest, args.output_dir, args.rerun_failed, workers=args.workers),
                    indent=2, sort_keys=True))
                return 0
            summary = run_manifest(
                manifest, args.output_dir, workers=args.workers, rerun_failed=args.rerun_failed,
                allow_cpu_overlap=args.allow_cpu_overlap, progress_seconds=args.progress_seconds,
                quiet=args.quiet)
            print(json.dumps(summary, indent=2, sort_keys=True))
            outcomes = summary["outcomes"]
            failures = sum(count for name, count in outcomes.items() if name not in ("success", "missing"))
            return 130 if summary["interrupted"] else (1 if failures or outcomes.get("missing", 0) else 0)
        if args.command == "merge":
            manifest = validate_manifest(_load_json(args.manifest))
            merged = merge_results(manifest, args.output_dir, args.output, allow_missing=args.allow_missing)
            print(json.dumps(merged["summary"], indent=2, sort_keys=True))
            return 0
        if args.command == "demo":
            output_dir = Path(args.output_dir).resolve()
            output_dir.mkdir(parents=True, exist_ok=True)
            manifest = build_manifest(_demo_spec(os.getcwd()), workers=min(2, len(_allowed_cpus()[0])))
            manifest_path = output_dir / "demo-manifest.json"
            _atomic_write_json(manifest_path, manifest)
            print("wrote demonstration manifest to {}".format(manifest_path))
            if args.dry_run:
                print(json.dumps(_dry_run_plan(manifest, output_dir / "results", False), indent=2, sort_keys=True))
            if args.run:
                summary = run_manifest(manifest, output_dir / "results")
                print(json.dumps(summary, indent=2, sort_keys=True))
                return 1 if any(name not in ("success", "missing") and count
                                for name, count in summary["outcomes"].items()) else 0
            return 0
        if args.command == "self-test":
            self_test()
            return 0
        parser.error("unknown command")
    except LabError as error:
        print("leopard2_lab: error: {}".format(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
