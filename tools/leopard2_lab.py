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

try:
    import resource
except ImportError:  # pragma: no cover - Unix/Linux is the production target.
    resource = None


MANIFEST_SCHEMA = "leopard2-lab-manifest/v2"
RESULT_SCHEMA = "leopard2-lab-result/v1"
MERGE_SCHEMA = "leopard2-lab-merged/v1"
JOB_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
TOKEN_REPLACEMENTS = ("seed", "job_id", "cpu_set")


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
    default_timeout = defaults.get("timeout_seconds", 300.0)
    _positive_number(default_timeout, "default timeout_seconds")
    default_memory = defaults.get("memory_mb", _default_memory_mb(topology, selected_workers))
    _positive_int(default_memory, "default memory_mb", allow_zero=True)
    default_cpu_count = defaults.get("cpu_count", 1)
    _positive_int(default_cpu_count, "default cpu_count")
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
        _positive_number(timeout_seconds, "job {} timeout_seconds".format(job_id))
        _positive_int(memory_mb, "job {} memory_mb".format(job_id), allow_zero=True)

        cpu_group = raw.get("cpu_group")
        if cpu_group is not None and (not isinstance(cpu_group, str) or not cpu_group):
            raise LabError("job {} cpu_group must be a non-empty string".format(job_id))
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

        cwd = raw.get("cwd", ".")
        if not isinstance(cwd, str) or not cwd:
            raise LabError("job {} cwd must be a non-empty string".format(job_id))
        minimum_memory_mb = raw.get("minimum_memory_mb", 0)
        _positive_int(
            minimum_memory_mb, "job {} minimum_memory_mb".format(job_id),
            allow_zero=True)
        executable = _resolve_executable(
            command, root_path, cwd, environment, executable_identities)
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
            "minimum_memory_mb": minimum_memory_mb,
            "cpu_set": cpu_set,
            "seed": seed,
            "executable": executable,
        }
        if raw.get("cpu_group") is not None:
            job["cpu_group"] = raw["cpu_group"]
        if benchmark_cell is not None:
            job["benchmark_cell"] = _json_copy(benchmark_cell)
        job["job_digest"] = _digest(job)
        jobs.append(job)

    manifest = {
        "schema": MANIFEST_SCHEMA,
        "root": str(root_path),
        "base_seed": selected_seed,
        "workers": selected_workers,
        "cpu_policy": cpu_policy,
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
        if not isinstance(job.get("cpu_set"), list) or not job["cpu_set"]:
            raise LabError("job {} has no CPU set".format(job["id"]))
        _positive_number(job.get("timeout_seconds"), "job {} timeout_seconds".format(job["id"]))
        _positive_int(job.get("memory_mb"), "job {} memory_mb".format(job["id"]), allow_zero=True)
        _positive_int(
            job.get("minimum_memory_mb"),
            "job {} minimum_memory_mb".format(job["id"]), allow_zero=True)
        executable = job.get("executable")
        if (not isinstance(executable, dict) or
                not isinstance(executable.get("path"), str) or
                not re.match(r"^[0-9a-f]{64}$", str(executable.get("sha256", ""))) or
                not isinstance(executable.get("size_bytes"), int) or
                executable["size_bytes"] < 0):
            raise LabError("job {} has an invalid executable identity".format(job["id"]))
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


def _validate_terminal_result(result_path, result, job):
    unsigned = dict(result)
    expected_result_digest = unsigned.pop("result_digest", None)
    if expected_result_digest != _digest(unsigned):
        raise LabError("terminal result digest is invalid for job {}".format(job["id"]))
    if result.get("job_id") != job["id"] or result.get("cpu_set") != job["cpu_set"]:
        raise LabError("terminal result identity is invalid for job {}".format(job["id"]))
    outputs = result.get("outputs")
    if not isinstance(outputs, dict):
        raise LabError("terminal result output identities are missing for job {}".format(job["id"]))
    job_dir = Path(result_path).parent
    for name in ("stdout", "stderr"):
        expected = outputs.get(name)
        if not isinstance(expected, dict) or _content_identity(
                job_dir / (name + ".txt")) != expected:
            raise LabError("{} changed after job {} completed".format(name, job["id"]))


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


def _expanded_environment(job):
    environment = os.environ.copy()
    environment.update({key: _replace_tokens(value, job) for key, value in job["env"].items()})
    environment["LEO2_LAB_SEED"] = str(job["seed"])
    environment["LEO2_LAB_JOB_ID"] = job["id"]
    environment["LEO2_LAB_CPUSET"] = format_cpu_list(job["cpu_set"])
    return environment


def _child_setup(cpu_set, memory_mb):
    if hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, set(cpu_set))
    elif cpu_set:
        raise RuntimeError("CPU affinity is unavailable on this platform")
    if memory_mb and resource is not None:
        limit = int(memory_mb) * 1024 * 1024
        resource.setrlimit(resource.RLIMIT_AS, (limit, limit))


def _base_result(job, command, duration_seconds, outcome, exit_code=None, detail=None):
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
        "minimum_memory_mb": job.get("minimum_memory_mb", 0),
        "timeout_seconds": job["timeout_seconds"],
        "command": command,
        "stdout": "stdout.txt",
        "stderr": "stderr.txt",
    }
    if detail:
        result["detail"] = detail
    return result


def _write_terminal_result(job_dir, result):
    job_dir = Path(job_dir)
    result["outputs"] = {
        "stdout": _content_identity(job_dir / "stdout.txt"),
        "stderr": _content_identity(job_dir / "stderr.txt"),
    }
    unsigned = dict(result)
    unsigned.pop("result_digest", None)
    result["result_digest"] = _digest(unsigned)
    _atomic_write_json(job_dir / "result.json", result)
    return result


def _job_executable_matches(job):
    current = _file_identity(job["executable"]["path"])
    expected = job["executable"]
    return (current["sha256"] == expected["sha256"] and
            current["size_bytes"] == expected["size_bytes"])


def _launch_job(job, manifest_root, output_dir):
    job_dir = _job_directory(output_dir, job["id"])
    job_dir.mkdir(parents=True, exist_ok=True)
    stdout_path = job_dir / "stdout.txt"
    stderr_path = job_dir / "stderr.txt"
    stdout_handle = stdout_path.open("wb")
    stderr_handle = stderr_path.open("wb")
    command = _expanded_command(job)
    cwd = Path(job["cwd"])
    if not cwd.is_absolute():
        cwd = Path(manifest_root) / cwd
    started = time.monotonic()
    try:
        if not _job_executable_matches(job):
            raise LabError("executable changed before launch: {}".format(
                job["executable"]["path"]))
        process = subprocess.Popen(
            command,
            cwd=str(cwd),
            env=_expanded_environment(job),
            stdin=subprocess.DEVNULL,
            stdout=stdout_handle,
            stderr=stderr_handle,
            start_new_session=True,
            preexec_fn=lambda: _child_setup(job["cpu_set"], job["memory_mb"]),
        )
    except BaseException as error:
        stdout_handle.close()
        stderr_handle.close()
        result = _base_result(job, command, time.monotonic() - started, "launch_error", detail=str(error))
        _write_terminal_result(job_dir, result)
        return None, result
    active = {
        "job": job,
        "command": command,
        "process": process,
        "started": started,
        "stdout_handle": stdout_handle,
        "stderr_handle": stderr_handle,
        "timed_out": False,
        "terminate_started": None,
    }
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


def _finish_active(active, output_dir, forced_outcome=None, detail=None):
    process = active["process"]
    try:
        exit_code = process.wait(timeout=1.0)
    except subprocess.TimeoutExpired:
        _signal_group(process, signal.SIGKILL)
        exit_code = process.wait()
    active["stdout_handle"].close()
    active["stderr_handle"].close()
    executable_detail = None
    try:
        if not _job_executable_matches(active["job"]):
            executable_detail = "executable changed during execution"
    except LabError as error:
        executable_detail = str(error)
    if forced_outcome:
        outcome = forced_outcome
    elif active["timed_out"]:
        outcome = "timeout"
    elif exit_code == 0:
        outcome = "success"
    else:
        outcome = "failed"
    if executable_detail:
        outcome = "evidence_invalid"
        detail = executable_detail if detail is None else detail + "; " + executable_detail
    result = _base_result(
        active["job"], active["command"], time.monotonic() - active["started"],
        outcome, exit_code=exit_code, detail=detail)
    result_path = _job_directory(output_dir, active["job"]["id"]) / "result.json"
    return _write_terminal_result(result_path.parent, result)


def _can_launch(job, active, allow_cpu_overlap, memory_budget_mb=None):
    if memory_budget_mb is not None:
        active_memory_mb = sum(
            current["job"].get("minimum_memory_mb", 0) for current in active)
        if active_memory_mb + job.get("minimum_memory_mb", 0) > memory_budget_mb:
            return False
    if allow_cpu_overlap:
        return True
    active_cpus = set()
    for current in active:
        active_cpus.update(current["job"]["cpu_set"])
    return not active_cpus.intersection(job["cpu_set"])


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


def _record_unavailable(job, output_dir, reason):
    job_dir = _job_directory(output_dir, job["id"])
    job_dir.mkdir(parents=True, exist_ok=True)
    (job_dir / "stdout.txt").write_text("", encoding="utf-8")
    (job_dir / "stderr.txt").write_text(reason + "\n", encoding="utf-8")
    result = _base_result(
        job, _expanded_command(job), 0.0, "unavailable", detail=reason)
    return _write_terminal_result(job_dir, result)


def _dry_run_plan(manifest, output_dir, rerun_failed, workers=None):
    _validate_runtime_cpus(manifest)
    _validate_runtime_executables(manifest)
    worker_count = manifest["workers"] if workers is None else workers
    _positive_int(worker_count, "workers")
    if worker_count > 128:
        raise LabError("workers may not exceed 128")
    memory_budget_mb = _manifest_memory_budget_mb(manifest)
    planned = []
    for job in manifest["jobs"]:
        completed = _read_completed_result(output_dir, job, rerun_failed)
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
            "minimum_memory_mb": job.get("minimum_memory_mb", 0),
            "command": _expanded_command(job),
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

    results = {}
    pending = []
    resumed = 0
    preflight_unavailable = 0
    memory_budget_mb = _manifest_memory_budget_mb(manifest)
    for job in manifest["jobs"]:
        completed = _read_completed_result(output_dir, job, rerun_failed)
        if completed is not None:
            results[job["id"]] = completed
            resumed += 1
        else:
            unavailable_reason = _memory_unavailable_reason(job, memory_budget_mb)
            if unavailable_reason:
                results[job["id"]] = _record_unavailable(
                    job, output_dir, unavailable_reason)
                preflight_unavailable += 1
            else:
                pending.append(job)

    active = []
    total = len(manifest["jobs"])
    run_started = time.monotonic()
    last_progress = 0.0
    interrupted = False
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
                    launched, launch_result = _launch_job(job, manifest["root"], output_dir)
                    if launched is not None:
                        active.append(launched)
                    else:
                        results[job["id"]] = launch_result
                    launched_any = True
                    break

            now = time.monotonic()
            for current in list(active):
                process = current["process"]
                exit_code = process.poll()
                elapsed = now - current["started"]
                if exit_code is None and not current["timed_out"] and elapsed >= current["job"]["timeout_seconds"]:
                    current["timed_out"] = True
                    current["terminate_started"] = now
                    _signal_group(process, signal.SIGTERM)
                elif (exit_code is None and current["timed_out"] and
                      now - current["terminate_started"] >= 0.5):
                    _signal_group(process, signal.SIGKILL)
                if process.poll() is not None:
                    result = _finish_active(current, output_dir)
                    results[current["job"]["id"]] = result
                    active.remove(current)

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
        "executed": total - resumed - preflight_unavailable - len(pending),
        "preflight_unavailable": preflight_unavailable,
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
                "command": [python, "-c", "import os; print(os.environ['LEO2_LAB_JOB_ID'], os.environ['LEO2_LAB_SEED'])"],
            },
            {
                "id": "demo.affinity",
                "command": [python, "-c", "import os; print(sorted(os.sched_getaffinity(0)))"],
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
                    "command": [python, "-c", "import os; print(os.environ['LEO2_LAB_SEED']); print(sorted(os.sched_getaffinity(0)))"],
                },
                {
                    "id": "memory-limit",
                    "command": [
                        python, "-c",
                        "try:\n bytearray(256 * 1024 * 1024)\nexcept MemoryError:\n print('limited')\nelse:\n raise SystemExit(9)",
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
        mutate_code = "from pathlib import Path; Path({!r}).write_text({!r})".format(
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
                    {"id": "aggregate-a", "command": [python, "-c", "pass"],
                     "minimum_memory_mb": per_job_mb},
                    {"id": "aggregate-b", "command": [python, "-c", "pass"],
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
