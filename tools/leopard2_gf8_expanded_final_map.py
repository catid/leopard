#!/usr/bin/env python3
"""Prepare and run the bounded final GF8 Leopard1/Leopard2 comparison.

The broad diagnostic stage is intentionally non-authoritative: it holds the
canonical lock exclusively and uses independent physical cores in parallel.
The confirmation stage holds the same lock exclusively and runs serial ABBA
only for near/slower cells and explicitly selected loss/reuse/batch/tiling
cells.  Both stages are resumable at one JSON document per cell.
"""

import argparse
import concurrent.futures
import contextlib
import fcntl
import hashlib
import json
import math
import os
import platform
import re
import statistics
import subprocess
import sys
import threading
import time
from pathlib import Path


SCHEMA = "leopard2-expanded-final-map/v2"
MATRIX_SCHEMA = "leopard2-expanded-final-matrix/v2"
MAIN_COMMIT = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
BASELINE_SHA256 = (
    "a43d7f43ff2e887ebcd47a1e94f806847a5d8b858a4e383e6c8d5e528a7dd910")
CANONICAL_MATRIX_CELL_COUNT = 341
CANONICAL_MATRIX_SHA256 = (
    "a2c7bc3873ead302e9254c03a7f0ae75bce8d0d48b6f1a565759f4df1a1d6a55")
COMMIT_PATTERN = re.compile(r"^[0-9a-f]{40}$")
ABBA_ORDER = ("baseline", "candidate", "candidate", "baseline")

# This is the frozen 144-cell final-map shape set (48 shapes x three sizes).
CORE_SHAPES = (
    ("t1-k1-r1", 1, 1),
    ("t1-k8-r1", 8, 1),
    ("t1-k129-r1", 129, 1),
    ("t1-k255-r1", 255, 1),
    ("t2-k2-r2", 2, 2),
    ("t2-k3-r2", 3, 2),
    ("t2-k8-r2", 8, 2),
    ("t2-k64-r2", 64, 2),
    ("t2-k254-r2", 254, 2),
    ("t4-k3-r3", 3, 3),
    ("t4-k4-r3", 4, 3),
    ("t4-k7-r3", 7, 3),
    ("t4-k8-r3", 8, 3),
    ("t4-k11-r3", 11, 3),
    ("t4-k12-r3", 12, 3),
    ("t4-k15-r3", 15, 3),
    ("t4-k16-r3", 16, 3),
    ("t4-k64-r3", 64, 3),
    ("t4-k252-r3", 252, 3),
    ("t4-k4-r4", 4, 4),
    ("t4-k7-r4", 7, 4),
    ("t4-k8-r4", 8, 4),
    ("t4-k11-r4", 11, 4),
    ("t4-k12-r4", 12, 4),
    ("t4-k15-r4", 15, 4),
    ("t4-k16-r4", 16, 4),
    ("t4-k64-r4", 64, 4),
    ("t4-k252-r4", 252, 4),
    ("t8-k8-r8", 8, 8),
    ("t8-k9-r8", 9, 8),
    ("t8-k16-r8", 16, 8),
    ("t8-k64-r8", 64, 8),
    ("t8-k248-r8", 248, 8),
    ("t16-k16-r16", 16, 16),
    ("t16-k17-r16", 17, 16),
    ("t16-k32-r16", 32, 16),
    ("t16-k64-r16", 64, 16),
    ("t16-k240-r16", 240, 16),
    ("t32-k32-r32", 32, 32),
    ("t32-k33-r32", 33, 32),
    ("t32-k64-r32", 64, 32),
    ("t32-k128-r32", 128, 32),
    ("t32-k224-r32", 224, 32),
    ("t64-k64-r64", 64, 64),
    ("t64-k65-r64", 65, 64),
    ("t64-k128-r64", 128, 64),
    ("t64-k192-r64", 192, 64),
    ("t128-k128-r128", 128, 128),
)

# Adds transform-boundary and non-power-R cells absent from the old map.
BOUNDARY_SHAPES = (
    ("t4-k5-r4", 5, 4),
    ("t8-k7-r7", 7, 7),
    ("t16-k15-r15", 15, 15),
    ("t32-k31-r31", 31, 31),
    ("t32-k34-r32", 34, 32),
    ("t64-k63-r63", 63, 63),
    ("t128-k127-r127", 127, 127),
)

TINY_SHAPES = (
    (1, 1), (129, 1), (3, 2), (64, 2), (4, 4), (16, 4),
    (7, 7), (9, 8), (15, 15), (17, 16), (31, 31), (34, 32),
    (65, 64), (128, 128),
)
LOSS_SHAPES = (
    (129, 1), (3, 2), (16, 4), (9, 8), (17, 16),
    (34, 32), (128, 32), (65, 64), (192, 64), (128, 128),
)
REUSE_SHAPES = LOSS_SHAPES
BATCH_SHAPES = LOSS_SHAPES
BATCH_TINY_SHAPES = (
    (3, 2), (16, 4), (9, 8), (17, 16), (34, 32), (65, 64),
)
LARGE_SHAPES = (
    (129, 1), (16, 4), (9, 8), (17, 16), (34, 32),
    (128, 32), (65, 64), (192, 64), (128, 128),
)
SELECTOR_ISOLATION_SHAPES = ((33, 32), (34, 32), (35, 32))

# The large-shard R=1 selector deliberately changes implementation at these
# four measured boundaries.  Retain the immediate lower-K fallback alongside
# each promoted cell, and exercise recovery rather than only no-loss decode.
R1_SELECTOR_BOUNDARIES = (
    (50, 51, 1048576),
    (30, 31, 2097152),
    (9, 10, 4194304),
    (7, 8, 8388608),
)

# The AVX2 source-major direct executor has specialized even-output kernels
# and composed odd-output paths.  The ordinary loss sweep contains only
# powers of two, so keep each non-power arity at the dispatch boundary, a
# cache-resident representative, and a large-shard tiling point.
DIRECT_ODD_LOSSES = (3, 5, 6, 7)

BYTE_LABELS = {
    64: "64b", 256: "256b", 1024: "1k", 2048: "2k",
    4096: "4k", 8192: "8k", 16384: "16k", 65536: "64k",
    1048576: "1m", 2097152: "2m", 4194304: "4m",
    8388608: "8m", 16777216: "16m",
}


def ceil_pow2(value):
    return 1 << (value - 1).bit_length()


def canonical_bytes(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def sha256_bytes(value):
    return hashlib.sha256(value).hexdigest()


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp.%d" % os.getpid())
    temporary.write_bytes(canonical_bytes(value) + b"\n")
    os.replace(str(temporary), str(path))


def read_json(path):
    with open(path, "r", encoding="utf-8") as stream:
        return json.load(stream)


def stable_seed(key):
    digest = hashlib.sha256((MATRIX_SCHEMA + ":" + key).encode()).digest()
    return int.from_bytes(digest[:8], "big") & ((1 << 63) - 1)


def iteration_policy(byte_count, confirmation):
    if byte_count <= 1024:
        return (21 if confirmation else 15, 2)
    if byte_count <= 4096:
        return (11 if confirmation else 7, 2)
    if byte_count <= 65536:
        return (9 if confirmation else 5, 2)
    if byte_count <= 1048576:
        return (3 if confirmation else 2, 1)
    return (1, 1)


def cell_key(k, r, byte_count, loss, reuse, batch):
    return "k%d-r%d-%s-l%d-u%d-b%d" % (
        k, r, BYTE_LABELS[byte_count], loss, reuse, batch)


def add_cell(cells, k, r, byte_count, loss, reuse, batch, tag):
    if (k <= 0 or r <= 0 or
            ceil_pow2(k + ceil_pow2(r)) > 256):
        raise RuntimeError("matrix contains an invalid GF8 public code")
    if loss < 0 or loss > min(k, r):
        raise RuntimeError("matrix contains an invalid loss count")
    key = (k, r, byte_count, loss, reuse, batch)
    if key not in cells:
        diagnostic_iterations, diagnostic_warmup = iteration_policy(
            byte_count, False)
        abba_iterations, abba_warmup = iteration_policy(byte_count, True)
        identifier = cell_key(k, r, byte_count, loss, reuse, batch)
        cells[key] = {
            "id": identifier,
            "K": k,
            "R": r,
            "T": ceil_pow2(r),
            "bytes": byte_count,
            "loss": loss,
            "reuse": reuse,
            "batch": batch,
            "seed": stable_seed(identifier),
            "diagnostic_iterations": diagnostic_iterations,
            "diagnostic_warmup": diagnostic_warmup,
            "abba_iterations": abba_iterations,
            "abba_warmup": abba_warmup,
            "estimated_peak_bytes": (
                6 * (k + r) * byte_count * batch + (64 << 20)),
            "tags": [],
        }
    if tag not in cells[key]["tags"]:
        cells[key]["tags"].append(tag)


def loss_values(k, r):
    values = {0, 1, 2, 4, r // 2, r}
    return sorted(value for value in values if value <= min(k, r))


def make_matrix():
    cells = {}
    all_shapes = CORE_SHAPES + BOUNDARY_SHAPES
    for unused_name, k, r in all_shapes:
        for byte_count in (2048, 4096, 65536):
            add_cell(cells, k, r, byte_count, min(k, r), 8, 1, "core")
    for k, r in TINY_SHAPES:
        for byte_count in (64, 256, 1024):
            add_cell(cells, k, r, byte_count, min(k, r), 8, 1,
                     "tiny_shard")
        add_cell(cells, k, r, 8192, min(k, r), 8, 1, "eight_kib")
    for k, r in LOSS_SHAPES:
        for loss in loss_values(k, r):
            add_cell(cells, k, r, 65536, loss, 8, 1, "loss_sweep")
    for k, r in REUSE_SHAPES:
        loss = r // 2
        for reuse in (1, 8, 64):
            add_cell(cells, k, r, 65536, loss, reuse, 1, "reuse_sweep")
    for k, r in BATCH_SHAPES:
        loss = r // 2
        for batch in (1, 8):
            add_cell(cells, k, r, 65536, loss, 8, batch, "batch_sweep")
    for k, r in BATCH_TINY_SHAPES:
        add_cell(cells, k, r, 256, r // 2, 8, 8, "batch_sweep")
    for k, r in LARGE_SHAPES:
        for byte_count in (1048576, 16777216):
            add_cell(cells, k, r, byte_count, r // 2, 8, 1,
                     "large_tiling")
    for k, r in SELECTOR_ISOLATION_SHAPES:
        for byte_count in (16384, 65536):
            add_cell(cells, k, r, byte_count, r, 8, 1,
                     "selector_isolation")
    for lower_k, promoted_k, byte_count in R1_SELECTOR_BOUNDARIES:
        add_cell(cells, lower_k, 1, byte_count, 1, 8, 1,
                 "r1_selector_fallback")
        add_cell(cells, promoted_k, 1, byte_count, 1, 8, 1,
                 "r1_selector_promoted")
    for k in (129, 255):
        add_cell(cells, k, 1, 1048576, 1, 8, 1,
                 "r1_large_recovery")
    for byte_count in (2048, 65536, 1048576):
        for loss in DIRECT_ODD_LOSSES:
            add_cell(cells, 65, 64, byte_count, loss, 8, 1,
                     "direct_odd_arity")
    result = sorted(cells.values(), key=lambda cell: cell["id"])
    for index, cell in enumerate(result):
        cell["index"] = index
        cell["tags"].sort()
    matrix = {
        "schema": MATRIX_SCHEMA,
        "description": (
            "Bounded broad diagnostic plus serial near/loss ABBA confirmation"),
        "baseline_commit": MAIN_COMMIT,
        "baseline_sha256": BASELINE_SHA256,
        "cell_count": len(result),
        "diagnostic_child_invocations": 2 * len(result),
        "maximum_abba_child_invocations_at_three_rounds": 12 * len(result),
        "sizes": sorted({cell["bytes"] for cell in result}),
        "transform_sides": sorted({cell["T"] for cell in result}),
        "losses": sorted({cell["loss"] for cell in result}),
        "reuse_counts": sorted({cell["reuse"] for cell in result}),
        "batch_counts": sorted({cell["batch"] for cell in result}),
        "cells": result,
    }
    return matrix


def matrix_digest(matrix):
    return sha256_bytes(canonical_bytes(matrix))


def validate_matrix(matrix):
    if matrix.get("schema") != MATRIX_SCHEMA:
        raise RuntimeError("matrix schema mismatch")
    regenerated = make_matrix()
    if matrix != regenerated:
        raise RuntimeError("matrix differs from the deterministic built-in matrix")
    if len(matrix["cells"]) != CANONICAL_MATRIX_CELL_COUNT:
        raise RuntimeError("canonical matrix cell count drifted")
    if matrix_digest(matrix) != CANONICAL_MATRIX_SHA256:
        raise RuntimeError("canonical matrix SHA-256 drifted")
    ids = [cell["id"] for cell in matrix["cells"]]
    seeds = [cell["seed"] for cell in matrix["cells"]]
    if len(ids) != len(set(ids)) or len(seeds) != len(set(seeds)):
        raise RuntimeError("matrix IDs or seeds are not unique")
    required_sizes = {64, 256, 1024, 2048, 4096, 8192, 16384, 65536,
                      1048576, 2097152, 4194304, 8388608, 16777216}
    if set(matrix["sizes"]) != required_sizes:
        raise RuntimeError("matrix size coverage drifted")
    if set(matrix["transform_sides"]) != {1, 2, 4, 8, 16, 32, 64, 128}:
        raise RuntimeError("matrix transform-side coverage drifted")
    if not {0, 1, 2, 4}.issubset(set(matrix["losses"])):
        raise RuntimeError("matrix small-loss coverage drifted")
    if matrix["reuse_counts"] != [1, 8, 64] or matrix["batch_counts"] != [1, 8]:
        raise RuntimeError("matrix reuse/batch coverage drifted")


def checked_output(command, cwd=None, timeout=60):
    completed = subprocess.run(
        command, cwd=cwd, stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        timeout=timeout, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            "%r failed rc=%d: %s" %
            (command, completed.returncode, completed.stderr.strip()))
    return completed.stdout.strip()


def validate_executable(path, expected_sha256=None):
    resolved = os.path.realpath(path)
    if not os.path.isfile(resolved) or not os.access(resolved, os.X_OK):
        raise RuntimeError("benchmark is not executable: " + resolved)
    digest = sha256_file(resolved)
    if expected_sha256 is not None and digest != expected_sha256:
        raise RuntimeError(
            "frozen baseline SHA-256 mismatch: expected %s, got %s" %
            (expected_sha256, digest))
    return {"path": resolved, "sha256": digest}


def inspect_isa(path):
    completed = subprocess.run(
        ["/usr/bin/objdump", "-d", "-M", "intel", path],
        stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, text=True, timeout=180, check=False)
    if completed.returncode != 0:
        raise RuntimeError("objdump failed: " + completed.stderr.strip())
    evex = 0
    ymm = 0
    instruction = re.compile(r"^\s*[0-9a-f]+:\s+((?:[0-9a-f]{2}\s+)+)")
    for line in completed.stdout.splitlines():
        match = instruction.match(line)
        if match and match.group(1).split()[0] == "62":
            evex += 1
        if re.search(r"\bymm[0-9]+\b", line):
            ymm += 1
    if evex:
        raise RuntimeError("%s contains %d EVEX instructions" % (path, evex))
    if not ymm:
        raise RuntimeError("%s contains no visible YMM instructions" % path)
    return {
        "disassembler": checked_output(
            ["/usr/bin/objdump", "--version"]).splitlines()[0],
        "evex_prefixed_instruction_count": evex,
        "ymm_operand_instruction_count": ymm,
    }


def source_provenance(source, expected_commit):
    source = os.path.realpath(source)
    if not COMMIT_PATTERN.match(expected_commit):
        raise RuntimeError("candidate commit must be a full lowercase SHA")
    head = checked_output(["git", "rev-parse", "HEAD"], cwd=source)
    tree = checked_output(["git", "rev-parse", "HEAD^{tree}"], cwd=source)
    expected_tree = checked_output(
        ["git", "rev-parse", expected_commit + "^{tree}"], cwd=source)
    status = checked_output(
        ["git", "status", "--porcelain", "--untracked-files=no"], cwd=source)
    if head != expected_commit or tree != expected_tree:
        raise RuntimeError(
            "candidate source does not exactly match expected commit")
    if status:
        raise RuntimeError("candidate source has tracked changes:\n" + status)
    return {"path": source, "head": head, "tree": tree,
            "tracked_clean": True}


def cmake_provenance(build, source):
    build = os.path.realpath(build)
    cache = os.path.join(build, "CMakeCache.txt")
    if not os.path.isfile(cache):
        raise RuntimeError("candidate build lacks CMakeCache.txt: " + cache)
    values = {}
    with open(cache, encoding="utf-8", errors="replace") as stream:
        for raw in stream:
            if "=" not in raw:
                continue
            left, value = raw.rstrip("\n").split("=", 1)
            if left in ("CMAKE_HOME_DIRECTORY:INTERNAL",
                        "CMAKE_GENERATOR:INTERNAL",
                        "CMAKE_BUILD_TYPE:STRING",
                        "CMAKE_CXX_COMPILER:FILEPATH",
                        "CMAKE_CXX_COMPILER_ID:STRING",
                        "CMAKE_CXX_COMPILER_VERSION:STRING"):
                values[left] = value
    home = values.get("CMAKE_HOME_DIRECTORY:INTERNAL")
    if os.path.realpath(home or "") != os.path.realpath(source):
        raise RuntimeError("candidate build points at another source tree")
    return {
        "path": build, "cache": cache,
        "cache_sha256": sha256_file(cache),
        "source": os.path.realpath(home),
        "generator": values.get("CMAKE_GENERATOR:INTERNAL"),
        "build_type": values.get("CMAKE_BUILD_TYPE:STRING"),
        "compiler": values.get("CMAKE_CXX_COMPILER:FILEPATH"),
        "compiler_id": values.get("CMAKE_CXX_COMPILER_ID:STRING"),
        "compiler_version": values.get("CMAKE_CXX_COMPILER_VERSION:STRING"),
    }


def candidate_source_attestation(executable, source):
    cpu = sorted(os.sched_getaffinity(0))[0]
    command = [
        "/usr/bin/taskset", "-c", str(cpu), executable,
        "--k", "3", "--r", "2", "--bytes", "64", "--loss", "1",
        "--batch", "1", "--reuse", "1", "--iterations", "1",
        "--warmup", "0", "--threads", "1", "--seed", "7",
        "--profile", "high", "--field", "gf8", "--backend", "avx2",
        "--skip-legacy", "--retain-samples", "--attest-source",
        "--json", "-",
    ]
    completed = subprocess.run(
        command, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, env=benchmark_environment(), timeout=120,
        check=False)
    if completed.returncode != 0:
        raise InvocationError(
            "candidate source-attestation preflight failed",
            {"command": command, "returncode": completed.returncode,
             "stdout": completed.stdout.decode(errors="replace"),
             "stderr": completed.stderr.decode(errors="replace")})
    try:
        document = json.loads(completed.stdout)
    except Exception as error:
        raise InvocationError(
            "candidate source-attestation JSON is invalid: %s" % error,
            {"command": command,
             "stdout": completed.stdout.decode(errors="replace"),
             "stderr": completed.stderr.decode(errors="replace")})
    build = document.get("build", {})
    if document.get("schema") != "leopard2-benchmark-v5":
        raise RuntimeError("candidate attestation did not emit schema v5")
    if build.get("source_commit") != source["head"]:
        raise RuntimeError("candidate binary attests another source commit")
    if build.get("source_tree") != source["tree"]:
        raise RuntimeError("candidate binary attests another source tree")
    if build.get("source_tracked_dirty") is not False:
        raise RuntimeError("candidate binary was configured from dirty source")
    parameters = document["parameters"]
    resolved = document["resolved"]
    if (parameters["requested_profile"], parameters["requested_field"],
            parameters["requested_backend"]) != (
                "legacy_high_v1", "gf8", "avx2"):
        raise RuntimeError("candidate attestation request identity changed")
    if not parameters.get("skip_legacy") or \
            not parameters.get("attest_source"):
        raise RuntimeError("candidate attestation flags were not recorded")
    if (resolved["profile"], resolved["field"], resolved["backend"],
            resolved["parent_count"], resolved["padded_side"],
            resolved["thread_count"]) != (
                "legacy_high_v1", "gf8", "avx2", 8, 2, 1):
        raise RuntimeError("candidate attestation resolved identity changed")
    correctness = document["correctness"]
    if not correctness.get("leopard2_round_trip") or \
            correctness.get("legacy_comparison") is not None:
        raise RuntimeError("candidate attestation round trip failed")
    stable = {
        "command": command,
        "schema": document["schema"], "build": build,
        "parameters": parameters, "resolved": resolved,
        "correctness": correctness,
        "workload_digests": document["workload_digests"],
    }
    stable["stable_identity_sha256"] = sha256_bytes(canonical_bytes(stable))
    return stable


def parse_cpu_list(text):
    result = set()
    for item in text.strip().split(","):
        if not item:
            continue
        bounds = item.split("-", 1)
        if len(bounds) == 1:
            result.add(int(bounds[0]))
        else:
            result.update(range(int(bounds[0]), int(bounds[1]) + 1))
    return result


def sibling_set(cpu):
    path = ("/sys/devices/system/cpu/cpu%d/topology/"
            "thread_siblings_list" % cpu)
    with open(path, encoding="ascii") as stream:
        return parse_cpu_list(stream.read())


def physical_cpu_pairs():
    allowed = set(os.sched_getaffinity(0))
    groups = {}
    for cpu in sorted(allowed):
        siblings = sibling_set(cpu) & allowed
        if len(siblings) < 2:
            continue
        key = tuple(sorted(siblings))
        groups[key] = key
    pairs = [(group[0], group[1]) for group in sorted(groups)]
    if not pairs:
        raise RuntimeError("no allowed SMT sibling pairs were detected")
    return pairs


def parse_pairs(text):
    if not text:
        return physical_cpu_pairs()
    pairs = []
    for item in text.split(","):
        values = item.split(":")
        if len(values) != 2:
            raise RuntimeError("CPU pairs use cpu:sibling syntax")
        pairs.append((int(values[0]), int(values[1])))
    for cpu, sibling in pairs:
        validate_cpu_pair(cpu, sibling)
    if len({value for pair in pairs for value in pair}) != 2 * len(pairs):
        raise RuntimeError("CPU pairs overlap")
    return pairs


def validate_cpu_pair(cpu, sibling):
    allowed = set(os.sched_getaffinity(0))
    if cpu == sibling or cpu not in allowed or sibling not in allowed:
        raise RuntimeError("benchmark CPU pair is not distinct and allowed")
    siblings = sibling_set(cpu) & allowed
    if cpu not in siblings or sibling not in siblings:
        raise RuntimeError("CPU%d and CPU%d are not SMT siblings" %
                           (cpu, sibling))
    if (sibling_set(sibling) & allowed) != siblings:
        raise RuntimeError("asymmetric sibling topology")
    return {"benchmark_cpu": cpu, "reserved_sibling": sibling,
            "thread_siblings": sorted(siblings),
            "process_allowed_cpus": sorted(allowed)}


def cpu_stat(cpu):
    prefix = "cpu%d " % cpu
    with open("/proc/stat", encoding="ascii") as stream:
        for line in stream:
            if line.startswith(prefix):
                values = [int(value) for value in line.split()[1:]]
                while len(values) < 8:
                    values.append(0)
                names = ("user", "nice", "system", "idle", "iowait",
                         "irq", "softirq", "steal")
                return dict(zip(names, values[:8]))
    raise RuntimeError("missing CPU%d in /proc/stat" % cpu)


def nonidle_delta(before, after):
    names = ("user", "nice", "system", "irq", "softirq", "steal")
    return sum(after[name] - before[name] for name in names)


@contextlib.contextmanager
def benchmark_lock(path, exclusive):
    descriptor = os.open(path, os.O_CREAT | os.O_RDWR, 0o600)
    try:
        mode = fcntl.LOCK_EX if exclusive else fcntl.LOCK_SH
        print("waiting for %s benchmark lock %s" %
              ("exclusive" if exclusive else "shared", path), flush=True)
        fcntl.flock(descriptor, mode)
        print("acquired benchmark lock", flush=True)
        yield
    finally:
        fcntl.flock(descriptor, fcntl.LOCK_UN)
        os.close(descriptor)


def mem_available_bytes():
    with open("/proc/meminfo", encoding="ascii") as stream:
        for line in stream:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) * 1024
    raise RuntimeError("MemAvailable is unavailable")


def mem_total_bytes():
    with open("/proc/meminfo", encoding="ascii") as stream:
        for line in stream:
            if line.startswith("MemTotal:"):
                return int(line.split()[1]) * 1024
    raise RuntimeError("MemTotal is unavailable")


class InvocationError(RuntimeError):
    def __init__(self, message, details):
        RuntimeError.__init__(self, message)
        self.details = details


def benchmark_environment():
    return {
        "PATH": "/usr/bin:/bin", "LANG": "C", "LC_ALL": "C",
        "TZ": "UTC", "OMP_DYNAMIC": "FALSE", "OMP_NUM_THREADS": "1",
        "OMP_PLACES": "cores", "OMP_PROC_BIND": "TRUE",
    }


def invoke_benchmark(kind, executable, cpu, sibling, cell, iterations,
                     warmup, maximum_attempts, candidate_commit,
                     require_idle_sibling):
    common = [
        "--k", str(cell["K"]), "--r", str(cell["R"]),
        "--bytes", str(cell["bytes"]), "--loss", str(cell["loss"]),
        "--batch", str(cell["batch"]), "--reuse", str(cell["reuse"]),
        "--iterations", str(iterations), "--warmup", str(warmup),
        "--threads", "1", "--seed", str(cell["seed"]), "--json", "-",
    ]
    if kind == "candidate":
        common += [
            "--profile", "high", "--field", "gf8", "--backend", "avx2",
            "--skip-legacy", "--retain-samples", "--report-decode-path",
        ]
    command = ["/usr/bin/taskset", "-c", str(cpu), executable] + common
    rejected = []
    timeout = 1200 if cell["bytes"] >= 1048576 else 300
    for attempt in range(maximum_attempts):
        cpu_before = cpu_stat(cpu)
        sibling_before = cpu_stat(sibling)
        started_ns = time.monotonic_ns()
        try:
            process = subprocess.run(
                command, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, env=benchmark_environment(),
                timeout=timeout, check=False)
        except subprocess.TimeoutExpired as error:
            details = {
                "kind": kind, "command": command, "timeout_seconds": timeout,
                "stdout": (error.stdout or b"").decode(errors="replace"),
                "stderr": (error.stderr or b"").decode(errors="replace"),
            }
            raise InvocationError("benchmark timed out", details)
        ended_ns = time.monotonic_ns()
        cpu_after = cpu_stat(cpu)
        sibling_after = cpu_stat(sibling)
        sibling_work = nonidle_delta(sibling_before, sibling_after)
        if process.returncode != 0:
            details = {
                "kind": kind, "command": command,
                "returncode": process.returncode,
                "stdout": process.stdout.decode(errors="replace"),
                "stderr": process.stderr.decode(errors="replace"),
            }
            raise InvocationError("benchmark child failed", details)
        if sibling_work and require_idle_sibling:
            rejected.append({
                "attempt": attempt,
                "reason": "reserved sibling performed non-idle work",
                "sibling_nonidle_jiffies": sibling_work,
                "started_ns": started_ns, "ended_ns": ended_ns,
                "stdout_sha256": hashlib.sha256(process.stdout).hexdigest(),
                "stderr_sha256": hashlib.sha256(process.stderr).hexdigest(),
            })
            time.sleep(0.05)
            continue
        try:
            document = json.loads(process.stdout)
        except Exception as error:
            details = {
                "kind": kind, "command": command,
                "stdout": process.stdout.decode(errors="replace"),
                "stderr": process.stderr.decode(errors="replace"),
            }
            raise InvocationError("invalid benchmark JSON: %s" % error,
                                  details)
        parameters = document["parameters"]
        observed = (
            parameters["K"], parameters["R"], parameters["shard_bytes"],
            parameters["loss_count"], parameters["batch"],
            parameters["reuse"], parameters["seed"])
        expected = (
            cell["K"], cell["R"], cell["bytes"], cell["loss"],
            cell["batch"], cell["reuse"], cell["seed"])
        if observed != expected:
            raise RuntimeError("benchmark parameter mismatch for " + cell["id"])
        resolved = document["resolved"]
        if (resolved["profile"], resolved["field"]) != (
                "legacy_high_v1", "gf8"):
            raise RuntimeError("wire profile or field mismatch")
        if parameters["thread_count"] != 1 or resolved["thread_count"] != 1:
            raise RuntimeError("benchmark did not remain single-threaded")
        if (resolved["parent_count"], resolved["padded_side"]) != (
                ceil_pow2(cell["K"] + cell["T"]), cell["T"]):
            raise RuntimeError("benchmark parent geometry mismatch")
        if kind == "baseline":
            if document["build"]["main_source_commit"] != MAIN_COMMIT:
                raise RuntimeError("baseline commit mismatch")
            if not document["correctness"]["round_trip"]:
                raise RuntimeError("baseline round trip failed")
            encode = document["metrics"]["encode_execution"][
                "median_us_per_batch_call"]
            decode = document["metrics"]["decode_including_setup"][
                "median_us_per_batch_call"]
            times = {"encode": encode, "decode_first": decode,
                     "decode_reuse": decode}
        else:
            if parameters["requested_backend"] != "avx2" or \
                    resolved["backend"] != "avx2":
                raise RuntimeError("candidate did not resolve explicit AVX2")
            runtime_commit = document.get("build", {}).get("source_commit")
            if runtime_commit is not None and runtime_commit != candidate_commit:
                raise RuntimeError("candidate runtime commit mismatch")
            if not document["correctness"]["leopard2_round_trip"]:
                raise RuntimeError("candidate round trip failed")
            encode = document["metrics"]["encode_execution"][
                "median_us_per_batch_call"]
            execution = document["metrics"]["decode_execution"][
                "median_us_per_batch_call"]
            setup = document["metrics"]["decode_plan_setup"]["median_us"]
            times = {
                "encode": encode,
                "decode_first": execution + setup,
                "decode_reuse": execution + setup / cell["reuse"],
            }
        if any(not isinstance(value, (int, float)) or value <= 0
               for value in times.values()):
            raise RuntimeError("non-positive benchmark timing")
        return {
            "kind": kind, "attempt": attempt, "command": command,
            "started_ns": started_ns, "ended_ns": ended_ns,
            "cpu_nonidle_jiffies": nonidle_delta(cpu_before, cpu_after),
            "sibling_nonidle_jiffies": sibling_work,
            "idle_sibling_required": require_idle_sibling,
            "stdout_sha256": hashlib.sha256(process.stdout).hexdigest(),
            "stderr_sha256": hashlib.sha256(process.stderr).hexdigest(),
            "rejected_attempts": rejected,
            "times_us": times, "document": document,
        }
    raise InvocationError(
        "reserved sibling contamination persisted",
        {"kind": kind, "command": command, "rejected_attempts": rejected})


def digest_tuple(document):
    value = document["workload_digests"]
    return (value["algorithm"], value["original_data"],
            value["transmitted_parity"], value["recovered_originals"])


def validate_call_identity(calls, cell):
    expected_digest = None
    expected_missing = None
    candidate_route = None
    for call in calls:
        document = call["document"]
        digest = digest_tuple(document)
        missing = tuple(document["parameters"]["missing_original_indices"])
        if expected_digest is None:
            expected_digest, expected_missing = digest, missing
        if digest != expected_digest or missing != expected_missing:
            raise RuntimeError(
                "wire/recovery digest or missing-set mismatch for " + cell["id"])
        if len(missing) != cell["loss"]:
            raise RuntimeError("missing-set cardinality mismatch")
        if call["kind"] == "candidate":
            route = (document["resolved"]["selected_decode_path"],
                     document["resolved"]["selected_decode_rule"])
            if candidate_route is None:
                candidate_route = route
            elif route != candidate_route:
                raise RuntimeError("candidate decode route changed")
    if expected_digest is None or candidate_route is None:
        raise RuntimeError("incomplete cell invocation set")
    return {
        "workload_digests": {
            "algorithm": expected_digest[0],
            "original_data": expected_digest[1],
            "transmitted_parity": expected_digest[2],
            "recovered_originals": expected_digest[3],
        },
        "missing_original_indices": list(expected_missing),
        "selected_decode_path": candidate_route[0],
        "selected_decode_rule": candidate_route[1],
    }


def point_ratios(calls):
    baseline = next(call for call in calls if call["kind"] == "baseline")
    candidate = next(call for call in calls if call["kind"] == "candidate")
    return {
        metric: baseline["times_us"][metric] / candidate["times_us"][metric]
        for metric in ("encode", "decode_first", "decode_reuse")
    }


def confidence(log_values):
    count = len(log_values)
    if count < 2:
        raise RuntimeError("confidence interval needs at least two rounds")
    critical_table = {2: 12.706204736, 3: 4.302652730, 4: 3.182446305,
                      5: 2.776445105, 6: 2.570581836, 7: 2.446911851,
                      8: 2.364624252, 9: 2.306004135, 10: 2.262157163}
    critical = critical_table.get(count, 1.96)
    mean = statistics.mean(log_values)
    half = critical * statistics.stdev(log_values) / math.sqrt(count)
    return {
        "ratio": math.exp(mean), "ci95_low": math.exp(mean - half),
        "ci95_high": math.exp(mean + half),
        "round_log_contrasts": log_values,
        "orientation": "leopard1_time_over_leopard2_time",
    }


def abba_metrics(calls, rounds):
    if len(calls) != rounds * len(ABBA_ORDER):
        raise RuntimeError("ABBA call count mismatch")
    result = {}
    for metric in ("encode", "decode_first", "decode_reuse"):
        contrasts = []
        for round_index in range(rounds):
            group = calls[round_index * 4:(round_index + 1) * 4]
            baseline = [math.log(call["times_us"][metric]) for call in group
                        if call["kind"] == "baseline"]
            candidate = [math.log(call["times_us"][metric]) for call in group
                         if call["kind"] == "candidate"]
            contrasts.append(statistics.mean(baseline) -
                             statistics.mean(candidate))
        result[metric] = confidence(contrasts)
    return result


def utc_now():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def machine_provenance():
    lscpu = json.loads(checked_output(["/usr/bin/lscpu", "-J"]))
    # Current CPU-frequency percentage is an observation, not platform
    # identity.  Excluding it makes an interrupted run safely resumable.
    lscpu["lscpu"] = [
        item for item in lscpu["lscpu"]
        if item.get("field", "").strip() != "CPU(s) scaling MHz:"]
    return {
        "platform": platform.platform(), "uname": list(platform.uname()),
        "python": sys.version, "python_executable": sys.executable,
        "allowed_cpus": sorted(os.sched_getaffinity(0)),
        "lscpu_json": lscpu,
        "mem_total_bytes": mem_total_bytes(),
    }


def common_provenance(args, matrix):
    baseline = validate_executable(args.baseline, BASELINE_SHA256)
    candidate = validate_executable(args.candidate)
    source = source_provenance(args.candidate_source, args.candidate_commit)
    build = cmake_provenance(args.candidate_build, args.candidate_source)
    attestation = candidate_source_attestation(candidate["path"], source)
    return {
        "schema": SCHEMA,
        "runner": {"path": os.path.realpath(__file__),
                   "sha256": sha256_file(__file__)},
        "matrix": {"path": os.path.realpath(args.matrix),
                   "sha256": matrix_digest(matrix),
                   "file_sha256": sha256_file(args.matrix),
                   "cell_count": len(matrix["cells"])},
        "frozen_main_commit": MAIN_COMMIT,
        "baseline": dict(baseline, isa=inspect_isa(baseline["path"])),
        "candidate": dict(candidate, isa=inspect_isa(candidate["path"]),
                          source_attestation=attestation),
        "candidate_source": source,
        "candidate_build": build,
        "machine": machine_provenance(),
    }


def stage_identity(stage, provenance, options):
    return {"schema": SCHEMA, "stage": stage, "provenance": provenance,
            "options": options}


def prepare_run_directory(run_dir, identity):
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    identity_path = run_dir / "identity.json"
    if identity_path.exists():
        existing = read_json(identity_path)
        if existing != identity:
            raise RuntimeError(
                "run directory identity differs; use a new run directory")
    else:
        atomic_json(identity_path, identity)
    (run_dir / "cells").mkdir(exist_ok=True)
    (run_dir / "failures").mkdir(exist_ok=True)
    return run_dir, sha256_bytes(canonical_bytes(identity))


def source_still_frozen(provenance):
    expected = provenance["candidate_source"]
    observed = source_provenance(expected["path"], expected["head"])
    if observed != expected:
        raise RuntimeError("candidate source provenance changed during stage")
    candidate_digest = sha256_file(provenance["candidate"]["path"])
    if candidate_digest != provenance["candidate"]["sha256"]:
        raise RuntimeError("candidate benchmark changed during stage")
    if sha256_file(provenance["baseline"]["path"]) != BASELINE_SHA256:
        raise RuntimeError("baseline benchmark changed during stage")


def completed_result(path, identity_digest, cell, stage):
    if not path.exists():
        return None
    value = read_json(path)
    if value.get("status") != "complete":
        return None
    if (value.get("identity_sha256") != identity_digest or
            value.get("stage") != stage or value.get("cell") != cell):
        raise RuntimeError("resumed cell identity mismatch: " + str(path))
    return value


def retain_failure(run_dir, stage, cell, pair, identity_digest, error):
    stamp = "%d-%d" % (time.time_ns(), os.getpid())
    details = getattr(error, "details", None)
    value = {
        "schema": SCHEMA, "status": "failed", "stage": stage,
        "identity_sha256": identity_digest, "cell": cell,
        "cpu_pair": list(pair), "error_type": type(error).__name__,
        "error": str(error), "details": details, "failed_at": utc_now(),
    }
    atomic_json(Path(run_dir) / "failures" /
                (cell["id"] + "-" + stamp + ".json"), value)


class MemoryGate:
    def __init__(self, budget):
        self.budget = budget
        self.used = 0
        self.condition = threading.Condition()

    @contextlib.contextmanager
    def hold(self, requested):
        amount = min(requested, self.budget)
        with self.condition:
            while self.used + amount > self.budget:
                self.condition.wait()
            self.used += amount
        try:
            yield
        finally:
            with self.condition:
                self.used -= amount
                self.condition.notify_all()


def diagnostic_cell(cell, pair, args, identity_digest, run_dir, memory_gate):
    result_path = Path(run_dir) / "cells" / (cell["id"] + ".json")
    resumed = completed_result(result_path, identity_digest, cell, "diagnostic")
    if resumed is not None:
        return resumed, True
    cpu, sibling = pair
    try:
        with memory_gate.hold(cell["estimated_peak_bytes"]):
            order = ("baseline", "candidate")
            if cell["seed"] & 1:
                order = tuple(reversed(order))
            calls = []
            for kind in order:
                executable = (args.baseline if kind == "baseline"
                              else args.candidate)
                calls.append(invoke_benchmark(
                    kind, executable, cpu, sibling, cell,
                    cell["diagnostic_iterations"],
                    cell["diagnostic_warmup"], args.maximum_attempts,
                    args.candidate_commit, False))
        identity = validate_call_identity(calls, cell)
        value = {
            "schema": SCHEMA, "status": "complete", "stage": "diagnostic",
            "authoritative": False, "identity_sha256": identity_digest,
            "cell": cell, "cpu_pair": list(pair), "calls": calls,
            "ratios": point_ratios(calls), **identity,
        }
        atomic_json(result_path, value)
        return value, False
    except Exception as error:
        retain_failure(run_dir, "diagnostic", cell, pair, identity_digest, error)
        raise


def diagnostic_summary(run_dir, identity, identity_digest, matrix, results,
                       resumed_count, elapsed_seconds):
    rows = sorted(results, key=lambda value: value["cell"]["id"])
    summary = {
        "schema": SCHEMA, "status": "complete", "stage": "diagnostic",
        "authoritative": False,
        "note": "Parallel screen only; ratios are not publication evidence.",
        "identity_sha256": identity_digest,
        "matrix_sha256": identity["provenance"]["matrix"]["sha256"],
        "cell_count": len(rows), "resumed_cell_count": resumed_count,
        "accepted_child_invocations": 2 * len(rows),
        "accepted_calls_with_sibling_activity": sum(
            int(call["sibling_nonidle_jiffies"] != 0)
            for row in results for call in row["calls"]),
        "elapsed_seconds": elapsed_seconds,
        "candidate_binary_sha256": identity["provenance"]["candidate"]["sha256"],
        "baseline_binary_sha256": identity["provenance"]["baseline"]["sha256"],
        "rows": [{
            "id": row["cell"]["id"], "cell": row["cell"],
            "ratios": row["ratios"], "cpu_pair": row["cpu_pair"],
            "selected_decode_path": row["selected_decode_path"],
            "selected_decode_rule": row["selected_decode_rule"],
            "workload_digests": row["workload_digests"],
            "missing_original_indices": row["missing_original_indices"],
        } for row in rows],
    }
    atomic_json(Path(run_dir) / "summary.json", summary)
    return summary


def run_diagnostic(args):
    matrix = read_json(args.matrix)
    validate_matrix(matrix)
    pairs = parse_pairs(args.cpu_pairs)
    if args.maximum_workers:
        pairs = pairs[:args.maximum_workers]
    if not pairs:
        raise RuntimeError("no CPU pairs selected")
    budget = args.memory_budget_gib * (1 << 30) if args.memory_budget_gib else \
        min(mem_total_bytes() // 4, 32 << 30)
    with benchmark_lock(args.lock, True):
        provenance = common_provenance(args, matrix)
        options = {
            "cpu_pairs": [list(pair) for pair in pairs],
            "maximum_attempts": args.maximum_attempts,
            "memory_budget_bytes": budget,
            "mode": "exclusive_parallel_non_authoritative_screen",
            "idle_sibling_required": False,
        }
        identity = stage_identity("diagnostic", provenance, options)
        run_dir, identity_digest = prepare_run_directory(args.run_dir, identity)
        memory_gate = MemoryGate(budget)
        started = time.monotonic()
        ordered_cells = sorted(
            matrix["cells"],
            key=lambda cell: (cell["estimated_peak_bytes"], cell["id"]))
        assignments = [[] for unused_pair in pairs]
        for index, cell in enumerate(ordered_cells):
            assignments[index % len(pairs)].append(cell)
        progress_lock = threading.Lock()
        progress = [0]

        def worker(pair, assigned_cells):
            local_results = []
            local_resumed = 0
            local_errors = []
            for cell in assigned_cells:
                failed = False
                try:
                    result, resumed = diagnostic_cell(
                        cell, pair, args, identity_digest, run_dir, memory_gate)
                    local_results.append(result)
                    local_resumed += int(resumed)
                except Exception as error:
                    failed = True
                    local_errors.append((cell["id"], error))
                with progress_lock:
                    progress[0] += 1
                    print("diagnostic %d/%d %s%s" % (
                        progress[0], len(matrix["cells"]), cell["id"],
                        " FAILED" if failed else ""), flush=True)
            return local_results, local_resumed, local_errors

        results = []
        resumed_count = 0
        errors = []
        with concurrent.futures.ThreadPoolExecutor(
                max_workers=len(pairs), thread_name_prefix="gf8-map") as executor:
            futures = [
                executor.submit(worker, pair, assignment)
                for pair, assignment in zip(pairs, assignments)]
            for future in concurrent.futures.as_completed(futures):
                worker_results, worker_resumed, worker_errors = future.result()
                results.extend(worker_results)
                resumed_count += worker_resumed
                errors.extend(worker_errors)
        if errors:
            raise RuntimeError(
                "%d diagnostic cells failed; first %s: %s" %
                (len(errors), errors[0][0], errors[0][1]))
        source_still_frozen(provenance)
        summary = diagnostic_summary(
            run_dir, identity, identity_digest, matrix, results,
            resumed_count, time.monotonic() - started)
        print(json.dumps({
            "stage": "diagnostic", "run_dir": str(run_dir),
            "cells": summary["cell_count"],
            "elapsed_seconds": summary["elapsed_seconds"],
            "summary_sha256": sha256_file(run_dir / "summary.json"),
        }, sort_keys=True))


def load_diagnostic_summary(path, matrix):
    summary = read_json(path)
    if (summary.get("schema") != SCHEMA or
            summary.get("stage") != "diagnostic" or
            summary.get("status") != "complete"):
        raise RuntimeError("invalid diagnostic summary")
    if summary.get("matrix_sha256") != matrix_digest(matrix):
        raise RuntimeError("diagnostic summary matrix mismatch")
    ids = {row["id"] for row in summary["rows"]}
    expected = {cell["id"] for cell in matrix["cells"]}
    if ids != expected or len(summary["rows"]) != len(expected):
        raise RuntimeError("diagnostic summary is incomplete")
    return summary


def selected_abba_cells(matrix, diagnostic, threshold, include, exclude):
    row_by_id = {row["id"]: row for row in diagnostic["rows"]}
    include = set(include)
    exclude = set(exclude)
    unknown = (include | exclude) - set(row_by_id)
    if unknown:
        raise RuntimeError("unknown selected cell IDs: " + ", ".join(sorted(unknown)))
    mandatory_tags = {
        "loss_sweep", "reuse_sweep", "batch_sweep", "large_tiling",
        "selector_isolation", "r1_selector_fallback",
        "r1_selector_promoted", "r1_large_recovery", "direct_odd_arity",
    }
    selected = []
    for cell in matrix["cells"]:
        reasons = []
        ratios = row_by_id[cell["id"]]["ratios"]
        near_metrics = sorted(
            metric for metric, ratio in ratios.items() if ratio <= threshold)
        if near_metrics:
            reasons.append("diagnostic_ratio_at_or_below_%g:%s" %
                           (threshold, ",".join(near_metrics)))
        tags = sorted(set(cell["tags"]) & mandatory_tags)
        if tags:
            reasons.append("mandatory_dimension:" + ",".join(tags))
        if cell["id"] in include:
            reasons.append("explicit_include")
        if reasons and cell["id"] not in exclude:
            value = dict(cell)
            value["abba_selection_reasons"] = reasons
            selected.append(value)
    return selected


def abba_cell(cell, pair, args, identity_digest, run_dir):
    result_path = Path(run_dir) / "cells" / (cell["id"] + ".json")
    resumed = completed_result(result_path, identity_digest, cell, "abba")
    if resumed is not None:
        return resumed, True
    cpu, sibling = pair
    calls = []
    try:
        for unused_round in range(args.rounds):
            for kind in ABBA_ORDER:
                executable = (args.baseline if kind == "baseline"
                              else args.candidate)
                calls.append(invoke_benchmark(
                    kind, executable, cpu, sibling, cell,
                    cell["abba_iterations"], cell["abba_warmup"],
                    args.maximum_attempts, args.candidate_commit, True))
        exact = validate_call_identity(calls, cell)
        value = {
            "schema": SCHEMA, "status": "complete", "stage": "abba",
            "authoritative": True, "identity_sha256": identity_digest,
            "cell": cell, "cpu_pair": list(pair), "rounds": args.rounds,
            "order_per_round": list(ABBA_ORDER), "calls": calls,
            "metrics": abba_metrics(calls, args.rounds), **exact,
        }
        atomic_json(result_path, value)
        return value, False
    except Exception as error:
        retain_failure(run_dir, "abba", cell, pair, identity_digest, error)
        raise


def run_abba(args):
    matrix = read_json(args.matrix)
    validate_matrix(matrix)
    diagnostic = load_diagnostic_summary(args.diagnostic_summary, matrix)
    selected = selected_abba_cells(
        matrix, diagnostic, args.near_ratio, args.include_cell,
        args.exclude_cell)
    if not selected:
        raise RuntimeError("ABBA selection is empty")
    if (args.cpu is None) != (args.sibling is None):
        raise RuntimeError("--cpu and --sibling must be supplied together")
    pair = ((args.cpu, args.sibling) if args.cpu is not None
            else physical_cpu_pairs()[0])
    topology = validate_cpu_pair(*pair)
    with benchmark_lock(args.lock, True):
        provenance = common_provenance(args, matrix)
        if diagnostic["candidate_binary_sha256"] != \
                provenance["candidate"]["sha256"]:
            raise RuntimeError("diagnostic used another candidate binary")
        if diagnostic["baseline_binary_sha256"] != \
                provenance["baseline"]["sha256"]:
            raise RuntimeError("diagnostic used another baseline binary")
        options = {
            "cpu_pair": list(pair), "topology": topology,
            "rounds": args.rounds, "order_per_round": list(ABBA_ORDER),
            "maximum_attempts": args.maximum_attempts,
            "idle_sibling_required": True,
            "near_ratio": args.near_ratio,
            "diagnostic_summary": os.path.realpath(args.diagnostic_summary),
            "diagnostic_summary_sha256": sha256_file(args.diagnostic_summary),
            "include_cell": sorted(args.include_cell),
            "exclude_cell": sorted(args.exclude_cell),
            "selected_cell_ids": [cell["id"] for cell in selected],
        }
        identity = stage_identity("abba", provenance, options)
        run_dir, identity_digest = prepare_run_directory(args.run_dir, identity)
        started = time.monotonic()
        results = []
        resumed_count = 0
        for index, cell in enumerate(selected, 1):
            result, resumed = abba_cell(
                cell, pair, args, identity_digest, run_dir)
            results.append(result)
            resumed_count += int(resumed)
            print("ABBA %d/%d %s%s" % (
                index, len(selected), cell["id"],
                " resumed" if resumed else ""), flush=True)
        source_still_frozen(provenance)
        rows = sorted(results, key=lambda value: value["cell"]["id"])
        summary = {
            "schema": SCHEMA, "status": "complete", "stage": "abba",
            "authoritative": True, "identity_sha256": identity_digest,
            "matrix_sha256": matrix_digest(matrix),
            "cell_count": len(rows), "resumed_cell_count": resumed_count,
            "accepted_child_invocations": len(rows) * args.rounds * 4,
            "elapsed_seconds": time.monotonic() - started,
            "candidate_binary_sha256": provenance["candidate"]["sha256"],
            "baseline_binary_sha256": provenance["baseline"]["sha256"],
            "rows": [{
                "id": row["cell"]["id"], "cell": row["cell"],
                "metrics": row["metrics"], "cpu_pair": row["cpu_pair"],
                "selected_decode_path": row["selected_decode_path"],
                "selected_decode_rule": row["selected_decode_rule"],
                "workload_digests": row["workload_digests"],
                "missing_original_indices": row["missing_original_indices"],
            } for row in rows],
        }
        atomic_json(run_dir / "summary.json", summary)
        print(json.dumps({
            "stage": "abba", "run_dir": str(run_dir),
            "cells": summary["cell_count"],
            "elapsed_seconds": summary["elapsed_seconds"],
            "summary_sha256": sha256_file(run_dir / "summary.json"),
        }, sort_keys=True))


def make_dry_run_manifest(matrix, diagnostic_summary_path=None, near_ratio=1.10):
    mandatory_tags = {"loss_sweep", "reuse_sweep", "batch_sweep",
                      "large_tiling", "selector_isolation",
                      "r1_selector_fallback", "r1_selector_promoted",
                      "r1_large_recovery", "direct_odd_arity"}
    mandatory = [cell for cell in matrix["cells"]
                 if set(cell["tags"]) & mandatory_tags]
    selected = None
    if diagnostic_summary_path:
        diagnostic = load_diagnostic_summary(diagnostic_summary_path, matrix)
        selected = selected_abba_cells(
            matrix, diagnostic, near_ratio, [], [])
    pairs = physical_cpu_pairs()
    ordered_cells = sorted(
        matrix["cells"],
        key=lambda cell: (cell["estimated_peak_bytes"], cell["id"]))
    assignments = [{"cell_id": cell["id"],
                    "cpu_pair": list(pairs[index % len(pairs)])}
                   for index, cell in enumerate(ordered_cells)]
    return {
        "schema": SCHEMA, "kind": "dry_run_manifest",
        "matrix_sha256": matrix_digest(matrix),
        "matrix_cell_count": len(matrix["cells"]),
        "diagnostic_child_invocations": len(matrix["cells"]) * 2,
        "detected_physical_cpu_pairs": [list(pair) for pair in pairs],
        "detected_pair_count": len(pairs),
        "deterministic_diagnostic_assignments": assignments,
        "default_memory_budget_bytes": min(mem_total_bytes() // 4, 32 << 30),
        "mandatory_abba_cell_count": len(mandatory),
        "mandatory_abba_child_invocations_three_rounds": len(mandatory) * 12,
        "selected_abba_cell_count": None if selected is None else len(selected),
        "selected_abba_child_invocations_three_rounds": (
            None if selected is None else len(selected) * 12),
        "maximum_abba_cell_count": len(matrix["cells"]),
        "maximum_abba_child_invocations_three_rounds":
            len(matrix["cells"]) * 12,
        "sizes": matrix["sizes"], "transform_sides": matrix["transform_sides"],
        "losses": matrix["losses"], "reuse_counts": matrix["reuse_counts"],
        "batch_counts": matrix["batch_counts"],
        "selector_isolation_cells": [
            cell["id"] for cell in matrix["cells"]
            if "selector_isolation" in cell["tags"]],
    }


def run_merge(args):
    diagnostic = read_json(args.diagnostic_summary)
    abba = read_json(args.abba_summary)
    if (diagnostic.get("schema") != SCHEMA or abba.get("schema") != SCHEMA or
            diagnostic.get("status") != "complete" or
            abba.get("status") != "complete" or
            diagnostic.get("stage") != "diagnostic" or
            abba.get("stage") != "abba" or
            diagnostic.get("authoritative") is not False or
            abba.get("authoritative") is not True):
        raise RuntimeError("merge inputs have wrong stages")
    if (diagnostic.get("cell_count") != len(diagnostic.get("rows", [])) or
            abba.get("cell_count") != len(abba.get("rows", []))):
        raise RuntimeError("merge input row count mismatch")
    if diagnostic.get("matrix_sha256") != abba.get("matrix_sha256"):
        raise RuntimeError("merge input matrix mismatch")
    if diagnostic.get("baseline_binary_sha256") != \
            abba.get("baseline_binary_sha256") or \
            diagnostic.get("candidate_binary_sha256") != \
            abba.get("candidate_binary_sha256"):
        raise RuntimeError("merge input binary mismatch")
    bundle = {
        "schema": SCHEMA, "kind": "merged_evidence_bundle",
        "diagnostic_authoritative": False,
        "abba_authoritative": True,
        "diagnostic_summary": diagnostic,
        "abba_summary": abba,
        "input_sha256": {
            "diagnostic": sha256_file(args.diagnostic_summary),
            "abba": sha256_file(args.abba_summary),
        },
    }
    atomic_json(args.output, bundle)
    print(json.dumps({"output": os.path.realpath(args.output),
                      "sha256": sha256_file(args.output)}, sort_keys=True))


def self_test():
    matrix = make_matrix()
    validate_matrix(matrix)
    assert parse_cpu_list("0-2,5,8-9") == {0, 1, 2, 5, 8, 9}
    selectors = {cell["id"] for cell in matrix["cells"]
                 if "selector_isolation" in cell["tags"]}
    expected_selectors = {
        cell_key(k, 32, size, 32, 8, 1)
        for k in (33, 34, 35) for size in (16384, 65536)}
    assert selectors == expected_selectors
    r1_boundaries = {
        (cell["K"], cell["bytes"], tuple(cell["tags"]))
        for cell in matrix["cells"]
        if "r1_selector_fallback" in cell["tags"] or
           "r1_selector_promoted" in cell["tags"]
    }
    expected_r1_boundaries = set()
    for lower_k, promoted_k, byte_count in R1_SELECTOR_BOUNDARIES:
        expected_r1_boundaries.add(
            (lower_k, byte_count, ("r1_selector_fallback",)))
        expected_r1_boundaries.add(
            (promoted_k, byte_count, ("r1_selector_promoted",)))
    assert r1_boundaries == expected_r1_boundaries
    direct_odd = {
        (cell["bytes"], cell["loss"])
        for cell in matrix["cells"]
        if "direct_odd_arity" in cell["tags"]
    }
    assert direct_odd == {
        (byte_count, loss)
        for byte_count in (2048, 65536, 1048576)
        for loss in DIRECT_ODD_LOSSES
    }
    calls = []
    for kind, value in zip(ABBA_ORDER * 3,
                           (10, 5, 5, 10) * 3):
        calls.append({"kind": kind, "times_us": {
            "encode": value, "decode_first": value,
            "decode_reuse": value}})
    metrics = abba_metrics(calls, 3)
    for value in metrics.values():
        assert abs(value["ratio"] - 2.0) < 1e-12
    encoded = canonical_bytes(matrix)
    assert json.loads(encoded) == matrix
    assert matrix_digest(matrix) == hashlib.sha256(encoded).hexdigest()
    print(json.dumps({
        "self_test": "passed", "matrix_cells": len(matrix["cells"]),
        "matrix_sha256": matrix_digest(matrix),
        "selector_isolation_cells": len(selectors),
        "r1_selector_boundary_cells": len(r1_boundaries),
        "direct_odd_arity_cells": len(direct_odd),
    }, sort_keys=True))


def add_benchmark_arguments(parser):
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--candidate-commit", required=True)
    parser.add_argument("--candidate-source", required=True)
    parser.add_argument("--candidate-build", required=True)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument(
        "--lock", default="/tmp/leopard-gf8-authoritative.lock",
        help="canonical lock shared by diagnostics and authoritative runs")
    parser.add_argument("--maximum-attempts", type=int, default=3)


def validate_benchmark_arguments(args):
    if args.maximum_attempts < 1:
        raise RuntimeError("--maximum-attempts must be positive")
    if not COMMIT_PATTERN.match(args.candidate_commit):
        raise RuntimeError("--candidate-commit must be a full lowercase SHA")


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    make = subparsers.add_parser("make-matrix")
    make.add_argument("--output", required=True)
    dry = subparsers.add_parser("dry-run")
    dry.add_argument("--matrix", required=True)
    dry.add_argument("--output", required=True)
    dry.add_argument("--diagnostic-summary")
    dry.add_argument("--near-ratio", type=float, default=1.10)
    subparsers.add_parser("self-test")
    diagnostic = subparsers.add_parser("diagnostic")
    add_benchmark_arguments(diagnostic)
    diagnostic.add_argument(
        "--cpu-pairs", help="optional comma-separated cpu:sibling pairs")
    diagnostic.add_argument("--maximum-workers", type=int, default=0)
    diagnostic.add_argument("--memory-budget-gib", type=int, default=0)
    abba = subparsers.add_parser("abba")
    add_benchmark_arguments(abba)
    abba.add_argument("--diagnostic-summary", required=True)
    abba.add_argument("--near-ratio", type=float, default=1.10)
    abba.add_argument("--rounds", type=int, default=3)
    abba.add_argument("--cpu", type=int)
    abba.add_argument("--sibling", type=int)
    abba.add_argument("--include-cell", action="append", default=[])
    abba.add_argument("--exclude-cell", action="append", default=[])
    merge = subparsers.add_parser("merge")
    merge.add_argument("--diagnostic-summary", required=True)
    merge.add_argument("--abba-summary", required=True)
    merge.add_argument("--output", required=True)
    return parser.parse_args()


def main():
    args = parse_arguments()
    if args.command == "self-test":
        self_test()
    elif args.command == "make-matrix":
        matrix = make_matrix()
        validate_matrix(matrix)
        atomic_json(args.output, matrix)
        print(json.dumps({"output": os.path.realpath(args.output),
                          "cells": len(matrix["cells"]),
                          "sha256": matrix_digest(matrix)}, sort_keys=True))
    elif args.command == "dry-run":
        matrix = read_json(args.matrix)
        validate_matrix(matrix)
        manifest = make_dry_run_manifest(
            matrix, args.diagnostic_summary, args.near_ratio)
        atomic_json(args.output, manifest)
        print(json.dumps({"output": os.path.realpath(args.output),
                          "sha256": sha256_file(args.output)}, sort_keys=True))
    elif args.command == "diagnostic":
        validate_benchmark_arguments(args)
        if args.maximum_workers < 0 or args.memory_budget_gib < 0:
            raise RuntimeError("worker and memory limits cannot be negative")
        run_diagnostic(args)
    elif args.command == "abba":
        validate_benchmark_arguments(args)
        if args.rounds < 2:
            raise RuntimeError("ABBA requires at least two rounds")
        if args.near_ratio <= 0:
            raise RuntimeError("--near-ratio must be positive")
        run_abba(args)
    elif args.command == "merge":
        run_merge(args)
    else:
        raise RuntimeError("unhandled command")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print("error: " + str(error), file=sys.stderr)
        sys.exit(1)
