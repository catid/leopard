#!/usr/bin/env python3
"""Broad diagnostic Leopard-main versus Leopard2 all-K gap map.

This is intentionally not promotion evidence.  It saturates all allowed CPUs
with independent single-thread cells to find algorithmic regions worth fixing.
Each cell uses an ABBA process order and retains exact JSON from independently
linked exact-main and Leopard2 executables.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import contextlib
import ctypes
import dataclasses
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import stat
import statistics
import subprocess
import sys
import time
from typing import Any, Mapping, Sequence


MAIN_COMMIT = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
PURE_AVX2_MAIN_SHA256 = (
    "a43d7f43ff2e887ebcd47a1e94f806847a5d8b858a4e383e6c8d5e528a7dd910")
ORDER = ("main", "leopard2", "leopard2", "main")
CHILD_ENV = {
    "LANG": "C", "LC_ALL": "C", "OMP_DYNAMIC": "FALSE",
    "OMP_NUM_THREADS": "1", "OMP_PROC_BIND": "TRUE",
    "OMP_PLACES": "cores", "PATH": "/usr/bin:/bin", "TZ": "UTC",
}
LINUX_F_ADD_SEALS = getattr(fcntl, "F_ADD_SEALS", 1033)
LINUX_F_GET_SEALS = getattr(fcntl, "F_GET_SEALS", 1034)
LINUX_F_SEAL_SEAL = getattr(fcntl, "F_SEAL_SEAL", 0x0001)
LINUX_F_SEAL_SHRINK = getattr(fcntl, "F_SEAL_SHRINK", 0x0002)
LINUX_F_SEAL_GROW = getattr(fcntl, "F_SEAL_GROW", 0x0004)
LINUX_F_SEAL_WRITE = getattr(fcntl, "F_SEAL_WRITE", 0x0008)
LINUX_MFD_CLOEXEC = getattr(os, "MFD_CLOEXEC", 0x0001)
LINUX_MFD_ALLOW_SEALING = getattr(os, "MFD_ALLOW_SEALING", 0x0002)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def linux_memfd_create(name: str) -> int:
    flags = LINUX_MFD_CLOEXEC | LINUX_MFD_ALLOW_SEALING
    if hasattr(os, "memfd_create"):
        return os.memfd_create(name, flags)
    libc = ctypes.CDLL(None, use_errno=True)
    creator = getattr(libc, "memfd_create", None)
    require(creator is not None,
            "all-K executable snapshots require Linux memfd_create support")
    creator.argtypes = (ctypes.c_char_p, ctypes.c_uint)
    creator.restype = ctypes.c_int
    descriptor = creator(name.encode("utf-8"), flags)
    if descriptor < 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error), name)
    return descriptor


def require_hex(value: str, label: str) -> str:
    require(isinstance(value, str) and len(value) == 40 and
            all(character in "0123456789abcdef" for character in value),
            f"{label} must be exactly 40 lowercase hexadecimal characters")
    return value


def require_sha256(value: str, label: str) -> str:
    require(isinstance(value, str) and len(value) == 64 and
            all(character in "0123456789abcdef" for character in value),
            f"{label} must be exactly 64 lowercase hexadecimal characters")
    return value


def inspect_isa_disassembly(disassembly: str) -> dict[str, int]:
    evex = 0
    ymm = 0
    for line in disassembly.splitlines():
        fields = line.lstrip().split()
        if fields and fields[0].endswith(":") and len(fields) > 1 and \
                fields[1] == "62":
            evex += 1
        if re.search(r"\bymm[0-9]+\b", line):
            ymm += 1
    return {"evex_prefixed_instruction_count": evex,
            "ymm_operand_instruction_count": ymm}


def inspect_pure_avx2(path: Path, label: str) -> dict[str, Any]:
    completed = subprocess.run(
        ["/usr/bin/objdump", "-d", "-M", "intel", str(path)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        timeout=180.0, check=False)
    require(completed.returncode == 0,
            f"objdump failed for {label}: {completed.stderr.strip()}")
    result = inspect_isa_disassembly(completed.stdout)
    require(result["evex_prefixed_instruction_count"] == 0,
            f"{label} contains EVEX-prefixed instructions")
    require(result["ymm_operand_instruction_count"] > 0,
            f"{label} contains no visible AVX2 YMM instructions")
    return {
        **result,
        "objdump_version": subprocess.check_output(
            ["/usr/bin/objdump", "--version"], text=True).splitlines()[0],
    }


@contextlib.contextmanager
def benchmark_lock(path: Path):
    resolved = path.resolve()
    resolved.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(resolved, os.O_CREAT | os.O_RDWR, 0o600)
    try:
        print(f"waiting for exclusive benchmark lock {resolved}", flush=True)
        fcntl.flock(descriptor, fcntl.LOCK_EX)
        print("acquired benchmark lock", flush=True)
        yield str(resolved)
    finally:
        fcntl.flock(descriptor, fcntl.LOCK_UN)
        os.close(descriptor)


def git_output(root: Path, *arguments: str) -> str:
    completed = subprocess.run(
        ["/usr/bin/git", "-C", str(root), *arguments],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        env={
            "GIT_NO_REPLACE_OBJECTS": "1",
            "LANG": "C", "LC_ALL": "C", "PATH": "/usr/bin:/bin",
        },
        check=False)
    require(completed.returncode == 0,
            f"git {' '.join(arguments)} failed for {root}: "
            f"{completed.stderr.strip()}")
    return completed.stdout.strip()


def git_identity(source_root: Path, requested_commit: str) -> dict[str, Any]:
    root = source_root.resolve(strict=True)
    require(root.is_dir(), f"current source root is not a directory: {root}")
    requested = require_hex(requested_commit, "current source commit")
    top = Path(git_output(root, "rev-parse", "--show-toplevel")).resolve(
        strict=True)
    require(top == root,
            f"current source root is not the Git top level: {root} != {top}")
    head = require_hex(git_output(root, "rev-parse", "HEAD"),
                       "current source HEAD")
    require(head == requested,
            f"current source HEAD mismatch: requested {requested}, got {head}")
    tree = require_hex(git_output(root, "rev-parse", "HEAD^{tree}"),
                       "current source tree")
    for flag in ("-v", "-f"):
        index_records = [record for record in
                         git_output(root, "ls-files", flag, "-z").split("\0")
                         if record]
        require(index_records and
                all(record.startswith("H ") for record in index_records),
                "current source index uses assume-unchanged, skip-worktree, "
                "fsmonitor-valid, or another non-default flag")
    status = git_output(root, "status", "--porcelain=v1",
                        "--untracked-files=normal")
    require(not status,
            "current source has tracked or untracked modifications; "
            "diagnostic identity requires a clean committed tree")
    return {
        "path": str(root), "head": head, "tree": tree,
        "tracked_status": "clean",
    }


def file_identity(path: Path, label: str) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
    descriptor = os.open(resolved, flags)
    try:
        before = os.fstat(descriptor)
        require(stat.S_ISREG(before.st_mode), f"{label} is not a regular file")
        digest = hashlib.sha256()
        while True:
            chunk = os.read(descriptor, 1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    stable_fields = ("st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns",
                     "st_ctime_ns")
    require(all(getattr(before, name) == getattr(after, name)
                for name in stable_fields),
            f"{label} changed while it was hashed")
    path_status = resolved.stat()
    require(all(getattr(after, name) == getattr(path_status, name)
                for name in stable_fields),
            f"{label} path changed while it was hashed")
    require(os.access(resolved, os.X_OK), f"{label} is not executable")
    return {
        "path": str(resolved), "sha256": digest.hexdigest(),
        "device": after.st_dev, "inode": after.st_ino,
        "mode": after.st_mode, "size": after.st_size,
        "mtime_ns": after.st_mtime_ns, "ctime_ns": after.st_ctime_ns,
    }


def sealed_snapshot_identity(descriptor: int, label: str) -> dict[str, Any]:
    status = os.fstat(descriptor)
    require(stat.S_ISREG(status.st_mode), f"{label} snapshot is not regular")
    digest = hashlib.sha256()
    offset = 0
    while offset < status.st_size:
        chunk = os.pread(descriptor, min(1024 * 1024, status.st_size - offset),
                         offset)
        require(bool(chunk), f"{label} snapshot ended before its recorded size")
        digest.update(chunk)
        offset += len(chunk)
    seals = fcntl.fcntl(descriptor, LINUX_F_GET_SEALS)
    required_seals = (LINUX_F_SEAL_SEAL | LINUX_F_SEAL_SHRINK |
                      LINUX_F_SEAL_GROW | LINUX_F_SEAL_WRITE)
    require(seals & required_seals == required_seals,
            f"{label} snapshot is not immutably sealed")
    return {
        "kind": "linux-sealed-memfd-v1",
        "sha256": digest.hexdigest(), "size": status.st_size,
        "mode": status.st_mode, "seals": seals,
    }


def snapshot_executable(path: Path, label: str) \
        -> tuple[dict[str, Any], int, dict[str, Any]]:
    require(sys.platform.startswith("linux") and hasattr(os, "pread"),
            "all-K executable snapshots require Linux sealed memfd support")
    source_identity = file_identity(path, label)
    source = os.open(source_identity["path"],
                     os.O_RDONLY | getattr(os, "O_CLOEXEC", 0))
    snapshot = -1
    try:
        snapshot = linux_memfd_create(
            "leopard2-allk-" + label.replace(" ", "-"))
        digest = hashlib.sha256()
        while True:
            chunk = os.read(source, 1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
            view = memoryview(chunk)
            while view:
                written = os.write(snapshot, view)
                require(written > 0, f"{label} snapshot write made no progress")
                view = view[written:]
        require(digest.hexdigest() == source_identity["sha256"],
                f"{label} changed between identity capture and snapshot copy")
        require(os.pread(snapshot, 4, 0) == b"\x7fELF",
                f"{label} is not an ELF executable")
        require(file_identity(Path(source_identity["path"]), label) ==
                source_identity,
                f"{label} changed while its executable snapshot was created")
        os.fchmod(snapshot, 0o500)
        required_seals = (LINUX_F_SEAL_SEAL | LINUX_F_SEAL_SHRINK |
                          LINUX_F_SEAL_GROW | LINUX_F_SEAL_WRITE)
        fcntl.fcntl(snapshot, LINUX_F_ADD_SEALS, required_seals)
        snapshot_identity = sealed_snapshot_identity(snapshot, label)
        require(snapshot_identity["sha256"] == source_identity["sha256"] and
                snapshot_identity["size"] == source_identity["size"],
                f"{label} sealed snapshot does not match its source identity")
        return source_identity, snapshot, snapshot_identity
    except BaseException:
        if snapshot >= 0:
            os.close(snapshot)
        raise
    finally:
        os.close(source)


def canonical_digest(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def write_json_atomic(path: Path, value: Any, *, pretty: bool = False) -> None:
    temporary = path.with_name(path.name + ".tmp")
    if pretty:
        encoded = json.dumps(value, indent=2, sort_keys=True) + "\n"
    else:
        encoded = json.dumps(value, sort_keys=True) + "\n"
    temporary.parent.mkdir(parents=True, exist_ok=True)
    temporary.write_text(encoded, encoding="utf-8")
    os.replace(temporary, path)


def validate_manifest(document: Mapping[str, Any], contract: Mapping[str, Any],
                      contract_digest: str,
                      cells: Sequence["Cell"]) -> None:
    require(set(document) == {
        "schema", "run_contract", "run_contract_sha256", "cells", "completion"
    }, "all-K manifest keys changed")
    require(document.get("schema") == "leopard2-all-k-gap-manifest/v2",
            "all-K manifest schema mismatch")
    require(document.get("run_contract") == contract,
            "existing all-K manifest run contract does not match this request")
    require(document.get("run_contract_sha256") == contract_digest,
            "existing all-K manifest contract digest mismatch")
    require(canonical_digest(document.get("run_contract")) == contract_digest,
            "existing all-K manifest contract bytes are internally inconsistent")
    require(document.get("cells") == [dataclasses.asdict(cell) for cell in cells],
            "existing all-K manifest cell matrix mismatch")
    completion = document.get("completion")
    require(completion is None or isinstance(completion, dict),
            "all-K manifest completion must be null or an object")


@dataclasses.dataclass(frozen=True)
class Cell:
    identifier: str
    family: str
    k: int
    r: int
    shard_bytes: int
    losses: int
    redundancy_band: str
    loss_band: str
    seed: int
    iterations: int
    reuse: int
    warmup: int


def ceil_pow2(value: int) -> int:
    return 1 << (value - 1).bit_length()


def parent_for(k: int, r: int) -> tuple[int, int]:
    padded = ceil_pow2(r)
    return ceil_pow2(k + padded), padded


def gf8_r_values(k: int) -> list[tuple[str, int]]:
    max_padded = 1 << ((256 - k).bit_length() - 1)
    maximum = min(k, max_padded)
    values = (("low-R", 1), ("mid-R", max(1, (maximum + 1) // 2)),
              ("max-GF8-R", maximum))
    result: list[tuple[str, int]] = []
    seen: set[int] = set()
    for name, value in values:
        if value not in seen:
            seen.add(value)
            result.append((name, value))
    return result


GF16_K = (
    129, 130, 191, 192, 193, 200, 224, 225, 239, 240, 241, 248,
    249, 255, 256, 257, 300, 511, 512, 513, 1000, 2048, 4096,
)


def gf16_r_values(k: int) -> list[tuple[str, int]]:
    if k <= 255:
        gap = 256 - k
        first_forced = (1 << (gap.bit_length() - 1)) + 1
        values = (("first-GF16-R", first_forced),
                  ("mid-R", min(k, max(first_forced, (k + 3) // 4))),
                  ("high-R", min(k, 512)))
    else:
        values = (("low-R", 1), ("mid-R", min(k, max(2, k // 8))),
                  ("high-R", min(k, 512)))
    result: list[tuple[str, int]] = []
    seen: set[int] = set()
    for name, value in values:
        parent, _ = parent_for(k, value)
        if value <= k and parent > 256 and value not in seen:
            seen.add(value)
            result.append((name, value))
    return result


def make_cells(*, gf8_only: bool = False) -> list[Cell]:
    cells: list[Cell] = []
    seed = 0x38000000
    for k in range(1, 256):
        for redundancy_band, r in gf8_r_values(k):
            parent, _ = parent_for(k, r)
            require(parent <= 256,
                    "GF8 all-K matrix generated an out-of-field parent")
            for shard_bytes in (4096, 65536):
                for loss_band, losses in (("one-loss", 1), ("max-loss", r)):
                    if loss_band == "max-loss" and losses == 1:
                        continue
                    seed += 1
                    cells.append(Cell(
                        f"gf8-k{k}-r{r}-b{shard_bytes}-l{losses}", "gf8-all-k",
                        k, r, shard_bytes, losses, redundancy_band, loss_band,
                        seed, 5, 16 if shard_bytes == 4096 else 4, 1))
    if not gf8_only:
        for k in GF16_K:
            for redundancy_band, r in gf16_r_values(k):
                parent, _ = parent_for(k, r)
                require(256 < parent <= 65536,
                        "GF16 representative matrix generated an invalid parent")
                for shard_bytes in (512, 4096):
                    for loss_band, losses in (("one-loss", 1),
                                              ("max-loss", r)):
                        if loss_band == "max-loss" and losses == 1:
                            continue
                        seed += 1
                        cells.append(Cell(
                            f"gf16-k{k}-r{r}-b{shard_bytes}-l{losses}",
                            "gf16-representative", k, r, shard_bytes, losses,
                            redundancy_band, loss_band, seed, 5,
                            16 if shard_bytes == 512 else 8, 1))
    require(len({cell.identifier for cell in cells}) == len(cells),
            "all-K cell identifiers are not unique")
    return cells


def direct_locator_cutoff(field: str, parent: int) -> int:
    if field == "gf8":
        if parent <= 8: return parent
        if parent == 16: return 8
        if parent in (32, 128): return 9
        if parent == 64: return 8
        return 7
    if parent <= 32: return parent
    if parent == 64: return 34
    if parent == 128: return 24
    if parent in (256, 512): return 16
    return 14


def classify_paths(cell: Cell, result: Mapping[str, Any]) -> dict[str, Any]:
    resolved = result["resolved"]
    profile = resolved["profile"]
    field = resolved["field"]
    backend = resolved["backend"]
    parent = int(resolved["parent_count"])
    padded = int(resolved["padded_side"])
    require(profile == "legacy_high_v1",
            "all-K current benchmark selected a non-legacy-high profile")
    expected_field = "gf8" if cell.family == "gf8-all-k" else "gf16"
    require(field == expected_field,
            "all-K current benchmark selected the wrong finite field")
    require(backend == "avx2",
            "all-K current benchmark did not select explicit AVX2")
    decode = resolved.get("selected_decode_path")
    decode_rule = resolved.get("selected_decode_rule")
    require(decode in ("direct", "generic", "materialized", "tiled"),
            "all-K current benchmark did not report its selected decode path")
    require(isinstance(decode_rule, str) and decode_rule != "unknown",
            "all-K current benchmark did not report its decode rule")
    if padded == 1:
        encode = "direct-xor-single-parity"
    else:
        encode = "specialized-high-transform"
    if decode == "direct":
        workspace = "direct-bounded"
        locator = "none"
    elif decode == "generic":
        workspace = "materialized-N"
        locator = "active-parent-" + (
            "sparse-direct" if padded <= direct_locator_cutoff(field, parent)
            else "walsh")
    else:
        workspace = ("materialized-N" if decode == "materialized" else
                     "tiled-side-plus-losses")
        permanent_cached = padded > cell.r and \
            cell.r <= direct_locator_cutoff(field, parent)
        effective = cell.r if permanent_cached else padded
        locator = "active-parent-" + (
            "sparse-direct" if effective <= direct_locator_cutoff(field, parent)
            else "walsh")
    return {
        "profile": profile, "field": field, "backend": backend,
        "parent_count": parent, "padded_side": padded,
        "parent_inflation": parent / float(cell.k + cell.r),
        "encode_path": encode, "decode_path": decode,
        "decode_rule": decode_rule,
        "decode_workspace": workspace, "locator_setup": locator,
        "decode_required_work_slots":
            int(resolved.get("decode_required_work_slots", 0)),
    }


def command(role: str, executable: Path, cell: Cell, cpu: int,
            with_current_legacy: bool) -> list[str]:
    common = [
        "/usr/bin/taskset", "-c", str(cpu), str(executable),
        "--k", str(cell.k), "--r", str(cell.r),
        "--bytes", str(cell.shard_bytes), "--loss", str(cell.losses),
        "--batch", "1", "--reuse", str(cell.reuse),
        "--iterations", str(cell.iterations), "--warmup", str(cell.warmup),
        "--threads", "1", "--seed", str(cell.seed),
    ]
    if role == "leopard2":
        field = "gf8" if cell.family == "gf8-all-k" else "gf16"
        common.extend(("--profile", "high", "--field", field,
                       "--backend", "avx2", "--retain-samples",
                       "--report-decode-path"))
        if not with_current_legacy:
            common.append("--skip-legacy")
    common.extend(("--json", "-"))
    return common


def source_attestation_command(executable: Path, cpu: int) -> list[str]:
    return [
        "/usr/bin/taskset", "-c", str(cpu), str(executable),
        "--k", "3", "--r", "2", "--bytes", "64", "--loss", "1",
        "--batch", "1", "--reuse", "1", "--iterations", "1",
        "--warmup", "0", "--threads", "1", "--seed", "7",
        "--profile", "high", "--field", "gf8", "--backend", "avx2",
        "--skip-legacy", "--retain-samples", "--attest-source",
        "--json", "-",
    ]


def run_one(role: str, command_value: Sequence[str], timeout: float,
            snapshot_descriptor: int,
            snapshot_identity: Mapping[str, Any]) -> dict[str, Any]:
    execution_command = list(command_value)
    require(len(execution_command) >= 4 and
            execution_command[0:2] == ["/usr/bin/taskset", "-c"],
            "all-K execution command has an unexpected launcher shape")
    execution_command[3] = f"/proc/self/fd/{snapshot_descriptor}"
    started = time.time_ns()
    completed = subprocess.run(
        execution_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        env=CHILD_ENV, timeout=timeout, pass_fds=(snapshot_descriptor,))
    finished = time.time_ns()
    record: dict[str, Any] = {
        "role": role, "command": list(command_value),
        "executable_snapshot_sha256": snapshot_identity["sha256"],
        "returncode": completed.returncode,
        "duration_ns": finished - started,
        "stdout_sha256": hashlib.sha256(completed.stdout).hexdigest(),
        "stderr": completed.stderr.decode("utf-8", errors="replace"),
    }
    if completed.returncode == 0:
        record["result"] = json.loads(completed.stdout.decode("utf-8"))
    else:
        record["stdout"] = completed.stdout.decode("utf-8", errors="replace")
    return record


def metric(record: Mapping[str, Any], path: Sequence[str]) -> float:
    value: Any = record["result"]
    for key in path:
        value = value[key]
    return float(value)


def geometric(values: Sequence[float]) -> float:
    require(bool(values) and
            all(value > 0 and math.isfinite(value) for value in values),
            "benchmark timing samples must be finite and positive")
    return math.exp(statistics.fmean(math.log(value) for value in values))


def validate_invocation_identities(
        invocations: Sequence[Mapping[str, Any]], main_commit: str,
        current_source: Mapping[str, Any],
        main_snapshot: Mapping[str, Any],
        current_snapshot: Mapping[str, Any]) -> None:
    require(len(invocations) == len(ORDER),
            "a valid all-K cell must contain exactly four invocations")
    require(tuple(record.get("role") for record in invocations) == ORDER,
            "all-K invocation role order mismatch")
    for index, record in enumerate(invocations):
        require(record.get("returncode") == 0,
                f"all-K invocation {index} did not exit successfully")
        result = record.get("result")
        require(isinstance(result, dict),
                f"all-K invocation {index} has no JSON result")
        build = result.get("build")
        require(isinstance(build, dict),
                f"all-K invocation {index} has no build identity")
        if record["role"] == "main":
            expected_snapshot = main_snapshot
            require(build.get("main_source_commit") == main_commit,
                    "exact-main benchmark embedded commit mismatch")
        else:
            expected_snapshot = current_snapshot
            require(result.get("schema") == "leopard2-benchmark-v3",
                    "Leopard2 benchmark decode-path schema mismatch")
            require(result.get("parameters", {}).get(
                "report_decode_path") is True,
                "Leopard2 benchmark did not record decode-path reporting")
        require(record.get("executable_snapshot_sha256") ==
                expected_snapshot.get("sha256"),
                f"all-K invocation {index} executable snapshot mismatch")


def validate_source_attestation(
        record: Mapping[str, Any], current_source: Mapping[str, Any],
        current_snapshot: Mapping[str, Any]) -> None:
    require(record.get("returncode") == 0,
            "Leopard2 source-attestation preflight failed")
    require(record.get("executable_snapshot_sha256") ==
            current_snapshot.get("sha256"),
            "Leopard2 source-attestation snapshot mismatch")
    result = record.get("result")
    require(isinstance(result, dict) and
            result.get("schema") == "leopard2-benchmark-v5",
            "Leopard2 benchmark source-attestation schema mismatch")
    require(result.get("parameters", {}).get("attest_source") is True,
            "Leopard2 benchmark did not record source attestation")
    build = result.get("build")
    require(isinstance(build, dict),
            "Leopard2 source attestation has no build identity")
    require(build.get("source_commit") == current_source.get("head"),
            "Leopard2 benchmark embedded commit mismatch")
    require(build.get("source_tree") == current_source.get("tree"),
            "Leopard2 benchmark embedded tree mismatch")
    require(build.get("source_tracked_dirty") is False,
            "Leopard2 benchmark was built from a tracked-dirty tree")
    resolved = result.get("resolved", {})
    require((resolved.get("profile"), resolved.get("field"),
             resolved.get("backend")) ==
            ("legacy_high_v1", "gf8", "avx2"),
            "Leopard2 source attestation resolved the wrong codec identity")
    require(result.get("correctness", {}).get(
        "leopard2_round_trip") is True,
        "Leopard2 source-attestation round trip failed")


def validate_correctness(cell: Cell,
                         invocations: Sequence[Mapping[str, Any]]) -> None:
    expected_parameters = {
        "K": cell.k, "R": cell.r, "shard_bytes": cell.shard_bytes,
        "loss_count": cell.losses, "batch": 1, "reuse": cell.reuse,
        "iterations": cell.iterations, "warmup": cell.warmup,
        "thread_count": 1, "seed": cell.seed,
    }
    fingerprints: list[Mapping[str, Any]] = []
    for index, record in enumerate(invocations):
        result = record["result"]
        parameters = result.get("parameters")
        require(isinstance(parameters, dict),
                f"all-K invocation {index} has no parameters")
        for name, value in expected_parameters.items():
            require(parameters.get(name) == value,
                    f"all-K invocation {index} parameter {name} mismatch")
        correctness = result.get("correctness")
        require(isinstance(correctness, dict),
                f"all-K invocation {index} has no correctness record")
        if record["role"] == "main":
            require(correctness.get("round_trip") is True,
                    "exact-main benchmark round trip failed")
        else:
            require(correctness.get("leopard2_round_trip") is True,
                    "Leopard2 benchmark round trip failed")
        fingerprint = result.get("workload_digests")
        expected_digest_keys = {
            "algorithm", "original_data", "transmitted_parity",
            "recovered_originals",
        }
        require(isinstance(fingerprint, dict) and
                set(fingerprint) == expected_digest_keys,
                f"all-K invocation {index} workload digest structure mismatch")
        require(fingerprint["algorithm"] == "fnv1a64",
                f"all-K invocation {index} workload digest algorithm mismatch")
        for name in ("original_data", "transmitted_parity",
                     "recovered_originals"):
            value = fingerprint[name]
            require(isinstance(value, str) and len(value) == 16 and
                    all(character in "0123456789abcdef"
                        for character in value),
                    f"all-K invocation {index} workload digest {name} "
                    "is not lowercase FNV-1a hex")
        fingerprints.append(fingerprint)
    require(all(fingerprint == fingerprints[0] for fingerprint in fingerprints),
            "Leopard1 and Leopard2 workload digests differ")


def gap_tags(cell: Cell, paths: Mapping[str, Any], encode_speedup: float,
             decode_first_speedup: float, decode_reuse_speedup: float) -> list[str]:
    tags: list[str] = []
    if encode_speedup < 1.05:
        if paths["encode_path"] == "direct-xor-single-parity":
            tags.append("encode:R1-xor/API-overhead")
        else:
            tags.append("encode:legacy-high-transform-not-faster")
        if cell.k % paths["padded_side"] != 0:
            tags.append("encode:partial-final-message-block")
    if decode_first_speedup < 1.05 or decode_reuse_speedup < 1.05:
        decode_path = paths["decode_path"]
        if decode_path == "direct":
            tags.append("decode:direct-path-overhead-or-kernel")
        elif decode_path == "generic":
            tags.append("decode:generic-fallback-crossover")
        else:
            tags.append("decode:specialized-high-crossover")
        if paths["decode_workspace"] == "materialized-N":
            tags.append("decode:materialized-N-workspace")
        if paths["locator_setup"].endswith("walsh"):
            tags.append("decode:active-parent-Walsh-setup")
    if paths["parent_inflation"] >= 1.5:
        tags.append("common:dyadic-parent-inflation>=1.5x")
    if cell.shard_bytes <= 4096:
        tags.append("common:small-shard-fixed-cost")
    return sorted(set(tags))


def analyze_cell(cell: Cell, invocations: Sequence[Mapping[str, Any]], cpu: int,
                 contract_digest: str, main_commit: str,
                 current_source: Mapping[str, Any],
                 main_snapshot: Mapping[str, Any],
                 current_snapshot: Mapping[str, Any]) -> dict[str, Any]:
    failures = [record for record in invocations if record["returncode"] != 0]
    result: dict[str, Any] = {
        "schema": "leopard2-all-k-gap-cell/v2",
        "run_contract_sha256": contract_digest,
        "cell": dataclasses.asdict(cell),
        "cpu": cpu, "order": list(ORDER), "invocations": list(invocations),
        "valid": not failures, "diagnostic_not_promotion_evidence": True,
    }
    if failures:
        result["failures"] = failures
        return result
    validate_invocation_identities(
        invocations, main_commit, current_source,
        main_snapshot, current_snapshot)
    validate_correctness(cell, invocations)
    main = [record for record in invocations if record["role"] == "main"]
    current = [record for record in invocations if record["role"] == "leopard2"]
    require(len(main) == len(current) == 2,
            "all-K ABBA role cardinality mismatch")
    paths = classify_paths(cell, current[0]["result"])
    require(all(classify_paths(cell, record["result"]) == paths
                for record in current),
            "Leopard2 selected-path identity changed within an all-K cell")
    main_encode = geometric([metric(record, ("metrics", "encode_execution",
        "median_us_per_batch_call")) for record in main])
    current_encode = geometric([metric(record, ("metrics", "encode_execution",
        "median_us_per_batch_call")) for record in current])
    main_decode = geometric([metric(record, ("metrics", "decode_including_setup",
        "median_us_per_batch_call")) for record in main])
    current_plan = geometric([metric(record, ("metrics", "decode_plan_setup",
        "median_us")) for record in current])
    current_decode = geometric([metric(record, ("metrics", "decode_execution",
        "median_us_per_batch_call")) for record in current])
    current_codec = geometric([metric(record, ("metrics", "codec_setup",
        "median_us")) for record in current])
    first = current_plan + current_decode
    amortized = current_decode + current_plan / cell.reuse
    encode_speedup = main_encode / current_encode
    first_speedup = main_decode / first
    reuse_speedup = main_decode / amortized
    paths["gap_tags"] = gap_tags(
        cell, paths, encode_speedup, first_speedup, reuse_speedup)
    result.update({
        "selected": paths,
        "timing_us": {
            "main_encode_execution": main_encode,
            "leopard2_codec_setup": current_codec,
            "leopard2_encode_execution": current_encode,
            "main_decode_including_setup": main_decode,
            "leopard2_decode_plan_setup": current_plan,
            "leopard2_decode_execution": current_decode,
            "leopard2_decode_first_use": first,
            "leopard2_decode_amortized_at_reuse": amortized,
        },
        "speedup_main_over_leopard2": {
            "encode": encode_speedup,
            "decode_first_use": first_speedup,
            "decode_at_reuse": reuse_speedup,
            "decode_execution_only_optimistic": main_decode / current_decode,
        },
        "significantly_beats_main_1_05": {
            "encode": encode_speedup >= 1.05,
            "decode_first_use": first_speedup >= 1.05,
            "decode_at_reuse": reuse_speedup >= 1.05,
        },
        "memory": current[0]["result"]["memory"],
    })
    if all(record["result"]["legacy"]["available"] for record in current):
        current_legacy_encode = geometric([metric(
            record, ("legacy", "encode_execution", "median_us_per_batch_call"))
            for record in current])
        current_legacy_decode = geometric([metric(
            record, ("legacy", "decode_including_setup",
                     "median_us_per_batch_call"))
            for record in current])
        result["timing_us"].update({
            "current_tree_legacy_encode_execution": current_legacy_encode,
            "current_tree_legacy_decode_including_setup": current_legacy_decode,
        })
        result["attribution_speedup"] = {
            "exact_main_over_current_tree_legacy": {
                "encode": main_encode / current_legacy_encode,
                "decode_including_setup": main_decode / current_legacy_decode,
            },
            "current_tree_legacy_over_leopard2": {
                "encode": current_legacy_encode / current_encode,
                "decode_first_use": current_legacy_decode / first,
                "decode_at_reuse": current_legacy_decode / amortized,
                "decode_execution_only_optimistic":
                    current_legacy_decode / current_decode,
            },
        }
    return result


def validate_cell_document(
        document: Mapping[str, Any], cell: Cell, cpu: int,
        contract_digest: str, main: Path, current: Path,
        with_current_legacy: bool, main_commit: str,
        current_source: Mapping[str, Any],
        main_snapshot: Mapping[str, Any],
        current_snapshot: Mapping[str, Any]) -> None:
    require(document.get("schema") == "leopard2-all-k-gap-cell/v2",
            f"stored cell {cell.identifier} schema mismatch")
    require(document.get("run_contract_sha256") == contract_digest,
            f"stored cell {cell.identifier} contract mismatch")
    require(document.get("cell") == dataclasses.asdict(cell),
            f"stored cell {cell.identifier} parameters mismatch")
    require(document.get("cpu") == cpu,
            f"stored cell {cell.identifier} CPU mismatch")
    require(document.get("order") == list(ORDER),
            f"stored cell {cell.identifier} order mismatch")
    invocations = document.get("invocations")
    require(isinstance(invocations, list) and 1 <= len(invocations) <= len(ORDER),
            f"stored cell {cell.identifier} invocation count is invalid")
    for slot, record in enumerate(invocations):
        require(isinstance(record, dict),
                f"stored cell {cell.identifier} invocation is not an object")
        role = ORDER[slot]
        executable = main if role == "main" else current
        require(record.get("role") == role,
                f"stored cell {cell.identifier} invocation role mismatch")
        require(record.get("command") == command(
            role, executable, cell, cpu, with_current_legacy),
            f"stored cell {cell.identifier} command mismatch")
    if document.get("valid") is True:
        require(len(invocations) == len(ORDER),
                f"stored valid cell {cell.identifier} is incomplete")
        validate_invocation_identities(
            invocations, main_commit, current_source,
            main_snapshot, current_snapshot)
        validate_correctness(cell, invocations)
    else:
        require(document.get("valid") is False,
                f"stored cell {cell.identifier} valid flag is not Boolean")


def run_cell(cell: Cell, index: int, cpus: Sequence[int], main: Path,
             current: Path, output: Path, timeout: float,
             with_current_legacy: bool, contract_digest: str,
             main_commit: str,
             current_source: Mapping[str, Any],
             main_snapshot_descriptor: int,
             main_snapshot: Mapping[str, Any],
             current_snapshot_descriptor: int,
             current_snapshot: Mapping[str, Any]) -> dict[str, Any]:
    path = output / "cells" / (cell.identifier + ".json")
    cpu = cpus[index % len(cpus)]
    if path.is_file():
        stored = json.loads(path.read_text(encoding="utf-8"))
        require(isinstance(stored, dict),
                f"stored cell {cell.identifier} is not a JSON object")
        validate_cell_document(
            stored, cell, cpu, contract_digest, main, current,
            with_current_legacy, main_commit, current_source,
            main_snapshot, current_snapshot)
        if stored["valid"] is True:
            return stored
    invocations = []
    for role in ORDER:
        executable = main if role == "main" else current
        snapshot_descriptor = (main_snapshot_descriptor if role == "main" else
                               current_snapshot_descriptor)
        snapshot_identity = (main_snapshot if role == "main" else
                             current_snapshot)
        logical_command = command(
            role, executable, cell, cpu, with_current_legacy)
        try:
            invocations.append(run_one(
                role, logical_command, timeout,
                snapshot_descriptor, snapshot_identity))
        except subprocess.TimeoutExpired as error:
            invocations.append({
                "role": role, "command": logical_command,
                "executable_snapshot_sha256": snapshot_identity["sha256"],
                "returncode": -999,
                "duration_ns": int(timeout * 1e9), "stderr": "timeout",
            })
            break
    result = analyze_cell(
        cell, invocations, cpu, contract_digest, main_commit, current_source,
        main_snapshot, current_snapshot)
    write_json_atomic(path, result)
    return result


def summarize(results: Sequence[Mapping[str, Any]], metadata: Mapping[str, Any]) \
        -> dict[str, Any]:
    valid = [result for result in results if result.get("valid") is True]
    failed = [result for result in results if result.get("valid") is not True]
    metrics = ("encode", "decode_first_use", "decode_at_reuse")
    summary: dict[str, Any] = {
        "schema": "leopard2-all-k-gap-summary/v2", "metadata": dict(metadata),
        "cell_count": len(results), "valid_cell_count": len(valid),
        "failed_cell_count": len(failed),
        "diagnostic_not_promotion_evidence": True,
        "threshold": "Leopard2 significant when main_time/Leopard2_time >= 1.05",
        "metrics": {}, "attribution_metrics": {}, "gap_tags": {},
        "worst_cells": {},
    }
    for name in metrics:
        values = [result["speedup_main_over_leopard2"][name] for result in valid]
        slow = [value for value in values if value < 1.05]
        summary["metrics"][name] = {
            "count": len(values), "gap_count": len(slow),
            "gap_fraction": len(slow) / len(values) if values else None,
            "median_speedup": statistics.median(values) if values else None,
            "p10_speedup": sorted(values)[max(0, len(values) // 10 - 1)]
            if values else None,
        }
        ranked = sorted(valid, key=lambda result:
                        result["speedup_main_over_leopard2"][name])[:30]
        summary["worst_cells"][name] = [{
            "cell": result["cell"], "selected": result["selected"],
            "speedup": result["speedup_main_over_leopard2"][name],
            "timing_us": result["timing_us"],
        } for result in ranked]
    attribution_paths = {
        "exact_main_over_current_tree_legacy_encode":
            ("exact_main_over_current_tree_legacy", "encode"),
        "exact_main_over_current_tree_legacy_decode":
            ("exact_main_over_current_tree_legacy", "decode_including_setup"),
        "current_tree_legacy_over_leopard2_encode":
            ("current_tree_legacy_over_leopard2", "encode"),
        "current_tree_legacy_over_leopard2_decode_first_use":
            ("current_tree_legacy_over_leopard2", "decode_first_use"),
        "current_tree_legacy_over_leopard2_decode_at_reuse":
            ("current_tree_legacy_over_leopard2", "decode_at_reuse"),
    }
    if valid and all("attribution_speedup" in result for result in valid):
        for name, (leg, metric_name) in attribution_paths.items():
            values = [result["attribution_speedup"][leg][metric_name]
                      for result in valid]
            ordered = sorted(values)
            summary["attribution_metrics"][name] = {
                "count": len(values),
                "median_ratio": statistics.median(values),
                "p10_ratio": ordered[max(0, len(ordered) // 10 - 1)],
                "p90_ratio": ordered[min(len(ordered) - 1,
                                         9 * len(ordered) // 10)],
                "count_below_1_0": sum(value < 1.0 for value in values),
                "count_below_1_05": sum(value < 1.05 for value in values),
            }
    tags: dict[str, list[Mapping[str, Any]]] = {}
    for result in valid:
        for tag in result["selected"]["gap_tags"]:
            tags.setdefault(tag, []).append(result)
    summary["gap_tags"] = {
        tag: {
            "cell_count": len(items),
            "median_encode_speedup": statistics.median(
                item["speedup_main_over_leopard2"]["encode"] for item in items),
            "median_decode_first_speedup": statistics.median(
                item["speedup_main_over_leopard2"]["decode_first_use"]
                for item in items),
            "median_decode_reuse_speedup": statistics.median(
                item["speedup_main_over_leopard2"]["decode_at_reuse"]
                for item in items),
        } for tag, items in sorted(tags.items(), key=lambda pair: -len(pair[1]))
    }
    summary["failed_cells"] = [result["cell"] for result in failed]
    return summary


def run(options: argparse.Namespace, lock_path: str) -> int:
    main_commit = require_hex(options.main_commit, "exact-main source commit")
    expected_main_sha256 = require_sha256(
        options.main_sha256, "exact-main executable SHA-256")
    current_commit = require_hex(options.current_commit,
                                 "current source commit")
    current_source_initial = git_identity(
        options.current_source_root, current_commit)
    (main_executable_initial, main_snapshot_descriptor,
     main_snapshot_initial) = snapshot_executable(
        options.main, "exact-main benchmark")
    (current_executable_initial, current_snapshot_descriptor,
     current_snapshot_initial) = snapshot_executable(
        options.current, "Leopard2 benchmark")
    require(main_executable_initial["sha256"] == expected_main_sha256,
            "exact-main executable SHA-256 does not match the frozen "
            "pure-AVX2 baseline")
    main_exe = Path(main_executable_initial["path"])
    current_exe = Path(current_executable_initial["path"])
    main_isa = inspect_pure_avx2(main_exe, "exact-main benchmark")
    current_isa = inspect_pure_avx2(current_exe, "Leopard2 benchmark")
    output = options.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    cpus = sorted(os.sched_getaffinity(0))
    require(options.workers > 0, "workers must be positive")
    require(options.timeout > 0 and math.isfinite(options.timeout),
            "timeout must be finite and positive")
    workers = min(options.workers, len(cpus), 128)
    require(workers > 0, "the process has no allowed CPUs")
    cpus = cpus[:workers]
    current_source_attestation = run_one(
        "leopard2",
        source_attestation_command(current_exe, cpus[0]),
        options.timeout, current_snapshot_descriptor,
        current_snapshot_initial)
    validate_source_attestation(
        current_source_attestation, current_source_initial,
        current_snapshot_initial)
    cells = make_cells(gf8_only=options.gf8_only)
    contract = {
        "schema": "leopard2-all-k-gap-contract/v2",
        "main_commit": main_commit, "current_commit": current_commit,
        "expected_main_sha256": expected_main_sha256,
        "current_source_initial": current_source_initial,
        "main_executable_initial": main_executable_initial,
        "current_executable_initial": current_executable_initial,
        "main_executable_snapshot": main_snapshot_initial,
        "current_executable_snapshot": current_snapshot_initial,
        "current_source_attestation": current_source_attestation,
        "main_isa": main_isa, "current_isa": current_isa,
        "benchmark_lock": lock_path,
        "allowed_cpus": sorted(os.sched_getaffinity(0)), "used_cpus": cpus,
        "workers": workers, "order": list(ORDER),
        "timeout_seconds": options.timeout,
        "with_current_legacy": options.with_current_legacy,
        "matrix": {
            "gf8_K": [1, 255], "gf8_shard_bytes": [4096, 65536],
            "gf16_K": [] if options.gf8_only else list(GF16_K),
            "gf16_shard_bytes": [] if options.gf8_only else [512, 4096],
            "gf8_only": options.gf8_only,
            "cell_count": len(cells),
        },
        "measurement_note": "all CPUs saturated; diagnostic crossover map, not isolated promotion evidence",
    }
    contract_digest = canonical_digest(contract)
    manifest_path = output / "manifest.json"
    manifest: dict[str, Any] = {
        "schema": "leopard2-all-k-gap-manifest/v2",
        "run_contract": contract,
        "run_contract_sha256": contract_digest,
        "cells": [dataclasses.asdict(cell) for cell in cells],
        "completion": None,
    }
    if manifest_path.is_file():
        existing_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        require(isinstance(existing_manifest, dict),
                "existing all-K manifest is not a JSON object")
        validate_manifest(existing_manifest, contract, contract_digest, cells)
        manifest = existing_manifest
    else:
        write_json_atomic(manifest_path, manifest)
    results: list[Mapping[str, Any]] = [None] * len(cells)  # type: ignore[list-item]
    completed = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        future_map = {
            executor.submit(run_cell, cell, index, cpus, main_exe, current_exe,
                            output, options.timeout,
                            options.with_current_legacy, contract_digest,
                            main_commit, current_source_initial,
                            main_snapshot_descriptor, main_snapshot_initial,
                            current_snapshot_descriptor,
                            current_snapshot_initial): index
            for index, cell in enumerate(cells)
        }
        for future in concurrent.futures.as_completed(future_map):
            index = future_map[future]
            results[index] = future.result()
            completed += 1
            if completed % 50 == 0 or completed == len(cells):
                print(f"{completed}/{len(cells)} cells", flush=True)
    current_source_final = git_identity(
        options.current_source_root, current_commit)
    main_executable_final = file_identity(main_exe, "exact-main benchmark")
    current_executable_final = file_identity(
        current_exe, "Leopard2 benchmark")
    main_snapshot_final = sealed_snapshot_identity(
        main_snapshot_descriptor, "exact-main benchmark")
    current_snapshot_final = sealed_snapshot_identity(
        current_snapshot_descriptor, "Leopard2 benchmark")
    require(current_source_final == current_source_initial,
            "current source identity changed during the all-K run")
    require(main_executable_final == main_executable_initial,
            "exact-main executable identity changed during the all-K run")
    require(current_executable_final == current_executable_initial,
            "Leopard2 executable identity changed during the all-K run")
    require(main_snapshot_final == main_snapshot_initial,
            "exact-main sealed executable snapshot changed")
    require(current_snapshot_final == current_snapshot_initial,
            "Leopard2 sealed executable snapshot changed")
    completion = {
        "current_source_final": current_source_final,
        "main_executable_final": main_executable_final,
        "current_executable_final": current_executable_final,
        "main_executable_snapshot_final": main_snapshot_final,
        "current_executable_snapshot_final": current_snapshot_final,
    }
    if manifest.get("completion") is not None:
        require(manifest["completion"] == completion,
                "existing all-K completion identity mismatch")
    manifest["completion"] = completion
    validate_manifest(manifest, contract, contract_digest, cells)
    write_json_atomic(manifest_path, manifest)
    cells_text = "".join(json.dumps(result, sort_keys=True) + "\n"
                         for result in results)
    cells_path = output / "cells.jsonl"
    cells_temporary = cells_path.with_name(cells_path.name + ".tmp")
    cells_temporary.write_text(cells_text, encoding="utf-8")
    os.replace(cells_temporary, cells_path)
    gaps = [result for result in results if result.get("valid") is True and
            not all(result["significantly_beats_main_1_05"].values())]
    gap_text = "".join(json.dumps(result, sort_keys=True) + "\n"
                       for result in gaps)
    gap_path = output / "gap_cells.jsonl"
    gap_temporary = gap_path.with_name(gap_path.name + ".tmp")
    gap_temporary.write_text(gap_text, encoding="utf-8")
    os.replace(gap_temporary, gap_path)
    metadata = {
        "run_contract": contract,
        "run_contract_sha256": contract_digest,
        "completion": completion,
    }
    summary = summarize(results, metadata)
    write_json_atomic(output / "summary.json", summary, pretty=True)
    print(output / "summary.json")
    exit_code = 0 if not summary["failed_cell_count"] else 1
    os.close(current_snapshot_descriptor)
    os.close(main_snapshot_descriptor)
    return exit_code


def main(arguments: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--main", type=Path, required=True)
    parser.add_argument("--main-commit", default=MAIN_COMMIT)
    parser.add_argument("--main-sha256", default=PURE_AVX2_MAIN_SHA256,
                        help="required frozen pure-AVX2 baseline digest")
    parser.add_argument("--current", type=Path, required=True)
    parser.add_argument("--current-source-root", type=Path, required=True)
    parser.add_argument("--current-commit", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int,
                        default=min(128, os.cpu_count() or 1))
    parser.add_argument("--timeout", type=float, default=120.0)
    parser.add_argument("--gf8-only", action="store_true",
                        help="run the 2,522-cell GF8 K=1..255 matrix only")
    parser.add_argument("--lock", type=Path,
                        default=Path("/tmp/leopard-gf8-authoritative.lock"),
                        help="exclusive lock shared by every build/test/timing lane")
    parser.add_argument("--with-current-legacy", action="store_true",
                        help="also time the current tree's retained legacy API")
    options = parser.parse_args(arguments)
    with benchmark_lock(options.lock) as lock_path:
        return run(options, lock_path)


if __name__ == "__main__":
    raise SystemExit(main())
