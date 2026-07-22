#!/usr/bin/env python3
"""Run and merge the exhaustive small-direct coefficient-matrix verifier."""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


SCHEMA = "leopard2-small-direct-exhaustive-campaign/v1"
SHARD_SCHEMA = "leopard2-small-direct-exhaustive/v1"
EXPECTED_PATTERNS = 1_982_812
DEFAULT_LOCK = Path("/tmp/leopard-gf8-authoritative.lock")


class EvidenceError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvidenceError(message)


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) +
            "\n").encode()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def file_identity(path: Path) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    payload = resolved.read_bytes()
    return {
        "path": str(resolved),
        "size": len(payload),
        "sha256": sha256_bytes(payload),
    }


def atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp-%d" % os.getpid())
    with temporary.open("wb") as stream:
        stream.write(json.dumps(value, sort_keys=True, indent=2).encode())
        stream.write(b"\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def acquire_lock(path: Path, timeout: float):
    require(path == DEFAULT_LOCK and path.is_absolute(),
            "exhaustive campaign requires the canonical GF8 lock")
    path.parent.mkdir(parents=True, exist_ok=True)
    stream = path.open("a+")
    deadline = time.monotonic() + timeout
    while True:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            return stream
        except BlockingIOError:
            if time.monotonic() >= deadline:
                stream.close()
                raise EvidenceError("timed out waiting for canonical GF8 lock")
            time.sleep(1.0)


def validate_shard(
        value: Any, shard_index: int, shard_count: int) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
        "schema", "mode", "shard_index", "shard_count", "total_patterns",
        "assigned_patterns", "recovered_shards", "verified_basis_symbols",
        "basis_seed", "assignment", "digest_fnv1a64",
    }, "shard result has the wrong shape")
    expected_assigned = EXPECTED_PATTERNS // shard_count + (
        1 if shard_index < EXPECTED_PATTERNS % shard_count else 0)
    require(value["schema"] == SHARD_SCHEMA and
            value["mode"] in (1, 2) and
            value["shard_index"] == shard_index and
            value["shard_count"] == shard_count and
            value["total_patterns"] == EXPECTED_PATTERNS and
            value["assigned_patterns"] == expected_assigned and
            type(value["recovered_shards"]) is int and
            value["recovered_shards"] >= expected_assigned * 5 and
            type(value["verified_basis_symbols"]) is int and
            value["verified_basis_symbols"] >=
                value["recovered_shards"] * 5 and
            value["basis_seed"] == 0 and
            value["assignment"] == "global_ordinal_mod_shard_count" and
            isinstance(value["digest_fnv1a64"], str) and
            re.fullmatch(r"[0-9a-f]{16}", value["digest_fnv1a64"]),
            "shard result identity/counts are invalid")
    return value


def retained_shard(
        root: Path, shard_index: int, shard_count: int) -> dict[str, Any] | None:
    directory = root / ("shard-%02d" % shard_index)
    result_path = directory / "result.json"
    stdout_path = directory / "stdout.json"
    stderr_path = directory / "stderr.txt"
    if not result_path.exists():
        return None
    require(stdout_path.is_file() and stderr_path.is_file(),
            "retained shard is missing stdout/stderr")
    record = json.loads(result_path.read_text())
    require(isinstance(record, dict) and set(record) == {
        "shard", "stdout", "stderr",
    }, "retained shard envelope has the wrong shape")
    stdout_identity = file_identity(stdout_path)
    stderr_identity = file_identity(stderr_path)
    require(record["stdout"] == stdout_identity and
            record["stderr"] == stderr_identity,
            "retained shard stream identity changed")
    parsed = json.loads(stdout_path.read_text())
    require(record["shard"] == parsed,
            "retained shard result differs from stdout")
    return validate_shard(parsed, shard_index, shard_count)


def launch_shards(
        binary: Path, output: Path, workers: int, timeout: float,
        shard_count: int) -> list[dict[str, Any]]:
    allowed = sorted(os.sched_getaffinity(0))
    require(workers > 0 and workers <= len(allowed),
            "workers exceed the process-visible CPU set")
    require(shard_count == workers,
            "one deterministic shard per worker is required")
    taskset = shutil.which("taskset")
    require(taskset is not None, "taskset is required for pinned shard jobs")

    results: list[dict[str, Any] | None] = [None] * shard_count
    pending: list[int] = []
    for shard_index in range(shard_count):
        retained = retained_shard(output, shard_index, shard_count)
        if retained is None:
            pending.append(shard_index)
        else:
            results[shard_index] = retained
    if not pending:
        return [value for value in results if value is not None]

    processes: dict[int, tuple[subprocess.Popen[bytes], Any, Any, float]] = {}
    try:
        for shard_index in pending:
            directory = output / ("shard-%02d" % shard_index)
            directory.mkdir(parents=True, exist_ok=True)
            stdout_stream = (directory / "stdout.json").open("wb")
            stderr_stream = (directory / "stderr.txt").open("wb")
            cpu = allowed[shard_index % workers]
            command = [
                taskset, "-c", str(cpu), str(binary),
                "--shard-index", str(shard_index),
                "--shard-count", str(shard_count),
            ]
            process = subprocess.Popen(
                command, stdout=stdout_stream, stderr=stderr_stream,
                start_new_session=True)
            processes[shard_index] = (
                process, stdout_stream, stderr_stream, time.monotonic())

        completed_count = shard_count - len(pending)
        while processes:
            for shard_index in list(processes):
                process, stdout_stream, stderr_stream, started = \
                    processes[shard_index]
                return_code = process.poll()
                if return_code is None and time.monotonic() - started > timeout:
                    process.kill()
                    return_code = process.wait()
                if return_code is None:
                    continue
                stdout_stream.close()
                stderr_stream.close()
                del processes[shard_index]
                directory = output / ("shard-%02d" % shard_index)
                stdout_path = directory / "stdout.json"
                stderr_path = directory / "stderr.txt"
                require(return_code == 0,
                        "shard %d failed with %d; see %s" % (
                            shard_index, return_code, stderr_path))
                try:
                    parsed = json.loads(stdout_path.read_text())
                except json.JSONDecodeError as error:
                    raise EvidenceError(
                        "shard %d emitted invalid JSON" % shard_index) from error
                value = validate_shard(parsed, shard_index, shard_count)
                envelope = {
                    "shard": value,
                    "stdout": file_identity(stdout_path),
                    "stderr": file_identity(stderr_path),
                }
                atomic_json(directory / "result.json", envelope)
                results[shard_index] = value
                completed_count += 1
                print("[%d/%d] shard %d verified %d matrices" % (
                    completed_count, shard_count, shard_index,
                    value["assigned_patterns"]), flush=True)
            if processes:
                time.sleep(0.05)
    finally:
        for process, stdout_stream, stderr_stream, unused_started in \
                processes.values():
            try:
                process.kill()
                process.wait(timeout=5)
            except (ProcessLookupError, subprocess.TimeoutExpired):
                pass
            stdout_stream.close()
            stderr_stream.close()
    require(all(value is not None for value in results),
            "exhaustive shard set is incomplete")
    return [value for value in results if value is not None]


def run(options: argparse.Namespace) -> int:
    binary = options.binary.resolve(strict=True)
    require(binary.is_file() and os.access(binary, os.X_OK),
            "exhaustive verifier is not executable")
    output = options.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    workers = options.workers
    request = {
        "schema": SCHEMA,
        "binary": file_identity(binary),
        "runner": file_identity(Path(__file__)),
        "workers": workers,
        "shard_count": workers,
        "expected_patterns": EXPECTED_PATTERNS,
        "timeout_seconds_per_shard": options.timeout,
        "allowed_cpus": sorted(os.sched_getaffinity(0)),
        "assignment": "global_ordinal_mod_shard_count",
        "basis_seed": 0,
    }
    request["digest"] = sha256_bytes(canonical_bytes(request))
    request_path = output / "request.json"
    if request_path.exists():
        require(json.loads(request_path.read_text()) == request,
                "resume request differs from retained campaign")
    else:
        atomic_json(request_path, request)

    lock = acquire_lock(options.lock, options.lock_timeout)
    try:
        shards = launch_shards(
            binary, output, workers, options.timeout, workers)
    finally:
        fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
        lock.close()
    modes = {value["mode"] for value in shards}
    require(len(modes) == 1,
            "exhaustive shards came from different experiment modes")
    require(sum(value["assigned_patterns"] for value in shards) ==
            EXPECTED_PATTERNS,
            "merged shard count does not cover the exact matrix")
    manifest = {
        "schema": SCHEMA,
        "request": request,
        "mode": next(iter(modes)),
        "shard_count": len(shards),
        "verified_patterns": sum(
            value["assigned_patterns"] for value in shards),
        "recovered_shards": sum(
            value["recovered_shards"] for value in shards),
        "verified_basis_symbols": sum(
            value["verified_basis_symbols"] for value in shards),
        "shards": shards,
    }
    manifest["digest"] = sha256_bytes(canonical_bytes(manifest))
    atomic_json(output / "manifest.json", manifest)
    print("verified %d coefficient matrices in %d pinned shards (mode %d)" % (
        manifest["verified_patterns"], len(shards), manifest["mode"]))
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--binary", required=True, type=Path)
    result.add_argument("--output", required=True, type=Path)
    result.add_argument("--workers", type=int, default=28)
    result.add_argument("--timeout", type=float, default=3600.0)
    result.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    result.add_argument("--lock-timeout", type=float, default=3600.0)
    return result


def main() -> int:
    try:
        return run(parser().parse_args())
    except (EvidenceError, OSError, subprocess.SubprocessError, ValueError) as error:
        print("error: %s" % error, file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
