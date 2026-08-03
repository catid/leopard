#!/usr/bin/env python3
"""Qualify the AVX2 GF8 Algorithm 4 P32/N64/B64 decode terminal.

The runner consumes a clean production build and an exact Leopard-main build,
copies every executable/archive into a lane-owned immutable artifact directory
under the canonical benchmark lock, runs the exhaustive focused correctness
gate, and executes nine mirrored ABBA rounds.  Target cells compare one shared
candidate/control executable and exact main; selector-boundary neighbors prove
that enabling the terminal is inert outside its qualified region.
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import statistics
import subprocess
import sys
import time
from typing import Any, Mapping, Sequence


SCHEMA = "leopard2-alg4-p32-terminal-abba/v1"
SUMMARY_SCHEMA = "leopard2-alg4-p32-terminal-summary/v1"
BENCHMARK_SCHEMA = "leopard2-benchmark-v18"
MAIN_SCHEMA = "leopard-main-benchmark-v1"
MAIN_COMMIT = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
LOCK_PATH = Path("/tmp/leopard-gf8-authoritative.lock")
T95_DF8 = 2.306004135204166
TARGET_CONTROL_FLOOR = 1.05
TARGET_MAIN_FLOOR = 1.05
NEIGHBOR_FLOOR = 1.0 / 1.02
CHILD_ENVIRONMENT = {
    "LANG": "C",
    "LC_ALL": "C",
    "OMP_DYNAMIC": "FALSE",
    "OMP_NUM_THREADS": "1",
    "OMP_PLACES": "cores",
    "OMP_PROC_BIND": "TRUE",
    "PATH": "/usr/bin:/bin",
    "TZ": "UTC",
}


class EvidenceError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvidenceError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_identity(path: Path, *, allow_empty: bool = False) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    status = resolved.stat()
    require(resolved.is_file() and (allow_empty or status.st_size > 0),
            f"required artifact is not a regular file with valid size: "
            f"{resolved}")
    return {
        "path": str(resolved),
        "size": status.st_size,
        "mode": status.st_mode & 0o777,
        "sha256": sha256(resolved),
    }


def write_exclusive(path: Path, value: object) -> None:
    payload = (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o444)
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(payload)
            output.flush()
            os.fsync(output.fileno())
    except BaseException:
        try:
            path.unlink()
        except OSError:
            pass
        raise


def write_bytes_exclusive(path: Path, payload: bytes) -> None:
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o444)
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(payload)
            output.flush()
            os.fsync(output.fileno())
    except BaseException:
        try:
            path.unlink()
        except OSError:
            pass
        raise


def copy_frozen(
    source: Path,
    destination: Path,
    executable: bool,
    expected_sha256: str,
) -> dict[str, Any]:
    source_identity = file_identity(source)
    require(source_identity["sha256"] == expected_sha256,
            f"source changed after build-closure validation: {source}")
    descriptor = os.open(
        destination, os.O_WRONLY | os.O_CREAT | os.O_EXCL,
        0o555 if executable else 0o444)
    try:
        with source.open("rb") as input_file, os.fdopen(descriptor, "wb") as output:
            shutil.copyfileobj(input_file, output, 1024 * 1024)
            output.flush()
            os.fsync(output.fileno())
        destination.chmod(0o555 if executable else 0o444)
    except BaseException:
        try:
            destination.unlink()
        except OSError:
            pass
        raise
    frozen_identity = file_identity(destination)
    require(frozen_identity["sha256"] == expected_sha256 and
            frozen_identity["sha256"] == source_identity["sha256"],
            f"frozen copy differs from source: {source}")
    return {"source": source_identity, "frozen": frozen_identity}


def run_text(arguments: Sequence[str], cwd: Path | None = None) -> str:
    completed = subprocess.run(
        list(arguments), cwd=str(cwd) if cwd else None,
        env=CHILD_ENVIRONMENT, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, timeout=20.0, check=False)
    require(completed.returncode == 0,
            f"command failed ({' '.join(arguments)}): " +
            completed.stderr.decode("utf-8", errors="replace")[-1000:])
    return completed.stdout.decode("utf-8", errors="strict").strip()


def git_identity(source: Path) -> dict[str, Any]:
    head = run_text(("/usr/bin/git", "rev-parse", "HEAD"), source)
    tree = run_text(("/usr/bin/git", "rev-parse", "HEAD^{tree}"), source)
    tracked_status = run_text(
        ("/usr/bin/git", "status", "--porcelain", "--untracked-files=no"),
        source)
    require(not tracked_status, f"tracked source is dirty: {source}")
    return {"path": str(source.resolve()), "head": head, "tree": tree,
            "tracked_status": "clean"}


def archive_members(archive: Path) -> list[str]:
    members = run_text(("/usr/bin/ar", "t", str(archive))).splitlines()
    require(members and len(members) == len(set(members)),
            f"archive is empty or has duplicate members: {archive}")
    return members


def candidate_build_closure(
    build: Path, archive: Path, executable: Path,
    correctness: Path,
) -> dict[str, Any]:
    metadata_paths = {
        "cmake_cache": build / "CMakeCache.txt",
        "compile_commands": build / "compile_commands.json",
        "archive_recipe": build / "CMakeFiles" / "leopard.dir" / "link.txt",
        "benchmark_recipe":
            build / "CMakeFiles" / "bench_leopard2.dir" / "link.txt",
    }
    metadata = {name: file_identity(path)
                for name, path in metadata_paths.items()}
    commands = json.loads(metadata_paths["compile_commands"].read_text())
    require(isinstance(commands, list), "compile_commands is not an array")
    p32 = [entry for entry in commands
           if Path(str(entry.get("file", ""))).name ==
           "Leopard2LowP32B64AVX2.cpp"]
    require(len(p32) == 1, "P32 terminal must compile exactly once")
    command = str(p32[0].get("command", ""))
    require(command.count("-mavx2") == 1 and
            command.count("-mno-avx512f") == 1 and
            "LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL=1" in command,
            "P32 terminal compile flags differ from the AVX2 contract")
    members = archive_members(archive)
    require(members.count("Leopard2LowP32B64AVX2.cpp.o") == 1,
            "production archive does not contain exactly one P32 terminal")
    for required in ("leopard2.cpp.o", "Leopard2BackendAVX2.cpp.o",
                     "Leopard2BackendAVX2Xor.cpp.o"):
        require(required in members,
                f"production archive lacks required member {required}")
    return {
        "build_directory": str(build.resolve()),
        "metadata": metadata,
        "archive": file_identity(archive),
        "benchmark": file_identity(executable),
        "correctness_test": file_identity(correctness),
        "archive_members": members,
        "p32_compile_command_sha256": hashlib.sha256(
            command.encode()).hexdigest(),
    }


def validate_main_manifest(
    path: Path, source: Mapping[str, Any], executable: Path, archive: Path,
) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise EvidenceError(f"exact-main manifest is unreadable: {error}") from error
    require(isinstance(value, dict) and
            value.get("schema") == "leopard2-p32-b64-exact-main-abba/v1",
            "exact-main manifest schema differs")
    baseline = value.get("baseline")
    require(isinstance(baseline, dict) and
            baseline.get("source_commit") == source["head"] and
            baseline.get("source_tree") == source["tree"] and
            baseline.get("pure_avx2") is True and
            baseline.get("executable_sha256") == sha256(executable) and
            baseline.get("archive_sha256") == sha256(archive),
            "exact-main manifest does not bind the supplied source/artifacts")
    return {"file": file_identity(path), "schema": value["schema"],
            "status": value.get("status"), "baseline": dict(baseline)}


def cpu_snapshot(cpu: int) -> dict[str, int]:
    prefix = f"cpu{cpu} "
    for line in Path("/proc/stat").read_text(encoding="ascii").splitlines():
        if line.startswith(prefix):
            fields = [int(value) for value in line.split()[1:]]
            return {"total": sum(fields), "idle": fields[3] + fields[4],
                    "nonidle": sum(fields) - fields[3] - fields[4]}
    raise EvidenceError(f"CPU {cpu} is absent from /proc/stat")


def cpu_delta(before: Mapping[str, int], after: Mapping[str, int]) -> dict[str, int]:
    return {name: after[name] - before[name]
            for name in ("total", "idle", "nonidle")}


def cells() -> list[dict[str, Any]]:
    return [
        {"id": "target-loss9", "K": 32, "R": 32, "bytes": 64,
         "loss": 9, "role": "target", "seed": 0xA4320901},
        {"id": "target-loss16", "K": 32, "R": 32, "bytes": 64,
         "loss": 16, "role": "target", "seed": 0xA4321001},
        {"id": "target-loss31", "K": 32, "R": 32, "bytes": 64,
         "loss": 31, "role": "target", "seed": 0xA4321F01},
        {"id": "byte-neighbor-63", "K": 32, "R": 32, "bytes": 63,
         "loss": 16, "role": "neighbor", "seed": 0xA4323F01},
        {"id": "byte-neighbor-65", "K": 32, "R": 32, "bytes": 65,
         "loss": 16, "role": "neighbor", "seed": 0xA4324101},
        {"id": "loss-neighbor-8", "K": 32, "R": 32, "bytes": 64,
         "loss": 8, "role": "neighbor", "seed": 0xA4320801},
        {"id": "loss-neighbor-32", "K": 32, "R": 32, "bytes": 64,
         "loss": 32, "role": "neighbor", "seed": 0xA4322001},
        {"id": "shape-neighbor-k31", "K": 31, "R": 32, "bytes": 64,
         "loss": 16, "role": "neighbor", "seed": 0xA4311001},
        {"id": "shape-neighbor-r31", "K": 32, "R": 31, "bytes": 64,
         "loss": 16, "role": "neighbor", "seed": 0xA4201001},
    ]


def benchmark_command(
    implementation: str, executable: Path, cell: Mapping[str, Any], cpu: int,
    iterations: int, warmup: int, output: Path,
) -> list[str]:
    common = [
        "/usr/bin/taskset", "-c", str(cpu), "/usr/bin/prlimit",
        "--as=268435456", "--", str(executable),
        "--k", str(cell["K"]), "--r", str(cell["R"]),
        "--bytes", str(cell["bytes"]), "--loss", str(cell["loss"]),
        "--batch", "1", "--reuse", "256", "--iterations", str(iterations),
        "--warmup", str(warmup), "--threads", "1",
        "--seed", str(cell["seed"]),
    ]
    if implementation == "main":
        return common + ["--json", str(output)]
    mode = "1" if implementation == "candidate" else "0"
    return common + [
        "--profile", "high", "--field", "gf8", "--backend", "avx2",
        "--skip-legacy", "--retain-samples", "--measure-one-shot-decode",
        "--low-p32-b64-terminal-mode", mode, "--attest-source",
        "--json", str(output),
    ]


def validate_result(
    implementation: str, result: object, cell: Mapping[str, Any],
    source: Mapping[str, Any], iterations: int, warmup: int,
) -> dict[str, Any]:
    require(isinstance(result, dict), "benchmark output is not a JSON object")
    parameters = result.get("parameters")
    require(isinstance(parameters, dict), "benchmark parameters are absent")
    expected = {
        "K": cell["K"], "R": cell["R"],
        "shard_bytes": cell["bytes"], "loss_count": cell["loss"],
        "batch": 1, "reuse": 256, "iterations": iterations,
        "warmup": warmup, "thread_count": 1, "seed": cell["seed"],
    }
    require(all(parameters.get(key) == value for key, value in expected.items()),
            "benchmark parameters differ from the frozen cell")
    resolved = result.get("resolved")
    digests = result.get("workload_digests")
    require(isinstance(resolved, dict) and
            resolved.get("profile") == "legacy_high_v1" and
            resolved.get("field") == "gf8" and
            resolved.get("thread_count") == 1,
            "benchmark resolved another profile, field, or thread count")
    require(isinstance(digests, dict) and digests.get("algorithm") == "fnv1a64",
            "workload digests are absent")
    for name in ("original_data", "transmitted_parity", "recovered_originals"):
        require(isinstance(digests.get(name), str) and len(digests[name]) == 16,
                f"workload digest is invalid: {name}")
    metrics = result.get("metrics")
    require(isinstance(metrics, dict), "metrics are absent")
    if implementation == "main":
        require(result.get("schema") == MAIN_SCHEMA and
                result.get("correctness", {}).get("round_trip") is True and
                result.get("build", {}).get("main_source_commit") == MAIN_COMMIT,
                "exact-main identity or correctness differs")
        metric = metrics.get("decode_including_setup")
    else:
        build = result.get("build")
        expected_selected = implementation == "candidate" and \
            cell["role"] == "target"
        require(result.get("schema") == BENCHMARK_SCHEMA and
                resolved.get("backend") == "avx2" and
                result.get("correctness", {}).get("leopard2_round_trip") is True,
                "Leopard2 schema, backend, or correctness differs")
        require(isinstance(build, dict) and
                build.get("source_commit") == source["head"] and
                build.get("source_tree") == source["tree"] and
                build.get("source_tracked_dirty") is False,
                "embedded candidate source identity differs")
        require(build.get("low_p32_b64_terminal_mode_word") ==
                (1 if implementation == "candidate" else 2),
                "runtime selector mode differs from the label")
        require(build.get(
            "low_p32_b64_terminal_one_shot_route_expected_selected") is
                expected_selected and
                build.get("low_p32_b64_terminal_one_shot_route_selected") is
                expected_selected,
                "actual terminal route differs from the frozen cell")
        metric = metrics.get("one_shot_decode_including_setup")
    require(isinstance(metric, dict), "decode timing summary is absent")
    samples = metric.get("samples_us_per_batch_call")
    require(isinstance(samples, list) and len(samples) == iterations and
            all(isinstance(value, (int, float)) and
                not isinstance(value, bool) and
                math.isfinite(float(value)) and float(value) > 0
                for value in samples),
            "decode raw samples are incomplete")
    numeric = [float(value) for value in samples]
    median = statistics.median(numeric)
    mad = statistics.median(abs(value - median) for value in numeric)
    reported = {
        "median": metric.get("median_us_per_batch_call"),
        "mad": metric.get("mad_us_per_batch_call"),
        "minimum": metric.get("minimum_us_per_batch_call"),
        "maximum": metric.get("maximum_us_per_batch_call"),
    }
    require(all(isinstance(value, (int, float)) and
                not isinstance(value, bool) and math.isfinite(float(value))
                for value in reported.values()) and
            math.isclose(float(reported["median"]), median,
                         rel_tol=0.0, abs_tol=2e-6) and
            math.isclose(float(reported["mad"]), mad,
                         rel_tol=0.0, abs_tol=2e-6) and
            math.isclose(float(reported["minimum"]), min(numeric),
                         rel_tol=0.0, abs_tol=2e-6) and
            math.isclose(float(reported["maximum"]), max(numeric),
                         rel_tol=0.0, abs_tol=2e-6),
            "decode summary disagrees with retained raw samples")
    return {"decode_us": median, "decode_samples_us": numeric,
            "digests": dict(digests)}


def run_one(
    implementation: str, executable: Path, identity: Mapping[str, Any],
    cell: Mapping[str, Any], cpu: int, source: Mapping[str, Any],
    iterations: int, warmup: int, directory: Path, label: str,
) -> dict[str, Any]:
    require(sha256(executable) == identity["sha256"],
            f"{implementation} executable changed before execution")
    result_path = directory / f"{label}.json"
    stdout_path = directory / f"{label}.stdout"
    stderr_path = directory / f"{label}.stderr"
    command = benchmark_command(
        implementation, executable, cell, cpu, iterations, warmup, result_path)
    started = time.monotonic_ns()
    try:
        completed = subprocess.run(
            command, env=CHILD_ENVIRONMENT, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, timeout=60.0, check=False)
    except subprocess.TimeoutExpired as error:
        write_bytes_exclusive(stdout_path, error.stdout or b"")
        write_bytes_exclusive(stderr_path, error.stderr or b"")
        raise EvidenceError(f"{implementation} timed out") from error
    write_bytes_exclusive(stdout_path, completed.stdout)
    write_bytes_exclusive(stderr_path, completed.stderr)
    require(completed.returncode == 0,
            f"{implementation} exited {completed.returncode}; see {stderr_path}")
    require(sha256(executable) == identity["sha256"],
            f"{implementation} executable changed after execution")
    try:
        result = json.loads(result_path.read_text(encoding="utf-8"))
        normalized = validate_result(
            implementation, result, cell, source, iterations, warmup)
    except (OSError, UnicodeError, json.JSONDecodeError, EvidenceError) as error:
        raise EvidenceError(
            f"{implementation} JSON/validation failure; retained {result_path}, "
            f"{stdout_path}, and {stderr_path}: {error}") from error
    return {
        "implementation": implementation,
        "command": command,
        "elapsed_ns": time.monotonic_ns() - started,
        "result": file_identity(result_path),
        "stdout": file_identity(stdout_path, allow_empty=True),
        "stderr": file_identity(stderr_path, allow_empty=True),
        "normalized": normalized,
    }


def confidence(round_logs: Sequence[float]) -> dict[str, Any]:
    require(len(round_logs) == 9, "nine independent rounds are required")
    center = statistics.mean(round_logs)
    half = T95_DF8 * statistics.stdev(round_logs) / math.sqrt(9)
    return {"speedup": math.exp(center),
            "ci95": [math.exp(center - half), math.exp(center + half)],
            "round_log_ratios": list(round_logs)}


def analyze(cell: Mapping[str, Any], rounds: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    reference = rounds[0]["invocations"][0]["normalized"]["digests"]
    labels = ("control", "main") if cell["role"] == "target" else ("control",)
    contrasts = {name: [] for name in labels}
    for round_value in rounds:
        invocations = round_value["invocations"]
        require(all(item["normalized"]["digests"] == reference
                    for item in invocations), "workload digests differ")
        candidates = [item["normalized"]["decode_us"] for item in invocations
                      if item["implementation"] == "candidate"]
        require(len(candidates) == 2, "round lacks two candidate observations")
        candidate_log = statistics.mean(math.log(value) for value in candidates)
        for label in labels:
            values = [item["normalized"]["decode_us"] for item in invocations
                      if item["implementation"] == label]
            require(len(values) == 2, f"round lacks two {label} observations")
            contrasts[label].append(
                statistics.mean(math.log(value) for value in values) -
                candidate_log)
    result = {"cell": dict(cell), "digests": reference,
              "control_over_candidate": confidence(contrasts["control"])}
    if "main" in contrasts:
        result["main_over_candidate"] = confidence(contrasts["main"])
    return result


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-build", required=True, type=Path)
    parser.add_argument("--candidate-source", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--candidate-archive", required=True, type=Path)
    parser.add_argument("--correctness-test", required=True, type=Path)
    parser.add_argument("--main-source", required=True, type=Path)
    parser.add_argument("--main", required=True, type=Path)
    parser.add_argument("--main-archive", required=True, type=Path)
    parser.add_argument("--main-manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cpu", required=True, type=int)
    parser.add_argument("--sibling", required=True, type=int)
    parser.add_argument("--iterations", type=int, default=21)
    parser.add_argument("--warmup", type=int, default=5)
    return parser.parse_args()


def main() -> int:
    options = parse_arguments()
    require(options.iterations >= 3 and options.warmup >= 1,
            "insufficient timing samples")
    require(not options.output.exists(), "output directory already exists")
    require(set(os.sched_getaffinity(0)) == {options.cpu},
            "runner must be singleton-pinned to the benchmark CPU")
    siblings = Path(
        f"/sys/devices/system/cpu/cpu{options.cpu}/topology/"
        "thread_siblings_list").read_text(encoding="ascii").strip()
    sibling_values = set()
    for part in siblings.split(","):
        if "-" in part:
            first, last = (int(value) for value in part.split("-", 1))
            sibling_values.update(range(first, last + 1))
        else:
            sibling_values.add(int(part))
    require(sibling_values == {options.cpu, options.sibling},
            "requested CPUs are not one SMT pair")

    lock = os.open(LOCK_PATH, os.O_RDWR | os.O_CREAT, 0o600)
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError as error:
        os.close(lock)
        raise EvidenceError(f"canonical benchmark lock is busy: {LOCK_PATH}") \
            from error

    raw: dict[str, Any] = {"schema": SCHEMA, "started_ns": time.time_ns()}
    try:
        options.output.mkdir(parents=True)
        artifacts = options.output / "artifacts"
        invocations = options.output / "invocations"
        artifacts.mkdir()
        invocations.mkdir()
        candidate_source = git_identity(options.candidate_source)
        main_source = git_identity(options.main_source)
        require(main_source["head"] == MAIN_COMMIT,
                "exact-main checkout is not the frozen main commit")
        main_manifest = validate_main_manifest(
            options.main_manifest, main_source,
            options.main, options.main_archive)
        build_closure = candidate_build_closure(
            options.candidate_build, options.candidate_archive,
            options.candidate, options.correctness_test)
        frozen = {
            "candidate": copy_frozen(
                options.candidate, artifacts / "bench_leopard2", True,
                build_closure["benchmark"]["sha256"]),
            "candidate_archive": copy_frozen(
                options.candidate_archive, artifacts / "libleopard.a", False,
                build_closure["archive"]["sha256"]),
            "correctness": copy_frozen(
                options.correctness_test,
                artifacts / "leopard2_low_p32_b64_terminal_test", True,
                build_closure["correctness_test"]["sha256"]),
            "main": copy_frozen(
                options.main, artifacts / "leopard_main_benchmark", True,
                main_manifest["baseline"]["executable_sha256"]),
            "main_archive": copy_frozen(
                options.main_archive, artifacts / "libleopard_main_exact.a", False,
                main_manifest["baseline"]["archive_sha256"]),
            "main_manifest": copy_frozen(
                options.main_manifest, artifacts / "main_manifest.json", False,
                main_manifest["file"]["sha256"]),
        }
        candidate = Path(frozen["candidate"]["frozen"]["path"])
        main_executable = Path(frozen["main"]["frozen"]["path"])
        correctness = Path(frozen["correctness"]["frozen"]["path"])
        identities = {
            "candidate": frozen["candidate"]["frozen"],
            "control": frozen["candidate"]["frozen"],
            "main": frozen["main"]["frozen"],
        }
        raw.update({
            "runner": file_identity(Path(__file__)),
            "candidate_source": candidate_source,
            "main_source": main_source,
            "validated_main_manifest": main_manifest,
            "candidate_build_closure": build_closure,
            "frozen_artifacts": frozen,
            "cpu": options.cpu,
            "reserved_sibling": options.sibling,
            "iterations": options.iterations,
            "warmup": options.warmup,
            "cells": [],
        })

        correctness_stdout = invocations / "correctness.stdout"
        correctness_stderr = invocations / "correctness.stderr"
        try:
            completed = subprocess.run(
                ["/usr/bin/taskset", "-c", str(options.cpu),
                 "/usr/bin/prlimit", "--as=268435456", "--", str(correctness)],
                env=CHILD_ENVIRONMENT, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, timeout=60.0, check=False)
        except subprocess.TimeoutExpired as error:
            write_bytes_exclusive(correctness_stdout, error.stdout or b"")
            write_bytes_exclusive(correctness_stderr, error.stderr or b"")
            raise EvidenceError(
                "focused correctness gate timed out; output retained") from error
        write_bytes_exclusive(correctness_stdout, completed.stdout)
        write_bytes_exclusive(correctness_stderr, completed.stderr)
        expected_pass = (b"PASS low_p32_b64_terminal payloads=2 patterns=67 "
                         b"parity_selections=2 routes=56\n")
        require(completed.returncode == 0 and completed.stdout == expected_pass,
                "focused correctness gate failed or changed its coverage")
        raw["correctness"] = {
            "stdout": file_identity(correctness_stdout, allow_empty=True),
            "stderr": file_identity(correctness_stderr, allow_empty=True),
            "coverage": {"payloads": 2, "patterns": 67,
                         "parity_selections": 2, "public_routes": 56,
                         "qualified_loss_counts": list(range(9, 32))},
        }

        for cell_index, cell in enumerate(cells()):
            cell_rounds = []
            for round_index in range(9):
                if cell["role"] == "target":
                    target_orders = (
                        ("main", "candidate", "control", "control", "candidate", "main"),
                        ("control", "main", "candidate", "candidate", "main", "control"),
                        ("candidate", "control", "main", "main", "control", "candidate"),
                    )
                    order = target_orders[round_index % len(target_orders)]
                elif round_index % 2 == 0:
                    order = ("control", "candidate", "candidate", "control")
                else:
                    order = ("candidate", "control", "control", "candidate")
                before_cpu = cpu_snapshot(options.cpu)
                before_sibling = cpu_snapshot(options.sibling)
                values = []
                for slot, implementation in enumerate(order):
                    label = (f"{cell['id']}-round{round_index}-slot{slot}-"
                             f"{implementation}")
                    executable = main_executable if implementation == "main" \
                        else candidate
                    values.append(run_one(
                        implementation, executable, identities[implementation],
                        cell, options.cpu, candidate_source,
                        options.iterations, options.warmup, invocations, label))
                isolation = {
                    "benchmark_cpu": cpu_delta(
                        before_cpu, cpu_snapshot(options.cpu)),
                    "reserved_sibling": cpu_delta(
                        before_sibling, cpu_snapshot(options.sibling)),
                }
                require(isolation["reserved_sibling"]["nonidle"] == 0,
                        f"SMT sibling was active in {cell['id']} round {round_index}")
                cell_rounds.append({"round": round_index, "order": list(order),
                                    "invocations": values,
                                    "isolation": isolation})
            raw["cells"].append({"cell": dict(cell), "rounds": cell_rounds})
            print(f"{cell_index + 1}/{len(cells())} {cell['id']}",
                  file=sys.stderr, flush=True)

        require(file_identity(Path(__file__)) == raw["runner"],
                "qualification runner changed during the campaign")
        require(candidate_build_closure(
                    options.candidate_build, options.candidate_archive,
                    options.candidate, options.correctness_test) ==
                build_closure,
                "candidate build closure changed during the campaign")
        require(validate_main_manifest(
                    options.main_manifest, main_source,
                    options.main, options.main_archive) == main_manifest,
                "exact-main build closure changed during the campaign")
        require(git_identity(options.candidate_source) == candidate_source and
                git_identity(options.main_source) == main_source,
                "source identity changed during the campaign")
        for name, item in frozen.items():
            identity = item["frozen"]
            require(file_identity(Path(identity["path"])) == identity,
                    f"frozen artifact changed during the campaign: {name}")
        analyses = [analyze(item["cell"], item["rounds"])
                    for item in raw["cells"]]
        target_control_failure = any(
            item["control_over_candidate"]["ci95"][0] < TARGET_CONTROL_FLOOR
            for item in analyses if item["cell"]["role"] == "target")
        target_main_failure = any(
            item["main_over_candidate"]["ci95"][0] < TARGET_MAIN_FLOOR
            for item in analyses if item["cell"]["role"] == "target")
        neighbor_regressions = [
            item["cell"]["id"] for item in analyses
            if item["cell"]["role"] == "neighbor" and
            item["control_over_candidate"]["ci95"][1] < NEIGHBOR_FLOOR]
        raw["completed_ns"] = time.time_ns()
        write_exclusive(options.output / "raw.json", raw)
        summary = {
            "schema": SUMMARY_SCHEMA,
            "status": "accepted" if not (
                target_control_failure or target_main_failure or
                neighbor_regressions) else "rejected",
            "candidate_source": candidate_source,
            "main_source": main_source,
            "binary_sha256": {name: value["sha256"]
                              for name, value in identities.items()},
            "correctness": raw["correctness"]["coverage"],
            "all_rounds_zero_sibling_nonidle": True,
            "target_control_failure": target_control_failure,
            "target_main_failure": target_main_failure,
            "credible_neighbor_regressions": neighbor_regressions,
            "cells": analyses,
            "raw_sha256": sha256(options.output / "raw.json"),
        }
        write_exclusive(options.output / "summary.json", summary)
        print(json.dumps({"status": summary["status"],
                          "targets": [item for item in analyses
                                      if item["cell"]["role"] == "target"],
                          "neighbor_regressions": neighbor_regressions},
                         sort_keys=True))
        return 0 if summary["status"] == "accepted" else 2
    except BaseException as error:
        raw["failure"] = {"type": type(error).__name__, "message": str(error)}
        if options.output.exists() and not (options.output / "failure.json").exists():
            write_exclusive(options.output / "failure.json", raw)
        print(f"evidence rejected: {error}", file=sys.stderr)
        return 1
    finally:
        os.close(lock)


if __name__ == "__main__":
    raise SystemExit(main())
