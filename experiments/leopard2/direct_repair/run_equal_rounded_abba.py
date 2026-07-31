#!/usr/bin/env python3
"""Authoritative equal-rounded GF8/AVX2 direct-repair comparison.

The runner consumes three already-built executables: a current-source
candidate with equal-rounded multi-loss repair enabled, a current-source
control with the same code and initialized-data selector disabled, and exact
Leopard main. It freezes each binary into a lane-owned directory, verifies
candidate/control executable sections are identical, executes every child on
one logical CPU, and rejects any invocation whose target or reserved-sibling
scheduler accounting exceeds the predeclared contamination bounds.

The coordinator itself must be pinned away from the measured CPU pair. No
foreign affinity is changed, so a killed coordinator cannot strand or
mis-restore another task's CPU mask.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Mapping, Sequence


SEED_SCHEMA = "leopard2-equal-rounded-abba/v1"
MANIFEST_SCHEMA = "leopard2-equal-rounded-manifest/v7"
CELL_SCHEMA = "leopard2-equal-rounded-cell/v7"
SUMMARY_SCHEMA = "leopard2-equal-rounded-summary/v7"
ONE_SHOT_MANIFEST_SCHEMA = "leopard2-equal-rounded-one-shot-manifest/v1"
ONE_SHOT_CELL_SCHEMA = "leopard2-equal-rounded-one-shot-cell/v1"
ONE_SHOT_SUMMARY_SCHEMA = "leopard2-equal-rounded-one-shot-summary/v1"
MAIN_COMMIT = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
LOCK_PATH = Path("/tmp/leopard-gf8-authoritative.lock")
TARGET_ORDER = (
    ("main", "candidate", "control", "control", "candidate", "main"),
    ("control", "main", "candidate", "candidate", "main", "control"),
    ("candidate", "control", "main", "main", "control", "candidate"),
)
NEIGHBOR_ORDER = (
    ("control", "candidate", "candidate", "control"),
    ("candidate", "control", "control", "candidate"),
    ("control", "candidate", "candidate", "control"),
)
T95 = {
    3: 4.302652729911275,
    9: 2.306004135204166,
}
TARGET_CONTROL_FLOOR = 1.05
TARGET_MAIN_FLOOR = 1.0
NEIGHBOR_FLOOR = 1.0 / 1.02
MAX_ATTEMPTS = 5
MAX_JSON_BYTES = 32 << 20
# CPU-wide and per-task schedstat clocks have both fixed and runtime-scaled
# endpoint disagreement on this host.  Across 645 otherwise-clean retained
# invocations, max(20 us, 50 ppm) accepts 630 while rejecting all 15 observed
# 116--456 ppm outliers.  The floor is at most 0.264% of the shortest observed
# child runtime and remains far below the 2% neighboring-cell regression gate.
TARGET_ENDPOINT_ABSOLUTE_TOLERANCE_NS = 20_000
TARGET_ENDPOINT_RELATIVE_TOLERANCE_PPM = 50
# The otherwise-clean sibling-runtime distribution is bimodal: ordinary
# periodic work reaches 175 ppm, while contention outliers start above
# 280 ppm.  A 50 us floor is at most 0.66% of the shortest retained runtime;
# the 200 ppm bound limits overlap to 0.02% on long invocations.
SIBLING_ABSOLUTE_TOLERANCE_NS = 50_000
SIBLING_RELATIVE_TOLERANCE_PPM = 200
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


def load_module(name: str, path: Path) -> Any:
    resolved = path.resolve(strict=True)
    specification = importlib.util.spec_from_file_location(name, resolved)
    require(specification is not None and specification.loader is not None,
            f"cannot import evidence support: {resolved}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


DIRECTORY = Path(__file__).resolve().parent
GATE = load_module(
    "leopard2_equal_rounded_gate_support",
    DIRECTORY / "run_small_direct_abba.py")
T8 = load_module(
    "leopard2_equal_rounded_file_support",
    DIRECTORY.parent / "gf8_high_encode" / "run_t8_two_block_abba.py")
PAIR = GATE.PROCESS_SUPPORT


def canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"),
        allow_nan=False).encode("utf-8")


def digest_object(value: object) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def stable_seed(identifier: str) -> int:
    material = (SEED_SCHEMA + ":" + identifier).encode("ascii")
    return int.from_bytes(hashlib.sha256(material).digest()[:8], "big")


def ceil64(value: int) -> int:
    require(type(value) is int and value > 0, "byte count is not positive")
    return (value + 63) & ~63


def cell(
    identifier: str,
    k: int,
    r: int,
    shard_bytes: int,
    losses: int,
    role: str,
    candidate_direct: bool,
    control_direct: bool,
) -> dict[str, Any]:
    require(
        role in ("target", "neighbor") and
        0 < losses <= min(k, r) and shard_bytes > 0,
        f"invalid cell: {identifier}")
    return {
        "id": identifier,
        "K": k,
        "R": r,
        "bytes": shard_bytes,
        "main_physical_bytes": ceil64(shard_bytes),
        "loss": losses,
        "role": role,
        "candidate_direct": candidate_direct,
        "control_direct": control_direct,
        "seed": stable_seed(identifier),
    }


def reusable_matrix() -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for k in (17, 31, 32, 33, 64, 66, 96, 127, 128):
        for losses in (2, 4, 8):
            result.append(cell(
                f"target-k{k}-r{k}-b65536-l{losses}",
                k, k, 65536, losses, "target", True, False))
    for k in (17, 128):
        for shard_bytes in (1, 31, 32, 33, 63, 64, 65):
            result.append(cell(
                f"tail-k{k}-r{k}-b{shard_bytes}-l8",
                k, k, shard_bytes, 8, "target", True, False))
    for k in (17, 66, 128):
        for losses in (2, 8):
            result.append(cell(
                f"small-k{k}-r{k}-b1024-l{losses}",
                k, k, 1024, losses, "target", True, False))
    result.extend((
        cell("neighbor-k65-r65-b65536-l8",
             65, 65, 65536, 8, "neighbor", True, True),
        cell("neighbor-k17-r16-b65536-l2",
             17, 16, 65536, 2, "neighbor", False, False),
        cell("neighbor-k17-r16-b65536-l8",
             17, 16, 65536, 8, "neighbor", False, False),
        cell("neighbor-k32-r16-b65536-l2",
             32, 16, 65536, 2, "neighbor", False, False),
        cell("neighbor-k32-r16-b65536-l8",
             32, 16, 65536, 8, "neighbor", False, False),
    ))
    require(len({item["id"] for item in result}) == len(result),
            "matrix contains duplicate identifiers")
    return result


def one_shot_cell(
    identifier: str,
    k: int,
    r: int,
    shard_bytes: int,
    losses: int,
    role: str,
    candidate_direct: bool,
    control_direct: bool,
) -> dict[str, Any]:
    result = cell(
        identifier, k, r, shard_bytes, losses, role,
        candidate_direct, control_direct)
    result["measurement_mode"] = "one-shot"
    return result


def one_shot_matrix() -> list[dict[str, Any]]:
    """Bound the promoted one-shot region and its immediate fallbacks."""
    result: list[dict[str, Any]] = []
    boundary_bytes = (
        1, 31, 32, 33, 63, 64, 65,
        127, 128, 129, 255, 256, 257,
        511, 512, 513, 514,
    )
    for k in (17, 128):
        for shard_bytes in boundary_bytes:
            result.append(one_shot_cell(
                f"oneshot-k{k}-r{k}-b{shard_bytes}-l8",
                k, k, shard_bytes, 8, "target", True, False))
        for shard_bytes in (1, 33, 65, 257, 514):
            result.append(one_shot_cell(
                f"oneshot-k{k}-r{k}-b{shard_bytes}-l1",
                k, k, shard_bytes, 1, "target", True, False))

    for k in (32, 66, 96, 127):
        for losses in (1, 8):
            result.append(one_shot_cell(
                f"oneshot-edge-k{k}-r{k}-b514-l{losses}",
                k, k, 514, losses, "target", True, False))

    # Exact Leopard1 rejects R > K, so keep the authoritative three-way
    # campaign to geometries both public APIs can represent.  R > K remains
    # covered by same-source correctness and directional screens.
    for k, r in ((32, 17), (128, 65)):
        for losses in (1, 8):
            for shard_bytes in (33, 514):
                result.append(one_shot_cell(
                    f"oneshot-asym-k{k}-r{r}-b{shard_bytes}-l{losses}",
                    k, r, shard_bytes, losses, "target", True, False))

    for k in (17, 32, 128):
        result.append(one_shot_cell(
            f"neighbor-cutoff-k{k}-r{k}-b515-l8",
            k, k, 515, 8, "neighbor", False, False))
    result.extend((
        one_shot_cell("neighbor-k65-r65-b33-l8",
                      65, 65, 33, 8, "neighbor", False, False),
        one_shot_cell("neighbor-k65-r65-b514-l8",
                      65, 65, 514, 8, "neighbor", False, False),
        one_shot_cell("neighbor-k16-r17-b33-l8",
                      16, 17, 33, 8, "neighbor", False, False),
        one_shot_cell("neighbor-k17-r17-b33-l9",
                      17, 17, 33, 9, "neighbor", False, False),
        one_shot_cell("neighbor-k128-r128-b514-l9",
                      128, 128, 514, 9, "neighbor", False, False),
    ))
    require(len({item["id"] for item in result}) == len(result),
            "one-shot matrix contains duplicate identifiers")
    return result


def matrix(mode: str = "reusable") -> list[dict[str, Any]]:
    require(mode in ("reusable", "one-shot"),
            f"unsupported measurement mode: {mode}")
    return one_shot_matrix() if mode == "one-shot" else reusable_matrix()


def manifest_schema(mode: str) -> str:
    return ONE_SHOT_MANIFEST_SCHEMA if mode == "one-shot" \
        else MANIFEST_SCHEMA


def cell_schema(item: Mapping[str, Any]) -> str:
    return ONE_SHOT_CELL_SCHEMA \
        if item.get("measurement_mode") == "one-shot" else CELL_SCHEMA


def summary_schema(mode: str) -> str:
    return ONE_SHOT_SUMMARY_SCHEMA if mode == "one-shot" \
        else SUMMARY_SCHEMA


def write_atomic_exclusive(path: Path, value: object) -> None:
    payload = json.dumps(
        value, indent=2, sort_keys=True, allow_nan=False
    ).encode("utf-8") + b"\n"
    temporary = path.with_name(
        f".{path.name}.tmp-{os.getpid()}-{time.monotonic_ns()}")
    descriptor = os.open(
        temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    try:
        view = memoryview(payload)
        while view:
            written = os.write(descriptor, view)
            require(written > 0, f"short write: {temporary}")
            view = view[written:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    try:
        os.link(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def load_json(path: Path) -> dict[str, Any]:
    status = path.stat()
    require(0 < status.st_size <= MAX_JSON_BYTES and path.is_file(),
            f"invalid JSON artifact size: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise EvidenceError(f"invalid JSON artifact {path}: {error}") from error
    require(isinstance(value, dict), f"JSON artifact is not an object: {path}")
    return value


def retained_file_identity(
    path: Path,
    maximum_bytes: int,
    require_nonempty: bool,
) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    status = resolved.stat()
    require(
        resolved.is_file() and status.st_size <= maximum_bytes and
        (status.st_size > 0 or not require_nonempty),
        f"invalid retained file: {resolved}")
    return {
        "path": str(resolved),
        "size": status.st_size,
        "mode": status.st_mode & 0o777,
        "sha256": T8.sha256(resolved),
    }


def allocate_raw_run_directory(
    raw_root: Path,
    identifier: str,
) -> tuple[Path, list[str]]:
    """Retain incomplete attempts while giving each resume fresh filenames."""
    cell_root = raw_root / identifier
    cell_root.mkdir(mode=0o700, exist_ok=True)
    prior = []
    for entry in sorted(cell_root.iterdir(), key=lambda value: value.name):
        require(entry.is_dir() and not entry.is_symlink(),
                f"unexpected raw cell entry: {entry}")
        prior.append(str(entry.relative_to(raw_root.parent)))
    for _ in range(100):
        directory = cell_root / (
            f"run-{os.getpid()}-{time.monotonic_ns()}")
        try:
            directory.mkdir(mode=0o700)
            return directory, prior
        except FileExistsError:
            continue
    raise EvidenceError(f"cannot allocate raw run directory: {identifier}")


def freeze_executable(origin: Path, directory: Path, name: str) -> dict[str, Any]:
    resolved = origin.resolve(strict=True)
    before = T8.file_identity(resolved)
    target = directory / name
    require(not target.exists(), f"frozen executable already exists: {target}")
    shutil.copyfile(resolved, target)
    os.chmod(target, 0o555)
    after = T8.file_identity(target)
    require(
        before["size"] == after["size"] and
        before["sha256"] == after["sha256"],
        f"frozen executable differs from origin: {name}")
    return {"origin": before, "frozen": after}


def metric(result: Mapping[str, Any], name: str, field: str) -> float:
    metrics = result.get("metrics")
    require(isinstance(metrics, dict), "benchmark metrics are absent")
    value = metrics.get(name)
    require(isinstance(value, dict), f"benchmark metric is absent: {name}")
    number = value.get(field)
    require(isinstance(number, (int, float)) and
            not isinstance(number, bool) and
            math.isfinite(float(number)) and float(number) >= 0,
            f"benchmark metric is invalid: {name}.{field}")
    return float(number)


def expected_parameters(
    implementation: str,
    item: Mapping[str, Any],
    iterations: int,
    warmup: int,
    reuse: int,
) -> dict[str, Any]:
    physical = int(item["main_physical_bytes"]) \
        if implementation == "main" else int(item["bytes"])
    expected = {
        "K": item["K"],
        "R": item["R"],
        "shard_bytes": physical,
        "loss_count": item["loss"],
        "batch": 1,
        "reuse": reuse,
        "iterations": iterations,
        "warmup": warmup,
        "thread_count": 1,
        "seed": item["seed"],
    }
    if implementation == "main" and physical != item["bytes"]:
        expected["logical_shard_bytes"] = item["bytes"]
    if (implementation != "main" and
            item.get("measurement_mode") == "one-shot"):
        expected["measure_one_shot_decode"] = True
    return expected


def validate_result(
    implementation: str,
    result: object,
    item: Mapping[str, Any],
    source_commit: str,
    source_tree: str,
    iterations: int,
    warmup: int,
    reuse: int,
) -> dict[str, Any]:
    require(isinstance(result, dict), "benchmark output is not an object")
    parameters = result.get("parameters")
    resolved = result.get("resolved")
    correctness = result.get("correctness")
    digests = result.get("workload_digests")
    require(all(isinstance(value, dict) for value in (
                parameters, resolved, correctness, digests)),
            "benchmark output is incomplete")
    for name, expected in expected_parameters(
            implementation, item, iterations, warmup, reuse).items():
        require(parameters.get(name) == expected,
                f"{implementation} parameter changed: {name}")
    require(
        resolved.get("profile") == "legacy_high_v1" and
        resolved.get("field") == "gf8" and
        resolved.get("thread_count") == 1,
        f"{implementation} resolved a different code identity")
    require(
        digests.get("algorithm") == "fnv1a64" and
        all(isinstance(digests.get(name), str) and len(digests[name]) == 16
            for name in (
                "original_data", "transmitted_parity",
                "recovered_originals")),
        f"{implementation} workload digests are incomplete")
    if implementation == "main":
        padded = int(item["main_physical_bytes"]) != int(item["bytes"])
        require(
            result.get("schema") ==
                ("leopard-main-benchmark-v2" if padded
                 else "leopard-main-benchmark-v1") and
            correctness.get("round_trip") is True,
            "exact-main identity or round trip failed")
        build = result.get("build")
        require(isinstance(build, dict) and
                build.get("main_source_commit") == MAIN_COMMIT,
                "exact-main source commit changed")
        if padded:
            require(
                resolved.get("padded_application_bytes") is True and
                correctness.get("logical_prefix_fingerprinted") is True,
                "exact-main logical-prefix semantics changed")
    elif item.get("measurement_mode") == "one-shot":
        expected_enabled = implementation == "candidate"
        build = result.get("build")
        require(
            result.get("schema") == "leopard2-benchmark-v8" and
            resolved.get("backend") == "avx2" and
            correctness.get("leopard2_round_trip") is True and
            parameters.get("skip_legacy") is True and
            parameters.get("retain_samples") is True and
            parameters.get("measure_one_shot_decode") is True and
            isinstance(build, dict) and
            build.get("one_shot_equal_rounded_direct_enabled")
                is expected_enabled,
            f"{implementation} one-shot build identity changed")
    else:
        expected_enabled = implementation == "candidate"
        expected_direct = bool(item[
            f"{implementation}_direct"])
        build = result.get("build")
        require(
            result.get("schema") == "leopard2-benchmark-v7" and
            resolved.get("backend") == "avx2" and
            correctness.get("leopard2_round_trip") is True and
            parameters.get("attest_source") is True and
            parameters.get("report_decode_path") is True and
            parameters.get("report_direct_executor") is True and
            isinstance(build, dict) and
            build.get("source_commit") == source_commit and
            build.get("source_tree") == source_tree and
            build.get("source_tracked_dirty") is False and
            build.get("equal_rounded_multi_loss_enabled") is expected_enabled,
            f"{implementation} build or path-attestation identity changed")
        if expected_direct:
            require(
                resolved.get("selected_decode_path") == "direct" and
                resolved.get("selected_direct_executor") == "source_major",
                f"{implementation} did not select source-major direct repair")
        else:
            require(
                resolved.get("selected_decode_path") != "direct" and
                resolved.get("selected_direct_executor") == "none",
                f"{implementation} unexpectedly selected direct repair")
    if implementation == "main":
        # Exact Leopard1 constructs decoder state inside every legacy API call
        # and reports one inseparable decode_including_setup value. Preserve
        # that public-call cost as both execution and first use; Leopard2
        # reports reusable plan setup separately below.
        setup = 0.0
        execution = metric(
            result, "decode_including_setup",
            "median_us_per_batch_call")
    elif item.get("measurement_mode") == "one-shot":
        setup = 0.0
        execution = metric(
            result, "one_shot_decode_including_setup",
            "median_us_per_batch_call")
    else:
        setup = metric(result, "decode_plan_setup", "median_us")
        execution = metric(
            result, "decode_execution", "median_us_per_batch_call")
    return {
        "decode_us": execution,
        "plan_setup_us": setup,
        "first_use_us": execution + setup,
        "amortized_us": execution + setup / reuse,
        "digests": {
            name: digests[name] for name in (
                "original_data", "transmitted_parity",
                "recovered_originals")
        },
    }


def benchmark_command(
    implementation: str,
    executable: Path,
    item: Mapping[str, Any],
    iterations: int,
    warmup: int,
    reuse: int,
) -> list[str]:
    physical = int(item["main_physical_bytes"]) \
        if implementation == "main" else int(item["bytes"])
    command = [
        "/usr/bin/prlimit", "--as=268435456",
        str(executable),
        "--k", str(item["K"]),
        "--r", str(item["R"]),
        "--bytes", str(physical),
        "--loss", str(item["loss"]),
        "--batch", "1",
        "--reuse", str(reuse),
        "--iterations", str(iterations),
        "--warmup", str(warmup),
        "--threads", "1",
        "--seed", str(item["seed"]),
    ]
    if implementation == "main":
        if physical != item["bytes"]:
            command.extend(("--logical-bytes", str(item["bytes"])))
    elif item.get("measurement_mode") == "one-shot":
        command.extend((
            "--profile", "high",
            "--field", "gf8",
            "--backend", "avx2",
            "--skip-legacy",
            "--retain-samples",
            "--measure-one-shot-decode",
        ))
    else:
        command.extend((
            "--profile", "high",
            "--field", "gf8",
            "--backend", "avx2",
            "--skip-legacy",
            "--retain-samples",
            "--report-direct-executor",
            "--attest-source",
        ))
    command.extend(("--json", "-"))
    return command


def attest_one_shot_source(
    executable: Path,
    identity: Mapping[str, Any],
    source_commit: str,
    source_tree: str,
) -> dict[str, Any]:
    """Bind a schema-v8 binary to its clean schema-v5 source identity."""
    require(T8.sha256(executable) == identity["sha256"],
            "one-shot executable changed before source attestation")
    command = [
        "/usr/bin/prlimit", "--as=268435456",
        str(executable),
        "--k", "32", "--r", "32",
        "--bytes", "64", "--loss", "8",
        "--batch", "1", "--reuse", "1",
        "--iterations", "1", "--warmup", "0",
        "--threads", "1", "--seed", "150",
        "--profile", "high", "--field", "gf8",
        "--backend", "avx2", "--skip-legacy",
        "--attest-source", "--json", "-",
    ]
    completed = subprocess.run(
        command, env=CHILD_ENVIRONMENT,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        timeout=30.0, check=False)
    require(
        completed.returncode == 0 and
        0 < len(completed.stdout) <= (8 << 20) and
        len(completed.stderr) <= (1 << 20),
        "one-shot source attestation process failed")
    try:
        result = json.loads(completed.stdout.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise EvidenceError(
            f"invalid one-shot source attestation: {error}") from error
    require(isinstance(result, dict),
            "one-shot source attestation is not an object")
    build = result.get("build")
    parameters = result.get("parameters")
    resolved = result.get("resolved")
    correctness = result.get("correctness")
    digests = result.get("workload_digests")
    require(
        result.get("schema") == "leopard2-benchmark-v5" and
        isinstance(build, dict) and
        build.get("source_commit") == source_commit and
        build.get("source_tree") == source_tree and
        build.get("source_tracked_dirty") is False and
        isinstance(parameters, dict) and
        parameters.get("attest_source") is True and
        isinstance(resolved, dict) and
        resolved.get("profile") == "legacy_high_v1" and
        resolved.get("field") == "gf8" and
        resolved.get("backend") == "avx2" and
        isinstance(correctness, dict) and
        correctness.get("leopard2_round_trip") is True and
        isinstance(digests, dict) and
        digests.get("algorithm") == "fnv1a64" and
        all(isinstance(digests.get(name), str) and
            len(digests[name]) == 16 for name in
            ("original_data", "transmitted_parity", "recovered_originals")),
        "one-shot source attestation identity changed")
    require(T8.sha256(executable) == identity["sha256"],
            "one-shot executable changed after source attestation")
    return {
        "command": command,
        "executable_sha256": identity["sha256"],
        "source_commit": build["source_commit"],
        "source_tree": build["source_tree"],
        "source_tracked_dirty": build["source_tracked_dirty"],
        "resolved": {
            name: resolved[name] for name in
            ("profile", "field", "backend")
        },
        "workload_digests": {
            name: digests[name] for name in
            ("original_data", "transmitted_parity", "recovered_originals")
        },
        "stdout_sha256": hashlib.sha256(completed.stdout).hexdigest(),
        "stderr_sha256": hashlib.sha256(completed.stderr).hexdigest(),
    }


def target_endpoint_tolerance_ns(child_runtime_ns: int) -> int:
    require(type(child_runtime_ns) is int and child_runtime_ns >= 0,
            "child runtime is invalid")
    relative = (
        child_runtime_ns * TARGET_ENDPOINT_RELATIVE_TOLERANCE_PPM + 999_999
    ) // 1_000_000
    return max(TARGET_ENDPOINT_ABSOLUTE_TOLERANCE_NS, relative)


def sibling_runtime_tolerance_ns(child_runtime_ns: int) -> int:
    require(type(child_runtime_ns) is int and child_runtime_ns >= 0,
            "child runtime is invalid")
    relative = (
        child_runtime_ns * SIBLING_RELATIVE_TOLERANCE_PPM + 999_999
    ) // 1_000_000
    return max(SIBLING_ABSOLUTE_TOLERANCE_NS, relative)


def isolation_accepted(gated: Mapping[str, Any]) -> bool:
    target = gated["target_runtime"]
    return (
        abs(int(target["signed_difference_ns"])) <=
            target_endpoint_tolerance_ns(int(target["child_cpu_time_ns"])) and
        gated["wait4_crosscheck"]["accepted"] and
        gated["target_interrupts"]["accepted"] and
        int(gated["sibling_runtime"]["scheduler_delta_ns"]) <=
            sibling_runtime_tolerance_ns(
                int(target["child_cpu_time_ns"])) and
        gated["sibling_interrupts"]["accepted"]
    )


def run_one(
    implementation: str,
    executable: Path,
    identity: Mapping[str, Any],
    item: Mapping[str, Any],
    source_commit: str,
    source_tree: str,
    iterations: int,
    warmup: int,
    reuse: int,
    cpu: int,
    sibling: int,
    raw_directory: Path,
    round_index: int,
    slot_index: int,
    lock_descriptor: int,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    rejected = []
    require(T8.sha256(executable) == identity["sha256"],
            f"{implementation} executable changed before invocation")
    command = benchmark_command(
        implementation, executable, item, iterations, warmup, reuse)
    for attempt in range(MAX_ATTEMPTS):
        stem = (
            f"round-{round_index:02d}-slot-{slot_index}-"
            f"{implementation}-attempt-{attempt:02d}"
        )
        stdout_path = raw_directory / f"{stem}.stdout.json"
        stderr_path = raw_directory / f"{stem}.stderr.txt"
        gated = GATE.run_gated_benchmark(
            command, cpu, sibling, 120.0, stdout_path, stderr_path,
            campaign_lock_descriptor=lock_descriptor)
        require(T8.sha256(executable) == identity["sha256"],
                f"{implementation} executable changed after invocation")
        record = {
            "implementation": implementation,
            "attempt": attempt,
            "command": command,
            "stdout": retained_file_identity(
                stdout_path, 8 << 20, True),
            "stderr": retained_file_identity(
                stderr_path, 1 << 20, False),
            "gated": gated,
            "accepted": False,
        }
        if gated["return_code"] == 0 and isolation_accepted(gated):
            result = load_json(stdout_path)
            normalized = validate_result(
                implementation, result, item, source_commit, source_tree,
                iterations, warmup, reuse)
            record["result"] = result
            record["normalized"] = normalized
            record["accepted"] = True
            envelope_path = raw_directory / f"{stem}.envelope.json"
            write_atomic_exclusive(envelope_path, record)
            record["envelope"] = retained_file_identity(
                envelope_path, MAX_JSON_BYTES, True)
            return record, rejected
        envelope_path = raw_directory / f"{stem}.envelope.json"
        write_atomic_exclusive(envelope_path, record)
        record["envelope"] = retained_file_identity(
            envelope_path, MAX_JSON_BYTES, True)
        rejected.append(record)
    raise EvidenceError(
        f"{item['id']} {implementation} failed isolation "
        f"for {MAX_ATTEMPTS} attempts")


def log_mean_ratio(
    invocations: Sequence[Mapping[str, Any]],
    numerator: str,
    denominator: str,
    metric_name: str,
) -> float:
    numerator_values = [
        float(item["normalized"][metric_name])
        for item in invocations if item["implementation"] == numerator
    ]
    denominator_values = [
        float(item["normalized"][metric_name])
        for item in invocations if item["implementation"] == denominator
    ]
    require(
        len(numerator_values) == len(denominator_values) and
        len(numerator_values) in (1, 2),
        "round does not contain a balanced comparison")
    return (
        statistics.mean(math.log(value) for value in numerator_values) -
        statistics.mean(math.log(value) for value in denominator_values)
    )


def confidence(log_ratios: Sequence[float]) -> dict[str, float]:
    require(len(log_ratios) in T95, "unsupported confidence sample count")
    center = statistics.mean(log_ratios)
    half = T95[len(log_ratios)] * statistics.stdev(log_ratios) / \
        math.sqrt(len(log_ratios))
    return {
        "ratio": math.exp(center),
        "lower": math.exp(center - half),
        "upper": math.exp(center + half),
    }


def analyze_cell(
    item: Mapping[str, Any],
    rounds: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    control_logs = [
        log_mean_ratio(
            round_value["invocations"], "control", "candidate", "decode_us")
        for round_value in rounds
    ]
    result = {
        "cell": dict(item),
        "control_over_candidate": confidence(control_logs),
    }
    if item["role"] == "target":
        main_logs = [
            log_mean_ratio(
                round_value["invocations"], "main", "candidate", "decode_us")
            for round_value in rounds
        ]
        first_use_logs = [
            log_mean_ratio(
                round_value["invocations"], "main", "candidate",
                "first_use_us")
            for round_value in rounds
        ]
        result["main_over_candidate"] = confidence(main_logs)
        result["main_over_candidate_first_use"] = confidence(first_use_logs)
    return result


def validate_cell_record(
    value: object,
    expected: Mapping[str, Any],
    manifest_sha256: str,
    target_rounds: int,
) -> dict[str, Any]:
    require(isinstance(value, dict) and
            value.get("schema") == cell_schema(expected),
            "cell record schema changed")
    require(
        value.get("manifest_sha256") == manifest_sha256 and
        value.get("cell") == expected and
        isinstance(value.get("rounds"), list) and
        isinstance(value.get("analysis"), dict),
        "cell record identity changed")
    raw_run_directory = value.get("raw_run_directory")
    prior_raw_directories = value.get("prior_incomplete_raw_directories")
    raw_prefix = f"raw/{expected['id']}/run-"
    require(
        isinstance(raw_run_directory, str) and
        raw_run_directory.startswith(raw_prefix) and
        isinstance(prior_raw_directories, list) and
        all(isinstance(path, str) and path.startswith(raw_prefix)
            for path in prior_raw_directories) and
        prior_raw_directories == sorted(prior_raw_directories) and
        len(set(prior_raw_directories)) == len(prior_raw_directories) and
        raw_run_directory not in prior_raw_directories,
        "cell raw-run identity changed")
    require(value.get("target_rounds") == target_rounds,
            "cell record target-round policy changed")
    expected_rounds = 3 if expected["role"] == "neighbor" else target_rounds
    require(len(value["rounds"]) == expected_rounds,
            "cell record has the wrong round count")
    for round_value in value["rounds"]:
        require(
            isinstance(round_value, dict) and
            isinstance(round_value.get("invocations"), list) and
            all(invocation.get("accepted") is True and
                invocation["gated"]["return_code"] == 0 and
                isolation_accepted(invocation["gated"])
                for invocation in round_value["invocations"]),
            "cell record contains an unaccepted invocation")
    require(
        isinstance(value.get("rejected_attempts"), list) and
        all(record.get("accepted") is False and
            (record["gated"]["return_code"] != 0 or
             not isolation_accepted(record["gated"]))
            for record in value["rejected_attempts"]),
        "cell record contains an accepted rejected-attempt")
    require(
        value["analysis"] == analyze_cell(expected, value["rounds"]),
        "cell analysis is not derived from retained rounds")
    return value


def aggregate(
    analyses: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    require(bool(analyses), "campaign has no analyses")
    targets = [item for item in analyses if item["cell"]["role"] == "target"]
    neighbors = [
        item for item in analyses if item["cell"]["role"] == "neighbor"
    ]
    target_failures = [
        item["cell"]["id"] for item in targets
        if item["control_over_candidate"]["lower"] < TARGET_CONTROL_FLOOR or
        item["main_over_candidate"]["lower"] < TARGET_MAIN_FLOOR
    ]
    neighbor_failures = [
        item["cell"]["id"] for item in neighbors
        if item["control_over_candidate"]["upper"] < NEIGHBOR_FLOOR
    ]
    target_control = [
        item["control_over_candidate"]["ratio"] for item in targets]
    target_main = [
        item["main_over_candidate"]["ratio"] for item in targets]
    target_main_first_use = [
        item["main_over_candidate_first_use"]["ratio"] for item in targets]
    target_control_lower = [
        item["control_over_candidate"]["lower"] for item in targets]
    target_main_lower = [
        item["main_over_candidate"]["lower"] for item in targets]
    target_main_first_use_lower = [
        item["main_over_candidate_first_use"]["lower"] for item in targets]
    neighbor_upper = [
        item["control_over_candidate"]["upper"] for item in neighbors]

    def geomean(values: Sequence[float]) -> float | None:
        return math.exp(statistics.mean(map(math.log, values))) \
            if values else None

    return {
        "target_count": len(targets),
        "neighbor_count": len(neighbors),
        "target_failures": target_failures,
        "neighbor_failures": neighbor_failures,
        "accepted": not target_failures and not neighbor_failures,
        "target_control_geomean": geomean(target_control),
        "target_main_geomean": geomean(target_main),
        "target_main_first_use_geomean": geomean(target_main_first_use),
        "minimum_target_control_lower":
            min(target_control_lower) if targets else None,
        "minimum_target_main_lower":
            min(target_main_lower) if targets else None,
        "minimum_target_main_first_use_lower":
            min(target_main_first_use_lower) if targets else None,
        "minimum_neighbor_upper":
            min(neighbor_upper) if neighbors else None,
    }


def run(options: argparse.Namespace) -> int:
    require(options.target_rounds in (3, 9),
            "--target-rounds must be 3 or 9")
    require(options.iterations > 0 and options.warmup > 0 and
            options.reuse > 0, "timing dimensions must be positive")
    allowed = set(os.sched_getaffinity(0))
    require(options.cpu not in allowed and options.sibling not in allowed,
            "coordinator affinity must exclude the measured CPU pair")
    sibling_text = Path(
        f"/sys/devices/system/cpu/cpu{options.cpu}/topology/"
        "thread_siblings_list").read_text(encoding="ascii")
    require(PAIR.parse_cpu_list(sibling_text) ==
            {options.cpu, options.sibling},
            "requested CPUs are not one SMT pair")
    cells = matrix(options.mode)
    if options.cell_id:
        by_id = {item["id"]: item for item in cells}
        unknown = sorted(set(options.cell_id) - set(by_id))
        require(not unknown, f"unknown cell IDs: {unknown}")
        cells = [by_id[identifier] for identifier in options.cell_id]
    require(cells, "campaign has no selected cells")
    require(options.output.is_absolute(),
            "output directory must be an absolute lane-owned path")
    if options.resume:
        require(options.output.is_dir(), "resume directory is absent")
    else:
        require(not options.output.exists(), "output already exists")
        options.output.mkdir(mode=0o700, parents=True)
    artifacts = options.output / "artifacts"
    if options.resume:
        require(artifacts.is_dir(), "resume artifacts are absent")
        identities = load_json(options.output / "identities.json")
    else:
        artifacts.mkdir(mode=0o700)
        frozen = {
            name: freeze_executable(path, artifacts, name)
            for name, path in (
                ("candidate", options.candidate),
                ("control", options.control),
                ("main", options.main),
            )
        }
        identities = {
            name: record["frozen"] for name, record in frozen.items()
        }
        require(len({value["sha256"] for value in identities.values()}) == 3,
                "candidate, control, and main binaries are not distinct")
        executable_sections = {
            name: T8.executable_sections_identity(
                Path(value["path"]))
            for name, value in identities.items()
            if name in ("candidate", "control")
        }
        require(
            executable_sections["candidate"]["sections"] ==
            executable_sections["control"]["sections"],
            "candidate/control executable instruction sections differ")
        identities = {
            "frozen": frozen,
            "executables": identities,
            "executable_sections": executable_sections,
        }
        write_atomic_exclusive(options.output / "identities.json", identities)
    executable_identities = identities["executables"]
    for identity in executable_identities.values():
        require(
            T8.sha256(Path(identity["path"])) == identity["sha256"] and
            (Path(identity["path"]).stat().st_mode & 0o777) == 0o555,
            "frozen executable changed before campaign")
    source_attestations: dict[str, Any] | None = None
    if options.mode == "one-shot":
        if options.resume:
            prior_manifest = load_json(options.output / "manifest.json")
            prior_attestations = prior_manifest.get("source_attestations")
            require(isinstance(prior_attestations, dict),
                    "resume source attestations are absent")
            source_attestations = prior_attestations
        else:
            source_attestations = {
                name: attest_one_shot_source(
                    Path(executable_identities[name]["path"]),
                    executable_identities[name],
                    options.source_commit, options.source_tree)
                for name in ("candidate", "control")
            }
            require(len({
                canonical_bytes(value["workload_digests"])
                for value in source_attestations.values()
            }) == 1, "one-shot source-attestation workloads differ")
    manifest = {
        "schema": manifest_schema(options.mode),
        "measurement_mode": options.mode,
        "source_commit": options.source_commit,
        "source_tree": options.source_tree,
        "main_commit": MAIN_COMMIT,
        "target_rounds": options.target_rounds,
        "neighbor_rounds": 3,
        "iterations": options.iterations,
        "warmup": options.warmup,
        "reuse": options.reuse,
        "cpu": options.cpu,
        "sibling": options.sibling,
        "coordinator_affinity": sorted(allowed),
        "runner": T8.regular_file_identity(Path(__file__)),
        "gate_runner": T8.regular_file_identity(
            DIRECTORY / "run_small_direct_abba.py"),
        "main_support": T8.regular_file_identity(
            DIRECTORY.parent / "main_compare" / "run_abba.py"),
        "identities": identities,
        "isolation_policy": {
            "target_scheduler_minus_child_runtime_effective_max":
                "max(absolute_floor_ns, child_runtime_ns * relative_ppm / 1e6)",
            "target_scheduler_minus_child_runtime_absolute_floor_ns":
                TARGET_ENDPOINT_ABSOLUTE_TOLERANCE_NS,
            "target_scheduler_minus_child_runtime_relative_ppm":
                TARGET_ENDPOINT_RELATIVE_TOLERANCE_PPM,
            "wait4_crosscheck_max_ns":
                GATE.RUSAGE_CROSSCHECK_TOLERANCE_NS,
            "target_rejected_interrupt_fields":
                ["irq", "softirq", "steal", "guest", "guest_nice"],
            "reserved_sibling_scheduler_runtime_effective_max":
                "max(absolute_floor_ns, child_runtime_ns * relative_ppm / 1e6)",
            "reserved_sibling_scheduler_runtime_absolute_floor_ns":
                SIBLING_ABSOLUTE_TOLERANCE_NS,
            "reserved_sibling_scheduler_runtime_relative_ppm":
                SIBLING_RELATIVE_TOLERANCE_PPM,
            "maximum_attempts": MAX_ATTEMPTS,
            "retry_trigger": "objective isolation failure only",
        },
        "cells": cells,
    }
    if source_attestations is not None:
        manifest["source_attestations"] = source_attestations
    manifest_sha256 = digest_object(manifest)
    manifest["digest"] = manifest_sha256
    manifest_path = options.output / "manifest.json"
    if options.resume:
        require(load_json(manifest_path) == manifest,
                "resume manifest changed")
    else:
        write_atomic_exclusive(manifest_path, manifest)
    cells_directory = options.output / "cells"
    raw_root = options.output / "raw"
    if options.resume:
        require(cells_directory.is_dir() and raw_root.is_dir(),
                "resume cell/raw directories are absent")
    else:
        cells_directory.mkdir(mode=0o700)
        raw_root.mkdir(mode=0o700)
    lock_descriptor = os.open(
        LOCK_PATH, os.O_RDWR | os.O_CREAT | os.O_CLOEXEC, 0o600)
    import fcntl
    fcntl.flock(lock_descriptor, fcntl.LOCK_EX)
    analyses = []
    rejected_attempt_count = 0
    try:
        with PAIR.StableLeaseAnchor(), \
                PAIR.PairLease(options.cpu, options.sibling) as pair_lease:
            for index, item in enumerate(cells):
                cell_path = cells_directory / f"{item['id']}.json"
                if cell_path.exists():
                    require(options.resume,
                            f"cell already exists: {item['id']}")
                    value = validate_cell_record(
                        load_json(cell_path), item, manifest_sha256,
                        options.target_rounds)
                    analyses.append(value["analysis"])
                    rejected_attempt_count += len(
                        value["rejected_attempts"])
                    print(
                        f"{index + 1}/{len(cells)} {item['id']} resumed",
                        file=sys.stderr, flush=True)
                    continue
                raw_directory, prior_raw_directories = \
                    allocate_raw_run_directory(raw_root, item["id"])
                orders = TARGET_ORDER * (
                    options.target_rounds // len(TARGET_ORDER)
                ) if item["role"] == "target" else NEIGHBOR_ORDER
                rounds = []
                rejected = []
                for round_index, order in enumerate(orders):
                    invocations = []
                    for slot_index, implementation in enumerate(order):
                        record, rejected_records = run_one(
                            implementation,
                            Path(executable_identities[
                                implementation]["path"]),
                            executable_identities[implementation],
                            item,
                            options.source_commit,
                            options.source_tree,
                            options.iterations,
                            options.warmup,
                            options.reuse,
                            options.cpu,
                            options.sibling,
                            raw_directory,
                            round_index,
                            slot_index,
                            lock_descriptor)
                        invocations.append(record)
                        rejected.extend(rejected_records)
                    rounds.append({
                        "round": round_index,
                        "order": list(order),
                        "invocations": invocations,
                    })
                analysis = analyze_cell(item, rounds)
                value = {
                    "schema": cell_schema(item),
                    "manifest_sha256": manifest_sha256,
                    "pair_lease": pair_lease,
                    "cell": item,
                    "target_rounds": options.target_rounds,
                    "raw_run_directory": str(
                        raw_directory.relative_to(options.output)),
                    "prior_incomplete_raw_directories":
                        prior_raw_directories,
                    "rounds": rounds,
                    "rejected_attempts": rejected,
                    "analysis": analysis,
                }
                write_atomic_exclusive(cell_path, value)
                validated = validate_cell_record(
                    load_json(cell_path), item, manifest_sha256,
                    options.target_rounds)
                analyses.append(validated["analysis"])
                rejected_attempt_count += len(rejected)
                print(
                    f"{index + 1}/{len(cells)} {item['id']}",
                    file=sys.stderr, flush=True)
    finally:
        os.close(lock_descriptor)
    for identity in executable_identities.values():
        require(T8.sha256(Path(identity["path"])) == identity["sha256"],
                "frozen executable changed after campaign")
    all_digests_matched = all(
        len({
            tuple(sorted(invocation["normalized"]["digests"].items()))
            for invocation in round_value["invocations"]
        }) == 1
        for item in cells
        for round_value in load_json(
            cells_directory / f"{item['id']}.json")["rounds"]
    )
    summary = {
        "schema": summary_schema(options.mode),
        "manifest_sha256": manifest_sha256,
        "cell_count": len(cells),
        "accepted_process_count": sum(
            len(round_value["invocations"])
            for item in cells
            for round_value in load_json(
                cells_directory / f"{item['id']}.json")["rounds"]
        ),
        "rejected_isolation_attempt_count": rejected_attempt_count,
        "all_digests_matched": all_digests_matched,
        "analysis": aggregate(analyses),
    }
    summary["digest"] = digest_object(summary)
    write_atomic_exclusive(options.output / "summary.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["analysis"]["accepted"] and \
        all_digests_matched else 2


def verify(options: argparse.Namespace) -> int:
    root = options.output.resolve(strict=True)
    manifest = load_json(root / "manifest.json")
    mode = manifest.get("measurement_mode", "reusable")
    require(mode in ("reusable", "one-shot") and
            manifest.get("schema") == manifest_schema(mode),
            "manifest schema changed")
    expected_digest = manifest.pop("digest", None)
    require(expected_digest == digest_object(manifest),
            "manifest digest changed")
    manifest["digest"] = expected_digest
    cells = manifest["cells"]
    analyses = []
    for item in cells:
        value = validate_cell_record(
            load_json(root / "cells" / f"{item['id']}.json"),
            item, expected_digest, int(manifest["target_rounds"]))
        analyses.append(value["analysis"])
    summary = load_json(root / "summary.json")
    require(summary.get("schema") == summary_schema(mode),
            "summary schema changed")
    summary_digest = summary.pop("digest", None)
    require(summary_digest == digest_object(summary),
            "summary digest changed")
    require(summary["analysis"] == aggregate(analyses),
            "summary analysis changed")
    require(summary.get("all_digests_matched") is True,
            "cross-implementation workload digests differ")
    identities = manifest["identities"]["executables"]
    for identity in identities.values():
        require(T8.sha256(Path(identity["path"])) == identity["sha256"],
                "frozen executable changed")
    if mode == "one-shot":
        require(
            T8.regular_file_identity(Path(__file__)) == manifest["runner"] and
            T8.regular_file_identity(
                DIRECTORY / "run_small_direct_abba.py") ==
                    manifest["gate_runner"] and
            T8.regular_file_identity(
                DIRECTORY.parent / "main_compare" / "run_abba.py") ==
                    manifest["main_support"],
            "one-shot evidence support identity changed")
        executable_sections = manifest["identities"].get(
            "executable_sections")
        require(isinstance(executable_sections, dict),
                "one-shot executable-section evidence is absent")
        recomputed_sections = {
            name: T8.executable_sections_identity(
                Path(identities[name]["path"]))
            for name in ("candidate", "control")
        }
        require(
            recomputed_sections == executable_sections and
            executable_sections["candidate"]["sections"] ==
                executable_sections["control"]["sections"],
            "one-shot executable sections changed")
        attestations = manifest.get("source_attestations")
        require(isinstance(attestations, dict) and
                set(attestations) == {"candidate", "control"},
                "one-shot source attestations changed")
        for name, attestation in attestations.items():
            require(
                isinstance(attestation, dict) and
                attestation.get("executable_sha256") ==
                    identities[name]["sha256"] and
                attestation.get("source_commit") ==
                    manifest["source_commit"] and
                attestation.get("source_tree") ==
                    manifest["source_tree"] and
                attestation.get("source_tracked_dirty") is False,
                f"{name} one-shot source attestation changed")
        require(len({
            canonical_bytes(attestation.get("workload_digests"))
            for attestation in attestations.values()
        }) == 1, "one-shot source-attestation workloads differ")
    print("equal-rounded evidence verified")
    return 0


def self_test() -> int:
    cells = matrix()
    require(len(cells) == 52, "matrix cell count changed")
    require(sum(item["role"] == "target" for item in cells) == 47,
            "target cell count changed")
    require(sum(item["role"] == "neighbor" for item in cells) == 5,
            "neighbor cell count changed")
    require(digest_object(cells) ==
            "96246fe43f13caf70bf615943782bcd4d53f6295b62308e04c09b43e7a7f8ff3",
            "matrix digest changed")
    one_shot_cells = matrix("one-shot")
    require(
        len(one_shot_cells) == 68 and
        sum(item["role"] == "target" for item in one_shot_cells) == 60 and
        sum(item["role"] == "neighbor" for item in one_shot_cells) == 8 and
        all(item.get("measurement_mode") == "one-shot"
            for item in one_shot_cells),
        "one-shot matrix shape changed")
    require(digest_object(one_shot_cells) ==
            "f261d82446d2920149c648e51e86beb60a3c25ae8ff4d419d92c0940d7b940ab",
            "one-shot matrix digest changed")
    sample = confidence((math.log(2.0), math.log(2.0), math.log(2.0)))
    require(all(abs(sample[name] - 2.0) < 1e-12
                for name in ("ratio", "lower", "upper")),
            "confidence calculation changed")
    require(ceil64(1) == 64 and ceil64(64) == 64 and ceil64(65) == 128,
            "padded main geometry changed")
    require(
        TARGET_ENDPOINT_ABSOLUTE_TOLERANCE_NS == 20_000 and
        TARGET_ENDPOINT_RELATIVE_TOLERANCE_PPM == 50 and
        target_endpoint_tolerance_ns(7_581_132) == 20_000 and
        target_endpoint_tolerance_ns(1_400_000_000) == 70_000,
        "target endpoint tolerance changed")
    require(
        SIBLING_ABSOLUTE_TOLERANCE_NS == 50_000 and
        SIBLING_RELATIVE_TOLERANCE_PPM == 200 and
        sibling_runtime_tolerance_ns(7_581_132) == 50_000 and
        sibling_runtime_tolerance_ns(1_400_000_000) == 280_000,
        "sibling runtime tolerance changed")
    require(MAX_ATTEMPTS == 5, "isolation retry policy changed")
    interval = {"ratio": 2.0, "lower": 1.5, "upper": 2.5}
    target_analysis = {
        "cell": {"id": "target", "role": "target"},
        "control_over_candidate": interval,
        "main_over_candidate": interval,
        "main_over_candidate_first_use": interval,
    }
    neighbor_analysis = {
        "cell": {"id": "neighbor", "role": "neighbor"},
        "control_over_candidate": interval,
    }
    target_only = aggregate((target_analysis,))
    neighbor_only = aggregate((neighbor_analysis,))
    require(
        target_only["target_count"] == 1 and
        target_only["neighbor_count"] == 0 and
        target_only["minimum_neighbor_upper"] is None and
        neighbor_only["target_count"] == 0 and
        neighbor_only["neighbor_count"] == 1 and
        neighbor_only["target_control_geomean"] is None,
        "single-role campaign aggregation changed")
    with tempfile.TemporaryDirectory(
            prefix="leopard2-equal-rounded-self-test-") as temporary:
        raw_root = Path(temporary) / "raw"
        raw_root.mkdir(mode=0o700)
        first, first_prior = allocate_raw_run_directory(raw_root, "cell")
        second, second_prior = allocate_raw_run_directory(raw_root, "cell")
        require(
            first_prior == [] and first != second and
            second_prior == [str(first.relative_to(raw_root.parent))],
            "resumable raw-run allocation changed")
    print("equal-rounded runner self-test passed")
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    run_parser = commands.add_parser("run")
    run_parser.add_argument("--candidate", required=True, type=Path)
    run_parser.add_argument("--control", required=True, type=Path)
    run_parser.add_argument("--main", required=True, type=Path)
    run_parser.add_argument("--source-commit", required=True)
    run_parser.add_argument("--source-tree", required=True)
    run_parser.add_argument("--output", required=True, type=Path)
    run_parser.add_argument("--cpu", required=True, type=int)
    run_parser.add_argument("--sibling", required=True, type=int)
    run_parser.add_argument("--target-rounds", type=int, default=3)
    run_parser.add_argument("--iterations", type=int, default=9)
    run_parser.add_argument("--warmup", type=int, default=2)
    run_parser.add_argument("--reuse", type=int, default=64)
    run_parser.add_argument("--cell-id", action="append", default=[])
    run_parser.add_argument("--resume", action="store_true")
    run_parser.add_argument(
        "--mode", choices=("reusable", "one-shot"), default="reusable")
    run_parser.set_defaults(function=run)
    verify_parser = commands.add_parser("verify")
    verify_parser.add_argument("--output", required=True, type=Path)
    verify_parser.set_defaults(function=verify)
    self_parser = commands.add_parser("self-test")
    self_parser.set_defaults(function=lambda unused: self_test())
    return result


def main() -> int:
    options = parser().parse_args()
    return int(options.function(options))


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except EvidenceError as error:
        print(f"evidence error: {error}", file=sys.stderr)
        raise SystemExit(1)
