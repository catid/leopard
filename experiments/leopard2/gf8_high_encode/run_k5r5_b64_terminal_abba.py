#!/usr/bin/env python3
"""Qualify the ordinary K=5/R=5/64-byte AVX2 encode terminal.

The runner consumes immutable binaries, proves that candidate and control have
identical executable sections, verifies their initialized diagnostic selector
words, pins every child to one logical CPU, reserves its SMT sibling, and
compares the target against both the same-source control and exact Leopard
main.  Neighbor cells compare only candidate and control because the terminal
must be inert there.
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import importlib.util
import json
import math
import os
import re
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA = "leopard2-gf8-k5r5-b64-terminal-abba/v1"
SUMMARY_SCHEMA = "leopard2-gf8-k5r5-b64-terminal-summary/v1"
MAIN_COMMIT = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
MODE_SYMBOL = "_ZN12_GLOBAL__N_1L24g_k5r5_b64_terminal_modeE"
ALLOW_IDENTICAL_CANDIDATE_CONTROL = False
ALLOW_MULTIPLE_TARGETS = False
CANDIDATE_SCHEMA = "leopard2-benchmark-v5"
CONTROL_SCHEMA = "leopard2-benchmark-v5"
CONTROL_EXTRA_ARGUMENTS: tuple[str, ...] = ()
CONTROL_BUILD_MARKER: str | None = None
REQUIRE_EXPECTED_IDENTITIES = False
REQUIRE_BUILD_CLOSURE = False
REQUIRE_FULL_ELF_IDENTITY = False
MAX_ISOLATION_ATTEMPTS = 3
LOCK_PATH = Path("/tmp/leopard-gf8-authoritative.lock")
T95_DF2 = 4.302652729911275
T95_DF3 = 3.182446305284263
T95_DF8 = 2.306004135204166
TARGET_CONTROL_FLOOR = 1.05
TARGET_MAIN_FLOOR = 1.0
NEIGHBOR_FLOOR = 1.0 / 1.02
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


def load_t8_support() -> Any:
    path = Path(__file__).resolve().with_name("run_t8_two_block_abba.py")
    specification = importlib.util.spec_from_file_location(
        "k5r5_b64_t8_evidence_support", path)
    require(specification is not None and specification.loader is not None,
            f"cannot load evidence support: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


T8_SUPPORT = load_t8_support()
MAIN_SUPPORT = T8_SUPPORT.SUPPORT
RUNNER_PATH = Path(__file__).resolve()
RUNNER_DEPENDENCIES = (
    RUNNER_PATH,
    Path(T8_SUPPORT.__file__).resolve(),
    Path(MAIN_SUPPORT.__file__).resolve(),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def support_file_identity(path: Path) -> dict[str, Any]:
    """Identify imported Python support without requiring executable mode."""
    resolved = path.resolve(strict=True)
    status = resolved.stat()
    require(resolved.is_file() and status.st_size > 0,
            f"runner support is not a nonempty regular file: {resolved}")
    return {
        "path": str(resolved),
        "size": status.st_size,
        "sha256": sha256(resolved),
    }


def run_tool(arguments: Sequence[str]) -> str:
    completed = subprocess.run(
        list(arguments), env=CHILD_ENVIRONMENT, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, timeout=10.0, check=False)
    require(completed.returncode == 0,
            f"{arguments[0]} failed: " +
            completed.stderr.decode("utf-8", errors="replace")[-1000:])
    return completed.stdout.decode("utf-8", errors="strict")


def mode_word_identity(executable: Path) -> dict[str, Any]:
    """Read the retained selector word by mapping its ELF symbol to the file."""
    symbols = run_tool(("/usr/bin/readelf", "-sW", str(executable)))
    matches = []
    for line in symbols.splitlines():
        tokens = line.split()
        if len(tokens) >= 8 and tokens[-1] == MODE_SYMBOL:
            matches.append(tokens)
    require(len(matches) == 1, "diagnostic selector symbol is missing or ambiguous")
    symbol = matches[0]
    require(symbol[2] == "4" and
            symbol[3:6] == ["OBJECT", "LOCAL", "DEFAULT"],
            "diagnostic selector symbol metadata changed")
    try:
        address = int(symbol[1], 16)
        section_index = int(symbol[6])
    except ValueError as error:
        raise EvidenceError("diagnostic selector symbol is not file-backed") from error

    sections = run_tool(("/usr/bin/readelf", "-SW", str(executable)))
    section_pattern = re.compile(
        r"^\s*\[\s*(\d+)\]\s+(\S+)\s+(\S+)\s+"
        r"([0-9A-Fa-f]+)\s+([0-9A-Fa-f]+)\s+([0-9A-Fa-f]+)\s+")
    selected: tuple[str, int, int, int] | None = None
    for line in sections.splitlines():
        match = section_pattern.match(line)
        if match and int(match.group(1)) == section_index:
            selected = (
                match.group(2), int(match.group(4), 16),
                int(match.group(5), 16), int(match.group(6), 16))
            break
    require(selected is not None, "diagnostic selector section is absent")
    section_name, section_address, section_offset, section_size = selected
    require(section_address <= address and
            address + 4 <= section_address + section_size,
            "diagnostic selector lies outside its section")
    file_offset = section_offset + address - section_address
    with executable.open("rb") as source:
        require(source.read(6) == b"\x7fELF\x02\x01",
                "benchmark is not little-endian ELF64")
        source.seek(file_offset)
        payload = source.read(4)
    require(len(payload) == 4, "diagnostic selector word is truncated")
    return {
        "symbol": MODE_SYMBOL,
        "virtual_address": address,
        "section_index": section_index,
        "section_name": section_name,
        "file_offset": file_offset,
        "bytes_hex": payload.hex(),
        "value": int.from_bytes(payload, "little"),
    }


def elf_sections(executable: Path) -> list[dict[str, Any]]:
    """Return the complete ELF64 section map used by full-file normalization."""
    sections = run_tool(("/usr/bin/readelf", "-SW", str(executable)))
    pattern = re.compile(
        r"^\s*\[\s*(\d+)\]\s+(\S+)\s+(\S+)\s+"
        r"([0-9A-Fa-f]+)\s+([0-9A-Fa-f]+)\s+([0-9A-Fa-f]+)\s+"
        r"\S+\s+(\S*)\s+")
    result = []
    for line in sections.splitlines():
        match = pattern.match(line)
        if match:
            result.append({
                "index": int(match.group(1)),
                "name": match.group(2),
                "type": match.group(3),
                "address": int(match.group(4), 16),
                "offset": int(match.group(5), 16),
                "size": int(match.group(6), 16),
                "flags": match.group(7),
            })
    require(any(row["name"] == ".text" for row in result) and
            any(row["name"] == ".note.gnu.build-id" for row in result),
            "ELF section map is incomplete")
    return result


def normalized_elf_identity(
    executable: Path,
    mode_word: Mapping[str, Any],
) -> dict[str, Any]:
    """Hash the full ELF after masking only build-id and selector bytes."""
    payload = bytearray(executable.read_bytes())
    sections = elf_sections(executable)
    build_id = next(
        row for row in sections if row["name"] == ".note.gnu.build-id")
    ranges = [
        (int(build_id["offset"]), int(build_id["size"]), "gnu_build_id"),
        (int(mode_word["file_offset"]), 4, "diagnostic_selector"),
    ]
    for offset, size, label in ranges:
        require(0 <= offset <= len(payload) and
                0 < size <= len(payload) - offset,
                f"normalized ELF range is outside the file: {label}")
        payload[offset:offset + size] = b"\0" * size
    section_records = []
    for row in sections:
        offset = int(row["offset"])
        size = int(row["size"])
        if row["type"] == "NOBITS" or size == 0:
            digest = None
        else:
            require(offset + size <= len(payload),
                    f"ELF section is outside the file: {row['name']}")
            digest = hashlib.sha256(payload[offset:offset + size]).hexdigest()
        section_records.append({
            "index": row["index"], "name": row["name"],
            "type": row["type"], "offset": offset, "size": size,
            "flags": row["flags"], "normalized_sha256": digest,
        })
    return {
        "size": len(payload),
        "normalized_sha256": hashlib.sha256(payload).hexdigest(),
        "normalized_ranges": [
            {"offset": offset, "size": size, "reason": label}
            for offset, size, label in ranges
        ],
        "sections": section_records,
    }


def normalized_elf_pair(
    candidate: Path,
    control: Path,
    mode_words: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    candidate_identity = normalized_elf_identity(
        candidate, mode_words["candidate"])
    control_identity = normalized_elf_identity(
        control, mode_words["control"])
    require(candidate_identity == control_identity,
            "candidate/control full ELF differs outside normalized fields")
    candidate_bytes = candidate.read_bytes()
    control_bytes = control.read_bytes()
    require(len(candidate_bytes) == len(control_bytes),
            "candidate/control ELF sizes differ")
    differing = [index for index, pair in enumerate(
        zip(candidate_bytes, control_bytes)) if pair[0] != pair[1]]
    allowed = set()
    for row in candidate_identity["normalized_ranges"]:
        allowed.update(range(row["offset"], row["offset"] + row["size"]))
    require(differing and set(differing) <= allowed,
            "candidate/control ELF has an unclassified byte difference")
    return {
        "normalized": candidate_identity,
        "differing_byte_count": len(differing),
        "differing_byte_offsets": differing,
    }


def freeze_executable(
    source: Path,
    expected_sha256: str | None,
    destination: Path,
) -> dict[str, Any]:
    source_identity = T8_SUPPORT.file_identity(source)
    if expected_sha256 is not None:
        require(source_identity["sha256"] == expected_sha256,
                f"input binary SHA-256 changed: {source}")
    descriptor = os.open(
        destination, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o555)
    try:
        with source.open("rb") as input_file, \
                os.fdopen(descriptor, "wb") as output_file:
            descriptor = -1
            for block in iter(lambda: input_file.read(1024 * 1024), b""):
                output_file.write(block)
            output_file.flush()
            os.fsync(output_file.fileno())
        destination.chmod(0o555)
    except BaseException:
        if descriptor >= 0:
            os.close(descriptor)
        try:
            destination.unlink()
        except OSError:
            pass
        raise
    frozen_identity = T8_SUPPORT.file_identity(destination)
    require(frozen_identity["sha256"] == source_identity["sha256"] and
            frozen_identity["mode"] == 0o555,
            f"frozen binary identity is invalid: {destination}")
    return {"input": source_identity, "frozen": frozen_identity}


def normalized_compile_commands_identity(
    candidate: Path,
    control: Path,
) -> dict[str, Any]:
    def normalize(path: Path, diagnostic_value: str) -> tuple[list[Any], str]:
        value = json.loads(path.read_text(encoding="utf-8"))
        require(isinstance(value, list) and value,
                "compile_commands is not a nonempty array")
        directories = {
            row.get("directory") for row in value if isinstance(row, dict)
        }
        require(len(directories) == 1 and
                all(isinstance(row, dict) for row in value),
                "compile_commands has ambiguous build roots")
        build_root = next(iter(directories))
        require(isinstance(build_root, str) and build_root,
                "compile_commands build root is invalid")
        normalized = []
        t32_rows = []
        for row in value:
            current = dict(row)
            for key, item in tuple(current.items()):
                if isinstance(item, str):
                    item = item.replace(build_root, "${BUILD}")
                    item = re.sub(
                        r"LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED="
                        + re.escape(diagnostic_value) + r"\b",
                        "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED=${MODE}",
                        item)
                    expected_route = "1" if diagnostic_value == "0" else "0"
                    item = re.sub(
                        r"LEO2_EXPECT_T32_B256_GENERATED=" + expected_route +
                        r"\b",
                        "LEO2_EXPECT_T32_B256_GENERATED=${EXPECT}", item)
                    item = re.sub(
                        r"LEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256=\\\""
                        r"[0-9a-f]{64}\\\"",
                        "LEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256=\\\"${CONFIG}\\\"",
                        item)
                    current[key] = item
            normalized.append(current)
            if Path(str(row.get("file", ""))).name == \
                    "Leopard2BackendAVX2T32B256.cpp":
                t32_rows.append(str(row.get("command", "")))
        require(len(t32_rows) == 1 and
                f"LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_GENERATED={diagnostic_value}" in
                    t32_rows[0] and
                "-mavx2" in t32_rows[0] and "-mno-avx512f" in t32_rows[0],
                "T32 compile-command contract changed")
        canonical = json.dumps(
            normalized, sort_keys=True, separators=(",", ":"))
        return normalized, hashlib.sha256(canonical.encode()).hexdigest()

    candidate_value, candidate_hash = normalize(candidate, "0")
    control_value, control_hash = normalize(control, "1")
    require(candidate_value == control_value and candidate_hash == control_hash,
            "candidate/control normalized compile commands differ")
    return {
        "normalized_sha256": candidate_hash,
        "entry_count": len(candidate_value),
    }


def build_closure_identity(options: argparse.Namespace) -> dict[str, Any]:
    values = {
        "candidate_archive": (
            options.candidate_archive, options.candidate_archive_sha256),
        "control_archive": (
            options.control_archive, options.control_archive_sha256),
        "candidate_compile_commands": (
            options.candidate_compile_commands,
            options.candidate_compile_commands_sha256),
        "control_compile_commands": (
            options.control_compile_commands,
            options.control_compile_commands_sha256),
    }
    identities = {}
    for name, (path, expected) in values.items():
        require(path is not None and expected is not None,
                f"build-closure input is required: {name}")
        identity = T8_SUPPORT.regular_file_identity(path)
        require(identity["sha256"] == expected,
                f"build-closure SHA-256 changed: {name}")
        identities[name] = identity
    candidate_archive = Path(str(identities["candidate_archive"]["path"]))
    control_archive = Path(str(identities["control_archive"]["path"]))
    candidate_payload = candidate_archive.read_bytes()
    control_payload = control_archive.read_bytes()
    require(len(candidate_payload) == len(control_payload),
            "candidate/control archive sizes differ")
    differences = [index for index, pair in enumerate(
        zip(candidate_payload, control_payload)) if pair[0] != pair[1]]
    require(len(differences) == 1 and
            candidate_payload[differences[0]] == 1 and
            control_payload[differences[0]] == 2,
            "candidate/control archives differ beyond the selector byte")
    candidate_members = run_tool(("/usr/bin/ar", "t", str(candidate_archive))).splitlines()
    control_members = run_tool(("/usr/bin/ar", "t", str(control_archive))).splitlines()
    require(candidate_members == control_members and
            candidate_members.count("Leopard2BackendAVX2T32B256.cpp.o") == 1,
            "candidate/control archive member closure changed")
    compile_identity = normalized_compile_commands_identity(
        Path(str(identities["candidate_compile_commands"]["path"])),
        Path(str(identities["control_compile_commands"]["path"])))
    return {
        "inputs": identities,
        "archive_members": candidate_members,
        "archive_selector_byte_offset": differences[0],
        "normalized_compile_commands": compile_identity,
    }


def cells() -> list[dict[str, Any]]:
    values = [
        ("target-k5-r5-b64-q1", 5, 5, 64, 1, "target"),
        ("byte-neighbor-k5-r5-b63-q1", 5, 5, 63, 1, "neighbor"),
        ("byte-neighbor-k5-r5-b65-q1", 5, 5, 65, 1, "neighbor"),
        ("byte-neighbor-k5-r5-b128-q1", 5, 5, 128, 1, "neighbor"),
        ("shape-neighbor-k4-r4-b64-q1", 4, 4, 64, 1, "neighbor"),
        ("shape-neighbor-k5-r4-b64-q1", 5, 4, 64, 1, "neighbor"),
        ("shape-neighbor-k6-r5-b64-q1", 6, 5, 64, 1, "neighbor"),
        ("shape-neighbor-k5-r6-b64-q1", 5, 6, 64, 1, "neighbor"),
        ("terminal-neighbor-k1-r1-b64-q1", 1, 1, 64, 1, "neighbor"),
        ("terminal-neighbor-k2-r1-b64-q1", 2, 1, 64, 1, "neighbor"),
        ("batch-neighbor-k5-r5-b64-q2", 5, 5, 64, 2, "neighbor"),
        ("batch-neighbor-k5-r5-b64-q8", 5, 5, 64, 8, "neighbor"),
    ]
    result = []
    for index, (name, k, r, shard_bytes, batch, role) in enumerate(values):
        result.append({
            "id": name,
            "K": k,
            "R": r,
            "bytes": shard_bytes,
            "batch": batch,
            "reuse": max(8192, 65536 // batch),
            "role": role,
            "seed": 0x5A5B6400 + index,
        })
    return result


def benchmark_command(
    implementation: str,
    executable: Path,
    cell: Mapping[str, Any],
    cpu: int,
    iterations: int,
    warmup: int,
) -> list[str]:
    logical_bytes = int(cell["bytes"])
    # The legacy API requires a 64-byte physical shard quantum.  Its
    # source-attested benchmark has an explicit zero-padded adapter so exact
    # main can still be compared with Leopard2 arbitrary-byte neighbors while
    # hashing only the common logical prefix.
    physical_bytes = logical_bytes if implementation != "main" else \
        (logical_bytes + 63) & ~63
    common = [
        "/usr/bin/prlimit", "--as=201326592",
        "/usr/bin/taskset", "-c", str(cpu), str(executable),
        "--k", str(cell["K"]), "--r", str(cell["R"]),
        "--bytes", str(physical_bytes),
        "--loss", str(cell.get("loss", 1)),
        "--batch", str(cell["batch"]), "--reuse", str(cell["reuse"]),
        "--iterations", str(iterations), "--warmup", str(warmup),
        "--threads", "1", "--seed", str(cell["seed"]), "--json", "-",
    ]
    if implementation == "main":
        if physical_bytes != logical_bytes:
            common[-2:-2] = ["--logical-bytes", str(logical_bytes)]
        return common
    command = common[:-2] + [
        "--profile", "high", "--field", "gf8", "--backend", "avx2",
        "--skip-legacy", "--retain-samples", "--attest-source",
    ]
    if cell.get("measure_one_shot") is True:
        command.append("--measure-one-shot-encode")
    if implementation == "control":
        command.extend(CONTROL_EXTRA_ARGUMENTS)
    return command + ["--json", "-"]


def positive_encode_metric(
    result: Mapping[str, Any],
    metric_name: str,
    iterations: int,
) -> dict[str, Any]:
    metrics = result.get("metrics")
    require(isinstance(metrics, dict), "benchmark metrics are absent")
    encode = metrics.get(metric_name)
    require(isinstance(encode, dict), f"{metric_name} metric is absent")
    median = encode.get("median_us_per_batch_call")
    mad = encode.get("mad_us_per_batch_call")
    minimum = encode.get("minimum_us_per_batch_call")
    maximum = encode.get("maximum_us_per_batch_call")
    samples = encode.get("samples_us_per_batch_call")
    require(isinstance(median, (int, float)) and
            not isinstance(median, bool) and
            math.isfinite(float(median)) and float(median) > 0,
            f"{metric_name} median is not finite and positive")
    require(isinstance(samples, list) and len(samples) == iterations and
            all(isinstance(value, (int, float)) and
                not isinstance(value, bool) and
                math.isfinite(float(value)) and float(value) > 0
                for value in samples),
            f"{metric_name} raw samples are incomplete")
    numeric = [float(value) for value in samples]
    recomputed_median = statistics.median(numeric)
    recomputed_mad = statistics.median(
        abs(value - recomputed_median) for value in numeric)
    reported = (mad, minimum, maximum)
    require(all(isinstance(value, (int, float)) and
                not isinstance(value, bool) and math.isfinite(float(value))
                for value in reported) and
            math.isclose(float(median), recomputed_median,
                         rel_tol=0.0, abs_tol=2e-6) and
            math.isclose(float(mad), recomputed_mad,
                         rel_tol=0.0, abs_tol=2e-6) and
            math.isclose(float(minimum), min(numeric),
                         rel_tol=0.0, abs_tol=2e-6) and
            math.isclose(float(maximum), max(numeric),
                         rel_tol=0.0, abs_tol=2e-6),
            f"{metric_name} summary disagrees with retained raw samples")
    return {"median_us": recomputed_median, "samples_us": numeric}


def validate_result(
    implementation: str,
    result: object,
    cell: Mapping[str, Any],
    source_commit: str,
    source_tree: str,
    iterations: int,
    warmup: int,
) -> dict[str, Any]:
    require(isinstance(result, dict), "benchmark output is not an object")
    logical_bytes = int(cell["bytes"])
    physical_bytes = logical_bytes if implementation != "main" else \
        (logical_bytes + 63) & ~63
    expected_parameters = {
        "K": cell["K"], "R": cell["R"],
        "shard_bytes": physical_bytes,
        "loss_count": cell.get("loss", 1),
        "batch": cell["batch"], "reuse": cell["reuse"],
        "iterations": iterations, "warmup": warmup,
        "thread_count": 1, "seed": cell["seed"],
    }
    if implementation == "main":
        expected_parameters["logical_shard_bytes"] = logical_bytes
    elif cell.get("measure_one_shot") is True:
        expected_parameters["measure_one_shot_encode"] = True
    parameters = result.get("parameters")
    require(isinstance(parameters, dict) and
            all(parameters.get(name) == value
                for name, value in expected_parameters.items()),
            "benchmark parameters differ from the frozen cell")
    resolved = result.get("resolved")
    correctness = result.get("correctness")
    digests = result.get("workload_digests")
    require(isinstance(resolved, dict) and
            resolved.get("profile") == "legacy_high_v1" and
            resolved.get("field") == "gf8" and
            resolved.get("thread_count") == 1,
            "benchmark resolved a different profile, field, or thread count")
    require(isinstance(digests, dict) and
            digests.get("algorithm") == "fnv1a64" and
            all(isinstance(digests.get(name), str) and len(digests[name]) == 16
                for name in (
                    "original_data", "transmitted_parity",
                    "recovered_originals")),
            "benchmark workload digests are incomplete")
    if implementation == "main":
        padded_application = physical_bytes != logical_bytes
        require(result.get("schema") == (
                    "leopard-main-benchmark-v2" if padded_application else
                    "leopard-main-benchmark-v1") and
                isinstance(correctness, dict) and
                correctness.get("round_trip") is True and
                correctness.get("logical_prefix_fingerprinted") is True and
                resolved.get("padded_application_bytes") is
                    padded_application and
                resolved.get("padding_policy") == "zero suffix per shard",
                "exact-main identity or round trip failed")
        build = result.get("build")
        require(isinstance(build, dict) and
                build.get("main_source_commit") == MAIN_COMMIT,
                "exact-main source identity changed")
    else:
        expected_schema = CONTROL_SCHEMA if implementation == "control" \
            else CANDIDATE_SCHEMA
        require(result.get("schema") == expected_schema and
                resolved.get("backend") == "avx2" and
                isinstance(correctness, dict) and
                correctness.get("leopard2_round_trip") is True,
                "Leopard2 identity, backend, or round trip failed")
        build = result.get("build")
        require(isinstance(build, dict) and
                build.get("source_commit") == source_commit and
                build.get("source_tree") == source_tree and
                build.get("source_tracked_dirty") is False,
                "Leopard2 embedded source identity changed")
        if CONTROL_BUILD_MARKER is not None:
            if implementation == "control":
                require(build.get(CONTROL_BUILD_MARKER) is True,
                        "runtime attribution marker differs from the label")
            else:
                require(CONTROL_BUILD_MARKER not in build or
                        build[CONTROL_BUILD_MARKER] is False,
                        "candidate unexpectedly reports a disabled terminal")
    encode_metric = positive_encode_metric(
        result, "encode_execution", iterations)
    one_shot_metric = None
    if cell.get("measure_one_shot") is True:
        one_shot_metric = encode_metric if implementation == "main" else \
            positive_encode_metric(result, "one_shot_encode", iterations)
    return {
        "encode_us": encode_metric["median_us"],
        "encode_samples_us": encode_metric["samples_us"],
        "one_shot_encode_us": None if one_shot_metric is None else
            one_shot_metric["median_us"],
        "one_shot_encode_samples_us": None if one_shot_metric is None else
            one_shot_metric["samples_us"],
        "digests": dict(digests),
        "schema": result["schema"],
    }


def retain_failure_streams(
    failure_output: Path,
    implementation: str,
    started_ns: int,
    stdout: bytes | str | None,
    stderr: bytes | str | None,
) -> dict[str, str]:
    """Retain child output for every rejected invocation, including timeouts."""
    def as_bytes(value: bytes | str | None) -> bytes:
        if value is None:
            return b""
        if isinstance(value, bytes):
            return value
        return value.encode("utf-8", errors="replace")

    prefix = failure_output / f"failure-{implementation}-{started_ns}"
    stdout_path = prefix.with_suffix(".stdout")
    stderr_path = prefix.with_suffix(".stderr")
    T8_SUPPORT.write_bytes_exclusive(stdout_path, as_bytes(stdout))
    T8_SUPPORT.write_bytes_exclusive(stderr_path, as_bytes(stderr))
    return {"stdout": str(stdout_path), "stderr": str(stderr_path)}


def run_one(
    implementation: str,
    identity: Mapping[str, Any],
    cell: Mapping[str, Any],
    cpu: int,
    source_commit: str,
    source_tree: str,
    iterations: int,
    warmup: int,
    failure_output: Path,
) -> dict[str, Any]:
    executable = Path(str(identity["path"]))
    require(sha256(executable) == identity["sha256"],
            f"{implementation} binary changed before execution")
    command = benchmark_command(
        implementation, executable, cell, cpu, iterations, warmup)
    started_ns = time.monotonic_ns()
    try:
        completed = subprocess.run(
            command, env=CHILD_ENVIRONMENT, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, timeout=60.0, check=False)
    except subprocess.TimeoutExpired as error:
        retained = retain_failure_streams(
            failure_output, implementation, started_ns,
            error.stdout, error.stderr)
        raise EvidenceError(
            f"{implementation} timed out; retained output: {retained}") \
            from error
    elapsed_ns = time.monotonic_ns() - started_ns
    if completed.returncode != 0:
        retained = retain_failure_streams(
            failure_output, implementation, started_ns,
            completed.stdout, completed.stderr)
        raise EvidenceError(
            f"{implementation} failed: " +
            completed.stderr.decode("utf-8", errors="replace")[-1000:] +
            f"; retained output: {retained}")
    try:
        require(sha256(executable) == identity["sha256"],
                f"{implementation} binary changed after execution")
        result = json.loads(completed.stdout.decode("utf-8"))
        normalized = validate_result(
            implementation, result, cell, source_commit, source_tree,
            iterations, warmup)
    except Exception as error:
        retained = retain_failure_streams(
            failure_output, implementation, started_ns,
            completed.stdout, completed.stderr)
        raise EvidenceError(
            f"{implementation} output was rejected: {error}; "
            f"retained output: {retained}") \
            from error
    return {
        "implementation": implementation,
        "command": command,
        "elapsed_ns": elapsed_ns,
        "stdout_sha256": hashlib.sha256(completed.stdout).hexdigest(),
        "stderr_sha256": hashlib.sha256(completed.stderr).hexdigest(),
        "normalized": normalized,
        "result": result,
    }


def confidence_interval(values: Sequence[float]) -> dict[str, Any]:
    require(len(values) in (3, 4, 9),
            "three, four, or nine independent round contrasts are required")
    center = statistics.mean(values)
    critical = {
        3: T95_DF2,
        4: T95_DF3,
        9: T95_DF8,
    }[len(values)]
    half_width = critical * statistics.stdev(values) / math.sqrt(len(values))
    return {
        "speedup": math.exp(center),
        "ci95": [
            math.exp(center - half_width),
            math.exp(center + half_width),
        ],
        "round_log_ratios": list(values),
    }


def analyze(
    cell: Mapping[str, Any],
    rounds: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    reference = rounds[0]["invocations"][0]["normalized"]["digests"]
    labels = ("control", "main") if cell.get(
        "compare_main", cell["role"] == "target") \
        else ("control",)
    metric_names = ["encode_us"]
    if cell.get("measure_one_shot") is True:
        metric_names.append("one_shot_encode_us")
    contrasts: dict[str, dict[str, list[float]]] = {
        metric: {label: [] for label in labels} for metric in metric_names
    }
    for round_value in rounds:
        require(round_value["isolation"]["accepted"] is True,
                "contaminated round cannot be analyzed")
        invocations = round_value["invocations"]
        require(all(item["normalized"]["digests"] == reference
                    for item in invocations),
                "candidate/control/main workload digests differ")
        for metric in metric_names:
            candidate = [
                item["normalized"][metric] for item in invocations
                if item["implementation"] == "candidate"]
            require(len(candidate) == 2 and
                    all(isinstance(value, (int, float)) and value > 0
                        for value in candidate),
                    f"round lacks two positive candidate {metric} observations")
            candidate_log = statistics.mean(
                math.log(value) for value in candidate)
            for label in labels:
                baseline = [
                    item["normalized"][metric] for item in invocations
                    if item["implementation"] == label]
                require(len(baseline) == 2 and
                        all(isinstance(value, (int, float)) and value > 0
                            for value in baseline),
                        f"round lacks two positive {label} {metric} observations")
                contrasts[metric][label].append(
                    statistics.mean(math.log(value) for value in baseline) -
                    candidate_log)
    output = {
        "cell": dict(cell),
        "digests": reference,
        "control_over_candidate": confidence_interval(
            contrasts["encode_us"]["control"]),
    }
    if "main" in labels:
        output["main_over_candidate"] = confidence_interval(
            contrasts["encode_us"]["main"])
    if "one_shot_encode_us" in contrasts:
        output["one_shot_control_over_candidate"] = confidence_interval(
            contrasts["one_shot_encode_us"]["control"])
        if "main" in labels:
            output["one_shot_main_over_candidate"] = confidence_interval(
                contrasts["one_shot_encode_us"]["main"])
    return output


def acquire_global_lock() -> int:
    descriptor = os.open(LOCK_PATH, os.O_RDWR | os.O_CREAT, 0o600)
    try:
        fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        os.close(descriptor)
        raise EvidenceError(f"authoritative lock is busy: {LOCK_PATH}")
    return descriptor


def select_round_orders(
    orders: Sequence[Sequence[str]],
    requested_rounds: int | None,
) -> tuple[tuple[str, ...], ...]:
    """Repeat a balanced order cycle to a supported independent-round count."""
    require(len(orders) > 0, "round-order cycle is empty")
    round_count = len(orders) if requested_rounds is None else requested_rounds
    require(round_count in (3, 4, 9),
            "--rounds must select 3, 4, or 9 independent contrasts")
    return tuple(tuple(orders[index % len(orders)])
                 for index in range(round_count))


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--candidate-sha256")
    parser.add_argument("--control", required=True, type=Path)
    parser.add_argument("--control-sha256")
    parser.add_argument("--main", required=True, type=Path)
    parser.add_argument("--main-sha256")
    parser.add_argument("--candidate-archive", type=Path)
    parser.add_argument("--candidate-archive-sha256")
    parser.add_argument("--control-archive", type=Path)
    parser.add_argument("--control-archive-sha256")
    parser.add_argument("--candidate-compile-commands", type=Path)
    parser.add_argument("--candidate-compile-commands-sha256")
    parser.add_argument("--control-compile-commands", type=Path)
    parser.add_argument("--control-compile-commands-sha256")
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--source-tree", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cpu", required=True, type=int)
    parser.add_argument("--sibling", required=True, type=int)
    parser.add_argument("--iterations", type=int, default=15)
    parser.add_argument("--warmup", type=int, default=64)
    parser.add_argument("--rounds", type=int)
    return parser.parse_args()


def main() -> int:
    options = parse_arguments()
    require(options.iterations >= 3 and options.warmup >= 1,
            "insufficient benchmark repetitions")
    require(not options.output.exists(), "output path already exists")
    expected_hashes = {
        "candidate": options.candidate_sha256,
        "control": options.control_sha256,
        "main": options.main_sha256,
    }
    if REQUIRE_EXPECTED_IDENTITIES:
        require(all(value is not None for value in expected_hashes.values()),
                "caller-supplied candidate/control/main SHA-256 is required")
    require(all(value is None or
                re.fullmatch(r"[0-9a-f]{64}", value) is not None
                for value in expected_hashes.values()),
            "an expected executable SHA-256 is malformed")
    raw: dict[str, Any] = {
        "schema": SCHEMA,
        "created_utc": MAIN_SUPPORT.utc_now(),
    }
    lock_descriptor: int | None = None
    try:
        options.output.mkdir(parents=True)
        raw.update({
            "source_commit": options.source_commit,
            "source_tree": options.source_tree,
            "main_commit": MAIN_COMMIT,
            "cpu": options.cpu,
            "reserved_sibling": options.sibling,
            "iterations": options.iterations,
            "warmup": options.warmup,
            "requested_rounds": options.rounds,
            "runner": T8_SUPPORT.file_identity(RUNNER_PATH),
            "runner_dependencies": [
                support_file_identity(path)
                for path in RUNNER_DEPENDENCIES
            ],
            "cells": [],
        })
        lock_descriptor = acquire_global_lock()
        frozen_directory = options.output / "frozen"
        frozen_directory.mkdir()
        source_paths = {
            "candidate": options.candidate,
            "control": options.control,
            "main": options.main,
        }
        freeze_records = {
            name: freeze_executable(
                source_paths[name], expected_hashes[name],
                frozen_directory / name)
            for name in ("candidate", "control", "main")
        }
        raw["frozen_executables"] = freeze_records
        input_identities = {
            name: record["input"] for name, record in freeze_records.items()}
        identities = {
            name: record["frozen"] for name, record in freeze_records.items()}
        raw["input_identities"] = input_identities
        raw["identities"] = identities
        closure_requested = REQUIRE_BUILD_CLOSURE or any(
            value is not None for value in (
                options.candidate_archive,
                options.candidate_archive_sha256,
                options.control_archive,
                options.control_archive_sha256,
                options.candidate_compile_commands,
                options.candidate_compile_commands_sha256,
                options.control_compile_commands,
                options.control_compile_commands_sha256))
        build_closure = build_closure_identity(options) \
            if closure_requested else None
        raw["build_closure"] = build_closure
        if ALLOW_IDENTICAL_CANDIDATE_CONTROL:
            require(
                identities["candidate"]["sha256"] ==
                    identities["control"]["sha256"] and
                identities["candidate"]["sha256"] !=
                    identities["main"]["sha256"],
                "runtime attribution requires one shared candidate/control binary")
        else:
            require(len({item["sha256"] for item in identities.values()}) == 3,
                    "candidate, control, and main binaries are not distinct")
        executable_sections = {
            name: T8_SUPPORT.executable_sections_identity(
                Path(str(identities[name]["path"])))
            for name in ("candidate", "control")
        }
        raw["executable_sections"] = executable_sections
        require(executable_sections["candidate"]["sections"] ==
                executable_sections["control"]["sections"],
                "candidate and control executable instruction sections differ")
        if ALLOW_IDENTICAL_CANDIDATE_CONTROL:
            shared_mode = mode_word_identity(
                Path(str(identities["candidate"]["path"])))
            mode_words = {"shared_binary_default": shared_mode}
        else:
            mode_words = {
                name: mode_word_identity(Path(str(identities[name]["path"])))
                for name in ("candidate", "control")
            }
        raw["mode_words"] = mode_words
        if ALLOW_IDENTICAL_CANDIDATE_CONTROL:
            require(mode_words["shared_binary_default"]["value"] == 1,
                    "shared binary does not default to the candidate route")
        else:
            require(mode_words["candidate"]["value"] == 1 and
                    mode_words["control"]["value"] == 2,
                    "candidate/control selectors are not enabled/disabled")
            if REQUIRE_FULL_ELF_IDENTITY:
                raw["normalized_candidate_control_elf"] = normalized_elf_pair(
                    Path(str(identities["candidate"]["path"])),
                    Path(str(identities["control"]["path"])), mode_words)
        require(set(os.sched_getaffinity(0)) == {options.cpu},
                "runner must be singleton-pinned to the benchmark CPU")
        sibling_text = Path(
            f"/sys/devices/system/cpu/cpu{options.cpu}/topology/"
            "thread_siblings_list").read_text(encoding="ascii")
        require(MAIN_SUPPORT.parse_cpu_list(sibling_text) ==
                {options.cpu, options.sibling},
                "requested CPUs are not one SMT pair")
        raw["host"] = {
            "runner_affinity": sorted(os.sched_getaffinity(0)),
            "benchmark_cpu": MAIN_SUPPORT.cpu_policy_identity(options.cpu),
            "reserved_sibling": MAIN_SUPPORT.cpu_policy_identity(options.sibling),
        }
        with MAIN_SUPPORT.StableLeaseAnchor(), \
                MAIN_SUPPORT.PairLease(
                    options.cpu, options.sibling) as pair_lease:
            raw["pair_lease"] = pair_lease
            before_cpu = MAIN_SUPPORT.cpu_stat_snapshot(options.cpu)
            before_sibling = MAIN_SUPPORT.cpu_stat_snapshot(options.sibling)
            before_ns = time.monotonic_ns()
            time.sleep(5.0)
            presample = MAIN_SUPPORT.isolation_record(
                options.cpu, options.sibling, pair_lease,
                before_ns, time.monotonic_ns(), before_cpu,
                MAIN_SUPPORT.cpu_stat_snapshot(options.cpu), before_sibling,
                MAIN_SUPPORT.cpu_stat_snapshot(options.sibling))
            raw["presample"] = presample
            require(presample["delta"]["benchmark_cpu"]["nonidle_jiffies"] <= 1 and
                    presample["delta"]["reserved_sibling"]["nonidle_jiffies"] == 0,
                    "CPU pair was not quiet during the presample")

            all_cells = cells()
            for cell_index, cell in enumerate(all_cells):
                cell_raw = {"cell": dict(cell), "rounds": []}
                raw["cells"].append(cell_raw)
                order_cycle = TARGET_ORDER if cell.get(
                    "compare_main", cell["role"] == "target") \
                    else NEIGHBOR_ORDER
                orders = select_round_orders(order_cycle, options.rounds)
                for round_index, order in enumerate(orders):
                    discarded_attempts = []
                    for attempt_index in range(MAX_ISOLATION_ATTEMPTS):
                        before_cpu = MAIN_SUPPORT.cpu_stat_snapshot(options.cpu)
                        before_sibling = MAIN_SUPPORT.cpu_stat_snapshot(
                            options.sibling)
                        before_ns = time.monotonic_ns()
                        invocations = [
                            run_one(
                                label, identities[label], cell, options.cpu,
                                options.source_commit, options.source_tree,
                                options.iterations, options.warmup,
                                options.output)
                            for label in order
                        ]
                        isolation = MAIN_SUPPORT.isolation_record(
                            options.cpu, options.sibling, pair_lease,
                            before_ns, time.monotonic_ns(), before_cpu,
                            MAIN_SUPPORT.cpu_stat_snapshot(options.cpu),
                            before_sibling,
                            MAIN_SUPPORT.cpu_stat_snapshot(options.sibling))
                        attempt = {
                            "attempt": attempt_index,
                            "order": list(order),
                            "invocations": invocations,
                            "isolation": isolation,
                        }
                        if isolation["accepted"]:
                            attempt["round"] = round_index
                            attempt["discarded_attempts"] = discarded_attempts
                            cell_raw["rounds"].append(attempt)
                            break
                        discarded_attempts.append(attempt)
                    else:
                        cell_raw["failed_round"] = {
                            "round": round_index,
                            "discarded_attempts": discarded_attempts,
                        }
                        require(False,
                                f"contaminated {cell['id']} round "
                                f"{round_index} for "
                                f"{MAX_ISOLATION_ATTEMPTS} attempts")
                print(f"{cell_index + 1}/{len(all_cells)} {cell['id']}",
                      file=sys.stderr, flush=True)

        post_identities = {
            name: T8_SUPPORT.file_identity(Path(str(identity["path"])))
            for name, identity in identities.items()
        }
        require(post_identities == identities,
                "a benchmark binary identity changed during the campaign")
        raw["identities_after"] = post_identities
        post_input_identities = {
            name: T8_SUPPORT.file_identity(source_paths[name])
            for name in source_paths
        }
        require(post_input_identities == input_identities,
                "an input executable identity changed during the campaign")
        raw["input_identities_after"] = post_input_identities
        post_runner = T8_SUPPORT.file_identity(RUNNER_PATH)
        post_runner_dependencies = [
            support_file_identity(path) for path in RUNNER_DEPENDENCIES]
        require(post_runner == raw["runner"] and
                post_runner_dependencies == raw["runner_dependencies"],
                "runner or imported support changed during the campaign")
        raw["runner_after"] = post_runner
        raw["runner_dependencies_after"] = post_runner_dependencies
        post_build_closure = build_closure_identity(options) \
            if closure_requested else None
        require(post_build_closure == build_closure,
                "build closure changed during the campaign")
        raw["build_closure_after"] = post_build_closure
        analyses = [
            analyze(item["cell"], item["rounds"]) for item in raw["cells"]]
        targets = [item for item in analyses
                   if item["cell"]["role"] == "target"]
        require(len(targets) >= 1, "campaign has no target cell")
        if not ALLOW_MULTIPLE_TARGETS:
            require(len(targets) == 1,
                    "campaign unexpectedly has multiple target cells")
        target_batch_control_failure = any(
            target["control_over_candidate"]["ci95"][0] <
                TARGET_CONTROL_FLOOR
            for target in targets)
        target_batch_main_failure = any(
            target["main_over_candidate"]["ci95"][0] < TARGET_MAIN_FLOOR
            for target in targets)
        one_shot_targets = [
            target for target in targets
            if target["cell"].get("measure_one_shot") is True]
        target_one_shot_control_failure = any(
            target["one_shot_control_over_candidate"]["ci95"][0] <
                TARGET_CONTROL_FLOOR
            for target in one_shot_targets)
        target_one_shot_main_failure = any(
            target["one_shot_main_over_candidate"]["ci95"][0] <
                TARGET_MAIN_FLOOR
            for target in one_shot_targets)
        target_control_failure = (
            target_batch_control_failure or target_one_shot_control_failure)
        target_main_failure = (
            target_batch_main_failure or target_one_shot_main_failure)
        credible_neighbor_regressions = [
            item["cell"]["id"] for item in analyses
            if item["cell"]["role"] == "neighbor" and
            item["control_over_candidate"]["ci95"][1] < NEIGHBOR_FLOOR
        ]
        credible_neighbor_one_shot_regressions = [
            item["cell"]["id"] for item in analyses
            if item["cell"]["role"] == "neighbor" and
            item["cell"].get("measure_one_shot") is True and
            item["one_shot_control_over_candidate"]["ci95"][1] <
                NEIGHBOR_FLOOR
        ]
        raw["completed_utc"] = MAIN_SUPPORT.utc_now()
        T8_SUPPORT.write_exclusive(options.output / "raw.json", raw)
        summary = {
            "schema": SUMMARY_SCHEMA,
            "status": "accepted" if not (
                target_control_failure or target_main_failure or
                credible_neighbor_regressions or
                credible_neighbor_one_shot_regressions) else "rejected",
            "source_commit": options.source_commit,
            "source_tree": options.source_tree,
            "main_commit": MAIN_COMMIT,
            "cell_count": len(analyses),
            "process_count": sum(
                len(round_value["invocations"])
                for cell_value in raw["cells"]
                for round_value in cell_value["rounds"]),
            "discarded_process_count": sum(
                len(attempt["invocations"])
                for cell_value in raw["cells"]
                for round_value in cell_value["rounds"]
                for attempt in round_value["discarded_attempts"]),
            "discarded_round_attempts": sum(
                len(round_value["discarded_attempts"])
                for cell_value in raw["cells"]
                for round_value in cell_value["rounds"]),
            "all_digests_matched": True,
            "all_rounds_zero_sibling_nonidle": all(
                round_value["isolation"]["delta"]
                    ["reserved_sibling"]["nonidle_jiffies"] == 0
                for cell_value in raw["cells"]
                for round_value in cell_value["rounds"]),
            "target_control_failure": target_control_failure,
            "target_main_failure": target_main_failure,
            "target_control_failure_by_metric": {
                "encode_execution": target_batch_control_failure,
                "one_shot_encode": target_one_shot_control_failure,
            },
            "target_main_failure_by_metric": {
                "encode_execution": target_batch_main_failure,
                "one_shot_encode": target_one_shot_main_failure,
            },
            "credible_neighbor_regressions": credible_neighbor_regressions,
            "credible_neighbor_one_shot_regressions":
                credible_neighbor_one_shot_regressions,
            "cells": analyses,
            "binary_sha256": {
                name: identity["sha256"] for name, identity in identities.items()},
            "candidate_control_executable_sections_sha256":
                executable_sections["candidate"]["combined_sha256"],
            "mode_words": mode_words,
            "raw_sha256": sha256(options.output / "raw.json"),
        }
        T8_SUPPORT.write_exclusive(options.output / "summary.json", summary)
        result_line = {
            "status": summary["status"],
            "cells": summary["cell_count"],
            "processes": summary["process_count"],
            "discarded_processes": summary["discarded_process_count"],
            "discarded_round_attempts": summary["discarded_round_attempts"],
            "credible_neighbor_regressions": credible_neighbor_regressions,
            "credible_neighbor_one_shot_regressions":
                credible_neighbor_one_shot_regressions,
        }
        if len(targets) == 1:
            result_line["target_control"] = \
                targets[0]["control_over_candidate"]
            result_line["target_main"] = targets[0]["main_over_candidate"]
            if targets[0]["cell"].get("measure_one_shot") is True:
                result_line["target_one_shot_control"] = \
                    targets[0]["one_shot_control_over_candidate"]
                result_line["target_one_shot_main"] = \
                    targets[0]["one_shot_main_over_candidate"]
        else:
            result_line["targets"] = {
                target["cell"]["id"]: {
                    "control": target["control_over_candidate"],
                    "main": target["main_over_candidate"],
                    "one_shot_control": target.get(
                        "one_shot_control_over_candidate"),
                    "one_shot_main": target.get(
                        "one_shot_main_over_candidate"),
                }
                for target in targets
            }
        print(json.dumps(result_line, sort_keys=True))
        return 0 if summary["status"] == "accepted" else 2
    except Exception as error:
        raw["failed_utc"] = MAIN_SUPPORT.utc_now()
        raw["failure"] = {"type": type(error).__name__, "message": str(error)}
        T8_SUPPORT.write_exclusive(options.output / "failure.json", raw)
        print(f"evidence rejected: {error}", file=sys.stderr)
        return 1
    finally:
        if lock_descriptor is not None:
            os.close(lock_descriptor)


if __name__ == "__main__":
    raise SystemExit(main())
