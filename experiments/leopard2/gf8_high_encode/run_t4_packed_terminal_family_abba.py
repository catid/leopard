#!/usr/bin/env python3
"""Measure the ordinary one-item GF8 T=4 packed-terminal family.

The runner consumes three immutable, already-built executables: a clean
candidate, a same-source compile-time control whose
LEO2_DIAGNOSTIC_DISABLE_HIGH_T4_BATCH_BINDING value is one, and the exact
Leopard main benchmark.  It verifies embedded source identity, the retained
candidate/control selector words, executable-section identity, workload
digests, singleton CPU pinning, and an exclusive SMT-pair lease before
reporting balanced ABBA confidence intervals.

This campaign is descriptive.  Integrity or isolation failures reject the
evidence, but a slower performance cell is recorded rather than changing the
process exit status.
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
from typing import Any, Callable, Mapping, Sequence


SCHEMA = "leopard2-gf8-t4-packed-terminal-family-abba/v1"
SUMMARY_SCHEMA = "leopard2-gf8-t4-packed-terminal-family-summary/v1"
MAIN_COMMIT = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
MODE_SYMBOL = "_ZN12_GLOBAL__N_1L28g_high_t4_batch_binding_modeE"
LOCK_PATH = Path("/tmp/leopard-gf8-authoritative.lock")
BENCHMARK_CPU = 14
RESERVED_SIBLING = 30
ADDRESS_SPACE_BYTES = 256 * 1024 * 1024
CELL_REUSE = 8192
T95 = {
    3: 4.302652729911275,
    5: 2.7764451051977987,
    9: 2.306004135204166,
}
ROUND_ORDERS = (
    ("main", "candidate", "control", "control", "candidate", "main"),
    ("control", "main", "candidate", "candidate", "main", "control"),
    ("candidate", "control", "main", "main", "control", "candidate"),
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


def load_support() -> Any:
    path = Path(__file__).resolve().with_name("run_t8_two_block_abba.py")
    specification = importlib.util.spec_from_file_location(
        "t4_packed_terminal_t8_evidence_support", path)
    require(specification is not None and specification.loader is not None,
            f"cannot load evidence support: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


T8_SUPPORT = load_support()
MAIN_SUPPORT = T8_SUPPORT.SUPPORT
RUNNER_PATH = Path(__file__).resolve()
# Opt-in provenance and binary-equivalence hooks for configured wrappers.  The
# original T=4 campaign retains its historical checks when these remain empty.
RUNNER_DEPENDENCIES: tuple[Path, ...] = ()
EXPECTED_BINARY_SHA256: Mapping[str, str] | None = None
REQUIRE_NORMALIZED_FULL_FILE_EQUIVALENCE = False
# Heap layout in the benchmark executable can depend on argv allocation.  A
# configured campaign may require equal-length executable paths so that the
# comparison does not accidentally benchmark different buffer alignments.
REQUIRE_EQUAL_EXECUTABLE_PATH_LENGTHS = False
# Runtime attribution avoids separate-ELF placement bias by invoking hard
# links to one immutable executable and changing only a setup-time diagnostic
# selector before codec creation.
ALLOW_IDENTICAL_CANDIDATE_CONTROL = False
CONTROL_EXTRA_ARGUMENTS: tuple[str, ...] = ()
CONTROL_BUILD_MARKER: str | None = None
CANDIDATE_SCHEMA = "leopard2-benchmark-v5"
CONTROL_SCHEMA = "leopard2-benchmark-v5"
# Configured wrappers may prove that additional retained selectors have the
# same expected value in both binaries.  The primary MODE_SYMBOL remains the
# only selector masked by normalized_full_file_comparison().
AUXILIARY_MODE_EXPECTATIONS: Mapping[
    str, Mapping[str, int]
] = {}
# Performance evidence remains descriptive by default.  A configured wrapper
# may attach a deterministic promotion policy without changing the validity or
# retention of a completed negative campaign.
PROMOTION_EVALUATOR: Callable[
    [Sequence[Mapping[str, Any]]], Mapping[str, Any]
] | None = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def support_file_identities(paths: Sequence[Path]) -> list[dict[str, Any]]:
    resolved = [path.resolve(strict=True) for path in paths]
    require(len(resolved) == len(set(resolved)),
            "runner dependency paths are not unique")
    return [T8_SUPPORT.regular_file_identity(path) for path in resolved]


def verify_support_file_identities(
    expected: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    current = [
        T8_SUPPORT.regular_file_identity(Path(str(identity["path"])))
        for identity in expected
    ]
    require(current == list(expected),
            "a runner dependency changed during the campaign")
    return current


def output_bytes(value: bytes | str | None) -> bytes:
    if value is None:
        return b""
    if isinstance(value, bytes):
        return value
    return value.encode("utf-8", errors="replace")


def run_tool(arguments: Sequence[str]) -> str:
    completed = subprocess.run(
        list(arguments), env=CHILD_ENVIRONMENT,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        timeout=10.0, check=False)
    require(completed.returncode == 0,
            f"{arguments[0]} failed: " +
            completed.stderr.decode("utf-8", errors="replace")[-1000:])
    return completed.stdout.decode("utf-8", errors="strict")


def mode_word_identity(
    executable: Path,
    symbol_name: str | None = None,
) -> dict[str, Any]:
    """Map the retained T=4 selector symbol to its initialized ELF word."""
    selected_symbol = MODE_SYMBOL if symbol_name is None else symbol_name
    symbols = run_tool(("/usr/bin/readelf", "-sW", str(executable)))
    matches = []
    for line in symbols.splitlines():
        tokens = line.split()
        if len(tokens) >= 8 and tokens[-1] == selected_symbol:
            matches.append(tokens)
    require(len(matches) == 1,
            "T=4 diagnostic selector symbol is missing or ambiguous")
    symbol = matches[0]
    require(symbol[2] == "4" and
            symbol[3:6] == ["OBJECT", "LOCAL", "DEFAULT"],
            "T=4 diagnostic selector symbol metadata changed")
    try:
        address = int(symbol[1], 16)
        section_index = int(symbol[6])
    except ValueError as error:
        raise EvidenceError(
            "T=4 diagnostic selector symbol is not file-backed") from error

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
    require(selected is not None,
            "T=4 diagnostic selector section is absent")
    section_name, section_address, section_offset, section_size = selected
    require(section_address <= address and
            address + 4 <= section_address + section_size,
            "T=4 diagnostic selector lies outside its section")
    file_offset = section_offset + address - section_address
    with executable.open("rb") as source:
        require(source.read(6) == b"\x7fELF\x02\x01",
                "benchmark is not little-endian ELF64")
        source.seek(file_offset)
        payload = source.read(4)
    require(len(payload) == 4,
            "T=4 diagnostic selector word is truncated")
    return {
        "symbol": selected_symbol,
        "virtual_address": address,
        "section_index": section_index,
        "section_name": section_name,
        "file_offset": file_offset,
        "bytes_hex": payload.hex(),
        "value": int.from_bytes(payload, "little"),
    }


def elf_section_file_range(
    executable: Path,
    section_name: str,
) -> dict[str, Any]:
    """Locate one nonempty ELF section in the immutable file image."""
    sections = run_tool(("/usr/bin/readelf", "-SW", str(executable)))
    section_pattern = re.compile(
        r"^\s*\[\s*(\d+)\]\s+(\S+)\s+(\S+)\s+"
        r"([0-9A-Fa-f]+)\s+([0-9A-Fa-f]+)\s+([0-9A-Fa-f]+)\s+")
    matches: list[dict[str, Any]] = []
    for line in sections.splitlines():
        match = section_pattern.match(line)
        if match and match.group(2) == section_name:
            matches.append({
                "name": match.group(2),
                "section_index": int(match.group(1)),
                "section_type": match.group(3),
                "virtual_address": int(match.group(4), 16),
                "file_offset": int(match.group(5), 16),
                "size": int(match.group(6), 16),
            })
    require(len(matches) == 1 and matches[0]["size"] > 0,
            f"ELF section {section_name} is absent, ambiguous, or empty")
    selected = matches[0]
    require(selected["file_offset"] + selected["size"] <=
            executable.stat().st_size,
            f"ELF section {section_name} lies outside the file")
    return selected


def difference_ranges(offsets: Sequence[int]) -> list[dict[str, int]]:
    if not offsets:
        return []
    output: list[dict[str, int]] = []
    start = offsets[0]
    previous = start
    for offset in offsets[1:]:
        if offset != previous + 1:
            output.append({"file_offset": start,
                           "size": previous - start + 1})
            start = offset
        previous = offset
    output.append({"file_offset": start, "size": previous - start + 1})
    return output


def normalized_full_file_comparison(
    candidate: Path,
    control: Path,
    allowed_ranges: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Prove two files differ only in named, nonoverlapping byte ranges."""
    candidate_bytes = candidate.read_bytes()
    control_bytes = control.read_bytes()
    candidate_payload = bytearray(candidate_bytes)
    control_payload = bytearray(control_bytes)
    require(len(candidate_payload) == len(control_payload),
            "candidate and control full ELF sizes differ")
    file_size = len(candidate_payload)
    normalized_ranges = sorted((
        {
            "name": str(item["name"]),
            "file_offset": int(item["file_offset"]),
            "size": int(item["size"]),
        }
        for item in allowed_ranges
    ), key=lambda item: item["file_offset"])
    require(normalized_ranges, "no full-file difference ranges were allowed")
    previous_end = 0
    allowed = bytearray(file_size)
    for item in normalized_ranges:
        start = item["file_offset"]
        end = start + item["size"]
        require(item["size"] > 0 and start >= previous_end and
                end <= file_size,
                "allowed full-file difference ranges overlap or are invalid")
        allowed[start:end] = b"\x01" * item["size"]
        candidate_payload[start:end] = b"\0" * item["size"]
        control_payload[start:end] = b"\0" * item["size"]
        previous_end = end

    differing_offsets = [
        offset for offset, (candidate_byte, control_byte) in enumerate(
            zip(candidate_bytes, control_bytes))
        if candidate_byte != control_byte
    ]
    unauthorized = [offset for offset in differing_offsets if not allowed[offset]]
    require(not unauthorized,
            "candidate/control differ outside selector and GNU build-id ranges")
    candidate_digest = hashlib.sha256(candidate_payload).hexdigest()
    control_digest = hashlib.sha256(control_payload).hexdigest()
    require(candidate_digest == control_digest,
            "candidate/control normalized full-file SHA-256 differs")
    require(differing_offsets,
            "candidate/control full ELF files are unexpectedly identical")
    return {
        "file_size": file_size,
        "allowed_ranges": normalized_ranges,
        "difference_count": len(differing_offsets),
        "difference_ranges": difference_ranges(differing_offsets),
        "normalized_sha256": candidate_digest,
    }


def candidate_control_full_file_equivalence(
    candidate: Path,
    control: Path,
    mode_words: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    """Bind the only permitted candidate/control full-ELF differences."""
    selector_objects = {
        label: {
            "symbol": mode_words[label]["symbol"],
            "file_offset": mode_words[label]["file_offset"],
            "size": 4,
        }
        for label in ("candidate", "control")
    }
    require(selector_objects["candidate"] == selector_objects["control"],
            "candidate/control selector object layouts differ")
    build_ids = {
        label: elf_section_file_range(path, ".note.gnu.build-id")
        for label, path in (("candidate", candidate), ("control", control))
    }
    require(
        all(build_ids["candidate"][key] == build_ids["control"][key]
            for key in (
                "name", "section_index", "section_type", "virtual_address",
                "file_offset", "size")),
        "candidate/control GNU build-id section layouts differ")
    selector = {
        "name": "diagnostic_selector",
        "file_offset": selector_objects["candidate"]["file_offset"],
        "size": 4,
    }
    build_id = {
        "name": ".note.gnu.build-id",
        "file_offset": build_ids["candidate"]["file_offset"],
        "size": build_ids["candidate"]["size"],
    }
    comparison = normalized_full_file_comparison(
        candidate, control, (selector, build_id))
    comparison["selector_objects"] = selector_objects
    comparison["build_id_sections"] = build_ids
    return comparison


def campaign_cells() -> list[dict[str, Any]]:
    cells: list[dict[str, Any]] = []

    for k in range(4, 8):
        for r in (3, 4):
            for shard_bytes in (64, 128, 256):
                role = "preexisting_reference" \
                    if (k, r, shard_bytes) == (5, 4, 64) \
                    else "selected_target"
                cells.append({
                    "id": f"{role}-k{k}-r{r}-b{shard_bytes}-q1",
                    "K": k,
                    "R": r,
                    "bytes": shard_bytes,
                    "batch": 1,
                    "reuse": CELL_REUSE,
                    "role": role,
                    "candidate_selected": True,
                    "control_selected": False,
                })
    for r in (3, 4):
        cells.append({
            "id": f"selected_target-k4-r{r}-b1024-q1",
            "K": 4,
            "R": r,
            "bytes": 1024,
            "batch": 1,
            "reuse": CELL_REUSE,
            "role": "selected_target",
            "candidate_selected": True,
            "control_selected": False,
        })

    for r in (3, 4):
        for shard_bytes in (64, 256):
            cells.append({
                "id": f"unselected_neighbor-k8-r{r}-b{shard_bytes}-q1",
                "K": 8,
                "R": r,
                "bytes": shard_bytes,
                "batch": 1,
                "reuse": CELL_REUSE,
                "role": "unselected_neighbor",
                "candidate_selected": False,
                "control_selected": False,
            })
    for k in range(4, 8):
        for r in (3, 4):
            cells.append({
                "id": f"unselected_neighbor-k{k}-r{r}-b512-q1",
                "K": k,
                "R": r,
                "bytes": 512,
                "batch": 1,
                "reuse": CELL_REUSE,
                "role": "unselected_neighbor",
                "candidate_selected": False,
                "control_selected": False,
            })

    for index, cell in enumerate(cells):
        cell["seed"] = 0x74350000 + index
    require(len(cells) == 38 and
            sum(cell["role"] == "selected_target" for cell in cells) == 25 and
            sum(cell["role"] == "preexisting_reference"
                for cell in cells) == 1 and
            sum(cell["role"] == "unselected_neighbor"
                for cell in cells) == 12 and
            len({cell["id"] for cell in cells}) == len(cells),
            "T=4 packed-terminal campaign matrix is incomplete")
    return cells


def benchmark_command(
    implementation: str,
    executable: Path,
    cell: Mapping[str, Any],
    cpu: int,
    iterations: int,
    warmup: int,
) -> list[str]:
    common = [
        "/usr/bin/prlimit", f"--as={ADDRESS_SPACE_BYTES}",
        "/usr/bin/taskset", "-c", str(cpu), str(executable),
        "--k", str(cell["K"]), "--r", str(cell["R"]),
        "--bytes", str(cell["bytes"]), "--loss", "1",
        "--batch", "1", "--reuse", str(cell["reuse"]),
        "--iterations", str(iterations), "--warmup", str(warmup),
        "--threads", "1", "--seed", str(cell["seed"]), "--json", "-",
    ]
    if implementation == "main":
        return common
    command = common[:-2] + [
        "--profile", "high", "--field", "gf8", "--backend", "avx2",
        "--skip-legacy", "--retain-samples", "--attest-source",
    ]
    if implementation == "control":
        command.extend(CONTROL_EXTRA_ARGUMENTS)
    return command + ["--json", "-"]


def positive_encode_metric(result: Mapping[str, Any]) -> float:
    metrics = result.get("metrics")
    require(isinstance(metrics, dict), "benchmark metrics are absent")
    encode = metrics.get("encode_execution")
    require(isinstance(encode, dict), "encode metric is absent")
    median = encode.get("median_us_per_batch_call")
    require(isinstance(median, (int, float)) and
            not isinstance(median, bool) and
            math.isfinite(float(median)) and float(median) > 0,
            "encode median is not finite and positive")
    return float(median)


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
    expected_parameters = {
        "K": cell["K"],
        "R": cell["R"],
        "shard_bytes": cell["bytes"],
        "loss_count": 1,
        "batch": 1,
        "reuse": cell["reuse"],
        "iterations": iterations,
        "warmup": warmup,
        "thread_count": 1,
        "seed": cell["seed"],
    }
    parameters = result.get("parameters")
    require(isinstance(parameters, dict) and
            all(parameters.get(name) == value
                for name, value in expected_parameters.items()),
            "benchmark parameters differ from the frozen one-item cell")
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
            all(isinstance(digests.get(name), str) and
                re.fullmatch(r"[0-9a-f]{16}", digests[name]) is not None
                for name in (
                    "original_data", "transmitted_parity",
                    "recovered_originals")),
            "benchmark workload digests are incomplete")

    if implementation == "main":
        require(result.get("schema") == "leopard-main-benchmark-v1" and
                isinstance(correctness, dict) and
                correctness.get("round_trip") is True,
                "exact-main benchmark identity or round trip failed")
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
                "Leopard2 benchmark identity, backend, or round trip failed")
        build = result.get("build")
        require(isinstance(build, dict) and
                build.get("source_commit") == source_commit and
                build.get("source_tree") == source_tree and
                build.get("source_tracked_dirty") is False,
                "Leopard2 embedded source identity changed")
        require("prevalidated_batch_experiment" not in build and
                "high_t4_batch_selected" not in build,
                "benchmark is a prevalidated binding binary, not ordinary batch")
        if CONTROL_BUILD_MARKER is not None:
            if implementation == "control":
                require(build.get(CONTROL_BUILD_MARKER) is True,
                        "runtime attribution marker differs from control")
            else:
                require(build.get(CONTROL_BUILD_MARKER) is False,
                        "candidate unexpectedly disabled the terminal")
    return {
        "encode_us": positive_encode_metric(result),
        "digests": dict(digests),
        "schema": result["schema"],
    }


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
    failure_prefix = failure_output / (
        f"failure-{cell['id']}-{implementation}-{started_ns}")

    def persist_failure(stdout: bytes, stderr: bytes) -> None:
        T8_SUPPORT.write_bytes_exclusive(
            failure_prefix.with_suffix(".stdout"), stdout)
        T8_SUPPORT.write_bytes_exclusive(
            failure_prefix.with_suffix(".stderr"), stderr)

    try:
        completed = subprocess.run(
            command, env=CHILD_ENVIRONMENT,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            timeout=60.0, check=False)
    except subprocess.TimeoutExpired as error:
        persist_failure(
            output_bytes(error.stdout), output_bytes(error.stderr))
        raise EvidenceError(
            f"{implementation} timed out after 60 seconds") from error
    elapsed_ns = time.monotonic_ns() - started_ns
    if completed.returncode != 0:
        persist_failure(completed.stdout, completed.stderr)
        raise EvidenceError(
            f"{implementation} failed: " +
            completed.stderr.decode(
                "utf-8", errors="replace")[-1000:])
    if sha256(executable) != identity["sha256"]:
        persist_failure(completed.stdout, completed.stderr)
        raise EvidenceError(
            f"{implementation} binary changed after execution")
    try:
        result = json.loads(completed.stdout.decode("utf-8"))
        normalized = validate_result(
            implementation, result, cell, source_commit, source_tree,
            iterations, warmup)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        persist_failure(completed.stdout, completed.stderr)
        raise EvidenceError(
            f"{implementation} output is not one JSON object: {error}") \
            from error
    except Exception:
        persist_failure(completed.stdout, completed.stderr)
        raise
    return {
        "implementation": implementation,
        "command": command,
        "elapsed_ns": elapsed_ns,
        "stdout_sha256": hashlib.sha256(completed.stdout).hexdigest(),
        "stderr_sha256": hashlib.sha256(completed.stderr).hexdigest(),
        "normalized": normalized,
        "result": result,
    }


def append_invocations(
    order: Sequence[str],
    destination: list[dict[str, Any]],
    invoke: Callable[[str], dict[str, Any]],
) -> None:
    """Retain each successful child before starting the next one."""
    for label in order:
        destination.append(invoke(label))


def confidence_interval(round_log_ratios: Sequence[float]) -> dict[str, Any]:
    require(len(round_log_ratios) in T95,
            "round count has no predeclared Student-t threshold")
    center = statistics.mean(round_log_ratios)
    half_width = T95[len(round_log_ratios)] * \
        statistics.stdev(round_log_ratios) / \
        math.sqrt(len(round_log_ratios))
    return {
        "speedup": math.exp(center),
        "speedup_definition": "baseline_encode_us / candidate_encode_us",
        "ci95": [
            math.exp(center - half_width),
            math.exp(center + half_width),
        ],
        "round_log_ratios": list(round_log_ratios),
    }


def analyze(
    cell: Mapping[str, Any],
    rounds: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    reference = rounds[0]["invocations"][0]["normalized"]["digests"]
    contrasts: dict[str, list[float]] = {
        "control": [],
        "main": [],
    }
    for round_value in rounds:
        require(round_value["isolation"]["accepted"] is True,
                "contaminated round cannot be analyzed")
        invocations = round_value["invocations"]
        require(all(item["normalized"]["digests"] == reference
                    for item in invocations),
                "candidate/control/main workload digests differ")
        candidate = [
            item["normalized"]["encode_us"] for item in invocations
            if item["implementation"] == "candidate"
        ]
        require(len(candidate) == 2,
                "round lacks two candidate observations")
        candidate_log = statistics.mean(math.log(value)
                                        for value in candidate)
        for label in contrasts:
            baseline = [
                item["normalized"]["encode_us"] for item in invocations
                if item["implementation"] == label
            ]
            require(len(baseline) == 2,
                    f"round lacks two {label} observations")
            contrasts[label].append(
                statistics.mean(math.log(value) for value in baseline) -
                candidate_log)
    return {
        "cell": dict(cell),
        "digests": reference,
        "candidate_vs_control": confidence_interval(contrasts["control"]),
        "candidate_vs_main": confidence_interval(contrasts["main"]),
    }


def acquire_global_lock() -> int:
    descriptor = os.open(LOCK_PATH, os.O_RDWR | os.O_CREAT, 0o600)
    try:
        fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        os.close(descriptor)
        raise EvidenceError(f"authoritative lock is busy: {LOCK_PATH}")
    return descriptor


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", required=True, type=Path,
                        help="clean ordinary bench_leopard2 candidate")
    parser.add_argument("--control", required=True, type=Path,
                        help="same-source compile-time-disabled control")
    parser.add_argument("--main", required=True, type=Path,
                        help="exact-main bench_leopard2 executable")
    parser.add_argument("--source-commit", required=True,
                        help="embedded candidate/control Git commit")
    parser.add_argument("--source-tree", required=True,
                        help="embedded candidate/control Git tree")
    parser.add_argument("--output", required=True, type=Path,
                        help="new exclusive evidence directory")
    parser.add_argument("--cpu", type=int, default=BENCHMARK_CPU)
    parser.add_argument("--sibling", type=int, default=RESERVED_SIBLING)
    parser.add_argument("--rounds", type=int, choices=tuple(T95), default=3)
    parser.add_argument("--iterations", type=int, default=15)
    parser.add_argument("--warmup", type=int, default=64)
    return parser.parse_args()


def verify_frozen_identity(
    label: str,
    identity: Mapping[str, Any],
) -> dict[str, Any]:
    current = T8_SUPPORT.file_identity(Path(str(identity["path"])))
    require(current == identity,
            f"{label} binary identity changed during the campaign")
    return current


def executable_path_lengths(
    identities: Mapping[str, Mapping[str, Any]],
) -> dict[str, int]:
    """Return path lengths and reject allocation-biased lane names."""
    lengths = {
        label: len(str(identity["path"]))
        for label, identity in identities.items()
    }
    require(len(set(lengths.values())) == 1,
            "candidate, control, and main executable path lengths must match")
    return lengths


def main() -> int:
    options = parse_arguments()
    require(options.cpu == BENCHMARK_CPU and
            options.sibling == RESERVED_SIBLING,
            "this campaign is frozen to CPU14 and its sibling CPU30")
    require(options.iterations >= 3 and options.warmup >= 1,
            "insufficient benchmark repetitions")
    require(re.fullmatch(r"[0-9a-f]{40}", options.source_commit) is not None and
            re.fullmatch(r"[0-9a-f]{40}", options.source_tree) is not None,
            "source commit and tree must be lowercase SHA-1 identities")
    require(not options.output.exists(), "output path already exists")
    options.output.mkdir(parents=True)
    cells = campaign_cells()
    orders = tuple(
        ROUND_ORDERS[index % len(ROUND_ORDERS)]
        for index in range(options.rounds))
    raw: dict[str, Any] = {
        "schema": SCHEMA,
        "created_utc": MAIN_SUPPORT.utc_now(),
        "source_commit": options.source_commit,
        "source_tree": options.source_tree,
        "main_commit": MAIN_COMMIT,
        "cpu": options.cpu,
        "reserved_sibling": options.sibling,
        "address_space_bytes": ADDRESS_SPACE_BYTES,
        "iterations": options.iterations,
        "warmup": options.warmup,
        "rounds": options.rounds,
        "round_orders": [list(order) for order in orders],
        "cells": [],
    }
    lock_descriptor: int | None = None
    try:
        raw["runner"] = T8_SUPPORT.file_identity(RUNNER_PATH)
        if RUNNER_DEPENDENCIES:
            raw["runner_dependencies_before"] = \
                support_file_identities(RUNNER_DEPENDENCIES)
        lock_descriptor = acquire_global_lock()
        identities = {
            "candidate": T8_SUPPORT.file_identity(options.candidate),
            "control": T8_SUPPORT.file_identity(options.control),
            "main": T8_SUPPORT.file_identity(options.main),
        }
        raw["identities_before"] = identities
        require(len({identity["path"] for identity in identities.values()}) == 3,
                "candidate, control, and main paths are not distinct")
        if ALLOW_IDENTICAL_CANDIDATE_CONTROL:
            require(
                identities["candidate"]["sha256"] ==
                    identities["control"]["sha256"] and
                identities["candidate"]["sha256"] !=
                    identities["main"]["sha256"],
                "runtime attribution requires one shared candidate/control "
                "binary")
            candidate_status = Path(
                str(identities["candidate"]["path"])).stat()
            control_status = Path(str(identities["control"]["path"])).stat()
            require(
                (candidate_status.st_dev, candidate_status.st_ino) ==
                (control_status.st_dev, control_status.st_ino),
                "runtime attribution requires candidate/control hard links "
                "to one inode")
            raw["candidate_control_shared_inode"] = {
                "device": candidate_status.st_dev,
                "inode": candidate_status.st_ino,
            }
        else:
            require(
                len({identity["sha256"] for identity in identities.values()}) == 3,
                "candidate, control, and main binaries are not distinct")
        if REQUIRE_EQUAL_EXECUTABLE_PATH_LENGTHS:
            raw["executable_path_lengths"] = \
                executable_path_lengths(identities)
        if EXPECTED_BINARY_SHA256 is not None:
            expected_hashes = dict(EXPECTED_BINARY_SHA256)
            require(set(expected_hashes) == set(identities) and
                    all(isinstance(value, str) and
                        re.fullmatch(r"[0-9a-f]{64}", value) is not None
                        for value in expected_hashes.values()),
                    "expected frozen binary SHA-256 mapping is invalid")
            raw["expected_binary_sha256"] = expected_hashes
            require(all(identities[label]["sha256"] == expected_hashes[label]
                        for label in identities),
                    "a benchmark binary does not match its frozen SHA-256")

        executable_sections = {
            label: T8_SUPPORT.executable_sections_identity(
                Path(str(identities[label]["path"])))
            for label in ("candidate", "control")
        }
        raw["candidate_control_executable_sections"] = executable_sections
        require(executable_sections["candidate"]["sections"] ==
                executable_sections["control"]["sections"],
                "candidate and control executable instruction sections differ")

        if ALLOW_IDENTICAL_CANDIDATE_CONTROL:
            mode_words = {
                "shared_binary_default": mode_word_identity(
                    Path(str(identities["candidate"]["path"])))
            }
        else:
            mode_words = {
                label: mode_word_identity(Path(str(identities[label]["path"])))
                for label in ("candidate", "control")
            }
        raw["mode_words"] = mode_words
        if ALLOW_IDENTICAL_CANDIDATE_CONTROL:
            require(mode_words["shared_binary_default"]["value"] == 1,
                    "shared binary does not default to the candidate route")
        else:
            require(mode_words["candidate"]["value"] == 1 and
                    mode_words["control"]["value"] == 2,
                    "candidate/control T=4 selector words are not "
                    "enabled/disabled")
            for key in (
                    "symbol", "virtual_address", "section_index",
                    "section_name", "file_offset"):
                require(
                    mode_words["candidate"][key] ==
                    mode_words["control"][key],
                    "candidate/control T=4 selector layouts differ")
        auxiliary_mode_words: dict[str, dict[str, Any]] = {}
        for symbol_name in sorted(AUXILIARY_MODE_EXPECTATIONS):
            expected_values = AUXILIARY_MODE_EXPECTATIONS[symbol_name]
            require(
                set(expected_values) == {"candidate", "control"} and
                all(isinstance(value, int) and not isinstance(value, bool)
                    for value in expected_values.values()),
                "auxiliary selector expectations are invalid")
            words = {
                label: mode_word_identity(
                    Path(str(identities[label]["path"])), symbol_name)
                for label in ("candidate", "control")
            }
            for label in ("candidate", "control"):
                require(words[label]["value"] == expected_values[label],
                        f"{label} auxiliary selector value changed")
            for key in (
                    "symbol", "virtual_address", "section_index",
                    "section_name", "file_offset"):
                require(words["candidate"][key] == words["control"][key],
                        "candidate/control auxiliary selector layouts differ")
            auxiliary_mode_words[symbol_name] = words
        if auxiliary_mode_words:
            raw["auxiliary_mode_words"] = auxiliary_mode_words
        if REQUIRE_NORMALIZED_FULL_FILE_EQUIVALENCE:
            require(not ALLOW_IDENTICAL_CANDIDATE_CONTROL,
                    "shared-inode attribution does not need file masking")
            raw["candidate_control_normalized_full_file"] = \
                candidate_control_full_file_equivalence(
                    Path(str(identities["candidate"]["path"])),
                    Path(str(identities["control"]["path"])), mode_words)

        require(set(os.sched_getaffinity(0)) == {options.cpu},
                "runner must be singleton-pinned to CPU14")
        sibling_text = Path(
            f"/sys/devices/system/cpu/cpu{options.cpu}/topology/"
            "thread_siblings_list").read_text(encoding="ascii")
        require(MAIN_SUPPORT.parse_cpu_list(sibling_text) ==
                {options.cpu, options.sibling},
                "CPU14 and CPU30 are not one SMT pair")
        raw["host"] = {
            "runner_affinity": sorted(os.sched_getaffinity(0)),
            "benchmark_cpu": MAIN_SUPPORT.cpu_policy_identity(options.cpu),
            "reserved_sibling":
                MAIN_SUPPORT.cpu_policy_identity(options.sibling),
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
            require(
                presample["delta"]["benchmark_cpu"]["nonidle_jiffies"] == 0 and
                presample["delta"]["reserved_sibling"]["nonidle_jiffies"] == 0,
                "CPU14/30 pair was not quiet during the presample")

            for cell_index, cell in enumerate(cells):
                cell_raw = {"cell": dict(cell), "rounds": []}
                raw["cells"].append(cell_raw)
                for round_index, order in enumerate(orders):
                    before_cpu = MAIN_SUPPORT.cpu_stat_snapshot(options.cpu)
                    before_sibling = \
                        MAIN_SUPPORT.cpu_stat_snapshot(options.sibling)
                    before_ns = time.monotonic_ns()
                    round_raw: dict[str, Any] = {
                        "round": round_index,
                        "order": list(order),
                        "invocations": [],
                        "isolation": None,
                    }
                    cell_raw["rounds"].append(round_raw)
                    append_invocations(
                        order, round_raw["invocations"], lambda label: run_one(
                            label, identities[label], cell, options.cpu,
                            options.source_commit, options.source_tree,
                            options.iterations, options.warmup,
                            options.output))
                    isolation = MAIN_SUPPORT.isolation_record(
                        options.cpu, options.sibling, pair_lease,
                        before_ns, time.monotonic_ns(), before_cpu,
                        MAIN_SUPPORT.cpu_stat_snapshot(options.cpu),
                        before_sibling,
                        MAIN_SUPPORT.cpu_stat_snapshot(options.sibling))
                    round_raw["isolation"] = isolation
                    require(isolation["accepted"],
                            f"contaminated {cell['id']} round {round_index}")
                print(f"{cell_index + 1}/{len(cells)} {cell['id']}",
                      file=sys.stderr, flush=True)

        post_identities = {
            label: verify_frozen_identity(label, identity)
            for label, identity in identities.items()
        }
        raw["identities_after"] = post_identities
        if RUNNER_DEPENDENCIES:
            raw["runner_dependencies_after"] = \
                verify_support_file_identities(
                    raw["runner_dependencies_before"])
        analyses = [
            analyze(item["cell"], item["rounds"])
            for item in raw["cells"]
        ]
        candidate_slower_than_control = [
            item["cell"]["id"] for item in analyses
            if item["candidate_vs_control"]["speedup"] < 1.0
        ]
        credible_candidate_control_regressions = [
            item["cell"]["id"] for item in analyses
            if item["candidate_vs_control"]["ci95"][1] < 1.0
        ]
        candidate_slower_than_main = [
            item["cell"]["id"] for item in analyses
            if item["candidate_vs_main"]["speedup"] < 1.0
        ]
        credible_candidate_main_regressions = [
            item["cell"]["id"] for item in analyses
            if item["candidate_vs_main"]["ci95"][1] < 1.0
        ]

        promotion_evaluation: dict[str, Any] | None = None
        if PROMOTION_EVALUATOR is not None:
            evaluated = PROMOTION_EVALUATOR(analyses)
            require(isinstance(evaluated, Mapping),
                    "promotion evaluator did not return a mapping")
            promotion_evaluation = dict(evaluated)
            require(
                isinstance(promotion_evaluation.get("passed"), bool),
                "promotion evaluator omitted its boolean result")
            raw["promotion_evaluation"] = promotion_evaluation

        raw["completed_utc"] = MAIN_SUPPORT.utc_now()
        T8_SUPPORT.write_exclusive(options.output / "raw.json", raw)
        summary = {
            "schema": SUMMARY_SCHEMA,
            "status": "completed",
            "performance_gate_applied":
                promotion_evaluation is not None,
            "source_commit": options.source_commit,
            "source_tree": options.source_tree,
            "main_commit": MAIN_COMMIT,
            "cell_count": len(analyses),
            "process_count": sum(
                len(round_value["invocations"])
                for cell_value in raw["cells"]
                for round_value in cell_value["rounds"]),
            "all_digests_matched": True,
            "all_rounds_accepted": True,
            "candidate_slower_than_control_cells":
                candidate_slower_than_control,
            "credible_candidate_control_regressions":
                credible_candidate_control_regressions,
            "candidate_slower_than_main_cells": candidate_slower_than_main,
            "credible_candidate_main_regressions":
                credible_candidate_main_regressions,
            "cells": analyses,
            "binary_sha256": {
                label: identity["sha256"]
                for label, identity in identities.items()
            },
            "candidate_control_executable_sections_sha256":
                executable_sections["candidate"]["combined_sha256"],
            "mode_words": mode_words,
            "raw_sha256": sha256(options.output / "raw.json"),
        }
        if promotion_evaluation is not None:
            summary["promotion_evaluation"] = promotion_evaluation
            summary["promotion_gate_passed"] = \
                promotion_evaluation["passed"]
        for optional_key in (
                "candidate_control_normalized_full_file",
                "candidate_control_shared_inode", "executable_path_lengths",
                "expected_binary_sha256", "runner_dependencies_before",
                "runner_dependencies_after", "auxiliary_mode_words"):
            if optional_key in raw:
                summary[optional_key] = raw[optional_key]
        T8_SUPPORT.write_exclusive(options.output / "summary.json", summary)
        print(json.dumps({
            "status": summary["status"],
            "cells": summary["cell_count"],
            "processes": summary["process_count"],
            "candidate_slower_than_control_cells":
                candidate_slower_than_control,
            "candidate_slower_than_main_cells": candidate_slower_than_main,
        }, sort_keys=True))
        return 0
    except Exception as error:
        if (RUNNER_DEPENDENCIES and
                "runner_dependencies_before" in raw and
                "runner_dependencies_after" not in raw):
            try:
                raw["runner_dependencies_after"] = \
                    verify_support_file_identities(
                        raw["runner_dependencies_before"])
            except Exception as dependency_error:
                raw["runner_dependencies_after_error"] = {
                    "type": type(dependency_error).__name__,
                    "message": str(dependency_error),
                }
        raw["failed_utc"] = MAIN_SUPPORT.utc_now()
        raw["failure"] = {
            "type": type(error).__name__,
            "message": str(error),
        }
        T8_SUPPORT.write_exclusive(options.output / "failure.json", raw)
        print(f"evidence rejected: {error}", file=sys.stderr)
        return 1
    finally:
        if lock_descriptor is not None:
            os.close(lock_descriptor)


if __name__ == "__main__":
    raise SystemExit(main())
