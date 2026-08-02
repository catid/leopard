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
from typing import Any, Mapping, Sequence


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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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


def mode_word_identity(executable: Path) -> dict[str, Any]:
    """Map the retained T=4 selector symbol to its initialized ELF word."""
    symbols = run_tool(("/usr/bin/readelf", "-sW", str(executable)))
    matches = []
    for line in symbols.splitlines():
        tokens = line.split()
        if len(tokens) >= 8 and tokens[-1] == MODE_SYMBOL:
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
        "symbol": MODE_SYMBOL,
        "virtual_address": address,
        "section_index": section_index,
        "section_name": section_name,
        "file_offset": file_offset,
        "bytes_hex": payload.hex(),
        "value": int.from_bytes(payload, "little"),
    }


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
    return common[:-2] + [
        "--profile", "high", "--field", "gf8", "--backend", "avx2",
        "--skip-legacy", "--retain-samples", "--attest-source",
        "--json", "-",
    ]


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
                len(digests[name]) == 16
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
        require(result.get("schema") == "leopard2-benchmark-v5" and
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
    require(sha256(executable) == identity["sha256"],
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
        "runner": T8_SUPPORT.file_identity(Path(__file__)),
        "cells": [],
    }
    lock_descriptor: int | None = None
    try:
        lock_descriptor = acquire_global_lock()
        identities = {
            "candidate": T8_SUPPORT.file_identity(options.candidate),
            "control": T8_SUPPORT.file_identity(options.control),
            "main": T8_SUPPORT.file_identity(options.main),
        }
        raw["identities_before"] = identities
        require(len({identity["path"] for identity in identities.values()}) == 3,
                "candidate, control, and main paths are not distinct")
        require(len({identity["sha256"] for identity in identities.values()}) == 3,
                "candidate, control, and main binaries are not distinct")

        executable_sections = {
            label: T8_SUPPORT.executable_sections_identity(
                Path(str(identities[label]["path"])))
            for label in ("candidate", "control")
        }
        raw["candidate_control_executable_sections"] = executable_sections
        require(executable_sections["candidate"]["sections"] ==
                executable_sections["control"]["sections"],
                "candidate and control executable instruction sections differ")

        mode_words = {
            label: mode_word_identity(Path(str(identities[label]["path"])))
            for label in ("candidate", "control")
        }
        raw["mode_words"] = mode_words
        require(mode_words["candidate"]["value"] == 1 and
                mode_words["control"]["value"] == 2,
                "candidate/control T=4 selector words are not enabled/disabled")
        for key in (
                "symbol", "virtual_address", "section_index",
                "section_name", "file_offset"):
            require(mode_words["candidate"][key] == mode_words["control"][key],
                    "candidate/control T=4 selector layouts differ")

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
                for round_index, order in enumerate(orders):
                    before_cpu = MAIN_SUPPORT.cpu_stat_snapshot(options.cpu)
                    before_sibling = \
                        MAIN_SUPPORT.cpu_stat_snapshot(options.sibling)
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
                    cell_raw["rounds"].append({
                        "round": round_index,
                        "order": list(order),
                        "invocations": invocations,
                        "isolation": isolation,
                    })
                    require(isolation["accepted"],
                            f"contaminated {cell['id']} round {round_index}")
                raw["cells"].append(cell_raw)
                print(f"{cell_index + 1}/{len(cells)} {cell['id']}",
                      file=sys.stderr, flush=True)

        post_identities = {
            label: verify_frozen_identity(label, identity)
            for label, identity in identities.items()
        }
        raw["identities_after"] = post_identities
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

        raw["completed_utc"] = MAIN_SUPPORT.utc_now()
        T8_SUPPORT.write_exclusive(options.output / "raw.json", raw)
        summary = {
            "schema": SUMMARY_SCHEMA,
            "status": "completed",
            "performance_gate_applied": False,
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
