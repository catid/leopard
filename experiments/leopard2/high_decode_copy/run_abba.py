#!/usr/bin/env python3
"""Collect fail-closed same-binary Algorithm 5 copy/no-copy ABBA evidence."""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import shlex
import statistics
import subprocess
import sys
import tempfile
import time
import traceback
from pathlib import Path
from types import ModuleType
from typing import Any, Mapping, Sequence


RAW_SCHEMA = "leopard2-high-decode-copy-raw/v1"
MANIFEST_SCHEMA = "leopard2-high-decode-copy-manifest/v1"
FAILURE_SCHEMA = "leopard2-high-decode-copy-failure/v1"
MATRIX_SCHEMA = "leopard2-high-decode-copy-matrix/v1"
RUNNER_RELATIVE = "experiments/leopard2/high_decode_copy/run_abba.py"
CONTRACT_RELATIVE = "experiments/leopard2/high_decode_copy/benchmark_contract.py"
MATRIX_RELATIVE = "experiments/leopard2/high_decode_copy/matrix.json"
TARGET_NAME = "bench_leopard2_high_decode_copy_attribution"


def load_module(name: str, path: Path) -> ModuleType:
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


ROOT = Path(__file__).resolve().parents[3]
SUPPORT = load_module(
    "leopard2_high_copy_exact_main_support",
    ROOT / "experiments/leopard2/main_compare/run_abba.py")
CONTRACT = load_module(
    "leopard2_high_copy_benchmark_contract",
    ROOT / CONTRACT_RELATIVE)
EvidenceError = SUPPORT.EvidenceError


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvidenceError(message)


def load_matrix(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, ValueError) as error:
        raise EvidenceError(f"cannot read matrix: {error}") from error
    require(isinstance(value, dict) and value.get("schema") == MATRIX_SCHEMA and
            set(value) == {"schema", "cells", "round_orders", "comparison_policy"},
            "high-decode copy matrix schema or shape changed")
    cells = value["cells"]
    require(isinstance(cells, list) and len(cells) == 16,
            "high-decode copy matrix must contain exactly 16 bounded cells")
    expected_keys = {"id", "backend", "field", "k", "r", "shard_bytes",
                     "losses", "seed", "workspace", "role",
                     "exact_main_eligible"}
    seen: set[str] = set()
    coverage: set[tuple[str, str, str]] = set()
    for cell in cells:
        require(isinstance(cell, dict) and set(cell) == expected_keys,
                "high-decode copy cell shape changed")
        identifier = cell["id"]
        require(isinstance(identifier, str) and identifier and identifier not in seen,
                "high-decode copy cell ID is empty or duplicated")
        seen.add(identifier)
        require(cell["backend"] in ("scalar", "ssse3", "avx2") and
                cell["field"] in ("gf8", "gf16") and
                cell["workspace"] in ("materialized", "tiled") and
                cell["role"] in ("target", "neighbor", "full-block", "tail") and
                all(type(cell[name]) is int and cell[name] > 0
                    for name in ("k", "r", "shard_bytes", "losses", "seed")) and
                cell["losses"] <= min(cell["k"], cell["r"]) and
                type(cell["exact_main_eligible"]) is bool,
                f"invalid high-decode copy cell {identifier}")
        padded_recovery = 1 << (cell["r"] - 1).bit_length()
        parent = 1 << (cell["k"] + padded_recovery - 1).bit_length()
        expected_field = "gf8" if parent <= 256 else "gf16"
        tail = cell["shard_bytes"] % 64
        full_dimensions = (128, 128, 128) if cell["field"] == "gf8" \
            else (256, 256, 256)
        require(parent <= 65536 and cell["field"] == expected_field and
                (cell["role"] == "tail") == (tail != 0) and
                (cell["role"] != "full-block" or
                 (cell["k"], cell["r"], cell["losses"]) == full_dimensions) and
                (cell["field"] != "gf16" or cell["shard_bytes"] % 2 == 0) and
                cell["exact_main_eligible"] == (tail == 0),
                f"cell {identifier} misclassifies its field or exact-main/tail eligibility")
        coverage.add((cell["field"], cell["workspace"], cell["role"]))
    require(coverage == {(field, workspace, role)
                         for field in ("gf8", "gf16")
                         for workspace in ("materialized", "tiled")
                         for role in ("target", "neighbor", "full-block", "tail")},
            "matrix lost field/workspace/target-neighbor-full-block-tail coverage")
    orders = value["round_orders"]
    require(orders == [
        ["copy-fallback", "no-copy", "no-copy", "copy-fallback"],
        ["no-copy", "copy-fallback", "copy-fallback", "no-copy"],
        ["copy-fallback", "no-copy", "no-copy", "copy-fallback"],
    ], "ABBA/BAAB/ABBA order changed")
    policy = value["comparison_policy"]
    require(isinstance(policy, dict) and
            policy.get("exact_main_commit") == SUPPORT.MAIN_COMMIT and
            policy.get("tail_classification") ==
                "same-source-only because exact Leopard main requires shard_bytes divisible by 64",
            "matrix exact-main policy changed")
    return value


def checked_output(arguments: list[str]) -> str:
    output = SUPPORT.run_checked(
        arguments, environment=SUPPORT.CHILD_ENVIRONMENT, timeout=30,
        max_stdout=16 * 1024 * 1024, max_stderr=1024 * 1024)
    try:
        return output.decode("utf-8", errors="strict").strip()
    except UnicodeError as error:
        raise EvidenceError(f"command output is not UTF-8: {error}") from error


def find_build_root(binary: Path) -> Path:
    for candidate in (binary.parent, *binary.parents):
        if (candidate / "CMakeCache.txt").is_file():
            return candidate
    raise EvidenceError("diagnostic benchmark has no enclosing CMake cache")


def build_identity(source_root: Path, binary: Path, hook_archive: Path) -> dict[str, Any]:
    source_root = source_root.resolve(strict=True)
    binary = binary.resolve(strict=True)
    hook_archive = hook_archive.resolve(strict=True)
    require(os.access(binary, os.X_OK), "diagnostic benchmark is not executable")
    build = find_build_root(binary)
    require(hook_archive.is_relative_to(build),
            "hook archive and diagnostic benchmark are not in one build tree")
    cache_path = build / "CMakeCache.txt"
    cache = SUPPORT.parse_cmake_cache(cache_path)
    require(Path(cache.get("CMAKE_HOME_DIRECTORY", "")).resolve() == source_root and
            cache.get("CMAKE_BUILD_TYPE") == "Release" and
            cache.get("LEO2_BUILD_TESTS") == "ON" and
            cache.get("LEO2_BUILD_BENCHMARKS") == "ON" and
            cache.get("LEO2_ENABLE_CUDA") == "OFF" and
            cache.get("LEOPARD_ENABLE_GF8") == "ON" and
            cache.get("LEOPARD_ENABLE_GF16") == "ON",
            "diagnostic CMake cache is not a Release tests/hooks build")
    link_recipe = build / "CMakeFiles" / f"{TARGET_NAME}.dir" / "link.txt"
    require(link_recipe.is_file(), "diagnostic benchmark link recipe is missing")
    link_text = link_recipe.read_text(encoding="utf-8", errors="strict")
    tokens = [token for line in link_text.splitlines() for token in shlex.split(line)]
    linked_paths: list[Path] = []
    for token in tokens:
        if token.startswith("-"):
            continue
        candidate = Path(token)
        candidate = candidate if candidate.is_absolute() else build / candidate
        try:
            linked_paths.append(candidate.resolve(strict=True))
        except OSError:
            continue
    require(linked_paths.count(hook_archive) == 1 and
            not any(path.name in ("libleopard.a", "leopard.lib")
                    for path in linked_paths),
            "diagnostic benchmark is not linked exactly once and exclusively "
            "to the hooks archive")
    commands_path = build / "compile_commands.json"
    require(commands_path.is_file(), "diagnostic build omitted compile_commands.json")
    commands = json.loads(commands_path.read_text(encoding="utf-8"))
    matching = [entry for entry in commands
                if isinstance(entry, dict) and
                TARGET_NAME in str(entry.get("output", "")) and
                str(entry.get("file", "")).endswith("/bench/leopard2/benchmark.cpp")]
    require(len(matching) == 1, "cannot identify the diagnostic benchmark compile action")
    compile_text = matching[0].get("command") or " ".join(matching[0].get("arguments", []))
    object_path = Path(str(matching[0].get("output", "")))
    object_path = object_path if object_path.is_absolute() else build / object_path
    require("LEO2_ENABLE_TEST_HOOKS=1" in compile_text and
            "LEO2_HIGH_DECODE_COPY_ATTRIBUTION=1" in compile_text and
            object_path.resolve(strict=True) in linked_paths,
            "diagnostic benchmark compile action lacks its two private gates")
    nm = Path("/usr/bin/nm").resolve(strict=True)
    required_symbols = tuple(
        f"leopard::{field}::{name}" for field in ("ff8", "ff16")
        for name in ("TestOnlySetHighDecodeCopyFallback",
                     "TestOnlyHighDecodeCopyFallbackEnabled"))
    hook_symbols = checked_output(
        [str(nm), "-C", "--defined-only", str(hook_archive)])
    binary_symbols = checked_output(
        [str(nm), "-C", "--defined-only", str(binary)])
    require(all(name in hook_symbols for name in required_symbols),
            "hook archive does not contain both fields' copy selectors")
    require(all(name in binary_symbols for name in required_symbols),
            "diagnostic executable did not retain both fields' copy selectors")
    return {
        "build_root": str(build),
        "cache": SUPPORT.artifact_identity(cache_path, "build_metadata"),
        "compile_commands": SUPPORT.artifact_identity(
            commands_path, "build_metadata"),
        "link_recipe": SUPPORT.artifact_identity(link_recipe, "build_metadata"),
        "binary": SUPPORT.artifact_identity(binary, "executable"),
        "hook_archive": SUPPORT.artifact_identity(hook_archive, "archive"),
        "nm": SUPPORT.artifact_identity(nm, "executable"),
        "selector_symbols_present_in_hook_archive_and_binary": True,
    }


def snapshot(source_root: Path, commit: str, binary: Path,
             hook_archive: Path, matrix_path: Path) -> dict[str, Any]:
    require(Path(__file__).resolve().relative_to(source_root.resolve()).as_posix() ==
            RUNNER_RELATIVE, "runner is not executing from its canonical path")
    return {
        "source": SUPPORT.git_identity(source_root, commit, detached=False),
        "build": build_identity(source_root, binary, hook_archive),
        "runner": SUPPORT.artifact_identity(Path(__file__), "source_file"),
        "contract": SUPPORT.artifact_identity(ROOT / CONTRACT_RELATIVE, "source_file"),
        "matrix": SUPPORT.artifact_identity(matrix_path, "source_file"),
        "support": SUPPORT.artifact_identity(
            ROOT / "experiments/leopard2/main_compare/run_abba.py", "source_file"),
        "runtime": SUPPORT.runtime_closure(Path("/usr/bin/ldd"), binary),
    }


def benchmark_command(binary: Path, cell: Mapping[str, Any], mode: str,
                      cpu: int, reuse: int, iterations: int,
                      warmup: int) -> list[str]:
    return [
        "/usr/bin/taskset", "-c", str(cpu), str(binary),
        "--k", str(cell["k"]), "--r", str(cell["r"]),
        "--profile", "high", "--field", str(cell["field"]),
        "--backend", str(cell["backend"]), "--bytes", str(cell["shard_bytes"]),
        "--loss", str(cell["losses"]), "--batch", "1", "--reuse", str(reuse),
        "--iterations", str(iterations), "--warmup", str(warmup),
        "--threads", "1", "--seed", str(cell["seed"]),
        "--force-specialized", "--force-" + str(cell["workspace"]),
        "--skip-legacy", "--retain-samples", "--report-decode-path",
        "--high-evaluator-mode", mode, "--json", "-",
    ]


def validate_result(document: object, cell: Mapping[str, Any], mode: str) -> dict[str, Any]:
    result = CONTRACT.validate_document(
        document, mode=mode, workspace=cell["workspace"], field=cell["field"],
        tail_bytes=cell["shard_bytes"] % 64,
        evaluator="mature" if cell["role"] == "full-block" else None)
    parameters = result["parameters"]
    resolved = result["resolved"]
    require(parameters.get("K") == cell["k"] and parameters.get("R") == cell["r"] and
            parameters.get("requested_backend") == cell["backend"] and
            parameters.get("shard_bytes") == cell["shard_bytes"] and
            parameters.get("loss_count") == cell["losses"] and
            parameters.get("seed") == cell["seed"] and
            resolved.get("backend") == cell["backend"],
            "benchmark output does not attest the matrix cell")
    return result


def run_child(binary: Path, cell: Mapping[str, Any], mode: str, cpu: int,
              sibling: int,
              reuse: int, iterations: int, warmup: int, timeout: float,
              round_index: int, slot: int) -> dict[str, Any]:
    command = benchmark_command(
        binary, cell, mode, cpu, reuse, iterations, warmup)
    before_cpu = SUPPORT.cpu_stat_snapshot(cpu)
    before_sibling = SUPPORT.cpu_stat_snapshot(sibling)
    started = time.monotonic_ns()
    completed = SUPPORT.run_process_bounded(
        command, environment=SUPPORT.CHILD_ENVIRONMENT, timeout=timeout,
        max_stdout=16 * 1024 * 1024, max_stderr=1024 * 1024)
    ended = time.monotonic_ns()
    after_cpu = SUPPORT.cpu_stat_snapshot(cpu)
    after_sibling = SUPPORT.cpu_stat_snapshot(sibling)
    cpu_delta = SUPPORT.cpu_stat_delta(before_cpu, after_cpu)
    sibling_delta = SUPPORT.cpu_stat_delta(before_sibling, after_sibling)
    require(completed.returncode == 0 and completed.stderr == b"" and
            cpu_delta["nonidle_jiffies"] > 0 and
            sibling_delta["nonidle_jiffies"] == 0,
            "child failed, emitted stderr, did no timed CPU work, or used the SMT sibling")
    require(len(completed.stdout) <= 16 * 1024 * 1024,
            "benchmark JSON exceeds its retained bound")
    try:
        document = json.loads(completed.stdout.decode("utf-8", errors="strict"))
    except (UnicodeError, ValueError) as error:
        raise EvidenceError(f"benchmark stdout is not JSON: {error}") from error
    result = validate_result(document, cell, mode)
    return {
        "cell": cell["id"], "mode": mode, "round": round_index, "slot": slot,
        "command": command, "command_sha256": SUPPORT.sha256_bytes(
            SUPPORT.canonical_bytes(command)),
        "started_monotonic_ns": started, "ended_monotonic_ns": ended,
        "cpu_before": before_cpu, "cpu_after": after_cpu, "cpu_delta": cpu_delta,
        "sibling_before": before_sibling, "sibling_after": after_sibling,
        "sibling_delta": sibling_delta,
        "stdout_sha256": SUPPORT.sha256_bytes(completed.stdout),
        "stdout_size": len(completed.stdout), "result": result,
    }


def grouped_records(records: Sequence[Mapping[str, Any]],
                    cells: Sequence[Mapping[str, Any]]) -> dict[str, list[Mapping[str, Any]]]:
    result = {cell["id"]: [] for cell in cells}
    for record in records:
        identifier = record.get("cell")
        require(identifier in result, "record refers to an unknown cell")
        result[identifier].append(record)
    return result


def analyze(records: Sequence[Mapping[str, Any]], matrix: Mapping[str, Any]) -> dict[str, Any]:
    grouped = grouped_records(records, matrix["cells"])
    output: dict[str, Any] = {}
    metrics = {
        "decode_execution": ("metrics", "decode_execution", "median_us"),
        "decode_amortized": (
            "metrics", "decode_amortized_at_reuse", "derived_median_us_per_batch_call"),
    }
    for cell in matrix["cells"]:
        cell_records = grouped[cell["id"]]
        require(len(cell_records) == 12, "cell does not have three four-slot rounds")
        documents = {mode: [record["result"] for record in cell_records
                            if record["mode"] == mode]
                     for mode in ("no-copy", "copy-fallback")}
        CONTRACT.validate_pair(
            documents["no-copy"][0], documents["copy-fallback"][0],
            workspace=cell["workspace"], field=cell["field"],
            tail_bytes=cell["shard_bytes"] % 64,
            evaluator="mature" if cell["role"] == "full-block" else None)
        digest = documents["no-copy"][0]["workload_digests"]
        missing = documents["no-copy"][0]["parameters"]["missing_original_indices"]
        require(all(document["workload_digests"] == digest and
                    document["parameters"]["missing_original_indices"] == missing
                    for values in documents.values() for document in values),
                "ABBA roles changed workload digest or loss set")
        cell_analysis: dict[str, Any] = {
            "workload_digests": digest,
            "missing_original_indices": missing,
            "exact_main_classification": "eligible" if cell["exact_main_eligible"]
                else "same-source-only-tail",
        }
        for metric, path in metrics.items():
            round_logs = []
            for round_index, order in enumerate(matrix["round_orders"]):
                round_records = sorted(
                    (record for record in cell_records if record["round"] == round_index),
                    key=lambda record: record["slot"])
                require([record["mode"] for record in round_records] == order,
                        "record order differs from the signed round order")
                values = []
                for record in round_records:
                    value: object = record["result"]
                    for name in path:
                        require(isinstance(value, dict), "timing metric path is malformed")
                        value = value[name]
                    require(isinstance(value, (int, float)) and not isinstance(value, bool) and
                            math.isfinite(float(value)) and float(value) > 0,
                            "timing metric is not finite and positive")
                    values.append(float(value))
                if order[0] == "copy-fallback":
                    contrasts = (math.log(values[0] / values[1]),
                                 math.log(values[3] / values[2]))
                else:
                    contrasts = (math.log(values[1] / values[0]),
                                 math.log(values[2] / values[3]))
                round_logs.append(statistics.fmean(contrasts))
            mean = statistics.fmean(round_logs)
            standard_error = statistics.stdev(round_logs) / math.sqrt(3.0)
            margin = 4.302652729911275 * standard_error
            cell_analysis[metric] = {
                "ratio_orientation": "copy_fallback_time_over_no_copy_time",
                "geometric_speedup": math.exp(mean),
                "ci95_lower": math.exp(mean - margin),
                "ci95_upper": math.exp(mean + margin),
                "independent_round_log_contrasts": round_logs,
                "independent_round_count": 3,
                "degrees_of_freedom": 2,
            }
        output[cell["id"]] = cell_analysis
    return output


def validate_campaign_isolation(
    isolation: object, campaign: Mapping[str, Any]
) -> dict[str, Any]:
    """Validate the canonical nested scheduler record used by main_compare."""
    cpu = campaign.get("cpu")
    sibling = campaign.get("reserved_sibling")
    require(isinstance(cpu, int) and not isinstance(cpu, bool) and
            isinstance(sibling, int) and not isinstance(sibling, bool),
            "campaign CPU identities are invalid")
    return SUPPORT.validate_isolation(isolation, cpu, sibling)


def validate_raw(raw: object, output: Path, *, current: bool) -> dict[str, Any]:
    SUPPORT.verify_signature(raw, "high-decode copy raw evidence")
    require(isinstance(raw, dict) and raw.get("schema") == RAW_SCHEMA and
            raw.get("validity_is_independent_of_speed") is True,
            "raw high-decode copy evidence schema changed")
    matrix = raw.get("matrix")
    require(matrix == load_matrix(ROOT / MATRIX_RELATIVE),
            "raw matrix differs from the canonical checked-in matrix")
    campaign = raw.get("campaign")
    require(isinstance(campaign, dict) and
            campaign.get("round_orders") == matrix["round_orders"] and
            campaign.get("child_environment") == SUPPORT.CHILD_ENVIRONMENT and
            campaign.get("thread_count") == 1,
            "raw campaign contract is incomplete")
    records = raw.get("records")
    invocation_count = sum(len(order) for order in matrix["round_orders"]) * \
        len(matrix["cells"])
    require(campaign.get("invocation_count") == invocation_count and
            isinstance(records, list) and len(records) == invocation_count,
            "raw campaign does not contain the exact bounded invocation count")
    binary = Path(raw["identities_initial"]["build"]["binary"]["path"])
    for cell in matrix["cells"]:
        matching = [record for record in records if record.get("cell") == cell["id"]]
        require(len(matching) == 12, "raw cell record count changed")
        for record in matching:
            command = benchmark_command(
                binary, cell, record["mode"], campaign["cpu"], campaign["reuse"],
                campaign["iterations"], campaign["warmup"])
            require(record.get("command") == command and
                    record.get("command_sha256") == SUPPORT.sha256_bytes(
                        SUPPORT.canonical_bytes(command)) and
                    record.get("cpu_delta", {}).get("cpu") == campaign["cpu"] and
                    record.get("sibling_delta", {}).get("cpu") ==
                        campaign["reserved_sibling"] and
                    record.get("sibling_delta", {}).get("nonidle_jiffies") == 0 and
                    record.get("cpu_delta", {}).get("nonidle_jiffies", 0) > 0,
                    "raw command or per-child CPU isolation attestation differs")
            validate_result(record.get("result"), cell, record["mode"])
    require(raw.get("analysis") == analyze(records, matrix),
            "raw analysis differs from retained invocations")
    validate_campaign_isolation(raw.get("isolation"), campaign)
    require(raw.get("identities_initial") == raw.get("identities_final"),
            "source/build identity changed during the campaign")
    if current:
        source = raw["identities_initial"]["source"]
        rebuilt = snapshot(
            Path(source["path"]), source["head"], binary,
            Path(raw["identities_initial"]["build"]["hook_archive"]["path"]),
            ROOT / MATRIX_RELATIVE)
        require(rebuilt == raw["identities_initial"],
                "current source/build closure differs from retained evidence")
    return raw["analysis"]


def run_campaign(options: argparse.Namespace) -> int:
    output = options.output.resolve()
    require(not output.exists(), "output directory already exists")
    output.mkdir(mode=0o700, parents=True)
    matrix_path = (options.source_root.resolve() / MATRIX_RELATIVE)
    matrix = load_matrix(matrix_path)
    require(options.iterations >= 3 and options.reuse >= 1 and options.warmup >= 1 and
            math.isfinite(options.timeout) and 0 < options.timeout <= 3600,
            "timing controls are outside their bounded policy")
    initial = None
    records: list[dict[str, Any]] = []
    isolation = None
    reservation = None
    pair_lease = None
    host_initial = None
    campaign = {
        "cpu": options.cpu, "reserved_sibling": options.reserved_sibling,
        "reuse": options.reuse, "iterations": options.iterations,
        "warmup": options.warmup, "thread_count": 1,
        "timeout_seconds": options.timeout,
        "round_orders": matrix["round_orders"],
        "child_environment": dict(SUPPORT.CHILD_ENVIRONMENT),
        "invocation_count": sum(len(order) for order in matrix["round_orders"]) *
            len(matrix["cells"]),
    }
    try:
        allowed, housekeeping = SUPPORT.validate_topology(
            options.cpu, options.reserved_sibling)
        host_initial = SUPPORT.host_identity(options.cpu, options.reserved_sibling, allowed)
        pair_guard = SUPPORT.PairLease(options.cpu, options.reserved_sibling)
        with SUPPORT.Reservation(
            options.reservation_file, options.cpu, options.reserved_sibling
        ) as reservation, pair_guard as pair_lease:
            os.sched_setaffinity(0, housekeeping)
            initial = snapshot(
                options.source_root, options.source_commit,
                options.binary, options.hook_archive, matrix_path)
            before_ns = time.monotonic_ns()
            before_cpu = SUPPORT.cpu_stat_snapshot(options.cpu)
            before_sibling = SUPPORT.cpu_stat_snapshot(options.reserved_sibling)
            for cell in matrix["cells"]:
                for round_index, order in enumerate(matrix["round_orders"]):
                    for slot, mode in enumerate(order):
                        record = run_child(
                            options.binary.resolve(strict=True), cell, mode,
                            options.cpu, options.reserved_sibling, options.reuse,
                            options.iterations, options.warmup, options.timeout,
                            round_index, slot)
                        pair_guard.validate_current()
                        SUPPORT.validate_reservation_current(reservation)
                        records.append(record)
            after_cpu = SUPPORT.cpu_stat_snapshot(options.cpu)
            after_sibling = SUPPORT.cpu_stat_snapshot(options.reserved_sibling)
            isolation = SUPPORT.isolation_record(
                options.cpu, options.reserved_sibling, pair_lease,
                before_ns, time.monotonic_ns(), before_cpu, after_cpu,
                before_sibling, after_sibling)
            require(isolation["accepted"] is True,
                    "reserved SMT sibling performed work during the campaign")
            final = snapshot(
                options.source_root, options.source_commit,
                options.binary, options.hook_archive, matrix_path)
            require(final == initial, "source/build closure changed during the campaign")
            host_final = SUPPORT.host_identity(
                options.cpu, options.reserved_sibling, allowed)
            require(host_final == host_initial,
                    "host topology/frequency policy changed during the campaign")
            raw = SUPPORT.signed({
                "schema": RAW_SCHEMA,
                "created_utc": SUPPORT.utc_now(),
                "validity_is_independent_of_speed": True,
                "matrix": matrix,
                "campaign": campaign,
                "host_initial": host_initial,
                "host_final": host_final,
                "reservation": reservation,
                "pair_lease": pair_lease,
                "isolation": isolation,
                "identities_initial": initial,
                "identities_final": final,
                "records": records,
                "analysis": analyze(records, matrix),
            })
            validate_raw(raw, output, current=True)
            raw_path = output / "raw.json"
            SUPPORT.write_json_exclusive(raw_path, raw)
            manifest = SUPPORT.signed({
                "schema": MANIFEST_SCHEMA,
                "created_utc": SUPPORT.utc_now(),
                "valid": True,
                "validity_is_independent_of_speed": True,
                "raw": {"path": "raw.json", "size": raw_path.stat().st_size,
                        "sha256": SUPPORT.sha256_file(raw_path),
                        "payload_digest": raw["digest"]},
                "matrix": matrix,
                "campaign": campaign,
                "identities": initial,
                "isolation": isolation,
                "analysis": raw["analysis"],
            })
            SUPPORT.write_json_exclusive(output / "manifest.json", manifest)
    except Exception as error:
        failure = SUPPORT.signed({
            "schema": FAILURE_SCHEMA, "created_utc": SUPPORT.utc_now(),
            "valid": False, "error_type": type(error).__name__, "error": str(error),
            "matrix": matrix, "campaign": campaign, "host_initial": host_initial,
            "reservation": reservation, "pair_lease": pair_lease,
            "isolation": isolation, "identities_initial": initial,
            "record_count": len(records), "traceback": traceback.format_exc(),
        })
        if not (output / "failure.json").exists():
            SUPPORT.write_json_exclusive(output / "failure.json", failure)
        raise
    print(output / "manifest.json")
    return 0


def verify_campaign(options: argparse.Namespace) -> int:
    path = options.manifest.resolve(strict=True)
    manifest = json.loads(path.read_text(encoding="utf-8"))
    SUPPORT.verify_signature(manifest, "high-decode copy manifest")
    require(manifest.get("schema") == MANIFEST_SCHEMA and manifest.get("valid") is True,
            "manifest is not valid high-decode copy evidence")
    raw_info = manifest.get("raw")
    require(isinstance(raw_info, dict), "manifest has no raw identity")
    raw_path = SUPPORT.safe_evidence_path(path.parent, raw_info.get("path"))
    require(raw_path.is_file() and raw_path.stat().st_size == raw_info.get("size") and
            SUPPORT.sha256_file(raw_path) == raw_info.get("sha256"),
            "raw high-decode copy bundle identity differs")
    raw = json.loads(raw_path.read_text(encoding="utf-8"))
    analysis = validate_raw(raw, path.parent, current=not options.no_current_input_check)
    require(raw.get("digest") == raw_info.get("payload_digest") and
            manifest.get("matrix") == raw.get("matrix") and
            manifest.get("campaign") == raw.get("campaign") and
            manifest.get("identities") == raw.get("identities_initial") and
            manifest.get("isolation") == raw.get("isolation") and
            manifest.get("analysis") == analysis,
            "manifest differs from its raw high-decode copy bundle")
    print("high-decode copy/no-copy ABBA evidence verified")
    return 0


def build_smoke(options: argparse.Namespace) -> int:
    identity = build_identity(
        options.source_root, options.binary, options.hook_archive)
    require(identity.get(
        "selector_symbols_present_in_hook_archive_and_binary") is True,
        "diagnostic selector symbol gate did not pass")
    print("high-decode copy diagnostic build identity verified")
    return 0


def self_test() -> None:
    matrix = load_matrix(ROOT / MATRIX_RELATIVE)
    wrong_field = json.loads(json.dumps(matrix))
    tail = next(cell for cell in wrong_field["cells"]
                if cell["id"] == "gf8-mat-tail")
    tail["k"] = 193
    wrong_full_block = json.loads(json.dumps(matrix))
    full_block = next(cell for cell in wrong_full_block["cells"]
                      if cell["id"] == "gf8-mat-full-block")
    full_block["losses"] = 127
    with tempfile.TemporaryDirectory(prefix="leo2-high-copy-matrix-") as root:
        for index, mutation in enumerate((wrong_field, wrong_full_block)):
            path = Path(root) / f"matrix-{index}.json"
            path.write_text(json.dumps(mutation), encoding="utf-8")
            try:
                load_matrix(path)
            except EvidenceError:
                pass
            else:
                raise AssertionError(
                    "invalid field/full-block matrix mutation was accepted")
    command_a = benchmark_command(
        Path("/tmp/bench"), matrix["cells"][0], "copy-fallback", 7, 8, 9, 2)
    command_b = benchmark_command(
        Path("/tmp/bench"), matrix["cells"][0], "no-copy", 7, 8, 9, 2)
    difference = [(left, right) for left, right in zip(command_a, command_b)
                  if left != right]
    require(difference == [("copy-fallback", "no-copy")],
            "same-binary roles differ outside the explicit selector")
    def cpu_stat(cpu: int, *, user: int, idle: int) -> dict[str, Any]:
        fields = {
            "user": user, "nice": 0, "system": 0, "idle": idle,
            "iowait": 0, "irq": 0, "softirq": 0, "steal": 0,
        }
        return {"cpu": cpu, "fields": fields,
                "total_jiffies": sum(fields.values())}

    payload = SUPPORT.pair_lease_payload(0, 1)
    lease = {
        "device": 1, "directory_device": 1, "directory_inode": 2,
        "inode": 3, "lock": "exclusive_nonblocking_pair_wide",
        "path": str(SUPPORT.pair_lease_directory() /
                    SUPPORT.pair_lease_name(0, 1)),
        "payload": payload,
        "sha256": SUPPORT.sha256_bytes(SUPPORT.canonical_bytes(payload)),
    }
    nested = SUPPORT.isolation_record(
        0, 1, lease, 1_000, 2_000,
        cpu_stat(0, user=100, idle=100),
        cpu_stat(0, user=110, idle=110),
        cpu_stat(1, user=100, idle=100),
        cpu_stat(1, user=100, idle=120))
    envelope = SUPPORT.signed({
        "schema": "high-decode-copy-isolation-fixture-v1",
        "isolation": nested,
    })
    SUPPORT.verify_signature(envelope, "nested isolation fixture")
    validate_campaign_isolation(
        envelope["isolation"], {"cpu": 0, "reserved_sibling": 1})
    flat = {"accepted": True, "reserved_sibling_nonidle_jiffies": 0}
    try:
        validate_campaign_isolation(flat, {"cpu": 0, "reserved_sibling": 1})
    except SUPPORT.EvidenceError:
        pass
    else:
        raise AssertionError("obsolete flat isolation fixture was accepted")
    CONTRACT.self_test()
    print("high-decode copy/no-copy runner self-test passed: "
          "16 cells, 192 invocations, canonical full-block, field, and nested isolation")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    commands.add_parser("self-test")
    build = commands.add_parser("build-smoke")
    build.add_argument("--source-root", required=True, type=Path)
    build.add_argument("--binary", required=True, type=Path)
    build.add_argument("--hook-archive", required=True, type=Path)
    run = commands.add_parser("run")
    run.add_argument("--source-root", required=True, type=Path)
    run.add_argument("--source-commit", required=True)
    run.add_argument("--binary", required=True, type=Path)
    run.add_argument("--hook-archive", required=True, type=Path)
    run.add_argument("--reservation-file", required=True, type=Path)
    run.add_argument("--output", required=True, type=Path)
    run.add_argument("--cpu", required=True, type=int)
    run.add_argument("--reserved-sibling", required=True, type=int)
    run.add_argument("--reuse", type=int, default=8)
    run.add_argument("--iterations", type=int, default=9)
    run.add_argument("--warmup", type=int, default=2)
    run.add_argument("--timeout", type=float, default=120.0)
    verify = commands.add_parser("verify")
    verify.add_argument("--manifest", required=True, type=Path)
    verify.add_argument("--no-current-input-check", action="store_true")
    return result


def main(arguments: Sequence[str] | None = None) -> int:
    options = parser().parse_args(arguments)
    if options.command == "self-test":
        self_test()
        return 0
    if options.command == "run":
        return run_campaign(options)
    if options.command == "build-smoke":
        return build_smoke(options)
    return verify_campaign(options)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (EvidenceError, CONTRACT.ContractError, OSError, ValueError,
            subprocess.SubprocessError) as error:
        print(f"high-decode copy evidence error: {error}", file=sys.stderr)
        raise SystemExit(1)
