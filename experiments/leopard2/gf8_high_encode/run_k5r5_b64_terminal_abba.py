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
LOCK_PATH = Path("/tmp/leopard-gf8-authoritative.lock")
T95_DF2 = 4.302652729911275
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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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
    common = [
        "/usr/bin/prlimit", "--as=201326592",
        "/usr/bin/taskset", "-c", str(cpu), str(executable),
        "--k", str(cell["K"]), "--r", str(cell["R"]),
        "--bytes", str(cell["bytes"]), "--loss", "1",
        "--batch", str(cell["batch"]), "--reuse", str(cell["reuse"]),
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
        "K": cell["K"], "R": cell["R"],
        "shard_bytes": cell["bytes"], "loss_count": 1,
        "batch": cell["batch"], "reuse": cell["reuse"],
        "iterations": iterations, "warmup": warmup,
        "thread_count": 1, "seed": cell["seed"],
    }
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
        require(result.get("schema") == "leopard-main-benchmark-v1" and
                isinstance(correctness, dict) and
                correctness.get("round_trip") is True,
                "exact-main identity or round trip failed")
        build = result.get("build")
        require(isinstance(build, dict) and
                build.get("main_source_commit") == MAIN_COMMIT,
                "exact-main source identity changed")
    else:
        require(result.get("schema") == "leopard2-benchmark-v5" and
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
    completed = subprocess.run(
        command, env=CHILD_ENVIRONMENT, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, timeout=60.0, check=False)
    elapsed_ns = time.monotonic_ns() - started_ns
    if completed.returncode != 0:
        prefix = failure_output / f"failure-{implementation}-{started_ns}"
        T8_SUPPORT.write_bytes_exclusive(
            prefix.with_suffix(".stdout"), completed.stdout)
        T8_SUPPORT.write_bytes_exclusive(
            prefix.with_suffix(".stderr"), completed.stderr)
        raise EvidenceError(
            f"{implementation} failed: " +
            completed.stderr.decode("utf-8", errors="replace")[-1000:])
    require(sha256(executable) == identity["sha256"],
            f"{implementation} binary changed after execution")
    try:
        result = json.loads(completed.stdout.decode("utf-8"))
        normalized = validate_result(
            implementation, result, cell, source_commit, source_tree,
            iterations, warmup)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise EvidenceError(
            f"{implementation} output is not one JSON object: {error}") \
            from error
    return {
        "implementation": implementation,
        "elapsed_ns": elapsed_ns,
        "stdout_sha256": hashlib.sha256(completed.stdout).hexdigest(),
        "stderr_sha256": hashlib.sha256(completed.stderr).hexdigest(),
        "normalized": normalized,
        "result": result,
    }


def confidence_interval(values: Sequence[float]) -> dict[str, Any]:
    require(len(values) == 3, "three independent round contrasts are required")
    center = statistics.mean(values)
    half_width = T95_DF2 * statistics.stdev(values) / math.sqrt(len(values))
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
    labels = ("control", "main") if cell["role"] == "target" \
        else ("control",)
    contrasts: dict[str, list[float]] = {label: [] for label in labels}
    for round_value in rounds:
        require(round_value["isolation"]["accepted"] is True,
                "contaminated round cannot be analyzed")
        invocations = round_value["invocations"]
        require(all(item["normalized"]["digests"] == reference
                    for item in invocations),
                "candidate/control/main workload digests differ")
        candidate = [
            item["normalized"]["encode_us"] for item in invocations
            if item["implementation"] == "candidate"]
        require(len(candidate) == 2, "round lacks two candidate observations")
        candidate_log = statistics.mean(math.log(value) for value in candidate)
        for label in labels:
            baseline = [
                item["normalized"]["encode_us"] for item in invocations
                if item["implementation"] == label]
            require(len(baseline) == 2,
                    f"round lacks two {label} observations")
            contrasts[label].append(
                statistics.mean(math.log(value) for value in baseline) -
                candidate_log)
    output = {
        "cell": dict(cell),
        "digests": reference,
        "control_over_candidate": confidence_interval(contrasts["control"]),
    }
    if "main" in contrasts:
        output["main_over_candidate"] = confidence_interval(contrasts["main"])
    return output


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
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--control", required=True, type=Path)
    parser.add_argument("--main", required=True, type=Path)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--source-tree", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cpu", required=True, type=int)
    parser.add_argument("--sibling", required=True, type=int)
    parser.add_argument("--iterations", type=int, default=15)
    parser.add_argument("--warmup", type=int, default=64)
    return parser.parse_args()


def main() -> int:
    options = parse_arguments()
    require(options.iterations >= 3 and options.warmup >= 1,
            "insufficient benchmark repetitions")
    require(not options.output.exists(), "output path already exists")
    options.output.mkdir(parents=True)
    raw: dict[str, Any] = {
        "schema": SCHEMA,
        "created_utc": MAIN_SUPPORT.utc_now(),
        "source_commit": options.source_commit,
        "source_tree": options.source_tree,
        "main_commit": MAIN_COMMIT,
        "cpu": options.cpu,
        "reserved_sibling": options.sibling,
        "iterations": options.iterations,
        "warmup": options.warmup,
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
        raw["identities"] = identities
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
        mode_words = {
            name: mode_word_identity(Path(str(identities[name]["path"])))
            for name in ("candidate", "control")
        }
        raw["mode_words"] = mode_words
        require(mode_words["candidate"]["value"] == 1 and
                mode_words["control"]["value"] == 2,
                "candidate/control diagnostic selectors are not enabled/disabled")
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
                orders = TARGET_ORDER if cell["role"] == "target" \
                    else NEIGHBOR_ORDER
                for round_index, order in enumerate(orders):
                    before_cpu = MAIN_SUPPORT.cpu_stat_snapshot(options.cpu)
                    before_sibling = MAIN_SUPPORT.cpu_stat_snapshot(options.sibling)
                    before_ns = time.monotonic_ns()
                    invocations = [
                        run_one(
                            label, identities[label], cell, options.cpu,
                            options.source_commit, options.source_tree,
                            options.iterations, options.warmup, options.output)
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
                print(f"{cell_index + 1}/{len(all_cells)} {cell['id']}",
                      file=sys.stderr, flush=True)

        analyses = [
            analyze(item["cell"], item["rounds"]) for item in raw["cells"]]
        target = next(item for item in analyses
                      if item["cell"]["role"] == "target")
        target_control_failure = \
            target["control_over_candidate"]["ci95"][0] < TARGET_CONTROL_FLOOR
        target_main_failure = \
            target["main_over_candidate"]["ci95"][0] < TARGET_MAIN_FLOOR
        credible_neighbor_regressions = [
            item["cell"]["id"] for item in analyses
            if item["cell"]["role"] == "neighbor" and
            item["control_over_candidate"]["ci95"][1] < NEIGHBOR_FLOOR
        ]
        raw["completed_utc"] = MAIN_SUPPORT.utc_now()
        T8_SUPPORT.write_exclusive(options.output / "raw.json", raw)
        summary = {
            "schema": SUMMARY_SCHEMA,
            "status": "accepted" if not (
                target_control_failure or target_main_failure or
                credible_neighbor_regressions) else "rejected",
            "source_commit": options.source_commit,
            "source_tree": options.source_tree,
            "main_commit": MAIN_COMMIT,
            "cell_count": len(analyses),
            "process_count": sum(
                len(round_value["invocations"])
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
            "credible_neighbor_regressions": credible_neighbor_regressions,
            "cells": analyses,
            "binary_sha256": {
                name: identity["sha256"] for name, identity in identities.items()},
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
            "target_control": target["control_over_candidate"],
            "target_main": target["main_over_candidate"],
            "credible_neighbor_regressions": credible_neighbor_regressions,
        }, sort_keys=True))
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
