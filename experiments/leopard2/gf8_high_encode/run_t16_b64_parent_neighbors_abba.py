#!/usr/bin/env python3
"""Screen T16/B64 layout neighbors against the exact parent executable."""

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


SCHEMA = "leopard2-gf8-t16-parent-neighbors-abba/v1"
SUMMARY_SCHEMA = "leopard2-gf8-t16-parent-neighbors-summary/v1"
LOCK_PATH = Path("/tmp/leopard-gf8-authoritative.lock")
NEW_COMMIT = "58f409ed42b298666a1b56b869db1beb9fc47ac4"
NEW_TREE = "065c09619e2c3b1f3e683704fed5d57f69241770"
PARENT_COMMIT = "068659e1dbcc8a2c75fef1348f1f916a8557599f"
PARENT_TREE = "e042f3d15f9e13e886d28ad1c18f161109950e13"
REGRESSION_FLOOR = 1.0 / 1.02
ENCODE_INTERNAL = "_ZL14EncodeInternalPK10leo2_codecmPKPKvPKPvS6_mS3_mb"
IFFT_ENCODER = (
    "_ZN7leopard3ff8L16IFFT_DIT_EncoderERKNS_7backend3OpsEmPKPKvjPPvSA_"
    "jPKhbNS0_19EncoderIFFTCallsiteEb")
T16_KERNEL = "_ZN7leopard7backendL23AVX2FF8HighEncodeT16B64EPKPKvPKPv"
ORDER = (
    ("parent", "new", "new", "parent"),
    ("new", "parent", "parent", "new"),
    ("parent", "new", "new", "parent"),
)
ENVIRONMENT = {
    "LANG": "C", "LC_ALL": "C", "PATH": "/usr/bin:/bin", "TZ": "UTC",
    "OMP_DYNAMIC": "FALSE", "OMP_NUM_THREADS": "1",
    "OMP_PLACES": "cores", "OMP_PROC_BIND": "TRUE",
}


class EvidenceError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvidenceError(message)


def load_base() -> Any:
    path = Path(__file__).resolve().with_name(
        "run_k5r5_b64_terminal_abba.py")
    spec = importlib.util.spec_from_file_location(
        "t16_parent_neighbor_evidence_base", path)
    require(spec is not None and spec.loader is not None,
            f"cannot load evidence support: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


BASE = load_base()
MAIN_SUPPORT = BASE.MAIN_SUPPORT
RUNNER = Path(__file__).resolve()
DEPENDENCIES = (
    RUNNER, Path(BASE.__file__).resolve(),
    Path(BASE.T8_SUPPORT.__file__).resolve(), Path(MAIN_SUPPORT.__file__).resolve())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_tool(arguments: Sequence[str]) -> str:
    completed = subprocess.run(
        list(arguments), env=ENVIRONMENT, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, timeout=10, check=False)
    require(completed.returncode == 0,
            f"{arguments[0]} failed: " +
            completed.stderr.decode(errors="replace")[-1000:])
    return completed.stdout.decode("utf-8")


def symbol_identity(executable: Path, symbol_name: str,
                    required: bool = True) -> dict[str, Any] | None:
    matches = []
    for line in run_tool(("/usr/bin/readelf", "-sW", str(executable))).splitlines():
        tokens = line.split()
        if len(tokens) >= 8 and tokens[-1] == symbol_name:
            matches.append(tokens)
    if not matches and not required:
        return None
    require(len(matches) == 1, f"hot symbol is absent or ambiguous: {symbol_name}")
    symbol = matches[0]
    try:
        address = int(symbol[1], 16)
        size = int(symbol[2])
        section_index = int(symbol[6])
    except ValueError as error:
        raise EvidenceError(f"hot symbol is not file-backed: {symbol_name}") from error
    require(size > 0, f"hot symbol has zero size: {symbol_name}")
    pattern = re.compile(
        r"^\s*\[\s*(\d+)\]\s+(\S+)\s+(\S+)\s+"
        r"([0-9A-Fa-f]+)\s+([0-9A-Fa-f]+)\s+([0-9A-Fa-f]+)\s+")
    section = None
    for line in run_tool(("/usr/bin/readelf", "-SW", str(executable))).splitlines():
        match = pattern.match(line)
        if match and int(match.group(1)) == section_index:
            section = {
                "name": match.group(2), "address": int(match.group(4), 16),
                "offset": int(match.group(5), 16), "size": int(match.group(6), 16)}
            break
    require(section is not None, f"hot-symbol section is absent: {symbol_name}")
    require(section["address"] <= address and
            address + size <= section["address"] + section["size"],
            f"hot symbol exceeds its section: {symbol_name}")
    file_offset = section["offset"] + address - section["address"]
    with executable.open("rb") as source:
        source.seek(file_offset)
        payload = source.read(size)
    require(len(payload) == size, f"hot symbol is truncated: {symbol_name}")
    return {
        "symbol": symbol_name, "section": section["name"], "size": size,
        "file_offset": file_offset, "sha256": hashlib.sha256(payload).hexdigest()}


def hot_identities(executable: Path, is_new: bool) -> dict[str, Any]:
    return {
        "encode_internal": symbol_identity(executable, ENCODE_INTERNAL),
        "ff8_ifft_encoder": symbol_identity(executable, IFFT_ENCODER),
        "t16_b64_kernel": symbol_identity(executable, T16_KERNEL, is_new),
    }


def cells() -> list[dict[str, Any]]:
    values = (
        ("byte-k16-r16-b63-q1", 16, 16, 63),
        ("byte-k16-r16-b65-q1", 16, 16, 65),
        ("byte-k16-r16-b128-q1", 16, 16, 128),
        ("shape-k15-r16-b64-q1", 15, 16, 64),
        ("shape-k16-r15-b64-q1", 16, 15, 64),
        ("transform-t8-k8-r8-b64-q1", 8, 8, 64),
        ("transform-t32-k32-r32-b64-q1", 32, 32, 64),
    )
    return [{
        "id": name, "K": k, "R": r, "bytes": shard_bytes,
        "batch": 1, "reuse": 8192, "seed": 0x50313600 + index,
    } for index, (name, k, r, shard_bytes) in enumerate(values)]


def command(executable: Path, cell: Mapping[str, Any], cpu: int,
            iterations: int, warmup: int) -> list[str]:
    return [
        "/usr/bin/prlimit", "--as=201326592", "/usr/bin/taskset", "-c",
        str(cpu), str(executable), "--k", str(cell["K"]), "--r", str(cell["R"]),
        "--bytes", str(cell["bytes"]), "--loss", "1", "--batch", "1",
        "--reuse", str(cell["reuse"]), "--iterations", str(iterations),
        "--warmup", str(warmup), "--threads", "1", "--seed", str(cell["seed"]),
        "--profile", "high", "--field", "gf8", "--backend", "avx2",
        "--skip-legacy", "--retain-samples", "--attest-source", "--json", "-"]


def run_one(label: str, identity: Mapping[str, Any], cell: Mapping[str, Any],
            cpu: int, iterations: int, warmup: int, failures: Path) -> dict[str, Any]:
    executable = Path(str(identity["path"]))
    require(sha256(executable) == identity["sha256"],
            f"{label} binary changed before execution")
    invocation = command(executable, cell, cpu, iterations, warmup)
    completed = subprocess.run(
        invocation, env=ENVIRONMENT, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, timeout=60, check=False)
    if completed.returncode != 0:
        stamp = time.monotonic_ns()
        BASE.T8_SUPPORT.write_bytes_exclusive(
            failures / f"failure-{label}-{stamp}.stdout", completed.stdout)
        BASE.T8_SUPPORT.write_bytes_exclusive(
            failures / f"failure-{label}-{stamp}.stderr", completed.stderr)
        raise EvidenceError(f"{label} failed: " +
                            completed.stderr.decode(errors="replace")[-1000:])
    require(sha256(executable) == identity["sha256"],
            f"{label} binary changed after execution")
    result = json.loads(completed.stdout)
    expected_commit, expected_tree = (
        (NEW_COMMIT, NEW_TREE) if label == "new" else (PARENT_COMMIT, PARENT_TREE))
    expected_parameters = {
        "K": cell["K"], "R": cell["R"], "shard_bytes": cell["bytes"],
        "loss_count": 1, "batch": 1, "reuse": cell["reuse"],
        "iterations": iterations, "warmup": warmup,
        "thread_count": 1, "seed": cell["seed"]}
    parameters = result.get("parameters")
    resolved = result.get("resolved")
    build = result.get("build")
    correctness = result.get("correctness")
    digests = result.get("workload_digests")
    require(result.get("schema") == "leopard2-benchmark-v5" and
            isinstance(parameters, dict) and
            all(parameters.get(key) == value for key, value in expected_parameters.items()),
            f"{label} schema or parameters changed")
    require(isinstance(resolved, dict) and
            resolved.get("profile") == "legacy_high_v1" and
            resolved.get("field") == "gf8" and resolved.get("backend") == "avx2" and
            resolved.get("thread_count") == 1,
            f"{label} resolved identity changed")
    require(isinstance(build, dict) and build.get("source_commit") == expected_commit and
            build.get("source_tree") == expected_tree and
            build.get("source_tracked_dirty") is False,
            f"{label} embedded source identity changed")
    require(isinstance(correctness, dict) and
            correctness.get("leopard2_round_trip") is True,
            f"{label} round trip failed")
    require(isinstance(digests, dict) and digests.get("algorithm") == "fnv1a64" and
            all(isinstance(digests.get(name), str) and len(digests[name]) == 16
                for name in ("original_data", "transmitted_parity", "recovered_originals")),
            f"{label} workload digests are incomplete")
    return {
        "implementation": label, "command": invocation,
        "stdout_sha256": hashlib.sha256(completed.stdout).hexdigest(),
        "stderr_sha256": hashlib.sha256(completed.stderr).hexdigest(),
        "encode_us": BASE.positive_encode_metric(result),
        "digests": dict(digests), "result": result,
    }


def analyze(cell: Mapping[str, Any], rounds: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    ratios = []
    reference = rounds[0]["invocations"][0]["digests"]
    for round_value in rounds:
        invocations = round_value["invocations"]
        require(round_value["isolation"]["accepted"] is True,
                "contaminated round cannot be analyzed")
        require(all(item["digests"] == reference for item in invocations),
                "parent/new workload digests differ")
        new_log = statistics.mean(math.log(item["encode_us"])
                                  for item in invocations
                                  if item["implementation"] == "new")
        parent_log = statistics.mean(math.log(item["encode_us"])
                                     for item in invocations
                                     if item["implementation"] == "parent")
        ratios.append(parent_log - new_log)
    return {
        "cell": dict(cell), "digests": reference,
        "parent_over_new": BASE.confidence_interval(ratios)}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--new", required=True, type=Path)
    parser.add_argument("--parent", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cpu", required=True, type=int)
    parser.add_argument("--sibling", required=True, type=int)
    parser.add_argument("--iterations", type=int, default=15)
    parser.add_argument("--warmup", type=int, default=64)
    return parser.parse_args()


def main() -> int:
    options = parse_arguments()
    require(options.iterations >= 3 and options.warmup >= 1,
            "insufficient repetitions")
    require(not options.output.exists(), "output path already exists")
    raw: dict[str, Any] = {"schema": SCHEMA, "created_utc": MAIN_SUPPORT.utc_now()}
    descriptor = None
    try:
        options.output.mkdir(parents=True)
        descriptor = os.open(LOCK_PATH, os.O_RDWR | os.O_CREAT, 0o600)
        fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        identities = {
            "new": BASE.T8_SUPPORT.file_identity(options.new),
            "parent": BASE.T8_SUPPORT.file_identity(options.parent)}
        require(identities["new"]["sha256"] != identities["parent"]["sha256"],
                "parent and new binaries are not distinct")
        executable_sections = {
            label: BASE.T8_SUPPORT.executable_sections_identity(Path(identity["path"]))
            for label, identity in identities.items()}
        hot_symbols = {
            label: hot_identities(Path(identity["path"]), label == "new")
            for label, identity in identities.items()}
        require(set(os.sched_getaffinity(0)) == {options.cpu},
                "runner must be singleton-pinned")
        siblings = MAIN_SUPPORT.parse_cpu_list(Path(
            f"/sys/devices/system/cpu/cpu{options.cpu}/topology/"
            "thread_siblings_list").read_text(encoding="ascii"))
        require(siblings == {options.cpu, options.sibling},
                "requested CPUs are not one SMT pair")
        raw.update({
            "new_commit": NEW_COMMIT, "new_tree": NEW_TREE,
            "parent_commit": PARENT_COMMIT, "parent_tree": PARENT_TREE,
            "cpu": options.cpu, "reserved_sibling": options.sibling,
            "iterations": options.iterations, "warmup": options.warmup,
            "runner": BASE.T8_SUPPORT.file_identity(RUNNER),
            "runner_dependencies": [BASE.support_file_identity(path) for path in DEPENDENCIES],
            "identities": identities, "executable_sections": executable_sections,
            "hot_symbols": hot_symbols, "cells": []})
        with MAIN_SUPPORT.StableLeaseAnchor(), \
                MAIN_SUPPORT.PairLease(options.cpu, options.sibling) as lease:
            before_cpu = MAIN_SUPPORT.cpu_stat_snapshot(options.cpu)
            before_sibling = MAIN_SUPPORT.cpu_stat_snapshot(options.sibling)
            before_ns = time.monotonic_ns()
            time.sleep(5)
            presample = MAIN_SUPPORT.isolation_record(
                options.cpu, options.sibling, lease, before_ns, time.monotonic_ns(),
                before_cpu, MAIN_SUPPORT.cpu_stat_snapshot(options.cpu),
                before_sibling, MAIN_SUPPORT.cpu_stat_snapshot(options.sibling))
            raw["presample"] = presample
            require(presample["delta"]["benchmark_cpu"]["nonidle_jiffies"] <= 1 and
                    presample["delta"]["reserved_sibling"]["nonidle_jiffies"] == 0,
                    "CPU pair was not quiet during presample")
            for cell_index, cell in enumerate(cells()):
                cell_raw = {"cell": dict(cell), "rounds": []}
                for round_index, order in enumerate(ORDER):
                    before_cpu = MAIN_SUPPORT.cpu_stat_snapshot(options.cpu)
                    before_sibling = MAIN_SUPPORT.cpu_stat_snapshot(options.sibling)
                    before_ns = time.monotonic_ns()
                    invocations = [run_one(
                        label, identities[label], cell, options.cpu,
                        options.iterations, options.warmup, options.output)
                        for label in order]
                    isolation = MAIN_SUPPORT.isolation_record(
                        options.cpu, options.sibling, lease, before_ns,
                        time.monotonic_ns(), before_cpu,
                        MAIN_SUPPORT.cpu_stat_snapshot(options.cpu), before_sibling,
                        MAIN_SUPPORT.cpu_stat_snapshot(options.sibling))
                    require(isolation["accepted"],
                            f"contaminated {cell['id']} round {round_index}")
                    cell_raw["rounds"].append({
                        "round": round_index, "order": list(order),
                        "invocations": invocations, "isolation": isolation})
                raw["cells"].append(cell_raw)
                print(f"{cell_index + 1}/{len(cells())} {cell['id']}",
                      file=sys.stderr, flush=True)
        post_identities = {
            label: BASE.T8_SUPPORT.file_identity(Path(identity["path"]))
            for label, identity in identities.items()}
        post_sections = {
            label: BASE.T8_SUPPORT.executable_sections_identity(Path(identity["path"]))
            for label, identity in identities.items()}
        post_hot = {
            label: hot_identities(Path(identity["path"]), label == "new")
            for label, identity in identities.items()}
        require(post_identities == identities and post_sections == executable_sections and
                post_hot == hot_symbols,
                "binary, executable-section, or hot-symbol identity changed")
        raw["identities_after"] = post_identities
        raw["executable_sections_after"] = post_sections
        raw["hot_symbols_after"] = post_hot
        analyses = [analyze(item["cell"], item["rounds"]) for item in raw["cells"]]
        regressions = [item["cell"]["id"] for item in analyses
                       if item["parent_over_new"]["ci95"][1] < REGRESSION_FLOOR]
        raw["completed_utc"] = MAIN_SUPPORT.utc_now()
        BASE.T8_SUPPORT.write_exclusive(options.output / "raw.json", raw)
        summary = {
            "schema": SUMMARY_SCHEMA,
            "status": "accepted" if not regressions else "rejected",
            "credible_regressions_over_2_percent": regressions,
            "all_digests_matched": True,
            "all_rounds_zero_sibling_nonidle": all(
                round_value["isolation"]["delta"]["reserved_sibling"]
                    ["nonidle_jiffies"] == 0
                for cell_value in raw["cells"]
                for round_value in cell_value["rounds"]),
            "cells": analyses,
            "binary_sha256": {label: value["sha256"]
                              for label, value in identities.items()},
            "executable_sections_sha256": {
                label: value["combined_sha256"]
                for label, value in executable_sections.items()},
            "hot_symbols": hot_symbols,
            "raw_sha256": sha256(options.output / "raw.json")}
        BASE.T8_SUPPORT.write_exclusive(options.output / "summary.json", summary)
        print(json.dumps({
            "status": summary["status"], "cells": len(analyses),
            "credible_regressions_over_2_percent": regressions,
            "ratios": {item["cell"]["id"]: item["parent_over_new"]
                       for item in analyses}}, sort_keys=True))
        return 0 if not regressions else 2
    except Exception as error:
        raw["failed_utc"] = MAIN_SUPPORT.utc_now()
        raw["failure"] = {"type": type(error).__name__, "message": str(error)}
        if options.output.exists():
            BASE.T8_SUPPORT.write_exclusive(options.output / "failure.json", raw)
        print(f"evidence rejected: {error}", file=sys.stderr)
        return 1
    finally:
        if descriptor is not None:
            os.close(descriptor)


if __name__ == "__main__":
    raise SystemExit(main())
