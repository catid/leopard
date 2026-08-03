#!/usr/bin/env python3
"""Authoritative low-P16 partial-output encode ABBA campaign.

One frozen, source-attested Leopard2 executable supplies both comparison arms:
runtime mode 0 retains the mature scratch-and-scatter path, while runtime mode
1 enables direct final-radix output.  Exact Leopard main is measured for the
three comparable target cells as descriptive context only.  The low-v1 and
legacy-high-v1 parity formats differ, so only original/recovered-data digests
are compared with main; the two same-source modes must match all digests.

The runner imports the project's audited executable-freezing, canonical-lock,
SMT-pair lease, descendant-containment, and schedstat-gating support.  Each
benchmark process retains 21 raw timing samples after five warmups.  Independent
observations are mirrored round-level log contrasts: nine for targets and three
for neighboring cells.
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
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA = "leopard2-low-p16-partial-output-abba/v1"
MANIFEST_SCHEMA = "leopard2-low-p16-partial-output-manifest/v1"
CELL_SCHEMA = "leopard2-low-p16-partial-output-cell/v1"
SUMMARY_SCHEMA = "leopard2-low-p16-partial-output-summary/v1"
SEED_SCHEMA = "leopard2-low-p16-partial-output-seed/v1"
SOURCE_COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")
TARGET_ROUNDS = 9
NEIGHBOR_ROUNDS = 3
ITERATIONS = 21
WARMUP = 5
REUSE = 4096
LOSS_COUNT = 1
TIMEOUT_SECONDS = 120.0
TARGET_FLOOR = 1.05
NEIGHBOR_FLOOR = 1.0 / 1.02
DEFAULT_CPU = 13
DEFAULT_SIBLING = 29
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
T95 = {3: 4.302652729911275, 9: 2.306004135204166}


class EvidenceError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvidenceError(message)


def load_support() -> Any:
    path = (Path(__file__).resolve().parent.parent / "direct_repair" /
            "run_equal_rounded_abba.py").resolve(strict=True)
    specification = importlib.util.spec_from_file_location(
        "leopard2_low_p16_equal_rounded_support", path)
    require(specification is not None and specification.loader is not None,
            f"cannot import authoritative evidence support: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


SUPPORT = load_support()
GATE = SUPPORT.GATE
PAIR = SUPPORT.PAIR
T8 = SUPPORT.T8
LOCK_PATH = SUPPORT.LOCK_PATH
MAX_ATTEMPTS = SUPPORT.MAX_ATTEMPTS
MAIN_COMMIT = SUPPORT.MAIN_COMMIT
DIRECTORY = Path(__file__).resolve().parent
SUPPORT_PATH = (DIRECTORY.parent / "direct_repair" /
                "run_equal_rounded_abba.py").resolve(strict=True)
GATE_PATH = (DIRECTORY.parent / "direct_repair" /
             "run_small_direct_abba.py").resolve(strict=True)
MAIN_SUPPORT_PATH = (DIRECTORY.parent / "main_compare" /
                     "run_abba.py").resolve(strict=True)
EVIDENCE_ERROR_TYPES = tuple(dict.fromkeys((
    EvidenceError, SUPPORT.EvidenceError, GATE.EvidenceError,
    PAIR.EvidenceError, T8.EvidenceError,
)))


def canonical_bytes(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      allow_nan=False).encode("utf-8")


def digest_object(value: object) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def stable_seed(identifier: str) -> int:
    payload = f"{SEED_SCHEMA}:{identifier}".encode("ascii")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")


def ceil_power_of_two(value: int) -> int:
    require(type(value) is int and value > 0, "count is not positive")
    return 1 << (value - 1).bit_length()


def route_eligible(k: int, r: int) -> bool:
    """Derive eligibility from code geometry, never from benchmark role."""
    p = ceil_power_of_two(k)
    return p == 16 and (r % p) in (4, 8, 12)


def make_cell(identifier: str, k: int, r: int, shard_bytes: int,
              role: str) -> dict[str, Any]:
    require(role in ("target", "neighbor") and k >= LOSS_COUNT and
            r >= LOSS_COUNT and shard_bytes > 0,
            f"invalid cell: {identifier}")
    p = ceil_power_of_two(k)
    parent = ceil_power_of_two(p + r)
    require(parent <= 256, f"cell is not a GF8 low-profile code: {identifier}")
    if role == "target":
        require(r <= k and shard_bytes % 64 == 0,
                f"target is not directly comparable with main: {identifier}")
    return {
        "id": identifier,
        "K": k,
        "R": r,
        "bytes": shard_bytes,
        "loss": LOSS_COUNT,
        "role": role,
        "P": p,
        "N": parent,
        "partial": r % p,
        "candidate_route_expected": route_eligible(k, r),
        "main_descriptive": role == "target",
        "seed": stable_seed(identifier),
    }


def matrix() -> list[dict[str, Any]]:
    cells = [
        make_cell("target-k16-r12-b64", 16, 12, 64, "target"),
        make_cell("target-k16-r12-b1024", 16, 12, 1024, "target"),
        make_cell("target-k16-r8-b64", 16, 8, 64, "target"),
        make_cell("neighbor-k16-r4-b128", 16, 4, 128, "neighbor"),
        make_cell("neighbor-k16-r4-b2048", 16, 4, 2048, "neighbor"),
        make_cell("neighbor-k9-r20-b65", 9, 20, 65, "neighbor"),
        make_cell("neighbor-k9-r228-b64", 9, 228, 64, "neighbor"),
        make_cell("neighbor-k16-r11-b64", 16, 11, 64, "neighbor"),
        make_cell("neighbor-k16-r13-b64", 16, 13, 64, "neighbor"),
        make_cell("neighbor-k8-r4-b64", 8, 4, 64, "neighbor"),
        make_cell("neighbor-k17-r20-b64", 17, 20, 64, "neighbor"),
        make_cell("neighbor-k9-r16-b256", 9, 16, 256, "neighbor"),
    ]
    require(len({item["id"] for item in cells}) == len(cells),
            "matrix contains duplicate identifiers")
    return cells


def load_json(path: Path) -> dict[str, Any]:
    return SUPPORT.load_json(path)


def write_json(path: Path, value: object) -> None:
    SUPPORT.write_atomic_exclusive(path, value)


def artifact_is_current(identity: Mapping[str, Any]) -> bool:
    path = Path(str(identity.get("path", "")))
    return (path.is_file() and (path.stat().st_mode & 0o777) == 0o555 and
            path.stat().st_size == identity.get("size") and
            T8.sha256(path) == identity.get("sha256"))


def validate_source_ids(commit: str, tree: str) -> None:
    require(SOURCE_COMMIT_RE.fullmatch(commit) is not None,
            "--source-commit must be a lowercase full Git object ID")
    require(SOURCE_COMMIT_RE.fullmatch(tree) is not None,
            "--source-tree must be a lowercase full Git object ID")


def expected_orders(item: Mapping[str, Any]) -> tuple[tuple[str, ...], ...]:
    if item["role"] == "target":
        return TARGET_ORDER * (TARGET_ROUNDS // len(TARGET_ORDER))
    return NEIGHBOR_ORDER


def expected_parameters(item: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "K": item["K"],
        "R": item["R"],
        "shard_bytes": item["bytes"],
        "loss_count": item["loss"],
        "batch": 1,
        "reuse": REUSE,
        "iterations": ITERATIONS,
        "warmup": WARMUP,
        "thread_count": 1,
        "seed": item["seed"],
    }


def benchmark_command(implementation: str, executable: Path,
                      item: Mapping[str, Any]) -> list[str]:
    require(implementation in ("main", "control", "candidate"),
            "unknown implementation")
    command = [
        "/usr/bin/prlimit", "--as=268435456", str(executable),
        "--k", str(item["K"]), "--r", str(item["R"]),
        "--bytes", str(item["bytes"]), "--loss", str(item["loss"]),
        "--batch", "1", "--reuse", str(REUSE),
        "--iterations", str(ITERATIONS), "--warmup", str(WARMUP),
        "--threads", "1", "--seed", str(item["seed"]),
    ]
    if implementation != "main":
        mode = 1 if implementation == "candidate" else 0
        command.extend((
            "--profile", "low", "--field", "gf8", "--backend", "avx2",
            "--skip-legacy", "--retain-samples", "--measure-one-shot-encode",
            "--attest-source", "--low-p16-partial-direct-output-mode",
            str(mode),
        ))
    command.extend(("--json", "-"))
    return command


def validate_common_parameters(value: Mapping[str, Any],
                               item: Mapping[str, Any]) -> Mapping[str, Any]:
    parameters = value.get("parameters")
    require(isinstance(parameters, dict), "benchmark parameters are absent")
    for name, expected in expected_parameters(item).items():
        require(parameters.get(name) == expected,
                f"benchmark parameter changed: {name}")
    missing = parameters.get("missing_original_indices")
    require(isinstance(missing, list) and len(missing) == item["loss"] and
            all(type(index) is int and 0 <= index < item["K"] for index in missing)
            and len(set(missing)) == len(missing) and missing == sorted(missing),
            "missing-original pattern is invalid")
    return parameters


def validate_result(implementation: str, value: object,
                    item: Mapping[str, Any], source_commit: str,
                    source_tree: str) -> dict[str, Any]:
    require(isinstance(value, dict), "benchmark output is not an object")
    parameters = validate_common_parameters(value, item)
    resolved = value.get("resolved")
    correctness = value.get("correctness")
    build = value.get("build")
    metrics = value.get("metrics")
    require(all(isinstance(part, dict) for part in
                (resolved, correctness, build, metrics)),
            "benchmark output is incomplete")
    digests = PAIR.validate_digest_object(value.get("workload_digests"))

    if implementation == "main":
        require(item["main_descriptive"] is True,
                "exact main was scheduled for an incomparable neighbor")
        t = ceil_power_of_two(int(item["R"]))
        n = ceil_power_of_two(int(item["K"]) + t)
        require(
            value.get("schema") == "leopard-main-benchmark-v1" and
            build.get("main_source_commit") == MAIN_COMMIT and
            correctness.get("round_trip") is True and
            resolved.get("profile") == "legacy_high_v1" and
            resolved.get("field") == "gf8" and
            resolved.get("thread_count") == 1 and
            resolved.get("parent_count") == n and
            resolved.get("padded_side") == t,
            "exact-main identity, geometry, or round trip changed")
        samples = PAIR.validate_summary(
            metrics.get("encode_execution"), ITERATIONS)
        return {
            "encode_execution_us": statistics.median(samples),
            "digests": digests,
            "missing_original_indices": parameters["missing_original_indices"],
            "route_expected": None,
            "route_selected": None,
        }

    mode = 1 if implementation == "candidate" else 0
    expected_selected = mode == 1 and route_eligible(
        int(item["K"]), int(item["R"]))
    require(
        value.get("schema") == "leopard2-benchmark-v21" and
        build.get("source_commit") == source_commit and
        build.get("source_tree") == source_tree and
        build.get("source_tracked_dirty") is False and
        build.get("low_p16_partial_direct_output_diagnostic_mode") == mode and
        build.get("low_p16_partial_direct_output_mode_word") ==
            (1 if mode == 1 else 2) and
        build.get("low_p16_partial_direct_output_enabled") is (mode == 1) and
        build.get("low_p16_partial_direct_output_route_expected_selected")
            is expected_selected and
        build.get("low_p16_partial_direct_output_route_selected")
            is expected_selected and
        build.get("low_p16_partial_direct_output_selector_contract") ==
            "LOW_V1,GF8,AVX2,P=16,d=R%16 in {4,8,12},"
            "dense transmitted prefix,transform encode",
        f"{implementation} v21 source or route attestation changed")
    require(
        parameters.get("requested_profile") == "low_v1" and
        parameters.get("requested_field") == "gf8" and
        parameters.get("requested_backend") == "avx2" and
        parameters.get("skip_legacy") is True and
        parameters.get("retain_samples") is True and
        parameters.get("measure_one_shot_encode") is True and
        parameters.get("attest_source") is True and
        parameters.get("low_p16_partial_direct_output_mode") == mode and
        resolved.get("profile") == "low_v1" and
        resolved.get("field") == "gf8" and
        resolved.get("backend") == "avx2" and
        resolved.get("thread_count") == 1 and
        resolved.get("parent_count") == item["N"] and
        resolved.get("padded_side") == item["P"] and
        correctness.get("leopard2_round_trip") is True and
        correctness.get("legacy_comparison") is None,
        f"{implementation} requested/resolved identity changed")
    legacy = value.get("legacy")
    require(isinstance(legacy, dict) and legacy.get("available") is False and
            legacy.get("unavailable_reason") == "disabled by --skip-legacy",
            f"{implementation} unexpectedly ran in-tree Leopard1")
    execution_samples = PAIR.validate_summary(
        metrics.get("encode_execution"), ITERATIONS)
    one_shot_samples = PAIR.validate_summary(
        metrics.get("one_shot_encode"), ITERATIONS)
    return {
        "encode_execution_us": statistics.median(execution_samples),
        "one_shot_encode_us": statistics.median(one_shot_samples),
        "digests": digests,
        "missing_original_indices": parameters["missing_original_indices"],
        "route_expected": expected_selected,
        "route_selected": expected_selected,
    }


def run_one(implementation: str, executable: Path,
            identity: Mapping[str, Any], item: Mapping[str, Any],
            source_commit: str, source_tree: str, cpu: int, sibling: int,
            raw_directory: Path, round_index: int, slot_index: int,
            lock_descriptor: int) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    require(artifact_is_current(identity),
            f"{implementation} frozen executable changed before invocation")
    command = benchmark_command(implementation, executable, item)
    rejected: list[dict[str, Any]] = []
    for attempt in range(MAX_ATTEMPTS):
        stem = (f"round-{round_index:02d}-slot-{slot_index}-"
                f"{implementation}-attempt-{attempt:02d}")
        stdout_path = raw_directory / f"{stem}.stdout.json"
        stderr_path = raw_directory / f"{stem}.stderr.txt"
        gated = GATE.run_gated_benchmark(
            command, cpu, sibling, TIMEOUT_SECONDS, stdout_path, stderr_path,
            campaign_lock_descriptor=lock_descriptor)
        require(artifact_is_current(identity),
                f"{implementation} frozen executable changed after invocation")
        record: dict[str, Any] = {
            "implementation": implementation,
            "attempt": attempt,
            "command": command,
            "executable_sha256": identity["sha256"],
            "stdout": SUPPORT.retained_file_identity(
                stdout_path, 8 << 20, False),
            "stderr": SUPPORT.retained_file_identity(
                stderr_path, 1 << 20, False),
            "gated": gated,
            "accepted": False,
        }
        if gated["return_code"] != 0:
            envelope = raw_directory / f"{stem}.envelope.json"
            write_json(envelope, record)
            record["envelope"] = SUPPORT.retained_file_identity(
                envelope, 32 << 20, True)
            raise EvidenceError(
                f"{item['id']} {implementation} benchmark failed; "
                f"see {stderr_path}")
        if not SUPPORT.isolation_accepted(gated):
            envelope = raw_directory / f"{stem}.envelope.json"
            write_json(envelope, record)
            record["envelope"] = SUPPORT.retained_file_identity(
                envelope, 32 << 20, True)
            rejected.append(record)
            continue
        result = load_json(stdout_path)
        normalized = validate_result(
            implementation, result, item, source_commit, source_tree)
        record["result"] = result
        record["normalized"] = normalized
        record["accepted"] = True
        envelope = raw_directory / f"{stem}.envelope.json"
        write_json(envelope, record)
        record["envelope"] = SUPPORT.retained_file_identity(
            envelope, 32 << 20, True)
        return record, rejected
    raise EvidenceError(
        f"{item['id']} {implementation} failed objective isolation "
        f"for {MAX_ATTEMPTS} attempts")


def round_log_ratio(invocations: Sequence[Mapping[str, Any]], numerator: str,
                    denominator: str, numerator_metric: str,
                    denominator_metric: str | None = None) -> float:
    denominator_metric = denominator_metric or numerator_metric
    numerator_values = [
        float(value["normalized"][numerator_metric])
        for value in invocations if value["implementation"] == numerator
    ]
    denominator_values = [
        float(value["normalized"][denominator_metric])
        for value in invocations if value["implementation"] == denominator
    ]
    require(len(numerator_values) == len(denominator_values) == 2 and
            all(value > 0 for value in (*numerator_values, *denominator_values)),
            "round is not a balanced positive timing comparison")
    return (statistics.mean(math.log(value) for value in numerator_values) -
            statistics.mean(math.log(value) for value in denominator_values))


def confidence(log_ratios: Sequence[float]) -> dict[str, Any]:
    require(len(log_ratios) in T95, "unsupported independent round count")
    center = statistics.mean(log_ratios)
    half = (T95[len(log_ratios)] * statistics.stdev(log_ratios) /
            math.sqrt(len(log_ratios)))
    return {
        "orientation": "numerator_time_over_denominator_time",
        "independent_round_count": len(log_ratios),
        "degrees_of_freedom": len(log_ratios) - 1,
        "round_log_contrasts": list(log_ratios),
        "ratio": math.exp(center),
        "lower": math.exp(center - half),
        "upper": math.exp(center + half),
    }


def digest_checks(item: Mapping[str, Any],
                  rounds: Sequence[Mapping[str, Any]]) -> dict[str, bool]:
    mode_records = [
        invocation for round_value in rounds
        for invocation in round_value["invocations"]
        if invocation["implementation"] in ("control", "candidate")
    ]
    mode_digests = [record["normalized"]["digests"] for record in mode_records]
    mode_patterns = [
        record["normalized"]["missing_original_indices"]
        for record in mode_records
    ]
    checks = {
        "mode_original_equal": len({value["original_data"]
                                     for value in mode_digests}) == 1,
        "mode_parity_equal": len({value["transmitted_parity"]
                                   for value in mode_digests}) == 1,
        "mode_recovered_equal": len({value["recovered_originals"]
                                      for value in mode_digests}) == 1,
        "mode_loss_pattern_equal": len({tuple(value)
                                         for value in mode_patterns}) == 1,
        "main_original_equal": True,
        "main_recovered_equal": True,
        "main_loss_pattern_equal": True,
    }
    if item["main_descriptive"]:
        main_records = [
            invocation for round_value in rounds
            for invocation in round_value["invocations"]
            if invocation["implementation"] == "main"
        ]
        require(main_records, "target has no exact-main records")
        reference = mode_digests[0]
        reference_pattern = tuple(mode_patterns[0])
        checks["main_original_equal"] = all(
            record["normalized"]["digests"]["original_data"] ==
            reference["original_data"] for record in main_records)
        checks["main_recovered_equal"] = all(
            record["normalized"]["digests"]["recovered_originals"] ==
            reference["recovered_originals"] for record in main_records)
        checks["main_loss_pattern_equal"] = all(
            tuple(record["normalized"]["missing_original_indices"]) ==
            reference_pattern for record in main_records)
    return checks


def analyze_cell(item: Mapping[str, Any],
                 rounds: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    control = {}
    for metric in ("encode_execution_us", "one_shot_encode_us"):
        control[metric] = confidence([
            round_log_ratio(round_value["invocations"], "control",
                            "candidate", metric)
            for round_value in rounds
        ])
    result: dict[str, Any] = {
        "cell": dict(item),
        "control_over_candidate": control,
        "digest_checks": digest_checks(item, rounds),
    }
    if item["main_descriptive"]:
        result["main_over_candidate_descriptive"] = {
            metric: confidence([
                round_log_ratio(round_value["invocations"], "main",
                                "candidate", "encode_execution_us", metric)
                for round_value in rounds
            ]) for metric in ("encode_execution_us", "one_shot_encode_us")
        }
    return result


def aggregate(analyses: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    require(analyses, "campaign has no analyses")
    target_failures: list[dict[str, str]] = []
    neighbor_failures: list[dict[str, str]] = []
    digest_failures: list[str] = []
    for analysis in analyses:
        identifier = analysis["cell"]["id"]
        if not all(analysis["digest_checks"].values()):
            digest_failures.append(identifier)
        for metric, interval in analysis["control_over_candidate"].items():
            if analysis["cell"]["role"] == "target":
                if interval["lower"] < TARGET_FLOOR:
                    target_failures.append({"cell": identifier,
                                            "metric": metric})
            elif interval["upper"] < NEIGHBOR_FLOOR:
                neighbor_failures.append({"cell": identifier,
                                          "metric": metric})
    return {
        "target_count": sum(value["cell"]["role"] == "target"
                            for value in analyses),
        "neighbor_count": sum(value["cell"]["role"] == "neighbor"
                              for value in analyses),
        "target_floor": TARGET_FLOOR,
        "neighbor_credible_regression_floor": NEIGHBOR_FLOOR,
        "target_failures": target_failures,
        "neighbor_failures": neighbor_failures,
        "digest_failures": digest_failures,
        "accepted": not target_failures and not neighbor_failures and
                    not digest_failures,
        "main_is_descriptive_only": True,
    }


def validate_retained_identity(identity: Mapping[str, Any], maximum: int,
                               nonempty: bool) -> None:
    path = Path(str(identity.get("path", ""))).resolve(strict=True)
    status = path.stat()
    require(path.is_file() and status.st_size == identity.get("size") and
            status.st_size <= maximum and (status.st_size > 0 or not nonempty)
            and (status.st_mode & 0o777) == identity.get("mode") and
            T8.sha256(path) == identity.get("sha256"),
            f"retained artifact changed: {path}")


def validate_cell_record(value: object, expected: Mapping[str, Any],
                         manifest_sha256: str, root: Path,
                         identities: Mapping[str, Any], source_commit: str,
                         source_tree: str) -> dict[str, Any]:
    require(isinstance(value, dict) and value.get("schema") == CELL_SCHEMA and
            value.get("manifest_sha256") == manifest_sha256 and
            value.get("cell") == expected,
            "cell record identity changed")
    rounds = value.get("rounds")
    orders = expected_orders(expected)
    require(isinstance(rounds, list) and len(rounds) == len(orders),
            "cell round count changed")
    for round_index, (round_value, order) in enumerate(zip(rounds, orders)):
        require(isinstance(round_value, dict) and
                round_value.get("round") == round_index and
                round_value.get("order") == list(order) and
                isinstance(round_value.get("invocations"), list) and
                len(round_value["invocations"]) == len(order),
                "cell round order changed")
        for invocation, implementation in zip(
                round_value["invocations"], order):
            identity_name = "main" if implementation == "main" else "leopard2"
            identity = identities[identity_name]
            require(
                invocation.get("implementation") == implementation and
                invocation.get("accepted") is True and
                invocation.get("executable_sha256") == identity["sha256"] and
                invocation.get("command") == benchmark_command(
                    implementation, Path(identity["path"]), expected) and
                invocation.get("gated", {}).get("return_code") == 0 and
                SUPPORT.isolation_accepted(invocation["gated"]),
                "cell contains an invalid accepted invocation")
            validate_retained_identity(invocation["stdout"], 8 << 20, True)
            validate_retained_identity(invocation["stderr"], 1 << 20, False)
            validate_retained_identity(invocation["envelope"], 32 << 20, True)
            envelope = load_json(Path(invocation["envelope"]["path"]))
            require(envelope == {
                name: retained for name, retained in invocation.items()
                if name != "envelope"
            }, "accepted invocation envelope changed")
            retained_result = load_json(Path(invocation["stdout"]["path"]))
            require(retained_result == invocation.get("result"),
                    "retained stdout differs from cell result")
            normalized = validate_result(
                implementation, retained_result, expected,
                source_commit, source_tree)
            require(normalized == invocation.get("normalized"),
                    "retained normalized result changed")
    rejected = value.get("rejected_attempts")
    require(isinstance(rejected, list) and all(
        record.get("accepted") is False and
        record.get("gated", {}).get("return_code") == 0 and
        not SUPPORT.isolation_accepted(record["gated"])
        for record in rejected), "rejected-attempt set is invalid")
    for record in rejected:
        validate_retained_identity(record["stdout"], 8 << 20, False)
        validate_retained_identity(record["stderr"], 1 << 20, False)
        validate_retained_identity(record["envelope"], 32 << 20, True)
        envelope = load_json(Path(record["envelope"]["path"]))
        require(envelope == {
            name: retained for name, retained in record.items()
            if name != "envelope"
        }, "rejected invocation envelope changed")
    raw_run = value.get("raw_run_directory")
    prior = value.get("prior_incomplete_raw_directories")
    prefix = f"raw/{expected['id']}/run-"
    require(isinstance(raw_run, str) and raw_run.startswith(prefix) and
            isinstance(prior, list) and prior == sorted(prior) and
            len(prior) == len(set(prior)) and
            all(isinstance(path, str) and path.startswith(prefix)
                for path in prior) and raw_run not in prior and
            (root / raw_run).resolve(strict=True).is_dir(),
            "cell raw-run provenance changed")
    analysis = analyze_cell(expected, rounds)
    require(value.get("analysis") == analysis,
            "cell analysis is not derived from retained evidence")
    return value


def support_identities() -> dict[str, Any]:
    return {
        "runner": T8.regular_file_identity(Path(__file__)),
        "equal_rounded_support": T8.regular_file_identity(SUPPORT_PATH),
        "gated_isolation_support": T8.regular_file_identity(GATE_PATH),
        "main_adapter_support": T8.regular_file_identity(MAIN_SUPPORT_PATH),
    }


def create_identities(output: Path, leopard2: Path,
                      main: Path) -> dict[str, Any]:
    artifacts = output / "artifacts"
    artifacts.mkdir(mode=0o700)
    frozen = {
        "leopard2": SUPPORT.freeze_executable(
            leopard2, artifacts, "bench_leopard2"),
        "main": SUPPORT.freeze_executable(
            main, artifacts, "leopard_main_benchmark"),
    }
    for record in frozen.values():
        require(T8.sha256(Path(record["origin"]["path"])) ==
                record["origin"]["sha256"],
                "source executable changed while being frozen")
    executables = {name: value["frozen"] for name, value in frozen.items()}
    require(executables["leopard2"]["sha256"] !=
            executables["main"]["sha256"],
            "Leopard2 and exact-main artifacts are unexpectedly identical")
    identities = {"frozen": frozen, "executables": executables}
    write_json(output / "identities.json", identities)
    return identities


def run(options: argparse.Namespace) -> int:
    validate_source_ids(options.source_commit, options.source_tree)
    allowed = set(os.sched_getaffinity(0))
    require(options.cpu not in allowed and options.sibling not in allowed,
            "coordinator affinity must exclude the measured CPU pair")
    sibling_text = Path(
        f"/sys/devices/system/cpu/cpu{options.cpu}/topology/"
        "thread_siblings_list").read_text(encoding="ascii")
    require(PAIR.parse_cpu_list(sibling_text) ==
            {options.cpu, options.sibling},
            "requested CPUs are not one SMT pair")
    cells = matrix()
    if options.cell_id:
        by_id = {item["id"]: item for item in cells}
        unknown = sorted(set(options.cell_id) - set(by_id))
        require(not unknown, f"unknown cell IDs: {unknown}")
        cells = [by_id[identifier] for identifier in options.cell_id]
    require(cells, "campaign has no cells")
    require(options.output.is_absolute(),
            "output must be an absolute lane-owned directory")
    if options.resume:
        require(options.output.is_dir(), "resume output directory is absent")
    else:
        require(options.leopard2 is not None and options.main is not None,
                "--leopard2 and --main are required for a new campaign")
        require(not options.output.exists(), "output directory already exists")
        options.output.mkdir(mode=0o700, parents=True)

    lock_descriptor = os.open(
        LOCK_PATH, os.O_RDWR | os.O_CREAT | os.O_CLOEXEC, 0o600)
    fcntl.flock(lock_descriptor, fcntl.LOCK_EX)
    try:
        if options.resume:
            identities = load_json(options.output / "identities.json")
        else:
            identities = create_identities(
                options.output, options.leopard2, options.main)
        executables = identities.get("executables")
        require(isinstance(executables, dict) and
                set(executables) == {"leopard2", "main"} and
                all(artifact_is_current(value)
                    for value in executables.values()),
                "frozen executable identities are invalid")

        manifest = {
            "schema": MANIFEST_SCHEMA,
            "runner_schema": SCHEMA,
            "source_commit": options.source_commit,
            "source_tree": options.source_tree,
            "source_clean_required": True,
            "main_commit": MAIN_COMMIT,
            "target_rounds": TARGET_ROUNDS,
            "neighbor_rounds": NEIGHBOR_ROUNDS,
            "iterations": ITERATIONS,
            "warmup": WARMUP,
            "reuse": REUSE,
            "loss_count": LOSS_COUNT,
            "cpu": options.cpu,
            "sibling": options.sibling,
            "coordinator_affinity": sorted(allowed),
            "support": support_identities(),
            "identities": identities,
            "promotion_policy": {
                "target": "both control/candidate CI95 lower >= 1.05",
                "neighbor":
                    "fail only when CI95 upper < 1/1.02 (credible >2% regression)",
                "main": "descriptive only; never a promotion gate",
            },
            "isolation_policy": {
                "canonical_lock": str(LOCK_PATH),
                "maximum_attempts": MAX_ATTEMPTS,
                "retry_trigger": "objective isolation failure only",
                "target_endpoint_absolute_ns":
                    SUPPORT.TARGET_ENDPOINT_ABSOLUTE_TOLERANCE_NS,
                "target_endpoint_relative_ppm":
                    SUPPORT.TARGET_ENDPOINT_RELATIVE_TOLERANCE_PPM,
                "sibling_absolute_ns":
                    SUPPORT.SIBLING_ABSOLUTE_TOLERANCE_NS,
                "sibling_relative_ppm":
                    SUPPORT.SIBLING_RELATIVE_TOLERANCE_PPM,
            },
            "cells": cells,
        }
        manifest_sha256 = digest_object(manifest)
        manifest["digest"] = manifest_sha256
        manifest_path = options.output / "manifest.json"
        if options.resume:
            require(load_json(manifest_path) == manifest,
                    "resume manifest changed")
        else:
            write_json(manifest_path, manifest)
            (options.output / "cells").mkdir(mode=0o700)
            (options.output / "raw").mkdir(mode=0o700)
        cells_directory = options.output / "cells"
        raw_root = options.output / "raw"
        require(cells_directory.is_dir() and raw_root.is_dir(),
                "cell/raw directories are absent")

        analyses: list[dict[str, Any]] = []
        rejected_count = 0
        with PAIR.StableLeaseAnchor(), \
                PAIR.PairLease(options.cpu, options.sibling) as pair_lease:
            for cell_index, item in enumerate(cells):
                cell_path = cells_directory / f"{item['id']}.json"
                if cell_path.exists():
                    require(options.resume,
                            f"cell already exists: {item['id']}")
                    record = validate_cell_record(
                        load_json(cell_path), item, manifest_sha256,
                        options.output, executables,
                        options.source_commit, options.source_tree)
                    analyses.append(record["analysis"])
                    rejected_count += len(record["rejected_attempts"])
                    print(f"{cell_index + 1}/{len(cells)} {item['id']} resumed",
                          file=sys.stderr, flush=True)
                    continue

                raw_directory, prior = SUPPORT.allocate_raw_run_directory(
                    raw_root, item["id"])
                rounds: list[dict[str, Any]] = []
                rejected: list[dict[str, Any]] = []
                for round_index, order in enumerate(expected_orders(item)):
                    invocations = []
                    for slot_index, implementation in enumerate(order):
                        identity_name = ("main" if implementation == "main"
                                         else "leopard2")
                        invocation, failed_isolation = run_one(
                            implementation,
                            Path(executables[identity_name]["path"]),
                            executables[identity_name], item,
                            options.source_commit, options.source_tree,
                            options.cpu, options.sibling, raw_directory,
                            round_index, slot_index, lock_descriptor)
                        invocations.append(invocation)
                        rejected.extend(failed_isolation)
                    rounds.append({"round": round_index,
                                   "order": list(order),
                                   "invocations": invocations})
                analysis = analyze_cell(item, rounds)
                value = {
                    "schema": CELL_SCHEMA,
                    "manifest_sha256": manifest_sha256,
                    "pair_lease": pair_lease,
                    "cell": item,
                    "raw_run_directory": str(
                        raw_directory.relative_to(options.output)),
                    "prior_incomplete_raw_directories": prior,
                    "rounds": rounds,
                    "rejected_attempts": rejected,
                    "analysis": analysis,
                }
                write_json(cell_path, value)
                validated = validate_cell_record(
                    load_json(cell_path), item, manifest_sha256,
                    options.output, executables,
                    options.source_commit, options.source_tree)
                analyses.append(validated["analysis"])
                rejected_count += len(rejected)
                print(f"{cell_index + 1}/{len(cells)} {item['id']}",
                      file=sys.stderr, flush=True)

        require(all(artifact_is_current(value)
                    for value in executables.values()),
                "frozen executable changed after campaign")
        summary = {
            "schema": SUMMARY_SCHEMA,
            "manifest_sha256": manifest_sha256,
            "cell_count": len(cells),
            "accepted_process_count": sum(
                len(round_value["invocations"])
                for item in cells
                for round_value in load_json(
                    cells_directory / f"{item['id']}.json")["rounds"]),
            "rejected_isolation_attempt_count": rejected_count,
            "analysis": aggregate(analyses),
        }
        summary["digest"] = digest_object(summary)
        summary_path = options.output / "summary.json"
        if summary_path.exists():
            require(options.resume and load_json(summary_path) == summary,
                    "existing summary differs from replayed campaign")
        else:
            write_json(summary_path, summary)
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 0 if summary["analysis"]["accepted"] else 2
    finally:
        os.close(lock_descriptor)


def verify(options: argparse.Namespace) -> int:
    root = options.output.resolve(strict=True)
    manifest = load_json(root / "manifest.json")
    require(manifest.get("schema") == MANIFEST_SCHEMA and
            manifest.get("runner_schema") == SCHEMA,
            "manifest schema changed")
    digest = manifest.pop("digest", None)
    require(digest == digest_object(manifest), "manifest digest changed")
    manifest["digest"] = digest
    require(manifest.get("support") == support_identities(),
            "evidence support files changed")
    identities_file = load_json(root / "identities.json")
    require(identities_file == manifest.get("identities"),
            "executable identity record changed")
    executables = identities_file.get("executables")
    require(isinstance(executables, dict) and
            all(artifact_is_current(value) for value in executables.values()),
            "frozen executable changed")
    cells = manifest.get("cells")
    require(isinstance(cells, list) and cells, "manifest cells are absent")
    analyses = []
    rejected_count = 0
    accepted_count = 0
    for item in cells:
        record = validate_cell_record(
            load_json(root / "cells" / f"{item['id']}.json"),
            item, digest, root, executables,
            manifest["source_commit"], manifest["source_tree"])
        analyses.append(record["analysis"])
        rejected_count += len(record["rejected_attempts"])
        accepted_count += sum(len(value["invocations"])
                              for value in record["rounds"])
    summary = load_json(root / "summary.json")
    summary_digest = summary.pop("digest", None)
    require(summary.get("schema") == SUMMARY_SCHEMA and
            summary_digest == digest_object(summary),
            "summary identity changed")
    expected_summary = {
        "schema": SUMMARY_SCHEMA,
        "manifest_sha256": digest,
        "cell_count": len(cells),
        "accepted_process_count": accepted_count,
        "rejected_isolation_attempt_count": rejected_count,
        "analysis": aggregate(analyses),
    }
    require(summary == expected_summary,
            "summary is not derived from retained cell evidence")
    print("low-P16 partial-output evidence verified")
    return 0


def self_test() -> int:
    cells = matrix()
    require(len(cells) == 12 and
            sum(item["role"] == "target" for item in cells) == 3 and
            sum(item["role"] == "neighbor" for item in cells) == 9,
            "matrix shape changed")
    require(digest_object(cells) ==
            "d030255fbf07484d516197341e25cc1f122d52ce73edc3cd16678f6ea5b96710",
            "matrix digest changed")
    require(route_eligible(16, 4) and route_eligible(16, 8) and
            route_eligible(16, 12) and route_eligible(9, 20) and
            route_eligible(9, 228) and
            not route_eligible(16, 11) and
            not route_eligible(16, 13) and
            not route_eligible(8, 4) and
            not route_eligible(17, 20) and
            not route_eligible(9, 16),
            "geometry-derived route eligibility changed")
    interval = confidence((math.log(1.1),) * 3)
    require(all(abs(interval[name] - 1.1) < 1e-12
                for name in ("ratio", "lower", "upper")),
            "Student-t log confidence changed")
    fake_target = {
        "cell": {"id": "target", "role": "target"},
        "control_over_candidate": {
            "encode_execution_us": {"lower": 1.051, "upper": 1.06},
            "one_shot_encode_us": {"lower": 1.052, "upper": 1.07},
        },
        "digest_checks": {"all": True},
    }
    fake_neighbor = {
        "cell": {"id": "neighbor", "role": "neighbor"},
        "control_over_candidate": {
            "encode_execution_us": {"lower": 0.96, "upper": 0.99},
            "one_shot_encode_us": {"lower": 0.95, "upper": 0.981},
        },
        "digest_checks": {"all": True},
    }
    require(aggregate((fake_target, fake_neighbor))["accepted"] is True,
            "promotion gates changed")
    fake_neighbor["control_over_candidate"]["one_shot_encode_us"]["upper"] = 0.97
    require(aggregate((fake_target, fake_neighbor))["accepted"] is False,
            "credible neighbor-regression gate changed")
    with tempfile.TemporaryDirectory(
            prefix="leopard2-low-p16-runner-test-") as temporary:
        raw = Path(temporary) / "raw"
        raw.mkdir(mode=0o700)
        first, first_prior = SUPPORT.allocate_raw_run_directory(raw, "cell")
        second, second_prior = SUPPORT.allocate_raw_run_directory(raw, "cell")
        require(first_prior == [] and first != second and
                second_prior == [str(first.relative_to(raw.parent))],
                "resumable raw-run allocation changed")
    print("low-P16 partial-output runner self-test passed")
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    run_parser = commands.add_parser("run")
    run_parser.add_argument("--leopard2", type=Path)
    run_parser.add_argument("--main", type=Path)
    run_parser.add_argument("--source-commit", required=True)
    run_parser.add_argument("--source-tree", required=True)
    run_parser.add_argument("--output", required=True, type=Path)
    run_parser.add_argument("--cpu", type=int, default=DEFAULT_CPU)
    run_parser.add_argument("--sibling", type=int, default=DEFAULT_SIBLING)
    run_parser.add_argument("--cell-id", action="append", default=[])
    run_parser.add_argument("--resume", action="store_true")
    run_parser.set_defaults(function=run)
    verify_parser = commands.add_parser("verify")
    verify_parser.add_argument("--output", required=True, type=Path)
    verify_parser.set_defaults(function=verify)
    test_parser = commands.add_parser("self-test")
    test_parser.set_defaults(function=lambda unused: self_test())
    return result


def main() -> int:
    options = parser().parse_args()
    return int(options.function(options))


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except EVIDENCE_ERROR_TYPES as error:
        print(f"evidence error: {error}", file=sys.stderr)
        raise SystemExit(1)
