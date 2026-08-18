#!/usr/bin/env python3
"""Authoritative 97-cell current-versus-main Leopard2 regression screen.

This is the final one-item-batch, one-shot boundary screen for 69 historical
current-atlas loss targets and 28 controls.  Candidate and control execute the
same default production route through one frozen inode: their contrast measures
replication noise, while exact Leopard main is the performance comparator.
Every process measures both encoding and setup-included decoding.  The
historically losing operation is normalized as the primary metric.
"""

from __future__ import annotations

import contextlib
import fcntl
import hashlib
import importlib.util
import io
import json
import math
import os
import re
import statistics
import sys
import time
from pathlib import Path
from typing import Any, Mapping, Sequence


PARENT_PATH = Path(__file__).resolve().with_name(
    "run_k66r16_b64_tail_abba.py")
SOURCE_ROOT = Path(__file__).resolve().parents[3]
ATLAS_MANIFEST_PATH = \
    SOURCE_ROOT / "docs/performance/leopard2_atlas/manifest.json"
ATLAS_SUMMARY_PATH = \
    SOURCE_ROOT / "docs/performance/leopard2_atlas/summary.json"
ATLAS_GENERATOR_PATH = \
    SOURCE_ROOT / "experiments/leopard2/performance_atlas/generate_atlas.py"
ATLAS_MANIFEST_SHA256 = \
    "dfd8db284b0b70d6dfdfad51a9105d69e568d45d26caea9a31f06a5e83783f55"
ATLAS_SUMMARY_SHA256 = \
    "0dc92139f524c46e8d3b317e60ca47b3529c82305bc93ad1f6a18e7f5c4e19f2"


def load_parent() -> Any:
    specification = importlib.util.spec_from_file_location(
        "current_atlas_regression_screen_parent", PARENT_PATH)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load qualification support: {PARENT_PATH}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


PARENT = load_parent()
BASE = PARENT.BASE


def load_atlas_generator() -> Any:
    specification = importlib.util.spec_from_file_location(
        "current_atlas_regression_screen_generator", ATLAS_GENERATOR_PATH)
    if specification is None or specification.loader is None:
        raise RuntimeError(
            f"cannot load atlas generator: {ATLAS_GENERATOR_PATH}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


ATLAS_GENERATOR = load_atlas_generator()
SUPPORT_DEPENDENCIES = tuple(BASE.RUNNER_DEPENDENCIES)
BASE.__doc__ = __doc__
PARENT.__doc__ = __doc__
PARENT.SUPPORT.__doc__ = __doc__
BASE.SCHEMA = "leopard2-current-atlas-regression-screen-abba/v2"
BASE.SUMMARY_SCHEMA = \
    "leopard2-current-atlas-regression-screen-preliminary-summary/v2"
FINAL_SUMMARY_SCHEMA = \
    "leopard2-current-atlas-regression-screen-summary/v2"
BASE.ALLOW_IDENTICAL_CANDIDATE_CONTROL = True
BASE.ALLOW_MULTIPLE_TARGETS = True
BASE.CANDIDATE_SCHEMA = "leopard2-benchmark-v9"
BASE.CONTROL_SCHEMA = "leopard2-benchmark-v9"
BASE.CONTROL_EXTRA_ARGUMENTS = ()
BASE.CONTROL_BUILD_MARKER = None
BASE.MODE_SYMBOL = None
BASE.TARGET_CONTROL_FLOOR = 1.0 / 1.02
BASE.TARGET_MAIN_FLOOR = 1.0
BASE.NEIGHBOR_FLOOR = 1.0 / 1.02
BASE.NEIGHBOR_EQUIVALENCE_LOWER = 1.0 / 1.02
BASE.NEIGHBOR_EQUIVALENCE_UPPER = 1.02
BASE.MAX_ISOLATION_ATTEMPTS = 8
BASE.RUNNER_PATH = Path(__file__).resolve()
BASE.RUNNER_DEPENDENCIES = (
    BASE.RUNNER_PATH,
    *SUPPORT_DEPENDENCIES,
    ATLAS_MANIFEST_PATH.resolve(),
    ATLAS_SUMMARY_PATH.resolve(),
    ATLAS_GENERATOR_PATH.resolve(),
)
BASE.CANONICAL_MAIN_SHA256 = \
    "78575e2fe1b9796d10dcb4c14d2bcc1b627a79e787be66d34dc7aa666899aa93"


# The inherited collector normally retains every full child payload in RAM.
# This campaign launches 14,550 children, so journal a completed round only
# after its isolation snapshot and retain a compact in-memory reference.
_SHARED_RUN_ONE = BASE.run_one
_BASE_ISOLATION_RECORD = BASE.MAIN_SUPPORT.isolation_record
_PENDING_INVOCATIONS: list[dict[str, Any]] = []
_PENDING_OUTPUT: Path | None = None
_ROUND_ARTIFACT_SEQUENCE = 0
JOURNAL_SETTLE_SECONDS = 0.02


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def run_one_with_round_journal(
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
    global _PENDING_OUTPUT
    BASE.require(not _PENDING_INVOCATIONS or
                 _PENDING_OUTPUT == failure_output,
                 "round journal mixed campaign outputs")
    try:
        record = _SHARED_RUN_ONE(
            implementation, identity, cell, cpu, source_commit, source_tree,
            iterations, warmup, failure_output)
    except Exception as error:
        _persist_pending_round({
            "accepted": False,
            "partial_attempt": True,
            "failure": {"type": type(error).__name__, "message": str(error)},
        })
        raise
    _PENDING_OUTPUT = failure_output
    _PENDING_INVOCATIONS.append(record)
    return record


def _persist_pending_round(isolation: Mapping[str, Any]) -> None:
    global _PENDING_OUTPUT, _ROUND_ARTIFACT_SEQUENCE
    if not _PENDING_INVOCATIONS:
        return
    BASE.require(_PENDING_OUTPUT is not None,
                 "round journal output is absent")
    directory = _PENDING_OUTPUT / "round_payloads"
    directory.mkdir(exist_ok=True)
    path = directory / f"round-{_ROUND_ARTIFACT_SEQUENCE:06d}.json"
    full_record = {
        "schema": "leopard2-current-atlas-round-payload/v1",
        "sequence": _ROUND_ARTIFACT_SEQUENCE,
        "isolation": dict(isolation),
        "invocations": _PENDING_INVOCATIONS,
    }
    BASE.T8_SUPPORT.write_atomic_exclusive(path, full_record)
    _fsync_directory(directory)
    artifact = BASE.support_file_identity(path)
    for index, invocation in enumerate(_PENDING_INVOCATIONS):
        normalized = invocation["normalized"]
        normalized.pop("encode_samples_us", None)
        normalized.pop("one_shot_encode_samples_us", None)
        invocation.pop("result", None)
        invocation.pop("command", None)
        invocation["round_payload_artifact"] = {
            **artifact,
            "invocation_index": index,
        }
    _PENDING_INVOCATIONS.clear()
    _PENDING_OUTPUT = None
    _ROUND_ARTIFACT_SEQUENCE += 1


def isolation_record_with_round_journal(*args: Any, **kwargs: Any) -> Any:
    try:
        record = _BASE_ISOLATION_RECORD(*args, **kwargs)
    except Exception as error:
        _persist_pending_round({
            "accepted": False,
            "partial_attempt": True,
            "failure": {"type": type(error).__name__, "message": str(error)},
        })
        raise
    _persist_pending_round(record)
    if _ROUND_ARTIFACT_SEQUENCE > 0:
        time.sleep(JOURNAL_SETTLE_SECONDS)
    return record


BASE.run_one = run_one_with_round_journal
BASE.MAIN_SUPPORT.isolation_record = isolation_record_with_round_journal


def _old_screen_classifications() -> list[tuple[str, str, str]]:
    """Return exact old-screen ID, operation, and historical role triples."""
    result: list[tuple[str, str, str]] = []

    def add(
        operation: str, k: int, shard_bytes: int, loss: int, role: str,
    ) -> None:
        result.append((
            f"{operation}-k{k}-r32-b{shard_bytes}-l{loss}", operation, role))

    for k in (33, 35, 41, 43, 45, 47, 51, 57, 59, 61, 63, 64, 79, 91):
        historical_role = "former_major" if k in (33, 63, 64) \
            else "former_marginal"
        add("encode", k, 64, 1, historical_role)
    add("encode", 35, 1024, 1, "former_marginal")

    for k in (65, 67, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95):
        add("decode", k, 64, 32, "former_decode")
    for k, loss in (
        (67, 7), (69, 7), (71, 8), (73, 8), (75, 8),
        (77, 8), (79, 8), (81, 9), (83, 9), (85, 9),
        (87, 9), (89, 9), (91, 10), (93, 10), (95, 10),
    ):
        add("decode", k, 1024, loss, "former_decode")
    for k in (81, 83, 85, 95):
        add("decode", k, 1024, 2, "former_decode")

    for k in (32, 34, 36, 62, 65, 80, 92, 96):
        add("encode", k, 64, 1, "encode_neighbor")
    for k in (34, 36):
        add("encode", k, 1024, 1, "encode_neighbor")
    for k in (64, 66, 96, 97):
        add("decode", k, 64, 32, "decode_neighbor")
    for k, loss in ((65, 7), (66, 7), (96, 10), (97, 10)):
        add("decode", k, 1024, loss, "decode_neighbor")
    for k in (79, 80, 86, 96, 97):
        add("decode", k, 1024, 2, "decode_neighbor")
    for k in (65, 81, 95, 96, 97):
        add("decode", k, 1024, 1, "one_loss_control")
    BASE.require(len(result) == 77, "old current-atlas screen changed")
    return result


SUPPLEMENTAL_TARGET_IDS = (
    "k033_b0000064_l02", "k033_b0000064_l04", "k033_b0000064_l32",
    "k035_b0000064_l04", "k035_b0000064_l32", "k037_b0000064_l02",
    "k039_b0000064_l02", "k043_b0000064_l05", "k045_b0000064_l05",
    "k051_b0000064_l02", "k053_b0000064_l02", "k053_b0000064_l32",
    "k059_b0000064_l02", "k059_b0000064_l06", "k059_b0000064_l32",
    "k061_b0000064_l07", "k063_b0000064_l07", "k064_b0000064_l07",
    "k081_b0000064_l02", "k209_b0000064_l32",
)


def _atlas_id_from_screen_id(screen_id: str) -> str:
    match = re.fullmatch(
        r"(?:encode|decode)-k([0-9]+)-r32-b([0-9]+)-l([0-9]+)",
        screen_id)
    BASE.require(match is not None, f"malformed historical cell ID: {screen_id}")
    k, shard_bytes, loss = (int(value) for value in match.groups())
    return f"k{k:03d}_b{shard_bytes:07d}_l{loss:02d}"


def _read_json(path: Path) -> Mapping[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as source:
            value = json.load(source)
    except (OSError, json.JSONDecodeError) as error:
        raise BASE.EvidenceError(f"cannot read atlas JSON {path}: {error}") \
            from error
    BASE.require(isinstance(value, Mapping), f"atlas JSON is not an object: {path}")
    return value


def load_atlas(
    manifest_path: Path = ATLAS_MANIFEST_PATH,
    summary_path: Path = ATLAS_SUMMARY_PATH,
) -> tuple[Mapping[str, Any], Mapping[str, Any]]:
    """Load, hash, and algebraically validate the committed atlas pair."""
    BASE.require(BASE.sha256(manifest_path) == ATLAS_MANIFEST_SHA256,
                 "committed atlas manifest bytes changed")
    BASE.require(BASE.sha256(summary_path) == ATLAS_SUMMARY_SHA256,
                 "committed atlas summary bytes changed")
    manifest = _read_json(manifest_path)
    summary = _read_json(summary_path)
    try:
        ATLAS_GENERATOR.validate_manifest(manifest)
        ATLAS_GENERATOR.validate_summary(summary, manifest,
                                         require_complete=True)
    except Exception as error:
        raise BASE.EvidenceError(f"committed atlas contract failed: {error}") \
            from error
    return manifest, summary


ATLAS_MANIFEST, ATLAS_SUMMARY = load_atlas()


def _classification() -> list[tuple[str, str, str, str]]:
    result = []
    for index, (screen_id, operation, historical_role) in enumerate(
            _old_screen_classifications()):
        role = "target" if index < 49 else "neighbor"
        result.append((screen_id, _atlas_id_from_screen_id(screen_id),
                       operation, role + ":" + historical_role))
    result.extend((atlas_id, atlas_id, "encode",
                   "target:supplemental_atlas_loss")
                  for atlas_id in SUPPLEMENTAL_TARGET_IDS)
    return result


def _derived_control(atlas_id: str) -> dict[str, Any]:
    match = re.fullmatch(r"k([0-9]{3})_b([0-9]{7})_l([0-9]{2})", atlas_id)
    BASE.require(match is not None, f"malformed derived atlas ID: {atlas_id}")
    k, shard_bytes, loss = (int(value) for value in match.groups())
    reuse, _iterations, _warmup = ATLAS_GENERATOR.repetition_policy(
        k, shard_bytes)
    return {
        "id": atlas_id, "K": k, "R": 32,
        "shard_bytes": shard_bytes, "loss_count": loss,
        "seed": ATLAS_GENERATOR.seed_for(k, shard_bytes, loss),
        "reuse": reuse,
    }


def _select_losses(k: int, loss_count: int, seed: int) -> list[int]:
    """Mirror benchmark.cpp's uint64_t xorshift/Fisher-Yates selector."""
    mask = (1 << 64) - 1
    state = (seed ^ 0xd1b54a32d192ed03) & mask
    if state == 0:
        state = 0x9e3779b97f4a7c15

    def next_value() -> int:
        nonlocal state
        value = state
        value = (value ^ ((value << 13) & mask)) & mask
        value = (value ^ (value >> 7)) & mask
        value = (value ^ ((value << 17) & mask)) & mask
        state = value
        return value

    order = list(range(k))
    for remaining in range(k, 1, -1):
        selected = next_value() % remaining
        order[remaining - 1], order[selected] = \
            order[selected], order[remaining - 1]
    return sorted(order[:loss_count])


def cells() -> list[dict[str, Any]]:
    """Return 69 historical loss targets plus 28 scientific controls."""
    manifest_cells = {item["id"]: item for item in ATLAS_MANIFEST["cells"]}
    summary_rows = {
        (item["cell_id"], item["codec"]): item
        for item in ATLAS_SUMMARY["rows"]
    }
    result: list[dict[str, Any]] = []
    for screen_id, atlas_id, operation, combined_role in _classification():
        role, historical_role = combined_role.split(":", 1)
        source = manifest_cells.get(atlas_id)
        origin = "committed_atlas_manifest"
        if source is None:
            BASE.require(role == "neighbor",
                         f"target is absent from committed atlas: {atlas_id}")
            source = _derived_control(atlas_id)
            origin = "canonical_policy_derived_control"
            expected_missing = _select_losses(
                source["K"], source["loss_count"], source["seed"])
        else:
            expected_missing = None
            for codec in ("leopard2", "leopard1"):
                row = summary_rows.get((atlas_id, codec))
                BASE.require(isinstance(row, Mapping) and
                             all(row.get(name) == source[name] for name in
                                 ("K", "R", "shard_bytes", "loss_count",
                                  "seed")),
                             f"atlas summary row differs: {atlas_id}/{codec}")
                missing = row.get("missing_original_indices")
                BASE.require(
                    isinstance(missing, list) and
                    all(type(index) is int for index in missing) and
                    (expected_missing is None or missing == expected_missing),
                    f"atlas erasure pattern differs: {atlas_id}/{codec}")
                expected_missing = list(missing)
            BASE.require(expected_missing == _select_losses(
                source["K"], source["loss_count"], source["seed"]),
                f"atlas erasure pattern is not reproducible: {atlas_id}")
        result.append({
            "id": screen_id,
            "atlas_cell_id": atlas_id,
            "K": source["K"], "R": source["R"],
            "bytes": source["shard_bytes"],
            "loss": source["loss_count"],
            "batch": 1, "reuse": source["reuse"],
            "role": role, "historical_role": historical_role,
            "operation": operation, "workload_origin": origin,
            "compare_main": True, "measure_one_shot": True,
            "seed": source["seed"],
            "missing_original_indices": expected_missing,
        })

    BASE.require(
        len(result) == 97 and
        sum(item["role"] == "target" for item in result) == 69 and
        sum(item["role"] == "neighbor" for item in result) == 28 and
        sum(item["role"] == "target" and item["operation"] == "encode"
            for item in result) == 35 and
        sum(item["role"] == "target" and item["operation"] == "decode"
            for item in result) == 34 and
        len({item["id"] for item in result}) == len(result) and
        sum(item["workload_origin"] == "committed_atlas_manifest"
            for item in result) == 81 and
        sum(item["workload_origin"] == "canonical_policy_derived_control"
            for item in result) == 16 and
        all(item["R"] == 32 and item["batch"] == 1 and item["reuse"] >= 1
            and item["compare_main"] is True and
            item["measure_one_shot"] is True and
            len(item["missing_original_indices"]) == item["loss"]
            for item in result),
        "historical current-atlas qualification matrix is incomplete")
    return result


BASE.cells = cells


def benchmark_command(
    implementation: str,
    executable: Path,
    cell: Mapping[str, Any],
    cpu: int,
    iterations: int,
    warmup: int,
) -> list[str]:
    """Measure encode and setup-included decode without diagnostic routes."""
    common = [
        "/usr/bin/prlimit", "--as=201326592",
        "/usr/bin/taskset", "-c", str(cpu), str(executable),
        "--k", str(cell["K"]), "--r", str(cell["R"]),
        "--bytes", str(cell["bytes"]), "--loss", str(cell["loss"]),
        "--batch", "1", "--reuse", str(cell["reuse"]),
        "--iterations", str(iterations), "--warmup", str(warmup),
        "--threads", "1", "--seed", str(cell["seed"]),
    ]
    if implementation == "main":
        return common + ["--json", "-"]
    BASE.require(implementation in {"candidate", "control"},
                 f"unknown benchmark implementation: {implementation}")
    return common + [
        "--profile", "high", "--field", "gf8", "--backend", "avx2",
        "--skip-legacy", "--retain-samples", "--measure-one-shot-decode",
        "--json", "-",
    ]


BASE.benchmark_command = benchmark_command


def _positive_metric(
    result: Mapping[str, Any], metric_name: str, iterations: int,
) -> dict[str, Any]:
    metrics = result.get("metrics")
    BASE.require(isinstance(metrics, Mapping), "benchmark metrics are absent")
    metric = metrics.get(metric_name)
    BASE.require(isinstance(metric, Mapping), f"{metric_name} metric is absent")
    median = metric.get("median_us_per_batch_call")
    mad = metric.get("mad_us_per_batch_call")
    minimum = metric.get("minimum_us_per_batch_call")
    maximum = metric.get("maximum_us_per_batch_call")
    samples = metric.get("samples_us_per_batch_call")
    BASE.require(
        isinstance(samples, list) and len(samples) == iterations and
        all(isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(float(value)) and float(value) > 0
            for value in samples),
        f"{metric_name} raw samples are incomplete")
    numeric = [float(value) for value in samples]
    expected_median = statistics.median(numeric)
    expected_mad = statistics.median(
        abs(value - expected_median) for value in numeric)
    BASE.require(
        all(isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(float(value)) for value in
            (median, mad, minimum, maximum)) and
        expected_median > 0 and float(mad) >= 0 and
        math.isclose(float(median), expected_median,
                     rel_tol=0.0, abs_tol=2e-6) and
        math.isclose(float(mad), expected_mad,
                     rel_tol=0.0, abs_tol=2e-6) and
        math.isclose(float(minimum), min(numeric),
                     rel_tol=0.0, abs_tol=2e-6) and
        math.isclose(float(maximum), max(numeric),
                     rel_tol=0.0, abs_tol=2e-6),
        f"{metric_name} summary disagrees with retained raw samples")
    return {"median_us": expected_median, "samples_us": numeric}


def _power_of_two_at_least(value: int) -> int:
    result = 1
    while result < value:
        result *= 2
    return result


def _require_integer_fields(
    payload: Mapping[str, Any], expected: Mapping[str, int], label: str,
) -> None:
    BASE.require(
        all(type(payload.get(name)) is int and payload[name] == value
            for name, value in expected.items()),
        f"{label} dimensions differ from the frozen cell")


def _require_digests(result: Mapping[str, Any]) -> dict[str, Any]:
    digests = result.get("workload_digests")
    names = ("original_data", "transmitted_parity", "recovered_originals")
    BASE.require(
        isinstance(digests, Mapping) and digests.get("algorithm") == "fnv1a64"
        and all(isinstance(digests.get(name), str) and
                re.fullmatch(r"[0-9a-f]{16}", digests[name]) is not None
                for name in names),
        "benchmark workload digests are incomplete or malformed")
    return dict(digests)


def validate_result(
    implementation: str,
    result: object,
    cell: Mapping[str, Any],
    source_commit: str,
    source_tree: str,
    iterations: int,
    warmup: int,
) -> dict[str, Any]:
    """Strictly validate v9/current and canonical exact-main payloads."""
    del source_commit, source_tree  # Bound by the production build closure.
    BASE.require(isinstance(result, Mapping),
                 "benchmark output is not an object")
    parameters = result.get("parameters")
    resolved = result.get("resolved")
    correctness = result.get("correctness")
    build = result.get("build")
    memory = result.get("memory")
    BASE.require(all(isinstance(value, Mapping) for value in
                     (parameters, resolved, correctness, build, memory)),
                 "benchmark identity sections are malformed")
    _require_integer_fields(parameters, {
        "K": int(cell["K"]), "R": int(cell["R"]),
        "shard_bytes": int(cell["bytes"]),
        "loss_count": int(cell["loss"]), "batch": 1,
        "reuse": int(cell["reuse"]),
        "iterations": iterations, "warmup": warmup,
        "thread_count": 1, "seed": int(cell["seed"]),
    }, "benchmark")
    missing = parameters.get("missing_original_indices")
    BASE.require(
        isinstance(missing, list) and
        missing == cell["missing_original_indices"],
        "missing-original indices differ from the frozen erasure pattern")
    expected_parent = _power_of_two_at_least(cell["K"] + 32)
    _require_integer_fields(resolved, {
        "thread_count": 1, "parent_count": expected_parent,
        "padded_side": 32,
    }, "resolved codec")
    BASE.require(resolved.get("profile") == "legacy_high_v1" and
                 resolved.get("field") == "gf8",
                 "benchmark resolved another wire profile or field")
    digests = _require_digests(result)

    if implementation == "main":
        BASE.require(
            result.get("schema") == "leopard-main-benchmark-v1" and
            parameters.get("logical_shard_bytes") == cell["bytes"] and
            correctness.get("round_trip") is True and
            correctness.get("logical_prefix_fingerprinted") is True and
            resolved.get("padded_application_bytes") is False and
            resolved.get("padding_policy") == "zero suffix per shard" and
            build.get("main_source_commit") == BASE.MAIN_COMMIT,
            "exact-main schema, source identity, or round trip changed")
        metrics = result.get("metrics")
        BASE.require(isinstance(metrics, Mapping) and
                     metrics.get("decode_timing_includes_setup") is True,
                     "exact-main decode does not include setup")
        _require_integer_fields(memory, {"alignment": 64}, "main memory")
        for name in ("encode_work_shards_per_stripe",
                     "decode_work_shards_per_stripe",
                     "encode_work_bytes_batch", "decode_work_bytes_batch"):
            BASE.require(type(memory.get(name)) is int and memory[name] > 0,
                         f"main {name} is malformed")
        BASE.require(
            memory["encode_work_bytes_batch"] ==
                memory["encode_work_shards_per_stripe"] * cell["bytes"] and
            memory["decode_work_bytes_batch"] ==
                memory["decode_work_shards_per_stripe"] * cell["bytes"],
            "main work-byte geometry is inconsistent")
        encode_metric_name = "encode_execution"
        decode_metric_name = "decode_including_setup"
    else:
        BASE.require(implementation in {"candidate", "control"},
                     f"unknown implementation: {implementation}")
        BASE.require(
            result.get("schema") == "leopard2-benchmark-v9" and
            parameters.get("requested_profile") == "legacy_high_v1" and
            parameters.get("requested_field") == "gf8" and
            parameters.get("requested_backend") == "avx2" and
            parameters.get("force_generic_decode") is False and
            parameters.get("force_specialized_decode") is False and
            parameters.get("force_tiled_decode") is False and
            parameters.get("force_materialized_decode") is False and
            parameters.get("skip_legacy") is True and
            parameters.get("retain_samples") is True and
            parameters.get("measure_one_shot_decode") is True and
            "measure_one_shot_encode" not in parameters and
            "attest_source" not in parameters and
            resolved.get("backend") == "avx2" and
            correctness.get("leopard2_round_trip") is True and
            correctness.get("legacy_comparison") is None,
            "Leopard2 requested/resolved identity or round trip changed")
        BASE.require(
            all(name not in build for name in
                ("source_commit", "source_tree", "source_tracked_dirty")),
            "standalone v9 payload unexpectedly contains source attestation")
        _require_integer_fields(memory, {"scratch_alignment": 64},
                                "Leopard2 memory")
        for prefix in ("encode_scratch_bytes", "decode_scratch_bytes",
                       "one_shot_decode_scratch_bytes"):
            per_stripe = memory.get(prefix + "_per_stripe")
            batch = memory.get(prefix + "_batch")
            BASE.require(type(per_stripe) is int and per_stripe > 0 and
                         type(batch) is int and batch == per_stripe,
                         f"Leopard2 {prefix} geometry is malformed")
        encode_metric_name = "encode_execution"
        decode_metric_name = "one_shot_decode_including_setup"

    encode = _positive_metric(result, encode_metric_name, iterations)
    decode = _positive_metric(result, decode_metric_name, iterations)
    if cell["operation"] == "encode":
        primary, companion = encode, decode
    else:
        primary, companion = decode, encode
    return {
        "encode_us": primary["median_us"],
        "encode_samples_us": primary["samples_us"],
        "one_shot_encode_us": companion["median_us"],
        "one_shot_encode_samples_us": companion["samples_us"],
        "digests": digests,
        "missing_original_indices": list(missing),
        "schema": result["schema"],
    }


BASE.validate_result = validate_result


# Bypass the older terminal wrapper's selector-specific equivalence assertions:
# this screen has no selector.  Performance acceptance is evaluated by the
# parent summary after raw.json has been durably written.
_BASE_ANALYZE = PARENT.SUPPORT._BASE_ANALYZE


def _mapping_for(cell: Mapping[str, Any]) -> dict[str, Any]:
    current_decode = "one_shot_decode_including_setup"
    main_decode = "decode_including_setup"
    if cell["operation"] == "encode":
        return {
            "primary_normalized_field": "encode_us",
            "primary_operation": "encode",
            "primary_current_metric": "encode_execution",
            "primary_main_metric": "encode_execution",
            "companion_normalized_field": "one_shot_encode_us",
            "companion_operation": "decode",
            "companion_current_metric": current_decode,
            "companion_main_metric": main_decode,
        }
    return {
        "primary_normalized_field": "encode_us",
        "primary_operation": "decode",
        "primary_current_metric": current_decode,
        "primary_main_metric": main_decode,
        "companion_normalized_field": "one_shot_encode_us",
        "companion_operation": "encode",
        "companion_current_metric": "encode_execution",
        "companion_main_metric": "encode_execution",
    }


def analyze(
    cell: Mapping[str, Any], rounds: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Annotate historical meaning without pre-raw performance exceptions."""
    expected_missing = list(cell["missing_original_indices"])
    BASE.require(all(
        invocation["normalized"].get("missing_original_indices") ==
            expected_missing
        for round_value in rounds
        for invocation in round_value["invocations"]),
        "candidate/control/main erasure patterns differ")
    result = _BASE_ANALYZE(cell, rounds)
    result["historical_role"] = cell["historical_role"]
    result["atlas_cell_id"] = cell["atlas_cell_id"]
    result["workload_origin"] = cell["workload_origin"]
    result["metric_mapping"] = _mapping_for(cell)
    result["missing_original_indices"] = expected_missing
    replicate = {
        "primary": result["control_over_candidate"],
        "companion": result["one_shot_control_over_candidate"],
    }
    main = {
        "primary": result["main_over_candidate"],
        "companion": result["one_shot_main_over_candidate"],
    }
    result["acceptance_observations"] = {
        "same_binary_ci95_equivalence_band": [
            BASE.NEIGHBOR_EQUIVALENCE_LOWER,
            BASE.NEIGHBOR_EQUIVALENCE_UPPER,
        ],
        "same_binary_equivalence_pass": {
            name: (ratio["ci95"][0] >=
                   BASE.NEIGHBOR_EQUIVALENCE_LOWER and
                   ratio["ci95"][1] <= BASE.NEIGHBOR_EQUIVALENCE_UPPER)
            for name, ratio in replicate.items()
        },
        "exact_main_lower_ci_floor": BASE.TARGET_MAIN_FLOOR,
        "exact_main_lower_ci_pass": {
            name: ratio["ci95"][0] >= BASE.TARGET_MAIN_FLOOR
            for name, ratio in main.items()
        },
    }
    return result


BASE.analyze = analyze


def _validated_ci(ratio: Mapping[str, Any], label: str) -> tuple[float, float]:
    interval = ratio.get("ci95")
    BASE.require(
        isinstance(interval, list) and len(interval) == 2 and
        all(isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(float(value)) and float(value) > 0
            for value in interval) and interval[0] <= interval[1],
        f"{label} confidence interval is malformed")
    return float(interval[0]), float(interval[1])


def validate_round_payload_artifacts(
    raw: Mapping[str, Any], summary: Mapping[str, Any], output: Path,
    expected_cell_count: int = 97,
) -> dict[str, Any]:
    """Reopen every full round payload and bind it to compact raw records."""
    raw_cells = raw.get("cells")
    strict_production_matrix = expected_cell_count == 97
    expected_cells = cells() if strict_production_matrix else None
    expected_orders = BASE.select_round_orders(
        BASE.TARGET_ORDER, 25) if strict_production_matrix else None
    BASE.require(raw.get("schema") == BASE.SCHEMA and
                 isinstance(raw_cells, list) and
                 len(raw_cells) == expected_cell_count,
                 "compact raw campaign matrix is incomplete")
    expected_sequence = 0
    accepted_processes = 0
    discarded_processes = 0
    artifact_hashes = []
    seen_paths = set()
    identities = raw.get("identities")
    BASE.require(isinstance(identities, Mapping),
                 "compact raw executable identities are absent")
    for cell_index, raw_cell in enumerate(raw_cells):
        cell = raw_cell.get("cell")
        rounds = raw_cell.get("rounds")
        BASE.require(isinstance(cell, Mapping) and
                     isinstance(rounds, list) and
                     (not strict_production_matrix or
                      (cell == expected_cells[cell_index] and
                       len(rounds) == 25)),
                     "compact raw cell rounds are malformed")
        for round_index, accepted_round in enumerate(rounds):
            discarded = accepted_round.get("discarded_attempts")
            expected_order = list(expected_orders[round_index]) \
                if strict_production_matrix else accepted_round.get("order")
            BASE.require(isinstance(discarded, list) and
                         (not strict_production_matrix or
                          (accepted_round.get("round") == round_index and
                           accepted_round.get("order") == expected_order and
                           accepted_round.get("isolation", {}).get(
                               "accepted") is True)),
                         "discarded-attempt list is malformed")
            for is_discarded, attempt in (
                [(True, item) for item in discarded] +
                [(False, accepted_round)]
            ):
                invocations = attempt.get("invocations")
                isolation = attempt.get("isolation")
                BASE.require(isinstance(invocations, list) and invocations and
                             isinstance(isolation, Mapping) and
                             attempt.get("order") == expected_order and
                             [item.get("implementation")
                              for item in invocations] == expected_order and
                             (not strict_production_matrix or
                              isolation.get("accepted") is
                                (not is_discarded)),
                             "compact round attempt is malformed")
                references = [item.get("round_payload_artifact")
                              for item in invocations]
                BASE.require(all(isinstance(item, Mapping)
                                 for item in references),
                             "compact invocation lacks its full payload")
                artifact_identity = {
                    name: references[0].get(name)
                    for name in ("path", "size", "sha256")
                }
                expected_path = (
                    output / "round_payloads" /
                    f"round-{expected_sequence:06d}.json").resolve()
                BASE.require(
                    all({name: reference.get(name) for name in
                         ("path", "size", "sha256")} == artifact_identity
                        for reference in references) and
                    [reference.get("invocation_index")
                     for reference in references] ==
                        list(range(len(invocations))) and
                    artifact_identity["path"] == str(expected_path) and
                    artifact_identity["path"] not in seen_paths,
                    "round payload references are duplicated or out of order")
                observed_identity = BASE.support_file_identity(expected_path)
                BASE.require(observed_identity == artifact_identity,
                             "round payload identity changed")
                payload = _read_json(expected_path)
                full_invocations = payload.get("invocations")
                BASE.require(
                    payload.get("schema") ==
                        "leopard2-current-atlas-round-payload/v1" and
                    payload.get("sequence") == expected_sequence and
                    payload.get("isolation") == isolation and
                    isinstance(full_invocations, list) and
                    len(full_invocations) == len(invocations),
                    "full round payload differs from compact evidence")
                for compact, full in zip(invocations, full_invocations):
                    BASE.require(isinstance(full.get("result"), Mapping) and
                                 isinstance(full.get("command"), list),
                                 "full child result or command is absent")
                    implementation = full.get("implementation")
                    identity = identities.get(implementation)
                    BASE.require(
                        implementation in {"candidate", "control", "main"} and
                        isinstance(identity, Mapping),
                        "full child implementation identity is absent")
                    expected_command = benchmark_command(
                        implementation, Path(str(identity["path"])), cell,
                        raw["cpu"], raw["iterations"], raw["warmup"])
                    expected_normalized_full = validate_result(
                        implementation, full["result"], cell,
                        raw["source_commit"], raw["source_tree"],
                        raw["iterations"], raw["warmup"])
                    BASE.require(full["command"] == expected_command and
                                 full.get("normalized") ==
                                    expected_normalized_full,
                                 "full child command or result was not reproducible")
                    expected_compact = {
                        name: value for name, value in full.items()
                        if name not in {"result", "command"}
                    }
                    expected_normalized = dict(expected_compact["normalized"])
                    expected_normalized.pop("encode_samples_us", None)
                    expected_normalized.pop(
                        "one_shot_encode_samples_us", None)
                    expected_compact["normalized"] = expected_normalized
                    observed_compact = {
                        name: value for name, value in compact.items()
                        if name != "round_payload_artifact"
                    }
                    BASE.require(observed_compact == expected_compact,
                                 "compact invocation differs from full payload")
                seen_paths.add(artifact_identity["path"])
                artifact_hashes.append(artifact_identity["sha256"])
                if is_discarded:
                    discarded_processes += len(invocations)
                else:
                    accepted_processes += len(invocations)
                expected_sequence += 1
    payload_directory = output / "round_payloads"
    directory_entries = list(payload_directory.iterdir()) \
        if payload_directory.is_dir() else []
    actual_paths = {
        str(path.resolve()) for path in directory_entries
    }
    BASE.require(
        summary.get("process_count") == accepted_processes and
        summary.get("discarded_process_count") == discarded_processes and
        all(path.is_file() and not path.is_symlink()
            for path in directory_entries) and
        len(directory_entries) == expected_sequence and
        actual_paths == seen_paths and len(seen_paths) == expected_sequence and
        (not strict_production_matrix or
         summary.get("cells") == [
             analyze(raw_cell["cell"], raw_cell["rounds"])
             for raw_cell in raw_cells]),
        "round payload counts differ from the summary")
    combined = hashlib.sha256(
        ("\n".join(artifact_hashes) + "\n").encode("ascii")).hexdigest()
    return {
        "schema": "leopard2-current-atlas-round-payload-closure/v1",
        "artifact_count": expected_sequence,
        "accepted_process_count": accepted_processes,
        "discarded_process_count": discarded_processes,
        "ordered_artifact_sha256": combined,
    }


def validate_final_live_identities(
    raw: Mapping[str, Any], options: Any,
) -> None:
    """Close the brief base-return/finalizer-reacquire identity window."""
    BASE.require(
        BASE.T8_SUPPORT.file_identity(BASE.RUNNER_PATH) ==
            raw.get("runner_after") and
        [BASE.support_file_identity(path)
         for path in BASE.RUNNER_DEPENDENCIES] ==
            raw.get("runner_dependencies_after"),
        "runner or support identity changed before finalization")
    for collection_name in ("identities_after", "input_identities_after"):
        collection = raw.get(collection_name)
        BASE.require(isinstance(collection, Mapping),
                     f"{collection_name} is absent")
        for identity in collection.values():
            BASE.require(isinstance(identity, Mapping) and
                         BASE.T8_SUPPORT.file_identity(
                             Path(str(identity["path"]))) == identity,
                         f"{collection_name} changed before finalization")
    BASE.require(BASE.build_closure_identity(options) ==
                 raw.get("build_closure_after"),
                 "production build closure changed before finalization")


def apply_final_acceptance(summary: dict[str, Any]) -> dict[str, Any]:
    """Apply all-cell main floors and a symmetric same-binary noise gate."""
    BASE.require(summary.get("schema") == BASE.SUMMARY_SCHEMA,
                 "final summary schema changed")
    analyses = summary.get("cells")
    expected_cells = cells()
    BASE.require(isinstance(analyses, list) and len(analyses) == 97 and
                 [analysis.get("cell") for analysis in analyses] ==
                    expected_cells and
                 len({analysis["cell"]["id"] for analysis in analyses}) == 97,
                 "final summary does not contain all 97 cells")
    exact_main_failures = []
    replicate_failures = []
    for analysis in analyses:
        cell = analysis.get("cell")
        mapping = analysis.get("metric_mapping")
        BASE.require(isinstance(cell, Mapping) and
                     isinstance(mapping, Mapping) and
                     mapping == _mapping_for(cell) and
                     analysis.get("historical_role") ==
                        cell.get("historical_role") and
                     analysis.get("atlas_cell_id") ==
                        cell.get("atlas_cell_id") and
                     analysis.get("workload_origin") ==
                        cell.get("workload_origin") and
                     analysis.get("missing_original_indices") ==
                        cell.get("missing_original_indices"),
                     "final cell analysis is malformed")
        for position, control_name, main_name in (
            ("primary", "control_over_candidate", "main_over_candidate"),
            ("companion", "one_shot_control_over_candidate",
             "one_shot_main_over_candidate"),
        ):
            operation = mapping.get(position + "_operation")
            BASE.require(operation in {"encode", "decode"},
                         "final metric operation is malformed")
            control_ci = _validated_ci(
                analysis.get(control_name, {}),
                f"{cell.get('id')} same-binary {operation}")
            main_ci = _validated_ci(
                analysis.get(main_name, {}),
                f"{cell.get('id')} exact-main {operation}")
            if not (BASE.NEIGHBOR_EQUIVALENCE_LOWER <= control_ci[0] and
                    control_ci[1] <= BASE.NEIGHBOR_EQUIVALENCE_UPPER):
                replicate_failures.append({
                    "cell_id": cell.get("id"), "role": cell.get("role"),
                    "operation": operation, "ci95": list(control_ci),
                })
            if main_ci[0] < BASE.TARGET_MAIN_FLOOR:
                exact_main_failures.append({
                    "cell_id": cell.get("id"), "role": cell.get("role"),
                    "operation": operation, "ci95": list(main_ci),
                })

    for name in (
        "target_control_failure", "target_main_failure",
        "target_control_failure_by_metric", "target_main_failure_by_metric",
        "credible_neighbor_regressions",
        "credible_neighbor_one_shot_regressions",
    ):
        summary.pop(name, None)
    summary["acceptance_contract"] = {
        "exact_main_scope": "all_97_cells_both_encode_and_decode",
        "exact_main_ci95_lower_floor": BASE.TARGET_MAIN_FLOOR,
        "same_binary_scope": "all_97_cells_both_encode_and_decode",
        "same_binary_ci95_equivalence_band": [
            BASE.NEIGHBOR_EQUIVALENCE_LOWER,
            BASE.NEIGHBOR_EQUIVALENCE_UPPER,
        ],
    }
    summary["exact_main_regressions"] = exact_main_failures
    summary["same_binary_equivalence_failures"] = replicate_failures
    summary["status"] = "accepted" if not (
        exact_main_failures or replicate_failures) else "rejected"
    summary["schema"] = FINAL_SUMMARY_SCHEMA
    return summary


def _replace_json(path: Path, value: Mapping[str, Any]) -> None:
    temporary = path.with_name(path.name + ".finalizing")
    BASE.require(not temporary.exists(),
                 f"stale summary finalization file exists: {temporary}")
    BASE.T8_SUPPORT.write_exclusive(temporary, value)
    os.replace(temporary, path)
    _fsync_directory(path.parent)


def _bind_failure_journal(output: Path) -> None:
    failure_path = output / "failure.json"
    if not failure_path.exists():
        return
    failure = dict(_read_json(failure_path))
    artifacts = sorted((output / "round_payloads").glob("round-*.json")) \
        if (output / "round_payloads").is_dir() else []
    failure["round_payload_artifacts"] = [
        BASE.support_file_identity(path) for path in artifacts]
    _replace_json(failure_path, failure)


def _begin_isolated_finalization(
    benchmark_cpu: int, sibling_cpu: int,
) -> tuple[int, int]:
    """Leave the timed pair, then serialize evidence finalization globally."""
    finalization_cpu = None
    for candidate in range(os.cpu_count() or 0):
        if candidate in {benchmark_cpu, sibling_cpu}:
            continue
        try:
            os.sched_setaffinity(0, {candidate})
        except OSError:
            continue
        if os.sched_getaffinity(0) == {candidate}:
            finalization_cpu = candidate
            break
    BASE.require(finalization_cpu is not None,
                 "no CPU outside the timed SMT pair is available for finalization")
    descriptor = os.open(BASE.LOCK_PATH, os.O_RDWR | os.O_CREAT, 0o600)
    try:
        fcntl.flock(descriptor, fcntl.LOCK_EX)
    except Exception:
        os.close(descriptor)
        raise
    return descriptor, finalization_cpu


def main() -> int:
    """Collect through the hardened base, then apply this screen's gates."""
    options = BASE.parse_arguments()
    captured_stdout = io.StringIO()
    with contextlib.redirect_stdout(captured_stdout):
        preliminary_status = BASE.main()
    finalization_lock, finalization_cpu = _begin_isolated_finalization(
        options.cpu, options.sibling)
    try:
        summary_path = options.output / "summary.json"
        if preliminary_status not in (0, 2):
            _bind_failure_journal(options.output)
            output = captured_stdout.getvalue()
            if output:
                print(output, end="")
            return preliminary_status
        if not summary_path.exists():
            _bind_failure_journal(options.output)
            print("evidence rejected: preliminary summary is absent",
                  file=sys.stderr)
            return 1
        summary = dict(_read_json(summary_path))
        summary["preliminary_generic_status"] = summary.get("status")
        summary["finalization_cpu"] = finalization_cpu
        finalization_failed = False
        try:
            raw_path = options.output / "raw.json"
            BASE.require(summary.get("raw_sha256") == BASE.sha256(raw_path),
                         "preliminary summary does not bind compact raw evidence")
            raw = _read_json(raw_path)
            validate_final_live_identities(raw, options)
            summary["round_payload_closure"] = \
                validate_round_payload_artifacts(raw, summary, options.output)
            apply_final_acceptance(summary)
        except Exception as error:
            finalization_failed = True
            summary["schema"] = FINAL_SUMMARY_SCHEMA
            summary["status"] = "rejected"
            summary["finalization_failure"] = {
                "type": type(error).__name__, "message": str(error),
            }
        _replace_json(summary_path, summary)
        print(json.dumps({
            "status": summary["status"],
            "cells": summary["cell_count"],
            "processes": summary["process_count"],
            "discarded_processes": summary["discarded_process_count"],
            "exact_main_regressions": len(
                summary.get("exact_main_regressions", [])),
            "same_binary_equivalence_failures": len(
                summary.get("same_binary_equivalence_failures", [])),
        }, sort_keys=True))
        if finalization_failed:
            return 1
        return 0 if summary["status"] == "accepted" else 2
    finally:
        os.close(finalization_lock)


if __name__ == "__main__":
    raise SystemExit(main())
