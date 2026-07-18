#!/usr/bin/env python3
"""Authenticate and summarize balanced same-binary forced-decoder evidence.

This analyzer deliberately does not accept caller-supplied source commits,
binary digests, or mode labels.  Those identities come only from the signed raw
manifest and are checked against the live clean source/build closure.  Every
timed role must be one exact forced tuple: generic, materialized Algorithm 5,
or tiled Algorithm 5.  An AUTO row cannot be reclassified as specialized.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys
import tempfile
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
import balanced_evidence_common as common  # noqa: E402
sys.path.pop(0)


T95 = {1: 12.706204736, 2: 4.302652730, 3: 3.182446305,
       4: 2.776445105, 5: 2.570581836, 6: 2.446911851,
       7: 2.364624252, 8: 2.306004135, 9: 2.262157163}


def close_six(actual: object, expected: float, label: str) -> None:
    value = common.finite(actual, label)
    common.require(abs(value - expected) <= 2.1e-6,
                   f"{label}={value!r}, recomputed {expected!r}")


def validate_signed(value: dict[str, Any], label: str) -> dict[str, Any]:
    digest = common.hex_digest(value.get("content_sha256"), label + " content digest")
    unsigned = dict(value)
    del unsigned["content_sha256"]
    common.require(common.canonical_sha256(unsigned) == digest,
                   f"{label} content digest mismatch")
    return unsigned


def validate_file_identity(value: object, label: str,
                           expected_relative: str | None = None) -> dict[str, Any]:
    identity = common.require_keys(value, {
        "path", "realpath", "relative_path", "sha256", "size", "mode",
        "device", "inode", "mtime_ns",
    }, label)
    common.require(isinstance(identity["path"], str) and
                   Path(identity["path"]).is_absolute(), f"{label} path is not absolute")
    common.require(isinstance(identity["realpath"], str) and
                   Path(identity["realpath"]).is_absolute(),
                   f"{label} realpath is not absolute")
    relative = identity["relative_path"]
    common.require(relative is None or (isinstance(relative, str) and relative and
                   not Path(relative).is_absolute() and ".." not in Path(relative).parts),
                   f"{label} relative path is unsafe")
    if expected_relative is not None:
        common.require(relative == expected_relative, f"{label} relative path changed")
    common.hex_digest(identity["sha256"], label + " SHA-256")
    common.strict_int(identity["size"], label + " size")
    common.strict_int(identity["mode"], label + " mode")
    common.strict_int(identity["device"], label + " device")
    common.strict_int(identity["inode"], label + " inode", 1)
    common.strict_int(identity["mtime_ns"], label + " mtime", 1)
    return identity


def validate_source(value: object) -> dict[str, Any]:
    source = common.require_keys(value, {
        "root", "head", "tree", "expected_commit", "status", "status_sha256",
    }, "source")
    common.require(isinstance(source["root"], str) and Path(source["root"]).is_absolute(),
                   "source root is not absolute")
    common.require(common.hex_digest(source["head"], "source head", 40) ==
                   common.hex_digest(source["expected_commit"],
                                     "source expected commit", 40),
                   "source head/expected commit differ")
    common.hex_digest(source["tree"], "source tree", 40)
    common.require(source["status"] == "clean" and
                   source["status_sha256"] == common.EMPTY_SHA256,
                   "source identity is not clean")
    return source


def validate_build(value: object, source_root: str) -> dict[str, Any]:
    build = common.require_keys(value, {
        "root", "cmake_home", "cache", "graph", "sources", "objects",
        "archive", "binary",
    }, "build")
    common.require(build["cmake_home"] == source_root and
                   isinstance(build["root"], str) and Path(build["root"]).is_absolute(),
                   "build root/CMake home differs from source")
    validate_file_identity(build["cache"], "CMake cache")
    validate_file_identity(build["graph"], "CMake build graph")
    sources = common.require_keys(build["sources"],
                                  {"benchmark", "decoder", "dispatch"}, "build sources")
    validate_file_identity(sources["benchmark"], "benchmark source",
                           "bench/leopard2/benchmark.cpp")
    validate_file_identity(sources["decoder"], "decoder source", "leopard2.cpp")
    validate_file_identity(sources["dispatch"], "dispatch source", "Leopard2Dispatch.h")
    objects = common.require_keys(build["objects"], {"benchmark", "decoder"},
                                  "build objects")
    validate_file_identity(objects["benchmark"], "benchmark object")
    validate_file_identity(objects["decoder"], "decoder object")
    validate_file_identity(build["archive"], "production archive")
    validate_file_identity(build["binary"], "benchmark executable")
    return build


def live_file(identity: dict[str, Any], label: str) -> None:
    path = Path(identity["path"])
    common.require(path.is_file(), f"live {label} is absent: {path}")
    descriptor = path.stat()
    common.require(str(path.resolve()) == identity["realpath"] and
                   descriptor.st_size == identity["size"] and
                   common.sha256_file(path) == identity["sha256"] and
                   (descriptor.st_mode & 0o7777) == identity["mode"] and
                   descriptor.st_dev == identity["device"] and
                   descriptor.st_ino == identity["inode"] and
                   descriptor.st_mtime_ns == identity["mtime_ns"],
                   f"live {label} metadata/bytes differ: {path}")


def validate_runtime(value: object) -> dict[str, Any]:
    runtime = common.require_keys(value, {
        "platform", "python", "affinity", "cpu", "sibling", "topology",
        "clock_ticks_per_second", "cpuinfo", "runtime_libraries",
        "ldd_normalized_sha256", "environment",
    }, "runtime")
    common.require(isinstance(runtime["platform"], dict) and
                   set(runtime["platform"]) == {"system", "node", "release", "version", "machine"},
                   "runtime platform identity changed")
    python = common.require_keys(runtime["python"],
                                 {"version", "implementation", "executable", "byteorder"},
                                 "runtime Python")
    validate_file_identity(python["executable"], "runtime Python executable")
    affinity = runtime["affinity"]
    common.require(isinstance(affinity, list) and affinity == sorted(set(affinity)) and
                   all(type(item) is int and item >= 0 for item in affinity),
                   "runtime affinity is invalid")
    cpu = common.strict_int(runtime["cpu"], "runtime CPU")
    sibling = common.strict_int(runtime["sibling"], "runtime sibling")
    common.require(cpu != sibling and cpu in affinity and sibling in affinity,
                   "runtime CPU pair is invalid")
    topology = runtime["topology"]
    common.require(isinstance(topology, list) and len(topology) == 2,
                   "runtime topology is not an exact pair")
    common.strict_int(runtime["clock_ticks_per_second"], "clock ticks", 1)
    common.require(isinstance(runtime["cpuinfo"], list) and
                   len(runtime["cpuinfo"]) == 2 and
                   all(isinstance(item, dict) and item for item in runtime["cpuinfo"]),
                   "runtime stable CPU information is invalid")
    libraries = runtime["runtime_libraries"]
    common.require(isinstance(libraries, list), "runtime library identity list is absent")
    for index, identity in enumerate(libraries):
        validate_file_identity(identity, f"runtime library {index}")
    common.hex_digest(runtime["ldd_normalized_sha256"], "normalized ldd digest")
    common.require(isinstance(runtime["environment"], dict) and
                   set(runtime["environment"]) == {
                       "LD_LIBRARY_PATH", "LD_PRELOAD", "LD_AUDIT", "OMP_NUM_THREADS",
                       "OMP_DYNAMIC", "OMP_PROC_BIND"},
                   "runtime environment identity changed")
    return runtime


def validate_live(source: dict[str, Any], build: dict[str, Any],
                  identities: dict[str, dict[str, Any]], campaign: dict[str, Any],
                  allow_self_test: bool = False) -> None:
    root = Path(source["root"])
    common.require(root.is_dir(), "live source root is absent")
    common.require(common.checked_output(["git", "-C", str(root), "rev-parse", "HEAD"]) ==
                   source["head"], "live source HEAD differs")
    common.require(common.checked_output(
        ["git", "-C", str(root), "rev-parse", "HEAD^{tree}"]) == source["tree"],
        "live source tree differs")
    common.require(not common.checked_output(
        ["git", "-C", str(root), "status", "--porcelain", "--untracked-files=all"]),
        "live source is dirty")
    if not allow_self_test:
        common.require(Path(__file__).resolve() == (root / common.ANALYZER_RELATIVE).resolve(),
                       "analyzer is not the canonical clean-tree source")
    for label in ("cache", "graph", "archive", "binary"):
        live_file(build[label], "build " + label)
    for group in ("sources", "objects"):
        for label, identity in build[group].items():
            live_file(identity, f"build {group} {label}")
    for label in ("runner", "common", "matrix", "lease_source"):
        live_file(identities[label], label)
    runtime = identities["runtime"]
    live_file(runtime["python"]["executable"], "Python executable")
    for index, identity in enumerate(runtime["runtime_libraries"]):
        live_file(identity, f"runtime library {index}")
    if not allow_self_test:
        recomputed = common.runtime_identity(
            Path(build["binary"]["path"]), campaign["cpu"], campaign["sibling"],
            set(os.sched_getaffinity(0)))
        common.require(recomputed == runtime, "live runtime identity differs")


def validate_stat(value: object, label: str, iterations: int,
                  execution: bool, rates: dict[str, int] | None = None) -> dict[str, Any]:
    suffix = "_us_per_batch_call" if execution else "_us"
    sample_key = "samples" + suffix
    keys = {"median" + suffix, "mad" + suffix, "minimum" + suffix,
            "maximum" + suffix, sample_key}
    if rates:
        keys.update(rates)
    summary = common.require_keys(value, keys, label)
    samples = summary[sample_key]
    common.require(isinstance(samples, list) and len(samples) == iterations,
                   f"{label} does not retain exactly {iterations} observations")
    numeric = [common.finite(item, label + " sample", 1e-300) for item in samples]
    median = statistics.median(numeric)
    mad = statistics.median(abs(item - median) for item in numeric)
    close_six(summary["median" + suffix], median, label + " median")
    close_six(summary["mad" + suffix], mad, label + " MAD")
    close_six(summary["minimum" + suffix], min(numeric), label + " minimum")
    close_six(summary["maximum" + suffix], max(numeric), label + " maximum")
    if rates:
        for key, byte_count in rates.items():
            close_six(summary[key], byte_count / (median * 1000.0), label + " " + key)
    return {"median": median, "mad": mad, "samples": numeric}


def validate_digests(value: object, label: str) -> dict[str, str]:
    digests = common.require_keys(value, {
        "algorithm", "original_data", "transmitted_parity", "recovered_originals",
    }, label)
    common.require(digests["algorithm"] == "fnv1a64", f"{label} algorithm changed")
    for key in ("original_data", "transmitted_parity", "recovered_originals"):
        common.hex_digest(digests[key], f"{label}.{key}", 16)
    return digests


def validate_benchmark(path: Path, case: dict[str, Any], mode: str,
                       expected_sha: str) -> dict[str, Any]:
    document, raw = common.load_json(path, "benchmark result")
    common.require(common.sha256_bytes(raw) == expected_sha,
                   f"benchmark bytes changed: {path}")
    common.require_keys(document, {
        "schema", "build", "parameters", "resolved", "correctness",
        "workload_digests", "memory", "metrics", "legacy",
    }, f"benchmark {path}")
    common.require(document["schema"] == common.BENCHMARK_SCHEMA,
                   f"{path}: benchmark schema is not v3")
    common.require_keys(document["build"], {"compiler", "compiler_version", "cplusplus"},
                        f"{path} compiler build")
    parameters = common.require_keys(document["parameters"], {
        "K", "R", "requested_profile", "requested_field", "requested_backend",
        "force_generic_decode", "force_specialized_decode", "force_tiled_decode",
        "force_materialized_decode", "skip_legacy", "retain_samples",
        "report_decode_path", "shard_bytes",
        "loss_count", "missing_original_indices", "batch", "reuse", "iterations",
        "warmup", "thread_count", "seed",
    }, f"{path} parameters")
    for key, expected in common.MODE_PARAMETERS[mode].items():
        common.require(type(parameters[key]) is bool and parameters[key] is expected,
                       f"{path}: {mode} force tuple is not exact at {key}")
    expected_parameters = {
        "K": case["K"], "R": case["R"], "requested_profile": "legacy_high_v1",
        "requested_field": "gf8", "requested_backend": case["backend"],
        "skip_legacy": True, "retain_samples": True,
        "report_decode_path": True,
        "shard_bytes": case["shard_bytes"], "loss_count": case["loss_count"],
        "batch": case["batch"], "reuse": case["reuse"],
        "iterations": case["iterations"], "warmup": case["warmup"],
        "thread_count": 1, "seed": case["seed"],
    }
    for key, expected in expected_parameters.items():
        common.require(type(parameters[key]) is type(expected) and parameters[key] == expected,
                       f"{path}: parameter {key} is not exact")
    common.require(parameters["missing_original_indices"] == list(range(case["K"])),
                   f"{path}: missing-original coordinates are not canonical full loss")
    resolved = common.require_keys(document["resolved"], {
        "profile", "field", "backend", "thread_count", "parent_count", "padded_side",
        "selected_decode_path", "selected_decode_rule",
        "decode_required_work_slots", "decode_aligned_prefix_bytes",
        "decode_tail_bytes", "decode_rounded_bytes", "decode_multi_item_batch",
    }, f"{path} resolved")
    expected_resolved = {
        "profile": "legacy_high_v1", "field": "gf8", "backend": case["backend"],
        "thread_count": 1, "parent_count": case["parent_count"],
        "padded_side": case["padded_side"],
        "selected_decode_path": mode,
        "selected_decode_rule": "forced_" + mode,
        "decode_required_work_slots": (
            case["parent_count"] if mode != "tiled" else
            2 * case["padded_side"] + case["loss_count"]),
        "decode_aligned_prefix_bytes": case["shard_bytes"] & ~63,
        "decode_tail_bytes": case["shard_bytes"] & 63,
        "decode_rounded_bytes": (case["shard_bytes"] + 63) & ~63,
        "decode_multi_item_batch": case["batch"] > 1,
    }
    for key, expected in expected_resolved.items():
        common.require(type(resolved[key]) is type(expected) and resolved[key] == expected,
                       f"{path}: resolved {key} is not exact")
    correctness = common.require_keys(document["correctness"],
                                      {"leopard2_round_trip", "legacy_comparison"},
                                      f"{path} correctness")
    common.require(correctness == {"leopard2_round_trip": True,
                                   "legacy_comparison": None},
                   f"{path}: correctness/skip-legacy contract failed")
    digests = validate_digests(document["workload_digests"], f"{path} digests")
    memory = common.require_keys(document["memory"], {
        "scratch_alignment", "encode_scratch_bytes_per_stripe",
        "decode_scratch_bytes_per_stripe", "encode_scratch_bytes_batch",
        "decode_scratch_bytes_batch",
    }, f"{path} memory")
    for key, value in memory.items():
        common.strict_int(value, f"{path} memory.{key}", 1)
    common.require(memory["encode_scratch_bytes_batch"] ==
                   memory["encode_scratch_bytes_per_stripe"] * case["batch"] and
                   memory["decode_scratch_bytes_batch"] ==
                   memory["decode_scratch_bytes_per_stripe"] * case["batch"],
                   f"{path}: batch scratch derivation failed")
    metrics = common.require_keys(document["metrics"], {
        "codec_setup", "encode_execution", "decode_plan_setup", "decode_execution",
        "decode_amortized_at_reuse", "rate_semantics",
    }, f"{path} metrics")
    encode_input = case["K"] * case["shard_bytes"] * case["batch"]
    encode_output = case["R"] * case["shard_bytes"] * case["batch"]
    decode_input = case["R"] * case["shard_bytes"] * case["batch"]
    decode_output = case["K"] * case["shard_bytes"] * case["batch"]
    validate_stat(metrics["codec_setup"], f"{path} codec setup",
                  case["iterations"], False)
    validate_stat(metrics["encode_execution"], f"{path} encode execution",
                  case["iterations"], True,
                  {"input_GB_per_s": encode_input,
                   "parity_output_GB_per_s": encode_output})
    setup = validate_stat(metrics["decode_plan_setup"], f"{path} plan setup",
                          case["iterations"], False)
    execution = validate_stat(metrics["decode_execution"], f"{path} decode execution",
                              case["iterations"], True,
                              {"offered_received_GB_per_s": decode_input,
                               "repaired_output_GB_per_s": decode_output})
    amortized = common.require_keys(metrics["decode_amortized_at_reuse"], {
        "reuse_count", "derived_median_us_per_batch_call",
        "offered_received_GB_per_s", "repaired_output_GB_per_s",
    }, f"{path} amortized")
    common.require(amortized["reuse_count"] == case["reuse"],
                   f"{path}: amortized reuse differs")
    amortized_us = execution["median"] + setup["median"] / case["reuse"]
    close_six(amortized["derived_median_us_per_batch_call"], amortized_us,
              f"{path} amortized median")
    close_six(amortized["offered_received_GB_per_s"],
              decode_input / (amortized_us * 1000.0), f"{path} amortized input")
    close_six(amortized["repaired_output_GB_per_s"],
              decode_output / (amortized_us * 1000.0), f"{path} amortized output")
    common.require(metrics["rate_semantics"] ==
                   "offered_received counts all non-null shard pointers supplied; "
                   "a plan may read a deterministic subset", f"{path}: rate semantics changed")
    legacy = common.require_keys(document["legacy"], {
        "available", "unavailable_reason", "codec_setup", "decode_timing_includes_setup",
        "encode_execution", "decode_including_setup",
    }, f"{path} legacy")
    common.require(legacy == {
        "available": False, "unavailable_reason": "disabled by --skip-legacy",
        "codec_setup": None, "decode_timing_includes_setup": True,
        "encode_execution": None, "decode_including_setup": None,
    }, f"{path}: legacy-disabled contract changed")
    return {
        "mode": mode, "decode_median_us": execution["median"],
        "decode_mad_us": execution["mad"], "decode_samples_us": execution["samples"],
        "setup_median_us": setup["median"], "setup_samples_us": setup["samples"],
        "scratch_bytes": memory["decode_scratch_bytes_per_stripe"],
        "digests": digests, "raw_sha256": expected_sha,
    }


def validate_isolation(value: object, label: str) -> dict[str, Any]:
    isolation = common.require_keys(value, {
        "policy", "accepted", "timed_before", "timed_after", "timed_delta",
        "sibling_before", "sibling_after", "sibling_delta",
    }, label)
    expected = common.isolation(
        isolation["timed_before"], isolation["timed_after"],
        isolation["sibling_before"], isolation["sibling_after"])
    common.require(isolation == expected and isolation["accepted"] is True,
                   f"{label}: CPU isolation was rejected or not recomputed")
    return isolation


def clustered(round_logs: list[float], label: str) -> dict[str, Any]:
    common.require(len(round_logs) == 3, f"{label}: exactly three rounds are required")
    mean = statistics.fmean(round_logs)
    deviation = statistics.stdev(round_logs)
    half_width = T95[2] * deviation / math.sqrt(3)
    lower = math.exp(mean - half_width)
    upper = math.exp(mean + half_width)
    ratio = math.exp(mean)
    return {
        "independent_round_count": 3, "degrees_of_freedom": 2,
        "independent_round_log_contrasts": round_logs,
        "geometric_control_over_candidate": ratio,
        "improvement_percent": 100.0 * (ratio - 1.0),
        "ci95_lower_percent": 100.0 * (lower - 1.0),
        "ci95_upper_percent": 100.0 * (upper - 1.0),
        "credible_candidate_gain_at_least_5_percent": lower >= 1.05,
        "credible_candidate_regression_over_2_percent": upper < 0.98,
    }


def validate_artifact(value: object, expected_name: str, path: Path,
                      label: str, empty: bool | None = None) -> str:
    identity = common.require_keys(value, {"name", "size", "sha256"}, label)
    common.require(identity["name"] == expected_name, f"{label} name differs")
    size = common.strict_int(identity["size"], label + " size")
    digest = common.hex_digest(identity["sha256"], label + " SHA-256")
    common.require(path.is_file() and path.stat().st_size == size and
                   common.sha256_file(path) == digest, f"{label} bytes changed")
    if empty is not None:
        common.require((size == 0) is empty, f"{label} emptiness differs")
    return digest


def analyze(manifest_path: Path, validate_live_files: bool = True,
            allow_self_test: bool = False) -> dict[str, Any]:
    manifest, manifest_raw = common.load_json(manifest_path, "balanced manifest")
    common.require_keys(manifest, {
        "schema", "status", "valid", "source", "build", "runner", "common",
        "matrix", "lease_source", "runtime", "campaign", "record_artifacts",
        "final", "content_sha256",
    }, "balanced manifest")
    validate_signed(manifest, "manifest")
    common.require(manifest["schema"] == common.RUN_SCHEMA and
                   manifest["status"] == "complete" and manifest["valid"] is True,
                   "manifest schema/status/validity differs")
    source = validate_source(manifest["source"])
    build = validate_build(manifest["build"], source["root"])
    runner = validate_file_identity(manifest["runner"], "runner", common.RUNNER_RELATIVE)
    common_id = validate_file_identity(manifest["common"], "common module",
                                       common.COMMON_RELATIVE)
    matrix_id = validate_file_identity(manifest["matrix"], "external matrix")
    lease_source = validate_file_identity(manifest["lease_source"], "lease source",
                                          common.LEASE_RELATIVE)
    runtime = validate_runtime(manifest["runtime"])
    identities = {"runner": runner, "common": common_id, "matrix": matrix_id,
                  "lease_source": lease_source, "runtime": runtime}
    final = common.require_keys(manifest["final"], {
        "source", "build", "runner", "common", "matrix", "lease_source", "runtime",
    }, "final identities")
    common.require(final == {
        "source": source, "build": build, "runner": runner, "common": common_id,
        "matrix": matrix_id, "lease_source": lease_source, "runtime": runtime,
    }, "source/object/executable/matrix/runtime identity changed during collection")
    campaign = common.require_keys(manifest["campaign"], {
        "matrix_name", "cpu", "sibling", "original_affinity", "housekeeping_affinity",
        "round_orders", "case_count", "record_count", "build_refresh_command",
        "lease", "lease_identity_sha256",
    }, "campaign")
    cpu = common.strict_int(campaign["cpu"], "campaign CPU")
    sibling = common.strict_int(campaign["sibling"], "campaign sibling")
    common.require(cpu != sibling, "campaign CPU and sibling are identical")
    original_affinity = campaign["original_affinity"]
    housekeeping = campaign["housekeeping_affinity"]
    common.require(isinstance(original_affinity, list) and
                   original_affinity == sorted(set(original_affinity)) and
                   isinstance(housekeeping, list) and housekeeping == sorted(set(housekeeping)) and
                   housekeeping == sorted(set(original_affinity) - {cpu, sibling}),
                   "campaign affinity/housekeeping identity differs")
    common.require(campaign["round_orders"] == [list(item) for item in common.ROUND_ORDERS],
                   "campaign is not ABBA/BAAB/ABBA")
    lease_digest = common.hex_digest(
        campaign["lease_identity_sha256"], "lease identity digest")
    common.require(common.canonical_sha256(campaign["lease"]) == lease_digest,
                   "lease identity digest mismatch")
    common.require(campaign["build_refresh_command"] == [
        "cmake", "--build", build["root"], "--target", "bench_leopard2", "-j", "1"],
        "build refresh command differs")
    matrix_document, _ = common.load_json(Path(matrix_id["path"]), "bound matrix")
    matrix_name, cases = common.normalize_matrix(matrix_document)
    common.require(matrix_name == campaign["matrix_name"] and
                   campaign["case_count"] == len(cases) and
                   campaign["record_count"] == len(cases) * 12,
                   "campaign/matrix cardinality differs")
    if validate_live_files:
        validate_live(source, build, identities, campaign, allow_self_test)
    record_artifacts = manifest["record_artifacts"]
    common.require(isinstance(record_artifacts, list) and
                   len(record_artifacts) == campaign["record_count"],
                   "record artifact list is incomplete")
    raw_root = manifest_path.parent / "raw"
    identity_digests = {
        "source": common.canonical_sha256(source),
        "build": common.canonical_sha256(build),
        "runner": common.canonical_sha256(runner),
        "common": common.canonical_sha256(common_id),
        "matrix": common.canonical_sha256(matrix_id),
        "lease_source": common.canonical_sha256(lease_source),
        "runtime": common.canonical_sha256(runtime),
    }
    grouped: dict[str, dict[int, dict[int, dict[str, Any]]]] = {}
    raw_hashes = []
    record_hashes = []
    sequence = 0
    for case in cases:
        for round_index, order in enumerate(common.ROUND_ORDERS):
            for slot, role in enumerate(order):
                prefix = f"{case['name']}-round{round_index}-slot{slot}-{role}"
                record_name = prefix + ".record.json"
                record_path = raw_root / record_name
                record_hash = validate_artifact(
                    record_artifacts[sequence], record_name, record_path,
                    f"record artifact {sequence}", empty=False)
                record_hashes.append(record_hash)
                envelope, _ = common.load_json(record_path, "record envelope")
                common.require_keys(envelope, {
                    "schema", "status", "sequence", "case", "round", "slot", "role",
                    "mode", "command", "command_sha256", "source_identity_sha256",
                    "build_identity_sha256", "runner_identity_sha256",
                    "common_identity_sha256", "matrix_identity_sha256",
                    "lease_source_identity_sha256", "runtime_identity_sha256",
                    "lease_identity_sha256", "start_time_ns", "end_time_ns",
                    "isolation", "artifacts", "content_sha256",
                }, f"record {sequence}")
                validate_signed(envelope, f"record {sequence}")
                mode = common.role_mode(case, role)
                common.require(envelope["schema"] == common.RECORD_SCHEMA and
                               envelope["status"] == "complete" and
                               envelope["sequence"] == sequence and
                               envelope["case"] == case["name"] and
                               envelope["round"] == round_index and
                               envelope["slot"] == slot and envelope["role"] == role and
                               envelope["mode"] == mode,
                               f"record {sequence} coordinate/forced-mode differs")
                benchmark_path = raw_root / (prefix + ".benchmark.json")
                expected_command = common.benchmark_command(
                    Path(build["binary"]["path"]), case, benchmark_path, cpu, role)
                common.require(envelope["command"] == expected_command and
                               envelope["command_sha256"] ==
                               common.canonical_sha256(expected_command),
                               f"record {sequence} exact argv differs")
                for key, expected in identity_digests.items():
                    common.require(envelope[key + "_identity_sha256"] == expected,
                                   f"record {sequence} {key} identity differs")
                common.require(envelope["lease_identity_sha256"] == lease_digest,
                               f"record {sequence} lease identity differs")
                start = common.strict_int(envelope["start_time_ns"],
                                          f"record {sequence} start", 1)
                end = common.strict_int(envelope["end_time_ns"],
                                        f"record {sequence} end", 1)
                common.require(end > start, f"record {sequence} timestamps do not increase")
                validate_isolation(envelope["isolation"], f"record {sequence} isolation")
                artifacts = common.require_keys(envelope["artifacts"],
                                                {"benchmark", "stdout", "stderr"},
                                                f"record {sequence} artifacts")
                benchmark_sha = validate_artifact(
                    artifacts["benchmark"], prefix + ".benchmark.json", benchmark_path,
                    f"record {sequence} benchmark", empty=False)
                validate_artifact(artifacts["stdout"], prefix + ".stdout",
                                  raw_root / (prefix + ".stdout"),
                                  f"record {sequence} stdout", empty=True)
                validate_artifact(artifacts["stderr"], prefix + ".stderr",
                                  raw_root / (prefix + ".stderr"),
                                  f"record {sequence} stderr", empty=True)
                observation = validate_benchmark(benchmark_path, case, mode, benchmark_sha)
                observation.update({"role": role, "round": round_index, "slot": slot,
                                    "record_sha256": record_hash})
                grouped.setdefault(case["name"], {}).setdefault(round_index, {})[slot] = observation
                raw_hashes.append(benchmark_sha)
                sequence += 1

    cells = []
    for case in cases:
        round_logs, setup_logs, retained_rounds = [], [], []
        scratch: dict[str, set[int]] = {case["control_mode"]: set(),
                                       case["candidate_mode"]: set()}
        digest_reference = None
        for round_index, order in enumerate(common.ROUND_ORDERS):
            slots = grouped[case["name"]][round_index]
            common.require(set(slots) == {0, 1, 2, 3}, "incomplete clustered round")
            observations = []
            pair_logs, pair_setup_logs = [], []
            for slot in range(4):
                item = slots[slot]
                scratch[item["mode"]].add(item["scratch_bytes"])
                if digest_reference is None:
                    digest_reference = item["digests"]
                common.require(item["digests"] == digest_reference,
                               f"{case['name']}: workload digests differ between modes")
                observations.append({
                    "slot": slot, "role": item["role"], "mode": item["mode"],
                    "decode_median_us": item["decode_median_us"],
                    "decode_mad_us": item["decode_mad_us"],
                    "decode_samples_us": item["decode_samples_us"],
                    "setup_median_us": item["setup_median_us"],
                    "setup_samples_us": item["setup_samples_us"],
                    "scratch_bytes": item["scratch_bytes"],
                    "raw_sha256": item["raw_sha256"],
                    "record_sha256": item["record_sha256"],
                })
            for left, right in ((0, 1), (2, 3)):
                pair = (slots[left], slots[right])
                control = pair[0] if pair[0]["role"] == "control" else pair[1]
                candidate = pair[1] if pair[1]["role"] == "candidate" else pair[0]
                pair_logs.append(math.log(
                    control["decode_median_us"] / candidate["decode_median_us"]))
                pair_setup_logs.append(math.log(
                    control["setup_median_us"] / candidate["setup_median_us"]))
            round_log = statistics.fmean(pair_logs)
            round_setup_log = statistics.fmean(pair_setup_logs)
            round_logs.append(round_log)
            setup_logs.append(round_setup_log)
            retained_rounds.append({
                "round": round_index, "order": list(order),
                "decode_log_contrast": round_log,
                "setup_log_contrast": round_setup_log,
                "observations": observations,
            })
        common.require(all(len(values) == 1 for values in scratch.values()),
                       f"{case['name']}: scratch changed within a forced mode")
        cell_case = {key: value for key, value in case.items()
                     if key not in {"padded_side", "parent_count"}}
        cells.append({
            **cell_case, "padded_side": case["padded_side"],
            "parent_count": case["parent_count"],
            "decode_execution": clustered(round_logs, case["name"] + " decode"),
            "decode_plan_setup": clustered(setup_logs, case["name"] + " setup"),
            "scratch_bytes_by_mode": {
                mode: next(iter(values)) for mode, values in scratch.items()},
            "workload_digests": digest_reference,
            "rounds": retained_rounds,
        })
    summary: dict[str, Any] = {
        "schema": common.SUMMARY_SCHEMA,
        "source_manifest_sha256": common.sha256_bytes(manifest_raw),
        "source": source, "build": build,
        "runner_sha256": runner["sha256"], "analyzer_sha256": common.sha256_file(Path(__file__)),
        "common_sha256": common_id["sha256"], "matrix": matrix_id,
        "lease_source_sha256": lease_source["sha256"],
        "runtime_identity_sha256": identity_digests["runtime"],
        "raw_result_set_sha256": common.canonical_sha256(raw_hashes),
        "record_set_sha256": common.canonical_sha256(record_hashes),
        "method": {
            "binary_relation": "same source/object/archive/executable closure",
            "mode_authentication": "exact benchmark argv plus exact four-boolean output tuple",
            "auto_rows_eligible": False,
            "cpu": cpu, "sibling": sibling,
            "round_orders": [list(item) for item in common.ROUND_ORDERS],
            "confidence_intervals": "clustered paired-log Student-t 95%",
            "ratio_orientation": "control_time_over_candidate_time",
            "raw_observations_retained": True,
            "strict_sibling_max_nonidle_jiffies_per_invocation": 0,
            "validity_is_independent_of_speed": True,
        },
        "cells": cells, "status": "pass",
    }
    summary["content_sha256"] = common.canonical_sha256(summary)
    return summary


def fixture_stats(samples: list[float], execution: bool,
                  rates: dict[str, int] | None = None) -> dict[str, Any]:
    suffix = "_us_per_batch_call" if execution else "_us"
    median = statistics.median(samples)
    mad = statistics.median(abs(item - median) for item in samples)
    result: dict[str, Any] = {
        "median" + suffix: median, "mad" + suffix: mad,
        "minimum" + suffix: min(samples), "maximum" + suffix: max(samples),
        "samples" + suffix: samples,
    }
    if rates:
        for key, count in rates.items():
            result[key] = count / (median * 1000.0)
    return result


def signed_write(path: Path, value: dict[str, Any]) -> None:
    document = dict(value)
    document.pop("content_sha256", None)
    document["content_sha256"] = common.canonical_sha256(document)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def make_fixture(root: Path) -> Path:
    source_root = root / "source"
    raw_root = root / "evidence/raw"
    for relative in (common.RUNNER_RELATIVE, common.COMMON_RELATIVE,
                     common.LEASE_RELATIVE, common.ANALYZER_RELATIVE,
                     "bench/leopard2/benchmark.cpp", "leopard2.cpp", "Leopard2Dispatch.h"):
        path = source_root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(relative + "\n", encoding="utf-8")
    matrix_path = source_root / "matrix.json"
    matrix = {"schema": common.MATRIX_SCHEMA, "name": "fixture", "cases": [{
        "name": "fixture", "K": 8, "R": 8, "profile": "legacy_high_v1",
        "field": "gf8", "backend": "scalar", "shard_bytes": 256,
        "loss_count": 8, "batch": 1, "reuse": 2, "iterations": 3,
        "warmup": 1, "seed": 77, "control_mode": "generic",
        "candidate_mode": "materialized",
    }]}
    matrix_path.write_text(json.dumps(matrix, sort_keys=True) + "\n", encoding="utf-8")
    build_root = source_root / "build"
    benchmark_object = build_root / "CMakeFiles/bench_leopard2.dir/bench/leopard2/benchmark.cpp.o"
    decoder_object = build_root / "CMakeFiles/leopard.dir/leopard2.cpp.o"
    for path in (benchmark_object, decoder_object):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(path.name.encode("ascii"))
    archive = build_root / "libleopard.a"
    archive.write_bytes(b"archive")
    binary = build_root / "bench_leopard2"
    binary.write_bytes(b"binary")
    binary.chmod(0o755)
    (build_root / "CMakeCache.txt").write_text(
        f"CMAKE_HOME_DIRECTORY:INTERNAL={source_root}\n", encoding="utf-8")
    (build_root / "build.ninja").write_text("fixture graph\n", encoding="utf-8")
    base_ns = 1700000000000000000
    for index, path in enumerate((benchmark_object, decoder_object, archive, binary)):
        os.utime(path, ns=(base_ns + index, base_ns + index))
    for command in (["git", "init", "-q"], ["git", "add", "."],
                    ["git", "-c", "user.name=Fixture", "-c",
                     "user.email=fixture@example.invalid", "commit", "-q", "-m", "fixture"]):
        subprocess.run(command, cwd=source_root, check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    head = common.checked_output(["git", "-C", str(source_root), "rev-parse", "HEAD"])
    source = common.git_identity(source_root, head)
    build = common.build_identity(source_root, "build/bench_leopard2")
    runner = common.file_identity(source_root / common.RUNNER_RELATIVE, source_root, "runner")
    common_id = common.file_identity(source_root / common.COMMON_RELATIVE, source_root, "common")
    matrix_id = common.file_identity(matrix_path, source_root, "matrix")
    lease_source = common.file_identity(source_root / common.LEASE_RELATIVE,
                                        source_root, "lease")
    python_id = common.file_identity(Path(sys.executable), None, "Python")
    cpuinfo_id = [{"processor": "1", "model name": "fixture"},
                  {"processor": "2", "model name": "fixture"}]
    runtime = {
        "platform": {"system": "fixture", "node": "fixture", "release": "1",
                     "version": "1", "machine": "fixture"},
        "python": {"version": sys.version, "implementation": "fixture",
                   "executable": python_id, "byteorder": sys.byteorder},
        "affinity": [0, 1, 2], "cpu": 1, "sibling": 2,
        "topology": [{"cpu": 1}, {"cpu": 2}], "clock_ticks_per_second": 100,
        "cpuinfo": cpuinfo_id, "runtime_libraries": [],
        "ldd_normalized_sha256": "1" * 64,
        "environment": {key: None for key in (
            "LD_LIBRARY_PATH", "LD_PRELOAD", "LD_AUDIT", "OMP_NUM_THREADS",
            "OMP_DYNAMIC", "OMP_PROC_BIND")},
    }
    identities = {"source": source, "build": build, "runner": runner,
                  "common": common_id, "matrix": matrix_id,
                  "lease_source": lease_source, "runtime": runtime}
    digests = {key: common.canonical_sha256(value) for key, value in identities.items()}
    lease = {"fixture": True}
    lease_digest = common.canonical_sha256(lease)
    raw_root.mkdir(parents=True)
    record_artifacts = []
    case = common.normalize_matrix(matrix)[1][0]
    sequence = 0
    for round_index, order in enumerate(common.ROUND_ORDERS):
        for slot, role in enumerate(order):
            mode = common.role_mode(case, role)
            prefix = f"fixture-round{round_index}-slot{slot}-{role}"
            benchmark_path = raw_root / (prefix + ".benchmark.json")
            stdout_path = raw_root / (prefix + ".stdout")
            stderr_path = raw_root / (prefix + ".stderr")
            control = role == "control"
            decode_samples = [20.0, 21.0, 22.0] if control else [18.0, 19.0, 20.0]
            setup_samples = [1.8, 2.0, 2.2]
            encode_input = case["K"] * case["shard_bytes"]
            encode_output = case["R"] * case["shard_bytes"]
            amortized = statistics.median(decode_samples) + 1.0
            raw = {
                "schema": common.BENCHMARK_SCHEMA,
                "build": {"compiler": "fixture", "compiler_version": "1", "cplusplus": 201103},
                "parameters": {"K": 8, "R": 8, "requested_profile": "legacy_high_v1",
                    "requested_field": "gf8", "requested_backend": "scalar",
                    **common.MODE_PARAMETERS[mode], "skip_legacy": True,
                    "retain_samples": True, "report_decode_path": True,
                    "shard_bytes": 256, "loss_count": 8,
                    "missing_original_indices": list(range(8)), "batch": 1, "reuse": 2,
                    "iterations": 3, "warmup": 1, "thread_count": 1, "seed": 77},
                "resolved": {"profile": "legacy_high_v1", "field": "gf8",
                    "backend": "scalar", "thread_count": 1, "parent_count": 16,
                    "padded_side": 8, "selected_decode_path": mode,
                    "selected_decode_rule": "forced_" + mode,
                    "decode_required_work_slots": 16 if mode != "tiled" else 24,
                    "decode_aligned_prefix_bytes": 256, "decode_tail_bytes": 0,
                    "decode_rounded_bytes": 256, "decode_multi_item_batch": False},
                "correctness": {"leopard2_round_trip": True, "legacy_comparison": None},
                "workload_digests": {"algorithm": "fnv1a64", "original_data": "1" * 16,
                    "transmitted_parity": "2" * 16, "recovered_originals": "3" * 16},
                "memory": {"scratch_alignment": 64, "encode_scratch_bytes_per_stripe": 1024,
                    "decode_scratch_bytes_per_stripe": 4096 if mode != "tiled" else 6144,
                    "encode_scratch_bytes_batch": 1024,
                    "decode_scratch_bytes_batch": 4096 if mode != "tiled" else 6144},
                "metrics": {
                    "codec_setup": fixture_stats([0.9, 1.0, 1.1], False),
                    "encode_execution": fixture_stats([9.0, 10.0, 11.0], True,
                        {"input_GB_per_s": encode_input,
                         "parity_output_GB_per_s": encode_output}),
                    "decode_plan_setup": fixture_stats(setup_samples, False),
                    "decode_execution": fixture_stats(decode_samples, True,
                        {"offered_received_GB_per_s": encode_output,
                         "repaired_output_GB_per_s": encode_input}),
                    "decode_amortized_at_reuse": {"reuse_count": 2,
                        "derived_median_us_per_batch_call": amortized,
                        "offered_received_GB_per_s": encode_output / (amortized * 1000.0),
                        "repaired_output_GB_per_s": encode_input / (amortized * 1000.0)},
                    "rate_semantics": "offered_received counts all non-null shard pointers supplied; "
                                      "a plan may read a deterministic subset"},
                "legacy": {"available": False,
                    "unavailable_reason": "disabled by --skip-legacy", "codec_setup": None,
                    "decode_timing_includes_setup": True, "encode_execution": None,
                    "decode_including_setup": None}}
            benchmark_path.write_text(json.dumps(raw, indent=2) + "\n", encoding="utf-8")
            stdout_path.write_bytes(b"")
            stderr_path.write_bytes(b"")
            command = common.benchmark_command(binary, case, benchmark_path, 1, role)
            before = {key: 100 + sequence for key in common.CPU_FIELDS}
            timed_after = dict(before); timed_after["user"] += 2
            sibling_after = dict(before); sibling_after["idle"] += 3
            record = {"schema": common.RECORD_SCHEMA, "status": "complete",
                "sequence": sequence, "case": "fixture", "round": round_index,
                "slot": slot, "role": role, "mode": mode, "command": command,
                "command_sha256": common.canonical_sha256(command),
                "source_identity_sha256": digests["source"],
                "build_identity_sha256": digests["build"],
                "runner_identity_sha256": digests["runner"],
                "common_identity_sha256": digests["common"],
                "matrix_identity_sha256": digests["matrix"],
                "lease_source_identity_sha256": digests["lease_source"],
                "runtime_identity_sha256": digests["runtime"],
                "lease_identity_sha256": lease_digest,
                "start_time_ns": 1000 + sequence * 100, "end_time_ns": 1050 + sequence * 100,
                "isolation": common.isolation(before, timed_after, before, sibling_after),
                "artifacts": {
                    "benchmark": {"name": benchmark_path.name, "size": benchmark_path.stat().st_size,
                                  "sha256": common.sha256_file(benchmark_path)},
                    "stdout": {"name": stdout_path.name, "size": 0,
                               "sha256": common.EMPTY_SHA256},
                    "stderr": {"name": stderr_path.name, "size": 0,
                               "sha256": common.EMPTY_SHA256}}}
            record_path = raw_root / (prefix + ".record.json")
            signed_write(record_path, record)
            record_artifacts.append({"name": record_path.name, "size": record_path.stat().st_size,
                                     "sha256": common.sha256_file(record_path)})
            sequence += 1
    campaign = {"matrix_name": "fixture", "cpu": 1, "sibling": 2,
        "original_affinity": [0, 1, 2], "housekeeping_affinity": [0],
        "round_orders": [list(item) for item in common.ROUND_ORDERS],
        "case_count": 1, "record_count": 12,
        "build_refresh_command": ["cmake", "--build", str(build_root), "--target",
                                  "bench_leopard2", "-j", "1"],
        "lease": lease, "lease_identity_sha256": lease_digest}
    manifest = {"schema": common.RUN_SCHEMA, "status": "complete", "valid": True,
                **identities, "campaign": campaign, "record_artifacts": record_artifacts,
                "final": copy.deepcopy(identities)}
    manifest_path = root / "evidence/manifest.json"
    signed_write(manifest_path, manifest)
    return manifest_path


def refresh_record(path: Path, manifest: dict[str, Any], index: int,
                   mutate_document) -> None:
    raw_root = path.parent / "raw"
    record_identity = manifest["record_artifacts"][index]
    record_path = raw_root / record_identity["name"]
    record, _ = common.load_json(record_path, "mutation record")
    benchmark_identity = record["artifacts"]["benchmark"]
    benchmark_path = raw_root / benchmark_identity["name"]
    document, _ = common.load_json(benchmark_path, "mutation benchmark")
    mutate_document(document)
    benchmark_path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    benchmark_identity["size"] = benchmark_path.stat().st_size
    benchmark_identity["sha256"] = common.sha256_file(benchmark_path)
    signed_write(record_path, record)
    record_identity["size"] = record_path.stat().st_size
    record_identity["sha256"] = common.sha256_file(record_path)


def self_test() -> None:
    with tempfile.TemporaryDirectory(prefix="balanced-evidence-self-test-") as temporary:
        root = Path(temporary)
        valid = make_fixture(root / "valid")
        summary = analyze(valid, validate_live_files=True, allow_self_test=True)
        common.require(len(summary["cells"]) == 1 and summary["status"] == "pass",
                       "valid fixture was rejected")

        def rejected(label: str, mutation, live: bool = False) -> None:
            path = make_fixture(root / label)
            manifest, _ = common.load_json(path, label + " manifest")
            mutation(path, manifest)
            signed_write(path, manifest)
            try:
                analyze(path, validate_live_files=live, allow_self_test=True)
            except common.EvidenceError:
                return
            raise common.EvidenceError("adversarial fixture was accepted: " + label)

        def auto_as_specialized(path, manifest):
            def mutate(document):
                for key in common.MODE_PARAMETERS["materialized"]:
                    document["parameters"][key] = False
            refresh_record(path, manifest, 1, mutate)
        rejected("auto-as-specialized", auto_as_specialized)

        def stale_source(_path, manifest):
            for target in (manifest["source"], manifest["final"]["source"]):
                target["head"] = "4" * 40
                target["expected_commit"] = "4" * 40
                target["tree"] = "5" * 40
            digest = common.canonical_sha256(manifest["source"])
            raw_root = Path(manifest["source"]["root"]).parent / "evidence/raw"
            for index, identity in enumerate(manifest["record_artifacts"]):
                record_path = raw_root / identity["name"]
                record, _ = common.load_json(record_path, "stale-source record")
                record["source_identity_sha256"] = digest
                signed_write(record_path, record)
                identity["size"] = record_path.stat().st_size
                identity["sha256"] = common.sha256_file(record_path)
        rejected("caller-asserted-stale-source", stale_source, live=True)

        def stale_binary(_path, manifest):
            for target in (manifest["build"], manifest["final"]["build"]):
                target["binary"]["sha256"] = "a" * 64
            digest = common.canonical_sha256(manifest["build"])
            raw_root = Path(manifest["source"]["root"]).parent / "evidence/raw"
            for identity in manifest["record_artifacts"]:
                record_path = raw_root / identity["name"]
                record, _ = common.load_json(record_path, "stale-binary record")
                record["build_identity_sha256"] = digest
                signed_write(record_path, record)
                identity["size"] = record_path.stat().st_size
                identity["sha256"] = common.sha256_file(record_path)
        rejected("caller-asserted-stale-binary", stale_binary, live=True)

        def sibling_work(_path, manifest):
            raw_root = Path(manifest["source"]["root"]).parent / "evidence/raw"
            identity = manifest["record_artifacts"][0]
            record_path = raw_root / identity["name"]
            record, _ = common.load_json(record_path, "sibling mutation record")
            record["isolation"]["sibling_after"]["system"] += 1
            record["isolation"] = common.isolation(
                record["isolation"]["timed_before"], record["isolation"]["timed_after"],
                record["isolation"]["sibling_before"], record["isolation"]["sibling_after"])
            signed_write(record_path, record)
            identity["size"] = record_path.stat().st_size
            identity["sha256"] = common.sha256_file(record_path)
        rejected("sibling-interference", sibling_work)

        def argv_rehash(_path, manifest):
            raw_root = Path(manifest["source"]["root"]).parent / "evidence/raw"
            identity = manifest["record_artifacts"][0]
            record_path = raw_root / identity["name"]
            record, _ = common.load_json(record_path, "argv mutation record")
            record["command"][record["command"].index("--seed") + 1] = "999"
            record["command_sha256"] = common.canonical_sha256(record["command"])
            signed_write(record_path, record)
            identity["size"] = record_path.stat().st_size
            identity["sha256"] = common.sha256_file(record_path)
        rejected("coherent-argv-rehash", argv_rehash)

        def raw_sample(path, manifest):
            refresh_record(path, manifest, 0, lambda document:
                document["metrics"]["decode_execution"][
                    "samples_us_per_batch_call"].__setitem__(0, 99.0))
        rejected("raw-sample-rehash", raw_sample)

        unsigned = make_fixture(root / "unsigned")
        manifest, _ = common.load_json(unsigned, "unsigned manifest")
        manifest["campaign"]["case_count"] = 9
        unsigned.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        try:
            analyze(unsigned, validate_live_files=False, allow_self_test=True)
        except common.EvidenceError:
            pass
        else:
            raise common.EvidenceError("unsigned manifest mutation was accepted")
    print("balanced forced-path analyzer self-test passed: valid + 7 adversarial fixtures")


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as output:
            json.dump(value, output, indent=2, sort_keys=True, allow_nan=False)
            output.write("\n")
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest")
    parser.add_argument("--output")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    if not args.manifest or not args.output:
        parser.error("--manifest and --output are required")
    summary = analyze(Path(args.manifest))
    write_json(Path(args.output), summary)
    print(json.dumps({"cells": len(summary["cells"]),
                      "content_sha256": summary["content_sha256"],
                      "schema": summary["schema"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except common.EvidenceError as error:
        print(f"analyze_balanced.py: {error}", file=sys.stderr)
        raise SystemExit(1)
