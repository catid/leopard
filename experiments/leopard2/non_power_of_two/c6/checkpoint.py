#!/usr/bin/env python3
"""Fail-closed merger for the C6 exact-GF256 experiment."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import statistics
import subprocess
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent
CORE_SHA = "48803c06fbd7a6802b4438af60e3104895938c9d"
BENCHMARK_EXECUTABLE_SHA = "11c1fd87f1c8470b19ad35de5fc33d659a91050ccda4c39c40eeb98e9f4f15aa"
BENCHMARK_LIBRARY_SHA = "69c195f4b6e6ce6a7713fed29e17f7171cb6f2e17783e8a3feb48b2fb3fcf925"
CPP_SCHEMA = "leopard2-c6-cpp/v2"
CHECKPOINT_SCHEMA = "leopard2-c6-checkpoint/v2"
EXPECTED_CORRECTNESS = {
    "cases": 12,
    "coefficients": 66920,
    "byte_comparisons": 278074,
    "digest": "0xc58c8188e359fdf2",
    "decode_cases": 320,
    "decode_byte_comparisons": 236500,
    "decode_digest": "0x29158fba4df259c1",
    "no_loss_calls": 96,
    "maximum_loss_cases": 96,
    "hot_path_allocations": 0,
}


class EvidenceError(ValueError):
    pass


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path: Path) -> Any:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise EvidenceError(f"cannot read {path}: {error}") from error
    return value


def load_algebra_module():
    path = HERE / "algebra.py"
    spec = importlib.util.spec_from_file_location("_leopard2_c6_algebra", path)
    if spec is None or spec.loader is None:
        raise EvidenceError("cannot load algebra reference")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def canonical(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=True).encode("utf-8")


def finite_number(value: Any, name: str, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise EvidenceError(f"{name} is not numeric")
    output = float(value)
    if not math.isfinite(output) or (positive and output <= 0) or output < 0:
        raise EvidenceError(f"{name} is outside its valid range")
    return output


def ceil_pow2(value: int) -> int:
    return 1 << (value - 1).bit_length()


def parent(profile: str, k: int, r: int) -> int:
    if profile == "legacy_high_v1":
        return ceil_pow2(k + ceil_pow2(r))
    if profile == "low_v1":
        return ceil_pow2(ceil_pow2(k) + r)
    raise EvidenceError("unknown benchmark profile")


def align_up(value: int, alignment: int) -> int:
    return (value + alignment - 1) & ~(alignment - 1)


def padded_scratch(profile: str, k: int, r: int, shard_bytes: int) -> int:
    rounded = align_up(shard_bytes, 64)
    side = ceil_pow2(r) if profile == "legacy_high_v1" else ceil_pow2(k)
    ranges = k + r
    pointers = k + 2 * side + (r if profile == "low_v1" else 0)
    slots = 2 * side
    pointer_offset = align_up(ranges * 16, 8)
    data_offset = align_up(pointer_offset + pointers * 8, 64)
    return data_offset + slots * rounded


def padded_decode_scratch(k: int, r: int, parent_count: int,
                          shard_bytes: int) -> int:
    rounded = align_up(shard_bytes, 64)
    ranges = 2 * k + r
    pointers = 2 * parent_count
    slots = k + r + parent_count
    pointer_offset = align_up(ranges * 16, 8)
    data_offset = align_up(pointer_offset + pointers * 8, 64)
    return data_offset + slots * rounded


def expected_specs() -> set[tuple[str, int, int, int, int, int]]:
    geometries = (
        ("legacy_high_v1", 2, 129),
        ("legacy_high_v1", 3, 129),
        ("legacy_high_v1", 3, 130),
        ("legacy_high_v1", 3, 131),
        ("legacy_high_v1", 4, 129),
        ("legacy_high_v1", 17, 129),
        ("legacy_high_v1", 17, 130),
        ("low_v1", 129, 2),
        ("low_v1", 129, 3),
        ("low_v1", 130, 3),
        ("low_v1", 131, 3),
        ("low_v1", 129, 4),
        ("low_v1", 129, 17),
        ("low_v1", 130, 17),
    )
    output = set()
    for profile, k, r in geometries:
        output.add((profile, k, r, 64, 8, 8))
        output.add((profile, k, r, 1024, 8, 4))
        output.add((profile, k, r, 65536, 1, 1))
        if k <= 3 or r <= 3:
            output.add((profile, k, r, 1048576, 1, 1))
    if len(output) != 50:
        raise AssertionError(f"internal benchmark matrix has {len(output)} cells")
    return output


def expected_decode_specs() -> set[tuple[str, int, int, int, int, int, int]]:
    patterns = (
        ("legacy_high_v1", 2, 129, 2),
        ("legacy_high_v1", 3, 129, 1),
        ("legacy_high_v1", 3, 129, 3),
        ("legacy_high_v1", 3, 130, 3),
        ("legacy_high_v1", 4, 129, 4),
        ("legacy_high_v1", 17, 129, 1),
        ("legacy_high_v1", 17, 129, 4),
        ("legacy_high_v1", 17, 129, 17),
        ("low_v1", 129, 2, 2),
        ("low_v1", 129, 3, 1),
        ("low_v1", 129, 3, 3),
        ("low_v1", 130, 3, 3),
        ("low_v1", 129, 4, 4),
        ("low_v1", 129, 17, 1),
        ("low_v1", 129, 17, 4),
        ("low_v1", 129, 17, 17),
    )
    output = set()
    for profile, k, r, losses in patterns:
        output.add((profile, k, r, losses, 64, 8, 8))
        output.add((profile, k, r, losses, 1024, 8, 4))
        output.add((profile, k, r, losses, 65536, 1, 1))
        if k <= 3 or r <= 3:
            output.add((profile, k, r, losses, 1048576, 1, 1))
    if len(output) != 56:
        raise AssertionError(f"internal decode matrix has {len(output)} cells")
    return output


def rounded_interval(value: float, places: int = 6) -> tuple[float, float]:
    error = 0.5 * 10.0 ** -places
    return value - error, value + error


def validate_derived(cell: dict[str, Any]) -> tuple[float, float]:
    exact = cell["exact_execution"]
    padded = cell["padded_execution"]
    em = finite_number(exact["median_us"], "exact median", positive=True)
    ea = finite_number(exact["mad_us"], "exact MAD")
    pm = finite_number(padded["median_us"], "padded median", positive=True)
    pa = finite_number(padded["mad_us"], "padded MAD")
    if ea >= em or pa >= pm:
        raise EvidenceError("MAD is not smaller than its positive median")
    for name in ("exact_setup", "padded_setup"):
        median = finite_number(cell[name]["median_us"], f"{name} median", positive=True)
        mad = finite_number(cell[name]["mad_us"], f"{name} MAD")
        if mad >= median:
            raise EvidenceError(f"{name} MAD is not smaller than its median")

    em_lo, em_hi = rounded_interval(em)
    ea_lo, ea_hi = rounded_interval(ea)
    pm_lo, pm_hi = rounded_interval(pm)
    pa_lo, pa_hi = rounded_interval(pa)
    ratio_lo, ratio_hi = pm_lo / em_hi, pm_hi / em_lo
    credible_lo = 100.0 * ((pm_lo - pa_hi) / (em_hi + ea_hi) - 1.0)
    credible_hi = 100.0 * ((pm_hi - max(0.0, pa_lo)) /
                           (em_lo + max(0.0, ea_lo)) - 1.0)
    ratio = finite_number(cell["padded_over_exact_ratio"], "reported ratio", positive=True)
    credible = cell["credible_gain_percent"]
    if isinstance(credible, bool) or not isinstance(credible, (int, float)) or not math.isfinite(credible):
        raise EvidenceError("reported credible gain is not finite")
    report_lo, report_hi = rounded_interval(ratio)
    if report_hi < ratio_lo or report_lo > ratio_hi:
        raise EvidenceError("reported ratio is not derived from raw medians")
    report_lo, report_hi = rounded_interval(float(credible))
    if report_hi < credible_lo or report_lo > credible_hi:
        raise EvidenceError("reported credible gain is not derived from raw timing")
    return ratio, float(credible)


def validate_cpp(path: Path, label: str, source_sha: str) -> dict[str, Any]:
    data = read_json(path)
    if not isinstance(data, dict) or data.get("schema") != CPP_SCHEMA:
        raise EvidenceError(f"{label}: wrong C++ artifact schema")
    if data.get("source_sha256") != source_sha:
        raise EvidenceError(f"{label}: artifact is not bound to current C++ source")
    if data.get("core_git_sha") != CORE_SHA:
        raise EvidenceError(f"{label}: wrong baseline core commit")
    library_sha = data.get("library_sha256")
    if not isinstance(library_sha, str) or len(library_sha) != 64 or any(
            character not in "0123456789abcdef" for character in library_sha):
        raise EvidenceError(f"{label}: malformed library SHA-256")
    if data.get("correctness") != EXPECTED_CORRECTNESS:
        raise EvidenceError(f"{label}: correctness evidence differs")
    if not isinstance(data.get("compiler"), str) or not data["compiler"]:
        raise EvidenceError(f"{label}: compiler is missing")
    if not isinstance(data.get("compiler_version"), str) or not data["compiler_version"]:
        raise EvidenceError(f"{label}: compiler version is missing")
    if not isinstance(data.get("cpu_model"), str) or not data["cpu_model"]:
        raise EvidenceError(f"{label}: CPU model is missing")
    return data


def verify_core_commit() -> None:
    result = subprocess.run(
        ["git", "cat-file", "-e", CORE_SHA + "^{commit}"], cwd=ROOT,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
    )
    if result.returncode != 0:
        raise EvidenceError("baseline core SHA is not a local commit")


def merge(paths: dict[str, Path]) -> dict[str, Any]:
    verify_core_commit()
    algebra_module = load_algebra_module()
    algebra_raw = read_json(paths["algebra"])
    algebra_expected = algebra_module.run()
    if canonical(algebra_raw) != canonical(algebra_expected):
        raise EvidenceError("algebra artifact does not reproduce byte-for-byte")

    cpp_source = HERE / "c6_gf256.cpp"
    source_sha = sha256(cpp_source)
    artifacts = {}
    for label in ("auto", "scalar", "ssse3", "avx2", "asan-ubsan", "benchmark"):
        artifacts[label] = validate_cpp(paths[label], label, source_sha)

    expected_runtime = {
        "scalar": "scalar", "ssse3": "ssse3", "avx2": "avx2",
    }
    reference_correctness = None
    for label in ("auto", "scalar", "ssse3", "avx2"):
        data = artifacts[label]
        if data.get("sanitizer_mode") != "none" or data.get("cells") != [] or \
                data.get("decode_cells") != []:
            raise EvidenceError(f"{label}: backend artifact is not correctness-only")
        if data.get("allocation_tracking") != "global-new":
            raise EvidenceError(f"{label}: hot-path allocation tracking is not enabled")
        if data.get("requested_backend") != label:
            raise EvidenceError(f"{label}: requested backend label differs")
        if label in expected_runtime and data.get("runtime_backend") != expected_runtime[label]:
            raise EvidenceError(f"{label}: forced runtime backend differs")
        if label == "auto" and data.get("runtime_backend") not in (
                "scalar", "ssse3", "avx2", "neon"):
            raise EvidenceError("auto: unknown runtime backend")
        current = canonical(data["correctness"])
        if reference_correctness is None:
            reference_correctness = current
        elif current != reference_correctness:
            raise EvidenceError("backend correctness digest is not deterministic")

    sanitizer = artifacts["asan-ubsan"]
    if sanitizer.get("sanitizer_mode") != "asan-ubsan" or \
            sanitizer.get("cells") != [] or sanitizer.get("decode_cells") != []:
        raise EvidenceError("sanitizer artifact does not identify correctness-only ASan+UBSan")
    if sanitizer.get("allocation_tracking") != "disabled-for-sanitizer":
        raise EvidenceError("sanitizer artifact does not declare disabled new interposition")

    benchmark = artifacts["benchmark"]
    if benchmark.get("requested_backend") != "avx2" or benchmark.get("runtime_backend") != "avx2":
        raise EvidenceError("benchmark did not run forced AVX2")
    if benchmark.get("sanitizer_mode") != "none":
        raise EvidenceError("benchmark is sanitizer-instrumented")
    if benchmark.get("allocation_tracking") != "global-new":
        raise EvidenceError("benchmark hot-path allocation tracking is not enabled")
    if benchmark.get("library_sha256") != BENCHMARK_LIBRARY_SHA:
        raise EvidenceError("benchmark linked archive differs from retained evidence")
    if benchmark.get("affinity") is None or len(benchmark["affinity"]) != 1:
        raise EvidenceError("benchmark was not pinned to exactly one CPU")
    if benchmark.get("omp_num_threads") != "1" or benchmark.get("openmp_max_threads") != 1:
        raise EvidenceError("benchmark was not OMP/OpenMP single-threaded")

    manifest = read_json(paths["benchmark-manifest"])
    if not isinstance(manifest, dict) or manifest.get("schema") != \
            "leopard2-c6-benchmark-run/v2" or manifest.get("status") != "pass":
        raise EvidenceError("benchmark run manifest has the wrong schema or status")
    expected_manifest = {
        "core_git_sha": CORE_SHA,
        "runner_sha256": sha256(HERE / "run_benchmark.py"),
        "source_sha256": source_sha,
        "library_sha256": benchmark["library_sha256"],
        "executable_sha256": BENCHMARK_EXECUTABLE_SHA,
        "result_sha256": sha256(paths["benchmark"]),
        "stdout_sha256": sha256(paths["benchmark-stdout"]),
        "stderr_sha256": sha256(paths["benchmark-stderr"]),
        "cpu": benchmark["affinity"][0],
        "encode_cell_count": 50,
        "decode_cell_count": 56,
        "environment": {"OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE"},
    }
    for name, value in expected_manifest.items():
        if manifest.get(name) != value:
            raise EvidenceError(f"benchmark manifest {name} differs")
    executable_sha = manifest.get("executable_sha256")
    if not isinstance(executable_sha, str) or len(executable_sha) != 64 or any(
            character not in "0123456789abcdef" for character in executable_sha):
        raise EvidenceError("benchmark manifest executable SHA-256 is malformed")
    command = manifest.get("command")
    if not isinstance(command, list) or len(command) != 4 or \
            command[1:3] != ["--backend", "avx2"]:
        raise EvidenceError("benchmark manifest command differs")

    wanted = expected_specs()
    seen = set()
    validated_cells = []
    for cell in benchmark.get("cells", []):
        if not isinstance(cell, dict):
            raise EvidenceError("benchmark cell is not an object")
        key = tuple(cell.get(name) for name in
                    ("profile", "K", "R", "bytes", "batch", "reuse"))
        if key not in wanted or key in seen:
            raise EvidenceError(f"unexpected or duplicate benchmark cell {key}")
        seen.add(key)
        profile, k, r, shard_bytes, batch, _reuse = key
        if cell.get("parent") != parent(profile, k, r) or cell["parent"] != 512:
            raise EvidenceError("benchmark cell is not a 512-parent field-boundary case")
        expected_table = k * r
        expected_terms = k * r * batch
        if cell.get("exact_table_bytes") != expected_table:
            raise EvidenceError("exact coefficient-table accounting differs")
        if cell.get("exact_execution_scratch_bytes") != 0:
            raise EvidenceError("exact executor scratch accounting differs")
        if cell.get("padded_execution_scratch_bytes") != padded_scratch(
                profile, k, r, shard_bytes):
            raise EvidenceError("padded scratch accounting differs")
        if cell.get("exact_fixed_multiply_terms") != expected_terms:
            raise EvidenceError("exact term accounting differs")
        if cell.get("exact_coefficient_payload_bytes") != expected_terms * shard_bytes:
            raise EvidenceError("exact coefficient payload accounting differs")
        if cell.get("input_bytes") != k * shard_bytes * batch or \
                cell.get("output_bytes") != r * shard_bytes * batch:
            raise EvidenceError("input/output byte accounting differs")
        for digest_name in ("exact_output_digest", "padded_output_digest"):
            digest = cell.get(digest_name)
            if not isinstance(digest, str) or len(digest) != 18 or not digest.startswith("0x"):
                raise EvidenceError("malformed output digest")
        if cell["exact_output_digest"] == cell["padded_output_digest"]:
            raise EvidenceError("benchmark lacks an output-level profile difference witness")
        ratio, credible = validate_derived(cell)
        validated_cells.append({**cell,
                                "padded_over_exact_ratio": ratio,
                                "credible_gain_percent": credible})
    if seen != wanted:
        raise EvidenceError(f"benchmark matrix is incomplete: missing {sorted(wanted - seen)}")

    plan_records = {
        (record["profile"], record["K"], record["R"], record["losses"]): record
        for record in algebra_raw.get("decoder_plan_records", [])
    }
    if len(plan_records) != 16 or algebra_raw.get("decoder") != {
        "plans": 16,
        "folded_nonzero_terms": 4933,
        "recovered_values": 70,
        "no_loss_execution_terms": 0,
    }:
        raise EvidenceError("independent decoder-plan algebra differs")
    wanted_decode = expected_decode_specs()
    wanted_plan_keys = {
        (profile, k, r, losses)
        for profile, k, r, losses, _bytes, _batch, _reuse in wanted_decode
    }
    if set(plan_records) != wanted_plan_keys:
        raise EvidenceError("independent decoder-plan geometry differs")
    for (profile, k, r, losses), record in plan_records.items():
        expected_missing = [(index * k) // losses for index in range(losses)]
        if record.get("missing_originals") != expected_missing or \
                record.get("selected_parities") != list(range(losses)):
            raise EvidenceError("independent decoder selection differs")
        digest = record.get("term_digest_sha256")
        if not isinstance(digest, str) or len(digest) != 64 or any(
                character not in "0123456789abcdef" for character in digest):
            raise EvidenceError("independent decoder term digest is malformed")
    seen_decode = set()
    validated_decode = []
    for cell in benchmark.get("decode_cells", []):
        if not isinstance(cell, dict):
            raise EvidenceError("decode benchmark cell is not an object")
        key = tuple(cell.get(name) for name in (
            "profile", "K", "R", "losses", "bytes", "batch", "reuse"
        ))
        if key not in wanted_decode or key in seen_decode:
            raise EvidenceError(f"unexpected or duplicate decode cell {key}")
        seen_decode.add(key)
        profile, k, r, losses, shard_bytes, batch, _reuse = key
        parent_count = parent(profile, k, r)
        if cell.get("parent") != parent_count or parent_count != 512:
            raise EvidenceError("decode cell is not a 512-parent boundary case")
        plan_record = plan_records.get((profile, k, r, losses))
        if plan_record is None:
            raise EvidenceError("decode cell has no independent algebra plan")
        term_count = plan_record["term_count"]
        if cell.get("selected_parity_count") != losses or \
                cell.get("selected_parity_prefix_end") != losses:
            raise EvidenceError("decode parity selection is not the declared prefix")
        if cell.get("exact_codec_table_bytes") != k * r:
            raise EvidenceError("decode codec-table accounting differs")
        if cell.get("exact_plan_terms") != term_count:
            raise EvidenceError("decode plan term count differs from independent algebra")
        if cell.get("exact_term_record_bytes") != 12 or \
                cell.get("exact_offset_record_bytes") != 8 or \
                cell.get("exact_coordinate_record_bytes") != 4:
            raise EvidenceError("decode record ABI accounting differs")
        plan_payload = term_count * 12 + (losses + 1) * 8 + losses * 2 * 4
        if cell.get("exact_plan_payload_bytes") != plan_payload:
            raise EvidenceError("decode plan payload accounting differs")
        execution_terms = term_count * batch
        if cell.get("exact_execution_terms") != execution_terms or \
                cell.get("exact_term_payload_bytes") != execution_terms * shard_bytes:
            raise EvidenceError("decode execution-term accounting differs")
        if cell.get("exact_execution_scratch_bytes") != 0:
            raise EvidenceError("exact decode scratch accounting differs")
        if cell.get("padded_execution_scratch_bytes") != padded_decode_scratch(
                k, r, parent_count, shard_bytes):
            raise EvidenceError("padded decode scratch accounting differs")
        if cell.get("offered_received_bytes") != \
                (k - losses + r) * shard_bytes * batch or \
                cell.get("repaired_output_bytes") != losses * shard_bytes * batch:
            raise EvidenceError("decode input/output byte accounting differs")
        if cell.get("baseline") != "padded_specialized_gf16":
            raise EvidenceError("decode baseline is not forced padded specialized GF16")
        for digest_name in ("exact_output_digest", "padded_output_digest"):
            digest = cell.get(digest_name)
            if not isinstance(digest, str) or len(digest) != 18 or \
                    not digest.startswith("0x"):
                raise EvidenceError("malformed decode output digest")
        if cell["exact_output_digest"] != cell["padded_output_digest"]:
            raise EvidenceError("exact and padded decoders did not restore identical originals")
        ratio, credible = validate_derived(cell)
        validated_decode.append({**cell,
                                 "padded_over_exact_ratio": ratio,
                                 "credible_gain_percent": credible})
    if seen_decode != wanted_decode:
        raise EvidenceError(
            f"decode benchmark matrix is incomplete: missing {sorted(wanted_decode - seen_decode)}"
        )

    target = [cell for cell in validated_cells if min(cell["K"], cell["R"]) == 3]
    neighbors = [cell for cell in validated_cells if min(cell["K"], cell["R"]) in (2, 4)]
    wide = [cell for cell in validated_cells if min(cell["K"], cell["R"]) == 17]
    if len(target) != 24 or len(neighbors) != 14 or len(wide) != 12:
        raise EvidenceError("benchmark region classification differs")
    target_gate = all(cell["credible_gain_percent"] >= 10.0 for cell in target)
    neighbor_gate = all(cell["credible_gain_percent"] >= -2.0 for cell in neighbors)
    measured_region = target_gate and neighbor_gate

    decode_target = [cell for cell in validated_decode
                     if min(cell["K"], cell["R"]) == 3]
    decode_neighbors = [cell for cell in validated_decode
                        if min(cell["K"], cell["R"]) in (2, 4)]
    decode_wide = [cell for cell in validated_decode
                   if min(cell["K"], cell["R"]) == 17]
    if len(decode_target) != 24 or len(decode_neighbors) != 14 or \
            len(decode_wide) != 18:
        raise EvidenceError("decode benchmark region classification differs")
    decode_target_gate = all(
        cell["credible_gain_percent"] >= 10.0 for cell in decode_target)
    decode_neighbor_gate = all(
        cell["credible_gain_percent"] >= -2.0 for cell in decode_neighbors)
    decode_measured_region = decode_target_gate and decode_neighbor_gate
    maximum_loss_cells = [cell for cell in validated_decode
                          if cell["losses"] == min(cell["K"], cell["R"])]
    if len(maximum_loss_cells) != 36:
        raise EvidenceError("decode maximum-loss coverage differs")
    decode_one_shot_ratios = [
        (cell["padded_setup"]["median_us"] +
         cell["padded_execution"]["median_us"]) /
        (cell["exact_setup"]["median_us"] +
         cell["exact_execution"]["median_us"])
        for cell in validated_decode
    ]
    decode_declared_reuse_ratios = [
        (cell["padded_setup"]["median_us"] / cell["reuse"] +
         cell["padded_execution"]["median_us"]) /
        (cell["exact_setup"]["median_us"] / cell["reuse"] +
         cell["exact_execution"]["median_us"])
        for cell in validated_decode
    ]
    decode_loss_groups = {}
    for losses in sorted({cell["losses"] for cell in validated_decode}):
        group = [cell for cell in validated_decode if cell["losses"] == losses]
        decode_loss_groups[str(losses)] = {
            "cells": len(group),
            "credible_gain_percent_min": min(
                cell["credible_gain_percent"] for cell in group),
            "credible_gain_percent_max": max(
                cell["credible_gain_percent"] for cell in group),
        }

    artifact_hashes = {label: sha256(path) for label, path in paths.items()}
    return {
        "schema": CHECKPOINT_SCHEMA,
        "status": "pass",
        "baseline_core_git_sha": CORE_SHA,
        "source_sha256": {
            "algebra.py": sha256(HERE / "algebra.py"),
            "c6_gf256.cpp": source_sha,
            "checkpoint.py": sha256(Path(__file__)),
            "run_benchmark.py": sha256(HERE / "run_benchmark.py"),
        },
        "artifacts_sha256": artifact_hashes,
        "correctness": {
            "gf4_exhaustive": algebra_raw["gf4_exhaustive"],
            "gf8_algebra": algebra_raw["gf8"],
            "gf8_decoder_algebra": algebra_raw["decoder"],
            "cpp_per_backend": EXPECTED_CORRECTNESS,
            "backends": {label: artifacts[label]["runtime_backend"]
                         for label in ("auto", "scalar", "ssse3", "avx2")},
            "asan_ubsan": "pass",
        },
        "geometry": {
            "affected_profile_cells": algebra_raw["affected_profile_cells"],
            "benchmark_cells": len(validated_cells),
            "decode_benchmark_cells": len(validated_decode),
            "benchmark_parent": 512,
        },
        "benchmark": {
            "cpu": benchmark["affinity"][0],
            "cpu_model": benchmark["cpu_model"],
            "compiler": benchmark["compiler"],
            "compiler_version": benchmark["compiler_version"],
            "runtime_backend": benchmark["runtime_backend"],
            "executable_sha256": executable_sha,
            "linked_library_sha256": benchmark["library_sha256"],
            "target_cells": len(target),
            "neighbor_cells": len(neighbors),
            "wide_cells": len(wide),
            "target_credible_gain_percent_min": min(
                cell["credible_gain_percent"] for cell in target),
            "target_credible_gain_percent_median": statistics.median(
                cell["credible_gain_percent"] for cell in target),
            "neighbor_credible_gain_percent_min": min(
                cell["credible_gain_percent"] for cell in neighbors),
            "wide_credible_gain_percent_min": min(
                cell["credible_gain_percent"] for cell in wide),
            "wide_credible_gain_percent_max": max(
                cell["credible_gain_percent"] for cell in wide),
            "target_at_least_10_percent": target_gate,
            "neighbors_no_regression_over_2_percent": neighbor_gate,
            "meaningful_region_gate": measured_region,
        },
        "decode_benchmark": {
            "baseline": "padded_specialized_gf16",
            "target_cells": len(decode_target),
            "neighbor_cells": len(decode_neighbors),
            "wide_cells": len(decode_wide),
            "maximum_loss_cells": len(maximum_loss_cells),
            "maximum_loss_credible_gain_percent_min": min(
                cell["credible_gain_percent"] for cell in maximum_loss_cells),
            "maximum_loss_credible_gain_percent_max": max(
                cell["credible_gain_percent"] for cell in maximum_loss_cells),
            "target_credible_gain_percent_min": min(
                cell["credible_gain_percent"] for cell in decode_target),
            "target_credible_gain_percent_median": statistics.median(
                cell["credible_gain_percent"] for cell in decode_target),
            "neighbor_credible_gain_percent_min": min(
                cell["credible_gain_percent"] for cell in decode_neighbors),
            "wide_credible_gain_percent_min": min(
                cell["credible_gain_percent"] for cell in decode_wide),
            "wide_credible_gain_percent_max": max(
                cell["credible_gain_percent"] for cell in decode_wide),
            "target_at_least_10_percent": decode_target_gate,
            "neighbors_no_regression_over_2_percent": decode_neighbor_gate,
            "meaningful_region_gate": decode_measured_region,
            "exact_setup_median_us_min": min(
                cell["exact_setup"]["median_us"] for cell in validated_decode),
            "exact_setup_median_us_max": max(
                cell["exact_setup"]["median_us"] for cell in validated_decode),
            "padded_setup_median_us_min": min(
                cell["padded_setup"]["median_us"] for cell in validated_decode),
            "padded_setup_median_us_max": max(
                cell["padded_setup"]["median_us"] for cell in validated_decode),
            "one_shot_setup_plus_execution_ratio_min": min(
                decode_one_shot_ratios),
            "one_shot_setup_plus_execution_ratio_median": statistics.median(
                decode_one_shot_ratios),
            "declared_reuse_setup_plus_execution_ratio_min": min(
                decode_declared_reuse_ratios),
            "declared_reuse_setup_plus_execution_ratio_median": statistics.median(
                decode_declared_reuse_ratios),
            "loss_count_groups": decode_loss_groups,
        },
        "memory": {
            "exact_execution_scratch_bytes": 0,
            "exact_table_bytes_min": min(cell["exact_table_bytes"] for cell in validated_cells),
            "exact_table_bytes_max": max(cell["exact_table_bytes"] for cell in validated_cells),
            "padded_execution_scratch_bytes_min": min(
                cell["padded_execution_scratch_bytes"] for cell in validated_cells),
            "padded_execution_scratch_bytes_max": max(
                cell["padded_execution_scratch_bytes"] for cell in validated_cells),
            "exact_decode_execution_scratch_bytes": 0,
            "exact_decode_plan_payload_bytes_min": min(
                cell["exact_plan_payload_bytes"] for cell in validated_decode),
            "exact_decode_plan_payload_bytes_max": max(
                cell["exact_plan_payload_bytes"] for cell in validated_decode),
            "padded_decode_execution_scratch_bytes_min": min(
                cell["padded_execution_scratch_bytes"] for cell in validated_decode),
            "padded_decode_execution_scratch_bytes_max": max(
                cell["padded_execution_scratch_bytes"] for cell in validated_decode),
        },
        "method_scope": algebra_raw["method_scope"],
        "disposition": {
            "measured_candidate_region": measured_region,
            "measured_decode_region": decode_measured_region,
            "production_promotion": "none",
            "reason": (
                "the exact GF8 prefix code has no frozen serialized profile identity, "
                "second-host timing, production API review, or calibrated dispatcher"
            ),
            "legacy_wire_compatible": False,
            "default_build_changed": False,
            "bead_recommendation": "close-experiment-result",
            "promotion_followup": "C7/C8/C10 and W",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    for name in (
        "algebra", "auto", "scalar", "ssse3", "avx2", "asan-ubsan",
        "benchmark", "benchmark-manifest", "benchmark-stdout", "benchmark-stderr",
    ):
        parser.add_argument("--" + name, dest=name.replace("-", "_"), type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    paths = {name: getattr(args, name.replace("-", "_")) for name in
             (
                 "algebra", "auto", "scalar", "ssse3", "avx2", "asan-ubsan",
                 "benchmark", "benchmark-manifest", "benchmark-stdout",
                 "benchmark-stderr",
             )}
    result = merge(paths)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except EvidenceError as error:
        print(f"c6 checkpoint: {error}", file=sys.stderr)
        raise SystemExit(1)
