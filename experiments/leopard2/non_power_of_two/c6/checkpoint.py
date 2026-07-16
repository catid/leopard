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
BENCHMARK_EXECUTABLE_SHA = "d4fe26a82e15099b4d5d096a186bfb0c03348146e802b1d84af3c5ea0f7049f0"
BENCHMARK_LIBRARY_SHA = "69c195f4b6e6ce6a7713fed29e17f7171cb6f2e17783e8a3feb48b2fb3fcf925"
CPP_SCHEMA = "leopard2-c6-cpp/v1"
CHECKPOINT_SCHEMA = "leopard2-c6-checkpoint/v1"
EXPECTED_CORRECTNESS = {
    "cases": 12,
    "coefficients": 66920,
    "byte_comparisons": 278074,
    "digest": "0xc58c8188e359fdf2",
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
        if data.get("sanitizer_mode") != "none" or data.get("cells") != []:
            raise EvidenceError(f"{label}: backend artifact is not correctness-only")
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
    if sanitizer.get("sanitizer_mode") != "asan-ubsan" or sanitizer.get("cells") != []:
        raise EvidenceError("sanitizer artifact does not identify correctness-only ASan+UBSan")

    benchmark = artifacts["benchmark"]
    if benchmark.get("requested_backend") != "avx2" or benchmark.get("runtime_backend") != "avx2":
        raise EvidenceError("benchmark did not run forced AVX2")
    if benchmark.get("sanitizer_mode") != "none":
        raise EvidenceError("benchmark is sanitizer-instrumented")
    if benchmark.get("library_sha256") != BENCHMARK_LIBRARY_SHA:
        raise EvidenceError("benchmark linked archive differs from retained evidence")
    if benchmark.get("affinity") is None or len(benchmark["affinity"]) != 1:
        raise EvidenceError("benchmark was not pinned to exactly one CPU")
    if benchmark.get("omp_num_threads") != "1" or benchmark.get("openmp_max_threads") != 1:
        raise EvidenceError("benchmark was not OMP/OpenMP single-threaded")

    manifest = read_json(paths["benchmark-manifest"])
    if not isinstance(manifest, dict) or manifest.get("schema") != \
            "leopard2-c6-benchmark-run/v1" or manifest.get("status") != "pass":
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
        "cell_count": 50,
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

    target = [cell for cell in validated_cells if min(cell["K"], cell["R"]) == 3]
    neighbors = [cell for cell in validated_cells if min(cell["K"], cell["R"]) in (2, 4)]
    wide = [cell for cell in validated_cells if min(cell["K"], cell["R"]) == 17]
    if len(target) != 24 or len(neighbors) != 14 or len(wide) != 12:
        raise EvidenceError("benchmark region classification differs")
    target_gate = all(cell["credible_gain_percent"] >= 10.0 for cell in target)
    neighbor_gate = all(cell["credible_gain_percent"] >= -2.0 for cell in neighbors)
    measured_region = target_gate and neighbor_gate

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
            "cpp_per_backend": EXPECTED_CORRECTNESS,
            "backends": {label: artifacts[label]["runtime_backend"]
                         for label in ("auto", "scalar", "ssse3", "avx2")},
            "asan_ubsan": "pass",
        },
        "geometry": {
            "affected_profile_cells": algebra_raw["affected_profile_cells"],
            "benchmark_cells": len(validated_cells),
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
        "memory": {
            "exact_execution_scratch_bytes": 0,
            "exact_table_bytes_min": min(cell["exact_table_bytes"] for cell in validated_cells),
            "exact_table_bytes_max": max(cell["exact_table_bytes"] for cell in validated_cells),
            "padded_execution_scratch_bytes_min": min(
                cell["padded_execution_scratch_bytes"] for cell in validated_cells),
            "padded_execution_scratch_bytes_max": max(
                cell["padded_execution_scratch_bytes"] for cell in validated_cells),
        },
        "method_scope": algebra_raw["method_scope"],
        "disposition": {
            "measured_candidate_region": measured_region,
            "production_promotion": "none",
            "reason": (
                "the exact GF8 prefix code has no frozen serialized profile identity, "
                "production decoder, second-host timing, or calibrated dispatcher"
            ),
            "legacy_wire_compatible": False,
            "default_build_changed": False,
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
