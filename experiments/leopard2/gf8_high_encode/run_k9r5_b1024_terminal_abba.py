#!/usr/bin/env python3
"""Qualify the ordinary GF8 K=9/R=5/1024-byte AVX2 terminal.

The target is compared with both an exact-main binary and the same frozen
Leopard2 executable with only the K=9/R=5/1024-byte terminal disabled.  Byte
and shape boundaries prove that the selector is inert, while the previously
promoted K=9/R=5/256-byte terminal remains an exact-main reference.
"""

from __future__ import annotations

import argparse
import importlib.util
import math
import re
import statistics
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence


BASE_PATH = Path(__file__).resolve().with_name(
    "run_k5r5_b64_terminal_abba.py")
SOURCE_ROOT = Path(__file__).resolve().parents[3]
PROVENANCE_PATH = SOURCE_ROOT / "tools/leopard2_build_provenance.py"


def load_base() -> Any:
    specification = importlib.util.spec_from_file_location(
        "k9r5_b1024_terminal_evidence_base", BASE_PATH)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load packed-terminal support: {BASE_PATH}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


BASE = load_base()


def load_build_provenance() -> Any:
    specification = importlib.util.spec_from_file_location(
        "k9r5_b1024_terminal_build_provenance", PROVENANCE_PATH)
    if specification is None or specification.loader is None:
        raise RuntimeError(
            f"cannot load production build provenance: {PROVENANCE_PATH}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


BUILD_PROVENANCE = load_build_provenance()
BASE.__doc__ = __doc__
BASE.SCHEMA = "leopard2-gf8-k9r5-b1024-terminal-abba/v2"
BASE.SUMMARY_SCHEMA = "leopard2-gf8-k9r5-b1024-terminal-summary/v2"
BASE.MODE_SYMBOL = \
    "_ZN12_GLOBAL__N_1L26g_k9r5_b1024_terminal_modeE"
BASE.ALLOW_IDENTICAL_CANDIDATE_CONTROL = True
BASE.CONTROL_EXTRA_ARGUMENTS = (
    "--disable-k9r5-b1024-terminal",
)
BASE.CONTROL_BUILD_MARKER = \
    "k9r5_b1024_terminal_diagnostic_disabled"
BASE.REQUIRE_EXPECTED_IDENTITIES = True
BASE.REQUIRE_BUILD_CLOSURE = True
BASE.REQUIRE_FULL_ELF_IDENTITY = True
BASE.TARGET_CONTROL_FLOOR = 1.05
BASE.TARGET_MAIN_FLOOR = 1.05
BASE.NEIGHBOR_EQUIVALENCE_LOWER = 1.0 / 1.02
BASE.NEIGHBOR_EQUIVALENCE_UPPER = 1.02
BASE.NEIGHBOR_ORDER = (
    ("control", "candidate", "candidate", "control"),
    ("candidate", "control", "control", "candidate"),
    ("control", "candidate", "candidate", "control"),
    ("candidate", "control", "control", "candidate"),
)
BASE.RUNNER_PATH = Path(__file__).resolve()
BASE.RUNNER_DEPENDENCIES = (
    BASE.RUNNER_PATH,
    BASE_PATH,
    Path(BASE.T8_SUPPORT.__file__).resolve(),
    Path(BASE.MAIN_SUPPORT.__file__).resolve(),
    PROVENANCE_PATH.resolve(),
)
BASE.RETAINED_MAIN_FLOOR = 1.0
BASE.CANONICAL_MAIN_SHA256 = \
    "78575e2fe1b9796d10dcb4c14d2bcc1b627a79e787be66d34dc7aa666899aa93"
CONFIRMATORY_ROUNDS = 25
CONFIDENCE_LEVEL = 0.95
T95_DF24 = 2.0638985616280205


_BASE_CONFIDENCE_INTERVAL = BASE.confidence_interval


def confidence_interval(values: Sequence[float]) -> dict[str, Any]:
    """Use a predeclared 25-round interval for this qualification only."""
    if len(values) != CONFIRMATORY_ROUNDS:
        return _BASE_CONFIDENCE_INTERVAL(values)
    BASE.require(all(isinstance(value, (int, float)) for value in values),
                 "contrast is not numeric")
    center = statistics.mean(values)
    half_width = T95_DF24 * statistics.stdev(values) / math.sqrt(len(values))
    return {
        "speedup": math.exp(center),
        "ci95": [
            math.exp(center - half_width),
            math.exp(center + half_width),
        ],
        "round_log_ratios": list(values),
        "confidence_level": CONFIDENCE_LEVEL,
        "degrees_of_freedom": CONFIRMATORY_ROUNDS - 1,
        "t_critical": T95_DF24,
    }


BASE.confidence_interval = confidence_interval


_BASE_SELECT_ROUND_ORDERS = BASE.select_round_orders


def select_round_orders(
    orders: Sequence[Sequence[str]], requested_rounds: int | None,
) -> tuple[tuple[str, ...], ...]:
    """Add one qualification-local 25-round confirmatory schedule."""
    if requested_rounds != CONFIRMATORY_ROUNDS:
        return _BASE_SELECT_ROUND_ORDERS(orders, requested_rounds)
    BASE.require(len(orders) > 0, "round-order cycle is empty")
    return tuple(
        tuple(orders[index % len(orders)])
        for index in range(CONFIRMATORY_ROUNDS)
    )


BASE.select_round_orders = select_round_orders


def parse_arguments() -> argparse.Namespace:
    """Expose only inputs that the shared-binary provenance gate consumes."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--candidate-sha256")
    parser.add_argument("--control", required=True, type=Path)
    parser.add_argument("--control-sha256")
    parser.add_argument("--main", required=True, type=Path)
    parser.add_argument(
        "--main-sha256", default=BASE.CANONICAL_MAIN_SHA256)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--source-tree", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cpu", required=True, type=int)
    parser.add_argument("--sibling", required=True, type=int)
    parser.add_argument("--iterations", type=int, default=31)
    parser.add_argument("--warmup", type=int, default=64)
    parser.add_argument(
        "--rounds", type=int, choices=(CONFIRMATORY_ROUNDS,),
        default=CONFIRMATORY_ROUNDS)
    # The imported runner still consults these generic compatibility fields.
    # This wrapper deliberately has no archive/compile-command CLI escape
    # hatch: the production validator discovers the entire graph from the
    # exact CMake build root.
    parser.set_defaults(
        candidate_archive=None,
        candidate_archive_sha256=None,
        control_archive=None,
        control_archive_sha256=None,
        candidate_compile_commands=None,
        candidate_compile_commands_sha256=None,
        control_compile_commands=None,
        control_compile_commands_sha256=None,
    )
    options = parser.parse_args()
    BASE.require(
        options.main_sha256 == BASE.CANONICAL_MAIN_SHA256,
        "exact-main executable SHA-256 differs from the canonical AVX2 "
        "Leopard-main benchmark")
    return options


BASE.parse_arguments = parse_arguments


def _require_identity(
    value: Any, expected_path: Path, label: str,
) -> Mapping[str, Any]:
    BASE.require(isinstance(value, Mapping), f"{label} identity is malformed")
    observed_path = Path(str(value.get("path", ""))).resolve(strict=True)
    BASE.require(observed_path == expected_path.resolve(strict=True),
                 f"{label} identity names another path")
    BASE.require(isinstance(value.get("sha256"), str) and
                 re.fullmatch(r"[0-9a-f]{64}", value["sha256"]) is not None
                 and type(value.get("size")) is int and value["size"] > 0,
                 f"{label} content identity is malformed")
    return value


def _source_dependency_identities(
    provenance: Mapping[str, Any],
) -> list[dict[str, Any]]:
    manifest = provenance.get("tracked_source_manifest")
    BASE.require(isinstance(manifest, Mapping) and
                 isinstance(manifest.get("files"), list),
                 "production provenance lacks its tracked source manifest")
    records = {
        record.get("path"): record
        for record in manifest["files"]
        if isinstance(record, Mapping) and isinstance(record.get("path"), str)
    }
    BASE.require(len(records) == len(manifest["files"]),
                 "tracked source manifest is malformed or duplicated")
    result = []
    for dependency in BASE.RUNNER_DEPENDENCIES:
        resolved = dependency.resolve(strict=True)
        BASE.require(resolved.is_relative_to(SOURCE_ROOT),
                     "runner dependency escapes the candidate source root")
        relative = resolved.relative_to(SOURCE_ROOT).as_posix()
        record = records.get(relative)
        identity = BASE.support_file_identity(resolved)
        BASE.require(isinstance(record, Mapping) and
                     record.get("sha256") == identity["sha256"] and
                     record.get("size") == identity["size"],
                     f"runner dependency is absent from tracked source: "
                     f"{relative}")
        result.append({
            "path": relative,
            "sha256": identity["sha256"],
            "size": identity["size"],
            "mode": record.get("mode"),
        })
    BASE.require(len(result) == len(set(item["path"] for item in result)),
                 "runner dependency closure contains duplicates")
    return result


def _validate_production_build_contract(
    provenance: Any, build: Path, executable: Path,
    options: argparse.Namespace,
) -> tuple[Mapping[str, Any], list[dict[str, Any]]]:
    BASE.require(isinstance(provenance, Mapping),
                 "production build provenance is malformed")
    BASE.require(
        provenance.get("schema") ==
            BUILD_PROVENANCE.PRODUCTION_BUILD_CLOSURE_SCHEMA and
        Path(str(provenance.get("build_root", ""))).resolve(strict=True) ==
            build and
        Path(str(provenance.get("physical_source_root", ""))).resolve(
            strict=True) == SOURCE_ROOT and
        Path(str(provenance.get("source_root", ""))).resolve(strict=True) ==
            SOURCE_ROOT and
        provenance.get("executable_target") == "bench_leopard2",
        "production provenance names another schema, source, build, or target")

    manifest = provenance.get("tracked_source_manifest")
    git = manifest.get("git") if isinstance(manifest, Mapping) else None
    BASE.require(isinstance(git, Mapping) and
                 git.get("commit") == options.source_commit and
                 git.get("tree") == options.source_tree and
                 git.get("dirty") is False,
                 "production source is dirty or differs from the requested "
                 "commit/tree")

    cache = provenance.get("validated_cache")
    BASE.require(isinstance(cache, Mapping) and
                 cache.get(
                     "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA") ==
                     BUILD_PROVENANCE.BENCHMARK_BUILD_CONFIGURATION_SCHEMA and
                 cache.get("LEO2_BUILD_BENCHMARKS") == "ON" and
                 cache.get("LEOPARD_ENABLE_GF8") == "ON" and
                 cache.get("LEO2_BACKEND_VARIANT") in {"auto", "avx2"} and
                 cache.get("LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR") ==
                     "OFF" and
                 cache.get("LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING") ==
                     "ON" and
                 isinstance(cache.get(
                     "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256"), str)
                 and re.fullmatch(
                     r"[0-9a-f]{64}", cache[
                         "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256"])
                 is not None,
                 "production CMake selector/schema contract differs")

    expected_paths = {
        "cmake_cache": build / "CMakeCache.txt",
        "compile_commands": build / "compile_commands.json",
        "archive_link_recipe": build / "CMakeFiles/leopard.dir/link.txt",
        "executable_link_recipe":
            build / "CMakeFiles/bench_leopard2.dir/link.txt",
        "archive": build / "libleopard.a",
        "executable": executable,
    }
    for key, path in expected_paths.items():
        _require_identity(provenance.get(key), path, key.replace("_", " "))

    closure = provenance.get("source_object_compile_closure")
    BASE.require(isinstance(closure, list) and closure,
                 "production compile closure is empty")
    archive_records = []
    benchmark_records = []
    for record in closure:
        BASE.require(isinstance(record, Mapping) and
                     record.get("role") in {"archive", "benchmark"} and
                     isinstance(record.get("source"), Mapping) and
                     isinstance(record.get("object"), Mapping) and
                     isinstance(record.get("compile_entry"), Mapping) and
                     isinstance(record["compile_entry"].get("arguments"),
                                list) and
                     bool(record["compile_entry"]["arguments"]) and
                     all(isinstance(argument, str) and argument for argument
                         in record["compile_entry"]["arguments"]) and
                     isinstance(record.get("flag_profile"), str) and
                     bool(record["flag_profile"]),
                     "production compile closure record is malformed")
        source = Path(str(record["source"].get("path", ""))).resolve(
            strict=True)
        obj = Path(str(record["object"].get("path", ""))).resolve(
            strict=True)
        BASE.require(source.is_relative_to(SOURCE_ROOT) and
                     obj.is_relative_to(build),
                     "production compile source/object escapes its root")
        (archive_records if record["role"] == "archive" else
         benchmark_records).append(record)
    BASE.require(archive_records and len(benchmark_records) == 1 and
                 Path(str(benchmark_records[0]["source"]["path"])).resolve(
                     strict=True) ==
                 (SOURCE_ROOT / "bench/leopard2/benchmark.cpp").resolve(
                     strict=True),
                 "production compile closure has the wrong target sources")

    members = provenance.get("archive_members")
    member_identities = provenance.get("archive_member_identities")
    BASE.require(isinstance(members, list) and members and
                 isinstance(member_identities, list) and
                 len(member_identities) == len(members) and
                 len(members) == len(archive_records),
                 "production archive member graph is malformed")
    object_by_name = {
        Path(str(record["object"]["path"])).name: record["object"]
        for record in archive_records
    }
    BASE.require(len(object_by_name) == len(archive_records) and
                 set(object_by_name) == set(members),
                 "production archive members differ from compile objects")
    for member in member_identities:
        BASE.require(isinstance(member, Mapping) and
                     member.get("member") in object_by_name and
                     member.get("sha256") ==
                         object_by_name[member["member"]].get("sha256") and
                     member.get("size") ==
                         object_by_name[member["member"]].get("size"),
                     "production archive member bytes differ from objects")

    archive_commands = provenance.get("archive_link_commands")
    executable_command = provenance.get("executable_link_command")
    BASE.require(isinstance(archive_commands, list) and
                 len(archive_commands) == 2 and
                 all(isinstance(command, list) and command and
                     all(isinstance(argument, str) and argument
                         for argument in command)
                     for command in archive_commands) and
                 isinstance(executable_command, list) and
                 executable_command and
                 all(isinstance(argument, str) and argument
                     for argument in executable_command),
                 "production archive/executable link graph is malformed")
    executable_identity = provenance["executable"]
    input_identity = BASE.T8_SUPPORT.regular_file_identity(executable)
    expected_candidate_sha256 = options.candidate_sha256
    expected_control_sha256 = options.control_sha256
    BASE.require(
        isinstance(expected_candidate_sha256, str) and
        re.fullmatch(r"[0-9a-f]{64}",
                     expected_candidate_sha256) is not None and
        expected_control_sha256 == expected_candidate_sha256,
        "candidate/control required executable hashes differ or are malformed")
    BASE.require(
        executable_identity.get("sha256") == input_identity["sha256"] and
        executable_identity.get("size") == input_identity["size"] and
        executable_identity.get("sha256") == expected_candidate_sha256,
        "production executable differs from the frozen benchmark input "
        "identity")
    return provenance, _source_dependency_identities(provenance)


_BUILD_CLOSURE_BASELINE: dict[str, Any] | None = None


def shared_build_closure_identity(options: Any) -> dict[str, Any]:
    """Bind both runtime labels to one reproducible production build."""
    global _BUILD_CLOSURE_BASELINE
    try:
        candidate = options.candidate.resolve(strict=True)
        control = options.control.resolve(strict=True)
        BASE.require(candidate == control,
                     "candidate/control must name the same production binary")
        BASE.require(candidate.name == "bench_leopard2",
                     "candidate is not the exact bench_leopard2 CMake target")
        build = candidate.parent.resolve(strict=True)
        captured = BUILD_PROVENANCE.candidate_build_provenance(
            build, SOURCE_ROOT, candidate, "bench_leopard2")
        captured, dependencies = _validate_production_build_contract(
            captured, build, candidate, options)
        if _BUILD_CLOSURE_BASELINE is None:
            replay = BUILD_PROVENANCE.verify_reproducible_candidate_build(
                captured, jobs=1)
            BUILD_PROVENANCE.validate_reproducible_build_proof(
                replay, captured, label="K9/R5/B1024 candidate")
            _BUILD_CLOSURE_BASELINE = {
                "candidate_build": captured,
                "reproducible_build": replay,
                "runner_source_dependencies": dependencies,
            }
        else:
            BASE.require(
                captured == _BUILD_CLOSURE_BASELINE["candidate_build"] and
                dependencies ==
                    _BUILD_CLOSURE_BASELINE["runner_source_dependencies"],
                "production build or runner source closure changed during "
                "the campaign")
            BUILD_PROVENANCE.validate_reproducible_build_proof(
                _BUILD_CLOSURE_BASELINE["reproducible_build"], captured,
                label="retained K9/R5/B1024 candidate")
        return _BUILD_CLOSURE_BASELINE
    except BASE.EvidenceError:
        raise
    except Exception as error:
        raise BASE.EvidenceError(
            f"production build provenance rejected the candidate: {error}") \
            from error


BASE.build_closure_identity = shared_build_closure_identity


def cells() -> list[dict[str, Any]]:
    values = (
        ("target-k9-r5-b1024-q1", 9, 5, 1024, "target", True),
        ("byte-control-k9-r5-b1023-q1", 9, 5, 1023, "neighbor", False),
        ("byte-control-k9-r5-b1025-q1", 9, 5, 1025, "neighbor", False),
        ("shape-control-k9-r6-b1024-q1", 9, 6, 1024, "neighbor", False),
        ("shape-control-k8-r5-b1024-q1", 8, 5, 1024, "neighbor", False),
        ("shape-control-k10-r5-b1024-q1", 10, 5, 1024, "neighbor", False),
        ("retained-k9-r5-b256-q1", 9, 5, 256, "neighbor", True),
    )
    result = []
    for index, (name, k, r, shard_bytes, role, compare_main) in enumerate(
            values):
        result.append({
            "id": name,
            "K": k,
            "R": r,
            "bytes": shard_bytes,
            "loss": 1,
            "batch": 1,
            "reuse": 8192,
            "role": role,
            "compare_main": compare_main,
            "measure_one_shot": True,
            "seed": 0x09541000 + index,
        })
    BASE.require(
        len(result) == 7 and
        sum(cell["role"] == "target" for cell in result) == 1 and
        len({cell["id"] for cell in result}) == len(result) and
        len({cell["seed"] for cell in result}) == len(result) and
        {cell["id"] for cell in result if cell["compare_main"]} == {
            "target-k9-r5-b1024-q1",
            "retained-k9-r5-b256-q1",
        } and
        all(cell["loss"] == 1 and cell["batch"] == 1 and
            cell["reuse"] == 8192 and cell["measure_one_shot"] is True
            for cell in result),
        "K9/R5/B1024 qualification matrix is incomplete")
    return result


BASE.cells = cells


_BASE_ANALYZE = BASE.analyze


def analyze(
    cell: Mapping[str, Any], rounds: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Prove selector inertness and retain the B256 exact-main reference."""
    result = _BASE_ANALYZE(cell, rounds)
    if cell.get("role") == "neighbor":
        contrasts = {
            "encode_execution": result.get("control_over_candidate"),
            "one_shot_encode":
                result.get("one_shot_control_over_candidate"),
        }
        confidence_intervals = {}
        for name, ratio in contrasts.items():
            BASE.require(isinstance(ratio, Mapping) and
                         isinstance(ratio.get("ci95"), list) and
                         len(ratio["ci95"]) == 2,
                         f"neighbor {name} selector contrast is missing")
            lower, upper = ratio["ci95"]
            BASE.require(
                isinstance(lower, (int, float)) and
                isinstance(upper, (int, float)) and
                lower >= BASE.NEIGHBOR_EQUIVALENCE_LOWER and
                upper <= BASE.NEIGHBOR_EQUIVALENCE_UPPER,
                f"neighbor {name} does not prove selector equivalence inside "
                f"[{BASE.NEIGHBOR_EQUIVALENCE_LOWER}, "
                f"{BASE.NEIGHBOR_EQUIVALENCE_UPPER}]: "
                f"CI [{lower}, {upper}]")
            confidence_intervals[name] = [lower, upper]
        result["neighbor_selector_inertness_validation"] = {
            "equivalence_band": [
                BASE.NEIGHBOR_EQUIVALENCE_LOWER,
                BASE.NEIGHBOR_EQUIVALENCE_UPPER,
            ],
            "ci95": confidence_intervals,
            "accepted": True,
        }
    if cell.get("id") != "retained-k9-r5-b256-q1":
        return result
    ratios = {
        "encode_execution": result.get("main_over_candidate"),
        "one_shot_encode": result.get("one_shot_main_over_candidate"),
    }
    lower_bounds = {}
    for name, ratio in ratios.items():
        BASE.require(isinstance(ratio, Mapping) and
                     isinstance(ratio.get("ci95"), list) and
                     len(ratio["ci95"]) == 2,
                     f"retained B256 {name} exact-main contrast is missing")
        lower = ratio["ci95"][0]
        BASE.require(isinstance(lower, (int, float)) and
                     lower >= BASE.RETAINED_MAIN_FLOOR,
                     f"retained B256 {name} regressed against exact main: "
                     f"lower CI {lower!r} < {BASE.RETAINED_MAIN_FLOOR}")
        lower_bounds[name] = lower
    result["retained_main_validation"] = {
        "floor": BASE.RETAINED_MAIN_FLOOR,
        "lower_ci95": lower_bounds,
        "accepted": True,
    }
    return result


BASE.analyze = analyze


if __name__ == "__main__":
    raise SystemExit(BASE.main())
