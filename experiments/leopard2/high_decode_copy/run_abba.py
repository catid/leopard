#!/usr/bin/env python3
"""Collect fail-closed same-binary Algorithm 5 copy/no-copy ABBA evidence."""

from __future__ import annotations

import argparse
import copy
import datetime as dt
import importlib.util
import json
import math
import os
import re
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


LEGACY_RAW_SCHEMA = "leopard2-high-decode-copy-raw/v2"
LEGACY_MANIFEST_SCHEMA = "leopard2-high-decode-copy-manifest/v2"
RAW_SCHEMA = "leopard2-high-decode-copy-raw/v3"
MANIFEST_SCHEMA = "leopard2-high-decode-copy-manifest/v3"
# Failure bundles are diagnostics rather than replay inputs, so only the
# current failure schema is emitted and no historical failure schema is
# advertised as verifiable.
FAILURE_SCHEMA = "leopard2-high-decode-copy-failure/v3"
RAW_SCHEMAS = frozenset({LEGACY_RAW_SCHEMA, RAW_SCHEMA})
MANIFEST_RAW_SCHEMAS = {
    LEGACY_MANIFEST_SCHEMA: LEGACY_RAW_SCHEMA,
    MANIFEST_SCHEMA: RAW_SCHEMA,
}
MATRIX_SCHEMA = "leopard2-high-decode-copy-matrix/v1"
RUNNER_RELATIVE = "experiments/leopard2/high_decode_copy/run_abba.py"
CONTRACT_RELATIVE = "experiments/leopard2/high_decode_copy/benchmark_contract.py"
MATRIX_RELATIVE = "experiments/leopard2/high_decode_copy/matrix.json"
TARGET_NAME = "bench_leopard2_high_decode_copy_attribution"
MASK64 = (1 << 64) - 1
RECORD_KEYS = {
    "cell", "mode", "round", "slot", "command", "command_sha256",
    "started_monotonic_ns", "ended_monotonic_ns",
    "cpu_before", "cpu_after", "cpu_delta",
    "sibling_before", "sibling_after", "sibling_delta", "result",
}
HEX40 = re.compile(r"^[0-9a-f]{40}$")
HEX64 = re.compile(r"^[0-9a-f]{64}$")
HOOK_ARCHIVE_DEFINITIONS = frozenset({
    "LEO2_DISABLE_AVX2_CODEGEN=1",
    "LEO2_DISABLE_SSSE3_CODEGEN=1",
    "LEO2_ENABLE_TEST_HOOKS=1",
    "LEO2_HAVE_AVX2_BACKEND=1",
    "LEO2_HAVE_AVX512_BACKEND=1",
    "LEO2_HAVE_SSSE3_BACKEND=1",
})


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


def parse_canonical_utc(value: object, label: str) -> dt.datetime:
    require(type(value) is str and value.endswith("Z"),
            f"{label} is not a canonical UTC timestamp")
    try:
        parsed = dt.datetime.fromisoformat(value[:-1] + "+00:00")
    except ValueError as error:
        raise EvidenceError(f"{label} is not an RFC3339 timestamp: {error}") from error
    require(parsed.tzinfo is not None and
            parsed.utcoffset() == dt.timedelta(0) and
            parsed.isoformat().replace("+00:00", "Z") == value,
            f"{label} is not canonical timezone-aware UTC")
    return parsed


def validate_manifest_raw_timestamps(
    manifest_created: object, raw_created: object
) -> None:
    manifest_time = parse_canonical_utc(manifest_created, "manifest created_utc")
    raw_time = parse_canonical_utc(raw_created, "raw created_utc")
    require(manifest_time >= raw_time,
            "manifest creation time predates its raw evidence")


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


class XorShift64:
    """Exact unsigned generator used by bench/leopard2/benchmark.cpp."""

    def __init__(self, seed: int):
        self.state = seed if seed else 0x9E3779B97F4A7C15

    def next(self) -> int:
        value = self.state
        value ^= (value << 13) & MASK64
        value ^= value >> 7
        value ^= (value << 17) & MASK64
        self.state = value & MASK64
        return self.state


def expected_missing_original_indices(cell: Mapping[str, Any]) -> list[int]:
    """Reproduce SelectLosses without trusting benchmark-emitted coordinates."""
    k = cell.get("k")
    losses = cell.get("losses")
    seed = cell.get("seed")
    require(all(type(value) is int for value in (k, losses, seed)) and
            0 < k and 0 < losses <= k and 0 <= seed <= MASK64,
            "loss-set inputs are not exact bounded integers")
    order = list(range(k))
    random = XorShift64(seed ^ 0xD1B54A32D192ED03)
    for remaining in range(len(order), 1, -1):
        selected = random.next() % remaining
        order[remaining - 1], order[selected] = \
            order[selected], order[remaining - 1]
    return sorted(order[:losses])


def validate_missing_original_indices(
    value: object, cell: Mapping[str, Any]
) -> list[int]:
    expected = expected_missing_original_indices(cell)
    require(type(value) is list and len(value) == cell["losses"] and
            all(type(index) is int and 0 <= index < cell["k"]
                for index in value) and
            value == sorted(set(value)) and value == expected,
            "benchmark loss set is not the deterministic matrix loss set")
    return value


def expected_record_sequence(
    matrix: Mapping[str, Any]
) -> list[tuple[str, int, int, str]]:
    return [
        (cell["id"], round_index, slot, mode)
        for cell in matrix["cells"]
        for round_index, order in enumerate(matrix["round_orders"])
        for slot, mode in enumerate(order)
    ]


def validate_record_sequence(
    records: object, matrix: Mapping[str, Any]
) -> list[Mapping[str, Any]]:
    expected = expected_record_sequence(matrix)
    require(type(records) is list and len(records) == len(expected) and
            all(isinstance(record, dict) for record in records),
            "raw campaign does not contain its exact record sequence")
    for record, coordinate in zip(records, expected):
        require(type(record.get("cell")) is str and
                type(record.get("round")) is int and
                type(record.get("slot")) is int and
                type(record.get("mode")) is str and
                (record["cell"], record["round"], record["slot"],
                 record["mode"]) == coordinate,
                "record order, round, slot, or ABBA role was relabeled")
    return records


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


def diagnostic_specifications(
    source_root: Path, build: Path, evidence_schema: str = RAW_SCHEMA,
    include_gfni: bool = False,
) -> list[tuple[Path, Path, set[str]]]:
    require(type(evidence_schema) is str and evidence_schema in RAW_SCHEMAS,
            "unknown high-decode copy evidence schema")
    ordinary = (
        "leopard.cpp", "leopard2.cpp", "Leopard2Backend.cpp",
        "Leopard2BackendScalar.cpp", "Leopard2CpuFeatures.cpp",
        "Leopard2Plan.cpp", "LeopardCommon.cpp", "LeopardFF8.cpp",
        "LeopardFF16.cpp",
    )
    backends = [
        ("Leopard2BackendSSSE3.cpp", "leopard2_backend_ssse3_test_hooks"),
        ("Leopard2BackendAVX2.cpp", "leopard2_backend_avx2_test_hooks"),
    ]
    if evidence_schema == RAW_SCHEMA:
        # V3 binds the split AVX2 XOR translation unit in its exact CMake
        # archive position.  V2 predates that split and remains replayable
        # against its historical twelve-member hooks archive.
        backends.append((
            "Leopard2BackendAVX2Xor.cpp",
            "leopard2_backend_avx2_test_hooks"))
    backends.append(
        ("Leopard2BackendAVX512.cpp", "leopard2_backend_avx512_test_hooks"),
    )
    if include_gfni:
        backends.append(
            ("Leopard2BackendGFNI.cpp", "leopard2_backend_gfni_test_hooks"))
    specifications: list[tuple[Path, Path, set[str]]] = [
        (source_root / name,
         build / "CMakeFiles/leopard_test_hooks.dir" / f"{name}.o",
         {"LEO2_ENABLE_TEST_HOOKS=1"})
        for name in ordinary
    ]
    specifications.extend(
        (source_root / name,
         build / "CMakeFiles" / f"{target}.dir" / f"{name}.o",
         {"LEO2_ENABLE_TEST_HOOKS=1"})
        for name, target in backends
    )
    specifications.append((
        source_root / "bench/leopard2/benchmark.cpp",
        build / "CMakeFiles" / f"{TARGET_NAME}.dir" /
            "bench/leopard2/benchmark.cpp.o",
        {"LEO2_ENABLE_TEST_HOOKS=1", "LEO2_HIGH_DECODE_COPY_ATTRIBUTION=1"}))
    return specifications


def validate_diagnostic_compile_flags(
    tokens: Sequence[str], label: str, source_root: Path, output: Path,
    required_definitions: set[str], include_gfni: bool = False,
) -> None:
    """Validate the experiment-specific preprocessor/include closure.

    The shared exact-main validator deliberately rejects arbitrary ``-D`` and
    ``-I`` options.  This diagnostic build has a small, source-attested set of
    private hook definitions and one source-root include, so validate those
    exactly before applying the shared compiler-flag policy to the remainder.
    """
    private_definitions = [
        token[2:] for token in tokens
        if token.startswith("-D") and token != "-DNDEBUG"
    ]
    expected_definitions = set(required_definitions)
    if output.parent.name == "leopard_test_hooks.dir":
        expected_definitions = set(HOOK_ARCHIVE_DEFINITIONS)
        if include_gfni:
            expected_definitions.add("LEO2_HAVE_GFNI_BACKEND=1")
    elif output.parent.name == "leopard2_backend_avx512_test_hooks.dir":
        expected_definitions.add("LEO2_HAVE_AVX2_BACKEND=1")
    elif output.parent.name == "leopard2_backend_gfni_test_hooks.dir":
        expected_definitions.update({
            "LEO2_HAVE_AVX2_BACKEND=1",
            "LEO2_HAVE_GFNI_BACKEND=1",
        })
    require(len(private_definitions) == len(set(private_definitions)) and
            set(private_definitions) == expected_definitions,
            f"{label} private definitions differ")

    include_flags = [token for token in tokens if token.startswith("-I")]
    require(include_flags == [f"-I{source_root}"],
            f"{label} include path differs")
    shared_tokens = [
        "-fopenmp" if token == "-fopenmp=libomp" else token
        for token in tokens
        if not (token.startswith("-D") and token != "-DNDEBUG") and
           not token.startswith("-I") and token != "-Werror" and
           not (output.parent.name ==
                "leopard2_backend_gfni_test_hooks.dir" and
                token == "-mgfni")
    ]
    if output.parent.name == "leopard2_backend_gfni_test_hooks.dir":
        require(tokens.count("-mgfni") == 1,
                f"{label} GFNI ISA flag differs")
    SUPPORT.validate_effective_flags(shared_tokens, label, "compile")


def compile_entries_include_gfni(
    entries: Sequence[Mapping[str, Any]], build: Path,
) -> bool:
    expected = (
        build / "CMakeFiles/leopard2_backend_gfni_test_hooks.dir" /
        "Leopard2BackendGFNI.cpp.o").resolve()
    matches = 0
    for entry in entries:
        output_value = entry.get("output")
        directory_value = entry.get("directory")
        if not isinstance(output_value, str) or not isinstance(
                directory_value, str):
            continue
        output = Path(output_value)
        if not output.is_absolute():
            output = Path(directory_value) / output
        if output.resolve() == expected:
            matches += 1
    require(matches <= 1,
            "diagnostic compile commands duplicate the optional GFNI object")
    return matches == 1


def diagnostic_compile_closure(
    source_root: Path, build: Path, commands_path: Path, compiler: Path,
    compiler_invocation: str, compiler_identity: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], Path, list[Path], str]:
    try:
        commands_text = commands_path.read_text(encoding="utf-8", errors="strict")
        entries = json.loads(commands_text)
    except json.JSONDecodeError as error:
        raise EvidenceError(f"invalid diagnostic compile_commands.json: {error}") from error
    require(isinstance(entries, list) and entries,
            "diagnostic compile_commands.json is empty")
    include_gfni = compile_entries_include_gfni(entries, build)
    specifications = diagnostic_specifications(
        source_root, build, include_gfni=include_gfni)
    benchmark_object = specifications[-1][1]
    expected_outputs = {output.resolve(): (source.resolve(), definitions)
                        for source, output, definitions in specifications}
    selected: dict[Path, tuple[Mapping[str, Any], list[str]]] = {}
    for entry in entries:
        require(isinstance(entry, dict), "diagnostic compile entry is not an object")
        output_value = entry.get("output")
        directory_value = entry.get("directory")
        if not isinstance(output_value, str) or not isinstance(directory_value, str):
            continue
        candidate = Path(output_value)
        candidate = candidate if candidate.is_absolute() \
            else Path(directory_value) / candidate
        if candidate.resolve() not in expected_outputs:
            continue
        tokens = SUPPORT.compile_command_tokens(entry)
        source_value = entry.get("file")
        require(isinstance(source_value, str),
                "diagnostic compile entry has no source")
        output = SUPPORT.validate_compile_entry_io(
            entry, tokens, Path(source_value).resolve(strict=True), build,
            compiler, compiler_invocation)
        require(output not in selected,
                f"duplicate diagnostic compile action for {output}")
        selected[output] = (entry, tokens)
    require(set(selected) == set(expected_outputs),
            "diagnostic compile commands omit a required hooks object")

    records: list[dict[str, Any]] = []
    for output, (expected_source, definitions) in expected_outputs.items():
        entry, tokens = selected[output]
        source_value = entry.get("file")
        directory_value = entry.get("directory")
        require(isinstance(directory_value, str) and
                Path(directory_value).resolve(strict=True) == build and
                isinstance(source_value, str) and
                Path(source_value).resolve(strict=True) == expected_source and
                tokens and tokens[0] == compiler_invocation and
                Path(tokens[0]).resolve(strict=True) == compiler and
                tokens.count("-c") == 1 and
                tokens.index("-c") + 1 < len(tokens) and
                Path(tokens[tokens.index("-c") + 1]).resolve(strict=True) ==
                    expected_source and
                not any(token.startswith("@") for token in tokens),
                f"diagnostic compile action source/compiler differs for {output}")
        validate_diagnostic_compile_flags(
            tokens, f"diagnostic compile action {output}", source_root,
            output, definitions, include_gfni)
        require("-fopenmp" in tokens or "-fopenmp=libomp" in tokens,
                f"diagnostic compile action lacks OpenMP: {output}")
        define_tokens = {token[2:] for token in tokens if token.startswith("-D")}
        require(definitions.issubset(define_tokens),
                f"diagnostic compile action lacks private gates: {output}")
        source_identity = SUPPORT.artifact_identity(expected_source, "source_file")
        object_identity = SUPPORT.artifact_identity(output, "object_file")
        require(object_identity["mtime_ns"] >= source_identity["mtime_ns"],
                f"diagnostic object predates source: {output}")
        records.append({
            "source": source_identity, "object": object_identity,
            "directory": str(build), "command_tokens": list(tokens),
            "compiler_invocation": compiler_invocation,
            "compiler_resolved_path": str(compiler),
            "compiler_invocation_identity": copy.deepcopy(compiler_identity),
            "required_definitions": sorted(definitions),
        })
    archive_objects = [output.resolve() for _, output, _ in specifications[:-1]]
    return (records, benchmark_object.resolve(strict=True), archive_objects,
            commands_text)


def validate_archive_closure(
    build: Path, archive: Path, recipe_path: Path, archive_objects: Sequence[Path],
    archiver: Path, ranlib: Path, archiver_invocation: str,
    ranlib_invocation: str,
) -> dict[str, Any]:
    recipe = recipe_path.read_text(encoding="utf-8", errors="strict")
    commands = [shlex.split(line) for line in recipe.splitlines() if line.strip()]
    require(len(commands) == 2,
            "hooks archive recipe must contain exactly ar and ranlib commands")
    archive_tokens, ranlib_tokens = commands
    require(len(archive_tokens) == 3 + len(archive_objects) and
            archive_tokens[0] == archiver_invocation and
            Path(archive_tokens[0]).resolve(strict=True) == archiver and
            archive_tokens[1] in ("qc", "rc", "rcs"),
            "hooks archive archiver recipe changed")
    archive_output = Path(archive_tokens[2])
    archive_output = archive_output if archive_output.is_absolute() \
        else build / archive_output
    recipe_objects = []
    for token in archive_tokens[3:]:
        require(not token.startswith(("@", "-")) and token.endswith(".o"),
                "hooks archive contains an indirect or non-object input")
        item = Path(token)
        recipe_objects.append((item if item.is_absolute() else build / item).resolve(strict=True))
    require(archive_output.resolve(strict=True) == archive and
            recipe_objects == list(archive_objects) and
            len(set(recipe_objects)) == len(recipe_objects),
            "hooks archive object order/closure differs from compile commands")
    require(len(ranlib_tokens) == 2 and
            ranlib_tokens[0] == ranlib_invocation and
            Path(ranlib_tokens[0]).resolve(strict=True) == ranlib and
            ranlib_tokens[1] == archive_tokens[2],
            "hooks archive ranlib recipe changed")
    members = checked_output([str(archiver), "t", str(archive)]).splitlines()
    expected_members = [path.name for path in recipe_objects]
    require(len(set(expected_members)) == len(expected_members) and
            members == expected_members,
            "hooks archive member list/order differs from its object closure")
    member_digests = []
    for member, object_path in zip(members, recipe_objects):
        member_bytes = SUPPORT.run_checked(
            [str(archiver), "p", str(archive), member],
            environment=SUPPORT.CHILD_ENVIRONMENT, timeout=30,
            max_stdout=16 * 1024 * 1024, max_stderr=1024 * 1024)
        require(SUPPORT.sha256_bytes(member_bytes) == SUPPORT.sha256_file(object_path),
                f"hooks archive member differs from retained object: {member}")
        member_digests.append({"member": member,
                               "sha256": SUPPORT.sha256_bytes(member_bytes)})
    archive_identity = SUPPORT.artifact_identity(archive, "archive")
    require(all(archive_identity["mtime_ns"] >=
                SUPPORT.artifact_identity(path, "object_file")["mtime_ns"]
                for path in recipe_objects),
            "hooks archive predates one of its retained objects")
    recipe_identity = SUPPORT.artifact_identity(recipe_path, "build_metadata")
    require(recipe_identity["sha256"] ==
            SUPPORT.sha256_bytes(recipe.encode("utf-8")),
            "hooks archive recipe changed while it was validated")
    return {
        "recipe": recipe_identity,
        "recipe_text": recipe,
        "recipe_tokens": [archive_tokens, ranlib_tokens],
        "archiver_invocation": archiver_invocation,
        "ranlib_invocation": ranlib_invocation,
        "archive": archive_identity,
        "members": member_digests,
        "archiver": SUPPORT.artifact_identity(archiver, "archiver"),
        "ranlib": SUPPORT.artifact_identity(ranlib, "ranlib"),
        "archiver_invocation_identity": SUPPORT.artifact_identity(
            Path(archiver_invocation), "archiver"),
        "ranlib_invocation_identity": SUPPORT.artifact_identity(
            Path(ranlib_invocation), "ranlib"),
    }


def validate_executable_link_inputs(
    tokens: Sequence[str], build: Path, source_root: Path, compiler: Path,
    compiler_invocation: str, compiler_identity: Mapping[str, Any],
    binary: Path, benchmark_object: Path, hook_archive: Path
) -> dict[str, Any]:
    require(tokens and tokens[0] == compiler_invocation and
            Path(tokens[0]).resolve(strict=True) == compiler and
            tokens.count("-o") == 1,
            "diagnostic executable link compiler/output shape changed")
    output_index = tokens.index("-o")
    require(output_index + 1 < len(tokens), "diagnostic link output is absent")
    output = Path(tokens[output_index + 1])
    output = output if output.is_absolute() else build / output
    require(output.resolve(strict=True) == binary,
            "diagnostic link recipe does not produce the declared executable")
    explicit_inputs: list[tuple[str, Path]] = []
    for index, token in enumerate(tokens[1:], start=1):
        if index == output_index or index == output_index + 1:
            continue
        require(not token.startswith("@"),
                "diagnostic link recipe uses an undeclared response file")
        if token.startswith("-"):
            require(not token.startswith(("-L", "-l", "-Xlinker")) and
                    (not token.startswith("-Wl,") or
                     token.startswith("-Wl,-rpath,")),
                    "diagnostic link recipe contains an undeclared linker input flag")
            if token.startswith("-Wl,-rpath,"):
                rpath = Path(token.split(",", 2)[2]).resolve(strict=True)
                require(rpath.is_relative_to(Path("/usr/lib")) or
                        rpath.is_relative_to(Path("/lib")),
                        "diagnostic link rpath is outside system library roots")
            continue
        item = Path(token)
        explicit_inputs.append((
            token, (item if item.is_absolute() else build / item).resolve(strict=True)))
    build_inputs = [(token, path) for token, path in explicit_inputs
                    if path.is_relative_to(build)]
    require([path for _token, path in build_inputs] ==
                [benchmark_object, hook_archive],
            "diagnostic executable has an extra/substituted build-tree input")
    require(not any(path.is_relative_to(source_root) and
                    not path.is_relative_to(build)
                    for _token, path in explicit_inputs),
            "diagnostic executable links a source-tree artifact")
    system_inputs = [(token, path) for token, path in explicit_inputs
                     if (token, path) not in build_inputs]
    allowed_system_names = re.compile(r"lib(?:gomp|omp|pthread)\.(?:a|so(?:\..*)?)$")
    require(all((path.is_relative_to(Path("/usr/lib")) or
                 path.is_relative_to(Path("/lib"))) and
                allowed_system_names.fullmatch(path.name) is not None
                for _token, path in system_inputs),
            "diagnostic executable has an undeclared external link input")
    return {
        "recipe_tokens": list(tokens),
        "compiler_invocation": compiler_invocation,
        "compiler_resolved_path": str(compiler),
        "compiler_invocation_identity": copy.deepcopy(compiler_identity),
        "output_token": tokens[output_index + 1],
        "output_path": str(binary),
        "build_inputs": [
            {"token": token, "resolved_path": str(path)}
            for token, path in build_inputs],
        "system_inputs": [
            {"token": token, "resolved_path": str(path)}
            for token, path in system_inputs],
    }


def validate_clean_rebuild(
    source_root: Path, cache: Mapping[str, str],
    binary_identity: Mapping[str, Any], archive_identity: Mapping[str, Any],
    executable_recipe: str, archive_recipe: str, cmake_invocation: str,
    cmake: Path,
) -> dict[str, Any]:
    require(Path(cmake_invocation).resolve(strict=True) == cmake,
            "diagnostic CMake invocation resolves to another executable")
    generator = cache.get("CMAKE_GENERATOR")
    require(generator == "Unix Makefiles",
            "diagnostic build identity requires the Unix Makefiles generator")
    retained_options = {
        "CMAKE_BUILD_TYPE": "Release",
        "CMAKE_AR": cache.get("CMAKE_AR", ""),
        "CMAKE_CXX_COMPILER": cache.get("CMAKE_CXX_COMPILER", ""),
        "CMAKE_CXX_FLAGS": cache.get("CMAKE_CXX_FLAGS", ""),
        "CMAKE_CXX_FLAGS_RELEASE": cache.get("CMAKE_CXX_FLAGS_RELEASE", ""),
        "CMAKE_EXPORT_COMPILE_COMMANDS": "ON",
        "CMAKE_RANLIB": cache.get("CMAKE_RANLIB", ""),
        "ENABLE_OPENMP": "ON",
        "LEO2_BACKEND_VARIANT": cache.get("LEO2_BACKEND_VARIANT", "auto"),
        "LEO2_BUILD_BENCHMARKS": "ON", "LEO2_BUILD_FUZZERS": "OFF",
        "LEO2_BUILD_TESTS": "ON", "LEO2_ENABLE_CUDA": "OFF",
        "LEOPARD_ENABLE_GF8": "ON", "LEOPARD_ENABLE_GF16": "ON",
    }
    require(all(retained_options[name]
                for name in ("CMAKE_AR", "CMAKE_CXX_COMPILER", "CMAKE_RANLIB")),
            "diagnostic build omitted a compiler/archive tool")
    with tempfile.TemporaryDirectory(prefix="leo2-high-copy-clean-build-") as root_text:
        clean = Path(root_text)
        configure = [str(cmake), "-S", str(source_root), "-B", str(clean),
                     "-G", generator]
        platform_value = cache.get("CMAKE_GENERATOR_PLATFORM", "")
        toolset_value = cache.get("CMAKE_GENERATOR_TOOLSET", "")
        if platform_value:
            configure.extend(("-A", platform_value))
        if toolset_value:
            configure.extend(("-T", toolset_value))
        configure.extend(f"-D{name}={value}"
                         for name, value in retained_options.items())
        SUPPORT.run_checked(
            configure, environment=SUPPORT.CHILD_ENVIRONMENT, timeout=120,
            max_stdout=4 * 1024 * 1024, max_stderr=4 * 1024 * 1024)
        SUPPORT.run_checked(
            [str(cmake), "--build", str(clean), "--target", TARGET_NAME,
             "--parallel", str(min(128, os.cpu_count() or 1))],
            environment=SUPPORT.CHILD_ENVIRONMENT, timeout=300,
            max_stdout=16 * 1024 * 1024, max_stderr=4 * 1024 * 1024)
        clean_binary = clean / TARGET_NAME
        clean_archive = clean / "libleopard_test_hooks.a"
        require(SUPPORT.sha256_file(clean_binary) == binary_identity["sha256"] and
                SUPPORT.sha256_file(clean_archive) == archive_identity["sha256"],
                "diagnostic archive/executable differ from a clean source rebuild")
        clean_executable_recipe = (clean / "CMakeFiles" / f"{TARGET_NAME}.dir" /
                                   "link.txt").read_text(
                                       encoding="utf-8", errors="strict")
        clean_archive_recipe = (clean / "CMakeFiles/leopard_test_hooks.dir/link.txt") \
            .read_text(encoding="utf-8", errors="strict")
        def normalized_recipe(text: str, root: Path) -> list[list[str]]:
            commands = []
            for line in text.splitlines():
                if not line.strip():
                    continue
                normalized = []
                for token in shlex.split(line):
                    if token.startswith("-Wl,-rpath,"):
                        path = Path(token.split(",", 2)[2]).resolve(strict=True)
                        normalized.append(f"-Wl,-rpath,{path}")
                        continue
                    if token.startswith("-") or token in ("qc", "rc", "rcs"):
                        normalized.append(token)
                        continue
                    item = Path(token)
                    candidate = item if item.is_absolute() else root / item
                    if not candidate.exists():
                        normalized.append(token)
                        continue
                    resolved = candidate.resolve(strict=True)
                    if resolved.is_relative_to(root):
                        normalized.append(
                            "BUILD/" + resolved.relative_to(root).as_posix())
                    else:
                        normalized.append(str(resolved))
                commands.append(normalized)
            return commands

        current_root = Path(str(binary_identity["path"])).parent
        require(normalized_recipe(clean_executable_recipe, clean) ==
                    normalized_recipe(executable_recipe, current_root) and
                normalized_recipe(clean_archive_recipe, clean) ==
                    normalized_recipe(archive_recipe, current_root),
                "diagnostic link/archive recipes differ from a clean source rebuild")
    return {
        "cmake": SUPPORT.artifact_identity(cmake, "build_tool"),
        "cmake_invocation": cmake_invocation,
        "cmake_resolved_path": str(cmake),
        "cmake_invocation_identity": SUPPORT.artifact_identity(
            Path(cmake_invocation), "build_tool"),
        "generator": generator,
        "retained_options": retained_options,
        "archive_sha256": archive_identity["sha256"],
        "binary_sha256": binary_identity["sha256"],
    }


def build_identity(source_root: Path, binary: Path, hook_archive: Path) -> dict[str, Any]:
    source_root = source_root.resolve(strict=True)
    binary = binary.resolve(strict=True)
    hook_archive = hook_archive.resolve(strict=True)
    require(os.access(binary, os.X_OK), "diagnostic benchmark is not executable")
    build = find_build_root(binary)
    require(binary == (build / TARGET_NAME).resolve(strict=True) and
            hook_archive == (build / "libleopard_test_hooks.a").resolve(strict=True),
            "diagnostic executable/archive are not the canonical target outputs")
    cache_path = build / "CMakeCache.txt"
    cache_bytes = cache_path.read_bytes()
    cache = SUPPORT.parse_cmake_cache(cache_path)
    cache_identity = SUPPORT.artifact_identity(cache_path, "build_metadata")
    require(cache_identity["sha256"] == SUPPORT.sha256_bytes(cache_bytes),
            "diagnostic CMake cache changed while it was validated")
    global_flags = shlex.split(cache.get("CMAKE_CXX_FLAGS", ""))
    require(cache.get("CMAKE_GENERATOR") == "Unix Makefiles",
            "diagnostic build identity requires the Unix Makefiles generator")
    require(Path(cache.get("CMAKE_HOME_DIRECTORY", "")).resolve() == source_root and
            cache.get("CMAKE_BUILD_TYPE") == "Release" and
            cache.get("CMAKE_CXX_FLAGS_RELEASE") == "-O3 -DNDEBUG" and
            all(flag in ("-Wall", "-Wextra", "-Werror")
                for flag in global_flags) and
            cache.get("ENABLE_OPENMP") == "ON" and
            cache.get("LEO2_BACKEND_VARIANT") == "auto" and
            cache.get("LEO2_BUILD_TESTS") == "ON" and
            cache.get("LEO2_BUILD_BENCHMARKS") == "ON" and
            cache.get("LEO2_BUILD_FUZZERS") == "OFF" and
            cache.get("LEO2_ENABLE_CUDA") == "OFF" and
            cache.get("LEOPARD_ENABLE_GF8") == "ON" and
            cache.get("LEOPARD_ENABLE_GF16") == "ON",
            "diagnostic CMake cache is not a Release tests/hooks build")
    commands_path = build / "compile_commands.json"
    require(commands_path.is_file(), "diagnostic build omitted compile_commands.json")
    compiler_invocation = cache.get("CMAKE_CXX_COMPILER", "")
    archiver_invocation = cache.get("CMAKE_AR", "")
    ranlib_invocation = cache.get("CMAKE_RANLIB", "")
    cmake_invocation = cache.get("CMAKE_COMMAND", "")
    require(all(type(item) is str and item for item in (
                compiler_invocation, archiver_invocation, ranlib_invocation,
                cmake_invocation)),
            "diagnostic CMake cache omitted a build-tool invocation")
    compiler = Path(compiler_invocation).resolve(strict=True)
    archiver = Path(archiver_invocation).resolve(strict=True)
    ranlib = Path(ranlib_invocation).resolve(strict=True)
    cmake = Path(cmake_invocation).resolve(strict=True)
    system_tool_roots = (Path("/usr/bin"), Path("/usr/lib"), Path("/lib"))
    require(all(any(tool.is_relative_to(root) for root in system_tool_roots)
                for tool in (compiler, archiver, ranlib, cmake)),
            "diagnostic compiler/archive tools are outside system roots")
    compiler_identity = SUPPORT.artifact_identity(compiler, "compiler")
    compile_records, benchmark_object, archive_objects, commands_text = \
        diagnostic_compile_closure(
            source_root, build, commands_path, compiler,
            compiler_invocation, compiler_identity)
    commands_identity = SUPPORT.artifact_identity(commands_path, "build_metadata")
    require(commands_identity["sha256"] ==
            SUPPORT.sha256_bytes(commands_text.encode("utf-8")),
            "diagnostic compile commands changed while they were validated")
    archive_recipe = build / "CMakeFiles/leopard_test_hooks.dir/link.txt"
    require(archive_recipe.is_file(), "hooks archive recipe is missing")
    archive_closure = validate_archive_closure(
        build, hook_archive, archive_recipe, archive_objects, archiver, ranlib,
        archiver_invocation, ranlib_invocation)
    link_recipe = build / "CMakeFiles" / f"{TARGET_NAME}.dir" / "link.txt"
    require(link_recipe.is_file(), "diagnostic benchmark link recipe is missing")
    link_text = link_recipe.read_text(encoding="utf-8", errors="strict")
    link_identity = SUPPORT.artifact_identity(link_recipe, "build_metadata")
    require(link_identity["sha256"] ==
            SUPPORT.sha256_bytes(link_text.encode("utf-8")),
            "diagnostic executable link recipe changed while it was validated")
    link_lines = [line for line in link_text.splitlines() if line.strip()]
    require(len(link_lines) == 1, "diagnostic benchmark has multiple link commands")
    tokens = shlex.split(link_lines[0])
    SUPPORT.validate_effective_flags(
        tokens, "diagnostic executable link recipe", "link")
    link_closure = validate_executable_link_inputs(
        tokens, build, source_root, compiler, compiler_invocation,
        compiler_identity, binary, benchmark_object, hook_archive)
    binary_identity = SUPPORT.artifact_identity(binary, "executable")
    require(binary_identity["mtime_ns"] >= archive_closure["archive"]["mtime_ns"] and
            binary_identity["mtime_ns"] >=
                SUPPORT.artifact_identity(benchmark_object, "object_file")["mtime_ns"],
            "diagnostic executable predates a retained link input")
    clean_rebuild = validate_clean_rebuild(
        source_root, cache, binary_identity, archive_closure["archive"],
        link_text, archive_closure["recipe_text"], cmake_invocation, cmake)
    with tempfile.TemporaryDirectory(prefix="leo2-high-copy-relink-") as temp_root:
        relinked = Path(temp_root) / binary.name
        relink_tokens = list(tokens)
        relink_tokens[relink_tokens.index("-o") + 1] = str(relinked)
        SUPPORT.run_checked(
            relink_tokens, cwd=build, environment=SUPPORT.CHILD_ENVIRONMENT,
            timeout=120, max_stdout=1024 * 1024, max_stderr=1024 * 1024)
        require(SUPPORT.sha256_file(relinked) == binary_identity["sha256"],
                "diagnostic executable bytes differ from a deterministic relink")
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
        "cache": cache_identity,
        "compile_commands": commands_identity,
        "compile_commands_text": commands_text,
        "validated_compile_closure": compile_records,
        "link_recipe": link_identity,
        "link_recipe_text": link_text,
        "validated_link_closure": link_closure,
        "validated_archive_closure": archive_closure,
        "validated_clean_rebuild": clean_rebuild,
        "binary": binary_identity,
        "hook_archive": archive_closure["archive"],
        "compiler_invocation": compiler_invocation,
        "compiler_invocation_identity": SUPPORT.artifact_identity(
            Path(compiler_invocation), "compiler"),
        "compiler": compiler_identity,
        "cmake_invocation": cmake_invocation,
        "cmake_invocation_identity": SUPPORT.artifact_identity(
            Path(cmake_invocation), "build_tool"),
        "cmake": SUPPORT.artifact_identity(cmake, "build_tool"),
        "system_link_inputs": [
            SUPPORT.artifact_identity(
                Path(item["resolved_path"]), "system_link_dependency")
            for item in link_closure["system_inputs"]],
        "deterministic_relink_sha256": binary_identity["sha256"],
        "nm": SUPPORT.artifact_identity(nm, "executable"),
        "selector_symbols_present_in_hook_archive_and_binary": True,
    }


def normalized_runtime_closure(ldd: Path, binary: Path) -> dict[str, Any]:
    closure = SUPPORT.runtime_closure(ldd, binary)
    for dependency in closure["dependencies"]:
        if "file" in dependency:
            dependency["loader_path"] = dependency["file"]["path"]
    return closure


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
        "runtime": normalized_runtime_closure(Path("/usr/bin/ldd"), binary),
    }


def canonical_absolute_path(value: object, label: str) -> Path:
    require(type(value) is str and value and os.path.isabs(value) and
            os.path.normpath(value) == value,
            f"{label} is not a canonical absolute path")
    return Path(value)


def lexical_path(value: str, directory: Path) -> Path:
    candidate = Path(value)
    if not candidate.is_absolute():
        candidate = directory / candidate
    return Path(os.path.normpath(str(candidate)))


def system_path(path: Path) -> bool:
    return any(path.is_relative_to(root)
               for root in (Path("/usr/bin"), Path("/usr/lib"),
                             Path("/usr/lib64"), Path("/lib"), Path("/lib64")))


def cxx_driver_path(path: Path) -> bool:
    return system_path(path) and re.fullmatch(
        r"(?:[^/]+-)?(?:g\+\+|c\+\+|clang(?:\+\+)?)(?:-\d+)?",
        path.name) is not None


def archive_tool_path(path: Path, name: str) -> bool:
    require(name in ("ar", "ranlib"), "unknown archive tool class")
    return system_path(path) and re.fullmatch(
        rf"(?:[^/]+-)?(?:llvm-)?{name}(?:-\d+)?", path.name) is not None


def cmake_tool_path(path: Path) -> bool:
    return system_path(path) and re.fullmatch(r"cmake(?:-\d+(?:\.\d+)*)?", path.name) \
        is not None


def validate_artifact_record(
    value: object, label: str, *, expected_path: Path | None = None,
    expected_kind: str | None = None,
) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
                "path", "kind", "size", "mode", "mtime_ns", "sha256"},
            f"{label} artifact shape changed")
    path = canonical_absolute_path(value.get("path"), f"{label} path")
    require(expected_path is None or path == expected_path,
            f"{label} artifact path changed")
    kind = value.get("kind")
    require(type(kind) is str and kind and
            (expected_kind is None or kind == expected_kind),
            f"{label} artifact kind changed")
    require(type(value.get("size")) is int and value["size"] > 0 and
            type(value.get("mode")) is int and 0 <= value["mode"] <= 0o7777 and
            type(value.get("mtime_ns")) is int and value["mtime_ns"] >= 0 and
            type(value.get("sha256")) is str and
            HEX64.fullmatch(value["sha256"]) is not None,
            f"{label} artifact metadata is invalid")
    if kind in ("executable", "compiler", "archiver", "ranlib", "build_tool"):
        require(value["mode"] & 0o111 != 0,
                f"{label} executable artifact has no execute bit")
    return value


def validate_source_identity(value: object) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
                "path", "head", "tree", "detached",
                "tracked_tree_listing_sha256", "tracked_status"},
            "retained source identity shape changed")
    canonical_absolute_path(value.get("path"), "source root")
    require(type(value.get("head")) is str and
            HEX40.fullmatch(value["head"]) is not None and
            type(value.get("tree")) is str and
            HEX40.fullmatch(value["tree"]) is not None and
            value.get("detached") is False and
            type(value.get("tracked_tree_listing_sha256")) is str and
            HEX64.fullmatch(value["tracked_tree_listing_sha256"]) is not None and
            value.get("tracked_status") == "clean",
            "retained source identity is not a clean attached commit")
    return value


def validate_compile_snapshot(
    build_value: Mapping[str, Any], source_root: Path, build: Path,
    compiler: Mapping[str, Any], compiler_invocation: str,
    evidence_schema: str,
) -> tuple[list[dict[str, Any]], Path]:
    commands_path = build / "compile_commands.json"
    commands_artifact = validate_artifact_record(
        build_value.get("compile_commands"), "compile commands",
        expected_path=commands_path, expected_kind="build_metadata")
    commands_text = build_value.get("compile_commands_text")
    require(type(commands_text) is str and commands_text and
            len(commands_text.encode("utf-8")) == commands_artifact["size"] and
            SUPPORT.sha256_bytes(commands_text.encode("utf-8")) ==
                commands_artifact["sha256"],
            "retained compile commands text differs from its artifact identity")
    try:
        entries = json.loads(commands_text)
    except json.JSONDecodeError as error:
        raise EvidenceError(
            f"retained compile commands are not JSON: {error}") from error
    require(type(entries) is list and all(isinstance(item, dict) for item in entries),
            "retained compile commands are not an object array")

    include_gfni = compile_entries_include_gfni(entries, build)
    specifications = diagnostic_specifications(
        source_root, build, evidence_schema, include_gfni)
    records = build_value.get("validated_compile_closure")
    require(type(records) is list and len(records) == len(specifications) and
            all(isinstance(record, dict) for record in records),
            "retained compile closure has the wrong action count")
    compiler_path = canonical_absolute_path(compiler["path"], "compiler path")
    seen_objects: set[Path] = set()
    for record, (source, output, definitions) in zip(records, specifications):
        require(set(record) == {
                    "source", "object", "directory", "command_tokens",
                    "compiler_invocation", "compiler_resolved_path",
                    "compiler_invocation_identity", "required_definitions"},
                "retained compile action shape changed")
        source_artifact = validate_artifact_record(
            record.get("source"), "compile source",
            expected_path=source, expected_kind="source_file")
        object_artifact = validate_artifact_record(
            record.get("object"), "compile object",
            expected_path=output, expected_kind="object_file")
        require(object_artifact["mtime_ns"] >= source_artifact["mtime_ns"] and
                output not in seen_objects,
                "retained compile object is stale or duplicated")
        seen_objects.add(output)
        require(record.get("directory") == str(build) and
                record.get("compiler_invocation") == compiler_invocation and
                record.get("compiler_resolved_path") == str(compiler_path) and
                record.get("compiler_invocation_identity") == compiler and
                record.get("required_definitions") == sorted(definitions),
                "retained compile action directory/compiler/definitions changed")
        tokens = record.get("command_tokens")
        require(type(tokens) is list and tokens and
                all(type(token) is str and token for token in tokens) and
                not any(token.startswith("@") for token in tokens) and
                tokens.count("-c") == 1 and tokens.count("-o") == 1,
                "retained compile action tokens are indirect or malformed")
        source_index = tokens.index("-c")
        output_index = tokens.index("-o")
        require(source_index + 1 < len(tokens) and output_index + 1 < len(tokens) and
                lexical_path(tokens[source_index + 1], build) == source and
                lexical_path(tokens[output_index + 1], build) == output and
                tokens[0] == compiler_invocation and
                cxx_driver_path(canonical_absolute_path(
                    tokens[0], "compile compiler invocation")),
                "retained compile action source/output/compiler changed")
        validate_diagnostic_compile_flags(
            tokens, f"retained compile action {output}", source_root,
            output, definitions, include_gfni)
        require("-fopenmp" in tokens or "-fopenmp=libomp" in tokens,
                "retained compile action lost OpenMP")
        define_tokens = {token[2:] for token in tokens if token.startswith("-D")}
        require(definitions.issubset(define_tokens),
                "retained compile action lost a private definition")

        matching = []
        for entry in entries:
            directory = entry.get("directory")
            output_value = entry.get("output")
            if type(directory) is not str or type(output_value) is not str:
                continue
            if lexical_path(output_value, Path(directory)) != output:
                continue
            matching.append(entry)
        require(len(matching) == 1,
                "retained compile_commands text omits or duplicates an action")
        entry = matching[0]
        require(entry.get("directory") == str(build) and
                type(entry.get("file")) is str and
                lexical_path(entry["file"], build) == source and
                SUPPORT.compile_command_tokens(entry) == tokens,
                "retained compile action differs from compile_commands text")
    return records, specifications[-1][1]


def validate_archive_snapshot(
    value: object, build: Path, compile_records: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    require(isinstance(value, dict) and set(value) == {
                "recipe", "recipe_text", "recipe_tokens",
                "archiver_invocation", "ranlib_invocation", "archive",
                "members", "archiver", "ranlib",
                "archiver_invocation_identity", "ranlib_invocation_identity"},
            "retained archive closure shape changed")
    recipe_path = build / "CMakeFiles/leopard_test_hooks.dir/link.txt"
    recipe = validate_artifact_record(
        value.get("recipe"), "archive recipe", expected_path=recipe_path,
        expected_kind="build_metadata")
    text = value.get("recipe_text")
    require(type(text) is str and text and
            len(text.encode("utf-8")) == recipe["size"] and
            SUPPORT.sha256_bytes(text.encode("utf-8")) == recipe["sha256"],
            "retained archive recipe text differs from its artifact")
    try:
        parsed = [shlex.split(line) for line in text.splitlines() if line.strip()]
    except ValueError as error:
        raise EvidenceError(f"retained archive recipe is malformed: {error}") from error
    tokens = value.get("recipe_tokens")
    require(type(tokens) is list and tokens == parsed and len(tokens) == 2 and
            all(type(command) is list and
                all(type(token) is str and token for token in command)
                for command in tokens),
            "retained archive tokenization differs from its recipe")
    archive_tokens, ranlib_tokens = tokens
    require(len(archive_tokens) == 3 + len(compile_records) - 1 and
            archive_tokens[0] == value.get("archiver_invocation") and
            archive_tokens[1] in ("qc", "rc", "rcs") and
            len(ranlib_tokens) == 2 and
            ranlib_tokens[0] == value.get("ranlib_invocation") and
            ranlib_tokens[1] == archive_tokens[2] and
            archive_tool_path(canonical_absolute_path(
                archive_tokens[0], "archiver invocation"), "ar") and
            archive_tool_path(canonical_absolute_path(
                ranlib_tokens[0], "ranlib invocation"), "ranlib"),
            "retained archive/ranlib command shape changed")
    archive_path = build / "libleopard_test_hooks.a"
    require(lexical_path(archive_tokens[2], build) == archive_path,
            "retained archive recipe output changed")
    expected_objects = [Path(record["object"]["path"])
                        for record in compile_records[:-1]]
    require([lexical_path(token, build) for token in archive_tokens[3:]] ==
                expected_objects and len(set(expected_objects)) == len(expected_objects),
            "retained archive recipe object closure changed")
    archive = validate_artifact_record(
        value.get("archive"), "hooks archive", expected_path=archive_path,
        expected_kind="archive")
    archiver = validate_artifact_record(
        value.get("archiver"), "archiver", expected_kind="archiver")
    ranlib = validate_artifact_record(
        value.get("ranlib"), "ranlib", expected_kind="ranlib")
    require(value.get("archiver_invocation_identity") == archiver and
            value.get("ranlib_invocation_identity") == ranlib and
            archive_tool_path(Path(archiver["path"]), "ar") and
            (archive_tool_path(Path(ranlib["path"]), "ranlib") or
             archive_tool_path(Path(ranlib["path"]), "ar")),
            "retained archive tools are outside system roots")
    members = value.get("members")
    require(type(members) is list and len(members) == len(expected_objects),
            "retained archive member count changed")
    for member, object_record, expected_object in zip(
            members, compile_records[:-1], expected_objects):
        require(isinstance(member, dict) and set(member) == {"member", "sha256"} and
                member.get("member") == expected_object.name and
                type(member.get("sha256")) is str and
                HEX64.fullmatch(member["sha256"]) is not None and
                member["sha256"] == object_record["object"]["sha256"],
                "retained archive member differs from its compiled object")
    require(archive["mtime_ns"] >=
            max(record["object"]["mtime_ns"] for record in compile_records[:-1]),
            "retained archive predates a compiled member")
    return value


def validate_link_snapshot(
    build_value: Mapping[str, Any], build: Path,
    benchmark_object: Path, archive: Mapping[str, Any],
    compiler: Mapping[str, Any], compiler_invocation: str,
) -> None:
    binary_path = build / TARGET_NAME
    link_path = build / "CMakeFiles" / f"{TARGET_NAME}.dir" / "link.txt"
    link = validate_artifact_record(
        build_value.get("link_recipe"), "executable link recipe",
        expected_path=link_path, expected_kind="build_metadata")
    text = build_value.get("link_recipe_text")
    require(type(text) is str and text and
            len(text.encode("utf-8")) == link["size"] and
            SUPPORT.sha256_bytes(text.encode("utf-8")) == link["sha256"],
            "retained executable link text differs from its artifact")
    lines = [line for line in text.splitlines() if line.strip()]
    require(len(lines) == 1, "retained executable has multiple link commands")
    try:
        parsed_tokens = shlex.split(lines[0])
    except ValueError as error:
        raise EvidenceError(f"retained executable link is malformed: {error}") from error
    closure = build_value.get("validated_link_closure")
    require(isinstance(closure, dict) and set(closure) == {
                "recipe_tokens", "compiler_invocation", "compiler_resolved_path",
                "compiler_invocation_identity",
                "output_token", "output_path", "build_inputs", "system_inputs"},
            "retained executable link closure shape changed")
    tokens = closure.get("recipe_tokens")
    require(type(tokens) is list and tokens == parsed_tokens and tokens and
            all(type(token) is str and token for token in tokens) and
            tokens.count("-o") == 1 and not any(token.startswith("@") for token in tokens),
            "retained executable link tokenization changed")
    output_index = tokens.index("-o")
    require(output_index + 1 < len(tokens) and
            tokens[0] == compiler_invocation and
            closure.get("compiler_invocation") == compiler_invocation and
            closure.get("compiler_resolved_path") == compiler["path"] and
            closure.get("compiler_invocation_identity") == compiler and
            cxx_driver_path(canonical_absolute_path(
                tokens[0], "link compiler invocation")) and
            tokens[output_index + 1] == closure.get("output_token") and
            lexical_path(tokens[output_index + 1], build) == binary_path and
            closure.get("output_path") == str(binary_path),
            "retained executable compiler/output closure changed")
    SUPPORT.validate_effective_flags(
        tokens, "retained executable link recipe", "link")
    explicit_tokens: list[str] = []
    for index, token in enumerate(tokens[1:], start=1):
        if index in (output_index, output_index + 1):
            continue
        if token.startswith("-"):
            require(not token.startswith(("-L", "-l", "-Xlinker")) and
                    (not token.startswith("-Wl,") or
                     token.startswith("-Wl,-rpath,")),
                    "retained executable contains an undeclared linker input flag")
            if token.startswith("-Wl,-rpath,"):
                require(system_path(canonical_absolute_path(
                    token.split(",", 2)[2], "link rpath")),
                    "retained executable rpath is outside system roots")
            continue
        explicit_tokens.append(token)
    build_inputs = closure.get("build_inputs")
    system_inputs = closure.get("system_inputs")
    require(type(build_inputs) is list and len(build_inputs) == 2 and
            type(system_inputs) is list and
            all(isinstance(item, dict) and
                set(item) == {"token", "resolved_path"}
                for item in [*build_inputs, *system_inputs]) and
            [item["token"] for item in [*build_inputs, *system_inputs]] ==
                explicit_tokens,
            "retained executable input token closure changed")
    expected_build_paths = [benchmark_object, Path(archive["archive"]["path"])]
    for item, expected in zip(build_inputs, expected_build_paths):
        require(type(item["token"]) is str and
                lexical_path(item["token"], build) == expected and
                item["resolved_path"] == str(expected),
                "retained executable build input changed")
    system_artifacts = build_value.get("system_link_inputs")
    require(type(system_artifacts) is list and
            len(system_artifacts) == len(system_inputs),
            "retained executable system input count changed")
    allowed_system_names = re.compile(r"lib(?:gomp|omp|pthread)\.(?:a|so(?:\..*)?)$")
    for item, artifact in zip(system_inputs, system_artifacts):
        resolved = canonical_absolute_path(
            item.get("resolved_path"), "system link input")
        token_path = canonical_absolute_path(item.get("token"), "system link token")
        validate_artifact_record(
            artifact, "system link input", expected_path=resolved,
            expected_kind="system_link_dependency")
        require(system_path(token_path) and system_path(resolved) and
                allowed_system_names.fullmatch(resolved.name) is not None,
                "retained executable has an undeclared system link input")
    binary = validate_artifact_record(
        build_value.get("binary"), "diagnostic binary", expected_path=binary_path,
        expected_kind="executable")
    require(binary["mtime_ns"] >= archive["archive"]["mtime_ns"] and
            binary["mtime_ns"] >=
                next(record for record in build_value["validated_compile_closure"]
                     if record["object"]["path"] == str(benchmark_object))
                ["object"]["mtime_ns"] and
            build_value.get("deterministic_relink_sha256") == binary["sha256"],
            "retained binary is stale or differs from its deterministic relink")


def validate_clean_rebuild_snapshot(
    value: object, build_value: Mapping[str, Any], archive: Mapping[str, Any],
    cmake: Mapping[str, Any], cmake_invocation: str,
) -> None:
    require(isinstance(value, dict) and set(value) == {
                "cmake", "cmake_invocation", "cmake_resolved_path",
                "cmake_invocation_identity", "generator", "retained_options",
                "archive_sha256", "binary_sha256"},
            "retained clean-rebuild identity shape changed")
    clean_cmake = validate_artifact_record(
        value.get("cmake"), "clean-rebuild CMake", expected_kind="build_tool")
    require(clean_cmake == cmake and
            value.get("cmake_invocation") == cmake_invocation and
            value.get("cmake_resolved_path") == cmake["path"] and
            value.get("cmake_invocation_identity") == cmake and
            cmake_tool_path(canonical_absolute_path(
                cmake_invocation, "CMake invocation")) and
            cmake_tool_path(Path(cmake["path"])) and
            value.get("generator") == "Unix Makefiles",
            "retained clean-rebuild tool/generator is invalid")
    options = value.get("retained_options")
    expected_keys = {
        "CMAKE_BUILD_TYPE", "CMAKE_AR", "CMAKE_CXX_COMPILER",
        "CMAKE_CXX_FLAGS", "CMAKE_CXX_FLAGS_RELEASE",
        "CMAKE_EXPORT_COMPILE_COMMANDS", "CMAKE_RANLIB", "ENABLE_OPENMP",
        "LEO2_BACKEND_VARIANT", "LEO2_BUILD_BENCHMARKS",
        "LEO2_BUILD_FUZZERS", "LEO2_BUILD_TESTS", "LEO2_ENABLE_CUDA",
        "LEOPARD_ENABLE_GF8", "LEOPARD_ENABLE_GF16",
    }
    require(isinstance(options, dict) and set(options) == expected_keys and
            all(type(item) is str for item in options.values()) and
            options["CMAKE_BUILD_TYPE"] == "Release" and
            options["CMAKE_CXX_FLAGS_RELEASE"] == "-O3 -DNDEBUG" and
            options["CMAKE_EXPORT_COMPILE_COMMANDS"] == "ON" and
            options["ENABLE_OPENMP"] == "ON" and
            options["LEO2_BACKEND_VARIANT"] == "auto" and
            options["LEO2_BUILD_BENCHMARKS"] == "ON" and
            options["LEO2_BUILD_FUZZERS"] == "OFF" and
            options["LEO2_BUILD_TESTS"] == "ON" and
            options["LEO2_ENABLE_CUDA"] == "OFF" and
            options["LEOPARD_ENABLE_GF8"] == "ON" and
            options["LEOPARD_ENABLE_GF16"] == "ON" and
            options["CMAKE_AR"] == archive["archiver_invocation"] and
            options["CMAKE_RANLIB"] == archive["ranlib_invocation"] and
            options["CMAKE_CXX_COMPILER"] == build_value["compiler_invocation"] and
            value.get("archive_sha256") == archive["archive"]["sha256"] and
            value.get("binary_sha256") == build_value["binary"]["sha256"],
            "retained clean-rebuild options or output hashes changed")


def validate_runtime_snapshot(value: object, binary: Path) -> None:
    require(isinstance(value, dict) and set(value) == {
                "executable", "dependencies"} and
            value.get("executable") == str(binary),
            "retained runtime closure shape/executable changed")
    dependencies = value.get("dependencies")
    require(type(dependencies) is list and dependencies and
            all(isinstance(item, dict) for item in dependencies) and
            dependencies == sorted(dependencies, key=lambda item: item.get("soname", "")) and
            len({item.get("soname") for item in dependencies}) == len(dependencies),
            "retained runtime dependency set is empty, duplicate, or unsorted")
    for item in dependencies:
        soname = item.get("soname")
        require(type(soname) is str and soname,
                "retained runtime dependency has no soname")
        if item.get("virtual") is True:
            require(set(item) == {"soname", "virtual"} and
                    soname == "linux-vdso.so.1",
                    "retained virtual runtime dependency changed")
            continue
        require(set(item) == {"soname", "loader_path", "file"},
                "retained runtime dependency shape changed")
        path = canonical_absolute_path(item.get("loader_path"),
                                       "runtime loader path")
        kind = "dynamic_loader" if Path(soname).name == path.name and \
            (soname.startswith("ld-linux") or "ld-" in soname) \
            else "shared_library"
        file_artifact = validate_artifact_record(
            item.get("file"), "runtime dependency", expected_path=path,
            expected_kind=kind)
        file_path = Path(file_artifact["path"])
        require(system_path(path) and system_path(file_path) and
                (file_path.name == soname or file_path.name.startswith(soname + ".")),
                "retained runtime dependency is outside system roots")


def validate_snapshot_identity(
    value: object, evidence_schema: str = RAW_SCHEMA,
) -> dict[str, Any]:
    require(type(evidence_schema) is str and evidence_schema in RAW_SCHEMAS,
            "unknown retained high-decode copy evidence schema")
    require(isinstance(value, dict) and set(value) == {
                "source", "build", "runner", "contract", "matrix",
                "support", "runtime"},
            "retained source/build snapshot shape changed")
    source = validate_source_identity(value.get("source"))
    source_root = Path(source["path"])
    build_value = value.get("build")
    expected_build_keys = {
        "build_root", "cache", "compile_commands", "compile_commands_text",
        "validated_compile_closure", "link_recipe", "link_recipe_text",
        "validated_link_closure", "validated_archive_closure",
        "validated_clean_rebuild", "binary", "hook_archive", "compiler",
        "compiler_invocation", "compiler_invocation_identity",
        "cmake", "cmake_invocation", "cmake_invocation_identity",
        "system_link_inputs", "deterministic_relink_sha256", "nm",
        "selector_symbols_present_in_hook_archive_and_binary",
    }
    require(isinstance(build_value, dict) and set(build_value) == expected_build_keys and
            build_value.get("selector_symbols_present_in_hook_archive_and_binary") is True,
            "retained diagnostic build snapshot shape changed")
    build = canonical_absolute_path(build_value.get("build_root"), "build root")
    validate_artifact_record(
        build_value.get("cache"), "CMake cache",
        expected_path=build / "CMakeCache.txt", expected_kind="build_metadata")
    compiler = validate_artifact_record(
        build_value.get("compiler"), "compiler", expected_kind="compiler")
    compiler_invocation = build_value.get("compiler_invocation")
    require(type(compiler_invocation) is str and
            cxx_driver_path(canonical_absolute_path(
                compiler_invocation, "primary compiler invocation")) and
            build_value.get("compiler_invocation_identity") == compiler,
            "retained primary compiler invocation is invalid")
    cmake = validate_artifact_record(
        build_value.get("cmake"), "primary CMake", expected_kind="build_tool")
    cmake_invocation = build_value.get("cmake_invocation")
    require(type(cmake_invocation) is str and
            cmake_tool_path(canonical_absolute_path(
                cmake_invocation, "primary CMake invocation")) and
            cmake_tool_path(Path(cmake["path"])) and
            build_value.get("cmake_invocation_identity") == cmake,
            "retained primary CMake invocation/identity is invalid")
    nm = validate_artifact_record(
        build_value.get("nm"), "nm", expected_kind="executable")
    nm_name = Path(nm["path"]).name
    require(cxx_driver_path(Path(compiler["path"])) and
            system_path(Path(nm["path"])) and
            re.fullmatch(r"(?:[^/]+-)?nm(?:-\d+)?", nm_name) is not None,
            "retained compiler/nm are outside system roots")
    compile_records, benchmark_object = validate_compile_snapshot(
        build_value, source_root, build, compiler, compiler_invocation,
        evidence_schema)
    archive = validate_archive_snapshot(
        build_value.get("validated_archive_closure"), build, compile_records)
    require(build_value.get("hook_archive") == archive["archive"],
            "top-level hooks archive differs from the validated archive closure")
    validate_link_snapshot(
        build_value, build, benchmark_object, archive, compiler,
        compiler_invocation)
    validate_clean_rebuild_snapshot(
        build_value.get("validated_clean_rebuild"), build_value, archive,
        cmake, cmake_invocation)
    binary = Path(build_value["binary"]["path"])
    for name, relative in (
        ("runner", RUNNER_RELATIVE), ("contract", CONTRACT_RELATIVE),
        ("matrix", MATRIX_RELATIVE),
        ("support", "experiments/leopard2/main_compare/run_abba.py")):
        validate_artifact_record(
            value.get(name), name, expected_path=source_root / relative,
            expected_kind="source_file")
    validate_runtime_snapshot(value.get("runtime"), binary)
    return value


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


def expected_rate(byte_count: int, microseconds: float) -> float:
    return byte_count / (microseconds * 1000.0)


def validate_rate(value: object, expected: float, name: str) -> None:
    actual = SUPPORT.finite_number(value, name)
    require(SUPPORT.close_enough(actual, expected),
            f"{name} is not derived from bytes and retained median")


def validate_timing_metrics(
    metrics: object, cell: Mapping[str, Any], controls: Mapping[str, Any]
) -> dict[str, float]:
    require(isinstance(metrics, dict) and set(metrics) == {
                "codec_setup", "encode_execution", "decode_plan_setup",
                "decode_execution", "decode_amortized_at_reuse", "rate_semantics"},
            "benchmark timing metrics are incomplete or contain unbound claims")
    iterations = controls.get("iterations")
    reuse = controls.get("reuse")
    require(type(iterations) is int and iterations >= 3 and
            type(reuse) is int and reuse >= 1,
            "timing controls are not exact bounded integers")
    SUPPORT.validate_summary(metrics["codec_setup"], iterations, setup=True)
    encode_samples = SUPPORT.validate_summary(
        metrics["encode_execution"], iterations)
    plan_samples = SUPPORT.validate_summary(
        metrics["decode_plan_setup"], iterations, setup=True)
    decode_samples = SUPPORT.validate_summary(
        metrics["decode_execution"], iterations)
    del encode_samples, plan_samples, decode_samples
    encode_median = SUPPORT.finite_number(
        metrics["encode_execution"].get("median_us_per_batch_call"),
        "encode execution median")
    plan_median = SUPPORT.finite_number(
        metrics["decode_plan_setup"].get("median_us"), "decode plan median")
    decode_median = SUPPORT.finite_number(
        metrics["decode_execution"].get("median_us_per_batch_call"),
        "decode execution median")
    encode_input_bytes = cell["k"] * cell["shard_bytes"]
    encode_output_bytes = cell["r"] * cell["shard_bytes"]
    decode_input_bytes = (cell["k"] - cell["losses"] + cell["r"]) * \
        cell["shard_bytes"]
    decode_output_bytes = cell["losses"] * cell["shard_bytes"]
    validate_rate(metrics["encode_execution"].get("input_GB_per_s"),
                  expected_rate(encode_input_bytes, encode_median),
                  "encode input rate")
    validate_rate(metrics["encode_execution"].get("parity_output_GB_per_s"),
                  expected_rate(encode_output_bytes, encode_median),
                  "encode parity rate")
    validate_rate(metrics["decode_execution"].get("offered_received_GB_per_s"),
                  expected_rate(decode_input_bytes, decode_median),
                  "decode offered-input rate")
    validate_rate(metrics["decode_execution"].get("repaired_output_GB_per_s"),
                  expected_rate(decode_output_bytes, decode_median),
                  "decode repaired-output rate")
    amortized = metrics["decode_amortized_at_reuse"]
    require(isinstance(amortized, dict) and set(amortized) == {
                "reuse_count", "derived_median_us_per_batch_call",
                "offered_received_GB_per_s", "repaired_output_GB_per_s"} and
            type(amortized.get("reuse_count")) is int and
            amortized["reuse_count"] == reuse,
            "amortized decode timing shape/reuse changed")
    derived_amortized = decode_median + plan_median / reuse
    actual_amortized = SUPPORT.finite_number(
        amortized.get("derived_median_us_per_batch_call"),
        "amortized decode median")
    require(SUPPORT.close_enough(actual_amortized, derived_amortized),
            "amortized decode median is not plan/reuse plus execution")
    validate_rate(amortized.get("offered_received_GB_per_s"),
                  expected_rate(decode_input_bytes, derived_amortized),
                  "amortized offered-input rate")
    validate_rate(amortized.get("repaired_output_GB_per_s"),
                  expected_rate(decode_output_bytes, derived_amortized),
                  "amortized repaired-output rate")
    require(metrics["rate_semantics"] ==
            "offered_received counts all non-null shard pointers supplied; "
            "a plan may read a deterministic subset",
            "decode rate semantics changed")
    return {
        "decode_execution_median_us": decode_median,
        "decode_amortized_median_us": actual_amortized,
    }


def validate_result(
    document: object, cell: Mapping[str, Any], mode: str,
    controls: Mapping[str, Any]
) -> dict[str, Any]:
    result = CONTRACT.validate_document(
        document, mode=mode, workspace=cell["workspace"], field=cell["field"],
        tail_bytes=cell["shard_bytes"] % 64,
        evaluator="mature" if cell["role"] == "full-block" else None)
    parameters = result["parameters"]
    resolved = result["resolved"]
    numeric_parameters = {
        "K": cell["k"], "R": cell["r"],
        "shard_bytes": cell["shard_bytes"], "loss_count": cell["losses"],
        "seed": cell["seed"], "batch": 1,
        "reuse": controls["reuse"], "iterations": controls["iterations"],
        "warmup": controls["warmup"], "thread_count": 1,
    }
    require(all(type(parameters.get(name)) is int and
                parameters[name] == expected
                for name, expected in numeric_parameters.items()) and
            parameters.get("requested_backend") == cell["backend"] and
            resolved.get("backend") == cell["backend"],
            "benchmark output does not attest the matrix cell")
    validate_missing_original_indices(
        parameters.get("missing_original_indices"), cell)
    validate_timing_metrics(result.get("metrics"), cell, controls)
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
    controls = {"reuse": reuse, "iterations": iterations, "warmup": warmup}
    result = validate_result(document, cell, mode, controls)
    return {
        "cell": cell["id"], "mode": mode, "round": round_index, "slot": slot,
        "command": command, "command_sha256": SUPPORT.sha256_bytes(
            SUPPORT.canonical_bytes(command)),
        "started_monotonic_ns": started, "ended_monotonic_ns": ended,
        "cpu_before": before_cpu, "cpu_after": after_cpu, "cpu_delta": cpu_delta,
        "sibling_before": before_sibling, "sibling_after": after_sibling,
        "sibling_delta": sibling_delta,
        "result": result,
    }


def grouped_records(records: Sequence[Mapping[str, Any]],
                    cells: Sequence[Mapping[str, Any]]) -> dict[str, list[Mapping[str, Any]]]:
    result = {cell["id"]: [] for cell in cells}
    for record in records:
        identifier = record.get("cell")
        require(identifier in result, "record refers to an unknown cell")
        result[identifier].append(record)
    return result


def analyze(
    records: object, matrix: Mapping[str, Any], campaign: Mapping[str, Any]
) -> dict[str, Any]:
    records = validate_record_sequence(records, matrix)
    grouped = grouped_records(records, matrix["cells"])
    output: dict[str, Any] = {}
    metrics = {
        "decode_execution": (
            "metrics", "decode_execution", "median_us_per_batch_call"),
        "decode_amortized": (
            "metrics", "decode_amortized_at_reuse", "derived_median_us_per_batch_call"),
    }
    for cell in matrix["cells"]:
        cell_records = grouped[cell["id"]]
        require(len(cell_records) == 12, "cell does not have three four-slot rounds")
        documents = {mode: [validate_result(
                                record.get("result"), cell, mode, campaign)
                            for record in cell_records if record["mode"] == mode]
                     for mode in ("no-copy", "copy-fallback")}
        CONTRACT.validate_pair(
            documents["no-copy"][0], documents["copy-fallback"][0],
            workspace=cell["workspace"], field=cell["field"],
            tail_bytes=cell["shard_bytes"] % 64,
            evaluator="mature" if cell["role"] == "full-block" else None)
        digest = documents["no-copy"][0]["workload_digests"]
        missing = expected_missing_original_indices(cell)
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


def validate_campaign_environment(
    raw: Mapping[str, Any], campaign: Mapping[str, Any],
    evidence_schema: str = RAW_SCHEMA,
) -> tuple[dict[str, Any], dict[str, Any]]:
    require(type(evidence_schema) is str and evidence_schema in RAW_SCHEMAS,
            "unknown campaign high-decode copy evidence schema")
    cpu = campaign["cpu"]
    sibling = campaign["reserved_sibling"]
    host_initial = raw.get("host_initial")
    host_final = raw.get("host_final")
    require(host_initial == host_final and isinstance(host_initial, dict),
            "host topology/frequency policy changed during campaign")
    allowed = host_initial.get("allowed_cpu_set_at_launch")
    require(type(allowed) is list and allowed == sorted(set(allowed)) and
            all(type(item) is int and 0 <= item <= SUPPORT.MAX_CPU_ID
                for item in allowed) and
            cpu in allowed and sibling in allowed,
            "host launch affinity record is invalid")
    SUPPORT.validate_host_record(
        host_initial, cpu, sibling, allowed, evidence_schema)
    isolation = validate_campaign_isolation(raw.get("isolation"), campaign)
    pair_lease = raw.get("pair_lease")
    SUPPORT.validate_pair_lease_identity(pair_lease, cpu, sibling)
    require(pair_lease == isolation["pair_lease"],
            "top-level pair lease differs from isolation evidence")
    reservation = raw.get("reservation")
    require(isinstance(reservation, dict) and set(reservation) == {
                "path", "sha256", "payload", "lock"} and
            reservation.get("lock") == "exclusive_nonblocking" and
            type(reservation.get("path")) is str and reservation["path"] and
            isinstance(reservation.get("payload"), dict) and
            SUPPORT.parse_reservation(
                SUPPORT.canonical_bytes(reservation["payload"]), cpu, sibling) ==
                reservation["payload"] and
            reservation.get("sha256") == SUPPORT.sha256_bytes(
                SUPPORT.canonical_bytes(reservation["payload"])),
            "coordinator reservation identity/semantics are invalid")
    return isolation, reservation


def validate_record_execution_attestation(
    record: Mapping[str, Any], *, interval_begin: int, interval_end: int,
    previous_end: int, cpu: int, sibling: int
) -> tuple[int, int]:
    started = record.get("started_monotonic_ns")
    ended = record.get("ended_monotonic_ns")
    require(type(started) is int and type(ended) is int and
            interval_begin <= started < ended <= interval_end and
            started >= previous_end,
            "child monotonic interval is invalid, overlapping, or out of order")
    cpu_delta = SUPPORT.cpu_stat_delta(
        record.get("cpu_before"), record.get("cpu_after"))
    sibling_delta = SUPPORT.cpu_stat_delta(
        record.get("sibling_before"), record.get("sibling_after"))
    require(record.get("cpu_delta") == cpu_delta and cpu_delta["cpu"] == cpu and
            cpu_delta["nonidle_jiffies"] > 0 and
            record.get("sibling_delta") == sibling_delta and
            sibling_delta["cpu"] == sibling and
            sibling_delta["nonidle_jiffies"] == 0,
            "per-child CPU/timestamp attestation differs")
    return ended, ended - started


def validate_raw(raw: object, output: Path, *, current: bool) -> dict[str, Any]:
    SUPPORT.verify_signature(raw, "high-decode copy raw evidence")
    raw_schema = raw.get("schema") if isinstance(raw, dict) else None
    require(isinstance(raw, dict) and type(raw_schema) is str and
            raw_schema in RAW_SCHEMAS and
            set(raw) == {
                "schema", "created_utc", "validity_is_independent_of_speed",
                "matrix", "campaign", "host_initial", "host_final",
                "reservation", "pair_lease", "isolation",
                "identities_initial", "identities_final", "records",
                "analysis", "digest"} and
            raw.get("validity_is_independent_of_speed") is True,
            "raw high-decode copy evidence schema changed")
    parse_canonical_utc(raw.get("created_utc"), "raw created_utc")
    matrix = raw.get("matrix")
    require(matrix == load_matrix(ROOT / MATRIX_RELATIVE),
            "raw matrix differs from the canonical checked-in matrix")
    campaign = raw.get("campaign")
    require(isinstance(campaign, dict) and set(campaign) == {
                "cpu", "reserved_sibling", "reuse", "iterations", "warmup",
                "thread_count", "timeout_seconds", "round_orders",
                "child_environment", "invocation_count"} and
            campaign.get("round_orders") == matrix["round_orders"] and
            campaign.get("child_environment") == SUPPORT.CHILD_ENVIRONMENT and
            type(campaign.get("cpu")) is int and
            0 <= campaign["cpu"] <= SUPPORT.MAX_CPU_ID and
            type(campaign.get("reserved_sibling")) is int and
            0 <= campaign["reserved_sibling"] <= SUPPORT.MAX_CPU_ID and
            campaign["cpu"] != campaign["reserved_sibling"] and
            type(campaign.get("reuse")) is int and campaign["reuse"] >= 1 and
            type(campaign.get("iterations")) is int and
            campaign["iterations"] >= 3 and
            type(campaign.get("warmup")) is int and campaign["warmup"] >= 1 and
            type(campaign.get("thread_count")) is int and
            campaign["thread_count"] == 1 and
            isinstance(campaign.get("timeout_seconds"), (int, float)) and
            not isinstance(campaign["timeout_seconds"], bool) and
            math.isfinite(campaign["timeout_seconds"]) and
            0 < campaign["timeout_seconds"] <= 3600,
            "raw campaign contract is incomplete")
    records = raw.get("records")
    invocation_count = sum(len(order) for order in matrix["round_orders"]) * \
        len(matrix["cells"])
    require(type(campaign.get("invocation_count")) is int and
            campaign["invocation_count"] == invocation_count,
            "raw campaign does not contain the exact bounded invocation count")
    records = validate_record_sequence(records, matrix)
    cpu = campaign["cpu"]
    sibling = campaign["reserved_sibling"]
    isolation, reservation = validate_campaign_environment(
        raw, campaign, raw_schema)
    initial = raw.get("identities_initial")
    require(isinstance(initial, dict) and initial == raw.get("identities_final"),
            "source/build identity changed during the campaign")
    validate_snapshot_identity(initial, raw_schema)
    binary_value = initial.get("build", {}).get("binary", {}).get("path")
    require(type(binary_value) is str and binary_value,
            "diagnostic binary identity path is absent")
    binary = Path(binary_value)
    cell_by_id = {cell["id"]: cell for cell in matrix["cells"]}
    interval_begin = isolation["before"]["monotonic_ns"]
    interval_end = isolation["after"]["monotonic_ns"]
    previous_end = interval_begin
    total_child_duration = 0
    for record in records:
        require(set(record) == RECORD_KEYS,
                "raw record contains an omitted or unverifiable claim")
        cell = cell_by_id[record["cell"]]
        command = benchmark_command(
            binary, cell, record["mode"], cpu, campaign["reuse"],
            campaign["iterations"], campaign["warmup"])
        previous_end, duration = validate_record_execution_attestation(
            record, interval_begin=interval_begin, interval_end=interval_end,
            previous_end=previous_end, cpu=cpu, sibling=sibling)
        total_child_duration += duration
        require(record.get("command") == command and
                record.get("command_sha256") == SUPPORT.sha256_bytes(
                    SUPPORT.canonical_bytes(command)),
                "raw benchmark command attestation differs")
        validate_result(record.get("result"), cell, record["mode"], campaign)
    require(interval_end - interval_begin >= total_child_duration,
            "campaign isolation interval does not cover child durations")
    require(raw.get("analysis") == analyze(records, matrix, campaign),
            "raw analysis differs from retained invocations")
    if current:
        source = initial["source"]
        rebuilt = snapshot(
            Path(source["path"]), source["head"], binary,
            Path(initial["build"]["hook_archive"]["path"]),
            ROOT / MATRIX_RELATIVE)
        require(rebuilt == initial,
                "current source/build closure differs from retained evidence")
        SUPPORT.validate_reservation_current(reservation)
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
                "analysis": analyze(records, matrix, campaign),
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
                "host_initial": host_initial,
                "host_final": host_final,
                "reservation": reservation,
                "pair_lease": pair_lease,
                "identities": initial,
                "isolation": isolation,
                "analysis": raw["analysis"],
            })
            validate_manifest_raw_timestamps(
                manifest["created_utc"], raw["created_utc"])
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
    manifest_schema = manifest.get("schema") if isinstance(manifest, dict) \
        else None
    require(type(manifest_schema) is str and
            manifest_schema in MANIFEST_RAW_SCHEMAS and
            set(manifest) == {
                "schema", "created_utc", "valid",
                "validity_is_independent_of_speed", "raw", "matrix", "campaign",
                "host_initial", "host_final", "reservation", "pair_lease",
                "identities", "isolation", "analysis", "digest"} and
            manifest.get("valid") is True and
            manifest.get("validity_is_independent_of_speed") is True,
            "manifest is not valid high-decode copy evidence")
    parse_canonical_utc(manifest.get("created_utc"), "manifest created_utc")
    raw_info = manifest.get("raw")
    require(isinstance(raw_info, dict), "manifest has no raw identity")
    raw_path = SUPPORT.safe_evidence_path(path.parent, raw_info.get("path"))
    require(raw_path.is_file() and raw_path.stat().st_size == raw_info.get("size") and
            SUPPORT.sha256_file(raw_path) == raw_info.get("sha256"),
            "raw high-decode copy bundle identity differs")
    raw = json.loads(raw_path.read_text(encoding="utf-8"))
    validate_manifest_raw_schema_pair(
        manifest_schema, raw.get("schema") if isinstance(raw, dict) else None)
    analysis = validate_raw(raw, path.parent, current=not options.no_current_input_check)
    validate_manifest_raw_timestamps(
        manifest.get("created_utc"), raw.get("created_utc"))
    require(raw.get("digest") == raw_info.get("payload_digest") and
            manifest.get("matrix") == raw.get("matrix") and
            manifest.get("campaign") == raw.get("campaign") and
            manifest.get("host_initial") == raw.get("host_initial") and
            manifest.get("host_final") == raw.get("host_final") and
            manifest.get("reservation") == raw.get("reservation") and
            manifest.get("pair_lease") == raw.get("pair_lease") and
            manifest.get("identities") == raw.get("identities_initial") and
            manifest.get("isolation") == raw.get("isolation") and
            manifest.get("analysis") == analysis,
            "manifest differs from its raw high-decode copy bundle")
    print("high-decode copy/no-copy ABBA evidence verified")
    return 0


def validate_manifest_raw_schema_pair(
    manifest_schema: object, raw_schema: object,
) -> tuple[str, str]:
    require(type(manifest_schema) is str and
            manifest_schema in MANIFEST_RAW_SCHEMAS,
            "unknown high-decode copy manifest schema")
    require(type(raw_schema) is str and
            raw_schema == MANIFEST_RAW_SCHEMAS[manifest_schema],
            "manifest/raw high-decode copy schema versions differ")
    return manifest_schema, raw_schema


def build_smoke(options: argparse.Namespace) -> int:
    identity = build_identity(
        options.source_root, options.binary, options.hook_archive)
    require(identity.get(
        "selector_symbols_present_in_hook_archive_and_binary") is True,
        "diagnostic selector symbol gate did not pass")
    print("high-decode copy diagnostic build identity verified")
    return 0


def synthetic_summary(samples: Sequence[float], *, setup: bool) -> dict[str, Any]:
    retained = list(samples)
    middle = statistics.median(retained)
    deviations = [abs(value - middle) for value in retained]
    suffix = "" if setup else "_per_batch_call"
    return {
        f"median_us{suffix}": middle,
        f"mad_us{suffix}": statistics.median(deviations),
        f"minimum_us{suffix}": min(retained),
        f"maximum_us{suffix}": max(retained),
        "samples_us" if setup else "samples_us_per_batch_call": retained,
    }


def synthetic_result(
    cell: Mapping[str, Any], controls: Mapping[str, Any], mode: str = "no-copy"
) -> dict[str, Any]:
    document = CONTRACT.synthetic_document(mode)
    document["parameters"].update({
        "K": cell["k"], "R": cell["r"],
        "requested_field": cell["field"],
        "requested_backend": cell["backend"],
        "force_materialized_decode": cell["workspace"] == "materialized",
        "force_tiled_decode": cell["workspace"] == "tiled",
        "shard_bytes": cell["shard_bytes"], "loss_count": cell["losses"],
        "missing_original_indices": expected_missing_original_indices(cell),
        "batch": 1, "reuse": controls["reuse"],
        "iterations": controls["iterations"], "warmup": controls["warmup"],
        "thread_count": 1, "seed": cell["seed"],
    })
    document["resolved"].update({
        "field": cell["field"], "backend": cell["backend"],
        "selected_decode_path": cell["workspace"],
        "selected_decode_rule": "forced_" + cell["workspace"],
        "decode_tail_bytes": cell["shard_bytes"] % 64,
    })
    encode = synthetic_summary((1.0, 2.0, 3.0), setup=False)
    plan = synthetic_summary((3.0, 4.0, 5.0), setup=True)
    decode = synthetic_summary((5.0, 6.0, 7.0), setup=False)
    codec = synthetic_summary((2.0, 3.0, 4.0), setup=True)
    encode_median = encode["median_us_per_batch_call"]
    decode_median = decode["median_us_per_batch_call"]
    plan_median = plan["median_us"]
    amortized = decode_median + plan_median / controls["reuse"]
    encode_input = cell["k"] * cell["shard_bytes"]
    encode_output = cell["r"] * cell["shard_bytes"]
    decode_input = (cell["k"] - cell["losses"] + cell["r"]) * \
        cell["shard_bytes"]
    decode_output = cell["losses"] * cell["shard_bytes"]
    encode.update({
        "input_GB_per_s": expected_rate(encode_input, encode_median),
        "parity_output_GB_per_s": expected_rate(encode_output, encode_median),
    })
    decode.update({
        "offered_received_GB_per_s": expected_rate(decode_input, decode_median),
        "repaired_output_GB_per_s": expected_rate(decode_output, decode_median),
    })
    document["metrics"] = {
        "codec_setup": codec, "encode_execution": encode,
        "decode_plan_setup": plan, "decode_execution": decode,
        "decode_amortized_at_reuse": {
            "reuse_count": controls["reuse"],
            "derived_median_us_per_batch_call": amortized,
            "offered_received_GB_per_s": expected_rate(decode_input, amortized),
            "repaired_output_GB_per_s": expected_rate(decode_output, amortized),
        },
        "rate_semantics":
            "offered_received counts all non-null shard pointers supplied; "
            "a plan may read a deterministic subset",
    }
    return document


def synthetic_snapshot_identity(
    evidence_schema: str = RAW_SCHEMA,
) -> dict[str, Any]:
    source_root = Path("/fixture/source")
    build = Path("/fixture/build")
    compiler_invocation = "/usr/bin/c++"
    archiver_invocation = "/usr/bin/ar"
    ranlib_invocation = "/usr/bin/ranlib"
    cmake_invocation = "/usr/bin/cmake"

    def artifact(
        path: Path, kind: str, *, payload: bytes | None = None,
        mtime_ns: int = 1, executable: bool = False,
    ) -> dict[str, Any]:
        content = payload if payload is not None else \
            f"{kind}:{path}".encode("utf-8")
        return {
            "path": str(path), "kind": kind, "size": len(content),
            "mode": 0o755 if executable else 0o644,
            "mtime_ns": mtime_ns, "sha256": SUPPORT.sha256_bytes(content),
        }

    compiler = artifact(Path(compiler_invocation), "compiler", executable=True)
    cmake = artifact(Path(cmake_invocation), "build_tool", executable=True)
    specifications = diagnostic_specifications(
        source_root, build, evidence_schema)
    compile_records: list[dict[str, Any]] = []
    compile_entries: list[dict[str, Any]] = []
    for source, output, definitions in specifications:
        compile_definitions = set(definitions)
        if output.parent.name == "leopard_test_hooks.dir":
            compile_definitions = set(HOOK_ARCHIVE_DEFINITIONS)
        elif output.parent.name == "leopard2_backend_avx512_test_hooks.dir":
            compile_definitions.add("LEO2_HAVE_AVX2_BACKEND=1")
        tokens = [compiler_invocation, "-O3", "-DNDEBUG", "-fopenmp",
                  *["-D" + item for item in sorted(compile_definitions)],
                  f"-I{source_root}",
                  "-c", str(source), "-o", str(output)]
        compile_records.append({
            "source": artifact(source, "source_file", mtime_ns=1),
            "object": artifact(output, "object_file", mtime_ns=2),
            "directory": str(build), "command_tokens": tokens,
            "compiler_invocation": compiler_invocation,
            "compiler_resolved_path": compiler["path"],
            "compiler_invocation_identity": copy.deepcopy(compiler),
            "required_definitions": sorted(definitions),
        })
        compile_entries.append({
            "directory": str(build), "file": str(source),
            "output": str(output), "arguments": tokens,
        })
    compile_text = json.dumps(
        compile_entries, sort_keys=True, separators=(",", ":"))
    compile_path = build / "compile_commands.json"

    archive_path = build / "libleopard_test_hooks.a"
    archive_recipe_path = build / "CMakeFiles/leopard_test_hooks.dir/link.txt"
    archive_tokens = [archiver_invocation, "qc", str(archive_path),
                      *[record["object"]["path"] for record in compile_records[:-1]]]
    ranlib_tokens = [ranlib_invocation, str(archive_path)]
    archive_text = " ".join(archive_tokens) + "\n" + \
        " ".join(ranlib_tokens) + "\n"
    archive_artifact = artifact(archive_path, "archive", mtime_ns=3)
    archive_closure = {
        "recipe": artifact(
            archive_recipe_path, "build_metadata",
            payload=archive_text.encode("utf-8")),
        "recipe_text": archive_text,
        "recipe_tokens": [archive_tokens, ranlib_tokens],
        "archiver_invocation": archiver_invocation,
        "ranlib_invocation": ranlib_invocation,
        "archive": archive_artifact,
        "members": [
            {"member": Path(record["object"]["path"]).name,
             "sha256": record["object"]["sha256"]}
            for record in compile_records[:-1]],
        "archiver": artifact(
            Path(archiver_invocation), "archiver", executable=True),
        "ranlib": artifact(
            Path(ranlib_invocation), "ranlib", executable=True),
    }
    archive_closure["archiver_invocation_identity"] = copy.deepcopy(
        archive_closure["archiver"])
    archive_closure["ranlib_invocation_identity"] = copy.deepcopy(
        archive_closure["ranlib"])

    binary_path = build / TARGET_NAME
    benchmark_object = specifications[-1][1]
    system_token = "/usr/lib/libgomp.so.1"
    link_tokens = [compiler_invocation, "-O3", str(benchmark_object),
                   "-o", str(binary_path), str(archive_path), system_token]
    link_text = " ".join(link_tokens) + "\n"
    binary = artifact(
        binary_path, "executable", mtime_ns=4, executable=True)
    system_link = artifact(
        Path(system_token), "system_link_dependency", mtime_ns=1)
    link_closure = {
        "recipe_tokens": link_tokens,
        "compiler_invocation": compiler_invocation,
        "compiler_resolved_path": compiler["path"],
        "compiler_invocation_identity": copy.deepcopy(compiler),
        "output_token": str(binary_path), "output_path": str(binary_path),
        "build_inputs": [
            {"token": str(benchmark_object),
             "resolved_path": str(benchmark_object)},
            {"token": str(archive_path), "resolved_path": str(archive_path)},
        ],
        "system_inputs": [
            {"token": system_token, "resolved_path": system_token}],
    }
    retained_options = {
        "CMAKE_BUILD_TYPE": "Release", "CMAKE_AR": archiver_invocation,
        "CMAKE_CXX_COMPILER": compiler_invocation, "CMAKE_CXX_FLAGS": "",
        "CMAKE_CXX_FLAGS_RELEASE": "-O3 -DNDEBUG",
        "CMAKE_EXPORT_COMPILE_COMMANDS": "ON",
        "CMAKE_RANLIB": ranlib_invocation, "ENABLE_OPENMP": "ON",
        "LEO2_BACKEND_VARIANT": "auto", "LEO2_BUILD_BENCHMARKS": "ON",
        "LEO2_BUILD_FUZZERS": "OFF", "LEO2_BUILD_TESTS": "ON",
        "LEO2_ENABLE_CUDA": "OFF", "LEOPARD_ENABLE_GF8": "ON",
        "LEOPARD_ENABLE_GF16": "ON",
    }
    build_identity = {
        "build_root": str(build),
        "cache": artifact(build / "CMakeCache.txt", "build_metadata"),
        "compile_commands": artifact(
            compile_path, "build_metadata", payload=compile_text.encode("utf-8")),
        "compile_commands_text": compile_text,
        "validated_compile_closure": compile_records,
        "link_recipe": artifact(
            build / "CMakeFiles" / f"{TARGET_NAME}.dir" / "link.txt",
            "build_metadata", payload=link_text.encode("utf-8")),
        "link_recipe_text": link_text,
        "validated_link_closure": link_closure,
        "validated_archive_closure": archive_closure,
        "validated_clean_rebuild": {
            "cmake": copy.deepcopy(cmake),
            "cmake_invocation": cmake_invocation,
            "cmake_resolved_path": cmake["path"],
            "cmake_invocation_identity": copy.deepcopy(cmake),
            "generator": "Unix Makefiles", "retained_options": retained_options,
            "archive_sha256": archive_artifact["sha256"],
            "binary_sha256": binary["sha256"],
        },
        "binary": binary, "hook_archive": copy.deepcopy(archive_artifact),
        "compiler_invocation": compiler_invocation,
        "compiler_invocation_identity": copy.deepcopy(compiler),
        "compiler": compiler,
        "cmake_invocation": cmake_invocation,
        "cmake_invocation_identity": copy.deepcopy(cmake),
        "cmake": cmake,
        "system_link_inputs": [system_link],
        "deterministic_relink_sha256": binary["sha256"],
        "nm": artifact(Path("/usr/bin/nm"), "executable", executable=True),
        "selector_symbols_present_in_hook_archive_and_binary": True,
    }
    runtime_dependencies = [
        {
            "soname": "ld-linux-x86-64.so.2",
            "loader_path": "/lib/ld-linux-x86-64.so.2",
            "file": artifact(
                Path("/lib/ld-linux-x86-64.so.2"), "dynamic_loader"),
        },
        {
            "soname": "libc.so.6", "loader_path": "/usr/lib/libc.so.6",
            "file": artifact(Path("/usr/lib/libc.so.6"), "shared_library"),
        },
        {"soname": "linux-vdso.so.1", "virtual": True},
    ]
    return {
        "source": {
            "path": str(source_root), "head": "a" * 40, "tree": "b" * 40,
            "detached": False, "tracked_tree_listing_sha256": "c" * 64,
            "tracked_status": "clean",
        },
        "build": build_identity,
        "runner": artifact(source_root / RUNNER_RELATIVE, "source_file"),
        "contract": artifact(source_root / CONTRACT_RELATIVE, "source_file"),
        "matrix": artifact(source_root / MATRIX_RELATIVE, "source_file"),
        "support": artifact(
            source_root / "experiments/leopard2/main_compare/run_abba.py",
            "source_file"),
        "runtime": {
            "executable": str(binary_path),
            "dependencies": runtime_dependencies,
        },
    }


def self_test() -> None:
    fixture_build = Path("/fixture/build")
    gfni_output = (
        fixture_build /
        "CMakeFiles/leopard2_backend_gfni_test_hooks.dir" /
        "Leopard2BackendGFNI.cpp.o")
    gfni_entry = {
        "directory": str(fixture_build),
        "output": str(gfni_output.relative_to(fixture_build)),
    }
    require(not compile_entries_include_gfni([], fixture_build) and
            compile_entries_include_gfni([gfni_entry], fixture_build),
            "optional GFNI compile-action detection differs")
    gfni_specs = diagnostic_specifications(
        Path("/fixture/source"), fixture_build, include_gfni=True)
    require(sum(
        output == gfni_output
        for unused_source, output, unused_definitions in gfni_specs) == 1,
        "optional GFNI compile specification differs")
    try:
        compile_entries_include_gfni(
            [gfni_entry, copy.deepcopy(gfni_entry)], fixture_build)
    except EvidenceError:
        pass
    else:
        raise AssertionError(
            "duplicate optional GFNI compile action was accepted")

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
    tail_cell = next(cell for cell in matrix["cells"]
                     if cell["id"] == "gf8-mat-tail")
    controls = {"reuse": 8, "iterations": 3, "warmup": 2}
    result = synthetic_result(tail_cell, controls)
    validate_result(result, tail_cell, "no-copy", controls)
    result_mutations = []
    duplicate_loss = copy.deepcopy(result)
    duplicate_loss["parameters"]["missing_original_indices"][1] = \
        duplicate_loss["parameters"]["missing_original_indices"][0]
    result_mutations.append(duplicate_loss)
    out_of_range_loss = copy.deepcopy(result)
    out_of_range_loss["parameters"]["missing_original_indices"][0] = tail_cell["k"]
    result_mutations.append(out_of_range_loss)
    reordered_losses = copy.deepcopy(result)
    reordered_losses["parameters"]["missing_original_indices"].reverse()
    result_mutations.append(reordered_losses)
    bool_loss = copy.deepcopy(result)
    bool_loss["parameters"]["missing_original_indices"][0] = True
    result_mutations.append(bool_loss)
    wrong_median = copy.deepcopy(result)
    wrong_median["metrics"]["decode_execution"][
        "median_us_per_batch_call"] += 1
    result_mutations.append(wrong_median)
    wrong_sample = copy.deepcopy(result)
    wrong_sample["metrics"]["decode_execution"][
        "samples_us_per_batch_call"][2] += 1
    result_mutations.append(wrong_sample)
    wrong_control = copy.deepcopy(result)
    wrong_control["parameters"]["reuse"] += 1
    result_mutations.append(wrong_control)
    wrong_amortized = copy.deepcopy(result)
    wrong_amortized["metrics"]["decode_amortized_at_reuse"][
        "derived_median_us_per_batch_call"] += 1
    result_mutations.append(wrong_amortized)
    for mutation in result_mutations:
        try:
            validate_result(mutation, tail_cell, "no-copy", controls)
        except (EvidenceError, CONTRACT.ContractError):
            continue
        raise AssertionError("adversarial loss/timing/control mutation was accepted")

    canonical_time = "2026-07-18T00:00:00Z"
    parse_canonical_utc(canonical_time, "fixture created_utc")
    validate_manifest_raw_timestamps(canonical_time, canonical_time)
    for invalid_time in (
            "fixture", "2026-07-18T00:00:00+00:00",
            "2026-07-18T00:00:00", "2026-07-18T00:00:00.000000Z"):
        try:
            parse_canonical_utc(invalid_time, "fixture created_utc")
        except EvidenceError:
            continue
        raise AssertionError("noncanonical raw/manifest timestamp was accepted")
    try:
        validate_manifest_raw_timestamps(
            "2026-07-18T00:00:00Z", "2026-07-18T00:00:01Z")
    except EvidenceError:
        pass
    else:
        raise AssertionError("manifest timestamp before raw evidence was accepted")

    sequence = [
        {"cell": cell, "round": round_index, "slot": slot, "mode": mode}
        for cell, round_index, slot, mode in expected_record_sequence(matrix)]
    validate_record_sequence(sequence, matrix)
    sequence_mutations = []
    all_zero_slots = copy.deepcopy(sequence)
    for record in all_zero_slots:
        record["slot"] = 0
    sequence_mutations.append(all_zero_slots)
    reordered = copy.deepcopy(sequence)
    reordered[0], reordered[1] = reordered[1], reordered[0]
    sequence_mutations.append(reordered)
    wrong_role = copy.deepcopy(sequence)
    wrong_role[0]["mode"] = "no-copy"
    sequence_mutations.append(wrong_role)
    bool_round = copy.deepcopy(sequence)
    bool_round[0]["round"] = False
    sequence_mutations.append(bool_round)
    for mutation in sequence_mutations:
        try:
            validate_record_sequence(mutation, matrix)
        except EvidenceError:
            continue
        raise AssertionError("adversarial ABBA round/slot mutation was accepted")

    with tempfile.TemporaryDirectory(prefix="leo2-high-copy-link-") as root_text:
        root = Path(root_text)
        source = root / "source"
        build = root / "build"
        source.mkdir()
        (build / "obj").mkdir(parents=True)
        benchmark_object = build / "obj/benchmark.cpp.o"
        hook_archive = build / "libleopard_test_hooks.a"
        binary = build / TARGET_NAME
        for path, payload in ((benchmark_object, b"object"),
                              (hook_archive, b"archive"), (binary, b"binary")):
            path.write_bytes(payload)
        compiler = Path("/usr/bin/true").resolve(strict=True)
        tokens = [str(compiler), "obj/benchmark.cpp.o", "-o", TARGET_NAME,
                  "libleopard_test_hooks.a"]
        compiler_identity = SUPPORT.artifact_identity(compiler, "compiler")
        validate_executable_link_inputs(
            tokens, build, source, compiler, str(compiler), compiler_identity,
            binary, benchmark_object, hook_archive)
        impostor = build / "obj/impostor.o"
        impostor.write_bytes(b"impostor")
        tampered = [*tokens, "obj/impostor.o"]
        try:
            validate_executable_link_inputs(
                tampered, build, source, compiler, str(compiler),
                compiler_identity, binary,
                benchmark_object, hook_archive)
        except EvidenceError:
            pass
        else:
            raise AssertionError("undeclared executable link object was accepted")
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
    def host_cpu(cpu: int) -> dict[str, Any]:
        return {
            "cpu": cpu, "cpuinfo": {"model name": "fixture"}, "online": "1",
            "topology": {"thread_siblings_list": "0-1"},
            "frequency_policy": {
                "scaling_driver": "fixture", "scaling_governor": "fixture",
                "energy_performance_preference": "fixture"},
        }

    host = {
        "system": {"release": "fixture"},
        "allowed_cpu_set_at_launch": [0, 1], "online_cpu_set": [0, 1],
        "benchmark_cpu": host_cpu(0), "reserved_sibling": host_cpu(1),
        "turbo_and_pstate": {},
    }
    reservation_payload = {
        "benchmark_cpu": 0, "nonce": "fixture-nonce", "owner": "fixture",
        "reserved_sibling": 1, "schema": SUPPORT.RESERVATION_SCHEMA,
        "status": "held",
    }
    reservation = {
        "path": "/tmp/fixture-reservation",
        "sha256": SUPPORT.sha256_bytes(
            SUPPORT.canonical_bytes(reservation_payload)),
        "payload": reservation_payload, "lock": "exclusive_nonblocking",
    }
    environment_fixture = {
        "host_initial": host, "host_final": copy.deepcopy(host),
        "isolation": nested, "pair_lease": lease, "reservation": reservation,
    }
    campaign_fixture = {"cpu": 0, "reserved_sibling": 1}
    validate_campaign_environment(environment_fixture, campaign_fixture)
    environment_mutations = []
    wrong_host = copy.deepcopy(environment_fixture)
    wrong_host["host_final"]["system"]["release"] = "changed"
    environment_mutations.append(wrong_host)
    wrong_lease = copy.deepcopy(environment_fixture)
    wrong_lease["pair_lease"]["inode"] += 1
    environment_mutations.append(wrong_lease)
    wrong_reservation = copy.deepcopy(environment_fixture)
    wrong_reservation["reservation"]["payload"]["benchmark_cpu"] = 1
    environment_mutations.append(wrong_reservation)
    for mutation in environment_mutations:
        try:
            validate_campaign_environment(mutation, campaign_fixture)
        except EvidenceError:
            continue
        raise AssertionError("adversarial host/lease/reservation mutation was accepted")

    execution_fixture = {
        "started_monotonic_ns": 1_100, "ended_monotonic_ns": 1_200,
        "cpu_before": cpu_stat(0, user=100, idle=100),
        "cpu_after": cpu_stat(0, user=110, idle=110),
        "sibling_before": cpu_stat(1, user=100, idle=100),
        "sibling_after": cpu_stat(1, user=100, idle=110),
    }
    execution_fixture["cpu_delta"] = SUPPORT.cpu_stat_delta(
        execution_fixture["cpu_before"], execution_fixture["cpu_after"])
    execution_fixture["sibling_delta"] = SUPPORT.cpu_stat_delta(
        execution_fixture["sibling_before"], execution_fixture["sibling_after"])
    validate_record_execution_attestation(
        execution_fixture, interval_begin=1_000, interval_end=2_000,
        previous_end=1_000, cpu=0, sibling=1)
    execution_mutations = []
    bool_timestamp = copy.deepcopy(execution_fixture)
    bool_timestamp["started_monotonic_ns"] = True
    execution_mutations.append(bool_timestamp)
    overlapping = copy.deepcopy(execution_fixture)
    overlapping["started_monotonic_ns"] = 1_050
    execution_mutations.append(overlapping)
    wrong_delta = copy.deepcopy(execution_fixture)
    wrong_delta["cpu_delta"]["nonidle_jiffies"] += 1
    execution_mutations.append(wrong_delta)
    for mutation in execution_mutations:
        try:
            validate_record_execution_attestation(
                mutation, interval_begin=1_000, interval_end=2_000,
                previous_end=1_075, cpu=0, sibling=1)
        except EvidenceError:
            continue
        raise AssertionError("adversarial timestamp/CPU mutation was accepted")

    raw_campaign = {
        "cpu": 0, "reserved_sibling": 1, "reuse": controls["reuse"],
        "iterations": controls["iterations"], "warmup": controls["warmup"],
        "thread_count": 1, "timeout_seconds": 120.0,
        "round_orders": matrix["round_orders"],
        "child_environment": dict(SUPPORT.CHILD_ENVIRONMENT),
        "invocation_count": len(expected_record_sequence(matrix)),
    }
    require(validate_manifest_raw_schema_pair(
                LEGACY_MANIFEST_SCHEMA, LEGACY_RAW_SCHEMA) ==
                (LEGACY_MANIFEST_SCHEMA, LEGACY_RAW_SCHEMA) and
            validate_manifest_raw_schema_pair(
                MANIFEST_SCHEMA, RAW_SCHEMA) ==
                (MANIFEST_SCHEMA, RAW_SCHEMA),
            "matching historical/current evidence schemas were rejected")
    for manifest_value, raw_value in (
        (LEGACY_MANIFEST_SCHEMA, RAW_SCHEMA),
        (MANIFEST_SCHEMA, LEGACY_RAW_SCHEMA),
        ("leopard2-high-decode-copy-manifest/unknown", RAW_SCHEMA),
        (MANIFEST_SCHEMA, "leopard2-high-decode-copy-raw/unknown"),
    ):
        try:
            validate_manifest_raw_schema_pair(manifest_value, raw_value)
        except EvidenceError:
            continue
        raise AssertionError(
            "mismatched/unknown manifest and raw schemas were accepted")
    fixture_identity = synthetic_snapshot_identity()
    validate_snapshot_identity(fixture_identity, RAW_SCHEMA)
    legacy_fixture_identity = synthetic_snapshot_identity(LEGACY_RAW_SCHEMA)
    validate_snapshot_identity(legacy_fixture_identity, LEGACY_RAW_SCHEMA)
    current_sources = [
        Path(record["source"]["path"]).name
        for record in fixture_identity["build"]["validated_compile_closure"]
    ]
    legacy_sources = [
        Path(record["source"]["path"]).name
        for record in legacy_fixture_identity["build"]
            ["validated_compile_closure"]
    ]
    expected_current_sources = [
        "leopard.cpp", "leopard2.cpp", "Leopard2Backend.cpp",
        "Leopard2BackendScalar.cpp", "Leopard2CpuFeatures.cpp",
        "Leopard2Plan.cpp", "LeopardCommon.cpp", "LeopardFF8.cpp",
        "LeopardFF16.cpp", "Leopard2BackendSSSE3.cpp",
        "Leopard2BackendAVX2.cpp", "Leopard2BackendAVX2Xor.cpp",
        "Leopard2BackendAVX512.cpp", "benchmark.cpp",
    ]
    require(current_sources == expected_current_sources and
            legacy_sources == [
                name for name in expected_current_sources
                if name != "Leopard2BackendAVX2Xor.cpp"],
            "current/legacy diagnostic source order changed")
    fixture_binary = Path(fixture_identity["build"]["binary"]["path"])
    raw_records = []
    for index, (identifier, round_index, slot, mode) in enumerate(
            expected_record_sequence(matrix)):
        cell = next(item for item in matrix["cells"] if item["id"] == identifier)
        command = benchmark_command(
            fixture_binary, cell, mode, 0, controls["reuse"],
            controls["iterations"], controls["warmup"])
        started = 1_001 + index * 3
        record = {
            "cell": identifier, "mode": mode, "round": round_index, "slot": slot,
            "command": command,
            "command_sha256": SUPPORT.sha256_bytes(
                SUPPORT.canonical_bytes(command)),
            "started_monotonic_ns": started, "ended_monotonic_ns": started + 1,
            "cpu_before": cpu_stat(0, user=100, idle=100),
            "cpu_after": cpu_stat(0, user=110, idle=110),
            "sibling_before": cpu_stat(1, user=100, idle=100),
            "sibling_after": cpu_stat(1, user=100, idle=110),
            "result": synthetic_result(cell, controls, mode),
        }
        record["cpu_delta"] = SUPPORT.cpu_stat_delta(
            record["cpu_before"], record["cpu_after"])
        record["sibling_delta"] = SUPPORT.cpu_stat_delta(
            record["sibling_before"], record["sibling_after"])
        raw_records.append(record)
    raw_fixture_payload = {
        "schema": RAW_SCHEMA, "created_utc": canonical_time,
        "validity_is_independent_of_speed": True, "matrix": matrix,
        "campaign": raw_campaign, "host_initial": host,
        "host_final": copy.deepcopy(host), "reservation": reservation,
        "pair_lease": lease, "isolation": nested,
        "identities_initial": fixture_identity,
        "identities_final": copy.deepcopy(fixture_identity),
        "records": raw_records,
        "analysis": analyze(raw_records, matrix, raw_campaign),
    }
    raw_fixture = SUPPORT.signed(raw_fixture_payload)
    validate_raw(raw_fixture, Path("."), current=False)
    legacy_raw_fixture_payload = copy.deepcopy(raw_fixture_payload)
    legacy_raw_fixture_payload["schema"] = LEGACY_RAW_SCHEMA
    legacy_raw_fixture_payload["identities_initial"] = legacy_fixture_identity
    legacy_raw_fixture_payload["identities_final"] = copy.deepcopy(
        legacy_fixture_identity)
    validate_raw(SUPPORT.signed(legacy_raw_fixture_payload), Path("."),
                 current=False)
    raw_mutations = []
    missing_host = copy.deepcopy(raw_fixture_payload)
    del missing_host["host_final"]
    raw_mutations.append(missing_host)
    extra_stdout_claim = copy.deepcopy(raw_fixture_payload)
    extra_stdout_claim["records"][0]["stdout_sha256"] = "0" * 64
    raw_mutations.append(extra_stdout_claim)
    overlapping_record = copy.deepcopy(raw_fixture_payload)
    overlapping_record["records"][1]["started_monotonic_ns"] = \
        overlapping_record["records"][0]["started_monotonic_ns"]
    raw_mutations.append(overlapping_record)
    bool_control = copy.deepcopy(raw_fixture_payload)
    bool_control["campaign"]["reuse"] = True
    raw_mutations.append(bool_control)
    bad_timestamp = copy.deepcopy(raw_fixture_payload)
    bad_timestamp["created_utc"] = "fixture"
    raw_mutations.append(bad_timestamp)
    oversized_allowed_cpu = copy.deepcopy(raw_fixture_payload)
    for name in ("host_initial", "host_final"):
        oversized_allowed_cpu[name]["allowed_cpu_set_at_launch"].append(
            SUPPORT.MAX_CPU_ID + 1)
    raw_mutations.append(oversized_allowed_cpu)

    def add_identity_mutation(editor: Any) -> None:
        mutation = copy.deepcopy(raw_fixture_payload)
        editor(mutation["identities_initial"])
        mutation["identities_final"] = copy.deepcopy(
            mutation["identities_initial"])
        raw_mutations.append(mutation)

    def substitute_artifact(
        original: Mapping[str, Any], path: str
    ) -> dict[str, Any]:
        result = copy.deepcopy(original)
        payload = (str(result["kind"]) + ":" + path).encode("utf-8")
        result.update({
            "path": path, "size": len(payload),
            "sha256": SUPPORT.sha256_bytes(payload),
        })
        return result

    def bind_text(identity: dict[str, Any], name: str, text: str) -> None:
        identity[name + "_text"] = text
        artifact_value = identity[name]
        encoded = text.encode("utf-8")
        artifact_value["size"] = len(encoded)
        artifact_value["sha256"] = SUPPORT.sha256_bytes(encoded)

    def bind_compile_entries(
        build_value: dict[str, Any], entries: Sequence[Mapping[str, Any]],
    ) -> None:
        text = json.dumps(entries, sort_keys=True, separators=(",", ":"))
        build_value["compile_commands_text"] = text
        encoded = text.encode("utf-8")
        build_value["compile_commands"]["size"] = len(encoded)
        build_value["compile_commands"]["sha256"] = \
            SUPPORT.sha256_bytes(encoded)

    def bind_archive_recipe(closure: dict[str, Any]) -> None:
        text = "\n".join(
            " ".join(command) for command in closure["recipe_tokens"]) + "\n"
        closure["recipe_text"] = text
        encoded = text.encode("utf-8")
        closure["recipe"]["size"] = len(encoded)
        closure["recipe"]["sha256"] = SUPPORT.sha256_bytes(encoded)

    def xor_record_index(identity: Mapping[str, Any]) -> int:
        records = identity["build"]["validated_compile_closure"]
        matches = [
            index for index, record in enumerate(records)
            if Path(record["source"]["path"]).name ==
                "Leopard2BackendAVX2Xor.cpp"
        ]
        require(len(matches) == 1,
                "synthetic current identity lost its AVX2 XOR action")
        return matches[0]

    def omit_xor_action_coherently(identity: dict[str, Any]) -> None:
        build_value = identity["build"]
        index = xor_record_index(identity)
        build_value["validated_compile_closure"].pop(index)
        entries = json.loads(build_value["compile_commands_text"])
        entries.pop(index)
        bind_compile_entries(build_value, entries)
        closure = build_value["validated_archive_closure"]
        closure["recipe_tokens"][0].pop(3 + index)
        closure["members"].pop(index)
        bind_archive_recipe(closure)

    def replace_xor_action_coherently(identity: dict[str, Any]) -> None:
        build_value = identity["build"]
        index = xor_record_index(identity)
        record = build_value["validated_compile_closure"][index]
        replacement_source = "/fixture/source/Leopard2BackendAVX2Replacement.cpp"
        replacement_object = (
            "/fixture/build/CMakeFiles/"
            "leopard2_backend_avx2_test_hooks.dir/"
            "Leopard2BackendAVX2Replacement.cpp.o")
        record["source"] = substitute_artifact(
            record["source"], replacement_source)
        record["object"] = substitute_artifact(
            record["object"], replacement_object)
        source_index = record["command_tokens"].index("-c") + 1
        output_index = record["command_tokens"].index("-o") + 1
        record["command_tokens"][source_index] = replacement_source
        record["command_tokens"][output_index] = replacement_object
        entries = json.loads(build_value["compile_commands_text"])
        entries[index]["file"] = replacement_source
        entries[index]["output"] = replacement_object
        entries[index]["arguments"] = list(record["command_tokens"])
        bind_compile_entries(build_value, entries)
        closure = build_value["validated_archive_closure"]
        closure["recipe_tokens"][0][3 + index] = replacement_object
        closure["members"][index] = {
            "member": Path(replacement_object).name,
            "sha256": record["object"]["sha256"],
        }
        bind_archive_recipe(closure)

    def reorder_xor_archive_recipe(identity: dict[str, Any]) -> None:
        index = xor_record_index(identity)
        require(index > 0, "synthetic AVX2 XOR action cannot be first")
        closure = identity["build"]["validated_archive_closure"]
        archive_tokens = closure["recipe_tokens"][0]
        archive_tokens[3 + index - 1], archive_tokens[3 + index] = \
            archive_tokens[3 + index], archive_tokens[3 + index - 1]
        closure["members"][index - 1], closure["members"][index] = \
            closure["members"][index], closure["members"][index - 1]
        bind_archive_recipe(closure)

    def replace_xor_archive_member(identity: dict[str, Any]) -> None:
        index = xor_record_index(identity)
        identity["build"]["validated_archive_closure"]["members"][index] \
            ["member"] = "Leopard2BackendAVX2Replacement.cpp.o"

    def substitute_compile_action(identity: dict[str, Any]) -> None:
        build_value = identity["build"]
        record = build_value["validated_compile_closure"][0]
        false_path = "/usr/bin/false"
        false_identity = substitute_artifact(
            record["compiler_invocation_identity"], false_path)
        record["command_tokens"][0] = false_path
        record["compiler_invocation"] = false_path
        record["compiler_resolved_path"] = false_path
        record["compiler_invocation_identity"] = false_identity
        entries = json.loads(build_value["compile_commands_text"])
        entries[0]["arguments"][0] = false_path
        text = json.dumps(entries, sort_keys=True, separators=(",", ":"))
        build_value["compile_commands_text"] = text
        encoded = text.encode("utf-8")
        build_value["compile_commands"]["size"] = len(encoded)
        build_value["compile_commands"]["sha256"] = \
            SUPPORT.sha256_bytes(encoded)

    def substitute_link_driver(identity: dict[str, Any]) -> None:
        build_value = identity["build"]
        false_path = "/usr/bin/false"
        false_identity = substitute_artifact(build_value["compiler"], false_path)
        closure = build_value["validated_link_closure"]
        closure["recipe_tokens"][0] = false_path
        closure["compiler_invocation"] = false_path
        closure["compiler_resolved_path"] = false_path
        closure["compiler_invocation_identity"] = false_identity
        bind_text(build_value, "link_recipe",
                  " ".join(closure["recipe_tokens"]) + "\n")
        build_value["validated_clean_rebuild"]["retained_options"] \
            ["CMAKE_CXX_COMPILER"] = false_path

    def substitute_complete_compiler(identity: dict[str, Any]) -> None:
        build_value = identity["build"]
        false_path = "/usr/bin/false"
        false_identity = substitute_artifact(build_value["compiler"], false_path)
        build_value["compiler"] = false_identity
        build_value["compiler_invocation"] = false_path
        build_value["compiler_invocation_identity"] = copy.deepcopy(false_identity)
        entries = json.loads(build_value["compile_commands_text"])
        for record, entry in zip(
                build_value["validated_compile_closure"], entries):
            record["command_tokens"][0] = false_path
            record["compiler_invocation"] = false_path
            record["compiler_resolved_path"] = false_path
            record["compiler_invocation_identity"] = copy.deepcopy(false_identity)
            entry["arguments"][0] = false_path
        compile_text = json.dumps(entries, sort_keys=True, separators=(",", ":"))
        build_value["compile_commands_text"] = compile_text
        encoded = compile_text.encode("utf-8")
        build_value["compile_commands"]["size"] = len(encoded)
        build_value["compile_commands"]["sha256"] = \
            SUPPORT.sha256_bytes(encoded)
        closure = build_value["validated_link_closure"]
        closure["recipe_tokens"][0] = false_path
        closure["compiler_invocation"] = false_path
        closure["compiler_resolved_path"] = false_path
        closure["compiler_invocation_identity"] = copy.deepcopy(false_identity)
        bind_text(build_value, "link_recipe",
                  " ".join(closure["recipe_tokens"]) + "\n")
        build_value["validated_clean_rebuild"]["retained_options"] \
            ["CMAKE_CXX_COMPILER"] = false_path

    def substitute_archiver(identity: dict[str, Any]) -> None:
        build_value = identity["build"]
        closure = build_value["validated_archive_closure"]
        true_path = "/usr/bin/true"
        true_identity = substitute_artifact(closure["archiver"], true_path)
        closure["archiver"] = true_identity
        closure["archiver_invocation_identity"] = copy.deepcopy(true_identity)
        closure["archiver_invocation"] = true_path
        closure["recipe_tokens"][0][0] = true_path
        recipe_text = "\n".join(
            " ".join(command) for command in closure["recipe_tokens"]) + "\n"
        closure["recipe_text"] = recipe_text
        encoded = recipe_text.encode("utf-8")
        closure["recipe"]["size"] = len(encoded)
        closure["recipe"]["sha256"] = SUPPORT.sha256_bytes(encoded)
        build_value["validated_clean_rebuild"]["retained_options"] \
            ["CMAKE_AR"] = true_path

    def substitute_cmake(identity: dict[str, Any]) -> None:
        build_value = identity["build"]
        true_path = "/usr/bin/true"
        true_identity = substitute_artifact(build_value["cmake"], true_path)
        build_value["cmake"] = true_identity
        build_value["cmake_invocation"] = true_path
        build_value["cmake_invocation_identity"] = copy.deepcopy(true_identity)
        clean = build_value["validated_clean_rebuild"]
        clean["cmake"] = copy.deepcopy(true_identity)
        clean["cmake_invocation"] = true_path
        clean["cmake_resolved_path"] = true_path
        clean["cmake_invocation_identity"] = copy.deepcopy(true_identity)

    add_identity_mutation(lambda identity: identity["source"].pop("tree"))
    add_identity_mutation(lambda identity: identity["build"].pop("cache"))
    add_identity_mutation(lambda identity: identity["build"]
                          ["validated_compile_closure"][0]["source"]
                          .__setitem__("path", "/fixture/source/impostor.cpp"))
    add_identity_mutation(lambda identity: identity["build"]
                          ["validated_compile_closure"][0]["object"]
                          .__setitem__("path", "/fixture/build/impostor.o"))
    add_identity_mutation(lambda identity: identity["build"]
                          ["validated_compile_closure"][0]["command_tokens"]
                          .append("-DUNRETAINED=1"))
    add_identity_mutation(omit_xor_action_coherently)
    add_identity_mutation(replace_xor_action_coherently)
    add_identity_mutation(substitute_compile_action)
    add_identity_mutation(substitute_complete_compiler)
    add_identity_mutation(lambda identity: identity["build"]
                          ["validated_archive_closure"]
                          .__setitem__("recipe_text", "substituted\n"))
    add_identity_mutation(lambda identity: identity["build"]
                          ["validated_archive_closure"]["members"][0]
                          .__setitem__("sha256", "e" * 64))
    add_identity_mutation(reorder_xor_archive_recipe)
    add_identity_mutation(replace_xor_archive_member)
    add_identity_mutation(lambda identity: identity["build"]
                          ["validated_archive_closure"]["archive"]
                          .__setitem__("sha256", "f" * 64))
    add_identity_mutation(lambda identity: identity["build"]
                          ["validated_archive_closure"]
                          .__setitem__("archiver_invocation", "/usr/bin/true"))
    add_identity_mutation(substitute_archiver)
    add_identity_mutation(lambda identity: identity["build"]
                          .__setitem__("link_recipe_text", "substituted\n"))
    add_identity_mutation(substitute_link_driver)
    add_identity_mutation(lambda identity: identity["build"]
                          ["validated_link_closure"]["build_inputs"][0]
                          .__setitem__("resolved_path", "/fixture/build/impostor.o"))
    add_identity_mutation(lambda identity: identity["build"]
                          ["system_link_inputs"][0]
                          .__setitem__("path", "/usr/lib/libomp.so.1"))
    add_identity_mutation(lambda identity: identity["build"]["binary"]
                          .__setitem__("sha256", "2" * 64))
    add_identity_mutation(lambda identity: identity["build"]
                          .__setitem__("deterministic_relink_sha256", "3" * 64))
    add_identity_mutation(lambda identity: identity["build"]
                          ["validated_clean_rebuild"]["retained_options"]
                          .__setitem__("CMAKE_BUILD_TYPE", "Debug"))
    add_identity_mutation(lambda identity: identity["build"]
                          ["validated_clean_rebuild"]
                          .__setitem__("generator", "Ninja"))
    add_identity_mutation(lambda identity: identity["build"]["compiler"]
                          .__setitem__("path", "/usr/bin/false"))
    add_identity_mutation(lambda identity: identity["build"]["nm"]
                          .__setitem__("path", "/usr/bin/false"))
    add_identity_mutation(substitute_cmake)
    for component in ("runner", "contract", "matrix", "support"):
        add_identity_mutation(
            lambda identity, name=component: identity[name].__setitem__(
                "path", f"/fixture/source/substituted-{name}"))
    add_identity_mutation(lambda identity: identity["runtime"]
                          ["dependencies"][1]["file"]
                          .__setitem__("path", "/usr/lib/libm.so.6"))
    add_identity_mutation(lambda identity: identity["runtime"]
                          ["dependencies"][1]
                          .__setitem__("loader_path", "/lib/libc.so.6"))
    for mutation in raw_mutations:
        try:
            validate_raw(SUPPORT.signed(mutation), Path("."), current=False)
        except (EvidenceError, CONTRACT.ContractError):
            continue
        raise AssertionError("adversarial integrated raw mutation was accepted")
    flat = {"accepted": True, "reserved_sibling_nonidle_jiffies": 0}
    try:
        validate_campaign_isolation(flat, {"cpu": 0, "reserved_sibling": 1})
    except SUPPORT.EvidenceError:
        pass
    else:
        raise AssertionError("obsolete flat isolation fixture was accepted")
    CONTRACT.self_test()
    print("high-decode copy/no-copy runner self-test passed: "
          "16 cells, 192 invocations, deterministic losses/timings, exact ABBA "
          "slots, current v3 AVX2-XOR closure, historical v2 replay, canonical "
          "build inputs, field, and nested isolation")


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
