#!/usr/bin/env python3
"""Run a bounded, fail-closed ABBA screen for small GF8 direct repair.

This runner compares two already-built Leopard2 benchmark executables whose
only semantic difference is an explicitly selected small-code direct-repair
mode: transform, output-major, or source-major.  Its separately versioned XMM
tail comparison keeps source-major mode 2 identical and permits only the AVX2
backend's explicit 0-versus-1 tail definition to change.  It is not for
comparison to the historical Leopard executable (use
../main_compare/run_abba.py for that).
It owns the canonical GF8 and physical-pair leases for the complete campaign,
moves all same-user threads off the held pair, pins every child to one logical
CPU, requires its SMT sibling to remain idle, retains raw output, and verifies
wire digests and the selected direct path before accepting timing.

SIGKILL cannot run userspace cleanup: kernel locks close, but a campaign killed
with SIGKILL may leave coordinator task affinities changed.  HUP, INT, and TERM
are deferred until child/session, affinity, pair-lease, and global-lock cleanup.
"""

from __future__ import annotations

import argparse
import datetime as dt
import fcntl
import hashlib
import importlib.util
import json
import math
import os
import platform
import re
import secrets
import signal
import shlex
import shutil
import stat
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Sequence


SCHEMA = "leopard2-small-direct-abba/v1"
XMM_TAIL_SCHEMA = "leopard2-small-direct-xmm-tail-abba/v1"
MATRIX_SCHEMA = "leopard2-small-direct-matrix/v1"
LARGE_MATRIX_SCHEMA = "leopard2-small-direct-full-matrix/v1"
TINY_MATRIX_SCHEMA = "leopard2-small-direct-tiny-matrix/v1"
XMM_TAIL_MATRIX_SCHEMA = "leopard2-small-direct-xmm-tail-matrix/v1"
CELL_SCHEMA = "leopard2-small-direct-cell/v1"
ENVELOPE_SCHEMA = "leopard2-small-direct-envelope/v1"
TILE_SPEC_SCHEMA = "leopard2-direct-tile-lab-spec/v2"
TILE_SCREEN_SCHEMA = "leopard2-direct-tile-screen/v2"
LAB_RESULT_SCHEMA = "leopard2-lab-result/v2"
BINARY_FILE_IDENTITY_KEYS = (
    "runner", "baseline_executable", "baseline_archive",
    "baseline_avx2_object", "baseline_compile_commands",
    "baseline_avx2_xor_object",
    "baseline_core_object", "baseline_avx512_object", "pair_lease_source",
    "candidate_executable", "candidate_archive", "candidate_avx2_object",
    "candidate_compile_commands", "candidate_core_object",
    "candidate_avx512_object", "candidate_avx2_xor_object",
)
ORDER = ("baseline", "candidate", "candidate", "baseline")
DEFAULT_LOCK = Path("/tmp/leopard-gf8-authoritative.lock")
DEFAULT_ROUNDS = 3
CONTROL_SIGNALS = frozenset((signal.SIGHUP, signal.SIGINT, signal.SIGTERM))
# Target CPU and stopped-child schedstat use the same nanosecond accounting
# domain.  A fixed five-microsecond endpoint allowance covers the scheduler's
# two stopped/zombie snapshot boundaries while rejecting any material foreign
# task.  The independent microsecond-resolution wait4 cross-check gets a
# separate fixed allowance and is never used to hide schedstat excess.
TARGET_RUNTIME_TOLERANCE_NS = 5_000
RUSAGE_CROSSCHECK_TOLERANCE_NS = 5_000
T_CRITICAL_95 = {2: 4.302652729911275, 4: 2.7764451051977987}
CHILD_ENV = {
    "LANG": "C",
    "LC_ALL": "C",
    "OMP_DYNAMIC": "FALSE",
    "OMP_NUM_THREADS": "1",
    "OMP_PLACES": "cores",
    "OMP_PROC_BIND": "TRUE",
    "PATH": "/usr/bin:/bin",
    "TZ": "UTC",
}
SOURCE_FILES = (
    "CMakeLists.txt",
    "leopard.cpp", "leopard.h", "leopard2.cpp", "leopard2.h",
    "Leopard2Backend.cpp", "Leopard2Backend.h",
    "Leopard2BackendScalar.cpp", "Leopard2BackendSSSE3.cpp",
    "Leopard2BackendAVX2.cpp", "Leopard2BackendAVX2Xor.cpp",
    "Leopard2BackendAVX512.cpp", "Leopard2CpuFeatures.cpp",
    "Leopard2Direct.h", "Leopard2Dispatch.h", "Leopard2Plan.cpp",
    "Leopard2Plan.h",
    "LeopardCommon.cpp", "LeopardCommon.h",
    "LeopardFF8.cpp", "LeopardFF8.h", "LeopardFF16.cpp",
    "LeopardFF16.h", "bench/leopard2/benchmark.cpp",
)
ARCHIVE_MEMBER_SPECS = (
    ("leopard.cpp.o", "leopard.cpp", "leopard.dir"),
    ("leopard2.cpp.o", "leopard2.cpp", "leopard.dir"),
    ("Leopard2Backend.cpp.o", "Leopard2Backend.cpp", "leopard.dir"),
    ("Leopard2BackendScalar.cpp.o", "Leopard2BackendScalar.cpp",
     "leopard.dir"),
    ("Leopard2CpuFeatures.cpp.o", "Leopard2CpuFeatures.cpp", "leopard.dir"),
    ("Leopard2Plan.cpp.o", "Leopard2Plan.cpp", "leopard.dir"),
    ("LeopardCommon.cpp.o", "LeopardCommon.cpp", "leopard.dir"),
    ("LeopardFF8.cpp.o", "LeopardFF8.cpp", "leopard.dir"),
    ("LeopardFF16.cpp.o", "LeopardFF16.cpp", "leopard.dir"),
    ("Leopard2BackendSSSE3.cpp.o", "Leopard2BackendSSSE3.cpp",
     "leopard2_backend_ssse3.dir"),
    ("Leopard2BackendAVX2.cpp.o", "Leopard2BackendAVX2.cpp",
     "leopard2_backend_avx2.dir"),
    ("Leopard2BackendAVX2Xor.cpp.o", "Leopard2BackendAVX2Xor.cpp",
     "leopard2_backend_avx2.dir"),
    ("Leopard2BackendAVX512.cpp.o", "Leopard2BackendAVX512.cpp",
     "leopard2_backend_avx512.dir"),
)
EXPECTED_ARCHIVE_MEMBERS = tuple(item[0] for item in ARCHIVE_MEMBER_SPECS)
MODE_COMPILE_DEFINITIONS = {
    "transform": "-DLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=0",
    "output": "-DLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=1",
    "source": "-DLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=2",
    "source-xmm0": "-DLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=2",
    "source-xmm1": "-DLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=2",
    # The production binary deliberately has no diagnostic definition.  Its
    # default is checked in leopard2.cpp before a campaign can start.
    "production": None,
}
XMM_TAIL_COMPILE_ARGUMENTS = {
    "source-xmm0": (
        "-DLEO2_EXPERIMENT_GF8_DIRECT_XMM_KAT=1",
        "-DLEO2_EXPERIMENT_GF8_DIRECT_XMM_TAIL=0",
    ),
    "source-xmm1": (
        "-DLEO2_EXPERIMENT_GF8_DIRECT_XMM_KAT=1",
        "-DLEO2_EXPERIMENT_GF8_DIRECT_XMM_TAIL=1",
    ),
}
DIAGNOSTIC_DEFINITION_NAMES = (
    "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE",
    "LEO2_EXPERIMENT_GF8_DIRECT_XMM_KAT",
    "LEO2_EXPERIMENT_GF8_DIRECT_XMM_TAIL",
)
XMM_TAIL_CONTROL_MODE = "source-xmm0"
XMM_TAIL_CANDIDATE_MODE = "source-xmm1"
XMM_TAIL_AUDIT_SCHEMA = "leopard2-avx2-xmm-tail-encoding-audit/v1"
XMM_TAIL_COMPARISON_AUDIT_SCHEMA = \
    "leopard2-avx2-xmm-tail-comparison-audit/v1"
PRODUCTION_TARGETS = {
    source: target for unused_member, source, target in ARCHIVE_MEMBER_SPECS
}
PRODUCTION_TARGETS["benchmark.cpp"] = "bench_leopard2.dir"


class EvidenceError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvidenceError(message)


def mode_compile_definition(mode: str) -> str | None:
    require(mode in MODE_COMPILE_DEFINITIONS,
            "unknown small-direct comparison mode")
    return MODE_COMPILE_DEFINITIONS[mode]


def mode_compile_arguments(mode: str) -> list[str]:
    definition = mode_compile_definition(mode)
    result = [] if definition is None else [definition]
    result.extend(XMM_TAIL_COMPILE_ARGUMENTS.get(mode, ()))
    return result


def validate_production_mode_source(source_root: Path) -> None:
    source = (source_root / "leopard2.cpp").read_text()
    pattern = re.compile(
        r"#ifndef\s+LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE\s+"
        r"#define\s+LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE\s+3\s+"
        r"#endif")
    require(len(pattern.findall(source)) == 1,
            "production source does not bind small-direct mode 3 exactly")


def comparison_modes(baseline: str, candidate: str) -> dict[str, str]:
    require(baseline in MODE_COMPILE_DEFINITIONS and
            candidate in MODE_COMPILE_DEFINITIONS and
            baseline != candidate,
            "baseline and candidate modes must be distinct and known")
    result = {"baseline": baseline, "candidate": candidate}
    if baseline in XMM_TAIL_COMPILE_ARGUMENTS or \
            candidate in XMM_TAIL_COMPILE_ARGUMENTS:
        require(baseline == XMM_TAIL_CONTROL_MODE and
                candidate == XMM_TAIL_CANDIDATE_MODE,
                "XMM-tail comparison requires source-xmm0 baseline and "
                "source-xmm1 candidate")
    return result


def is_xmm_tail_comparison(modes: dict[str, str]) -> bool:
    return modes == {
        "baseline": XMM_TAIL_CONTROL_MODE,
        "candidate": XMM_TAIL_CANDIDATE_MODE,
    }


def comparison_schema(modes: dict[str, str]) -> str:
    comparison_modes(modes["baseline"], modes["candidate"])
    return XMM_TAIL_SCHEMA if is_xmm_tail_comparison(modes) else SCHEMA


def require_comparison_schema(
        schema: Any, modes: dict[str, str], label: str) -> None:
    require(schema == comparison_schema(modes),
            "%s schema does not match comparison modes" % label)


def strip_mode_definition(arguments: list[str], mode: str,
                          label: str) -> list[str]:
    expected = mode_compile_arguments(mode)
    result = list(arguments)
    def is_diagnostic_definition(argument: str) -> bool:
        return any(argument == "-D" + name or
                   argument.startswith("-D" + name + "=")
                   for name in DIAGNOSTIC_DEFINITION_NAMES)
    if not expected:
        require(not any(is_diagnostic_definition(value) for value in result),
                "%s production build contains a diagnostic mode" % label)
        return result
    for definition in expected:
        require(result.count(definition) == 1,
                "%s must contain its exact mode definition once" % label)
        result.remove(definition)
    require(not any(is_diagnostic_definition(value) for value in result),
            "%s contains another mode definition" % label)
    return result


def comparison_changed_members(modes: dict[str, str]) -> set[str]:
    comparison_modes(modes["baseline"], modes["candidate"])
    if is_xmm_tail_comparison(modes):
        return {"Leopard2BackendAVX2.cpp.o"}
    return {"leopard2.cpp.o"}


def validate_comparison_object_deltas(
        modes: dict[str, str],
        baseline: dict[str, str],
        candidate: dict[str, str]) -> None:
    require(set(baseline) == set(EXPECTED_ARCHIVE_MEMBERS) and
            set(candidate) == set(EXPECTED_ARCHIVE_MEMBERS),
            "comparison object map does not cover the exact archive")
    changed = comparison_changed_members(modes)
    for member in EXPECTED_ARCHIVE_MEMBERS:
        if member in changed:
            require(baseline[member] != candidate[member],
                    "comparison-changing object is byte-identical: %s" %
                    member)
        else:
            require(baseline[member] == candidate[member],
                    "unexpected cross-build object difference: %s" % member)


def isolation_policy() -> dict[str, Any]:
    return {
        "same_user_pair_affinity": "empty before and after invocation",
        "target_runtime_accounting":
            "target CPU schedstat minus stopped-child schedstat delta",
        "target_runtime_tolerance_ns": TARGET_RUNTIME_TOLERANCE_NS,
        "wait4_crosscheck_tolerance_ns":
            RUSAGE_CROSSCHECK_TOLERANCE_NS,
        "target_interrupt_fields_rejected": [
            "irq", "softirq", "steal", "guest", "guest_nice"],
        "child_scope": (
            "forked session leader pins and SIGSTOPs before snapshots, "
            "then directly execs benchmark; waitid WNOWAIT preserves "
            "zombie schedstat and wait4 aggregates all child threads; "
            "HUP/INT/TERM remain blocked in the parent until exact cleanup"),
        "sigkill_limitation": (
            "SIGKILL cannot run userspace affinity restoration; kernel-held "
            "locks still close on process death"),
    }


def require_no_pending_control_signal() -> None:
    pending = set(signal.sigpending()) & CONTROL_SIGNALS
    require(not pending,
            "campaign interrupted by pending control signal(s): %s" %
            sorted(item.name for item in pending))


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def file_identity(path: Path) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    require(resolved.is_file(), "%s is not a regular file" % resolved)
    return {
        "path": str(resolved),
        "size": resolved.stat().st_size,
        "sha256": sha256_file(resolved),
    }


def summarize_avx2_disassembly(output: bytes) -> dict[str, Any]:
    instruction_count = 0
    vex_instruction_count = 0
    evex_instruction_count = 0
    target_instruction_count = 0
    target_stack_vector_access_count = 0
    target_xmm_instruction_count = 0
    target_xmm_vpshufb_count = 0
    target_xmm_memory_instruction_count = 0
    target_legacy_xmm_instruction_count = 0
    per_function = {
        "pair": {"instruction": 0, "stack_vector_access": 0,
                 "xmm_instruction": 0, "xmm_vpshufb": 0,
                 "xmm_memory_instruction": 0, "legacy_xmm_instruction": 0},
        "group4": {"instruction": 0, "stack_vector_access": 0,
                   "xmm_instruction": 0, "xmm_vpshufb": 0,
                   "xmm_memory_instruction": 0,
                   "legacy_xmm_instruction": 0},
    }
    target_function = None
    pattern = re.compile(
        r"^\s*[0-9a-fA-F]+:\s+"
        r"((?:[0-9a-fA-F]{2}(?:\s+|$))+)\s*(.*)$")
    for line in output.decode("utf-8", errors="replace").splitlines():
        if re.match(r"^\s*[0-9a-fA-F]+\s+<.*>:$", line):
            if "AVX2FF8MultiplyAddOutputPair" in line:
                target_function = "pair"
            elif "AVX2FF8MultiplyAddOutputGroup4" in line:
                target_function = "group4"
            else:
                target_function = None
            continue
        match = pattern.match(line)
        if match is None:
            continue
        encoding = match.group(1).split()
        if not encoding:
            continue
        first = int(encoding[0], 16)
        instruction_count += 1
        if first in (0xc4, 0xc5):
            vex_instruction_count += 1
        elif first == 0x62:
            evex_instruction_count += 1
        if target_function:
            counts = per_function[target_function]
            target_instruction_count += 1
            counts["instruction"] += 1
            assembly = match.group(2)
            has_xmm = "%xmm" in assembly
            if has_xmm:
                target_xmm_instruction_count += 1
                counts["xmm_instruction"] += 1
                if first not in (0xc4, 0xc5, 0x62):
                    target_legacy_xmm_instruction_count += 1
                    counts["legacy_xmm_instruction"] += 1
                if assembly.lstrip().startswith("vpshufb"):
                    target_xmm_vpshufb_count += 1
                    counts["xmm_vpshufb"] += 1
                if assembly.lstrip().startswith("vmov") and "(" in assembly:
                    target_xmm_memory_instruction_count += 1
                    counts["xmm_memory_instruction"] += 1
            if assembly.lstrip().startswith("v") and re.search(
                    r"\(%r(?:sp|bp)\)", assembly):
                target_stack_vector_access_count += 1
                counts["stack_vector_access"] += 1
    require(instruction_count > 0,
            "AVX2 object disassembly contains no instructions")
    require(evex_instruction_count == 0,
            "AVX2 object disassembly contains EVEX instructions")
    require(vex_instruction_count > 0,
            "AVX2 object disassembly contains no VEX instructions")
    return {
        "output_sha256": sha256_bytes(output),
        "instruction_count": instruction_count,
        "vex_instruction_count": vex_instruction_count,
        "evex_instruction_count": evex_instruction_count,
        "target_instruction_count": target_instruction_count,
        "target_stack_vector_access_count":
            target_stack_vector_access_count,
        "target_xmm_instruction_count": target_xmm_instruction_count,
        "target_xmm_vpshufb_count": target_xmm_vpshufb_count,
        "target_xmm_memory_instruction_count":
            target_xmm_memory_instruction_count,
        "target_legacy_xmm_instruction_count":
            target_legacy_xmm_instruction_count,
        **{
            "target_%s_%s_count" % (function, counter): value
            for function, counters in per_function.items()
            for counter, value in counters.items()
        },
    }


def summarize_text_sections(output: bytes) -> dict[str, Any]:
    text_bytes = 0
    pattern = re.compile(
        r"^\s*\d+\s+(\S+)\s+([0-9a-fA-F]+)(?:\s+|$)")
    for line in output.decode("utf-8", errors="replace").splitlines():
        match = pattern.match(line)
        if match is not None and match.group(1).startswith(".text"):
            text_bytes += int(match.group(2), 16)
    require(text_bytes > 0, "AVX2 object has no executable text section")
    return {
        "section_output_sha256": sha256_bytes(output),
        "text_bytes": text_bytes,
    }


def avx2_encoding_audit(path: Path) -> dict[str, Any]:
    executable = shutil.which("objdump", path="/usr/bin:/bin")
    require(executable is not None,
            "objdump is required for the XMM-tail encoding audit")
    tool = Path(executable).resolve(strict=True)
    object_path = path.resolve(strict=True)
    disassembly = subprocess.run(
        [str(tool), "-d", "-w", "-C", str(object_path)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=CHILD_ENV,
        check=False)
    require(disassembly.returncode == 0 and not disassembly.stderr,
            "cannot disassemble AVX2 object for XMM-tail audit")
    sections = subprocess.run(
        [str(tool), "-h", "-w", str(object_path)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=CHILD_ENV,
        check=False)
    require(sections.returncode == 0 and not sections.stderr,
            "cannot inspect AVX2 text sections for XMM-tail audit")
    result = {
        "schema": XMM_TAIL_AUDIT_SCHEMA,
        "tool": file_identity(tool),
        "object": file_identity(object_path),
        **summarize_avx2_disassembly(disassembly.stdout),
        **summarize_text_sections(sections.stdout),
    }
    result["digest"] = sha256_bytes(canonical_bytes(result))
    return result


def validate_avx2_encoding_audit(
        value: Any, expected_object: dict[str, Any], label: str) -> None:
    per_function_keys = {
        "target_%s_%s_count" % (function, counter)
        for function in ("pair", "group4")
        for counter in (
            "instruction", "stack_vector_access", "xmm_instruction",
            "xmm_vpshufb", "xmm_memory_instruction",
            "legacy_xmm_instruction")
    }
    require_exact_keys(value, {
        "schema", "tool", "object", "output_sha256", "instruction_count",
        "vex_instruction_count", "evex_instruction_count", "digest",
        "target_instruction_count", "target_stack_vector_access_count",
        "target_xmm_instruction_count", "target_xmm_vpshufb_count",
        "target_xmm_memory_instruction_count",
        "target_legacy_xmm_instruction_count",
        "section_output_sha256", "text_bytes",
    } | per_function_keys, label)
    validate_file_identity_record(value["tool"], label + " tool")
    validate_file_identity_record(value["object"], label + " object")
    require(value["schema"] == XMM_TAIL_AUDIT_SCHEMA and
            value["object"] == expected_object and
            isinstance(value["output_sha256"], str) and
            re.fullmatch(r"[0-9a-f]{64}", value["output_sha256"]) is not None and
            isinstance(value["section_output_sha256"], str) and
            re.fullmatch(r"[0-9a-f]{64}",
                         value["section_output_sha256"]) is not None and
            type(value["instruction_count"]) is int and
            type(value["vex_instruction_count"]) is int and
            type(value["evex_instruction_count"]) is int and
            type(value["target_instruction_count"]) is int and
            type(value["target_stack_vector_access_count"]) is int and
            type(value["target_xmm_instruction_count"]) is int and
            type(value["target_xmm_vpshufb_count"]) is int and
            type(value["target_xmm_memory_instruction_count"]) is int and
            type(value["target_legacy_xmm_instruction_count"]) is int and
            all(type(value[key]) is int and value[key] >= 0
                for key in per_function_keys) and
            type(value["text_bytes"]) is int and
            value["instruction_count"] >= value["vex_instruction_count"] > 0 and
            value["evex_instruction_count"] == 0,
            label + " contents are invalid")
    require(value["target_instruction_count"] > 0 and
            value["target_stack_vector_access_count"] == 0 and
            value["target_legacy_xmm_instruction_count"] == 0 and
            all(value["target_%s_instruction_count" % function] > 0 and
                value["target_%s_stack_vector_access_count" % function] == 0
                and value["target_%s_legacy_xmm_instruction_count" %
                          function] == 0
                for function in ("pair", "group4")) and
            value["text_bytes"] > 0,
            label + " target kernel has spills or missing text")
    copy = dict(value)
    digest = copy.pop("digest")
    require(digest == sha256_bytes(canonical_bytes(copy)),
            label + " digest is invalid")


def validate_xmm_tail_audit_pair(
        baseline: dict[str, Any], candidate: dict[str, Any]) -> None:
    delta = candidate["text_bytes"] - baseline["text_bytes"]
    require(0 < delta <= 4096,
            "XMM-tail AVX2 text growth is outside the 1..4096-byte cap")
    for function in ("pair", "group4"):
        prefix = "target_%s_" % function
        require(candidate[prefix + "xmm_instruction_count"] >
                    baseline[prefix + "xmm_instruction_count"] and
                candidate[prefix + "xmm_vpshufb_count"] >
                    baseline[prefix + "xmm_vpshufb_count"] and
                candidate[prefix + "xmm_memory_instruction_count"] >
                    baseline[prefix + "xmm_memory_instruction_count"],
                "candidate lacks an attributable %s XMM shuffle/load/store "
                "delta" % function)


def xmm_tail_comparison_audit(
        baseline: dict[str, Any], candidate: dict[str, Any]) -> dict[str, Any]:
    """Bind the exact code-size and per-kernel instruction deltas."""
    validate_xmm_tail_audit_pair(baseline, candidate)
    result = {
        "schema": XMM_TAIL_COMPARISON_AUDIT_SCHEMA,
        "baseline_audit_digest": baseline["digest"],
        "candidate_audit_digest": candidate["digest"],
        "text_delta_bytes": candidate["text_bytes"] - baseline["text_bytes"],
    }
    for function in ("pair", "group4"):
        for counter in (
                "instruction", "xmm_instruction", "xmm_vpshufb",
                "xmm_memory_instruction"):
            key = "target_%s_%s_count" % (function, counter)
            result["target_%s_%s_delta" % (function, counter)] = \
                candidate[key] - baseline[key]
    require(all(result["target_%s_%s_delta" % (function, counter)] > 0
                for function in ("pair", "group4")
                for counter in (
                    "instruction", "xmm_instruction", "xmm_vpshufb",
                    "xmm_memory_instruction")),
            "candidate per-kernel XMM code deltas are not all positive")
    result["digest"] = sha256_bytes(canonical_bytes(result))
    return result


def validate_xmm_tail_comparison_audit(
        value: Any, baseline: dict[str, Any],
        candidate: dict[str, Any]) -> None:
    expected = xmm_tail_comparison_audit(baseline, candidate)
    require(value == expected,
            "XMM-tail comparison audit does not reconstruct exactly")


def archive_members_identity(
        path: Path, archiver: Path = Path("/usr/bin/ar")
) -> dict[str, Any]:
    archive = path.resolve(strict=True)
    ar = archiver.resolve(strict=True)
    completed = subprocess.run(
        [str(ar), "t", str(archive)], stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, check=False, env=CHILD_ENV)
    require(completed.returncode == 0,
            "cannot enumerate archive members: %s" %
            completed.stderr.decode("utf-8", errors="replace"))
    members = completed.stdout.decode("utf-8", errors="strict").splitlines()
    require(members == list(EXPECTED_ARCHIVE_MEMBERS),
            "archive member closure or exact order is invalid")
    member_files = {}
    for member in EXPECTED_ARCHIVE_MEMBERS:
        extracted = subprocess.run(
            [str(ar), "p", str(archive), member], stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False, env=CHILD_ENV)
        require(extracted.returncode == 0,
                "cannot extract archive member %s: %s" % (
                    member, extracted.stderr.decode(
                        "utf-8", errors="replace")))
        member_files[member] = {
            "size": len(extracted.stdout),
            "sha256": sha256_bytes(extracted.stdout),
        }
    result = {
        "schema": "leopard2-archive-members/v2",
        "path": str(archive), "members": members,
        "member_files": member_files,
        "archiver": file_identity(ar),
    }
    result["sha256"] = sha256_bytes(canonical_bytes(result))
    return result


def compile_entry_for(path: Path, filename: str) -> dict[str, Any]:
    compile_commands = path.resolve(strict=True)
    value = json.loads(compile_commands.read_text())
    require(isinstance(value, list), "compile_commands.json is not an array")
    entries = [entry for entry in value if isinstance(entry, dict) and
               Path(str(entry.get("file", ""))).name ==
               filename]
    target = PRODUCTION_TARGETS.get(filename)
    if target:
        marker = "CMakeFiles/%s/" % target
        entries = [entry for entry in entries if
                   marker in str(entry.get("output", "")) or
                   marker in str(entry.get("command", "")) or
                   any(marker in argument for argument in
                       entry.get("arguments", []) if
                       isinstance(argument, str))]
    require(len(entries) == 1,
            ("compile_commands.json must identify exactly one production "
             "%s entry") % filename)
    entry = entries[0]
    require(isinstance(entry.get("command") or entry.get("arguments"),
                       (str, list)),
            "%s compile entry has no command" % filename)
    return entry


def avx2_compile_entry(path: Path) -> dict[str, Any]:
    return compile_entry_for(path, "Leopard2BackendAVX2.cpp")


def compile_entry_source(entry: dict[str, Any]) -> Path:
    source = Path(entry["file"])
    if not source.is_absolute():
        source = Path(entry.get("directory", ".")) / source
    return source.resolve()


def compile_entry_arguments(entry: dict[str, Any]) -> list[str]:
    arguments = entry.get("arguments")
    if isinstance(arguments, list) and all(
            isinstance(argument, str) for argument in arguments):
        return list(arguments)
    command = entry.get("command")
    require(isinstance(command, str) and command,
            "compile entry has no usable argv")
    return shlex.split(command)


def compile_entry_output(entry: dict[str, Any]) -> Path:
    raw_output = entry.get("output")
    if not isinstance(raw_output, str) or not raw_output:
        arguments = compile_entry_arguments(entry)
        try:
            raw_output = arguments[arguments.index("-o") + 1]
        except (ValueError, IndexError):
            raise EvidenceError("compile entry has no identifiable output")
    output = Path(raw_output)
    if not output.is_absolute():
        directory = entry.get("directory")
        require(isinstance(directory, str) and directory,
                "relative compile output has no directory")
        output = Path(directory) / output
    return output.resolve(strict=True)


def validate_compile_profile(arguments: Sequence[str], target: str) -> None:
    optimization_flags = [argument for argument in arguments
                          if re.fullmatch(r"-O(?:[0-3]|fast|g|s|z)",
                                          argument)]
    standard_flags = [argument for argument in arguments
                      if argument.startswith("-std=")]
    forbidden_prefixes = (
        "-flto", "-fprofile", "-fsanitize", "-march=", "-mcpu=",
        "-mtune=",
    )
    require(optimization_flags and
            all(argument == "-O3" for argument in optimization_flags) and
            standard_flags == ["-std=gnu++11"] and
            not any(argument.startswith(forbidden_prefixes)
                    for argument in arguments) and
            not any(argument in ("--coverage", "-fno-openmp")
                    for argument in arguments) and
            "-O3" in arguments and "-DNDEBUG" in arguments and
            "-Wall" in arguments and "-Wextra" in arguments and
            "-std=gnu++11" in arguments and "-fopenmp" in arguments and
            arguments.count("-DNDEBUG") == 1,
            "production compile command is outside the exact Release profile")
    isa_flags = {argument for argument in arguments
                 if argument.startswith("-m")}
    if target == "leopard2_backend_ssse3.dir":
        require(isa_flags == {"-mssse3", "-mno-avx"},
                "SSSE3 production target has the wrong ISA profile")
    elif target == "leopard2_backend_avx2.dir":
        require(isa_flags == {"-mavx2", "-mno-avx512f"} and
                "-falign-functions=64" in arguments,
                "AVX2 production target has the wrong ISA profile")
    elif target == "leopard2_backend_avx512.dir":
        require(isa_flags == {
                    "-mavx2", "-mavx512f", "-mavx512bw", "-mavx512vl",
                    "-mprefer-vector-width=256"} and
                "-falign-functions=64" in arguments,
                "AVX-512 production target has the wrong ISA profile")
    else:
        require(not isa_flags,
                "portable production target contains an isolated ISA flag")


def parse_cmake_cache_text(text_value: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in text_value.splitlines():
        if not line or line.startswith(("#", "//")) or "=" not in line:
            continue
        declaration, value = line.split("=", 1)
        name, separator, unused_type = declaration.partition(":")
        if separator and name:
            values[name] = value
    require(values.get("CMAKE_BUILD_TYPE") == "Release",
            "timing build is not CMake Release")
    require(values.get("CMAKE_GENERATOR") == "Unix Makefiles",
            "timing provenance requires the explicit Unix Makefiles generator")
    return values


def parse_cmake_cache(path: Path) -> tuple[dict[str, str], dict[str, Any]]:
    cache = path.resolve(strict=True)
    values = parse_cmake_cache_text(cache.read_text(encoding="utf-8"))
    return values, file_identity(cache)


def exact_text_identity(path: Path, schema: str) -> dict[str, Any]:
    recipe = path.resolve(strict=True)
    payload = recipe.read_bytes()
    require(0 < len(payload) <= (1 << 20) and b"\0" not in payload,
            "%s is outside the retained recipe byte bound" % recipe)
    result = {
        "schema": schema,
        "file": file_identity(recipe),
        "encoding": "utf-8",
        "text": payload.decode("utf-8", errors="strict"),
    }
    result["sha256"] = sha256_bytes(canonical_bytes(result))
    return result


def validate_compile_entry_binding(
        entry: dict[str, Any], source: Path, build: Path, compiler: Path,
        target: str) -> dict[str, Any]:
    require(isinstance(entry, dict) and
            set(entry) in ({"directory", "command", "file", "output"},
                           {"arguments", "directory", "file", "output"}),
            "production compile entry has unexpected fields")
    expected_source = source.resolve(strict=True)
    require(Path(str(entry["directory"])).resolve(strict=True) == build and
            Path(str(entry["file"])).resolve(strict=True) == expected_source,
            "production compile entry source/build binding is invalid")
    arguments = compile_entry_arguments(entry)
    require(arguments and all("@" not in argument for argument in arguments) and
            Path(arguments[0]).resolve(strict=True) == compiler,
            "production compile entry compiler/response-file binding is invalid")
    validate_compile_profile(arguments, target)
    require(arguments.count("-c") == 1 and
            arguments.index("-c") + 1 < len(arguments) and
            Path(arguments[arguments.index("-c") + 1]).resolve(strict=True) ==
                expected_source and
            [Path(argument).resolve(strict=False) for argument in arguments
             if argument.endswith(".cpp")] == [expected_source],
            "production compile entry does not compile its exact source")
    output_positions = [
        index for index, argument in enumerate(arguments) if argument == "-o"]
    require(len(output_positions) == 1 and
            output_positions[0] + 1 < len(arguments),
            "production compile entry has invalid output semantics")
    output = compile_entry_output(entry)
    token_output = Path(arguments[output_positions[0] + 1])
    if not token_output.is_absolute():
        token_output = build / token_output
    marker = "CMakeFiles/%s/" % target
    require(token_output.resolve(strict=True) == output and
            output.is_relative_to(build) and marker in output.as_posix(),
            "production compile output is not in its exact CMake target")
    result = {
        "schema": "leopard2-production-compile-entry/v1",
        "source": str(expected_source),
        "target": target,
        "entry": entry,
        "arguments": arguments,
        "output": file_identity(output),
    }
    result["sha256"] = sha256_bytes(canonical_bytes(result))
    return result


def object_closure_identity(
        compile_commands: Path, archive_members: dict[str, Any],
        source_root: Path, compiler: Path) -> dict[str, Any]:
    commands = compile_commands.resolve(strict=True)
    commands_content = exact_text_identity(
        commands, "leopard2-compile-commands-content/v1")
    build = commands.parent
    root = source_root.resolve(strict=True)
    members = {}
    for member, filename, target in ARCHIVE_MEMBER_SPECS:
        binding = validate_compile_entry_binding(
            compile_entry_for(commands, filename), root / filename,
            build, compiler, target)
        archive_member = archive_members["member_files"][member]
        output = binding["output"]
        require(output["size"] == archive_member["size"] and
                output["sha256"] == archive_member["sha256"],
                "archive member %s differs from its production object" % member)
        members[member] = {
            "source": filename,
            "target": target,
            "compile": binding,
            "archive_member": archive_member,
        }
    benchmark_binding = validate_compile_entry_binding(
        compile_entry_for(commands, "benchmark.cpp"),
        root / "bench/leopard2/benchmark.cpp", build, compiler,
        "bench_leopard2.dir")
    require(exact_text_identity(
                commands, "leopard2-compile-commands-content/v1") ==
            commands_content,
            "compile_commands.json changed while build closure was captured")
    result = {
        "schema": "leopard2-production-object-closure/v1",
        "compile_commands": file_identity(commands),
        "compile_commands_content": commands_content,
        "source_root": str(root),
        "members": members,
        "benchmark": benchmark_binding,
    }
    result["sha256"] = sha256_bytes(canonical_bytes(result))
    return result


def validate_link_recipes(
        build: Path, executable: Path, archive: Path,
        closure: dict[str, Any], compiler: Path, archiver: Path,
        ranlib: Path) -> dict[str, Any]:
    archive_recipe = exact_text_identity(
        build / "CMakeFiles/leopard.dir/link.txt",
        "leopard2-archive-link-recipe/v1")
    archive_commands = [
        shlex.split(line) for line in archive_recipe["text"].splitlines()
        if line.strip()
    ]
    require(len(archive_commands) == 2,
            "archive link recipe must contain exactly ar and ranlib commands")
    archive_command, ranlib_command = archive_commands
    expected_objects = [
        Path(closure["members"][member]["compile"]["output"]["path"])
            .relative_to(build).as_posix()
        for member in EXPECTED_ARCHIVE_MEMBERS
    ]
    require(len(archive_command) == 3 + len(expected_objects) and
            Path(archive_command[0]).resolve(strict=True) == archiver and
            archive_command[1] in ("qc", "rc", "rcs") and
            archive_command[2] == archive.name and
            archive_command[3:] == expected_objects and
            len(ranlib_command) == 2 and
            Path(ranlib_command[0]).resolve(strict=True) == ranlib and
            ranlib_command[1] == archive.name,
            "archive CMake recipe/tool/object order differs from provenance")

    executable_recipe = exact_text_identity(
        build / "CMakeFiles/bench_leopard2.dir/link.txt",
        "leopard2-executable-link-recipe/v1")
    executable_lines = [
        line for line in executable_recipe["text"].splitlines() if line.strip()]
    require(len(executable_lines) == 1,
            "benchmark link recipe must contain exactly one command")
    link_arguments = shlex.split(executable_lines[0])
    benchmark_object = Path(
        closure["benchmark"]["output"]["path"]).relative_to(build).as_posix()
    require(link_arguments and
            Path(link_arguments[0]).resolve(strict=True) == compiler and
            link_arguments.count(benchmark_object) == 1 and
            link_arguments.count(archive.name) == 1 and
            link_arguments.index(benchmark_object) <
                link_arguments.index(archive.name) and
            link_arguments.count("-o") == 1 and
            link_arguments.index("-o") + 1 < len(link_arguments) and
            link_arguments[link_arguments.index("-o") + 1] == executable.name and
            all("@" not in argument for argument in link_arguments),
            "benchmark CMake link recipe/object/archive order is invalid")
    require([argument for argument in link_arguments
             if argument.endswith(".o")] == [benchmark_object],
            "benchmark link recipe contains an unproven object")
    result = {
        "schema": "leopard2-production-link-recipes/v1",
        "archive": archive_recipe,
        "executable": executable_recipe,
        "archive_objects": expected_objects,
        "benchmark_object": benchmark_object,
        "archive_tool": file_identity(archiver),
        "ranlib_tool": file_identity(ranlib),
        "compiler": file_identity(compiler),
    }
    result["sha256"] = sha256_bytes(canonical_bytes(result))
    return result


def build_provenance_identity(
        executable: Path, archive: Path, compile_commands: Path,
        source_root: Path) -> dict[str, Any]:
    binary = executable.resolve(strict=True)
    library = archive.resolve(strict=True)
    commands = compile_commands.resolve(strict=True)
    build = commands.parent
    root = source_root.resolve(strict=True)
    require(binary.parent == build and library.parent == build,
            "binary/archive do not belong to compile_commands build directory")
    cache_path = build / "CMakeCache.txt"
    cache_content = exact_text_identity(
        cache_path, "leopard2-cmake-cache-content/v1")
    cache, cache_file = parse_cmake_cache(cache_path)
    require(exact_text_identity(
                cache_path, "leopard2-cmake-cache-content/v1") ==
            cache_content,
            "CMakeCache.txt changed while build provenance was captured")
    require(Path(cache.get("CMAKE_HOME_DIRECTORY", "")).resolve(strict=True) ==
                root,
            "CMake source root differs from the attested source root")
    compiler = Path(cache.get("CMAKE_CXX_COMPILER", "")).resolve(strict=True)
    archiver = Path(cache.get("CMAKE_AR", "")).resolve(strict=True)
    ranlib = Path(cache.get("CMAKE_RANLIB", "")).resolve(strict=True)
    archive_members = archive_members_identity(library, archiver)
    closure = object_closure_identity(
        commands, archive_members, root, compiler)
    recipes = validate_link_recipes(
        build, binary, library, closure, compiler, archiver, ranlib)
    result = {
        "schema": "leopard2-production-build-provenance/v1",
        "build": str(build),
        "executable": file_identity(binary),
        "archive": file_identity(library),
        "cache": cache_file,
        "cache_content": cache_content,
        "cache_values": {
            key: cache.get(key) for key in (
                "CMAKE_AR", "CMAKE_BUILD_TYPE", "CMAKE_CXX_COMPILER",
                "CMAKE_CXX_FLAGS", "CMAKE_CXX_FLAGS_RELEASE",
                "CMAKE_EXE_LINKER_FLAGS", "CMAKE_EXE_LINKER_FLAGS_RELEASE",
                "CMAKE_EXPORT_COMPILE_COMMANDS", "CMAKE_GENERATOR",
                "CMAKE_HOME_DIRECTORY", "CMAKE_RANLIB",
                "CMAKE_STATIC_LINKER_FLAGS",
                "CMAKE_STATIC_LINKER_FLAGS_RELEASE",
                "ENABLE_OPENMP", "LEO2_BACKEND_VARIANT",
                "LEO2_BUILD_BENCHMARKS", "LEO2_BUILD_TESTS",
                "LEO2_ENABLE_CUDA", "LEOPARD_ENABLE_GF8",
                "LEOPARD_ENABLE_GF16")
        },
        "archive_members": archive_members,
        "object_closure": closure,
        "link_recipes": recipes,
    }
    result["sha256"] = sha256_bytes(canonical_bytes(result))
    return result


def build_variant_identity(
        executable: Path, compile_commands: Path | None = None) -> dict[str, Any]:
    binary = executable.resolve(strict=True)
    build = binary.parent
    commands = (compile_commands.resolve(strict=True) if compile_commands
                else (build / "compile_commands.json").resolve(strict=True))
    avx2_entry = compile_entry_for(commands, "Leopard2BackendAVX2.cpp")
    avx2_xor_entry = compile_entry_for(
        commands, "Leopard2BackendAVX2Xor.cpp")
    avx512_entry = compile_entry_for(commands, "Leopard2BackendAVX512.cpp")
    core_entry = compile_entry_for(commands, "leopard2.cpp")
    archive = (build / "libleopard.a").resolve(strict=True)
    source_root = compile_entry_source(core_entry).parent
    provenance = build_provenance_identity(
        binary, archive, commands, source_root)
    return {
        "executable": file_identity(binary),
        "archive": file_identity(archive),
        "archive_members": provenance["archive_members"],
        "compile_commands": file_identity(commands),
        "avx2_compile_entry": avx2_entry,
        "avx2_object": file_identity(compile_entry_output(avx2_entry)),
        "avx2_xor_compile_entry": avx2_xor_entry,
        "avx2_xor_object": file_identity(
            compile_entry_output(avx2_xor_entry)),
        "avx512_compile_entry": avx512_entry,
        "avx512_object": file_identity(compile_entry_output(avx512_entry)),
        "core_compile_entry": core_entry,
        "core_object": file_identity(compile_entry_output(core_entry)),
        "build_provenance": provenance,
    }


def git_output(source_root: Path, *arguments: str) -> str:
    completed = subprocess.run(
        ["/usr/bin/git", "-C", str(source_root), *arguments],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
        env=CHILD_ENV)
    require(completed.returncode == 0,
            "git identity command failed: %s" % completed.stderr.decode(
                "utf-8", errors="replace"))
    return completed.stdout.decode("utf-8", errors="strict").strip()


def source_identity(source_root: Path) -> dict[str, Any]:
    root = source_root.resolve(strict=True)
    files = {}
    for relative in SOURCE_FILES:
        path = root / relative
        require(path.is_file(), "source snapshot is missing %s" % relative)
        files[relative] = {
            "size": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    status = git_output(root, "status", "--short", "--untracked-files=all")
    result = {
        "root": str(root),
        "head": git_output(root, "rev-parse", "HEAD"),
        "head_tree": git_output(root, "rev-parse", "HEAD^{tree}"),
        "status_short": status.splitlines() if status else [],
        "files": files,
    }
    result["snapshot_sha256"] = sha256_bytes(canonical_bytes(result))
    return result


def atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp.%d" % os.getpid())
    payload = canonical_bytes(value) + b"\n"
    descriptor = os.open(
        temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC,
        0o600)
    try:
        offset = 0
        while offset < len(payload):
            written = os.write(descriptor, payload[offset:])
            require(written > 0, "short atomic JSON write")
            offset += written
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    os.replace(str(temporary), str(path))
    directory = os.open(
        path.parent, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)


def stable_seed(identifier: str, schema: str = MATRIX_SCHEMA) -> int:
    digest = hashlib.sha256((schema + ":" + identifier).encode()).digest()
    return int.from_bytes(digest[:8], "big") & ((1 << 63) - 1)


def iteration_policy(byte_count: int) -> tuple[int, int]:
    if byte_count <= 4096:
        return 15, 4
    return 9, 3


def ceil_power_of_two(value: int) -> int:
    require(type(value) is int and value > 0,
            "power-of-two input must be positive")
    return 1 << (value - 1).bit_length()


def expected_direct_executor(mode: str, cell: dict[str, Any]) -> str:
    mode_compile_definition(mode)
    if (cell.get("K") == 65 and cell.get("R") == 65 and
            2 <= cell.get("loss", 0) <= 8):
        return "source_major" if cell.get("bytes", 0) >= 2048 \
            else "output_major"
    if cell["loss"] <= 4:
        return "output_major"
    if mode == "transform":
        return "none"
    if mode in ("source", "production",
                XMM_TAIL_CONTROL_MODE, XMM_TAIL_CANDIDATE_MODE):
        return "source_major"
    return mode + "_major"


def make_matrix() -> dict[str, Any]:
    cells = []
    index = 0
    for original_count in (5, 8, 9, 12, 16):
        shapes = [(5, 4), (8, 4), (5, 5), (8, 5)]
        if original_count >= 6:
            shapes.extend(((6, 6), (8, 6)))
        if original_count >= 7:
            shapes.extend(((7, 7), (8, 7)))
        if original_count >= 8:
            shapes.append((8, 8))
        for recovery_count, loss_count in shapes:
            for byte_count in (64, 2048, 4096, 65536):
                identifier = "k%d-r%d-b%d-l%d" % (
                    original_count, recovery_count, byte_count, loss_count)
                iterations, warmup = iteration_policy(byte_count)
                cells.append({
                    "index": index,
                    "id": identifier,
                    "K": original_count,
                    "R": recovery_count,
                    "bytes": byte_count,
                    "loss": loss_count,
                    "batch": 1,
                    "reuse": 8,
                    "iterations": iterations,
                    "warmup": warmup,
                    "seed": stable_seed(identifier),
                    "estimated_peak_bytes": 6 *
                        (original_count + recovery_count) * byte_count +
                        (64 << 20),
                    "role": "loss4_neighbor" if loss_count == 4
                        else "loss5_to_8_target",
                    "exact_main_required": (
                        recovery_count <= original_count and loss_count >= 5
                        and byte_count in (2048, 4096, 65536)),
                })
                index += 1
    result = {
        "schema": MATRIX_SCHEMA,
        "cell_count": len(cells),
        "cells": cells,
        "directional_scope": {
            "baseline_candidate_pairs": [
                ["transform", "output"],
                ["transform", "source"]
            ],
            "direct_head_to_head": (
                "Run output versus source only for overlapping winners; do "
                "not infer it by dividing separate noisy ratios."),
        },
        "same_source_promotion": {
            "target": "whole decode_execution",
            "candidate_ci95_low_at_least": 1.05,
            "neighbor_ci95_high_at_least": 0.98,
            "orientation": "baseline_time_over_candidate_time",
        },
        "exact_main_note": (
            "Run the cells marked exact_main_required with "
            "experiments/leopard2/main_compare/run_abba.py against exact "
            "Leopard main commit 6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198. "
            "Cells with R>K are outside Leopard1's R<=K API contract."),
    }
    result["matrix_sha256"] = sha256_bytes(canonical_bytes(result))
    return result


def make_large_matrix() -> dict[str, Any]:
    cells = []
    index = 0
    for original_count in range(5, 17):
        for recovery_count in range(5, 9):
            for loss_count in range(4, min(original_count, recovery_count) + 1):
                for byte_count in (
                        64, 65, 256, 1024, 2048, 2049, 4096, 65536):
                    identifier = "full-k%d-r%d-b%d-l%d" % (
                        original_count, recovery_count, byte_count, loss_count)
                    iterations, warmup = iteration_policy(byte_count)
                    cells.append({
                        "index": index,
                        "id": identifier,
                        "K": original_count,
                        "R": recovery_count,
                        "bytes": byte_count,
                        "loss": loss_count,
                        "batch": 1,
                        "reuse": 8,
                        "iterations": iterations,
                        "warmup": warmup,
                        "seed": stable_seed(identifier, LARGE_MATRIX_SCHEMA),
                        "estimated_peak_bytes": 6 *
                            (original_count + recovery_count) * byte_count +
                            (64 << 20),
                        "role": "loss4_neighbor" if loss_count == 4
                            else "loss5_to_8_target",
                        "requested_tile_screen_loss": False,
                        "exact_main_required": (
                            recovery_count <= original_count and
                            loss_count >= 5 and
                            byte_count in (2048, 4096, 65536)),
                    })
                    index += 1
    result = {
        "schema": LARGE_MATRIX_SCHEMA,
        "cell_count": len(cells),
        "cells": cells,
        "same_source_promotion": {
            "target": "whole decode_execution",
            "candidate_ci95_low_at_least": 1.05,
            "neighbor_ci95_high_at_least": 0.98,
            "orientation": "baseline_time_over_candidate_time",
        },
        "coverage": (
            "Every K=5..16, R=5..8, L=4..min(K,R), and requested "
            "64-byte through 64-KiB execution/tail boundary."),
    }
    result["matrix_sha256"] = sha256_bytes(canonical_bytes(result))
    return result


def make_tiny_matrix() -> dict[str, Any]:
    """Cover every public byte-tail boundary below the 64-byte screen.

    Source-major direct repair transposes fixed multipliers into a stack-local
    table on each execution.  That setup is amortized at ordinary shard sizes,
    but cannot be assumed free for tiny shards.  Keep this matrix separate so
    retained 64-byte-and-larger evidence remains stable.
    """
    cells = []
    index = 0
    shapes = (
        (5, 5, 4), (5, 5, 5),
        (8, 8, 4), (8, 8, 5), (8, 8, 8),
        (16, 8, 4), (16, 8, 5), (16, 8, 8),
    )
    for original_count, recovery_count, loss_count in shapes:
        for byte_count in (1, 2, 3, 7, 8, 15, 16, 17, 31, 32, 33, 63):
            identifier = "tiny-k%d-r%d-b%d-l%d" % (
                original_count, recovery_count, byte_count, loss_count)
            iterations, warmup = iteration_policy(byte_count)
            cells.append({
                "index": index,
                "id": identifier,
                "K": original_count,
                "R": recovery_count,
                "bytes": byte_count,
                "loss": loss_count,
                "batch": 1,
                "reuse": 8,
                "iterations": iterations,
                "warmup": warmup,
                "seed": stable_seed(identifier, TINY_MATRIX_SCHEMA),
                "estimated_peak_bytes": 6 *
                    (original_count + recovery_count) * byte_count +
                    (64 << 20),
                "role": "loss4_neighbor" if loss_count == 4
                    else "loss5_to_8_target",
                "exact_main_required": False,
            })
            index += 1
    result = {
        "schema": TINY_MATRIX_SCHEMA,
        "cell_count": len(cells),
        "cells": cells,
        "same_source_promotion": {
            "target": "whole decode_execution",
            "candidate_ci95_low_at_least": 1.05,
            "neighbor_ci95_high_at_least": 0.98,
            "orientation": "baseline_time_over_candidate_time",
        },
        "coverage": (
            "All required positive byte-tail boundaries below 64 bytes for "
            "representative K=5,8,16 and loss=4,5,8 shapes."),
    }
    result["matrix_sha256"] = sha256_bytes(canonical_bytes(result))
    return result


def xmm_tail_segments(byte_count: int) -> list[int]:
    """Return the VEX-XMM tail widths reached after complete YMM chunks."""
    remainder = byte_count % 32
    result = []
    if remainder >= 16:
        result.append(16)
        remainder -= 16
    if remainder >= 8:
        result.append(8)
    return result


def direct_kernel_shape(loss_count: int) -> str:
    shapes = {
        2: "pair",
        4: "group4",
        5: "group4_single",
        6: "group4_pair",
        7: "group4_pair_single",
        8: "group4_group4",
    }
    require(loss_count in shapes,
            "XMM-tail matrix has an unsupported direct kernel shape")
    return shapes[loss_count]


def make_xmm_tail_matrix() -> dict[str, Any]:
    """Attribute both direct kernels at every 8-byte tail transition.

    The compact K=8 block exercises Group4, Pair, Pair+single, and composed
    Group4 paths.  The five K=65 boundaries bind the experiment to the
    existing production-size source-major regime and include exact no-tail
    neighbors.  This is separately versioned so the historical 96-cell tiny
    matrix and all retained evidence remain byte-stable.
    """
    cells = []
    index = 0
    residue_boundaries = tuple(sorted({
        value + delta
        for value in (8, 16, 24, 32, 40, 48, 56, 64)
        for delta in (-1, 0, 1)
    }))
    shapes = [
        (8, 8, loss_count, byte_count)
        for loss_count in (4, 5, 6, 7, 8)
        for byte_count in residue_boundaries
    ]
    shapes.extend(
        (65, 65, loss_count, byte_count)
        for loss_count in (2, 4, 8)
        for byte_count in (2048, 2063, 2064, 2079, 2080)
    )
    for original_count, recovery_count, loss_count, byte_count in shapes:
        kernel_shape = direct_kernel_shape(loss_count)
        executor = expected_direct_executor(XMM_TAIL_CONTROL_MODE, {
            "K": original_count, "R": recovery_count,
            "bytes": byte_count, "loss": loss_count,
        })
        segments = xmm_tail_segments(byte_count) \
            if executor == "source_major" else []
        active = bool(segments)
        identifier = "xmm-tail-k%d-r%d-b%d-l%d" % (
            original_count, recovery_count, byte_count, loss_count)
        iterations, warmup = iteration_policy(byte_count)
        cells.append({
            "index": index,
            "id": identifier,
            "K": original_count,
            "R": recovery_count,
            "bytes": byte_count,
            "loss": loss_count,
            "batch": 1,
            "reuse": 8,
            "iterations": iterations,
            "warmup": warmup,
            "seed": stable_seed(identifier, XMM_TAIL_MATRIX_SCHEMA),
            "estimated_peak_bytes": 6 *
                (original_count + recovery_count) * byte_count +
                (64 << 20),
            "role": ("xmm_tail_target_" if active else
                     "no_xmm_regression_neighbor_") + kernel_shape,
            "kernel_shape": kernel_shape,
            "xmm_tail_executes": active,
            "xmm_tail_segments": segments,
            "expected_executor": executor,
            "expected_decode_path": "direct",
            "exact_main_required": False,
        })
        index += 1
    result = {
        "schema": XMM_TAIL_MATRIX_SCHEMA,
        "cell_count": len(cells),
        "cells": cells,
        "same_source_promotion": {
            "target": "whole decode_execution",
            "candidate_ci95_low_at_least": 1.05,
            "neighbor_ci95_high_at_least": 0.98,
            "orientation": "baseline_time_over_candidate_time",
        },
        "coverage": (
            "K=8/R=8 losses 4..8 at n-1/n/n+1 around every 8-byte "
            "boundary from 8 through 64, plus K=65/R=65 losses 2,4,8 "
            "at 2048/2063/2064/2079/2080 bytes.  Complete selection is "
            "required so Group4, Pair, composed, active-tail, and no-tail "
            "neighbor behavior are all retained together."),
    }
    result["matrix_sha256"] = sha256_bytes(canonical_bytes(result))
    return result


def matrix_for_preset(preset: str) -> dict[str, Any]:
    if preset == "core":
        return make_matrix()
    if preset == "large":
        return make_large_matrix()
    if preset == "tiny":
        return make_tiny_matrix()
    if preset == "xmm-tail":
        return make_xmm_tail_matrix()
    raise EvidenceError("unknown matrix preset")


def validate_comparison_matrix_selection(
        preset: str, modes: dict[str, str], matrix: dict[str, Any],
        selected: list[dict[str, Any]]) -> None:
    """Keep XMM attribution inseparable from its complete control matrix."""
    if not is_xmm_tail_comparison(modes):
        require(preset != "xmm-tail",
                "the XMM-tail preset requires the XMM-tail comparison")
        return
    require(preset == "xmm-tail" and
            matrix.get("schema") == XMM_TAIL_MATRIX_SCHEMA,
            "XMM-tail comparison requires the exact XMM-tail preset")
    require([cell.get("id") for cell in selected] ==
                [cell.get("id") for cell in matrix.get("cells", [])],
            "XMM-tail comparison requires every frozen matrix cell")
    active = set()
    inactive = set()
    executors = set()
    for cell in selected:
        shape = direct_kernel_shape(cell["loss"])
        executor = expected_direct_executor(XMM_TAIL_CONTROL_MODE, cell)
        segments = xmm_tail_segments(cell["bytes"]) \
            if executor == "source_major" else []
        require(cell.get("kernel_shape") == shape and
                cell.get("xmm_tail_segments") == segments and
                cell.get("xmm_tail_executes") == bool(segments) and
                cell.get("expected_executor") == executor and
                cell.get("expected_decode_path") == "direct" and
                cell.get("role") == (("xmm_tail_target_" if segments else
                    "no_xmm_regression_neighbor_") + shape),
                "XMM-tail cell role or path contract is invalid")
        (active if segments else inactive).add(shape)
        executors.add(executor)
    required_shapes = {
        "pair", "group4", "group4_single", "group4_pair",
        "group4_pair_single", "group4_group4",
    }
    require(active == required_shapes and
            {"pair", "group4", "group4_group4"}.issubset(inactive) and
            executors == {"output_major", "source_major"},
            "XMM-tail selection does not attribute every required path")


def parse_cpu_list(text: str) -> set[int]:
    result: set[int] = set()
    for component in text.strip().split(","):
        if not component:
            continue
        bounds = component.split("-", 1)
        first = int(bounds[0])
        last = int(bounds[1]) if len(bounds) == 2 else first
        require(0 <= first <= last, "invalid CPU list")
        result.update(range(first, last + 1))
    return result


def cgroup_effective_cpus() -> set[int]:
    """Return the inherited cgroup-v2 CPU allowance, not process affinity."""
    unified = None
    for line in Path("/proc/self/cgroup").read_text().splitlines():
        fields = line.split(":", 2)
        if len(fields) == 3 and fields[0] == "0" and fields[1] == "":
            unified = fields[2]
            break
    require(unified is not None, "unified cgroup-v2 membership is unavailable")
    root = Path("/sys/fs/cgroup").resolve(strict=True)
    current = (root / unified.lstrip("/")).resolve(strict=True)
    require(current == root or root in current.parents,
            "process cgroup escapes the cgroup-v2 root")
    while True:
        effective = current / "cpuset.cpus.effective"
        if effective.is_file():
            text = effective.read_text().strip()
            if text:
                return parse_cpu_list(text)
        if current == root:
            break
        current = current.parent
    raise EvidenceError("cgroup-v2 effective CPU set is unavailable")


def validate_smt_pair_identity(
        cpu: int, sibling: int, cgroup_allowed: set[int],
        launch_affinity: set[int], reported_siblings: set[int]
) -> dict[str, Any]:
    require(cpu != sibling, "benchmark CPU and reserved sibling must differ")
    require(cpu in cgroup_allowed and sibling in cgroup_allowed,
            "benchmark CPU pair is outside the cgroup-effective CPU set")
    require(launch_affinity - {cpu, sibling},
            "no housekeeping CPU remains outside the benchmark pair")
    require(reported_siblings == {cpu, sibling},
            "requested CPUs are not the complete two-thread SMT sibling pair")
    return {
        "cgroup_effective_cpus": sorted(cgroup_allowed),
        "launch_affinity": sorted(launch_affinity),
        "housekeeping_affinity": sorted(launch_affinity - {cpu, sibling}),
    }


def require_smt_pair(cpu: int, sibling: int) -> dict[str, Any]:
    cgroup_allowed = cgroup_effective_cpus()
    launch_affinity = set(os.sched_getaffinity(0))
    topology = Path("/sys/devices/system/cpu/cpu%d/topology/thread_siblings_list" % cpu)
    require(topology.is_file(), "CPU topology file is unavailable")
    siblings = parse_cpu_list(topology.read_text())
    return validate_smt_pair_identity(
        cpu, sibling, cgroup_allowed, launch_affinity, siblings)


def load_pair_lease(source_root: Path):
    path = (source_root / "tools/leopard2_jerasure_compare.py").resolve(
        strict=True)
    specification = importlib.util.spec_from_file_location(
        "leopard2_direct_repair_pair_lease", path)
    require(specification is not None and specification.loader is not None,
            "cannot import canonical CPU-pair lease")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    try:
        specification.loader.exec_module(module)
    except Exception:
        sys.modules.pop(specification.name, None)
        raise
    require(hasattr(module, "PairLease"),
            "canonical lease module does not export PairLease")
    return module.PairLease, file_identity(path)


def same_user_pair_affinity(cpu: int, sibling: int) -> list[dict[str, Any]]:
    pair = {cpu, sibling}
    uid = os.getuid()
    found = []
    for process in Path("/proc").iterdir():
        if not process.name.isdigit():
            continue
        try:
            if process.stat().st_uid != uid:
                continue
            command = (process / "comm").read_text().strip()
            tasks = process / "task"
            for task in tasks.iterdir():
                if not task.name.isdigit():
                    continue
                affinity = set(os.sched_getaffinity(int(task.name)))
                if affinity & pair:
                    stat_text = (task / "stat").read_text()
                    stat_end = stat_text.rfind(")")
                    require(stat_end >= 0, "malformed /proc task identity")
                    stat_fields = stat_text[stat_end + 2:].split()
                    require(len(stat_fields) > 19,
                            "short /proc task identity")
                    found.append({
                        "pid": int(process.name), "tid": int(task.name),
                        "command": command, "affinity": sorted(affinity),
                        "start_time_ticks": int(stat_fields[19]),
                        "task_inode": task.stat(
                            follow_symlinks=False).st_ino,
                    })
        except (FileNotFoundError, ProcessLookupError):
            continue
    return sorted(found, key=lambda value: (value["pid"], value["tid"]))


def current_task_identity(pid: int, tid: int) -> tuple[int, int]:
    task = Path("/proc") / str(pid) / "task" / str(tid)
    inode = task.stat(follow_symlinks=False).st_ino
    text = (task / "stat").read_text()
    end = text.rfind(")")
    require(end >= 0, "malformed current /proc task identity")
    fields = text[end + 2:].split()
    require(len(fields) > 19, "short current /proc task identity")
    return int(fields[19]), inode


def current_process_state(pid: int) -> dict[str, Any]:
    text = (Path("/proc") / str(pid) / "stat").read_text()
    end = text.rfind(")")
    require(end >= 0, "malformed current process identity")
    fields = text[end + 2:].split()
    require(len(fields) > 36, "short current process identity")
    return {
        "state": fields[0], "pgrp": int(fields[2]),
        "session": int(fields[3]), "processor": int(fields[36]),
    }


def require_current_task(record: dict[str, Any]) -> None:
    try:
        start_time, inode = current_task_identity(
            record["pid"], record["tid"])
    except (FileNotFoundError, ProcessLookupError) as error:
        raise EvidenceError(
            "recorded affinity task vanished before mutation: %s" % record) \
            from error
    require(start_time == record["start_time_ticks"],
            "recorded affinity task identity was reused or replaced")
    require(inode == record["task_inode"],
            "recorded affinity task inode was reused or replaced")


def runner_affinity_record() -> dict[str, Any]:
    pid = os.getpid()
    start_time, inode = current_task_identity(pid, pid)
    return {
        "pid": pid, "tid": pid,
        "affinity": sorted(os.sched_getaffinity(pid)),
        "start_time_ticks": start_time, "task_inode": inode,
    }


def exclude_same_user_from_pair(cpu: int, sibling: int) -> dict[str, Any]:
    pair = {cpu, sibling}
    started_ns = time.time_ns()
    runner_before = runner_affinity_record()
    fallback_housekeeping = cgroup_effective_cpus() - pair
    require(fallback_housekeeping,
            "cgroup has no housekeeping CPU outside reserved pair")
    changed = []
    quiescence_passes = 0
    try:
        for quiescence_passes in range(1, 9):
            pending = same_user_pair_affinity(cpu, sibling)
            if not pending:
                break
            for record in pending:
                require_current_task(record)
                before = set(record["affinity"])
                current = set(os.sched_getaffinity(record["tid"]))
                require(current == before,
                        "recorded affinity changed before exclusion mutation")
                after = before - pair
                if not after:
                    after = fallback_housekeeping
                try:
                    os.sched_setaffinity(record["tid"], after)
                except OSError as error:
                    raise EvidenceError(
                        "failed to exclude recorded task from CPU pair: %s" %
                        record) from error
                require_current_task(record)
                require(set(os.sched_getaffinity(record["tid"])) == after,
                        "affinity exclusion mutation did not take effect")
                changed.append({**record, "after": sorted(after)})
        remaining = same_user_pair_affinity(cpu, sibling)
        require(not remaining,
                "same-user affinity did not quiesce after eight passes: %s" %
                remaining)
    except Exception as error:
        rollback_error = None
        try:
            restore_same_user_affinity({"changed": changed})
        except Exception as rollback:
            rollback_error = rollback
        if rollback_error is not None:
            raise EvidenceError(
                "affinity exclusion failed and rollback also failed: %s; %s" %
                (error, rollback_error)) from error
        raise
    return {
        "schema": "leopard2-user-affinity-exclusion/v1",
        "uid": os.getuid(), "pair": sorted(pair), "changed": changed,
        "after_snapshot": remaining, "started_ns": started_ns,
        "completed_ns": time.time_ns(),
        "quiescence_passes": quiescence_passes,
        "runner_before": runner_before,
        "runner_after": runner_affinity_record(),
    }


def restore_same_user_affinity(exclusion: dict[str, Any]) -> dict[str, Any]:
    restored = []
    failures = []
    for record in reversed(exclusion.get("changed", [])):
        tid = record["tid"]
        try:
            require_current_task(record)
            current = sorted(os.sched_getaffinity(tid))
            if current != record["after"]:
                failures.append({
                    "reason": "changed_elsewhere",
                    "pid": record["pid"], "tid": tid,
                    "expected": record["after"], "actual": current,
                })
                continue
            os.sched_setaffinity(tid, set(record["affinity"]))
            require_current_task(record)
            actual = sorted(os.sched_getaffinity(tid))
            if actual != record["affinity"]:
                failures.append({
                    "reason": "restore_not_applied",
                    "pid": record["pid"], "tid": tid,
                    "expected": record["affinity"], "actual": actual,
                })
                continue
            restored.append({"pid": record["pid"], "tid": tid,
                             "affinity": record["affinity"],
                             "start_time_ticks": record["start_time_ticks"],
                             "task_inode": record["task_inode"]})
        except (EvidenceError, OSError, ValueError) as error:
            failures.append({
                "reason": "restore_exception", "pid": record["pid"],
                "tid": tid, "error": str(error),
            })
    require(not failures,
            "excluded thread affinity restoration failed: %s" % failures)
    return {
        "schema": "leopard2-user-affinity-restoration/v1",
        "restored": restored, "failures": failures,
        "runner_restored": runner_affinity_record(),
    }


def reservation_identity(path: Path, cpu: int, sibling: int) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    value = json.loads(resolved.read_text())
    require(isinstance(value, dict) and
            value.get("schema") == "leopard2-cpu-reservation/v1" and
            value.get("status") == "held" and
            value.get("benchmark_cpu") == cpu and
            value.get("reserved_sibling") == sibling and
            isinstance(value.get("owner"), str) and value["owner"] and
            isinstance(value.get("nonce"), str) and value["nonce"],
            "CPU reservation does not match the requested held pair")
    return {"file": file_identity(resolved), "value": value}


def cpu_ticks(cpu: int) -> tuple[int, ...]:
    prefix = "cpu%d " % cpu
    with Path("/proc/stat").open("r", encoding="ascii") as stream:
        for line in stream:
            if line.startswith(prefix):
                fields = tuple(int(value) for value in line.split()[1:])
                require(len(fields) >= 5, "short /proc/stat CPU record")
                return fields
    raise EvidenceError("reserved sibling is absent from /proc/stat")


def cpu_scheduler_runtime_ns(cpu: int) -> int:
    lines = Path("/proc/schedstat").read_text(encoding="ascii").splitlines()
    require(lines and lines[0].startswith("version "),
            "scheduler runtime version is unavailable")
    version = int(lines[0].split()[1])
    require(version == 15,
            "scheduler runtime accounting requires exactly schedstat v15")
    prefix = "cpu%d " % cpu
    for line in lines:
        if line.startswith(prefix):
            fields = line.split()
            require(len(fields) == 10,
                    "unexpected schedstat v15 CPU record")
            return int(fields[7])
    raise EvidenceError("benchmark CPU is absent from /proc/schedstat")


def child_cpu_time_ns(before: Any, after: Any) -> int:
    seconds = (after.ru_utime - before.ru_utime) + \
        (after.ru_stime - before.ru_stime)
    require(math.isfinite(seconds) and seconds >= 0,
            "child CPU accounting moved backwards")
    return int(round(seconds * 1_000_000_000.0))


def rusage_cpu_time_ns(value: Any) -> int:
    seconds = value.ru_utime + value.ru_stime
    require(math.isfinite(seconds) and seconds >= 0,
            "child CPU accounting is invalid")
    return int(round(seconds * 1_000_000_000.0))


def task_scheduler_runtime_ns(pid: int) -> int:
    fields = (Path("/proc") / str(pid) / "schedstat").read_text().split()
    require(len(fields) == 3,
            "unexpected per-task schedstat record")
    return require_int(int(fields[0]), "task scheduler runtime")


def target_runtime_evidence(
        scheduler_before_ns: int, scheduler_after_ns: int,
        child_runtime_ns: int,
        tolerance_ns: int = TARGET_RUNTIME_TOLERANCE_NS) -> dict[str, Any]:
    require(scheduler_after_ns >= scheduler_before_ns and
            child_runtime_ns >= 0,
            "target CPU runtime accounting moved backwards")
    scheduler_delta = scheduler_after_ns - scheduler_before_ns
    signed_difference = scheduler_delta - child_runtime_ns
    return {
        "scheduler_before_ns": scheduler_before_ns,
        "scheduler_after_ns": scheduler_after_ns,
        "scheduler_delta_ns": scheduler_delta,
        "child_cpu_time_ns": child_runtime_ns,
        "signed_difference_ns": signed_difference,
        "unexplained_runtime_ns": max(0, signed_difference),
        "tolerance_ns": tolerance_ns,
        "accepted": abs(signed_difference) <= tolerance_ns,
    }


def nonidle_delta(before: tuple[int, ...], after: tuple[int, ...]) -> int:
    require(len(before) == len(after), "/proc/stat field count changed")
    total = sum(after) - sum(before)
    idle = (after[3] - before[3]) + (after[4] - before[4])
    return total - idle


def target_interrupt_evidence(
        before: tuple[int, ...], after: tuple[int, ...]) -> dict[str, Any]:
    require(len(before) == len(after), "/proc/stat field count changed")
    deltas = [later - earlier for earlier, later in zip(before, after)]
    require(all(value >= 0 for value in deltas),
            "/proc/stat target CPU fields moved backwards")
    names = (
        "user", "nice", "system", "idle", "iowait", "irq", "softirq",
        "steal", "guest", "guest_nice",
    )
    fields = {
        names[index] if index < len(names) else "field_%d" % index: value
        for index, value in enumerate(deltas)
    }
    rejected = {
        name: fields.get(name, 0)
        for name in ("irq", "softirq", "steal", "guest", "guest_nice")
    }
    return {
        "before": before, "after": after, "deltas": fields,
        "rejected_fields": rejected,
        "accepted": all(value == 0 for value in rejected.values()),
    }


def acquire_lock(path: Path, timeout: float):
    require(path == DEFAULT_LOCK and path.is_absolute(),
            "authoritative campaign requires the canonical GF8 lock path")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(
        path, os.O_RDWR | os.O_CREAT | os.O_CLOEXEC | os.O_NOFOLLOW, 0o600)
    stream = os.fdopen(descriptor, "r+")
    deadline = time.monotonic() + timeout
    while True:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            os.fchmod(stream.fileno(), 0o600)
            return stream
        except BlockingIOError:
            if time.monotonic() >= deadline:
                stream.close()
                raise EvidenceError("timed out waiting for canonical GF8 lock")
            time.sleep(1.0)


def global_lock_identity(stream: Any, path: Path) -> dict[str, Any]:
    require(path == DEFAULT_LOCK and path.is_absolute(),
            "noncanonical GF8 lock path")
    opened = os.fstat(stream.fileno())
    linked = path.lstat()
    require(stat.S_ISREG(opened.st_mode) and opened.st_uid == os.getuid() and
            opened.st_nlink == 1 and
            stat.S_IMODE(opened.st_mode) == 0o600 and
            opened.st_dev == linked.st_dev and
            opened.st_ino == linked.st_ino,
            "canonical GF8 lock path was replaced or is not regular")
    return {
        "schema": "leopard2-global-benchmark-lock/v1",
        "path": str(path), "device": opened.st_dev,
        "inode": opened.st_ino, "uid": opened.st_uid,
        "mode": stat.S_IMODE(opened.st_mode), "nlink": opened.st_nlink,
        "mechanism": "fcntl-flock-exclusive",
    }


def validate_result(result: Any, cell: dict[str, Any], mode: str) -> dict[str, Any]:
    require(isinstance(result, dict), "benchmark output is not an object")
    require(result.get("schema") == "leopard2-benchmark-v3",
            "benchmark did not emit the path-attested v3 schema")
    parameters = result.get("parameters")
    resolved = result.get("resolved")
    correctness = result.get("correctness")
    digests = result.get("workload_digests")
    metrics = result.get("metrics")
    for value, name in ((parameters, "parameters"), (resolved, "resolved"),
                        (correctness, "correctness"), (digests, "digests"),
                        (metrics, "metrics")):
        require(isinstance(value, dict), "benchmark %s is missing" % name)
    expected_parameters = {
        "K": cell["K"], "R": cell["R"],
        "shard_bytes": cell["bytes"], "loss_count": cell["loss"],
        "batch": cell["batch"], "reuse": cell["reuse"],
        "iterations": cell["iterations"], "warmup": cell["warmup"],
        "thread_count": 1, "seed": cell["seed"],
        "requested_profile": "legacy_high_v1",
        "requested_field": "gf8", "requested_backend": "avx2",
        "force_generic_decode": False,
        "force_specialized_decode": False,
        "force_tiled_decode": False,
        "force_materialized_decode": False,
        "skip_legacy": True, "retain_samples": True,
        "report_decode_path": True,
    }
    for key, expected in expected_parameters.items():
        require(parameters.get(key) == expected,
                "benchmark parameter %s does not match request" % key)
    require(resolved.get("profile") == "legacy_high_v1",
            "benchmark resolved the wrong profile")
    require(resolved.get("field") == "gf8", "benchmark resolved the wrong field")
    require(resolved.get("backend") == "avx2",
            "benchmark did not resolve the AVX2 backend")
    require(resolved.get("thread_count") == 1,
            "benchmark did not resolve one execution thread")
    padded_side = ceil_power_of_two(cell["R"])
    parent_count = ceil_power_of_two(cell["K"] + padded_side)
    require(resolved.get("parent_count") == parent_count and
            resolved.get("padded_side") == padded_side,
            "benchmark resolved the wrong direct-repair parent geometry")
    expected_executor = expected_direct_executor(mode, cell)
    selected_path = resolved.get("selected_decode_path")
    selected_rule = resolved.get("selected_decode_rule")
    if "expected_executor" in cell or "expected_decode_path" in cell:
        require(cell.get("expected_executor") == expected_executor and
                cell.get("expected_decode_path") == "direct",
                "frozen cell has the wrong selector/path assertion")
    if expected_executor == "none":
        require(selected_path in ("generic", "materialized", "tiled") and
                selected_rule not in ("no_op", "direct", "unsupported_profile"),
                "benchmark did not select transform repair")
    else:
        require(selected_path == "direct" and selected_rule == "direct",
                "benchmark did not select direct repair")
    require(correctness.get("leopard2_round_trip") is True,
            "benchmark round trip failed")
    require(digests.get("algorithm") == "fnv1a64",
            "benchmark emitted the wrong digest algorithm")
    for key in ("original_data", "transmitted_parity", "recovered_originals"):
        value = digests.get(key)
        require(isinstance(value, str) and
                re.fullmatch(r"[0-9a-f]{16}", value) is not None,
                "benchmark emitted an invalid %s digest" % key)
    missing = parameters.get("missing_original_indices")
    require(isinstance(missing, list) and len(missing) == cell["loss"] and
            all(type(index) is int and 0 <= index < cell["K"]
                for index in missing) and
            missing == sorted(set(missing)),
            "benchmark emitted an invalid missing-index list")
    execution = metrics.get("decode_execution")
    setup = metrics.get("decode_plan_setup")
    require(isinstance(execution, dict), "decode execution metric is missing")
    require(isinstance(setup, dict), "decode plan-setup metric is missing")
    median = execution.get("median_us_per_batch_call")
    setup_median = setup.get("median_us")
    samples = execution.get("samples_us_per_batch_call")
    require(type(median) in (int, float) and
            math.isfinite(median) and median > 0,
            "decode execution median is invalid")
    require(isinstance(samples, list) and len(samples) == cell["iterations"],
            "decode execution samples are missing")
    require(all(type(sample) in (int, float) and
                math.isfinite(sample) and sample > 0 for sample in samples),
            "decode execution samples contain an invalid value")
    require(type(setup_median) in (int, float) and
            math.isfinite(setup_median) and setup_median >= 0,
            "decode plan-setup median is invalid")
    return {
        "median_us": float(median),
        "plan_setup_us": float(setup_median),
        "digests": digests,
        "missing_original_indices": missing,
        "decode_path": selected_path,
        "decode_rule": selected_rule,
        # This identity is bound by the exact retained compiler definition and
        # core-object/source closure.  The C++ regression suite separately
        # queries the shared private selector.  Do not extend benchmark-v3's
        # exact JSON map merely for an experiment-only attribution field.
        "build_bound_executor": expected_executor,
    }


def ratio_summary(invocations: list[dict[str, Any]], rounds: int,
                  metric) -> dict[str, Any]:
    contrasts = []
    round_ratios = []
    for round_index in range(rounds):
        group = invocations[round_index * len(ORDER):(round_index + 1) * len(ORDER)]
        baseline = [metric(item["normalized"]) for item in group
                    if item["implementation"] == "baseline"]
        candidate = [metric(item["normalized"]) for item in group
                     if item["implementation"] == "candidate"]
        require(len(baseline) == 2 and len(candidate) == 2 and
                all(value > 0 for value in baseline + candidate),
                "ABBA round has an invalid metric")
        contrast = (sum(math.log(value) for value in baseline) -
                    sum(math.log(value) for value in candidate)) / 2.0
        contrasts.append(contrast)
        round_ratios.append(math.exp(contrast))
    mean_log = statistics.mean(contrasts)
    if rounds == 1:
        interval = [None, None]
    else:
        standard_error = statistics.stdev(contrasts) / math.sqrt(rounds)
        critical = T_CRITICAL_95[rounds - 1]
        interval = [math.exp(mean_log - critical * standard_error),
                    math.exp(mean_log + critical * standard_error)]
    return {
        "round_ratios": round_ratios,
        "geometric_mean_ratio": math.exp(mean_log),
        "ci95": interval,
        "directional_only": rounds == 1,
    }


def benchmark_arguments(cell: dict[str, Any]) -> list[str]:
    return [
        "--k", str(cell["K"]), "--r", str(cell["R"]),
        "--loss", str(cell["loss"]), "--bytes", str(cell["bytes"]),
        "--batch", str(cell["batch"]), "--profile", "high",
        "--field", "gf8", "--backend", "avx2",
        "--reuse", str(cell["reuse"]),
        "--iterations", str(cell["iterations"]),
        "--warmup", str(cell["warmup"]), "--threads", "1",
        "--seed", str(cell["seed"]), "--skip-legacy",
        "--retain-samples", "--report-decode-path", "--json", "-",
    ]


def benchmark_command(binary: Path, cell: dict[str, Any], cpu: int) -> list[str]:
    return [
        str(binary),
        *benchmark_arguments(cell),
    ]


def wait4_until(pid: int, flags: int, deadline: float):
    while True:
        result = os.wait4(pid, flags | os.WNOHANG)
        if result[0] == pid:
            return result
        require(result[0] == 0, "wait4 returned an unexpected child")
        if time.monotonic() >= deadline:
            raise TimeoutError("timed out waiting for benchmark child state")
        time.sleep(0.001)


def waitid_exit_without_reap(pid: int, deadline: float):
    flags = os.WEXITED | os.WNOHANG | os.WNOWAIT
    while True:
        result = os.waitid(os.P_PID, pid, flags)
        if result is not None:
            require(result.si_pid == pid,
                    "waitid observed an unexpected benchmark child")
            return result
        if time.monotonic() >= deadline:
            raise TimeoutError("timed out waiting for benchmark child exit")
        time.sleep(0.001)


def process_group_or_session_members(identifier: int) -> list[dict[str, Any]]:
    """Return all processes whose pgrp or session is *identifier*.

    The caller retains an unreaped session-leader zombie while using this
    function, so the numeric PID/PGID/SID cannot be recycled during the scan.
    Only normal process-disappearance races are ignored; unreadable proc state
    fails the evidence run closed.
    """
    members = []
    for process in Path("/proc").iterdir():
        if not process.name.isdigit():
            continue
        try:
            text = (process / "stat").read_text()
            end = text.rfind(")")
            require(end >= 0, "malformed /proc process identity")
            fields = text[end + 2:].split()
            require(len(fields) > 19, "short /proc process identity")
            process_group = int(fields[2])
            session = int(fields[3])
            if process_group == identifier or session == identifier:
                members.append({
                    "pid": int(process.name), "state": fields[0],
                    "pgrp": process_group, "session": session,
                    "start_time_ticks": int(fields[19]),
                })
        except (FileNotFoundError, ProcessLookupError):
            continue
    return sorted(members, key=lambda value: value["pid"])


def signal_retained_child_session(
        pid: int, session_proven: bool) -> list[str]:
    """Best-effort SIGKILL while the unreaped leader pins its namespace."""
    failures = []
    operations = [("leader", os.kill, pid)]
    if session_proven:
        operations.append(("process_group", os.killpg, pid))
    for label, operation, identifier in operations:
        try:
            operation(identifier, signal.SIGKILL)
        except ProcessLookupError:
            pass
        except OSError as error:
            failures.append("signal_%s: %s" % (label, error))
    return failures


def require_child_identity_gone_after_reap(pid: int) -> None:
    """Observe only after wait4; never signal a possibly recycled PID/PGID."""
    require(not (Path("/proc") / str(pid)).exists(),
            "benchmark child PID was unexpectedly reused after wait4")


def cleanup_gated_child(
        pid: int, reaped: bool, session_proven: bool,
        timeout: float = 5.0) -> list[str]:
    """Kill and reap a retained child without signaling after PID release."""
    failures = []
    if not reaped:
        for attempt in range(2):
            failures.extend(signal_retained_child_session(
                pid, session_proven))
            try:
                waitid_exit_without_reap(
                    pid, time.monotonic() + timeout)
            except TimeoutError as error:
                failures.append("reap_attempt_%d: %s" % (attempt + 1, error))
                continue
            except (ChildProcessError, EvidenceError, OSError) as error:
                # Child identity ownership is no longer provable after an
                # unexpected wait error, so a numeric PID must not be signaled
                # again.
                failures.append("reap_attempt_%d: %s" % (attempt + 1, error))
                break
            if session_proven:
                quiescence_deadline = time.monotonic() + timeout
                clean = False
                try:
                    while True:
                        members = process_group_or_session_members(pid)
                        clean = len(members) == 1 and \
                            members[0]["pid"] == pid and \
                            members[0]["state"] == "Z" and \
                            members[0]["pgrp"] == pid and \
                            members[0]["session"] == pid
                        if clean:
                            break
                        failures.extend(signal_retained_child_session(
                            pid, True))
                        if time.monotonic() >= quiescence_deadline:
                            failures.append(
                                "retained session did not quiesce: %s" %
                                members)
                            break
                        time.sleep(0.001)
                except Exception as error:
                    failures.append("retained_session_scan: %s" % error)
            try:
                waited_pid, _, _ = os.wait4(pid, 0)
                require(waited_pid == pid,
                        "cleanup wait4 reaped an unexpected child")
                reaped = True
            except (ChildProcessError, OSError) as error:
                failures.append("final_reap: %s" % error)
            break
        if not reaped:
            try:
                waited_pid, _, _ = wait4_until(
                    pid, 0, time.monotonic() + timeout)
                require(waited_pid == pid,
                        "fallback wait4 reaped an unexpected child")
                reaped = True
            except (ChildProcessError, EvidenceError, OSError,
                    TimeoutError) as error:
                failures.append("fallback_reap: %s" % error)
    if reaped:
        try:
            require_child_identity_gone_after_reap(pid)
        except Exception as error:
            failures.append("post_reap_identity: %s" % error)
    else:
        failures.append("leader was not reaped")
    return failures


def run_gated_benchmark(
        command: list[str], cpu: int, sibling: int, timeout: float,
        stdout_path: Path, stderr_path: Path,
        child_signal_mask: set[signal.Signals] | None = None
) -> dict[str, Any]:
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    output_flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
    stdout_fd = os.open(stdout_path, output_flags, 0o600)
    try:
        stderr_fd = os.open(stderr_path, output_flags, 0o600)
    except Exception:
        os.close(stdout_fd)
        raise
    previous_mask = signal.pthread_sigmask(
        signal.SIG_BLOCK, CONTROL_SIGNALS)
    effective_child_mask = previous_mask if child_signal_mask is None \
        else child_signal_mask
    try:
        pid = os.fork()
    except Exception:
        signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)
        os.close(stdout_fd)
        os.close(stderr_fd)
        raise
    if pid == 0:
        try:
            os.dup2(stdout_fd, 1)
            os.dup2(stderr_fd, 2)
            os.close(stdout_fd)
            os.close(stderr_fd)
            os.setsid()
            os.sched_setaffinity(0, {cpu})
            os.kill(os.getpid(), signal.SIGSTOP)
            # Restore the caller's mask only after the parent has observed the
            # stopped, pinned session leader and explicitly continued it.  A
            # control signal received before then remains pending instead of
            # escaping the parent's owned-child cleanup interval.
            signal.pthread_sigmask(signal.SIG_SETMASK, effective_child_mask)
            os.execve(command[0], command, CHILD_ENV)
        except BaseException as error:
            try:
                os.write(2, ("gated benchmark launch failed: %s\n" % error).encode(
                    "utf-8", errors="replace"))
            finally:
                os._exit(126)
    reaped = False
    session_proven = False
    result = None
    try:
        os.close(stdout_fd)
        os.close(stderr_fd)
        stopped_pid, stopped_status, gate_rusage = wait4_until(
            pid, os.WUNTRACED,
            time.monotonic() + min(timeout, 5.0))
        if os.WIFEXITED(stopped_status) or os.WIFSIGNALED(stopped_status):
            reaped = True
        require(not reaped and os.WIFSTOPPED(stopped_status) and
                os.WSTOPSIG(stopped_status) == signal.SIGSTOP,
                "benchmark child failed before the timing gate")
        require(stopped_pid == pid and os.getsid(pid) == pid and
                set(os.sched_getaffinity(pid)) == {cpu},
                "benchmark child gate identity/session/affinity is invalid")
        child_start_time, child_inode = current_task_identity(pid, pid)
        gate_state = current_process_state(pid)
        require(gate_state["state"] in ("T", "t") and
                gate_state["pgrp"] == pid and
                gate_state["session"] == pid and
                gate_state["processor"] == cpu,
                "benchmark child stopped-state identity is invalid")
        session_proven = True
        gate_cpu_time_ns = rusage_cpu_time_ns(gate_rusage)
        child_sched_before = task_scheduler_runtime_ns(pid)

        target_stat_before = cpu_ticks(cpu)
        target_sched_before = cpu_scheduler_runtime_ns(cpu)
        topology = Path(
            "/sys/devices/system/cpu/cpu%d/topology/thread_siblings_list" % cpu)
        sibling_set = parse_cpu_list(topology.read_text()) - {cpu}
        require(sibling_set == {sibling},
                "gated child CPU has a different SMT sibling")
        sibling_stat_before = cpu_ticks(sibling)
        sibling_sched_before = cpu_scheduler_runtime_ns(sibling)
        started_ns = time.time_ns()
        os.kill(pid, signal.SIGCONT)
        timed_out = False
        try:
            waitid_exit_without_reap(pid, time.monotonic() + timeout)
        except TimeoutError:
            timed_out = True
            signal_failures = signal_retained_child_session(pid, True)
            waitid_exit_without_reap(pid, time.monotonic() + 5.0)
            require(not signal_failures,
                    "timed-out benchmark signaling failed: %s" %
                    signal_failures)
        ended_ns = time.time_ns()
        zombie_start_time, zombie_inode = current_task_identity(pid, pid)
        zombie_state = current_process_state(pid)
        require(zombie_start_time == child_start_time and
                zombie_inode == child_inode and
                zombie_state["state"] == "Z" and
                zombie_state["pgrp"] == pid and
                zombie_state["session"] == pid and
                zombie_state["processor"] == cpu and
                set(os.sched_getaffinity(pid)) == {cpu},
                "benchmark zombie identity/session/affinity is invalid")
        retained_session_members = process_group_or_session_members(pid)
        require(retained_session_members == [{
                    "pid": pid, "state": "Z", "pgrp": pid,
                    "session": pid, "start_time_ticks": child_start_time,
                }],
                "benchmark session retained unexpected descendants: %s" %
                retained_session_members)
        child_sched_after = task_scheduler_runtime_ns(pid)
        target_sched_after = cpu_scheduler_runtime_ns(cpu)
        target_stat_after = cpu_ticks(cpu)
        sibling_sched_after = cpu_scheduler_runtime_ns(sibling)
        sibling_stat_after = cpu_ticks(sibling)
        exited_pid, status, exit_rusage = os.wait4(pid, 0)
        reaped = True
        require(exited_pid == pid,
                "wait4 reaped an unexpected benchmark child")
        require_child_identity_gone_after_reap(pid)
        if os.WIFEXITED(status):
            return_code = os.WEXITSTATUS(status)
        elif os.WIFSIGNALED(status):
            return_code = -os.WTERMSIG(status)
        else:
            raise EvidenceError("benchmark child exited in an unexpected state")
        child_runtime_ns = child_sched_after - child_sched_before
        require(child_runtime_ns >= 0,
                "benchmark child schedstat runtime moved backwards")
        wait4_runtime_ns = child_cpu_time_ns(gate_rusage, exit_rusage)
        rusage_difference_ns = abs(wait4_runtime_ns - child_runtime_ns)
        result = {
            "return_code": -1 if timed_out else return_code,
            "timed_out": timed_out,
            "started_ns": started_ns,
            "ended_ns": ended_ns,
            "gate": {
                "pid": pid, "session": pid, "cpu": cpu,
                "start_time_ticks": child_start_time,
                "task_inode": child_inode,
                "stopped_process_state": gate_state,
                "zombie_process_state": zombie_state,
                "retained_session_members": retained_session_members,
                "pre_target_child_cpu_time_ns": gate_cpu_time_ns,
                "child_schedstat_before_ns": child_sched_before,
                "child_schedstat_after_ns": child_sched_after,
                "affinity": [cpu],
            },
            "target_runtime": target_runtime_evidence(
                target_sched_before, target_sched_after, child_runtime_ns,
                TARGET_RUNTIME_TOLERANCE_NS),
            "wait4_crosscheck": {
                "child_cpu_time_ns": wait4_runtime_ns,
                "child_schedstat_delta_ns": child_runtime_ns,
                "absolute_difference_ns": rusage_difference_ns,
                "tolerance_ns": RUSAGE_CROSSCHECK_TOLERANCE_NS,
                "accepted":
                    rusage_difference_ns <= RUSAGE_CROSSCHECK_TOLERANCE_NS,
            },
            "target_interrupts": target_interrupt_evidence(
                target_stat_before, target_stat_after),
            "sibling_runtime": {
                "scheduler_before_ns": sibling_sched_before,
                "scheduler_after_ns": sibling_sched_after,
                "scheduler_delta_ns":
                    sibling_sched_after - sibling_sched_before,
                "accepted": sibling_sched_after == sibling_sched_before,
            },
            "sibling_interrupts": target_interrupt_evidence(
                sibling_stat_before, sibling_stat_after),
        }
    finally:
        try:
            cleanup_failures = cleanup_gated_child(
                pid, reaped, session_proven)
        finally:
            # This is the first parent unmask after fork.  Any pending control
            # signal is delivered only after the retained child/session has
            # been killed if necessary, reaped, and observed gone.
            signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)
        if cleanup_failures:
            raise EvidenceError(
                "benchmark child cleanup failed: %s" % cleanup_failures)
    require(result is not None, "gated benchmark produced no result")
    return result


def run_invocation(binary: Path, implementation: str, cell: dict[str, Any],
                   mode: str,
                   cpu: int, sibling: int, timeout: float,
                   raw_directory: Path, round_index: int,
                   slot_index: int, max_retries: int,
                   reservation_path: Path,
                   reservation: dict[str, Any],
                   isolation_epoch_digest: str,
                   child_signal_mask: set[signal.Signals]) -> dict[str, Any]:
    command = benchmark_command(binary, cell, cpu)
    for attempt in range(max_retries + 1):
        stem = "round%02d-slot%d-%s-attempt%02d" % (
            round_index, slot_index, implementation, attempt)
        reservation_before = reservation_identity(
            reservation_path, cpu, sibling)
        require(reservation_before == reservation,
                "CPU reservation changed before invocation")
        pair_affinity_before = same_user_pair_affinity(cpu, sibling)
        require(not pair_affinity_before,
                "same-user thread can run on reserved pair before invocation")
        stdout_path = raw_directory / (stem + ".stdout.json")
        stderr_path = raw_directory / (stem + ".stderr.txt")
        gated = run_gated_benchmark(
            command, cpu, sibling, timeout, stdout_path, stderr_path,
            child_signal_mask)
        require_no_pending_control_signal()
        stdout = stdout_path.read_bytes()
        stderr = stderr_path.read_bytes()
        reservation_after = reservation_identity(
            reservation_path, cpu, sibling)
        require(reservation_after == reservation,
                "CPU reservation changed during invocation")
        pair_affinity_after = same_user_pair_affinity(cpu, sibling)
        require(not pair_affinity_after,
                "same-user thread can run on reserved pair after invocation")
        sibling_nonidle = nonidle_delta(
            tuple(gated["sibling_interrupts"]["before"]),
            tuple(gated["sibling_interrupts"]["after"]))
        envelope = {
            "schema": ENVELOPE_SCHEMA,
            "implementation": implementation,
            "command": command,
            "started_ns": gated["started_ns"],
            "ended_ns": gated["ended_ns"],
            "return_code": gated["return_code"],
            "stdout": {"path": str(stdout_path), "size": len(stdout),
                       "sha256": sha256_bytes(stdout)},
            "stderr": {"path": str(stderr_path), "size": len(stderr),
                       "sha256": sha256_bytes(stderr)},
            "reserved_sibling_before": gated["sibling_interrupts"]["before"],
            "reserved_sibling_after": gated["sibling_interrupts"]["after"],
            "reserved_sibling_nonidle_jiffies": sibling_nonidle,
            "reservation_before": reservation_before,
            "reservation_after": reservation_after,
            "same_user_pair_affinity_before": pair_affinity_before,
            "same_user_pair_affinity_after": pair_affinity_after,
            "gate": gated["gate"],
            "target_runtime": gated["target_runtime"],
            "target_interrupts": gated["target_interrupts"],
            "sibling_runtime": gated["sibling_runtime"],
            "sibling_interrupts": gated["sibling_interrupts"],
            "wait4_crosscheck": gated["wait4_crosscheck"],
            "isolation_epoch_digest": isolation_epoch_digest,
            "accepted": False,
        }
        isolation_accepted = sibling_nonidle == 0 and \
            gated["target_runtime"]["accepted"] and \
            gated["target_interrupts"]["accepted"] and \
            gated["sibling_runtime"]["accepted"] and \
            gated["sibling_interrupts"]["accepted"] and \
            gated["wait4_crosscheck"]["accepted"]
        if gated["return_code"] == 0 and isolation_accepted:
            try:
                parsed = json.loads(stdout)
                normalized = validate_result(parsed, cell, mode)
                envelope["result"] = parsed
                envelope["normalized"] = normalized
                envelope["accepted"] = True
            except (EvidenceError, json.JSONDecodeError) as error:
                envelope["validation_error"] = str(error)
        envelope_path = raw_directory / (stem + ".envelope.json")
        envelope["envelope_path"] = str(envelope_path)
        atomic_json(envelope_path, envelope)
        if envelope["accepted"]:
            return envelope
        if isolation_accepted:
            reason = envelope.get(
                "validation_error",
                "benchmark exited with %d" % gated["return_code"])
            raise EvidenceError("invalid invocation for %s: %s" % (
                cell["id"], reason))
    raise EvidenceError("no uncontaminated valid invocation for %s" % cell["id"])


def summarize_cell(cell: dict[str, Any], invocations: list[dict[str, Any]],
                   rounds: int, binary_identity: dict[str, Any]) -> dict[str, Any]:
    require(len(invocations) == rounds * len(ORDER),
            "wrong invocation count for cell")
    reference = invocations[0]["normalized"]
    for invocation in invocations[1:]:
        normalized = invocation["normalized"]
        require(normalized["digests"] == reference["digests"],
                "wire/workload digest mismatch between binaries")
        require(normalized["missing_original_indices"] ==
                reference["missing_original_indices"],
                "missing-coordinate mismatch between binaries")
    execution_ratio = ratio_summary(
        invocations, rounds, lambda value: value["median_us"])
    ratios = {"execution": execution_ratio}
    for reuse in (1, 8, 64):
        ratios["amortized_reuse_%d" % reuse] = ratio_summary(
            invocations, rounds,
            lambda value, reuse_count=reuse: value["median_us"] +
                value["plan_setup_us"] / reuse_count)
    result = {
        "schema": CELL_SCHEMA,
        "cell": cell,
        "binary_identity": binary_identity,
        "orientation": "baseline_time_over_candidate_time",
        "round_ratios": execution_ratio["round_ratios"],
        "geometric_mean_ratio": execution_ratio["geometric_mean_ratio"],
        "ci95": execution_ratio["ci95"],
        "metric_ratios": ratios,
        "digests": reference["digests"],
        "missing_original_indices": reference["missing_original_indices"],
        "selected_paths": {
            implementation: {
                "decode_path": next(
                    item["normalized"]["decode_path"] for item in invocations
                    if item["implementation"] == implementation),
                "decode_rule": next(
                    item["normalized"]["decode_rule"] for item in invocations
                    if item["implementation"] == implementation),
                "build_bound_executor": next(
                    item["normalized"]["build_bound_executor"]
                    for item in invocations
                    if item["implementation"] == implementation),
            } for implementation in ("baseline", "candidate")
        },
        "invocations": [
            {
                "implementation": item["implementation"],
                "median_us": item["normalized"]["median_us"],
                "plan_setup_us": item["normalized"]["plan_setup_us"],
                "envelope_path": item["envelope_path"],
                "envelope_sha256": sha256_file(Path(item["envelope_path"])),
                "reserved_sibling_nonidle_jiffies":
                    item["reserved_sibling_nonidle_jiffies"],
                "same_user_pair_affinity_before":
                    item["same_user_pair_affinity_before"],
                "same_user_pair_affinity_after":
                    item["same_user_pair_affinity_after"],
                "target_runtime": item["target_runtime"],
                "target_interrupts": item["target_interrupts"],
                "isolation_epoch_digest": item["isolation_epoch_digest"],
            } for item in invocations
        ],
    }
    result["digest"] = sha256_bytes(canonical_bytes(result))
    return result


def validate_current_binary_identity(identity: dict[str, Any]) -> None:
    for key in BINARY_FILE_IDENTITY_KEYS:
        record = identity.get(key)
        require(isinstance(record, dict) and
                file_identity(Path(record["path"])) == record,
                "current %s identity changed" % key)
    recorded_source = identity.get("source")
    require(isinstance(recorded_source, dict) and
            source_identity(Path(recorded_source["root"])) == recorded_source,
            "current source snapshot changed")
    modes = identity.get("comparison_modes")
    require(isinstance(modes, dict),
            "current comparison modes are missing")
    if is_xmm_tail_comparison(modes):
        for label in ("baseline", "candidate"):
            require(avx2_encoding_audit(Path(
                        identity[label + "_avx2_object"]["path"])) ==
                    identity[label + "_avx2_encoding_audit"],
                    "current %s AVX2 encoding audit changed" % label)
        validate_xmm_tail_comparison_audit(
            identity["xmm_tail_comparison_audit"],
            identity["baseline_avx2_encoding_audit"],
            identity["candidate_avx2_encoding_audit"])
    for label in ("baseline", "candidate"):
        require(build_provenance_identity(
                    Path(identity[label + "_executable"]["path"]),
                    Path(identity[label + "_archive"]["path"]),
                    Path(identity[label + "_compile_commands"]["path"]),
                    Path(recorded_source["root"])) ==
                identity[label + "_build_provenance"],
                "current %s production build provenance changed" % label)


def run_campaign(options: argparse.Namespace) -> int:
    require(options.lock_timeout >= 0,
            "lock timeout must not be negative")
    lock_stream = acquire_lock(options.lock, options.lock_timeout)
    previous_mask = signal.pthread_sigmask(
        signal.SIG_BLOCK, CONTROL_SIGNALS)
    try:
        require_no_pending_control_signal()
        lock_identity = global_lock_identity(lock_stream, options.lock)
        return run_campaign_locked(
            options, lock_stream, lock_identity, previous_mask)
    finally:
        try:
            try:
                fcntl.flock(lock_stream.fileno(), fcntl.LOCK_UN)
            finally:
                lock_stream.close()
        finally:
            # Pending control signals are delivered only after every campaign
            # resource, including the canonical global lock, is released.
            signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)


def run_campaign_locked(
        options: argparse.Namespace, lock_stream: Any,
        lock_identity: dict[str, Any],
        child_signal_mask: set[signal.Signals]) -> int:
    require(options.rounds in (1, 3, 5), "--rounds must be 1, 3, or 5")
    require(0 <= options.max_retries <= 99,
            "--max-retries must be between 0 and 99")
    require(options.timeout > 0 and options.lock_timeout >= 0,
            "benchmark timeout must be positive and lock timeout nonnegative")
    modes = comparison_modes(options.baseline_mode, options.candidate_mode)
    if is_xmm_tail_comparison(modes):
        require(options.preset == "xmm-tail" and not options.cell,
                "XMM-tail comparison requires the complete exact XMM-tail "
                "preset")
    execution_nonce = secrets.token_hex(16)
    topology_identity = require_smt_pair(
        options.cpu, options.reserved_sibling)
    source_root = options.source_root.resolve(strict=True)
    if "production" in modes.values():
        validate_production_mode_source(source_root)
    PairLease, pair_lease_source = load_pair_lease(source_root)
    reservation = reservation_identity(
        options.reservation_file, options.cpu, options.reserved_sibling)
    baseline = options.baseline.resolve(strict=True)
    candidate = options.candidate.resolve(strict=True)
    require(os.access(baseline, os.X_OK) and os.access(candidate, os.X_OK),
            "both benchmark paths must be executable")
    matrix = matrix_for_preset(options.preset)
    selected = matrix["cells"]
    if options.cell:
        identifiers = set(options.cell)
        selected = [cell for cell in selected if cell["id"] in identifiers]
        require(len(selected) == len(identifiers), "unknown --cell identifier")
    validate_comparison_matrix_selection(
        options.preset, modes, matrix, selected)
    baseline_provenance = build_provenance_identity(
        baseline, options.baseline_archive,
        options.baseline_compile_commands, source_root)
    candidate_provenance = build_provenance_identity(
        candidate, options.candidate_archive,
        options.candidate_compile_commands, source_root)
    binary_identity = {
        "comparison_modes": modes,
        "runner": file_identity(Path(__file__)),
        "baseline_executable": file_identity(baseline),
        "baseline_archive": file_identity(options.baseline_archive),
        "baseline_avx2_object": file_identity(options.baseline_avx2_object),
        "baseline_avx2_xor_object": file_identity(
            options.baseline_avx2_xor_object),
        "baseline_compile_commands": file_identity(
            options.baseline_compile_commands),
        "baseline_avx2_compile_entry": avx2_compile_entry(
            options.baseline_compile_commands),
        "baseline_avx2_xor_compile_entry": compile_entry_for(
            options.baseline_compile_commands,
            "Leopard2BackendAVX2Xor.cpp"),
        "baseline_core_compile_entry": compile_entry_for(
            options.baseline_compile_commands, "leopard2.cpp"),
        "baseline_avx512_compile_entry": compile_entry_for(
            options.baseline_compile_commands, "Leopard2BackendAVX512.cpp"),
        "baseline_build_id": options.baseline_build_id,
        "candidate_executable": file_identity(candidate),
        "candidate_archive": file_identity(options.candidate_archive),
        "candidate_avx2_object": file_identity(options.candidate_avx2_object),
        "candidate_avx2_xor_object": file_identity(
            options.candidate_avx2_xor_object),
        "candidate_compile_commands": file_identity(
            options.candidate_compile_commands),
        "candidate_avx2_compile_entry": avx2_compile_entry(
            options.candidate_compile_commands),
        "candidate_avx2_xor_compile_entry": compile_entry_for(
            options.candidate_compile_commands,
            "Leopard2BackendAVX2Xor.cpp"),
        "candidate_core_compile_entry": compile_entry_for(
            options.candidate_compile_commands, "leopard2.cpp"),
        "candidate_avx512_compile_entry": compile_entry_for(
            options.candidate_compile_commands, "Leopard2BackendAVX512.cpp"),
        "candidate_build_id": options.candidate_build_id,
        "source": source_identity(source_root),
        "pair_lease_source": pair_lease_source,
        "baseline_archive_members": baseline_provenance["archive_members"],
        "candidate_archive_members": candidate_provenance["archive_members"],
        "baseline_build_provenance": baseline_provenance,
        "candidate_build_provenance": candidate_provenance,
    }
    binary_identity["baseline_core_object"] = file_identity(
        compile_entry_output(binary_identity["baseline_core_compile_entry"]))
    binary_identity["candidate_core_object"] = file_identity(
        compile_entry_output(binary_identity["candidate_core_compile_entry"]))
    binary_identity["baseline_avx512_object"] = file_identity(
        compile_entry_output(
            binary_identity["baseline_avx512_compile_entry"]))
    binary_identity["candidate_avx512_object"] = file_identity(
        compile_entry_output(
            binary_identity["candidate_avx512_compile_entry"]))
    if is_xmm_tail_comparison(modes):
        binary_identity["baseline_avx2_encoding_audit"] = \
            avx2_encoding_audit(options.baseline_avx2_object)
        binary_identity["candidate_avx2_encoding_audit"] = \
            avx2_encoding_audit(options.candidate_avx2_object)
        binary_identity["xmm_tail_comparison_audit"] = \
            xmm_tail_comparison_audit(
                binary_identity["baseline_avx2_encoding_audit"],
                binary_identity["candidate_avx2_encoding_audit"])
    for label in ("baseline", "candidate"):
        require(Path(binary_identity[label + "_avx2_object"]["path"]) ==
                compile_entry_output(
                    binary_identity[label + "_avx2_compile_entry"]),
                "%s AVX2 object is not the production compile output" % label)
        require(Path(binary_identity[label + "_avx2_xor_object"]["path"]) ==
                compile_entry_output(
                    binary_identity[label + "_avx2_xor_compile_entry"]),
                "%s AVX2 Xor object is not the production compile output" %
                label)
    require(binary_identity["baseline_executable"]["sha256"] !=
            binary_identity["candidate_executable"]["sha256"] and
            binary_identity["baseline_archive"]["sha256"] !=
            binary_identity["candidate_archive"]["sha256"],
            "baseline and candidate archive/executable are byte-identical")
    expected_avx2_source = (
        Path(binary_identity["source"]["root"]) /
        "Leopard2BackendAVX2.cpp").resolve()
    expected_core_source = (
        Path(binary_identity["source"]["root"]) / "leopard2.cpp").resolve()
    expected_avx512_source = (
        Path(binary_identity["source"]["root"]) /
        "Leopard2BackendAVX512.cpp").resolve()
    expected_avx2_xor_source = (
        Path(binary_identity["source"]["root"]) /
        "Leopard2BackendAVX2Xor.cpp").resolve()
    for label in ("baseline", "candidate"):
        avx2_entry = binary_identity[label + "_avx2_compile_entry"]
        avx2_xor_entry = binary_identity[
            label + "_avx2_xor_compile_entry"]
        core_entry = binary_identity[label + "_core_compile_entry"]
        avx512_entry = binary_identity[label + "_avx512_compile_entry"]
        require(compile_entry_source(avx2_entry) == expected_avx2_source,
                "%s AVX2 compile entry uses a different source" % label)
        require(compile_entry_source(avx2_xor_entry) ==
                expected_avx2_xor_source,
                "%s AVX2 Xor compile entry uses a different source" % label)
        require(compile_entry_source(core_entry) == expected_core_source,
                "%s core compile entry uses a different source" % label)
        require(compile_entry_source(avx512_entry) == expected_avx512_source,
                "%s AVX-512 compile entry uses a different source" % label)
    baseline_members = baseline_provenance["object_closure"]["members"]
    candidate_members = candidate_provenance["object_closure"]["members"]
    validate_comparison_object_deltas(
        modes,
        {member: baseline_members[member]["compile"]["output"]["sha256"]
         for member in EXPECTED_ARCHIVE_MEMBERS},
        {member: candidate_members[member]["compile"]["output"]["sha256"]
         for member in EXPECTED_ARCHIVE_MEMBERS})
    require(baseline_provenance["object_closure"]["benchmark"]["output"]
                ["sha256"] ==
            candidate_provenance["object_closure"]["benchmark"]["output"]
                ["sha256"],
            "diagnostic unexpectedly changed the benchmark object")
    baseline_cache = baseline_provenance["cache_values"]
    candidate_cache = candidate_provenance["cache_values"]
    require(shlex.split(baseline_cache["CMAKE_CXX_FLAGS"]) ==
                mode_compile_arguments(modes["baseline"]) and
            shlex.split(candidate_cache["CMAKE_CXX_FLAGS"]) ==
                mode_compile_arguments(modes["candidate"]) and
            baseline_cache["CMAKE_CXX_FLAGS_RELEASE"] ==
                candidate_cache["CMAKE_CXX_FLAGS_RELEASE"] ==
                "-O3 -DNDEBUG",
            "diagnostic build flags do not bind the requested modes")
    for key in baseline_cache:
        if key != "CMAKE_CXX_FLAGS":
            require(baseline_cache[key] == candidate_cache[key],
                    "diagnostic CMake cache differs at %s" % key)
    require(baseline_provenance["link_recipes"]["archive"]["text"] ==
                candidate_provenance["link_recipes"]["archive"]["text"],
            "diagnostic archive link recipe differs")
    baseline_link = shlex.split(
        baseline_provenance["link_recipes"]["executable"]["text"])
    candidate_link = shlex.split(
        candidate_provenance["link_recipes"]["executable"]["text"])
    baseline_link = strip_mode_definition(
        baseline_link, modes["baseline"], "baseline executable link")
    candidate_link = strip_mode_definition(
        candidate_link, modes["candidate"], "candidate executable link")
    require(baseline_link == candidate_link,
            "diagnostic executable link recipe differs beyond the modes")
    for tool in ("archive_tool", "ranlib_tool", "compiler"):
        require(baseline_provenance["link_recipes"][tool]["sha256"] ==
                    candidate_provenance["link_recipes"][tool]["sha256"],
                "diagnostic build tool differs at %s" % tool)
    for member in EXPECTED_ARCHIVE_MEMBERS:
        baseline_arguments = list(
            baseline_members[member]["compile"]["arguments"])
        candidate_arguments = list(
            candidate_members[member]["compile"]["arguments"])
        baseline_arguments = strip_mode_definition(
            baseline_arguments, modes["baseline"],
            "baseline compile %s" % member)
        candidate_arguments = strip_mode_definition(
            candidate_arguments, modes["candidate"],
            "candidate compile %s" % member)
        require(baseline_arguments == candidate_arguments,
                "production compile commands differ beyond the modes: %s" %
                    member)
    baseline_benchmark_arguments = list(
        baseline_provenance["object_closure"]["benchmark"]["arguments"])
    candidate_benchmark_arguments = list(
        candidate_provenance["object_closure"]["benchmark"]["arguments"])
    baseline_benchmark_arguments = strip_mode_definition(
        baseline_benchmark_arguments, modes["baseline"],
        "baseline benchmark compile")
    candidate_benchmark_arguments = strip_mode_definition(
        candidate_benchmark_arguments, modes["candidate"],
        "candidate benchmark compile")
    require(baseline_benchmark_arguments == candidate_benchmark_arguments,
            "benchmark compile commands differ beyond the modes")
    # Replay the complete retained structure before launching the first timed
    # child.  Offline verification repeats this check from manifest bytes, but
    # an invalid build must fail before consuming an isolation lease or
    # producing any timing observations.
    validate_binary_identity_structure(binary_identity)
    request = {
        "schema": comparison_schema(modes),
        "matrix_sha256": matrix["matrix_sha256"],
        "matrix_preset": options.preset,
        "selected_cell_ids": [cell["id"] for cell in selected],
        "rounds": options.rounds,
        "order": ORDER,
        "binary_identity": binary_identity,
        "cpu": options.cpu,
        "reserved_sibling": options.reserved_sibling,
        "reservation": reservation,
        "topology_identity": topology_identity,
        "isolation_policy": isolation_policy(),
        "timeout_seconds": options.timeout,
        "max_retries": options.max_retries,
        "execution_nonce": execution_nonce,
        "resume_policy": "fresh complete authoritative campaign only",
        "lock": lock_identity,
    }
    request["request_digest"] = sha256_bytes(canonical_bytes(request))
    output = options.output.resolve()
    try:
        output.relative_to(Path(binary_identity["source"]["root"]))
    except ValueError:
        pass
    else:
        raise EvidenceError(
            "campaign output must be outside the source worktree")
    output.mkdir(parents=True, exist_ok=True)
    request_path = output / "request.json"
    require(not request_path.exists() and
            not (output / "manifest.json").exists() and
            not (output / "cells").exists() and
            not (output / "raw").exists(),
            "authoritative ABBA output is not empty; use a fresh directory")
    atomic_json(request_path, request)
    atomic_json(output / "matrix.json", matrix)
    pair_guard = None
    affinity_exclusion = None
    affinity_restoration = None
    try:
        try:
            pair_guard = PairLease(options.cpu, options.reserved_sibling)
            pair_guard.__enter__()
            pair_lease_identity = pair_guard.revalidate()
        except Exception as error:
            raise EvidenceError(
                "cannot acquire canonical CPU-pair lease: %s" % error) from error
        affinity_exclusion = exclude_same_user_from_pair(
            options.cpu, options.reserved_sibling)
        isolation_epoch_start = {
            "schema": "leopard2-direct-isolation-epoch/v4",
            "execution_nonce": execution_nonce,
            "global_lock": lock_identity,
            "pair_lease": pair_lease_identity,
            "reservation": reservation,
            "topology_identity": topology_identity,
            "affinity_exclusion": affinity_exclusion,
        }
        isolation_epoch_digest = sha256_bytes(
            canonical_bytes(isolation_epoch_start))
        summaries = []
        for cell_index, cell in enumerate(selected):
            require_no_pending_control_signal()
            cell_path = output / "cells" / (cell["id"] + ".json")
            require(not cell_path.exists(),
                    "fresh authoritative campaign contains a retained cell")
            invocations = []
            raw_directory = output / "raw" / cell["id"]
            for round_index in range(options.rounds):
                for slot_index, implementation in enumerate(ORDER):
                    binary = baseline if implementation == "baseline" else candidate
                    require(global_lock_identity(
                                lock_stream, options.lock) == lock_identity,
                            "canonical GF8 lock identity changed")
                    require(pair_guard.revalidate() == pair_lease_identity,
                            "canonical CPU-pair lease changed before invocation")
                    invocations.append(run_invocation(
                        binary, implementation, cell, modes[implementation],
                        options.cpu,
                        options.reserved_sibling, options.timeout, raw_directory,
                        round_index, slot_index, options.max_retries,
                        options.reservation_file, reservation,
                        isolation_epoch_digest, child_signal_mask))
                    require_no_pending_control_signal()
                    require(pair_guard.revalidate() == pair_lease_identity,
                            "canonical CPU-pair lease changed during invocation")
            summary = summarize_cell(
                cell, invocations, options.rounds, binary_identity)
            atomic_json(cell_path, summary)
            summaries.append(summary)
            if options.rounds == 1:
                print("[%d/%d] %s %.6fx directional" % (
                    cell_index + 1, len(selected), cell["id"],
                    summary["geometric_mean_ratio"]), flush=True)
            else:
                print("[%d/%d] %s %.6fx [%.6f, %.6f]" % (
                    cell_index + 1, len(selected), cell["id"],
                    summary["geometric_mean_ratio"], summary["ci95"][0],
                    summary["ci95"][1]), flush=True)
        validate_current_binary_identity(binary_identity)
        require(global_lock_identity(lock_stream, options.lock) ==
                lock_identity,
                "canonical GF8 lock identity changed before manifest")
        require(reservation_identity(
                    options.reservation_file, options.cpu,
                    options.reserved_sibling) == reservation,
                "CPU reservation changed before manifest")
        require(pair_guard.revalidate() == pair_lease_identity,
                "canonical CPU-pair lease changed before manifest")
        affinity_restoration = restore_same_user_affinity(affinity_exclusion)
        isolation_epoch = {
            "start": isolation_epoch_start,
            "start_digest": isolation_epoch_digest,
            "restoration": affinity_restoration,
        }
        manifest = {
            "schema": comparison_schema(modes),
            "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "request": request,
            "global_lock": lock_identity,
            "pair_lease": pair_lease_identity,
            "affinity_exclusion": affinity_exclusion,
            "isolation_epoch": isolation_epoch,
            "host": {"node": platform.node(), "platform": platform.platform(),
                     "machine": platform.machine(),
                     "allowed_cpus": sorted(os.sched_getaffinity(0))},
            "statistics": {
                "method": "one mean log contrast per independent ABBA round",
                "confidence": 0.95,
                "rounds": options.rounds,
                "degrees_of_freedom": options.rounds - 1,
                "directional_only": options.rounds == 1,
                "child_estimator": "median of retained per-invocation samples",
            },
            "cell_summaries": summaries,
        }
        manifest["digest"] = sha256_bytes(canonical_bytes(manifest))
        validate_current_binary_identity(binary_identity)
        require(global_lock_identity(lock_stream, options.lock) ==
                lock_identity,
                "canonical GF8 lock identity changed at manifest commit")
        require(reservation_identity(
                    options.reservation_file, options.cpu,
                    options.reserved_sibling) == reservation,
                "CPU reservation changed at manifest commit")
        require(pair_guard.revalidate() == pair_lease_identity,
                "canonical CPU-pair lease changed at manifest commit")
        atomic_json(output / "manifest.json", manifest)
        require(global_lock_identity(lock_stream, options.lock) ==
                lock_identity,
                "canonical GF8 lock identity changed after manifest commit")
    finally:
        cleanup_error = None
        if affinity_exclusion is not None and affinity_restoration is None:
            try:
                restore_same_user_affinity(affinity_exclusion)
            except Exception as error:
                cleanup_error = error
        if pair_guard is not None:
            try:
                pair_guard.__exit__(None, None, None)
            except Exception as error:
                if cleanup_error is None:
                    cleanup_error = error
        if cleanup_error is not None:
            raise EvidenceError(
                "benchmark isolation cleanup failed: %s" % cleanup_error) \
                from cleanup_error
    return 0


def require_exact_keys(value: Any, keys: set[str], label: str) -> None:
    require(isinstance(value, dict) and set(value) == keys,
            "%s has unexpected structure" % label)


def require_int(value: Any, label: str, minimum: int = 0) -> int:
    require(type(value) is int and value >= minimum,
            "%s is not a valid integer" % label)
    return value


def validate_topology_identity(
        value: Any, cpu: int, sibling: int) -> None:
    require_exact_keys(value, {
        "cgroup_effective_cpus", "launch_affinity", "housekeeping_affinity",
    }, "topology identity")
    cgroup = value["cgroup_effective_cpus"]
    launch = value["launch_affinity"]
    housekeeping = value["housekeeping_affinity"]
    for name, cpus in (("cgroup", cgroup), ("launch", launch),
                       ("housekeeping", housekeeping)):
        require(isinstance(cpus, list) and
                all(type(item) is int and item >= 0 for item in cpus) and
                cpus == sorted(set(cpus)),
                "%s CPU identity is not canonical" % name)
    require(cpu in cgroup and sibling in cgroup,
            "pair is outside recorded cgroup eligibility")
    require(housekeeping == sorted(set(launch) - {cpu, sibling}) and
            housekeeping,
            "recorded inherited/housekeeping affinity is inconsistent")


def validate_pair_lease_record(value: Any, cpu: int, sibling: int) -> None:
    require_exact_keys(value, {
        "lock", "path", "payload", "sha256", "device", "inode",
        "directory_device", "directory_inode", "kernel_lease",
    }, "CPU-pair lease")
    expected_lock = (
        "dual_linux_abstract_af_unix_bind_and_fcntl_flock_"
        "exclusive_nonblocking_pair_wide")
    require(value["lock"] == expected_lock and
            isinstance(value["path"], str) and
            Path(value["path"]).is_absolute() and
            isinstance(value["sha256"], str) and
            re.fullmatch(r"[0-9a-f]{64}", value["sha256"]) is not None,
            "CPU-pair lease strings are invalid")
    for key in ("device", "inode", "directory_device", "directory_inode"):
        require_int(value[key], "CPU-pair lease %s" % key)
    payload = value["payload"]
    require_exact_keys(payload, {"cpus", "schema", "uid"},
                       "CPU-pair lease payload")
    require(payload["cpus"] == sorted((cpu, sibling)) and
            payload["schema"] == "leopard2-cpu-pair-lease/v1" and
            type(payload["uid"]) is int and payload["uid"] >= 0,
            "CPU-pair lease payload is invalid")
    expected_path = Path("/run/user") / str(payload["uid"]) / \
        "leopard2-cpu-leases" / (
            "leopard2-cpu-pair-%d-%d-%d.lock" %
            (payload["uid"], min(cpu, sibling), max(cpu, sibling)))
    require(Path(value["path"]) == expected_path and
            value["sha256"] == sha256_bytes(canonical_bytes(payload)),
            "CPU-pair lease path/payload digest is invalid")
    kernel = value["kernel_lease"]
    require_exact_keys(kernel, {"schema", "mechanism", "name_sha256"},
                       "kernel CPU-pair lease")
    lease_root = expected_path.parent.parent
    material = canonical_bytes({
        "cpus": sorted((cpu, sibling)),
        "root": os.path.abspath(expected_path.parent),
        "schema": "leopard2-cpu-pair-lease/v1",
        "uid": payload["uid"],
    })
    kernel_name = b"\0leopard2-pair-v1-" + \
        hashlib.sha256(material).hexdigest()[:40].encode("ascii")
    require(kernel["schema"] == "leopard2-kernel-cpu-pair-lease/v1" and
            kernel["mechanism"] ==
                "linux-abstract-af-unix-bind-exclusive" and
            kernel["name_sha256"] == hashlib.sha256(kernel_name).hexdigest() and
            lease_root == Path("/run/user") / str(payload["uid"]),
            "kernel CPU-pair lease is invalid")


def validate_affinity_exclusion(
        exclusion: Any, restoration: Any, cpu: int, sibling: int) -> None:
    require_exact_keys(exclusion, {
        "schema", "uid", "pair", "changed", "after_snapshot",
        "started_ns", "completed_ns", "quiescence_passes",
        "runner_before", "runner_after",
    }, "affinity exclusion")
    require(exclusion["schema"] == "leopard2-user-affinity-exclusion/v1" and
            exclusion["pair"] == sorted((cpu, sibling)) and
            type(exclusion["uid"]) is int and exclusion["uid"] >= 0 and
            exclusion["after_snapshot"] == [],
            "affinity exclusion header is invalid")
    started = require_int(exclusion["started_ns"], "exclusion start", 1)
    completed = require_int(exclusion["completed_ns"], "exclusion end", 1)
    require(completed >= started and
            1 <= require_int(exclusion["quiescence_passes"],
                             "quiescence passes", 1) <= 8,
            "affinity exclusion timing/pass count is invalid")
    changed = exclusion["changed"]
    require(isinstance(changed, list), "affinity changed list is invalid")
    retained = {}
    pair = {cpu, sibling}
    for record in changed:
        require_exact_keys(record, {
            "pid", "tid", "command", "affinity", "start_time_ticks",
            "task_inode", "after",
        }, "affinity exclusion task")
        key = (require_int(record["pid"], "excluded PID", 1),
               require_int(record["tid"], "excluded TID", 1))
        require(key not in retained, "duplicate excluded task identity")
        before = record["affinity"]
        after = record["after"]
        require(isinstance(record["command"], str) and record["command"] and
                isinstance(before, list) and before == sorted(set(before)) and
                isinstance(after, list) and after == sorted(set(after)) and
                bool(set(before) & pair) and not (set(after) & pair) and after,
                "excluded task masks are invalid")
        require_int(record["start_time_ticks"], "task start ticks", 1)
        require_int(record["task_inode"], "task inode", 1)
        retained[key] = before
    require_exact_keys(restoration, {
        "schema", "restored", "failures", "runner_restored",
    },
                       "affinity restoration")
    require(restoration["schema"] ==
                "leopard2-user-affinity-restoration/v1" and
            restoration["failures"] == [] and
            isinstance(restoration["restored"], list),
            "affinity restoration header is invalid")
    restored = {}
    for record in restoration["restored"]:
        require_exact_keys(record, {
            "pid", "tid", "affinity", "start_time_ticks", "task_inode",
        },
                           "restored task")
        key = (require_int(record["pid"], "restored PID", 1),
               require_int(record["tid"], "restored TID", 1))
        source = next((item for item in changed
                       if (item["pid"], item["tid"]) == key), None)
        require(key not in restored and key in retained and source is not None and
                record["affinity"] == retained[key] and
                record["start_time_ticks"] == source["start_time_ticks"] and
                record["task_inode"] == source["task_inode"],
                "restored task does not match excluded mask")
        restored[key] = record["affinity"]
    require(restored == retained,
            "not every excluded task has an exact restoration record")
    require_exact_keys(restoration["runner_restored"], {
        "pid", "tid", "affinity", "start_time_ticks", "task_inode",
    }, "restored runner")
    for name in ("runner_before", "runner_after"):
        require_exact_keys(exclusion[name], {
            "pid", "tid", "affinity", "start_time_ticks", "task_inode",
        }, name)
    runner_before = exclusion["runner_before"]
    runner_after = exclusion["runner_after"]
    runner_restored = restoration["runner_restored"]
    require((runner_before["pid"], runner_before["tid"],
             runner_before["start_time_ticks"], runner_before["task_inode"]) ==
            (runner_after["pid"], runner_after["tid"],
             runner_after["start_time_ticks"], runner_after["task_inode"]) ==
            (runner_restored["pid"], runner_restored["tid"],
             runner_restored["start_time_ticks"],
             runner_restored["task_inode"]) and
            not (set(runner_after["affinity"]) & pair) and
            runner_restored["affinity"] == runner_before["affinity"],
            "runner affinity identity/masks were not exactly restored")


def validate_file_identity_record(value: Any, label: str) -> None:
    require_exact_keys(value, {"path", "size", "sha256"}, label)
    require(isinstance(value["path"], str) and
            Path(value["path"]).is_absolute() and
            type(value["size"]) is int and value["size"] >= 0 and
            isinstance(value["sha256"], str) and
            re.fullmatch(r"[0-9a-f]{64}", value["sha256"]) is not None,
            "%s identity fields are invalid" % label)


def validate_global_lock_record(value: Any) -> None:
    require_exact_keys(value, {
        "schema", "path", "device", "inode", "uid", "mode", "nlink",
        "mechanism",
    }, "global benchmark lock")
    require(value["schema"] == "leopard2-global-benchmark-lock/v1" and
            value["path"] == str(DEFAULT_LOCK) and
            value["mechanism"] == "fcntl-flock-exclusive" and
            value["mode"] == 0o600 and value["nlink"] == 1,
            "global benchmark lock identity is not canonical")
    for key in ("device", "inode", "uid", "mode", "nlink"):
        require_int(value[key], "global benchmark lock %s" % key)


def validate_reservation_record(
        value: Any, cpu: int, sibling: int) -> None:
    require_exact_keys(value, {"file", "value"}, "CPU reservation")
    validate_file_identity_record(value["file"], "CPU reservation file")
    reservation = value["value"]
    require_exact_keys(reservation, {
        "benchmark_cpu", "nonce", "owner", "reserved_sibling",
        "schema", "status",
    }, "CPU reservation payload")
    require(reservation["schema"] == "leopard2-cpu-reservation/v1" and
            reservation["status"] == "held" and
            reservation["benchmark_cpu"] == cpu and
            reservation["reserved_sibling"] == sibling and
            isinstance(reservation["nonce"], str) and reservation["nonce"] and
            isinstance(reservation["owner"], str) and reservation["owner"],
            "CPU reservation payload is invalid")


def validate_interrupt_record(value: Any, label: str) -> None:
    require_exact_keys(value, {
        "before", "after", "deltas", "rejected_fields", "accepted",
    }, label)
    before = value["before"]
    after = value["after"]
    require(isinstance(before, list) and isinstance(after, list) and
            len(before) == len(after) and len(before) >= 5 and
            all(type(item) is int and item >= 0 for item in before + after),
            "%s /proc/stat snapshots are invalid" % label)
    expected = json.loads(canonical_bytes(target_interrupt_evidence(
        tuple(before), tuple(after))))
    require(value == expected and value["accepted"] is True,
            "%s arithmetic or acceptance is invalid" % label)


def validate_gate_record(value: Any, cpu: int) -> int:
    require_exact_keys(value, {
        "pid", "session", "cpu", "start_time_ticks", "task_inode",
        "stopped_process_state", "zombie_process_state",
        "retained_session_members", "pre_target_child_cpu_time_ns",
        "child_schedstat_before_ns", "child_schedstat_after_ns", "affinity",
    }, "gated child")
    pid = require_int(value["pid"], "gated child PID", 1)
    require(value["session"] == pid and value["cpu"] == cpu and
            value["affinity"] == [cpu],
            "gated child identity/affinity is invalid")
    start_time = require_int(
        value["start_time_ticks"], "gated child start time", 1)
    require_int(value["task_inode"], "gated child task inode", 1)
    require_int(value["pre_target_child_cpu_time_ns"],
                "gated child pre-target CPU time")
    before = require_int(value["child_schedstat_before_ns"],
                         "gated child schedstat before")
    after = require_int(value["child_schedstat_after_ns"],
                        "gated child schedstat after")
    require(after >= before, "gated child schedstat moved backwards")
    for key, expected_state in (("stopped_process_state", ("T", "t")),
                                ("zombie_process_state", ("Z",))):
        state = value[key]
        require_exact_keys(state, {"state", "pgrp", "session", "processor"},
                           key)
        require(state["state"] in expected_state and
                state["pgrp"] == pid and state["session"] == pid and
                state["processor"] == cpu,
                "%s is invalid" % key)
    members = value["retained_session_members"]
    require(members == [{
                "pid": pid, "state": "Z", "pgrp": pid, "session": pid,
                "start_time_ticks": start_time,
            }],
            "retained child session was not leader-only")
    return after - before


def validate_runtime_records(
        envelope: dict[str, Any], child_runtime_ns: int) -> None:
    target = envelope["target_runtime"]
    require_exact_keys(target, {
        "scheduler_before_ns", "scheduler_after_ns", "scheduler_delta_ns",
        "child_cpu_time_ns", "signed_difference_ns",
        "unexplained_runtime_ns", "tolerance_ns", "accepted",
    }, "target runtime")
    before = require_int(target["scheduler_before_ns"],
                         "target scheduler before")
    after = require_int(target["scheduler_after_ns"],
                        "target scheduler after")
    expected_target = target_runtime_evidence(
        before, after, child_runtime_ns, TARGET_RUNTIME_TOLERANCE_NS)
    require(target == expected_target and target["accepted"] is True,
            "target runtime arithmetic or acceptance is invalid")

    sibling = envelope["sibling_runtime"]
    require_exact_keys(sibling, {
        "scheduler_before_ns", "scheduler_after_ns", "scheduler_delta_ns",
        "accepted",
    }, "sibling runtime")
    sibling_before = require_int(
        sibling["scheduler_before_ns"], "sibling scheduler before")
    sibling_after = require_int(
        sibling["scheduler_after_ns"], "sibling scheduler after")
    require(sibling == {
                "scheduler_before_ns": sibling_before,
                "scheduler_after_ns": sibling_after,
                "scheduler_delta_ns": sibling_after - sibling_before,
                "accepted": sibling_after == sibling_before,
            } and sibling["accepted"] is True,
            "sibling runtime arithmetic or acceptance is invalid")

    crosscheck = envelope["wait4_crosscheck"]
    require_exact_keys(crosscheck, {
        "child_cpu_time_ns", "child_schedstat_delta_ns",
        "absolute_difference_ns", "tolerance_ns", "accepted",
    }, "wait4 cross-check")
    wait4_runtime = require_int(
        crosscheck["child_cpu_time_ns"], "wait4 child CPU time")
    expected_difference = abs(wait4_runtime - child_runtime_ns)
    require(crosscheck == {
                "child_cpu_time_ns": wait4_runtime,
                "child_schedstat_delta_ns": child_runtime_ns,
                "absolute_difference_ns": expected_difference,
                "tolerance_ns": RUSAGE_CROSSCHECK_TOLERANCE_NS,
                "accepted":
                    expected_difference <= RUSAGE_CROSSCHECK_TOLERANCE_NS,
            } and crosscheck["accepted"] is True,
            "wait4 cross-check arithmetic or acceptance is invalid")


def validate_retained_stream(
        value: Any, expected_path: Path, label: str) -> bytes:
    validate_file_identity_record(value, label)
    require(Path(value["path"]) == expected_path and
            expected_path.resolve(strict=True) == expected_path and
            file_identity(expected_path) == value,
            "%s retained identity changed" % label)
    return expected_path.read_bytes()


def validate_envelope(
        value: Any, envelope_path: Path, cell: dict[str, Any],
        implementation: str, round_index: int, slot_index: int,
        request: dict[str, Any], isolation_epoch_digest: str,
        campaign_root: Path) -> dict[str, Any]:
    require_exact_keys(value, {
        "schema", "implementation", "command", "started_ns", "ended_ns",
        "return_code", "stdout", "stderr", "reserved_sibling_before",
        "reserved_sibling_after", "reserved_sibling_nonidle_jiffies",
        "reservation_before", "reservation_after",
        "same_user_pair_affinity_before", "same_user_pair_affinity_after",
        "gate", "target_runtime", "target_interrupts", "sibling_runtime",
        "sibling_interrupts", "wait4_crosscheck",
        "isolation_epoch_digest", "accepted", "result", "normalized",
        "envelope_path",
    }, "retained invocation envelope")
    require(value["schema"] == ENVELOPE_SCHEMA and
            value["implementation"] == implementation and
            value["return_code"] == 0 and value["accepted"] is True and
            value["envelope_path"] == str(envelope_path),
            "retained invocation envelope header is invalid")
    started = require_int(value["started_ns"], "invocation start", 1)
    ended = require_int(value["ended_ns"], "invocation end", 1)
    require(ended >= started, "invocation timing moved backwards")
    expected_binary = request["binary_identity"][
        implementation + "_executable"]["path"]
    require(value["command"] == benchmark_command(
                Path(expected_binary), cell, request["cpu"]),
            "retained invocation command does not match the frozen cell")
    require(value["reservation_before"] == request["reservation"] and
            value["reservation_after"] == request["reservation"] and
            value["same_user_pair_affinity_before"] == [] and
            value["same_user_pair_affinity_after"] == [] and
            value["isolation_epoch_digest"] == isolation_epoch_digest,
            "retained invocation is outside the recorded isolation epoch")

    expected_directory = campaign_root / "raw" / cell["id"]
    require(envelope_path.parent == expected_directory and
            envelope_path.resolve(strict=True) == envelope_path,
            "retained invocation envelope escaped its cell directory")
    expected_prefix = "round%02d-slot%d-%s-attempt" % (
        round_index, slot_index, implementation)
    match = re.fullmatch(
        re.escape(expected_prefix) + r"([0-9]{2})\.envelope\.json",
        envelope_path.name)
    require(match is not None and
            int(match.group(1)) <= request["max_retries"],
            "retained invocation envelope has the wrong ABBA slot/attempt")
    stem = envelope_path.name[:-len(".envelope.json")]
    stdout_path = expected_directory / (stem + ".stdout.json")
    stderr_path = expected_directory / (stem + ".stderr.txt")
    stdout = validate_retained_stream(value["stdout"], stdout_path, "stdout")
    validate_retained_stream(value["stderr"], stderr_path, "stderr")
    try:
        parsed = json.loads(stdout)
    except json.JSONDecodeError as error:
        raise EvidenceError("retained stdout is not JSON") from error
    require(parsed == value["result"],
            "retained stdout and embedded result differ")
    mode = request["binary_identity"]["comparison_modes"][implementation]
    normalized = validate_result(parsed, cell, mode)
    require(normalized == value["normalized"],
            "retained normalized result was forged")

    child_runtime = validate_gate_record(value["gate"], request["cpu"])
    validate_runtime_records(value, child_runtime)
    validate_interrupt_record(value["target_interrupts"],
                              "target interrupts")
    validate_interrupt_record(value["sibling_interrupts"],
                              "sibling interrupts")
    require(value["reserved_sibling_before"] ==
                value["sibling_interrupts"]["before"] and
            value["reserved_sibling_after"] ==
                value["sibling_interrupts"]["after"] and
            value["reserved_sibling_nonidle_jiffies"] == nonidle_delta(
                tuple(value["reserved_sibling_before"]),
                tuple(value["reserved_sibling_after"])) == 0,
            "reserved sibling jiffy evidence is invalid")
    return value


def validate_archive_members_structure(
        value: Any, archive: dict[str, Any], label: str) -> None:
    require_exact_keys(value, {
        "schema", "path", "members", "member_files", "archiver", "sha256",
    }, label)
    validate_file_identity_record(value["archiver"], label + " archiver")
    require(value["schema"] == "leopard2-archive-members/v2" and
            value["path"] == archive["path"] and
            value["members"] == list(EXPECTED_ARCHIVE_MEMBERS) and
            isinstance(value["member_files"], dict) and
            set(value["member_files"]) == set(EXPECTED_ARCHIVE_MEMBERS),
            label + " closure/order is invalid")
    for member in EXPECTED_ARCHIVE_MEMBERS:
        record = value["member_files"][member]
        require_exact_keys(record, {"size", "sha256"},
                           label + " " + member)
        require(type(record["size"]) is int and record["size"] > 0 and
                isinstance(record["sha256"], str) and
                re.fullmatch(r"[0-9a-f]{64}", record["sha256"]) is not None,
                label + " member bytes are invalid")
    copy = dict(value)
    digest = copy.pop("sha256")
    require(digest == sha256_bytes(canonical_bytes(copy)),
            label + " digest is invalid")


def validate_compile_binding_structure(
        value: Any, build: Path, root: Path, source: str, target: str,
        compiler: dict[str, Any], label: str) -> None:
    require_exact_keys(value, {
        "schema", "source", "target", "entry", "arguments", "output",
        "sha256",
    }, label)
    require(value["schema"] == "leopard2-production-compile-entry/v1" and
            value["source"] == str(root / source) and
            value["target"] == target and
            isinstance(value["entry"], dict) and
            isinstance(value["arguments"], list) and value["arguments"] and
            all(isinstance(argument, str) and "@" not in argument
                for argument in value["arguments"]),
            label + " header is invalid")
    validate_file_identity_record(value["output"], label + " object")
    output = Path(value["output"]["path"])
    require(output.is_absolute() and output.is_relative_to(build) and
            "CMakeFiles/%s/" % target in output.as_posix() and
            Path(value["arguments"][0]).resolve(strict=True) ==
                Path(compiler["path"]),
            label + " compiler/output target is invalid")
    entry = value["entry"]
    require(set(entry) in ({"directory", "command", "file", "output"},
                           {"arguments", "directory", "file", "output"}) and
            Path(str(entry["directory"])) == build and
            Path(str(entry["file"])) == root / source,
            label + " raw compile entry is invalid")
    arguments = value["arguments"]
    require(compile_entry_arguments(entry) == arguments,
            label + " retained argv differs from raw compile entry")
    validate_compile_profile(arguments, target)
    require(arguments.count("-c") == 1 and
            arguments.index("-c") + 1 < len(arguments) and
            arguments[arguments.index("-c") + 1] == str(root / source) and
            arguments.count(str(root / source)) == 1 and
            [Path(argument) for argument in arguments
             if argument.endswith(".cpp")] == [root / source] and
            arguments.count("-o") == 1 and
            arguments.index("-o") + 1 < len(arguments),
            label + " compile/source operands are invalid")
    raw_output = Path(arguments[arguments.index("-o") + 1])
    if not raw_output.is_absolute():
        raw_output = build / raw_output
    entry_output = Path(str(entry["output"]))
    if not entry_output.is_absolute():
        entry_output = build / entry_output
    require(raw_output == output and entry_output == output,
            label + " compile output does not match the object identity")
    copy = dict(value)
    digest = copy.pop("sha256")
    require(digest == sha256_bytes(canonical_bytes(copy)),
            label + " digest is invalid")


def validate_exact_text_structure(value: Any, schema: str, label: str) -> None:
    require_exact_keys(value, {
        "schema", "file", "encoding", "text", "sha256",
    }, label)
    validate_file_identity_record(value["file"], label + " file")
    require(value["schema"] == schema and value["encoding"] == "utf-8" and
            isinstance(value["text"], str), label + " header is invalid")
    payload = value["text"].encode("utf-8")
    require(value["file"]["size"] == len(payload) and
            value["file"]["sha256"] == sha256_bytes(payload),
            label + " retained bytes differ from file identity")
    copy = dict(value)
    digest = copy.pop("sha256")
    require(digest == sha256_bytes(canonical_bytes(copy)),
            label + " digest is invalid")


def validate_build_provenance_structure(
        value: Any, source: dict[str, Any], label: str) -> None:
    require_exact_keys(value, {
        "schema", "build", "executable", "archive", "cache",
        "cache_content",
        "cache_values", "archive_members", "object_closure",
        "link_recipes", "sha256",
    }, label)
    require(value["schema"] == "leopard2-production-build-provenance/v1" and
            isinstance(value["build"], str) and
            Path(value["build"]).is_absolute() and
            isinstance(value["cache_values"], dict),
            label + " header is invalid")
    for name in ("executable", "archive", "cache"):
        validate_file_identity_record(value[name], label + " " + name)
    build = Path(value["build"])
    root = Path(source["root"])
    require(Path(value["executable"]["path"]) == build / "bench_leopard2" and
            Path(value["archive"]["path"]) == build / "libleopard.a" and
            Path(value["cache"]["path"]) == build / "CMakeCache.txt",
            label + " canonical build artifact paths are invalid")
    validate_exact_text_structure(
        value["cache_content"], "leopard2-cmake-cache-content/v1",
        label + " CMake cache content")
    require(value["cache_content"]["file"] == value["cache"],
            label + " retained CMake cache bytes are not cross-bound")
    cache = value["cache_values"]
    require(set(cache) == {
                "CMAKE_AR", "CMAKE_BUILD_TYPE", "CMAKE_CXX_COMPILER",
                "CMAKE_CXX_FLAGS", "CMAKE_CXX_FLAGS_RELEASE",
                "CMAKE_EXE_LINKER_FLAGS", "CMAKE_EXE_LINKER_FLAGS_RELEASE",
                "CMAKE_EXPORT_COMPILE_COMMANDS", "CMAKE_GENERATOR",
                "CMAKE_HOME_DIRECTORY", "CMAKE_RANLIB",
                "CMAKE_STATIC_LINKER_FLAGS",
                "CMAKE_STATIC_LINKER_FLAGS_RELEASE",
                "ENABLE_OPENMP", "LEO2_BACKEND_VARIANT",
                "LEO2_BUILD_BENCHMARKS", "LEO2_BUILD_TESTS",
                "LEO2_ENABLE_CUDA", "LEOPARD_ENABLE_GF8",
                "LEOPARD_ENABLE_GF16"} and
            all(isinstance(cache[key], str) and cache[key]
                for key in ("CMAKE_AR", "CMAKE_BUILD_TYPE",
                            "CMAKE_CXX_COMPILER", "CMAKE_GENERATOR",
                            "CMAKE_HOME_DIRECTORY", "CMAKE_RANLIB")) and
            cache["CMAKE_BUILD_TYPE"] == "Release" and
            cache["CMAKE_GENERATOR"] == "Unix Makefiles" and
            Path(cache["CMAKE_HOME_DIRECTORY"]) == root and
            cache["CMAKE_EXE_LINKER_FLAGS"] == "" and
            cache["CMAKE_EXE_LINKER_FLAGS_RELEASE"] == "" and
            cache["CMAKE_STATIC_LINKER_FLAGS"] == "" and
            cache["CMAKE_STATIC_LINKER_FLAGS_RELEASE"] == "" and
            cache["CMAKE_EXPORT_COMPILE_COMMANDS"] == "ON" and
            cache["ENABLE_OPENMP"] == "ON" and
            cache["LEO2_BACKEND_VARIANT"] == "avx2" and
            cache["LEO2_BUILD_BENCHMARKS"] == "ON" and
            cache["LEO2_BUILD_TESTS"] == "ON" and
            cache["LEO2_ENABLE_CUDA"] == "OFF" and
            cache["LEOPARD_ENABLE_GF8"] == "ON" and
            cache["LEOPARD_ENABLE_GF16"] == "ON",
            label + " CMake cache semantics are invalid")
    parsed_cache = parse_cmake_cache_text(value["cache_content"]["text"])
    require({key: parsed_cache.get(key) for key in cache} == cache,
            label + " retained CMake cache values differ")
    validate_archive_members_structure(
        value["archive_members"], value["archive"],
        label + " archive members")

    closure = value["object_closure"]
    require_exact_keys(closure, {
        "schema", "compile_commands", "compile_commands_content",
        "source_root", "members", "benchmark", "sha256",
    }, label + " object closure")
    validate_file_identity_record(
        closure["compile_commands"], label + " compile commands")
    validate_exact_text_structure(
        closure["compile_commands_content"],
        "leopard2-compile-commands-content/v1",
        label + " compile commands content")
    require(closure["schema"] == "leopard2-production-object-closure/v1" and
            closure["compile_commands_content"]["file"] ==
                closure["compile_commands"] and
            closure["source_root"] == str(root) and
            Path(closure["compile_commands"]["path"]) ==
                build / "compile_commands.json" and
            isinstance(closure["members"], dict) and
            set(closure["members"]) == set(EXPECTED_ARCHIVE_MEMBERS),
            label + " object closure header is invalid")
    try:
        retained_entries = json.loads(
            closure["compile_commands_content"]["text"])
    except json.JSONDecodeError as error:
        raise EvidenceError(
            label + " retained compile_commands JSON is invalid: %s" % error)
    require(isinstance(retained_entries, list) and retained_entries,
            label + " retained compile_commands closure is empty")
    require(isinstance(value["link_recipes"], dict) and
            isinstance(value["link_recipes"].get("compiler"), dict),
            label + " compiler record is invalid")
    compiler_record = value["link_recipes"]["compiler"]
    validate_file_identity_record(compiler_record, label + " compiler")
    for member, filename, target in ARCHIVE_MEMBER_SPECS:
        record = closure["members"][member]
        require_exact_keys(record, {
            "source", "target", "compile", "archive_member",
        }, label + " " + member)
        require(record["source"] == filename and record["target"] == target and
                record["archive_member"] ==
                    value["archive_members"]["member_files"][member],
                label + " archive/object cross-binding is invalid")
        require(retained_entries.count(record["compile"]["entry"]) == 1,
                label + " production compile entry is absent or duplicated")
        validate_compile_binding_structure(
            record["compile"], build, root, filename, target,
            compiler_record, label + " " + member + " compile")
        require(record["compile"]["output"]["size"] ==
                    record["archive_member"]["size"] and
                record["compile"]["output"]["sha256"] ==
                    record["archive_member"]["sha256"],
                label + " archive member bytes differ from object bytes")
    validate_compile_binding_structure(
        closure["benchmark"], build, root, "bench/leopard2/benchmark.cpp",
        "bench_leopard2.dir", compiler_record,
        label + " benchmark compile")
    require(retained_entries.count(closure["benchmark"]["entry"]) == 1,
            label + " benchmark compile entry is absent or duplicated")
    closure_copy = dict(closure)
    closure_digest = closure_copy.pop("sha256")
    require(closure_digest == sha256_bytes(canonical_bytes(closure_copy)),
            label + " object closure digest is invalid")

    recipes = value["link_recipes"]
    require_exact_keys(recipes, {
        "schema", "archive", "executable", "archive_objects",
        "benchmark_object", "archive_tool", "ranlib_tool", "compiler",
        "sha256",
    }, label + " link recipes")
    for tool in ("archive_tool", "ranlib_tool", "compiler"):
        validate_file_identity_record(recipes[tool], label + " " + tool)
    require(recipes["schema"] == "leopard2-production-link-recipes/v1" and
            recipes["archive_tool"] == value["archive_members"]["archiver"] and
            Path(cache["CMAKE_AR"]).resolve(strict=True) ==
                Path(recipes["archive_tool"]["path"]) and
            Path(cache["CMAKE_RANLIB"]).resolve(strict=True) ==
                Path(recipes["ranlib_tool"]["path"]) and
            Path(cache["CMAKE_CXX_COMPILER"]).resolve(strict=True) ==
                Path(recipes["compiler"]["path"]),
            label + " CMake tool cross-binding is invalid")
    validate_exact_text_structure(
        recipes["archive"], "leopard2-archive-link-recipe/v1",
        label + " archive recipe")
    validate_exact_text_structure(
        recipes["executable"], "leopard2-executable-link-recipe/v1",
        label + " executable recipe")
    require(Path(recipes["archive"]["file"]["path"]) ==
                build / "CMakeFiles/leopard.dir/link.txt" and
            Path(recipes["executable"]["file"]["path"]) ==
                build / "CMakeFiles/bench_leopard2.dir/link.txt",
            label + " CMake recipe paths are invalid")
    expected_objects = [
        Path(closure["members"][member]["compile"]["output"]["path"])
            .relative_to(build).as_posix()
        for member in EXPECTED_ARCHIVE_MEMBERS
    ]
    benchmark_object = Path(
        closure["benchmark"]["output"]["path"]).relative_to(build).as_posix()
    require(recipes["archive_objects"] == expected_objects and
            recipes["benchmark_object"] == benchmark_object,
            label + " recipe object closure is invalid")
    archive_commands = [
        shlex.split(line) for line in recipes["archive"]["text"].splitlines()
        if line.strip()]
    require(len(archive_commands) == 2 and
            len(archive_commands[0]) == 3 + len(expected_objects) and
            archive_commands[0][1] in ("qc", "rc", "rcs") and
            archive_commands[0][3:] == expected_objects and
            Path(archive_commands[0][0]).resolve(strict=True) ==
                Path(recipes["archive_tool"]["path"]) and
            archive_commands[0][2] == Path(value["archive"]["path"]).name and
            len(archive_commands[1]) == 2 and
            Path(archive_commands[1][0]).resolve(strict=True) ==
                Path(recipes["ranlib_tool"]["path"]) and
            archive_commands[1][1] == Path(value["archive"]["path"]).name,
            label + " retained archive recipe semantics are invalid")
    executable_lines = [line for line in recipes["executable"]["text"].splitlines()
                        if line.strip()]
    require(len(executable_lines) == 1,
            label + " retained executable recipe line count is invalid")
    link_arguments = shlex.split(executable_lines[0])
    require(link_arguments and
            Path(link_arguments[0]).resolve(strict=True) ==
                Path(recipes["compiler"]["path"]) and
            link_arguments.count(benchmark_object) == 1 and
            link_arguments.count(Path(value["archive"]["path"]).name) == 1 and
            link_arguments.index(benchmark_object) <
                link_arguments.index(Path(value["archive"]["path"]).name) and
            link_arguments.count("-o") == 1 and
            link_arguments.index("-o") + 1 < len(link_arguments) and
            link_arguments[link_arguments.index("-o") + 1] ==
                Path(value["executable"]["path"]).name and
            all("@" not in argument for argument in link_arguments) and
            [argument for argument in link_arguments
             if argument.endswith(".o")] == [benchmark_object],
            label + " retained executable recipe semantics are invalid")
    recipes_copy = dict(recipes)
    recipes_digest = recipes_copy.pop("sha256")
    require(recipes_digest == sha256_bytes(canonical_bytes(recipes_copy)),
            label + " link recipe digest is invalid")
    copy = dict(value)
    digest = copy.pop("sha256")
    require(digest == sha256_bytes(canonical_bytes(copy)),
            label + " digest is invalid")


def validate_binary_identity_structure(value: Any) -> None:
    require(isinstance(value, dict), "binary identity is not an object")
    modes = value.get("comparison_modes")
    require(isinstance(modes, dict) and
            set(modes) == {"baseline", "candidate"},
            "binary identity comparison modes are invalid")
    comparison_modes(modes["baseline"], modes["candidate"])
    expected_keys = set(BINARY_FILE_IDENTITY_KEYS) | {
        "comparison_modes",
        "baseline_avx2_compile_entry", "baseline_avx2_xor_compile_entry",
        "baseline_core_compile_entry",
        "baseline_avx512_compile_entry", "baseline_build_id",
        "candidate_avx2_compile_entry", "candidate_avx2_xor_compile_entry",
        "candidate_core_compile_entry",
        "candidate_avx512_compile_entry", "candidate_build_id", "source",
        "baseline_archive_members", "candidate_archive_members",
        "baseline_build_provenance", "candidate_build_provenance",
    }
    if is_xmm_tail_comparison(modes):
        expected_keys.update({
            "baseline_avx2_encoding_audit",
            "candidate_avx2_encoding_audit",
            "xmm_tail_comparison_audit",
        })
    require_exact_keys(value, expected_keys, "binary identity")
    require(isinstance(value.get("source"), dict) and
            isinstance(value["source"].get("root"), str),
            "binary identity source is invalid")
    for key in BINARY_FILE_IDENTITY_KEYS:
        validate_file_identity_record(value[key], key)
    for label in ("baseline", "candidate"):
        require(isinstance(value[label + "_build_id"], str) and
                value[label + "_build_id"],
                "%s build ID is invalid" % label)
        for component in ("avx2", "avx2_xor", "core", "avx512"):
            require(isinstance(
                        value[label + "_" + component + "_compile_entry"],
                        dict),
                    "%s %s compile entry is invalid" % (label, component))
        archive_members = value[label + "_archive_members"]
        validate_archive_members_structure(
            archive_members, value[label + "_archive"],
            "%s archive membership" % label)
        provenance = value[label + "_build_provenance"]
        require(isinstance(provenance, dict) and
                isinstance(provenance.get("object_closure"), dict),
                "%s build provenance is invalid" % label)
        require(provenance.get("executable") ==
                    value[label + "_executable"] and
                provenance.get("archive") == value[label + "_archive"] and
                provenance.get("archive_members") == archive_members and
                provenance.get("object_closure", {}).get(
                    "compile_commands") ==
                    value[label + "_compile_commands"],
                "%s build provenance is not cross-bound" % label)
        validate_build_provenance_structure(
            provenance, value["source"], "%s build provenance" % label)
        closure_members = provenance["object_closure"]["members"]
        explicit_members = {
            "avx2": "Leopard2BackendAVX2.cpp.o",
            "avx2_xor": "Leopard2BackendAVX2Xor.cpp.o",
            "avx512": "Leopard2BackendAVX512.cpp.o",
            "core": "leopard2.cpp.o",
        }
        for component, member in explicit_members.items():
            binding = closure_members[member]["compile"]
            require(value[label + "_" + component + "_object"] ==
                        binding["output"] and
                    value[label + "_" + component + "_compile_entry"] ==
                        binding["entry"],
                    "%s %s record is not bound to production closure" %
                        (label, component))
        if is_xmm_tail_comparison(modes):
            validate_avx2_encoding_audit(
                value[label + "_avx2_encoding_audit"],
                value[label + "_avx2_object"],
                label + " AVX2 encoding audit")
    if is_xmm_tail_comparison(modes):
        validate_xmm_tail_audit_pair(
            value["baseline_avx2_encoding_audit"],
            value["candidate_avx2_encoding_audit"])
        validate_xmm_tail_comparison_audit(
            value["xmm_tail_comparison_audit"],
            value["baseline_avx2_encoding_audit"],
            value["candidate_avx2_encoding_audit"])
    require(value["baseline_executable"]["sha256"] !=
            value["candidate_executable"]["sha256"] and
            value["baseline_archive"]["sha256"] !=
            value["candidate_archive"]["sha256"] and
            value["baseline_avx2_xor_object"]["sha256"] ==
            value["candidate_avx2_xor_object"]["sha256"] and
            value["baseline_avx512_object"]["sha256"] ==
            value["candidate_avx512_object"]["sha256"],
            "binary comparison identity is invalid")
    baseline_closure = value["baseline_build_provenance"]["object_closure"]
    candidate_closure = value["candidate_build_provenance"]["object_closure"]
    validate_comparison_object_deltas(
        modes,
        {member: baseline_closure["members"][member]["compile"]["output"]
            ["sha256"] for member in EXPECTED_ARCHIVE_MEMBERS},
        {member: candidate_closure["members"][member]["compile"]["output"]
            ["sha256"] for member in EXPECTED_ARCHIVE_MEMBERS})
    require(baseline_closure["benchmark"]["output"]["sha256"] ==
            candidate_closure["benchmark"]["output"]["sha256"],
            "benchmark object differs between diagnostic builds")
    baseline_cache = value["baseline_build_provenance"]["cache_values"]
    candidate_cache = value["candidate_build_provenance"]["cache_values"]
    require(shlex.split(baseline_cache["CMAKE_CXX_FLAGS"]) ==
                mode_compile_arguments(modes["baseline"]) and
            shlex.split(candidate_cache["CMAKE_CXX_FLAGS"]) ==
                mode_compile_arguments(modes["candidate"]) and
            baseline_cache["CMAKE_CXX_FLAGS_RELEASE"] ==
                candidate_cache["CMAKE_CXX_FLAGS_RELEASE"] ==
                "-O3 -DNDEBUG",
            "retained diagnostic build flags are invalid")
    for key in baseline_cache:
        if key != "CMAKE_CXX_FLAGS":
            require(baseline_cache[key] == candidate_cache[key],
                    "retained diagnostic CMake cache differs at %s" % key)
    require(value["baseline_build_provenance"]["link_recipes"]["archive"]
                ["text"] ==
            value["candidate_build_provenance"]["link_recipes"]["archive"]
                ["text"],
            "retained diagnostic archive link recipe differs")
    baseline_link = shlex.split(
        value["baseline_build_provenance"]["link_recipes"]["executable"]
            ["text"])
    candidate_link = shlex.split(
        value["candidate_build_provenance"]["link_recipes"]["executable"]
            ["text"])
    baseline_link = strip_mode_definition(
        baseline_link, modes["baseline"], "retained baseline link")
    candidate_link = strip_mode_definition(
        candidate_link, modes["candidate"], "retained candidate link")
    require(baseline_link == candidate_link,
            "retained executable link differs beyond the modes")
    for tool in ("archive_tool", "ranlib_tool", "compiler"):
        require(value["baseline_build_provenance"]["link_recipes"][tool]
                    ["sha256"] ==
                value["candidate_build_provenance"]["link_recipes"][tool]
                    ["sha256"],
                "retained diagnostic tool differs at %s" % tool)
    for member in EXPECTED_ARCHIVE_MEMBERS:
        baseline_arguments = list(
            baseline_closure["members"][member]["compile"]["arguments"])
        candidate_arguments = list(
            candidate_closure["members"][member]["compile"]["arguments"])
        baseline_arguments = strip_mode_definition(
            baseline_arguments, modes["baseline"],
            "retained baseline compile %s" % member)
        candidate_arguments = strip_mode_definition(
            candidate_arguments, modes["candidate"],
            "retained candidate compile %s" % member)
        require(baseline_arguments == candidate_arguments,
                "retained compile commands differ beyond modes: %s" %
                    member)
    baseline_arguments = list(
        baseline_closure["benchmark"]["arguments"])
    candidate_arguments = list(
        candidate_closure["benchmark"]["arguments"])
    baseline_arguments = strip_mode_definition(
        baseline_arguments, modes["baseline"],
        "retained baseline benchmark compile")
    candidate_arguments = strip_mode_definition(
        candidate_arguments, modes["candidate"],
        "retained candidate benchmark compile")
    require(baseline_arguments == candidate_arguments,
            "retained benchmark commands differ beyond modes")
    source = value["source"]
    require_exact_keys(source, {
        "root", "head", "head_tree", "status_short", "files",
        "snapshot_sha256",
    }, "source identity")
    require(isinstance(source["root"], str) and
            Path(source["root"]).is_absolute() and
            isinstance(source["head"], str) and source["head"] and
            isinstance(source["head_tree"], str) and source["head_tree"] and
            isinstance(source["status_short"], list) and
            all(isinstance(line, str) for line in source["status_short"]),
            "source identity header is invalid")
    files = source["files"]
    require(isinstance(files, dict) and set(files) == set(SOURCE_FILES),
            "source identity file set is invalid")
    for relative, record in files.items():
        require_exact_keys(record, {"size", "sha256"},
                           "source file %s" % relative)
        require(type(record["size"]) is int and record["size"] >= 0 and
                isinstance(record["sha256"], str) and
                re.fullmatch(r"[0-9a-f]{64}", record["sha256"]) is not None,
                "source file identity is invalid")
    source_copy = dict(source)
    snapshot_digest = source_copy.pop("snapshot_sha256")
    require(snapshot_digest == sha256_bytes(canonical_bytes(source_copy)),
            "source snapshot digest is invalid")


def verify_campaign(options: argparse.Namespace) -> int:
    manifest_path = options.manifest.resolve(strict=True)
    campaign_root = manifest_path.parent
    require(manifest_path == campaign_root / "manifest.json",
            "campaign manifest must use its canonical filename")
    manifest = json.loads(manifest_path.read_text())
    require_exact_keys(manifest, {
        "schema", "created_utc", "request", "global_lock", "pair_lease",
        "affinity_exclusion", "isolation_epoch", "host", "statistics",
        "cell_summaries", "digest",
    }, "campaign manifest")
    require(manifest["schema"] in (SCHEMA, XMM_TAIL_SCHEMA),
            "wrong manifest schema")
    manifest_copy = dict(manifest)
    expected_digest = manifest_copy.pop("digest")
    require(expected_digest == sha256_bytes(canonical_bytes(manifest_copy)),
            "manifest digest mismatch")
    require(isinstance(manifest["created_utc"], str),
            "manifest timestamp is invalid")
    timestamp = dt.datetime.fromisoformat(manifest["created_utc"])
    require(timestamp.tzinfo is not None,
            "manifest timestamp is not timezone-aware")

    request = manifest["request"]
    require_exact_keys(request, {
        "schema", "matrix_sha256", "matrix_preset", "selected_cell_ids",
        "rounds", "order", "binary_identity", "cpu", "reserved_sibling",
        "reservation", "topology_identity", "isolation_policy",
        "timeout_seconds", "max_retries", "execution_nonce",
        "resume_policy", "lock", "request_digest",
    }, "campaign request")
    request_copy = dict(request)
    request_digest = request_copy.pop("request_digest")
    require(request_digest == sha256_bytes(canonical_bytes(request_copy)),
            "request digest mismatch")
    binary_identity = request.get("binary_identity")
    modes = (binary_identity.get("comparison_modes")
             if isinstance(binary_identity, dict) else None)
    require(isinstance(modes, dict) and
            set(modes) == {"baseline", "candidate"},
            "request comparison modes are invalid")
    comparison_modes(modes["baseline"], modes["candidate"])
    require_comparison_schema(
        request["schema"], modes, "request")
    require_comparison_schema(
        manifest["schema"], modes, "manifest")
    require(request["resume_policy"] ==
                "fresh complete authoritative campaign only" and
            request["isolation_policy"] == isolation_policy() and
            request["order"] == list(ORDER) and
            request["rounds"] in (1, 3, 5) and
            type(request["timeout_seconds"]) in (int, float) and
            math.isfinite(request["timeout_seconds"]) and
            request["timeout_seconds"] > 0 and
            type(request["max_retries"]) is int and
            0 <= request["max_retries"] <= 99 and
            isinstance(request["execution_nonce"], str) and
            re.fullmatch(r"[0-9a-f]{32}",
                         request["execution_nonce"]) is not None,
            "campaign request policy is invalid")
    cpu = require_int(request["cpu"], "benchmark CPU")
    sibling = require_int(request["reserved_sibling"], "reserved sibling")
    require(cpu != sibling, "campaign CPU pair is invalid")

    matrix = matrix_for_preset(request["matrix_preset"])
    require(request["matrix_sha256"] == matrix["matrix_sha256"],
            "campaign does not use the frozen matrix")
    matrix_path = campaign_root / "matrix.json"
    require(matrix_path.is_file() and json.loads(matrix_path.read_text()) ==
            matrix, "retained frozen matrix changed")
    selected_ids = request["selected_cell_ids"]
    require(isinstance(selected_ids, list) and selected_ids and
            all(isinstance(identifier, str) for identifier in selected_ids) and
            len(selected_ids) == len(set(selected_ids)),
            "selected cell IDs are invalid")
    selected = [cell for cell in matrix["cells"]
                if cell["id"] in set(selected_ids)]
    require([cell["id"] for cell in selected] == selected_ids,
            "selected cells are unknown or not in frozen order")
    validate_comparison_matrix_selection(
        request["matrix_preset"], modes, matrix, selected)
    request_path = campaign_root / "request.json"
    require(request_path.is_file() and
            json.loads(request_path.read_text()) == request,
            "retained request changed")

    validate_binary_identity_structure(request["binary_identity"])
    validate_global_lock_record(request["lock"])
    validate_topology_identity(request["topology_identity"], cpu, sibling)
    validate_reservation_record(request["reservation"], cpu, sibling)
    validate_pair_lease_record(manifest["pair_lease"], cpu, sibling)
    validate_global_lock_record(manifest["global_lock"])

    epoch = manifest["isolation_epoch"]
    require_exact_keys(epoch, {"start", "start_digest", "restoration"},
                       "isolation epoch")
    start = epoch["start"]
    require_exact_keys(start, {
        "schema", "execution_nonce", "global_lock", "pair_lease",
        "reservation", "topology_identity", "affinity_exclusion",
    }, "isolation epoch start")
    require(start["schema"] == "leopard2-direct-isolation-epoch/v4" and
            epoch["start_digest"] == sha256_bytes(canonical_bytes(start)) and
            start["execution_nonce"] == request["execution_nonce"] and
            start["global_lock"] == request["lock"] ==
                manifest["global_lock"] and
            start["pair_lease"] == manifest["pair_lease"] and
            start["reservation"] == request["reservation"] and
            start["topology_identity"] == request["topology_identity"] and
            start["affinity_exclusion"] == manifest["affinity_exclusion"],
            "isolation epoch records are not cross-bound")
    validate_affinity_exclusion(
        manifest["affinity_exclusion"], epoch["restoration"], cpu, sibling)
    lease_uid = manifest["pair_lease"]["payload"]["uid"]
    require(lease_uid == manifest["affinity_exclusion"]["uid"] ==
            manifest["global_lock"]["uid"],
            "isolation records belong to different users")
    runner_before = manifest["affinity_exclusion"]["runner_before"]
    runner_after = manifest["affinity_exclusion"]["runner_after"]
    require(runner_before["affinity"] ==
                request["topology_identity"]["launch_affinity"] and
            runner_after["affinity"] ==
                request["topology_identity"]["housekeeping_affinity"],
            "runner affinity is not bound to topology identity")

    host = manifest["host"]
    require_exact_keys(host, {
        "node", "platform", "machine", "allowed_cpus",
    }, "host identity")
    require(all(isinstance(host[key], str)
                for key in ("node", "platform", "machine")) and
            host["allowed_cpus"] == runner_before["affinity"],
            "host/runner CPU identity is inconsistent")
    expected_statistics = {
        "method": "one mean log contrast per independent ABBA round",
        "confidence": 0.95,
        "rounds": request["rounds"],
        "degrees_of_freedom": request["rounds"] - 1,
        "directional_only": request["rounds"] == 1,
        "child_estimator": "median of retained per-invocation samples",
    }
    require(manifest["statistics"] == expected_statistics,
            "campaign statistics declaration is invalid")

    if not options.no_current_input_check:
        validate_current_binary_identity(request["binary_identity"])
        require(reservation_identity(
                    Path(request["reservation"]["file"]["path"]),
                    cpu, sibling) == request["reservation"],
                "current CPU reservation changed")
        linked_lock = DEFAULT_LOCK.lstat()
        require(linked_lock.st_dev == request["lock"]["device"] and
                linked_lock.st_ino == request["lock"]["inode"] and
                linked_lock.st_uid == request["lock"]["uid"] and
                linked_lock.st_nlink == request["lock"]["nlink"] and
                stat.S_IMODE(linked_lock.st_mode) == request["lock"]["mode"],
                "current canonical GF8 lock identity changed")

    cells = manifest["cell_summaries"]
    require(isinstance(cells, list) and len(cells) == len(selected),
            "manifest cell set is incomplete")
    seen_envelopes = set()
    for summary, expected_cell in zip(cells, selected):
        cell_path = campaign_root / "cells" / (expected_cell["id"] + ".json")
        require(cell_path.is_file() and
                json.loads(cell_path.read_text()) == summary,
                "retained cell summary changed")
        require(isinstance(summary, dict) and
                summary.get("cell") == expected_cell and
                summary.get("binary_identity") == request["binary_identity"],
                "cell summary is not bound to request identity")
        invocation_records = summary.get("invocations")
        expected_count = request["rounds"] * len(ORDER)
        require(isinstance(invocation_records, list) and
                len(invocation_records) == expected_count,
                "cell summary has the wrong invocation count")
        retained = []
        for invocation_index, record in enumerate(invocation_records):
            require_exact_keys(record, {
                "implementation", "median_us", "plan_setup_us",
                "envelope_path", "envelope_sha256",
                "reserved_sibling_nonidle_jiffies",
                "same_user_pair_affinity_before",
                "same_user_pair_affinity_after", "target_runtime",
                "target_interrupts", "isolation_epoch_digest",
            }, "cell invocation record")
            round_index, slot_index = divmod(invocation_index, len(ORDER))
            implementation = ORDER[slot_index]
            envelope_path = Path(record["envelope_path"])
            require(envelope_path not in seen_envelopes and
                    envelope_path.is_file() and
                    sha256_file(envelope_path) == record["envelope_sha256"],
                    "retained invocation envelope identity changed")
            seen_envelopes.add(envelope_path)
            value = json.loads(envelope_path.read_text())
            retained.append(validate_envelope(
                value, envelope_path, expected_cell, implementation,
                round_index, slot_index, request, epoch["start_digest"],
                campaign_root))
        reconstructed = summarize_cell(
            expected_cell, retained, request["rounds"],
            request["binary_identity"])
        require(summary == reconstructed,
                "cell summary does not reconstruct from retained envelopes")
    print("verified %d uncontaminated ABBA cells" % len(cells))
    return 0


def available_physical_cpus() -> list[int]:
    allowed = sorted(os.sched_getaffinity(0))
    physical = []
    seen = set()
    for cpu in allowed:
        topology = Path("/sys/devices/system/cpu/cpu%d/topology" % cpu)
        try:
            package = int((topology / "physical_package_id").read_text())
            core = int((topology / "core_id").read_text())
            key = ("core", package, core)
        except (OSError, ValueError):
            siblings_path = topology / "thread_siblings_list"
            try:
                siblings = tuple(sorted(
                    parse_cpu_list(siblings_path.read_text()) & set(allowed)))
            except (OSError, ValueError):
                siblings = (cpu,)
            key = ("siblings", *siblings)
        if key in seen:
            continue
        seen.add(key)
        physical.append(cpu)
    require(physical, "no process-visible physical CPU is available")
    return physical


def parse_variants(values: Sequence[str]) -> dict[str, Path]:
    variants: dict[str, Path] = {}
    for value in values:
        label, separator, raw_path = value.partition("=")
        require(separator == "=" and re.fullmatch(r"[a-z][a-z0-9_-]*", label),
                "--variant must be LABEL=EXECUTABLE with a safe label")
        require(label not in variants, "duplicate variant label %s" % label)
        path = Path(raw_path).resolve(strict=True)
        require(path.is_file() and os.access(path, os.X_OK),
                "variant %s is not executable" % label)
        variants[label] = path
    required = {
        "full", "tile4096", "tile8192", "tile16384", "tile32768",
        "tile65536",
    }
    require(set(variants) == required,
            "tile screen requires exactly: %s" % ", ".join(sorted(required)))
    identities = {label: file_identity(path) for label, path in variants.items()}
    require(len({value["sha256"] for value in identities.values()}) ==
            len(identities), "tile-screen executables must all differ")
    return variants


def write_tile_screen_spec(options: argparse.Namespace) -> int:
    require(options.timeout > 0, "tile-screen timeout must be positive")
    variants = parse_variants(options.variant)
    source_root = options.source_root.resolve(strict=True)
    variant_builds = {
        label: build_variant_identity(executable)
        for label, executable in variants.items()
    }
    expected_avx2_source = (source_root / "Leopard2BackendAVX2.cpp").resolve()
    expected_avx2_xor_source = (
        source_root / "Leopard2BackendAVX2Xor.cpp").resolve()
    expected_avx512_source = (
        source_root / "Leopard2BackendAVX512.cpp").resolve()
    expected_core_source = (source_root / "leopard2.cpp").resolve()
    for label, identity in variant_builds.items():
        require(compile_entry_source(identity["avx2_compile_entry"]) ==
                expected_avx2_source,
                "tile variant %s uses a different AVX2 source" % label)
        require(compile_entry_source(identity["avx2_xor_compile_entry"]) ==
                expected_avx2_xor_source,
                "tile variant %s uses a different AVX2 Xor source" % label)
        require(compile_entry_source(identity["core_compile_entry"]) ==
                expected_core_source,
                "tile variant %s uses a different core source" % label)
        require(compile_entry_source(identity["avx512_compile_entry"]) ==
                expected_avx512_source,
                "tile variant %s uses a different AVX-512 source" % label)
    require(len({
                identity["avx2_xor_object"]["sha256"]
                for identity in variant_builds.values()}) == 1,
            "payload-tiling diagnostics unexpectedly changed AVX2 Xor code")
    require(len({
                identity["avx512_object"]["sha256"]
                for identity in variant_builds.values()}) == 1,
            "payload-tiling diagnostics unexpectedly changed AVX-512 code")
    output = options.output.resolve()
    try:
        output.relative_to(source_root)
    except ValueError:
        pass
    else:
        raise EvidenceError(
            "tile-screen spec must be written outside the source worktree")
    physical_cpus = available_physical_cpus()
    matrix = make_large_matrix()
    cells = [cell for cell in matrix["cells"]
             if cell["requested_tile_screen_loss"]]
    jobs = []
    tile_labels = (
        "tile4096", "tile8192", "tile16384", "tile32768", "tile65536",
    )
    variant_labels = ("full", *tile_labels)
    for screen_index, cell in enumerate(cells):
        group = "r%03d-b%09d-l%02d" % (
            cell["R"], cell["bytes"], cell["loss"])
        cpu = physical_cpus[screen_index % len(physical_cpus)]
        minimum_memory_mb = (
            cell["estimated_peak_bytes"] + (1 << 20) - 1) >> 20
        memory_limit_mb = minimum_memory_mb + 256
        rotation = screen_index % len(tile_labels)
        rotated_tiles = tile_labels[rotation:] + tile_labels[:rotation]
        variant_order = (("full", *rotated_tiles) if screen_index % 2 == 0
                         else (*rotated_tiles, "full"))
        for variant_index, label in enumerate(variant_order):
            benchmark_cell = dict(cell)
            benchmark_cell["variant"] = label
            # The retained tile-screen utility predates this experiment and
            # is not part of its promotion evidence.  If it is used for an
            # exploratory source-major build, retain the execution-mode
            # identity in each job so result validation is still explicit.
            benchmark_cell["mode"] = "source"
            benchmark_cell["directional_order_index"] = variant_index
            jobs.append({
                "id": "tile.%s.v%02d-%s" % (
                    group, variant_index, label),
                "command": [str(variants[label]),
                            *benchmark_arguments(cell)],
                "cwd": str(source_root),
                "cpu_set": [cpu],
                "cpu_group": group,
                "resume_group": group,
                "timeout_seconds": options.timeout,
                "memory_mb": memory_limit_mb,
                "rss_limit_mb": memory_limit_mb,
                "minimum_memory_mb": minimum_memory_mb,
                "benchmark_cell": benchmark_cell,
            })
    spec = {
        "schema": TILE_SPEC_SCHEMA,
        "root": str(source_root),
        "workers": len(physical_cpus),
        "base_seed": 0xD1EC7A11,
        "metadata": {
            "campaign": "non-authoritative GF8 direct payload-tile screen",
            "directional_only": True,
            "matrix_sha256": matrix["matrix_sha256"],
            "physical_cpus": physical_cpus,
            "variant_labels": list(variant_labels),
            "ordering": (
                "Tile order rotates by cell; full alternates first/last to "
                "reduce monotonic drift in this directional-only screen."),
            "variant_build_identities": {
                label: variant_builds[label] for label in variant_labels
            },
            "runner": file_identity(Path(__file__)),
            "source": source_identity(source_root),
            "follow_up": (
                "Use serial sibling-idle ABBA only for per-loss/size winners "
                "and adjacent unchanged controls."),
        },
        "defaults": {
            "cpu_count": 1,
            "cpu_policy": "physical-first",
            "memory_mb": 0,
            "rss_limit_mb": 0,
            "minimum_memory_mb": 0,
            "timeout_seconds": options.timeout,
        },
        "jobs": jobs,
    }
    atomic_json(output, spec)
    print("wrote %d directional jobs over %d cells and %d physical CPUs" % (
        len(jobs), len(cells), len(physical_cpus)))
    return 0


def run_tile_screen_lab(options: argparse.Namespace) -> int:
    require(options.lock_timeout >= 0 and options.progress_seconds >= 0,
            "tile-screen lock timeout and progress interval must be nonnegative")
    manifest = options.manifest.resolve(strict=True)
    value = json.loads(manifest.read_text())
    require(value.get("source_spec", {}).get("schema") == TILE_SPEC_SCHEMA,
            "manifest was not generated from the frozen tile-screen spec")
    lab = options.lab.resolve(strict=True)
    require(lab.is_file(), "leopard2 lab runner is missing")
    command = [
        sys.executable, str(lab), "run", "--manifest", str(manifest),
        "--output-dir", str(options.output_dir.resolve()),
        "--progress-seconds", str(options.progress_seconds),
    ]
    if options.rerun_failed:
        command.append("--rerun-failed")
    lock_stream = acquire_lock(options.lock, options.lock_timeout)
    try:
        completed = subprocess.run(command, check=False)
    finally:
        fcntl.flock(lock_stream.fileno(), fcntl.LOCK_UN)
        lock_stream.close()
    require(completed.returncode == 0,
            "leopard2 lab screen failed with exit code %d" %
            completed.returncode)
    return 0


def analyze_tile_screen(options: argparse.Namespace) -> int:
    manifest_path = options.manifest.resolve(strict=True)
    manifest = json.loads(manifest_path.read_text())
    manifest_unsigned = dict(manifest)
    manifest_digest = manifest_unsigned.pop("manifest_digest", None)
    require(manifest_digest == sha256_bytes(canonical_bytes(manifest_unsigned)),
            "lab manifest digest is invalid")
    require(manifest.get("source_spec", {}).get("schema") == TILE_SPEC_SCHEMA,
            "manifest was not generated from the frozen tile-screen spec")
    output_dir = options.output_dir.resolve(strict=True)
    by_cell: dict[tuple[int, int, int, int], dict[str, Any]] = {}
    variant_set = set(manifest["source_spec"]["metadata"]["variant_labels"])
    for job in manifest.get("jobs", []):
        benchmark_cell = job.get("benchmark_cell")
        require(isinstance(benchmark_cell, dict),
                "manifest job lacks benchmark-cell identity")
        variant = benchmark_cell.get("variant")
        require(variant in variant_set, "manifest job has unknown variant")
        result_path = output_dir / "jobs" / job["id"] / "result.json"
        result = json.loads(result_path.read_text())
        result_unsigned = dict(result)
        result_digest = result_unsigned.pop("result_digest", None)
        require(result_digest == sha256_bytes(canonical_bytes(result_unsigned)),
                "lab job %s result digest is invalid" % job["id"])
        require(result.get("schema") == LAB_RESULT_SCHEMA and
                result.get("outcome") == "success" and
                result.get("job_digest") == job["job_digest"],
                "lab job %s did not complete successfully" % job["id"])
        stdout_path = result_path.parent / "stdout.txt"
        stdout = stdout_path.read_bytes()
        identity = result.get("outputs", {}).get("stdout", {})
        require(identity.get("size_bytes") == len(stdout) and
                identity.get("sha256") == sha256_bytes(stdout),
                "lab job %s stdout identity changed" % job["id"])
        try:
            parsed = json.loads(stdout)
        except json.JSONDecodeError as error:
            raise EvidenceError(
                "lab job %s emitted invalid JSON: %s" % (job["id"], error))
        normalized = validate_result(
            parsed, benchmark_cell, benchmark_cell["mode"])
        key = (benchmark_cell["K"], benchmark_cell["R"],
               benchmark_cell["bytes"], benchmark_cell["loss"])
        cell = by_cell.setdefault(key, {
            "K": key[0], "R": key[1], "bytes": key[2], "loss": key[3],
            "variants": {},
        })
        require(variant not in cell["variants"],
                "duplicate variant in a tile-screen cell")
        cell["variants"][variant] = normalized

    cells = []
    aggregate_logs: dict[tuple[int, int, str], list[float]] = {}
    for key in sorted(by_cell):
        cell = by_cell[key]
        require(set(cell["variants"]) == variant_set,
                "tile-screen cell is incomplete")
        reference = cell["variants"]["full"]
        for value in cell["variants"].values():
            require(value["digests"] == reference["digests"] and
                    value["missing_original_indices"] ==
                    reference["missing_original_indices"],
                    "tile-screen variants changed wire data or erasures")
        ratios = {}
        for variant in sorted(variant_set - {"full"}):
            value = cell["variants"][variant]
            metrics = {
                "execution": reference["median_us"] / value["median_us"],
            }
            for reuse in (1, 8, 64):
                metrics["amortized_reuse_%d" % reuse] = (
                    reference["median_us"] +
                    reference["plan_setup_us"] / reuse) / (
                    value["median_us"] + value["plan_setup_us"] / reuse)
            ratios[variant] = metrics
            aggregate_logs.setdefault(
                (cell["bytes"], cell["loss"], variant), []).append(
                    math.log(metrics["execution"]))
        winner = max(ratios, key=lambda label: ratios[label]["execution"])
        cells.append({
            "K": cell["K"], "R": cell["R"], "bytes": cell["bytes"],
            "loss": cell["loss"], "orientation": "full_over_tiled",
            "ratios": ratios, "directional_winner": winner,
            "directional_winner_ratio": ratios[winner]["execution"],
            "digests": reference["digests"],
            "missing_original_indices": reference["missing_original_indices"],
        })
    aggregates = []
    for (byte_count, loss_count, variant), logs in sorted(
            aggregate_logs.items()):
        require(len(logs) == 3,
                "tile-screen aggregate does not contain all R values")
        aggregates.append({
            "bytes": byte_count, "loss": loss_count, "variant": variant,
            "R_values": [65, 96, 128],
            "geometric_mean_directional_ratio": math.exp(
                statistics.mean(logs)),
        })
    result = {
        "schema": TILE_SCREEN_SCHEMA,
        "directional_only": True,
        "warning": (
            "Cells ran concurrently across physical cores and are only a "
            "shortlist. Promotion requires serial sibling-idle ABBA."),
        "manifest": file_identity(manifest_path),
        "cell_count": len(cells),
        "cells": cells,
        "aggregates": aggregates,
    }
    result["digest"] = sha256_bytes(canonical_bytes(result))
    if options.output == Path("-"):
        print(json.dumps(result, sort_keys=True, indent=2))
    else:
        atomic_json(options.output, result)
        print("wrote %d directional tile-screen cells" % len(cells))
    return 0


def write_matrix(options: argparse.Namespace) -> int:
    value = matrix_for_preset(options.preset)
    if options.output == Path("-"):
        print(json.dumps(value, sort_keys=True, indent=2))
    else:
        atomic_json(options.output, value)
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    matrix = commands.add_parser("matrix", help="write a frozen matrix")
    matrix.add_argument("--output", type=Path, default=Path("-"))
    matrix.add_argument(
        "--preset", choices=("core", "large", "tiny", "xmm-tail"),
        default="core")
    matrix.set_defaults(function=write_matrix)
    run = commands.add_parser(
        "run", help="run a fresh, complete same-source campaign")
    run.add_argument("--baseline", required=True, type=Path)
    run.add_argument("--candidate", required=True, type=Path)
    run.add_argument("--baseline-archive", required=True, type=Path)
    run.add_argument("--candidate-archive", required=True, type=Path)
    run.add_argument("--baseline-avx2-object", required=True, type=Path)
    run.add_argument("--candidate-avx2-object", required=True, type=Path)
    run.add_argument("--baseline-avx2-xor-object", required=True, type=Path)
    run.add_argument("--candidate-avx2-xor-object", required=True, type=Path)
    run.add_argument("--baseline-compile-commands", required=True, type=Path)
    run.add_argument("--candidate-compile-commands", required=True, type=Path)
    run.add_argument("--source-root", required=True, type=Path)
    run.add_argument("--baseline-build-id", required=True)
    run.add_argument("--candidate-build-id", required=True)
    run.add_argument(
        "--baseline-mode", required=True,
        choices=tuple(MODE_COMPILE_DEFINITIONS))
    run.add_argument(
        "--candidate-mode", required=True,
        choices=tuple(MODE_COMPILE_DEFINITIONS))
    run.add_argument("--output", required=True, type=Path)
    run.add_argument(
        "--preset", choices=("core", "large", "tiny", "xmm-tail"),
        default="core")
    run.add_argument("--cpu", required=True, type=int)
    run.add_argument("--reserved-sibling", required=True, type=int)
    run.add_argument("--reservation-file", required=True, type=Path)
    run.add_argument("--cell", action="append",
                     help="run only this frozen matrix cell ID; repeat as needed")
    run.add_argument("--rounds", type=int, default=DEFAULT_ROUNDS)
    run.add_argument("--timeout", type=float, default=120.0)
    run.add_argument("--max-retries", type=int, default=10)
    run.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    run.add_argument("--lock-timeout", type=float, default=3600.0)
    run.set_defaults(function=run_campaign)
    verify = commands.add_parser("verify", help="verify a retained campaign")
    verify.add_argument("--manifest", required=True, type=Path)
    verify.add_argument("--no-current-input-check", action="store_true",
                        help="verify retained evidence after build paths moved")
    verify.set_defaults(function=verify_campaign)
    tile_spec = commands.add_parser(
        "tile-spec", help="write a 15-core memory-budgeted tile-screen spec")
    tile_spec.add_argument(
        "--variant", action="append", required=True,
        help=("one LABEL=EXECUTABLE; require full and tile4096 through "
              "tile65536"))
    tile_spec.add_argument("--source-root", required=True, type=Path)
    tile_spec.add_argument("--output", required=True, type=Path)
    tile_spec.add_argument("--timeout", type=float, default=600.0)
    tile_spec.set_defaults(function=write_tile_screen_spec)
    tile_run = commands.add_parser(
        "tile-run", help="run a tile-screen lab manifest under the GF8 lock")
    tile_run.add_argument("--manifest", required=True, type=Path)
    tile_run.add_argument("--output-dir", required=True, type=Path)
    tile_run.add_argument("--lab", required=True, type=Path)
    tile_run.add_argument("--rerun-failed", action="store_true")
    tile_run.add_argument("--progress-seconds", type=float, default=1.0)
    tile_run.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    tile_run.add_argument("--lock-timeout", type=float, default=3600.0)
    tile_run.set_defaults(function=run_tile_screen_lab)
    tile_analyze = commands.add_parser(
        "tile-analyze", help="analyze a completed directional tile screen")
    tile_analyze.add_argument("--manifest", required=True, type=Path)
    tile_analyze.add_argument("--output-dir", required=True, type=Path)
    tile_analyze.add_argument("--output", type=Path, default=Path("-"))
    tile_analyze.set_defaults(function=analyze_tile_screen)
    return result


def main(arguments: Sequence[str] | None = None) -> int:
    options = parser().parse_args(arguments)
    try:
        return options.function(options)
    except (EvidenceError, IndexError, KeyError, OSError,
            subprocess.SubprocessError, TypeError, ValueError) as error:
        print("error: %s" % error, file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
