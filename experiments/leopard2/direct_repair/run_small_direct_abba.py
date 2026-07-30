#!/usr/bin/env python3
"""Run a bounded, fail-closed ABBA screen for small GF8 direct repair.

This runner compares two already-built Leopard2 benchmark executables whose
only semantic difference is an explicitly selected small-code direct-repair
mode: transform, output-major, or source-major.  It is not for comparison to
the historical Leopard executable (use ../main_compare/run_abba.py for that).
It owns the canonical GF8 and physical-pair leases for the complete campaign,
moves all same-user threads off the held pair, pins every child to one logical
CPU, requires its SMT sibling to remain idle, retains raw output, and verifies
wire digests and the selected direct path before accepting timing.

SIGKILL cannot run userspace affinity cleanup, but every timed child inherits
the canonical global-lock open-file description so overlapping authoritative
work remains excluded until its deepest surviving descendant exits.  HUP, INT,
and TERM are deferred until child/session, affinity, pair-lease, and global-lock
cleanup.
"""

from __future__ import annotations

import argparse
import ctypes
import datetime as dt
import errno
import fcntl
import hashlib
import importlib.util
import json
import math
import os
import platform
import re
import secrets
import selectors
import signal
import stat
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Sequence


SCHEMA = "leopard2-small-direct-abba/v5"
MATRIX_SCHEMA = "leopard2-small-direct-matrix/v1"
LARGE_MATRIX_SCHEMA = "leopard2-small-direct-full-matrix/v1"
TINY_MATRIX_SCHEMA = "leopard2-small-direct-tiny-matrix/v1"
CELL_SCHEMA = "leopard2-small-direct-cell/v3"
ENVELOPE_SCHEMA = "leopard2-small-direct-envelope/v2"
EXECUTION_SCHEMA = "leopard2-small-direct-immutable-execution/v2"
BUILD_SCHEMA_V4 = "leopard2-small-direct-build/v4"
BUILD_SCHEMA = "leopard2-small-direct-build/v5"
BINARY_IDENTITY_SCHEMA_V4 = "leopard2-small-direct-binary-identity/v4"
BINARY_IDENTITY_SCHEMA = "leopard2-small-direct-binary-identity/v5"
PROMOTION_SCHEMA = "leopard2-small-direct-promotion-analysis/v3"
EXACT_LEOPARD1_COMMIT = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
ORDER = ("baseline", "candidate", "candidate", "baseline")
DEFAULT_LOCK = Path("/tmp/leopard-gf8-authoritative.lock")
DEFAULT_ROUNDS = 3
CONTROL_SIGNALS = frozenset((signal.SIGHUP, signal.SIGINT, signal.SIGTERM))
MAX_GENERAL_FILE_BYTES = 256 << 20
MAX_SOURCE_FILE_BYTES = 256 << 20
MAX_EXECUTABLE_BYTES = 256 << 20
MAX_RESERVATION_JSON_BYTES = 64 << 10
MAX_BENCHMARK_STDOUT_BYTES = 8 << 20
MAX_BENCHMARK_STDERR_BYTES = 1 << 20
MAX_MATRIX_JSON_BYTES = 16 << 20
MAX_REQUEST_JSON_BYTES = 32 << 20
MAX_CELL_JSON_BYTES = 8 << 20
MAX_ENVELOPE_JSON_BYTES = 16 << 20
MAX_MANIFEST_JSON_BYTES = 256 << 20
LINUX_F_ADD_SEALS = getattr(fcntl, "F_ADD_SEALS", 1033)
LINUX_F_GET_SEALS = getattr(fcntl, "F_GET_SEALS", 1034)
LINUX_REQUIRED_EXECUTABLE_SEALS = (
    getattr(fcntl, "F_SEAL_SEAL", 0x0001) |
    getattr(fcntl, "F_SEAL_SHRINK", 0x0002) |
    getattr(fcntl, "F_SEAL_GROW", 0x0004) |
    getattr(fcntl, "F_SEAL_WRITE", 0x0008)
)
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
PROCESS_SUPPORT_RELATIVE = "experiments/leopard2/main_compare/run_abba.py"
SOURCE_FILES = (
    "CMakeLists.txt",
    "cmake/GenerateBenchmarkSourceAttestation.cmake",
    "cmake/Leopard2BenchmarkAttestation.cmake",
    "leopard.cpp", "leopard.h", "leopard2.cpp", "leopard2.h",
    "Leopard2Backend.cpp", "Leopard2Backend.h",
    "Leopard2BackendScalar.cpp", "Leopard2BackendSSSE3.cpp",
    "Leopard2BackendAVX2.cpp", "Leopard2BackendAVX2Xor.cpp",
    "Leopard2CpuFeatures.cpp",
    "Leopard2Direct.h", "Leopard2Dispatch.h", "Leopard2Plan.cpp",
    "Leopard2Plan.h",
    "LeopardCommon.cpp", "LeopardCommon.h",
    "LeopardFF8.cpp", "LeopardFF8.h", "LeopardFF16.cpp",
    "LeopardFF16.h", "Leopard2BackendGFNI.cpp",
    "bench/leopard2/benchmark.cpp",
    PROCESS_SUPPORT_RELATIVE,
    "tools/leopard2_build_provenance.py",
    "tools/leopard2_direct_encode_crossover.py",
)
MODE_COMPILE_DEFINITIONS = {
    # Mode zero is the source default and must remain absent from production
    # compile commands so the diagnostic does not perturb ordinary artifacts.
    "transform": None,
    "output": "-DLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=1",
    "source": "-DLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=2",
}
MODE_CACHE_VALUES = {"transform": "0", "output": "1", "source": "2"}
MODE_CACHE_VARIABLE = "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE"
BUILD_CONFIGURATION_ATTESTATION_SCHEMA_V2 = (
    "leopard2-benchmark-build-configuration-attestation/v2")
BUILD_CONFIGURATION_ATTESTATION_SCHEMA = (
    "leopard2-benchmark-build-configuration-attestation/v3")
BUILD_CONFIGURATION_VARIABLES_V2 = (
    "CMAKE_BUILD_TYPE",
    "CMAKE_GENERATOR",
    "CMAKE_CONFIGURATION_TYPES",
    "CMAKE_CXX_COMPILER",
    "CMAKE_CXX_FLAGS",
    "CMAKE_CXX_FLAGS_DEBUG",
    "CMAKE_CXX_FLAGS_RELEASE",
    "CMAKE_CXX_FLAGS_RELWITHDEBINFO",
    "CMAKE_CXX_FLAGS_MINSIZEREL",
    "ENABLE_OPENMP",
    "LEOPARD_ENABLE_GF8",
    "LEOPARD_ENABLE_GF16",
    "LEO2_BACKEND_VARIANT",
    "LEO2_BENCHMARK_GIT_EXECUTABLE",
    "LEO2_BUILD_BENCHMARKS",
    "LEO2_BUILD_TESTS",
    "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN",
    "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE",
    "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE",
)
BUILD_CONFIGURATION_VARIABLES = (
    *BUILD_CONFIGURATION_VARIABLES_V2[:-1],
    "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT",
    BUILD_CONFIGURATION_VARIABLES_V2[-1],
)
REQUIRED_DISABLED_EXPERIMENTS_V2 = {
    "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN": "OFF",
    "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE": "OFF",
}
REQUIRED_DISABLED_EXPERIMENTS = {
    **REQUIRED_DISABLED_EXPERIMENTS_V2,
    "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "OFF",
}
BENCHMARK_CONFIGURATION_DEFINITION = (
    "LEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256")
BENCHMARK_SOURCE_HEADER_DEFINITION = (
    "LEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER")
BENCHMARK_SOURCE_HEADER_RELATIVE = (
    "generated/leopard2-benchmark-attestation/"
    "leopard2_benchmark_source_attestation.h")


class EvidenceError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvidenceError(message)


def validated_inherited_descriptors(
        descriptors: Sequence[int], label: str) -> tuple[int, ...]:
    retained = tuple(descriptors)
    require(all(type(descriptor) is int and descriptor >= 0
                for descriptor in retained),
            "%s inherited descriptor set is invalid" % label)
    return tuple(sorted(set(retained)))


def load_process_containment_support() -> tuple[Any, Path]:
    """Load the repository's audited subreaper/procfs/pidfd implementation."""
    path = Path(__file__).resolve().parents[3] / PROCESS_SUPPORT_RELATIVE
    resolved = path.resolve(strict=True)
    specification = importlib.util.spec_from_file_location(
        "leopard2_small_direct_process_containment", resolved)
    require(specification is not None and specification.loader is not None,
            "cannot import process-containment support")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    try:
        specification.loader.exec_module(module)
    except Exception:
        sys.modules.pop(specification.name, None)
        raise
    require(callable(getattr(module, "LinuxDescendantContainment", None)) and
            callable(getattr(module, "_get_child_subreaper", None)),
            "process-containment support contract changed")
    return module, resolved


PROCESS_SUPPORT, PROCESS_SUPPORT_PATH = load_process_containment_support()


def mode_compile_definition(mode: str) -> str | None:
    require(mode in MODE_COMPILE_DEFINITIONS,
            "unknown small-direct comparison mode")
    return MODE_COMPILE_DEFINITIONS[mode]


def mode_compile_arguments(mode: str) -> list[str]:
    definition = mode_compile_definition(mode)
    return [] if definition is None else [definition]


def cpp_directive_view(source: str) -> str:
    """Expose logical C++ directives while blanking comments and literals."""
    if source.startswith("\ufeff"):
        source = source[1:]
    # Translation phase one maps implementation line endings to newlines.
    source = source.replace("\r\n", "\n").replace("\r", "\n")
    # C++11/14 translation phase one recognizes trigraphs.  Normalize all
    # spellings, not only ??=, because ??/ can create the phase-two splice
    # that joins a hidden directive.
    trigraphs = {
        "=": "#", "/": "\\", "'": "^", "(": "[", ")": "]",
        "!": "|", "<": "{", ">": "}", "-": "~",
    }
    source = re.sub(
        r"\?\?([=/\'()!<>-])",
        lambda match: trigraphs[match.group(1)], source)
    require(re.search(r"\\[ \t\v\f]+\n", source) is None,
            "leopard2.cpp uses unsupported whitespace after a line splice")
    # Translation phase two removes escaped physical newlines before comments
    # and preprocessing directives are recognized.  Apply it first so a
    # line-spliced comment or define/undef cannot evade the source contract.
    source = re.sub(r"\\(?:\r\n|\n|\r)", "", source)
    output = list(source)
    size = len(source)

    def blank(first: int, last: int) -> None:
        for index in range(first, last):
            if output[index] not in ("\n", "\r"):
                output[index] = " "

    index = 0
    while index < size:
        if source.startswith("//", index):
            end = source.find("\n", index + 2)
            end = size if end < 0 else end
            blank(index, end)
            index = end
            continue
        if source.startswith("/*", index):
            end = source.find("*/", index + 2)
            require(end >= 0, "leopard2.cpp contains an unterminated comment")
            end += 2
            blank(index, end)
            index = end
            continue
        if source.startswith('R"', index):
            opening = source.find(
                "(", index + 2, min(size, index + 19))
            require(opening >= 0,
                    "leopard2.cpp contains an invalid raw string")
            delimiter = source[index + 2:opening]
            require(len(delimiter) <= 16 and
                    all(character not in " ()\\\t\v\f\r\n"
                        for character in delimiter),
                    "leopard2.cpp contains an invalid raw-string delimiter")
            terminator = ")" + delimiter + '"'
            end = source.find(terminator, opening + 1)
            require(end >= 0,
                    "leopard2.cpp contains an unterminated raw string")
            end += len(terminator)
            blank(index, end)
            index = end
            continue
        if source[index] in ('"', "'"):
            quote = source[index]
            end = index + 1
            while end < size:
                if source[end] == "\\":
                    end += 2
                    continue
                if source[end] == quote:
                    end += 1
                    break
                require(source[end] not in "\r\n",
                        "leopard2.cpp contains an unterminated literal")
                end += 1
            require(end <= size and source[end - 1] == quote,
                    "leopard2.cpp contains an unterminated literal")
            blank(index, end)
            index = end
            continue
        if source.startswith("%:", index):
            # The C++ digraph %: is the preprocessing-token spelling of '#'.
            output[index] = "#"
            output[index + 1] = " "
            index += 2
            continue
        index += 1
    return "".join(output)


def validate_default_mode_source(source_root: Path) -> None:
    unused_identity, payload = stable_file_snapshot(
        source_root / "leopard2.cpp", "leopard2.cpp default-mode source",
        MAX_SOURCE_FILE_BYTES, require_canonical=True)
    try:
        source = payload.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise EvidenceError(
            "leopard2.cpp default-mode source is not UTF-8") from error
    require(all(
                character in "\r\n\t\v\f" or
                (ord(character) >= 32 and
                 not 127 <= ord(character) <= 159)
                for character in source),
            "leopard2.cpp contains unsupported control characters")
    directives = cpp_directive_view(source)
    pattern = re.compile(
        r"(?m)^[ \t\v\f]*#[ \t\v\f]*ifndef[ \t\v\f]+"
        r"LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE[ \t\v\f]*\r?\n"
        r"^[ \t\v\f]*#[ \t\v\f]*define[ \t\v\f]+"
        r"LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE[ \t\v\f]+0"
        r"[ \t\v\f]*\r?\n"
        r"^[ \t\v\f]*#[ \t\v\f]*endif[ \t\v\f]*$")
    mutations = re.findall(
        r"(?m)^[ \t\v\f]*#[ \t\v\f]*(?:define|undef)[ \t\v\f]+"
        r"LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE\b[^\r\n]*$",
        directives)
    matches = list(pattern.finditer(directives))
    conditional_depth = 0
    if matches:
        for conditional in re.finditer(
                r"(?m)^[ \t\v\f]*#[ \t\v\f]*"
                r"(if|ifdef|ifndef|endif|elif|elifdef|elifndef|else)\b",
                directives[:matches[0].start()]):
            directive = conditional.group(1)
            if directive in ("if", "ifdef", "ifndef"):
                conditional_depth += 1
            elif directive == "endif":
                conditional_depth -= 1
                require(conditional_depth >= 0,
                        "leopard2.cpp has unbalanced conditional directives")
    require(len(matches) == 1 and conditional_depth == 0 and
            mutations == [
                "#define LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE 0"],
            "production source does not default small-direct mode to zero")


def comparison_modes(baseline: str, candidate: str) -> dict[str, str]:
    require(baseline in MODE_COMPILE_DEFINITIONS and
            candidate in MODE_COMPILE_DEFINITIONS and
            baseline != candidate,
            "baseline and candidate modes must be distinct and known")
    return {"baseline": baseline, "candidate": candidate}


def strip_mode_definition(arguments: list[str], mode: str,
                          label: str) -> list[str]:
    definition = mode_compile_definition(mode)
    result = list(arguments)
    mentions = [
        (index, argument) for index, argument in enumerate(result)
        if MODE_CACHE_VARIABLE in argument
    ]
    if definition is None:
        require(not mentions,
                "%s production build contains a diagnostic mode" % label)
        return result
    require(len(mentions) == 1 and mentions[0][1] == definition and
            result.count(definition) == 1,
            "%s must contain its exact mode definition once" % label)
    result.remove(definition)
    return result


def build_record_contract(
        build_schema: str) -> tuple[
            str, tuple[str, ...], dict[str, str], str]:
    if build_schema == BUILD_SCHEMA:
        return (
            BUILD_CONFIGURATION_ATTESTATION_SCHEMA,
            BUILD_CONFIGURATION_VARIABLES,
            REQUIRED_DISABLED_EXPERIMENTS,
            BINARY_IDENTITY_SCHEMA,
        )
    if build_schema == BUILD_SCHEMA_V4:
        return (
            BUILD_CONFIGURATION_ATTESTATION_SCHEMA_V2,
            BUILD_CONFIGURATION_VARIABLES_V2,
            REQUIRED_DISABLED_EXPERIMENTS_V2,
            BINARY_IDENTITY_SCHEMA_V4,
        )
    raise EvidenceError("small-direct build schema is unsupported")


def build_configuration_digest(
        entries: dict[str, str], variables: Sequence[str]) -> str:
    require(tuple(variables) in (
                BUILD_CONFIGURATION_VARIABLES,
                BUILD_CONFIGURATION_VARIABLES_V2) and
            isinstance(entries, dict) and
            set(entries) == set(variables) and
            all(isinstance(entries[name], str)
                for name in variables),
            "effective CMake configuration entries have the wrong shape")
    material = "".join(
        "%s=%s\n" % (name, entries[name]) for name in variables)
    return sha256_bytes(material.encode("utf-8"))


def require_disabled_experiments(
        entries: dict[str, str], label: str,
        disabled: dict[str, str] = REQUIRED_DISABLED_EXPERIMENTS) -> None:
    require(all(entries.get(name) == value
                for name, value in disabled.items()),
            "%s does not disable unrelated direct-path experiments" % label)


def validate_effective_configuration_attestation(
        attestation: Any, mode: str, build_schema: str,
        label: str) -> dict[str, str]:
    expected_schema, variables, disabled, unused_identity_schema = \
        build_record_contract(build_schema)
    del unused_identity_schema
    require(isinstance(attestation, dict) and set(attestation) == {
                "entries", "path", "schema", "sha256"} and
            attestation.get("schema") == expected_schema and
            isinstance(attestation.get("path"), str) and
            bool(attestation["path"]) and
            isinstance(attestation.get("sha256"), str) and
            re.fullmatch(r"[0-9a-f]{64}", attestation["sha256"]) is not None,
            "%s attestation framing is invalid" % label)
    entries = attestation["entries"]
    expected_mode = MODE_CACHE_VALUES.get(mode)
    require(isinstance(entries, dict) and
            set(entries) == set(variables) and
            entries.get(MODE_CACHE_VARIABLE) == expected_mode,
            "%s attestation has the wrong diagnostic mode or selector "
            "keys differ" % label)
    require_disabled_experiments(entries, label, disabled)
    require(attestation["sha256"] ==
                build_configuration_digest(entries, variables),
            "%s attestation digest differs" % label)
    return entries


def require_exact_macro_argument(
        arguments: Sequence[str], variable: str,
        expected: str | None, label: str) -> None:
    mentions = [argument for argument in arguments if variable in argument]
    if expected is None:
        require(not mentions, "%s unexpectedly mentions %s" %
                (label, variable))
        return
    require(mentions == [expected] and arguments.count(expected) == 1,
            "%s does not bind exact %s definition" % (label, variable))


def benchmark_attestation_definitions(
        build_root: str | Path, configuration_sha256: str) -> tuple[str, str]:
    require(isinstance(configuration_sha256, str) and
            re.fullmatch(r"[0-9a-f]{64}", configuration_sha256) is not None,
            "effective configuration digest is malformed")
    header = Path(build_root) / BENCHMARK_SOURCE_HEADER_RELATIVE
    return (
        '-D%s="%s"' %
            (BENCHMARK_CONFIGURATION_DEFINITION, configuration_sha256),
        '-D%s="%s"' %
            (BENCHMARK_SOURCE_HEADER_DEFINITION, header),
    )


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
            "forked session leader pins and SIGSTOPs before snapshots, then "
            "execs a read-only descriptor for a write/grow/shrink-sealed "
            "memfd; bounded pipes "
            "retain stdout/stderr; waitid WNOWAIT preserves zombie schedstat "
            "and wait4 aggregates all child threads; Linux child-subreaper, "
            "procfs start-time identities, and pidfd signals contain setsid/"
            "double-fork descendants; HUP/INT/TERM remain blocked until the "
            "complete tree is killed, reaped, and proven empty"),
        "sigkill_limitation": (
            "SIGKILL cannot run userspace affinity restoration; kernel-held "
            "global flock remains inherited by surviving timed descendants"),
    }


def require_no_pending_control_signal() -> None:
    pending = set(signal.sigpending()) & CONTROL_SIGNALS
    require(not pending,
            "campaign interrupted by pending control signal(s): %s" %
            sorted(item.name for item in pending))


def canonical_bytes(value: Any) -> bytes:
    try:
        return json.dumps(
            value, sort_keys=True, separators=(",", ":"),
            allow_nan=False).encode("utf-8", errors="strict")
    except (TypeError, ValueError, UnicodeError) as error:
        raise EvidenceError(
            "evidence value is not canonical finite JSON") from error


def strict_json_bytes(payload: bytes, label: str) -> Any:
    """Decode evidence JSON without accepting ambiguous/non-standard values."""
    require(isinstance(payload, bytes), "%s JSON payload is not bytes" % label)

    def reject_constant(value: str) -> Any:
        raise EvidenceError(
            "%s JSON contains non-finite constant %s" % (label, value))

    def unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result = {}
        for key, value in pairs:
            if key in result:
                raise EvidenceError(
                    "%s JSON contains duplicate key %r" % (label, key))
            result[key] = value
        return result

    def finite_float(value: str) -> float:
        result = float(value)
        if not math.isfinite(result):
            raise EvidenceError(
                "%s JSON contains non-finite number %s" % (label, value))
        return result

    try:
        text = payload.decode("utf-8", errors="strict")
        result = json.loads(
            text, object_pairs_hook=unique_object,
            parse_constant=reject_constant, parse_float=finite_float)
        pending = [result]
        while pending:
            value = pending.pop()
            require(not (type(value) is float and not math.isfinite(value)),
                    "%s JSON contains a non-finite number" % label)
            if isinstance(value, dict):
                pending.extend(value.values())
            elif isinstance(value, list):
                pending.extend(value)
        return result
    except EvidenceError:
        raise
    except (UnicodeDecodeError, ValueError, OverflowError,
            RecursionError) as error:
        raise EvidenceError("%s is not strict JSON" % label) from error


def _stable_stat_fields(metadata: os.stat_result) -> tuple[int, ...]:
    return (
        metadata.st_dev, metadata.st_ino, metadata.st_mode, metadata.st_uid,
        metadata.st_nlink, metadata.st_size, metadata.st_mtime_ns,
        metadata.st_ctime_ns,
    )


def stable_file_snapshot(
        path: Path, label: str, maximum_bytes: int,
        *, require_canonical: bool = False) -> tuple[dict[str, Any], bytes]:
    """Read one regular file once, with a size bound and stable pathname.

    Every digest and parser must consume the returned bytes.  Opening the final
    component with O_NOFOLLOW and comparing both descriptor metadata and the
    post-read directory entry rejects symlinks, replacement, truncation, and
    in-place mutation during the snapshot.
    """
    require(type(maximum_bytes) is int and maximum_bytes >= 0,
            "%s has an invalid snapshot bound" % label)
    supplied = Path(path)
    try:
        resolved = supplied.resolve(strict=True)
    except (OSError, RuntimeError, ValueError) as error:
        raise EvidenceError("%s path cannot be resolved" % label) from error
    if require_canonical:
        require(supplied.is_absolute() and supplied == resolved,
                "%s path is not canonical" % label)
    flags = os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW
    flags |= getattr(os, "O_NONBLOCK", 0)
    try:
        descriptor = os.open(resolved, flags)
    except OSError as error:
        raise EvidenceError("%s cannot be opened safely" % label) from error
    try:
        before = os.fstat(descriptor)
        require(stat.S_ISREG(before.st_mode),
                "%s is not a regular file" % label)
        require(0 <= before.st_size <= maximum_bytes,
                "%s exceeds its retained byte bound" % label)
        remaining = before.st_size
        chunks = []
        while remaining:
            block = os.read(descriptor, min(1 << 20, remaining))
            require(block, "%s was truncated while being read" % label)
            chunks.append(block)
            remaining -= len(block)
        require(os.read(descriptor, 1) == b"",
                "%s grew while being read" % label)
        after = os.fstat(descriptor)
        require(_stable_stat_fields(before) == _stable_stat_fields(after),
                "%s changed while being read" % label)
        try:
            linked = os.lstat(resolved)
        except OSError as error:
            raise EvidenceError(
                "%s path changed while being read" % label) from error
        require(stat.S_ISREG(linked.st_mode) and
                linked.st_dev == after.st_dev and
                linked.st_ino == after.st_ino,
                "%s path was replaced while being read" % label)
        payload = b"".join(chunks)
        require(len(payload) == before.st_size,
                "%s stable snapshot length is inconsistent" % label)
        identity = {
            "path": str(resolved),
            "size": len(payload),
            "sha256": sha256_bytes(payload),
        }
        return identity, payload
    finally:
        os.close(descriptor)


def strict_json_file(
        path: Path, label: str,
        maximum_bytes: int = MAX_GENERAL_FILE_BYTES) -> Any:
    unused_identity, payload = stable_file_snapshot(
        path, label, maximum_bytes, require_canonical=True)
    return strict_json_bytes(payload, label)


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(
        path: Path, maximum_bytes: int = MAX_GENERAL_FILE_BYTES) -> str:
    identity, unused_payload = stable_file_snapshot(
        path, "file digest", maximum_bytes, require_canonical=True)
    return identity["sha256"]


def file_identity(
        path: Path, maximum_bytes: int = MAX_GENERAL_FILE_BYTES) -> dict[str, Any]:
    identity, unused_payload = stable_file_snapshot(
        path, "file identity", maximum_bytes, require_canonical=True)
    return identity


def git_output(
        source_root: Path, *arguments: str,
        inherited_descriptors: Sequence[int] = ()) -> str:
    completed = subprocess.run(
        ["/usr/bin/git", "-C", str(source_root), *arguments],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
        env=CHILD_ENV,
        pass_fds=validated_inherited_descriptors(
            inherited_descriptors, "Git identity helper"))
    require(completed.returncode == 0,
            "git identity command failed: %s" % completed.stderr.decode(
                "utf-8", errors="replace"))
    return completed.stdout.decode("utf-8", errors="strict").strip()


def source_identity(
        source_root: Path, *,
        inherited_descriptors: Sequence[int] = ()) -> dict[str, Any]:
    root = source_root.resolve(strict=True)
    files = {}
    for relative in SOURCE_FILES:
        path = root / relative
        identity, unused_payload = stable_file_snapshot(
            path, "source file %s" % relative, MAX_SOURCE_FILE_BYTES,
            require_canonical=True)
        files[relative] = {
            "size": identity["size"],
            "sha256": identity["sha256"],
        }
    status = git_output(
        root, "status", "--short", "--untracked-files=all",
        inherited_descriptors=inherited_descriptors)
    result = {
        "root": str(root),
        "head": git_output(
            root, "rev-parse", "HEAD",
            inherited_descriptors=inherited_descriptors),
        "head_tree": git_output(
            root, "rev-parse", "HEAD^{tree}",
            inherited_descriptors=inherited_descriptors),
        "status_short": status.splitlines() if status else [],
        "files": files,
    }
    result["snapshot_sha256"] = sha256_bytes(canonical_bytes(result))
    return result


def load_current_provenance_modules(
        source_root: Path) -> tuple[Any, Any]:
    root = source_root.resolve(strict=True)

    def load(relative: str, name: str) -> Any:
        path = root / relative
        spec = importlib.util.spec_from_file_location(name, path)
        require(spec is not None and spec.loader is not None,
                "cannot load current provenance module %s" % relative)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
        return module

    return (
        load("tools/leopard2_build_provenance.py",
             "leopard2_small_direct_build_provenance"),
        load("tools/leopard2_direct_encode_crossover.py",
             "leopard2_small_direct_configuration_attestation"),
    )


def validate_reproducible_build_proof(
        proof: Any, closure: dict[str, Any], label: str,
        provenance_module: Any | None = None) -> None:
    """Delegate to the shared complete immutable-replay proof validator."""
    try:
        if provenance_module is None:
            provenance_module, unused_attestation = \
                load_current_provenance_modules(
                    Path(closure["source_root"]))
            del unused_attestation
        provenance_module.validate_reproducible_build_proof(
            proof, closure, label=label)
    except Exception as error:
        raise EvidenceError(
            "%s clean reproducible-build proof is invalid: %s" %
            (label, error)) from error


def capture_current_build(
        executable: Path, source_root: Path, mode: str, *,
        inherited_descriptors: Sequence[int] = (),
        retained_reproducible_build: dict[str, Any] | None = None,
        ) -> dict[str, Any]:
    expected_mode = MODE_CACHE_VALUES.get(mode)
    require(expected_mode is not None, "unknown small-direct build mode")
    source = source_root.resolve(strict=True)
    binary = executable.resolve(strict=True)
    build = binary.parent
    provenance_module, attestation_module = \
        load_current_provenance_modules(source)
    try:
        closure = provenance_module.candidate_build_provenance(
            build, source, binary, "bench_leopard2",
            inherited_descriptors=inherited_descriptors)
        metadata = attestation_module.cmake_build_metadata(binary)
        cache_path = Path(closure["cmake_cache"]["path"])
        cache_identity, cache_bytes = provenance_module.file_snapshot(
            cache_path, "small-direct CMake cache")
        cache = provenance_module.parse_cmake_cache(cache_bytes)
        attestation_module.validate_build_configuration_attestation(
            metadata["effective_configuration_attestation"],
            build / attestation_module.BUILD_CONFIGURATION_RELATIVE_PATH)
        if retained_reproducible_build is None:
            reproducible_build = \
                provenance_module.verify_reproducible_candidate_build(
                    closure,
                    inherited_descriptors=inherited_descriptors)
        else:
            reproducible_build = retained_reproducible_build
    except Exception as error:
        raise EvidenceError(
            "current build provenance rejected %s: %s" %
            (binary, error)) from error

    require(cache_identity == closure["cmake_cache"],
            "current provenance modules disagree on CMake cache identity")
    require(cache.get(MODE_CACHE_VARIABLE) == expected_mode,
            "CMake cache does not bind requested small-direct mode")
    require_disabled_experiments(cache, "CMake cache")
    attested_configuration = metadata[
        "effective_configuration_attestation"]
    validate_effective_configuration_attestation(
        attested_configuration, mode, BUILD_SCHEMA,
        "effective-configuration sidecar")
    validated_cache = closure["validated_cache"]
    require(validated_cache.get("LEO2_BACKEND_VARIANT") == "avx2",
            "small-direct timing requires an explicit AVX2 build")
    require(validated_cache.get(MODE_CACHE_VARIABLE) == expected_mode,
            "validated build closure binds the wrong small-direct mode")
    require_disabled_experiments(
        validated_cache, "validated build closure")
    require(metadata["build_root"] == str(build) and
            metadata["executable"]["sha256"] ==
                closure["executable"]["sha256"] and
            metadata["executable"]["size"] ==
                closure["executable"]["size"],
            "effective-configuration metadata is not bound to benchmark")

    benchmark_definitions = benchmark_attestation_definitions(
        build, attested_configuration["sha256"])
    for record in closure["source_object_compile_closure"]:
        arguments = record["compile_entry"]["arguments"]
        if record["role"] == "archive" and \
                mode_bearing_archive_record(record):
            strip_mode_definition(
                arguments, mode,
                "archive compile command")
        else:
            strip_mode_definition(
                arguments, "transform",
                "backend-object or benchmark translation unit")
        if record["role"] == "benchmark":
            require_exact_macro_argument(
                arguments, BENCHMARK_CONFIGURATION_DEFINITION,
                benchmark_definitions[0], "benchmark compile command")
            require_exact_macro_argument(
                arguments, BENCHMARK_SOURCE_HEADER_DEFINITION,
                benchmark_definitions[1], "benchmark compile command")
        else:
            require_exact_macro_argument(
                arguments, BENCHMARK_CONFIGURATION_DEFINITION, None,
                "archive compile command")
            require_exact_macro_argument(
                arguments, BENCHMARK_SOURCE_HEADER_DEFINITION, None,
                "archive compile command")

    result = {
        "schema": BUILD_SCHEMA,
        "mode": mode,
        "cache_mode": expected_mode,
        "closure": closure,
        "reproducible_build": reproducible_build,
        "effective_configuration": metadata,
        "attestation_module": file_identity(
            source / "cmake/Leopard2BenchmarkAttestation.cmake"),
        "attestation_generator": file_identity(
            source / "cmake/GenerateBenchmarkSourceAttestation.cmake"),
        "provenance_module": file_identity(
            source / "tools/leopard2_build_provenance.py"),
        "configuration_reader": file_identity(
            source / "tools/leopard2_direct_encode_crossover.py"),
    }
    validate_reproducible_build_proof(
        result["reproducible_build"], closure, mode,
        provenance_module=provenance_module)
    result["digest"] = sha256_bytes(canonical_bytes(result))
    return result


def normalized_compile_arguments(
        record: dict[str, Any], mode: str,
        benchmark_definitions: Sequence[str] = ()) -> list[str]:
    arguments = list(record["compile_entry"]["arguments"])
    role = record.get("role")
    if role == "archive":
        require_exact_macro_argument(
            arguments, BENCHMARK_CONFIGURATION_DEFINITION, None,
            "archive build closure")
        require_exact_macro_argument(
            arguments, BENCHMARK_SOURCE_HEADER_DEFINITION, None,
            "archive build closure")
        if mode_bearing_archive_record(record):
            return strip_mode_definition(
                arguments, mode, "archive build closure")
        return strip_mode_definition(
            arguments, "transform", "backend-object build closure")
    require(role in ("archive", "benchmark"),
            "build closure contains an unknown compile role")
    strip_mode_definition(
        arguments, "transform",
        "backend-object or benchmark build closure")
    if role == "benchmark":
        require(len(benchmark_definitions) == 2,
                "benchmark attestation derivatives are missing")
        for variable, definition in zip((
                BENCHMARK_CONFIGURATION_DEFINITION,
                BENCHMARK_SOURCE_HEADER_DEFINITION),
                benchmark_definitions):
            require_exact_macro_argument(
                arguments, variable, definition,
                "benchmark compile closure")
            arguments.remove(definition)
    return arguments


def mode_bearing_archive_record(record: dict[str, Any]) -> bool:
    if record.get("role") != "archive":
        return False
    source = record.get("source", {}).get("path")
    require(isinstance(source, str) and source,
            "archive compile closure omits its source path")
    return Path(source).name == "leopard2.cpp"


def normalized_effective_configuration(
        attestation: dict[str, Any], mode: str,
        build_schema: str = BUILD_SCHEMA) -> dict[str, str]:
    entries = validate_effective_configuration_attestation(
        attestation, mode, build_schema,
        "effective CMake configuration")
    normalized = dict(entries)
    del normalized[MODE_CACHE_VARIABLE]
    return normalized


def normalized_validated_cache(
        cache: dict[str, Any], mode: str,
        build_schema: str = BUILD_SCHEMA) -> dict[str, Any]:
    require(isinstance(cache, dict),
            "validated CMake semantics are not an object")
    unused_schema, unused_variables, disabled, unused_identity_schema = \
        build_record_contract(build_schema)
    del unused_schema, unused_variables, unused_identity_schema
    expected_mode = MODE_CACHE_VALUES.get(mode)
    require(cache.get(MODE_CACHE_VARIABLE) == expected_mode,
            "validated CMake semantics bind the wrong diagnostic mode")
    general = "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT"
    require((build_schema == BUILD_SCHEMA and general in cache) or
            (build_schema == BUILD_SCHEMA_V4 and general not in cache),
            "validated CMake semantics use the wrong versioned selector keys")
    require_disabled_experiments(
        cache, "validated CMake semantics", disabled)
    normalized = dict(cache)
    del normalized[MODE_CACHE_VARIABLE]
    return normalized


def compare_current_builds(
        baseline: dict[str, Any], candidate: dict[str, Any],
        modes: dict[str, str]) -> dict[str, Any]:
    require(baseline["mode"] == modes["baseline"] and
            candidate["mode"] == modes["candidate"] and
            baseline.get("schema") == candidate.get("schema") and
            baseline.get("schema") in (BUILD_SCHEMA_V4, BUILD_SCHEMA),
            "build records use wrong comparison modes")
    build_schema = baseline["schema"]
    for key in ("attestation_module", "attestation_generator",
                "provenance_module", "configuration_reader"):
        require(baseline[key]["sha256"] == candidate[key]["sha256"] and
                baseline[key]["size"] == candidate[key]["size"],
                "builds use different %s bytes" % key)
    baseline_attestation = baseline["effective_configuration"][
        "effective_configuration_attestation"]
    candidate_attestation = candidate["effective_configuration"][
        "effective_configuration_attestation"]
    baseline_normalized = normalized_effective_configuration(
        baseline_attestation, modes["baseline"], build_schema)
    candidate_normalized = normalized_effective_configuration(
        candidate_attestation, modes["candidate"], build_schema)
    require(baseline_normalized == candidate_normalized,
            "effective CMake configurations differ beyond diagnostic mode")
    require(baseline_attestation["sha256"] !=
                candidate_attestation["sha256"],
            "mode-distinct effective configurations share one digest")
    baseline_benchmark_definitions = benchmark_attestation_definitions(
        baseline["closure"]["build_root"],
        baseline_attestation["sha256"])
    candidate_benchmark_definitions = benchmark_attestation_definitions(
        candidate["closure"]["build_root"],
        candidate_attestation["sha256"])

    def closure_map(build: dict[str, Any]) -> dict[tuple[str, str], Any]:
        source = Path(build["closure"]["source_root"])
        result = {}
        for record in build["closure"]["source_object_compile_closure"]:
            key = (
                record["role"],
                str(Path(record["source"]["path"]).relative_to(source)),
            )
            require(key not in result, "build closure contains duplicate source")
            result[key] = record
        return result

    baseline_map = closure_map(baseline)
    candidate_map = closure_map(candidate)
    require(set(baseline_map) == set(candidate_map),
            "diagnostic build source/object closures differ")
    differing_objects = []
    for key in sorted(baseline_map):
        left = baseline_map[key]
        right = candidate_map[key]
        require(left["source"]["sha256"] == right["source"]["sha256"] and
                left["source"]["size"] == right["source"]["size"] and
                left["flag_profile"] == right["flag_profile"],
                "diagnostic build source or ISA profile differs at %s" %
                key[1])
        require(normalized_compile_arguments(
                    left, modes["baseline"],
                    baseline_benchmark_definitions) ==
                normalized_compile_arguments(
                    right, modes["candidate"],
                    candidate_benchmark_definitions),
                "compile commands differ beyond diagnostic mode at %s" %
                key[1])
        bytes_equal = (
            left["object"]["sha256"] == right["object"]["sha256"] and
            left["object"]["size"] == right["object"]["size"])
        if key == ("archive", "leopard2.cpp"):
            require(not bytes_equal,
                    "diagnostic selector did not change leopard2.cpp object")
            differing_objects.append(key[1])
        else:
            require(bytes_equal,
                    "diagnostic unexpectedly changed object %s" % key[1])

    left_closure = baseline["closure"]
    right_closure = candidate["closure"]
    require(normalized_validated_cache(
                left_closure["validated_cache"], modes["baseline"],
                build_schema) ==
            normalized_validated_cache(
                right_closure["validated_cache"], modes["candidate"],
                build_schema),
            "validated CMake semantics differ beyond diagnostic mode")
    require(left_closure["compiler"]["sha256"] ==
                right_closure["compiler"]["sha256"] and
            left_closure["compiler_version_sha256"] ==
                right_closure["compiler_version_sha256"] and
            left_closure["archiver"]["sha256"] ==
                right_closure["archiver"]["sha256"],
            "diagnostic builds use different build tools")
    require(left_closure["archive"]["sha256"] !=
                right_closure["archive"]["sha256"] and
            left_closure["executable"]["sha256"] !=
                right_closure["executable"]["sha256"],
            "diagnostic archive/executable bytes are identical")
    result = {
        "schema": "leopard2-small-direct-build-comparison/v2",
        "modes": modes,
        "effective_configuration_sha256": {
            "baseline": baseline_attestation["sha256"],
            "candidate": candidate_attestation["sha256"],
        },
        "normalized_effective_configuration_sha256":
            sha256_bytes(canonical_bytes(baseline_normalized)),
        "differing_objects": differing_objects,
        "all_other_linked_objects_identical": True,
    }
    result["digest"] = sha256_bytes(canonical_bytes(result))
    return result


def _atomic_write_new(path: Path, payload: bytes, mode: int, label: str) -> None:
    """Durably publish one new regular file without replacing caller data."""
    require(isinstance(payload, bytes), "%s payload is not bytes" % label)
    parent = path.parent
    require(parent.is_absolute() and parent.resolve(strict=True) == parent,
            "%s parent directory is not canonical" % label)
    temporary = parent / (
        ".%s.tmp-%d-%s" % (path.name, os.getpid(), secrets.token_hex(8)))
    descriptor = None
    published = False
    final_created = False
    owned_identity = None
    try:
        descriptor = os.open(
            temporary,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC |
                os.O_NOFOLLOW,
            mode)
        owned_metadata = os.fstat(descriptor)
        owned_identity = (owned_metadata.st_dev, owned_metadata.st_ino)
        offset = 0
        while offset < len(payload):
            written = os.write(descriptor, payload[offset:])
            require(written > 0, "short %s write" % label)
            offset += written
        os.fchmod(descriptor, mode)
        os.fsync(descriptor)
        opened = os.fstat(descriptor)
        linked_temporary = os.lstat(temporary)
        require(stat.S_ISREG(opened.st_mode) and
                opened.st_dev == linked_temporary.st_dev and
                opened.st_ino == linked_temporary.st_ino,
                "%s temporary path changed" % label)
        try:
            os.link(temporary, path, follow_symlinks=False)
        except FileExistsError as error:
            raise EvidenceError(
                "%s destination already exists" % label) from error
        final_created = True
        linked_final = os.lstat(path)
        require(stat.S_ISREG(linked_final.st_mode) and
                linked_final.st_dev == opened.st_dev and
                linked_final.st_ino == opened.st_ino,
                "%s publication path changed" % label)
        published = True
        os.unlink(temporary)
        os.close(descriptor)
        descriptor = None
        directory = os.open(
            parent, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC |
                os.O_NOFOLLOW)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        if descriptor is not None:
            os.close(descriptor)
        try:
            temporary_metadata = os.lstat(temporary)
        except FileNotFoundError:
            pass
        else:
            if owned_identity == (
                    temporary_metadata.st_dev, temporary_metadata.st_ino):
                os.unlink(temporary)
        if final_created and not published:
            try:
                final_metadata = os.lstat(path)
            except FileNotFoundError:
                pass
            else:
                if owned_identity == (
                        final_metadata.st_dev, final_metadata.st_ino):
                    os.unlink(path)


def freeze_executable(
        origin: Path, artifact_root: Path, label: str) -> dict[str, Any]:
    source = origin.resolve(strict=True)
    directory = artifact_root / label
    directory.mkdir(parents=True, mode=0o700, exist_ok=False)
    os.chmod(directory, 0o700)
    frozen = directory / "bench_leopard2"
    origin_identity, origin_bytes = stable_file_snapshot(
        source, "%s origin benchmark" % label, MAX_EXECUTABLE_BYTES,
        require_canonical=True)
    _atomic_write_new(
        frozen, origin_bytes, 0o555, "%s frozen benchmark" % label)
    frozen_identity = file_identity(frozen, MAX_EXECUTABLE_BYTES)
    frozen_metadata = frozen.stat()
    require(frozen_metadata.st_uid == os.getuid() and
            frozen_metadata.st_nlink == 1 and
            stat.S_ISREG(frozen_metadata.st_mode) and
            origin_identity["sha256"] == frozen_identity["sha256"] and
            origin_identity["size"] == frozen_identity["size"],
            "frozen benchmark copy has unsafe metadata or changed bytes")
    return {
        "schema": "leopard2-small-direct-frozen-executable/v1",
        "origin": origin_identity,
        "frozen": frozen_identity,
    }


def descriptor_file_identity(
        descriptor: int, maximum_bytes: int = MAX_EXECUTABLE_BYTES
        ) -> dict[str, Any]:
    require(type(descriptor) is int and descriptor >= 0,
            "immutable executable descriptor is invalid")
    metadata = os.fstat(descriptor)
    require(stat.S_ISREG(metadata.st_mode) and
            0 <= metadata.st_size <= maximum_bytes,
            "immutable executable descriptor is not a bounded regular file")
    digest = hashlib.sha256()
    offset = 0
    while offset < metadata.st_size:
        block = os.pread(
            descriptor, min(1 << 20, metadata.st_size - offset), offset)
        require(block, "immutable executable descriptor was truncated")
        digest.update(block)
        offset += len(block)
    after = os.fstat(descriptor)
    require(_stable_stat_fields(metadata) == _stable_stat_fields(after),
            "immutable executable descriptor changed while hashing")
    try:
        seals = fcntl.fcntl(descriptor, LINUX_F_GET_SEALS)
        descriptor_flags = fcntl.fcntl(descriptor, fcntl.F_GETFL)
    except OSError as error:
        raise EvidenceError(
            "immutable executable descriptor is not a sealed memfd") \
            from error
    return {
        "size": metadata.st_size,
        "sha256": digest.hexdigest(),
        "mode": stat.S_IMODE(metadata.st_mode),
        "device": metadata.st_dev,
        "inode": metadata.st_ino,
        "nlink": metadata.st_nlink,
        "seals": seals,
        "access_mode": descriptor_flags & os.O_ACCMODE,
    }


def validate_execution_identity(
        value: Any, binary_identity: dict[str, Any]) -> dict[str, Any]:
    require(isinstance(value, dict) and
            set(value) == {"baseline", "candidate"},
            "immutable execution identity is incomplete")
    descriptors = set()
    for label in ("baseline", "candidate"):
        item = value.get(label)
        frozen = binary_identity[label]["executable"]["frozen"]
        require(isinstance(item, dict) and set(item) == {
                    "schema", "strategy", "descriptor", "command_path",
                    "logical_path", "sha256", "size", "mode", "seals",
                },
                "%s immutable execution identity is incomplete" % label)
        descriptor = item.get("descriptor")
        require(item.get("schema") == EXECUTION_SCHEMA and
                item.get("strategy") == "linux_sealed_memfd" and
                type(descriptor) is int and 3 <= descriptor <= 1_048_575 and
                descriptor not in descriptors and
                item.get("command_path") ==
                    "/proc/self/fd/%d" % descriptor and
                item.get("logical_path") == frozen["path"] and
                item.get("sha256") == frozen["sha256"] and
                item.get("size") == frozen["size"] and
                item.get("mode") == 0o500 and
                type(item.get("seals")) is int and
                item["seals"] & LINUX_REQUIRED_EXECUTABLE_SEALS ==
                    LINUX_REQUIRED_EXECUTABLE_SEALS,
                "%s immutable execution identity differs from frozen "
                "attestation" % label)
        descriptors.add(descriptor)
    return value


def validate_execution_descriptor(
        descriptor: int, identity: dict[str, Any]) -> None:
    current = descriptor_file_identity(descriptor)
    require(current["sha256"] == identity["sha256"] and
            current["size"] == identity["size"] and
            current["mode"] == identity["mode"] == 0o500 and
            current["nlink"] == 0 and
            current["access_mode"] == os.O_RDONLY and
            current["seals"] == identity["seals"] and
            current["seals"] & LINUX_REQUIRED_EXECUTABLE_SEALS ==
                LINUX_REQUIRED_EXECUTABLE_SEALS and
            identity["descriptor"] == descriptor and
            identity["command_path"] == "/proc/self/fd/%d" % descriptor,
            "immutable execution descriptor changed")


def create_executable_memfd(name: str) -> int:
    """Create a sealable executable memfd across old and new Linux kernels."""
    base_flags = getattr(os, "MFD_ALLOW_SEALING", 0x0002)
    executable_flag = getattr(os, "MFD_EXEC", 0x0010)

    def create(flags: int) -> int:
        wrapper = getattr(os, "memfd_create", None)
        if callable(wrapper):
            return wrapper(name, flags)
        require(sys.platform.startswith("linux"),
                "immutable executable staging requires Linux memfd_create")
        try:
            function = ctypes.CDLL(None, use_errno=True).memfd_create
        except (AttributeError, OSError) as error:
            raise EvidenceError(
                "immutable executable staging requires Linux memfd_create") \
                from error
        function.argtypes = (ctypes.c_char_p, ctypes.c_uint)
        function.restype = ctypes.c_int
        ctypes.set_errno(0)
        descriptor = function(name.encode("utf-8"), ctypes.c_uint(flags))
        if descriptor >= 0:
            return descriptor
        number = ctypes.get_errno()
        raise OSError(number or errno.EPERM,
                      os.strerror(number or errno.EPERM))

    try:
        return create(base_flags | executable_flag)
    except OSError as error:
        # MFD_EXEC was added after memfd_create.  Older kernels reject the
        # unknown bit but create executable memfds by default.
        if error.errno != errno.EINVAL:
            raise
    return create(base_flags)


class ImmutableFrozenExecutables:
    """Stage attested frozen bytes into kernel-immutable sealed memfds."""

    def __init__(
            self, campaign_root: Path,
            binary_identity: dict[str, Any]) -> None:
        self.campaign_root = campaign_root
        self.binary_identity = binary_identity
        self.descriptors: dict[str, int] = {}
        self.identities: dict[str, dict[str, Any]] = {}

    def __enter__(self) -> tuple[dict[str, Any], dict[str, int]]:
        require(sys.platform.startswith("linux"),
                "immutable executable staging requires Linux memfd_create")
        try:
            for label in ("baseline", "candidate"):
                frozen = self.binary_identity[label]["executable"]["frozen"]
                path = Path(frozen["path"])
                current, payload = stable_file_snapshot(
                    path, "%s frozen execution source" % label,
                    MAX_EXECUTABLE_BYTES, require_canonical=True)
                metadata = path.stat()
                require(current == frozen and
                        stat.S_IMODE(metadata.st_mode) == 0o555 and
                        metadata.st_uid == os.getuid() and
                        metadata.st_nlink == 1,
                        "%s frozen executable changed before immutable "
                        "staging" % label)
                writable = create_executable_memfd(
                    "leopard2-small-direct-" + label)
                descriptor = -1
                try:
                    offset = 0
                    while offset < len(payload):
                        written = os.write(writable, payload[offset:])
                        require(written > 0,
                                "short immutable executable write")
                        offset += written
                    os.fchmod(writable, 0o500)
                    os.fsync(writable)
                    fcntl.fcntl(
                        writable, LINUX_F_ADD_SEALS,
                        LINUX_REQUIRED_EXECUTABLE_SEALS)
                    descriptor = os.open(
                        "/proc/self/fd/%d" % writable, os.O_RDONLY)
                finally:
                    os.close(writable)
                try:
                    os.set_inheritable(descriptor, True)
                    descriptor_identity = descriptor_file_identity(descriptor)
                    require(
                        descriptor_identity["sha256"] == frozen["sha256"] and
                        descriptor_identity["size"] == frozen["size"] and
                        descriptor_identity["mode"] == 0o500 and
                        descriptor_identity["nlink"] == 0 and
                        descriptor_identity["access_mode"] == os.O_RDONLY and
                        descriptor_identity["seals"] &
                            LINUX_REQUIRED_EXECUTABLE_SEALS ==
                            LINUX_REQUIRED_EXECUTABLE_SEALS,
                        "%s immutable executable differs from frozen bytes" %
                        label)
                except BaseException:
                    os.close(descriptor)
                    raise
                self.descriptors[label] = descriptor
                self.identities[label] = {
                    "schema": EXECUTION_SCHEMA,
                    "strategy": "linux_sealed_memfd",
                    "descriptor": descriptor,
                    "command_path": "/proc/self/fd/%d" % descriptor,
                    "logical_path": frozen["path"],
                    "sha256": frozen["sha256"],
                    "size": frozen["size"],
                    "mode": 0o500,
                    "seals": descriptor_identity["seals"],
                }
            validate_execution_identity(
                self.identities, self.binary_identity)
            for label in ("baseline", "candidate"):
                validate_execution_descriptor(
                    self.descriptors[label], self.identities[label])
            return (
                json.loads(json.dumps(self.identities)),
                dict(self.descriptors),
            )
        except BaseException:
            for descriptor in self.descriptors.values():
                try:
                    os.close(descriptor)
                except OSError:
                    pass
            self.descriptors.clear()
            raise

    def __exit__(self, exc_type: Any, exc: Any, traceback: Any) -> None:
        del exc_type, exc, traceback
        failures = []
        descriptors = list(self.descriptors.items())
        self.descriptors.clear()
        for label, descriptor in descriptors:
            try:
                os.close(descriptor)
            except OSError as error:
                failures.append("%s fd %d: %s" % (
                    label, descriptor, error))
        require(not failures,
                "cannot close immutable executable descriptors: %s" %
                "; ".join(failures))


def atomic_json(
        path: Path, value: Any,
        maximum_bytes: int = MAX_GENERAL_FILE_BYTES) -> None:
    payload = canonical_bytes(value) + b"\n"
    require(len(payload) <= maximum_bytes,
            "atomic JSON exceeds its retained byte bound")
    try:
        destination = Path(os.path.abspath(os.fspath(path)))
    except (OSError, TypeError, ValueError) as error:
        raise EvidenceError("atomic JSON destination is invalid") from error
    _atomic_write_new(destination, payload, 0o600, "atomic JSON")


def validate_output_root_candidate(
        requested: Path, source_root: Path) -> tuple[Path, tuple[int, int]]:
    """Validate a fresh campaign root without creating or changing anything."""
    try:
        output = Path(os.path.abspath(os.fspath(requested)))
        parent = output.parent.resolve(strict=True)
    except (OSError, RuntimeError, TypeError, ValueError) as error:
        raise EvidenceError("campaign output path is invalid") from error
    require(output.parent == parent and output.name not in ("", ".", ".."),
            "campaign output parent is not canonical")
    try:
        output.relative_to(source_root)
    except ValueError:
        pass
    else:
        raise EvidenceError(
            "campaign output must be outside the source worktree")
    try:
        os.lstat(output)
    except FileNotFoundError:
        pass
    except OSError as error:
        raise EvidenceError("campaign output cannot be inspected") from error
    else:
        raise EvidenceError(
            "authoritative ABBA output must not already exist")
    metadata = os.stat(parent, follow_symlinks=False)
    require(stat.S_ISDIR(metadata.st_mode),
            "campaign output parent is not a directory")
    return output, (metadata.st_dev, metadata.st_ino)


def create_output_root(
        output: Path, expected_parent: tuple[int, int]) -> Path:
    """Atomically create the previously validated owned 0700 empty root."""
    parent = output.parent
    flags = os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW
    try:
        parent_descriptor = os.open(parent, flags)
    except OSError as error:
        raise EvidenceError(
            "campaign output parent cannot be opened safely") from error
    created = False
    try:
        parent_metadata = os.fstat(parent_descriptor)
        require((parent_metadata.st_dev, parent_metadata.st_ino) ==
                    expected_parent and stat.S_ISDIR(parent_metadata.st_mode),
                "campaign output parent changed before creation")
        try:
            os.stat(output.name, dir_fd=parent_descriptor,
                    follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            raise EvidenceError(
                "authoritative ABBA output appeared before creation")
        os.mkdir(output.name, mode=0o700, dir_fd=parent_descriptor)
        created = True
        root_descriptor = os.open(
            output.name, flags, dir_fd=parent_descriptor)
        try:
            os.fchmod(root_descriptor, 0o700)
            root_metadata = os.fstat(root_descriptor)
            linked_metadata = os.stat(
                output.name, dir_fd=parent_descriptor, follow_symlinks=False)
            require(stat.S_ISDIR(root_metadata.st_mode) and
                    root_metadata.st_uid == os.getuid() and
                    stat.S_IMODE(root_metadata.st_mode) == 0o700 and
                    root_metadata.st_dev == linked_metadata.st_dev and
                    root_metadata.st_ino == linked_metadata.st_ino,
                    "campaign output root is not an owned canonical 0700 "
                    "directory")
        finally:
            os.close(root_descriptor)
        require(output.resolve(strict=True) == output and
                not any(output.iterdir()),
                "new campaign output root is not canonical and empty")
        return output
    except Exception:
        if created:
            try:
                os.rmdir(output.name, dir_fd=parent_descriptor)
            except OSError:
                pass
        raise
    finally:
        os.close(parent_descriptor)


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
    if cell["loss"] <= 4:
        return "output_major"
    if mode == "transform":
        return "none"
    return mode + "_major"


def frozen_directional_scope() -> dict[str, Any]:
    return {
        "baseline_candidate_pairs": [
            ["transform", "output"],
            ["transform", "source"],
        ],
        "direct_head_to_head": (
            "Run output versus source only for overlapping winners; do "
            "not infer it by dividing separate noisy ratios."),
    }


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
        "directional_scope": frozen_directional_scope(),
        "same_source_promotion": {
            "target": "whole decode_execution",
            "candidate_ci95_low_at_least": 1.05,
            "neighbor_ci95_low_at_least": 0.98,
            "orientation": "baseline_time_over_candidate_time",
        },
        "exact_main_note": (
            "Run the cells marked exact_main_required with "
            "experiments/leopard2/main_compare/run_abba.py against exact "
            "Leopard main commit " + EXACT_LEOPARD1_COMMIT + ". "
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
                    for reuse_count in (1, 8, 64):
                        identifier = "full-k%d-r%d-b%d-l%d-u%d" % (
                            original_count, recovery_count, byte_count,
                            loss_count, reuse_count)
                        iterations, warmup = iteration_policy(byte_count)
                        cells.append({
                            "index": index,
                            "id": identifier,
                            "K": original_count,
                            "R": recovery_count,
                            "bytes": byte_count,
                            "loss": loss_count,
                            "batch": 1,
                            "reuse": reuse_count,
                            "iterations": iterations,
                            "warmup": warmup,
                            "seed": stable_seed(
                                identifier, LARGE_MATRIX_SCHEMA),
                            "estimated_peak_bytes": 6 *
                                (original_count + recovery_count) * byte_count +
                                (64 << 20),
                            "role": "loss4_neighbor" if loss_count == 4
                                else "loss5_to_8_target",
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
        "directional_scope": frozen_directional_scope(),
        "same_source_promotion": {
            "target": "whole decode_execution",
            "candidate_ci95_low_at_least": 1.05,
            "neighbor_ci95_low_at_least": 0.98,
            "orientation": "baseline_time_over_candidate_time",
        },
        "coverage": (
            "Every K=5..16, R=5..8, L=4..min(K,R), and requested "
            "64-byte through 64-KiB execution/tail boundary at plan reuse "
            "counts 1, 8, and 64."),
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
        "directional_scope": frozen_directional_scope(),
        "same_source_promotion": {
            "target": "whole decode_execution",
            "candidate_ci95_low_at_least": 1.05,
            "neighbor_ci95_low_at_least": 0.98,
            "orientation": "baseline_time_over_candidate_time",
        },
        "coverage": (
            "All required positive byte-tail boundaries below 64 bytes for "
            "representative K=5,8,16 and loss=4,5,8 shapes."),
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
    raise EvidenceError("unknown matrix preset")


def analyze_same_source_promotion(
        matrix: dict[str, Any],
        summaries: Sequence[dict[str, Any]],
        modes: dict[str, str]) -> dict[str, Any]:
    """Apply the frozen lower-confidence-bound same-source screen.

    Target cells must clear the 5% lower-bound win and loss-four neighbors
    must keep their lower bound at or above 0.98.  This runner cannot bind an
    exact historical Leopard1 executable, so even a complete passing screen is
    explicitly non-promotional until main-compare evidence is integrated.
    """
    require(isinstance(modes, dict) and
            set(modes) == {"baseline", "candidate"},
            "promotion analysis has no exact requested mode pair")
    requested_modes = comparison_modes(
        modes.get("baseline"), modes.get("candidate"))
    directional_scope = matrix.get("directional_scope")
    require(isinstance(directional_scope, dict) and
            set(directional_scope) == {
                "baseline_candidate_pairs", "direct_head_to_head"} and
            directional_scope["baseline_candidate_pairs"] == [
                ["transform", "output"], ["transform", "source"]] and
            isinstance(directional_scope["direct_head_to_head"], str) and
            directional_scope["direct_head_to_head"],
            "matrix has an invalid directional promotion scope")
    requested_pair = [
        requested_modes["baseline"], requested_modes["candidate"]]
    same_source_direction_authorized = requested_pair in \
        directional_scope["baseline_candidate_pairs"]
    policy = matrix.get("same_source_promotion")
    require(isinstance(policy, dict) and set(policy) == {
                "target", "candidate_ci95_low_at_least",
                "neighbor_ci95_low_at_least", "orientation",
            } and
            policy["target"] == "whole decode_execution" and
            policy["orientation"] == "baseline_time_over_candidate_time",
            "matrix has an invalid same-source promotion contract")
    target_threshold = policy["candidate_ci95_low_at_least"]
    neighbor_threshold = policy["neighbor_ci95_low_at_least"]
    require(type(target_threshold) in (int, float) and
            math.isfinite(target_threshold) and target_threshold > 1 and
            type(neighbor_threshold) in (int, float) and
            math.isfinite(neighbor_threshold) and
            0 < neighbor_threshold <= 1,
            "matrix has invalid same-source promotion thresholds")
    cells = matrix.get("cells")
    require(isinstance(cells, list) and cells,
            "matrix has no promotion cells")
    matrix_ids = [cell.get("id") for cell in cells]
    require(all(isinstance(identifier, str) and identifier
                for identifier in matrix_ids) and
            len(matrix_ids) == len(set(matrix_ids)),
            "matrix has invalid promotion cell IDs")
    known_cells = {
        cell["id"]: cell for cell in cells
    }

    checks = []
    seen = set()
    for summary in summaries:
        require(isinstance(summary, dict) and
                isinstance(summary.get("cell"), dict),
                "promotion summary is malformed")
        cell = summary["cell"]
        identifier = cell.get("id")
        require(identifier in known_cells and
                cell == known_cells[identifier] and
                identifier not in seen,
                "promotion summary is not bound to one frozen matrix cell")
        seen.add(identifier)
        role = cell.get("role")
        require(role in ("loss4_neighbor", "loss5_to_8_target"),
                "promotion cell has an unknown role")
        threshold = neighbor_threshold if role == "loss4_neighbor" \
            else target_threshold
        ratios = summary.get("metric_ratios")
        execution = ratios.get("execution") \
            if isinstance(ratios, dict) else None
        interval = execution.get("ci95") \
            if isinstance(execution, dict) else None
        require(isinstance(interval, list) and len(interval) == 2,
                "promotion summary has no execution confidence interval")
        low, high = interval
        has_interval = (
            type(low) in (int, float) and math.isfinite(low) and low > 0 and
            type(high) in (int, float) and math.isfinite(high) and
            high >= low)
        checks.append({
            "cell_id": identifier,
            "role": role,
            "required_ci95_low": threshold,
            "observed_ci95": interval,
            "passes": has_interval and low >= threshold,
        })

    selected_ids = [check["cell_id"] for check in checks]
    complete = selected_ids == matrix_ids
    all_pass = bool(checks) and all(check["passes"] for check in checks)
    roles = {check["role"] for check in checks}
    same_source_screen_eligible = (
        same_source_direction_authorized and complete and all_pass and
        roles == {"loss4_neighbor", "loss5_to_8_target"})
    result = {
        "schema": PROMOTION_SCHEMA,
        "policy": policy,
        "directional_scope": directional_scope,
        "requested_mode_pair": requested_pair,
        "same_source_direction_authorized":
            same_source_direction_authorized,
        "same_source_screen_eligible": same_source_screen_eligible,
        "exact_leopard1_evidence": {
            "required": True,
            "status": "not_bound",
            "required_commit": EXACT_LEOPARD1_COMMIT,
            "note": (
                "This runner compares same-source Leopard2 modes only. "
                "Exact Leopard1 evidence must be produced and bound by the "
                "main-compare runner before promotion."),
        },
        "promotion_authorized": False,
        "promotion_ineligibility_reason": (
            "exact Leopard1 evidence is not bound to this same-source "
            "screen"),
        "checks": checks,
        "complete_frozen_matrix": complete,
        "all_evaluated_cells_pass": all_pass,
        "promotion_eligible": False,
    }
    result["digest"] = sha256_bytes(canonical_bytes(result))
    return result


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


def validate_reservation_payload(
        value: Any, cpu: int, sibling: int) -> None:
    require(isinstance(value, dict) and set(value) == {
                "benchmark_cpu", "nonce", "owner", "reserved_sibling",
                "schema", "status",
            },
            "CPU reservation payload has unexpected or missing fields")
    require(value["schema"] == "leopard2-cpu-reservation/v1" and
            value["status"] == "held" and
            type(value["benchmark_cpu"]) is int and
            value["benchmark_cpu"] == cpu and
            type(value["reserved_sibling"]) is int and
            value["reserved_sibling"] == sibling and
            isinstance(value["owner"], str) and value["owner"] and
            isinstance(value["nonce"], str) and value["nonce"],
            "CPU reservation does not match the requested held pair")


def reservation_identity(path: Path, cpu: int, sibling: int) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    identity, payload = stable_file_snapshot(
        resolved, "CPU reservation", MAX_RESERVATION_JSON_BYTES,
        require_canonical=True)
    value = strict_json_bytes(payload, "CPU reservation")
    validate_reservation_payload(value, cpu, sibling)
    return {"file": identity, "value": value}


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


def finite_timing_number(
        value: Any, label: str, *, positive: bool) -> float:
    require(type(value) in (int, float),
            "%s is not numeric" % label)
    try:
        result = float(value)
    except (OverflowError, TypeError, ValueError) as error:
        raise EvidenceError("%s is not bounded" % label) from error
    require(math.isfinite(result) and
            (result > 0 if positive else result >= 0),
            "%s is not a finite %s value" % (
                label, "positive" if positive else "nonnegative"))
    return result


def timing_close(actual: float, expected: float) -> bool:
    # The benchmark serializes doubles with decimal stream formatting.  This
    # tolerance covers that round trip only; every statistic is still derived
    # independently from the retained arrays.
    return math.isfinite(actual) and math.isfinite(expected) and \
        abs(actual - expected) <= max(
        0.000002, abs(expected) * 0.000002)


def validate_timing_summary(
        value: Any, iterations: int, *, setup: bool) -> list[float]:
    require(isinstance(value, dict), "timing summary is not an object")
    suffix = "" if setup else "_per_batch_call"
    names = {
        "median": "median_us" + suffix,
        "mad": "mad_us" + suffix,
        "minimum": "minimum_us" + suffix,
        "maximum": "maximum_us" + suffix,
        "samples": "samples_us" if setup
            else "samples_us_per_batch_call",
    }
    require(set(value) == set(names.values()),
            "timing summary has unexpected or missing fields")
    retained = value[names["samples"]]
    require(isinstance(retained, list) and len(retained) == iterations,
            "%s does not contain exactly %d samples" % (
                names["samples"], iterations))
    samples = [
        finite_timing_number(
            item, names["samples"], positive=not setup)
        for item in retained
    ]
    median = statistics.median(samples)
    derived = {
        names["median"]: median,
        names["mad"]: statistics.median(
            [abs(item - median) for item in samples]),
        names["minimum"]: min(samples),
        names["maximum"]: max(samples),
    }
    for name, expected in derived.items():
        actual = finite_timing_number(
            value.get(name), name,
            positive=(name != names["mad"] and (not setup or expected > 0)))
        require(timing_close(actual, expected),
                "%s is not derived from retained samples" % name)
    return samples


def validate_execution_timing(
        value: Any, iterations: int,
        input_name: str, output_name: str,
        input_bytes: int, output_bytes: int) -> list[float]:
    require(isinstance(value, dict) and set(value) == {
                "median_us_per_batch_call", "mad_us_per_batch_call",
                "minimum_us_per_batch_call", "maximum_us_per_batch_call",
                "samples_us_per_batch_call", input_name, output_name,
            },
            "execution timing summary has unexpected or missing fields")
    summary_only = {
        key: value[key] for key in (
            "median_us_per_batch_call", "mad_us_per_batch_call",
            "minimum_us_per_batch_call", "maximum_us_per_batch_call",
            "samples_us_per_batch_call",
        )
    }
    samples = validate_timing_summary(
        summary_only, iterations, setup=False)
    median = statistics.median(samples)
    for name, byte_count in (
            (input_name, input_bytes), (output_name, output_bytes)):
        actual = finite_timing_number(
            value.get(name), name, positive=True)
        expected = byte_count / (median * 1000.0)
        require(math.isfinite(expected) and expected > 0 and
                timing_close(actual, expected),
                "%s is not derived from retained median and bytes" % name)
    return samples


def validate_result(result: Any, cell: dict[str, Any], mode: str) -> dict[str, Any]:
    require(isinstance(result, dict), "benchmark output is not an object")
    require(result.get("schema") == "leopard2-benchmark-v6",
            "benchmark did not emit the direct-executor v6 schema")
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
        "report_direct_executor": True,
    }
    require(set(parameters) ==
            set(expected_parameters) | {"missing_original_indices"},
            "benchmark parameters have unexpected or missing fields")
    for key, expected in expected_parameters.items():
        actual = parameters.get(key)
        if type(expected) is bool:
            exact = type(actual) is bool and actual is expected
        elif type(expected) is int:
            exact = type(actual) is int and actual == expected
        else:
            exact = type(actual) is type(expected) and actual == expected
        require(exact,
                "benchmark parameter %s does not match request" % key)
    expected_resolved_keys = {
        "profile", "field", "backend", "thread_count", "parent_count",
        "padded_side", "selected_decode_path", "selected_decode_rule",
        "decode_required_work_slots", "decode_aligned_prefix_bytes",
        "decode_tail_bytes", "decode_rounded_bytes",
        "decode_multi_item_batch", "selected_direct_executor",
    }
    require(set(resolved) == expected_resolved_keys,
            "benchmark resolved fields are unexpected or incomplete")
    require(resolved.get("profile") == "legacy_high_v1",
            "benchmark resolved the wrong profile")
    require(resolved.get("field") == "gf8", "benchmark resolved the wrong field")
    require(resolved.get("backend") == "avx2",
            "benchmark did not resolve the AVX2 backend")
    require(type(resolved.get("thread_count")) is int and
            resolved["thread_count"] == 1,
            "benchmark did not resolve one execution thread")
    padded_side = ceil_power_of_two(cell["R"])
    parent_count = ceil_power_of_two(cell["K"] + padded_side)
    require(type(resolved.get("parent_count")) is int and
            type(resolved.get("padded_side")) is int and
            resolved["parent_count"] == parent_count and
            resolved["padded_side"] == padded_side,
            "benchmark resolved the wrong direct-repair parent geometry")
    for key in (
            "decode_required_work_slots", "decode_aligned_prefix_bytes",
            "decode_tail_bytes", "decode_rounded_bytes"):
        require(type(resolved.get(key)) is int and resolved[key] >= 0,
                "benchmark resolved %s is not an exact nonnegative integer" %
                key)
    aligned = resolved["decode_aligned_prefix_bytes"]
    tail = resolved["decode_tail_bytes"]
    rounded = resolved["decode_rounded_bytes"]
    require(
        tail < 64 and aligned % 64 == 0 and
        aligned + tail == cell["bytes"] and
        rounded == (aligned if tail == 0 else aligned + 64) and
        type(resolved.get("decode_multi_item_batch")) is bool and
        resolved["decode_multi_item_batch"] is False,
        "benchmark resolved invalid byte/batch geometry")
    expected_executor = expected_direct_executor(mode, cell)
    selected_path = resolved.get("selected_decode_path")
    selected_rule = resolved.get("selected_decode_rule")
    require(resolved.get("selected_direct_executor") == expected_executor,
            "benchmark reported the wrong direct-repair loop order")
    if expected_executor == "none":
        require(selected_path in ("generic", "materialized", "tiled") and
                selected_rule not in ("no_op", "direct", "unsupported_profile"),
                "benchmark did not select transform repair")
    else:
        require(selected_path == "direct" and selected_rule == "direct" and
                resolved["decode_required_work_slots"] == 0,
                "benchmark did not select direct repair")
    require(set(correctness) == {
                "leopard2_round_trip", "legacy_comparison"} and
            correctness.get("leopard2_round_trip") is True and
            correctness.get("legacy_comparison") is None,
            "benchmark round trip failed")
    require(set(digests) == {
                "algorithm", "original_data", "transmitted_parity",
                "recovered_originals"} and
            digests.get("algorithm") == "fnv1a64",
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
    require(set(metrics) == {
                "codec_setup", "encode_execution", "decode_plan_setup",
                "decode_execution", "decode_amortized_at_reuse",
                "rate_semantics",
            } and
            metrics.get("rate_semantics") ==
                "offered_received counts all non-null shard pointers supplied; "
                "a plan may read a deterministic subset",
            "benchmark timing metrics are missing or changed")
    iterations = cell["iterations"]
    codec_setup_samples = validate_timing_summary(
        metrics["codec_setup"], iterations, setup=True)
    encode_samples = validate_execution_timing(
        metrics["encode_execution"], iterations,
        "input_GB_per_s", "parity_output_GB_per_s",
        cell["K"] * cell["bytes"], cell["R"] * cell["bytes"])
    setup_samples = validate_timing_summary(
        metrics["decode_plan_setup"], iterations, setup=True)
    execution_samples = validate_execution_timing(
        metrics["decode_execution"], iterations,
        "offered_received_GB_per_s", "repaired_output_GB_per_s",
        (cell["K"] - cell["loss"] + cell["R"]) * cell["bytes"],
        cell["loss"] * cell["bytes"])
    execution_median = statistics.median(execution_samples)
    setup_median = statistics.median(setup_samples)
    amortized = metrics["decode_amortized_at_reuse"]
    require(isinstance(amortized, dict) and set(amortized) == {
                "reuse_count", "derived_median_us_per_batch_call",
                "offered_received_GB_per_s", "repaired_output_GB_per_s",
            } and type(amortized.get("reuse_count")) is int and
            amortized["reuse_count"] == cell["reuse"],
            "amortized decode summary has unexpected fields or reuse")
    expected_amortized = execution_median + \
        setup_median / cell["reuse"]
    require(timing_close(
                finite_timing_number(
                    amortized.get("derived_median_us_per_batch_call"),
                    "amortized decode median", positive=True),
                expected_amortized),
            "amortized decode median is not derived from retained samples")
    for name, byte_count in (
            ("offered_received_GB_per_s",
             (cell["K"] - cell["loss"] + cell["R"]) * cell["bytes"]),
            ("repaired_output_GB_per_s",
             cell["loss"] * cell["bytes"])):
        expected_rate = byte_count / (expected_amortized * 1000.0)
        require(timing_close(
                    finite_timing_number(
                        amortized.get(name), "amortized " + name,
                        positive=True),
                    expected_rate) and
                math.isfinite(expected_rate) and expected_rate > 0,
                "amortized %s is not derived from retained samples" % name)
    return {
        "median_us": execution_median,
        "plan_setup_us": setup_median,
        "decode_execution_samples_us_per_batch_call": execution_samples,
        "decode_plan_setup_samples_us": setup_samples,
        "codec_setup_samples_us": codec_setup_samples,
        "encode_execution_samples_us_per_batch_call": encode_samples,
        "digests": digests,
        "missing_original_indices": missing,
        "decode_path": selected_path,
        "decode_rule": selected_rule,
        "build_bound_executor": resolved["selected_direct_executor"],
    }


def ratio_summary(invocations: list[dict[str, Any]], rounds: int,
                  metric) -> dict[str, Any]:
    def finite_ratio(logarithm: float) -> float:
        require(math.isfinite(logarithm),
                "ABBA log contrast is not finite")
        try:
            result = math.exp(logarithm)
        except OverflowError as error:
            raise EvidenceError(
                "ABBA ratio is outside the finite numeric range") from error
        require(math.isfinite(result) and result > 0,
                "ABBA ratio is outside the finite numeric range")
        return result

    contrasts = []
    round_ratios = []
    for round_index in range(rounds):
        group = invocations[round_index * len(ORDER):(round_index + 1) * len(ORDER)]
        baseline = [metric(item["normalized"]) for item in group
                    if item["implementation"] == "baseline"]
        candidate = [metric(item["normalized"]) for item in group
                     if item["implementation"] == "candidate"]
        require(len(baseline) == 2 and len(candidate) == 2 and
                all(type(value) in (int, float) and
                    math.isfinite(value) and value > 0
                    for value in baseline + candidate),
                "ABBA round has an invalid metric")
        contrast = (sum(math.log(value) for value in baseline) -
                    sum(math.log(value) for value in candidate)) / 2.0
        contrasts.append(contrast)
        round_ratios.append(finite_ratio(contrast))
    mean_log = statistics.mean(contrasts)
    if rounds == 1:
        interval = [None, None]
    else:
        standard_error = statistics.stdev(contrasts) / math.sqrt(rounds)
        critical = T_CRITICAL_95[rounds - 1]
        interval = [
            finite_ratio(mean_log - critical * standard_error),
            finite_ratio(mean_log + critical * standard_error),
        ]
    return {
        "round_ratios": round_ratios,
        "geometric_mean_ratio": finite_ratio(mean_log),
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
        "--retain-samples", "--report-decode-path",
        "--report-direct-executor", "--json", "-",
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


def _wait_status_return_code(status: int) -> int:
    if os.WIFEXITED(status):
        return os.WEXITSTATUS(status)
    if os.WIFSIGNALED(status):
        return -os.WTERMSIG(status)
    raise EvidenceError("benchmark child exited in an unexpected state")


class ForkedProcessHandle:
    """Minimal Popen-compatible ownership handle for audited containment."""

    def __init__(self, pid: int) -> None:
        self.pid = pid
        self.returncode: int | None = None

    def mark_reaped(self, status: int) -> None:
        require(self.returncode is None,
                "benchmark child was marked reaped twice")
        self.returncode = _wait_status_return_code(status)

    def poll(self) -> int | None:
        if self.returncode is not None:
            return self.returncode
        try:
            waited, status = os.waitpid(self.pid, os.WNOHANG)
        except ChildProcessError:
            # A normal path may already have consumed wait4 before assigning
            # the adapter result.  Treat that identity as gone; containment
            # still proves all adopted descendants empty.
            self.returncode = -int(signal.SIGKILL)
            return self.returncode
        if waited == 0:
            return None
        require(waited == self.pid,
                "containment reaped an unexpected child")
        self.returncode = _wait_status_return_code(status)
        return self.returncode

    def wait(self, timeout: float | None = None) -> int:
        deadline = None if timeout is None else time.monotonic() + timeout
        while True:
            result = self.poll()
            if result is not None:
                return result
            if deadline is not None and time.monotonic() >= deadline:
                raise subprocess.TimeoutExpired(
                    ["/proc/self/fd/immutable"], timeout)
            time.sleep(0.001)


def _write_all(descriptor: int, payload: bytes) -> None:
    view = memoryview(payload)
    while view:
        written = os.write(descriptor, view)
        require(written > 0, "short retained benchmark output write")
        view = view[written:]


def prepare_gated_resources(
        stdout_path: Path, stderr_path: Path,
) -> tuple[int, int, int, int, int, int, set[signal.Signals], Any]:
    """Create signal-safe bounded-output resources with rollback on failure."""
    stdout_fd = stderr_fd = -1
    stdout_read = stdout_write = -1
    stderr_read = stderr_write = -1
    previous_mask: set[signal.Signals] | None = None
    selector = None
    try:
        previous_mask = signal.pthread_sigmask(
            signal.SIG_BLOCK, set())
        blocked_from = signal.pthread_sigmask(
            signal.SIG_BLOCK, CONTROL_SIGNALS)
        require(blocked_from == previous_mask,
                "signal mask changed during gated resource setup")
        selector = selectors.DefaultSelector()
        output_flags = \
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC | \
            os.O_NOFOLLOW
        stdout_fd = os.open(stdout_path, output_flags, 0o600)
        stderr_fd = os.open(stderr_path, output_flags, 0o600)
        stdout_read, stdout_write = os.pipe2(os.O_CLOEXEC)
        stderr_read, stderr_write = os.pipe2(os.O_CLOEXEC)
        return (
            stdout_fd, stderr_fd, stdout_read, stdout_write,
            stderr_read, stderr_write, previous_mask, selector,
        )
    except BaseException as primary:
        cleanup_failures = []
        if selector is not None:
            try:
                selector.close()
            except BaseException as error:
                cleanup_failures.append("selector: %s" % error)
        for descriptor in (
                stdout_read, stdout_write, stderr_read, stderr_write,
                stdout_fd, stderr_fd):
            if descriptor < 0:
                continue
            try:
                os.close(descriptor)
            except OSError as error:
                cleanup_failures.append(
                    "fd %d: %s" % (descriptor, error))
        if previous_mask is not None:
            try:
                signal.pthread_sigmask(
                    signal.SIG_SETMASK, previous_mask)
            except BaseException as error:
                cleanup_failures.append("signal mask: %s" % error)
        if cleanup_failures:
            raise EvidenceError(
                "gated resource setup failed and rollback was incomplete: "
                "%s; primary failure: %s: %s" % (
                    "; ".join(cleanup_failures),
                    type(primary).__name__, primary)) from primary
        raise


def run_gated_benchmark(
        command: list[str], cpu: int, sibling: int, timeout: float,
        stdout_path: Path, stderr_path: Path,
        child_signal_mask: set[signal.Signals] | None = None,
        executable_descriptor: int | None = None,
        execution_identity: dict[str, Any] | None = None,
        campaign_lock_descriptor: int | None = None,
) -> dict[str, Any]:
    # Reject an invalid executable binding before creating retained output
    # files or pipes.  Besides being cheaper, this keeps a malformed identity
    # from leaking descriptors on the pre-containment failure path.
    require((executable_descriptor is None) ==
            (execution_identity is None),
            "immutable execution descriptor and identity must be supplied "
            "together")
    require(campaign_lock_descriptor is None or
            (type(campaign_lock_descriptor) is int and
             campaign_lock_descriptor >= 0),
            "campaign lock descriptor is invalid")
    if campaign_lock_descriptor is not None:
        lock_metadata = os.fstat(campaign_lock_descriptor)
        require(stat.S_ISREG(lock_metadata.st_mode),
                "campaign lock descriptor is not a regular file")
    if executable_descriptor is not None:
        validate_execution_descriptor(
            executable_descriptor, execution_identity)
        require(command[0] == execution_identity["command_path"],
                "benchmark command does not execute immutable descriptor")
    output_directory = stdout_path.parent
    directory_metadata = os.stat(output_directory, follow_symlinks=False)
    require(output_directory.is_absolute() and
            output_directory.resolve(strict=True) == output_directory and
            stat.S_ISDIR(directory_metadata.st_mode) and
            directory_metadata.st_uid == os.getuid() and
            stat.S_IMODE(directory_metadata.st_mode) == 0o700,
            "benchmark raw-output directory is not canonical and owned")
    (stdout_fd, stderr_fd, stdout_read, stdout_write,
     stderr_read, stderr_write, previous_mask, selector) = \
        prepare_gated_resources(stdout_path, stderr_path)
    effective_child_mask = previous_mask if child_signal_mask is None \
        else child_signal_mask
    pid = None
    handle = None
    result = None
    try:
        try:
            with PROCESS_SUPPORT.LinuxDescendantContainment() as containment:
                try:
                    pid = os.fork()
                except Exception:
                    raise
                if pid == 0:
                    try:
                        os.close(stdout_read)
                        os.close(stderr_read)
                        os.dup2(stdout_write, 1)
                        os.dup2(stderr_write, 2)
                        os.close(stdout_write)
                        os.close(stderr_write)
                        os.close(stdout_fd)
                        os.close(stderr_fd)
                        os.setsid()
                        os.sched_setaffinity(0, {cpu})
                        if campaign_lock_descriptor is not None:
                            # The campaign owner closes rather than explicitly
                            # unlocking this shared open-file description.
                            # Keep it across exec so coordinator SIGKILL cannot
                            # release the authoritative lease while timed work
                            # or one of its descendants survives.
                            os.set_inheritable(
                                campaign_lock_descriptor, True)
                        os.kill(os.getpid(), signal.SIGSTOP)
                        signal.pthread_sigmask(
                            signal.SIG_SETMASK, effective_child_mask)
                        os.execve(command[0], command, CHILD_ENV)
                    except BaseException as error:
                        try:
                            os.write(
                                2,
                                ("gated benchmark launch failed: %s\n" %
                                 error).encode("utf-8", errors="replace"))
                        finally:
                            os._exit(126)
                handle = ForkedProcessHandle(pid)
                containment.spawned_process = handle
                containment.observe_spawn(handle)
                containment.attach(pid)
                os.close(stdout_write)
                stdout_write = -1
                os.close(stderr_write)
                stderr_write = -1
                for descriptor in (stdout_read, stderr_read):
                    os.set_blocking(descriptor, False)
                    selector.register(descriptor, selectors.EVENT_READ)
                retained = {
                    stdout_read: {
                        "output": stdout_fd,
                        "limit": MAX_BENCHMARK_STDOUT_BYTES,
                        "count": 0,
                        "label": "stdout",
                    },
                    stderr_read: {
                        "output": stderr_fd,
                        "limit": MAX_BENCHMARK_STDERR_BYTES,
                        "count": 0,
                        "label": "stderr",
                    },
                }

                stopped_pid, stopped_status, gate_rusage = wait4_until(
                    pid, os.WUNTRACED,
                    time.monotonic() + min(timeout, 5.0))
                require(not (os.WIFEXITED(stopped_status) or
                             os.WIFSIGNALED(stopped_status)) and
                        os.WIFSTOPPED(stopped_status) and
                        os.WSTOPSIG(stopped_status) == signal.SIGSTOP,
                        "benchmark child failed before the timing gate")
                require(stopped_pid == pid and os.getsid(pid) == pid and
                        set(os.sched_getaffinity(pid)) == {cpu},
                        "benchmark child gate identity/session/affinity is "
                        "invalid")
                child_start_time, child_inode = current_task_identity(pid, pid)
                gate_state = current_process_state(pid)
                require(gate_state["state"] in ("T", "t") and
                        gate_state["pgrp"] == pid and
                        gate_state["session"] == pid and
                        gate_state["processor"] == cpu,
                        "benchmark child stopped-state identity is invalid")
                gate_cpu_time_ns = rusage_cpu_time_ns(gate_rusage)
                child_sched_before = task_scheduler_runtime_ns(pid)

                target_stat_before = cpu_ticks(cpu)
                target_sched_before = cpu_scheduler_runtime_ns(cpu)
                topology = Path(
                    "/sys/devices/system/cpu/cpu%d/topology/"
                    "thread_siblings_list" % cpu)
                sibling_set = parse_cpu_list(topology.read_text()) - {cpu}
                require(sibling_set == {sibling},
                        "gated child CPU has a different SMT sibling")
                sibling_stat_before = cpu_ticks(sibling)
                sibling_sched_before = cpu_scheduler_runtime_ns(sibling)
                started_ns = time.time_ns()
                os.kill(pid, signal.SIGCONT)
                deadline = time.monotonic() + timeout
                exit_observed = False
                output_failure = None
                while not exit_observed or selector.get_map():
                    observed = os.waitid(
                        os.P_PID, pid,
                        os.WEXITED | os.WNOHANG | os.WNOWAIT)
                    if observed is not None:
                        require(observed.si_pid == pid,
                                "waitid observed an unexpected benchmark child")
                        exit_observed = True
                    remaining = deadline - time.monotonic()
                    if remaining <= 0:
                        output_failure = (
                            "benchmark exceeded %.3f seconds" % timeout)
                        break
                    events = selector.select(min(remaining, 0.05)) \
                        if selector.get_map() else ()
                    if not events and not selector.get_map():
                        time.sleep(min(0.001, remaining))
                    for key, unused_mask in events:
                        del unused_mask
                        descriptor = key.fd
                        try:
                            block = os.read(descriptor, 65536)
                        except BlockingIOError:
                            continue
                        if not block:
                            selector.unregister(descriptor)
                            continue
                        stream = retained[descriptor]
                        available = stream["limit"] - stream["count"]
                        prefix = block[:available]
                        if prefix:
                            _write_all(stream["output"], prefix)
                            stream["count"] += len(prefix)
                        if len(block) > available:
                            output_failure = (
                                "benchmark %s exceeded retained byte limit" %
                                stream["label"])
                            break
                    if output_failure is not None:
                        break
                ended_ns = time.time_ns()
                if output_failure is not None:
                    raise EvidenceError(output_failure)

                zombie_start_time, zombie_inode = current_task_identity(
                    pid, pid)
                zombie_state = current_process_state(pid)
                require(zombie_start_time == child_start_time and
                        zombie_inode == child_inode and
                        zombie_state["state"] == "Z" and
                        zombie_state["pgrp"] == pid and
                        zombie_state["session"] == pid and
                        zombie_state["processor"] == cpu and
                        set(os.sched_getaffinity(pid)) == {cpu},
                        "benchmark zombie identity/session/affinity is invalid")
                retained_session_members = \
                    process_group_or_session_members(pid)
                require(retained_session_members == [{
                            "pid": pid, "state": "Z", "pgrp": pid,
                            "session": pid,
                            "start_time_ticks": child_start_time,
                        }],
                        "benchmark session retained unexpected descendants: %s"
                        % retained_session_members)
                child_sched_after = task_scheduler_runtime_ns(pid)
                target_sched_after = cpu_scheduler_runtime_ns(cpu)
                target_stat_after = cpu_ticks(cpu)
                sibling_sched_after = cpu_scheduler_runtime_ns(sibling)
                sibling_stat_after = cpu_ticks(sibling)
                exited_pid, status, exit_rusage = os.wait4(pid, 0)
                require(exited_pid == pid,
                        "wait4 reaped an unexpected benchmark child")
                handle.mark_reaped(status)
                require_child_identity_gone_after_reap(pid)
                child_runtime_ns = child_sched_after - child_sched_before
                require(child_runtime_ns >= 0,
                        "benchmark child schedstat runtime moved backwards")
                wait4_runtime_ns = child_cpu_time_ns(
                    gate_rusage, exit_rusage)
                rusage_difference_ns = abs(
                    wait4_runtime_ns - child_runtime_ns)
                result = {
                    "return_code": handle.returncode,
                    "timed_out": False,
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
                        target_sched_before, target_sched_after,
                        child_runtime_ns, TARGET_RUNTIME_TOLERANCE_NS),
                    "wait4_crosscheck": {
                        "child_cpu_time_ns": wait4_runtime_ns,
                        "child_schedstat_delta_ns": child_runtime_ns,
                        "absolute_difference_ns": rusage_difference_ns,
                        "tolerance_ns": RUSAGE_CROSSCHECK_TOLERANCE_NS,
                        "accepted": rusage_difference_ns <=
                            RUSAGE_CROSSCHECK_TOLERANCE_NS,
                    },
                    "target_interrupts": target_interrupt_evidence(
                        target_stat_before, target_stat_after),
                    "sibling_runtime": {
                        "scheduler_before_ns": sibling_sched_before,
                        "scheduler_after_ns": sibling_sched_after,
                        "scheduler_delta_ns":
                            sibling_sched_after - sibling_sched_before,
                        "accepted":
                            sibling_sched_after == sibling_sched_before,
                    },
                    "sibling_interrupts": target_interrupt_evidence(
                        sibling_stat_before, sibling_stat_after),
                }
                containment.terminate_and_reap(handle)
        except PROCESS_SUPPORT.EvidenceError as error:
            raise EvidenceError(
                "benchmark descendant containment failed: %s" % error) \
                from error
    finally:
        primary = sys.exc_info()[1]
        cleanup_failures = []
        try:
            selector.close()
        except BaseException as error:
            cleanup_failures.append("selector: %s" % error)
        for descriptor in (
                stdout_read, stderr_read, stdout_write, stderr_write):
            if descriptor < 0:
                continue
            try:
                os.close(descriptor)
            except OSError as error:
                cleanup_failures.append(
                    "pipe fd %d: %s" % (descriptor, error))
        for descriptor in (stdout_fd, stderr_fd):
            if descriptor < 0:
                continue
            try:
                os.fsync(descriptor)
            except OSError as error:
                cleanup_failures.append(
                    "output fsync fd %d: %s" % (descriptor, error))
            try:
                os.close(descriptor)
            except OSError as error:
                cleanup_failures.append(
                    "output close fd %d: %s" % (descriptor, error))
        try:
            signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)
        except BaseException as error:
            cleanup_failures.append("signal mask: %s" % error)
        if cleanup_failures:
            message = "gated benchmark cleanup failed: %s" % \
                "; ".join(cleanup_failures)
            if primary is not None:
                message += "; primary failure: %s: %s" % (
                    type(primary).__name__, primary)
            raise EvidenceError(message) from primary
    if executable_descriptor is not None:
        validate_execution_descriptor(
            executable_descriptor, execution_identity)
    require(result is not None, "gated benchmark produced no result")
    return result


def invocation_isolation_accepted(envelope: dict[str, Any]) -> bool:
    return envelope["reserved_sibling_nonidle_jiffies"] == 0 and \
        envelope["target_runtime"]["accepted"] and \
        envelope["target_interrupts"]["accepted"] and \
        envelope["sibling_runtime"]["accepted"] and \
        envelope["sibling_interrupts"]["accepted"] and \
        envelope["wait4_crosscheck"]["accepted"]


def run_invocation(binary: Path, implementation: str, cell: dict[str, Any],
                   mode: str,
                   cpu: int, sibling: int, timeout: float,
                   raw_directory: Path, round_index: int,
                   slot_index: int, max_retries: int,
                   reservation_path: Path,
                   reservation: dict[str, Any],
                   isolation_epoch_digest: str,
                   child_signal_mask: set[signal.Signals],
                   campaign_lock_descriptor: int,
                   execution_identity: dict[str, Any] | None = None,
                   executable_descriptor: int | None = None,
                   ) -> dict[str, Any]:
    require((executable_descriptor is None) ==
            (execution_identity is None),
            "invocation immutable descriptor and identity must be supplied "
            "together")
    command = benchmark_command(binary, cell, cpu)
    if executable_descriptor is not None:
        validate_execution_descriptor(
            executable_descriptor, execution_identity)
        require(binary == Path(execution_identity["command_path"]),
                "invocation command path is not immutable")
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
            child_signal_mask, executable_descriptor, execution_identity,
            campaign_lock_descriptor)
        if executable_descriptor is not None:
            validate_execution_descriptor(
                executable_descriptor, execution_identity)
        require_no_pending_control_signal()
        stdout_identity, stdout = stable_file_snapshot(
            stdout_path, "benchmark stdout", MAX_BENCHMARK_STDOUT_BYTES,
            require_canonical=True)
        stderr_identity, stderr = stable_file_snapshot(
            stderr_path, "benchmark stderr", MAX_BENCHMARK_STDERR_BYTES,
            require_canonical=True)
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
            "execution_identity": execution_identity,
            "started_ns": gated["started_ns"],
            "ended_ns": gated["ended_ns"],
            "return_code": gated["return_code"],
            "stdout": stdout_identity,
            "stderr": stderr_identity,
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
        isolation_accepted = invocation_isolation_accepted(envelope)
        if gated["return_code"] == 0 and isolation_accepted:
            try:
                parsed = strict_json_bytes(stdout, "benchmark stdout")
                normalized = validate_result(parsed, cell, mode)
                envelope["result"] = parsed
                envelope["normalized"] = normalized
                envelope["accepted"] = True
            except EvidenceError as error:
                envelope["validation_error"] = str(error)
        envelope_path = raw_directory / (stem + ".envelope.json")
        envelope["envelope_path"] = str(envelope_path)
        atomic_json(envelope_path, envelope, MAX_ENVELOPE_JSON_BYTES)
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
    binary_identity_digest = binary_identity.get("digest")
    require(isinstance(binary_identity_digest, str) and
            re.fullmatch(r"[0-9a-f]{64}", binary_identity_digest) is not None,
            "cell has no valid binary-identity digest")
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
    first_use_ratio = ratio_summary(
        invocations, rounds,
        lambda value: value["median_us"] + value["plan_setup_us"])
    ratios = {
        "execution": execution_ratio,
        "first_use_plan_plus_execution": first_use_ratio,
    }
    for reuse in (1, 8, 64):
        ratios["amortized_reuse_%d" % reuse] = ratio_summary(
            invocations, rounds,
            lambda value, reuse_count=reuse: value["median_us"] +
                value["plan_setup_us"] / reuse_count)
    result = {
        "schema": CELL_SCHEMA,
        "cell": cell,
        "binary_identity_digest": binary_identity_digest,
        "orientation": "baseline_time_over_candidate_time",
        "round_ratios": execution_ratio["round_ratios"],
        "geometric_mean_ratio": execution_ratio["geometric_mean_ratio"],
        "ci95": execution_ratio["ci95"],
        "metric_ratios": ratios,
        "plan_setup_us": {
            implementation: {
                "samples": [
                    item["normalized"]["plan_setup_us"]
                    for item in invocations
                    if item["implementation"] == implementation
                ],
                "median": statistics.median(
                    item["normalized"]["plan_setup_us"]
                    for item in invocations
                    if item["implementation"] == implementation),
            } for implementation in ("baseline", "candidate")
        },
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
                "execution_identity": item["execution_identity"],
                "median_us": item["normalized"]["median_us"],
                "plan_setup_us": item["normalized"]["plan_setup_us"],
                "envelope_path": item["envelope_path"],
                "envelope_sha256": sha256_file(
                    Path(item["envelope_path"]), MAX_ENVELOPE_JSON_BYTES),
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


def validate_current_binary_identity(
        identity: dict[str, Any], *,
        inherited_descriptors: Sequence[int] = ()) -> None:
    validate_binary_identity_structure(identity)
    recorded_source = identity["source"]
    require(source_identity(
                Path(recorded_source["root"]),
                inherited_descriptors=inherited_descriptors) ==
            recorded_source,
            "current source snapshot changed")
    for key in (
            "runner", "process_containment_source", "pair_lease_source"):
        record = identity[key]
        require(file_identity(Path(record["path"])) == record,
                "current %s identity changed" % key)
    modes = identity["comparison_modes"]
    for label in ("baseline", "candidate"):
        executable = identity[label]["executable"]
        for role in ("origin", "frozen"):
            record = executable[role]
            require(file_identity(Path(record["path"])) == record,
                    "current %s %s executable changed" % (label, role))
        current_build = capture_current_build(
            Path(executable["origin"]["path"]),
            Path(recorded_source["root"]), modes[label],
            inherited_descriptors=inherited_descriptors,
            retained_reproducible_build=
                identity[label]["build"]["reproducible_build"])
        require(current_build == identity[label]["build"],
                "current %s build provenance changed" % label)
    require(compare_current_builds(
                identity["baseline"]["build"],
                identity["candidate"]["build"], modes) ==
            identity["build_comparison"],
            "current cross-build comparison changed")

def run_campaign(options: argparse.Namespace) -> int:
    require(type(options.lock_timeout) in (int, float) and
            math.isfinite(options.lock_timeout) and
            options.lock_timeout >= 0,
            "lock timeout must be finite and nonnegative")
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
            # Never issue LOCK_UN on this shared open-file description.
            # A timed child may retain a duplicate after coordinator failure;
            # close-only release keeps the flock held until the last survivor
            # exits.
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
    require(type(options.timeout) in (int, float) and
            math.isfinite(options.timeout) and options.timeout > 0 and
            type(options.lock_timeout) in (int, float) and
            math.isfinite(options.lock_timeout) and
            options.lock_timeout >= 0,
            "benchmark timeout must be finite and positive and lock timeout "
            "finite and nonnegative")
    modes = comparison_modes(options.baseline_mode, options.candidate_mode)
    lock_descriptor = lock_stream.fileno()
    inherited_descriptors = (lock_descriptor,)
    execution_nonce = secrets.token_hex(16)
    topology_identity = require_smt_pair(
        options.cpu, options.reserved_sibling)
    source_root = options.source_root.resolve(strict=True)
    validate_default_mode_source(source_root)
    PairLease, pair_lease_source = load_pair_lease(source_root)
    reservation = reservation_identity(
        options.reservation_file, options.cpu, options.reserved_sibling)
    baseline_origin = options.baseline.resolve(strict=True)
    candidate_origin = options.candidate.resolve(strict=True)
    require(os.access(baseline_origin, os.X_OK) and
            os.access(candidate_origin, os.X_OK),
            "both benchmark paths must be executable")
    matrix = matrix_for_preset(options.preset)
    selected = matrix["cells"]
    if options.cell:
        identifiers = set(options.cell)
        selected = [cell for cell in selected if cell["id"] in identifiers]
        require(len(selected) == len(identifiers), "unknown --cell identifier")

    output, output_parent_identity = validate_output_root_candidate(
        options.output, source_root)

    source_record = source_identity(
        source_root, inherited_descriptors=inherited_descriptors)
    require(not source_record["status_short"],
            "authoritative campaign requires a clean source worktree")
    baseline_build = capture_current_build(
        baseline_origin, source_root, modes["baseline"],
        inherited_descriptors=inherited_descriptors)
    candidate_build = capture_current_build(
        candidate_origin, source_root, modes["candidate"],
        inherited_descriptors=inherited_descriptors)
    build_comparison = compare_current_builds(
        baseline_build, candidate_build, modes)
    output = create_output_root(output, output_parent_identity)
    artifacts = output / "artifacts"
    artifacts.mkdir(mode=0o700)
    os.chmod(artifacts, 0o700)
    baseline_artifact = freeze_executable(
        baseline_origin, artifacts, "baseline")
    candidate_artifact = freeze_executable(
        candidate_origin, artifacts, "candidate")
    baseline = Path(baseline_artifact["frozen"]["path"])
    candidate = Path(candidate_artifact["frozen"]["path"])
    binary_identity = {
        "schema": BINARY_IDENTITY_SCHEMA,
        "comparison_modes": modes,
        "runner": file_identity(Path(__file__).resolve(strict=True)),
        "process_containment_source": file_identity(PROCESS_SUPPORT_PATH),
        "pair_lease_source": pair_lease_source,
        "source": source_record,
        "baseline": {
            "build": baseline_build,
            "executable": baseline_artifact,
        },
        "candidate": {
            "build": candidate_build,
            "executable": candidate_artifact,
        },
        "build_comparison": build_comparison,
    }
    binary_identity["digest"] = sha256_bytes(
        canonical_bytes(binary_identity))
    validate_binary_identity_structure(binary_identity)
    execution_guard = ImmutableFrozenExecutables(output, binary_identity)
    execution_guard_active = False
    try:
        execution_identity, execution_descriptors = \
            execution_guard.__enter__()
        execution_guard_active = True
        request = {
            "schema": SCHEMA,
            "matrix_sha256": matrix["matrix_sha256"],
            "matrix_preset": options.preset,
            "selected_cell_ids": [cell["id"] for cell in selected],
            "rounds": options.rounds,
            "order": ORDER,
            "binary_identity": binary_identity,
            "execution_identity": execution_identity,
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
        request_path = output / "request.json"
        atomic_json(request_path, request, MAX_REQUEST_JSON_BYTES)
        atomic_json(output / "matrix.json", matrix, MAX_MATRIX_JSON_BYTES)
        (output / "cells").mkdir(mode=0o700)
        (output / "raw").mkdir(mode=0o700)
    except Exception:
        if execution_guard_active:
            execution_guard.__exit__(None, None, None)
        raise
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
            raw_directory.mkdir(mode=0o700)
            for round_index in range(options.rounds):
                for slot_index, implementation in enumerate(ORDER):
                    execution = execution_identity[implementation]
                    binary = Path(execution["command_path"])
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
                        isolation_epoch_digest, child_signal_mask,
                        lock_descriptor,
                        execution, execution_descriptors[implementation]))
                    require_no_pending_control_signal()
                    require(pair_guard.revalidate() == pair_lease_identity,
                            "canonical CPU-pair lease changed during invocation")
            summary = summarize_cell(
                cell, invocations, options.rounds, binary_identity)
            atomic_json(cell_path, summary, MAX_CELL_JSON_BYTES)
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
        validate_current_binary_identity(
            binary_identity,
            inherited_descriptors=inherited_descriptors)
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
        validate_campaign_inventory(
            output, [cell["id"] for cell in selected],
            manifest_present=False)
        isolation_epoch = {
            "start": isolation_epoch_start,
            "start_digest": isolation_epoch_digest,
            "restoration": affinity_restoration,
        }
        manifest = {
            "schema": SCHEMA,
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
            "promotion_analysis":
                analyze_same_source_promotion(matrix, summaries, modes),
        }
        manifest["digest"] = sha256_bytes(canonical_bytes(manifest))
        validate_current_binary_identity(
            binary_identity,
            inherited_descriptors=inherited_descriptors)
        require(global_lock_identity(lock_stream, options.lock) ==
                lock_identity,
                "canonical GF8 lock identity changed at manifest commit")
        require(reservation_identity(
                    options.reservation_file, options.cpu,
                    options.reserved_sibling) == reservation,
                "CPU reservation changed at manifest commit")
        require(pair_guard.revalidate() == pair_lease_identity,
                "canonical CPU-pair lease changed at manifest commit")
        atomic_json(
            output / "manifest.json", manifest, MAX_MANIFEST_JSON_BYTES)
        require(global_lock_identity(lock_stream, options.lock) ==
                lock_identity,
                "canonical GF8 lock identity changed after manifest commit")
        validate_campaign_inventory(
            output, [cell["id"] for cell in selected])
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
        if execution_guard_active:
            try:
                execution_guard.__exit__(None, None, None)
                execution_guard_active = False
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
    validate_reservation_payload(reservation, cpu, sibling)


def validate_interrupt_record(
        value: Any, label: str, require_accepted: bool = True) -> None:
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
    require(value == expected and
            (not require_accepted or value["accepted"] is True),
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
        envelope: dict[str, Any], child_runtime_ns: int,
        require_accepted: bool = True) -> None:
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
    require(target == expected_target and
            (not require_accepted or target["accepted"] is True),
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
            } and
            (not require_accepted or sibling["accepted"] is True),
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
            } and
            (not require_accepted or crosscheck["accepted"] is True),
            "wait4 cross-check arithmetic or acceptance is invalid")


def validate_retained_stream(
        value: Any, expected_path: Path, label: str,
        maximum_bytes: int) -> bytes:
    validate_file_identity_record(value, label)
    identity, payload = stable_file_snapshot(
        expected_path, label, maximum_bytes, require_canonical=True)
    require(Path(value["path"]) == expected_path and identity == value,
            "%s retained identity changed" % label)
    return payload


def retained_attempt_paths(
        accepted_path: Path, cell: dict[str, Any], implementation: str,
        round_index: int, slot_index: int, max_retries: int,
        campaign_root: Path) -> list[tuple[Path, Path, Path]]:
    expected_directory = campaign_root / "raw" / cell["id"]
    prefix = "round%02d-slot%d-%s-attempt" % (
        round_index, slot_index, implementation)
    require(accepted_path.parent == expected_directory and
            accepted_path.is_file() and
            accepted_path.resolve(strict=True) == accepted_path,
            "accepted invocation envelope escaped its cell directory")
    match = re.fullmatch(
        re.escape(prefix) + r"([0-9]{2})\.envelope\.json",
        accepted_path.name)
    require(match is not None,
            "accepted invocation envelope has the wrong ABBA slot")
    accepted_attempt = int(match.group(1))
    require(accepted_attempt <= max_retries,
            "accepted invocation exceeds the retry limit")
    attempts = []
    expected_names = set()
    for attempt in range(accepted_attempt + 1):
        stem = prefix + ("%02d" % attempt)
        envelope = expected_directory / (stem + ".envelope.json")
        stdout = expected_directory / (stem + ".stdout.json")
        stderr = expected_directory / (stem + ".stderr.txt")
        attempts.append((envelope, stdout, stderr))
        expected_names.update((envelope.name, stdout.name, stderr.name))
    require(expected_directory.is_dir(),
            "retained invocation cell directory is missing")
    actual_names = {
        path.name for path in expected_directory.iterdir()
        if path.name.startswith(prefix)
    }
    require(actual_names == expected_names,
            "retained invocation retry chain has extra or missing attempts")
    return attempts


def validate_retry_envelope(
        value: Any, envelope_path: Path, cell: dict[str, Any],
        implementation: str, round_index: int, slot_index: int,
        attempt_index: int, request: dict[str, Any],
        isolation_epoch_digest: str, campaign_root: Path) -> dict[str, Any]:
    require_exact_keys(value, {
        "schema", "implementation", "command", "execution_identity",
        "started_ns", "ended_ns",
        "return_code", "stdout", "stderr", "reserved_sibling_before",
        "reserved_sibling_after", "reserved_sibling_nonidle_jiffies",
        "reservation_before", "reservation_after",
        "same_user_pair_affinity_before", "same_user_pair_affinity_after",
        "gate", "target_runtime", "target_interrupts", "sibling_runtime",
        "sibling_interrupts", "wait4_crosscheck",
        "isolation_epoch_digest", "accepted", "envelope_path",
    }, "retained rejected retry envelope")
    expected_directory = campaign_root / "raw" / cell["id"]
    expected_stem = "round%02d-slot%d-%s-attempt%02d" % (
        round_index, slot_index, implementation, attempt_index)
    require(envelope_path ==
                expected_directory / (expected_stem + ".envelope.json") and
            envelope_path.resolve(strict=True) == envelope_path and
            value["schema"] == ENVELOPE_SCHEMA and
            value["implementation"] == implementation and
            type(value["return_code"]) is int and
            value["accepted"] is False and
            value["envelope_path"] == str(envelope_path),
            "retained rejected retry header is invalid")
    started = require_int(value["started_ns"], "retry start", 1)
    ended = require_int(value["ended_ns"], "retry end", 1)
    require(ended >= started, "retry timing moved backwards")
    expected_execution = request["execution_identity"][implementation]
    expected_binary = expected_execution["command_path"]
    require(value["command"] == benchmark_command(
                Path(expected_binary), cell, request["cpu"]),
            "retained retry command does not match the frozen cell")
    require(value["execution_identity"] == expected_execution,
            "retained retry used a different immutable executable")
    require(value["reservation_before"] == request["reservation"] and
            value["reservation_after"] == request["reservation"] and
            value["same_user_pair_affinity_before"] == [] and
            value["same_user_pair_affinity_after"] == [] and
            value["isolation_epoch_digest"] == isolation_epoch_digest,
            "retained retry is outside the recorded isolation epoch")
    validate_retained_stream(
        value["stdout"], expected_directory / (expected_stem + ".stdout.json"),
        "retry stdout", MAX_BENCHMARK_STDOUT_BYTES)
    validate_retained_stream(
        value["stderr"], expected_directory / (expected_stem + ".stderr.txt"),
        "retry stderr", MAX_BENCHMARK_STDERR_BYTES)
    child_runtime = validate_gate_record(value["gate"], request["cpu"])
    validate_runtime_records(
        value, child_runtime, require_accepted=False)
    validate_interrupt_record(
        value["target_interrupts"], "retry target interrupts",
        require_accepted=False)
    validate_interrupt_record(
        value["sibling_interrupts"], "retry sibling interrupts",
        require_accepted=False)
    require(value["reserved_sibling_before"] ==
                value["sibling_interrupts"]["before"] and
            value["reserved_sibling_after"] ==
                value["sibling_interrupts"]["after"] and
            value["reserved_sibling_nonidle_jiffies"] == nonidle_delta(
                tuple(value["reserved_sibling_before"]),
                tuple(value["reserved_sibling_after"])),
            "rejected retry sibling evidence is invalid")
    require(not invocation_isolation_accepted(value),
            "rejected retry has no declared isolation failure")
    return value


def validate_envelope(
        value: Any, envelope_path: Path, cell: dict[str, Any],
        implementation: str, round_index: int, slot_index: int,
        request: dict[str, Any], isolation_epoch_digest: str,
        campaign_root: Path, attempt_index: int | None = None
        ) -> dict[str, Any]:
    require_exact_keys(value, {
        "schema", "implementation", "command", "execution_identity",
        "started_ns", "ended_ns",
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
            type(value["return_code"]) is int and
            value["return_code"] == 0 and value["accepted"] is True and
            value["envelope_path"] == str(envelope_path),
            "retained invocation envelope header is invalid")
    started = require_int(value["started_ns"], "invocation start", 1)
    ended = require_int(value["ended_ns"], "invocation end", 1)
    require(ended >= started, "invocation timing moved backwards")
    expected_execution = request["execution_identity"][implementation]
    expected_binary = expected_execution["command_path"]
    require(value["command"] == benchmark_command(
                Path(expected_binary), cell, request["cpu"]),
            "retained invocation command does not match the frozen cell")
    require(value["execution_identity"] == expected_execution,
            "retained invocation used a different immutable executable")
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
            int(match.group(1)) <= request["max_retries"] and
            (attempt_index is None or
             int(match.group(1)) == attempt_index),
            "retained invocation envelope has the wrong ABBA slot/attempt")
    stem = envelope_path.name[:-len(".envelope.json")]
    stdout_path = expected_directory / (stem + ".stdout.json")
    stderr_path = expected_directory / (stem + ".stderr.txt")
    stdout = validate_retained_stream(
        value["stdout"], stdout_path, "stdout", MAX_BENCHMARK_STDOUT_BYTES)
    validate_retained_stream(
        value["stderr"], stderr_path, "stderr", MAX_BENCHMARK_STDERR_BYTES)
    parsed = strict_json_bytes(stdout, "retained stdout")
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


def validate_source_identity_structure(source: Any) -> None:
    require_exact_keys(source, {
        "root", "head", "head_tree", "status_short", "files",
        "snapshot_sha256",
    }, "source identity")
    root_text = source["root"]
    root = Path(root_text) if isinstance(root_text, str) else Path(".")
    require(isinstance(root_text, str) and root.is_absolute() and
            str(root) == root_text and
            os.path.normpath(root_text) == root_text and
            all(part not in (".", "..") for part in root.parts) and
            isinstance(source["head"], str) and
            re.fullmatch(r"[0-9a-f]{40}", source["head"]) is not None and
            isinstance(source["head_tree"], str) and
            re.fullmatch(r"[0-9a-f]{40}",
                         source["head_tree"]) is not None and
            source["status_short"] == [],
            "source identity header is invalid")
    require(isinstance(source["files"], dict) and
            set(source["files"]) == set(SOURCE_FILES),
            "source identity file set is invalid")
    for relative, record in source["files"].items():
        path = Path(relative)
        require(not path.is_absolute() and str(path) == relative and
                all(part not in ("", ".", "..") for part in path.parts),
                "source file key is not a safe canonical relative path")
        require_exact_keys(record, {"size", "sha256"},
                           "source file %s" % relative)
        require(type(record["size"]) is int and record["size"] >= 0 and
                isinstance(record["sha256"], str) and
                re.fullmatch(r"[0-9a-f]{64}", record["sha256"]),
                "source file identity is invalid")
    source_copy = dict(source)
    source_digest = source_copy.pop("snapshot_sha256")
    require(source_digest == sha256_bytes(canonical_bytes(source_copy)),
            "source snapshot digest is invalid")


def validate_binary_identity_structure(value: Any) -> None:
    require_exact_keys(value, {
        "schema", "comparison_modes", "runner",
        "process_containment_source", "pair_lease_source", "source",
        "baseline", "candidate", "build_comparison", "digest",
    }, "binary identity")
    require(value["schema"] in (
                BINARY_IDENTITY_SCHEMA_V4, BINARY_IDENTITY_SCHEMA),
            "binary identity has the wrong schema")
    copy = dict(value)
    digest = copy.pop("digest")
    require(isinstance(digest, str) and
            digest == sha256_bytes(canonical_bytes(copy)),
            "binary identity digest is invalid")
    modes = value["comparison_modes"]
    require(isinstance(modes, dict) and
            set(modes) == {"baseline", "candidate"},
            "binary identity comparison modes are invalid")
    comparison_modes(modes["baseline"], modes["candidate"])
    expected_build_schema = (
        BUILD_SCHEMA if value["schema"] == BINARY_IDENTITY_SCHEMA
        else BUILD_SCHEMA_V4)
    validate_file_identity_record(value["runner"], "runner")
    validate_file_identity_record(
        value["process_containment_source"],
        "process containment source")
    validate_file_identity_record(
        value["pair_lease_source"], "pair lease source")
    validate_source_identity_structure(value["source"])
    support = value["process_containment_source"]
    source_support = value["source"]["files"][PROCESS_SUPPORT_RELATIVE]
    require(support["size"] == source_support["size"] and
            support["sha256"] == source_support["sha256"],
            "loaded process containment source differs from the attested "
            "source snapshot")

    for label in ("baseline", "candidate"):
        role = value[label]
        require_exact_keys(role, {"build", "executable"},
                           "%s binary role" % label)
        executable = role["executable"]
        require_exact_keys(executable, {"schema", "origin", "frozen"},
                           "%s frozen executable" % label)
        require(executable["schema"] ==
                    "leopard2-small-direct-frozen-executable/v1",
                "%s frozen executable schema is invalid" % label)
        validate_file_identity_record(
            executable["origin"], "%s origin executable" % label)
        validate_file_identity_record(
            executable["frozen"], "%s frozen executable" % label)
        require(executable["origin"]["sha256"] ==
                    executable["frozen"]["sha256"] and
                executable["origin"]["size"] ==
                    executable["frozen"]["size"] and
                executable["origin"]["path"] !=
                    executable["frozen"]["path"],
                "%s frozen executable is not an immutable copy" % label)

        build = role["build"]
        require_exact_keys(build, {
            "schema", "mode", "cache_mode", "closure",
            "reproducible_build",
            "effective_configuration", "attestation_module",
            "attestation_generator", "provenance_module",
            "configuration_reader", "digest",
        }, "%s build identity" % label)
        require(build["schema"] == expected_build_schema and
                build["mode"] == modes[label] and
                build["cache_mode"] == MODE_CACHE_VALUES[modes[label]],
                "%s build mode identity is invalid" % label)
        build_copy = dict(build)
        build_digest = build_copy.pop("digest")
        require(build_digest == sha256_bytes(canonical_bytes(build_copy)),
                "%s build digest is invalid" % label)
        require(isinstance(build["closure"], dict) and
                build["closure"].get("schema") ==
                    "leopard2-production-build-closure/v1" and
                isinstance(build["effective_configuration"], dict),
                "%s current provenance structure is invalid" % label)
        validate_effective_configuration_attestation(
            build["effective_configuration"].get(
                "effective_configuration_attestation"),
            build["mode"], expected_build_schema,
            "%s retained effective configuration" % label)
        validate_reproducible_build_proof(
            build["reproducible_build"], build["closure"], label)
        for key in ("attestation_module", "attestation_generator",
                    "provenance_module", "configuration_reader"):
            validate_file_identity_record(
                build[key], "%s %s" % (label, key))
        require(build["closure"].get("executable", {}).get("sha256") ==
                    executable["origin"]["sha256"] and
                build["closure"].get("executable", {}).get("size") ==
                    executable["origin"]["size"],
                "%s build closure is not bound to origin binary" % label)

    expected_comparison = compare_current_builds(
        value["baseline"]["build"], value["candidate"]["build"], modes)
    require(value["build_comparison"] == expected_comparison,
            "retained cross-build comparison is invalid")
    require(value["baseline"]["executable"]["frozen"]["sha256"] !=
            value["candidate"]["executable"]["frozen"]["sha256"],
            "frozen comparison executables are byte-identical")


def validate_retained_frozen_executables(
        identity: dict[str, Any], campaign_root: Path) -> None:
    for label in ("baseline", "candidate"):
        expected = campaign_root / "artifacts" / label / "bench_leopard2"
        record = identity[label]["executable"]["frozen"]
        require(Path(record["path"]) == expected,
                "%s frozen executable escaped its artifact lane" % label)
        require(expected.is_file() and
                expected.resolve(strict=True) == expected,
                "%s retained frozen executable changed" % label)
        metadata = expected.stat()
        require(
                stat.S_ISREG(metadata.st_mode) and
                metadata.st_uid == os.getuid() and
                metadata.st_nlink == 1 and
                stat.S_IMODE(metadata.st_mode) == 0o555 and
                file_identity(expected, MAX_EXECUTABLE_BYTES) == record,
                "%s retained frozen executable changed" % label)


def validate_campaign_inventory(
        campaign_root: Path, selected_cell_ids: Sequence[str],
        *, manifest_present: bool = True) -> None:
    root_metadata = os.stat(campaign_root, follow_symlinks=False)
    require(campaign_root.is_dir() and
            campaign_root.resolve(strict=True) == campaign_root and
            stat.S_ISDIR(root_metadata.st_mode) and
            root_metadata.st_uid == os.getuid() and
            stat.S_IMODE(root_metadata.st_mode) == 0o700,
            "campaign root is not a canonical owned 0700 directory")
    expected_top = {
        "request.json", "matrix.json",
        "artifacts", "cells", "raw",
    }
    if manifest_present:
        expected_top.add("manifest.json")
    require({path.name for path in campaign_root.iterdir()} == expected_top,
            "campaign root contains extra or missing artifacts")
    for name in ("artifacts", "cells", "raw"):
        path = campaign_root / name
        metadata = os.stat(path, follow_symlinks=False)
        require(path.is_dir() and path.resolve(strict=True) == path and
                stat.S_ISDIR(metadata.st_mode) and
                metadata.st_uid == os.getuid() and
                stat.S_IMODE(metadata.st_mode) == 0o700,
                "campaign %s directory is not canonical and owned" % name)
    artifacts = campaign_root / "artifacts"
    require({path.name for path in artifacts.iterdir()} ==
                {"baseline", "candidate"},
            "campaign artifact lanes are incomplete or contain extras")
    for label in ("baseline", "candidate"):
        lane = artifacts / label
        metadata = os.stat(lane, follow_symlinks=False)
        require(lane.is_dir() and lane.resolve(strict=True) == lane and
                stat.S_ISDIR(metadata.st_mode) and
                metadata.st_uid == os.getuid() and
                stat.S_IMODE(metadata.st_mode) == 0o700 and
                {path.name for path in lane.iterdir()} == {"bench_leopard2"},
                "%s artifact lane contains extra or missing files" % label)
    expected_cells = {
        identifier + ".json" for identifier in selected_cell_ids
    }
    require({path.name for path in (campaign_root / "cells").iterdir()} ==
                expected_cells,
            "campaign cells directory contains extra or missing summaries")
    raw = campaign_root / "raw"
    require({path.name for path in raw.iterdir()} ==
                set(selected_cell_ids),
            "campaign raw directory contains extra or missing cell lanes")
    for identifier in selected_cell_ids:
        lane = raw / identifier
        metadata = os.stat(lane, follow_symlinks=False)
        require(lane.is_dir() and lane.resolve(strict=True) == lane and
                stat.S_ISDIR(metadata.st_mode) and
                metadata.st_uid == os.getuid() and
                stat.S_IMODE(metadata.st_mode) == 0o700,
                "campaign raw cell lane is not canonical and owned")


def verify_campaign(options: argparse.Namespace) -> int:
    manifest_path = options.manifest.resolve(strict=True)
    campaign_root = manifest_path.parent
    require(manifest_path == campaign_root / "manifest.json",
            "campaign manifest must use its canonical filename")
    manifest = strict_json_file(
        manifest_path, "campaign manifest", MAX_MANIFEST_JSON_BYTES)
    require_exact_keys(manifest, {
        "schema", "created_utc", "request", "global_lock", "pair_lease",
        "affinity_exclusion", "isolation_epoch", "host", "statistics",
        "cell_summaries", "promotion_analysis", "digest",
    }, "campaign manifest")
    require(manifest["schema"] == SCHEMA, "wrong manifest schema")
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
        "rounds", "order", "binary_identity", "execution_identity",
        "cpu", "reserved_sibling",
        "reservation", "topology_identity", "isolation_policy",
        "timeout_seconds", "max_retries", "execution_nonce",
        "resume_policy", "lock", "request_digest",
    }, "campaign request")
    request_copy = dict(request)
    request_digest = request_copy.pop("request_digest")
    require(request_digest == sha256_bytes(canonical_bytes(request_copy)),
            "request digest mismatch")
    require(request["schema"] == SCHEMA and
            request["resume_policy"] ==
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
    require(matrix_path.is_file() and
            strict_json_file(
                matrix_path, "frozen matrix", MAX_MATRIX_JSON_BYTES) ==
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
    validate_campaign_inventory(campaign_root, selected_ids)
    request_path = campaign_root / "request.json"
    require(request_path.is_file() and
            strict_json_file(
                request_path, "retained request", MAX_REQUEST_JSON_BYTES) ==
            request,
            "retained request changed")

    validate_binary_identity_structure(request["binary_identity"])
    validate_execution_identity(
        request["execution_identity"], request["binary_identity"])
    validate_retained_frozen_executables(
        request["binary_identity"], campaign_root)
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
                strict_json_file(
                    cell_path, "retained cell summary", MAX_CELL_JSON_BYTES) ==
                summary,
                "retained cell summary changed")
        require(isinstance(summary, dict) and
                summary.get("cell") == expected_cell and
                summary.get("binary_identity_digest") ==
                    request["binary_identity"]["digest"],
                "cell summary is not bound to request identity")
        invocation_records = summary.get("invocations")
        expected_count = request["rounds"] * len(ORDER)
        require(isinstance(invocation_records, list) and
                len(invocation_records) == expected_count,
                "cell summary has the wrong invocation count")
        retained = []
        retained_attempt_files = set()
        for invocation_index, record in enumerate(invocation_records):
            require_exact_keys(record, {
                "implementation", "execution_identity",
                "median_us", "plan_setup_us",
                "envelope_path", "envelope_sha256",
                "reserved_sibling_nonidle_jiffies",
                "same_user_pair_affinity_before",
                "same_user_pair_affinity_after", "target_runtime",
                "target_interrupts", "isolation_epoch_digest",
            }, "cell invocation record")
            round_index, slot_index = divmod(invocation_index, len(ORDER))
            implementation = ORDER[slot_index]
            require(record["execution_identity"] ==
                    request["execution_identity"][implementation],
                    "cell invocation record used a different immutable "
                    "executable")
            envelope_path = Path(record["envelope_path"])
            attempts = retained_attempt_paths(
                envelope_path, expected_cell, implementation, round_index,
                slot_index, request["max_retries"], campaign_root)
            require(envelope_path == attempts[-1][0],
                    "cell summary does not name the accepted retry")
            for attempt_index, paths in enumerate(attempts[:-1]):
                require(not any(path in retained_attempt_files
                                for path in paths),
                        "retained retry artifact was reused")
                retained_attempt_files.update(paths)
                retry_path = paths[0]
                retry_value = strict_json_file(
                    retry_path, "retained rejected retry envelope",
                    MAX_ENVELOPE_JSON_BYTES)
                validate_retry_envelope(
                    retry_value, retry_path, expected_cell, implementation,
                    round_index, slot_index, attempt_index, request,
                    epoch["start_digest"], campaign_root)
            require(envelope_path not in seen_envelopes and
                    envelope_path.is_file() and
                    sha256_file(
                        envelope_path, MAX_ENVELOPE_JSON_BYTES) ==
                    record["envelope_sha256"],
                    "retained invocation envelope identity changed")
            seen_envelopes.add(envelope_path)
            require(not any(path in retained_attempt_files
                            for path in attempts[-1]),
                    "retained accepted artifact was reused")
            retained_attempt_files.update(attempts[-1])
            value = strict_json_file(
                envelope_path, "retained accepted envelope",
                MAX_ENVELOPE_JSON_BYTES)
            retained.append(validate_envelope(
                value, envelope_path, expected_cell, implementation,
                round_index, slot_index, request, epoch["start_digest"],
                campaign_root, len(attempts) - 1))
        raw_directory = campaign_root / "raw" / expected_cell["id"]
        require(set(raw_directory.iterdir()) == retained_attempt_files,
                "cell raw directory contains unaudited retry artifacts")
        reconstructed = summarize_cell(
            expected_cell, retained, request["rounds"],
            request["binary_identity"])
        require(summary == reconstructed,
                "cell summary does not reconstruct from retained envelopes")
    require(manifest["promotion_analysis"] ==
            analyze_same_source_promotion(
                matrix, cells,
                request["binary_identity"]["comparison_modes"]),
            "retained promotion analysis does not reconstruct from cells")
    print("verified %d uncontaminated ABBA cells" % len(cells))
    return 0


def write_matrix(options: argparse.Namespace) -> int:
    value = matrix_for_preset(options.preset)
    if options.output == Path("-"):
        print(json.dumps(value, sort_keys=True, indent=2))
    else:
        atomic_json(options.output, value, MAX_MATRIX_JSON_BYTES)
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    matrix = commands.add_parser("matrix", help="write a frozen matrix")
    matrix.add_argument("--output", type=Path, default=Path("-"))
    matrix.add_argument(
        "--preset", choices=("core", "large", "tiny"), default="core")
    matrix.set_defaults(function=write_matrix)
    run = commands.add_parser(
        "run", help="run a fresh, complete same-source campaign")
    run.add_argument("--baseline", required=True, type=Path)
    run.add_argument("--candidate", required=True, type=Path)
    run.add_argument("--source-root", required=True, type=Path)
    run.add_argument(
        "--baseline-mode", required=True,
        choices=tuple(MODE_COMPILE_DEFINITIONS))
    run.add_argument(
        "--candidate-mode", required=True,
        choices=tuple(MODE_COMPILE_DEFINITIONS))
    run.add_argument("--output", required=True, type=Path)
    run.add_argument(
        "--preset", choices=("core", "large", "tiny"), default="core")
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
