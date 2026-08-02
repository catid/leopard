#!/usr/bin/env python3
"""Reproduce the frozen pre-promotion GF8 AVX2 R=1 campaign.

The runner never builds code.  Candidate and control are two equal-length hard
links to one frozen ``bench_leopard2`` inode; the only behavioral difference is
``--r1-small-reduction-mode 1`` versus ``0``.  Exact Leopard main is measured in
alternating MAAM/AMMA segments.  Every completed cell is an atomic resume unit.

The path model is intentionally bound to the clean candidate commit below.
Later production selectors may make mode zero resolve to a promoted path, so a
post-promotion binary is rejected rather than silently changing this campaign.

API names in the evidence are deliberately literal.  The Leopard2 benchmark
uses the ordinary one-item batch APIs for encode and reusable-plan
execution, not prevalidated bindings or ``leo2_decode_plan_execute``.  Its
separately timed one-shot lanes really invoke ``leo2_encode`` and
``leo2_decode``.  Exact main invokes the old ``leo_encode`` and ``leo_decode``
APIs.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import re
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA = "leopard2-r1-small-reduction-abba/v1"
MANIFEST_SCHEMA = "leopard2-r1-small-reduction-manifest/v1"
CELL_SCHEMA = "leopard2-r1-small-reduction-cell/v1"
SUMMARY_SCHEMA = "leopard2-r1-small-reduction-summary/v1"
FAILURE_SCHEMA = "leopard2-r1-small-reduction-failure/v1"
MAIN_COMMIT = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
MAIN_SHA256 = \
    "e252aa2c03c1efdda9f7de256ea0d5bf459310d8cef5a9ba69c8cf2619cd3048"
CANDIDATE_COMMIT = "8ff4ed9f7e02343763ab133ff50b70f76d0f4327"
CANDIDATE_TREE = "af404318034483d15148eedb868d9d9eb7f53bdb"
CANDIDATE_SHA256 = \
    "18c9ad5bc664c8a339ab63e6e6c391237b4d4ccf808d9aa7993dfed7d1d086d5"
MODE_SYMBOL = "_ZN12_GLOBAL__N_1L25g_r1_small_reduction_modeE"
LOCK_PATH = Path("/tmp/leopard-gf8-authoritative.lock")
BENCHMARK_CPU = 14
RESERVED_SIBLING = 30
ADDRESS_SPACE_BYTES = 256 * 1024 * 1024
TARGET_INPUT_BYTES_PER_METRIC_SAMPLE = 16 * 1024 * 1024
MIN_REUSE = 64
MAX_REUSE = 262144
MASK64 = (1 << 64) - 1
T95 = {
    3: 4.302652729911275,
    5: 2.7764451051977987,
    9: 2.306004135204166,
}
CHILD_ENVIRONMENT = {
    "LANG": "C",
    "LC_ALL": "C",
    "OMP_DYNAMIC": "FALSE",
    "OMP_NUM_THREADS": "1",
    "OMP_PLACES": "cores",
    "OMP_PROC_BIND": "TRUE",
    "PATH": "/usr/bin:/bin",
    "TZ": "UTC",
}

# A means mode-one Leopard2, B means its same-inode mode-zero control, and M
# means exact Leopard main.  Each round contains one exact-main MAAM/AMMA
# segment and one same-source BAAB/ABBA segment.
ROUND_SEGMENTS = (
    (
        ("main", "candidate", "candidate", "main"),
        ("control", "candidate", "candidate", "control"),
    ),
    (
        ("candidate", "main", "main", "candidate"),
        ("candidate", "control", "control", "candidate"),
    ),
)

API_LANES = {
    "main_encode": {
        "implementation": "exact_leopard_main",
        "api": "leo_encode",
        "metric": "encode_execution",
    },
    "main_decode_one_call": {
        "implementation": "exact_leopard_main",
        "api": "leo_decode",
        "metric": "decode_including_setup",
        "semantics": "old API call includes its erasure setup",
    },
    "leopard2_encode_execution": {
        "implementation": "leopard2",
        "api": "leo2_encode_batch",
        "item_count": 1,
        "batch_preflight_scratch_bytes": 0,
        "metric": "encode_execution",
        "not_api": "leo2_encode",
    },
    "leopard2_one_shot_encode": {
        "implementation": "leopard2",
        "api": "leo2_encode",
        "item_count": 1,
        "metric": "one_shot_encode",
        "not_api": "batch helper or prevalidated binding",
    },
    "leopard2_reused_plan_execution": {
        "implementation": "leopard2",
        "api": "leo2_decode_plan_execute_batch",
        "item_count": 1,
        "batch_preflight_scratch_bytes": 0,
        "execution_scope": "one_loss_direct_xor",
        "metric": "decode_execution",
        "not_api": "leo2_decode_plan_execute",
    },
    "leopard2_decode_plan_setup": {
        "implementation": "leopard2",
        "api": "leo2_decode_plan_create",
        "metric": "decode_plan_setup",
    },
    "leopard2_one_shot_decode": {
        "implementation": "leopard2",
        "api": "leo2_decode",
        "metric": "one_shot_decode_including_setup",
        "not_api": "batch or plan-create-plus-batch execution",
    },
}

ORDINARY_BATCH_ROUTE_PROOF = {
    "required_schema": "leopard2-benchmark-v12",
    "required_batch_item_count": 1,
    "required_loss_count": 1,
    "forbidden_build_fields": ["prevalidated_batch_experiment"],
    "forbidden_metric_fields": [
        "encode_binding_setup", "decode_binding_setup"],
    "required_build_fields": {
        "r1_timed_encode_api":
            "leo2_encode_batch:item_count=1:no_preflight_scratch",
        "r1_timed_one_shot_encode_api": "leo2_encode",
        "r1_timed_reused_decode_api":
            "leo2_decode_plan_execute_batch:item_count=1:"
            "no_preflight_scratch:one_loss_direct_xor",
        "r1_timed_one_shot_decode_api": "leo2_decode:one_loss",
    },
    "required_memory_fields": {
        "encode_batch_preflight_scratch_bytes": 0,
        "decode_batch_preflight_scratch_bytes": 0,
    },
    "path_attestation_scope": (
        "The schema-v12 codec path fields are accepted as execution-path "
        "evidence only when the serialized API scope and zero-preflight "
        "fields prove the ordinary one-item batch APIs and actual leo2_encode "
        "and leo2_decode one-shot lanes. Prevalidated bindings are excluded "
        "to preserve this exact timed API scope even when they resolve the "
        "same reduction path."),
}


class EvidenceError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvidenceError(message)


def validate_frozen_candidate_identity(
    commit: str,
    tree: str,
    executable_sha256: str,
) -> None:
    require(re.fullmatch(r"[0-9a-f]{40}", commit or "") is not None and
            re.fullmatch(r"[0-9a-f]{40}", tree or "") is not None and
            re.fullmatch(r"[0-9a-f]{64}", executable_sha256 or "") is not
                None,
            "candidate identities must be full lowercase hashes")
    require(commit == CANDIDATE_COMMIT and tree == CANDIDATE_TREE and
            executable_sha256 == CANDIDATE_SHA256,
            "candidate is not the frozen pre-promotion source and binary")


def load_module(path: Path, name: str) -> Any:
    specification = importlib.util.spec_from_file_location(name, path)
    require(specification is not None and specification.loader is not None,
            f"cannot load evidence dependency: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


REPOSITORY = Path(__file__).resolve().parents[3]
T4_BASE = load_module(
    REPOSITORY / "experiments/leopard2/gf8_high_encode/"
    "run_t4_packed_terminal_family_abba.py",
    "r1_small_reduction_t4_evidence_base")
BUILD_PROVENANCE = load_module(
    REPOSITORY / "tools/leopard2_build_provenance.py",
    "r1_small_reduction_build_provenance")
MAIN_SUPPORT = T4_BASE.MAIN_SUPPORT
T8_SUPPORT = T4_BASE.T8_SUPPORT
RUNNER_PATH = Path(__file__).resolve()
RUNNER_DEPENDENCIES = (
    RUNNER_PATH,
    Path(T4_BASE.__file__).resolve(),
    Path(T8_SUPPORT.__file__).resolve(),
    Path(MAIN_SUPPORT.__file__).resolve(),
    Path(MAIN_SUPPORT.git_capture.__file__).resolve(),
    Path(BUILD_PROVENANCE.__file__).resolve(),
    Path(MAIN_SUPPORT.link_common.__file__).resolve(),
)


def canonical_bytes(value: object) -> bytes:
    try:
        return json.dumps(
            value, sort_keys=True, separators=(",", ":"),
            ensure_ascii=False, allow_nan=False).encode("utf-8")
    except (TypeError, ValueError, UnicodeEncodeError, RecursionError) as error:
        raise EvidenceError(f"value is not canonical JSON: {error}") from error


def digest_object(value: object) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_identity(path: Path) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    before = resolved.stat()
    require(resolved.is_file() and before.st_size > 0,
            f"identity path is not a nonempty regular file: {resolved}")
    content_hash = sha256(resolved)
    after = resolved.stat()
    fields = (
        "st_dev", "st_ino", "st_mode", "st_uid", "st_gid", "st_nlink",
        "st_size", "st_mtime_ns", "st_ctime_ns")
    require(all(getattr(before, name) == getattr(after, name)
                for name in fields),
            f"file changed while it was hashed: {resolved}")
    return {
        "path": str(resolved),
        "sha256": content_hash,
        "device": after.st_dev,
        "inode": after.st_ino,
        "mode": after.st_mode,
        "uid": after.st_uid,
        "gid": after.st_gid,
        "links": after.st_nlink,
        "size": after.st_size,
        "mtime_ns": after.st_mtime_ns,
        "ctime_ns": after.st_ctime_ns,
    }


def dependency_identities() -> list[dict[str, Any]]:
    resolved = [path.resolve(strict=True) for path in RUNNER_DEPENDENCIES]
    require(len(resolved) == len(set(resolved)),
            "runner dependencies are not unique")
    return [file_identity(path) for path in resolved]


def write_atomic(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    require(path.parent.is_dir(), f"artifact parent is not a directory: {path}")
    temporary = path.parent / (
        f".{path.name}.tmp-{os.getpid()}-{time.monotonic_ns()}")
    descriptor = os.open(
        temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    try:
        offset = 0
        while offset < len(payload):
            offset += os.write(descriptor, payload[offset:])
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    try:
        os.replace(temporary, path)
        directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def write_json_atomic(path: Path, value: object) -> None:
    write_atomic(path, json.dumps(
        value, sort_keys=True, indent=2, allow_nan=False,
        ensure_ascii=False).encode("utf-8") + b"\n")


def read_json(path: Path, label: str) -> object:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise EvidenceError(f"cannot read {label}: {error}") from error


def xorshift64(value: int) -> int:
    value ^= (value << 13) & MASK64
    value ^= value >> 7
    value ^= (value << 17) & MASK64
    return value & MASK64


def selected_loss(k: int, seed: int) -> int:
    order = list(range(k))
    state = (seed ^ 0xD1B54A32D192ED03) & MASK64
    if state == 0:
        state = 0x9E3779B97F4A7C15
    for remaining in range(len(order), 1, -1):
        state = xorshift64(state)
        selected = state % remaining
        order[remaining - 1], order[selected] = \
            order[selected], order[remaining - 1]
    # Both benchmark adapters keep the prefix after the Fisher-Yates shuffle
    # (`order.resize(losses)`), so the one-loss selection is order[0].
    return order[0]


def desired_loss(k: int, position: str) -> int:
    if position in {"first", "only"}:
        return 0
    if position == "middle":
        return k // 2
    if position == "last":
        return k - 1
    raise EvidenceError(f"unknown missing-position label: {position}")


def seed_for_loss(k: int, position: str, start: int) -> int:
    desired = desired_loss(k, position)
    for offset in range(1_000_000):
        seed = (start + offset) & MASK64
        if seed != 0 and selected_loss(k, seed) == desired:
            return seed
    raise EvidenceError("could not derive a deterministic positional seed")


def expected_reduction_path(k: int, shard_bytes: int, mode: int) -> str:
    require(mode in (0, 1), "R=1 diagnostic mode is invalid")
    if k == 1:
        return "k1_copy"
    if mode == 1:
        if k == 2 and shard_bytes in (64, 256, 1024):
            return "k2_terminal"
        if k == 3 and shard_bytes in (64, 256, 1024):
            return "fused_final"
        if k == 4 and shard_bytes in (64, 256, 1024):
            return "dense"
        if k >= 8 and shard_bytes in (64, 256):
            return "coarse"
    if k == 2 and shard_bytes >= 2048:
        return "dense"
    if k == 4 and shard_bytes >= 2048:
        return "dense"
    if shard_bytes >= 4096 and k in (3, 5, 6):
        return "fused_final"
    if shard_bytes >= 1024 and k >= 7:
        return "coarse"
    return "pairwise"


def reuse_for_cell(k: int, shard_bytes: int) -> int:
    """Scale calls so each retained timing sample covers about 16 MiB."""
    require(k > 0 and shard_bytes > 0,
            "reuse geometry must be positive")
    input_bytes_per_call = k * shard_bytes
    nearest = (
        TARGET_INPUT_BYTES_PER_METRIC_SAMPLE + input_bytes_per_call // 2
    ) // input_bytes_per_call
    return max(MIN_REUSE, min(MAX_REUSE, nearest))


def campaign_cells() -> list[dict[str, Any]]:
    specifications: list[tuple[str, int, int, str]] = []

    for k in (2, 3, 4):
        for index, shard_bytes in enumerate((64, 256, 1024)):
            specifications.append((
                "target_small_exact", k, shard_bytes,
                ("first", "middle", "last")[index]))

    # The extra near-boundary K values make the diagnostic useful even when
    # the coarse-reduction crossover is nonmonotonic.  "Target" means that
    # mode one changes the resolved path, not that a speedup is presumed.
    wide_k = (
        8, 9, 10, 11, 12, 14, 15, 16, 17, 23, 24,
        25, 30, 31, 32, 63, 64, 127, 128, 191, 224, 255,
    )
    for index, k in enumerate(wide_k):
        for byte_index, shard_bytes in enumerate((64, 256)):
            specifications.append((
                "target_wide_exact", k, shard_bytes,
                ("first", "middle", "last")[(2 * index + byte_index) % 3]))

    for index, k in enumerate((1, 5, 6, 7)):
        for byte_index, shard_bytes in enumerate((64, 256, 1024)):
            position = "only" if k == 1 else \
                ("first", "middle", "last")[(index + byte_index) % 3]
            specifications.append(("control_k", k, shard_bytes, position))

    # Exact Leopard main accepts only positive shard sizes divisible by 64.
    # These aligned sizes remain inert in both diagnostic modes while probing
    # byte-count sensitivity away from the selected 64/256/1024 cells.
    aligned_controls = (128, 192, 320, 512, 768, 1088)
    for k_index, k in enumerate((2, 3, 4)):
        for byte_index, shard_bytes in enumerate(aligned_controls):
            specifications.append((
                "control_small_aligned_bytes", k, shard_bytes,
                ("first", "middle", "last")[(k_index + byte_index) % 3]))

    for byte_index, shard_bytes in enumerate(aligned_controls):
        specifications.append((
            "control_wide_aligned_bytes", 8, shard_bytes,
            ("first", "middle", "last")[byte_index % 3]))

    for index, k in enumerate(wide_k):
        specifications.append((
            "control_wide_b1024", k, 1024,
            ("first", "middle", "last")[index % 3]))

    cells: list[dict[str, Any]] = []
    for index, (role, k, shard_bytes, position) in enumerate(specifications):
        seed = seed_for_loss(k, position, 0x523100000000 + index * 0x1000)
        expected_missing = desired_loss(k, position)
        candidate_path = expected_reduction_path(k, shard_bytes, 1)
        control_path = expected_reduction_path(k, shard_bytes, 0)
        selected = candidate_path != control_path
        require(selected == role.startswith("target_"),
                "directional role disagrees with the R=1 policy")
        cells.append({
            "id": (f"{role}-k{k}-r1-b{shard_bytes}-"
                   f"m{position}-q1"),
            "K": k,
            "R": 1,
            "bytes": shard_bytes,
            "batch": 1,
            "reuse": reuse_for_cell(k, shard_bytes),
            "loss": 1,
            "missing_position": position,
            "expected_missing_original_indices": [expected_missing],
            "role": role,
            "candidate_selected": selected,
            "candidate_reduction_path": candidate_path,
            "control_reduction_path": control_path,
            "seed": seed,
        })

    expected_roles = {
        "target_small_exact": 9,
        "target_wide_exact": 44,
        "control_k": 12,
        "control_small_aligned_bytes": 18,
        "control_wide_aligned_bytes": 6,
        "control_wide_b1024": 22,
    }
    role_counts = {
        role: sum(cell["role"] == role for cell in cells)
        for role in expected_roles
    }
    require(
        len(cells) == 111 and role_counts == expected_roles and
        sum(bool(cell["candidate_selected"]) for cell in cells) == 53 and
        len({cell["id"] for cell in cells}) == len(cells) and
        len({cell["seed"] for cell in cells}) == len(cells) and
        all(cell["bytes"] > 0 and cell["bytes"] % 64 == 0
            for cell in cells) and
        all(cell["reuse"] == reuse_for_cell(cell["K"], cell["bytes"])
            for cell in cells) and
        all(selected_loss(cell["K"], cell["seed"]) ==
            cell["expected_missing_original_indices"][0] for cell in cells),
        "R=1 small-reduction directional matrix is incomplete")
    return cells


def inspect_isa_disassembly(disassembly: str) -> dict[str, int]:
    evex = 0
    ymm = 0
    for line in disassembly.splitlines():
        fields = line.lstrip().split()
        if fields and fields[0].endswith(":") and len(fields) > 1 and \
                fields[1].lower() == "62":
            evex += 1
        if re.search(r"\bymm[0-9]+\b", line):
            ymm += 1
    return {
        "evex_prefixed_instruction_count": evex,
        "ymm_operand_instruction_count": ymm,
    }


def inspect_pure_avx2(executable: Path, label: str) -> dict[str, Any]:
    completed = subprocess.run(
        ["/usr/bin/objdump", "-d", "-M", "intel", str(executable)],
        env=CHILD_ENVIRONMENT, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, timeout=180.0, check=False)
    require(completed.returncode == 0,
            f"objdump failed for {label}: " +
            completed.stderr.decode("utf-8", errors="replace")[-1000:])
    try:
        disassembly = completed.stdout.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise EvidenceError(f"objdump output for {label} is not UTF-8") from error
    result = inspect_isa_disassembly(disassembly)
    require(result["evex_prefixed_instruction_count"] == 0,
            f"{label} contains EVEX-prefixed instructions")
    require(result["ymm_operand_instruction_count"] > 0,
            f"{label} contains no visible AVX2 YMM instructions")
    version = subprocess.run(
        ["/usr/bin/objdump", "--version"], env=CHILD_ENVIRONMENT,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        timeout=10.0, check=False)
    require(version.returncode == 0 and version.stdout.splitlines(),
            "objdump --version failed")
    result["objdump_version"] = version.stdout.splitlines()[0].decode(
        "utf-8", errors="strict")
    return result


def mode_word_identity(executable: Path) -> dict[str, Any]:
    """Map the retained selector, including its ELF NOBITS zero initializer."""
    symbols = T4_BASE.run_tool((
        "/usr/bin/readelf", "-sW", str(executable)))
    matches = []
    for line in symbols.splitlines():
        tokens = line.split()
        if len(tokens) >= 8 and tokens[-1] == MODE_SYMBOL:
            matches.append(tokens)
    require(len(matches) == 1,
            "R=1 diagnostic selector symbol is missing or ambiguous")
    symbol = matches[0]
    require(symbol[2] == "4" and
            symbol[3:6] == ["OBJECT", "LOCAL", "DEFAULT"],
            "R=1 diagnostic selector symbol metadata changed")
    try:
        address = int(symbol[1], 16)
        section_index = int(symbol[6])
    except ValueError as error:
        raise EvidenceError(
            "R=1 diagnostic selector symbol is not section-backed") from error

    sections = T4_BASE.run_tool((
        "/usr/bin/readelf", "-SW", str(executable)))
    pattern = re.compile(
        r"^\s*\[\s*(\d+)\]\s+(\S+)\s+(\S+)\s+"
        r"([0-9A-Fa-f]+)\s+([0-9A-Fa-f]+)\s+([0-9A-Fa-f]+)\s+")
    selected: tuple[str, str, int, int, int] | None = None
    for line in sections.splitlines():
        match = pattern.match(line)
        if match and int(match.group(1)) == section_index:
            selected = (
                match.group(2), match.group(3), int(match.group(4), 16),
                int(match.group(5), 16), int(match.group(6), 16))
            break
    require(selected is not None,
            "R=1 diagnostic selector section is absent")
    section_name, section_type, section_address, section_offset, \
        section_size = selected
    require(section_address <= address and
            address + 4 <= section_address + section_size,
            "R=1 diagnostic selector lies outside its section")
    if section_type == "NOBITS":
        # ELF guarantees that the loader supplies zero bytes for SHT_NOBITS.
        # The production selector intentionally defaults to zero, so unlike
        # older nonzero T=4 selectors it has no file payload to read.
        require(section_name in {".bss", ".sbss"},
                "zero-initialized selector is in an unexpected NOBITS section")
        payload = b"\0" * 4
        file_offset: int | None = None
        storage = "elf_nobits_zero_fill"
    else:
        require(section_type == "PROGBITS",
                "R=1 selector section is neither PROGBITS nor NOBITS")
        file_offset = section_offset + address - section_address
        with executable.open("rb") as source:
            require(source.read(6) == b"\x7fELF\x02\x01",
                    "benchmark is not little-endian ELF64")
            source.seek(file_offset)
            payload = source.read(4)
        require(len(payload) == 4,
                "R=1 diagnostic selector payload is truncated")
        storage = "elf_progbits_payload"
    return {
        "symbol": MODE_SYMBOL,
        "virtual_address": address,
        "section_index": section_index,
        "section_name": section_name,
        "section_type": section_type,
        "file_offset": file_offset,
        "storage": storage,
        "bytes_hex": payload.hex(),
        "value": int.from_bytes(payload, "little"),
    }


def validate_shared_binary_identities(
    identities: Mapping[str, Mapping[str, Any]],
    expected_candidate_sha256: str,
    expected_main_sha256: str,
) -> dict[str, Any]:
    require(set(identities) == {"candidate", "control", "main"},
            "binary identity set is incomplete")
    candidate = identities["candidate"]
    control = identities["control"]
    main = identities["main"]
    require(candidate["sha256"] == control["sha256"] ==
            expected_candidate_sha256 and main["sha256"] ==
            expected_main_sha256 and candidate["sha256"] != main["sha256"],
            "frozen benchmark SHA-256 identity differs")
    require((candidate["device"], candidate["inode"]) ==
            (control["device"], control["inode"]) and
            candidate["links"] >= 2 and control["links"] >= 2,
            "candidate/control are not hard links to one inode")
    require(len({candidate["path"], control["path"], main["path"]}) == 3,
            "benchmark pathnames are not distinct")
    lengths = {name: len(str(identity["path"]))
               for name, identity in identities.items()}
    require(len(set(lengths.values())) == 1,
            "candidate/control/main executable paths have unequal lengths")
    return {
        "candidate_control_device": candidate["device"],
        "candidate_control_inode": candidate["inode"],
        "candidate_control_link_count": candidate["links"],
        "executable_path_lengths": lengths,
    }


def validate_lane_owned_copy_identities(
    identities: Mapping[str, Mapping[str, Any]],
    provenance: Mapping[str, Any],
) -> dict[str, Any]:
    candidate_build = provenance.get("candidate_build")
    main_build = provenance.get("main_build")
    require(isinstance(candidate_build, Mapping) and
            isinstance(main_build, Mapping),
            "build provenance is incomplete")
    candidate_build_executable = candidate_build.get("executable")
    main_build_executable = main_build.get("validated_executable")
    require(isinstance(candidate_build_executable, Mapping) and
            isinstance(main_build_executable, Mapping),
            "build executable provenance is incomplete")
    candidate = identities.get("candidate")
    control = identities.get("control")
    main = identities.get("main")
    require(isinstance(candidate, Mapping) and
            isinstance(control, Mapping) and isinstance(main, Mapping),
            "frozen lane identity set is incomplete")
    require(candidate_build_executable.get("sha256") ==
                candidate.get("sha256") == control.get("sha256") and
            main_build_executable.get("sha256") == main.get("sha256"),
            "frozen executable bytes differ from exact build provenance")
    candidate_build_inode = (
        candidate_build_executable.get("device"),
        candidate_build_executable.get("inode"))
    main_build_inode = (
        main_build_executable.get("device"),
        main_build_executable.get("inode"))
    candidate_lane_inode = (candidate.get("device"), candidate.get("inode"))
    control_lane_inode = (control.get("device"), control.get("inode"))
    main_lane_inode = (main.get("device"), main.get("inode"))
    require(candidate_lane_inode == control_lane_inode and
            candidate_lane_inode != candidate_build_inode and
            main_lane_inode != main_build_inode,
            "benchmark lanes are hard links to mutable build artifacts instead "
            "of lane-owned copies")
    return {
        "candidate_build_device_inode": list(candidate_build_inode),
        "candidate_lane_device_inode": list(candidate_lane_inode),
        "main_build_device_inode": list(main_build_inode),
        "main_lane_device_inode": list(main_lane_inode),
        "lane_inodes_distinct_from_build_artifacts": True,
    }


def benchmark_command(
    implementation: str,
    executable: Path,
    cell: Mapping[str, Any],
    cpu: int,
    iterations: int,
    warmup: int,
) -> list[str]:
    common = [
        "/usr/bin/prlimit", f"--as={ADDRESS_SPACE_BYTES}",
        "/usr/bin/taskset", "-c", str(cpu), str(executable),
        "--k", str(cell["K"]), "--r", "1",
        "--bytes", str(cell["bytes"]), "--loss", "1",
        "--batch", "1", "--reuse", str(cell["reuse"]),
        "--iterations", str(iterations), "--warmup", str(warmup),
        "--threads", "1", "--seed", str(cell["seed"]), "--json", "-",
    ]
    if implementation == "main":
        return common
    require(implementation in {"candidate", "control"},
            "unknown Leopard2 implementation lane")
    mode = "1" if implementation == "candidate" else "0"
    return common[:-2] + [
        "--profile", "high", "--field", "gf8", "--backend", "avx2",
        "--skip-legacy", "--retain-samples", "--measure-one-shot-decode",
        "--r1-small-reduction-mode", mode, "--json", "-",
    ]


def positive_number(value: object, label: str) -> float:
    require(isinstance(value, (int, float)) and not isinstance(value, bool) and
            math.isfinite(float(value)) and float(value) > 0,
            f"{label} is not finite and positive")
    return float(value)


def metric_median(
    metrics: Mapping[str, Any], name: str, key: str,
    iterations: int, samples_key: str | None,
) -> float:
    summary = metrics.get(name)
    require(isinstance(summary, Mapping), f"metric {name} is absent")
    median = positive_number(summary.get(key), f"metric {name}/{key}")
    if samples_key is not None:
        samples = summary.get(samples_key)
        require(isinstance(samples, list) and len(samples) == iterations and
                all(isinstance(item, (int, float)) and
                    not isinstance(item, bool) and
                    math.isfinite(float(item)) and float(item) >= 0
                    for item in samples),
                f"metric {name} retained samples are invalid")
    return median


def validate_result(
    implementation: str,
    result: object,
    cell: Mapping[str, Any],
    iterations: int,
    warmup: int,
) -> dict[str, Any]:
    require(isinstance(result, Mapping), "benchmark result is not an object")
    parameters = result.get("parameters")
    require(isinstance(parameters, Mapping), "benchmark parameters are absent")
    expected = {
        "K": cell["K"], "R": 1, "shard_bytes": cell["bytes"],
        "loss_count": 1, "missing_original_indices":
            cell["expected_missing_original_indices"],
        "batch": 1, "reuse": cell["reuse"], "iterations": iterations,
        "warmup": warmup, "thread_count": 1, "seed": cell["seed"],
    }
    require(all(parameters.get(name) == value for name, value in expected.items()),
            "benchmark parameters or missing-position selection differ")
    resolved = result.get("resolved")
    correctness = result.get("correctness")
    digests = result.get("workload_digests")
    require(isinstance(resolved, Mapping) and
            resolved.get("profile") == "legacy_high_v1" and
            resolved.get("field") == "gf8" and
            resolved.get("thread_count") == 1 and
            resolved.get("padded_side") == 1,
            "benchmark resolved a different code geometry")
    require(isinstance(digests, Mapping) and
            digests.get("algorithm") == "fnv1a64" and
            all(isinstance(digests.get(name), str) and
                re.fullmatch(r"[0-9a-f]{16}", str(digests[name])) is not None
                for name in ("original_data", "transmitted_parity",
                             "recovered_originals")),
            "benchmark workload digests are incomplete")
    metrics = result.get("metrics")
    require(isinstance(metrics, Mapping), "benchmark metrics are absent")

    if implementation == "main":
        require(result.get("schema") == "leopard-main-benchmark-v1" and
                isinstance(correctness, Mapping) and
                correctness.get("round_trip") is True,
                "exact-main schema or round trip differs")
        build = result.get("build")
        require(isinstance(build, Mapping) and
                build.get("main_source_commit") == MAIN_COMMIT and
                metrics.get("codec_setup") is None and
                metrics.get("decode_timing_includes_setup") is True,
                "exact-main source or API timing semantics differ")
        return {
            "metrics_us": {
                "legacy_encode": metric_median(
                    metrics, "encode_execution", "median_us_per_batch_call",
                    iterations, None),
                "legacy_decode": metric_median(
                    metrics, "decode_including_setup",
                    "median_us_per_batch_call", iterations, None),
            },
            "digests": dict(digests),
            "api_lanes": {
                "legacy_encode": API_LANES["main_encode"],
                "legacy_decode": API_LANES["main_decode_one_call"],
            },
        }

    require(implementation in {"candidate", "control"},
            "unknown Leopard2 result lane")
    require(result.get("schema") == "leopard2-benchmark-v12" and
            resolved.get("backend") == "avx2" and
            isinstance(correctness, Mapping) and
            correctness.get("leopard2_round_trip") is True and
            parameters.get("requested_profile") == "legacy_high_v1" and
            parameters.get("requested_field") == "gf8" and
            parameters.get("requested_backend") == "avx2" and
            parameters.get("force_generic_decode") is False and
            parameters.get("force_specialized_decode") is False and
            parameters.get("force_tiled_decode") is False and
            parameters.get("force_materialized_decode") is False and
            parameters.get("measure_one_shot_decode") is True and
            parameters.get("skip_legacy") is True and
            parameters.get("retain_samples") is True,
            "Leopard2 schema, backend, one-shot route, or round trip differs")
    build = result.get("build")
    memory = result.get("memory")
    require(isinstance(build, Mapping) and
            isinstance(memory, Mapping) and
            all(name not in build for name in
                ORDINARY_BATCH_ROUTE_PROOF["forbidden_build_fields"]) and
            all(name not in metrics for name in
                ORDINARY_BATCH_ROUTE_PROOF["forbidden_metric_fields"]) and
            all(build.get(name) == value for name, value in
                ORDINARY_BATCH_ROUTE_PROOF["required_build_fields"].items()) and
            all(memory.get(name) == value for name, value in
                ORDINARY_BATCH_ROUTE_PROOF["required_memory_fields"].items()),
            "benchmark selected prevalidated bindings instead of ordinary "
            "zero-preflight one-item batch APIs, or API scope differs")
    mode = 1 if implementation == "candidate" else 0
    expected_path = str(cell[
        "candidate_reduction_path" if mode == 1 else
        "control_reduction_path"])
    require(build.get("r1_small_reduction_diagnostic_mode") == mode and
            build.get("r1_small_reduction_codec_enabled") is (mode == 1) and
            build.get("r1_encode_reduction_path") == expected_path and
            build.get("r1_decode_reduction_path") == expected_path,
            "resolved R=1 reduction mode or execution path differs")
    codec_setup = metric_median(
        metrics, "codec_setup", "median_us", iterations, "samples_us")
    encode = metric_median(
        metrics, "encode_execution", "median_us_per_batch_call",
        iterations, "samples_us_per_batch_call")
    one_shot_encode = metric_median(
        metrics, "one_shot_encode", "median_us_per_batch_call",
        iterations, "samples_us_per_batch_call")
    plan_setup = metric_median(
        metrics, "decode_plan_setup", "median_us", iterations, "samples_us")
    decode_execution = metric_median(
        metrics, "decode_execution", "median_us_per_batch_call",
        iterations, "samples_us_per_batch_call")
    one_shot = metric_median(
        metrics, "one_shot_decode_including_setup",
        "median_us_per_batch_call", iterations,
        "samples_us_per_batch_call")
    amortized = metrics.get("decode_amortized_at_reuse")
    require(isinstance(amortized, Mapping) and
            amortized.get("reuse_count") == cell["reuse"],
            "reused-plan amortization semantics differ")
    derived = positive_number(
        amortized.get("derived_median_us_per_batch_call"),
        "decode amortized median")
    expected_derived = decode_execution + plan_setup / float(cell["reuse"])
    require(abs(derived - expected_derived) <=
            max(1e-6, 5e-6 * expected_derived),
            "decode amortized result is not execution plus setup/reuse")
    return {
        "metrics_us": {
            "codec_setup": codec_setup,
            "batch_encode": encode,
            "one_shot_encode": one_shot_encode,
            "decode_plan_setup": plan_setup,
            "decode_execution": decode_execution,
            "decode_first_use": decode_execution + plan_setup,
            "decode_reuse_amortized": derived,
            "one_shot_decode": one_shot,
        },
        "digests": dict(digests),
        "mode": mode,
        "reduction_path": expected_path,
        "api_lanes": {
            "batch_encode": API_LANES["leopard2_encode_execution"],
            "one_shot_encode": API_LANES["leopard2_one_shot_encode"],
            "decode_execution":
                API_LANES["leopard2_reused_plan_execution"],
            "decode_plan_setup": API_LANES["leopard2_decode_plan_setup"],
            "one_shot_decode": API_LANES["leopard2_one_shot_decode"],
        },
    }


def run_one(
    implementation: str,
    identity: Mapping[str, Any],
    cell: Mapping[str, Any],
    cpu: int,
    iterations: int,
    warmup: int,
    timeout: float,
) -> dict[str, Any]:
    executable = Path(str(identity["path"]))
    require(file_identity(executable) == identity,
            f"{implementation} binary changed before execution")
    command = benchmark_command(
        implementation, executable, cell, cpu, iterations, warmup)
    started_ns = time.monotonic_ns()
    try:
        completed = subprocess.run(
            command, env=CHILD_ENVIRONMENT, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, timeout=timeout, check=False)
    except subprocess.TimeoutExpired as error:
        raise EvidenceError(
            f"{implementation} benchmark timed out after {timeout}s") from error
    elapsed_ns = time.monotonic_ns() - started_ns
    require(completed.returncode == 0,
            f"{implementation} benchmark failed: " +
            completed.stderr.decode("utf-8", errors="replace")[-1000:])
    require(file_identity(executable) == identity,
            f"{implementation} binary changed after execution")
    try:
        result = json.loads(completed.stdout.decode("utf-8", errors="strict"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise EvidenceError(
            f"{implementation} benchmark did not emit one JSON object") from error
    normalized = validate_result(
        implementation, result, cell, iterations, warmup)
    return {
        "implementation": implementation,
        "command": command,
        "elapsed_ns": elapsed_ns,
        "returncode": completed.returncode,
        "stdout_sha256": hashlib.sha256(completed.stdout).hexdigest(),
        "stderr_sha256": hashlib.sha256(completed.stderr).hexdigest(),
        "normalized": normalized,
        "result": result,
    }


def confidence_interval(log_ratios: Sequence[float]) -> dict[str, Any]:
    require(len(log_ratios) in T95,
            "round count has no predeclared Student-t threshold")
    center = statistics.mean(log_ratios)
    half_width = T95[len(log_ratios)] * statistics.stdev(log_ratios) / \
        math.sqrt(len(log_ratios))
    return {
        "speedup": math.exp(center),
        "speedup_definition": "baseline_time/candidate_time",
        "ci95": [math.exp(center - half_width),
                 math.exp(center + half_width)],
        "round_log_ratios": list(log_ratios),
    }


def analyze_cell(
    cell: Mapping[str, Any],
    rounds: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    require(bool(rounds), "cell has no completed rounds")
    reference = rounds[0]["segments"][0]["invocations"][0][
        "normalized"]["digests"]
    same_source_names = (
        "codec_setup", "batch_encode", "one_shot_encode",
        "decode_plan_setup", "decode_execution", "decode_first_use",
        "decode_reuse_amortized", "one_shot_decode")
    exact_main_comparisons = {
        "legacy_encode_vs_batch_encode":
            ("legacy_encode", "batch_encode"),
        "legacy_encode_vs_one_shot_encode":
            ("legacy_encode", "one_shot_encode"),
        "legacy_decode_vs_one_shot_decode":
            ("legacy_decode", "one_shot_decode"),
    }
    logs = {
        "candidate_vs_control": {name: [] for name in same_source_names},
        "candidate_vs_main": {
            name: [] for name in exact_main_comparisons},
    }
    for round_index, round_value in enumerate(rounds):
        require(round_value.get("round") == round_index,
                "round sequence is not contiguous and zero-based")
        isolation = round_value.get("isolation")
        require(isinstance(isolation, Mapping) and
                isolation.get("accepted") is True and
                isolation.get("delta", {}).get(
                    "reserved_sibling", {}).get("nonidle_jiffies") == 0,
                "contaminated round cannot be analyzed")
        segments = round_value.get("segments")
        require(isinstance(segments, list) and len(segments) == 2,
                "round does not contain two balanced segments")
        expected_segments = ROUND_SEGMENTS[round_index % 2]
        for segment_index, segment in enumerate(segments):
            require(isinstance(segment, Mapping),
                    "round segment is not an object")
            expected_order = expected_segments[segment_index]
            invocations = segment.get("invocations")
            require(isinstance(invocations, list) and
                    segment.get("order") == list(expected_order) and
                    [item.get("implementation")
                     if isinstance(item, Mapping) else None
                     for item in invocations] == list(expected_order) and
                    len(invocations) == 4 and
                    all(item["normalized"]["digests"] == reference
                        for item in invocations),
                    "balanced segment order or workload digests differ")
            baseline = segment.get("baseline")
            expected_baseline = "main" if "main" in expected_order \
                else "control"
            require(baseline == expected_baseline,
                    "segment baseline differs from its declared order")
            candidates = [item for item in invocations
                          if item["implementation"] == "candidate"]
            baselines = [item for item in invocations
                         if item["implementation"] == baseline]
            require(len(candidates) == len(baselines) == 2,
                    "balanced segment lacks two observations per lane")
            contrast = "candidate_vs_main" if baseline == "main" \
                else "candidate_vs_control"
            comparisons = exact_main_comparisons if baseline == "main" else {
                name: (name, name) for name in same_source_names}
            for name, (baseline_name, candidate_name) in comparisons.items():
                candidate_log = statistics.mean(math.log(
                    item["normalized"]["metrics_us"][candidate_name])
                    for item in candidates)
                baseline_log = statistics.mean(math.log(
                    item["normalized"]["metrics_us"][baseline_name])
                    for item in baselines)
                logs[contrast][name].append(baseline_log - candidate_log)
    return {
        "cell": dict(cell),
        "digests": dict(reference),
        "candidate_vs_control": {
            name: confidence_interval(values)
            for name, values in logs["candidate_vs_control"].items()
        },
        "candidate_vs_main": {
            name: confidence_interval(values)
            for name, values in logs["candidate_vs_main"].items()
        },
    }


def cell_artifact_path(output: Path, cell: Mapping[str, Any]) -> Path:
    identifier = str(cell["id"])
    require(re.fullmatch(r"[a-z0-9_-]+", identifier) is not None,
            "cell identifier is not pathname-safe")
    return output / "cells" / f"{identifier}.json"


def validate_cell_artifact(
    value: object,
    cell: Mapping[str, Any],
    contract_sha256: str,
    rounds: int,
) -> dict[str, Any]:
    require(isinstance(value, Mapping) and value.get("schema") == CELL_SCHEMA and
            value.get("contract_sha256") == contract_sha256 and
            value.get("cell") == cell and
            isinstance(value.get("rounds"), list) and
            len(value["rounds"]) == rounds and
            isinstance(value.get("analysis"), Mapping),
            "resumed cell artifact does not match the campaign contract")
    expected = analyze_cell(cell, value["rounds"])
    require(digest_object(expected) == digest_object(value["analysis"]),
            "resumed cell analysis differs from its retained rounds")
    result = dict(value)
    require(result.get("artifact_sha256") == digest_object({
                key: item for key, item in result.items()
                if key != "artifact_sha256"}),
            "resumed cell artifact digest differs")
    return result


def capture_sources_and_builds(options: argparse.Namespace) -> dict[str, Any]:
    candidate_source = options.source_root.resolve(strict=True)
    main_source = options.main_source_root.resolve(strict=True)
    try:
        sources = MAIN_SUPPORT.git_capture.capture_git_identities((
            (candidate_source, options.source_commit, False),
            (main_source, MAIN_COMMIT, True),
        ))
    except Exception as error:
        raise EvidenceError(f"source capture failed: {error}") from error
    require(sources[0].get("tree") == options.source_tree and
            sources[0].get("tracked_status") == "clean" and
            sources[1].get("tree") is not None and
            sources[1].get("tracked_status") == "clean",
            "source tree or tracked cleanliness differs")

    candidate_build = options.build_dir.resolve(strict=True)
    candidate_executable = candidate_build / "bench_leopard2"
    try:
        candidate_provenance = BUILD_PROVENANCE.candidate_build_provenance(
            candidate_build, candidate_source, candidate_executable,
            "bench_leopard2")
    except Exception as error:
        raise EvidenceError(f"candidate build provenance failed: {error}") \
            from error

    main_build = options.main_build_dir.resolve(strict=True)
    main_build_executable = main_build / "leopard_main_benchmark"
    main_archive = main_build / "libleopard_main_exact.a"
    main_specification = {
        "baseline_build_dir": str(main_build),
        "baseline_executable": str(main_build_executable),
        "baseline_archive": str(main_archive),
        "baseline_source_root": str(main_source),
        "candidate_source_root": str(candidate_source),
        "baseline_pure_avx2": True,
    }
    try:
        main_provenance = MAIN_SUPPORT.build_provenance(
            "baseline", main_specification)
    except Exception as error:
        raise EvidenceError(f"exact-main build provenance failed: {error}") \
            from error
    return {
        "candidate_source": sources[0],
        "main_source": sources[1],
        "candidate_build": candidate_provenance,
        "main_build": main_provenance,
    }


def make_contract(
    options: argparse.Namespace,
    cells: Sequence[Mapping[str, Any]],
    identities: Mapping[str, Mapping[str, Any]],
    shared: Mapping[str, Any],
    provenance: Mapping[str, Any],
    isa: Mapping[str, Any],
    dependencies: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    lane_owned_copy_proof = validate_lane_owned_copy_identities(
        identities, provenance)
    mode_word = mode_word_identity(
        Path(str(identities["candidate"]["path"])))
    require(mode_word["value"] == 0,
            "shared production binary does not default R=1 mode to zero")
    candidate_example = benchmark_command(
        "candidate", Path(str(identities["candidate"]["path"])),
        cells[0], options.cpu, options.iterations, options.warmup)
    control_example = benchmark_command(
        "control", Path(str(identities["control"]["path"])),
        cells[0], options.cpu, options.iterations, options.warmup)
    require(len(candidate_example) == len(control_example) and
            all(len(left) == len(right)
                for left, right in zip(candidate_example, control_example)),
            "candidate/control argv allocation shape differs")
    return {
        "schema": SCHEMA,
        "main_commit": MAIN_COMMIT,
        "source_commit": options.source_commit,
        "source_tree": options.source_tree,
        "binary_identities": dict(identities),
        "shared_inode_and_equal_path_proof": dict(shared),
        "lane_owned_copy_proof": lane_owned_copy_proof,
        "source_archive_build_provenance": dict(provenance),
        "pure_avx2_no_evex": dict(isa),
        "mode_word_default": mode_word,
        "runtime_attribution": {
            "candidate_arguments": ["--r1-small-reduction-mode", "1"],
            "control_arguments": ["--r1-small-reduction-mode", "0"],
            "only_runtime_difference": True,
            "candidate_control_equal_argument_lengths": True,
        },
        "api_lane_contract": API_LANES,
        "ordinary_batch_route_proof": ORDINARY_BATCH_ROUTE_PROOF,
        "benchmark_api_scope": (
            "Schema v12 serializes the exact zero-preflight APIs timed by "
            "this lane: leo2_encode_batch(item_count=1), "
            "leo2_decode_plan_execute_batch(item_count=1), plus the distinct "
            "actual leo2_encode and leo2_decode one-shot calls. The runner "
            "reports batch and one-shot encode separately and rejects binding "
            "and "
            "with-preflight-scratch variants."),
        "isolation": {
            "canonical_lock": str(LOCK_PATH),
            "benchmark_cpu": options.cpu,
            "reserved_sibling": options.sibling,
            "child_address_space_bytes": ADDRESS_SPACE_BYTES,
            "child_environment": CHILD_ENVIRONMENT,
        },
        "timing": {
            "rounds": options.rounds,
            "iterations": options.iterations,
            "warmup": options.warmup,
            "reuse_policy": {
                "formula": (
                    "clamp(nearest_integer(16777216 / (K * shard_bytes)), "
                    "64, 262144)"),
                "target_input_bytes_per_metric_sample":
                    TARGET_INPUT_BYTES_PER_METRIC_SAMPLE,
                "minimum": MIN_REUSE,
                "maximum": MAX_REUSE,
                "serialized_per_cell": True,
                "same_value_for_main_candidate_control": True,
            },
            "reuse_values": sorted({cell["reuse"] for cell in cells}),
            "round_segments": [
                [list(segment) for segment in ROUND_SEGMENTS[index % 2]]
                for index in range(options.rounds)],
            "method": (
                "alternating exact-main MAAM/AMMA and same-source "
                "BAAB/ABBA geometric contrasts; Student-t 95% interval over "
                "independent round log contrasts"),
        },
        "matrix": [dict(cell) for cell in cells],
        "matrix_role_semantics": (
            "target means mode one resolves a different reduction path; it "
            "does not presume a performance win. Controls require the same "
            "resolved path in both modes."),
        "runner_dependencies": list(dependencies),
    }


def verify_contract_stability(
    options: argparse.Namespace,
    identities: Mapping[str, Mapping[str, Any]],
    initial_provenance: Mapping[str, Any],
    dependencies: Sequence[Mapping[str, Any]],
) -> None:
    current_identities = {
        name: file_identity(Path(str(identity["path"])))
        for name, identity in identities.items()
    }
    require(current_identities == dict(identities),
            "frozen binary identities changed during the campaign")
    require(dependency_identities() == list(dependencies),
            "runner dependency changed during the campaign")
    final_provenance = capture_sources_and_builds(options)
    require(digest_object(final_provenance) == digest_object(initial_provenance),
            "source/archive/build provenance changed during the campaign")


def run_campaign(options: argparse.Namespace) -> int:
    require(options.cpu == BENCHMARK_CPU and
            options.sibling == RESERVED_SIBLING,
            "authoritative campaign is frozen to CPU14 and sibling CPU30")
    require(options.rounds in T95 and options.iterations >= 3 and
            options.warmup >= 1 and options.timeout > 0 and
            math.isfinite(options.timeout),
            "benchmark timing controls are invalid")
    validate_frozen_candidate_identity(
        options.source_commit or "", options.source_tree or "",
        options.candidate_sha256 or "")
    require(re.fullmatch(r"[0-9a-f]{64}", options.main_sha256 or ""),
            "exact-main identity must be a full lowercase hash")
    require(options.main_sha256 == MAIN_SHA256,
            "exact-main SHA-256 is not the frozen pure-AVX2 baseline")
    output = options.output.resolve()
    output.mkdir(parents=True, exist_ok=True, mode=0o700)
    require(output.is_dir(), "output is not a directory")
    lock_descriptor = T4_BASE.acquire_global_lock()
    failure_context: dict[str, Any] = {"schema": FAILURE_SCHEMA}
    try:
        require(set(os.sched_getaffinity(0)) == {options.cpu},
                "runner must be singleton-pinned to the benchmark CPU")
        cells = campaign_cells()
        identities = {
            "candidate": file_identity(options.candidate),
            "control": file_identity(options.control),
            "main": file_identity(options.main),
        }
        shared = validate_shared_binary_identities(
            identities, options.candidate_sha256, options.main_sha256)
        dependencies = dependency_identities()
        provenance = capture_sources_and_builds(options)
        isa = {
            "candidate_control": inspect_pure_avx2(
                options.candidate, "shared candidate/control benchmark"),
            "main": inspect_pure_avx2(options.main, "exact-main benchmark"),
        }
        contract = make_contract(
            options, cells, identities, shared, provenance, isa, dependencies)
        contract_sha256 = digest_object(contract)
        manifest_path = output / "manifest.json"
        manifest = {
            "schema": MANIFEST_SCHEMA,
            "contract_sha256": contract_sha256,
            "contract": contract,
        }
        if manifest_path.exists():
            existing = read_json(manifest_path, "existing manifest")
            require(digest_object(existing) == digest_object(manifest),
                    "existing manifest belongs to a different campaign")
        else:
            write_json_atomic(manifest_path, manifest)
        failure_context.update({
            "contract_sha256": contract_sha256,
            "manifest": str(manifest_path),
        })

        sibling_text = Path(
            f"/sys/devices/system/cpu/cpu{options.cpu}/topology/"
            "thread_siblings_list").read_text(encoding="ascii")
        require(MAIN_SUPPORT.parse_cpu_list(sibling_text) ==
                {options.cpu, options.sibling},
                "requested CPUs are not one SMT pair")
        # Retain only compact analyses in memory.  Raw invocation JSON is
        # durably stored per cell and can be tens of megabytes over nine
        # rounds; accumulating it would undermine the campaign's OOM bounds.
        analyses: list[Mapping[str, Any]] = []
        with MAIN_SUPPORT.StableLeaseAnchor(), \
                MAIN_SUPPORT.PairLease(
                    options.cpu, options.sibling) as pair_lease:
            before_cpu = MAIN_SUPPORT.cpu_stat_snapshot(options.cpu)
            before_sibling = MAIN_SUPPORT.cpu_stat_snapshot(options.sibling)
            before_ns = time.monotonic_ns()
            time.sleep(5.0)
            presample = MAIN_SUPPORT.isolation_record(
                options.cpu, options.sibling, pair_lease,
                before_ns, time.monotonic_ns(), before_cpu,
                MAIN_SUPPORT.cpu_stat_snapshot(options.cpu), before_sibling,
                MAIN_SUPPORT.cpu_stat_snapshot(options.sibling))
            require(presample["delta"]["benchmark_cpu"]["nonidle_jiffies"] == 0 and
                    presample["delta"]["reserved_sibling"]["nonidle_jiffies"] == 0,
                    "CPU pair was not idle during the presample")

            for cell_index, cell in enumerate(cells):
                path = cell_artifact_path(output, cell)
                if path.exists():
                    resumed = validate_cell_artifact(
                        read_json(path, "resumed cell artifact"), cell,
                        contract_sha256, options.rounds)
                    analyses.append(resumed["analysis"])
                    print(f"{cell_index + 1}/{len(cells)} {cell['id']} resumed",
                          file=sys.stderr, flush=True)
                    continue
                rounds: list[dict[str, Any]] = []
                for round_index in range(options.rounds):
                    before_cpu = MAIN_SUPPORT.cpu_stat_snapshot(options.cpu)
                    before_sibling = \
                        MAIN_SUPPORT.cpu_stat_snapshot(options.sibling)
                    before_ns = time.monotonic_ns()
                    segments = []
                    for order in ROUND_SEGMENTS[round_index % 2]:
                        baseline = "main" if "main" in order else "control"
                        invocations = [run_one(
                            label, identities[label], cell, options.cpu,
                            options.iterations, options.warmup,
                            options.timeout) for label in order]
                        segments.append({
                            "baseline": baseline,
                            "order": list(order),
                            "invocations": invocations,
                        })
                    isolation = MAIN_SUPPORT.isolation_record(
                        options.cpu, options.sibling, pair_lease,
                        before_ns, time.monotonic_ns(), before_cpu,
                        MAIN_SUPPORT.cpu_stat_snapshot(options.cpu),
                        before_sibling,
                        MAIN_SUPPORT.cpu_stat_snapshot(options.sibling))
                    require(isolation["accepted"] and
                            isolation["delta"]["reserved_sibling"][
                                "nonidle_jiffies"] == 0,
                            f"contaminated {cell['id']} round {round_index}")
                    rounds.append({
                        "round": round_index,
                        "segments": segments,
                        "isolation": isolation,
                    })
                artifact: dict[str, Any] = {
                    "schema": CELL_SCHEMA,
                    "contract_sha256": contract_sha256,
                    "cell": dict(cell),
                    "rounds": rounds,
                    "analysis": analyze_cell(cell, rounds),
                }
                artifact["artifact_sha256"] = digest_object(artifact)
                write_json_atomic(path, artifact)
                validated = validate_cell_artifact(
                    artifact, cell, contract_sha256, options.rounds)
                analyses.append(validated["analysis"])
                print(f"{cell_index + 1}/{len(cells)} {cell['id']}",
                      file=sys.stderr, flush=True)

        verify_contract_stability(
            options, identities, provenance, dependencies)
        summary = {
            "schema": SUMMARY_SCHEMA,
            "status": "completed",
            "contract_sha256": contract_sha256,
            "cell_count": len(analyses),
            "target_cell_count": sum(
                bool(item["cell"]["candidate_selected"])
                for item in analyses),
            "control_cell_count": sum(
                not bool(item["cell"]["candidate_selected"])
                for item in analyses),
            "process_count": len(analyses) * options.rounds * 8,
            "all_workload_digests_matched": True,
            "all_rounds_accepted": True,
            "all_reserved_sibling_nonidle_jiffies_zero": True,
            "api_lane_contract": API_LANES,
            "cells": analyses,
        }
        write_json_atomic(output / "summary.json", summary)
        try:
            (output / "failure.json").unlink()
        except FileNotFoundError:
            pass
        print(json.dumps({
            "status": "completed",
            "cells": len(analyses),
            "processes": summary["process_count"],
            "summary": str(output / "summary.json"),
        }, sort_keys=True))
        return 0
    except Exception as error:
        failure_context.update({
            "failed_utc": MAIN_SUPPORT.utc_now(),
            "error_type": type(error).__name__,
            "error": str(error),
        })
        write_json_atomic(output / "failure.json", failure_context)
        print(f"evidence rejected: {error}", file=sys.stderr)
        return 1
    finally:
        os.close(lock_descriptor)


def self_test() -> int:
    validate_frozen_candidate_identity(
        CANDIDATE_COMMIT, CANDIDATE_TREE, CANDIDATE_SHA256)
    rejected_identities = (
        ("0" * 40, CANDIDATE_TREE, CANDIDATE_SHA256),
        (CANDIDATE_COMMIT, "0" * 40, CANDIDATE_SHA256),
        (CANDIDATE_COMMIT, CANDIDATE_TREE, "0" * 64),
        ("not-a-commit", CANDIDATE_TREE, CANDIDATE_SHA256),
    )
    for identity in rejected_identities:
        try:
            validate_frozen_candidate_identity(*identity)
        except EvidenceError:
            pass
        else:
            raise EvidenceError(
                "self-test accepted a non-frozen candidate identity")
    cells = campaign_cells()
    require(len(cells) == 111 and
            sum(cell["candidate_selected"] for cell in cells) == 53 and
            all(cell["bytes"] > 0 and cell["bytes"] % 64 == 0
                for cell in cells) and
            reuse_for_cell(1, 64) == MAX_REUSE and
            reuse_for_cell(2, 64) == 131072 and
            reuse_for_cell(255, 1024) == MIN_REUSE and
            all(MIN_REUSE <= cell["reuse"] <= MAX_REUSE and
                abs(cell["reuse"] * cell["K"] * cell["bytes"] -
                    TARGET_INPUT_BYTES_PER_METRIC_SAMPLE) <=
                    TARGET_INPUT_BYTES_PER_METRIC_SAMPLE // 50
                for cell in cells),
            "self-test matrix accounting failed")
    require(expected_reduction_path(2, 64, 1) == "k2_terminal" and
            expected_reduction_path(3, 256, 1) == "fused_final" and
            expected_reduction_path(4, 1024, 1) == "dense" and
            expected_reduction_path(8, 64, 1) == "coarse" and
            expected_reduction_path(8, 1024, 1) == "coarse" and
            expected_reduction_path(8, 64, 0) == "pairwise" and
            expected_reduction_path(3, 1025, 1) == "pairwise",
            "self-test reduction boundary failed")
    # Fixed vector from both C++ benchmark adapters.  This independently
    # guards their prefix-after-shuffle loss-selection convention.
    require(selected_loss(2, 90370406875137) == 1,
            "self-test benchmark loss-selection vector changed")
    require(inspect_isa_disassembly(
        "  10: c5 fd ef c0 vpxor ymm0,ymm0,ymm0") == {
            "evex_prefixed_instruction_count": 0,
            "ymm_operand_instruction_count": 1,
        } and inspect_isa_disassembly(
            "  20: 62 f1 7d 48 ef c0 vpxord zmm0,zmm0,zmm0")[
                "evex_prefixed_instruction_count"] == 1,
            "self-test ISA parser failed open")
    sample = cells[0]
    candidate = benchmark_command(
        "candidate", Path("/frozen/a/benchmark"), sample, 14, 15, 6)
    control = benchmark_command(
        "control", Path("/frozen/b/benchmark"), sample, 14, 15, 6)
    require(candidate[-4:] == [
                "--r1-small-reduction-mode", "1", "--json", "-"] and
            control[-4:] == [
                "--r1-small-reduction-mode", "0", "--json", "-"] and
            len(candidate) == len(control) and
            all(len(left) == len(right)
                for left, right in zip(candidate, control)),
            "self-test runtime attribution argv failed")
    require(API_LANES["leopard2_one_shot_encode"]["api"] == "leo2_encode" and
            API_LANES["leopard2_one_shot_decode"]["api"] == "leo2_decode" and
            API_LANES["leopard2_encode_execution"].get("not_api") ==
                "leo2_encode" and
            API_LANES["leopard2_encode_execution"]["api"] ==
                "leo2_encode_batch" and
            API_LANES["leopard2_reused_plan_execution"]["api"] ==
                "leo2_decode_plan_execute_batch",
            "self-test API lane labels became ambiguous")
    require(
        ORDINARY_BATCH_ROUTE_PROOF["required_schema"] ==
            "leopard2-benchmark-v12" and
        ORDINARY_BATCH_ROUTE_PROOF["required_batch_item_count"] == 1 and
        ORDINARY_BATCH_ROUTE_PROOF["required_loss_count"] == 1 and
        ORDINARY_BATCH_ROUTE_PROOF["forbidden_build_fields"] ==
            ["prevalidated_batch_experiment"] and
        ORDINARY_BATCH_ROUTE_PROOF["forbidden_metric_fields"] ==
            ["encode_binding_setup", "decode_binding_setup"] and
        ORDINARY_BATCH_ROUTE_PROOF["required_build_fields"] == {
            "r1_timed_encode_api":
                "leo2_encode_batch:item_count=1:no_preflight_scratch",
            "r1_timed_one_shot_encode_api": "leo2_encode",
            "r1_timed_reused_decode_api":
                "leo2_decode_plan_execute_batch:item_count=1:"
                "no_preflight_scratch:one_loss_direct_xor",
            "r1_timed_one_shot_decode_api": "leo2_decode:one_loss",
        } and
        ORDINARY_BATCH_ROUTE_PROOF["required_memory_fields"] == {
            "encode_batch_preflight_scratch_bytes": 0,
            "decode_batch_preflight_scratch_bytes": 0,
        } and
        "preserve this exact timed API scope" in
            ORDINARY_BATCH_ROUTE_PROOF["path_attestation_scope"],
        "self-test prevalidated-binding rejection contract changed")

    def setup_summary(value: float) -> dict[str, Any]:
        return {
            "median_us": value,
            "samples_us": [value] * 15,
        }

    def execution_summary(value: float) -> dict[str, Any]:
        return {
            "median_us_per_batch_call": value,
            "samples_us_per_batch_call": [value] * 15,
        }

    candidate_fixture: dict[str, Any] = {
        "schema": "leopard2-benchmark-v12",
        "parameters": {
            "K": sample["K"], "R": 1,
            "shard_bytes": sample["bytes"], "loss_count": 1,
            "missing_original_indices":
                sample["expected_missing_original_indices"],
            "batch": 1, "reuse": sample["reuse"],
            "iterations": 15, "warmup": 6, "thread_count": 1,
            "seed": sample["seed"], "measure_one_shot_decode": True,
            "skip_legacy": True, "retain_samples": True,
            "requested_profile": "legacy_high_v1",
            "requested_field": "gf8", "requested_backend": "avx2",
            "force_generic_decode": False,
            "force_specialized_decode": False,
            "force_tiled_decode": False,
            "force_materialized_decode": False,
        },
        "resolved": {
            "profile": "legacy_high_v1", "field": "gf8",
            "backend": "avx2", "thread_count": 1, "padded_side": 1,
        },
        "correctness": {"leopard2_round_trip": True},
        "workload_digests": {
            "algorithm": "fnv1a64", "original_data": "0" * 16,
            "transmitted_parity": "1" * 16,
            "recovered_originals": "2" * 16,
        },
        "build": {
            "r1_small_reduction_diagnostic_mode": 1,
            "r1_small_reduction_codec_enabled": True,
            "r1_encode_reduction_path":
                sample["candidate_reduction_path"],
            "r1_decode_reduction_path":
                sample["candidate_reduction_path"],
            **ORDINARY_BATCH_ROUTE_PROOF["required_build_fields"],
        },
        "memory": dict(ORDINARY_BATCH_ROUTE_PROOF[
            "required_memory_fields"]),
        "metrics": {
            "codec_setup": setup_summary(1.0),
            "encode_execution": execution_summary(2.0),
            "one_shot_encode": execution_summary(2.5),
            "decode_plan_setup": setup_summary(3.0),
            "decode_execution": execution_summary(4.0),
            "decode_amortized_at_reuse": {
                "reuse_count": sample["reuse"],
                "derived_median_us_per_batch_call":
                    4.0 + 3.0 / sample["reuse"],
            },
            "one_shot_decode_including_setup": execution_summary(7.0),
        },
    }
    normalized_fixture = validate_result(
        "candidate", candidate_fixture, sample, 15, 6)
    require(normalized_fixture["api_lanes"]["one_shot_encode"]["api"] ==
            "leo2_encode" and
            normalized_fixture["api_lanes"]["one_shot_decode"]["api"] ==
                "leo2_decode" and
            normalized_fixture["reduction_path"] == "k2_terminal",
            "self-test ordinary schema-v12 fixture was mislabeled")
    main_fixture: dict[str, Any] = {
        "schema": "leopard-main-benchmark-v1",
        "parameters": {
            "K": sample["K"], "R": 1,
            "shard_bytes": sample["bytes"], "loss_count": 1,
            "missing_original_indices":
                sample["expected_missing_original_indices"],
            "batch": 1, "reuse": sample["reuse"],
            "iterations": 15, "warmup": 6, "thread_count": 1,
            "seed": sample["seed"],
        },
        "resolved": {
            "profile": "legacy_high_v1", "field": "gf8",
            "thread_count": 1, "padded_side": 1,
        },
        "correctness": {"round_trip": True},
        "workload_digests": dict(candidate_fixture["workload_digests"]),
        "build": {"main_source_commit": MAIN_COMMIT},
        "metrics": {
            "codec_setup": None,
            "decode_timing_includes_setup": True,
            "encode_execution": execution_summary(4.0),
            "decode_including_setup": execution_summary(14.0),
        },
    }
    normalized_main = validate_result(
        "main", main_fixture, sample, 15, 6)
    require(normalized_main["api_lanes"]["legacy_encode"]["api"] ==
                "leo_encode" and
            normalized_main["api_lanes"]["legacy_decode"]["api"] ==
                "leo_decode",
            "self-test exact-main API lanes were mislabeled")

    normalized_control = json.loads(json.dumps(normalized_fixture))
    normalized_by_implementation = {
        "candidate": normalized_fixture,
        "control": normalized_control,
        "main": normalized_main,
    }
    synthetic_rounds: list[dict[str, Any]] = []
    for round_index in range(3):
        synthetic_segments = []
        for order in ROUND_SEGMENTS[round_index % 2]:
            synthetic_segments.append({
                "baseline": "main" if "main" in order else "control",
                "order": list(order),
                "invocations": [
                    {
                        "implementation": implementation,
                        "normalized": normalized_by_implementation[
                            implementation],
                    }
                    for implementation in order
                ],
            })
        synthetic_rounds.append({
            "round": round_index,
            "isolation": {
                "accepted": True,
                "delta": {"reserved_sibling": {"nonidle_jiffies": 0}},
            },
            "segments": synthetic_segments,
        })
    synthetic_analysis = analyze_cell(sample, synthetic_rounds)
    require(
        set(synthetic_analysis["candidate_vs_main"]) == {
            "legacy_encode_vs_batch_encode",
            "legacy_encode_vs_one_shot_encode",
            "legacy_decode_vs_one_shot_decode",
        } and
        math.isclose(synthetic_analysis["candidate_vs_main"][
            "legacy_encode_vs_batch_encode"]["speedup"], 2.0) and
        math.isclose(synthetic_analysis["candidate_vs_main"][
            "legacy_encode_vs_one_shot_encode"]["speedup"], 1.6) and
        math.isclose(synthetic_analysis["candidate_vs_main"][
            "legacy_decode_vs_one_shot_decode"]["speedup"], 2.0),
        "self-test exact-main comparison mapping changed")
    for mutation in ("round_number", "serialized_order", "invocation_order"):
        changed_rounds = json.loads(json.dumps(synthetic_rounds))
        if mutation == "round_number":
            changed_rounds[1]["round"] = 9
        elif mutation == "serialized_order":
            changed_rounds[0]["segments"][0]["order"][0:2] = \
                reversed(changed_rounds[0]["segments"][0]["order"][0:2])
        else:
            invocations = changed_rounds[0]["segments"][0]["invocations"]
            invocations[0], invocations[1] = invocations[1], invocations[0]
        rejected = False
        try:
            analyze_cell(sample, changed_rounds)
        except EvidenceError:
            rejected = True
        require(rejected,
                f"self-test accepted malformed resume order: {mutation}")
    for mutation in (
            "prevalidated_marker", "binding_metric",
            "missing_one_shot_encode", "missing_one_shot_decode",
            "encode_api_scope", "one_shot_encode_api_scope",
            "decode_api_scope", "one_shot_decode_api_scope",
            "encode_preflight", "decode_preflight"):
        changed = json.loads(json.dumps(candidate_fixture))
        if mutation == "prevalidated_marker":
            changed["build"]["prevalidated_batch_experiment"] = True
        elif mutation == "binding_metric":
            changed["metrics"]["encode_binding_setup"] = setup_summary(1.0)
        elif mutation == "missing_one_shot_encode":
            changed["metrics"].pop("one_shot_encode")
        elif mutation == "missing_one_shot_decode":
            changed["metrics"].pop("one_shot_decode_including_setup")
        elif mutation == "encode_api_scope":
            changed["build"]["r1_timed_encode_api"] = "leo2_encode"
        elif mutation == "one_shot_encode_api_scope":
            changed["build"]["r1_timed_one_shot_encode_api"] = \
                "leo2_encode_batch:item_count=1"
        elif mutation == "decode_api_scope":
            changed["build"]["r1_timed_reused_decode_api"] = \
                "leo2_decode_plan_execute"
        elif mutation == "one_shot_decode_api_scope":
            changed["build"]["r1_timed_one_shot_decode_api"] = \
                "leo2_decode_plan_execute_batch:item_count=1"
        elif mutation == "encode_preflight":
            changed["memory"]["encode_batch_preflight_scratch_bytes"] = 64
        else:
            changed["memory"]["decode_batch_preflight_scratch_bytes"] = 64
        rejected = False
        try:
            validate_result("candidate", changed, sample, 15, 6)
        except EvidenceError:
            rejected = True
        require(rejected,
                f"self-test accepted forbidden route mutation: {mutation}")

    with tempfile.TemporaryDirectory(prefix="leo2-r1-runner-selftest-") as raw:
        directory = Path(raw)
        shared = directory / "aaaaaaaa"
        control_link = directory / "bbbbbbbb"
        main = directory / "mmmmmmmm"
        shared.write_bytes(b"candidate")
        os.link(shared, control_link)
        main.write_bytes(b"main")
        identities = {name: file_identity(path) for name, path in (
            ("candidate", shared), ("control", control_link),
            ("main", main))}
        proof = validate_shared_binary_identities(
            identities, sha256(shared), sha256(main))
        require(proof["candidate_control_inode"] == shared.stat().st_ino,
                "self-test shared-inode proof failed")
        candidate_build = directory / "candidate-build"
        main_build = directory / "main-build"
        candidate_build.write_bytes(shared.read_bytes())
        main_build.write_bytes(main.read_bytes())
        provenance_fixture = {
            "candidate_build": {
                "executable": file_identity(candidate_build),
            },
            "main_build": {
                "validated_executable": file_identity(main_build),
            },
        }
        copy_proof = validate_lane_owned_copy_identities(
            identities, provenance_fixture)
        require(copy_proof["lane_inodes_distinct_from_build_artifacts"] is True,
                "self-test lane-owned copy proof failed")
        for direct_lane in ("candidate", "main"):
            direct_provenance = json.loads(json.dumps(provenance_fixture))
            if direct_lane == "candidate":
                direct_provenance["candidate_build"]["executable"] = \
                    dict(identities["candidate"])
            else:
                direct_provenance["main_build"]["validated_executable"] = \
                    dict(identities["main"])
            rejected = False
            try:
                validate_lane_owned_copy_identities(
                    identities, direct_provenance)
            except EvidenceError:
                rejected = True
            require(rejected,
                    f"self-test accepted direct {direct_lane} build hardlink")
        copied = directory / "cccccccc"
        copied.write_bytes(shared.read_bytes())
        bad = dict(identities)
        bad["control"] = file_identity(copied)
        rejected = False
        try:
            validate_shared_binary_identities(
                bad, sha256(shared), sha256(main))
        except EvidenceError:
            rejected = True
        require(rejected, "self-test accepted equal bytes on distinct inodes")

        artifact = {
            "schema": CELL_SCHEMA,
            "contract_sha256": "0" * 64,
            "cell": sample,
            "rounds": [],
            "analysis": {},
        }
        artifact["artifact_sha256"] = digest_object(artifact)
        rejected = False
        try:
            validate_cell_artifact(artifact, sample, "1" * 64, 3)
        except EvidenceError:
            rejected = True
        require(rejected, "self-test accepted a foreign resume contract")
    print(json.dumps({"status": "self-test-passed", "cells": len(cells)}))
    return 0


def parse_arguments(arguments: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--control", type=Path)
    parser.add_argument("--main", type=Path)
    parser.add_argument("--candidate-sha256")
    parser.add_argument("--main-sha256", default=MAIN_SHA256)
    parser.add_argument("--build-dir", type=Path,
                        help="provenance root containing bench_leopard2")
    parser.add_argument("--source-root", type=Path)
    parser.add_argument("--source-commit")
    parser.add_argument("--source-tree")
    parser.add_argument("--main-build-dir", type=Path,
                        help="provenance root containing leopard_main_benchmark")
    parser.add_argument("--main-source-root", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--cpu", type=int, default=BENCHMARK_CPU)
    parser.add_argument("--sibling", type=int, default=RESERVED_SIBLING)
    parser.add_argument("--rounds", type=int, choices=tuple(T95), default=3)
    parser.add_argument("--iterations", type=int, default=15)
    parser.add_argument("--warmup", type=int, default=6)
    parser.add_argument("--timeout", type=float, default=60.0)
    options = parser.parse_args(arguments)
    if not options.self_test:
        required = {
            "candidate": options.candidate,
            "control": options.control,
            "main": options.main,
            "candidate-sha256": options.candidate_sha256,
            "build-dir": options.build_dir,
            "source-root": options.source_root,
            "source-commit": options.source_commit,
            "source-tree": options.source_tree,
            "main-build-dir": options.main_build_dir,
            "main-source-root": options.main_source_root,
            "output": options.output,
        }
        missing = sorted(name for name, value in required.items()
                         if value is None)
        require(not missing, "missing required arguments: " + ", ".join(missing))
    return options


def main(arguments: Sequence[str] | None = None) -> int:
    options = parse_arguments(arguments)
    return self_test() if options.self_test else run_campaign(options)


if __name__ == "__main__":
    raise SystemExit(main())
