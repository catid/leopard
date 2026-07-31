#!/usr/bin/env python3
"""Measure the two-block T=8 encoder candidate against control and main.

This runner is deliberately narrow.  It consumes immutable, already-built
executables, runs the reusable prevalidated-batch API on one logical CPU, and
rejects a round if its reserved SMT sibling performs non-idle work.  The
candidate and control must report the same clean source commit and tree; exact
Leopard main must report the pinned historical commit.
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
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA = "leopard2-gf8-t8-two-block-abba/v1"
SUMMARY_SCHEMA = "leopard2-gf8-t8-two-block-summary/v1"
ONE_KIB_SCHEMA = "leopard2-gf8-t8-one-kib-extension-abba/v1"
ONE_KIB_SUMMARY_SCHEMA = \
    "leopard2-gf8-t8-one-kib-extension-summary/v1"
TINY_SCHEMA = "leopard2-gf8-t8-tiny-extension-abba/v4"
TINY_SUMMARY_SCHEMA = "leopard2-gf8-t8-tiny-extension-summary/v4"
RAGGED_SCHEMA = "leopard2-gf8-t8-ragged-extension-abba/v1"
RAGGED_SUMMARY_SCHEMA = "leopard2-gf8-t8-ragged-extension-summary/v1"
MAIN_COMMIT = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
T95_DF2 = 4.302652729911275
T95_DF8 = 2.306004135204166
TARGET_CONTROL_FLOOR = 1.05
TARGET_MAIN_FLOOR = 1.0
NEIGHBOR_FLOOR = 1.0 / 1.02
LOCK_PATH = Path("/tmp/leopard-gf8-authoritative.lock")
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
EXTENDED_PRODUCTION_MASKS = {
    384: 0xFFFFFFF6,
    448: 0xFFFF5FF4,
    512: 0xFFFFEFF0,
    576: 0xFFFF3FF0,
    640: 0xFFFF1FF0,
    704: 0xFFFF2F60,
    768: 0x6FFF0E70,
    832: 0x5FFF0D80,
    896: 0xFFFF0FD0,
    960: 0x5FFF0D40,
    1024: 0x6FF70C00,
}
ONE_KIB_ONE_BLOCK_PRIOR_MASK = 0x4FCC
ONE_KIB_ONE_BLOCK_EXTENSION_MASK = 0x0030
ONE_KIB_TWO_BLOCK_PRIOR_MASK = EXTENDED_PRODUCTION_MASKS[1024]
ONE_KIB_TWO_BLOCK_EXTENSION_MASK = 0x10000080
TINY_BYTE_COUNTS = (1, 2, 3, 7, 8, 15, 16, 17, 31, 32, 33, 63)
TINY_NEIGHBOR_BYTE_COUNTS = (64, 65)
RAGGED_TARGET_BYTE_COUNTS = (
    65, 66, 95, 96, 97, 127, 129, 191, 193, 224,
    257, 319, 321, 352, 416, 449, 480, 513, 544,
    577, 608, 641, 672, 736, 769, 800, 864, 897, 928,
)
RAGGED_NEIGHBOR_BYTE_COUNTS = (
    225, 255, 353, 385, 417, 481, 511, 545, 609, 673,
    705, 737, 801, 833, 865, 929, 961, 992, 993, 1023,
)
PADDED_MAX_ISOLATION_ATTEMPTS = 3


class EvidenceError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvidenceError(message)


def load_support() -> Any:
    path = Path(__file__).resolve().parents[1] / "main_compare" / "run_abba.py"
    specification = importlib.util.spec_from_file_location(
        "t8_two_block_main_compare_support", path)
    require(specification is not None and specification.loader is not None,
            f"cannot load evidence support: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


SUPPORT = load_support()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_exclusive(path: Path, value: object) -> None:
    payload = json.dumps(
        value, indent=2, sort_keys=True, allow_nan=False
    ).encode("utf-8") + b"\n"
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    try:
        require(os.write(descriptor, payload) == len(payload),
                f"short write: {path}")
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def write_bytes_exclusive(path: Path, payload: bytes) -> None:
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    try:
        require(os.write(descriptor, payload) == len(payload),
                f"short write: {path}")
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def file_identity(path: Path) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    status = resolved.stat()
    require(status.st_size > 0 and os.access(resolved, os.X_OK),
            f"benchmark is not executable: {resolved}")
    return {
        "path": str(resolved),
        "size": status.st_size,
        "mode": status.st_mode & 0o777,
        "sha256": sha256(resolved),
    }


def regular_file_identity(path: Path) -> dict[str, Any]:
    resolved = path.resolve(strict=True)
    status = resolved.stat()
    require(status.st_size > 0 and resolved.is_file(),
            f"evidence input is not a nonempty regular file: {resolved}")
    return {
        "path": str(resolved),
        "size": status.st_size,
        "mode": status.st_mode & 0o777,
        "sha256": sha256(resolved),
    }


def executable_sections_identity(executable: Path) -> dict[str, Any]:
    """Hash every ELF section carrying executable instructions."""
    readelf_name = shutil.which("readelf")
    objcopy_name = shutil.which("objcopy")
    require(readelf_name is not None and objcopy_name is not None,
            "readelf and objcopy are required for executable provenance")
    readelf = Path(readelf_name).resolve(strict=True)
    objcopy = Path(objcopy_name).resolve(strict=True)
    section_table = subprocess.run(
        [str(readelf), "-W", "-S", str(executable)],
        env=CHILD_ENVIRONMENT, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, timeout=10.0, check=False)
    require(section_table.returncode == 0,
            "readelf failed: " +
            section_table.stderr.decode(
                "utf-8", errors="replace")[-1000:])
    section_pattern = re.compile(
        r"^\s*\[\s*\d+\]\s+(\.\S+)\s+\S+\s+"
        r"[0-9A-Fa-f]+\s+[0-9A-Fa-f]+\s+[0-9A-Fa-f]+\s+"
        r"\S+\s+([A-Z]*)\s+")
    section_names: list[str] = []
    for line in section_table.stdout.decode(
            "utf-8", errors="strict").splitlines():
        match = section_pattern.match(line)
        if match and "X" in match.group(2):
            section_names.append(match.group(1))
    require(".text" in section_names and section_names,
            "ELF executable sections are incomplete")

    records: list[dict[str, Any]] = []
    combined = hashlib.sha256()
    with tempfile.TemporaryDirectory(
            prefix="leopard-t8-executable-sections-") as directory:
        root = Path(directory)
        for index, section_name in enumerate(section_names):
            output = root / f"section-{index}.bin"
            copied_elf = root / f"copy-{index}.elf"
            completed = subprocess.run(
                [str(objcopy), "--dump-section",
                 f"{section_name}={output}", str(executable),
                 str(copied_elf)],
                env=CHILD_ENVIRONMENT, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, timeout=10.0, check=False)
            require(completed.returncode == 0 and output.is_file() and
                    copied_elf.is_file(),
                    f"objcopy failed for {section_name}: " +
                    completed.stderr.decode(
                        "utf-8", errors="replace")[-1000:])
            payload = output.read_bytes()
            records.append({
                "name": section_name,
                "size": len(payload),
                "sha256": hashlib.sha256(payload).hexdigest(),
            })
            combined.update(section_name.encode("utf-8"))
            combined.update(b"\0")
            combined.update(payload)
    return {
        "sections": records,
        "combined_sha256": combined.hexdigest(),
        "readelf": file_identity(readelf),
        "objcopy": file_identity(objcopy),
    }


def production_selected(k: int, r: int, shard_bytes: int) -> bool:
    if not (9 <= k <= 16 and 5 <= r <= 8):
        return False
    if shard_bytes in (64, 128, 192, 256, 320):
        return True
    mask = EXTENDED_PRODUCTION_MASKS.get(shard_bytes)
    if mask is None:
        return False
    bit = 4 * (k - 9) + (r - 5)
    return (mask & (1 << bit)) != 0


def one_kib_selection(
    k: int,
    r: int,
    extension_enabled: bool,
) -> bool:
    if 5 <= k <= 8 and 5 <= r <= min(k, 8):
        bit = 4 * (k - 5) + (r - 5)
        mask = ONE_KIB_ONE_BLOCK_PRIOR_MASK
        if extension_enabled:
            mask |= ONE_KIB_ONE_BLOCK_EXTENSION_MASK
        return (mask & (1 << bit)) != 0
    if 9 <= k <= 16 and 5 <= r <= 8:
        bit = 4 * (k - 9) + (r - 5)
        mask = ONE_KIB_TWO_BLOCK_PRIOR_MASK
        if extension_enabled:
            mask |= ONE_KIB_TWO_BLOCK_EXTENSION_MASK
        return (mask & (1 << bit)) != 0
    return False


def one_kib_extension_cells() -> list[dict[str, Any]]:
    cells: list[dict[str, Any]] = []
    index = 0
    for k in range(5, 17):
        for r in range(5, min(k, 8) + 1):
            control_selected = one_kib_selection(k, r, False)
            candidate_selected = one_kib_selection(k, r, True)
            cells.append({
                "id": f"one-kib-k{k}-r{r}",
                "K": k,
                "R": r,
                "bytes": 1024,
                "role": (
                    "target"
                    if candidate_selected and not control_selected
                    else "neighbor"
                ),
                "seed": 0x172E000 + index,
                "candidate_selected": candidate_selected,
                "control_selected": control_selected,
            })
            index += 1
    require(len(cells) == 42, "one-kib shape matrix is incomplete")
    require(sum(cell["role"] == "target" for cell in cells) == 4,
            "one-kib target intersection changed")
    return cells


def tiny_extension_cells() -> list[dict[str, Any]]:
    cells: list[dict[str, Any]] = []
    index = 0
    for k in range(5, 17):
        for r in range(5, min(k, 8) + 1):
            for shard_bytes in TINY_BYTE_COUNTS:
                cells.append({
                    "id": f"tiny-k{k}-r{r}-b{shard_bytes}",
                    "K": k,
                    "R": r,
                    "bytes": shard_bytes,
                    "role": "target",
                    "seed": 0x182E000 + index,
                    "candidate_selected": True,
                    "control_selected": False,
                    "main_physical_shard_bytes": 64,
                })
                index += 1
            for shard_bytes in TINY_NEIGHBOR_BYTE_COUNTS:
                selected = shard_bytes == 64
                cells.append({
                    "id": f"tiny-neighbor-k{k}-r{r}-b{shard_bytes}",
                    "K": k,
                    "R": r,
                    "bytes": shard_bytes,
                    "role": "neighbor",
                    "seed": 0x182E000 + index,
                    "candidate_selected": selected,
                    "control_selected": selected,
                })
                index += 1
    require(len(cells) == 588, "tiny T=8 matrix is incomplete")
    require(sum(cell["role"] == "target" for cell in cells) == 504,
            "tiny T=8 target matrix is incomplete")
    return cells


def ragged_selected(shard_bytes: int) -> bool:
    if shard_bytes < 65 or shard_bytes > 1024 or shard_bytes % 64 == 0:
        return False
    return (
        65 <= shard_bytes <= 191 or
        193 <= shard_bytes <= 224 or
        257 <= shard_bytes <= 352 or
        shard_bytes == 416 or
        449 <= shard_bytes <= 480 or
        513 <= shard_bytes <= 544 or
        577 <= shard_bytes <= 608 or
        641 <= shard_bytes <= 672 or
        shard_bytes == 736 or
        769 <= shard_bytes <= 800 or
        shard_bytes == 864 or
        897 <= shard_bytes <= 928
    )


def ragged_extension_cells() -> list[dict[str, Any]]:
    cells: list[dict[str, Any]] = []
    index = 0
    for k in range(5, 17):
        for r in range(5, min(k, 8) + 1):
            for role, byte_counts in (
                ("target", RAGGED_TARGET_BYTE_COUNTS),
                ("neighbor", RAGGED_NEIGHBOR_BYTE_COUNTS),
            ):
                for shard_bytes in byte_counts:
                    selected = ragged_selected(shard_bytes)
                    require(selected is (role == "target"),
                            "ragged target/neighbor table is inconsistent")
                    cells.append({
                        "id": f"ragged-{role}-k{k}-r{r}-b{shard_bytes}",
                        "K": k,
                        "R": r,
                        "bytes": shard_bytes,
                        "role": role,
                        "seed": 0x192E000 + index,
                        "candidate_selected": selected,
                        "control_selected": False,
                        "main_physical_shard_bytes":
                            ((shard_bytes + 63) // 64) * 64,
                    })
                    index += 1
    expected_targets = 42 * len(RAGGED_TARGET_BYTE_COUNTS)
    expected_neighbors = 42 * len(RAGGED_NEIGHBOR_BYTE_COUNTS)
    require(len(cells) == expected_targets + expected_neighbors,
            "ragged T=8 matrix is incomplete")
    require(sum(cell["role"] == "target" for cell in cells) ==
            expected_targets,
            "ragged T=8 target matrix is incomplete")
    return cells


def padded_campaign(campaign: str) -> bool:
    return campaign in ("tiny-extension", "ragged-extension")


def target_cells(
    target_bytes: int = 64,
    final_selector: bool = False,
) -> list[dict[str, Any]]:
    cells = []
    index = 0
    for k in range(9, 17):
        for r in range(5, 9):
            role = "target"
            if final_selector and not production_selected(
                    k, r, target_bytes):
                role = "excluded_neighbor"
            cells.append({
                "id": f"target-k{k}-r{r}-b{target_bytes}",
                "K": k,
                "R": r,
                "bytes": target_bytes,
                "role": role,
                "seed": (
                    0x152E000 if final_selector else 0x142E000
                ) + index,
            })
            index += 1
    return cells


def neighbor_cells(
    target_bytes: int = 64,
    final_selector: bool = False,
) -> list[dict[str, Any]]:
    cells = []
    byte_neighbors_by_target = {
        64: (32, 33, 63, 65),
        128: (64, 96, 127, 129, 160, 256, 320),
        192: (
            64, 128, 160, 191, 193, 224,
            255, 256, 257, 320, 384, 512,
        ),
        256: (64, 65, 255, 257, 1024),
        320: (256, 257, 288, 319, 321, 352, 384, 512),
    }
    if target_bytes >= 384 and target_bytes <= 1024 and \
            target_bytes % 64 == 0:
        byte_neighbors = tuple(sorted({
            320, target_bytes - 64, target_bytes - 1,
            target_bytes + 1, target_bytes + 64,
        }))
    else:
        require(target_bytes in byte_neighbors_by_target,
                f"unsupported target byte count: {target_bytes}")
        byte_neighbors = byte_neighbors_by_target[target_bytes]
    for k, r in ((9, 5), (16, 8)):
        for shard_bytes in byte_neighbors:
            cells.append({
                "id": f"bytes-k{k}-r{r}-b{shard_bytes}",
                "K": k,
                "R": r,
                "bytes": shard_bytes,
                "role": "byte_neighbor",
            })
    for k, r in (
        (8, 5), (8, 8), (17, 5), (17, 8),
        (9, 4), (16, 4), (9, 9), (16, 9),
    ):
        cells.append({
            "id": f"shape-k{k}-r{r}-b{target_bytes}",
            "K": k,
            "R": r,
            "bytes": target_bytes,
            "role": "shape_neighbor",
        })
    for index, cell in enumerate(cells):
        cell["seed"] = (
            0x152F000 if final_selector else 0x142F000
        ) + index
    return cells


def benchmark_command(
    implementation: str,
    executable: Path,
    cell: Mapping[str, Any],
    cpu: int,
    iterations: int,
    warmup: int,
    campaign: str = "two-block",
) -> list[str]:
    shard_bytes = (
        int(cell["main_physical_shard_bytes"])
        if implementation == "main" and padded_campaign(campaign)
        else int(cell["bytes"])
    )
    common = [
        "/usr/bin/prlimit", "--as=201326592",
        "/usr/bin/taskset", "-c", str(cpu), str(executable),
        "--k", str(cell["K"]), "--r", str(cell["R"]),
        "--bytes", str(shard_bytes), "--loss", "1",
        "--batch", "64", "--reuse", "64",
        "--iterations", str(iterations), "--warmup", str(warmup),
        "--threads", "1", "--seed", str(cell["seed"]), "--json", "-",
    ]
    if implementation == "main":
        if padded_campaign(campaign):
            common[-2:-2] = [
                "--logical-bytes", str(cell["bytes"]),
            ]
        return common
    return common[:-2] + [
        "--profile", "high", "--field", "gf8", "--backend", "avx2",
        "--skip-legacy", "--retain-samples", "--attest-source",
        "--json", "-",
    ]


def positive_metric(result: Mapping[str, Any], name: str) -> float:
    metrics = result.get("metrics")
    require(isinstance(metrics, dict), "benchmark metrics are absent")
    value = metrics.get(name)
    require(isinstance(value, dict), f"missing benchmark metric: {name}")
    median = value.get("median_us_per_batch_call")
    require(isinstance(median, (int, float)) and
            not isinstance(median, bool) and
            math.isfinite(float(median)) and float(median) > 0,
            f"invalid benchmark metric: {name}")
    return float(median)


def validate_result(
    implementation: str,
    result: object,
    cell: Mapping[str, Any],
    source_commit: str,
    source_tree: str,
    iterations: int,
    warmup: int,
    target_bytes: int,
    campaign: str = "two-block",
    final_selector: bool = False,
) -> dict[str, Any]:
    require(isinstance(result, dict), "benchmark output is not an object")
    expected_shard_bytes = (
        int(cell["main_physical_shard_bytes"])
        if implementation == "main" and padded_campaign(campaign)
        else int(cell["bytes"])
    )
    expected_parameters = {
        "K": cell["K"], "R": cell["R"],
        "shard_bytes": expected_shard_bytes,
        "loss_count": 1, "batch": 64, "reuse": 64,
        "iterations": iterations, "warmup": warmup,
        "thread_count": 1, "seed": cell["seed"],
    }
    if implementation == "main" and padded_campaign(campaign):
        expected_parameters["logical_shard_bytes"] = int(cell["bytes"])
    parameters = result.get("parameters")
    require(isinstance(parameters, dict) and
            all(parameters.get(name) == value
                for name, value in expected_parameters.items()),
            "benchmark parameters differ from the frozen cell")
    resolved = result.get("resolved")
    correctness = result.get("correctness")
    digests = result.get("workload_digests")
    require(isinstance(resolved, dict) and
            resolved.get("profile") == "legacy_high_v1" and
            resolved.get("field") == "gf8" and
            resolved.get("thread_count") == 1,
            "benchmark resolved a different profile, field, or thread count")
    require(isinstance(digests, dict) and digests.get("algorithm") == "fnv1a64" and
            all(isinstance(digests.get(name), str) and len(digests[name]) == 16
                for name in (
                    "original_data", "transmitted_parity",
                    "recovered_originals")),
            "benchmark digests are incomplete")
    if implementation == "main":
        expected_schema = (
            "leopard-main-benchmark-v2"
            if padded_campaign(campaign)
            else "leopard-main-benchmark-v1"
        )
        require(result.get("schema") == expected_schema and
                isinstance(correctness, dict) and
                correctness.get("round_trip") is True,
                "exact-main benchmark identity or round trip failed")
        build = result.get("build")
        require(isinstance(build, dict) and
                build.get("main_source_commit") == MAIN_COMMIT,
                "exact-main source identity changed")
        if padded_campaign(campaign):
            require(
                resolved.get("padded_application_bytes") is True and
                resolved.get("padding_policy") == "zero suffix per shard" and
                correctness.get("logical_prefix_fingerprinted") is True,
                "exact-main padded comparison semantics changed")
    else:
        require(result.get("schema") == "leopard2-benchmark-v5" and
                resolved.get("backend") == "avx2" and
                isinstance(correctness, dict) and
                correctness.get("leopard2_round_trip") is True,
                "Leopard2 benchmark identity or round trip failed")
        build = result.get("build")
        require(campaign in (
                    "two-block",
                    "one-block-extended",
                    "one-block-beyond512",
                    "one-kib-extension",
                    "tiny-extension",
                    "ragged-extension",
                ),
                f"unsupported campaign identity: {campaign}")
        require(isinstance(build, dict),
                "Leopard2 build identity is absent")
        if campaign == "tiny-extension":
            expected_marker = implementation == "candidate"
            expected_selected = cell.get(
                f"{implementation}_selected") is True
            one_block_shape = int(cell["K"]) <= 8
            marker_valid = (
                build.get("high_t8_tiny_binding_enabled") is expected_marker and
                build.get("high_t8_one_block_extended_enabled") is True and
                build.get("high_t8_one_block_beyond_512_enabled") is True and
                build.get("high_t8_one_kilobyte_extension_enabled") is True and
                build.get("high_t8_two_block_128_192_enabled") is True and
                build.get("high_t8_two_block_320_enabled") is True and
                build.get("high_t8_two_block_extended_enabled") is True and
                build.get("high_t8_one_block_selected") is
                    (expected_selected and one_block_shape) and
                build.get("high_t8_two_block_selected") is
                    (expected_selected and not one_block_shape)
            )
        elif campaign == "ragged-extension":
            expected_marker = implementation == "candidate"
            expected_selected = cell.get(
                f"{implementation}_selected") is True
            one_block_shape = int(cell["K"]) <= 8
            marker_valid = (
                build.get("high_t8_ragged_binding_enabled") is
                    expected_marker and
                build.get("high_t8_tiny_binding_enabled") is True and
                build.get("high_t8_one_block_extended_enabled") is True and
                build.get("high_t8_one_block_beyond_512_enabled") is True and
                build.get("high_t8_one_kilobyte_extension_enabled") is True and
                build.get("high_t8_two_block_128_192_enabled") is True and
                build.get("high_t8_two_block_320_enabled") is True and
                build.get("high_t8_two_block_extended_enabled") is True and
                build.get("high_t8_one_block_selected") is
                    (expected_selected and one_block_shape) and
                build.get("high_t8_two_block_selected") is
                    (expected_selected and not one_block_shape)
            )
        elif campaign == "one-kib-extension":
            expected_extension = implementation == "candidate"
            expected_selected = cell.get(
                f"{implementation}_selected") is True
            one_block_shape = int(cell["K"]) <= 8
            marker_valid = (
                build.get("high_t8_one_block_extended_enabled") is True and
                build.get("high_t8_one_block_beyond_512_enabled") is True and
                build.get("high_t8_one_kilobyte_extension_enabled") is
                    expected_extension and
                build.get("high_t8_two_block_128_192_enabled") is True and
                build.get("high_t8_two_block_320_enabled") is True and
                build.get("high_t8_two_block_extended_enabled") is True and
                build.get("high_t8_one_block_selected") is
                    (expected_selected and one_block_shape) and
                build.get("high_t8_two_block_selected") is
                    (expected_selected and not one_block_shape)
            )
        elif campaign == "two-block":
            expected_128_192 = (
                implementation == "candidate" or
                target_bytes not in (128, 192)
            )
            expected_320 = (
                implementation == "candidate" or target_bytes != 320
            )
            expected_extended = (
                implementation == "candidate" or target_bytes < 384
            )
            eligible_shape = (
                9 <= cell["K"] <= 16 and 5 <= cell["R"] <= 8
            )
            byte_count = int(cell["bytes"])
            if final_selector and implementation == "candidate":
                expected_selected = production_selected(
                    int(cell["K"]), int(cell["R"]), byte_count)
            else:
                expected_selected = eligible_shape and (
                    byte_count == 64 or
                    (byte_count in (128, 192) and expected_128_192) or
                    byte_count == 256 or
                    (byte_count == 320 and expected_320) or
                    (384 <= byte_count <= 1024 and
                     byte_count % 64 == 0 and expected_extended)
                )
            marker_valid = (
                build.get("high_t8_two_block_128_192_enabled") is
                    expected_128_192 and
                build.get("high_t8_two_block_320_enabled") is expected_320 and
                build.get("high_t8_two_block_extended_enabled") is
                    expected_extended and
                build.get("high_t8_two_block_selected") is expected_selected
            )
        elif campaign == "one-block-extended":
            expected_selected = (
                (str(cell["role"]).startswith("target") and
                 implementation == "candidate") or
                (cell["role"] == "neighbor" and cell["bytes"] == 64 and
                 5 <= cell["K"] <= 8 and 5 <= cell["R"] <= 8)
            )
            marker_valid = (
                build.get("high_t8_one_block_extended_enabled") is
                    (implementation == "candidate") and
                build.get("high_t8_one_block_selected") is
                    expected_selected and
                build.get("high_t8_two_block_128_192_enabled") is True and
                build.get("high_t8_two_block_320_enabled") is True
            )
        else:
            byte_count = int(cell["bytes"])
            eligible_shape = (
                5 <= cell["K"] <= 8 and 5 <= cell["R"] <= 8
            )
            expected_beyond = implementation == "candidate"
            if final_selector and implementation == "candidate":
                expected_selected = cell.get(
                    "candidate_selected") is True
            else:
                expected_selected = eligible_shape and (
                    byte_count == 64 or
                    (128 <= byte_count <= 512 and
                     byte_count % 64 == 0) or
                    (576 <= byte_count <= 1024 and
                     byte_count % 64 == 0 and expected_beyond)
                )
            marker_valid = (
                build.get("high_t8_one_block_extended_enabled") is True and
                build.get("high_t8_one_block_beyond_512_enabled") is
                    expected_beyond and
                build.get("high_t8_one_block_selected") is
                    expected_selected and
                build.get("high_t8_two_block_128_192_enabled") is True and
                build.get("high_t8_two_block_320_enabled") is True
            )
        require(build.get("prevalidated_batch_experiment") is True and
                marker_valid and
                build.get("source_commit") == source_commit and
                build.get("source_tree") == source_tree and
                build.get("source_tracked_dirty") is False,
                "Leopard2 build mode or embedded source identity changed")
    return {
        "encode_us": positive_metric(result, "encode_execution"),
        "digests": dict(digests),
    }


def run_one(
    implementation: str,
    identity: Mapping[str, Any],
    cell: Mapping[str, Any],
    cpu: int,
    source_commit: str,
    source_tree: str,
    iterations: int,
    warmup: int,
    target_bytes: int,
    campaign: str = "two-block",
    final_selector: bool = False,
    failure_output: Path | None = None,
) -> dict[str, Any]:
    executable = Path(str(identity["path"]))
    require(sha256(executable) == identity["sha256"],
            f"{implementation} binary changed before execution")
    command = benchmark_command(
        implementation, executable, cell, cpu, iterations, warmup, campaign)
    start = time.monotonic_ns()
    failure_prefix = (
        failure_output /
        f"failure-{implementation}-{time.monotonic_ns()}"
        if failure_output is not None else None
    )

    def persist_failure(stdout: bytes, stderr: bytes) -> None:
        if failure_prefix is None:
            return
        write_bytes_exclusive(failure_prefix.with_suffix(".stdout"), stdout)
        write_bytes_exclusive(failure_prefix.with_suffix(".stderr"), stderr)

    try:
        completed = subprocess.run(
            command, env=CHILD_ENVIRONMENT, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, timeout=60.0, check=False)
    except subprocess.TimeoutExpired as error:
        persist_failure(error.stdout or b"", error.stderr or b"")
        raise EvidenceError(
            f"{implementation} timed out after 60 seconds") from error
    elapsed = time.monotonic_ns() - start
    try:
        require(completed.returncode == 0,
                f"{implementation} failed: " +
                completed.stderr.decode(
                    "utf-8", errors="replace")[-1000:])
        require(sha256(executable) == identity["sha256"],
                f"{implementation} binary changed after execution")
        result = json.loads(completed.stdout.decode("utf-8"))
        normalized = validate_result(
            implementation, result, cell, source_commit, source_tree,
            iterations, warmup, target_bytes, campaign, final_selector)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        persist_failure(completed.stdout, completed.stderr)
        raise EvidenceError(
            f"{implementation} output is not one JSON object: {error}") from error
    except Exception:
        persist_failure(completed.stdout, completed.stderr)
        raise
    return {
        "implementation": implementation,
        "elapsed_ns": elapsed,
        "stdout_sha256": hashlib.sha256(completed.stdout).hexdigest(),
        "stderr_sha256": hashlib.sha256(completed.stderr).hexdigest(),
        "normalized": normalized,
        "result": result,
    }


def compact_invocation(
    invocation: Mapping[str, Any],
    artifact: Path,
) -> dict[str, Any]:
    """
    Retain only analysis/provenance fields in the in-memory campaign record.

    The complete benchmark JSON remains in the exclusive per-invocation
    artifact.  This avoids growing resident memory with thousands of retained
    sample arrays while preserving a hash-bound audit trail.
    """
    status = artifact.stat()
    return {
        "implementation": invocation["implementation"],
        "elapsed_ns": invocation["elapsed_ns"],
        "stdout_sha256": invocation["stdout_sha256"],
        "stderr_sha256": invocation["stderr_sha256"],
        "normalized": invocation["normalized"],
        "result_artifact": {
            "name": artifact.name,
            "size": status.st_size,
            "sha256": sha256(artifact),
        },
    }


def confidence_interval(round_log_ratios: Sequence[float]) -> dict[str, Any]:
    round_count = len(round_log_ratios)
    require(round_count in (3, 9),
            "three or nine independent round contrasts are required")
    critical_value = T95_DF2 if round_count == 3 else T95_DF8
    center = statistics.mean(round_log_ratios)
    half = critical_value * statistics.stdev(
        round_log_ratios) / math.sqrt(round_count)
    return {
        "speedup": math.exp(center),
        "ci95": [math.exp(center - half), math.exp(center + half)],
        "round_log_ratios": list(round_log_ratios),
    }


def analyze(
    cell: Mapping[str, Any],
    rounds: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    reference = rounds[0]["invocations"][0]["normalized"]["digests"]
    labels = ("control", "main") if cell["role"] == "target" else ("control",)
    contrasts: dict[str, list[float]] = {label: [] for label in labels}
    for round_value in rounds:
        require(round_value["isolation"]["accepted"] is True,
                "contaminated round cannot be analyzed")
        invocations = round_value["invocations"]
        require(all(item["normalized"]["digests"] == reference
                    for item in invocations),
                "candidate/control/main workload digests differ")
        candidate = [
            item["normalized"]["encode_us"] for item in invocations
            if item["implementation"] == "candidate"
        ]
        require(len(candidate) == 2, "round lacks two candidate observations")
        candidate_log = statistics.mean(math.log(value) for value in candidate)
        for label in labels:
            baseline = [
                item["normalized"]["encode_us"] for item in invocations
                if item["implementation"] == label
            ]
            require(len(baseline) == 2,
                    f"round lacks two {label} observations")
            contrasts[label].append(
                statistics.mean(math.log(value) for value in baseline) -
                candidate_log)
    output = {
        "cell": dict(cell),
        "digests": reference,
        "control_over_candidate": confidence_interval(contrasts["control"]),
    }
    if "main" in contrasts:
        output["main_over_candidate"] = confidence_interval(contrasts["main"])
    return output


def acquire_global_lock() -> int:
    descriptor = os.open(LOCK_PATH, os.O_RDWR | os.O_CREAT, 0o600)
    try:
        fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        os.close(descriptor)
        raise EvidenceError(f"authoritative lock is busy: {LOCK_PATH}")
    return descriptor


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--control", required=True, type=Path)
    parser.add_argument("--main", required=True, type=Path)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--source-tree", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cpu", required=True, type=int)
    parser.add_argument("--sibling", required=True, type=int)
    parser.add_argument("--iterations", type=int, default=9)
    parser.add_argument("--warmup", type=int, default=2)
    parser.add_argument(
        "--target-rounds", type=int, choices=(3, 9), default=3,
        help=(
            "independent balanced target contrasts; neighbors retain three "
            "rounds"
        ))
    parser.add_argument(
        "--target-bytes", type=int,
        choices=(
            64, 128, 192, 256, 320, 384, 448, 512,
            576, 640, 704, 768, 832, 896, 960, 1024,
        ), default=64,
        help="dense target shard size; default preserves the original campaign")
    parser.add_argument(
        "--final-selector", action="store_true",
        help=(
            "validate the narrowed production selector: selected shapes are "
            "targets and excluded shapes are same-source regression neighbors"
        ))
    parser.add_argument(
        "--one-kib-extension", action="store_true",
        help=(
            "qualify only the four new 1024-byte T=8 selector shapes while "
            "treating every other legal K=5..16/R=5..8 shape as a neighbor"
        ))
    parser.add_argument(
        "--tiny-extension", action="store_true",
        help=(
            "qualify the 1..63-byte T=8 selector across every legal "
            "K=5..16/R=5..8 shape; exact main performs an explicitly "
            "labeled zero-padded 64-byte call with matching logical digests"
        ))
    parser.add_argument(
        "--ragged-extension", action="store_true",
        help=(
            "qualify the narrowed 65..928-byte ragged T=8 selector across "
            "every legal K=5..16/R=5..8 shape; exact main processes the "
            "smallest legal zero-padded multiple of 64"
        ))
    parser.add_argument(
        "--cell-id", action="append", default=[],
        help=(
            "run only the named generated cell; repeat for a predeclared "
            "holdout subset (tiny/ragged extension campaigns only)"
        ))
    return parser.parse_args()


def main() -> int:
    options = parse_arguments()
    require(options.iterations >= 3 and options.warmup >= 1,
            "insufficient benchmark repetitions")
    require(
        sum((
            bool(options.one_kib_extension),
            bool(options.tiny_extension),
            bool(options.ragged_extension),
            bool(options.final_selector),
        )) <= 1,
        "--one-kib-extension, --tiny-extension, --ragged-extension, and "
        "--final-selector are distinct campaigns")
    require(not options.cell_id or
            options.tiny_extension or options.ragged_extension,
            "--cell-id is supported only with a padded extension campaign")
    require(len(options.cell_id) == len(set(options.cell_id)),
            "--cell-id values must be unique")
    require(not options.output.exists(), "output path already exists")
    options.output.mkdir(parents=True)
    if options.tiny_extension:
        target_bytes = 0
        campaign = "tiny-extension"
        schema = TINY_SCHEMA
        summary_schema = TINY_SUMMARY_SCHEMA
        cells = tiny_extension_cells()
    elif options.ragged_extension:
        target_bytes = 0
        campaign = "ragged-extension"
        schema = RAGGED_SCHEMA
        summary_schema = RAGGED_SUMMARY_SCHEMA
        cells = ragged_extension_cells()
    elif options.one_kib_extension:
        target_bytes = 1024
        campaign = "one-kib-extension"
        schema = ONE_KIB_SCHEMA
        summary_schema = ONE_KIB_SUMMARY_SCHEMA
        cells = one_kib_extension_cells()
    else:
        target_bytes = options.target_bytes
        campaign = "two-block"
        schema = SCHEMA
        summary_schema = SUMMARY_SCHEMA
        cells = target_cells(
            target_bytes, options.final_selector) + \
            neighbor_cells(target_bytes, options.final_selector)
    if options.cell_id:
        cells_by_id = {str(cell["id"]): cell for cell in cells}
        unknown = sorted(set(options.cell_id) - set(cells_by_id))
        require(not unknown, f"unknown --cell-id values: {unknown}")
        cells = [cells_by_id[cell_id] for cell_id in options.cell_id]
    raw: dict[str, Any] = {
        "schema": schema,
        "created_utc": SUPPORT.utc_now(),
        "source_commit": options.source_commit,
        "source_tree": options.source_tree,
        "main_commit": MAIN_COMMIT,
        "cpu": options.cpu,
        "reserved_sibling": options.sibling,
        "iterations": options.iterations,
        "warmup": options.warmup,
        "target_rounds": options.target_rounds,
        "neighbor_rounds": len(NEIGHBOR_ORDER),
        "target_bytes": target_bytes,
        "target_byte_counts": (
            list(TINY_BYTE_COUNTS) if campaign == "tiny-extension" else
            list(RAGGED_TARGET_BYTE_COUNTS)
            if campaign == "ragged-extension" else [target_bytes]
        ),
        "main_physical_shard_bytes": (
            64 if campaign == "tiny-extension" else
            "ceil(logical_shard_bytes/64)*64"
            if campaign == "ragged-extension" else target_bytes
        ),
        "main_comparison_semantics": (
            "Leopard main processes a zero-padded legal physical shard; "
            "input, parity, and recovery digests cover the matching logical "
            "prefix only"
            if padded_campaign(campaign)
            else "equal physical and logical shard bytes"
        ),
        "campaign": campaign,
        "cell_filter": list(options.cell_id),
        "final_selector": options.final_selector,
        "batch": 64,
        "reuse": 64,
        "runner": regular_file_identity(Path(__file__)),
        "invocation_storage": (
            "hash-bound per-invocation artifacts with compact in-memory index"
            if padded_campaign(campaign)
            else "embedded full invocation records"
        ),
        "isolation_retry_policy": {
            "maximum_attempts_per_round": (
                PADDED_MAX_ISOLATION_ATTEMPTS
                if padded_campaign(campaign) else 1
            ),
            "retry_trigger": (
                "objective CPU-pair isolation rejection only; timing values "
                "are never consulted"
            ),
            "rejected_attempts_retained": padded_campaign(campaign),
        },
        "cells": [],
    }
    lock_descriptor: int | None = None
    try:
        lock_descriptor = acquire_global_lock()
        identities = {
            "candidate": file_identity(options.candidate),
            "control": file_identity(options.control),
            "main": file_identity(options.main),
        }
        raw["identities"] = identities
        require(
            len({identity["sha256"] for identity in identities.values()}) == 3,
            "candidate, control, and main binaries are not distinct")
        executable_sections = {
            name: executable_sections_identity(Path(str(identity["path"])))
            for name, identity in identities.items()
            if name in ("candidate", "control")
        }
        raw["executable_sections"] = executable_sections
        require(
            executable_sections["candidate"]["sections"] ==
                executable_sections["control"]["sections"],
            "candidate and control executable instruction sections differ")
        require(set(os.sched_getaffinity(0)) == {options.cpu},
                "runner must be singleton-pinned to the benchmark CPU")
        sibling_text = Path(
            f"/sys/devices/system/cpu/cpu{options.cpu}/topology/"
            "thread_siblings_list").read_text(encoding="ascii")
        require(SUPPORT.parse_cpu_list(sibling_text) ==
                {options.cpu, options.sibling},
                "requested CPUs are not one SMT pair")
        raw["host"] = {
            "runner_affinity": sorted(os.sched_getaffinity(0)),
            "benchmark_cpu": SUPPORT.cpu_policy_identity(options.cpu),
            "reserved_sibling": SUPPORT.cpu_policy_identity(options.sibling),
        }
        with SUPPORT.StableLeaseAnchor(), \
                SUPPORT.PairLease(options.cpu, options.sibling) as pair_lease:
            raw["pair_lease"] = pair_lease
            before_cpu = SUPPORT.cpu_stat_snapshot(options.cpu)
            before_sibling = SUPPORT.cpu_stat_snapshot(options.sibling)
            before_ns = time.monotonic_ns()
            time.sleep(5.0)
            presample = SUPPORT.isolation_record(
                options.cpu, options.sibling, pair_lease,
                before_ns, time.monotonic_ns(),
                before_cpu, SUPPORT.cpu_stat_snapshot(options.cpu),
                before_sibling, SUPPORT.cpu_stat_snapshot(options.sibling))
            raw["presample"] = presample
            require(
                presample["delta"]["benchmark_cpu"]["nonidle_jiffies"] == 0 and
                presample["delta"]["reserved_sibling"]["nonidle_jiffies"] == 0,
                "CPU pair was not quiet during the presample")
            for cell_index, cell in enumerate(cells):
                orders = (
                    TARGET_ORDER * (options.target_rounds // len(TARGET_ORDER))
                    if cell["role"] == "target"
                    else NEIGHBOR_ORDER
                )
                cell_raw: dict[str, Any] = {
                    "cell": dict(cell),
                    "rounds": [],
                }
                if padded_campaign(campaign):
                    cell_raw["rejected_isolation_attempts"] = []
                for round_index, order in enumerate(orders):
                    maximum_attempts = (
                        PADDED_MAX_ISOLATION_ATTEMPTS
                        if padded_campaign(campaign) else 1
                    )
                    accepted = False
                    for attempt in range(maximum_attempts):
                        before_cpu = SUPPORT.cpu_stat_snapshot(options.cpu)
                        before_sibling = SUPPORT.cpu_stat_snapshot(
                            options.sibling)
                        before_ns = time.monotonic_ns()
                        invocations = []
                        for slot_index, label in enumerate(order):
                            invocation = run_one(
                                label, identities[label], cell, options.cpu,
                                options.source_commit, options.source_tree,
                                options.iterations, options.warmup,
                                target_bytes,
                                campaign=campaign,
                                final_selector=options.final_selector,
                                failure_output=options.output)
                            attempt_suffix = (
                                f"-attempt-{attempt}"
                                if padded_campaign(campaign) else ""
                            )
                            artifact = options.output / (
                                f"partial-{cell['id']}-round-{round_index}"
                                f"{attempt_suffix}-slot-{slot_index}.json"
                            )
                            write_exclusive(artifact, invocation)
                            invocations.append(
                                compact_invocation(invocation, artifact)
                                if padded_campaign(campaign)
                                else invocation
                            )
                        isolation = SUPPORT.isolation_record(
                            options.cpu, options.sibling, pair_lease,
                            before_ns, time.monotonic_ns(),
                            before_cpu, SUPPORT.cpu_stat_snapshot(options.cpu),
                            before_sibling,
                            SUPPORT.cpu_stat_snapshot(options.sibling))
                        round_record = {
                            "round": round_index,
                            "order": list(order),
                            "invocations": invocations,
                            "isolation": isolation,
                        }
                        if padded_campaign(campaign):
                            round_record["attempt"] = attempt
                        if isolation["accepted"]:
                            cell_raw["rounds"].append(round_record)
                            accepted = True
                            break
                        if padded_campaign(campaign):
                            cell_raw[
                                "rejected_isolation_attempts"
                            ].append(round_record)
                        else:
                            cell_raw["rounds"].append(round_record)
                    if not accepted:
                        raw["cells"].append(cell_raw)
                        raise EvidenceError(
                            f"contaminated {cell['id']} round {round_index} "
                            f"for all {maximum_attempts} attempts")
                raw["cells"].append(cell_raw)
                print(
                    f"{cell_index + 1}/{len(cells)} {cell['id']}",
                    file=sys.stderr, flush=True)
        analyses = [
            analyze(item["cell"], item["rounds"]) for item in raw["cells"]
        ]
        target_failures = []
        neighbor_failures = []
        for result in analyses:
            if result["cell"]["role"] == "target":
                if (result["control_over_candidate"]["ci95"][0] <
                        TARGET_CONTROL_FLOOR or
                        result["main_over_candidate"]["ci95"][0] <
                        TARGET_MAIN_FLOOR):
                    target_failures.append(result["cell"]["id"])
            elif result["control_over_candidate"]["ci95"][1] < NEIGHBOR_FLOOR:
                neighbor_failures.append(result["cell"]["id"])
        raw["completed_utc"] = SUPPORT.utc_now()
        write_exclusive(options.output / "raw.json", raw)
        accepted_process_count = sum(
            len(round_value["invocations"])
            for item in raw["cells"] for round_value in item["rounds"])
        rejected_isolation_attempt_count = sum(
            len(item.get("rejected_isolation_attempts", []))
            for item in raw["cells"])
        rejected_process_count = sum(
            len(round_value["invocations"])
            for item in raw["cells"]
            for round_value in item.get("rejected_isolation_attempts", []))
        summary = {
            "schema": summary_schema,
            "status": "accepted" if not target_failures and
                not neighbor_failures else "rejected",
            "source_commit": options.source_commit,
            "source_tree": options.source_tree,
            "main_commit": MAIN_COMMIT,
            "target_bytes": target_bytes,
            "target_byte_counts": (
                list(TINY_BYTE_COUNTS) if campaign == "tiny-extension" else
                list(RAGGED_TARGET_BYTE_COUNTS)
                if campaign == "ragged-extension" else [target_bytes]
            ),
            "main_physical_shard_bytes": (
                64 if campaign == "tiny-extension" else
                "ceil(logical_shard_bytes/64)*64"
                if campaign == "ragged-extension" else target_bytes
            ),
            "main_comparison_semantics": raw["main_comparison_semantics"],
            "runner": raw["runner"],
            "invocation_storage": raw["invocation_storage"],
            "isolation_retry_policy": raw["isolation_retry_policy"],
            "campaign": campaign,
            "cell_filter": raw["cell_filter"],
            "target_rounds": options.target_rounds,
            "neighbor_rounds": len(NEIGHBOR_ORDER),
            "cell_count": len(analyses),
            "target_count": sum(
                result["cell"]["role"] == "target"
                for result in analyses),
            "neighbor_count": sum(
                result["cell"]["role"] != "target"
                for result in analyses),
            "final_selector": options.final_selector,
            "process_count": (
                accepted_process_count + rejected_process_count
            ),
            "accepted_process_count": accepted_process_count,
            "rejected_isolation_process_count": rejected_process_count,
            "rejected_isolation_attempt_count":
                rejected_isolation_attempt_count,
            "all_digests_matched": True,
            "all_accepted_rounds_isolated": all(
                round_value["isolation"]["accepted"] is True
                for item in raw["cells"] for round_value in item["rounds"]
            ),
            "all_rounds_zero_sibling_nonidle": all(
                round_value["isolation"]["delta"][
                    "reserved_sibling"]["nonidle_jiffies"] == 0
                for item in raw["cells"] for round_value in item["rounds"]
            ),
            "target_failures": target_failures,
            "neighbor_failures": neighbor_failures,
            "cells": analyses,
            "binary_sha256": {
                name: identity["sha256"]
                for name, identity in identities.items()
            },
            "candidate_control_executable_sections_sha256":
                executable_sections["candidate"]["combined_sha256"],
            "raw_sha256": sha256(options.output / "raw.json"),
        }
        write_exclusive(options.output / "summary.json", summary)
        print(json.dumps({
            "status": summary["status"],
            "cells": summary["cell_count"],
            "processes": summary["process_count"],
            "target_failures": target_failures,
            "neighbor_failures": neighbor_failures,
        }, sort_keys=True))
        return 0 if summary["status"] == "accepted" else 2
    except Exception as error:
        raw["failed_utc"] = SUPPORT.utc_now()
        raw["failure"] = {
            "type": type(error).__name__,
            "message": str(error),
        }
        write_exclusive(options.output / "failure.json", raw)
        print(f"evidence rejected: {error}", file=sys.stderr)
        return 1
    finally:
        if lock_descriptor is not None:
            os.close(lock_descriptor)


if __name__ == "__main__":
    raise SystemExit(main())
