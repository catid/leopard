#!/usr/bin/env python3
"""Prepare and run the bounded final GF8 Leopard1/Leopard2 comparison.

The broad diagnostic stage is intentionally non-authoritative: it holds the
canonical lock exclusively and uses independent physical cores in parallel.
The confirmation stage holds the same lock exclusively and runs serial ABBA
only for near/slower cells and explicitly selected loss/reuse/batch/tiling
cells.  Both stages are resumable at one JSON document per cell.
"""

import argparse
import ast
import concurrent.futures
import contextlib
import ctypes
import fcntl
import hashlib
import json
import math
import os
import platform
import re
import stat
import statistics
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

from leopard2_build_provenance import (
    candidate_build_provenance,
    verify_reproducible_candidate_build,
)


SCHEMA = "leopard2-expanded-final-map/v2"
MATRIX_SCHEMA = "leopard2-expanded-final-matrix/v2"
MAIN_COMMIT = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
BASELINE_SHA256 = (
    "a43d7f43ff2e887ebcd47a1e94f806847a5d8b858a4e383e6c8d5e528a7dd910")
PADDED_MAIN_BASELINE_SHA256 = (
    "34a73d4481c989cf8e5b881b24c1a066b6bf356a707e01998e9e19d008f2d973")
BASELINE_CONTRACT_LEGACY_V1 = "legacy_exact_v1"
BASELINE_CONTRACT_LOGICAL_ZERO_PAD_V2 = "logical_zero_pad_v2"
CANONICAL_MATRIX_CELL_COUNT = 341
CANONICAL_MATRIX_SHA256 = (
    "3d5d423fe35760727ac1e9751036d9e208be1e69886e27a18c95009204d85d6e")
CANONICAL_MATRIX_PROFILE = "canonical"
R1_TINY_EXACT_MAIN_PROFILE = "r1-tiny-exact-main"
R1_TINY_EXACT_MAIN_K = (3, 4, 8, 24, 32)
R1_TINY_EXACT_MAIN_BYTES = (64, 256)
R1_TINY_EXACT_MAIN_REUSE = 87381
R1_TINY_MISSING_POSITIONS = ("first", "middle", "last")
R1_TINY_EXACT_MAIN_MATRIX_SHA256 = (
    "7a51f3bc4f4a9bffafef1f29b962b1176b04a128a4ef64f8a729e7a077b42695")
K4_R2_PRODUCTION_AVX2_PROFILE = "k4-r2-production-avx2-exact-main"
K4_R2_PRODUCTION_AVX2_BYTES = (63, 64, 65, 256, 1024, 2048, 2049)
K4_R2_PRODUCTION_AVX2_REUSE = 8192
K4_R2_PRODUCTION_AVX2_ITERATIONS = 21
K4_R2_PRODUCTION_AVX2_WARMUP = 32
K4_R2_PRODUCTION_AVX2_MATRIX_SHA256 = (
    "f085476ffe395f7104a7111cfa4d6ec6cfa5d5bf9646f544f6477f47fd187bfb")
K4_R2_TARGET_LOSS_PATTERNS = (
    (63, (0, 1), "adjacent_front"),
    (64, (0, 3), "separated_edges"),
    (65, (2, 3), "adjacent_tail"),
    (256, (1, 2), "adjacent_middle"),
    (1024, (0, 3), "separated_edges"),
    (2048, (0, 1), "adjacent_front"),
    (2049, (1, 3), "separated_interior_tail"),
)
# Keep R=2 neighbors at their maximum two-original loss.  R=1 necessarily
# uses its sole correctable loss, while R=3 uses maximum loss to exercise the
# neighboring T=4 decoder rather than timing a disguised R=2 pattern.
K4_R2_NEIGHBOR_LOSS_PATTERNS = (
    (2, 2, (0, 1), "k_lower_all_originals"),
    (3, 2, (0, 2), "k_lower_separated"),
    (5, 2, (1, 4), "k_upper_separated"),
    (4, 1, (2,), "r_lower_single"),
    (4, 3, (0, 1, 3), "r_upper_maximum"),
)
COMMIT_PATTERN = re.compile(r"^[0-9a-f]{40}$")
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
ABBA_ORDER = ("baseline", "candidate", "candidate", "baseline")
RATIO_METRICS = ("encode", "decode_first", "decode_reuse")
LINUX_F_ADD_SEALS = getattr(fcntl, "F_ADD_SEALS", 1033)
LINUX_F_GET_SEALS = getattr(fcntl, "F_GET_SEALS", 1034)
LINUX_F_SEAL_SEAL = getattr(fcntl, "F_SEAL_SEAL", 0x0001)
LINUX_F_SEAL_SHRINK = getattr(fcntl, "F_SEAL_SHRINK", 0x0002)
LINUX_F_SEAL_GROW = getattr(fcntl, "F_SEAL_GROW", 0x0004)
LINUX_F_SEAL_WRITE = getattr(fcntl, "F_SEAL_WRITE", 0x0008)
LINUX_MFD_CLOEXEC = getattr(os, "MFD_CLOEXEC", 0x0001)
LINUX_MFD_ALLOW_SEALING = getattr(os, "MFD_ALLOW_SEALING", 0x0002)

# This is the frozen 144-cell final-map shape set (48 shapes x three sizes).
CORE_SHAPES = (
    ("t1-k1-r1", 1, 1),
    ("t1-k8-r1", 8, 1),
    ("t1-k129-r1", 129, 1),
    ("t1-k255-r1", 255, 1),
    ("t2-k2-r2", 2, 2),
    ("t2-k3-r2", 3, 2),
    ("t2-k8-r2", 8, 2),
    ("t2-k64-r2", 64, 2),
    ("t2-k254-r2", 254, 2),
    ("t4-k3-r3", 3, 3),
    ("t4-k4-r3", 4, 3),
    ("t4-k7-r3", 7, 3),
    ("t4-k8-r3", 8, 3),
    ("t4-k11-r3", 11, 3),
    ("t4-k12-r3", 12, 3),
    ("t4-k15-r3", 15, 3),
    ("t4-k16-r3", 16, 3),
    ("t4-k64-r3", 64, 3),
    ("t4-k252-r3", 252, 3),
    ("t4-k4-r4", 4, 4),
    ("t4-k7-r4", 7, 4),
    ("t4-k8-r4", 8, 4),
    ("t4-k11-r4", 11, 4),
    ("t4-k12-r4", 12, 4),
    ("t4-k15-r4", 15, 4),
    ("t4-k16-r4", 16, 4),
    ("t4-k64-r4", 64, 4),
    ("t4-k252-r4", 252, 4),
    ("t8-k8-r8", 8, 8),
    ("t8-k9-r8", 9, 8),
    ("t8-k16-r8", 16, 8),
    ("t8-k64-r8", 64, 8),
    ("t8-k248-r8", 248, 8),
    ("t16-k16-r16", 16, 16),
    ("t16-k17-r16", 17, 16),
    ("t16-k32-r16", 32, 16),
    ("t16-k64-r16", 64, 16),
    ("t16-k240-r16", 240, 16),
    ("t32-k32-r32", 32, 32),
    ("t32-k33-r32", 33, 32),
    ("t32-k64-r32", 64, 32),
    ("t32-k128-r32", 128, 32),
    ("t32-k224-r32", 224, 32),
    ("t64-k64-r64", 64, 64),
    ("t64-k65-r64", 65, 64),
    ("t64-k128-r64", 128, 64),
    ("t64-k192-r64", 192, 64),
    ("t128-k128-r128", 128, 128),
)

# Adds transform-boundary and non-power-R cells absent from the old map.
BOUNDARY_SHAPES = (
    ("t4-k5-r4", 5, 4),
    ("t8-k7-r7", 7, 7),
    ("t16-k15-r15", 15, 15),
    ("t32-k31-r31", 31, 31),
    ("t32-k34-r32", 34, 32),
    ("t64-k63-r63", 63, 63),
    ("t128-k127-r127", 127, 127),
)

TINY_SHAPES = (
    (1, 1), (129, 1), (3, 2), (64, 2), (4, 4), (16, 4),
    (7, 7), (9, 8), (15, 15), (17, 16), (31, 31), (34, 32),
    (65, 64), (128, 128),
)
LOSS_SHAPES = (
    (129, 1), (3, 2), (16, 4), (9, 8), (17, 16),
    (34, 32), (128, 32), (65, 64), (192, 64), (128, 128),
)
REUSE_SHAPES = LOSS_SHAPES
BATCH_SHAPES = LOSS_SHAPES
BATCH_TINY_SHAPES = (
    (3, 2), (16, 4), (9, 8), (17, 16), (34, 32), (65, 64),
)
LARGE_SHAPES = (
    (129, 1), (16, 4), (9, 8), (17, 16), (34, 32),
    (128, 32), (65, 64), (192, 64), (128, 128),
)
SELECTOR_ISOLATION_SHAPES = ((33, 32), (34, 32), (35, 32))

# The large-shard R=1 selector deliberately changes implementation at these
# four measured boundaries.  Retain the immediate lower-K fallback alongside
# each promoted cell, and exercise recovery rather than only no-loss decode.
R1_SELECTOR_BOUNDARIES = (
    (50, 51, 1048576),
    (30, 31, 2097152),
    (9, 10, 4194304),
    (7, 8, 8388608),
)

# The AVX2 source-major direct executor has specialized even-output kernels
# and composed odd-output paths.  The ordinary loss sweep contains only
# powers of two, so keep each non-power arity at the dispatch boundary, a
# cache-resident representative, and a large-shard tiling point.
DIRECT_ODD_LOSSES = (3, 5, 6, 7)

BYTE_LABELS = {
    63: "63b", 64: "64b", 65: "65b", 256: "256b", 1024: "1k",
    2048: "2k", 2049: "2049b",
    4096: "4k", 8192: "8k", 16384: "16k", 65536: "64k",
    1048576: "1m", 2097152: "2m", 4194304: "4m",
    8388608: "8m", 16777216: "16m",
}


def ceil_pow2(value):
    return 1 << (value - 1).bit_length()


def canonical_bytes(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def sha256_bytes(value):
    return hashlib.sha256(value).hexdigest()


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp.%d" % os.getpid())
    temporary.write_bytes(canonical_bytes(value) + b"\n")
    os.replace(str(temporary), str(path))


def read_json(path):
    with open(path, "r", encoding="utf-8") as stream:
        return json.load(stream)


def read_json_snapshot(path, label, maximum_bytes=64 << 20):
    """Parse exactly the bytes whose identity is retained by the caller."""
    resolved = Path(path).resolve(strict=True)
    descriptor = os.open(
        resolved, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0))
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise RuntimeError(label + " is not a regular file")
        if before.st_size < 1 or before.st_size > maximum_bytes:
            raise RuntimeError(label + " exceeds its retained byte bound")
        chunks = []
        remaining = before.st_size
        while remaining:
            block = os.read(descriptor, min(1 << 20, remaining))
            if not block:
                raise RuntimeError(label + " ended before its recorded size")
            chunks.append(block)
            remaining -= len(block)
        if os.read(descriptor, 1):
            raise RuntimeError(label + " grew while it was read")
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    stable_fields = (
        "st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns",
        "st_ctime_ns")
    if any(getattr(before, field) != getattr(after, field)
           for field in stable_fields):
        raise RuntimeError(label + " changed while it was read")
    path_status = resolved.stat()
    if any(getattr(after, field) != getattr(path_status, field)
           for field in stable_fields):
        raise RuntimeError(label + " pathname changed while it was read")
    content = b"".join(chunks)
    try:
        document = json.loads(content.decode("utf-8", errors="strict"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise RuntimeError(label + " is not valid strict UTF-8 JSON: " +
                           str(error)) from error
    identity = {
        "path": str(resolved), "sha256": sha256_bytes(content),
        "device": after.st_dev, "inode": after.st_ino,
        "mode": after.st_mode, "size": after.st_size,
        "mtime_ns": after.st_mtime_ns, "ctime_ns": after.st_ctime_ns,
    }
    return document, identity


def require_json_snapshot_unchanged(identity, label):
    unused_document, observed = read_json_snapshot(identity["path"], label)
    if observed != identity:
        raise RuntimeError(label + " changed after its exact bytes were used")


def stable_seed(key):
    digest = hashlib.sha256((MATRIX_SCHEMA + ":" + key).encode()).digest()
    return int.from_bytes(digest[:8], "big") & ((1 << 63) - 1)


def xorshift64_next(state):
    mask = (1 << 64) - 1
    value = state
    value = (value ^ (value << 13)) & mask
    value = (value ^ (value >> 7)) & mask
    value = (value ^ (value << 17)) & mask
    return value


def selected_missing_indices(k, loss, seed):
    if k <= 0 or loss < 0 or loss > k:
        raise RuntimeError("loss selection requires 0 <= loss <= K")
    state = (seed ^ 0xd1b54a32d192ed03) & ((1 << 64) - 1)
    if state == 0:
        state = 0x9e3779b97f4a7c15
    order = list(range(k))
    for remaining in range(k, 1, -1):
        state = xorshift64_next(state)
        selected = state % remaining
        order[remaining - 1], order[selected] = (
            order[selected], order[remaining - 1])
    return tuple(sorted(order[:loss]))


def r1_missing_index(k, seed):
    return selected_missing_indices(k, 1, seed)[0]


def r1_missing_position_index(k, position):
    if position == "first":
        return 0
    if position == "middle":
        return k // 2
    if position == "last":
        return k - 1
    raise RuntimeError("unknown R=1 missing-position label: " + position)


def r1_seed_for_missing_index(identifier, k, missing_index, used_seeds):
    return seed_for_missing_indices(
        identifier, k, (missing_index,), used_seeds)


def seed_for_missing_indices(identifier, k, missing_indices, used_seeds):
    missing_indices = tuple(sorted(missing_indices))
    if (not missing_indices or len(missing_indices) != len(set(missing_indices))
            or missing_indices[0] < 0 or missing_indices[-1] >= k):
        raise RuntimeError("expected missing indices are invalid")
    modulus = 1 << 63
    initial = stable_seed(identifier)
    for offset in range(1 << 20):
        seed = (initial + offset) % modulus
        if seed in used_seeds:
            continue
        if selected_missing_indices(
                k, len(missing_indices), seed) == missing_indices:
            used_seeds.add(seed)
            return seed
    raise RuntimeError(
        "could not find a bounded unique loss seed for " + identifier)


def iteration_policy(byte_count, confirmation):
    if byte_count <= 1024:
        return (21 if confirmation else 15, 2)
    if byte_count <= 4096:
        return (11 if confirmation else 7, 2)
    if byte_count <= 65536:
        return (9 if confirmation else 5, 2)
    if byte_count <= 1048576:
        return (3 if confirmation else 2, 1)
    return (1, 1)


def cell_key(k, r, byte_count, loss, reuse, batch):
    return "k%d-r%d-%s-l%d-u%d-b%d" % (
        k, r, BYTE_LABELS[byte_count], loss, reuse, batch)


def add_cell(cells, k, r, byte_count, loss, reuse, batch, tag):
    if (k <= 0 or r <= 0 or
            ceil_pow2(k + ceil_pow2(r)) > 256):
        raise RuntimeError("matrix contains an invalid GF8 public code")
    if loss < 0 or loss > min(k, r):
        raise RuntimeError("matrix contains an invalid loss count")
    key = (k, r, byte_count, loss, reuse, batch)
    if key not in cells:
        diagnostic_iterations, diagnostic_warmup = iteration_policy(
            byte_count, False)
        abba_iterations, abba_warmup = iteration_policy(byte_count, True)
        identifier = cell_key(k, r, byte_count, loss, reuse, batch)
        cells[key] = {
            "id": identifier,
            "K": k,
            "R": r,
            "T": ceil_pow2(r),
            "bytes": byte_count,
            "loss": loss,
            "reuse": reuse,
            "batch": batch,
            "seed": stable_seed(identifier),
            "diagnostic_iterations": diagnostic_iterations,
            "diagnostic_warmup": diagnostic_warmup,
            "abba_iterations": abba_iterations,
            "abba_warmup": abba_warmup,
            "estimated_peak_bytes": (
                6 * (k + r) * byte_count * batch + (64 << 20)),
            "tags": [],
        }
    if tag not in cells[key]["tags"]:
        cells[key]["tags"].append(tag)


def loss_values(k, r):
    values = {0, 1, 2, 4, r // 2, r}
    return sorted(value for value in values if value <= min(k, r))


def make_canonical_matrix():
    cells = {}
    all_shapes = CORE_SHAPES + BOUNDARY_SHAPES
    for unused_name, k, r in all_shapes:
        for byte_count in (2048, 4096, 65536):
            add_cell(cells, k, r, byte_count, min(k, r), 8, 1, "core")
    for k, r in TINY_SHAPES:
        for byte_count in (64, 256, 1024):
            add_cell(cells, k, r, byte_count, min(k, r), 8, 1,
                     "tiny_shard")
        add_cell(cells, k, r, 8192, min(k, r), 8, 1, "eight_kib")
    for k, r in LOSS_SHAPES:
        for loss in loss_values(k, r):
            add_cell(cells, k, r, 65536, loss, 8, 1, "loss_sweep")
    for k, r in REUSE_SHAPES:
        loss = r // 2
        for reuse in (1, 8, 64):
            add_cell(cells, k, r, 65536, loss, reuse, 1, "reuse_sweep")
    for k, r in BATCH_SHAPES:
        loss = r // 2
        for batch in (1, 8):
            add_cell(cells, k, r, 65536, loss, 8, batch, "batch_sweep")
    for k, r in BATCH_TINY_SHAPES:
        add_cell(cells, k, r, 256, r // 2, 8, 8, "batch_sweep")
    for k, r in LARGE_SHAPES:
        for byte_count in (1048576, 16777216):
            add_cell(cells, k, r, byte_count, r // 2, 8, 1,
                     "large_tiling")
    for k, r in SELECTOR_ISOLATION_SHAPES:
        for byte_count in (16384, 65536):
            add_cell(cells, k, r, byte_count, r, 8, 1,
                     "selector_isolation")
    for lower_k, promoted_k, byte_count in R1_SELECTOR_BOUNDARIES:
        add_cell(cells, lower_k, 1, byte_count, 1, 8, 1,
                 "r1_selector_fallback")
        add_cell(cells, promoted_k, 1, byte_count, 1, 8, 1,
                 "r1_selector_promoted")
    for k in (129, 255):
        add_cell(cells, k, 1, 1048576, 1, 8, 1,
                 "r1_large_recovery")
    for byte_count in (2048, 65536, 1048576):
        for loss in DIRECT_ODD_LOSSES:
            # R=65 rounds to the same P=T=128 construction as K=65 and is
            # therefore the measured expanded direct-repair profile.  R=64
            # takes the transform decoder and would not exercise this path.
            add_cell(cells, 65, 65, byte_count, loss, 8, 1,
                     "direct_odd_arity")
    result = sorted(cells.values(), key=lambda cell: cell["id"])
    for index, cell in enumerate(result):
        cell["index"] = index
        cell["tags"].sort()
    matrix = {
        "schema": MATRIX_SCHEMA,
        "description": (
            "Bounded broad diagnostic plus serial near/loss ABBA confirmation"),
        "baseline_commit": MAIN_COMMIT,
        "baseline_sha256": BASELINE_SHA256,
        "cell_count": len(result),
        "diagnostic_child_invocations": 2 * len(result),
        "maximum_abba_child_invocations_at_three_rounds": 12 * len(result),
        "sizes": sorted({cell["bytes"] for cell in result}),
        "transform_sides": sorted({cell["T"] for cell in result}),
        "losses": sorted({cell["loss"] for cell in result}),
        "reuse_counts": sorted({cell["reuse"] for cell in result}),
        "batch_counts": sorted({cell["batch"] for cell in result}),
        "cells": result,
    }
    return matrix


def make_r1_tiny_exact_main_matrix():
    cells = []
    used_seeds = set()
    for k in R1_TINY_EXACT_MAIN_K:
        for byte_count in R1_TINY_EXACT_MAIN_BYTES:
            for position in R1_TINY_MISSING_POSITIONS:
                missing_index = r1_missing_position_index(k, position)
                identifier = "%s-missing-%s" % (
                    cell_key(
                        k, 1, byte_count, 1,
                        R1_TINY_EXACT_MAIN_REUSE, 1),
                    position,
                )
                seed = r1_seed_for_missing_index(
                    identifier, k, missing_index, used_seeds)
                diagnostic_iterations, diagnostic_warmup = iteration_policy(
                    byte_count, False)
                abba_iterations, abba_warmup = iteration_policy(
                    byte_count, True)
                cells.append({
                    "id": identifier,
                    "K": k,
                    "R": 1,
                    "T": 1,
                    "bytes": byte_count,
                    "loss": 1,
                    "reuse": R1_TINY_EXACT_MAIN_REUSE,
                    "batch": 1,
                    "seed": seed,
                    "expected_missing_original_index": missing_index,
                    "diagnostic_iterations": diagnostic_iterations,
                    "diagnostic_warmup": diagnostic_warmup,
                    "abba_iterations": abba_iterations,
                    "abba_warmup": abba_warmup,
                    "estimated_peak_bytes": (
                        6 * (k + 1) * byte_count + (64 << 20)),
                    "tags": [
                        "r1_tiny_exact_main", "missing_" + position],
                })
    cells.sort(key=lambda cell: cell["id"])
    for index, cell in enumerate(cells):
        cell["index"] = index
        cell["tags"].sort()
    return {
        "schema": MATRIX_SCHEMA,
        "profile": R1_TINY_EXACT_MAIN_PROFILE,
        "description": (
            "Tiny R=1 exact-main first/middle/last loss confirmation"),
        "baseline_commit": MAIN_COMMIT,
        "baseline_sha256": BASELINE_SHA256,
        "cell_count": len(cells),
        "diagnostic_child_invocations": 2 * len(cells),
        "maximum_abba_child_invocations_at_three_rounds": 12 * len(cells),
        "sizes": list(R1_TINY_EXACT_MAIN_BYTES),
        "transform_sides": [1],
        "losses": [1],
        "reuse_counts": [R1_TINY_EXACT_MAIN_REUSE],
        "batch_counts": [1],
        "cells": cells,
    }


def make_k4_r2_production_avx2_matrix():
    cells = []
    used_seeds = set()

    def append_cell(k, r, byte_count, missing_indices, pattern, role):
        loss = len(missing_indices)
        if loss < 1 or loss > min(k, r):
            raise RuntimeError("K4/R2 gate contains an invalid loss pattern")
        identifier = "%s-missing-%s" % (
            cell_key(
                k, r, byte_count, loss,
                K4_R2_PRODUCTION_AVX2_REUSE, 1),
            "-".join(str(index) for index in missing_indices),
        )
        seed = seed_for_missing_indices(
            identifier, k, missing_indices, used_seeds)
        cells.append({
            "id": identifier,
            "K": k,
            "R": r,
            "T": ceil_pow2(r),
            "bytes": byte_count,
            "loss": loss,
            "loss_pattern": pattern,
            "reuse": K4_R2_PRODUCTION_AVX2_REUSE,
            "batch": 1,
            "seed": seed,
            "expected_missing_original_indices": list(missing_indices),
            "diagnostic_iterations": K4_R2_PRODUCTION_AVX2_ITERATIONS,
            "diagnostic_warmup": K4_R2_PRODUCTION_AVX2_WARMUP,
            "abba_iterations": K4_R2_PRODUCTION_AVX2_ITERATIONS,
            "abba_warmup": K4_R2_PRODUCTION_AVX2_WARMUP,
            "estimated_peak_bytes": (
                6 * (k + r) * byte_count + (64 << 20)),
            "tags": [
                "k4_r2_production_avx2_gate", role,
                "loss_" + pattern,
            ],
        })

    for byte_count, missing_indices, pattern in K4_R2_TARGET_LOSS_PATTERNS:
        append_cell(
            4, 2, byte_count, missing_indices, pattern, "k4_r2_target")
    for k, r, missing_indices, pattern in K4_R2_NEIGHBOR_LOSS_PATTERNS:
        append_cell(
            k, r, 64, missing_indices, pattern, "k4_r2_neighbor")

    cells.sort(key=lambda cell: cell["id"])
    for index, cell in enumerate(cells):
        cell["index"] = index
        cell["tags"].sort()
    return {
        "schema": MATRIX_SCHEMA,
        "profile": K4_R2_PRODUCTION_AVX2_PROFILE,
        "description": (
            "Strict exact-main K=4,R=2 AVX2 production gate with adjacent "
            "and separated maximum-loss target patterns plus K/R neighbors"),
        "loss_pattern_policy": {
            "target": (
                "K=4,R=2 always loses the maximum two originals; adjacent "
                "pairs represent burst losses and separated pairs represent "
                "independent losses"),
            "neighbors": (
                "R=2 neighbors lose two originals, R=1 loses one, and the "
                "R=3 neighbor loses three to exercise its maximum-loss T=4 "
                "decoder"),
        },
        "baseline_commit": MAIN_COMMIT,
        "baseline_sha256": PADDED_MAIN_BASELINE_SHA256,
        "cell_count": len(cells),
        "diagnostic_child_invocations": 2 * len(cells),
        "maximum_abba_child_invocations_at_three_rounds": 12 * len(cells),
        "sizes": sorted({cell["bytes"] for cell in cells}),
        "transform_sides": sorted({cell["T"] for cell in cells}),
        "losses": sorted({cell["loss"] for cell in cells}),
        "reuse_counts": [K4_R2_PRODUCTION_AVX2_REUSE],
        "batch_counts": [1],
        "cells": cells,
    }


def make_matrix(profile=CANONICAL_MATRIX_PROFILE):
    if profile == CANONICAL_MATRIX_PROFILE:
        return make_canonical_matrix()
    if profile == R1_TINY_EXACT_MAIN_PROFILE:
        return make_r1_tiny_exact_main_matrix()
    if profile == K4_R2_PRODUCTION_AVX2_PROFILE:
        return make_k4_r2_production_avx2_matrix()
    raise RuntimeError("unknown matrix profile: " + profile)


def matrix_digest(matrix):
    return sha256_bytes(canonical_bytes(matrix))


def validate_matrix(matrix):
    if matrix.get("schema") != MATRIX_SCHEMA:
        raise RuntimeError("matrix schema mismatch")
    profile = matrix.get("profile", CANONICAL_MATRIX_PROFILE)
    regenerated = make_matrix(profile)
    if matrix != regenerated:
        raise RuntimeError("matrix differs from the deterministic built-in matrix")
    ids = [cell["id"] for cell in matrix["cells"]]
    seeds = [cell["seed"] for cell in matrix["cells"]]
    if len(ids) != len(set(ids)) or len(seeds) != len(set(seeds)):
        raise RuntimeError("matrix IDs or seeds are not unique")
    if profile == R1_TINY_EXACT_MAIN_PROFILE:
        if len(matrix["cells"]) != (
                len(R1_TINY_EXACT_MAIN_K) *
                len(R1_TINY_EXACT_MAIN_BYTES) *
                len(R1_TINY_MISSING_POSITIONS)):
            raise RuntimeError("R=1 tiny exact-main cell count drifted")
        if matrix["baseline_sha256"] != BASELINE_SHA256:
            raise RuntimeError("R=1 tiny exact-main baseline identity drifted")
        if matrix_digest(matrix) != R1_TINY_EXACT_MAIN_MATRIX_SHA256:
            raise RuntimeError("R=1 tiny exact-main matrix SHA-256 drifted")
        for cell in matrix["cells"]:
            if r1_missing_index(cell["K"], cell["seed"]) != \
                    cell["expected_missing_original_index"]:
                raise RuntimeError("R=1 tiny exact-main seed coverage drifted")
        return
    if profile == K4_R2_PRODUCTION_AVX2_PROFILE:
        expected_count = (
            len(K4_R2_TARGET_LOSS_PATTERNS) +
            len(K4_R2_NEIGHBOR_LOSS_PATTERNS))
        if len(matrix["cells"]) != expected_count:
            raise RuntimeError("K4/R2 production AVX2 cell count drifted")
        if matrix["baseline_sha256"] != PADDED_MAIN_BASELINE_SHA256:
            raise RuntimeError(
                "K4/R2 production AVX2 baseline identity drifted")
        if matrix_digest(matrix) != K4_R2_PRODUCTION_AVX2_MATRIX_SHA256:
            raise RuntimeError(
                "K4/R2 production AVX2 matrix SHA-256 drifted")
        target_cells = [
            cell for cell in matrix["cells"]
            if "k4_r2_target" in cell["tags"]]
        neighbor_cells = [
            cell for cell in matrix["cells"]
            if "k4_r2_neighbor" in cell["tags"]]
        if (tuple(sorted(cell["bytes"] for cell in target_cells)) !=
                K4_R2_PRODUCTION_AVX2_BYTES or
                any((cell["K"], cell["R"], cell["loss"]) != (4, 2, 2)
                    for cell in target_cells) or
                any(cell["bytes"] != 64 for cell in neighbor_cells)):
            raise RuntimeError(
                "K4/R2 production AVX2 target/neighbor grid drifted")
        for cell in matrix["cells"]:
            expected_missing = tuple(
                cell["expected_missing_original_indices"])
            if selected_missing_indices(
                    cell["K"], cell["loss"], cell["seed"]) != \
                    expected_missing:
                raise RuntimeError(
                    "K4/R2 production AVX2 seed coverage drifted")
            if (cell["reuse"] != K4_R2_PRODUCTION_AVX2_REUSE or
                    cell["batch"] != 1 or
                    cell["diagnostic_iterations"] !=
                    K4_R2_PRODUCTION_AVX2_ITERATIONS or
                    cell["diagnostic_warmup"] !=
                    K4_R2_PRODUCTION_AVX2_WARMUP or
                    cell["abba_iterations"] !=
                    K4_R2_PRODUCTION_AVX2_ITERATIONS or
                    cell["abba_warmup"] !=
                    K4_R2_PRODUCTION_AVX2_WARMUP):
                raise RuntimeError(
                    "K4/R2 production AVX2 timing policy drifted")
            if "k4_r2_production_avx2_gate" not in cell["tags"]:
                raise RuntimeError(
                    "K4/R2 production AVX2 mandatory tag drifted")
        return
    if profile != CANONICAL_MATRIX_PROFILE:
        raise RuntimeError("unknown matrix profile: " + str(profile))
    if len(matrix["cells"]) != CANONICAL_MATRIX_CELL_COUNT:
        raise RuntimeError("canonical matrix cell count drifted")
    if matrix["baseline_sha256"] != BASELINE_SHA256:
        raise RuntimeError("canonical matrix baseline identity drifted")
    if matrix_digest(matrix) != CANONICAL_MATRIX_SHA256:
        raise RuntimeError("canonical matrix SHA-256 drifted")
    required_sizes = {64, 256, 1024, 2048, 4096, 8192, 16384, 65536,
                      1048576, 2097152, 4194304, 8388608, 16777216}
    if set(matrix["sizes"]) != required_sizes:
        raise RuntimeError("matrix size coverage drifted")
    if set(matrix["transform_sides"]) != {1, 2, 4, 8, 16, 32, 64, 128}:
        raise RuntimeError("matrix transform-side coverage drifted")
    if not {0, 1, 2, 4}.issubset(set(matrix["losses"])):
        raise RuntimeError("matrix small-loss coverage drifted")
    if matrix["reuse_counts"] != [1, 8, 64] or matrix["batch_counts"] != [1, 8]:
        raise RuntimeError("matrix reuse/batch coverage drifted")


def checked_output(command, cwd=None, timeout=60):
    completed = subprocess.run(
        command, cwd=cwd, stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        timeout=timeout, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            "%r failed rc=%d: %s" %
            (command, completed.returncode, completed.stderr.strip()))
    return completed.stdout.strip()


def git_output(source, *arguments):
    git = Path("/usr/bin/git").resolve(strict=True)
    environment = {
        "GIT_CONFIG_GLOBAL": "/dev/null",
        "GIT_CONFIG_NOSYSTEM": "1",
        "GIT_NO_REPLACE_OBJECTS": "1",
        "GIT_OPTIONAL_LOCKS": "0",
        "LANG": "C", "LC_ALL": "C", "PATH": "/usr/bin:/bin",
    }
    completed = subprocess.run(
        [str(git), "-C", str(source), *arguments],
        stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, text=True, env=environment,
        timeout=60, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            "git %s failed rc=%d: %s" % (
                " ".join(arguments), completed.returncode,
                completed.stderr.strip()))
    return completed.stdout.strip()


def linux_memfd_create(name):
    flags = LINUX_MFD_CLOEXEC | LINUX_MFD_ALLOW_SEALING
    if hasattr(os, "memfd_create"):
        return os.memfd_create(name, flags)
    libc = ctypes.CDLL(None, use_errno=True)
    creator = getattr(libc, "memfd_create", None)
    if creator is None:
        raise RuntimeError(
            "benchmark snapshots require Linux memfd_create support")
    creator.argtypes = (ctypes.c_char_p, ctypes.c_uint)
    creator.restype = ctypes.c_int
    descriptor = creator(name.encode("utf-8"), flags)
    if descriptor < 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error), name)
    return descriptor


def file_identity(path, label):
    resolved = Path(path).resolve(strict=True)
    descriptor = os.open(
        resolved, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0))
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise RuntimeError(label + " is not a regular file")
        digest = hashlib.sha256()
        while True:
            block = os.read(descriptor, 1 << 20)
            if not block:
                break
            digest.update(block)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    stable_fields = (
        "st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns",
        "st_ctime_ns")
    if any(getattr(before, name) != getattr(after, name)
           for name in stable_fields):
        raise RuntimeError(label + " changed while it was hashed")
    path_status = resolved.stat()
    if any(getattr(after, name) != getattr(path_status, name)
           for name in stable_fields):
        raise RuntimeError(label + " path changed while it was hashed")
    if not os.access(resolved, os.X_OK):
        raise RuntimeError(label + " is not executable")
    return {
        "path": str(resolved), "sha256": digest.hexdigest(),
        "device": after.st_dev, "inode": after.st_ino,
        "mode": after.st_mode, "size": after.st_size,
        "mtime_ns": after.st_mtime_ns, "ctime_ns": after.st_ctime_ns,
    }


def equivalent_executable_identities(build_path, artifact_path):
    build = file_identity(build_path, "candidate build-tree executable")
    artifact = file_identity(artifact_path, "candidate timing artifact")
    if (build["size"], build["sha256"]) != \
            (artifact["size"], artifact["sha256"]):
        raise RuntimeError(
            "candidate build-tree executable and timing artifact differ")
    build_descriptor = os.open(
        build["path"], os.O_RDONLY | getattr(os, "O_CLOEXEC", 0))
    artifact_descriptor = os.open(
        artifact["path"], os.O_RDONLY | getattr(os, "O_CLOEXEC", 0))
    try:
        offset = 0
        while offset < build["size"]:
            length = min(1 << 20, build["size"] - offset)
            build_block = os.pread(build_descriptor, length, offset)
            artifact_block = os.pread(artifact_descriptor, length, offset)
            if (len(build_block) != length or
                    artifact_block != build_block):
                raise RuntimeError(
                    "candidate build-tree executable and timing artifact "
                    "differ byte-for-byte")
            offset += length
    finally:
        os.close(artifact_descriptor)
        os.close(build_descriptor)
    if file_identity(build["path"], "candidate build-tree executable") != \
            build or file_identity(
                artifact["path"], "candidate timing artifact") != artifact:
        raise RuntimeError(
            "candidate executable changed during byte-equivalence validation")
    return {
        "schema": "leopard2-executable-equivalence/v1",
        "build_tree_executable": build,
        "selected_timing_artifact": artifact,
        "equal_size": build["size"],
        "equal_sha256": build["sha256"],
        "byte_for_byte_equal": True,
    }


def sealed_snapshot_identity(descriptor, label):
    status = os.fstat(descriptor)
    if not stat.S_ISREG(status.st_mode):
        raise RuntimeError(label + " sealed snapshot is not regular")
    digest = hashlib.sha256()
    offset = 0
    while offset < status.st_size:
        block = os.pread(
            descriptor, min(1 << 20, status.st_size - offset), offset)
        if not block:
            raise RuntimeError(label + " sealed snapshot ended early")
        digest.update(block)
        offset += len(block)
    seals = fcntl.fcntl(descriptor, LINUX_F_GET_SEALS)
    required = (LINUX_F_SEAL_SEAL | LINUX_F_SEAL_SHRINK |
                LINUX_F_SEAL_GROW | LINUX_F_SEAL_WRITE)
    if seals & required != required:
        raise RuntimeError(label + " snapshot is not immutably sealed")
    return {
        "kind": "linux-sealed-memfd-v1", "sha256": digest.hexdigest(),
        "size": status.st_size, "mode": status.st_mode, "seals": seals,
    }


def snapshot_executable(path, label, expected_sha256=None):
    if not sys.platform.startswith("linux") or not hasattr(os, "pread"):
        raise RuntimeError(
            "benchmark snapshots require Linux sealed memfd support")
    source_identity = file_identity(path, label)
    if expected_sha256 is not None and \
            source_identity["sha256"] != expected_sha256:
        raise RuntimeError(
            "frozen baseline SHA-256 mismatch: expected %s, got %s" %
            (expected_sha256, source_identity["sha256"]))
    source = os.open(
        source_identity["path"], os.O_RDONLY | getattr(os, "O_CLOEXEC", 0))
    snapshot = -1
    try:
        snapshot = linux_memfd_create(
            "leopard2-expanded-" + label.replace(" ", "-"))
        copied_digest = hashlib.sha256()
        while True:
            block = os.read(source, 1 << 20)
            if not block:
                break
            copied_digest.update(block)
            view = memoryview(block)
            while view:
                written = os.write(snapshot, view)
                if written <= 0:
                    raise RuntimeError(label + " snapshot copy made no progress")
                view = view[written:]
        if copied_digest.hexdigest() != source_identity["sha256"]:
            raise RuntimeError(label + " changed before snapshot completion")
        if os.pread(snapshot, 4, 0) != b"\x7fELF":
            raise RuntimeError(label + " is not an ELF executable")
        if file_identity(source_identity["path"], label) != source_identity:
            raise RuntimeError(label + " changed while snapshotting")
        os.fchmod(snapshot, 0o500)
        required = (LINUX_F_SEAL_SEAL | LINUX_F_SEAL_SHRINK |
                    LINUX_F_SEAL_GROW | LINUX_F_SEAL_WRITE)
        fcntl.fcntl(snapshot, LINUX_F_ADD_SEALS, required)
        snapshot_identity = sealed_snapshot_identity(snapshot, label)
        if (snapshot_identity["sha256"] != source_identity["sha256"] or
                snapshot_identity["size"] != source_identity["size"]):
            raise RuntimeError(label + " sealed snapshot identity mismatch")
        return source_identity, snapshot, snapshot_identity
    except BaseException:
        if snapshot >= 0:
            os.close(snapshot)
        raise
    finally:
        os.close(source)


def snapshot_path(descriptor):
    return "/proc/self/fd/%d" % descriptor


def inspect_isa(path, pass_fds=()):
    completed = subprocess.run(
        ["/usr/bin/objdump", "-d", "-M", "intel", path],
        stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, text=True, timeout=180, check=False,
        pass_fds=pass_fds)
    if completed.returncode != 0:
        raise RuntimeError("objdump failed: " + completed.stderr.strip())
    evex = 0
    ymm = 0
    instruction = re.compile(r"^\s*[0-9a-f]+:\s+((?:[0-9a-f]{2}\s+)+)")
    for line in completed.stdout.splitlines():
        match = instruction.match(line)
        if match and match.group(1).split()[0] == "62":
            evex += 1
        if re.search(r"\bymm[0-9]+\b", line):
            ymm += 1
    if evex:
        raise RuntimeError("%s contains %d EVEX instructions" % (path, evex))
    if not ymm:
        raise RuntimeError("%s contains no visible YMM instructions" % path)
    return {
        "disassembler": checked_output(
            ["/usr/bin/objdump", "--version"]).splitlines()[0],
        "evex_prefixed_instruction_count": evex,
        "ymm_operand_instruction_count": ymm,
    }


def source_provenance(source, expected_commit):
    source = os.path.realpath(source)
    if not os.path.isdir(source):
        raise RuntimeError("candidate source is not a directory: " + source)
    if not COMMIT_PATTERN.match(expected_commit):
        raise RuntimeError("candidate commit must be a full lowercase SHA")
    replace_refs = git_output(
        source, "for-each-ref", "--format=%(refname)", "refs/replace")
    if replace_refs:
        raise RuntimeError(
            "candidate repository contains replace refs; benchmark identity "
            "requires an unmodified object graph")
    top = git_output(source, "rev-parse", "--show-toplevel")
    if os.path.realpath(top) != source:
        raise RuntimeError("candidate source is not the Git top level")
    head = git_output(source, "rev-parse", "HEAD")
    tree = git_output(source, "rev-parse", "HEAD^{tree}")
    expected_tree = git_output(source, "rev-parse", expected_commit + "^{tree}")
    for flag in ("-v", "-f"):
        records = [record for record in git_output(
            source, "ls-files", flag, "-z").split("\0")
                   if record]
        if not records or any(not record.startswith("H ")
                              for record in records):
            raise RuntimeError(
                "candidate index uses assume-unchanged, skip-worktree, "
                "fsmonitor-valid, or another non-default flag")
    status = git_output(
        source, "status", "--porcelain=v1", "--untracked-files=normal")
    if head != expected_commit or tree != expected_tree:
        raise RuntimeError(
            "candidate source does not exactly match expected commit")
    if status:
        raise RuntimeError("candidate source has tracked changes:\n" + status)
    return {"path": source, "head": head, "tree": tree,
            "status": "clean"}


def cmake_provenance(build, source, executable):
    return candidate_build_provenance(
        build, source, executable, "bench_leopard2")


def candidate_source_attestation(executable, source, snapshot_descriptor):
    cpu = sorted(os.sched_getaffinity(0))[0]
    command = [
        "/usr/bin/taskset", "-c", str(cpu), executable,
        "--k", "3", "--r", "2", "--bytes", "64", "--loss", "1",
        "--batch", "1", "--reuse", "1", "--iterations", "1",
        "--warmup", "0", "--threads", "1", "--seed", "7",
        "--profile", "high", "--field", "gf8", "--backend", "avx2",
        "--skip-legacy", "--retain-samples", "--attest-source",
        "--json", "-",
    ]
    completed = subprocess.run(
        command, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, env=benchmark_environment(), timeout=120,
        check=False, pass_fds=(snapshot_descriptor,))
    if completed.returncode != 0:
        raise InvocationError(
            "candidate source-attestation preflight failed",
            {"command": command, "returncode": completed.returncode,
             "stdout": completed.stdout.decode(errors="replace"),
             "stderr": completed.stderr.decode(errors="replace")})
    try:
        document = json.loads(completed.stdout)
    except Exception as error:
        raise InvocationError(
            "candidate source-attestation JSON is invalid: %s" % error,
            {"command": command,
             "stdout": completed.stdout.decode(errors="replace"),
             "stderr": completed.stderr.decode(errors="replace")})
    build = document.get("build", {})
    if document.get("schema") != "leopard2-benchmark-v5":
        raise RuntimeError("candidate attestation did not emit schema v5")
    if build.get("source_commit") != source["head"]:
        raise RuntimeError("candidate binary attests another source commit")
    if build.get("source_tree") != source["tree"]:
        raise RuntimeError("candidate binary attests another source tree")
    if build.get("source_tracked_dirty") is not False:
        raise RuntimeError("candidate binary was configured from dirty source")
    parameters = document["parameters"]
    resolved = document["resolved"]
    if (parameters["requested_profile"], parameters["requested_field"],
            parameters["requested_backend"]) != (
                "legacy_high_v1", "gf8", "avx2"):
        raise RuntimeError("candidate attestation request identity changed")
    if not parameters.get("skip_legacy") or \
            not parameters.get("attest_source"):
        raise RuntimeError("candidate attestation flags were not recorded")
    if (resolved["profile"], resolved["field"], resolved["backend"],
            resolved["parent_count"], resolved["padded_side"],
            resolved["thread_count"]) != (
                "legacy_high_v1", "gf8", "avx2", 8, 2, 1):
        raise RuntimeError("candidate attestation resolved identity changed")
    correctness = document["correctness"]
    if not correctness.get("leopard2_round_trip") or \
            correctness.get("legacy_comparison") is not None:
        raise RuntimeError("candidate attestation round trip failed")
    recorded_command = list(command)
    recorded_command[3] = "<sealed-candidate-snapshot>"
    stable = {
        "command": recorded_command,
        "schema": document["schema"], "build": build,
        "parameters": parameters, "resolved": resolved,
        "correctness": correctness,
        "workload_digests": document["workload_digests"],
    }
    stable["stable_identity_sha256"] = sha256_bytes(canonical_bytes(stable))
    return stable


def parse_cpu_list(text):
    result = set()
    for item in text.strip().split(","):
        if not item:
            continue
        bounds = item.split("-", 1)
        if len(bounds) == 1:
            result.add(int(bounds[0]))
        else:
            result.update(range(int(bounds[0]), int(bounds[1]) + 1))
    return result


def sibling_set(cpu):
    path = ("/sys/devices/system/cpu/cpu%d/topology/"
            "thread_siblings_list" % cpu)
    with open(path, encoding="ascii") as stream:
        return parse_cpu_list(stream.read())


def physical_cpu_pairs():
    allowed = set(os.sched_getaffinity(0))
    groups = {}
    for cpu in sorted(allowed):
        siblings = sibling_set(cpu) & allowed
        if len(siblings) < 2:
            continue
        key = tuple(sorted(siblings))
        groups[key] = key
    pairs = [(group[0], group[1]) for group in sorted(groups)]
    if not pairs:
        raise RuntimeError("no allowed SMT sibling pairs were detected")
    return pairs


def parse_pairs(text):
    if not text:
        return physical_cpu_pairs()
    pairs = []
    for item in text.split(","):
        values = item.split(":")
        if len(values) != 2:
            raise RuntimeError("CPU pairs use cpu:sibling syntax")
        pairs.append((int(values[0]), int(values[1])))
    for cpu, sibling in pairs:
        validate_cpu_pair(cpu, sibling)
    if len({value for pair in pairs for value in pair}) != 2 * len(pairs):
        raise RuntimeError("CPU pairs overlap")
    return pairs


def validate_cpu_pair(cpu, sibling):
    allowed = set(os.sched_getaffinity(0))
    if cpu == sibling or cpu not in allowed or sibling not in allowed:
        raise RuntimeError("benchmark CPU pair is not distinct and allowed")
    siblings = sibling_set(cpu) & allowed
    if cpu not in siblings or sibling not in siblings:
        raise RuntimeError("CPU%d and CPU%d are not SMT siblings" %
                           (cpu, sibling))
    if (sibling_set(sibling) & allowed) != siblings:
        raise RuntimeError("asymmetric sibling topology")
    return {"benchmark_cpu": cpu, "reserved_sibling": sibling,
            "thread_siblings": sorted(siblings),
            "process_allowed_cpus": sorted(allowed)}


def cpu_stat(cpu):
    prefix = "cpu%d " % cpu
    with open("/proc/stat", encoding="ascii") as stream:
        for line in stream:
            if line.startswith(prefix):
                values = [int(value) for value in line.split()[1:]]
                while len(values) < 8:
                    values.append(0)
                names = ("user", "nice", "system", "idle", "iowait",
                         "irq", "softirq", "steal")
                return dict(zip(names, values[:8]))
    raise RuntimeError("missing CPU%d in /proc/stat" % cpu)


def nonidle_delta(before, after):
    names = ("user", "nice", "system", "irq", "softirq", "steal")
    return sum(after[name] - before[name] for name in names)


@contextlib.contextmanager
def benchmark_lock(path, exclusive):
    descriptor = os.open(path, os.O_CREAT | os.O_RDWR, 0o600)
    try:
        mode = fcntl.LOCK_EX if exclusive else fcntl.LOCK_SH
        print("waiting for %s benchmark lock %s" %
              ("exclusive" if exclusive else "shared", path), flush=True)
        fcntl.flock(descriptor, mode)
        print("acquired benchmark lock", flush=True)
        yield
    finally:
        fcntl.flock(descriptor, fcntl.LOCK_UN)
        os.close(descriptor)


def mem_available_bytes():
    with open("/proc/meminfo", encoding="ascii") as stream:
        for line in stream:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) * 1024
    raise RuntimeError("MemAvailable is unavailable")


def mem_total_bytes():
    with open("/proc/meminfo", encoding="ascii") as stream:
        for line in stream:
            if line.startswith("MemTotal:"):
                return int(line.split()[1]) * 1024
    raise RuntimeError("MemTotal is unavailable")


class InvocationError(RuntimeError):
    def __init__(self, message, details):
        RuntimeError.__init__(self, message)
        self.details = details


def benchmark_environment():
    return {
        "PATH": "/usr/bin:/bin", "LANG": "C", "LC_ALL": "C",
        "TZ": "UTC", "OMP_DYNAMIC": "FALSE", "OMP_NUM_THREADS": "1",
        "OMP_PLACES": "cores", "OMP_PROC_BIND": "TRUE",
    }


def validate_baseline_contract(baseline_contract):
    if baseline_contract not in (
            BASELINE_CONTRACT_LEGACY_V1,
            BASELINE_CONTRACT_LOGICAL_ZERO_PAD_V2):
        raise RuntimeError("unknown exact-main baseline contract")


def physical_shard_bytes(kind, logical_bytes, baseline_contract):
    if kind not in ("baseline", "candidate"):
        raise RuntimeError("unknown benchmark implementation kind")
    validate_baseline_contract(baseline_contract)
    if (isinstance(logical_bytes, bool) or
            not isinstance(logical_bytes, int) or logical_bytes < 1):
        raise RuntimeError("logical shard bytes must be a positive integer")
    if kind == "baseline":
        if baseline_contract == BASELINE_CONTRACT_LOGICAL_ZERO_PAD_V2:
            return (logical_bytes + 63) & ~63
        if logical_bytes & 63:
            raise RuntimeError(
                "legacy exact-main v1 requires 64-byte-aligned shards")
    return logical_bytes


def benchmark_byte_arguments(kind, logical_bytes, baseline_contract):
    physical_bytes = physical_shard_bytes(
        kind, logical_bytes, baseline_contract)
    result = ["--bytes", str(physical_bytes)]
    if (kind == "baseline" and
            baseline_contract == BASELINE_CONTRACT_LOGICAL_ZERO_PAD_V2 and
            physical_bytes != logical_bytes):
        result += ["--logical-bytes", str(logical_bytes)]
    return result


def expected_benchmark_byte_identity(
        kind, logical_bytes, baseline_contract):
    physical_bytes = physical_shard_bytes(
        kind, logical_bytes, baseline_contract)
    return {
        "baseline_contract": baseline_contract,
        "logical_shard_bytes": logical_bytes,
        "physical_shard_bytes": physical_bytes,
        "padded_application_bytes": physical_bytes != logical_bytes,
        "ratio_scope": "time_per_call",
        "digest_scope": "logical_prefix",
    }


def validate_benchmark_byte_identity(
        kind, document, cell, baseline_contract):
    logical_bytes = cell["bytes"]
    physical_bytes = physical_shard_bytes(
        kind, logical_bytes, baseline_contract)
    parameters = document.get("parameters")
    resolved = document.get("resolved")
    if not isinstance(parameters, dict) or not isinstance(resolved, dict):
        raise RuntimeError("benchmark byte identity is missing")
    if parameters.get("shard_bytes") != physical_bytes:
        raise RuntimeError("benchmark physical shard bytes mismatch")
    if kind == "baseline":
        if baseline_contract == BASELINE_CONTRACT_LEGACY_V1:
            correctness = document.get("correctness")
            if (document.get("schema") != "leopard-main-benchmark-v1" or
                    "logical_shard_bytes" in parameters or
                    "padded_application_bytes" in resolved or
                    "padding_policy" in resolved or
                    not isinstance(correctness, dict) or
                    "logical_prefix_fingerprinted" in correctness):
                raise RuntimeError(
                    "legacy v1 baseline byte identity mismatch")
        else:
            padded = physical_bytes != logical_bytes
            expected_schema = (
                "leopard-main-benchmark-v2" if padded else
                "leopard-main-benchmark-v1")
            correctness = document.get("correctness")
            if (document.get("schema") != expected_schema or
                    parameters.get("logical_shard_bytes") != logical_bytes or
                    resolved.get("padded_application_bytes") is not padded or
                    resolved.get("padding_policy") !=
                        "zero suffix per shard" or
                    not isinstance(correctness, dict) or
                    correctness.get(
                        "logical_prefix_fingerprinted") is not True):
                raise RuntimeError(
                    "baseline logical-byte padding identity mismatch")
    elif "logical_shard_bytes" in parameters:
        raise RuntimeError(
            "candidate unexpectedly reported logical-byte adaptation")
    return expected_benchmark_byte_identity(
        kind, logical_bytes, baseline_contract)


def invoke_benchmark(kind, executable, cpu, sibling, cell, iterations,
                     warmup, maximum_attempts, candidate_commit,
                     require_idle_sibling, snapshot_descriptor,
                     baseline_contract):
    common = [
        "--k", str(cell["K"]), "--r", str(cell["R"]),
    ] + benchmark_byte_arguments(
        kind, cell["bytes"], baseline_contract) + [
        "--loss", str(cell["loss"]),
        "--batch", str(cell["batch"]), "--reuse", str(cell["reuse"]),
        "--iterations", str(iterations), "--warmup", str(warmup),
        "--threads", "1", "--seed", str(cell["seed"]), "--json", "-",
    ]
    if kind == "candidate":
        common += [
            "--profile", "high", "--field", "gf8", "--backend", "avx2",
            "--skip-legacy", "--retain-samples", "--report-decode-path",
        ]
    command = ["/usr/bin/taskset", "-c", str(cpu), executable] + common
    rejected = []
    timeout = 1200 if cell["bytes"] >= 1048576 else 300
    for attempt in range(maximum_attempts):
        cpu_before = cpu_stat(cpu)
        sibling_before = cpu_stat(sibling)
        started_ns = time.monotonic_ns()
        try:
            process = subprocess.run(
                command, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, env=benchmark_environment(),
                timeout=timeout, check=False,
                pass_fds=(snapshot_descriptor,))
        except subprocess.TimeoutExpired as error:
            details = {
                "kind": kind, "command": command, "timeout_seconds": timeout,
                "stdout": (error.stdout or b"").decode(errors="replace"),
                "stderr": (error.stderr or b"").decode(errors="replace"),
            }
            raise InvocationError("benchmark timed out", details)
        ended_ns = time.monotonic_ns()
        cpu_after = cpu_stat(cpu)
        sibling_after = cpu_stat(sibling)
        sibling_work = nonidle_delta(sibling_before, sibling_after)
        if process.returncode != 0:
            details = {
                "kind": kind, "command": command,
                "returncode": process.returncode,
                "stdout": process.stdout.decode(errors="replace"),
                "stderr": process.stderr.decode(errors="replace"),
            }
            raise InvocationError("benchmark child failed", details)
        if sibling_work and require_idle_sibling:
            rejected.append({
                "attempt": attempt,
                "reason": "reserved sibling performed non-idle work",
                "sibling_nonidle_jiffies": sibling_work,
                "started_ns": started_ns, "ended_ns": ended_ns,
                "stdout_sha256": hashlib.sha256(process.stdout).hexdigest(),
                "stderr_sha256": hashlib.sha256(process.stderr).hexdigest(),
            })
            time.sleep(0.05)
            continue
        try:
            document = json.loads(process.stdout)
        except Exception as error:
            details = {
                "kind": kind, "command": command,
                "stdout": process.stdout.decode(errors="replace"),
                "stderr": process.stderr.decode(errors="replace"),
            }
            raise InvocationError("invalid benchmark JSON: %s" % error,
                                  details)
        parameters = document["parameters"]
        byte_identity = validate_benchmark_byte_identity(
            kind, document, cell, baseline_contract)
        observed = (
            parameters["K"], parameters["R"], parameters["shard_bytes"],
            parameters["loss_count"], parameters["batch"],
            parameters["reuse"], parameters["seed"])
        expected = (
            cell["K"], cell["R"],
            physical_shard_bytes(
                kind, cell["bytes"], baseline_contract), cell["loss"],
            cell["batch"], cell["reuse"], cell["seed"])
        if observed != expected:
            raise RuntimeError("benchmark parameter mismatch for " + cell["id"])
        resolved = document["resolved"]
        if (resolved["profile"], resolved["field"]) != (
                "legacy_high_v1", "gf8"):
            raise RuntimeError("wire profile or field mismatch")
        if parameters["thread_count"] != 1 or resolved["thread_count"] != 1:
            raise RuntimeError("benchmark did not remain single-threaded")
        if (resolved["parent_count"], resolved["padded_side"]) != (
                ceil_pow2(cell["K"] + cell["T"]), cell["T"]):
            raise RuntimeError("benchmark parent geometry mismatch")
        if kind == "baseline":
            if document["build"]["main_source_commit"] != MAIN_COMMIT:
                raise RuntimeError("baseline commit mismatch")
            if not document["correctness"]["round_trip"]:
                raise RuntimeError("baseline round trip failed")
            encode = document["metrics"]["encode_execution"][
                "median_us_per_batch_call"]
            decode = document["metrics"]["decode_including_setup"][
                "median_us_per_batch_call"]
            times = {"encode": encode, "decode_first": decode,
                     "decode_reuse": decode}
        else:
            if parameters["requested_backend"] != "avx2" or \
                    resolved["backend"] != "avx2":
                raise RuntimeError("candidate did not resolve explicit AVX2")
            runtime_commit = document.get("build", {}).get("source_commit")
            if runtime_commit is not None and runtime_commit != candidate_commit:
                raise RuntimeError("candidate runtime commit mismatch")
            if not document["correctness"]["leopard2_round_trip"]:
                raise RuntimeError("candidate round trip failed")
            encode = document["metrics"]["encode_execution"][
                "median_us_per_batch_call"]
            execution = document["metrics"]["decode_execution"][
                "median_us_per_batch_call"]
            setup = document["metrics"]["decode_plan_setup"]["median_us"]
            times = {
                "encode": encode,
                "decode_first": execution + setup,
                "decode_reuse": execution + setup / cell["reuse"],
            }
        if any(isinstance(value, bool) or
               not isinstance(value, (int, float)) or
               not math.isfinite(value) or value <= 0
               for value in times.values()):
            raise RuntimeError("benchmark timing is not positive and finite")
        recorded_command = list(command)
        recorded_command[3] = "<sealed-%s-snapshot>" % kind
        return {
            "kind": kind, "attempt": attempt, "command": recorded_command,
            "started_ns": started_ns, "ended_ns": ended_ns,
            "cpu_nonidle_jiffies": nonidle_delta(cpu_before, cpu_after),
            "sibling_nonidle_jiffies": sibling_work,
            "idle_sibling_required": require_idle_sibling,
            "stdout_sha256": hashlib.sha256(process.stdout).hexdigest(),
            "stderr_sha256": hashlib.sha256(process.stderr).hexdigest(),
            "rejected_attempts": rejected,
            "times_us": times, "byte_identity": byte_identity,
            "document": document,
        }
    raise InvocationError(
        "reserved sibling contamination persisted",
        {"kind": kind, "command": command, "rejected_attempts": rejected})


def digest_tuple(document):
    value = document["workload_digests"]
    return (value["algorithm"], value["original_data"],
            value["transmitted_parity"], value["recovered_originals"])


def validate_call_identity(calls, cell):
    expected_digest = None
    expected_missing = None
    required_missing_index = cell.get("expected_missing_original_index")
    required_missing_indices = cell.get("expected_missing_original_indices")
    if (required_missing_index is not None and
            required_missing_indices is not None):
        raise RuntimeError("cell declares two expected missing-set formats")
    if required_missing_indices is not None:
        if (not isinstance(required_missing_indices, list) or
                any(isinstance(index, bool) or not isinstance(index, int)
                    for index in required_missing_indices)):
            raise RuntimeError("expected missing original set is malformed")
        required_missing = tuple(required_missing_indices)
        if (required_missing != tuple(sorted(set(required_missing))) or
                len(required_missing) != cell["loss"] or
                any(index < 0 or index >= cell["K"]
                    for index in required_missing)):
            raise RuntimeError("expected missing original set is invalid")
    elif required_missing_index is not None:
        required_missing = (required_missing_index,)
    else:
        required_missing = None
    candidate_route = None
    for call in calls:
        document = call["document"]
        digest = digest_tuple(document)
        missing = tuple(document["parameters"]["missing_original_indices"])
        if required_missing is not None and missing != required_missing:
            message = "expected missing original set mismatch for "
            if required_missing_index is not None:
                message = "expected missing original index mismatch for "
            raise RuntimeError(message + cell["id"])
        if expected_digest is None:
            expected_digest, expected_missing = digest, missing
        if digest != expected_digest or missing != expected_missing:
            raise RuntimeError(
                "wire/recovery digest or missing-set mismatch for " + cell["id"])
        if len(missing) != cell["loss"]:
            raise RuntimeError("missing-set cardinality mismatch")
        if call["kind"] == "candidate":
            route = (document["resolved"]["selected_decode_path"],
                     document["resolved"]["selected_decode_rule"])
            if candidate_route is None:
                candidate_route = route
            elif route != candidate_route:
                raise RuntimeError("candidate decode route changed")
    if expected_digest is None or candidate_route is None:
        raise RuntimeError("incomplete cell invocation set")
    return {
        "workload_digests": {
            "algorithm": expected_digest[0],
            "original_data": expected_digest[1],
            "transmitted_parity": expected_digest[2],
            "recovered_originals": expected_digest[3],
        },
        "missing_original_indices": list(expected_missing),
        "selected_decode_path": candidate_route[0],
        "selected_decode_rule": candidate_route[1],
    }


def point_ratios(calls):
    by_kind = {}
    for call in calls:
        kind = call.get("kind")
        if kind not in ("baseline", "candidate") or kind in by_kind:
            raise RuntimeError(
                "point ratios require one baseline and one candidate call")
        by_kind[kind] = call
    if set(by_kind) != {"baseline", "candidate"}:
        raise RuntimeError(
            "point ratios require one baseline and one candidate call")
    baseline = by_kind["baseline"]
    candidate = by_kind["candidate"]
    for kind, call in by_kind.items():
        times = call.get("times_us")
        if not isinstance(times, dict) or set(times) != set(RATIO_METRICS):
            raise RuntimeError(kind + " point timings have wrong metrics")
        if any(isinstance(value, bool) or
               not isinstance(value, (int, float)) or
               not math.isfinite(value) or value <= 0
               for value in times.values()):
            raise RuntimeError(
                kind + " point timings are not positive and finite")
    return {
        metric: baseline["times_us"][metric] / candidate["times_us"][metric]
        for metric in RATIO_METRICS
    }


def confidence(log_values):
    count = len(log_values)
    if count < 2:
        raise RuntimeError("confidence interval needs at least two rounds")
    if any(isinstance(value, bool) or
           not isinstance(value, (int, float)) or
           not math.isfinite(value) for value in log_values):
        raise RuntimeError("confidence inputs must be finite real numbers")
    critical_table = {2: 12.706204736, 3: 4.302652730, 4: 3.182446305,
                      5: 2.776445105, 6: 2.570581836, 7: 2.446911851,
                      8: 2.364624252, 9: 2.306004135, 10: 2.262157163}
    critical = critical_table.get(count, 1.96)
    mean = statistics.mean(log_values)
    half = critical * statistics.stdev(log_values) / math.sqrt(count)
    return {
        "ratio": math.exp(mean), "ci95_low": math.exp(mean - half),
        "ci95_high": math.exp(mean + half),
        "round_log_contrasts": log_values,
        "orientation": "leopard1_time_over_leopard2_time",
    }


def abba_metrics(calls, rounds):
    if len(calls) != rounds * len(ABBA_ORDER):
        raise RuntimeError("ABBA call count mismatch")
    result = {}
    for metric in RATIO_METRICS:
        contrasts = []
        for round_index in range(rounds):
            group = calls[round_index * 4:(round_index + 1) * 4]
            baseline = [math.log(call["times_us"][metric]) for call in group
                        if call["kind"] == "baseline"]
            candidate = [math.log(call["times_us"][metric]) for call in group
                         if call["kind"] == "candidate"]
            contrasts.append(statistics.mean(baseline) -
                             statistics.mean(candidate))
        result[metric] = confidence(contrasts)
    return result


def utc_now():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def machine_provenance():
    lscpu = json.loads(checked_output(["/usr/bin/lscpu", "-J"]))
    # Current CPU-frequency percentage is an observation, not platform
    # identity.  Excluding it makes an interrupted run safely resumable.
    lscpu["lscpu"] = [
        item for item in lscpu["lscpu"]
        if item.get("field", "").strip() != "CPU(s) scaling MHz:"]
    return {
        "platform": platform.platform(), "uname": list(platform.uname()),
        "python": sys.version, "python_executable": sys.executable,
        "allowed_cpus": sorted(os.sched_getaffinity(0)),
        "lscpu_json": lscpu,
        "mem_total_bytes": mem_total_bytes(),
    }


def matrix_baseline_sha256(matrix):
    observed = matrix.get("baseline_sha256")
    if (not isinstance(observed, str) or
            SHA256_PATTERN.fullmatch(observed) is None):
        raise RuntimeError("matrix baseline SHA-256 is invalid")
    profile = matrix.get("profile", CANONICAL_MATRIX_PROFILE)
    if profile in (CANONICAL_MATRIX_PROFILE, R1_TINY_EXACT_MAIN_PROFILE):
        required = BASELINE_SHA256
    elif profile == K4_R2_PRODUCTION_AVX2_PROFILE:
        required = PADDED_MAIN_BASELINE_SHA256
    else:
        raise RuntimeError("unknown matrix profile baseline identity")
    if observed != required:
        raise RuntimeError("matrix baseline SHA-256 differs from its profile")
    return required


def matrix_baseline_contract(matrix):
    baseline_sha256 = matrix_baseline_sha256(matrix)
    if baseline_sha256 == BASELINE_SHA256:
        return BASELINE_CONTRACT_LEGACY_V1
    if baseline_sha256 == PADDED_MAIN_BASELINE_SHA256:
        return BASELINE_CONTRACT_LOGICAL_ZERO_PAD_V2
    raise RuntimeError("matrix baseline has no byte contract")


def common_provenance(args, matrix):
    expected_baseline_sha256 = matrix_baseline_sha256(matrix)
    baseline_contract = matrix_baseline_contract(matrix)
    source = source_provenance(args.candidate_source, args.candidate_commit)
    build_executable = (
        args.candidate_build_executable or args.candidate)
    build = cmake_provenance(
        args.candidate_build, args.candidate_source, build_executable)
    reproducible_build = verify_reproducible_candidate_build(build)
    executable_equivalence = equivalent_executable_identities(
        build_executable, args.candidate)
    descriptors = []
    try:
        baseline_source, baseline_descriptor, baseline_snapshot = \
            snapshot_executable(
                args.baseline, "baseline benchmark",
                expected_baseline_sha256)
        descriptors.append(baseline_descriptor)
        candidate_source, candidate_descriptor, candidate_snapshot = \
            snapshot_executable(args.candidate, "candidate benchmark")
        descriptors.append(candidate_descriptor)
        if candidate_source != \
                executable_equivalence["selected_timing_artifact"]:
            raise RuntimeError(
                "candidate timing artifact changed before snapshot")
        baseline_path = snapshot_path(baseline_descriptor)
        candidate_path = snapshot_path(candidate_descriptor)
        attestation = candidate_source_attestation(
            candidate_path, source, candidate_descriptor)
        provenance = {
            "schema": SCHEMA,
            "runner": {"path": os.path.realpath(__file__),
                       "sha256": sha256_file(__file__)},
            "matrix": {"path": os.path.realpath(args.matrix),
                       "sha256": matrix_digest(matrix),
                       "file_sha256": sha256_file(args.matrix),
                       "cell_count": len(matrix["cells"])},
            "frozen_main_commit": MAIN_COMMIT,
            "baseline_contract": baseline_contract,
            "baseline": {
                "path": baseline_source["path"],
                "sha256": baseline_snapshot["sha256"],
                "expected_sha256": expected_baseline_sha256,
                "source_file": baseline_source,
                "snapshot": baseline_snapshot,
                "isa": inspect_isa(
                    baseline_path, (baseline_descriptor,)),
            },
            "candidate": {
                "path": candidate_source["path"],
                "sha256": candidate_snapshot["sha256"],
                "source_file": candidate_source,
                "snapshot": candidate_snapshot,
                "isa": inspect_isa(
                    candidate_path, (candidate_descriptor,)),
                "source_attestation": attestation,
            },
            "candidate_source": source,
            "candidate_build": build,
            "candidate_executable_equivalence": executable_equivalence,
            "candidate_reproducible_build": reproducible_build,
            "machine": machine_provenance(),
        }
        runtime = {
            "baseline_contract": baseline_contract,
            "baseline": {
                "path": baseline_path, "descriptor": baseline_descriptor},
            "candidate": {
                "path": candidate_path, "descriptor": candidate_descriptor},
        }
        return provenance, runtime
    except BaseException:
        for descriptor in descriptors:
            os.close(descriptor)
        raise


def stage_identity(stage, provenance, options):
    return {"schema": SCHEMA, "stage": stage, "provenance": provenance,
            "options": options}


def prepare_run_directory(run_dir, identity):
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    identity_path = run_dir / "identity.json"
    if identity_path.exists():
        existing = read_json(identity_path)
        if existing != identity:
            raise RuntimeError(
                "run directory identity differs; use a new run directory")
    else:
        atomic_json(identity_path, identity)
    (run_dir / "cells").mkdir(exist_ok=True)
    (run_dir / "failures").mkdir(exist_ok=True)
    return run_dir, sha256_bytes(canonical_bytes(identity))


def close_runtime_snapshots(runtime):
    for role in ("baseline", "candidate"):
        descriptor = runtime.get(role, {}).get("descriptor")
        if isinstance(descriptor, int) and descriptor >= 0:
            os.close(descriptor)


def source_still_frozen(provenance, runtime):
    expected = provenance["candidate_source"]
    observed = source_provenance(expected["path"], expected["head"])
    if observed != expected:
        raise RuntimeError("candidate source provenance changed during stage")
    expected_build = provenance["candidate_build"]
    observed_build = candidate_build_provenance(
        expected_build["build_root"], expected_build["source_root"],
        expected_build["executable"]["path"],
        expected_build["executable_target"])
    if observed_build != expected_build:
        raise RuntimeError(
            "candidate source/object/archive/link closure changed during stage")
    expected_equivalence = provenance["candidate_executable_equivalence"]
    observed_equivalence = equivalent_executable_identities(
        expected_equivalence["build_tree_executable"]["path"],
        expected_equivalence["selected_timing_artifact"]["path"])
    if observed_equivalence != expected_equivalence:
        raise RuntimeError(
            "candidate build executable or timing artifact changed during stage")
    for role in ("baseline", "candidate"):
        recorded = provenance[role]
        if file_identity(recorded["path"], role + " benchmark") != \
                recorded["source_file"]:
            raise RuntimeError(role + " benchmark source changed during stage")
        observed_snapshot = sealed_snapshot_identity(
            runtime[role]["descriptor"], role + " benchmark")
        if observed_snapshot != recorded["snapshot"]:
            raise RuntimeError(role + " sealed snapshot changed during stage")


def completed_result(path, identity_digest, cell, stage):
    if not path.exists():
        return None
    value = read_json(path)
    if value.get("status") != "complete":
        return None
    if (value.get("identity_sha256") != identity_digest or
            value.get("stage") != stage or value.get("cell") != cell):
        raise RuntimeError("resumed cell identity mismatch: " + str(path))
    return value


def retain_failure(run_dir, stage, cell, pair, identity_digest, error):
    stamp = "%d-%d" % (time.time_ns(), os.getpid())
    details = getattr(error, "details", None)
    value = {
        "schema": SCHEMA, "status": "failed", "stage": stage,
        "identity_sha256": identity_digest, "cell": cell,
        "cpu_pair": list(pair), "error_type": type(error).__name__,
        "error": str(error), "details": details, "failed_at": utc_now(),
    }
    atomic_json(Path(run_dir) / "failures" /
                (cell["id"] + "-" + stamp + ".json"), value)


class MemoryGate:
    def __init__(self, budget):
        self.budget = budget
        self.used = 0
        self.condition = threading.Condition()

    @contextlib.contextmanager
    def hold(self, requested):
        amount = min(requested, self.budget)
        with self.condition:
            while self.used + amount > self.budget:
                self.condition.wait()
            self.used += amount
        try:
            yield
        finally:
            with self.condition:
                self.used -= amount
                self.condition.notify_all()


def diagnostic_cell(cell, pair, args, identity_digest, run_dir, memory_gate,
                    runtime):
    result_path = Path(run_dir) / "cells" / (cell["id"] + ".json")
    resumed = completed_result(result_path, identity_digest, cell, "diagnostic")
    if resumed is not None:
        return resumed, True
    cpu, sibling = pair
    try:
        with memory_gate.hold(cell["estimated_peak_bytes"]):
            order = ("baseline", "candidate")
            if cell["seed"] & 1:
                order = tuple(reversed(order))
            calls = []
            for kind in order:
                executable = runtime[kind]["path"]
                calls.append(invoke_benchmark(
                    kind, executable, cpu, sibling, cell,
                    cell["diagnostic_iterations"],
                    cell["diagnostic_warmup"], args.maximum_attempts,
                    args.candidate_commit, False,
                    runtime[kind]["descriptor"],
                    runtime["baseline_contract"]))
        identity = validate_call_identity(calls, cell)
        value = {
            "schema": SCHEMA, "status": "complete", "stage": "diagnostic",
            "authoritative": False, "identity_sha256": identity_digest,
            "cell": cell, "cpu_pair": list(pair), "calls": calls,
            "ratios": point_ratios(calls), **identity,
        }
        atomic_json(result_path, value)
        return value, False
    except Exception as error:
        retain_failure(run_dir, "diagnostic", cell, pair, identity_digest, error)
        raise


def diagnostic_summary(run_dir, identity, identity_digest, matrix, results,
                       resumed_count, elapsed_seconds):
    rows = sorted(results, key=lambda value: value["cell"]["id"])
    summary = {
        "schema": SCHEMA, "status": "complete", "stage": "diagnostic",
        "authoritative": False,
        "note": "Parallel screen only; ratios are not publication evidence.",
        "identity_sha256": identity_digest,
        "matrix_sha256": identity["provenance"]["matrix"]["sha256"],
        "cell_count": len(rows), "resumed_cell_count": resumed_count,
        "accepted_child_invocations": 2 * len(rows),
        "accepted_calls_with_sibling_activity": sum(
            int(call["sibling_nonidle_jiffies"] != 0)
            for row in results for call in row["calls"]),
        "elapsed_seconds": elapsed_seconds,
        "candidate_binary_sha256": identity["provenance"]["candidate"]["sha256"],
        "baseline_binary_sha256": identity["provenance"]["baseline"]["sha256"],
        "rows": [{
            "id": row["cell"]["id"], "cell": row["cell"],
            "ratios": row["ratios"], "cpu_pair": row["cpu_pair"],
            "selected_decode_path": row["selected_decode_path"],
            "selected_decode_rule": row["selected_decode_rule"],
            "workload_digests": row["workload_digests"],
            "benchmark_byte_identity": {
                kind: expected_benchmark_byte_identity(
                    kind, row["cell"]["bytes"],
                    matrix_baseline_contract(matrix))
                for kind in ("baseline", "candidate")
            },
            "missing_original_indices": row["missing_original_indices"],
        } for row in rows],
    }
    atomic_json(Path(run_dir) / "summary.json", summary)
    return summary


def run_diagnostic(args):
    matrix = read_json(args.matrix)
    validate_matrix(matrix)
    pairs = parse_pairs(args.cpu_pairs)
    if args.maximum_workers:
        pairs = pairs[:args.maximum_workers]
    if not pairs:
        raise RuntimeError("no CPU pairs selected")
    budget = args.memory_budget_gib * (1 << 30) if args.memory_budget_gib else \
        min(mem_total_bytes() // 4, 32 << 30)
    with benchmark_lock(args.lock, True):
        provenance, runtime = common_provenance(args, matrix)
        try:
            options = {
                "cpu_pairs": [list(pair) for pair in pairs],
                "maximum_attempts": args.maximum_attempts,
                "memory_budget_bytes": budget,
                "mode": "exclusive_parallel_non_authoritative_screen",
                "idle_sibling_required": False,
            }
            identity = stage_identity("diagnostic", provenance, options)
            run_dir, identity_digest = prepare_run_directory(
                args.run_dir, identity)
            memory_gate = MemoryGate(budget)
            started = time.monotonic()
            ordered_cells = sorted(
                matrix["cells"],
                key=lambda cell: (cell["estimated_peak_bytes"], cell["id"]))
            assignments = [[] for unused_pair in pairs]
            for index, cell in enumerate(ordered_cells):
                assignments[index % len(pairs)].append(cell)
            progress_lock = threading.Lock()
            progress = [0]

            def worker(pair, assigned_cells):
                local_results = []
                local_resumed = 0
                local_errors = []
                for cell in assigned_cells:
                    failed = False
                    try:
                        result, resumed = diagnostic_cell(
                            cell, pair, args, identity_digest, run_dir,
                            memory_gate, runtime)
                        local_results.append(result)
                        local_resumed += int(resumed)
                    except Exception as error:
                        failed = True
                        local_errors.append((cell["id"], error))
                    with progress_lock:
                        progress[0] += 1
                        print("diagnostic %d/%d %s%s" % (
                            progress[0], len(matrix["cells"]), cell["id"],
                            " FAILED" if failed else ""), flush=True)
                return local_results, local_resumed, local_errors

            results = []
            resumed_count = 0
            errors = []
            with concurrent.futures.ThreadPoolExecutor(
                    max_workers=len(pairs),
                    thread_name_prefix="gf8-map") as executor:
                futures = [
                    executor.submit(worker, pair, assignment)
                    for pair, assignment in zip(pairs, assignments)]
                for future in concurrent.futures.as_completed(futures):
                    worker_results, worker_resumed, worker_errors = \
                        future.result()
                    results.extend(worker_results)
                    resumed_count += worker_resumed
                    errors.extend(worker_errors)
            if errors:
                raise RuntimeError(
                    "%d diagnostic cells failed; first %s: %s" %
                    (len(errors), errors[0][0], errors[0][1]))
            source_still_frozen(provenance, runtime)
            summary = diagnostic_summary(
                run_dir, identity, identity_digest, matrix, results,
                resumed_count, time.monotonic() - started)
            print(json.dumps({
                "stage": "diagnostic", "run_dir": str(run_dir),
                "cells": summary["cell_count"],
                "elapsed_seconds": summary["elapsed_seconds"],
                "summary_sha256": sha256_file(run_dir / "summary.json"),
            }, sort_keys=True))
        finally:
            close_runtime_snapshots(runtime)


def validate_ratio_map(ratios, label):
    if not isinstance(ratios, dict) or set(ratios) != set(RATIO_METRICS):
        raise RuntimeError(label + " has the wrong ratio metrics")
    if any(isinstance(value, bool) or
           not isinstance(value, (int, float)) or
           not math.isfinite(value) or value <= 0
           for value in ratios.values()):
        raise RuntimeError(label + " ratios are not positive and finite")


def validate_diagnostic_summary(summary, matrix):
    if (summary.get("schema") != SCHEMA or
            summary.get("stage") != "diagnostic" or
            summary.get("status") != "complete" or
            summary.get("authoritative") is not False):
        raise RuntimeError("invalid diagnostic summary")
    if summary.get("matrix_sha256") != matrix_digest(matrix):
        raise RuntimeError("diagnostic summary matrix mismatch")
    rows = summary.get("rows")
    if not isinstance(rows, list):
        raise RuntimeError("diagnostic summary rows are missing")
    expected = {cell["id"]: cell for cell in matrix["cells"]}
    if summary.get("cell_count") != len(rows):
        raise RuntimeError("diagnostic summary row count mismatch")
    ids = [row.get("id") for row in rows if isinstance(row, dict)]
    if (len(ids) != len(rows) or len(ids) != len(expected) or
            set(ids) != set(expected)):
        raise RuntimeError("diagnostic summary is incomplete")
    for row in rows:
        identifier = row["id"]
        if row.get("cell") != expected[identifier]:
            raise RuntimeError(
                "diagnostic summary cell identity mismatch: " + identifier)
        validate_ratio_map(
            row.get("ratios"), "diagnostic summary row " + identifier)
    for field in ("baseline_binary_sha256", "candidate_binary_sha256"):
        if not isinstance(summary.get(field), str) or \
                not SHA256_PATTERN.fullmatch(summary[field]):
            raise RuntimeError(
                "diagnostic summary has invalid " + field)
    return summary


def load_diagnostic_summary(path, matrix):
    return validate_diagnostic_summary(read_json(path), matrix)


def selected_abba_cells(matrix, diagnostic, threshold, include, exclude):
    row_by_id = {row["id"]: row for row in diagnostic["rows"]}
    include = set(include)
    exclude = set(exclude)
    unknown = (include | exclude) - set(row_by_id)
    if unknown:
        raise RuntimeError("unknown selected cell IDs: " + ", ".join(sorted(unknown)))
    mandatory_tags = {
        "loss_sweep", "reuse_sweep", "batch_sweep", "large_tiling",
        "selector_isolation", "r1_selector_fallback",
        "r1_selector_promoted", "r1_large_recovery", "direct_odd_arity",
        "r1_tiny_exact_main", "k4_r2_production_avx2_gate",
    }
    selected = []
    for cell in matrix["cells"]:
        reasons = []
        ratios = row_by_id[cell["id"]]["ratios"]
        near_metrics = sorted(
            metric for metric, ratio in ratios.items() if ratio <= threshold)
        if near_metrics:
            reasons.append("diagnostic_ratio_at_or_below_%g:%s" %
                           (threshold, ",".join(near_metrics)))
        tags = sorted(set(cell["tags"]) & mandatory_tags)
        if tags:
            reasons.append("mandatory_dimension:" + ",".join(tags))
        if cell["id"] in include:
            reasons.append("explicit_include")
        if reasons and cell["id"] not in exclude:
            value = dict(cell)
            value["abba_selection_reasons"] = reasons
            selected.append(value)
    return selected


def abba_cell(cell, pair, args, identity_digest, run_dir, runtime):
    result_path = Path(run_dir) / "cells" / (cell["id"] + ".json")
    resumed = completed_result(result_path, identity_digest, cell, "abba")
    if resumed is not None:
        return resumed, True
    cpu, sibling = pair
    calls = []
    try:
        for unused_round in range(args.rounds):
            for kind in ABBA_ORDER:
                executable = runtime[kind]["path"]
                calls.append(invoke_benchmark(
                    kind, executable, cpu, sibling, cell,
                    cell["abba_iterations"], cell["abba_warmup"],
                    args.maximum_attempts, args.candidate_commit, True,
                    runtime[kind]["descriptor"],
                    runtime["baseline_contract"]))
        exact = validate_call_identity(calls, cell)
        value = {
            "schema": SCHEMA, "status": "complete", "stage": "abba",
            "authoritative": True, "identity_sha256": identity_digest,
            "cell": cell, "cpu_pair": list(pair), "rounds": args.rounds,
            "order_per_round": list(ABBA_ORDER), "calls": calls,
            "metrics": abba_metrics(calls, args.rounds), **exact,
        }
        atomic_json(result_path, value)
        return value, False
    except Exception as error:
        retain_failure(run_dir, "abba", cell, pair, identity_digest, error)
        raise


def run_abba(args):
    matrix = read_json(args.matrix)
    validate_matrix(matrix)
    if (args.cpu is None) != (args.sibling is None):
        raise RuntimeError("--cpu and --sibling must be supplied together")
    pair = ((args.cpu, args.sibling) if args.cpu is not None
            else physical_cpu_pairs()[0])
    topology = validate_cpu_pair(*pair)
    with benchmark_lock(args.lock, True):
        diagnostic_document, diagnostic_identity = read_json_snapshot(
            args.diagnostic_summary, "diagnostic summary")
        diagnostic = validate_diagnostic_summary(
            diagnostic_document, matrix)
        selected = selected_abba_cells(
            matrix, diagnostic, args.near_ratio, args.include_cell,
            args.exclude_cell)
        if not selected:
            raise RuntimeError("ABBA selection is empty")
        provenance, runtime = common_provenance(args, matrix)
        try:
            if diagnostic["candidate_binary_sha256"] != \
                    provenance["candidate"]["sha256"]:
                raise RuntimeError("diagnostic used another candidate binary")
            if diagnostic["baseline_binary_sha256"] != \
                    provenance["baseline"]["sha256"]:
                raise RuntimeError("diagnostic used another baseline binary")
            options = {
                "cpu_pair": list(pair), "topology": topology,
                "rounds": args.rounds, "order_per_round": list(ABBA_ORDER),
                "maximum_attempts": args.maximum_attempts,
                "idle_sibling_required": True,
                "near_ratio": args.near_ratio,
                "diagnostic_summary": diagnostic_identity,
                "include_cell": sorted(args.include_cell),
                "exclude_cell": sorted(args.exclude_cell),
                "selected_cell_ids": [cell["id"] for cell in selected],
            }
            identity = stage_identity("abba", provenance, options)
            run_dir, identity_digest = prepare_run_directory(
                args.run_dir, identity)
            started = time.monotonic()
            results = []
            resumed_count = 0
            for index, cell in enumerate(selected, 1):
                result, resumed = abba_cell(
                    cell, pair, args, identity_digest, run_dir, runtime)
                results.append(result)
                resumed_count += int(resumed)
                print("ABBA %d/%d %s%s" % (
                    index, len(selected), cell["id"],
                    " resumed" if resumed else ""), flush=True)
            source_still_frozen(provenance, runtime)
            require_json_snapshot_unchanged(
                diagnostic_identity, "diagnostic summary")
            rows = sorted(results, key=lambda value: value["cell"]["id"])
            summary = {
                "schema": SCHEMA, "status": "complete", "stage": "abba",
                "authoritative": True, "identity_sha256": identity_digest,
                "matrix_sha256": matrix_digest(matrix),
                "cell_count": len(rows), "resumed_cell_count": resumed_count,
                "accepted_child_invocations": len(rows) * args.rounds * 4,
                "elapsed_seconds": time.monotonic() - started,
                "candidate_binary_sha256": provenance["candidate"]["sha256"],
                "baseline_binary_sha256": provenance["baseline"]["sha256"],
                "rows": [{
                    "id": row["cell"]["id"], "cell": row["cell"],
                    "metrics": row["metrics"], "cpu_pair": row["cpu_pair"],
                    "selected_decode_path": row["selected_decode_path"],
                    "selected_decode_rule": row["selected_decode_rule"],
                    "workload_digests": row["workload_digests"],
                    "benchmark_byte_identity": {
                        kind: expected_benchmark_byte_identity(
                            kind, row["cell"]["bytes"],
                            matrix_baseline_contract(matrix))
                        for kind in ("baseline", "candidate")
                    },
                    "missing_original_indices":
                        row["missing_original_indices"],
                } for row in rows],
            }
            atomic_json(run_dir / "summary.json", summary)
            print(json.dumps({
                "stage": "abba", "run_dir": str(run_dir),
                "cells": summary["cell_count"],
                "elapsed_seconds": summary["elapsed_seconds"],
                "summary_sha256": sha256_file(run_dir / "summary.json"),
            }, sort_keys=True))
        finally:
            close_runtime_snapshots(runtime)


def make_dry_run_manifest(matrix, diagnostic_summary_path=None, near_ratio=1.10):
    mandatory_tags = {"loss_sweep", "reuse_sweep", "batch_sweep",
                      "large_tiling", "selector_isolation",
                      "r1_selector_fallback", "r1_selector_promoted",
                      "r1_large_recovery", "direct_odd_arity",
                      "r1_tiny_exact_main",
                      "k4_r2_production_avx2_gate"}
    mandatory = [cell for cell in matrix["cells"]
                 if set(cell["tags"]) & mandatory_tags]
    selected = None
    if diagnostic_summary_path:
        diagnostic = load_diagnostic_summary(diagnostic_summary_path, matrix)
        selected = selected_abba_cells(
            matrix, diagnostic, near_ratio, [], [])
    pairs = physical_cpu_pairs()
    ordered_cells = sorted(
        matrix["cells"],
        key=lambda cell: (cell["estimated_peak_bytes"], cell["id"]))
    assignments = [{"cell_id": cell["id"],
                    "cpu_pair": list(pairs[index % len(pairs)])}
                   for index, cell in enumerate(ordered_cells)]
    return {
        "schema": SCHEMA, "kind": "dry_run_manifest",
        "matrix_sha256": matrix_digest(matrix),
        "matrix_cell_count": len(matrix["cells"]),
        "diagnostic_child_invocations": len(matrix["cells"]) * 2,
        "detected_physical_cpu_pairs": [list(pair) for pair in pairs],
        "detected_pair_count": len(pairs),
        "deterministic_diagnostic_assignments": assignments,
        "default_memory_budget_bytes": min(mem_total_bytes() // 4, 32 << 30),
        "mandatory_abba_cell_count": len(mandatory),
        "mandatory_abba_child_invocations_three_rounds": len(mandatory) * 12,
        "selected_abba_cell_count": None if selected is None else len(selected),
        "selected_abba_child_invocations_three_rounds": (
            None if selected is None else len(selected) * 12),
        "maximum_abba_cell_count": len(matrix["cells"]),
        "maximum_abba_child_invocations_three_rounds":
            len(matrix["cells"]) * 12,
        "sizes": matrix["sizes"], "transform_sides": matrix["transform_sides"],
        "losses": matrix["losses"], "reuse_counts": matrix["reuse_counts"],
        "batch_counts": matrix["batch_counts"],
        "selector_isolation_cells": [
            cell["id"] for cell in matrix["cells"]
            if "selector_isolation" in cell["tags"]],
    }


def run_merge(args):
    diagnostic, diagnostic_identity = read_json_snapshot(
        args.diagnostic_summary, "diagnostic merge summary")
    abba, abba_identity = read_json_snapshot(
        args.abba_summary, "ABBA merge summary")
    if (diagnostic.get("schema") != SCHEMA or abba.get("schema") != SCHEMA or
            diagnostic.get("status") != "complete" or
            abba.get("status") != "complete" or
            diagnostic.get("stage") != "diagnostic" or
            abba.get("stage") != "abba" or
            diagnostic.get("authoritative") is not False or
            abba.get("authoritative") is not True):
        raise RuntimeError("merge inputs have wrong stages")
    if (diagnostic.get("cell_count") != len(diagnostic.get("rows", [])) or
            abba.get("cell_count") != len(abba.get("rows", []))):
        raise RuntimeError("merge input row count mismatch")
    if diagnostic.get("matrix_sha256") != abba.get("matrix_sha256"):
        raise RuntimeError("merge input matrix mismatch")
    if diagnostic.get("baseline_binary_sha256") != \
            abba.get("baseline_binary_sha256") or \
            diagnostic.get("candidate_binary_sha256") != \
            abba.get("candidate_binary_sha256"):
        raise RuntimeError("merge input binary mismatch")
    bundle = {
        "schema": SCHEMA, "kind": "merged_evidence_bundle",
        "diagnostic_authoritative": False,
        "abba_authoritative": True,
        "diagnostic_summary": diagnostic,
        "abba_summary": abba,
        "input_identity": {
            "diagnostic": diagnostic_identity,
            "abba": abba_identity,
        },
    }
    require_json_snapshot_unchanged(
        diagnostic_identity, "diagnostic merge summary")
    require_json_snapshot_unchanged(abba_identity, "ABBA merge summary")
    atomic_json(args.output, bundle)
    print(json.dumps({"output": os.path.realpath(args.output),
                      "sha256": sha256_file(args.output)}, sort_keys=True))


def self_test():
    check_count = 0

    def check(condition, label):
        nonlocal check_count
        check_count += 1
        if not condition:
            raise RuntimeError("expanded final-map self-test failed: " + label)

    def expect_runtime_error(action, message, label):
        try:
            action()
        except RuntimeError as error:
            check(message in str(error), label + " error contract")
        else:
            check(False, label + " rejection")

    def assert_line_numbers(source):
        return tuple(
            node.lineno for node in ast.walk(ast.parse(source))
            if isinstance(node, ast.Assert)
        )

    matrix = make_matrix()
    validate_matrix(matrix)
    canonical_digest = matrix_digest(matrix)
    check("profile" not in matrix, "canonical matrix has no new profile field")
    check(
        len(matrix["cells"]) == CANONICAL_MATRIX_CELL_COUNT,
        "canonical matrix retained cell count",
    )
    check(
        canonical_digest == CANONICAL_MATRIX_SHA256,
        "canonical matrix retained frozen digest",
    )
    check(
        matrix_baseline_sha256(matrix) == BASELINE_SHA256 and
        matrix_baseline_contract(matrix) == BASELINE_CONTRACT_LEGACY_V1,
        "canonical matrix retains exact-main v1 baseline",
    )
    r1_matrix = make_matrix(R1_TINY_EXACT_MAIN_PROFILE)
    validate_matrix(r1_matrix)
    check(
        r1_matrix == make_matrix(R1_TINY_EXACT_MAIN_PROFILE),
        "R=1 tiny exact-main profile determinism",
    )
    check(
        matrix_digest(r1_matrix) == R1_TINY_EXACT_MAIN_MATRIX_SHA256,
        "R=1 tiny exact-main retained frozen digest",
    )
    check(
        matrix_baseline_sha256(r1_matrix) == BASELINE_SHA256 and
        matrix_baseline_contract(r1_matrix) ==
            BASELINE_CONTRACT_LEGACY_V1,
        "R=1 tiny exact-main retains exact-main v1 baseline",
    )
    r1_ids = [cell["id"] for cell in r1_matrix["cells"]]
    r1_seeds = [cell["seed"] for cell in r1_matrix["cells"]]
    check(
        len(r1_ids) == len(set(r1_ids)) and
        len(r1_seeds) == len(set(r1_seeds)),
        "R=1 tiny exact-main unique IDs and seeds",
    )
    observed_r1_grid = {
        (cell["K"], cell["bytes"],
         next(tag[len("missing_"):] for tag in cell["tags"]
              if tag.startswith("missing_")),
         cell["expected_missing_original_index"])
        for cell in r1_matrix["cells"]
    }
    expected_r1_grid = {
        (k, byte_count, position,
         r1_missing_position_index(k, position))
        for k in R1_TINY_EXACT_MAIN_K
        for byte_count in R1_TINY_EXACT_MAIN_BYTES
        for position in R1_TINY_MISSING_POSITIONS
    }
    check(
        observed_r1_grid == expected_r1_grid,
        "R=1 tiny exact-main K/bytes/loss-position coverage",
    )
    check(
        all(r1_missing_index(cell["K"], cell["seed"]) ==
            cell["expected_missing_original_index"]
            for cell in r1_matrix["cells"]),
        "R=1 tiny exact-main seeds select declared losses",
    )
    check(
        all(cell["reuse"] == R1_TINY_EXACT_MAIN_REUSE and
            cell["loss"] == 1 and cell["batch"] == 1
            for cell in r1_matrix["cells"]),
        "R=1 tiny exact-main bounded call shape",
    )
    k4_r2_matrix = make_matrix(K4_R2_PRODUCTION_AVX2_PROFILE)
    validate_matrix(k4_r2_matrix)
    check(
        k4_r2_matrix == make_matrix(K4_R2_PRODUCTION_AVX2_PROFILE),
        "K4/R2 production AVX2 profile determinism",
    )
    check(
        matrix_digest(k4_r2_matrix) ==
        K4_R2_PRODUCTION_AVX2_MATRIX_SHA256,
        "K4/R2 production AVX2 frozen digest",
    )
    check(
        matrix_baseline_sha256(k4_r2_matrix) ==
        PADDED_MAIN_BASELINE_SHA256 and
        matrix_baseline_contract(k4_r2_matrix) ==
            BASELINE_CONTRACT_LOGICAL_ZERO_PAD_V2 and
        matrix_baseline_sha256(k4_r2_matrix) != BASELINE_SHA256,
        "K4/R2 production AVX2 pins exact-main v2 baseline",
    )
    k4_r2_ids = [cell["id"] for cell in k4_r2_matrix["cells"]]
    k4_r2_seeds = [cell["seed"] for cell in k4_r2_matrix["cells"]]
    check(
        len(k4_r2_ids) == len(set(k4_r2_ids)) and
        len(k4_r2_seeds) == len(set(k4_r2_seeds)),
        "K4/R2 production AVX2 unique IDs and seeds",
    )
    observed_k4_r2_targets = {
        (cell["bytes"],
         tuple(cell["expected_missing_original_indices"]),
         cell["loss_pattern"])
        for cell in k4_r2_matrix["cells"]
        if "k4_r2_target" in cell["tags"]
    }
    check(
        observed_k4_r2_targets == set(K4_R2_TARGET_LOSS_PATTERNS),
        "K4/R2 production AVX2 target size/loss-pattern grid",
    )
    observed_k4_r2_neighbors = {
        (cell["K"], cell["R"],
         tuple(cell["expected_missing_original_indices"]),
         cell["loss_pattern"])
        for cell in k4_r2_matrix["cells"]
        if "k4_r2_neighbor" in cell["tags"]
    }
    check(
        observed_k4_r2_neighbors == set(K4_R2_NEIGHBOR_LOSS_PATTERNS),
        "K4/R2 production AVX2 neighbor loss-pattern grid",
    )
    check(
        all(selected_missing_indices(
                cell["K"], cell["loss"], cell["seed"]) ==
            tuple(cell["expected_missing_original_indices"])
            for cell in k4_r2_matrix["cells"]),
        "K4/R2 production AVX2 seeds select declared losses",
    )
    check(
        all(cell["reuse"] == K4_R2_PRODUCTION_AVX2_REUSE and
            cell["batch"] == 1 and
            cell["diagnostic_iterations"] ==
                K4_R2_PRODUCTION_AVX2_ITERATIONS and
            cell["diagnostic_warmup"] ==
                K4_R2_PRODUCTION_AVX2_WARMUP and
            cell["abba_iterations"] ==
                K4_R2_PRODUCTION_AVX2_ITERATIONS and
            cell["abba_warmup"] == K4_R2_PRODUCTION_AVX2_WARMUP
            for cell in k4_r2_matrix["cells"]),
        "K4/R2 production AVX2 timing and reuse policy",
    )

    def byte_identity_document(kind, logical_bytes, baseline_contract):
        physical_bytes = physical_shard_bytes(
            kind, logical_bytes, baseline_contract)
        document = {
            "schema": "candidate-self-test",
            "parameters": {"shard_bytes": physical_bytes},
            "resolved": {},
            "correctness": {},
        }
        if (kind == "baseline" and baseline_contract ==
                BASELINE_CONTRACT_LOGICAL_ZERO_PAD_V2):
            padded = physical_bytes != logical_bytes
            document["schema"] = (
                "leopard-main-benchmark-v2" if padded else
                "leopard-main-benchmark-v1")
            document["parameters"]["logical_shard_bytes"] = logical_bytes
            document["resolved"].update({
                "padded_application_bytes": padded,
                "padding_policy": "zero suffix per shard",
            })
            document["correctness"][
                "logical_prefix_fingerprinted"] = True
        elif kind == "baseline":
            document["schema"] = "leopard-main-benchmark-v1"
        return document

    for logical_bytes, baseline_bytes in (
            (63, 64), (64, 64), (65, 128), (2049, 2112)):
        baseline_arguments = benchmark_byte_arguments(
            "baseline", logical_bytes,
            BASELINE_CONTRACT_LOGICAL_ZERO_PAD_V2)
        candidate_arguments = benchmark_byte_arguments(
            "candidate", logical_bytes,
            BASELINE_CONTRACT_LOGICAL_ZERO_PAD_V2)
        check(
            baseline_arguments[:2] == ["--bytes", str(baseline_bytes)] and
            candidate_arguments == ["--bytes", str(logical_bytes)],
            "exact-main tail physical-byte command %d" % logical_bytes,
        )
        if baseline_bytes == logical_bytes:
            check(
                "--logical-bytes" not in baseline_arguments,
                "aligned exact-main command remains unchanged",
            )
        else:
            logical_position = baseline_arguments.index("--logical-bytes")
            check(
                baseline_arguments[logical_position + 1] ==
                str(logical_bytes),
                "exact-main tail command declares logical bytes %d" %
                logical_bytes,
            )
        baseline_document = byte_identity_document(
            "baseline", logical_bytes,
            BASELINE_CONTRACT_LOGICAL_ZERO_PAD_V2)
        candidate_document = byte_identity_document(
            "candidate", logical_bytes,
            BASELINE_CONTRACT_LOGICAL_ZERO_PAD_V2)
        cell = {"bytes": logical_bytes}
        baseline_identity = validate_benchmark_byte_identity(
            "baseline", baseline_document, cell,
            BASELINE_CONTRACT_LOGICAL_ZERO_PAD_V2)
        candidate_identity = validate_benchmark_byte_identity(
            "candidate", candidate_document, cell,
            BASELINE_CONTRACT_LOGICAL_ZERO_PAD_V2)
        check(
            baseline_identity["physical_shard_bytes"] == baseline_bytes and
            baseline_identity["logical_shard_bytes"] == logical_bytes and
            baseline_identity["padded_application_bytes"] is
                (baseline_bytes != logical_bytes) and
            baseline_identity["ratio_scope"] == "time_per_call" and
            baseline_identity["digest_scope"] == "logical_prefix" and
            candidate_identity["physical_shard_bytes"] == logical_bytes and
            candidate_identity["logical_shard_bytes"] == logical_bytes and
            candidate_identity["padded_application_bytes"] is False,
            "exact-main/candidate byte identity %d" % logical_bytes,
        )
        wrong_physical = json.loads(json.dumps(baseline_document))
        wrong_physical["parameters"]["shard_bytes"] = logical_bytes
        if logical_bytes == baseline_bytes:
            wrong_physical["parameters"]["shard_bytes"] += 64
        expect_runtime_error(
            lambda value=wrong_physical, selected_cell=cell:
                validate_benchmark_byte_identity(
                    "baseline", value, selected_cell,
                    BASELINE_CONTRACT_LOGICAL_ZERO_PAD_V2),
            "physical shard bytes mismatch",
            "exact-main rejects wrong physical bytes %d" % logical_bytes,
        )
        wrong_logical = json.loads(json.dumps(baseline_document))
        wrong_logical["parameters"]["logical_shard_bytes"] += 1
        expect_runtime_error(
            lambda value=wrong_logical, selected_cell=cell:
                validate_benchmark_byte_identity(
                    "baseline", value, selected_cell,
                    BASELINE_CONTRACT_LOGICAL_ZERO_PAD_V2),
            "logical-byte padding identity mismatch",
            "exact-main rejects wrong logical bytes %d" % logical_bytes,
        )
        wrong_padding = json.loads(json.dumps(baseline_document))
        wrong_padding["resolved"]["padded_application_bytes"] = not (
            baseline_bytes != logical_bytes)
        expect_runtime_error(
            lambda value=wrong_padding, selected_cell=cell:
                validate_benchmark_byte_identity(
                    "baseline", value, selected_cell,
                    BASELINE_CONTRACT_LOGICAL_ZERO_PAD_V2),
            "logical-byte padding identity mismatch",
            "exact-main rejects wrong padding flag %d" % logical_bytes,
        )
        wrong_policy = json.loads(json.dumps(baseline_document))
        wrong_policy["resolved"]["padding_policy"] = "unspecified"
        expect_runtime_error(
            lambda value=wrong_policy, selected_cell=cell:
                validate_benchmark_byte_identity(
                    "baseline", value, selected_cell,
                    BASELINE_CONTRACT_LOGICAL_ZERO_PAD_V2),
            "logical-byte padding identity mismatch",
            "exact-main rejects wrong padding policy %d" % logical_bytes,
        )
        wrong_candidate_physical = json.loads(json.dumps(candidate_document))
        wrong_candidate_physical["parameters"]["shard_bytes"] += 64
        expect_runtime_error(
            lambda value=wrong_candidate_physical, selected_cell=cell:
                validate_benchmark_byte_identity(
                    "candidate", value, selected_cell,
                    BASELINE_CONTRACT_LOGICAL_ZERO_PAD_V2),
            "physical shard bytes mismatch",
            "candidate rejects padded physical bytes %d" % logical_bytes,
        )
        adapted_candidate = json.loads(json.dumps(candidate_document))
        adapted_candidate["parameters"][
            "logical_shard_bytes"] = logical_bytes
        expect_runtime_error(
            lambda value=adapted_candidate, selected_cell=cell:
                validate_benchmark_byte_identity(
                    "candidate", value, selected_cell,
                    BASELINE_CONTRACT_LOGICAL_ZERO_PAD_V2),
            "unexpectedly reported logical-byte adaptation",
            "candidate rejects logical-byte adaptation %d" % logical_bytes,
        )
    legacy_v1_cell = {"bytes": 64}
    legacy_v1_document = byte_identity_document(
        "baseline", 64, BASELINE_CONTRACT_LEGACY_V1)
    legacy_v1_document["correctness"]["round_trip"] = True
    check(
        benchmark_byte_arguments(
            "baseline", 64, BASELINE_CONTRACT_LEGACY_V1) ==
        ["--bytes", "64"] and
        validate_benchmark_byte_identity(
            "baseline", legacy_v1_document, legacy_v1_cell,
            BASELINE_CONTRACT_LEGACY_V1)["baseline_contract"] ==
        BASELINE_CONTRACT_LEGACY_V1,
        "sparse exact-main v1 byte identity remains accepted",
    )
    expect_runtime_error(
        lambda: validate_benchmark_byte_identity(
            "baseline", legacy_v1_document, legacy_v1_cell,
            BASELINE_CONTRACT_LOGICAL_ZERO_PAD_V2),
        "logical-byte padding identity mismatch",
        "K4 v2 contract rejects sparse exact-main v1 output",
    )
    aligned_v2_document = byte_identity_document(
        "baseline", 64, BASELINE_CONTRACT_LOGICAL_ZERO_PAD_V2)
    expect_runtime_error(
        lambda: validate_benchmark_byte_identity(
            "baseline", aligned_v2_document, legacy_v1_cell,
            BASELINE_CONTRACT_LEGACY_V1),
        "legacy v1 baseline byte identity mismatch",
        "legacy contract rejects enhanced exact-main v2 output",
    )
    expect_runtime_error(
        lambda: benchmark_byte_arguments(
            "baseline", 63, BASELINE_CONTRACT_LEGACY_V1),
        "requires 64-byte-aligned shards",
        "legacy exact-main v1 tail command",
    )
    identity_cell = r1_matrix["cells"][0]

    def identity_call(kind, missing_indices):
        if isinstance(missing_indices, int):
            missing_indices = [missing_indices]
        return {
            "kind": kind,
            "document": {
                "parameters": {
                    "missing_original_indices": list(missing_indices)},
                "resolved": {
                    "selected_decode_path": "direct_xor",
                    "selected_decode_rule": "r1_self_test"},
                "workload_digests": {
                    "algorithm": "sha256", "original_data": "a",
                    "transmitted_parity": "b",
                    "recovered_originals": "c"},
            },
        }

    required_missing = identity_cell["expected_missing_original_index"]
    identity_calls = [
        identity_call("baseline", required_missing),
        identity_call("candidate", required_missing),
    ]
    check(
        validate_call_identity(identity_calls, identity_cell)[
            "missing_original_indices"] == [required_missing],
        "R=1 exact expected loss accepts every matching invocation",
    )
    for invocation_index, kind in enumerate(("baseline", "candidate")):
        wrong_calls = json.loads(json.dumps(identity_calls))
        wrong_calls[invocation_index]["document"]["parameters"][
            "missing_original_indices"] = [
                (required_missing + 1) % identity_cell["K"]]
        expect_runtime_error(
            lambda calls=wrong_calls: validate_call_identity(
                calls, identity_cell),
            "expected missing original index mismatch",
            "R=1 exact rejects wrong %s loss" % kind,
        )
    exact_identity_cell = next(
        cell for cell in k4_r2_matrix["cells"]
        if cell["K"] == 4 and cell["R"] == 2 and cell["bytes"] == 64)
    exact_missing = exact_identity_cell[
        "expected_missing_original_indices"]
    exact_identity_calls = [
        identity_call("baseline", exact_missing),
        identity_call("candidate", exact_missing),
    ]
    check(
        validate_call_identity(exact_identity_calls, exact_identity_cell)[
            "missing_original_indices"] == exact_missing,
        "K4/R2 exact expected loss set accepts matching invocations",
    )
    for invocation_index, kind in enumerate(("baseline", "candidate")):
        wrong_calls = json.loads(json.dumps(exact_identity_calls))
        wrong_calls[invocation_index]["document"]["parameters"][
            "missing_original_indices"] = [0, 1]
        expect_runtime_error(
            lambda calls=wrong_calls: validate_call_identity(
                calls, exact_identity_cell),
            "expected missing original set mismatch",
            "K4/R2 exact rejects wrong %s loss set" % kind,
        )
    digest_mismatch_calls = json.loads(json.dumps(exact_identity_calls))
    digest_mismatch_calls[1]["document"]["workload_digests"][
        "transmitted_parity"] = "different"
    expect_runtime_error(
        lambda: validate_call_identity(
            digest_mismatch_calls, exact_identity_cell),
        "wire/recovery digest or missing-set mismatch",
        "K4/R2 exact rejects baseline/candidate digest mismatch",
    )
    check(
        parse_cpu_list("0-2,5,8-9") == {0, 1, 2, 5, 8, 9},
        "CPU-list parsing",
    )
    selectors = {cell["id"] for cell in matrix["cells"]
                 if "selector_isolation" in cell["tags"]}
    expected_selectors = {
        cell_key(k, 32, size, 32, 8, 1)
        for k in (33, 34, 35) for size in (16384, 65536)}
    check(selectors == expected_selectors, "selector-isolation grid")
    r1_boundaries = {
        (cell["K"], cell["bytes"], tuple(cell["tags"]))
        for cell in matrix["cells"]
        if "r1_selector_fallback" in cell["tags"] or
           "r1_selector_promoted" in cell["tags"]
    }
    expected_r1_boundaries = set()
    for lower_k, promoted_k, byte_count in R1_SELECTOR_BOUNDARIES:
        expected_r1_boundaries.add(
            (lower_k, byte_count, ("r1_selector_fallback",)))
        expected_r1_boundaries.add(
            (promoted_k, byte_count, ("r1_selector_promoted",)))
    check(
        r1_boundaries == expected_r1_boundaries,
        "R=1 selector-boundary grid",
    )
    direct_odd = {
        (cell["K"], cell["R"], cell["bytes"], cell["loss"])
        for cell in matrix["cells"]
        if "direct_odd_arity" in cell["tags"]
    }
    check(
        direct_odd == {
            (65, 65, byte_count, loss)
            for byte_count in (2048, 65536, 1048576)
            for loss in DIRECT_ODD_LOSSES
        },
        "direct odd-arity grid",
    )
    encoded = canonical_bytes(matrix)
    check(json.loads(encoded) == matrix, "canonical matrix round trip")
    check(
        matrix_digest(matrix) == hashlib.sha256(encoded).hexdigest(),
        "canonical matrix digest",
    )
    selector_mutation = json.loads(encoded)
    selector_cell = next(
        cell for cell in selector_mutation["cells"]
        if "selector_isolation" in cell["tags"])
    selector_cell["tags"].remove("selector_isolation")
    expect_runtime_error(
        lambda: validate_matrix(selector_mutation),
        "deterministic built-in matrix",
        "selector-grid mutation",
    )
    source_mutation = json.loads(encoded)
    source_mutation["baseline_commit"] = "0" * 40
    expect_runtime_error(
        lambda: validate_matrix(source_mutation),
        "deterministic built-in matrix",
        "matrix provenance mutation",
    )
    k4_r2_mutation = json.loads(canonical_bytes(k4_r2_matrix))
    k4_r2_mutation["cells"][0]["expected_missing_original_indices"][0] += 1
    expect_runtime_error(
        lambda: validate_matrix(k4_r2_mutation),
        "deterministic built-in matrix",
        "K4/R2 exact loss-pattern mutation",
    )
    k4_r2_v1_baseline_mutation = json.loads(canonical_bytes(k4_r2_matrix))
    k4_r2_v1_baseline_mutation["baseline_sha256"] = BASELINE_SHA256
    expect_runtime_error(
        lambda: matrix_baseline_sha256(k4_r2_v1_baseline_mutation),
        "differs from its profile",
        "K4/R2 rejects exact-main v1 baseline identity",
    )
    canonical_v2_baseline_mutation = json.loads(encoded)
    canonical_v2_baseline_mutation[
        "baseline_sha256"] = PADDED_MAIN_BASELINE_SHA256
    expect_runtime_error(
        lambda: matrix_baseline_sha256(canonical_v2_baseline_mutation),
        "differs from its profile",
        "canonical profile rejects exact-main v2 baseline identity",
    )
    malformed_baseline_mutation = json.loads(encoded)
    malformed_baseline_mutation["baseline_sha256"] = "invalid"
    expect_runtime_error(
        lambda: matrix_baseline_sha256(malformed_baseline_mutation),
        "is invalid",
        "malformed matrix baseline identity",
    )

    calls = []
    for kind, value in zip(ABBA_ORDER * 3,
                           (10, 5, 5, 10) * 3):
        calls.append({"kind": kind, "times_us": {
            "encode": value, "decode_first": value,
            "decode_reuse": value}})
    metrics = abba_metrics(calls, 3)
    for metric, value in sorted(metrics.items()):
        check(
            abs(value["ratio"] - 2.0) < 1e-12,
            metric + " ABBA ratio",
        )
        check(
            value["ci95_low"] <= value["ratio"] <= value["ci95_high"],
            metric + " ABBA confidence interval",
        )
    point_calls = [
        {"kind": "baseline", "times_us": {
            "encode": 12.0, "decode_first": 18.0, "decode_reuse": 24.0}},
        {"kind": "candidate", "times_us": {
            "encode": 6.0, "decode_first": 9.0, "decode_reuse": 12.0}},
    ]
    check(
        point_ratios(point_calls) == {
            "encode": 2.0, "decode_first": 2.0, "decode_reuse": 2.0},
        "diagnostic point ratios",
    )
    mutated_calls = json.loads(json.dumps(calls))
    mutated_calls[0]["times_us"]["encode"] = 20
    check(
        abba_metrics(mutated_calls, 3)["encode"]["ratio"] !=
        metrics["encode"]["ratio"],
        "ratio mutation changes evidence",
    )
    expect_runtime_error(
        lambda: abba_metrics(calls[:-1], 3),
        "ABBA call count mismatch",
        "truncated ABBA evidence",
    )
    expect_runtime_error(
        lambda: point_ratios(point_calls + [point_calls[0]]),
        "one baseline and one candidate",
        "duplicate point-ratio call",
    )
    nonfinite_point_calls = json.loads(json.dumps(point_calls))
    nonfinite_point_calls[1]["times_us"]["encode"] = float("nan")
    expect_runtime_error(
        lambda: point_ratios(nonfinite_point_calls),
        "positive and finite",
        "non-finite point-ratio timing",
    )
    expect_runtime_error(
        lambda: confidence([0.0, float("nan")]),
        "finite real numbers",
        "non-finite confidence input",
    )

    diagnostic_fixture = {
        "schema": SCHEMA,
        "stage": "diagnostic",
        "status": "complete",
        "authoritative": False,
        "matrix_sha256": matrix_digest(matrix),
        "cell_count": len(matrix["cells"]),
        "baseline_binary_sha256": "1" * 64,
        "candidate_binary_sha256": "2" * 64,
        "rows": [
            {
                "id": cell["id"],
                "cell": cell,
                "ratios": {
                    "encode": 2.0,
                    "decode_first": 2.0,
                    "decode_reuse": 2.0,
                },
            }
            for cell in matrix["cells"]
        ],
    }
    check(
        validate_diagnostic_summary(diagnostic_fixture, matrix) is
        diagnostic_fixture,
        "complete diagnostic summary",
    )
    r1_diagnostic_fixture = json.loads(json.dumps(diagnostic_fixture))
    r1_diagnostic_fixture["matrix_sha256"] = matrix_digest(r1_matrix)
    r1_diagnostic_fixture["cell_count"] = len(r1_matrix["cells"])
    r1_diagnostic_fixture["rows"] = [
        {
            "id": cell["id"], "cell": cell,
            "ratios": {metric: 2.0 for metric in RATIO_METRICS},
        }
        for cell in r1_matrix["cells"]
    ]
    validate_diagnostic_summary(r1_diagnostic_fixture, r1_matrix)
    check(
        {cell["id"] for cell in selected_abba_cells(
            r1_matrix, r1_diagnostic_fixture, 1.10, [], [])} ==
        {cell["id"] for cell in r1_matrix["cells"]},
        "R=1 tiny exact-main profile makes every cell ABBA-mandatory",
    )
    k4_r2_diagnostic_fixture = json.loads(json.dumps(diagnostic_fixture))
    k4_r2_diagnostic_fixture["matrix_sha256"] = matrix_digest(k4_r2_matrix)
    k4_r2_diagnostic_fixture["cell_count"] = len(k4_r2_matrix["cells"])
    k4_r2_diagnostic_fixture["rows"] = [
        {
            "id": cell["id"], "cell": cell,
            "ratios": {metric: 2.0 for metric in RATIO_METRICS},
        }
        for cell in k4_r2_matrix["cells"]
    ]
    validate_diagnostic_summary(k4_r2_diagnostic_fixture, k4_r2_matrix)
    check(
        {cell["id"] for cell in selected_abba_cells(
            k4_r2_matrix, k4_r2_diagnostic_fixture, 1.10, [], [])} ==
        {cell["id"] for cell in k4_r2_matrix["cells"]},
        "K4/R2 production AVX2 profile makes every cell ABBA-mandatory",
    )
    check(
        make_dry_run_manifest(k4_r2_matrix)[
            "mandatory_abba_cell_count"] == len(k4_r2_matrix["cells"]),
        "K4/R2 dry run declares every cell ABBA-mandatory",
    )
    target_cell = next(
        cell for cell in matrix["cells"] if set(cell["tags"]) == {"core"})
    selected_before = selected_abba_cells(
        matrix, diagnostic_fixture, 1.10, [], [])
    check(
        target_cell["id"] not in {
            cell["id"] for cell in selected_before},
        "fast core-only ratio stays out of ABBA confirmation",
    )
    ratio_selection_fixture = json.loads(json.dumps(diagnostic_fixture))
    ratio_selection_row = next(
        row for row in ratio_selection_fixture["rows"]
        if row["id"] == target_cell["id"])
    ratio_selection_row["ratios"]["encode"] = 1.05
    validate_diagnostic_summary(ratio_selection_fixture, matrix)
    selected_after = selected_abba_cells(
        matrix, ratio_selection_fixture, 1.10, [], [])
    check(
        target_cell["id"] in {cell["id"] for cell in selected_after},
        "near diagnostic ratio enters ABBA confirmation",
    )
    ratio_mutation = json.loads(json.dumps(diagnostic_fixture))
    ratio_mutation["rows"][0]["ratios"]["encode"] = float("nan")
    expect_runtime_error(
        lambda: validate_diagnostic_summary(ratio_mutation, matrix),
        "positive and finite",
        "non-finite diagnostic ratio",
    )
    ratio_shape_mutation = json.loads(json.dumps(diagnostic_fixture))
    del ratio_shape_mutation["rows"][0]["ratios"]["decode_reuse"]
    expect_runtime_error(
        lambda: validate_diagnostic_summary(ratio_shape_mutation, matrix),
        "wrong ratio metrics",
        "missing diagnostic ratio",
    )
    diagnostic_cell_mutation = json.loads(json.dumps(diagnostic_fixture))
    diagnostic_cell_mutation["rows"][0]["cell"]["K"] += 1
    expect_runtime_error(
        lambda: validate_diagnostic_summary(
            diagnostic_cell_mutation, matrix),
        "cell identity mismatch",
        "diagnostic selector-grid mutation",
    )
    diagnostic_provenance_mutation = json.loads(
        json.dumps(diagnostic_fixture))
    diagnostic_provenance_mutation["candidate_binary_sha256"] = "invalid"
    expect_runtime_error(
        lambda: validate_diagnostic_summary(
            diagnostic_provenance_mutation, matrix),
        "invalid candidate_binary_sha256",
        "diagnostic binary provenance mutation",
    )

    synthetic_provenance = {
        "candidate": {"sha256": "a" * 64},
        "baseline": {"sha256": "b" * 64},
        "candidate_source": {"head": "c" * 40, "tree": "d" * 40},
    }
    identity = stage_identity(
        "diagnostic", synthetic_provenance, {"maximum_attempts": 3})
    identity_digest = sha256_bytes(canonical_bytes(identity))
    provenance_mutation = json.loads(json.dumps(identity))
    provenance_mutation["provenance"]["candidate"]["sha256"] = "e" * 64
    check(
        sha256_bytes(canonical_bytes(provenance_mutation)) != identity_digest,
        "stage identity binds candidate provenance",
    )
    with tempfile.TemporaryDirectory(
            prefix="leopard2-expanded-identity-") as directory:
        run_dir, observed_digest = prepare_run_directory(directory, identity)
        check(observed_digest == identity_digest, "run identity digest")
        resumed_dir, resumed_digest = prepare_run_directory(
            directory, identity)
        check(
            resumed_dir == run_dir and resumed_digest == identity_digest,
            "identical run identity resumes",
        )
        expect_runtime_error(
            lambda: prepare_run_directory(directory, provenance_mutation),
            "run directory identity differs",
            "provenance mutation",
        )

    with tempfile.TemporaryDirectory(prefix="leopard2-expanded-snapshot-") \
            as directory:
        executable = Path(directory) / "fixture"
        executable.write_bytes(Path("/bin/true").read_bytes())
        executable.chmod(0o500)
        build_executable = Path(directory) / "build-fixture"
        build_executable.write_bytes(executable.read_bytes())
        build_executable.chmod(0o500)
        split_equivalence = equivalent_executable_identities(
            build_executable, executable)
        check(
            split_equivalence["byte_for_byte_equal"] is True and
            split_equivalence["build_tree_executable"]["path"] !=
            split_equivalence["selected_timing_artifact"]["path"],
            "separate build executable and timing artifact equivalence",
        )
        check(
            equivalent_executable_identities(
                executable, executable)["equal_sha256"] ==
            split_equivalence["equal_sha256"],
            "candidate build executable defaults to selected artifact",
        )
        mismatched_executable = Path(directory) / "mismatched-fixture"
        mismatched_executable.write_bytes(executable.read_bytes() + b"X")
        mismatched_executable.chmod(0o500)
        expect_runtime_error(
            lambda: equivalent_executable_identities(
                build_executable, mismatched_executable),
            "timing artifact differ",
            "mismatched build executable and timing artifact",
        )
        source_identity, descriptor, snapshot_identity = snapshot_executable(
            executable, "self-test fixture")
        try:
            check(
                source_identity["sha256"] == snapshot_identity["sha256"],
                "snapshot content identity",
            )
            executable.chmod(0o700)
            executable.write_bytes(b"replaced after immutable snapshot\n")
            executable.chmod(0o500)
            check(
                file_identity(
                    executable, "mutated self-test fixture") !=
                source_identity,
                "source artifact replacement changes identity",
            )
            completed = subprocess.run(
                [snapshot_path(descriptor)], stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                pass_fds=(descriptor,), timeout=5, check=False)
            check(
                completed.returncode == 0,
                "sealed snapshot survives source replacement",
            )
            try:
                os.pwrite(descriptor, b"X", 0)
            except OSError:
                sealed_write_rejected = True
            else:
                sealed_write_rejected = False
            check(sealed_write_rejected, "sealed snapshot rejects writes")
            check(
                sealed_snapshot_identity(
                    descriptor, "self-test fixture") == snapshot_identity,
                "sealed snapshot identity remains stable",
            )
        finally:
            os.close(descriptor)
        expect_runtime_error(
            lambda: snapshot_executable(
                executable, "wrong-hash fixture", "0" * 64),
            "SHA-256 mismatch",
            "artifact hash mutation",
        )
        unsealed_descriptor = linux_memfd_create(
            "leopard2-expanded-unsealed-self-test")
        try:
            os.write(unsealed_descriptor, Path("/bin/true").read_bytes())
            expect_runtime_error(
                lambda: sealed_snapshot_identity(
                    unsealed_descriptor, "unsealed self-test fixture"),
                "not immutably sealed",
                "unsealed artifact",
            )
        finally:
            os.close(unsealed_descriptor)

    with tempfile.TemporaryDirectory(prefix="leopard2-expanded-json-") \
            as directory:
        summary_path = Path(directory) / "summary.json"
        summary_path.write_text('{"selected":"old"}\n', encoding="utf-8")
        document, identity = read_json_snapshot(
            summary_path, "self-test diagnostic summary")
        check(document == {"selected": "old"}, "snapshot JSON document")
        check(
            identity["sha256"] == sha256_bytes(b'{"selected":"old"}\n'),
            "snapshot JSON byte identity",
        )
        replacement = Path(directory) / "replacement.json"
        replacement.write_text('{"selected":"new"}\n', encoding="utf-8")
        os.replace(replacement, summary_path)
        check(
            document == {"selected": "old"},
            "parsed JSON is isolated from pathname replacement",
        )
        expect_runtime_error(
            lambda: require_json_snapshot_unchanged(
                identity, "self-test diagnostic summary"),
            "changed after",
            "diagnostic pathname replacement",
        )

    with tempfile.TemporaryDirectory(prefix="leopard2-expanded-source-") \
            as directory:
        root = Path(directory)
        checked_output(["git", "init", "-q"], cwd=root)
        tracked = root / "tracked.txt"
        tracked.write_text("tracked\n", encoding="utf-8")
        checked_output(["git", "add", "tracked.txt"], cwd=root)
        checked_output([
            "git", "-c", "user.name=Leopard2 Self Test", "-c",
            "user.email=leopard2-self-test.invalid", "commit", "-qm",
            "fixture"], cwd=root)
        commit = checked_output(["git", "rev-parse", "HEAD"], cwd=root)
        clean_provenance = source_provenance(root, commit)
        check(clean_provenance["status"] == "clean", "clean source status")
        check(clean_provenance["head"] == commit, "clean source commit")
        check(
            clean_provenance["tree"] ==
            checked_output(["/usr/bin/git", "rev-parse", "HEAD^{tree}"],
                           cwd=root),
            "clean source tree",
        )
        rogue = root / "untracked.txt"
        rogue.write_text("untracked\n", encoding="utf-8")
        expect_runtime_error(
            lambda: source_provenance(root, commit),
            "tracked changes",
            "untracked source mutation",
        )
        rogue.unlink()

        tracked.write_text("dirty\n", encoding="utf-8")
        expect_runtime_error(
            lambda: source_provenance(root, commit),
            "tracked changes",
            "tracked source mutation",
        )
        checked_output(
            ["/usr/bin/git", "checkout", "--", "tracked.txt"], cwd=root)

        checked_output(
            ["git", "update-index", "--assume-unchanged", "tracked.txt"],
            cwd=root)
        expect_runtime_error(
            lambda: source_provenance(root, commit),
            "non-default flag",
            "assume-unchanged index flag",
        )
        checked_output(
            ["git", "update-index", "--no-assume-unchanged", "tracked.txt"],
            cwd=root)
        checked_output(
            ["git", "update-index", "--skip-worktree", "tracked.txt"],
            cwd=root)
        expect_runtime_error(
            lambda: source_provenance(root, commit),
            "non-default flag",
            "skip-worktree index flag",
        )
        checked_output(
            ["git", "update-index", "--no-skip-worktree", "tracked.txt"],
            cwd=root)

        with tempfile.TemporaryDirectory(
                prefix="leopard2-expanded-fake-git-") as fake_directory:
            fake_bin = Path(fake_directory)
            fake_git = fake_bin / "git"
            fake_git.write_text("#!/bin/sh\nexit 97\n", encoding="utf-8")
            fake_git.chmod(0o700)
            saved_path = os.environ.get("PATH")
            os.environ["PATH"] = str(fake_bin)
            try:
                check(
                    source_provenance(root, commit)["head"] == commit,
                    "source provenance uses absolute Git",
                )
            finally:
                if saved_path is None:
                    del os.environ["PATH"]
                else:
                    os.environ["PATH"] = saved_path
        tracked.write_text("second\n", encoding="utf-8")
        checked_output(["/usr/bin/git", "add", "tracked.txt"], cwd=root)
        checked_output([
            "/usr/bin/git", "-c", "user.name=Leopard2 Self Test", "-c",
            "user.email=leopard2-self-test.invalid", "commit", "-qm",
            "second"], cwd=root)
        second = checked_output(
            ["/usr/bin/git", "rev-parse", "HEAD"], cwd=root)
        checked_output(
            ["/usr/bin/git", "reset", "--hard", "-q", commit], cwd=root)
        expect_runtime_error(
            lambda: source_provenance(root, second),
            "does not exactly match expected commit",
            "wrong expected source commit",
        )
        checked_output(
            ["/usr/bin/git", "replace", commit, second], cwd=root)
        expect_runtime_error(
            lambda: source_provenance(root, commit),
            "replace refs",
            "Git replacement ref",
        )

    check(
        assert_line_numbers("value = 1\n") == (),
        "assert scanner accepts assert-free source",
    )
    check(
        assert_line_numbers("value = 1\nassert value\n") == (2,),
        "assert scanner detects optimized-away checks",
    )
    check(
        assert_line_numbers(Path(__file__).read_text(encoding="utf-8")) == (),
        "expanded final-map runner contains no Python assert statements",
    )
    check(check_count > 1, "self-test executed a nontrivial check set")
    print(json.dumps({
        "self_test": "passed", "matrix_cells": len(matrix["cells"]),
        "matrix_sha256": canonical_digest,
        "r1_tiny_exact_main_cells": len(r1_matrix["cells"]),
        "r1_tiny_exact_main_sha256": matrix_digest(r1_matrix),
        "k4_r2_production_avx2_cells": len(k4_r2_matrix["cells"]),
        "k4_r2_production_avx2_sha256": matrix_digest(k4_r2_matrix),
        "selector_isolation_cells": len(selectors),
        "r1_selector_boundary_cells": len(r1_boundaries),
        "direct_odd_arity_cells": len(direct_odd),
        "sealed_snapshot_execution": "passed",
        "diagnostic_byte_identity": "passed",
        "absolute_git_and_replace_ref_gate": "passed",
        "clean_source_identity": "passed",
        "checks": check_count,
    }, sort_keys=True))


def add_benchmark_arguments(parser):
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument(
        "--candidate-build-executable",
        help=("build-tree executable used to validate CMake/object/archive/"
              "link provenance; defaults to --candidate"))
    parser.add_argument("--candidate-commit", required=True)
    parser.add_argument("--candidate-source", required=True)
    parser.add_argument("--candidate-build", required=True)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument(
        "--lock", default="/tmp/leopard-gf8-authoritative.lock",
        help="canonical lock shared by diagnostics and authoritative runs")
    parser.add_argument("--maximum-attempts", type=int, default=3)


def validate_benchmark_arguments(args):
    if args.maximum_attempts < 1:
        raise RuntimeError("--maximum-attempts must be positive")
    if not COMMIT_PATTERN.match(args.candidate_commit):
        raise RuntimeError("--candidate-commit must be a full lowercase SHA")


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    make = subparsers.add_parser("make-matrix")
    make.add_argument("--output", required=True)
    make.add_argument(
        "--profile",
        choices=(
            CANONICAL_MATRIX_PROFILE,
            R1_TINY_EXACT_MAIN_PROFILE,
            K4_R2_PRODUCTION_AVX2_PROFILE,
        ),
        default=CANONICAL_MATRIX_PROFILE,
        help="deterministic matrix profile (default: canonical)")
    dry = subparsers.add_parser("dry-run")
    dry.add_argument("--matrix", required=True)
    dry.add_argument("--output", required=True)
    dry.add_argument("--diagnostic-summary")
    dry.add_argument("--near-ratio", type=float, default=1.10)
    subparsers.add_parser("self-test")
    diagnostic = subparsers.add_parser("diagnostic")
    add_benchmark_arguments(diagnostic)
    diagnostic.add_argument(
        "--cpu-pairs", help="optional comma-separated cpu:sibling pairs")
    diagnostic.add_argument("--maximum-workers", type=int, default=0)
    diagnostic.add_argument("--memory-budget-gib", type=int, default=0)
    abba = subparsers.add_parser("abba")
    add_benchmark_arguments(abba)
    abba.add_argument("--diagnostic-summary", required=True)
    abba.add_argument("--near-ratio", type=float, default=1.10)
    abba.add_argument("--rounds", type=int, default=3)
    abba.add_argument("--cpu", type=int)
    abba.add_argument("--sibling", type=int)
    abba.add_argument("--include-cell", action="append", default=[])
    abba.add_argument("--exclude-cell", action="append", default=[])
    merge = subparsers.add_parser("merge")
    merge.add_argument("--diagnostic-summary", required=True)
    merge.add_argument("--abba-summary", required=True)
    merge.add_argument("--output", required=True)
    return parser.parse_args()


def main():
    args = parse_arguments()
    if args.command == "self-test":
        self_test()
    elif args.command == "make-matrix":
        matrix = make_matrix(args.profile)
        validate_matrix(matrix)
        atomic_json(args.output, matrix)
        print(json.dumps({"output": os.path.realpath(args.output),
                          "cells": len(matrix["cells"]),
                          "sha256": matrix_digest(matrix)}, sort_keys=True))
    elif args.command == "dry-run":
        matrix = read_json(args.matrix)
        validate_matrix(matrix)
        manifest = make_dry_run_manifest(
            matrix, args.diagnostic_summary, args.near_ratio)
        atomic_json(args.output, manifest)
        print(json.dumps({"output": os.path.realpath(args.output),
                          "sha256": sha256_file(args.output)}, sort_keys=True))
    elif args.command == "diagnostic":
        validate_benchmark_arguments(args)
        if args.maximum_workers < 0 or args.memory_budget_gib < 0:
            raise RuntimeError("worker and memory limits cannot be negative")
        run_diagnostic(args)
    elif args.command == "abba":
        validate_benchmark_arguments(args)
        if args.rounds < 2:
            raise RuntimeError("ABBA requires at least two rounds")
        if args.near_ratio <= 0:
            raise RuntimeError("--near-ratio must be positive")
        run_abba(args)
    elif args.command == "merge":
        run_merge(args)
    else:
        raise RuntimeError("unhandled command")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print("error: " + str(error), file=sys.stderr)
        sys.exit(1)
