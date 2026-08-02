#!/usr/bin/env python3
"""Screen one K range of the dense full-parity GF8 AVX2 terminal.

Build two executable-section-identical Leopard2 benchmark binaries with
LEO2_DIAGNOSTIC_DENSE_FULL_PARITY_ONE_BLOCK_MODE set to one and two,
respectively.  Freeze and hash them before invoking this runner.  The runner
reuses the hardened pinned ABBA machinery and compares every selected K with
its same-source control.  Exact Leopard main is additionally invoked only for
the public K >= R cells accepted by its API; K < R targets use candidate/control
orders and analysis only.

The side, byte size, and K interval are explicit so the 480-cell B64/B256
matrix can be sharded into independently retained result directories.  This
screen deliberately applies no promotion threshold: it records crossover
data from which contiguous production calibration intervals can be derived.
Authoritative confirmation of any proposed interval remains a separate
nine-round campaign with the ordinary five-percent promotion gate.
"""

from __future__ import annotations

import argparse
import fcntl
import importlib.util
import os
import sys
import tempfile
from pathlib import Path
from typing import Any, Sequence


def load_base() -> Any:
    path = Path(__file__).resolve().with_name(
        "run_k5r5_b64_terminal_abba.py")
    specification = importlib.util.spec_from_file_location(
        "dense_full_parity_screen_base", path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load evidence base: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


BASE = load_base()
BASE.__doc__ = __doc__
BASE.SCHEMA = "leopard2-gf8-dense-full-parity-one-block-screen/v1"
BASE.SUMMARY_SCHEMA = \
    "leopard2-gf8-dense-full-parity-one-block-screen-summary/v1"
BASE.MODE_SYMBOL = \
    "_ZN12_GLOBAL__N_1L34g_dense_full_parity_one_block_modeE"
BASE.ALLOW_MULTIPLE_TARGETS = True
BASE.TARGET_ROLES = ("target_main", "target_control")
BASE.MAIN_TARGET_ROLES = ("target_main",)
BASE.TARGET_CONTROL_FLOOR = 0.0
BASE.TARGET_MAIN_FLOOR = 0.0
BASE.REQUIRE_NORMALIZED_FULL_FILE_EQUIVALENCE = True
# This is a saturated all-physical-core crossover diagnostic, not promotion
# evidence.  Retain every isolation record but do not reject a shard merely
# because unrelated activity reached a measured CPU or its SMT sibling.
BASE.REQUIRE_STRICT_ISOLATION = False
BASE.REQUIRE_STABLE_LEASE_ANCHOR = False
BALANCED_MODE_SYMBOL = \
    "_ZN12_GLOBAL__N_1L28g_balanced_b64_terminal_modeE"
BASE.AUXILIARY_MODE_EXPECTATIONS = {
    BALANCED_MODE_SYMBOL: {"candidate": 1, "control": 1},
}
BASE.RUNNER_PATH = Path(__file__).resolve()
BASE.RUNNER_DEPENDENCIES = (
    BASE.RUNNER_PATH,
    Path(BASE.__file__).resolve(),
    Path(BASE.T8_SUPPORT.__file__).resolve(),
    Path(BASE.MAIN_SUPPORT.__file__).resolve(),
)


SCREEN: argparse.Namespace | None = None


def parse_screen_arguments(arguments: Sequence[str]) -> tuple[
        argparse.Namespace, list[str]]:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--screen-side", type=int,
                        choices=(16, 32, 64, 128))
    parser.add_argument("--screen-bytes", type=int,
                        choices=(64, 256))
    parser.add_argument("--screen-k-start", type=int, default=1)
    parser.add_argument("--screen-k-end", type=int)
    return parser.parse_known_args(list(arguments))


def target_cells(
    side: int,
    shard_bytes: int,
    first: int = 1,
    last: int | None = None,
) -> list[dict[str, Any]]:
    last = side if last is None else last
    BASE.require(1 <= first <= last <= side,
                 "screen K interval must lie within 1..T")

    result: list[dict[str, Any]] = []
    for original_count in range(first, last + 1):
        role = "target_main" if original_count >= side else "target_control"
        result.append({
            "id": (f"{role}-k{original_count}-r{side}-"
                   f"b{shard_bytes}-q1"),
            "K": original_count,
            "R": side,
            "bytes": shard_bytes,
            "batch": 1,
            "reuse": 8192,
            "role": role,
            "seed": (0x44504600 + side * 0x10000 +
                     shard_bytes * 0x100 + original_count),
        })

    return result


def cells() -> list[dict[str, Any]]:
    BASE.require(SCREEN is not None, "screen arguments were not initialized")
    BASE.require(SCREEN.screen_side is not None and
                 SCREEN.screen_bytes is not None,
                 "--screen-side and --screen-bytes are required")
    side = int(SCREEN.screen_side)
    shard_bytes = int(SCREEN.screen_bytes)
    first = int(SCREEN.screen_k_start)
    last = side if SCREEN.screen_k_end is None else int(SCREEN.screen_k_end)
    result = target_cells(side, shard_bytes, first, last)

    # A bounded set of selector neighbors travels with every shard.  Duplicate
    # IDs are impossible because target byte/count tuples differ from these.
    sentinels = sorted({1, max(1, side // 2), max(1, side - 1), side})
    for original_count in sentinels:
        for neighbor_bytes in (shard_bytes - 1, shard_bytes + 1):
            result.append({
                "id": (f"byte-neighbor-k{original_count}-r{side}-"
                       f"b{neighbor_bytes}-q1"),
                "K": original_count,
                "R": side,
                "bytes": neighbor_bytes,
                "batch": 1,
                "reuse": 8192,
                "role": "neighbor",
                "seed": (0x44504E00 + side * 0x10000 +
                         neighbor_bytes * 0x100 + original_count),
            })
        result.append({
            "id": (f"r-neighbor-k{original_count}-r{side - 1}-"
                   f"b{shard_bytes}-q1"),
            "K": original_count,
            "R": side - 1,
            "bytes": shard_bytes,
            "batch": 1,
            "reuse": 8192,
            "role": "neighbor",
            "seed": (0x44505200 + side * 0x10000 +
                     shard_bytes * 0x100 + original_count),
        })
    if side < 128:
        result.append({
            "id": f"k-neighbor-k{side + 1}-r{side}-b{shard_bytes}-q1",
            "K": side + 1,
            "R": side,
            "bytes": shard_bytes,
            "batch": 1,
            "reuse": 8192,
            "role": "neighbor",
            "seed": 0x44504B00 + side * 0x10000 + shard_bytes * 0x100,
        })
    return result


BASE.cells = cells


def synthetic_rounds(order: Sequence[str]) -> list[dict[str, Any]]:
    digests = {
        "algorithm": "fnv1a64",
        "original_data": "0000000000000001",
        "transmitted_parity": "0000000000000002",
        "recovered_originals": "0000000000000003",
    }
    timings = {"candidate": 1.0, "control": 1.1, "main": 1.2}
    return [{
        "round": round_index,
        "order": list(order),
        "invocations": [{
            "implementation": label,
            "normalized": {
                "digests": dict(digests),
                "encode_us": timings[label],
            },
        } for label in order],
        "isolation": {"accepted": True},
    } for round_index in range(3)]


def self_test() -> int:
    targets = [
        cell
        for side in (16, 32, 64, 128)
        for shard_bytes in (64, 256)
        for cell in target_cells(side, shard_bytes)
    ]
    BASE.require(len(targets) == 480,
                 "full dense screen does not contain 480 target cells")
    BASE.require(sum(cell["role"] == "target_main" for cell in targets) == 8,
                 "exact-main-compatible dense target count changed")
    BASE.require(sum(
        cell["role"] == "target_control" for cell in targets) == 472,
        "same-source-only dense target count changed")
    control_cell = next(
        cell for cell in targets if cell["role"] == "target_control")
    main_cell = next(cell for cell in targets if cell["role"] == "target_main")
    BASE.require(BASE.cell_is_target(control_cell) and
                 not BASE.cell_uses_main(control_cell),
                 "K<R target was not classified same-source-only")
    BASE.require(BASE.cell_is_target(main_cell) and
                 BASE.cell_uses_main(main_cell),
                 "K>=R target was not classified exact-main-compatible")
    control_analysis = BASE.analyze(
        control_cell, synthetic_rounds(BASE.NEIGHBOR_ORDER[0]))
    BASE.require("main_over_candidate" not in control_analysis,
                 "K<R analysis unexpectedly requested exact main")
    main_analysis = BASE.analyze(
        main_cell, synthetic_rounds(BASE.TARGET_ORDER[0]))
    BASE.require("main_over_candidate" in main_analysis,
                 "K>=R analysis omitted exact main")
    contaminated = synthetic_rounds(BASE.NEIGHBOR_ORDER[0])
    for round_value in contaminated:
        round_value["isolation"]["accepted"] = False
    BASE.analyze(control_cell, contaminated)
    BASE.require(BASE.REQUIRE_NORMALIZED_FULL_FILE_EQUIVALENCE,
                 "dense screen does not require full-file provenance")
    BASE.require(not BASE.REQUIRE_STRICT_ISOLATION,
                 "saturated dense diagnostic unexpectedly claims isolation")
    BASE.require(not BASE.REQUIRE_STABLE_LEASE_ANCHOR,
                 "saturated dense diagnostic unexpectedly serializes shards")
    BASE.require(BASE.AUXILIARY_MODE_EXPECTATIONS == {
        BALANCED_MODE_SYMBOL: {"candidate": 1, "control": 1}},
        "balanced B64 selector provenance is not pinned enabled/equal")
    with tempfile.TemporaryDirectory(
            prefix="leopard-dense-full-file-self-test-") as directory:
        root = Path(directory)
        candidate = root / "candidate"
        control = root / "control"
        candidate.write_bytes(b"sameAtext-build-id-X")
        control.write_bytes(b"sameBtext-build-id-Y")
        comparison = BASE.normalized_full_file_comparison(
            candidate, control, (
                {"name": "diagnostic_selector", "file_offset": 4, "size": 1},
                {"name": ".note.gnu.build-id", "file_offset": 19, "size": 1},
            ))
        BASE.require(comparison["difference_count"] == 2,
                     "allowed full-file differences were not retained")
        rejected = False
        try:
            BASE.normalized_full_file_comparison(
                candidate, control, (
                    {"name": "diagnostic_selector",
                     "file_offset": 4, "size": 1},
                ))
        except BASE.EvidenceError:
            rejected = True
        BASE.require(rejected,
                     "unauthorized full-file difference was accepted")
    with tempfile.TemporaryDirectory(
            prefix="leopard-dense-lock-self-test-") as directory:
        original_lock_path = BASE.LOCK_PATH
        lock_path = Path(directory) / "campaign.lock"
        BASE.LOCK_PATH = lock_path
        owner = os.open(lock_path, os.O_RDWR | os.O_CREAT, 0o600)
        adopted = None
        try:
            fcntl.flock(owner, fcntl.LOCK_EX | fcntl.LOCK_NB)
            adopted = BASE.acquire_global_lock(owner)
            os.close(owner)
            owner = -1
            contender = os.open(lock_path, os.O_RDWR)
            blocked = False
            try:
                try:
                    fcntl.flock(contender,
                                fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError:
                    blocked = True
            finally:
                os.close(contender)
            BASE.require(blocked,
                         "adopted descriptor did not retain umbrella lease")
            os.close(adopted)
            adopted = None
            contender = os.open(lock_path, os.O_RDWR)
            try:
                fcntl.flock(contender, fcntl.LOCK_EX | fcntl.LOCK_NB)
                fcntl.flock(contender, fcntl.LOCK_UN)
            finally:
                os.close(contender)
        finally:
            if adopted is not None:
                os.close(adopted)
            if owner >= 0:
                os.close(owner)
            BASE.LOCK_PATH = original_lock_path
    print("dense full-parity runner self-test passed")
    return 0


def main() -> int:
    global SCREEN
    SCREEN, remaining = parse_screen_arguments(sys.argv[1:])
    if SCREEN.self_test:
        BASE.require(not remaining,
                     "--self-test does not accept benchmark arguments")
        return self_test()
    BASE.require(SCREEN.screen_side is not None and
                 SCREEN.screen_bytes is not None,
                 "--screen-side and --screen-bytes are required")
    sys.argv = [sys.argv[0]] + remaining
    return BASE.main()


if __name__ == "__main__":
    raise SystemExit(main())
