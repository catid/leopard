#!/usr/bin/env python3
"""Qualify the packed T=R=32/B=256 multi-block AVX2 encoder.

The campaign compares one source-attested candidate with a same-source
selector-disabled control and exact Leopard main.  Exact K=64, 96, and 128
cells must clear both batch-call and one-shot five-percent floors.  Immediate
byte, K, and R neighbors must remain inert.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import sys
from pathlib import Path
from typing import Any


def load_base() -> Any:
    path = Path(__file__).resolve().with_name(
        "run_t32_b256_generated_abba.py")
    specification = importlib.util.spec_from_file_location(
        "t32_b256_multiblock_evidence_base", path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load evidence base: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


PARENT = load_base()
BASE = PARENT.BASE
BASE.__doc__ = __doc__
BASE.SCHEMA = "leopard2-gf8-t32-b256-multiblock-abba/v1"
BASE.SUMMARY_SCHEMA = "leopard2-gf8-t32-b256-multiblock-summary/v1"
BASE.MODE_SYMBOL = \
    "_ZN7leopard7backend12_GLOBAL__N_1L26g_t32_b256_multiblock_modeE"
BASE.TARGET_CONTROL_FLOOR = 1.05
BASE.TARGET_MAIN_FLOOR = 1.05
BASE.ALLOW_MULTIPLE_TARGETS = True
BASE.REQUIRE_EXPECTED_IDENTITIES = True
BASE.REQUIRE_BUILD_CLOSURE = True
BASE.REQUIRE_FULL_ELF_IDENTITY = True
BASE.RUNNER_PATH = Path(__file__).resolve()
BASE.RUNNER_DEPENDENCIES = (
    BASE.RUNNER_PATH,
    Path(PARENT.__file__).resolve(),
    Path(BASE.__file__).resolve(),
    Path(BASE.T8_SUPPORT.__file__).resolve(),
    Path(BASE.MAIN_SUPPORT.__file__).resolve(),
)


def cells() -> list[dict[str, Any]]:
    values = [
        ("target-k64-r32-b256-l1", 64, 32, 256, 1, "target"),
        ("target-k96-r32-b256-l8", 96, 32, 256, 8, "target"),
        ("target-k128-r32-b256-l16", 128, 32, 256, 16, "target"),
        ("byte-neighbor-k64-r32-b255", 64, 32, 255, 1, "neighbor"),
        ("byte-neighbor-k64-r32-b257", 64, 32, 257, 1, "neighbor"),
        ("byte-neighbor-k64-r32-b64", 64, 32, 64, 1, "neighbor"),
        ("shape-neighbor-k63-r32-b256", 63, 32, 256, 1, "neighbor"),
        ("shape-neighbor-k65-r32-b256", 65, 32, 256, 1, "neighbor"),
        ("shape-neighbor-k95-r32-b256", 95, 32, 256, 1, "neighbor"),
        ("shape-neighbor-k97-r32-b256", 97, 32, 256, 1, "neighbor"),
        ("shape-neighbor-k127-r32-b256", 127, 32, 256, 1, "neighbor"),
        ("shape-neighbor-k129-r32-b256", 129, 32, 256, 1, "neighbor"),
        ("shape-neighbor-k64-r31-b256", 64, 31, 256, 1, "neighbor"),
        ("shape-neighbor-k64-r33-b256", 64, 33, 256, 1, "neighbor"),
    ]
    result = []
    for index, (name, k, r, shard_bytes, loss, role) in enumerate(values):
        result.append({
            "id": name,
            "K": k,
            "R": r,
            "bytes": shard_bytes,
            "loss": loss,
            "batch": 1,
            "measure_one_shot": True,
            "reuse": 8192,
            "role": role,
            "compare_main": True,
            "seed": 0x54324D00 + index,
        })
    BASE.require(
        len(result) == 14 and
        sum(cell["role"] == "target" for cell in result) == 3 and
        len({cell["id"] for cell in result}) == len(result) and
        len({cell["seed"] for cell in result}) == len(result) and
        all(1 <= cell["loss"] <= cell["R"] for cell in result),
        "T32 multi-block qualification matrix is incomplete")
    return result


def normalized_compile_commands_identity(
    candidate: Path,
    control: Path,
) -> dict[str, Any]:
    def normalize(path: Path, diagnostic_value: str) -> tuple[list[Any], str]:
        value = json.loads(path.read_text(encoding="utf-8"))
        BASE.require(isinstance(value, list) and value,
                     "compile_commands is not a nonempty array")
        directories = {
            row.get("directory") for row in value if isinstance(row, dict)
        }
        BASE.require(len(directories) == 1 and
                     all(isinstance(row, dict) for row in value),
                     "compile_commands has ambiguous build roots")
        build_root = next(iter(directories))
        BASE.require(isinstance(build_root, str) and build_root,
                     "compile_commands build root is invalid")
        normalized = []
        t32_rows = []
        for row in value:
            current = dict(row)
            for key, item in tuple(current.items()):
                if isinstance(item, str):
                    item = item.replace(build_root, "${BUILD}")
                    item = re.sub(
                        r"LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_MULTIBLOCK="
                        + re.escape(diagnostic_value) + r"\b",
                        "LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_MULTIBLOCK=${MODE}",
                        item)
                    expected_route = "1" if diagnostic_value == "0" else "0"
                    item = re.sub(
                        r"LEO2_EXPECT_T32_B256_MULTIBLOCK=" +
                        expected_route + r"\b",
                        "LEO2_EXPECT_T32_B256_MULTIBLOCK=${EXPECT}", item)
                    item = re.sub(
                        r"LEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256=\\\""
                        r"[0-9a-f]{64}\\\"",
                        "LEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256=\\\"${CONFIG}\\\"",
                        item)
                    current[key] = item
            normalized.append(current)
            if Path(str(row.get("file", ""))).name == \
                    "Leopard2BackendAVX2T32B256.cpp":
                t32_rows.append(str(row.get("command", "")))
        BASE.require(
            len(t32_rows) == 1 and
            ("LEO2_DIAGNOSTIC_DISABLE_HIGH_T32_B256_MULTIBLOCK=" +
             diagnostic_value) in t32_rows[0] and
            "-mavx2" in t32_rows[0] and "-mno-avx512f" in t32_rows[0],
            "T32 multi-block compile-command contract changed")
        canonical = json.dumps(
            normalized, sort_keys=True, separators=(",", ":"))
        return normalized, hashlib.sha256(canonical.encode()).hexdigest()

    candidate_value, candidate_hash = normalize(candidate, "0")
    control_value, control_hash = normalize(control, "1")
    BASE.require(
        candidate_value == control_value and candidate_hash == control_hash,
        "candidate/control normalized compile commands differ")
    return {
        "normalized_sha256": candidate_hash,
        "entry_count": len(candidate_value),
    }


BASE.cells = cells
BASE.normalized_compile_commands_identity = \
    normalized_compile_commands_identity


if __name__ == "__main__":
    raise SystemExit(BASE.main())
