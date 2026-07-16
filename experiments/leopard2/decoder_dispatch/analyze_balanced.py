#!/usr/bin/env python3
"""Validate and merge the isolated balanced-decoder crossover runs."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import tempfile
from pathlib import Path


BACKENDS = ("scalar", "ssse3", "avx2")
SIZES = (64, 256, 4096, 65536, 1048576)
MODES = ("specialized", "generic")


class EvidenceError(ValueError):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_result(path: Path, backend: str, shard_bytes: int, mode: str) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise EvidenceError(f"cannot load {path}: {error}") from error
    if value.get("schema") != "leopard2-benchmark-v1":
        raise EvidenceError(f"unexpected schema in {path}")
    parameters = value.get("parameters", {})
    resolved = value.get("resolved", {})
    expected_generic = mode == "generic"
    checks = {
        "K": 128,
        "R": 128,
        "requested_profile": "legacy_high_v1",
        "requested_field": "gf8",
        "requested_backend": backend,
        "force_generic_decode": expected_generic,
        "shard_bytes": shard_bytes,
        "loss_count": 128,
        "batch": 1,
        "reuse": 4,
        "iterations": 15,
        "warmup": 3,
        "thread_count": 1,
        "seed": 2129325312,
    }
    for key, expected in checks.items():
        if parameters.get(key) != expected:
            raise EvidenceError(
                f"{path}: parameter {key}={parameters.get(key)!r}, expected {expected!r}")
    if parameters.get("missing_original_indices") != list(range(128)):
        raise EvidenceError(f"{path}: full-loss coordinate list is not canonical")
    if (resolved.get("profile"), resolved.get("field"), resolved.get("backend"),
            resolved.get("parent_count"), resolved.get("padded_side")) != (
                "legacy_high_v1", "gf8", backend, 256, 128):
        raise EvidenceError(f"{path}: resolved codec identity differs")
    if value.get("correctness", {}).get("leopard2_round_trip") is not True:
        raise EvidenceError(f"{path}: round trip did not pass")
    metric = value.get("metrics", {}).get("decode_execution", {})
    median = float(metric.get("median_us_per_batch_call", 0.0))
    mad = float(metric.get("mad_us_per_batch_call", -1.0))
    if not math.isfinite(median) or not math.isfinite(mad) or median <= 0 or mad < 0:
        raise EvidenceError(f"{path}: invalid timing statistic")
    setup_metric = value.get("metrics", {}).get("decode_plan_setup", {})
    setup_median = float(setup_metric.get("median_us", 0.0))
    setup_mad = float(setup_metric.get("mad_us", -1.0))
    if (not math.isfinite(setup_median) or not math.isfinite(setup_mad) or
            setup_median <= 0 or setup_mad < 0):
        raise EvidenceError(f"{path}: invalid decode-plan setup statistic")
    return {
        "median_us": median,
        "mad_us": mad,
        "decode_plan_setup_median_us": setup_median,
        "decode_plan_setup_mad_us": setup_mad,
        "sha256": sha256_file(path),
    }


def parse_binary(value: str) -> tuple[str, str]:
    try:
        backend, digest = value.split("=", 1)
    except ValueError as error:
        raise argparse.ArgumentTypeError("binary identity must be BACKEND=SHA256") from error
    if backend not in BACKENDS or len(digest) != 64:
        raise argparse.ArgumentTypeError("invalid backend or SHA-256")
    try:
        int(digest, 16)
    except ValueError as error:
        raise argparse.ArgumentTypeError("invalid SHA-256") from error
    return backend, digest.lower()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as output:
            json.dump(value, output, indent=2, sort_keys=True)
            output.write("\n")
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--binary", action="append", type=parse_binary, required=True)
    args = parser.parse_args()
    binaries = dict(args.binary)
    if set(binaries) != set(BACKENDS) or len(args.binary) != len(BACKENDS):
        parser.error("provide exactly one --binary for scalar, ssse3, and avx2")

    root = Path(args.input_root)
    cells = []
    for backend in BACKENDS:
        for shard_bytes in SIZES:
            runs = []
            for run in (1, 2):
                values = {}
                for mode in MODES:
                    path = root / f"run{run}" / f"{backend}.{shard_bytes}.{mode}.json"
                    values[mode] = load_result(path, backend, shard_bytes, mode)
                ratio = values["generic"]["median_us"] / values["specialized"]["median_us"]
                values["generic_over_specialized"] = ratio
                runs.append(values)
            mean_ratio = sum(run["generic_over_specialized"] for run in runs) / len(runs)
            maximum_mad = max(
                100.0 * run[mode]["mad_us"] / run[mode]["median_us"]
                for run in runs for mode in MODES)
            cells.append({
                "backend": backend,
                "shard_bytes": shard_bytes,
                "runs": runs,
                "mean_generic_gain_percent": 100.0 * (1.0 - mean_ratio),
                "maximum_mad_percent": maximum_mad,
            })

    output = {
        "schema": "leopard2-balanced-dispatch-evidence/v1",
        "source_commit": args.source_commit,
        "binary_sha256": binaries,
        "method": {
            "cpu_set": [0],
            "runs": 2,
            "reversed_order": True,
            "run_order": {
                "run1": {
                    "backend_order": list(BACKENDS),
                    "size_order": list(SIZES),
                    "mode_order": ["specialized", "generic"],
                },
                "run2": {
                    "backend_order": list(BACKENDS),
                    "size_order": list(reversed(SIZES)),
                    "mode_order": ["generic", "specialized"],
                },
            },
            "samples_per_cell_per_run": 15,
            "warmups": 3,
            "batch": 1,
            "reuse": 4,
            "threads": 1,
            "seed": 2129325312,
            "profile": "legacy_high_v1",
            "field": "gf8",
            "K": 128,
            "R": 128,
            "missing_originals": 128,
            "shard_bytes": list(SIZES),
        },
        "cells": cells,
        "status": "pass",
    }
    canonical = json.dumps(output, sort_keys=True, separators=(",", ":")).encode("ascii")
    output["content_sha256"] = hashlib.sha256(canonical).hexdigest()
    write_json(Path(args.output), output)
    print(json.dumps({"cells": len(cells), "content_sha256": output["content_sha256"]},
                     sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
