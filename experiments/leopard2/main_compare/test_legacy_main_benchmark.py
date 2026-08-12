#!/usr/bin/env python3
"""Smoke and schema checks for the exact Leopard-main benchmark adapter."""

from __future__ import annotations

import argparse
import json
import math
import shlex
import subprocess
from pathlib import Path


COMMIT = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
MASK64 = (1 << 64) - 1


class XorShift64:
    def __init__(self, seed: int):
        self.state = seed or 0x9E3779B97F4A7C15

    def next(self) -> int:
        value = self.state
        value ^= (value << 13) & MASK64
        value ^= value >> 7
        value ^= (value << 17) & MASK64
        self.state = value & MASK64
        return self.state


def fnv1a64(chunks) -> str:
    value = 14695981039346656037
    for chunk in chunks:
        for byte in chunk:
            value ^= byte
            value = (value * 1099511628211) & MASK64
    return f"{value:016x}"


def original_shards(k: int, batch: int, byte_count: int, seed: int):
    stripes = []
    for stripe in range(batch):
        random = XorShift64(
            seed ^ ((0x9E3779B97F4A7C15 * (stripe + 1)) & MASK64))
        data = bytes(random.next() >> 56 for _ in range(k * byte_count))
        stripes.append([
            data[index * byte_count:(index + 1) * byte_count]
            for index in range(k)
        ])
    return stripes


def run(executable: Path, arguments: list[str], success: bool = True):
    completed = subprocess.run(
        [str(executable), *arguments], text=True, capture_output=True,
        timeout=60, check=False)
    if (completed.returncode == 0) != success:
        raise RuntimeError(
            f"unexpected exit {completed.returncode}: {completed.stderr}")
    return completed


def positive_summary(value: object) -> None:
    if not isinstance(value, dict):
        raise RuntimeError("timing summary is not an object")
    for name in (
        "median_us_per_batch_call", "mad_us_per_batch_call",
        "minimum_us_per_batch_call", "maximum_us_per_batch_call",
    ):
        number = value.get(name)
        if not isinstance(number, (int, float)) or not math.isfinite(number):
            raise RuntimeError(f"invalid {name}")
        if name != "mad_us_per_batch_call" and number <= 0:
            raise RuntimeError(f"non-positive {name}")
    samples = value.get("samples_us_per_batch_call")
    if (not isinstance(samples, list) or not samples or
            any(not isinstance(sample, (int, float)) or
                not math.isfinite(sample) or sample <= 0 for sample in samples)):
        raise RuntimeError("invalid raw timing samples")
    if len(samples) != 2:
        raise RuntimeError("raw timing sample count does not match iterations")
    ordered = sorted(samples)
    median = (ordered[0] + ordered[1]) * 0.5
    deviations = sorted(abs(sample - median) for sample in samples)
    mad = (deviations[0] + deviations[1]) * 0.5
    expected = {
        "median_us_per_batch_call": median,
        "mad_us_per_batch_call": mad,
        "minimum_us_per_batch_call": min(samples),
        "maximum_us_per_batch_call": max(samples),
    }
    for name, derived in expected.items():
        if abs(value[name] - derived) > 0.000002:
            raise RuntimeError(f"{name} is not derived from retained samples")


def check_case(executable: Path, k: int, r: int, loss: int, seed: int,
               pure_avx2: bool) -> None:
    arguments = [
        "--k", str(k), "--r", str(r), "--bytes", "64",
        "--loss", str(loss), "--batch", "2", "--reuse", "1",
        "--iterations", "2", "--warmup", "1", "--threads", "1",
        "--seed", str(seed), "--json", "-",
    ]
    first = run(executable, arguments)
    second = run(executable, arguments)
    one = json.loads(first.stdout)
    two = json.loads(second.stdout)
    if one.get("schema") != "leopard-main-benchmark-v1":
        raise RuntimeError("wrong schema")
    if one.get("build", {}).get("main_source_commit") != COMMIT:
        raise RuntimeError("wrong source commit")
    if one.get("build", {}).get("pure_avx2") is not pure_avx2:
        raise RuntimeError("wrong exact-main ISA-profile attestation")
    parameters = one.get("parameters", {})
    if (parameters.get("K"), parameters.get("R"),
            parameters.get("loss_count"),
            parameters.get("logical_shard_bytes")) != (k, r, loss, 64):
        raise RuntimeError("parameter echo mismatch")
    missing = parameters.get("missing_original_indices")
    if (not isinstance(missing, list) or len(missing) != loss or
            missing != sorted(set(missing)) or
            missing != two.get("parameters", {}).get("missing_original_indices")):
        raise RuntimeError("missing-set determinism failed")
    if one.get("correctness", {}).get("round_trip") is not True:
        raise RuntimeError("round trip was not reported")
    if one.get("resolved", {}).get("thread_count") != 1:
        raise RuntimeError("resolved thread count was not reported truthfully")
    digests = one.get("workload_digests", {})
    if digests.get("algorithm") != "fnv1a64":
        raise RuntimeError("wrong workload digest algorithm")
    for name in ("original_data", "transmitted_parity", "recovered_originals"):
        digest = digests.get(name)
        if (not isinstance(digest, str) or len(digest) != 16 or
                any(character not in "0123456789abcdef" for character in digest)):
            raise RuntimeError(f"invalid workload digest {name}")
        if digest != two.get("workload_digests", {}).get(name):
            raise RuntimeError(f"nondeterministic workload digest {name}")
    originals = original_shards(k, 2, 64, seed)
    if digests["original_data"] != fnv1a64(
            shard for stripe in originals for shard in stripe):
        raise RuntimeError("original-data digest disagrees with independent generator")
    if digests["recovered_originals"] != fnv1a64(
            originals[stripe][index]
            for stripe in range(2) for index in missing):
        raise RuntimeError("recovered digest disagrees with missing originals")
    metrics = one.get("metrics", {})
    if metrics.get("decode_timing_includes_setup") is not True:
        raise RuntimeError("decode setup semantics are missing")
    positive_summary(metrics.get("encode_execution"))
    positive_summary(metrics.get("decode_including_setup"))


def check_padded_case(
        executable: Path, k: int, r: int, logical_bytes: int,
        loss: int, seed: int) -> None:
    arguments = [
        "--k", str(k), "--r", str(r), "--bytes", "64",
        "--logical-bytes", str(logical_bytes),
        "--loss", str(loss), "--batch", "2", "--reuse", "1",
        "--iterations", "2", "--warmup", "1", "--threads", "1",
        "--seed", str(seed), "--json", "-",
    ]
    first = json.loads(run(executable, arguments).stdout)
    second = json.loads(run(executable, arguments).stdout)
    if first.get("schema") != "leopard-main-benchmark-v2":
        raise RuntimeError("wrong padded-comparison schema")
    parameters = first.get("parameters", {})
    if (
        parameters.get("shard_bytes"),
        parameters.get("logical_shard_bytes"),
    ) != (64, logical_bytes):
        raise RuntimeError("padded byte-count echo mismatch")
    resolved = first.get("resolved", {})
    if (
        resolved.get("padded_application_bytes") is not True or
        resolved.get("padding_policy") != "zero suffix per shard"
    ):
        raise RuntimeError("padded comparison semantics are absent")
    if first.get("correctness", {}).get("round_trip") is not True:
        raise RuntimeError("padded round trip failed")
    if (
        first.get("correctness", {}).get(
            "logical_prefix_fingerprinted") is not True
    ):
        raise RuntimeError("logical-prefix correctness evidence is absent")

    digests = first.get("workload_digests", {})
    if digests != second.get("workload_digests"):
        raise RuntimeError("padded workload digests are nondeterministic")
    originals = original_shards(k, 2, logical_bytes, seed)
    if digests.get("original_data") != fnv1a64(
            shard for stripe in originals for shard in stripe):
        raise RuntimeError(
            "padded original digest disagrees with independent generator")
    missing = parameters.get("missing_original_indices")
    if digests.get("recovered_originals") != fnv1a64(
            originals[stripe][index]
            for stripe in range(2) for index in missing):
        raise RuntimeError(
            "padded recovered digest disagrees with logical originals")
    positive_summary(first.get("metrics", {}).get("encode_execution"))


def check_compile_commands(path: Path, pure_avx2: bool) -> None:
    commands = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(commands, list):
        raise RuntimeError("compile_commands is not an array")
    required_sources = {
        "leopard.cpp", "LeopardCommon.cpp", "LeopardFF8.cpp",
        "LeopardFF16.cpp", "legacy_main_benchmark.cpp",
    }
    selected = {}
    for row in commands:
        source = Path(row.get("file", "")).name
        if source in required_sources:
            selected[source] = shlex.split(row.get("command", ""))
    if set(selected) != required_sources:
        raise RuntimeError("compile_commands omits an exact-main source")
    isa_flags = ({
        "-march=x86-64", "-mtune=generic", "-mavx2", "-mno-avx512f",
    } if pure_avx2 else {"-march=native"})
    required_flags = isa_flags | {
        "-Wall", "-Wextra", "-fopenmp", "-g", "-O0", "-O3",
    }
    for source, command in selected.items():
        missing = required_flags - set(command)
        if missing:
            raise RuntimeError(
                f"{source} is missing canonical flags: {sorted(missing)}")
        if "-DNDEBUG" in command:
            raise RuntimeError(f"{source} unexpectedly defines NDEBUG")
        if pure_avx2 and any(flag.startswith("-mavx512") for flag in command):
            raise RuntimeError(f"{source} unexpectedly enables AVX-512")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", required=True, type=Path)
    parser.add_argument("--compile-commands", required=True, type=Path)
    parser.add_argument("--pure-avx2", action="store_true")
    options = parser.parse_args()
    executable = options.benchmark.resolve()
    check_compile_commands(
        options.compile_commands.resolve(), options.pure_avx2)
    check_case(executable, 8, 4, 0, 1, options.pure_avx2)
    check_case(executable, 8, 4, 4, 7, options.pure_avx2)
    check_case(executable, 129, 1, 1, 11, options.pure_avx2)
    check_case(executable, 240, 16, 3, 19, options.pure_avx2)
    check_padded_case(executable, 8, 4, 1, 1, 23)
    check_padded_case(executable, 8, 4, 33, 4, 29)
    check_padded_case(executable, 16, 8, 63, 1, 31)
    run(executable, ["--k", "8", "--r", "9"], success=False)
    run(executable, ["--k", "8", "--r", "4", "--bytes", "65"],
        success=False)
    run(executable, [
        "--k", "8", "--r", "4", "--bytes", "64",
        "--logical-bytes", "0",
    ], success=False)
    run(executable, [
        "--k", "8", "--r", "4", "--bytes", "64",
        "--logical-bytes", "65",
    ], success=False)
    run(executable, ["--k", "8", "--r", "4", "--loss", "5"],
        success=False)
    run(executable, ["--k", str((1 << 32) - 1), "--r", "1"],
        success=False)
    print("exact Leopard-main benchmark smoke passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
