#!/usr/bin/env python3
"""Deterministic benchmark-v2 emitter used only by run_abba.py self-test."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path


MASK64 = (1 << 64) - 1


class XorShift64:
    def __init__(self, seed: int):
        self.state = seed & MASK64
        if self.state == 0:
            self.state = 0x9E3779B97F4A7C15

    def next(self) -> int:
        value = self.state
        value ^= (value << 13) & MASK64
        value ^= value >> 7
        value ^= (value << 17) & MASK64
        self.state = value & MASK64
        return self.state


def losses(k: int, count: int, seed: int) -> list[int]:
    order = list(range(k))
    random = XorShift64(seed ^ 0xD1B54A32D192ED03)
    for remaining in range(k, 1, -1):
        selected = random.next() % remaining
        order[remaining - 1], order[selected] = order[selected], order[remaining - 1]
    return sorted(order[:count])


def summary(
    value: float,
    iterations: int,
    setup: bool,
    input_name: str = "",
    output_name: str = "",
    input_bytes: int = 0,
    output_bytes: int = 0,
) -> dict[str, object]:
    samples = [value] * iterations
    if setup:
        return {
            "median_us": value,
            "mad_us": 0.0,
            "minimum_us": value,
            "maximum_us": value,
            "samples_us": samples,
        }
    return {
        "median_us_per_batch_call": value,
        "mad_us_per_batch_call": 0.0,
        "minimum_us_per_batch_call": value,
        "maximum_us_per_batch_call": value,
        "samples_us_per_batch_call": samples,
        input_name: input_bytes / (value * 1000.0),
        output_name: output_bytes / (value * 1000.0),
    }


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument("--k", type=int, required=True)
    result.add_argument("--r", type=int, required=True)
    result.add_argument("--profile", required=True)
    result.add_argument("--field", required=True)
    result.add_argument("--backend", required=True)
    result.add_argument("--bytes", type=int, required=True)
    result.add_argument("--loss", type=int, required=True)
    result.add_argument("--batch", type=int, required=True)
    result.add_argument("--reuse", type=int, required=True)
    result.add_argument("--iterations", type=int, required=True)
    result.add_argument("--warmup", type=int, required=True)
    result.add_argument("--threads", type=int, required=True)
    result.add_argument("--seed", type=int, required=True)
    result.add_argument("--skip-legacy", action="store_true")
    result.add_argument("--retain-samples", action="store_true")
    result.add_argument("--json", required=True)
    return result


def main() -> int:
    options = parser().parse_args()
    role_name = Path(sys.argv[0]).name
    try:
        contents = Path(__file__).read_text(encoding="utf-8")
    except OSError:
        contents = ""
    for line in contents.splitlines():
        if line.startswith("# MOCK_ROLE="):
            role_name = line.split("=", 1)[1]
            break
    candidate = "candidate" in role_name
    slow = "slow" in role_name
    encode = 11.0 if candidate and slow else (8.0 if candidate else 10.0)
    setup = 5.0 if candidate else 5.0
    padded = 1 << (options.k - 1).bit_length()
    parent = 1 << (padded + options.r - 1).bit_length()
    identity = (
        f"{options.k}:{options.r}:{options.bytes}:{options.seed}:"
        f"{options.field}:{options.backend}"
    ).encode("ascii")
    original_digest = hashlib.sha256(b"original:" + identity).hexdigest()[:16]
    parity_digest = hashlib.sha256(b"parity:" + identity).hexdigest()[:16]
    recovered_digest = hashlib.sha256(b"recovered:" + identity).hexdigest()[:16]
    encode_input_bytes = options.k * options.bytes
    encode_output_bytes = options.r * options.bytes
    decode_input_bytes = (options.k - options.loss + options.r) * options.bytes
    decode_output_bytes = options.loss * options.bytes
    amortized_decode = 7.0 + 3.0 / options.reuse
    document = {
        "schema": "leopard2-benchmark-v2",
        "build": {
            "compiler": "mock",
            "compiler_version": "self-test",
            "cplusplus": 202002,
        },
        "parameters": {
            "K": options.k,
            "R": options.r,
            "requested_profile": "low_v1",
            "requested_field": options.field,
            "requested_backend": options.backend,
            "force_generic_decode": False,
            "force_specialized_decode": False,
            "force_tiled_decode": False,
            "force_materialized_decode": False,
            "skip_legacy": options.skip_legacy,
            "retain_samples": options.retain_samples,
            "shard_bytes": options.bytes,
            "loss_count": options.loss,
            "missing_original_indices": losses(options.k, options.loss, options.seed),
            "batch": options.batch,
            "reuse": options.reuse,
            "iterations": options.iterations,
            "warmup": options.warmup,
            "thread_count": options.threads,
            "seed": options.seed,
        },
        "resolved": {
            "profile": "low_v1",
            "field": options.field,
            "backend": options.backend,
            "thread_count": 1,
            "parent_count": parent,
            "padded_side": padded,
        },
        "correctness": {
            "leopard2_round_trip": True,
            "legacy_comparison": None,
        },
        "workload_digests": {
            "algorithm": "fnv1a64",
            "original_data": original_digest,
            "transmitted_parity": parity_digest,
            "recovered_originals": recovered_digest,
        },
        "memory": {
            "scratch_alignment": 64,
            "encode_scratch_bytes_per_stripe": padded * options.bytes,
            "decode_scratch_bytes_per_stripe": parent * options.bytes,
            "encode_scratch_bytes_batch": padded * options.bytes,
            "decode_scratch_bytes_batch": parent * options.bytes,
        },
        "metrics": {
            "codec_setup": summary(setup, options.iterations, True),
            "encode_execution": summary(
                encode, options.iterations, False,
                "input_GB_per_s", "parity_output_GB_per_s",
                encode_input_bytes, encode_output_bytes),
            "decode_plan_setup": summary(3.0, options.iterations, True),
            "decode_execution": summary(
                7.0, options.iterations, False,
                "offered_received_GB_per_s", "repaired_output_GB_per_s",
                decode_input_bytes, decode_output_bytes),
            "decode_amortized_at_reuse": {
                "reuse_count": options.reuse,
                "derived_median_us_per_batch_call": amortized_decode,
                "offered_received_GB_per_s":
                    decode_input_bytes / (amortized_decode * 1000.0),
                "repaired_output_GB_per_s":
                    decode_output_bytes / (amortized_decode * 1000.0),
            },
            "rate_semantics":
                "offered_received counts all non-null shard pointers supplied; "
                "a plan may read a deterministic subset",
        },
        "legacy": {
            "available": False,
            "unavailable_reason": "disabled by --skip-legacy",
            "codec_setup": None,
            "decode_timing_includes_setup": True,
            "encode_execution": None,
            "decode_including_setup": None,
        },
    }
    encoded = json.dumps(document, sort_keys=True, separators=(",", ":"))
    if options.json == "-":
        print(encoded)
    else:
        Path(options.json).write_text(encoded + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
