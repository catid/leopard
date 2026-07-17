#!/usr/bin/env python3
"""Validate and summarize pinned materialized-versus-tiled high decode evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import tempfile


SCHEMA = "leopard2-tiling-followup/v1"
BENCHMARK_SCHEMA = "leopard2-benchmark-v2"
SUMMARY_SCHEMA = "leopard2-tiled-high-dispatch-evidence/v1"
BACKENDS = {"scalar", "ssse3", "avx2"}
# Two-sided Student-t 95% critical values.  Production evidence currently uses
# three independent ABBA rounds (df=2), but the small table keeps replay clear.
T95 = {
    1: 12.706204736,
    2: 4.302652730,
    3: 3.182446305,
    4: 2.776445105,
    5: 2.570581836,
    6: 2.446911851,
    7: 2.364624252,
    8: 2.306004135,
    9: 2.262157163,
    10: 2.228138852,
    11: 2.200985160,
    12: 2.178812830,
    13: 2.160368656,
    14: 2.144786688,
    15: 2.131449546,
}


class EvidenceError(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvidenceError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path, label: str) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise EvidenceError(f"cannot read {label} {path}: {error}") from error
    require(isinstance(value, dict), f"{label} is not a JSON object")
    return value


def hex_digest(value: object, label: str) -> str:
    require(isinstance(value, str) and len(value) == 64,
            f"{label} is not a SHA-256")
    try:
        int(value, 16)
    except ValueError as error:
        raise EvidenceError(f"{label} is not hexadecimal") from error
    return value.lower()


def git_oid(value: object, label: str) -> str:
    require(isinstance(value, str) and len(value) == 40,
            f"{label} is not a full Git object ID")
    try:
        int(value, 16)
    except ValueError as error:
        raise EvidenceError(f"{label} is not hexadecimal") from error
    return value.lower()


def command_options(command: object) -> dict[str, str]:
    require(isinstance(command, list) and all(isinstance(x, str) for x in command),
            "record command is not a string array")
    options: dict[str, str] = {}
    index = 0
    while index < len(command):
        item = command[index]
        if item.startswith("--"):
            require(item not in options, f"duplicate command option {item}")
            if item in {"--force-specialized", "--skip-legacy", "--retain-samples"}:
                options[item] = "true"
            else:
                require(index + 1 < len(command), f"missing value after {item}")
                options[item] = command[index + 1]
                index += 1
        index += 1
    return options


def normalize_case(raw: object) -> dict:
    require(isinstance(raw, list), "case is not an array")
    if len(raw) == 6:
        name, k, r, byte_count, warmup, reuse = raw
        loss, backend, batch = 8, "avx2", 1
    elif len(raw) == 9:
        name, k, r, byte_count, warmup, reuse, loss, backend, batch = raw
    else:
        raise EvidenceError("case must contain six legacy or nine calibrated fields")
    require(isinstance(name, str) and name, "case name is empty")
    integers = (k, r, byte_count, warmup, reuse, loss, batch)
    require(all(isinstance(value, int) and value > 0 for value in integers),
            f"case {name} contains a non-positive integer")
    require(isinstance(backend, str) and backend in BACKENDS,
            f"case {name} has an unsupported backend")
    require(k + r <= 256 and loss <= min(k, r), f"case {name} has invalid GF8 counts")
    padded_side = 1 << (r - 1).bit_length()
    parent_count = 1 << (k + padded_side - 1).bit_length()
    return {
        "name": name, "K": k, "R": r, "shard_bytes": byte_count,
        "warmup": warmup, "reuse": reuse, "loss": loss,
        "backend": backend, "batch": batch, "padded_side": padded_side,
        "parent_count": parent_count,
    }


def finite_positive(value: object, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise EvidenceError(f"{label} is not numeric") from error
    require(math.isfinite(result) and result > 0.0, f"{label} is not positive finite")
    return result


def validate_result(path: Path, case: dict, variant: str, expected_sha: str) -> dict:
    require(path.is_file(), f"missing benchmark result {path}")
    require(sha256_file(path) == expected_sha, f"benchmark result hash changed: {path}")
    document = load_json(path, "benchmark result")
    require(document.get("schema") == BENCHMARK_SCHEMA,
            f"unexpected benchmark schema in {path}")
    parameters = document.get("parameters")
    resolved = document.get("resolved")
    require(isinstance(parameters, dict) and isinstance(resolved, dict),
            f"missing benchmark identity in {path}")
    expected = {
        "K": case["K"], "R": case["R"],
        "requested_profile": "legacy_high_v1", "requested_field": "gf8",
        "requested_backend": case["backend"], "force_generic_decode": False,
        "force_specialized_decode": True, "skip_legacy": True,
        "retain_samples": True, "shard_bytes": case["shard_bytes"],
        "loss_count": case["loss"], "batch": case["batch"],
        "reuse": case["reuse"], "warmup": case["warmup"], "thread_count": 1,
    }
    for key, value in expected.items():
        require(parameters.get(key) == value,
                f"{path}: {key}={parameters.get(key)!r}, expected {value!r}")
    require(resolved.get("profile") == "legacy_high_v1" and
            resolved.get("field") == "gf8" and
            resolved.get("backend") == case["backend"] and
            resolved.get("parent_count") == case["parent_count"] and
            resolved.get("padded_side") == case["padded_side"],
            f"{path}: resolved codec is not the declared high-rate parent shape")
    require(document.get("correctness", {}).get("leopard2_round_trip") is True,
            f"{path}: round trip failed")
    metrics = document.get("metrics", {})
    execution = metrics.get("decode_execution", {})
    setup = metrics.get("decode_plan_setup", {})
    memory = document.get("memory", {})
    digest = document.get("workload_digests")
    require(isinstance(digest, dict), f"{path}: missing workload digests")
    return {
        "variant": variant,
        "decode_us": finite_positive(
            execution.get("median_us_per_batch_call"), f"{path} decode median"),
        "decode_mad_us": finite_positive(
            execution.get("mad_us_per_batch_call", 0.0) + 1e-300,
            f"{path} decode MAD"),
        "setup_us": finite_positive(setup.get("median_us"), f"{path} setup median"),
        "scratch_bytes": int(memory.get("decode_scratch_bytes_per_stripe", -1)),
        "digests": digest,
    }


def clustered_summary(round_logs: list[float]) -> dict:
    require(len(round_logs) >= 2, "at least two independent rounds are required")
    degrees = len(round_logs) - 1
    require(degrees in T95, "unsupported independent-round count")
    mean = statistics.fmean(round_logs)
    deviation = statistics.stdev(round_logs)
    half_width = T95[degrees] * deviation / math.sqrt(len(round_logs))
    lower = math.exp(mean - half_width)
    upper = math.exp(mean + half_width)
    speedup = math.exp(mean)
    return {
        "independent_round_count": len(round_logs),
        "degrees_of_freedom": degrees,
        "independent_round_log_contrasts": round_logs,
        "geometric_control_over_tiled": speedup,
        "improvement_percent": 100.0 * (speedup - 1.0),
        "ci95_lower_percent": 100.0 * (lower - 1.0),
        "ci95_upper_percent": 100.0 * (upper - 1.0),
        "credible_regression_over_2_percent": upper < 0.98,
        "credible_gain_at_least_5_percent": lower >= 1.05,
    }


def analyze(manifest_path: Path, require_rounds: int) -> dict:
    manifest = load_json(manifest_path, "manifest")
    require(manifest.get("schema") == SCHEMA and manifest.get("valid") is True,
            "manifest schema/validity is wrong")
    git_oid(manifest.get("control_commit"), "control commit")
    git_oid(manifest.get("candidate_commit"), "candidate commit")
    identities = manifest.get("identities")
    require(isinstance(identities, dict) and identities, "missing build identities")
    for key, value in identities.items():
        hex_digest(value, f"identity {key}")

    cases = [normalize_case(value) for value in manifest.get("cases", [])]
    require(cases and len({case["name"] for case in cases}) == len(cases),
            "case names are empty or duplicated")
    case_by_name = {case["name"]: case for case in cases}
    orders = manifest.get("round_orders")
    require(isinstance(orders, list) and len(orders) == require_rounds,
            f"manifest must contain exactly {require_rounds} independent rounds")
    for order in orders:
        require(isinstance(order, list) and len(order) == 4 and
                order.count("control") == 2 and order.count("candidate") == 2,
                "each round must contain two control and two candidate invocations")
        require(order in (["control", "candidate", "candidate", "control"],
                          ["candidate", "control", "control", "candidate"]),
                "round is not ABBA or BAAB")

    raw_root = manifest_path.parent / "raw"
    records = manifest.get("records")
    require(isinstance(records, list), "manifest records are missing")
    expected_count = len(cases) * len(orders) * 4
    require(len(records) == expected_count,
            f"manifest has {len(records)} records, expected {expected_count}")
    grouped: dict[str, dict[int, dict[int, dict]]] = {}
    seen = set()
    for record in records:
        require(isinstance(record, dict), "manifest record is not an object")
        name = record.get("case")
        round_index = record.get("round")
        slot = record.get("slot")
        variant = record.get("variant")
        require(name in case_by_name and isinstance(round_index, int) and
                isinstance(slot, int) and 0 <= round_index < len(orders) and 0 <= slot < 4,
                "manifest record coordinate is invalid")
        require(variant == orders[round_index][slot], "record variant/order mismatch")
        coordinate = (name, round_index, slot)
        require(coordinate not in seen, "duplicate manifest record")
        seen.add(coordinate)
        options = command_options(record.get("command"))
        case = case_by_name[name]
        command_expected = {
            "--k": str(case["K"]), "--r": str(case["R"]),
            "--profile": "high", "--field": "gf8", "--backend": case["backend"],
            "--bytes": str(case["shard_bytes"]), "--loss": str(case["loss"]),
            "--batch": str(case["batch"]), "--reuse": str(case["reuse"]),
            "--warmup": str(case["warmup"]), "--threads": "1",
            "--force-specialized": "true", "--skip-legacy": "true",
            "--retain-samples": "true",
        }
        for key, value in command_expected.items():
            require(options.get(key) == value, f"record command {key} mismatch")
        relative = record.get("json")
        require(isinstance(relative, str) and Path(relative).name == relative,
                "record JSON path is not a basename")
        json_path = raw_root / relative
        for suffix in (".stdout", ".stderr"):
            stream_path = json_path.with_suffix(suffix)
            require(stream_path.is_file() and stream_path.stat().st_size == 0,
                    f"benchmark stream is missing or nonempty: {stream_path}")
        result = validate_result(
            json_path, case, variant,
            hex_digest(record.get("json_sha256"), "record result digest"))
        grouped.setdefault(name, {}).setdefault(round_index, {})[slot] = result

    cells = []
    for case in cases:
        round_logs = []
        setup_logs = []
        scratch = {"control": set(), "candidate": set()}
        digest_reference = None
        for round_index, order in enumerate(orders):
            slots = grouped[case["name"]][round_index]
            require(set(slots) == {0, 1, 2, 3}, "incomplete ABBA round")
            for result in slots.values():
                scratch[result["variant"]].add(result["scratch_bytes"])
                if digest_reference is None:
                    digest_reference = result["digests"]
                require(result["digests"] == digest_reference,
                        "control/candidate workload digests differ")
            pair_logs = []
            pair_setup_logs = []
            for left, right in ((0, 1), (2, 3)):
                pair = (slots[left], slots[right])
                control = pair[0] if pair[0]["variant"] == "control" else pair[1]
                candidate = pair[1] if pair[1]["variant"] == "candidate" else pair[0]
                pair_logs.append(math.log(control["decode_us"] / candidate["decode_us"]))
                pair_setup_logs.append(math.log(control["setup_us"] / candidate["setup_us"]))
            round_logs.append(statistics.fmean(pair_logs))
            setup_logs.append(statistics.fmean(pair_setup_logs))
        require(all(len(values) == 1 for values in scratch.values()),
                "scratch size changed within a variant")
        control_scratch = next(iter(scratch["control"]))
        candidate_scratch = next(iter(scratch["candidate"]))
        require(control_scratch > 0 and 0 < candidate_scratch <= control_scratch,
                "candidate scratch is invalid")
        cells.append({
            **case,
            "decode_execution": clustered_summary(round_logs),
            "decode_plan_setup": clustered_summary(setup_logs),
            "scratch": {
                "materialized_bytes": control_scratch,
                "tiled_bytes": candidate_scratch,
                "reduction_percent": 100.0 *
                    (control_scratch - candidate_scratch) / control_scratch,
            },
            "workload_digests": digest_reference,
        })

    summary = {
        "schema": SUMMARY_SCHEMA,
        "source_manifest_sha256": sha256_file(manifest_path),
        "control_commit": manifest["control_commit"],
        "candidate_commit": manifest["candidate_commit"],
        "identities": identities,
        "method": {
            "cpu": manifest.get("cpu"), "sibling": manifest.get("sibling"),
            "lease": manifest.get("lease"), "round_orders": orders,
            "confidence_intervals": "clustered paired-log Student-t 95%",
            "ratio_orientation": "materialized_time_over_tiled_time",
            "validity_is_independent_of_speed": True,
        },
        "cells": cells,
        "status": "pass",
    }
    canonical = json.dumps(summary, sort_keys=True, separators=(",", ":")).encode("ascii")
    summary["content_sha256"] = hashlib.sha256(canonical).hexdigest()
    return summary


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


def self_test() -> None:
    summary = clustered_summary([math.log(1.10), math.log(1.11), math.log(1.09)])
    require(summary["credible_gain_at_least_5_percent"], "gain classification failed")
    regression = clustered_summary([math.log(0.95), math.log(0.96), math.log(0.95)])
    require(regression["credible_regression_over_2_percent"],
            "regression classification failed")
    options = command_options(["bench", "--k", "224", "--force-specialized"])
    require(options == {"--k": "224", "--force-specialized": "true"},
            "command parser self-test failed")
    print("tiled-high evidence analyzer self-test passed")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest")
    parser.add_argument("--output")
    parser.add_argument("--require-rounds", type=int, default=3)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    if not args.manifest or not args.output:
        parser.error("--manifest and --output are required unless --self-test is used")
    summary = analyze(Path(args.manifest), args.require_rounds)
    write_json(Path(args.output), summary)
    print(json.dumps({"cells": len(summary["cells"]),
                      "content_sha256": summary["content_sha256"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except EvidenceError as error:
        print(f"analyze_tiled_high.py: {error}", file=os.sys.stderr)
        raise SystemExit(1)
