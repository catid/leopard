#!/usr/bin/env python3
"""Validate and summarize the pinned K=1/R=1 binding holdout."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import random
import re
import statistics
import tempfile
from typing import Any


FILE_RE = re.compile(
    r"^holdout-b(?P<bytes>[0-9]+)-r(?P<round>[0-9]+)-"
    r"s(?P<sequence>[0-9]+)-(?P<lane>candidate|control|main)\.json$")
BYTE_COUNTS = (64, 1024, 4096, 65536)
LANES = ("candidate", "control", "main")
ROUNDS = tuple(range(1, 10))
EXPECTED_PER_ROUND = 2
BOOTSTRAP_SEED = 0x4B31523142494E44


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while True:
            block = stream.read(1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as stream:
        value = json.load(stream)
    require(isinstance(value, dict), f"{path}: root is not an object")
    return value


def parse_checksum_manifest(artifact: Path) -> dict[str, str]:
    manifest = artifact / "sha256.txt"
    checksums: dict[str, str] = {}
    for line in manifest.read_text(encoding="utf-8").splitlines():
        fields = line.split(maxsplit=1)
        require(len(fields) == 2, "malformed sha256 manifest line")
        name = Path(fields[1].lstrip("* ")).name
        require(name not in checksums, f"duplicate checksum entry: {name}")
        checksums[name] = fields[0]
    expected_names = {
        "candidate", "control", "main", "candidate.text", "control.text"
    }
    require(set(checksums) == expected_names,
            "sha256 manifest membership differs")
    for name, expected in checksums.items():
        require(sha256(artifact / name) == expected,
                f"frozen artifact hash differs: {name}")
    require(checksums["candidate.text"] == checksums["control.text"],
            "candidate/control executable text differs")
    return checksums


def parse_cpu_snapshot(path: Path, prefix: str) -> list[int]:
    fields = path.read_text(encoding="utf-8").strip().split()
    require(len(fields) == 12 and fields[0] == prefix and fields[1] == "cpu20",
            f"malformed {prefix} snapshot")
    return [int(value) for value in fields[2:]]


def validate_cpu_isolation(artifact: Path) -> dict[str, Any]:
    before = parse_cpu_snapshot(artifact / "cpu20-before.txt", "cpu20_before")
    after = parse_cpu_snapshot(artifact / "cpu20-after.txt", "cpu20_after")
    non_idle_indices = (0, 1, 2, 4, 5, 6, 7, 8, 9)
    deltas = [right - left for left, right in zip(before, after)]
    require(all(deltas[index] == 0 for index in non_idle_indices),
            "SMT sibling accumulated non-idle ticks")
    require(deltas[3] > 0, "SMT sibling idle counter did not advance")
    return {
        "timed_cpu": 4,
        "reserved_smt_sibling": 20,
        "before": before,
        "after": after,
        "delta": deltas,
        "sibling_non_idle_delta": 0,
    }


def metric(document: dict[str, Any], lane: str, operation: str) -> float:
    metrics = document["metrics"]
    if operation == "encode":
        summary = metrics["encode_execution"]
    elif lane == "main":
        summary = metrics["decode_including_setup"]
    else:
        summary = metrics["decode_execution"]
    value = summary["median_us_per_batch_call"]
    require(type(value) in (int, float) and math.isfinite(value) and value > 0,
            f"invalid {lane} {operation} process median")
    samples = summary["samples_us_per_batch_call"]
    require(isinstance(samples, list) and len(samples) == 11 and
            all(type(item) in (int, float) and math.isfinite(item) and item > 0
                for item in samples),
            f"invalid {lane} {operation} retained samples")
    return float(value)


def setup_metric(document: dict[str, Any], name: str) -> float:
    value = document["metrics"][name]["median_us"]
    require(type(value) in (int, float) and math.isfinite(value) and value >= 0,
            f"invalid {name} median")
    return float(value)


def geometric_mean(values: list[float]) -> float:
    require(values and all(value > 0 for value in values),
            "geometric mean requires positive values")
    return math.exp(sum(math.log(value) for value in values) / len(values))


def ratio_summary(
    baseline: dict[int, list[float]],
    candidate: dict[int, list[float]],
    repetitions: int,
    seed: int,
) -> dict[str, float]:
    round_ratios: list[float] = []
    for round_number in ROUNDS:
        require(len(baseline[round_number]) == EXPECTED_PER_ROUND and
                len(candidate[round_number]) == EXPECTED_PER_ROUND,
                "round pairing is incomplete")
        round_ratios.append(
            geometric_mean(baseline[round_number]) /
            geometric_mean(candidate[round_number]))
    center = geometric_mean(round_ratios)
    rng = random.Random(seed)
    samples = []
    for _ in range(repetitions):
        samples.append(geometric_mean([
            round_ratios[rng.randrange(len(round_ratios))]
            for _ in round_ratios
        ]))
    samples.sort()
    low = samples[int(0.025 * (repetitions - 1))]
    high = samples[int(0.975 * (repetitions - 1))]
    return {
        "geometric_mean_ratio": round(center, 9),
        "bootstrap_95_low": round(low, 9),
        "bootstrap_95_high": round(high, 9),
    }


def validate_document(
    document: dict[str, Any],
    lane: str,
    byte_count: int,
    expected_source_commit: str,
    expected_source_tree: str,
    expected_main_commit: str,
) -> None:
    parameters = document["parameters"]
    expected_reuse = 500000 if byte_count in (64, 1024) else (
        100000 if byte_count == 4096 else 10000)
    common = {
        "K": 1,
        "R": 1,
        "shard_bytes": byte_count,
        "loss_count": 1,
        "missing_original_indices": [0],
        "batch": 1,
        "reuse": expected_reuse,
        "iterations": 11,
        "warmup": 100,
        "thread_count": 1,
        "seed": 1,
    }
    for name, expected in common.items():
        require(parameters.get(name) == expected,
                f"{lane}/{byte_count}: parameter {name} differs")
    if lane == "main":
        require(document.get("schema") == "leopard-main-benchmark-v1",
                "main schema differs")
        require(document["correctness"].get("round_trip") is True,
                "main round trip failed")
        require(document["build"].get("main_source_commit") ==
                expected_main_commit, "main source commit differs")
    else:
        require(document.get("schema") == "leopard2-benchmark-v5",
                f"{lane} schema differs")
        require(parameters.get("requested_profile") == "legacy_high_v1" and
                parameters.get("requested_field") == "gf8" and
                parameters.get("requested_backend") == "avx2" and
                parameters.get("skip_legacy") is True and
                parameters.get("retain_samples") is True and
                parameters.get("attest_source") is True,
                f"{lane} public execution request differs")
        require(document["resolved"].get("profile") == "legacy_high_v1" and
                document["resolved"].get("field") == "gf8" and
                document["resolved"].get("backend") == "avx2" and
                document["resolved"].get("thread_count") == 1,
                f"{lane} resolved identity differs")
        require(document["correctness"].get("leopard2_round_trip") is True,
                f"{lane} round trip failed")
        build = document["build"]
        require(build.get("prevalidated_batch_experiment") is True and
                build.get("source_commit") == expected_source_commit and
                build.get("source_tree") == expected_source_tree and
                build.get("source_tracked_dirty") is False,
                f"{lane} build/source identity differs")
        setup_metric(document, "encode_binding_setup")
        setup_metric(document, "decode_binding_setup")
    metric(document, lane, "encode")
    metric(document, lane, "decode")


def atomic_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=str(path.parent))
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(value, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-source-commit", required=True)
    parser.add_argument("--expected-source-tree", required=True)
    parser.add_argument("--expected-main-commit", required=True)
    parser.add_argument("--bootstrap-repetitions", type=int, default=100000)
    args = parser.parse_args()
    require(args.bootstrap_repetitions >= 1000,
            "bootstrap repetitions must be at least 1000")

    artifact = args.artifact.resolve()
    checksums = parse_checksum_manifest(artifact)
    isolation = validate_cpu_isolation(artifact)
    records: dict[int, dict[str, dict[int, list[dict[str, Any]]]]] = {
        byte_count: {
            lane: {round_number: [] for round_number in ROUNDS}
            for lane in LANES
        } for byte_count in BYTE_COUNTS
    }
    sequences: dict[int, set[int]] = {value: set() for value in BYTE_COUNTS}
    paths = sorted(artifact.glob("holdout-b*-r*-s*-*.json"))
    require(len(paths) == 216, f"expected 216 raw files, found {len(paths)}")
    for path in paths:
        match = FILE_RE.fullmatch(path.name)
        require(match is not None, f"unexpected raw filename: {path.name}")
        byte_count = int(match.group("bytes"))
        round_number = int(match.group("round"))
        sequence = int(match.group("sequence"))
        lane = match.group("lane")
        require(byte_count in records and round_number in ROUNDS,
                f"raw file domain differs: {path.name}")
        require(sequence not in sequences[byte_count],
                f"duplicate sequence for {byte_count}: {sequence}")
        sequences[byte_count].add(sequence)
        document = load_json(path)
        validate_document(document, lane, byte_count,
                          args.expected_source_commit,
                          args.expected_source_tree,
                          args.expected_main_commit)
        records[byte_count][lane][round_number].append(document)

    for byte_count in BYTE_COUNTS:
        require(sequences[byte_count] == set(range(1, 55)),
                f"sequence domain differs for {byte_count}")
        expected_digest = None
        for lane in LANES:
            for round_number in ROUNDS:
                documents = records[byte_count][lane][round_number]
                require(len(documents) == EXPECTED_PER_ROUND,
                        f"{byte_count}/{lane}/{round_number}: count differs")
                for document in documents:
                    digest = document.get("workload_digests")
                    require(isinstance(digest, dict) and
                            digest.get("algorithm") == "fnv1a64",
                            "workload digest is missing")
                    if expected_digest is None:
                        expected_digest = digest
                    require(digest == expected_digest,
                            f"workload digest differs for {byte_count}")

    summaries = []
    for byte_index, byte_count in enumerate(BYTE_COUNTS):
        operations: dict[str, Any] = {}
        for operation_index, operation in enumerate(("encode", "decode")):
            values: dict[str, dict[int, list[float]]] = {}
            for lane in LANES:
                values[lane] = {
                    round_number: [
                        metric(document, lane, operation)
                        for document in records[byte_count][lane][round_number]
                    ] for round_number in ROUNDS
                }
            process_medians = {
                lane: round(statistics.median([
                    value for round_values in values[lane].values()
                    for value in round_values
                ]), 9) for lane in LANES
            }
            seed_base = BOOTSTRAP_SEED + byte_index * 17 + operation_index
            operations[operation] = {
                "process_median_us": process_medians,
                "main_over_candidate": ratio_summary(
                    values["main"], values["candidate"],
                    args.bootstrap_repetitions, seed_base),
                "control_over_candidate": ratio_summary(
                    values["control"], values["candidate"],
                    args.bootstrap_repetitions, seed_base + 1000),
            }

        candidate_documents = [
            document
            for round_values in records[byte_count]["candidate"].values()
            for document in round_values
        ]
        encode_setup = statistics.median([
            setup_metric(document, "encode_binding_setup")
            for document in candidate_documents
        ])
        decode_setup = statistics.median([
            setup_metric(document, "decode_binding_setup")
            for document in candidate_documents
        ])
        for operation, setup in (("encode", encode_setup),
                                 ("decode", decode_setup)):
            medians = operations[operation]["process_median_us"]
            saved = medians["main"] - medians["candidate"]
            operations[operation]["binding_setup_median_us"] = round(setup, 9)
            operations[operation]["main_break_even_calls"] = (
                math.ceil(setup / saved) if saved > 0 else None)
        summaries.append({"shard_bytes": byte_count, **operations})

    result = {
        "schema": "leopard2-k1-binding-holdout-v1",
        "source": {
            "commit": args.expected_source_commit,
            "tree": args.expected_source_tree,
            "tracked_dirty": False,
            "main_commit": args.expected_main_commit,
        },
        "frozen_sha256": checksums,
        "isolation": isolation,
        "campaign": {
            "processes": 216,
            "rounds": 9,
            "processes_per_lane_per_cell": 18,
            "samples_per_process": 11,
            "warmups_per_process": 100,
            "bootstrap_repetitions": args.bootstrap_repetitions,
            "bootstrap_unit": "counterbalanced round",
            "correctness_all": True,
            "workload_identity_count": len(BYTE_COUNTS),
        },
        "summaries": summaries,
    }
    atomic_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
