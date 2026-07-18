#!/usr/bin/env python3
"""Validate the private Algorithm 5 copy/no-copy attribution benchmark."""

from __future__ import annotations

import argparse
import copy
import json
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Any


HEX64 = re.compile(r"^[0-9a-f]{16}$")


class ContractError(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def validate_document(
    document: object,
    *,
    mode: str,
    workspace: str,
    field: str,
    tail_bytes: int,
) -> dict[str, Any]:
    require(isinstance(document, dict), "benchmark result is not an object")
    result = document
    require(result.get("schema") == "leopard2-benchmark-v4",
            "attribution benchmark schema changed")
    parameters = result.get("parameters")
    resolved = result.get("resolved")
    attribution = result.get("high_evaluator_attribution")
    correctness = result.get("correctness")
    digests = result.get("workload_digests")
    require(isinstance(parameters, dict) and isinstance(resolved, dict) and
            isinstance(attribution, dict) and isinstance(correctness, dict) and
            isinstance(digests, dict), "benchmark attestation is incomplete")
    require(parameters.get("high_evaluator_mode") == mode and
            parameters.get("requested_profile") == "legacy_high_v1" and
            parameters.get("requested_field") == field and
            parameters.get("force_specialized_decode") is True and
            parameters.get("force_generic_decode") is False and
            parameters.get("force_materialized_decode") is
                (workspace == "materialized") and
            parameters.get("force_tiled_decode") is (workspace == "tiled") and
            parameters.get("skip_legacy") is True and
            parameters.get("retain_samples") is True and
            parameters.get("report_decode_path") is True,
            "benchmark request attestation differs from the signed role")
    require(resolved.get("profile") == "legacy_high_v1" and
            resolved.get("field") == field and
            resolved.get("selected_decode_path") == workspace and
            resolved.get("selected_decode_rule") == "forced_" + workspace and
            resolved.get("high_evaluator_mode") == mode and
            resolved.get("decode_tail_bytes") == tail_bytes,
            "resolved field/path/tail/mode attestation differs")
    require(correctness.get("leopard2_round_trip") is True and
            correctness.get("legacy_comparison") is None,
            "benchmark did not attest a private round trip")
    require(digests.get("algorithm") == "fnv1a64" and
            all(isinstance(digests.get(name), str) and
                HEX64.fullmatch(digests[name]) is not None
                for name in ("original_data", "transmitted_parity",
                             "recovered_originals")),
            "workload digest attestation is malformed")
    output_blocks = attribution.get("output_blocks")
    two_way = attribution.get("fft_butterfly2_out_of_place")
    four_way = attribution.get("fft_butterfly4_out_of_place")
    fallbacks = attribution.get("compatibility_copy_fallbacks")
    require(attribution.get("mode") == mode and
            attribution.get("invariant_passed") is True and
            all(type(value) is int and value >= 0 for value in
                (output_blocks, two_way, four_way, fallbacks)) and
            output_blocks > 0,
            "high-evaluator counter attestation is malformed")
    if mode == "no-copy":
        require(fallbacks == 0 and two_way + four_way > 0,
                "no-copy role entered a copy fallback or no out-of-place layer")
    elif mode == "copy-fallback":
        require(fallbacks == output_blocks and two_way + four_way == 0,
                "copy role did not force exactly one fallback per output block")
    else:
        raise ContractError("unknown signed high-evaluator mode")
    return result


def validate_pair(
    no_copy: object,
    copy_fallback: object,
    *,
    workspace: str,
    field: str,
    tail_bytes: int,
) -> None:
    first = validate_document(
        no_copy, mode="no-copy", workspace=workspace, field=field,
        tail_bytes=tail_bytes)
    second = validate_document(
        copy_fallback, mode="copy-fallback", workspace=workspace, field=field,
        tail_bytes=tail_bytes)
    require(first["workload_digests"] == second["workload_digests"] and
            first["parameters"]["missing_original_indices"] ==
                second["parameters"]["missing_original_indices"] and
            first["memory"] == second["memory"],
            "copy/no-copy roles do not describe one identical workload")
    first_parameters = dict(first["parameters"])
    second_parameters = dict(second["parameters"])
    first_parameters.pop("high_evaluator_mode")
    second_parameters.pop("high_evaluator_mode")
    require(first_parameters == second_parameters,
            "copy/no-copy parameters differ outside the mode selector")
    first_resolved = dict(first["resolved"])
    second_resolved = dict(second["resolved"])
    first_resolved.pop("high_evaluator_mode")
    second_resolved.pop("high_evaluator_mode")
    require(first_resolved == second_resolved,
            "copy/no-copy resolution differs outside the mode selector")


def benchmark_command(
    executable: Path,
    output: Path,
    *,
    k: int,
    r: int,
    shard_bytes: int,
    losses: int,
    seed: int,
    field: str,
    workspace: str,
    mode: str,
) -> list[str]:
    return [
        str(executable), "--k", str(k), "--r", str(r),
        "--profile", "high", "--field", field, "--backend", "scalar",
        "--bytes", str(shard_bytes), "--loss", str(losses), "--batch", "1",
        "--reuse", "1", "--iterations", "1", "--warmup", "0",
        "--threads", "1", "--seed", str(seed), "--force-specialized",
        "--force-" + workspace, "--skip-legacy", "--retain-samples",
        "--report-decode-path", "--high-evaluator-mode", mode,
        "--json", str(output),
    ]


def run_checked(arguments: list[str], *, expect_success: bool = True) -> None:
    completed = subprocess.run(
        arguments, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, timeout=30, check=False,
        env={
            "LANG": "C", "LC_ALL": "C", "OMP_DYNAMIC": "FALSE",
            "OMP_NUM_THREADS": "1", "OMP_PLACES": "cores",
            "OMP_PROC_BIND": "TRUE", "PATH": "/usr/bin:/bin", "TZ": "UTC",
        })
    if expect_success:
        require(completed.returncode == 0 and not completed.stdout and
                not completed.stderr,
                "attribution benchmark failed or emitted an unexpected stream")
    else:
        require(completed.returncode != 0,
                "attribution benchmark accepted a malformed selector contract")


def smoke(executable: Path) -> None:
    executable = executable.resolve(strict=True)
    cases = (
        ("gf8", "materialized", 192, 64, 64, 4, 101, 0),
        ("gf8", "materialized", 192, 64, 65, 4, 103, 1),
        ("gf8", "tiled", 240, 16, 256, 8, 107, 0),
        ("gf8", "tiled", 240, 16, 257, 8, 109, 1),
        ("gf16", "materialized", 257, 63, 64, 4, 113, 0),
        ("gf16", "materialized", 257, 63, 66, 4, 127, 2),
        ("gf16", "tiled", 1000, 200, 256, 8, 131, 0),
        ("gf16", "tiled", 1000, 200, 258, 8, 137, 2),
    )
    with tempfile.TemporaryDirectory(prefix="leo2-high-copy-contract-") as root_text:
        root = Path(root_text)
        for index, (field, workspace, k, r, byte_count, losses, seed,
                    tail_bytes) in enumerate(cases):
            documents: dict[str, object] = {}
            for mode in ("no-copy", "copy-fallback"):
                output = root / f"{index}-{mode}.json"
                run_checked(benchmark_command(
                    executable, output, k=k, r=r, shard_bytes=byte_count,
                    losses=losses, seed=seed, field=field,
                    workspace=workspace, mode=mode))
                documents[mode] = json.loads(output.read_text(encoding="utf-8"))
            validate_pair(
                documents["no-copy"], documents["copy-fallback"],
                workspace=workspace, field=field, tail_bytes=tail_bytes)

        base = benchmark_command(
            executable, root / "invalid.json", k=192, r=64, shard_bytes=64,
            losses=4, seed=149, field="gf8", workspace="materialized",
            mode="no-copy")
        missing_mode = base[:]
        mode_index = missing_mode.index("--high-evaluator-mode")
        del missing_mode[mode_index:mode_index + 2]
        run_checked(missing_mode, expect_success=False)
        invalid_mode = base[:]
        invalid_mode[invalid_mode.index("no-copy")] = "unknown"
        run_checked(invalid_mode, expect_success=False)
        missing_path = [item for item in base if item != "--force-materialized"]
        run_checked(missing_path, expect_success=False)
    print("high-decode copy attribution benchmark contract passed: 8 paired cells")


def synthetic_document(mode: str) -> dict[str, Any]:
    copy_mode = mode == "copy-fallback"
    return {
        "schema": "leopard2-benchmark-v4",
        "parameters": {
            "requested_profile": "legacy_high_v1", "requested_field": "gf8",
            "force_specialized_decode": True, "force_generic_decode": False,
            "force_materialized_decode": True, "force_tiled_decode": False,
            "skip_legacy": True, "retain_samples": True,
            "report_decode_path": True, "high_evaluator_mode": mode,
            "missing_original_indices": [7],
        },
        "resolved": {
            "profile": "legacy_high_v1", "field": "gf8",
            "selected_decode_path": "materialized",
            "selected_decode_rule": "forced_materialized",
            "decode_tail_bytes": 1, "high_evaluator_mode": mode,
        },
        "high_evaluator_attribution": {
            "mode": mode, "output_blocks": 2,
            "fft_butterfly2_out_of_place": 0 if copy_mode else 16,
            "fft_butterfly4_out_of_place": 0,
            "compatibility_copy_fallbacks": 2 if copy_mode else 0,
            "invariant_passed": True,
        },
        "correctness": {"leopard2_round_trip": True,
                        "legacy_comparison": None},
        "workload_digests": {
            "algorithm": "fnv1a64", "original_data": "0" * 16,
            "transmitted_parity": "1" * 16,
            "recovered_originals": "2" * 16,
        },
        "memory": {"decode_scratch_bytes_per_stripe": 64},
    }


def self_test() -> None:
    no_copy = synthetic_document("no-copy")
    fallback = synthetic_document("copy-fallback")
    validate_pair(no_copy, fallback, workspace="materialized", field="gf8",
                  tail_bytes=1)
    mutations = []
    wrong_mode = copy.deepcopy(no_copy)
    wrong_mode["resolved"]["high_evaluator_mode"] = "copy-fallback"
    mutations.append(wrong_mode)
    wrong_counter = copy.deepcopy(no_copy)
    wrong_counter["high_evaluator_attribution"]["compatibility_copy_fallbacks"] = 1
    mutations.append(wrong_counter)
    wrong_path = copy.deepcopy(no_copy)
    wrong_path["resolved"]["selected_decode_path"] = "tiled"
    mutations.append(wrong_path)
    wrong_digest = copy.deepcopy(no_copy)
    wrong_digest["workload_digests"]["recovered_originals"] = "g" * 16
    mutations.append(wrong_digest)
    for mutation in mutations:
        try:
            validate_pair(mutation, fallback, workspace="materialized", field="gf8",
                          tail_bytes=1)
        except ContractError:
            continue
        raise ContractError("adversarial mode/counter/path/digest mutation passed")
    print("high-decode copy attribution contract self-test passed: 4 mutations rejected")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("self-test")
    smoke_parser = subparsers.add_parser("smoke")
    smoke_parser.add_argument("--benchmark", required=True, type=Path)
    options = parser.parse_args()
    if options.command == "self-test":
        self_test()
    else:
        smoke(options.benchmark)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ContractError, OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"benchmark_contract.py: {error}", file=__import__("sys").stderr)
        raise SystemExit(1)
