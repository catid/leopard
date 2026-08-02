#!/usr/bin/env python3
"""Focused fail-closed checks for the T=32 final-IFFT ABBA specialization."""

from __future__ import annotations

import importlib.util
import math
import os
import sys
import tempfile
from pathlib import Path
from typing import Any, Callable


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_runner() -> Any:
    path = Path(__file__).resolve().with_name(
        "run_t32_final_ifft_abba.py")
    specification = importlib.util.spec_from_file_location(
        "leopard2_t32_final_ifft_runner_test_target", path)
    require(specification is not None and specification.loader is not None,
            "cannot load T=32 runner")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


RUNNER = load_runner()
BASE = RUNNER.BASE


def expect_failure(action: Callable[[], object], message: str) -> None:
    try:
        action()
    except Exception:
        return
    raise RuntimeError(message)


def normalized(
    implementation: str,
    encode_us: float,
    digests: dict[str, str],
) -> dict[str, Any]:
    return {
        "implementation": implementation,
        "normalized": {
            "encode_us": encode_us,
            "digests": dict(digests),
        },
    }


def synthetic_round(order: tuple[str, ...]) -> dict[str, Any]:
    digests = {
        "original_data": "1" * 16,
        "transmitted_parity": "2" * 16,
        "recovered_originals": "1" * 16,
    }
    times = {"candidate": 1.0, "control": 1.25, "main": 1.5}
    return {
        "order": list(order),
        "invocations": [
            normalized(label, times[label], digests) for label in order
        ],
        "isolation": {"accepted": True},
    }


def benchmark_result(
    cell: dict[str, Any],
    source_commit: str,
    source_tree: str,
) -> dict[str, Any]:
    return {
        "schema": "leopard2-benchmark-v5",
        "parameters": {
            "K": cell["K"],
            "R": cell["R"],
            "shard_bytes": cell["bytes"],
            "loss_count": 1,
            "batch": 1,
            "reuse": cell["reuse"],
            "iterations": 15,
            "warmup": 64,
            "thread_count": 1,
            "seed": cell["seed"],
        },
        "resolved": {
            "profile": "legacy_high_v1",
            "field": "gf8",
            "backend": "avx2",
            "thread_count": 1,
        },
        "correctness": {"leopard2_round_trip": True},
        "workload_digests": {
            "algorithm": "fnv1a64",
            "original_data": "1" * 16,
            "transmitted_parity": "2" * 16,
            "recovered_originals": "1" * 16,
        },
        "build": {
            "source_commit": source_commit,
            "source_tree": source_tree,
            "source_tracked_dirty": False,
        },
        "metrics": {
            "encode_execution": {"median_us_per_batch_call": 1.0},
        },
    }


def main() -> int:
    cells = RUNNER.campaign_cells()
    require(len(cells) == 72, "T=32 matrix size changed")
    require(
        {(cell["K"], cell["R"], cell["bytes"]) for cell in cells} == {
            (k, r, shard_bytes)
            for k in RUNNER.ORIGINAL_COUNTS
            for r in RUNNER.RECOVERY_COUNTS
            for shard_bytes in RUNNER.BYTE_COUNTS
        },
        "T=32 Cartesian acceptance matrix changed")
    primary = [cell for cell in cells if cell["primary_target"]]
    require(
        {(cell["K"], cell["R"], cell["bytes"]) for cell in primary} ==
            set(RUNNER.PRIMARY_TARGETS) and
        all(cell["role"] == "primary_target" for cell in primary),
        "T=32 primary targets changed")
    require(
        sum(cell["candidate_selected"] for cell in cells) == 24 and
        all(cell["candidate_selected"] == (cell["bytes"] == 64)
            for cell in cells) and
        all(cell["control_selected"] is False for cell in cells),
        "T=32 selector annotations changed")
    require(len({cell["id"] for cell in cells}) == 72 and
            len({cell["seed"] for cell in cells}) == 72,
            "T=32 IDs or seeds are not unique")

    require(BASE.BENCHMARK_CPU == 14 and BASE.RESERVED_SIBLING == 30,
            "authoritative CPU pair changed")
    require(BASE.MODE_SYMBOL ==
            "_ZN7leopard3ff8L29g_high_final_ifft2_range_modeE",
            "T=32 selector symbol changed")
    require(BASE.SCHEMA == "leopard2-gf8-t32-final-ifft-abba/v1" and
            BASE.SUMMARY_SCHEMA ==
                "leopard2-gf8-t32-final-ifft-summary/v1",
            "T=32 evidence schema changed")
    runner_path = Path(RUNNER.__file__).resolve()
    require(BASE.RUNNER_PATH == runner_path and
            os.access(runner_path, os.X_OK),
            "T=32 runner is not its executable provenance root")
    expected_dependencies = (
        runner_path,
        Path(BASE.__file__).resolve(),
        Path(BASE.T8_SUPPORT.__file__).resolve(),
        Path(BASE.MAIN_SUPPORT.__file__).resolve(),
        Path(BASE.MAIN_SUPPORT.git_capture.__file__).resolve(),
        Path(BASE.MAIN_SUPPORT.git_capture._build_provenance.__file__).resolve(),
        Path(BASE.MAIN_SUPPORT.link_common.__file__).resolve(),
    )
    require(BASE.RUNNER_DEPENDENCIES == expected_dependencies and
            len(set(expected_dependencies)) == 7,
            "T=32 runner dependency chain is incomplete")
    actual_dependency_identities = BASE.support_file_identities(
        BASE.RUNNER_DEPENDENCIES)
    require(len(actual_dependency_identities) == 7 and
            BASE.verify_support_file_identities(
                actual_dependency_identities) == actual_dependency_identities,
            "T=32 runner dependency preflight failed")
    require(BASE.REQUIRE_NORMALIZED_FULL_FILE_EQUIVALENCE is True and
            BASE.EXPECTED_BINARY_SHA256 == {
                "candidate":
                    "6737262ea4b206690678c01a742c6e8e56f99cfdec2ead968f0c26aa9db69d42",
                "control":
                    "b88bd5c7c8a57915ebb4cfd09e74e39ee9dcd1691f7fea1d4d16d614be691a64",
                "main":
                    "a43d7f43ff2e887ebcd47a1e94f806847a5d8b858a4e383e6c8d5e528a7dd910",
            },
            "T=32 frozen binary policy changed")

    cell = primary[0]
    candidate = BASE.benchmark_command(
        "candidate", Path("/frozen/candidate"), cell,
        BASE.BENCHMARK_CPU, 15, 64)
    control = BASE.benchmark_command(
        "control", Path("/frozen/control"), cell,
        BASE.BENCHMARK_CPU, 15, 64)
    exact_main = BASE.benchmark_command(
        "main", Path("/frozen/main"), cell,
        BASE.BENCHMARK_CPU, 15, 64)
    require(candidate[0:4] == [
                "/usr/bin/prlimit", f"--as={BASE.ADDRESS_SPACE_BYTES}",
                "/usr/bin/taskset", "-c"] and
            str(BASE.BENCHMARK_CPU) in candidate and
            "--profile" in candidate and "high" in candidate and
            "--field" in candidate and "gf8" in candidate and
            "--backend" in candidate and "avx2" in candidate and
            "--attest-source" in candidate,
            "candidate command lost pinning or source attestation")
    require("--attest-source" in control and
            "--profile" not in exact_main and
            "--attest-source" not in exact_main,
            "control or exact-main command identity changed")

    rounds = [synthetic_round(order) for order in BASE.ROUND_ORDERS]
    analysis = BASE.analyze(cell, rounds)
    require(math.isclose(
                analysis["candidate_vs_control"]["speedup"], 1.25) and
            math.isclose(
                analysis["candidate_vs_main"]["speedup"], 1.5),
            "T=32 ABBA contrast definition changed")
    bad_rounds = [synthetic_round(order) for order in BASE.ROUND_ORDERS]
    bad_rounds[1]["invocations"][0]["normalized"]["digests"][
        "transmitted_parity"] = "f" * 16
    expect_failure(lambda: BASE.analyze(cell, bad_rounds),
                   "mismatched process digests were accepted")

    source_commit = "a" * 40
    source_tree = "b" * 40
    valid_result = benchmark_result(cell, source_commit, source_tree)
    BASE.validate_result(
        "candidate", valid_result, cell, source_commit, source_tree, 15, 64)
    invalid_result = benchmark_result(cell, source_commit, source_tree)
    invalid_result["workload_digests"]["original_data"] = "A" * 16
    expect_failure(
        lambda: BASE.validate_result(
            "candidate", invalid_result, cell, source_commit, source_tree,
            15, 64),
        "uppercase workload digest was accepted")

    retained: list[dict[str, Any]] = []

    def fail_third(label: str) -> dict[str, Any]:
        if label == "control":
            raise BASE.EvidenceError("synthetic third-child failure")
        return {"implementation": label}

    expect_failure(
        lambda: BASE.append_invocations(
            ("main", "candidate", "control"), retained, fail_third),
        "synthetic child failure was accepted")
    require(retained == [
                {"implementation": "main"},
                {"implementation": "candidate"},
            ],
            "successful children were not retained before a later failure")

    with tempfile.TemporaryDirectory(
            prefix="leopard-t32-evidence-self-test-") as directory:
        root = Path(directory)
        dependency_paths = tuple(root / f"support-{index}.py"
                                 for index in range(4))
        for index, path in enumerate(dependency_paths):
            path.write_text(f"support {index}\n", encoding="ascii")
        dependency_identities = BASE.support_file_identities(dependency_paths)
        require(BASE.verify_support_file_identities(dependency_identities) ==
                dependency_identities,
                "stable dependency identities were rejected")
        dependency_paths[2].write_text(
            "mutated support dependency\n", encoding="ascii")
        expect_failure(
            lambda: BASE.verify_support_file_identities(
                dependency_identities),
            "mutated runner dependency was accepted")

        candidate_path = root / "candidate.elf"
        control_path = root / "control.elf"
        candidate_payload = bytearray(b"x" * 128)
        control_payload = bytearray(candidate_payload)
        candidate_payload[16:20] = (1).to_bytes(4, "little")
        control_payload[16:20] = (2).to_bytes(4, "little")
        candidate_payload[64:84] = b"c" * 20
        control_payload[64:84] = b"d" * 20
        candidate_path.write_bytes(candidate_payload)
        control_path.write_bytes(control_payload)
        allowed_ranges = (
            {"name": "diagnostic_selector", "file_offset": 16, "size": 4},
            {"name": ".note.gnu.build-id", "file_offset": 64, "size": 20},
        )
        full_file = BASE.normalized_full_file_comparison(
            candidate_path, control_path, allowed_ranges)
        require(full_file["difference_count"] == 21 and
                full_file["difference_ranges"] == [
                    {"file_offset": 16, "size": 1},
                    {"file_offset": 64, "size": 20},
                ] and
                len(full_file["normalized_sha256"]) == 64,
                "normalized full-file evidence changed")
        control_payload[40] ^= 1
        control_path.write_bytes(control_payload)
        expect_failure(
            lambda: BASE.normalized_full_file_comparison(
                candidate_path, control_path, allowed_ranges),
            "unauthorized candidate/control data-byte difference was accepted")

        executable = root / "candidate"
        executable.write_bytes(b"synthetic executable")
        failure_output = root / "failure-output"
        failure_output.mkdir()
        expected_hash = "e" * 64
        identity = {"path": str(executable), "sha256": expected_hash}
        original_sha256 = BASE.sha256
        original_run = BASE.subprocess.run
        completed_type = BASE.subprocess.CompletedProcess
        hash_values = iter((expected_hash, "f" * 64))
        try:
            BASE.sha256 = lambda unused_path: next(hash_values)
            BASE.subprocess.run = lambda *args, **kwargs: completed_type(
                args[0] if args else [], 0,
                stdout=b"captured stdout", stderr=b"captured stderr")
            expect_failure(
                lambda: BASE.run_one(
                    "candidate", identity, cell, BASE.BENCHMARK_CPU,
                    source_commit, source_tree, 15, 64, failure_output),
                "post-child binary drift was accepted")
        finally:
            BASE.sha256 = original_sha256
            BASE.subprocess.run = original_run
        stdout_files = list(failure_output.glob("*.stdout"))
        stderr_files = list(failure_output.glob("*.stderr"))
        require(len(stdout_files) == 1 and len(stderr_files) == 1 and
                stdout_files[0].read_bytes() == b"captured stdout" and
                stderr_files[0].read_bytes() == b"captured stderr",
                "post-child binary drift lost captured process output")

    print("T32 final-IFFT ABBA runner self-test passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
