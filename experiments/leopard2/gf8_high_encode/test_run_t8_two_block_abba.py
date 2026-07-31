#!/usr/bin/env python3
"""Fail-closed unit tests for streamed ragged-T8 evidence checkpoints."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import sys
import tempfile
from pathlib import Path
from typing import Any, Callable


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_runner() -> Any:
    path = Path(__file__).resolve().with_name("run_t8_two_block_abba.py")
    specification = importlib.util.spec_from_file_location(
        "leopard2_t8_ragged_runner_test_target", path)
    require(specification is not None and specification.loader is not None,
            "cannot load ragged-T8 runner")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


RUNNER = load_runner()


def expect_failure(action: Callable[[], object], message: str) -> None:
    try:
        action()
    except Exception:
        return
    raise RuntimeError(message)


def invocation(
    implementation: str,
    encode_us: float,
    digests: dict[str, str],
) -> dict[str, Any]:
    return {
        "implementation": implementation,
        "elapsed_ns": 1000,
        "stdout_sha256": hashlib.sha256(
            f"{implementation}-{encode_us}".encode("ascii")).hexdigest(),
        "stderr_sha256": hashlib.sha256(b"").hexdigest(),
        "normalized": {
            "encode_us": encode_us,
            "digests": dict(digests),
        },
        "result": {
            "schema": "synthetic-runner-test/v1",
            "implementation": implementation,
            "encode_us": encode_us,
        },
    }


def build_record(
    root: Path,
    cell: dict[str, Any],
    manifest_sha256: str,
) -> tuple[dict[str, Any], Path]:
    pair_lease = {
        "payload": {
            "schema": "leopard2-cpu-pair-lease/v1",
            "cpus": [4, 20],
            "uid": 1000,
        },
        "sha256": "1" * 64,
    }
    presample = {
        "pair_lease": pair_lease,
        "delta": {
            "benchmark_cpu": {"nonidle_jiffies": 0},
            "reserved_sibling": {"nonidle_jiffies": 0},
        },
    }
    digests = {
        "original_data": "1" * 16,
        "transmitted_parity": "2" * 16,
        "recovered_originals": "1" * 16,
    }
    rounds = []
    first_artifact: Path | None = None
    for round_index, order in enumerate(RUNNER.TARGET_ORDER):
        compact = []
        for slot_index, implementation in enumerate(order):
            value = invocation(
                implementation,
                {
                    "candidate": 1.0,
                    "control": 1.2,
                    "main": 1.3,
                }[implementation],
                digests)
            artifact = root / (
                f"partial-{cell['id']}-round-{round_index}"
                f"-slot-{slot_index}.json")
            RUNNER.write_exclusive(artifact, value)
            if first_artifact is None:
                first_artifact = artifact
            compact.append(RUNNER.compact_invocation(value, artifact))
        rounds.append({
            "round": round_index,
            "order": list(order),
            "attempt": 0,
            "invocations": compact,
            "isolation": {
                "accepted": True,
                "pair_lease": pair_lease,
                "delta": {
                    "benchmark_cpu": {"nonidle_jiffies": 1},
                    "reserved_sibling": {"nonidle_jiffies": 0},
                },
            },
        })
    require(first_artifact is not None, "synthetic invocation set is empty")
    return {
        "schema": RUNNER.RAGGED_CELL_SCHEMA,
        "manifest_sha256": manifest_sha256,
        "completed_utc": "2026-07-31T00:00:00Z",
        "run_nonce": "synthetic",
        "pair_lease": pair_lease,
        "presample": presample,
        "cell": dict(cell),
        "rounds": rounds,
        "rejected_isolation_attempts": [],
        "analysis": RUNNER.analyze(cell, rounds),
    }, first_artifact


def main() -> int:
    cells = RUNNER.ragged_extension_cells()
    require(len(cells) == 2058, "ragged matrix size changed")
    cell = dict(cells[0])
    require(cell["role"] == "target", "first ragged cell is not a target")
    cells_by_id = {item["id"]: item for item in cells}
    excluded = {
        "ragged-target-k5-r5-b191",
        "ragged-target-k6-r5-b319",
        "ragged-target-k6-r6-b319",
        "ragged-target-k7-r5-b319",
        "ragged-target-k7-r6-b319",
    }
    require(sum(item["role"] == "target" for item in cells) == 1213,
            "ragged final target count changed")
    require(sum(item["role"] == "neighbor" for item in cells) == 845,
            "ragged final neighbor count changed")
    require(all(
            cells_by_id[cell_id]["role"] == "neighbor" and
            cells_by_id[cell_id]["candidate_selected"] is False and
            cells_by_id[cell_id]["control_selected"] is False
            for cell_id in excluded),
        "ragged final selector lost a holdout exclusion")
    require(
        cells_by_id["ragged-target-k6-r5-b191"]["role"] == "target" and
        cells_by_id["ragged-target-k5-r5-b319"]["role"] == "target" and
        cells_by_id["ragged-target-k7-r7-b319"]["role"] == "target",
        "ragged final selector widened a holdout exclusion")
    manifest_sha256 = "a" * 64
    with tempfile.TemporaryDirectory(
            prefix="leopard2-t8-ragged-runner-test-") as directory:
        root = Path(directory)
        atomic = root / "atomic.json"
        RUNNER.write_atomic_exclusive(atomic, {"value": 1})
        require(json.loads(atomic.read_text(encoding="utf-8")) == {"value": 1},
                "atomic evidence write changed its payload")
        expect_failure(
            lambda: RUNNER.write_atomic_exclusive(atomic, {"value": 2}),
            "atomic evidence writer replaced an existing artifact")
        require(not list(root.glob(".atomic.json.tmp-*")),
                "atomic evidence writer leaked a temporary file")

        record, first_artifact = build_record(
            root, cell, manifest_sha256)
        validated = RUNNER.validate_ragged_cell_record(
            record, cell, manifest_sha256, root, 3)
        require(validated["accepted_process_count"] == 18 and
                validated["rejected_process_count"] == 0 and
                validated["rejected_attempt_count"] == 0 and
                validated["all_rounds_zero_sibling_nonidle"] is True,
                "valid streamed cell accounting changed")

        runner_wakeup = copy.deepcopy(record)
        runner_wakeup["presample"]["delta"]["benchmark_cpu"][
            "nonidle_jiffies"] = 1
        RUNNER.validate_ragged_cell_record(
            runner_wakeup, cell, manifest_sha256, root, 3)

        contaminated_presample = copy.deepcopy(record)
        contaminated_presample["presample"]["delta"]["benchmark_cpu"][
            "nonidle_jiffies"] = 2
        expect_failure(
            lambda: RUNNER.validate_ragged_cell_record(
                contaminated_presample, cell, manifest_sha256, root, 3),
            "cell with excess benchmark-CPU presample work was accepted")

        wrong_manifest = copy.deepcopy(record)
        wrong_manifest["manifest_sha256"] = "b" * 64
        expect_failure(
            lambda: RUNNER.validate_ragged_cell_record(
                wrong_manifest, cell, manifest_sha256, root, 3),
            "cell from a different manifest was accepted")

        wrong_lease = copy.deepcopy(record)
        wrong_lease["rounds"][0]["isolation"]["pair_lease"] = {
            "sha256": "2" * 64
        }
        expect_failure(
            lambda: RUNNER.validate_ragged_cell_record(
                wrong_lease, cell, manifest_sha256, root, 3),
            "cell with a cross-wired CPU lease was accepted")

        first_artifact.write_bytes(first_artifact.read_bytes()[:-1])
        expect_failure(
            lambda: RUNNER.validate_ragged_cell_record(
                record, cell, manifest_sha256, root, 3),
            "cell with a truncated invocation artifact was accepted")
    print("ragged T8 streamed-runner self-test passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
