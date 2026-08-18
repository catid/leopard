#!/usr/bin/env python3
"""Isolated contract checks for the 97-cell current-atlas ABBA screen."""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
import tempfile
from copy import deepcopy
from pathlib import Path
from typing import Any, Mapping


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_runner() -> Any:
    path = Path(__file__).resolve().with_name(
        "run_current_atlas_regression_screen_abba.py")
    specification = importlib.util.spec_from_file_location(
        "current_atlas_regression_screen_test_target", path)
    require(specification is not None and specification.loader is not None,
            "cannot load current-atlas regression runner")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


RUNNER = load_runner()
BASE = RUNNER.BASE
PROVENANCE_RUNNER = RUNNER.PARENT.SUPPORT
BUILD_PROVENANCE = PROVENANCE_RUNNER.BUILD_PROVENANCE


def require_rejected(callback: Any, message: str) -> None:
    try:
        callback()
    except BASE.EvidenceError:
        return
    raise RuntimeError(message)


def expected_old_ids() -> list[str]:
    result: list[str] = []

    def add(operation: str, k: int, byte_count: int, loss: int) -> None:
        result.append(f"{operation}-k{k}-r32-b{byte_count}-l{loss}")

    for k in (33, 35, 41, 43, 45, 47, 51, 57, 59, 61, 63, 64, 79, 91):
        add("encode", k, 64, 1)
    add("encode", 35, 1024, 1)
    for k in (65, 67, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95):
        add("decode", k, 64, 32)
    for k, loss in (
        (67, 7), (69, 7), (71, 8), (73, 8), (75, 8),
        (77, 8), (79, 8), (81, 9), (83, 9), (85, 9),
        (87, 9), (89, 9), (91, 10), (93, 10), (95, 10),
    ):
        add("decode", k, 1024, loss)
    for k in (81, 83, 85, 95):
        add("decode", k, 1024, 2)
    for k in (32, 34, 36, 62, 65, 80, 92, 96):
        add("encode", k, 64, 1)
    for k in (34, 36):
        add("encode", k, 1024, 1)
    for k in (64, 66, 96, 97):
        add("decode", k, 64, 32)
    for k, loss in ((65, 7), (66, 7), (96, 10), (97, 10)):
        add("decode", k, 1024, loss)
    for k in (79, 80, 86, 96, 97):
        add("decode", k, 1024, 2)
    for k in (65, 81, 95, 96, 97):
        add("decode", k, 1024, 1)
    return result


def metric(value: float, iterations: int = 31) -> dict[str, Any]:
    samples = [value] * iterations
    return {
        "median_us_per_batch_call": value,
        "mad_us_per_batch_call": 0.0,
        "minimum_us_per_batch_call": value,
        "maximum_us_per_batch_call": value,
        "samples_us_per_batch_call": samples,
    }


def payload(
    implementation: str, cell: Mapping[str, Any],
    *, encode: float = 2.0, decode: float = 3.0,
    iterations: int = 31, warmup: int = 64,
) -> dict[str, Any]:
    parent = 1
    while parent < cell["K"] + 32:
        parent *= 2
    common_parameters = {
        "K": cell["K"], "R": cell["R"],
        "shard_bytes": cell["bytes"], "loss_count": cell["loss"],
        "missing_original_indices": list(cell["missing_original_indices"]),
        "batch": 1, "reuse": cell["reuse"], "iterations": iterations,
        "warmup": warmup, "thread_count": 1, "seed": cell["seed"],
    }
    digests = {
        "algorithm": "fnv1a64", "original_data": "1" * 16,
        "transmitted_parity": "2" * 16,
        "recovered_originals": "3" * 16,
    }
    if implementation == "main":
        return {
            "schema": "leopard-main-benchmark-v1",
            "build": {"main_source_commit": BASE.MAIN_COMMIT,
                      "cplusplus": 201103},
            "parameters": {
                **common_parameters,
                "logical_shard_bytes": cell["bytes"],
            },
            "resolved": {
                "profile": "legacy_high_v1", "field": "gf8",
                "thread_count": 1, "parent_count": parent,
                "padded_side": 32, "padded_application_bytes": False,
                "padding_policy": "zero suffix per shard",
            },
            "correctness": {"round_trip": True,
                            "logical_prefix_fingerprinted": True},
            "workload_digests": digests,
            "memory": {
                "alignment": 64,
                "encode_work_shards_per_stripe": 32,
                "decode_work_shards_per_stripe": parent,
                "encode_work_bytes_batch": 32 * cell["bytes"],
                "decode_work_bytes_batch": parent * cell["bytes"],
            },
            "metrics": {
                "decode_timing_includes_setup": True,
                "encode_execution": metric(encode, iterations),
                "decode_including_setup": metric(decode, iterations),
            },
        }
    return {
        "schema": "leopard2-benchmark-v9",
        "build": {"compiler": "gcc", "compiler_version": "13.3.0",
                  "cplusplus": 201103},
        "parameters": {
            **common_parameters,
            "requested_profile": "legacy_high_v1",
            "requested_field": "gf8", "requested_backend": "avx2",
            "force_generic_decode": False,
            "force_specialized_decode": False,
            "force_tiled_decode": False,
            "force_materialized_decode": False,
            "skip_legacy": True, "retain_samples": True,
            "measure_one_shot_decode": True,
        },
        "resolved": {
            "profile": "legacy_high_v1", "field": "gf8",
            "backend": "avx2", "thread_count": 1,
            "parent_count": parent, "padded_side": 32,
        },
        "correctness": {"leopard2_round_trip": True,
                        "legacy_comparison": None},
        "workload_digests": digests,
        "memory": {
            "scratch_alignment": 64,
            "encode_scratch_bytes_per_stripe": 128,
            "encode_scratch_bytes_batch": 128,
            "decode_scratch_bytes_per_stripe": 256,
            "decode_scratch_bytes_batch": 256,
            "one_shot_decode_scratch_bytes_per_stripe": 384,
            "one_shot_decode_scratch_bytes_batch": 384,
        },
        "metrics": {
            "encode_execution": metric(encode, iterations),
            "one_shot_decode_including_setup": metric(decode, iterations),
        },
    }


def synthetic_round(
    order: tuple[str, ...], *, control_primary: float = 1.0,
    control_companion: float = 1.0, main_primary: float = 1.05,
    main_companion: float = 1.05, digest: str = "2" * 16,
    missing: list[int] | None = None,
) -> dict[str, Any]:
    values = {
        "candidate": (1.0, 1.0),
        "control": (control_primary, control_companion),
        "main": (main_primary, main_companion),
    }
    digests = {
        "original_data": "1" * 16, "transmitted_parity": digest,
        "recovered_originals": "3" * 16,
    }
    return {
        "invocations": [{
            "implementation": label,
            "normalized": {
                "encode_us": values[label][0],
                "one_shot_encode_us": values[label][1],
                "digests": dict(digests),
                "missing_original_indices": list(
                    [0] if missing is None else missing),
            },
        } for label in order],
        "isolation": {"accepted": True},
    }


def identity(path: Path) -> dict[str, Any]:
    return BASE.T8_SUPPORT.regular_file_identity(path)


def reset_shared_freeze() -> None:
    RUNNER.PARENT._SHARED_FROZEN_CANDIDATE = None
    RUNNER.PARENT._SHARED_FROZEN_PHYSICAL_IDENTITY = None


def exercise_shared_inode_contract() -> None:
    with tempfile.TemporaryDirectory(prefix="leo2-atlas-shared-") as directory:
        root = Path(directory)
        executable = root / "bench_leopard2"
        executable.write_bytes(b"one production executable\n")
        executable.chmod(0o755)
        digest = identity(executable)["sha256"]
        frozen = root / "frozen"
        frozen.mkdir()
        reset_shared_freeze()
        candidate = BASE.freeze_executable(
            executable, digest, frozen / "candidate")
        control = BASE.freeze_executable(
            executable, digest, frozen / "control")
        require(candidate["input"] == control["input"] and
                candidate["frozen"] == control["frozen"] and
                candidate["frozen_physical_identity"]["inode"] ==
                    control["frozen_physical_identity"]["inode"] and
                not (frozen / "control").exists(),
                "candidate/control did not share one frozen inode")

        other = root / "other"
        other.write_bytes(executable.read_bytes())
        other.chmod(0o755)
        second = root / "second"
        second.mkdir()
        reset_shared_freeze()
        BASE.freeze_executable(executable, digest, second / "candidate")
        require_rejected(
            lambda: BASE.freeze_executable(other, digest, second / "control"),
            "same bytes on distinct input inodes were accepted")


def exercise_round_journal() -> None:
    with tempfile.TemporaryDirectory(prefix="leo2-atlas-journal-") as directory:
        output = Path(directory)
        saved_pending = list(RUNNER._PENDING_INVOCATIONS)
        saved_output = RUNNER._PENDING_OUTPUT
        saved_sequence = RUNNER._ROUND_ARTIFACT_SEQUENCE
        try:
            RUNNER._PENDING_INVOCATIONS.clear()
            RUNNER._PENDING_OUTPUT = output
            RUNNER._ROUND_ARTIFACT_SEQUENCE = 0
            cell = RUNNER.cells()[0]
            source_commit, source_tree = "e" * 40, "f" * 40
            result = payload("candidate", cell)
            normalized = RUNNER.validate_result(
                "candidate", result, cell, source_commit, source_tree, 31, 64)
            executable = Path("/frozen/candidate")
            invocation = {
                "implementation": "candidate",
                "command": RUNNER.benchmark_command(
                    "candidate", executable, cell, 0, 31, 64),
                "elapsed_ns": 123,
                "stdout_sha256": "a" * 64,
                "stderr_sha256": "b" * 64,
                "normalized": normalized,
                "result": result,
            }
            RUNNER._PENDING_INVOCATIONS.append(invocation)
            RUNNER._persist_pending_round({"accepted": True})
            artifact_path = output / "round_payloads/round-000000.json"
            full = json.loads(artifact_path.read_text())
            require(full["invocations"][0]["result"]["schema"] ==
                        "leopard2-benchmark-v9" and
                    full["invocations"][0]["normalized"][
                        "encode_samples_us"] == [2.0] * 31 and
                    "result" not in invocation and
                    "command" not in invocation and
                    "encode_samples_us" not in invocation["normalized"] and
                    invocation["round_payload_artifact"]["sha256"] ==
                        BASE.sha256(artifact_path),
                    "round journal lost full evidence or failed to compact RAM")
            raw = {
                "schema": BASE.SCHEMA,
                "source_commit": source_commit, "source_tree": source_tree,
                "cpu": 0, "iterations": 31, "warmup": 64,
                "identities": {
                    "candidate": {"path": str(executable)},
                },
                "cells": [{"cell": cell, "rounds": [{
                    "order": ["candidate"],
                    "invocations": [invocation],
                    "isolation": {"accepted": True},
                    "discarded_attempts": [],
                }]}],
            }
            summary = {"process_count": 1, "discarded_process_count": 0}
            closure = RUNNER.validate_round_payload_artifacts(
                raw, summary, output, expected_cell_count=1)
            require(closure["artifact_count"] == 1 and
                    closure["accepted_process_count"] == 1,
                    "round payload closure rejected complete evidence")
            original_bytes = artifact_path.read_bytes()
            artifact_path.write_bytes(original_bytes + b" ")
            require_rejected(
                lambda: RUNNER.validate_round_payload_artifacts(
                    raw, summary, output, expected_cell_count=1),
                "mutated round payload survived closure validation")
            artifact_path.write_bytes(original_bytes)
            original_reference = deepcopy(
                invocation["round_payload_artifact"])

            invocation["round_payload_artifact"]["invocation_index"] = 1
            require_rejected(
                lambda: RUNNER.validate_round_payload_artifacts(
                    raw, summary, output, expected_cell_count=1),
                "wrong invocation index survived closure validation")
            invocation["round_payload_artifact"] = deepcopy(original_reference)

            orphan = output / "round_payloads/round-999999.json"
            orphan.write_text("{}\n", encoding="utf-8")
            require_rejected(
                lambda: RUNNER.validate_round_payload_artifacts(
                    raw, summary, output, expected_cell_count=1),
                "orphan round payload survived closure validation")
            orphan.unlink()

            def publish_mutation(value: dict[str, Any]) -> None:
                artifact_path.write_text(json.dumps(
                    value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
                observed = BASE.support_file_identity(artifact_path)
                invocation["round_payload_artifact"] = {
                    **observed, "invocation_index": 0,
                }

            changed = deepcopy(full)
            changed["invocations"][0]["result"]["parameters"]["reuse"] += 1
            publish_mutation(changed)
            require_rejected(
                lambda: RUNNER.validate_round_payload_artifacts(
                    raw, summary, output, expected_cell_count=1),
                "full child-result mutation survived semantic validation")
            changed = deepcopy(full)
            changed["invocations"][0]["command"].append("--unknown")
            publish_mutation(changed)
            require_rejected(
                lambda: RUNNER.validate_round_payload_artifacts(
                    raw, summary, output, expected_cell_count=1),
                "full child-command mutation survived semantic validation")
            changed = deepcopy(full)
            changed["isolation"] = {"accepted": False}
            publish_mutation(changed)
            require_rejected(
                lambda: RUNNER.validate_round_payload_artifacts(
                    raw, summary, output, expected_cell_count=1),
                "full isolation mutation survived semantic validation")
            artifact_path.write_bytes(original_bytes)
            invocation["round_payload_artifact"] = deepcopy(original_reference)

            pending = deepcopy(full["invocations"][0])
            RUNNER._PENDING_INVOCATIONS.append(pending)
            RUNNER._PENDING_OUTPUT = output
            saved_run_one = RUNNER._SHARED_RUN_ONE

            def fail_child(*_args: Any, **_kwargs: Any) -> Any:
                raise BASE.EvidenceError("synthetic mid-round child failure")

            RUNNER._SHARED_RUN_ONE = fail_child
            try:
                require_rejected(lambda: RUNNER.run_one_with_round_journal(
                    "control", {}, {}, 0, "", "", 3, 1, output),
                    "mid-round child failure was accepted")
            finally:
                RUNNER._SHARED_RUN_ONE = saved_run_one
            require(not RUNNER._PENDING_INVOCATIONS and
                    (output / "round_payloads/round-000001.json").is_file(),
                    "mid-round successful children were not journaled")
            BASE.T8_SUPPORT.write_exclusive(
                output / "failure.json", {"failure": {"type": "synthetic"}})
            RUNNER._bind_failure_journal(output)
            failure = json.loads((output / "failure.json").read_text())
            require(len(failure["round_payload_artifacts"]) == 2,
                    "failure evidence did not bind partial round artifacts")
        finally:
            RUNNER._PENDING_INVOCATIONS.clear()
            RUNNER._PENDING_INVOCATIONS.extend(saved_pending)
            RUNNER._PENDING_OUTPUT = saved_output
            RUNNER._ROUND_ARTIFACT_SEQUENCE = saved_sequence


def fake_production_provenance(
    build: Path, executable: Path, commit: str, tree: str,
) -> dict[str, Any]:
    files = {
        "cmake_cache": build / "CMakeCache.txt",
        "compile_commands": build / "compile_commands.json",
        "archive_link_recipe": build / "CMakeFiles/leopard.dir/link.txt",
        "executable_link_recipe":
            build / "CMakeFiles/bench_leopard2.dir/link.txt",
        "archive": build / "libleopard.a",
    }
    for index, path in enumerate(files.values()):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(f"fixture {index}\n".encode())
    archive_object = build / "CMakeFiles/leopard.dir/leopard2.cpp.o"
    benchmark_object = (
        build / "CMakeFiles/bench_leopard2.dir/bench/leopard2/benchmark.cpp.o")
    archive_object.parent.mkdir(parents=True, exist_ok=True)
    benchmark_object.parent.mkdir(parents=True, exist_ok=True)
    archive_object.write_bytes(b"archive object\n")
    benchmark_object.write_bytes(b"benchmark object\n")
    archive_identity = identity(archive_object)
    benchmark_identity = identity(benchmark_object)
    dependencies = []
    for dependency in BASE.RUNNER_DEPENDENCIES:
        source_identity = BASE.support_file_identity(dependency)
        dependencies.append({
            "path": dependency.resolve().relative_to(
                PROVENANCE_RUNNER.SOURCE_ROOT).as_posix(),
            "sha256": source_identity["sha256"],
            "size": source_identity["size"], "mode": 0o644,
        })
    return {
        "schema": BUILD_PROVENANCE.PRODUCTION_BUILD_CLOSURE_SCHEMA,
        "build_root": str(build),
        "physical_source_root": str(PROVENANCE_RUNNER.SOURCE_ROOT),
        "source_root": str(PROVENANCE_RUNNER.SOURCE_ROOT),
        "executable_target": "bench_leopard2",
        "tracked_source_manifest": {
            "files": dependencies,
            "git": {"commit": commit, "tree": tree, "dirty": False},
        },
        "validated_cache": {
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                BUILD_PROVENANCE.BENCHMARK_BUILD_CONFIGURATION_SCHEMA,
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256": "9" * 64,
            "LEO2_BUILD_BENCHMARKS": "ON", "LEOPARD_ENABLE_GF8": "ON",
            "LEO2_BACKEND_VARIANT": "avx2",
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR": "OFF",
            "LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING": "ON",
        },
        **{name: identity(path) for name, path in files.items()},
        "executable": identity(executable),
        "source_object_compile_closure": [
            {
                "role": "archive",
                "source": identity(PROVENANCE_RUNNER.SOURCE_ROOT /
                                   "leopard2.cpp"),
                "object": archive_identity,
                "compile_entry": {"arguments": ["c++"]},
                "flag_profile": "portable-core",
            },
            {
                "role": "benchmark",
                "source": identity(PROVENANCE_RUNNER.SOURCE_ROOT /
                                   "bench/leopard2/benchmark.cpp"),
                "object": benchmark_identity,
                "compile_entry": {"arguments": ["c++"]},
                "flag_profile": "portable-core",
            },
        ],
        "archive_members": [archive_object.name],
        "archive_member_identities": [{
            "member": archive_object.name,
            "sha256": archive_identity["sha256"],
            "size": archive_identity["size"],
        }],
        "archive_link_commands": [["ar"], ["ranlib"]],
        "executable_link_command": ["c++"],
    }


def exercise_provenance_contract() -> None:
    with tempfile.TemporaryDirectory(prefix="leo2-atlas-provenance-") as directory:
        build = Path(directory) / "build"
        build.mkdir()
        executable = build / "bench_leopard2"
        executable.write_bytes(b"fixture executable\n")
        executable.chmod(0o755)
        digest = identity(executable)["sha256"]
        commit, tree = "e" * 40, "f" * 40
        pristine = fake_production_provenance(build, executable, commit, tree)
        options = argparse.Namespace(
            candidate=executable, control=executable,
            candidate_sha256=digest, control_sha256=digest,
            source_commit=commit, source_tree=tree,
        )
        frozen = Path(directory) / "frozen"
        frozen.mkdir()
        reset_shared_freeze()
        BASE.freeze_executable(executable, digest, frozen / "candidate")
        BASE.freeze_executable(executable, digest, frozen / "control")

        saved_capture = BUILD_PROVENANCE.candidate_build_provenance
        saved_replay = BUILD_PROVENANCE.verify_reproducible_candidate_build
        saved_validate = BUILD_PROVENANCE.validate_reproducible_build_proof
        current = {"value": pristine}
        calls: list[str] = []

        def capture(*_args: Any, **_kwargs: Any) -> dict[str, Any]:
            calls.append("capture")
            return deepcopy(current["value"])

        def replay(*_args: Any, **kwargs: Any) -> dict[str, Any]:
            require(kwargs == {"jobs": 1}, "production replay lost jobs=1")
            calls.append("replay")
            return {"schema": "mock-reproducible-build", "valid": True}

        def validate(proof: Any, _captured: Any, **_kwargs: Any) -> None:
            require(proof == {"schema": "mock-reproducible-build",
                              "valid": True}, "wrong replay proof")
            calls.append("validate")

        def reset(value: dict[str, Any]) -> None:
            PROVENANCE_RUNNER._BUILD_CLOSURE_BASELINE = None
            current["value"] = value
            calls.clear()

        try:
            BUILD_PROVENANCE.candidate_build_provenance = capture
            BUILD_PROVENANCE.verify_reproducible_candidate_build = replay
            BUILD_PROVENANCE.validate_reproducible_build_proof = validate
            reset(pristine)
            closure = BASE.build_closure_identity(options)
            require(len(closure["runner_source_dependencies"]) ==
                    len(BASE.RUNNER_DEPENDENCIES) and
                    calls == ["capture", "replay", "validate"],
                    "production replay/dependency closure was bypassed")
            BASE.build_closure_identity(options)
            require(calls == ["capture", "replay", "validate",
                              "capture", "validate"],
                    "post-run build closure did not recapture/revalidate")

            malformed = deepcopy(pristine)
            malformed["tracked_source_manifest"]["files"].pop()
            reset(malformed)
            require_rejected(
                lambda: BASE.build_closure_identity(options),
                "missing runner dependency was accepted")
        finally:
            BUILD_PROVENANCE.candidate_build_provenance = saved_capture
            BUILD_PROVENANCE.verify_reproducible_candidate_build = saved_replay
            BUILD_PROVENANCE.validate_reproducible_build_proof = saved_validate
            PROVENANCE_RUNNER._BUILD_CLOSURE_BASELINE = None


def exercise_atlas_contract(cells: list[dict[str, Any]]) -> None:
    manifest_by_id = {
        item["id"]: item for item in RUNNER.ATLAS_MANIFEST["cells"]}
    rows_by_cell: dict[str, dict[str, Mapping[str, Any]]] = {}
    for row in RUNNER.ATLAS_SUMMARY["rows"]:
        rows_by_cell.setdefault(row["cell_id"], {})[row["codec"]] = row
    losing = set()
    for atlas_id, rows in rows_by_cell.items():
        if "leopard1" not in rows:
            continue
        leopard2, main = rows["leopard2"], rows["leopard1"]
        if (main["encode_execution_us"] / leopard2["encode_execution_us"] <
                0.98 or
                main["decode_first_us"] / leopard2["decode_first_us"] <
                0.98):
            losing.add(atlas_id)
    curated_target_atlas_ids = {
        cell["atlas_cell_id"] for cell in cells[:49]}
    require(losing - curated_target_atlas_ids ==
            set(RUNNER.SUPPLEMENTAL_TARGET_IDS) and len(losing) == 56,
            "supplemental IDs differ from exact atlas <0.98 union-minus-curated")
    for cell in cells:
        if cell["workload_origin"] == "committed_atlas_manifest":
            source = manifest_by_id[cell["atlas_cell_id"]]
            require(
                (cell["K"], cell["R"], cell["bytes"], cell["loss"],
                 cell["seed"], cell["reuse"]) ==
                (source["K"], source["R"], source["shard_bytes"],
                 source["loss_count"], source["seed"], source["reuse"]),
                f"manifest cell drifted: {cell['id']}")
        else:
            require(cell["seed"] == RUNNER.ATLAS_GENERATOR.seed_for(
                        cell["K"], cell["bytes"], cell["loss"]) and
                    cell["reuse"] ==
                        RUNNER.ATLAS_GENERATOR.repetition_policy(
                            cell["K"], cell["bytes"])[0],
                    f"canonical derived control drifted: {cell['id']}")
        require(cell["missing_original_indices"] ==
                RUNNER._select_losses(
                    cell["K"], cell["loss"], cell["seed"]),
                f"erasure-pattern derivation drifted: {cell['id']}")

    with tempfile.TemporaryDirectory(prefix="leo2-atlas-json-") as directory:
        root = Path(directory)
        manifest = root / "manifest.json"
        summary = root / "summary.json"
        manifest.write_bytes(RUNNER.ATLAS_MANIFEST_PATH.read_bytes())
        summary.write_bytes(RUNNER.ATLAS_SUMMARY_PATH.read_bytes())
        loaded, loaded_summary = RUNNER.load_atlas(manifest, summary)
        require(loaded == RUNNER.ATLAS_MANIFEST and
                loaded_summary == RUNNER.ATLAS_SUMMARY,
                "exact atlas copies did not load identically")

        mutated = bytearray(manifest.read_bytes())
        mutated[-2] = ord(" ") if mutated[-2] != ord(" ") else ord("\t")
        manifest.write_bytes(mutated)
        require_rejected(lambda: RUNNER.load_atlas(manifest, summary),
                         "manifest byte mutation was accepted")
        manifest.write_bytes(RUNNER.ATLAS_MANIFEST_PATH.read_bytes())

        value = json.loads(summary.read_text())
        value["rows"][0]["seed"] += 1
        summary.write_text(json.dumps(value), encoding="utf-8")
        saved_hash = RUNNER.ATLAS_SUMMARY_SHA256
        RUNNER.ATLAS_SUMMARY_SHA256 = BASE.sha256(summary)
        try:
            require_rejected(lambda: RUNNER.load_atlas(manifest, summary),
                             "summary row mutation survived semantic checks")
        finally:
            RUNNER.ATLAS_SUMMARY_SHA256 = saved_hash


def exercise_payload_contract(cells: list[dict[str, Any]]) -> None:
    encode_cell = next(cell for cell in cells
                       if cell["operation"] == "encode")
    decode_cell = next(cell for cell in cells
                       if cell["operation"] == "decode")
    for implementation in ("candidate", "control", "main"):
        normalized = RUNNER.validate_result(
            implementation, payload(implementation, encode_cell,
                                    encode=2.0, decode=3.0),
            encode_cell, "e" * 40, "f" * 40, 31, 64)
        require(normalized["encode_us"] == 2.0 and
                normalized["one_shot_encode_us"] == 3.0,
                "encode-primary metric mapping changed")
        normalized = RUNNER.validate_result(
            implementation, payload(implementation, decode_cell,
                                    encode=2.0, decode=3.0),
            decode_cell, "e" * 40, "f" * 40, 31, 64)
        require(normalized["encode_us"] == 3.0 and
                normalized["one_shot_encode_us"] == 2.0,
                "decode-primary metric mapping changed")

    malformed = payload("candidate", encode_cell)
    malformed["metrics"]["encode_execution"][
        "median_us_per_batch_call"] += 0.1
    require_rejected(
        lambda: RUNNER.validate_result(
            "candidate", malformed, encode_cell, "e" * 40, "f" * 40,
            31, 64),
        "forged metric median was accepted")
    malformed = payload("candidate", encode_cell)
    malformed["metrics"]["encode_execution"][
        "samples_us_per_batch_call"].pop()
    require_rejected(
        lambda: RUNNER.validate_result(
            "candidate", malformed, encode_cell, "e" * 40, "f" * 40,
            31, 64),
        "incomplete metric samples were accepted")
    malformed = payload("candidate", encode_cell)
    malformed["schema"] = "leopard2-benchmark-v8"
    require_rejected(
        lambda: RUNNER.validate_result(
            "candidate", malformed, encode_cell, "e" * 40, "f" * 40,
            31, 64), "stale current schema was accepted")
    malformed = payload("candidate", encode_cell)
    malformed["parameters"]["reuse"] += 1
    require_rejected(
        lambda: RUNNER.validate_result(
            "candidate", malformed, encode_cell, "e" * 40, "f" * 40,
            31, 64), "wrong workload geometry was accepted")
    malformed = payload("candidate", encode_cell)
    malformed["parameters"]["missing_original_indices"] = [
        (encode_cell["missing_original_indices"][0] + 1) % encode_cell["K"]]
    require_rejected(
        lambda: RUNNER.validate_result(
            "candidate", malformed, encode_cell, "e" * 40, "f" * 40,
            31, 64), "wrong erasure pattern was accepted")
    malformed = payload("main", encode_cell)
    malformed["build"]["main_source_commit"] = "0" * 40
    require_rejected(
        lambda: RUNNER.validate_result(
            "main", malformed, encode_cell, "e" * 40, "f" * 40, 31, 64),
        "noncanonical main source was accepted")


def exercise_argument_contract() -> None:
    digest = "a" * 64
    arguments = [
        str(BASE.RUNNER_PATH), "--candidate", "/build/bench_leopard2",
        "--candidate-sha256", digest,
        "--control", "/build/bench_leopard2",
        "--control-sha256", digest, "--main", "/frozen/main",
        "--main-sha256", BASE.CANONICAL_MAIN_SHA256,
        "--source-commit", "e" * 40, "--source-tree", "f" * 40,
        "--output", "/results/current-atlas", "--cpu", "30",
        "--sibling", "31",
    ]
    saved_argv = sys.argv
    try:
        sys.argv = arguments
        parsed = BASE.parse_arguments()
        require(parsed.candidate == parsed.control and parsed.rounds == 25 and
                parsed.iterations == 31 and parsed.warmup == 64,
                "authoritative CLI defaults changed")
        sys.argv = arguments[:]
        index = sys.argv.index("--main-sha256") + 1
        sys.argv[index] = "b" * 64
        require_rejected(BASE.parse_arguments,
                         "noncanonical main SHA-256 was accepted")
    finally:
        sys.argv = saved_argv


def main() -> int:
    cells = RUNNER.cells()
    expected_ids = expected_old_ids() + list(RUNNER.SUPPLEMENTAL_TARGET_IDS)
    require([cell["id"] for cell in cells] == expected_ids and
            len(cells) == 97 and
            sum(cell["role"] == "target" for cell in cells) == 69 and
            sum(cell["role"] == "neighbor" for cell in cells) == 28 and
            sum(cell["workload_origin"] == "committed_atlas_manifest"
                for cell in cells) == 81 and
            sum(cell["workload_origin"] ==
                "canonical_policy_derived_control" for cell in cells) == 16,
            "exact 97-cell classification changed")
    require(
        BASE.SCHEMA == "leopard2-current-atlas-regression-screen-abba/v2" and
        BASE.SUMMARY_SCHEMA ==
            "leopard2-current-atlas-regression-screen-preliminary-summary/v2" and
        RUNNER.FINAL_SUMMARY_SCHEMA ==
            "leopard2-current-atlas-regression-screen-summary/v2" and
        BASE.CANDIDATE_SCHEMA == BASE.CONTROL_SCHEMA ==
            "leopard2-benchmark-v9" and
        BASE.CONTROL_EXTRA_ARGUMENTS == () and
        BASE.CONTROL_BUILD_MARKER is None and
        BASE.MODE_SYMBOL is None and
        BASE.ALLOW_IDENTICAL_CANDIDATE_CONTROL is True and
        BASE.ALLOW_MULTIPLE_TARGETS is True and
        BASE.TARGET_CONTROL_FLOOR == 1.0 / 1.02 and
        BASE.TARGET_MAIN_FLOOR == 1.0 and
        BASE.NEIGHBOR_FLOOR == 1.0 / 1.02 and
        BASE.MAX_ISOLATION_ATTEMPTS == 8 and
        BASE.CANONICAL_MAIN_SHA256 ==
            "78575e2fe1b9796d10dcb4c14d2bcc1b627a79e787be66d34dc7aa666899aa93" and
        RUNNER.PARENT.SUPPORT.__doc__ == RUNNER.__doc__,
        "schema, immutable labels, or acceptance gates changed")
    expected_dependencies = (
        Path(RUNNER.__file__).resolve(), *RUNNER.SUPPORT_DEPENDENCIES,
        RUNNER.ATLAS_MANIFEST_PATH.resolve(),
        RUNNER.ATLAS_SUMMARY_PATH.resolve(),
        RUNNER.ATLAS_GENERATOR_PATH.resolve(),
    )
    require(BASE.RUNNER_DEPENDENCIES == expected_dependencies and
            len(set(expected_dependencies)) == len(expected_dependencies),
            "source dependency attestation did not inherit the full closure")

    representative = cells[0]
    candidate = BASE.benchmark_command(
        "candidate", Path("/frozen/shared"), representative, 30, 31, 64)
    control = BASE.benchmark_command(
        "control", Path("/frozen/shared"), representative, 30, 31, 64)
    main_command = BASE.benchmark_command(
        "main", Path("/frozen/main"), representative, 30, 31, 64)
    forbidden = {
        "--attest-source", "--measure-one-shot-encode",
        "--balanced-b64-terminal-mode", "--k66r16-b64-tail-mode",
    }
    require(candidate == control and
            "--measure-one-shot-decode" in candidate and
            forbidden.isdisjoint(candidate) and
            "--measure-one-shot-decode" not in main_command and
            main_command[main_command.index("--reuse") + 1] ==
                str(representative["reuse"]),
            "same-route command or standalone v9 restriction changed")

    orders = BASE.select_round_orders(BASE.TARGET_ORDER, 25)
    require(len(orders) == 25 and orders[:3] == BASE.TARGET_ORDER,
            "25-round ABBA schedule changed")
    analysis = BASE.analyze(
        representative, [synthetic_round(
            order, missing=representative["missing_original_indices"])
            for order in orders])
    require(analysis["historical_role"] == "former_major" and
            analysis["metric_mapping"]["primary_operation"] == "encode" and
            analysis["acceptance_observations"][
                "same_binary_equivalence_pass"] ==
                    {"primary": True, "companion": True} and
            analysis["acceptance_observations"][
                "exact_main_lower_ci_pass"] ==
                    {"primary": True, "companion": True},
            "25-round ratio annotation changed")
    below = BASE.analyze(
        representative,
        [synthetic_round(
            order, main_primary=0.99,
            missing=representative["missing_original_indices"])
            for order in orders])
    require(below["acceptance_observations"][
                "exact_main_lower_ci_pass"]["primary"] is False,
            "pre-raw analysis hid an exact-main regression")
    corrupted = [synthetic_round(
        order, missing=representative["missing_original_indices"])
        for order in orders]
    corrupted[0]["invocations"][0]["normalized"]["digests"][
        "transmitted_parity"] = "f" * 16
    require_rejected(lambda: BASE.analyze(representative, corrupted),
                     "cross-label digest mismatch was accepted")

    wrong_pattern = [synthetic_round(
        order, missing=representative["missing_original_indices"])
        for order in orders]
    wrong_pattern[0]["invocations"][0]["normalized"][
        "missing_original_indices"] = [
            (representative["missing_original_indices"][0] + 1) %
            representative["K"]]
    require_rejected(lambda: BASE.analyze(representative, wrong_pattern),
                     "cross-label erasure-pattern mismatch was accepted")

    analyses = [BASE.analyze(cell, [synthetic_round(
        order, missing=cell["missing_original_indices"])
        for order in orders]) for cell in cells]
    preliminary = {
        "schema": BASE.SUMMARY_SCHEMA, "cells": analyses,
    }
    accepted = RUNNER.apply_final_acceptance(deepcopy(preliminary))
    require(accepted["status"] == "accepted" and
            accepted["schema"] == RUNNER.FINAL_SUMMARY_SCHEMA and
            not accepted["exact_main_regressions"] and
            not accepted["same_binary_equivalence_failures"],
            "valid all-cell acceptance summary was rejected")
    require(preliminary["schema"] != RUNNER.FINAL_SUMMARY_SCHEMA,
            "preliminary generic output uses the final evidence schema")
    substituted = deepcopy(preliminary)
    substituted["cells"][-1] = deepcopy(substituted["cells"][0])
    require_rejected(
        lambda: RUNNER.apply_final_acceptance(substituted),
        "duplicate/substituted cell matrix was accepted")
    biased = deepcopy(preliminary)
    biased["cells"][0]["control_over_candidate"]["ci95"] = [1.04, 1.06]
    RUNNER.apply_final_acceptance(biased)
    require(biased["status"] == "rejected" and
            biased["same_binary_equivalence_failures"][0]["operation"] ==
                "encode",
            "candidate-favoring same-binary label bias was accepted")
    neighbor_regression = deepcopy(preliminary)
    neighbor = next(item for item in neighbor_regression["cells"]
                    if item["cell"]["role"] == "neighbor")
    neighbor["one_shot_main_over_candidate"]["ci95"] = [0.97, 0.99]
    RUNNER.apply_final_acceptance(neighbor_regression)
    require(neighbor_regression["status"] == "rejected" and
            any(item["cell_id"] == neighbor["cell"]["id"]
                for item in neighbor_regression["exact_main_regressions"]),
            "neighbor exact-main regression was accepted")

    exercise_atlas_contract(cells)
    exercise_payload_contract(cells)
    exercise_argument_contract()
    exercise_shared_inode_contract()
    exercise_round_journal()
    exercise_provenance_contract()
    print("Current-atlas 97-cell hardened runner checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
