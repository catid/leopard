#!/usr/bin/env python3
"""Deterministic fail-closed checks for the K65/R65/B64 ABBA runner."""

from __future__ import annotations

import argparse
import contextlib
import importlib.util
import io
import math
import statistics
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
        "run_k65r65_b64_packed_terminal_abba.py")
    specification = importlib.util.spec_from_file_location(
        "leopard2_k65r65_b64_runner_test_target", path)
    require(specification is not None and specification.loader is not None,
            "cannot load K65/R65/B64 runner")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


RUNNER = load_runner()
BASE = RUNNER.BASE
K66 = RUNNER.SUPPORT
PRODUCTION = RUNNER.PRODUCTION_SUPPORT
BUILD_PROVENANCE = PRODUCTION.BUILD_PROVENANCE


def require_rejected(callback: Any, message: str) -> None:
    try:
        callback()
    except BASE.EvidenceError:
        return
    raise RuntimeError(message)


def require_argument_rejected(callback: Any, message: str) -> None:
    try:
        with contextlib.redirect_stderr(io.StringIO()):
            callback()
    except SystemExit as error:
        require(error.code != 0, "malformed arguments exited successfully")
        return
    raise RuntimeError(message)


def identity(path: Path) -> dict[str, Any]:
    return BASE.T8_SUPPORT.regular_file_identity(path)


def reset_shared_freezer() -> None:
    K66._SHARED_FROZEN_CANDIDATE = None
    K66._SHARED_FROZEN_PHYSICAL_IDENTITY = None


def exercise_shared_inode_freezer() -> None:
    with tempfile.TemporaryDirectory(
            prefix="leo2-k65r65-shared-freeze-") as directory:
        root = Path(directory)
        executable = root / "bench_leopard2"
        executable.write_bytes(b"shared K65/R65 executable fixture\n")
        executable.chmod(0o755)
        main = root / "bench_main"
        main.write_bytes(b"canonical-main fixture\n")
        main.chmod(0o755)
        executable_sha256 = identity(executable)["sha256"]
        main_sha256 = identity(main)["sha256"]
        frozen = root / "frozen"
        frozen.mkdir()

        reset_shared_freezer()
        candidate = BASE.freeze_executable(
            executable, executable_sha256, frozen / "candidate")
        control = BASE.freeze_executable(
            executable, executable_sha256, frozen / "control")
        frozen_main = BASE.freeze_executable(
            main, main_sha256, frozen / "main")
        require(
            candidate["input"] == control["input"] and
            candidate["frozen"] == control["frozen"] and
            candidate["frozen"]["path"] ==
                str((frozen / "candidate").resolve()) and
            candidate["frozen_physical_identity"] ==
                control["frozen_physical_identity"] and
            candidate["frozen_physical_identity"]["device"] ==
                (frozen / "candidate").stat().st_dev and
            candidate["frozen_physical_identity"]["inode"] ==
                (frozen / "candidate").stat().st_ino and
            control["shared_physical_executable"] == "candidate" and
            not (frozen / "control").exists() and
            Path(frozen_main["frozen"]["path"]) ==
                (frozen / "main").resolve(),
            "candidate/control were not frozen to one physical inode")
        K66.require_shared_frozen_physical_identity()

        reset_shared_freezer()
        control_first = root / "control-first"
        control_first.mkdir()
        require_rejected(
            lambda: BASE.freeze_executable(
                executable, executable_sha256, control_first / "control"),
            "control was accepted before the shared candidate")

        different = root / "same-bytes-different-inode"
        different.write_bytes(executable.read_bytes())
        different.chmod(0o755)
        mismatch = root / "mismatch"
        mismatch.mkdir()
        reset_shared_freezer()
        BASE.freeze_executable(
            executable, executable_sha256, mismatch / "candidate")
        require_rejected(
            lambda: BASE.freeze_executable(
                different, executable_sha256, mismatch / "control"),
            "distinct same-byte input inodes were accepted")

        replaced = root / "replaced"
        replaced.mkdir()
        reset_shared_freezer()
        record = BASE.freeze_executable(
            executable, executable_sha256, replaced / "candidate")
        frozen_path = Path(record["frozen"]["path"])
        replacement = replaced / "replacement"
        replacement.write_bytes(frozen_path.read_bytes())
        replacement.chmod(0o555)
        replacement.replace(frozen_path)
        require_rejected(
            K66.require_shared_frozen_physical_identity,
            "same-byte replacement of the shared inode was accepted")
    reset_shared_freezer()


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
        path.write_bytes(f"fixture {index}\n".encode("ascii"))

    archive_object = build / "CMakeFiles/leopard.dir/leopard2.cpp.o"
    benchmark_object = (
        build / "CMakeFiles/bench_leopard2.dir/bench/leopard2/benchmark.cpp.o")
    archive_object.parent.mkdir(parents=True, exist_ok=True)
    benchmark_object.parent.mkdir(parents=True, exist_ok=True)
    archive_object.write_bytes(b"archive object\n")
    benchmark_object.write_bytes(b"benchmark object\n")
    archive_object_identity = identity(archive_object)
    benchmark_object_identity = identity(benchmark_object)

    dependencies = []
    for dependency in BASE.RUNNER_DEPENDENCIES:
        source_identity = BASE.support_file_identity(dependency)
        dependencies.append({
            "path": dependency.resolve().relative_to(
                PRODUCTION.SOURCE_ROOT).as_posix(),
            "sha256": source_identity["sha256"],
            "size": source_identity["size"],
            "mode": 0o644,
        })
    return {
        "schema": BUILD_PROVENANCE.PRODUCTION_BUILD_CLOSURE_SCHEMA,
        "build_root": str(build),
        "physical_source_root": str(PRODUCTION.SOURCE_ROOT),
        "source_root": str(PRODUCTION.SOURCE_ROOT),
        "executable_target": "bench_leopard2",
        "tracked_source_manifest": {
            "files": dependencies,
            "git": {"commit": commit, "tree": tree, "dirty": False},
        },
        "validated_cache": {
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                BUILD_PROVENANCE.BENCHMARK_BUILD_CONFIGURATION_SCHEMA,
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256": "9" * 64,
            "LEO2_BUILD_BENCHMARKS": "ON",
            "LEOPARD_ENABLE_GF8": "ON",
            "LEO2_BACKEND_VARIANT": "avx2",
            "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR": "OFF",
            "LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING": "ON",
        },
        **{name: identity(path) for name, path in files.items()},
        "executable": identity(executable),
        "source_object_compile_closure": [
            {
                "role": "archive",
                "source": identity(PRODUCTION.SOURCE_ROOT / "leopard2.cpp"),
                "object": archive_object_identity,
                "compile_entry": {"arguments": ["c++"]},
                "flag_profile": "portable-core",
            },
            {
                "role": "benchmark",
                "source": identity(
                    PRODUCTION.SOURCE_ROOT / "bench/leopard2/benchmark.cpp"),
                "object": benchmark_object_identity,
                "compile_entry": {"arguments": ["c++"]},
                "flag_profile": "portable-core",
            },
        ],
        "archive_members": [archive_object.name],
        "archive_member_identities": [{
            "member": archive_object.name,
            "sha256": archive_object_identity["sha256"],
            "size": archive_object_identity["size"],
        }],
        "archive_link_commands": [["ar"], ["ranlib"]],
        "executable_link_command": ["c++"],
    }


def exercise_production_build_closure() -> None:
    with tempfile.TemporaryDirectory(
            prefix="leo2-k65r65-provenance-") as directory:
        build = Path(directory)
        executable = build / "bench_leopard2"
        executable.write_bytes(b"production executable fixture\n")
        executable.chmod(0o755)
        commit = "e" * 40
        tree = "f" * 40
        pristine = fake_production_provenance(
            build, executable, commit, tree)
        executable_sha256 = identity(executable)["sha256"]
        options = argparse.Namespace(
            candidate=executable, control=executable,
            candidate_sha256=executable_sha256,
            control_sha256=executable_sha256,
            source_commit=commit, source_tree=tree)

        frozen = build / "frozen-fixture"
        frozen.mkdir()
        reset_shared_freezer()
        BASE.freeze_executable(
            executable, executable_sha256, frozen / "candidate")
        BASE.freeze_executable(
            executable, executable_sha256, frozen / "control")

        saved_capture = BUILD_PROVENANCE.candidate_build_provenance
        saved_replay = BUILD_PROVENANCE.verify_reproducible_candidate_build
        saved_validate = BUILD_PROVENANCE.validate_reproducible_build_proof
        calls: list[tuple[Any, ...]] = []
        current = {"value": pristine}

        def capture(*arguments: Any, **keywords: Any) -> dict[str, Any]:
            calls.append(("capture", arguments, keywords))
            return deepcopy(current["value"])

        def replay(*arguments: Any, **keywords: Any) -> dict[str, Any]:
            calls.append(("replay", arguments, keywords))
            return {"schema": "mock-replay", "valid": True}

        def validate(
            proof: Any, candidate: Any, **keywords: Any,
        ) -> None:
            calls.append(("validate", proof, candidate, keywords))
            if proof != {"schema": "mock-replay", "valid": True}:
                raise RuntimeError("mock replay proof rejected")

        def reset(value: dict[str, Any] | None = None) -> None:
            PRODUCTION._BUILD_CLOSURE_BASELINE = None
            current["value"] = pristine if value is None else value
            calls.clear()

        try:
            BUILD_PROVENANCE.candidate_build_provenance = capture
            BUILD_PROVENANCE.verify_reproducible_candidate_build = replay
            BUILD_PROVENANCE.validate_reproducible_build_proof = validate

            reset()
            closure = BASE.build_closure_identity(options)
            require(
                closure["candidate_build"] == pristine and
                closure["reproducible_build"] ==
                    {"schema": "mock-replay", "valid": True} and
                len(closure["runner_source_dependencies"]) ==
                    len(BASE.RUNNER_DEPENDENCIES),
                "production provenance/replay result changed")
            require(
                calls[0] == (
                    "capture",
                    (build.resolve(), PRODUCTION.SOURCE_ROOT,
                     executable.resolve(), "bench_leopard2"),
                    {}) and
                any(call[0] == "replay" and call[2] == {"jobs": 1}
                    for call in calls) and
                any(call[0] == "validate" for call in calls),
                "production target or single-job replay was not bound")
            BASE.build_closure_identity(options)
            require(sum(call[0] == "capture" for call in calls) == 2 and
                    sum(call[0] == "replay" for call in calls) == 1 and
                    sum(call[0] == "validate" for call in calls) == 2,
                    "retained replay proof was not recaptured/revalidated")

            dirty = deepcopy(pristine)
            dirty["tracked_source_manifest"]["git"]["dirty"] = True
            reset(dirty)
            require_rejected(
                lambda: BASE.build_closure_identity(options),
                "dirty production source was accepted")

            missing_dependency = deepcopy(pristine)
            missing_dependency["tracked_source_manifest"]["files"].pop()
            reset(missing_dependency)
            require_rejected(
                lambda: BASE.build_closure_identity(options),
                "incomplete runner source closure was accepted")

            wrong_control = argparse.Namespace(
                candidate=executable, control=build / "another",
                candidate_sha256=executable_sha256,
                control_sha256=executable_sha256,
                source_commit=commit, source_tree=tree)
            wrong_control.control.write_bytes(b"other executable\n")
            wrong_control.control.chmod(0o755)
            reset()
            require_rejected(
                lambda: BASE.build_closure_identity(wrong_control),
                "distinct production candidate/control paths were accepted")

            def bad_replay(*arguments: Any, **keywords: Any) -> dict[str, Any]:
                del arguments, keywords
                return {"schema": "unvalidated-replay"}

            BUILD_PROVENANCE.verify_reproducible_candidate_build = bad_replay
            reset()
            require_rejected(
                lambda: BASE.build_closure_identity(options),
                "invalid reproducible-build proof was accepted")
        finally:
            BUILD_PROVENANCE.candidate_build_provenance = saved_capture
            BUILD_PROVENANCE.verify_reproducible_candidate_build = \
                saved_replay
            BUILD_PROVENANCE.validate_reproducible_build_proof = \
                saved_validate
            PRODUCTION._BUILD_CLOSURE_BASELINE = None
            reset_shared_freezer()


def metric(value: float, iterations: int) -> dict[str, Any]:
    samples = [value] * iterations
    return {
        "median_us_per_batch_call": value,
        "mad_us_per_batch_call": 0.0,
        "minimum_us_per_batch_call": value,
        "maximum_us_per_batch_call": value,
        "samples_us_per_batch_call": samples,
    }


def synthetic_v31_result(
    implementation: str, cell: Mapping[str, Any], commit: str, tree: str,
    iterations: int = RUNNER.FIXED_ITERATIONS,
    warmup: int = RUNNER.FIXED_WARMUP,
) -> dict[str, Any]:
    mode = 0 if implementation == "control" else 1
    selected = (mode == 1 and cell["K"] == 65 and cell["R"] == 65 and
                cell["bytes"] == 64)
    return {
        "schema": "leopard2-benchmark-v31",
        "parameters": {
            "K": cell["K"], "R": cell["R"],
            "requested_profile": "legacy_high_v1",
            "requested_field": "gf8",
            "requested_backend": "avx2",
            "force_generic_decode": False,
            "force_specialized_decode": False,
            "force_tiled_decode": False,
            "force_materialized_decode": False,
            "skip_legacy": True,
            "retain_samples": True,
            "attest_source": True,
            "measure_one_shot_encode": True,
            "shard_bytes": cell["bytes"], "loss_count": cell["loss"],
            "missing_original_indices":
                RUNNER.expected_missing_original_indices(cell),
            "batch": cell["batch"], "reuse": cell["reuse"],
            "iterations": iterations, "warmup": warmup,
            "thread_count": 1, "seed": cell["seed"],
            "measure_one_shot_encode": True,
            "k65r65_b64_packed_terminal_mode": mode,
        },
        "resolved": {
            "profile": "legacy_high_v1", "field": "gf8",
            "backend": "avx2", "thread_count": 1,
            "parent_count": RUNNER.expected_codec_geometry(cell)[0],
            "padded_side": RUNNER.expected_codec_geometry(cell)[1],
        },
        "correctness": {
            "leopard2_round_trip": True, "legacy_comparison": None,
        },
        "workload_digests": {
            "algorithm": "fnv1a64",
            "original_data": "1" * 16,
            "transmitted_parity": "2" * 16,
            "recovered_originals": "1" * 16,
        },
        "build": {
            "compiler": "GNU",
            "compiler_version": "fixture",
            "cplusplus": 201703,
            "k8r3r4_t4_terminal_diagnostic_disabled": False,
            "t8_full_parity_terminal_diagnostic_disabled": False,
            "k16r8_b256_terminal_diagnostic_disabled": False,
            "k9r5_b256_terminal_diagnostic_disabled": False,
            "k9r5_b1024_terminal_diagnostic_disabled": False,
            "high_t16_q2_b64_fused_diagnostic_disabled": False,
            "source_commit": commit,
            "source_tree": tree,
            "source_tracked_dirty": False,
            "k65r65_b64_packed_terminal_diagnostic_mode": mode,
            "k65r65_b64_packed_terminal_diagnostic_disabled": mode == 0,
            "k65r65_b64_packed_terminal_mode_latched": mode,
            "k65r65_b64_packed_terminal_selector_expected_selected":
                selected,
            "k65r65_b64_packed_terminal_selector_selected": selected,
            "k65r65_b64_packed_terminal_selector_contract":
                RUNNER.SELECTOR_CONTRACT,
            "k65r65_b64_packed_terminal_timed_ordinary_encode_api":
                RUNNER.ORDINARY_API_MARKER,
            "k65r65_b64_packed_terminal_timed_one_shot_encode_api":
                RUNNER.ONE_SHOT_API_MARKER,
        },
        "metrics": {
            "codec_setup": {},
            "encode_execution": metric(1.0, iterations),
            "one_shot_encode": metric(1.1, iterations),
            "decode_plan_setup": {},
            "decode_execution": {},
            "decode_amortized_at_reuse": {},
            "rate_semantics": "fixture",
        },
        "memory": {
            "scratch_alignment": 64,
            "encode_scratch_bytes_per_stripe": 1,
            "decode_scratch_bytes_per_stripe": 1,
            "encode_scratch_bytes_batch": 1,
            "decode_scratch_bytes_batch": 1,
        },
        "legacy": {
            "available": False,
            "unavailable_reason": "skipped",
            "codec_setup": None,
            "decode_timing_includes_setup": True,
            "encode_execution": None,
            "decode_including_setup": None,
        },
    }


def exercise_v31_payload_validation(cells: list[dict[str, Any]]) -> None:
    commit = "e" * 40
    tree = "f" * 40
    target = cells[0]
    neighbor = cells[1]
    for implementation, cell, expected_selected in (
        ("candidate", target, True),
        ("control", target, False),
        ("candidate", neighbor, False),
        ("control", neighbor, False),
    ):
        document = synthetic_v31_result(
            implementation, cell, commit, tree)
        normalized = BASE.validate_result(
            implementation, document, cell, commit, tree,
            RUNNER.FIXED_ITERATIONS, RUNNER.FIXED_WARMUP)
        attribution = normalized[
            "k65r65_b64_packed_terminal_attribution"]
        expected_mode = 0 if implementation == "control" else 1
        require(
            normalized["schema"] == "leopard2-benchmark-v31" and
            attribution["requested_mode"] == expected_mode and
            attribution["latched_mode"] == expected_mode and
            attribution["selector_expected_selected"] is expected_selected and
            attribution["selector_selected"] is expected_selected and
            attribution["ordinary_api"] == RUNNER.ORDINARY_API_MARKER and
            attribution["one_shot_api"] == RUNNER.ONE_SHOT_API_MARKER,
            "v31 normalized selector/API attribution changed")

    pristine = synthetic_v31_result("candidate", target, commit, tree)
    mutations = (
        (lambda value: value.update({"schema": "leopard2-benchmark-v30"}),
         "non-v31 benchmark schema was accepted"),
        (lambda value: value["parameters"].update(
            {"k65r65_b64_packed_terminal_mode": 0}),
         "mismatched requested selector mode was accepted"),
        (lambda value: value["build"].update(
            {"k65r65_b64_packed_terminal_diagnostic_mode": False}),
         "boolean diagnostic mode was accepted as integer mode"),
        (lambda value: value["build"].update(
            {"k65r65_b64_packed_terminal_diagnostic_disabled": True}),
         "mismatched disabled marker was accepted"),
        (lambda value: value["build"].update(
            {"k65r65_b64_packed_terminal_mode_latched": 0}),
         "mismatched latched mode was accepted"),
        (lambda value: value["build"].update({
            "k65r65_b64_packed_terminal_selector_expected_selected": False}),
         "mismatched expected-selected marker was accepted"),
        (lambda value: value["build"].update(
            {"k65r65_b64_packed_terminal_selector_selected": False}),
         "mismatched selected marker was accepted"),
        (lambda value: value["build"].update(
            {"k65r65_b64_packed_terminal_selector_contract": "other"}),
         "mismatched selector contract was accepted"),
        (lambda value: value["build"].pop(
            "k65r65_b64_packed_terminal_timed_ordinary_encode_api"),
         "missing ordinary API marker was accepted"),
        (lambda value: value["build"].update({
            "k65r65_b64_packed_terminal_timed_one_shot_encode_api":
                "other"}),
         "mismatched one-shot API marker was accepted"),
        (lambda value: value["build"].update(
            {"source_tracked_dirty": True}),
         "dirty source attestation was accepted"),
        (lambda value: value["build"].update(
            {"source_commit": "d" * 40}),
         "mismatched embedded source commit was accepted"),
        (lambda value: value["workload_digests"].pop(
            "transmitted_parity"),
         "incomplete workload digests were accepted"),
        (lambda value: value["parameters"].update({"extra": False}),
         "extra v31 parameter key was accepted"),
        (lambda value: value["build"].pop("compiler_version"),
         "missing v31 build key was accepted"),
        (lambda value: value["metrics"].update({"extra": {}}),
         "extra v31 metric key was accepted"),
    )
    for mutate, message in mutations:
        malformed = deepcopy(pristine)
        mutate(malformed)
        require_rejected(
            lambda malformed=malformed: BASE.validate_result(
                "candidate", malformed, target, commit, tree,
                RUNNER.FIXED_ITERATIONS, RUNNER.FIXED_WARMUP),
            message)

    exact_parameter_values = {
        "K": target["K"],
        "R": target["R"],
        "shard_bytes": target["bytes"],
        "loss_count": target["loss"],
        "batch": target["batch"],
        "reuse": target["reuse"],
        "iterations": RUNNER.FIXED_ITERATIONS,
        "warmup": RUNNER.FIXED_WARMUP,
        "thread_count": 1,
        "seed": target["seed"],
    }
    for name, expected in exact_parameter_values.items():
        malformed = deepcopy(pristine)
        malformed["parameters"][name] = expected + 1
        require_rejected(
            lambda malformed=malformed: BASE.validate_result(
                "candidate", malformed, target, commit, tree,
                RUNNER.FIXED_ITERATIONS, RUNNER.FIXED_WARMUP),
            f"wrong exact {name} value was accepted")
    for name in ("loss_count", "batch", "thread_count"):
        malformed = deepcopy(pristine)
        malformed["parameters"][name] = True
        require_rejected(
            lambda malformed=malformed: BASE.validate_result(
                "candidate", malformed, target, commit, tree,
                RUNNER.FIXED_ITERATIONS, RUNNER.FIXED_WARMUP),
            f"boolean {name} was accepted as integer")

    missing = RUNNER.expected_missing_original_indices(target)
    malformed = deepcopy(pristine)
    malformed["parameters"]["missing_original_indices"] = [
        (missing[0] + 1) % target["K"]]
    require_rejected(
        lambda: BASE.validate_result(
            "candidate", malformed, target, commit, tree,
            RUNNER.FIXED_ITERATIONS, RUNNER.FIXED_WARMUP),
        "wrong deterministic missing-original selection was accepted")
    float_missing = deepcopy(pristine)
    float_missing["parameters"]["missing_original_indices"] = [
        float(missing[0])]
    require_rejected(
        lambda: BASE.validate_result(
            "candidate", float_missing, target, commit, tree,
            RUNNER.FIXED_ITERATIONS, RUNNER.FIXED_WARMUP),
        "float missing-original index was accepted as an exact integer")
    zero_missing_cell = cells[3]
    require(RUNNER.expected_missing_original_indices(zero_missing_cell) == [0],
            "Boolean missing-index fixture no longer aliases integer zero")
    bool_missing = synthetic_v31_result(
        "candidate", zero_missing_cell, commit, tree)
    bool_missing["parameters"]["missing_original_indices"] = [False]
    require_rejected(
        lambda: BASE.validate_result(
            "candidate", bool_missing, zero_missing_cell, commit, tree,
            RUNNER.FIXED_ITERATIONS, RUNNER.FIXED_WARMUP),
        "Boolean missing-original index was accepted as integer zero")

    parent_count, padded_side = RUNNER.expected_codec_geometry(target)
    for name, expected in (
            ("parent_count", parent_count),
            ("padded_side", padded_side)):
        malformed = deepcopy(pristine)
        malformed["resolved"][name] = expected + 1
        require_rejected(
            lambda malformed=malformed: BASE.validate_result(
                "candidate", malformed, target, commit, tree,
                RUNNER.FIXED_ITERATIONS, RUNNER.FIXED_WARMUP),
            f"wrong derived resolved {name} was accepted")

    unrelated_disables = (
        "k8r3r4_t4_terminal_diagnostic_disabled",
        "t8_full_parity_terminal_diagnostic_disabled",
        "k16r8_b256_terminal_diagnostic_disabled",
        "k9r5_b256_terminal_diagnostic_disabled",
        "k9r5_b1024_terminal_diagnostic_disabled",
        "high_t16_q2_b64_fused_diagnostic_disabled",
    )
    for name in unrelated_disables:
        malformed = deepcopy(pristine)
        malformed["build"][name] = True
        require_rejected(
            lambda malformed=malformed: BASE.validate_result(
                "candidate", malformed, target, commit, tree,
                RUNNER.FIXED_ITERATIONS, RUNNER.FIXED_WARMUP),
            f"active unrelated disable marker {name} was accepted")


def synthetic_round(
    order: tuple[str, ...], *, control: float = 1.0,
    main: float = 1.0, one_shot_control: float | None = None,
    one_shot_main: float | None = None,
) -> dict[str, Any]:
    execution = {"candidate": 1.0, "control": control, "main": main}
    one_shot = {
        "candidate": 1.0,
        "control": control if one_shot_control is None else one_shot_control,
        "main": main if one_shot_main is None else one_shot_main,
    }
    digests = {
        "original_data": "1" * 16,
        "transmitted_parity": "2" * 16,
        "recovered_originals": "1" * 16,
    }
    return {
        "invocations": [{
            "implementation": label,
            "normalized": {
                "encode_us": execution[label],
                "one_shot_encode_us": one_shot[label],
                "digests": dict(digests),
            },
        } for label in order],
        "isolation": {"accepted": True},
    }


def exercise_analysis(cells: list[dict[str, Any]]) -> None:
    target_orders = BASE.select_round_orders(
        BASE.TARGET_ORDER, RUNNER.CONFIRMATORY_ROUNDS)
    neighbor_orders = BASE.select_round_orders(
        BASE.NEIGHBOR_ORDER, RUNNER.CONFIRMATORY_ROUNDS)
    require(len(target_orders) == 25 and len(neighbor_orders) == 25 and
            target_orders[:3] == BASE.TARGET_ORDER and
            neighbor_orders[:4] == BASE.NEIGHBOR_ORDER,
            "fixed 25-round ABBA order changed")
    require_rejected(
        lambda: BASE.select_round_orders(BASE.TARGET_ORDER, 9),
        "non-confirmatory round count was accepted")

    log_ratios = [(index - 12) / 1000.0 for index in range(25)]
    interval = BASE.confidence_interval(log_ratios)
    expected_half_width = PRODUCTION.T95_DF24 * statistics.stdev(
        log_ratios) / math.sqrt(25)
    require(
        math.isclose(interval["speedup"], 1.0) and
        math.isclose(interval["ci95"][0], math.exp(-expected_half_width)) and
        math.isclose(interval["ci95"][1], math.exp(expected_half_width)) and
        interval["confidence_level"] == 0.95 and
        interval["degrees_of_freedom"] == 24 and
        interval["t_critical"] == PRODUCTION.T95_DF24,
        "25-round Student-t df=24 interval changed")
    require_rejected(
        lambda: BASE.confidence_interval([0.0] * 9),
        "non-df24 confidence interval was accepted")

    target_rounds = [
        synthetic_round(order, control=1.06, main=1.06)
        for order in target_orders]
    target = BASE.analyze(cells[0], target_rounds)
    require(
        target["target_acceptance_validation"]["accepted"] is True and
        target["target_acceptance_validation"]["control_floor"] == 1.05 and
        target["target_acceptance_validation"]["main_floor"] == 1.05,
        "target batch/one-shot acceptance gate changed")
    target_regression = [
        synthetic_round(order, control=1.04, main=1.06)
        for order in target_orders]
    require_rejected(
        lambda: BASE.analyze(cells[0], target_regression),
        "target below the mature-control floor was accepted")
    main_regression = [
        synthetic_round(order, control=1.06, main=1.04)
        for order in target_orders]
    require_rejected(
        lambda: BASE.analyze(cells[0], main_regression),
        "target below the exact-main floor was accepted")

    immediate = cells[1]
    inert_rounds = [synthetic_round(order) for order in neighbor_orders]
    inert = BASE.analyze(immediate, inert_rounds)
    require(
        inert["neighbor_selector_inertness_validation"]["accepted"] is True
        and "main_over_candidate" not in inert,
        "immediate neighbor selector-inertness classification changed")
    mismatched_digests = deepcopy(inert_rounds)
    mismatched_digests[0]["invocations"][0]["normalized"]["digests"][
        "transmitted_parity"] = "9" * 16
    require_rejected(
        lambda: BASE.analyze(immediate, mismatched_digests),
        "cross-mode workload digest mismatch was accepted")
    leaked = [
        synthetic_round(order, control=1.021)
        for order in neighbor_orders]
    require_rejected(
        lambda: BASE.analyze(immediate, leaked),
        "full CI outside the two-percent equivalence band was accepted")
    leaked_one_shot = [
        synthetic_round(order, one_shot_control=0.979)
        for order in neighbor_orders]
    require_rejected(
        lambda: BASE.analyze(immediate, leaked_one_shot),
        "one-shot CI outside the equivalence band was accepted")

    retained = cells[7]
    retained_rounds = [
        synthetic_round(order, main=0.99)
        for order in target_orders]
    retained_analysis = BASE.analyze(retained, retained_rounds)
    require(
        retained_analysis["exact_main_context"]["policy"] == "gated" and
        retained_analysis["exact_main_context"]["floor"] == 0.98 and
        retained_analysis["exact_main_context"]["accepted"] is True,
        "retained B8192 exact-main gate changed")
    retained_regression = [
        synthetic_round(order, main=0.97)
        for order in target_orders]
    require_rejected(
        lambda: BASE.analyze(retained, retained_regression),
        "retained B8192 exact-main regression was accepted")

    for balanced in cells[8:10]:
        accepted = BASE.analyze(
            balanced,
            [synthetic_round(order, main=0.99) for order in target_orders])
        require(accepted["exact_main_context"]["policy"] == "gated",
                "balanced layout canary lost its established main gate")
    context_cell = cells[10]
    context = BASE.analyze(
        context_cell,
        [synthetic_round(order, main=0.50) for order in target_orders])
    require(
        context["exact_main_context"]["policy"] == "context_only" and
        context["exact_main_context"]["floor"] is None and
        context["exact_main_context"]["metrics"]
            ["encode_execution"]["ci95"][0] == 0.5,
        "unestablished layout main context became a release gate")


def exercise_arguments_and_commands(cells: list[dict[str, Any]]) -> None:
    runner_path = Path(RUNNER.__file__).resolve()
    digest = "a" * 64
    arguments = [
        str(runner_path),
        "--candidate", "/build/bench_leopard2",
        "--candidate-sha256", digest,
        "--control", "/build/bench_leopard2",
        "--control-sha256", digest,
        "--main", "/frozen/main",
        "--main-sha256", BASE.CANONICAL_MAIN_SHA256,
        "--source-commit", "e" * 40,
        "--source-tree", "f" * 40,
        "--output", "/results/k65r65-b64",
        "--cpu", "13", "--sibling", "29",
    ]
    saved_argv = sys.argv
    try:
        sys.argv = arguments
        parsed = BASE.parse_arguments()
        require(
            parsed.candidate == Path("/build/bench_leopard2") and
            parsed.control == parsed.candidate and
            parsed.iterations == RUNNER.FIXED_ITERATIONS and
            parsed.warmup == RUNNER.FIXED_WARMUP and
            parsed.rounds == RUNNER.CONFIRMATORY_ROUNDS and
            parsed.candidate_archive is None and
            parsed.candidate_compile_commands is None,
            "fixed authoritative command-line policy changed")

        sys.argv = arguments + ["--iterations", "15"]
        require_argument_rejected(
            BASE.parse_arguments,
            "timing iteration escape hatch was accepted")
        wrong_main = list(arguments)
        wrong_main[wrong_main.index("--main-sha256") + 1] = "b" * 64
        sys.argv = wrong_main
        require_rejected(
            BASE.parse_arguments,
            "noncanonical exact-main SHA-256 was accepted")
        wrong_control_hash = list(arguments)
        wrong_control_hash[
            wrong_control_hash.index("--control-sha256") + 1] = "b" * 64
        sys.argv = wrong_control_hash
        require_rejected(
            BASE.parse_arguments,
            "different candidate/control SHA-256 values were accepted")
    finally:
        sys.argv = saved_argv

    target = cells[0]
    candidate = BASE.benchmark_command(
        "candidate", Path("/frozen/candidate"), target, 13,
        RUNNER.FIXED_ITERATIONS, RUNNER.FIXED_WARMUP)
    control = BASE.benchmark_command(
        "control", Path("/frozen/candidate"), target, 13,
        RUNNER.FIXED_ITERATIONS, RUNNER.FIXED_WARMUP)
    expected_prefix = [
        "/usr/bin/prlimit", "--as=201326592",
        "/usr/bin/taskset", "-c", "13", "/frozen/candidate",
        "--k", "65", "--r", "65", "--bytes", "64",
        "--loss", "1", "--batch", "1", "--reuse", "8192",
        "--iterations", "31", "--warmup", "64", "--threads", "1",
        "--seed", str(target["seed"]),
        "--profile", "high", "--field", "gf8", "--backend", "avx2",
        "--skip-legacy", "--retain-samples", "--attest-source",
        "--measure-one-shot-encode",
    ]
    require(
        candidate == expected_prefix + [
            RUNNER.SELECTOR_ARGUMENT, "1", "--json", "-"] and
        control == expected_prefix + [
            RUNNER.SELECTOR_ARGUMENT, "0", "--json", "-"] and
        candidate[5] == control[5],
        "candidate/control exact benchmark commands changed")

    byte_neighbor = cells[5]
    main = BASE.benchmark_command(
        "main", Path("/frozen/main"), byte_neighbor, 13,
        RUNNER.FIXED_ITERATIONS, RUNNER.FIXED_WARMUP)
    require(
        main[main.index("--bytes") + 1] == "64" and
        main[main.index("--logical-bytes") + 1] == "63" and
        RUNNER.SELECTOR_ARGUMENT not in main,
        "canonical-main logical-byte command adapter changed")


def exercise_workload_sized_child_timeout(
    cells: list[dict[str, Any]],
) -> None:
    target = cells[0]
    retained = cells[7]
    candidate_budget = RUNNER.child_timeout_budget(
        "candidate", retained, RUNNER.FIXED_ITERATIONS,
        RUNNER.FIXED_WARMUP)
    main_budget = RUNNER.child_timeout_budget(
        "main", retained, RUNNER.FIXED_ITERATIONS,
        RUNNER.FIXED_WARMUP)
    target_budget = RUNNER.child_timeout_budget(
        "candidate", target, RUNNER.FIXED_ITERATIONS,
        RUNNER.FIXED_WARMUP)
    expected_calls = (
        RUNNER.FIXED_ITERATIONS * RUNNER.FIXED_REUSE +
        RUNNER.FIXED_WARMUP)
    expected_bytes_per_call = (retained["K"] + retained["R"]) * 8192
    expected_visits = 3 * expected_calls * expected_bytes_per_call
    expected_seconds = RUNNER.CHILD_TIMEOUT_SETUP_SECONDS + math.ceil(
        expected_visits / RUNNER.CHILD_TIMEOUT_LOGICAL_BYTES_PER_SECOND)
    require(
        candidate_budget == {
            "timeout_seconds": expected_seconds,
            "setup_seconds": RUNNER.CHILD_TIMEOUT_SETUP_SECONDS,
            "logical_bytes_per_second_floor":
                RUNNER.CHILD_TIMEOUT_LOGICAL_BYTES_PER_SECOND,
            "measured_metric_count": 3,
            "calls_per_metric": expected_calls,
            "logical_bytes_per_call": expected_bytes_per_call,
            "logical_byte_visits": expected_visits,
        } and
        candidate_budget["timeout_seconds"] > 60 and
        candidate_budget["timeout_seconds"] >
            target_budget["timeout_seconds"] and
        main_budget["measured_metric_count"] == 2 and
        main_budget["timeout_seconds"] > 60,
        "B8192 fixed-reuse workload timeout budget changed")
    require(
        K66._BASE_RUN_ONE is RUNNER.run_one_with_workload_timeout,
        "shared-inode runner did not install the workload-sized inner child")

    with tempfile.TemporaryDirectory(
            prefix="leo2-k65r65-timeout-launch-") as directory:
        root = Path(directory)
        executable = root / "candidate"
        executable.write_bytes(b"timeout launcher fixture\n")
        executable.chmod(0o555)
        executable_identity = identity(executable)
        captured: dict[str, Any] = {}

        class Completed:
            returncode = 0
            stdout = b"{}"
            stderr = b""

        def fake_run(command: Any, **keywords: Any) -> Completed:
            captured["command"] = list(command)
            captured["keywords"] = dict(keywords)
            return Completed()

        saved_run = RUNNER.subprocess.run
        saved_validate = BASE.validate_result
        try:
            RUNNER.subprocess.run = fake_run
            BASE.validate_result = lambda *args, **kwargs: {
                "schema": "mock-v31", "digests": {}}
            invocation = RUNNER.run_one_with_workload_timeout(
                "candidate", executable_identity, retained, 13,
                "e" * 40, "f" * 40, RUNNER.FIXED_ITERATIONS,
                RUNNER.FIXED_WARMUP, root)
        finally:
            RUNNER.subprocess.run = saved_run
            BASE.validate_result = saved_validate
        require(
            captured["keywords"]["timeout"] ==
                candidate_budget["timeout_seconds"] and
            captured["keywords"]["check"] is False and
            captured["command"] == invocation["command"] and
            invocation["timeout_budget"] == candidate_budget and
            invocation["command"][
                invocation["command"].index(RUNNER.SELECTOR_ARGUMENT) + 1]
                == "1",
            "inner child launch did not use its workload-sized timeout")


def main() -> int:
    cells = RUNNER.cells()
    require(
        [(cell["K"], cell["R"], cell["bytes"]) for cell in cells] == [
            (65, 65, 64),
            (64, 65, 64), (66, 65, 64),
            (65, 64, 64), (65, 66, 64),
            (65, 65, 63), (65, 65, 65),
            (65, 65, 8192),
            (64, 64, 64), (128, 128, 64),
            (79, 32, 64), (62, 8, 64), (66, 16, 64),
        ],
        "K65/R65 target, boundaries, retained path, or canaries changed")
    require(
        [cell["seed"] for cell in cells] ==
            list(range(0x6565B640, 0x6565B64D)) and
        all(cell["batch"] == 1 and cell["reuse"] == 8192 and
            cell["loss"] == 1 and cell["measure_one_shot"] is True
            for cell in cells),
        "cell seed or fixed timing payload changed")
    require(
        [RUNNER.expected_codec_geometry(cell) for cell in cells] == [
            (256, 128), (256, 128), (256, 128), (256, 64),
            (256, 128), (256, 128), (256, 128), (256, 128),
            (128, 64), (256, 128), (128, 32), (128, 8), (128, 16),
        ] and
        [RUNNER.expected_missing_original_indices(cell) for cell in cells] == [
            [28], [47], [11], [0], [2], [41], [24],
            [17], [43], [66], [18], [22], [31],
        ],
        "derived codec geometry or deterministic loss selection changed")
    require(
        BASE.SCHEMA ==
            "leopard2-gf8-k65r65-b64-packed-terminal-abba/v1" and
        BASE.SUMMARY_SCHEMA ==
            "leopard2-gf8-k65r65-b64-packed-terminal-summary/v1" and
        BASE.MODE_SYMBOL ==
            "_ZN12_GLOBAL__N_1L33g_k65r65_b64_packed_terminal_modeE" and
        BASE.CANDIDATE_SCHEMA == "leopard2-benchmark-v31" and
        BASE.CONTROL_SCHEMA == "leopard2-benchmark-v31" and
        BASE.CONTROL_BUILD_MARKER ==
            "k65r65_b64_packed_terminal_diagnostic_disabled" and
        BASE.ALLOW_IDENTICAL_CANDIDATE_CONTROL is True and
        BASE.REQUIRE_EXPECTED_IDENTITIES is True and
        BASE.REQUIRE_BUILD_CLOSURE is True and
        BASE.REQUIRE_FULL_ELF_IDENTITY is True and
        BASE.TARGET_CONTROL_FLOOR == 1.05 and
        BASE.TARGET_MAIN_FLOOR == 1.05 and
        BASE.NEIGHBOR_EQUIVALENCE_LOWER == 1.0 / 1.02 and
        BASE.NEIGHBOR_EQUIVALENCE_UPPER == 1.02 and
        BASE.RETAINED_MAIN_FLOOR == 0.98 and
        BASE.MAX_ISOLATION_ATTEMPTS == 3,
        "attribution, schema, retry, or promotion policy changed")
    require(
        BASE.CANONICAL_MAIN_SHA256 ==
            "78575e2fe1b9796d10dcb4c14d2bcc1b627a79e787be66d34dc7aa666899aa93",
        "canonical exact-main executable identity changed")
    require(
        BASE.RUNNER_DEPENDENCIES[0] == Path(RUNNER.__file__).resolve() and
        BASE.RUNNER_DEPENDENCIES[1:] ==
            RUNNER._INHERITED_RUNNER_DEPENDENCIES and
        len(BASE.RUNNER_DEPENDENCIES) ==
            len(set(BASE.RUNNER_DEPENDENCIES)),
        "runner/source support attestation closure changed")

    exercise_arguments_and_commands(cells)
    exercise_workload_sized_child_timeout(cells)
    exercise_v31_payload_validation(cells)
    exercise_analysis(cells)
    exercise_shared_inode_freezer()
    exercise_production_build_closure()
    print("K65/R65/B64 hardened ABBA runner checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
