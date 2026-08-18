#!/usr/bin/env python3
"""Deterministic contract checks for the T32/Q3/B64 boundary runner."""

from __future__ import annotations

import argparse
import importlib.util
import math
import sys
import tempfile
from copy import deepcopy
from pathlib import Path
from typing import Any


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_runner() -> Any:
    path = Path(__file__).resolve().with_name(
        "run_t32_q3_b64_boundary_abba.py")
    specification = importlib.util.spec_from_file_location(
        "leopard2_t32_q3_b64_boundary_runner_test_target", path)
    require(specification is not None and specification.loader is not None,
            "cannot load T32/Q3/B64 boundary runner")
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


def identity(path: Path) -> dict[str, Any]:
    return BASE.T8_SUPPORT.regular_file_identity(path)


def synthetic_round(
    order: tuple[str, ...], *, control_execution: float = 1.10,
    control_one_shot: float = 1.10, main_execution: float = 1.01,
    main_one_shot: float = 1.01, digest: str = "2" * 16,
) -> dict[str, Any]:
    digests = {
        "original_data": "1" * 16,
        "transmitted_parity": digest,
        "recovered_originals": "1" * 16,
    }
    execution = {
        "candidate": 1.0,
        "control": control_execution,
        "main": main_execution,
    }
    one_shot = {
        "candidate": 1.0,
        "control": control_one_shot,
        "main": main_one_shot,
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
    archive_identity = identity(archive_object)
    benchmark_identity = identity(benchmark_object)

    dependencies = []
    for dependency in BASE.RUNNER_DEPENDENCIES:
        source_identity = BASE.support_file_identity(dependency)
        dependencies.append({
            "path": dependency.resolve().relative_to(
                PROVENANCE_RUNNER.SOURCE_ROOT).as_posix(),
            "sha256": source_identity["sha256"],
            "size": source_identity["size"],
            "mode": 0o644,
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
                "source": identity(
                    PROVENANCE_RUNNER.SOURCE_ROOT / "leopard2.cpp"),
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


def reset_shared_freeze() -> None:
    RUNNER.PARENT._SHARED_FROZEN_CANDIDATE = None
    RUNNER.PARENT._SHARED_FROZEN_PHYSICAL_IDENTITY = None


def exercise_shared_inode_contract() -> None:
    with tempfile.TemporaryDirectory(
            prefix="leo2-t32-q3-b64-shared-") as directory:
        root = Path(directory)
        executable = root / "bench_leopard2"
        executable.write_bytes(b"shared candidate and control\n")
        executable.chmod(0o755)
        digest = identity(executable)["sha256"]
        frozen = root / "frozen"
        frozen.mkdir()
        reset_shared_freeze()
        candidate = BASE.freeze_executable(
            executable, digest, frozen / "candidate")
        control = BASE.freeze_executable(
            executable, digest, frozen / "control")
        require(
            candidate["input"] == control["input"] and
            candidate["frozen"] == control["frozen"] and
            candidate["frozen_physical_identity"] ==
                control["frozen_physical_identity"] and
            candidate["frozen_physical_identity"]["inode"] ==
                (frozen / "candidate").stat().st_ino and
            not (frozen / "control").exists(),
            "candidate/control did not retain one frozen inode")

        other = root / "same-bytes-other-inode"
        other.write_bytes(executable.read_bytes())
        other.chmod(0o755)
        second = root / "second"
        second.mkdir()
        reset_shared_freeze()
        BASE.freeze_executable(executable, digest, second / "candidate")
        require_rejected(
            lambda: BASE.freeze_executable(other, digest, second / "control"),
            "same-byte distinct input inodes were accepted")

        replaced = root / "replaced"
        replaced.mkdir()
        reset_shared_freeze()
        record = BASE.freeze_executable(
            executable, digest, replaced / "candidate")
        path = Path(record["frozen"]["path"])
        replacement = replaced / "replacement"
        replacement.write_bytes(path.read_bytes())
        replacement.chmod(0o555)
        replacement.replace(path)
        require_rejected(
            RUNNER.PARENT.require_shared_frozen_physical_identity,
            "same-byte frozen executable inode replacement was accepted")


def exercise_production_replay_contract() -> None:
    with tempfile.TemporaryDirectory(
            prefix="leo2-t32-q3-b64-provenance-") as directory:
        build = Path(directory) / "build"
        build.mkdir()
        executable = build / "bench_leopard2"
        executable.write_bytes(b"fixture executable\n")
        executable.chmod(0o755)
        digest = identity(executable)["sha256"]
        commit = "e" * 40
        tree = "f" * 40
        pristine = fake_production_provenance(
            build, executable, commit, tree)
        options = argparse.Namespace(
            candidate=executable,
            control=executable,
            candidate_sha256=digest,
            control_sha256=digest,
            source_commit=commit,
            source_tree=tree,
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
        calls: list[tuple[Any, ...]] = []

        def capture(*arguments: Any, **keywords: Any) -> dict[str, Any]:
            calls.append(("capture", arguments, keywords))
            return deepcopy(current["value"])

        def replay(*arguments: Any, **keywords: Any) -> dict[str, Any]:
            calls.append(("replay", arguments, keywords))
            return {"schema": "mock-production-replay", "valid": True}

        def validate(
            proof: Any, candidate: Any, **keywords: Any,
        ) -> None:
            calls.append(("validate", proof, candidate, keywords))
            if proof != {"schema": "mock-production-replay", "valid": True}:
                raise RuntimeError("mock replay proof rejected")

        def reset(value: dict[str, Any] | None = None) -> None:
            PROVENANCE_RUNNER._BUILD_CLOSURE_BASELINE = None
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
                    {"schema": "mock-production-replay", "valid": True} and
                len(closure["runner_source_dependencies"]) ==
                    len(BASE.RUNNER_DEPENDENCIES),
                "production capture/replay/source closure changed")
            require(
                calls[0] == (
                    "capture",
                    (build.resolve(), PROVENANCE_RUNNER.SOURCE_ROOT,
                     executable.resolve(), "bench_leopard2"),
                    {}) and
                any(call[0] == "replay" and call[2] == {"jobs": 1}
                    for call in calls) and
                any(call[0] == "validate" for call in calls),
                "exact production target was not captured and replayed")
            BASE.build_closure_identity(options)
            require(sum(call[0] == "capture" for call in calls) == 2 and
                    sum(call[0] == "replay" for call in calls) == 1 and
                    sum(call[0] == "validate" for call in calls) == 2,
                    "post-campaign closure did not recapture and revalidate")

            malformed = deepcopy(pristine)
            malformed["tracked_source_manifest"]["files"].pop()
            reset(malformed)
            require_rejected(
                lambda: BASE.build_closure_identity(options),
                "missing runner dependency was accepted")

            malformed = deepcopy(pristine)
            malformed["source_object_compile_closure"] = []
            reset(malformed)
            require_rejected(
                lambda: BASE.build_closure_identity(options),
                "empty production compile closure was accepted")

            malformed = deepcopy(pristine)
            malformed["validated_cache"][
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA"] = "old"
            reset(malformed)
            require_rejected(
                lambda: BASE.build_closure_identity(options),
                "wrong production configuration schema was accepted")

            wrong_hash = argparse.Namespace(**vars(options))
            wrong_hash.candidate_sha256 = "0" * 64
            wrong_hash.control_sha256 = "0" * 64
            reset()
            require_rejected(
                lambda: BASE.build_closure_identity(wrong_hash),
                "wrong required production executable hash was accepted")

            other = build / "other-bench_leopard2"
            other.write_bytes(executable.read_bytes())
            other.chmod(0o755)
            distinct = argparse.Namespace(**vars(options))
            distinct.control = other
            reset()
            require_rejected(
                lambda: BASE.build_closure_identity(distinct),
                "distinct candidate/control paths were accepted")
        finally:
            BUILD_PROVENANCE.candidate_build_provenance = saved_capture
            BUILD_PROVENANCE.verify_reproducible_candidate_build = saved_replay
            BUILD_PROVENANCE.validate_reproducible_build_proof = saved_validate
            PROVENANCE_RUNNER._BUILD_CLOSURE_BASELINE = None


def exercise_result_mode_contract() -> None:
    saved = RUNNER._BASE_VALIDATE_RESULT
    fixture = {"normalized": True}

    def base_validate(*arguments: Any, **keywords: Any) -> dict[str, Any]:
        del arguments, keywords
        return fixture

    RUNNER._BASE_VALIDATE_RESULT = base_validate
    cell = RUNNER.cells()[0]
    try:
        for implementation, mode in (("candidate", 1), ("control", 0)):
            result = {
                "build": {
                    "balanced_b64_terminal_diagnostic_mode": mode,
                    "balanced_b64_terminal_enabled": mode == 1,
                },
            }
            require(
                RUNNER.validate_result(
                    implementation, result, cell, "e" * 40, "f" * 40,
                    31, 64) is fixture,
                f"{implementation} correct diagnostic mode was rejected")
            wrong = deepcopy(result)
            wrong["build"]["balanced_b64_terminal_diagnostic_mode"] = 1 - mode
            require_rejected(
                lambda implementation=implementation, wrong=wrong:
                    RUNNER.validate_result(
                        implementation, wrong, cell, "e" * 40, "f" * 40,
                        31, 64),
                f"{implementation} mislabeled diagnostic mode was accepted")
            wrong = deepcopy(result)
            wrong["build"]["balanced_b64_terminal_enabled"] = mode != 1
            require_rejected(
                lambda implementation=implementation, wrong=wrong:
                    RUNNER.validate_result(
                        implementation, wrong, cell, "e" * 40, "f" * 40,
                        31, 64),
                f"{implementation} inconsistent enabled marker was accepted")
    finally:
        RUNNER._BASE_VALIDATE_RESULT = saved


def main() -> int:
    cells = RUNNER.cells()
    require(
        [(cell["id"], cell["K"], cell["R"], cell["bytes"],
          cell["qualification"])
         for cell in cells] == [
            ("target-k77-r32-b64-q3", 77, 32, 64, "focus"),
            ("target-k78-r32-b64-q3", 78, 32, 64, "focus"),
            ("target-k79-r32-b64-q3", 79, 32, 64, "focus"),
            ("context-k65-r32-b64-q3", 65, 32, 64, "context"),
            ("context-k76-r32-b64-q3", 76, 32, 64, "context"),
            ("context-k80-r32-b64-q3", 80, 32, 64, "context"),
            ("context-k81-r32-b64-q3", 81, 32, 64, "context"),
            ("context-k95-r32-b64-q3", 95, 32, 64, "context"),
            ("context-k96-r32-b64-q3", 96, 32, 64, "context"),
            ("r-neighbor-k79-r31-b64-q3", 79, 31, 64, "inactive"),
            ("r-neighbor-k79-r33-b64-q2", 79, 33, 64, "inactive"),
            ("byte-neighbor-k79-r32-b63-q3", 79, 32, 63, "inactive"),
            ("byte-neighbor-k79-r32-b65-q3", 79, 32, 65, "inactive"),
            ("k-neighbor-k97-r32-b64-q4", 97, 32, 64, "inactive"),
        ] and
        [cell["seed"] for cell in cells] ==
            list(range(0x54334200, 0x5433420E)) and
        all(cell["batch"] == 1 and cell["reuse"] == 8192 and
            cell["loss"] == 1 and cell["measure_one_shot"] is True
            for cell in cells),
        "qualification dimensions, order, seeds, or call scope changed")
    require(
        BASE.SCHEMA == "leopard2-gf8-t32-q3-b64-boundary-abba/v1" and
        BASE.SUMMARY_SCHEMA ==
            "leopard2-gf8-t32-q3-b64-boundary-summary/v1" and
        BASE.MODE_SYMBOL ==
            "_ZN12_GLOBAL__N_1L28g_balanced_b64_terminal_modeE" and
        BASE.ALLOW_IDENTICAL_CANDIDATE_CONTROL is True and
        BASE.ALLOW_MULTIPLE_TARGETS is True and
        BASE.CANDIDATE_SCHEMA == BASE.CONTROL_SCHEMA ==
            "leopard2-benchmark-v22" and
        BASE.CONTROL_EXTRA_ARGUMENTS ==
            ("--balanced-b64-terminal-mode", "0") and
        BASE.CONTROL_BUILD_MARKER is None and
        BASE.TARGET_CONTROL_FLOOR == 1.05 and
        BASE.TARGET_MAIN_FLOOR == 1.0 and
        BASE.NEIGHBOR_EQUIVALENCE_LOWER == 1.0 / 1.02 and
        BASE.NEIGHBOR_EQUIVALENCE_UPPER == 1.02 and
        BASE.REQUIRE_EXPECTED_IDENTITIES is True and
        BASE.REQUIRE_BUILD_CLOSURE is True and
        BASE.REQUIRE_FULL_ELF_IDENTITY is True and
        BASE.CANONICAL_MAIN_SHA256 ==
            "78575e2fe1b9796d10dcb4c14d2bcc1b627a79e787be66d34dc7aa666899aa93",
        "qualification identity, schema, or acceptance policy changed")
    expected_dependencies = (
        Path(RUNNER.__file__).resolve(),
        *RUNNER.SUPPORT_DEPENDENCIES,
    )
    require(BASE.RUNNER_PATH == Path(RUNNER.__file__).resolve() and
            BASE.RUNNER_DEPENDENCIES == expected_dependencies and
            len(set(expected_dependencies)) == len(expected_dependencies),
            "runner source dependency attestation changed")

    digest = "a" * 64
    saved_argv = sys.argv
    try:
        sys.argv = [
            str(BASE.RUNNER_PATH),
            "--candidate", "/build/bench_leopard2",
            "--candidate-sha256", digest,
            "--control", "/build/bench_leopard2",
            "--control-sha256", digest,
            "--main", "/frozen/main",
            "--main-sha256", BASE.CANONICAL_MAIN_SHA256,
            "--source-commit", "e" * 40,
            "--source-tree", "f" * 40,
            "--output", "/results/t32-q3-b64",
            "--cpu", "30", "--sibling", "31",
        ]
        parsed = BASE.parse_arguments()
    finally:
        sys.argv = saved_argv
    require(
        parsed.candidate == parsed.control == Path("/build/bench_leopard2") and
        parsed.rounds == 25 and parsed.iterations == 31 and
        parsed.warmup == 64 and parsed.candidate_archive is None and
        parsed.control_archive is None and
        parsed.candidate_compile_commands is None and
        parsed.control_compile_commands is None,
        "25-round production command-line contract changed")
    saved_argv = sys.argv
    try:
        sys.argv = [
            str(BASE.RUNNER_PATH),
            "--candidate", "/build/bench_leopard2",
            "--candidate-sha256", digest,
            "--control", "/build/bench_leopard2",
            "--control-sha256", digest,
            "--main", "/frozen/main",
            "--main-sha256", "b" * 64,
            "--source-commit", "e" * 40,
            "--source-tree", "f" * 40,
            "--output", "/results/t32-q3-b64",
            "--cpu", "30", "--sibling", "31",
        ]
        require_rejected(
            BASE.parse_arguments,
            "noncanonical exact-main executable hash was accepted")
    finally:
        sys.argv = saved_argv

    target = cells[0]
    candidate_command = BASE.benchmark_command(
        "candidate", Path("/frozen/shared"), target, 30, 31, 64)
    control_command = BASE.benchmark_command(
        "control", Path("/frozen/shared"), target, 30, 31, 64)
    main_command = BASE.benchmark_command(
        "main", Path("/frozen/main"), target, 30, 31, 64)
    require(
        candidate_command.count("--balanced-b64-terminal-mode") == 1 and
        candidate_command[
            candidate_command.index("--balanced-b64-terminal-mode") + 1] ==
            "1" and
        control_command.count("--balanced-b64-terminal-mode") == 1 and
        control_command[
            control_command.index("--balanced-b64-terminal-mode") + 1] ==
            "0" and
        "--balanced-b64-terminal-mode" not in main_command and
        "--measure-one-shot-encode" in candidate_command and
        candidate_command[candidate_command.index("--batch") + 1] == "1" and
        candidate_command[candidate_command.index("--iterations") + 1] ==
            "31" and
        candidate_command[candidate_command.index("--warmup") + 1] == "64",
        "same-inode mode or one-item/one-shot benchmark command changed")

    orders = BASE.select_round_orders(BASE.TARGET_ORDER, 25)
    neighbor_orders = BASE.select_round_orders(BASE.NEIGHBOR_ORDER, 25)
    require(len(orders) == len(neighbor_orders) == 25 and
            orders[:3] == BASE.TARGET_ORDER and
            neighbor_orders[:4] == BASE.NEIGHBOR_ORDER,
            "25-round ABBA schedules changed")
    target_analysis = BASE.analyze(
        target, [synthetic_round(order) for order in orders])
    require(target_analysis["active_boundary_validation"]["accepted"] is True,
            "valid active boundary contrast was rejected")
    require_rejected(
        lambda: BASE.analyze(
            target,
            [synthetic_round(order, control_execution=1.049,
                             control_one_shot=1.10)
             for order in orders]),
        "active execution contrast below five percent was accepted")
    require_rejected(
        lambda: BASE.analyze(
            target,
            [synthetic_round(order, main_execution=0.999,
                             main_one_shot=1.01)
             for order in orders]),
        "active exact-main regression was accepted")

    neighbor = cells[-1]
    neighbor_rounds = [
        synthetic_round(order, control_execution=1.0,
                        control_one_shot=1.0)
        for order in neighbor_orders
    ]
    neighbor_analysis = BASE.analyze(neighbor, neighbor_rounds)
    require(
        neighbor_analysis["inactive_boundary_validation"]["accepted"] is True,
        "valid inactive selector equivalence was rejected")
    require_rejected(
        lambda: BASE.analyze(
            neighbor,
            [synthetic_round(order, control_execution=1.021,
                             control_one_shot=1.0)
             for order in neighbor_orders]),
        "inactive execution effect outside two percent was accepted")
    require_rejected(
        lambda: BASE.analyze(
            neighbor,
            [synthetic_round(order, control_execution=1.0,
                             control_one_shot=0.979)
             for order in neighbor_orders]),
        "inactive one-shot effect outside two percent was accepted")
    corrupted = deepcopy(neighbor_rounds)
    corrupted[0]["invocations"][0]["normalized"]["digests"][
        "transmitted_parity"] = "f" * 16
    require_rejected(
        lambda: BASE.analyze(neighbor, corrupted),
        "mismatched parity digest was accepted")

    exercise_result_mode_contract()
    exercise_shared_inode_contract()
    exercise_production_replay_contract()
    print("T32/Q3/B64 boundary hardened runner checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
