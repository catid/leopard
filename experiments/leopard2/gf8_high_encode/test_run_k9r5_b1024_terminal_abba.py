#!/usr/bin/env python3
"""Deterministic contract checks for the K9/R5/B1024 ABBA wrapper."""

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


def require_rejected(callback: Any, message: str) -> None:
    try:
        callback()
    except BASE.EvidenceError:
        return
    raise RuntimeError(message)


def load_runner() -> Any:
    path = Path(__file__).resolve().with_name(
        "run_k9r5_b1024_terminal_abba.py")
    specification = importlib.util.spec_from_file_location(
        "leopard2_k9r5_b1024_runner_test_target", path)
    require(specification is not None and specification.loader is not None,
            "cannot load K9/R5/B1024 runner")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


RUNNER = load_runner()
BASE = RUNNER.BASE


def synthetic_round(
    order: tuple[str, ...], *, main_execution: float = 1.3,
    main_one_shot: float = 1.25, control_execution: float = 1.4,
    control_one_shot: float = 1.35,
) -> dict[str, Any]:
    digests = {
        "original_data": "1" * 16,
        "transmitted_parity": "2" * 16,
        "recovered_originals": "1" * 16,
    }
    execution = {
        "candidate": 1.0, "control": control_execution,
        "main": main_execution}
    one_shot = {
        "candidate": 1.0, "control": control_one_shot,
        "main": main_one_shot}
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


def identity(path: Path) -> dict[str, Any]:
    return BASE.T8_SUPPORT.regular_file_identity(path)


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
                RUNNER.SOURCE_ROOT).as_posix(),
            "sha256": source_identity["sha256"],
            "size": source_identity["size"],
            "mode": 0o644,
        })
    return {
        "schema": RUNNER.BUILD_PROVENANCE.PRODUCTION_BUILD_CLOSURE_SCHEMA,
        "build_root": str(build),
        "physical_source_root": str(RUNNER.SOURCE_ROOT),
        "source_root": str(RUNNER.SOURCE_ROOT),
        "executable_target": "bench_leopard2",
        "tracked_source_manifest": {
            "files": dependencies,
            "git": {"commit": commit, "tree": tree, "dirty": False},
        },
        "validated_cache": {
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                RUNNER.BUILD_PROVENANCE.
                    BENCHMARK_BUILD_CONFIGURATION_SCHEMA,
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
                "source": identity(RUNNER.SOURCE_ROOT / "leopard2.cpp"),
                "object": archive_object_identity,
                "compile_entry": {"arguments": ["c++"]},
                "flag_profile": "portable-core",
            },
            {
                "role": "benchmark",
                "source": identity(
                    RUNNER.SOURCE_ROOT / "bench/leopard2/benchmark.cpp"),
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


def exercise_mocked_provenance_contract() -> None:
    with tempfile.TemporaryDirectory(
            prefix="leo2-k9r5-b1024-provenance-") as directory:
        build = Path(directory)
        executable = build / "bench_leopard2"
        executable.write_bytes(b"fixture executable\n")
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

        saved_capture = RUNNER.BUILD_PROVENANCE.candidate_build_provenance
        saved_replay = \
            RUNNER.BUILD_PROVENANCE.verify_reproducible_candidate_build
        saved_validate = \
            RUNNER.BUILD_PROVENANCE.validate_reproducible_build_proof
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
            RUNNER._BUILD_CLOSURE_BASELINE = None
            current["value"] = pristine if value is None else value
            calls.clear()

        try:
            RUNNER.BUILD_PROVENANCE.candidate_build_provenance = capture
            RUNNER.BUILD_PROVENANCE.verify_reproducible_candidate_build = \
                replay
            RUNNER.BUILD_PROVENANCE.validate_reproducible_build_proof = \
                validate

            reset()
            closure = BASE.build_closure_identity(options)
            require(
                closure["candidate_build"] == pristine and
                closure["reproducible_build"] ==
                    {"schema": "mock-replay", "valid": True} and
                len(closure["runner_source_dependencies"]) ==
                    len(BASE.RUNNER_DEPENDENCIES),
                "production provenance result changed")
            capture_call = calls[0]
            require(
                capture_call == (
                    "capture",
                    (build.resolve(), RUNNER.SOURCE_ROOT,
                     executable.resolve(), "bench_leopard2"),
                    {}) and
                any(call[0] == "replay" and call[2] == {"jobs": 1}
                    for call in calls) and
                any(call[0] == "validate" for call in calls),
                "production provenance was not called with the exact build "
                "target or replay contract")
            BASE.build_closure_identity(options)
            require(sum(call[0] == "capture" for call in calls) == 2 and
                    sum(call[0] == "replay" for call in calls) == 1 and
                    sum(call[0] == "validate" for call in calls) == 2,
                    "post-campaign provenance did not recapture the graph or "
                    "revalidate the retained replay proof")

            different = argparse.Namespace(
                candidate=executable,
                control=build / "another-bench_leopard2",
                candidate_sha256=executable_sha256,
                control_sha256=executable_sha256,
                source_commit=commit, source_tree=tree)
            different.control.write_bytes(b"other\n")
            different.control.chmod(0o755)
            reset()
            require_rejected(
                lambda: BASE.build_closure_identity(different),
                "distinct candidate/control executables were accepted")

            mismatched_hash = argparse.Namespace(
                candidate=executable, control=executable,
                candidate_sha256="0" * 64,
                control_sha256="0" * 64,
                source_commit=commit, source_tree=tree)
            reset()
            require_rejected(
                lambda: BASE.build_closure_identity(mismatched_hash),
                "production closure accepted a hash other than the frozen "
                "benchmark identity")
            mismatched_labels = argparse.Namespace(
                candidate=executable, control=executable,
                candidate_sha256=executable_sha256,
                control_sha256="0" * 64,
                source_commit=commit, source_tree=tree)
            reset()
            require_rejected(
                lambda: BASE.build_closure_identity(mismatched_labels),
                "candidate/control required hashes were allowed to differ")

            for mutate, message in (
                (lambda value: value.update(
                    {"build_root": str(build / "other")}),
                 "wrong build root was accepted"),
                (lambda value: value["validated_cache"].update({
                    "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                        "old-schema"}),
                 "wrong CMake configuration schema was accepted"),
                (lambda value: value.update(
                    {"source_object_compile_closure": []}),
                 "empty compile/object closure was accepted"),
                (lambda value: value["source_object_compile_closure"][0].
                    update({"flag_profile": {"avx2": True}}),
                 "malformed compile flag profile was accepted"),
                (lambda value: value.update({"archive_members": ["other.o"]}),
                 "inconsistent archive graph was accepted"),
                (lambda value: value["tracked_source_manifest"]["files"].pop(),
                 "missing runner dependency was accepted"),
                (lambda value: value.update({"executable_link_command": []}),
                 "empty executable link graph was accepted"),
            ):
                malformed = deepcopy(pristine)
                mutate(malformed)
                reset(malformed)
                require_rejected(
                    lambda: BASE.build_closure_identity(options), message)

            def bad_replay(*arguments: Any, **keywords: Any) -> dict[str, Any]:
                del arguments, keywords
                return {"schema": "unvalidated-replay"}

            RUNNER.BUILD_PROVENANCE.verify_reproducible_candidate_build = \
                bad_replay
            reset()
            require_rejected(
                lambda: BASE.build_closure_identity(options),
                "invalid reproducible-build proof was accepted")
        finally:
            RUNNER.BUILD_PROVENANCE.candidate_build_provenance = saved_capture
            RUNNER.BUILD_PROVENANCE.verify_reproducible_candidate_build = \
                saved_replay
            RUNNER.BUILD_PROVENANCE.validate_reproducible_build_proof = \
                saved_validate
            RUNNER._BUILD_CLOSURE_BASELINE = None


def main() -> int:
    cells = RUNNER.cells()
    require(
        [(cell["K"], cell["R"], cell["bytes"]) for cell in cells] == [
            (9, 5, 1024),
            (9, 5, 1023),
            (9, 5, 1025),
            (9, 6, 1024),
            (8, 5, 1024),
            (10, 5, 1024),
            (9, 5, 256),
        ],
        "K9/R5/B1024 dimensions or order changed")
    require(
        [cell["seed"] for cell in cells] ==
            list(range(0x09541000, 0x09541007)) and
        all(cell["batch"] == 1 and cell["reuse"] == 8192 and
            cell["loss"] == 1 and cell["measure_one_shot"] is True
            for cell in cells),
        "K9/R5/B1024 timing or seed policy changed")
    require(
        {cell["id"] for cell in cells if cell["compare_main"]} == {
            "target-k9-r5-b1024-q1",
            "retained-k9-r5-b256-q1",
        } and
        [cell["role"] for cell in cells].count("target") == 1,
        "exact-main or target classification changed")

    runner_path = Path(RUNNER.__file__).resolve()
    expected_dependencies = (
        runner_path,
        RUNNER.BASE_PATH,
        Path(BASE.T8_SUPPORT.__file__).resolve(),
        Path(BASE.MAIN_SUPPORT.__file__).resolve(),
        RUNNER.PROVENANCE_PATH.resolve(),
    )
    require(BASE.RUNNER_PATH == runner_path and
            BASE.RUNNER_DEPENDENCIES == expected_dependencies and
            len(set(expected_dependencies)) == 5,
            "runner dependency attestation changed")
    require(
        BASE.SCHEMA ==
            "leopard2-gf8-k9r5-b1024-terminal-abba/v1" and
        BASE.SUMMARY_SCHEMA ==
            "leopard2-gf8-k9r5-b1024-terminal-summary/v1" and
        BASE.MODE_SYMBOL ==
            "_ZN12_GLOBAL__N_1L26g_k9r5_b1024_terminal_modeE" and
        BASE.ALLOW_IDENTICAL_CANDIDATE_CONTROL is True and
        BASE.CONTROL_EXTRA_ARGUMENTS ==
            ("--disable-k9r5-b1024-terminal",) and
        BASE.CONTROL_BUILD_MARKER ==
            "k9r5_b1024_terminal_diagnostic_disabled" and
        BASE.REQUIRE_EXPECTED_IDENTITIES is True and
        BASE.REQUIRE_BUILD_CLOSURE is True and
        BASE.REQUIRE_FULL_ELF_IDENTITY is True and
        BASE.TARGET_CONTROL_FLOOR == 1.05 and
        BASE.TARGET_MAIN_FLOOR == 1.05 and
        BASE.NEIGHBOR_EQUIVALENCE_LOWER == 1.0 / 1.02 and
        BASE.NEIGHBOR_EQUIVALENCE_UPPER == 1.02 and
        BASE.RETAINED_MAIN_FLOOR == 1.0,
        "K9/R5/B1024 attribution or promotion policy changed")
    require(
        BASE.CANONICAL_MAIN_SHA256 ==
            "78575e2fe1b9796d10dcb4c14d2bcc1b627a79e787be66d34dc7aa666899aa93",
        "canonical exact-main executable identity changed")

    digest = "a" * 64
    saved_argv = sys.argv
    try:
        sys.argv = [
            str(runner_path),
            "--candidate", "/build/bench_leopard2",
            "--candidate-sha256", digest,
            "--control", "/build/bench_leopard2",
            "--control-sha256", digest,
            "--main", "/frozen/main",
            "--main-sha256", BASE.CANONICAL_MAIN_SHA256,
            "--source-commit", "e" * 40,
            "--source-tree", "f" * 40,
            "--output", "/results/k9r5-b1024",
            "--cpu", "30", "--sibling", "31",
        ]
        parsed = BASE.parse_arguments()
    finally:
        sys.argv = saved_argv
    require(
        parsed.candidate == Path("/build/bench_leopard2") and
        parsed.control == parsed.candidate and
        parsed.candidate_archive is None and
        parsed.control_archive is None and
        parsed.candidate_compile_commands is None and
        parsed.control_compile_commands is None,
        "production build-root command-line contract changed")
    saved_argv = sys.argv
    try:
        sys.argv = [
            str(runner_path),
            "--candidate", "/build/bench_leopard2",
            "--candidate-sha256", digest,
            "--control", "/build/bench_leopard2",
            "--control-sha256", digest,
            "--main", "/frozen/main",
            "--main-sha256", "b" * 64,
            "--source-commit", "e" * 40,
            "--source-tree", "f" * 40,
            "--output", "/results/k9r5-b1024",
            "--cpu", "30", "--sibling", "31",
        ]
        require_rejected(
            BASE.parse_arguments,
            "noncanonical exact-main executable SHA-256 was accepted")
    finally:
        sys.argv = saved_argv
    exercise_mocked_provenance_contract()

    orders = BASE.select_round_orders(BASE.TARGET_ORDER, 9)
    neighbor_orders = BASE.select_round_orders(BASE.NEIGHBOR_ORDER, 9)
    require(len(orders) == 9 and len(neighbor_orders) == 9 and
            orders[:3] == BASE.TARGET_ORDER and
            neighbor_orders[:4] == BASE.NEIGHBOR_ORDER,
            "nine-round balanced ordering changed")

    target = cells[0]
    candidate_command = BASE.benchmark_command(
        "candidate", Path("/frozen/candidate"), target, 30, 31, 16)
    control_command = BASE.benchmark_command(
        "control", Path("/frozen/control"), target, 30, 31, 16)
    require("--disable-k9r5-b1024-terminal" not in candidate_command and
            control_command.count("--disable-k9r5-b1024-terminal") == 1 and
            candidate_command[candidate_command.index("--reuse") + 1] ==
                "8192" and
            candidate_command[candidate_command.index("--iterations") + 1] ==
                "31" and
            candidate_command[candidate_command.index("--warmup") + 1] ==
                "16" and
            candidate_command[candidate_command.index("-c") + 1] == "30" and
            "--measure-one-shot-encode" in candidate_command,
            "candidate/control benchmark command changed")

    byte_control = cells[1]
    main_command = BASE.benchmark_command(
        "main", Path("/frozen/main"), byte_control, 30, 31, 16)
    candidate_byte_command = BASE.benchmark_command(
        "candidate", Path("/frozen/candidate"), byte_control, 30, 31, 16)
    require(main_command[main_command.index("--bytes") + 1] == "1024" and
            main_command[main_command.index("--logical-bytes") + 1] ==
                "1023" and
            candidate_byte_command[
                candidate_byte_command.index("--bytes") + 1] == "1023" and
            "--logical-bytes" not in candidate_byte_command,
            "exact-main logical-byte adapter changed")

    target_rounds = [synthetic_round(order) for order in orders]
    analysis = BASE.analyze(target, target_rounds)
    require(math.isclose(
                analysis["control_over_candidate"]["speedup"], 1.4) and
            math.isclose(
                analysis["main_over_candidate"]["speedup"], 1.3) and
            math.isclose(
                analysis["one_shot_control_over_candidate"]["speedup"],
                1.35) and
            math.isclose(
                analysis["one_shot_main_over_candidate"]["speedup"], 1.25),
            "target execution or one-shot contrasts changed")
    control_rounds = [
        synthetic_round(
            order, control_execution=1.0, control_one_shot=1.0)
        for order in neighbor_orders]
    control_analysis = BASE.analyze(byte_control, control_rounds)
    require("main_over_candidate" not in control_analysis and
            "one_shot_main_over_candidate" not in control_analysis and
            math.isclose(
                control_analysis["control_over_candidate"]["speedup"], 1.0) and
            control_analysis[
                "neighbor_selector_inertness_validation"]["accepted"] is True,
            "same-binary-only control classification changed")

    leaked_speedup = [
        synthetic_round(
            order, control_execution=1.10, control_one_shot=1.10)
        for order in neighbor_orders]
    require_rejected(
        lambda: BASE.analyze(byte_control, leaked_speedup),
        "credible selector speedup leaked into a neighbor")
    leaked_regression = [
        synthetic_round(
            order, control_execution=0.90, control_one_shot=0.90)
        for order in neighbor_orders]
    require_rejected(
        lambda: BASE.analyze(byte_control, leaked_regression),
        "credible selector regression leaked into a neighbor")
    wide_uncertain_effect = [
        synthetic_round(
            order,
            control_execution=0.5 if index % 2 == 0 else 2.0,
            control_one_shot=2.0 if index % 2 == 0 else 0.5)
        for index, order in enumerate(neighbor_orders)]
    require_rejected(
        lambda: BASE.analyze(byte_control, wide_uncertain_effect),
        "wide selector-effect confidence interval was treated as equivalent")

    retained = cells[-1]
    retained_rounds = [
        synthetic_round(
            order, control_execution=1.0, control_one_shot=1.0)
        for order in orders]
    retained_analysis = BASE.analyze(retained, retained_rounds)
    require(
        retained_analysis["retained_main_validation"]["accepted"] is True and
        retained_analysis["retained_main_validation"]["floor"] == 1.0 and
        retained_analysis["retained_main_validation"]["lower_ci95"] == {
            "encode_execution": 1.3,
            "one_shot_encode": 1.25,
        },
        "retained B256 exact-main acceptance gate changed")
    retained_regression = [
        synthetic_round(
            order, main_execution=0.90, main_one_shot=0.92,
            control_execution=1.0, control_one_shot=1.0)
        for order in orders
    ]
    require_rejected(
        lambda: BASE.analyze(retained, retained_regression),
        "credible retained B256 exact-main regression was accepted")

    print("K9/R5/B1024 hardened ABBA wrapper checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
