#!/usr/bin/env python3
"""Portable adversarial tests for the exact-main authority record."""

from __future__ import annotations

import ast
import builtins
import copy
import hashlib
import importlib.util
import io
import os
import pathlib
import subprocess
import sys
import unittest
from unittest import mock


def load_module(name: str, path: pathlib.Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load test dependency {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


TOOLS = pathlib.Path(__file__).resolve().parent
identity_contract = load_module(
    "leopard2_exact_main_baseline",
    TOOLS / "leopard2_exact_main_baseline.py",
)
contract = load_module(
    "leopard2_exact_main_baseline_record_under_test",
    TOOLS / "leopard2_exact_main_baseline_record.py",
)


def digest(label: str) -> str:
    return hashlib.sha256(label.encode("ascii")).hexdigest()


def oid(label: str) -> str:
    return hashlib.sha1(label.encode("ascii")).hexdigest()


def text_identity(path: str, label: str, *, empty: bool = False) -> dict:
    size = 0 if empty else len(label)
    return {
        "relative_path": path,
        "size": size,
        "sha256": (
            identity_contract.EMPTY_CONTENT_SHA256 if empty else digest(label)),
    }


def file_identity(path: str, label: str) -> dict:
    return {"relative_path": path, "size": len(label) + 1,
            "sha256": digest(label)}


def lane_fixture(*, successful: bool = True, failed_index: int = 2) -> dict:
    names = list(contract.STAGE_NAMES)
    if not successful:
        names = names[:failed_index + 1]
    return {
        "root": "/tmp/leopard-exact-main-baseline-a1",
        "attempt": 1,
        "attempt_budget": 3,
        "record_relative_path": (
            "baseline-authority.json" if successful else "FAILED.json"),
        "seal_protocol": "owner-only-tree-sha256sums/v1",
        "stages": [{
            "name": name,
            "status": (
                "complete" if successful or index < len(names) - 1
                else "failed"),
            "log": file_identity(
                f"logs/{index:02d}-{name}.log", f"stage-{name}"),
        } for index, name in enumerate(names)],
    }


def host_fixture() -> dict:
    return {
        "hostname": "fixture-host",
        "kernel": "Linux 6.8.0 fixture",
        "architecture": "x86_64",
        "cpu_model": "Fixture AMD Family 1Ah Model 08h",
        "online_cpus": [0, 1, 2, 3],
        "sc_clk_tck": 100,
    }


def archive_fixture(path: str, prefix: str, label: str) -> dict:
    sha256 = digest(label)
    return {
        "relative_path": path,
        "prefix": prefix,
        "size": 8192 + len(label),
        "sha256": sha256,
        "replay_sha256": sha256,
        "replay_identical": True,
    }


def source_fixture() -> dict:
    return {
        "baseline": {
            "commit": contract.BASELINE_COMMIT,
            "tree": contract.BASELINE_TREE,
            "submodule": {
                "path": "sse2neon",
                "commit": contract.BASELINE_SSE2NEON_COMMIT,
            },
            "git_capture": file_identity(
                "source/leopard-main-git-capture.json", "baseline-git"),
            "archive": archive_fixture(
                "source/leopard-main-source.tar", "leopard-main-source/",
                "baseline-source-archive"),
        },
        "adapter_repository": {
            "commit": "1" * 40,
            "tree": "2" * 40,
            "clean": True,
            "git_capture": file_identity(
                "source/adapter-git-capture.json", "adapter-git"),
            "archive": archive_fixture(
                "source/leopard2-adapter-source.tar",
                "leopard2-adapter-source/", "adapter-source-archive"),
        },
    }


def adapter_fixture() -> dict:
    files = [{
        "path": path,
        "git_blob_sha1": oid(path),
        "size": 1000 + index,
        "sha256": digest(path),
    } for index, path in enumerate(contract.ADAPTER_PATHS)]
    return {
        "minimum_harness_commit": contract.MINIMUM_HARNESS_COMMIT,
        "files": files,
        "files_sha256": identity_contract.canonical_sha256(files),
    }


TOOL_PATHS = {
    "archiver": "/usr/bin/ar",
    "cmake": "/usr/bin/cmake",
    "compiler": "/usr/bin/x86_64-linux-gnu-g++-13",
    "ctest": "/usr/bin/ctest",
    "git": "/usr/bin/git",
    "ldd": "/usr/bin/ldd",
    "linker": "/usr/bin/x86_64-linux-gnu-ld",
    "make": "/usr/bin/make",
    "python": "/usr/bin/python3.12",
    "ranlib": "/usr/bin/ranlib",
    "assembler": "/usr/bin/x86_64-linux-gnu-as",
    "cc1plus": "/usr/lib/gcc/x86_64-linux-gnu/13/cc1plus",
    "collect2": "/usr/lib/gcc/x86_64-linux-gnu/13/collect2",
}


def tool_fixture(role: str) -> dict:
    return {
        "role": role,
        "path": TOOL_PATHS[role],
        "resolved_path": TOOL_PATHS[role],
        "size": 100_000 + len(role),
        "mode": 0o755,
        "sha256": digest("tool-" + role),
    }


def toolchain_fixture() -> dict:
    tools = [tool_fixture(role) for role in contract.TOOL_ROLES]
    subtools = [tool_fixture(role) for role in contract.SUBTOOL_ROLES]
    return {
        "tools": tools,
        "subtools": subtools,
        "versions": [{
            "role": role,
            "argv": [TOOL_PATHS[role], "--version"],
            "stdout": text_identity(
                f"toolchain/versions/{role}.stdout", f"{role}-version"),
            "stderr": text_identity(
                f"toolchain/versions/{role}.stderr",
                f"{role}-version-stderr", empty=True),
            "exit_status": 0,
        } for role in contract.VERSION_ROLES],
    }


def roots_fixture(role: str) -> dict:
    family = "canonical" if role != "path_variant" else "path-variant"
    base = f"/tmp/leopard-exact-main-{family}"
    return {
        "adapter_source_root": base + "/adapter",
        "baseline_source_root": base + "/baseline",
        "build_root": base + "/build",
    }


def build_fixture(role: str, *, executable_sha256: str,
                  archive_sha256: str) -> dict:
    roots = roots_fixture(role)
    role_path = role.replace("_", "-")
    return {
        "role": role,
        "roots": roots,
        "configure_argv": [
            TOOL_PATHS["cmake"],
            "-S", roots["adapter_source_root"] +
            "/experiments/leopard2/main_compare",
            "-B", roots["build_root"],
            "-G", "Unix Makefiles",
            "-DCMAKE_BUILD_TYPE=Release",
            "-DLEO_MAIN_PURE_AVX2=ON",
            "-DLEOPARD_MAIN_SOURCE_DIR=" + roots["baseline_source_root"],
            "-DCMAKE_CXX_COMPILER=" + TOOL_PATHS["compiler"],
        ],
        "build_argv": [
            TOOL_PATHS["cmake"], "--build", roots["build_root"],
            "--target", "leopard_main_benchmark", "--", "-j1",
        ],
        "configure_log": file_identity(
            f"builds/{role_path}/configure.log", role + "-configure"),
        "build_log": file_identity(
            f"builds/{role_path}/build.log", role + "-build"),
        "cmake_cache": file_identity(
            f"builds/{role_path}/CMakeCache.txt", role + "-cache"),
        "compile_commands": file_identity(
            f"builds/{role_path}/compile_commands.json",
            role + "-compile-commands"),
        "executable": {
            "name": "leopard_main_benchmark",
            "build_relative_path": "leopard_main_benchmark",
            "retained_relative_path": (
                f"artifacts/{role_path}/leopard_main_benchmark"),
            "size": 1_175_200,
            "sha256": executable_sha256,
        },
        "archive": {
            "name": "libleopard_main_exact.a",
            "build_relative_path": "libleopard_main_exact.a",
            "retained_relative_path": (
                f"artifacts/{role_path}/libleopard_main_exact.a"),
            "size": 2_345_678,
            "sha256": archive_sha256,
        },
        "closure": {
            "relative_path": f"builds/{role_path}/build-closure.json",
            "size": 12_345,
            "sha256": digest(role + "-closure"),
            "file_count": 23,
        },
    }


def builds_fixture() -> dict:
    canonical_executable = digest("canonical-executable")
    canonical_archive = digest("canonical-archive")
    return {
        "canonical_first": build_fixture(
            "canonical_first", executable_sha256=canonical_executable,
            archive_sha256=canonical_archive),
        "canonical_second": build_fixture(
            "canonical_second", executable_sha256=canonical_executable,
            archive_sha256=canonical_archive),
        "path_variant": build_fixture(
            "path_variant", executable_sha256=digest("variant-executable"),
            archive_sha256=digest("variant-archive")),
    }


def normalized_sections() -> list[dict]:
    return [
        {
            "index": 2, "name": ".text", "type": "PROGBITS",
            "flags": 0x6, "address": 0x401000, "size": 4096,
            "alignment": 64, "content_sha256": digest("text"),
        },
        {
            "index": 3, "name": ".rodata", "type": "PROGBITS",
            "flags": 0x2, "address": 0x403000, "size": 2048,
            "alignment": 32, "content_sha256": digest("rodata"),
        },
        {
            "index": 4, "name": ".bss", "type": "NOBITS",
            "flags": 0x3, "address": 0x405000, "size": 1024,
            "alignment": 32, "content_sha256": None,
        },
    ]


def zero_census(roots: dict) -> dict:
    sections = normalized_sections()
    return {
        "match_rule": "non-overlapping-exact-utf8-byte-substring/v1",
        "roots": [{
            "role": role,
            "path": roots[role],
            "occurrences": 0,
            "sections": [{"name": section["name"], "occurrences": 0}
                         for section in sections],
        } for role in (
            "adapter_source_root", "baseline_source_root", "build_root")],
    }


def identity_fixture(builds: dict) -> dict:
    identities = {}
    for role in contract.BUILD_ROLES:
        build = builds[role]
        identities[role] = identity_contract.normalized_code_identity_record(
            artifact={
                "size": build["executable"]["size"],
                "sha256": build["executable"]["sha256"],
            },
            sections=normalized_sections(),
            path_string_census=zero_census(build["roots"]),
        )
    return {
        **identities,
        "combined_sha256": identities["canonical_first"]["combined_sha256"],
        "canonical_raw_bytes_identical": True,
        "path_variant_raw_bytes_differ": True,
        "normalized_match": True,
    }


def runtime_fixture(builds: dict) -> dict:
    return {
        "schema": contract.RUNTIME_CLOSURE_SCHEMA,
        "normalization": "canonical-ldd-C-v1",
        "records": [{
            "role": role,
            "executable_sha256": builds[role]["executable"]["sha256"],
            "canonical_ldd_output": text_identity(
                f"runtime/{role.replace('_', '-')}/ldd.txt",
                "canonical-ldd"),
            "dependencies": [
                {
                    "soname": "ld-linux-x86-64.so.2", "kind": "file",
                    "path": "/usr/lib/x86_64-linux-gnu/ld-linux-x86-64.so.2",
                    "size": 234_000, "sha256": digest("ld-linux"),
                },
                {
                    "soname": "libc.so.6", "kind": "file",
                    "path": "/usr/lib/x86_64-linux-gnu/libc.so.6",
                    "size": 2_000_000, "sha256": digest("libc"),
                },
                {
                    "soname": "linux-vdso.so.1", "kind": "virtual",
                    "path": None, "size": None, "sha256": None,
                },
            ],
        } for role in contract.BUILD_ROLES],
    }


def attestation_fixture(builds: dict, adapter: dict) -> dict:
    records = []
    for role in contract.BUILD_ROLES:
        build = builds[role]
        executable = (build["roots"]["build_root"] + "/" +
                      build["executable"]["build_relative_path"])
        records.append({
            "role": role,
            "argv": [
                executable,
                "--k", "8", "--r", "4", "--bytes", "64", "--loss", "1",
                "--batch", "2", "--reuse", "1", "--iterations", "2",
                "--warmup", "1", "--threads", "1", "--seed", "7",
                "--json", "-",
            ],
            "stdout": text_identity(
                f"attestations/{role.replace('_', '-')}/benchmark.stdout.json",
                role + "-probe-stdout"),
            "stderr": text_identity(
                f"attestations/{role.replace('_', '-')}/benchmark.stderr",
                role + "-probe-stderr", empty=True),
            "exit_status": 0,
            "reported_schema": "leopard-main-benchmark-v1",
            "main_source_commit": contract.BASELINE_COMMIT,
            "pure_avx2": True,
            "round_trip": True,
            "ctest": {
                "argv": [
                    TOOL_PATHS["ctest"], "--test-dir",
                    build["roots"]["build_root"], "--output-on-failure",
                    "-R", "^leopard_main_benchmark_smoke$",
                ],
                "stdout": text_identity(
                    f"attestations/{role.replace('_', '-')}/ctest.stdout.log",
                    role + "-ctest-stdout"),
                "stderr": text_identity(
                    f"attestations/{role.replace('_', '-')}/ctest.stderr.log",
                    role + "-ctest-stderr", empty=True),
                "exit_status": 0,
                "passed": 1,
                "failed": 0,
            },
        })
    return {
        "schema": contract.ATTESTATION_SCHEMA,
        "test_controller": {
            "relative_path": "controllers/test_legacy_main_benchmark.py",
            "size": adapter["files"][2]["size"],
            "sha256": adapter["files"][2]["sha256"],
        },
        "records": records,
    }


def authority_fixture() -> dict:
    builds = builds_fixture()
    adapter = adapter_fixture()
    return contract.baseline_authority_record(
        created_utc="2026-08-29T20:30:00Z",
        host=host_fixture(),
        lane=lane_fixture(),
        source=source_fixture(),
        adapter=adapter,
        toolchain=toolchain_fixture(),
        builds=builds,
        runtime_closure=runtime_fixture(builds),
        attestation=attestation_fixture(builds, adapter),
        identity=identity_fixture(builds),
    )


def resign(value: dict) -> dict:
    changed = copy.deepcopy(value)
    changed["record_sha256"] = identity_contract.canonical_sha256({
        key: copy.deepcopy(field) for key, field in changed.items()
        if key != "record_sha256"
    })
    return changed


def walk_dicts(value: object, path: tuple[object, ...] = ()):
    if type(value) is dict:
        yield path, value
        for key, child in value.items():
            yield from walk_dicts(child, path + (key,))
    elif type(value) is list:
        for index, child in enumerate(value):
            yield from walk_dicts(child, path + (index,))


def at_path(value: object, path: tuple[object, ...]):
    current = value
    for component in path:
        current = current[component]  # type: ignore[index]
    return current


class ExactMainBaselineAuthorityRecordTest(unittest.TestCase):
    def assertAuthorityRejected(self, value: object) -> None:  # noqa: N802
        with self.assertRaises(contract.ExactMainBaselineRecordError):
            contract.validate_baseline_authority_record(value)

    def assertFailureRejected(self, value: object) -> None:  # noqa: N802
        with self.assertRaises(contract.ExactMainBaselineRecordError):
            contract.validate_baseline_failure_record(value)

    def test_literal_golden_roundtrip_and_detachment(self) -> None:
        record = authority_fixture()
        self.assertEqual(record["schema"],
                         "leopard2-gf8-exact-main-pure-avx2-baseline/v1")
        self.assertEqual(record["source"]["baseline"]["commit"],
                         "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198")
        self.assertEqual(
            record["record_sha256"],
            "949dd158451daa32e1d99ad531ceda4b9a33d13ec315a241b24822471088c998",
        )
        encoded = identity_contract.canonical_json_bytes(record)
        self.assertEqual(contract.load_baseline_authority_record(encoded),
                         record)
        loaded = contract.validate_baseline_authority_record(record)
        loaded["host"]["hostname"] = "changed"
        self.assertEqual(record["host"]["hostname"], "fixture-host")

        host = host_fixture()
        builds = builds_fixture()
        adapter = adapter_fixture()
        created = contract.baseline_authority_record(
            created_utc="2026-08-29T20:30:00Z",
            host=host,
            lane=lane_fixture(),
            source=source_fixture(),
            adapter=adapter,
            toolchain=toolchain_fixture(),
            builds=builds,
            runtime_closure=runtime_fixture(builds),
            attestation=attestation_fixture(builds, adapter),
            identity=identity_fixture(builds),
        )
        host["hostname"] = "mutated"
        builds["canonical_first"]["executable"]["sha256"] = "f" * 64
        self.assertEqual(created, record)

    def test_every_nested_object_rejects_missing_and_extra_keys(self) -> None:
        record = authority_fixture()
        paths = [path for path, item in walk_dicts(record) if item]
        self.assertEqual(len(paths), 220)
        for path in paths:
            original = at_path(record, path)
            key = next(iter(original))
            with self.subTest(path=path, mutation="missing"):
                changed = copy.deepcopy(record)
                del at_path(changed, path)[key]
                self.assertAuthorityRejected(changed)
            with self.subTest(path=path, mutation="extra"):
                changed = copy.deepcopy(record)
                at_path(changed, path)["unexpected_key"] = None
                self.assertAuthorityRejected(changed)

    def test_strict_json_and_exact_scalar_types(self) -> None:
        record = authority_fixture()
        encoded = identity_contract.canonical_json_bytes(record)
        hostile = (
            encoded + b"{}\n",
            b'{"schema":"one","schema":"two"}\n',
            b'{"value":NaN}\n',
            b'{"value":Infinity}\n',
            b'{"value":1e309}\n',
            b"[" * 10_000 + b"]" * 10_000,
        )
        for index, payload in enumerate(hostile):
            with self.subTest(index=index):
                with self.assertRaises(contract.ExactMainBaselineRecordError):
                    contract.load_baseline_authority_record(payload)

        mutations = []
        changed = copy.deepcopy(record)
        changed["host"]["sc_clk_tck"] = True
        mutations.append(changed)
        changed = copy.deepcopy(record)
        changed["lane"]["attempt"] = 1.0
        mutations.append(changed)
        changed = copy.deepcopy(record)
        changed["toolchain"]["tools"][0]["mode"] = True
        mutations.append(changed)
        changed = copy.deepcopy(record)
        changed["attestation"]["records"][0]["ctest"]["passed"] = True
        mutations.append(changed)
        for index, changed in enumerate(mutations):
            with self.subTest(type_index=index):
                self.assertAuthorityRejected(resign(changed))

    def test_public_producers_normalize_errors_and_submodule_types(self) -> None:
        builds = builds_fixture()
        adapter = adapter_fixture()
        malformed_host = host_fixture()
        malformed_host["cpu_model"] = object()
        with self.assertRaises(
                contract.ExactMainBaselineRecordError) as authority_error:
            contract.baseline_authority_record(
                created_utc="2026-08-29T20:30:00Z",
                host=malformed_host,
                lane=lane_fixture(),
                source=source_fixture(),
                adapter=adapter,
                toolchain=toolchain_fixture(),
                builds=builds,
                runtime_closure=runtime_fixture(builds),
                attestation=attestation_fixture(builds, adapter),
                identity=identity_fixture(builds),
            )
        self.assertIs(type(authority_error.exception),
                      contract.ExactMainBaselineRecordError)

        failure_lane = lane_fixture(successful=False)
        with self.assertRaises(
                contract.ExactMainBaselineRecordError) as failure_error:
            contract.baseline_acquisition_failure_record(
                created_utc="2026-08-29T20:31:00Z",
                lane=failure_lane,
                stage=failure_lane["stages"][-1]["name"],
                error={"kind": "build_error", "message": object(),
                       "exit_status": 1},
                retained_files=[],
            )
        self.assertIs(type(failure_error.exception),
                      contract.ExactMainBaselineRecordError)

        class EqualToAnything:
            def __eq__(self, other: object) -> bool:
                return True

        record = authority_fixture()
        for key in ("path", "commit"):
            with self.subTest(submodule_key=key):
                changed = copy.deepcopy(record)
                changed["source"]["baseline"]["submodule"][key] = \
                    EqualToAnything()
                with self.assertRaises(
                        contract.ExactMainBaselineRecordError) as error:
                    contract.validate_baseline_authority_record(changed)
                self.assertIs(type(error.exception),
                              contract.ExactMainBaselineRecordError)

    def test_source_adapter_and_toolchain_drift_fail_closed(self) -> None:
        record = authority_fixture()
        changes = []
        changed = copy.deepcopy(record)
        changed["source"]["baseline"]["commit"] = "0" * 40
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["source"]["baseline"]["tree"] = "0" * 40
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["source"]["baseline"]["submodule"]["commit"] = "0" * 40
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["source"]["baseline"]["archive"]["replay_identical"] = False
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["adapter"]["minimum_harness_commit"] = "0" * 40
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["adapter"]["files"][0]["sha256"] = "f" * 64
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["toolchain"]["tools"][0], changed["toolchain"]["tools"][1] = (
            changed["toolchain"]["tools"][1],
            changed["toolchain"]["tools"][0])
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["toolchain"]["versions"][0]["argv"][0] = "/usr/bin/false"
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["configure"]["isa_flags"].append("-march=native")
        changes.append(changed)
        for index, changed in enumerate(changes):
            with self.subTest(index=index):
                self.assertAuthorityRejected(resign(changed))

    def test_build_identity_census_and_promotion_joins(self) -> None:
        record = authority_fixture()
        changes = []
        changed = copy.deepcopy(record)
        changed["builds"]["canonical_second"]["executable"]["sha256"] = \
            digest("different-second")
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["builds"]["canonical_second"]["roots"]["build_root"] += "-x"
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["builds"]["path_variant"]["configure_argv"][0] = "/usr/bin/false"
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["builds"]["canonical_second"]["archive"]["sha256"] = \
            digest("different-second-archive")
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["identity"]["path_variant"]["artifact"]["sha256"] = "f" * 64
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["identity"]["normalized_match"] = False
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["promotion"]["selected_section_census_zero"] = False
        changed["promotion"]["promoted"] = False
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["superseded_references"][0]["authority"] = True
        changes.append(changed)
        for index, changed in enumerate(changes):
            with self.subTest(index=index):
                self.assertAuthorityRejected(resign(changed))

        nested_resigned = copy.deepcopy(record)
        old_identity = nested_resigned["identity"]["path_variant"]
        nested_resigned["identity"]["path_variant"] = \
            identity_contract.normalized_code_identity_record(
                artifact={
                    "size": old_identity["artifact"]["size"],
                    "sha256": "f" * 64,
                },
                sections=copy.deepcopy(old_identity["sections"]),
                path_string_census=copy.deepcopy(
                    old_identity["path_string_census"]),
            )
        with self.assertRaisesRegex(
                contract.ExactMainBaselineRecordError,
                "not bound to its executable"):
            contract.validate_baseline_authority_record(
                resign(nested_resigned))

        archive_failure = copy.deepcopy(record)
        archive_failure["builds"]["canonical_second"]["archive"][
            "sha256"] = digest("honest-different-second-archive")
        archive_failure["promotion"][
            "same_path_archive_bytes_identical"] = False
        archive_failure["promotion"]["promoted"] = False
        with self.assertRaisesRegex(
                contract.ExactMainBaselineRecordError,
                "did not satisfy every promotion gate"):
            contract.validate_baseline_authority_record(
                resign(archive_failure))

        canonical_identity_failure = copy.deepcopy(record)
        second_build = canonical_identity_failure["builds"][
            "canonical_second"]
        second_sha256 = digest("honest-different-second-executable")
        second_build["executable"]["sha256"] = second_sha256
        canonical_identity_failure["runtime_closure"]["records"][1][
            "executable_sha256"] = second_sha256
        old_second = canonical_identity_failure["identity"][
            "canonical_second"]
        canonical_identity_failure["identity"]["canonical_second"] = \
            identity_contract.normalized_code_identity_record(
                artifact={
                    "size": second_build["executable"]["size"],
                    "sha256": second_sha256,
                },
                sections=copy.deepcopy(old_second["sections"]),
                path_string_census=copy.deepcopy(
                    old_second["path_string_census"]),
            )
        canonical_identity_failure["identity"][
            "canonical_raw_bytes_identical"] = False
        with self.assertRaisesRegex(
                contract.ExactMainBaselineRecordError,
                "same-path normalized identities were not byte-identical"):
            contract.validate_baseline_authority_record(
                resign(canonical_identity_failure))

        variant_raw_failure = copy.deepcopy(record)
        variant_build = variant_raw_failure["builds"]["path_variant"]
        canonical_sha256 = variant_raw_failure["builds"]["canonical_first"][
            "executable"]["sha256"]
        variant_build["executable"]["sha256"] = canonical_sha256
        variant_raw_failure["runtime_closure"]["records"][2][
            "executable_sha256"] = canonical_sha256
        old_variant = variant_raw_failure["identity"]["path_variant"]
        variant_raw_failure["identity"]["path_variant"] = \
            identity_contract.normalized_code_identity_record(
                artifact={
                    "size": variant_build["executable"]["size"],
                    "sha256": canonical_sha256,
                },
                sections=copy.deepcopy(old_variant["sections"]),
                path_string_census=copy.deepcopy(
                    old_variant["path_string_census"]),
            )
        variant_raw_failure["identity"][
            "path_variant_raw_bytes_differ"] = False
        with self.assertRaisesRegex(
                contract.ExactMainBaselineRecordError,
                "path-variant raw executable did not differ"):
            contract.validate_baseline_authority_record(
                resign(variant_raw_failure))

        normalized_failure = copy.deepcopy(record)
        old_variant = normalized_failure["identity"]["path_variant"]
        changed_sections = copy.deepcopy(old_variant["sections"])
        changed_sections[1]["content_sha256"] = digest(
            "honest-different-path-variant-rodata")
        normalized_failure["identity"]["path_variant"] = \
            identity_contract.normalized_code_identity_record(
                artifact=copy.deepcopy(old_variant["artifact"]),
                sections=changed_sections,
                path_string_census=copy.deepcopy(
                    old_variant["path_string_census"]),
            )
        normalized_failure["identity"]["normalized_match"] = False
        with self.assertRaisesRegex(
                contract.ExactMainBaselineRecordError,
                "path-variant normalized identities did not match"):
            contract.validate_baseline_authority_record(
                resign(normalized_failure))

        contaminated = copy.deepcopy(record)
        variant = contaminated["identity"]["path_variant"]
        census = copy.deepcopy(variant["path_string_census"])
        census["roots"][2]["occurrences"] = 1
        census["roots"][2]["sections"][1]["occurrences"] = 1
        rebuilt = identity_contract.normalized_code_identity_record(
            artifact=copy.deepcopy(variant["artifact"]),
            sections=copy.deepcopy(variant["sections"]),
            path_string_census=census,
        )
        contaminated["identity"]["path_variant"] = rebuilt
        contaminated["promotion"]["selected_section_census_zero"] = False
        contaminated["promotion"]["promoted"] = False
        self.assertEqual(rebuilt["combined_sha256"],
                         contaminated["identity"]["combined_sha256"])
        with self.assertRaisesRegex(
                contract.ExactMainBaselineRecordError,
                "did not satisfy every promotion gate"):
            contract.validate_baseline_authority_record(resign(contaminated))

    def test_global_paths_and_independent_roots_fail_closed(self) -> None:
        record = authority_fixture()
        changes = []
        changed = copy.deepcopy(record)
        changed["lane"]["stages"][0]["log"]["relative_path"] = \
            "logs/00-wrong-source-acquisition.log"
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["source"]["baseline"]["git_capture"]["relative_path"] = \
            "source/wrong-leopard-main-git-capture.json"
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["source"]["adapter_repository"]["git_capture"][
            "relative_path"] = "source/wrong-adapter-git-capture.json"
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["builds"]["canonical_second"]["configure_log"][
            "relative_path"] = changed["builds"]["canonical_first"][
                "configure_log"]["relative_path"]
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["builds"]["canonical_second"]["cmake_cache"] = \
            copy.deepcopy(changed["builds"]["canonical_first"]["cmake_cache"])
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["lane"]["root"] = changed["builds"]["canonical_first"][
            "roots"]["build_root"]
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["lane"]["stages"][0]["log"]["relative_path"] = \
            changed["source"]["baseline"]["archive"]["relative_path"]
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["lane"]["stages"][0]["log"]["relative_path"] = "builds"
        changes.append(changed)

        changed = copy.deepcopy(record)
        variant = changed["builds"]["path_variant"]
        canonical_roots = changed["builds"]["canonical_first"]["roots"]
        variant["roots"] = {
            key: canonical_roots[key] + "/path-variant"
            for key in canonical_roots
        }
        roots = variant["roots"]
        variant["configure_argv"] = [
            TOOL_PATHS["cmake"],
            "-S", roots["adapter_source_root"] +
            "/experiments/leopard2/main_compare",
            "-B", roots["build_root"],
            "-G", "Unix Makefiles",
            "-DCMAKE_BUILD_TYPE=Release",
            "-DLEO_MAIN_PURE_AVX2=ON",
            "-DLEOPARD_MAIN_SOURCE_DIR=" + roots["baseline_source_root"],
            "-DCMAKE_CXX_COMPILER=" + TOOL_PATHS["compiler"],
        ]
        variant["build_argv"] = [
            TOOL_PATHS["cmake"], "--build", roots["build_root"],
            "--target", "leopard_main_benchmark", "--", "-j1",
        ]
        changes.append(changed)

        for index, changed in enumerate(changes):
            with self.subTest(index=index):
                self.assertAuthorityRejected(resign(changed))

    def test_byte_identity_and_host_path_functions_fail_closed(self) -> None:
        record = authority_fixture()
        changes = []
        derived_claim_changes = []

        changed = copy.deepcopy(record)
        changed["builds"]["path_variant"]["build_log"]["sha256"] = \
            changed["builds"]["canonical_first"]["configure_log"]["sha256"]
        changes.append(changed)

        changed = copy.deepcopy(record)
        changed["lane"]["stages"][5]["log"]["sha256"] = \
            changed["lane"]["stages"][0]["log"]["sha256"]
        changes.append(changed)

        changed = copy.deepcopy(record)
        baseline_archive = changed["source"]["baseline"]["archive"]
        adapter_archive = changed["source"]["adapter_repository"]["archive"]
        adapter_archive["sha256"] = baseline_archive["sha256"]
        adapter_archive["replay_sha256"] = baseline_archive["sha256"]
        changes.append(changed)

        changed = copy.deepcopy(record)
        baseline_archive = changed["source"]["baseline"]["archive"]
        adapter_archive = changed["source"]["adapter_repository"]["archive"]
        adapter_archive["size"] = baseline_archive["size"]
        adapter_archive["sha256"] = baseline_archive["sha256"]
        adapter_archive["replay_sha256"] = baseline_archive["sha256"]
        derived_claim_changes.append((
            changed, "one SHA-256 with two content-derived identities"))

        changed = copy.deepcopy(record)
        first_closure = changed["builds"]["canonical_first"]["closure"]
        second_closure = changed["builds"]["canonical_second"]["closure"]
        second_closure["size"] = first_closure["size"]
        second_closure["sha256"] = first_closure["sha256"]
        second_closure["file_count"] = first_closure["file_count"] + 7
        derived_claim_changes.append((
            changed, "one SHA-256 with two content-derived identities"))

        changed = copy.deepcopy(record)
        adapter_files = changed["adapter"]["files"]
        adapter_files[0]["git_blob_sha1"] = adapter_files[1]["git_blob_sha1"]
        changed["adapter"]["files_sha256"] = \
            identity_contract.canonical_sha256(adapter_files)
        derived_claim_changes.append((
            changed, "reuse one Git blob object ID at two identities"))

        changed = copy.deepcopy(record)
        adapter_files = changed["adapter"]["files"]
        adapter_files[1]["size"] = adapter_files[0]["size"]
        adapter_files[1]["sha256"] = adapter_files[0]["sha256"]
        changed["adapter"]["files_sha256"] = \
            identity_contract.canonical_sha256(adapter_files)
        derived_claim_changes.append((
            changed, "reuse one byte identity with two Git blob object IDs"))

        changed = copy.deepcopy(record)
        changed["source"]["adapter_repository"]["commit"] = \
            changed["source"]["baseline"]["commit"]
        derived_claim_changes.append((
            changed, "reuse one commit object ID with two trees"))

        changed = copy.deepcopy(record)
        changed["source"]["adapter_repository"]["tree"] = \
            changed["source"]["baseline"]["commit"]
        derived_claim_changes.append((
            changed, "reuse one Git object ID as two object types"))

        changed = copy.deepcopy(record)
        compiler = next(item for item in changed["toolchain"]["tools"]
                        if item["role"] == "compiler")
        cc1plus = next(item for item in changed["toolchain"]["subtools"]
                       if item["role"] == "cc1plus")
        compiler["sha256"] = cc1plus["sha256"]
        changes.append(changed)

        changed = copy.deepcopy(record)
        archiver = next(item for item in changed["toolchain"]["tools"]
                        if item["role"] == "archiver")
        ranlib = next(item for item in changed["toolchain"]["tools"]
                      if item["role"] == "ranlib")
        ranlib["resolved_path"] = archiver["resolved_path"]
        version = next(item for item in changed["toolchain"]["versions"]
                       if item["role"] == "ranlib")
        version["argv"][0] = archiver["resolved_path"]
        changes.append(changed)

        changed = copy.deepcopy(record)
        compiler = next(item for item in changed["toolchain"]["tools"]
                        if item["role"] == "compiler")
        cc1plus = next(item for item in changed["toolchain"]["subtools"]
                       if item["role"] == "cc1plus")
        cc1plus["resolved_path"] = compiler["resolved_path"]
        version = next(item for item in changed["toolchain"]["versions"]
                       if item["role"] == "cc1plus")
        version["argv"][0] = compiler["resolved_path"]
        changes.append(changed)

        changed = copy.deepcopy(record)
        archiver = next(item for item in changed["toolchain"]["tools"]
                        if item["role"] == "archiver")
        ranlib = next(item for item in changed["toolchain"]["tools"]
                      if item["role"] == "ranlib")
        ranlib["path"] = archiver["path"]
        derived_claim_changes.append((
            changed, "resolve one tool path two ways"))

        changed = copy.deepcopy(record)
        archiver = next(item for item in changed["toolchain"]["tools"]
                        if item["role"] == "archiver")
        ranlib = next(item for item in changed["toolchain"]["tools"]
                      if item["role"] == "ranlib")
        ranlib["resolved_path"] = archiver["resolved_path"]
        ranlib["size"] = archiver["size"]
        ranlib["sha256"] = archiver["sha256"]
        ranlib["mode"] = 0o700
        version = next(item for item in changed["toolchain"]["versions"]
                       if item["role"] == "ranlib")
        version["argv"][0] = archiver["resolved_path"]
        derived_claim_changes.append((
            changed, "give one resolved tool path two modes"))

        changed = copy.deepcopy(record)
        python_tool = next(item for item in changed["toolchain"]["tools"]
                           if item["role"] == "python")
        for runtime_record in changed["runtime_closure"]["records"]:
            runtime_record["dependencies"][1]["path"] = \
                python_tool["resolved_path"]
        changes.append(changed)

        for index, changed in enumerate(changes):
            with self.subTest(index=index):
                self.assertAuthorityRejected(resign(changed))
        for index, (changed, message) in enumerate(derived_claim_changes):
            with self.subTest(derived_claim_index=index, message=message):
                with self.assertRaisesRegex(
                        contract.ExactMainBaselineRecordError, message):
                    contract.validate_baseline_authority_record(
                        resign(changed))

        shared_tool = copy.deepcopy(record)
        archiver = next(item for item in shared_tool["toolchain"]["tools"]
                        if item["role"] == "archiver")
        ranlib = next(item for item in shared_tool["toolchain"]["tools"]
                      if item["role"] == "ranlib")
        ranlib["resolved_path"] = archiver["resolved_path"]
        ranlib["size"] = archiver["size"]
        ranlib["sha256"] = archiver["sha256"]
        version = next(item for item in shared_tool["toolchain"]["versions"]
                       if item["role"] == "ranlib")
        version["argv"][0] = archiver["resolved_path"]
        shared_tool = resign(shared_tool)
        self.assertEqual(
            contract.validate_baseline_authority_record(shared_tool),
            shared_tool,
        )

        shared_blob = copy.deepcopy(record)
        adapter_files = shared_blob["adapter"]["files"]
        adapter_files[1]["git_blob_sha1"] = adapter_files[0]["git_blob_sha1"]
        adapter_files[1]["size"] = adapter_files[0]["size"]
        adapter_files[1]["sha256"] = adapter_files[0]["sha256"]
        shared_blob["adapter"]["files_sha256"] = \
            identity_contract.canonical_sha256(adapter_files)
        shared_blob = resign(shared_blob)
        self.assertEqual(
            contract.validate_baseline_authority_record(shared_blob),
            shared_blob,
        )

    def test_authority_evidence_inventory_is_complete_and_reached(self) -> None:
        expected_paths = [
            "artifacts/canonical-first/leopard_main_benchmark",
            "artifacts/canonical-first/libleopard_main_exact.a",
            "artifacts/canonical-second/leopard_main_benchmark",
            "artifacts/canonical-second/libleopard_main_exact.a",
            "artifacts/path-variant/leopard_main_benchmark",
            "artifacts/path-variant/libleopard_main_exact.a",
            "attestations/canonical-first/benchmark.stderr",
            "attestations/canonical-first/benchmark.stdout.json",
            "attestations/canonical-first/ctest.stderr.log",
            "attestations/canonical-first/ctest.stdout.log",
            "attestations/canonical-second/benchmark.stderr",
            "attestations/canonical-second/benchmark.stdout.json",
            "attestations/canonical-second/ctest.stderr.log",
            "attestations/canonical-second/ctest.stdout.log",
            "attestations/path-variant/benchmark.stderr",
            "attestations/path-variant/benchmark.stdout.json",
            "attestations/path-variant/ctest.stderr.log",
            "attestations/path-variant/ctest.stdout.log",
            "baseline-authority.json",
            "builds/canonical-first/CMakeCache.txt",
            "builds/canonical-first/build-closure.json",
            "builds/canonical-first/build.log",
            "builds/canonical-first/compile_commands.json",
            "builds/canonical-first/configure.log",
            "builds/canonical-second/CMakeCache.txt",
            "builds/canonical-second/build-closure.json",
            "builds/canonical-second/build.log",
            "builds/canonical-second/compile_commands.json",
            "builds/canonical-second/configure.log",
            "builds/path-variant/CMakeCache.txt",
            "builds/path-variant/build-closure.json",
            "builds/path-variant/build.log",
            "builds/path-variant/compile_commands.json",
            "builds/path-variant/configure.log",
            "controllers/test_legacy_main_benchmark.py",
            "logs/00-source_acquisition.log",
            "logs/01-canonical_first_build.log",
            "logs/02-canonical_second_build.log",
            "logs/03-path_variant_build.log",
            "logs/04-independent_verification.log",
            "logs/05-seal.log",
            "runtime/canonical-first/ldd.txt",
            "runtime/canonical-second/ldd.txt",
            "runtime/path-variant/ldd.txt",
            "source/adapter-git-capture.json",
            "source/leopard-main-git-capture.json",
            "source/leopard-main-source.tar",
            "source/leopard2-adapter-source.tar",
            "toolchain/versions/archiver.stderr",
            "toolchain/versions/archiver.stdout",
            "toolchain/versions/assembler.stderr",
            "toolchain/versions/assembler.stdout",
            "toolchain/versions/cc1plus.stderr",
            "toolchain/versions/cc1plus.stdout",
            "toolchain/versions/cmake.stderr",
            "toolchain/versions/cmake.stdout",
            "toolchain/versions/collect2.stderr",
            "toolchain/versions/collect2.stdout",
            "toolchain/versions/compiler.stderr",
            "toolchain/versions/compiler.stdout",
            "toolchain/versions/ctest.stderr",
            "toolchain/versions/ctest.stdout",
            "toolchain/versions/git.stderr",
            "toolchain/versions/git.stdout",
            "toolchain/versions/ldd.stderr",
            "toolchain/versions/ldd.stdout",
            "toolchain/versions/linker.stderr",
            "toolchain/versions/linker.stdout",
            "toolchain/versions/make.stderr",
            "toolchain/versions/make.stdout",
            "toolchain/versions/python.stderr",
            "toolchain/versions/python.stdout",
            "toolchain/versions/ranlib.stderr",
            "toolchain/versions/ranlib.stdout",
        ]
        record = authority_fixture()
        with mock.patch.object(
                contract, "_require_unique_retained_paths",
                wraps=contract._require_unique_retained_paths) as path_guard:
            self.assertEqual(
                contract.validate_baseline_authority_record(record), record)
        self.assertEqual(path_guard.call_count, 1)
        paths, label = path_guard.call_args.args
        self.assertEqual(label, "baseline authority record")
        self.assertEqual(len(paths), 74)
        self.assertEqual(sorted(paths), expected_paths)
        inventory = contract.authority_retained_inventory(record)
        self.assertIs(type(inventory), tuple)
        self.assertEqual(
            [item["relative_path"] for item in inventory], expected_paths)
        self.assertEqual(inventory[18], {
            "relative_path": "baseline-authority.json",
            "size": None,
            "sha256": None,
        })
        detached = contract.authority_retained_inventory(record)
        detached[0]["relative_path"] = "changed"
        self.assertEqual(
            contract.authority_retained_inventory(record)[0]["relative_path"],
            expected_paths[0],
        )

    def test_canonical_ldd_grammar_and_cache_expansion(self) -> None:
        payload = (
            b"ld-linux-x86-64.so.2\tfile\t"
            b"/usr/lib/x86_64-linux-gnu/ld-linux-x86-64.so.2\n"
            b"libc.so.6\tfile\t/usr/lib/x86_64-linux-gnu/libc.so.6\n"
            b"linux-vdso.so.1\tvirtual\n"
        )
        self.assertEqual(contract.parse_canonical_ldd_output(payload), (
            {
                "soname": "ld-linux-x86-64.so.2",
                "kind": "file",
                "path": "/usr/lib/x86_64-linux-gnu/ld-linux-x86-64.so.2",
            },
            {
                "soname": "libc.so.6",
                "kind": "file",
                "path": "/usr/lib/x86_64-linux-gnu/libc.so.6",
            },
            {"soname": "linux-vdso.so.1", "kind": "virtual", "path": None},
        ))
        hostile = (
            payload.replace(b"\n", b"\r\n"),
            payload[:-1],
            payload + payload.splitlines(keepends=True)[0],
            b"libc.so.6\tfile\trelative/path\n",
            b"libc.so.6\tother\n",
            b"libc.so.6\tvirtual\textra\n",
            b"z.so\tvirtual\na.so\tvirtual\n",
            b"libc.so.6\tvirtual\nlibc.so.6\tvirtual\n",
            b"libc.so.6\tvirtual\n\xff",
        )
        for index, changed in enumerate(hostile):
            with self.subTest(ldd_index=index):
                with self.assertRaises(contract.ExactMainBaselineRecordError):
                    contract.parse_canonical_ldd_output(changed)

        roots = roots_fixture("canonical_first")
        requirements = contract.exact_main_build_cache_requirements(roots)
        self.assertIs(type(requirements), tuple)
        self.assertEqual(requirements, (
            {"name": "CMAKE_BUILD_TYPE", "type": "STRING",
             "value": "Release"},
            {"name": "CMAKE_CXX_FLAGS_RELEASE", "type": "STRING",
             "value": "-g -O0 -O3"},
            {"name": "LEOPARD_MAIN_SOURCE_DIR", "type": "PATH",
             "value": roots["baseline_source_root"]},
            {"name": "LEO_MAIN_PURE_AVX2", "type": "BOOL", "value": "ON"},
        ))
        requirements[0]["value"] = "changed"
        self.assertEqual(
            contract.exact_main_build_cache_requirements(roots)[0]["value"],
            "Release")

    def test_runtime_and_attestation_joins(self) -> None:
        record = authority_fixture()
        changes = []
        changed = copy.deepcopy(record)
        changed["runtime_closure"]["records"][0][
            "executable_sha256"] = "f" * 64
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["runtime_closure"]["records"][0]["dependencies"].reverse()
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["runtime_closure"]["records"][0][
            "canonical_ldd_output"] = {
                "relative_path": "runtime/canonical-first/ldd.txt",
                "size": 0,
                "sha256": identity_contract.EMPTY_CONTENT_SHA256,
            }
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["runtime_closure"]["records"][0][
            "canonical_ldd_output"] = text_identity(
                "runtime/canonical-first/ldd.txt", "different-ldd")
        changes.append(changed)
        changed = copy.deepcopy(record)
        dependency = changed["runtime_closure"]["records"][0][
            "dependencies"][1]
        dependency["path"] = changed["runtime_closure"]["records"][0][
            "dependencies"][0]["path"]
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["attestation"]["test_controller"] = {
            "relative_path": "controllers/evil.py",
            "size": 7,
            "sha256": digest("evil-controller"),
        }
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["attestation"]["records"][0]["stdout"] = {
            "relative_path": (
                "attestations/canonical-first/benchmark.stdout.json"),
            "size": 0,
            "sha256": identity_contract.EMPTY_CONTENT_SHA256,
        }
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["attestation"]["records"][0]["ctest"]["stdout"] = {
            "relative_path": (
                "attestations/canonical-first/ctest.stdout.log"),
            "size": 0,
            "sha256": identity_contract.EMPTY_CONTENT_SHA256,
        }
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["attestation"]["records"][0]["pure_avx2"] = False
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["attestation"]["records"][1]["main_source_commit"] = "0" * 40
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["attestation"]["records"][2]["ctest"]["failed"] = 1
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["attestation"]["records"][2]["argv"][-2:] = ["--json", "/tmp/x"]
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["attestation"]["records"][0]["stderr"] = {
            "relative_path": (
                "attestations/canonical-first/benchmark.stderr"),
            "size": 0,
            "sha256": "b" * 64,
        }
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["builds"]["canonical_first"]["build_log"] = {
            "relative_path": "builds/canonical-first/build.log",
            "size": 99,
            "sha256": identity_contract.EMPTY_CONTENT_SHA256,
        }
        changes.append(changed)
        for index, changed in enumerate(changes):
            with self.subTest(index=index):
                self.assertAuthorityRejected(resign(changed))

    def test_failure_variants_are_exact_nonpromotion_records(self) -> None:
        lane = lane_fixture(successful=False, failed_index=2)
        retained = [
            file_identity("logs/00-source.log", "retained-source"),
            file_identity("source/leopard-main-source.tar", "retained-tar"),
        ]
        retained.sort(key=lambda item: item["relative_path"])
        error = {"kind": "build_error", "message": "compiler failed",
                 "exit_status": 2}
        acquisition = contract.baseline_acquisition_failure_record(
            created_utc="2026-08-29T20:31:00Z",
            lane=lane,
            stage=lane["stages"][-1]["name"],
            error=error,
            retained_files=retained,
        )
        self.assertEqual(
            acquisition["schema"],
            "leopard2-gf8-exact-main-pure-avx2-baseline-"
            "acquisition-failure/v1",
        )
        self.assertEqual(
            acquisition["record_sha256"],
            "b457482e5a3fda4ca302bca79b4197f0cacf29080c96450ecb7affd2caa888dd",
        )
        self.assertIsNone(acquisition["authority_record_sha256"])
        self.assertFalse(acquisition["promoted"])
        self.assertEqual(contract.load_baseline_failure_record(
            identity_contract.canonical_json_bytes(acquisition)), acquisition)

        verification_lane = lane_fixture(successful=False, failed_index=4)
        verification = contract.baseline_verification_failure_record(
            created_utc="2026-08-29T20:32:00Z",
            lane=verification_lane,
            stage=verification_lane["stages"][-1]["name"],
            error={"kind": "verification_error",
                   "message": "identity replay failed", "exit_status": 1},
            retained_files=retained,
            authority_record_sha256=authority_fixture()["record_sha256"],
        )
        self.assertEqual(
            verification["schema"],
            "leopard2-gf8-exact-main-pure-avx2-baseline-"
            "verification-failure/v1",
        )
        self.assertEqual(
            verification["record_sha256"],
            "2c86af5066a3611307f8b0133996220726133667a3dae86d9f428f2ee9ad8b15",
        )
        self.assertFalse(verification["promoted"])
        acquisition_inventory = contract.failure_retained_inventory(acquisition)
        expected_failure_paths = sorted(
            ["FAILED.json"] +
            [stage["log"]["relative_path"] for stage in lane["stages"]] +
            [item["relative_path"] for item in retained])
        self.assertEqual(
            [item["relative_path"] for item in acquisition_inventory],
            expected_failure_paths,
        )
        terminal = next(item for item in acquisition_inventory
                        if item["relative_path"] == "FAILED.json")
        self.assertEqual(terminal, {
            "relative_path": "FAILED.json", "size": None, "sha256": None})
        self.assertEqual(
            len(contract.failure_retained_inventory(verification)),
            1 + len(verification_lane["stages"]) + len(retained),
        )

        changes = []
        changed = copy.deepcopy(acquisition)
        changed["promoted"] = True
        changes.append(changed)
        changed = copy.deepcopy(acquisition)
        changed["authority_record_sha256"] = "a" * 64
        changes.append(changed)
        changed = copy.deepcopy(verification)
        changed["schema"] = contract.ACQUISITION_FAILURE_SCHEMA
        changes.append(changed)
        changed = copy.deepcopy(verification)
        changed["stage"] = "source_acquisition"
        changes.append(changed)
        changed = copy.deepcopy(verification)
        changed["retained_files"].reverse()
        changes.append(changed)
        changed = copy.deepcopy(verification)
        changed["schema"] = contract.AUTHORITY_SCHEMA
        changes.append(changed)
        changed = copy.deepcopy(acquisition)
        changed["lane"] = lane_fixture(successful=False, failed_index=4)
        changed["stage"] = "independent_verification"
        changes.append(changed)
        changed = copy.deepcopy(verification)
        changed["lane"] = lane_fixture(successful=False, failed_index=0)
        changed["stage"] = "source_acquisition"
        changes.append(changed)
        changed = copy.deepcopy(acquisition)
        changed["retained_files"][0] = {
            "relative_path": changed["lane"]["stages"][2]["log"][
                "relative_path"],
            "size": 101,
            "sha256": digest("contradictory-stage-log"),
        }
        changed["retained_files"].sort(
            key=lambda item: item["relative_path"])
        changes.append(changed)
        changed = copy.deepcopy(acquisition)
        changed["retained_files"] = [{
            "relative_path": "logs",
            "size": 9,
            "sha256": digest("nested-retained-file"),
        }]
        changes.append(changed)
        for index, changed in enumerate(changes):
            with self.subTest(index=index):
                self.assertFailureRejected(resign(changed))

        unsigned_authority = authority_fixture()
        unsigned_authority["host"]["hostname"] = "unsigned-change"
        self.assertAuthorityRejected(unsigned_authority)
        unsigned_failure = copy.deepcopy(acquisition)
        unsigned_failure["error"]["message"] = "unsigned change"
        self.assertFailureRejected(unsigned_failure)

    def test_static_backstops_are_reached_directly(self) -> None:
        record = authority_fixture()
        changes = []
        changed = copy.deepcopy(record)
        changed["lane"]["stages"] = changed["lane"]["stages"][:-1]
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["lane"]["record_relative_path"] = "FAILED.json"
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["host"]["online_cpus"] = [0, 2, 1, 3]
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["toolchain"]["tools"][0]["mode"] = 0o644
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["builds"]["canonical_first"]["build_argv"][-1] = "-j2"
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["builds"]["canonical_first"]["configure_log"][
            "relative_path"] = "logs/configure.log"
        changes.append(changed)
        changed = copy.deepcopy(record)
        second = changed["builds"]["canonical_second"]
        second["roots"]["build_root"] += "-other"
        second["configure_argv"][4] = second["roots"]["build_root"]
        second["build_argv"][2] = second["roots"]["build_root"]
        changes.append(changed)
        changed = copy.deepcopy(record)
        changed["attestation"]["records"][0]["ctest"]["argv"][-1] = \
            "^different$"
        changes.append(changed)

        changed = copy.deepcopy(record)
        original = changed["identity"]["path_variant"]
        census = copy.deepcopy(original["path_string_census"])
        census["roots"][2]["path"] += "-different"
        changed["identity"]["path_variant"] = \
            identity_contract.normalized_code_identity_record(
                artifact=copy.deepcopy(original["artifact"]),
                sections=copy.deepcopy(original["sections"]),
                path_string_census=census,
            )
        changes.append(changed)

        for index, changed in enumerate(changes):
            with self.subTest(index=index):
                self.assertAuthorityRejected(resign(changed))

        larger_closure = copy.deepcopy(record)
        larger_closure["builds"]["canonical_first"]["closure"][
            "file_count"] = 257
        larger_closure = resign(larger_closure)
        self.assertEqual(
            contract.validate_baseline_authority_record(larger_closure),
            larger_closure,
        )
        oversized_closure = copy.deepcopy(record)
        oversized_closure["builds"]["canonical_first"]["closure"][
            "file_count"] = contract.MAX_CLOSURE_FILES + 1
        self.assertAuthorityRejected(resign(oversized_closure))

    def test_reachable_authority_and_failure_backstops(self) -> None:
        record = authority_fixture()
        authority_changes = []

        changed = copy.deepcopy(record)
        changed["created_utc"] = "2026-02-30T00:00:00Z"
        authority_changes.append((changed, "not a real UTC timestamp"))

        changed = copy.deepcopy(record)
        changed["lane"]["stages"][0]["status"] = "failed"
        authority_changes.append((changed, "lane stage status changed"))

        changed = copy.deepcopy(record)
        changed["lane"]["attempt_budget"] = 4
        authority_changes.append((changed, "attempt budget changed"))

        changed = copy.deepcopy(record)
        changed["host"]["online_cpus"] = []
        authority_changes.append((changed, "online CPU list is invalid"))

        changed = copy.deepcopy(record)
        changed["source"]["adapter_repository"]["clean"] = False
        authority_changes.append((changed, "repository was not clean"))

        changed = copy.deepcopy(record)
        changed["source"]["baseline"]["archive"]["prefix"] = "wrong/"
        authority_changes.append((changed, "source archive prefix changed"))

        changed = copy.deepcopy(record)
        changed["source"]["baseline"]["archive"][
            "relative_path"] = "source/wrong-main-source.tar"
        authority_changes.append((changed, "source archive path changed"))

        changed = copy.deepcopy(record)
        changed["adapter"]["files"] = changed["adapter"]["files"][:-1]
        authority_changes.append((changed, "adapter file inventory changed"))

        changed = copy.deepcopy(record)
        changed["toolchain"]["tools"] = changed["toolchain"]["tools"][:-1]
        authority_changes.append((changed, "tool inventory changed"))

        changed = copy.deepcopy(record)
        changed["toolchain"]["subtools"] = \
            changed["toolchain"]["subtools"][:-1]
        authority_changes.append((changed, "compiler-subtool inventory changed"))

        changed = copy.deepcopy(record)
        changed["toolchain"]["versions"] = \
            changed["toolchain"]["versions"][:-1]
        authority_changes.append((changed, "tool version inventory changed"))

        changed = copy.deepcopy(record)
        changed["toolchain"]["versions"][0]["exit_status"] = 1
        authority_changes.append((changed, "version command failed"))

        changed = copy.deepcopy(record)
        changed["builds"]["canonical_first"]["executable"][
            "retained_relative_path"] = \
            "artifacts/canonical-first/wrong-main"
        authority_changes.append((changed, "executable retained path changed"))

        changed = copy.deepcopy(record)
        changed["builds"]["canonical_first"]["closure"][
            "relative_path"] = "builds/canonical-first/wrong-closure.json"
        authority_changes.append((changed, "closure path changed"))

        changed = copy.deepcopy(record)
        changed["runtime_closure"]["records"] = \
            changed["runtime_closure"]["records"][:-1]
        authority_changes.append((changed, "runtime closure build inventory"))

        changed = copy.deepcopy(record)
        changed["runtime_closure"]["records"][0]["dependencies"][0][
            "kind"] = "unknown"
        authority_changes.append((changed, "dependency 0 kind changed"))

        changed = copy.deepcopy(record)
        virtual = changed["runtime_closure"]["records"][0][
            "dependencies"][2]
        virtual["path"] = "/tmp/not-virtual"
        authority_changes.append((changed, "virtual dependency has file claims"))

        changed = copy.deepcopy(record)
        changed["attestation"]["records"] = \
            changed["attestation"]["records"][:-1]
        authority_changes.append((changed, "attestation build inventory changed"))

        changed = copy.deepcopy(record)
        changed["attestation"]["records"][0]["exit_status"] = 1
        authority_changes.append((changed, "attestation command failed"))

        changed = copy.deepcopy(record)
        changed["attestation"]["records"][0][
            "reported_schema"] = "leopard-main-benchmark-v2"
        authority_changes.append((changed, "benchmark schema changed"))

        changed = copy.deepcopy(record)
        changed["attestation"]["records"][0]["round_trip"] = False
        authority_changes.append((changed, "did not attest correctness"))

        changed = copy.deepcopy(record)
        changed["superseded_references"] = \
            changed["superseded_references"][:-1]
        authority_changes.append((changed, "reference inventory changed"))

        for index, (changed, message) in enumerate(authority_changes):
            with self.subTest(authority_index=index, message=message):
                with self.assertRaisesRegex(
                        contract.ExactMainBaselineRecordError, message):
                    contract.validate_baseline_authority_record(
                        resign(changed))

        failure_lane = lane_fixture(successful=False, failed_index=2)
        failure = contract.baseline_acquisition_failure_record(
            created_utc="2026-08-29T20:31:00Z",
            lane=failure_lane,
            stage=failure_lane["stages"][-1]["name"],
            error={"kind": "build_error", "message": "compiler failed",
                   "exit_status": 2},
            retained_files=[],
        )
        failure_changes = []

        changed = copy.deepcopy(failure)
        changed["created_utc"] = "2026-02-30T00:00:00Z"
        failure_changes.append((changed, "not a real UTC timestamp"))

        changed = copy.deepcopy(failure)
        changed["lane"]["stages"][-1]["status"] = "complete"
        failure_changes.append((changed, "lane stage status changed"))

        changed = copy.deepcopy(failure)
        changed["error"]["kind"] = "bad kind"
        failure_changes.append((changed, "failure kind is not a canonical token"))

        changed = copy.deepcopy(failure)
        changed["error"]["message"] = "first line\nsecond line"
        failure_changes.append((changed, "contains a non-display-safe character"))

        changed = copy.deepcopy(failure)
        changed["error"]["exit_status"] = 0
        failure_changes.append((changed, "exit status is outside"))

        changed = copy.deepcopy(failure)
        changed["retained_files"] = [file_identity(
            f"retained/{index:03d}.bin", f"retained-{index}")
            for index in range(257)]
        failure_changes.append((changed, "file inventory is invalid"))

        changed = copy.deepcopy(failure)
        stage_log = changed["lane"]["stages"][0]["log"]
        changed["retained_files"] = [{
            "relative_path": "retained/conflicting-size.bin",
            "size": stage_log["size"] + 1,
            "sha256": stage_log["sha256"],
        }]
        failure_changes.append((changed, "one SHA-256 at two byte lengths"))

        changed = copy.deepcopy(failure)
        changed["retained_files"] = [file_identity(
            "baseline-authority.json", "failed-authority-draft")]
        failure_changes.append((changed, "retained a promoted terminal name"))

        changed = copy.deepcopy(failure)
        changed["retained_files"] = [file_identity(
            "logs/05-seal.log", "unreachable-seal-log")]
        failure_changes.append((
            changed, "retained a log for a stage that never ran"))

        for index, (changed, message) in enumerate(failure_changes):
            with self.subTest(failure_index=index, message=message):
                with self.assertRaisesRegex(
                        contract.ExactMainBaselineRecordError, message):
                    contract.validate_baseline_failure_record(resign(changed))

        with self.assertRaises(contract.ExactMainBaselineRecordError):
            contract.validate_baseline_failure_record([])
        with self.assertRaises(
                contract.ExactMainBaselineRecordError) as loader_error:
            contract.load_baseline_failure_record(b"{}\n{}\n")
        self.assertIs(type(loader_error.exception),
                      contract.ExactMainBaselineRecordError)

    def test_public_operations_are_pure_under_host_access_canaries(self) -> None:
        tree = ast.parse(
            (TOOLS / "leopard2_exact_main_baseline_record.py").read_text(
                encoding="utf-8"))
        imports = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.update(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                imports.add(node.module)
        self.assertEqual(imports, {
            "__future__", "copy", "datetime", "re", "typing",
            "leopard2_exact_main_baseline",
        })
        with mock.patch.object(
                builtins, "open", side_effect=AssertionError("file I/O")), \
             mock.patch.object(
                 io, "open", side_effect=AssertionError("file I/O")), \
             mock.patch.object(
                 os, "open", side_effect=AssertionError("file I/O")), \
             mock.patch.object(
                 os, "popen", side_effect=AssertionError("process")), \
             mock.patch.object(
                 os, "system", side_effect=AssertionError("host command")), \
             mock.patch.object(
                 subprocess, "run", side_effect=AssertionError("process")), \
             mock.patch.object(
                 subprocess, "Popen", side_effect=AssertionError("process")):
            record = authority_fixture()
            self.assertEqual(
                contract.validate_baseline_authority_record(record), record)
            authority_bytes = contract.canonical_json_bytes(record)
            self.assertEqual(contract.load_baseline_authority_record(
                authority_bytes), record)
            self.assertEqual(contract.exact_main_build_profile(),
                             record["configure"])
            self.assertEqual(len(contract.authority_retained_inventory(record)),
                             74)
            self.assertEqual(
                contract.exact_main_build_cache_requirements(
                    record["builds"]["canonical_first"]["roots"])[2]["value"],
                record["builds"]["canonical_first"]["roots"][
                    "baseline_source_root"],
            )
            self.assertEqual(contract.superseded_historical_references(),
                             record["superseded_references"])

            acquisition_lane = lane_fixture(successful=False)
            acquisition = contract.baseline_acquisition_failure_record(
                created_utc="2026-08-29T20:31:00Z",
                lane=acquisition_lane,
                stage=acquisition_lane["stages"][-1]["name"],
                error={"kind": "build_error", "message": "failed",
                       "exit_status": 1},
                retained_files=[],
            )
            self.assertEqual(
                contract.validate_baseline_failure_record(acquisition),
                acquisition)
            self.assertEqual(
                len(contract.failure_retained_inventory(acquisition)),
                1 + len(acquisition_lane["stages"]),
            )
            self.assertEqual(contract.load_baseline_failure_record(
                contract.canonical_json_bytes(acquisition)), acquisition)

            verification_lane = lane_fixture(
                successful=False, failed_index=4)
            verification = contract.baseline_verification_failure_record(
                created_utc="2026-08-29T20:32:00Z",
                lane=verification_lane,
                stage=verification_lane["stages"][-1]["name"],
                error={"kind": "verification_error",
                       "message": "failed", "exit_status": 1},
                retained_files=[],
                authority_record_sha256=record["record_sha256"],
            )
            self.assertEqual(
                contract.validate_baseline_failure_record(verification),
                verification)
            self.assertEqual(contract.load_baseline_failure_record(
                contract.canonical_json_bytes(verification)), verification)


if __name__ == "__main__":
    unittest.main()
