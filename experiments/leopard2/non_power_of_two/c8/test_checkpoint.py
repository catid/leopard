#!/usr/bin/env python3
"""Fail-closed unit tests for the Leopard2 C8 evidence merger."""

from __future__ import annotations

import copy
import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SPEC = importlib.util.spec_from_file_location("c8_checkpoint", HERE / "checkpoint.py")
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load the C8 checkpoint module")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def core_revision() -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value), encoding="utf-8")


def benchmark_payload(source_hash: str, library_hash: str) -> dict:
    rows = []
    for field, k, r, byte_count, batch, reuse in MODULE.EXPECTED_BENCHMARKS:
        exact = 10.
        exact_mad = 0.
        padded = 12.
        padded_mad = 0.
        rows.append({
            "field": field, "K": k, "R": r, "bytes": byte_count,
            "batch": batch, "reuse": reuse,
            "exact_median_us": exact, "exact_mad_us": exact_mad,
            "padded_median_us": padded, "padded_mad_us": padded_mad,
            "ratio": round(padded / exact, 6),
            "credible_gain": round((padded - padded_mad) /
                                    (exact + exact_mad) - 1., 6),
            "exact_setup_median_us": 20., "exact_setup_mad_us": 0.,
            "padded_setup_median_us": 2., "padded_setup_mad_us": 0.,
            "exact_table_bytes": k * r * (1 if field == "gf8" else 2),
            "exact_scratch_bytes": 0, "padded_scratch_bytes": 64,
            "exact_logical_bytes": r * (1 + 3 * (k - 1)) * byte_count,
            "exact_samples_us": [exact] * 11,
            "padded_samples_us": [padded] * 11,
            "exact_setup_samples_us": [20.] * 7,
            "padded_setup_samples_us": [2.] * 7,
            "digest": "0x1",
        })
    return {
        "schema": "leopard2-c8-executable-v1", "mode": "benchmark",
        "backend_label": "avx2", "runtime_backend": "avx2",
        "source_sha256": source_hash, "core_git_sha": core_revision(),
        "library_sha256": library_hash, "sanitizer": "none",
        "allocation_tracking": "global-new",
        "profile": "exact_high_prefix_v1_candidate", "default_off": True,
        "openmp_max_threads": 1, "affinity": [15], "benchmarks": rows,
    }


def executable_payload(source_hash: str, library_hash: str,
                       label: str = "scalar") -> dict:
    return {
        "schema": "leopard2-c8-executable-v1", "mode": "correctness",
        "backend_label": label, "runtime_backend": label,
        "source_sha256": source_hash, "core_git_sha": core_revision(),
        "library_sha256": library_hash, "sanitizer": "none",
        "allocation_tracking": "global-new",
        "profile": "exact_high_prefix_v1_candidate", "default_off": True,
        "openmp_max_threads": 1, "affinity": [0, 1],
        "correctness": copy.deepcopy(MODULE.EXPECTED_CORRECTNESS),
        "benchmarks": [],
    }


def isolation_payload() -> dict:
    fields = {name: 0 for name in
              ("user", "nice", "system", "idle", "iowait", "irq",
               "softirq", "steal", "guest", "guest_nice")}
    timing = dict(fields)
    timing["user"] = 1
    return {
        "schema": "leopard2-c8-isolation-v1", "exit_code": 0,
        "timing_cpu": 15, "sibling_cpu": 31, "sibling_idle": True,
        "sibling_nonidle_jiffies": 0, "elapsed_seconds": .1,
        "environment": {"OMP_DYNAMIC": "FALSE", "OMP_NUM_THREADS": "1"},
        "command": ["c8", "--mode", "benchmark", "--backend-label", "avx2"],
        "delta": {"15": timing, "31": fields},
    }


class CheckpointTests(unittest.TestCase):
    def test_parse_mapping_rejects_duplicate(self) -> None:
        with self.assertRaises(ValueError):
            MODULE.parse_mapping(["auto=a", "auto=b"], "--backend")

    def test_benchmark_accepts_complete_matrix(self) -> None:
        with tempfile.TemporaryDirectory() as directory_name:
            directory = Path(directory_name)
            library = directory / "library.a"
            library.write_bytes(b"library")
            source_hash = "a" * 64
            path = directory / "benchmark.json"
            write_json(path, benchmark_payload(source_hash, MODULE.sha256(library)))
            _, rows = MODULE.validate_benchmark(path, source_hash, library, ROOT)
            self.assertEqual(len(rows), len(MODULE.EXPECTED_BENCHMARKS))

    def test_benchmark_rejects_forged_fields(self) -> None:
        mutations = (
            lambda value: value["benchmarks"][0].update(ratio=99.),
            lambda value: value["benchmarks"].pop(),
            lambda value: value["benchmarks"][0].update(exact_table_bytes=0),
            lambda value: value["benchmarks"][0].update(
                exact_samples_us=[99.] * 11),
            lambda value: value["benchmarks"][0][
                "padded_setup_samples_us"].pop(),
            lambda value: value.update(affinity=[15, 31]),
            lambda value: value.update(runtime_backend="scalar"),
        )
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                with tempfile.TemporaryDirectory() as directory_name:
                    directory = Path(directory_name)
                    library = directory / "library.a"
                    library.write_bytes(b"library")
                    source_hash = "a" * 64
                    value = benchmark_payload(source_hash, MODULE.sha256(library))
                    mutation(value)
                    path = directory / "benchmark.json"
                    write_json(path, value)
                    with self.assertRaises(ValueError):
                        MODULE.validate_benchmark(path, source_hash, library, ROOT)

    def test_executable_accepts_projection(self) -> None:
        with tempfile.TemporaryDirectory() as directory_name:
            directory = Path(directory_name)
            library = directory / "library.a"
            library.write_bytes(b"library")
            source_hash = "b" * 64
            path = directory / "scalar.json"
            write_json(path, executable_payload(source_hash, MODULE.sha256(library)))
            MODULE.validate_executable(path, source_hash, library, "scalar", ROOT)

    def test_executable_rejects_digest_backend_and_source(self) -> None:
        mutations = (
            lambda value: value["correctness"].update(digest="0x0"),
            lambda value: value.update(runtime_backend="avx2"),
            lambda value: value.update(source_sha256="c" * 64),
        )
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                with tempfile.TemporaryDirectory() as directory_name:
                    directory = Path(directory_name)
                    library = directory / "library.a"
                    library.write_bytes(b"library")
                    source_hash = "b" * 64
                    value = executable_payload(source_hash, MODULE.sha256(library))
                    mutation(value)
                    path = directory / "scalar.json"
                    write_json(path, value)
                    with self.assertRaises(ValueError):
                        MODULE.validate_executable(
                            path, source_hash, library, "scalar", ROOT)

    def test_isolation_accepts_idle_sibling(self) -> None:
        with tempfile.TemporaryDirectory() as directory_name:
            path = Path(directory_name) / "isolation.json"
            write_json(path, isolation_payload())
            MODULE.validate_isolation(path, 15)

    def test_isolation_rejects_sibling_activity(self) -> None:
        with tempfile.TemporaryDirectory() as directory_name:
            value = isolation_payload()
            value["sibling_idle"] = False
            value["sibling_nonidle_jiffies"] = 1
            path = Path(directory_name) / "isolation.json"
            write_json(path, value)
            with self.assertRaises(ValueError):
                MODULE.validate_isolation(path, 15)

    def test_default_build_is_unmodified(self) -> None:
        MODULE.validate_default_off(ROOT)

    def test_sha256_is_content_bound(self) -> None:
        with tempfile.TemporaryDirectory() as directory_name:
            path = Path(directory_name) / "data"
            path.write_bytes(b"c8\n")
            first = MODULE.sha256(path)
            path.write_bytes(b"C8\n")
            self.assertNotEqual(first, MODULE.sha256(path))


if __name__ == "__main__":
    unittest.main()
