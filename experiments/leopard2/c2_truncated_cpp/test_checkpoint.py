#!/usr/bin/env python3
"""Fail-closed and adversarial tests for the C2 checkpoint result merger."""

from __future__ import annotations

import importlib.util
import json
import pathlib
import tempfile
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[3]
MODULE_PATH = ROOT / "tools" / "leopard2_c2_checkpoint.py"
SPEC = importlib.util.spec_from_file_location("leopard2_c2_checkpoint", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
checkpoint = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(checkpoint)


def cases() -> list[dict]:
    result = []
    for name, field, parent, shift, direction, active, requested, compared in (
        checkpoint.EXPECTED_CASE_MATRIX
    ):
        slots = max(1, active)
        result.append(
            {
                "name": name,
                "field_bits": field,
                "parent_size": parent,
                "shift": shift,
                "direction": direction,
                "active": active,
                "requested": requested,
                "compared_bytes": compared,
                "pair_operations": 1,
                "complete_operations": 0,
                "maximum_complete_block": 0,
                "scratch_slots": slots,
                "peak_live_slots": slots,
                "execution_scratch_bytes": slots * 64 + 4 * 64,
                "padded_scratch_bytes": parent * (64 + 8),
                "serialized_schedule_bytes": 1024,
                "resident_plan_bytes_excluding_allocator_headers": 2048,
            }
        )
    return result


def result(
    label: str,
    runtime: str,
    source_sha: str,
    core_sha: str,
    core_matrix_sha: str,
    library_sha: str,
    sanitizer: str,
    mode: str = "correctness",
) -> dict:
    return {
        "schema_version": checkpoint.RESULT_SCHEMA,
        "status": "pass",
        "wire_identity": "existing padded dyadic parent",
        "requested_backend": label,
        "runtime_backend": runtime,
        "mode": mode,
        "context_thread_count": 1,
        "shard_bytes": 64,
        "pointer_bytes": 8,
        "pair_staging_shards": 4,
        "runtime_environment": {
            "affinity_api": "sched_getaffinity",
            "process_affinity_cpus": [0, 1],
            "omp_num_threads_env": "1",
            "openmp_max_threads": 1,
            "hostname": "fixture-host",
            "machine": "Linux fixture x86_64",
            "cpu_model": "Fixture CPU",
        },
        "build_binding": {
            "source_sha256": source_sha,
            "core_git_sha": core_sha,
            "core_matrix_sha256": core_matrix_sha,
            "linked_library_sha256": library_sha,
            "sanitizer_mode": sanitizer,
            "compiler": "fixture compiler 1.0",
        },
        "production_multiplier_checks": 65790,
        "correctness_digest_fnv1a64": checkpoint.EXPECTED_CORRECTNESS_DIGEST,
        "cases": cases(),
        "benchmarks": [],
    }


def benchmark_rows(case_by_name: dict[str, dict]) -> list[dict]:
    rows = []
    for name, field, parent, direction in checkpoint.EXPECTED_BENCHMARK_MATRIX:
        case = case_by_name[name]
        rows.append(
            {
                "name": name,
                "field_bits": field,
                "parent_size": parent,
                "direction": direction,
                "repetitions": 11,
                "setup_median_us": 10.0,
                "setup_mad_us": 0.1,
                "candidate_median_us": 100.0,
                "candidate_mad_us": 1.0,
                "padded_median_us": 1.0,
                "padded_mad_us": 0.01,
                "padded_over_candidate": 0.01,
                "candidate_scratch_bytes": case["execution_scratch_bytes"],
                "padded_scratch_bytes": case["padded_scratch_bytes"],
            }
        )
    return rows


class MergeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.directory = tempfile.TemporaryDirectory(prefix="leo2-c2-merge-")
        self.root = pathlib.Path(self.directory.name)
        self.source = ROOT / checkpoint.EXPECTED_SOURCE_RELATIVE
        self.source_sha = checkpoint.sha256(self.source)
        self.core_sha = checkpoint.repository_head(ROOT)
        self.core_matrix_sha = "b" * 64
        self.libraries = {}
        self.library_hashes = {}
        for label in checkpoint.EXPECTED_BACKENDS:
            path = self.root / f"lib-{label}.a"
            path.write_bytes(f"library:{label}\n".encode())
            self.libraries[label] = path
            self.library_hashes[label] = checkpoint.sha256(path)
        self.sanitizer_library = self.root / "lib-asan.a"
        self.sanitizer_library.write_bytes(b"library:asan-ubsan\n")
        sanitizer_library_sha = checkpoint.sha256(self.sanitizer_library)

        self.backend_paths = []
        for label, runtime in checkpoint.EXPECTED_BACKENDS.items():
            path = self.root / f"{label}.json"
            path.write_text(
                json.dumps(
                    result(
                        label,
                        runtime,
                        self.source_sha,
                        self.core_sha,
                        self.core_matrix_sha,
                        self.library_hashes[label],
                        "none",
                    )
                )
                + "\n",
                encoding="utf-8",
            )
            self.backend_paths.append(path)

        self.sanitizer = self.root / "sanitizer.json"
        self.sanitizer.write_text(
            json.dumps(
                result(
                    "asan-ubsan",
                    "avx2",
                    self.source_sha,
                    self.core_sha,
                    self.core_matrix_sha,
                    sanitizer_library_sha,
                    "asan-ubsan",
                )
            )
            + "\n",
            encoding="utf-8",
        )
        benchmark = result(
            "auto",
            "avx2",
            self.source_sha,
            self.core_sha,
            self.core_matrix_sha,
            self.library_hashes["auto"],
            "none",
            "all",
        )
        benchmark["benchmarks"] = benchmark_rows(
            {item["name"]: item for item in benchmark["cases"]}
        )
        benchmark["runtime_environment"]["process_affinity_cpus"] = [3]
        self.benchmark = self.root / "benchmark.json"
        self.benchmark.write_text(json.dumps(benchmark) + "\n", encoding="utf-8")

    def tearDown(self) -> None:
        self.directory.cleanup()

    def merge(self, source: pathlib.Path | None = None):
        return checkpoint.merge(
            self.backend_paths,
            self.sanitizer,
            self.benchmark,
            self.source if source is None else source,
            self.libraries,
            self.sanitizer_library,
            self.core_sha,
            self.core_matrix_sha,
            ROOT,
        )

    def rewrite(self, path: pathlib.Path, update) -> None:
        payload = json.loads(path.read_text(encoding="utf-8"))
        update(payload)
        path.write_text(json.dumps(payload) + "\n", encoding="utf-8")

    def test_valid_matrix_merges(self) -> None:
        payload = self.merge()
        self.assertEqual(payload["status"], "pass")
        self.assertEqual(payload["disposition"], "not_promoted")
        self.assertEqual(payload["correctness"]["cases_per_backend"], 22)
        self.assertEqual(payload["benchmark"]["slowdown_max"], 100.0)

    def test_backend_digest_drift_is_rejected(self) -> None:
        self.rewrite(
            self.backend_paths[-1],
            lambda payload: payload.update(
                correctness_digest_fnv1a64="0xffffffffffffffff"
            ),
        )
        with self.assertRaisesRegex(checkpoint.CheckpointError, "digest differs"):
            self.merge()

    def test_runtime_backend_drift_is_rejected(self) -> None:
        self.rewrite(
            self.backend_paths[1],
            lambda payload: payload.update(runtime_backend="avx2"),
        )
        with self.assertRaisesRegex(checkpoint.CheckpointError, "does not match"):
            self.merge()

    def test_wire_identity_drift_is_rejected(self) -> None:
        self.rewrite(
            self.sanitizer,
            lambda payload: payload.update(wire_identity="new exact profile"),
        )
        with self.assertRaisesRegex(checkpoint.CheckpointError, "wire identity"):
            self.merge()

    def test_missing_benchmark_is_rejected(self) -> None:
        self.rewrite(
            self.benchmark, lambda payload: payload.update(benchmarks=[])
        )
        with self.assertRaisesRegex(checkpoint.CheckpointError, "row count"):
            self.merge()

    def test_unrelated_source_is_rejected(self) -> None:
        unrelated = self.root / "unrelated.cpp"
        unrelated.write_text("int unrelated;\n", encoding="utf-8")
        with self.assertRaisesRegex(checkpoint.CheckpointError, "selected source"):
            self.merge(unrelated)

    def test_changed_avx2_scratch_is_rejected(self) -> None:
        self.rewrite(
            self.backend_paths[-1],
            lambda payload: payload["cases"][0].update(
                execution_scratch_bytes=(
                    payload["cases"][0]["execution_scratch_bytes"] + 64
                )
            ),
        )
        with self.assertRaisesRegex(checkpoint.CheckpointError, "scratch bound"):
            self.merge()

    def test_fake_thousand_x_row_is_rejected(self) -> None:
        def fake(payload):
            payload["benchmarks"][0].update(
                candidate_median_us=1000.0,
                padded_median_us=1.0,
                padded_over_candidate=0.001,
            )

        self.rewrite(self.benchmark, fake)
        with self.assertRaisesRegex(checkpoint.CheckpointError, "runtime ratio"):
            self.merge()

    def test_scalar_sanitizer_is_rejected(self) -> None:
        self.rewrite(
            self.sanitizer,
            lambda payload: payload.update(runtime_backend="scalar"),
        )
        with self.assertRaisesRegex(checkpoint.CheckpointError, "sanitizer runtime"):
            self.merge()

    def test_unpinned_benchmark_is_rejected(self) -> None:
        self.rewrite(
            self.benchmark,
            lambda payload: payload["runtime_environment"].update(
                process_affinity_cpus=[2, 3]
            ),
        )
        with self.assertRaisesRegex(checkpoint.CheckpointError, "not actually pinned"):
            self.merge()

    def test_unbound_openmp_environment_is_rejected(self) -> None:
        self.rewrite(
            self.backend_paths[0],
            lambda payload: payload["runtime_environment"].update(
                omp_num_threads_env="unset"
            ),
        )
        with self.assertRaisesRegex(checkpoint.CheckpointError, "OMP_NUM_THREADS"):
            self.merge()

    def test_changed_library_binding_is_rejected(self) -> None:
        self.rewrite(
            self.backend_paths[0],
            lambda payload: payload["build_binding"].update(
                linked_library_sha256="0" * 64
            ),
        )
        with self.assertRaisesRegex(checkpoint.CheckpointError, "selected library"):
            self.merge()


if __name__ == "__main__":
    unittest.main()
