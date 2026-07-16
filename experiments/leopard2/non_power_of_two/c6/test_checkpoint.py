#!/usr/bin/env python3
"""Adversarial tests for the C6 fail-closed evidence merger."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"


def load_checkpoint():
    spec = importlib.util.spec_from_file_location("_leopard2_c6_checkpoint",
                                                  HERE / "checkpoint.py")
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load checkpoint")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


checkpoint = load_checkpoint()


class CheckpointTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.names = (
            "algebra", "auto", "scalar", "ssse3", "avx2", "asan-ubsan",
            "benchmark", "benchmark-manifest", "benchmark-stdout",
            "benchmark-stderr",
        )
        filenames = {
            "algebra": "algebra.json",
            "auto": "auto.json",
            "scalar": "scalar.json",
            "ssse3": "ssse3.json",
            "avx2": "avx2.json",
            "asan-ubsan": "asan-ubsan.json",
            "benchmark": "benchmark.json",
            "benchmark-manifest": "benchmark-manifest.json",
            "benchmark-stdout": "benchmark.stdout.txt",
            "benchmark-stderr": "benchmark.stderr.txt",
        }
        cls.contents = {name: (RESULTS / filename).read_bytes()
                        for name, filename in filenames.items()}
        cls.algebra = json.loads(cls.contents["algebra"])

    def materialize(self, changes=None):
        changes = changes or {}
        temporary = tempfile.TemporaryDirectory(prefix="leopard2-c6-test-")
        root = Path(temporary.name)
        paths = {}
        for name in self.names:
            path = root / name
            payload = changes.get(name, self.contents[name])
            if not isinstance(payload, bytes):
                payload = (json.dumps(payload, indent=2, sort_keys=True,
                                      allow_nan=True) + "\n").encode("utf-8")
            path.write_bytes(payload)
            paths[name] = path
        return temporary, paths

    def coherent_benchmark_change(self, mutate):
        benchmark = json.loads(self.contents["benchmark"])
        mutate(benchmark)
        benchmark_bytes = (json.dumps(benchmark, indent=2, sort_keys=True,
                                      allow_nan=True) + "\n").encode("utf-8")
        manifest = json.loads(self.contents["benchmark-manifest"])
        manifest["result_sha256"] = hashlib.sha256(benchmark_bytes).hexdigest()
        return {"benchmark": benchmark_bytes, "benchmark-manifest": manifest}

    def expect_reject(self, changes):
        temporary, paths = self.materialize(changes)
        try:
            fake = types.SimpleNamespace(run=lambda: copy.deepcopy(self.algebra))
            with mock.patch.object(checkpoint, "load_algebra_module", return_value=fake):
                with self.assertRaises(checkpoint.EvidenceError):
                    checkpoint.merge(paths)
        finally:
            temporary.cleanup()

    def test_happy_path_recomputes_real_algebra(self):
        temporary, paths = self.materialize()
        try:
            result = checkpoint.merge(paths)
            self.assertEqual(result["status"], "pass")
            self.assertEqual(result["geometry"]["benchmark_cells"], 50)
            self.assertEqual(result["geometry"]["decode_benchmark_cells"], 56)
            self.assertTrue(result["benchmark"]["meaningful_region_gate"])
            self.assertEqual(result["decode_benchmark"]["maximum_loss_cells"], 36)
            self.assertGreater(
                result["decode_benchmark"][
                    "one_shot_setup_plus_execution_ratio_min"], 1.0)
            self.assertGreater(
                result["decode_benchmark"][
                    "declared_reuse_setup_plus_execution_ratio_min"], 1.0)
            self.assertEqual(result["memory"]["exact_decode_execution_scratch_bytes"], 0)
            self.assertEqual(result["correctness"]["gf8_decoder_algebra"], {
                "plans": 16,
                "folded_nonzero_terms": 4933,
                "recovered_values": 70,
                "no_loss_execution_terms": 0,
            })
            self.assertEqual(
                result["correctness"]["cpp_per_backend"]["hot_path_allocations"], 0)
            self.assertEqual(
                result["correctness"]["cpp_per_backend"]["no_loss_calls"], 96)
            self.assertEqual(
                result["correctness"]["cpp_per_backend"]["maximum_loss_cases"], 96)
            self.assertEqual(result["disposition"]["production_promotion"], "none")
            self.assertEqual(
                result["disposition"]["bead_recommendation"],
                "close-experiment-result",
            )
            raw = json.loads(self.contents["benchmark"])
            gains = sorted(
                cell["credible_gain_percent"] for cell in raw["cells"]
                if min(cell["K"], cell["R"]) == 3
            )
            even_median = (gains[11] + gains[12]) * 0.5
            self.assertEqual(len(gains), 24)
            self.assertAlmostEqual(
                result["benchmark"]["target_credible_gain_percent_median"],
                even_median,
            )
            self.assertNotEqual(even_median, gains[12])
        finally:
            temporary.cleanup()

    def test_algebra_mutation_rejected(self):
        algebra = copy.deepcopy(self.algebra)
        algebra["gf8"]["boundary_cases"] -= 1
        self.expect_reject({"algebra": algebra})

    def test_decoder_algebra_mutation_rejected(self):
        algebra = copy.deepcopy(self.algebra)
        algebra["decoder_plan_records"][0]["term_count"] += 1
        self.expect_reject({"algebra": algebra})

    def test_missing_duplicate_extra_and_relabelled_cells_rejected(self):
        mutations = (
            lambda data: data["cells"].pop(),
            lambda data: data["cells"].append(copy.deepcopy(data["cells"][0])),
            lambda data: data["cells"][0].__setitem__("bytes", 128),
            lambda data: data["cells"][0].__setitem__("profile", "low_v1"),
        )
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                self.expect_reject(self.coherent_benchmark_change(mutation))

    def test_coordinated_derived_timing_forgery_rejected(self):
        def mutate(data):
            data["cells"][0]["padded_over_exact_ratio"] += 0.01
            data["cells"][0]["credible_gain_percent"] += 1.0
        self.expect_reject(self.coherent_benchmark_change(mutate))

    def test_missing_duplicate_extra_and_relabelled_decode_cells_rejected(self):
        mutations = (
            lambda data: data["decode_cells"].pop(),
            lambda data: data["decode_cells"].append(
                copy.deepcopy(data["decode_cells"][0])),
            lambda data: data["decode_cells"][0].__setitem__("bytes", 128),
            lambda data: data["decode_cells"][0].__setitem__("profile", "low_v1"),
            lambda data: data["decode_cells"][0].__setitem__("losses", 1),
        )
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                self.expect_reject(self.coherent_benchmark_change(mutation))

    def test_coordinated_decode_derived_timing_forgery_rejected(self):
        def mutate(data):
            data["decode_cells"][0]["padded_over_exact_ratio"] += 0.01
            data["decode_cells"][0]["credible_gain_percent"] += 1.0
        self.expect_reject(self.coherent_benchmark_change(mutate))

    def test_nonfinite_or_negative_uncertainty_rejected(self):
        for value in (float("nan"), -1.0):
            def mutate(data, replacement=value):
                data["cells"][0]["exact_execution"]["mad_us"] = replacement
            with self.subTest(value=value):
                self.expect_reject(self.coherent_benchmark_change(mutate))

    def test_decode_nonfinite_or_negative_uncertainty_rejected(self):
        for value in (float("nan"), -1.0):
            def mutate(data, replacement=value):
                data["decode_cells"][0]["exact_execution"]["mad_us"] = replacement
            with self.subTest(value=value):
                self.expect_reject(self.coherent_benchmark_change(mutate))

    def test_accounting_forgery_rejected(self):
        fields = (
            "exact_table_bytes", "exact_execution_scratch_bytes",
            "padded_execution_scratch_bytes", "exact_fixed_multiply_terms",
            "exact_coefficient_payload_bytes", "input_bytes", "output_bytes",
        )
        for field in fields:
            def mutate(data, name=field):
                data["cells"][0][name] += 1
            with self.subTest(field=field):
                self.expect_reject(self.coherent_benchmark_change(mutate))

    def test_decode_accounting_forgery_rejected(self):
        numeric_fields = (
            "selected_parity_count", "selected_parity_prefix_end",
            "exact_codec_table_bytes", "exact_plan_terms",
            "exact_term_record_bytes", "exact_offset_record_bytes",
            "exact_coordinate_record_bytes", "exact_plan_payload_bytes",
            "exact_execution_terms", "exact_term_payload_bytes",
            "exact_execution_scratch_bytes", "padded_execution_scratch_bytes",
            "offered_received_bytes", "repaired_output_bytes",
        )
        for field in numeric_fields:
            def mutate(data, name=field):
                data["decode_cells"][0][name] += 1
            with self.subTest(field=field):
                self.expect_reject(self.coherent_benchmark_change(mutate))

        def mutate_baseline(data):
            data["decode_cells"][0]["baseline"] = "padded_generic_gf16"
        self.expect_reject(self.coherent_benchmark_change(mutate_baseline))

    def test_profile_output_identity_forgery_rejected(self):
        def mutate(data):
            data["cells"][0]["padded_output_digest"] = \
                data["cells"][0]["exact_output_digest"]
        self.expect_reject(self.coherent_benchmark_change(mutate))

    def test_decode_output_difference_forgery_rejected(self):
        def mutate(data):
            data["decode_cells"][0]["padded_output_digest"] = \
                "0x0000000000000000"
        self.expect_reject(self.coherent_benchmark_change(mutate))

    def test_affinity_backend_and_threading_forgery_rejected(self):
        mutations = (
            lambda data: data.__setitem__("affinity", [15, 31]),
            lambda data: data.__setitem__("runtime_backend", "scalar"),
            lambda data: data.__setitem__("requested_backend", "auto"),
            lambda data: data.__setitem__(
                "allocation_tracking", "disabled-for-sanitizer"),
            lambda data: data.__setitem__("omp_num_threads", "2"),
            lambda data: data.__setitem__("openmp_max_threads", 2),
        )
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                self.expect_reject(self.coherent_benchmark_change(mutation))

    def test_source_core_and_correctness_binding_rejected(self):
        for field in ("source_sha256", "core_git_sha", "library_sha256"):
            def mutate(data, name=field):
                data[name] = "0" * 64
            with self.subTest(field=field):
                self.expect_reject(self.coherent_benchmark_change(mutate))
        def mutate_correctness(data):
            data["correctness"]["digest"] = "0x0000000000000000"
        self.expect_reject(self.coherent_benchmark_change(mutate_correctness))
        for field, replacement in (
            ("decode_digest", "0x0000000000000000"),
            ("hot_path_allocations", 1),
        ):
            def mutate_decode_correctness(data, name=field, value=replacement):
                data["correctness"][name] = value
            with self.subTest(field=field):
                self.expect_reject(
                    self.coherent_benchmark_change(mutate_decode_correctness))

    def test_manifest_executable_result_and_log_binding_rejected(self):
        for field in (
            "executable_sha256", "result_sha256", "runner_sha256",
            "source_sha256", "library_sha256", "stdout_sha256", "stderr_sha256",
        ):
            manifest = json.loads(self.contents["benchmark-manifest"])
            manifest[field] = "0" * 64
            with self.subTest(field=field):
                self.expect_reject({"benchmark-manifest": manifest})

        for field in ("encode_cell_count", "decode_cell_count"):
            manifest = json.loads(self.contents["benchmark-manifest"])
            manifest[field] -= 1
            with self.subTest(field=field):
                self.expect_reject({"benchmark-manifest": manifest})

    def test_backend_matrix_and_sanitizer_labels_rejected(self):
        scalar = json.loads(self.contents["scalar"])
        scalar["runtime_backend"] = "avx2"
        self.expect_reject({"scalar": scalar})
        scalar = json.loads(self.contents["scalar"])
        scalar["allocation_tracking"] = "disabled-for-sanitizer"
        self.expect_reject({"scalar": scalar})
        sanitizer = json.loads(self.contents["asan-ubsan"])
        sanitizer["sanitizer_mode"] = "none"
        self.expect_reject({"asan-ubsan": sanitizer})
        sanitizer = json.loads(self.contents["asan-ubsan"])
        sanitizer["allocation_tracking"] = "global-new"
        self.expect_reject({"asan-ubsan": sanitizer})

    def test_log_tamper_rejected_even_with_unchanged_manifest(self):
        self.expect_reject({"benchmark-stderr": self.contents["benchmark-stderr"] + b"tamper\n"})


if __name__ == "__main__":
    unittest.main()
