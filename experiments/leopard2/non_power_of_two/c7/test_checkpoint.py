#!/usr/bin/env python3
"""Portable validation for the C7 exact-low correctness checkpoint."""

from __future__ import annotations

import copy
import hashlib
import json
import pathlib
import unittest


HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[3]
RESULTS = HERE / "results"
BACKENDS = ("scalar", "ssse3", "avx2", "auto")
EXPECTED_RUNTIME = {
    "scalar": "scalar", "ssse3": "ssse3", "avx2": "avx2", "auto": "avx2"
}
EXPECTED_CORRECTNESS = {
    "gf8_cases": 9,
    "gf16_cases": 5,
    "coefficients": 118717,
    "encode_executions": 117,
    "encode_symbol_comparisons": 1030423,
    "subset_encode_executions": 117,
    "decode_plans": 403,
    "decode_executions": 403,
    "decode_symbol_comparisons": 272487,
    "maximum_loss_plans": 117,
    "unavailable_parity_plans": 175,
    "no_loss_null_calls": 14,
    "parity_rebuilds": 403,
    "odd_gf16_rejections": 10,
    "overlap_rejections": 28,
    "read_only_input_alias_calls": 13,
    "hot_path_allocations": 0,
    "digest_fnv64": "0xec4179e9f2776a58",
}


def load(name: str) -> dict:
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def validate_algebra(data: dict) -> None:
    if set(data) != {
        "coordinate_comparison", "derivation", "gf4_exhaustive",
        "large_field", "profile", "schema", "source_sha256", "status"
    }:
        raise ValueError("unexpected algebra keys")
    if data["schema"] != "leopard2-c7-algebra/v1" or data["status"] != "pass":
        raise ValueError("bad algebra status")
    profile = data["profile"]
    if profile != {
        "identity_family": 3,
        "profile_version": 1,
        "coordinate_map_version": 1,
        "systematic_coordinates": "omega_0 .. omega_(K-1)",
        "parity_coordinates": "omega_K .. omega_(K+R-1)",
        "parent_count": "K+R",
        "padded_side": "K",
        "shortening": "none",
        "puncturing": "none",
        "future_exact_high_identity_family": 4,
    }:
        raise ValueError("algebra profile changed")
    source_hash = hashlib.sha256((HERE / "algebra.py").read_bytes()).hexdigest()
    if data["source_sha256"] != source_hash:
        raise ValueError("algebra source hash mismatch")
    small = data["gf4_exhaustive"]
    expected = {
        "geometries": 120,
        "mds_subsets": 131038,
        "dense_coefficients": 3060,
        "vandermonde_oracle_coefficients": 3060,
        "affine_maps": 28800,
        "affine_coefficients": 734400,
        "searched_full_field_partitions": 65534,
        "transform_search_candidates": 65519,
        "transform_search_wins": 57,
    }
    for key, value in expected.items():
        if small.get(key) != value:
            raise ValueError(f"algebra counter changed: {key}")
    aligned = small.get("aligned_union", {})
    if aligned != {
        "available_geometries": 85,
        "globally_affine_geometries": 49,
        "nonaffine_changed_coefficients": 758,
        "dense_coefficients": 1812,
        "prefix_one_coefficients": 297,
        "aligned_one_coefficients": 198,
        "prefix_symbolic_fragments": 358,
        "aligned_symbolic_fragments": 328,
    }:
        raise ValueError("aligned-union comparison changed")
    if len(small.get("search", [])) != 15 or len(
            small.get("transform_search_focus", [])) != 5:
        raise ValueError("coordinate search records are incomplete")
    large = data["large_field"]
    if (large.get("cases"), large.get("dense_coefficients"),
            large.get("affine_coefficients")) != (9, 83281, 177600):
        raise ValueError("large-field counters changed")
    if len(large.get("records", [])) != 9 or len(
            large.get("aligned_union_records", [])) != 9:
        raise ValueError("large-field records are incomplete")
    gf8_boundary = large["aligned_union_records"][0]
    if gf8_boundary.get("available") is not False:
        raise ValueError("aligned map must not cross the GF8 field boundary")
    gf16_nonaffine = large["aligned_union_records"][6]
    if (gf16_nonaffine.get("globally_affine") is not False or
            gf16_nonaffine.get("changed_coefficients") != 12900):
        raise ValueError("non-affine aligned GF16 witness changed")
    comparison = data["coordinate_comparison"]
    if comparison.get("decision") != (
            "freeze the simple prefix map; do not add search tables to V1"):
        raise ValueError("coordinate-map decision changed")


def validate_correctness(data: dict, backend: str, sanitizer: bool = False) -> None:
    required = {
        "allocation_tracking", "benchmarks", "core_git_sha", "correctness",
        "library_sha256", "production_constructor_rejected", "profile",
        "requested_backend", "runtime_backend", "sanitizer", "schema",
        "source_sha256", "status", "timing_scope"
    }
    if set(data) != required:
        raise ValueError("unexpected C++ result keys")
    if data["schema"] != "leopard2-c7-exact-low/v1" or data["status"] != "pass":
        raise ValueError("bad C++ status")
    if data["profile"] != {
        "family": 3, "version": 1, "coordinate_map": 1,
        "systematic": "0..K-1", "parity": "K..K+R-1",
        "production_enabled": False,
    }:
        raise ValueError("C++ profile changed")
    if data["production_constructor_rejected"] is not True:
        raise ValueError("production unexpectedly enabled C7")
    if data["timing_scope"] != "none-correctness-only" or data["benchmarks"]:
        raise ValueError("correctness result contains timing evidence")
    if data["correctness"] != EXPECTED_CORRECTNESS:
        raise ValueError("C++ correctness counters changed")
    source_hash = hashlib.sha256((HERE / "c7_exact_low.cpp").read_bytes()).hexdigest()
    if data["source_sha256"] != source_hash:
        raise ValueError("C++ source hash mismatch")
    if data["core_git_sha"] != "48803c06fbd7a6802b4438af60e3104895938c9d":
        raise ValueError("unexpected linked core commit")
    if len(data["library_sha256"]) != 64:
        raise ValueError("missing linked-library hash")
    if sanitizer:
        if (data["requested_backend"] != "auto" or
                data["runtime_backend"] != "avx2" or
                data["sanitizer"] != "asan-ubsan" or
                data["allocation_tracking"] != "disabled-for-sanitizer"):
            raise ValueError("bad sanitizer provenance")
    elif (data["requested_backend"] != backend or
          data["runtime_backend"] != EXPECTED_RUNTIME[backend] or
          data["sanitizer"] != "none" or
          data["allocation_tracking"] != "global-new"):
        raise ValueError("bad backend provenance")


def validate_smoke(data: dict) -> None:
    if data.get("timing_scope") != "non-authoritative-smoke":
        raise ValueError("smoke timing was not labelled non-authoritative")
    if data.get("correctness") != EXPECTED_CORRECTNESS:
        raise ValueError("smoke correctness differs")
    cells = data.get("benchmarks")
    if not isinstance(cells, list) or len(cells) != 1:
        raise ValueError("smoke must contain exactly one cell")
    cell = cells[0]
    if [cell.get(key) for key in ("K", "R", "bytes", "batch", "losses")] != [
            3, 253, 64, 8, 3]:
        raise ValueError("smoke geometry changed")
    for key in ("exact_encode_samples_us", "padded_encode_samples_us",
                "exact_decode_samples_us", "padded_decode_samples_us"):
        samples = cell.get(key)
        if not isinstance(samples, list) or len(samples) != 7 or not all(
                isinstance(value, (int, float)) and value > 0 for value in samples):
            raise ValueError("smoke raw samples are incomplete")
    if cell.get("exact_coefficients") != 759 or cell.get("exact_decode_terms") != 9:
        raise ValueError("smoke exact plan accounting changed")


class CheckpointTests(unittest.TestCase):
    def test_algebra(self) -> None:
        validate_algebra(load("algebra.json"))

    def test_backends(self) -> None:
        hashes = set()
        for backend in BACKENDS:
            data = load(f"{backend}.json")
            validate_correctness(data, backend)
            hashes.add(data["correctness"]["digest_fnv64"])
        self.assertEqual(hashes, {"0xec4179e9f2776a58"})

    def test_sanitizer(self) -> None:
        validate_correctness(load("asan-ubsan.json"), "auto", sanitizer=True)

    def test_smoke_is_non_authoritative(self) -> None:
        validate_smoke(load("smoke-nonauthoritative.json"))

    def test_algebra_mutations_rejected(self) -> None:
        original = load("algebra.json")
        mutations = []
        for path, value in (
            (("profile", "identity_family"), 4),
            (("gf4_exhaustive", "mds_subsets"), 1),
            (("gf4_exhaustive", "aligned_union", "globally_affine_geometries"), 85),
            (("large_field", "dense_coefficients"), 0),
            (("coordinate_comparison", "decision"), "aligned"),
            (("source_sha256",), "0" * 64),
        ):
            candidate = copy.deepcopy(original)
            cursor = candidate
            for key in path[:-1]:
                cursor = cursor[key]
            cursor[path[-1]] = value
            mutations.append(candidate)
        for candidate in mutations:
            with self.assertRaises(ValueError):
                validate_algebra(candidate)

    def test_cpp_mutations_rejected(self) -> None:
        original = load("auto.json")
        mutations = []
        for path, value in (
            (("profile", "coordinate_map"), 2),
            (("production_constructor_rejected",), False),
            (("correctness", "coefficients"), 0),
            (("correctness", "digest_fnv64"), "0x0"),
            (("correctness", "hot_path_allocations"), 1),
            (("timing_scope",), "authoritative"),
            (("source_sha256",), "0" * 64),
            (("runtime_backend",), "scalar"),
        ):
            candidate = copy.deepcopy(original)
            cursor = candidate
            for key in path[:-1]:
                cursor = cursor[key]
            cursor[path[-1]] = value
            mutations.append(candidate)
        for candidate in mutations:
            with self.assertRaises(ValueError):
                validate_correctness(candidate, "auto")

    def test_smoke_mutations_rejected(self) -> None:
        original = load("smoke-nonauthoritative.json")
        for mutation in ("scope", "samples", "geometry"):
            candidate = copy.deepcopy(original)
            if mutation == "scope":
                candidate["timing_scope"] = "authoritative"
            elif mutation == "samples":
                candidate["benchmarks"][0]["exact_encode_samples_us"] = [1]
            else:
                candidate["benchmarks"][0]["K"] = 4
            with self.assertRaises(ValueError):
                validate_smoke(candidate)


if __name__ == "__main__":
    unittest.main(verbosity=2)
