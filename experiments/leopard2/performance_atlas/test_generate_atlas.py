#!/usr/bin/env python3
"""Focused fail-closed tests for the performance-atlas runner and plotter."""

from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


MODULE_PATH = Path(__file__).with_name("generate_atlas.py")
SPEC = importlib.util.spec_from_file_location("leopard2_performance_atlas", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
ATLAS = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = ATLAS
SPEC.loader.exec_module(ATLAS)


def samples(iterations: int, batch: bool = False) -> dict:
    key = "median_us_per_batch_call" if batch else "median_us"
    sample_key = "samples_us_per_batch_call" if batch else "samples_us"
    return {key: 2.0, sample_key: [2.0] * iterations}


def base_parameters(cell: dict) -> dict:
    return {
        "K": cell["K"], "R": cell["R"],
        "shard_bytes": cell["shard_bytes"],
        "loss_count": cell["loss_count"],
        "missing_original_indices": list(range(cell["loss_count"])),
        "reuse": cell["reuse"], "iterations": cell["iterations"],
        "warmup": cell["warmup"], "seed": cell["seed"],
        "batch": 1, "thread_count": 1,
    }


def payload(codec: str, cell: dict) -> dict:
    parameters = base_parameters(cell)
    digest = {
        "algorithm": "fnv1a64", "original_data": "1" * 16,
        "recovered_originals": "2" * 16,
    }
    if codec == "leopard2":
        profile = "low_v1" if cell["K"] < 32 else "legacy_high_v1"
        parent, padded = ATLAS.expected_geometry(cell["K"], profile)
        parameters.update({
            "requested_profile": "auto", "requested_field": "auto",
            "requested_backend": "auto", "force_generic_decode": False,
            "force_specialized_decode": False, "force_tiled_decode": False,
            "force_materialized_decode": False, "skip_legacy": True,
            "retain_samples": True, "measure_one_shot_encode": True,
            "measure_one_shot_decode": True,
        })
        digest["transmitted_parity"] = "3" * 16
        return {
            "schema": "leopard2-benchmark-v9", "parameters": parameters,
            "correctness": {"leopard2_round_trip": True,
                            "legacy_comparison": None},
            "resolved": {
                "profile": profile, "field": "gf8", "backend": "avx2",
                "thread_count": 1, "parent_count": parent,
                "padded_side": padded,
            },
            "workload_digests": digest,
            "memory": {
                "encode_scratch_bytes_per_stripe": 128,
                "decode_scratch_bytes_per_stripe": 64,
                "encode_scratch_bytes_batch": 128,
                "decode_scratch_bytes_batch": 64,
            },
            "metrics": {
                "codec_setup": samples(cell["iterations"]),
                "encode_execution": samples(cell["iterations"], True),
                "one_shot_encode": samples(cell["iterations"], True),
                "decode_plan_setup": samples(cell["iterations"]),
                "decode_execution": samples(cell["iterations"], True),
                "one_shot_decode_including_setup": samples(
                    cell["iterations"], True),
                "rate_semantics": "offered_received counts all non-null "
                    "shard pointers supplied; a plan may read a deterministic "
                    "subset",
            },
        }
    if codec == "leopard1":
        parent, padded = ATLAS.expected_geometry(
            cell["K"], "legacy_high_v1")
        parameters["logical_shard_bytes"] = cell["shard_bytes"]
        digest["transmitted_parity"] = "3" * 16
        return {
            "schema": "leopard-main-benchmark-v1", "parameters": parameters,
            "build": {"main_source_commit": ATLAS.EXPECTED_MAIN_COMMIT,
                      "pure_avx2": True},
            "correctness": {"round_trip": True},
            "resolved": {
                "profile": "legacy_high_v1", "field": "gf8",
                "thread_count": 1, "parent_count": parent,
                "padded_side": padded, "padded_application_bytes": False,
                "padding_policy": "zero suffix per shard",
            },
            "workload_digests": digest,
            "memory": {"alignment": 64,
                       "encode_work_shards_per_stripe": 2 * padded,
                       "decode_work_shards_per_stripe": parent,
                       "encode_work_bytes_batch":
                           2 * padded * cell["shard_bytes"],
                       "decode_work_bytes_batch":
                           parent * cell["shard_bytes"]},
            "metrics": {
                "codec_setup": None,
                "decode_timing_includes_setup": True,
                "encode_execution": samples(cell["iterations"], True),
                "decode_including_setup": samples(cell["iterations"], True),
                "rate_semantics": "offered_received counts every non-null "
                    "original and all R supplied parity shards",
            },
        }
    if codec == "wirehair":
        digest["generated_repair"] = "4" * 16
        return {
            "schema": "leopard2-performance-atlas-wirehair-v1/v2",
            "parameters": parameters,
            "build": {"wirehair_source_commit": ATLAS.EXPECTED_WIREHAIR_COMMIT,
                      "pure_avx2": False, "wirehair_abi_version": 2,
                      "wire_profile_version": 1,
                      "wire_profile_id": ATLAS.EXPECTED_WIREHAIR_PROFILE_ID,
                      "isa_claim": "wirehair_v1_compact_path_avx2",
                      "target_qualified_avx512_helpers_present": True,
                      "wirehair_v1_wide_xor_forced_off": True,
                      "runtime_wide_xor_enabled": False,
                      "measured_path_avx512": False,
                      "active_x86_features": {
                          "ssse3": True, "avx2": True,
                          "gfni": False, "avx512": True}},
            "correctness": {"round_trip": True},
            "workload_digests": digest,
            "decode_input": {
                "surviving_systematic_blocks":
                    cell["K"] - cell["loss_count"],
                "repair_blocks_consumed": cell["loss_count"],
                "extra_repair_blocks": 0,
                "arrival_order":
                    "surviving_systematic_ascending_then_repair_ascending",
            },
            "path_semantics": {
                "codec": "shipping_wirehair_v1",
                "threading": "single_caller_thread",
                "wide_xor": "forced_off_on_benchmark_thread",
                "avx512_target_helpers":
                    "present_but_unreachable_from_measured_v1_compact_path",
            },
            "timing_semantics": {
                "message_precode_setup":
                    "fresh wirehair_encoder_create_ex; no repair emission",
                "encode_execution": "reuse one message-precode encoder; emit "
                    "exactly R repairs; exclude encoder creation",
                "encode_including_setup": "fresh wirehair_encoder_create_ex "
                    "then emit exactly R repairs",
                "decode_including_setup": "fresh decoder; ingest surviving "
                    "systematic blocks then ascending repair blocks through "
                    "solve; recover missing originals",
            },
            "metrics": {
                "message_precode_setup": samples(cell["iterations"]),
                "encode_execution": samples(cell["iterations"]),
                "encode_including_setup": samples(cell["iterations"]),
                "decode_including_setup": samples(cell["iterations"]),
            },
        }
    raise AssertionError(codec)


def external_build_fixture(root: Path, *, nested_archive: bool = False
                           ) -> tuple[dict, Path, Path, Path]:
    compiler = shutil.which("c++")
    archiver = shutil.which("ar")
    ranlib = shutil.which("ranlib")
    if not compiler or not archiver or not ranlib:
        raise unittest.SkipTest("C++/ar/ranlib tools are required")
    source_root = root / "source"
    build = root / "build"
    archive_root = build / "codec-subdir" if nested_archive else build
    archive_object = archive_root / "CMakeFiles/archive.dir/codec.cpp.o"
    executable_object = build / "CMakeFiles/bench.dir/adapter.cpp.o"
    archive_link = archive_root / "CMakeFiles/archive.dir/link.txt"
    executable_link = build / "CMakeFiles/bench.dir/link.txt"
    for directory in (source_root, archive_object.parent,
                      executable_object.parent):
        directory.mkdir(parents=True, exist_ok=True)
    codec_source = source_root / "codec.cpp"
    adapter_source = source_root / "adapter.cpp"
    codec_source.write_text("int codec_source = 1;\n", encoding="utf-8")
    adapter_source.write_text("int main() { return 0; }\n", encoding="utf-8")
    archive_object.write_bytes(b"compiled codec object\n")
    executable_object.write_bytes(b"compiled adapter object\n")
    archive = archive_root / "libcodec.a"
    subprocess.run([archiver, "qc", str(archive), str(archive_object)],
                   check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)
    subprocess.run([ranlib, str(archive)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    executable = build / "bench"
    executable.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    executable.chmod(0o755)
    flags = ["-march=x86-64", "-mtune=generic", "-mavx2",
             "-mno-avx512f"]
    rows = []
    for source, output in ((codec_source, archive_object),
                           (adapter_source, executable_object)):
        relative_output = output.relative_to(build).as_posix()
        command_directory = (archive_root if source == codec_source else
                             build)
        command_output = output.relative_to(command_directory).as_posix()
        command = [compiler, *flags, "-o", command_output,
                   "-c", str(source)]
        rows.append({
            "directory": str(command_directory),
            "command": " ".join(command),
            "file": str(source), "output": relative_output,
        })
    (build / "compile_commands.json").write_text(
        json.dumps(rows), encoding="utf-8")
    (build / "CMakeCache.txt").write_text(
        "CMAKE_BUILD_TYPE:STRING=Release\n"
        "PURE_MODE:BOOL=ON\n"
        f"CMAKE_CXX_COMPILER:FILEPATH={compiler}\n"
        f"CMAKE_AR:FILEPATH={archiver}\n"
        f"CMAKE_RANLIB:FILEPATH={ranlib}\n", encoding="utf-8")
    archive_link.write_text(
        f"{archiver} qc libcodec.a "
        f"{archive_object.relative_to(archive_root).as_posix()}\n"
        f"{ranlib} libcodec.a\n", encoding="utf-8")
    executable_link.write_text(
        f"{compiler} {executable_object.relative_to(build).as_posix()} "
        f"-o bench {archive.relative_to(build).as_posix()}\n",
        encoding="utf-8")
    arguments = {
        "build_root": build,
        "executable": executable,
        "archive_relative": archive.relative_to(build),
        "archive_link_relative": archive_link.relative_to(build),
        "archive_sources": {"codec.cpp": codec_source},
        "executable_sources": {"adapter.cpp": adapter_source},
        "cache_contract": {"PURE_MODE": "ON",
                           "CMAKE_BUILD_TYPE": "Release"},
        "archive_target_fragment": (
            "codec-subdir/CMakeFiles/archive.dir/" if nested_archive else
            "CMakeFiles/archive.dir/"),
        "executable_target_fragment": "CMakeFiles/bench.dir/",
        "label": "synthetic external",
    }
    return arguments, archive_object, archive, codec_source


class AtlasTests(unittest.TestCase):
    def test_grid_and_loss_contract(self) -> None:
        self.assertEqual(len(ATLAS.k_values()), 120)
        self.assertEqual(ATLAS.k_values()[-1], 224)
        self.assertEqual(ATLAS.conceptual_loss(1, "two"), None)
        self.assertEqual(ATLAS.conceptual_loss(101, "ten_percent"), 11)
        self.assertEqual(ATLAS.conceptual_loss(224, "full"), 32)
        self.assertEqual(ATLAS.seed_for(33, 64, 2),
                         ATLAS.seed_for(33, 1024 * 1024, 2))
        self.assertGreaterEqual(len(ATLAS.nice_log_ticks(2.52, 3.97)), 2)
        low_cells = [cell for cell in ATLAS.build_manifest()["cells"]
                     if cell["K"] < 32 and
                     len(cell["available_codecs"]) == 2]
        low_order = [ATLAS.codec_order(cell["id"])[0] for cell in low_cells]
        self.assertEqual(low_order.count("leopard2"),
                         low_order.count("wirehair"))

    def test_manifest_rejects_deletion_and_mutation(self) -> None:
        manifest = ATLAS.build_manifest()
        missing = copy.deepcopy(manifest)
        missing["cells"].pop()
        missing["matrix"]["unique_cell_count"] -= 1
        with self.assertRaises(ATLAS.AtlasError):
            ATLAS.validate_manifest(missing)
        changed = copy.deepcopy(manifest)
        changed["cells"][100]["available_codecs"] = ["leopard2"]
        with self.assertRaises(ATLAS.AtlasError):
            ATLAS.validate_manifest(changed)

    def test_payloads_bind_shape_samples_and_digests(self) -> None:
        cell = next(cell for cell in ATLAS.build_manifest()["cells"]
                    if cell["K"] == 32 and cell["shard_bytes"] == 64 and
                    cell["loss_count"] == 2)
        payloads = {codec: payload(codec, cell) for codec in ATLAS.CODECS}
        for codec, value in payloads.items():
            ATLAS.validate_payload(codec, value, cell)
        ATLAS.compare_cell_payloads(payloads, cell)
        broken = copy.deepcopy(payloads["wirehair"])
        broken["parameters"]["missing_original_indices"][-1] += 1
        with self.assertRaises(ATLAS.AtlasError):
            ATLAS.compare_cell_payloads(
                {**payloads, "wirehair": broken}, cell)
        broken = copy.deepcopy(payloads["leopard1"])
        broken["workload_digests"]["transmitted_parity"] = "f" * 16
        with self.assertRaises(ATLAS.AtlasError):
            ATLAS.compare_cell_payloads(
                {**payloads, "leopard1": broken}, cell)
        for codec, mutate in (
                ("leopard2", lambda value: value["parameters"].update(
                    {"thread_count": 2})),
                ("leopard2", lambda value: value["resolved"].update(
                    {"parent_count": 128})),
                ("leopard1", lambda value: value["metrics"].update(
                    {"decode_timing_includes_setup": False})),
                ("wirehair", lambda value: value["build"].update(
                    {"runtime_wide_xor_enabled": True})),
                ("wirehair", lambda value: value["decode_input"].update(
                    {"arrival_order": "random"}))):
            changed = copy.deepcopy(payloads[codec])
            mutate(changed)
            with self.assertRaises(ATLAS.AtlasError):
                ATLAS.validate_payload(codec, changed, cell)

    def test_short_samples_and_nonfinite_metrics_reject(self) -> None:
        cell = next(cell for cell in ATLAS.build_manifest()["cells"]
                    if cell["K"] == 32)
        candidate = payload("leopard2", cell)
        candidate["metrics"]["decode_execution"][
            "samples_us_per_batch_call"].pop()
        with self.assertRaises(ATLAS.AtlasError):
            ATLAS.validate_payload("leopard2", candidate, cell)
        candidate = payload("leopard2", cell)
        candidate["metrics"]["decode_execution"][
            "median_us_per_batch_call"] = 999999.0
        with self.assertRaises(ATLAS.AtlasError):
            ATLAS.validate_payload("leopard2", candidate, cell)

    def test_summary_binds_complete_cell_and_derived_rates(self) -> None:
        manifest = ATLAS.build_manifest()
        cell = next(cell for cell in manifest["cells"] if cell["K"] == 32)
        rows = []
        for codec in cell["available_codecs"]:
            value = payload(codec, cell)
            row = {
                "cell_id": cell["id"], "codec": codec,
                "K": cell["K"], "R": cell["R"],
                "shard_bytes": cell["shard_bytes"],
                "loss_count": cell["loss_count"],
                "loss_labels": cell["loss_labels"], "seed": cell["seed"],
            }
            row.update(ATLAS.normalize_result(codec, value, cell))
            rows.append(row)
        summary = {
            "schema": ATLAS.SUMMARY_SCHEMA,
            "manifest_sha256": ATLAS.sha256_bytes(
                ATLAS.canonical_json(manifest).encode()),
            "run_metadata_sha256": "0" * 64,
            "complete": False, "complete_cell_count": 1,
            "expected_cell_count": len(manifest["cells"]), "rows": rows,
        }
        ATLAS.validate_summary(summary, manifest, require_complete=False)
        partial = copy.deepcopy(summary)
        partial["rows"].pop()
        with self.assertRaises(ATLAS.AtlasError):
            ATLAS.validate_summary(partial, manifest, require_complete=False)
        forged = copy.deepcopy(summary)
        forged["rows"][0]["encode_execution_message_GBps"] *= 2
        with self.assertRaises(ATLAS.AtlasError):
            ATLAS.validate_summary(forged, manifest, require_complete=False)
        candidate = payload("leopard2", cell)
        candidate["metrics"]["decode_execution"][
            "median_us_per_batch_call"] = float("nan")
        with self.assertRaises(ATLAS.AtlasError):
            ATLAS.validate_payload("leopard2", candidate, cell)

    def test_raw_bundle_binds_metadata_and_exact_record_set(self) -> None:
        manifest = ATLAS.build_manifest()
        records = []
        for cell in manifest["cells"]:
            for codec in cell["available_codecs"]:
                records.append({"cell_id": cell["id"], "codec": codec,
                                "payload": payload(codec, cell)})
        bundle = {
            "schema": ATLAS.RAW_BUNDLE_SCHEMA,
            "manifest_sha256": ATLAS.sha256_bytes(
                ATLAS.canonical_json(manifest).encode()),
            "run_metadata_sha256": "a" * 64,
            "record_count": len(records), "records": records,
        }
        ATLAS.validate_raw_bundle(bundle, manifest, "a" * 64)
        missing = copy.deepcopy(bundle)
        missing["records"].pop()
        missing["record_count"] -= 1
        with self.assertRaises(ATLAS.AtlasError):
            ATLAS.validate_raw_bundle(missing, manifest, "a" * 64)
        with self.assertRaises(ATLAS.AtlasError):
            ATLAS.validate_raw_bundle(bundle, manifest, "b" * 64)

    def test_external_build_closure_binds_sources_objects_and_archive(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            arguments, _, _, _ = external_build_fixture(Path(directory))
            closure = ATLAS.capture_external_avx2_compile_build(**arguments)
            graph = closure["source_object_graph"]
            self.assertEqual(set(graph["archive"]), {"codec.cpp"})
            self.assertEqual(set(graph["executable"]), {"adapter.cpp"})
            self.assertEqual(
                closure["archive_members"][0]["sha256"],
                graph["archive"]["codec.cpp"]["object"]["sha256"])

        with tempfile.TemporaryDirectory() as directory:
            arguments, _, _, _ = external_build_fixture(
                Path(directory), nested_archive=True)
            closure = ATLAS.capture_external_avx2_compile_build(**arguments)
            self.assertIn("codec-subdir/CMakeFiles/archive.dir/codec.cpp.o",
                          closure["source_object_graph"]["archive"]
                          ["codec.cpp"]["object"]["path"])

        with tempfile.TemporaryDirectory() as directory:
            arguments, _, _, source = external_build_fixture(Path(directory))
            wrong_source = source.with_name("wrong.cpp")
            wrong_source.write_text(source.read_text(encoding="utf-8"),
                                    encoding="utf-8")
            arguments["archive_sources"] = {"codec.cpp": wrong_source}
            with self.assertRaisesRegex(
                    ATLAS.AtlasError, "expected the declared source"):
                ATLAS.capture_external_avx2_compile_build(**arguments)

        with tempfile.TemporaryDirectory() as directory:
            arguments, object_path, _, _ = external_build_fixture(
                Path(directory))
            with object_path.open("ab") as output:
                output.write(b"object changed after archive\n")
            with self.assertRaisesRegex(
                    ATLAS.AtlasError, "differs from its compiled object"):
                ATLAS.capture_external_avx2_compile_build(**arguments)

        with tempfile.TemporaryDirectory() as directory:
            arguments, _, archive, _ = external_build_fixture(Path(directory))
            extra = Path(directory) / "extra.cpp.o"
            extra.write_bytes(b"unexpected archive object\n")
            subprocess.run([shutil.which("ar"), "q", str(archive), str(extra)],
                           check=True, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL)
            with self.assertRaisesRegex(ATLAS.AtlasError,
                                        "archive member set differs"):
                ATLAS.capture_external_avx2_compile_build(**arguments)

    def test_plotter_emits_all_expected_dependency_free_svgs(self) -> None:
        rows = []
        for block_bytes in ATLAS.BLOCK_BYTES:
            for k in (1, 2, 32, 224):
                for loss_label in ATLAS.LOSS_LABELS:
                    losses = ATLAS.conceptual_loss(k, loss_label)
                    if losses is None:
                        continue
                    for codec in ATLAS.available_codecs(k):
                        base = 1.0 + k / 100.0 + block_bytes / 1048576.0
                        scale = {"leopard2": 1.0, "leopard1": 1.2,
                                 "wirehair": 1.5}[codec]
                        usec = base * scale * (1 + losses / 64.0)
                        rows.append({
                            "cell_id": f"synthetic-{block_bytes}-{k}-{loss_label}",
                            "codec": codec, "K": k, "R": 32,
                            "shard_bytes": block_bytes,
                            "loss_count": losses, "loss_labels": [loss_label],
                            "encode_execution_us": usec,
                            "decode_first_us": usec * 1.3,
                            "decode_reused_us": (None if codec == "wirehair"
                                                  else usec),
                            "encode_execution_message_GBps": k * block_bytes /
                                (usec * 1000.0),
                            "encode_execution_output_GBps": 32 * block_bytes /
                                (usec * 1000.0),
                            "wirehair_repair_emission_us": (usec / 3
                                if codec == "wirehair" else None),
                            "wirehair_repair_emission_output_GBps": (
                                32 * block_bytes / (usec / 3 * 1000.0)
                                if codec == "wirehair" else None),
                            "decode_first_message_GBps": k * block_bytes /
                                (usec * 1.3 * 1000.0),
                            "decode_first_repaired_GBps": losses * block_bytes /
                                (usec * 1.3 * 1000.0),
                            "decode_first_received_GBps": (k - losses + 32) *
                                block_bytes / (usec * 1.3 * 1000.0),
                            "decode_reused_message_GBps": (None if codec == "wirehair"
                                else k * block_bytes / (usec * 1000.0)),
                            "extra_repair_blocks": 0,
                        })
        summary = {"rows": rows}
        with tempfile.TemporaryDirectory() as directory:
            paths = ATLAS.generate_plots(summary, Path(directory))
            self.assertEqual(len(paths), 39)
            for path in paths:
                text = path.read_text(encoding="utf-8")
                self.assertTrue(text.startswith("<svg"))
                self.assertIn("K (log", text)

    def test_command_contract(self) -> None:
        cell = ATLAS.build_manifest()["cells"][0]
        command = ATLAS.command_for("leopard2", Path("/candidate"), cell)
        self.assertIn("--measure-one-shot-decode", command)
        self.assertIn("--skip-legacy", command)
        self.assertEqual(command[-1], "--measure-one-shot-decode")
        main = ATLAS.command_for("leopard1", Path("/main"), {
            **cell, "K": 32})
        self.assertNotIn("--skip-legacy", main)


if __name__ == "__main__":
    unittest.main()
