#!/usr/bin/env python3
"""Independent tests for tools/leopard2_operation_counts.py."""

from __future__ import annotations

import csv
import importlib.util
import io
import json
import pathlib
import subprocess
import sys
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[3]
TOOL = ROOT / "tools" / "leopard2_operation_counts.py"
SPEC = importlib.util.spec_from_file_location("leopard2_operation_counts", TOOL)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load operation-count tool")
COUNTS = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = COUNTS
SPEC.loader.exec_module(COUNTS)


def independent_full_transform_count(size: int) -> int:
    operations = 0
    span = 2
    while span <= size:
        for _block in range(0, size, span):
            for _lane in range(span // 2):
                operations += 1
        span *= 2
    return operations


class OperationCountTests(unittest.TestCase):
    def test_full_schedules_against_explicit_layer_enumeration(self) -> None:
        for size in (1, 2, 4, 8, 16, 32, 64, 128, 256):
            expected = independent_full_transform_count(size)
            self.assertEqual(COUNTS.ifft_prefix_butterflies(size, size), expected)
            self.assertEqual(COUNTS.fft_prefix_butterflies(size, size), expected)
            self.assertEqual(
                COUNTS.fft_mask_butterflies(size, set(range(size))), expected
            )

    def test_mask_schedule_is_prefix_schedule_for_prefix_masks(self) -> None:
        for size in (2, 4, 8, 16, 32, 64):
            for prefix in range(size + 1):
                self.assertEqual(
                    COUNTS.fft_mask_butterflies(size, set(range(prefix))),
                    COUNTS.fft_prefix_butterflies(size, prefix),
                )

    def test_high_and_low_representative_cells(self) -> None:
        high = COUNTS.model_high_encode(240, 16, 16, 1024, set(range(16)))
        self.assertEqual(high.butterflies, 512)
        self.assertEqual(high.nontransform_xor_vectors, 224)
        low = COUNTS.model_low_encode(8, 248, 8, 1024, set(range(248)))
        self.assertEqual(low.butterflies, 384)
        self.assertEqual(low.details["active_parity_blocks"], 31)
        self.assertEqual(low.copies, 8 + 248)
        self.assertTrue(low.details["out_of_place_first_fft_layer"])

        high_decode = COUNTS.model_high_decode(
            240, 16, 256, 16, 1024, {0}
        )
        self.assertEqual(
            high_decode.details["receive_source_fused_radix4_groups"], 56
        )
        self.assertEqual(high_decode.details["receive_copy_vectors"], 16)
        self.assertEqual(high_decode.details["receive_zero_vectors"], 16)
        self.assertEqual(
            high_decode.details["receive_copy_vectors_removed"], 224
        )
        self.assertEqual(
            high_decode.details["receive_exact_pruned_staged_blocks"], 2
        )
        # These are exact deltas at the new receive boundary.  The broader
        # model's pre-existing absolute copy total is tracked separately and
        # is not an oracle for aligned coordinate-pointer retention.
        high_decode_gf16 = COUNTS.model_high_decode(
            240, 16, 256, 16, 1024, {0}, "gf16"
        )
        self.assertEqual(
            high_decode_gf16.details["receive_source_fused_radix4_groups"], 0
        )
        self.assertEqual(high_decode_gf16.details["receive_copy_vectors"], 240)
        self.assertEqual(
            high_decode_gf16.details["receive_copy_vectors_removed"], 0
        )
        # Selected input coordinates are retained as pointers on aligned
        # passes; they are not K shard copies.  High Algorithm 5 separately
        # stages K rows into its transform workspace.
        self.assertEqual(high_decode.details["coordinate_pointer_mappings"], 240)
        self.assertEqual(high_decode.details["aligned_input_staging_copies"], 0)
        self.assertEqual(high_decode.copies, 240)
        self.assertEqual(high_decode.decode_output_gather_payload_bytes, 1024)

        low_decode = COUNTS.model_low_decode(8, 248, 256, 8, 65, {0, 1})
        self.assertEqual(low_decode.copies, 0)
        self.assertEqual(low_decode.decode_output_gather_payload_bytes, 2)
        self.assertEqual(low_decode.decode_coordinate_pointer_mappings, 8)

    def test_low_encode_source_guard_rejects_old_copy_loop(self) -> None:
        old_copy = """
            for (unsigned i = 0; i < p; ++i)
                memcpy(work[p + i], work[i], buffer_bytes);
        """
        for filename in ("LeopardFF8.cpp", "LeopardFF16.cpp"):
            source = (ROOT / filename).read_text(encoding="utf-8")
            COUNTS.verify_low_encode_no_copy_source(source, filename)
            with self.assertRaises(COUNTS.ModelError):
                COUNTS.verify_low_encode_no_copy_source(
                    source + old_copy, filename + " mutation"
                )

    def test_high_decode_source_boundary_guard(self) -> None:
        filename = "LeopardFF8.cpp"
        source = (ROOT / filename).read_text(encoding="utf-8")
        COUNTS.verify_high_decode_receive_fusion_source(source, filename)
        decoder_call = (
            "ops.ff8_ifft_butterfly4_out(\n"
            "                sources[r], sources[r + 1], sources[r + 2], "
            "sources[r + 3],"
        )
        mutated_decoder = source.replace(
            decoder_call,
            decoder_call.replace(
                "ops.ff8_ifft_butterfly4_out(",
                "ops.ff8_ifft_butterfly4_out_removed(",
            ),
            1,
        )
        self.assertNotEqual(mutated_decoder, source)
        with self.assertRaises(COUNTS.ModelError):
            COUNTS.verify_high_decode_receive_fusion_source(
                mutated_decoder, filename + " decoder-call mutation"
            )
        mutated_sink = source.replace(
            "ExecutePrunedInverseTransformPlanAccumulate(",
            "ExecutePrunedInverseTransformPlanAccumulateRemoved(",
            1,
        )
        self.assertNotEqual(mutated_sink, source)
        with self.assertRaises(COUNTS.ModelError):
            COUNTS.verify_high_decode_receive_fusion_source(
                mutated_sink, filename + " pruned-sink mutation"
            )
        mutated = source.replace(
            "IFFT_DIT_DecoderFromSources(",
            "IFFT_DIT_DecoderCopyFirst(",
        )
        with self.assertRaises(COUNTS.ModelError):
            COUNTS.verify_high_decode_receive_fusion_source(
                mutated, filename + " mutation"
            )

        filename = "LeopardFF16.cpp"
        source = (ROOT / filename).read_text(encoding="utf-8")
        COUNTS.verify_high_decode_gf16_copy_first_source(source, filename)
        mutated = source.replace(
            "StageHighDecodeSources(buffer_bytes, coordinate_data, work, n);",
            "StageHighDecodeSources(buffer_bytes, coordinate_data, work, t);",
        )
        self.assertNotEqual(mutated, source)
        with self.assertRaises(COUNTS.ModelError):
            COUNTS.verify_high_decode_gf16_copy_first_source(
                mutated, filename + " per-block mutation"
            )

    def test_high_decode_integration_keeps_both_counter_families(self) -> None:
        tokens = (
            "syndrome_pruned_accumulated_blocks",
            "syndrome_pruned_fallback_blocks",
            "receive_ifft_butterfly4_out_of_place",
            "receive_copy_shards",
            "receive_zero_shards",
        )
        for stem in ("LeopardFF8", "LeopardFF16"):
            source = (ROOT / (stem + ".cpp")).read_text(encoding="utf-8")
            header = (ROOT / (stem + ".h")).read_text(encoding="utf-8")
            for token in tokens:
                self.assertIn(token, header, stem + " header lost " + token)
            self.assertIn(
                "TestHighSyndromePrunedAccumulatedBlocks", source,
                stem + " source lost exact-pruned accounting",
            )
            self.assertIn(
                "TestHighReceiveCopyShards", source,
                stem + " source lost receive-stage accounting",
            )

    def test_high_pruned_stage_hook_target_cannot_compile_out(self) -> None:
        filename = "CMakeLists.txt"
        source = (ROOT / filename).read_text(encoding="utf-8")
        COUNTS.verify_high_pruned_hook_registration(source, filename)
        mutated = source.replace(
            "        LEO2_ENABLE_TEST_HOOKS=1\n"
            "        LEO2_REQUIRE_HIGH_PRUNED_STAGE_HOOKS=1)",
            "        LEO2_REQUIRE_HIGH_PRUNED_STAGE_HOOKS=1)",
            1,
        )
        self.assertNotEqual(mutated, source)
        with self.assertRaises(COUNTS.ModelError):
            COUNTS.verify_high_pruned_hook_registration(
                mutated, filename + " hook-definition mutation"
            )

    def test_sparse_requested_parity_reduces_low_work(self) -> None:
        all_outputs = COUNTS.model_low_encode(8, 248, 8, 1024, set(range(248)))
        sparse = COUNTS.model_low_encode(8, 248, 8, 1024, {0, 247})
        self.assertLess(sparse.butterflies, all_outputs.butterflies)
        self.assertEqual(sparse.details["active_parity_blocks"], 2)
        self.assertEqual(sparse.copies, 8 + 2)

    def test_deterministic_received_subset(self) -> None:
        high_data, parity = COUNTS.deterministic_decode_coordinates(
            9, 7, "high", 8, {0, 3}
        )
        self.assertEqual(parity, {0, 1})
        self.assertEqual(len(high_data), 9)
        self.assertIn(8 + 1, high_data)
        self.assertNotIn(8 + 3, high_data)

    def test_direct_bound_and_dispatch_gates(self) -> None:
        model = COUNTS.model_direct_repair(16, 8, 16, 16, "low", set(range(4)))
        self.assertEqual(model.fixed_multiply_vectors, 64)
        self.assertEqual(model.nontransform_xor_vectors, 60)
        self.assertEqual(model.api_scratch_slots, 0)
        with self.assertRaises(COUNTS.ModelError):
            COUNTS.model_direct_repair(17, 8, 32, 32, "low", {0})

    def test_decode_scratch_exact_components(self) -> None:
        tiled = COUNTS.decode_scratch_accounting(
            240, 16, 256, 16, "high", 65, 1, "tiled", 8,
            codec_workspace="tiled",
        )
        self.assertEqual(tiled.range_count, 496)
        self.assertEqual(tiled.range_metadata_bytes, 7936)
        self.assertEqual(tiled.plan_pointer_count, 256 + 2 * 16 + 1)
        self.assertEqual(tiled.tail_reserved_slots, 256)
        self.assertEqual(tiled.tail_selected_staged_slots, 240)
        self.assertEqual(tiled.tail_staged_payload_bytes, 240)
        self.assertEqual(tiled.tail_staged_zero_padding_bytes, 240 * 63)
        self.assertEqual(tiled.work_slot_bytes, 64)
        self.assertEqual(tiled.plan_work_slots, 33)
        self.assertEqual(tiled.codec_work_slots, 48)
        self.assertEqual(tiled.plan_total_bytes, 28800)
        self.assertEqual(tiled.codec_total_bytes, 29824)

        materialized = COUNTS.decode_scratch_accounting(
            240, 16, 256, 16, "high", 128, 1, "materialized", 8,
        )
        self.assertEqual(materialized.tail_reserved_slots, 0)
        self.assertEqual(materialized.plan_work_slots, 256)
        self.assertEqual(materialized.work_slot_bytes, 128)
        direct = COUNTS.decode_scratch_accounting(
            16, 8, 32, 16, "low", 65, 4, "specialized", 8,
            codec_workspace="specialized", direct=True,
        )
        self.assertEqual(direct.plan_work_slots, 0)
        self.assertEqual(direct.plan_pointer_count, 0)
        self.assertEqual(direct.plan_total_bytes, 640)
        self.assertGreater(direct.codec_total_bytes, direct.plan_total_bytes)
        no_op = COUNTS.decode_scratch_accounting(
            240, 16, 256, 16, "high", 65, 0, "tiled", 8,
            codec_workspace="tiled", no_op=True,
        )
        self.assertEqual(no_op.plan_total_bytes, 0)
        self.assertGreater(no_op.codec_total_bytes, 0)

    def test_decode_scratch_source_guard_rejects_layout_mutations(self) -> None:
        filename = "leopard2.cpp"
        source = (ROOT / filename).read_text(encoding="utf-8")
        COUNTS.verify_decode_scratch_source(source, filename)
        mutations = (
            (
                "static_cast<size_t>(codec->original_count) * 2 + "
                "codec->recovery_count",
                "static_cast<size_t>(codec->original_count) + "
                "codec->recovery_count",
            ),
            (
                "ComputeScratchLayout(range_count, 0, 0, rounded_bytes, layout)",
                "ComputeScratchLayout(range_count, 1, 0, rounded_bytes, layout)",
            ),
        )
        for old, new in mutations:
            mutated = source.replace(old, new, 1)
            self.assertNotEqual(mutated, source)
            with self.assertRaises(COUNTS.ModelError):
                COUNTS.verify_decode_scratch_source(mutated, filename + " mutation")

    def test_json_and_csv_are_deterministic(self) -> None:
        command = [
            sys.executable, str(TOOL), "report",
            "--path", "generic_decode", "--profile", "high",
            "--k", "9", "--r", "7", "--field", "gf8",
            "--shard-bytes", "65", "--loss-mask", "0,3",
        ]
        first = subprocess.check_output(command, text=True)
        second = subprocess.check_output(command, text=True)
        self.assertEqual(first, second)
        document = json.loads(first)
        self.assertEqual(document["schema_version"], 2)
        self.assertEqual(
            document["metrics"]["decode_plan_scratch_bytes"]["classification"],
            "exact_schedule",
        )
        self.assertEqual(
            document["metrics"]["transform_butterflies"]["classification"],
            "exact_schedule",
        )
        csv_text = subprocess.check_output(command + ["--format", "csv"], text=True)
        rows = list(csv.DictReader(io.StringIO(csv_text)))
        self.assertGreater(len(rows), 10)
        self.assertEqual([row["metric"] for row in rows],
                         sorted(row["metric"] for row in rows))

    def test_gf16_odd_length_is_rejected(self) -> None:
        completed = subprocess.run(
            [sys.executable, str(TOOL), "report", "--path", "low_encode",
             "--k", "8", "--r", "8", "--field", "gf16",
             "--shard-bytes", "3"],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            check=False,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("even shard byte count", completed.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
