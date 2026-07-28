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
import types
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


def build_decode_report(**overrides: object) -> dict:
    """Build one report with explicit stable defaults for selector tests."""
    values = {
        "path": "legacy_high_decode",
        "k": 240,
        "r": 16,
        "profile": None,
        "field": "gf8",
        "backend": "scalar",
        "shard_bytes": 4096,
        "requested_parity": None,
        "loss_count": 1,
        "loss_mask": None,
        "direction": "forward",
        "active_input_count": None,
        "transform_output_mask": None,
        "pointer_bytes": 8,
        "decode_workspace": "materialized",
        "decode_selection": "auto",
        "multi_item_batch": False,
        "context_auto_requested": False,
        "gf16_detected_l3_bytes": 0,
    }
    values.update(overrides)
    return COUNTS.build_report(types.SimpleNamespace(**values))


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
        self.assertEqual(low.copies, 8)
        self.assertEqual(low.details["direct_output_blocks"], 31)
        self.assertTrue(low.details["out_of_place_first_fft_layer"])

        high_decode = COUNTS.model_high_decode(
            240, 16, 256, 16, 1024, {0}
        )
        self.assertEqual(
            high_decode.details["receive_source_fused_radix4_groups"], 56
        )
        self.assertEqual(high_decode.details["receive_copy_vectors"], 15)
        self.assertEqual(high_decode.details["receive_zero_vectors"], 1)
        self.assertEqual(
            high_decode.details["receive_copy_vectors_removed"], 224
        )
        self.assertEqual(
            high_decode.details["receive_exact_pruned_staged_blocks"], 1
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
        self.assertEqual(high_decode_gf16.details["receive_copy_vectors"], 239)
        self.assertEqual(
            high_decode_gf16.details["receive_copy_vectors_removed"], 0
        )
        # Selected input coordinates are retained as pointers on aligned
        # passes; they are not K shard copies.  High Algorithm 5 separately
        # stages K rows into its transform workspace.
        self.assertEqual(high_decode.details["coordinate_pointer_mappings"], 240)
        self.assertEqual(high_decode.details["aligned_input_staging_copies"], 0)
        self.assertEqual(high_decode.copies, 239)
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
        mutated_backend = source.replace(
            "        ops.kind == LEO2_BACKEND_AVX2 ||\n"
            "        ops.kind == LEO2_BACKEND_AVX512 ||\n"
            "        ops.kind == LEO2_BACKEND_GFNI;",
            "        ops.kind == LEO2_BACKEND_AVX2 ||\n"
            "        ops.kind == LEO2_BACKEND_GFNI;",
            1,
        )
        self.assertNotEqual(mutated_backend, source)
        with self.assertRaises(COUNTS.ModelError):
            COUNTS.verify_high_decode_receive_fusion_source(
                mutated_backend, filename + " backend mutation"
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
            "buffer_bytes, coordinate_data + offset, block_work, t);",
            "buffer_bytes, coordinate_data + offset, block_work, n);",
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
            "syndrome_block_zero_ifft_elisions",
            "syndrome_forward_transforms",
            "syndrome_forward_transform_elisions",
            "syndrome_block_zero_xor_shards",
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

    def test_decode_fusion_source_guards_reject_policy_mutations(self) -> None:
        core = (ROOT / "leopard2.cpp").read_text(encoding="utf-8")
        ff8 = (ROOT / "LeopardFF8.cpp").read_text(encoding="utf-8")
        ff16 = (ROOT / "LeopardFF16.cpp").read_text(encoding="utf-8")
        COUNTS.verify_decode_fusion_sources(core, ff8, ff16)
        mutations = (
            (
                core.replace("aligned_prefix_bytes >= 4096",
                             "aligned_prefix_bytes >= 2048", 1),
                ff8,
                ff16,
                "generic reveal crossover",
            ),
            (
                core,
                ff8.replace("buffer_bytes >= 65536",
                            "buffer_bytes >= 32768", 1),
                ff16,
                "GF8 pruned sink crossover",
            ),
            (
                core,
                ff8,
                ff16.replace("bytes == 64", "bytes == 32", 1),
                "GF16 compact fused-four policy",
            ),
            (
                core.replace(
                    "         codec->context->backend == "
                    "LEO2_BACKEND_AVX2 ||\n"
                    "         codec->context->backend == "
                    "LEO2_BACKEND_AVX512",
                    "         codec->context->backend == "
                    "LEO2_BACKEND_AVX2",
                    1,
                ),
                ff8,
                ff16,
                "generic reveal AVX512 qualification",
            ),
            (
                core,
                ff8.replace(
                    "    if (ops.kind == LEO2_BACKEND_AVX2 ||\n"
                    "        ops.kind == LEO2_BACKEND_AVX512 ||\n"
                    "        ops.kind == LEO2_BACKEND_GFNI)\n"
                    "        return buffer_bytes >= 1024;",
                    "    if (ops.kind == LEO2_BACKEND_AVX2 ||\n"
                    "        ops.kind == LEO2_BACKEND_GFNI)\n"
                    "        return buffer_bytes >= 1024;",
                    1,
                ),
                ff16,
                "GF8 pruned sink AVX512 qualification",
            ),
            (
                core,
                ff8,
                ff16.replace(
                    "         (ops.kind == LEO2_BACKEND_AVX2 ||\n"
                    "          ops.kind == LEO2_BACKEND_AVX512 ||\n"
                    "          ops.kind == LEO2_BACKEND_GFNI));",
                    "         (ops.kind == LEO2_BACKEND_AVX2 ||\n"
                    "          ops.kind == LEO2_BACKEND_GFNI));",
                    1,
                ),
                "GF16 compact AVX512 qualification",
            ),
            (
                core.replace(
                    "const bool reveal_aligned_outputs_in_place =\n"
                    "        !(fuse_generic_reveal_scatter || "
                    "fuse_low_reveal_scatter);",
                    "const bool reveal_aligned_outputs_in_place = true;",
                    1,
                ),
                ff8,
                ff16,
                "aligned reveal execution wiring",
            ),
            (
                core.replace(
                    "ExecuteTransformDecodePass(\n"
                    "            plan, kScratchAlignment, coordinate_input, "
                    "work,\n"
                    "            use_generic, use_tiled, true, NULL);",
                    "ExecuteTransformDecodePass(\n"
                    "            plan, kScratchAlignment, coordinate_input, "
                    "work,\n"
                    "            use_generic, use_tiled, "
                    "reveal_aligned_outputs_in_place, NULL);",
                    1,
                ),
                ff8,
                ff16,
                "tail reveal execution wiring",
            ),
            (
                core,
                ff8.replace(
                    "const bool use_accumulating_sink =\n"
                    "                ShouldUsePrunedHighSyndromeSink("
                    "ops, buffer_bytes);",
                    "const bool use_accumulating_sink = false;",
                    1,
                ),
                ff16,
                "GF8 pruned sink wiring",
            ),
            (
                core,
                ff8,
                ff16.replace(
                    "const bool use_accumulating_sink =\n"
                    "                ShouldUsePrunedHighSyndromeSink("
                    "ops, buffer_bytes);",
                    "const bool use_accumulating_sink = false;",
                    1,
                ),
                "GF16 pruned sink wiring",
            ),
        )
        for mutated_core, mutated_ff8, mutated_ff16, label in mutations:
            with self.subTest(label=label):
                self.assertTrue(
                    mutated_core != core or mutated_ff8 != ff8 or
                    mutated_ff16 != ff16,
                    label + " mutation did not apply",
                )
                with self.assertRaises(COUNTS.ModelError):
                    COUNTS.verify_decode_fusion_sources(
                        mutated_core, mutated_ff8, mutated_ff16
                    )

    def test_decode_selection_matches_promoted_boundaries(self) -> None:
        for backend in ("scalar", "ssse3", "avx2", "avx512"):
            with self.subTest(balanced_backend=backend):
                balanced = build_decode_report(
                    k=128, r=128, loss_count=128, shard_bytes=256,
                    backend=backend,
                )
                self.assertEqual(
                    balanced["selection"],
                    {
                        "path": "materialized",
                        "rule": "translated_low",
                        "matching_auto_rules": 0,
                        "required_work_slots": 256,
                    },
                )
        below_balanced = build_decode_report(
            k=128, r=128, loss_count=128, shard_bytes=192,
        )
        self.assertEqual(
            (below_balanced["selection"]["path"],
             below_balanced["selection"]["rule"]),
            ("materialized", "translated_low"),
        )
        ragged_balanced = build_decode_report(
            k=128, r=128, loss_count=128, shard_bytes=193,
        )
        self.assertEqual(ragged_balanced["selection"]["rule"],
                         "translated_low")
        above_balanced = build_decode_report(
            k=128, r=128, loss_count=128, shard_bytes=1024 * 1024 + 1,
        )
        self.assertEqual(above_balanced["selection"]["rule"],
                         "translated_low")
        forced_specialized = build_decode_report(
            k=128, r=128, loss_count=128, shard_bytes=256,
            decode_selection="specialized",
        )
        self.assertEqual(
            (forced_specialized["selection"]["path"],
             forced_specialized["selection"]["rule"]),
            ("materialized", "translated_low"),
        )
        self.assertEqual(
            forced_specialized["selection"]["matching_auto_rules"], 0
        )
        forced_generic = build_decode_report(
            k=128, r=128, loss_count=128, shard_bytes=256,
            decode_selection="generic",
        )
        self.assertEqual(
            forced_generic["selection"],
            {
                "path": "generic",
                "rule": "forced_generic",
                "matching_auto_rules": 1,
                "required_work_slots": 256,
            },
        )

        materialized = build_decode_report(
            k=224, r=32, loss_count=8, shard_bytes=32 * 1024,
            backend="avx2",
        )
        self.assertEqual(
            (materialized["selection"]["path"],
             materialized["selection"]["rule"]),
            ("materialized", "measured_materialized"),
        )
        batch = build_decode_report(
            k=224, r=32, loss_count=8, shard_bytes=32 * 1024,
            backend="avx2", multi_item_batch=True,
        )
        self.assertEqual(
            (batch["selection"]["path"], batch["selection"]["rule"]),
            ("tiled", "measured_batch_tiled"),
        )
        self.assertEqual(
            batch["metrics"]["execution_working_slots"]["value"], 72
        )
        self.assertEqual(
            batch["metrics"]["decode_plan_work_slots"]["value"], 256
        )
        avx512_batch = build_decode_report(
            k=224, r=32, loss_count=8, shard_bytes=32 * 1024,
            backend="avx512", multi_item_batch=True,
        )
        self.assertEqual(
            (avx512_batch["selection"]["path"],
             avx512_batch["selection"]["rule"]),
            ("tiled", "measured_batch_tiled"),
        )
        below_materialized = build_decode_report(
            k=224, r=32, loss_count=8, shard_bytes=16 * 1024,
            backend="avx2",
        )
        self.assertEqual(
            (below_materialized["selection"]["path"],
             below_materialized["selection"]["rule"]),
            ("tiled", "workspace_tiled"),
        )
        for backend, threshold in (("avx2", 24 * 1024),
                                   ("avx512", 24 * 1024),
                                   ("ssse3", 32 * 1024)):
            with self.subTest(backend=backend):
                below = build_decode_report(
                    k=224, r=32, loss_count=1,
                    shard_bytes=threshold - 64, backend=backend,
                )
                ragged_boundary = build_decode_report(
                    k=224, r=32, loss_count=1,
                    shard_bytes=threshold - 63, backend=backend,
                )
                self.assertEqual(below["selection"]["rule"],
                                 "workspace_tiled")
                self.assertEqual(ragged_boundary["selection"]["rule"],
                                 "measured_materialized")
        above_materialized = build_decode_report(
            k=224, r=32, loss_count=8, shard_bytes=64 * 1024 + 1,
            backend="avx2",
        )
        self.assertEqual(above_materialized["selection"]["rule"],
                         "workspace_tiled")
        too_many_losses = build_decode_report(
            k=224, r=32, loss_count=9, shard_bytes=32 * 1024,
            backend="avx2",
        )
        self.assertEqual(too_many_losses["selection"]["rule"],
                         "workspace_tiled")

        no_op = build_decode_report(loss_count=0)
        self.assertEqual(no_op["selection"]["path"], "no_op")
        geometry = COUNTS.parent_geometry(240, 16, "high")
        for unvalidated_bytes in (0, COUNTS.UINT64_MAX):
            no_op_selection = COUNTS.select_decode_execution(
                240, 16, "high", "gf8", "scalar", geometry,
                unvalidated_bytes, 0,
            )
            self.assertEqual(no_op_selection.path, "no_op")
        with self.assertRaises(COUNTS.ModelError):
            build_decode_report(loss_count=0, shard_bytes=COUNTS.UINT64_MAX)
        xor_terminal = build_decode_report(k=9, r=1, loss_count=1)
        self.assertEqual(xor_terminal["selection"]["path"], "direct")
        self.assertEqual(xor_terminal["details"]["direct_mode"], "r1_xor")
        self.assertEqual(
            xor_terminal["metrics"]["nontransform_xor_vectors"]["value"], 8
        )
        xor_terminal_ragged = build_decode_report(
            k=9, r=1, loss_count=1, shard_bytes=65
        )
        self.assertEqual(
            xor_terminal_ragged["metrics"]["kernel_shard_bytes"]["value"], 65
        )
        copy_terminal = build_decode_report(
            path="low_decode", k=1, r=8, loss_count=1,
        )
        self.assertEqual(copy_terminal["selection"]["path"], "direct")
        self.assertEqual(copy_terminal["details"]["direct_mode"], "k1_copy")
        self.assertEqual(
            copy_terminal["metrics"]["logical_copy_vectors"]["value"], 1
        )
        direct = build_decode_report(
            path="low_decode", profile=None, k=16, r=8, loss_count=4,
        )
        self.assertEqual(direct["selection"]["path"], "direct")
        direct_disabled = build_decode_report(
            path="low_decode", profile=None, k=16, r=8, loss_count=4,
            decode_selection="specialized",
        )
        self.assertNotEqual(direct_disabled["selection"]["path"], "direct")

    def test_backend_qualified_reveal_traffic(self) -> None:
        fused = build_decode_report(
            path="generic_decode", profile="high", k=240, r=16,
            loss_count=2, shard_bytes=4096, backend="avx2",
            decode_selection="path",
        )["metrics"]
        self.assertEqual(
            fused["decode_reveal_direct_payload_bytes"]["value"], 8192
        )
        self.assertEqual(
            fused["decode_reveal_inplace_temporary_payload_bytes"]["value"], 0
        )
        self.assertEqual(
            fused["decode_reveal_scatter_payload_bytes"]["value"], 0
        )
        avx512_fused = build_decode_report(
            path="generic_decode", profile="high", k=240, r=16,
            loss_count=2, shard_bytes=4096, backend="avx512",
            decode_selection="path",
        )["metrics"]
        for metric in (
            "decode_reveal_direct_payload_bytes",
            "decode_reveal_inplace_temporary_payload_bytes",
            "decode_reveal_scatter_payload_bytes",
        ):
            self.assertEqual(avx512_fused[metric], fused[metric])

        ragged = build_decode_report(
            path="generic_decode", profile="high", k=240, r=16,
            loss_count=2, shard_bytes=4097, backend="avx2",
            decode_selection="path",
        )["metrics"]
        self.assertEqual(
            ragged["decode_reveal_direct_payload_bytes"]["value"], 8192
        )
        self.assertEqual(
            ragged["decode_reveal_inplace_temporary_payload_bytes"]["value"],
            128,
        )
        self.assertEqual(
            ragged["decode_reveal_scatter_payload_bytes"]["value"], 2
        )
        scalar = build_decode_report(
            path="generic_decode", profile="high", k=240, r=16,
            loss_count=2, shard_bytes=4096, backend="scalar",
            decode_selection="path",
        )["metrics"]
        self.assertEqual(
            scalar["decode_reveal_direct_payload_bytes"]["value"], 0
        )
        self.assertEqual(
            scalar["decode_reveal_inplace_temporary_payload_bytes"]["value"],
            8192,
        )
        self.assertEqual(
            scalar["decode_reveal_scatter_payload_bytes"]["value"], 8192
        )

        low_tail = build_decode_report(
            path="low_decode", k=8, r=248, field="gf16", backend="scalar",
            shard_bytes=130, loss_count=2, decode_workspace="tiled",
            decode_selection="path",
        )["metrics"]
        self.assertEqual(
            low_tail["decode_reveal_direct_payload_bytes"]["value"], 256
        )
        self.assertEqual(
            low_tail["decode_reveal_inplace_temporary_payload_bytes"]["value"],
            128,
        )
        self.assertEqual(
            low_tail["decode_reveal_scatter_payload_bytes"]["value"], 4
        )

        high_materialized = build_decode_report(
            loss_count=2, shard_bytes=128, decode_workspace="materialized",
            decode_selection="path",
        )["metrics"]
        high_tiled = build_decode_report(
            loss_count=2, shard_bytes=128, decode_workspace="tiled",
            decode_selection="path",
        )["metrics"]
        self.assertEqual(
            high_materialized[
                "decode_reveal_inplace_temporary_payload_bytes"
            ]["value"],
            0,
        )
        self.assertEqual(
            high_materialized["decode_reveal_scatter_payload_bytes"]["value"],
            0,
        )
        self.assertEqual(
            high_materialized["decode_reveal_direct_payload_bytes"]["value"],
            256,
        )
        self.assertEqual(
            high_tiled["decode_reveal_inplace_temporary_payload_bytes"]["value"],
            0,
        )
        self.assertEqual(
            high_tiled["decode_reveal_scatter_payload_bytes"]["value"], 0
        )
        self.assertEqual(
            high_tiled["decode_reveal_direct_payload_bytes"]["value"], 256
        )

        high_materialized_tail = build_decode_report(
            loss_count=2, shard_bytes=129, decode_workspace="materialized",
            decode_selection="path",
        )["metrics"]
        high_tiled_tail = build_decode_report(
            loss_count=2, shard_bytes=129, decode_workspace="tiled",
            decode_selection="path",
        )["metrics"]
        self.assertEqual(
            high_materialized_tail[
                "decode_reveal_inplace_temporary_payload_bytes"
            ]["value"],
            128,
        )
        self.assertEqual(
            high_tiled_tail[
                "decode_reveal_inplace_temporary_payload_bytes"
            ]["value"],
            0,
        )
        for metrics in (high_materialized_tail, high_tiled_tail):
            self.assertEqual(
                metrics["decode_reveal_direct_payload_bytes"]["value"], 256
            )
            self.assertEqual(
                metrics["decode_reveal_scatter_payload_bytes"]["value"], 2
            )

    def test_backend_qualified_high_syndrome_traffic(self) -> None:
        def metrics(shard_bytes: int, backend: str = "avx2") -> dict:
            return build_decode_report(
                k=1000, r=200, field="gf16", backend=backend,
                shard_bytes=shard_bytes, loss_count=8,
                decode_workspace="tiled", decode_selection="path",
            )["metrics"]

        compact = metrics(64)
        self.assertEqual(
            compact["decode_syndrome_fused_accumulation_payload_bytes"]["value"],
            0,
        )
        self.assertEqual(
            compact["decode_syndrome_materialized_xor_payload_bytes"]["value"],
            4 * 256 * 64,
        )
        split = metrics(192)
        self.assertEqual(
            split["decode_syndrome_fused_accumulation_payload_bytes"]["value"],
            2 * 256 * 192,
        )
        self.assertEqual(
            split["decode_syndrome_materialized_xor_payload_bytes"]["value"],
            2 * 256 * 192,
        )
        all_fused = metrics(1024)
        self.assertEqual(
            all_fused["decode_syndrome_fused_accumulation_payload_bytes"]["value"],
            4 * 256 * 1024,
        )
        self.assertEqual(
            all_fused["decode_syndrome_materialized_xor_payload_bytes"]["value"],
            0,
        )
        scalar = metrics(192, "scalar")
        self.assertEqual(
            scalar["decode_syndrome_fused_accumulation_payload_bytes"]["value"],
            0,
        )
        self.assertEqual(
            scalar["decode_syndrome_materialized_xor_payload_bytes"]["value"],
            4 * 256 * 192,
        )
        avx512 = metrics(192, "avx512")
        self.assertEqual(
            avx512["decode_syndrome_fused_accumulation_payload_bytes"],
            split["decode_syndrome_fused_accumulation_payload_bytes"],
        )
        self.assertEqual(
            avx512["decode_syndrome_materialized_xor_payload_bytes"],
            split["decode_syndrome_materialized_xor_payload_bytes"],
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
        self.assertEqual(direct.work_slot_bytes, 0)
        self.assertEqual(direct.plan_total_bytes, 640)
        self.assertGreater(direct.codec_total_bytes, direct.plan_total_bytes)
        no_op = COUNTS.decode_scratch_accounting(
            240, 16, 256, 16, "high", 65, 0, "tiled", 8,
            codec_workspace="tiled", no_op=True,
        )
        self.assertEqual(no_op.plan_total_bytes, 0)
        self.assertEqual(no_op.work_slot_bytes, 0)
        self.assertGreater(no_op.codec_total_bytes, 0)
        k1 = COUNTS.decode_scratch_accounting(
            1, 1, 2, 1, "high", 65, 1, "specialized", 8,
            codec_workspace="specialized", direct=True,
        )
        self.assertEqual(k1.plan_work_slots, 0)
        self.assertEqual(k1.codec_work_slots, 0)
        self.assertEqual(k1.work_slot_bytes, 0)
        self.assertEqual(k1.codec_work_slot_bytes, 0)

    def test_gf16_cache_policy_and_balanced_pass_geometry(self) -> None:
        mib = COUNTS.MIB
        expected = {
            0: (32, 16, 64),
            8: (32, 16, 64),
            32: (32, 16, 64),
            64: (64, 16, 64),
            96: (96, 16, 96),
            256: (96, 16, 96),
        }
        for detected_mib, policy_mib in expected.items():
            with self.subTest(detected_mib=detected_mib):
                policy = COUNTS.derive_gf16_cache_policy(
                    detected_mib * mib
                )
                self.assertEqual(
                    (
                        policy.effective_l3_bytes // mib,
                        policy.live_set_target_bytes // mib,
                        policy.tile_threshold_bytes // mib,
                    ),
                    policy_mib,
                )
        quantized = COUNTS.derive_gf16_cache_policy(64 * mib + 12345)
        self.assertEqual(quantized.effective_l3_bytes, 64 * mib)
        with self.assertRaises(COUNTS.ModelError):
            COUNTS.derive_gf16_cache_policy(-1)

        self.assertEqual(
            COUNTS.balanced_execution_tiles(128 * 1024, 48 * 1024),
            COUNTS.ExecutionTiles(3, 43712),
        )
        self.assertEqual(
            COUNTS.balanced_execution_tiles(1024 * 1024, 384 * 1024),
            COUNTS.ExecutionTiles(3, 349568),
        )
        self.assertEqual(
            COUNTS.balanced_execution_tiles(64 * 1024, 64 * 1024),
            COUNTS.ExecutionTiles(1, 64 * 1024),
        )

    def test_gf16_decode_cache_crossovers_and_backend_scope(self) -> None:
        mib = COUNTS.MIB
        # LOW_V1 side 64 has 128 retained rows.  These three payloads land
        # exactly on the 32/64/96-MiB policy crossovers.
        expected_passes = {
            (32, 512 * 1024): (4, 128 * 1024),
            (64, 512 * 1024): (4, 128 * 1024),
            (96, 512 * 1024): (1, 512 * 1024),
            (32, 768 * 1024): (6, 128 * 1024),
            (64, 768 * 1024): (6, 128 * 1024),
            (96, 768 * 1024): (6, 128 * 1024),
            (32, 1024 * 1024): (8, 128 * 1024),
            (64, 1024 * 1024): (8, 128 * 1024),
            (96, 1024 * 1024): (8, 128 * 1024),
        }
        for backend in ("avx2",):
            for (cache_mib, shard_bytes), expected in expected_passes.items():
                with self.subTest(
                    backend=backend, cache_mib=cache_mib,
                    shard_bytes=shard_bytes,
                ):
                    report = build_decode_report(
                        path="low_decode", k=64, r=193, field="gf16",
                        backend=backend, shard_bytes=shard_bytes,
                        loss_count=9, decode_selection="specialized",
                        gf16_detected_l3_bytes=cache_mib * mib,
                    )
                    metrics = report["metrics"]
                    self.assertEqual(
                        (
                            metrics["decode_plan_execution_tile_count"]["value"],
                            metrics["decode_plan_execution_tile_bytes"]["value"],
                        ),
                        expected,
                    )
                    self.assertEqual(
                        metrics["decode_work_slot_bytes"]["value"],
                        expected[1],
                    )
        gfni = build_decode_report(
            path="low_decode", k=64, r=193, field="gf16",
            backend="gfni", shard_bytes=768 * 1024,
            loss_count=9, decode_selection="specialized",
            gf16_detected_l3_bytes=96 * mib,
        )
        self.assertEqual(
            (
                gfni["metrics"][
                    "decode_plan_execution_tile_count"
                ]["value"],
                gfni["metrics"][
                    "decode_plan_execution_tile_bytes"
                ]["value"],
            ),
            (6, 128 * 1024),
        )
        self.assertEqual(
            gfni["inputs"]["gf16_cache_policy"],
            {
                "detected_l3_bytes": 0,
                "effective_l3_bytes": 32 * mib,
                "live_set_target_bytes": 16 * mib,
                "tile_threshold_bytes": 64 * mib,
            },
        )
        auto_avx2 = build_decode_report(
            path="low_decode", k=64, r=193, field="gf16",
            backend="avx2", shard_bytes=512 * 1024,
            loss_count=9, decode_selection="specialized",
            gf16_detected_l3_bytes=96 * mib,
            context_auto_requested=True,
        )
        self.assertEqual(
            (
                auto_avx2["metrics"][
                    "decode_plan_execution_tile_count"
                ]["value"],
                auto_avx2["metrics"][
                    "decode_plan_execution_tile_bytes"
                ]["value"],
            ),
            (4, 128 * 1024),
        )
        self.assertEqual(
            auto_avx2["inputs"]["gf16_cache_policy"][
                "effective_l3_bytes"
            ],
            32 * mib,
        )

        for backend in ("scalar", "ssse3", "avx512"):
            with self.subTest(excluded_backend=backend):
                report = build_decode_report(
                    path="low_decode", k=64, r=193, field="gf16",
                    backend=backend, shard_bytes=1024 * 1024,
                    loss_count=9, decode_selection="specialized",
                    gf16_detected_l3_bytes=32 * mib,
                )
                self.assertEqual(
                    report["metrics"][
                        "decode_plan_execution_tile_bytes"
                    ]["value"],
                    1024 * 1024,
                )

        low_side256 = {}
        for cache_mib in (32, 64, 96):
            report = build_decode_report(
                path="low_decode", k=200, r=800, field="gf16",
                backend="avx2", shard_bytes=64 * 1024,
                loss_count=9, decode_selection="specialized",
                gf16_detected_l3_bytes=cache_mib * mib,
            )
            low_side256[cache_mib] = (
                report["metrics"][
                    "decode_plan_execution_tile_count"
                ]["value"],
                report["metrics"][
                    "decode_plan_execution_tile_bytes"
                ]["value"],
            )
        self.assertEqual(
            low_side256,
            {
                32: (1, 64 * 1024),
                64: (1, 64 * 1024),
                96: (1, 64 * 1024),
            },
        )
        low_neighbors = {
            "loss": build_decode_report(
                path="low_decode", k=200, r=800, field="gf16",
                backend="avx2", shard_bytes=64 * 1024,
                loss_count=8, decode_selection="specialized",
                gf16_detected_l3_bytes=32 * mib,
            ),
            "bytes": build_decode_report(
                path="low_decode", k=200, r=800, field="gf16",
                backend="avx2", shard_bytes=96 * 1024,
                loss_count=9, decode_selection="specialized",
                gf16_detected_l3_bytes=32 * mib,
            ),
            "ragged": build_decode_report(
                path="low_decode", k=200, r=800, field="gf16",
                backend="avx2", shard_bytes=64 * 1024 + 2,
                loss_count=9, decode_selection="specialized",
                gf16_detected_l3_bytes=32 * mib,
            ),
            "auto": build_decode_report(
                path="low_decode", k=200, r=800, field="gf16",
                backend="avx2", shard_bytes=64 * 1024,
                loss_count=9, decode_selection="specialized",
                gf16_detected_l3_bytes=32 * mib,
                context_auto_requested=True,
            ),
            "count": build_decode_report(
                path="low_decode", k=201, r=799, field="gf16",
                backend="avx2", shard_bytes=64 * 1024,
                loss_count=9, decode_selection="specialized",
                gf16_detected_l3_bytes=32 * mib,
            ),
        }
        for neighbor, report in low_neighbors.items():
            with self.subTest(low_side256_neighbor=neighbor):
                prefix = (
                    report["inputs"]["shard_bytes"] //
                    COUNTS.SCRATCH_ALIGNMENT * COUNTS.SCRATCH_ALIGNMENT
                )
                self.assertEqual(
                    (
                        report["metrics"][
                            "decode_plan_execution_tile_count"
                        ]["value"],
                        report["metrics"][
                            "decode_plan_execution_tile_bytes"
                        ]["value"],
                    ),
                    (1, prefix),
                )

        gf8 = build_decode_report(
            path="low_decode", k=8, r=248, field="gf8",
            backend="avx2", shard_bytes=1024 * 1024,
            loss_count=8, decode_selection="specialized",
            gf16_detected_l3_bytes=96 * mib,
        )
        self.assertEqual(
            gf8["metrics"]["decode_plan_execution_tile_bytes"]["value"],
            1024 * 1024,
        )
        gf8_tiled = build_decode_report(
            k=192, r=64, field="gf8", backend="avx2",
            shard_bytes=256 * 1024, loss_count=9,
            decode_selection="specialized",
        )
        self.assertEqual(
            (
                gf8_tiled["metrics"][
                    "decode_plan_execution_tile_count"
                ]["value"],
                gf8_tiled["metrics"][
                    "decode_plan_execution_tile_bytes"
                ]["value"],
                gf8_tiled["metrics"][
                    "decode_codec_execution_tile_count"
                ]["value"],
                gf8_tiled["metrics"][
                    "decode_codec_execution_tile_bytes"
                ]["value"],
            ),
            (8, 32 * 1024, 8, 32 * 1024),
        )

    def test_gf16_plan_rows_and_conservative_codec_query(self) -> None:
        mib = COUNTS.MIB
        expected = {
            32: ((8, 57344), (7, 65536)),
            64: ((8, 57344), (7, 65536)),
            96: ((8, 57344), (7, 65536)),
        }
        for cache_mib, (plan_geometry, codec_geometry) in expected.items():
            with self.subTest(cache_mib=cache_mib):
                report = build_decode_report(
                    k=300, r=100, field="gf16", backend="avx2",
                    shard_bytes=448 * 1024, loss_count=9,
                    decode_selection="specialized",
                    gf16_detected_l3_bytes=cache_mib * mib,
                )
                metrics = report["metrics"]
                self.assertEqual(
                    metrics["decode_plan_work_slots"]["value"], 265
                )
                self.assertEqual(
                    metrics["decode_codec_work_slots"]["value"], 356
                )
                self.assertEqual(
                    (
                        metrics[
                            "decode_plan_execution_tile_count"
                        ]["value"],
                        metrics[
                            "decode_plan_execution_tile_bytes"
                        ]["value"],
                    ),
                    plan_geometry,
                )
                self.assertEqual(
                    (
                        metrics[
                            "decode_codec_execution_tile_count"
                        ]["value"],
                        metrics[
                            "decode_codec_execution_tile_bytes"
                        ]["value"],
                    ),
                    codec_geometry,
                )
                self.assertGreaterEqual(
                    metrics["decode_codec_scratch_bytes"]["value"],
                    metrics["decode_plan_scratch_bytes"]["value"],
                )

        # GFNI retains the conservative 32-MiB policy because only explicit
        # AVX2 contexts were qualified for affinity-derived cache sizing.
        # Both the known plan and plan-null query therefore take the existing
        # side-512 16-KiB override and the codec remains an upper bound.
        maximum_loss = build_decode_report(
            k=2000, r=512, field="gf16", backend="gfni",
            shard_bytes=64 * 1024, loss_count=512,
            decode_selection="specialized",
            gf16_detected_l3_bytes=96 * mib,
        )["metrics"]
        self.assertEqual(
            (
                maximum_loss["decode_plan_execution_tile_count"]["value"],
                maximum_loss["decode_plan_execution_tile_bytes"]["value"],
            ),
            (4, 16 * 1024),
        )
        self.assertEqual(
            (
                maximum_loss["decode_codec_execution_tile_count"]["value"],
                maximum_loss["decode_codec_execution_tile_bytes"]["value"],
            ),
            (4, 16 * 1024),
        )
        self.assertGreaterEqual(
            maximum_loss["decode_codec_scratch_bytes"]["value"],
            maximum_loss["decode_plan_scratch_bytes"]["value"],
        )

    def test_decode_scratch_declared_abi_boundaries(self) -> None:
        for pointer_bytes in (4, 8):
            maximum = COUNTS.size_t_max(pointer_bytes)
            largest_roundable = maximum - 63
            self.assertEqual(
                COUNTS.round_shard_bytes(largest_roundable, pointer_bytes),
                largest_roundable,
            )
            self.assertEqual(
                COUNTS.checked_size_add(maximum - 1, 1, maximum), maximum
            )
            self.assertEqual(
                COUNTS.checked_size_multiply(maximum, 1, maximum), maximum
            )
            with self.assertRaises(COUNTS.ModelError):
                COUNTS.checked_size_add(maximum, 1, maximum)
            with self.assertRaises(COUNTS.ModelError):
                COUNTS.checked_size_multiply(maximum, 2, maximum)
            for rejected in (maximum - 62, maximum - 1, maximum):
                with self.assertRaises(COUNTS.ModelError):
                    COUNTS.round_shard_bytes(rejected, pointer_bytes)

            # K=R=1 is a direct plan, but the pattern-independent codec query
            # still owns two transform work slots.  Pin the largest aligned
            # shard whose complete query fits size_t and reject its +64 neighbor.
            control = COUNTS.decode_scratch_accounting(
                1, 1, 2, 1, "low", 64, 1, "specialized",
                pointer_bytes, codec_workspace="specialized", direct=True,
            )
            largest_layout_shard = (
                (maximum - control.codec_data_offset) // control.codec_work_slots
            ) & ~63
            boundary = COUNTS.decode_scratch_accounting(
                1, 1, 2, 1, "low", largest_layout_shard, 1, "specialized",
                pointer_bytes, codec_workspace="specialized", direct=True,
            )
            self.assertLessEqual(boundary.codec_total_bytes, maximum)
            with self.assertRaises(COUNTS.ModelError):
                COUNTS.decode_scratch_accounting(
                    1, 1, 2, 1, "low", largest_layout_shard + 64, 1,
                    "specialized", pointer_bytes,
                    codec_workspace="specialized", direct=True,
                )

            with self.assertRaises(COUNTS.ModelError):
                COUNTS.decode_scratch_accounting(
                    1, 1, 2, 1, "low", COUNTS.UINT64_MAX, 1,
                    "specialized", pointer_bytes,
                    codec_workspace="specialized", direct=True,
                )

        with self.assertRaises(COUNTS.ModelError):
            COUNTS.round_shard_bytes(COUNTS.UINT64_MAX + 1, 8)
        with self.assertRaises(COUNTS.ModelError):
            COUNTS.decode_scratch_accounting(
                1, 1, 4, 1, "low", 64, 1, "specialized", 8,
                codec_workspace="specialized", direct=True,
            )
        with self.assertRaises(COUNTS.ModelError):
            COUNTS.decode_scratch_accounting(
                300, 100, 512, 128, "high", 64, 1, "tiled", 8,
                codec_workspace="tiled", field_name="gf8",
            )
        with self.assertRaises(COUNTS.ModelError):
            COUNTS.decode_scratch_accounting(
                1, 1, 2, 1, "low", 65, 1, "specialized", 8,
                codec_workspace="specialized", direct=True,
                field_name="gf16",
            )

        completed = subprocess.run(
            [sys.executable, str(TOOL), "report", "--path", "low_decode",
             "--k", "1", "--r", "1", "--field", "gf8",
             "--shard-bytes", str(COUNTS.UINT64_MAX), "--loss-count", "1",
             "--pointer-bytes", "8"],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            check=False,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("declared ABI size_t", completed.stderr)

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
            (
                "codec->context->gf16_live_set_target_bytes / work_rows",
                "codec->context->gf16_live_set_target_bytes / "
                "(work_rows + 1)",
            ),
            (
                "policy.live_set_target_bytes =\n"
                "        (kFallbackGF16L3Bytes / 2U)",
                "policy.live_set_target_bytes =\n"
                "        (kFallbackGF16L3Bytes / 4U)",
            ),
            (
                "std::max(effective, kFallbackGF16L3Bytes * 2U)",
                "std::max(effective, kFallbackGF16L3Bytes * 3U)",
            ),
            (
                "scaled /= kFallbackGF16L3Bytes / kMiB;",
                "scaled /= kMaximumCalibratedGF16L3Bytes / kMiB;",
            ),
            (
                "(aligned_bytes - 1U) / requested_tile_bytes",
                "(aligned_bytes + 1U) / requested_tile_bytes",
            ),
            (
                "case 64:\n"
                "            return aligned_prefix_bytes >= 192U * 1024U",
                "case 64:\n"
                "            return aligned_prefix_bytes >= 256U * 1024U",
            ),
            (
                "const bool cache_sensitive_gf16_backend =\n"
                "        requested == LEO2_BACKEND_AVX2;",
                "const bool cache_sensitive_gf16_backend =\n"
                "        effective_backend == LEO2_BACKEND_AVX2;",
            ),
        )
        for old, new in mutations:
            mutated = source.replace(old, new, 1)
            self.assertNotEqual(mutated, source)
            with self.assertRaises(COUNTS.ModelError):
                COUNTS.verify_decode_scratch_source(mutated, filename + " mutation")

    def test_decode_policy_source_guard_rejects_backend_mutations(self) -> None:
        filename = "Leopard2Dispatch.h"
        source = (ROOT / filename).read_text(encoding="utf-8")
        COUNTS.verify_decode_policy_source(source, filename)
        mutations = (
            (
                "backend == LEO2_BACKEND_AVX2 ||\n"
                "         backend == LEO2_BACKEND_AVX512",
                "backend == LEO2_BACKEND_AVX2",
            ),
            (
                "input.backend == LEO2_BACKEND_AVX2 ||\n"
                "             input.backend == LEO2_BACKEND_AVX512",
                "input.backend == LEO2_BACKEND_AVX2",
            ),
        )
        for old, new in mutations:
            mutated = source.replace(old, new, 1)
            self.assertNotEqual(mutated, source)
            with self.assertRaises(COUNTS.ModelError):
                COUNTS.verify_decode_policy_source(
                    mutated, filename + " backend mutation"
                )

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
        self.assertEqual(document["schema_version"], 4)
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
