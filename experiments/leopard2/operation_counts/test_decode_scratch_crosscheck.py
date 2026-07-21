#!/usr/bin/env python3
"""Cross-check schema-v3 decode policy/scratch accounting against production."""

from __future__ import annotations

import argparse
import importlib.util
import json
import pathlib
import subprocess
import sys
from typing import Dict, List


ROOT = pathlib.Path(__file__).resolve().parents[3]
TOOL = ROOT / "tools" / "leopard2_operation_counts.py"
SPEC = importlib.util.spec_from_file_location("leopard2_operation_counts", TOOL)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load operation-count tool")
COUNTS = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = COUNTS
SPEC.loader.exec_module(COUNTS)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def load_rows(probe: pathlib.Path) -> List[Dict[str, object]]:
    completed = subprocess.run(
        [str(probe)], text=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "decode scratch probe failed ({}):\n{}".format(
                completed.returncode, completed.stderr
            )
        )
    rows: List[Dict[str, object]] = []
    for line_number, line in enumerate(completed.stdout.splitlines(), 1):
        try:
            row = json.loads(line)
        except json.JSONDecodeError as error:
            raise AssertionError(
                "probe line {} is not JSON: {}".format(line_number, error)
            ) from error
        require(isinstance(row, dict), "probe row must be an object")
        rows.append(row)
    require(rows, "probe emitted no rows")
    return rows


def normalized_workspace(path: str) -> str:
    if path == "generic":
        return "generic"
    if path == "materialized":
        return "materialized"
    if path == "tiled":
        return "tiled"
    raise AssertionError("unexpected transform path: " + path)


def decode_force_from_flags(flags: int) -> str:
    if flags & 1:
        return "generic"
    if flags & 8:
        return "materialized"
    if flags & 4:
        return "tiled"
    if flags & 2:
        return "specialized"
    return "auto"


def cross_check_rows(rows: List[Dict[str, object]]) -> int:
    names = set()
    checks = 0
    auto_rows = 0
    for row in rows:
        name = str(row["name"])
        require(name not in names, "duplicate probe row " + name)
        names.add(name)
        require(
            row.get("schema") == "leopard2-decode-scratch-probe/v2",
            name + " probe schema changed",
        )
        selected_path = str(row["selected_path"])
        scratch_path = str(row["scratch_path"])
        codec_path = normalized_workspace(str(row["codec_path"]))
        direct = selected_path == "direct"
        no_op = selected_path == "no_op"
        profile = str(row["profile"])
        padded = int(row["padded"])
        parent = int(row["parent"])
        geometry = {
            "padded_side": padded,
            "parent_count": parent,
            "parent_dimension": (
                parent - padded if profile == "high" else padded
            ),
        }
        force = decode_force_from_flags(int(row["flags"]))
        selected_model = COUNTS.select_decode_execution(
            int(row["K"]), int(row["R"]), profile, str(row["field"]),
            str(row["backend"]), geometry, int(row["shard_bytes"]),
            int(row["losses"]), force,
            multi_item_batch=bool(row["multi_item_batch"]),
        )
        require(selected_model.path == selected_path,
                name + " selected policy path mismatch")
        require(selected_model.rule == str(row["selected_rule"]),
                name + " selected policy rule mismatch")
        require(
            selected_model.matching_auto_rules ==
            int(row["selected_matching_auto_rules"]),
            name + " selected matching-rule mask mismatch",
        )
        require(
            selected_model.required_work_slots ==
            int(row["selected_required_work_slots"]),
            name + " selected policy work-slot mismatch",
        )
        codec_model = COUNTS.select_codec_transform_execution(
            int(row["K"]), int(row["R"]), profile, str(row["field"]),
            str(row["backend"]), geometry, int(row["shard_bytes"]), force,
        )
        require(codec_model.path == codec_path,
                name + " codec-query policy path mismatch")
        require(codec_model.rule == str(row["codec_rule"]),
                name + " codec-query policy rule mismatch")
        require(
            codec_model.matching_auto_rules ==
            int(row["codec_matching_auto_rules"]),
            name + " codec-query matching-rule mask mismatch",
        )
        require(codec_model.required_work_slots == int(row["codec_work_slots"]),
                name + " codec-query policy work-slot mismatch")
        require(
            scratch_path == selected_path or bool(row["multi_item_batch"]),
            name + " single-item selected/scratch path differs",
        )
        plan_workspace = codec_path if (direct or no_op) else normalized_workspace(
            scratch_path
        )
        accounting = COUNTS.decode_scratch_accounting(
            int(row["K"]), int(row["R"]), int(row["parent"]),
            int(row["padded"]), str(row["profile"]),
            int(row["shard_bytes"]), int(row["losses"]),
            plan_workspace, int(row["pointer_bytes"]),
            codec_workspace=codec_path, no_op=no_op, direct=direct,
        )
        require(
            accounting.plan_total_bytes == int(row["plan_scratch_bytes"]),
            "{} plan scratch mismatch: model={} public={}".format(
                name, accounting.plan_total_bytes, row["plan_scratch_bytes"]
            ),
        )
        require(
            accounting.codec_total_bytes == int(row["codec_scratch_bytes"]),
            "{} codec scratch mismatch: model={} public={}".format(
                name, accounting.codec_total_bytes, row["codec_scratch_bytes"]
            ),
        )
        require(
            accounting.plan_work_slots == int(row["plan_work_slots"]),
            name + " plan work-slot mismatch",
        )
        require(
            accounting.codec_work_slots == int(row["codec_work_slots"]),
            name + " codec work-slot mismatch",
        )
        selected_accounting = COUNTS.decode_scratch_accounting(
            int(row["K"]), int(row["R"]), int(row["parent"]),
            int(row["padded"]), str(row["profile"]),
            int(row["shard_bytes"]), int(row["losses"]),
            (codec_path if (direct or no_op) else
             normalized_workspace(selected_path)),
            int(row["pointer_bytes"]), codec_workspace=codec_path,
            no_op=no_op, direct=direct,
        )
        require(
            selected_accounting.plan_work_slots ==
            int(row["selected_required_work_slots"]),
            name + " selected-path work-slot mismatch",
        )
        maximum = COUNTS.size_t_max(int(row["pointer_bytes"]))
        require(
            accounting.plan_total_bytes <= maximum and
            accounting.codec_total_bytes <= maximum,
            name + " public query escaped declared size_t",
        )
        require(
            COUNTS.round_shard_bytes(
                int(row["shard_bytes"]), int(row["pointer_bytes"])
            ) == ((int(row["shard_bytes"]) + 63) & ~63),
            name + " valid probe control did not round like production",
        )
        tail = int(row["shard_bytes"]) & 63
        prefix = int(row["shard_bytes"]) - tail
        require(int(row["tail_bytes"]) == tail,
                name + " introspected tail differs")
        require(int(row["aligned_prefix_bytes"]) == prefix,
                name + " introspected aligned prefix differs")
        require(int(row["rounded_bytes"]) == ((int(row["shard_bytes"]) + 63) & ~63),
                name + " introspected rounded bytes differ")
        if tail and not (direct or no_op):
            require(
                accounting.tail_reserved_slots == int(row["K"]) + int(row["R"]),
                name + " must reserve K+R ragged slots",
            )
            require(
                accounting.tail_selected_staged_slots == int(row["K"]),
                name + " must stage exactly K selected inputs",
            )
        else:
            require(
                accounting.tail_selected_staged_slots == 0,
                name + " must not stage transform inputs",
            )
        observed = bool(row["auto_host_backend_observed"])
        require(observed == name.startswith("auto_"), name + " AUTO label mismatch")
        if observed:
            auto_rows += 1
            require(
                str(row["backend"]) in
                ("scalar", "ssse3", "avx2", "avx512", "neon"),
                name + " did not record the selected production backend",
            )
        require(str(row["selected_rule"]) and str(row["scratch_rule"]) and
                str(row["codec_rule"]), name + " missing selected rule")
        checks += 24

    required = {
        "noop_ragged", "forced_high_materialized_aligned",
        "forced_high_materialized_ragged", "forced_high_tiled_aligned",
        "forced_high_tiled_ragged", "forced_generic_ragged",
        "forced_low_materialized_ragged", "forced_low_tiled_aligned",
        "direct_xor", "direct_copy", "direct_repair",
        "auto_ordinary_high", "auto_balanced_64", "auto_balanced_256",
        "auto_high_16k", "auto_high_32k", "auto_high_128k", "auto_low",
        "auto_high_32k_batch",
    }
    require(required.issubset(names), "probe is missing required rows")
    require(auto_rows >= 7, "probe did not cover AUTO size transitions")
    return checks + 2


def mutation_checks() -> int:
    source_path = ROOT / "leopard2.cpp"
    source = source_path.read_text(encoding="utf-8")
    COUNTS.verify_decode_scratch_source(source, "leopard2.cpp")
    mutations = (
        (
            "static_cast<size_t>(codec->original_count) * 2 + codec->recovery_count",
            "static_cast<size_t>(codec->original_count) + codec->recovery_count",
            "range-count mutation",
        ),
        (
            "static_cast<size_t>(codec->original_count) + codec->recovery_count\n        : 0;",
            "static_cast<size_t>(codec->original_count)\n        : 0;",
            "ragged-reservation mutation",
        ),
        (
            "plan\n            ? plan->requested_coordinates.size()\n            : codec->recovery_count",
            "plan\n            ? codec->recovery_count\n            : codec->recovery_count",
            "plan-output-slot mutation",
        ),
        (
            "ComputeScratchLayout(range_count, 0, 0, rounded_bytes, layout)",
            "ComputeScratchLayout(range_count, 1, 0, rounded_bytes, layout)",
            "direct-pointer mutation",
        ),
        (
            "static_cast<size_t>(codec->parent_count),\n            geometry.work_slot_count,\n            pointer_count",
            "static_cast<size_t>(codec->original_count),\n            geometry.work_slot_count,\n            pointer_count",
            "coordinate-pointer mutation",
        ),
    )
    checks = 1
    for old, new, label in mutations:
        mutated = source.replace(old, new, 1)
        require(mutated != source, label + " did not apply")
        try:
            COUNTS.verify_decode_scratch_source(mutated, label)
        except COUNTS.ModelError:
            checks += 1
        else:
            raise AssertionError(label + " escaped the source guard")

    # Formula-level sentinels catch plausible but wrong rewrites that can look
    # reasonable in source review: rounded work slots, K-only tail reservation,
    # and R output slots in a pattern-specific high tiled plan.
    exact = COUNTS.decode_scratch_accounting(
        240, 16, 256, 16, "high", 65, 1, "tiled", 8,
        codec_workspace="tiled",
    )
    require(exact.work_slot_bytes == 64, "ragged work-slot sentinel")
    require(exact.tail_reserved_slots == 256, "K+R tail sentinel")
    require(exact.plan_work_slots == 33, "high plan L-slot sentinel")
    require(exact.codec_work_slots == 48, "high codec R-slot sentinel")
    ff8 = (ROOT / "LeopardFF8.cpp").read_text(encoding="utf-8")
    ff16 = (ROOT / "LeopardFF16.cpp").read_text(encoding="utf-8")
    COUNTS.verify_decode_fusion_sources(source, ff8, ff16)
    fusion_mutations = (
        (
            source.replace("aligned_prefix_bytes >= 4096",
                           "aligned_prefix_bytes >= 2048", 1),
            ff8,
            ff16,
            "reveal crossover mutation",
        ),
        (
            source,
            ff8.replace("buffer_bytes >= 65536",
                        "buffer_bytes >= 32768", 1),
            ff16,
            "GF8 sink crossover mutation",
        ),
        (
            source,
            ff8,
            ff16.replace("bytes == 64", "bytes == 32", 1),
            "GF16 fused-four mutation",
        ),
        (
            source.replace(
                "const bool reveal_aligned_outputs_in_place =\n"
                "        !(fuse_generic_reveal_scatter || "
                "fuse_low_reveal_scatter);",
                "const bool reveal_aligned_outputs_in_place = true;",
                1,
            ),
            ff8,
            ff16,
            "aligned reveal wiring mutation",
        ),
        (
            source.replace(
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
            "tail reveal wiring mutation",
        ),
        (
            source,
            ff8.replace(
                "const bool use_accumulating_sink =\n"
                "                ShouldUsePrunedHighSyndromeSink("
                "ops, buffer_bytes);",
                "const bool use_accumulating_sink = false;",
                1,
            ),
            ff16,
            "GF8 pruned sink wiring mutation",
        ),
        (
            source,
            ff8,
            ff16.replace(
                "const bool use_accumulating_sink =\n"
                "                ShouldUsePrunedHighSyndromeSink("
                "ops, buffer_bytes);",
                "const bool use_accumulating_sink = false;",
                1,
            ),
            "GF16 pruned sink wiring mutation",
        ),
    )
    checks += 1
    for mutated_core, mutated_ff8, mutated_ff16, label in fusion_mutations:
        require(
            mutated_core != source or mutated_ff8 != ff8 or
            mutated_ff16 != ff16,
            label + " did not apply",
        )
        try:
            COUNTS.verify_decode_fusion_sources(
                mutated_core, mutated_ff8, mutated_ff16
            )
        except COUNTS.ModelError:
            checks += 1
        else:
            raise AssertionError(label + " escaped the source guard")
    return checks + 4


def abi_boundary_checks() -> int:
    checks = 0
    for pointer_bytes in (4, 8):
        maximum = COUNTS.size_t_max(pointer_bytes)
        largest_roundable = maximum - 63
        require(
            COUNTS.round_shard_bytes(largest_roundable, pointer_bytes) ==
            largest_roundable,
            "{}-bit largest roundable shard mismatch".format(pointer_bytes * 8),
        )
        checks += 1
        rejected_boundaries = dict.fromkeys(
            (maximum - 62, maximum - 1, maximum, COUNTS.UINT64_MAX)
        )
        for rejected in rejected_boundaries:
            try:
                COUNTS.round_shard_bytes(rejected, pointer_bytes)
            except COUNTS.ModelError:
                checks += 1
            else:
                raise AssertionError(
                    "{}-bit shard boundary was accepted: {}".format(
                        pointer_bytes * 8, rejected
                    )
                )

        control = COUNTS.decode_scratch_accounting(
            1, 1, 2, 1, "low", 64, 1, "specialized", pointer_bytes,
            codec_workspace="specialized", direct=True,
        )
        largest_layout_shard = (
            (maximum - control.codec_data_offset) // control.codec_work_slots
        ) & ~63
        boundary = COUNTS.decode_scratch_accounting(
            1, 1, 2, 1, "low", largest_layout_shard, 1, "specialized",
            pointer_bytes, codec_workspace="specialized", direct=True,
        )
        require(
            boundary.codec_total_bytes <= maximum,
            "{}-bit valid layout boundary overflowed".format(pointer_bytes * 8),
        )
        try:
            COUNTS.decode_scratch_accounting(
                1, 1, 2, 1, "low", largest_layout_shard + 64, 1,
                "specialized", pointer_bytes,
                codec_workspace="specialized", direct=True,
            )
        except COUNTS.ModelError:
            checks += 2
        else:
            raise AssertionError(
                "{}-bit layout overflow neighbor was accepted".format(
                    pointer_bytes * 8
                )
            )

    try:
        COUNTS.round_shard_bytes(COUNTS.UINT64_MAX + 1, 8)
    except COUNTS.ModelError:
        return checks + 1
    raise AssertionError("value above uint64_t was accepted")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", type=pathlib.Path, required=True)
    args = parser.parse_args()
    checks = (cross_check_rows(load_rows(args.probe)) + mutation_checks() +
              abi_boundary_checks())
    print("decode scratch cross-check passed: {} checks".format(checks))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (AssertionError, OSError, RuntimeError) as error:
        print("decode scratch cross-check failed: " + str(error), file=sys.stderr)
        sys.exit(1)
