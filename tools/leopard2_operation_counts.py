#!/usr/bin/env python3
"""Deterministic operation and logical-traffic model for Leopard2 paths.

This utility models the transform schedules implemented by Leopard's radix-4
wrappers in terms of their equivalent radix-2 butterflies.  It deliberately
does not pretend that a scheduled butterfly always executes a multiplication:
zero skew factors are specialized by the production kernels.  Arithmetic
totals involving those factors are therefore emitted as bounds.  The
structural schedule remains ISA-independent; separately classified traffic
metrics apply the production backend's decode-fusion predicates.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import pathlib
import re
import struct
import sys
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, TextIO, Tuple


SCHEMA_VERSION = 3
UINT32_MAX = (1 << 32) - 1
UINT64_MAX = (1 << 64) - 1
SCRATCH_ALIGNMENT = 64
PATHS = (
    "legacy_high_encode",
    "legacy_high_decode",
    "low_encode",
    "low_decode",
    "generic_transform",
    "generic_decode",
    "direct_repair",
)
FIELDS = ("gf8", "gf16")
PROFILES = ("high", "low")
DECODE_WORKSPACES = ("materialized", "tiled")
BACKENDS = ("scalar", "ssse3", "avx2", "avx512", "neon")
DECODE_SELECTIONS = ("path", "auto", "specialized")


_LOW_COEFFICIENT_COPY_PATTERNS = (
    re.compile(
        r"(?:std::)?memcpy\(work\[p\+[^\]]+\],work\[[^\]]+\],buffer_bytes\)"
    ),
    re.compile(
        r"(?:std::)?memcpy\(evaluation_work\[[^\]]+\],"
        r"coefficients\[[^\]]+\],bytes\)"
    ),
)


_PREPROCESSOR_DIRECTIVE = re.compile(
    r"^[ \t]*#[ \t]*(if|ifdef|ifndef|elif|else|endif)\b([^\r\n]*)",
    re.MULTILINE,
)


def _compact_source_with_offsets(source: str) -> Tuple[str, List[int]]:
    """Return whitespace-free source and each retained character's offset."""
    characters: List[str] = []
    offsets: List[int] = []
    for offset, character in enumerate(source):
        if character.isspace():
            continue
        characters.append(character)
        offsets.append(offset)
    return "".join(characters), offsets


def _canonical_test_hook_guard_start(source: str, offset: int) -> Optional[int]:
    """Find an enclosing positive, canonical test-hook preprocessor branch.

    This is intentionally a small structural preprocessor parser rather than a
    nearby-text heuristic.  In particular, an unrelated nested ``#endif`` must
    not hide the still-active outer test-hook guard, while a copy after the
    outer ``#endif`` (or in its ``#else`` arm) must never be exempted.
    """
    stack: List[Tuple[bool, bool, int]] = []
    for directive in _PREPROCESSOR_DIRECTIVE.finditer(source):
        if directive.start() >= offset:
            break
        operation = directive.group(1)
        expression = directive.group(2).strip()
        if operation in ("if", "ifdef", "ifndef"):
            canonical = (
                operation == "if" and
                re.fullmatch(
                    r"defined[ \t]*\([ \t]*LEO2_ENABLE_TEST_HOOKS[ \t]*\)",
                    expression,
                ) is not None
            )
            stack.append((canonical, canonical, directive.start()))
        elif operation in ("elif", "else"):
            if stack:
                canonical, _, start = stack[-1]
                stack[-1] = (canonical, False, start)
        elif operation == "endif" and stack:
            stack.pop()
    starts = [start for canonical, active, start in stack
              if canonical and active]
    return max(starts) if starts else None


def _is_forced_high_decode_test_copy(
    source: str, original_offset: int
) -> bool:
    guard_start = _canonical_test_hook_guard_start(source, original_offset)
    if guard_start is None:
        return False
    guarded_prefix = re.sub(r"\s+", "", source[guard_start:original_offset])
    return (
        "callsite==SourceEvaluationHighDecode&&" in guarded_prefix and
        "TestForceHighDecodeCopyFallback" in guarded_prefix
    )


class ModelError(ValueError):
    """Invalid or unsupported model input."""


def verify_low_encode_no_copy_source(source: str, label: str) -> None:
    """Reject the former whole-P coefficient-copy pass.

    The executable backend-call tests prove which first-layer implementation
    ran.  This independent source guard covers the otherwise unobservable
    regression where a whole coefficient block is copied before that call.
    Whitespace is removed so the known loop cannot evade the check through
    formatting changes.
    """
    compact, offsets = _compact_source_with_offsets(source)
    for pattern_index, pattern in enumerate(_LOW_COEFFICIENT_COPY_PATTERNS):
        for match in pattern.finditer(compact):
            # The private Algorithm-5 attribution build deliberately restores
            # this copy behind the test-hooks preprocessor boundary.  It is a
            # high-decode control, not the low encoder regression guarded
            # here; test_hook_isolation.cmake separately proves it is absent
            # from the production archive.
            forced_high_decode_test_copy = (
                pattern_index == 1 and
                _is_forced_high_decode_test_copy(
                    source, offsets[match.start()])
            )
            if forced_high_decode_test_copy:
                continue
            raise ModelError(
                "{} contains a whole-P low-encode coefficient copy".format(label)
            )
    call = (
        "FFT_DIT_FromCoefficients(ops,buffer_bytes,work,work+p,"
        "requested_count,p,"
    )
    if call not in compact:
        raise ModelError(
            "{} no longer routes low parity evaluation through the "
            "out-of-place first layer".format(label)
        )


_HIGH_EVALUATOR_COPY_PATTERNS = (
    re.compile(r"memcpy\(work\[offset\+i\],work\[i\],buffer_bytes\)"),
    re.compile(r"memcpy\(tile\[i\],accumulator\[i\],buffer_bytes\)"),
)


def verify_high_decode_no_copy_source(source: str, label: str) -> None:
    """Reject an unconditional whole-T Algorithm 5 evaluator copy."""
    compact = re.sub(r"\s+", "", source)
    # The compatibility fallback is deliberately retained after a failed
    # trusted-plan execution.  Only copies before selection of the immutable-
    # source executor would restore the former production hot path.
    for pattern in _HIGH_EVALUATOR_COPY_PATTERNS:
        for match in pattern.finditer(compact):
            nearby = compact[max(0, match.start() - 256):match.start()]
            if "if(!executed){" not in nearby:
                raise ModelError(
                    "{} contains an unconditional whole-T high-decode "
                    "evaluator copy".format(label)
                )
    if "ExecutePrunedForwardTransformPlanFromSources(" not in compact:
        raise ModelError(
            "{} no longer routes sparse Algorithm 5 evaluation through the "
            "immutable-source executor".format(label)
        )
    if "FFT_DIT_FromCoefficients(" not in compact:
        raise ModelError(
            "{} no longer routes mature Algorithm 5 evaluation through the "
            "out-of-place first layer".format(label)
        )


def verify_high_decode_receive_fusion_source(source: str, label: str) -> None:
    """Require mature Algorithm 5 input blocks to use source staging."""
    compact = re.sub(r"\s+", "", source)
    begin = compact.find("staticvoidIFFT_DIT_DecoderFromSources(")
    end = compact.find("/*DecimationintimeFFT:", begin)
    if begin < 0 or end < 0:
        raise ModelError(
            "{} no longer defines the Algorithm 5 source-boundary inverse "
            "transform".format(label)
        )
    helper = compact[begin:end]
    if "ops.ff8_ifft_butterfly4_out(" not in helper:
        raise ModelError(
            "{} no longer enters the GF8 out-of-place inverse boundary".format(
                label
            )
        )
    qualified_backends = (
        "constboolqualified_source_backend="
        "ops.kind==LEO2_BACKEND_SSSE3||"
        "ops.kind==LEO2_BACKEND_AVX2||"
        "ops.kind==LEO2_BACKEND_AVX512;"
    )
    if qualified_backends not in helper or \
            "!qualified_source_backend" not in helper:
        raise ModelError(
            "{} no longer confines GF8 receive fusion to qualified SIMD "
            "backends".format(label)
        )
    materialized = compact.find(
        "voidReedSolomonDecodeHighPrunedPlanned(constbackend::Ops&"
    )
    tiled = compact.find(
        "voidReedSolomonDecodeHighTiledPrunedPlanned(constbackend::Ops&"
    )
    if materialized < 0 or tiled < 0:
        raise ModelError("{} is missing a planned Algorithm 5 path".format(label))
    materialized_end = compact.find(
        "voidReedSolomonDecodeHighPlanned(", materialized
    )
    tiled_end = compact.find("voidReedSolomonDecodeHighTiledPlanned(", tiled)
    if materialized_end < 0 or tiled_end < 0:
        raise ModelError("{} has an unbounded Algorithm 5 source check".format(label))
    for name, segment in (
        ("materialized", compact[materialized:materialized_end]),
        ("tiled", compact[tiled:tiled_end]),
    ):
        if "IFFT_DIT_DecoderFromSources(" not in segment:
            raise ModelError(
                "{} {} Algorithm 5 path restored copy-first receive staging".format(
                    label, name
                )
            )
        if "ExecutePrunedInverseTransformPlanAccumulate(" not in segment or \
                "StageHighDecodeSources(" not in segment:
            raise ModelError(
                "{} {} Algorithm 5 path lost the exact-pruned sink union".format(
                    label, name
                )
            )


def verify_high_decode_weighted_locator_source(source: str, label: str) -> None:
    """Bind the Algorithm 5 locator delta to its measured dispatch region."""
    compact = re.sub(r"\s+", "", source)
    begin = compact.find("staticvoidIFFT_DIT_DecoderWeightedLocator(")
    end = compact.find("staticvoidStageHighDecodeSources(", begin)
    if begin < 0 or end < 0:
        raise ModelError(
            "{} is missing the Algorithm 5 weighted locator boundary".format(
                label
            )
        )
    helper = compact[begin:end]
    required = (
        "(ops.kind==LEO2_BACKEND_AVX2||"
        "ops.kind==LEO2_BACKEND_AVX512)",
        "ops.ff8_weighted_ifft_butterfly4!=NULL",
        "m>=64",
        "bytes>=16U*1024U",
        "bytes<=256U*1024U",
        "live_count>=(m+1U)/2U",
        "mul_mem_inplace(ops,work[i],locator_logs[i],bytes)",
        "IFFT_DIT_Decoder(ops,bytes,m,work,m,skewLUT)",
        "ops.ff8_weighted_ifft_butterfly4(",
    )
    for token in required:
        if token not in helper:
            raise ModelError(
                "{} widened or obscured the qualified weighted locator "
                "boundary ({})".format(label, token)
            )
    for callsite in (
        "voidReedSolomonDecodeHighPrepared(constbackend::Ops&",
        "voidReedSolomonDecodeHighPrunedPlanned(constbackend::Ops&",
        "voidReedSolomonDecodeHighTiledPrunedPlanned(constbackend::Ops&",
    ):
        start = compact.find(callsite)
        if start < 0:
            raise ModelError("{} lost {}".format(label, callsite))
        body = compact.find("{", start + len(callsite))
        if body < 0:
            raise ModelError("{} has no body".format(callsite))
        depth = 0
        end = body
        for end in range(body, len(compact)):
            if compact[end] == "{":
                depth += 1
            elif compact[end] == "}":
                depth -= 1
                if depth == 0:
                    break
        if depth != 0:
            raise ModelError("{} has an unbalanced body".format(callsite))
        segment = compact[start:end + 1]
        if "IFFT_DIT_DecoderWeightedLocator(" not in segment:
            raise ModelError(
                "{} callsite bypasses the weighted locator dispatcher".format(
                    callsite
                )
            )


def verify_high_decode_gf16_copy_first_source(source: str, label: str) -> None:
    """Reject promotion of the measured-below-threshold GF16 candidate."""
    compact = re.sub(r"\s+", "", source)
    begin = compact.find("staticvoidStageHighDecodeSources(")
    end = compact.find("/*DecimationintimeFFT:", begin)
    if begin < 0 or end < 0:
        raise ModelError("{} is missing the GF16 staging helper".format(label))
    segment = compact[begin:end]
    if "ff16_ifft_butterfly4_out(" in segment:
        raise ModelError(
            "{} promoted the below-threshold GF16 source candidate".format(label)
        )

    materialized = compact.find(
        "voidReedSolomonDecodeHighPrunedPlanned(constbackend::Ops&"
    )
    materialized_end = compact.find(
        "voidReedSolomonDecodeHighPlanned(", materialized
    )
    if materialized < 0 or materialized_end < 0:
        raise ModelError("{} is missing materialized Algorithm 5".format(label))
    materialized_segment = compact[materialized:materialized_end]
    stage = "StageHighDecodeSources(buffer_bytes,coordinate_data,work,n);"
    loop = materialized_segment.find("for(unsignedblock=0;")
    stage_at = materialized_segment.find(stage)
    if stage_at < 0 or loop < 0 or stage_at > loop or \
            materialized_segment.count("StageHighDecodeSources(") != 1:
        raise ModelError(
            "{} no longer stages materialized GF16 in one parent-wide pass".format(
                label
            )
        )
    if "IFFT_DIT_DecoderFromSources(" in materialized_segment:
        raise ModelError(
            "{} reopened per-block staging in materialized GF16".format(label)
        )

    tiled = compact.find(
        "voidReedSolomonDecodeHighTiledPrunedPlanned(constbackend::Ops&"
    )
    tiled_end = compact.find("voidReedSolomonDecodeHighTiledPlanned(", tiled)
    if tiled < 0 or tiled_end < 0:
        raise ModelError("{} is missing tiled Algorithm 5".format(label))
    tiled_segment = compact[tiled:tiled_end]
    if "StageHighDecodeSources(buffer_bytes,coordinate_data,accumulator,t);" \
            not in tiled_segment or \
            "StageHighDecodeSources(buffer_bytes,coordinate_data+offset,tile,t);" \
            not in tiled_segment:
        raise ModelError(
            "{} no longer retains per-tile GF16 copy-first staging".format(label)
        )
    if "ExecutePrunedInverseTransformPlanAccumulate(" not in tiled_segment:
        raise ModelError(
            "{} tiled GF16 path lost the exact-pruned sink".format(label)
        )


def verify_high_pruned_hook_registration(source: str, label: str) -> None:
    """Require the GF16 parent-wide counter gate to execute with hooks."""
    compact = re.sub(r"\s+", "", source)
    begin = compact.find(
        "add_executable(leopard2_high_pruned_stage_test"
    )
    end = compact.find(
        "add_executable(leopard2_decode_low_acceptance_test", begin
    )
    if begin < 0 or end < 0:
        raise ModelError(
            "{} no longer registers the high-pruned stage counter target".format(
                label
            )
        )
    segment = compact[begin:end]
    required = (
        "target_compile_definitions(leopard2_high_pruned_stage_testPRIVATE"
        "LEO2_ENABLE_TEST_HOOKS=1",
        "LEO2_REQUIRE_HIGH_PRUNED_STAGE_HOOKS=1)",
        "target_link_libraries(leopard2_high_pruned_stage_test"
        "leopard_test_hooks)",
        "add_test(NAMEleopard2_high_pruned_stageCOMMAND"
        "leopard2_high_pruned_stage_test)",
    )
    if any(token not in segment for token in required):
        raise ModelError(
            "{} can compile out or bypass the high-pruned stage counters".format(
                label
            )
        )
    if "install(TARGETSleopard2_high_pruned_stage_test" in compact:
        raise ModelError("{} installs a private hook executable".format(label))


def verify_low_decode_weighted_fusion_source(source: str, label: str) -> None:
    """Require Algorithm 4's weighted reduction to use one mul-add pass."""
    compact = re.sub(r"\s+", "", source)
    materialized = (
        "muladd_mem(ops,work[i],work[offset+i],"
        "block_factors[block-1],buffer_bytes);"
    )
    tiled = (
        "muladd_mem(ops,accumulator[i],tile[i],"
        "block_factors[block-1],buffer_bytes);"
    )
    if compact.count(materialized) < 2 or tiled not in compact:
        raise ModelError(
            "{} no longer fuses every Algorithm 4 weighted block reduction "
            "into multiply-add".format(label)
        )


def verify_decode_scratch_source(source: str, label: str) -> None:
    """Guard the production layout boundaries mirrored by schema v2.

    Runtime cross-checks compare public query bytes.  This scoped source guard
    additionally makes the intended distinctions reviewable: N coordinate
    pointers are metadata, ragged K+R input slots are fixed 64-byte staging,
    and direct execution retains range metadata without shard-data slots.
    """
    compact = re.sub(r"\s+", "", source)
    selector_begin = compact.find("staticboolSelectTransformDecodePath(")
    begin = compact.find("staticleo2_resultDecodeLayout(", selector_begin)
    end = compact.find("staticleo2_resultDirectDecodeLayout(", begin)
    direct_end = compact.find("staticboolUseSingleSideEncodeLayout(", end)
    populate = compact.find("staticvoidPopulateDecodeCoordinates(")
    populate_end = compact.find("staticLEO_FORCE_INLINEvoidGatherTransformDecodeOne(", populate)
    execute = compact.find("staticleo2_resultDecodePlanExecuteInternal(", direct_end)
    if min(selector_begin, begin, end, direct_end, populate, populate_end,
           execute) < 0:
        raise ModelError("{} is missing a decode scratch boundary".format(label))
    selector = compact[selector_begin:begin]
    layout = compact[begin:end]
    direct = compact[end:direct_end]
    staging = compact[populate:populate_end]
    execution = compact[execute:]
    required_selector = (
        "plan?plan->requested_coordinates.size():"
        "codec->recovery_count",
        "input.actual_shard_bytes=shard_bytes;",
        "input.aligned_prefix_bytes=geometry.aligned_prefix_bytes;",
        "input.tail_bytes=geometry.tail_bytes;",
        "input.rounded_shard_bytes=geometry.rounded_bytes;",
        "input.multi_item_batch=multi_item_batch;",
        "returnleopard2_internal::SelectDecodePath(input,selection);",
    )
    if any(token not in selector for token in required_selector):
        raise ModelError(
            "{} decode selector no longer owns complete byte/batch state".format(
                label
            )
        )
    required_layout = (
        "geometry.work_slot_bytes=geometry.tail_bytes==0?"
        "static_cast<size_t>(shard_bytes):"
        "std::max(geometry.aligned_prefix_bytes,kScratchAlignment);",
        "constsize_trange_count=static_cast<size_t>(codec->original_count)*2+"
        "codec->recovery_count;",
        "static_cast<size_t>(codec->parent_count),"
        "geometry.work_slot_count,pointer_count",
        "constsize_tinput_slot_count=geometry.tail_bytes!=0?"
        "static_cast<size_t>(codec->original_count)+codec->recovery_count:0;",
        "ComputeSplitScratchLayout(range_count,pointer_count,input_slot_count,"
        "geometry.work_slot_count,geometry.work_slot_bytes,geometry.layout,"
        "geometry.work_data_offset)",
        "SelectTransformDecodePath(codec,plan,shard_bytes,multi_item_batch,"
        "geometry,geometry.selection)",
        "geometry.work_slot_count=geometry.selection.required_work_slots;",
    )
    if any(token not in layout for token in required_layout):
        raise ModelError(
            "{} decode transform scratch formula diverged from schema v2".format(
                label
            )
        )
    required_direct = (
        "constsize_trange_count=static_cast<size_t>(codec->original_count)*2+"
        "codec->recovery_count;",
        "ComputeScratchLayout(range_count,0,0,rounded_bytes,layout)",
    )
    if any(token not in direct for token in required_direct):
        raise ModelError(
            "{} direct decode no longer uses range-only scratch".format(label)
        )
    required_staging = (
        "constboolstage_inputs=staging_slots!=NULL;",
        "static_cast<size_t>(i)*kScratchAlignment",
        "(static_cast<size_t>(codec->original_count)+i)*kScratchAlignment",
        "StageShardForKernel(codec,slot,source,pass_bytes,kScratchAlignment);",
    )
    if any(token not in staging for token in required_staging):
        raise ModelError(
            "{} ragged decode staging map diverged from schema v2".format(label)
        )
    required_execution = (
        "constbooluse_generic=geometry.selection.path=="
        "leopard2_internal::kDecodePathGeneric;",
        "constbooluse_tiled=geometry.selection.path=="
        "leopard2_internal::kDecodePathTiled;",
    )
    if any(token not in execution for token in required_execution):
        raise ModelError(
            "{} execution duplicated or bypassed the selected decode path".format(
                label
            )
        )


def verify_decode_policy_source(source: str, label: str) -> None:
    """Bind the selector mirror to backend-qualified production predicates."""
    compact = re.sub(r"\s+", "", source)
    balanced_begin = compact.find(
        "staticinlineboolShouldUseBalancedGenericDecode("
    )
    materialized_begin = compact.find(
        "staticinlineboolShouldUseMaterializedHighDecode("
    )
    geometry_begin = compact.find(
        "staticinlineboolIsDecodeByteGeometryConsistent("
    )
    selector_begin = compact.find("staticinlineboolSelectDecodePath(")
    selector_end = compact.find(
        "staticinlineconstchar*DecodePathName(", selector_begin
    )
    if min(balanced_begin, materialized_begin, geometry_begin,
           selector_begin, selector_end) < 0 or \
            not (balanced_begin < materialized_begin < geometry_begin):
        raise ModelError("{} is missing a decode policy boundary".format(label))

    balanced = compact[balanced_begin:materialized_begin]
    materialized = compact[materialized_begin:geometry_begin]
    selector = compact[selector_begin:selector_end]
    balanced_backend = (
        "backend==LEO2_BACKEND_SCALAR||"
        "backend==LEO2_BACKEND_SSSE3||"
        "backend==LEO2_BACKEND_AVX2||"
        "backend==LEO2_BACKEND_AVX512);"
    )
    materialized_backend = (
        "if(backend==LEO2_BACKEND_AVX2||"
        "backend==LEO2_BACKEND_AVX512)"
        "returnrounded_shard_bytes>=24*1024;"
    )
    batch_backend = (
        "if(input.multi_item_batch&&"
        "(input.backend==LEO2_BACKEND_AVX2||"
        "input.backend==LEO2_BACKEND_AVX512)){"
    )
    if balanced_backend not in balanced:
        raise ModelError(
            "{} balanced decode backends diverged from schema v3".format(label)
        )
    if materialized_backend not in materialized:
        raise ModelError(
            "{} materialized decode backends diverged from schema v3".format(
                label
            )
        )
    if batch_backend not in selector:
        raise ModelError(
            "{} batch decode backends diverged from schema v3".format(label)
        )


def verify_decode_fusion_sources(
    core_source: str,
    ff8_source: str,
    ff16_source: str,
) -> None:
    """Bind schema-v3 traffic predicates to the production source.

    The operation model is intentionally independent code, but backend traffic
    would become misleading if a production crossover changed without a schema
    update.  These narrow guards reject that drift and are mutation-tested.
    """
    core = re.sub(r"\s+", "", core_source)
    ff8 = re.sub(r"\s+", "", ff8_source)
    ff16 = re.sub(r"\s+", "", ff16_source)
    core_tokens = (
        "returncodec->field==LEO2_FIELD_GF8&&"
        "aligned_prefix_bytes>=4096&&"
        "(codec->context->backend==LEO2_BACKEND_SSSE3||"
        "codec->context->backend==LEO2_BACKEND_AVX2||"
        "codec->context->backend==LEO2_BACKEND_AVX512);",
        "returncodec->profile==LEO2_PROFILE_LOW_V1&&"
        "aligned_prefix_bytes!=0;",
        "constboolfuse_generic_reveal_scatter=use_generic&&"
        "UseFusedGenericRevealScatter(codec,geometry.aligned_prefix_bytes);",
        "constboolfuse_low_reveal_scatter=!use_generic&&"
        "UseFusedLowRevealScatter(codec,geometry.aligned_prefix_bytes);",
        "constboolreveal_aligned_outputs_in_place="
        "!(fuse_generic_reveal_scatter||fuse_low_reveal_scatter);",
        "ExecuteTransformDecodePass(plan,geometry.aligned_prefix_bytes,"
        "coordinate_input,work,use_generic,use_tiled,"
        "reveal_aligned_outputs_in_place);",
        "GatherTransformDecodePass(plan,restored_original,0,"
        "geometry.aligned_prefix_bytes,work,use_generic,use_tiled,"
        "reveal_aligned_outputs_in_place);",
        "ExecuteTransformDecodePass(plan,kScratchAlignment,coordinate_input,"
        "work,use_generic,use_tiled,true);",
        "GatherTransformDecodePass(plan,restored_original,"
        "geometry.aligned_prefix_bytes,geometry.tail_bytes,work,use_generic,"
        "use_tiled,true);",
    )
    if any(token not in core for token in core_tokens):
        raise ModelError(
            "leopard2.cpp decode reveal-fusion policy diverged from schema v3"
        )
    ff8_tokens = (
        "if(ops.kind==LEO2_BACKEND_AVX2||"
        "ops.kind==LEO2_BACKEND_AVX512)returnbuffer_bytes>=1024;",
        "if(ops.kind==LEO2_BACKEND_SSSE3)returnbuffer_bytes>=65536;",
        "IFFT_DIT_DecoderImpl(ops,bytes,m_truncated,work,m,skewLUT,"
        "xor_result,1);",
    )
    if any(token not in ff8 for token in ff8_tokens):
        raise ModelError(
            "LeopardFF8.cpp syndrome-fusion policy diverged from schema v3"
        )
    ff8_sink_assignment = (
        "constbooluse_accumulating_sink="
        "ShouldUsePrunedHighSyndromeSink(ops,buffer_bytes);"
    )
    if ff8.count(ff8_sink_assignment) != 2:
        raise ModelError(
            "LeopardFF8.cpp pruned syndrome-sink wiring diverged from schema v3"
        )
    ff16_tokens = (
        "returnbytes==64||(bytes==128&&"
        "(ops.kind==LEO2_BACKEND_AVX2||"
        "ops.kind==LEO2_BACKEND_AVX512));",
        "return(ops.kind==LEO2_BACKEND_AVX2||"
        "ops.kind==LEO2_BACKEND_AVX512)&&buffer_bytes>=1024;",
        "if(ops.kind==LEO2_BACKEND_SSSE3||"
        "ops.kind==LEO2_BACKEND_AVX2||"
        "ops.kind==LEO2_BACKEND_AVX512){"
        "returnIFFT_DIT_DecoderImpl<true>(",
    )
    if any(token not in ff16 for token in ff16_tokens):
        raise ModelError(
            "LeopardFF16.cpp syndrome-fusion policy diverged from schema v3"
        )
    ff16_sink_assignment = (
        "constbooluse_accumulating_sink="
        "ShouldUsePrunedHighSyndromeSink(ops,buffer_bytes);"
    )
    if ff16.count(ff16_sink_assignment) != 2:
        raise ModelError(
            "LeopardFF16.cpp pruned syndrome-sink wiring diverged from schema v3"
        )


@dataclass(frozen=True)
class Metric:
    value: object
    unit: str
    classification: str
    note: str = ""

    def as_dict(self) -> Dict[str, object]:
        result: Dict[str, object] = {
            "value": self.value,
            "unit": self.unit,
            "classification": self.classification,
        }
        if self.note:
            result["note"] = self.note
        return result


@dataclass
class Schedule:
    butterflies: int = 0
    transform_invocations: int = 0
    radix2_layers: int = 0
    nontransform_xor_vectors: int = 0
    fixed_multiply_vectors: int = 0
    fused_multiply_add_vectors: int = 0
    copies: int = 0
    zero_fills: int = 0
    api_scratch_slots: int = 0
    execution_slots: int = 0
    details: Dict[str, object] = field(default_factory=dict)
    arithmetic_classification: str = "derived_upper_bound"
    memory_override: Optional[Tuple[int, int, int, int]] = None
    execution_bytes_override: Optional[int] = None
    decode_scratch: Optional["DecodeScratchAccounting"] = None
    decode_coordinate_pointer_mappings: Optional[int] = None
    decode_output_gather_payload_bytes: Optional[int] = None
    decode_output_gather_classification: str = "exact_schedule"
    decode_output_gather_note: str = ""
    decode_reveal_temporary_payload_bytes: Optional[int] = None
    decode_reveal_scatter_payload_bytes: Optional[int] = None
    decode_reveal_direct_payload_bytes: Optional[int] = None
    decode_syndrome_fused_payload_bytes: Optional[int] = None
    decode_syndrome_materialized_payload_bytes: Optional[int] = None
    selected_decode_path: Optional[str] = None
    selected_decode_rule: Optional[str] = None
    selected_decode_matching_auto_rules: int = 0

    def add_transform(self, count: int, layers: int) -> None:
        self.butterflies += count
        self.transform_invocations += 1
        self.radix2_layers += layers


@dataclass(frozen=True)
class DecodeScratchAccounting:
    """Exact public decode scratch layout for one declared ABI and path.

    The production layout has three different kinds of storage that must not
    be conflated: AddressRange validation records, coordinate/work pointer
    maps, and shard-data slots.  A ragged final tile additionally reserves one
    fixed 64-byte input slot per public coordinate, although exactly K selected
    coordinates are populated by a valid decode plan.
    """

    workspace: str
    pointer_bytes: int
    range_count: int
    range_metadata_bytes: int
    plan_pointer_count: int
    plan_pointer_map_bytes: int
    plan_metadata_alignment_bytes: int
    plan_data_offset: int
    tail_reserved_slots: int
    tail_reserved_bytes: int
    tail_selected_staged_slots: int
    tail_staged_payload_bytes: int
    tail_staged_zero_padding_bytes: int
    work_slot_bytes: int
    plan_work_slots: int
    plan_work_bytes: int
    plan_total_bytes: int
    codec_pointer_count: int
    codec_range_count: int
    codec_range_metadata_bytes: int
    codec_pointer_map_bytes: int
    codec_metadata_alignment_bytes: int
    codec_data_offset: int
    codec_tail_reserved_slots: int
    codec_tail_reserved_bytes: int
    codec_work_slots: int
    codec_work_bytes: int
    codec_total_bytes: int


def size_t_max(pointer_bytes: int) -> int:
    if pointer_bytes not in (4, 8):
        raise ModelError("pointer bytes must be 4 or 8")
    return (1 << (pointer_bytes * 8)) - 1


def checked_size_add(a: int, b: int, maximum: int) -> int:
    if a < 0 or b < 0 or maximum < 0 or a > maximum or b > maximum - a:
        raise ModelError("decode scratch arithmetic exceeds declared ABI size_t")
    return a + b


def checked_size_multiply(a: int, b: int, maximum: int) -> int:
    if a < 0 or b < 0 or maximum < 0 or a > maximum or b > maximum:
        raise ModelError("decode scratch arithmetic exceeds declared ABI size_t")
    if a != 0 and b > maximum // a:
        raise ModelError("decode scratch arithmetic exceeds declared ABI size_t")
    return a * b


def checked_size_align_up(value: int, alignment: int, maximum: int) -> int:
    if alignment <= 0 or not is_power_of_two(alignment):
        raise ModelError("invalid alignment input")
    summed = checked_size_add(value, alignment - 1, maximum)
    return summed & ~(alignment - 1)


def round_shard_bytes(shard_bytes: int, pointer_bytes: int) -> int:
    """Mirror production RoundShardBytes for a declared size_t width."""
    maximum = size_t_max(pointer_bytes)
    if shard_bytes <= 0 or shard_bytes > UINT64_MAX or shard_bytes > maximum:
        raise ModelError("shard bytes do not fit the declared ABI size_t")
    try:
        return checked_size_align_up(
            shard_bytes, SCRATCH_ALIGNMENT, maximum
        )
    except ModelError as error:
        raise ModelError(
            "shard bytes cannot be rounded within the declared ABI size_t"
        ) from error


def _decode_transform_work_slots(
    parent: int,
    padded: int,
    profile: str,
    output_slots: int,
    workspace: str,
    maximum: int,
) -> int:
    if workspace == "generic" or workspace == "materialized":
        return parent
    if workspace not in ("tiled", "specialized"):
        raise ModelError("decode workspace must be materialized or tiled")
    slots = checked_size_multiply(2, padded, maximum)
    if profile == "high":
        slots = checked_size_add(slots, output_slots, maximum)
    # AUTO specialized execution retains N slots whenever that is no larger
    # than the side-sized form.  A forced tiled diagnostic intentionally does
    # not apply this cap, so callers asking for `tiled` get the forced oracle.
    return min(parent, slots) if workspace == "specialized" else slots


def _scratch_components(
    range_count: int,
    pointer_count: int,
    tail_reserved_slots: int,
    work_slots: int,
    work_slot_bytes: int,
    pointer_bytes: int,
) -> Tuple[int, int, int, int, int, int, int]:
    maximum = size_t_max(pointer_bytes)
    address_range_bytes = checked_size_multiply(2, pointer_bytes, maximum)
    range_bytes = checked_size_multiply(
        range_count, address_range_bytes, maximum
    )
    pointer_offset = checked_size_align_up(
        range_bytes, pointer_bytes, maximum
    )
    pointer_bytes_total = checked_size_multiply(
        pointer_count, pointer_bytes, maximum
    )
    end_pointers = checked_size_add(
        pointer_offset, pointer_bytes_total, maximum
    )
    data_offset = checked_size_align_up(
        end_pointers, SCRATCH_ALIGNMENT, maximum
    )
    alignment_bytes = data_offset - range_bytes - pointer_bytes_total
    tail_bytes = checked_size_multiply(
        tail_reserved_slots, SCRATCH_ALIGNMENT, maximum
    )
    work_data_offset = checked_size_add(data_offset, tail_bytes, maximum)
    work_bytes = checked_size_multiply(work_slots, work_slot_bytes, maximum)
    total = checked_size_add(work_data_offset, work_bytes, maximum)
    return (
        range_bytes,
        pointer_bytes_total,
        alignment_bytes,
        data_offset,
        tail_bytes,
        work_bytes,
        total,
    )


def decode_scratch_accounting(
    k: int,
    r: int,
    parent: int,
    padded: int,
    profile: str,
    shard_bytes: int,
    loss_count: int,
    plan_workspace: str,
    pointer_bytes: int = struct.calcsize("P"),
    codec_workspace: Optional[str] = None,
    no_op: bool = False,
    direct: bool = False,
) -> DecodeScratchAccounting:
    """Mirror DecodeLayout/DirectDecodeLayout and both public queries.

    `plan_workspace` is the path selected by an already-created plan.
    `codec_workspace` is the conservative one-shot query path; for tiled high
    decode it reserves R output slots because no erasure pattern is known yet.
    Forced-path probes pass the same workspace for both.  Direct/no-op plans
    still have an ordinary transform-capable codec query.
    """
    maximum = size_t_max(pointer_bytes)
    if not (0 < k <= UINT32_MAX and 0 < r <= UINT32_MAX):
        raise ModelError("decode counts must fit positive uint32_t values")
    if loss_count < 0 or loss_count > min(k, r):
        raise ModelError("invalid decode loss count")
    if profile not in PROFILES:
        raise ModelError("profile must be high or low")
    expected_geometry = parent_geometry(k, r, profile)
    if (parent != expected_geometry["parent_count"] or
            padded != expected_geometry["padded_side"] or
            parent > 65536 or parent > UINT32_MAX or padded > UINT32_MAX):
        raise ModelError("decode geometry does not match a public codec")
    if no_op and loss_count != 0:
        raise ModelError("a no-op decode plan cannot contain losses")
    if direct and (no_op or loss_count == 0):
        raise ModelError("a direct decode plan requires a loss")

    # Production rejects the public uint64_t input if it cannot first be
    # represented and rounded in size_t, even when a direct layout owns no
    # shard-data slots.  The rounded value is a validation result here; the
    # split layout below deliberately sizes ragged work from the aligned prefix.
    round_shard_bytes(shard_bytes, pointer_bytes)

    range_count = checked_size_add(
        checked_size_multiply(2, k, maximum), r, maximum
    )
    tail = shard_bytes & (SCRATCH_ALIGNMENT - 1)
    prefix = shard_bytes - tail
    work_slot_bytes = shard_bytes if tail == 0 else max(
        prefix, SCRATCH_ALIGNMENT
    )
    tail_reserved_slots = checked_size_add(k, r, maximum) if tail else 0
    tail_selected_staged_slots = k if tail and not no_op and not direct else 0

    if no_op:
        plan_pointer_count = 0
        plan_work_slots = 0
        plan_values = (0, 0, 0, 0, 0, 0, 0)
        plan_range_count = 0
        plan_tail_slots = 0
    elif direct:
        plan_pointer_count = 0
        plan_work_slots = 0
        plan_values = _scratch_components(
            range_count, 0, 0, 0, work_slot_bytes, pointer_bytes
        )
        plan_range_count = range_count
        plan_tail_slots = 0
    else:
        plan_work_slots = _decode_transform_work_slots(
            parent, padded, profile, loss_count, plan_workspace, maximum
        )
        plan_pointer_count = checked_size_add(
            parent, plan_work_slots, maximum
        )
        plan_values = _scratch_components(
            range_count, plan_pointer_count, tail_reserved_slots,
            plan_work_slots, work_slot_bytes, pointer_bytes
        )
        plan_range_count = range_count
        plan_tail_slots = tail_reserved_slots

    if codec_workspace is None:
        codec_workspace = plan_workspace
    codec_work_slots = _decode_transform_work_slots(
        parent, padded, profile,
        r if profile == "high" else loss_count,
        codec_workspace,
        maximum,
    )
    codec_pointer_count = checked_size_add(parent, codec_work_slots, maximum)
    codec_values = _scratch_components(
        range_count, codec_pointer_count, tail_reserved_slots,
        codec_work_slots, work_slot_bytes, pointer_bytes
    )

    (plan_range_bytes, plan_pointer_bytes, plan_alignment, plan_data_offset,
     plan_tail_bytes, plan_work_bytes, plan_total) = plan_values
    (codec_range_bytes, codec_pointer_bytes, codec_alignment,
     codec_data_offset, codec_tail_bytes, codec_work_bytes,
     codec_total) = codec_values
    expected_range_bytes = checked_size_multiply(
        range_count, checked_size_multiply(2, pointer_bytes, maximum), maximum
    )
    if codec_range_bytes != expected_range_bytes:
        raise AssertionError("codec range accounting mismatch")

    tail_payload_bytes = checked_size_multiply(
        tail_selected_staged_slots, tail, maximum
    )
    tail_zero_padding_bytes = checked_size_multiply(
        tail_selected_staged_slots, SCRATCH_ALIGNMENT - tail, maximum
    )

    return DecodeScratchAccounting(
        workspace=("no_op" if no_op else "direct" if direct else plan_workspace),
        pointer_bytes=pointer_bytes,
        range_count=plan_range_count,
        range_metadata_bytes=plan_range_bytes,
        plan_pointer_count=plan_pointer_count,
        plan_pointer_map_bytes=plan_pointer_bytes,
        plan_metadata_alignment_bytes=plan_alignment,
        plan_data_offset=plan_data_offset,
        tail_reserved_slots=plan_tail_slots,
        tail_reserved_bytes=plan_tail_bytes,
        tail_selected_staged_slots=tail_selected_staged_slots,
        tail_staged_payload_bytes=tail_payload_bytes,
        tail_staged_zero_padding_bytes=tail_zero_padding_bytes,
        work_slot_bytes=work_slot_bytes,
        plan_work_slots=plan_work_slots,
        plan_work_bytes=plan_work_bytes,
        plan_total_bytes=plan_total,
        codec_pointer_count=codec_pointer_count,
        codec_range_count=range_count,
        codec_range_metadata_bytes=codec_range_bytes,
        codec_pointer_map_bytes=codec_pointer_bytes,
        codec_metadata_alignment_bytes=codec_alignment,
        codec_data_offset=codec_data_offset,
        codec_tail_reserved_slots=tail_reserved_slots,
        codec_tail_reserved_bytes=codec_tail_bytes,
        codec_work_slots=codec_work_slots,
        codec_work_bytes=codec_work_bytes,
        codec_total_bytes=codec_total,
    )


def is_power_of_two(value: int) -> bool:
    return value > 0 and (value & (value - 1)) == 0


def ceil_power_of_two(value: int) -> int:
    if value <= 0:
        raise ModelError("counts must be positive")
    return 1 << (value - 1).bit_length()


def ceil_div(numerator: int, denominator: int) -> int:
    return (numerator + denominator - 1) // denominator


def rounded_kernel_bytes(shard_bytes: int) -> int:
    # Operation reports model a uint64_t public shard length.  Decode scratch
    # accounting applies the stricter declared-ABI size_t limit separately.
    return round_shard_bytes(shard_bytes, 8)


def parent_geometry(k: int, r: int, profile: str) -> Dict[str, int]:
    if k <= 0 or r <= 0:
        raise ModelError("K and R must be positive")
    if profile == "high":
        padded = ceil_power_of_two(r)
        parent = ceil_power_of_two(k + padded)
        dimension = parent - padded
        shortened = dimension - k
        punctured = padded - r
    elif profile == "low":
        padded = ceil_power_of_two(k)
        parent = ceil_power_of_two(padded + r)
        dimension = padded
        shortened = padded - k
        punctured = parent - padded - r
    else:
        raise ModelError("profile must be high or low")
    return {
        "padded_side": padded,
        "parent_count": parent,
        "parent_dimension": dimension,
        "shortened_count": shortened,
        "punctured_count": punctured,
    }


def validate_field(parent: int, field_name: str) -> None:
    if field_name not in FIELDS:
        raise ModelError("field must be gf8 or gf16")
    order = 256 if field_name == "gf8" else 65536
    if parent > order:
        raise ModelError(
            "parent length {} exceeds {} order {}".format(parent, field_name, order)
        )


def parse_mask(spec: Optional[str], limit: int, default_all: bool) -> Set[int]:
    if limit < 0:
        raise ModelError("mask limit cannot be negative")
    if spec is None:
        return set(range(limit)) if default_all else set()
    text = spec.strip().lower()
    if text in ("", "none"):
        return set()
    if text == "all":
        return set(range(limit))
    result: Set[int] = set()
    for token in text.split(","):
        token = token.strip()
        if not token:
            raise ModelError("empty mask token")
        if "-" in token:
            parts = token.split("-")
            if len(parts) != 2 or not parts[0] or not parts[1]:
                raise ModelError("invalid mask range: {}".format(token))
            begin = int(parts[0], 10)
            end = int(parts[1], 10)
            if begin > end:
                raise ModelError("descending mask range: {}".format(token))
            result.update(range(begin, end + 1))
        else:
            result.add(int(token, 10))
    if any(index < 0 or index >= limit for index in result):
        raise ModelError("mask index outside [0, {})".format(limit))
    return result


def mask_to_ranges(indices: Iterable[int]) -> str:
    ordered = sorted(set(indices))
    if not ordered:
        return "none"
    ranges: List[str] = []
    begin = ordered[0]
    end = begin
    for value in ordered[1:]:
        if value == end + 1:
            end = value
            continue
        ranges.append(str(begin) if begin == end else "{}-{}".format(begin, end))
        begin = end = value
    ranges.append(str(begin) if begin == end else "{}-{}".format(begin, end))
    return ",".join(ranges)


def _transform_layers(size: int) -> int:
    if not is_power_of_two(size):
        raise ModelError("transform size must be a power of two")
    return size.bit_length() - 1


def ifft_prefix_butterflies(size: int, active_prefix: int) -> int:
    """Exact current DIT-IFFT schedule, expressed as radix-2 butterflies."""
    layers = _transform_layers(size)
    if active_prefix < 0 or active_prefix > size:
        raise ModelError("active prefix outside transform")
    if size == 1:
        return 0
    count = 0
    dist = 1
    dist4 = 4
    while dist4 <= size:
        count += ceil_div(active_prefix, dist4) * dist4
        dist = dist4
        dist4 <<= 2
    # This final layer is unconditional in IFFT_DIT_Decoder once present.
    # Encoder calls always have a non-empty prefix; retaining the production
    # behavior for an empty decoder block makes this a schedule model, not an
    # idealized dependency count.
    if dist < size:
        count += size // 2
    expected_full = size * layers // 2
    if active_prefix == size and count != expected_full:
        raise AssertionError("internal IFFT full-schedule mismatch")
    return count


def fft_prefix_butterflies(size: int, requested_prefix: int) -> int:
    """Exact current prefix-pruned DIT-FFT schedule in radix-2 butterflies."""
    layers = _transform_layers(size)
    if requested_prefix < 0 or requested_prefix > size:
        raise ModelError("requested prefix outside transform")
    if size == 1 or requested_prefix == 0:
        return 0
    count = 0
    dist4 = size
    dist = size >> 2
    while dist != 0:
        count += ceil_div(requested_prefix, dist4) * dist4
        dist4 = dist
        dist >>= 2
    if dist4 == 2:
        count += ceil_div(requested_prefix, 2)
    expected_full = size * layers // 2
    if requested_prefix == size and count != expected_full:
        raise AssertionError("internal FFT full-schedule mismatch")
    return count


def fft_mask_butterflies(size: int, requested: Set[int]) -> int:
    """Exact ErrorBitfield block schedule for an arbitrary output mask."""
    layers = _transform_layers(size)
    if any(index < 0 or index >= size for index in requested):
        raise ModelError("requested output outside transform")
    if size == 1 or not requested:
        return 0
    count = 0
    dist4 = size
    dist = size >> 2
    while dist != 0:
        blocks = {index // dist4 for index in requested}
        count += len(blocks) * dist4
        dist4 = dist
        dist >>= 2
    if dist4 == 2:
        count += len({index // 2 for index in requested})
    expected_full = size * layers // 2
    if len(requested) == size and count != expected_full:
        raise AssertionError("internal masked FFT full-schedule mismatch")
    return count


def derivative_xor_vectors(size: int) -> int:
    """Exact VectorXOR shard count in Leopard's normalized derivative loop."""
    if not is_power_of_two(size):
        raise ModelError("derivative size must be a power of two")
    return sum(index & -index for index in range(1, size))


def _profile_for_path(path: str, requested_profile: Optional[str]) -> str:
    if path.startswith("legacy_high"):
        implied = "high"
    elif path.startswith("low_"):
        implied = "low"
    else:
        implied = requested_profile or "high"
    if requested_profile and requested_profile != implied:
        raise ModelError("path {} requires profile {}".format(path, implied))
    return implied


def _coordinate_for_original(profile: str, padded: int, original: int) -> int:
    return padded + original if profile == "high" else original


def _coordinate_for_recovery(profile: str, padded: int, recovery: int) -> int:
    return recovery if profile == "high" else padded + recovery


def deterministic_decode_coordinates(
    k: int, r: int, profile: str, padded: int, losses: Set[int]
) -> Tuple[Set[int], Set[int]]:
    if len(losses) > r:
        raise ModelError("loss count exceeds public recovery count")
    data = {
        _coordinate_for_original(profile, padded, original)
        for original in range(k)
        if original not in losses
    }
    selected_parities = set(range(len(losses)))
    data.update(
        _coordinate_for_recovery(profile, padded, parity)
        for parity in selected_parities
    )
    if len(data) != k:
        raise AssertionError("deterministic received subset is not K coordinates")
    return data, selected_parities


def _block_prefix(data_coordinates: Set[int], offset: int, width: int) -> int:
    local = [coordinate - offset for coordinate in data_coordinates
             if offset <= coordinate < offset + width]
    return max(local) + 1 if local else 0


@dataclass(frozen=True)
class DecodeSelection:
    path: str
    rule: str
    matching_auto_rules: int
    required_work_slots: int


def select_decode_execution(
    k: int,
    r: int,
    profile: str,
    field_name: str,
    backend: str,
    geometry: Dict[str, int],
    shard_bytes: int,
    loss_count: int,
    force: str = "auto",
    multi_item_batch: bool = False,
) -> DecodeSelection:
    """Independent mirror of terminal-plan and transform path selection."""
    if backend not in BACKENDS:
        raise ModelError("unsupported backend")
    if force not in ("auto", "generic", "specialized", "materialized", "tiled"):
        raise ModelError("unsupported decode selection")
    if loss_count < 0 or loss_count > min(k, r):
        raise ModelError("invalid decode loss count")
    padded = geometry["padded_side"]
    parent = geometry["parent_count"]
    # Production terminal plan introspection and execution return no-op before
    # constructing byte geometry.  Public scratch queries are stricter; the
    # report's later scratch accounting still validates their shard length.
    if loss_count == 0:
        return DecodeSelection("no_op", "no_op", 0, 0)
    rounded = rounded_kernel_bytes(shard_bytes)
    # R=1 legacy-high and K=1 low profiles are unconditional terminal plans,
    # including codecs with forced transform flags.
    if padded == 1:
        return DecodeSelection("direct", "direct", 0, 0)
    # Ordinary AUTO prepares the bounded LxL direct plan before transform
    # metadata.  Any explicit transform selection disables this path.
    if (force == "auto" and 2 <= k <= 16 and
            geometry["parent_dimension"] <= 256 and loss_count <= 4):
        return DecodeSelection("direct", "direct", 0, 0)

    tiled_slots = 2 * padded + (loss_count if profile == "high" else 0)
    materialized_slots = parent
    balanced = (
        profile == "high" and field_name == "gf8" and
        k == 128 and r == 128 and padded == 128 and parent == 256 and
        loss_count == 128 and 256 <= rounded <= 1024 * 1024 and
        backend in ("scalar", "ssse3", "avx2", "avx512")
    )
    measured_materialized = (
        profile == "high" and field_name == "gf8" and
        k == 224 and r == 32 and padded == 32 and parent == 256 and
        0 < loss_count <= 8 and rounded <= 64 * 1024 and
        ((backend in ("avx2", "avx512") and rounded >= 24 * 1024) or
         (backend == "ssse3" and rounded >= 32 * 1024))
    )
    matching = (1 if balanced else 0) | (2 if measured_materialized else 0)
    if force == "generic":
        return DecodeSelection(
            "generic", "forced_generic", matching, materialized_slots
        )
    if force == "materialized":
        return DecodeSelection(
            "materialized", "forced_materialized", matching,
            materialized_slots,
        )
    if force == "tiled":
        return DecodeSelection("tiled", "forced_tiled", matching, tiled_slots)

    if force != "specialized" and balanced:
        return DecodeSelection(
            "generic", "balanced_generic", matching, materialized_slots
        )
    if measured_materialized:
        if multi_item_batch and backend in ("avx2", "avx512"):
            return DecodeSelection(
                "tiled", "measured_batch_tiled", matching, tiled_slots
            )
        return DecodeSelection(
            "materialized", "measured_materialized", matching,
            materialized_slots,
        )
    if tiled_slots < materialized_slots:
        return DecodeSelection("tiled", "workspace_tiled", matching, tiled_slots)
    return DecodeSelection(
        "materialized", "workspace_materialized", matching,
        materialized_slots,
    )


def select_codec_transform_execution(
    k: int,
    r: int,
    profile: str,
    field_name: str,
    backend: str,
    geometry: Dict[str, int],
    shard_bytes: int,
    force: str = "auto",
) -> DecodeSelection:
    """Mirror the pattern-independent one-shot codec scratch query."""
    padded = geometry["padded_side"]
    parent = geometry["parent_count"]
    rounded = rounded_kernel_bytes(shard_bytes)
    tiled_slots = 2 * padded + (r if profile == "high" else 0)
    if force not in ("auto", "specialized"):
        if force not in ("generic", "materialized", "tiled"):
            raise ModelError("unsupported codec decode selection")
    measured_materialized = (
        profile == "high" and field_name == "gf8" and
        k == 224 and r == 32 and padded == 32 and parent == 256 and
        rounded <= 64 * 1024 and
        ((backend in ("avx2", "avx512") and rounded >= 24 * 1024) or
         (backend == "ssse3" and rounded >= 32 * 1024))
    )
    matching = 2 if measured_materialized else 0
    if force == "generic":
        return DecodeSelection("generic", "forced_generic", matching, parent)
    if force == "materialized":
        return DecodeSelection(
            "materialized", "forced_materialized", matching, parent
        )
    if force == "tiled":
        return DecodeSelection("tiled", "forced_tiled", matching, tiled_slots)
    if measured_materialized:
        return DecodeSelection(
            "materialized", "measured_materialized", matching, parent
        )
    if tiled_slots < parent:
        return DecodeSelection("tiled", "workspace_tiled", 0, tiled_slots)
    return DecodeSelection("materialized", "workspace_materialized", 0, parent)


def _gf16_mature_syndrome_fuses(
    padded: int, backend: str, pass_bytes: int
) -> bool:
    if backend not in ("ssse3", "avx2", "avx512"):
        return False
    # An odd transform depth ends in a two-way layer, which the qualified
    # vector decoder always accumulates.  Even depths end in a four-way layer;
    # compact fused-four kernels retain the measured materialize-then-XOR path.
    if _transform_layers(padded) & 1:
        return True
    if pass_bytes == 64:
        return False
    if pass_bytes == 128 and backend in ("avx2", "avx512"):
        return False
    return True


def _pruned_syndrome_fuses(
    field_name: str, backend: str, pass_bytes: int
) -> bool:
    if backend in ("avx2", "avx512"):
        return pass_bytes >= 1024
    return field_name == "gf8" and backend == "ssse3" and pass_bytes >= 65536


def apply_decode_backend_accounting(
    schedule: Schedule,
    k: int,
    r: int,
    profile: str,
    field_name: str,
    backend: str,
    geometry: Dict[str, int],
    shard_bytes: int,
    losses: Set[int],
    selection: DecodeSelection,
    multi_item_batch: bool,
) -> None:
    """Attach exact selected-rule and backend-qualified boundary traffic."""
    schedule.selected_decode_path = selection.path
    schedule.selected_decode_rule = selection.rule
    schedule.selected_decode_matching_auto_rules = selection.matching_auto_rules
    schedule.details.update({
        "selected_decode_path": selection.path,
        "selected_decode_rule": selection.rule,
        "selected_decode_matching_auto_rules": selection.matching_auto_rules,
        "selected_backend": backend,
        "multi_item_batch": multi_item_batch,
    })
    if selection.path in ("no_op", "direct"):
        return

    loss_count = len(losses)
    tail = shard_bytes & (SCRATCH_ALIGNMENT - 1)
    aligned = shard_bytes - tail
    kernel_bytes = rounded_kernel_bytes(shard_bytes)
    use_generic = selection.path == "generic"
    use_low_specialized = not use_generic and profile == "low"
    fuse_aligned_reveal = (
        aligned != 0 and
        (use_low_specialized or
         (use_generic and field_name == "gf8" and aligned >= 4096 and
          backend in ("ssse3", "avx2", "avx512")))
    )

    if use_generic or use_low_specialized:
        unfused_kernel = (0 if fuse_aligned_reveal else aligned) + (64 if tail else 0)
        scatter = (0 if fuse_aligned_reveal else aligned) + tail
        direct = aligned if fuse_aligned_reveal else 0
        temporary = unfused_kernel
    elif selection.path == "materialized":
        # Algorithm 5 reveals in-place in the materialized output block.
        temporary = kernel_bytes
        scatter = shard_bytes
        direct = 0
    else:
        # Tiled Algorithm 5 multiplies out of place into a compact requested-
        # output workspace, then scatters the public payload.
        temporary = 0
        scatter = shard_bytes
        direct = 0
    schedule.decode_reveal_temporary_payload_bytes = temporary * loss_count
    schedule.decode_reveal_scatter_payload_bytes = scatter * loss_count
    schedule.decode_reveal_direct_payload_bytes = direct * loss_count
    schedule.decode_output_gather_payload_bytes = scatter * loss_count
    schedule.decode_output_gather_classification = "exact_schedule"
    schedule.decode_output_gather_note = (
        "Backend-qualified execution: direct reveal owns the aligned prefix; "
        "only the reported residual payload is copied from scratch."
        if fuse_aligned_reveal else
        "Selected execution retains scratch reveal and the reported scatter."
    )

    fused_syndrome = 0
    materialized_syndrome = 0
    if profile == "high" and not use_generic:
        padded = geometry["padded_side"]
        parent = geometry["parent_count"]
        data, _ = deterministic_decode_coordinates(k, r, profile, padded, losses)
        pass_sizes = ([aligned] if aligned else []) + ([64] if tail else [])
        for offset in range(padded, parent, padded):
            live_count = sum(
                1 for coordinate in data
                if offset <= coordinate < offset + padded
            )
            if live_count == 0:
                continue
            pruned = live_count != padded
            for pass_bytes in pass_sizes:
                if pruned:
                    fused = _pruned_syndrome_fuses(
                        field_name, backend, pass_bytes
                    )
                elif field_name == "gf8":
                    # The mature GF8 decoder always folds the final inverse
                    # layer into the accumulator; backend qualification affects
                    # source staging, not this boundary.
                    fused = True
                else:
                    fused = _gf16_mature_syndrome_fuses(
                        padded, backend, pass_bytes
                    )
                payload = padded * pass_bytes
                if fused:
                    fused_syndrome += payload
                else:
                    materialized_syndrome += payload
    schedule.decode_syndrome_fused_payload_bytes = fused_syndrome
    schedule.decode_syndrome_materialized_payload_bytes = materialized_syndrome
    schedule.details.update({
        "reveal_aligned_fused": fuse_aligned_reveal,
        "reveal_inplace_temporary_payload_bytes":
            schedule.decode_reveal_temporary_payload_bytes,
        "reveal_scatter_payload_bytes": schedule.decode_reveal_scatter_payload_bytes,
        "reveal_direct_payload_bytes": schedule.decode_reveal_direct_payload_bytes,
        "syndrome_fused_accumulation_payload_bytes": fused_syndrome,
        "syndrome_materialized_xor_payload_bytes": materialized_syndrome,
        "backend_traffic_classification": "exact_production_predicates",
    })


def model_high_encode(
    k: int, r: int, padded: int, shard_bytes: int, requested: Set[int]
) -> Schedule:
    schedule = Schedule()
    partial = rounded_kernel_bytes(shard_bytes) != shard_bytes
    schedule.api_scratch_slots = 2 * padded + (k if partial else 0)
    schedule.execution_slots = 2 * padded
    schedule.details["requested_parity"] = mask_to_ranges(requested)
    if not requested:
        schedule.api_scratch_slots = 2 * padded + (k if partial else 0)
        schedule.execution_slots = 0
        schedule.details["no_op"] = True
        return schedule
    requested_prefix = max(requested) + 1
    if padded == 1:
        schedule.nontransform_xor_vectors = k - 1
        schedule.copies = 1 + (k if partial else 0) + (len(requested) if partial else 0)
        schedule.zero_fills = 0
        schedule.details["direct_r1_xor"] = True
        return schedule
    block_count = ceil_div(k, padded)
    remaining = k
    for _ in range(block_count):
        active = min(remaining, padded)
        schedule.copies += active
        schedule.zero_fills += padded - active
        schedule.add_transform(
            ifft_prefix_butterflies(padded, active), _transform_layers(padded)
        )
        remaining -= active
    schedule.nontransform_xor_vectors += (block_count - 1) * padded
    schedule.add_transform(
        fft_prefix_butterflies(padded, requested_prefix),
        _transform_layers(padded),
    )
    if partial:
        schedule.copies += k + len(requested)
    schedule.details.update({
        "message_block_count": block_count,
        "requested_prefix": requested_prefix,
        "fused_accumulation_xor_vectors": (block_count - 1) * padded,
    })
    return schedule


def model_low_encode(
    k: int, r: int, padded: int, shard_bytes: int, requested: Set[int]
) -> Schedule:
    schedule = Schedule()
    partial = rounded_kernel_bytes(shard_bytes) != shard_bytes
    schedule.api_scratch_slots = 2 * padded + (k + r if partial else 0)
    schedule.execution_slots = 2 * padded
    schedule.details["requested_parity"] = mask_to_ranges(requested)
    if not requested:
        schedule.execution_slots = 0
        schedule.details["no_op"] = True
        return schedule
    if padded == 1:
        schedule.copies = len(requested)
        schedule.details["direct_k1_copy"] = True
        return schedule
    schedule.copies += k
    schedule.zero_fills += padded - k
    schedule.add_transform(
        ifft_prefix_butterflies(padded, k), _transform_layers(padded)
    )
    active_blocks = 0
    direct_output_blocks = 0
    block_prefixes: List[int] = []
    for offset in range(0, r, padded):
        block_mask = {index - offset for index in requested
                      if offset <= index < min(offset + padded, r)}
        if not block_mask:
            continue
        prefix = max(block_mask) + 1
        active_blocks += 1
        block_prefixes.append(prefix)
        # Production backends load the immutable coefficient block directly
        # into an out-of-place first FFT layer.  The transform traffic below
        # therefore accounts for those reads and evaluation-work writes; there
        # is no separate padded-shard copy to charge here.
        schedule.add_transform(
            fft_prefix_butterflies(padded, prefix), _transform_layers(padded)
        )
        # A complete dense parity block executes its final forward layer
        # out-of-place into the disjoint caller destinations.  Sparse and
        # partial blocks retain scratch evaluation plus the public scatter.
        if len(block_mask) == padded:
            direct_output_blocks += 1
        else:
            schedule.copies += len(block_mask)
    if partial:
        schedule.copies += k + len(requested)
    schedule.details.update({
        "active_parity_blocks": active_blocks,
        "direct_output_blocks": direct_output_blocks,
        "parity_block_prefixes": block_prefixes,
        "out_of_place_first_fft_layer": True,
        "out_of_place_final_fft_layer_for_dense_blocks": True,
    })
    return schedule


def _decode_base(
    k: int,
    r: int,
    parent: int,
    padded: int,
    profile: str,
    shard_bytes: int,
    losses: Set[int],
) -> Tuple[Schedule, Set[int], Set[int]]:
    schedule = Schedule()
    data, selected_parities = deterministic_decode_coordinates(
        k, r, profile, padded, losses
    )
    # The N coordinate entries are a pointer map, not shard-data storage.
    # Complete tiles retain the selected K caller pointers directly.  Output
    # gathers can be shorter than one rounded kernel vector, so decode paths
    # report those boundary payload bytes separately.
    schedule.decode_coordinate_pointer_mappings = len(data)
    schedule.details.update({
        "missing_originals": mask_to_ranges(losses),
        "selected_parity": mask_to_ranges(selected_parities),
        "received_coordinate_count": len(data),
        "coordinate_pointer_mappings": len(data),
        "aligned_input_staging_copies": 0,
        "kernel_rounding_bytes": rounded_kernel_bytes(shard_bytes) - shard_bytes,
    })
    return schedule, data, selected_parities


def model_generic_decode(
    k: int,
    r: int,
    parent: int,
    padded: int,
    profile: str,
    shard_bytes: int,
    losses: Set[int],
) -> Schedule:
    if not losses:
        result = Schedule(details={"no_op": True, "missing_originals": "none"})
        return result
    schedule, data, _ = _decode_base(
        k, r, parent, padded, profile, shard_bytes, losses
    )
    requested = {
        _coordinate_for_original(profile, padded, original) for original in losses
    }
    schedule.fixed_multiply_vectors += k
    schedule.zero_fills += parent - k
    input_prefix = max(data) + 1
    schedule.add_transform(
        ifft_prefix_butterflies(parent, input_prefix), _transform_layers(parent)
    )
    derivative = derivative_xor_vectors(parent)
    schedule.nontransform_xor_vectors += derivative
    schedule.add_transform(
        fft_mask_butterflies(parent, requested), _transform_layers(parent)
    )
    schedule.fixed_multiply_vectors += len(losses)
    schedule.decode_output_gather_payload_bytes = len(losses) * shard_bytes
    schedule.decode_output_gather_classification = "estimated_upper_bound"
    schedule.decode_output_gather_note = (
        "Copy-first backend-neutral baseline. Qualified GF8 "
        "SSSE3/AVX2/AVX512 "
        "generic reveal fusion can remove the aligned-prefix gather."
    )
    schedule.details.update({
        "input_prefix": input_prefix,
        "requested_coordinates": mask_to_ranges(requested),
        "derivative_xor_vectors": derivative,
    })
    return schedule


def model_low_decode(
    k: int,
    r: int,
    parent: int,
    padded: int,
    shard_bytes: int,
    losses: Set[int],
) -> Schedule:
    if not losses:
        return Schedule(details={"no_op": True, "missing_originals": "none"})
    if padded == 1:
        return Schedule(
            decode_output_gather_payload_bytes=shard_bytes,
            details={"direct_k1_copy": True, "missing_originals": "0"},
        )
    schedule, data, _ = _decode_base(
        k, r, parent, padded, "low", shard_bytes, losses
    )
    schedule.fixed_multiply_vectors += k
    schedule.zero_fills += parent - k
    prefixes: List[int] = []
    for offset in range(0, parent, padded):
        prefix = _block_prefix(data, offset, padded)
        prefixes.append(prefix)
        schedule.add_transform(
            ifft_prefix_butterflies(padded, prefix), _transform_layers(padded)
        )
    derivative = derivative_xor_vectors(padded)
    schedule.nontransform_xor_vectors += derivative
    reduction = (parent // padded - 1) * padded
    schedule.fixed_multiply_vectors += reduction
    schedule.nontransform_xor_vectors += reduction
    schedule.fused_multiply_add_vectors += reduction
    schedule.add_transform(
        fft_mask_butterflies(padded, set(losses)), _transform_layers(padded)
    )
    schedule.fixed_multiply_vectors += len(losses)
    schedule.decode_output_gather_payload_bytes = (
        shard_bytes & 63
    ) * len(losses)
    schedule.decode_output_gather_note = (
        "Algorithm 4 reveals aligned output directly; only a ragged tail is gathered."
    )
    schedule.details.update({
        "block_input_prefixes": prefixes,
        "derivative_xor_vectors": derivative,
        "weighted_reduction_vectors": reduction,
        "fused_weighted_multiply_add_vectors": reduction,
    })
    return schedule


def model_high_decode(
    k: int,
    r: int,
    parent: int,
    padded: int,
    shard_bytes: int,
    losses: Set[int],
    field_name: str = "gf8",
) -> Schedule:
    if not losses:
        return Schedule(details={"no_op": True, "missing_originals": "none"})
    if padded == 1:
        if len(losses) != 1:
            raise ModelError("R=1 cannot recover more than one original")
        return Schedule(
            nontransform_xor_vectors=k - 1,
            decode_output_gather_payload_bytes=shard_bytes,
            details={"direct_r1_xor": True, "missing_originals": mask_to_ranges(losses)},
        )
    schedule, data, selected_parities = _decode_base(
        k, r, parent, padded, "high", shard_bytes, losses
    )
    # The mature Algorithm 5 input schedule consumes complete live blocks in
    # four-row groups through the first two inverse layers.  A partial block
    # owns an exact pruned input plan and is staged as a whole before that plan
    # executes; source fusion must not be credited inside it.  An empty later
    # block is skipped.  This is an execution-boundary change: butterfly
    # counts are identical to the former stage-then-transform schedule.
    receive_copy_vectors = 0
    receive_zero_vectors = 0
    receive_fused_groups = 0
    receive_exact_pruned_blocks = 0
    receive_skipped_blocks = 0
    active_later_blocks = 0
    prefixes: List[int] = []
    for offset in range(0, parent, padded):
        prefix = _block_prefix(data, offset, padded)
        prefixes.append(prefix)
        live = {coordinate - offset for coordinate in data
                if offset <= coordinate < offset + padded}
        if not live and offset != 0:
            receive_skipped_blocks += 1
            continue
        if offset != 0:
            active_later_blocks += 1
        if padded <= 4 or field_name != "gf8":
            receive_copy_vectors += len(live)
            receive_zero_vectors += padded - len(live)
        elif len(live) != padded:
            receive_copy_vectors += len(live)
            receive_zero_vectors += padded - len(live)
            if live:
                receive_exact_pruned_blocks += 1
        else:
            receive_fused_groups += padded // 4
        schedule.add_transform(
            ifft_prefix_butterflies(padded, prefix), _transform_layers(padded)
        )
    # Keep the ISA-independent schedule totals on the original copy-first
    # contract.  The source-boundary values below are a qualified
    # SSSE3/AVX2/AVX512 delta, not a backend-neutral replacement for those
    # totals.
    # A future backend-aware model can apply the delta to absolute traffic.
    schedule.copies += k
    schedule.zero_fills += parent - k
    reduction = active_later_blocks * padded
    schedule.nontransform_xor_vectors += reduction
    schedule.add_transform(
        fft_prefix_butterflies(padded, padded), _transform_layers(padded)
    )
    schedule.fixed_multiply_vectors += len(selected_parities)
    schedule.add_transform(
        ifft_prefix_butterflies(padded, padded), _transform_layers(padded)
    )

    evaluation_prefixes: List[Tuple[int, int]] = []
    requested_coordinates = {
        _coordinate_for_original("high", padded, original) for original in losses
    }
    for offset in range(padded, parent, padded):
        local = {coordinate - offset for coordinate in requested_coordinates
                 if offset <= coordinate < offset + padded}
        if not local:
            continue
        prefix = max(local) + 1
        evaluation_prefixes.append((offset // padded, prefix))
        schedule.add_transform(
            fft_prefix_butterflies(padded, prefix), _transform_layers(padded)
        )
    schedule.fixed_multiply_vectors += len(losses)
    schedule.decode_output_gather_payload_bytes = len(losses) * shard_bytes
    schedule.decode_output_gather_note = (
        "Algorithm 5 gathers each requested original from its evaluator workspace."
    )
    schedule.details.update({
        "block_input_prefixes": prefixes,
        "receive_source_fused_radix4_groups": receive_fused_groups,
        "receive_copy_vectors": receive_copy_vectors,
        "receive_zero_vectors": receive_zero_vectors,
        "receive_copy_vectors_removed": 4 * receive_fused_groups,
        "receive_exact_pruned_staged_blocks": receive_exact_pruned_blocks,
        "receive_skipped_empty_blocks": receive_skipped_blocks,
        "active_later_syndrome_blocks": active_later_blocks,
        "receive_source_fusion_scope": (
            "GF8 qualified SSSE3/AVX2/AVX512 mature unpruned schedule delta"
            if field_name == "gf8" else
            "GF16 deterministic copy-first policy"
        ),
        "syndrome_reduction_vectors": reduction,
        "locator_scale_vectors": len(selected_parities),
        "locator_weighted_fusion_scope": (
            "qualified GF8 AVX2/AVX512 delta only: T>=64, "
            "live_count>=ceil(T/2), "
            "and each kernel pass in [16 KiB, 256 KiB]"
        ),
        "locator_weighted_live_scan_scope": (
            "statically qualified GF8 AVX2/AVX512 passes scan receive pointers "
            "until ceil(T/2) live rows are found; sparse fallback scans all T"
        ),
        "locator_weighted_fusion_applied_to_isa_independent_totals": False,
        "evaluation_block_prefixes": evaluation_prefixes,
        "out_of_place_evaluator_first_layer": True,
        "requested_coordinates": mask_to_ranges(requested_coordinates),
    })
    return schedule


def model_direct_repair(
    k: int,
    r: int,
    parent_dimension: int,
    padded: int,
    profile: str,
    losses: Set[int],
    shard_bytes: Optional[int] = None,
) -> Schedule:
    loss_count = len(losses)
    if loss_count == 0:
        return Schedule(details={"no_op": True, "missing_originals": "none"})
    if not (2 <= k <= 16):
        raise ModelError("bounded direct repair requires 2 <= K <= 16")
    if loss_count > 4:
        raise ModelError("bounded direct repair supports at most four losses")
    if loss_count > r:
        raise ModelError("loss count exceeds public recovery count")
    if parent_dimension > 256 or padded < 2:
        raise ModelError("bounded direct repair parent is outside dispatch limits")
    # The execution term list is coefficient-dependent.  K dense terms per
    # recovered output is a safe bound: L selected parity equations plus all
    # K-L surviving originals.  Cancellation can make the actual list shorter,
    # and a coefficient equal to one removes its multiplication.
    dense_terms = loss_count * k
    additions = loss_count * (k - 1)
    reads = loss_count * (1 + 2 * (k - 1))
    writes = dense_terms
    schedule = Schedule(
        nontransform_xor_vectors=additions,
        fixed_multiply_vectors=dense_terms,
        api_scratch_slots=0,
        execution_slots=0,
        arithmetic_classification="estimated_upper_bound",
        # Both memory bounds are the dense fused multiply/add traversal.  It is
        # an upper bound, not an exact term-list count.
        memory_override=(reads, writes, reads, writes),
        execution_bytes_override=shard_bytes,
        details={
            "missing_originals": mask_to_ranges(losses),
            "dense_term_upper_bound": dense_terms,
            "coefficient_dependent_terms": True,
            "profile": profile,
        },
    )
    return schedule


def model_direct_terminal(
    k: int, profile: str, losses: Set[int], shard_bytes: int
) -> Schedule:
    """Model the unconditional R=1 XOR or K=1 copy terminal plan."""
    if len(losses) != 1:
        raise ModelError("terminal direct path requires exactly one loss")
    schedule = Schedule(
        execution_slots=0,
        api_scratch_slots=0,
        execution_bytes_override=shard_bytes,
    )
    schedule.copies = 1
    if profile == "high":
        # Copy the XOR parity, then XOR every surviving original into it.
        schedule.nontransform_xor_vectors = k - 1
        schedule.details["direct_mode"] = "r1_xor"
    else:
        if k != 1:
            raise ModelError("low terminal direct-copy path requires K=1")
        schedule.details["direct_mode"] = "k1_copy"
    schedule.details.update({
        "missing_originals": mask_to_ranges(losses),
        "plan_setup": "excluded",
    })
    return schedule


def model_generic_transform(
    parent: int,
    direction: str,
    active_input_count: Optional[int],
    requested: Set[int],
) -> Schedule:
    schedule = Schedule(api_scratch_slots=parent, execution_slots=parent)
    layers = _transform_layers(parent)
    if direction == "inverse":
        active = parent if active_input_count is None else active_input_count
        schedule.add_transform(ifft_prefix_butterflies(parent, active), layers)
        schedule.details["active_input_prefix"] = active
    elif direction == "forward":
        schedule.add_transform(fft_mask_butterflies(parent, requested), layers)
        schedule.details["requested_outputs"] = mask_to_ranges(requested)
    else:
        raise ModelError("direction must be forward or inverse")
    schedule.details["direction"] = direction
    return schedule


def schedule_metrics(
    schedule: Schedule, field_name: str, shard_bytes: int
) -> Dict[str, Metric]:
    kernel_bytes = (
        schedule.execution_bytes_override
        if schedule.execution_bytes_override is not None
        else rounded_kernel_bytes(shard_bytes)
    )
    symbol_bytes = 1 if field_name == "gf8" else 2
    symbols = kernel_bytes // symbol_bytes
    butterfly_xors_lower = schedule.butterflies
    butterfly_xors_upper = 2 * schedule.butterflies
    xor_vectors_lower = butterfly_xors_lower + schedule.nontransform_xor_vectors
    xor_vectors_upper = butterfly_xors_upper + schedule.nontransform_xor_vectors
    multiply_vectors_upper = schedule.butterflies + schedule.fixed_multiply_vectors

    if schedule.memory_override is None:
        # One fused butterfly reads and writes both operands.  A scalar-style
        # two-primitive implementation can reread both operands for the second
        # update.  These bounds exclude cache-line fetch amplification.
        reads_lower = (
            2 * schedule.butterflies
            + 2 * schedule.nontransform_xor_vectors
            + schedule.fixed_multiply_vectors
            + schedule.copies
            - schedule.fused_multiply_add_vectors
        )
        writes_lower = (
            2 * schedule.butterflies
            + schedule.nontransform_xor_vectors
            + schedule.fixed_multiply_vectors
            + schedule.copies
            + schedule.zero_fills
            - schedule.fused_multiply_add_vectors
        )
        reads_upper = (
            4 * schedule.butterflies
            + 2 * schedule.nontransform_xor_vectors
            + schedule.fixed_multiply_vectors
            + schedule.copies
            - schedule.fused_multiply_add_vectors
        )
        writes_upper = writes_lower
    else:
        reads_lower, writes_lower, reads_upper, writes_upper = schedule.memory_override

    arithmetic_note = (
        "Transform butterflies with zero skew execute one XOR and no fixed "
        "multiply; nonzero skew executes two XORs and one multiply."
    )
    memory_note = (
        "Logical full-kernel shard traffic only; excludes decode ragged/output "
        "boundary payload metrics, metadata, allocator effects, cache-line "
        "amplification, and backend table loads."
    )
    full_pass_denominator = max(schedule.execution_slots, 1)
    metrics = {
        "transform_butterflies": Metric(
            schedule.butterflies, "radix2_butterfly_equivalents",
            "exact_schedule", "Fused radix-4 work expanded into radix-2 equivalents."
        ),
        "transform_invocations": Metric(
            schedule.transform_invocations, "calls", "exact_schedule"
        ),
        "radix2_layers_entered": Metric(
            schedule.radix2_layers, "layers", "exact_schedule",
            "Invocation layer counts; a pruned layer may touch only selected blocks."
        ),
        "nontransform_xor_vectors": Metric(
            schedule.nontransform_xor_vectors, "shard_vectors",
            schedule.arithmetic_classification
        ),
        "fused_multiply_add_vectors": Metric(
            schedule.fused_multiply_add_vectors, "shard_vectors",
            "exact_schedule",
            "Each vector combines a fixed multiplication and accumulator XOR "
            "in one source/destination memory pass."
        ),
        "logical_copy_vectors": Metric(
            schedule.copies, "shard_vectors", "exact_schedule",
            "Full kernel-work copies only. Decode output/tail boundary copies "
            "are reported as exact payload bytes because they may be ragged."
        ),
        "logical_zero_fill_vectors": Metric(
            schedule.zero_fills, "shard_vectors", "exact_schedule"
        ),
        "field_additions_lower": Metric(
            xor_vectors_lower * symbols, "field_additions",
            "derived_lower_bound", arithmetic_note
        ),
        "field_additions_upper": Metric(
            xor_vectors_upper * symbols, "field_additions",
            schedule.arithmetic_classification, arithmetic_note
        ),
        "fixed_multiplications_upper": Metric(
            multiply_vectors_upper * symbols, "fixed_field_multiplications",
            schedule.arithmetic_classification, arithmetic_note
        ),
        "xor_bytes_lower": Metric(
            xor_vectors_lower * kernel_bytes, "bytes", "derived_lower_bound",
            "Counts byte lanes participating in field additions."
        ),
        "xor_bytes_upper": Metric(
            xor_vectors_upper * kernel_bytes, "bytes",
            schedule.arithmetic_classification,
            "Counts byte lanes participating in field additions."
        ),
        "estimated_bytes_read_lower": Metric(
            reads_lower * kernel_bytes, "bytes", "estimated_lower_bound", memory_note
        ),
        "estimated_bytes_read_upper": Metric(
            reads_upper * kernel_bytes, "bytes", "estimated_upper_bound", memory_note
        ),
        "estimated_bytes_written_lower": Metric(
            writes_lower * kernel_bytes, "bytes", "estimated_lower_bound", memory_note
        ),
        "estimated_bytes_written_upper": Metric(
            writes_upper * kernel_bytes, "bytes", "estimated_upper_bound", memory_note
        ),
        "logical_shard_accesses_lower": Metric(
            reads_lower + writes_lower, "shard_accesses", "estimated_lower_bound", memory_note
        ),
        "logical_shard_accesses_upper": Metric(
            reads_upper + writes_upper, "shard_accesses", "estimated_upper_bound", memory_note
        ),
        "estimated_full_workspace_passes_lower": Metric(
            round((reads_lower + writes_lower) / full_pass_denominator, 9),
            "workspace_equivalents", "estimated_lower_bound", memory_note
        ),
        "estimated_full_workspace_passes_upper": Metric(
            round((reads_upper + writes_upper) / full_pass_denominator, 9),
            "workspace_equivalents", "estimated_upper_bound", memory_note
        ),
        "execution_working_slots": Metric(
            schedule.execution_slots, "rounded_shard_slots", "exact_schedule"
        ),
        "kernel_shard_bytes": Metric(
            kernel_bytes, "bytes", "exact_schedule",
            (
                "Direct execution consumes the exact public payload length."
                if schedule.execution_bytes_override is not None else
                "Transform kernels round execution to complete 64-byte tiles."
            )
        ),
    }
    if schedule.decode_coordinate_pointer_mappings is not None:
        metrics["decode_coordinate_pointer_mappings"] = Metric(
            schedule.decode_coordinate_pointer_mappings, "pointer_entries",
            "exact_schedule",
            "Selected caller shards retained by address; this is metadata, not data traffic.",
        )
    if schedule.decode_output_gather_payload_bytes is not None:
        metrics["decode_output_gather_payload_bytes"] = Metric(
            schedule.decode_output_gather_payload_bytes, "bytes",
            schedule.decode_output_gather_classification,
            schedule.decode_output_gather_note,
        )
    if schedule.decode_reveal_temporary_payload_bytes is not None:
        reveal_temporary = schedule.decode_reveal_temporary_payload_bytes
        reveal_scatter = schedule.decode_reveal_scatter_payload_bytes or 0
        reveal_direct = schedule.decode_reveal_direct_payload_bytes or 0
        syndrome_fused = schedule.decode_syndrome_fused_payload_bytes or 0
        syndrome_materialized = (
            schedule.decode_syndrome_materialized_payload_bytes or 0
        )
        auxiliary = reveal_temporary + reveal_scatter
        adjusted_note = (
            "Backend-qualified logical traffic adds the alias-safe in-place "
            "reveal temporary and public scatter, and subtracts one read plus "
            "one write for each syndrome byte accumulated by the final inverse "
            "layer. Schema v3 applies only these reveal/syndrome deltas; "
            "receive-source and weighted-locator effects remain separately "
            "reported structural details. This is logical, not cache traffic."
        )
        metrics.update({
            "decode_reveal_inplace_temporary_payload_bytes": Metric(
                reveal_temporary, "bytes", "exact_schedule",
                "Payload copied to a disjoint 64-byte temporary before a "
                "restrict-qualified in-place reveal multiply; contributes "
                "the same number of logical read and write bytes.",
            ),
            "decode_reveal_scatter_payload_bytes": Metric(
                reveal_scatter, "bytes", "exact_schedule",
                "Public payload copied from transform/requested-output scratch; "
                "contributes the same number of logical read and write bytes.",
            ),
            "decode_reveal_direct_payload_bytes": Metric(
                reveal_direct, "bytes", "exact_schedule",
                "Aligned payload multiplied directly from transform scratch "
                "into caller output by the selected backend policy.",
            ),
            "decode_syndrome_fused_accumulation_payload_bytes": Metric(
                syndrome_fused, "bytes", "exact_schedule",
                "Algorithm 5 output-row bytes whose final inverse layer "
                "accumulates directly into the syndrome workspace.",
            ),
            "decode_syndrome_materialized_xor_payload_bytes": Metric(
                syndrome_materialized, "bytes", "exact_schedule",
                "Algorithm 5 output-row bytes retaining a separate "
                "materialized XOR reduction.",
            ),
            "estimated_backend_bytes_read_lower": Metric(
                reads_lower * kernel_bytes - syndrome_fused + auxiliary,
                "bytes", "estimated_lower_bound", adjusted_note,
            ),
            "estimated_backend_bytes_read_upper": Metric(
                reads_upper * kernel_bytes - syndrome_fused + auxiliary,
                "bytes", "estimated_upper_bound", adjusted_note,
            ),
            "estimated_backend_bytes_written_lower": Metric(
                writes_lower * kernel_bytes - syndrome_fused + auxiliary,
                "bytes", "estimated_lower_bound", adjusted_note,
            ),
            "estimated_backend_bytes_written_upper": Metric(
                writes_upper * kernel_bytes - syndrome_fused + auxiliary,
                "bytes", "estimated_upper_bound", adjusted_note,
            ),
        })
    scratch = schedule.decode_scratch
    if scratch is not None:
        scratch_note = (
            "Exact for the declared pointer width; includes AddressRange and "
            "pointer-map metadata, 64-byte alignment, ragged input staging, "
            "and the selected decode workspace."
        )
        metrics.update({
            "decode_plan_scratch_bytes": Metric(
                scratch.plan_total_bytes, "bytes", "exact_schedule", scratch_note
            ),
            "decode_codec_scratch_bytes": Metric(
                scratch.codec_total_bytes, "bytes", "exact_schedule",
                "One-shot codec query; high tiled layout reserves R output slots "
                "because an erasure pattern is not yet known. " + scratch_note,
            ),
            "decode_plan_range_metadata_bytes": Metric(
                scratch.range_metadata_bytes, "bytes", "exact_schedule"
            ),
            "decode_plan_pointer_map_bytes": Metric(
                scratch.plan_pointer_map_bytes, "bytes", "exact_schedule",
                "Pointer mappings are metadata, not shard-data copies."
            ),
            "decode_plan_metadata_alignment_bytes": Metric(
                scratch.plan_metadata_alignment_bytes, "bytes", "exact_schedule"
            ),
            "decode_plan_tail_reserved_bytes": Metric(
                scratch.tail_reserved_bytes, "bytes", "exact_schedule",
                "Ragged transform plans reserve K+R fixed 64-byte public-coordinate slots."
            ),
            "decode_tail_staged_payload_bytes": Metric(
                scratch.tail_staged_payload_bytes, "bytes", "exact_schedule",
                "Exactly K selected input coordinates are copied for a valid ragged pass."
            ),
            "decode_tail_staged_zero_padding_bytes": Metric(
                scratch.tail_staged_zero_padding_bytes, "bytes", "exact_schedule"
            ),
            "decode_work_slot_bytes": Metric(
                scratch.work_slot_bytes, "bytes", "exact_schedule"
            ),
            "decode_plan_work_slots": Metric(
                scratch.plan_work_slots, "shard_slots", "exact_schedule"
            ),
            "decode_plan_work_bytes": Metric(
                scratch.plan_work_bytes, "bytes", "exact_schedule"
            ),
            "decode_codec_work_slots": Metric(
                scratch.codec_work_slots, "shard_slots", "exact_schedule"
            ),
            "decode_codec_work_bytes": Metric(
                scratch.codec_work_bytes, "bytes", "exact_schedule"
            ),
            "decode_codec_range_metadata_bytes": Metric(
                scratch.codec_range_metadata_bytes, "bytes", "exact_schedule"
            ),
            "decode_codec_pointer_map_bytes": Metric(
                scratch.codec_pointer_map_bytes, "bytes", "exact_schedule"
            ),
            "decode_codec_metadata_alignment_bytes": Metric(
                scratch.codec_metadata_alignment_bytes, "bytes", "exact_schedule"
            ),
            "decode_codec_tail_reserved_bytes": Metric(
                scratch.codec_tail_reserved_bytes, "bytes", "exact_schedule",
                "Pattern-independent ragged query reserves K+R fixed 64-byte slots."
            ),
        })
    return metrics


def attach_decode_scratch(
    schedule: Schedule,
    k: int,
    r: int,
    parent: int,
    padded: int,
    profile: str,
    shard_bytes: int,
    losses: Set[int],
    plan_workspace: str,
    pointer_bytes: int,
    codec_workspace: Optional[str] = None,
    direct: bool = False,
) -> None:
    scratch = decode_scratch_accounting(
        k, r, parent, padded, profile, shard_bytes, len(losses),
        plan_workspace, pointer_bytes, codec_workspace=codec_workspace,
        no_op=not losses, direct=direct and bool(losses),
    )
    schedule.decode_scratch = scratch
    schedule.execution_slots = scratch.plan_work_slots
    schedule.api_scratch_slots = 0
    schedule.details.update({
        "decode_workspace": scratch.workspace,
        "pointer_bytes": pointer_bytes,
        "coordinate_pointer_map_entries": (
            parent if scratch.plan_work_slots != 0 else 0
        ),
        "work_pointer_entries": scratch.plan_work_slots,
        "tail_reserved_input_slots": scratch.tail_reserved_slots,
        "tail_selected_staged_slots": scratch.tail_selected_staged_slots,
        "scratch_accounting": "exact_public_api_bytes",
    })


def build_report(args: argparse.Namespace) -> Dict[str, object]:
    path = args.path
    profile = _profile_for_path(path, args.profile)
    geometry = parent_geometry(args.k, args.r, profile)
    parent = geometry["parent_count"]
    padded = geometry["padded_side"]
    validate_field(parent, args.field)
    if args.field == "gf16" and args.shard_bytes % 2:
        raise ModelError("current GF16 profile requires an even shard byte count")

    requested = parse_mask(args.requested_parity, args.r, default_all=True)
    if args.loss_mask is not None and args.loss_count is not None:
        raise ModelError("use either --loss-count or --loss-mask, not both")
    if args.loss_mask is not None:
        losses = parse_mask(args.loss_mask, args.k, default_all=False)
    else:
        loss_count = 1 if args.loss_count is None else args.loss_count
        if loss_count < 0 or loss_count > args.k:
            raise ModelError("loss count outside [0, K]")
        losses = set(range(loss_count))

    decode_paths = {
        "legacy_high_decode", "low_decode", "generic_decode", "direct_repair"
    }
    selection: Optional[DecodeSelection] = None
    if path in decode_paths:
        if path == "direct_repair":
            selection = DecodeSelection("direct", "direct", 0, 0)
        else:
            if path == "generic_decode":
                if args.decode_selection != "path":
                    raise ModelError(
                        "generic_decode already names its selected path"
                    )
                force = "generic"
            elif args.decode_selection == "path":
                force = args.decode_workspace
            else:
                force = args.decode_selection
            selection = select_decode_execution(
                args.k, args.r, profile, args.field, args.backend, geometry,
                args.shard_bytes, len(losses), force,
                multi_item_batch=args.multi_item_batch,
            )

    if path == "legacy_high_encode":
        schedule = model_high_encode(args.k, args.r, padded, args.shard_bytes, requested)
    elif path == "low_encode":
        schedule = model_low_encode(args.k, args.r, padded, args.shard_bytes, requested)
    elif path in ("legacy_high_decode", "low_decode", "generic_decode"):
        assert selection is not None
        if selection.path == "no_op":
            schedule = Schedule(details={
                "no_op": True, "missing_originals": "none"
            })
        elif selection.path == "direct":
            if padded == 1:
                schedule = model_direct_terminal(
                    args.k, profile, losses, args.shard_bytes
                )
            else:
                schedule = model_direct_repair(
                    args.k, args.r, geometry["parent_dimension"], padded,
                    profile, losses, args.shard_bytes,
                )
        elif selection.path == "generic":
            schedule = model_generic_decode(
                args.k, args.r, parent, padded, profile,
                args.shard_bytes, losses,
            )
        elif profile == "high":
            schedule = model_high_decode(
                args.k, args.r, parent, padded, args.shard_bytes,
                losses, args.field,
            )
        else:
            schedule = model_low_decode(
                args.k, args.r, parent, padded, args.shard_bytes, losses
            )
    elif path == "direct_repair":
        schedule = model_direct_repair(
            args.k, args.r, geometry["parent_dimension"], padded, profile,
            losses, args.shard_bytes,
        )
    elif path == "generic_transform":
        transform_requested = parse_mask(
            args.transform_output_mask, parent, default_all=True
        )
        schedule = model_generic_transform(
            parent, args.direction, args.active_input_count, transform_requested
        )
    else:
        raise ModelError("unsupported path: {}".format(path))

    if path in decode_paths:
        assert selection is not None
        apply_decode_backend_accounting(
            schedule, args.k, args.r, profile, args.field, args.backend,
            geometry, args.shard_bytes, losses, selection,
            args.multi_item_batch,
        )
        if selection.path == "generic":
            plan_workspace = "generic"
        elif selection.path in ("materialized", "tiled"):
            plan_workspace = selection.path
        else:
            # Direct/no-op plan layout ignores this value.  The one-shot codec
            # query remains transform-capable and uses the ordinary specialized
            # capacity policy for the same immutable codec.
            plan_workspace = "specialized"
        if path == "generic_decode":
            codec_workspace = "generic"
        elif args.decode_selection == "path" and path != "direct_repair":
            codec_workspace = args.decode_workspace
        else:
            codec_force = (
                args.decode_selection
                if args.decode_selection in ("auto", "specialized") else "auto"
            )
            codec_workspace = select_codec_transform_execution(
                args.k, args.r, profile, args.field, args.backend,
                geometry, args.shard_bytes, codec_force,
            ).path
        direct = selection.path == "direct"
        attach_decode_scratch(
            schedule, args.k, args.r, parent, padded, profile,
            args.shard_bytes, losses, plan_workspace, args.pointer_bytes,
            codec_workspace=codec_workspace, direct=direct,
        )

    metrics = schedule_metrics(schedule, args.field, args.shard_bytes)
    return {
        "schema_version": SCHEMA_VERSION,
        "model": "leopard2_operation_counts_with_backend_traffic",
        "path": path,
        "inputs": {
            "K": args.k,
            "R": args.r,
            "profile": profile,
            "field": args.field,
            "backend": args.backend,
            "shard_bytes": args.shard_bytes,
            "pointer_bytes": args.pointer_bytes,
            "decode_workspace": (
                schedule.decode_scratch.workspace
                if schedule.decode_scratch is not None else None
            ),
            "decode_selection": args.decode_selection,
            "multi_item_batch": args.multi_item_batch,
            "requested_parity": mask_to_ranges(requested),
            "missing_originals": mask_to_ranges(losses),
        },
        "geometry": geometry,
        "scope": {
            "plan_setup_included": False,
            "structural_schedule_isa_independent": True,
            "backend_traffic_predicates_applied": path in decode_paths,
            "schedule_source": "Leopard DIT radix-4 loop topology",
        },
        "selection": (
            {
                "path": selection.path,
                "rule": selection.rule,
                "matching_auto_rules": selection.matching_auto_rules,
                "required_work_slots": selection.required_work_slots,
            } if selection is not None else None
        ),
        "details": schedule.details,
        "metrics": {name: metric.as_dict() for name, metric in sorted(metrics.items())},
    }


def write_json(report: Dict[str, object], output: TextIO) -> None:
    json.dump(report, output, sort_keys=True, indent=2)
    output.write("\n")


def write_csv(report: Dict[str, object], output: TextIO) -> None:
    inputs = report["inputs"]
    geometry = report["geometry"]
    assert isinstance(inputs, dict)
    assert isinstance(geometry, dict)
    metrics = report["metrics"]
    assert isinstance(metrics, dict)
    columns = [
        "schema_version", "path", "K", "R", "profile", "field",
        "backend", "shard_bytes", "pointer_bytes", "decode_workspace",
        "decode_selection", "multi_item_batch", "selected_decode_path",
        "selected_decode_rule", "selected_decode_matching_auto_rules",
        "selected_decode_required_work_slots", "parent_count",
        "padded_side", "metric", "value", "unit", "classification", "note",
    ]
    writer = csv.DictWriter(output, fieldnames=columns, lineterminator="\n")
    writer.writeheader()
    for name in sorted(metrics):
        metric = metrics[name]
        assert isinstance(metric, dict)
        writer.writerow({
            "schema_version": report["schema_version"],
            "path": report["path"],
            "K": inputs["K"],
            "R": inputs["R"],
            "profile": inputs["profile"],
            "field": inputs["field"],
            "backend": inputs["backend"],
            "shard_bytes": inputs["shard_bytes"],
            "pointer_bytes": inputs["pointer_bytes"],
            "decode_workspace": inputs["decode_workspace"],
            "decode_selection": inputs["decode_selection"],
            "multi_item_batch": inputs["multi_item_batch"],
            "selected_decode_path": (
                report["selection"]["path"] if report["selection"] else ""
            ),
            "selected_decode_rule": (
                report["selection"]["rule"] if report["selection"] else ""
            ),
            "selected_decode_matching_auto_rules": (
                report["selection"]["matching_auto_rules"]
                if report["selection"] else ""
            ),
            "selected_decode_required_work_slots": (
                report["selection"]["required_work_slots"]
                if report["selection"] else ""
            ),
            "parent_count": geometry["parent_count"],
            "padded_side": geometry["padded_side"],
            "metric": name,
            "value": metric["value"],
            "unit": metric["unit"],
            "classification": metric["classification"],
            "note": metric.get("note", ""),
        })


def _brute_full_butterflies(size: int) -> int:
    count = 0
    width = 2
    while width <= size:
        count += size // 2
        width *= 2
    return count


def run_self_test(verbose: bool = True) -> None:
    checks = 0
    for exponent in range(0, 17):
        size = 1 << exponent
        expected = _brute_full_butterflies(size)
        assert ifft_prefix_butterflies(size, size) == expected
        assert fft_prefix_butterflies(size, size) == expected
        assert fft_mask_butterflies(size, set(range(size))) == expected
        checks += 3
    assert parent_geometry(240, 16, "high") == {
        "padded_side": 16,
        "parent_count": 256,
        "parent_dimension": 240,
        "shortened_count": 0,
        "punctured_count": 0,
    }
    assert parent_geometry(8, 248, "low")["parent_count"] == 256
    checks += 2

    high = model_high_encode(240, 16, 16, 1024, set(range(16)))
    assert high.butterflies == 512
    assert high.nontransform_xor_vectors == 224
    assert high.execution_slots == 32
    checks += 3

    high_decode = model_high_decode(240, 16, 256, 16, 1024, {0})
    assert high_decode.decode_coordinate_pointer_mappings == 240
    assert high_decode.copies == 240
    assert high_decode.decode_output_gather_payload_bytes == 1024
    assert high_decode.details["locator_weighted_live_scan_scope"] == (
        "statically qualified GF8 AVX2/AVX512 passes scan receive pointers "
        "until ceil(T/2) live rows are found; sparse fallback scans all T"
    )
    low_decode_ragged = model_low_decode(8, 248, 256, 8, 65, {0, 1})
    assert low_decode_ragged.decode_coordinate_pointer_mappings == 8
    assert low_decode_ragged.copies == 0
    assert low_decode_ragged.decode_output_gather_payload_bytes == 2
    scratch = decode_scratch_accounting(
        240, 16, 256, 16, "high", 65, 1, "tiled", 8,
        codec_workspace="tiled",
    )
    assert scratch.plan_total_bytes == 28800
    assert scratch.codec_total_bytes == 29824
    checks += 9

    for pointer_bytes in (4, 8):
        maximum = size_t_max(pointer_bytes)
        largest_roundable = maximum - (SCRATCH_ALIGNMENT - 1)
        assert round_shard_bytes(largest_roundable, pointer_bytes) == \
            largest_roundable
        try:
            round_shard_bytes(maximum, pointer_bytes)
        except ModelError:
            pass
        else:
            raise AssertionError("SIZE_MAX shard escaped round-up rejection")
        try:
            round_shard_bytes(UINT64_MAX, pointer_bytes)
        except ModelError:
            pass
        else:
            raise AssertionError("UINT64_MAX shard escaped round-up rejection")

        control = decode_scratch_accounting(
            1, 1, 2, 1, "low", 64, 1, "specialized", pointer_bytes,
            codec_workspace="specialized", direct=True,
        )
        largest_layout_shard = (
            (maximum - control.codec_data_offset) // control.codec_work_slots
        ) & ~(SCRATCH_ALIGNMENT - 1)
        boundary = decode_scratch_accounting(
            1, 1, 2, 1, "low", largest_layout_shard, 1, "specialized",
            pointer_bytes, codec_workspace="specialized", direct=True,
        )
        assert boundary.codec_total_bytes <= maximum
        try:
            decode_scratch_accounting(
                1, 1, 2, 1, "low",
                largest_layout_shard + SCRATCH_ALIGNMENT, 1, "specialized",
                pointer_bytes, codec_workspace="specialized", direct=True,
            )
        except ModelError:
            pass
        else:
            raise AssertionError("decode scratch overflow neighbor was accepted")
        checks += 5
    try:
        round_shard_bytes(UINT64_MAX + 1, 8)
    except ModelError:
        pass
    else:
        raise AssertionError("shard length above uint64_t was accepted")
    checks += 1

    low = model_low_encode(8, 248, 8, 1024, set(range(248)))
    assert low.butterflies == 384
    assert low.details["active_parity_blocks"] == 31
    assert low.details["direct_output_blocks"] == 31
    assert low.copies == 8
    checks += 4

    low_decode = model_low_decode(8, 248, 256, 8, 1024, {0, 1})
    assert low_decode.fused_multiply_add_vectors == 248
    assert low_decode.details["fused_weighted_multiply_add_vectors"] == 248
    low_decode_metrics = schedule_metrics(low_decode, "gf8", 1024)
    unfused_reads = (
        2 * low_decode.butterflies
        + 2 * low_decode.nontransform_xor_vectors
        + low_decode.fixed_multiply_vectors
        + low_decode.copies
    ) * 1024
    assert (low_decode_metrics["estimated_bytes_read_lower"].value ==
            unfused_reads - 248 * 1024)
    checks += 3

    root = pathlib.Path(__file__).resolve().parent.parent
    cmake_source = (root / "CMakeLists.txt").read_text(encoding="utf-8")
    verify_high_pruned_hook_registration(cmake_source, "CMakeLists.txt")
    checks += 1
    core_source = (root / "leopard2.cpp").read_text(encoding="utf-8")
    verify_decode_scratch_source(core_source, "leopard2.cpp")
    try:
        verify_decode_scratch_source(
            core_source.replace(
                "ComputeScratchLayout(range_count, 0, 0, rounded_bytes, layout)",
                "ComputeScratchLayout(range_count, 1, 0, rounded_bytes, layout)",
                1,
            ),
            "leopard2.cpp direct-pointer mutation",
        )
    except ModelError:
        pass
    else:
        raise AssertionError("direct decode scratch mutation escaped source guard")
    checks += 2
    policy_source = (root / "Leopard2Dispatch.h").read_text(encoding="utf-8")
    verify_decode_policy_source(policy_source, "Leopard2Dispatch.h")
    policy_mutations = (
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
    for old, new in policy_mutations:
        mutated = policy_source.replace(old, new, 1)
        if mutated == policy_source:
            raise AssertionError("decode policy mutation did not apply")
        try:
            verify_decode_policy_source(
                mutated, "Leopard2Dispatch.h backend mutation"
            )
        except ModelError:
            pass
        else:
            raise AssertionError("decode policy mutation escaped source guard")
    checks += 1 + len(policy_mutations)
    ff8_source = (root / "LeopardFF8.cpp").read_text(encoding="utf-8")
    ff16_source = (root / "LeopardFF16.cpp").read_text(encoding="utf-8")
    verify_decode_fusion_sources(core_source, ff8_source, ff16_source)
    fusion_mutations = (
        (
            core_source.replace(
                "aligned_prefix_bytes >= 4096",
                "aligned_prefix_bytes >= 2048", 1,
            ),
            ff8_source,
            ff16_source,
        ),
        (
            core_source,
            ff8_source.replace(
                "buffer_bytes >= 65536", "buffer_bytes >= 32768", 1
            ),
            ff16_source,
        ),
        (
            core_source,
            ff8_source,
            ff16_source.replace("bytes == 64", "bytes == 32", 1),
        ),
        (
            core_source.replace(
                "const bool reveal_aligned_outputs_in_place =\n"
                "        !(fuse_generic_reveal_scatter || "
                "fuse_low_reveal_scatter);",
                "const bool reveal_aligned_outputs_in_place = true;",
                1,
            ),
            ff8_source,
            ff16_source,
        ),
        (
            core_source.replace(
                "ExecuteTransformDecodePass(\n"
                "            plan, kScratchAlignment, coordinate_input, "
                "work,\n"
                "            use_generic, use_tiled, true);",
                "ExecuteTransformDecodePass(\n"
                "            plan, kScratchAlignment, coordinate_input, "
                "work,\n"
                "            use_generic, use_tiled, "
                "reveal_aligned_outputs_in_place);",
                1,
            ),
            ff8_source,
            ff16_source,
        ),
        (
            core_source,
            ff8_source.replace(
                "const bool use_accumulating_sink =\n"
                "                ShouldUsePrunedHighSyndromeSink("
                "ops, buffer_bytes);",
                "const bool use_accumulating_sink = false;",
                1,
            ),
            ff16_source,
        ),
        (
            core_source,
            ff8_source,
            ff16_source.replace(
                "const bool use_accumulating_sink =\n"
                "                ShouldUsePrunedHighSyndromeSink("
                "ops, buffer_bytes);",
                "const bool use_accumulating_sink = false;",
                1,
            ),
        ),
    )
    for mutated_core, mutated_ff8, mutated_ff16 in fusion_mutations:
        if (mutated_core == core_source and mutated_ff8 == ff8_source and
                mutated_ff16 == ff16_source):
            raise AssertionError("decode fusion mutation did not apply")
        try:
            verify_decode_fusion_sources(
                mutated_core, mutated_ff8, mutated_ff16
            )
        except ModelError:
            pass
        else:
            raise AssertionError("decode fusion mutation escaped source guard")
    checks += 1 + len(fusion_mutations)

    low_call_fixture = """
        FFT_DIT_FromCoefficients(
            ops, buffer_bytes, work, work + p, requested_count, p,
            FFTSkewStorage + p, SourceEvaluationLowEncode);
    """
    nested_hook_fixture = low_call_fixture + """
        #if defined(LEO2_ENABLE_TEST_HOOKS)
        if (callsite == SourceEvaluationHighDecode &&
            TestForceHighDecodeCopyFallback.load())
        {
        #if defined(LEO2_NESTED_TEST_FIXTURE)
            (void)callsite;
        #endif
            memcpy(evaluation_work[i], coefficients[i], bytes);
        }
        #endif
    """
    verify_low_encode_no_copy_source(
        nested_hook_fixture, "nested test-hook fixture")
    escaped_hook_fixtures = (
        low_call_fixture + """
            #if defined(LEO2_ENABLE_TEST_HOOKS)
            if (callsite == SourceEvaluationHighDecode &&
                TestForceHighDecodeCopyFallback.load()) {}
            #endif
            memcpy(evaluation_work[i], coefficients[i], bytes);
        """,
        low_call_fixture + """
            #if defined(LEO2_ENABLE_TEST_HOOKS)
            (void)callsite;
            #else
            if (callsite == SourceEvaluationHighDecode &&
                TestForceHighDecodeCopyFallback.load())
                memcpy(evaluation_work[i], coefficients[i], bytes);
            #endif
        """,
    )
    for fixture in escaped_hook_fixtures:
        try:
            verify_low_encode_no_copy_source(fixture, "escaped hook fixture")
        except ModelError:
            pass
        else:
            raise AssertionError(
                "copy outside the canonical positive hook branch escaped")
    checks += 3

    for filename in ("LeopardFF8.cpp", "LeopardFF16.cpp"):
        source = ff8_source if filename == "LeopardFF8.cpp" else ff16_source
        verify_low_encode_no_copy_source(source, filename)
        verify_high_decode_no_copy_source(source, filename)
        if filename == "LeopardFF8.cpp":
            verify_high_decode_receive_fusion_source(source, filename)
            verify_high_decode_weighted_locator_source(source, filename)
        else:
            verify_high_decode_gf16_copy_first_source(source, filename)
        verify_low_decode_weighted_fusion_source(source, filename)
        try:
            verify_low_encode_no_copy_source(
                source +
                "\nfor (unsigned i = 0; i < p; ++i) "
                "memcpy(work[p + i], work[i], buffer_bytes);\n",
                filename + " mutation",
            )
        except ModelError:
            pass
        else:
            raise AssertionError("whole-P copy mutation escaped source guard")
        try:
            verify_high_decode_no_copy_source(
                source + "\nfor (unsigned i = 0; i < t; ++i) "
                "memcpy(tile[i], accumulator[i], buffer_bytes);\n",
                filename + " high mutation",
            )
        except ModelError:
            pass
        else:
            raise AssertionError("whole-T copy mutation escaped source guard")
        checks += 6 if filename == "LeopardFF8.cpp" else 5

    direct = model_direct_repair(16, 8, 16, 16, "low", set(range(4)))
    assert direct.nontransform_xor_vectors == 60
    assert direct.fixed_multiply_vectors == 64
    assert direct.api_scratch_slots == 0
    checks += 3

    for size in (2, 4, 8, 16, 32, 64):
        previous = -1
        for prefix in range(size + 1):
            value = fft_prefix_butterflies(size, prefix)
            assert value >= previous
            previous = value
            assert value == fft_mask_butterflies(size, set(range(prefix)))
            checks += 2
    assert mask_to_ranges(parse_mask("0-3,7,9-10", 16, False)) == "0-3,7,9-10"
    checks += 1

    namespace = argparse.Namespace(
        path="legacy_high_encode", k=9, r=7, profile=None, field="gf8",
        backend="scalar",
        shard_bytes=65, requested_parity="0-3,6", loss_count=None,
        loss_mask=None, direction="forward", active_input_count=None,
        transform_output_mask=None, pointer_bytes=struct.calcsize("P"),
        decode_workspace="materialized", decode_selection="path",
        multi_item_batch=False,
    )
    first = json.dumps(build_report(namespace), sort_keys=True)
    second = json.dumps(build_report(namespace), sort_keys=True)
    assert first == second
    checks += 1
    if verbose:
        print("operation-count self-test passed: {} checks".format(checks))


def _add_model_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--path", choices=PATHS, required=True)
    parser.add_argument("--k", type=int, required=True, help="public original count")
    parser.add_argument("--r", type=int, required=True, help="public recovery count")
    parser.add_argument("--profile", choices=PROFILES,
                        help="required for generic/direct paths; inferred otherwise")
    parser.add_argument("--field", choices=FIELDS, required=True)
    parser.add_argument(
        "--backend", choices=BACKENDS, default="scalar",
        help="production backend used for decode fusion traffic (default: scalar)",
    )
    parser.add_argument("--shard-bytes", type=int, required=True)
    parser.add_argument(
        "--pointer-bytes", type=int, choices=(4, 8),
        default=struct.calcsize("P"),
        help="ABI pointer width used for exact decode metadata accounting",
    )
    parser.add_argument(
        "--decode-workspace", choices=DECODE_WORKSPACES,
        default="materialized",
        help="forced specialized decode workspace (default: materialized)",
    )
    parser.add_argument(
        "--decode-selection", choices=DECODE_SELECTIONS, default="path",
        help=(
            "path forces the named generic/specialized path; auto mirrors "
            "production AUTO; specialized mirrors FORCE_SPECIALIZED"
        ),
    )
    parser.add_argument(
        "--multi-item-batch", action="store_true",
        help="apply the production multi-item batch decode selector rule",
    )
    parser.add_argument(
        "--requested-parity", default=None,
        help="all, none, or comma-separated indices/ranges (default: all)",
    )
    parser.add_argument("--loss-count", type=int,
                        help="missing originals 0..count-1 (decode default: 1)")
    parser.add_argument("--loss-mask",
                        help="missing-original indices/ranges; overrides no defaults")
    parser.add_argument("--direction", choices=("forward", "inverse"),
                        default="forward", help="generic_transform direction")
    parser.add_argument("--active-input-count", type=int,
                        help="inverse generic_transform input prefix")
    parser.add_argument("--transform-output-mask",
                        help="forward generic_transform parent-coordinate mask")


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    report = subparsers.add_parser("report", help="emit one deterministic report")
    _add_model_arguments(report)
    report.add_argument("--format", choices=("json", "csv"), default="json")
    report.add_argument("--output", type=pathlib.Path,
                        help="write to a file instead of stdout")
    subparsers.add_parser("self-test", help="run deterministic model invariants")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = make_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "self-test":
            run_self_test()
            return 0
        report = build_report(args)
        if args.output:
            with args.output.open("w", encoding="utf-8", newline="") as output:
                if args.format == "json":
                    write_json(report, output)
                else:
                    write_csv(report, output)
        elif args.format == "json":
            write_json(report, sys.stdout)
        else:
            write_csv(report, sys.stdout)
        return 0
    except (ModelError, OSError, ValueError) as error:
        parser.error(str(error))
    return 2


if __name__ == "__main__":
    sys.exit(main())
