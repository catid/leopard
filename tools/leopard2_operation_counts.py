#!/usr/bin/env python3
"""Deterministic, ISA-independent operation model for Leopard2 paths.

This utility models the transform schedules implemented by Leopard's radix-4
wrappers in terms of their equivalent radix-2 butterflies.  It deliberately
does not pretend that a scheduled butterfly always executes a multiplication:
zero skew factors are specialized by the production kernels.  Arithmetic
totals involving those factors are therefore emitted as bounds.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import pathlib
import re
import sys
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, TextIO, Tuple


SCHEMA_VERSION = 1
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


_LOW_COEFFICIENT_COPY_PATTERNS = (
    re.compile(
        r"(?:std::)?memcpy\(work\[p\+[^\]]+\],work\[[^\]]+\],buffer_bytes\)"
    ),
    re.compile(
        r"(?:std::)?memcpy\(evaluation_work\[[^\]]+\],"
        r"coefficients\[[^\]]+\],bytes\)"
    ),
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
    compact = re.sub(r"\s+", "", source)
    for pattern in _LOW_COEFFICIENT_COPY_PATTERNS:
        if pattern.search(compact):
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
    if "ops.kind==LEO2_BACKEND_SSSE3" not in helper or \
            "ops.kind==LEO2_BACKEND_AVX2" not in helper or \
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

    def add_transform(self, count: int, layers: int) -> None:
        self.butterflies += count
        self.transform_invocations += 1
        self.radix2_layers += layers


def is_power_of_two(value: int) -> bool:
    return value > 0 and (value & (value - 1)) == 0


def ceil_power_of_two(value: int) -> int:
    if value <= 0:
        raise ModelError("counts must be positive")
    return 1 << (value - 1).bit_length()


def ceil_div(numerator: int, denominator: int) -> int:
    return (numerator + denominator - 1) // denominator


def rounded_kernel_bytes(shard_bytes: int) -> int:
    if shard_bytes <= 0:
        raise ModelError("shard bytes must be positive")
    return ceil_div(shard_bytes, 64) * 64


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
        schedule.copies += len(block_mask)
    if partial:
        schedule.copies += k + len(requested)
    schedule.details.update({
        "active_parity_blocks": active_blocks,
        "parity_block_prefixes": block_prefixes,
        "out_of_place_first_fft_layer": True,
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
    schedule.api_scratch_slots = k + r + parent
    schedule.execution_slots = parent
    # leo2_decode_plan_execute stages exactly the selected K coordinates and
    # gathers only the missing originals after transform execution.
    schedule.copies = k + len(losses)
    schedule.details.update({
        "missing_originals": mask_to_ranges(losses),
        "selected_parity": mask_to_ranges(selected_parities),
        "received_coordinate_count": len(data),
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
            copies=1,
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
            copies=1,
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
    prefixes: List[int] = []
    for offset in range(0, parent, padded):
        prefix = _block_prefix(data, offset, padded)
        prefixes.append(prefix)
        live = {coordinate - offset for coordinate in data
                if offset <= coordinate < offset + padded}
        if not live and offset != 0:
            receive_skipped_blocks += 1
            continue
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
    # SSSE3/AVX2 delta, not a backend-neutral replacement for those totals.
    # A future backend-aware model can apply the delta to absolute traffic.
    schedule.copies += k
    schedule.zero_fills += parent - k
    reduction = (parent // padded - 1) * padded
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
    schedule.details.update({
        "block_input_prefixes": prefixes,
        "receive_source_fused_radix4_groups": receive_fused_groups,
        "receive_copy_vectors": receive_copy_vectors,
        "receive_zero_vectors": receive_zero_vectors,
        "receive_copy_vectors_removed": 4 * receive_fused_groups,
        "receive_exact_pruned_staged_blocks": receive_exact_pruned_blocks,
        "receive_skipped_empty_blocks": receive_skipped_blocks,
        "receive_source_fusion_scope": (
            "GF8 qualified SSSE3/AVX2 mature unpruned schedule delta"
            if field_name == "gf8" else
            "GF16 deterministic copy-first policy"
        ),
        "syndrome_reduction_vectors": reduction,
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
        details={
            "missing_originals": mask_to_ranges(losses),
            "dense_term_upper_bound": dense_terms,
            "coefficient_dependent_terms": True,
            "profile": profile,
        },
    )
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
    kernel_bytes = rounded_kernel_bytes(shard_bytes)
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
        "Logical shard traffic only; excludes metadata, allocator effects, "
        "cache-line amplification, and backend table loads."
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
            schedule.copies, "shard_vectors", "exact_schedule"
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
        "api_scratch_data_slots": Metric(
            schedule.api_scratch_slots, "rounded_shard_slots", "exact_schedule",
            "Excludes pointer/range metadata and alignment padding."
        ),
        "execution_working_slots": Metric(
            schedule.execution_slots, "rounded_shard_slots", "exact_schedule"
        ),
        "kernel_shard_bytes": Metric(
            kernel_bytes, "bytes", "exact_schedule",
            "Current kernels round execution to complete 64-byte tiles."
        ),
    }
    return metrics


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

    if path == "legacy_high_encode":
        schedule = model_high_encode(args.k, args.r, padded, args.shard_bytes, requested)
    elif path == "low_encode":
        schedule = model_low_encode(args.k, args.r, padded, args.shard_bytes, requested)
    elif path == "legacy_high_decode":
        schedule = model_high_decode(
            args.k, args.r, parent, padded, args.shard_bytes, losses, args.field
        )
    elif path == "low_decode":
        schedule = model_low_decode(
            args.k, args.r, parent, padded, args.shard_bytes, losses
        )
    elif path == "generic_decode":
        schedule = model_generic_decode(
            args.k, args.r, parent, padded, profile, args.shard_bytes, losses
        )
    elif path == "direct_repair":
        schedule = model_direct_repair(
            args.k, args.r, geometry["parent_dimension"], padded, profile, losses
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

    metrics = schedule_metrics(schedule, args.field, args.shard_bytes)
    return {
        "schema_version": SCHEMA_VERSION,
        "model": "leopard2_isa_independent_operation_counts",
        "path": path,
        "inputs": {
            "K": args.k,
            "R": args.r,
            "profile": profile,
            "field": args.field,
            "shard_bytes": args.shard_bytes,
            "requested_parity": mask_to_ranges(requested),
            "missing_originals": mask_to_ranges(losses),
        },
        "geometry": geometry,
        "scope": {
            "plan_setup_included": False,
            "isa_independent": True,
            "schedule_source": "Leopard DIT radix-4 loop topology",
        },
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
        "shard_bytes", "parent_count", "padded_side", "metric", "value",
        "unit", "classification", "note",
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
            "shard_bytes": inputs["shard_bytes"],
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

    low = model_low_encode(8, 248, 8, 1024, set(range(248)))
    assert low.butterflies == 384
    assert low.details["active_parity_blocks"] == 31
    checks += 2

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
    for filename in ("LeopardFF8.cpp", "LeopardFF16.cpp"):
        source = (root / filename).read_text(encoding="utf-8")
        verify_low_encode_no_copy_source(source, filename)
        verify_high_decode_no_copy_source(source, filename)
        if filename == "LeopardFF8.cpp":
            verify_high_decode_receive_fusion_source(source, filename)
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
        checks += 5

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
        shard_bytes=65, requested_parity="0-3,6", loss_count=None,
        loss_mask=None, direction="forward", active_input_count=None,
        transform_output_mask=None,
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
    parser.add_argument("--shard-bytes", type=int, required=True)
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
