/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the name of Leopard-RS nor the names of its contributors may be
      used to endorse or promote products derived from this software without
      specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once

#include "leopard2.h"

#include <stddef.h>
#include <stdint.h>

namespace leopard2_internal {

enum DecodePath
{
    kDecodePathNoOp = 0,
    kDecodePathDirect,
    kDecodePathGeneric,
    kDecodePathMaterialized,
    kDecodePathTiled
};

enum DecodePathRule
{
    kDecodeRuleNoOp = 0,
    kDecodeRuleDirect,
    kDecodeRuleForcedGeneric,
    kDecodeRuleForcedMaterialized,
    kDecodeRuleForcedTiled,
    kDecodeRuleBalancedGeneric,
    kDecodeRuleMeasuredBatchTiled,
    kDecodeRuleMeasuredMaterialized,
    kDecodeRuleWorkspaceTiled,
    kDecodeRuleWorkspaceMaterialized,
    kDecodeRuleTranslatedLowDiagnostic,
    kDecodeRuleUnsupportedProfile
};

enum DecodeAutoRule
{
    kDecodeAutoRuleBalancedGeneric = 1u << 0,
    kDecodeAutoRuleMeasuredMaterialized = 1u << 1
};

struct DecodePathInput
{
    leo2_profile profile;
    leo2_field field;
    leo2_backend backend;
    uint32_t original_count;
    uint32_t recovery_count;
    uint32_t padded_side;
    uint32_t parent_count;
    uint32_t missing_original_count;
    uint32_t requested_output_count;
    uint32_t codec_flags;
    uint64_t actual_shard_bytes;
    size_t aligned_prefix_bytes;
    size_t tail_bytes;
    size_t rounded_shard_bytes;
    bool plan_known;
    /*
        Execution intentionally exposes only the stable single-item versus
        multi-item class to this selector.  The batch entry point already
        knows its exact item count, but callers size each item's decode
        scratch through the ordinary plan query.  A batch-only rule must
        therefore retain or reduce that conservative one-item allocation.

        Expected plan reuse is not a selector input.  Setup is complete before
        this byte-heavy decision, immutable plans may be shared by callers
        with different reuse behavior, and benchmarks report setup separately
        at several amortization counts.  A future caller hint, if evidence
        justifies one, requires an additive versioned API rather than changing
        this V1 default or consuming a reserved public field.
    */
    bool multi_item_batch;
};

struct DecodePathInfo
{
    DecodePath path;
    DecodePathRule rule;
    uint32_t matching_auto_rules;
    size_t required_work_slots;
    size_t aligned_prefix_bytes;
    size_t tail_bytes;
    size_t rounded_shard_bytes;
    bool multi_item_batch;
};

/*
    Pure policy predicate kept separate from execution so every boundary of the
    evidence-scoped automatic crossover can be tested without timing inference.
*/
static inline bool ShouldUseBalancedGenericDecode(
    leo2_profile profile,
    leo2_field field,
    uint32_t original_count,
    uint32_t recovery_count,
    uint32_t padded_side,
    uint32_t parent_count,
    uint32_t missing_original_count,
    size_t rounded_shard_bytes,
    leo2_backend backend)
{
    return profile == LEO2_PROFILE_LEGACY_HIGH_V1 &&
        field == LEO2_FIELD_GF8 &&
        original_count == 128 && recovery_count == 128 &&
        padded_side == 128 && parent_count == 256 &&
        missing_original_count == 128 &&
        rounded_shard_bytes >= 256 && rounded_shard_bytes <= 1024 * 1024 &&
        (backend == LEO2_BACKEND_SCALAR ||
         backend == LEO2_BACKEND_SSSE3 ||
         backend == LEO2_BACKEND_AVX2 ||
         backend == LEO2_BACKEND_AVX512);
}

/*
    Offline-calibrated exception to the side-sized high decoder.  Both kernels
    implement the same legacy-high mathematics and wire profile; this policy
    selects the retained regular N-slot traversal only where three-round pinned
    evidence found a credible tiled regression.  Batch calls are handled by the
    caller because two or more stripes restore the tiled cache advantage.
*/
static inline bool ShouldUseMaterializedHighDecode(
    leo2_profile profile,
    leo2_field field,
    uint32_t original_count,
    uint32_t recovery_count,
    uint32_t padded_side,
    uint32_t parent_count,
    uint32_t missing_original_count,
    size_t rounded_shard_bytes,
    leo2_backend backend)
{
    if (profile != LEO2_PROFILE_LEGACY_HIGH_V1 ||
        field != LEO2_FIELD_GF8 ||
        original_count != 224 || recovery_count != 32 ||
        padded_side != 32 || parent_count != 256 ||
        missing_original_count == 0 || missing_original_count > 8 ||
        rounded_shard_bytes > 64 * 1024)
        return false;

    if (backend == LEO2_BACKEND_AVX2 || backend == LEO2_BACKEND_AVX512)
        return rounded_shard_bytes >= 24 * 1024;
    if (backend == LEO2_BACKEND_SSSE3)
        return rounded_shard_bytes >= 32 * 1024;
    return false;
}

static inline bool IsDecodeByteGeometryConsistent(
    const DecodePathInput& input)
{
    if (input.actual_shard_bytes == 0 ||
        input.actual_shard_bytes > static_cast<uint64_t>(SIZE_MAX) ||
        input.tail_bytes >= 64 ||
        (input.aligned_prefix_bytes & 63u) != 0 ||
        input.aligned_prefix_bytes > SIZE_MAX - input.tail_bytes ||
        input.aligned_prefix_bytes + input.tail_bytes !=
            static_cast<size_t>(input.actual_shard_bytes))
        return false;
    if (input.tail_bytes != 0 && input.aligned_prefix_bytes > SIZE_MAX - 64)
        return false;
    const size_t expected_rounded = input.tail_bytes == 0
        ? input.aligned_prefix_bytes : input.aligned_prefix_bytes + 64;
    return input.rounded_shard_bytes == expected_rounded;
}

/*
    One ordered transform-decode selector owns both workspace sizing and
    transform execution.  No-op and direct repair are immutable terminal plan
    states that bypass transform workspace and are reported by the private
    introspection wrapper before this selector is called.  The actual byte
    split is deliberately part of the input even though the currently promoted
    crossover rules use the 64-byte-rounded size: future rules cannot
    accidentally receive only a rounded surrogate and ignore a ragged pass.
    matching_auto_rules excludes the capacity fallback and lets tests prove
    that evidence-scoped rules do not overlap.
*/
static inline bool SelectDecodePath(
    const DecodePathInput& input,
    DecodePathInfo& selection)
{
    selection.path = kDecodePathMaterialized;
    selection.rule = kDecodeRuleUnsupportedProfile;
    selection.matching_auto_rules = 0;
    selection.required_work_slots = 0;
    selection.aligned_prefix_bytes = input.aligned_prefix_bytes;
    selection.tail_bytes = input.tail_bytes;
    selection.rounded_shard_bytes = input.rounded_shard_bytes;
    selection.multi_item_batch = input.multi_item_batch;

    const uint32_t supported_flags =
        LEO2_CODEC_FORCE_GENERIC_DECODE |
        LEO2_CODEC_FORCE_SPECIALIZED_DECODE |
        LEO2_CODEC_FORCE_TILED_DECODE |
        LEO2_CODEC_FORCE_MATERIALIZED_DECODE;
    const uint32_t algorithm_flags = input.codec_flags &
        (LEO2_CODEC_FORCE_GENERIC_DECODE |
         LEO2_CODEC_FORCE_SPECIALIZED_DECODE);
    const uint32_t workspace_flags = input.codec_flags &
        (LEO2_CODEC_FORCE_TILED_DECODE |
         LEO2_CODEC_FORCE_MATERIALIZED_DECODE);
    if (!IsDecodeByteGeometryConsistent(input) ||
        input.original_count == 0 || input.recovery_count == 0 ||
        input.padded_side == 0 || input.parent_count == 0 ||
        input.missing_original_count > input.original_count ||
        input.missing_original_count > input.recovery_count ||
        input.requested_output_count > input.recovery_count ||
        (input.codec_flags & ~supported_flags) != 0 ||
        algorithm_flags == (LEO2_CODEC_FORCE_GENERIC_DECODE |
                            LEO2_CODEC_FORCE_SPECIALIZED_DECODE) ||
        workspace_flags == (LEO2_CODEC_FORCE_TILED_DECODE |
                            LEO2_CODEC_FORCE_MATERIALIZED_DECODE) ||
        ((algorithm_flags & LEO2_CODEC_FORCE_GENERIC_DECODE) != 0 &&
         workspace_flags != 0))
        return false;

    if (static_cast<size_t>(input.padded_side) > SIZE_MAX / 2)
        return false;
    size_t tiled_work_slots = static_cast<size_t>(input.padded_side) * 2;
    if (input.profile == LEO2_PROFILE_LEGACY_HIGH_V1)
    {
        if (tiled_work_slots > SIZE_MAX - input.requested_output_count)
            return false;
        tiled_work_slots += input.requested_output_count;
    }
    const size_t materialized_work_slots = input.parent_count;

    const bool balanced_generic = input.plan_known &&
        ShouldUseBalancedGenericDecode(
            input.profile, input.field, input.original_count,
            input.recovery_count, input.padded_side, input.parent_count,
            input.missing_original_count, input.rounded_shard_bytes,
            input.backend);
    const uint32_t materialized_missing = input.plan_known
        ? input.missing_original_count : 1;
    const bool measured_materialized = ShouldUseMaterializedHighDecode(
        input.profile, input.field, input.original_count,
        input.recovery_count, input.padded_side, input.parent_count,
        materialized_missing, input.rounded_shard_bytes, input.backend);
    if (balanced_generic)
        selection.matching_auto_rules |= kDecodeAutoRuleBalancedGeneric;
    if (measured_materialized)
        selection.matching_auto_rules |= kDecodeAutoRuleMeasuredMaterialized;

    if ((input.codec_flags & LEO2_CODEC_FORCE_GENERIC_DECODE) != 0)
    {
        selection.path = kDecodePathGeneric;
        selection.rule = kDecodeRuleForcedGeneric;
        selection.required_work_slots = materialized_work_slots;
        return true;
    }
    if ((input.codec_flags & LEO2_CODEC_FORCE_MATERIALIZED_DECODE) != 0)
    {
        selection.path = kDecodePathMaterialized;
        selection.rule = kDecodeRuleForcedMaterialized;
        selection.required_work_slots = materialized_work_slots;
        return true;
    }
    if ((input.codec_flags & LEO2_CODEC_FORCE_TILED_DECODE) != 0)
    {
        selection.path = kDecodePathTiled;
        selection.rule = kDecodeRuleForcedTiled;
        /* Codec construction currently rejects unsupported profiles, but the
           retained fallback historically reserved N slots before the
           forced execution predicate selected tiled traversal.  Preserve that
           conservative contract for private selector probes and future
           versioned profiles until their workspace is defined explicitly. */
        selection.required_work_slots =
            input.profile == LEO2_PROFILE_LOW_V1 ||
            input.profile == LEO2_PROFILE_LEGACY_HIGH_V1
                ? tiled_work_slots : materialized_work_slots;
        return true;
    }

    const bool force_specialized =
        (input.codec_flags & LEO2_CODEC_FORCE_SPECIALIZED_DECODE) != 0;
    if (!force_specialized && balanced_generic)
    {
        selection.path = kDecodePathGeneric;
        selection.rule = kDecodeRuleBalancedGeneric;
        selection.required_work_slots = materialized_work_slots;
        return true;
    }
    if (measured_materialized)
    {
        if (input.multi_item_batch &&
            (input.backend == LEO2_BACKEND_AVX2 ||
             input.backend == LEO2_BACKEND_AVX512))
        {
            selection.path = kDecodePathTiled;
            selection.rule = kDecodeRuleMeasuredBatchTiled;
            selection.required_work_slots = tiled_work_slots;
        }
        else
        {
            selection.path = kDecodePathMaterialized;
            selection.rule = kDecodeRuleMeasuredMaterialized;
            selection.required_work_slots = materialized_work_slots;
        }
        return true;
    }

    if (input.profile != LEO2_PROFILE_LOW_V1 &&
        input.profile != LEO2_PROFILE_LEGACY_HIGH_V1)
    {
        selection.path = kDecodePathMaterialized;
        selection.rule = kDecodeRuleUnsupportedProfile;
        selection.required_work_slots = materialized_work_slots;
        return true;
    }
    if (tiled_work_slots < materialized_work_slots)
    {
        selection.path = kDecodePathTiled;
        selection.rule = kDecodeRuleWorkspaceTiled;
        selection.required_work_slots = tiled_work_slots;
    }
    else
    {
        selection.path = kDecodePathMaterialized;
        selection.rule = kDecodeRuleWorkspaceMaterialized;
        selection.required_work_slots = materialized_work_slots;
    }
    return true;
}

static inline const char* DecodePathName(DecodePath path)
{
    switch (path)
    {
    case kDecodePathNoOp: return "no_op";
    case kDecodePathDirect: return "direct";
    case kDecodePathGeneric: return "generic";
    case kDecodePathMaterialized: return "materialized";
    case kDecodePathTiled: return "tiled";
    }
    return "unknown";
}

static inline const char* DecodePathRuleName(DecodePathRule rule)
{
    switch (rule)
    {
    case kDecodeRuleNoOp: return "no_op";
    case kDecodeRuleDirect: return "direct";
    case kDecodeRuleForcedGeneric: return "forced_generic";
    case kDecodeRuleForcedMaterialized: return "forced_materialized";
    case kDecodeRuleForcedTiled: return "forced_tiled";
    case kDecodeRuleBalancedGeneric: return "balanced_generic";
    case kDecodeRuleMeasuredBatchTiled: return "measured_batch_tiled";
    case kDecodeRuleMeasuredMaterialized: return "measured_materialized";
    case kDecodeRuleWorkspaceTiled: return "workspace_tiled";
    case kDecodeRuleWorkspaceMaterialized: return "workspace_materialized";
    case kDecodeRuleTranslatedLowDiagnostic:
        return "translated_low_diagnostic";
    case kDecodeRuleUnsupportedProfile: return "unsupported_profile";
    }
    return "unknown";
}

/* Private build-tree introspection, intentionally absent from leopard2.h. */
leo2_result GetDecodePlanPathInfo(
    const leo2_decode_plan* plan,
    uint64_t shard_bytes,
    bool multi_item_batch,
    DecodePathInfo* info_out);

leo2_result GetDecodeCodecScratchPathInfo(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    DecodePathInfo* info_out);

} // namespace leopard2_internal
