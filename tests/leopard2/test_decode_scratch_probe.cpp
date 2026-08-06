#include "leopard2.h"
#include "Leopard2Backend.h"
#include "Leopard2Dispatch.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct ProbeCase
{
    const char* name;
    uint32_t k;
    uint32_t r;
    leo2_profile profile;
    leo2_field field;
    uint32_t flags;
    uint64_t shard_bytes;
    uint32_t losses;
    bool auto_observed;
};

void require_result(leo2_result actual, const char* operation)
{
    if (actual != LEO2_SUCCESS)
        throw std::runtime_error(
            std::string(operation) + ": " + leo2_result_string(actual));
}

void require_result(
    leo2_result actual,
    leo2_result expected,
    const char* operation)
{
    if (actual != expected)
        throw std::runtime_error(
            std::string(operation) + ": expected " +
            leo2_result_string(expected) + ", got " +
            leo2_result_string(actual));
}

size_t align_up(size_t value, size_t alignment)
{
    return (value + alignment - 1) & ~(alignment - 1);
}

size_t transform_scratch_bytes(
    const ProbeCase& test,
    uint32_t parent,
    size_t work_slots,
    size_t work_slot_bytes,
    size_t input_slot_bytes = 64)
{
    const size_t range_count = static_cast<size_t>(test.k) * 2 + test.r;
    const size_t ranges = range_count * sizeof(uintptr_t) * 2;
    const size_t pointer_count = static_cast<size_t>(parent) + work_slots;
    const size_t pointer_offset = align_up(ranges, alignof(void*));
    const size_t data_offset = align_up(
        pointer_offset + pointer_count * sizeof(void*), 64);
    const size_t tail = static_cast<size_t>(test.shard_bytes) & 63;
    const size_t input_bytes = tail == 0
        ? 0
        : static_cast<size_t>(test.k) * input_slot_bytes;
    return data_offset + input_bytes + work_slots * work_slot_bytes;
}

size_t direct_scratch_bytes(const ProbeCase& test)
{
    const size_t range_count =
        test.profile == LEO2_PROFILE_LEGACY_HIGH_V1 &&
                test.k == 1 && test.r == 1
            ? 0
            : static_cast<size_t>(test.k) * 2 + test.r;
    return align_up(range_count * sizeof(uintptr_t) * 2, 64);
}

size_t observe_work_slot_bytes(
    const ProbeCase& test,
    uint32_t parent,
    size_t scratch_bytes,
    size_t work_slots,
    bool codec_query,
    size_t input_slot_bytes = 64)
{
    if (work_slots == 0)
        throw std::runtime_error("cannot observe zero transform work slots");
    const size_t base = transform_scratch_bytes(
        test, parent, work_slots, 0, input_slot_bytes);
    if (scratch_bytes < base ||
        (scratch_bytes - base) % work_slots != 0)
        throw std::runtime_error(codec_query
            ? "codec scratch query has a non-integral work-slot width"
            : "plan scratch query has a non-integral work-slot width");
    const size_t bytes = (scratch_bytes - base) / work_slots;
    if (bytes == 0 || (bytes & 63U) != 0)
        throw std::runtime_error(codec_query
            ? "codec scratch query has an invalid work-slot width"
            : "plan scratch query has an invalid work-slot width");
    return bytes;
}

const char* backend_name(leo2_backend backend)
{
    switch (backend)
    {
    case LEO2_BACKEND_SCALAR: return "scalar";
    case LEO2_BACKEND_SSSE3: return "ssse3";
    case LEO2_BACKEND_AVX2: return "avx2";
    case LEO2_BACKEND_NEON: return "neon";
    case LEO2_BACKEND_AVX512: return "avx512";
    case LEO2_BACKEND_GFNI: return "gfni";
    case LEO2_BACKEND_AUTO: return "auto";
    default: return "unknown";
    }
}

const char* profile_name(leo2_profile profile)
{
    return profile == LEO2_PROFILE_LOW_V1 ? "low" : "high";
}

const char* field_name(leo2_field field)
{
    return field == LEO2_FIELD_GF16 ? "gf16" : "gf8";
}

void emit_case(
    leo2_context* context,
    const ProbeCase& test,
    uint64_t detected_l3_bytes,
    bool multi_item_batch = false,
    bool fused_ragged = false,
    bool context_auto_requested = true)
{
    leo2_codec_options options;
    std::memset(&options, 0, sizeof(options));
    options.struct_size = sizeof(options);
    options.flags = test.flags;

    leo2_codec* codec = NULL;
    require_result(leo2_codec_create(
        context, test.k, test.r, test.profile, test.field, &options, &codec),
        "codec create");
    const uint32_t parent = leo2_codec_parent_count(codec);
    const uint32_t padded = leo2_codec_padded_side(codec);

    std::vector<uint8_t> originals(test.k, 1);
    std::vector<uint8_t> recovery(test.r, 1);
    for (uint32_t i = 0; i < test.losses; ++i)
        originals[i] = 0;

    leo2_decode_plan* plan = NULL;
    require_result(leo2_decode_plan_create(
        codec, originals.data(), recovery.data(), &plan), "plan create");
    size_t plan_scratch = 0;
    size_t codec_scratch = 0;
    require_result(leo2_decode_plan_scratch_size(
        plan, test.shard_bytes, &plan_scratch), "plan scratch query");
    require_result(leo2_decode_scratch_size(
        codec, test.shard_bytes, &codec_scratch), "codec scratch query");

    const leo2_backend backend = leo2_context_backend(context);
    const size_t direct_bytes = direct_scratch_bytes(test);
    leopard2_internal::DecodePathInfo selected;
    leopard2_internal::DecodePathInfo scratch_selection;
    leopard2_internal::DecodePathInfo codec_selection;
    require_result(leopard2_internal::GetDecodePlanPathInfo(
        plan, test.shard_bytes, multi_item_batch, &selected),
        "selected path introspection");
    require_result(leopard2_internal::GetDecodePlanPathInfo(
        plan, test.shard_bytes, false, &scratch_selection),
        "scratch path introspection");
    require_result(leopard2_internal::GetDecodeCodecScratchPathInfo(
        codec, test.shard_bytes, &codec_selection),
        "codec path introspection");
    const bool no_op = selected.path == leopard2_internal::kDecodePathNoOp;
    const bool direct = selected.path == leopard2_internal::kDecodePathDirect;
    const size_t input_slot_bytes = fused_ragged
        ? align_up(static_cast<size_t>(test.shard_bytes), 64)
        : 64;
    if (direct != (!no_op && plan_scratch == direct_bytes))
        throw std::runtime_error("direct introspection/scratch mismatch");
    const size_t plan_work_slots = direct || no_op
        ? 0 : scratch_selection.required_work_slots;
    const size_t plan_work_slot_bytes = direct || no_op
        ? 0
        : observe_work_slot_bytes(
            test, parent, plan_scratch, plan_work_slots, false,
            input_slot_bytes);
    const bool k1_legacy_high =
        test.profile == LEO2_PROFILE_LEGACY_HIGH_V1 &&
        test.k == 1 && test.r == 1;
    const size_t codec_work_slots = k1_legacy_high
        ? 0 : codec_selection.required_work_slots;
    const size_t codec_work_slot_bytes = k1_legacy_high
        ? 0
        : observe_work_slot_bytes(
            test, parent, codec_scratch, codec_work_slots, true,
            input_slot_bytes);
    size_t plan_execution_tile_count = 0;
    size_t plan_execution_tile_bytes = 0;
    require_result(leopard2_internal::GetDecodePlanExecutionTiles(
        plan, test.shard_bytes, &plan_execution_tile_count,
        &plan_execution_tile_bytes), "plan execution-tile query");
    if (direct || no_op)
    {
        if (plan_execution_tile_count != 0 ||
            plan_execution_tile_bytes != 0)
            throw std::runtime_error(
                "terminal plan reported transform execution tiles");
    }
    else
    {
        const size_t expected_slot_bytes = std::max(
            plan_execution_tile_bytes,
            (static_cast<size_t>(test.shard_bytes) & 63U) ? size_t(64) : 0);
        if (expected_slot_bytes != plan_work_slot_bytes)
            throw std::runtime_error(
                "plan execution-tile/work-slot mismatch");
    }
    size_t codec_execution_tile_count = 0;
    size_t codec_execution_tile_bytes = 0;
    require_result(leopard2_internal::GetDecodeCodecExecutionTiles(
        codec, test.shard_bytes, &codec_execution_tile_count,
        &codec_execution_tile_bytes), "codec execution-tile query");
    if (k1_legacy_high)
    {
        if (codec_execution_tile_count != 0 ||
            codec_execution_tile_bytes != 0)
            throw std::runtime_error(
                "terminal codec reported transform execution tiles");
    }
    else
    {
        const size_t expected_slot_bytes = std::max(
            codec_execution_tile_bytes,
            (static_cast<size_t>(test.shard_bytes) & 63U) ? size_t(64) : 0);
        if (expected_slot_bytes != codec_work_slot_bytes)
            throw std::runtime_error(
                "codec execution-tile/work-slot mismatch");
    }
    uint64_t effective_l3_bytes = 0;
    uint64_t live_set_target_bytes = 0;
    uint64_t tile_threshold_bytes = 0;
    require_result(leopard2_internal::GetContextGF16CachePolicy(
        context, &effective_l3_bytes, &live_set_target_bytes,
        &tile_threshold_bytes), "context GF16 cache-policy query");

    std::cout
        << "{\"schema\":\"leopard2-decode-scratch-probe/v3\""
        << ",\"name\":\"" << test.name
        << "\",\"K\":" << test.k
        << ",\"R\":" << test.r
        << ",\"profile\":\"" << profile_name(test.profile)
        << "\",\"field\":\"" << field_name(test.field)
        << "\",\"flags\":" << test.flags
        << ",\"shard_bytes\":" << test.shard_bytes
        << ",\"losses\":" << test.losses
        << ",\"parent\":" << parent
        << ",\"padded\":" << padded
        << ",\"pointer_bytes\":" << sizeof(void*)
        << ",\"detected_l3_bytes\":" << detected_l3_bytes
        << ",\"effective_l3_bytes\":" << effective_l3_bytes
        << ",\"live_set_target_bytes\":" << live_set_target_bytes
        << ",\"tile_threshold_bytes\":" << tile_threshold_bytes
        << ",\"backend\":\"" << backend_name(backend)
        << "\",\"selected_path\":\""
        << leopard2_internal::DecodePathName(selected.path)
        << "\",\"selected_rule\":\""
        << leopard2_internal::DecodePathRuleName(selected.rule)
        << "\",\"selected_matching_auto_rules\":"
        << selected.matching_auto_rules
        << ",\"selected_required_work_slots\":"
        << selected.required_work_slots
        << ",\"scratch_path\":\""
        << leopard2_internal::DecodePathName(scratch_selection.path)
        << "\",\"scratch_rule\":\""
        << leopard2_internal::DecodePathRuleName(scratch_selection.rule)
        << "\",\"codec_path\":\""
        << leopard2_internal::DecodePathName(codec_selection.path)
        << "\",\"codec_rule\":\""
        << leopard2_internal::DecodePathRuleName(codec_selection.rule)
        << "\",\"codec_matching_auto_rules\":"
        << codec_selection.matching_auto_rules
        << ",\"multi_item_batch\":"
        << (multi_item_batch ? "true" : "false")
        << ",\"aligned_prefix_bytes\":" << selected.aligned_prefix_bytes
        << ",\"tail_bytes\":" << selected.tail_bytes
        << ",\"rounded_bytes\":" << selected.rounded_shard_bytes
        << ",\"context_auto_requested\":"
        << (context_auto_requested ? "true" : "false")
        << ",\"auto_host_backend_observed\":"
        << (test.auto_observed ? "true" : "false")
        << ",\"plan_work_slots\":" << plan_work_slots
        << ",\"plan_work_slot_bytes\":" << plan_work_slot_bytes
        << ",\"plan_execution_tile_count\":"
        << plan_execution_tile_count
        << ",\"plan_execution_tile_bytes\":"
        << plan_execution_tile_bytes
        << ",\"codec_work_slots\":" << codec_work_slots
        << ",\"codec_work_slot_bytes\":" << codec_work_slot_bytes
        << ",\"codec_execution_tile_count\":"
        << codec_execution_tile_count
        << ",\"codec_execution_tile_bytes\":"
        << codec_execution_tile_bytes
        << ",\"plan_scratch_bytes\":" << plan_scratch
        << ",\"codec_scratch_bytes\":" << codec_scratch
        << "}\n";

    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
}

void verify_query_boundaries(leo2_context* context, leo2_field field)
{
    ProbeCase test = {
        "abi_boundary_control", 1, 1, LEO2_PROFILE_LOW_V1,
        field, 0, 64, 1, false
    };
    leo2_codec* codec = NULL;
    require_result(leo2_codec_create(
        context, test.k, test.r, test.profile, test.field, NULL, &codec),
        "boundary codec create");
    const uint32_t parent = leo2_codec_parent_count(codec);
    const uint32_t padded = leo2_codec_padded_side(codec);
    if (parent != 2 || padded != 1)
        throw std::runtime_error("unexpected boundary-control geometry");

    const size_t maximum = std::numeric_limits<size_t>::max();
    const size_t work_slots = 2;
    const size_t data_offset = transform_scratch_bytes(
        test, parent, work_slots, static_cast<size_t>(test.shard_bytes)) -
        work_slots * test.shard_bytes;
    const size_t largest_layout_shard =
        ((maximum - data_offset) / work_slots) & ~static_cast<size_t>(63);
    test.shard_bytes = static_cast<uint64_t>(largest_layout_shard);
    size_t scratch_bytes = 0;
    require_result(leo2_decode_scratch_size(
        codec, test.shard_bytes, &scratch_bytes), "boundary codec scratch");
    if (scratch_bytes != transform_scratch_bytes(
            test, parent, work_slots,
            static_cast<size_t>(test.shard_bytes)))
        throw std::runtime_error("boundary codec scratch byte mismatch");

    require_result(leo2_decode_scratch_size(
        codec, static_cast<uint64_t>(largest_layout_shard) + 64,
        &scratch_bytes), LEO2_INVALID_COUNTS,
        "overflow-neighbor codec scratch");
    require_result(leo2_decode_scratch_size(
        codec, static_cast<uint64_t>(maximum) - 63, &scratch_bytes),
        LEO2_INVALID_COUNTS, "largest-roundable codec scratch");
    require_result(leo2_decode_scratch_size(
        codec, static_cast<uint64_t>(maximum), &scratch_bytes),
        LEO2_INVALID_ARGUMENT, "SIZE_MAX codec scratch");
    require_result(leo2_decode_scratch_size(
        codec, UINT64_MAX, &scratch_bytes),
        LEO2_INVALID_ARGUMENT, "UINT64_MAX codec scratch");
    leo2_codec_destroy(codec);
}

} // namespace

int main()
{
    try
    {
        leo2_context_options options;
        const uint64_t detected_l3_bytes =
            leopard::backend::DetectConservativeL3Bytes();
        std::memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AVX2;
        options.thread_count = 1;
        leo2_context* context = NULL;
        const leo2_result avx2_result = leo2_context_create(&options, &context);
        if (avx2_result == LEO2_SUCCESS)
        {
            const ProbeCase fused_case = {
                "avx2_fused_ragged_65", 57, 29,
                LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
                0, 65, 2, false
            };
            emit_case(context, fused_case, detected_l3_bytes,
                false, true, false);
            leo2_context_destroy(context);
        }
        else if (avx2_result != LEO2_UNSUPPORTED)
            require_result(avx2_result, "explicit AVX2 context create");

        std::memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AUTO;
        options.thread_count = 1;
        context = NULL;
        require_result(leo2_context_create(&options, &context), "context create");

        const ProbeCase cases[] = {
            { "noop_ragged", 240, 16, LEO2_PROFILE_LEGACY_HIGH_V1,
              LEO2_FIELD_GF8, LEO2_CODEC_FORCE_TILED_DECODE, 65, 0, false },
            { "forced_high_materialized_aligned", 240, 16,
              LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
              LEO2_CODEC_FORCE_MATERIALIZED_DECODE, 128, 1, false },
            { "forced_high_materialized_ragged", 240, 16,
              LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
              LEO2_CODEC_FORCE_MATERIALIZED_DECODE, 129, 1, false },
            { "forced_high_tiled_aligned", 240, 16,
              LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
              LEO2_CODEC_FORCE_TILED_DECODE, 64, 1, false },
            { "forced_high_tiled_ragged", 240, 16,
              LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
              LEO2_CODEC_FORCE_TILED_DECODE, 65, 1, false },
            { "forced_gf8_high_tiled_side64", 192, 64,
              LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
              LEO2_CODEC_FORCE_SPECIALIZED_DECODE, 256 * 1024, 9, false },
            { "forced_generic_ragged", 240, 16,
              LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
              LEO2_CODEC_FORCE_GENERIC_DECODE, 129, 2, false },
            { "forced_low_materialized_ragged", 8, 248,
              LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8,
              LEO2_CODEC_FORCE_MATERIALIZED_DECODE, 65, 8, false },
            { "forced_low_tiled_aligned", 8, 248,
              LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8,
              LEO2_CODEC_FORCE_TILED_DECODE, 128, 8, false },
            { "forced_gf16_high_tiled_ragged", 240, 16,
              LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16,
              LEO2_CODEC_FORCE_TILED_DECODE, 66, 1, false },
            { "gf16_low_cache_512k", 64, 193,
              LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16,
              LEO2_CODEC_FORCE_SPECIALIZED_DECODE, 512 * 1024, 9, false },
            { "gf16_low_cache_768k", 64, 193,
              LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16,
              LEO2_CODEC_FORCE_SPECIALIZED_DECODE, 768 * 1024, 9, false },
            { "gf16_low_cache_1m", 64, 193,
              LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16,
              LEO2_CODEC_FORCE_SPECIALIZED_DECODE, 1024 * 1024, 9, false },
            { "gf16_low_side256", 200, 800,
              LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16,
              LEO2_CODEC_FORCE_SPECIALIZED_DECODE, 64 * 1024, 9, false },
            { "gf16_high_codec_bound", 300, 100,
              LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16,
              LEO2_CODEC_FORCE_SPECIALIZED_DECODE, 448 * 1024, 9, false },
            { "gf16_high_max_loss", 2000, 512,
              LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16,
              LEO2_CODEC_FORCE_SPECIALIZED_DECODE, 64 * 1024, 512, false },
            { "direct_xor", 9, 1, LEO2_PROFILE_LEGACY_HIGH_V1,
              LEO2_FIELD_GF8, 0, 65, 1, false },
            { "direct_k1_r1", 1, 1, LEO2_PROFILE_LEGACY_HIGH_V1,
              LEO2_FIELD_GF8, 0, 17, 1, false },
            { "direct_copy", 1, 8, LEO2_PROFILE_LOW_V1,
              LEO2_FIELD_GF8, 0, 65, 1, false },
            { "direct_repair", 16, 8, LEO2_PROFILE_LOW_V1,
              LEO2_FIELD_GF8, 0, 65, 4, false },
            { "auto_ordinary_high", 240, 16,
              LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 0, 4096, 8, true },
            { "auto_balanced_64", 128, 128,
              LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 0, 64, 128, true },
            { "auto_balanced_256", 128, 128,
              LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 0, 256, 128, true },
            { "auto_high_16k", 224, 32,
              LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 0, 16 * 1024, 8, true },
            { "auto_high_32k", 224, 32,
              LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 0, 32 * 1024, 8, true },
            { "auto_high_128k", 224, 32,
              LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 0, 128 * 1024, 8, true },
            { "auto_low", 8, 248, LEO2_PROFILE_LOW_V1,
              LEO2_FIELD_GF8, 0, 4096, 8, true },
        };
        const uint32_t fields = leo2_context_field_mask(context);
        verify_query_boundaries(context,
            (fields & LEO2_FIELD_MASK_GF8) != 0
                ? LEO2_FIELD_GF8 : LEO2_FIELD_GF16);
        for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i)
        {
            const uint32_t bit = cases[i].field == LEO2_FIELD_GF16
                ? LEO2_FIELD_MASK_GF16 : LEO2_FIELD_MASK_GF8;
            if ((fields & bit) != 0)
                emit_case(context, cases[i], detected_l3_bytes);
        }
        if ((fields & LEO2_FIELD_MASK_GF8) != 0)
        {
            const ProbeCase batch_case = {
                "auto_high_32k_batch", 224, 32,
                LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
                0, 32 * 1024, 8, true
            };
            emit_case(context, batch_case, detected_l3_bytes, true);
        }
        leo2_context_destroy(context);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "decode scratch probe failed: " << error.what() << '\n';
        return 1;
    }
}
