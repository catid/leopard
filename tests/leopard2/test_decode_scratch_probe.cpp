#include "leopard2.h"
#include "Leopard2Dispatch.h"

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
    size_t work_slots)
{
    const size_t range_count = static_cast<size_t>(test.k) * 2 + test.r;
    const size_t ranges = range_count * sizeof(uintptr_t) * 2;
    const size_t pointer_count = static_cast<size_t>(parent) + work_slots;
    const size_t pointer_offset = align_up(ranges, alignof(void*));
    const size_t data_offset = align_up(
        pointer_offset + pointer_count * sizeof(void*), 64);
    const size_t tail = static_cast<size_t>(test.shard_bytes) & 63;
    const size_t prefix = static_cast<size_t>(test.shard_bytes) - tail;
    const size_t input_bytes = tail == 0
        ? 0
        : (static_cast<size_t>(test.k) + test.r) * 64;
    const size_t work_slot_bytes = tail == 0
        ? static_cast<size_t>(test.shard_bytes)
        : (prefix > 64 ? prefix : 64);
    return data_offset + input_bytes + work_slots * work_slot_bytes;
}

size_t direct_scratch_bytes(const ProbeCase& test)
{
    const size_t range_count = static_cast<size_t>(test.k) * 2 + test.r;
    return align_up(range_count * sizeof(uintptr_t) * 2, 64);
}

size_t observe_work_slots(
    const ProbeCase& test,
    uint32_t parent,
    uint32_t padded,
    size_t scratch_bytes,
    bool codec_query)
{
    const size_t upper = static_cast<size_t>(parent) + 2 * padded + test.r + 1;
    size_t match = static_cast<size_t>(-1);
    for (size_t slots = 0; slots <= upper; ++slots)
    {
        if (transform_scratch_bytes(test, parent, slots) != scratch_bytes)
            continue;
        if (match != static_cast<size_t>(-1))
            throw std::runtime_error("scratch query has ambiguous work-slot count");
        match = slots;
    }
    if (match == static_cast<size_t>(-1))
        throw std::runtime_error(codec_query
            ? "codec scratch query does not match any production layout"
            : "plan scratch query does not match any production layout");
    return match;
}

const char* backend_name(leo2_backend backend)
{
    switch (backend)
    {
    case LEO2_BACKEND_SCALAR: return "scalar";
    case LEO2_BACKEND_SSSE3: return "ssse3";
    case LEO2_BACKEND_AVX2: return "avx2";
    case LEO2_BACKEND_NEON: return "neon";
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

const char* classify_transform_path(
    const ProbeCase& test,
    leo2_backend backend,
    uint32_t parent,
    uint32_t padded,
    size_t rounded_bytes,
    uint32_t output_slots,
    bool codec_query)
{
    if ((test.flags & LEO2_CODEC_FORCE_GENERIC_DECODE) != 0)
        return "generic";
    if ((test.flags & LEO2_CODEC_FORCE_MATERIALIZED_DECODE) != 0)
        return "materialized";
    if ((test.flags & LEO2_CODEC_FORCE_TILED_DECODE) != 0)
        return "tiled";

    if (!codec_query &&
        leopard2_internal::ShouldUseBalancedGenericDecode(
            test.profile, test.field, test.k, test.r, padded, parent,
            output_slots, rounded_bytes, backend))
        return "generic";

    const uint32_t materialized_losses = codec_query ? 1 : output_slots;
    if (leopard2_internal::ShouldUseMaterializedHighDecode(
            test.profile, test.field, test.k, test.r, padded, parent,
            materialized_losses, rounded_bytes, backend))
        return "materialized";

    const size_t tiled = static_cast<size_t>(padded) * 2 +
        (test.profile == LEO2_PROFILE_LEGACY_HIGH_V1
            ? (codec_query ? test.r : output_slots)
            : 0);
    return tiled < parent ? "tiled" : "materialized";
}

void emit_case(leo2_context* context, const ProbeCase& test)
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
    const size_t rounded = align_up(static_cast<size_t>(test.shard_bytes), 64);
    const size_t direct_bytes = direct_scratch_bytes(test);
    const bool no_op = test.losses == 0;
    const bool direct = !no_op && plan_scratch == direct_bytes;
    const char* selected_path = NULL;
    if (no_op)
        selected_path = "no_op";
    else if (direct)
    {
        if (test.profile == LEO2_PROFILE_LEGACY_HIGH_V1 && padded == 1)
            selected_path = "direct_xor";
        else if (test.profile == LEO2_PROFILE_LOW_V1 && padded == 1)
            selected_path = "direct_copy";
        else
            selected_path = "direct_repair";
    }
    else
    {
        selected_path = classify_transform_path(
            test, backend, parent, padded, rounded, test.losses, false);
    }
    const char* codec_path = classify_transform_path(
        test, backend, parent, padded, rounded, test.r, true);
    const size_t plan_work_slots = direct || no_op
        ? 0
        : observe_work_slots(test, parent, padded, plan_scratch, false);
    const size_t codec_work_slots = observe_work_slots(
        test, parent, padded, codec_scratch, true);

    std::cout
        << "{\"name\":\"" << test.name
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
        << ",\"backend\":\"" << backend_name(backend)
        << "\",\"selected_path\":\"" << selected_path
        << "\",\"codec_path\":\"" << codec_path
        << "\",\"auto_host_backend_observed\":"
        << (test.auto_observed ? "true" : "false")
        << ",\"plan_work_slots\":" << plan_work_slots
        << ",\"codec_work_slots\":" << codec_work_slots
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
        test, parent, work_slots) - work_slots * test.shard_bytes;
    const size_t largest_layout_shard =
        ((maximum - data_offset) / work_slots) & ~static_cast<size_t>(63);
    test.shard_bytes = static_cast<uint64_t>(largest_layout_shard);
    size_t scratch_bytes = 0;
    require_result(leo2_decode_scratch_size(
        codec, test.shard_bytes, &scratch_bytes), "boundary codec scratch");
    if (scratch_bytes != transform_scratch_bytes(test, parent, work_slots))
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
        std::memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AUTO;
        options.thread_count = 1;
        leo2_context* context = NULL;
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
            { "direct_xor", 9, 1, LEO2_PROFILE_LEGACY_HIGH_V1,
              LEO2_FIELD_GF8, 0, 65, 1, false },
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
                emit_case(context, cases[i]);
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
