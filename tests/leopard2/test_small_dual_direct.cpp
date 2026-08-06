/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in the project
    LICENSE file are met.
*/

#include "direct_oracle.h"
#include "direct_repair.h"
#include "Leopard2Direct.h"
#include "Leopard2Dispatch.h"
#include "leopard2.h"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <vector>

#ifndef LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT
#define LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT 1
#endif

#ifndef LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE
#define LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE 0
#endif

#ifndef LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS
#define LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS 1
#endif

#ifndef LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK
#define LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK 0
#endif

namespace {

using leopard2_test::BinaryField;
using leopard2_test::Matrix;
using leopard2_test::ProfileLayout;

size_t g_max_direct_retained_bytes = 0;

const uint8_t kOriginalGuard = 0x6d;
const uint8_t kParityGuard = 0x7e;
const uint8_t kRestoredGuard = 0xa5;

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void require_result(leo2_result result, const char* operation)
{
    if (result != LEO2_SUCCESS)
    {
        std::ostringstream stream;
        stream << operation << ": " << leo2_result_string(result)
               << " (" << static_cast<int>(result) << ")";
        throw std::runtime_error(stream.str());
    }
}

struct Scratch
{
    std::vector<uint8_t> storage;
    void* data;
    size_t bytes;

    explicit Scratch(size_t requested)
        : storage(requested == 0 ? 0 : requested + 63, 0)
        , data(NULL)
        , bytes(requested)
    {
        if (requested != 0)
        {
            const uintptr_t raw =
                reinterpret_cast<uintptr_t>(storage.data());
            data = reinterpret_cast<void*>((raw + 63U) & ~uintptr_t(63U));
        }
    }
};

struct Buffers
{
    size_t shard_bytes;
    std::vector<std::vector<uint8_t> > originals;
    std::vector<std::vector<uint8_t> > parities;
    std::vector<std::vector<uint8_t> > restored;
    std::vector<const void*> original_input;
    std::vector<const void*> recovery_input;
    std::vector<void*> restored_output;
    Scratch scratch;

    Buffers(unsigned k, unsigned r, unsigned losses, size_t bytes,
            size_t scratch_bytes, const BinaryField& field,
            const Matrix& generator, uint64_t seed)
        : shard_bytes(bytes)
        , originals(k, std::vector<uint8_t>(bytes + 2, kOriginalGuard))
        , parities(r, std::vector<uint8_t>(bytes + 2, kParityGuard))
        , restored(k, std::vector<uint8_t>(bytes + 2, kRestoredGuard))
        , original_input(k, NULL)
        , recovery_input(r, NULL)
        , restored_output(k, NULL)
        , scratch(scratch_bytes)
    {
        for (unsigned original = 0; original < k; ++original)
        {
            for (size_t byte = 0; byte < bytes; ++byte)
            {
                seed = seed * UINT64_C(6364136223846793005) +
                    UINT64_C(1442695040888963407);
                originals[original][byte + 1] = static_cast<uint8_t>(
                    seed >> 56);
            }
        }
        for (unsigned parity = 0; parity < r; ++parity)
        {
            for (size_t byte = 0; byte < bytes; ++byte)
            {
                uint8_t value = 0;
                for (unsigned original = 0; original < k; ++original)
                {
                    value ^= static_cast<uint8_t>(field.multiply(
                        generator[k + parity][original],
                        originals[original][byte + 1]));
                }
                parities[parity][byte + 1] = value;
            }
        }
        for (unsigned original = 0; original < k; ++original)
        {
            if (original < losses)
                restored_output[original] = restored[original].data() + 1;
            else
                original_input[original] = originals[original].data() + 1;
        }
        for (unsigned parity = 0; parity < r; ++parity)
            recovery_input[parity] = parities[parity].data() + 1;
    }

    void clear_restored(unsigned losses)
    {
        for (unsigned original = 0; original < losses; ++original)
            std::fill(restored[original].begin(), restored[original].end(),
                kRestoredGuard);
    }

    void apply_presence(
        const std::vector<uint8_t>& original_present,
        const std::vector<uint8_t>& recovery_present)
    {
        require(original_present.size() == originals.size(),
            "original presence size mismatch");
        require(recovery_present.size() == parities.size(),
            "recovery presence size mismatch");
        std::fill(original_input.begin(), original_input.end(),
            static_cast<const void*>(NULL));
        std::fill(recovery_input.begin(), recovery_input.end(),
            static_cast<const void*>(NULL));
        std::fill(restored_output.begin(), restored_output.end(),
            static_cast<void*>(NULL));
        for (size_t original = 0; original < originals.size(); ++original)
        {
            if (original_present[original] != 0)
                original_input[original] = originals[original].data() + 1;
            else
                restored_output[original] = restored[original].data() + 1;
        }
        for (size_t parity = 0; parity < parities.size(); ++parity)
        {
            if (recovery_present[parity] != 0)
                recovery_input[parity] = parities[parity].data() + 1;
        }
    }

    void verify_guards() const
    {
        for (size_t original = 0; original < originals.size(); ++original)
        {
            require(originals[original].front() == kOriginalGuard &&
                    originals[original].back() == kOriginalGuard,
                "original input guard byte was modified");
            require(restored[original].front() == kRestoredGuard &&
                    restored[original].back() == kRestoredGuard,
                "restored output guard byte was modified");
        }
        for (size_t parity = 0; parity < parities.size(); ++parity)
        {
            require(parities[parity].front() == kParityGuard &&
                    parities[parity].back() == kParityGuard,
                "recovery input guard byte was modified");
        }
    }

    void verify(unsigned losses) const
    {
        for (unsigned original = 0; original < losses; ++original)
        {
            require(std::memcmp(restored[original].data() + 1,
                        originals[original].data() + 1, shard_bytes) == 0,
                "recovered original differs from independent oracle input");
        }
        verify_guards();
    }

    void verify_missing(const std::vector<unsigned>& missing) const
    {
        for (size_t missing_i = 0; missing_i < missing.size(); ++missing_i)
        {
            const unsigned original = missing[missing_i];
            require(original < originals.size(),
                "missing original index is out of range");
            require(std::memcmp(restored[original].data() + 1,
                        originals[original].data() + 1, shard_bytes) == 0,
                "non-prefix recovered original differs from independent oracle input");
        }
        verify_guards();
    }

    leo2_decode_batch_item item()
    {
        leo2_decode_batch_item result = {};
        result.shard_bytes = shard_bytes;
        result.original = original_input.data();
        result.recovery = recovery_input.data();
        result.restored_original = restored_output.data();
        result.scratch = scratch.data;
        result.scratch_bytes = scratch.bytes;
        return result;
    }
};

leo2_codec* make_codec(
    leo2_context* context,
    unsigned k,
    unsigned r,
    uint32_t flags)
{
    leo2_codec_options options = {};
    options.struct_size = sizeof(options);
    options.flags = flags;
    leo2_codec* codec = NULL;
    require_result(leo2_codec_create(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, &options, &codec),
        "codec create");
    return codec;
}

leo2_decode_plan* make_plan(
    const leo2_codec* codec,
    const std::vector<uint8_t>& original_present,
    const std::vector<uint8_t>& recovery_present)
{
    leo2_decode_plan* plan = NULL;
    require_result(leo2_decode_plan_create(codec, original_present.data(),
        recovery_present.data(), &plan), "decode plan create");
    return plan;
}

size_t direct_retained_bytes(
    const leo2_decode_plan* plan,
    bool expect_direct_capability,
    const char* operation)
{
    size_t bytes = 0;
    size_t term_count = 0;
    size_t source_row_count = 0;
    require(leopard2_internal::GetDecodePlanDirectStorageInfo(
            plan, &bytes, &term_count, &source_row_count),
        operation);
    require((term_count != 0) == expect_direct_capability,
        "direct retained-term count disagrees with plan capability");
    require((bytes != 0) == expect_direct_capability,
        "direct retained bytes disagree with plan capability");
    require(bytes <= 4096,
        "small dual plan retained more than 4096 direct bytes");
    require(source_row_count <= 16,
        "small dual plan retained too many source-major rows");
    g_max_direct_retained_bytes = std::max(
        g_max_direct_retained_bytes, bytes);
    return bytes;
}

leo2_decode_plan* make_plan(
    const leo2_codec* codec,
    unsigned k,
    unsigned r,
    unsigned losses)
{
    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 1);
    for (unsigned original = 0; original < losses; ++original)
        original_present[original] = 0;
    return make_plan(codec, original_present, recovery_present);
}

void verify_transform_only(
    leo2_context* context,
    unsigned k,
    unsigned r,
    unsigned losses,
    leo2_profile profile,
    leo2_field field)
{
    leo2_codec_options options = {};
    options.struct_size = sizeof(options);
    leo2_codec* codec = NULL;
    require_result(leo2_codec_create(context, k, r, profile, field,
        &options, &codec), "negative codec create");
    leo2_decode_plan* plan = make_plan(codec, k, r, losses);
    leopard2_internal::DecodePathInfo path = {};
    require_result(leopard2_internal::GetDecodePlanPathInfo(
        plan, 4096, false, &path), "negative path introspection");
    require(path.path != leopard2_internal::kDecodePathDirect &&
            path.direct_executor ==
                leopard2_internal::kDirectRepairExecutorNone,
        "ineligible codec selected the small dual direct path");
    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
}

#if LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE == 0 && \
    LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT && \
    LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS
void verify_locator_coefficient_map(
    leo2_context* context,
    unsigned k,
    unsigned r,
    const std::vector<uint8_t>& original_present,
    const std::vector<uint8_t>& recovery_present)
{
    require(original_present.size() == k && recovery_present.size() == r,
        "coefficient-map presence size mismatch");
    std::vector<unsigned> missing;
    for (unsigned original = 0; original < k; ++original)
        if (!original_present[original])
            missing.push_back(original);
    require(missing.size() >= 5 && missing.size() <= 8,
        "coefficient-map loss count is outside the dual region");

    const BinaryField field = leopard2_test::make_legacy_gf8();
    const ProfileLayout layout = leopard2_test::make_profile_layout(
        leopard2_test::kLegacyHigh, k, r);
    const Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);
    const leopard2_test::DirectRepairPlan oracle =
        leopard2_test::make_direct_repair_plan(
            field, generator, missing, recovery_present);

    leo2_codec* codec = make_codec(context, k, r, 0);
    leo2_decode_plan* plan = make_plan(
        codec, original_present, recovery_present);
    require(leopard2_internal::DecodePlanUsesLocatorDirectTerms(plan),
        "coefficient-map plan did not use locator-derived terms");

    static const size_t kShardBytes = 1024;
    size_t scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(
        plan, kShardBytes, &scratch_bytes),
        "coefficient-map scratch query");
    Scratch scratch(scratch_bytes);
    std::vector<std::vector<uint8_t> > originals(
        k, std::vector<uint8_t>(kShardBytes, 0));
    std::vector<std::vector<uint8_t> > recovery(
        r, std::vector<uint8_t>(kShardBytes, 0));
    std::vector<std::vector<uint8_t> > restored(
        k, std::vector<uint8_t>(kShardBytes, 0));
    std::vector<const void*> original_input(k, NULL);
    std::vector<const void*> recovery_input(r, NULL);
    std::vector<void*> restored_output(k, NULL);
    for (unsigned original = 0; original < k; ++original)
    {
        if (original_present[original])
            original_input[original] = originals[original].data();
        else
            restored_output[original] = restored[original].data();
    }
    for (unsigned parity = 0; parity < r; ++parity)
        if (recovery_present[parity])
            recovery_input[parity] = recovery[parity].data();

    struct BasisSource
    {
        bool parity;
        unsigned index;
    };
    std::vector<BasisSource> sources;
    for (unsigned parity_i = 0;
         parity_i < oracle.selected_parities.size(); ++parity_i)
    {
        const BasisSource source = {
            true, oracle.selected_parities[parity_i]
        };
        sources.push_back(source);
    }
    for (unsigned original = 0; original < k; ++original)
    {
        if (original_present[original])
        {
            const BasisSource source = { false, original };
            sources.push_back(source);
        }
    }
    require(sources.size() == k,
        "coefficient-map selected source count is not K");

    for (size_t source_i = 0; source_i < sources.size(); ++source_i)
    {
        for (unsigned original = 0; original < k; ++original)
            std::fill(originals[original].begin(), originals[original].end(), 0);
        for (unsigned parity = 0; parity < r; ++parity)
            std::fill(recovery[parity].begin(), recovery[parity].end(), 0);
        for (unsigned original = 0; original < k; ++original)
            std::fill(restored[original].begin(), restored[original].end(), 0xa5);
        const BasisSource& source = sources[source_i];
        (source.parity ? recovery[source.index] : originals[source.index])[0] = 1;

        require_result(leo2_decode_plan_execute(plan, kShardBytes,
            original_input.data(), recovery_input.data(),
            restored_output.data(), scratch.data, scratch.bytes),
            "coefficient-map plan execute");
        for (size_t output = 0; output < missing.size(); ++output)
        {
            uint8_t expected = 0;
            const std::vector<leopard2_test::DirectRepairTerm>& terms =
                oracle.output_terms[output];
            for (size_t term_i = 0; term_i < terms.size(); ++term_i)
            {
                if (terms[term_i].parity == source.parity &&
                    terms[term_i].index == source.index)
                {
                    expected = static_cast<uint8_t>(terms[term_i].coefficient);
                    break;
                }
            }
            const std::vector<uint8_t>& output_shard =
                restored[missing[output]];
            require(output_shard[0] == expected,
                "locator-derived coefficient differs from independent matrix oracle");
            require(std::find_if(output_shard.begin() + 1,
                        output_shard.end(), [](uint8_t value) {
                            return value != 0;
                        }) == output_shard.end(),
                "coefficient-map zero tail was modified incorrectly");
        }
    }

    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
}

void verify_locator_pattern_sweep(leo2_context* context)
{
    size_t plan_count = 0;
    for (unsigned k = 5; k <= 16; ++k)
    {
        /* Explicit LEGACY_HIGH is public even when R > K.  AUTO would choose
           LOW for those counts, but the high-profile dual-plan algebra and
           xor-P execution view still require direct coverage. */
        for (unsigned r = 5; r <= 8; ++r)
        {
            leo2_codec* codec = make_codec(context, k, r, 0);
            for (unsigned losses = 5;
                 losses <= std::min(k, r); ++losses)
            {
                if (losses == k && k >= 7)
                    continue;
                for (unsigned pattern = 0; pattern < 2; ++pattern)
                {
                    std::vector<uint8_t> original_present(k, 1);
                    for (unsigned i = 0; i < losses; ++i)
                    {
                        const unsigned original = pattern == 0
                            ? i : k - 1U - i;
                        original_present[original] = 0;
                    }
                    std::vector<uint8_t> recovery_present(r, 1);
                    if (pattern != 0)
                    {
                        for (unsigned parity = 0;
                             parity < r - losses; ++parity)
                            recovery_present[parity] = 0;
                    }
                    leo2_decode_plan* plan = make_plan(
                        codec, original_present, recovery_present);
                    require(leopard2_internal::
                            DecodePlanUsesLocatorDirectTerms(plan),
                        "eligible sweep plan did not use locator-derived terms");
                    leo2_decode_plan_destroy(plan);
                    ++plan_count;
                }
            }
            leo2_codec_destroy(codec);
        }
    }
    require(plan_count == 214,
        "locator coefficient sweep did not cover the expected plan count");
}
#endif

size_t verify_raw_native_high_oracle_matrix(leo2_context* context)
{
    require(leopard2_internal::SetOneShotPlanSetupModeForDiagnostics(3),
        "raw native-high setup mode select failed");
    const BinaryField field = leopard2_test::make_legacy_gf8();
    const unsigned original_counts[] = { 9, 12, 16 };
    std::vector<size_t> byte_sizes;
    for (size_t bytes = 1; bytes <= 256; ++bytes)
        byte_sizes.push_back(bytes);
    const size_t extended_bytes[] = {
        257,
        319, 320, 321,
        511, 512, 513,
        1023, 1024, 1025,
        2047, 2048, 2049,
        4095, 4096, 4097,
        6143, 6144, 6145,
        7103, 7104, 7105,
        7135, 7136, 7137,
        7167, 7168
    };
    byte_sizes.insert(byte_sizes.end(), extended_bytes,
        extended_bytes +
            sizeof(extended_bytes) / sizeof(extended_bytes[0]));
    const size_t source_major_and_fallback_bytes[] = {
        12287, 12288, 12289,
        12319, 12320, 12321,
        12350, 12351,
        16383, 16384, 16385,
        32768, 65536
    };
    byte_sizes.insert(byte_sizes.end(), source_major_and_fallback_bytes,
        source_major_and_fallback_bytes +
            sizeof(source_major_and_fallback_bytes) /
                sizeof(source_major_and_fallback_bytes[0]));
    size_t cases = 0;
    for (size_t k_i = 0;
         k_i < sizeof(original_counts) / sizeof(original_counts[0]); ++k_i)
    {
        const unsigned k = original_counts[k_i];
        for (unsigned r = 5; r <= 8; ++r)
        {
            const ProfileLayout layout = leopard2_test::make_profile_layout(
                leopard2_test::kLegacyHigh, k, r);
            const Matrix generator =
                leopard2_test::direct_systematic_generator(field, layout);
            leo2_codec* codec = make_codec(context, k, r, 0);
            leo2_codec* reusable_codec = make_codec(
                context, k, r, LEO2_CODEC_FORCE_SPECIALIZED_DECODE);
            const unsigned loss_counts[] = { 4, 5, r };
            for (size_t loss_i = 0;
                 loss_i < sizeof(loss_counts) / sizeof(loss_counts[0]);
                 ++loss_i)
            {
                const unsigned losses = loss_counts[loss_i];
                bool duplicate_loss_count = false;
                for (size_t previous = 0; previous < loss_i; ++previous)
                    duplicate_loss_count |= losses == loss_counts[previous];
                if (duplicate_loss_count)
                    continue;
                for (unsigned shape = 0; shape < 2; ++shape)
                {
                    std::vector<uint8_t> original_present(k, 1);
                    std::vector<unsigned> missing;
                    for (unsigned i = 0; i < losses; ++i)
                    {
                        const unsigned original = shape == 0
                            ? i : k - 1U - i;
                        original_present[original] = 0;
                        missing.push_back(original);
                    }
                    std::sort(missing.begin(), missing.end());
                    std::vector<uint8_t> recovery_present(r, 1);
                    if (shape != 0 && losses < r)
                    {
                        std::fill(recovery_present.begin(),
                            recovery_present.end(), 0);
                        for (unsigned i = 0; i < losses; ++i)
                            recovery_present[r - 1U - i] = 1;
                    }
                    leo2_decode_plan* reusable_plan = make_plan(
                        reusable_codec, original_present, recovery_present);

                    for (size_t byte_index = 0;
                         byte_index < byte_sizes.size(); ++byte_index)
                    {
                        const size_t bytes = byte_sizes[byte_index];
                        size_t raw_scratch_bytes = 0;
                        size_t reusable_scratch_bytes = 0;
                        require_result(leo2_decode_scratch_size(
                            codec, bytes, &raw_scratch_bytes),
                            "raw native-high scratch query");
                        require_result(leo2_decode_plan_scratch_size(
                            reusable_plan, bytes, &reusable_scratch_bytes),
                            "raw native-high reusable scratch query");
                        const uint64_t seed =
                            UINT64_C(0x5241574849474800) +
                            k * 65537U + r * 257U + losses * 17U +
                            shape * 7U + bytes;
                        Buffers raw(k, r, 0, bytes, raw_scratch_bytes,
                            field, generator, seed);
                        Buffers reusable(k, r, 0, bytes,
                            reusable_scratch_bytes, field, generator, seed);
                        raw.apply_presence(
                            original_present, recovery_present);
                        reusable.apply_presence(
                            original_present, recovery_present);
                        require_result(leo2_decode(codec, bytes,
                            original_present.data(), recovery_present.data(),
                            raw.original_input.data(),
                            raw.recovery_input.data(),
                            raw.restored_output.data(), raw.scratch.data,
                            raw.scratch.bytes),
                            "raw native-high one-shot decode");
                        require_result(leo2_decode_plan_execute(
                            reusable_plan, bytes,
                            reusable.original_input.data(),
                            reusable.recovery_input.data(),
                            reusable.restored_output.data(),
                            reusable.scratch.data, reusable.scratch.bytes),
                            "raw native-high reusable Algorithm 5 decode");
                        raw.verify_missing(missing);
                        reusable.verify_missing(missing);
                        for (size_t missing_i = 0;
                             missing_i < missing.size(); ++missing_i)
                        {
                            const unsigned original = missing[missing_i];
                            require(std::memcmp(
                                    raw.restored[original].data() + 1,
                                    reusable.restored[original].data() + 1,
                                    bytes) == 0,
                                "raw native-high output differs from reusable Algorithm 5");
                        }
                        ++cases;
                    }
                    leo2_decode_plan_destroy(reusable_plan);
                }
            }
            leo2_codec_destroy(reusable_codec);
            leo2_codec_destroy(codec);
        }
    }
    return cases;
}

void verify_routes_and_execution(
    leo2_context* context,
    unsigned k,
    unsigned r,
    unsigned losses)
{
    const BinaryField field = leopard2_test::make_legacy_gf8();
    const ProfileLayout layout = leopard2_test::make_profile_layout(
        leopard2_test::kLegacyHigh, k, r);
    const Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);
    leo2_codec* codec = make_codec(context, k, r, 0);
    leo2_codec* transform_codec = make_codec(context, k, r,
        LEO2_CODEC_FORCE_SPECIALIZED_DECODE);
    leo2_decode_plan* plan = make_plan(codec, k, r, losses);
    leo2_decode_plan* transform_plan =
        make_plan(transform_codec, k, r, losses);
    const bool production_direct_shape =
        !(losses == k && k >= 7);
    const bool diagnostic_direct =
        LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE != 0;
    const bool expect_direct_capability = diagnostic_direct ||
        (LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT != 0 &&
         production_direct_shape);
    direct_retained_bytes(plan, expect_direct_capability,
        "dual direct storage introspection failed");
    direct_retained_bytes(transform_plan, false,
        "forced transform storage introspection failed");
    const bool expect_locator_terms =
        LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE == 0 &&
        LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT != 0 &&
        LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS != 0 &&
        production_direct_shape;
    require(leopard2_internal::DecodePlanUsesLocatorDirectTerms(plan) ==
            expect_locator_terms,
        "dual plan direct-term provenance disagrees with the selector");
    require(!leopard2_internal::DecodePlanUsesLocatorDirectTerms(
            transform_plan),
        "forced transform plan retained locator-derived direct terms");
    if (LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE == 0 &&
        LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT != 0 &&
        production_direct_shape &&
        LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK != 0 &&
        (leo2_test_decode_plan_uses_translated_low(plan) || losses == r))
    {
        leopard2_internal::DecodePlanPrunedScheduleInfo schedules = {};
        require(leopard2_internal::GetDecodePlanPrunedScheduleInfo(
                plan, &schedules),
            "dual pruned-schedule introspection failed");
        require(schedules.low_input_plan_count == 0 &&
                schedules.low_output_plan_count == 0 &&
                schedules.high_input_plan_count == 0 &&
                schedules.high_output_plan_count == 0,
            "regular-fallback dual plan retained a pruned schedule");
    }

    const size_t byte_counts[] = {
        1, 63, 64, 1023, 1024, 1025, 4095, 4096, 4097
    };
    for (size_t byte_i = 0;
         byte_i < sizeof(byte_counts) / sizeof(byte_counts[0]); ++byte_i)
    {
        const size_t bytes = byte_counts[byte_i];
        leopard2_internal::DecodePathInfo path = {};
        require_result(leopard2_internal::GetDecodePlanPathInfo(
            plan, bytes, false, &path), "path introspection");
        const bool expect_direct =
            diagnostic_direct ||
            (LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT != 0 &&
             production_direct_shape &&
             (bytes <= 63 || bytes >= 1024));
        const leopard2_internal::DirectRepairExecutor expected_executor =
            expect_direct
                ? (LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE == 1
                    ? leopard2_internal::kDirectRepairExecutorOutputMajor
                    : leopard2_internal::kDirectRepairExecutorSourceMajor)
                : leopard2_internal::kDirectRepairExecutorNone;
        if ((path.path == leopard2_internal::kDecodePathDirect) !=
            expect_direct)
        {
            std::ostringstream stream;
            stream << "dual plan selected the wrong side of the 1024-byte "
                   << "boundary: K=" << k << " R=" << r
                   << " losses=" << losses << " bytes=" << bytes
                   << " expected_direct=" << expect_direct
                   << " actual_path=" << static_cast<int>(path.path);
            throw std::runtime_error(stream.str());
        }
        require(path.direct_executor == expected_executor,
            "dual plan reported the wrong byte executor");

        size_t scratch_bytes = 0;
        size_t transform_scratch_bytes = 0;
        require_result(leo2_decode_plan_scratch_size(
            plan, bytes, &scratch_bytes), "dual scratch query");
        require_result(leo2_decode_plan_scratch_size(
            transform_plan, bytes, &transform_scratch_bytes),
            "transform scratch query");
        require(expect_direct
                ? scratch_bytes < transform_scratch_bytes
                : scratch_bytes == transform_scratch_bytes,
            "scratch query disagrees with the selected execution route");

        Buffers buffers(k, r, losses, bytes, scratch_bytes, field,
            generator, UINT64_C(0x4455414c00000000) + bytes + k * 257U);
        require_result(leo2_decode_plan_execute(plan, bytes,
            buffers.original_input.data(), buffers.recovery_input.data(),
            buffers.restored_output.data(), buffers.scratch.data,
            buffers.scratch.bytes), "dual plan execute");
        buffers.verify(losses);
    }

    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 1);
    for (unsigned original = 0; original < losses; ++original)
        original_present[original] = 0;
    const size_t one_shot_bytes[] = { 64, 1023, 4095, 4096, 4097 };
    for (size_t byte_i = 0;
         byte_i < sizeof(one_shot_bytes) / sizeof(one_shot_bytes[0]); ++byte_i)
    {
        const size_t bytes = one_shot_bytes[byte_i];
        const bool expect_direct =
            diagnostic_direct ||
            (LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT != 0 &&
             production_direct_shape && bytes >= 4096);
        const unsigned setup_modes[] = { 0, 3 };
        for (size_t mode_i = 0;
             mode_i < sizeof(setup_modes) / sizeof(setup_modes[0]); ++mode_i)
        {
            require(leopard2_internal::SetOneShotPlanSetupModeForDiagnostics(
                    setup_modes[mode_i]),
                "one-shot setup mode select failed");
            leo2_decode_plan* diagnostic_plan = NULL;
            const leo2_result diagnostic_result =
                leopard2_internal::CreateOneShotTransformPlanForDiagnostics(
                    codec, bytes, original_present.data(),
                    recovery_present.data(), &diagnostic_plan);
            require(diagnostic_result ==
                    (expect_direct ? LEO2_UNSUPPORTED : LEO2_SUCCESS),
                "one-shot factory selected the wrong side of the cold threshold");
            require((diagnostic_plan != NULL) == !expect_direct,
                "one-shot diagnostic factory returned an inconsistent handle");
            leo2_decode_plan_destroy(diagnostic_plan);
        }
        require(leopard2_internal::SetOneShotPlanSetupModeForDiagnostics(3),
            "one-shot setup mode restore failed");

        size_t scratch_bytes = 0;
        require_result(leo2_decode_scratch_size(
            codec, bytes, &scratch_bytes), "one-shot scratch query");
        Buffers buffers(k, r, losses, bytes, scratch_bytes, field,
            generator, UINT64_C(0x4f4e4553484f5400) + bytes + k * 257U);
        require_result(leo2_decode(codec, bytes, original_present.data(),
            recovery_present.data(), buffers.original_input.data(),
            buffers.recovery_input.data(), buffers.restored_output.data(),
            buffers.scratch.data, buffers.scratch.bytes),
            "one-shot dual decode");
        buffers.verify(losses);
    }

    size_t below_scratch = 0;
    size_t above_scratch = 0;
    require_result(leo2_decode_plan_scratch_size(
        plan, 1023, &below_scratch), "mixed batch below scratch");
    require_result(leo2_decode_plan_scratch_size(
        plan, 1024, &above_scratch), "mixed batch above scratch");
    Buffers below(k, r, losses, 1023, below_scratch, field, generator,
        UINT64_C(0x4d49584544000001) + k);
    Buffers above(k, r, losses, 1024, above_scratch, field, generator,
        UINT64_C(0x4d49584544000002) + k);
    leo2_decode_batch_item items[2] = { below.item(), above.item() };
    require_result(leo2_decode_plan_execute_batch(plan, items, 2),
        "mixed dual batch execute");
    below.verify(losses);
    above.verify(losses);

    below.clear_restored(losses);
    above.clear_restored(losses);
    leo2_decode_batch_binding* binding = NULL;
    require_result(leo2_decode_batch_binding_create(
        plan, items, 2, &binding), "mixed dual binding create");
    require_result(leo2_decode_batch_binding_execute(binding),
        "mixed dual binding execute");
    below.verify(losses);
    above.verify(losses);
    leo2_decode_batch_binding_destroy(binding);

    below.clear_restored(losses);
    above.clear_restored(losses);
    std::atomic<int> below_result(LEO2_INTERNAL_ERROR);
    std::atomic<int> above_result(LEO2_INTERNAL_ERROR);
    std::thread below_thread([&]() {
        below_result.store(leo2_decode_plan_execute(plan, 1023,
            below.original_input.data(), below.recovery_input.data(),
            below.restored_output.data(), below.scratch.data,
            below.scratch.bytes));
    });
    std::thread above_thread([&]() {
        above_result.store(leo2_decode_plan_execute(plan, 1024,
            above.original_input.data(), above.recovery_input.data(),
            above.restored_output.data(), above.scratch.data,
            above.scratch.bytes));
    });
    below_thread.join();
    above_thread.join();
    require(below_result.load() == LEO2_SUCCESS &&
            above_result.load() == LEO2_SUCCESS,
        "concurrent mixed-route execution failed");
    below.verify(losses);
    above.verify(losses);

    leo2_decode_plan_destroy(transform_plan);
    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(transform_codec);
    leo2_codec_destroy(codec);
}

void verify_irregular_presence_boundary(leo2_context* context)
{
    const unsigned k = 16;
    const unsigned r = 8;
    const unsigned missing_indices[] = { 1, 4, 7, 12, 15 };
    const std::vector<unsigned> missing(
        missing_indices,
        missing_indices + sizeof(missing_indices) / sizeof(missing_indices[0]));
    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 1);
    for (size_t missing_i = 0; missing_i < missing.size(); ++missing_i)
        original_present[missing[missing_i]] = 0;
    recovery_present[0] = 0;
    recovery_present[6] = 0;

    const BinaryField field = leopard2_test::make_legacy_gf8();
    const ProfileLayout layout = leopard2_test::make_profile_layout(
        leopard2_test::kLegacyHigh, k, r);
    const Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);
    leo2_codec* codec = make_codec(context, k, r, 0);
    leo2_codec* transform_codec = make_codec(context, k, r,
        LEO2_CODEC_FORCE_SPECIALIZED_DECODE);
    leo2_decode_plan* plan = make_plan(
        codec, original_present, recovery_present);
    leo2_decode_plan* transform_plan = make_plan(
        transform_codec, original_present, recovery_present);
    const bool diagnostic_direct =
        LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE != 0;
    const bool expect_direct_capability = diagnostic_direct ||
        LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT != 0;
    direct_retained_bytes(plan, expect_direct_capability,
        "irregular direct storage introspection failed");
    direct_retained_bytes(transform_plan, false,
        "irregular forced-transform storage introspection failed");
    const bool expect_locator_terms =
        LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE == 0 &&
        LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT != 0 &&
        LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS != 0;
    require(leopard2_internal::DecodePlanUsesLocatorDirectTerms(plan) ==
            expect_locator_terms,
        "irregular direct-term provenance disagrees with the selector");
    require(!leopard2_internal::DecodePlanUsesLocatorDirectTerms(
            transform_plan),
        "irregular forced transform retained locator-derived terms");

    const size_t byte_counts[] = { 1023, 1024, 1025 };
    for (size_t byte_i = 0;
         byte_i < sizeof(byte_counts) / sizeof(byte_counts[0]); ++byte_i)
    {
        const size_t bytes = byte_counts[byte_i];
        const bool expect_direct = diagnostic_direct ||
            (LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT != 0 && bytes >= 1024);
        leopard2_internal::DecodePathInfo path = {};
        require_result(leopard2_internal::GetDecodePlanPathInfo(
            plan, bytes, false, &path), "irregular path introspection");
        require((path.path == leopard2_internal::kDecodePathDirect) ==
                expect_direct,
            "irregular dual plan selected the wrong boundary route");
        const leopard2_internal::DirectRepairExecutor expected_executor =
            expect_direct
                ? (LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE == 1
                    ? leopard2_internal::kDirectRepairExecutorOutputMajor
                    : leopard2_internal::kDirectRepairExecutorSourceMajor)
                : leopard2_internal::kDirectRepairExecutorNone;
        require(path.direct_executor == expected_executor,
            "irregular dual plan reported the wrong byte executor");

        size_t scratch_bytes = 0;
        size_t transform_scratch_bytes = 0;
        require_result(leo2_decode_plan_scratch_size(
            plan, bytes, &scratch_bytes), "irregular dual scratch query");
        require_result(leo2_decode_plan_scratch_size(transform_plan, bytes,
            &transform_scratch_bytes), "irregular transform scratch query");
        require(expect_direct
                ? scratch_bytes < transform_scratch_bytes
                : scratch_bytes == transform_scratch_bytes,
            "irregular scratch query disagrees with selected route");

        Buffers buffers(k, r, 0, bytes, scratch_bytes, field, generator,
            UINT64_C(0x4952524547554c00) + bytes);
        buffers.apply_presence(original_present, recovery_present);
        require_result(leo2_decode_plan_execute(plan, bytes,
            buffers.original_input.data(), buffers.recovery_input.data(),
            buffers.restored_output.data(), buffers.scratch.data,
            buffers.scratch.bytes), "irregular dual plan execute");
        buffers.verify_missing(missing);

        Buffers transform_buffers(k, r, 0, bytes, transform_scratch_bytes,
            field, generator, UINT64_C(0x5452414e53464f00) + bytes);
        transform_buffers.apply_presence(original_present, recovery_present);
        require_result(leo2_decode_plan_execute(transform_plan, bytes,
            transform_buffers.original_input.data(),
            transform_buffers.recovery_input.data(),
            transform_buffers.restored_output.data(),
            transform_buffers.scratch.data,
            transform_buffers.scratch.bytes),
            "irregular forced-transform plan execute");
        transform_buffers.verify_missing(missing);
    }

    leo2_decode_plan_destroy(transform_plan);
    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(transform_codec);
    leo2_codec_destroy(codec);
}

} // namespace

int main()
{
    try
    {
        leo2_context_options options = {};
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AVX2;
        options.thread_count = 1;
        leo2_context* context = NULL;
        const leo2_result context_result =
            leo2_context_create(&options, &context);
        if (context_result == LEO2_UNSUPPORTED)
        {
            std::cout << "small dual-direct test skipped: AVX2 unavailable\n";
            return 0;
        }
        require_result(context_result, "AVX2 context create");
        require(leo2_context_backend(context) == LEO2_BACKEND_AVX2,
            "explicit AVX2 context selected another backend");
        if ((leo2_context_field_mask(context) & LEO2_FIELD_MASK_GF8) == 0)
        {
            leo2_context_destroy(context);
            std::cout << "small dual-direct test skipped: GF8 disabled\n";
            return 0;
        }

#if LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE == 0 && \
    LEO2_ENABLE_GF8_SMALL_DUAL_DIRECT && \
    LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS
        verify_locator_pattern_sweep(context);
        verify_locator_coefficient_map(context, 7, 5,
            std::vector<uint8_t>{ 0, 1, 0, 0, 1, 0, 0 },
            std::vector<uint8_t>{ 1, 1, 1, 1, 1 });
        verify_locator_coefficient_map(context, 8, 8,
            std::vector<uint8_t>{ 0, 1, 0, 1, 0, 0, 1, 0 },
            std::vector<uint8_t>{ 1, 0, 1, 1, 1, 1, 0, 1 });
        verify_locator_coefficient_map(context, 8, 8,
            std::vector<uint8_t>{ 0, 0, 0, 1, 0, 0, 0, 0 },
            std::vector<uint8_t>{ 1, 1, 0, 1, 1, 1, 1, 1 });
        verify_locator_coefficient_map(context, 5, 8,
            std::vector<uint8_t>{ 0, 0, 0, 0, 0 },
            std::vector<uint8_t>{ 0, 1, 0, 1, 1, 1, 1, 1 });
        {
            std::vector<uint8_t> original_present(16, 1);
            const unsigned missing[] = { 0, 3, 7, 11, 15 };
            for (size_t i = 0; i < sizeof(missing) / sizeof(missing[0]); ++i)
                original_present[missing[i]] = 0;
            verify_locator_coefficient_map(context, 16, 5,
                original_present, std::vector<uint8_t>(5, 1));
        }
        {
            std::vector<uint8_t> original_present(16, 1);
            const unsigned missing[] = { 1, 4, 7, 12, 15 };
            for (size_t i = 0; i < sizeof(missing) / sizeof(missing[0]); ++i)
                original_present[missing[i]] = 0;
            std::vector<uint8_t> recovery_present(8, 1);
            recovery_present[0] = 0;
            recovery_present[6] = 0;
            verify_locator_coefficient_map(context, 16, 8,
                original_present, recovery_present);
        }
        {
            std::vector<uint8_t> original_present(16, 1);
            const unsigned missing[] = { 0, 2, 4, 6, 9, 11, 13, 15 };
            for (size_t i = 0; i < sizeof(missing) / sizeof(missing[0]); ++i)
                original_present[missing[i]] = 0;
            verify_locator_coefficient_map(context, 16, 8,
                original_present, std::vector<uint8_t>(8, 1));
        }
#endif

        verify_routes_and_execution(context, 5, 5, 5);
        verify_routes_and_execution(context, 5, 8, 5);
        verify_routes_and_execution(context, 6, 8, 6);
        verify_routes_and_execution(context, 7, 7, 7);
        verify_routes_and_execution(context, 8, 8, 5);
        verify_routes_and_execution(context, 8, 8, 8);
        verify_routes_and_execution(context, 16, 8, 5);
        verify_routes_and_execution(context, 16, 8, 8);
        const size_t raw_native_high_oracle_cases =
            verify_raw_native_high_oracle_matrix(context);
        require(raw_native_high_oracle_cases == 19536,
            "raw native-high oracle matrix case count drifted");
        verify_irregular_presence_boundary(context);
        verify_transform_only(context, 16, 9, 8,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        verify_transform_only(context, 17, 8, 8,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        verify_transform_only(context, 8, 8, 8,
            LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8);
        if ((leo2_context_field_mask(context) & LEO2_FIELD_MASK_GF16) != 0)
        {
            verify_transform_only(context, 8, 8, 8,
                LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16);
        }

        leo2_context_destroy(context);

        const leo2_backend fallback_backends[] = {
            LEO2_BACKEND_SCALAR, LEO2_BACKEND_SSSE3
        };
        for (size_t backend_i = 0;
             backend_i < sizeof(fallback_backends) /
                 sizeof(fallback_backends[0]); ++backend_i)
        {
            options.backend = fallback_backends[backend_i];
            context = NULL;
            const leo2_result result = leo2_context_create(&options, &context);
            if (result == LEO2_UNSUPPORTED)
                continue;
            require_result(result, "fallback context create");
            verify_transform_only(context, 8, 8, 8,
                LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
            leo2_context_destroy(context);
        }
        std::cout << "small dual transform/direct checks passed"
                  << " direct_retained_max_bytes="
                  << g_max_direct_retained_bytes << "\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "small dual transform/direct failure: "
                  << error.what() << std::endl;
        return 1;
    }
}
