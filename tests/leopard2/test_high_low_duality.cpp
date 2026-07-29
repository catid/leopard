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
    ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
*/

/*
    Paper-duality diagnostic.

    If P = ceil_pow2(K) = ceil_pow2(R) = T, both public profiles have the same
    parent [2P,P] code.  Translation by omega_P is coordinate xor P, swapping
    the two P-point cosets.  Consequently LOW_V1 Algorithm 4 can recover a
    LEGACY_HIGH_V1 codeword after translating its immutable plan and pointer
    views, without changing the public wire profile or parity ordinal bytes.
*/

#include "Leopard2Direct.h"
#include "Leopard2Dispatch.h"
#include "leopard.h"
#include "leopard2.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <new>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

#ifndef LEO2_ENABLE_TEST_HOOKS
#error "high/low duality test requires LEO2_ENABLE_TEST_HOOKS"
#endif

#ifndef LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE
#define LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE 0
#endif

namespace {

typedef std::vector<uint8_t> Bytes;
typedef std::vector<Bytes> Shards;

struct Counts
{
    uint64_t parity_cases;
    uint64_t decode_cases;
    uint64_t restored_shards;
    uint64_t exhaustive_patterns;
    uint64_t random_patterns;
    uint64_t backends;
    uint64_t parity_rebuilds;

    Counts()
        : parity_cases(0)
        , decode_cases(0)
        , restored_shards(0)
        , exhaustive_patterns(0)
        , random_patterns(0)
        , backends(0)
        , parity_rebuilds(0)
    {}
};

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void require_result(
    leo2_result actual,
    leo2_result expected,
    const std::string& operation)
{
    if (actual == expected)
        return;
    std::ostringstream stream;
    stream << operation << ": expected " << leo2_result_string(expected)
           << " (" << static_cast<int>(expected) << "), received "
           << leo2_result_string(actual) << " (" << static_cast<int>(actual)
           << ')';
    throw std::runtime_error(stream.str());
}

void require_success(leo2_result result, const std::string& operation)
{
    require_result(result, LEO2_SUCCESS, operation);
}

class AlignedBuffer
{
public:
    explicit AlignedBuffer(size_t bytes)
        : data_(NULL)
    {
        if (bytes == 0)
            return;
#if defined(_MSC_VER)
        data_ = _aligned_malloc(bytes, leo2_scratch_alignment());
        if (!data_)
            throw std::bad_alloc();
#else
        if (posix_memalign(&data_, leo2_scratch_alignment(), bytes) != 0)
            throw std::bad_alloc();
#endif
        memset(data_, 0, bytes);
    }

    ~AlignedBuffer()
    {
#if defined(_MSC_VER)
        _aligned_free(data_);
#else
        free(data_);
#endif
    }

    void* data() { return data_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);

    void* data_;
};

uint64_t splitmix64(uint64_t* state)
{
    uint64_t value = (*state += UINT64_C(0x9e3779b97f4a7c15));
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

Shards make_originals(unsigned count, size_t bytes, uint64_t seed)
{
    Shards shards(count, Bytes(bytes, 0));
    for (unsigned shard = 0; shard < count; ++shard)
    {
        for (size_t offset = 0; offset < bytes; ++offset)
        {
            shards[shard][offset] = static_cast<uint8_t>(
                splitmix64(&seed) + shard * 37u + offset * 149u);
        }
    }
    return shards;
}

std::vector<const void*> const_pointers(const Shards& shards)
{
    std::vector<const void*> result(shards.size(), NULL);
    for (size_t i = 0; i < shards.size(); ++i)
        result[i] = shards[i].empty() ? NULL : &shards[i][0];
    return result;
}

std::vector<void*> mutable_pointers(Shards* shards)
{
    std::vector<void*> result(shards->size(), NULL);
    for (size_t i = 0; i < shards->size(); ++i)
        result[i] = (*shards)[i].empty() ? NULL : &(*shards)[i][0];
    return result;
}

leo2_codec* create_codec(
    leo2_context* context,
    unsigned k,
    unsigned r,
    leo2_profile profile,
    leo2_field field,
    uint32_t flags)
{
    leo2_codec_options options;
    memset(&options, 0, sizeof(options));
    options.struct_size = sizeof(options);
    options.flags = flags;
    leo2_codec* codec = NULL;
    require_success(leo2_codec_create(
        context, k, r, profile, field, &options, &codec), "codec create");
    require(codec != NULL, "codec create returned null");
    require(leo2_codec_profile(codec) == profile,
        "explicit profile changed during codec creation");
    require(leo2_codec_field(codec) == field,
        "explicit field changed during codec creation");
    return codec;
}

Shards encode(
    const leo2_codec* codec,
    const Shards& originals,
    unsigned recovery_count,
    size_t bytes)
{
    Shards recovery(recovery_count, Bytes(bytes, 0));
    std::vector<const void*> input = const_pointers(originals);
    std::vector<void*> output = mutable_pointers(&recovery);
    size_t scratch_bytes = 0;
    require_success(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        "encode scratch query");
    AlignedBuffer scratch(scratch_bytes);
    require_success(leo2_encode(codec, bytes, &input[0], &output[0],
        scratch.data(), scratch_bytes), "encode");
    return recovery;
}

Shards decode(
    leo2_codec* codec,
    const Shards& originals,
    const Shards& recovery,
    const std::vector<uint8_t>& original_present,
    const std::vector<uint8_t>& recovery_present,
    size_t bytes,
    leo2_test_decode_mode decode_mode,
    bool test_plan_immutability,
    bool expect_tiled)
{
    const bool translated =
        decode_mode == LEO2_TEST_DECODE_FORCE_TRANSLATED_LOW;
    require_success(leo2_test_codec_set_decode_mode(codec,
        decode_mode),
        "set decode mode");
    leo2_decode_plan* plan = NULL;
    require_success(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), "decode plan create");
    require(plan != NULL, "decode plan create returned null");
    const bool has_missing_original = std::find(
        original_present.begin(), original_present.end(), 0) !=
            original_present.end();
    if (decode_mode != LEO2_TEST_DECODE_AUTO)
    {
        require(leo2_test_decode_plan_uses_translated_low(plan) ==
                (translated && has_missing_original ? 1 : 0),
            "decode plan did not capture the requested translation mode");
    }
    if (test_plan_immutability)
    {
        require_success(leo2_test_codec_set_decode_mode(
            codec, LEO2_TEST_DECODE_AUTO), "reset codec decode mode");
        require(leo2_test_decode_plan_uses_translated_low(plan) ==
                (translated && has_missing_original ? 1 : 0),
            "codec mutation changed immutable plan translation metadata");
    }

    std::vector<const void*> original_input = const_pointers(originals);
    std::vector<const void*> recovery_input = const_pointers(recovery);
    Shards restored(originals.size(), Bytes(bytes, 0xa5));
    std::vector<void*> output(originals.size(), NULL);
    for (size_t i = 0; i < originals.size(); ++i)
    {
        if (!original_present[i])
        {
            original_input[i] = NULL;
            output[i] = &restored[i][0];
        }
    }
    for (size_t i = 0; i < recovery.size(); ++i)
        if (!recovery_present[i])
            recovery_input[i] = NULL;

    size_t scratch_bytes = 0;
    require_success(leo2_decode_plan_scratch_size(plan, bytes, &scratch_bytes),
        "decode scratch query");
    unsigned missing_count = 0;
    for (size_t i = 0; i < original_present.size(); ++i)
        missing_count += original_present[i] ? 0u : 1u;
    if (translated && missing_count != 0)
    {
        leopard2_internal::DecodePathInfo path;
        require_success(leopard2_internal::GetDecodePlanPathInfo(
            plan, bytes, false, &path), "translated path introspection");
        require(path.path == (expect_tiled
                ? leopard2_internal::kDecodePathTiled
                : leopard2_internal::kDecodePathMaterialized),
            "translated plan selected the wrong forced workspace traversal");
        require(path.rule ==
                leopard2_internal::kDecodeRuleTranslatedLow,
            "translated plan did not report its production dispatch rule");
        require(path.required_work_slots == leo2_codec_parent_count(codec),
            "translated P=T plan did not report its exact N-slot workspace");
    }
    AlignedBuffer scratch(scratch_bytes);
    if (translated)
    {
        leo2_test_reset_low_reveal_counts();
        leo2_test_reset_high_reveal_counts();
    }
    const unsigned execution_count = translated ? 2u : 1u;
    for (unsigned execution = 0; execution < execution_count; ++execution)
    {
        for (size_t i = 0; i < restored.size(); ++i)
            if (!original_present[i])
                std::fill(restored[i].begin(), restored[i].end(),
                    static_cast<uint8_t>(0xa5 + execution));
        require_success(leo2_decode_plan_execute(plan, bytes,
            &original_input[0], &recovery_input[0], &output[0],
            scratch.data(), scratch_bytes), "decode execute");
    }
    if (translated)
    {
        // A translated Algorithm 4 plan retains a legacy-high profile ID, but
        // must never enter Algorithm 5's direct-output suppression.  This is
        // the regression gate for the materialized-high reveal integration.
        require(leo2_test_high_direct_reveal_shards() == 0 &&
                leo2_test_high_materialized_direct_reveal_shards() == 0 &&
                leo2_test_high_scratch_reveal_shards() == 0,
            "translated Algorithm 4 plan entered Algorithm 5 reveal handling");
        const uint64_t expected_direct = bytes >= 64
            ? static_cast<uint64_t>(missing_count) * execution_count : 0;
        const uint64_t expected_scratch = (bytes & 63u) != 0
            ? static_cast<uint64_t>(missing_count) * execution_count : 0;
        require(leo2_test_low_direct_reveal_shards() == expected_direct &&
                leo2_test_low_scratch_reveal_shards() == expected_scratch,
            "translated Algorithm 4 reveal/gather accounting is incomplete");
    }
    leo2_decode_plan_destroy(plan);
    return restored;
}

void require_selected_path(
    leo2_codec* codec,
    const std::vector<uint8_t>& original_present,
    const std::vector<uint8_t>& recovery_present,
    size_t bytes,
    leopard2_internal::DecodePath expected_path,
    leopard2_internal::DecodePathRule expected_rule)
{
    leo2_decode_plan* plan = NULL;
    require_success(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), "dispatch plan create");
    leopard2_internal::DecodePathInfo path;
    require_success(leopard2_internal::GetDecodePlanPathInfo(
        plan, bytes, false, &path), "dispatch path introspection");
    require(path.path == expected_path && path.rule == expected_rule,
        "decode dispatcher selected an unexpected path/rule");
    leo2_decode_plan_destroy(plan);
}

void require_recovery(
    const Shards& originals,
    const Shards& translated,
    const Shards& translated_tiled,
    const Shards& ordinary_high,
    const Shards& ordinary_low,
    const std::vector<uint8_t>& original_present,
    Counts* counts)
{
    for (size_t i = 0; i < originals.size(); ++i)
    {
        if (original_present[i])
            continue;
        require(translated[i] == originals[i],
            "translated Algorithm 4 restored the wrong original");
        require(translated_tiled[i] == originals[i],
            "translated tiled Algorithm 4 restored the wrong original");
        require(ordinary_high[i] == originals[i],
            "ordinary Algorithm 5 restored the wrong original");
        require(ordinary_low[i] == originals[i],
            "ordinary low-profile decoder restored the wrong original");
        require(translated[i] == translated_tiled[i] &&
                translated[i] == ordinary_high[i] &&
                translated[i] == ordinary_low[i],
            "high/low recovery results differ");
        ++counts->restored_shards;
    }
}

void run_pattern(
    leo2_codec* translated_high,
    leo2_codec* translated_high_tiled,
    leo2_codec* ordinary_high,
    leo2_codec* ordinary_low,
    const Shards& originals,
    const Shards& recovery,
    const std::vector<uint8_t>& original_present,
    const std::vector<uint8_t>& recovery_present,
    size_t bytes,
    bool test_plan_immutability,
    Counts* counts)
{
    const Shards translated = decode(translated_high, originals, recovery,
        original_present, recovery_present, bytes,
        LEO2_TEST_DECODE_FORCE_TRANSLATED_LOW,
        test_plan_immutability, false);
    const Shards translated_tiled = decode(translated_high_tiled, originals,
        recovery, original_present, recovery_present, bytes,
        LEO2_TEST_DECODE_FORCE_TRANSLATED_LOW, false, true);
    const Shards high = decode(ordinary_high, originals, recovery,
        original_present, recovery_present, bytes,
        LEO2_TEST_DECODE_FORCE_NATIVE_HIGH, false, false);
    const Shards low = decode(ordinary_low, originals, recovery,
        original_present, recovery_present, bytes,
        LEO2_TEST_DECODE_AUTO, false, false);
    require_recovery(originals, translated, translated_tiled, high, low,
        original_present, counts);
    Shards complete = originals;
    for (size_t i = 0; i < complete.size(); ++i)
        if (!original_present[i])
            complete[i] = translated[i];
    require(encode(translated_high, complete,
            static_cast<unsigned>(recovery.size()), bytes) == recovery,
        "parity rebuild after translated recovery changed legacy-high parity");
    ++counts->parity_rebuilds;
    ++counts->decode_cases;
}

unsigned popcount(unsigned value)
{
    unsigned count = 0;
    while (value != 0)
    {
        count += value & 1u;
        value >>= 1;
    }
    return count;
}

void run_case(
    leo2_context* context,
    unsigned k,
    unsigned r,
    leo2_field field,
    size_t bytes,
    bool exhaustive,
    unsigned random_pattern_count,
    uint64_t seed,
    Counts* counts)
{
    leo2_codec* translated_high = create_codec(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, field,
        LEO2_CODEC_FORCE_MATERIALIZED_DECODE);
    leo2_codec* translated_high_tiled = create_codec(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, field,
        LEO2_CODEC_FORCE_TILED_DECODE);
    leo2_codec* ordinary_high = create_codec(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, field,
        LEO2_CODEC_FORCE_MATERIALIZED_DECODE);
    leo2_codec* ordinary_low = create_codec(context, k, r,
        LEO2_PROFILE_LOW_V1, field,
        LEO2_CODEC_FORCE_MATERIALIZED_DECODE);

    require(leo2_test_codec_translated_low_capable(translated_high) == 1,
        "P=T high codec was not translation-capable");
    require_success(leo2_test_codec_set_decode_mode(translated_high,
        LEO2_TEST_DECODE_FORCE_TRANSLATED_LOW),
        "enable translated low decoder");
    require_success(leo2_test_codec_set_decode_mode(translated_high_tiled,
        LEO2_TEST_DECODE_FORCE_TRANSLATED_LOW),
        "enable translated tiled low decoder");
    require_result(leo2_test_codec_set_decode_mode(ordinary_low,
        LEO2_TEST_DECODE_FORCE_TRANSLATED_LOW), LEO2_UNSUPPORTED,
        "reject translation on low profile");

    const Shards originals = make_originals(k, bytes, seed);
    const Shards high_recovery = encode(translated_high, originals, r, bytes);
    const Shards low_recovery = encode(ordinary_low, originals, r, bytes);
    require(high_recovery == low_recovery,
        "translated high and low profiles emitted different parity ordinals");
    ++counts->parity_cases;

    bool tested_immutability = false;
    // A true no-loss call stays a no-op even when every parity is supplied.
    std::vector<uint8_t> all_original_present(k, 1);
    std::vector<uint8_t> all_recovery_present(r, 1);
    run_pattern(translated_high, translated_high_tiled, ordinary_high,
        ordinary_low, originals, high_recovery, all_original_present,
        all_recovery_present, bytes, false, counts);

    // Keep surplus parity present so deterministic virtual erasures, fixed
    // shortening, and puncturing all participate in the translated plan.
    std::vector<uint8_t> surplus_original_present(k, 1);
    const unsigned surplus_losses = std::min<unsigned>(5, std::min(k, r));
    for (unsigned i = 0; i < surplus_losses; ++i)
        surplus_original_present[(i * 7u + 1u) % k] = 0;
    run_pattern(translated_high, translated_high_tiled, ordinary_high,
        ordinary_low, originals, high_recovery, surplus_original_present,
        all_recovery_present, bytes, true, counts);
    tested_immutability = true;

    // This is the paper's balanced full-loss corner and the production
    // selector's existing generic-special-case boundary.
    if (k == 128 && r == 128)
    {
        std::vector<uint8_t> none_original_present(k, 0);
        run_pattern(translated_high, translated_high_tiled, ordinary_high,
            ordinary_low, originals, high_recovery, none_original_present,
            all_recovery_present, bytes, false, counts);

        leo2_codec* automatic_high = create_codec(context, k, r,
            LEO2_PROFILE_LEGACY_HIGH_V1, field, 0);
        require_selected_path(automatic_high, none_original_present,
            all_recovery_present, bytes,
            leopard2_internal::kDecodePathMaterialized,
            leopard2_internal::kDecodeRuleTranslatedLow);
        const Shards automatic = decode(automatic_high, originals, high_recovery,
            none_original_present, all_recovery_present, bytes,
            LEO2_TEST_DECODE_AUTO, false, false);
        for (unsigned i = 0; i < k; ++i)
            require(automatic[i] == originals[i],
                "production translated Algorithm 4 disagrees with source");
        leo2_codec_destroy(automatic_high);
    }

    if (exhaustive)
    {
        const unsigned public_count = k + r;
        require(public_count < 31, "exhaustive mask is too wide");
        const unsigned mask_limit = 1u << public_count;
        for (unsigned mask = 0; mask < mask_limit; ++mask)
        {
            if (popcount(mask) != k)
                continue;
            std::vector<uint8_t> original_present(k, 0);
            std::vector<uint8_t> recovery_present(r, 0);
            for (unsigned i = 0; i < k; ++i)
                original_present[i] = (mask >> i) & 1u;
            for (unsigned i = 0; i < r; ++i)
                recovery_present[i] = (mask >> (k + i)) & 1u;
            run_pattern(translated_high, translated_high_tiled,
                ordinary_high, ordinary_low,
                originals, high_recovery, original_present, recovery_present,
                bytes, !tested_immutability, counts);
            tested_immutability = true;
            ++counts->exhaustive_patterns;
        }
    }
    else
    {
        std::vector<unsigned> coordinates(k + r, 0);
        for (unsigned i = 0; i < k + r; ++i)
            coordinates[i] = i;
        for (unsigned pattern = 0; pattern < random_pattern_count; ++pattern)
        {
            for (unsigned i = k + r; i > 1; --i)
            {
                const unsigned j = static_cast<unsigned>(splitmix64(&seed) % i);
                std::swap(coordinates[i - 1], coordinates[j]);
            }
            std::vector<uint8_t> original_present(k, 0);
            std::vector<uint8_t> recovery_present(r, 0);
            for (unsigned i = 0; i < k; ++i)
            {
                const unsigned coordinate = coordinates[i];
                if (coordinate < k)
                    original_present[coordinate] = 1;
                else
                    recovery_present[coordinate - k] = 1;
            }
            run_pattern(translated_high, translated_high_tiled,
                ordinary_high, ordinary_low,
                originals, high_recovery, original_present, recovery_present,
                bytes, !tested_immutability, counts);
            tested_immutability = true;
            ++counts->random_patterns;
        }
    }

    leo2_codec_destroy(ordinary_low);
    leo2_codec_destroy(ordinary_high);
    leo2_codec_destroy(translated_high_tiled);
    leo2_codec_destroy(translated_high);
}

void test_rejection(leo2_context* context)
{
    leo2_codec* unequal = create_codec(context, 9, 8,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 0);
    require(leo2_test_codec_translated_low_capable(unequal) == 0,
        "unequal rounded sides were reported translation-capable");
    require_result(leo2_test_codec_set_decode_mode(unequal,
        LEO2_TEST_DECODE_FORCE_TRANSLATED_LOW), LEO2_UNSUPPORTED,
        "reject unequal rounded sides");
    require_result(leo2_test_codec_set_decode_mode(unequal,
        LEO2_TEST_DECODE_FORCE_NATIVE_HIGH), LEO2_UNSUPPORTED,
        "reject native-high control outside the duality region");
    require_result(leo2_test_codec_set_decode_mode(unequal,
        -1), LEO2_INVALID_ARGUMENT,
        "reject negative decode mode");
    require_result(leo2_test_codec_set_decode_mode(unequal,
        99), LEO2_INVALID_ARGUMENT,
        "reject unknown decode mode");
    leo2_codec_destroy(unequal);

    leo2_codec* forced_generic = create_codec(context, 8, 8,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        LEO2_CODEC_FORCE_GENERIC_DECODE);
    require_result(leo2_test_codec_set_decode_mode(forced_generic,
        LEO2_TEST_DECODE_FORCE_TRANSLATED_LOW), LEO2_INVALID_ARGUMENT,
        "reject translation with forced generic decoder");
    require_result(leo2_test_codec_set_decode_mode(forced_generic,
        LEO2_TEST_DECODE_FORCE_NATIVE_HIGH), LEO2_INVALID_ARGUMENT,
        "reject native-high control with forced generic decoder");
    leo2_codec_destroy(forced_generic);
}

void test_direct_dispatch_bypass(leo2_context* context)
{
    const unsigned k = 8;
    const unsigned r = 8;
    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 1);
    original_present[3] = 0;

    leo2_codec* automatic = create_codec(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 0);
    require_selected_path(automatic, original_present, recovery_present, 64,
        leopard2_internal::kDecodePathDirect,
        leopard2_internal::kDecodeRuleDirect);

    // The native Algorithm 5 control exists only in hook builds and must
    // bypass production direct repair so attribution cannot silently time the
    // small-loss solver instead of the requested transform.
    require_success(leo2_test_codec_set_decode_mode(automatic,
        LEO2_TEST_DECODE_FORCE_NATIVE_HIGH),
        "force native Algorithm 5 direct bypass");
    require_selected_path(automatic, original_present, recovery_present, 64,
        leopard2_internal::kDecodePathMaterialized,
        leopard2_internal::kDecodeRuleWorkspaceMaterialized);
    require_success(leo2_test_codec_set_decode_mode(automatic,
        LEO2_TEST_DECODE_AUTO), "restore automatic direct repair");
    require_selected_path(automatic, original_present, recovery_present, 64,
        leopard2_internal::kDecodePathDirect,
        leopard2_internal::kDecodeRuleDirect);

    // Once the measured direct-repair region is exceeded, production AUTO
    // deterministically selects the translated Algorithm 4 plan.
    for (unsigned i = 0; i < 5; ++i)
        original_present[i] = 0;
    leo2_decode_plan* automatic_transform_plan = NULL;
    require_success(leo2_decode_plan_create(automatic, &original_present[0],
        &recovery_present[0], &automatic_transform_plan),
        "automatic translated plan create");
    leopard2_internal::DecodePathInfo automatic_path;
    require_success(leopard2_internal::GetDecodePlanPathInfo(
        automatic_transform_plan, 64, false, &automatic_path),
        "automatic translated path introspection");
#if LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE != 0
    if (leo2_context_backend(context) == LEO2_BACKEND_AVX2)
    {
        require(leo2_test_decode_plan_uses_translated_low(
                    automatic_transform_plan) == 0,
            "experimental AUTO direct plan retained translated Algorithm 4");
        require(automatic_path.path == leopard2_internal::kDecodePathDirect &&
                automatic_path.rule == leopard2_internal::kDecodeRuleDirect,
            "experimental AUTO reported the wrong direct dispatch rule");
        const leopard2_internal::DirectRepairExecutor expected_executor =
            LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE == 2
                ? leopard2_internal::kDirectRepairExecutorSourceMajor
                : leopard2_internal::kDirectRepairExecutorOutputMajor;
        require(automatic_path.direct_executor == expected_executor,
            "experimental AUTO reported the wrong direct executor");
    }
    else
#endif
    {
        require(leo2_test_decode_plan_uses_translated_low(
                    automatic_transform_plan) == 1,
            "production AUTO did not select translated Algorithm 4");
        require(automatic_path.path ==
                    leopard2_internal::kDecodePathMaterialized &&
                automatic_path.rule ==
                    leopard2_internal::kDecodeRuleTranslatedLow,
            "production AUTO reported the wrong translated dispatch rule");
    }
    leo2_decode_plan_destroy(automatic_transform_plan);
    leo2_codec_destroy(automatic);

    leo2_codec* translated = create_codec(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        LEO2_CODEC_FORCE_MATERIALIZED_DECODE);
    require_success(leo2_test_codec_set_decode_mode(translated,
        LEO2_TEST_DECODE_FORCE_TRANSLATED_LOW),
        "force translated direct-bypass plan");
    require_selected_path(translated, original_present, recovery_present, 64,
        leopard2_internal::kDecodePathMaterialized,
        leopard2_internal::kDecodeRuleTranslatedLow);

    require_success(leo2_test_codec_set_decode_mode(translated,
        LEO2_TEST_DECODE_FORCE_NATIVE_HIGH),
        "force native Algorithm 5 control");
    require_selected_path(translated, original_present, recovery_present, 64,
        leopard2_internal::kDecodePathMaterialized,
        leopard2_internal::kDecodeRuleForcedMaterialized);
    leo2_codec_destroy(translated);
}

void test_legacy_wire_and_recovery(leo2_context* context)
{
    const unsigned k = 8;
    const unsigned r = 8;
    const size_t bytes = 64;
    const Shards originals = make_originals(
        k, bytes, UINT64_C(0x4c454741435932));

    leo2_codec* codec = create_codec(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        LEO2_CODEC_FORCE_MATERIALIZED_DECODE);
    require(leo2_codec_profile(codec) == LEO2_PROFILE_LEGACY_HIGH_V1,
        "translated-capable codec changed its persistent profile identity");
    const Shards parity = encode(codec, originals, r, bytes);

    const unsigned encode_work_count = leo_encode_work_count(k, r);
    require(encode_work_count >= r,
        "old Leopard rejected a translated-capable legacy profile");
    Shards old_encode_work(encode_work_count, Bytes(bytes, 0));
    std::vector<const void*> original_input = const_pointers(originals);
    std::vector<void*> old_encode_output = mutable_pointers(&old_encode_work);
    require(leo_encode(bytes, k, r, encode_work_count, &original_input[0],
            &old_encode_output[0]) == Leopard_Success,
        "old Leopard encode failed");
    for (unsigned i = 0; i < r; ++i)
        require(old_encode_work[i] == parity[i],
            "legacy-high parity differs from old Leopard");

    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 1);
    original_present[1] = original_present[3] = original_present[6] = 0;
    require_success(leo2_test_codec_set_decode_mode(codec,
        LEO2_TEST_DECODE_FORCE_TRANSLATED_LOW),
        "enable translated old-data decode");
    const Shards translated = decode(codec, originals, parity,
        original_present, recovery_present, bytes,
        LEO2_TEST_DECODE_FORCE_TRANSLATED_LOW, false, false);

    const unsigned decode_work_count = leo_decode_work_count(k, r);
    require(decode_work_count >= k, "old Leopard decode work count failed");
    Shards old_decode_work(decode_work_count, Bytes(bytes, 0));
    std::vector<void*> old_decode_output = mutable_pointers(&old_decode_work);
    std::vector<const void*> old_original_input = const_pointers(originals);
    std::vector<const void*> old_recovery_input = const_pointers(parity);
    for (unsigned i = 0; i < k; ++i)
        if (!original_present[i])
            old_original_input[i] = NULL;
    require(leo_decode(bytes, k, r, decode_work_count,
            &old_original_input[0], &old_recovery_input[0],
            &old_decode_output[0]) == Leopard_Success,
        "old Leopard recovery failed");
    for (unsigned i = 0; i < k; ++i)
    {
        if (!original_present[i])
        {
            require(old_decode_work[i] == originals[i] &&
                    translated[i] == old_decode_work[i],
                "translated recovery differs from old Leopard");
        }
    }
    leo2_codec_destroy(codec);
}

void test_unaligned_tail(
    leo2_context* context,
    leo2_field field,
    unsigned k,
    unsigned r,
    size_t bytes)
{
    leo2_codec* codec = create_codec(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, field,
        LEO2_CODEC_FORCE_MATERIALIZED_DECODE);
    require_success(leo2_test_codec_set_decode_mode(codec,
        LEO2_TEST_DECODE_FORCE_TRANSLATED_LOW),
        "enable unaligned translated decoder");

    const Shards originals = make_originals(
        k, bytes, UINT64_C(0x554e414c49474e45) + field);
    const Shards expected_parity = encode(codec, originals, r, bytes);
    Shards original_storage(k, Bytes(bytes + 2, 0xd1));
    Shards parity_storage(r, Bytes(bytes + 2, 0xe2));
    std::vector<const void*> original_input(k, NULL);
    std::vector<void*> parity_output(r, NULL);
    for (unsigned i = 0; i < k; ++i)
    {
        memcpy(&original_storage[i][1], &originals[i][0], bytes);
        original_input[i] = &original_storage[i][1];
    }
    for (unsigned i = 0; i < r; ++i)
        parity_output[i] = &parity_storage[i][1];
    size_t encode_scratch_bytes = 0;
    require_success(leo2_encode_scratch_size(
        codec, bytes, &encode_scratch_bytes), "unaligned encode scratch query");
    AlignedBuffer encode_scratch(encode_scratch_bytes);
    require_success(leo2_encode(codec, bytes, &original_input[0],
        &parity_output[0], encode_scratch.data(), encode_scratch_bytes),
        "unaligned encode");
    for (unsigned i = 0; i < r; ++i)
    {
        require(memcmp(&parity_storage[i][1], &expected_parity[i][0],
                bytes) == 0,
            "unaligned parity differs from aligned legacy-high parity");
        require(parity_storage[i][0] == 0xe2 &&
                parity_storage[i][bytes + 1] == 0xe2,
            "unaligned encode changed a parity guard");
    }

    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 1);
    const unsigned losses = std::min<unsigned>(5, std::min(k, r));
    for (unsigned i = 0; i < losses; ++i)
        original_present[(i * 3u + 1u) % k] = 0;
    leo2_decode_plan* plan = NULL;
    require_success(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), "unaligned decode plan create");
    std::vector<const void*> received_original(k, NULL);
    std::vector<const void*> received_recovery(r, NULL);
    for (unsigned i = 0; i < k; ++i)
        if (original_present[i])
            received_original[i] = &original_storage[i][1];
    for (unsigned i = 0; i < r; ++i)
        received_recovery[i] = &parity_storage[i][1];
    Shards restored(k, Bytes(bytes + 2, 0xa5));
    std::vector<void*> restored_output(k, NULL);
    for (unsigned i = 0; i < k; ++i)
        if (!original_present[i])
            restored_output[i] = &restored[i][1];
    size_t decode_scratch_bytes = 0;
    require_success(leo2_decode_plan_scratch_size(
        plan, bytes, &decode_scratch_bytes), "unaligned decode scratch query");
    AlignedBuffer decode_scratch(decode_scratch_bytes);
    require_success(leo2_decode_plan_execute(plan, bytes,
        &received_original[0], &received_recovery[0], &restored_output[0],
        decode_scratch.data(), decode_scratch_bytes), "unaligned decode");
    for (unsigned i = 0; i < k; ++i)
    {
        if (!original_present[i])
        {
            require(memcmp(&restored[i][1], &originals[i][0], bytes) == 0,
                "unaligned translated decode restored the wrong original");
            require(restored[i][0] == 0xa5 &&
                    restored[i][bytes + 1] == 0xa5,
                "unaligned translated decode changed an output guard");
        }
    }
    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
}

void test_batch_and_concurrent_plan(leo2_context* context)
{
    const unsigned k = 16;
    const unsigned r = 9;
    const size_t bytes = 65;
    leo2_codec* codec = create_codec(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        LEO2_CODEC_FORCE_MATERIALIZED_DECODE);
    require_success(leo2_test_codec_set_decode_mode(codec,
        LEO2_TEST_DECODE_FORCE_TRANSLATED_LOW),
        "enable shared translated plan");
    const Shards originals = make_originals(
        k, bytes, UINT64_C(0x4241544348434f4e));
    const Shards parity = encode(codec, originals, r, bytes);
    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 1);
    for (unsigned i = 0; i < 5; ++i)
        original_present[i * 3u] = 0;
    leo2_decode_plan* plan = NULL;
    require_success(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), "shared decode plan create");
    std::vector<const void*> original_input = const_pointers(originals);
    std::vector<const void*> recovery_input = const_pointers(parity);
    for (unsigned i = 0; i < k; ++i)
        if (!original_present[i])
            original_input[i] = NULL;

    size_t scratch_bytes = 0;
    require_success(leo2_decode_plan_scratch_size(plan, bytes, &scratch_bytes),
        "shared plan scratch query");
    Shards restored_a(k, Bytes(bytes, 0xa5));
    Shards restored_b(k, Bytes(bytes, 0x5a));
    std::vector<void*> output_a(k, NULL);
    std::vector<void*> output_b(k, NULL);
    for (unsigned i = 0; i < k; ++i)
    {
        if (!original_present[i])
        {
            output_a[i] = &restored_a[i][0];
            output_b[i] = &restored_b[i][0];
        }
    }
    AlignedBuffer scratch_a(scratch_bytes);
    AlignedBuffer scratch_b(scratch_bytes);
    leo2_decode_batch_item items[2];
    memset(items, 0, sizeof(items));
    items[0].shard_bytes = items[1].shard_bytes = bytes;
    items[0].original = items[1].original = &original_input[0];
    items[0].recovery = items[1].recovery = &recovery_input[0];
    items[0].restored_original = &output_a[0];
    items[1].restored_original = &output_b[0];
    items[0].scratch = scratch_a.data();
    items[1].scratch = scratch_b.data();
    items[0].scratch_bytes = items[1].scratch_bytes = scratch_bytes;
    require_success(leo2_decode_plan_execute_batch(plan, items, 2),
        "translated decode batch");
    for (unsigned i = 0; i < k; ++i)
        if (!original_present[i])
            require(restored_a[i] == originals[i] &&
                    restored_b[i] == originals[i],
                "translated batch decode differs from source");

    for (unsigned i = 0; i < k; ++i)
    {
        std::fill(restored_a[i].begin(), restored_a[i].end(), 0x33);
        std::fill(restored_b[i].begin(), restored_b[i].end(), 0x44);
    }
    leo2_result result_a = LEO2_INTERNAL_ERROR;
    leo2_result result_b = LEO2_INTERNAL_ERROR;
    std::thread first([&]() {
        result_a = leo2_decode_plan_execute(plan, bytes, &original_input[0],
            &recovery_input[0], &output_a[0], scratch_a.data(), scratch_bytes);
    });
    std::thread second([&]() {
        result_b = leo2_decode_plan_execute(plan, bytes, &original_input[0],
            &recovery_input[0], &output_b[0], scratch_b.data(), scratch_bytes);
    });
    first.join();
    second.join();
    require_success(result_a, "first concurrent translated decode");
    require_success(result_b, "second concurrent translated decode");
    for (unsigned i = 0; i < k; ++i)
        if (!original_present[i])
            require(restored_a[i] == originals[i] &&
                    restored_b[i] == originals[i],
                "concurrent immutable-plan decode differs from source");
    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
}

void run_backend_suite(leo2_context* context, Counts* counts)
{
    test_rejection(context);
    test_direct_dispatch_bypass(context);
    test_legacy_wire_and_recovery(context);
    test_unaligned_tail(context, LEO2_FIELD_GF8, 17, 17, 65);
    test_unaligned_tail(context, LEO2_FIELD_GF16, 17, 17, 66);
    test_batch_and_concurrent_plan(context);
    run_case(context, 2, 2, LEO2_FIELD_GF8, 1, true, 0,
        UINT64_C(0x02020001), counts);
    run_case(context, 3, 3, LEO2_FIELD_GF8, 17, true, 0,
        UINT64_C(0x03030011), counts);
    run_case(context, 3, 4, LEO2_FIELD_GF8, 65, true, 0,
        UINT64_C(0x03040041), counts);
    run_case(context, 4, 3, LEO2_FIELD_GF8, 193, true, 0,
        UINT64_C(0x040300c1), counts);

    const unsigned gf8_cases[][2] = {
            {5, 8}, {8, 5}, {9, 16}, {16, 9}, {17, 32}, {31, 17},
            {33, 64}, {63, 33}, {65, 128}, {127, 65}, {128, 127},
            {128, 128}
    };
    const size_t gf8_bytes[] = { 1, 7, 31, 63, 64, 65, 193 };
    for (size_t i = 0; i < sizeof(gf8_cases) / sizeof(gf8_cases[0]); ++i)
    {
        const size_t bytes = gf8_cases[i][0] == 128 &&
                gf8_cases[i][1] == 128
            ? 256
            : gf8_bytes[i %
                (sizeof(gf8_bytes) / sizeof(gf8_bytes[0]))];
        run_case(context, gf8_cases[i][0], gf8_cases[i][1],
            LEO2_FIELD_GF8, bytes,
            false, 8, UINT64_C(0x8f000000) + i, counts);
    }

    const unsigned gf16_cases[][2] = {
        {129, 256}, {255, 129}, {257, 512}, {511, 257}
    };
    const size_t gf16_bytes[] = { 2, 62, 66, 130 };
    for (size_t i = 0;
         i < sizeof(gf16_cases) / sizeof(gf16_cases[0]); ++i)
    {
        run_case(context, gf16_cases[i][0], gf16_cases[i][1],
            LEO2_FIELD_GF16, gf16_bytes[i], false, 6,
            UINT64_C(0x16000000) + i, counts);
    }
}

} // namespace

int main()
{
    try
    {
        Counts counts;
        const leo2_backend backends[] = {
            LEO2_BACKEND_SCALAR,
            LEO2_BACKEND_SSSE3,
            LEO2_BACKEND_AVX2,
            LEO2_BACKEND_AVX512,
            LEO2_BACKEND_NEON,
            LEO2_BACKEND_GFNI
        };
        for (size_t i = 0; i < sizeof(backends) / sizeof(backends[0]); ++i)
        {
            leo2_context_options options;
            memset(&options, 0, sizeof(options));
            options.struct_size = sizeof(options);
            options.backend = backends[i];
            options.thread_count = 1;
            leo2_context* context = NULL;
            const leo2_result result = leo2_context_create(&options, &context);
            if (result == LEO2_UNSUPPORTED)
                continue;
            require_success(result, "explicit-backend context create");
            require(context != NULL, "context create returned null");
            require(leo2_context_backend(context) == backends[i],
                "explicit context selected a different backend");
            run_backend_suite(context, &counts);
            leo2_context_destroy(context);
            ++counts.backends;
        }
        require(counts.backends != 0, "no runtime backend was available");
        std::cout << "high/low duality test passed: parity_cases="
                  << counts.parity_cases
                  << " decode_cases=" << counts.decode_cases
                  << " restored_shards=" << counts.restored_shards
                  << " exhaustive_patterns=" << counts.exhaustive_patterns
                  << " random_patterns=" << counts.random_patterns
                  << " parity_rebuilds=" << counts.parity_rebuilds
                  << " backends=" << counts.backends << '\n';
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "high/low duality test failed: " << error.what() << '\n';
        return 1;
    }
}
