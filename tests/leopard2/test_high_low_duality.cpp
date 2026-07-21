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
#include "leopard2.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <new>
#include <sstream>
#include <stdexcept>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

#ifndef LEO2_ENABLE_TEST_HOOKS
#error "high/low duality test requires LEO2_ENABLE_TEST_HOOKS"
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

    Counts()
        : parity_cases(0)
        , decode_cases(0)
        , restored_shards(0)
        , exhaustive_patterns(0)
        , random_patterns(0)
        , backends(0)
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
    bool translated,
    bool test_plan_immutability)
{
    require_success(leo2_test_codec_set_decode_mode(codec,
        translated ? LEO2_TEST_DECODE_FORCE_TRANSLATED_LOW
                   : LEO2_TEST_DECODE_AUTO),
        "set decode mode");
    leo2_decode_plan* plan = NULL;
    require_success(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), "decode plan create");
    require(plan != NULL, "decode plan create returned null");
    require(leo2_test_decode_plan_uses_translated_low(plan) ==
            (translated ? 1 : 0),
        "decode plan did not capture the requested translation mode");
    if (test_plan_immutability)
    {
        require_success(leo2_test_codec_set_decode_mode(
            codec, LEO2_TEST_DECODE_AUTO), "reset codec decode mode");
        require(leo2_test_decode_plan_uses_translated_low(plan) == 1,
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
    AlignedBuffer scratch(scratch_bytes);
    require_success(leo2_decode_plan_execute(plan, bytes, &original_input[0],
        &recovery_input[0], &output[0], scratch.data(), scratch_bytes),
        "decode execute");
    leo2_decode_plan_destroy(plan);
    return restored;
}

void require_recovery(
    const Shards& originals,
    const Shards& translated,
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
        require(ordinary_high[i] == originals[i],
            "ordinary Algorithm 5 restored the wrong original");
        require(ordinary_low[i] == originals[i],
            "ordinary low-profile decoder restored the wrong original");
        require(translated[i] == ordinary_high[i] &&
                translated[i] == ordinary_low[i],
            "high/low recovery results differ");
        ++counts->restored_shards;
    }
}

void run_pattern(
    leo2_codec* translated_high,
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
        original_present, recovery_present, bytes, true,
        test_plan_immutability);
    const Shards high = decode(ordinary_high, originals, recovery,
        original_present, recovery_present, bytes, false, false);
    const Shards low = decode(ordinary_low, originals, recovery,
        original_present, recovery_present, bytes, false, false);
    require_recovery(originals, translated, high, low, original_present,
        counts);
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
            run_pattern(translated_high, ordinary_high, ordinary_low,
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
            run_pattern(translated_high, ordinary_high, ordinary_low,
                originals, high_recovery, original_present, recovery_present,
                bytes, !tested_immutability, counts);
            tested_immutability = true;
            ++counts->random_patterns;
        }
    }

    leo2_codec_destroy(ordinary_low);
    leo2_codec_destroy(ordinary_high);
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
        static_cast<leo2_test_decode_mode>(99)), LEO2_INVALID_ARGUMENT,
        "reject unknown decode mode");
    leo2_codec_destroy(unequal);

    leo2_codec* forced_generic = create_codec(context, 8, 8,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        LEO2_CODEC_FORCE_GENERIC_DECODE);
    require_result(leo2_test_codec_set_decode_mode(forced_generic,
        LEO2_TEST_DECODE_FORCE_TRANSLATED_LOW), LEO2_INVALID_ARGUMENT,
        "reject translation with forced generic decoder");
    leo2_codec_destroy(forced_generic);
}

void run_backend_suite(leo2_context* context, Counts* counts)
{
    test_rejection(context);
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
        {33, 64}, {63, 33}, {65, 128}, {127, 65}, {128, 127}
    };
    const size_t gf8_bytes[] = { 1, 7, 31, 63, 64, 65, 193 };
    for (size_t i = 0; i < sizeof(gf8_cases) / sizeof(gf8_cases[0]); ++i)
    {
        run_case(context, gf8_cases[i][0], gf8_cases[i][1],
            LEO2_FIELD_GF8, gf8_bytes[i %
                (sizeof(gf8_bytes) / sizeof(gf8_bytes[0]))],
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
            LEO2_BACKEND_NEON
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
                  << " backends=" << counts.backends << '\n';
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "high/low duality test failed: " << error.what() << '\n';
        return 1;
    }
}
