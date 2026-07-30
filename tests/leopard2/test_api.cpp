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

#include "direct_oracle.h"
#include "LeopardFF8.h"
#include "LeopardFF16.h"
#include "Leopard2Dispatch.h"
#include "Leopard2Direct.h"
#include "leopard.h"
#include "leopard2.h"

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

#ifndef LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE
#define LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE 0
#endif

#ifndef LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT
#define LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT 0
#endif

namespace {

using leopard2_test::BinaryField;
using leopard2_test::Element;
using leopard2_test::ProfileLayout;

struct TestCounts
{
    uint64_t high_compatibility;
    uint64_t low_oracle_symbols;
    uint64_t recovered_shards;
    uint64_t tail_cases;
    uint64_t plan_executions;

    TestCounts()
        : high_compatibility(0)
        , low_oracle_symbols(0)
        , recovered_shards(0)
        , tail_cases(0)
        , plan_executions(0)
    {}
};

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void require_result(leo2_result result, const std::string& operation)
{
    if (result != LEO2_SUCCESS)
    {
        std::ostringstream stream;
        stream << operation << ": " << leo2_result_string(result)
               << " (" << static_cast<int>(result) << ")";
        throw std::runtime_error(stream.str());
    }
}

struct AlignedBuffer
{
    void* data;
    size_t bytes;

    explicit AlignedBuffer(size_t size)
        : data(NULL)
        , bytes(size)
    {
        if (size != 0)
        {
#if defined(_MSC_VER)
            data = _aligned_malloc(size, leo2_scratch_alignment());
            if (data == NULL)
                throw std::bad_alloc();
#else
            const int error = posix_memalign(&data, leo2_scratch_alignment(), size);
            if (error != 0)
                throw std::bad_alloc();
#endif
            memset(data, 0, size);
        }
    }

    ~AlignedBuffer()
    {
#if defined(_MSC_VER)
        _aligned_free(data);
#else
        free(data);
#endif
    }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);
};

typedef std::vector<std::vector<uint8_t> > Shards;

std::vector<uint8_t> compact_pack_gf16(const std::vector<uint8_t>& input)
{
    require(!input.empty() && (input.size() & 1u) == 0,
        "GF16 compact pack requires a positive even byte count");
    const size_t rounded = (input.size() + 63u) & ~static_cast<size_t>(63u);
    std::vector<uint8_t> output(rounded, 0);
    const size_t complete = input.size() & ~static_cast<size_t>(63u);
    std::copy(input.begin(), input.begin() + complete, output.begin());
    const size_t symbols = (input.size() - complete) / 2;
    if (symbols != 0)
    {
        std::copy(input.begin() + complete, input.begin() + complete + symbols,
            output.begin() + complete);
        std::copy(input.begin() + complete + symbols, input.end(),
            output.begin() + complete + 32);
    }
    return output;
}

std::vector<uint8_t> compact_gather_gf16(
    const std::vector<uint8_t>& input,
    size_t bytes)
{
    require(bytes != 0 && (bytes & 1u) == 0 && input.size() >= bytes,
        "GF16 compact gather requires a positive even byte count");
    std::vector<uint8_t> output(bytes, 0);
    const size_t complete = bytes & ~static_cast<size_t>(63u);
    std::copy(input.begin(), input.begin() + complete, output.begin());
    const size_t symbols = (bytes - complete) / 2;
    if (symbols != 0)
    {
        std::copy(input.begin() + complete, input.begin() + complete + symbols,
            output.begin() + complete);
        std::copy(input.begin() + complete + 32,
            input.begin() + complete + 32 + symbols,
            output.begin() + complete + symbols);
    }
    return output;
}

Shards make_originals(unsigned count, size_t bytes, uint64_t seed)
{
    Shards shards(count, std::vector<uint8_t>(bytes));
    uint64_t state = seed;
    for (unsigned i = 0; i < count; ++i)
        for (size_t j = 0; j < bytes; ++j)
        {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            shards[i][j] = static_cast<uint8_t>(state + i * 29u + j * 131u);
        }
    return shards;
}

#if LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE != 0
uint64_t shards_digest(const Shards& shards)
{
    uint64_t digest = UINT64_C(1469598103934665603);
    for (size_t shard = 0; shard < shards.size(); ++shard)
        for (size_t byte_i = 0; byte_i < shards[shard].size(); ++byte_i)
        {
            digest ^= shards[shard][byte_i];
            digest *= UINT64_C(1099511628211);
        }
    return digest;
}
#endif

std::vector<const void*> const_pointers(const Shards& shards)
{
    std::vector<const void*> pointers(shards.size());
    for (size_t i = 0; i < shards.size(); ++i)
        pointers[i] = &shards[i][0];
    return pointers;
}

std::vector<void*> mutable_pointers(Shards& shards)
{
    std::vector<void*> pointers(shards.size());
    for (size_t i = 0; i < shards.size(); ++i)
        pointers[i] = &shards[i][0];
    return pointers;
}

leo2_codec* make_codec(
    leo2_context* context,
    unsigned k,
    unsigned r,
    leo2_profile profile,
    leo2_field field)
{
    leo2_codec* codec = NULL;
    require_result(leo2_codec_create(
        context, k, r, profile, field, NULL, &codec), "codec create");
    require(codec != NULL, "codec create returned a null codec");
    return codec;
}

Shards encode_new(const leo2_codec* codec, const Shards& original, size_t bytes)
{
    const unsigned recovery_count = leo2_codec_recovery_count(codec);
    Shards recovery(recovery_count, std::vector<uint8_t>(bytes));
    std::vector<const void*> original_ptrs = const_pointers(original);
    std::vector<void*> recovery_ptrs = mutable_pointers(recovery);
    size_t scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        "encode scratch query");
    AlignedBuffer scratch(scratch_bytes);
    require_result(leo2_encode(
        codec, bytes, &original_ptrs[0], &recovery_ptrs[0],
        scratch.data, scratch.bytes), "encode");
    return recovery;
}

Shards encode_legacy(const Shards& original, unsigned recovery_count, size_t bytes)
{
    const size_t rounded = (bytes + 63u) & ~static_cast<size_t>(63u);
    Shards padded(original.size(), std::vector<uint8_t>(rounded, 0));
    for (size_t i = 0; i < original.size(); ++i)
        memcpy(&padded[i][0], &original[i][0], bytes);
    std::vector<const void*> original_ptrs = const_pointers(padded);
    const unsigned work_count = leo_encode_work_count(
        static_cast<unsigned>(original.size()), recovery_count);
    Shards work(work_count, std::vector<uint8_t>(rounded, 0));
    std::vector<void*> work_ptrs = mutable_pointers(work);
    require(leo_encode(
        rounded, static_cast<unsigned>(original.size()), recovery_count,
        work_count, &original_ptrs[0], &work_ptrs[0]) == Leopard_Success,
        "legacy encode failed");
    Shards recovery(recovery_count, std::vector<uint8_t>(bytes));
    for (unsigned i = 0; i < recovery_count; ++i)
        memcpy(&recovery[i][0], &work[i][0], bytes);
    return recovery;
}

Shards encode_legacy_gf16_compact(
    const Shards& original,
    unsigned recovery_count,
    size_t bytes)
{
    const size_t rounded = (bytes + 63u) & ~static_cast<size_t>(63u);
    Shards packed(original.size());
    for (size_t i = 0; i < original.size(); ++i)
        packed[i] = compact_pack_gf16(original[i]);
    std::vector<const void*> original_ptrs = const_pointers(packed);
    const unsigned work_count = leo_encode_work_count(
        static_cast<unsigned>(original.size()), recovery_count);
    Shards work(work_count, std::vector<uint8_t>(rounded, 0));
    std::vector<void*> work_ptrs = mutable_pointers(work);
    require(leo_encode(rounded, static_cast<unsigned>(original.size()),
        recovery_count, work_count, &original_ptrs[0], &work_ptrs[0]) ==
        Leopard_Success, "legacy compact GF16 encode failed");
    Shards recovery(recovery_count);
    for (unsigned i = 0; i < recovery_count; ++i)
        recovery[i] = compact_gather_gf16(work[i], bytes);
    return recovery;
}

void test_profile_metadata(leo2_context* context)
{
    leo2_codec* high = make_codec(
        context, 193, 64, LEO2_PROFILE_AUTO, LEO2_FIELD_AUTO);
    require(leo2_codec_profile(high) == LEO2_PROFILE_LEGACY_HIGH_V1,
        "AUTO did not choose high for R <= K");
    require(leo2_codec_field(high) == LEO2_FIELD_GF16,
        "high parent field boundary was incorrect");
    require(leo2_codec_padded_side(high) == 64 && leo2_codec_parent_count(high) == 512,
        "high parent layout was incorrect");
    leo2_codec_destroy(high);

    leo2_codec* low = make_codec(
        context, 17, 40, LEO2_PROFILE_AUTO, LEO2_FIELD_AUTO);
    require(leo2_codec_profile(low) == LEO2_PROFILE_LOW_V1,
        "AUTO did not choose low for R > K");
    require(leo2_codec_padded_side(low) == 32 && leo2_codec_parent_count(low) == 128,
        "low parent layout was incorrect");
    leo2_codec_destroy(low);
}

void compare_high_with_legacy(
    leo2_context* context,
    unsigned k,
    unsigned r,
    const std::vector<size_t>& byte_counts,
    TestCounts* counts)
{
    leo2_codec* codec = make_codec(
        context, k, r, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_AUTO);
    for (size_t case_i = 0; case_i < byte_counts.size(); ++case_i)
    {
        const size_t bytes = byte_counts[case_i];
        const Shards original = make_originals(k, bytes, 0x12345678ULL + k * 97u + r * 13u);
        const Shards expected = encode_legacy(original, r, bytes);
        const Shards actual = encode_new(codec, original, bytes);
        require(actual == expected, "legacy-high parity mismatch");
        counts->high_compatibility += static_cast<uint64_t>(r) * bytes;
        if ((bytes & 63u) != 0)
            ++counts->tail_cases;

        if (r > 2)
        {
            std::vector<const void*> original_ptrs = const_pointers(original);
            Shards subset(r, std::vector<uint8_t>(bytes, 0xa5));
            std::vector<void*> subset_ptrs(r, NULL);
            subset_ptrs[1] = &subset[1][0];
            subset_ptrs[r - 1] = &subset[r - 1][0];
            size_t scratch_bytes = 0;
            require_result(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
                "subset scratch query");
            AlignedBuffer scratch(scratch_bytes);
            require_result(leo2_encode(
                codec, bytes, &original_ptrs[0], &subset_ptrs[0],
                scratch.data, scratch.bytes), "subset encode");
            require(subset[1] == expected[1] && subset[r - 1] == expected[r - 1],
                "requested parity subset mismatch");
        }
    }
    leo2_codec_destroy(codec);
}

void compare_avx2_fused_high_with_legacy(TestCounts* counts)
{
    leo2_context_options options;
    memset(&options, 0, sizeof(options));
    options.struct_size = sizeof(options);
    options.backend = LEO2_BACKEND_AVX2;
    options.thread_count = 1;
    leo2_context* context = NULL;
    const leo2_result result = leo2_context_create(&options, &context);
    require(result == LEO2_SUCCESS || result == LEO2_UNSUPPORTED,
        "explicit AVX2 compatibility context returned an unexpected result");
    if (result == LEO2_UNSUPPORTED)
    {
        require(context == NULL,
            "unsupported AVX2 compatibility context was not cleared");
        return;
    }
    require(context != NULL &&
            leo2_context_backend(context) == LEO2_BACKEND_AVX2,
        "explicit AVX2 compatibility context selected the wrong backend");

    // Mask-specialized inverse input boundaries: three live rows in the
    // shortened first T=8 block, one live row in the second block, and the
    // promoted large-shard T=16 boundary.  The +1 byte cases prove that the
    // aligned prefix and arbitrary public tail both retain legacy parity.
    compare_high_with_legacy(
        context, 7, 7, std::vector<size_t>{4096, 4097, 65536}, counts);
    compare_high_with_legacy(
        context, 9, 8, std::vector<size_t>{2048, 2049, 4097}, counts);
    compare_high_with_legacy(
        context, 15, 15, std::vector<size_t>{65536, 65537}, counts);

    // Exercise the lower and upper recovery-count edges of the promoted T
    // bins.  Ragged lengths cover the aligned fused pass followed by an
    // independently padded tail.  The T=64 high-rate case also covers a
    // partial final message block.
    compare_high_with_legacy(
        context, 16, 9, std::vector<size_t>{1089, 4097}, counts);
    compare_high_with_legacy(
        context, 16, 16, std::vector<size_t>{1089}, counts);
    compare_high_with_legacy(
        context, 32, 17, std::vector<size_t>{1089, 4097}, counts);
    compare_high_with_legacy(
        context, 32, 32, std::vector<size_t>{1089}, counts);
    compare_high_with_legacy(
        context, 64, 33, std::vector<size_t>{1089, 4097}, counts);
    compare_high_with_legacy(
        context, 190, 64, std::vector<size_t>{4096}, counts);
    compare_high_with_legacy(
        context, 127, 65, std::vector<size_t>{1089, 4097}, counts);
    compare_high_with_legacy(
        context, 128, 128, std::vector<size_t>{1089}, counts);
    leo2_context_destroy(context);
}

void compare_high_gf16_compact_with_legacy(
    leo2_context* context,
    unsigned k,
    unsigned r,
    size_t bytes,
    TestCounts* counts)
{
    leo2_codec* codec = make_codec(
        context, k, r, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16);
    const Shards original = make_originals(
        k, bytes, 0x16a17d5ULL + bytes * 19u);
    const Shards actual = encode_new(codec, original, bytes);
    const Shards expected = encode_legacy_gf16_compact(original, r, bytes);
    require(actual == expected,
        "compact GF16 high parity differs from packed legacy ALTMAP");
    counts->high_compatibility += static_cast<uint64_t>(r) * bytes;
    leo2_codec_destroy(codec);
}

void compare_low_gf8_with_oracle(
    leo2_context* context,
    unsigned k,
    unsigned r,
    size_t bytes,
    const BinaryField& field,
    TestCounts* counts)
{
    leo2_codec* codec = make_codec(
        context, k, r, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8);
    const Shards original = make_originals(k, bytes, 0x91e10da5ULL + k * 41u + r);
    const Shards actual = encode_new(codec, original, bytes);
    const ProfileLayout layout = leopard2_test::make_profile_layout(
        leopard2_test::kLow, k, r);
    for (size_t offset = 0; offset < bytes; ++offset)
    {
        std::vector<Element> message(k);
        for (unsigned i = 0; i < k; ++i)
            message[i] = original[i][offset];
        const std::vector<Element> expected = leopard2_test::direct_encode(
            field, layout, message);
        for (unsigned i = 0; i < r; ++i)
        {
            require(actual[i][offset] == expected[k + i],
                "low GF8 parity differs from direct interpolation");
            ++counts->low_oracle_symbols;
        }
    }
    leo2_codec_destroy(codec);
}

void compare_low_gf16_with_oracle(
    leo2_context* context,
    const BinaryField& field,
    size_t bytes,
    TestCounts* counts)
{
    const unsigned k = 3;
    const unsigned r = 5;
    leo2_codec* codec = make_codec(
        context, k, r, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16);
    const Shards original = make_originals(k, bytes, 0x6789abcdefULL);
    const Shards actual = encode_new(codec, original, bytes);
    const ProfileLayout layout = leopard2_test::make_profile_layout(
        leopard2_test::kLow, k, r);
    const size_t rounded = (bytes + 63u) & ~static_cast<size_t>(63u);
    Shards padded(k, std::vector<uint8_t>(rounded, 0));
    for (unsigned i = 0; i < k; ++i)
        padded[i] = compact_pack_gf16(original[i]);
    Shards expected(r, std::vector<uint8_t>(rounded, 0));
    for (size_t tile = 0; tile < rounded; tile += 64)
    {
        for (size_t lane = 0; lane < 32; ++lane)
        {
            std::vector<Element> message(k);
            for (unsigned i = 0; i < k; ++i)
                message[i] = static_cast<Element>(padded[i][tile + lane] |
                    (static_cast<unsigned>(padded[i][tile + 32 + lane]) << 8));
            const std::vector<Element> codeword = leopard2_test::direct_encode(
                field, layout, message);
            for (unsigned i = 0; i < r; ++i)
            {
                expected[i][tile + lane] = static_cast<uint8_t>(codeword[k + i]);
                expected[i][tile + 32 + lane] = static_cast<uint8_t>(codeword[k + i] >> 8);
            }
        }
    }
    for (unsigned i = 0; i < r; ++i)
    {
        require(actual[i] == compact_gather_gf16(expected[i], bytes),
            "low GF16 parity differs from direct interpolation");
        counts->low_oracle_symbols += bytes;
    }
    leo2_codec_destroy(codec);
}

void run_decode_case(
    leo2_context* context,
    unsigned k,
    unsigned r,
    leo2_profile profile,
    leo2_field field,
    size_t bytes,
    const std::vector<unsigned>& missing_originals,
    const std::vector<unsigned>& missing_recovery,
    TestCounts* counts)
{
    leo2_codec* codec = make_codec(context, k, r, profile, field);
    const Shards source = make_originals(k, bytes, 0xc001d00dULL + k * 11u + r);
    const Shards parity = encode_new(codec, source, bytes);
    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 1);
    for (size_t i = 0; i < missing_originals.size(); ++i)
        original_present[missing_originals[i]] = 0;
    for (size_t i = 0; i < missing_recovery.size(); ++i)
        recovery_present[missing_recovery[i]] = 0;

    leo2_decode_plan* plan = NULL;
    require_result(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), "decode plan create");
    require(leo2_decode_plan_missing_original_count(plan) == missing_originals.size(),
        "decode plan missing-original count mismatch");

    size_t scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(plan, bytes, &scratch_bytes),
        "decode scratch query");
    AlignedBuffer scratch(scratch_bytes);
    Shards restored(k, std::vector<uint8_t>(bytes, 0));
    std::vector<const void*> original_ptrs(k, NULL);
    std::vector<const void*> recovery_ptrs(r, NULL);
    std::vector<void*> restored_ptrs(k, NULL);
    for (unsigned i = 0; i < k; ++i)
    {
        if (original_present[i])
            original_ptrs[i] = &source[i][0];
        else
            restored_ptrs[i] = &restored[i][0];
    }
    for (unsigned i = 0; i < r; ++i)
        if (recovery_present[i])
            recovery_ptrs[i] = &parity[i][0];

    for (unsigned repeat = 0; repeat < 2; ++repeat)
    {
        require_result(leo2_decode_plan_execute(
            plan, bytes, &original_ptrs[0], &recovery_ptrs[0], &restored_ptrs[0],
            scratch.data, scratch.bytes), "decode plan execute");
        ++counts->plan_executions;
        for (size_t i = 0; i < missing_originals.size(); ++i)
        {
            const unsigned index = missing_originals[i];
            if (restored[index] != source[index])
            {
                size_t first_difference = 0;
                while (first_difference < bytes &&
                       restored[index][first_difference] == source[index][first_difference])
                    ++first_difference;
                std::ostringstream stream;
                stream << "restored original mismatch: K=" << k << " R=" << r
                       << " profile=" << static_cast<int>(profile)
                       << " field=" << static_cast<int>(field)
                       << " bytes=" << bytes << " original=" << index
                       << " repeat=" << repeat << " first_difference=" << first_difference;
                throw std::runtime_error(stream.str());
            }
            ++counts->recovered_shards;
        }
    }

    if (!missing_originals.empty())
    {
        leo2_codec_options generic_options;
        memset(&generic_options, 0, sizeof(generic_options));
        generic_options.struct_size = sizeof(generic_options);
        generic_options.flags = LEO2_CODEC_FORCE_GENERIC_DECODE;
        leo2_codec* generic_codec = NULL;
        require_result(leo2_codec_create(context, k, r, profile, field,
            &generic_options, &generic_codec), "generic codec create");
        leo2_decode_plan* generic_plan = NULL;
        require_result(leo2_decode_plan_create(generic_codec, &original_present[0],
            &recovery_present[0], &generic_plan), "generic plan create");
        size_t generic_scratch_bytes = 0;
        require_result(leo2_decode_plan_scratch_size(
            generic_plan, bytes, &generic_scratch_bytes),
            "generic scratch query");
        AlignedBuffer generic_scratch(generic_scratch_bytes);
        Shards generic_restored(k, std::vector<uint8_t>(bytes, 0));
        std::vector<void*> generic_restored_ptrs(k, NULL);
        for (size_t i = 0; i < missing_originals.size(); ++i)
            generic_restored_ptrs[missing_originals[i]] =
                &generic_restored[missing_originals[i]][0];
        require_result(leo2_decode_plan_execute(generic_plan, bytes,
            &original_ptrs[0], &recovery_ptrs[0], &generic_restored_ptrs[0],
            generic_scratch.data, generic_scratch.bytes),
            "generic plan execute");
        for (size_t i = 0; i < missing_originals.size(); ++i)
        {
            const unsigned index = missing_originals[i];
            require(generic_restored[index] == restored[index],
                "selected decoder and generic recovery differ");
        }
        leo2_decode_plan_destroy(generic_plan);
        leo2_codec_destroy(generic_codec);

        if (profile == LEO2_PROFILE_LEGACY_HIGH_V1 &&
            field == LEO2_FIELD_GF8 && k == 128 && r == 128 &&
            missing_originals.size() == 128 && bytes >= 256 && bytes <= 1024 * 1024)
        {
            leo2_codec_options specialized_options;
            memset(&specialized_options, 0, sizeof(specialized_options));
            specialized_options.struct_size = sizeof(specialized_options);
            specialized_options.flags = LEO2_CODEC_FORCE_SPECIALIZED_DECODE;
            leo2_codec* specialized_codec = NULL;
            require_result(leo2_codec_create(context, k, r, profile, field,
                &specialized_options, &specialized_codec),
                "specialized codec create");
            leo2_decode_plan* specialized_plan = NULL;
            require_result(leo2_decode_plan_create(specialized_codec,
                &original_present[0], &recovery_present[0], &specialized_plan),
                "specialized plan create");
            size_t specialized_scratch_bytes = 0;
            require_result(leo2_decode_plan_scratch_size(
                specialized_plan, bytes, &specialized_scratch_bytes),
                "specialized scratch query");
            AlignedBuffer specialized_scratch(specialized_scratch_bytes);
            Shards specialized_restored(k, std::vector<uint8_t>(bytes, 0));
            std::vector<void*> specialized_restored_ptrs(k, NULL);
            for (size_t i = 0; i < missing_originals.size(); ++i)
                specialized_restored_ptrs[missing_originals[i]] =
                    &specialized_restored[missing_originals[i]][0];
            require_result(leo2_decode_plan_execute(specialized_plan, bytes,
                &original_ptrs[0], &recovery_ptrs[0],
                &specialized_restored_ptrs[0], specialized_scratch.data,
                specialized_scratch.bytes), "specialized plan execute");
            for (size_t i = 0; i < missing_originals.size(); ++i)
            {
                const unsigned index = missing_originals[i];
                require(specialized_restored[index] == restored[index],
                    "automatic and forced-specialized recovery differ");
            }
            leo2_decode_plan_destroy(specialized_plan);
            leo2_codec_destroy(specialized_codec);
        }
    }
    if ((bytes & 63u) != 0)
        ++counts->tail_cases;

    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
}

void test_concurrent_preformatted_source_major_plan(
    leo2_context* context,
    TestCounts* counts)
{
#if defined(LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN)
    static const unsigned k = 65;
    static const unsigned r = 65;
    static const size_t bytes = 2049;
    static const unsigned worker_count = 4;
    static const unsigned repetitions = 4;
    const std::vector<unsigned> missing_originals = {
        0, 1, 7, 16, 32, 48, 63, 64
    };
    const std::vector<unsigned> missing_recovery = { 0, 3, 64 };
    const Shards source = make_originals(k, bytes, 0x5a17c0deULL);
    leo2_codec* codec = make_codec(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    const Shards parity = encode_new(codec, source, bytes);

    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 1);
    for (size_t i = 0; i < missing_originals.size(); ++i)
        original_present[missing_originals[i]] = 0;
    for (size_t i = 0; i < missing_recovery.size(); ++i)
        recovery_present[missing_recovery[i]] = 0;

    leo2_decode_plan* plan = NULL;
    require_result(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan),
        "concurrent source-major plan create");
    require(leo2_test_decode_plan_direct_source_rows(plan) != 0,
        "concurrent source-major plan was not preformatted");
    size_t scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(
        plan, bytes, &scratch_bytes),
        "concurrent source-major scratch query");

    std::atomic<unsigned> ready(0);
    std::atomic<bool> go(false);
    std::atomic<unsigned> failures(0);
    std::vector<std::thread> workers;
    workers.reserve(worker_count);
    for (unsigned worker = 0; worker < worker_count; ++worker)
    {
        workers.push_back(std::thread([&]() {
            bool announced_ready = false;
            try
            {
                AlignedBuffer scratch(scratch_bytes);
                Shards restored(k, std::vector<uint8_t>(bytes, 0));
                std::vector<const void*> original_ptrs(k, NULL);
                std::vector<const void*> recovery_ptrs(r, NULL);
                std::vector<void*> restored_ptrs(k, NULL);
                for (unsigned original = 0; original < k; ++original)
                {
                    if (original_present[original])
                        original_ptrs[original] = &source[original][0];
                    else
                        restored_ptrs[original] = &restored[original][0];
                }
                for (unsigned recovery_index = 0;
                     recovery_index < r; ++recovery_index)
                {
                    if (recovery_present[recovery_index])
                    {
                        recovery_ptrs[recovery_index] =
                            &parity[recovery_index][0];
                    }
                }

                ready.fetch_add(1, std::memory_order_release);
                announced_ready = true;
                while (!go.load(std::memory_order_acquire))
                    std::this_thread::yield();
                for (unsigned repeat = 0; repeat < repetitions; ++repeat)
                {
                    require_result(leo2_decode_plan_execute(plan, bytes,
                        &original_ptrs[0], &recovery_ptrs[0],
                        &restored_ptrs[0], scratch.data, scratch.bytes),
                        "concurrent source-major plan execute");
                    for (size_t missing = 0;
                         missing < missing_originals.size(); ++missing)
                    {
                        const unsigned index = missing_originals[missing];
                        if (restored[index] != source[index])
                            throw std::runtime_error(
                                "concurrent source-major recovery mismatch");
                    }
                }
            }
            catch (...)
            {
                if (!announced_ready)
                    ready.fetch_add(1, std::memory_order_release);
                failures.fetch_add(1, std::memory_order_relaxed);
            }
        }));
    }
    while (ready.load(std::memory_order_acquire) != worker_count)
        std::this_thread::yield();
    go.store(true, std::memory_order_release);
    for (size_t worker = 0; worker < workers.size(); ++worker)
        workers[worker].join();
    require(failures.load(std::memory_order_relaxed) == 0,
        "concurrent immutable source-major plan execution failed");
    counts->plan_executions += worker_count * repetitions;
    counts->recovered_shards += worker_count * repetitions *
        missing_originals.size();
    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
#else
    (void)context;
    (void)counts;
#endif
}

void test_no_loss_no_op(leo2_context* context)
{
    leo2_codec* codec = make_codec(
        context, 9, 7, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    std::vector<uint8_t> original_present(9, 1);
    std::vector<uint8_t> recovery_present(7, 0);
    leo2_decode_plan* plan = NULL;
    require_result(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), "no-loss plan create");
    size_t scratch_bytes = 99;
    require_result(leo2_decode_plan_scratch_size(plan, 17, &scratch_bytes),
        "no-loss scratch query");
    require(scratch_bytes == 0, "no-loss plan unexpectedly requires scratch");
    scratch_bytes = 99;
    require_result(leo2_decode_plan_scratch_size(
        plan, UINT64_MAX, &scratch_bytes),
        "no-loss maximum-byte scratch query");
    require(scratch_bytes == 0,
        "no-loss maximum-byte query unexpectedly requires scratch");
    scratch_bytes = 99;
    require(leo2_decode_plan_scratch_size(plan, 0, &scratch_bytes) ==
            LEO2_INVALID_ARGUMENT,
        "no-loss zero-byte scratch query did not reject the length");
    require(scratch_bytes == 0,
        "rejected no-loss zero-byte query did not clear output");
    require_result(leo2_decode_plan_execute(
        plan, 17, NULL, NULL, NULL, NULL, 0), "no-loss no-op execute");
    require_result(leo2_decode_plan_execute(
        plan, 0, NULL, NULL, NULL, NULL, 0),
        "no-loss zero-byte no-op execute");
    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
}

void test_direct_repair_dispatch_bounds(leo2_context* context)
{
    static const bool experimental_small_direct =
        LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE != 0;
    static const bool experimental_general_one_loss =
#if defined(LEO2_GFNI_VARIANT)
        false;
#else
        LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT != 0;
#endif
    struct Case
    {
        unsigned k;
        unsigned r;
        leo2_profile profile;
        leo2_field field;
        size_t bytes;
        unsigned losses;
        bool expect_direct;
        bool avx2_only;
    };
    const Case cases[] = {
        { 16, 8, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          33, 4, true, false },
        { 5, 8, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          64, 5, experimental_small_direct, true },
        { 8, 8, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          65, 8, experimental_small_direct, true },
        { 16, 8, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          4096, 8, experimental_small_direct, true },
        { 16, 9, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          4096, 8, false, false },
        { 17, 8, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          33, 4, false, false },
        { 17, 17, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          1, 1, true, true },
        { 17, 17, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          1, 4, false, false },
        { 17, 32, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          64, 1, true, true },
        { 17, 33, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          64, 1, true, true },
        { 17, 128, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          64, 1, true, true },
        { 32, 17, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          63, 1, true, true },
        { 32, 16, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          64, 1, true, true },
        { 33, 33, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          65, 1, true, true },
        { 64, 64, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          1025, 1, true, true },
        { 16, 31, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16,
          66, 4, true, false },
        { 17, 31, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16,
          66, 4, false, false },
        { 65, 65, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          63, 8, true, true },
        { 65, 65, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          2047, 2, true, true },
        { 65, 65, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          2048, 2, true, true },
        { 65, 65, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          2047, 8, true, true },
        { 65, 65, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          2048, 8, true, true },
        { 65, 66, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          64, 1, true, true },
        { 65, 66, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          64, 4, true, true },
        { 65, 66, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          64, 8, true, true },
        { 65, 66, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          2047, 8, true, true },
        { 65, 66, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          2048, 8, true, true },
        { 65, 96, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          64, 1, true, true },
        { 65, 96, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          64, 4, true, true },
        { 65, 96, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          64, 8, true, true },
        { 65, 127, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          65, 8, true, true },
        { 65, 128, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          64, 1, true, true },
        { 65, 128, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          64, 4, true, true },
        { 65, 128, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          64, 8, true, true },
        { 65, 128, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          2047, 8, true, true },
        { 65, 128, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          2048, 8, true, true },
        { 65, 65, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          64, 9, false, false },
        { 66, 65, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          64, 8, false, false },
        { 66, 65, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          64, 1, true, true },
        { 128, 65, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          65, 1, true, true },
        { 128, 128, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          1024, 1, true, true },
        { 128, 128, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          1024, 4, false, false },
        { 128, 64, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          64, 1, true, true },
        { 65, 64, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          64, 1, true, true },
        { 191, 64, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          63, 1, true, true },
        { 192, 63, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          64, 1, true, true },
        { 192, 64, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          65, 1, true, true },
        { 192, 64, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          65, 2, false, false },
        { 31, 200, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8,
          33, 1, true, true },
        { 32, 32, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8,
          33, 1, false, false },
        { 65, 64, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          64, 8, false, false },
        { 65, 65, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16,
          64, 8, false, false },
        { 65, 65, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8,
          64, 8, false, false }
    };
    for (size_t case_i = 0; case_i < sizeof(cases) / sizeof(cases[0]); ++case_i)
    {
        const Case& test = cases[case_i];
        std::vector<uint8_t> original_present(test.k, 1);
        std::vector<uint8_t> recovery_present(test.r, 1);
        for (unsigned i = 0; i < test.losses; ++i)
            original_present[i] = 0;
        // Force deterministic selection to skip parity zero when enough
        // parity remains. A full-loss R=L case requires every parity.
        if (test.losses < test.r)
            recovery_present[0] = 0;

        leo2_codec* codec = make_codec(
            context, test.k, test.r, test.profile, test.field);
        leo2_decode_plan* plan = NULL;
        require_result(leo2_decode_plan_create(codec, &original_present[0],
            &recovery_present[0], &plan), "direct-dispatch plan create");
        size_t scratch_bytes = 0;
        require_result(leo2_decode_plan_scratch_size(plan, test.bytes, &scratch_bytes),
            "direct-dispatch scratch query");

        leo2_codec_options reference_options;
        memset(&reference_options, 0, sizeof(reference_options));
        reference_options.struct_size = sizeof(reference_options);
        reference_options.flags = LEO2_CODEC_FORCE_SPECIALIZED_DECODE;
        leo2_codec* reference_codec = NULL;
        require_result(leo2_codec_create(context, test.k, test.r, test.profile,
            test.field, &reference_options, &reference_codec),
            "direct-dispatch reference codec create");
        leo2_decode_plan* reference_plan = NULL;
        require_result(leo2_decode_plan_create(reference_codec,
            &original_present[0], &recovery_present[0], &reference_plan),
            "direct-dispatch reference plan create");
        size_t reference_scratch_bytes = 0;
        require_result(leo2_decode_plan_scratch_size(
            reference_plan, test.bytes, &reference_scratch_bytes),
            "direct-dispatch reference scratch query");
        uint32_t padded_k = 1;
        while (padded_k < test.k)
            padded_k <<= 1;
        uint32_t padded_r = 1;
        while (padded_r < test.r)
            padded_r <<= 1;
        const bool general_one_loss_shape =
            test.field == LEO2_FIELD_GF8 && test.k > 16 &&
            ((test.profile == LEO2_PROFILE_LEGACY_HIGH_V1 &&
              test.r > 1 && padded_k != padded_r) ||
             (test.profile == LEO2_PROFILE_LOW_V1 && test.r > test.k));
        const bool expect_direct = test.expect_direct &&
            (!general_one_loss_shape || experimental_general_one_loss) &&
            (!test.avx2_only ||
             leo2_context_backend(context) == LEO2_BACKEND_AVX2);
        leopard2_internal::DecodePathInfo path_info;
        require_result(leopard2_internal::GetDecodePlanPathInfo(
            plan, test.bytes, false, &path_info),
            "direct-dispatch path introspection");
        if ((path_info.path == leopard2_internal::kDecodePathDirect) !=
            expect_direct)
        {
            std::ostringstream message;
            message << "direct-repair path introspection disagrees for K="
                    << test.k << " R=" << test.r
                    << " profile=" << static_cast<unsigned>(test.profile)
                    << " field=" << static_cast<unsigned>(test.field)
                    << " losses=" << test.losses
                    << " expected=" << expect_direct
                    << " actual_path="
                    << static_cast<unsigned>(path_info.path);
            require(false, message.str());
        }
        if (expect_direct)
        {
            const bool expanded_k65_source_major =
                test.k == 65 && test.r >= 65 && test.r <= 128 &&
                test.losses >= 2 &&
                test.bytes >= 2048;
            const bool experimental_small_source_major =
                test.k >= 5 && test.k <= 16 &&
                test.r >= 5 && test.r <= 8 && test.losses >= 5 &&
                LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE == 2;
            const leopard2_internal::DirectRepairExecutor expected_executor =
                expanded_k65_source_major || experimental_small_source_major
                    ? leopard2_internal::kDirectRepairExecutorSourceMajor
                    : leopard2_internal::kDirectRepairExecutorOutputMajor;
            require(path_info.direct_executor == expected_executor,
                "direct-repair executor introspection selected wrong loop order");
        }
        else
        {
            require(path_info.direct_executor ==
                    leopard2_internal::kDirectRepairExecutorNone,
                "transform plan reported a direct-repair executor");
        }
#if defined(LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN)
        const bool expect_source_rows = expect_direct && test.k == 65 &&
            test.losses >= 2 && test.losses <= 8 &&
            test.profile == LEO2_PROFILE_LEGACY_HIGH_V1 &&
            test.field == LEO2_FIELD_GF8 &&
            leo2_context_backend(context) == LEO2_BACKEND_AVX2;
        require((leo2_test_decode_plan_direct_source_rows(plan) != 0) ==
                expect_source_rows,
            "direct-repair source-plan formatting boundary mismatch");
#else
        require(leo2_test_decode_plan_direct_source_rows(plan) == 0,
            "default build unexpectedly retained a source-major plan");
#endif
        require(expect_direct
                ? scratch_bytes < reference_scratch_bytes
                : scratch_bytes == reference_scratch_bytes,
            "direct-repair dispatch boundary selected the wrong scratch shape");
        leo2_decode_plan_destroy(reference_plan);
        leo2_codec_destroy(reference_codec);
        leo2_decode_plan_destroy(plan);
        leo2_codec_destroy(codec);
    }

    if (leo2_context_backend(context) == LEO2_BACKEND_AVX2)
    {
        std::vector<uint8_t> original_present(65, 1);
        std::vector<uint8_t> recovery_present(65, 1);
        for (unsigned i = 0; i < 8; ++i)
            original_present[i] = 0;
        recovery_present[0] = 0;

        leo2_codec* automatic = make_codec(context, 65, 65,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        leo2_decode_plan* automatic_plan = NULL;
        require_result(leo2_decode_plan_create(automatic,
            &original_present[0], &recovery_present[0], &automatic_plan),
            "forced-exclusion automatic plan create");
        size_t automatic_scratch_bytes = 0;
        require_result(leo2_decode_plan_scratch_size(
            automatic_plan, 64, &automatic_scratch_bytes),
            "forced-exclusion automatic scratch query");

        const uint32_t forced_flags[] = {
            LEO2_CODEC_FORCE_GENERIC_DECODE,
            LEO2_CODEC_FORCE_SPECIALIZED_DECODE
        };
        for (size_t i = 0;
             i < sizeof(forced_flags) / sizeof(forced_flags[0]); ++i)
        {
            leo2_codec_options options;
            memset(&options, 0, sizeof(options));
            options.struct_size = sizeof(options);
            options.flags = forced_flags[i];
            leo2_codec* forced = NULL;
            require_result(leo2_codec_create(context, 65, 65,
                LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
                &options, &forced), "forced-exclusion codec create");
            leo2_decode_plan* forced_plan = NULL;
            require_result(leo2_decode_plan_create(forced,
                &original_present[0], &recovery_present[0], &forced_plan),
                "forced-exclusion plan create");
            size_t forced_scratch_bytes = 0;
            require_result(leo2_decode_plan_scratch_size(
                forced_plan, 64, &forced_scratch_bytes),
                "forced-exclusion scratch query");
            require(forced_scratch_bytes > automatic_scratch_bytes,
                "forced decoder mode unexpectedly selected direct repair");
            leo2_decode_plan_destroy(forced_plan);
            leo2_codec_destroy(forced);
        }
        leo2_decode_plan_destroy(automatic_plan);
        leo2_codec_destroy(automatic);
    }
}

#if LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT
void test_generalized_one_loss_direct_repair_execution(
    leo2_context* context,
    TestCounts* counts)
{
    struct Case
    {
        unsigned k;
        unsigned r;
        unsigned missing_original;
        size_t bytes;
    };
    const Case cases[] = {
        { 17, 17, 0, 1 },
        { 17, 32, 16, 65 },
        { 32, 17, 31, 1025 },
        { 32, 32, 15, 64 },
        { 33, 33, 32, 63 },
        { 33, 64, 0, 65 },
        { 64, 33, 63, 1025 },
        { 64, 64, 31, 64 },
        { 66, 65, 65, 1 },
        { 66, 128, 0, 65 },
        { 128, 65, 127, 1025 },
        { 128, 128, 64, 1024 * 1024 }
    };
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i)
    {
        const Case& test = cases[i];
        run_decode_case(context, test.k, test.r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
            test.bytes, std::vector<unsigned>{test.missing_original},
            std::vector<unsigned>{0, test.r / 2, test.r - 1}, counts);
    }

    struct GeneralCase
    {
        unsigned k;
        unsigned r;
        leo2_profile profile;
        unsigned missing_original;
        size_t bytes;
    };
    const GeneralCase general_cases[] = {
        { 32, 16, LEO2_PROFILE_LEGACY_HIGH_V1, 0, 17 },
        { 65, 64, LEO2_PROFILE_LEGACY_HIGH_V1, 64, 33 },
        { 128, 64, LEO2_PROFILE_LEGACY_HIGH_V1, 63, 257 },
        { 191, 64, LEO2_PROFILE_LEGACY_HIGH_V1, 0, 1 },
        { 192, 63, LEO2_PROFILE_LEGACY_HIGH_V1, 191, 63 },
        { 192, 64, LEO2_PROFILE_LEGACY_HIGH_V1, 96, 65 },
        { 224, 32, LEO2_PROFILE_LEGACY_HIGH_V1, 223, 1025 },
        { 240, 16, LEO2_PROFILE_LEGACY_HIGH_V1, 0, 64 * 1024 },
        { 192, 64, LEO2_PROFILE_LEGACY_HIGH_V1, 191, 64 * 1024 + 1 },
        { 17, 31, LEO2_PROFILE_LOW_V1, 16, 7 },
        { 31, 200, LEO2_PROFILE_LOW_V1, 0, 33 },
        { 32, 224, LEO2_PROFILE_LOW_V1, 31, 1025 },
        { 127, 128, LEO2_PROFILE_LOW_V1, 63, 64 * 1024 },
        { 31, 200, LEO2_PROFILE_LOW_V1, 30, 64 * 1024 + 1 }
    };
    for (size_t i = 0;
         i < sizeof(general_cases) / sizeof(general_cases[0]); ++i)
    {
        const GeneralCase& test = general_cases[i];
        run_decode_case(context, test.k, test.r, test.profile,
            LEO2_FIELD_GF8, test.bytes,
            std::vector<unsigned>{test.missing_original},
            std::vector<unsigned>{0, test.r / 2, test.r - 1}, counts);
    }

    /*
        The ordinary cases leave many parity rows available, so deterministic
        received-subset selection tends to choose the same low row.  Leave
        exactly the first or last transmitted parity row available to bind
        both ends of the promoted high and low generator-coordinate maps.
    */
    struct ParityRowCase
    {
        unsigned k;
        unsigned r;
        leo2_profile profile;
        unsigned missing_original;
        unsigned available_recovery;
    };
    const ParityRowCase parity_row_cases[] = {
        { 254, 2, LEO2_PROFILE_LEGACY_HIGH_V1, 0, 0 },
        { 254, 2, LEO2_PROFILE_LEGACY_HIGH_V1, 253, 1 },
        { 192, 64, LEO2_PROFILE_LEGACY_HIGH_V1, 0, 0 },
        { 192, 64, LEO2_PROFILE_LEGACY_HIGH_V1, 191, 63 },
        { 31, 200, LEO2_PROFILE_LOW_V1, 0, 0 },
        { 31, 200, LEO2_PROFILE_LOW_V1, 30, 199 }
    };
    for (size_t i = 0;
         i < sizeof(parity_row_cases) / sizeof(parity_row_cases[0]); ++i)
    {
        const ParityRowCase& test = parity_row_cases[i];
        std::vector<unsigned> missing_recovery;
        missing_recovery.reserve(test.r - 1);
        for (unsigned recovery = 0; recovery < test.r; ++recovery)
        {
            if (recovery != test.available_recovery)
                missing_recovery.push_back(recovery);
        }
        run_decode_case(context, test.k, test.r, test.profile,
            LEO2_FIELD_GF8, 65,
            std::vector<unsigned>{test.missing_original},
            missing_recovery, counts);
    }
}

#if !defined(LEO2_GFNI_VARIANT)
void test_generalized_four_tiny_unaligned(
    leo2_context* context,
    TestCounts* counts)
{
    static const size_t byte_counts[] = {
        1, 2, 3, 4, 5, 6, 7, 8, 15, 16, 17, 23, 24, 25,
        31, 32, 33, 47, 48, 49, 55, 56, 57, 63, 64
    };
    struct Shape
    {
        unsigned k;
        unsigned r;
        leo2_profile profile;
    };
    static const Shape shapes[] = {
        { 240, 16, LEO2_PROFILE_LEGACY_HIGH_V1 },
        { 18, 16, LEO2_PROFILE_LEGACY_HIGH_V1 },
        { 17, 31, LEO2_PROFILE_LOW_V1 },
        { 31, 200, LEO2_PROFILE_LOW_V1 },
        { 32, 224, LEO2_PROFILE_LOW_V1 }
    };

    for (size_t shape_i = 0;
         shape_i < sizeof(shapes) / sizeof(shapes[0]); ++shape_i)
    {
        const Shape& shape = shapes[shape_i];
        const unsigned missing = shape.k - 1;
        leo2_codec* codec = make_codec(context, shape.k, shape.r,
            shape.profile, LEO2_FIELD_GF8);
        std::vector<uint8_t> original_present(shape.k, 1);
        std::vector<uint8_t> recovery_present(shape.r, 1);
        original_present[missing] = 0;
        leo2_decode_plan* plan = NULL;
        require_result(leo2_decode_plan_create(codec,
            original_present.data(), recovery_present.data(), &plan),
            "general four-tiny plan create");

        for (size_t byte_i = 0;
             byte_i < sizeof(byte_counts) / sizeof(byte_counts[0]); ++byte_i)
        {
            const size_t bytes = byte_counts[byte_i];
            leopard2_internal::DecodePathInfo path;
            require_result(leopard2_internal::GetDecodePlanPathInfo(
                plan, bytes, false, &path),
                "general four-tiny path query");
            require(path.path == leopard2_internal::kDecodePathDirect,
                "general four-tiny case did not retain direct repair");

            const Shards source = make_originals(shape.k, bytes,
                UINT64_C(0x4f54110000000000) +
                    shape_i * UINT64_C(0x100000) + bytes);
            const Shards parity = encode_new(codec, source, bytes);
            Shards original_storage(
                shape.k, std::vector<uint8_t>(bytes + 2, 0x5a));
            Shards recovery_storage(
                shape.r, std::vector<uint8_t>(bytes + 2, 0x6b));
            Shards output_storage(
                shape.k, std::vector<uint8_t>(bytes + 2, 0xa5));
            std::vector<const void*> original_input(shape.k, NULL);
            std::vector<const void*> recovery_input(shape.r, NULL);
            std::vector<void*> output(shape.k, NULL);
            for (unsigned original = 0; original < shape.k; ++original)
            {
                std::memcpy(&original_storage[original][1],
                    source[original].data(), bytes);
                if (original != missing)
                    original_input[original] =
                        &original_storage[original][1];
            }
            for (unsigned recovery = 0; recovery < shape.r; ++recovery)
            {
                std::memcpy(&recovery_storage[recovery][1],
                    parity[recovery].data(), bytes);
                recovery_input[recovery] =
                    &recovery_storage[recovery][1];
            }
            output[missing] = &output_storage[missing][1];
            size_t scratch_bytes = 0;
            require_result(leo2_decode_plan_scratch_size(
                plan, bytes, &scratch_bytes),
                "general four-tiny scratch query");
            AlignedBuffer scratch(scratch_bytes);

            leo2_test_reset_direct_four_tiny_calls();
            require_result(leo2_decode_plan_execute(plan, bytes,
                original_input.data(), recovery_input.data(), output.data(),
                scratch.data, scratch.bytes),
                "general four-tiny execute");
            const bool selected =
#if defined(LEO2_GFNI_VARIANT)
                false;
#else
                bytes != 0 && bytes <= 63 && bytes != 7;
#endif
            const uint64_t expected_calls =
                selected ? shape.k / 4U : 0;
            require(leo2_test_direct_four_tiny_calls() == expected_calls,
                "general four-tiny callback count differs from selector");
            require(std::memcmp(&output_storage[missing][1],
                    source[missing].data(), bytes) == 0,
                "general four-tiny output differs from source");
            require(output_storage[missing][0] == 0xa5 &&
                    output_storage[missing][bytes + 1] == 0xa5,
                "general four-tiny output changed a guard byte");
            for (unsigned original = 0; original < shape.k; ++original)
            {
                require(original_storage[original][0] == 0x5a &&
                        original_storage[original][bytes + 1] == 0x5a,
                    "general four-tiny changed an original guard byte");
            }
            for (unsigned recovery = 0; recovery < shape.r; ++recovery)
            {
                require(recovery_storage[recovery][0] == 0x6b &&
                        recovery_storage[recovery][bytes + 1] == 0x6b,
                    "general four-tiny changed a recovery guard byte");
            }
            ++counts->recovered_shards;
            ++counts->plan_executions;
            if ((bytes & 63U) != 0)
                ++counts->tail_cases;
        }
        leo2_decode_plan_destroy(plan);
        leo2_codec_destroy(codec);
    }
}
#endif
#endif

#if LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE != 0
void test_experimental_small_direct_unaligned(
    leo2_context* context,
    size_t bytes,
    TestCounts* counts)
{
    const unsigned k = 16;
    const unsigned r = 8;
    const unsigned losses = 5;
    leo2_codec* codec = make_codec(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    const Shards source = make_originals(
        k, bytes, UINT64_C(0x51a11d0000000000) + bytes);
    const Shards parity = encode_new(codec, source, bytes);
    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 1);
    for (unsigned i = 0; i < losses; ++i)
        original_present[i] = 0;
    leo2_decode_plan* plan = NULL;
    require_result(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), "unaligned small-direct plan create");
    size_t scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(
        plan, bytes, &scratch_bytes), "unaligned small-direct scratch query");
    AlignedBuffer scratch(scratch_bytes);

    Shards original_storage(k, std::vector<uint8_t>(bytes + 2, 0x5a));
    Shards recovery_storage(r, std::vector<uint8_t>(bytes + 2, 0x5a));
    Shards output_storage(k, std::vector<uint8_t>(bytes + 2, 0xa5));
    std::vector<const void*> original_input(k, NULL);
    std::vector<const void*> recovery_input(r, NULL);
    std::vector<void*> output(k, NULL);
    for (unsigned i = 0; i < k; ++i)
    {
        std::memcpy(&original_storage[i][1], source[i].data(), bytes);
        if (original_present[i])
            original_input[i] = &original_storage[i][1];
        else
            output[i] = &output_storage[i][1];
    }
    for (unsigned i = 0; i < r; ++i)
    {
        std::memcpy(&recovery_storage[i][1], parity[i].data(), bytes);
        recovery_input[i] = &recovery_storage[i][1];
    }
    require_result(leo2_decode_plan_execute(plan, bytes,
        original_input.data(), recovery_input.data(), output.data(),
        scratch.data, scratch.bytes), "unaligned small-direct execute");
    for (unsigned i = 0; i < losses; ++i)
    {
        require(std::memcmp(&output_storage[i][1], source[i].data(), bytes) == 0,
            "unaligned small-direct output differs from source");
        require(output_storage[i][0] == 0xa5 &&
                output_storage[i][bytes + 1] == 0xa5,
            "unaligned small-direct output changed a guard byte");
    }
    for (unsigned i = losses; i < k; ++i)
        require(original_storage[i][0] == 0x5a &&
                original_storage[i][bytes + 1] == 0x5a,
            "unaligned small-direct execution changed an input guard");
    for (unsigned i = 0; i < r; ++i)
        require(recovery_storage[i][0] == 0x5a &&
                recovery_storage[i][bytes + 1] == 0x5a,
            "unaligned small-direct execution changed a parity guard");
    counts->recovered_shards += losses;
    ++counts->plan_executions;
    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
}

void test_experimental_small_direct_batch_and_concurrent(
    leo2_context* context,
    TestCounts* counts)
{
    static const unsigned k = 8;
    static const unsigned r = 8;
    static const size_t bytes = 65;
    static const unsigned lanes = 2;
    const Shards source[lanes] = {
        make_originals(k, bytes, UINT64_C(0x82a1000000000001)),
        make_originals(k, bytes, UINT64_C(0x82b2000000000002))
    };
    require(shards_digest(source[0]) != shards_digest(source[1]),
        "small-direct batch sources did not have distinct digests");

    leo2_codec* codec = make_codec(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    const Shards parity[lanes] = {
        encode_new(codec, source[0], bytes),
        encode_new(codec, source[1], bytes)
    };
    require(shards_digest(parity[0]) != shards_digest(parity[1]),
        "small-direct batch parity did not have distinct digests");

    std::vector<uint8_t> original_present(k, 0);
    std::vector<uint8_t> recovery_present(r, 1);
    leo2_decode_plan* plan = NULL;
    require_result(leo2_decode_plan_create(codec, original_present.data(),
        recovery_present.data(), &plan),
        "small-direct shared plan create");
    size_t scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(
        plan, bytes, &scratch_bytes), "small-direct shared scratch query");
    leopard2_internal::DecodePathInfo path;
    require_result(leopard2_internal::GetDecodePlanPathInfo(
        plan, bytes, true, &path), "small-direct batch path introspection");
    const leopard2_internal::DirectRepairExecutor expected_executor =
        LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE == 2
            ? leopard2_internal::kDirectRepairExecutorSourceMajor
            : leopard2_internal::kDirectRepairExecutorOutputMajor;
    require(path.path == leopard2_internal::kDecodePathDirect &&
            path.direct_executor == expected_executor,
        "small-direct batch selected the wrong executor");

    std::vector<const void*> original_input(k, NULL);
    Shards recovery_storage[lanes] = {
        Shards(r, std::vector<uint8_t>(bytes + 2, 0x51)),
        Shards(r, std::vector<uint8_t>(bytes + 2, 0x62))
    };
    Shards batch_output[lanes] = {
        Shards(k, std::vector<uint8_t>(bytes + 2, 0xa1)),
        Shards(k, std::vector<uint8_t>(bytes + 2, 0xb2))
    };
    std::vector<const void*> recovery_input[lanes] = {
        std::vector<const void*>(r, NULL),
        std::vector<const void*>(r, NULL)
    };
    std::vector<void*> output[lanes] = {
        std::vector<void*>(k, NULL),
        std::vector<void*>(k, NULL)
    };
    for (unsigned lane = 0; lane < lanes; ++lane)
    {
        for (unsigned recovery = 0; recovery < r; ++recovery)
        {
            std::memcpy(&recovery_storage[lane][recovery][1],
                parity[lane][recovery].data(), bytes);
            recovery_input[lane][recovery] =
                &recovery_storage[lane][recovery][1];
        }
        for (unsigned original = 0; original < k; ++original)
            output[lane][original] = &batch_output[lane][original][1];
    }
    AlignedBuffer batch_scratch0(scratch_bytes);
    AlignedBuffer batch_scratch1(scratch_bytes);
    AlignedBuffer* batch_scratch[lanes] = {
        &batch_scratch0, &batch_scratch1
    };
    leo2_decode_batch_item items[lanes];
    std::memset(items, 0, sizeof(items));
    for (unsigned lane = 0; lane < lanes; ++lane)
    {
        items[lane].shard_bytes = bytes;
        items[lane].original = original_input.data();
        items[lane].recovery = recovery_input[lane].data();
        items[lane].restored_original = output[lane].data();
        items[lane].scratch = batch_scratch[lane]->data;
        items[lane].scratch_bytes = batch_scratch[lane]->bytes;
    }
    require_result(leo2_decode_plan_execute_batch(plan, items, lanes),
        "small-direct two-item batch execute");
    for (unsigned lane = 0; lane < lanes; ++lane)
    {
        const uint8_t recovery_guard =
            lane == 0 ? static_cast<uint8_t>(0x51)
                      : static_cast<uint8_t>(0x62);
        const uint8_t output_guard =
            lane == 0 ? static_cast<uint8_t>(0xa1)
                      : static_cast<uint8_t>(0xb2);
        for (unsigned recovery = 0; recovery < r; ++recovery)
            require(recovery_storage[lane][recovery][0] == recovery_guard &&
                    recovery_storage[lane][recovery][bytes + 1] ==
                        recovery_guard,
                "small-direct batch changed a parity guard");
        for (unsigned original = 0; original < k; ++original)
        {
            require(std::memcmp(&batch_output[lane][original][1],
                        source[lane][original].data(), bytes) == 0,
                "small-direct batch recovery mismatch");
            require(batch_output[lane][original][0] == output_guard &&
                    batch_output[lane][original][bytes + 1] == output_guard,
                "small-direct batch changed an output guard");
        }
    }

    Shards concurrent_output[lanes] = {
        Shards(k, std::vector<uint8_t>(bytes + 2, 0xc3)),
        Shards(k, std::vector<uint8_t>(bytes + 2, 0xd4))
    };
    std::vector<void*> concurrent_ptrs[lanes] = {
        std::vector<void*>(k, NULL), std::vector<void*>(k, NULL)
    };
    for (unsigned lane = 0; lane < lanes; ++lane)
        for (unsigned original = 0; original < k; ++original)
            concurrent_ptrs[lane][original] =
                &concurrent_output[lane][original][1];
    std::atomic<unsigned> ready(0);
    std::atomic<bool> go(false);
    std::vector<std::string> errors(lanes);
    std::vector<std::thread> workers;
    AlignedBuffer concurrent_scratch0(scratch_bytes);
    AlignedBuffer concurrent_scratch1(scratch_bytes);
    AlignedBuffer* concurrent_scratch[lanes] = {
        &concurrent_scratch0, &concurrent_scratch1
    };
    for (unsigned lane = 0; lane < lanes; ++lane)
    {
        workers.push_back(std::thread([&, lane]() {
            try
            {
                ready.fetch_add(1, std::memory_order_release);
                while (!go.load(std::memory_order_acquire))
                    std::this_thread::yield();
                for (unsigned repeat = 0; repeat < 4; ++repeat)
                    require_result(leo2_decode_plan_execute(
                        plan, bytes, original_input.data(),
                        recovery_input[lane].data(),
                        concurrent_ptrs[lane].data(),
                        concurrent_scratch[lane]->data,
                        concurrent_scratch[lane]->bytes),
                        "small-direct concurrent shared-plan execute");
            }
            catch (const std::exception& error)
            {
                errors[lane] = error.what();
            }
        }));
    }
    while (ready.load(std::memory_order_acquire) != lanes)
        std::this_thread::yield();
    go.store(true, std::memory_order_release);
    for (size_t lane = 0; lane < workers.size(); ++lane)
        workers[lane].join();
    for (unsigned lane = 0; lane < lanes; ++lane)
    {
        require(errors[lane].empty(),
            std::string("small-direct shared-plan lane failed: ") +
                errors[lane]);
        const uint8_t guard =
            lane == 0 ? static_cast<uint8_t>(0xc3)
                      : static_cast<uint8_t>(0xd4);
        for (unsigned original = 0; original < k; ++original)
        {
            require(std::memcmp(&concurrent_output[lane][original][1],
                        source[lane][original].data(), bytes) == 0,
                "small-direct concurrent recovery mismatch");
            require(concurrent_output[lane][original][0] == guard &&
                    concurrent_output[lane][original][bytes + 1] == guard,
                "small-direct concurrent execution changed an output guard");
        }
    }
    counts->recovered_shards += lanes * k * 5;
    counts->plan_executions += lanes * 5;
    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
}

void test_experimental_small_direct_repair_execution(
    leo2_context* context,
    TestCounts* counts)
{
    struct Shape
    {
        unsigned k;
        unsigned r;
        unsigned losses;
    };
    static const Shape shapes[] = {
        { 5, 5, 5 }, { 8, 8, 5 }, { 8, 8, 8 },
        { 9, 6, 6 }, { 12, 7, 7 }, { 14, 8, 8 },
        { 16, 6, 6 }, { 16, 8, 5 }, { 16, 8, 7 }, { 16, 8, 8 }
    };
    static const size_t boundary_bytes[] = {
        31, 32, 33, 63, 64, 65,
        2047, 2048, 2049, 4095, 4096, 4097,
        65535, 65536, 65537
    };
    for (size_t shape_i = 0;
         shape_i < sizeof(shapes) / sizeof(shapes[0]); ++shape_i)
    {
        const Shape& shape = shapes[shape_i];
        std::vector<unsigned> missing(shape.losses, 0);
        for (unsigned i = 0; i < shape.losses; ++i)
            missing[i] = i;
        for (size_t byte_i = 0;
             byte_i < sizeof(boundary_bytes) / sizeof(boundary_bytes[0]);
             ++byte_i)
        {
            run_decode_case(context, shape.k, shape.r,
                LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
                boundary_bytes[byte_i], missing, std::vector<unsigned>(),
                counts);
        }
    }

    std::vector<unsigned> all_missing(8, 0);
    for (unsigned i = 0; i < all_missing.size(); ++i)
        all_missing[i] = i;
    for (size_t bytes = 1; bytes <= 257; ++bytes)
    {
        run_decode_case(context, 8, 8,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
            bytes, all_missing, std::vector<unsigned>(), counts);
    }
    test_experimental_small_direct_unaligned(context, 65, counts);
    test_experimental_small_direct_unaligned(context, 2049, counts);
    test_experimental_small_direct_batch_and_concurrent(context, counts);
}
#endif

void test_concurrent_expanded_direct_repair_cache(
    leo2_context* context,
    TestCounts* counts)
{
    /*
        This must run before any other K=65 codec is created on this context:
        all lanes race the context-owned call_once cache, then compare every
        promoted loss count with an independently forced transform plan.
    */
    static const unsigned lane_count = 8;
    static const unsigned recovery_counts[lane_count] = {
        65, 66, 96, 127, 128, 65, 96, 128
    };
    static const unsigned missing_candidates[lane_count] = {
        0, 1, 7, 16, 32, 48, 63, 64
    };
    std::atomic<unsigned> ready(0);
    std::atomic<bool> go(false);
    std::vector<TestCounts> lane_counts(lane_count);
    std::vector<std::string> errors(lane_count);
    std::vector<std::thread> threads;
    threads.reserve(lane_count);

    for (unsigned lane = 0; lane < lane_count; ++lane)
    {
        threads.push_back(std::thread(
            [context, lane, &ready, &go, &lane_counts, &errors]() {
            bool announced_ready = false;
            try
            {
                const unsigned recovery_count = recovery_counts[lane];
                const unsigned loss_count = lane + 1;
                std::vector<unsigned> missing_originals(
                    missing_candidates,
                    missing_candidates + loss_count);
                const std::vector<unsigned> missing_recovery = {
                    0, recovery_count / 2, recovery_count - 1
                };
                std::vector<uint8_t> original_present(65, 1);
                std::vector<uint8_t> recovery_present(recovery_count, 1);
                for (size_t i = 0; i < missing_originals.size(); ++i)
                    original_present[missing_originals[i]] = 0;
                for (size_t i = 0; i < missing_recovery.size(); ++i)
                    recovery_present[missing_recovery[i]] = 0;

                ready.fetch_add(1, std::memory_order_release);
                announced_ready = true;
                while (!go.load(std::memory_order_acquire))
                    std::this_thread::yield();

                leo2_codec* automatic = make_codec(context, 65,
                    recovery_count, LEO2_PROFILE_LEGACY_HIGH_V1,
                    LEO2_FIELD_GF8);
                leo2_decode_plan* automatic_plan = NULL;
                require_result(leo2_decode_plan_create(automatic,
                    &original_present[0], &recovery_present[0],
                    &automatic_plan), "concurrent direct plan create");
                size_t automatic_scratch = 0;
                require_result(leo2_decode_plan_scratch_size(
                    automatic_plan, 64, &automatic_scratch),
                    "concurrent direct scratch query");

                leo2_codec_options transform_options;
                memset(&transform_options, 0, sizeof(transform_options));
                transform_options.struct_size = sizeof(transform_options);
                transform_options.flags = LEO2_CODEC_FORCE_GENERIC_DECODE;
                leo2_codec* transform = NULL;
                require_result(leo2_codec_create(context, 65,
                    recovery_count, LEO2_PROFILE_LEGACY_HIGH_V1,
                    LEO2_FIELD_GF8, &transform_options, &transform),
                    "concurrent transform codec create");
                leo2_decode_plan* transform_plan = NULL;
                require_result(leo2_decode_plan_create(transform,
                    &original_present[0], &recovery_present[0],
                    &transform_plan), "concurrent transform plan create");
                size_t transform_scratch = 0;
                require_result(leo2_decode_plan_scratch_size(
                    transform_plan, 64, &transform_scratch),
                    "concurrent transform scratch query");
                require(automatic_scratch < transform_scratch,
                    "concurrent K=65 plan did not select direct repair");
                leo2_decode_plan_destroy(transform_plan);
                leo2_codec_destroy(transform);
                leo2_decode_plan_destroy(automatic_plan);
                leo2_codec_destroy(automatic);

                run_decode_case(context, 65, recovery_count,
                    LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 64,
                    missing_originals, missing_recovery,
                    &lane_counts[lane]);
            }
            catch (const std::exception& error)
            {
                if (!announced_ready)
                    ready.fetch_add(1, std::memory_order_release);
                errors[lane] = error.what();
            }
            catch (...)
            {
                if (!announced_ready)
                    ready.fetch_add(1, std::memory_order_release);
                errors[lane] = "unknown exception before or during lane";
            }
        }));
    }

    while (ready.load(std::memory_order_acquire) != lane_count)
        std::this_thread::yield();
    go.store(true, std::memory_order_release);
    for (size_t lane = 0; lane < threads.size(); ++lane)
        threads[lane].join();
    for (unsigned lane = 0; lane < lane_count; ++lane)
    {
        require(errors[lane].empty(),
            std::string("concurrent direct cache lane failed: ") +
                errors[lane]);
        counts->high_compatibility += lane_counts[lane].high_compatibility;
        counts->low_oracle_symbols += lane_counts[lane].low_oracle_symbols;
        counts->recovered_shards += lane_counts[lane].recovered_shards;
        counts->tail_cases += lane_counts[lane].tail_cases;
        counts->plan_executions += lane_counts[lane].plan_executions;
    }
}

void test_expanded_direct_repair_execution(TestCounts* counts)
{
    leo2_context_options options;
    memset(&options, 0, sizeof(options));
    options.struct_size = sizeof(options);
    options.backend = LEO2_BACKEND_AVX2;
    options.thread_count = 1;
    leo2_context* context = NULL;
    const leo2_result result = leo2_context_create(&options, &context);
    require(result == LEO2_SUCCESS || result == LEO2_UNSUPPORTED,
        "explicit AVX2 direct-repair context returned an unexpected result");
    if (result == LEO2_UNSUPPORTED)
    {
        require(context == NULL,
            "unsupported AVX2 direct-repair context was not cleared");
        return;
    }
    require(context != NULL &&
            leo2_context_backend(context) == LEO2_BACKEND_AVX2,
        "explicit AVX2 direct-repair context selected the wrong backend");

    test_concurrent_expanded_direct_repair_cache(context, counts);
    test_direct_repair_dispatch_bounds(context);
#if LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT
    test_generalized_one_loss_direct_repair_execution(context, counts);
#if !defined(LEO2_GFNI_VARIANT)
    test_generalized_four_tiny_unaligned(context, counts);
#endif
#endif
#if LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE != 0
    test_experimental_small_direct_repair_execution(context, counts);
#endif
    const std::vector<unsigned> missing_originals = {
        0, 1, 7, 16, 32, 48, 63, 64
    };
    const std::vector<unsigned> missing_recovery = { 0, 3, 64 };
    const size_t byte_counts[] = {
        1, 7, 31, 63, 64, 65, 127, 128, 1024, 1025,
        2048, 2049, 4099, 65536
    };
    for (size_t i = 0;
         i < sizeof(byte_counts) / sizeof(byte_counts[0]); ++i)
    {
        run_decode_case(context, 65, 65,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
            byte_counts[i], missing_originals, missing_recovery, counts);
    }
    /* The AVX2 source-major executor composes its established two- and
       four-output kernels for every loss count through eight.  Exercise each
       composition at the exact 2 KiB dispatch boundary and with an unaligned
       vector tail, comparing each result with the independent forced-generic
       decoder inside run_decode_case. */
    static const size_t source_major_bytes[] = { 2048, 2049 };
    static const size_t source_major_losses[] = { 2, 3, 4, 5, 6, 7, 8 };
    for (size_t loss_i = 0;
         loss_i < sizeof(source_major_losses) /
             sizeof(source_major_losses[0]); ++loss_i)
    {
        const std::vector<unsigned> missing(
            missing_originals.begin(),
            missing_originals.begin() + source_major_losses[loss_i]);
        for (size_t byte_i = 0;
             byte_i < sizeof(source_major_bytes) /
                 sizeof(source_major_bytes[0]); ++byte_i)
        {
            run_decode_case(context, 65, 65,
                LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
                source_major_bytes[byte_i], missing, missing_recovery,
                counts);
        }
    }
    test_concurrent_preformatted_source_major_plan(context, counts);
    run_decode_case(context, 65, 127,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        63, missing_originals, std::vector<unsigned>{0, 64, 126}, counts);
    const unsigned recovery_counts[] = { 65, 66, 96, 128 };
    const unsigned loss_counts[] = { 1, 2, 3, 4, 5, 6, 7, 8 };
    for (size_t recovery_i = 0;
         recovery_i < sizeof(recovery_counts) / sizeof(recovery_counts[0]);
         ++recovery_i)
    {
        const unsigned recovery_count = recovery_counts[recovery_i];
        for (size_t loss_i = 0;
             loss_i < sizeof(loss_counts) / sizeof(loss_counts[0]); ++loss_i)
        {
            const unsigned loss_count = loss_counts[loss_i];
            run_decode_case(context, 65, recovery_count,
                LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
                64,
                std::vector<unsigned>(
                    missing_originals.begin(),
                    missing_originals.begin() + loss_count),
                std::vector<unsigned>{
                    0, recovery_count / 2, recovery_count - 1
                },
                counts);
            if (loss_count > 1)
            {
                run_decode_case(context, 65, recovery_count,
                    LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
                    2049,
                    std::vector<unsigned>(
                        missing_originals.begin(),
                        missing_originals.begin() + loss_count),
                    std::vector<unsigned>{
                        0, recovery_count / 2, recovery_count - 1
                    },
                    counts);
            }
        }
    }
    leo2_context_destroy(context);
}

void test_direct_repair_field_helpers()
{
    const BinaryField gf8 = leopard2_test::make_legacy_gf8();
    for (unsigned a = 0; a < 256; ++a)
        for (unsigned b = 0; b < 256; ++b)
            require(leopard::ff8::MultiplyElements(
                        static_cast<uint8_t>(a), static_cast<uint8_t>(b)) ==
                    gf8.multiply(static_cast<Element>(a), static_cast<Element>(b)),
                "GF8 production element multiply differs from direct oracle");
    for (unsigned value = 1; value < 256; ++value)
        require(leopard::ff8::InverseElement(static_cast<uint8_t>(value)) ==
                gf8.inverse(static_cast<Element>(value)),
            "GF8 production inverse differs from direct oracle");

    const uint8_t gf8_coefficient = 173;
    const uint8_t gf8_log = leopard::ff8::ElementLog(gf8_coefficient);
    const size_t gf8_sizes[] = { 1, 63, 64, 65, 257 };
    for (size_t size_i = 0; size_i < sizeof(gf8_sizes) / sizeof(gf8_sizes[0]); ++size_i)
    {
        const size_t bytes = gf8_sizes[size_i];
        std::vector<uint8_t> source(bytes);
        std::vector<uint8_t> product(bytes, 0);
        std::vector<uint8_t> sum(bytes, 0xa5);
        for (size_t i = 0; i < bytes; ++i)
            source[i] = static_cast<uint8_t>(i * 37u + 11u);
        leopard::ff8::MultiplyBytes(&product[0], &source[0], gf8_log, bytes);
        leopard::ff8::MultiplyAddBytes(&sum[0], &source[0], gf8_log, bytes);
        for (size_t i = 0; i < bytes; ++i)
        {
            const uint8_t expected = static_cast<uint8_t>(gf8.multiply(
                source[i], gf8_coefficient));
            require(product[i] == expected && sum[i] == (0xa5 ^ expected),
                "GF8 fixed multiplier helper mishandled a SIMD tail");
        }
    }

    const BinaryField gf16 = leopard2_test::make_legacy_gf16();
    uint32_t state = 0x91e10da5u;
    for (unsigned i = 0; i < 4096; ++i)
    {
        state = state * 1664525u + 1013904223u;
        const uint16_t a = static_cast<uint16_t>(state);
        state = state * 1664525u + 1013904223u;
        const uint16_t b = static_cast<uint16_t>(state);
        require(leopard::ff16::MultiplyElements(a, b) == gf16.multiply(a, b),
            "GF16 production element multiply differs from direct oracle");
        if (a != 0)
            require(leopard::ff16::InverseElement(a) == gf16.inverse(a),
                "GF16 production inverse differs from direct oracle");
    }

    const uint16_t gf16_coefficient = 0x93a7;
    const uint16_t gf16_log = leopard::ff16::ElementLog(gf16_coefficient);
    const size_t gf16_sizes[] = { 2, 34, 64, 66, 1026 };
    for (size_t size_i = 0;
         size_i < sizeof(gf16_sizes) / sizeof(gf16_sizes[0]); ++size_i)
    {
        const size_t bytes = gf16_sizes[size_i];
        const size_t symbols = bytes / 2;
        std::vector<uint16_t> values(symbols);
        std::vector<uint8_t> source(bytes, 0);
        std::vector<uint8_t> product(bytes, 0);
        std::vector<uint8_t> sum(bytes, 0x5a);
        for (size_t symbol = 0; symbol < symbols; ++symbol)
            values[symbol] = static_cast<uint16_t>(symbol * 4051u + 0x1234u);
        for (size_t tile = 0; tile < bytes; tile += 64)
        {
            const size_t tile_symbols = std::min<size_t>(32, (bytes - tile) / 2);
            const size_t first_symbol = tile / 2;
            for (size_t lane = 0; lane < tile_symbols; ++lane)
            {
                const uint16_t value = values[first_symbol + lane];
                source[tile + lane] = static_cast<uint8_t>(value);
                source[tile + tile_symbols + lane] = static_cast<uint8_t>(value >> 8);
            }
        }
        leopard::ff16::MultiplyBytes(&product[0], &source[0], gf16_log, bytes);
        leopard::ff16::MultiplyAddBytes(&sum[0], &source[0], gf16_log, bytes);
        for (size_t tile = 0; tile < bytes; tile += 64)
        {
            const size_t tile_symbols = std::min<size_t>(32, (bytes - tile) / 2);
            const size_t first_symbol = tile / 2;
            for (size_t lane = 0; lane < tile_symbols; ++lane)
            {
                const uint16_t expected = static_cast<uint16_t>(gf16.multiply(
                    values[first_symbol + lane], gf16_coefficient));
                require(product[tile + lane] == static_cast<uint8_t>(expected) &&
                        product[tile + tile_symbols + lane] ==
                            static_cast<uint8_t>(expected >> 8) &&
                        sum[tile + lane] == (0x5a ^ static_cast<uint8_t>(expected)) &&
                        sum[tile + tile_symbols + lane] ==
                            (0x5a ^ static_cast<uint8_t>(expected >> 8)),
                    "GF16 fixed multiplier helper mishandled a compact tail");
            }
        }
    }
}

void test_overlap_rejection(leo2_context* context)
{
    leo2_codec* codec = make_codec(
        context, 3, 2, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    Shards original = make_originals(3, 17, 42);
    std::vector<const void*> original_ptrs = const_pointers(original);
    std::vector<void*> recovery_ptrs(2, NULL);
    recovery_ptrs[0] = &original[0][0];
    Shards output(1, std::vector<uint8_t>(17));
    recovery_ptrs[1] = &output[0][0];
    size_t scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(codec, 17, &scratch_bytes),
        "overlap scratch query");
    AlignedBuffer scratch(scratch_bytes);
    require(leo2_encode(codec, 17, &original_ptrs[0], &recovery_ptrs[0],
        scratch.data, scratch.bytes) == LEO2_OVERLAP,
        "unsupported encode overlap was not rejected");
    leo2_codec_destroy(codec);
}

void test_decode_late_scratch_alias_failure_atomicity(leo2_context* context)
{
    const unsigned k = 8;
    const unsigned r = 4;
    const size_t bytes = 64;
    leo2_codec_options options;
    memset(&options, 0, sizeof(options));
    options.struct_size = sizeof(options);
    options.flags = LEO2_CODEC_FORCE_SPECIALIZED_DECODE;
    leo2_codec* codec = NULL;
    require_result(leo2_codec_create(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, &options, &codec),
        "late-alias codec create");
    const Shards original = make_originals(k, bytes, 0x1a7ea11a5ULL);
    const Shards recovery = encode_new(codec, original, bytes);

    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 1);
    original_present[0] = 0;
    original_present[k - 1] = 0;
    leo2_decode_plan* plan = NULL;
    require_result(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), "late-alias plan create");
    size_t scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(plan, bytes, &scratch_bytes),
        "late-alias scratch query");
    require(scratch_bytes >= bytes,
        "late-alias fixture scratch is smaller than one shard");
    AlignedBuffer scratch(scratch_bytes);
    Shards restored(k, std::vector<uint8_t>(bytes, 0x5a));
    std::vector<const void*> original_ptrs = const_pointers(original);
    std::vector<const void*> recovery_ptrs = const_pointers(recovery);
    std::vector<void*> restored_ptrs(k, NULL);
    original_ptrs[0] = NULL;
    original_ptrs[k - 1] = NULL;
    restored_ptrs[0] = &restored[0][0];
    restored_ptrs[k - 1] = &restored[k - 1][0];

    const std::vector<const void*> valid_original_ptrs = original_ptrs;
    const std::vector<const void*> valid_recovery_ptrs = recovery_ptrs;
    const std::vector<void*> valid_restored_ptrs = restored_ptrs;

    require_result(leo2_decode_plan_execute(plan, bytes, &original_ptrs[0],
        &recovery_ptrs[0], &restored_ptrs[0], scratch.data, scratch.bytes),
        "late-alias valid control execute");
    require(restored[0] == original[0] && restored[k - 1] == original[k - 1],
        "late-alias valid control recovered incorrect data");
    std::fill(restored[0].begin(), restored[0].end(), 0x5a);
    std::fill(restored[k - 1].begin(), restored[k - 1].end(), 0x5a);
    const Shards original_before = original;
    const Shards recovery_before = recovery;
    const Shards restored_before = restored;
    const auto require_atomic_overlap = [&](const char* label) {
        memset(scratch.data, 0xa5, scratch.bytes);
        std::vector<uint8_t> scratch_before(scratch.bytes);
        memcpy(&scratch_before[0], scratch.data, scratch.bytes);
        require(leo2_decode_plan_execute(plan, bytes, &original_ptrs[0],
            &recovery_ptrs[0], &restored_ptrs[0], scratch.data,
            scratch.bytes) == LEO2_OVERLAP, label);
        require(memcmp(scratch.data, &scratch_before[0], scratch.bytes) == 0,
            std::string(label) + " changed scratch before rejecting overlap");
        require(original == original_before && recovery == recovery_before,
            std::string(label) + " changed an input before rejecting overlap");
        require(restored == restored_before,
            std::string(label) + " changed an output before rejecting overlap");
        original_ptrs = valid_original_ptrs;
        recovery_ptrs = valid_recovery_ptrs;
        restored_ptrs = valid_restored_ptrs;
    };

    original_ptrs[k - 2] = scratch.data;
    require_atomic_overlap("late original/scratch overlap was not rejected");
    recovery_ptrs[r - 1] = scratch.data;
    require_atomic_overlap("late recovery/scratch overlap was not rejected");
    restored_ptrs[k - 1] = scratch.data;
    require_atomic_overlap("late output/scratch overlap was not rejected");

    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
}

void test_gf16_byte_granularity(leo2_context* context)
{
    leo2_codec* codec = make_codec(
        context, 257, 33, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16);
    const size_t odd_counts[] = { 1, 3, 33, 63, 65, 1023, 1025 };
    for (size_t count_i = 0;
         count_i < sizeof(odd_counts) / sizeof(odd_counts[0]); ++count_i)
    {
        const size_t bytes = odd_counts[count_i];
        size_t scratch_bytes = 0;
        require(leo2_encode_scratch_size(codec, bytes, &scratch_bytes) ==
            LEO2_UNSUPPORTED, "odd GF16 encode length was not rejected");
        require(leo2_encode(codec, bytes, NULL, NULL, NULL, 0) ==
            LEO2_UNSUPPORTED, "odd GF16 encode execution disagrees with query");
        require(leo2_decode_scratch_size(codec, bytes, &scratch_bytes) ==
            LEO2_UNSUPPORTED, "odd GF16 decode length was not rejected");
    }

    std::vector<uint8_t> original_present(257, 1);
    std::vector<uint8_t> recovery_present(33, 1);
    original_present[0] = 0;
    leo2_decode_plan* loss_plan = NULL;
    require_result(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &loss_plan), "odd GF16 loss plan create");
    size_t scratch_bytes = 0;
    require(leo2_decode_plan_scratch_size(loss_plan, 65, &scratch_bytes) ==
        LEO2_UNSUPPORTED, "odd GF16 loss-plan scratch was not rejected");
    require(leo2_decode_plan_execute(loss_plan, 65, NULL, NULL, NULL, NULL, 0) ==
        LEO2_UNSUPPORTED, "odd GF16 loss-plan execution was not rejected");
    leo2_decode_plan_destroy(loss_plan);

    original_present[0] = 1;
    std::fill(recovery_present.begin(), recovery_present.end(), 0);
    leo2_decode_plan* no_loss_plan = NULL;
    require_result(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &no_loss_plan), "odd GF16 no-loss plan create");
    scratch_bytes = 99;
    require_result(leo2_decode_plan_scratch_size(no_loss_plan, 65, &scratch_bytes),
        "odd GF16 no-loss scratch query");
    require(scratch_bytes == 0, "odd GF16 no-loss plan requires scratch");
    require_result(leo2_decode_plan_execute(
        no_loss_plan, 65, NULL, NULL, NULL, NULL, 0),
        "odd GF16 no-loss execution");
    leo2_decode_plan_destroy(no_loss_plan);
    leo2_codec_destroy(codec);
}

void test_forced_backend(leo2_context* context)
{
    const char* expected = std::getenv("LEO2_EXPECT_BACKEND");
    if (!expected || expected[0] == '\0')
        return;
    leo2_backend backend = LEO2_BACKEND_AUTO;
    if (std::strcmp(expected, "scalar") == 0)
        backend = LEO2_BACKEND_SCALAR;
    else if (std::strcmp(expected, "ssse3") == 0)
        backend = LEO2_BACKEND_SSSE3;
    else if (std::strcmp(expected, "avx2") == 0)
        backend = LEO2_BACKEND_AVX2;
    else if (std::strcmp(expected, "avx512") == 0)
        backend = LEO2_BACKEND_AVX512;
    else if (std::strcmp(expected, "gfni") == 0 ||
             std::strcmp(expected, "avx2-gfni") == 0)
        backend = LEO2_BACKEND_GFNI;
    else
        throw std::runtime_error("invalid LEO2_EXPECT_BACKEND test value");
    require(leo2_context_backend(context) == backend,
        "forced build selected the wrong runtime backend");
}

void test_codec_flag_validation(leo2_context* context)
{
    const int invalid_profiles[] = { -1, 99 };
    for (size_t i = 0;
         i < sizeof(invalid_profiles) / sizeof(invalid_profiles[0]); ++i)
    {
        leo2_codec* invalid =
            reinterpret_cast<leo2_codec*>(static_cast<uintptr_t>(1));
        require(leo2_codec_create(context, 9, 7, invalid_profiles[i],
                LEO2_FIELD_GF8, NULL, &invalid) == LEO2_INVALID_ARGUMENT,
            "invalid integer profile was accepted");
        require(invalid == NULL,
            "invalid integer profile did not clear codec output");
    }
    const int invalid_fields[] = { -1, 99 };
    for (size_t i = 0;
         i < sizeof(invalid_fields) / sizeof(invalid_fields[0]); ++i)
    {
        leo2_codec* invalid =
            reinterpret_cast<leo2_codec*>(static_cast<uintptr_t>(1));
        require(leo2_codec_create(context, 9, 7,
                LEO2_PROFILE_LEGACY_HIGH_V1, invalid_fields[i], NULL,
                &invalid) == LEO2_INVALID_ARGUMENT,
            "invalid integer field was accepted");
        require(invalid == NULL,
            "invalid integer field did not clear codec output");
    }

    leo2_codec_options options;
    memset(&options, 0, sizeof(options));
    options.struct_size = sizeof(options);
    options.flags = LEO2_CODEC_FORCE_GENERIC_DECODE |
        LEO2_CODEC_FORCE_SPECIALIZED_DECODE;
    leo2_codec* codec = reinterpret_cast<leo2_codec*>(static_cast<uintptr_t>(1));
    require(leo2_codec_create(context, 9, 7, LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF8, &options, &codec) == LEO2_INVALID_ARGUMENT,
        "mutually exclusive decoder flags were accepted");
    require(codec == NULL, "failed codec creation did not clear its output");

    options.flags = LEO2_CODEC_FORCE_TILED_DECODE |
        LEO2_CODEC_FORCE_MATERIALIZED_DECODE;
    codec = reinterpret_cast<leo2_codec*>(static_cast<uintptr_t>(1));
    require(leo2_codec_create(context, 9, 7, LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF8, &options, &codec) == LEO2_INVALID_ARGUMENT,
        "mutually exclusive specialized workspace flags were accepted");
    require(codec == NULL, "workspace-flag failure did not clear its output");

    options.flags = LEO2_CODEC_FORCE_GENERIC_DECODE |
        LEO2_CODEC_FORCE_TILED_DECODE;
    codec = reinterpret_cast<leo2_codec*>(static_cast<uintptr_t>(1));
    require(leo2_codec_create(context, 9, 7, LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF8, &options, &codec) == LEO2_INVALID_ARGUMENT,
        "generic decode accepted a specialized workspace flag");
    require(codec == NULL, "generic/workspace failure did not clear its output");

    options.flags = LEO2_CODEC_FORCE_SPECIALIZED_DECODE |
        LEO2_CODEC_FORCE_TILED_DECODE;
    codec = NULL;
    require_result(leo2_codec_create(context, 9, 7,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, &options, &codec),
        "specialized tiled codec create");
    leo2_codec_destroy(codec);

    options.flags = LEO2_CODEC_FORCE_SPECIALIZED_DECODE |
        LEO2_CODEC_FORCE_MATERIALIZED_DECODE;
    codec = NULL;
    require_result(leo2_codec_create(context, 9, 7,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, &options, &codec),
        "specialized materialized codec create");
    leo2_codec_destroy(codec);

    options.flags = 0x80000000u;
    codec = reinterpret_cast<leo2_codec*>(static_cast<uintptr_t>(1));
    require(leo2_codec_create(context, 9, 7, LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF8, &options, &codec) == LEO2_INVALID_ARGUMENT,
        "unknown codec flag was accepted");
    require(codec == NULL, "unknown-flag failure did not clear its output");
}

void test_tiled_high_dispatch_policy()
{
    using leopard2_internal::ShouldUseMaterializedHighDecode;
    const leo2_profile profile = LEO2_PROFILE_LEGACY_HIGH_V1;
    const leo2_field field = LEO2_FIELD_GF8;
    const uint32_t k = 224;
    const uint32_t r = 32;
    const uint32_t padded = 32;
    const uint32_t parent = 256;

    require(ShouldUseMaterializedHighDecode(profile, field, k, r, padded,
        parent, 1, 24 * 1024, LEO2_BACKEND_AVX2),
        "AVX2 materialized policy rejected its lower byte boundary");
    require(ShouldUseMaterializedHighDecode(profile, field, k, r, padded,
        parent, 8, 64 * 1024, LEO2_BACKEND_AVX2),
        "AVX2 materialized policy rejected its upper byte/loss boundary");
    require(ShouldUseMaterializedHighDecode(profile, field, k, r, padded,
        parent, 1, 24 * 1024, LEO2_BACKEND_AVX512),
        "AVX-512 materialized policy rejected its lower byte boundary");
    require(ShouldUseMaterializedHighDecode(profile, field, k, r, padded,
        parent, 1, 32 * 1024, LEO2_BACKEND_SSSE3),
        "SSSE3 materialized policy rejected its lower byte boundary");
    require(ShouldUseMaterializedHighDecode(profile, field, k, r, padded,
        parent, 8, 64 * 1024, LEO2_BACKEND_SSSE3),
        "SSSE3 materialized policy rejected its upper byte/loss boundary");

    require(!ShouldUseMaterializedHighDecode(profile, field, k, r, padded,
        parent, 1, 16 * 1024, LEO2_BACKEND_AVX2),
        "AVX2 materialized policy crossed its lower byte boundary");
    require(!ShouldUseMaterializedHighDecode(profile, field, k, r, padded,
        parent, 1, 24 * 1024, LEO2_BACKEND_SSSE3),
        "SSSE3 materialized policy crossed its lower byte boundary");
    require(!ShouldUseMaterializedHighDecode(profile, field, k, r, padded,
        parent, 1, 64 * 1024 + 64, LEO2_BACKEND_AVX2),
        "materialized policy crossed its upper byte boundary");
    require(!ShouldUseMaterializedHighDecode(profile, field, k, r, padded,
        parent, 9, 32 * 1024, LEO2_BACKEND_AVX2),
        "materialized policy accepted too many missing originals");
    require(!ShouldUseMaterializedHighDecode(profile, field, k, r, padded,
        parent, 1, 32 * 1024, LEO2_BACKEND_SCALAR),
        "materialized policy accepted the scalar backend");
    require(!ShouldUseMaterializedHighDecode(LEO2_PROFILE_LOW_V1, field,
        k, r, padded, parent, 1, 32 * 1024, LEO2_BACKEND_AVX2),
        "materialized policy accepted the low profile");
    require(!ShouldUseMaterializedHighDecode(profile, LEO2_FIELD_GF16,
        k, r, padded, parent, 1, 32 * 1024, LEO2_BACKEND_AVX2),
        "materialized policy accepted GF16");
    require(!ShouldUseMaterializedHighDecode(profile, field, k - 1, r,
        padded, parent, 1, 32 * 1024, LEO2_BACKEND_AVX2),
        "materialized policy accepted a neighboring K");
    require(!ShouldUseMaterializedHighDecode(profile, field, k, r - 1,
        padded, parent, 1, 32 * 1024, LEO2_BACKEND_AVX2),
        "materialized policy accepted a neighboring R");
    require(!ShouldUseMaterializedHighDecode(profile, field, k, r,
        padded / 2, parent, 1, 32 * 1024, LEO2_BACKEND_AVX2),
        "materialized policy accepted a neighboring padded side");
    require(!ShouldUseMaterializedHighDecode(profile, field, k, r, padded,
        parent / 2, 1, 32 * 1024, LEO2_BACKEND_AVX2),
        "materialized policy accepted a neighboring parent");
}

leo2_codec* make_flagged_codec(
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
    require_result(leo2_codec_create(
        context, k, r, profile, field, &options, &codec),
        "flagged codec create");
    require(codec != NULL, "flagged codec create returned a null codec");
    return codec;
}

void run_low_reveal_fusion_case(
    leo2_context* context,
    leo2_field field,
    uint32_t codec_flags,
    size_t bytes,
    unsigned missing_count)
{
    // K is the complete padded side, so the dense case covers the full final
    // transform while the one-loss case covers the C1-pruned final-output
    // schedule.  This is an Algorithm-4 attribution test, so its otherwise
    // AUTO control explicitly stays on the specialized transform path now
    // that redundancy-dominant one-loss AUTO uses direct repair on AVX2.
    const unsigned k = 32;
    const unsigned r = 33;
    uint32_t effective_codec_flags = codec_flags;
    if (effective_codec_flags == 0)
        effective_codec_flags = LEO2_CODEC_FORCE_SPECIALIZED_DECODE;
    leo2_codec* codec = make_flagged_codec(context, k, r,
        LEO2_PROFILE_LOW_V1, field, effective_codec_flags);
    const Shards source = make_originals(k, bytes,
        UINT64_C(0x89ef342100000000) + bytes + field + missing_count);
    const Shards parity = encode_new(codec, source, bytes);
    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 1);
    std::vector<const void*> original_input = const_pointers(source);
    const std::vector<const void*> recovery_input = const_pointers(parity);
    for (unsigned i = 0; i < missing_count; ++i)
    {
        original_present[i] = 0;
        original_input[i] = NULL;
    }

    leo2_decode_plan* plan = NULL;
    require_result(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), "low reveal plan create");
    size_t scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(plan, bytes, &scratch_bytes),
        "low reveal scratch query");
    AlignedBuffer scratch(scratch_bytes);
    Shards restored(k, std::vector<uint8_t>(bytes + 1, 0));
    std::vector<void*> output(k, NULL);
    for (unsigned i = 0; i < missing_count; ++i)
        output[i] = &restored[i][1]; // Exercise unaligned caller outputs.

    leo2_test_reset_low_reveal_counts();
    require_result(leo2_decode_plan_execute(plan, bytes, &original_input[0],
        &recovery_input[0], &output[0], scratch.data, scratch.bytes),
        "low reveal decode execute");
    for (unsigned i = 0; i < missing_count; ++i)
        require(memcmp(&restored[i][1], &source[i][0], bytes) == 0,
            "low reveal decode output mismatch");

    const size_t aligned_prefix = bytes & ~static_cast<size_t>(63u);
    const uint64_t expected_direct =
        aligned_prefix != 0 ? missing_count : 0;
    const uint64_t expected_scratch =
        (bytes & 63u) != 0 ? missing_count : 0;
    require(leo2_test_low_direct_reveal_shards() == expected_direct,
        "low direct reveal counter disagrees with complete-tile policy");
    require(leo2_test_low_scratch_reveal_shards() == expected_scratch,
        "low scratch reveal counter disagrees with tail policy");

    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
}

void test_low_reveal_fusion(leo2_context* context)
{
    const uint32_t codec_flags[] = {
        0,
        LEO2_CODEC_FORCE_SPECIALIZED_DECODE |
            LEO2_CODEC_FORCE_TILED_DECODE,
        LEO2_CODEC_FORCE_SPECIALIZED_DECODE |
            LEO2_CODEC_FORCE_MATERIALIZED_DECODE
    };
    for (size_t mode = 0;
         mode < sizeof(codec_flags) / sizeof(codec_flags[0]); ++mode)
    {
        const unsigned losses[] = { 1, 32 };
        for (size_t loss = 0; loss < sizeof(losses) / sizeof(losses[0]); ++loss)
        {
            run_low_reveal_fusion_case(context, LEO2_FIELD_GF8,
                codec_flags[mode], 17, losses[loss]);
            run_low_reveal_fusion_case(context, LEO2_FIELD_GF8,
                codec_flags[mode], 64, losses[loss]);
            run_low_reveal_fusion_case(context, LEO2_FIELD_GF8,
                codec_flags[mode], 65, losses[loss]);
#ifdef LEO_HAS_FF16
            run_low_reveal_fusion_case(context, LEO2_FIELD_GF16,
                codec_flags[mode], 34, losses[loss]);
            run_low_reveal_fusion_case(context, LEO2_FIELD_GF16,
                codec_flags[mode], 64, losses[loss]);
            run_low_reveal_fusion_case(context, LEO2_FIELD_GF16,
                codec_flags[mode], 66, losses[loss]);
#endif
        }
    }
}

void run_high_reveal_fusion_case(
    leo2_context* context,
    leo2_field field,
    size_t bytes,
    const std::vector<unsigned>& missing,
    uint32_t execution_flag)
{
    const unsigned k = 64;
    const unsigned r = 16;
    require(!missing.empty() && missing.size() <= r,
        "invalid high reveal loss pattern");
    leo2_codec* codec = make_flagged_codec(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, field,
        LEO2_CODEC_FORCE_SPECIALIZED_DECODE |
            execution_flag);
    const Shards source = make_originals(k, bytes,
        UINT64_C(0x62af134900000000) + bytes + field + missing.size());
    const Shards parity = encode_new(codec, source, bytes);
    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 1);
    std::vector<const void*> original_input = const_pointers(source);
    const std::vector<const void*> recovery_input = const_pointers(parity);
    for (size_t i = 0; i < missing.size(); ++i)
    {
        require(missing[i] < k &&
                (i == 0 || missing[i - 1] < missing[i]),
            "high reveal loss pattern is not sorted and unique");
        original_present[missing[i]] = 0;
        original_input[missing[i]] = NULL;
    }

    leo2_decode_plan* plan = NULL;
    require_result(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), "high reveal plan create");
    leopard2_internal::DecodePathInfo path_info;
    require_result(leopard2_internal::GetDecodePlanPathInfo(
        plan, bytes, false, &path_info), "high reveal path introspection");
    require(path_info.path == leopard2_internal::kDecodePathTiled ||
            path_info.path == leopard2_internal::kDecodePathMaterialized,
        "high reveal test did not select a specialized workspace");
    size_t scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(plan, bytes, &scratch_bytes),
        "high reveal scratch query");
    AlignedBuffer scratch(scratch_bytes);
    Shards restored(k, std::vector<uint8_t>(bytes + 2, 0xa5));
    std::vector<void*> output(k, NULL);
    for (size_t i = 0; i < missing.size(); ++i)
        output[missing[i]] = &restored[missing[i]][1];

    leo2_test_reset_high_reveal_counts();
    require_result(leo2_decode_plan_execute(plan, bytes, &original_input[0],
        &recovery_input[0], &output[0], scratch.data, scratch.bytes),
        "high reveal decode execute");
    for (size_t i = 0; i < missing.size(); ++i)
    {
        const unsigned original_index = missing[i];
        require(memcmp(&restored[original_index][1],
                    &source[original_index][0], bytes) == 0,
            "high reveal decode output mismatch");
        require(restored[original_index][0] == 0xa5 &&
                restored[original_index][bytes + 1] == 0xa5,
            "high reveal decode wrote outside its output shard");
    }
    for (unsigned i = 0; i < k; ++i)
    {
        if (!original_present[i])
            continue;
        require(std::find_if(restored[i].begin(), restored[i].end(),
                    [](uint8_t value) { return value != 0xa5; }) ==
                restored[i].end(),
            "high reveal decode touched an unrequested output");
    }

    const size_t aligned_prefix = bytes & ~static_cast<size_t>(63u);
    const uint64_t expected_direct =
        aligned_prefix != 0 ? missing.size() : 0;
    const uint64_t expected_scratch =
        (bytes & 63u) != 0 ? missing.size() : 0;
    require(leo2_test_high_direct_reveal_shards() == expected_direct,
        "high direct reveal counter disagrees with complete-tile policy");
    const uint64_t expected_materialized_direct =
        path_info.path == leopard2_internal::kDecodePathMaterialized &&
            aligned_prefix != 0
        ? missing.size() : 0;
    require(leo2_test_high_materialized_direct_reveal_shards() ==
            expected_materialized_direct,
        "materialized high reveal counter disagrees with selected path");
    require(leo2_test_high_scratch_reveal_shards() == expected_scratch,
        "high scratch reveal counter disagrees with tail policy");

    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
}

void test_high_reveal_fusion(leo2_context* context)
{
    const uint32_t execution_flags[] = {
        0,
        LEO2_CODEC_FORCE_TILED_DECODE,
        LEO2_CODEC_FORCE_MATERIALIZED_DECODE
    };
    const unsigned losses[] = { 1, 16 };
    for (size_t mode = 0;
         mode < sizeof(execution_flags) / sizeof(execution_flags[0]); ++mode)
    {
        for (size_t loss = 0;
             loss < sizeof(losses) / sizeof(losses[0]); ++loss)
        {
            std::vector<unsigned> missing(losses[loss]);
            for (unsigned i = 0; i < losses[loss]; ++i)
                missing[i] = i;
            run_high_reveal_fusion_case(
                context, LEO2_FIELD_GF8, 17, missing,
                execution_flags[mode]);
            run_high_reveal_fusion_case(
                context, LEO2_FIELD_GF8, 64, missing,
                execution_flags[mode]);
            run_high_reveal_fusion_case(
                context, LEO2_FIELD_GF8, 65, missing,
                execution_flags[mode]);
#ifdef LEO_HAS_FF16
            run_high_reveal_fusion_case(
                context, LEO2_FIELD_GF16, 34, missing,
                execution_flags[mode]);
            run_high_reveal_fusion_case(
                context, LEO2_FIELD_GF16, 64, missing,
                execution_flags[mode]);
            run_high_reveal_fusion_case(
                context, LEO2_FIELD_GF16, 66, missing,
                execution_flags[mode]);
#endif
        }
    }

    /* Exercise requested-output indexing across every T-sized message block,
       including aligned-plus-tail restoration of the scratch destinations. */
    const unsigned scattered_indices[] = { 0, 15, 16, 31, 32, 47, 48, 63 };
    const std::vector<unsigned> scattered(scattered_indices,
        scattered_indices + sizeof(scattered_indices) /
            sizeof(scattered_indices[0]));
    for (size_t mode = 0;
         mode < sizeof(execution_flags) / sizeof(execution_flags[0]); ++mode)
    {
        run_high_reveal_fusion_case(
            context, LEO2_FIELD_GF8, 65, scattered,
            execution_flags[mode]);
#ifdef LEO_HAS_FF16
        run_high_reveal_fusion_case(
            context, LEO2_FIELD_GF16, 66, scattered,
            execution_flags[mode]);
#endif
    }
}

void test_tiled_materialized_execution(leo2_context* context)
{
    const unsigned k = 224;
    const unsigned r = 32;
    const size_t bytes = 32 * 1024;
    const unsigned missing_count = 8;
    leo2_codec* auto_codec = make_codec(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    leo2_codec* tiled_codec = make_flagged_codec(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        LEO2_CODEC_FORCE_SPECIALIZED_DECODE |
            LEO2_CODEC_FORCE_TILED_DECODE);
    leo2_codec* materialized_codec = make_flagged_codec(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        LEO2_CODEC_FORCE_SPECIALIZED_DECODE |
            LEO2_CODEC_FORCE_MATERIALIZED_DECODE);

    const Shards source = make_originals(k, bytes, UINT64_C(0x71e2d15a));
    const Shards parity = encode_new(auto_codec, source, bytes);
    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 1);
    for (unsigned i = 0; i < missing_count; ++i)
        original_present[i] = 0;

    leo2_decode_plan* plans[3] = { NULL, NULL, NULL };
    leo2_codec* codecs[3] = { auto_codec, tiled_codec, materialized_codec };
    size_t scratch_bytes[3] = { 0, 0, 0 };
    for (unsigned i = 0; i < 3; ++i)
    {
        require_result(leo2_decode_plan_create(codecs[i], &original_present[0],
            &recovery_present[0], &plans[i]), "workspace plan create");
        require_result(leo2_decode_plan_scratch_size(
            plans[i], bytes, &scratch_bytes[i]), "workspace scratch query");
    }
    require(scratch_bytes[2] > scratch_bytes[1],
        "materialized workspace is not larger than tiled workspace");
    const leo2_backend backend = leo2_context_backend(context);
    if (backend == LEO2_BACKEND_AVX2 ||
        backend == LEO2_BACKEND_AVX512 ||
        backend == LEO2_BACKEND_SSSE3)
        require(scratch_bytes[0] == scratch_bytes[2],
            "AUTO did not reserve its calibrated materialized workspace");
    else
        require(scratch_bytes[0] == scratch_bytes[1],
            "uncalibrated backend did not retain the tiled workspace");
    leopard2_internal::DecodePathInfo path_info[3];
    for (unsigned i = 0; i < 3; ++i)
        require_result(leopard2_internal::GetDecodePlanPathInfo(
            plans[i], bytes, false, &path_info[i]),
            "workspace path introspection");
    const leopard2_internal::DecodePath expected_auto =
        backend == LEO2_BACKEND_AVX2 ||
        backend == LEO2_BACKEND_AVX512 ||
        backend == LEO2_BACKEND_SSSE3
            ? leopard2_internal::kDecodePathMaterialized
            : leopard2_internal::kDecodePathTiled;
    require(path_info[0].path == expected_auto &&
            path_info[1].path == leopard2_internal::kDecodePathTiled &&
            path_info[2].path == leopard2_internal::kDecodePathMaterialized,
        "workspace introspection differs from selected single-item paths");
    require(path_info[0].required_work_slots ==
                (expected_auto == leopard2_internal::kDecodePathMaterialized
                    ? 256u : 72u) &&
            path_info[1].required_work_slots == 72 &&
            path_info[2].required_work_slots == 256,
        "workspace introspection slot accounting differs");

    leopard2_internal::DecodePathInfo batch_path;
    require_result(leopard2_internal::GetDecodePlanPathInfo(
        plans[0], bytes, true, &batch_path),
        "AUTO batch path introspection");
    const leopard2_internal::DecodePath expected_batch =
        backend == LEO2_BACKEND_SSSE3
            ? leopard2_internal::kDecodePathMaterialized
            : leopard2_internal::kDecodePathTiled;
    require(batch_path.path == expected_batch &&
            batch_path.required_work_slots ==
                (expected_batch == leopard2_internal::kDecodePathMaterialized
                    ? 256u : 72u),
        "AUTO batch introspection differs from execution workspace");

    std::vector<const void*> original_inputs = const_pointers(source);
    std::vector<const void*> recovery_inputs = const_pointers(parity);
    for (unsigned i = 0; i < missing_count; ++i)
        original_inputs[i] = NULL;
    Shards restored[3] = {
        Shards(k, std::vector<uint8_t>(bytes, 0)),
        Shards(k, std::vector<uint8_t>(bytes, 0)),
        Shards(k, std::vector<uint8_t>(bytes, 0))
    };
    for (unsigned path = 0; path < 3; ++path)
    {
        std::vector<void*> output(k, NULL);
        for (unsigned i = 0; i < missing_count; ++i)
            output[i] = &restored[path][i][0];
        AlignedBuffer scratch(scratch_bytes[path]);
        leo2_test_reset_high_reveal_counts();
        require_result(leo2_decode_plan_execute(plans[path], bytes,
            &original_inputs[0], &recovery_inputs[0], &output[0],
            scratch.data, scratch.bytes), "workspace decode execute");
        const uint64_t expected_direct = missing_count;
        const uint64_t expected_materialized_direct =
            path_info[path].path == leopard2_internal::kDecodePathMaterialized
                ? missing_count : 0;
        require(leo2_test_high_direct_reveal_shards() == expected_direct &&
                leo2_test_high_materialized_direct_reveal_shards() ==
                    expected_materialized_direct &&
                leo2_test_high_scratch_reveal_shards() == 0,
            "high reveal fusion escaped the selected specialized path");
        for (unsigned i = 0; i < missing_count; ++i)
            require(restored[path][i] == source[i],
                "workspace decode restored the wrong original");
    }
    for (unsigned i = 0; i < missing_count; ++i)
    {
        require(restored[0][i] == restored[1][i],
            "AUTO and tiled decode differ");
        require(restored[0][i] == restored[2][i],
            "AUTO and materialized decode differ");
    }

    Shards batch_restored[2] = {
        Shards(k, std::vector<uint8_t>(bytes, 0)),
        Shards(k, std::vector<uint8_t>(bytes, 0))
    };
    std::vector<void*> batch_output[2] = {
        std::vector<void*>(k, NULL), std::vector<void*>(k, NULL)
    };
    AlignedBuffer batch_scratch0(scratch_bytes[0]);
    AlignedBuffer batch_scratch1(scratch_bytes[0]);
    leo2_decode_batch_item items[2];
    memset(items, 0, sizeof(items));
    for (unsigned item = 0; item < 2; ++item)
    {
        for (unsigned i = 0; i < missing_count; ++i)
            batch_output[item][i] = &batch_restored[item][i][0];
        items[item].shard_bytes = bytes;
        items[item].original = &original_inputs[0];
        items[item].recovery = &recovery_inputs[0];
        items[item].restored_original = &batch_output[item][0];
        items[item].scratch = item == 0
            ? batch_scratch0.data : batch_scratch1.data;
        items[item].scratch_bytes = scratch_bytes[0];
    }
    leo2_test_reset_high_reveal_counts();
    require_result(leo2_decode_plan_execute_batch(plans[0], items, 2),
        "AUTO tiled batch execute");
    const uint64_t expected_batch_direct = 2u * missing_count;
    const uint64_t expected_batch_materialized_direct =
        batch_path.path == leopard2_internal::kDecodePathMaterialized
            ? expected_batch_direct : 0;
    require(leo2_test_high_direct_reveal_shards() == expected_batch_direct &&
            leo2_test_high_materialized_direct_reveal_shards() ==
                expected_batch_materialized_direct &&
            leo2_test_high_scratch_reveal_shards() == 0,
        "high reveal fusion escaped the selected AUTO batch path");
    for (unsigned item = 0; item < 2; ++item)
        for (unsigned i = 0; i < missing_count; ++i)
            require(batch_restored[item][i] == source[i],
                "AUTO tiled batch restored the wrong original");

    for (unsigned i = 0; i < 3; ++i)
        leo2_decode_plan_destroy(plans[i]);
    leo2_codec_destroy(materialized_codec);
    leo2_codec_destroy(tiled_codec);
    leo2_codec_destroy(auto_codec);
}

void test_available_avx512_batch_dispatch()
{
    leo2_context_options options;
    memset(&options, 0, sizeof(options));
    options.struct_size = sizeof(options);
    options.backend = LEO2_BACKEND_AVX512;
    options.thread_count = 1;
    leo2_context* context = NULL;
    const leo2_result result = leo2_context_create(&options, &context);
    if (result == LEO2_UNSUPPORTED)
    {
        require(context == NULL,
            "unsupported AVX-512 context did not clear its output");
        return;
    }
    require_result(result, "AVX-512 batch-dispatch context create");
    require(leo2_context_backend(context) == LEO2_BACKEND_AVX512,
        "explicit AVX-512 context selected a different backend");
    test_tiled_materialized_execution(context);
    leo2_context_destroy(context);
}

void test_batch_materialized_capacity(leo2_context* context)
{
    const unsigned k = 128;
    const unsigned r = 128;
    const size_t bytes = 257;
    leo2_codec* codec = make_flagged_codec(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        LEO2_CODEC_FORCE_SPECIALIZED_DECODE);
    leo2_codec* tiled_codec = make_flagged_codec(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        LEO2_CODEC_FORCE_SPECIALIZED_DECODE |
            LEO2_CODEC_FORCE_TILED_DECODE);
    const Shards source = make_originals(k, bytes, UINT64_C(0xb47c4d15));
    const Shards parity = encode_new(codec, source, bytes);
    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 1);
    original_present[0] = 0;
    leo2_decode_plan* plan = NULL;
    require_result(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), "materialized-capacity plan create");
    size_t scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(plan, bytes, &scratch_bytes),
        "materialized-capacity scratch query");
    leo2_decode_plan* tiled_plan = NULL;
    require_result(leo2_decode_plan_create(tiled_codec, &original_present[0],
        &recovery_present[0], &tiled_plan), "forced-tiled-capacity plan create");
    size_t tiled_scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(
        tiled_plan, bytes, &tiled_scratch_bytes),
        "forced-tiled-capacity scratch query");
    require(tiled_scratch_bytes == scratch_bytes,
        "translated tiled/materialized N-slot queries differ");

    std::vector<const void*> original_inputs = const_pointers(source);
    std::vector<const void*> recovery_inputs = const_pointers(parity);
    original_inputs[0] = NULL;
    Shards restored[2] = {
        Shards(k, std::vector<uint8_t>(bytes, 0)),
        Shards(k, std::vector<uint8_t>(bytes, 0))
    };
    std::vector<void*> outputs[2] = {
        std::vector<void*>(k, NULL), std::vector<void*>(k, NULL)
    };
    outputs[0][0] = &restored[0][0][0];
    outputs[1][0] = &restored[1][0][0];
    AlignedBuffer scratch0(scratch_bytes);
    AlignedBuffer scratch1(scratch_bytes);
    leo2_decode_batch_item items[2];
    memset(items, 0, sizeof(items));
    for (unsigned item = 0; item < 2; ++item)
    {
        items[item].shard_bytes = bytes;
        items[item].original = &original_inputs[0];
        items[item].recovery = &recovery_inputs[0];
        items[item].restored_original = &outputs[item][0];
        items[item].scratch = item == 0 ? scratch0.data : scratch1.data;
        items[item].scratch_bytes = scratch_bytes;
    }
    require_result(leo2_decode_plan_execute_batch(plan, items, 2),
        "materialized-capacity batch execute");
    require(restored[0][0] == source[0] && restored[1][0] == source[0],
        "materialized-capacity batch restored the wrong original");

    Shards tiled_restored(k, std::vector<uint8_t>(bytes, 0));
    std::vector<void*> tiled_output(k, NULL);
    tiled_output[0] = &tiled_restored[0][0];
    AlignedBuffer tiled_scratch(tiled_scratch_bytes);
    require_result(leo2_decode_plan_execute(tiled_plan, bytes,
        &original_inputs[0], &recovery_inputs[0], &tiled_output[0],
        tiled_scratch.data, tiled_scratch.bytes),
        "forced-tiled-capacity execute");
    require(tiled_restored[0] == source[0],
        "forced tiled capacity restored the wrong original");

    leo2_decode_plan_destroy(tiled_plan);
    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(tiled_codec);
    leo2_codec_destroy(codec);
}

void test_balanced_dispatch_policy()
{
    using leopard2_internal::ShouldUseBalancedGenericDecode;
    const leo2_profile profile = LEO2_PROFILE_LEGACY_HIGH_V1;
    const leo2_field field = LEO2_FIELD_GF8;
    const uint32_t k = 128;
    const uint32_t r = 128;
    const uint32_t padded = 128;
    const uint32_t parent = 256;
    const uint32_t missing = 128;

    require(ShouldUseBalancedGenericDecode(profile, field, k, r, padded,
        parent, missing, 256, LEO2_BACKEND_SCALAR),
        "balanced dispatch rejected its lower byte boundary");
    require(ShouldUseBalancedGenericDecode(profile, field, k, r, padded,
        parent, missing, 1024 * 1024, LEO2_BACKEND_SSSE3),
        "balanced dispatch rejected its upper byte boundary");
    require(ShouldUseBalancedGenericDecode(profile, field, k, r, padded,
        parent, missing, 4096, LEO2_BACKEND_AVX2),
        "balanced dispatch rejected AVX2");
    require(ShouldUseBalancedGenericDecode(profile, field, k, r, padded,
        parent, missing, 4096, LEO2_BACKEND_AVX512),
        "balanced dispatch rejected AVX-512");

    require(!ShouldUseBalancedGenericDecode(profile, field, k, r, padded,
        parent, missing, 255, LEO2_BACKEND_SCALAR),
        "balanced dispatch crossed its lower byte boundary");
    require(!ShouldUseBalancedGenericDecode(profile, field, k, r, padded,
        parent, missing, 1024 * 1024 + 64, LEO2_BACKEND_SCALAR),
        "balanced dispatch crossed its upper byte boundary");
    require(!ShouldUseBalancedGenericDecode(LEO2_PROFILE_LOW_V1, field, k, r,
        padded, parent, missing, 4096, LEO2_BACKEND_SCALAR),
        "balanced dispatch accepted the low profile");
    require(!ShouldUseBalancedGenericDecode(profile, LEO2_FIELD_GF16, k, r,
        padded, parent, missing, 4096, LEO2_BACKEND_SCALAR),
        "balanced dispatch accepted GF16");
    require(!ShouldUseBalancedGenericDecode(profile, field, k - 1, r, padded,
        parent, missing, 4096, LEO2_BACKEND_SCALAR),
        "balanced dispatch accepted a neighboring K");
    require(!ShouldUseBalancedGenericDecode(profile, field, k, r - 1, padded,
        parent, missing, 4096, LEO2_BACKEND_SCALAR),
        "balanced dispatch accepted a neighboring R");
    require(!ShouldUseBalancedGenericDecode(profile, field, k, r, padded / 2,
        parent, missing, 4096, LEO2_BACKEND_SCALAR),
        "balanced dispatch accepted a different padded side");
    require(!ShouldUseBalancedGenericDecode(profile, field, k, r, padded,
        parent / 2, missing, 4096, LEO2_BACKEND_SCALAR),
        "balanced dispatch accepted a different parent");
    require(!ShouldUseBalancedGenericDecode(profile, field, k, r, padded,
        parent, missing - 1, 4096, LEO2_BACKEND_SCALAR),
        "balanced dispatch accepted partial recovery");
    require(!ShouldUseBalancedGenericDecode(profile, field, k, r, padded,
        parent, missing, 4096, LEO2_BACKEND_NEON),
        "balanced dispatch accepted an unmeasured backend");
}

uint32_t test_ceil_pow2(uint32_t value)
{
    uint32_t result = 1;
    while (result < value)
        result <<= 1;
    return result;
}

leopard2_internal::DecodePathInput make_decode_path_input(
    leo2_profile profile,
    leo2_field field,
    leo2_backend backend,
    uint32_t k,
    uint32_t r,
    uint32_t padded,
    uint32_t parent,
    uint32_t missing,
    uint64_t bytes,
    uint32_t flags = 0,
    bool multi_item_batch = false,
    bool plan_known = true)
{
    leopard2_internal::DecodePathInput input;
    input.profile = profile;
    input.field = field;
    input.backend = backend;
    input.original_count = k;
    input.recovery_count = r;
    input.padded_side = padded;
    input.parent_count = parent;
    input.missing_original_count = plan_known ? missing : 0;
    input.requested_output_count = profile == LEO2_PROFILE_LEGACY_HIGH_V1
        ? (plan_known ? missing : r) : 0;
    input.codec_flags = flags;
    input.actual_shard_bytes = bytes;
    input.tail_bytes = static_cast<size_t>(bytes) & 63u;
    input.aligned_prefix_bytes =
        static_cast<size_t>(bytes) - input.tail_bytes;
    input.rounded_shard_bytes = input.aligned_prefix_bytes +
        (input.tail_bytes == 0 ? 0 : 64);
    input.plan_known = plan_known;
    input.multi_item_batch = multi_item_batch;
    return input;
}

leopard2_internal::DecodePathInfo require_decode_path(
    const leopard2_internal::DecodePathInput& input,
    leopard2_internal::DecodePath expected,
    leopard2_internal::DecodePathRule expected_rule,
    const char* operation)
{
    leopard2_internal::DecodePathInfo selection;
    require(leopard2_internal::SelectDecodePath(input, selection),
        std::string(operation) + " rejected a valid selector input");
    require(selection.path == expected,
        std::string(operation) + " selected " +
        leopard2_internal::DecodePathName(selection.path));
    require(selection.rule == expected_rule,
        std::string(operation) + " selected rule " +
        leopard2_internal::DecodePathRuleName(selection.rule));
    require(selection.direct_executor ==
            leopard2_internal::kDirectRepairExecutorNone,
        std::string(operation) +
        " transform selector reported a direct-repair executor");
    require(selection.aligned_prefix_bytes == input.aligned_prefix_bytes &&
            selection.tail_bytes == input.tail_bytes &&
            selection.rounded_shard_bytes == input.rounded_shard_bytes &&
            selection.multi_item_batch == input.multi_item_batch,
        std::string(operation) + " lost byte or batch state");
    return selection;
}

void require_batch_decode_scratch_bound(
    leopard2_internal::DecodePathInput input,
    const char* operation)
{
    leopard2_internal::DecodePathInfo single;
    leopard2_internal::DecodePathInfo batch;
    input.multi_item_batch = false;
    require(leopard2_internal::SelectDecodePath(input, single),
        std::string(operation) + " rejected the single-item input");
    input.multi_item_batch = true;
    require(leopard2_internal::SelectDecodePath(input, batch),
        std::string(operation) + " rejected the multi-item input");
    require(batch.required_work_slots <= single.required_work_slots,
        std::string(operation) +
        " made batch execution exceed one-item scratch");
}

void test_unified_decode_path_selector()
{
    using namespace leopard2_internal;
    DecodePathInput balanced = make_decode_path_input(
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        LEO2_BACKEND_AVX2, 128, 128, 128, 256, 128, 193);
    DecodePathInfo selection = require_decode_path(balanced,
        kDecodePathGeneric, kDecodeRuleBalancedGeneric,
        "balanced ragged lower boundary");
    require(selection.required_work_slots == 256 &&
            selection.matching_auto_rules == kDecodeAutoRuleBalancedGeneric,
        "balanced selector accounting differs");

    balanced.actual_shard_bytes = 192;
    balanced.aligned_prefix_bytes = 192;
    balanced.tail_bytes = 0;
    balanced.rounded_shard_bytes = 192;
    require_decode_path(balanced, kDecodePathMaterialized,
        kDecodeRuleWorkspaceMaterialized,
        "balanced immediate lower neighbor");
    balanced = make_decode_path_input(LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF8, LEO2_BACKEND_SCALAR,
        128, 128, 128, 256, 128, 1024 * 1024);
    require_decode_path(balanced, kDecodePathGeneric,
        kDecodeRuleBalancedGeneric, "balanced upper boundary");
    balanced = make_decode_path_input(LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF8, LEO2_BACKEND_SCALAR,
        128, 128, 128, 256, 128, 1024 * 1024 + 1);
    require_decode_path(balanced, kDecodePathMaterialized,
        kDecodeRuleWorkspaceMaterialized,
        "balanced immediate upper neighbor");

    balanced = make_decode_path_input(LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF8, LEO2_BACKEND_AVX2,
        128, 128, 128, 256, 128, 4097,
        LEO2_CODEC_FORCE_SPECIALIZED_DECODE);
    require_decode_path(balanced, kDecodePathMaterialized,
        kDecodeRuleWorkspaceMaterialized,
        "force-specialized precedence");
    balanced.codec_flags = LEO2_CODEC_FORCE_GENERIC_DECODE;
    require_decode_path(balanced, kDecodePathGeneric,
        kDecodeRuleForcedGeneric, "force-generic precedence");
    balanced.codec_flags = LEO2_CODEC_FORCE_MATERIALIZED_DECODE;
    require_decode_path(balanced, kDecodePathMaterialized,
        kDecodeRuleForcedMaterialized, "force-materialized precedence");
    balanced.codec_flags = LEO2_CODEC_FORCE_TILED_DECODE;
    selection = require_decode_path(balanced, kDecodePathTiled,
        kDecodeRuleForcedTiled, "force-tiled precedence");
    require(selection.required_work_slots == 384,
        "forced balanced tiled work slots differ");

    DecodePathInput high = make_decode_path_input(
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        LEO2_BACKEND_AVX2, 224, 32, 32, 256, 1,
        24 * 1024 - 64);
    require_decode_path(high, kDecodePathTiled, kDecodeRuleWorkspaceTiled,
        "AVX2 materialized immediate lower neighbor");
    high = make_decode_path_input(LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF8, LEO2_BACKEND_AVX2,
        224, 32, 32, 256, 1, 24 * 1024 - 63);
    require_decode_path(high, kDecodePathMaterialized,
        kDecodeRuleMeasuredMaterialized,
        "AVX2 materialized ragged lower boundary");
    high = make_decode_path_input(LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF8, LEO2_BACKEND_SSSE3,
        224, 32, 32, 256, 1, 32 * 1024 - 64);
    require_decode_path(high, kDecodePathTiled, kDecodeRuleWorkspaceTiled,
        "SSSE3 materialized immediate lower neighbor");
    high = make_decode_path_input(LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF8, LEO2_BACKEND_SSSE3,
        224, 32, 32, 256, 1, 32 * 1024 - 63);
    require_decode_path(high, kDecodePathMaterialized,
        kDecodeRuleMeasuredMaterialized,
        "SSSE3 materialized ragged lower boundary");
    high = make_decode_path_input(LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF8, LEO2_BACKEND_AVX2,
        224, 32, 32, 256, 8, 64 * 1024);
    require_decode_path(high, kDecodePathMaterialized,
        kDecodeRuleMeasuredMaterialized,
        "materialized byte/loss upper boundary");
    high = make_decode_path_input(LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF8, LEO2_BACKEND_AVX2,
        224, 32, 32, 256, 8, 64 * 1024 + 1);
    require_decode_path(high, kDecodePathTiled, kDecodeRuleWorkspaceTiled,
        "materialized byte upper neighbor");
    high = make_decode_path_input(LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF8, LEO2_BACKEND_AVX2,
        224, 32, 32, 256, 9, 32 * 1024);
    require_decode_path(high, kDecodePathTiled, kDecodeRuleWorkspaceTiled,
        "materialized loss upper neighbor");

    high = make_decode_path_input(LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF8, LEO2_BACKEND_AVX2,
        224, 32, 32, 256, 8, 32 * 1024, 0, true);
    selection = require_decode_path(high, kDecodePathTiled,
        kDecodeRuleMeasuredBatchTiled, "AVX2 multi-item exception");
    require(selection.required_work_slots == 72,
        "batch tiled selector did not account for L outputs");
    high.backend = LEO2_BACKEND_AVX512;
    selection = require_decode_path(high, kDecodePathTiled,
        kDecodeRuleMeasuredBatchTiled, "AVX-512 multi-item exception");
    require(selection.required_work_slots == 72 &&
            selection.matching_auto_rules ==
                kDecodeAutoRuleMeasuredMaterialized,
        "AVX-512 batch exception lost workspace or measured-rule state");
    high.backend = LEO2_BACKEND_SSSE3;
    require_decode_path(high, kDecodePathMaterialized,
        kDecodeRuleMeasuredMaterialized,
        "SSSE3 multi-item materialized control");
    high.backend = LEO2_BACKEND_AVX2;
    high.codec_flags = LEO2_CODEC_FORCE_MATERIALIZED_DECODE;
    require_decode_path(high, kDecodePathMaterialized,
        kDecodeRuleForcedMaterialized,
        "forced materialized batch precedence");

    DecodePathInput low = make_decode_path_input(
        LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8, LEO2_BACKEND_AVX2,
        8, 248, 8, 256, 8, 65);
    selection = require_decode_path(low, kDecodePathTiled,
        kDecodeRuleWorkspaceTiled, "low tiled fallback");
    require(selection.required_work_slots == 16,
        "low tiled work slots differ");
    low = make_decode_path_input(LEO2_PROFILE_LOW_V1,
        LEO2_FIELD_GF8, LEO2_BACKEND_AVX2,
        8, 8, 8, 16, 8, 64);
    require_decode_path(low, kDecodePathMaterialized,
        kDecodeRuleWorkspaceMaterialized,
        "workspace equality must retain materialized");
    low.codec_flags = LEO2_CODEC_FORCE_TILED_DECODE;
    require_decode_path(low, kDecodePathTiled,
        kDecodeRuleForcedTiled, "forced workspace equality");

    DecodePathInput fallback = make_decode_path_input(
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16,
        LEO2_BACKEND_AVX2, 224, 32, 32, 256, 1, 32 * 1024);
    require_decode_path(fallback, kDecodePathTiled,
        kDecodeRuleWorkspaceTiled, "GF16 measured-rule fallback");
    fallback.field = LEO2_FIELD_GF8;
    fallback.backend = LEO2_BACKEND_NEON;
    require_decode_path(fallback, kDecodePathTiled,
        kDecodeRuleWorkspaceTiled, "NEON measured-rule fallback");

    DecodePathInput unsupported = make_decode_path_input(
        LEO2_PROFILE_EXACT_EXPERIMENTAL_V1, LEO2_FIELD_GF8,
        LEO2_BACKEND_AVX2, 224, 32, 32, 256, 8, 32 * 1024);
    selection = require_decode_path(unsupported, kDecodePathMaterialized,
        kDecodeRuleUnsupportedProfile,
        "unsupported-profile retained fallback");
    require(selection.required_work_slots == 256,
        "unsupported-profile fallback did not retain N slots");
    unsupported.codec_flags = LEO2_CODEC_FORCE_TILED_DECODE;
    selection = require_decode_path(unsupported, kDecodePathTiled,
        kDecodeRuleForcedTiled,
        "unsupported-profile forced tiled traversal");
    require(selection.required_work_slots == 256,
        "unsupported forced tiled path did not retain fallback N slots");

    DecodePathInput codec_query = make_decode_path_input(
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        LEO2_BACKEND_AVX2, 224, 32, 32, 256, 0, 32 * 1024,
        0, false, false);
    require_decode_path(codec_query, kDecodePathMaterialized,
        kDecodeRuleMeasuredMaterialized,
        "one-shot conservative materialized query");

    DecodePathInput invalid = balanced;
    invalid.tail_bytes = 0;
    DecodePathInfo unused;
    require(!SelectDecodePath(invalid, unused),
        "inconsistent tail geometry was accepted");
    invalid = balanced;
    invalid.codec_flags = LEO2_CODEC_FORCE_GENERIC_DECODE |
        LEO2_CODEC_FORCE_TILED_DECODE;
    require(!SelectDecodePath(invalid, unused),
        "conflicting force precedence was accepted");

    const uint32_t counts[] = { 31, 32, 33, 127, 128, 129, 223, 224, 225 };
    const uint64_t byte_counts[] = {
        192, 193, 256, 24 * 1024 - 64, 24 * 1024 - 63,
        32 * 1024, 64 * 1024, 64 * 1024 + 1,
        1024 * 1024, 1024 * 1024 + 1
    };
    const leo2_backend backends[] = {
        LEO2_BACKEND_SCALAR, LEO2_BACKEND_SSSE3,
        LEO2_BACKEND_AVX2, LEO2_BACKEND_AVX512, LEO2_BACKEND_NEON,
        LEO2_BACKEND_GFNI
    };
    for (size_t k_i = 0; k_i < sizeof(counts) / sizeof(counts[0]); ++k_i)
    for (size_t r_i = 0; r_i < sizeof(counts) / sizeof(counts[0]); ++r_i)
    {
        const uint32_t k = counts[k_i];
        const uint32_t r = counts[r_i];
        const uint32_t padded = test_ceil_pow2(r);
        const uint32_t parent = test_ceil_pow2(k + padded);
        if (parent > 256)
            continue;
        const uint32_t missing = std::min(k, r);
        for (size_t byte_i = 0;
             byte_i < sizeof(byte_counts) / sizeof(byte_counts[0]); ++byte_i)
        for (size_t backend_i = 0;
             backend_i < sizeof(backends) / sizeof(backends[0]); ++backend_i)
        {
            DecodePathInput input = make_decode_path_input(
                LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
                backends[backend_i], k, r, padded, parent, missing,
                byte_counts[byte_i]);
            DecodePathInfo result;
            require(SelectDecodePath(input, result),
                "overlap sweep rejected a valid cell");
            require(result.matching_auto_rules !=
                    (kDecodeAutoRuleBalancedGeneric |
                     kDecodeAutoRuleMeasuredMaterialized),
                "evidence-scoped AUTO rules overlap");
            require_batch_decode_scratch_bound(
                input, "high-profile batch scratch invariant");

            const uint32_t low_padded = test_ceil_pow2(k);
            const uint32_t low_parent = test_ceil_pow2(low_padded + r);
            if (low_parent <= 256)
            {
                DecodePathInput low_input = make_decode_path_input(
                    LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8,
                    backends[backend_i], k, r, low_padded, low_parent,
                    missing, byte_counts[byte_i]);
                require_batch_decode_scratch_bound(
                    low_input, "low-profile batch scratch invariant");
            }
        }
    }
}

void require_balanced_family(
    bool condition,
    unsigned k,
    size_t bytes,
    const char* message)
{
    if (condition)
        return;
    std::ostringstream stream;
    stream << "balanced family K=R=" << k;
    if (bytes != 0)
        stream << " bytes=" << bytes;
    stream << ": " << message;
    throw std::runtime_error(stream.str());
}

void test_balanced_family_forced_equivalence(leo2_context* context)
{
    const size_t execution_bytes = 193;
    const size_t scratch_sizes[] = { execution_bytes, 256, 4097 };
    const uint32_t flags[] = {
        0,
        LEO2_CODEC_FORCE_GENERIC_DECODE,
        LEO2_CODEC_FORCE_SPECIALIZED_DECODE |
            LEO2_CODEC_FORCE_MATERIALIZED_DECODE,
        LEO2_CODEC_FORCE_SPECIALIZED_DECODE |
            LEO2_CODEC_FORCE_TILED_DECODE
    };
    const size_t path_count = sizeof(flags) / sizeof(flags[0]);

    for (unsigned k = 5; k <= 128; ++k)
    {
        leo2_codec* codecs[path_count] = { NULL, NULL, NULL, NULL };
        leo2_decode_plan* plans[path_count] = { NULL, NULL, NULL, NULL };
        std::vector<uint8_t> original_present(k, 0);
        std::vector<uint8_t> recovery_present(k, 1);

        for (size_t path = 0; path < path_count; ++path)
        {
            codecs[path] = flags[path] == 0
                ? make_codec(context, k, k,
                    LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8)
                : make_flagged_codec(context, k, k,
                    LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
                    flags[path]);
            require_balanced_family(
                leo2_codec_field(codecs[path]) == LEO2_FIELD_GF8,
                k, 0, "resolved outside GF8");
            require_balanced_family(
                leo2_codec_parent_count(codecs[path]) ==
                    leo2_codec_padded_side(codecs[path]) * 2,
                k, 0, "did not resolve to N=2T");
            require_result(leo2_decode_plan_create(
                codecs[path], &original_present[0], &recovery_present[0],
                &plans[path]), "balanced family plan create");
        }

        size_t execution_scratch[path_count] = { 0, 0, 0, 0 };
        for (size_t size_i = 0;
             size_i < sizeof(scratch_sizes) / sizeof(scratch_sizes[0]);
             ++size_i)
        {
            const size_t bytes = scratch_sizes[size_i];
            size_t scratch[path_count] = { 0, 0, 0, 0 };
            for (size_t path = 0; path < path_count; ++path)
            {
                require_result(leo2_decode_plan_scratch_size(
                    plans[path], bytes, &scratch[path]),
                    "balanced family scratch query");
            }
            const bool experimental_auto_direct =
                LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE != 0 &&
                leo2_context_backend(context) == LEO2_BACKEND_AVX2 &&
                k >= 5 && k <= 8;
            require_balanced_family(
                scratch[1] == scratch[2] &&
                (experimental_auto_direct
                    ? scratch[0] < scratch[1]
                    : scratch[0] == scratch[1]),
                k, bytes,
                experimental_auto_direct
                    ? "experimental AUTO did not use bounded direct scratch"
                    : "AUTO, generic, and materialized scratch differ");
            require_balanced_family(
                scratch[3] == scratch[2], k, bytes,
                "translated tiled and materialized scratch differ at N=2P");

            size_t one_shot_scratch = 0;
            require_result(leo2_decode_scratch_size(
                codecs[0], bytes, &one_shot_scratch),
                "balanced family one-shot scratch query");
            require_balanced_family(
                experimental_auto_direct
                    ? one_shot_scratch == scratch[1] &&
                        one_shot_scratch > scratch[0]
                    : one_shot_scratch == scratch[0],
                k, bytes,
                experimental_auto_direct
                    ? "one-shot query did not retain conservative scratch"
                    : "one-shot query does not cover the full-loss AUTO plan");
            if (bytes == execution_bytes)
                std::copy(scratch, scratch + path_count, execution_scratch);
        }

        const Shards source = make_originals(
            k, execution_bytes, UINT64_C(0xb41a4ced00000000) + k);
        const Shards parity = encode_new(codecs[0], source, execution_bytes);
        std::vector<const void*> original_input(k, NULL);
        const std::vector<const void*> recovery_input = const_pointers(parity);
        for (size_t path = 0; path < path_count; ++path)
        {
            Shards restored(k, std::vector<uint8_t>(execution_bytes, 0));
            std::vector<void*> output = mutable_pointers(restored);
            AlignedBuffer scratch(execution_scratch[path]);
            require_result(leo2_decode_plan_execute(
                plans[path], execution_bytes, &original_input[0],
                &recovery_input[0], &output[0], scratch.data, scratch.bytes),
                "balanced family decode execute");
            require_balanced_family(restored == source, k, execution_bytes,
                "AUTO/generic/materialized/tiled recovery mismatch");
        }

        for (size_t path = 0; path < path_count; ++path)
        {
            leo2_decode_plan_destroy(plans[path]);
            leo2_codec_destroy(codecs[path]);
        }
    }
}

} // namespace

int main()
{
    try
    {
        TestCounts counts;
        leo2_context* context = NULL;
        leo2_context_options options;
        memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AUTO;
        require_result(leo2_context_create(&options, &context), "context create");
        require(context != NULL && leo2_context_backend(context) != LEO2_BACKEND_AUTO,
            "context backend introspection failed");
        test_forced_backend(context);
        test_codec_flag_validation(context);
        test_balanced_dispatch_policy();
        test_unified_decode_path_selector();
        test_balanced_family_forced_equivalence(context);
        test_tiled_high_dispatch_policy();
        test_tiled_materialized_execution(context);
        test_available_avx512_batch_dispatch();
        test_low_reveal_fusion(context);
        test_high_reveal_fusion(context);
        test_batch_materialized_capacity(context);

        test_profile_metadata(context);
        compare_high_with_legacy(context, 3, 2,
            std::vector<size_t>{1, 17, 64, 65, 257}, &counts);
        compare_high_with_legacy(context, 9, 7,
            std::vector<size_t>{3, 64, 129}, &counts);
        compare_high_with_legacy(context, 129, 2,
            std::vector<size_t>{7, 64}, &counts);
        compare_high_with_legacy(context, 240, 16,
            std::vector<size_t>{64}, &counts);
        compare_high_with_legacy(context, 1000, 200,
            std::vector<size_t>{64}, &counts);
        compare_avx2_fused_high_with_legacy(&counts);

        const BinaryField gf8 = leopard2_test::make_legacy_gf8();
        compare_low_gf8_with_oracle(context, 2, 5, 65, gf8, &counts);
        compare_low_gf8_with_oracle(context, 3, 5, 129, gf8, &counts);
        compare_low_gf8_with_oracle(context, 5, 11, 17, gf8, &counts);
        const BinaryField gf16 = leopard2_test::make_legacy_gf16();
        for (size_t bytes = 2; bytes <= 64; bytes += 2)
            compare_low_gf16_with_oracle(context, gf16, bytes, &counts);
        compare_low_gf16_with_oracle(context, gf16, 66, &counts);
        compare_low_gf16_with_oracle(context, gf16, 1024, &counts);
        compare_low_gf16_with_oracle(context, gf16, 1026, &counts);

        run_decode_case(context, 9, 7, LEO2_PROFILE_LEGACY_HIGH_V1,
            LEO2_FIELD_GF8, 17, std::vector<unsigned>{0, 4, 8},
            std::vector<unsigned>{5, 6}, &counts);
        run_decode_case(context, 3, 5, LEO2_PROFILE_LOW_V1,
            LEO2_FIELD_GF8, 65, std::vector<unsigned>{0, 2},
            std::vector<unsigned>{1, 3, 4}, &counts);
        run_decode_case(context, 3, 5, LEO2_PROFILE_LOW_V1,
            LEO2_FIELD_GF8, 129, std::vector<unsigned>{0, 1, 2},
            std::vector<unsigned>{3, 4}, &counts);
        run_decode_case(context, 257, 33, LEO2_PROFILE_LEGACY_HIGH_V1,
            LEO2_FIELD_GF16, 64, std::vector<unsigned>{1, 64, 128, 256},
            std::vector<unsigned>{29, 30, 31, 32}, &counts);
        const size_t gf16_boundaries[] = { 2, 32, 34, 62, 64, 66, 1024, 1026 };
        for (size_t count_i = 0;
             count_i < sizeof(gf16_boundaries) / sizeof(gf16_boundaries[0]);
             ++count_i)
        {
            const size_t bytes = gf16_boundaries[count_i];
            compare_high_gf16_compact_with_legacy(
                context, 257, 33, bytes, &counts);
            run_decode_case(context, 257, 33, LEO2_PROFILE_LEGACY_HIGH_V1,
                LEO2_FIELD_GF16, bytes, std::vector<unsigned>{1, 128, 256},
                std::vector<unsigned>{31, 32}, &counts);
            run_decode_case(context, 100, 156, LEO2_PROFILE_LOW_V1,
                LEO2_FIELD_GF16, bytes, std::vector<unsigned>{0, 37, 99},
                std::vector<unsigned>{151, 152, 155}, &counts);
        }
        run_decode_case(context, 257, 1, LEO2_PROFILE_LEGACY_HIGH_V1,
            LEO2_FIELD_GF16, 66, std::vector<unsigned>{128},
            std::vector<unsigned>(), &counts);
        run_decode_case(context, 1, 300, LEO2_PROFILE_LOW_V1,
            LEO2_FIELD_GF16, 66, std::vector<unsigned>{0},
            std::vector<unsigned>{0, 299}, &counts);
        run_decode_case(context, 8, 1, LEO2_PROFILE_LEGACY_HIGH_V1,
            LEO2_FIELD_GF8, 33, std::vector<unsigned>{3},
            std::vector<unsigned>(), &counts);
        run_decode_case(context, 1, 5, LEO2_PROFILE_LOW_V1,
            LEO2_FIELD_GF8, 31, std::vector<unsigned>{0},
            std::vector<unsigned>{0, 1, 3, 4}, &counts);

        std::vector<unsigned> balanced_full_recovery(128);
        for (size_t i = 0; i < balanced_full_recovery.size(); ++i)
            balanced_full_recovery[i] = static_cast<unsigned>(i);
        leo2_test_reset_generic_reveal_counts();
        leo2_test_reset_low_reveal_counts();
        run_decode_case(context, 128, 128, LEO2_PROFILE_LEGACY_HIGH_V1,
            LEO2_FIELD_GF8, 257, balanced_full_recovery,
            std::vector<unsigned>(), &counts);
        require(leo2_test_generic_direct_reveal_shards() == 0,
            "small balanced decode unexpectedly fused reveal/scatter");
        require(leo2_test_low_direct_reveal_shards() == 128 * 3,
            "small balanced AUTO did not use translated Algorithm 4 reveal");

        leo2_test_reset_generic_reveal_counts();
        leo2_test_reset_low_reveal_counts();
        run_decode_case(context, 128, 128, LEO2_PROFILE_LEGACY_HIGH_V1,
            LEO2_FIELD_GF8, 4097, balanced_full_recovery,
            std::vector<unsigned>(), &counts);
        const leo2_backend balanced_backend = leo2_context_backend(context);
        const uint64_t expected_direct_reveals =
            balanced_backend == LEO2_BACKEND_SSSE3 ||
            balanced_backend == LEO2_BACKEND_AVX2 ||
            balanced_backend == LEO2_BACKEND_AVX512
                ? 128
                : 0;
        require(leo2_test_generic_direct_reveal_shards() ==
                expected_direct_reveals,
            "forced generic reveal/scatter did not match backend policy");
        require(leo2_test_low_direct_reveal_shards() == 128 * 3,
            "balanced AUTO did not use translated Algorithm 4 reveal");

        const size_t direct_gf16_boundaries[] = { 2, 34, 64, 66, 1026 };
        for (size_t count_i = 0;
             count_i < sizeof(direct_gf16_boundaries) /
                 sizeof(direct_gf16_boundaries[0]);
             ++count_i)
        {
            const size_t bytes = direct_gf16_boundaries[count_i];
            run_decode_case(context, 9, 7, LEO2_PROFILE_LEGACY_HIGH_V1,
                LEO2_FIELD_GF16, bytes, std::vector<unsigned>{0, 4, 8},
                std::vector<unsigned>{0, 2, 6}, &counts);
            run_decode_case(context, 5, 11, LEO2_PROFILE_LOW_V1,
                LEO2_FIELD_GF16, bytes, std::vector<unsigned>{0, 1, 3, 4},
                std::vector<unsigned>{0, 2, 10}, &counts);
        }

        test_no_loss_no_op(context);
        test_direct_repair_dispatch_bounds(context);
        test_expanded_direct_repair_execution(&counts);
        test_direct_repair_field_helpers();
        test_overlap_rejection(context);
        test_decode_late_scratch_alias_failure_atomicity(context);
        test_gf16_byte_granularity(context);
        leo2_context_destroy(context);

        std::cout << "high_compatibility_bytes=" << counts.high_compatibility
                  << " low_oracle_symbols=" << counts.low_oracle_symbols
                  << " recovered_shards=" << counts.recovered_shards
                  << " plan_executions=" << counts.plan_executions
                  << " tail_cases=" << counts.tail_cases << std::endl;
        return 0;
    }
    catch (const std::exception& exception)
    {
        std::cerr << "leopard2 API test failed: " << exception.what() << std::endl;
        return 1;
    }
}
