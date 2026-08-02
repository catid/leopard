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
#include "leopard2.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

namespace {

static const unsigned kOriginalCount = 32;
static const unsigned kRecoveryCount = 32;
#ifndef LEO2_BALANCED_PRODUCTION_SHARD_BYTES
#define LEO2_BALANCED_PRODUCTION_SHARD_BYTES 64
#endif
#ifndef LEO2_EXPECT_T32_B256_GENERATED
#define LEO2_EXPECT_T32_B256_GENERATED 0
#endif
static const size_t kShardBytes = LEO2_BALANCED_PRODUCTION_SHARD_BYTES;
static const size_t kGuardBytes = 64;
static const uint8_t kGuardValue = 0x6d;

void Require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void RequireResult(
    leo2_result actual,
    leo2_result expected,
    const char* message)
{
    if (actual != expected)
    {
        throw std::runtime_error(std::string(message) + ": expected " +
            leo2_result_string(expected) + ", received " +
            leo2_result_string(actual));
    }
}

class AlignedBuffer
{
public:
    explicit AlignedBuffer(size_t bytes)
        : value_(NULL)
        , bytes_(bytes)
    {
#if defined(_MSC_VER)
        value_ = _aligned_malloc(bytes, leo2_scratch_alignment());
#else
        if (posix_memalign(&value_, leo2_scratch_alignment(), bytes) != 0)
            value_ = NULL;
#endif
        if (!value_)
            throw std::bad_alloc();
        std::memset(value_, 0, bytes);
    }

    ~AlignedBuffer()
    {
#if defined(_MSC_VER)
        _aligned_free(value_);
#else
        std::free(value_);
#endif
    }

    uint8_t* bytes() const { return static_cast<uint8_t*>(value_); }
    void* data() const { return value_; }
    size_t size() const { return bytes_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);
    void* value_;
    size_t bytes_;
};

void FillInput(uint8_t* input)
{
    uint64_t state = UINT64_C(0x54333250524f4436);
    for (size_t i = 0; i < kOriginalCount * kShardBytes; ++i)
    {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        input[i] = static_cast<uint8_t>(state >> 24);
    }
}

void SetPackedPointers(
    uint8_t* input,
    uint8_t* output,
    const void** original,
    void** recovery)
{
    for (unsigned i = 0; i < kOriginalCount; ++i)
        original[i] = input + static_cast<size_t>(i) * kShardBytes;
    for (unsigned i = 0; i < kRecoveryCount; ++i)
        recovery[i] = output + static_cast<size_t>(i) * kShardBytes;
}

void CheckParity(
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator,
    const void* const* original,
    void* const* recovery)
{
    for (unsigned parity = 0; parity < kRecoveryCount; ++parity)
    {
        const std::vector<leopard2_test::Element>& row =
            generator[kOriginalCount + parity];
        const uint8_t* const output =
            static_cast<const uint8_t*>(recovery[parity]);
        for (size_t offset = 0; offset < kShardBytes; ++offset)
        {
            leopard2_test::Element expected = 0;
            for (unsigned source = 0; source < kOriginalCount; ++source)
            {
                const uint8_t value =
                    static_cast<const uint8_t*>(original[source])[offset];
                expected = field.add(expected,
                    field.multiply(row[source], value));
            }
            Require(output[offset] == expected,
                "production parity differs from the independent oracle");
        }
    }
}

size_t ExpectedProductionScratch(unsigned side)
{
    const size_t alignment = leo2_scratch_alignment();
    const size_t metadata_bytes =
        4U * side * sizeof(uintptr_t) + 3U * side * sizeof(void*);
    const size_t data_offset =
        (metadata_bytes + alignment - 1U) & ~(alignment - 1U);
    return data_offset + 2U * side * kShardBytes;
}

void CheckGuards(
    const AlignedBuffer& buffer,
    size_t payload_offset,
    size_t payload_bytes,
    const char* message)
{
    for (size_t i = 0; i < payload_offset; ++i)
        Require(buffer.bytes()[i] == kGuardValue, message);
    for (size_t i = payload_offset + payload_bytes;
         i < buffer.size(); ++i)
        Require(buffer.bytes()[i] == kGuardValue, message);
}

void ExerciseProduction(leo2_context* context)
{
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context,
        kOriginalCount, kRecoveryCount,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        NULL, &codec), LEO2_SUCCESS, "create production K32/R32 codec");
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(codec, kShardBytes,
        &scratch_bytes), LEO2_SUCCESS,
        "query production K32/R32 scratch");
    Require(scratch_bytes == ExpectedProductionScratch(kOriginalCount),
        "production scratch differs from the portable fixed geometry");
    AlignedBuffer scratch(scratch_bytes + leo2_scratch_alignment());
#if LEO2_EXPECT_T32_B256_GENERATED
    std::memset(scratch.bytes(), 0x39, scratch.size());
    const std::vector<uint8_t> scratch_before(
        scratch.bytes(), scratch.bytes() + scratch.size());
#endif

    const size_t input_bytes = kOriginalCount * kShardBytes;
    const size_t output_bytes = kRecoveryCount * kShardBytes;
    AlignedBuffer input(input_bytes + 2U * kGuardBytes + 8U);
    AlignedBuffer output(output_bytes + 2U * kGuardBytes + 8U);
    std::memset(input.bytes(), kGuardValue, input.size());
    std::memset(output.bytes(), kGuardValue, output.size());
    const size_t input_offset = kGuardBytes + 1U;
    const size_t output_offset = kGuardBytes + 3U;
    uint8_t* const input_base = input.bytes() + input_offset;
    uint8_t* const output_base = output.bytes() + output_offset;
    FillInput(input_base);
    const std::vector<uint8_t> input_before(
        input.bytes(), input.bytes() + input.size());
    const void* original[kOriginalCount];
    void* recovery[kRecoveryCount];
    SetPackedPointers(input_base, output_base, original, recovery);

    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh,
            kOriginalCount, kRecoveryCount);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);

    RequireResult(leo2_encode(codec, kShardBytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute production packed terminal");
#if LEO2_EXPECT_T32_B256_GENERATED
    Require(std::memcmp(scratch.bytes(), &scratch_before[0],
            scratch_before.size()) == 0,
        "generated T32/B256 route unexpectedly staged scratch data");
#endif
    Require(std::memcmp(input.bytes(), &input_before[0],
            input_before.size()) == 0,
        "production encode modified source data or guards");
    CheckGuards(output, output_offset, output_bytes,
        "production encode modified a destination guard");
    CheckParity(field, generator, original, recovery);

    std::memset(output.bytes(), kGuardValue, output.size());
    leo2_encode_batch_item item = {
        kShardBytes, original, recovery, scratch.data(), scratch_bytes
    };
    RequireResult(leo2_encode_batch(codec, &item, 1), LEO2_SUCCESS,
        "execute production one-item batch terminal");
#if LEO2_EXPECT_T32_B256_GENERATED
    Require(std::memcmp(scratch.bytes(), &scratch_before[0],
            scratch_before.size()) == 0,
        "generated T32/B256 batch unexpectedly staged scratch data");
#endif
    Require(std::memcmp(input.bytes(), &input_before[0],
            input_before.size()) == 0,
        "production batch modified source data or guards");
    CheckGuards(output, output_offset, output_bytes,
        "production batch modified a destination guard");
    CheckParity(field, generator, original, recovery);

    std::memset(output.bytes(), kGuardValue, output.size());
    const std::vector<uint8_t> output_before(
        output.bytes(), output.bytes() + output.size());
    RequireResult(leo2_encode(codec, kShardBytes, original, recovery,
        scratch.data(), scratch_bytes - 1U), LEO2_SCRATCH_TOO_SMALL,
        "reject undersized production scratch");
    Require(std::memcmp(output.bytes(), &output_before[0],
            output_before.size()) == 0,
        "undersized production scratch modified output");
    RequireResult(leo2_encode(codec, kShardBytes, original, recovery,
        scratch.bytes() + 1U, scratch_bytes), LEO2_BAD_ALIGNMENT,
        "reject misaligned production scratch");
    Require(std::memcmp(output.bytes(), &output_before[0],
            output_before.size()) == 0,
        "misaligned production scratch modified output");

    leo2_codec_destroy(codec);
}

void ExerciseProductionBalancedSide(leo2_context* context, unsigned side)
{
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, side, side,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        LEO2_SUCCESS, "create production balanced codec");
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(codec, kShardBytes,
        &scratch_bytes), LEO2_SUCCESS,
        "query production balanced scratch");
    Require(scratch_bytes == ExpectedProductionScratch(side),
        "balanced production scratch differs from fixed geometry");
    AlignedBuffer scratch(scratch_bytes);

    const size_t payload_bytes = static_cast<size_t>(side) * kShardBytes;
    AlignedBuffer input(payload_bytes + 2U * kGuardBytes + 8U);
    AlignedBuffer output(payload_bytes + 2U * kGuardBytes + 8U);
    std::memset(input.bytes(), kGuardValue, input.size());
    std::memset(output.bytes(), kGuardValue, output.size());
    const size_t input_offset = kGuardBytes + 1U;
    const size_t output_offset = kGuardBytes + 3U;
    uint8_t* const input_base = input.bytes() + input_offset;
    uint8_t* const output_base = output.bytes() + output_offset;
    uint64_t state = UINT64_C(0x50524f4442414c36) ^ side;
    for (size_t i = 0; i < payload_bytes; ++i)
    {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        input_base[i] = static_cast<uint8_t>(state >> 24);
    }
    const std::vector<uint8_t> input_before(
        input.bytes(), input.bytes() + input.size());
    std::vector<const void*> original(side);
    std::vector<void*> recovery(side);
    for (unsigned i = 0; i < side; ++i)
    {
        original[i] = input_base + static_cast<size_t>(i) * kShardBytes;
        recovery[i] = output_base + static_cast<size_t>(i) * kShardBytes;
    }

    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh, side, side);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);

    RequireResult(leo2_encode(codec, kShardBytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute production balanced terminal");
    Require(std::memcmp(input.bytes(), &input_before[0], input.size()) == 0,
        "balanced production encode modified source or guards");
    CheckGuards(output, output_offset, payload_bytes,
        "balanced production encode modified output guard");

    for (unsigned parity = 0; parity < side; ++parity)
    {
        const std::vector<leopard2_test::Element>& row =
            generator[side + parity];
        const uint8_t* const encoded =
            static_cast<const uint8_t*>(recovery[parity]);
        for (size_t offset = 0; offset < kShardBytes; ++offset)
        {
            leopard2_test::Element expected = 0;
            for (unsigned source = 0; source < side; ++source)
            {
                expected = field.add(expected, field.multiply(row[source],
                    static_cast<const uint8_t*>(original[source])[offset]));
            }
            Require(encoded[offset] == expected,
                "balanced production parity differs from direct oracle");
        }
    }

    std::memset(output.bytes(), kGuardValue, output.size());
    leo2_encode_batch_item item = {
        kShardBytes, &original[0], &recovery[0], scratch.data(), scratch.size()
    };
    RequireResult(leo2_encode_batch(codec, &item, 1), LEO2_SUCCESS,
        "execute production balanced one-item batch");
    CheckGuards(output, output_offset, payload_bytes,
        "balanced production batch modified output guard");
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
        const leo2_result result = leo2_context_create(&options, &context);
        if (result == LEO2_UNSUPPORTED)
        {
            std::printf(
                "production balanced 64-byte terminal test skipped: AVX2 unavailable\n");
            return 0;
        }
        RequireResult(result, LEO2_SUCCESS,
            "create production AVX2 context");
        ExerciseProduction(context);
        ExerciseProductionBalancedSide(context, 64);
        ExerciseProductionBalancedSide(context, 128);
        leo2_context_destroy(context);
        std::printf("production balanced 64-byte terminal checks passed\n");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr,
            "production balanced 64-byte terminal failure: %s\n",
            error.what());
        return 1;
    }
}
