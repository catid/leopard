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
#if defined(LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK)
#include "Leopard2Backend.h"
#endif

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
#ifndef LEO2_EXPECT_T32_B256_TWO_BLOCK
#define LEO2_EXPECT_T32_B256_TWO_BLOCK 0
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

#if defined(LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK)
void CheckRegionGuards(
    const AlignedBuffer& buffer,
    size_t region_offset,
    size_t region_bytes,
    size_t payload_offset,
    size_t payload_bytes,
    const char* message)
{
    Require(payload_offset >= region_offset &&
        payload_offset + payload_bytes <= region_offset + region_bytes,
        "invalid guarded-region geometry");
    for (size_t i = region_offset; i < payload_offset; ++i)
        Require(buffer.bytes()[i] == kGuardValue, message);
    for (size_t i = payload_offset + payload_bytes;
         i < region_offset + region_bytes; ++i)
        Require(buffer.bytes()[i] == kGuardValue, message);
}
#endif

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

#if defined(LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK)
void CheckVariableParity(
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator,
    unsigned original_count,
    unsigned recovery_count,
    const std::vector<const void*>& original,
    const std::vector<void*>& recovery)
{
    for (unsigned parity = 0; parity < recovery_count; ++parity)
    {
        const std::vector<leopard2_test::Element>& row =
            generator[original_count + parity];
        const uint8_t* const encoded =
            static_cast<const uint8_t*>(recovery[parity]);
        for (size_t offset = 0; offset < kShardBytes; ++offset)
        {
            leopard2_test::Element expected = 0;
            for (unsigned source = 0; source < original_count; ++source)
            {
                expected = field.add(expected, field.multiply(row[source],
                    static_cast<const uint8_t*>(
                        original[source])[offset]));
            }
            Require(encoded[offset] == expected,
                "T32 multi-block parity differs from direct oracle");
        }
    }
}

void ExerciseProductionT32TwoBlockFamily(
    leo2_context* context,
    unsigned original_count)
{
    static const unsigned recovery_count = 32;
    Require(kShardBytes == 256,
        "T32 multi-block production test requires 256-byte shards");
    Require(original_count == 64 || original_count == 96 ||
            original_count == 128,
        "invalid T32 multi-block production shape");

    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, original_count, recovery_count,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        LEO2_SUCCESS, "create T32 multi-block codec");
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(codec, kShardBytes,
        &scratch_bytes), LEO2_SUCCESS,
        "query T32 multi-block scratch");
    const size_t alignment = leo2_scratch_alignment();
    AlignedBuffer scratch(scratch_bytes + 2U * alignment);
    std::memset(scratch.bytes(), kGuardValue, scratch.size());
    void* const scratch_payload = scratch.bytes() + alignment;

    const size_t input_bytes =
        static_cast<size_t>(original_count) * kShardBytes;
    const size_t output_bytes =
        static_cast<size_t>(recovery_count) * kShardBytes;
    AlignedBuffer input(input_bytes + 2U * kGuardBytes + 8U);
    AlignedBuffer output(output_bytes + 2U * kGuardBytes + 8U);
    std::memset(input.bytes(), kGuardValue, input.size());
    std::memset(output.bytes(), kGuardValue, output.size());
    const size_t input_offset = kGuardBytes + 1U;
    const size_t output_offset = kGuardBytes + 3U;
    uint8_t* const input_base = input.bytes() + input_offset;
    uint8_t* const output_base = output.bytes() + output_offset;
    uint64_t state = UINT64_C(0x5433324d554c5449) ^ original_count;
    for (size_t i = 0; i < input_bytes; ++i)
    {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        input_base[i] = static_cast<uint8_t>(state >> 24);
    }
    const std::vector<uint8_t> input_before(
        input.bytes(), input.bytes() + input.size());
    std::vector<const void*> original(original_count);
    std::vector<void*> recovery(recovery_count);
    for (unsigned i = 0; i < original_count; ++i)
        original[i] = input_base + static_cast<size_t>(i) * kShardBytes;
    for (unsigned i = 0; i < recovery_count; ++i)
        recovery[i] = output_base + static_cast<size_t>(i) * kShardBytes;

    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh,
            original_count, recovery_count);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);

    RequireResult(leo2_encode(codec, kShardBytes, &original[0], &recovery[0],
        scratch_payload, scratch_bytes), LEO2_SUCCESS,
        "execute public T32 multi-block encode");
    Require(std::memcmp(input.bytes(), &input_before[0], input.size()) == 0,
        "T32 multi-block encode modified source or guards");
    CheckGuards(output, output_offset, output_bytes,
        "T32 multi-block encode modified output guard");
    CheckGuards(scratch, alignment, scratch_bytes,
        "T32 multi-block encode modified scratch guard");
    CheckVariableParity(field, generator, original_count, recovery_count,
        original, recovery);

    std::memset(output.bytes(), kGuardValue, output.size());
    leo2_encode_batch_item item = {
        kShardBytes, &original[0], &recovery[0],
        scratch_payload, scratch_bytes
    };
    RequireResult(leo2_encode_batch(codec, &item, 1), LEO2_SUCCESS,
        "execute public T32 multi-block batch encode");
    CheckGuards(output, output_offset, output_bytes,
        "T32 multi-block batch modified output guard");
    CheckVariableParity(field, generator, original_count, recovery_count,
        original, recovery);

    AlignedBuffer direct_output(output_bytes + 2U * kGuardBytes + 8U);
    AlignedBuffer temporary(output_bytes + 2U * kGuardBytes + 8U);
    std::memset(direct_output.bytes(), kGuardValue, direct_output.size());
    std::memset(temporary.bytes(), kGuardValue, temporary.size());
    const size_t direct_offset = kGuardBytes + 3U;
    const size_t temporary_offset = kGuardBytes + 5U;
    uint8_t* const direct_base = direct_output.bytes() + direct_offset;
    uint8_t* const temporary_base = temporary.bytes() + temporary_offset;
    std::vector<void*> direct_recovery(recovery_count);
    for (unsigned i = 0; i < recovery_count; ++i)
    {
        direct_recovery[i] =
            direct_base + static_cast<size_t>(i) * kShardBytes;
    }
    const bool selected =
        original_count == 64U &&
        leopard::backend::TryAVX2FF8HighEncodeT32B256TwoBlockPacked(
            input_base, direct_base, temporary_base);
#if LEO2_EXPECT_T32_B256_TWO_BLOCK
    Require(selected == (original_count == 64U),
        "T32 two-block callback selector disagrees with its exact shape");
    if (selected)
    {
        CheckGuards(direct_output, direct_offset, output_bytes,
            "T32 two-block callback modified output guard");
        CheckGuards(temporary, temporary_offset, output_bytes,
            "T32 two-block callback modified temporary guard");
        CheckVariableParity(field, generator, original_count, recovery_count,
            original, direct_recovery);
    }
#else
    Require(!selected,
        "disabled T32 two-block callback accepted its exact shape");
#endif

    std::memset(output.bytes(), kGuardValue, output.size());
    const std::vector<uint8_t> output_before(
        output.bytes(), output.bytes() + output.size());
    RequireResult(leo2_encode(codec, kShardBytes, &original[0], &recovery[0],
        scratch_payload, scratch_bytes - 1U), LEO2_SCRATCH_TOO_SMALL,
        "reject undersized T32 multi-block scratch");
    RequireResult(leo2_encode(codec, kShardBytes, &original[0], &recovery[0],
        static_cast<uint8_t*>(scratch_payload) + 1U, scratch_bytes),
        LEO2_BAD_ALIGNMENT, "reject misaligned T32 multi-block scratch");
    Require(std::memcmp(output.bytes(), &output_before[0], output.size()) == 0,
        "rejected T32 multi-block scratch modified output");

    void* const saved_recovery = recovery[0];
    recovery[0] = input_base;
    RequireResult(leo2_encode(codec, kShardBytes, &original[0], &recovery[0],
        scratch_payload, scratch_bytes), LEO2_OVERLAP,
        "reject T32 multi-block source/output overlap");
    recovery[0] = saved_recovery;
    Require(std::memcmp(input.bytes(), &input_before[0], input.size()) == 0,
        "rejected T32 multi-block overlap modified input");
    Require(std::memcmp(output.bytes(), &output_before[0], output.size()) == 0,
        "rejected T32 multi-block overlap modified output");

    leo2_codec_destroy(codec);
}

void ExerciseProductionT32TwoBlockBatch(leo2_context* context)
{
    static const size_t kBatchCount = 8;
    static const unsigned original_count = 64;
    static const unsigned recovery_count = 32;
    const size_t alignment = leo2_scratch_alignment();
    const size_t input_bytes = original_count * kShardBytes;
    const size_t output_bytes = recovery_count * kShardBytes;
    const size_t input_offset_in_region = kGuardBytes + 1U;
    const size_t output_offset_in_region = kGuardBytes + 3U;
    const size_t input_region_bytes =
        (input_offset_in_region + input_bytes + kGuardBytes + 8U +
            alignment - 1U) & ~(alignment - 1U);
    const size_t output_region_bytes =
        (output_offset_in_region + output_bytes + kGuardBytes + 8U +
            alignment - 1U) & ~(alignment - 1U);

    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, original_count, recovery_count,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        LEO2_SUCCESS, "create T32 two-block batch codec");
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(codec, kShardBytes,
        &scratch_bytes), LEO2_SUCCESS,
        "query T32 two-block batch scratch");
    const size_t scratch_region_bytes =
        (alignment + scratch_bytes + alignment + alignment - 1U) &
            ~(alignment - 1U);

    AlignedBuffer input(kBatchCount * input_region_bytes);
    AlignedBuffer output(kBatchCount * output_region_bytes);
    AlignedBuffer scratch(kBatchCount * scratch_region_bytes);
    std::memset(input.bytes(), kGuardValue, input.size());
    std::memset(output.bytes(), kGuardValue, output.size());
    std::memset(scratch.bytes(), kGuardValue, scratch.size());

    std::vector<std::vector<const void*> > original(
        kBatchCount, std::vector<const void*>(original_count));
    std::vector<std::vector<void*> > recovery(
        kBatchCount, std::vector<void*>(recovery_count));
    std::vector<leo2_encode_batch_item> items(kBatchCount);
    for (size_t item_i = 0; item_i < kBatchCount; ++item_i)
    {
        uint8_t* const input_base = input.bytes() +
            item_i * input_region_bytes + input_offset_in_region;
        uint8_t* const output_base = output.bytes() +
            item_i * output_region_bytes + output_offset_in_region;
        uint64_t state = UINT64_C(0x5433324241544348) ^ item_i;
        for (size_t byte_i = 0; byte_i < input_bytes; ++byte_i)
        {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            input_base[byte_i] = static_cast<uint8_t>(state >> 24);
        }
        for (unsigned i = 0; i < original_count; ++i)
            original[item_i][i] =
                input_base + static_cast<size_t>(i) * kShardBytes;
        for (unsigned i = 0; i < recovery_count; ++i)
            recovery[item_i][i] =
                output_base + static_cast<size_t>(i) * kShardBytes;
        items[item_i].shard_bytes = kShardBytes;
        items[item_i].original = &original[item_i][0];
        items[item_i].recovery = &recovery[item_i][0];
        items[item_i].scratch = scratch.bytes() +
            item_i * scratch_region_bytes + alignment;
        items[item_i].scratch_bytes = scratch_bytes;
    }

    const std::vector<uint8_t> input_before(
        input.bytes(), input.bytes() + input.size());
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh,
            original_count, recovery_count);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);

    RequireResult(leo2_encode_batch(codec, &items[0], items.size()),
        LEO2_SUCCESS, "execute public T32 two-block multi-item batch");
    Require(std::memcmp(input.bytes(), &input_before[0], input.size()) == 0,
        "T32 two-block multi-item batch modified source or guards");
    for (size_t item_i = 0; item_i < kBatchCount; ++item_i)
    {
        CheckRegionGuards(output, item_i * output_region_bytes,
            output_region_bytes,
            item_i * output_region_bytes + output_offset_in_region,
            output_bytes,
            "T32 two-block multi-item batch modified output guard");
        CheckRegionGuards(scratch, item_i * scratch_region_bytes,
            scratch_region_bytes,
            item_i * scratch_region_bytes + alignment,
            scratch_bytes,
            "T32 two-block multi-item batch modified scratch guard");
        CheckVariableParity(field, generator, original_count, recovery_count,
            original[item_i], recovery[item_i]);
    }

    leo2_encode_batch_binding* binding = NULL;
    RequireResult(leo2_encode_batch_binding_create(
        codec, &items[0], items.size(), &binding), LEO2_SUCCESS,
        "create packed T32 two-block batch binding");
    std::memset(output.bytes(), kGuardValue, output.size());
    RequireResult(leo2_encode_batch_binding_execute(binding),
        LEO2_SUCCESS, "execute packed T32 two-block batch binding");
    for (size_t item_i = 0; item_i < kBatchCount; ++item_i)
    {
        CheckRegionGuards(output, item_i * output_region_bytes,
            output_region_bytes,
            item_i * output_region_bytes + output_offset_in_region,
            output_bytes,
            "T32 two-block batch binding modified output guard");
        CheckRegionGuards(scratch, item_i * scratch_region_bytes,
            scratch_region_bytes,
            item_i * scratch_region_bytes + alignment,
            scratch_bytes,
            "T32 two-block batch binding modified scratch guard");
        CheckVariableParity(field, generator, original_count, recovery_count,
            original[item_i], recovery[item_i]);
    }
    leo2_encode_batch_binding_destroy(binding);

    /* One valid aliased-input stripe must fall back without disabling the
       packed terminal for its seven ordinary neighbors. */
    std::memset(output.bytes(), kGuardValue, output.size());
    const void* const saved_original = original.back()[1];
    original.back()[1] = original.back()[0];
    RequireResult(leo2_encode_batch(codec, &items[0], items.size()),
        LEO2_SUCCESS, "execute mixed packed/fallback T32 batch");
    for (size_t item_i = 0; item_i < kBatchCount; ++item_i)
    {
        CheckRegionGuards(output, item_i * output_region_bytes,
            output_region_bytes,
            item_i * output_region_bytes + output_offset_in_region,
            output_bytes,
            "mixed T32 public batch modified output guard");
        CheckRegionGuards(scratch, item_i * scratch_region_bytes,
            scratch_region_bytes,
            item_i * scratch_region_bytes + alignment,
            scratch_bytes,
            "mixed T32 public batch modified scratch guard");
        CheckVariableParity(field, generator, original_count, recovery_count,
            original[item_i], recovery[item_i]);
    }

    binding = NULL;
    RequireResult(leo2_encode_batch_binding_create(
        codec, &items[0], items.size(), &binding), LEO2_SUCCESS,
        "create mixed packed/fallback T32 batch binding");
    std::memset(output.bytes(), kGuardValue, output.size());
    RequireResult(leo2_encode_batch_binding_execute(binding),
        LEO2_SUCCESS, "execute mixed packed/fallback T32 batch binding");
    for (size_t item_i = 0; item_i < kBatchCount; ++item_i)
    {
        CheckRegionGuards(output, item_i * output_region_bytes,
            output_region_bytes,
            item_i * output_region_bytes + output_offset_in_region,
            output_bytes,
            "mixed T32 batch binding modified output guard");
        CheckRegionGuards(scratch, item_i * scratch_region_bytes,
            scratch_region_bytes,
            item_i * scratch_region_bytes + alignment,
            scratch_bytes,
            "mixed T32 batch binding modified scratch guard");
        CheckVariableParity(field, generator, original_count, recovery_count,
            original[item_i], recovery[item_i]);
    }
    leo2_encode_batch_binding_destroy(binding);
    original.back()[1] = saved_original;

    /* Every failing check must complete before an earlier worker can enter
       the terminal and write parity. */
    std::memset(output.bytes(), kGuardValue, output.size());
    const std::vector<uint8_t> output_before(
        output.bytes(), output.bytes() + output.size());
    --items.back().scratch_bytes;
    RequireResult(leo2_encode_batch(codec, &items[0], items.size()),
        LEO2_SCRATCH_TOO_SMALL,
        "reject late undersized scratch in T32 multi-item batch");
    ++items.back().scratch_bytes;
    Require(std::memcmp(output.bytes(), &output_before[0], output.size()) == 0,
        "failed T32 multi-item preflight modified an earlier output");

    void* const* const saved_recovery = items.back().recovery;
    items.back().recovery = &recovery.front()[0];
    RequireResult(leo2_encode_batch(codec, &items[0], items.size()),
        LEO2_OVERLAP,
        "reject cross-item T32 multi-item output overlap");
    items.back().recovery = saved_recovery;
    Require(std::memcmp(output.bytes(), &output_before[0], output.size()) == 0,
        "overlapping T32 multi-item batch modified output");
    Require(std::memcmp(input.bytes(), &input_before[0], input.size()) == 0,
        "T32 two-block batch or binding modified source or guards");

    leo2_codec_destroy(codec);
}
#endif

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
#if defined(LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK)
        ExerciseProductionT32TwoBlockFamily(context, 64);
        ExerciseProductionT32TwoBlockFamily(context, 96);
        ExerciseProductionT32TwoBlockFamily(context, 128);
        ExerciseProductionT32TwoBlockBatch(context);
#endif
        leo2_context_destroy(context);
#if defined(LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK)
        options.thread_count = 2;
        context = NULL;
        RequireResult(leo2_context_create(&options, &context),
            LEO2_SUCCESS, "create production pooled AVX2 context");
        ExerciseProductionT32TwoBlockBatch(context);
        leo2_context_destroy(context);
#endif
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
