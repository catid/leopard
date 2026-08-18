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
#include "Leopard2Direct.h"
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

void CheckPackedParity(
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator,
    unsigned original_count,
    unsigned recovery_count,
    const void* const* original,
    void* const* recovery)
{
    for (unsigned parity = 0; parity < recovery_count; ++parity)
    {
        const std::vector<leopard2_test::Element>& row =
            generator[original_count + parity];
        const uint8_t* const output =
            static_cast<const uint8_t*>(recovery[parity]);
        for (size_t offset = 0; offset < kShardBytes; ++offset)
        {
            leopard2_test::Element expected = 0;
            for (unsigned source = 0; source < original_count; ++source)
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

void CheckParity(
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator,
    const void* const* original,
    void* const* recovery)
{
    CheckPackedParity(field, generator, kOriginalCount, kRecoveryCount,
        original, recovery);
}

size_t ExpectedProductionWorkDataOffset(
    unsigned original_count,
    unsigned recovery_count,
    unsigned side)
{
    const size_t alignment = leo2_scratch_alignment();
    const size_t metadata_bytes =
        2U * (original_count + recovery_count) * sizeof(uintptr_t) +
        (original_count + 2U * side) * sizeof(void*);
    return (metadata_bytes + alignment - 1U) & ~(alignment - 1U);
}

size_t ExpectedProductionScratch(
    unsigned original_count,
    unsigned recovery_count,
    unsigned side)
{
    return ExpectedProductionWorkDataOffset(
        original_count, recovery_count, side) +
        2U * side * kShardBytes;
}

size_t ExpectedProductionScratch(unsigned side)
{
    return ExpectedProductionScratch(side, side, side);
}

void CheckProductionPackedScratchMetadata(
    const AlignedBuffer& scratch,
    const std::vector<uint8_t>& before,
    size_t work_data_offset,
    unsigned original_count,
    unsigned recovery_count,
    unsigned transform_side,
    bool avx512_gfni_selected,
    void* const* recovery,
    const char* message)
{
    const bool k65r65_terminal =
        original_count == 65 && recovery_count == 65 &&
        transform_side == 128 && kShardBytes == 64;
    if (!k65r65_terminal)
    {
        Require(std::memcmp(scratch.bytes(), &before[0],
                work_data_offset) == 0, message);
        return;
    }

    if (avx512_gfni_selected)
    {
        Require(std::memcmp(scratch.bytes(), &before[0],
                work_data_offset) == 0,
            "K65/R65 AVX-512/GFNI leaf modified scratch metadata");
        const size_t live_state_bytes =
            static_cast<size_t>(transform_side) * kShardBytes;
        Require(work_data_offset + live_state_bytes <= scratch.size(),
            "K65/R65 AVX-512/GFNI state exceeds public scratch");
        Require(std::memcmp(scratch.bytes() + work_data_offset +
                    live_state_bytes,
                &before[work_data_offset + live_state_bytes],
                scratch.size() - work_data_offset - live_state_bytes) == 0,
            "K65/R65 AVX-512/GFNI leaf modified unused scratch rows");
        return;
    }

    const size_t pointer_offset =
        2U * (original_count + recovery_count) * sizeof(uintptr_t);
    const size_t work_pointer_offset =
        pointer_offset + original_count * sizeof(void*);
    const size_t work_pointer_bytes =
        static_cast<size_t>(transform_side) * sizeof(void*);
    Require(work_pointer_offset + work_pointer_bytes <= work_data_offset,
        "K65/R65 staged-pointer geometry exceeds scratch metadata");
    Require(std::memcmp(scratch.bytes(), &before[0],
            work_pointer_offset) == 0,
        "K65/R65 terminal modified metadata before its work pointers");
    Require(std::memcmp(scratch.bytes() + work_pointer_offset +
                work_pointer_bytes,
            &before[work_pointer_offset + work_pointer_bytes],
            work_data_offset - work_pointer_offset - work_pointer_bytes) == 0,
        "K65/R65 terminal modified metadata after its work pointers");

    void* const* const work = reinterpret_cast<void* const*>(
        scratch.bytes() + work_pointer_offset);
    for (unsigned parity = 0; parity < recovery_count; ++parity)
    {
        Require(work[parity] == recovery[parity],
            "K65/R65 terminal staged the wrong public output pointer");
    }
    for (unsigned slot = recovery_count; slot < transform_side; ++slot)
    {
        Require(work[slot] == scratch.bytes() + work_data_offset +
                static_cast<size_t>(slot - recovery_count) * kShardBytes,
            "K65/R65 terminal staged the wrong private work pointer");
    }

    const size_t live_private_rows = transform_side - recovery_count;
    const size_t unused_data_offset =
        work_data_offset + live_private_rows * kShardBytes;
    Require(unused_data_offset <= scratch.size(),
        "K65/R65 live private rows exceed public scratch");
    Require(std::memcmp(scratch.bytes() + unused_data_offset,
            &before[unused_data_offset],
            scratch.size() - unused_data_offset) == 0,
        "K65/R65 terminal modified unused public work-data rows");
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

void ExerciseProductionPackedSide(
    leo2_context* context,
    unsigned original_count,
    unsigned side)
{
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, original_count, side,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        LEO2_SUCCESS, "create production packed codec");
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(codec, kShardBytes,
        &scratch_bytes), LEO2_SUCCESS,
        "query production balanced scratch");
    unsigned transform_side = 1;
    while (transform_side < side)
        transform_side <<= 1;
    Require(scratch_bytes ==
            ExpectedProductionScratch(
                original_count, side, transform_side),
        "packed production scratch differs from fixed geometry");
    AlignedBuffer scratch(scratch_bytes);
    const size_t work_data_offset = ExpectedProductionWorkDataOffset(
        original_count, side, transform_side);
    std::memset(scratch.bytes(), 0x5c, scratch.size());
    const std::vector<uint8_t> scratch_prefix_before(
        scratch.bytes(), scratch.bytes() + scratch.size());

    const size_t input_bytes =
        static_cast<size_t>(original_count) * kShardBytes;
    const size_t output_bytes = static_cast<size_t>(side) * kShardBytes;
    const size_t input_allocation_bytes = input_bytes > scratch_bytes
        ? input_bytes : scratch_bytes;
    const size_t output_allocation_bytes = output_bytes > scratch_bytes
        ? output_bytes : scratch_bytes;
    AlignedBuffer input(
        input_allocation_bytes + 2U * kGuardBytes + 8U);
    AlignedBuffer output(
        output_allocation_bytes + 2U * kGuardBytes + 8U);
    std::memset(input.bytes(), kGuardValue, input.size());
    std::memset(output.bytes(), kGuardValue, output.size());
    const size_t input_offset = kGuardBytes + 1U;
    const size_t output_offset = kGuardBytes + 3U;
    uint8_t* const input_base = input.bytes() + input_offset;
    uint8_t* const output_base = output.bytes() + output_offset;
    uint64_t state = UINT64_C(0x50524f4442414c36) ^
        original_count ^ (static_cast<uint64_t>(side) << 32);
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
    std::vector<void*> recovery(side);
    for (unsigned i = 0; i < original_count; ++i)
    {
        original[i] = input_base + static_cast<size_t>(i) * kShardBytes;
    }
    for (unsigned i = 0; i < side; ++i)
    {
        recovery[i] = output_base + static_cast<size_t>(i) * kShardBytes;
    }

    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh, original_count, side);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);
    const bool avx512_gfni_selected = leopard2_internal::
        K65R65B64AVX512GFNISelectedForDiagnostics(codec, kShardBytes);

    RequireResult(leo2_encode(codec, kShardBytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute production packed terminal");
    Require(std::memcmp(input.bytes(), &input_before[0], input.size()) == 0,
        "packed production encode modified source or guards");
    CheckProductionPackedScratchMetadata(scratch, scratch_prefix_before,
        work_data_offset, original_count, side, transform_side,
        avx512_gfni_selected,
        &recovery[0],
        "packed production encode staged general-path scratch metadata");
    CheckGuards(output, output_offset, output_bytes,
        "packed production encode modified output guard");

    CheckPackedParity(field, generator, original_count, side,
        &original[0], &recovery[0]);

    std::memset(output.bytes(), kGuardValue, output.size());
    leo2_encode_batch_item item = {
        kShardBytes, &original[0], &recovery[0], scratch.data(), scratch.size()
    };
    RequireResult(leo2_encode_batch(codec, &item, 1), LEO2_SUCCESS,
        "execute production packed one-item batch");
    CheckProductionPackedScratchMetadata(scratch, scratch_prefix_before,
        work_data_offset, original_count, side, transform_side,
        avx512_gfni_selected,
        &recovery[0],
        "packed production batch staged general-path scratch metadata");
    CheckGuards(output, output_offset, output_bytes,
        "packed production batch modified output guard");
    CheckPackedParity(field, generator, original_count, side,
        &original[0], &recovery[0]);

    std::memset(output.bytes(), kGuardValue, output.size());
    const std::vector<uint8_t> output_before(
        output.bytes(), output.bytes() + output.size());
    RequireResult(leo2_encode(codec, kShardBytes,
        &original[0], &recovery[0], scratch.data(), scratch.size() - 1U),
        LEO2_SCRATCH_TOO_SMALL,
        "reject undersized packed production scratch");
    RequireResult(leo2_encode(codec, kShardBytes,
        &original[0], &recovery[0], scratch.bytes() + 1U, scratch.size()),
        LEO2_BAD_ALIGNMENT,
        "reject misaligned packed production scratch");
    Require(std::memcmp(output.bytes(), &output_before[0], output.size()) == 0,
        "rejected packed production scratch modified output");

    void* const saved_recovery = recovery[0];
    recovery[0] = input_base;
    RequireResult(leo2_encode(codec, kShardBytes,
        &original[0], &recovery[0], scratch.data(), scratch.size()),
        LEO2_OVERLAP, "reject packed production source/output overlap");
    recovery[0] = saved_recovery;
    Require(std::memcmp(input.bytes(), &input_before[0], input.size()) == 0,
        "rejected packed production overlap modified input");
    Require(std::memcmp(output.bytes(), &output_before[0], output.size()) == 0,
        "rejected packed production overlap modified output");

    const bool aggregate_overlap_shape =
        (original_count == 62 && side == 8)
        || (original_count == 79 && side == 32)
        || (original_count == 65 && side == 65)
#if LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED
        || (original_count == 65 && side == 9)
        || (original_count == 66 && side == 16)
#endif
        ;
    if (aggregate_overlap_shape)
    {
        std::vector<void*> packed_overlap(side);
        for (unsigned i = 0; i < side; ++i)
        {
            packed_overlap[i] =
                input_base + static_cast<size_t>(i) * kShardBytes;
        }
        RequireResult(leo2_encode(codec, kShardBytes,
            &original[0], &packed_overlap[0], scratch.data(), scratch.size()),
            LEO2_OVERLAP,
            "reject aggregate packed source/output overlap");
        Require(std::memcmp(input.bytes(), &input_before[0], input.size()) == 0,
            "aggregate packed overlap modified input");
        Require(std::memcmp(output.bytes(), &output_before[0], output.size()) == 0,
            "aggregate packed overlap modified output");

        Require(input.size() >= scratch.size(),
            "aggregate input allocation cannot cover scratch-overlap probe");
        RequireResult(leo2_encode(codec, kShardBytes,
            &original[0], &recovery[0], input.bytes(), scratch.size()),
            LEO2_OVERLAP,
            "reject aggregate scratch/data overlap");
        Require(std::memcmp(input.bytes(), &input_before[0], input.size()) == 0,
            "aggregate scratch/data overlap modified input");
        Require(std::memcmp(output.bytes(), &output_before[0], output.size()) == 0,
            "aggregate scratch/data overlap modified output");

        Require(output.size() >= scratch.size(),
            "aggregate output allocation cannot cover scratch-overlap probe");
        RequireResult(leo2_encode(codec, kShardBytes,
            &original[0], &recovery[0], output.data(), scratch.size()),
            LEO2_OVERLAP,
            "reject aggregate scratch/output overlap");
        Require(std::memcmp(input.bytes(), &input_before[0], input.size()) == 0,
            "aggregate scratch/output overlap modified input");
        Require(std::memcmp(output.bytes(), &output_before[0], output.size()) == 0,
            "aggregate scratch/output overlap modified output");

        AlignedBuffer metadata_scratch(scratch.size());
        const void** const overlapping_original =
            reinterpret_cast<const void**>(metadata_scratch.data());
        for (unsigned i = 0; i < original_count; ++i)
            overlapping_original[i] = original[i];
        RequireResult(leo2_encode(codec, kShardBytes,
            overlapping_original, &recovery[0], metadata_scratch.data(),
            metadata_scratch.size()), LEO2_OVERLAP,
            "reject aggregate scratch/pointer-array overlap");
        Require(std::memcmp(output.bytes(), &output_before[0], output.size()) == 0,
            "aggregate pointer overlap modified output");

        AlignedBuffer batch_scratch(scratch.size());
        leo2_encode_batch_item* const overlapping_item =
            new (batch_scratch.data()) leo2_encode_batch_item;
        overlapping_item->shard_bytes = kShardBytes;
        overlapping_item->original = &original[0];
        overlapping_item->recovery = &recovery[0];
        overlapping_item->scratch = batch_scratch.data();
        overlapping_item->scratch_bytes = batch_scratch.size();
        RequireResult(leo2_encode_batch(codec, overlapping_item, 1),
            LEO2_OVERLAP,
            "reject aggregate scratch/batch-descriptor overlap");
        Require(std::memcmp(output.bytes(), &output_before[0], output.size()) == 0,
            "aggregate batch overlap modified output");
    }
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
        ExerciseProductionPackedSide(context, 16, 16);
        ExerciseProductionPackedSide(context, 62, 8);
        // Exercise both three- and four-block T16 endpoints in the ordinary
        // production archive.  The diagnostic test performs the bounded
        // all-K sweep; these cases retain guard, scratch, batch, overlap, and
        // direct-generator checks without duplicating that full matrix here.
        ExerciseProductionPackedSide(context, 33, 9);
        ExerciseProductionPackedSide(context, 64, 9);
#if LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED
        ExerciseProductionPackedSide(context, 65, 9);
#endif
        ExerciseProductionPackedSide(context, 33, 16);
        ExerciseProductionPackedSide(context, 62, 16);
        ExerciseProductionPackedSide(context, 64, 16);
#if LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED
        ExerciseProductionPackedSide(context, 65, 16);
        ExerciseProductionPackedSide(context, 66, 16);
#endif
        ExerciseProductionPackedSide(context, 33, 17);
        ExerciseProductionPackedSide(context, 62, 17);
        ExerciseProductionPackedSide(context, 62, 30);
        ExerciseProductionPackedSide(context, 33, 31);
        ExerciseProductionPackedSide(context, 62, 31);
        ExerciseProductionPackedSide(context, 33, 32);
        ExerciseProductionPackedSide(context, 34, 32);
        ExerciseProductionPackedSide(context, 62, 32);
        ExerciseProductionPackedSide(context, 63, 32);
        ExerciseProductionPackedSide(context, 64, 32);
        ExerciseProductionPackedSide(context, 65, 32);
        ExerciseProductionPackedSide(context, 79, 32);
        ExerciseProductionPackedSide(context, 96, 32);
        ExerciseProductionPackedSide(context, 33, 33);
        ExerciseProductionPackedSide(context, 34, 33);
        ExerciseProductionPackedSide(context, 62, 33);
        ExerciseProductionPackedSide(context, 64, 33);
        ExerciseProductionPackedSide(context, 33, 34);
        ExerciseProductionPackedSide(context, 62, 34);
        ExerciseProductionPackedSide(context, 64, 34);
        ExerciseProductionPackedSide(context, 62, 35);
        ExerciseProductionPackedSide(context, 33, 64);
        ExerciseProductionPackedSide(context, 62, 64);
        ExerciseProductionPackedSide(context, 64, 64);
        ExerciseProductionPackedSide(context, 65, 65);
        ExerciseProductionPackedSide(context, 128, 128);
#if defined(LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK)
        ExerciseProductionT32TwoBlockFamily(context, 64);
        ExerciseProductionT32TwoBlockFamily(context, 96);
        ExerciseProductionT32TwoBlockFamily(context, 128);
        ExerciseProductionT32TwoBlockBatch(context);
#endif
        leo2_context_destroy(context);

        if (kShardBytes == 64)
        {
            // The operation-specific raised-ISA route is deliberately
            // AUTO-only and byte-exact.  Recreate only its B=64 production
            // geometry; this source also builds the independent B=256 suite.
            options.backend = LEO2_BACKEND_AUTO;
            options.thread_count = 1;
            context = NULL;
            RequireResult(leo2_context_create(&options, &context),
                LEO2_SUCCESS, "create production AUTO T128 context");
            ExerciseProductionPackedSide(context, 65, 65);
            leo2_context_destroy(context);
        }
#if defined(LEO2_EXPERIMENT_HIGH_T32_B256_TWO_BLOCK)
        options.backend = LEO2_BACKEND_AVX2;
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
