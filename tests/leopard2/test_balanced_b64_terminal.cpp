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
#include "LeopardFF8.h"
#include "leopard.h"
#include "leopard2.h"

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
static const size_t kShardBytes = 64;

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

void FillInput(uint8_t* input, size_t shard_bytes)
{
    uint64_t state = UINT64_C(0x4b33325233324236) ^ shard_bytes;
    for (size_t i = 0; i < kOriginalCount * shard_bytes; ++i)
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
    size_t shard_bytes,
    const void** original,
    void** recovery)
{
    for (unsigned i = 0; i < kOriginalCount; ++i)
        original[i] = input + static_cast<size_t>(i) * shard_bytes;
    for (unsigned i = 0; i < kRecoveryCount; ++i)
        recovery[i] = output + static_cast<size_t>(i) * shard_bytes;
}

void CheckParity(
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator,
    const void* const* original,
    void* const* recovery,
    size_t shard_bytes)
{
    for (unsigned parity = 0; parity < kRecoveryCount; ++parity)
    {
        if (!recovery[parity])
            continue;
        const std::vector<leopard2_test::Element>& row =
            generator[kOriginalCount + parity];
        const uint8_t* const output =
            static_cast<const uint8_t*>(recovery[parity]);
        for (size_t offset = 0; offset < shard_bytes; ++offset)
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
                "K32/R32 parity differs from the independent oracle");
        }
    }
}

void RequireTerminalCalls(uint64_t expected, const char* message)
{
    const leopard::ff8::TestOnlyHighEncodeCounts counts =
        leopard::ff8::TestOnlyGetHighEncodeCounts();
    Require(counts.balanced_b64_packed_calls == expected, message);
}

bool CreateOptionalGF16Codec(
    leo2_context* context,
    unsigned original_count,
    unsigned recovery_count,
    leo2_codec** codec_out,
    const char* message)
{
    const bool supported =
        (leo2_context_field_mask(context) & LEO2_FIELD_MASK_GF16) != 0;
    const leo2_result result = leo2_codec_create(
        context, original_count, recovery_count,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, NULL, codec_out);
    RequireResult(result,
        supported ? LEO2_SUCCESS : LEO2_UNSUPPORTED, message);
    Require((*codec_out != NULL) == supported,
        "GF16 codec result disagrees with the context field mask");
    return supported;
}

leo2_codec* CreateCodec(leo2_context* context, leo2_profile profile)
{
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context,
        kOriginalCount, kRecoveryCount, profile, LEO2_FIELD_GF8,
        NULL, &codec), LEO2_SUCCESS, "create K32/R32 codec");
    return codec;
}

void FillVariableInput(uint8_t* input, unsigned original_count,
    size_t shard_bytes, uint64_t seed)
{
    uint64_t state = seed ^ original_count ^
        (static_cast<uint64_t>(shard_bytes) << 32);
    for (size_t i = 0;
         i < static_cast<size_t>(original_count) * shard_bytes; ++i)
    {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        input[i] = static_cast<uint8_t>(state >> 24);
    }
}

void SetVariablePackedPointers(
    uint8_t* input,
    uint8_t* output,
    unsigned original_count,
    unsigned recovery_count,
    size_t shard_bytes,
    std::vector<const void*>& original,
    std::vector<void*>& recovery)
{
    original.resize(original_count);
    recovery.resize(recovery_count);
    for (unsigned i = 0; i < original_count; ++i)
        original[i] = input + static_cast<size_t>(i) * shard_bytes;
    for (unsigned i = 0; i < recovery_count; ++i)
        recovery[i] = output + static_cast<size_t>(i) * shard_bytes;
}

void CheckVariableParity(
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator,
    unsigned original_count,
    unsigned recovery_count,
    const std::vector<const void*>& original,
    const std::vector<void*>& recovery,
    size_t shard_bytes,
    const char* message)
{
    for (unsigned parity = 0; parity < recovery_count; ++parity)
    {
        if (!recovery[parity])
            continue;
        const std::vector<leopard2_test::Element>& row =
            generator[original_count + parity];
        const uint8_t* const output =
            static_cast<const uint8_t*>(recovery[parity]);
        for (size_t offset = 0; offset < shard_bytes; ++offset)
        {
            leopard2_test::Element expected = 0;
            for (unsigned source = 0; source < original_count; ++source)
            {
                const uint8_t value = static_cast<const uint8_t*>(
                    original[source])[offset];
                expected = field.add(expected,
                    field.multiply(row[source], value));
            }
            Require(output[offset] == expected, message);
        }
    }
}

void CheckVariableParityGF16(
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator,
    unsigned original_count,
    unsigned recovery_count,
    const std::vector<const void*>& original,
    const std::vector<void*>& recovery,
    size_t shard_bytes,
    const char* message);

void CheckLegacyParity(
    unsigned original_count,
    unsigned recovery_count,
    const std::vector<const void*>& original,
    const std::vector<void*>& recovery,
    size_t shard_bytes,
    const char* message);

void CheckVariableBasisParity(
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator,
    unsigned original_count,
    unsigned recovery_count,
    unsigned live_source,
    const std::vector<const void*>& original,
    const std::vector<void*>& recovery,
    size_t shard_bytes,
    const char* message)
{
    Require(live_source < original_count,
        "basis source lies outside the systematic dimension");
    const uint8_t* const source =
        static_cast<const uint8_t*>(original[live_source]);
    for (unsigned parity = 0; parity < recovery_count; ++parity)
    {
        const leopard2_test::Element coefficient =
            generator[original_count + parity][live_source];
        const uint8_t* const output =
            static_cast<const uint8_t*>(recovery[parity]);
        for (size_t offset = 0; offset < shard_bytes; ++offset)
        {
            Require(output[offset] ==
                    field.multiply(coefficient, source[offset]),
                message);
        }
    }
}

void ExerciseK65R65PackedTerminal(leo2_context* context)
{
    static const unsigned original_count = 65;
    static const unsigned recovery_count = 65;
    static const size_t shard_bytes = 64;

    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, original_count, recovery_count,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        LEO2_SUCCESS, "create K65/R65/B64 packed-terminal codec");
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(codec, shard_bytes,
        &scratch_bytes), LEO2_SUCCESS,
        "query K65/R65/B64 packed-terminal scratch");
    AlignedBuffer scratch(scratch_bytes);
    AlignedBuffer second_scratch(scratch_bytes);
    const size_t input_bytes =
        static_cast<size_t>(original_count) * shard_bytes;
    const size_t output_bytes =
        static_cast<size_t>(recovery_count) * shard_bytes;
    const size_t data_allocation_bytes =
        scratch_bytes > input_bytes ? scratch_bytes : input_bytes;
    AlignedBuffer input(data_allocation_bytes + 8U);
    AlignedBuffer output(data_allocation_bytes + 8U);
    AlignedBuffer second_input(input_bytes);
    AlignedBuffer second_output(output_bytes);
    uint8_t* const input_base = input.bytes() + 1U;
    uint8_t* const output_base = output.bytes() + 3U;
    FillVariableInput(input_base, original_count, shard_bytes,
        UINT64_C(0x4b36355236354236));
    FillVariableInput(second_input.bytes(), original_count, shard_bytes,
        UINT64_C(0x4b36355236354232));
    std::vector<const void*> original;
    std::vector<void*> recovery;
    std::vector<const void*> second_original;
    std::vector<void*> second_recovery;
    SetVariablePackedPointers(input_base, output_base, original_count,
        recovery_count, shard_bytes, original, recovery);
    SetVariablePackedPointers(second_input.bytes(), second_output.bytes(),
        original_count, recovery_count, shard_bytes,
        second_original, second_recovery);

    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(leopard2_test::kLegacyHigh,
            original_count, recovery_count);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);

    std::memset(output_base, 0xa5, output_bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute unaligned K65/R65/B64 packed terminal");
    RequireTerminalCalls(1,
        "K65/R65/B64 one-shot missed the exact packed terminal");
    CheckVariableParity(field, generator, original_count, recovery_count,
        original, recovery, shard_bytes,
        "K65/R65/B64 one-shot parity differs from direct oracle");
    CheckLegacyParity(original_count, recovery_count, original, recovery,
        shard_bytes,
        "K65/R65/B64 one-shot parity differs from legacy leo_encode");

    std::memset(output_base, 0xa5, output_bytes);
    leo2_encode_batch_item item = {
        shard_bytes, &original[0], &recovery[0], scratch.data(), scratch.size()
    };
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch(codec, &item, 1), LEO2_SUCCESS,
        "execute one-item K65/R65/B64 packed batch");
    RequireTerminalCalls(1,
        "K65/R65/B64 one-item batch missed the exact packed terminal");
    CheckVariableParity(field, generator, original_count, recovery_count,
        original, recovery, shard_bytes,
        "K65/R65/B64 one-item batch parity differs from direct oracle");

    /* Exercise every shortened systematic column independently.  Keeping
       source 64 last makes the Leopard1 comparison below an explicit guard
       for the first coordinate beyond the complete 64-shard half block. */
    for (unsigned source = 0; source < original_count; ++source)
    {
        std::memset(input_base, 0, input_bytes);
        for (size_t offset = 0; offset < shard_bytes; ++offset)
        {
            input_base[static_cast<size_t>(source) * shard_bytes + offset] =
                static_cast<uint8_t>(offset + 1U);
        }
        std::memset(output_base, 0xa5, output_bytes);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode(codec, shard_bytes,
            &original[0], &recovery[0], scratch.data(), scratch.size()),
            LEO2_SUCCESS, "execute K65/R65 systematic basis encode");
        RequireTerminalCalls(1,
            "K65/R65 systematic basis encode missed the packed terminal");
        CheckVariableBasisParity(field, generator, original_count,
            recovery_count, source, original, recovery, shard_bytes,
            "K65/R65 systematic basis parity differs from direct oracle");
    }
    CheckLegacyParity(original_count, recovery_count, original, recovery,
        shard_bytes,
        "K65/R65 source-64 basis parity differs from legacy leo_encode");
    Require(static_cast<const uint8_t*>(recovery[64])[0] != 0,
        "K65/R65 source-64 contribution unexpectedly vanished");

    FillVariableInput(input_base, original_count, shard_bytes,
        UINT64_C(0x4b36355236354236));
    Require(leopard2_internal::
            SetK65R65B64PackedTerminalEnabledForDiagnostics(false),
        "disable exact K65/R65/B64 packed terminal");
    std::memset(output_base, 0xa5, output_bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute same-binary K65/R65/B64 mature control");
    RequireTerminalCalls(0,
        "disabled K65/R65/B64 terminal was still selected");
    CheckVariableParity(field, generator, original_count, recovery_count,
        original, recovery, shard_bytes,
        "disabled K65/R65/B64 terminal changed parity");
    Require(leopard2_internal::
            SetK65R65B64PackedTerminalEnabledForDiagnostics(true),
        "restore exact K65/R65/B64 packed terminal");

    RequireResult(leo2_test_codec_set_encode_mode(
        codec, LEO2_TEST_ENCODE_FORCE_TRANSFORM), LEO2_SUCCESS,
        "force mature K65/R65/B64 transform");
    std::memset(output_base, 0xa5, output_bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute forced mature K65/R65/B64 transform");
    RequireTerminalCalls(0,
        "forced K65/R65 transform entered the packed terminal");
    CheckVariableParity(field, generator, original_count, recovery_count,
        original, recovery, shard_bytes,
        "forced mature K65/R65 parity differs from direct oracle");
    RequireResult(leo2_test_codec_set_encode_mode(
        codec, LEO2_TEST_ENCODE_AUTO), LEO2_SUCCESS,
        "restore K65/R65/B64 AUTO encode mode");

    std::vector<uint8_t> detached_source(shard_bytes);
    std::memcpy(&detached_source[0], original[64], shard_bytes);
    original[64] = &detached_source[0];
    std::memset(output_base, 0xa5, output_bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute detached-source K65/R65 fallback");
    RequireTerminalCalls(0,
        "detached K65/R65 source entered the packed terminal");
    CheckVariableParity(field, generator, original_count, recovery_count,
        original, recovery, shard_bytes,
        "detached-source K65/R65 parity differs from direct oracle");
    original[64] = input_base + 64U * shard_bytes;

    const void* const saved_original3 = original[3];
    original[3] = original[51];
    original[51] = saved_original3;
    std::memset(output_base, 0xa5, output_bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute reordered-source K65/R65 fallback");
    RequireTerminalCalls(0,
        "reordered K65/R65 sources entered the packed terminal");
    CheckVariableParity(field, generator, original_count, recovery_count,
        original, recovery, shard_bytes,
        "reordered-source K65/R65 parity differs from direct oracle");
    original[51] = original[3];
    original[3] = saved_original3;

    std::vector<uint8_t> detached_output(shard_bytes);
    void* const saved_recovery27 = recovery[27];
    recovery[27] = &detached_output[0];
    std::memset(output_base, 0xa5, output_bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute detached-output K65/R65 fallback");
    RequireTerminalCalls(0,
        "detached K65/R65 output entered the packed terminal");
    CheckVariableParity(field, generator, original_count, recovery_count,
        original, recovery, shard_bytes,
        "detached-output K65/R65 parity differs from direct oracle");
    recovery[27] = saved_recovery27;

    void* const saved_recovery17 = recovery[17];
    recovery[17] = NULL;
    std::memset(output_base, 0xa5, output_bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute sparse-output K65/R65 fallback");
    RequireTerminalCalls(0,
        "sparse-output K65/R65 encode entered the packed terminal");
    CheckVariableParity(field, generator, original_count, recovery_count,
        original, recovery, shard_bytes,
        "sparse-output K65/R65 parity differs from direct oracle");
    for (size_t offset = 0; offset < shard_bytes; ++offset)
    {
        Require(output_base[17U * shard_bytes + offset] == 0xa5,
            "sparse K65/R65 encode modified the omitted output");
    }
    recovery[17] = saved_recovery17;

    std::memset(output_base, 0xa5, output_bytes);
    const std::vector<uint8_t> invalid_output_before(
        output_base, output_base + output_bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, NULL, &recovery[0],
        scratch.data(), scratch.size()), LEO2_INVALID_ARGUMENT,
        "reject null K65/R65 original array");
    RequireTerminalCalls(0,
        "null K65/R65 original array reached the packed terminal");
    RequireResult(leo2_encode(codec, shard_bytes, &original[0], NULL,
        scratch.data(), scratch.size()), LEO2_INVALID_ARGUMENT,
        "reject null K65/R65 recovery array");
    RequireTerminalCalls(0,
        "null K65/R65 recovery array reached the packed terminal");
    const void* const saved_original13 = original[13];
    original[13] = NULL;
    RequireResult(leo2_encode(codec, shard_bytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()), LEO2_INVALID_ARGUMENT,
        "reject interior null K65/R65 source");
    RequireTerminalCalls(0,
        "interior null K65/R65 source reached the packed terminal");
    original[13] = saved_original13;
    Require(std::memcmp(output_base, &invalid_output_before[0],
            output_bytes) == 0,
        "invalid K65/R65 pointer input modified output");

    void* const saved_recovery0 = recovery[0];
    recovery[0] = NULL;
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute null-recovery-zero K65/R65 fallback");
    RequireTerminalCalls(0,
        "null recovery[0] entered the K65/R65 packed terminal");
    CheckVariableParity(field, generator, original_count, recovery_count,
        original, recovery, shard_bytes,
        "null-recovery-zero K65/R65 parity differs from direct oracle");
    for (size_t offset = 0; offset < shard_bytes; ++offset)
    {
        Require(output_base[offset] == 0xa5,
            "null recovery[0] fallback modified the omitted output");
    }
    recovery[0] = saved_recovery0;

    const size_t byte_neighbors[] = { 63U, 65U };
    for (size_t neighbor_i = 0;
         neighbor_i < sizeof(byte_neighbors) / sizeof(byte_neighbors[0]);
         ++neighbor_i)
    {
        const size_t neighbor_bytes = byte_neighbors[neighbor_i];
        size_t neighbor_scratch_bytes = 0;
        RequireResult(leo2_encode_scratch_size(codec, neighbor_bytes,
            &neighbor_scratch_bytes), LEO2_SUCCESS,
            "query K65/R65 neighboring-byte scratch");
        AlignedBuffer neighbor_scratch(neighbor_scratch_bytes);
        AlignedBuffer neighbor_input(
            static_cast<size_t>(original_count) * neighbor_bytes);
        AlignedBuffer neighbor_output(
            static_cast<size_t>(recovery_count) * neighbor_bytes);
        FillVariableInput(neighbor_input.bytes(), original_count,
            neighbor_bytes, UINT64_C(0x4b36355236354e45) + neighbor_i);
        std::vector<const void*> neighbor_original;
        std::vector<void*> neighbor_recovery;
        SetVariablePackedPointers(neighbor_input.bytes(),
            neighbor_output.bytes(), original_count, recovery_count,
            neighbor_bytes, neighbor_original, neighbor_recovery);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode(codec, neighbor_bytes,
            &neighbor_original[0], &neighbor_recovery[0],
            neighbor_scratch.data(), neighbor_scratch.size()),
            LEO2_SUCCESS, "execute K65/R65 neighboring-byte mature path");
        RequireTerminalCalls(0,
            "neighboring byte count entered the K65/R65/B64 terminal");
        CheckVariableParity(field, generator, original_count, recovery_count,
            neighbor_original, neighbor_recovery, neighbor_bytes,
            "K65/R65 neighboring-byte parity differs from direct oracle");
    }

    std::memset(output_base, 0xa5, output_bytes);
    const std::vector<uint8_t> output_before(
        output_base, output_base + output_bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, &original[0], &recovery[0],
        scratch.data(), scratch.size() - 1U), LEO2_SCRATCH_TOO_SMALL,
        "reject undersized K65/R65 packed-terminal scratch");
    RequireTerminalCalls(0,
        "undersized K65/R65 scratch reached the packed terminal");
    Require(std::memcmp(output_base, &output_before[0], output_bytes) == 0,
        "undersized K65/R65 scratch modified output");
    RequireResult(leo2_encode(codec, shard_bytes, &original[0], &recovery[0],
        scratch.bytes() + 1U, scratch.size()), LEO2_BAD_ALIGNMENT,
        "reject misaligned K65/R65 packed-terminal scratch");
    RequireTerminalCalls(0,
        "misaligned K65/R65 scratch reached the packed terminal");
    Require(std::memcmp(output_base, &output_before[0], output_bytes) == 0,
        "misaligned K65/R65 scratch modified output");

    std::vector<void*> overlapping_recovery(recovery_count);
    for (unsigned parity = 0; parity < recovery_count; ++parity)
    {
        overlapping_recovery[parity] =
            input_base + static_cast<size_t>(parity) * shard_bytes;
    }
    const std::vector<uint8_t> input_before(input_base,
        input_base + input_bytes);
    RequireResult(leo2_encode(codec, shard_bytes, &original[0],
        &overlapping_recovery[0], scratch.data(), scratch.size()),
        LEO2_OVERLAP, "reject aggregate K65/R65 source/output overlap");
    Require(std::memcmp(input_base, &input_before[0], input_bytes) == 0,
        "K65/R65 overlap rejection modified input");
    Require(std::memcmp(output_base, &output_before[0], output_bytes) == 0,
        "K65/R65 overlap rejection modified output");

    RequireResult(leo2_encode(codec, shard_bytes, &original[0], &recovery[0],
        input.data(), scratch.size()), LEO2_OVERLAP,
        "reject K65/R65 scratch/source overlap");
    Require(std::memcmp(input_base, &input_before[0], input_bytes) == 0,
        "K65/R65 scratch/source rejection modified input");
    RequireResult(leo2_encode(codec, shard_bytes, &original[0], &recovery[0],
        output.data(), scratch.size()), LEO2_OVERLAP,
        "reject K65/R65 scratch/output overlap");
    Require(std::memcmp(output_base, &output_before[0], output_bytes) == 0,
        "K65/R65 scratch/output rejection modified output");

    AlignedBuffer metadata_scratch(scratch.size());
    const void** const overlapping_original =
        reinterpret_cast<const void**>(metadata_scratch.data());
    for (unsigned source = 0; source < original_count; ++source)
        overlapping_original[source] = original[source];
    RequireResult(leo2_encode(codec, shard_bytes, overlapping_original,
        &recovery[0], metadata_scratch.data(), metadata_scratch.size()),
        LEO2_OVERLAP, "reject K65/R65 scratch/pointer-array overlap");
    Require(std::memcmp(output_base, &output_before[0], output_bytes) == 0,
        "K65/R65 pointer-array rejection modified output");

    AlignedBuffer protected_storage(output_bytes);
    leo2_encode_batch_item* const protected_item =
        new (protected_storage.data()) leo2_encode_batch_item;
    protected_item->shard_bytes = shard_bytes;
    protected_item->original = &original[0];
    protected_item->recovery = &recovery[0];
    protected_item->scratch = scratch.data();
    protected_item->scratch_bytes = scratch.size();
    for (unsigned parity = 0; parity < recovery_count; ++parity)
    {
        recovery[parity] = protected_storage.bytes() +
            static_cast<size_t>(parity) * shard_bytes;
    }
    const leo2_encode_batch_item protected_before = *protected_item;
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch(codec, protected_item, 1),
        LEO2_OVERLAP, "reject K65/R65 output/batch-metadata overlap");
    RequireTerminalCalls(0,
        "K65/R65 batch-metadata overlap reached the packed terminal");
    Require(std::memcmp(protected_item, &protected_before,
            sizeof(*protected_item)) == 0,
        "K65/R65 metadata rejection modified the batch descriptor");
    SetVariablePackedPointers(input_base, output_base, original_count,
        recovery_count, shard_bytes, original, recovery);

    leo2_encode_batch_item two_items[2] = {
        { shard_bytes, &original[0], &recovery[0],
            scratch.data(), scratch.size() },
        { shard_bytes, &second_original[0], &second_recovery[0],
            second_scratch.data(), second_scratch.size() }
    };
    std::memset(output_base, 0xa5, output_bytes);
    std::memset(second_output.bytes(), 0xa5, output_bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch(codec, two_items, 2), LEO2_SUCCESS,
        "execute two-item K65/R65 mature batch");
    RequireTerminalCalls(0,
        "two-item K65/R65 batch entered the one-item packed terminal");
    CheckVariableParity(field, generator, original_count, recovery_count,
        original, recovery, shard_bytes,
        "first two-item K65/R65 parity differs from direct oracle");
    CheckVariableParity(field, generator, original_count, recovery_count,
        second_original, second_recovery, shard_bytes,
        "second two-item K65/R65 parity differs from direct oracle");

    item.original = &original[0];
    item.recovery = &recovery[0];
    leo2_encode_batch_binding* binding = NULL;
    RequireResult(leo2_encode_batch_binding_create(codec, &item, 1, &binding),
        LEO2_SUCCESS, "create reusable K65/R65 batch binding");
    std::memset(output_base, 0xa5, output_bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch_binding_execute(binding), LEO2_SUCCESS,
        "execute reusable K65/R65 batch binding");
    RequireTerminalCalls(0,
        "reusable K65/R65 binding entered the one-shot packed terminal");
    CheckVariableParity(field, generator, original_count, recovery_count,
        original, recovery, shard_bytes,
        "reusable K65/R65 binding parity differs from direct oracle");
    leo2_encode_batch_binding_destroy(binding);

    leo2_codec_destroy(codec);

    leo2_codec* low_codec = NULL;
    RequireResult(leo2_codec_create(context, original_count, recovery_count,
        LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8, NULL, &low_codec),
        LEO2_SUCCESS, "create low-profile K65/R65 exclusion codec");
    size_t low_scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(low_codec, shard_bytes,
        &low_scratch_bytes), LEO2_SUCCESS,
        "query low-profile K65/R65 scratch");
    AlignedBuffer low_scratch(low_scratch_bytes);
    const leopard2_test::ProfileLayout low_layout =
        leopard2_test::make_profile_layout(leopard2_test::kLow,
            original_count, recovery_count);
    const leopard2_test::Matrix low_generator =
        leopard2_test::direct_systematic_generator(field, low_layout);
    FillVariableInput(input_base, original_count, shard_bytes,
        UINT64_C(0x4b36354c4f575052));
    std::memset(output_base, 0xa5, output_bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(low_codec, shard_bytes,
        &original[0], &recovery[0], low_scratch.data(), low_scratch.size()),
        LEO2_SUCCESS, "execute low-profile K65/R65 exclusion");
    RequireTerminalCalls(0,
        "low-profile K65/R65 entered the legacy-high packed terminal");
    CheckVariableParity(field, low_generator, original_count, recovery_count,
        original, recovery, shard_bytes,
        "low-profile K65/R65 parity differs from direct oracle");
    leo2_codec_destroy(low_codec);

    leo2_codec* gf16_codec = NULL;
    if (CreateOptionalGF16Codec(context, original_count, recovery_count,
            &gf16_codec, "create GF16 K65/R65 exclusion codec"))
    {
        size_t gf16_scratch_bytes = 0;
        RequireResult(leo2_encode_scratch_size(gf16_codec, shard_bytes,
            &gf16_scratch_bytes), LEO2_SUCCESS,
            "query GF16 K65/R65 exclusion scratch");
        AlignedBuffer gf16_scratch(gf16_scratch_bytes);
        const leopard2_test::BinaryField gf16 =
            leopard2_test::make_legacy_gf16();
        const leopard2_test::Matrix gf16_generator =
            leopard2_test::direct_systematic_generator(gf16, layout);
        FillVariableInput(input_base, original_count, shard_bytes,
            UINT64_C(0x4b36354746313650));
        std::memset(output_base, 0xa5, output_bytes);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode(gf16_codec, shard_bytes,
            &original[0], &recovery[0], gf16_scratch.data(),
            gf16_scratch.size()), LEO2_SUCCESS,
            "execute GF16 K65/R65 exclusion");
        RequireTerminalCalls(0,
            "GF16 K65/R65 entered the GF8 packed terminal");
        CheckVariableParityGF16(gf16, gf16_generator, original_count,
            recovery_count, original, recovery, shard_bytes,
            "GF16 K65/R65 parity differs from direct oracle");
        leo2_codec_destroy(gf16_codec);
    }
}

leopard2_test::Element ReadNativeGF16Symbol(
    const void* shard,
    size_t shard_bytes,
    size_t symbol)
{
    const uint8_t* const bytes = static_cast<const uint8_t*>(shard);
    const size_t total_symbols = shard_bytes / 2U;
    const size_t tile_symbol = symbol & 31U;
    const size_t tile_first = symbol - tile_symbol;
    const size_t tile_symbols =
        total_symbols - tile_first < 32U
            ? total_symbols - tile_first : 32U;
    const size_t tile_byte = tile_first * 2U;
    return static_cast<leopard2_test::Element>(
        bytes[tile_byte + tile_symbol] |
        (static_cast<unsigned>(
            bytes[tile_byte + tile_symbols + tile_symbol]) << 8));
}

void CheckVariableParityGF16(
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator,
    unsigned original_count,
    unsigned recovery_count,
    const std::vector<const void*>& original,
    const std::vector<void*>& recovery,
    size_t shard_bytes,
    const char* message)
{
    Require((shard_bytes & 1U) == 0,
        "native GF16 parity check requires an even byte count");
    const size_t symbols = shard_bytes / 2U;
    for (unsigned parity = 0; parity < recovery_count; ++parity)
    {
        const std::vector<leopard2_test::Element>& row =
            generator[original_count + parity];
        for (size_t symbol = 0; symbol < symbols; ++symbol)
        {
            leopard2_test::Element expected = 0;
            for (unsigned source = 0; source < original_count; ++source)
            {
                expected = field.add(expected, field.multiply(row[source],
                    ReadNativeGF16Symbol(
                        original[source], shard_bytes, symbol)));
            }
            Require(ReadNativeGF16Symbol(
                    recovery[parity], shard_bytes, symbol) == expected,
                message);
        }
    }
}
void CheckLegacyParity(
    unsigned original_count,
    unsigned recovery_count,
    const std::vector<const void*>& original,
    const std::vector<void*>& recovery,
    size_t shard_bytes,
    const char* message)
{
    const unsigned work_count =
        leo_encode_work_count(original_count, recovery_count);
    Require(work_count >= recovery_count,
        "legacy work-count query rejected selected packed shape");
    AlignedBuffer legacy_storage(
        static_cast<size_t>(work_count) * shard_bytes);
    std::vector<void*> legacy_work(work_count);
    for (unsigned i = 0; i < work_count; ++i)
    {
        legacy_work[i] = legacy_storage.bytes() +
            static_cast<size_t>(i) * shard_bytes;
    }
    Require(leo_encode(shard_bytes, original_count, recovery_count,
            work_count, &original[0], &legacy_work[0]) == Leopard_Success,
        "legacy selected packed encode failed");
    for (unsigned parity = 0; parity < recovery_count; ++parity)
    {
        Require(std::memcmp(recovery[parity], legacy_work[parity],
                shard_bytes) == 0,
            message);
    }
}

void ExerciseSelectedPackedShape(
    leo2_context* context,
    unsigned original_count,
    unsigned side)
{
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, original_count, side,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        LEO2_SUCCESS, "create selected packed codec");
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, kShardBytes, &scratch_bytes), LEO2_SUCCESS,
        "query selected balanced scratch");
    AlignedBuffer scratch(scratch_bytes);
    AlignedBuffer input(
        static_cast<size_t>(original_count) * kShardBytes);
    AlignedBuffer output(static_cast<size_t>(side) * kShardBytes);
    FillVariableInput(input.bytes(), original_count, kShardBytes,
        UINT64_C(0x42414c414e434536));
    std::vector<const void*> original;
    std::vector<void*> recovery;
    SetVariablePackedPointers(input.bytes(), output.bytes(), original_count,
        side, kShardBytes, original, recovery);

    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh, original_count, side);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);

    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute selected balanced terminal");
    const bool k62_fused_shape = original_count == 62 && side == 8;
    RequireTerminalCalls(k62_fused_shape ? 0U : 1U,
        "selected shape used the wrong public terminal route");
    const leopard::ff8::TestOnlyHighEncodeCounts one_shot_counts =
        leopard::ff8::TestOnlyGetHighEncodeCounts();
    Require(one_shot_counts.t8_k62_b64_fused_calls ==
            (k62_fused_shape ? 1U : 0U),
        "selected one-shot K62/R8 fused-leaf attribution is wrong");
    if (side >= 9 && side <= 16)
    {
        Require(one_shot_counts.tail_column_calls ==
                (original_count > 64 ? 1U : 0U),
            "selected T16 one-shot tail attribution is wrong");
    }
    CheckVariableParity(field, generator, original_count, side,
        original, recovery, kShardBytes,
        "selected packed parity differs from direct oracle");
    const bool legacy_boundary =
        (original_count == 62 && side == 8) ||
        (side == 32 && original_count >= 65 && original_count <= 96) ||
        ((side == 9 || side == 16) &&
         (original_count == 33 || original_count == 48 ||
          original_count == 49 || original_count == 64 ||
          original_count == 65 || original_count == 66));
    if (legacy_boundary)
    {
        CheckLegacyParity(original_count, side, original, recovery,
            kShardBytes,
            "selected packed parity differs from legacy leo_encode");
    }

    std::memset(output.bytes(), 0xa5,
        static_cast<size_t>(side) * kShardBytes);
    leo2_encode_batch_item item = {
        kShardBytes, &original[0], &recovery[0], scratch.data(), scratch.size()
    };
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch(codec, &item, 1), LEO2_SUCCESS,
        "execute selected balanced one-item batch terminal");
    RequireTerminalCalls(k62_fused_shape ? 0U : 1U,
        "selected batch used the wrong public terminal route");
    const leopard::ff8::TestOnlyHighEncodeCounts batch_counts =
        leopard::ff8::TestOnlyGetHighEncodeCounts();
    Require(batch_counts.t8_k62_b64_fused_calls ==
            (k62_fused_shape ? 1U : 0U),
        "selected batch K62/R8 fused-leaf attribution is wrong");
    if (side >= 9 && side <= 16)
    {
        Require(batch_counts.tail_column_calls ==
                (original_count > 64 ? 1U : 0U),
            "selected T16 batch tail attribution is wrong");
    }
    CheckVariableParity(field, generator, original_count, side,
        original, recovery, kShardBytes,
        "selected packed batch parity differs from oracle");

    if (side == 32 && original_count >= 65 && original_count <= 96)
    {
        /* Each K changes the last live group in the third T=32 message
           block.  Isolate that final systematic coordinate so shortened
           zero handling cannot be hidden by a favorable dense random
           message. */
        std::memset(input.bytes(), 0,
            static_cast<size_t>(original_count) * kShardBytes);
        std::memset(input.bytes() +
                static_cast<size_t>(original_count - 1U) * kShardBytes,
            1, kShardBytes);
        std::memset(output.bytes(), 0xa5,
            static_cast<size_t>(side) * kShardBytes);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode(codec, kShardBytes,
            &original[0], &recovery[0], scratch.data(), scratch.size()),
            LEO2_SUCCESS, "execute T32/Q3 final-source basis encode");
        RequireTerminalCalls(1,
            "T32/Q3 final-source basis encode missed the packed terminal");
        CheckVariableParity(field, generator, original_count, side,
            original, recovery, kShardBytes,
            "T32/Q3 final-source basis parity differs from oracle");
        CheckLegacyParity(original_count, side, original, recovery,
            kShardBytes,
            "T32/Q3 final-source basis parity differs from legacy leo_encode");
        FillVariableInput(input.bytes(), original_count, kShardBytes,
            UINT64_C(0x42414c414e434536));
    }

#if LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING && \
    !defined(LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR)
    if (original_count == 62 && side == 9)
    {
        Require(leopard2_internal::
                SetK62R8B64FusedEnabledForDiagnostics(false),
            "disable only the K62/R8 fused arithmetic leaf");
        std::memset(output.bytes(), 0xa5,
            static_cast<size_t>(side) * kShardBytes);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode(codec, kShardBytes,
            &original[0], &recovery[0], scratch.data(), scratch.size()),
            LEO2_SUCCESS, "execute R9 with the K62/R8 selector disabled");
        RequireTerminalCalls(1,
            "K62/R8-only selector disabled the adjacent R9 terminal");
        Require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                    t8_k62_b64_fused_calls == 0,
            "adjacent R9 encode entered the K62/R8 fused leaf");
        CheckVariableParity(field, generator, original_count, side,
            original, recovery, kShardBytes,
            "K62/R8-only selector changed adjacent R9 parity");
        Require(leopard2_internal::
                SetK62R8B64FusedEnabledForDiagnostics(true),
            "restore the K62/R8 fused arithmetic leaf");
    }
#endif

    if ((original_count == 62 ||
         (original_count == 79 && side == 32)) &&
        (side == 8 || side == 16 || side == 32 || side == 33))
    {
        RequireResult(leo2_test_codec_set_encode_mode(
            codec, LEO2_TEST_ENCODE_FORCE_TRANSFORM), LEO2_SUCCESS,
            "force mature selected-shape transform");
        std::memset(output.bytes(), 0xa5,
            static_cast<size_t>(side) * kShardBytes);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode(codec, kShardBytes,
            &original[0], &recovery[0], scratch.data(), scratch.size()),
            LEO2_SUCCESS, "execute forced mature selected-shape transform");
        RequireTerminalCalls(0,
            "forced transform entered generalized packed terminal");
        CheckVariableParity(field, generator, original_count, side,
            original, recovery, kShardBytes,
            "forced mature selected-shape parity differs from oracle");
        RequireResult(leo2_test_codec_set_encode_mode(
            codec, LEO2_TEST_ENCODE_AUTO), LEO2_SUCCESS,
            "restore selected-shape AUTO encode mode");
    }

    if (original_count == 65 && side == 9)
    {
        std::vector<uint8_t> detached(kShardBytes);
        std::memcpy(&detached[0], original[64], kShardBytes);
        original[64] = &detached[0];
        std::memset(output.bytes(), 0xa5,
            static_cast<size_t>(side) * kShardBytes);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode(codec, kShardBytes,
            &original[0], &recovery[0], scratch.data(), scratch.size()),
            LEO2_SUCCESS, "execute detached K65/R9 fallback");
        RequireTerminalCalls(0,
            "detached K65/R9 source entered the packed terminal");
        CheckVariableParity(field, generator, original_count, side,
            original, recovery, kShardBytes,
            "detached K65/R9 parity differs from oracle");
        original[64] = input.bytes() + 64U * kShardBytes;

        void* const saved_recovery = recovery[8];
        recovery[8] = NULL;
        std::memset(output.bytes(), 0xa5,
            static_cast<size_t>(side) * kShardBytes);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode(codec, kShardBytes,
            &original[0], &recovery[0], scratch.data(), scratch.size()),
            LEO2_SUCCESS, "execute sparse-output K65/R9 fallback");
        RequireTerminalCalls(0,
            "sparse-output K65/R9 entered the packed terminal");
        CheckVariableParity(field, generator, original_count, side - 1U,
            original, recovery, kShardBytes,
            "sparse-output K65/R9 parity differs from oracle");
        recovery[8] = saved_recovery;
    }

    if ((original_count == 33 || original_count == 79) && side == 32)
    {
        /* A detached final source must retain the mature path: the terminal
           is a packed-layout optimization, not a new API restriction. */
        std::vector<uint8_t> detached(kShardBytes);
        const unsigned detached_source = original_count - 1U;
        std::memcpy(&detached[0], original[detached_source], kShardBytes);
        original[detached_source] = &detached[0];
        std::memset(output.bytes(), 0xa5,
            static_cast<size_t>(side) * kShardBytes);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode(codec, kShardBytes,
            &original[0], &recovery[0], scratch.data(), scratch.size()),
            LEO2_SUCCESS, "execute detached K33/R32 fallback");
        RequireTerminalCalls(0,
            "detached K33/R32 source entered the packed terminal");
        CheckVariableParity(field, generator, original_count, side,
            original, recovery, kShardBytes,
            "detached K33/R32 parity differs from oracle");

        original[detached_source] = input.bytes() +
            static_cast<size_t>(detached_source) * kShardBytes;
        if (original_count == 79)
        {
            void* const saved_recovery = recovery[31];
            recovery[31] = NULL;
            std::memset(output.bytes(), 0xa5,
                static_cast<size_t>(side) * kShardBytes);
            leopard::ff8::TestOnlyResetHighEncodeCounts();
            RequireResult(leo2_encode(codec, kShardBytes,
                &original[0], &recovery[0], scratch.data(), scratch.size()),
                LEO2_SUCCESS, "execute sparse-output K79/R32 fallback");
            RequireTerminalCalls(0,
                "sparse-output K79/R32 entered the packed terminal");
            CheckVariableParity(field, generator, original_count, side - 1U,
                original, recovery, kShardBytes,
                "sparse-output K79/R32 parity differs from oracle");
            recovery[31] = saved_recovery;
        }
        Require(leopard2_internal::
                SetBalancedB64TerminalEnabledForDiagnostics(false),
            "disable the K33/R32 terminal");
        std::memset(output.bytes(), 0xa5,
            static_cast<size_t>(side) * kShardBytes);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode(codec, kShardBytes,
            &original[0], &recovery[0], scratch.data(), scratch.size()),
            LEO2_SUCCESS, "execute same-binary K33/R32 control route");
        RequireTerminalCalls(0,
            "same-binary K33/R32 control entered the packed terminal");
        CheckVariableParity(field, generator, original_count, side,
            original, recovery, kShardBytes,
            "same-binary K33/R32 control differs from oracle");
        Require(leopard2_internal::
                SetBalancedB64TerminalEnabledForDiagnostics(true),
            "restore the K33/R32 terminal");

        const size_t byte_neighbors[] = { 63U, 65U };
        for (size_t neighbor_index = 0;
             neighbor_index < sizeof(byte_neighbors) /
                 sizeof(byte_neighbors[0]); ++neighbor_index)
        {
            const size_t neighbor_bytes = byte_neighbors[neighbor_index];
            size_t neighbor_scratch_bytes = 0;
            RequireResult(leo2_encode_scratch_size(codec, neighbor_bytes,
                &neighbor_scratch_bytes), LEO2_SUCCESS,
                "query K33/R32 neighboring-byte scratch");
            AlignedBuffer neighbor_scratch(neighbor_scratch_bytes);
            AlignedBuffer neighbor_input(
                static_cast<size_t>(original_count) * neighbor_bytes);
            AlignedBuffer neighbor_output(
                static_cast<size_t>(side) * neighbor_bytes);
            FillVariableInput(neighbor_input.bytes(), original_count,
                neighbor_bytes, UINT64_C(0x4b33334e45494748));
            std::vector<const void*> neighbor_original;
            std::vector<void*> neighbor_recovery;
            SetVariablePackedPointers(neighbor_input.bytes(),
                neighbor_output.bytes(), original_count, side,
                neighbor_bytes, neighbor_original, neighbor_recovery);
            leopard::ff8::TestOnlyResetHighEncodeCounts();
            RequireResult(leo2_encode(codec, neighbor_bytes,
                &neighbor_original[0], &neighbor_recovery[0],
                neighbor_scratch.data(), neighbor_scratch.size()),
                LEO2_SUCCESS,
                "execute K33/R32 neighboring-byte mature path");
            RequireTerminalCalls(0,
                "neighboring byte count entered the K33/R32 terminal");
            CheckVariableParity(field, generator, original_count, side,
                neighbor_original, neighbor_recovery, neighbor_bytes,
                "K33/R32 neighboring-byte parity differs from oracle");
        }
    }

    if ((original_count == 62 || original_count == 65 ||
         original_count == 66) &&
        side >= 8 && side <= 64)
    {
        const bool k66_target = original_count == 66 && side == 16;
        Require(k62_fused_shape
                ? leopard2_internal::
                    SetK62R8B64FusedEnabledForDiagnostics(false)
                : k66_target
                    ? leopard2_internal::
                        SetK66R16B64TailEnabledForDiagnostics(false)
                    : leopard2_internal::
                        SetBalancedB64TerminalEnabledForDiagnostics(false),
            "disable the selected K62/K65/K66 arithmetic route");
        std::memset(output.bytes(), 0xa5,
            static_cast<size_t>(side) * kShardBytes);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode(codec, kShardBytes,
            &original[0], &recovery[0], scratch.data(), scratch.size()),
            LEO2_SUCCESS, "execute same-binary K62 control route");
        RequireTerminalCalls(k66_target ? 1U : 0U,
            "same-binary control used the wrong packed terminal route");
        if (k66_target)
        {
            Require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                        tail_column_calls == 0,
                "same-binary K66 control entered the fused tail leaf");
        }
        CheckVariableParity(field, generator, original_count, side,
            original, recovery, kShardBytes,
            "same-binary K62 control differs from oracle");
        Require(k62_fused_shape
                ? leopard2_internal::
                    SetK62R8B64FusedEnabledForDiagnostics(true)
                : k66_target
                    ? leopard2_internal::
                        SetK66R16B64TailEnabledForDiagnostics(true)
                    : leopard2_internal::
                        SetBalancedB64TerminalEnabledForDiagnostics(true),
            "restore the selected K62/K65/K66 arithmetic route");

        std::memset(output.bytes(), 0xa5,
            static_cast<size_t>(side) * kShardBytes);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode(codec, kShardBytes,
            &original[0], &recovery[0], scratch.data(), scratch.size()),
            LEO2_SUCCESS, "reexecute restored K62 terminal");
        RequireTerminalCalls(k62_fused_shape ? 0U : 1U,
            "restored K62 route used the wrong public terminal");
        Require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                    t8_k62_b64_fused_calls ==
                (k62_fused_shape ? 1U : 0U),
            "restored K62 route used the wrong fused arithmetic leaf");
        CheckVariableParity(field, generator, original_count, side,
            original, recovery, kShardBytes,
            "restored K62 terminal differs from oracle");
    }

    if (original_count == 16 && side == 16)
    {
        /* The codec is GF(2)-linear.  Exercising every systematic basis row
           at both 32-byte slices proves the generated transform matrix rather
           than relying only on a favorable random vector. */
        for (unsigned source = 0; source < side; ++source)
        {
            std::memset(input.bytes(), 0,
                static_cast<size_t>(original_count) * kShardBytes);
            std::memset(output.bytes(), 0xa5,
                static_cast<size_t>(side) * kShardBytes);
            input.bytes()[static_cast<size_t>(source) * kShardBytes] = 1;
            input.bytes()[static_cast<size_t>(source) * kShardBytes + 63] =
                0xa7;
            leopard::ff8::TestOnlyResetHighEncodeCounts();
            RequireResult(leo2_encode(codec, kShardBytes,
                &original[0], &recovery[0], scratch.data(), scratch.size()),
                LEO2_SUCCESS, "execute T16 systematic basis message");
            RequireTerminalCalls(1,
                "T16 basis message missed the generated terminal");
            CheckVariableParity(field, generator, side, side,
                original, recovery, kShardBytes,
                "T16 basis parity differs from direct oracle");
        }
    }

    if (original_count == 62 && side == 8)
    {
        /* Exercise every systematic generator column at both AVX2 slices.
           This includes all six live rows in the shortened final block. */
        for (unsigned source = 0; source < original_count; ++source)
        {
            std::memset(input.bytes(), 0,
                static_cast<size_t>(original_count) * kShardBytes);
            std::memset(output.bytes(), 0xa5,
                static_cast<size_t>(side) * kShardBytes);
            input.bytes()[static_cast<size_t>(source) * kShardBytes] = 1;
            input.bytes()[static_cast<size_t>(source) * kShardBytes + 63] =
                0xa7;
            leopard::ff8::TestOnlyResetHighEncodeCounts();
            RequireResult(leo2_encode(codec, kShardBytes,
                &original[0], &recovery[0], scratch.data(), scratch.size()),
                LEO2_SUCCESS, "execute K62/R8 systematic basis message");
            RequireTerminalCalls(0,
                "K62/R8 basis message entered a public packed terminal");
            const leopard::ff8::TestOnlyHighEncodeCounts basis_counts =
                leopard::ff8::TestOnlyGetHighEncodeCounts();
            Require(basis_counts.t8_k62_b64_fused_calls == 1,
                "K62/R8 basis message missed the fused leaf");
            CheckVariableParity(field, generator, original_count, side,
                original, recovery, kShardBytes,
                "K62/R8 basis parity differs from direct oracle");
        }
    }

    if (original_count == 62 && side == 8)
    {
        FillVariableInput(input.bytes(), original_count, kShardBytes,
            UINT64_C(0x4b36324445544143));
        std::vector<uint8_t> detached(kShardBytes);
        std::memcpy(&detached[0], original[61], kShardBytes);
        original[61] = &detached[0];
        std::memset(output.bytes(), 0xa5,
            static_cast<size_t>(side) * kShardBytes);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode(codec, kShardBytes,
            &original[0], &recovery[0], scratch.data(), scratch.size()),
            LEO2_SUCCESS, "execute detached K62/R8 fallback");
        RequireTerminalCalls(0,
            "detached K62/R8 source entered the packed terminal");
        CheckVariableParity(field, generator, original_count, side,
            original, recovery, kShardBytes,
            "detached K62/R8 parity differs from direct oracle");
        original[61] = input.bytes() + 61U * kShardBytes;

        void* const saved_recovery = recovery[7];
        recovery[7] = NULL;
        std::memset(output.bytes(), 0xa5,
            static_cast<size_t>(side) * kShardBytes);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode(codec, kShardBytes,
            &original[0], &recovery[0], scratch.data(), scratch.size()),
            LEO2_SUCCESS, "execute sparse-output K62/R8 fallback");
        RequireTerminalCalls(0,
            "sparse-output K62/R8 entered the packed terminal");
        CheckVariableParity(field, generator, original_count, side - 1U,
            original, recovery, kShardBytes,
            "sparse-output K62/R8 parity differs from direct oracle");
        recovery[7] = saved_recovery;

        const size_t byte_neighbors[] = { 63U, 65U };
        for (size_t neighbor_index = 0;
             neighbor_index < sizeof(byte_neighbors) /
                 sizeof(byte_neighbors[0]); ++neighbor_index)
        {
            const size_t neighbor_bytes = byte_neighbors[neighbor_index];
            size_t neighbor_scratch_bytes = 0;
            RequireResult(leo2_encode_scratch_size(codec, neighbor_bytes,
                &neighbor_scratch_bytes), LEO2_SUCCESS,
                "query K62/R8 neighboring-byte scratch");
            AlignedBuffer neighbor_scratch(neighbor_scratch_bytes);
            AlignedBuffer neighbor_input(
                static_cast<size_t>(original_count) * neighbor_bytes);
            AlignedBuffer neighbor_output(
                static_cast<size_t>(side) * neighbor_bytes);
            FillVariableInput(neighbor_input.bytes(), original_count,
                neighbor_bytes, UINT64_C(0x4b363252384e4549));
            std::vector<const void*> neighbor_original;
            std::vector<void*> neighbor_recovery;
            SetVariablePackedPointers(neighbor_input.bytes(),
                neighbor_output.bytes(), original_count, side,
                neighbor_bytes, neighbor_original, neighbor_recovery);
            leopard::ff8::TestOnlyResetHighEncodeCounts();
            RequireResult(leo2_encode(codec, neighbor_bytes,
                &neighbor_original[0], &neighbor_recovery[0],
                neighbor_scratch.data(), neighbor_scratch.size()),
                LEO2_SUCCESS,
                "execute K62/R8 neighboring-byte mature path");
            RequireTerminalCalls(0,
                "neighboring byte count entered the K62/R8 terminal");
            CheckVariableParity(field, generator, original_count, side,
                neighbor_original, neighbor_recovery, neighbor_bytes,
                "K62/R8 neighboring-byte parity differs from oracle");
        }
    }
    leo2_codec_destroy(codec);
}

void ExerciseK62MatureScope(leo2_context* context)
{
    static const unsigned original_count = 62;
    static const unsigned recovery_count = 8;
    static const size_t batch_count = 2;

    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, original_count, recovery_count,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        LEO2_SUCCESS, "create K62/R8 mature-scope codec");
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(codec, kShardBytes,
        &scratch_bytes), LEO2_SUCCESS,
        "query K62/R8 mature-scope scratch");
    AlignedBuffer input0(
        static_cast<size_t>(original_count) * kShardBytes);
    AlignedBuffer input1(
        static_cast<size_t>(original_count) * kShardBytes);
    AlignedBuffer output0(
        static_cast<size_t>(recovery_count) * kShardBytes);
    AlignedBuffer output1(
        static_cast<size_t>(recovery_count) * kShardBytes);
    AlignedBuffer scratch0(scratch_bytes);
    AlignedBuffer scratch1(scratch_bytes);
    std::vector<std::vector<const void*> > original(
        batch_count, std::vector<const void*>(original_count));
    std::vector<std::vector<void*> > recovery(
        batch_count, std::vector<void*>(recovery_count));
    SetVariablePackedPointers(input0.bytes(), output0.bytes(),
        original_count, recovery_count, kShardBytes,
        original[0], recovery[0]);
    SetVariablePackedPointers(input1.bytes(), output1.bytes(),
        original_count, recovery_count, kShardBytes,
        original[1], recovery[1]);
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh,
            original_count, recovery_count);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);
    leo2_encode_batch_item items[batch_count] = {
        { kShardBytes, &original[0][0], &recovery[0][0],
            scratch0.data(), scratch0.size() },
        { kShardBytes, &original[1][0], &recovery[1][0],
            scratch1.data(), scratch1.size() }
    };

    FillVariableInput(input0.bytes(), original_count, kShardBytes,
        UINT64_C(0x4b36324d554c5449));
    FillVariableInput(input1.bytes(), original_count, kShardBytes,
        UINT64_C(0x4b36324d554c5432));
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch(codec, items, batch_count),
        LEO2_SUCCESS, "execute K62/R8 multi-item mature batch");
    leopard::ff8::TestOnlyHighEncodeCounts counts =
        leopard::ff8::TestOnlyGetHighEncodeCounts();
    Require(counts.t8_k62_b64_fused_calls == 0,
        "K62/R8 multi-item batch entered the one-item fused leaf");
    Require(counts.balanced_b64_packed_calls == 0,
        "K62/R8 multi-item batch entered the public packed terminal");
    for (size_t item_i = 0; item_i < batch_count; ++item_i)
    {
        CheckVariableParity(field, generator, original_count,
            recovery_count, original[item_i], recovery[item_i],
            kShardBytes,
            "K62/R8 multi-item mature parity differs from oracle");
    }

    leo2_encode_batch_binding* binding = NULL;
    RequireResult(leo2_encode_batch_binding_create(
        codec, items, batch_count, &binding), LEO2_SUCCESS,
        "create reusable K62/R8 mature binding");
    Require(leo2_encode_batch_binding_item_count(binding) == batch_count,
        "K62/R8 reusable binding lost an item");
    for (unsigned execution = 0; execution < 2; ++execution)
    {
        FillVariableInput(input0.bytes(), original_count, kShardBytes,
            UINT64_C(0x4b363242494e4430) + execution);
        FillVariableInput(input1.bytes(), original_count, kShardBytes,
            UINT64_C(0x4b363242494e4431) + execution);
        std::memset(output0.bytes(), 0xa5, output0.size());
        std::memset(output1.bytes(), 0xa5, output1.size());
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode_batch_binding_execute(binding),
            LEO2_SUCCESS, "execute reusable K62/R8 mature binding");
        counts = leopard::ff8::TestOnlyGetHighEncodeCounts();
        Require(counts.t8_k62_b64_fused_calls == 0,
            "K62/R8 reusable binding entered the one-item fused leaf");
        Require(counts.balanced_b64_packed_calls == 0,
            "K62/R8 reusable binding entered the public packed terminal");
        for (size_t item_i = 0; item_i < batch_count; ++item_i)
        {
            CheckVariableParity(field, generator, original_count,
                recovery_count, original[item_i], recovery[item_i],
                kShardBytes,
                "K62/R8 reusable binding parity differs from oracle");
        }
    }
    leo2_encode_batch_binding_destroy(binding);

    /* A non-null but detached recovery row is legal and must retain the
       mature per-pointer validation path. */
    std::vector<uint8_t> detached_recovery(kShardBytes, 0xa5);
    void* const saved_recovery = recovery[0][7];
    recovery[0][7] = &detached_recovery[0];
    std::memset(output0.bytes(), 0xa5, output0.size());
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes,
        &original[0][0], &recovery[0][0], scratch0.data(), scratch0.size()),
        LEO2_SUCCESS, "execute detached-recovery K62/R8 fallback");
    counts = leopard::ff8::TestOnlyGetHighEncodeCounts();
    Require(counts.t8_k62_b64_fused_calls == 0,
        "detached K62/R8 recovery entered the fused leaf");
    Require(counts.balanced_b64_packed_calls == 0,
        "detached K62/R8 recovery entered the packed terminal");
    CheckVariableParity(field, generator, original_count, recovery_count,
        original[0], recovery[0], kShardBytes,
        "detached K62/R8 recovery parity differs from oracle");
    recovery[0][7] = saved_recovery;
    leo2_codec_destroy(codec);

    leo2_codec* low_codec = NULL;
    RequireResult(leo2_codec_create(context, original_count, recovery_count,
        LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8, NULL, &low_codec),
        LEO2_SUCCESS, "create low-profile K62/R8 exclusion codec");
    size_t low_scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(low_codec, kShardBytes,
        &low_scratch_bytes), LEO2_SUCCESS,
        "query low-profile K62/R8 exclusion scratch");
    AlignedBuffer low_scratch(low_scratch_bytes);
    const leopard2_test::ProfileLayout low_layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLow, original_count, recovery_count);
    const leopard2_test::Matrix low_generator =
        leopard2_test::direct_systematic_generator(field, low_layout);
    FillVariableInput(input0.bytes(), original_count, kShardBytes,
        UINT64_C(0x4b36324c4f575052));
    std::memset(output0.bytes(), 0xa5, output0.size());
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(low_codec, kShardBytes,
        &original[0][0], &recovery[0][0],
        low_scratch.data(), low_scratch.size()), LEO2_SUCCESS,
        "execute low-profile K62/R8 exclusion");
    counts = leopard::ff8::TestOnlyGetHighEncodeCounts();
    Require(counts.t8_k62_b64_fused_calls == 0,
        "low-profile K62/R8 entered the legacy-high fused leaf");
    Require(counts.balanced_b64_packed_calls == 0,
        "low-profile K62/R8 entered the legacy-high packed terminal");
    CheckVariableParity(field, low_generator, original_count,
        recovery_count, original[0], recovery[0], kShardBytes,
        "low-profile K62/R8 parity differs from oracle");
    leo2_codec_destroy(low_codec);

    leo2_codec* gf16_codec = NULL;
    if (!CreateOptionalGF16Codec(context, original_count, recovery_count,
            &gf16_codec, "create GF16 K62/R8 exclusion codec"))
        return;
    size_t gf16_scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(gf16_codec, kShardBytes,
        &gf16_scratch_bytes), LEO2_SUCCESS,
        "query GF16 K62/R8 exclusion scratch");
    AlignedBuffer gf16_scratch(gf16_scratch_bytes);
    const leopard2_test::BinaryField gf16 =
        leopard2_test::make_legacy_gf16();
    const leopard2_test::Matrix gf16_generator =
        leopard2_test::direct_systematic_generator(gf16, layout);
    FillVariableInput(input0.bytes(), original_count, kShardBytes,
        UINT64_C(0x4b36324746313650));
    std::memset(output0.bytes(), 0xa5, output0.size());
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(gf16_codec, kShardBytes,
        &original[0][0], &recovery[0][0],
        gf16_scratch.data(), gf16_scratch.size()), LEO2_SUCCESS,
        "execute GF16 K62/R8 exclusion");
    counts = leopard::ff8::TestOnlyGetHighEncodeCounts();
    Require(counts.t8_k62_b64_fused_calls == 0,
        "GF16 K62/R8 entered the GF8 fused leaf");
    Require(counts.balanced_b64_packed_calls == 0,
        "GF16 K62/R8 entered the GF8 packed terminal");
    CheckVariableParityGF16(gf16, gf16_generator, original_count,
        recovery_count, original[0], recovery[0], kShardBytes,
        "GF16 K62/R8 parity differs from oracle");
    leo2_codec_destroy(gf16_codec);
}

void ExerciseK79ProfileFieldExclusions(leo2_context* context)
{
    static const unsigned original_count = 79;
    static const unsigned recovery_count = 32;
    AlignedBuffer input(
        static_cast<size_t>(original_count) * kShardBytes);
    AlignedBuffer output(
        static_cast<size_t>(recovery_count) * kShardBytes);
    FillVariableInput(input.bytes(), original_count, kShardBytes,
        UINT64_C(0x4b373950524f4649));
    std::vector<const void*> original;
    std::vector<void*> recovery;
    SetVariablePackedPointers(input.bytes(), output.bytes(), original_count,
        recovery_count, kShardBytes, original, recovery);

    const leopard2_test::BinaryField gf8 =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout high_gf8_layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh, original_count, recovery_count);
    const leopard2_test::Matrix high_gf8_generator =
        leopard2_test::direct_systematic_generator(gf8, high_gf8_layout);
    leo2_codec* auto_codec = NULL;
    RequireResult(leo2_codec_create(context, original_count, recovery_count,
        LEO2_PROFILE_AUTO, LEO2_FIELD_AUTO, NULL, &auto_codec),
        LEO2_SUCCESS, "create AUTO K79/R32 codec");
    Require(leo2_codec_profile(auto_codec) ==
            LEO2_PROFILE_LEGACY_HIGH_V1 &&
            leo2_codec_field(auto_codec) == LEO2_FIELD_GF8,
        "AUTO K79/R32 resolved to the wrong wire profile or field");
    size_t auto_scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        auto_codec, kShardBytes, &auto_scratch_bytes), LEO2_SUCCESS,
        "query AUTO K79/R32 scratch");
    AlignedBuffer auto_scratch(auto_scratch_bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(auto_codec, kShardBytes,
        &original[0], &recovery[0], auto_scratch.data(), auto_scratch.size()),
        LEO2_SUCCESS, "execute AUTO K79/R32 encode");
    RequireTerminalCalls(1,
        "AUTO K79/R32 encode missed the packed terminal");
    CheckVariableParity(gf8, high_gf8_generator, original_count,
        recovery_count, original, recovery, kShardBytes,
        "AUTO K79/R32 parity differs from oracle");
    leo2_codec_destroy(auto_codec);

    const leopard2_test::ProfileLayout low_layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLow, original_count, recovery_count);
    const leopard2_test::Matrix low_generator =
        leopard2_test::direct_systematic_generator(gf8, low_layout);
    leo2_codec* low_codec = NULL;
    RequireResult(leo2_codec_create(context, original_count, recovery_count,
        LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8, NULL, &low_codec),
        LEO2_SUCCESS, "create low-profile K79/R32 exclusion codec");
    size_t low_scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        low_codec, kShardBytes, &low_scratch_bytes), LEO2_SUCCESS,
        "query low-profile K79/R32 exclusion scratch");
    AlignedBuffer low_scratch(low_scratch_bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(low_codec, kShardBytes,
        &original[0], &recovery[0], low_scratch.data(), low_scratch.size()),
        LEO2_SUCCESS, "execute low-profile K79/R32 exclusion");
    RequireTerminalCalls(0,
        "low-profile K79/R32 entered the legacy-high packed terminal");
    CheckVariableParity(gf8, low_generator, original_count, recovery_count,
        original, recovery, kShardBytes,
        "low-profile K79/R32 parity differs from oracle");
    leo2_codec_destroy(low_codec);

    leo2_codec* gf16_codec = NULL;
    if (!CreateOptionalGF16Codec(context, original_count, recovery_count,
            &gf16_codec, "create GF16 K79/R32 exclusion codec"))
        return;
    const leopard2_test::BinaryField gf16 =
        leopard2_test::make_legacy_gf16();
    const leopard2_test::ProfileLayout high_layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh, original_count, recovery_count);
    const leopard2_test::Matrix high_generator =
        leopard2_test::direct_systematic_generator(gf16, high_layout);
    size_t gf16_scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        gf16_codec, kShardBytes, &gf16_scratch_bytes), LEO2_SUCCESS,
        "query GF16 K79/R32 exclusion scratch");
    AlignedBuffer gf16_scratch(gf16_scratch_bytes);
    std::memset(output.bytes(), 0xa5, output.size());
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(gf16_codec, kShardBytes,
        &original[0], &recovery[0], gf16_scratch.data(),
        gf16_scratch.size()), LEO2_SUCCESS,
        "execute GF16 K79/R32 exclusion");
    RequireTerminalCalls(0,
        "GF16 K79/R32 entered the GF8 packed terminal");
    CheckVariableParityGF16(gf16, high_generator, original_count,
        recovery_count, original, recovery, kShardBytes,
        "GF16 K79/R32 parity differs from oracle");
    leo2_codec_destroy(gf16_codec);
}

#if LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED
void ExerciseK65R9TailBasis(leo2_context* context)
{
    static const unsigned original_count = 65;
    static const unsigned recovery_count = 9;

    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, original_count, recovery_count,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        LEO2_SUCCESS, "create K65/R9 tail-basis codec");
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(codec, kShardBytes,
        &scratch_bytes), LEO2_SUCCESS,
        "query K65/R9 tail-basis scratch");
    AlignedBuffer scratch(scratch_bytes);
    AlignedBuffer input(
        static_cast<size_t>(original_count) * kShardBytes);
    AlignedBuffer output(
        static_cast<size_t>(recovery_count) * kShardBytes);
    // Isolate source 64 so the ninth output detects a lost one-row remainder
    // in the K65 tail multiply-add path.
    for (size_t offset = 0; offset < kShardBytes; ++offset)
    {
        input.bytes()[64U * kShardBytes + offset] =
            static_cast<uint8_t>(offset + 1U);
    }

    std::vector<const void*> original;
    std::vector<void*> recovery;
    SetVariablePackedPointers(input.bytes(), output.bytes(), original_count,
        recovery_count, kShardBytes, original, recovery);
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(leopard2_test::kLegacyHigh,
            original_count, recovery_count);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);

    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute K65/R9 tail-basis terminal");
    RequireTerminalCalls(1,
        "K65/R9 tail-basis message missed the packed terminal");
    CheckVariableParity(field, generator, original_count, recovery_count,
        original, recovery, kShardBytes,
        "K65/R9 tail-basis parity differs from direct oracle");

    const unsigned work_count =
        leo_encode_work_count(original_count, recovery_count);
    Require(work_count >= recovery_count,
        "legacy work-count query rejected K65/R9");
    AlignedBuffer legacy_storage(
        static_cast<size_t>(work_count) * kShardBytes);
    std::vector<void*> legacy_work(work_count);
    for (unsigned i = 0; i < work_count; ++i)
        legacy_work[i] = legacy_storage.bytes() +
            static_cast<size_t>(i) * kShardBytes;
    Require(leo_encode(kShardBytes, original_count, recovery_count,
            work_count, &original[0], &legacy_work[0]) == Leopard_Success,
        "legacy K65/R9 tail-basis encode failed");
    for (unsigned parity = 0; parity < recovery_count; ++parity)
    {
        Require(std::memcmp(recovery[parity], legacy_work[parity],
                kShardBytes) == 0,
            "K65/R9 tail-basis parity differs from legacy leo_encode");
    }
    Require(static_cast<const uint8_t*>(recovery[8])[0] != 0,
        "K65/R9 tail contribution to parity 8 unexpectedly vanished");
    Require(static_cast<const uint8_t*>(recovery[8])[0] ==
            static_cast<const uint8_t*>(legacy_work[8])[0],
        "K65/R9 parity 8 tail contribution differs from legacy leo_encode");
    leo2_codec_destroy(codec);
}

void ExerciseK66R16TailBasis(leo2_context* context)
{
    static const unsigned original_count = 66;
    static const unsigned recovery_count = 16;
    static const uint8_t expected_tail_elements[recovery_count] = {
        127, 85, 96, 251, 80, 120, 244, 102,
        173, 197, 134, 68, 201, 167, 64, 142
    };

    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, original_count, recovery_count,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        LEO2_SUCCESS, "create K66/R16 tail-basis codec");
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(codec, kShardBytes,
        &scratch_bytes), LEO2_SUCCESS,
        "query K66/R16 tail-basis scratch");
    AlignedBuffer scratch(scratch_bytes);
    AlignedBuffer input(
        static_cast<size_t>(original_count) * kShardBytes);
    AlignedBuffer output(
        static_cast<size_t>(recovery_count) * kShardBytes);
    std::memset(input.bytes() + 65U * kShardBytes, 1, kShardBytes);

    std::vector<const void*> original;
    std::vector<void*> recovery;
    SetVariablePackedPointers(input.bytes(), output.bytes(), original_count,
        recovery_count, kShardBytes, original, recovery);
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(leopard2_test::kLegacyHigh,
            original_count, recovery_count);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);

    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute K66/R16 tail-basis terminal");
    RequireTerminalCalls(1,
        "K66/R16 tail-basis message missed the packed terminal");
    Require(leopard::ff8::TestOnlyGetHighEncodeCounts().tail_column_calls == 1,
        "K66/R16 tail-basis message missed the second-tail leaf");
    CheckVariableParity(field, generator, original_count, recovery_count,
        original, recovery, kShardBytes,
        "K66/R16 tail-basis parity differs from direct oracle");
    for (unsigned parity = 0; parity < recovery_count; ++parity)
    {
        const uint8_t* const bytes =
            static_cast<const uint8_t*>(recovery[parity]);
        for (size_t offset = 0; offset < kShardBytes; ++offset)
        {
            Require(bytes[offset] == expected_tail_elements[parity],
                "K66/R16 coordinate-81 golden column differs");
        }
    }

    const unsigned work_count =
        leo_encode_work_count(original_count, recovery_count);
    Require(work_count >= recovery_count,
        "legacy work-count query rejected K66/R16");
    AlignedBuffer legacy_storage(
        static_cast<size_t>(work_count) * kShardBytes);
    std::vector<void*> legacy_work(work_count);
    for (unsigned i = 0; i < work_count; ++i)
    {
        legacy_work[i] = legacy_storage.bytes() +
            static_cast<size_t>(i) * kShardBytes;
    }
    Require(leo_encode(kShardBytes, original_count, recovery_count,
            work_count, &original[0], &legacy_work[0]) == Leopard_Success,
        "legacy K66/R16 tail-basis encode failed");
    for (unsigned parity = 0; parity < recovery_count; ++parity)
    {
        Require(std::memcmp(recovery[parity], legacy_work[parity],
                kShardBytes) == 0,
            "K66/R16 tail-basis parity differs from legacy leo_encode");
    }

    RequireResult(leo2_test_codec_set_encode_mode(
        codec, LEO2_TEST_ENCODE_FORCE_TRANSFORM), LEO2_SUCCESS,
        "force mature K66/R16 transform");
    std::memset(output.bytes(), 0xa5, output.size());
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute forced mature K66/R16 transform");
    RequireTerminalCalls(0,
        "forced K66/R16 transform entered the packed terminal");
    Require(leopard::ff8::TestOnlyGetHighEncodeCounts().tail_column_calls == 0,
        "forced K66/R16 transform entered the fused tail leaf");
    CheckVariableParity(field, generator, original_count, recovery_count,
        original, recovery, kShardBytes,
        "forced mature K66/R16 parity differs from direct oracle");
    RequireResult(leo2_test_codec_set_encode_mode(
        codec, LEO2_TEST_ENCODE_AUTO), LEO2_SUCCESS,
        "restore K66/R16 AUTO encode mode");

    leo2_codec* low_codec = NULL;
    RequireResult(leo2_codec_create(context, original_count, recovery_count,
        LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8, NULL, &low_codec),
        LEO2_SUCCESS, "create low-profile K66/R16 exclusion codec");
    size_t low_scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(low_codec, kShardBytes,
        &low_scratch_bytes), LEO2_SUCCESS,
        "query low-profile K66/R16 exclusion scratch");
    AlignedBuffer low_scratch(low_scratch_bytes);
    const leopard2_test::ProfileLayout low_layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLow, original_count, recovery_count);
    const leopard2_test::Matrix low_generator =
        leopard2_test::direct_systematic_generator(field, low_layout);
    std::memset(output.bytes(), 0xa5, output.size());
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(low_codec, kShardBytes,
        &original[0], &recovery[0], low_scratch.data(), low_scratch.size()),
        LEO2_SUCCESS, "execute low-profile K66/R16 exclusion");
    RequireTerminalCalls(0,
        "low-profile K66/R16 entered the legacy-high packed terminal");
    Require(leopard::ff8::TestOnlyGetHighEncodeCounts().tail_column_calls == 0,
        "low-profile K66/R16 entered the legacy-high tail leaf");
    CheckVariableParity(field, low_generator, original_count, recovery_count,
        original, recovery, kShardBytes,
        "low-profile K66/R16 parity differs from direct oracle");
    leo2_codec_destroy(low_codec);

    leo2_codec* gf16_codec = NULL;
    if (CreateOptionalGF16Codec(context, original_count, recovery_count,
            &gf16_codec, "create GF16 K66/R16 exclusion codec"))
    {
        size_t gf16_scratch_bytes = 0;
        RequireResult(leo2_encode_scratch_size(gf16_codec, kShardBytes,
            &gf16_scratch_bytes), LEO2_SUCCESS,
            "query GF16 K66/R16 exclusion scratch");
        AlignedBuffer gf16_scratch(gf16_scratch_bytes);
        const leopard2_test::BinaryField gf16 =
            leopard2_test::make_legacy_gf16();
        const leopard2_test::Matrix gf16_generator =
            leopard2_test::direct_systematic_generator(gf16, layout);
        std::memset(output.bytes(), 0xa5, output.size());
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode(gf16_codec, kShardBytes,
            &original[0], &recovery[0],
            gf16_scratch.data(), gf16_scratch.size()), LEO2_SUCCESS,
            "execute GF16 K66/R16 exclusion");
        RequireTerminalCalls(0,
            "GF16 K66/R16 entered the GF8 packed terminal");
        Require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                    tail_column_calls == 0,
            "GF16 K66/R16 entered the GF8 tail leaf");
        CheckVariableParityGF16(gf16, gf16_generator, original_count,
            recovery_count, original, recovery, kShardBytes,
            "GF16 K66/R16 parity differs from direct oracle");
        leo2_codec_destroy(gf16_codec);
    }

    std::vector<uint8_t> detached(kShardBytes);
    std::memcpy(&detached[0], original[65], kShardBytes);
    original[65] = &detached[0];
    std::memset(output.bytes(), 0xa5, output.size());
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute detached K66/R16 fallback");
    RequireTerminalCalls(0,
        "detached K66/R16 source entered the packed terminal");
    Require(leopard::ff8::TestOnlyGetHighEncodeCounts().tail_column_calls == 0,
        "detached K66/R16 source entered the fused tail leaf");
    CheckVariableParity(field, generator, original_count, recovery_count,
        original, recovery, kShardBytes,
        "detached K66/R16 parity differs from direct oracle");
    original[65] = input.bytes() + 65U * kShardBytes;

    void* const saved_recovery = recovery[15];
    recovery[15] = NULL;
    std::memset(output.bytes(), 0xa5, output.size());
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute sparse-output K66/R16 fallback");
    RequireTerminalCalls(0,
        "sparse-output K66/R16 entered the packed terminal");
    Require(leopard::ff8::TestOnlyGetHighEncodeCounts().tail_column_calls == 0,
        "sparse-output K66/R16 entered the fused tail leaf");
    CheckVariableParity(field, generator, original_count, recovery_count - 1U,
        original, recovery, kShardBytes,
        "sparse-output K66/R16 parity differs from direct oracle");
    recovery[15] = saved_recovery;

    const size_t byte_neighbors[] = { 63U, 65U };
    for (size_t neighbor_i = 0;
         neighbor_i < sizeof(byte_neighbors) / sizeof(byte_neighbors[0]);
         ++neighbor_i)
    {
        const size_t neighbor_bytes = byte_neighbors[neighbor_i];
        size_t neighbor_scratch_bytes = 0;
        RequireResult(leo2_encode_scratch_size(codec, neighbor_bytes,
            &neighbor_scratch_bytes), LEO2_SUCCESS,
            "query K66/R16 neighboring-byte scratch");
        AlignedBuffer neighbor_scratch(neighbor_scratch_bytes);
        AlignedBuffer neighbor_input(
            static_cast<size_t>(original_count) * neighbor_bytes);
        AlignedBuffer neighbor_output(
            static_cast<size_t>(recovery_count) * neighbor_bytes);
        FillVariableInput(neighbor_input.bytes(), original_count,
            neighbor_bytes, UINT64_C(0x4b36365231364e45) + neighbor_i);
        std::vector<const void*> neighbor_original;
        std::vector<void*> neighbor_recovery;
        SetVariablePackedPointers(neighbor_input.bytes(),
            neighbor_output.bytes(), original_count, recovery_count,
            neighbor_bytes, neighbor_original, neighbor_recovery);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode(codec, neighbor_bytes,
            &neighbor_original[0], &neighbor_recovery[0],
            neighbor_scratch.data(), neighbor_scratch.size()),
            LEO2_SUCCESS, "execute K66/R16 neighboring-byte mature path");
        RequireTerminalCalls(0,
            "neighboring byte count entered the K66/R16 terminal");
        Require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                    tail_column_calls == 0,
            "neighboring byte count entered the K66/R16 tail leaf");
        CheckVariableParity(field, generator, original_count, recovery_count,
            neighbor_original, neighbor_recovery, neighbor_bytes,
            "K66/R16 neighboring-byte parity differs from oracle");
    }

    AlignedBuffer second_scratch(scratch_bytes);
    AlignedBuffer second_input(
        static_cast<size_t>(original_count) * kShardBytes);
    AlignedBuffer second_output(
        static_cast<size_t>(recovery_count) * kShardBytes);
    FillVariableInput(second_input.bytes(), original_count, kShardBytes,
        UINT64_C(0x4b36365231364232));
    std::vector<const void*> second_original;
    std::vector<void*> second_recovery;
    SetVariablePackedPointers(second_input.bytes(), second_output.bytes(),
        original_count, recovery_count, kShardBytes,
        second_original, second_recovery);
    leo2_encode_batch_item two_items[2] = {
        { kShardBytes, &original[0], &recovery[0],
            scratch.data(), scratch.size() },
        { kShardBytes, &second_original[0], &second_recovery[0],
            second_scratch.data(), second_scratch.size() }
    };
    std::memset(output.bytes(), 0xa5, output.size());
    std::memset(second_output.bytes(), 0xa5, second_output.size());
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch(codec, two_items, 2), LEO2_SUCCESS,
        "execute two-item K66/R16 mature batch");
    RequireTerminalCalls(0,
        "two-item K66/R16 batch entered the one-item terminal");
    Require(leopard::ff8::TestOnlyGetHighEncodeCounts().tail_column_calls == 0,
        "two-item K66/R16 batch entered the one-item tail leaf");
    CheckVariableParity(field, generator, original_count, recovery_count,
        original, recovery, kShardBytes,
        "first two-item K66/R16 parity differs from oracle");
    CheckVariableParity(field, generator, original_count, recovery_count,
        second_original, second_recovery, kShardBytes,
        "second two-item K66/R16 parity differs from oracle");

    leo2_encode_batch_item item = {
        kShardBytes, &original[0], &recovery[0], scratch.data(), scratch.size()
    };
    leo2_encode_batch_binding* binding = NULL;
    RequireResult(leo2_encode_batch_binding_create(
        codec, &item, 1, &binding), LEO2_SUCCESS,
        "create one-item reusable K66/R16 binding");
    std::memset(output.bytes(), 0xa5, output.size());
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch_binding_execute(binding), LEO2_SUCCESS,
        "execute one-item reusable K66/R16 binding");
    RequireTerminalCalls(0,
        "reusable K66/R16 binding entered the one-shot terminal");
    Require(leopard::ff8::TestOnlyGetHighEncodeCounts().tail_column_calls == 0,
        "reusable K66/R16 binding entered the one-shot tail leaf");
    CheckVariableParity(field, generator, original_count, recovery_count,
        original, recovery, kShardBytes,
        "reusable K66/R16 binding parity differs from oracle");
    leo2_encode_batch_binding_destroy(binding);

    leo2_codec_destroy(codec);
}
#endif

void ExerciseExcludedShape(leo2_context* context,
    unsigned original_count, unsigned recovery_count)
{
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, original_count, recovery_count,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        LEO2_SUCCESS, "create excluded balanced-neighbor codec");
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, kShardBytes, &scratch_bytes), LEO2_SUCCESS,
        "query excluded balanced-neighbor scratch");
    AlignedBuffer scratch(scratch_bytes);
    AlignedBuffer input(
        static_cast<size_t>(original_count) * kShardBytes);
    AlignedBuffer output(
        static_cast<size_t>(recovery_count) * kShardBytes);
    FillVariableInput(input.bytes(), original_count, kShardBytes,
        UINT64_C(0x4e45494748424f52) ^ recovery_count);
    std::vector<const void*> original;
    std::vector<void*> recovery;
    SetVariablePackedPointers(input.bytes(), output.bytes(), original_count,
        recovery_count, kShardBytes, original, recovery);

    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(leopard2_test::kLegacyHigh,
            original_count, recovery_count);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute excluded balanced-neighbor shape");
    const leopard::ff8::TestOnlyHighEncodeCounts counts =
        leopard::ff8::TestOnlyGetHighEncodeCounts();
    if (counts.balanced_b64_packed_calls != 0)
    {
        throw std::runtime_error(
            "excluded balanced-neighbor K=" +
            std::to_string(original_count) + "/R=" +
            std::to_string(recovery_count) + " entered packed terminal");
    }
    CheckVariableParity(field, generator, original_count, recovery_count,
        original, recovery, kShardBytes,
        "excluded balanced-neighbor parity differs from oracle");
    leo2_codec_destroy(codec);
}

void ExerciseNeighborBytes(
    leo2_codec* codec,
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator,
    size_t shard_bytes)
{
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, shard_bytes, &scratch_bytes), LEO2_SUCCESS,
        "query neighboring-byte scratch");
    AlignedBuffer scratch(scratch_bytes);
    AlignedBuffer input(kOriginalCount * shard_bytes);
    AlignedBuffer output(kRecoveryCount * shard_bytes);
    FillInput(input.bytes(), shard_bytes);
    const void* original[kOriginalCount];
    void* recovery[kRecoveryCount];
    SetPackedPointers(input.bytes(), output.bytes(), shard_bytes,
        original, recovery);

    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, original, recovery,
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute neighboring-byte mature path");
    RequireTerminalCalls(0,
        "neighboring byte count entered the exact-64 terminal");
    CheckParity(field, generator, original, recovery, shard_bytes);
}

void ExerciseAVX2Terminal(leo2_context* context)
{
    leo2_codec* codec = CreateCodec(
        context, LEO2_PROFILE_LEGACY_HIGH_V1);
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, kShardBytes, &scratch_bytes), LEO2_SUCCESS,
        "query K32/R32 terminal scratch");
    Require(scratch_bytes >= 5888,
        "K32/R32 scratch query lost the production upper bound");

    AlignedBuffer scratch(scratch_bytes + leo2_scratch_alignment());
    AlignedBuffer input(kOriginalCount * kShardBytes + 8);
    AlignedBuffer output(kRecoveryCount * kShardBytes + 8);
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh,
            kOriginalCount, kRecoveryCount);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);
    const void* original[kOriginalCount];
    void* recovery[kRecoveryCount];

    FillInput(input.bytes(), kShardBytes);
    SetPackedPointers(input.bytes(), output.bytes(), kShardBytes,
        original, recovery);
    std::memset(output.bytes(), 0xa5,
        kRecoveryCount * kShardBytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute packed public K32/R32 terminal");
    RequireTerminalCalls(1,
        "packed public encode did not select the T32 terminal");
    CheckParity(field, generator, original, recovery, kShardBytes);

    std::memset(output.bytes(), 0xa5,
        kRecoveryCount * kShardBytes);
    leo2_encode_batch_item item = {
        kShardBytes, original, recovery, scratch.data(), scratch_bytes
    };
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch(codec, &item, 1), LEO2_SUCCESS,
        "execute packed one-item K32/R32 batch terminal");
    RequireTerminalCalls(1,
        "one-item batch did not select the T32 terminal");
    CheckParity(field, generator, original, recovery, kShardBytes);

    FillInput(input.bytes() + 1, kShardBytes);
    SetPackedPointers(input.bytes() + 1, output.bytes() + 3, kShardBytes,
        original, recovery);
    std::memset(output.bytes() + 3, 0xa5,
        kRecoveryCount * kShardBytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute unaligned packed K32/R32 terminal");
    RequireTerminalCalls(1,
        "unaligned packed encode missed the T32 terminal");
    CheckParity(field, generator, original, recovery, kShardBytes);

    FillInput(input.bytes(), kShardBytes);
    SetPackedPointers(input.bytes(), output.bytes(), kShardBytes,
        original, recovery);
    RequireResult(leo2_test_codec_set_encode_mode(
        codec, LEO2_TEST_ENCODE_FORCE_TRANSFORM), LEO2_SUCCESS,
        "force mature K32/R32 transform");
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute forced mature K32/R32 transform");
    RequireTerminalCalls(0,
        "forced transform unexpectedly entered the T32 terminal");
    CheckParity(field, generator, original, recovery, kShardBytes);
    RequireResult(leo2_test_codec_set_encode_mode(
        codec, LEO2_TEST_ENCODE_AUTO), LEO2_SUCCESS,
        "restore K32/R32 AUTO encode mode");

    Require(leopard2_internal::
            SetBalancedB64TerminalEnabledForDiagnostics(false),
        "disable the K32/R32 terminal");
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute same-binary K32/R32 control route");
    RequireTerminalCalls(0,
        "same-binary control unexpectedly entered the T32 terminal");
    CheckParity(field, generator, original, recovery, kShardBytes);
    Require(leopard2_internal::
            SetBalancedB64TerminalEnabledForDiagnostics(true),
        "restore the K32/R32 terminal");

    std::memset(output.bytes(), 0xa5,
        kRecoveryCount * kShardBytes);
    recovery[11] = NULL;
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute sparse K32/R32 fallback");
    RequireTerminalCalls(0,
        "sparse output entered the dense T32 terminal");
    CheckParity(field, generator, original, recovery, kShardBytes);
    for (size_t i = 0; i < kShardBytes; ++i)
        Require(output.bytes()[11 * kShardBytes + i] == 0xa5,
            "sparse fallback modified a null parity shard");

    SetPackedPointers(input.bytes(), output.bytes(), kShardBytes,
        original, recovery);
    std::vector<uint8_t> detached_input(kShardBytes);
    std::memcpy(&detached_input[0], original[13], kShardBytes);
    original[13] = &detached_input[0];
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute detached-input K32/R32 fallback");
    RequireTerminalCalls(0,
        "detached input entered the packed T32 terminal");
    CheckParity(field, generator, original, recovery, kShardBytes);

    SetPackedPointers(input.bytes(), output.bytes(), kShardBytes,
        original, recovery);
    std::vector<uint8_t> detached_output(kShardBytes);
    recovery[19] = &detached_output[0];
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute detached-output K32/R32 fallback");
    RequireTerminalCalls(0,
        "detached output entered the packed T32 terminal");
    CheckParity(field, generator, original, recovery, kShardBytes);

    SetPackedPointers(input.bytes(), output.bytes(), kShardBytes,
        original, recovery);
    std::memset(output.bytes(), 0xa5,
        kRecoveryCount * kShardBytes);
    const std::vector<uint8_t> output_before(
        output.bytes(), output.bytes() + kRecoveryCount * kShardBytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, original, recovery,
        scratch.data(), scratch_bytes - 1), LEO2_SCRATCH_TOO_SMALL,
        "reject undersized K32/R32 terminal scratch");
    RequireTerminalCalls(0,
        "undersized scratch reached the T32 terminal");
    Require(std::memcmp(output.bytes(), &output_before[0],
            output_before.size()) == 0,
        "undersized scratch modified K32/R32 output");

    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, original, recovery,
        scratch.bytes() + 1, scratch_bytes), LEO2_BAD_ALIGNMENT,
        "reject misaligned K32/R32 terminal scratch");
    RequireTerminalCalls(0,
        "misaligned scratch reached the T32 terminal");
    Require(std::memcmp(output.bytes(), &output_before[0],
            output_before.size()) == 0,
        "misaligned scratch modified K32/R32 output");

    for (unsigned i = 0; i < kRecoveryCount; ++i)
        recovery[i] = input.bytes() + i * kShardBytes;
    const std::vector<uint8_t> input_before(
        input.bytes(), input.bytes() + kOriginalCount * kShardBytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_OVERLAP,
        "reject overlapping packed K32/R32 slabs");
    RequireTerminalCalls(0,
        "overlapping slabs reached the T32 terminal");
    Require(std::memcmp(input.bytes(), &input_before[0],
            input_before.size()) == 0,
        "overlap rejection modified K32/R32 input");

    SetPackedPointers(input.bytes(), output.bytes(), kShardBytes,
        original, recovery);
    alignas(64) uint8_t protected_storage[
        kRecoveryCount * kShardBytes] = {};
    leo2_encode_batch_item* const protected_item =
        new (protected_storage) leo2_encode_batch_item;
    protected_item->shard_bytes = kShardBytes;
    protected_item->original = original;
    protected_item->recovery = recovery;
    protected_item->scratch = scratch.data();
    protected_item->scratch_bytes = scratch_bytes;
    for (unsigned i = 0; i < kRecoveryCount; ++i)
        recovery[i] = protected_storage + i * kShardBytes;
    const leo2_encode_batch_item protected_before = *protected_item;
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch(codec, protected_item, 1),
        LEO2_OVERLAP, "reject K32/R32 output/batch-metadata overlap");
    RequireTerminalCalls(0,
        "metadata overlap reached the T32 terminal");
    Require(std::memcmp(protected_item, &protected_before,
            sizeof(*protected_item)) == 0,
        "metadata rejection modified the batch descriptor");

    SetPackedPointers(input.bytes(), output.bytes(), kShardBytes,
        original, recovery);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, NULL, recovery,
        scratch.data(), scratch_bytes), LEO2_INVALID_ARGUMENT,
        "reject null K32/R32 source array");
    RequireTerminalCalls(0,
        "null source array reached the T32 terminal");
    const void* const saved_original0 = original[0];
    original[0] = NULL;
    RequireResult(leo2_encode(codec, kShardBytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_INVALID_ARGUMENT,
        "reject null K32/R32 source pointer");
    original[0] = saved_original0;

    ExerciseNeighborBytes(codec, field, generator, 63);
    ExerciseNeighborBytes(codec, field, generator, 65);
    ExerciseNeighborBytes(codec, field, generator, 256);
    leo2_codec_destroy(codec);

    leo2_codec* const low_codec = CreateCodec(context, LEO2_PROFILE_LOW_V1);
    size_t low_scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        low_codec, kShardBytes, &low_scratch_bytes), LEO2_SUCCESS,
        "query low-profile K32/R32 scratch");
    AlignedBuffer low_scratch(low_scratch_bytes);
    SetPackedPointers(input.bytes(), output.bytes(), kShardBytes,
        original, recovery);
    const leopard2_test::ProfileLayout low_layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLow, kOriginalCount, kRecoveryCount);
    const leopard2_test::Matrix low_generator =
        leopard2_test::direct_systematic_generator(field, low_layout);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(low_codec, kShardBytes, original, recovery,
        low_scratch.data(), low_scratch.size()), LEO2_SUCCESS,
        "execute low-profile K32/R32 exclusion");
    RequireTerminalCalls(0,
        "low profile entered the legacy-high T32 terminal");
    CheckParity(field, low_generator, original, recovery, kShardBytes);
    leo2_codec_destroy(low_codec);
}

void ExerciseBackendShape(
    leo2_context* context,
    unsigned original_count,
    unsigned recovery_count)
{
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, original_count, recovery_count,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        LEO2_SUCCESS, "create backend-route packed codec");
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, kShardBytes, &scratch_bytes), LEO2_SUCCESS,
        "query backend-route packed scratch");
    AlignedBuffer scratch(scratch_bytes);
    AlignedBuffer input(
        static_cast<size_t>(original_count) * kShardBytes);
    AlignedBuffer output(
        static_cast<size_t>(recovery_count) * kShardBytes);
    FillVariableInput(input.bytes(), original_count, kShardBytes,
        UINT64_C(0x4241434b454e4452) ^ recovery_count);
    std::vector<const void*> original;
    std::vector<void*> recovery;
    SetVariablePackedPointers(input.bytes(), output.bytes(), original_count,
        recovery_count, kShardBytes, original, recovery);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes,
        &original[0], &recovery[0], scratch.data(), scratch.size()),
        LEO2_SUCCESS, "execute backend-route packed encode");
    bool terminal_shape_enabled = true;
#if !LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING || \
    defined(LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR)
    if (original_count == 62 && recovery_count == 8)
        terminal_shape_enabled = false;
#endif
    const bool k62_fused_shape =
        original_count == 62 && recovery_count == 8;
    const bool expected_avx2_route =
        leo2_context_backend(context) == LEO2_BACKEND_AVX2 &&
        terminal_shape_enabled;
    RequireTerminalCalls(
        expected_avx2_route && !k62_fused_shape ? 1U : 0U,
        "packed terminal selection disagrees with effective backend");
    const leopard::ff8::TestOnlyHighEncodeCounts counts =
        leopard::ff8::TestOnlyGetHighEncodeCounts();
    Require(counts.t8_k62_b64_fused_calls ==
            (expected_avx2_route && k62_fused_shape ? 1U : 0U),
        "K62/R8 fused-leaf selection disagrees with effective backend");

    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh, original_count, recovery_count);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);
    CheckVariableParity(field, generator, original_count, recovery_count,
        original, recovery, kShardBytes,
        "backend-route packed parity differs from oracle");
    leo2_codec_destroy(codec);
}

void ExerciseBackendRoute(leo2_backend requested_backend)
{
    leo2_context_options options = {};
    options.struct_size = sizeof(options);
    options.backend = requested_backend;
    options.thread_count = 1;
    leo2_context* context = NULL;
    const leo2_result context_result = leo2_context_create(&options, &context);
    if (context_result == LEO2_UNSUPPORTED)
        return;
    RequireResult(context_result, LEO2_SUCCESS,
        "create backend-route K32/R32 context");
    if (requested_backend != LEO2_BACKEND_AUTO)
    {
        leo2_codec* probe = NULL;
        RequireResult(leo2_codec_create(context, 65, 65,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &probe),
            LEO2_SUCCESS, "create explicit-backend K65/R65 probe");
        Require(!leopard2_internal::
                K65R65B64AVX512GFNISelectedForDiagnostics(probe, 64),
            "explicit backend request widened to AVX-512/GFNI T128");
        Require(!leopard2_internal::
                K65R65T128AVX512GFNILargerSelectedForDiagnostics(probe, 512),
            "explicit backend request widened to larger AVX-512/GFNI T128");
        leo2_codec_destroy(probe);
    }
    ExerciseBackendShape(context, 32, 32);
    ExerciseBackendShape(context, 62, 8);
    ExerciseBackendShape(context, 62, 16);
#if LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED
    ExerciseBackendShape(context, 66, 16);
#endif
    ExerciseBackendShape(context, 62, 32);
    ExerciseBackendShape(context, 79, 32);
    ExerciseBackendShape(context, 62, 33);
    ExerciseBackendShape(context, 65, 65);
    leo2_context_destroy(context);
}

void ExerciseK65R65AVX512GFNIAutoRoute()
{
    static const unsigned original_count = 65;
    static const unsigned recovery_count = 65;
    static const size_t shard_bytes = 64;

    leo2_context_options options = {};
    options.struct_size = sizeof(options);
    options.backend = LEO2_BACKEND_AUTO;
    options.thread_count = 1;
    leo2_context* context = NULL;
    RequireResult(leo2_context_create(&options, &context), LEO2_SUCCESS,
        "create AUTO K65/R65 AVX-512/GFNI context");
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, original_count, recovery_count,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        LEO2_SUCCESS, "create AUTO K65/R65 AVX-512/GFNI codec");

    if (!leopard2_internal::K65R65B64AVX512GFNISelectedForDiagnostics(
            codec, shard_bytes))
    {
        leo2_codec_destroy(codec);
        leo2_context_destroy(context);
        return;
    }

    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(codec, shard_bytes,
        &scratch_bytes), LEO2_SUCCESS,
        "query AUTO K65/R65 AVX-512/GFNI scratch");
    AlignedBuffer scratch(scratch_bytes);
    AlignedBuffer control_scratch(scratch_bytes);
    AlignedBuffer input_allocation(original_count * shard_bytes + 8U);
    AlignedBuffer output_allocation(recovery_count * shard_bytes + 8U);
    AlignedBuffer control_allocation(recovery_count * shard_bytes + 8U);
    uint8_t* const input = input_allocation.bytes() + 1U;
    uint8_t* const output = output_allocation.bytes() + 3U;
    uint8_t* const control = control_allocation.bytes() + 5U;
    std::vector<const void*> original;
    std::vector<void*> recovery;
    std::vector<void*> control_recovery;
    SetVariablePackedPointers(input, output, original_count, recovery_count,
        shard_bytes, original, recovery);
    std::vector<const void*> unused_original;
    SetVariablePackedPointers(input, control, original_count, recovery_count,
        shard_bytes, unused_original, control_recovery);

    FillVariableInput(input, original_count, shard_bytes,
        UINT64_C(0x4155544f36354746));
    Require(leopard2_internal::
            SetK65R65B64AVX512GFNIEnabledForDiagnostics(false),
        "disable K65/R65 AVX-512/GFNI leaf");
    RequireResult(leo2_encode(codec, shard_bytes, &original[0],
        &control_recovery[0], control_scratch.data(), control_scratch.size()),
        LEO2_SUCCESS, "execute mature K65/R65 same-binary control");
    Require(leopard2_internal::
            FinishK65R65B64AVX512GFNIRouteProbeForDiagnostics(),
        "finish disabled K65/R65 AVX-512/GFNI route probe");
    Require(!leopard2_internal::
            K65R65B64AVX512GFNISelectedForDiagnostics(codec, shard_bytes),
        "disabled K65/R65 AVX-512/GFNI leaf remained selected");
    Require(leopard2_internal::
            SetK65R65B64AVX512GFNIEnabledForDiagnostics(true),
        "restore K65/R65 AVX-512/GFNI leaf");
    Require(leopard2_internal::
            K65R65B64AVX512GFNISelectedForDiagnostics(codec, shard_bytes),
        "restored K65/R65 AVX-512/GFNI leaf was not selected");
    RequireResult(leo2_encode(codec, shard_bytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute K65/R65 AVX-512/GFNI leaf");
    Require(std::memcmp(output, control,
            recovery_count * shard_bytes) == 0,
        "K65/R65 AVX-512/GFNI parity differs from mature control");

    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(leopard2_test::kLegacyHigh,
            original_count, recovery_count);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);
    for (unsigned source = 0; source < original_count; ++source)
    {
        std::memset(input, 0, original_count * shard_bytes);
        for (size_t offset = 0; offset < shard_bytes; ++offset)
        {
            input[static_cast<size_t>(source) * shard_bytes + offset] =
                static_cast<uint8_t>(offset * 29U + source * 17U + 1U);
        }
        std::memset(output, 0xa5, recovery_count * shard_bytes);
        RequireResult(leo2_encode(codec, shard_bytes,
            &original[0], &recovery[0], scratch.data(), scratch.size()),
            LEO2_SUCCESS,
            "execute K65/R65 AVX-512/GFNI systematic basis encode");
        CheckVariableBasisParity(field, generator, original_count,
            recovery_count, source, original, recovery, shard_bytes,
            "K65/R65 AVX-512/GFNI basis differs from direct oracle");
    }

    FillVariableInput(input, original_count, shard_bytes,
        UINT64_C(0x4155544f42415443));
    leo2_encode_batch_item item = {
        shard_bytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()
    };
    RequireResult(leo2_encode_batch(codec, &item, 1), LEO2_SUCCESS,
        "execute K65/R65 AVX-512/GFNI one-item batch");
    CheckVariableParity(field, generator, original_count, recovery_count,
        original, recovery, shard_bytes,
        "K65/R65 AVX-512/GFNI batch differs from direct oracle");
    Require(leopard2_internal::
            FinishK65R65B64AVX512GFNIRouteProbeForDiagnostics(),
        "finish enabled K65/R65 AVX-512/GFNI route probe");

    leo2_codec_destroy(codec);
    leo2_context_destroy(context);
}

void ExerciseK65R65AVX512GFNILargerAutoRoute()
{
    static const unsigned original_count = 65;
    static const unsigned recovery_count = 65;
    static const size_t target_bytes[] = {
        128, 256, 512, 1024, 2048, 4096
    };
    static const size_t neighbor_bytes[] = {
        127, 129, 192, 768, 4095, 4097, 8192
    };

    leo2_context_options options = {};
    options.struct_size = sizeof(options);
    options.backend = LEO2_BACKEND_AUTO;
    options.thread_count = 1;
    leo2_context* context = NULL;
    RequireResult(leo2_context_create(&options, &context), LEO2_SUCCESS,
        "create AUTO K65/R65 larger AVX-512/GFNI context");
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, original_count, recovery_count,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        LEO2_SUCCESS, "create AUTO K65/R65 larger AVX-512/GFNI codec");

    if (!leopard2_internal::
            K65R65T128AVX512GFNILargerAvailableForDiagnostics(codec))
    {
        leo2_codec_destroy(codec);
        leo2_context_destroy(context);
        return;
    }

    for (size_t i = 0;
        i < sizeof(neighbor_bytes) / sizeof(neighbor_bytes[0]); ++i)
    {
        Require(!leopard2_internal::
                K65R65T128AVX512GFNILargerSelectedForDiagnostics(
                    codec, neighbor_bytes[i]),
            "inactive K65/R65 larger byte neighbor selected GFNI leaf");
    }

    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(leopard2_test::kLegacyHigh,
            original_count, recovery_count);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);

    for (size_t byte_index = 0;
        byte_index < sizeof(target_bytes) / sizeof(target_bytes[0]);
        ++byte_index)
    {
        const size_t shard_bytes = target_bytes[byte_index];
        Require(leopard2_internal::
                K65R65T128AVX512GFNILargerSelectedForDiagnostics(
                    codec, shard_bytes),
            "target K65/R65 larger byte cell missed GFNI leaf");

        size_t scratch_bytes = 0;
        RequireResult(leo2_encode_scratch_size(codec, shard_bytes,
            &scratch_bytes), LEO2_SUCCESS,
            "query larger K65/R65 AVX-512/GFNI scratch");
        AlignedBuffer scratch(scratch_bytes);
        AlignedBuffer control_scratch(scratch_bytes);
        AlignedBuffer input_allocation(
            original_count * shard_bytes + 8U);
        AlignedBuffer output_allocation(
            recovery_count * shard_bytes + 8U);
        AlignedBuffer control_allocation(
            recovery_count * shard_bytes + 8U);
        uint8_t* const input = input_allocation.bytes() + 1U;
        uint8_t* const output = output_allocation.bytes() + 3U;
        uint8_t* const control = control_allocation.bytes() + 5U;
        std::vector<const void*> original;
        std::vector<void*> recovery;
        std::vector<void*> control_recovery;
        SetVariablePackedPointers(input, output,
            original_count, recovery_count, shard_bytes,
            original, recovery);
        std::vector<const void*> unused_original;
        SetVariablePackedPointers(input, control,
            original_count, recovery_count, shard_bytes,
            unused_original, control_recovery);

        FillVariableInput(input, original_count, shard_bytes,
            UINT64_C(0x4c41524745523634) + shard_bytes);
        Require(leopard2_internal::
                SetK65R65T128AVX512GFNIEnabledForDiagnostics(false),
            "disable larger K65/R65 AVX-512/GFNI leaf");
        RequireResult(leo2_encode(codec, shard_bytes, &original[0],
            &control_recovery[0], control_scratch.data(),
            control_scratch.size()), LEO2_SUCCESS,
            "execute larger K65/R65 mature control");
        Require(leopard2_internal::
                K65R65T128AVX512GFNICallCountForDiagnostics() == 0U &&
                leopard2_internal::
                    K65R65T128AVX512GFNITileCountForDiagnostics() == 0U,
            "disabled larger K65/R65 GFNI route recorded a call");
        Require(leopard2_internal::
                FinishK65R65T128AVX512GFNIRouteProbeForDiagnostics(),
            "finish disabled larger K65/R65 GFNI route probe");

        Require(leopard2_internal::
                SetK65R65T128AVX512GFNIEnabledForDiagnostics(true),
            "enable larger K65/R65 AVX-512/GFNI leaf");
        RequireResult(leo2_encode(codec, shard_bytes, &original[0],
            &recovery[0], scratch.data(), scratch.size()), LEO2_SUCCESS,
            "execute larger K65/R65 AVX-512/GFNI leaf");
        Require(leopard2_internal::
                K65R65T128AVX512GFNICallCountForDiagnostics() == 1U &&
                leopard2_internal::
                    K65R65T128AVX512GFNITileCountForDiagnostics() ==
                        shard_bytes / 64U,
            "larger K65/R65 GFNI route counts differ");
        Require(std::memcmp(output, control,
                recovery_count * shard_bytes) == 0,
            "larger K65/R65 GFNI parity differs from mature control");
        CheckVariableParity(field, generator,
            original_count, recovery_count, original, recovery, shard_bytes,
            "larger K65/R65 GFNI parity differs from direct oracle");
        Require(leopard2_internal::
                FinishK65R65T128AVX512GFNIRouteProbeForDiagnostics(),
            "finish enabled larger K65/R65 GFNI route probe");

        FillVariableInput(input, original_count, shard_bytes,
            UINT64_C(0x4c41524745424154) + shard_bytes);
        leo2_encode_batch_item item = {
            shard_bytes, &original[0], &recovery[0],
            scratch.data(), scratch.size()
        };
        Require(leopard2_internal::
                SetK65R65T128AVX512GFNIEnabledForDiagnostics(true),
            "enable larger K65/R65 GFNI one-item batch probe");
        RequireResult(leo2_encode_batch(codec, &item, 1), LEO2_SUCCESS,
            "execute larger K65/R65 GFNI one-item batch");
        Require(leopard2_internal::
                K65R65T128AVX512GFNICallCountForDiagnostics() == 1U &&
                leopard2_internal::
                    K65R65T128AVX512GFNITileCountForDiagnostics() ==
                        shard_bytes / 64U,
            "larger K65/R65 GFNI batch route counts differ");
        CheckVariableParity(field, generator,
            original_count, recovery_count, original, recovery, shard_bytes,
            "larger K65/R65 GFNI batch differs from direct oracle");
        Require(leopard2_internal::
                FinishK65R65T128AVX512GFNIRouteProbeForDiagnostics(),
            "finish larger K65/R65 GFNI batch route probe");

        if (shard_bytes == 128)
        {
            for (unsigned source = 0; source < original_count; ++source)
            {
                std::memset(input, 0, original_count * shard_bytes);
                for (size_t offset = 0; offset < shard_bytes; ++offset)
                {
                    input[static_cast<size_t>(source) * shard_bytes +
                        offset] = static_cast<uint8_t>(
                            offset * 29U + source * 17U + 1U);
                }
                std::memset(output, 0xa5,
                    recovery_count * shard_bytes);
                RequireResult(leo2_encode(codec, shard_bytes,
                    &original[0], &recovery[0],
                    scratch.data(), scratch.size()), LEO2_SUCCESS,
                    "execute larger K65/R65 GFNI systematic basis");
                CheckVariableBasisParity(field, generator,
                    original_count, recovery_count, source,
                    original, recovery, shard_bytes,
                    "larger K65/R65 GFNI basis differs from direct oracle");
            }
        }
    }

    /*
        The larger leaf is intentionally limited to ordinary dense-packed
        one-item execution.  Exercise the nearby public shapes that must keep
        using the mature implementation while a route probe is armed.
    */
    {
        static const size_t shard_bytes = 512;
        size_t scratch_bytes = 0;
        RequireResult(leo2_encode_scratch_size(codec, shard_bytes,
            &scratch_bytes), LEO2_SUCCESS,
            "query larger GFNI fallback scratch");
        AlignedBuffer scratch(scratch_bytes);
        AlignedBuffer second_scratch(scratch_bytes);
        AlignedBuffer input(original_count * shard_bytes);
        AlignedBuffer output(recovery_count * shard_bytes);
        AlignedBuffer second_input(original_count * shard_bytes);
        AlignedBuffer second_output(recovery_count * shard_bytes);
        std::vector<const void*> original;
        std::vector<void*> recovery;
        std::vector<const void*> second_original;
        std::vector<void*> second_recovery;
        SetVariablePackedPointers(input.bytes(), output.bytes(),
            original_count, recovery_count, shard_bytes, original, recovery);
        SetVariablePackedPointers(second_input.bytes(), second_output.bytes(),
            original_count, recovery_count, shard_bytes,
            second_original, second_recovery);
        FillVariableInput(input.bytes(), original_count, shard_bytes,
            UINT64_C(0x4c4152474546414c));
        FillVariableInput(second_input.bytes(), original_count, shard_bytes,
            UINT64_C(0x4c4152474546414d));

        std::vector<void*> overlap_recovery(recovery_count);
        for (unsigned row = 0; row < recovery_count; ++row)
        {
            overlap_recovery[row] = input.bytes() +
                static_cast<size_t>(row) * shard_bytes;
        }
        Require(leopard2_internal::
                SetK65R65T128AVX512GFNIEnabledForDiagnostics(true),
            "arm larger GFNI overlap probe");
        RequireResult(leo2_encode(codec, shard_bytes,
            &original[0], &overlap_recovery[0],
            scratch.data(), scratch.size()), LEO2_OVERLAP,
            "reject larger GFNI input/output overlap");
        Require(leopard2_internal::
                K65R65T128AVX512GFNICallCountForDiagnostics() == 0U &&
                leopard2_internal::
                    K65R65T128AVX512GFNITileCountForDiagnostics() == 0U,
            "larger GFNI overlap entered candidate leaf");
        Require(leopard2_internal::
                FinishK65R65T128AVX512GFNIRouteProbeForDiagnostics(),
            "finish larger GFNI overlap probe");

        void* const saved_recovery0 = recovery[0];
        recovery[0] = NULL;
        std::memset(output.bytes(), 0xa5, output.size());
        Require(leopard2_internal::
                SetK65R65T128AVX512GFNIEnabledForDiagnostics(true),
            "arm larger GFNI partial-output probe");
        RequireResult(leo2_encode(codec, shard_bytes,
            &original[0], &recovery[0], scratch.data(), scratch.size()),
            LEO2_SUCCESS, "execute larger GFNI partial-output fallback");
        Require(leopard2_internal::
                K65R65T128AVX512GFNICallCountForDiagnostics() == 0U &&
                leopard2_internal::
                    K65R65T128AVX512GFNITileCountForDiagnostics() == 0U,
            "larger GFNI partial output entered candidate leaf");
        CheckVariableParity(field, generator,
            original_count, recovery_count, original, recovery, shard_bytes,
            "larger GFNI partial-output fallback differs from oracle");
        for (size_t offset = 0; offset < shard_bytes; ++offset)
        {
            Require(output.bytes()[offset] == 0xa5,
                "larger GFNI partial-output fallback wrote omitted parity");
        }
        Require(leopard2_internal::
                FinishK65R65T128AVX512GFNIRouteProbeForDiagnostics(),
            "finish larger GFNI partial-output probe");
        recovery[0] = saved_recovery0;

        leo2_encode_batch_item items[2] = {
            { shard_bytes, &original[0], &recovery[0],
                scratch.data(), scratch.size() },
            { shard_bytes, &second_original[0], &second_recovery[0],
                second_scratch.data(), second_scratch.size() }
        };
        Require(leopard2_internal::
                SetK65R65T128AVX512GFNIEnabledForDiagnostics(true),
            "arm larger GFNI multi-item probe");
        RequireResult(leo2_encode_batch(codec, items, 2), LEO2_SUCCESS,
            "execute larger GFNI multi-item mature batch");
        Require(leopard2_internal::
                K65R65T128AVX512GFNICallCountForDiagnostics() == 0U &&
                leopard2_internal::
                    K65R65T128AVX512GFNITileCountForDiagnostics() == 0U,
            "larger GFNI multi-item batch entered candidate leaf");
        CheckVariableParity(field, generator,
            original_count, recovery_count, original, recovery, shard_bytes,
            "first larger GFNI multi-item parity differs from oracle");
        CheckVariableParity(field, generator,
            original_count, recovery_count,
            second_original, second_recovery, shard_bytes,
            "second larger GFNI multi-item parity differs from oracle");
        Require(leopard2_internal::
                FinishK65R65T128AVX512GFNIRouteProbeForDiagnostics(),
            "finish larger GFNI multi-item probe");

        leo2_encode_batch_item binding_item = {
            shard_bytes, &original[0], &recovery[0],
            scratch.data(), scratch.size()
        };
        leo2_encode_batch_binding* binding = NULL;
        RequireResult(leo2_encode_batch_binding_create(
            codec, &binding_item, 1, &binding), LEO2_SUCCESS,
            "create larger GFNI reusable binding fallback");
        Require(leopard2_internal::
                SetK65R65T128AVX512GFNIEnabledForDiagnostics(true),
            "arm larger GFNI reusable-binding probe");
        RequireResult(leo2_encode_batch_binding_execute(binding),
            LEO2_SUCCESS,
            "execute larger GFNI reusable-binding mature fallback");
        Require(leopard2_internal::
                K65R65T128AVX512GFNICallCountForDiagnostics() == 0U &&
                leopard2_internal::
                    K65R65T128AVX512GFNITileCountForDiagnostics() == 0U,
            "larger GFNI reusable binding entered candidate leaf");
        CheckVariableParity(field, generator,
            original_count, recovery_count, original, recovery, shard_bytes,
            "larger GFNI reusable-binding parity differs from oracle");
        Require(leopard2_internal::
                FinishK65R65T128AVX512GFNIRouteProbeForDiagnostics(),
            "finish larger GFNI reusable-binding probe");
        leo2_encode_batch_binding_destroy(binding);

        const size_t row_stride = shard_bytes + 64U;
        AlignedBuffer strided_input(original_count * row_stride);
        AlignedBuffer strided_output(recovery_count * row_stride);
        std::vector<const void*> strided_original(original_count);
        std::vector<void*> strided_recovery(recovery_count);
        FillVariableInput(strided_input.bytes(), original_count, row_stride,
            UINT64_C(0x4c41524745535452));
        for (unsigned row = 0; row < original_count; ++row)
        {
            strided_original[row] = strided_input.bytes() +
                static_cast<size_t>(row) * row_stride;
        }
        for (unsigned row = 0; row < recovery_count; ++row)
        {
            strided_recovery[row] = strided_output.bytes() +
                static_cast<size_t>(row) * row_stride;
        }
        Require(leopard2_internal::
                SetK65R65T128AVX512GFNIEnabledForDiagnostics(true),
            "arm larger GFNI strided-layout probe");
        RequireResult(leo2_encode(codec, shard_bytes,
            &strided_original[0], &strided_recovery[0],
            scratch.data(), scratch.size()), LEO2_SUCCESS,
            "execute larger GFNI strided-layout mature fallback");
        Require(leopard2_internal::
                K65R65T128AVX512GFNICallCountForDiagnostics() == 0U &&
                leopard2_internal::
                    K65R65T128AVX512GFNITileCountForDiagnostics() == 0U,
            "larger GFNI strided layout entered candidate leaf");
        CheckVariableParity(field, generator,
            original_count, recovery_count,
            strided_original, strided_recovery, shard_bytes,
            "larger GFNI strided-layout fallback differs from oracle");
        Require(leopard2_internal::
                FinishK65R65T128AVX512GFNIRouteProbeForDiagnostics(),
            "finish larger GFNI strided-layout probe");

        std::memset(output.bytes(), 0xa5, output.size());
        Require(leopard2_internal::
                SetK65R65T128AVX512GFNIEnabledForDiagnostics(true),
            "arm larger GFNI undersized-scratch probe");
        RequireResult(leo2_encode(codec, shard_bytes,
            &original[0], &recovery[0], scratch.data(), scratch.size() - 1U),
            LEO2_SCRATCH_TOO_SMALL,
            "reject larger GFNI undersized scratch");
        Require(leopard2_internal::
                K65R65T128AVX512GFNICallCountForDiagnostics() == 0U &&
                leopard2_internal::
                    K65R65T128AVX512GFNITileCountForDiagnostics() == 0U,
            "larger GFNI invalid scratch entered candidate leaf");
        Require(leopard2_internal::
                FinishK65R65T128AVX512GFNIRouteProbeForDiagnostics(),
            "finish larger GFNI undersized-scratch probe");
    }

    static const unsigned shape_neighbors[][2] = {
        { 64, 65 }, { 66, 65 }, { 65, 64 }, { 65, 66 }
    };
    for (size_t i = 0;
        i < sizeof(shape_neighbors) / sizeof(shape_neighbors[0]); ++i)
    {
        leo2_codec* neighbor = NULL;
        RequireResult(leo2_codec_create(context,
            shape_neighbors[i][0], shape_neighbors[i][1],
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
            NULL, &neighbor), LEO2_SUCCESS,
            "create larger GFNI shape neighbor");
        Require(!leopard2_internal::
                K65R65T128AVX512GFNILargerSelectedForDiagnostics(
                    neighbor, 512),
            "larger GFNI shape neighbor selected candidate leaf");
        leo2_codec_destroy(neighbor);
    }

    leo2_codec_destroy(codec);
    leo2_context_destroy(context);
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
            std::printf(
                "balanced 64-byte AVX2 terminal test skipped: AVX2 unavailable\n");
            return 0;
        }
        RequireResult(context_result, LEO2_SUCCESS,
            "create AVX2 K32/R32 terminal context");
        ExerciseAVX2Terminal(context);
#if defined(LEO2_EXPERIMENT_HIGH_T16_B64_GENERATED)
        ExerciseSelectedPackedShape(context, 16, 16);
#else
        ExerciseExcludedShape(context, 16, 16);
#endif
        // The T16 aggregate terminal covers K=33..64 and R=9..16.  The
        // generated AVX2 circuit extends that band through K=65 and adds the
        // exact K=66/R=16 second-tail case.  Sweep
        // every supported K at the full T16 output width against the direct
        // systematic generator, then cover both Q3/Q4 block counts at the
        // smallest transmitted prefix.  ExerciseSelectedPackedShape also
        // proves that both one-shot and one-item batch calls select the route.
        const unsigned t16_maximum_original_count =
#if LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED
            65;
#else
            64;
#endif
        for (unsigned original_count = 33;
             original_count <= t16_maximum_original_count; ++original_count)
        {
            ExerciseSelectedPackedShape(context, original_count, 16);
        }
        ExerciseSelectedPackedShape(context, 33, 9);
        ExerciseSelectedPackedShape(context, 62, 9);
        ExerciseSelectedPackedShape(context, 64, 9);
#if LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING && \
    !defined(LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_VECTOR)
        ExerciseSelectedPackedShape(context, 62, 8);
#else
        ExerciseExcludedShape(context, 62, 8);
#endif
        ExerciseK62MatureScope(context);
        ExerciseK79ProfileFieldExclusions(context);
#if LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED
        ExerciseSelectedPackedShape(context, 65, 9);
        ExerciseK65R9TailBasis(context);
        ExerciseSelectedPackedShape(context, 66, 16);
        ExerciseK66R16TailBasis(context);
#else
        ExerciseExcludedShape(context, 65, 9);
        ExerciseExcludedShape(context, 66, 16);
#endif
        ExerciseSelectedPackedShape(context, 33, 17);
        ExerciseSelectedPackedShape(context, 62, 17);
        ExerciseSelectedPackedShape(context, 62, 30);
        ExerciseSelectedPackedShape(context, 33, 31);
        ExerciseSelectedPackedShape(context, 62, 31);
        // Every K changes the ragged inverse-transform suffix.  Sweep the
        // complete two- and three-block T32/R32 band, and retain the one-block
        // T64/R33 sweep separately.
        for (unsigned original_count = 33; original_count <= 96;
             ++original_count)
        {
            ExerciseSelectedPackedShape(context, original_count, 32);
            if (original_count <= 64)
                ExerciseSelectedPackedShape(context, original_count, 33);
        }
        ExerciseSelectedPackedShape(context, 33, 34);
        ExerciseSelectedPackedShape(context, 62, 34);
        ExerciseSelectedPackedShape(context, 64, 34);
        ExerciseSelectedPackedShape(context, 62, 35);
        ExerciseSelectedPackedShape(context, 33, 64);
        ExerciseSelectedPackedShape(context, 62, 64);
        ExerciseSelectedPackedShape(context, 63, 64);
        ExerciseSelectedPackedShape(context, 64, 63);
        ExerciseSelectedPackedShape(context, 64, 64);
        ExerciseK65R65PackedTerminal(context);
        ExerciseK65R65AVX512GFNIAutoRoute();
        ExerciseK65R65AVX512GFNILargerAutoRoute();
        ExerciseSelectedPackedShape(context, 128, 128);

        // The selector is deliberately exact.  Exercise both immediate count
        // neighbors at each promoted side and unrelated earlier tiny
        // terminals as layout controls.
        ExerciseExcludedShape(context, 15, 16);
        ExerciseExcludedShape(context, 16, 15);
        ExerciseExcludedShape(context, 31, 32);
        ExerciseExcludedShape(context, 32, 31);
        ExerciseExcludedShape(context, 79, 31);
        ExerciseExcludedShape(context, 79, 33);
        ExerciseExcludedShape(context, 97, 32);
        ExerciseExcludedShape(context, 32, 33);
        ExerciseExcludedShape(context, 65, 33);
        ExerciseExcludedShape(context, 65, 34);
        ExerciseExcludedShape(context, 61, 8);
        ExerciseExcludedShape(context, 63, 8);
        ExerciseExcludedShape(context, 62, 7);
        ExerciseExcludedShape(context, 66, 15);
        ExerciseExcludedShape(context, 66, 17);
        ExerciseExcludedShape(context, 67, 16);
        ExerciseExcludedShape(context, 62, 65);
        ExerciseExcludedShape(context, 64, 65);
        ExerciseExcludedShape(context, 66, 65);
        ExerciseExcludedShape(context, 65, 64);
        ExerciseExcludedShape(context, 65, 66);
        ExerciseExcludedShape(context, 127, 128);
        ExerciseExcludedShape(context, 128, 127);
        ExerciseExcludedShape(context, 8, 8);
        ExerciseExcludedShape(context, 5, 5);
        ExerciseExcludedShape(context, 16, 8);
        leo2_context_destroy(context);
        ExerciseBackendRoute(LEO2_BACKEND_AUTO);
        ExerciseBackendRoute(LEO2_BACKEND_SCALAR);
        ExerciseBackendRoute(LEO2_BACKEND_SSSE3);
        ExerciseBackendRoute(LEO2_BACKEND_AVX512);
        ExerciseBackendRoute(LEO2_BACKEND_GFNI);
        std::printf("balanced 64-byte packed AVX2 terminal checks passed\n");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "balanced 64-byte terminal failure: %s\n",
            error.what());
        return 1;
    }
}
