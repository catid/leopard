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
    RequireTerminalCalls(1,
        "selected shape missed the packed terminal");
    CheckVariableParity(field, generator, original_count, side,
        original, recovery, kShardBytes,
        "selected packed parity differs from direct oracle");

    std::memset(output.bytes(), 0xa5,
        static_cast<size_t>(side) * kShardBytes);
    leo2_encode_batch_item item = {
        kShardBytes, &original[0], &recovery[0], scratch.data(), scratch.size()
    };
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch(codec, &item, 1), LEO2_SUCCESS,
        "execute selected balanced one-item batch terminal");
    RequireTerminalCalls(1,
        "selected batch missed the packed terminal");
    CheckVariableParity(field, generator, original_count, side,
        original, recovery, kShardBytes,
        "selected packed batch parity differs from oracle");

    if (original_count == 33 && side == 32)
    {
        /* A detached final source must retain the mature path: the terminal
           is a packed-layout optimization, not a new API restriction. */
        std::vector<uint8_t> detached(kShardBytes);
        std::memcpy(&detached[0], original[32], kShardBytes);
        original[32] = &detached[0];
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

        original[32] = input.bytes() + 32U * kShardBytes;
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

    if (original_count == 62 && side == 32)
    {
        Require(leopard2_internal::
                SetBalancedB64TerminalEnabledForDiagnostics(false),
            "disable the K62/R32 terminal");
        std::memset(output.bytes(), 0xa5,
            static_cast<size_t>(side) * kShardBytes);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode(codec, kShardBytes,
            &original[0], &recovery[0], scratch.data(), scratch.size()),
            LEO2_SUCCESS, "execute same-binary K62/R32 control route");
        RequireTerminalCalls(0,
            "same-binary K62/R32 control entered the packed terminal");
        CheckVariableParity(field, generator, original_count, side,
            original, recovery, kShardBytes,
            "same-binary K62/R32 control differs from oracle");
        Require(leopard2_internal::
                SetBalancedB64TerminalEnabledForDiagnostics(true),
            "restore the K62/R32 terminal");

        std::memset(output.bytes(), 0xa5,
            static_cast<size_t>(side) * kShardBytes);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode(codec, kShardBytes,
            &original[0], &recovery[0], scratch.data(), scratch.size()),
            LEO2_SUCCESS, "reexecute restored K62/R32 terminal");
        RequireTerminalCalls(1,
            "restored K62/R32 route missed the packed terminal");
        CheckVariableParity(field, generator, original_count, side,
            original, recovery, kShardBytes,
            "restored K62/R32 terminal differs from oracle");
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
    leo2_codec_destroy(codec);
}

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
    RequireTerminalCalls(0,
        "excluded balanced-neighbor shape entered packed terminal");
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
    leo2_codec* codec = CreateCodec(
        context, LEO2_PROFILE_LEGACY_HIGH_V1);
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, kShardBytes, &scratch_bytes), LEO2_SUCCESS,
        "query backend-route K32/R32 scratch");
    AlignedBuffer scratch(scratch_bytes);
    AlignedBuffer input(kOriginalCount * kShardBytes);
    AlignedBuffer output(kRecoveryCount * kShardBytes);
    FillInput(input.bytes(), kShardBytes);
    const void* original[kOriginalCount];
    void* recovery[kRecoveryCount];
    SetPackedPointers(input.bytes(), output.bytes(), kShardBytes,
        original, recovery);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, original, recovery,
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute backend-route K32/R32 encode");
    const uint64_t expected_calls =
        leo2_context_backend(context) == LEO2_BACKEND_AVX2 ? 1U : 0U;
    RequireTerminalCalls(expected_calls,
        "balanced terminal selection disagrees with effective backend");

    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh,
            kOriginalCount, kRecoveryCount);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);
    CheckParity(field, generator, original, recovery, kShardBytes);
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
        ExerciseSelectedPackedShape(context, 33, 32);
        ExerciseSelectedPackedShape(context, 34, 32);
        ExerciseSelectedPackedShape(context, 62, 32);
        ExerciseSelectedPackedShape(context, 63, 32);
        ExerciseSelectedPackedShape(context, 64, 32);
        ExerciseSelectedPackedShape(context, 64, 64);
        ExerciseSelectedPackedShape(context, 128, 128);

        // The selector is deliberately exact.  Exercise both immediate count
        // neighbors at each promoted side and unrelated earlier tiny
        // terminals as layout controls.
        ExerciseExcludedShape(context, 15, 16);
        ExerciseExcludedShape(context, 16, 15);
        ExerciseExcludedShape(context, 31, 32);
        ExerciseExcludedShape(context, 32, 31);
        ExerciseExcludedShape(context, 33, 31);
        ExerciseExcludedShape(context, 65, 32);
        ExerciseExcludedShape(context, 63, 64);
        ExerciseExcludedShape(context, 64, 63);
        ExerciseExcludedShape(context, 127, 128);
        ExerciseExcludedShape(context, 128, 127);
        ExerciseExcludedShape(context, 8, 8);
        ExerciseExcludedShape(context, 5, 5);
        ExerciseExcludedShape(context, 16, 8);
        leo2_context_destroy(context);
        ExerciseBackendRoute(LEO2_BACKEND_AUTO);
        ExerciseBackendRoute(LEO2_BACKEND_SCALAR);
        ExerciseBackendRoute(LEO2_BACKEND_SSSE3);
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
