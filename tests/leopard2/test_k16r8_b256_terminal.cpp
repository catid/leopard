/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
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

static const unsigned kOriginalCount = 16;
static const unsigned kRecoveryCount = 8;
static const size_t kShardBytes = 256;

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
    uint64_t state = UINT64_C(0x4b31365242323536);
    for (size_t i = 0; i < kOriginalCount * kShardBytes; ++i)
    {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        input[i] = static_cast<uint8_t>(state >> 24);
    }
}

void SetPackedPointers(
    uint8_t* input_base,
    uint8_t* output_base,
    const void** original,
    void** recovery)
{
    for (unsigned i = 0; i < kOriginalCount; ++i)
        original[i] = input_base + i * kShardBytes;
    for (unsigned i = 0; i < kRecoveryCount; ++i)
        recovery[i] = output_base + i * kShardBytes;
}

void ResetOutput(uint8_t* output)
{
    std::memset(output, 0xa5, kRecoveryCount * kShardBytes);
}

void CheckParity(
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator,
    const void* const* original,
    void* const* recovery)
{
    for (unsigned parity = 0; parity < kRecoveryCount; ++parity)
    {
        if (!recovery[parity])
            continue;
        uint8_t* const output = static_cast<uint8_t*>(recovery[parity]);
        const std::vector<leopard2_test::Element>& row =
            generator[kOriginalCount + parity];
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
                "K16/R8/256 parity differs from the independent oracle");
        }
    }
}

void RequireTwoBlockCalls(uint64_t expected, const char* message)
{
    const leopard::ff8::TestOnlyHighEncodeCounts counts =
        leopard::ff8::TestOnlyGetHighEncodeCounts();
    Require(counts.two_block_calls == expected, message);
}

void ExerciseAVX2Terminal(leo2_context* context)
{
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context,
        kOriginalCount, kRecoveryCount,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        NULL, &codec), LEO2_SUCCESS, "create K16/R8 terminal codec");

    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, kShardBytes, &scratch_bytes), LEO2_SUCCESS,
        "query K16/R8 terminal scratch");
    Require(scratch_bytes >= 4736,
        "K16/R8 scratch query lost the established upper bound");
    AlignedBuffer scratch(scratch_bytes + leo2_scratch_alignment());
    AlignedBuffer input(kOriginalCount * kShardBytes + 8);
    AlignedBuffer output(kRecoveryCount * kShardBytes + 8);
    FillInput(input.bytes());

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
    SetPackedPointers(input.bytes(), output.bytes(), original, recovery);

    ResetOutput(output.bytes());
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute packed public K16/R8 terminal");
    RequireTwoBlockCalls(1,
        "packed public encode did not select the two-block terminal");
    CheckParity(field, generator, original, recovery);

    ResetOutput(output.bytes());
    leo2_encode_batch_item item = {
        kShardBytes, original, recovery, scratch.data(), scratch_bytes
    };
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch(codec, &item, 1), LEO2_SUCCESS,
        "execute packed one-item K16/R8 batch terminal");
    RequireTwoBlockCalls(1,
        "one-item batch did not select the two-block terminal");
    CheckParity(field, generator, original, recovery);

    SetPackedPointers(
        input.bytes() + 1, output.bytes() + 3, original, recovery);
    FillInput(input.bytes() + 1);
    ResetOutput(output.bytes() + 3);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute unaligned packed K16/R8 terminal");
    RequireTwoBlockCalls(1,
        "unaligned packed encode missed the two-block terminal");
    CheckParity(field, generator, original, recovery);

    SetPackedPointers(input.bytes(), output.bytes(), original, recovery);
    FillInput(input.bytes());
    ResetOutput(output.bytes());
    RequireResult(leo2_test_codec_set_encode_mode(
        codec, LEO2_TEST_ENCODE_FORCE_TRANSFORM), LEO2_SUCCESS,
        "force mature K16/R8 transform");
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute forced mature K16/R8 transform");
    RequireTwoBlockCalls(0,
        "forced transform unexpectedly entered the two-block terminal");
    CheckParity(field, generator, original, recovery);
    RequireResult(leo2_test_codec_set_encode_mode(
        codec, LEO2_TEST_ENCODE_AUTO), LEO2_SUCCESS,
        "restore K16/R8 AUTO encode mode");

    ResetOutput(output.bytes());
    recovery[3] = NULL;
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute sparse K16/R8 fallback");
    RequireTwoBlockCalls(0,
        "sparse output entered the dense two-block terminal");
    CheckParity(field, generator, original, recovery);
    for (size_t i = 0; i < kShardBytes; ++i)
        Require(output.bytes()[3 * kShardBytes + i] == 0xa5,
            "sparse fallback modified a null parity shard");

    SetPackedPointers(input.bytes(), output.bytes(), original, recovery);
    std::vector<uint8_t> detached(kShardBytes);
    std::memcpy(&detached[0], original[7], kShardBytes);
    original[7] = &detached[0];
    ResetOutput(output.bytes());
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute non-packed K16/R8 fallback");
    RequireTwoBlockCalls(0,
        "non-packed source entered the packed two-block terminal");
    CheckParity(field, generator, original, recovery);

    SetPackedPointers(input.bytes(), output.bytes(), original, recovery);
    ResetOutput(output.bytes());
    const std::vector<uint8_t> output_before(
        output.bytes(), output.bytes() + kRecoveryCount * kShardBytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, original, recovery,
        scratch.data(), scratch_bytes - 1), LEO2_SCRATCH_TOO_SMALL,
        "reject undersized K16/R8 terminal scratch");
    RequireTwoBlockCalls(0,
        "undersized scratch reached the two-block terminal");
    Require(std::memcmp(output.bytes(), &output_before[0],
            output_before.size()) == 0,
        "undersized scratch modified K16/R8 output");

    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, original, recovery,
        scratch.bytes() + 1, scratch_bytes), LEO2_BAD_ALIGNMENT,
        "reject misaligned K16/R8 terminal scratch");
    RequireTwoBlockCalls(0,
        "misaligned scratch reached the two-block terminal");

    for (unsigned i = 0; i < kRecoveryCount; ++i)
        recovery[i] = input.bytes() + i * kShardBytes;
    const std::vector<uint8_t> input_before(
        input.bytes(), input.bytes() + kOriginalCount * kShardBytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_OVERLAP,
        "reject overlapping packed K16/R8 slabs");
    RequireTwoBlockCalls(0,
        "overlapping slabs reached the two-block terminal");
    Require(std::memcmp(input.bytes(), &input_before[0],
            input_before.size()) == 0,
        "overlap rejection modified K16/R8 input");

    SetPackedPointers(input.bytes(), output.bytes(), original, recovery);
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
        LEO2_OVERLAP, "reject K16/R8 output/batch-metadata overlap");
    RequireTwoBlockCalls(0,
        "metadata overlap reached the two-block terminal");
    Require(std::memcmp(protected_item, &protected_before,
            sizeof(*protected_item)) == 0,
        "metadata rejection modified the batch descriptor");

    leo2_codec_destroy(codec);
}

void ExerciseScalarFallback()
{
    leo2_context_options options = {};
    options.struct_size = sizeof(options);
    options.backend = LEO2_BACKEND_SCALAR;
    options.thread_count = 1;
    leo2_context* context = NULL;
    RequireResult(leo2_context_create(&options, &context), LEO2_SUCCESS,
        "create scalar K16/R8 fallback context");
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context,
        kOriginalCount, kRecoveryCount,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        NULL, &codec), LEO2_SUCCESS, "create scalar K16/R8 codec");

    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, kShardBytes, &scratch_bytes), LEO2_SUCCESS,
        "query scalar K16/R8 scratch");
    AlignedBuffer scratch(scratch_bytes);
    AlignedBuffer input(kOriginalCount * kShardBytes);
    AlignedBuffer output(kRecoveryCount * kShardBytes);
    FillInput(input.bytes());
    const void* original[kOriginalCount];
    void* recovery[kRecoveryCount];
    SetPackedPointers(input.bytes(), output.bytes(), original, recovery);

    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, original, recovery,
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute scalar K16/R8 fallback");
    RequireTwoBlockCalls(0,
        "scalar context entered the AVX2 K16/R8 terminal");

    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh,
            kOriginalCount, kRecoveryCount);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);
    CheckParity(field, generator, original, recovery);
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
                "K16/R8/256 AVX2 terminal test skipped: AVX2 unavailable\n");
            return 0;
        }
        RequireResult(context_result, LEO2_SUCCESS,
            "create AVX2 K16/R8 terminal context");
        ExerciseAVX2Terminal(context);
        leo2_context_destroy(context);
        ExerciseScalarFallback();
        std::printf("K16/R8/256 packed AVX2 terminal checks passed\n");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "K16/R8/256 terminal failure: %s\n",
            error.what());
        return 1;
    }
}
