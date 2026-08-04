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

static const unsigned kOriginalCount = 4;
static const unsigned kRecoveryCount = 2;
static const uint8_t kGuardValue = 0x6d;
static const size_t kGuardBytes = 64;

void Require(bool condition, size_t shard_bytes, const char* message)
{
    if (!condition)
    {
        char prefix[48];
        std::snprintf(prefix, sizeof(prefix), "K=4/R=2/B=%zu: ",
            shard_bytes);
        throw std::runtime_error(std::string(prefix) + message);
    }
}

void RequireResult(
    leo2_result actual,
    leo2_result expected,
    size_t shard_bytes,
    const char* message)
{
    if (actual != expected)
    {
        throw std::runtime_error(std::string(message) + " at B=" +
            std::to_string(shard_bytes) + ": expected " +
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
    uint64_t state = UINT64_C(0x54324b3450524f44) ^ shard_bytes;
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

const leopard2_test::Matrix& Generator()
{
    static const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    static const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh, kOriginalCount, kRecoveryCount);
    static const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);
    return generator;
}

void CheckParity(
    size_t shard_bytes,
    const void* const* original,
    void* const* recovery)
{
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::Matrix& generator = Generator();
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
            Require(output[offset] == expected, shard_bytes,
                "production parity differs from direct generator");
        }
    }
}

leo2_codec* CreateCodec(leo2_context* context, size_t shard_bytes)
{
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, kOriginalCount, kRecoveryCount,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        LEO2_SUCCESS, shard_bytes, "create production K4/R2 codec");
    return codec;
}

void ExerciseCell(leo2_context* context, size_t shard_bytes)
{
    leo2_codec* codec = CreateCodec(context, shard_bytes);
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, shard_bytes, &scratch_bytes),
        LEO2_SUCCESS, shard_bytes, "query production scratch");
    Require(scratch_bytes != 0, shard_bytes,
        "production scratch query returned zero");
    AlignedBuffer scratch(scratch_bytes + leo2_scratch_alignment());

    const size_t input_payload = kOriginalCount * shard_bytes;
    const size_t output_payload = kRecoveryCount * shard_bytes;
    AlignedBuffer input(input_payload + 2U * kGuardBytes + 8U);
    AlignedBuffer output(output_payload + 2U * kGuardBytes + 8U);
    std::memset(input.bytes(), kGuardValue, input.size());
    std::memset(output.bytes(), kGuardValue, output.size());
    uint8_t* const input_base = input.bytes() + kGuardBytes + 1U;
    uint8_t* const output_base = output.bytes() + kGuardBytes + 3U;
    FillInput(input_base, shard_bytes);
    const std::vector<uint8_t> input_before(
        input.bytes(), input.bytes() + input.size());
    const void* original[kOriginalCount] = { NULL, NULL, NULL, NULL };
    void* recovery[kRecoveryCount] = { NULL, NULL };
    SetPackedPointers(
        input_base, output_base, shard_bytes, original, recovery);

    RequireResult(leo2_encode(codec, shard_bytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS, shard_bytes,
        "execute production K4/R2 encode");
    Require(std::memcmp(input.bytes(), &input_before[0],
            input_before.size()) == 0,
        shard_bytes, "production encode modified input or guards");
    CheckParity(shard_bytes, original, recovery);

    std::memset(output_base, kGuardValue, output_payload);
    leo2_encode_batch_item item = {
        shard_bytes, original, recovery, scratch.data(), scratch_bytes
    };
    RequireResult(leo2_encode_batch(codec, &item, 1),
        LEO2_SUCCESS, shard_bytes,
        "execute production K4/R2 one-item batch");
    CheckParity(shard_bytes, original, recovery);

    std::memset(output_base, kGuardValue, output_payload);
    leo2_encode_batch_binding* binding = NULL;
    RequireResult(leo2_encode_batch_binding_create(
        codec, &item, 1, &binding), LEO2_SUCCESS, shard_bytes,
        "create production K4/R2 encode binding");
    RequireResult(leo2_encode_batch_binding_execute(binding),
        LEO2_SUCCESS, shard_bytes,
        "execute production K4/R2 encode binding");
    CheckParity(shard_bytes, original, recovery);
    leo2_encode_batch_binding_destroy(binding);

    for (size_t i = 0; i < kGuardBytes + 3U; ++i)
        Require(output.bytes()[i] == kGuardValue, shard_bytes,
            "production encode modified leading output guard");
    for (size_t i = kGuardBytes + 3U + output_payload;
         i < output.size(); ++i)
    {
        Require(output.bytes()[i] == kGuardValue, shard_bytes,
            "production encode modified trailing output guard");
    }
    leo2_codec_destroy(codec);
}

void ExerciseProductionContracts(leo2_context* context)
{
    const size_t shard_bytes = 65;
    leo2_codec* codec = CreateCodec(context, shard_bytes);
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(codec, shard_bytes,
        &scratch_bytes), LEO2_SUCCESS, shard_bytes,
        "query production contract scratch");
    AlignedBuffer scratch(scratch_bytes + leo2_scratch_alignment());
    AlignedBuffer input(kOriginalCount * shard_bytes);
    AlignedBuffer output(kRecoveryCount * shard_bytes);
    FillInput(input.bytes(), shard_bytes);
    const void* original[kOriginalCount] = { NULL, NULL, NULL, NULL };
    void* recovery[kRecoveryCount] = { NULL, NULL };
    SetPackedPointers(input.bytes(), output.bytes(), shard_bytes,
        original, recovery);

    std::memset(output.bytes(), kGuardValue, output.size());
    const std::vector<uint8_t> output_before(
        output.bytes(), output.bytes() + output.size());
    RequireResult(leo2_encode(codec, shard_bytes, original, recovery,
        scratch.data(), scratch_bytes - 1U), LEO2_SCRATCH_TOO_SMALL,
        shard_bytes, "reject undersized production scratch");
    Require(std::memcmp(output.bytes(), &output_before[0],
            output_before.size()) == 0,
        shard_bytes, "undersized scratch modified output");

    RequireResult(leo2_encode(codec, shard_bytes, original, recovery,
        scratch.bytes() + 1U, scratch_bytes), LEO2_BAD_ALIGNMENT,
        shard_bytes, "reject misaligned production scratch");
    Require(std::memcmp(output.bytes(), &output_before[0],
            output_before.size()) == 0,
        shard_bytes, "misaligned scratch modified output");

    recovery[0] = input.bytes();
    recovery[1] = input.bytes() + shard_bytes;
    const std::vector<uint8_t> input_before(
        input.bytes(), input.bytes() + input.size());
    RequireResult(leo2_encode(codec, shard_bytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_OVERLAP,
        shard_bytes, "reject production input/output overlap");
    Require(std::memcmp(input.bytes(), &input_before[0],
            input_before.size()) == 0,
        shard_bytes, "overlap rejection modified input");

    SetPackedPointers(input.bytes(), output.bytes(), shard_bytes,
        original, recovery);
    std::vector<uint8_t> detached(shard_bytes);
    std::memcpy(&detached[0], original[3], shard_bytes);
    original[3] = &detached[0];
    std::memset(output.bytes(), kGuardValue, output.size());
    RequireResult(leo2_encode(codec, shard_bytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS, shard_bytes,
        "execute production nonpacked fallback");
    CheckParity(shard_bytes, original, recovery);

    leo2_codec_destroy(codec);
}

void ExerciseThreadedBatch(leo2_context* context)
{
    static const size_t kItems = 2;
    const size_t shard_bytes = 65;
    leo2_codec* codec = CreateCodec(context, shard_bytes);
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(codec, shard_bytes,
        &scratch_bytes), LEO2_SUCCESS, shard_bytes,
        "query threaded production scratch");
    const size_t alignment = leo2_scratch_alignment();
    const size_t scratch_stride =
        (scratch_bytes + alignment - 1U) & ~(alignment - 1U);
    AlignedBuffer scratch(kItems * scratch_stride);
    AlignedBuffer input(kItems * kOriginalCount * shard_bytes);
    AlignedBuffer output(kItems * kRecoveryCount * shard_bytes);
    const void* original[kItems][kOriginalCount] = {};
    void* recovery[kItems][kRecoveryCount] = {};
    leo2_encode_batch_item items[kItems] = {};

    for (size_t item = 0; item < kItems; ++item)
    {
        uint8_t* const input_base = input.bytes() +
            item * kOriginalCount * shard_bytes;
        uint8_t* const output_base = output.bytes() +
            item * kRecoveryCount * shard_bytes;
        FillInput(input_base, shard_bytes);
        for (size_t byte = 0; byte < kOriginalCount * shard_bytes; ++byte)
            input_base[byte] ^= static_cast<uint8_t>(item * 0x5bU);
        SetPackedPointers(input_base, output_base, shard_bytes,
            original[item], recovery[item]);
        items[item].shard_bytes = shard_bytes;
        items[item].original = original[item];
        items[item].recovery = recovery[item];
        items[item].scratch = scratch.bytes() + item * scratch_stride;
        items[item].scratch_bytes = scratch_bytes;
    }

    RequireResult(leo2_encode_batch(codec, items, kItems),
        LEO2_SUCCESS, shard_bytes,
        "execute threaded production K4/R2 batch");
    for (size_t item = 0; item < kItems; ++item)
        CheckParity(shard_bytes, original[item], recovery[item]);

    /* Whole-batch preflight must reject an item-0 write into item-1 input
       before either worker starts. */
    const std::vector<uint8_t> input_before(
        input.bytes(), input.bytes() + input.size());
    const std::vector<uint8_t> output_before(
        output.bytes(), output.bytes() + output.size());
    void* const saved_recovery = recovery[0][0];
    recovery[0][0] = const_cast<void*>(original[1][0]);
    RequireResult(leo2_encode_batch(codec, items, kItems),
        LEO2_OVERLAP, shard_bytes,
        "reject threaded output/other-item-input alias");
    recovery[0][0] = saved_recovery;
    Require(std::memcmp(input.bytes(), &input_before[0],
            input_before.size()) == 0,
        shard_bytes, "threaded alias rejection modified input");
    Require(std::memcmp(output.bytes(), &output_before[0],
            output_before.size()) == 0,
        shard_bytes, "threaded alias rejection modified output");

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
        options.thread_count = 2;
        leo2_context* context = NULL;
        const leo2_result result = leo2_context_create(&options, &context);
        if (result == LEO2_UNSUPPORTED)
        {
            std::printf(
                "production K4/R2 AVX2 test skipped: AVX2 unavailable\n");
            return 0;
        }
        if (result != LEO2_SUCCESS)
            throw std::runtime_error("create production AVX2 context");
        Require(leo2_context_thread_count(context) == 2U, 0,
            "production context did not retain two worker threads");

        static const size_t sizes[] = {
            1, 2, 3, 7, 15, 16, 17, 31, 32, 33, 63, 64, 65,
            127, 128, 129, 255, 256, 257, 1023, 1024, 1025,
            1984, 2016, 2047, 2048, 2049, 2111, 2112, 2113,
            4095, 4096, 4097, 8192, 16384, 16385
        };
        for (size_t i = 0; i < sizeof(sizes) / sizeof(sizes[0]); ++i)
            ExerciseCell(context, sizes[i]);
        ExerciseProductionContracts(context);
        ExerciseThreadedBatch(context);
        leo2_context_destroy(context);
        std::printf("production K4/R2 AVX2 checks passed\n");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "production K4/R2 failure: %s\n",
            error.what());
        return 1;
    }
}
