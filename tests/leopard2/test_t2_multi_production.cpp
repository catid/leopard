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

static const unsigned kRecoveryCount = 2;
static const size_t kGuardBytes = 64;
static const uint8_t kGuardValue = 0x6d;

void Require(bool condition, unsigned k, size_t bytes, const char* message)
{
    if (!condition)
    {
        throw std::runtime_error("K=" + std::to_string(k) + "/R=2/B=" +
            std::to_string(bytes) + ": " + message);
    }
}

void RequireResult(
    leo2_result actual,
    leo2_result expected,
    unsigned k,
    size_t bytes,
    const char* message)
{
    if (actual != expected)
    {
        throw std::runtime_error("K=" + std::to_string(k) + "/R=2/B=" +
            std::to_string(bytes) + ": " + message + ": expected " +
            leo2_result_string(expected) + ", received " +
            leo2_result_string(actual));
    }
}

class AlignedBuffer
{
public:
    explicit AlignedBuffer(size_t bytes) : value_(NULL), bytes_(bytes)
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

void FillInput(uint8_t* input, unsigned k, size_t bytes)
{
    uint64_t state = UINT64_C(0x54324d554c544950) ^
        (static_cast<uint64_t>(k) << 32) ^ bytes;
    for (size_t i = 0; i < static_cast<size_t>(k) * bytes; ++i)
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
    unsigned k,
    size_t bytes,
    const void** original,
    void** recovery)
{
    for (unsigned i = 0; i < k; ++i)
        original[i] = input + static_cast<size_t>(i) * bytes;
    for (unsigned i = 0; i < kRecoveryCount; ++i)
        recovery[i] = output + static_cast<size_t>(i) * bytes;
}

void CheckParity(
    unsigned k,
    size_t bytes,
    const void* const* original,
    void* const* recovery,
    const leopard2_test::Matrix& generator)
{
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    for (unsigned parity = 0; parity < kRecoveryCount; ++parity)
    {
        const std::vector<leopard2_test::Element>& row =
            generator[k + parity];
        const uint8_t* output = static_cast<const uint8_t*>(recovery[parity]);
        for (size_t offset = 0; offset < bytes; ++offset)
        {
            leopard2_test::Element expected = 0;
            for (unsigned source = 0; source < k; ++source)
            {
                expected = field.add(expected, field.multiply(row[source],
                    static_cast<const uint8_t*>(original[source])[offset]));
            }
            Require(output[offset] == expected, k, bytes,
                "production parity differs from direct generator");
        }
    }
}

void ExerciseCell(leo2_context* context, unsigned k, size_t bytes)
{
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, k, kRecoveryCount,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        LEO2_SUCCESS, k, bytes, "create production codec");
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        LEO2_SUCCESS, k, bytes, "query production scratch");
    Require(scratch_bytes != 0, k, bytes,
        "production scratch query returned zero");

    AlignedBuffer scratch(scratch_bytes + leo2_scratch_alignment());
    const size_t input_payload = static_cast<size_t>(k) * bytes;
    const size_t output_payload = kRecoveryCount * bytes;
    AlignedBuffer input(input_payload + 2U * kGuardBytes + 8U);
    AlignedBuffer output(output_payload + 2U * kGuardBytes + 8U);
    std::memset(input.bytes(), kGuardValue, input.size());
    std::memset(output.bytes(), kGuardValue, output.size());
    uint8_t* const input_base = input.bytes() + kGuardBytes + 1U;
    uint8_t* const output_base = output.bytes() + kGuardBytes + 3U;
    FillInput(input_base, k, bytes);
    const std::vector<uint8_t> input_before(
        input.bytes(), input.bytes() + input.size());
    const void* original[16] = {};
    void* recovery[kRecoveryCount] = {};
    SetPackedPointers(input_base, output_base, k, bytes, original, recovery);

    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh, k, kRecoveryCount);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);

    RequireResult(leo2_encode(codec, bytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS, k, bytes,
        "execute production public encode");
    Require(std::memcmp(input.bytes(), &input_before[0], input.size()) == 0,
        k, bytes, "production encode modified input or guards");
    CheckParity(k, bytes, original, recovery, generator);

    std::memset(output_base, kGuardValue, output_payload);
    leo2_encode_batch_item item = {
        bytes, original, recovery, scratch.data(), scratch_bytes
    };
    RequireResult(leo2_encode_batch(codec, &item, 1), LEO2_SUCCESS,
        k, bytes, "execute production batch encode");
    CheckParity(k, bytes, original, recovery, generator);

    std::memset(output_base, kGuardValue, output_payload);
    leo2_encode_batch_binding* binding = NULL;
    RequireResult(leo2_encode_batch_binding_create(codec, &item, 1, &binding),
        LEO2_SUCCESS, k, bytes, "create production encode binding");
    RequireResult(leo2_encode_batch_binding_execute(binding), LEO2_SUCCESS,
        k, bytes, "execute production encode binding");
    CheckParity(k, bytes, original, recovery, generator);
    leo2_encode_batch_binding_destroy(binding);

    for (size_t i = 0; i < kGuardBytes + 3U; ++i)
        Require(output.bytes()[i] == kGuardValue, k, bytes,
            "production encode modified leading output guard");
    for (size_t i = kGuardBytes + 3U + output_payload;
         i < output.size(); ++i)
    {
        Require(output.bytes()[i] == kGuardValue, k, bytes,
            "production encode modified trailing output guard");
    }
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
            std::printf("production T=2 multi AVX2 test skipped\n");
            return 0;
        }
        if (result != LEO2_SUCCESS)
            throw std::runtime_error("create production AVX2 context");

        static const size_t sizes[] = {
            64, 128, 1024, 1984, 2048, 2112, 3072, 3136,
            4096, 65536, 65600
        };
        for (unsigned k = 5; k <= 16; ++k)
            for (size_t i = 0; i < sizeof(sizes) / sizeof(sizes[0]); ++i)
                ExerciseCell(context, k, sizes[i]);

        leo2_context_destroy(context);
        std::printf("production T=2 K=5..16 AVX2 checks passed\n");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "production T=2 multi failure: %s\n",
            error.what());
        return 1;
    }
}
