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

#if !defined(LEO2_EXPECT_HIGH_DIRECT_PRODUCTION)
#error "production high-direct expectation must be explicit"
#endif

namespace {

typedef std::vector<uint8_t> Bytes;
typedef std::vector<Bytes> Shards;

void Require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void RequireResult(leo2_result result, const char* message)
{
    if (result != LEO2_SUCCESS)
        throw std::runtime_error(
            std::string(message) + ": " + leo2_result_string(result));
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

    void* data() const { return value_; }
    size_t size() const { return bytes_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);
    void* value_;
    size_t bytes_;
};

Shards MakeOriginal(unsigned k, size_t bytes)
{
    Shards original(k, Bytes(bytes, 0));
    uint64_t state = UINT64_C(0x484947484b325231);
    for (unsigned shard = 0; shard < k; ++shard)
        for (size_t offset = 0; offset < bytes; ++offset)
        {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            original[shard][offset] = static_cast<uint8_t>(
                state + shard * 29u + offset * 131u);
        }
    return original;
}

Bytes OracleParity(
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator,
    const Shards& original,
    unsigned recovery_index)
{
    const unsigned k = static_cast<unsigned>(original.size());
    Bytes parity(original[0].size(), 0);
    const std::vector<leopard2_test::Element>& row =
        generator[k + recovery_index];
    for (size_t offset = 0; offset < parity.size(); ++offset)
    {
        leopard2_test::Element value = 0;
        for (unsigned source = 0; source < k; ++source)
        {
            value = field.add(value, field.multiply(
                row[source], original[source][offset]));
        }
        parity[offset] = static_cast<uint8_t>(value);
    }
    return parity;
}

void EncodeAndCheck(
    const leo2_codec* codec,
    const Shards& original,
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator,
    unsigned recovery_index)
{
    const uint8_t sentinel = 0xa5;
    const size_t bytes = original[0].size();
    std::vector<const void*> input(original.size(), NULL);
    for (size_t i = 0; i < original.size(); ++i)
        input[i] = &original[i][0];
    Shards recovery(16, Bytes(bytes, sentinel));
    std::vector<void*> output(recovery.size(), NULL);
    output[recovery_index] = &recovery[recovery_index][0];

    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, bytes, &scratch_bytes), "production scratch query");
    AlignedBuffer scratch(scratch_bytes);
    RequireResult(leo2_encode(codec, bytes, &input[0], &output[0],
        scratch.data(), scratch.size()), "production AUTO encode");

    Require(recovery[recovery_index] == OracleParity(
        field, generator, original, recovery_index),
        "production AUTO parity differs from the independent oracle");
    for (unsigned i = 0; i < recovery.size(); ++i)
    {
        if (i == recovery_index)
            continue;
        for (size_t offset = 0; offset < bytes; ++offset)
            Require(recovery[i][offset] == sentinel,
                "production AUTO modified an unrequested parity shard");
    }
}

} // namespace

int main()
{
    try
    {
        const unsigned k = 2;
        const unsigned r = 16;
        const size_t bytes = 4096;

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
                "Leopard2 production high-direct smoke skipped: AVX2 unavailable\n");
            return 0;
        }
        RequireResult(context_result, "explicit AVX2 context");
        Require(context != NULL &&
            leo2_context_backend(context) == LEO2_BACKEND_AVX2,
            "explicit AVX2 context selected another backend");

        leo2_codec* codec = NULL;
        RequireResult(leo2_codec_create(context, k, r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
            "production high codec");

        leopard2_internal::CodecEncodePathInfo path = {};
        Require(leopard2_internal::GetCodecEncodePathInfo(
            codec, bytes, 1, &path),
            "production encode-path introspection rejected the campaign cell");
#if LEO2_EXPECT_HIGH_DIRECT_PRODUCTION
        Require(path.direct_generator_rows == r,
            "option-ON production codec omitted direct generator rows");
        Require(path.auto_direct_selected,
            "option-ON production AUTO did not select the direct encoder");
#else
        Require(path.direct_generator_rows == 0,
            "default-OFF production codec retained direct generator rows");
        Require(!path.auto_direct_selected,
            "default-OFF production AUTO escaped the transform encoder");
#endif

        const Shards original = MakeOriginal(k, bytes);
        const Shards original_before = original;
        const leopard2_test::BinaryField field =
            leopard2_test::make_legacy_gf8();
        const leopard2_test::ProfileLayout layout =
            leopard2_test::make_profile_layout(
                leopard2_test::kLegacyHigh, k, r);
        const leopard2_test::Matrix generator =
            leopard2_test::direct_systematic_generator(field, layout);
        EncodeAndCheck(codec, original, field, generator, 0);
        EncodeAndCheck(codec, original, field, generator, r - 1);
        Require(original == original_before,
            "production AUTO modified an original shard");

        leo2_codec_destroy(codec);
        leo2_context_destroy(context);
        std::printf(
            "Leopard2 production high-direct smoke passed: option=%s "
            "K=2 R=16 bytes=4096 Q=1 parity=0,15\n",
#if LEO2_EXPECT_HIGH_DIRECT_PRODUCTION
            "ON"
#else
            "OFF"
#endif
        );
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr,
            "Leopard2 production high-direct smoke failed: %s\n",
            error.what());
        return 1;
    }
}
