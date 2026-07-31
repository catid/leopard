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

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
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
#if !defined(LEO2_EXPECT_HIGH_DIRECT_AUTO)
#error "production high-direct AUTO expectation must be explicit"
#endif
#if !defined(LEO2_EXPECT_HIGH_T8_VECTOR)
#error "production high-T8 expectation must be explicit"
#endif
#if !defined(LEO2_EXPECT_HIGH_T8_PARTIAL_BINDING)
#error "production high-T8 partial-binding expectation must be explicit"
#endif
#if !defined(LEO2_EXPECT_HIGH_T8_TWO_BLOCK_BINDING)
#error "production high-T8 two-block-binding expectation must be explicit"
#endif
#ifndef LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_ONE_BLOCK_EXTENDED
#define LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_ONE_BLOCK_EXTENDED 0
#endif
#if LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_ONE_BLOCK_EXTENDED < 0 || \
    LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_ONE_BLOCK_EXTENDED > 1
#error "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_ONE_BLOCK_EXTENDED must be 0 or 1"
#endif
#ifndef LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_ONE_BLOCK_BEYOND_512
#define LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_ONE_BLOCK_BEYOND_512 0
#endif
#if LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_ONE_BLOCK_BEYOND_512 < 0 || \
    LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_ONE_BLOCK_BEYOND_512 > 1
#error "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_ONE_BLOCK_BEYOND_512 must be 0 or 1"
#endif
#ifndef LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_TWO_BLOCK_256
#define LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_TWO_BLOCK_256 0
#endif
#if LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_TWO_BLOCK_256 < 0 || \
    LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_TWO_BLOCK_256 > 1
#error "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_TWO_BLOCK_256 must be 0 or 1"
#endif
#ifndef LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_TWO_BLOCK_128_192
#define LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_TWO_BLOCK_128_192 0
#endif
#if LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_TWO_BLOCK_128_192 < 0 || \
    LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_TWO_BLOCK_128_192 > 1
#error "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_TWO_BLOCK_128_192 must be 0 or 1"
#endif
#ifndef LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_TWO_BLOCK_320
#define LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_TWO_BLOCK_320 0
#endif
#if LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_TWO_BLOCK_320 < 0 || \
    LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_TWO_BLOCK_320 > 1
#error "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_TWO_BLOCK_320 must be 0 or 1"
#endif
#ifndef LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_TWO_BLOCK_EXTENDED
#define LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_TWO_BLOCK_EXTENDED 0
#endif
#if LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_TWO_BLOCK_EXTENDED < 0 || \
    LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_TWO_BLOCK_EXTENDED > 1
#error "LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_TWO_BLOCK_EXTENDED must be 0 or 1"
#endif

namespace {

typedef std::vector<uint8_t> Bytes;
typedef std::vector<Bytes> Shards;

void Require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

bool IsExpectedT8OneBlockBeyond512ShapeByteCount(
    unsigned k,
    unsigned r,
    size_t bytes)
{
    static const uint16_t shape_masks[] = {
        UINT16_C(0xffff), UINT16_C(0xffff),
        UINT16_C(0x4ffd), UINT16_C(0x4fdd),
        UINT16_C(0x4fdd), UINT16_C(0x0fd4),
        UINT16_C(0x4efd), UINT16_C(0x4fcc)
    };
    if (LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_ONE_BLOCK_BEYOND_512 != 0 ||
        k < 5 || k > 8 || r < 5 || r > 8 ||
        bytes < 576 || bytes > 1088 || (bytes & 63U) != 0)
        return false;
    if (bytes == 1088)
        return k == 5 && r == 5;
    const size_t byte_index = (bytes - 576) / 64;
    const unsigned shape_bit = 4U * (k - 5U) + (r - 5U);
    return (shape_masks[byte_index] & (UINT16_C(1) << shape_bit)) != 0;
}

bool IsExpectedT8OneBlockByteCount(
    unsigned k,
    unsigned r,
    size_t bytes)
{
    return bytes == 64 ||
        (LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_ONE_BLOCK_EXTENDED == 0 &&
         bytes >= 128 && bytes <= 512 && (bytes & 63U) == 0) ||
        IsExpectedT8OneBlockBeyond512ShapeByteCount(k, r, bytes);
}

bool IsExpectedT8TwoBlockExtendedShapeByteCount(
    unsigned k,
    unsigned r,
    size_t bytes)
{
    static const uint32_t shape_masks[] = {
        UINT32_C(0xfffffff6), UINT32_C(0xffff5ff4),
        UINT32_C(0xffffeff0), UINT32_C(0xffff3ff0),
        UINT32_C(0xffff1ff0), UINT32_C(0xffff2f60),
        UINT32_C(0x6fff0e70), UINT32_C(0x5fff0d80),
        UINT32_C(0xffff0fd0), UINT32_C(0x5fff0d40),
        UINT32_C(0x6ff70c00)
    };
    if (LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_TWO_BLOCK_EXTENDED != 0 ||
        k < 9 || k > 16 || r < 5 || r > 8 ||
        bytes < 384 || bytes > 1024 || (bytes & 63U) != 0)
        return false;
    const size_t byte_index = (bytes - 384) / 64;
    const unsigned shape_bit = 4U * (k - 9U) + (r - 5U);
    return (shape_masks[byte_index] & (UINT32_C(1) << shape_bit)) != 0;
}

bool IsExpectedT8TwoBlockByteCount(
    unsigned k,
    unsigned r,
    size_t bytes)
{
    return bytes == 64 ||
        ((bytes == 128 || bytes == 192) &&
         LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_TWO_BLOCK_128_192 == 0) ||
        (bytes == 256 &&
         LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_TWO_BLOCK_256 == 0) ||
        (bytes == 320 &&
         LEO2_DIAGNOSTIC_DISABLE_HIGH_T8_TWO_BLOCK_320 == 0) ||
        IsExpectedT8TwoBlockExtendedShapeByteCount(k, r, bytes);
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

void EncodeAllAndCheck(
    const leo2_codec* codec,
    const Shards& original,
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator,
    unsigned recovery_count)
{
    const size_t bytes = original[0].size();
    std::vector<const void*> input(original.size(), NULL);
    for (size_t i = 0; i < original.size(); ++i)
        input[i] = &original[i][0];
    Shards recovery(recovery_count, Bytes(bytes, 0));
    std::vector<void*> output(recovery_count, NULL);
    for (unsigned i = 0; i < recovery_count; ++i)
        output[i] = &recovery[i][0];

    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, bytes, &scratch_bytes), "tiny full-output scratch query");
    AlignedBuffer scratch(scratch_bytes);
    RequireResult(leo2_encode(codec, bytes, &input[0], &output[0],
        scratch.data(), scratch.size()), "tiny full-output AUTO encode");

    for (unsigned recovery_index = 0;
         recovery_index < recovery_count;
         ++recovery_index)
    {
        Require(recovery[recovery_index] == OracleParity(
            field, generator, original, recovery_index),
            "tiny full-output AUTO parity differs from the independent oracle");
    }
}

uint64_t ExerciseT4BatchBindings(leo2_context* context)
{
    struct Shape
    {
        unsigned k;
        unsigned r;
        size_t maximum_bytes;
    };
    static const Shape shapes[] = {
        { 3, 3, 16384 }, { 3, 4, 2048 },
        { 4, 3, 12288 }, { 4, 4, 6144 },
        { 5, 3, 6144 }, { 5, 4, 4096 },
        { 6, 3, 16384 }, { 6, 4, 8192 },
        { 7, 3, 16384 }, { 7, 4, 4096 },
        { 9, 3, 8192 }, { 9, 4, 6144 },
        { 10, 3, 8192 }, { 10, 4, 6144 },
        { 11, 3, 6144 }, { 11, 4, 4096 }
    };
    static const size_t item_count = 2;
    static const uint8_t sentinel = 0xa5;
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    uint64_t checks = 0;

    for (size_t shape_i = 0;
         shape_i < sizeof(shapes) / sizeof(shapes[0]); ++shape_i)
    {
        const unsigned k = shapes[shape_i].k;
        const unsigned r = shapes[shape_i].r;
        const size_t maximum_bytes = shapes[shape_i].maximum_bytes;
        Require(leopard2_internal::HighT4BatchMaximumBytes(k, r) ==
                maximum_bytes,
            "T4 binding maximum-byte policy differs from test matrix");
        leo2_codec* codec = NULL;
        RequireResult(leo2_codec_create(context, k, r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
            "T4 binding codec");
        const leopard2_test::ProfileLayout layout =
            leopard2_test::make_profile_layout(
                leopard2_test::kLegacyHigh, k, r);
        const leopard2_test::Matrix generator =
            leopard2_test::direct_systematic_generator(field, layout);

        std::vector<size_t> byte_counts;
        byte_counts.push_back(32);
        byte_counts.push_back(64);
        byte_counts.push_back(96);
        byte_counts.push_back(2016);
        byte_counts.push_back(2048);
        byte_counts.push_back(2049);
        byte_counts.push_back(2080);
        byte_counts.push_back(3072);
        byte_counts.push_back(maximum_bytes - 32);
        byte_counts.push_back(maximum_bytes);
        byte_counts.push_back(maximum_bytes + 32);
        std::sort(byte_counts.begin(), byte_counts.end());
        byte_counts.erase(
            std::unique(byte_counts.begin(), byte_counts.end()),
            byte_counts.end());

        for (size_t byte_i = 0; byte_i < byte_counts.size(); ++byte_i)
        {
            const size_t bytes = byte_counts[byte_i];
            leopard2_internal::CodecEncodePathInfo path = {};
            Require(leopard2_internal::GetCodecEncodePathInfo(
                    codec, bytes, r, &path),
                "T4 binding path introspection");
            const bool expected_selected =
                bytes >= 32 && bytes <= maximum_bytes &&
                (bytes & 31U) == 0;
            Require(path.high_t4_batch_binding_selected ==
                    expected_selected,
                "T4 binding selector differs from its byte predicate");

            std::vector<Shards> original(
                item_count, MakeOriginal(k, bytes));
            original[1][0][0] ^= 0x6du;
            std::vector<Shards> recovery(
                item_count, Shards(r, Bytes(bytes, sentinel)));
            std::vector<std::vector<const void*> > input(
                item_count, std::vector<const void*>(k, NULL));
            std::vector<std::vector<void*> > output(
                item_count, std::vector<void*>(r, NULL));
            std::vector<std::unique_ptr<AlignedBuffer> > scratch;
            std::vector<leo2_encode_batch_item> items(item_count);
            size_t scratch_bytes = 0;
            RequireResult(leo2_encode_scratch_size(
                codec, bytes, &scratch_bytes),
                "T4 binding scratch query");
            for (size_t item = 0; item < item_count; ++item)
            {
                for (unsigned source = 0; source < k; ++source)
                    input[item][source] = &original[item][source][0];
                for (unsigned parity = 0; parity < r; ++parity)
                    output[item][parity] = &recovery[item][parity][0];
                scratch.push_back(std::unique_ptr<AlignedBuffer>(
                    new AlignedBuffer(scratch_bytes)));
                items[item].shard_bytes = bytes;
                items[item].original = &input[item][0];
                items[item].recovery = &output[item][0];
                items[item].scratch = scratch[item]->data();
                items[item].scratch_bytes = scratch[item]->size();
            }

            leo2_encode_batch_binding* binding = NULL;
            RequireResult(leo2_encode_batch_binding_create(
                codec, &items[0], items.size(), &binding),
                "T4 binding create");
            RequireResult(leo2_encode_batch_binding_execute(binding),
                "T4 binding execute");
            for (size_t item = 0; item < item_count; ++item)
                for (unsigned parity = 0; parity < r; ++parity)
                {
                    Require(recovery[item][parity] == OracleParity(
                        field, generator, original[item], parity),
                        "T4 binding parity differs from oracle");
                    ++checks;
                }

            if (bytes == 64 && k == 4 && r == 4)
            {
                original[1][3][17] ^= 0x93u;
                RequireResult(leo2_encode_batch_binding_execute(binding),
                    "T4 changed-source binding execute");
                for (unsigned parity = 0; parity < r; ++parity)
                {
                    Require(recovery[1][parity] == OracleParity(
                        field, generator, original[1], parity),
                        "T4 changed-source parity differs from oracle");
                    ++checks;
                }
                leo2_encode_batch_binding_destroy(binding);
                binding = NULL;

                std::fill(recovery[0][1].begin(),
                    recovery[0][1].end(), sentinel);
                output[0][1] = NULL;
                path = leopard2_internal::CodecEncodePathInfo();
                Require(leopard2_internal::GetCodecEncodePathInfo(
                        codec, bytes, r - 1, &path),
                    "sparse T4 path introspection");
                Require(!path.high_t4_batch_binding_selected,
                    "sparse T4 request selected dense batch shortcut");
                RequireResult(leo2_encode_batch_binding_create(
                    codec, &items[0], items.size(), &binding),
                    "sparse T4 binding create");
                RequireResult(leo2_encode_batch_binding_execute(binding),
                    "sparse T4 binding execute");
                Require(recovery[0][1] == Bytes(bytes, sentinel),
                    "sparse T4 binding modified a null output");
                ++checks;
            }
            leo2_encode_batch_binding_destroy(binding);
        }
        leo2_codec_destroy(codec);
    }
    return checks;
}

uint64_t ExerciseT8BatchBinding(
    leo2_context* context,
    size_t bytes)
{
    static const unsigned k = 8;
    static const unsigned r = 8;
    static const size_t batch_count = 8;
    static const uint8_t sentinel = 0xa5;

    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        "T8 binding codec");
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh, k, r);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);

    Shards original = MakeOriginal(k, bytes);
    std::vector<const void*> input(k, NULL);
    for (unsigned i = 0; i < k; ++i)
        input[i] = &original[i][0];
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, bytes, &scratch_bytes), "T8 binding scratch query");
    leopard2_internal::CodecEncodePathInfo t8_path = {};
    Require(leopard2_internal::GetCodecEncodePathInfo(
            codec, bytes, r, &t8_path),
        "T8 binding path introspection");
    Require(t8_path.high_t8_vector_selected ==
            (LEO2_EXPECT_HIGH_T8_VECTOR != 0 &&
             IsExpectedT8OneBlockByteCount(k, r, bytes)),
        "T8 binding selector differs from the production expectation");

    std::vector<Shards> recovery(
        batch_count, Shards(r, Bytes(bytes, sentinel)));
    std::vector<std::vector<void*> > output(
        batch_count, std::vector<void*>(r, NULL));
    std::vector<std::unique_ptr<AlignedBuffer> > scratch(batch_count);
    std::vector<leo2_encode_batch_item> items(batch_count);
    for (size_t batch = 0; batch < batch_count; ++batch)
    {
        for (unsigned parity = 0; parity < r; ++parity)
            output[batch][parity] = &recovery[batch][parity][0];
        scratch[batch].reset(new AlignedBuffer(scratch_bytes));
        items[batch].shard_bytes = bytes;
        items[batch].original = &input[0];
        items[batch].recovery = &output[batch][0];
        items[batch].scratch = scratch[batch]->data();
        items[batch].scratch_bytes = scratch[batch]->size();
    }

    leo2_encode_batch_binding* binding = NULL;
    RequireResult(leo2_encode_batch_binding_create(
        codec, &items[0], items.size(), &binding),
        "T8 binding create");
    Require(binding != NULL &&
            leo2_encode_batch_binding_item_count(binding) == batch_count,
        "T8 binding item count");
    RequireResult(leo2_encode_batch_binding_execute(binding),
        "T8 binding execute");

    uint64_t checks = 0;
    for (size_t batch = 0; batch < batch_count; ++batch)
        for (unsigned parity = 0; parity < r; ++parity)
        {
            Require(recovery[batch][parity] == OracleParity(
                field, generator, original, parity),
                "T8 binding parity differs from the independent oracle");
            ++checks;
        }

    /*
        Execution must read current bytes from captured addresses rather than
        caching a value-dependent transform during setup.
    */
    original[3][17] ^= 0x6du;
    RequireResult(leo2_encode_batch_binding_execute(binding),
        "T8 binding changed-source execute");
    for (size_t batch = 0; batch < batch_count; ++batch)
        for (unsigned parity = 0; parity < r; ++parity)
        {
            Require(recovery[batch][parity] == OracleParity(
                field, generator, original, parity),
                "T8 binding changed-source parity differs from oracle");
            ++checks;
        }
    original[3][17] ^= 0x6du;
    leo2_encode_batch_binding_destroy(binding);

    /*
        A sparse output mask must bypass the dense shortcut and retain the
        ordinary no-write guarantee for its null entry.
    */
    for (unsigned parity = 0; parity < r; ++parity)
        std::fill(recovery[0][parity].begin(),
            recovery[0][parity].end(), sentinel);
    output[0][3] = NULL;
    binding = NULL;
    RequireResult(leo2_encode_batch_binding_create(
        codec, &items[0], 1, &binding),
        "T8 sparse binding create");
    RequireResult(leo2_encode_batch_binding_execute(binding),
        "T8 sparse binding execute");
    for (unsigned parity = 0; parity < r; ++parity)
    {
        if (parity == 3)
        {
            Require(recovery[0][parity] == Bytes(bytes, sentinel),
                "T8 sparse binding modified a null output");
        }
        else
        {
            Require(recovery[0][parity] == OracleParity(
                field, generator, original, parity),
                "T8 sparse binding parity differs from oracle");
        }
        ++checks;
    }
    leo2_encode_batch_binding_destroy(binding);
    leo2_codec_destroy(codec);
    Require(checks == 136,
        "T8 binding check count changed unexpectedly");
    return checks;
}

uint64_t ExerciseT8PartialBindings(leo2_context* context)
{
    static const size_t byte_counts[] = {
        64, 128, 192, 256, 320, 384, 448, 512,
        576, 640, 704, 768, 832, 896, 960, 1024, 1088
    };
    static const uint8_t sentinel = 0xa5;
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    uint64_t checks = 0;

    for (size_t byte_index = 0;
         byte_index < sizeof(byte_counts) / sizeof(byte_counts[0]);
         ++byte_index)
    {
        const size_t bytes = byte_counts[byte_index];
        for (unsigned k = 5; k <= 8; ++k)
            for (unsigned r = 5; r <= 8; ++r)
            {
                leo2_codec* codec = NULL;
                RequireResult(leo2_codec_create(context, k, r,
                    LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
                    "partial T8 binding codec");
                const leopard2_test::ProfileLayout layout =
                    leopard2_test::make_profile_layout(
                        leopard2_test::kLegacyHigh, k, r);
                const leopard2_test::Matrix generator =
                    leopard2_test::direct_systematic_generator(field, layout);
                Shards original = MakeOriginal(k, bytes);
                std::vector<const void*> input(k, NULL);
                for (unsigned i = 0; i < k; ++i)
                    input[i] = &original[i][0];
                Shards recovery(r, Bytes(bytes, sentinel));
                std::vector<void*> output(r, NULL);
                for (unsigned parity = 0; parity < r; ++parity)
                    output[parity] = &recovery[parity][0];

                size_t scratch_bytes = 0;
                RequireResult(leo2_encode_scratch_size(
                    codec, bytes, &scratch_bytes),
                    "partial T8 binding scratch query");
                AlignedBuffer scratch(scratch_bytes);
                leopard2_internal::CodecEncodePathInfo path = {};
                Require(leopard2_internal::GetCodecEncodePathInfo(
                        codec, bytes, r, &path),
                    "partial T8 binding path introspection");
                const bool expected =
                    LEO2_EXPECT_HIGH_T8_PARTIAL_BINDING != 0 &&
                    IsExpectedT8OneBlockByteCount(k, r, bytes) &&
                    (k != 8 || r != 8);
                Require(path.high_t8_partial_binding_selected == expected,
                    "partial T8 binding selector differs from expectation");
                Require(path.high_t8_vector_selected ==
                        (LEO2_EXPECT_HIGH_T8_VECTOR != 0 &&
                         IsExpectedT8OneBlockByteCount(k, r, bytes) &&
                         k == 8 && r == 8),
                    "full T8 vector selector differs from expectation");

                leo2_encode_batch_item item = {};
                item.shard_bytes = bytes;
                item.original = &input[0];
                item.recovery = &output[0];
                item.scratch = scratch.data();
                item.scratch_bytes = scratch.size();
                leo2_encode_batch_binding* binding = NULL;
                RequireResult(leo2_encode_batch_binding_create(
                    codec, &item, 1, &binding),
                    "partial T8 binding create");
                RequireResult(leo2_encode_batch_binding_execute(binding),
                    "partial T8 binding execute");
                for (unsigned parity = 0; parity < r; ++parity)
                {
                    Require(recovery[parity] == OracleParity(
                        field, generator, original, parity),
                        "partial T8 binding parity differs from oracle");
                    ++checks;
                }

                if (k == 5 && r == 5)
                {
                    original[1][13] ^= 0x6du;
                    RequireResult(leo2_encode_batch_binding_execute(binding),
                        "partial T8 changed-source binding execute");
                    for (unsigned parity = 0; parity < r; ++parity)
                    {
                        Require(recovery[parity] == OracleParity(
                            field, generator, original, parity),
                            "partial T8 changed-source parity differs from oracle");
                        ++checks;
                    }
                    original[1][13] ^= 0x6du;
                    leo2_encode_batch_binding_destroy(binding);
                    binding = NULL;
                    for (unsigned parity = 0; parity < r; ++parity)
                        std::fill(recovery[parity].begin(),
                            recovery[parity].end(), sentinel);
                    output[2] = NULL;
                    path = leopard2_internal::CodecEncodePathInfo();
                    Require(leopard2_internal::GetCodecEncodePathInfo(
                            codec, bytes, r - 1, &path),
                        "sparse partial T8 path introspection");
                    Require(!path.high_t8_partial_binding_selected,
                        "sparse partial T8 binding selected dense shortcut");
                    binding = NULL;
                    RequireResult(leo2_encode_batch_binding_create(
                        codec, &item, 1, &binding),
                        "sparse partial T8 binding create");
                    RequireResult(leo2_encode_batch_binding_execute(binding),
                        "sparse partial T8 binding execute");
                    for (unsigned parity = 0; parity < r; ++parity)
                    {
                        if (parity == 2)
                        {
                            Require(recovery[parity] == Bytes(bytes, sentinel),
                                "sparse partial T8 modified a null output");
                        }
                        else
                        {
                            Require(recovery[parity] == OracleParity(
                                field, generator, original, parity),
                                "sparse partial T8 parity differs from oracle");
                        }
                        ++checks;
                    }
                    leo2_encode_batch_binding_destroy(binding);
                }
                else
                    leo2_encode_batch_binding_destroy(binding);
                leo2_codec_destroy(codec);
            }
    }

    Require(checks == 114 *
            (sizeof(byte_counts) / sizeof(byte_counts[0])),
        "partial T8 binding check count changed unexpectedly");
    return checks;
}

uint64_t ExerciseT8TwoBlockBindings(leo2_context* context)
{
    static const size_t byte_counts[] = {
        64, 128, 192, 256, 320, 384, 448, 512,
        576, 640, 704, 768, 832, 896, 960, 1024
    };
    static const uint8_t sentinel = 0xa5;
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    uint64_t checks = 0;

    for (size_t byte_index = 0;
         byte_index < sizeof(byte_counts) / sizeof(byte_counts[0]);
         ++byte_index)
    {
        const size_t bytes = byte_counts[byte_index];
        for (unsigned k = 9; k <= 16; ++k)
            for (unsigned r = 5; r <= 8; ++r)
            {
                leo2_codec* codec = NULL;
                RequireResult(leo2_codec_create(context, k, r,
                    LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
                    NULL, &codec),
                    "two-block T8 binding codec");
                const leopard2_test::ProfileLayout layout =
                    leopard2_test::make_profile_layout(
                        leopard2_test::kLegacyHigh, k, r);
                const leopard2_test::Matrix generator =
                    leopard2_test::direct_systematic_generator(
                        field, layout);
                Shards original = MakeOriginal(k, bytes);
                std::vector<const void*> input(k, NULL);
                for (unsigned i = 0; i < k; ++i)
                    input[i] = &original[i][0];
                Shards recovery(r, Bytes(bytes, sentinel));
                std::vector<void*> output(r, NULL);
                for (unsigned i = 0; i < r; ++i)
                    output[i] = &recovery[i][0];

                size_t scratch_bytes = 0;
                RequireResult(leo2_encode_scratch_size(
                    codec, bytes, &scratch_bytes),
                    "two-block T8 binding scratch query");
                AlignedBuffer scratch(scratch_bytes);
                leopard2_internal::CodecEncodePathInfo path = {};
                Require(leopard2_internal::GetCodecEncodePathInfo(
                        codec, bytes, r, &path),
                    "two-block T8 binding path introspection");
                const bool expect_two_block =
                    LEO2_EXPECT_HIGH_T8_TWO_BLOCK_BINDING != 0 &&
                    IsExpectedT8TwoBlockByteCount(k, r, bytes);
                Require(path.high_t8_two_block_binding_selected ==
                        expect_two_block,
                    "two-block T8 selector differs from expectation");
                Require(!path.high_t8_partial_binding_selected,
                    "two-block T8 shape selected the one-block binding");

                leo2_encode_batch_item item = {};
                item.shard_bytes = bytes;
                item.original = &input[0];
                item.recovery = &output[0];
                item.scratch = scratch.data();
                item.scratch_bytes = scratch.size();
                leo2_encode_batch_binding* binding = NULL;
                RequireResult(leo2_encode_batch_binding_create(
                    codec, &item, 1, &binding),
                    "two-block T8 binding create");
                RequireResult(leo2_encode_batch_binding_execute(binding),
                    "two-block T8 binding execute");
                for (unsigned parity = 0; parity < r; ++parity)
                {
                    Require(recovery[parity] == OracleParity(
                        field, generator, original, parity),
                        "two-block T8 parity differs from oracle");
                    ++checks;
                }

                const bool exercise_mutation_and_sparse =
                    (bytes < 384 && k == 9 && r == 5) ||
                    (bytes >= 384 && k == 13 && r == 5);
                if (exercise_mutation_and_sparse)
                {
                    original[8][bytes - 1] ^= 0x6du;
                    RequireResult(
                        leo2_encode_batch_binding_execute(binding),
                        "two-block T8 changed-source execute");
                    for (unsigned parity = 0; parity < r; ++parity)
                    {
                        Require(recovery[parity] == OracleParity(
                            field, generator, original, parity),
                            "two-block T8 changed-source parity differs from oracle");
                        ++checks;
                    }
                    original[8][bytes - 1] ^= 0x6du;
                    leo2_encode_batch_binding_destroy(binding);
                    binding = NULL;

                    for (unsigned parity = 0; parity < r; ++parity)
                        std::fill(recovery[parity].begin(),
                            recovery[parity].end(), sentinel);
                    output[2] = NULL;
                    path = leopard2_internal::CodecEncodePathInfo();
                    Require(leopard2_internal::GetCodecEncodePathInfo(
                            codec, bytes, r - 1, &path),
                        "sparse two-block T8 path introspection");
                    Require(!path.high_t8_two_block_binding_selected,
                        "sparse two-block T8 selected the dense shortcut");
                    RequireResult(leo2_encode_batch_binding_create(
                        codec, &item, 1, &binding),
                        "sparse two-block T8 binding create");
                    RequireResult(
                        leo2_encode_batch_binding_execute(binding),
                        "sparse two-block T8 binding execute");
                    for (unsigned parity = 0; parity < r; ++parity)
                    {
                        if (parity == 2)
                        {
                            Require(recovery[parity] ==
                                    Bytes(bytes, sentinel),
                                "sparse two-block T8 modified a null output");
                        }
                        else
                        {
                            Require(recovery[parity] == OracleParity(
                                field, generator, original, parity),
                                "sparse two-block T8 parity differs from oracle");
                        }
                        ++checks;
                    }
                }
                leo2_encode_batch_binding_destroy(binding);
                leo2_codec_destroy(codec);
            }
    }

    Require(checks == 218 *
            (sizeof(byte_counts) / sizeof(byte_counts[0])),
        "two-block T8 binding check count changed unexpectedly");
    return checks;
}

uint64_t ExerciseT8PartialThreadPool(size_t bytes)
{
    const unsigned k = bytes == 1088 ? 5 : (bytes > 512 ? 7 : 5);
    const unsigned r = bytes == 1088 ? 5 : (bytes > 512 ? 7 : 5);
    static const size_t batch_count = 8;

    leo2_context_options options = {};
    options.struct_size = sizeof(options);
    options.backend = LEO2_BACKEND_AVX2;
    options.thread_count = 4;
    leo2_context* context = NULL;
    RequireResult(leo2_context_create(&options, &context),
        "partial T8 thread-pool context");
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        "partial T8 thread-pool codec");
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh, k, r);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);
    Shards original = MakeOriginal(k, bytes);
    std::vector<const void*> input(k, NULL);
    for (unsigned i = 0; i < k; ++i)
        input[i] = &original[i][0];

    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, bytes, &scratch_bytes),
        "partial T8 thread-pool scratch query");
    std::vector<Shards> recovery(
        batch_count, Shards(r, Bytes(bytes, 0)));
    std::vector<std::vector<void*> > output(
        batch_count, std::vector<void*>(r, NULL));
    std::vector<std::unique_ptr<AlignedBuffer> > scratch(batch_count);
    std::vector<leo2_encode_batch_item> items(batch_count);
    for (size_t batch = 0; batch < batch_count; ++batch)
    {
        for (unsigned parity = 0; parity < r; ++parity)
            output[batch][parity] = &recovery[batch][parity][0];
        scratch[batch].reset(new AlignedBuffer(scratch_bytes));
        items[batch].shard_bytes = bytes;
        items[batch].original = &input[0];
        items[batch].recovery = &output[batch][0];
        items[batch].scratch = scratch[batch]->data();
        items[batch].scratch_bytes = scratch[batch]->size();
    }

    leo2_encode_batch_binding* binding = NULL;
    RequireResult(leo2_encode_batch_binding_create(
        codec, &items[0], items.size(), &binding),
        "partial T8 thread-pool binding create");
    RequireResult(leo2_encode_batch_binding_execute(binding),
        "partial T8 thread-pool binding execute");
    uint64_t checks = 0;
    for (size_t batch = 0; batch < batch_count; ++batch)
        for (unsigned parity = 0; parity < r; ++parity)
        {
            Require(recovery[batch][parity] == OracleParity(
                field, generator, original, parity),
                "partial T8 thread-pool parity differs from oracle");
            ++checks;
        }

    leo2_encode_batch_binding_destroy(binding);
    leo2_codec_destroy(codec);
    leo2_context_destroy(context);
    Require(checks == batch_count * r,
        "partial T8 thread-pool check count changed unexpectedly");
    return checks;
}

uint64_t ExerciseT8TwoBlockThreadPool(size_t bytes)
{
    static const unsigned k = 13;
    static const unsigned r = 5;
    static const size_t batch_count = 8;

    leo2_context_options options = {};
    options.struct_size = sizeof(options);
    options.backend = LEO2_BACKEND_AVX2;
    options.thread_count = 4;
    leo2_context* context = NULL;
    RequireResult(leo2_context_create(&options, &context),
        "two-block T8 thread-pool context");
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        "two-block T8 thread-pool codec");
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh, k, r);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);
    Shards original = MakeOriginal(k, bytes);
    std::vector<const void*> input(k, NULL);
    for (unsigned i = 0; i < k; ++i)
        input[i] = &original[i][0];

    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, bytes, &scratch_bytes),
        "two-block T8 thread-pool scratch query");
    std::vector<Shards> recovery(
        batch_count, Shards(r, Bytes(bytes, 0)));
    std::vector<std::vector<void*> > output(
        batch_count, std::vector<void*>(r, NULL));
    std::vector<std::unique_ptr<AlignedBuffer> > scratch(batch_count);
    std::vector<leo2_encode_batch_item> items(batch_count);
    for (size_t batch = 0; batch < batch_count; ++batch)
    {
        for (unsigned parity = 0; parity < r; ++parity)
            output[batch][parity] = &recovery[batch][parity][0];
        scratch[batch].reset(new AlignedBuffer(scratch_bytes));
        items[batch].shard_bytes = bytes;
        items[batch].original = &input[0];
        items[batch].recovery = &output[batch][0];
        items[batch].scratch = scratch[batch]->data();
        items[batch].scratch_bytes = scratch[batch]->size();
    }

    leo2_encode_batch_binding* binding = NULL;
    RequireResult(leo2_encode_batch_binding_create(
        codec, &items[0], items.size(), &binding),
        "two-block T8 thread-pool binding create");
    RequireResult(leo2_encode_batch_binding_execute(binding),
        "two-block T8 thread-pool binding execute");
    uint64_t checks = 0;
    for (size_t batch = 0; batch < batch_count; ++batch)
        for (unsigned parity = 0; parity < r; ++parity)
        {
            Require(recovery[batch][parity] == OracleParity(
                field, generator, original, parity),
                "two-block T8 thread-pool parity differs from oracle");
            ++checks;
        }

    leo2_encode_batch_binding_destroy(binding);
    leo2_codec_destroy(codec);
    leo2_context_destroy(context);
    Require(checks == batch_count * r,
        "two-block T8 thread-pool check count changed unexpectedly");
    return checks;
}

uint64_t ExerciseT8PartialUnaligned(
    leo2_context* context,
    size_t bytes)
{
    const unsigned k = bytes == 1088 ? 5 : (bytes > 512 ? 7 : 5);
    const unsigned r = bytes == 1088 ? 5 : (bytes > 512 ? 7 : 5);
    static const uint8_t sentinel = 0xa5;

    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        "partial T8 unaligned codec");
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh, k, r);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);
    const Shards original = MakeOriginal(k, bytes);
    Shards original_storage(k, Bytes(bytes + 2, sentinel));
    std::vector<const void*> input(k, NULL);
    for (unsigned source = 0; source < k; ++source)
    {
        std::copy(original[source].begin(), original[source].end(),
            original_storage[source].begin() + 1);
        input[source] = &original_storage[source][1];
    }
    Shards recovery_storage(r, Bytes(bytes + 2, sentinel));
    std::vector<void*> output(r, NULL);
    for (unsigned parity = 0; parity < r; ++parity)
        output[parity] = &recovery_storage[parity][1];

    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, bytes, &scratch_bytes),
        "partial T8 unaligned scratch query");
    AlignedBuffer scratch(scratch_bytes);
    leo2_encode_batch_item item = {};
    item.shard_bytes = bytes;
    item.original = &input[0];
    item.recovery = &output[0];
    item.scratch = scratch.data();
    item.scratch_bytes = scratch.size();
    leo2_encode_batch_binding* binding = NULL;
    RequireResult(leo2_encode_batch_binding_create(
        codec, &item, 1, &binding),
        "partial T8 unaligned binding create");
    RequireResult(leo2_encode_batch_binding_execute(binding),
        "partial T8 unaligned binding execute");

    uint64_t checks = 0;
    for (unsigned parity = 0; parity < r; ++parity)
    {
        const Bytes actual(
            recovery_storage[parity].begin() + 1,
            recovery_storage[parity].end() - 1);
        Require(actual == OracleParity(
            field, generator, original, parity),
            "partial T8 unaligned parity differs from oracle");
        Require(recovery_storage[parity].front() == sentinel &&
                recovery_storage[parity].back() == sentinel,
            "partial T8 unaligned encode modified a guard byte");
        checks += 2;
    }

    leo2_encode_batch_binding_destroy(binding);
    leo2_codec_destroy(codec);
    Require(checks == 2 * r,
        "partial T8 unaligned check count changed unexpectedly");
    return checks;
}

uint64_t ExerciseT8TwoBlockUnaligned(
    leo2_context* context,
    size_t bytes)
{
    static const unsigned k = 13;
    static const unsigned r = 5;
    static const uint8_t sentinel = 0xa5;

    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        "two-block T8 unaligned codec");
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh, k, r);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);
    const Shards original = MakeOriginal(k, bytes);
    Shards original_storage(k, Bytes(bytes + 2, sentinel));
    std::vector<const void*> input(k, NULL);
    for (unsigned source = 0; source < k; ++source)
    {
        std::copy(original[source].begin(), original[source].end(),
            original_storage[source].begin() + 1);
        input[source] = &original_storage[source][1];
    }
    Shards recovery_storage(r, Bytes(bytes + 2, sentinel));
    std::vector<void*> output(r, NULL);
    for (unsigned parity = 0; parity < r; ++parity)
        output[parity] = &recovery_storage[parity][1];

    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, bytes, &scratch_bytes),
        "two-block T8 unaligned scratch query");
    AlignedBuffer scratch(scratch_bytes);
    leo2_encode_batch_item item = {};
    item.shard_bytes = bytes;
    item.original = &input[0];
    item.recovery = &output[0];
    item.scratch = scratch.data();
    item.scratch_bytes = scratch.size();
    leo2_encode_batch_binding* binding = NULL;
    RequireResult(leo2_encode_batch_binding_create(
        codec, &item, 1, &binding),
        "two-block T8 unaligned binding create");
    RequireResult(leo2_encode_batch_binding_execute(binding),
        "two-block T8 unaligned binding execute");

    uint64_t checks = 0;
    for (unsigned parity = 0; parity < r; ++parity)
    {
        const Bytes actual(
            recovery_storage[parity].begin() + 1,
            recovery_storage[parity].end() - 1);
        Require(actual == OracleParity(
            field, generator, original, parity),
            "two-block T8 unaligned parity differs from oracle");
        Require(recovery_storage[parity].front() == sentinel &&
                recovery_storage[parity].back() == sentinel,
            "two-block T8 unaligned encode modified a guard byte");
        checks += 2;
    }

    leo2_encode_batch_binding_destroy(binding);
    leo2_codec_destroy(codec);
    Require(checks == 10,
        "two-block T8 unaligned check count changed unexpectedly");
    return checks;
}

void ExerciseTinyFullOutputRegion(
    leo2_context* context,
    uint64_t& codec_checks,
    uint64_t& encode_checks,
    uint64_t& direct_checks,
    uint64_t& transform_checks)
{
    const size_t byte_counts[] = {
        1, 2, 3, 4, 7, 8, 15, 16,
        17, 31, 32, 33, 63, 64, 65, 70,
        128, 192, 256, 320, 384, 448, 512, 513,
        576, 640, 704, 768, 832, 896, 960, 1024, 1025, 1088
    };
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();

    for (unsigned k = 5; k <= 16; ++k)
    {
        const unsigned maximum_r = k < 8 ? k : 8;
        for (unsigned r = 5; r <= maximum_r; ++r)
        {
            leo2_codec* codec = NULL;
            RequireResult(leo2_codec_create(context, k, r,
                LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
                "tiny full-output codec");
            const leopard2_test::ProfileLayout layout =
                leopard2_test::make_profile_layout(
                    leopard2_test::kLegacyHigh, k, r);
            const leopard2_test::Matrix generator =
                leopard2_test::direct_systematic_generator(field, layout);

            for (size_t byte_index = 0;
                 byte_index < sizeof(byte_counts) / sizeof(byte_counts[0]);
                 ++byte_index)
            {
                const size_t bytes = byte_counts[byte_index];
                leopard2_internal::CodecEncodePathInfo path = {};
                Require(leopard2_internal::GetCodecEncodePathInfo(
                    codec, bytes, r, &path),
                    "tiny full-output path introspection");
                const bool expected_partial_binding =
                    LEO2_EXPECT_HIGH_T8_PARTIAL_BINDING != 0 &&
                    IsExpectedT8OneBlockByteCount(k, r, bytes) && k <= 8 &&
                    (k != 8 || r != 8);
                Require(path.high_t8_partial_binding_selected ==
                        expected_partial_binding,
                    "tiny partial T8 selector differs from expectation");
                const bool expected_two_block_binding =
                    LEO2_EXPECT_HIGH_T8_TWO_BLOCK_BINDING != 0 &&
                    IsExpectedT8TwoBlockByteCount(k, r, bytes) &&
                    k >= 9;
                Require(path.high_t8_two_block_binding_selected ==
                        expected_two_block_binding,
                    "tiny two-block T8 selector differs from expectation");
#if LEO2_EXPECT_HIGH_DIRECT_PRODUCTION
                Require(path.direct_generator_rows == r,
                    "tiny production codec omitted direct generator rows");
#else
                Require(path.direct_generator_rows == 0,
                    "tiny default production codec retained direct rows");
#endif
                if (path.auto_direct_selected)
                    ++direct_checks;
                else
                    ++transform_checks;

                const Shards original = MakeOriginal(k, bytes);
                const Shards original_before = original;
                EncodeAllAndCheck(
                    codec, original, field, generator, r);
                Require(original == original_before,
                    "tiny full-output AUTO modified an original shard");
                ++encode_checks;
            }
            leo2_codec_destroy(codec);
            ++codec_checks;
        }
    }

    Require(codec_checks == 42,
        "tiny full-output shape count changed unexpectedly");
    Require(encode_checks == codec_checks *
        (sizeof(byte_counts) / sizeof(byte_counts[0])),
        "tiny full-output execution count changed unexpectedly");
#if LEO2_EXPECT_HIGH_DIRECT_AUTO
    Require(direct_checks != 0 && transform_checks != 0,
        "AUTO-on tiny coverage did not exercise both encoder paths");
#else
    Require(direct_checks == 0,
        "AUTO-off tiny coverage unexpectedly selected direct encoding");
#endif
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
#else
        Require(path.direct_generator_rows == 0,
            "default-OFF production codec retained direct generator rows");
#endif
#if LEO2_EXPECT_HIGH_DIRECT_AUTO
        Require(path.auto_direct_selected,
            "AUTO-on production codec did not select the direct encoder");
#else
        Require(!path.auto_direct_selected,
            "AUTO-off production codec escaped the transform encoder");
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
        uint64_t tiny_codec_checks = 0;
        uint64_t tiny_encode_checks = 0;
        uint64_t tiny_direct_checks = 0;
        uint64_t tiny_transform_checks = 0;
        ExerciseTinyFullOutputRegion(
            context, tiny_codec_checks, tiny_encode_checks,
            tiny_direct_checks, tiny_transform_checks);
        const uint64_t t4_binding_checks =
            ExerciseT4BatchBindings(context);
        const uint64_t t8_binding_checks =
            ExerciseT8BatchBinding(context, 64) +
            ExerciseT8BatchBinding(context, 128) +
            ExerciseT8BatchBinding(context, 192) +
            ExerciseT8BatchBinding(context, 256) +
            ExerciseT8BatchBinding(context, 320) +
            ExerciseT8BatchBinding(context, 384) +
            ExerciseT8BatchBinding(context, 448) +
            ExerciseT8BatchBinding(context, 512) +
            ExerciseT8BatchBinding(context, 576) +
            ExerciseT8BatchBinding(context, 640) +
            ExerciseT8BatchBinding(context, 704) +
            ExerciseT8BatchBinding(context, 768) +
            ExerciseT8BatchBinding(context, 832) +
            ExerciseT8BatchBinding(context, 896) +
            ExerciseT8BatchBinding(context, 960) +
            ExerciseT8BatchBinding(context, 1024);
        const uint64_t t8_partial_binding_checks =
            ExerciseT8PartialBindings(context);
        const uint64_t t8_two_block_binding_checks =
            ExerciseT8TwoBlockBindings(context);
        const uint64_t t8_partial_unaligned_checks =
            ExerciseT8PartialUnaligned(context, 64) +
            ExerciseT8PartialUnaligned(context, 128) +
            ExerciseT8PartialUnaligned(context, 192) +
            ExerciseT8PartialUnaligned(context, 256) +
            ExerciseT8PartialUnaligned(context, 320) +
            ExerciseT8PartialUnaligned(context, 384) +
            ExerciseT8PartialUnaligned(context, 448) +
            ExerciseT8PartialUnaligned(context, 512) +
            ExerciseT8PartialUnaligned(context, 576) +
            ExerciseT8PartialUnaligned(context, 640) +
            ExerciseT8PartialUnaligned(context, 704) +
            ExerciseT8PartialUnaligned(context, 768) +
            ExerciseT8PartialUnaligned(context, 832) +
            ExerciseT8PartialUnaligned(context, 896) +
            ExerciseT8PartialUnaligned(context, 960) +
            ExerciseT8PartialUnaligned(context, 1024) +
            ExerciseT8PartialUnaligned(context, 1088);
        const uint64_t t8_two_block_unaligned_checks =
            ExerciseT8TwoBlockUnaligned(context, 64) +
            ExerciseT8TwoBlockUnaligned(context, 128) +
            ExerciseT8TwoBlockUnaligned(context, 192) +
            ExerciseT8TwoBlockUnaligned(context, 256) +
            ExerciseT8TwoBlockUnaligned(context, 320) +
            ExerciseT8TwoBlockUnaligned(context, 384) +
            ExerciseT8TwoBlockUnaligned(context, 448) +
            ExerciseT8TwoBlockUnaligned(context, 512) +
            ExerciseT8TwoBlockUnaligned(context, 576) +
            ExerciseT8TwoBlockUnaligned(context, 640) +
            ExerciseT8TwoBlockUnaligned(context, 704) +
            ExerciseT8TwoBlockUnaligned(context, 768) +
            ExerciseT8TwoBlockUnaligned(context, 832) +
            ExerciseT8TwoBlockUnaligned(context, 896) +
            ExerciseT8TwoBlockUnaligned(context, 960) +
            ExerciseT8TwoBlockUnaligned(context, 1024);
        leo2_context_destroy(context);
        const uint64_t t8_partial_thread_pool_checks =
            ExerciseT8PartialThreadPool(64) +
            ExerciseT8PartialThreadPool(128) +
            ExerciseT8PartialThreadPool(192) +
            ExerciseT8PartialThreadPool(256) +
            ExerciseT8PartialThreadPool(320) +
            ExerciseT8PartialThreadPool(384) +
            ExerciseT8PartialThreadPool(448) +
            ExerciseT8PartialThreadPool(512) +
            ExerciseT8PartialThreadPool(576) +
            ExerciseT8PartialThreadPool(640) +
            ExerciseT8PartialThreadPool(704) +
            ExerciseT8PartialThreadPool(768) +
            ExerciseT8PartialThreadPool(832) +
            ExerciseT8PartialThreadPool(896) +
            ExerciseT8PartialThreadPool(960) +
            ExerciseT8PartialThreadPool(1024) +
            ExerciseT8PartialThreadPool(1088);
        const uint64_t t8_two_block_thread_pool_checks =
            ExerciseT8TwoBlockThreadPool(128) +
            ExerciseT8TwoBlockThreadPool(192) +
            ExerciseT8TwoBlockThreadPool(256) +
            ExerciseT8TwoBlockThreadPool(320) +
            ExerciseT8TwoBlockThreadPool(384) +
            ExerciseT8TwoBlockThreadPool(448) +
            ExerciseT8TwoBlockThreadPool(512) +
            ExerciseT8TwoBlockThreadPool(576) +
            ExerciseT8TwoBlockThreadPool(640) +
            ExerciseT8TwoBlockThreadPool(704) +
            ExerciseT8TwoBlockThreadPool(768) +
            ExerciseT8TwoBlockThreadPool(832) +
            ExerciseT8TwoBlockThreadPool(896) +
            ExerciseT8TwoBlockThreadPool(960) +
            ExerciseT8TwoBlockThreadPool(1024);
#if LEO2_EXPECT_HIGH_DIRECT_PRODUCTION
        const char* const table_state = "ON";
#else
        const char* const table_state = "OFF";
#endif
#if LEO2_EXPECT_HIGH_DIRECT_AUTO
        const char* const auto_state = "ON";
#else
        const char* const auto_state = "OFF";
#endif
#if LEO2_EXPECT_HIGH_T8_VECTOR
        const char* const t8_state = "ON";
#else
        const char* const t8_state = "OFF";
#endif
#if LEO2_EXPECT_HIGH_T8_PARTIAL_BINDING
        const char* const t8_partial_state = "ON";
#else
        const char* const t8_partial_state = "OFF";
#endif
#if LEO2_EXPECT_HIGH_T8_TWO_BLOCK_BINDING
        const char* const t8_two_block_state = "ON";
#else
        const char* const t8_two_block_state = "OFF";
#endif
        std::printf(
            "Leopard2 production high-direct smoke passed: "
            "tables=%s auto=%s t8_vector=%s t8_partial_binding=%s "
            "t8_two_block_binding=%s "
            "K=2 R=16 bytes=4096 Q=1 parity=0,15 "
            "tiny_codecs=%llu tiny_encodes=%llu direct=%llu transform=%llu "
            "t4_binding_checks=%llu t8_binding_checks=%llu "
            "t8_partial_binding_checks=%llu "
            "t8_two_block_binding_checks=%llu "
            "t8_partial_unaligned_checks=%llu "
            "t8_two_block_unaligned_checks=%llu "
            "t8_partial_thread_pool_checks=%llu "
            "t8_two_block_thread_pool_checks=%llu\n",
            table_state, auto_state, t8_state, t8_partial_state,
            t8_two_block_state,
            static_cast<unsigned long long>(tiny_codec_checks),
            static_cast<unsigned long long>(tiny_encode_checks),
            static_cast<unsigned long long>(tiny_direct_checks),
            static_cast<unsigned long long>(tiny_transform_checks),
            static_cast<unsigned long long>(t4_binding_checks),
            static_cast<unsigned long long>(t8_binding_checks),
            static_cast<unsigned long long>(t8_partial_binding_checks),
            static_cast<unsigned long long>(
                t8_two_block_binding_checks),
            static_cast<unsigned long long>(
                t8_partial_unaligned_checks),
            static_cast<unsigned long long>(
                t8_two_block_unaligned_checks),
            static_cast<unsigned long long>(
                t8_partial_thread_pool_checks),
            static_cast<unsigned long long>(
                t8_two_block_thread_pool_checks));
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
