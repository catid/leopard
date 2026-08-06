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

#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

namespace {

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

void FillInput(uint8_t* data, size_t bytes, uint64_t seed)
{
    uint64_t state = seed;
    for (size_t i = 0; i < bytes; ++i)
    {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        data[i] = static_cast<uint8_t>(state >> 24);
    }
}

void CheckParity(
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator,
    const std::vector<const void*>& original,
    const std::vector<void*>& recovery,
    size_t shard_bytes)
{
    const unsigned count = static_cast<unsigned>(original.size());
    for (unsigned parity = 0; parity < count; ++parity)
    {
        const std::vector<leopard2_test::Element>& row =
            generator[count + parity];
        const uint8_t* const output =
            static_cast<const uint8_t*>(recovery[parity]);
        for (size_t offset = 0; offset < shard_bytes; ++offset)
        {
            leopard2_test::Element expected = 0;
            for (unsigned source = 0; source < count; ++source)
            {
                expected = field.add(expected, field.multiply(row[source],
                    static_cast<const uint8_t*>(
                        original[source])[offset]));
            }
            Require(output[offset] == expected,
                "production prepared T16 parity differs from direct oracle");
        }
    }
}

void ExerciseCell(
    leo2_context* context,
    unsigned count,
    size_t shard_bytes)
{
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, count, count,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        LEO2_SUCCESS, "create production prepared T16 codec");
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, shard_bytes, &scratch_bytes), LEO2_SUCCESS,
        "query production prepared T16 scratch");
    const size_t input_bytes = static_cast<size_t>(count) * shard_bytes;
    const size_t output_bytes = static_cast<size_t>(count) * shard_bytes;
    AlignedBuffer input(input_bytes + 1U);
    AlignedBuffer output(output_bytes + 3U);
    AlignedBuffer scratch(scratch_bytes);
    FillInput(input.bytes() + 1U, input_bytes,
        UINT64_C(0x54313650524f4400) ^ count ^ shard_bytes);
    std::vector<const void*> original(count);
    std::vector<void*> recovery(count);
    for (unsigned i = 0; i < count; ++i)
    {
        original[i] = input.bytes() + 1U +
            static_cast<size_t>(i) * shard_bytes;
        recovery[i] = output.bytes() + 3U +
            static_cast<size_t>(i) * shard_bytes;
    }

    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh, count, count);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);

    RequireResult(leo2_encode(codec, shard_bytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute production prepared T16 encode");
    CheckParity(field, generator, original, recovery, shard_bytes);

    leo2_encode_batch_item item = {
        shard_bytes, &original[0], &recovery[0], scratch.data(), scratch.size()
    };
    leo2_encode_batch_binding* binding = NULL;
    RequireResult(leo2_encode_batch_binding_create(
        codec, &item, 1, &binding), LEO2_SUCCESS,
        "create production prepared T16 binding");
    std::memset(output.bytes() + 3U, 0xa5, output_bytes);
    RequireResult(leo2_encode_batch_binding_execute(binding), LEO2_SUCCESS,
        "execute production prepared T16 binding");
    CheckParity(field, generator, original, recovery, shard_bytes);
    leo2_encode_batch_binding_destroy(binding);
    leo2_codec_destroy(codec);
}

void ExerciseConcurrent(leo2_context* context)
{
    static const unsigned kCount = 12;
    static const size_t kShardBytes = 1024;
    static const unsigned kThreadCount = 8;
    static const unsigned kRepetitions = 8;
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, kCount, kCount,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        LEO2_SUCCESS, "create concurrent prepared T16 codec");
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, kShardBytes, &scratch_bytes), LEO2_SUCCESS,
        "query concurrent prepared T16 scratch");
    const size_t slab_bytes = kCount * kShardBytes;
    AlignedBuffer input(slab_bytes);
    AlignedBuffer expected(slab_bytes);
    AlignedBuffer reference_scratch(scratch_bytes);
    FillInput(input.bytes(), slab_bytes, UINT64_C(0x5431365448524541));
    std::vector<const void*> original(kCount);
    std::vector<void*> expected_recovery(kCount);
    for (unsigned i = 0; i < kCount; ++i)
    {
        original[i] = input.bytes() + i * kShardBytes;
        expected_recovery[i] = expected.bytes() + i * kShardBytes;
    }
    RequireResult(leo2_encode(codec, kShardBytes,
        &original[0], &expected_recovery[0], reference_scratch.data(),
        reference_scratch.size()), LEO2_SUCCESS,
        "create concurrent prepared T16 reference");

    std::atomic<unsigned> ready(0);
    std::atomic<bool> go(false);
    std::atomic<bool> failed(false);
    std::vector<std::thread> threads;
    for (unsigned thread = 0; thread < kThreadCount; ++thread)
    {
        threads.push_back(std::thread([&]() {
            ready.fetch_add(1, std::memory_order_release);
            try
            {
                AlignedBuffer output(slab_bytes);
                AlignedBuffer scratch(scratch_bytes);
                std::vector<void*> recovery(kCount);
                for (unsigned i = 0; i < kCount; ++i)
                    recovery[i] = output.bytes() + i * kShardBytes;
                while (!go.load(std::memory_order_acquire))
                    std::this_thread::yield();
                for (unsigned repetition = 0;
                     repetition < kRepetitions; ++repetition)
                {
                    if (leo2_encode(codec, kShardBytes,
                            &original[0], &recovery[0], scratch.data(),
                            scratch.size()) != LEO2_SUCCESS ||
                        std::memcmp(output.bytes(), expected.bytes(),
                            slab_bytes) != 0)
                    {
                        failed.store(true, std::memory_order_relaxed);
                        return;
                    }
                }
            }
            catch (...)
            {
                failed.store(true, std::memory_order_relaxed);
            }
        }));
    }
    while (ready.load(std::memory_order_acquire) != kThreadCount)
        std::this_thread::yield();
    go.store(true, std::memory_order_release);
    for (size_t i = 0; i < threads.size(); ++i)
        threads[i].join();
    Require(!failed.load(std::memory_order_relaxed),
        "concurrent prepared T16 encode differed from reference");
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
                "production prepared T16 test skipped: AVX2 unavailable\n");
            return 0;
        }
        RequireResult(result, LEO2_SUCCESS,
            "create production prepared T16 context");
        static const size_t shard_bytes[] = {
            64, 128, 192, 256, 320, 384, 448, 512, 1024, 2048
        };
        for (unsigned count = 9; count <= 16; ++count)
        {
            for (size_t i = 0;
                 i < sizeof(shard_bytes) / sizeof(shard_bytes[0]); ++i)
                ExerciseCell(context, count, shard_bytes[i]);
        }
        ExerciseConcurrent(context);
        leo2_context_destroy(context);
        std::printf("production prepared T16 checks passed\n");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "production prepared T16 failure: %s\n",
            error.what());
        return 1;
    }
}
