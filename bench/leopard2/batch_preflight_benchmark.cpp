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

#include "leopard2.h"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

namespace {

typedef std::vector<uint8_t> Bytes;

void RequireResult(leo2_result result, const char* operation)
{
    if (result == LEO2_SUCCESS)
        return;
    throw std::runtime_error(std::string(operation) + ": " +
        leo2_result_string(result));
}

class AlignedBuffer
{
public:
    explicit AlignedBuffer(size_t bytes)
        : data_(NULL)
        , bytes_(bytes)
    {
        if (bytes_ == 0)
            return;
#if defined(_MSC_VER)
        data_ = _aligned_malloc(bytes_, leo2_scratch_alignment());
#else
        if (posix_memalign(&data_, leo2_scratch_alignment(), bytes_) != 0)
            data_ = NULL;
#endif
        if (!data_)
            throw std::bad_alloc();
        memset(data_, 0, bytes_);
    }

    ~AlignedBuffer()
    {
#if defined(_MSC_VER)
        _aligned_free(data_);
#else
        free(data_);
#endif
    }

    void* data() { return data_; }
    size_t size() const { return bytes_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);
    void* data_;
    size_t bytes_;
};

struct Fixture
{
    Fixture(
        size_t item_count,
        size_t shard_bytes,
        uint32_t original_count,
        uint32_t recovery_count)
        : context(NULL)
        , codec(NULL)
        , plan(NULL)
        , shard_bytes(shard_bytes)
        , original(original_count, Bytes(shard_bytes, 0))
        , original_pointers(original_count, NULL)
        , reference_parity(recovery_count, Bytes(shard_bytes, 0))
        , reference_parity_pointers(recovery_count, NULL)
        , encode_outputs(item_count,
              std::vector<Bytes>(recovery_count, Bytes(shard_bytes, 0)))
        , encode_output_pointers(item_count,
              std::vector<void*>(recovery_count, NULL))
        , encode_scratches(item_count)
        , encode_items(item_count)
        , restored(item_count, Bytes(shard_bytes, 0))
        , decode_original_pointers(item_count,
              std::vector<const void*>(original_count, NULL))
        , decode_recovery_pointers(item_count,
              std::vector<const void*>(recovery_count, NULL))
        , decode_output_pointers(item_count,
              std::vector<void*>(original_count, NULL))
        , decode_scratches(item_count)
        , decode_items(item_count)
        , encode_preflight_bytes(0)
        , decode_preflight_bytes(0)
    {
        leo2_context_options options;
        memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_SCALAR;
        options.thread_count = 1;
        RequireResult(leo2_context_create(&options, &context),
            "context create");
        RequireResult(leo2_codec_create(context, original_count, recovery_count,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
            "codec create");

        uint32_t state = 1;
        for (size_t shard = 0; shard < original.size(); ++shard)
        {
            for (size_t i = 0; i < shard_bytes; ++i)
            {
                state = state * 1664525u + 1013904223u;
                original[shard][i] = static_cast<uint8_t>(state >> 24);
            }
            original_pointers[shard] = &original[shard][0];
        }
        for (size_t i = 0; i < reference_parity.size(); ++i)
            reference_parity_pointers[i] = &reference_parity[i][0];

        size_t encode_scratch_bytes = 0;
        RequireResult(leo2_encode_scratch_size(
            codec, shard_bytes, &encode_scratch_bytes),
            "encode scratch query");
        AlignedBuffer reference_scratch(encode_scratch_bytes);
        RequireResult(leo2_encode(codec, shard_bytes,
            &original_pointers[0], &reference_parity_pointers[0],
            reference_scratch.data(), reference_scratch.size()),
            "reference encode");

        std::vector<uint8_t> original_present(original_count, 1);
        std::vector<uint8_t> recovery_present(recovery_count, 0);
        original_present[0] = 0;
        recovery_present[0] = 1;
        RequireResult(leo2_decode_plan_create(codec, &original_present[0],
            &recovery_present[0], &plan), "decode plan create");
        size_t decode_scratch_bytes = 0;
        RequireResult(leo2_decode_plan_scratch_size(
            plan, shard_bytes, &decode_scratch_bytes),
            "decode scratch query");

        for (size_t item = 0; item < item_count; ++item)
        {
            for (size_t parity = 0; parity < recovery_count; ++parity)
                encode_output_pointers[item][parity] =
                    &encode_outputs[item][parity][0];
            encode_scratches[item].reset(
                new AlignedBuffer(encode_scratch_bytes));
            encode_items[item].shard_bytes = shard_bytes;
            encode_items[item].original = &original_pointers[0];
            encode_items[item].recovery = &encode_output_pointers[item][0];
            encode_items[item].scratch = encode_scratches[item]->data();
            encode_items[item].scratch_bytes =
                encode_scratches[item]->size();

            decode_original_pointers[item][0] = NULL;
            for (size_t original_index = 1;
                 original_index < original_count;
                 ++original_index)
            {
                decode_original_pointers[item][original_index] =
                    original_pointers[original_index];
            }
            decode_recovery_pointers[item][0] =
                &reference_parity[0][0];
            decode_output_pointers[item][0] = &restored[item][0];
            decode_scratches[item].reset(
                new AlignedBuffer(decode_scratch_bytes));
            decode_items[item].shard_bytes = shard_bytes;
            decode_items[item].original = &decode_original_pointers[item][0];
            decode_items[item].recovery = &decode_recovery_pointers[item][0];
            decode_items[item].restored_original =
                &decode_output_pointers[item][0];
            decode_items[item].scratch = decode_scratches[item]->data();
            decode_items[item].scratch_bytes =
                decode_scratches[item]->size();
        }

        RequireResult(leo2_encode_batch_preflight_scratch_size(
            codec, item_count, &encode_preflight_bytes),
            "encode preflight scratch query");
        RequireResult(leo2_decode_plan_batch_preflight_scratch_size(
            plan, item_count, &decode_preflight_bytes),
            "decode preflight scratch query");
        encode_preflight.reset(new AlignedBuffer(encode_preflight_bytes));
        decode_preflight.reset(new AlignedBuffer(decode_preflight_bytes));
    }

    ~Fixture()
    {
        leo2_decode_plan_destroy(plan);
        leo2_codec_destroy(codec);
        leo2_context_destroy(context);
    }

    leo2_result EncodeCompatibility()
    {
        return leo2_encode_batch(codec, &encode_items[0], encode_items.size());
    }

    leo2_result EncodeScalable()
    {
        return leo2_encode_batch_with_preflight_scratch(codec,
            &encode_items[0], encode_items.size(), encode_preflight->data(),
            encode_preflight->size());
    }

    leo2_result DecodeCompatibility()
    {
        return leo2_decode_plan_execute_batch(
            plan, &decode_items[0], decode_items.size());
    }

    leo2_result DecodeScalable()
    {
        return leo2_decode_plan_execute_batch_with_preflight_scratch(plan,
            &decode_items[0], decode_items.size(), decode_preflight->data(),
            decode_preflight->size());
    }

    leo2_context* context;
    leo2_codec* codec;
    leo2_decode_plan* plan;
    size_t shard_bytes;
    std::vector<Bytes> original;
    std::vector<const void*> original_pointers;
    std::vector<Bytes> reference_parity;
    std::vector<void*> reference_parity_pointers;
    std::vector<std::vector<Bytes> > encode_outputs;
    std::vector<std::vector<void*> > encode_output_pointers;
    std::vector<std::unique_ptr<AlignedBuffer> > encode_scratches;
    std::vector<leo2_encode_batch_item> encode_items;
    std::vector<Bytes> restored;
    std::vector<std::vector<const void*> > decode_original_pointers;
    std::vector<std::vector<const void*> > decode_recovery_pointers;
    std::vector<std::vector<void*> > decode_output_pointers;
    std::vector<std::unique_ptr<AlignedBuffer> > decode_scratches;
    std::vector<leo2_decode_batch_item> decode_items;
    size_t encode_preflight_bytes;
    size_t decode_preflight_bytes;
    std::unique_ptr<AlignedBuffer> encode_preflight;
    std::unique_ptr<AlignedBuffer> decode_preflight;

private:
    Fixture(const Fixture&);
    Fixture& operator=(const Fixture&);
};

typedef leo2_result (Fixture::*Operation)();

double Measure(Fixture& fixture, Operation operation, size_t iterations)
{
    const std::chrono::steady_clock::time_point begin =
        std::chrono::steady_clock::now();
    for (size_t i = 0; i < iterations; ++i)
        RequireResult((fixture.*operation)(), "timed batch call");
    const std::chrono::steady_clock::time_point end =
        std::chrono::steady_clock::now();
    const double nanoseconds =
        std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
    return nanoseconds / static_cast<double>(iterations);
}

double Median(std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
}

struct Comparison
{
    double compatibility;
    double scalable;
};

Comparison Compare(
    Fixture& fixture,
    Operation compatibility,
    Operation scalable,
    size_t iterations)
{
    for (unsigned i = 0; i < 8; ++i)
    {
        RequireResult((fixture.*compatibility)(), "compatibility warmup");
        RequireResult((fixture.*scalable)(), "scalable warmup");
    }
    std::vector<double> compatibility_samples;
    std::vector<double> scalable_samples;
    for (unsigned round = 0; round < 7; ++round)
    {
        double compatibility_time = 0;
        double scalable_time = 0;
        if ((round & 1u) == 0)
        {
            compatibility_time = Measure(
                fixture, compatibility, iterations);
            scalable_time = Measure(fixture, scalable, iterations);
        }
        else
        {
            scalable_time = Measure(fixture, scalable, iterations);
            compatibility_time = Measure(
                fixture, compatibility, iterations);
        }
        compatibility_samples.push_back(compatibility_time);
        scalable_samples.push_back(scalable_time);
    }
    Comparison comparison;
    comparison.compatibility = Median(compatibility_samples);
    comparison.scalable = Median(scalable_samples);
    return comparison;
}

void RunCell(
    size_t item_count,
    uint32_t original_count,
    uint32_t recovery_count)
{
    const bool large_code = original_count > 16;
    const size_t iterations = large_code
        ? (item_count == 1 ? 2000 :
            (item_count == 8 ? 300 : (item_count == 64 ? 20 : 1)))
        : (item_count == 1 ? 50000 :
            (item_count == 8 ? 10000 : (item_count == 64 ? 200 : 2)));
    Fixture fixture(item_count, 1, original_count, recovery_count);
    const Operation encode_scalable_operation =
        fixture.encode_preflight_bytes == 0
            ? &Fixture::EncodeCompatibility
            : &Fixture::EncodeScalable;
    const Operation decode_scalable_operation =
        fixture.decode_preflight_bytes == 0
            ? &Fixture::DecodeCompatibility
            : &Fixture::DecodeScalable;
    const Comparison encode = Compare(fixture,
        &Fixture::EncodeCompatibility, encode_scalable_operation, iterations);
    const Comparison decode = Compare(fixture,
        &Fixture::DecodeCompatibility, decode_scalable_operation, iterations);
    std::cout << std::fixed << std::setprecision(3)
              << "{\"k\":" << original_count
              << ",\"r\":" << recovery_count
              << ",\"batch\":" << item_count
              << ",\"iterations\":" << iterations
              << ",\"encode_compatibility_ns\":" << encode.compatibility
              << ",\"encode_scalable_ns\":" << encode.scalable
              << ",\"encode_speedup\":"
              << encode.compatibility / encode.scalable
              << ",\"decode_compatibility_ns\":" << decode.compatibility
              << ",\"decode_scalable_ns\":" << decode.scalable
              << ",\"decode_speedup\":"
              << decode.compatibility / decode.scalable
              << ",\"encode_preflight_bytes\":"
              << fixture.encode_preflight_bytes
              << ",\"decode_preflight_bytes\":"
              << fixture.decode_preflight_bytes << "}" << std::endl;
}

} // namespace

int main()
{
    try
    {
        const size_t batches[] = { 1, 8, 64, 1024 };
        const uint32_t codes[][2] = { { 3, 2 }, { 100, 28 } };
        for (size_t code = 0; code < sizeof(codes) / sizeof(codes[0]); ++code)
            for (size_t i = 0; i < sizeof(batches) / sizeof(batches[0]); ++i)
                RunCell(batches[i], codes[code][0], codes[code][1]);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "batch preflight benchmark failed: "
                  << error.what() << std::endl;
        return 1;
    }
}
