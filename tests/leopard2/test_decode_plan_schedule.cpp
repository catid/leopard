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

#include "Leopard2Plan.h"
#include "LeopardFF8.h"
#include "LeopardFF16.h"
#include "leopard.h"
#include "leopard2.h"

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <new>
#include <thread>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

static std::atomic<bool> gTrackAllocations(false);
static std::atomic<uint64_t> gTrackedAllocations(0);

#if defined(_MSC_VER)
#define LEO2_TEST_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define LEO2_TEST_NOINLINE __attribute__((noinline))
#else
#define LEO2_TEST_NOINLINE
#endif

// Keep GCC from inlining across this legal replacement new/delete pair and
// misdiagnosing its intentional malloc/free implementation as a caller error.
LEO2_TEST_NOINLINE void* operator new(size_t bytes)
{
    if (gTrackAllocations.load(std::memory_order_relaxed))
        gTrackedAllocations.fetch_add(1, std::memory_order_relaxed);
    void* result = malloc(bytes == 0 ? 1 : bytes);
    if (!result)
        throw std::bad_alloc();
    return result;
}

LEO2_TEST_NOINLINE void* operator new[](size_t bytes)
{
    return ::operator new(bytes);
}

LEO2_TEST_NOINLINE void* operator new(
    size_t bytes, const std::nothrow_t&) noexcept
{
    try
    {
        return ::operator new(bytes);
    }
    catch (...)
    {
        return NULL;
    }
}

LEO2_TEST_NOINLINE void* operator new[](
    size_t bytes, const std::nothrow_t&) noexcept
{
    return ::operator new(bytes, std::nothrow);
}

LEO2_TEST_NOINLINE void operator delete(void* pointer) noexcept
{
    free(pointer);
}

LEO2_TEST_NOINLINE void operator delete[](void* pointer) noexcept
{
    free(pointer);
}

LEO2_TEST_NOINLINE void operator delete(
    void* pointer, const std::nothrow_t&) noexcept
{
    ::operator delete(pointer);
}

LEO2_TEST_NOINLINE void operator delete[](
    void* pointer, const std::nothrow_t&) noexcept
{
    ::operator delete[](pointer);
}

#undef LEO2_TEST_NOINLINE

namespace {

static const size_t kBytes = 64;

static void Require(bool condition, const char* message)
{
    if (!condition)
    {
        std::cerr << "decode-plan schedule test failed: " << message << std::endl;
        std::exit(1);
    }
}

static uint32_t Mix(uint32_t value)
{
    value ^= value >> 16;
    value *= 0x7feb352du;
    value ^= value >> 15;
    value *= 0x846ca68bu;
    return value ^ (value >> 16);
}

class AlignedBytes
{
public:
    explicit AlignedBytes(size_t bytes) : data_(NULL), bytes_(bytes)
    {
        if (bytes_ == 0)
            return;
#if defined(_MSC_VER)
        data_ = _aligned_malloc(bytes_, 64);
#else
        if (posix_memalign(&data_, 64, bytes) != 0)
            data_ = NULL;
#endif
        if (data_)
            memset(data_, 0, bytes);
    }
    ~AlignedBytes()
    {
#if defined(_MSC_VER)
        _aligned_free(data_);
#else
        free(data_);
#endif
    }
    uint8_t* data() { return static_cast<uint8_t*>(data_); }
    const uint8_t* data() const { return static_cast<const uint8_t*>(data_); }
    bool valid() const { return data_ != NULL; }
    size_t size() const { return bytes_; }

private:
    AlignedBytes(const AlignedBytes&);
    AlignedBytes& operator=(const AlignedBytes&);
    void* data_;
    size_t bytes_;
};

static bool AnyRequested(
    const std::vector<uint8_t>& requested,
    unsigned begin,
    unsigned end)
{
    for (unsigned i = begin; i < end; ++i)
        if (requested[i])
            return true;
    return false;
}

static uint64_t VerifyDependencyQueries(
    unsigned n,
    const std::vector<uint8_t>& requested,
    const leopard2_internal::OutputDependencyView& view)
{
    uint64_t comparisons = 0;
    unsigned mip_level = leopard2_internal::Log2PowerOfTwo(n);
    while (mip_level >= 2)
    {
        const unsigned group_size = 1u << mip_level;
        for (unsigned begin = 0; begin < n; begin += group_size)
        {
            Require(view.IsNeeded(mip_level, begin) ==
                AnyRequested(requested, begin, begin + group_size),
                "fused-group dependency mismatch");
            ++comparisons;
        }
        mip_level -= 2;
    }
    if (mip_level == 1)
    {
        for (unsigned begin = 0; begin < n; begin += 2)
        {
            Require(view.IsNeeded(1, begin) ==
                AnyRequested(requested, begin, begin + 2),
                "final-pair dependency mismatch");
            ++comparisons;
        }
    }
    return comparisons;
}

static uint64_t CheckDependencyBitmap()
{
    Require(leopard2_internal::OutputDependencyBitCount(256) == 85,
        "GF8 dependency bit count");
    Require(leopard2_internal::OutputDependencyWordCount(256) * 8 == 16,
        "GF8 dependency byte count");
    Require(leopard2_internal::OutputDependencyBitCount(65536) == 21845,
        "GF16 dependency bit count");
    Require(leopard2_internal::OutputDependencyWordCount(65536) * 8 == 2736,
        "GF16 dependency byte count");
    Require(leopard2_internal::OutputDependencyWordCount(3) == 0,
        "invalid non-power-of-two size accepted");

    uint64_t comparisons = 0;
    for (unsigned n = 2; n <= 16; n <<= 1)
    {
        const uint32_t mask_count = 1u << n;
        for (uint32_t mask = 0; mask < mask_count; ++mask)
        {
            std::vector<uint8_t> requested(n, 0);
            std::vector<uint32_t> coordinates;
            for (unsigned i = 0; i < n; ++i)
            {
                requested[i] = static_cast<uint8_t>((mask >> i) & 1u);
                if (requested[i])
                    coordinates.push_back(i);
            }
            std::vector<uint64_t> words(
                leopard2_internal::OutputDependencyWordCount(n));
            Require(leopard2_internal::BuildOutputDependencies(
                n, coordinates.empty() ? NULL : &coordinates[0], coordinates.size(),
                &words[0], words.size()), "exhaustive dependency build");
            const leopard2_internal::OutputDependencyView view =
                leopard2_internal::MakeOutputDependencyView(n, &words[0], words.size());
            comparisons += VerifyDependencyQueries(n, requested, view);
            Require(!view.IsNeeded(0, 0), "zero mip accepted");
            Require(!view.IsNeeded(view.log2_size + 1, 0), "oversized mip accepted");
            Require(!view.IsNeeded(view.log2_size, n),
                "out-of-range query coordinate accepted");
            if (view.log2_size >= 2)
                Require(!view.IsNeeded(view.log2_size - 1, 0),
                    "wrong-parity mip accepted");
        }
    }

    const unsigned sizes[] = { 32, 128, 256, 1024, 4096, 65536 };
    for (unsigned s = 0; s < sizeof(sizes) / sizeof(sizes[0]); ++s)
    {
        const unsigned n = sizes[s];
        for (unsigned pattern = 0; pattern < 24; ++pattern)
        {
            std::vector<uint8_t> requested(n, 0);
            std::vector<uint32_t> coordinates;
            const unsigned count = pattern == 0 ? 0 :
                std::min(n, 1u + pattern * pattern);
            for (unsigned i = 0; i < count; ++i)
            {
                const uint32_t coordinate = Mix(
                    0x51c9u + n * 13u + pattern * 991u + i) & (n - 1);
                if (!requested[coordinate])
                {
                    requested[coordinate] = 1;
                    coordinates.push_back(coordinate);
                }
            }
            std::sort(coordinates.begin(), coordinates.end());
            std::vector<uint64_t> words(
                leopard2_internal::OutputDependencyWordCount(n));
            Require(leopard2_internal::BuildOutputDependencies(
                n, coordinates.empty() ? NULL : &coordinates[0], coordinates.size(),
                &words[0], words.size()), "random dependency build");
            const leopard2_internal::OutputDependencyView view =
                leopard2_internal::MakeOutputDependencyView(n, &words[0], words.size());
            comparisons += VerifyDependencyQueries(n, requested, view);
        }
    }

    uint64_t invalid_word = 0;
    const uint32_t out_of_range = 8;
    Require(!leopard2_internal::BuildOutputDependencies(
        8, &out_of_range, 1, &invalid_word, 1),
        "out-of-range coordinate accepted");
    Require(!leopard2_internal::BuildOutputDependencies(
        8, NULL, 1, &invalid_word, 1), "null coordinate list accepted");
    Require(!leopard2_internal::BuildOutputDependencies(
        8, NULL, 0, &invalid_word, 2), "wrong word count accepted");
    return comparisons;
}

static void RequireResult(leo2_result actual, const char* operation)
{
    if (actual != LEO2_SUCCESS)
    {
        std::cerr << operation << ": " << leo2_result_string(actual) << std::endl;
        Require(false, operation);
    }
}

struct PublicObjects
{
    PublicObjects() : context(NULL), codec(NULL), plan(NULL) {}
    ~PublicObjects()
    {
        leo2_decode_plan_destroy(plan);
        leo2_codec_destroy(codec);
        leo2_context_destroy(context);
    }
    leo2_context* context;
    leo2_codec* codec;
    leo2_decode_plan* plan;
};

struct ConcurrentInvocation
{
    explicit ConcurrentInvocation(size_t scratch_bytes)
        : outputs(8 * kBytes), scratch(scratch_bytes), restored(32, NULL)
    {
        Require(outputs.valid() && scratch.valid(), "concurrent buffer allocation");
        for (unsigned i = 0; i < 8; ++i)
            restored[i] = outputs.data() + static_cast<size_t>(i) * kBytes;
    }
    AlignedBytes outputs;
    AlignedBytes scratch;
    std::vector<void*> restored;
};

static void CheckPublicPlanReuseAndAllocation()
{
    const uint32_t k = 32;
    const uint32_t r = 224;
    PublicObjects objects;
    leo2_context_options context_options = {};
    context_options.struct_size = sizeof(context_options);
    context_options.thread_count = 1;
    RequireResult(leo2_context_create(&context_options, &objects.context),
        "context create");
    leo2_codec_options codec_options = {};
    codec_options.struct_size = sizeof(codec_options);
    codec_options.flags = LEO2_CODEC_FORCE_SPECIALIZED_DECODE;
    RequireResult(leo2_codec_create(objects.context, k, r, LEO2_PROFILE_LOW_V1,
        LEO2_FIELD_GF8, &codec_options, &objects.codec), "codec create");

    AlignedBytes originals(static_cast<size_t>(k) * kBytes);
    AlignedBytes recovery(static_cast<size_t>(r) * kBytes);
    Require(originals.valid() && recovery.valid(), "public shard allocation");
    std::vector<const void*> encode_original(k);
    std::vector<void*> encode_recovery(r);
    for (uint32_t i = 0; i < k; ++i)
    {
        uint8_t* shard = originals.data() + static_cast<size_t>(i) * kBytes;
        for (size_t b = 0; b < kBytes; ++b)
            shard[b] = static_cast<uint8_t>(Mix(i * 173u + b));
        encode_original[i] = shard;
    }
    for (uint32_t i = 0; i < r; ++i)
        encode_recovery[i] = recovery.data() + static_cast<size_t>(i) * kBytes;
    size_t encode_scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        objects.codec, kBytes, &encode_scratch_bytes), "encode scratch query");
    AlignedBytes encode_scratch(encode_scratch_bytes);
    Require(encode_scratch.valid(), "encode scratch allocation");
    RequireResult(leo2_encode(objects.codec, kBytes, &encode_original[0],
        &encode_recovery[0], encode_scratch.data(), encode_scratch.size()),
        "encode");

    std::vector<uint8_t> original_present(k, 1), recovery_present(r, 1);
    std::vector<const void*> decode_original = encode_original;
    std::vector<const void*> decode_recovery(r);
    for (uint32_t i = 0; i < r; ++i)
        decode_recovery[i] = encode_recovery[i];
    for (uint32_t i = 0; i < 8; ++i)
    {
        original_present[i] = 0;
        decode_original[i] = NULL;
    }
    RequireResult(leo2_decode_plan_create(objects.codec, &original_present[0],
        &recovery_present[0], &objects.plan), "decode plan create");
    size_t decode_scratch_bytes = 0;
    RequireResult(leo2_decode_plan_scratch_size(
        objects.plan, kBytes, &decode_scratch_bytes), "decode scratch query");

    ConcurrentInvocation warmup(decode_scratch_bytes);
    RequireResult(leo2_decode_plan_execute(objects.plan, kBytes,
        &decode_original[0], &decode_recovery[0], &warmup.restored[0],
        warmup.scratch.data(), warmup.scratch.size()), "warm decode execute");
    for (uint32_t i = 0; i < 8; ++i)
        Require(memcmp(warmup.restored[i],
            originals.data() + static_cast<size_t>(i) * kBytes, kBytes) == 0,
            "warm decode output mismatch");

    gTrackedAllocations.store(0, std::memory_order_relaxed);
    gTrackAllocations.store(true, std::memory_order_release);
    const leo2_result allocation_result = leo2_decode_plan_execute(
        objects.plan, kBytes, &decode_original[0], &decode_recovery[0],
        &warmup.restored[0], warmup.scratch.data(), warmup.scratch.size());
    gTrackAllocations.store(false, std::memory_order_release);
    RequireResult(allocation_result, "allocation-trapped decode execute");
    Require(gTrackedAllocations.load(std::memory_order_relaxed) == 0,
        "decode execution allocated C++ storage");

    const unsigned worker_count = 8;
    std::vector<ConcurrentInvocation*> invocations(worker_count, NULL);
    for (unsigned i = 0; i < worker_count; ++i)
        invocations[i] = new ConcurrentInvocation(decode_scratch_bytes);
    std::atomic<unsigned> failures(0);
    std::vector<std::thread> threads;
    for (unsigned worker = 0; worker < worker_count; ++worker)
    {
        threads.push_back(std::thread([&, worker]() {
            ConcurrentInvocation* invocation = invocations[worker];
            for (unsigned repeat = 0; repeat < 16; ++repeat)
            {
                if (leo2_decode_plan_execute(objects.plan, kBytes,
                    &decode_original[0], &decode_recovery[0],
                    &invocation->restored[0], invocation->scratch.data(),
                    invocation->scratch.size()) != LEO2_SUCCESS)
                {
                    ++failures;
                    continue;
                }
                for (uint32_t i = 0; i < 8; ++i)
                    if (memcmp(invocation->restored[i],
                        originals.data() + static_cast<size_t>(i) * kBytes,
                        kBytes) != 0)
                        ++failures;
            }
        }));
    }
    for (size_t i = 0; i < threads.size(); ++i)
        threads[i].join();
    for (unsigned i = 0; i < worker_count; ++i)
        delete invocations[i];
    Require(failures.load() == 0, "concurrent immutable plan reuse mismatch");
}

template<class Ffe>
static std::vector<Ffe> MakeLocator(unsigned n, unsigned modulus, uint32_t seed)
{
    std::vector<Ffe> result(n);
    for (unsigned i = 0; i < n; ++i)
        result[i] = static_cast<Ffe>(Mix(seed + i * 17u) % modulus);
    return result;
}

static void PrepareWorkPointers(
    unsigned n,
    AlignedBytes& storage,
    std::vector<void*>& work)
{
    Require(storage.valid(), "aligned work allocation");
    work.resize(n);
    for (unsigned i = 0; i < n; ++i)
        work[i] = storage.data() + static_cast<size_t>(i) * kBytes;
}

static std::vector<uint32_t> CoordinatesFromMask(
    const std::vector<uint8_t>& requested)
{
    std::vector<uint32_t> result;
    for (unsigned i = 0; i < requested.size(); ++i)
        if (requested[i])
            result.push_back(i);
    return result;
}

static std::vector<uint16_t> BlockInputCounts(
    const std::vector<const void*>& coordinate_data,
    unsigned block_size)
{
    std::vector<uint16_t> counts(coordinate_data.size() / block_size, 0);
    for (unsigned i = 0; i < coordinate_data.size(); ++i)
        if (coordinate_data[i])
            counts[i / block_size] = static_cast<uint16_t>(i % block_size + 1);
    return counts;
}

static std::vector<leopard2_internal::DecodeOutputBlock> OutputBlocks(
    const std::vector<uint32_t>& coordinates,
    unsigned block_size)
{
    std::vector<leopard2_internal::DecodeOutputBlock> blocks;
    size_t begin = 0;
    while (begin < coordinates.size())
    {
        const uint32_t block = coordinates[begin] / block_size;
        size_t end = begin + 1;
        while (end < coordinates.size() && coordinates[end] / block_size == block)
            ++end;
        leopard2_internal::DecodeOutputBlock descriptor = {
            block,
            coordinates[end - 1] % block_size + 1,
            static_cast<uint32_t>(begin),
            static_cast<uint32_t>(end)
        };
        blocks.push_back(descriptor);
        begin = end;
    }
    return blocks;
}

template<class Ffe>
static uint64_t CheckFieldKernels(
    unsigned n,
    unsigned side,
    unsigned modulus,
    void (*generic_prepared)(uint64_t, unsigned, const void* const*,
        const uint8_t*, const Ffe*, void**),
    void (*generic_planned)(uint64_t, unsigned, const void* const*, unsigned,
        const uint32_t*, unsigned,
        const leopard2_internal::OutputDependencyView&, const Ffe*, void**),
    void (*prepare_low)(unsigned, unsigned, Ffe*),
    void (*low_prepared)(uint64_t, unsigned, unsigned, const void* const*,
        const uint8_t*, const Ffe*, const Ffe*, void**),
    void (*low_planned)(uint64_t, unsigned, unsigned, const void* const*,
        const uint16_t*, const uint32_t*, unsigned,
        const leopard2_internal::OutputDependencyView&, const Ffe*, const Ffe*,
        void**),
    void (*prepare_high)(unsigned, unsigned, Ffe*),
    void (*high_prepared)(uint64_t, unsigned, unsigned, const void* const*,
        const uint8_t*, const Ffe*, const Ffe*, void**),
    void (*high_planned)(uint64_t, unsigned, unsigned, const void* const*,
        const uint16_t*, const uint32_t*,
        const leopard2_internal::DecodeOutputBlock*, unsigned, const Ffe*,
        const Ffe*, void**))
{
    AlignedBytes sources(static_cast<size_t>(n) * kBytes);
    Require(sources.valid(), "aligned source allocation");
    std::vector<const void*> coordinate_data(n, NULL);
    for (unsigned i = 0; i < n; ++i)
    {
        uint8_t* shard = sources.data() + static_cast<size_t>(i) * kBytes;
        for (unsigned b = 0; b < kBytes; ++b)
            shard[b] = static_cast<uint8_t>(Mix(i * 131u + b * 7u + n));
        // Leave complete blocks 1 and 5 empty, with sparse holes elsewhere.
        const unsigned block = i / side;
        if (block != 1 && block != 5 && (Mix(i + n) & 3u) != 0)
            coordinate_data[i] = shard;
    }
    unsigned generic_input_count = n;
    while (generic_input_count > 0 && !coordinate_data[generic_input_count - 1])
        --generic_input_count;
    const std::vector<uint16_t> block_inputs =
        BlockInputCounts(coordinate_data, side);
    const std::vector<Ffe> locator = MakeLocator<Ffe>(n, modulus, 0x1181u + n);

    AlignedBytes expected_storage(static_cast<size_t>(n) * kBytes);
    AlignedBytes actual_storage(static_cast<size_t>(n) * kBytes);
    std::vector<void*> expected_work, actual_work;
    PrepareWorkPointers(n, expected_storage, expected_work);
    PrepareWorkPointers(n, actual_storage, actual_work);

    std::vector<uint8_t> generic_requested(n, 0);
    generic_requested[0] = generic_requested[n / 3] =
        generic_requested[n - 1] = 1;
    const std::vector<uint32_t> generic_coordinates =
        CoordinatesFromMask(generic_requested);
    std::vector<uint64_t> generic_words(
        leopard2_internal::OutputDependencyWordCount(n));
    Require(leopard2_internal::BuildOutputDependencies(
        n, &generic_coordinates[0], generic_coordinates.size(),
        &generic_words[0], generic_words.size()), "generic schedule build");
    const leopard2_internal::OutputDependencyView generic_view =
        leopard2_internal::MakeOutputDependencyView(
            n, &generic_words[0], generic_words.size());
    generic_prepared(kBytes, n, &coordinate_data[0], &generic_requested[0],
        &locator[0], &expected_work[0]);
    generic_planned(kBytes, n, &coordinate_data[0], generic_input_count,
        &generic_coordinates[0], generic_coordinates.size(), generic_view,
        &locator[0], &actual_work[0]);
    Require(memcmp(expected_storage.data(), actual_storage.data(),
        expected_storage.size()) == 0, "generic planned/prepared mismatch");

    std::vector<uint8_t> low_requested(n, 0);
    low_requested[0] = low_requested[side / 2] = low_requested[side - 1] = 1;
    const std::vector<uint32_t> low_coordinates = CoordinatesFromMask(low_requested);
    std::vector<uint64_t> low_words(
        leopard2_internal::OutputDependencyWordCount(side));
    Require(leopard2_internal::BuildOutputDependencies(
        side, &low_coordinates[0], low_coordinates.size(),
        &low_words[0], low_words.size()), "low schedule build");
    const leopard2_internal::OutputDependencyView low_view =
        leopard2_internal::MakeOutputDependencyView(
            side, &low_words[0], low_words.size());
    std::vector<Ffe> low_factors(n / side - 1);
    prepare_low(n, side, &low_factors[0]);
    low_prepared(kBytes, n, side, &coordinate_data[0], &low_requested[0],
        &locator[0], &low_factors[0], &expected_work[0]);
    low_planned(kBytes, n, side, &coordinate_data[0], &block_inputs[0],
        &low_coordinates[0], low_coordinates.size(), low_view, &locator[0],
        &low_factors[0], &actual_work[0]);
    Require(memcmp(expected_storage.data(), actual_storage.data(),
        expected_storage.size()) == 0, "low planned/prepared mismatch");

    std::vector<uint8_t> high_requested(n, 0);
    high_requested[side + 1] = high_requested[side * 3 + side / 2] =
        high_requested[n - 2] = 1;
    const std::vector<uint32_t> high_coordinates =
        CoordinatesFromMask(high_requested);
    const std::vector<leopard2_internal::DecodeOutputBlock> high_blocks =
        OutputBlocks(high_coordinates, side);
    std::vector<Ffe> high_factors(n);
    prepare_high(n, side, &high_factors[0]);
    high_prepared(kBytes, n, side, &coordinate_data[0], &high_requested[0],
        &locator[0], &high_factors[0], &expected_work[0]);
    high_planned(kBytes, n, side, &coordinate_data[0], &block_inputs[0],
        &high_coordinates[0], &high_blocks[0], high_blocks.size(), &locator[0],
        &high_factors[0], &actual_work[0]);
    Require(memcmp(expected_storage.data(), actual_storage.data(),
        expected_storage.size()) == 0, "high planned/prepared mismatch");

    return static_cast<uint64_t>(n) * 3;
}

} // namespace

int main()
{
    // The public initializer probes the CPU before selecting and allocating
    // field tables.  Calling an internal field initializer first would freeze
    // scalar tables and then make a later SIMD probe observe missing tables.
    Require(leo_init() == Leopard_Success, "Leopard initialization");

    const uint64_t dependency_comparisons = CheckDependencyBitmap();
    uint64_t kernel_slots = 0;
    kernel_slots += CheckFieldKernels<leopard::ff8::ffe_t>(
        256, 32, leopard::ff8::kModulus,
        leopard::ff8::ReedSolomonDecodePrepared,
        leopard::ff8::ReedSolomonDecodePlanned,
        leopard::ff8::PrepareLowDecode,
        leopard::ff8::ReedSolomonDecodeLowPrepared,
        leopard::ff8::ReedSolomonDecodeLowPlanned,
        leopard::ff8::PrepareHighDecode,
        leopard::ff8::ReedSolomonDecodeHighPrepared,
        leopard::ff8::ReedSolomonDecodeHighPlanned);
    kernel_slots += CheckFieldKernels<leopard::ff16::ffe_t>(
        1024, 128, leopard::ff16::kModulus,
        leopard::ff16::ReedSolomonDecodePrepared,
        leopard::ff16::ReedSolomonDecodePlanned,
        leopard::ff16::PrepareLowDecode,
        leopard::ff16::ReedSolomonDecodeLowPrepared,
        leopard::ff16::ReedSolomonDecodeLowPlanned,
        leopard::ff16::PrepareHighDecode,
        leopard::ff16::ReedSolomonDecodeHighPrepared,
        leopard::ff16::ReedSolomonDecodeHighPlanned);
    CheckPublicPlanReuseAndAllocation();

    std::cout << "leopard2 decode-plan schedule tests passed: dependency_queries="
              << dependency_comparisons << " differential_work_slots="
              << kernel_slots << " concurrent_plan_executions=128"
              << " execution_cpp_allocations=0 gf16_max_schedule_bytes=2736"
              << std::endl;
    return 0;
}
