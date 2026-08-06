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
#include "allocation_audit_config.h"

#if defined(LEO2_TEST_ONE_SHOT_NO_LOSS_HOOK_BUILD)
#include "Leopard2Direct.h"
#endif

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <new>
#include <stdexcept>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

#ifndef LEO2_EXPERIMENT_ONE_SHOT_NO_LOSS_SHORT_CIRCUIT
#define LEO2_EXPERIMENT_ONE_SHOT_NO_LOSS_SHORT_CIRCUIT 1
#endif

#ifndef LEO2_TEST_EXPECT_GF8
#define LEO2_TEST_EXPECT_GF8 0
#endif

#ifndef LEO2_TEST_EXPECT_GF16
#define LEO2_TEST_EXPECT_GF16 0
#endif

#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
static std::atomic<bool> g_track_allocations(false);
static std::atomic<uint64_t> g_tracked_allocations(0);

#if defined(_MSC_VER)
#define LEO2_TEST_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define LEO2_TEST_NOINLINE __attribute__((noinline))
#else
#define LEO2_TEST_NOINLINE
#endif

LEO2_TEST_NOINLINE void* operator new(size_t bytes)
{
    if (g_track_allocations.load(std::memory_order_relaxed))
        g_tracked_allocations.fetch_add(1, std::memory_order_relaxed);
    void* const result = malloc(bytes == 0 ? 1 : bytes);
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
#endif

namespace {

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

struct DecodeObservation
{
    leo2_result result;
    uint64_t allocations;
};

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void require_result(
    leo2_result actual,
    leo2_result expected,
    const char* operation)
{
    if (actual != expected)
    {
        std::cerr << operation << ": expected " << expected
                  << ", got " << actual << std::endl;
        throw std::runtime_error(operation);
    }
}

void begin_allocation_audit()
{
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
    g_tracked_allocations.store(0, std::memory_order_relaxed);
    g_track_allocations.store(true, std::memory_order_release);
#endif
}

uint64_t end_allocation_audit()
{
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
    g_track_allocations.store(false, std::memory_order_release);
    return g_tracked_allocations.load(std::memory_order_relaxed);
#else
    return 0;
#endif
}

DecodeObservation observe_decode(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    const uint8_t* original_present,
    const uint8_t* recovery_present,
    const void* const* original,
    const void* const* recovery,
    void* const* restored,
    void* scratch,
    size_t scratch_bytes)
{
    begin_allocation_audit();
    const leo2_result result = leo2_decode(codec, shard_bytes,
        original_present, recovery_present, original, recovery, restored,
        scratch, scratch_bytes);
    const uint64_t allocations = end_allocation_audit();
    const DecodeObservation observation = { result, allocations };
    return observation;
}

void require_no_loss_observation(
    const DecodeObservation& observation,
    const char* operation)
{
    require_result(observation.result, LEO2_SUCCESS, operation);
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
# if LEO2_EXPERIMENT_ONE_SHOT_NO_LOSS_SHORT_CIRCUIT
    require(observation.allocations == 0,
        "enabled one-shot no-loss path allocated");
# else
    require(observation.allocations > 0,
        "disabled one-shot no-loss control avoided plan allocation");
# endif
#else
    (void)observation;
#endif
}

leo2_codec* create_codec(
    leo2_context* context,
    uint32_t k,
    uint32_t r,
    leo2_profile profile,
    leo2_field field,
    uint32_t flags = 0)
{
    leo2_codec_options options = {};
    options.struct_size = sizeof(options);
    options.flags = flags;
    leo2_codec* codec = NULL;
    require_result(leo2_codec_create(
        context, k, r, profile, field, &options, &codec),
        LEO2_SUCCESS, "codec create");
    return codec;
}

#if LEO2_TEST_EXPECT_GF8
class RawTransientFixture
{
public:
    RawTransientFixture(
        leo2_context* context,
        uint32_t original_count,
        uint32_t recovery_count,
        size_t shard_bytes,
        leo2_profile profile = LEO2_PROFILE_LEGACY_HIGH_V1,
        leo2_field field = LEO2_FIELD_GF8,
        uint32_t flags = 0)
        : k(original_count)
        , r(recovery_count)
        , bytes(shard_bytes)
        , codec(create_codec(context, k, r, profile, field, flags))
        , originals(k, std::vector<uint8_t>(bytes))
        , parity(r, std::vector<uint8_t>(bytes))
        , restored_storage(k, std::vector<uint8_t>(bytes, 0xa5))
        , original_present(k, 1)
        , recovery_present(r, 1)
        , original(k, NULL)
        , recovery(r, NULL)
        , restored(k, NULL)
        , encode_scratch(NULL)
        , decode_scratch(NULL)
        , decode_scratch_bytes(0)
    {
        for (uint32_t shard = 0; shard < k; ++shard)
        {
            for (size_t i = 0; i < bytes; ++i)
            {
                originals[shard][i] = static_cast<uint8_t>(
                    shard * 53u + i * 29u + (shard ^ i) * 7u + 11u);
            }
            original[shard] = originals[shard].data();
        }
        std::vector<void*> parity_output(r, NULL);
        for (uint32_t i = 0; i < r; ++i)
            parity_output[i] = parity[i].data();
        size_t encode_scratch_bytes = 0;
        require_result(leo2_encode_scratch_size(
            codec, bytes, &encode_scratch_bytes),
            LEO2_SUCCESS, "raw transient encode scratch query");
        encode_scratch = new AlignedBuffer(encode_scratch_bytes);
        require_result(leo2_encode(codec, bytes, original.data(),
            parity_output.data(), encode_scratch->data(),
            encode_scratch->size()), LEO2_SUCCESS,
            "raw transient fixture encode");
        require_result(leo2_decode_scratch_size(
            codec, bytes, &decode_scratch_bytes),
            LEO2_SUCCESS, "raw transient decode scratch query");
        decode_scratch = new AlignedBuffer(decode_scratch_bytes);
    }

    ~RawTransientFixture()
    {
        delete decode_scratch;
        delete encode_scratch;
        leo2_codec_destroy(codec);
    }

    void Configure(uint32_t losses, unsigned shape, bool mixed_recovery)
    {
        require(losses > 0 && losses <= k && losses <= r,
            "invalid raw transient fixture loss count");
        std::fill(original_present.begin(), original_present.end(), 1);
        std::fill(recovery_present.begin(), recovery_present.end(), 1);
        if (shape == 0)
        {
            for (uint32_t i = 0; i < losses; ++i)
                original_present[i] = 0;
        }
        else if (shape == 1)
        {
            for (uint32_t i = 0; i < losses; ++i)
                original_present[k - 1U - i] = 0;
        }
        else
        {
            // A k-1 stride is coprime to every k, so this remains a scattered
            // full permutation even at K=7/14 boundaries where the former
            // fixed stride of seven could cycle forever.
            for (uint32_t marked = 0; marked < losses; ++marked)
            {
                const uint32_t candidate = static_cast<uint32_t>(
                    (3U + static_cast<uint64_t>(marked) * (k - 1U)) % k);
                original_present[candidate] = 0;
            }
        }

        if (mixed_recovery && losses < r)
        {
            const uint32_t target = std::min<uint32_t>(r, losses + 3U);
            uint32_t present = r;
            for (uint32_t i = 1; i < r && present > target; i += 3)
            {
                recovery_present[i] = 0;
                --present;
            }
            for (uint32_t i = 0; i < r && present > target; ++i)
            {
                if (recovery_present[i])
                {
                    recovery_present[i] = 0;
                    --present;
                }
            }
        }

        for (uint32_t i = 0; i < k; ++i)
        {
            original[i] = original_present[i]
                ? static_cast<const void*>(originals[i].data()) : NULL;
            restored[i] = original_present[i]
                ? NULL : static_cast<void*>(restored_storage[i].data());
        }
        for (uint32_t i = 0; i < r; ++i)
        {
            recovery[i] = recovery_present[i]
                ? static_cast<const void*>(parity[i].data()) : NULL;
        }
        ResetOutputs();
    }

    void ResetOutputs()
    {
        for (uint32_t i = 0; i < k; ++i)
            std::fill(restored_storage[i].begin(),
                restored_storage[i].end(), 0xa5);
    }

    void SelectRecoverySubset(uint32_t retained, unsigned shape)
    {
        require(retained >= Missing().size() && retained <= r,
            "invalid raw transient recovery subset size");
        std::fill(recovery_present.begin(), recovery_present.end(), 0);
        if (shape == 0)
        {
            for (uint32_t i = 0; i < retained; ++i)
                recovery_present[i] = 1;
        }
        else if (shape == 1)
        {
            for (uint32_t i = 0; i < retained; ++i)
                recovery_present[r - 1U - i] = 1;
        }
        else
        {
            uint32_t selected = 0;
            for (uint32_t distance = 0;
                 distance < r && selected < retained; ++distance)
            {
                const uint32_t index = (distance & 1U) == 0
                    ? distance / 2U : r - 1U - distance / 2U;
                if (!recovery_present[index])
                {
                    recovery_present[index] = 1;
                    ++selected;
                }
            }
        }
        for (uint32_t i = 0; i < r; ++i)
        {
            recovery[i] = recovery_present[i]
                ? static_cast<const void*>(parity[i].data()) : NULL;
        }
        ResetOutputs();
    }

    bool OutputsMatch() const
    {
        for (uint32_t i = 0; i < k; ++i)
        {
            if (!original_present[i] &&
                restored_storage[i] != originals[i])
                return false;
        }
        return true;
    }

    bool OutputsUntouched() const
    {
        for (uint32_t i = 0; i < k; ++i)
        {
            if (!original_present[i])
            {
                for (size_t j = 0; j < bytes; ++j)
                    if (restored_storage[i][j] != 0xa5)
                        return false;
            }
        }
        return true;
    }

    std::vector<uint32_t> Missing() const
    {
        std::vector<uint32_t> result;
        for (uint32_t i = 0; i < k; ++i)
            if (!original_present[i])
                result.push_back(i);
        return result;
    }

    DecodeObservation Observe(
        void* scratch_override = NULL,
        size_t scratch_bytes_override = std::numeric_limits<size_t>::max())
    {
        void* selected_scratch = scratch_override
            ? scratch_override : decode_scratch->data();
        const size_t selected_bytes =
            scratch_bytes_override == std::numeric_limits<size_t>::max()
                ? decode_scratch_bytes : scratch_bytes_override;
        return observe_decode(codec, bytes, original_present.data(),
            recovery_present.data(), original.data(), recovery.data(),
            restored.data(), selected_scratch, selected_bytes);
    }

    uint32_t k;
    uint32_t r;
    size_t bytes;
    leo2_codec* codec;
    std::vector<std::vector<uint8_t> > originals;
    std::vector<std::vector<uint8_t> > parity;
    std::vector<std::vector<uint8_t> > restored_storage;
    std::vector<uint8_t> original_present;
    std::vector<uint8_t> recovery_present;
    std::vector<const void*> original;
    std::vector<const void*> recovery;
    std::vector<void*> restored;
    AlignedBuffer* encode_scratch;
    AlignedBuffer* decode_scratch;
    size_t decode_scratch_bytes;

private:
    RawTransientFixture(const RawTransientFixture&);
    RawTransientFixture& operator=(const RawTransientFixture&);
};

void require_raw_transient_success(
    const DecodeObservation& observation,
    const RawTransientFixture& fixture,
    bool expect_allocation_free,
    const char* operation)
{
    require_result(observation.result, LEO2_SUCCESS, operation);
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
    const std::string allocation_message = expect_allocation_free
        ? std::string(operation) +
            " eligible raw transient decode allocated"
        : std::string(operation) +
            " fallback unexpectedly avoided plan allocation";
    if (expect_allocation_free)
        require(observation.allocations == 0,
            allocation_message.c_str());
    else
        require(observation.allocations > 0,
            allocation_message.c_str());
#else
    (void)expect_allocation_free;
#endif
    const std::string output_message =
        std::string(operation) + " restored incorrect bytes";
    require(fixture.OutputsMatch(), output_message.c_str());
}

void run_raw_transient_case(
    leo2_context* context,
    uint32_t k,
    uint32_t r,
    size_t bytes,
    uint32_t losses,
    unsigned shape,
    bool mixed_recovery,
    bool expect_allocation_free)
{
    RawTransientFixture fixture(context, k, r, bytes);
    fixture.Configure(losses, shape, mixed_recovery);
    const std::string label =
        "raw transient one-shot decode K=" + std::to_string(k) +
        " R=" + std::to_string(r) + " B=" + std::to_string(bytes) +
        " L=" + std::to_string(losses) +
        " shape=" + std::to_string(shape) +
        " mixed=" + std::to_string(mixed_recovery ? 1 : 0);
    require_raw_transient_success(fixture.Observe(), fixture,
        expect_allocation_free, label.c_str());
}

void require_raw_transient_failure(
    const DecodeObservation& observation,
    leo2_result expected,
    const RawTransientFixture& fixture,
    const char* operation)
{
    require_result(observation.result, expected, operation);
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
    require(observation.allocations == 0,
        "raw transient validation failure allocated");
#endif
    require(fixture.OutputsUntouched(),
        "raw transient validation failure modified an output");
}

void test_raw_transient_failure_atomicity(
    leo2_context* context,
    uint32_t k,
    uint32_t r,
    uint32_t losses,
    size_t bytes)
{
    RawTransientFixture fixture(context, k, r, bytes);
    fixture.Configure(losses, 2, true);
    require(fixture.decode_scratch_bytes > 1,
        "raw transient fixture has no decode scratch");

    require_raw_transient_failure(
        fixture.Observe(fixture.decode_scratch->data(),
            fixture.decode_scratch_bytes - 1U),
        LEO2_SCRATCH_TOO_SMALL, fixture,
        "raw transient short scratch");

    AlignedBuffer misaligned(fixture.decode_scratch_bytes + 64U);
    fixture.ResetOutputs();
    require_raw_transient_failure(fixture.Observe(
        static_cast<uint8_t*>(misaligned.data()) + 1,
        fixture.decode_scratch_bytes), LEO2_BAD_ALIGNMENT, fixture,
        "raw transient misaligned scratch");

    const std::vector<uint32_t> missing = fixture.Missing();
    require(missing.size() == losses,
        "raw transient missing fixture mismatch");
    const uint32_t missing0 = missing[0];
    const uint32_t missing1 = missing[1];
    uint32_t survivor = 0;
    while (!fixture.original_present[survivor])
        ++survivor;

    fixture.ResetOutputs();
    fixture.original[missing0] = fixture.originals[missing0].data();
    require_raw_transient_failure(fixture.Observe(),
        LEO2_INVALID_ARGUMENT, fixture,
        "raw transient presence pointer mismatch");
    fixture.original[missing0] = NULL;

    fixture.ResetOutputs();
    fixture.restored[missing0] = NULL;
    require_raw_transient_failure(fixture.Observe(),
        LEO2_INVALID_ARGUMENT, fixture,
        "raw transient null restored output");
    fixture.restored[missing0] = fixture.restored_storage[missing0].data();

    fixture.ResetOutputs();
    fixture.restored[missing0] = const_cast<void*>(fixture.original[survivor]);
    require_raw_transient_failure(fixture.Observe(),
        LEO2_OVERLAP, fixture, "raw transient output input overlap");
    fixture.restored[missing0] = fixture.restored_storage[missing0].data();

    fixture.ResetOutputs();
    fixture.restored[missing1] = fixture.restored[missing0];
    require_raw_transient_failure(fixture.Observe(),
        LEO2_OVERLAP, fixture, "raw transient output output overlap");
    fixture.restored[missing1] = fixture.restored_storage[missing1].data();

    // Every immutable metadata array is protected independently.  Exercise
    // the raw validator's complete five-span contract on this eligible route.
    fixture.ResetOutputs();
    fixture.restored[missing0] = fixture.original.data();
    require_raw_transient_failure(fixture.Observe(),
        LEO2_OVERLAP, fixture,
        "raw transient output original-pointer metadata overlap");
    fixture.restored[missing0] = fixture.restored_storage[missing0].data();

    fixture.ResetOutputs();
    fixture.restored[missing0] = fixture.recovery.data();
    require_raw_transient_failure(fixture.Observe(),
        LEO2_OVERLAP, fixture,
        "raw transient output recovery-pointer metadata overlap");
    fixture.restored[missing0] = fixture.restored_storage[missing0].data();

    fixture.ResetOutputs();
    fixture.restored[missing0] = fixture.restored.data();
    require_raw_transient_failure(fixture.Observe(),
        LEO2_OVERLAP, fixture,
        "raw transient output restored-pointer metadata overlap");
    fixture.restored[missing0] = fixture.restored_storage[missing0].data();

    fixture.ResetOutputs();
    fixture.restored[missing0] = fixture.original_present.data();
    require_raw_transient_failure(fixture.Observe(),
        LEO2_OVERLAP, fixture,
        "raw transient output original-presence metadata overlap");
    fixture.restored[missing0] = fixture.restored_storage[missing0].data();

    fixture.ResetOutputs();
    fixture.restored[missing0] = fixture.recovery_present.data();
    require_raw_transient_failure(fixture.Observe(),
        LEO2_OVERLAP, fixture,
        "raw transient output recovery-presence metadata overlap");
    fixture.restored[missing0] = fixture.restored_storage[missing0].data();

    fixture.ResetOutputs();
    const void* const saved_original = fixture.original[survivor];
    fixture.original[survivor] = reinterpret_cast<const void*>(
        std::numeric_limits<uintptr_t>::max() - fixture.bytes / 2U);
    require_raw_transient_failure(fixture.Observe(),
        LEO2_INVALID_ARGUMENT, fixture,
        "raw transient overflowing original shard range");
    fixture.original[survivor] = saved_original;

    uint32_t recovery_survivor = 0;
    while (!fixture.recovery_present[recovery_survivor])
        ++recovery_survivor;
    fixture.ResetOutputs();
    const void* const saved_recovery = fixture.recovery[recovery_survivor];
    fixture.recovery[recovery_survivor] = reinterpret_cast<const void*>(
        std::numeric_limits<uintptr_t>::max() - fixture.bytes / 2U);
    require_raw_transient_failure(fixture.Observe(),
        LEO2_INVALID_ARGUMENT, fixture,
        "raw transient overflowing recovery shard range");
    fixture.recovery[recovery_survivor] = saved_recovery;

    fixture.ResetOutputs();
    const uint8_t saved_recovery_presence = fixture.recovery_present.back();
    fixture.recovery_present.back() = 2;
    require_raw_transient_failure(fixture.Observe(),
        LEO2_INVALID_ARGUMENT, fixture,
        "raw transient invalid late presence");
    fixture.recovery_present.back() = saved_recovery_presence;

    fixture.Configure(losses, 2, false);
    std::fill(fixture.recovery_present.begin(),
        fixture.recovery_present.end(), 0);
    std::fill(fixture.recovery.begin(), fixture.recovery.end(),
        static_cast<const void*>(NULL));
    require_raw_transient_failure(fixture.Observe(),
        LEO2_NEED_MORE_DATA, fixture,
        "raw transient insufficient data");

    fixture.Configure(losses, 2, true);
    AlignedBuffer presence_scratch(fixture.decode_scratch_bytes);
    memcpy(presence_scratch.data(), fixture.original_present.data(), fixture.k);
    fixture.ResetOutputs();
    begin_allocation_audit();
    const leo2_result overlap_result = leo2_decode(fixture.codec,
        fixture.bytes, static_cast<const uint8_t*>(presence_scratch.data()),
        fixture.recovery_present.data(), fixture.original.data(),
        fixture.recovery.data(), fixture.restored.data(),
        presence_scratch.data(), presence_scratch.size());
    const uint64_t overlap_allocations = end_allocation_audit();
    require_result(overlap_result, LEO2_OVERLAP,
        "raw transient scratch presence overlap");
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
    require(overlap_allocations == 0,
        "raw transient scratch metadata overlap allocated");
#endif
    require(fixture.OutputsUntouched(),
        "raw transient scratch metadata overlap modified output");
}

void test_raw_transient_reusable_plan(
    leo2_context* context,
    uint32_t k,
    uint32_t r,
    uint32_t losses)
{
    RawTransientFixture fixture(context, k, r, 255);
    fixture.Configure(losses, 2, true);
    leo2_decode_plan* plan = NULL;
    require_result(leo2_decode_plan_create(fixture.codec,
        fixture.original_present.data(), fixture.recovery_present.data(),
        &plan), LEO2_SUCCESS, "raw transient reusable plan create");
    size_t scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(
        plan, fixture.bytes, &scratch_bytes),
        LEO2_SUCCESS, "raw transient reusable plan scratch query");
    AlignedBuffer scratch(scratch_bytes);
    fixture.ResetOutputs();
    begin_allocation_audit();
    const leo2_result result = leo2_decode_plan_execute(plan, fixture.bytes,
        fixture.original.data(), fixture.recovery.data(),
        fixture.restored.data(), scratch.data(), scratch.size());
    const uint64_t allocations = end_allocation_audit();
    require_result(result, LEO2_SUCCESS,
        "raw transient reusable plan execute");
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
    require(allocations == 0,
        "reusable plan execution allocated after raw transient refactor");
#endif
    require(fixture.OutputsMatch(),
        "reusable plan diverged from raw transient one-shot decode");
    leo2_decode_plan_destroy(plan);
}

uint64_t test_raw_native_high_matrix(leo2_context* avx2)
{
    std::vector<size_t> boundary_bytes;
    boundary_bytes.push_back(1);
    boundary_bytes.push_back(63);
    boundary_bytes.push_back(64);
    for (size_t tail = 1; tail < 64; ++tail)
        boundary_bytes.push_back(64 + tail);
    boundary_bytes.push_back(128);
    boundary_bytes.push_back(129);
    boundary_bytes.push_back(193);
    boundary_bytes.push_back(255);
    boundary_bytes.push_back(256);
    const size_t extended_bytes[] = {
        257,
        319, 320, 321,
        511, 512, 513,
        1023, 1024, 1025,
        2047, 2048, 2049,
        4095, 4096, 4097,
        6143, 6144, 6145,
        7103, 7104, 7105,
        7135, 7136, 7137,
        7167, 7168
    };
    boundary_bytes.insert(boundary_bytes.end(),
        extended_bytes,
        extended_bytes +
            sizeof(extended_bytes) / sizeof(extended_bytes[0]));
    // The allocation-free source-major interval begins on a vector boundary.
    // Exercise every 64-byte residue without letting a padded-tail bug hide
    // behind naturally aligned benchmark sizes.
    for (size_t tail = 0; tail < 64; ++tail)
        boundary_bytes.push_back(12288 + tail);
    boundary_bytes.push_back(16383);
    boundary_bytes.push_back(16384);
    uint64_t cases = 0;
    for (uint32_t k = 9; k <= 16; ++k)
    {
        for (uint32_t r = 5; r <= 8; ++r)
        {
            for (size_t byte_i = 0;
                 byte_i < boundary_bytes.size(); ++byte_i)
            {
                RawTransientFixture fixture(
                    avx2, k, r, boundary_bytes[byte_i]);
                for (uint32_t losses = 3; losses <= r; ++losses)
                {
                    const unsigned original_shape = static_cast<unsigned>(
                        (k + r + losses + byte_i) % 3U);
                    fixture.Configure(
                        losses, original_shape, false);
                    require_raw_transient_success(
                        fixture.Observe(), fixture, true,
                        "raw native-high boundary decode");
                    ++cases;

                    if (losses < r)
                    {
                        fixture.Configure(
                            losses, (original_shape + 1U) % 3U, false);
                        fixture.SelectRecoverySubset(
                            losses, (original_shape + 2U) % 3U);
                        require_raw_transient_success(
                            fixture.Observe(), fixture, true,
                            "raw native-high parity-loss decode");
                        ++cases;
                    }
                }
            }
        }
    }
    return cases;
}

void test_raw_transient_decode(leo2_context* automatic_context)
{
    leo2_context_options options;
    memset(&options, 0, sizeof(options));
    options.struct_size = sizeof(options);
    options.backend = LEO2_BACKEND_AVX2;
    options.thread_count = 1;
    leo2_context* avx2 = NULL;
    const leo2_result context_result = leo2_context_create(&options, &avx2);
    if (context_result == LEO2_UNSUPPORTED)
        return;
    require_result(context_result, LEO2_SUCCESS,
        "raw transient AVX2 context create");
    require(leo2_context_backend(avx2) == LEO2_BACKEND_AVX2,
        "raw transient explicit context did not select AVX2");

    const size_t boundary_bytes[] = { 1, 63, 64, 65, 255, 256 };
    for (size_t i = 0;
         i < sizeof(boundary_bytes) / sizeof(boundary_bytes[0]); ++i)
    {
        run_raw_transient_case(avx2, 32, 32, boundary_bytes[i],
            9, static_cast<unsigned>(i % 3), true, true);
    }
    // Exercise the smallest reachable translated-low parent above the
    // eight-loss direct-repair limit, then both sides of the P=32 count
    // boundary.  These are eligibility boundaries, not merely larger-shape
    // repetitions of the P=32/P=128 cases below.
    run_raw_transient_case(avx2, 9, 9, 1, 9, 0, true, true);
    run_raw_transient_case(avx2, 16, 16, 256, 9, 1, true, true);
    run_raw_transient_case(avx2, 17, 17, 64, 9, 2, true, true);
    run_raw_transient_case(avx2, 17, 17, 255, 17, 0, false, true);
    run_raw_transient_case(avx2, 95, 95, 65, 47, 2, true, true);
    run_raw_transient_case(avx2, 95, 95, 256, 95, 1, false, true);
    run_raw_transient_case(avx2, 17, 31, 255, 17, 0, true, true);
    run_raw_transient_case(avx2, 31, 17, 64, 17, 2, false, true);

    // Translated-low raw setup remains independently bounded to 256 bytes.
    run_raw_transient_case(avx2, 32, 32, 257, 9, 2, true, false);
    const uint64_t raw_native_high_cases =
        test_raw_native_high_matrix(avx2);
    require(raw_native_high_cases == 41984,
        "raw native-high matrix case count drifted");

    /*
        Dense partial-loss T=64 calls use a separate scratch-owned Algorithm 5
        view at exactly 64 bytes.  Cover the K/R/loss eligibility boundaries,
        all three systematic-loss shapes, and surplus/non-prefix recovery
        subsets.  Allocation audit is the route oracle: a canonical transient
        plan still allocates, while the promoted terminal does not.
    */
    run_raw_transient_case(avx2, 65, 33, 64, 16, 0, true, true);
    run_raw_transient_case(avx2, 99, 50, 64, 25, 2, true, true);
    run_raw_transient_case(avx2, 99, 50, 64, 49, 1, true, true);
    run_raw_transient_case(avx2, 124, 63, 64, 31, 0, true, true);
    run_raw_transient_case(avx2, 124, 64, 64, 63, 2, true, true);
    run_raw_transient_case(avx2, 191, 62, 64, 61, 1, true, true);
    run_raw_transient_case(avx2, 192, 33, 64, 16, 0, true, true);
    run_raw_transient_case(avx2, 192, 34, 64, 17, 1, true, true);
    run_raw_transient_case(avx2, 192, 62, 64, 31, 2, true, true);
    run_raw_transient_case(avx2, 192, 64, 64, 63, 0, true, true);

    /* K=64 is outside this terminal but remains allocation-free through the
       independently qualified translated-low raw path. */
    run_raw_transient_case(avx2, 64, 33, 64, 16, 0, true, true);
    run_raw_transient_case(avx2, 99, 50, 63, 25, 2, true, false);
    run_raw_transient_case(avx2, 99, 50, 65, 25, 0, true, false);
    run_raw_transient_case(avx2, 99, 50, 64, 24, 1, true, false);
    run_raw_transient_case(avx2, 99, 50, 64, 50, 2, false, false);
    run_raw_transient_case(avx2, 192, 62, 64, 30, 1, true, false);
    run_raw_transient_case(avx2, 192, 62, 64, 62, 2, false, false);

    // Native-high direct execution has separately measured output-major and
    // source-major intervals.  Their gap and the first byte above the latter
    // retain the heap-owned direct-plan fallback.
    run_raw_transient_case(avx2, 16, 8, 257, 8, 2, false, true);
    run_raw_transient_case(avx2, 16, 8, 7168, 8, 1, true, true);
    run_raw_transient_case(avx2, 16, 8, 7169, 8, 2, false, false);
    run_raw_transient_case(avx2, 16, 8, 12287, 8, 0, true, false);
    run_raw_transient_case(avx2, 16, 8, 12288, 8, 1, true, true);
    run_raw_transient_case(avx2, 16, 8, 12351, 8, 2, true, true);
    run_raw_transient_case(avx2, 16, 8, 16384, 8, 0, true, true);
    run_raw_transient_case(avx2, 16, 8, 16385, 8, 1, true, false);
    run_raw_transient_case(avx2, 16, 8, 32768, 8, 2, true, false);
    run_raw_transient_case(avx2, 16, 8, 65536, 8, 0, true, false);
    run_raw_transient_case(avx2, 16, 8, 1024U * 1024U,
        8, 1, true, false);
    test_raw_transient_failure_atomicity(avx2, 32, 32, 9, 65);
    test_raw_transient_failure_atomicity(avx2, 16, 8, 8, 7168);
    test_raw_transient_failure_atomicity(avx2, 16, 8, 8, 16384);
    test_raw_transient_failure_atomicity(avx2, 99, 50, 49, 64);
    test_raw_transient_failure_atomicity(avx2, 192, 62, 31, 64);
    test_raw_transient_reusable_plan(avx2, 32, 32, 9);
    test_raw_transient_reusable_plan(avx2, 16, 8, 8);

    {
        RawTransientFixture fixture(avx2, 16, 8, 64,
            LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8);
        fixture.Configure(8, 1, false);
        require_raw_transient_success(fixture.Observe(), fixture, false,
            "raw native-high low-profile fallback");
    }
#if LEO2_TEST_EXPECT_GF16
    {
        RawTransientFixture fixture(avx2, 16, 8, 64,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16);
        fixture.Configure(8, 1, false);
        require_raw_transient_success(fixture.Observe(), fixture, false,
            "raw native-high GF16 fallback");
    }
#endif
    {
        RawTransientFixture fixture(avx2, 16, 8, 64,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
            LEO2_CODEC_FORCE_TILED_DECODE);
        fixture.Configure(8, 1, false);
        require_raw_transient_success(fixture.Observe(), fixture, false,
            "raw native-high forced-workspace fallback");
    }

#if defined(LEO2_TEST_ONE_SHOT_NO_LOSS_HOOK_BUILD)
    {
        RawTransientFixture fixture(avx2, 32, 32, 65);
        fixture.Configure(9, 2, true);
        require_result(leo2_test_codec_set_decode_mode(fixture.codec,
            LEO2_TEST_DECODE_FORCE_TRANSLATED_LOW),
            LEO2_SUCCESS, "raw transient force translated-low mode");
        require_raw_transient_success(fixture.Observe(), fixture, false,
            "forced mode raw transient fallback");
        require_result(leo2_test_codec_set_decode_mode(fixture.codec,
            LEO2_TEST_DECODE_AUTO), LEO2_SUCCESS,
            "raw transient restore automatic mode");
    }
    {
        require(leopard2_internal::SetOneShotPlanSetupModeForDiagnostics(2),
            "raw native-high disable mode failed");
        run_raw_transient_case(avx2, 16, 8, 64, 8, 2, false, false);
        require(leopard2_internal::SetOneShotPlanSetupModeForDiagnostics(3),
            "raw native-high production mode restore failed");
    }
#endif

    if (leo2_context_backend(automatic_context) == LEO2_BACKEND_AVX2)
    {
        run_raw_transient_case(automatic_context,
            32, 32, 65, 9, 2, true, true);
        run_raw_transient_case(automatic_context,
            16, 8, 65, 8, 2, false, true);
    }

    options.backend = LEO2_BACKEND_SCALAR;
    leo2_context* scalar = NULL;
    require_result(leo2_context_create(&options, &scalar), LEO2_SUCCESS,
        "raw native-high scalar context create");
    run_raw_transient_case(
        scalar, 16, 8, 64, 8, 2, false, false);
    leo2_context_destroy(scalar);
    leo2_context_destroy(avx2);
}
#endif

void test_no_loss_shape(
    leo2_context* context,
    uint32_t k,
    uint32_t r,
    leo2_profile profile,
    leo2_field field,
    uint64_t odd_or_zero_bytes)
{
    leo2_codec* const codec = create_codec(
        context, k, r, profile, field);
    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 0);

    DecodeObservation observation = observe_decode(codec, 0,
        original_present.data(), recovery_present.data(), NULL, NULL, NULL,
        NULL, 0);
    require_no_loss_observation(
        observation, "zero-byte missing-parity no-loss decode");
    if (odd_or_zero_bytes != 0)
    {
        observation = observe_decode(codec, odd_or_zero_bytes,
            original_present.data(), recovery_present.data(), NULL, NULL,
            NULL, NULL, 0);
        require_no_loss_observation(
            observation, "odd-byte missing-parity no-loss decode");
    }

    std::fill(recovery_present.begin(), recovery_present.end(), 1);
    observation = observe_decode(codec,
        std::numeric_limits<uint64_t>::max(), original_present.data(),
        recovery_present.data(), NULL, NULL, NULL, NULL, 0);
    require_no_loss_observation(
        observation, "surplus-parity no-loss decode");

    for (uint32_t i = 0; i < r; ++i)
        recovery_present[i] = static_cast<uint8_t>(i & 1u);
    observation = observe_decode(codec, odd_or_zero_bytes,
        original_present.data(), recovery_present.data(), NULL, NULL, NULL,
        NULL, 0);
    require_no_loss_observation(
        observation, "mixed-parity no-loss decode");

    /*
        No-loss execution is a true no-op after presence validation.  Supply
        data that would fail ordinary alias and scratch checks and prove it is
        neither inspected nor modified.
    */
    uint8_t payload = 0x5a;
    std::vector<const void*> original(k, &payload);
    std::vector<const void*> recovery(r, &payload);
    std::vector<void*> restored(k, original_present.data());
    const std::vector<uint8_t> original_before = original_present;
    const std::vector<uint8_t> recovery_before = recovery_present;
    observation = observe_decode(codec,
        std::numeric_limits<uint64_t>::max(), original_present.data(),
        recovery_present.data(), original.data(), recovery.data(),
        restored.data(), original_present.data(),
        std::numeric_limits<size_t>::max());
    require_no_loss_observation(
        observation, "alias-heavy no-loss decode");
    require(original_present == original_before &&
            recovery_present == recovery_before && payload == 0x5a,
        "one-shot no-loss decode modified caller storage");

    leo2_codec_destroy(codec);
}

void require_preplan_failure(
    const DecodeObservation& observation,
    leo2_result expected,
    const char* operation)
{
    require_result(observation.result, expected, operation);
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
    require(observation.allocations == 0,
        "presence validation allocated before rejecting the call");
#else
    (void)observation;
#endif
}

void test_validation_and_failure_atomicity(
    leo2_context* context,
    leo2_field field)
{
    const uint32_t k = 9;
    const uint32_t r = 7;
    leo2_codec* const codec = create_codec(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, field);
    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 0);
    uint8_t output[17];
    memset(output, 0xa5, sizeof(output));
    void* restored[9] = { output, output, output, output, output,
                          output, output, output, output };

    original_present.back() = 2;
    require_preplan_failure(observe_decode(codec, sizeof(output),
        original_present.data(), recovery_present.data(), NULL, NULL,
        restored, output, sizeof(output)), LEO2_INVALID_ARGUMENT,
        "late invalid original presence");
    original_present.back() = 1;

    recovery_present.back() = 2;
    require_preplan_failure(observe_decode(codec, sizeof(output),
        original_present.data(), recovery_present.data(), NULL, NULL,
        restored, output, sizeof(output)), LEO2_INVALID_ARGUMENT,
        "late invalid recovery presence");
    recovery_present.back() = 0;

    /*
        The speculative scan stops at the first loss.  The canonical fallback
        must still validate every later presence byte.
    */
    original_present.front() = 0;
    original_present.back() = 2;
    require_preplan_failure(observe_decode(codec, sizeof(output),
        original_present.data(), recovery_present.data(), NULL, NULL,
        restored, output, sizeof(output)), LEO2_INVALID_ARGUMENT,
        "loss-prefix late invalid original presence");
    original_present.back() = 1;
    recovery_present.back() = 2;
    require_preplan_failure(observe_decode(codec, sizeof(output),
        original_present.data(), recovery_present.data(), NULL, NULL,
        restored, output, sizeof(output)), LEO2_INVALID_ARGUMENT,
        "loss-prefix late invalid recovery presence");
    recovery_present.back() = 0;

    require_preplan_failure(observe_decode(codec, sizeof(output),
        original_present.data(), recovery_present.data(), NULL, NULL,
        restored, output, sizeof(output)), LEO2_NEED_MORE_DATA,
        "loss-prefix insufficient data");
    original_present.front() = 1;

    const uintptr_t address_end =
        std::numeric_limits<uintptr_t>::max() - 1;
    const uint8_t* const overflowing =
        reinterpret_cast<const uint8_t*>(address_end);
    require_preplan_failure(observe_decode(codec, sizeof(output), overflowing,
        recovery_present.data(), NULL, NULL, restored, output,
        sizeof(output)), LEO2_INVALID_ARGUMENT,
        "overflowing original-presence span");
    require_preplan_failure(observe_decode(codec, sizeof(output),
        original_present.data(), overflowing, NULL, NULL, restored, output,
        sizeof(output)), LEO2_INVALID_ARGUMENT,
        "overflowing recovery-presence span");
    require_preplan_failure(observe_decode(NULL, sizeof(output),
        original_present.data(), recovery_present.data(), NULL, NULL,
        restored, output, sizeof(output)), LEO2_INVALID_ARGUMENT,
        "null codec");
    require_preplan_failure(observe_decode(codec, sizeof(output), NULL,
        recovery_present.data(), NULL, NULL, restored, output,
        sizeof(output)), LEO2_INVALID_ARGUMENT,
        "null original presence");
    require_preplan_failure(observe_decode(codec, sizeof(output),
        original_present.data(), NULL, NULL, NULL, restored, output,
        sizeof(output)), LEO2_INVALID_ARGUMENT,
        "null recovery presence");

    for (size_t i = 0; i < sizeof(output); ++i)
        require(output[i] == 0xa5,
            "rejected one-shot decode modified output or scratch");
    leo2_codec_destroy(codec);
}

void test_r1_presence_validation(
    leo2_context* context,
    leo2_field field)
{
    const uint32_t k = 9;
    leo2_codec* const codec = create_codec(context, k, 1,
        LEO2_PROFILE_LEGACY_HIGH_V1, field);
    std::vector<uint8_t> original_present(k, 1);
    uint8_t recovery_present[1] = { 0 };

    recovery_present[0] = 2;
    require_preplan_failure(observe_decode(codec, 0,
        original_present.data(), recovery_present, NULL, NULL, NULL, NULL, 0),
        LEO2_INVALID_ARGUMENT, "R=1 invalid recovery presence after no loss");

    recovery_present[0] = 0;
    original_present.back() = 2;
    require_preplan_failure(observe_decode(codec, 0,
        original_present.data(), recovery_present, NULL, NULL, NULL, NULL, 0),
        LEO2_INVALID_ARGUMENT, "R=1 late invalid original presence");

    recovery_present[0] = 1;
    original_present.front() = 0;
    require_preplan_failure(observe_decode(codec, 0,
        original_present.data(), recovery_present, NULL, NULL, NULL, NULL, 0),
        LEO2_INVALID_ARGUMENT,
        "R=1 invalid original presence after an earlier loss");

    original_present.back() = 1;
    recovery_present[0] = 2;
    require_preplan_failure(observe_decode(codec, 0,
        original_present.data(), recovery_present, NULL, NULL, NULL, NULL, 0),
        LEO2_INVALID_ARGUMENT,
        "R=1 invalid recovery presence after an earlier loss");

    leo2_codec_destroy(codec);
}

void test_loss_fallback(
    leo2_context* context,
    leo2_field field)
{
    const uint32_t k = 4;
    const uint32_t r = 3;
    const size_t bytes = 64;
    leo2_codec* const codec = create_codec(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, field);
    std::vector<std::vector<uint8_t> > originals(
        k, std::vector<uint8_t>(bytes));
    std::vector<std::vector<uint8_t> > parity(
        r, std::vector<uint8_t>(bytes));
    for (uint32_t shard = 0; shard < k; ++shard)
        for (size_t i = 0; i < bytes; ++i)
            originals[shard][i] = static_cast<uint8_t>(
                shard * 53u + i * 29u + 7u);
    const void* original_input[k];
    void* parity_output[r];
    for (uint32_t i = 0; i < k; ++i)
        original_input[i] = originals[i].data();
    for (uint32_t i = 0; i < r; ++i)
        parity_output[i] = parity[i].data();
    size_t encode_scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(
        codec, bytes, &encode_scratch_bytes),
        LEO2_SUCCESS, "encode scratch query");
    AlignedBuffer encode_scratch(encode_scratch_bytes);
    require_result(leo2_encode(codec, bytes, original_input, parity_output,
        encode_scratch.data(), encode_scratch.size()),
        LEO2_SUCCESS, "loss fallback fixture encode");

    size_t decode_scratch_bytes = 0;
    require_result(leo2_decode_scratch_size(
        codec, bytes, &decode_scratch_bytes),
        LEO2_SUCCESS, "decode scratch query");
    AlignedBuffer decode_scratch(decode_scratch_bytes);
    std::vector<uint8_t> restored_bytes(bytes, 0xa5);
    for (uint32_t missing_case = 0; missing_case < 2; ++missing_case)
    {
        const uint32_t missing = missing_case == 0 ? 0 : k - 1;
        uint8_t original_present[k];
        uint8_t recovery_present[r];
        const void* decode_original[k];
        const void* decode_recovery[r];
        void* restored[k];
        for (uint32_t i = 0; i < k; ++i)
        {
            original_present[i] = i == missing ? 0 : 1;
            decode_original[i] =
                i == missing ? NULL : originals[i].data();
            restored[i] =
                i == missing ? restored_bytes.data() : NULL;
        }
        for (uint32_t i = 0; i < r; ++i)
        {
            recovery_present[i] = 1;
            decode_recovery[i] = parity[i].data();
        }
        std::fill(restored_bytes.begin(), restored_bytes.end(), 0xa5);
        const DecodeObservation observation = observe_decode(
            codec, bytes, original_present, recovery_present,
            decode_original, decode_recovery, restored,
            decode_scratch.data(), decode_scratch.size());
        require_result(observation.result, LEO2_SUCCESS,
            "one-loss canonical fallback");
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
        require(observation.allocations > 0,
            "one-loss call bypassed canonical plan construction");
#endif
        require(restored_bytes == originals[missing],
            "one-loss canonical fallback restored wrong bytes");
    }
    leo2_codec_destroy(codec);
}

#if defined(LEO2_TEST_ONE_SHOT_NO_LOSS_HOOK_BUILD)
void test_forced_hook_modes(
    leo2_context* context,
    leo2_field field)
{
    leo2_codec* const codec = create_codec(context, 8, 8,
        LEO2_PROFILE_LEGACY_HIGH_V1, field);
    uint8_t original_present[8];
    uint8_t recovery_present[8];
    memset(original_present, 1, sizeof(original_present));
    memset(recovery_present, 0, sizeof(recovery_present));
    const leo2_test_decode_mode modes[] = {
        LEO2_TEST_DECODE_FORCE_TRANSLATED_LOW,
        LEO2_TEST_DECODE_FORCE_NATIVE_HIGH
    };
    for (size_t i = 0; i < sizeof(modes) / sizeof(modes[0]); ++i)
    {
        require_result(leo2_test_codec_set_decode_mode(codec, modes[i]),
            LEO2_SUCCESS, "select forced decode mode");
        const DecodeObservation observation = observe_decode(codec, 17,
            original_present, recovery_present, NULL, NULL, NULL, NULL, 0);
        require_result(observation.result, LEO2_SUCCESS,
            "forced-mode one-shot no-loss decode");
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
        require(observation.allocations > 0,
            "forced hook mode bypassed canonical plan construction");
#endif
    }
    leo2_codec_destroy(codec);
}
#endif

} // namespace

int main()
{
    try
    {
        leo2_context_options options;
        memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AUTO;
        options.thread_count = 1;
        leo2_context* context = NULL;
        require_result(leo2_context_create(&options, &context),
            LEO2_SUCCESS, "context create");

#if LEO2_TEST_EXPECT_GF8
        test_no_loss_shape(context, 9, 7,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 0);
        test_no_loss_shape(context, 7, 9,
            LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8, 17);
        test_no_loss_shape(context, 65, 1,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 17);
        test_validation_and_failure_atomicity(
            context, LEO2_FIELD_GF8);
        test_r1_presence_validation(context, LEO2_FIELD_GF8);
        test_loss_fallback(context, LEO2_FIELD_GF8);
        test_raw_transient_decode(context);
#endif
#if LEO2_TEST_EXPECT_GF16
        test_no_loss_shape(context, 9, 7,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, 65);
        test_no_loss_shape(context, 7, 9,
            LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16, 65);
        test_validation_and_failure_atomicity(
            context, LEO2_FIELD_GF16);
        test_r1_presence_validation(context, LEO2_FIELD_GF16);
        test_loss_fallback(context, LEO2_FIELD_GF16);
#endif

#if defined(LEO2_TEST_ONE_SHOT_NO_LOSS_HOOK_BUILD)
# if LEO2_TEST_EXPECT_GF8
        test_forced_hook_modes(context, LEO2_FIELD_GF8);
# elif LEO2_TEST_EXPECT_GF16
        test_forced_hook_modes(context, LEO2_FIELD_GF16);
# endif
#endif
        leo2_context_destroy(context);

        std::cout << "{\"one_shot_no_loss_short_circuit\":"
                  << LEO2_EXPERIMENT_ONE_SHOT_NO_LOSS_SHORT_CIRCUIT
                  << ",\"hook_build\":"
#if defined(LEO2_TEST_ONE_SHOT_NO_LOSS_HOOK_BUILD)
                  << 1
#else
                  << 0
#endif
                  << ",\"allocation_audit\":"
                  << LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
                  << ",\"status\":\"ok\"}" << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << std::endl;
        return 1;
    }
}
