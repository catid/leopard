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
#include "Leopard2Direct.h"
#include "allocation_audit_config.h"
#include "direct_oracle.h"

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

namespace {

static const unsigned kOriginalCount = 8;
static const unsigned kRecoveryCount = 4;
static const uint8_t kGuard = 0xa5;
static std::atomic<bool> g_track_allocations(false);
static std::atomic<uint64_t> g_tracked_allocations(0);

#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
#if defined(_MSC_VER)
#define LEO2_TEST_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define LEO2_TEST_NOINLINE __attribute__((noinline))
#else
#define LEO2_TEST_NOINLINE
#endif

} // namespace

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

#if defined(__cpp_sized_deallocation)
LEO2_TEST_NOINLINE void operator delete(void* pointer, size_t) noexcept
{
    free(pointer);
}

LEO2_TEST_NOINLINE void operator delete[](void* pointer, size_t) noexcept
{
    free(pointer);
}
#endif

#undef LEO2_TEST_NOINLINE

namespace {
#endif

void Require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void RequireResult(
    leo2_result actual,
    leo2_result expected,
    const std::string& label)
{
    if (actual == expected)
        return;
    throw std::runtime_error(label + ": expected " +
        leo2_result_string(expected) + ", got " +
        leo2_result_string(actual));
}

class AlignedBuffer
{
public:
    explicit AlignedBuffer(size_t bytes)
        : pointer_(NULL)
        , bytes_(bytes)
    {
        const size_t allocated = std::max(bytes, leo2_scratch_alignment());
#if defined(_MSC_VER)
        pointer_ = _aligned_malloc(allocated, leo2_scratch_alignment());
#else
        if (posix_memalign(
                &pointer_, leo2_scratch_alignment(), allocated) != 0)
            pointer_ = NULL;
#endif
        if (!pointer_)
            throw std::bad_alloc();
        std::memset(pointer_, 0, allocated);
    }

    ~AlignedBuffer()
    {
#if defined(_MSC_VER)
        _aligned_free(pointer_);
#else
        free(pointer_);
#endif
    }

    void* data() { return pointer_; }
    const void* data() const { return pointer_; }
    size_t size() const { return bytes_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);

    void* pointer_;
    size_t bytes_;
};

class Context
{
public:
    explicit Context(leo2_backend backend)
        : pointer_(NULL)
        , result_(LEO2_INTERNAL_ERROR)
    {
        leo2_context_options options = {};
        options.struct_size = sizeof(options);
        options.backend = backend;
        options.thread_count = 1;
        result_ = leo2_context_create(&options, &pointer_);
    }

    ~Context() { leo2_context_destroy(pointer_); }
    leo2_context* get() const { return pointer_; }
    leo2_result result() const { return result_; }

private:
    Context(const Context&);
    Context& operator=(const Context&);

    leo2_context* pointer_;
    leo2_result result_;
};

class Codec
{
public:
    Codec(
        leo2_context* context,
        unsigned original_count,
        unsigned recovery_count,
        leo2_profile profile = LEO2_PROFILE_LEGACY_HIGH_V1,
        uint32_t flags = 0)
        : pointer_(NULL)
    {
        leo2_codec_options options = {};
        options.struct_size = sizeof(options);
        options.flags = flags;
        RequireResult(leo2_codec_create(context,
            original_count, recovery_count, profile, LEO2_FIELD_GF8,
            &options, &pointer_), LEO2_SUCCESS, "create codec");
    }

    ~Codec() { leo2_codec_destroy(pointer_); }
    leo2_codec* get() const { return pointer_; }

private:
    Codec(const Codec&);
    Codec& operator=(const Codec&);

    leo2_codec* pointer_;
};

unsigned CountBits(unsigned value)
{
    unsigned count = 0;
    while (value != 0)
    {
        count += value & 1U;
        value >>= 1;
    }
    return count;
}

struct Fixture
{
    explicit Fixture(leo2_codec* codec, size_t bytes, uint64_t seed)
        : shard_bytes(bytes)
        , stride(bytes + 11)
        , original_base(1)
        , recovery_base(3)
        , restored_base(5)
        , original_storage(original_base + kOriginalCount * stride + 13,
              kGuard)
        , recovery_storage(recovery_base + kRecoveryCount * stride + 13,
              kGuard)
        , restored_storage(restored_base + kOriginalCount * stride + 13,
              kGuard)
        , original(kOriginalCount)
        , recovery(kRecoveryCount)
        , recovery_output(kRecoveryCount)
        , restored(kOriginalCount, NULL)
        , original_present(kOriginalCount, 1)
        , recovery_present(kRecoveryCount, 1)
        , encode_scratch_bytes(0)
        , decode_scratch_bytes(0)
        , encode_scratch(QueryEncodeScratch(codec, bytes))
        , decode_scratch(QueryDecodeScratch(codec, bytes))
    {
        encode_scratch_bytes = encode_scratch.size();
        decode_scratch_bytes = decode_scratch.size();
        uint64_t state = seed ^ bytes;
        for (unsigned row = 0; row < kOriginalCount; ++row)
        {
            uint8_t* const shard = OriginalData(row);
            original[row] = shard;
            for (size_t offset = 0; offset < shard_bytes; ++offset)
            {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                shard[offset] = static_cast<uint8_t>(state >> 24);
            }
        }
        for (unsigned row = 0; row < kRecoveryCount; ++row)
        {
            recovery[row] = RecoveryData(row);
            recovery_output[row] = RecoveryData(row);
        }
    }

    static size_t QueryEncodeScratch(leo2_codec* codec, size_t bytes)
    {
        size_t result = 0;
        RequireResult(leo2_encode_scratch_size(codec, bytes, &result),
            LEO2_SUCCESS, "query encode scratch");
        return result;
    }

    static size_t QueryDecodeScratch(leo2_codec* codec, size_t bytes)
    {
        size_t result = 0;
        RequireResult(leo2_decode_scratch_size(codec, bytes, &result),
            LEO2_SUCCESS, "query decode scratch");
        return result;
    }

    uint8_t* OriginalData(unsigned row)
    {
        return &original_storage[original_base + row * stride];
    }

    uint8_t* RecoveryData(unsigned row)
    {
        return &recovery_storage[recovery_base + row * stride];
    }

    uint8_t* RestoredData(unsigned row)
    {
        return &restored_storage[restored_base + row * stride];
    }

    void Encode(leo2_codec* codec)
    {
        RequireResult(leo2_encode(codec, shard_bytes,
            &original[0], &recovery_output[0], encode_scratch.data(),
            encode_scratch_bytes), LEO2_SUCCESS, "encode fixture");
    }

    void SelectMissing(unsigned mask, bool in_place = false)
    {
        std::fill(original_present.begin(), original_present.end(), 1);
        std::fill(recovery_present.begin(), recovery_present.end(), 1);
        std::fill(restored.begin(), restored.end(), static_cast<void*>(NULL));
        std::fill(restored_storage.begin(), restored_storage.end(), kGuard);
        for (unsigned row = 0; row < kOriginalCount; ++row)
        {
            const bool missing = (mask & (1U << row)) != 0;
            original_present[row] = missing ? 0 : 1;
            original[row] = missing ? NULL : OriginalData(row);
            if (missing)
            {
                restored[row] = in_place
                    ? static_cast<void*>(OriginalData(row))
                    : static_cast<void*>(RestoredData(row));
            }
        }
    }

    size_t shard_bytes;
    size_t stride;
    size_t original_base;
    size_t recovery_base;
    size_t restored_base;
    std::vector<uint8_t> original_storage;
    std::vector<uint8_t> recovery_storage;
    std::vector<uint8_t> restored_storage;
    std::vector<const void*> original;
    std::vector<const void*> recovery;
    std::vector<void*> recovery_output;
    std::vector<void*> restored;
    std::vector<uint8_t> original_present;
    std::vector<uint8_t> recovery_present;
    size_t encode_scratch_bytes;
    size_t decode_scratch_bytes;
    AlignedBuffer encode_scratch;
    AlignedBuffer decode_scratch;
};

void CheckParity(
    const Fixture& fixture,
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator)
{
    std::vector<leopard2_test::Element> message(kOriginalCount);
    for (size_t offset = 0; offset < fixture.shard_bytes; ++offset)
    {
        for (unsigned original = 0;
             original < kOriginalCount; ++original)
        {
            message[original] = fixture.original_storage[
                fixture.original_base + original * fixture.stride + offset];
        }
        for (unsigned parity = 0; parity < kRecoveryCount; ++parity)
        {
            leopard2_test::Element expected = 0;
            for (unsigned original = 0;
                 original < kOriginalCount; ++original)
            {
                expected ^= field.multiply(
                    generator[kOriginalCount + parity][original],
                    message[original]);
            }
            const uint8_t actual = fixture.recovery_storage[
                fixture.recovery_base + parity * fixture.stride + offset];
            Require(actual == expected,
                "production parity differs from independent direct oracle");
        }
    }
}

uint64_t DecodeAndCountAllocations(leo2_codec* codec, Fixture& fixture)
{
    g_tracked_allocations.store(0, std::memory_order_relaxed);
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
    g_track_allocations.store(true, std::memory_order_release);
#endif
    const leo2_result result = leo2_decode(codec, fixture.shard_bytes,
        &fixture.original_present[0], &fixture.recovery_present[0],
        &fixture.original[0], &fixture.recovery[0], &fixture.restored[0],
        fixture.decode_scratch.data(), fixture.decode_scratch_bytes);
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
    g_track_allocations.store(false, std::memory_order_release);
#endif
    RequireResult(result, LEO2_SUCCESS, "one-shot decode");
    return g_tracked_allocations.load(std::memory_order_relaxed);
}

void CheckRestored(const Fixture& fixture, unsigned mask, bool in_place = false)
{
    for (unsigned row = 0; row < kOriginalCount; ++row)
    {
        if ((mask & (1U << row)) == 0)
            continue;
        const uint8_t* const output = in_place
            ? &fixture.original_storage[fixture.original_base + row * fixture.stride]
            : &fixture.restored_storage[fixture.restored_base + row * fixture.stride];
        const uint8_t* const expected =
            &fixture.original_storage[fixture.original_base + row * fixture.stride];
        Require(std::memcmp(output, expected, fixture.shard_bytes) == 0,
            "restored original differs from source");
    }
}

void CheckRestoredGuards(const Fixture& fixture, unsigned mask)
{
    for (size_t offset = 0; offset < fixture.restored_base; ++offset)
        Require(fixture.restored_storage[offset] == kGuard,
            "restored prefix guard changed");
    for (unsigned row = 0; row < kOriginalCount; ++row)
    {
        const size_t begin = fixture.restored_base + row * fixture.stride;
        if ((mask & (1U << row)) == 0)
        {
            for (size_t offset = 0; offset < fixture.shard_bytes; ++offset)
                Require(fixture.restored_storage[begin + offset] == kGuard,
                    "unrequested restored payload changed");
        }
        for (size_t offset = fixture.shard_bytes;
             offset < fixture.stride; ++offset)
        {
            Require(fixture.restored_storage[begin + offset] == kGuard,
                "restored row guard changed");
        }
    }
    const size_t used = fixture.restored_base +
        kOriginalCount * fixture.stride;
    for (size_t offset = used;
         offset < fixture.restored_storage.size(); ++offset)
    {
        Require(fixture.restored_storage[offset] == kGuard,
            "restored suffix guard changed");
    }
}

leopard2_internal::CodecDecodeMetadataInfo Metadata(leo2_codec* codec)
{
    leopard2_internal::CodecDecodeMetadataInfo info = {};
    Require(leopard2_internal::GetCodecDecodeMetadataInfo(codec, &info),
        "query codec metadata");
    return info;
}

void CheckCachePredicate(leo2_context* avx2_context)
{
    Codec target(avx2_context, 8, 4);
    Require(Metadata(target.get()).codec_direct_repair_generator_bytes == 0,
        "K=8/R=4 AVX2 codec duplicated context-owned generator state");
    Require(Metadata(target.get()).codec_k8r4_terminal_cache_bytes ==
            256U + 70U * 4U + 70U * 8U * 4U,
        "K=8/R=4 AVX2 codec did not cache its maximum-loss terminal");

    struct Neighbor
    {
        unsigned k;
        unsigned r;
        leo2_profile profile;
        uint32_t flags;
    };
    const Neighbor neighbors[] = {
        { 7, 4, LEO2_PROFILE_LEGACY_HIGH_V1, 0 },
        { 8, 3, LEO2_PROFILE_LEGACY_HIGH_V1, 0 },
        { 9, 4, LEO2_PROFILE_LEGACY_HIGH_V1, 0 },
        { 8, 4, LEO2_PROFILE_LOW_V1, 0 },
        { 8, 4, LEO2_PROFILE_LEGACY_HIGH_V1,
          LEO2_CODEC_FORCE_SPECIALIZED_DECODE }
    };
    for (size_t i = 0; i < sizeof(neighbors) / sizeof(neighbors[0]); ++i)
    {
        Codec codec(avx2_context, neighbors[i].k, neighbors[i].r,
            neighbors[i].profile, neighbors[i].flags);
        const leopard2_internal::CodecDecodeMetadataInfo info =
            Metadata(codec.get());
        Require(info.codec_direct_repair_generator_bytes == 0,
            "neighbor codec unexpectedly cached K=8/R=4 repair rows");
        Require(info.codec_k8r4_terminal_cache_bytes == 0,
            "neighbor codec unexpectedly exposed the K=8/R=4 cache");
    }
}

void ExerciseExhaustiveAVX2(
    leo2_context* context,
    uint64_t* decode_count_out,
    uint64_t* oracle_byte_count_out)
{
    Codec codec(context, kOriginalCount, kRecoveryCount);
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh,
            kOriginalCount, kRecoveryCount);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);
    const size_t byte_counts[] = {
        1, 2, 3, 7, 8, 15, 16, 17, 31, 32, 33,
        63, 64, 65, 127, 128, 129, 255, 256, 257,
        511, 512, 513, 1023, 1024, 1025, 2047, 2048
    };
    for (size_t byte_i = 0;
         byte_i < sizeof(byte_counts) / sizeof(byte_counts[0]); ++byte_i)
    {
        Fixture fixture(codec.get(), byte_counts[byte_i],
            UINT64_C(0x4b3852344f4e4553) + byte_i);
        fixture.Encode(codec.get());
        CheckParity(fixture, field, generator);
        *oracle_byte_count_out +=
            fixture.shard_bytes * kRecoveryCount;
        const std::vector<uint8_t> original_before = fixture.original_storage;
        const std::vector<uint8_t> recovery_before = fixture.recovery_storage;
        for (unsigned mask = 0; mask < (1U << kOriginalCount); ++mask)
        {
            if (CountBits(mask) != kRecoveryCount)
                continue;
            fixture.SelectMissing(mask);
            const uint64_t allocations =
                DecodeAndCountAllocations(codec.get(), fixture);
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
            Require(allocations == 0,
                "bounded K=8/R=4 one-shot route allocated");
#else
            (void)allocations;
#endif
            CheckRestored(fixture, mask);
            CheckRestoredGuards(fixture, mask);
            Require(fixture.original_storage == original_before,
                "one-shot route modified an input original or guard");
            Require(fixture.recovery_storage == recovery_before,
                "one-shot route modified parity or guard bytes");
            ++*decode_count_out;
        }
    }
}

void ExerciseFallbackBoundaries(leo2_context* context)
{
    Codec codec(context, kOriginalCount, kRecoveryCount);
    const unsigned full_mask = 0x0f;
    Fixture long_fixture(codec.get(), 2049,
        UINT64_C(0x424f554e44415259));
    long_fixture.Encode(codec.get());
    long_fixture.SelectMissing(full_mask);
    const uint64_t long_allocations =
        DecodeAndCountAllocations(codec.get(), long_fixture);
    CheckRestored(long_fixture, full_mask);
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
    Require(long_allocations != 0,
        "2049-byte neighbor unexpectedly entered bounded one-shot route");
#endif

    Fixture three_loss(codec.get(), 64,
        UINT64_C(0x54485245454c4f53));
    three_loss.Encode(codec.get());
    const unsigned three_mask = 0x0b;
    three_loss.SelectMissing(three_mask);
    const uint64_t three_allocations =
        DecodeAndCountAllocations(codec.get(), three_loss);
    CheckRestored(three_loss, three_mask);
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
    Require(three_allocations != 0,
        "three-loss neighbor unexpectedly entered maximum-loss route");
#endif
}

void ExerciseFailureAtomicity(leo2_context* context)
{
    Codec codec(context, kOriginalCount, kRecoveryCount);
    Fixture fixture(codec.get(), 65, UINT64_C(0x4641494c41544f4d));
    fixture.Encode(codec.get());
    const unsigned mask = 0xa5;
    fixture.SelectMissing(mask);
    const std::vector<uint8_t> output_before = fixture.restored_storage;

    fixture.recovery_present[3] = 0;
    fixture.recovery[3] = NULL;
    RequireResult(leo2_decode(codec.get(), fixture.shard_bytes,
        &fixture.original_present[0], &fixture.recovery_present[0],
        &fixture.original[0], &fixture.recovery[0], &fixture.restored[0],
        fixture.decode_scratch.data(), fixture.decode_scratch_bytes),
        LEO2_NEED_MORE_DATA, "reject insufficient data");
    Require(fixture.restored_storage == output_before,
        "insufficient-data failure modified output");
    fixture.recovery_present[3] = 1;
    fixture.recovery[3] = fixture.RecoveryData(3);

    RequireResult(leo2_decode(codec.get(), fixture.shard_bytes,
        &fixture.original_present[0], &fixture.recovery_present[0],
        &fixture.original[0], &fixture.recovery[0], &fixture.restored[0],
        fixture.decode_scratch.data(), 0),
        LEO2_SCRATCH_TOO_SMALL, "reject short scratch");
    Require(fixture.restored_storage == output_before,
        "short-scratch failure modified output");

    RequireResult(leo2_decode(codec.get(), fixture.shard_bytes,
        &fixture.original_present[0], &fixture.recovery_present[0],
        &fixture.original[0], &fixture.recovery[0], &fixture.restored[0],
        static_cast<uint8_t*>(fixture.decode_scratch.data()) + 1,
        fixture.decode_scratch_bytes), LEO2_BAD_ALIGNMENT,
        "reject misaligned scratch");
    Require(fixture.restored_storage == output_before,
        "misaligned-scratch failure modified output");

    unsigned present_row = 0;
    while (!fixture.original_present[present_row])
        ++present_row;
    fixture.original_present[present_row] = 2;
    RequireResult(leo2_decode(codec.get(), fixture.shard_bytes,
        &fixture.original_present[0], &fixture.recovery_present[0],
        &fixture.original[0], &fixture.recovery[0], &fixture.restored[0],
        fixture.decode_scratch.data(), fixture.decode_scratch_bytes),
        LEO2_INVALID_ARGUMENT, "reject invalid presence marker");
    Require(fixture.restored_storage == output_before,
        "invalid-presence failure modified output");
    fixture.original_present[present_row] = 1;

    unsigned missing_first = 0;
    while (fixture.original_present[missing_first])
        ++missing_first;
    unsigned missing_second = missing_first + 1;
    while (fixture.original_present[missing_second])
        ++missing_second;
    fixture.restored[missing_second] = fixture.restored[missing_first];
    RequireResult(leo2_decode(codec.get(), fixture.shard_bytes,
        &fixture.original_present[0], &fixture.recovery_present[0],
        &fixture.original[0], &fixture.recovery[0], &fixture.restored[0],
        fixture.decode_scratch.data(), fixture.decode_scratch_bytes),
        LEO2_OVERLAP, "reject overlapping restored outputs");
    Require(fixture.restored_storage == output_before,
        "overlap failure modified output");
    fixture.restored[missing_second] = fixture.RestoredData(missing_second);

    fixture.restored[missing_first] = NULL;
    RequireResult(leo2_decode(codec.get(), fixture.shard_bytes,
        &fixture.original_present[0], &fixture.recovery_present[0],
        &fixture.original[0], &fixture.recovery[0], &fixture.restored[0],
        fixture.decode_scratch.data(), fixture.decode_scratch_bytes),
        LEO2_INVALID_ARGUMENT, "reject null restored output");
    Require(fixture.restored_storage == output_before,
        "null-output failure modified output");
}

void ExerciseAllowedAliases(leo2_context* context)
{
    Codec codec(context, kOriginalCount, kRecoveryCount);
    Fixture fixture(codec.get(), 129, UINT64_C(0x414c494153494e50));
    std::memcpy(fixture.OriginalData(1), fixture.OriginalData(0),
        fixture.shard_bytes);
    fixture.original[1] = fixture.original[0];
    fixture.Encode(codec.get());
    const unsigned missing_mask = 0xf0;
    fixture.SelectMissing(missing_mask);
    fixture.original[1] = fixture.original[0];
    const uint64_t allocations = DecodeAndCountAllocations(codec.get(), fixture);
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
    Require(allocations == 0, "aliased-input terminal decode allocated");
#endif
    CheckRestored(fixture, missing_mask);

    Fixture in_place(codec.get(), 33, UINT64_C(0x494e504c41434544));
    in_place.Encode(codec.get());
    const std::vector<uint8_t> expected = in_place.original_storage;
    in_place.SelectMissing(missing_mask, true);
    for (unsigned row = 0; row < kOriginalCount; ++row)
    {
        if ((missing_mask & (1U << row)) != 0)
        {
            std::memset(in_place.OriginalData(row), 0, in_place.shard_bytes);
        }
    }
    const uint64_t in_place_allocations =
        DecodeAndCountAllocations(codec.get(), in_place);
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
    Require(in_place_allocations == 0, "in-place terminal decode allocated");
#endif
    for (unsigned row = 0; row < kOriginalCount; ++row)
    {
        if ((missing_mask & (1U << row)) == 0)
            continue;
        const size_t begin = in_place.original_base + row * in_place.stride;
        Require(std::memcmp(&in_place.original_storage[begin],
                    &expected[begin], in_place.shard_bytes) == 0,
            "in-place restored original differs from source");
    }
}

void ExerciseLowerBackend(leo2_backend backend, const char* name)
{
    Context context(backend);
    if (context.result() == LEO2_UNSUPPORTED)
    {
        std::cout << "SKIP " << name << " backend unavailable\n";
        return;
    }
    RequireResult(context.result(), LEO2_SUCCESS,
        std::string("create ") + name + " context");
    Codec codec(context.get(), kOriginalCount, kRecoveryCount);
    Require(Metadata(codec.get()).codec_direct_repair_generator_bytes == 0,
        std::string(name) + " codec retained AVX2-only row cache");
    Require(Metadata(codec.get()).codec_k8r4_terminal_cache_bytes == 0,
        std::string(name) + " codec retained AVX2-only terminal cache");
    Fixture fixture(codec.get(), 65,
        UINT64_C(0x4c4f574552424143) ^ static_cast<uint64_t>(backend));
    fixture.Encode(codec.get());
    const unsigned mask = 0xa5;
    fixture.SelectMissing(mask);
    (void)DecodeAndCountAllocations(codec.get(), fixture);
    CheckRestored(fixture, mask);
}

} // namespace

int main()
{
    try
    {
        Context avx2(LEO2_BACKEND_AVX2);
        if (avx2.result() == LEO2_UNSUPPORTED)
        {
            std::cout << "K=8/R=4 one-shot direct test skipped: AVX2 unavailable\n";
            return 0;
        }
        RequireResult(avx2.result(), LEO2_SUCCESS, "create AVX2 context");

        CheckCachePredicate(avx2.get());
        uint64_t decode_count = 0;
        uint64_t oracle_byte_count = 0;
        ExerciseExhaustiveAVX2(
            avx2.get(), &decode_count, &oracle_byte_count);
        ExerciseFallbackBoundaries(avx2.get());
        ExerciseFailureAtomicity(avx2.get());
        ExerciseAllowedAliases(avx2.get());
        ExerciseLowerBackend(LEO2_BACKEND_SCALAR, "scalar");
        ExerciseLowerBackend(LEO2_BACKEND_SSSE3, "SSSE3");

        std::cout << "K=8/R=4 one-shot direct checks passed: decodes="
                  << decode_count << " oracle_parity_bytes="
                  << oracle_byte_count
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
                  << " allocation_audit=enabled"
#else
                  << " allocation_audit=disabled"
#endif
                  << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
        g_track_allocations.store(false, std::memory_order_release);
#endif
        std::cerr << "K=8/R=4 one-shot direct failure: "
                  << error.what() << std::endl;
        return 1;
    }
}
