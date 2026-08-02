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
#include "Leopard2Backend.h"
#include "Leopard2Direct.h"
#include "LeopardFF8.h"
#include "leopard.h"
#include "leopard2.h"

#include <algorithm>
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

static const unsigned kOriginalCount = 32;
static const unsigned kRecoveryCount = 32;
static const unsigned kParentCount = 64;
static const size_t kShardBytes = 64;

void Require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void RequireResult(leo2_result result, const char* message)
{
    if (result != LEO2_SUCCESS)
        throw std::runtime_error(std::string(message) + ": " +
            leo2_result_string(result));
}

class AlignedBuffer
{
public:
    explicit AlignedBuffer(size_t bytes)
        : value_(NULL)
        , bytes_(bytes)
    {
#if defined(_MSC_VER)
        value_ = _aligned_malloc(bytes, 64);
#else
        if (posix_memalign(&value_, 64, bytes) != 0)
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

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);
    void* value_;
    size_t bytes_;
};

void FillMessage(std::vector<std::vector<uint8_t> >& message, uint64_t seed)
{
    message.assign(kOriginalCount, std::vector<uint8_t>(kShardBytes));
    uint64_t state = seed;
    for (unsigned row = 0; row < kOriginalCount; ++row)
    {
        for (size_t offset = 0; offset < kShardBytes; ++offset)
        {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            message[row][offset] = static_cast<uint8_t>(
                (state >> 24) ^ row ^ offset);
        }
    }
}

void EncodeIndependentParity(
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator,
    const std::vector<std::vector<uint8_t> >& message,
    std::vector<std::vector<uint8_t> >& parity)
{
    parity.assign(kRecoveryCount, std::vector<uint8_t>(kShardBytes));
    for (unsigned output = 0; output < kRecoveryCount; ++output)
    {
        const std::vector<leopard2_test::Element>& row =
            generator[kOriginalCount + output];
        for (size_t offset = 0; offset < kShardBytes; ++offset)
        {
            leopard2_test::Element value = 0;
            for (unsigned source = 0; source < kOriginalCount; ++source)
            {
                value = field.add(value, field.multiply(
                    row[source], message[source][offset]));
            }
            parity[output][offset] = static_cast<uint8_t>(value);
        }
    }
}

std::vector<std::vector<unsigned> > BuildPatterns()
{
    std::vector<std::vector<unsigned> > patterns;
    std::vector<unsigned> missing;
    for (unsigned i = 0; i < 16; ++i)
        missing.push_back(i);
    patterns.push_back(missing);
    missing.clear();
    for (unsigned i = 16; i < 32; ++i)
        missing.push_back(i);
    patterns.push_back(missing);
    missing.clear();
    for (unsigned i = 0; i < 32; i += 2)
        missing.push_back(i);
    patterns.push_back(missing);
    missing.clear();
    for (unsigned i = 1; i < 32; i += 2)
        missing.push_back(i);
    patterns.push_back(missing);

    uint64_t state = UINT64_C(0x414c473450333236);
    for (unsigned sample = 0; sample < 8; ++sample)
    {
        std::vector<unsigned> order(kOriginalCount);
        for (unsigned i = 0; i < kOriginalCount; ++i)
            order[i] = i;
        for (unsigned i = kOriginalCount - 1; i != 0; --i)
        {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            const unsigned selected = static_cast<unsigned>(state % (i + 1));
            std::swap(order[i], order[selected]);
        }
        order.resize(sample < 4 ? 16 : 31);
        std::sort(order.begin(), order.end());
        patterns.push_back(order);
    }

    // Every possible R-1 survivor location covers the most severe practical
    // pattern while varying every requested-output mask bit.
    for (unsigned survivor = 0; survivor < kOriginalCount; ++survivor)
    {
        missing.clear();
        for (unsigned i = 0; i < kOriginalCount; ++i)
            if (i != survivor)
                missing.push_back(i);
        patterns.push_back(missing);
    }
    return patterns;
}

void SetWorkPointers(AlignedBuffer& storage, void** pointers)
{
    for (unsigned i = 0; i < kParentCount; ++i)
        pointers[i] = storage.bytes() + static_cast<size_t>(i) * kShardBytes;
}

void ExercisePattern(
    const leopard::backend::Ops& avx2,
    const std::vector<std::vector<uint8_t> >& message,
    const std::vector<std::vector<uint8_t> >& parity,
    const std::vector<unsigned>& missing)
{
    uint8_t erased[kParentCount] = {};
    uint8_t requested[kParentCount] = {};
    const void* coordinates[kParentCount] = {};
    uint32_t requested_coordinates[kOriginalCount] = {};
    void* restored[kOriginalCount] = {};
    uint8_t locator[kParentCount] = {};
    uint8_t factor = 0;

    for (unsigned i = 0; i < kOriginalCount; ++i)
        coordinates[i] = &message[i][0];
    for (size_t i = 0; i < missing.size(); ++i)
    {
        const unsigned coordinate = missing[i];
        coordinates[coordinate] = NULL;
        erased[coordinate] = 1;
        requested[coordinate] = 1;
        requested_coordinates[i] = coordinate;
    }
    // The deterministic translated-low selector consumes the lowest L parity
    // rows and marks every surplus row as a virtual erasure.
    for (unsigned i = 0; i < kRecoveryCount; ++i)
    {
        const unsigned coordinate = kOriginalCount + i;
        if (i < missing.size())
            coordinates[coordinate] = &parity[i][0];
        else
            erased[coordinate] = 1;
    }
    leopard::ff8::PrepareDecodeKnownCount(
        kParentCount, erased, kRecoveryCount, locator);
    leopard::ff8::PrepareLowDecode(kParentCount, kOriginalCount, &factor);

    AlignedBuffer candidate_work_storage(kParentCount * kShardBytes);
    AlignedBuffer mature_work_storage(kParentCount * kShardBytes);
    AlignedBuffer generic_work_storage(kParentCount * kShardBytes);
    AlignedBuffer restored_storage(kOriginalCount * kShardBytes);
    void* candidate_work[kParentCount];
    void* mature_work[kParentCount];
    void* generic_work[kParentCount];
    SetWorkPointers(candidate_work_storage, candidate_work);
    SetWorkPointers(mature_work_storage, mature_work);
    SetWorkPointers(generic_work_storage, generic_work);
    std::memset(restored_storage.bytes(), 0xa5,
        kOriginalCount * kShardBytes);
    for (size_t i = 0; i < missing.size(); ++i)
    {
        const unsigned coordinate = missing[i];
        restored[coordinate] = restored_storage.bytes() +
            static_cast<size_t>(coordinate) * kShardBytes;
    }

    Require(leopard::ff8::ReedSolomonDecodeLowP32B64TerminalExperimental(
            coordinates, requested_coordinates,
            static_cast<unsigned>(missing.size()), locator, factor,
            restored, candidate_work),
        "P32/B64 AVX2 terminal rejected a valid pattern");
    leopard::ff8::ReedSolomonDecodeLowPrepared(
        avx2, kShardBytes, kParentCount, kOriginalCount,
        coordinates, requested, locator, &factor, mature_work);
    leopard::ff8::ReedSolomonDecodePrepared(
        avx2, kShardBytes, kParentCount, coordinates,
        requested, locator, generic_work);

    for (unsigned coordinate = 0; coordinate < kOriginalCount; ++coordinate)
    {
        if (requested[coordinate])
        {
            Require(std::memcmp(restored[coordinate],
                    &message[coordinate][0], kShardBytes) == 0,
                "P32/B64 terminal differs from direct systematic data");
            Require(std::memcmp(restored[coordinate],
                    mature_work[coordinate], kShardBytes) == 0,
                "P32/B64 terminal differs from mature Algorithm 4");
            Require(std::memcmp(restored[coordinate],
                    generic_work[coordinate], kShardBytes) == 0,
                "P32/B64 terminal differs from generic decode");
        }
        else
        {
            const uint8_t* untouched = restored_storage.bytes() +
                static_cast<size_t>(coordinate) * kShardBytes;
            for (size_t offset = 0; offset < kShardBytes; ++offset)
                Require(untouched[offset] == 0xa5,
                    "P32/B64 terminal modified an unrequested output");
        }
    }
}

void ExercisePublicRoute(
    leo2_context* context,
    size_t shard_bytes,
    unsigned missing_count,
    bool reusable_plan,
    uint64_t expected_terminal_calls)
{
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context,
        kOriginalCount, kRecoveryCount, LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF8, NULL, &codec), "create public-route codec");

    std::vector<uint8_t> message(kOriginalCount * shard_bytes);
    std::vector<uint8_t> parity(kRecoveryCount * shard_bytes);
    uint64_t state = UINT64_C(0x5055424c49435032) ^ shard_bytes ^
        (static_cast<uint64_t>(missing_count) << 32);
    for (size_t i = 0; i < message.size(); ++i)
    {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        message[i] = static_cast<uint8_t>(state >> 24);
    }
    std::vector<const void*> original(kOriginalCount);
    std::vector<void*> recovery_output(kRecoveryCount);
    for (unsigned i = 0; i < kOriginalCount; ++i)
        original[i] = &message[static_cast<size_t>(i) * shard_bytes];
    for (unsigned i = 0; i < kRecoveryCount; ++i)
        recovery_output[i] = &parity[static_cast<size_t>(i) * shard_bytes];

    size_t encode_scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, shard_bytes, &encode_scratch_bytes),
        "query public-route encode scratch");
    AlignedBuffer encode_scratch(encode_scratch_bytes);
    RequireResult(leo2_encode(codec, shard_bytes, &original[0],
        &recovery_output[0], encode_scratch.bytes(), encode_scratch_bytes),
        "encode public-route parity");

    std::vector<uint8_t> original_present(kOriginalCount, 1);
    std::vector<uint8_t> recovery_present(kRecoveryCount, 1);
    std::vector<const void*> decode_original(original);
    std::vector<const void*> decode_recovery(kRecoveryCount);
    std::vector<void*> restored(kOriginalCount, static_cast<void*>(NULL));
    std::vector<uint8_t> restored_storage(
        static_cast<size_t>(kOriginalCount) * shard_bytes, 0xa5);
    for (unsigned i = 0; i < kRecoveryCount; ++i)
        decode_recovery[i] = recovery_output[i];
    for (unsigned i = 0; i < missing_count; ++i)
    {
        const unsigned coordinate =
            missing_count == 16 ? i * 2U : i;
        original_present[coordinate] = 0;
        decode_original[coordinate] = NULL;
        restored[coordinate] = &restored_storage[
            static_cast<size_t>(coordinate) * shard_bytes];
    }

    leo2_decode_plan* plan = NULL;
    size_t decode_scratch_bytes = 0;
    if (reusable_plan)
    {
        RequireResult(leo2_decode_plan_create(codec,
            &original_present[0], &recovery_present[0], &plan),
            "create reusable public-route plan");
        RequireResult(leo2_decode_plan_scratch_size(
            plan, shard_bytes, &decode_scratch_bytes),
            "query reusable public-route scratch");
    }
    else
    {
        RequireResult(leo2_decode_scratch_size(
            codec, shard_bytes, &decode_scratch_bytes),
            "query one-shot public-route scratch");
    }
    AlignedBuffer decode_scratch(decode_scratch_bytes);
    leo2_test_reset_low_p32_b64_terminal_calls();
    if (reusable_plan)
    {
        RequireResult(leo2_decode_plan_execute(plan, shard_bytes,
            &decode_original[0], &decode_recovery[0], &restored[0],
            decode_scratch.bytes(), decode_scratch_bytes),
            "execute reusable public-route plan");
    }
    else
    {
        RequireResult(leo2_decode(codec, shard_bytes,
            &original_present[0], &recovery_present[0],
            &decode_original[0], &decode_recovery[0], &restored[0],
            decode_scratch.bytes(), decode_scratch_bytes),
            "execute one-shot public route");
    }
    Require(leo2_test_low_p32_b64_terminal_calls() ==
            expected_terminal_calls,
        "public-route terminal call count differs from exact predicate");
    for (unsigned i = 0; i < kOriginalCount; ++i)
    {
        if (original_present[i])
            continue;
        Require(std::memcmp(restored[i],
                &message[static_cast<size_t>(i) * shard_bytes],
                shard_bytes) == 0,
            "public-route terminal recovered incorrect bytes");
    }
    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
}

} // namespace

int main()
{
    try
    {
        Require(leo_init() == Leopard_Success, "Leopard initialization");
        leopard::backend::QualificationStatus status =
            leopard::backend::QualificationAvailable;
        const leopard::backend::Ops* avx2 =
            leopard::backend::GetQualifiedOps(LEO2_BACKEND_AVX2, &status);
        if (!avx2)
        {
            std::printf("SKIP low_p32_b64_terminal: AVX2 unavailable\n");
            return 0;
        }

        const leopard2_test::BinaryField field =
            leopard2_test::make_legacy_gf8();
        const leopard2_test::ProfileLayout layout =
            leopard2_test::make_profile_layout(
                leopard2_test::kLegacyHigh,
                kOriginalCount, kRecoveryCount);
        const leopard2_test::Matrix generator =
            leopard2_test::direct_systematic_generator(field, layout);
        const std::vector<std::vector<unsigned> > patterns = BuildPatterns();
        for (unsigned payload = 0; payload < 2; ++payload)
        {
            std::vector<std::vector<uint8_t> > message;
            std::vector<std::vector<uint8_t> > parity;
            FillMessage(message,
                UINT64_C(0x6c6f775033324236) + payload);
            EncodeIndependentParity(field, generator, message, parity);
            for (size_t i = 0; i < patterns.size(); ++i)
                ExercisePattern(*avx2, message, parity, patterns[i]);
        }
        leo2_context_options options = {};
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AVX2;
        options.thread_count = 1;
        leo2_context* context = NULL;
        RequireResult(leo2_context_create(&options, &context),
            "create public-route AVX2 context");
        ExercisePublicRoute(context, 64, 16, false, 1);
        ExercisePublicRoute(context, 64, 31, false, 1);
        ExercisePublicRoute(context, 63, 16, false, 0);
        ExercisePublicRoute(context, 65, 31, false, 0);
        ExercisePublicRoute(context, 64, 16, true, 0);
        leo2_context_destroy(context);
        std::printf(
            "PASS low_p32_b64_terminal payloads=2 patterns=%zu routes=5\n",
            patterns.size());
        return 0;
    }
    catch (const std::exception& exception)
    {
        std::fprintf(stderr, "FAIL low_p32_b64_terminal: %s\n",
            exception.what());
        return 1;
    }
}
