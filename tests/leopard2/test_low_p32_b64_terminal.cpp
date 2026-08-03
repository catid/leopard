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

void RequireFill(
    const uint8_t* data, size_t bytes, uint8_t expected,
    const char* message)
{
    for (size_t i = 0; i < bytes; ++i)
        Require(data[i] == expected, message);
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

    // Exercise every loss count admitted by the production selector.  The
    // stride is coprime to 32, so each prefix is a distinct, non-contiguous
    // requested-output mask rather than another lowest-index special case.
    for (unsigned loss_count = 9; loss_count < kOriginalCount; ++loss_count)
    {
        missing.clear();
        for (unsigned i = 0; i < loss_count; ++i)
            missing.push_back((3U + i * 13U) & 31U);
        std::sort(missing.begin(), missing.end());
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
    const std::vector<unsigned>& missing,
    bool select_high_parity)
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
    // Cover both ends of the parity block.  Selecting only the upper L rows
    // makes the unselected lower rows genuine/virtual parity erasures and
    // exercises locator normalization independently of the usual prefix.
    for (unsigned i = 0; i < kRecoveryCount; ++i)
    {
        const unsigned coordinate = kOriginalCount + i;
        const bool selected = select_high_parity
            ? i >= kRecoveryCount - missing.size()
            : i < missing.size();
        if (selected)
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
            avx2, coordinates, requested_coordinates,
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
    unsigned original_count,
    unsigned recovery_count,
    size_t shard_bytes,
    unsigned missing_count,
    bool reusable_plan,
    bool select_high_parity,
    uint64_t expected_p32_terminal_calls,
    uint64_t expected_p128_terminal_calls = 0,
    unsigned missing_stride = 1,
    unsigned missing_offset = 0,
    bool unaligned = false,
    unsigned parity_stride = 0,
    unsigned parity_offset = 0,
    bool alias_surviving_inputs = false,
    bool batch_binding = false)
{
    Require(!batch_binding || reusable_plan,
        "decode binding requires a reusable public-route plan");
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context,
        original_count, recovery_count, LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF8, NULL, &codec), "create public-route codec");

    Require(missing_stride != 0,
        "public-route missing stride must be nonzero");
    const size_t row_stride = shard_bytes + (unaligned ? 7 : 0);
    const size_t message_base = unaligned ? 1 : 0;
    const size_t parity_base = unaligned ? 3 : 0;
    const size_t restored_base = unaligned ? 5 : 0;
    std::vector<uint8_t> message(
        message_base + original_count * row_stride + 8, 0xcc);
    std::vector<uint8_t> parity(
        parity_base + recovery_count * row_stride + 8, 0xdd);
    uint64_t state = UINT64_C(0x5055424c49435032) ^ shard_bytes ^
        (static_cast<uint64_t>(missing_count) << 32) ^
        (static_cast<uint64_t>(original_count) << 16) ^ recovery_count ^
        (select_high_parity ? UINT64_C(0x8000000000000000) : 0);
    for (unsigned row = 0; row < original_count; ++row)
    {
        uint8_t* output = &message[
            message_base + static_cast<size_t>(row) * row_stride];
        for (size_t offset = 0; offset < shard_bytes; ++offset)
        {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            output[offset] = static_cast<uint8_t>(state >> 24);
        }
    }
    if (alias_surviving_inputs)
    {
        Require(original_count >= 2,
            "public-route input-alias case needs two originals");
        std::memcpy(&message[message_base +
                static_cast<size_t>(original_count - 1) * row_stride],
            &message[message_base +
                static_cast<size_t>(original_count - 2) * row_stride],
            shard_bytes);
    }
    const std::vector<uint8_t> original_message = message;
    std::vector<const void*> original(original_count);
    std::vector<void*> recovery_output(recovery_count);
    for (unsigned i = 0; i < original_count; ++i)
        original[i] = &message[
            message_base + static_cast<size_t>(i) * row_stride];
    for (unsigned i = 0; i < recovery_count; ++i)
        recovery_output[i] = &parity[
            parity_base + static_cast<size_t>(i) * row_stride];

    size_t encode_scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, shard_bytes, &encode_scratch_bytes),
        "query public-route encode scratch");
    AlignedBuffer encode_scratch(encode_scratch_bytes + 128);
    std::memset(encode_scratch.bytes(), 0xe1, encode_scratch_bytes + 128);
    RequireResult(leo2_encode(codec, shard_bytes, &original[0],
        &recovery_output[0], encode_scratch.bytes() + 64,
        encode_scratch_bytes),
        "encode public-route parity");
    RequireFill(encode_scratch.bytes(), 64, 0xe1,
        "public-route encode scratched before its range");
    RequireFill(encode_scratch.bytes() + 64 + encode_scratch_bytes,
        64, 0xe1, "public-route encode scratched after its range");
    Require(message == original_message,
        "public-route encoder modified systematic input or its guards");
    const std::vector<uint8_t> encoded_parity = parity;

    std::vector<uint8_t> original_present(original_count, 1);
    std::vector<uint8_t> recovery_present(recovery_count, 1);
    std::vector<const void*> decode_original(original);
    std::vector<const void*> decode_recovery(recovery_count);
    std::vector<void*> restored(original_count, static_cast<void*>(NULL));
    std::vector<uint8_t> restored_storage(
        restored_base + static_cast<size_t>(original_count) * row_stride + 8,
        0xa5);
    for (unsigned i = 0; i < recovery_count; ++i)
        decode_recovery[i] = recovery_output[i];
    for (unsigned i = 0; i < original_count; ++i)
        restored[i] = &restored_storage[restored_base +
            static_cast<size_t>(i) * row_stride];
    for (unsigned i = 0; i < missing_count; ++i)
    {
        const unsigned coordinate = static_cast<unsigned>(
            (missing_offset + static_cast<uint64_t>(i) * missing_stride) %
            original_count);
        Require(original_present[coordinate] != 0,
            "public-route missing stride generated a duplicate coordinate");
        original_present[coordinate] = 0;
        decode_original[coordinate] = NULL;
    }
    if (alias_surviving_inputs)
    {
        const unsigned first = original_count - 2;
        const unsigned second = original_count - 1;
        Require(original_present[first] != 0 && original_present[second] != 0,
            "public-route input-alias rows must both survive");
        Require(std::memcmp(decode_original[first], decode_original[second],
                shard_bytes) == 0,
            "public-route aliased rows must contain identical bytes");
        decode_original[second] = decode_original[first];
    }
    if (parity_stride != 0)
    {
        std::fill(recovery_present.begin(), recovery_present.end(), 0);
        std::fill(decode_recovery.begin(), decode_recovery.end(),
            static_cast<const void*>(NULL));
        for (unsigned i = 0; i < missing_count; ++i)
        {
            const unsigned coordinate = static_cast<unsigned>(
                (parity_offset + static_cast<uint64_t>(i) * parity_stride) %
                recovery_count);
            Require(recovery_present[coordinate] == 0,
                "public-route parity stride generated a duplicate coordinate");
            recovery_present[coordinate] = 1;
            decode_recovery[coordinate] = recovery_output[coordinate];
        }
    }
    else if (select_high_parity)
    {
        Require(missing_count <= recovery_count,
            "high-parity route does not have enough recovery rows");
        for (unsigned i = 0; i < recovery_count - missing_count; ++i)
        {
            recovery_present[i] = 0;
            decode_recovery[i] = NULL;
        }
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
    AlignedBuffer decode_scratch(decode_scratch_bytes + 128);
    std::memset(decode_scratch.bytes(), 0xe2, decode_scratch_bytes + 128);
    leo2_test_reset_low_p32_b64_terminal_calls();
    leo2_test_reset_low_p128_b64_terminal_calls();
    if (batch_binding)
    {
        std::vector<uint8_t> second_restored_storage(
            restored_storage.size(), 0xa5);
        std::vector<void*> second_restored(
            original_count, static_cast<void*>(NULL));
        for (unsigned i = 0; i < original_count; ++i)
        {
            second_restored[i] = &second_restored_storage[restored_base +
                static_cast<size_t>(i) * row_stride];
        }
        AlignedBuffer second_scratch(decode_scratch_bytes + 128);
        std::memset(
            second_scratch.bytes(), 0xe2, decode_scratch_bytes + 128);
        leo2_decode_batch_item items[2] = {};
        for (unsigned item = 0; item < 2; ++item)
        {
            items[item].shard_bytes = shard_bytes;
            items[item].original = &decode_original[0];
            items[item].recovery = &decode_recovery[0];
            items[item].restored_original = item == 0
                ? &restored[0] : &second_restored[0];
            items[item].scratch = item == 0
                ? decode_scratch.bytes() + 64
                : second_scratch.bytes() + 64;
            items[item].scratch_bytes = decode_scratch_bytes;
        }
        leo2_decode_batch_binding* binding = NULL;
        RequireResult(leo2_decode_batch_binding_create(
            plan, items, 2, &binding),
            "create reusable public-route binding");
        Require(leo2_decode_batch_binding_item_count(binding) == 2,
            "public-route binding item count differs");
        RequireResult(leo2_decode_batch_binding_execute(binding),
            "execute reusable public-route binding");
        leo2_decode_batch_binding_destroy(binding);
        RequireFill(second_scratch.bytes(), 64, 0xe2,
            "binding decode scratched before its second range");
        RequireFill(second_scratch.bytes() + 64 + decode_scratch_bytes,
            64, 0xe2, "binding decode scratched after its second range");
        std::vector<uint8_t> expected_second(
            second_restored_storage.size(), 0xa5);
        for (unsigned i = 0; i < original_count; ++i)
        {
            if (original_present[i])
                continue;
            Require(std::memcmp(second_restored[i], original[i],
                    shard_bytes) == 0,
                "binding terminal recovered incorrect second-item bytes");
            std::memcpy(&expected_second[restored_base +
                    static_cast<size_t>(i) * row_stride],
                original[i], shard_bytes);
        }
        Require(second_restored_storage == expected_second,
            "binding terminal wrote outside second-item output ranges");
    }
    else if (reusable_plan)
    {
        RequireResult(leo2_decode_plan_execute(plan, shard_bytes,
            &decode_original[0], &decode_recovery[0], &restored[0],
            decode_scratch.bytes() + 64, decode_scratch_bytes),
            "execute reusable public-route plan");
    }
    else
    {
        RequireResult(leo2_decode(codec, shard_bytes,
            &original_present[0], &recovery_present[0],
            &decode_original[0], &decode_recovery[0], &restored[0],
            decode_scratch.bytes() + 64, decode_scratch_bytes),
            "execute one-shot public route");
    }
    RequireFill(decode_scratch.bytes(), 64, 0xe2,
        "public-route decode scratched before its range");
    RequireFill(decode_scratch.bytes() + 64 + decode_scratch_bytes,
        64, 0xe2, "public-route decode scratched after its range");
    const uint64_t actual_p32_calls =
        leo2_test_low_p32_b64_terminal_calls();
    const uint64_t actual_p128_calls =
        leo2_test_low_p128_b64_terminal_calls();
    if (actual_p32_calls != expected_p32_terminal_calls ||
        actual_p128_calls != expected_p128_terminal_calls)
    {
        std::fprintf(stderr,
            "route mismatch K=%u R=%u B=%zu L=%u reusable=%u "
            "binding=%u P32=%llu/%llu P128=%llu/%llu\n",
            original_count, recovery_count, shard_bytes, missing_count,
            reusable_plan ? 1U : 0U, batch_binding ? 1U : 0U,
            static_cast<unsigned long long>(actual_p32_calls),
            static_cast<unsigned long long>(expected_p32_terminal_calls),
            static_cast<unsigned long long>(actual_p128_calls),
            static_cast<unsigned long long>(expected_p128_terminal_calls));
    }
    Require(actual_p32_calls == expected_p32_terminal_calls,
        "public-route P32 terminal call count differs from exact predicate");
    Require(actual_p128_calls == expected_p128_terminal_calls,
        "public-route P128 terminal call count differs from exact predicate");
    std::vector<uint8_t> expected_restored(restored_storage.size(), 0xa5);
    for (unsigned i = 0; i < original_count; ++i)
    {
        if (original_present[i])
            continue;
        Require(std::memcmp(restored[i],
                original[i],
                shard_bytes) == 0,
            "public-route terminal recovered incorrect bytes");
        std::memcpy(&expected_restored[restored_base +
                static_cast<size_t>(i) * row_stride],
            original[i], shard_bytes);
    }
    Require(restored_storage == expected_restored,
        "public-route terminal wrote outside requested output ranges");
    Require(message == original_message,
        "public-route decode modified systematic input or its guards");
    Require(parity == encoded_parity,
        "public-route decode modified recovery input or its guards");
    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
}

void EncodeIndependentParityForCount(
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator,
    const std::vector<std::vector<uint8_t> >& message,
    std::vector<std::vector<uint8_t> >& parity,
    unsigned count)
{
    parity.assign(count, std::vector<uint8_t>(kShardBytes));
    for (unsigned output = 0; output < count; ++output)
    {
        const std::vector<leopard2_test::Element>& row =
            generator[count + output];
        for (size_t offset = 0; offset < kShardBytes; ++offset)
        {
            leopard2_test::Element value = 0;
            for (unsigned source = 0; source < count; ++source)
            {
                value = field.add(value, field.multiply(
                    row[source], message[source][offset]));
            }
            parity[output][offset] = static_cast<uint8_t>(value);
        }
    }
}

std::vector<std::vector<unsigned> > BuildP128Patterns(unsigned count)
{
    Require(count == 95 || count == 128,
        "P128 pattern count is outside the qualified shapes");
    std::vector<std::vector<unsigned> > patterns;
    std::vector<unsigned> missing;
    const unsigned half_loss = count == 95 ? 47 : 64;
    const unsigned stride = count == 95 ? 37 : 29;
    for (unsigned sample = 0; sample < 8; ++sample)
    {
        missing.clear();
        for (unsigned i = 0; i < half_loss; ++i)
            missing.push_back((sample + i * stride) % count);
        std::sort(missing.begin(), missing.end());
        patterns.push_back(missing);
    }
    // Every possible single-survivor position stresses the L=K-1 locator,
    // weighted live-mask lane, requested-mask boundary, and shortened tail.
    for (unsigned survivor = 0; survivor < count; ++survivor)
    {
        missing.clear();
        for (unsigned i = 0; i < count; ++i)
            if (i != survivor)
                missing.push_back(i);
        patterns.push_back(missing);
    }
    return patterns;
}

void ExerciseP128Pattern(
    const leopard::backend::Ops& avx2,
    const std::vector<std::vector<uint8_t> >& message,
    const std::vector<std::vector<uint8_t> >& parity,
    const std::vector<unsigned>& missing,
    bool select_high_parity,
    unsigned count)
{
    static const unsigned kPadded = 128;
    static const unsigned kParent = 256;
    uint8_t erased[kParent] = {};
    uint8_t requested[kParent] = {};
    const void* coordinates[kParent] = {};
    uint32_t requested_coordinates[kPadded] = {};
    void* restored[kPadded] = {};
    uint8_t locator[kParent] = {};
    uint8_t factor = 0;

    for (unsigned i = 0; i < count; ++i)
        coordinates[i] = &message[i][0];
    for (size_t i = 0; i < missing.size(); ++i)
    {
        const unsigned coordinate = missing[i];
        coordinates[coordinate] = NULL;
        erased[coordinate] = 1;
        requested[coordinate] = 1;
        requested_coordinates[i] = coordinate;
    }
    for (unsigned i = 0; i < count; ++i)
    {
        const unsigned coordinate = kPadded + i;
        const bool selected = select_high_parity
            ? i >= count - missing.size()
            : i < missing.size();
        if (selected)
            coordinates[coordinate] = &parity[i][0];
        else
            erased[coordinate] = 1;
    }
    for (unsigned coordinate = kPadded + count;
         coordinate < kParent; ++coordinate)
        erased[coordinate] = 1;
    leopard::ff8::PrepareDecodeKnownCount(
        kParent, erased, kPadded, locator);
    leopard::ff8::PrepareLowDecode(kParent, kPadded, &factor);

    AlignedBuffer candidate_work_storage(kParent * kShardBytes + 128);
    AlignedBuffer mature_work_storage(kParent * kShardBytes);
    AlignedBuffer generic_work_storage(kParent * kShardBytes);
    AlignedBuffer restored_storage(kPadded * kShardBytes + 128);
    void* candidate_work[kParent];
    void* mature_work[kParent];
    void* generic_work[kParent];
    for (unsigned i = 0; i < kParent; ++i)
    {
        candidate_work[i] = candidate_work_storage.bytes() + 64 +
            static_cast<size_t>(i) * kShardBytes;
        mature_work[i] = mature_work_storage.bytes() +
            static_cast<size_t>(i) * kShardBytes;
        generic_work[i] = generic_work_storage.bytes() +
            static_cast<size_t>(i) * kShardBytes;
    }
    std::memset(candidate_work_storage.bytes(), 0xd1,
        kParent * kShardBytes + 128);
    std::memset(restored_storage.bytes(), 0xa5,
        kPadded * kShardBytes + 128);
    for (unsigned coordinate = 0; coordinate < kPadded; ++coordinate)
        restored[coordinate] = restored_storage.bytes() + 64 +
            static_cast<size_t>(coordinate) * kShardBytes;
    for (size_t i = 0; i < missing.size(); ++i)
    {
        const unsigned coordinate = missing[i];
        Require(coordinate < count,
            "P128 pattern requested a shortened coordinate");
    }

    Require(leopard::ff8::ReedSolomonDecodeLowP128B64TerminalExperimental(
            avx2, coordinates, requested_coordinates,
            static_cast<unsigned>(missing.size()), locator, factor,
            restored, candidate_work),
        "P128/B64 AVX2 terminal rejected a valid pattern");
    leopard::ff8::ReedSolomonDecodeLowPrepared(
        avx2, kShardBytes, kParent, kPadded,
        coordinates, requested, locator, &factor, mature_work);
    leopard::ff8::ReedSolomonDecodePrepared(
        avx2, kShardBytes, kParent, coordinates,
        requested, locator, generic_work);

    RequireFill(candidate_work_storage.bytes(), 64, 0xd1,
        "P128/B64 terminal scratched before its work range");
    RequireFill(candidate_work_storage.bytes() + 64 +
            kParent * kShardBytes,
        64, 0xd1, "P128/B64 terminal scratched after its work range");
    RequireFill(restored_storage.bytes(), 64, 0xa5,
        "P128/B64 terminal scratched before its output range");
    RequireFill(restored_storage.bytes() + 64 +
            kPadded * kShardBytes,
        64, 0xa5, "P128/B64 terminal scratched after its output range");

    for (unsigned coordinate = 0; coordinate < kPadded; ++coordinate)
    {
        if (requested[coordinate])
        {
            Require(coordinate < count,
                "P128 terminal requested a shortened coordinate");
            Require(std::memcmp(restored[coordinate],
                    &message[coordinate][0], kShardBytes) == 0,
                "P128/B64 terminal differs from direct systematic data");
            Require(std::memcmp(restored[coordinate],
                    mature_work[coordinate], kShardBytes) == 0,
                "P128/B64 terminal differs from mature Algorithm 4");
            Require(std::memcmp(restored[coordinate],
                    generic_work[coordinate], kShardBytes) == 0,
                "P128/B64 terminal differs from generic decode");
        }
        else
        {
            const uint8_t* untouched = restored_storage.bytes() + 64 +
                static_cast<size_t>(coordinate) * kShardBytes;
            for (size_t offset = 0; offset < kShardBytes; ++offset)
                Require(untouched[offset] == 0xa5,
                    "P128/B64 terminal modified an unrequested output");
        }
    }
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
            {
                ExercisePattern(
                    *avx2, message, parity, patterns[i], false);
                ExercisePattern(
                    *avx2, message, parity, patterns[i], true);
            }
        }

        const leopard2_test::ProfileLayout layout128 =
            leopard2_test::make_profile_layout(
                leopard2_test::kLegacyHigh, 128, 128);
        const leopard2_test::Matrix generator128 =
            leopard2_test::direct_systematic_generator(field, layout128);
        const std::vector<std::vector<unsigned> > patterns128 =
            BuildP128Patterns(128);
        for (unsigned payload = 0; payload < 2; ++payload)
        {
            std::vector<std::vector<uint8_t> > message128(
                128, std::vector<uint8_t>(kShardBytes));
            std::vector<std::vector<uint8_t> > parity128;
            uint64_t state = UINT64_C(0x6c6f775031323836) + payload;
            for (unsigned row = 0; row < 128; ++row)
            {
                for (size_t offset = 0; offset < kShardBytes; ++offset)
                {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    message128[row][offset] = static_cast<uint8_t>(
                        (state >> 24) ^ row ^ offset);
                }
            }
            EncodeIndependentParityForCount(
                field, generator128, message128, parity128, 128);
            for (size_t i = 0; i < patterns128.size(); ++i)
            {
                ExerciseP128Pattern(
                    *avx2, message128, parity128, patterns128[i], false, 128);
                ExerciseP128Pattern(
                    *avx2, message128, parity128, patterns128[i], true, 128);
            }
        }

        // Shortened/punctured K=R=95 receives independently generated
        // legacy-high parity, then feeds the exact translated P=128 layout to
        // the candidate, mature Algorithm 4, and generic decoder.
        const leopard2_test::ProfileLayout layout95 =
            leopard2_test::make_profile_layout(
                leopard2_test::kLegacyHigh, 95, 95);
        const leopard2_test::Matrix generator95 =
            leopard2_test::direct_systematic_generator(field, layout95);
        const std::vector<std::vector<unsigned> > patterns95 =
            BuildP128Patterns(95);
        std::vector<std::vector<uint8_t> > message95(
            95, std::vector<uint8_t>(kShardBytes));
        std::vector<std::vector<uint8_t> > parity95;
        uint64_t state95 = UINT64_C(0x6c6f775039353634);
        for (unsigned row = 0; row < 95; ++row)
        {
            for (size_t offset = 0; offset < kShardBytes; ++offset)
            {
                state95 ^= state95 << 13;
                state95 ^= state95 >> 7;
                state95 ^= state95 << 17;
                message95[row][offset] = static_cast<uint8_t>(
                    (state95 >> 24) ^ row ^ offset);
            }
        }
        EncodeIndependentParityForCount(
            field, generator95, message95, parity95, 95);
        for (size_t i = 0; i < patterns95.size(); ++i)
        {
            ExerciseP128Pattern(
                *avx2, message95, parity95, patterns95[i], false, 95);
            ExerciseP128Pattern(
                *avx2, message95, parity95, patterns95[i], true, 95);
        }
        leo2_context_options options = {};
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AVX2;
        // Two-item bindings exercise one immutable plan concurrently through
        // the persistent worker pool; ordinary plan calls remain single-item.
        options.thread_count = 2;
        leo2_context* context = NULL;
        RequireResult(leo2_context_create(&options, &context),
            "create public-route AVX2 context");
        Require(leopard2_internal::LowP32B64TerminalModeForDiagnostics() == 1,
            "production terminal mode did not start enabled");
        Require(leopard2_internal::LowP128B64TerminalModeForDiagnostics() == 1,
            "qualified P128 terminal did not start enabled");
        unsigned route_count = 0;
        for (unsigned loss_count = 9; loss_count < kOriginalCount;
             ++loss_count)
        {
            ExercisePublicRoute(context, 32, 32, 64, loss_count,
                false, false, 1);
            ++route_count;
            ExercisePublicRoute(context, 32, 32, 64, loss_count,
                false, true, 1);
            ++route_count;
        }
        // Explicit route-off neighbors prove every selector boundary: direct
        // repair below L=9, all-loss at L=32, ragged byte tails, and adjacent
        // public K/R shapes.  The exact reusable plan now shares the promoted
        // terminal with the one-shot route.
        ExercisePublicRoute(context, 32, 32, 64, 8, false, false, 0);
        ExercisePublicRoute(context, 32, 32, 64, 32, false, false, 0);
        ExercisePublicRoute(context, 32, 32, 63, 16, false, false, 0);
        ExercisePublicRoute(context, 32, 32, 65, 16, false, false, 0);
        ExercisePublicRoute(context, 32, 32, 64, 16, true, false, 1);
        ExercisePublicRoute(context, 32, 32, 63, 16, true, false, 0);
        ExercisePublicRoute(context, 32, 32, 65, 16, true, false, 0);
        ExercisePublicRoute(context, 31, 32, 64, 16, false, false, 0);
        ExercisePublicRoute(context, 32, 31, 64, 16, false, false, 0);
        route_count += 9;
        ExercisePublicRoute(context, 32, 32, 64, 16,
            true, false, 2, 0, 1, 0, false, 0, 0, false, true);
        ++route_count;

        Require(leopard2_internal::SetLowP32B64TerminalEnabledForDiagnostics(
                true),
            "arm enabled terminal route probe");
        ExercisePublicRoute(
            context, 32, 32, 64, 16, false, true, 1);
        ++route_count;
        Require(leopard2_internal::
                LowP32B64TerminalRouteSelectedForDiagnostics(),
            "enabled route probe did not observe the exact terminal");
        Require(leopard2_internal::
                FinishLowP32B64TerminalRouteProbeForDiagnostics(),
            "finish enabled terminal route probe");

        Require(leopard2_internal::SetLowP32B64TerminalEnabledForDiagnostics(
                true),
            "arm route-off loss-neighbor probe");
        ExercisePublicRoute(context, 32, 32, 64, 8, false, false, 0);
        ++route_count;
        Require(!leopard2_internal::
                LowP32B64TerminalRouteSelectedForDiagnostics(),
            "L=8 direct-repair neighbor reached the exact terminal");
        Require(leopard2_internal::
                FinishLowP32B64TerminalRouteProbeForDiagnostics(),
            "finish route-off loss-neighbor probe");
        Require(leopard2_internal::SetLowP32B64TerminalEnabledForDiagnostics(
                false),
            "disable terminal for same-executable control");
        Require(leopard2_internal::LowP32B64TerminalModeForDiagnostics() == 2,
            "disabled terminal did not select mode word two");
        ExercisePublicRoute(
            context, 32, 32, 64, 31, false, false, 0);
        ++route_count;
        Require(!leopard2_internal::
                LowP32B64TerminalRouteSelectedForDiagnostics(),
            "disabled route probe observed the exact terminal");
        Require(leopard2_internal::
                FinishLowP32B64TerminalRouteProbeForDiagnostics(),
            "finish disabled terminal route probe");
        Require(leopard2_internal::SetLowP32B64TerminalEnabledForDiagnostics(
                true),
            "restore terminal after same-executable control");
        Require(leopard2_internal::LowP32B64TerminalModeForDiagnostics() == 1,
            "enabled terminal did not restore mode word one");
        Require(leopard2_internal::
                FinishLowP32B64TerminalRouteProbeForDiagnostics(),
            "finish final enabled terminal route probe");

        Require(leopard2_internal::SetLowP128B64TerminalEnabledForDiagnostics(
                true),
            "arm enabled P128 terminal route probe");
        Require(leopard2_internal::LowP128B64TerminalModeForDiagnostics() == 1,
            "armed enabled P128 terminal did not normalize to mode one");
        const unsigned p128_shapes[][2] = {
            { 95, 47 }, { 95, 94 }, { 128, 64 }, { 128, 127 }
        };
        for (size_t shape = 0;
             shape < sizeof(p128_shapes) / sizeof(p128_shapes[0]); ++shape)
        {
            for (unsigned high = 0; high < 2; ++high)
            {
                ExercisePublicRoute(context,
                    p128_shapes[shape][0], p128_shapes[shape][0], 64,
                    p128_shapes[shape][1], false, high != 0, 0, 1);
                ++route_count;
            }
            ExercisePublicRoute(context,
                p128_shapes[shape][0], p128_shapes[shape][0], 64,
                p128_shapes[shape][1], true, false, 0, 1);
            ++route_count;
        }
        ExercisePublicRoute(context, 128, 128, 64, 127,
            true, false, 0, 2, 1, 0, false, 0, 0, false, true);
        ++route_count;
        // Non-contiguous K=95 patterns exercise shortening/puncturing and
        // deliberately unaligned live inputs/outputs through both parity
        // selection ends.  Stride 37 is coprime to 95.
        ExercisePublicRoute(context, 95, 95, 64, 47,
            false, false, 0, 1, 37, 11, true, 41, 13);
        ExercisePublicRoute(context, 95, 95, 64, 94,
            false, true, 0, 1, 37, 7, true, 41, 3);
        // K=128 uses an odd stride to cover both requested-mask words while
        // retaining unique coordinates; all AVX loads/stores must tolerate
        // the distinct source/parity/output misalignments above.
        ExercisePublicRoute(context, 128, 128, 64, 64,
            false, false, 0, 1, 29, 9, true);
        ExercisePublicRoute(context, 128, 128, 64, 127,
            false, true, 0, 1, 29, 5, true);
        // Input/input aliases are permitted when the bytes agree.  The last
        // two systematic rows survive this prefix-loss case and deliberately
        // share one decode pointer after their equal data was encoded.
        ExercisePublicRoute(context, 128, 128, 64, 64,
            false, false, 0, 1, 1, 0, false, 0, 0, true);
        route_count += 5;
        Require(leopard2_internal::
                LowP128B64TerminalRouteSelectedForDiagnostics(),
            "enabled P128 route probe did not observe the exact terminal");
        Require(leopard2_internal::
                FinishLowP128B64TerminalRouteProbeForDiagnostics(),
            "finish enabled P128 terminal route probe");
        Require(leopard2_internal::LowP128B64TerminalModeForDiagnostics() == 1,
            "finished enabled P128 terminal did not retain mode one");
        ExercisePublicRoute(context, 128, 128, 64, 64,
            false, false, 0, 1);
        ++route_count;

        Require(leopard2_internal::SetLowP128B64TerminalEnabledForDiagnostics(
                true),
            "arm off-target P128 route probe");
        ExercisePublicRoute(context, 128, 128, 64, 65,
            false, false, 0, 0);
        ++route_count;
        Require(!leopard2_internal::
                LowP128B64TerminalRouteSelectedForDiagnostics(),
            "off-target P128 route probe observed the exact terminal");
        Require(leopard2_internal::
                FinishLowP128B64TerminalRouteProbeForDiagnostics(),
            "finish off-target P128 route probe");

        // Every public selector edge remains on the mature path.
        ExercisePublicRoute(context, 95, 95, 64, 46,
            false, false, 0, 0);
        ExercisePublicRoute(context, 95, 95, 64, 48,
            false, false, 0, 0);
        ExercisePublicRoute(context, 95, 95, 64, 93,
            false, false, 0, 0);
        ExercisePublicRoute(context, 95, 95, 64, 95,
            false, false, 0, 0);
        ExercisePublicRoute(context, 128, 128, 64, 63,
            false, false, 0, 0);
        ExercisePublicRoute(context, 128, 128, 64, 65,
            false, false, 0, 0);
        ExercisePublicRoute(context, 128, 128, 64, 126,
            false, false, 0, 0);
        ExercisePublicRoute(context, 128, 128, 64, 128,
            false, false, 0, 0);
        ExercisePublicRoute(context, 128, 128, 63, 127,
            false, false, 0, 0);
        ExercisePublicRoute(context, 128, 128, 65, 127,
            false, false, 0, 0);
        ExercisePublicRoute(context, 128, 128, 63, 127,
            true, false, 0, 0);
        ExercisePublicRoute(context, 128, 128, 65, 127,
            true, false, 0, 0);
        ExercisePublicRoute(context, 127, 128, 64, 127,
            false, false, 0, 0);
        ExercisePublicRoute(context, 128, 127, 64, 126,
            false, false, 0, 0);
        ExercisePublicRoute(context, 94, 95, 64, 47,
            false, false, 0, 0);
        ExercisePublicRoute(context, 95, 94, 64, 47,
            false, false, 0, 0);
        ExercisePublicRoute(context, 127, 128, 64, 64,
            false, false, 0, 0);
        ExercisePublicRoute(context, 128, 127, 64, 64,
            false, false, 0, 0);
        route_count += 18;

        Require(leopard2_internal::SetLowP128B64TerminalEnabledForDiagnostics(
                false),
            "disable P128 terminal for same-executable control");
        Require(leopard2_internal::LowP128B64TerminalModeForDiagnostics() == 2,
            "disabled P128 terminal did not normalize to mode two");
        ExercisePublicRoute(context, 128, 128, 64, 127,
            false, false, 0, 0);
        ++route_count;
        Require(!leopard2_internal::
                LowP128B64TerminalRouteSelectedForDiagnostics(),
            "disabled P128 route probe observed the exact terminal");
        Require(leopard2_internal::
                FinishLowP128B64TerminalRouteProbeForDiagnostics(),
            "finish disabled P128 terminal route probe");
        Require(leopard2_internal::LowP128B64TerminalModeForDiagnostics() == 2,
            "finished disabled P128 terminal did not retain mode two");
        Require(leopard2_internal::SetLowP128B64TerminalEnabledForDiagnostics(
                true),
            "restore qualified P128 terminal after disabled control");
        Require(leopard2_internal::LowP128B64TerminalModeForDiagnostics() == 1,
            "restored P128 terminal did not normalize to mode one");
        Require(leopard2_internal::
                FinishLowP128B64TerminalRouteProbeForDiagnostics(),
            "finish restored P128 terminal route probe");
        leo2_context_destroy(context);
        std::printf(
            "PASS low_p32_p128_b64_terminal p32_payloads=2 "
            "p32_patterns=%zu p128_payloads=2 p128_patterns=%zu "
            "p95_payloads=1 p95_patterns=%zu "
            "parity_selections=2 routes=%u\n",
            patterns.size(), patterns128.size(), patterns95.size(),
            route_count);
        return 0;
    }
    catch (const std::exception& exception)
    {
        std::fprintf(stderr, "FAIL low_p32_b64_terminal: %s\n",
            exception.what());
        return 1;
    }
}
