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
#include "leopard.h"
#include "leopard2.h"

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

namespace {

using leopard2_test::BinaryField;
using leopard2_test::Element;
using leopard2_test::Matrix;
using leopard2_test::ProfileLayout;

static const uint8_t kGuardByte = 0xd3;

struct TestCounts
{
    uint64_t high_cases;
    uint64_t high_bytes_compared;
    uint64_t low_cases;
    uint64_t low_symbols_compared;
    uint64_t randomized_patterns;
    uint64_t recovered_shards;
    uint64_t repeated_executions;
    uint64_t generic_comparisons;
    uint64_t subset_cases;
    uint64_t batch_stripes;
    uint64_t no_op_cases;
    uint64_t invalid_checks;
    uint64_t guard_checks;

    TestCounts()
        : high_cases(0)
        , high_bytes_compared(0)
        , low_cases(0)
        , low_symbols_compared(0)
        , randomized_patterns(0)
        , recovered_shards(0)
        , repeated_executions(0)
        , generic_comparisons(0)
        , subset_cases(0)
        , batch_stripes(0)
        , no_op_cases(0)
        , invalid_checks(0)
        , guard_checks(0)
    {}
};

class Random
{
public:
    explicit Random(uint64_t seed)
        : state_(seed ? seed : UINT64_C(0x9e3779b97f4a7c15))
    {}

    uint64_t next()
    {
        uint64_t value = state_;
        value ^= value >> 12;
        value ^= value << 25;
        value ^= value >> 27;
        state_ = value;
        return value * UINT64_C(2685821657736338717);
    }

    unsigned below(unsigned limit)
    {
        return limit == 0 ? 0 : static_cast<unsigned>(next() % limit);
    }

private:
    uint64_t state_;
};

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void require_result(leo2_result actual, leo2_result expected, const std::string& operation)
{
    if (actual == expected)
        return;
    std::ostringstream stream;
    stream << operation << ": expected " << leo2_result_string(expected)
           << " (" << static_cast<int>(expected) << "), received "
           << leo2_result_string(actual) << " (" << static_cast<int>(actual) << ')';
    throw std::runtime_error(stream.str());
}

void require_success(leo2_result result, const std::string& operation)
{
    require_result(result, LEO2_SUCCESS, operation);
}

class AlignedBuffer
{
public:
    explicit AlignedBuffer(size_t bytes)
        : data_(NULL)
        , bytes_(bytes)
    {
        if (bytes != 0)
        {
#if defined(_MSC_VER)
            data_ = _aligned_malloc(bytes, leo2_scratch_alignment());
            if (data_ == NULL)
                throw std::bad_alloc();
#else
            if (posix_memalign(&data_, leo2_scratch_alignment(), bytes) != 0)
                throw std::bad_alloc();
#endif
            memset(data_, 0, bytes);
        }
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
    size_t bytes() const { return bytes_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);

    void* data_;
    size_t bytes_;
};

class GuardedShard
{
public:
    GuardedShard(size_t bytes, unsigned requested_residue)
        : storage_(bytes + 192, kGuardByte)
        , offset_(0)
        , bytes_(bytes)
    {
        const uintptr_t base = reinterpret_cast<uintptr_t>(&storage_[0]);
        for (size_t candidate = 64; candidate < 128; ++candidate)
        {
            if (((base + candidate) & 63u) == (requested_residue & 63u))
            {
                offset_ = candidate;
                break;
            }
        }
        if (offset_ == 0)
            throw std::runtime_error("could not construct guarded shard alignment");
    }

    uint8_t* data() { return &storage_[offset_]; }
    const uint8_t* data() const { return &storage_[offset_]; }
    size_t bytes() const { return bytes_; }

    void fill(uint8_t value)
    {
        memset(data(), value, bytes_);
    }

    bool guards_intact() const
    {
        for (size_t i = 0; i < offset_; ++i)
            if (storage_[i] != kGuardByte)
                return false;
        for (size_t i = offset_ + bytes_; i < storage_.size(); ++i)
            if (storage_[i] != kGuardByte)
                return false;
        return true;
    }

private:
    std::vector<uint8_t> storage_;
    size_t offset_;
    size_t bytes_;
};

class Shards
{
public:
    Shards(unsigned count, size_t bytes, unsigned residue_base)
    {
        shards_.reserve(count);
        for (unsigned i = 0; i < count; ++i)
        {
            unsigned residue = (residue_base + i * 11u) & 63u;
            if (residue == 0)
                residue = 1;
            shards_.push_back(GuardedShard(bytes, residue));
        }
    }

    unsigned size() const { return static_cast<unsigned>(shards_.size()); }
    size_t bytes() const { return shards_.empty() ? 0 : shards_[0].bytes(); }
    GuardedShard& operator[](unsigned i) { return shards_[i]; }
    const GuardedShard& operator[](unsigned i) const { return shards_[i]; }

    std::vector<const void*> const_pointers() const
    {
        std::vector<const void*> result(size());
        for (unsigned i = 0; i < size(); ++i)
            result[i] = shards_[i].data();
        return result;
    }

    std::vector<void*> mutable_pointers()
    {
        std::vector<void*> result(size());
        for (unsigned i = 0; i < size(); ++i)
            result[i] = shards_[i].data();
        return result;
    }

    bool equals(unsigned i, const Shards& other, unsigned j) const
    {
        return bytes() == other.bytes() &&
            memcmp(shards_[i].data(), other.shards_[j].data(), bytes()) == 0;
    }

    bool all_guards_intact() const
    {
        for (unsigned i = 0; i < size(); ++i)
            if (!shards_[i].guards_intact())
                return false;
        return true;
    }

private:
    std::vector<GuardedShard> shards_;
};

void fill_random(Shards& shards, Random& random)
{
    for (unsigned i = 0; i < shards.size(); ++i)
        for (size_t j = 0; j < shards.bytes(); ++j)
            shards[i].data()[j] = static_cast<uint8_t>(random.next());
}

leo2_codec* create_codec(
    leo2_context* context,
    unsigned k,
    unsigned r,
    leo2_profile profile,
    leo2_field field)
{
    leo2_codec* codec = NULL;
    require_success(leo2_codec_create(context, k, r, profile, field, NULL, &codec),
        "codec create");
    require(codec != NULL, "codec create returned a null codec");
    return codec;
}

void encode(
    const leo2_codec* codec,
    const Shards& original,
    Shards& recovery,
    std::vector<void*>* recovery_pointers_override = NULL)
{
    std::vector<const void*> original_pointers = original.const_pointers();
    std::vector<void*> recovery_pointers = recovery.mutable_pointers();
    if (recovery_pointers_override)
        recovery_pointers = *recovery_pointers_override;
    size_t scratch_bytes = 0;
    require_success(leo2_encode_scratch_size(codec, original.bytes(), &scratch_bytes),
        "encode scratch query");
    AlignedBuffer scratch(scratch_bytes);
    require_success(leo2_encode(codec, original.bytes(), &original_pointers[0],
        &recovery_pointers[0], scratch.data(), scratch.bytes()), "encode");
}

void verify_guards(const Shards& first, const Shards& second, TestCounts* counts)
{
    require(first.all_guards_intact() && second.all_guards_intact(),
        "shard guard region was modified");
    counts->guard_checks += first.size() + second.size();
}

std::vector<std::vector<uint8_t> > legacy_encode(
    const Shards& original,
    unsigned recovery_count)
{
    const size_t bytes = original.bytes();
    const size_t rounded = (bytes + 63u) & ~static_cast<size_t>(63u);
    std::vector<std::vector<uint8_t> > padded(
        original.size(), std::vector<uint8_t>(rounded, 0));
    std::vector<const void*> original_pointers(original.size());
    for (unsigned i = 0; i < original.size(); ++i)
    {
        memcpy(&padded[i][0], original[i].data(), bytes);
        original_pointers[i] = &padded[i][0];
    }
    const unsigned work_count = leo_encode_work_count(original.size(), recovery_count);
    std::vector<std::vector<uint8_t> > work(
        work_count, std::vector<uint8_t>(rounded, 0));
    std::vector<void*> work_pointers(work_count);
    for (unsigned i = 0; i < work_count; ++i)
        work_pointers[i] = &work[i][0];
    require(leo_encode(rounded, original.size(), recovery_count, work_count,
        &original_pointers[0], &work_pointers[0]) == Leopard_Success,
        "legacy encoder failed");
    work.resize(recovery_count);
    for (unsigned i = 0; i < recovery_count; ++i)
        work[i].resize(bytes);
    return work;
}

void test_high_legacy_compatibility(
    leo2_context* context,
    Random& random,
    TestCounts* counts)
{
    struct Pair { unsigned k; unsigned r; };
    static const Pair pairs[] = {
        { 1, 1 }, { 2, 1 }, { 3, 2 }, { 4, 4 }, { 7, 3 }, { 8, 8 },
        { 9, 7 }, { 15, 1 }, { 16, 16 }, { 17, 8 }, { 31, 16 },
        { 32, 32 }, { 33, 31 }, { 63, 1 }, { 64, 64 }, { 65, 32 },
        { 127, 64 }, { 128, 128 }, { 129, 1 }, { 191, 64 },
        { 192, 64 }, { 193, 32 }, { 223, 32 }, { 224, 32 }, { 255, 1 }
    };
    static const size_t byte_counts[] = { 1, 2, 7, 17, 33, 64, 65, 129, 257 };

    for (unsigned case_i = 0; case_i < sizeof(pairs) / sizeof(pairs[0]); ++case_i)
    {
        const Pair pair = pairs[case_i];
        const size_t bytes = byte_counts[case_i %
            (sizeof(byte_counts) / sizeof(byte_counts[0]))];
        leo2_codec* codec = create_codec(context, pair.k, pair.r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        require(leo2_codec_field(codec) == LEO2_FIELD_GF8,
            "legacy boundary case unexpectedly left GF8");

        Shards original(pair.k, bytes, 1u + case_i);
        Shards recovery(pair.r, bytes, 3u + case_i);
        fill_random(original, random);
        encode(codec, original, recovery);
        const std::vector<std::vector<uint8_t> > expected =
            legacy_encode(original, pair.r);
        for (unsigned i = 0; i < pair.r; ++i)
            require(memcmp(recovery[i].data(), &expected[i][0], bytes) == 0,
                "GF8 legacy-high parity mismatch at boundary case");

        if (pair.r > 1)
        {
            Shards subset(pair.r, bytes, 5u + case_i);
            for (unsigned i = 0; i < pair.r; ++i)
                subset[i].fill(0xa7);
            std::vector<void*> subset_pointers(pair.r, NULL);
            const unsigned first = random.below(pair.r);
            unsigned second = random.below(pair.r);
            if (second == first)
                second = (second + 1) % pair.r;
            subset_pointers[first] = subset[first].data();
            subset_pointers[second] = subset[second].data();
            encode(codec, original, subset, &subset_pointers);
            require(memcmp(subset[first].data(), &expected[first][0], bytes) == 0 &&
                    memcmp(subset[second].data(), &expected[second][0], bytes) == 0,
                "GF8 requested high parity subset mismatch");
            for (unsigned i = 0; i < pair.r; ++i)
                if (i != first && i != second)
                    for (size_t j = 0; j < bytes; ++j)
                        require(subset[i].data()[j] == 0xa7,
                            "unrequested high parity output was modified");
            verify_guards(original, subset, counts);
            ++counts->subset_cases;
        }
        verify_guards(original, recovery, counts);
        ++counts->high_cases;
        counts->high_bytes_compared += static_cast<uint64_t>(pair.r) * bytes;
        leo2_codec_destroy(codec);
    }
}

void test_low_direct_parity(
    leo2_context* context,
    const BinaryField& field,
    Random& random,
    TestCounts* counts)
{
    struct Pair { unsigned k; unsigned r; size_t bytes; };
    static const Pair pairs[] = {
        { 1, 3, 1 }, { 2, 5, 2 }, { 3, 5, 7 }, { 5, 11, 17 },
        { 7, 9, 3 }, { 9, 17, 7 }, { 17, 33, 2 }
    };
    for (unsigned case_i = 0; case_i < sizeof(pairs) / sizeof(pairs[0]); ++case_i)
    {
        const Pair pair = pairs[case_i];
        leo2_codec* codec = create_codec(context, pair.k, pair.r,
            LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8);
        Shards original(pair.k, pair.bytes, 7u + case_i);
        Shards recovery(pair.r, pair.bytes, 13u + case_i);
        fill_random(original, random);
        encode(codec, original, recovery);

        const ProfileLayout layout = leopard2_test::make_profile_layout(
            leopard2_test::kLow, pair.k, pair.r);
        const Matrix generator = leopard2_test::direct_systematic_generator(field, layout);
        for (size_t byte_i = 0; byte_i < pair.bytes; ++byte_i)
        {
            std::vector<Element> message(pair.k);
            for (unsigned i = 0; i < pair.k; ++i)
                message[i] = original[i].data()[byte_i];
            const std::vector<Element> expected =
                leopard2_test::matrix_vector_multiply(field, generator, message);
            for (unsigned i = 0; i < pair.r; ++i)
            {
                require(recovery[i].data()[byte_i] == expected[pair.k + i],
                    "GF8 low parity differs from direct generator oracle");
                ++counts->low_symbols_compared;
            }
        }

        if (pair.r > 2)
        {
            Shards subset(pair.r, pair.bytes, 17u + case_i);
            for (unsigned i = 0; i < pair.r; ++i)
                subset[i].fill(0x4e);
            std::vector<void*> pointers(pair.r, NULL);
            pointers[0] = subset[0].data();
            pointers[pair.r - 1] = subset[pair.r - 1].data();
            encode(codec, original, subset, &pointers);
            require(subset.equals(0, recovery, 0) &&
                    subset.equals(pair.r - 1, recovery, pair.r - 1),
                "GF8 requested low parity subset mismatch");
            for (unsigned i = 1; i + 1 < pair.r; ++i)
                for (size_t j = 0; j < pair.bytes; ++j)
                    require(subset[i].data()[j] == 0x4e,
                        "unrequested low parity output was modified");
            ++counts->subset_cases;
            verify_guards(original, subset, counts);
        }

        verify_guards(original, recovery, counts);
        ++counts->low_cases;
        leo2_codec_destroy(codec);
    }
}

void shuffle(std::vector<unsigned>& values, Random& random)
{
    for (size_t i = values.size(); i > 1; --i)
        std::swap(values[i - 1], values[random.below(static_cast<unsigned>(i))]);
}

void run_random_decode_case(
    leo2_context* context,
    Random& random,
    unsigned iteration,
    TestCounts* counts)
{
    struct Shape { unsigned k; unsigned r; leo2_profile profile; };
    static const Shape shapes[] = {
        { 3, 2, LEO2_PROFILE_LEGACY_HIGH_V1 },
        { 7, 3, LEO2_PROFILE_LEGACY_HIGH_V1 },
        { 9, 7, LEO2_PROFILE_LEGACY_HIGH_V1 },
        { 17, 8, LEO2_PROFILE_LEGACY_HIGH_V1 },
        { 31, 16, LEO2_PROFILE_LEGACY_HIGH_V1 },
        { 65, 32, LEO2_PROFILE_LEGACY_HIGH_V1 },
        { 129, 1, LEO2_PROFILE_LEGACY_HIGH_V1 },
        { 191, 64, LEO2_PROFILE_LEGACY_HIGH_V1 },
        { 224, 32, LEO2_PROFILE_LEGACY_HIGH_V1 },
        { 1, 5, LEO2_PROFILE_LOW_V1 },
        { 2, 5, LEO2_PROFILE_LOW_V1 },
        { 3, 7, LEO2_PROFILE_LOW_V1 },
        { 5, 11, LEO2_PROFILE_LOW_V1 },
        { 9, 17, LEO2_PROFILE_LOW_V1 },
        { 17, 33, LEO2_PROFILE_LOW_V1 },
        { 33, 65, LEO2_PROFILE_LOW_V1 }
    };
    static const size_t byte_counts[] = {
        1, 2, 3, 7, 8, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 129, 257
    };

    const Shape shape = shapes[random.below(
        static_cast<unsigned>(sizeof(shapes) / sizeof(shapes[0])))];
    const size_t bytes = byte_counts[random.below(
        static_cast<unsigned>(sizeof(byte_counts) / sizeof(byte_counts[0])))];
    leo2_codec* codec = create_codec(context, shape.k, shape.r,
        shape.profile, LEO2_FIELD_GF8);
    Shards original(shape.k, bytes, 1u + iteration);
    Shards recovery(shape.r, bytes, 23u + iteration);
    fill_random(original, random);
    encode(codec, original, recovery);

    std::vector<unsigned> coordinates(shape.k + shape.r);
    for (unsigned i = 0; i < coordinates.size(); ++i)
        coordinates[i] = i;
    shuffle(coordinates, random);
    unsigned loss_count = 1 + random.below(shape.r);
    if (iteration % 8 == 0)
        loss_count = shape.r;
    else if (iteration % 8 == 1)
        loss_count = 1;
    if (iteration % 4 != 0)
    {
        const unsigned forced_original = random.below(shape.k);
        std::vector<unsigned>::iterator location = std::find(
            coordinates.begin(), coordinates.end(), forced_original);
        std::swap(coordinates[0], *location);
    }

    std::vector<uint8_t> original_present(shape.k, 1);
    std::vector<uint8_t> recovery_present(shape.r, 1);
    for (unsigned i = 0; i < loss_count; ++i)
    {
        if (coordinates[i] < shape.k)
            original_present[coordinates[i]] = 0;
        else
            recovery_present[coordinates[i] - shape.k] = 0;
    }

    unsigned missing_original_count = 0;
    for (unsigned i = 0; i < shape.k; ++i)
        if (!original_present[i])
            ++missing_original_count;

    leo2_decode_plan* plan = NULL;
    require_success(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), "random decode plan create");
    require(leo2_decode_plan_missing_original_count(plan) == missing_original_count,
        "random plan reports the wrong missing-original count");

    std::vector<const void*> original_pointers = original.const_pointers();
    std::vector<const void*> recovery_pointers = recovery.const_pointers();
    Shards restored(shape.k, bytes, 37u + iteration);
    std::vector<void*> restored_pointers(shape.k, NULL);
    for (unsigned i = 0; i < shape.k; ++i)
    {
        if (!original_present[i])
        {
            original_pointers[i] = NULL;
            restored_pointers[i] = restored[i].data();
        }
    }
    for (unsigned i = 0; i < shape.r; ++i)
        if (!recovery_present[i])
            recovery_pointers[i] = NULL;

    size_t scratch_bytes = 0;
    require_success(leo2_decode_plan_scratch_size(plan, bytes, &scratch_bytes),
        "random decode scratch query");
    AlignedBuffer scratch(scratch_bytes);
    for (unsigned repeat = 0; repeat < 2; ++repeat)
    {
        require_success(leo2_decode_plan_execute(plan, bytes, &original_pointers[0],
            &recovery_pointers[0], &restored_pointers[0], scratch.data(), scratch.bytes()),
            "random decode plan execute");
        ++counts->repeated_executions;
        for (unsigned i = 0; i < shape.k; ++i)
        {
            if (!original_present[i])
            {
                if (!restored.equals(i, original, i))
                {
                    std::ostringstream stream;
                    stream << "random recovery mismatch: iteration=" << iteration
                           << " K=" << shape.k << " R=" << shape.r
                           << " profile=" << static_cast<int>(shape.profile)
                           << " bytes=" << bytes << " original=" << i;
                    throw std::runtime_error(stream.str());
                }
                ++counts->recovered_shards;
            }
        }
    }

    if (missing_original_count != 0 && iteration % 7 == 0)
    {
        Shards one_shot_restored(shape.k, bytes, 41u + iteration);
        std::vector<void*> one_shot_pointers(shape.k, NULL);
        for (unsigned i = 0; i < shape.k; ++i)
            if (!original_present[i])
                one_shot_pointers[i] = one_shot_restored[i].data();
        size_t one_shot_scratch_bytes = 0;
        require_success(leo2_decode_scratch_size(codec, bytes, &one_shot_scratch_bytes),
            "one-shot scratch query");
        AlignedBuffer one_shot_scratch(one_shot_scratch_bytes);
        require_success(leo2_decode(codec, bytes, &original_present[0],
            &recovery_present[0], &original_pointers[0], &recovery_pointers[0],
            &one_shot_pointers[0], one_shot_scratch.data(), one_shot_scratch.bytes()),
            "one-shot decode");
        for (unsigned i = 0; i < shape.k; ++i)
            if (!original_present[i])
                require(one_shot_restored.equals(i, original, i),
                    "one-shot recovery mismatch");
        require(one_shot_restored.all_guards_intact(),
            "one-shot output guard was modified");
        counts->guard_checks += shape.k;
    }

    if (missing_original_count != 0 && iteration % 4 == 0)
    {
        leo2_codec_options generic_options;
        memset(&generic_options, 0, sizeof(generic_options));
        generic_options.struct_size = sizeof(generic_options);
        generic_options.flags = LEO2_CODEC_FORCE_GENERIC_DECODE;
        leo2_codec* generic_codec = NULL;
        require_success(leo2_codec_create(context, shape.k, shape.r,
            shape.profile, LEO2_FIELD_GF8, &generic_options, &generic_codec),
            "generic-fallback codec create");
        leo2_decode_plan* generic_plan = NULL;
        require_success(leo2_decode_plan_create(generic_codec, &original_present[0],
            &recovery_present[0], &generic_plan), "generic-fallback plan create");
        size_t generic_scratch_bytes = 0;
        require_success(leo2_decode_plan_scratch_size(
            generic_plan, bytes, &generic_scratch_bytes),
            "generic-fallback scratch query");
        AlignedBuffer generic_scratch(generic_scratch_bytes);
        Shards generic_restored(shape.k, bytes, 43u + iteration);
        std::vector<void*> generic_restored_pointers(shape.k, NULL);
        for (unsigned i = 0; i < shape.k; ++i)
            if (!original_present[i])
                generic_restored_pointers[i] = generic_restored[i].data();
        require_success(leo2_decode_plan_execute(generic_plan, bytes,
            &original_pointers[0], &recovery_pointers[0],
            &generic_restored_pointers[0], generic_scratch.data(),
            generic_scratch.bytes()), "generic-fallback execute");
        for (unsigned i = 0; i < shape.k; ++i)
        {
            if (!original_present[i])
            {
                require(generic_restored.equals(i, restored, i),
                    "specialized and generic fallback recovery differ");
                ++counts->generic_comparisons;
            }
        }
        require(generic_restored.all_guards_intact(),
            "generic-fallback output guard was modified");
        counts->guard_checks += shape.k;
        leo2_decode_plan_destroy(generic_plan);
        leo2_codec_destroy(generic_codec);
    }

    verify_guards(original, recovery, counts);
    require(restored.all_guards_intact(), "restored shard guard region was modified");
    counts->guard_checks += shape.k;
    ++counts->randomized_patterns;
    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
}

void test_no_loss_and_parity_only(leo2_context* context, TestCounts* counts)
{
    struct Shape { unsigned k; unsigned r; leo2_profile profile; };
    static const Shape shapes[] = {
        { 9, 7, LEO2_PROFILE_LEGACY_HIGH_V1 },
        { 3, 7, LEO2_PROFILE_LOW_V1 }
    };
    for (unsigned shape_i = 0; shape_i < sizeof(shapes) / sizeof(shapes[0]); ++shape_i)
    {
        const Shape shape = shapes[shape_i];
        leo2_codec* codec = create_codec(context, shape.k, shape.r,
            shape.profile, LEO2_FIELD_GF8);
        std::vector<uint8_t> original_present(shape.k, 1);
        std::vector<uint8_t> recovery_present(shape.r, 1);
        for (unsigned mode = 0; mode < 2; ++mode)
        {
            if (mode == 0)
                std::fill(recovery_present.begin(), recovery_present.end(), 0);
            else
            {
                std::fill(recovery_present.begin(), recovery_present.end(), 1);
                recovery_present[0] = 0;
                recovery_present[shape.r - 1] = 0;
            }
            leo2_decode_plan* plan = NULL;
            require_success(leo2_decode_plan_create(codec, &original_present[0],
                &recovery_present[0], &plan), "no-op plan create");
            require(leo2_decode_plan_missing_original_count(plan) == 0,
                "no-op plan reports missing originals");
            size_t scratch_bytes = 99;
            require_success(leo2_decode_plan_scratch_size(plan, 17, &scratch_bytes),
                "no-op scratch query");
            require(scratch_bytes == 0, "no-op plan unexpectedly requires scratch");
            require_success(leo2_decode_plan_execute(plan, 17, NULL, NULL, NULL, NULL, 0),
                "no-op plan execute");
            leo2_decode_plan_destroy(plan);
            ++counts->no_op_cases;
        }
        leo2_codec_destroy(codec);
    }
}

struct BatchStripe
{
    std::unique_ptr<Shards> original;
    std::unique_ptr<Shards> recovery;
    std::unique_ptr<Shards> restored;
    std::vector<const void*> original_pointers;
    std::vector<void*> recovery_output_pointers;
    std::vector<const void*> recovery_input_pointers;
    std::vector<void*> restored_pointers;
    std::unique_ptr<AlignedBuffer> encode_scratch;
    std::unique_ptr<AlignedBuffer> decode_scratch;
};

void test_batch_apis(leo2_context* context, Random& random, TestCounts* counts)
{
    const unsigned k = 7;
    const unsigned r = 3;
    const size_t byte_counts[] = { 1, 17, 65, 129 };
    leo2_codec* codec = create_codec(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    std::vector<BatchStripe> stripes(sizeof(byte_counts) / sizeof(byte_counts[0]));
    std::vector<leo2_encode_batch_item> encode_items(stripes.size());
    for (unsigned stripe_i = 0; stripe_i < stripes.size(); ++stripe_i)
    {
        BatchStripe& stripe = stripes[stripe_i];
        const size_t bytes = byte_counts[stripe_i];
        stripe.original.reset(new Shards(k, bytes, 3u + stripe_i));
        stripe.recovery.reset(new Shards(r, bytes, 19u + stripe_i));
        stripe.restored.reset(new Shards(k, bytes, 31u + stripe_i));
        fill_random(*stripe.original, random);
        stripe.original_pointers = stripe.original->const_pointers();
        stripe.recovery_output_pointers = stripe.recovery->mutable_pointers();
        size_t scratch_bytes = 0;
        require_success(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
            "batch encode scratch query");
        stripe.encode_scratch.reset(new AlignedBuffer(scratch_bytes));
        leo2_encode_batch_item& item = encode_items[stripe_i];
        item.shard_bytes = bytes;
        item.original = &stripe.original_pointers[0];
        item.recovery = &stripe.recovery_output_pointers[0];
        item.scratch = stripe.encode_scratch->data();
        item.scratch_bytes = stripe.encode_scratch->bytes();
    }
    require_success(leo2_encode_batch(codec, &encode_items[0], encode_items.size()),
        "batch encode");

    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 1);
    original_present[1] = 0;
    original_present[5] = 0;
    recovery_present[2] = 0;
    leo2_decode_plan* plan = NULL;
    require_success(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), "batch decode plan create");

    std::vector<leo2_decode_batch_item> decode_items(stripes.size());
    for (unsigned stripe_i = 0; stripe_i < stripes.size(); ++stripe_i)
    {
        BatchStripe& stripe = stripes[stripe_i];
        stripe.original_pointers = stripe.original->const_pointers();
        stripe.recovery_input_pointers = stripe.recovery->const_pointers();
        stripe.restored_pointers.assign(k, NULL);
        stripe.original_pointers[1] = NULL;
        stripe.original_pointers[5] = NULL;
        stripe.recovery_input_pointers[2] = NULL;
        stripe.restored_pointers[1] = (*stripe.restored)[1].data();
        stripe.restored_pointers[5] = (*stripe.restored)[5].data();
        size_t scratch_bytes = 0;
        require_success(leo2_decode_plan_scratch_size(
            plan, byte_counts[stripe_i], &scratch_bytes), "batch decode scratch query");
        stripe.decode_scratch.reset(new AlignedBuffer(scratch_bytes));
        leo2_decode_batch_item& item = decode_items[stripe_i];
        item.shard_bytes = byte_counts[stripe_i];
        item.original = &stripe.original_pointers[0];
        item.recovery = &stripe.recovery_input_pointers[0];
        item.restored_original = &stripe.restored_pointers[0];
        item.scratch = stripe.decode_scratch->data();
        item.scratch_bytes = stripe.decode_scratch->bytes();
    }
    require_success(leo2_decode_plan_execute_batch(plan, &decode_items[0],
        decode_items.size()), "batch decode");
    for (unsigned stripe_i = 0; stripe_i < stripes.size(); ++stripe_i)
    {
        BatchStripe& stripe = stripes[stripe_i];
        require(stripe.restored->equals(1, *stripe.original, 1) &&
                stripe.restored->equals(5, *stripe.original, 5),
            "batch recovered data mismatch");
        verify_guards(*stripe.original, *stripe.recovery, counts);
        require(stripe.restored->all_guards_intact(),
            "batch restored output guard was modified");
        counts->guard_checks += k;
        ++counts->batch_stripes;
    }

    require_success(leo2_encode_batch(codec, NULL, 0), "empty encode batch");
    require_success(leo2_decode_plan_execute_batch(plan, NULL, 0), "empty decode batch");
    require_result(leo2_encode_batch(codec, NULL, 1), LEO2_INVALID_ARGUMENT,
        "null encode batch");
    ++counts->invalid_checks;
    require_result(leo2_decode_plan_execute_batch(plan, NULL, 1), LEO2_INVALID_ARGUMENT,
        "null decode batch");
    ++counts->invalid_checks;
    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
}

void test_invalid_contracts(leo2_context* context, TestCounts* counts)
{
    leo2_codec* codec = reinterpret_cast<leo2_codec*>(static_cast<uintptr_t>(1));
    require_result(leo2_codec_create(context, 0, 1, LEO2_PROFILE_AUTO,
        LEO2_FIELD_AUTO, NULL, &codec), LEO2_INVALID_ARGUMENT, "zero K");
    ++counts->invalid_checks;
    require_result(leo2_codec_create(context, 1, 0, LEO2_PROFILE_AUTO,
        LEO2_FIELD_AUTO, NULL, &codec), LEO2_INVALID_ARGUMENT, "zero R");
    ++counts->invalid_checks;
    require_result(leo2_codec_create(context, UINT32_MAX, UINT32_MAX,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, NULL, &codec),
        LEO2_INVALID_COUNTS, "overflowing counts");
    ++counts->invalid_checks;
    require_result(leo2_codec_create(context, 256, 1,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        LEO2_INVALID_COUNTS, "GF8 parent overflow");
    ++counts->invalid_checks;
    require_result(leo2_codec_create(context, 3, 2,
        LEO2_PROFILE_EXACT_EXPERIMENTAL_V1, LEO2_FIELD_GF8, NULL, &codec),
        LEO2_UNSUPPORTED, "experimental profile");
    ++counts->invalid_checks;

    codec = create_codec(context, 3, 2,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    size_t scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(codec, 0, &scratch_bytes),
        LEO2_INVALID_ARGUMENT, "zero-byte encode scratch");
    ++counts->invalid_checks;
    require_result(leo2_encode_scratch_size(codec, UINT64_MAX, &scratch_bytes),
        LEO2_INVALID_ARGUMENT, "overflowing encode scratch");
    ++counts->invalid_checks;
    require_result(leo2_decode_scratch_size(codec, UINT64_MAX, &scratch_bytes),
        LEO2_INVALID_ARGUMENT, "overflowing decode scratch");
    ++counts->invalid_checks;

    const size_t bytes = 17;
    Shards original(3, bytes, 1);
    Shards recovery(2, bytes, 5);
    Random random(UINT64_C(0x1234));
    fill_random(original, random);
    std::vector<const void*> original_pointers = original.const_pointers();
    std::vector<void*> recovery_pointers = recovery.mutable_pointers();
    require_success(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        "invalid-case scratch query");
    AlignedBuffer scratch(scratch_bytes + leo2_scratch_alignment());
    require_result(leo2_encode(codec, bytes, &original_pointers[0],
        &recovery_pointers[0], scratch.data(), scratch_bytes - 1),
        LEO2_SCRATCH_TOO_SMALL, "short encode scratch");
    ++counts->invalid_checks;
    require_result(leo2_encode(codec, bytes, &original_pointers[0],
        &recovery_pointers[0], static_cast<uint8_t*>(scratch.data()) + 1, scratch_bytes),
        LEO2_BAD_ALIGNMENT, "unaligned encode scratch");
    ++counts->invalid_checks;
    recovery_pointers[0] = const_cast<uint8_t*>(original[0].data());
    require_result(leo2_encode(codec, bytes, &original_pointers[0],
        &recovery_pointers[0], scratch.data(), scratch_bytes),
        LEO2_OVERLAP, "encode input/output overlap");
    ++counts->invalid_checks;
    recovery_pointers = recovery.mutable_pointers();
    recovery_pointers[1] = recovery_pointers[0];
    require_result(leo2_encode(codec, bytes, &original_pointers[0],
        &recovery_pointers[0], scratch.data(), scratch_bytes),
        LEO2_OVERLAP, "encode output/output overlap");
    ++counts->invalid_checks;
    recovery_pointers = recovery.mutable_pointers();
    require_success(leo2_encode(codec, bytes, &original_pointers[0],
        &recovery_pointers[0], scratch.data(), scratch_bytes),
        "invalid-case initial encode");

    std::vector<uint8_t> original_present(3, 1);
    std::vector<uint8_t> recovery_present(2, 0);
    original_present[0] = 0;
    leo2_decode_plan* plan = NULL;
    require_result(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), LEO2_NEED_MORE_DATA, "insufficient plan");
    ++counts->invalid_checks;
    original_present[0] = 2;
    recovery_present[0] = 1;
    require_result(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), LEO2_INVALID_ARGUMENT, "non-binary presence");
    ++counts->invalid_checks;

    original_present.assign(3, 1);
    recovery_present.assign(2, 1);
    original_present[0] = 0;
    recovery_present[1] = 0;
    require_success(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), "invalid-case plan create");
    size_t decode_scratch_bytes = 0;
    require_success(leo2_decode_plan_scratch_size(plan, bytes, &decode_scratch_bytes),
        "invalid-case decode scratch query");
    AlignedBuffer decode_scratch(decode_scratch_bytes);
    std::vector<const void*> decode_original = original.const_pointers();
    std::vector<const void*> decode_recovery = recovery.const_pointers();
    std::vector<void*> restored(3, NULL);
    decode_original[0] = NULL;
    decode_recovery[1] = NULL;
    restored[0] = const_cast<uint8_t*>(recovery[0].data());
    require_result(leo2_decode_plan_execute(plan, bytes, &decode_original[0],
        &decode_recovery[0], &restored[0], decode_scratch.data(), decode_scratch.bytes()),
        LEO2_OVERLAP, "decode input/output overlap");
    ++counts->invalid_checks;
    restored[0] = original[0].data();
    decode_original[1] = NULL;
    require_result(leo2_decode_plan_execute(plan, bytes, &decode_original[0],
        &decode_recovery[0], &restored[0], decode_scratch.data(), decode_scratch.bytes()),
        LEO2_INVALID_ARGUMENT, "decode presence/pointer mismatch");
    ++counts->invalid_checks;

    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
}

uint64_t parse_unsigned(const char* text, const char* label)
{
    if (!text || !*text)
        throw std::invalid_argument(std::string("empty ") + label);
    errno = 0;
    char* end = NULL;
    const unsigned long long value = strtoull(text, &end, 0);
    if (errno != 0 || !end || *end != '\0')
        throw std::invalid_argument(std::string("invalid ") + label + ": " + text);
    return static_cast<uint64_t>(value);
}

struct Configuration
{
    uint64_t seed;
    unsigned random_cases;
    uint32_t thread_count;
};

Configuration configuration_from_args(int argc, char** argv)
{
    Configuration configuration;
    configuration.seed = UINT64_C(0x4c656f7061726432);
    configuration.random_cases = 48;
    configuration.thread_count = 1;
    const char* seed_environment = getenv("LEO2_RANDOM_SEED");
    const char* cases_environment = getenv("LEO2_RANDOM_CASES");
    if (seed_environment)
        configuration.seed = parse_unsigned(seed_environment, "LEO2_RANDOM_SEED");
    if (cases_environment)
    {
        const uint64_t count = parse_unsigned(cases_environment, "LEO2_RANDOM_CASES");
        if (count > 1000000u)
            throw std::invalid_argument("random case count is unreasonably large");
        configuration.random_cases = static_cast<unsigned>(count);
    }

    for (int i = 1; i < argc; ++i)
    {
        const std::string argument(argv[i]);
        if (argument == "--seed" && i + 1 < argc)
            configuration.seed = parse_unsigned(argv[++i], "seed");
        else if (argument == "--cases" && i + 1 < argc)
        {
            const uint64_t count = parse_unsigned(argv[++i], "cases");
            if (count > 1000000u)
                throw std::invalid_argument("random case count is unreasonably large");
            configuration.random_cases = static_cast<unsigned>(count);
        }
        else if (argument == "--threads" && i + 1 < argc)
        {
            const uint64_t count = parse_unsigned(argv[++i], "threads");
            if (count == 0 || count > 128)
                throw std::invalid_argument("thread count must be in 1..128");
            configuration.thread_count = static_cast<uint32_t>(count);
        }
        else
            throw std::invalid_argument(
                "usage: test_random [--seed N] [--cases N] [--threads N]");
    }
    return configuration;
}

} // namespace

int main(int argc, char** argv)
{
    uint64_t seed_for_failure = 0;
    try
    {
        const Configuration configuration = configuration_from_args(argc, argv);
        seed_for_failure = configuration.seed;
        std::cout << "leopard2_random seed=" << configuration.seed
                  << " random_cases=" << configuration.random_cases
                  << " threads=" << configuration.thread_count << std::endl;
        Random random(configuration.seed);
        TestCounts counts;

        leo2_context_options options;
        memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AUTO;
        options.thread_count = configuration.thread_count;
        leo2_context* context = NULL;
        require_success(leo2_context_create(&options, &context), "context create");
        require(context != NULL, "context create returned null");

        test_high_legacy_compatibility(context, random, &counts);
        const BinaryField gf8 = leopard2_test::make_legacy_gf8();
        test_low_direct_parity(context, gf8, random, &counts);
        test_no_loss_and_parity_only(context, &counts);
        test_batch_apis(context, random, &counts);
        for (unsigned i = 0; i < configuration.random_cases; ++i)
            run_random_decode_case(context, random, i, &counts);
        test_invalid_contracts(context, &counts);
        leo2_context_destroy(context);

        std::cout << "leopard2_random passed"
                  << " seed=" << configuration.seed
                  << " high_cases=" << counts.high_cases
                  << " high_bytes_compared=" << counts.high_bytes_compared
                  << " low_cases=" << counts.low_cases
                  << " low_symbols_compared=" << counts.low_symbols_compared
                  << " randomized_patterns=" << counts.randomized_patterns
                  << " recovered_shards=" << counts.recovered_shards
                  << " repeated_executions=" << counts.repeated_executions
                  << " generic_comparisons=" << counts.generic_comparisons
                  << " subset_cases=" << counts.subset_cases
                  << " batch_stripes=" << counts.batch_stripes
                  << " no_op_cases=" << counts.no_op_cases
                  << " invalid_checks=" << counts.invalid_checks
                  << " guard_checks=" << counts.guard_checks << std::endl;
        return 0;
    }
    catch (const std::exception& exception)
    {
        std::cerr << "leopard2_random failed seed=" << seed_for_failure
                  << ": " << exception.what() << std::endl;
        return 1;
    }
}
