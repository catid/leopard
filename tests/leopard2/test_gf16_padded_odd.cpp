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

#include <algorithm>
#include <stddef.h>
#include <stdint.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

namespace {

using leopard2_test::BinaryField;
using leopard2_test::Element;
using leopard2_test::Matrix;
using leopard2_test::ProfileLayout;

typedef std::vector<uint8_t> Bytes;
typedef std::vector<Bytes> Shards;

struct Counts
{
    uint64_t payload_cases;
    uint64_t direct_symbols;
    uint64_t restored_shards;
    uint64_t mds_basis_recoveries;
    uint64_t exhaustive_patterns;
    uint64_t guard_checks;

    Counts()
        : payload_cases(0)
        , direct_symbols(0)
        , restored_shards(0)
        , mds_basis_recoveries(0)
        , exhaustive_patterns(0)
        , guard_checks(0)
    {}
};

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void require_result(leo2_result result, const std::string& operation)
{
    if (result != LEO2_SUCCESS)
    {
        std::ostringstream stream;
        stream << operation << ": " << leo2_result_string(result)
               << " (" << static_cast<int>(result) << ')';
        throw std::runtime_error(stream.str());
    }
}

struct AlignedBuffer
{
    void* data;
    size_t bytes;

    explicit AlignedBuffer(size_t size)
        : data(NULL)
        , bytes(size)
    {
        if (size != 0)
        {
#if defined(_MSC_VER)
            data = _aligned_malloc(size, leo2_scratch_alignment());
            if (!data)
                throw std::bad_alloc();
#else
            if (posix_memalign(&data, leo2_scratch_alignment(), size) != 0)
                throw std::bad_alloc();
#endif
            memset(data, 0, size);
        }
    }

    ~AlignedBuffer()
    {
#if defined(_MSC_VER)
        _aligned_free(data);
#else
        free(data);
#endif
    }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);
};

struct CodecOwner
{
    leo2_codec* codec;

    CodecOwner()
        : codec(NULL)
    {}

    ~CodecOwner()
    {
        leo2_codec_destroy(codec);
    }

private:
    CodecOwner(const CodecOwner&);
    CodecOwner& operator=(const CodecOwner&);
};

class GuardedShards
{
public:
    GuardedShards(size_t count, size_t bytes, uint8_t guard, uint8_t fill)
        : bytes_(bytes)
        , guard_(guard)
        , storage_(count, Bytes(kPrefix + bytes + kSuffix, guard))
    {
        for (size_t i = 0; i < storage_.size(); ++i)
            std::fill(storage_[i].begin() + kPrefix,
                storage_[i].begin() + kPrefix + bytes_, fill);
    }

    GuardedShards(const Shards& shards, uint8_t guard)
        : bytes_(shards.empty() ? 0 : shards[0].size())
        , guard_(guard)
        , storage_(shards.size(), Bytes(kPrefix + bytes_ + kSuffix, guard))
    {
        for (size_t i = 0; i < shards.size(); ++i)
        {
            require(shards[i].size() == bytes_,
                "guarded shards have inconsistent lengths");
            std::copy(shards[i].begin(), shards[i].end(),
                storage_[i].begin() + kPrefix);
        }
    }

    const void* data(size_t index) const
    {
        return &storage_[index][kPrefix];
    }

    void* mutable_data(size_t index)
    {
        return &storage_[index][kPrefix];
    }

    Bytes shard(size_t index) const
    {
        return Bytes(storage_[index].begin() + kPrefix,
            storage_[index].begin() + kPrefix + bytes_);
    }

    Shards shards() const
    {
        Shards result(storage_.size());
        for (size_t i = 0; i < result.size(); ++i)
            result[i] = shard(i);
        return result;
    }

    bool guards_intact() const
    {
        for (size_t i = 0; i < storage_.size(); ++i)
        {
            for (size_t j = 0; j < kPrefix; ++j)
                if (storage_[i][j] != guard_)
                    return false;
            for (size_t j = kPrefix + bytes_; j < storage_[i].size(); ++j)
                if (storage_[i][j] != guard_)
                    return false;
        }
        return true;
    }

private:
    enum { kPrefix = 3, kSuffix = 5 };

    size_t bytes_;
    uint8_t guard_;
    std::vector<Bytes> storage_;
};

std::vector<const void*> const_pointers(const Shards& shards)
{
    std::vector<const void*> result(shards.size(), NULL);
    for (size_t i = 0; i < shards.size(); ++i)
        result[i] = &shards[i][0];
    return result;
}

std::vector<void*> mutable_pointers(Shards& shards)
{
    std::vector<void*> result(shards.size(), NULL);
    for (size_t i = 0; i < shards.size(); ++i)
        result[i] = &shards[i][0];
    return result;
}

leo2_codec* create_codec(
    leo2_context* context,
    unsigned k,
    unsigned r,
    leo2_profile profile,
    leo2_shard_layout layout,
    uint32_t flags)
{
    leo2_codec_options options;
    memset(&options, 0, sizeof(options));
    options.struct_size = sizeof(options);
    options.flags = flags;
    options.shard_layout = layout;
    leo2_codec* codec = NULL;
    require_result(leo2_codec_create(context, k, r, profile, LEO2_FIELD_GF16,
        &options, &codec), "codec create");
    require(codec != NULL, "codec create returned null");
    require(leo2_codec_shard_layout(codec) == layout,
        "codec shard-layout introspection mismatch");
    return codec;
}

Bytes deterministic_bytes(size_t count, uint64_t seed)
{
    Bytes result(count, 0);
    uint64_t state = seed;
    for (size_t i = 0; i < count; ++i)
    {
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
        result[i] = static_cast<uint8_t>(
            (state * UINT64_C(2685821657736338717)) >> 56);
    }
    return result;
}

Shards deterministic_payloads(unsigned count, size_t bytes, uint64_t seed)
{
    Shards result(count);
    for (unsigned i = 0; i < count; ++i)
        result[i] = deterministic_bytes(bytes, seed ^
            (UINT64_C(0x9e3779b97f4a7c15) * (i + 1)));
    return result;
}

size_t rounded_bytes(size_t bytes)
{
    return (bytes + 63u) & ~static_cast<size_t>(63u);
}

Bytes compact_pack(const Bytes& input)
{
    require(!input.empty() && (input.size() & 1u) == 0,
        "compact GF16 input must contain complete symbols");
    Bytes result(rounded_bytes(input.size()), 0);
    const size_t complete = input.size() & ~static_cast<size_t>(63u);
    std::copy(input.begin(), input.begin() + complete, result.begin());
    const size_t symbols = (input.size() - complete) / 2;
    if (symbols != 0)
    {
        std::copy(input.begin() + complete,
            input.begin() + complete + symbols,
            result.begin() + complete);
        std::copy(input.begin() + complete + symbols, input.end(),
            result.begin() + complete + 32);
    }
    return result;
}

Bytes compact_gather(const Bytes& input, size_t bytes)
{
    require(bytes != 0 && (bytes & 1u) == 0 &&
        input.size() == rounded_bytes(bytes), "invalid compact GF16 gather");
    Bytes result(bytes, 0);
    const size_t complete = bytes & ~static_cast<size_t>(63u);
    std::copy(input.begin(), input.begin() + complete, result.begin());
    const size_t symbols = (bytes - complete) / 2;
    if (symbols != 0)
    {
        std::copy(input.begin() + complete,
            input.begin() + complete + symbols,
            result.begin() + complete);
        std::copy(input.begin() + complete + 32,
            input.begin() + complete + 32 + symbols,
            result.begin() + complete + symbols);
    }
    return result;
}

Element load_altmap(const Bytes& bytes, size_t tile, size_t lane)
{
    return static_cast<Element>(bytes[tile + lane] |
        (static_cast<unsigned>(bytes[tile + 32 + lane]) << 8));
}

void store_altmap(Bytes& bytes, size_t tile, size_t lane, Element value)
{
    bytes[tile + lane] = static_cast<uint8_t>(value);
    bytes[tile + 32 + lane] = static_cast<uint8_t>(value >> 8);
}

Shards direct_parity(
    const BinaryField& field,
    const Matrix& generator,
    unsigned k,
    unsigned r,
    const Shards& wire)
{
    const size_t bytes = wire[0].size();
    const size_t rounded = rounded_bytes(bytes);
    Shards packed(k);
    for (unsigned i = 0; i < k; ++i)
        packed[i] = compact_pack(wire[i]);
    Shards parity_altmap(r, Bytes(rounded, 0));
    for (size_t tile = 0; tile < rounded; tile += 64)
        for (size_t lane = 0; lane < 32; ++lane)
        {
            std::vector<Element> message(k, 0);
            for (unsigned i = 0; i < k; ++i)
                message[i] = load_altmap(packed[i], tile, lane);
            const std::vector<Element> codeword =
                leopard2_test::matrix_vector_multiply(field, generator, message);
            for (unsigned i = 0; i < r; ++i)
                store_altmap(parity_altmap[i], tile, lane, codeword[k + i]);
        }

    Shards result(r);
    for (unsigned i = 0; i < r; ++i)
        result[i] = compact_gather(parity_altmap[i], bytes);
    return result;
}

Shards direct_gf8_parity(
    const BinaryField& field,
    const Matrix& generator,
    unsigned k,
    unsigned r,
    const Shards& original)
{
    const size_t bytes = original[0].size();
    Shards result(r, Bytes(bytes, 0));
    for (size_t byte = 0; byte < bytes; ++byte)
    {
        std::vector<Element> message(k, 0);
        for (unsigned i = 0; i < k; ++i)
            message[i] = original[i][byte];
        const std::vector<Element> codeword =
            leopard2_test::matrix_vector_multiply(field, generator, message);
        for (unsigned i = 0; i < r; ++i)
            result[i][byte] = static_cast<uint8_t>(codeword[k + i]);
    }
    return result;
}

Shards encode(const leo2_codec* codec, const Shards& wire)
{
    const size_t bytes = wire[0].size();
    const unsigned r = leo2_codec_recovery_count(codec);
    Shards parity(r, Bytes(bytes, 0));
    std::vector<const void*> originals = const_pointers(wire);
    std::vector<void*> recoveries = mutable_pointers(parity);
    size_t scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        "encode scratch query");
    AlignedBuffer scratch(scratch_bytes);
    require_result(leo2_encode(codec, bytes, &originals[0], &recoveries[0],
        scratch.data, scratch.bytes), "encode");
    return parity;
}

Shards encode_guarded(
    const leo2_codec* codec,
    const Shards& wire,
    Counts* counts)
{
    const size_t bytes = wire[0].size();
    const unsigned r = leo2_codec_recovery_count(codec);
    GuardedShards original(wire, 0xa7);
    GuardedShards parity(r, bytes, 0xd3, 0x5c);
    std::vector<const void*> original_pointers(wire.size(), NULL);
    std::vector<void*> recovery_pointers(r, NULL);
    for (size_t i = 0; i < wire.size(); ++i)
        original_pointers[i] = original.data(i);
    for (unsigned i = 0; i < r; ++i)
        recovery_pointers[i] = parity.mutable_data(i);

    size_t scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        "guarded encode scratch query");
    AlignedBuffer scratch(scratch_bytes);
    require_result(leo2_encode(codec, bytes, &original_pointers[0],
        &recovery_pointers[0], scratch.data, scratch.bytes),
        "guarded encode");
    require(original.guards_intact(),
        "guarded encode changed an input guard");
    require(original.shards() == wire,
        "guarded encode changed an input payload");
    require(parity.guards_intact(),
        "guarded encode changed an output guard");
    if (counts)
        counts->guard_checks += wire.size() + r;
    return parity.shards();
}

unsigned popcount32(uint32_t value)
{
    unsigned count = 0;
    while (value != 0)
    {
        count += value & 1u;
        value >>= 1;
    }
    return count;
}

std::vector<std::vector<unsigned> > coordinate_subsets(unsigned n, unsigned k)
{
    std::vector<std::vector<unsigned> > result;
    std::vector<unsigned> subset(k, 0);
    for (unsigned i = 0; i < k; ++i)
        subset[i] = i;
    for (;;)
    {
        result.push_back(subset);
        int i = static_cast<int>(k) - 1;
        while (i >= 0 && subset[static_cast<unsigned>(i)] ==
            n - k + static_cast<unsigned>(i))
            --i;
        if (i < 0)
            break;
        ++subset[static_cast<unsigned>(i)];
        for (unsigned j = static_cast<unsigned>(i) + 1; j < k; ++j)
            subset[j] = subset[j - 1] + 1;
    }
    return result;
}

void verify_generator_mds(
    const BinaryField& field,
    const Matrix& generator,
    unsigned k,
    Counts& counts)
{
    const std::vector<std::vector<unsigned> > subsets =
        coordinate_subsets(static_cast<unsigned>(generator.size()), k);
    for (size_t subset_i = 0; subset_i < subsets.size(); ++subset_i)
    {
        Matrix selected(k, std::vector<Element>(k, 0));
        for (unsigned row = 0; row < k; ++row)
            selected[row] = generator[subsets[subset_i][row]];
        Matrix inverse;
        require(leopard2_test::invert_matrix(field, selected, &inverse),
            "direct generator contains a singular public K-subset");
        for (unsigned basis = 0; basis < k; ++basis)
        {
            std::vector<Element> message(k, 0);
            message[basis] = 1;
            const std::vector<Element> encoded =
                leopard2_test::matrix_vector_multiply(field, generator, message);
            std::vector<Element> received(k, 0);
            for (unsigned row = 0; row < k; ++row)
                received[row] = encoded[subsets[subset_i][row]];
            require(leopard2_test::matrix_vector_multiply(
                field, inverse, received) == message,
                "direct generator basis recovery failed");
            ++counts.mds_basis_recoveries;
        }
    }
}

Shards pack_payloads(
    const leo2_codec* codec,
    const Shards& payload,
    size_t wire_bytes)
{
    Shards wire(payload.size(), Bytes(wire_bytes, 0xa5));
    for (size_t i = 0; i < payload.size(); ++i)
        require_result(leo2_pack_systematic_shard(codec, payload[i].size(),
            &payload[i][0], &wire[i][0], wire[i].size()), "systematic pack");
    return wire;
}

void recover_missing(
    const leo2_codec* codec,
    const Shards& payload,
    const Shards& wire,
    const Shards& parity,
    const std::vector<unsigned>& missing,
    Counts& counts)
{
    const unsigned k = static_cast<unsigned>(wire.size());
    const unsigned r = static_cast<unsigned>(parity.size());
    const size_t bytes = wire[0].size();
    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 1);
    for (size_t i = 0; i < missing.size(); ++i)
        original_present[missing[i]] = 0;

    leo2_decode_plan* plan = NULL;
    require_result(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), "decode plan create");
    require(plan != NULL, "decode plan create returned null");
    size_t scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(plan, bytes, &scratch_bytes),
        "decode scratch query");
    AlignedBuffer scratch(scratch_bytes);

    Shards restored(k, Bytes(bytes, 0xcc));
    std::vector<const void*> originals(k, NULL);
    std::vector<const void*> recoveries = const_pointers(parity);
    std::vector<void*> outputs(k, NULL);
    for (unsigned i = 0; i < k; ++i)
    {
        if (original_present[i])
            originals[i] = &wire[i][0];
        else
            outputs[i] = &restored[i][0];
    }
    require_result(leo2_decode_plan_execute(plan, bytes, &originals[0],
        &recoveries[0], &outputs[0], scratch.data, scratch.bytes),
        "decode execute");

    for (size_t i = 0; i < missing.size(); ++i)
    {
        const unsigned index = missing[i];
        require(restored[index] == wire[index],
            "restored physical systematic shard differs");
        Bytes unpacked(payload[index].size(), 0);
        require_result(leo2_unpack_systematic_shard(codec,
            payload[index].size(), &restored[index][0], restored[index].size(),
            &unpacked[0]), "systematic unpack");
        require(unpacked == payload[index], "restored payload differs");
        ++counts.restored_shards;
    }

    Shards complete_wire = wire;
    for (size_t i = 0; i < missing.size(); ++i)
        complete_wire[missing[i]] = restored[missing[i]];
    require(encode(codec, complete_wire) == parity,
        "parity rebuild after message recovery differs");
    leo2_decode_plan_destroy(plan);
}

void recover_erasure_pattern_guarded(
    const leo2_codec* codec,
    const Shards& payload,
    const Shards& wire,
    const Shards& parity,
    uint32_t loss_mask,
    Counts* counts)
{
    const unsigned k = static_cast<unsigned>(wire.size());
    const unsigned r = static_cast<unsigned>(parity.size());
    const size_t bytes = wire[0].size();
    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 1);
    unsigned missing_originals = 0;
    for (unsigned i = 0; i < k; ++i)
        if ((loss_mask & (UINT32_C(1) << i)) != 0)
        {
            original_present[i] = 0;
            ++missing_originals;
        }
    for (unsigned i = 0; i < r; ++i)
        if ((loss_mask & (UINT32_C(1) << (k + i))) != 0)
            recovery_present[i] = 0;

    leo2_decode_plan* plan = NULL;
    require_result(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), "guarded erasure plan create");
    require(plan != NULL, "guarded erasure plan create returned null");
    size_t scratch_bytes = 99;
    require_result(leo2_decode_plan_scratch_size(plan, bytes, &scratch_bytes),
        "guarded erasure scratch query");

    if (missing_originals == 0)
    {
        require(scratch_bytes == 0,
            "parity-only loss unexpectedly requires decode scratch");
        require_result(leo2_decode_plan_execute(
            plan, bytes, NULL, NULL, NULL, NULL, 0),
            "parity-only loss no-op decode");
        leo2_decode_plan_destroy(plan);
        if (counts)
            ++counts->exhaustive_patterns;
        return;
    }

    GuardedShards guarded_original(wire, 0xb5);
    GuardedShards guarded_parity(parity, 0x6d);
    GuardedShards restored(k, bytes, 0xe7, 0xcc);
    std::vector<const void*> original_input(k, NULL);
    std::vector<const void*> parity_input(r, NULL);
    std::vector<void*> restored_output(k, NULL);
    for (unsigned i = 0; i < k; ++i)
    {
        if (original_present[i])
            original_input[i] = guarded_original.data(i);
        else
            restored_output[i] = restored.mutable_data(i);
    }
    for (unsigned i = 0; i < r; ++i)
        if (recovery_present[i])
            parity_input[i] = guarded_parity.data(i);

    AlignedBuffer scratch(scratch_bytes);
    require_result(leo2_decode_plan_execute(plan, bytes,
        &original_input[0], &parity_input[0], &restored_output[0],
        scratch.data, scratch.bytes), "guarded erasure decode");
    require(guarded_original.guards_intact(),
        "guarded decode changed an original input guard");
    require(guarded_parity.guards_intact(),
        "guarded decode changed a parity input guard");
    require(guarded_original.shards() == wire,
        "guarded decode changed an original input payload");
    require(guarded_parity.shards() == parity,
        "guarded decode changed a parity input payload");
    require(restored.guards_intact(),
        "guarded decode changed a restored-output guard");

    Shards complete_wire = wire;
    for (unsigned i = 0; i < k; ++i)
        if (!original_present[i])
        {
            complete_wire[i] = restored.shard(i);
            require(complete_wire[i] == wire[i],
                "guarded erasure recovery differs from the physical systematic shard");
            Bytes unpacked(payload[i].size(), 0);
            require_result(leo2_unpack_systematic_shard(codec,
                payload[i].size(), &complete_wire[i][0], bytes,
                &unpacked[0]), "guarded recovered systematic unpack");
            require(unpacked == payload[i],
                "guarded erasure recovery differs from the payload");
        }
    require(encode_guarded(codec, complete_wire, NULL) == parity,
        "guarded parity rebuild differs after erasure recovery");

    if (counts)
    {
        ++counts->exhaustive_patterns;
        counts->restored_shards += missing_originals;
        counts->guard_checks += static_cast<uint64_t>(k) * 2 + r;
    }
    leo2_decode_plan_destroy(plan);
}

void test_requested_odd_matrix(
    leo2_context* context,
    const BinaryField& field,
    Counts* counts)
{
    const size_t payload_sizes[] = {
        1, 3, 17, 33, 65, 129, 257, 1025
    };
    struct ProfileCase
    {
        leo2_profile profile;
        unsigned k;
        unsigned r;
    };
    const ProfileCase cases[] = {
        { LEO2_PROFILE_LEGACY_HIGH_V1, 3, 2 },
        { LEO2_PROFILE_LOW_V1, 2, 3 }
    };

    for (size_t case_i = 0; case_i < sizeof(cases) / sizeof(cases[0]);
         ++case_i)
    {
        const ProfileCase& test = cases[case_i];
        CodecOwner codec;
        codec.codec = create_codec(context, test.k, test.r, test.profile,
            LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1, 0);
        const ProfileLayout layout = leopard2_test::make_profile_layout(
            test.profile == LEO2_PROFILE_LEGACY_HIGH_V1
                ? leopard2_test::kLegacyHigh : leopard2_test::kLow,
            test.k, test.r);
        const Matrix generator =
            leopard2_test::direct_systematic_generator(field, layout);
        if (counts)
            verify_generator_mds(field, generator, test.k, *counts);

        for (size_t size_i = 0;
             size_i < sizeof(payload_sizes) / sizeof(payload_sizes[0]);
             ++size_i)
        {
            const size_t payload_bytes = payload_sizes[size_i];
            uint64_t wire_bytes_u64 = 0;
            require_result(leo2_codec_wire_shard_bytes(codec.codec,
                payload_bytes, &wire_bytes_u64),
                "requested odd-matrix wire-size query");
            require(wire_bytes_u64 == payload_bytes + 1,
                "requested odd-matrix physical length is not B+1");
            const size_t wire_bytes = static_cast<size_t>(wire_bytes_u64);
            const Shards payload = deterministic_payloads(test.k,
                payload_bytes, UINT64_C(0x6f64646d61747269) ^
                (static_cast<uint64_t>(case_i) << 56) ^ payload_bytes);
            const Shards wire = pack_payloads(
                codec.codec, payload, wire_bytes);
            const Shards expected = direct_parity(
                field, generator, test.k, test.r, wire);
            const Shards actual = encode_guarded(codec.codec, wire, counts);
            require(actual == expected,
                "requested odd-matrix parity differs from the direct GF16 oracle");
            if (counts)
            {
                ++counts->payload_cases;
                counts->direct_symbols += static_cast<uint64_t>(test.r) *
                    (wire_bytes / 2);
            }

            const unsigned transmitted = test.k + test.r;
            const uint32_t mask_limit = UINT32_C(1) << transmitted;
            for (uint32_t loss_mask = 0; loss_mask < mask_limit; ++loss_mask)
                if (popcount32(loss_mask) <= test.r)
                    recover_erasure_pattern_guarded(codec.codec, payload,
                        wire, actual, loss_mask, counts);
        }
    }
}

void test_profile_payloads(
    leo2_context* context,
    const BinaryField& field,
    leo2_profile profile,
    unsigned k,
    unsigned r,
    const std::vector<unsigned>& missing,
    Counts& counts)
{
    CodecOwner native;
    CodecOwner padded;
    native.codec = create_codec(context, k, r, profile,
        LEO2_SHARD_LAYOUT_NATIVE_V1, LEO2_CODEC_FORCE_SPECIALIZED_DECODE);
    padded.codec = create_codec(context, k, r, profile,
        LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1,
        LEO2_CODEC_FORCE_SPECIALIZED_DECODE);

    const ProfileLayout oracle_layout = leopard2_test::make_profile_layout(
        profile == LEO2_PROFILE_LEGACY_HIGH_V1
            ? leopard2_test::kLegacyHigh : leopard2_test::kLow,
        k, r);
    const Matrix generator =
        leopard2_test::direct_systematic_generator(field, oracle_layout);
    verify_generator_mds(field, generator, k, counts);

    std::vector<size_t> payload_sizes;
    for (size_t bytes = 1; bytes <= 65; ++bytes)
        payload_sizes.push_back(bytes);
    payload_sizes.push_back(129);
    payload_sizes.push_back(257);
    payload_sizes.push_back(1023);
    payload_sizes.push_back(1024);
    payload_sizes.push_back(1025);

    for (size_t size_i = 0; size_i < payload_sizes.size(); ++size_i)
    {
        const size_t payload_bytes = payload_sizes[size_i];
        const leo2_codec* codec = (payload_bytes & 1u) != 0
            ? padded.codec : native.codec;
        uint64_t wire_bytes_u64 = 0;
        require_result(leo2_codec_wire_shard_bytes(
            codec, payload_bytes, &wire_bytes_u64), "wire-size query");
        const size_t wire_bytes = static_cast<size_t>(wire_bytes_u64);
        require(wire_bytes == payload_bytes + (payload_bytes & 1u),
            "wire-size query returned an incorrect physical size");

        const Shards payload = deterministic_payloads(k, payload_bytes,
            UINT64_C(0x7061646465640000) ^ (payload_bytes * 131u) ^
            (static_cast<unsigned>(profile) * 17u));
        const Shards wire = pack_payloads(codec, payload, wire_bytes);
        if ((payload_bytes & 1u) != 0)
            for (unsigned i = 0; i < k; ++i)
            {
                require(std::equal(payload[i].begin(), payload[i].end(), wire[i].begin()),
                    "padded systematic wire did not retain the payload prefix");
                require(wire[i][payload_bytes] == 0,
                    "padded systematic wire has a nonzero pad");
            }

        const Shards actual = encode(codec, wire);
        const Shards expected = direct_parity(field, generator, k, r, wire);
        require(actual == expected,
            "padded/native parity differs from independent direct GF16 generator");
        counts.direct_symbols += static_cast<uint64_t>(r) *
            (wire_bytes / 2);
        recover_missing(codec, payload, wire, actual, missing, counts);
        ++counts.payload_cases;
    }
}

void test_decode_dispatch_variants(
    leo2_context* context,
    const BinaryField& field,
    Counts& counts)
{
    const uint32_t flags[] = {
        0,
        LEO2_CODEC_FORCE_GENERIC_DECODE,
        LEO2_CODEC_FORCE_SPECIALIZED_DECODE
    };
    for (unsigned profile_i = 0; profile_i < 2; ++profile_i)
    {
        const leo2_profile profile = profile_i == 0
            ? LEO2_PROFILE_LEGACY_HIGH_V1 : LEO2_PROFILE_LOW_V1;
        const unsigned k = profile_i == 0 ? 17 : 5;
        const unsigned r = profile_i == 0 ? 5 : 11;
        const size_t payload_bytes = 33;
        const size_t wire_bytes = payload_bytes + 1;
        const ProfileLayout oracle_layout = leopard2_test::make_profile_layout(
            profile == LEO2_PROFILE_LEGACY_HIGH_V1
                ? leopard2_test::kLegacyHigh : leopard2_test::kLow,
            k, r);
        const Matrix generator =
            leopard2_test::direct_systematic_generator(field, oracle_layout);
        const Shards payload = deterministic_payloads(k, payload_bytes,
            UINT64_C(0x6469737061746368) ^ static_cast<unsigned>(profile));
        std::vector<unsigned> missing;
        for (unsigned i = 0; i < 5; ++i)
            missing.push_back(i);

        for (size_t flag_i = 0; flag_i < sizeof(flags) / sizeof(flags[0]);
             ++flag_i)
        {
            CodecOwner codec;
            codec.codec = create_codec(context, k, r, profile,
                LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1, flags[flag_i]);
            const Shards wire = pack_payloads(codec.codec, payload, wire_bytes);
            const Shards parity = encode(codec.codec, wire);
            require(parity == direct_parity(field, generator, k, r, wire),
                "decode-dispatch parity differs from direct generator");
            recover_missing(codec.codec, payload, wire, parity, missing, counts);
            counts.direct_symbols += static_cast<uint64_t>(r) * (wire_bytes / 2);
            ++counts.payload_cases;
        }
    }
}

void test_default_direct_repair(
    leo2_context* context,
    const BinaryField& field,
    Counts& counts)
{
    for (unsigned profile_i = 0; profile_i < 2; ++profile_i)
    {
        const leo2_profile profile = profile_i == 0
            ? LEO2_PROFILE_LEGACY_HIGH_V1 : LEO2_PROFILE_LOW_V1;
        const unsigned k = profile_i == 0 ? 9 : 5;
        const unsigned r = profile_i == 0 ? 5 : 11;
        CodecOwner codec;
        codec.codec = create_codec(context, k, r, profile,
            LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1, 0);
        const ProfileLayout oracle_layout = leopard2_test::make_profile_layout(
            profile == LEO2_PROFILE_LEGACY_HIGH_V1
                ? leopard2_test::kLegacyHigh : leopard2_test::kLow,
            k, r);
        const Matrix generator =
            leopard2_test::direct_systematic_generator(field, oracle_layout);
        const Shards payload = deterministic_payloads(k, 33,
            UINT64_C(0x6469726563740000) ^ static_cast<unsigned>(profile));
        const Shards wire = pack_payloads(codec.codec, payload, 34);
        const Shards parity = encode(codec.codec, wire);
        require(parity == direct_parity(field, generator, k, r, wire),
            "direct-repair parity differs from direct generator");
        recover_missing(codec.codec, payload, wire, parity,
            std::vector<unsigned>{ 0, 1, 2, 3 }, counts);
        counts.direct_symbols += static_cast<uint64_t>(r) * 17;
        ++counts.payload_cases;
    }
}

void test_parity_pad_is_physical(
    leo2_context* context,
    const BinaryField& field,
    leo2_profile profile,
    unsigned k,
    unsigned r)
{
    CodecOwner codec;
    codec.codec = create_codec(context, k, r, profile,
        LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1, 0);
    const ProfileLayout oracle_layout = leopard2_test::make_profile_layout(
        profile == LEO2_PROFILE_LEGACY_HIGH_V1
            ? leopard2_test::kLegacyHigh : leopard2_test::kLow,
        k, r);
    const Matrix generator =
        leopard2_test::direct_systematic_generator(field, oracle_layout);

    unsigned selected_original = 0;
    unsigned selected_parity = 0;
    Element selected_value = 0;
    Element expected_value = 0;
    bool found = false;
    for (unsigned parity = 0; parity < r && !found; ++parity)
        for (unsigned original = 0; original < k && !found; ++original)
            for (unsigned value = 1; value <= 255; ++value)
            {
                const Element product = field.multiply(
                    generator[k + parity][original], static_cast<Element>(value));
                if ((product >> 8) != 0)
                {
                    selected_original = original;
                    selected_parity = parity;
                    selected_value = static_cast<Element>(value);
                    expected_value = product;
                    found = true;
                    break;
                }
            }
    require(found, "could not construct a nonzero parity-pad witness");

    Shards payload(k, Bytes(1, 0));
    payload[selected_original][0] = static_cast<uint8_t>(selected_value);
    const Shards wire = pack_payloads(codec.codec, payload, 2);
    const Shards parity = encode(codec.codec, wire);
    require(parity[selected_parity][0] == static_cast<uint8_t>(expected_value) &&
        parity[selected_parity][1] == static_cast<uint8_t>(expected_value >> 8),
        "parity-pad witness differs from direct GF16 multiplication");
    require(parity[selected_parity][1] != 0,
        "the required transmitted parity pad was unexpectedly zero");
}

void test_pack_unpack_contract(leo2_context* context)
{
    CodecOwner native;
    CodecOwner padded;
    native.codec = create_codec(context, 5, 3, LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_SHARD_LAYOUT_NATIVE_V1, 0);
    padded.codec = create_codec(context, 5, 3, LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1, 0);

    const size_t sizes[] = {
        1, 3, 17, 33, 63, 65, 129, 257, 1023, 1025
    };
    for (size_t size_i = 0; size_i < sizeof(sizes) / sizeof(sizes[0]); ++size_i)
    {
        const size_t payload_bytes = sizes[size_i];
        const size_t wire_bytes = payload_bytes + 1;
        const Bytes payload = deterministic_bytes(payload_bytes,
            UINT64_C(0x6f7665726c617000) ^ payload_bytes);

        Bytes guarded_wire(wire_bytes + 2, 0xd3);
        require_result(leo2_pack_systematic_shard(padded.codec, payload_bytes,
            &payload[0], &guarded_wire[1], wire_bytes), "guarded pack");
        require(guarded_wire.front() == 0xd3 && guarded_wire.back() == 0xd3,
            "pack wrote outside the physical wire range");
        require(std::equal(payload.begin(), payload.end(), guarded_wire.begin() + 1) &&
            guarded_wire[wire_bytes] == 0, "guarded pack contents differ");

        Bytes guarded_payload(payload_bytes + 2, 0x7b);
        require_result(leo2_unpack_systematic_shard(padded.codec, payload_bytes,
            &guarded_wire[1], wire_bytes, &guarded_payload[1]), "guarded unpack");
        require(guarded_payload.front() == 0x7b && guarded_payload.back() == 0x7b,
            "unpack wrote outside the payload range");
        require(std::equal(payload.begin(), payload.end(), guarded_payload.begin() + 1),
            "guarded unpack contents differ");

        Bytes overlap(wire_bytes + 1, 0x55);
        std::copy(payload.begin(), payload.end(), overlap.begin());
        require_result(leo2_pack_systematic_shard(padded.codec, payload_bytes,
            &overlap[0], &overlap[1], wire_bytes), "right-overlap pack");
        require(std::equal(payload.begin(), payload.end(), overlap.begin() + 1) &&
            overlap[wire_bytes] == 0, "right-overlap pack failed");
        require_result(leo2_unpack_systematic_shard(padded.codec, payload_bytes,
            &overlap[1], wire_bytes, &overlap[0]), "left-overlap unpack");
        require(std::equal(payload.begin(), payload.end(), overlap.begin()),
            "left-overlap unpack failed");

        Bytes left_pack(wire_bytes + 1, 0x66);
        std::copy(payload.begin(), payload.end(), left_pack.begin() + 1);
        require_result(leo2_pack_systematic_shard(padded.codec, payload_bytes,
            &left_pack[1], &left_pack[0], wire_bytes), "left-overlap pack");
        require(std::equal(payload.begin(), payload.end(), left_pack.begin()) &&
            left_pack[payload_bytes] == 0, "left-overlap pack failed");

        Bytes right_unpack(wire_bytes + 1, 0x77);
        std::copy(payload.begin(), payload.end(), right_unpack.begin());
        right_unpack[payload_bytes] = 0;
        require_result(leo2_unpack_systematic_shard(padded.codec, payload_bytes,
            &right_unpack[0], wire_bytes, &right_unpack[1]),
            "right-overlap unpack");
        require(std::equal(payload.begin(), payload.end(), right_unpack.begin() + 1),
            "right-overlap unpack failed");

        Bytes in_place(wire_bytes, 0);
        std::copy(payload.begin(), payload.end(), in_place.begin());
        require_result(leo2_pack_systematic_shard(padded.codec, payload_bytes,
            &in_place[0], &in_place[0], wire_bytes), "in-place pack");
        require_result(leo2_unpack_systematic_shard(padded.codec, payload_bytes,
            &in_place[0], wire_bytes, &in_place[0]), "in-place unpack");
        require(std::equal(payload.begin(), payload.end(), in_place.begin()),
            "in-place pack/unpack failed");
    }

    uint64_t wire_bytes = 77;
    require(leo2_codec_wire_shard_bytes(NULL, 1, &wire_bytes) ==
        LEO2_INVALID_ARGUMENT && wire_bytes == 0,
        "null-codec wire query did not fail cleanly");
    require(leo2_codec_wire_shard_bytes(padded.codec, 0, &wire_bytes) ==
        LEO2_INVALID_ARGUMENT && wire_bytes == 0,
        "zero payload wire query did not fail cleanly");
    require(leo2_codec_wire_shard_bytes(padded.codec, 2, &wire_bytes) ==
        LEO2_UNSUPPORTED && wire_bytes == 0,
        "padded-odd layout accepted an even payload");
    require(leo2_codec_wire_shard_bytes(native.codec, 1, &wire_bytes) ==
        LEO2_UNSUPPORTED && wire_bytes == 0,
        "native GF16 layout accepted an odd payload");
    require(leo2_codec_wire_shard_bytes(padded.codec,
        std::numeric_limits<uint64_t>::max(), &wire_bytes) ==
        LEO2_INVALID_ARGUMENT && wire_bytes == 0,
        "wire-size overflow was not rejected");
    require(leo2_codec_wire_shard_bytes(padded.codec, 1, NULL) ==
        LEO2_INVALID_ARGUMENT, "null wire-size output was accepted");

    const uint8_t payload = 0x42;
    uint8_t wire[2] = { 0, 0 };
    require(leo2_pack_systematic_shard(padded.codec, 1, &payload, wire, 1) ==
        LEO2_INVALID_ARGUMENT, "undersized pack was accepted");
    require(leo2_pack_systematic_shard(padded.codec, 1, NULL, wire, 2) ==
        LEO2_INVALID_ARGUMENT, "null pack payload was accepted");
    require(leo2_pack_systematic_shard(padded.codec, 1, &payload, NULL, 2) ==
        LEO2_INVALID_ARGUMENT, "null pack output was accepted");
    wire[0] = payload;
    wire[1] = 1;
    uint8_t unpacked = 0;
    require(leo2_unpack_systematic_shard(padded.codec, 1, wire, 2, &unpacked) ==
        LEO2_INVALID_ARGUMENT, "nonzero systematic pad was accepted");
    require(leo2_unpack_systematic_shard(padded.codec, 1, wire, 1, &unpacked) ==
        LEO2_INVALID_ARGUMENT, "undersized unpack was accepted");

    const uintptr_t invalid_address =
        std::numeric_limits<uintptr_t>::max() - 31u;
    void* const invalid_mutable = reinterpret_cast<void*>(invalid_address);
    const void* const invalid_const =
        reinterpret_cast<const void*>(invalid_address);
    Bytes span_payload(65, 0x42);
    Bytes span_wire(66, 0x5a);
    const Bytes untouched_wire = span_wire;
    require(leo2_pack_systematic_shard(padded.codec, span_payload.size(),
        invalid_const, &span_wire[0], span_wire.size()) ==
        LEO2_INVALID_ARGUMENT,
        "unrepresentable pack source span was accepted");
    require(span_wire == untouched_wire,
        "rejected pack source span modified its valid destination");
    require(leo2_pack_systematic_shard(padded.codec, span_payload.size(),
        &span_payload[0], invalid_mutable, span_wire.size()) ==
        LEO2_INVALID_ARGUMENT,
        "unrepresentable pack destination span was accepted");

    span_wire.assign(66, 0x24);
    span_wire.back() = 0;
    Bytes span_output(65, 0xa5);
    const Bytes untouched_output = span_output;
    require(leo2_unpack_systematic_shard(padded.codec, span_output.size(),
        invalid_const, span_wire.size(), &span_output[0]) ==
        LEO2_INVALID_ARGUMENT,
        "unrepresentable unpack source span was accepted");
    require(span_output == untouched_output,
        "rejected unpack source span modified its valid destination");
    require(leo2_unpack_systematic_shard(padded.codec, span_output.size(),
        &span_wire[0], span_wire.size(), invalid_mutable) ==
        LEO2_INVALID_ARGUMENT,
        "unrepresentable unpack destination span was accepted");

    Bytes native_payload(64, 0x39);
    Bytes native_wire(64, 0);
    Bytes native_output(64, 0);
    require_result(leo2_pack_systematic_shard(native.codec,
        native_payload.size(), &native_payload[0], &native_wire[0],
        native_wire.size()), "native representable-span pack control");
    require_result(leo2_unpack_systematic_shard(native.codec,
        native_output.size(), &native_wire[0], native_wire.size(),
        &native_output[0]), "native representable-span unpack control");
    require(native_output == native_payload,
        "native representable-span control did not round-trip");
}

void test_native_odd_rejection_matrix(leo2_context* context)
{
    const size_t odd_sizes[] = { 1, 3, 17, 33, 65, 129, 257, 1025 };
    const leo2_profile profiles[] = {
        LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_PROFILE_LOW_V1
    };
    for (size_t profile_i = 0;
         profile_i < sizeof(profiles) / sizeof(profiles[0]); ++profile_i)
    {
        CodecOwner native;
        CodecOwner padded;
        native.codec = create_codec(context, 3, 2, profiles[profile_i],
            LEO2_SHARD_LAYOUT_NATIVE_V1, 0);
        padded.codec = create_codec(context, 3, 2, profiles[profile_i],
            LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1, 0);
        std::vector<uint8_t> original_present(3, 1);
        std::vector<uint8_t> recovery_present(2, 1);
        original_present[0] = 0;
        leo2_decode_plan* loss_plan = NULL;
        leo2_decode_plan* padded_loss_plan = NULL;
        require_result(leo2_decode_plan_create(native.codec,
            &original_present[0], &recovery_present[0], &loss_plan),
            "native odd-rejection loss plan create");
        require_result(leo2_decode_plan_create(padded.codec,
            &original_present[0], &recovery_present[0], &padded_loss_plan),
            "padded odd-physical rejection loss plan create");

        for (size_t size_i = 0;
             size_i < sizeof(odd_sizes) / sizeof(odd_sizes[0]); ++size_i)
        {
            const size_t bytes = odd_sizes[size_i];
            uint64_t wire_bytes = 77;
            require(leo2_codec_wire_shard_bytes(native.codec, bytes,
                &wire_bytes) == LEO2_UNSUPPORTED && wire_bytes == 0,
                "native GF16 odd wire query was not rejected deterministically");
            size_t scratch_bytes = 77;
            require(leo2_encode_scratch_size(native.codec, bytes,
                &scratch_bytes) == LEO2_UNSUPPORTED && scratch_bytes == 0,
                "native GF16 odd encode query was not rejected deterministically");
            require(leo2_encode(native.codec, bytes, NULL, NULL, NULL, 0) ==
                LEO2_UNSUPPORTED,
                "native GF16 odd encode execution disagrees with its query");
            scratch_bytes = 77;
            require(leo2_decode_scratch_size(native.codec, bytes,
                &scratch_bytes) == LEO2_UNSUPPORTED && scratch_bytes == 0,
                "native GF16 odd one-shot decode query was not rejected deterministically");
            require(leo2_decode(native.codec, bytes,
                &original_present[0], &recovery_present[0],
                NULL, NULL, NULL, NULL, 0) == LEO2_UNSUPPORTED,
                "native GF16 odd one-shot decode disagrees with its query");
            scratch_bytes = 77;
            require(leo2_decode_plan_scratch_size(loss_plan, bytes,
                &scratch_bytes) == LEO2_UNSUPPORTED && scratch_bytes == 0,
                "native GF16 odd plan query was not rejected deterministically");
            require(leo2_decode_plan_execute(loss_plan, bytes,
                NULL, NULL, NULL, NULL, 0) == LEO2_UNSUPPORTED,
                "native GF16 odd plan execution disagrees with its query");

            wire_bytes = 0;
            require_result(leo2_codec_wire_shard_bytes(padded.codec, bytes,
                &wire_bytes), "padded odd payload wire query");
            require(wire_bytes == bytes + 1,
                "padded odd payload query did not return physical B+1");
            scratch_bytes = 77;
            require(leo2_encode_scratch_size(padded.codec, bytes,
                &scratch_bytes) == LEO2_UNSUPPORTED && scratch_bytes == 0,
                "padded GF16 accepted an odd physical encode length");
            require(leo2_encode(padded.codec, bytes,
                NULL, NULL, NULL, 0) == LEO2_UNSUPPORTED,
                "padded GF16 odd physical encode disagrees with its query");
            scratch_bytes = 77;
            require(leo2_decode_plan_scratch_size(padded_loss_plan, bytes,
                &scratch_bytes) == LEO2_UNSUPPORTED && scratch_bytes == 0,
                "padded GF16 accepted an odd physical decode length");
            require(leo2_decode_plan_execute(padded_loss_plan, bytes,
                NULL, NULL, NULL, NULL, 0) == LEO2_UNSUPPORTED,
                "padded GF16 odd physical decode disagrees with its query");
        }
        leo2_decode_plan_destroy(padded_loss_plan);
        leo2_decode_plan_destroy(loss_plan);

        original_present[0] = 1;
        std::fill(recovery_present.begin(), recovery_present.end(), 0);
        leo2_decode_plan* no_loss_plan = NULL;
        leo2_decode_plan* padded_no_loss_plan = NULL;
        require_result(leo2_decode_plan_create(native.codec,
            &original_present[0], &recovery_present[0], &no_loss_plan),
            "native odd-rejection no-loss plan create");
        require_result(leo2_decode_plan_create(padded.codec,
            &original_present[0], &recovery_present[0],
            &padded_no_loss_plan),
            "padded odd-physical no-loss plan create");
        for (size_t size_i = 0;
             size_i < sizeof(odd_sizes) / sizeof(odd_sizes[0]); ++size_i)
        {
            size_t scratch_bytes = 77;
            require_result(leo2_decode_plan_scratch_size(no_loss_plan,
                odd_sizes[size_i], &scratch_bytes),
                "native odd no-loss scratch query");
            require(scratch_bytes == 0,
                "native odd no-loss plan unexpectedly requires scratch");
            require_result(leo2_decode_plan_execute(no_loss_plan,
                odd_sizes[size_i], NULL, NULL, NULL, NULL, 0),
                "native odd no-loss execution");
            scratch_bytes = 77;
            require_result(leo2_decode_plan_scratch_size(padded_no_loss_plan,
                odd_sizes[size_i], &scratch_bytes),
                "padded odd-physical no-loss scratch query");
            require(scratch_bytes == 0,
                "padded odd-physical no-loss plan requires scratch");
            require_result(leo2_decode_plan_execute(padded_no_loss_plan,
                odd_sizes[size_i], NULL, NULL, NULL, NULL, 0),
                "padded odd-physical no-loss execution");
        }
        leo2_decode_plan_destroy(padded_no_loss_plan);
        leo2_decode_plan_destroy(no_loss_plan);
    }
}

void test_auto_uses_gf8_for_legal_odd_code(
    leo2_context* context,
    Counts& counts)
{
    const unsigned k = 5;
    const unsigned r = 3;
    CodecOwner codec;
    require_result(leo2_codec_create(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_AUTO, NULL, &codec.codec),
        "AUTO odd codec create");
    require(leo2_codec_field(codec.codec) == LEO2_FIELD_GF8,
        "AUTO did not select GF8 for a legal odd-byte code");
    require(leo2_codec_shard_layout(codec.codec) ==
            LEO2_SHARD_LAYOUT_NATIVE_V1,
        "AUTO changed shard layout for a legal odd-byte code");

    const size_t bytes = 17;
    uint64_t wire_bytes = 0;
    require_result(leo2_codec_wire_shard_bytes(
        codec.codec, bytes, &wire_bytes), "AUTO GF8 wire-size query");
    require(wire_bytes == bytes,
        "AUTO GF8 changed the native odd physical length");
    const Shards payload = deterministic_payloads(
        k, bytes, UINT64_C(0x6175746f67663831));
    const Shards wire = pack_payloads(codec.codec, payload, bytes);
    const Shards parity = encode_guarded(codec.codec, wire, &counts);
    const BinaryField field = leopard2_test::make_legacy_gf8();
    const ProfileLayout layout = leopard2_test::make_profile_layout(
        leopard2_test::kLegacyHigh, k, r);
    const Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);
    require(parity == direct_gf8_parity(field, generator, k, r, wire),
        "AUTO GF8 odd parity differs from the direct oracle");
    recover_missing(codec.codec, payload, wire, parity,
        std::vector<unsigned>{ 0, 2, 4 }, counts);
    ++counts.payload_cases;
    counts.direct_symbols += static_cast<uint64_t>(r) * bytes;
}

void test_layout_creation_and_core_enforcement(leo2_context* context)
{
    leo2_codec_options options;
    memset(&options, 0, sizeof(options));
    options.struct_size = sizeof(options);
    options.shard_layout = LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1;
    leo2_codec* codec = reinterpret_cast<leo2_codec*>(static_cast<uintptr_t>(1));
    require(leo2_codec_create(context, 5, 3, LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_AUTO, &options, &codec) == LEO2_INVALID_ARGUMENT && codec == NULL,
        "padded-odd layout was AUTO field-selected");
    codec = reinterpret_cast<leo2_codec*>(static_cast<uintptr_t>(1));
    require(leo2_codec_create(context, 5, 3, LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF8, &options, &codec) == LEO2_INVALID_ARGUMENT && codec == NULL,
        "padded-odd layout accepted GF8");
    options.shard_layout = static_cast<leo2_shard_layout>(99);
    codec = reinterpret_cast<leo2_codec*>(static_cast<uintptr_t>(1));
    require(leo2_codec_create(context, 5, 3, LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF16, &options, &codec) == LEO2_INVALID_ARGUMENT && codec == NULL,
        "unknown shard layout was accepted");

    memset(&options, 0, sizeof(options));
    options.struct_size = offsetof(leo2_codec_options, shard_layout);
    options.shard_layout = LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1;
    require_result(leo2_codec_create(context, 5, 3,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, &options, &codec),
        "version-1 codec options");
    require(leo2_codec_shard_layout(codec) == LEO2_SHARD_LAYOUT_NATIVE_V1,
        "version-1 codec options did not default to native layout");
    leo2_codec_destroy(codec);
    codec = NULL;

    options.struct_size = offsetof(leo2_codec_options, shard_layout) +
        sizeof(options.shard_layout) - 1;
    require(leo2_codec_create(context, 5, 3, LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF16, &options, &codec) == LEO2_INVALID_ARGUMENT,
        "partial shard-layout field was accepted");

    CodecOwner padded;
    padded.codec = create_codec(context, 5, 3, LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1, 0);
    size_t scratch_bytes = 0;
    require(leo2_encode_scratch_size(padded.codec, 1, &scratch_bytes) ==
        LEO2_UNSUPPORTED, "padded codec accepted odd physical bytes");

    Shards wire(5, Bytes(2, 0));
    wire[0][1] = 1;
    Shards parity(3, Bytes(2, 0));
    std::vector<const void*> originals = const_pointers(wire);
    std::vector<void*> recoveries = mutable_pointers(parity);
    require_result(leo2_encode_scratch_size(padded.codec, 2, &scratch_bytes),
        "padded enforcement scratch query");
    AlignedBuffer scratch(scratch_bytes);
    require(leo2_encode(padded.codec, 2, &originals[0], &recoveries[0],
        scratch.data, scratch.bytes) == LEO2_INVALID_ARGUMENT,
        "padded codec accepted a nonzero systematic pad");

    wire[0][1] = 0;
    parity = encode(padded.codec, wire);
    std::vector<uint8_t> original_present(5, 1);
    std::vector<uint8_t> recovery_present(3, 1);
    original_present[0] = 0;
    leo2_decode_plan* loss_plan = NULL;
    require_result(leo2_decode_plan_create(padded.codec, &original_present[0],
        &recovery_present[0], &loss_plan), "pad-enforcement loss plan");
    require_result(leo2_decode_plan_scratch_size(loss_plan, 2, &scratch_bytes),
        "pad-enforcement loss scratch query");
    AlignedBuffer decode_scratch(scratch_bytes);
    wire[1][1] = 1;
    originals = const_pointers(wire);
    originals[0] = NULL;
    std::vector<const void*> parity_inputs = const_pointers(parity);
    Shards restored(5, Bytes(2, 0));
    std::vector<void*> restored_outputs(5, NULL);
    restored_outputs[0] = &restored[0][0];
    require(leo2_decode_plan_execute(loss_plan, 2, &originals[0],
        &parity_inputs[0], &restored_outputs[0], decode_scratch.data,
        decode_scratch.bytes) == LEO2_INVALID_ARGUMENT,
        "loss decode accepted a nonzero surviving systematic pad");
    leo2_decode_plan_destroy(loss_plan);

    std::fill(original_present.begin(), original_present.end(), 1);
    std::fill(recovery_present.begin(), recovery_present.end(), 0);
    leo2_decode_plan* no_loss_plan = NULL;
    require_result(leo2_decode_plan_create(padded.codec, &original_present[0],
        &recovery_present[0], &no_loss_plan), "pad-enforcement no-loss plan");
    scratch_bytes = 99;
    require_result(leo2_decode_plan_scratch_size(no_loss_plan, 2, &scratch_bytes),
        "pad-enforcement no-loss scratch query");
    require(scratch_bytes == 0, "padded no-loss plan requires scratch");
    require_result(leo2_decode_plan_execute(no_loss_plan, 2, &originals[0],
        &parity_inputs[0], &restored_outputs[0], NULL, 0),
        "padded no-loss execute");
    leo2_decode_plan_destroy(no_loss_plan);
}

void test_native_even_identity(leo2_context* context)
{
    CodecOwner implicit_native;
    CodecOwner explicit_native;
    require_result(leo2_codec_create(context, 5, 3,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, NULL,
        &implicit_native.codec), "implicit native codec create");
    explicit_native.codec = create_codec(context, 5, 3,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_SHARD_LAYOUT_NATIVE_V1, 0);
    require(leo2_codec_shard_layout(implicit_native.codec) ==
        LEO2_SHARD_LAYOUT_NATIVE_V1, "null options did not select native layout");

    const Shards original = deterministic_payloads(5, 64,
        UINT64_C(0x6e61746976653634));
    require(encode(implicit_native.codec, original) ==
        encode(explicit_native.codec, original),
        "explicit native layout changed existing even-byte parity");
}

void test_requested_matrix_on_other_available_backends(
    leo2_backend already_tested,
    const BinaryField& field)
{
    const leo2_backend backends[] = {
        LEO2_BACKEND_SCALAR,
        LEO2_BACKEND_SSSE3,
        LEO2_BACKEND_AVX2,
        LEO2_BACKEND_AVX512
    };
    for (size_t backend_i = 0;
         backend_i < sizeof(backends) / sizeof(backends[0]); ++backend_i)
    {
        if (backends[backend_i] == already_tested)
            continue;
        leo2_context_options options;
        memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = backends[backend_i];
        options.thread_count = 1;
        leo2_context* context = NULL;
        const leo2_result result = leo2_context_create(&options, &context);
        require(result == LEO2_SUCCESS || result == LEO2_UNSUPPORTED,
            "explicit backend returned an unexpected qualification result");
        if (result == LEO2_UNSUPPORTED)
        {
            require(context == NULL,
                "unsupported explicit backend returned a context");
            continue;
        }
        require(context != NULL &&
                leo2_context_backend(context) == backends[backend_i],
            "explicit backend context reports the wrong backend");
        test_requested_odd_matrix(context, field, NULL);
        leo2_context_destroy(context);
    }
}

} // namespace

int main()
{
    try
    {
        leo2_context_options options;
        memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.thread_count = 1;
        leo2_context* context = NULL;
        require_result(leo2_context_create(&options, &context), "context create");
        require(context != NULL, "context create returned null");

        Counts counts;
        const BinaryField field = leopard2_test::make_legacy_gf16();
        test_layout_creation_and_core_enforcement(context);
        test_pack_unpack_contract(context);
        test_native_odd_rejection_matrix(context);
        test_native_even_identity(context);
        test_auto_uses_gf8_for_legal_odd_code(context, counts);
        test_requested_odd_matrix(context, field, &counts);
        test_profile_payloads(context, field,
            LEO2_PROFILE_LEGACY_HIGH_V1, 5, 3,
            std::vector<unsigned>{ 0, 2, 4 }, counts);
        test_profile_payloads(context, field,
            LEO2_PROFILE_LOW_V1, 3, 5,
            std::vector<unsigned>{ 0, 1, 2 }, counts);
        test_decode_dispatch_variants(context, field, counts);
        test_default_direct_repair(context, field, counts);
        test_parity_pad_is_physical(context, field,
            LEO2_PROFILE_LOW_V1, 2, 256);
        test_requested_matrix_on_other_available_backends(
            leo2_context_backend(context), field);
        leo2_context_destroy(context);

        std::cout << "GF16 padded-odd passed: payload_cases="
                  << counts.payload_cases
                  << " direct_symbols=" << counts.direct_symbols
                  << " restored_shards=" << counts.restored_shards
                  << " mds_basis_recoveries=" << counts.mds_basis_recoveries
                  << " exhaustive_patterns=" << counts.exhaustive_patterns
                  << " guard_checks=" << counts.guard_checks
                  << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "GF16 padded-odd failed: " << error.what() << std::endl;
        return 1;
    }
}
