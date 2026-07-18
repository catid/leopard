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
    ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
*/

/*
    Bounded acceptance gate for the legacy-high decoder.  Forced-specialized
    execution exercises the Algorithm 5 path, forced-generic execution remains
    an independent production fallback, and direct_oracle.cpp supplies a slow
    polynomial-basis generator and interpolation oracle for tractable profiles.
*/

#include "direct_oracle.h"
#include "leopard2.h"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <new>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

namespace {

using leopard2_test::BinaryField;
using leopard2_test::Element;
using leopard2_test::Matrix;
using leopard2_test::ProfileLayout;

static const uint8_t kGuard = 0xd7;

typedef std::vector<uint8_t> Bytes;

struct Counts
{
    uint64_t profiles;
    uint64_t specialized_executions;
    uint64_t generic_executions;
    uint64_t concurrent_executions;
    uint64_t restored_shards;
    uint64_t direct_parity_symbols;
    uint64_t direct_recovered_symbols;
    uint64_t parity_rebuilds;
    uint64_t untouched_parity_bytes;
    uint64_t guard_checks;
    uint64_t overlap_checks;

    Counts()
        : profiles(0)
        , specialized_executions(0)
        , generic_executions(0)
        , concurrent_executions(0)
        , restored_shards(0)
        , direct_parity_symbols(0)
        , direct_recovered_symbols(0)
        , parity_rebuilds(0)
        , untouched_parity_bytes(0)
        , guard_checks(0)
        , overlap_checks(0)
    {}
};

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void require_result(
    leo2_result actual,
    leo2_result expected,
    const std::string& operation)
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
        if (bytes_ == 0)
            return;
#if defined(_MSC_VER)
        data_ = _aligned_malloc(bytes_, leo2_scratch_alignment());
        if (!data_)
            throw std::bad_alloc();
#else
        if (posix_memalign(&data_, leo2_scratch_alignment(), bytes_) != 0)
            throw std::bad_alloc();
#endif
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
    GuardedShard(size_t bytes, unsigned residue)
        : storage_(bytes + 256, kGuard)
        , offset_(0)
        , bytes_(bytes)
    {
        const uintptr_t base = reinterpret_cast<uintptr_t>(&storage_[0]);
        for (size_t candidate = 64; candidate < 192; ++candidate)
        {
            if (((base + candidate) & 63u) == (residue & 63u))
            {
                offset_ = candidate;
                break;
            }
        }
        require(offset_ != 0, "could not construct guarded shard alignment");
    }

    uint8_t* data() { return &storage_[offset_]; }
    const uint8_t* data() const { return &storage_[offset_]; }
    size_t bytes() const { return bytes_; }

    void fill(uint8_t value)
    {
        memset(data(), value, bytes_);
    }

    bool payload_equals(const GuardedShard& other) const
    {
        return bytes_ == other.bytes_ && memcmp(data(), other.data(), bytes_) == 0;
    }

    bool payload_is(uint8_t value) const
    {
        for (size_t i = 0; i < bytes_; ++i)
            if (data()[i] != value)
                return false;
        return true;
    }

    Bytes snapshot() const
    {
        return Bytes(data(), data() + bytes_);
    }

    bool payload_equals(const Bytes& expected) const
    {
        return expected.size() == bytes_ &&
            memcmp(data(), &expected[0], bytes_) == 0;
    }

    bool guards_intact() const
    {
        for (size_t i = 0; i < offset_; ++i)
            if (storage_[i] != kGuard)
                return false;
        for (size_t i = offset_ + bytes_; i < storage_.size(); ++i)
            if (storage_[i] != kGuard)
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
    Shards(unsigned count, size_t bytes, unsigned residue_seed)
        : bytes_(bytes)
    {
        shards_.reserve(count);
        for (unsigned i = 0; i < count; ++i)
        {
            unsigned residue = (residue_seed + i * 13u) & 63u;
            if (residue == 0)
                residue = 1;
            shards_.push_back(GuardedShard(bytes, residue));
        }
    }

    unsigned size() const { return static_cast<unsigned>(shards_.size()); }
    size_t bytes() const { return bytes_; }
    GuardedShard& operator[](unsigned index) { return shards_[index]; }
    const GuardedShard& operator[](unsigned index) const { return shards_[index]; }

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

    void fill(uint8_t value)
    {
        for (unsigned i = 0; i < size(); ++i)
            shards_[i].fill(value);
    }

    bool guards_intact() const
    {
        for (unsigned i = 0; i < size(); ++i)
            if (!shards_[i].guards_intact())
                return false;
        return true;
    }

private:
    std::vector<GuardedShard> shards_;
    size_t bytes_;
};

uint64_t next_random(uint64_t* state)
{
    uint64_t value = *state;
    value ^= value >> 12;
    value ^= value << 25;
    value ^= value >> 27;
    *state = value;
    return value * UINT64_C(2685821657736338717);
}

void fill_originals(Shards& shards, uint64_t seed)
{
    uint64_t state = seed ? seed : UINT64_C(0x9e3779b97f4a7c15);
    for (unsigned shard = 0; shard < shards.size(); ++shard)
    {
        for (size_t offset = 0; offset < shards.bytes(); ++offset)
        {
            shards[shard].data()[offset] = static_cast<uint8_t>(
                next_random(&state) + shard * 29u + offset * 131u);
        }
    }
}

leo2_codec* create_codec(
    leo2_context* context,
    unsigned k,
    unsigned r,
    leo2_field field,
    uint32_t flags)
{
    leo2_codec_options options;
    memset(&options, 0, sizeof(options));
    options.struct_size = sizeof(options);
    options.flags = flags;
    leo2_codec* codec = NULL;
    require_success(leo2_codec_create(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, field, &options, &codec), "codec create");
    require(codec != NULL, "codec create returned null");
    require(leo2_codec_profile(codec) == LEO2_PROFILE_LEGACY_HIGH_V1,
        "codec selected a non-high profile");
    require(leo2_codec_field(codec) == field, "codec selected the wrong field");
    return codec;
}

void encode(const leo2_codec* codec, const Shards& original, Shards& recovery)
{
    std::vector<const void*> original_pointers = original.const_pointers();
    std::vector<void*> recovery_pointers = recovery.mutable_pointers();
    size_t scratch_bytes = 0;
    require_success(leo2_encode_scratch_size(codec, original.bytes(), &scratch_bytes),
        "encode scratch query");
    AlignedBuffer scratch(scratch_bytes);
    require_success(leo2_encode(codec, original.bytes(), &original_pointers[0],
        &recovery_pointers[0], scratch.data(), scratch.bytes()), "encode");
}

Element shard_symbol(
    const GuardedShard& shard,
    leo2_field field,
    size_t symbol)
{
    if (field == LEO2_FIELD_GF8)
        return shard.data()[symbol];

    const size_t bytes = shard.bytes();
    const size_t complete_bytes = bytes & ~static_cast<size_t>(63u);
    const size_t complete_symbols = complete_bytes / 2;
    size_t low = 0;
    size_t high = 0;
    if (symbol < complete_symbols)
    {
        const size_t tile = symbol / 32;
        const size_t lane = symbol % 32;
        low = tile * 64 + lane;
        high = low + 32;
    }
    else
    {
        const size_t tail = symbol - complete_symbols;
        const size_t tail_symbols = (bytes - complete_bytes) / 2;
        low = complete_bytes + tail;
        high = complete_bytes + tail_symbols + tail;
    }
    require(high < bytes, "GF16 compact symbol index is out of range");
    return static_cast<Element>(shard.data()[low] |
        (static_cast<unsigned>(shard.data()[high]) << 8));
}

struct DirectOracle
{
    Matrix generator;
    Matrix inverse_received;
    std::vector<unsigned> received_rows;
};

DirectOracle build_direct_oracle(
    const BinaryField& field,
    unsigned k,
    unsigned r,
    const std::vector<uint8_t>& original_present,
    const std::vector<uint8_t>& recovery_present)
{
    const ProfileLayout layout = leopard2_test::make_profile_layout(
        leopard2_test::kLegacyHigh, k, r);
    DirectOracle oracle;
    oracle.generator = leopard2_test::direct_systematic_generator(field, layout);
    for (unsigned i = 0; i < k; ++i)
        if (original_present[i])
            oracle.received_rows.push_back(i);
    for (unsigned i = 0;
         i < r && oracle.received_rows.size() < static_cast<size_t>(k);
         ++i)
    {
        if (recovery_present[i])
            oracle.received_rows.push_back(k + i);
    }
    require(oracle.received_rows.size() == static_cast<size_t>(k),
        "direct oracle did not select K received rows");

    Matrix selected(k, std::vector<Element>(k, 0));
    for (unsigned row = 0; row < k; ++row)
        selected[row] = oracle.generator[oracle.received_rows[row]];
    require(leopard2_test::invert_matrix(field, selected, &oracle.inverse_received),
        "direct oracle received matrix is singular");
    return oracle;
}

void verify_direct_parity(
    const BinaryField& field,
    leo2_field public_field,
    const DirectOracle& oracle,
    const Shards& original,
    const Shards& recovery,
    Counts* counts)
{
    const size_t symbols = public_field == LEO2_FIELD_GF8
        ? original.bytes() : original.bytes() / 2;
    std::vector<Element> message(original.size(), 0);
    for (size_t symbol = 0; symbol < symbols; ++symbol)
    {
        for (unsigned i = 0; i < original.size(); ++i)
            message[i] = shard_symbol(original[i], public_field, symbol);
        const std::vector<Element> codeword = leopard2_test::matrix_vector_multiply(
            field, oracle.generator, message);
        for (unsigned i = 0; i < recovery.size(); ++i)
        {
            require(shard_symbol(recovery[i], public_field, symbol) ==
                    codeword[original.size() + i],
                "encoded parity differs from the independent direct generator");
            ++counts->direct_parity_symbols;
        }
    }
}

void verify_direct_recovery(
    const BinaryField& field,
    leo2_field public_field,
    const DirectOracle& oracle,
    const Shards& original,
    const Shards& recovery,
    const Shards& restored,
    const std::vector<uint8_t>& original_present,
    Counts* counts)
{
    const size_t symbols = public_field == LEO2_FIELD_GF8
        ? original.bytes() : original.bytes() / 2;
    std::vector<Element> received(oracle.received_rows.size(), 0);
    for (size_t symbol = 0; symbol < symbols; ++symbol)
    {
        for (size_t i = 0; i < oracle.received_rows.size(); ++i)
        {
            const unsigned row = oracle.received_rows[i];
            received[i] = row < original.size()
                ? shard_symbol(original[row], public_field, symbol)
                : shard_symbol(recovery[row - original.size()], public_field, symbol);
        }
        const std::vector<Element> message = leopard2_test::matrix_vector_multiply(
            field, oracle.inverse_received, received);
        for (unsigned i = 0; i < original.size(); ++i)
        {
            if (!original_present[i])
            {
                require(shard_symbol(restored[i], public_field, symbol) == message[i],
                    "specialized recovery differs from direct interpolation");
                ++counts->direct_recovered_symbols;
            }
        }
    }
}

struct Scenario
{
    unsigned k;
    unsigned r;
    leo2_field field;
    size_t bytes;
    const unsigned* missing_originals;
    size_t missing_original_count;
    const unsigned* missing_recovery;
    size_t missing_recovery_count;
    bool use_direct_oracle;
};

std::vector<Bytes> snapshot_shards(const Shards& shards)
{
    std::vector<Bytes> result(shards.size());
    for (unsigned i = 0; i < shards.size(); ++i)
        result[i] = shards[i].snapshot();
    return result;
}

void require_snapshots(
    const Shards& shards,
    const std::vector<Bytes>& expected,
    const std::string& message)
{
    require(shards.size() == expected.size(), "snapshot shard count differs");
    for (unsigned i = 0; i < shards.size(); ++i)
        require(shards[i].payload_equals(expected[i]), message);
}

void run_scenario(
    leo2_context* context,
    const Scenario& scenario,
    const BinaryField& gf8,
    const BinaryField& gf16,
    unsigned scenario_index,
    Counts* counts)
{
    leo2_codec* specialized_codec = create_codec(context, scenario.k, scenario.r,
        scenario.field, LEO2_CODEC_FORCE_SPECIALIZED_DECODE);
    leo2_codec* generic_codec = create_codec(context, scenario.k, scenario.r,
        scenario.field, LEO2_CODEC_FORCE_GENERIC_DECODE);

    Shards original(scenario.k, scenario.bytes, 1u + scenario_index * 5u);
    Shards recovery(scenario.r, scenario.bytes, 7u + scenario_index * 9u);
    fill_originals(original, UINT64_C(0x6a09e667f3bcc909) + scenario_index * 101u);
    encode(specialized_codec, original, recovery);

    std::vector<uint8_t> original_present(scenario.k, 1);
    std::vector<uint8_t> recovery_present(scenario.r, 1);
    for (size_t i = 0; i < scenario.missing_original_count; ++i)
        original_present[scenario.missing_originals[i]] = 0;
    for (size_t i = 0; i < scenario.missing_recovery_count; ++i)
        recovery_present[scenario.missing_recovery[i]] = 0;

    std::unique_ptr<DirectOracle> direct;
    const BinaryField& direct_field = scenario.field == LEO2_FIELD_GF8 ? gf8 : gf16;
    if (scenario.use_direct_oracle)
    {
        direct.reset(new DirectOracle(build_direct_oracle(direct_field,
            scenario.k, scenario.r, original_present, recovery_present)));
        verify_direct_parity(direct_field, scenario.field, *direct,
            original, recovery, counts);
    }

    leo2_decode_plan* specialized_plan = NULL;
    leo2_decode_plan* generic_plan = NULL;
    require_success(leo2_decode_plan_create(specialized_codec, &original_present[0],
        &recovery_present[0], &specialized_plan), "specialized plan create");
    require_success(leo2_decode_plan_create(generic_codec, &original_present[0],
        &recovery_present[0], &generic_plan), "generic plan create");
    require(leo2_decode_plan_missing_original_count(specialized_plan) ==
            scenario.missing_original_count,
        "specialized plan reports the wrong missing-original count");

    std::vector<const void*> original_inputs = original.const_pointers();
    std::vector<const void*> recovery_inputs = recovery.const_pointers();
    for (unsigned i = 0; i < scenario.k; ++i)
        if (!original_present[i])
            original_inputs[i] = NULL;
    for (unsigned i = 0; i < scenario.r; ++i)
        if (!recovery_present[i])
            recovery_inputs[i] = NULL;

    const std::vector<Bytes> original_before = snapshot_shards(original);
    const std::vector<Bytes> recovery_before = snapshot_shards(recovery);

    Shards specialized_restored(scenario.k, scenario.bytes, 19u + scenario_index);
    std::vector<void*> specialized_outputs(scenario.k, NULL);
    for (unsigned i = 0; i < scenario.k; ++i)
        if (!original_present[i])
            specialized_outputs[i] = specialized_restored[i].data();
    size_t specialized_scratch_bytes = 0;
    require_success(leo2_decode_plan_scratch_size(
        specialized_plan, scenario.bytes, &specialized_scratch_bytes),
        "specialized scratch query");
    AlignedBuffer specialized_scratch(specialized_scratch_bytes);
    for (unsigned repeat = 0; repeat < 2; ++repeat)
    {
        specialized_restored.fill(0xa5);
        require_success(leo2_decode_plan_execute(specialized_plan, scenario.bytes,
            &original_inputs[0], &recovery_inputs[0], &specialized_outputs[0],
            specialized_scratch.data(), specialized_scratch.bytes()),
            "specialized plan execute");
        ++counts->specialized_executions;
        for (unsigned i = 0; i < scenario.k; ++i)
        {
            if (original_present[i])
                require(specialized_restored[i].payload_is(0xa5),
                    "message-only decode touched an unrequested original output");
            else
            {
                require(specialized_restored[i].payload_equals(original[i]),
                    "specialized decoder restored the wrong original");
                ++counts->restored_shards;
            }
        }
        require_snapshots(original, original_before,
            "specialized decode modified an original input");
        require_snapshots(recovery, recovery_before,
            "message-only decode modified or rebuilt parity");
        counts->untouched_parity_bytes +=
            static_cast<uint64_t>(scenario.r) * scenario.bytes;
    }

    if (direct.get())
    {
        verify_direct_recovery(direct_field, scenario.field, *direct,
            original, recovery, specialized_restored, original_present, counts);
    }

    Shards generic_restored(scenario.k, scenario.bytes, 31u + scenario_index);
    generic_restored.fill(0x5a);
    std::vector<void*> generic_outputs(scenario.k, NULL);
    for (unsigned i = 0; i < scenario.k; ++i)
        if (!original_present[i])
            generic_outputs[i] = generic_restored[i].data();
    size_t generic_scratch_bytes = 0;
    require_success(leo2_decode_plan_scratch_size(
        generic_plan, scenario.bytes, &generic_scratch_bytes),
        "generic scratch query");
    AlignedBuffer generic_scratch(generic_scratch_bytes);
    require_success(leo2_decode_plan_execute(generic_plan, scenario.bytes,
        &original_inputs[0], &recovery_inputs[0], &generic_outputs[0],
        generic_scratch.data(), generic_scratch.bytes()), "generic plan execute");
    ++counts->generic_executions;
    for (unsigned i = 0; i < scenario.k; ++i)
    {
        if (original_present[i])
            require(generic_restored[i].payload_is(0x5a),
                "generic fallback touched an unrequested output");
        else
            require(generic_restored[i].payload_equals(specialized_restored[i]),
                "specialized and generic recovery differ");
    }

    Shards rebuilt_original(scenario.k, scenario.bytes, 41u + scenario_index);
    for (unsigned i = 0; i < scenario.k; ++i)
    {
        const GuardedShard& source = original_present[i]
            ? original[i] : specialized_restored[i];
        memcpy(rebuilt_original[i].data(), source.data(), scenario.bytes);
    }
    Shards rebuilt_recovery(scenario.r, scenario.bytes, 47u + scenario_index);
    encode(specialized_codec, rebuilt_original, rebuilt_recovery);
    for (unsigned i = 0; i < scenario.r; ++i)
        require(rebuilt_recovery[i].payload_equals(recovery[i]),
            "explicit parity rebuild differs from the original parity");
    ++counts->parity_rebuilds;

    unsigned present_parity = 0;
    while (present_parity < scenario.r && !recovery_present[present_parity])
        ++present_parity;
    require(present_parity < scenario.r, "alias test has no surviving parity");
    std::vector<void*> alias_outputs(scenario.k, NULL);
    for (unsigned i = 0; i < scenario.k; ++i)
        if (!original_present[i])
            alias_outputs[i] = specialized_restored[i].data();
    alias_outputs[scenario.missing_originals[0]] =
        const_cast<uint8_t*>(recovery[present_parity].data());
    require_result(leo2_decode_plan_execute(specialized_plan, scenario.bytes,
        &original_inputs[0], &recovery_inputs[0], &alias_outputs[0],
        specialized_scratch.data(), specialized_scratch.bytes()),
        LEO2_OVERLAP, "decode input/output alias");
    ++counts->overlap_checks;
    require_snapshots(recovery, recovery_before,
        "rejected alias call modified parity");

    require(original.guards_intact() && recovery.guards_intact() &&
            specialized_restored.guards_intact() && generic_restored.guards_intact() &&
            rebuilt_original.guards_intact() && rebuilt_recovery.guards_intact(),
        "scenario guard region was modified");
    counts->guard_checks += static_cast<uint64_t>(scenario.k) * 4u +
        static_cast<uint64_t>(scenario.r) * 2u;
    ++counts->profiles;

    leo2_decode_plan_destroy(generic_plan);
    leo2_decode_plan_destroy(specialized_plan);
    leo2_codec_destroy(generic_codec);
    leo2_codec_destroy(specialized_codec);
}

struct ConcurrentStripe
{
    std::unique_ptr<Shards> original;
    std::unique_ptr<Shards> recovery;
    std::unique_ptr<Shards> restored;
    std::vector<const void*> original_inputs;
    std::vector<const void*> recovery_inputs;
    std::vector<void*> outputs;
    std::unique_ptr<AlignedBuffer> scratch;
};

void test_concurrent_reuse(leo2_context* context, Counts* counts)
{
    const unsigned k = 65;
    const unsigned r = 32;
    const size_t bytes = 65;
    const unsigned missing[] = { 0, 16, 32, 48, 64 };
    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 1);
    for (size_t i = 0; i < sizeof(missing) / sizeof(missing[0]); ++i)
        original_present[missing[i]] = 0;
    recovery_present[30] = 0;
    recovery_present[31] = 0;

    leo2_codec* codec = create_codec(context, k, r, LEO2_FIELD_GF8,
        LEO2_CODEC_FORCE_SPECIALIZED_DECODE);
    leo2_decode_plan* plan = NULL;
    require_success(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), "concurrent plan create");
    size_t scratch_bytes = 0;
    require_success(leo2_decode_plan_scratch_size(plan, bytes, &scratch_bytes),
        "concurrent scratch query");

    const unsigned thread_count = 4;
    const unsigned repeats = 8;
    std::vector<std::unique_ptr<ConcurrentStripe> > stripes;
    for (unsigned thread_index = 0; thread_index < thread_count; ++thread_index)
    {
        std::unique_ptr<ConcurrentStripe> stripe(new ConcurrentStripe);
        stripe->original.reset(new Shards(k, bytes, 3u + thread_index));
        stripe->recovery.reset(new Shards(r, bytes, 17u + thread_index));
        stripe->restored.reset(new Shards(k, bytes, 29u + thread_index));
        fill_originals(*stripe->original,
            UINT64_C(0xbb67ae8584caa73b) + thread_index * 103u);
        encode(codec, *stripe->original, *stripe->recovery);
        stripe->restored->fill(0xcc);
        stripe->original_inputs = stripe->original->const_pointers();
        stripe->recovery_inputs = stripe->recovery->const_pointers();
        stripe->outputs.assign(k, NULL);
        for (unsigned i = 0; i < k; ++i)
        {
            if (!original_present[i])
            {
                stripe->original_inputs[i] = NULL;
                stripe->outputs[i] = (*stripe->restored)[i].data();
            }
        }
        stripe->recovery_inputs[30] = NULL;
        stripe->recovery_inputs[31] = NULL;
        stripe->scratch.reset(new AlignedBuffer(scratch_bytes));
        stripes.push_back(std::move(stripe));
    }

    std::atomic<unsigned> failures(0);
    std::vector<std::thread> threads;
    for (unsigned thread_index = 0; thread_index < thread_count; ++thread_index)
    {
        ConcurrentStripe* stripe = stripes[thread_index].get();
        threads.push_back(std::thread([plan, stripe, &failures]() {
            for (unsigned repeat = 0; repeat < repeats; ++repeat)
            {
                const leo2_result result = leo2_decode_plan_execute(plan, bytes,
                    &stripe->original_inputs[0], &stripe->recovery_inputs[0],
                    &stripe->outputs[0], stripe->scratch->data(),
                    stripe->scratch->bytes());
                if (result != LEO2_SUCCESS)
                    ++failures;
            }
        }));
    }
    for (size_t i = 0; i < threads.size(); ++i)
        threads[i].join();
    require(failures.load() == 0, "concurrent plan execution failed");

    for (unsigned thread_index = 0; thread_index < thread_count; ++thread_index)
    {
        const ConcurrentStripe& stripe = *stripes[thread_index];
        for (size_t i = 0; i < sizeof(missing) / sizeof(missing[0]); ++i)
        {
            require((*stripe.restored)[missing[i]].payload_equals(
                        (*stripe.original)[missing[i]]),
                "concurrent plan restored the wrong original");
        }
        require(stripe.original->guards_intact() &&
                stripe.recovery->guards_intact() &&
                stripe.restored->guards_intact(),
            "concurrent execution modified a guard region");
    }
    counts->concurrent_executions += thread_count * repeats;
    counts->guard_checks += static_cast<uint64_t>(thread_count) * (k + r + k);
    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
}

void verify_expected_backend(const leo2_context* context)
{
    const char* expected = std::getenv("LEO2_EXPECT_BACKEND");
    if (!expected || expected[0] == '\0')
        return;
    leo2_backend backend = LEO2_BACKEND_AUTO;
    if (std::strcmp(expected, "scalar") == 0)
        backend = LEO2_BACKEND_SCALAR;
    else if (std::strcmp(expected, "ssse3") == 0)
        backend = LEO2_BACKEND_SSSE3;
    else if (std::strcmp(expected, "avx2") == 0)
        backend = LEO2_BACKEND_AVX2;
    else if (std::strcmp(expected, "avx512") == 0)
        backend = LEO2_BACKEND_AVX512;
    else
        throw std::runtime_error("invalid LEO2_EXPECT_BACKEND value");
    require(leo2_context_backend(context) == backend,
        "forced build selected the wrong runtime backend");
}

} // namespace

int main()
{
    try
    {
        static const unsigned missing_a[] = { 0, 2, 4, 6, 8, 10, 12, 16 };
        static const unsigned missing_b[] = { 0, 16, 32, 48, 64 };
        static const unsigned missing_b_parity[] = { 29, 30, 31 };
        static const unsigned missing_c[] = { 0, 8, 16 };
        static const unsigned missing_c_parity[] = { 6, 7 };
        static const unsigned missing_d[] = { 0, 128, 256 };
        static const unsigned missing_d_parity[] = { 31, 32 };
        static const Scenario scenarios[] = {
            { 17, 8, LEO2_FIELD_GF8, 17, missing_a,
              sizeof(missing_a) / sizeof(missing_a[0]), NULL, 0, true },
            { 65, 32, LEO2_FIELD_GF8, 65, missing_b,
              sizeof(missing_b) / sizeof(missing_b[0]), missing_b_parity,
              sizeof(missing_b_parity) / sizeof(missing_b_parity[0]), false },
            { 17, 8, LEO2_FIELD_GF16, 34, missing_c,
              sizeof(missing_c) / sizeof(missing_c[0]), missing_c_parity,
              sizeof(missing_c_parity) / sizeof(missing_c_parity[0]), true },
            { 257, 33, LEO2_FIELD_GF16, 66, missing_d,
              sizeof(missing_d) / sizeof(missing_d[0]), missing_d_parity,
              sizeof(missing_d_parity) / sizeof(missing_d_parity[0]), false }
        };

        leo2_context_options options;
        memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AUTO;
        options.thread_count = 1;
        leo2_context* context = NULL;
        require_success(leo2_context_create(&options, &context), "context create");
        require(context != NULL, "context create returned null");
        verify_expected_backend(context);

        const BinaryField gf8 = leopard2_test::make_legacy_gf8();
        const BinaryField gf16 = leopard2_test::make_legacy_gf16();
        Counts counts;
        for (size_t i = 0; i < sizeof(scenarios) / sizeof(scenarios[0]); ++i)
            run_scenario(context, scenarios[i], gf8, gf16,
                static_cast<unsigned>(i), &counts);
        test_concurrent_reuse(context, &counts);
        leo2_context_destroy(context);

        std::cout << "high decoder acceptance passed: profiles=" << counts.profiles
                  << " specialized_executions=" << counts.specialized_executions
                  << " generic_executions=" << counts.generic_executions
                  << " concurrent_executions=" << counts.concurrent_executions
                  << " restored_shards=" << counts.restored_shards
                  << " direct_parity_symbols=" << counts.direct_parity_symbols
                  << " direct_recovered_symbols=" << counts.direct_recovered_symbols
                  << " parity_rebuilds=" << counts.parity_rebuilds
                  << " untouched_parity_bytes=" << counts.untouched_parity_bytes
                  << " guard_checks=" << counts.guard_checks
                  << " overlap_checks=" << counts.overlap_checks
                  << std::endl;
        return 0;
    }
    catch (const std::exception& exception)
    {
        std::cerr << "high decoder acceptance failed: "
                  << exception.what() << std::endl;
        return 1;
    }
}
