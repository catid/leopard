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
    Bounded low-decoder acceptance gate.

    Unlike codec round-trip tests, this executable constructs every parity
    input and every expected recovery with the slow direct generator in
    direct_oracle.cpp.  Production encode is invoked only after recovery, as
    the explicit parity-rebuild operation whose output is checked against that
    independent codeword.
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
#include <string>
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

typedef std::vector<uint8_t> Bytes;
typedef std::vector<Bytes> Shards;

const size_t kGuardBytes = 32;

struct Counts
{
    uint64_t profile_variants;
    uint64_t maximum_erasure_patterns;
    uint64_t plan_executions;
    uint64_t one_shot_executions;
    uint64_t batch_executions;
    uint64_t concurrent_executions;
    uint64_t direct_interpolation_symbols;
    uint64_t recovered_shards;
    uint64_t parity_rebuilds;
    uint64_t rebuilt_parity_bytes;
    uint64_t immutable_input_checks;
    uint64_t guard_checks;
    uint64_t overlap_rejections;
    uint64_t allowed_input_aliases;
    uint64_t digest;

    Counts()
        : profile_variants(0)
        , maximum_erasure_patterns(0)
        , plan_executions(0)
        , one_shot_executions(0)
        , batch_executions(0)
        , concurrent_executions(0)
        , direct_interpolation_symbols(0)
        , recovered_shards(0)
        , parity_rebuilds(0)
        , rebuilt_parity_bytes(0)
        , immutable_input_checks(0)
        , guard_checks(0)
        , overlap_rejections(0)
        , allowed_input_aliases(0)
        , digest(UINT64_C(14695981039346656037))
    {
    }
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

void digest_bytes(Counts* counts, const uint8_t* data, size_t bytes)
{
    for (size_t i = 0; i < bytes; ++i)
    {
        counts->digest ^= data[i];
        counts->digest *= UINT64_C(1099511628211);
    }
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

class GuardedShard
{
public:
    GuardedShard(size_t bytes, uint8_t guard)
        : storage_(bytes + kGuardBytes * 2 + 1, guard)
        , offset_(kGuardBytes)
        , bytes_(bytes)
        , guard_(guard)
    {
        if ((reinterpret_cast<uintptr_t>(&storage_[offset_]) &
             (leo2_scratch_alignment() - 1)) == 0)
        {
            ++offset_;
        }
        require((reinterpret_cast<uintptr_t>(&storage_[offset_]) &
                 (leo2_scratch_alignment() - 1)) != 0,
            "guarded application shard unexpectedly has scratch alignment");
    }

    uint8_t* data() { return &storage_[offset_]; }
    const uint8_t* data() const { return &storage_[offset_]; }
    size_t size() const { return bytes_; }
    const Bytes& raw() const { return storage_; }

    void assign(const Bytes& value)
    {
        require(value.size() == bytes_, "guarded shard assignment size mismatch");
        if (bytes_ != 0)
            memcpy(data(), &value[0], bytes_);
    }

    void fill(uint8_t value)
    {
        if (bytes_ != 0)
            memset(data(), value, bytes_);
    }

    Bytes payload() const
    {
        return Bytes(data(), data() + bytes_);
    }

    bool guards_intact() const
    {
        for (size_t i = 0; i < offset_; ++i)
            if (storage_[i] != guard_)
                return false;
        for (size_t i = offset_ + bytes_; i < storage_.size(); ++i)
            if (storage_[i] != guard_)
                return false;
        return true;
    }

private:
    Bytes storage_;
    size_t offset_;
    size_t bytes_;
    uint8_t guard_;
};

class CodecOwner
{
public:
    CodecOwner() : codec(NULL) {}
    ~CodecOwner() { leo2_codec_destroy(codec); }
    leo2_codec* codec;

private:
    CodecOwner(const CodecOwner&);
    CodecOwner& operator=(const CodecOwner&);
};

class PlanOwner
{
public:
    PlanOwner() : plan(NULL) {}
    ~PlanOwner() { leo2_decode_plan_destroy(plan); }
    leo2_decode_plan* plan;

private:
    PlanOwner(const PlanOwner&);
    PlanOwner& operator=(const PlanOwner&);
};

uint64_t splitmix64(uint64_t* state)
{
    uint64_t value = (*state += UINT64_C(0x9e3779b97f4a7c15));
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

struct CaseSpec
{
    unsigned k;
    unsigned r;
    leo2_field field;
    size_t bytes;
    std::vector<unsigned> missing_originals;
    std::vector<unsigned> selected_parities;
    const char* label;
};

size_t symbol_count(const CaseSpec& spec)
{
    if (spec.field == LEO2_FIELD_GF8)
        return spec.bytes;
    require((spec.bytes & 1u) == 0, "GF16 acceptance case has odd bytes");
    return spec.bytes / 2;
}

Element load_symbol(const CaseSpec& spec, const uint8_t* shard, size_t symbol)
{
    if (spec.field == LEO2_FIELD_GF8)
        return shard[symbol];
    const size_t symbols = symbol_count(spec);
    const size_t tile_symbol = symbol & 31u;
    const size_t tile_first = symbol - tile_symbol;
    const size_t tile_symbols = std::min<size_t>(32, symbols - tile_first);
    const size_t tile_byte = tile_first * 2;
    return static_cast<Element>(shard[tile_byte + tile_symbol] |
        (static_cast<unsigned>(
            shard[tile_byte + tile_symbols + tile_symbol]) << 8));
}

void store_symbol(
    const CaseSpec& spec,
    uint8_t* shard,
    size_t symbol,
    Element value)
{
    if (spec.field == LEO2_FIELD_GF8)
    {
        shard[symbol] = static_cast<uint8_t>(value);
        return;
    }
    const size_t symbols = symbol_count(spec);
    const size_t tile_symbol = symbol & 31u;
    const size_t tile_first = symbol - tile_symbol;
    const size_t tile_symbols = std::min<size_t>(32, symbols - tile_first);
    const size_t tile_byte = tile_first * 2;
    shard[tile_byte + tile_symbol] = static_cast<uint8_t>(value);
    shard[tile_byte + tile_symbols + tile_symbol] =
        static_cast<uint8_t>(value >> 8);
}

std::vector<uint8_t> presence(unsigned count, const std::vector<unsigned>& present)
{
    std::vector<uint8_t> result(count, 0);
    for (size_t i = 0; i < present.size(); ++i)
    {
        require(present[i] < count, "presence index out of range");
        require(result[present[i]] == 0, "duplicate presence index");
        result[present[i]] = 1;
    }
    return result;
}

std::vector<unsigned> surviving_originals(const CaseSpec& spec)
{
    std::vector<uint8_t> missing = presence(spec.k, spec.missing_originals);
    std::vector<unsigned> result;
    for (unsigned i = 0; i < spec.k; ++i)
        if (!missing[i])
            result.push_back(i);
    return result;
}

class DirectOracle
{
public:
    DirectOracle(const BinaryField& field, const CaseSpec& spec)
        : field_(field)
        , spec_(spec)
        , layout_(leopard2_test::make_profile_layout(
              leopard2_test::kLow, spec.k, spec.r))
        , generator_(leopard2_test::direct_systematic_generator(field, layout_))
    {
        require(spec.r > spec.k, "low acceptance case does not have R > K");
        require(spec.missing_originals.size() == spec.selected_parities.size(),
            "low acceptance case does not provide exactly one parity per loss");
        require(layout_.parent_size <= field.order(),
            "low acceptance case exceeds its field");
        const std::vector<unsigned> surviving = surviving_originals(spec);
        received_rows_ = surviving;
        for (size_t i = 0; i < spec.selected_parities.size(); ++i)
        {
            require(spec.selected_parities[i] < spec.r,
                "selected parity index out of range");
            received_rows_.push_back(spec.k + spec.selected_parities[i]);
        }
        require(received_rows_.size() == spec.k,
            "direct interpolation does not have exactly K coordinates");

        Matrix selected(spec.k, std::vector<Element>(spec.k, 0));
        for (unsigned row = 0; row < spec.k; ++row)
            selected[row] = generator_[received_rows_[row]];
        require(leopard2_test::invert_matrix(field_, selected, &inverse_),
            "selected direct-interpolation coordinates are singular");

        const uint64_t transmitted_erasures =
            static_cast<uint64_t>(spec.missing_originals.size()) +
            (spec.r - spec.selected_parities.size());
        require(transmitted_erasures == spec.r,
            "acceptance case is not a maximum-erasure pattern");
    }

    unsigned parent_count() const { return layout_.parent_size; }

    Shards make_originals(uint64_t seed) const
    {
        Shards result(spec_.k, Bytes(spec_.bytes, 0));
        uint64_t state = seed;
        const size_t symbols = symbol_count(spec_);
        for (unsigned original = 0; original < spec_.k; ++original)
            for (size_t symbol = 0; symbol < symbols; ++symbol)
                store_symbol(spec_, &result[original][0], symbol,
                    static_cast<Element>(splitmix64(&state)));
        return result;
    }

    Shards encode(const Shards& originals) const
    {
        require(originals.size() == spec_.k, "oracle original count mismatch");
        Shards parity(spec_.r, Bytes(spec_.bytes, 0));
        const size_t symbols = symbol_count(spec_);
        for (size_t symbol = 0; symbol < symbols; ++symbol)
        {
            std::vector<Element> message(spec_.k, 0);
            for (unsigned original = 0; original < spec_.k; ++original)
                message[original] = load_symbol(
                    spec_, &originals[original][0], symbol);
            const std::vector<Element> codeword =
                leopard2_test::matrix_vector_multiply(
                    field_, generator_, message);
            for (unsigned recovery = 0; recovery < spec_.r; ++recovery)
                store_symbol(spec_, &parity[recovery][0], symbol,
                    codeword[spec_.k + recovery]);
        }
        return parity;
    }

    Shards interpolate(const Shards& originals, const Shards& parity) const
    {
        Shards recovered(spec_.k, Bytes(spec_.bytes, 0));
        const size_t symbols = symbol_count(spec_);
        for (size_t symbol = 0; symbol < symbols; ++symbol)
        {
            std::vector<Element> received(spec_.k, 0);
            for (unsigned row = 0; row < spec_.k; ++row)
            {
                const unsigned public_row = received_rows_[row];
                received[row] = public_row < spec_.k
                    ? load_symbol(spec_, &originals[public_row][0], symbol)
                    : load_symbol(spec_, &parity[public_row - spec_.k][0], symbol);
            }
            const std::vector<Element> message =
                leopard2_test::matrix_vector_multiply(
                    field_, inverse_, received);
            for (unsigned original = 0; original < spec_.k; ++original)
                store_symbol(spec_, &recovered[original][0], symbol,
                    message[original]);
        }
        require(recovered == originals,
            "direct interpolation did not recover the source codeword");
        return recovered;
    }

private:
    const BinaryField& field_;
    const CaseSpec& spec_;
    ProfileLayout layout_;
    Matrix generator_;
    std::vector<unsigned> received_rows_;
    Matrix inverse_;
};

class Stripe
{
public:
    Stripe(
        const CaseSpec& spec,
        const DirectOracle& oracle,
        uint64_t seed,
        size_t scratch_bytes)
        : spec_(spec)
        , source_(oracle.make_originals(seed))
        , parity_value_(oracle.encode(source_))
        , expected_(oracle.interpolate(source_, parity_value_))
        , original_present_(spec.k, 1)
        , recovery_present_(presence(spec.r, spec.selected_parities))
        , original_input_(spec.k, NULL)
        , recovery_input_(spec.r, NULL)
        , restored_output_(spec.k, NULL)
        , scratch_(new AlignedBuffer(scratch_bytes))
    {
        const std::vector<uint8_t> missing = presence(
            spec.k, spec.missing_originals);
        originals_.reserve(spec.k);
        restored_.reserve(spec.k);
        for (unsigned i = 0; i < spec.k; ++i)
        {
            originals_.push_back(GuardedShard(
                spec.bytes, static_cast<uint8_t>(0x41u + i * 13u)));
            originals_.back().assign(source_[i]);
            restored_.push_back(GuardedShard(
                spec.bytes, static_cast<uint8_t>(0xb3u + i * 7u)));
            restored_.back().fill(0xcc);
            original_present_[i] = missing[i] ? 0 : 1;
        }
        parity_.reserve(spec.r);
        for (unsigned i = 0; i < spec.r; ++i)
        {
            parity_.push_back(GuardedShard(
                spec.bytes, static_cast<uint8_t>(0x69u + i * 5u)));
            parity_.back().assign(parity_value_[i]);
        }
        refresh_pointers();
        snapshot_inputs();
        present_restored_snapshot_.reserve(spec.k);
        for (unsigned i = 0; i < spec.k; ++i)
            present_restored_snapshot_.push_back(restored_[i].raw());
    }

    void reset_outputs()
    {
        for (size_t i = 0; i < spec_.missing_originals.size(); ++i)
            restored_[spec_.missing_originals[i]].fill(0xcc);
        refresh_pointers();
    }

    const void* const* original_void() const
    {
        return original_input_.empty() ? NULL : &original_input_[0];
    }

    const void* const* recovery_void() const
    {
        return recovery_input_.empty() ? NULL : &recovery_input_[0];
    }

    void* const* restored_void()
    {
        return restored_output_.empty() ? NULL : &restored_output_[0];
    }

    const uint8_t* original_presence() const { return &original_present_[0]; }
    const uint8_t* recovery_presence() const { return &recovery_present_[0]; }
    void* scratch_data() { return scratch_->data(); }
    size_t scratch_size() const { return scratch_->size(); }

    const std::vector<void*>& restored_pointers() const
    {
        return restored_output_;
    }

    const std::vector<const void*>& recovery_pointers() const
    {
        return recovery_input_;
    }

    void verify_contents() const
    {
        for (unsigned i = 0; i < spec_.k; ++i)
        {
            require(originals_[i].guards_intact(),
                "decode modified an original guard");
            require(originals_[i].raw() == original_snapshot_[i],
                "decode modified an original input");
            require(restored_[i].guards_intact(),
                "decode modified a restored-output guard");
            if (original_present_[i])
            {
                require(restored_[i].raw() == present_restored_snapshot_[i],
                    "decode wrote an output for a present original");
            }
            else
            {
                require(restored_[i].payload() == expected_[i],
                    "low decoder differs from direct interpolation");
            }
        }
        for (unsigned i = 0; i < spec_.r; ++i)
        {
            require(parity_[i].guards_intact(),
                "decode modified a recovery-input guard");
            require(parity_[i].raw() == parity_snapshot_[i],
                "message-only decode wrote a parity input");
        }
    }

    void verify(Counts* counts, bool add_digest) const
    {
        verify_contents();
        for (unsigned i = 0; i < spec_.k; ++i)
        {
            counts->guard_checks += 2;
            ++counts->immutable_input_checks;
            if (!original_present_[i])
            {
                ++counts->recovered_shards;
                if (add_digest)
                    digest_bytes(counts, restored_[i].data(), spec_.bytes);
            }
        }
        for (unsigned i = 0; i < spec_.r; ++i)
        {
            ++counts->guard_checks;
            ++counts->immutable_input_checks;
        }
    }

    void verify_inputs_only() const
    {
        for (unsigned i = 0; i < spec_.k; ++i)
            require(originals_[i].raw() == original_snapshot_[i],
                "failed decode modified an original input");
        for (unsigned i = 0; i < spec_.r; ++i)
            require(parity_[i].raw() == parity_snapshot_[i],
                "failed decode modified a parity input");
    }

    void rebuild_parity(const leo2_codec* codec, Counts* counts)
    {
        std::vector<const void*> complete(spec_.k, NULL);
        for (unsigned i = 0; i < spec_.k; ++i)
            complete[i] = original_present_[i]
                ? static_cast<const void*>(originals_[i].data())
                : static_cast<const void*>(restored_[i].data());

        std::vector<GuardedShard> rebuilt;
        rebuilt.reserve(spec_.r);
        std::vector<void*> output(spec_.r, NULL);
        for (unsigned i = 0; i < spec_.r; ++i)
        {
            rebuilt.push_back(GuardedShard(
                spec_.bytes, static_cast<uint8_t>(0x27u + i * 3u)));
            rebuilt.back().fill(0x5d);
            output[i] = rebuilt.back().data();
        }
        size_t scratch_bytes = 0;
        require_result(leo2_encode_scratch_size(
            codec, spec_.bytes, &scratch_bytes), "parity-rebuild scratch query");
        AlignedBuffer scratch(scratch_bytes);
        require_result(leo2_encode(codec, spec_.bytes, &complete[0], &output[0],
            scratch.data(), scratch.size()), "explicit parity rebuild");
        for (unsigned i = 0; i < spec_.r; ++i)
        {
            require(rebuilt[i].guards_intact(),
                "explicit parity rebuild modified an output guard");
            require(rebuilt[i].payload() == parity_value_[i],
                "explicit parity rebuild differs from direct encode");
            digest_bytes(counts, rebuilt[i].data(), spec_.bytes);
            ++counts->guard_checks;
        }
        verify_inputs_only();
        ++counts->parity_rebuilds;
        counts->rebuilt_parity_bytes +=
            static_cast<uint64_t>(spec_.r) * spec_.bytes;
    }

private:
    void refresh_pointers()
    {
        for (unsigned i = 0; i < spec_.k; ++i)
        {
            original_input_[i] = original_present_[i]
                ? static_cast<const void*>(originals_[i].data()) : NULL;
            restored_output_[i] = original_present_[i]
                ? NULL : static_cast<void*>(restored_[i].data());
        }
        for (unsigned i = 0; i < spec_.r; ++i)
            recovery_input_[i] = recovery_present_[i]
                ? static_cast<const void*>(parity_[i].data()) : NULL;
    }

    void snapshot_inputs()
    {
        original_snapshot_.clear();
        parity_snapshot_.clear();
        for (unsigned i = 0; i < spec_.k; ++i)
            original_snapshot_.push_back(originals_[i].raw());
        for (unsigned i = 0; i < spec_.r; ++i)
            parity_snapshot_.push_back(parity_[i].raw());
    }

    const CaseSpec& spec_;
    Shards source_;
    Shards parity_value_;
    Shards expected_;
    std::vector<GuardedShard> originals_;
    std::vector<GuardedShard> parity_;
    std::vector<GuardedShard> restored_;
    std::vector<uint8_t> original_present_;
    std::vector<uint8_t> recovery_present_;
    std::vector<const void*> original_input_;
    std::vector<const void*> recovery_input_;
    std::vector<void*> restored_output_;
    std::unique_ptr<AlignedBuffer> scratch_;
    std::vector<Bytes> original_snapshot_;
    std::vector<Bytes> parity_snapshot_;
    std::vector<Bytes> present_restored_snapshot_;
};

leo2_codec* create_codec(
    leo2_context* context,
    const CaseSpec& spec,
    uint32_t flags)
{
    leo2_codec_options options;
    memset(&options, 0, sizeof(options));
    options.struct_size = sizeof(options);
    options.flags = flags;
    options.shard_layout = LEO2_SHARD_LAYOUT_NATIVE_V1;
    leo2_codec* codec = NULL;
    require_result(leo2_codec_create(context, spec.k, spec.r,
        LEO2_PROFILE_LOW_V1, spec.field, &options, &codec),
        std::string(spec.label) + " codec create");
    require(codec != NULL, "codec create returned null");
    require(leo2_codec_profile(codec) == LEO2_PROFILE_LOW_V1,
        "codec selected a non-low profile");
    require(leo2_codec_field(codec) == spec.field,
        "codec selected the wrong field");
    return codec;
}

void create_plan(
    const CaseSpec& spec,
    const leo2_codec* codec,
    PlanOwner* owner,
    std::vector<uint8_t>* original_present,
    std::vector<uint8_t>* recovery_present)
{
    const std::vector<uint8_t> missing = presence(
        spec.k, spec.missing_originals);
    original_present->assign(spec.k, 1);
    for (unsigned i = 0; i < spec.k; ++i)
        if (missing[i])
            (*original_present)[i] = 0;
    *recovery_present = presence(spec.r, spec.selected_parities);
    require_result(leo2_decode_plan_create(codec, &(*original_present)[0],
        &(*recovery_present)[0], &owner->plan),
        std::string(spec.label) + " plan create");
    require(leo2_decode_plan_missing_original_count(owner->plan) ==
            spec.missing_originals.size(),
        "plan reports the wrong missing-original count");
}

void test_overlap_contract(
    const CaseSpec& spec,
    const leo2_decode_plan* plan,
    Stripe* stripe,
    Counts* counts)
{
    require(spec.missing_originals.size() >= 2,
        "overlap test requires at least two missing originals");
    stripe->reset_outputs();
    std::vector<void*> aliased = stripe->restored_pointers();
    aliased[spec.missing_originals[1]] =
        aliased[spec.missing_originals[0]];
    require(leo2_decode_plan_execute(plan, spec.bytes,
        stripe->original_void(), stripe->recovery_void(), &aliased[0],
        stripe->scratch_data(), stripe->scratch_size()) == LEO2_OVERLAP,
        "overlapping restored outputs were not rejected");
    stripe->verify_inputs_only();
    ++counts->overlap_rejections;

    stripe->reset_outputs();
    aliased = stripe->restored_pointers();
    const std::vector<const void*>& recovery = stripe->recovery_pointers();
    const void* received_parity = NULL;
    for (unsigned i = 0; i < spec.r; ++i)
        if (recovery[i])
        {
            received_parity = recovery[i];
            break;
        }
    require(received_parity != NULL, "overlap test has no received parity");
    aliased[spec.missing_originals[0]] = const_cast<void*>(received_parity);
    require(leo2_decode_plan_execute(plan, spec.bytes,
        stripe->original_void(), stripe->recovery_void(), &aliased[0],
        stripe->scratch_data(), stripe->scratch_size()) == LEO2_OVERLAP,
        "restored-output/input overlap was not rejected");
    stripe->verify_inputs_only();
    ++counts->overlap_rejections;
}

void run_profile_variant(
    leo2_context* context,
    const CaseSpec& spec,
    const DirectOracle& oracle,
    uint32_t flags,
    Counts* counts)
{
    CodecOwner codec;
    codec.codec = create_codec(context, spec, flags);
    require(leo2_codec_parent_count(codec.codec) == oracle.parent_count(),
        "codec parent differs from the direct profile");
    PlanOwner plan;
    std::vector<uint8_t> original_present;
    std::vector<uint8_t> recovery_present;
    create_plan(spec, codec.codec, &plan, &original_present, &recovery_present);
    size_t scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(
        plan.plan, spec.bytes, &scratch_bytes), "plan scratch query");
    require(scratch_bytes != 0, "lossy low plan unexpectedly has zero scratch");

    for (unsigned stripe_index = 0; stripe_index < 2; ++stripe_index)
    {
        Stripe stripe(spec, oracle,
            UINT64_C(0x4c4f574445434f44) ^
                (static_cast<uint64_t>(flags) << 48) ^ stripe_index,
            scratch_bytes);
        for (unsigned repeat = 0; repeat < 2; ++repeat)
        {
            stripe.reset_outputs();
            require_result(leo2_decode_plan_execute(plan.plan, spec.bytes,
                stripe.original_void(), stripe.recovery_void(),
                stripe.restored_void(), stripe.scratch_data(),
                stripe.scratch_size()), "reusable low-plan execute");
            stripe.verify(counts, true);
            ++counts->plan_executions;
        }
        if (stripe_index == 0)
        {
            stripe.rebuild_parity(codec.codec, counts);
            test_overlap_contract(spec, plan.plan, &stripe, counts);
        }
    }

    Stripe one_shot(spec, oracle,
        UINT64_C(0x6f6e6573686f7400) ^ flags, scratch_bytes);
    size_t one_shot_scratch_bytes = 0;
    require_result(leo2_decode_scratch_size(
        codec.codec, spec.bytes, &one_shot_scratch_bytes),
        "one-shot scratch query");
    AlignedBuffer one_shot_scratch(one_shot_scratch_bytes);
    one_shot.reset_outputs();
    require_result(leo2_decode(codec.codec, spec.bytes,
        one_shot.original_presence(), one_shot.recovery_presence(),
        one_shot.original_void(), one_shot.recovery_void(),
        one_shot.restored_void(), one_shot_scratch.data(),
        one_shot_scratch.size()), "one-shot low decode");
    one_shot.verify(counts, true);
    ++counts->one_shot_executions;

    ++counts->profile_variants;
    ++counts->maximum_erasure_patterns;
    counts->direct_interpolation_symbols +=
        static_cast<uint64_t>(symbol_count(spec)) * spec.k * 3;
}

void test_allowed_input_alias(leo2_context* context, Counts* counts)
{
    const CaseSpec spec = {
        1, 7, LEO2_FIELD_GF8, 17,
        std::vector<unsigned>(1, 0),
        std::vector<unsigned>{ 0, 1, 2, 3, 4, 5, 6 },
        "K=1 input-alias"
    };
    CodecOwner codec;
    codec.codec = create_codec(context, spec, 0);
    std::vector<uint8_t> original_present(1, 0);
    std::vector<uint8_t> recovery_present(7, 1);
    PlanOwner plan;
    require_result(leo2_decode_plan_create(codec.codec, &original_present[0],
        &recovery_present[0], &plan.plan), "input-alias plan create");
    size_t scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(plan.plan, spec.bytes,
        &scratch_bytes), "input-alias scratch query");
    AlignedBuffer scratch(scratch_bytes);
    GuardedShard shared(spec.bytes, 0x73);
    Bytes value(spec.bytes, 0);
    for (size_t i = 0; i < value.size(); ++i)
        value[i] = static_cast<uint8_t>(i * 19u + 7u);
    shared.assign(value);
    const Bytes snapshot = shared.raw();
    std::vector<const void*> originals(1, NULL);
    std::vector<const void*> recoveries(7, shared.data());
    GuardedShard restored(spec.bytes, 0xa9);
    restored.fill(0xcc);
    std::vector<void*> output(1, restored.data());
    require_result(leo2_decode_plan_execute(plan.plan, spec.bytes,
        &originals[0], &recoveries[0], &output[0], scratch.data(),
        scratch.size()), "allowed input-alias decode");
    require(restored.payload() == value,
        "K=1 aliased parity inputs restored the wrong message");
    require(shared.raw() == snapshot && shared.guards_intact() &&
            restored.guards_intact(),
        "allowed input-alias decode modified a guard or input");
    digest_bytes(counts, restored.data(), restored.size());
    ++counts->allowed_input_aliases;
    counts->guard_checks += 2;
    ++counts->immutable_input_checks;
    ++counts->recovered_shards;
}

void test_concurrent_and_batch_plan(
    leo2_context* context,
    const BinaryField& field,
    Counts* counts)
{
    const CaseSpec spec = {
        7, 17, LEO2_FIELD_GF8, 33,
        std::vector<unsigned>{ 0, 2, 3, 5, 6 },
        std::vector<unsigned>{ 0, 3, 7, 12, 16 },
        "concurrent low plan"
    };
    const DirectOracle oracle(field, spec);
    CodecOwner codec;
    codec.codec = create_codec(
        context, spec, LEO2_CODEC_FORCE_SPECIALIZED_DECODE);
    PlanOwner plan;
    std::vector<uint8_t> original_present;
    std::vector<uint8_t> recovery_present;
    create_plan(spec, codec.codec, &plan, &original_present, &recovery_present);
    size_t scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(plan.plan, spec.bytes,
        &scratch_bytes), "concurrent-plan scratch query");

    const unsigned thread_count = 8;
    const unsigned repetitions = 8;
    std::vector<std::unique_ptr<Stripe> > stripes;
    for (unsigned i = 0; i < thread_count; ++i)
        stripes.push_back(std::unique_ptr<Stripe>(new Stripe(
            spec, oracle, UINT64_C(0x7468726561640000) ^ i, scratch_bytes)));
    std::vector<std::string> errors(thread_count);
    std::atomic<unsigned> ready(0);
    std::atomic<bool> start(false);
    std::vector<std::thread> workers;
    for (unsigned thread = 0; thread < thread_count; ++thread)
    {
        workers.push_back(std::thread([&, thread]() {
            ready.fetch_add(1, std::memory_order_release);
            while (!start.load(std::memory_order_acquire))
                std::this_thread::yield();
            try
            {
                for (unsigned repeat = 0; repeat < repetitions; ++repeat)
                {
                    stripes[thread]->reset_outputs();
                    require_result(leo2_decode_plan_execute(
                        plan.plan, spec.bytes,
                        stripes[thread]->original_void(),
                        stripes[thread]->recovery_void(),
                        stripes[thread]->restored_void(),
                        stripes[thread]->scratch_data(),
                        stripes[thread]->scratch_size()),
                        "concurrent immutable-plan execute");
                    stripes[thread]->verify_contents();
                }
            }
            catch (const std::exception& error)
            {
                errors[thread] = error.what();
            }
        }));
    }
    while (ready.load(std::memory_order_acquire) != thread_count)
        std::this_thread::yield();
    start.store(true, std::memory_order_release);
    for (size_t i = 0; i < workers.size(); ++i)
        workers[i].join();
    for (unsigned i = 0; i < thread_count; ++i)
    {
        require(errors[i].empty(), errors[i]);
        stripes[i]->verify(counts, true);
    }
    counts->concurrent_executions +=
        static_cast<uint64_t>(thread_count) * repetitions;

    const size_t batch_size = 8;
    stripes.clear();
    std::vector<leo2_decode_batch_item> items(batch_size);
    for (size_t i = 0; i < batch_size; ++i)
    {
        stripes.push_back(std::unique_ptr<Stripe>(new Stripe(
            spec, oracle, UINT64_C(0x6261746368000000) ^ i, scratch_bytes)));
        stripes.back()->reset_outputs();
        leo2_decode_batch_item item;
        item.shard_bytes = spec.bytes;
        item.original = stripes.back()->original_void();
        item.recovery = stripes.back()->recovery_void();
        item.restored_original = stripes.back()->restored_void();
        item.scratch = stripes.back()->scratch_data();
        item.scratch_bytes = stripes.back()->scratch_size();
        items[i] = item;
    }
    require_result(leo2_decode_plan_execute_batch(
        plan.plan, &items[0], items.size()), "low-plan batch execute");
    for (size_t i = 0; i < stripes.size(); ++i)
        stripes[i]->verify(counts, true);
    counts->batch_executions += batch_size;
    counts->direct_interpolation_symbols +=
        static_cast<uint64_t>(symbol_count(spec)) * spec.k *
        (thread_count + batch_size);
}

} // namespace

int main()
{
    try
    {
        const BinaryField gf8 = leopard2_test::make_legacy_gf8();
        const BinaryField gf16 = leopard2_test::make_legacy_gf16();
        const CaseSpec cases[] = {
            {
                31, 65, LEO2_FIELD_GF8, 17,
                std::vector<unsigned>{
                    0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                    10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                    20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30
                },
                std::vector<unsigned>{
                    0, 2, 4, 6, 8, 10, 12, 14, 16, 18,
                    20, 22, 24, 26, 28, 30, 32, 34, 36, 38,
                    40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 64
                },
                "large GF8 maximum erasures"
            },
            {
                33, 65, LEO2_FIELD_GF16, 66,
                std::vector<unsigned>{
                    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                    11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                    21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
                    31, 32
                },
                std::vector<unsigned>{
                    0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
                    22, 24, 26, 28, 30, 32, 34, 36, 38, 40,
                    42, 44, 46, 48, 50, 52, 54, 56, 58, 60,
                    62, 64
                },
                "large GF16 maximum erasures"
            }
        };

        leo2_context_options options;
        memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AUTO;
        options.thread_count = 4;
        leo2_context* context = NULL;
        require_result(leo2_context_create(&options, &context),
            "low acceptance context create");
        require(context != NULL, "context create returned null");

        Counts counts;
        const uint32_t flags[] = {
            0,
            LEO2_CODEC_FORCE_SPECIALIZED_DECODE,
            LEO2_CODEC_FORCE_GENERIC_DECODE
        };
        for (size_t case_i = 0; case_i < sizeof(cases) / sizeof(cases[0]);
             ++case_i)
        {
            const BinaryField& field = cases[case_i].field == LEO2_FIELD_GF8
                ? gf8 : gf16;
            const DirectOracle oracle(field, cases[case_i]);
            for (size_t flag_i = 0; flag_i < sizeof(flags) / sizeof(flags[0]);
                 ++flag_i)
                run_profile_variant(
                    context, cases[case_i], oracle, flags[flag_i], &counts);
        }

        test_concurrent_and_batch_plan(context, gf8, &counts);
        test_allowed_input_alias(context, &counts);
        leo2_context_destroy(context);

        std::cout << "leopard2 low-decoder acceptance passed: variants="
                  << counts.profile_variants
                  << " max_patterns=" << counts.maximum_erasure_patterns
                  << " plan_exec=" << counts.plan_executions
                  << " one_shot=" << counts.one_shot_executions
                  << " batch=" << counts.batch_executions
                  << " concurrent=" << counts.concurrent_executions
                  << " direct_symbols=" << counts.direct_interpolation_symbols
                  << " recovered=" << counts.recovered_shards
                  << " parity_rebuilds=" << counts.parity_rebuilds
                  << " rebuilt_bytes=" << counts.rebuilt_parity_bytes
                  << " immutable_inputs=" << counts.immutable_input_checks
                  << " guards=" << counts.guard_checks
                  << " overlap_rejections=" << counts.overlap_rejections
                  << " allowed_input_aliases=" << counts.allowed_input_aliases
                  << " digest=0x" << std::hex << counts.digest << std::dec
                  << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "leopard2 low-decoder acceptance failed: "
                  << error.what() << std::endl;
        return 1;
    }
}
