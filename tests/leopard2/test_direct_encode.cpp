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

#include "Leopard2Direct.h"
#include "Leopard2Backend.h"
#include "LeopardFF16.h"
#include "LeopardFF8.h"
#include "direct_oracle.h"
#include "leopard.h"
#include "leopard2.h"
#include "allocation_audit_config.h"

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <new>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

#ifndef LEO2_ENABLE_TEST_HOOKS
#error "The direct-encode test requires LEO2_ENABLE_TEST_HOOKS"
#endif

#ifndef LEO2_EXPERIMENT_HIGH_HALF_TAIL_COLUMN
#define LEO2_EXPERIMENT_HIGH_HALF_TAIL_COLUMN 1
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
    void* result = malloc(bytes == 0 ? 1 : bytes);
    if (!result)
        throw std::bad_alloc();
    return result;
}

LEO2_TEST_NOINLINE void* operator new[](size_t bytes)
{
    return ::operator new(bytes);
}

LEO2_TEST_NOINLINE void* operator new(size_t bytes, const std::nothrow_t&) noexcept
{
    try { return ::operator new(bytes); }
    catch (...) { return NULL; }
}

LEO2_TEST_NOINLINE void* operator new[](size_t bytes, const std::nothrow_t&) noexcept
{
    return ::operator new(bytes, std::nothrow);
}

LEO2_TEST_NOINLINE void operator delete(void* pointer) noexcept { free(pointer); }
LEO2_TEST_NOINLINE void operator delete[](void* pointer) noexcept { free(pointer); }
LEO2_TEST_NOINLINE void operator delete(void* pointer, const std::nothrow_t&) noexcept
{
    ::operator delete(pointer);
}
LEO2_TEST_NOINLINE void operator delete[](void* pointer, const std::nothrow_t&) noexcept
{
    ::operator delete[](pointer);
}

#undef LEO2_TEST_NOINLINE
#endif

static void begin_allocation_audit()
{
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
    g_tracked_allocations.store(0, std::memory_order_relaxed);
    g_track_allocations.store(true, std::memory_order_release);
#endif
}

static uint64_t end_allocation_audit()
{
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
    g_track_allocations.store(false, std::memory_order_release);
    return g_tracked_allocations.load(std::memory_order_relaxed);
#else
    return 0;
#endif
}

namespace {

using leopard2_test::BinaryField;
using leopard2_test::Element;
using leopard2_test::Matrix;
using leopard2_test::ProfileLayout;

typedef std::vector<uint8_t> Bytes;
typedef std::vector<Bytes> Shards;

struct Counts
{
    uint64_t profiles;
    uint64_t basis_messages;
    uint64_t random_messages;
    uint64_t parity_symbols;
    uint64_t mask_executions;
    uint64_t boundary_profiles;
    uint64_t allocation_checks;
    uint64_t concurrent_executions;
    uint64_t contract_checks;
    uint64_t dispatch_checks;
    uint64_t unaligned_checks;
    uint64_t batch_executions;
    uint64_t no_copy_checks;
    uint64_t low_partial_output_checks;
    uint64_t high_source_staging_checks;
    uint64_t high_byte_tiling_checks;
    uint64_t gf8_coarse_oracle_checks;
    uint64_t high_small_transform_checks;
    uint64_t high_tail_column_checks;
    uint64_t sparse_auto_promotion_checks;

    Counts()
        : profiles(0), basis_messages(0), random_messages(0), parity_symbols(0)
        , mask_executions(0), boundary_profiles(0), allocation_checks(0)
        , concurrent_executions(0), contract_checks(0), dispatch_checks(0)
        , unaligned_checks(0), batch_executions(0), no_copy_checks(0)
        , low_partial_output_checks(0)
        , high_source_staging_checks(0), high_byte_tiling_checks(0)
        , gf8_coarse_oracle_checks(0)
        , high_small_transform_checks(0)
        , high_tail_column_checks(0)
        , sparse_auto_promotion_checks(0)
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

class AlignedBuffer
{
public:
    explicit AlignedBuffer(size_t bytes)
        : data_(NULL), bytes_(bytes)
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
    const void* data() const { return data_; }
    size_t size() const { return bytes_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);
    void* data_;
    size_t bytes_;
};

struct CodecOwner
{
    CodecOwner() : codec(NULL) {}
    ~CodecOwner() { leo2_codec_destroy(codec); }
    leo2_codec* codec;

private:
    CodecOwner(const CodecOwner&);
    CodecOwner& operator=(const CodecOwner&);
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

Shards random_shards(unsigned count, size_t bytes, uint64_t seed)
{
    Shards result(count, Bytes(bytes, 0));
    for (unsigned shard = 0; shard < count; ++shard)
        for (size_t offset = 0; offset < bytes; ++offset)
            result[shard][offset] = static_cast<uint8_t>(next_random(&seed) >> 56);
    return result;
}

std::vector<const void*> const_pointers(const Shards& shards)
{
    std::vector<const void*> result(shards.size(), NULL);
    for (size_t i = 0; i < shards.size(); ++i)
        result[i] = &shards[i][0];
    return result;
}

uint64_t hash_shards(const Shards& shards)
{
    uint64_t hash = UINT64_C(1469598103934665603);
    for (size_t shard_i = 0; shard_i < shards.size(); ++shard_i)
    {
        const Bytes& shard = shards[shard_i];
        for (size_t byte_i = 0; byte_i < shard.size(); ++byte_i)
        {
            hash ^= shard[byte_i];
            hash *= UINT64_C(1099511628211);
        }
    }
    return hash;
}

leopard2_test::ProfileKind oracle_profile(leo2_profile profile)
{
    return profile == LEO2_PROFILE_LEGACY_HIGH_V1
        ? leopard2_test::kLegacyHigh : leopard2_test::kLow;
}

CodecOwner* make_codec(
    leo2_context* context,
    unsigned k,
    unsigned r,
    leo2_profile profile,
    leo2_field field,
    leo2_shard_layout layout = LEO2_SHARD_LAYOUT_NATIVE_V1)
{
    CodecOwner* owner = new CodecOwner;
    leo2_codec_options options = {};
    options.struct_size = sizeof(options);
    options.shard_layout = layout;
    const leo2_result result = leo2_codec_create(
        context, k, r, profile, field, &options, &owner->codec);
    if (result != LEO2_SUCCESS)
    {
        delete owner;
        require_result(result, "codec create");
    }
    return owner;
}

Shards oracle_parity(
    const BinaryField& field,
    const Matrix& generator,
    const Shards& original,
    unsigned recovery_count,
    leo2_field public_field)
{
    const unsigned k = static_cast<unsigned>(original.size());
    const size_t bytes = original[0].size();
    Shards result(recovery_count, Bytes(bytes, 0));
    if (public_field == LEO2_FIELD_GF8)
    {
        for (size_t offset = 0; offset < bytes; ++offset)
        {
            std::vector<Element> message(k, 0);
            for (unsigned i = 0; i < k; ++i)
                message[i] = original[i][offset];
            const std::vector<Element> codeword =
                leopard2_test::matrix_vector_multiply(field, generator, message);
            for (unsigned parity = 0; parity < recovery_count; ++parity)
                result[parity][offset] =
                    static_cast<uint8_t>(codeword[k + parity]);
        }
    }
    else
    {
        require((bytes & 1u) == 0, "GF16 oracle received an odd physical size");
        for (size_t tile = 0; tile < bytes; tile += 64)
        {
            const size_t tile_bytes = std::min<size_t>(64, bytes - tile);
            const size_t symbols = tile_bytes / 2;
            for (size_t lane = 0; lane < symbols; ++lane)
            {
                std::vector<Element> message(k, 0);
                for (unsigned i = 0; i < k; ++i)
                    message[i] = static_cast<Element>(original[i][tile + lane] |
                        (static_cast<unsigned>(original[i][tile + symbols + lane]) << 8));
                const std::vector<Element> codeword =
                    leopard2_test::matrix_vector_multiply(field, generator, message);
                for (unsigned parity = 0; parity < recovery_count; ++parity)
                {
                    const Element value = codeword[k + parity];
                    result[parity][tile + lane] = static_cast<uint8_t>(value);
                    result[parity][tile + symbols + lane] =
                        static_cast<uint8_t>(value >> 8);
                }
            }
        }
    }
    return result;
}

struct EncodeResult
{
    Shards recovery;
    leo2_result result;
};

EncodeResult encode(
    leo2_codec* codec,
    leo2_test_encode_mode mode,
    const Shards& original,
    const std::vector<uint8_t>& requested,
    uint8_t sentinel = 0xa5)
{
    require_result(leo2_test_codec_set_encode_mode(codec, mode), "set encode mode");
    const size_t bytes = original[0].size();
    const unsigned r = leo2_codec_recovery_count(codec);
    require(requested.size() == r, "encode request bitmap size mismatch");
    size_t scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        "encode scratch query");
    AlignedBuffer scratch(scratch_bytes);
    std::vector<const void*> input = const_pointers(original);
    EncodeResult output;
    output.recovery.assign(r, Bytes(bytes, sentinel));
    std::vector<void*> recovery(r, NULL);
    for (unsigned i = 0; i < r; ++i)
        if (requested[i])
            recovery[i] = &output.recovery[i][0];
    output.result = leo2_encode(codec, bytes, &input[0], &recovery[0],
        scratch.data(), scratch.size());
    return output;
}

void compare_requested(
    const Shards& actual,
    const Shards& expected,
    const std::vector<uint8_t>& requested,
    uint8_t sentinel,
    const std::string& label,
    Counts* counts)
{
    for (unsigned parity = 0; parity < requested.size(); ++parity)
    {
        if (requested[parity])
        {
            require(actual[parity] == expected[parity], label + " parity mismatch");
            counts->parity_symbols += actual[parity].size();
        }
        else
        {
            require(std::find(actual[parity].begin(), actual[parity].end(),
                    static_cast<uint8_t>(sentinel ^ 0xffu)) == actual[parity].end(),
                label + " invalid sentinel check");
            for (size_t i = 0; i < actual[parity].size(); ++i)
                require(actual[parity][i] == sentinel,
                    label + " unrequested output was modified");
        }
    }
}

std::vector<std::vector<uint8_t> > request_masks(unsigned r)
{
    std::vector<std::vector<uint8_t> > masks;
    masks.push_back(std::vector<uint8_t>(r, 1));
    std::vector<uint8_t> prefix(r, 0);
    for (unsigned i = 0; i < std::max(1u, r / 2); ++i)
        prefix[i] = 1;
    masks.push_back(prefix);
    std::vector<uint8_t> holey(r, 0);
    for (unsigned i = 0; i < r; i += 2)
        holey[i] = 1;
    holey[r - 1] = 1;
    masks.push_back(holey);
    std::vector<uint8_t> first(r, 0);
    first[0] = 1;
    masks.push_back(first);
    std::vector<uint8_t> middle(r, 0);
    middle[r / 2] = 1;
    masks.push_back(middle);
    std::vector<uint8_t> last(r, 0);
    last[r - 1] = 1;
    masks.push_back(last);
    std::vector<uint8_t> suffix(r, 0);
    for (unsigned i = r / 2; i < r; ++i)
        suffix[i] = 1;
    masks.push_back(suffix);
    masks.push_back(std::vector<uint8_t>(r, 0));
    return masks;
}

void test_profile_matrix(
    leo2_context* context,
    const BinaryField& field,
    leo2_profile profile,
    leo2_field public_field,
    Counts* counts)
{
    for (unsigned k = 1; k <= 16; ++k)
        for (unsigned r = 1; r <= 16; ++r)
        {
            CodecOwner* owner = make_codec(context, k, r, profile, public_field);
            require(leo2_test_codec_direct_encode_capable(owner->codec) == 1,
                "in-range tiny codec was not direct-capable");
            const ProfileLayout layout = leopard2_test::make_profile_layout(
                oracle_profile(profile), k, r);
            const Matrix generator =
                leopard2_test::direct_systematic_generator(field, layout);
            const size_t symbol_bytes = public_field == LEO2_FIELD_GF8 ? 1 : 2;
            const std::vector<uint8_t> all(r, 1);

            for (unsigned basis = 0; basis < k; ++basis)
            {
                Shards original(k, Bytes(symbol_bytes, 0));
                original[basis][0] = 1;
                const Shards expected = oracle_parity(
                    field, generator, original, r, public_field);
                const EncodeResult direct = encode(owner->codec,
                    LEO2_TEST_ENCODE_FORCE_DIRECT, original, all);
                const EncodeResult transform = encode(owner->codec,
                    LEO2_TEST_ENCODE_FORCE_TRANSFORM, original, all);
                require_result(direct.result, "basis direct encode");
                require_result(transform.result, "basis transform encode");
                compare_requested(direct.recovery, expected, all, 0xa5,
                    "basis direct/oracle", counts);
                require(transform.recovery == direct.recovery,
                    "basis direct and transform encoders differ: K=" +
                        std::to_string(k) + " R=" + std::to_string(r) +
                        " basis=" + std::to_string(basis) +
                        " profile=" +
                            std::to_string(static_cast<unsigned>(profile)) +
                        " field=" +
                            std::to_string(static_cast<unsigned>(public_field)));
                ++counts->basis_messages;
            }

            const size_t random_bytes = public_field == LEO2_FIELD_GF8 ? 7 : 6;
            Shards original = random_shards(k, random_bytes,
                UINT64_C(0x78656e636f6465) ^ (static_cast<uint64_t>(k) << 32) ^
                (static_cast<uint64_t>(r) << 16) ^
                (static_cast<uint64_t>(profile) << 8) ^ public_field);
            const Shards expected = oracle_parity(
                field, generator, original, r, public_field);
            const std::vector<std::vector<uint8_t> > masks = request_masks(r);
            for (size_t mask_i = 0; mask_i < masks.size(); ++mask_i)
            {
                const EncodeResult direct = encode(owner->codec,
                    LEO2_TEST_ENCODE_FORCE_DIRECT, original, masks[mask_i]);
                const EncodeResult transform = encode(owner->codec,
                    LEO2_TEST_ENCODE_FORCE_TRANSFORM, original, masks[mask_i]);
                require_result(direct.result, "random direct encode");
                require_result(transform.result, "random transform encode");
                compare_requested(direct.recovery, expected, masks[mask_i], 0xa5,
                    "random direct/oracle", counts);
                require(transform.recovery == direct.recovery,
                    "random direct and transform output masks differ");
                ++counts->mask_executions;
            }

            require_result(leo2_test_codec_set_encode_mode(
                owner->codec, LEO2_TEST_ENCODE_AUTO), "set AUTO encode mode");
            int direct_path = -1;
            require_result(leo2_test_codec_encode_path(owner->codec,
                random_bytes, r, &direct_path), "AUTO path query");
            int expected_auto_direct = 0;
#if defined(LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE) && \
    LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE_AUTO
            if (profile == LEO2_PROFILE_LEGACY_HIGH_V1 &&
                public_field == LEO2_FIELD_GF8 &&
                leo2_context_backend(context) == LEO2_BACKEND_AVX2 &&
                k == 5 && r == 5)
                expected_auto_direct = 1;
#endif
            require(direct_path == expected_auto_direct,
                "AUTO selected an uncalibrated tiny direct-encode threshold");
            ++counts->random_messages;
            ++counts->profiles;
            delete owner;
        }
}

void test_capability_boundaries(leo2_context* context, Counts* counts)
{
    CodecOwner* mode_validation = make_codec(context, 3, 3,
        LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8);
    require_result(leo2_test_codec_set_encode_mode(mode_validation->codec,
        LEO2_TEST_ENCODE_FORCE_TRANSFORM), "set validation transform mode");
    require(leo2_test_codec_set_encode_mode(mode_validation->codec, -1) ==
            LEO2_INVALID_ARGUMENT &&
            leo2_test_codec_set_encode_mode(mode_validation->codec, 99) ==
            LEO2_INVALID_ARGUMENT,
        "invalid integer encode mode was accepted");
    int validation_direct = -1;
    require_result(leo2_test_codec_encode_path(mode_validation->codec,
        1024, 1, &validation_direct), "post-rejection encode path query");
    require(validation_direct == 0,
        "invalid encode mode changed the prior codec mode");
    delete mode_validation;

    const unsigned pairs[][2] = {
        { 17, 1 }, { 17, 16 }, { 1, 17 }, { 16, 17 }, { 17, 17 }
    };
    const leo2_profile profiles[] = {
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_PROFILE_LOW_V1
    };
    const leo2_field fields[] = { LEO2_FIELD_GF8, LEO2_FIELD_GF16 };
    for (unsigned profile_i = 0; profile_i < 2; ++profile_i)
        for (unsigned field_i = 0; field_i < 2; ++field_i)
            for (unsigned pair_i = 0; pair_i < sizeof(pairs) / sizeof(pairs[0]); ++pair_i)
            {
                CodecOwner* owner = make_codec(context, pairs[pair_i][0],
                    pairs[pair_i][1], profiles[profile_i], fields[field_i]);
                require(leo2_test_codec_direct_encode_capable(owner->codec) == 0,
                    "out-of-range tiny codec was direct-capable");
                require(leo2_test_codec_set_encode_mode(owner->codec,
                    LEO2_TEST_ENCODE_FORCE_DIRECT) == LEO2_UNSUPPORTED,
                    "force-direct accepted K or R equal to 17");
                require_result(leo2_test_codec_set_encode_mode(owner->codec,
                    LEO2_TEST_ENCODE_FORCE_TRANSFORM), "force transform boundary");
                int direct_path = -1;
                const uint64_t bytes = fields[field_i] == LEO2_FIELD_GF8 ? 1 : 2;
                require_result(leo2_test_codec_encode_path(owner->codec, bytes,
                    pairs[pair_i][1], &direct_path), "boundary path query");
                require(direct_path == 0, "forced transform boundary chose direct");
                ++counts->boundary_profiles;
                delete owner;
            }
}

void test_sparse_schedule_budget_fallback(
    leo2_context* context,
    Counts* counts)
{
    const unsigned k = 513;
    const unsigned r = 53000;
    CodecOwner* owner = make_codec(
        context, k, r, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16);
    const Shards original = random_shards(
        k, 2, UINT64_C(0x5350415253454244));
    std::vector<uint8_t> requested(r, 0);
    requested[0] = 1;
    requested[r - 1] = 1;

    const EncodeResult prefix = encode(
        owner->codec, LEO2_TEST_ENCODE_AUTO, original, requested);
    const EncodeResult bounded = encode(
        owner->codec, LEO2_TEST_ENCODE_FORCE_TRANSFORM,
        original, requested);
    require_result(prefix.result, "schedule-budget prefix encode");
    require_result(bounded.result, "schedule-budget forced encode");
    require(prefix.recovery == bounded.recovery,
        "bounded schedule fallback differs from mature prefix encode");
    ++counts->boundary_profiles;
    delete owner;
}

void test_low_transform_no_coefficient_copy(
    leo2_context* context,
    const BinaryField& gf8,
    const BinaryField& gf16,
    Counts* counts)
{
    struct Case
    {
        leo2_field field;
        unsigned k;
        unsigned r;
        size_t bytes;
    };
    const Case cases[] = {
        { LEO2_FIELD_GF8, 1, 5, 65 },
        { LEO2_FIELD_GF8, 2, 5, 65 },
        { LEO2_FIELD_GF8, 3, 8, 64 },
        { LEO2_FIELD_GF8, 5, 19, 129 },
        { LEO2_FIELD_GF8, 9, 32, 64 },
        { LEO2_FIELD_GF8, 17, 64, 64 },
        { LEO2_FIELD_GF16, 1, 5, 66 },
        { LEO2_FIELD_GF16, 2, 5, 66 },
        { LEO2_FIELD_GF16, 3, 8, 1024 },
        { LEO2_FIELD_GF16, 5, 19, 130 },
        { LEO2_FIELD_GF16, 9, 32, 1024 },
        { LEO2_FIELD_GF16, 17, 64, 1024 }
    };
    for (unsigned case_i = 0;
         case_i < sizeof(cases) / sizeof(cases[0]); ++case_i)
    {
        const Case& c = cases[case_i];
        CodecOwner* owner = make_codec(context, c.k, c.r,
            LEO2_PROFILE_LOW_V1, c.field);
        const Shards original = random_shards(c.k, c.bytes,
            UINT64_C(0x4e4f434f50590000) + case_i);
        const Shards original_before = original;
        const BinaryField& field = c.field == LEO2_FIELD_GF8 ? gf8 : gf16;
        const ProfileLayout layout = leopard2_test::make_profile_layout(
            leopard2_test::kLow, c.k, c.r);
        const Matrix generator =
            leopard2_test::direct_systematic_generator(field, layout);
        const Shards expected = oracle_parity(
            field, generator, original, c.r, c.field);

        std::vector<uint8_t> requested(c.r, 0);
        requested[0] = 1;
        requested[c.k > 2 ? 2 : 1] = 1;
        const unsigned p = leo2_codec_padded_side(owner->codec);
        requested[p - 1] = 1;
        if (p < c.r)
            requested[p] = 1;
        requested[c.r - 1] = 1;

        if (c.field == LEO2_FIELD_GF8)
            leopard::ff8::TestOnlyResetLowEncodeCounts();
        else
            leopard::ff16::TestOnlyResetLowEncodeCounts();
        const EncodeResult transformed = encode(owner->codec,
            LEO2_TEST_ENCODE_AUTO, original, requested);
        require_result(transformed.result, "no-copy low transform encode");
        require(original == original_before,
            "low transform modified a systematic source shard");
        compare_requested(transformed.recovery, expected, requested, 0xa5,
            "no-copy low transform/oracle", counts);

        uint64_t nonempty_blocks = 0;
        for (unsigned offset = 0; offset < c.r; offset += p)
        {
            const unsigned end = std::min(c.r, offset + p);
            for (unsigned i = offset; i < end; ++i)
                if (requested[i])
                {
                    ++nonempty_blocks;
                    break;
                }
        }
        // Ragged public shards execute the aligned prefix and one padded tail
        // tile as independent transform passes.  Count both without treating
        // the padding bytes as a second public shard.
        const uint64_t transform_passes =
            (c.bytes >= 64 ? 1U : 0U) + ((c.bytes & 63U) != 0 ? 1U : 0U);
        const uint64_t executed_blocks = nonempty_blocks * transform_passes;
        if (c.field == LEO2_FIELD_GF8)
        {
            const leopard::ff8::TestOnlyLowEncodeCounts actual =
                leopard::ff8::TestOnlyGetLowEncodeCounts();
            /*
                The property this protects is the no-copy contract: every
                coefficient enters the transform exactly once through an
                out-of-place first layer, with no staging copy.  Assert that
                directly in coefficients rather than pinning one operation's
                call count, so a backend that publishes a wider out-of-place
                first stage (radix-eight covers eight coefficients per call
                instead of four) still has to satisfy it exactly.
            */
            const uint64_t staged_coefficients =
                actual.fft_butterfly2_out_of_place * 2 +
                actual.fft_butterfly4_out_of_place * 4 +
                actual.fft_butterfly8_out_of_place * 8;
            if (p == 1)
                require(staged_coefficients == 0,
                    "GF8 P=1 entered a transform first layer");
            else
                require(staged_coefficients == executed_blocks * p,
                    "GF8 out-of-place first layer did not stage every "
                    "coefficient exactly once");
            if (p == 2)
                require(actual.fft_butterfly2_out_of_place == executed_blocks &&
                        actual.fft_butterfly4_out_of_place == 0 &&
                        actual.fft_butterfly8_out_of_place == 0,
                    "GF8 P=2 out-of-place call count mismatch");
            else if (p >= 4)
                require(actual.fft_butterfly2_out_of_place == 0,
                    "GF8 fused first layer fell back to radix two");
        }
        else
        {
            const leopard::ff16::TestOnlyLowEncodeCounts actual =
                leopard::ff16::TestOnlyGetLowEncodeCounts();
            if (p == 1)
                require(actual.fft_butterfly2_out_of_place == 0 &&
                        actual.fft_butterfly4_out_of_place == 0,
                    "GF16 P=1 entered a transform first layer");
            else if (p == 2)
                require(actual.fft_butterfly2_out_of_place == executed_blocks &&
                        actual.fft_butterfly4_out_of_place == 0,
                    "GF16 P=2 out-of-place call count mismatch");
            else
            {
                const size_t prefix_bytes = c.bytes & ~size_t(63);
                const leo2_backend backend = leo2_context_backend(context);
                const bool prefix_fused = prefix_bytes == 64 ||
                    (prefix_bytes == 128 &&
                     (backend == LEO2_BACKEND_AVX2 ||
                      backend == LEO2_BACKEND_AVX512));
                const uint64_t fused_passes =
                    ((c.bytes & 63U) != 0 ? 1U : 0U) +
                    (prefix_bytes != 0 && prefix_fused ? 1U : 0U);
                const uint64_t split_passes =
                    prefix_bytes != 0 && !prefix_fused ? 1U : 0U;
                require(actual.fft_butterfly2_out_of_place ==
                            nonempty_blocks * split_passes * (p / 2) &&
                        actual.fft_butterfly4_out_of_place ==
                            nonempty_blocks * fused_passes * (p / 4),
                    "GF16 out-of-place first-layer call count mismatch");
            }
        }

        if (p >= 2)
        {
            std::vector<uint8_t> dense_requested(c.r, 1);
            if (c.field == LEO2_FIELD_GF8)
                leopard::ff8::TestOnlyResetLowEncodeCounts();
            else
                leopard::ff16::TestOnlyResetLowEncodeCounts();
            const EncodeResult dense = encode(owner->codec,
                LEO2_TEST_ENCODE_AUTO, original, dense_requested);
            require_result(dense.result, "dense low transform encode");
            compare_requested(dense.recovery, expected, dense_requested, 0xa5,
                "dense low direct-output/oracle", counts);

            const uint64_t expected_direct_blocks =
                (c.r / p) * transform_passes;
            const uint64_t actual_direct_blocks = c.field == LEO2_FIELD_GF8
                ? leopard::ff8::TestOnlyGetLowEncodeCounts()
                    .direct_output_blocks
                : leopard::ff16::TestOnlyGetLowEncodeCounts()
                    .direct_output_blocks;
            require(actual_direct_blocks == expected_direct_blocks,
                "dense low direct-output block count mismatch");
        }
        ++counts->no_copy_checks;
        delete owner;
    }
}

void test_low_p16_partial_direct_output(
    const BinaryField& gf8,
    Counts* counts)
{
    leo2_context_options avx2_options = {};
    avx2_options.struct_size = sizeof(avx2_options);
    avx2_options.backend = LEO2_BACKEND_AVX2;
    avx2_options.thread_count = 1;
    leo2_context* avx2_context = NULL;
    const leo2_result avx2_result = leo2_context_create(
        &avx2_options, &avx2_context);
    if (avx2_result == LEO2_UNSUPPORTED)
        return;
    require_result(avx2_result, "P16 partial-output AVX2 context create");
    require(leo2_context_backend(avx2_context) == LEO2_BACKEND_AVX2,
        "P16 partial-output context did not resolve to AVX2");

    leo2_context_options scalar_options = {};
    scalar_options.struct_size = sizeof(scalar_options);
    scalar_options.backend = LEO2_BACKEND_SCALAR;
    scalar_options.thread_count = 1;
    leo2_context* scalar_context = NULL;
    require_result(leo2_context_create(&scalar_options, &scalar_context),
        "P16 partial-output scalar context create");
    require(leo2_context_backend(scalar_context) == LEO2_BACKEND_SCALAR,
        "P16 partial-output context did not resolve to scalar");

    struct ModeGuard
    {
        explicit ModeGuard(bool initial) : initial_(initial) {}
        ~ModeGuard()
        {
            leopard::ff8::TestOnlySetLowP16PartialDirectOutputEnabled(initial_);
        }
        bool initial_;
    } mode_guard(
        leopard::ff8::TestOnlyLowP16PartialDirectOutputEnabled());
    (void)mode_guard;

    struct Shape
    {
        unsigned k;
        unsigned r;
    };
    const Shape shapes[] = {
        { 16, 4 }, { 16, 8 }, { 16, 12 },
        { 9, 20 }, { 9, 24 }, { 9, 28 }
    };
    const size_t byte_counts[] = {
        1, 3, 31, 32, 33, 63, 64, 65,
        127, 128, 129, 255, 256, 257, 1024
    };
    for (size_t shape_i = 0;
         shape_i < sizeof(shapes) / sizeof(shapes[0]); ++shape_i)
    {
        const Shape& shape = shapes[shape_i];
        const unsigned partial = shape.r & 15U;
        require(partial == 4 || partial == 8 || partial == 12,
            "P16 partial-output test shape is not aligned");
        CodecOwner* avx2 = make_codec(avx2_context, shape.k, shape.r,
            LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8);
        CodecOwner* scalar = make_codec(scalar_context, shape.k, shape.r,
            LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8);
        require(leo2_codec_padded_side(avx2->codec) == 16,
            "P16 partial-output test selected another side");
        const ProfileLayout layout = leopard2_test::make_profile_layout(
            leopard2_test::kLow, shape.k, shape.r);
        const Matrix generator =
            leopard2_test::direct_systematic_generator(gf8, layout);
        const std::vector<uint8_t> dense(shape.r, 1);

        for (size_t bytes_i = 0;
             bytes_i < sizeof(byte_counts) / sizeof(byte_counts[0]); ++bytes_i)
        {
            const size_t bytes = byte_counts[bytes_i];
            const uint64_t seed = UINT64_C(0x5031365041525400) +
                static_cast<uint64_t>(shape_i) * 0x10000U + bytes;
            const Shards original = random_shards(shape.k, bytes, seed);
            const Shards original_before = original;
            const Shards expected = oracle_parity(
                gf8, generator, original, shape.r, LEO2_FIELD_GF8);

            leopard::ff8::TestOnlySetLowP16PartialDirectOutputEnabled(false);
            leopard::ff8::TestOnlyResetLowEncodeCounts();
            const EncodeResult mature = encode(avx2->codec,
                LEO2_TEST_ENCODE_AUTO, original, dense);
            require_result(mature.result, "P16 mature transform encode");
            compare_requested(mature.recovery, expected, dense, 0xa5,
                "P16 mature transform/oracle", counts);
            const leopard::ff8::TestOnlyLowEncodeCounts mature_counts =
                leopard::ff8::TestOnlyGetLowEncodeCounts();

            leopard::ff8::TestOnlySetLowP16PartialDirectOutputEnabled(true);
            leopard::ff8::TestOnlyResetLowEncodeCounts();
            const EncodeResult direct = encode(avx2->codec,
                LEO2_TEST_ENCODE_AUTO, original, dense);
            require_result(direct.result, "P16 direct-output encode");
            compare_requested(direct.recovery, expected, dense, 0xa5,
                "P16 direct-output/oracle", counts);
            require(direct.recovery == mature.recovery,
                "P16 direct output differs from mature transform");
            require(original == original_before,
                "P16 direct output modified a systematic source");
            const leopard::ff8::TestOnlyLowEncodeCounts direct_counts =
                leopard::ff8::TestOnlyGetLowEncodeCounts();

            const uint64_t passes =
                (bytes >= 64 ? 1U : 0U) + ((bytes & 63U) != 0 ? 1U : 0U);
            const uint64_t rounded = (bytes + 63U) & ~UINT64_C(63);
            const uint64_t full_blocks = shape.r / 16U;
            if (mature_counts.direct_partial_output_blocks != 0 ||
                direct_counts.direct_partial_output_blocks != passes)
            {
                std::ostringstream message;
                message << "P16 partial-output block route count mismatch: K="
                        << shape.k << " R=" << shape.r << " bytes=" << bytes
                        << " mature="
                        << mature_counts.direct_partial_output_blocks
                        << " direct="
                        << direct_counts.direct_partial_output_blocks
                        << " expected=" << passes;
                require(false, message.str());
            }
            require(mature_counts.direct_output_blocks ==
                        full_blocks * passes &&
                    direct_counts.direct_output_blocks ==
                        (full_blocks + 1U) * passes,
                "P16 direct-output total block count mismatch");
            require(direct_counts.direct_output_shards -
                        mature_counts.direct_output_shards ==
                        static_cast<uint64_t>(partial) * passes,
                "P16 direct-output shard count mismatch");
            require(direct_counts.avoided_scatter_bytes -
                        mature_counts.avoided_scatter_bytes ==
                        2U * static_cast<uint64_t>(partial) * rounded,
                "P16 direct-output traffic accounting mismatch");
            require(mature_counts.scatter_bytes ==
                        2U * static_cast<uint64_t>(partial) * rounded &&
                    direct_counts.scatter_bytes == 0,
                "P16 actual scatter traffic did not match the eliminated pass");
            require(direct_counts.direct_output_butterfly4_out_of_place -
                        mature_counts.direct_output_butterfly4_out_of_place ==
                        static_cast<uint64_t>(partial / 4U) * passes,
                "P16 direct-output radix-four count mismatch");

            leopard::ff8::TestOnlyResetLowEncodeCounts();
            const EncodeResult scalar_result = encode(scalar->codec,
                LEO2_TEST_ENCODE_AUTO, original, dense);
            require_result(scalar_result.result,
                "P16 scalar fallback encode");
            compare_requested(scalar_result.recovery, expected, dense, 0xa5,
                "P16 scalar fallback/oracle", counts);
            require(leopard::ff8::TestOnlyGetLowEncodeCounts()
                        .direct_partial_output_blocks == 0,
                "P16 scalar fallback entered the AVX2 partial route");

            std::vector<uint8_t> holey = dense;
            holey[shape.r - 2U] = 0;
            leopard::ff8::TestOnlyResetLowEncodeCounts();
            leopard::ff8::TestOnlyResetSparseEncodeCounts();
            const EncodeResult sparse = encode(avx2->codec,
                LEO2_TEST_ENCODE_FORCE_TRANSFORM, original, holey);
            require_result(sparse.result, "P16 sparse fallback encode");
            compare_requested(sparse.recovery, expected, holey, 0xa5,
                "P16 sparse fallback/oracle", counts);
            require(leopard::ff8::TestOnlyGetLowEncodeCounts()
                        .direct_partial_output_blocks == 0,
                "P16 sparse block entered the dense partial route");
            require(leopard::ff8::TestOnlyGetSparseEncodeCounts()
                        .exact_blocks != 0,
                "P16 forced sparse fallback did not execute a sparse plan");
            ++counts->low_partial_output_checks;
        }
        delete scalar;
        delete avx2;
    }

    // Every aligned GF8 parity coset has a distinct transform shift.  Sweep
    // the complete legal P=16 range so the partial final callback is checked
    // at offsets 0,16,...,224 rather than only at the first two blocks.
    for (unsigned recovery_offset = 0; recovery_offset <= 224;
         recovery_offset += 16)
    {
        const unsigned partials[] = { 4, 8, 12 };
        for (unsigned partial_i = 0; partial_i < 3; ++partial_i)
        {
            const unsigned k = 9;
            const unsigned r = recovery_offset + partials[partial_i];
            const size_t bytes = 64;
            CodecOwner* owner = make_codec(avx2_context, k, r,
                LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8);
            require(leo2_codec_padded_side(owner->codec) == 16,
                "P16 coset sweep selected another side");
            const ProfileLayout layout = leopard2_test::make_profile_layout(
                leopard2_test::kLow, k, r);
            const Matrix generator =
                leopard2_test::direct_systematic_generator(gf8, layout);
            const Shards original = random_shards(k, bytes,
                UINT64_C(0x503136434f534554) + recovery_offset +
                    partials[partial_i]);
            const Shards expected = oracle_parity(
                gf8, generator, original, r, LEO2_FIELD_GF8);
            const std::vector<uint8_t> dense(r, 1);

            leopard::ff8::TestOnlySetLowP16PartialDirectOutputEnabled(true);
            leopard::ff8::TestOnlyResetLowEncodeCounts();
            const EncodeResult actual = encode(owner->codec,
                LEO2_TEST_ENCODE_AUTO, original, dense);
            require_result(actual.result, "P16 parity-coset sweep encode");
            compare_requested(actual.recovery, expected, dense, 0xa5,
                "P16 parity-coset sweep/oracle", counts);
            const leopard::ff8::TestOnlyLowEncodeCounts route =
                leopard::ff8::TestOnlyGetLowEncodeCounts();
            require(route.direct_partial_output_blocks == 1 &&
                    route.direct_output_blocks == recovery_offset / 16U + 1U &&
                    route.scatter_bytes == 0,
                "P16 parity-coset sweep selected the wrong output route");
            ++counts->low_partial_output_checks;
            delete owner;
        }
    }

    // The ordinary benchmark switches candidate/control inside one frozen
    // executable.  Exercise its armed route probe and normalized timing modes
    // here so a benchmark cannot silently measure the same path twice.
    {
        const unsigned k = 16;
        const unsigned r = 4;
        const size_t bytes = 64;
        CodecOwner* owner = make_codec(avx2_context, k, r,
            LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8);
        const ProfileLayout layout = leopard2_test::make_profile_layout(
            leopard2_test::kLow, k, r);
        const Matrix generator =
            leopard2_test::direct_systematic_generator(gf8, layout);
        const Shards original = random_shards(
            k, bytes, UINT64_C(0x50313650524f4245));
        const Shards expected = oracle_parity(
            gf8, generator, original, r, LEO2_FIELD_GF8);
        const std::vector<uint8_t> dense(r, 1);
        for (unsigned enabled = 0; enabled <= 1; ++enabled)
        {
            require(leopard2_internal::
                    SetLowP16PartialDirectOutputEnabledForDiagnostics(
                        enabled != 0),
                "P16 diagnostic route probe could not be armed");
            const EncodeResult actual = encode(owner->codec,
                LEO2_TEST_ENCODE_AUTO, original, dense);
            require_result(actual.result, "P16 diagnostic probe encode");
            compare_requested(actual.recovery, expected, dense, 0xa5,
                "P16 diagnostic probe/oracle", counts);
            require(leopard2_internal::
                    LowP16PartialDirectOutputRouteSelectedForDiagnostics() ==
                        (enabled != 0),
                "P16 diagnostic probe reported the wrong route");
            require(leopard2_internal::
                    FinishLowP16PartialDirectOutputRouteProbeForDiagnostics(),
                "P16 diagnostic route probe could not finish");
            require(leopard2_internal::
                    LowP16PartialDirectOutputModeForDiagnostics() ==
                        (enabled != 0 ? 1U : 2U),
                "P16 diagnostic timing mode did not normalize");
            ++counts->low_partial_output_checks;
        }
        delete owner;
    }

    // Exercise the qualified path with unaligned inputs and outputs, a padded
    // tail pass, guards, allowed source aliasing, failure atomicity, and no
    // execution allocation.
    {
        const unsigned k = 9;
        const unsigned r = 28;
        const size_t bytes = 65;
        const size_t input_prefix = 3;
        const size_t output_prefix = 5;
        const size_t suffix = 11;
        const uint8_t canary = 0xd7;
        CodecOwner* owner = make_codec(avx2_context, k, r,
            LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8);
        require_result(leo2_test_codec_set_encode_mode(owner->codec,
            LEO2_TEST_ENCODE_AUTO), "select P16 guarded AUTO transform");
        const ProfileLayout layout = leopard2_test::make_profile_layout(
            leopard2_test::kLow, k, r);
        const Matrix generator =
            leopard2_test::direct_systematic_generator(gf8, layout);
        const Shards original = random_shards(
            k, bytes, UINT64_C(0x5031364755415244));
        const Shards expected = oracle_parity(
            gf8, generator, original, r, LEO2_FIELD_GF8);
        Shards input_storage(k, Bytes(input_prefix + bytes + suffix, canary));
        Shards output_storage(r, Bytes(output_prefix + bytes + suffix, canary));
        std::vector<const void*> input(k, NULL);
        std::vector<void*> output(r, NULL);
        for (unsigned i = 0; i < k; ++i)
        {
            memcpy(&input_storage[i][input_prefix], &original[i][0], bytes);
            input[i] = &input_storage[i][input_prefix];
        }
        for (unsigned i = 0; i < r; ++i)
            output[i] = &output_storage[i][output_prefix];
        size_t scratch_bytes = 0;
        require_result(leo2_encode_scratch_size(
            owner->codec, bytes, &scratch_bytes),
            "P16 guarded scratch query");
        const size_t scratch_guard = leo2_scratch_alignment();
        AlignedBuffer scratch_storage(
            scratch_bytes + scratch_guard * 2U);
        uint8_t* const scratch =
            static_cast<uint8_t*>(scratch_storage.data()) + scratch_guard;
        memset(scratch_storage.data(), canary, scratch_storage.size());
        leopard::ff8::TestOnlyResetLowEncodeCounts();
        begin_allocation_audit();
        const leo2_result encoded = leo2_encode(owner->codec, bytes,
            &input[0], &output[0], scratch, scratch_bytes);
        const uint64_t allocations = end_allocation_audit();
        require_result(encoded, "P16 guarded direct-output encode");
        require(allocations == 0,
            "P16 direct-output execution allocated C++ storage");
        require(leopard::ff8::TestOnlyGetLowEncodeCounts()
                    .direct_partial_output_blocks == 2,
            "P16 guarded tail did not use both partial-output passes");
        for (unsigned i = 0; i < r; ++i)
            require(memcmp(&output_storage[i][output_prefix],
                        &expected[i][0], bytes) == 0,
                "P16 guarded direct output differs from oracle");

        Shards aliased_original = original;
        aliased_original[1] = aliased_original[0];
        const Shards aliased_expected = oracle_parity(
            gf8, generator, aliased_original, r, LEO2_FIELD_GF8);
        input[1] = input[0];
        for (unsigned i = 0; i < r; ++i)
            std::fill(output_storage[i].begin() + output_prefix,
                output_storage[i].begin() + output_prefix + bytes, canary);
        require_result(leo2_encode(owner->codec, bytes, &input[0], &output[0],
            scratch, scratch_bytes), "P16 aliased-source encode");
        for (unsigned i = 0; i < r; ++i)
            require(memcmp(&output_storage[i][output_prefix],
                        &aliased_expected[i][0], bytes) == 0,
                "P16 allowed source alias changed parity");

        input[1] = &input_storage[1][input_prefix];
        const Shards input_before = input_storage;
        const Shards output_before = output_storage;
        output[0] = const_cast<void*>(input[0]);
        require(leo2_encode(owner->codec, bytes, &input[0], &output[0],
                    scratch, scratch_bytes) == LEO2_OVERLAP,
            "P16 direct output accepted an input/output overlap");
        require(input_storage == input_before &&
                output_storage == output_before,
            "P16 overlap failure modified caller storage");
        output[0] = &output_storage[0][output_prefix];

        const uint8_t* const scratch_storage_bytes =
            static_cast<const uint8_t*>(scratch_storage.data());
        require(std::all_of(scratch_storage_bytes,
                    scratch_storage_bytes + scratch_guard,
                    [](uint8_t value) { return value == canary; }) &&
                std::all_of(scratch_storage_bytes + scratch_guard +
                        scratch_bytes,
                    scratch_storage_bytes + scratch_storage.size(),
                    [](uint8_t value) { return value == canary; }),
            "P16 direct output changed a scratch guard");

        for (unsigned i = 0; i < k; ++i)
            require(std::all_of(input_storage[i].begin(),
                        input_storage[i].begin() + input_prefix,
                        [](uint8_t value) { return value == canary; }) &&
                    std::all_of(input_storage[i].begin() +
                            input_prefix + bytes, input_storage[i].end(),
                        [](uint8_t value) { return value == canary; }),
                "P16 direct output changed an input guard");
        for (unsigned i = 0; i < r; ++i)
            require(std::all_of(output_storage[i].begin(),
                        output_storage[i].begin() + output_prefix,
                        [](uint8_t value) { return value == canary; }) &&
                    std::all_of(output_storage[i].begin() +
                            output_prefix + bytes, output_storage[i].end(),
                        [](uint8_t value) { return value == canary; }),
                "P16 direct output changed an output guard");
        counts->allocation_checks += 1;
        counts->contract_checks += 2;
        counts->unaligned_checks += k + r;
        counts->low_partial_output_checks += 6;

        // One immutable codec may execute the direct-output schedule
        // concurrently when callers own disjoint scratch and destinations.
        const std::vector<const void*> aligned_input = const_pointers(original);
        struct Invocation
        {
            Invocation(unsigned count, size_t shard_bytes, size_t scratch_size)
                : recovery(count, Bytes(shard_bytes, 0))
                , pointers(count, NULL), scratch(scratch_size)
            {
                for (unsigned i = 0; i < count; ++i)
                    pointers[i] = &recovery[i][0];
            }
            Shards recovery;
            std::vector<void*> pointers;
            AlignedBuffer scratch;
        };
        const unsigned worker_count = 2;
        const unsigned repeats = 4;
        std::vector<Invocation*> invocations(worker_count, NULL);
        for (unsigned worker = 0; worker < worker_count; ++worker)
            invocations[worker] = new Invocation(r, bytes, scratch_bytes);
        std::atomic<unsigned> failures(0);
        std::atomic<unsigned> ready(0);
        std::atomic<bool> start(false);
        std::vector<std::thread> workers;
        for (unsigned worker = 0; worker < worker_count; ++worker)
        {
            workers.push_back(std::thread([&, worker]() {
                ready.fetch_add(1, std::memory_order_release);
                while (!start.load(std::memory_order_acquire))
                    std::this_thread::yield();
                for (unsigned repeat = 0; repeat < repeats; ++repeat)
                {
                    Invocation* call = invocations[worker];
                    if (leo2_encode(owner->codec, bytes, &aligned_input[0],
                            &call->pointers[0], call->scratch.data(),
                            call->scratch.size()) != LEO2_SUCCESS ||
                        call->recovery != expected)
                    {
                        failures.fetch_add(1, std::memory_order_relaxed);
                    }
                }
            }));
        }
        while (ready.load(std::memory_order_acquire) != worker_count)
            std::this_thread::yield();
        start.store(true, std::memory_order_release);
        for (size_t worker = 0; worker < workers.size(); ++worker)
            workers[worker].join();
        require(failures.load(std::memory_order_relaxed) == 0,
            "concurrent P16 direct output was nondeterministic");
        for (unsigned worker = 0; worker < worker_count; ++worker)
            delete invocations[worker];
        counts->concurrent_executions += worker_count * repeats;
        delete owner;
    }

    leo2_context_destroy(scalar_context);
    leo2_context_destroy(avx2_context);
}

void test_high_transform_source_staging(
    leo2_context* context,
    const BinaryField& gf8,
    const BinaryField& gf16,
    Counts* counts)
{
    struct Case
    {
        leo2_field field;
        unsigned k;
        uint64_t expected_out_of_place_calls;
        uint64_t expected_copied_input_shards;
    };
    const Case cases[] = {
        { LEO2_FIELD_GF8, 17, 4, 0 },
        { LEO2_FIELD_GF8, 32, 8, 0 },
        { LEO2_FIELD_GF8, 33, 8, 1 },
        { LEO2_FIELD_GF16, 32, 8, 0 },
        { LEO2_FIELD_GF16, 33, 8, 1 }
    };
    const unsigned r = 16;
    const size_t bytes = 64;
    for (unsigned case_i = 0;
         case_i < sizeof(cases) / sizeof(cases[0]); ++case_i)
    {
        const Case& c = cases[case_i];
        CodecOwner* owner = make_codec(context, c.k, r,
            LEO2_PROFILE_LEGACY_HIGH_V1, c.field);
        const Shards original = random_shards(c.k, bytes,
            UINT64_C(0x4849474853544147) + case_i);
        const Shards original_before = original;
        const BinaryField& field = c.field == LEO2_FIELD_GF8 ? gf8 : gf16;
        const ProfileLayout layout = leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh, c.k, r);
        const Matrix generator =
            leopard2_test::direct_systematic_generator(field, layout);
        const Shards expected = oracle_parity(
            field, generator, original, r, c.field);
        const std::vector<uint8_t> requested(r, 1);

        if (c.field == LEO2_FIELD_GF8)
            leopard::ff8::TestOnlyResetHighEncodeCounts();
        else
            leopard::ff16::TestOnlyResetHighEncodeCounts();
        const EncodeResult transformed = encode(owner->codec,
            LEO2_TEST_ENCODE_FORCE_TRANSFORM, original, requested);
        require_result(transformed.result, "high source-staging encode");
        require(original == original_before,
            "high transform modified a caller source shard");
        compare_requested(transformed.recovery, expected, requested, 0xa5,
            "high source-staging/oracle", counts);

        uint64_t out_of_place_calls = 0;
        uint64_t copied_input_shards = 0;
        if (c.field == LEO2_FIELD_GF8)
        {
            const leopard::ff8::TestOnlyHighEncodeCounts actual =
                leopard::ff8::TestOnlyGetHighEncodeCounts();
            out_of_place_calls = actual.ifft_butterfly4_out_of_place;
            copied_input_shards = actual.input_copy_shards;
        }
        else
        {
            const leopard::ff16::TestOnlyHighEncodeCounts actual =
                leopard::ff16::TestOnlyGetHighEncodeCounts();
            out_of_place_calls = actual.ifft_butterfly4_out_of_place;
            copied_input_shards = actual.input_copy_shards;
        }
        require(out_of_place_calls == c.expected_out_of_place_calls,
            "high source-staging out-of-place call count mismatch");
        require(copied_input_shards == c.expected_copied_input_shards,
            "high source-staging input-copy count mismatch: field=" +
                std::to_string(static_cast<unsigned>(c.field)) +
                " K=" + std::to_string(c.k) +
                " actual=" + std::to_string(copied_input_shards) +
                " expected=" +
                    std::to_string(c.expected_copied_input_shards));
        ++counts->high_source_staging_checks;
        delete owner;
    }

    // The measured GF16 AVX2 crossover returns every legacy-high T>=256 block
    // to copy-first above 16 KiB.  Compare the complete staged prefix and
    // independently oracle-check the next complete 64-byte tile.  This route
    // test also proves that the crossover cannot change parity.
    {
        const unsigned large_k = 512;
        const unsigned large_r = 256;
        const size_t boundary_bytes = 16U * 1024U;
        const size_t tail_bytes = 64U;
        const size_t fallback_bytes = boundary_bytes + tail_bytes;
        CodecOwner* owner = make_codec(context, large_k, large_r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16);
        const Shards large_original = random_shards(
            large_k, fallback_bytes, UINT64_C(0x4746313647415445));
        Shards boundary_original(large_k, Bytes(boundary_bytes));
        Shards tail_original(large_k, Bytes(tail_bytes));
        for (unsigned i = 0; i < large_k; ++i)
        {
            std::copy(large_original[i].begin(),
                large_original[i].begin() + boundary_bytes,
                boundary_original[i].begin());
            std::copy(large_original[i].begin() + boundary_bytes,
                large_original[i].end(), tail_original[i].begin());
        }
        const std::vector<uint8_t> requested(large_r, 1);

        leopard::ff16::TestOnlyResetHighEncodeCounts();
        const EncodeResult boundary = encode(owner->codec,
            LEO2_TEST_ENCODE_FORCE_TRANSFORM, boundary_original, requested);
        require_result(boundary.result, "GF16 source-staging boundary encode");
        const leopard::ff16::TestOnlyHighEncodeCounts boundary_counts =
            leopard::ff16::TestOnlyGetHighEncodeCounts();
        require(boundary_counts.ifft_butterfly4_out_of_place == 128 &&
                boundary_counts.input_copy_shards == 0,
            "GF16 source-staging boundary route mismatch");

        leopard::ff16::TestOnlyResetHighEncodeCounts();
        const EncodeResult fallback = encode(owner->codec,
            LEO2_TEST_ENCODE_FORCE_TRANSFORM, large_original, requested);
        require_result(fallback.result, "GF16 source-staging fallback encode");
        const leopard::ff16::TestOnlyHighEncodeCounts fallback_counts =
            leopard::ff16::TestOnlyGetHighEncodeCounts();
        const leo2_backend backend = leo2_context_backend(context);
        const bool avx2_family = backend == LEO2_BACKEND_AVX2 ||
            backend == LEO2_BACKEND_AVX512;
        require(fallback_counts.ifft_butterfly4_out_of_place ==
                    (avx2_family ? 0U : 128U) &&
                fallback_counts.input_copy_shards ==
                    (avx2_family ? 512U : 0U),
            "GF16 source-staging fallback route mismatch: backend=" +
                std::to_string(static_cast<unsigned>(
                    leo2_context_backend(context))) +
                " out=" + std::to_string(
                    fallback_counts.ifft_butterfly4_out_of_place) +
                " copied=" + std::to_string(
                    fallback_counts.input_copy_shards));

        const unsigned legacy_work_count =
            leo_encode_work_count(large_k, large_r);
        require(legacy_work_count >= large_r,
            "GF16 fallback legacy work-count query failed");
        Shards legacy_work(legacy_work_count, Bytes(tail_bytes, 0));
        std::vector<const void*> legacy_original(large_k, NULL);
        std::vector<void*> legacy_pointers(legacy_work_count, NULL);
        for (unsigned i = 0; i < large_k; ++i)
            legacy_original[i] = tail_original[i].data();
        for (unsigned i = 0; i < legacy_work_count; ++i)
            legacy_pointers[i] = legacy_work[i].data();
        require(leo_encode(tail_bytes, large_k, large_r,
                    legacy_work_count, legacy_original.data(),
                    legacy_pointers.data()) == Leopard_Success,
            "GF16 fallback legacy tail encode failed");
        for (unsigned i = 0; i < large_r; ++i)
        {
            require(std::equal(boundary.recovery[i].begin(),
                    boundary.recovery[i].end(), fallback.recovery[i].begin()),
                "GF16 source-staging fallback changed prefix parity");
            require(std::equal(legacy_work[i].begin(),
                    legacy_work[i].end(),
                    fallback.recovery[i].begin() + boundary_bytes),
                "GF16 source-staging fallback tail differs from legacy oracle");
        }
        ++counts->high_source_staging_checks;
        delete owner;
    }

    // The large-shard fallback is specific to the legacy-high block
    // accumulator.  A low-profile interpolation with the same GF16 side and
    // byte count must retain direct source staging.
    {
        const unsigned low_k = 129;
        const unsigned low_r = 1;
        const size_t low_bytes = 16U * 1024U + 64U;
        CodecOwner* owner = make_codec(context, low_k, low_r,
            LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16);
        const Shards original = random_shards(
            low_k, low_bytes, UINT64_C(0x4c4f575354414745));
        const Shards original_before = original;
        const size_t oracle_bytes = 64;
        Shards oracle_original(low_k, Bytes(oracle_bytes));
        for (unsigned i = 0; i < low_k; ++i)
            std::copy(original[i].begin(),
                original[i].begin() + oracle_bytes,
                oracle_original[i].begin());
        const BinaryField& field = gf16;
        const ProfileLayout layout = leopard2_test::make_profile_layout(
            leopard2_test::kLow, low_k, low_r);
        const Matrix generator =
            leopard2_test::direct_systematic_generator(field, layout);
        const Shards expected = oracle_parity(
            field, generator, oracle_original, low_r, LEO2_FIELD_GF16);
        const std::vector<uint8_t> requested(low_r, 1);

        leopard::ff16::TestOnlyResetHighEncodeCounts();
        const EncodeResult transformed = encode(owner->codec,
            LEO2_TEST_ENCODE_FORCE_TRANSFORM, original, requested);
        require_result(transformed.result,
            "GF16 low-profile large source-staging encode");
        const leopard::ff16::TestOnlyHighEncodeCounts actual =
            leopard::ff16::TestOnlyGetHighEncodeCounts();
        require(actual.ifft_butterfly4_out_of_place == 32 &&
                actual.input_copy_shards == 1,
            "GF16 low-profile large source-staging route mismatch");
        require(original == original_before,
            "GF16 low-profile large source-staging modified caller data");
        require(std::equal(expected[0].begin(), expected[0].end(),
                transformed.recovery[0].begin()),
            "GF16 low-profile large source-staging/oracle prefix mismatch");
        ++counts->parity_symbols;
        ++counts->high_source_staging_checks;
        delete owner;
    }
}

void test_high_small_coarse_kernel(
    leo2_context* context,
    const BinaryField& gf8,
    Counts* counts)
{
    struct Case
    {
        unsigned k;
        unsigned r;
        size_t bytes;
    };
    const Case cases[] = {
        { 2, 2, 2048 },
        { 254, 2, 2049 },
        // Independently check every register-fused T=4 specialization against
        // the direct systematic generator.  R=3 is the contiguous punctured
        // prefix of the same four-coordinate parent transform.
        { 3, 3, 4096 },
        { 3, 4, 4096 },
        { 5, 3, 2048 },
        { 5, 4, 2048 },
        { 6, 3, 2048 },
        { 6, 4, 2048 },
        { 7, 3, 4096 },
        { 7, 4, 4096 },
        { 10, 3, 2048 },
        { 10, 4, 2048 },
        { 11, 3, 2048 },
        { 11, 4, 2048 },
        { 16, 4, 2048 },
        { 8, 4, 4096 },
        { 12, 3, 4097 },
        { 9, 4, 8192 },
        { 4, 4, 65536 },
        { 64, 2, 4096 },
        { 64, 3, 4097 },
        { 240, 4, 65536 },
        { 252, 4, 65536 }
    };
    for (unsigned case_i = 0;
         case_i < sizeof(cases) / sizeof(cases[0]); ++case_i)
    {
        const Case& c = cases[case_i];
        CodecOwner* owner = make_codec(context, c.k, c.r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        const Shards original = random_shards(c.k, c.bytes,
            UINT64_C(0x534d414c4c540000) + case_i);
        const Shards original_before = original;
        const ProfileLayout layout = leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh, c.k, c.r);
        const Matrix generator =
            leopard2_test::direct_systematic_generator(gf8, layout);
        const Shards expected = oracle_parity(
            gf8, generator, original, c.r, LEO2_FIELD_GF8);
        const std::vector<uint8_t> all(c.r, 1);

        int direct = -1;
        require_result(leo2_test_codec_encode_path(
            owner->codec, c.bytes, c.r, &direct),
            "small-T AUTO path query");
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        const EncodeResult dense = encode(owner->codec,
            LEO2_TEST_ENCODE_AUTO, original, all);
        require_result(dense.result, "small-T dense AUTO encode");
        compare_requested(dense.recovery, expected, all, 0xa5,
            "small-T dense AUTO/oracle", counts);
        require(original == original_before,
            "small-T dense AUTO modified caller input");
        const bool avx2 = leo2_context_backend(context) == LEO2_BACKEND_AVX2;
        require((leopard::ff8::TestOnlyGetHighEncodeCounts().
                    small_transform_calls != 0) == (avx2 && direct == 0),
            "small-T coarse-kernel route disagrees with the backend: K=" +
                std::to_string(c.k) + " R=" + std::to_string(c.r) +
                " bytes=" + std::to_string(c.bytes) + " backend=" +
                std::to_string(static_cast<unsigned>(
                    leo2_context_backend(context))) + " direct=" +
                std::to_string(direct) + " calls=" +
                std::to_string(leopard::ff8::TestOnlyGetHighEncodeCounts().
                    small_transform_calls));

        std::vector<uint8_t> sparse(c.r, 0);
        sparse[0] = 1;
        if (c.r > 2)
            sparse[c.r - 1] = 1;
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        const EncodeResult partial = encode(owner->codec,
            LEO2_TEST_ENCODE_FORCE_TRANSFORM, original, sparse);
        require_result(partial.result, "small-T sparse transform encode");
        compare_requested(partial.recovery, expected, sparse, 0xa5,
            "small-T sparse transform/oracle", counts);
        require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                    small_transform_calls == 0,
            "small-T sparse transform used the dense coarse kernel");

        ++counts->high_small_transform_checks;
        delete owner;
    }
}

void test_high_t2_t4_tiny_multi_block(
    const BinaryField& gf8,
    Counts* counts)
{
    leo2_context_options options = {};
    options.struct_size = sizeof(options);
    options.backend = LEO2_BACKEND_AVX2;
    options.thread_count = 1;
    leo2_context* context = NULL;
    const leo2_result created = leo2_context_create(&options, &context);
    if (created == LEO2_UNSUPPORTED)
        return;
    require_result(created, "tiny multi-block AVX2 context create");

    struct Case
    {
        unsigned k;
        unsigned r;
        size_t bytes;
        bool expect_small_transform;
    };
    const Case cases[] = {
        // A padded tail is deliberately kept off the direct-input callback.
        { 5, 2, 32, false },
        { 5, 2, 63, false },
        // T=2 begins immediately after the K=2..4 packed terminals.
        { 5, 2, 64, true },
        { 11, 2, 65, true },
        { 15, 2, 1024, true },
        { 15, 2, 2047, true },
        { 253, 2, 64, true },
        { 254, 2, 65, true },
        // T=4 begins after the K=3..11 terminal/register families.
        { 12, 3, 64, true },
        { 13, 4, 65, true },
        { 15, 4, 1024, true },
        { 15, 4, 2047, true },
        { 249, 4, 1024, true },
        { 252, 4, 2047, true },
        // The established large-shard callback remains selected at 2 KiB.
        { 5, 2, 2048, true }
    };

    for (unsigned case_i = 0;
         case_i < sizeof(cases) / sizeof(cases[0]); ++case_i)
    {
        const Case& c = cases[case_i];
        CodecOwner* owner = make_codec(context, c.k, c.r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        const Shards original = random_shards(c.k, c.bytes,
            UINT64_C(0x543254344d420000) + case_i);
        const Shards original_before = original;
        const ProfileLayout layout = leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh, c.k, c.r);
        const Matrix generator =
            leopard2_test::direct_systematic_generator(gf8, layout);
        const Shards expected = oracle_parity(
            gf8, generator, original, c.r, LEO2_FIELD_GF8);
        const std::vector<uint8_t> requested(c.r, 1);

        leopard::ff8::TestOnlyResetHighEncodeCounts();
        const EncodeResult actual = encode(owner->codec,
            LEO2_TEST_ENCODE_AUTO, original, requested);
        require_result(actual.result, "tiny multi-block encode");
        compare_requested(actual.recovery, expected, requested, 0xa5,
            "tiny multi-block/direct oracle", counts);
        require(original == original_before,
            "tiny multi-block encode modified caller input");
        const leopard::ff8::TestOnlyHighEncodeCounts route =
            leopard::ff8::TestOnlyGetHighEncodeCounts();
        require((route.small_transform_calls != 0) ==
                c.expect_small_transform,
            "tiny multi-block route mismatch: K=" + std::to_string(c.k) +
                " R=" + std::to_string(c.r) +
                " bytes=" + std::to_string(c.bytes) +
                " calls=" + std::to_string(route.small_transform_calls));
        if (c.expect_small_transform && (c.bytes & 63U) == 0)
        {
            const uint64_t expected_tail_copies = c.r == 2
                ? c.k % 2U
                : (c.k % 4U == 3U ? 0U : c.k % 4U);
            require(route.input_copy_shards == expected_tail_copies,
                "tiny multi-block source-copy count mismatch: K=" +
                    std::to_string(c.k) + " R=" + std::to_string(c.r) +
                    " bytes=" + std::to_string(c.bytes) +
                    " copies=" +
                        std::to_string(route.input_copy_shards));
        }
        ++counts->high_small_transform_checks;
        delete owner;
    }
    leo2_context_destroy(context);
}

// Mirrors ComputeBalancedExecutionTiles in leopard2.cpp.  A qualifying aligned
// prefix is split into ceil(prefix / requested) passes carrying equal 64-byte
// tile counts, and the work rows are sized to the largest pass.  Keeping the
// formula here rather than hand-computed constants lets the size assertions
// below stay exact as the calibrated size range widens.
static size_t balanced_pass_bytes(size_t aligned_prefix, size_t requested)
{
    if (requested == 0 || requested >= aligned_prefix)
        return aligned_prefix;
    const size_t passes = 1U + (aligned_prefix - 1U) / requested;
    const size_t total_tiles = aligned_prefix / 64U;
    const size_t pass_tiles =
        total_tiles / passes + (total_tiles % passes != 0);
    return pass_tiles * 64U;
}

void test_high_gf16_byte_tiling(
    leo2_context* context,
    const BinaryField& gf16,
    Counts* counts)
{
    // T=256 enters the qualified two-pass 32-KiB execution layout at 64 KiB.
    // Exercise the exact multi-message-block target used by the same-source
    // crossover screen and compare every transmitted byte with the old API.
    const unsigned k = 1000;
    const unsigned r = 200;
    const size_t tile_bytes = 32U * 1024U;
    const size_t bytes = tile_bytes * 2U;
    CodecOwner* owner = make_codec(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16);

    size_t one_pass_scratch = 0;
    size_t tiled_scratch = 0;
    size_t neighboring_scratch = 0;
    size_t ragged_scratch = 0;
    require_result(leo2_encode_scratch_size(
        owner->codec, tile_bytes, &one_pass_scratch),
        "GF16 byte-tile one-pass scratch query");
    require_result(leo2_encode_scratch_size(
        owner->codec, bytes, &tiled_scratch),
        "GF16 byte-tile scratch query");
    require_result(leo2_encode_scratch_size(
        owner->codec, bytes + 64U, &neighboring_scratch),
        "GF16 byte-tile neighboring-size scratch query");
    require_result(leo2_encode_scratch_size(
        owner->codec, bytes + 2U, &ragged_scratch),
        "GF16 byte-tile ragged scratch query");
    uint64_t effective_l3_bytes = 0;
    uint64_t live_set_target_bytes = 0;
    uint64_t tile_threshold_bytes = 0;
    require(leopard::backend::TestOnlyGetContextGF16CachePolicy(
            context, &effective_l3_bytes, &live_set_target_bytes,
            &tile_threshold_bytes),
        "GF16 byte-tile cache-policy query");
    (void)live_set_target_bytes;
    (void)tile_threshold_bytes;
    const bool byte_tiling_active = tiled_scratch == one_pass_scratch;
    require(byte_tiling_active ||
            tiled_scratch == one_pass_scratch + 512U * tile_bytes,
        "GF16 byte tiling/fallback scratch geometry mismatch");
    if (byte_tiling_active)
    {
        /*
            Byte tiling is calibrated across shard sizes for GF16 legacy-high
            AVX2 at padded_side 256 and 512, no longer at one exact size.  The
            measured gain rises with the shard rather than peaking at the
            original 64-KiB cell -- K=1000 R=200 is 1.14x at 64 KiB, 1.92x at
            128 KiB, 2.52x at 256 KiB and 2.92x at 1 MiB, with encode scratch
            falling 536.9 MB -> 16.8 MB at 1 MiB.  Neighboring and larger sizes
            therefore stay tiled, each sized to its own balanced pass, and the
            work rows never exceed the requested tile.
        */
        const size_t work_base = one_pass_scratch - 512U * tile_bytes;
        require(tiled_scratch ==
                work_base + 512U * balanced_pass_bytes(bytes, tile_bytes),
            "GF16 byte tiling lost its balanced pass geometry");
        /*
            The exact 64-KiB AUTO/AVX-512 cell is independently calibrated and
            intentionally ignores the general cache budget.  On a context
            created inside the 96-MiB affinity domain, that exact cell remains
            two-pass while its 64-byte neighbor correctly stays one-pass until
            the larger-cache threshold.  Every other active case retains the
            general neighboring-size coverage.
        */
        if (effective_l3_bytes == UINT64_C(32) * 1024 * 1024)
        {
            require(neighboring_scratch ==
                    work_base + 512U *
                        balanced_pass_bytes(bytes + 64U, tile_bytes),
                "GF16 byte tiling stopped covering the neighboring size");
        }
        else
        {
            require(effective_l3_bytes >
                    UINT64_C(32) * 1024 * 1024,
                "GF16 neighboring tile was lost at the fallback cache policy");
            require(neighboring_scratch ==
                    work_base + 512U * (bytes + 64U),
                "exact GF16 AUTO tile has invalid one-pass neighbor");
        }
        require(ragged_scratch > tiled_scratch,
            "GF16 byte tiling lost its ragged staging slots");

        leo2_backend dense_backend = LEO2_BACKEND_AUTO;
        leo2_backend partial_backend = LEO2_BACKEND_AUTO;
        require_result(leo2_test_codec_transform_encode_backend(
            owner->codec, bytes, r, r, &dense_backend),
            "GF16 byte-tile dense backend query");
        require_result(leo2_test_codec_transform_encode_backend(
            owner->codec, bytes, r - 1U, r, &partial_backend),
            "GF16 byte-tile partial backend query");
        const leo2_backend context_backend = leo2_context_backend(context);
        leo2_backend expected_dense_backend = context_backend;
#if defined(LEO2_TEST_AUTO_MAY_USE_AVX512_ENCODE)
        if (context_backend == LEO2_BACKEND_AVX2 &&
            leopard::backend::IsCalibratedAutoAVX512EncodeHost())
            expected_dense_backend = LEO2_BACKEND_AVX512;
#endif
        require(dense_backend == expected_dense_backend &&
                partial_backend == context_backend,
            "GF16 byte-tile/backend dispatch boundary mismatch: dense=" +
                std::to_string(static_cast<unsigned>(dense_backend)) +
                " expected_dense=" + std::to_string(
                    static_cast<unsigned>(expected_dense_backend)) +
                " partial=" +
                std::to_string(static_cast<unsigned>(partial_backend)));
    }
    else
    {
        require(neighboring_scratch == tiled_scratch + 512U * 64U,
            "GF16 one-pass neighboring-size scratch geometry mismatch");
        require(ragged_scratch > tiled_scratch,
            "GF16 one-pass ragged scratch geometry mismatch");
    }

    // Promotion covers the calibrated AVX2 transform path, which the AUTO host
    // context also selects as its baseline ops, and the GFNI member that
    // executes the same AVX2-tier schedule.  Every other explicit backend is
    // uncalibrated and retains one complete byte pass.
    const leo2_backend fixed_backends[] = {
        LEO2_BACKEND_SCALAR,
        LEO2_BACKEND_SSSE3,
        LEO2_BACKEND_AVX2,
        LEO2_BACKEND_AVX512,
        LEO2_BACKEND_NEON,
        LEO2_BACKEND_GFNI
    };
    for (unsigned backend_i = 0;
        backend_i < sizeof(fixed_backends) / sizeof(fixed_backends[0]);
        ++backend_i)
    {
        leo2_context_options fixed_options = {};
        fixed_options.struct_size = sizeof(fixed_options);
        fixed_options.backend = fixed_backends[backend_i];
        fixed_options.thread_count = 1;
        leo2_context* fixed_context = NULL;
        const leo2_result created = leo2_context_create(
            &fixed_options, &fixed_context);
        if (created == LEO2_UNSUPPORTED)
            continue;
        require_result(created, "fixed-backend byte-tile context create");
        // This is a deterministic regression test for the established
        // 32-MiB calibration, independent of the CI worker's host affinity.
        leopard::backend::TestOnlySetContextL3Bytes(
            fixed_context, UINT64_C(32) * 1024 * 1024);
        CodecOwner* fixed_owner = make_codec(fixed_context, k, r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16);
        size_t fixed_small_scratch = 0;
        size_t fixed_large_scratch = 0;
        require_result(leo2_encode_scratch_size(
            fixed_owner->codec, tile_bytes, &fixed_small_scratch),
            "fixed-backend small scratch query");
        require_result(leo2_encode_scratch_size(
            fixed_owner->codec, bytes, &fixed_large_scratch),
            "fixed-backend large scratch query");
        const bool backend_tiles =
            fixed_backends[backend_i] == LEO2_BACKEND_AVX2 ||
            fixed_backends[backend_i] == LEO2_BACKEND_GFNI;
        const size_t expected_work_bytes = backend_tiles
            ? balanced_pass_bytes(bytes, tile_bytes) : bytes;
        require(fixed_large_scratch == fixed_small_scratch +
                    512U * (expected_work_bytes - tile_bytes),
            backend_tiles
                ? "explicit AVX2-tier backend lost the calibrated GF16 byte "
                  "tiling"
                : "explicit backend unexpectedly selected GF16 byte tiling");
        ++counts->dispatch_checks;
        delete fixed_owner;
        leo2_context_destroy(fixed_context);
    }

    // The large differential below exists only to validate the promoted
    // schedule.  Generic GF16 encoding and the one-pass fallback are covered
    // by the ordinary profile matrix and boundary tests.  Avoid imposing this
    // roughly 300-MiB vector on uncalibrated and non-x86 CI hosts that cannot
    // execute the candidate.
    if (!byte_tiling_active)
    {
        ++counts->high_byte_tiling_checks;
        delete owner;
        return;
    }

    const Shards original = random_shards(
        k, bytes, UINT64_C(0x4746313642595445));
    const uint64_t original_hash = hash_shards(original);
    const std::vector<const void*> original_pointers =
        const_pointers(original);

    const unsigned legacy_work_count = leo_encode_work_count(k, r);
    require(legacy_work_count == 512,
        "GF16 byte-tile legacy work-count geometry mismatch");
    Shards expected;
    {
        Shards legacy_work(legacy_work_count, Bytes(bytes, 0));
        std::vector<void*> legacy_pointers(legacy_work_count, NULL);
        for (unsigned i = 0; i < legacy_work_count; ++i)
            legacy_pointers[i] = &legacy_work[i][0];
        require(leo_encode(bytes, k, r, legacy_work_count,
                    &original_pointers[0], &legacy_pointers[0]) ==
                Leopard_Success,
            "GF16 byte-tile legacy compatibility encode failed");

        // Retain only the transmitted parity buffers.  Moving the vector
        // storage avoids copying another R complete 64-KiB shards and frees
        // the unused transform half before the Leopard2 executions begin.
        legacy_work.resize(r);
        expected.swap(legacy_work);
    }

    const std::vector<uint8_t> dense_requested(r, 1);
    leopard::ff16::TestOnlyResetHighEncodeCounts();
    leopard::ff16::TestOnlyHighEncodeCounts route = {};
    {
        const EncodeResult dense = encode(owner->codec,
            LEO2_TEST_ENCODE_FORCE_TRANSFORM, original, dense_requested);
        require_result(dense.result, "GF16 byte-tile dense encode");
        require(hash_shards(original) == original_hash,
            "GF16 byte tiling modified a caller source shard");
        compare_requested(dense.recovery, expected, dense_requested, 0xa5,
            "GF16 byte-tile/legacy", counts);

        route = leopard::ff16::TestOnlyGetHighEncodeCounts();

        // The large compatibility vector is complemented by a mathematically
        // independent direct interpolation of its first GF16 symbol.
        const ProfileLayout layout = leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh, k, r);
        std::vector<Element> points(layout.parent_dimension, 0);
        std::vector<Element> values(layout.parent_dimension, 0);
        for (unsigned i = 0; i < layout.parent_dimension; ++i)
        {
            points[i] = static_cast<Element>(
                layout.systematic_coordinates[i]);
        }
        for (unsigned i = 0; i < k; ++i)
        {
            values[i] = static_cast<Element>(original[i][0] |
                (static_cast<unsigned>(original[i][32]) << 8));
        }
        const Element direct_symbol = leopard2_test::lagrange_evaluate(
            gf16, points, values,
            static_cast<Element>(layout.parity_coordinates[0]));
        const Element encoded_symbol = static_cast<Element>(
            dense.recovery[0][0] |
            (static_cast<unsigned>(dense.recovery[0][32]) << 8));
        require(encoded_symbol == direct_symbol,
            "GF16 byte-tile direct-symbol mismatch");
        ++counts->parity_symbols;
    }

    const leo2_backend backend = leo2_context_backend(context);
    const bool copy_first = backend == LEO2_BACKEND_AVX2 ||
        backend == LEO2_BACKEND_AVX512;
    const uint64_t execution_passes = byte_tiling_active ? 2U : 1U;
    require(route.ifft_butterfly4_out_of_place ==
                (copy_first ? 0U : 250U * execution_passes) &&
            route.input_copy_shards ==
                (copy_first ? 1000U * execution_passes : 0U),
        "GF16 byte tiling changed the full-call source-staging policy");

    {
        std::vector<uint8_t> sparse_requested(r, 0);
        sparse_requested[0] = 1;
        sparse_requested[r / 2] = 1;
        sparse_requested[r - 1] = 1;
        const EncodeResult sparse = encode(owner->codec,
            LEO2_TEST_ENCODE_FORCE_TRANSFORM, original, sparse_requested);
        require_result(sparse.result, "GF16 byte-tile sparse encode");
        compare_requested(sparse.recovery, expected, sparse_requested, 0xa5,
            "GF16 byte-tile sparse/legacy", counts);
    }

    {
        std::vector<uint8_t> prefix_requested(r, 1);
        prefix_requested[r - 1U] = 0;
        const EncodeResult prefix = encode(owner->codec,
            LEO2_TEST_ENCODE_FORCE_TRANSFORM, original, prefix_requested);
        require_result(prefix.result, "GF16 byte-tile prefix encode");
        compare_requested(prefix.recovery, expected, prefix_requested, 0xa5,
            "GF16 byte-tile prefix/legacy", counts);
    }

    {
        AlignedBuffer alias_scratch(tiled_scratch);
        std::vector<void*> aliased_output(r, NULL);
        aliased_output[0] = const_cast<void*>(original_pointers[0]);
        require(leo2_encode(owner->codec, bytes, &original_pointers[0],
                    &aliased_output[0], alias_scratch.data(),
                    alias_scratch.size()) == LEO2_OVERLAP,
            "GF16 byte tiling accepted output/input overlap");
        ++counts->contract_checks;
    }

    // One immutable codec may execute the tiled schedule concurrently when
    // each caller supplies disjoint outputs and scratch.
    struct Invocation
    {
        Invocation(unsigned output_count, size_t shard_bytes,
            size_t scratch_bytes)
            : output(output_count, Bytes(shard_bytes, 0))
            , pointers(output_count, NULL), scratch(scratch_bytes)
        {
            for (unsigned i = 0; i < output_count; ++i)
                pointers[i] = &output[i][0];
        }
        Shards output;
        std::vector<void*> pointers;
        AlignedBuffer scratch;
    };
    const unsigned worker_count = 2;
    const unsigned repeats = 2;
    std::vector<Invocation*> invocation(worker_count, NULL);
    for (unsigned i = 0; i < worker_count; ++i)
        invocation[i] = new Invocation(r, bytes, tiled_scratch);
    std::atomic<unsigned> failures(0);
    std::vector<std::thread> workers;
    for (unsigned worker = 0; worker < worker_count; ++worker)
        workers.push_back(std::thread([&, worker]() {
            for (unsigned repeat = 0; repeat < repeats; ++repeat)
            {
                Invocation* call = invocation[worker];
                if (leo2_encode(owner->codec, bytes, &original_pointers[0],
                        &call->pointers[0], call->scratch.data(),
                        call->scratch.size()) != LEO2_SUCCESS ||
                    call->output != expected)
                {
                    failures.fetch_add(1, std::memory_order_relaxed);
                }
            }
        }));
    for (size_t i = 0; i < workers.size(); ++i)
        workers[i].join();
    require(failures.load(std::memory_order_relaxed) == 0,
        "concurrent GF16 byte tiling was nondeterministic");
    for (unsigned i = 0; i < worker_count; ++i)
        delete invocation[i];
    counts->concurrent_executions += worker_count * repeats;
    require(hash_shards(original) == original_hash,
        "concurrent GF16 byte tiling modified a caller source shard");

    ++counts->high_byte_tiling_checks;
    delete owner;
}

void test_gf8_high_tail_generator_column(
    const BinaryField& gf8,
    Counts* counts)
{
    leo2_context_options context_options = {};
    context_options.struct_size = sizeof(context_options);
    context_options.backend = LEO2_BACKEND_AVX2;
    context_options.thread_count = 1;
    leo2_context* context = NULL;
    const leo2_result context_result = leo2_context_create(
        &context_options, &context);
    if (context_result == LEO2_UNSUPPORTED)
        return;
    require_result(context_result, "explicit AVX2 tail-column context");

    static const unsigned transform_sizes[] = { 8, 16, 32, 64 };
    const size_t bytes = 2051;
    for (unsigned size_i = 0;
         size_i < sizeof(transform_sizes) / sizeof(transform_sizes[0]);
         ++size_i)
    {
        const unsigned t = transform_sizes[size_i];
        const unsigned k = t + 1U;
        CodecOwner* owner = make_codec(context, k, t,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        Shards original(k, Bytes(bytes, 0));
        original[k - 1] = random_shards(
            1, bytes, UINT64_C(0x5441494c434f4c00) + t)[0];
        const Shards original_before = original;
        const ProfileLayout layout = leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh, k, t);
        std::vector<Element> systematic_points(
            layout.systematic_coordinates.begin(),
            layout.systematic_coordinates.end());
        const Element source_coordinate = systematic_points[k - 1];
        Element derivative = 1;
        for (size_t point = 0; point < systematic_points.size(); ++point)
        {
            if (systematic_points[point] != source_coordinate)
            {
                derivative = gf8.multiply(derivative,
                    gf8.add(systematic_points[point], source_coordinate));
            }
        }
        require(derivative != 0,
            "GF8 tail generator derivative was zero");
        Shards expected(t, Bytes(bytes, 0));
        for (unsigned parity = 0; parity < t; ++parity)
        {
            const Element parity_coordinate = static_cast<Element>(
                layout.parity_coordinates[parity]);
            Element vanishing = 1;
            for (size_t point = 0; point < systematic_points.size(); ++point)
            {
                vanishing = gf8.multiply(vanishing,
                    gf8.add(parity_coordinate, systematic_points[point]));
            }
            const Element denominator = gf8.multiply(derivative,
                gf8.add(parity_coordinate, source_coordinate));
            const Element coefficient = gf8.divide(vanishing, denominator);
            require(coefficient != 0,
                "GF8 tail generator column contained zero");
            for (size_t offset = 0; offset < bytes; ++offset)
            {
                expected[parity][offset] = static_cast<uint8_t>(gf8.multiply(
                    coefficient, original[k - 1][offset]));
            }
        }

        const std::vector<uint8_t> all(t, 1);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        const EncodeResult dense = encode(owner->codec,
            LEO2_TEST_ENCODE_FORCE_TRANSFORM, original, all);
        require_result(dense.result, "GF8 tail-column dense encode");
        compare_requested(dense.recovery, expected, all, 0xa5,
            "GF8 tail-column/direct generator", counts);
        require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                    tail_column_calls != 0,
            "GF8 dense tail-column route was not exercised");
        require(original == original_before,
            "GF8 tail-column encode modified caller input");

        std::vector<uint8_t> sparse(t, 0);
        sparse[0] = 1;
        sparse[t / 2] = 1;
        sparse[t - 1] = 1;
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        const EncodeResult partial = encode(owner->codec,
            LEO2_TEST_ENCODE_FORCE_TRANSFORM, original, sparse);
        require_result(partial.result, "GF8 tail-column sparse encode");
        compare_requested(partial.recovery, expected, sparse, 0xa5,
            "GF8 sparse tail-column/direct generator", counts);
        require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                    tail_column_calls == 0,
            "GF8 sparse output used the dense tail-column path");

        ++counts->high_tail_column_checks;
        delete owner;
    }
    leo2_context_destroy(context);
}

#if LEO2_EXPERIMENT_HIGH_HALF_TAIL_COLUMN
void test_gf8_high_half_tail_generator_column(
    const BinaryField& gf8,
    Counts* counts)
{
    leo2_context_options context_options = {};
    context_options.struct_size = sizeof(context_options);
    context_options.backend = LEO2_BACKEND_AVX2;
    context_options.thread_count = 1;
    leo2_context* context = NULL;
    const leo2_result context_result = leo2_context_create(
        &context_options, &context);
    if (context_result == LEO2_UNSUPPORTED)
        return;
    require_result(context_result,
        "explicit AVX2 half-tail-column context");

    static const unsigned k = 65;
    static const unsigned r = 65;
    static const size_t bytes = 8193;
    CodecOwner* owner = make_codec(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    Shards original(k, Bytes(bytes, 0));
    original[k - 1] = random_shards(
        1, bytes, UINT64_C(0x48414c465441494c))[0];
    const Shards original_before = original;
    const ProfileLayout layout = leopard2_test::make_profile_layout(
        leopard2_test::kLegacyHigh, k, r);
    const Matrix generator =
        leopard2_test::direct_systematic_generator(gf8, layout);
    const Shards expected =
        oracle_parity(gf8, generator, original, r, LEO2_FIELD_GF8);

    const std::vector<uint8_t> all(r, 1);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    const EncodeResult dense = encode(owner->codec,
        LEO2_TEST_ENCODE_FORCE_TRANSFORM, original, all);
    require_result(dense.result, "GF8 half-tail-column dense encode");
    compare_requested(dense.recovery, expected, all, 0xa5,
        "GF8 half-tail-column/direct generator", counts);
    require(original == original_before,
        "GF8 half-tail-column encode modified caller input");
    require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                half_tail_column_calls != 0,
        "GF8 dense half-tail-column route was not exercised");

    // A sparse public output mask deliberately falls back to the mature
    // pruned transform; it must remain wire-identical to the same oracle.
    std::vector<uint8_t> sparse(r, 0);
    sparse[0] = 1;
    sparse[r / 2] = 1;
    sparse[r - 1] = 1;
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    const EncodeResult partial = encode(owner->codec,
        LEO2_TEST_ENCODE_FORCE_TRANSFORM, original, sparse);
    require_result(partial.result, "GF8 half-tail-column sparse encode");
    compare_requested(partial.recovery, expected, sparse, 0xa5,
        "GF8 sparse half-tail-column/direct generator", counts);
    require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                half_tail_column_calls == 0,
        "GF8 sparse output used the dense half-tail-column path");

    static const size_t input_prefix = 3;
    static const size_t output_prefix = 5;
    static const size_t suffix = 11;
    static const uint8_t canary = 0xd7;
    Shards input_storage(
        k, Bytes(input_prefix + bytes + suffix, canary));
    Shards output_storage(
        r, Bytes(output_prefix + bytes + suffix, canary));
    std::vector<const void*> input(k, NULL);
    std::vector<void*> output(r, NULL);
    for (unsigned i = 0; i < k; ++i)
    {
        memcpy(&input_storage[i][input_prefix], &original[i][0], bytes);
        input[i] = &input_storage[i][input_prefix];
    }
    for (unsigned i = 0; i < r; ++i)
        output[i] = &output_storage[i][output_prefix];
    size_t scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(
        owner->codec, bytes, &scratch_bytes),
        "GF8 half-tail-column unaligned scratch query");
    AlignedBuffer scratch(scratch_bytes);
    require_result(leo2_test_codec_set_encode_mode(owner->codec,
        LEO2_TEST_ENCODE_FORCE_TRANSFORM),
        "force GF8 half-tail-column unaligned transform");
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    require_result(leo2_encode(owner->codec, bytes, &input[0], &output[0],
        scratch.data(), scratch.size()),
        "GF8 half-tail-column unaligned encode");
    require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                half_tail_column_calls != 0,
        "unaligned GF8 half-tail-column route was not exercised");
    for (unsigned i = 0; i < r; ++i)
    {
        require(memcmp(&output_storage[i][output_prefix],
                    &expected[i][0], bytes) == 0,
            "unaligned GF8 half-tail-column differs from oracle");
        require(std::all_of(output_storage[i].begin(),
                    output_storage[i].begin() + output_prefix,
                    [](uint8_t value) { return value == canary; }) &&
                std::all_of(
                    output_storage[i].begin() + output_prefix + bytes,
                    output_storage[i].end(),
                    [](uint8_t value) { return value == canary; }),
            "unaligned GF8 half-tail-column changed an output guard");
    }
    for (unsigned i = 0; i < k; ++i)
    {
        require(std::all_of(input_storage[i].begin(),
                    input_storage[i].begin() + input_prefix,
                    [](uint8_t value) { return value == canary; }) &&
                std::all_of(
                    input_storage[i].begin() + input_prefix + bytes,
                    input_storage[i].end(),
                    [](uint8_t value) { return value == canary; }),
            "unaligned GF8 half-tail-column changed an input guard");
    }
    counts->unaligned_checks += k + r;

    ++counts->high_tail_column_checks;
    delete owner;
    leo2_context_destroy(context);
}
#endif

void test_high_gf8_byte_tiling(
    const BinaryField& gf8,
    Counts* counts)
{
    leo2_context_options avx2_options = {};
    avx2_options.struct_size = sizeof(avx2_options);
    avx2_options.backend = LEO2_BACKEND_AVX2;
    avx2_options.thread_count = 1;
    leo2_context* avx2_context = NULL;
    const leo2_result avx2_result = leo2_context_create(
        &avx2_options, &avx2_context);
    if (avx2_result == LEO2_UNSUPPORTED)
        return;
    require_result(avx2_result, "GF8 byte-tile AVX2 context create");

    leo2_context_options scalar_options = {};
    scalar_options.struct_size = sizeof(scalar_options);
    scalar_options.backend = LEO2_BACKEND_SCALAR;
    scalar_options.thread_count = 1;
    leo2_context* scalar_context = NULL;
    require_result(leo2_context_create(&scalar_options, &scalar_context),
        "GF8 byte-tile scalar context create");

    const unsigned k = 65;
    const unsigned r = 64;
    const size_t threshold_bytes = 256U * 1024U;
    const size_t below_bytes = threshold_bytes - 64U;
    const size_t tail_bytes = threshold_bytes + 17U;
    CodecOwner* target = make_codec(avx2_context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    CodecOwner* scalar = make_codec(scalar_context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);

    size_t target_below_scratch = 0;
    size_t scalar_below_scratch = 0;
    size_t target_scratch = 0;
    size_t scalar_scratch = 0;
    size_t target_tail_scratch = 0;
    size_t scalar_tail_scratch = 0;
    require_result(leo2_encode_scratch_size(
        target->codec, below_bytes, &target_below_scratch),
        "GF8 byte-tile below scratch query");
    require_result(leo2_encode_scratch_size(
        scalar->codec, below_bytes, &scalar_below_scratch),
        "GF8 byte-tile scalar below scratch query");
    require_result(leo2_encode_scratch_size(
        target->codec, threshold_bytes, &target_scratch),
        "GF8 byte-tile threshold scratch query");
    require_result(leo2_encode_scratch_size(
        scalar->codec, threshold_bytes, &scalar_scratch),
        "GF8 byte-tile scalar threshold scratch query");
    require_result(leo2_encode_scratch_size(
        target->codec, tail_bytes, &target_tail_scratch),
        "GF8 byte-tile tail scratch query");
    require_result(leo2_encode_scratch_size(
        scalar->codec, tail_bytes, &scalar_tail_scratch),
        "GF8 byte-tile scalar tail scratch query");
    require(target_below_scratch == scalar_below_scratch,
        "GF8 byte tiling crossed its measured lower boundary");
    require(target_scratch < scalar_scratch &&
            target_tail_scratch < scalar_tail_scratch,
        "GF8 byte tiling did not reduce transform scratch");

    // A one-block transform has a higher measured crossover.  Lock both sides
    // of the K>T policy boundary without allocating either large scratch area.
    CodecOwner* one_block = make_codec(avx2_context, 64, 64,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    CodecOwner* one_block_scalar = make_codec(scalar_context, 64, 64,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    size_t one_block_small = 0;
    size_t one_block_small_scalar = 0;
    size_t one_block_large = 0;
    size_t one_block_large_scalar = 0;
    require_result(leo2_encode_scratch_size(one_block->codec,
        threshold_bytes, &one_block_small),
        "GF8 byte-tile one-block small query");
    require_result(leo2_encode_scratch_size(one_block_scalar->codec,
        threshold_bytes, &one_block_small_scalar),
        "GF8 byte-tile scalar one-block small query");
    require_result(leo2_encode_scratch_size(one_block->codec,
        1024U * 1024U, &one_block_large),
        "GF8 byte-tile one-block large query");
    require_result(leo2_encode_scratch_size(one_block_scalar->codec,
        1024U * 1024U, &one_block_large_scalar),
        "GF8 byte-tile scalar one-block large query");
    require(one_block_small == one_block_small_scalar &&
            one_block_large < one_block_large_scalar,
        "GF8 byte-tile one-block crossover mismatch");
    delete one_block;
    delete one_block_scalar;

    // T=8 is a measured negative region even for large shards.
    CodecOwner* tiny_side = make_codec(avx2_context, 8, 8,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    CodecOwner* tiny_side_scalar = make_codec(scalar_context, 8, 8,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    size_t tiny_side_scratch = 0;
    size_t tiny_side_scalar_scratch = 0;
    require_result(leo2_encode_scratch_size(tiny_side->codec,
        1024U * 1024U, &tiny_side_scratch),
        "GF8 byte-tile T8 query");
    require_result(leo2_encode_scratch_size(tiny_side_scalar->codec,
        1024U * 1024U, &tiny_side_scalar_scratch),
        "GF8 byte-tile scalar T8 query");
    require(tiny_side_scratch == tiny_side_scalar_scratch,
        "GF8 byte tiling entered the losing T8 region");
    delete tiny_side;
    delete tiny_side_scalar;

    const Shards original = random_shards(
        k, tail_bytes, UINT64_C(0x474638425954494c));
    const Shards original_before = original;
    std::vector<uint8_t> dense_requested(r, 1);
    const EncodeResult reference = encode(scalar->codec,
        LEO2_TEST_ENCODE_FORCE_TRANSFORM, original, dense_requested);
    require_result(reference.result, "GF8 byte-tile scalar reference");
    const EncodeResult actual = encode(target->codec,
        LEO2_TEST_ENCODE_FORCE_TRANSFORM, original, dense_requested);
    require_result(actual.result, "GF8 byte-tile tail encode");
    compare_requested(actual.recovery, reference.recovery,
        dense_requested, 0xa5, "GF8 byte-tile/scalar", counts);
    require(original == original_before,
        "GF8 byte tiling modified a caller source shard");

    const ProfileLayout layout = leopard2_test::make_profile_layout(
        leopard2_test::kLegacyHigh, k, r);
    const Matrix generator =
        leopard2_test::direct_systematic_generator(gf8, layout);
    Shards first_byte(k, Bytes(1, 0));
    for (unsigned i = 0; i < k; ++i)
        first_byte[i][0] = original[i][0];
    const Shards direct = oracle_parity(
        gf8, generator, first_byte, r, LEO2_FIELD_GF8);
    for (unsigned parity = 0; parity < r; ++parity)
    {
        require(actual.recovery[parity][0] == direct[parity][0],
            "GF8 byte-tile output differs from direct algebra");
        ++counts->parity_symbols;
    }

    std::vector<uint8_t> sparse_requested(r, 0);
    sparse_requested[0] = 1;
    sparse_requested[r / 2] = 1;
    sparse_requested[r - 1] = 1;
    const EncodeResult sparse_reference = encode(scalar->codec,
        LEO2_TEST_ENCODE_FORCE_TRANSFORM, original, sparse_requested);
    const EncodeResult sparse_actual = encode(target->codec,
        LEO2_TEST_ENCODE_FORCE_TRANSFORM, original, sparse_requested);
    require_result(sparse_reference.result,
        "GF8 byte-tile sparse scalar encode");
    require_result(sparse_actual.result, "GF8 byte-tile sparse encode");
    compare_requested(sparse_actual.recovery, sparse_reference.recovery,
        sparse_requested, 0xa5, "GF8 byte-tile sparse/scalar", counts);

    const std::vector<const void*> input = const_pointers(original);
    Shards audited_output(r, Bytes(tail_bytes, 0));
    std::vector<void*> audited_recovery(r, NULL);
    for (unsigned i = 0; i < r; ++i)
        audited_recovery[i] = &audited_output[i][0];
    AlignedBuffer audited_scratch(target_tail_scratch);
    begin_allocation_audit();
    const leo2_result audited_result = leo2_encode(target->codec, tail_bytes,
        &input[0], &audited_recovery[0], audited_scratch.data(),
        audited_scratch.size());
    const uint64_t hot_allocations = end_allocation_audit();
    require_result(audited_result, "allocation-trapped GF8 byte tiling");
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
    require(hot_allocations == 0,
        "GF8 byte tiling allocated C++ storage");
    ++counts->allocation_checks;
#else
    (void)hot_allocations;
#endif
    require(audited_output == actual.recovery,
        "allocation-trapped GF8 byte tiling changed parity");

    ++counts->high_byte_tiling_checks;
    delete target;
    delete scalar;
    leo2_context_destroy(avx2_context);
    leo2_context_destroy(scalar_context);
}

void test_gf8_high_two_tail_generator_columns(
    const BinaryField& gf8,
    Counts* counts)
{
    leo2_context_options context_options = {};
    context_options.struct_size = sizeof(context_options);
    context_options.backend = LEO2_BACKEND_AVX2;
    context_options.thread_count = 1;
    leo2_context* context = NULL;
    const leo2_result context_result = leo2_context_create(
        &context_options, &context);
    if (context_result == LEO2_UNSUPPORTED)
        return;
    require_result(context_result, "explicit AVX2 two-tail context");

    static const size_t byte_counts[] = {
        65, 4097, 16383, 16384, 16385
    };
    static const unsigned t = 32;
    static const unsigned k = t + 2;
    CodecOwner* owner = make_codec(context, k, t,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    const ProfileLayout layout = leopard2_test::make_profile_layout(
        leopard2_test::kLegacyHigh, k, t);
    std::vector<Element> systematic_points(
        layout.systematic_coordinates.begin(),
        layout.systematic_coordinates.end());
    const std::vector<uint8_t> all(t, 1);

    for (size_t count_i = 0;
         count_i < sizeof(byte_counts) / sizeof(byte_counts[0]); ++count_i)
    {
        const size_t bytes = byte_counts[count_i];
        Shards original(k, Bytes(bytes, 0));
        const Shards tails = random_shards(
            2, bytes, UINT64_C(0x325441494c000000) + bytes);
        original[t] = tails[0];
        original[t + 1] = tails[1];
        const Shards original_before = original;
        Shards expected(t, Bytes(bytes, 0));
        for (unsigned tail = 0; tail < 2; ++tail)
        {
            const Element source_coordinate = systematic_points[t + tail];
            Element derivative = 1;
            for (size_t point = 0; point < systematic_points.size(); ++point)
            {
                if (systematic_points[point] != source_coordinate)
                {
                    derivative = gf8.multiply(derivative,
                        gf8.add(systematic_points[point], source_coordinate));
                }
            }
            require(derivative != 0,
                "GF8 two-tail generator derivative was zero");
            for (unsigned parity = 0; parity < t; ++parity)
            {
                const Element parity_coordinate = static_cast<Element>(
                    layout.parity_coordinates[parity]);
                Element vanishing = 1;
                for (size_t point = 0; point < systematic_points.size();
                     ++point)
                {
                    vanishing = gf8.multiply(vanishing,
                        gf8.add(parity_coordinate,
                            systematic_points[point]));
                }
                const Element denominator = gf8.multiply(derivative,
                    gf8.add(parity_coordinate, source_coordinate));
                const Element coefficient = gf8.divide(
                    vanishing, denominator);
                require(coefficient != 0,
                    "GF8 two-tail generator column contained zero");
                require(coefficient != 1,
                    "GF8 two-tail generator column contained one");
                for (size_t offset = 0; offset < bytes; ++offset)
                {
                    expected[parity][offset] ^= static_cast<uint8_t>(
                        gf8.multiply(coefficient,
                            original[t + tail][offset]));
                }
            }
        }

        leopard::ff8::TestOnlyResetHighEncodeCounts();
        const EncodeResult dense = encode(owner->codec,
            LEO2_TEST_ENCODE_FORCE_TRANSFORM, original, all);
        require_result(dense.result, "GF8 two-tail dense encode");
        compare_requested(dense.recovery, expected, all, 0xa5,
            "GF8 two-tail/direct generator", counts);
        const bool expected_direct = bytes >= 16384;
        require((leopard::ff8::TestOnlyGetHighEncodeCounts().
                    tail_column_calls != 0) == expected_direct,
            "GF8 two-tail route selected at the wrong byte boundary");
        require(original == original_before,
            "GF8 two-tail encode modified caller input");

        if (bytes == 16385)
        {
            std::vector<uint8_t> sparse(t, 0);
            sparse[0] = 1;
            sparse[15] = 1;
            sparse[31] = 1;
            leopard::ff8::TestOnlyResetHighEncodeCounts();
            const EncodeResult partial = encode(owner->codec,
                LEO2_TEST_ENCODE_FORCE_TRANSFORM, original, sparse);
            require_result(partial.result, "GF8 two-tail sparse encode");
            compare_requested(partial.recovery, expected, sparse, 0xa5,
                "GF8 sparse two-tail/direct generator", counts);
            require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                        tail_column_calls == 0,
                "GF8 sparse output used the dense two-tail path");
        }
        ++counts->high_tail_column_checks;
    }
    delete owner;
    leo2_context_destroy(context);
}

void test_gf8_high_coarse_direct_oracle(
    const BinaryField& gf8,
    Counts* counts)
{
    leo2_context_options context_options = {};
    context_options.struct_size = sizeof(context_options);
    context_options.backend = LEO2_BACKEND_AVX512;
    context_options.thread_count = 1;
    leo2_context* context = NULL;
    const leo2_result context_result = leo2_context_create(
        &context_options, &context);
    if (context_result == LEO2_UNSUPPORTED)
        return;
    require_result(context_result, "explicit AVX-512 oracle context");

    struct Case
    {
        unsigned k;
        unsigned r;
    };
    const Case cases[] = {
        { 8, 8 }, { 7, 8 }, { 8, 7 }, { 7, 7 },
        { 15, 16 }, { 16, 15 }, { 15, 15 },
        { 31, 32 }, { 32, 31 }, { 31, 31 },
        { 63, 64 }, { 64, 63 }, { 63, 63 }
    };
    const size_t bytes = 2049;
    const size_t oracle_bytes = 65;
    for (unsigned case_i = 0;
         case_i < sizeof(cases) / sizeof(cases[0]); ++case_i)
    {
        const Case& c = cases[case_i];
        CodecOwner* owner = make_codec(context, c.k, c.r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        const Shards original = random_shards(
            c.k, bytes, UINT64_C(0x434f415253450000) + case_i);
        const Shards original_before = original;
        const std::vector<uint8_t> requested(c.r, 1);

        leopard::ff8::TestOnlyResetHighEncodeCounts();
        const EncodeResult candidate = encode(owner->codec,
            LEO2_TEST_ENCODE_AUTO, original, requested);
        require_result(candidate.result, "GF8 coarse direct-oracle encode");
        require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                    whole_transform_calls == 1,
            "GF8 coarse direct-oracle encode missed the aligned callback");
        require(original == original_before,
            "GF8 coarse direct-oracle encode modified a source shard");

        const ProfileLayout layout = leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh, c.k, c.r);
        const Matrix generator =
            leopard2_test::direct_systematic_generator(gf8, layout);
        Shards head_original(c.k, Bytes(oracle_bytes));
        Shards tail_original(c.k, Bytes(oracle_bytes));
        for (unsigned original_i = 0; original_i < c.k; ++original_i)
        {
            std::copy(original[original_i].begin(),
                original[original_i].begin() + oracle_bytes,
                head_original[original_i].begin());
            std::copy(original[original_i].end() - oracle_bytes,
                original[original_i].end(), tail_original[original_i].begin());
        }
        const Shards head_expected = oracle_parity(
            gf8, generator, head_original, c.r, LEO2_FIELD_GF8);
        const Shards tail_expected = oracle_parity(
            gf8, generator, tail_original, c.r, LEO2_FIELD_GF8);
        Shards head_actual(c.r, Bytes(oracle_bytes));
        Shards tail_actual(c.r, Bytes(oracle_bytes));
        for (unsigned recovery = 0; recovery < c.r; ++recovery)
        {
            std::copy(candidate.recovery[recovery].begin(),
                candidate.recovery[recovery].begin() + oracle_bytes,
                head_actual[recovery].begin());
            std::copy(candidate.recovery[recovery].end() - oracle_bytes,
                candidate.recovery[recovery].end(),
                tail_actual[recovery].begin());
        }
        compare_requested(head_actual, head_expected, requested, 0xa5,
            "GF8 coarse head/direct oracle", counts);
        compare_requested(tail_actual, tail_expected, requested, 0xa5,
            "GF8 coarse tail/direct oracle", counts);
        ++counts->gf8_coarse_oracle_checks;
        delete owner;
    }
    leo2_context_destroy(context);
}

void test_auto_dispatch_threshold(leo2_context* context, Counts* counts)
{
    const leo2_backend backend = leo2_context_backend(context);
    const bool qualified_simd = backend == LEO2_BACKEND_SSSE3 ||
        backend == LEO2_BACKEND_AVX2 ||
        backend == LEO2_BACKEND_AVX512 ||
        backend == LEO2_BACKEND_GFNI;
    const bool qualified_k7 = qualified_simd || backend == LEO2_BACKEND_SCALAR;

    CodecOwner* low = make_codec(context, 7, 7,
        LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8);
    int direct = -1;
    require_result(leo2_test_codec_encode_path(low->codec, 1024, 1, &direct),
        "AUTO qualified path query");
    require(direct == (qualified_k7 ? 1 : 0),
        "AUTO qualified region disagrees with the runtime backend");
    ++counts->dispatch_checks;
    require_result(leo2_test_codec_encode_path(low->codec, 1023, 1, &direct),
        "AUTO byte-neighbor path query");
    require(direct == 0, "AUTO selected a ragged byte neighbor");
    ++counts->dispatch_checks;
    require_result(leo2_test_codec_encode_path(low->codec, 64, 1, &direct),
        "AUTO size-threshold path query");
    require(direct == 0, "AUTO selected below its measured size threshold");
    ++counts->dispatch_checks;
    require_result(leo2_test_codec_encode_path(low->codec, 1024, 2, &direct),
        "AUTO Q-neighbor path query");
    require(direct == 0, "AUTO selected an unpromoted Q neighbor");
    ++counts->dispatch_checks;
    const Shards low_original = random_shards(7, 1024, UINT64_C(0xa07010ff));
    std::vector<uint8_t> singleton(7, 0);
    singleton[6] = 1;
    const EncodeResult low_automatic = encode(low->codec,
        LEO2_TEST_ENCODE_AUTO, low_original, singleton);
    const EncodeResult low_direct = encode(low->codec,
        LEO2_TEST_ENCODE_FORCE_DIRECT, low_original, singleton);
    const EncodeResult low_transform = encode(low->codec,
        LEO2_TEST_ENCODE_FORCE_TRANSFORM, low_original, singleton);
    require_result(low_automatic.result, "AUTO qualified encode");
    require(low_automatic.recovery == low_direct.recovery &&
            low_automatic.recovery == low_transform.recovery,
        "AUTO qualified output differs between direct and transform");
    ++counts->dispatch_checks;
    delete low;

    CodecOwner* k2 = make_codec(context, 2, 7,
        LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8);
    require_result(leo2_test_codec_encode_path(k2->codec, 1024, 1, &direct),
        "AUTO K=2 boundary path query");
    require(direct == (qualified_simd ? 1 : 0),
        "AUTO K=2 boundary disagrees with the runtime backend");
    ++counts->dispatch_checks;
    delete k2;

    CodecOwner* k3 = make_codec(context, 3, 7,
        LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8);
    require_result(leo2_test_codec_encode_path(k3->codec, 960, 1, &direct),
        "AUTO scalar size boundary path query");
    require(direct == 0, "AUTO selected below its 1-KiB threshold");
    ++counts->dispatch_checks;
    require_result(leo2_test_codec_encode_path(k3->codec, 1024, 1, &direct),
        "AUTO scalar K=3 boundary path query");
    require(direct == (qualified_k7 ? 1 : 0),
        "AUTO scalar K=3 boundary disagrees with the runtime backend");
    ++counts->dispatch_checks;
    require_result(leo2_test_codec_encode_path(k3->codec, 1088, 1, &direct),
        "AUTO regular-byte neighbor path query");
    require(direct == (qualified_k7 ? 1 : 0),
        "AUTO rejected the measured regular-byte neighbor");
    ++counts->dispatch_checks;
    delete k3;

    CodecOwner* high = make_codec(context, 7, 7,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    require_result(leo2_test_codec_encode_path(high->codec, 1024, 1, &direct),
        "AUTO high-profile path query");
#if defined(LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE) && \
    LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE_AUTO
    const bool measured_high_backend = backend == LEO2_BACKEND_AVX2;
#else
    const bool measured_high_backend = false;
#endif
    require(direct == (measured_high_backend ? 1 : 0),
        "AUTO high-profile selection escaped its measured AVX2 region");
    ++counts->dispatch_checks;
    delete high;

    CodecOwner* historical = make_codec(context, 2, 16,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    require_result(leo2_test_codec_encode_path(
        historical->codec, 4096, 1, &direct),
        "historical high-profile path query");
    require(direct == (measured_high_backend ? 1 : 0),
        "historical high-profile path escaped its default-off AVX2 gate");
    ++counts->dispatch_checks;
    const Shards historical_original = random_shards(
        2, 4096, UINT64_C(0x484947484b325231));
    for (unsigned recovery_index = 0; recovery_index < 16;
         recovery_index += 15)
    {
        std::vector<uint8_t> requested(16, 0);
        requested[recovery_index] = 1;
        const EncodeResult automatic = encode(historical->codec,
            LEO2_TEST_ENCODE_AUTO, historical_original, requested);
        const EncodeResult forced_direct = encode(historical->codec,
            LEO2_TEST_ENCODE_FORCE_DIRECT, historical_original, requested);
        const EncodeResult forced_transform = encode(historical->codec,
            LEO2_TEST_ENCODE_FORCE_TRANSFORM, historical_original, requested);
        require_result(automatic.result, "historical high AUTO encode");
        require_result(forced_direct.result,
            "historical high force-direct encode");
        require_result(forced_transform.result,
            "historical high force-transform encode");
        require(automatic.recovery == forced_direct.recovery &&
                automatic.recovery == forced_transform.recovery,
            "historical high AUTO/direct/transform parity mismatch");
        ++counts->dispatch_checks;
    }
    delete historical;

#if defined(LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE) && \
    LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE_AUTO
    leo2_context_options avx2_options = {};
    avx2_options.struct_size = sizeof(avx2_options);
    avx2_options.backend = LEO2_BACKEND_AVX2;
    avx2_options.thread_count = 1;
    leo2_context* avx2_context = NULL;
    const leo2_result avx2_created =
        leo2_context_create(&avx2_options, &avx2_context);
    if (avx2_created != LEO2_UNSUPPORTED)
    {
        require_result(avx2_created,
            "high-selector bounded AVX2 context create");
        for (unsigned k = 2; k <= 16; ++k)
            for (unsigned r = 2; r <= 16; ++r)
            {
                CodecOwner* bounded = make_codec(avx2_context, k, r,
                    LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
                require_result(leo2_test_codec_encode_path(
                    bounded->codec, 4096, 1, &direct),
                    "high-selector bounded path query");
                require(direct == 1,
                    "high-selector rejected a bounded diagnostic AVX2 shape");
                ++counts->dispatch_checks;
                delete bounded;
            }
        leo2_context_destroy(avx2_context);
    }
#endif

    CodecOwner* single_side = make_codec(context, 2, 1,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    require_result(leo2_test_codec_encode_path(
        single_side->codec, 4096, 1, &direct),
        "AUTO high-profile R=1 path query");
    require(direct == 0,
        "AUTO path query bypassed the existing single-side encoder");
    ++counts->dispatch_checks;
    direct = -1;
    require(leo2_test_codec_encode_path(
            single_side->codec, 0, 1, &direct) == LEO2_INVALID_ARGUMENT &&
            direct == 0,
        "AUTO high-profile R=1 path query accepted zero shard bytes");
    ++counts->dispatch_checks;
    direct = -1;
    require(leo2_test_codec_encode_path(
            single_side->codec, UINT64_MAX, 1, &direct) ==
                LEO2_INVALID_ARGUMENT &&
            direct == 0,
        "AUTO high-profile R=1 path query accepted overflowing shard bytes");
    ++counts->dispatch_checks;
    require_result(leo2_test_codec_set_encode_mode(
        single_side->codec, LEO2_TEST_ENCODE_FORCE_DIRECT),
        "force direct for high-profile R=1 control");
    require_result(leo2_test_codec_encode_path(
        single_side->codec, 4096, 1, &direct),
        "forced high-profile R=1 path query");
    require(direct == 1,
        "forced high-profile R=1 control did not select generator-row direct");
    ++counts->dispatch_checks;
    delete single_side;

    CodecOwner* single_side_gf16 = make_codec(context, 2, 1,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16);
    direct = -1;
    require(leo2_test_codec_encode_path(
            single_side_gf16->codec, 1, 1, &direct) == LEO2_UNSUPPORTED &&
            direct == 0,
        "AUTO high-profile GF16 R=1 path query accepted an odd shard");
    ++counts->dispatch_checks;
    delete single_side_gf16;

    CodecOwner* high_gf16 = make_codec(context, 7, 7,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16);
    require_result(leo2_test_codec_encode_path(
        high_gf16->codec, 1024, 1, &direct),
        "AUTO high-profile GF16 path query");
    require(direct == 0,
        "experimental high selector admitted unmeasured GF16");
    ++counts->dispatch_checks;
    delete high_gf16;

    const leo2_backend high_control_backends[] = {
        LEO2_BACKEND_SCALAR,
        LEO2_BACKEND_SSSE3,
        LEO2_BACKEND_AVX2,
        LEO2_BACKEND_AVX512,
        LEO2_BACKEND_GFNI,
        LEO2_BACKEND_NEON
    };
    for (size_t i = 0;
         i < sizeof(high_control_backends) /
             sizeof(high_control_backends[0]); ++i)
    {
        leo2_context_options fixed_options = {};
        fixed_options.struct_size = sizeof(fixed_options);
        fixed_options.backend = high_control_backends[i];
        fixed_options.thread_count = 1;
        leo2_context* fixed_context = NULL;
        const leo2_result created =
            leo2_context_create(&fixed_options, &fixed_context);
        if (created == LEO2_UNSUPPORTED)
            continue;
        require_result(created, "high-selector control context create");
        CodecOwner* control = make_codec(fixed_context, 7, 7,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        require_result(leo2_test_codec_encode_path(
            control->codec, 1024, 1, &direct),
            "high-selector backend control path query");
#if defined(LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE) && \
    LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE_AUTO
        const bool expected_high =
            high_control_backends[i] == LEO2_BACKEND_AVX2;
#else
        const bool expected_high = false;
#endif
        require(direct == (expected_high ? 1 : 0),
            "experimental high selector admitted an unmeasured backend");
        ++counts->dispatch_checks;
        delete control;
        leo2_context_destroy(fixed_context);
    }

    CodecOwner* identity = make_codec(context, 1, 7,
        LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8);
    require_result(leo2_test_codec_encode_path(identity->codec, 1024, 1, &direct),
        "AUTO identity path query");
    require(direct == 0, "AUTO replaced the existing K=1 copy path");
    ++counts->dispatch_checks;
    delete identity;

    CodecOwner* gf16 = make_codec(context, 7, 7,
        LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16);
    require_result(leo2_test_codec_encode_path(gf16->codec, 1024, 1, &direct),
        "AUTO GF16 qualified path query");
    require(direct == (qualified_k7 ? 1 : 0),
        "AUTO GF16 region disagrees with the runtime backend");
    ++counts->dispatch_checks;
    const Shards gf16_original = random_shards(7, 1024,
        UINT64_C(0xa07016ff));
    const EncodeResult gf16_automatic = encode(gf16->codec,
        LEO2_TEST_ENCODE_AUTO, gf16_original, singleton);
    const EncodeResult gf16_direct = encode(gf16->codec,
        LEO2_TEST_ENCODE_FORCE_DIRECT, gf16_original, singleton);
    const EncodeResult gf16_transform = encode(gf16->codec,
        LEO2_TEST_ENCODE_FORCE_TRANSFORM, gf16_original, singleton);
    require_result(gf16_automatic.result, "AUTO GF16 qualified encode");
    require(gf16_automatic.recovery == gf16_direct.recovery &&
            gf16_automatic.recovery == gf16_transform.recovery,
        "AUTO GF16 output differs between direct and transform");
    ++counts->dispatch_checks;
    delete gf16;
}

void test_k1_binding_respects_encode_mode(Counts* counts)
{
    leo2_context_options options = {};
    options.struct_size = sizeof(options);
    options.backend = LEO2_BACKEND_AVX2;
    options.thread_count = 1;
    leo2_context* context = NULL;
    const leo2_result created = leo2_context_create(&options, &context);
    if (created == LEO2_UNSUPPORTED)
        return;
    require_result(created, "K1 binding AVX2 context create");

    CodecOwner* owner = make_codec(context, 1, 1,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    static const uint64_t bytes = 1024;
    const Shards original = random_shards(
        1, static_cast<size_t>(bytes), UINT64_C(0x4b3142494e44494e));
    const void* original_pointer[1] = { &original[0][0] };
    Bytes recovery(static_cast<size_t>(bytes), 0);
    void* recovery_pointer[1] = { &recovery[0] };

    require_result(leo2_test_codec_set_encode_mode(
        owner->codec, LEO2_TEST_ENCODE_FORCE_TRANSFORM),
        "force K1 binding transform mode");
    size_t scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(
        owner->codec, bytes, &scratch_bytes),
        "forced K1 binding scratch query");
    AlignedBuffer scratch(scratch_bytes);
    leo2_encode_batch_item item = {
        bytes, original_pointer, recovery_pointer,
        scratch.data(), scratch_bytes
    };
    leo2_encode_batch_binding* binding = NULL;
    require_result(leo2_encode_batch_binding_create(
        owner->codec, &item, 1, &binding),
        "forced K1 transform binding create");
    require(leo2_test_encode_batch_binding_uses_k1_copy(binding) == 0,
        "forced K1 transform binding captured the copy terminal");
    require_result(leo2_encode_batch_binding_execute(binding),
        "forced K1 transform binding execute");
    require(recovery == original[0],
        "forced K1 transform binding parity mismatch");
    leo2_encode_batch_binding_destroy(binding);
    ++counts->dispatch_checks;

    require_result(leo2_test_codec_set_encode_mode(
        owner->codec, LEO2_TEST_ENCODE_AUTO),
        "restore K1 binding AUTO mode");
    item.scratch = NULL;
    item.scratch_bytes = 0;
    binding = NULL;
    require_result(leo2_encode_batch_binding_create(
        owner->codec, &item, 1, &binding),
        "AUTO K1 copy binding create");
    require(leo2_test_encode_batch_binding_uses_k1_copy(binding) == 1,
        "AUTO K1 binding did not capture the copy terminal");
    leo2_encode_batch_binding_destroy(binding);
    ++counts->dispatch_checks;

    delete owner;
    leo2_context_destroy(context);
}

void test_tail_allocation_and_contracts(
    leo2_context* context,
    const BinaryField& gf8,
    const BinaryField& gf16,
    Counts* counts)
{
    struct Case {
        leo2_profile profile;
        leo2_field field;
        leo2_shard_layout layout;
        unsigned k;
        unsigned r;
        size_t bytes;
    };
    const Case cases[] = {
        { LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          LEO2_SHARD_LAYOUT_NATIVE_V1, 9, 7, 65 },
        { LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8,
          LEO2_SHARD_LAYOUT_NATIVE_V1, 5, 11, 129 },
        { LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16,
          LEO2_SHARD_LAYOUT_NATIVE_V1, 9, 7, 66 },
        { LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16,
          LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1, 5, 11, 66 },
        { LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          LEO2_SHARD_LAYOUT_NATIVE_V1, 7, 5, 4095 },
        { LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8,
          LEO2_SHARD_LAYOUT_NATIVE_V1, 5, 7, 4096 },
        { LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8,
          LEO2_SHARD_LAYOUT_NATIVE_V1, 3, 5, 4097 },
        { LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16,
          LEO2_SHARD_LAYOUT_NATIVE_V1, 7, 5, 4094 },
        { LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16,
          LEO2_SHARD_LAYOUT_NATIVE_V1, 5, 7, 4096 },
        { LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16,
          LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1, 3, 5, 4098 }
    };
    for (unsigned case_i = 0; case_i < sizeof(cases) / sizeof(cases[0]); ++case_i)
    {
        const Case& c = cases[case_i];
        CodecOwner* owner = make_codec(
            context, c.k, c.r, c.profile, c.field, c.layout);
        require_result(leo2_test_codec_set_encode_mode(owner->codec,
            LEO2_TEST_ENCODE_FORCE_DIRECT), "force direct tail");
        Shards original = random_shards(c.k, c.bytes,
            UINT64_C(0x7461696c5f636173) + case_i);
        if (c.layout == LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1)
            for (unsigned i = 0; i < c.k; ++i)
                original[i][c.bytes - 1] = 0;
        const BinaryField& oracle_field = c.field == LEO2_FIELD_GF8 ? gf8 : gf16;
        const ProfileLayout layout = leopard2_test::make_profile_layout(
            oracle_profile(c.profile), c.k, c.r);
        const Matrix generator =
            leopard2_test::direct_systematic_generator(oracle_field, layout);
        const Shards expected = oracle_parity(
            oracle_field, generator, original, c.r, c.field);
        std::vector<const void*> input = const_pointers(original);
        Shards output(c.r, Bytes(c.bytes, 0));
        std::vector<void*> recovery(c.r, NULL);
        for (unsigned i = 0; i < c.r; ++i)
            recovery[i] = &output[i][0];
        size_t scratch_bytes = 0;
        require_result(leo2_encode_scratch_size(
            owner->codec, c.bytes, &scratch_bytes), "tail scratch query");
        AlignedBuffer scratch(scratch_bytes);

        begin_allocation_audit();
        const leo2_result result = leo2_encode(owner->codec, c.bytes, &input[0],
            &recovery[0], scratch.data(), scratch.size());
        const uint64_t direct_allocations = end_allocation_audit();
        require_result(result, "allocation-trapped tail encode");
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
        require(direct_allocations == 0,
            "direct encode allocated C++ storage");
        ++counts->allocation_checks;
#else
        (void)direct_allocations;
#endif
        require(output == expected, "tail direct encode differs from oracle");

        require_result(leo2_test_codec_set_encode_mode(owner->codec,
            LEO2_TEST_ENCODE_FORCE_TRANSFORM), "force transform tail");
        output.assign(c.r, Bytes(c.bytes, 0xa5));
        for (unsigned i = 0; i < c.r; ++i)
            recovery[i] = &output[i][0];
        begin_allocation_audit();
        const leo2_result transform_result = leo2_encode(
            owner->codec, c.bytes, &input[0], &recovery[0],
            scratch.data(), scratch.size());
        const uint64_t transform_allocations = end_allocation_audit();
        require_result(transform_result, "allocation-trapped transform tail encode");
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
        require(transform_allocations == 0,
            "transform tail encode allocated C++ storage");
        ++counts->allocation_checks;
#else
        (void)transform_allocations;
#endif
        require(output == expected,
            "tail transform encode differs from direct generator oracle");

        std::vector<uint8_t> subset(c.r, 0);
        subset[0] = subset[c.r / 2] = subset[c.r - 1] = 1;
        output.assign(c.r, Bytes(c.bytes, 0xa5));
        recovery.assign(c.r, NULL);
        for (unsigned i = 0; i < c.r; ++i)
            if (subset[i])
                recovery[i] = &output[i][0];
        begin_allocation_audit();
        const leo2_result sparse_transform_result = leo2_encode(
            owner->codec, c.bytes, &input[0], &recovery[0],
            scratch.data(), scratch.size());
        const uint64_t sparse_transform_allocations = end_allocation_audit();
        require_result(sparse_transform_result,
            "allocation-trapped sparse transform encode");
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
        require(sparse_transform_allocations == 0,
            "sparse transform encode allocated C++ storage");
        ++counts->allocation_checks;
#else
        (void)sparse_transform_allocations;
#endif
        compare_requested(output, expected, subset, 0xa5,
            "allocation-trapped sparse transform/oracle", counts);

        const EncodeResult transform_subset = encode(
            owner->codec, LEO2_TEST_ENCODE_FORCE_TRANSFORM,
            original, subset);
        require_result(transform_subset.result,
            "requested transform tail encode");
        compare_requested(
            transform_subset.recovery, expected, subset, 0xa5,
            "requested transform tail/oracle", counts);

        std::vector<void*> one_output(c.r, NULL);
        one_output[0] = const_cast<void*>(input[0]);
        require(leo2_encode(owner->codec, c.bytes, &input[0], &one_output[0],
            scratch.data(), scratch.size()) == LEO2_OVERLAP,
            "direct encode accepted output/input aliasing");
        ++counts->contract_checks;
        if (c.r > 1)
        {
            one_output.assign(c.r, NULL);
            one_output[0] = &output[0][0];
            one_output[1] = &output[0][0];
            require(leo2_encode(owner->codec, c.bytes, &input[0], &one_output[0],
                scratch.data(), scratch.size()) == LEO2_OVERLAP,
                "direct encode accepted overlapping outputs");
            ++counts->contract_checks;
        }
        require(scratch.size() != 0, "encode scratch unexpectedly empty");
        require(leo2_encode(owner->codec, c.bytes, &input[0], &recovery[0],
            scratch.data(), scratch.size() - 1) == LEO2_SCRATCH_TOO_SMALL,
            "direct encode accepted undersized scratch");
        ++counts->contract_checks;

        if (c.layout == LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1)
        {
            original[0][c.bytes - 1] = 1;
            input = const_pointers(original);
            require(leo2_encode(owner->codec, c.bytes, &input[0], &recovery[0],
                scratch.data(), scratch.size()) == LEO2_INVALID_ARGUMENT,
                "padded-odd direct encode accepted a nonzero systematic pad");
            ++counts->contract_checks;
        }
        delete owner;
    }
}

void test_unaligned_guarded_buffers(
    leo2_context* context,
    const BinaryField& gf8,
    const BinaryField& gf16,
    Counts* counts)
{
    struct Case {
        leo2_profile profile;
        leo2_field field;
        unsigned k;
        unsigned r;
        size_t bytes;
        bool direct_first;
    };
    const Case cases[] = {
        { LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 5, 3, 65, true },
        { LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 252, 4, 65, false },
        { LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16, 5, 3, 66, true }
    };
    const size_t input_prefix = 3;
    const size_t output_prefix = 5;
    const size_t suffix = 11;
    const uint8_t canary = 0xd7;
    for (unsigned case_i = 0; case_i < sizeof(cases) / sizeof(cases[0]); ++case_i)
    {
        const Case& c = cases[case_i];
        const BinaryField& field = c.field == LEO2_FIELD_GF8 ? gf8 : gf16;
        CodecOwner* owner = make_codec(
            context, c.k, c.r, c.profile, c.field);
        if (c.direct_first)
        {
            require_result(leo2_test_codec_set_encode_mode(owner->codec,
                LEO2_TEST_ENCODE_FORCE_DIRECT), "force direct unaligned");
        }

        const Shards original = random_shards(c.k, c.bytes,
            UINT64_C(0x756e616c69676e) + case_i);
        const ProfileLayout layout = leopard2_test::make_profile_layout(
            oracle_profile(c.profile), c.k, c.r);
        const Matrix generator =
            leopard2_test::direct_systematic_generator(field, layout);
        const Shards expected = oracle_parity(
            field, generator, original, c.r, c.field);

        Shards input_storage(
            c.k, Bytes(input_prefix + c.bytes + suffix, canary));
        Shards output_storage(
            c.r, Bytes(output_prefix + c.bytes + suffix, canary));
        std::vector<const void*> input(c.k, NULL);
        std::vector<void*> output(c.r, NULL);
        for (unsigned i = 0; i < c.k; ++i)
        {
            memcpy(&input_storage[i][input_prefix], &original[i][0], c.bytes);
            input[i] = &input_storage[i][input_prefix];
            require((reinterpret_cast<uintptr_t>(input[i]) & 15u) != 0,
                "test input unexpectedly remained SIMD-aligned");
        }
        for (unsigned i = 0; i < c.r; ++i)
        {
            output[i] = &output_storage[i][output_prefix];
            require((reinterpret_cast<uintptr_t>(output[i]) & 15u) != 0,
                "test output unexpectedly remained SIMD-aligned");
        }

        size_t scratch_bytes = 0;
        require_result(leo2_encode_scratch_size(
            owner->codec, c.bytes, &scratch_bytes), "unaligned scratch query");
        AlignedBuffer scratch(scratch_bytes);
        require_result(leo2_encode(owner->codec, c.bytes, &input[0], &output[0],
            scratch.data(), scratch.size()), "unaligned initial encode");
        for (unsigned i = 0; i < c.r; ++i)
            require(memcmp(&output_storage[i][output_prefix],
                        &expected[i][0], c.bytes) == 0,
                "unaligned initial encode differs from the oracle");

        require_result(leo2_test_codec_set_encode_mode(owner->codec,
            LEO2_TEST_ENCODE_FORCE_TRANSFORM), "force transform unaligned");
        for (unsigned i = 0; i < c.r; ++i)
            std::fill(output_storage[i].begin() + output_prefix,
                output_storage[i].begin() + output_prefix + c.bytes, canary);
        require_result(leo2_encode(owner->codec, c.bytes, &input[0], &output[0],
            scratch.data(), scratch.size()), "unaligned transform encode");
        for (unsigned i = 0; i < c.r; ++i)
            require(memcmp(&output_storage[i][output_prefix],
                        &expected[i][0], c.bytes) == 0,
                "unaligned transform encode differs from the oracle");

        if (case_i == 0)
        {
            Shards aliased_original = original;
            aliased_original[1] = aliased_original[0];
            const Shards aliased_expected = oracle_parity(
                field, generator, aliased_original, c.r, c.field);
            input[1] = input[0];
            for (unsigned i = 0; i < c.r; ++i)
                std::fill(output_storage[i].begin() + output_prefix,
                    output_storage[i].begin() + output_prefix + c.bytes, canary);
            require_result(leo2_encode(owner->codec, c.bytes, &input[0],
                &output[0], scratch.data(), scratch.size()),
                "aliased-input transform encode");
            for (unsigned i = 0; i < c.r; ++i)
                require(memcmp(&output_storage[i][output_prefix],
                            &aliased_expected[i][0], c.bytes) == 0,
                    "allowed input alias changed transform parity");
            ++counts->contract_checks;
        }

        for (unsigned i = 0; i < c.k; ++i)
        {
            require(std::all_of(input_storage[i].begin(),
                        input_storage[i].begin() + input_prefix,
                        [](uint8_t value) { return value == canary; }) &&
                    std::all_of(input_storage[i].begin() + input_prefix + c.bytes,
                        input_storage[i].end(),
                        [](uint8_t value) { return value == canary; }),
                "unaligned encode changed an input guard");
        }
        for (unsigned i = 0; i < c.r; ++i)
        {
            require(std::all_of(output_storage[i].begin(),
                        output_storage[i].begin() + output_prefix,
                        [](uint8_t value) { return value == canary; }) &&
                    std::all_of(output_storage[i].begin() + output_prefix + c.bytes,
                        output_storage[i].end(),
                        [](uint8_t value) { return value == canary; }),
                "unaligned encode changed an output guard");
        }
        counts->unaligned_checks += c.k + c.r;
        delete owner;
    }
}

void test_auto_encode_batch(leo2_context* context, Counts* counts)
{
    const unsigned k = 3;
    const unsigned r = 3;
    const size_t bytes = 1024;
    const size_t batch_size = 4;
    const uint8_t sentinel = 0xa5;
    CodecOwner* owner = make_codec(context, k, r,
        LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8);
    size_t scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(owner->codec, bytes, &scratch_bytes),
        "batch scratch query");

    std::vector<Shards> originals(batch_size);
    std::vector<Shards> outputs(batch_size);
    std::vector<Shards> expected(batch_size);
    std::vector<std::vector<const void*> > inputs(batch_size);
    std::vector<std::vector<void*> > recoveries(batch_size);
    std::vector<AlignedBuffer*> scratches(batch_size, NULL);
    std::vector<leo2_encode_batch_item> items(batch_size);
    for (size_t item = 0; item < batch_size; ++item)
    {
        originals[item] = random_shards(k, bytes,
            UINT64_C(0xb47c0000) + item);
        std::vector<uint8_t> requested(r, 0);
        requested[item % r] = 1;
        const EncodeResult reference = encode(owner->codec,
            LEO2_TEST_ENCODE_FORCE_TRANSFORM, originals[item], requested,
            sentinel);
        require_result(reference.result, "batch transform reference");
        expected[item] = reference.recovery;

        outputs[item].assign(r, Bytes(bytes, sentinel));
        inputs[item] = const_pointers(originals[item]);
        recoveries[item].assign(r, NULL);
        recoveries[item][item % r] = &outputs[item][item % r][0];
        scratches[item] = new AlignedBuffer(scratch_bytes);
        items[item].shard_bytes = bytes;
        items[item].original = &inputs[item][0];
        items[item].recovery = &recoveries[item][0];
        items[item].scratch = scratches[item]->data();
        items[item].scratch_bytes = scratches[item]->size();
    }

    require_result(leo2_test_codec_set_encode_mode(
        owner->codec, LEO2_TEST_ENCODE_AUTO), "set batch AUTO mode");
    require_result(leo2_encode_batch(owner->codec, &items[0], items.size()),
        "AUTO encode batch");
    for (size_t item = 0; item < batch_size; ++item)
    {
        require(outputs[item] == expected[item],
            "AUTO batch output differs from transform reference");
        outputs[item].assign(r, Bytes(bytes, sentinel));
        recoveries[item][item % r] = &outputs[item][item % r][0];
    }
    require_result(leo2_test_codec_set_encode_mode(
        owner->codec, LEO2_TEST_ENCODE_FORCE_DIRECT),
        "set batch forced-direct mode");
    require_result(leo2_encode_batch(owner->codec, &items[0], items.size()),
        "forced-direct encode batch");
    for (size_t item = 0; item < batch_size; ++item)
    {
        require(outputs[item] == expected[item],
            "forced-direct batch output differs from transform reference");
        delete scratches[item];
    }
    counts->batch_executions += batch_size * 2;
    delete owner;
}

void test_sparse_low_gf16_auto_promotion(
    const BinaryField& gf16,
    Counts* counts)
{
    const unsigned k = 128;
    const unsigned r = 896;
    const size_t bytes = 1024;
    const size_t tail_bytes = 1026;
    static const unsigned edge_indices[] = {
        0, 127, 128, 383, 895
    };
    static const unsigned scattered_indices[] = {
        7, 63, 135, 255, 519, 895
    };

    std::vector<uint8_t> edge(r, 0);
    std::vector<uint8_t> scattered(r, 0);
    for (size_t i = 0;
         i < sizeof(edge_indices) / sizeof(edge_indices[0]); ++i)
        edge[edge_indices[i]] = 1;
    for (size_t i = 0;
         i < sizeof(scattered_indices) / sizeof(scattered_indices[0]); ++i)
        scattered[scattered_indices[i]] = 1;

    // Request AVX2 explicitly so an AVX-512 AUTO policy cannot accidentally
    // turn this into a skipped or misclassified test.  Unsupported hosts still
    // cover exact schedules through FORCE_TRANSFORM in the generic matrix.
    leo2_context_options target_options = {};
    target_options.struct_size = sizeof(target_options);
    target_options.backend = LEO2_BACKEND_AVX2;
    target_options.thread_count = 1;
    leo2_context* target_context = NULL;
    const leo2_result target_context_result = leo2_context_create(
        &target_options, &target_context);
    if (target_context_result == LEO2_UNSUPPORTED)
        return;
    require_result(target_context_result,
        "sparse AUTO AVX2 context create");

    leo2_context_options scalar_options = {};
    scalar_options.struct_size = sizeof(scalar_options);
    scalar_options.backend = LEO2_BACKEND_SCALAR;
    scalar_options.thread_count = 1;
    leo2_context* scalar_context = NULL;
    require_result(leo2_context_create(&scalar_options, &scalar_context),
        "sparse AUTO scalar context create");
    CodecOwner* target = make_codec(target_context, k, r,
        LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16);
    CodecOwner* scalar = make_codec(scalar_context, k, r,
        LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16);

    require_result(leo2_test_codec_set_encode_mode(
        target->codec, LEO2_TEST_ENCODE_AUTO), "sparse AUTO target mode");
    require_result(leo2_test_codec_set_encode_mode(
        scalar->codec, LEO2_TEST_ENCODE_AUTO), "sparse AUTO scalar mode");
    size_t target_scratch_bytes = 0;
    size_t scalar_scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(
        target->codec, bytes, &target_scratch_bytes),
        "sparse AUTO target scratch query");
    require_result(leo2_encode_scratch_size(
        scalar->codec, bytes, &scalar_scratch_bytes),
        "sparse AUTO scalar scratch query");
    require(target_scratch_bytes == scalar_scratch_bytes,
        "hook-build sparse scratch upper bound changed across backends");
    size_t target_small_scratch = 0;
    size_t scalar_small_scratch = 0;
    require_result(leo2_encode_scratch_size(
        target->codec, 64, &target_small_scratch),
        "sparse AUTO 64-byte target scratch query");
    require_result(leo2_encode_scratch_size(
        scalar->codec, 64, &scalar_small_scratch),
        "sparse AUTO 64-byte scalar scratch query");
    require(target_small_scratch == scalar_small_scratch,
        "64-byte negative cell reserved an exact schedule");

    const Shards original = random_shards(
        k, tail_bytes, UINT64_C(0x5350415253454155));
    Shards aligned_original(k, Bytes(bytes, 0));
    Shards small_original(k, Bytes(64, 0));
    for (unsigned i = 0; i < k; ++i)
    {
        std::copy(original[i].begin(), original[i].begin() + bytes,
            aligned_original[i].begin());
        std::copy(original[i].begin(), original[i].begin() + 64,
            small_original[i].begin());
    }

    const EncodeResult edge_reference = encode(
        scalar->codec, LEO2_TEST_ENCODE_AUTO, aligned_original, edge);
    const EncodeResult scattered_reference = encode(
        scalar->codec, LEO2_TEST_ENCODE_AUTO, aligned_original, scattered);
    require_result(edge_reference.result, "sparse AUTO edge scalar encode");
    require_result(scattered_reference.result,
        "sparse AUTO scattered scalar encode");

    leopard::ff16::TestOnlyResetSparseEncodeCounts();
    const EncodeResult edge_actual = encode(
        target->codec, LEO2_TEST_ENCODE_AUTO, aligned_original, edge);
    require_result(edge_actual.result, "sparse AUTO edge encode");
    compare_requested(edge_actual.recovery, edge_reference.recovery, edge,
        0xa5, "sparse AUTO edge/scalar", counts);
    leopard::ff16::TestOnlySparseEncodeCounts route =
        leopard::ff16::TestOnlyGetSparseEncodeCounts();
    require(route.exact_blocks == 4 &&
            route.retained_butterflies < route.prefix_butterflies &&
            route.requested_output_copies == 5,
        "edge mask did not select the exact sparse schedule");

    leopard::ff16::TestOnlyResetSparseEncodeCounts();
    const EncodeResult scattered_actual = encode(
        target->codec, LEO2_TEST_ENCODE_AUTO, aligned_original, scattered);
    require_result(scattered_actual.result, "sparse AUTO scattered encode");
    compare_requested(scattered_actual.recovery,
        scattered_reference.recovery, scattered, 0xa5,
        "sparse AUTO scattered/scalar", counts);
    route = leopard::ff16::TestOnlyGetSparseEncodeCounts();
    require(route.exact_blocks == 4 &&
            route.retained_butterflies < route.prefix_butterflies &&
            route.requested_output_copies == 6,
        "scattered mask did not select the exact sparse schedule");

    // Independently evaluate one native GF16 symbol at every requested edge
    // coordinate.  This checks the promoted path against direct algebra rather
    // than only against another transform implementation.
    const ProfileLayout layout = leopard2_test::make_profile_layout(
        leopard2_test::kLow, k, r);
    std::vector<Element> points(layout.parent_dimension, 0);
    std::vector<Element> values(layout.parent_dimension, 0);
    for (unsigned i = 0; i < layout.parent_dimension; ++i)
    {
        points[i] = static_cast<Element>(layout.systematic_coordinates[i]);
        values[i] = static_cast<Element>(aligned_original[i][0] |
            (static_cast<unsigned>(aligned_original[i][32]) << 8));
    }
    for (size_t i = 0;
         i < sizeof(edge_indices) / sizeof(edge_indices[0]); ++i)
    {
        const unsigned parity = edge_indices[i];
        const Element expected = leopard2_test::lagrange_evaluate(
            gf16, points, values,
            static_cast<Element>(layout.parity_coordinates[parity]));
        const Element actual = static_cast<Element>(
            edge_actual.recovery[parity][0] |
            (static_cast<unsigned>(edge_actual.recovery[parity][32]) << 8));
        require(actual == expected,
            "sparse AUTO edge output differs from direct GF16 algebra");
        ++counts->parity_symbols;
    }

    // 64-byte and near-mask controls must retain the mature prefix evaluator.
    const EncodeResult small_reference = encode(
        scalar->codec, LEO2_TEST_ENCODE_AUTO, small_original, edge);
    leopard::ff16::TestOnlyResetSparseEncodeCounts();
    const EncodeResult small_actual = encode(
        target->codec, LEO2_TEST_ENCODE_AUTO, small_original, edge);
    require_result(small_actual.result, "sparse AUTO 64-byte fallback encode");
    compare_requested(small_actual.recovery, small_reference.recovery, edge,
        0xa5, "sparse AUTO 64-byte fallback/scalar", counts);
    require(leopard::ff16::TestOnlyGetSparseEncodeCounts().exact_blocks == 0,
        "64-byte negative cell selected the exact sparse schedule");

    std::vector<uint8_t> near_edge = edge;
    near_edge[1] = 1;
    const EncodeResult near_reference = encode(
        scalar->codec, LEO2_TEST_ENCODE_AUTO, aligned_original, near_edge);
    leopard::ff16::TestOnlyResetSparseEncodeCounts();
    const EncodeResult near_actual = encode(
        target->codec, LEO2_TEST_ENCODE_AUTO, aligned_original, near_edge);
    require_result(near_actual.result, "sparse AUTO near-mask fallback encode");
    compare_requested(near_actual.recovery, near_reference.recovery, near_edge,
        0xa5, "sparse AUTO near-mask fallback/scalar", counts);
    require(leopard::ff16::TestOnlyGetSparseEncodeCounts().exact_blocks == 0,
        "unmeasured near-mask selected the exact sparse schedule");

    std::vector<uint8_t> dense_prefix(r, 0);
    std::fill(dense_prefix.begin(), dense_prefix.begin() + 128, 1);
    const EncodeResult dense_reference = encode(
        scalar->codec, LEO2_TEST_ENCODE_AUTO, aligned_original, dense_prefix);
    leopard::ff16::TestOnlyResetSparseEncodeCounts();
    const EncodeResult dense_actual = encode(
        target->codec, LEO2_TEST_ENCODE_AUTO, aligned_original, dense_prefix);
    require_result(dense_actual.result, "sparse AUTO dense fallback encode");
    compare_requested(dense_actual.recovery, dense_reference.recovery,
        dense_prefix, 0xa5, "sparse AUTO dense fallback/scalar", counts);
    require(leopard::ff16::TestOnlyGetSparseEncodeCounts().exact_blocks == 0,
        "dense-prefix control selected the exact sparse schedule");

    // A compact GF16 tail reuses the schedule compiled before either pass.
    const EncodeResult tail_reference = encode(
        scalar->codec, LEO2_TEST_ENCODE_AUTO, original, edge);
    leopard::ff16::TestOnlyResetSparseEncodeCounts();
    const EncodeResult tail_actual = encode(
        target->codec, LEO2_TEST_ENCODE_AUTO, original, edge);
    require_result(tail_actual.result, "sparse AUTO tail encode");
    compare_requested(tail_actual.recovery, tail_reference.recovery, edge,
        0xa5, "sparse AUTO tail/scalar", counts);
    route = leopard::ff16::TestOnlyGetSparseEncodeCounts();
    require(route.exact_blocks == 8 && route.requested_output_copies == 10,
        "GF16 compact tail did not execute the compiled schedule twice");
    for (unsigned i = 0; i < k; ++i)
        values[i] = static_cast<Element>(original[i][1024] |
            (static_cast<unsigned>(original[i][1025]) << 8));
    for (size_t i = 0;
         i < sizeof(edge_indices) / sizeof(edge_indices[0]); ++i)
    {
        const unsigned parity = edge_indices[i];
        const Element expected = leopard2_test::lagrange_evaluate(
            gf16, points, values,
            static_cast<Element>(layout.parity_coordinates[parity]));
        const Element actual = static_cast<Element>(
            tail_actual.recovery[parity][1024] |
            (static_cast<unsigned>(tail_actual.recovery[parity][1025]) << 8));
        require(actual == expected,
            "sparse AUTO compact tail differs from direct GF16 algebra");
        ++counts->parity_symbols;
    }

    // Prepare the complete invocation before enabling the allocation trap.
    const std::vector<const void*> input = const_pointers(original);
    Shards audited_output(r, Bytes(tail_bytes, 0xa5));
    std::vector<void*> audited_recovery(r, NULL);
    for (unsigned i = 0; i < r; ++i)
        if (edge[i])
            audited_recovery[i] = &audited_output[i][0];
    size_t tail_scratch_bytes = 0;
    size_t scalar_tail_scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(
        target->codec, tail_bytes, &tail_scratch_bytes),
        "sparse AUTO tail scratch query");
    require_result(leo2_encode_scratch_size(
        scalar->codec, tail_bytes, &scalar_tail_scratch_bytes),
        "sparse AUTO scalar tail scratch query");
    require(tail_scratch_bytes == scalar_tail_scratch_bytes,
        "hook-build tail scratch upper bound changed across backends");
    AlignedBuffer tail_scratch(tail_scratch_bytes);
    begin_allocation_audit();
    const leo2_result audited_result = leo2_encode(target->codec, tail_bytes,
        &input[0], &audited_recovery[0], tail_scratch.data(),
        tail_scratch.size());
    const uint64_t hot_allocations = end_allocation_audit();
    require_result(audited_result, "allocation-trapped sparse AUTO encode");
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
    require(hot_allocations == 0,
        "sparse AUTO encode allocated C++ storage");
    ++counts->allocation_checks;
#else
    (void)hot_allocations;
#endif
    compare_requested(audited_output, tail_reference.recovery, edge, 0xa5,
        "allocation-trapped sparse AUTO/scalar", counts);

    // Validation and scratch failures occur before schedule compilation or any
    // output write.  The UINT64_MAX query also locks overflow behavior.
    Shards rejected_output(r, Bytes(tail_bytes, 0x6d));
    const Shards rejected_before = rejected_output;
    std::vector<void*> rejected_recovery(r, NULL);
    for (unsigned i = 0; i < r; ++i)
        if (edge[i])
            rejected_recovery[i] = &rejected_output[i][0];
    leopard::ff16::TestOnlyResetSparseEncodeCounts();
    require(leo2_encode(target->codec, tail_bytes, &input[0],
                &rejected_recovery[0], tail_scratch.data(),
                tail_scratch.size() - 1) == LEO2_SCRATCH_TOO_SMALL,
        "sparse AUTO accepted one-byte-short scratch");
    require(rejected_output == rejected_before &&
            leopard::ff16::TestOnlyGetSparseEncodeCounts().exact_blocks == 0,
        "short scratch compiled or partially wrote sparse AUTO output");
    rejected_recovery[edge_indices[0]] = const_cast<void*>(input[0]);
    const uint64_t input_hash = hash_shards(original);
    require(leo2_encode(target->codec, tail_bytes, &input[0],
                &rejected_recovery[0], tail_scratch.data(),
                tail_scratch.size()) == LEO2_OVERLAP,
        "sparse AUTO accepted output/input aliasing");
    require(hash_shards(original) == input_hash &&
            rejected_output == rejected_before &&
            leopard::ff16::TestOnlyGetSparseEncodeCounts().exact_blocks == 0,
        "alias rejection compiled or modified sparse AUTO buffers");
    size_t invalid_scratch = 123;
    require(leo2_encode_scratch_size(target->codec, UINT64_MAX,
                &invalid_scratch) == LEO2_INVALID_ARGUMENT &&
            invalid_scratch == 0,
        "sparse AUTO scratch overflow was not rejected atomically");
    invalid_scratch = 123;
    require(leo2_encode_scratch_size(target->codec, 1025,
                &invalid_scratch) == LEO2_UNSUPPORTED &&
            invalid_scratch == 0,
        "sparse AUTO accepted an odd native GF16 shard length");
    counts->contract_checks += 4;

    // Recover one missing original from the sparse parity set, then re-encode
    // the repaired message and require exact parity reproduction.
    const unsigned missing = 37;
    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present = edge;
    original_present[missing] = 0;
    leo2_decode_plan* plan = NULL;
    require_result(leo2_decode_plan_create(target->codec,
        &original_present[0], &recovery_present[0], &plan),
        "sparse AUTO decode plan create");
    std::vector<const void*> decode_original = input;
    decode_original[missing] = NULL;
    std::vector<const void*> decode_recovery(r, NULL);
    for (unsigned i = 0; i < r; ++i)
        if (edge[i])
            decode_recovery[i] = &tail_actual.recovery[i][0];
    Shards restored(k, Bytes(tail_bytes, 0));
    std::vector<void*> restored_output(k, NULL);
    restored_output[missing] = &restored[missing][0];
    size_t decode_scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(
        plan, tail_bytes, &decode_scratch_bytes),
        "sparse AUTO decode scratch query");
    AlignedBuffer decode_scratch(decode_scratch_bytes);
    require_result(leo2_decode_plan_execute(plan, tail_bytes,
        &decode_original[0], &decode_recovery[0], &restored_output[0],
        decode_scratch.data(), decode_scratch.size()),
        "sparse AUTO decode execute");
    require(restored[missing] == original[missing],
        "sparse AUTO parity did not recover the missing original");
    Shards repaired = original;
    repaired[missing] = restored[missing];
    const EncodeResult rebuilt = encode(
        target->codec, LEO2_TEST_ENCODE_AUTO, repaired, edge);
    require_result(rebuilt.result, "sparse AUTO parity rebuild");
    compare_requested(rebuilt.recovery, tail_actual.recovery, edge, 0xa5,
        "sparse AUTO rebuilt parity", counts);
    leo2_decode_plan_destroy(plan);

    // Batch and concurrent callers share only immutable codec state; each
    // execution owns its schedule scratch and outputs.
    struct Invocation
    {
        Invocation(unsigned output_count, size_t shard_bytes,
            size_t scratch_bytes)
            : output(output_count, Bytes(shard_bytes, 0xa5))
            , pointers(output_count, NULL), scratch(scratch_bytes)
        {}
        Shards output;
        std::vector<void*> pointers;
        AlignedBuffer scratch;
    };
    Invocation edge_batch(r, bytes, target_scratch_bytes);
    Invocation scattered_batch(r, bytes, target_scratch_bytes);
    for (unsigned i = 0; i < r; ++i)
    {
        if (edge[i])
            edge_batch.pointers[i] = &edge_batch.output[i][0];
        if (scattered[i])
            scattered_batch.pointers[i] = &scattered_batch.output[i][0];
    }
    const std::vector<const void*> aligned_input =
        const_pointers(aligned_original);
    leo2_encode_batch_item batch_items[2] = {
        { bytes, &aligned_input[0], &edge_batch.pointers[0],
            edge_batch.scratch.data(), edge_batch.scratch.size() },
        { bytes, &aligned_input[0], &scattered_batch.pointers[0],
            scattered_batch.scratch.data(), scattered_batch.scratch.size() }
    };
    require_result(leo2_encode_batch(target->codec, batch_items, 2),
        "sparse AUTO batch encode");
    compare_requested(edge_batch.output, edge_reference.recovery, edge, 0xa5,
        "sparse AUTO edge batch/scalar", counts);
    compare_requested(scattered_batch.output, scattered_reference.recovery,
        scattered, 0xa5, "sparse AUTO scattered batch/scalar", counts);
    counts->batch_executions += 2;

    std::atomic<unsigned> failures(0);
    std::vector<std::thread> threads;
    for (unsigned worker = 0; worker < 2; ++worker)
    {
        threads.push_back(std::thread([&]() {
            try
            {
                Invocation call(r, bytes, target_scratch_bytes);
                for (unsigned i = 0; i < r; ++i)
                    if (edge[i])
                        call.pointers[i] = &call.output[i][0];
                if (leo2_encode(target->codec, bytes, &aligned_input[0],
                        &call.pointers[0], call.scratch.data(),
                        call.scratch.size()) != LEO2_SUCCESS)
                {
                    failures.fetch_add(1, std::memory_order_relaxed);
                    return;
                }
                for (unsigned i = 0; i < r; ++i)
                    if (call.output[i] != edge_reference.recovery[i])
                    {
                        failures.fetch_add(1, std::memory_order_relaxed);
                        return;
                    }
            }
            catch (...)
            {
                failures.fetch_add(1, std::memory_order_relaxed);
            }
        }));
    }
    for (size_t i = 0; i < threads.size(); ++i)
        threads[i].join();
    require(failures.load(std::memory_order_relaxed) == 0,
        "concurrent sparse AUTO encoding was nondeterministic");
    counts->concurrent_executions += 2;
    counts->sparse_auto_promotion_checks += 18;

    delete scalar;
    delete target;
    leo2_context_destroy(scalar_context);
    leo2_context_destroy(target_context);
}

void test_concurrent_codec(leo2_context* context, Counts* counts)
{
    const unsigned k = 9;
    const unsigned r = 7;
    const size_t bytes = 65;
    CodecOwner* owner = make_codec(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    require_result(leo2_test_codec_set_encode_mode(owner->codec,
        LEO2_TEST_ENCODE_FORCE_DIRECT), "force direct concurrent");
    const Shards original = random_shards(k, bytes, UINT64_C(0xc011ab1e));
    const std::vector<const void*> input = const_pointers(original);
    const std::vector<uint8_t> all(r, 1);
    const EncodeResult reference = encode(owner->codec,
        LEO2_TEST_ENCODE_FORCE_DIRECT, original, all);
    require_result(reference.result, "concurrent reference encode");
    size_t scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(owner->codec, bytes, &scratch_bytes),
        "concurrent scratch query");

    struct Invocation
    {
        Invocation(unsigned r_, size_t bytes_, size_t scratch_bytes_)
            : output(r_, Bytes(bytes_, 0)), pointers(r_, NULL), scratch(scratch_bytes_)
        {
            for (unsigned i = 0; i < r_; ++i)
                pointers[i] = &output[i][0];
        }
        Shards output;
        std::vector<void*> pointers;
        AlignedBuffer scratch;
    };

    const unsigned workers = 4;
    const unsigned repeats = 8;
    std::vector<Invocation*> invocation(workers, NULL);
    for (unsigned i = 0; i < workers; ++i)
        invocation[i] = new Invocation(r, bytes, scratch_bytes);
    std::atomic<unsigned> failures(0);
    std::vector<std::thread> threads;
    for (unsigned worker = 0; worker < workers; ++worker)
        threads.push_back(std::thread([&, worker]() {
            for (unsigned repeat = 0; repeat < repeats; ++repeat)
            {
                Invocation* call = invocation[worker];
                if (leo2_encode(owner->codec, bytes, &input[0], &call->pointers[0],
                        call->scratch.data(), call->scratch.size()) != LEO2_SUCCESS ||
                    call->output != reference.recovery)
                    failures.fetch_add(1, std::memory_order_relaxed);
            }
        }));
    for (size_t i = 0; i < threads.size(); ++i)
        threads[i].join();
    require(failures.load(std::memory_order_relaxed) == 0,
        "concurrent direct encoding was nondeterministic");
    for (unsigned i = 0; i < workers; ++i)
        delete invocation[i];
    counts->concurrent_executions += workers * repeats;
    delete owner;
}

} // namespace

int main()
{
    try
    {
        leo2_context_options options = {};
        options.struct_size = sizeof(options);
        options.thread_count = 1;
        leo2_context* context = NULL;
        require_result(leo2_context_create(&options, &context), "context create");

        const BinaryField gf8 = leopard2_test::make_legacy_gf8();
        const BinaryField gf16 = leopard2_test::make_legacy_gf16();
        Counts counts;
        test_profile_matrix(context, gf8,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, &counts);
        test_profile_matrix(context, gf8,
            LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8, &counts);
        test_profile_matrix(context, gf16,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, &counts);
        test_profile_matrix(context, gf16,
            LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16, &counts);
        test_capability_boundaries(context, &counts);
        test_sparse_schedule_budget_fallback(context, &counts);
        test_low_transform_no_coefficient_copy(
            context, gf8, gf16, &counts);
        test_low_p16_partial_direct_output(gf8, &counts);
        test_high_transform_source_staging(
            context, gf8, gf16, &counts);
        test_high_small_coarse_kernel(context, gf8, &counts);
        test_high_t2_t4_tiny_multi_block(gf8, &counts);
        test_high_gf16_byte_tiling(context, gf16, &counts);
        test_gf8_high_tail_generator_column(gf8, &counts);
#if LEO2_EXPERIMENT_HIGH_HALF_TAIL_COLUMN
        test_gf8_high_half_tail_generator_column(gf8, &counts);
#endif
        test_gf8_high_two_tail_generator_columns(gf8, &counts);
        test_high_gf8_byte_tiling(gf8, &counts);
        test_gf8_high_coarse_direct_oracle(gf8, &counts);
        test_auto_dispatch_threshold(context, &counts);
        test_k1_binding_respects_encode_mode(&counts);
        test_tail_allocation_and_contracts(context, gf8, gf16, &counts);
        test_unaligned_guarded_buffers(context, gf8, gf16, &counts);
        test_auto_encode_batch(context, &counts);
        test_sparse_low_gf16_auto_promotion(gf16, &counts);
        test_concurrent_codec(context, &counts);
        leo2_context_destroy(context);

        std::cout << "leopard2 direct encode passed: profiles=" << counts.profiles
                  << " basis_messages=" << counts.basis_messages
                  << " random_messages=" << counts.random_messages
                  << " parity_symbols=" << counts.parity_symbols
                  << " mask_executions=" << counts.mask_executions
                  << " boundary_profiles=" << counts.boundary_profiles
                  << " allocation_checks=" << counts.allocation_checks
                  << " concurrent_executions=" << counts.concurrent_executions
                  << " contract_checks=" << counts.contract_checks
                  << " dispatch_checks=" << counts.dispatch_checks
                  << " unaligned_checks=" << counts.unaligned_checks
                  << " batch_executions=" << counts.batch_executions
                  << " no_copy_checks=" << counts.no_copy_checks
                  << " low_partial_output_checks="
                  << counts.low_partial_output_checks
                  << " high_source_staging_checks="
                  << counts.high_source_staging_checks
                  << " high_byte_tiling_checks="
                  << counts.high_byte_tiling_checks
                  << " high_tail_column_checks="
                  << counts.high_tail_column_checks
                  << " gf8_coarse_oracle_checks="
                  << counts.gf8_coarse_oracle_checks
                  << " high_small_transform_checks="
                  << counts.high_small_transform_checks
                  << " sparse_auto_promotion_checks="
                  << counts.sparse_auto_promotion_checks
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
                  << " allocation_audit=enabled"
#else
                  << " allocation_audit=disabled-thread-sanitizer"
#endif
                  << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "leopard2 direct encode failed: " << error.what() << std::endl;
        return 1;
    }
}
