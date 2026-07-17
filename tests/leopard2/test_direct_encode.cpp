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
#include "direct_oracle.h"
#include "leopard2.h"

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

    Counts()
        : profiles(0), basis_messages(0), random_messages(0), parity_symbols(0)
        , mask_executions(0), boundary_profiles(0), allocation_checks(0)
        , concurrent_executions(0), contract_checks(0), dispatch_checks(0)
        , unaligned_checks(0), batch_executions(0)
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
                    "basis direct and transform encoders differ");
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
            require(direct_path == 0,
                "AUTO selected an uncalibrated tiny direct-encode threshold");
            ++counts->random_messages;
            ++counts->profiles;
            delete owner;
        }
}

void test_capability_boundaries(leo2_context* context, Counts* counts)
{
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

void test_auto_dispatch_threshold(leo2_context* context, Counts* counts)
{
    const leo2_backend backend = leo2_context_backend(context);
    const bool qualified_simd = backend == LEO2_BACKEND_SSSE3 ||
        backend == LEO2_BACKEND_AVX2;
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
    require(direct == 0, "AUTO selected the unpromoted high profile");
    ++counts->dispatch_checks;
    delete high;

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

        g_tracked_allocations.store(0, std::memory_order_relaxed);
        g_track_allocations.store(true, std::memory_order_release);
        const leo2_result result = leo2_encode(owner->codec, c.bytes, &input[0],
            &recovery[0], scratch.data(), scratch.size());
        g_track_allocations.store(false, std::memory_order_release);
        require_result(result, "allocation-trapped tail encode");
        require(g_tracked_allocations.load(std::memory_order_relaxed) == 0,
            "direct encode allocated C++ storage");
        require(output == expected, "tail direct encode differs from oracle");
        ++counts->allocation_checks;

        require_result(leo2_test_codec_set_encode_mode(owner->codec,
            LEO2_TEST_ENCODE_FORCE_TRANSFORM), "force transform tail");
        output.assign(c.r, Bytes(c.bytes, 0xa5));
        for (unsigned i = 0; i < c.r; ++i)
            recovery[i] = &output[i][0];
        g_tracked_allocations.store(0, std::memory_order_relaxed);
        g_track_allocations.store(true, std::memory_order_release);
        const leo2_result transform_result = leo2_encode(
            owner->codec, c.bytes, &input[0], &recovery[0],
            scratch.data(), scratch.size());
        g_track_allocations.store(false, std::memory_order_release);
        require_result(transform_result, "allocation-trapped transform tail encode");
        require(g_tracked_allocations.load(std::memory_order_relaxed) == 0,
            "transform tail encode allocated C++ storage");
        require(output == expected,
            "tail transform encode differs from direct generator oracle");
        ++counts->allocation_checks;

        std::vector<uint8_t> subset(c.r, 0);
        subset[0] = subset[c.r / 2] = subset[c.r - 1] = 1;
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
    };
    const Case cases[] = {
        { LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 5, 3, 65 },
        { LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16, 5, 3, 66 }
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
        require_result(leo2_test_codec_set_encode_mode(owner->codec,
            LEO2_TEST_ENCODE_FORCE_DIRECT), "force direct unaligned");

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
            scratch.data(), scratch.size()), "unaligned direct encode");
        for (unsigned i = 0; i < c.r; ++i)
            require(memcmp(&output_storage[i][output_prefix],
                        &expected[i][0], c.bytes) == 0,
                "unaligned direct encode differs from the oracle");

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
        test_auto_dispatch_threshold(context, &counts);
        test_tail_allocation_and_contracts(context, gf8, gf16, &counts);
        test_unaligned_guarded_buffers(context, gf8, gf16, &counts);
        test_auto_encode_batch(context, &counts);
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
                  << " batch_executions=" << counts.batch_executions << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "leopard2 direct encode failed: " << error.what() << std::endl;
        return 1;
    }
}
