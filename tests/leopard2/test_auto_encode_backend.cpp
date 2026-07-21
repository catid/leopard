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

#include "Leopard2Backend.h"
#include "Leopard2Direct.h"
#include "LeopardFF8.h"
#include "leopard.h"
#include "leopard2.h"

#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

#ifndef LEO2_ENABLE_TEST_HOOKS
#error "The AUTO encode-backend test requires LEO2_ENABLE_TEST_HOOKS"
#endif

namespace {

typedef std::vector<uint8_t> Bytes;
typedef std::vector<Bytes> Shards;

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void require_result(
    leo2_result actual, leo2_result expected, const char* message)
{
    if (actual != expected)
        throw std::runtime_error(std::string(message) + ": " +
            leo2_result_string(actual));
}

class Context
{
public:
    explicit Context(leo2_backend requested)
        : value_(NULL), result_(LEO2_INTERNAL_ERROR)
    {
        leo2_context_options options;
        std::memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = requested;
        options.thread_count = 1;
        result_ = leo2_context_create(&options, &value_);
    }

    ~Context() { leo2_context_destroy(value_); }
    leo2_context* get() const { return value_; }
    leo2_result result() const { return result_; }

private:
    Context(const Context&);
    Context& operator=(const Context&);
    leo2_context* value_;
    leo2_result result_;
};

class Codec
{
public:
    Codec(
        leo2_context* context,
        uint32_t k,
        uint32_t r,
        leo2_profile profile = LEO2_PROFILE_LEGACY_HIGH_V1,
        leo2_field field = LEO2_FIELD_GF16)
        : value_(NULL)
    {
        require_result(leo2_codec_create(context, k, r,
            profile, field, NULL, &value_),
            LEO2_SUCCESS, "codec create");
    }

    ~Codec() { leo2_codec_destroy(value_); }
    leo2_codec* get() const { return value_; }

private:
    Codec(const Codec&);
    Codec& operator=(const Codec&);
    leo2_codec* value_;
};

class AlignedBuffer
{
public:
    explicit AlignedBuffer(size_t bytes)
        : value_(NULL)
    {
#if defined(_MSC_VER)
        value_ = _aligned_malloc(bytes, leo2_scratch_alignment());
#else
        if (posix_memalign(&value_, leo2_scratch_alignment(), bytes) != 0)
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

    void* get() const { return value_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);
    void* value_;
};

leo2_backend selected_backend(
    const leo2_codec* codec,
    uint64_t bytes,
    uint32_t count,
    uint32_t prefix)
{
    leo2_backend backend = LEO2_BACKEND_AUTO;
    require_result(leo2_test_codec_transform_encode_backend(
        codec, bytes, count, prefix, &backend), LEO2_SUCCESS,
        "encode backend query");
    return backend;
}

Shards make_original(uint32_t k, size_t bytes)
{
    Shards original(k, Bytes(bytes));
    uint32_t state = 0x243f6a88U;
    for (uint32_t shard = 0; shard < k; ++shard)
        for (size_t i = 0; i < bytes; ++i)
        {
            state = state * 1664525U + 1013904223U;
            original[shard][i] = static_cast<uint8_t>(
                (state >> 24) ^ shard ^ static_cast<uint32_t>(i));
        }
    return original;
}

Shards encode(
    const leo2_codec* codec,
    const Shards& original,
    uint32_t r,
    bool sparse)
{
    const size_t bytes = original[0].size();
    std::vector<const void*> original_ptrs(original.size());
    for (size_t i = 0; i < original.size(); ++i)
        original_ptrs[i] = original[i].data();
    Shards recovery(r, Bytes(bytes));
    std::vector<void*> recovery_ptrs(r);
    for (uint32_t i = 0; i < r; ++i)
        recovery_ptrs[i] = sparse && (i & 1U) != 0
            ? NULL : recovery[i].data();
    size_t scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        LEO2_SUCCESS, "scratch query");
    AlignedBuffer scratch(scratch_bytes);
    require_result(leo2_encode(codec, bytes, original_ptrs.data(),
        recovery_ptrs.data(), scratch.get(), scratch_bytes),
        LEO2_SUCCESS, "encode");
    return recovery;
}

Shards encode_unaligned(
    const leo2_codec* codec,
    const Shards& original,
    uint32_t r)
{
    const size_t bytes = original[0].size();
    Shards original_storage(original.size(), Bytes(bytes + 3));
    std::vector<const void*> original_ptrs(original.size());
    for (size_t i = 0; i < original.size(); ++i)
    {
        std::memcpy(original_storage[i].data() + 1,
            original[i].data(), bytes);
        original_ptrs[i] = original_storage[i].data() + 1;
    }

    Shards recovery_storage(r, Bytes(bytes + 5));
    std::vector<void*> recovery_ptrs(r);
    for (uint32_t i = 0; i < r; ++i)
        recovery_ptrs[i] = recovery_storage[i].data() + 3;

    size_t scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        LEO2_SUCCESS, "unaligned scratch query");
    AlignedBuffer scratch(scratch_bytes);
    require_result(leo2_encode(codec, bytes, original_ptrs.data(),
        recovery_ptrs.data(), scratch.get(), scratch_bytes),
        LEO2_SUCCESS, "unaligned encode");

    Shards recovery(r, Bytes(bytes));
    for (uint32_t i = 0; i < r; ++i)
        std::memcpy(recovery[i].data(),
            recovery_storage[i].data() + 3, bytes);
    return recovery;
}

Shards encode_legacy(const Shards& original, uint32_t r)
{
    const size_t bytes = original[0].size();
    const uint32_t k = static_cast<uint32_t>(original.size());
    std::vector<const void*> original_ptrs(k);
    for (uint32_t i = 0; i < k; ++i)
        original_ptrs[i] = original[i].data();
    const unsigned work_count = leo_encode_work_count(k, r);
    require(work_count >= r, "legacy encode work count");
    Shards work(work_count, Bytes(bytes));
    std::vector<void*> work_ptrs(work_count);
    for (unsigned i = 0; i < work_count; ++i)
        work_ptrs[i] = work[i].data();
    require(leo_encode(bytes, k, r, work_count,
                original_ptrs.data(), work_ptrs.data()) == Leopard_Success,
        "legacy encode");
    work.resize(r);
    return work;
}

void require_sparse_matches_full(
    const Shards& sparse, const Shards& full)
{
    require(sparse.size() == full.size(), "sparse parity size mismatch");
    for (size_t i = 0; i < full.size(); i += 2)
        require(sparse[i] == full[i], "sparse parity differs from full parity");
}

void require_explicit_backend(Context& context, leo2_backend expected)
{
    if (context.result() != LEO2_SUCCESS)
        return;
    require(leo2_context_backend(context.get()) == expected,
        "explicit context reported the wrong backend");
    Codec codec(context.get(), 1000, 200);
    require(selected_backend(codec.get(), 4U * 1024U * 1024U + 64U,
                200, 200) == expected,
        "explicit backend was changed by the AUTO calibration bounds");
}

void test_selection_and_bytes(
    Context& automatic,
    Context& scalar,
    Context& ssse3,
    Context& avx2,
    Context& avx512)
{
    Codec auto_codec(automatic.get(), 1000, 200);
    Codec avx2_codec(avx2.get(), 1000, 200);
    Codec avx512_codec(avx512.get(), 1000, 200);

    require(selected_backend(auto_codec.get(), 32, 200, 200) ==
            LEO2_BACKEND_AVX2,
        "AUTO widened a partial transform tile");
    require(selected_backend(auto_codec.get(), 63, 200, 200) ==
            LEO2_BACKEND_AVX2,
        "AUTO widened immediately below the minimum shard length");
    require(selected_backend(auto_codec.get(), 64, 200, 200) ==
            LEO2_BACKEND_AVX512,
        "AUTO did not widen at the complete-tile boundary");
    require(selected_backend(auto_codec.get(), 65, 200, 200) ==
            LEO2_BACKEND_AVX2,
        "AUTO widened an immediate ragged-tail boundary");
    require(selected_backend(auto_codec.get(), 4098, 200, 200) ==
            LEO2_BACKEND_AVX2,
        "AUTO widened an uncalibrated tail");
    require(selected_backend(auto_codec.get(), 4U * 1024U * 1024U,
                200, 200) ==
            LEO2_BACKEND_AVX512,
        "AUTO did not widen a calibrated large complete tile");
    require(selected_backend(auto_codec.get(),
                4U * 1024U * 1024U + 64U, 200, 200) ==
            LEO2_BACKEND_AVX2,
        "AUTO widened immediately above the calibrated shard-length range");
    require(selected_backend(auto_codec.get(), 4096, 199, 200) ==
            LEO2_BACKEND_AVX2,
        "AUTO widened a partial-output encode");
    require(selected_backend(avx2_codec.get(), 4096, 200, 200) ==
            LEO2_BACKEND_AVX2,
        "explicit AVX2 widened");
    require(selected_backend(avx512_codec.get(), 4032, 200, 200) ==
            LEO2_BACKEND_AVX512,
        "explicit AVX512 did not remain exact");
    require(selected_backend(avx512_codec.get(),
                4U * 1024U * 1024U + 64U, 200, 200) ==
            LEO2_BACKEND_AVX512,
        "explicit AVX512 was constrained by AUTO-only byte bounds");

    Codec large(automatic.get(), 4096, 512);
    require(selected_backend(large.get(), 4032, 512, 512) ==
            LEO2_BACKEND_AVX512,
        "large AUTO did not widen a complete tile");
    require(selected_backend(large.get(), 4096, 512, 512) ==
            LEO2_BACKEND_AVX512,
        "large AUTO did not widen in the calibrated cell");
    require(selected_backend(large.get(), 4160, 512, 512) ==
            LEO2_BACKEND_AVX512,
        "large AUTO did not widen a neighboring complete tile");

    Codec minimum_shape(automatic.get(), 8, 2);
    require(leo2_codec_parent_count(minimum_shape.get()) == 16,
        "minimum calibrated parent changed");
    require(selected_backend(minimum_shape.get(), 64, 2, 2) ==
            LEO2_BACKEND_AVX512,
        "minimum K/N transform shape did not widen");
    Codec below_k(automatic.get(), 7, 2);
    require(leo2_codec_parent_count(below_k.get()) == 16,
        "below-K parent does not isolate the K boundary");
    require(selected_backend(below_k.get(), 64, 2, 2) ==
            LEO2_BACKEND_AVX2,
        "AUTO widened immediately below K=8");
    Codec below_parent(automatic.get(), 7, 1);
    require(leo2_codec_parent_count(below_parent.get()) == 8,
        "immediate lower parent is not N=8");
    require(selected_backend(below_parent.get(), 64, 1, 1) ==
            LEO2_BACKEND_AVX2,
        "AUTO widened below N=16");
    Codec direct_xor(automatic.get(), 8, 1);
    require(leo2_codec_parent_count(direct_xor.get()) == 16,
        "R=1 boundary did not retain N=16");
    require(selected_backend(direct_xor.get(), 4096, 1, 1) ==
            LEO2_BACKEND_AVX2,
        "T=1 direct codec unnecessarily qualified transform widening");

    Codec maximum_r(automatic.get(), 8, 4096);
    require(selected_backend(maximum_r.get(), 64, 4096, 4096) ==
            LEO2_BACKEND_AVX512,
        "AUTO did not widen at R=4096");
    Codec above_r(automatic.get(), 8, 4097);
    require(selected_backend(above_r.get(), 64, 4097, 4097) ==
            LEO2_BACKEND_AVX2,
        "AUTO widened immediately above R=4096");

    Codec low_profile(automatic.get(), 8, 8, LEO2_PROFILE_LOW_V1);
    require(selected_backend(low_profile.get(), 64, 8, 8) ==
            LEO2_BACKEND_AVX2,
        "AUTO widened a non-legacy-high profile");
    Codec gf8_unbalanced(automatic.get(), 8, 2,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    require(selected_backend(gf8_unbalanced.get(), 4096, 2, 2) ==
            LEO2_BACKEND_AVX2,
        "AUTO widened an unqualified GF8 codec");

    static const uint32_t gf8_sides[] = { 16, 32, 64 };
    for (size_t side_i = 0;
         side_i < sizeof(gf8_sides) / sizeof(gf8_sides[0]); ++side_i)
    {
        const uint32_t side = gf8_sides[side_i];
        Codec balanced(automatic.get(), side, side,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        const uint64_t minimum_bytes = side == 16 ? 4096 : 2048;
        require(selected_backend(
                    balanced.get(), minimum_bytes - 64, side, side) ==
                LEO2_BACKEND_AVX2,
            "balanced GF8 AUTO widened below the calibrated byte range");
        require(selected_backend(
                    balanced.get(), minimum_bytes, side, side) ==
                LEO2_BACKEND_AVX512,
            "balanced GF8 AUTO did not widen at the lower byte boundary");
        if (side == 16)
            require(selected_backend(balanced.get(), 2048, side, side) ==
                    LEO2_BACKEND_AVX2,
                "balanced GF8 T=16 widened at the inconclusive 2 KiB cell");
        require(selected_backend(balanced.get(), 4096, side, side) ==
                LEO2_BACKEND_AVX512,
            "balanced GF8 AUTO did not widen in the calibrated region");
        require(selected_backend(balanced.get(), 65536, side, side) ==
                LEO2_BACKEND_AVX512,
            "balanced GF8 AUTO did not widen at the upper byte boundary");
        require(selected_backend(balanced.get(), 65600, side, side) ==
                LEO2_BACKEND_AVX2,
            "balanced GF8 AUTO widened above the calibrated byte range");
        require(selected_backend(balanced.get(), 4097, side, side) ==
                LEO2_BACKEND_AVX2,
            "balanced GF8 AUTO widened a ragged tail");
        require(selected_backend(balanced.get(), 4096, side - 1, side) ==
                LEO2_BACKEND_AVX2,
            "balanced GF8 AUTO widened a partial-output encode");

        Codec explicit_avx2(avx2.get(), side, side,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        Codec explicit_avx512(avx512.get(), side, side,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        Codec explicit_scalar(scalar.get(), side, side,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        Codec explicit_ssse3(ssse3.get(), side, side,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        require(selected_backend(explicit_avx2.get(), 4096, side, side) ==
                LEO2_BACKEND_AVX2,
            "balanced GF8 explicit AVX2 widened");
        require(selected_backend(explicit_avx512.get(), 4097, side, side) ==
                LEO2_BACKEND_AVX512,
            "balanced GF8 explicit AVX512 was constrained by AUTO bounds");

        const Shards balanced_original = make_original(side, 4096);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        const Shards balanced_auto = encode(
            balanced.get(), balanced_original, side, false);
        const leopard::ff8::TestOnlyHighEncodeCounts whole_counts =
            leopard::ff8::TestOnlyGetHighEncodeCounts();
        require(whole_counts.whole_transform_calls == 1,
            "balanced GF8 AUTO did not execute the coarse transform");
        const Shards balanced_avx2 = encode(
            explicit_avx2.get(), balanced_original, side, false);
        const Shards balanced_avx512 = encode(
            explicit_avx512.get(), balanced_original, side, false);
        require(balanced_auto == balanced_avx2 &&
                balanced_auto == balanced_avx512 &&
                balanced_auto == encode(
                    explicit_scalar.get(), balanced_original, side, false) &&
                balanced_auto == encode(
                    explicit_ssse3.get(), balanced_original, side, false) &&
                balanced_auto == encode_legacy(balanced_original, side),
            "balanced GF8 coarse transform changed legacy parity bytes");
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        require_sparse_matches_full(
            encode(balanced.get(), balanced_original, side, true),
            balanced_auto);
        require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                    whole_transform_calls == 0,
            "balanced GF8 sparse encode used the dense coarse transform");

        static const size_t tail_bytes[] = { 1025, 2049, 4097 };
        for (size_t tail_i = 0;
             tail_i < sizeof(tail_bytes) / sizeof(tail_bytes[0]); ++tail_i)
        {
            const Shards tail_original =
                make_original(side, tail_bytes[tail_i]);
            const Shards tail_reference = encode(
                explicit_scalar.get(), tail_original, side, false);
            const Shards tail_ssse3 = encode(
                explicit_ssse3.get(), tail_original, side, false);
            const Shards tail_avx2 = encode(
                explicit_avx2.get(), tail_original, side, false);
            leopard::ff8::TestOnlyResetHighEncodeCounts();
            const Shards tail_avx512 = encode(
                explicit_avx512.get(), tail_original, side, false);
            const uint64_t expected_whole_calls = tail_bytes[tail_i] >= 2049
                ? 1 : 0;
            require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                    whole_transform_calls == expected_whole_calls,
                "balanced GF8 tail used the wrong coarse-transform split");
            require(tail_reference == tail_ssse3 &&
                    tail_reference == tail_avx2 &&
                    tail_reference == tail_avx512,
                "balanced GF8 coarse transform changed a byte tail");
            require(tail_reference == encode_unaligned(
                        explicit_avx2.get(), tail_original, side) &&
                    tail_reference == encode_unaligned(
                        explicit_avx512.get(), tail_original, side),
                "balanced GF8 coarse transform mishandled unaligned buffers");
        }
    }

    Codec gf8_below_side(automatic.get(), 8, 8,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    Codec gf8_above_side(automatic.get(), 128, 128,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    Codec gf8_shortened(automatic.get(), 15, 16,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    Codec gf8_punctured(automatic.get(), 16, 15,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    require(selected_backend(gf8_below_side.get(), 4096, 8, 8) ==
            LEO2_BACKEND_AVX2,
        "balanced GF8 AUTO widened T=8");
    require(selected_backend(gf8_above_side.get(), 4096, 128, 128) ==
            LEO2_BACKEND_AVX2,
        "balanced GF8 AUTO widened T=128");
    require(selected_backend(gf8_shortened.get(), 4096, 16, 16) ==
            LEO2_BACKEND_AVX2,
        "balanced GF8 AUTO widened a shortened input block");
    require(selected_backend(gf8_punctured.get(), 4096, 15, 15) ==
            LEO2_BACKEND_AVX2,
        "balanced GF8 AUTO widened a punctured parity block");

    const Shards shortened_original = make_original(15, 4096);
    Codec shortened_scalar(scalar.get(), 15, 16,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    Codec shortened_avx2(avx2.get(), 15, 16,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    Codec shortened_avx512(avx512.get(), 15, 16,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    const Shards shortened_reference = encode(
        shortened_scalar.get(), shortened_original, 16, false);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    const Shards shortened_auto = encode(
        gf8_shortened.get(), shortened_original, 16, false);
    require(shortened_reference == shortened_auto &&
            shortened_reference == encode(
                shortened_avx2.get(), shortened_original, 16, false) &&
            shortened_reference == encode(
                shortened_avx512.get(), shortened_original, 16, false),
        "shortened GF8 neighbor changed legacy parity bytes");
    require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                whole_transform_calls == 0,
        "shortened GF8 neighbor used the exact balanced transform");

    const Shards punctured_original = make_original(16, 4096);
    const Shards punctured_reference = encode_legacy(punctured_original, 15);
    Codec punctured_scalar(scalar.get(), 16, 15,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    Codec punctured_avx2(avx2.get(), 16, 15,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    Codec punctured_avx512(avx512.get(), 16, 15,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    const Shards punctured_auto = encode(
        gf8_punctured.get(), punctured_original, 15, false);
    require(punctured_reference == punctured_auto &&
            punctured_reference == encode(
                punctured_scalar.get(), punctured_original, 15, false) &&
            punctured_reference == encode(
                punctured_avx2.get(), punctured_original, 15, false) &&
            punctured_reference == encode(
                punctured_avx512.get(), punctured_original, 15, false),
        "punctured GF8 neighbor changed legacy parity bytes");
    require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                whole_transform_calls == 0,
        "punctured GF8 neighbor used the exact balanced transform");

    Codec maximum_parent(automatic.get(), 60000, 1024);
    require(selected_backend(maximum_parent.get(), 1024, 1024, 1024) ==
            LEO2_BACKEND_AVX512,
        "maximum-parent GF16 high transform did not widen");

    const Shards original = make_original(1000, 4096);
    const Shards automatic_parity = encode(
        auto_codec.get(), original, 200, false);
    const Shards avx2_parity = encode(
        avx2_codec.get(), original, 200, false);
    const Shards avx512_parity = encode(
        avx512_codec.get(), original, 200, false);
    require(automatic_parity == avx2_parity &&
            automatic_parity == avx512_parity,
        "AUTO widening changed GF16 legacy parity bytes");
    require_sparse_matches_full(
        encode(auto_codec.get(), original, 200, true), automatic_parity);

    std::atomic<unsigned> failures(0);
    std::vector<std::thread> threads;
    for (unsigned lane = 0; lane < 4; ++lane)
    {
        threads.push_back(std::thread([&]() {
            try
            {
                require(encode(auto_codec.get(), original, 200, false) ==
                        automatic_parity,
                    "concurrent AUTO parity mismatch");
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
        "concurrent AUTO widening failed");
}

} // namespace

int main()
{
    try
    {
        require(leo_init() == Leopard_Success, "Leopard initialization");
        Context automatic(LEO2_BACKEND_AUTO);
        require_result(automatic.result(), LEO2_SUCCESS, "AUTO context");
        Context scalar(LEO2_BACKEND_SCALAR);
        Context ssse3(LEO2_BACKEND_SSSE3);
        Context avx2(LEO2_BACKEND_AVX2);
        Context avx512(LEO2_BACKEND_AVX512);

        require_explicit_backend(scalar, LEO2_BACKEND_SCALAR);
        require_explicit_backend(ssse3, LEO2_BACKEND_SSSE3);
        require_explicit_backend(avx2, LEO2_BACKEND_AVX2);
        require_explicit_backend(avx512, LEO2_BACKEND_AVX512);

        if (leopard::backend::IsCalibratedAutoAVX512EncodeHost() &&
            leo2_context_backend(automatic.get()) == LEO2_BACKEND_AVX2 &&
            scalar.result() == LEO2_SUCCESS &&
            ssse3.result() == LEO2_SUCCESS &&
            avx2.result() == LEO2_SUCCESS &&
            avx512.result() == LEO2_SUCCESS)
        {
            test_selection_and_bytes(
                automatic, scalar, ssse3, avx2, avx512);
        }
        else
        {
            Codec codec(automatic.get(), 1000, 200);
            require(selected_backend(codec.get(), 4096, 200, 200) ==
                    leo2_context_backend(automatic.get()),
                "unsupported host widened away from its AUTO baseline");
        }

        std::printf("Leopard2 AUTO encode backend passed\n");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "Leopard2 AUTO encode backend failed: %s\n",
            error.what());
        return 1;
    }
}
