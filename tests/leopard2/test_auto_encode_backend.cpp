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

#include <algorithm>
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

void test_balanced_execution_tile_geometry()
{
    static const size_t requested_tile_bytes = 32U * 1024U;
    static const size_t requested_tile_count =
        requested_tile_bytes / leo2_scratch_alignment();
    static const size_t pass_counts[] = { 1, 2, 63, 64, 65, 129, 513 };
    for (size_t count_i = 0;
         count_i < sizeof(pass_counts) / sizeof(pass_counts[0]); ++count_i)
    {
        const size_t expected_count = pass_counts[count_i];
        for (size_t final_pass_tiles = 1;
             final_pass_tiles <= requested_tile_count;
             ++final_pass_tiles)
        {
            const size_t total_tiles =
                (expected_count - 1U) * requested_tile_count +
                final_pass_tiles;
            const size_t aligned_bytes =
                total_tiles * leo2_scratch_alignment();
            size_t execution_tile_count = 0;
            size_t maximum_pass_bytes = 0;
            require_result(leo2_test_balanced_execution_tiles(
                    aligned_bytes, requested_tile_bytes,
                    &execution_tile_count, &maximum_pass_bytes),
                LEO2_SUCCESS, "balanced execution-tile geometry");
            require(execution_tile_count == expected_count,
                "balanced execution-tile count differs from ceiling division");

            size_t remaining_tiles = total_tiles;
            size_t reference_maximum_tiles = 0;
            for (size_t pass = 0; pass < expected_count; ++pass)
            {
                const size_t passes_left = expected_count - pass;
                const size_t pass_tiles =
                    remaining_tiles / passes_left +
                    (remaining_tiles % passes_left != 0);
                reference_maximum_tiles = std::max(
                    reference_maximum_tiles, pass_tiles);
                remaining_tiles -= pass_tiles;
            }
            require(remaining_tiles == 0,
                "balanced execution reference did not consume every tile");
            require(maximum_pass_bytes ==
                    reference_maximum_tiles * leo2_scratch_alignment(),
                "balanced execution scratch is smaller than a distributed pass");
        }
    }

    size_t count = 99;
    size_t bytes = 99;
    require_result(leo2_test_balanced_execution_tiles(
            0, 0, &count, &bytes),
        LEO2_SUCCESS, "empty balanced execution geometry");
    require(count == 0 && bytes == 0,
        "empty balanced execution geometry is not empty");
    require_result(leo2_test_balanced_execution_tiles(
            64, 0, &count, &bytes),
        LEO2_INVALID_ARGUMENT, "zero requested execution tile");
    require_result(leo2_test_balanced_execution_tiles(
            65, 64, &count, &bytes),
        LEO2_INVALID_ARGUMENT, "unaligned execution byte count");
    require_result(leo2_test_balanced_execution_tiles(
            128, 65, &count, &bytes),
        LEO2_INVALID_ARGUMENT, "unaligned requested execution tile");
    require_result(leo2_test_balanced_execution_tiles(
            128, 64, NULL, &bytes),
        LEO2_INVALID_ARGUMENT, "null execution tile-count output");
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

void require_encode_overlap_rejected(
    const leo2_codec* codec,
    Shards& original,
    uint32_t r)
{
    const size_t bytes = original[0].size();
    const Shards original_before = original;
    std::vector<const void*> original_ptrs(original.size());
    for (size_t i = 0; i < original.size(); ++i)
        original_ptrs[i] = original[i].data();
    Shards recovery(r, Bytes(bytes));
    std::vector<void*> recovery_ptrs(r);
    for (uint32_t i = 0; i < r; ++i)
        recovery_ptrs[i] = recovery[i].data();
    recovery_ptrs[0] = original[0].data();
    size_t scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        LEO2_SUCCESS, "overlap scratch query");
    AlignedBuffer scratch(scratch_bytes);
    require_result(leo2_encode(codec, bytes, original_ptrs.data(),
        recovery_ptrs.data(), scratch.get(), scratch_bytes),
        LEO2_OVERLAP, "encode input/output overlap");
    require(original == original_before,
        "rejected encode overlap modified a source shard");
}

Shards encode_prefix(
    const leo2_codec* codec,
    const Shards& original,
    uint32_t r,
    uint32_t prefix)
{
    const size_t bytes = original[0].size();
    std::vector<const void*> original_ptrs(original.size());
    for (size_t i = 0; i < original.size(); ++i)
        original_ptrs[i] = original[i].data();
    Shards recovery(r, Bytes(bytes));
    std::vector<void*> recovery_ptrs(r, NULL);
    for (uint32_t i = 0; i < prefix; ++i)
        recovery_ptrs[i] = recovery[i].data();
    size_t scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        LEO2_SUCCESS, "prefix scratch query");
    AlignedBuffer scratch(scratch_bytes);
    require_result(leo2_encode(codec, bytes, original_ptrs.data(),
        recovery_ptrs.data(), scratch.get(), scratch_bytes),
        LEO2_SUCCESS, "prefix encode");
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

void test_small_high_encode(
    Context& scalar,
    Context& ssse3,
    Context& avx2,
    Context& avx512)
{
    if (avx2.result() != LEO2_SUCCESS)
        return;

    struct TestCase
    {
        uint32_t k;
        uint32_t r;
        size_t bytes;
    };
    static const TestCase cases[] = {
        // T=2 crosses at 2 KiB for every valid K.  Cover the minimum K, an
        // odd message count and ragged byte tail, and the largest GF8 K.
        { 2, 2, 2048 },
        { 3, 2, 2049 },
        { 254, 2, 2048 },
        // T=4 uses the evidence-derived K/byte staircase.  The new lower cells
        // fuse only the final forward butterfly; the existing thresholds keep
        // their register-fused or mature coarse implementations.
        { 3, 4, 2048 },
        { 3, 3, 2049 },
        { 3, 4, 4095 },
        { 3, 4, 4096 },
        { 3, 3, 4097 },
        { 4, 4, 2048 },
        { 4, 3, 4097 },
        { 4, 4, 8191 },
        { 4, 4, 8192 },
        { 4, 3, 8193 },
        { 5, 4, 2048 },
        { 6, 3, 2049 },
        { 7, 4, 2048 },
        { 7, 4, 4095 },
        { 7, 4, 4096 },
        { 16, 4, 2048 },
        { 16, 3, 2049 },
        { 8, 4, 2048 },
        { 8, 4, 4095 },
        { 8, 4, 4096 },
        { 8, 3, 4097 },
        { 9, 4, 2048 },
        { 10, 3, 2049 },
        { 11, 4, 2048 },
        { 12, 3, 2049 },
        // Retain the previously qualified interior and upper-bound cases.
        { 64, 2, 4096 },
        { 64, 3, 4097 },
        { 240, 4, 65536 },
        { 251, 3, 65536 },
        { 252, 4, 65536 }
    };

    for (size_t case_i = 0;
         case_i < sizeof(cases) / sizeof(cases[0]); ++case_i)
    {
        const TestCase& test_case = cases[case_i];
        const Shards original = make_original(test_case.k, test_case.bytes);
        const Shards original_before = original;
        Codec avx2_codec(avx2.get(), test_case.k, test_case.r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);

        leopard::ff8::TestOnlyResetHighEncodeCounts();
        const Shards actual = encode(
            avx2_codec.get(), original, test_case.r, false);
        const leopard::ff8::TestOnlyHighEncodeCounts route =
            leopard::ff8::TestOnlyGetHighEncodeCounts();
        // The aligned prefix is one execution pass.  A ragged 64-byte staging
        // pass deliberately remains below the production byte threshold.
        const uint64_t expected_passes = 1;
        const uint32_t side = test_case.r == 2 ? 2U : 4U;
        const bool fused_t4 = side == 4 &&
            (((test_case.k == 5 || test_case.k == 6 ||
                test_case.k >= 9) && test_case.bytes >= 2U * 1024U) ||
             ((test_case.k == 3 || test_case.k == 7) &&
                test_case.bytes >= 4U * 1024U) ||
             (test_case.k == 4 && test_case.bytes >= 8U * 1024U)) &&
            ((test_case.k >= 3 && test_case.k <= 7) ||
             (test_case.k >= 9 && test_case.k <= 11));
        const uint64_t expected_input_copies =
            (fused_t4 || test_case.k <= side ||
                    (side == 4 && test_case.k % side == 3)
                ? 0U : test_case.k % side) * expected_passes +
            ((test_case.bytes & 63U) != 0 ? test_case.k : 0U);
        require(route.small_transform_calls == expected_passes,
            "dense T=2/T=4 AVX2 encode missed the coarse kernel");
        require(route.input_copy_shards ==
                expected_input_copies,
            "dense T=2/T=4 AVX2 encode retained avoidable input copies");
        require(original == original_before,
            "dense T=2/T=4 AVX2 encode modified caller input");
        if ((test_case.bytes & 63U) == 0 &&
            leo_encode_work_count(test_case.k, test_case.r) != 0)
        {
            require(actual == encode_legacy(original, test_case.r),
                "dense T=2/T=4 AVX2 encode changed legacy parity bytes");
        }
        require(actual == encode_unaligned(
                    avx2_codec.get(), original, test_case.r),
            "dense T=2/T=4 AVX2 encode mishandled unaligned buffers");

        leopard::ff8::TestOnlyResetHighEncodeCounts();
        require_sparse_matches_full(
            encode(avx2_codec.get(), original, test_case.r, true), actual);
        const uint64_t sparse_calls = leopard::ff8::
            TestOnlyGetHighEncodeCounts().small_transform_calls;
        require((sparse_calls != 0) == (test_case.r == 2),
            "T=2/T=4 coarse kernel mishandled a prefix/holey output mask");

        if (test_case.k == 64 && test_case.r == 3)
        {
            leopard::ff8::TestOnlyResetHighEncodeCounts();
            const Shards prefix = encode_prefix(
                avx2_codec.get(), original, test_case.r, 2);
            require(prefix[0] == actual[0] && prefix[1] == actual[1],
                "T=4 coarse kernel changed a requested parity prefix");
            require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                        small_transform_calls == 1,
                "T=4 parity prefix missed the coarse kernel");
        }

        if (scalar.result() == LEO2_SUCCESS)
        {
            Codec scalar_codec(scalar.get(), test_case.k, test_case.r,
                LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
            require(actual == encode(
                        scalar_codec.get(), original, test_case.r, false),
                "dense T=2/T=4 AVX2 encode differs from scalar");
        }
        if (ssse3.result() == LEO2_SUCCESS)
        {
            Codec ssse3_codec(ssse3.get(), test_case.k, test_case.r,
                LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
            require(actual == encode(
                        ssse3_codec.get(), original, test_case.r, false),
                "dense T=2/T=4 AVX2 encode differs from SSSE3");
        }
        if (avx512.result() == LEO2_SUCCESS)
        {
            Codec avx512_codec(avx512.get(), test_case.k, test_case.r,
                LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
            require(actual == encode(
                        avx512_codec.get(), original, test_case.r, false),
                "dense T=2/T=4 AVX2 encode differs from AVX-512");
        }
    }

    // Lock every immediate lower K/byte boundary in the measured staircase.
    // The K<side and T=8 controls also protect the callback's preconditions.
    const TestCase controls[] = {
        { 1, 2, 65536 },
        { 2, 2, 2047 },
        { 3, 4, 2047 },
        { 4, 4, 2047 },
        { 5, 4, 2047 },
        { 6, 3, 2047 },
        { 7, 4, 2047 },
        { 8, 4, 2047 },
        { 9, 4, 2047 },
        { 11, 4, 2047 },
        { 12, 4, 2047 },
        { 16, 4, 2047 },
        { 64, 5, 65536 }
    };
    for (size_t case_i = 0;
         case_i < sizeof(controls) / sizeof(controls[0]); ++case_i)
    {
        const TestCase& test_case = controls[case_i];
        const Shards original = make_original(test_case.k, test_case.bytes);
        Codec codec(avx2.get(), test_case.k, test_case.r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        (void)encode(codec.get(), original, test_case.r, false);
        require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                    small_transform_calls == 0,
            "T=2/T=4 coarse-kernel policy escaped its shape bounds");
    }
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

    static const uint32_t gf8_sides[] = { 8, 16, 32, 64 };
    for (size_t side_i = 0;
         side_i < sizeof(gf8_sides) / sizeof(gf8_sides[0]); ++side_i)
    {
        const uint32_t side = gf8_sides[side_i];
        Codec balanced(automatic.get(), side, side,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        const size_t execution_bytes = side == 8 ? 8192 : 4096;
        if (side == 8)
        {
            require(selected_backend(balanced.get(), 1984, side, side) ==
                    LEO2_BACKEND_AVX2,
                "balanced GF8 T=8 widened below the 2 KiB singleton");
            require(selected_backend(balanced.get(), 2048, side, side) ==
                    LEO2_BACKEND_AVX512,
                "balanced GF8 T=8 did not widen at 2 KiB");
            require(selected_backend(balanced.get(), 2112, side, side) ==
                    LEO2_BACKEND_AVX2,
                "balanced GF8 T=8 widened above the 2 KiB singleton");
            require(selected_backend(balanced.get(), 4096, side, side) ==
                    LEO2_BACKEND_AVX2,
                "balanced GF8 T=8 ignored the 4 KiB exact-main veto");
            require(selected_backend(balanced.get(), 8128, side, side) ==
                    LEO2_BACKEND_AVX2,
                "balanced GF8 T=8 widened below the upper interval");
            require(selected_backend(balanced.get(), 8192, side, side) ==
                    LEO2_BACKEND_AVX512,
                "balanced GF8 T=8 did not widen at 8 KiB");
        }
        else
        {
            const uint64_t minimum_bytes = side == 16 ? 4096 : 2048;
            require(selected_backend(
                        balanced.get(), minimum_bytes - 64, side, side) ==
                    LEO2_BACKEND_AVX2,
                "balanced GF8 AUTO widened below the calibrated byte range");
            require(selected_backend(
                        balanced.get(), minimum_bytes, side, side) ==
                    LEO2_BACKEND_AVX512,
                "balanced GF8 AUTO did not widen at the lower byte boundary");
        }
        if (side == 16)
            require(selected_backend(balanced.get(), 2048, side, side) ==
                    LEO2_BACKEND_AVX2,
                "balanced GF8 T=16 widened at the inconclusive 2 KiB cell");
        require(selected_backend(
                    balanced.get(), execution_bytes, side, side) ==
                LEO2_BACKEND_AVX512,
            "balanced GF8 AUTO did not widen in a calibrated cell");
        require(selected_backend(balanced.get(), 65536, side, side) ==
                LEO2_BACKEND_AVX512,
            "balanced GF8 AUTO did not widen at the upper byte boundary");
        require(selected_backend(balanced.get(), 65600, side, side) ==
                LEO2_BACKEND_AVX2,
            "balanced GF8 AUTO widened above the calibrated byte range");
        require(selected_backend(balanced.get(), 4097, side, side) ==
                LEO2_BACKEND_AVX2,
            "balanced GF8 AUTO widened a ragged tail");
        require(selected_backend(
                    balanced.get(), execution_bytes, side - 1, side) ==
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

        Shards balanced_original = make_original(side, execution_bytes);
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
        if (side == 8)
            require_encode_overlap_rejected(
                balanced.get(), balanced_original, side);
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

    // Exercise concurrent immutable-codec use through the exact T=8 AUTO
    // route independently of the larger neighboring callback below.
    Codec concurrent_t8(automatic.get(), 8, 8,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    const Shards concurrent_t8_original = make_original(8, 8192);
    const Shards concurrent_t8_reference = encode(
        concurrent_t8.get(), concurrent_t8_original, 8, false);
    std::atomic<unsigned> concurrent_t8_failures(0);
    std::vector<std::thread> concurrent_t8_threads;
    for (unsigned lane = 0; lane < 4; ++lane)
    {
        concurrent_t8_threads.push_back(std::thread([&]() {
            try
            {
                if (encode(concurrent_t8.get(), concurrent_t8_original,
                        8, false) != concurrent_t8_reference)
                    concurrent_t8_failures.fetch_add(
                        1, std::memory_order_relaxed);
            }
            catch (...)
            {
                concurrent_t8_failures.fetch_add(
                    1, std::memory_order_relaxed);
            }
        }));
    }
    for (size_t i = 0; i < concurrent_t8_threads.size(); ++i)
        concurrent_t8_threads[i].join();
    require(concurrent_t8_failures.load(std::memory_order_relaxed) == 0,
        "concurrent exact T=8 AUTO encode failed");

    Codec gf8_t8(automatic.get(), 8, 8,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    Codec gf8_above_side(automatic.get(), 128, 128,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    Codec gf8_shortened(automatic.get(), 15, 16,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    Codec gf8_punctured(automatic.get(), 16, 15,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    require(selected_backend(gf8_t8.get(), 4096, 8, 8) ==
            LEO2_BACKEND_AVX2,
        "balanced GF8 AUTO ignored the T=8 4 KiB veto");
    require(selected_backend(gf8_t8.get(), 8192, 8, 8) ==
            LEO2_BACKEND_AVX512,
        "balanced GF8 AUTO did not widen qualified T=8");
    require(selected_backend(gf8_above_side.get(), 4096, 128, 128) ==
            LEO2_BACKEND_AVX2,
        "balanced GF8 AUTO widened T=128");
    require(selected_backend(gf8_shortened.get(), 4096, 16, 16) ==
            LEO2_BACKEND_AVX2,
        "balanced GF8 T=16 widened an unqualified shortened input block");
    require(selected_backend(gf8_punctured.get(), 4096, 15, 15) ==
            LEO2_BACKEND_AVX2,
        "balanced GF8 T=16 widened an unqualified punctured parity block");

    struct GF8CoarseCase
    {
        uint32_t k;
        uint32_t r;
    };
    static const GF8CoarseCase coarse_cases[] = {
        { 8, 8 }, { 7, 8 }, { 8, 7 }, { 7, 7 },
        { 15, 16 }, { 16, 15 }, { 15, 15 },
        { 31, 32 }, { 32, 31 }, { 31, 31 },
        { 63, 64 }, { 64, 63 }, { 63, 63 }
    };
    static const size_t coarse_bytes[] = {
        2048, 2049, 4096, 4097, 8192, 65536, 65537
    };
    for (size_t case_i = 0;
         case_i < sizeof(coarse_cases) / sizeof(coarse_cases[0]); ++case_i)
    {
        const GF8CoarseCase& current = coarse_cases[case_i];
        Codec automatic_codec(automatic.get(), current.k, current.r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        Codec scalar_codec(scalar.get(), current.k, current.r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        Codec ssse3_codec(ssse3.get(), current.k, current.r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        Codec avx2_codec(avx2.get(), current.k, current.r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
        Codec avx512_codec(avx512.get(), current.k, current.r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);

        const uint32_t side = current.r <= 8 ? 8U :
            (current.r <= 16 ? 16U : (current.r <= 32 ? 32U : 64U));
        const bool qualified_shape =
            (side == 8 && current.k == 8 && current.r == 8) ||
            (side == 64 &&
                (current.k == 64 || current.k == 63) &&
                (current.r == 64 || current.r == 63));

        // AUTO widens only the exact aligned cells that passed the isolated
        // crossover gate.  Explicit AVX-512 exercises every neighboring
        // candidate independently of the production promotion decision.
        for (size_t byte_i = 0;
             byte_i < sizeof(coarse_bytes) / sizeof(coarse_bytes[0]); ++byte_i)
        {
            const size_t bytes = coarse_bytes[byte_i];
            const bool qualified_bytes = side == 8
                ? bytes == 2U * 1024U ||
                    (bytes >= 8U * 1024U && bytes <= 64U * 1024U)
                : current.r == 64 || bytes == 64U * 1024U;
            const leo2_backend expected_auto = qualified_shape &&
                    qualified_bytes && (bytes & 63U) == 0
                ? LEO2_BACKEND_AVX512 : LEO2_BACKEND_AVX2;
            require(selected_backend(automatic_codec.get(), bytes,
                        current.r, current.r) == expected_auto,
                "GF8 coarse neighbor selected the wrong AUTO backend");
            require(selected_backend(avx512_codec.get(), bytes,
                        current.r, current.r) == LEO2_BACKEND_AVX512,
                "explicit GF8 coarse neighbor lost AVX-512");

            const Shards original = make_original(current.k, bytes);
            const Shards reference = encode(
                scalar_codec.get(), original, current.r, false);
            require(reference == encode(
                        ssse3_codec.get(), original, current.r, false) &&
                    reference == encode(
                        avx2_codec.get(), original, current.r, false),
                "GF8 coarse neighbor changed a mature backend result");

            leopard::ff8::TestOnlyResetHighEncodeCounts();
            const Shards widened = encode(
                avx512_codec.get(), original, current.r, false);
            require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                        whole_transform_calls == 1,
                "GF8 coarse neighbor did not execute one aligned-prefix callback");
            const Shards automatic_result = encode(
                automatic_codec.get(), original, current.r, false);
            require(reference == widened && reference == automatic_result,
                "GF8 coarse neighbor changed parity bytes");
            require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                        whole_transform_calls ==
                    (expected_auto == LEO2_BACKEND_AVX512 ? 2U : 1U),
                "GF8 coarse neighbor executed the wrong AUTO callback route");
            if ((bytes & 63U) == 0 && current.r <= current.k)
                require(reference == encode_legacy(original, current.r),
                    "GF8 coarse neighbor differs from legacy Leopard");

            if ((bytes & 63U) != 0)
            {
                require(reference == encode_unaligned(
                            avx2_codec.get(), original, current.r) &&
                        reference == encode_unaligned(
                            avx512_codec.get(), original, current.r),
                    "GF8 coarse neighbor mishandled an unaligned byte tail");
            }
        }

        const Shards partial_original = make_original(current.k, 4097);
        const Shards partial_reference = encode(
            scalar_codec.get(), partial_original, current.r, false);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        require_sparse_matches_full(encode(
                avx512_codec.get(), partial_original, current.r, true),
            partial_reference);
        require(leopard::ff8::TestOnlyGetHighEncodeCounts().
                    whole_transform_calls == 0,
            "partial-output GF8 neighbor used the dense coarse callback");
    }

    // One ragged, shortened-and-punctured codec is also executed concurrently;
    // each call owns its output and scratch while sharing the immutable codec.
    Codec concurrent_codec(avx512.get(), 63, 63,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8);
    const Shards concurrent_original = make_original(63, 4097);
    const Shards concurrent_reference = encode(
        concurrent_codec.get(), concurrent_original, 63, false);
    std::atomic<unsigned> coarse_failures(0);
    std::vector<std::thread> coarse_threads;
    for (unsigned lane = 0; lane < 4; ++lane)
    {
        coarse_threads.push_back(std::thread([&]() {
            try
            {
                if (encode(concurrent_codec.get(), concurrent_original,
                        63, false) != concurrent_reference)
                    coarse_failures.fetch_add(1, std::memory_order_relaxed);
            }
            catch (...)
            {
                coarse_failures.fetch_add(1, std::memory_order_relaxed);
            }
        }));
    }
    for (size_t i = 0; i < coarse_threads.size(); ++i)
        coarse_threads[i].join();
    require(coarse_failures.load(std::memory_order_relaxed) == 0,
        "concurrent GF8 coarse neighbor encode failed");

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
        test_balanced_execution_tile_geometry();

        test_small_high_encode(scalar, ssse3, avx2, avx512);

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
