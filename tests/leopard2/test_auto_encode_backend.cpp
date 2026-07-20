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
    Codec(leo2_context* context, uint32_t k, uint32_t r)
        : value_(NULL)
    {
        require_result(leo2_codec_create(context, k, r,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, NULL, &value_),
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

void require_sparse_matches_full(
    const Shards& sparse, const Shards& full)
{
    require(sparse.size() == full.size(), "sparse parity size mismatch");
    for (size_t i = 0; i < full.size(); i += 2)
        require(sparse[i] == full[i], "sparse parity differs from full parity");
}

void test_selection_and_bytes(
    Context& automatic, Context& avx2, Context& avx512)
{
    Codec auto_codec(automatic.get(), 1000, 200);
    Codec avx2_codec(avx2.get(), 1000, 200);
    Codec avx512_codec(avx512.get(), 1000, 200);

    require(selected_backend(auto_codec.get(), 32, 200, 200) ==
            LEO2_BACKEND_AVX2,
        "AUTO widened a partial transform tile");
    require(selected_backend(auto_codec.get(), 64, 200, 200) ==
            LEO2_BACKEND_AVX512,
        "AUTO did not widen at the complete-tile boundary");
    require(selected_backend(auto_codec.get(), 4098, 200, 200) ==
            LEO2_BACKEND_AVX2,
        "AUTO widened an uncalibrated tail");
    require(selected_backend(auto_codec.get(), 4U * 1024U * 1024U,
                200, 200) ==
            LEO2_BACKEND_AVX512,
        "AUTO did not widen a calibrated large complete tile");
    require(selected_backend(auto_codec.get(), 4096, 199, 200) ==
            LEO2_BACKEND_AVX2,
        "AUTO widened a partial-output encode");
    require(selected_backend(avx2_codec.get(), 4096, 200, 200) ==
            LEO2_BACKEND_AVX2,
        "explicit AVX2 widened");
    require(selected_backend(avx512_codec.get(), 4032, 200, 200) ==
            LEO2_BACKEND_AVX512,
        "explicit AVX512 did not remain exact");

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

    Codec tiny(automatic.get(), 8, 8);
    require(selected_backend(tiny.get(), 64, 8, 8) ==
            LEO2_BACKEND_AVX512,
        "tiny GF16 high transform did not widen");
    Codec maximum_parent(automatic.get(), 60000, 1024);
    require(selected_backend(maximum_parent.get(), 1024, 1024, 1024) ==
            LEO2_BACKEND_AVX512,
        "maximum-parent GF16 high transform did not widen");
    Codec direct_xor(automatic.get(), 129, 1);
    require(selected_backend(direct_xor.get(), 4096, 1, 1) ==
            LEO2_BACKEND_AVX2,
        "T=1 direct codec unnecessarily qualified transform widening");

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
        Context avx2(LEO2_BACKEND_AVX2);
        Context avx512(LEO2_BACKEND_AVX512);

        if (leopard::backend::IsCalibratedAutoAVX512EncodeHost() &&
            leo2_context_backend(automatic.get()) == LEO2_BACKEND_AVX2 &&
            avx2.result() == LEO2_SUCCESS &&
            avx512.result() == LEO2_SUCCESS)
        {
            test_selection_and_bytes(automatic, avx2, avx512);
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
