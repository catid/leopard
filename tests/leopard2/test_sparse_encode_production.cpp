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

typedef std::vector<uint8_t> Bytes;
typedef std::vector<Bytes> Shards;

void Require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void RequireResult(leo2_result actual, leo2_result expected, const char* message)
{
    if (actual != expected)
        throw std::runtime_error(std::string(message) + ": " +
            leo2_result_string(actual));
}

class Context
{
public:
    explicit Context(leo2_backend backend)
        : value_(NULL), result_(LEO2_INTERNAL_ERROR)
    {
        leo2_context_options options = {};
        options.struct_size = sizeof(options);
        options.backend = backend;
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
        RequireResult(leo2_codec_create(context, k, r,
            LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16, NULL, &value_),
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
        : value_(NULL), bytes_(bytes)
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
    size_t size() const { return bytes_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);
    void* value_;
    size_t bytes_;
};

size_t Scratch(const leo2_codec* codec, uint64_t bytes)
{
    size_t result = 0;
    RequireResult(leo2_encode_scratch_size(codec, bytes, &result),
        LEO2_SUCCESS, "scratch query");
    return result;
}

Shards MakeOriginal(uint32_t k, size_t bytes)
{
    Shards result(k, Bytes(bytes, 0));
    uint32_t state = 0x9e3779b9U;
    for (uint32_t shard = 0; shard < k; ++shard)
        for (size_t i = 0; i < bytes; ++i)
        {
            state = state * 1664525U + 1013904223U;
            result[shard][i] = static_cast<uint8_t>(
                (state >> 24) ^ shard ^ static_cast<uint32_t>(i));
        }
    return result;
}

Shards Encode(
    const leo2_codec* codec,
    const Shards& original,
    const std::vector<uint8_t>& requested)
{
    const size_t bytes = original[0].size();
    std::vector<const void*> input(original.size(), NULL);
    for (size_t i = 0; i < original.size(); ++i)
        input[i] = &original[i][0];
    Shards output(requested.size(), Bytes(bytes, 0xa5));
    std::vector<void*> recovery(requested.size(), NULL);
    for (size_t i = 0; i < requested.size(); ++i)
        if (requested[i])
            recovery[i] = &output[i][0];
    AlignedBuffer scratch(Scratch(codec, bytes));
    RequireResult(leo2_encode(codec, bytes, &input[0], &recovery[0],
        scratch.get(), scratch.size()), LEO2_SUCCESS, "encode");
    return output;
}

void RequireRequestedEqual(
    const Shards& actual,
    const Shards& expected,
    const std::vector<uint8_t>& requested,
    const char* message)
{
    for (size_t i = 0; i < requested.size(); ++i)
    {
        if (requested[i])
            Require(actual[i] == expected[i], message);
        else
            Require(std::count(actual[i].begin(), actual[i].end(),
                        static_cast<uint8_t>(0xa5)) ==
                    static_cast<ptrdiff_t>(actual[i].size()),
                "unrequested parity sentinel mismatch");
    }
}

void RequireScratchDelta(
    leo2_context* candidate_context,
    leo2_context* scalar_context,
    uint32_t k,
    uint32_t r,
    uint64_t bytes,
    size_t expected_delta)
{
    Codec candidate(candidate_context, k, r);
    Codec scalar(scalar_context, k, r);
    const size_t candidate_bytes = Scratch(candidate.get(), bytes);
    const size_t scalar_bytes = Scratch(scalar.get(), bytes);
    Require(candidate_bytes >= scalar_bytes,
        "candidate scratch is unexpectedly smaller than scalar scratch");
    Require(candidate_bytes - scalar_bytes == expected_delta,
        "production sparse scratch delta mismatch");
}

} // namespace

int main()
{
    try
    {
        Context scalar(LEO2_BACKEND_SCALAR);
        Context automatic(LEO2_BACKEND_AUTO);
        RequireResult(scalar.result(), LEO2_SUCCESS, "scalar context");
        RequireResult(automatic.result(), LEO2_SUCCESS, "AUTO context");

        // Seven parity cosets each retain 448 two-bit butterfly decisions
        // (112 bytes), plus a 128-bit dependency workspace: 7*112+16 = 800.
        const size_t expected_schedule_bytes = 800;
        Require(expected_schedule_bytes <= 65536,
            "production sparse schedule exceeds its 64-KiB cap");

        Context explicit_avx2(LEO2_BACKEND_AVX2);
        if (explicit_avx2.result() == LEO2_SUCCESS)
        {
            RequireScratchDelta(explicit_avx2.get(), scalar.get(),
                128, 896, 1024, expected_schedule_bytes);
            RequireScratchDelta(explicit_avx2.get(), scalar.get(),
                128, 896, 1026, expected_schedule_bytes);
            RequireScratchDelta(explicit_avx2.get(), scalar.get(),
                128, 896, 65536, expected_schedule_bytes);
            RequireScratchDelta(explicit_avx2.get(), scalar.get(),
                128, 896, 64, 0);
            RequireScratchDelta(explicit_avx2.get(), scalar.get(),
                128, 896, 1022, 0);
            RequireScratchDelta(explicit_avx2.get(), scalar.get(),
                127, 896, 1024, 0);
            RequireScratchDelta(explicit_avx2.get(), scalar.get(),
                129, 896, 1024, 0);
            RequireScratchDelta(explicit_avx2.get(), scalar.get(),
                128, 895, 1024, 0);
            RequireScratchDelta(explicit_avx2.get(), scalar.get(),
                128, 897, 1024, 0);
        }
        else
            RequireResult(explicit_avx2.result(), LEO2_UNSUPPORTED,
                "explicit AVX2 context");

        const leo2_backend auto_backend =
            leo2_context_backend(automatic.get());
        RequireScratchDelta(automatic.get(), scalar.get(), 128, 896, 1024,
            auto_backend == LEO2_BACKEND_AVX2 ? expected_schedule_bytes : 0);

        const uint32_t k = 128;
        const uint32_t r = 896;
        const size_t bytes = 1026;
        const unsigned edge_indices[] = { 0, 127, 128, 383, 895 };
        const unsigned scattered_indices[] = { 7, 63, 135, 255, 519, 895 };
        std::vector<uint8_t> edge(r, 0);
        std::vector<uint8_t> scattered(r, 0);
        for (size_t i = 0;
             i < sizeof(edge_indices) / sizeof(edge_indices[0]); ++i)
            edge[edge_indices[i]] = 1;
        for (size_t i = 0;
             i < sizeof(scattered_indices) / sizeof(scattered_indices[0]); ++i)
            scattered[scattered_indices[i]] = 1;
        const Shards original = MakeOriginal(k, bytes);
        Codec scalar_codec(scalar.get(), k, r);
        Codec auto_codec(automatic.get(), k, r);
        const Shards edge_reference = Encode(
            scalar_codec.get(), original, edge);
        const Shards scattered_reference = Encode(
            scalar_codec.get(), original, scattered);
        RequireRequestedEqual(Encode(auto_codec.get(), original, edge),
            edge_reference, edge, "default AUTO edge parity mismatch");
        RequireRequestedEqual(Encode(auto_codec.get(), original, scattered),
            scattered_reference, scattered,
            "default AUTO scattered parity mismatch");

        // A one-byte-short span is rejected before call-local plan compilation
        // and before any output write in the ordinary, no-hook archive.
        std::vector<const void*> input(k, NULL);
        for (uint32_t i = 0; i < k; ++i)
            input[i] = &original[i][0];
        Shards rejected(r, Bytes(bytes, 0x6d));
        const Shards rejected_before = rejected;
        std::vector<void*> recovery(r, NULL);
        for (uint32_t i = 0; i < r; ++i)
            if (edge[i])
                recovery[i] = &rejected[i][0];
        const size_t scratch_bytes = Scratch(auto_codec.get(), bytes);
        AlignedBuffer scratch(scratch_bytes);
        RequireResult(leo2_encode(auto_codec.get(), bytes, &input[0],
            &recovery[0], scratch.get(), scratch.size() - 1),
            LEO2_SCRATCH_TOO_SMALL, "one-byte-short production encode");
        Require(rejected == rejected_before,
            "one-byte-short production encode modified output");

        std::printf(
            "Leopard2 production sparse encode passed: auto_backend=%u "
            "schedule_delta=%zu schedule_cap=65536 scratch_cells=10\n",
            static_cast<unsigned>(auto_backend), expected_schedule_bytes);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr,
            "Leopard2 production sparse encode failed: %s\n", error.what());
        return 1;
    }
}
