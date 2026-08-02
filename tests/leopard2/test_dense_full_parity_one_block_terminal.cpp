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
#include "Leopard2Direct.h"
#include "LeopardFF8.h"
#include "leopard2.h"

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

void Require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void RequireResult(
    leo2_result actual,
    leo2_result expected,
    const char* message)
{
    if (actual != expected)
    {
        throw std::runtime_error(std::string(message) + ": expected " +
            leo2_result_string(expected) + ", received " +
            leo2_result_string(actual));
    }
}

class AlignedBuffer
{
public:
    explicit AlignedBuffer(size_t bytes)
        : value_(NULL)
        , bytes_(bytes)
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

    uint8_t* bytes() const { return static_cast<uint8_t*>(value_); }
    void* data() const { return value_; }
    size_t size() const { return bytes_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);
    void* value_;
    size_t bytes_;
};

void FillInput(uint8_t* input, unsigned original_count,
    size_t shard_bytes, uint64_t seed)
{
    uint64_t state = seed ^ original_count ^
        (static_cast<uint64_t>(shard_bytes) << 32);
    for (size_t i = 0;
         i < static_cast<size_t>(original_count) * shard_bytes; ++i)
    {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        input[i] = static_cast<uint8_t>(state >> 24);
    }
}

void SetPackedPointers(
    uint8_t* input,
    uint8_t* output,
    unsigned original_count,
    unsigned recovery_count,
    size_t shard_bytes,
    std::vector<const void*>& original,
    std::vector<void*>& recovery)
{
    original.resize(original_count);
    recovery.resize(recovery_count);
    for (unsigned i = 0; i < original_count; ++i)
        original[i] = input + static_cast<size_t>(i) * shard_bytes;
    for (unsigned i = 0; i < recovery_count; ++i)
        recovery[i] = output + static_cast<size_t>(i) * shard_bytes;
}

void CheckParity(
    unsigned original_count,
    unsigned recovery_count,
    size_t shard_bytes,
    const std::vector<const void*>& original,
    const std::vector<void*>& recovery)
{
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh, original_count, recovery_count);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);
    for (unsigned parity = 0; parity < recovery_count; ++parity)
    {
        if (!recovery[parity])
            continue;
        const std::vector<leopard2_test::Element>& row =
            generator[original_count + parity];
        const uint8_t* const encoded =
            static_cast<const uint8_t*>(recovery[parity]);
        for (size_t offset = 0; offset < shard_bytes; ++offset)
        {
            leopard2_test::Element expected = 0;
            for (unsigned source = 0; source < original_count; ++source)
            {
                expected = field.add(expected, field.multiply(row[source],
                    static_cast<const uint8_t*>(
                        original[source])[offset]));
            }
            Require(encoded[offset] == expected,
                "dense full-parity result differs from direct oracle");
        }
    }
}

void ResetCandidateCalls()
{
#ifdef LEO2_ENABLE_TEST_HOOKS
    leopard::ff8::TestOnlyResetHighEncodeCounts();
#endif
}

void RequireCandidateCalls(uint64_t expected)
{
#ifdef LEO2_ENABLE_TEST_HOOKS
    const uint64_t actual = leopard::ff8::TestOnlyGetHighEncodeCounts().
        dense_full_parity_one_block_calls;
    Require(actual == expected,
        "dense full-parity route count disagrees with selector");
#else
    (void)expected;
#endif
}

void RequireBalancedCalls(uint64_t expected)
{
#ifdef LEO2_ENABLE_TEST_HOOKS
    const uint64_t actual = leopard::ff8::TestOnlyGetHighEncodeCounts().
        balanced_b64_packed_calls;
    Require(actual == expected,
        "balanced B64 route count disagrees with calibration precedence");
#else
    (void)expected;
#endif
}

#ifndef LEO2_ENABLE_TEST_HOOKS
size_t ExpectedProductionScratch(
    unsigned original_count,
    unsigned side,
    size_t shard_bytes)
{
    const size_t alignment = leo2_scratch_alignment();
    const size_t metadata_bytes =
        2U * (original_count + side) * sizeof(uintptr_t) +
        (original_count + 2U * side) * sizeof(void*);
    const size_t data_offset =
        (metadata_bytes + alignment - 1U) & ~(alignment - 1U);
    return data_offset + 2U * side * shard_bytes;
}
#endif

enum LayoutMutation
{
    kPacked,
    kDetachedInput,
    kDetachedOutput,
    kSparseOutput
};

void ExerciseSuccess(
    leo2_context* context,
    unsigned original_count,
    unsigned recovery_count,
    size_t shard_bytes,
    uint64_t expected_calls,
    LayoutMutation mutation = kPacked,
    bool batch = false,
    size_t input_offset = 0,
    size_t output_offset = 0)
{
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context,
        original_count, recovery_count, LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF8, NULL, &codec), LEO2_SUCCESS,
        "create dense full-parity codec");
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, shard_bytes, &scratch_bytes), LEO2_SUCCESS,
        "query dense full-parity scratch");
    AlignedBuffer scratch(scratch_bytes);
    AlignedBuffer input(
        static_cast<size_t>(original_count) * shard_bytes + input_offset + 8U);
    AlignedBuffer output(
        static_cast<size_t>(recovery_count) * shard_bytes + output_offset + 8U);
    FillInput(input.bytes() + input_offset, original_count, shard_bytes,
        UINT64_C(0x44454e534546554c) ^ recovery_count);
    std::vector<const void*> original;
    std::vector<void*> recovery;
    SetPackedPointers(input.bytes() + input_offset,
        output.bytes() + output_offset, original_count, recovery_count,
        shard_bytes, original, recovery);

    std::vector<uint8_t> detached_input;
    std::vector<uint8_t> detached_output;
    if (mutation == kDetachedInput)
    {
        detached_input.resize(shard_bytes);
        std::memcpy(&detached_input[0], original[original_count / 2U],
            shard_bytes);
        original[original_count / 2U] = &detached_input[0];
    }
    else if (mutation == kDetachedOutput)
    {
        detached_output.resize(shard_bytes);
        recovery[recovery_count / 2U] = &detached_output[0];
    }
    else if (mutation == kSparseOutput)
        recovery[recovery_count / 2U] = NULL;

    ResetCandidateCalls();
    leo2_result result;
    if (batch)
    {
        leo2_encode_batch_item item = {
            shard_bytes, &original[0], &recovery[0],
            scratch.data(), scratch.size()
        };
        result = leo2_encode_batch(codec, &item, 1);
    }
    else
    {
        result = leo2_encode(codec, shard_bytes, &original[0], &recovery[0],
            scratch.data(), scratch.size());
    }
    RequireResult(result, LEO2_SUCCESS, "execute dense full-parity case");
    RequireCandidateCalls(expected_calls);
    CheckParity(original_count, recovery_count, shard_bytes,
        original, recovery);
    leo2_codec_destroy(codec);
}

void ExerciseOverlap(leo2_context* context, bool scratch_overlap)
{
    static const unsigned k = 31;
    static const unsigned r = 32;
    static const size_t bytes = 64;
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        LEO2_SUCCESS, "create overlap codec");
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        LEO2_SUCCESS, "query overlap scratch");
    AlignedBuffer scratch(scratch_bytes);
    AlignedBuffer data(static_cast<size_t>(r) * bytes);
    FillInput(data.bytes(), k, bytes, UINT64_C(0x4f5645524c415031));
    std::vector<const void*> original;
    std::vector<void*> recovery;
    SetPackedPointers(data.bytes(),
        scratch_overlap ? scratch.bytes() : data.bytes(),
        k, r, bytes, original, recovery);

    ResetCandidateCalls();
    const leo2_result result = leo2_encode(codec, bytes,
        &original[0], &recovery[0], scratch.data(), scratch.size());
    RequireResult(result, LEO2_OVERLAP,
        scratch_overlap ? "reject scratch/output overlap" :
            "reject input/output overlap");
    RequireCandidateCalls(0);
    leo2_codec_destroy(codec);
}

void ExerciseScratchContract(leo2_context* context)
{
    static const unsigned k = 31;
    static const unsigned r = 32;
    static const size_t bytes = 64;
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        LEO2_SUCCESS, "create scratch-contract codec");
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        LEO2_SUCCESS, "query scratch-contract workspace");
#ifndef LEO2_ENABLE_TEST_HOOKS
    Require(scratch_bytes == ExpectedProductionScratch(k, r, bytes),
        "dense full-parity fixed scratch geometry differs from query");
#endif
    AlignedBuffer scratch(scratch_bytes + leo2_scratch_alignment());
    AlignedBuffer input(static_cast<size_t>(k) * bytes);
    AlignedBuffer output(static_cast<size_t>(r) * bytes);
    FillInput(input.bytes(), k, bytes, UINT64_C(0x5343524154434831));
    std::memset(output.bytes(), 0xa5, output.size());
    std::vector<const void*> original;
    std::vector<void*> recovery;
    SetPackedPointers(input.bytes(), output.bytes(), k, r, bytes,
        original, recovery);
    const std::vector<uint8_t> output_before(
        output.bytes(), output.bytes() + output.size());

    ResetCandidateCalls();
    RequireResult(leo2_encode(codec, bytes, &original[0], &recovery[0],
        scratch.data(), scratch_bytes - 1U), LEO2_SCRATCH_TOO_SMALL,
        "reject undersized dense full-parity scratch");
    RequireCandidateCalls(0);
    Require(std::memcmp(output.bytes(), &output_before[0],
            output_before.size()) == 0,
        "undersized scratch modified dense full-parity output");

    ResetCandidateCalls();
    RequireResult(leo2_encode(codec, bytes, &original[0], &recovery[0],
        scratch.bytes() + 1U, scratch_bytes), LEO2_BAD_ALIGNMENT,
        "reject misaligned dense full-parity scratch");
    RequireCandidateCalls(0);
    Require(std::memcmp(output.bytes(), &output_before[0],
            output_before.size()) == 0,
        "misaligned scratch modified dense full-parity output");
    leo2_codec_destroy(codec);
}

leo2_context* CreateContext(leo2_backend backend)
{
    leo2_context_options options = {};
    options.struct_size = sizeof(options);
    options.backend = backend;
    options.thread_count = 1;
    leo2_context* context = NULL;
    RequireResult(leo2_context_create(&options, &context), LEO2_SUCCESS,
        "create dense full-parity context");
    return context;
}

} // namespace

int main()
{
    try
    {
        leo2_context_options options = {};
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AVX2;
        options.thread_count = 1;
        leo2_context* context = NULL;
        const leo2_result context_result =
            leo2_context_create(&options, &context);
        if (context_result == LEO2_UNSUPPORTED)
        {
            std::printf(
                "dense full-parity AVX2 terminal test skipped: AVX2 unavailable\n");
            return 0;
        }
        RequireResult(context_result, LEO2_SUCCESS,
            "create AVX2 dense full-parity context");

        Require(!leopard2_internal::
                SetDenseFullParityOneBlockModeForDiagnostics(3),
            "dense full-parity diagnostic selector accepted invalid mode");

        /* Production starts fail-closed until the all-K table is calibrated. */
        Require(leopard2_internal::
                SetDenseFullParityOneBlockModeForDiagnostics(0),
            "select dense full-parity production policy");
        ExerciseSuccess(context, 15, 16, 64, 0);
        RequireBalancedCalls(0);
        ExerciseSuccess(context, 32, 32, 64, 0);
        RequireBalancedCalls(1);

        Require(leopard2_internal::
                SetDenseFullParityOneBlockModeForDiagnostics(1),
            "force dense full-parity candidate envelope");
        ExerciseSuccess(context, 32, 32, 64, 1);
        RequireBalancedCalls(0);
        const unsigned boundary[][2] = {
            { 1, 16 }, { 15, 16 }, { 16, 16 },
            { 1, 32 }, { 31, 32 }, { 32, 32 },
            { 63, 64 }, { 64, 64 },
            { 127, 128 }, { 128, 128 }
        };
        for (size_t i = 0; i < sizeof(boundary) / sizeof(boundary[0]); ++i)
        {
            ExerciseSuccess(context, boundary[i][0], boundary[i][1],
                (i & 1U) == 0 ? 64U : 256U, 1);
        }
        ExerciseSuccess(context, 31, 32, 64, 1, kPacked, true);
        ExerciseSuccess(context, 31, 32, 64, 1, kPacked, false, 1, 3);

        /* Immediate count/byte and non-packed neighbors must fall through. */
        ExerciseSuccess(context, 31, 31, 64, 0);
        const size_t byte_neighbors[] = { 63, 65, 255, 257, 1024 };
        for (size_t i = 0;
             i < sizeof(byte_neighbors) / sizeof(byte_neighbors[0]); ++i)
            ExerciseSuccess(context, 31, 32, byte_neighbors[i], 0);
        ExerciseSuccess(context, 33, 32, 64, 0);
        ExerciseSuccess(context, 31, 32, 64, 0, kDetachedInput);
        ExerciseSuccess(context, 31, 32, 64, 0, kDetachedOutput);
        ExerciseSuccess(context, 31, 32, 64, 0, kSparseOutput);
        ExerciseOverlap(context, false);
        ExerciseOverlap(context, true);
        ExerciseScratchContract(context);

        /* Same executable, explicit mature route. */
        Require(leopard2_internal::
                SetDenseFullParityOneBlockModeForDiagnostics(2),
            "force dense full-parity control route");
        ExerciseSuccess(context, 31, 32, 64, 0);
        RequireBalancedCalls(0);
        ExerciseSuccess(context, 32, 32, 64, 0);
        RequireBalancedCalls(1);
        leo2_context_destroy(context);

        /* A forced candidate mode never overrides the selected backend. */
        Require(leopard2_internal::
                SetDenseFullParityOneBlockModeForDiagnostics(1),
            "restore dense full-parity candidate envelope");
        leo2_context* scalar = CreateContext(LEO2_BACKEND_SCALAR);
        ExerciseSuccess(scalar, 31, 32, 64, 0);
        leo2_context_destroy(scalar);

        Require(leopard2_internal::
                SetDenseFullParityOneBlockModeForDiagnostics(0),
            "restore dense full-parity production policy");
        std::printf("dense full-parity one-block terminal tests passed\n");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "dense full-parity terminal test failed: %s\n",
            error.what());
        return 1;
    }
}
