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
#include "leopard.h"
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
        const size_t allocated = bytes == 0 ? leo2_scratch_alignment() : bytes;
#if defined(_MSC_VER)
        value_ = _aligned_malloc(allocated, leo2_scratch_alignment());
#else
        if (posix_memalign(
                &value_, leo2_scratch_alignment(), allocated) != 0)
            value_ = NULL;
#endif
        if (!value_)
            throw std::bad_alloc();
        std::memset(value_, 0, allocated);
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

void FillInput(uint8_t* data, size_t bytes, uint64_t seed)
{
    uint64_t state = seed;
    for (size_t i = 0; i < bytes; ++i)
    {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        data[i] = static_cast<uint8_t>(state >> 24);
    }
}

void SetPackedPointers(
    uint8_t* input,
    uint8_t* output,
    unsigned count,
    size_t shard_bytes,
    std::vector<const void*>& original,
    std::vector<void*>& recovery)
{
    original.resize(count);
    recovery.resize(count);
    for (unsigned i = 0; i < count; ++i)
    {
        original[i] = input + static_cast<size_t>(i) * shard_bytes;
        recovery[i] = output + static_cast<size_t>(i) * shard_bytes;
    }
}

void CheckParity(
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator,
    const std::vector<const void*>& original,
    const std::vector<void*>& recovery,
    size_t shard_bytes)
{
    const unsigned count = static_cast<unsigned>(original.size());
    for (unsigned parity = 0; parity < count; ++parity)
    {
        if (!recovery[parity])
            continue;
        const std::vector<leopard2_test::Element>& row =
            generator[count + parity];
        const uint8_t* const output =
            static_cast<const uint8_t*>(recovery[parity]);
        for (size_t offset = 0; offset < shard_bytes; ++offset)
        {
            leopard2_test::Element expected = 0;
            for (unsigned source = 0; source < count; ++source)
            {
                expected = field.add(expected, field.multiply(row[source],
                    static_cast<const uint8_t*>(
                        original[source])[offset]));
            }
            Require(output[offset] == expected,
                "prepared T16 parity differs from direct oracle");
        }
    }
}

uint64_t PreparedCalls()
{
    return leopard::ff8::TestOnlyGetHighEncodeCounts().t16_prepared_calls;
}

leo2_codec* CreateCodec(leo2_context* context, unsigned count)
{
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, count, count,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        LEO2_SUCCESS, "create balanced T16 codec");
    return codec;
}

void CompareLegacy(
    unsigned count,
    size_t shard_bytes,
    const std::vector<const void*>& original,
    const std::vector<void*>& recovery)
{
    const unsigned work_count = leo_encode_work_count(count, count);
    Require(work_count >= count, "legacy T16 work count is too small");
    AlignedBuffer legacy_storage(
        static_cast<size_t>(work_count) * shard_bytes);
    std::vector<void*> legacy(work_count);
    for (unsigned i = 0; i < work_count; ++i)
        legacy[i] = legacy_storage.bytes() +
            static_cast<size_t>(i) * shard_bytes;
    Require(leo_encode(shard_bytes, count, count, work_count,
            &original[0], &legacy[0]) == Leopard_Success,
        "legacy T16 encode failed");
    for (unsigned parity = 0; parity < count; ++parity)
    {
        Require(std::memcmp(
                recovery[parity], legacy[parity], shard_bytes) == 0,
            "prepared T16 parity differs from legacy Leopard");
    }
}

void ExerciseSelectedCell(
    leo2_context* context,
    unsigned count,
    size_t shard_bytes,
    size_t input_offset,
    size_t output_offset)
{
    leo2_codec* codec = CreateCodec(context, count);
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, shard_bytes, &scratch_bytes), LEO2_SUCCESS,
        "query prepared T16 scratch");
    AlignedBuffer scratch(scratch_bytes);
    AlignedBuffer input(
        static_cast<size_t>(count) * shard_bytes + input_offset + 8U);
    AlignedBuffer output(
        static_cast<size_t>(count) * shard_bytes + output_offset + 8U);
    FillInput(input.bytes() + input_offset,
        static_cast<size_t>(count) * shard_bytes,
        UINT64_C(0x5431365052455000) ^ count ^ shard_bytes);
    std::vector<const void*> original;
    std::vector<void*> recovery;
    SetPackedPointers(input.bytes() + input_offset,
        output.bytes() + output_offset, count, shard_bytes,
        original, recovery);

    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh, count, count);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);

    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute prepared T16 encode");
    Require(PreparedCalls() == 1,
        "qualified packed T16 encode missed prepared terminal");
    CheckParity(field, generator, original, recovery, shard_bytes);
    CompareLegacy(count, shard_bytes, original, recovery);

    std::memset(output.bytes() + output_offset, 0xa5,
        static_cast<size_t>(count) * shard_bytes);
    leo2_encode_batch_item item = {
        shard_bytes, &original[0], &recovery[0], scratch.data(), scratch.size()
    };
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch(codec, &item, 1), LEO2_SUCCESS,
        "execute one-item prepared T16 batch");
    Require(PreparedCalls() == 1,
        "qualified T16 one-item batch missed prepared terminal");
    CheckParity(field, generator, original, recovery, shard_bytes);

    leo2_encode_batch_binding* binding = NULL;
    RequireResult(leo2_encode_batch_binding_create(
        codec, &item, 1, &binding), LEO2_SUCCESS,
        "create prepared T16 binding");
    std::memset(output.bytes() + output_offset, 0xa5,
        static_cast<size_t>(count) * shard_bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch_binding_execute(binding),
        LEO2_SUCCESS, "execute prepared T16 binding");
    Require(PreparedCalls() == 1,
        "qualified T16 binding missed prepared terminal");
    CheckParity(field, generator, original, recovery, shard_bytes);
    leo2_encode_batch_binding_destroy(binding);
    leo2_codec_destroy(codec);
}

void ExerciseExcludedCell(
    leo2_context* context,
    unsigned original_count,
    unsigned recovery_count,
    size_t shard_bytes,
    bool detach_input,
    bool sparse_output)
{
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, original_count, recovery_count,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        LEO2_SUCCESS, "create excluded T16 codec");
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, shard_bytes, &scratch_bytes), LEO2_SUCCESS,
        "query excluded T16 scratch");
    AlignedBuffer scratch(scratch_bytes);
    AlignedBuffer input(
        static_cast<size_t>(original_count) * shard_bytes);
    AlignedBuffer output(
        static_cast<size_t>(recovery_count) * shard_bytes);
    FillInput(input.bytes(), input.size(),
        UINT64_C(0x5431364558434c00) ^ original_count ^ recovery_count ^
            shard_bytes);
    std::vector<const void*> original(original_count);
    std::vector<void*> recovery(recovery_count);
    for (unsigned i = 0; i < original_count; ++i)
        original[i] = input.bytes() + static_cast<size_t>(i) * shard_bytes;
    for (unsigned i = 0; i < recovery_count; ++i)
        recovery[i] = output.bytes() + static_cast<size_t>(i) * shard_bytes;
    std::vector<uint8_t> detached;
    if (detach_input)
    {
        detached.resize(shard_bytes);
        std::memcpy(&detached[0], original[original_count / 2], shard_bytes);
        original[original_count / 2] = &detached[0];
    }
    if (sparse_output)
        recovery[recovery_count / 2] = NULL;

    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute excluded T16 cell");
    Require(PreparedCalls() == 0,
        "excluded T16 cell entered prepared terminal");

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
        const uint8_t* const parity_bytes =
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
            Require(parity_bytes[offset] == expected,
                "excluded T16 parity differs from direct oracle");
        }
    }
    leo2_codec_destroy(codec);
}

void ExerciseDiagnosticControl(leo2_context* context)
{
    Require(leopard2_internal::
            SetHighT16PreparedTerminalEnabledForDiagnostics(false),
        "disable prepared T16 terminal");
    ExerciseExcludedCell(context, 12, 12, 1024, false, false);
    Require(leopard2_internal::
            SetHighT16PreparedTerminalEnabledForDiagnostics(true),
        "restore prepared T16 terminal");
}

void ExerciseBackendExclusion(leo2_backend requested_backend)
{
    leo2_context_options options = {};
    options.struct_size = sizeof(options);
    options.backend = requested_backend;
    options.thread_count = 1;
    leo2_context* context = NULL;
    const leo2_result result = leo2_context_create(&options, &context);
    if (result == LEO2_UNSUPPORTED)
        return;
    RequireResult(result, LEO2_SUCCESS, "create excluded backend context");
    ExerciseExcludedCell(context, 12, 12, 1024, false, false);
    leo2_context_destroy(context);
}

} // namespace

int main()
{
    try
    {
        Require(leo_init() == Leopard_Success,
            "legacy initialization failed");
        leo2_context_options options = {};
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AVX2;
        options.thread_count = 1;
        leo2_context* context = NULL;
        const leo2_result result = leo2_context_create(&options, &context);
        if (result == LEO2_UNSUPPORTED)
        {
            std::printf("prepared T16 AVX2 test skipped: AVX2 unavailable\n");
            return 0;
        }
        RequireResult(result, LEO2_SUCCESS,
            "create prepared T16 AVX2 context");

        static const size_t selected_bytes[] = {
            64, 128, 192, 256, 320, 384, 448, 512, 1024, 2048
        };
        for (unsigned count = 9; count <= 16; ++count)
        {
            for (size_t byte_index = 0;
                 byte_index < sizeof(selected_bytes) /
                     sizeof(selected_bytes[0]); ++byte_index)
            {
                if (count == 16 && selected_bytes[byte_index] == 64)
                    continue;
                ExerciseSelectedCell(context, count,
                    selected_bytes[byte_index], 0, 0);
            }
        }
        ExerciseSelectedCell(context, 9, 64, 1, 3);
        ExerciseSelectedCell(context, 9, 512, 1, 3);

        ExerciseExcludedCell(context, 12, 12, 63, false, false);
        ExerciseExcludedCell(context, 12, 12, 65, false, false);
        ExerciseExcludedCell(context, 12, 12, 511, false, false);
        ExerciseExcludedCell(context, 12, 12, 513, false, false);
        ExerciseExcludedCell(context, 16, 16, 64, false, false);
        ExerciseExcludedCell(context, 12, 12, 2112, false, false);
        ExerciseExcludedCell(context, 8, 8, 1024, false, false);
        ExerciseExcludedCell(context, 9, 8, 1024, false, false);
        ExerciseExcludedCell(context, 9, 10, 1024, false, false);
        ExerciseExcludedCell(context, 12, 12, 1024, true, false);
        ExerciseExcludedCell(context, 12, 12, 1024, false, true);
        ExerciseDiagnosticControl(context);

        leo2_context_destroy(context);
        ExerciseBackendExclusion(LEO2_BACKEND_SCALAR);
        ExerciseBackendExclusion(LEO2_BACKEND_SSSE3);
        std::printf("prepared T16 AVX2 terminal checks passed\n");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "prepared T16 terminal failure: %s\n",
            error.what());
        return 1;
    }
}
