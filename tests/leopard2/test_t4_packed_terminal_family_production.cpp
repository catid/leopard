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
#include "leopard2.h"

#include <cstdint>
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

struct Cell
{
    unsigned original_count;
    unsigned recovery_count;
    size_t shard_bytes;
};

static const size_t kGuardBytes = 64;
static const uint8_t kGuardValue = 0x6d;

std::string CellLabel(const Cell& cell)
{
    char text[96];
    std::snprintf(text, sizeof(text), "K=%u/R=%u/B=%zu",
        cell.original_count, cell.recovery_count, cell.shard_bytes);
    return text;
}

void RequireCell(bool condition, const Cell& cell, const char* message)
{
    if (!condition)
        throw std::runtime_error(CellLabel(cell) + ": " + message);
}

void RequireResult(
    leo2_result actual,
    leo2_result expected,
    const Cell& cell,
    const char* message)
{
    if (actual != expected)
    {
        throw std::runtime_error(CellLabel(cell) + ": " + message +
            ": expected " + leo2_result_string(expected) + ", received " +
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

void FillInput(uint8_t* input, const Cell& cell)
{
    uint64_t state = UINT64_C(0x543450524f445543) ^
        (static_cast<uint64_t>(cell.original_count) << 48) ^
        (static_cast<uint64_t>(cell.recovery_count) << 40) ^
        cell.shard_bytes;
    const size_t total = static_cast<size_t>(cell.original_count) *
        cell.shard_bytes;
    for (size_t i = 0; i < total; ++i)
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
    const Cell& cell,
    std::vector<const void*>& original,
    std::vector<void*>& recovery)
{
    for (unsigned i = 0; i < cell.original_count; ++i)
        original[i] = input + static_cast<size_t>(i) * cell.shard_bytes;
    for (unsigned i = 0; i < cell.recovery_count; ++i)
        recovery[i] = output + static_cast<size_t>(i) * cell.shard_bytes;
}

leopard2_test::Matrix MakeGenerator(const Cell& cell)
{
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh,
            cell.original_count, cell.recovery_count);
    return leopard2_test::direct_systematic_generator(field, layout);
}

void CheckParity(
    const Cell& cell,
    const leopard2_test::Matrix& generator,
    const std::vector<const void*>& original,
    const std::vector<void*>& recovery)
{
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    for (unsigned parity = 0; parity < cell.recovery_count; ++parity)
    {
        const uint8_t* const output =
            static_cast<const uint8_t*>(recovery[parity]);
        const std::vector<leopard2_test::Element>& row =
            generator[cell.original_count + parity];
        for (size_t offset = 0; offset < cell.shard_bytes; ++offset)
        {
            leopard2_test::Element expected = 0;
            for (unsigned source = 0; source < cell.original_count; ++source)
            {
                const uint8_t value =
                    static_cast<const uint8_t*>(original[source])[offset];
                expected = field.add(expected,
                    field.multiply(row[source], value));
            }
            RequireCell(output[offset] == expected, cell,
                "production parity differs from the independent oracle");
        }
    }
}

void CheckOutputGuards(
    const AlignedBuffer& output,
    size_t payload_offset,
    const Cell& cell)
{
    const size_t payload_bytes =
        static_cast<size_t>(cell.recovery_count) * cell.shard_bytes;
    for (size_t i = 0; i < payload_offset; ++i)
        RequireCell(output.bytes()[i] == kGuardValue, cell,
            "production encode modified a leading destination guard");
    for (size_t i = payload_offset + payload_bytes;
         i < output.size(); ++i)
    {
        RequireCell(output.bytes()[i] == kGuardValue, cell,
            "production encode modified a trailing destination guard");
    }
}

void ExerciseCell(leo2_context* context, const Cell& cell)
{
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context,
        cell.original_count, cell.recovery_count,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        NULL, &codec), LEO2_SUCCESS, cell, "create production codec");

    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, cell.shard_bytes, &scratch_bytes),
        LEO2_SUCCESS, cell, "query exact production scratch");
    RequireCell(scratch_bytes != 0, cell,
        "production scratch query returned an empty extent");
    AlignedBuffer scratch(scratch_bytes);
    RequireCell(reinterpret_cast<uintptr_t>(scratch.data()) %
            leo2_scratch_alignment() == 0,
        cell, "production scratch allocation is misaligned");

    const size_t input_payload_bytes =
        static_cast<size_t>(cell.original_count) * cell.shard_bytes;
    const size_t output_payload_bytes =
        static_cast<size_t>(cell.recovery_count) * cell.shard_bytes;
    AlignedBuffer input(input_payload_bytes + 2U * kGuardBytes + 8U);
    AlignedBuffer output(output_payload_bytes + 2U * kGuardBytes + 8U);
    std::memset(input.bytes(), kGuardValue, input.size());
    std::memset(output.bytes(), kGuardValue, output.size());

    const size_t input_offset = kGuardBytes + 1U;
    const size_t output_offset = kGuardBytes + 3U;
    uint8_t* const input_base = input.bytes() + input_offset;
    uint8_t* const output_base = output.bytes() + output_offset;
    FillInput(input_base, cell);
    const std::vector<uint8_t> input_before(
        input.bytes(), input.bytes() + input.size());
    std::vector<const void*> original(cell.original_count);
    std::vector<void*> recovery(cell.recovery_count);
    SetPackedPointers(input_base, output_base, cell, original, recovery);
    const leopard2_test::Matrix generator = MakeGenerator(cell);

    RequireResult(leo2_encode(codec, cell.shard_bytes,
        &original[0], &recovery[0], scratch.data(), scratch.size()),
        LEO2_SUCCESS, cell, "execute production packed terminal");
    RequireCell(std::memcmp(input.bytes(), &input_before[0],
            input_before.size()) == 0,
        cell, "production encode modified source data or guards");
    CheckOutputGuards(output, output_offset, cell);
    CheckParity(cell, generator, original, recovery);

    std::memset(output.bytes(), kGuardValue, output.size());
    leo2_encode_batch_item item = {
        cell.shard_bytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()
    };
    RequireResult(leo2_encode_batch(codec, &item, 1),
        LEO2_SUCCESS, cell,
        "execute production packed one-item batch terminal");
    RequireCell(std::memcmp(input.bytes(), &input_before[0],
            input_before.size()) == 0,
        cell, "production batch encode modified source data or guards");
    CheckOutputGuards(output, output_offset, cell);
    CheckParity(cell, generator, original, recovery);

    leo2_codec_destroy(codec);
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
        const leo2_result result = leo2_context_create(&options, &context);
        if (result == LEO2_UNSUPPORTED)
        {
            std::printf(
                "production T=4 packed-terminal smoke skipped: AVX2 unavailable\n");
            return 0;
        }
        if (result != LEO2_SUCCESS)
            throw std::runtime_error("create production AVX2 context");

        static const size_t bytes[] = { 64, 128, 256, 512 };
        for (unsigned k = 4; k <= 7; ++k)
        {
            for (unsigned r = 3; r <= 4; ++r)
            {
                for (size_t i = 0; i < sizeof(bytes) / sizeof(bytes[0]); ++i)
                    ExerciseCell(context, Cell{ k, r, bytes[i] });
            }
        }
        static const size_t k8_bytes[] = { 64, 128, 256, 512 };
        for (unsigned r = 3; r <= 4; ++r)
        {
            for (size_t i = 0;
                 i < sizeof(k8_bytes) / sizeof(k8_bytes[0]); ++i)
                ExerciseCell(context, Cell{ 8, r, k8_bytes[i] });
        }
        ExerciseCell(context, Cell{ 4, 3, 1024 });
        ExerciseCell(context, Cell{ 4, 4, 1024 });

        leo2_context_destroy(context);
        std::printf("production T=4 packed-terminal smoke checks passed\n");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr,
            "production T=4 packed-terminal smoke failure: %s\n",
            error.what());
        return 1;
    }
}
