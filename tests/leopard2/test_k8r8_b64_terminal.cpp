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
static const uint8_t kGuardValue = 0xd7;

void Require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

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

void FillGuards(AlignedBuffer& buffer)
{
    std::memset(buffer.bytes(), kGuardValue, buffer.size());
}

std::vector<uint8_t> Snapshot(const AlignedBuffer& buffer)
{
    return std::vector<uint8_t>(
        buffer.bytes(), buffer.bytes() + buffer.size());
}

void RequireSourceUnchanged(
    const AlignedBuffer& input,
    const std::vector<uint8_t>& before,
    const Cell& cell)
{
    RequireCell(before.size() == input.size() &&
            std::memcmp(input.bytes(), &before[0], before.size()) == 0,
        cell, "encode modified a source byte or source guard");
}

void RequireOutputGuards(
    const AlignedBuffer& output,
    size_t payload_offset,
    const Cell& cell)
{
    const size_t payload_bytes =
        static_cast<size_t>(cell.recovery_count) * cell.shard_bytes;
    RequireCell(payload_offset <= output.size() &&
            payload_bytes <= output.size() - payload_offset,
        cell, "test output guard geometry overflowed");
    for (size_t i = 0; i < payload_offset; ++i)
        RequireCell(output.bytes()[i] == kGuardValue, cell,
            "encode modified a leading destination guard");
    for (size_t i = payload_offset + payload_bytes;
         i < output.size(); ++i)
    {
        RequireCell(output.bytes()[i] == kGuardValue, cell,
            "encode modified a trailing destination guard");
    }
}

void FillInput(uint8_t* input, const Cell& cell, uint64_t salt)
{
    uint64_t state = UINT64_C(0x54345041434b4544) ^ salt ^
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

void ResetOutput(uint8_t* output, const Cell& cell)
{
    std::memset(output, 0xa5,
        static_cast<size_t>(cell.recovery_count) * cell.shard_bytes);
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
        if (!recovery[parity])
            continue;
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
                "parity differs from the independent generator oracle");
        }
    }
}

uint64_t T8PackedCalls()
{
    return leopard::ff8::TestOnlyGetHighEncodeCounts().t8_packed_calls;
}

uint64_t T8K8B1024DirectCalls()
{
    return leopard::ff8::TestOnlyGetHighEncodeCounts().
        t8_k8_b1024_direct_calls;
}

uint64_t T8K7B1024DirectCalls()
{
    return leopard::ff8::TestOnlyGetHighEncodeCounts().
        t8_k7_b1024_direct_calls;
}

leo2_codec* CreateCodec(leo2_context* context, const Cell& cell)
{
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context,
        cell.original_count, cell.recovery_count,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        NULL, &codec), LEO2_SUCCESS, cell, "create codec");
    return codec;
}

size_t QueryScratch(leo2_codec* codec, const Cell& cell)
{
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, cell.shard_bytes, &scratch_bytes),
        LEO2_SUCCESS, cell, "query encode scratch");
    const size_t alignment = leo2_scratch_alignment();
    RequireCell(alignment != 0 && (alignment & (alignment - 1U)) == 0,
        cell, "scratch alignment is not a nonzero power of two");
    RequireCell(scratch_bytes != 0, cell,
        "scratch query unexpectedly returned an empty extent");
    return scratch_bytes;
}

size_t ExpectedK8R8ProductionScratch()
{
    const size_t alignment = leo2_scratch_alignment();
    const size_t metadata_bytes =
        16U * 2U * sizeof(uintptr_t) + 24U * sizeof(void*);
    const size_t data_offset =
        (metadata_bytes + alignment - 1U) & ~(alignment - 1U);
    return data_offset + 16U * 64U;
}

void ExerciseCell(
    leo2_context* context,
    const Cell& cell,
    bool expect_terminal)
{
    leo2_codec* codec = CreateCodec(context, cell);
    const size_t scratch_bytes = QueryScratch(codec, cell);
    if (cell.original_count == 8 && cell.recovery_count == 8 &&
        cell.shard_bytes == 64)
    {
        RequireCell(scratch_bytes >= ExpectedK8R8ProductionScratch(), cell,
            "hook scratch query lost the portable production geometry");
    }
    AlignedBuffer scratch(scratch_bytes);
    RequireCell(reinterpret_cast<uintptr_t>(scratch.data()) %
            leo2_scratch_alignment() == 0,
        cell, "exact scratch allocation is misaligned");
    AlignedBuffer input(static_cast<size_t>(cell.original_count) *
        cell.shard_bytes + 2U * kGuardBytes + 8U);
    AlignedBuffer output(static_cast<size_t>(cell.recovery_count) *
        cell.shard_bytes + 2U * kGuardBytes + 8U);
    std::vector<const void*> original(cell.original_count);
    std::vector<void*> recovery(cell.recovery_count);
    const leopard2_test::Matrix generator = MakeGenerator(cell);
    const bool expect_k7_direct = expect_terminal &&
        cell.original_count == 7 && cell.recovery_count == 7 &&
        cell.shard_bytes == 1024;
    const bool expect_k8_direct = expect_terminal &&
        cell.original_count == 8 && cell.recovery_count == 8 &&
        cell.shard_bytes == 1024;

    FillGuards(input);
    FillGuards(output);
    uint8_t* input_base = input.bytes() + kGuardBytes;
    uint8_t* output_base = output.bytes() + kGuardBytes;
    SetPackedPointers(input_base, output_base, cell, original, recovery);
    FillInput(input_base, cell, 0);
    const std::vector<uint8_t> aligned_input_before = Snapshot(input);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, cell.shard_bytes,
        &original[0], &recovery[0], scratch.data(), scratch.size()),
        LEO2_SUCCESS, cell, "execute packed public encode");
    RequireCell(T8PackedCalls() == (expect_terminal ? 1U : 0U), cell,
        "packed public terminal selection mismatch");
    RequireCell(T8K7B1024DirectCalls() ==
            (expect_k7_direct ? 1U : 0U), cell,
        "packed public K7 direct-leaf selection mismatch");
    RequireCell(T8K8B1024DirectCalls() ==
            (expect_k8_direct ? 1U : 0U), cell,
        "packed public direct-leaf selection mismatch");
    CheckParity(cell, generator, original, recovery);
    RequireSourceUnchanged(input, aligned_input_before, cell);
    RequireOutputGuards(output, kGuardBytes, cell);

    FillGuards(output);
    leo2_encode_batch_item item = {
        cell.shard_bytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()
    };
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch(codec, &item, 1),
        LEO2_SUCCESS, cell, "execute packed one-item batch encode");
    RequireCell(T8PackedCalls() == (expect_terminal ? 1U : 0U), cell,
        "packed one-item batch terminal selection mismatch");
    RequireCell(T8K7B1024DirectCalls() ==
            (expect_k7_direct ? 1U : 0U), cell,
        "packed one-item batch K7 direct-leaf selection mismatch");
    RequireCell(T8K8B1024DirectCalls() ==
            (expect_k8_direct ? 1U : 0U), cell,
        "packed one-item batch direct-leaf selection mismatch");
    CheckParity(cell, generator, original, recovery);
    RequireSourceUnchanged(input, aligned_input_before, cell);
    RequireOutputGuards(output, kGuardBytes, cell);

    /* Packed source and destination slabs need not be vector aligned. */
    FillGuards(input);
    FillGuards(output);
    input_base = input.bytes() + kGuardBytes + 1U;
    output_base = output.bytes() + kGuardBytes + 3U;
    SetPackedPointers(input_base, output_base, cell, original, recovery);
    FillInput(input_base, cell, UINT64_C(0x554e414c49474e45));
    const std::vector<uint8_t> unaligned_input_before = Snapshot(input);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, cell.shard_bytes,
        &original[0], &recovery[0], scratch.data(), scratch.size()),
        LEO2_SUCCESS, cell, "execute unaligned packed public encode");
    RequireCell(T8PackedCalls() == (expect_terminal ? 1U : 0U), cell,
        "unaligned packed terminal selection mismatch");
    RequireCell(T8K7B1024DirectCalls() ==
            (expect_k7_direct ? 1U : 0U), cell,
        "unaligned packed K7 direct-leaf selection mismatch");
    RequireCell(T8K8B1024DirectCalls() ==
            (expect_k8_direct ? 1U : 0U), cell,
        "unaligned packed direct-leaf selection mismatch");
    CheckParity(cell, generator, original, recovery);
    RequireSourceUnchanged(input, unaligned_input_before, cell);
    RequireOutputGuards(output, kGuardBytes + 3U, cell);

    leo2_codec_destroy(codec);
}

void ExercisePromotedMatrix(
    leo2_context* context,
    bool full_parity_terminal_available)
{
    static const Cell cells[] = {
        { 5, 5, 256 },
        { 6, 6, 256 },
        { 7, 7, 256 },
        { 8, 8, 256 },
        { 5, 5, 1024 },
        { 6, 6, 1024 },
        { 7, 7, 1024 },
        { 8, 8, 1024 },
        /* Retain the separately qualified balanced 64-byte terminal. */
        { 8, 8, 64 }
    };
    for (size_t i = 0; i < sizeof(cells) / sizeof(cells[0]); ++i)
    {
        const bool is_established_b64 =
            cells[i].original_count == 8 &&
            cells[i].recovery_count == 8 && cells[i].shard_bytes == 64;
        ExerciseCell(context, cells[i],
            is_established_b64 || full_parity_terminal_available);
    }
}

void ExerciseNonPromotedCells(leo2_context* context)
{
    static const size_t boundary_bytes[] = {
        255, 257, 1023, 1025
    };
    for (unsigned count = 5; count <= 8; ++count)
    {
        for (size_t byte_i = 0;
             byte_i < sizeof(boundary_bytes) / sizeof(boundary_bytes[0]);
             ++byte_i)
        {
            ExerciseCell(context,
                Cell{ count, count, boundary_bytes[byte_i] }, false);
        }
    }

    /* Every punctured T=8 neighbor stays on the mature general path. */
    static const size_t exact_bytes[] = { 256, 1024 };
    for (unsigned original_count = 5; original_count <= 8; ++original_count)
        for (unsigned recovery_count = 5; recovery_count <= 8;
             ++recovery_count)
        {
            if (original_count == recovery_count)
                continue;
            for (size_t byte_i = 0;
                 byte_i < sizeof(exact_bytes) / sizeof(exact_bytes[0]);
                 ++byte_i)
            {
                ExerciseCell(context,
                    Cell{ original_count, recovery_count,
                        exact_bytes[byte_i] }, false);
            }
        }

    /* Preserve the established B64 terminal's immediate neighbors. */
    ExerciseCell(context, Cell{ 8, 8, 63 }, false);
    ExerciseCell(context, Cell{ 8, 8, 128 }, false);
    ExerciseCell(context, Cell{ 8, 7, 64 }, false);
    ExerciseCell(context, Cell{ 7, 8, 64 }, false);
}

void ExerciseForcedTransform(leo2_context* context, const Cell& cell)
{
    leo2_codec* codec = CreateCodec(context, cell);
    const size_t scratch_bytes = QueryScratch(codec, cell);
    AlignedBuffer scratch(scratch_bytes);
    AlignedBuffer input(
        static_cast<size_t>(cell.original_count) * cell.shard_bytes);
    AlignedBuffer output(
        static_cast<size_t>(cell.recovery_count) * cell.shard_bytes);
    std::vector<const void*> original(cell.original_count);
    std::vector<void*> recovery(cell.recovery_count);
    SetPackedPointers(
        input.bytes(), output.bytes(), cell, original, recovery);
    FillInput(input.bytes(), cell, UINT64_C(0x464f524345543454));
    ResetOutput(output.bytes(), cell);
    RequireResult(leo2_test_codec_set_encode_mode(
        codec, LEO2_TEST_ENCODE_FORCE_TRANSFORM),
        LEO2_SUCCESS, cell, "force mature transform");
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, cell.shard_bytes,
        &original[0], &recovery[0], scratch.data(), scratch.size()),
        LEO2_SUCCESS, cell, "execute forced mature transform");
    RequireCell(T8PackedCalls() == 0, cell,
        "forced transform entered the packed K=8/R=8/T=8 terminal");
    RequireCell(T8K7B1024DirectCalls() == 0, cell,
        "forced transform entered the K=7/R=7/B=1024 direct leaf");
    RequireCell(T8K8B1024DirectCalls() == 0, cell,
        "forced transform entered the K=8/R=8/B=1024 direct leaf");
    CheckParity(cell, MakeGenerator(cell), original, recovery);
    RequireResult(leo2_test_codec_set_encode_mode(
        codec, LEO2_TEST_ENCODE_AUTO),
        LEO2_SUCCESS, cell, "restore AUTO encode mode");
    leo2_codec_destroy(codec);
}

void ExerciseFallbackLayouts(leo2_context* context, const Cell& cell)
{
    leo2_codec* codec = CreateCodec(context, cell);
    const size_t scratch_bytes = QueryScratch(codec, cell);
    AlignedBuffer scratch(scratch_bytes);
    AlignedBuffer input(
        static_cast<size_t>(cell.original_count) * cell.shard_bytes);
    AlignedBuffer output(
        static_cast<size_t>(cell.recovery_count) * cell.shard_bytes);
    std::vector<const void*> original(cell.original_count);
    std::vector<void*> recovery(cell.recovery_count);
    SetPackedPointers(
        input.bytes(), output.bytes(), cell, original, recovery);
    FillInput(input.bytes(), cell, UINT64_C(0x46414c4c4241434b));
    const leopard2_test::Matrix generator = MakeGenerator(cell);

    /* An all-null output set is a valid no-op through both public forms. */
    ResetOutput(output.bytes(), cell);
    const std::vector<uint8_t> no_output_before(
        output.bytes(), output.bytes() + output.size());
    for (unsigned i = 0; i < cell.recovery_count; ++i)
        recovery[i] = NULL;
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, cell.shard_bytes,
        &original[0], &recovery[0], scratch.data(), scratch.size()),
        LEO2_SUCCESS, cell, "execute all-null public no-op");
    RequireCell(T8PackedCalls() == 0, cell,
        "all-null public output entered the dense terminal");
    RequireCell(std::memcmp(output.bytes(), &no_output_before[0],
            no_output_before.size()) == 0,
        cell, "all-null public encode modified output");

    leo2_encode_batch_item no_output_item = {
        cell.shard_bytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()
    };
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch(codec, &no_output_item, 1),
        LEO2_SUCCESS, cell, "execute all-null one-item batch no-op");
    RequireCell(T8PackedCalls() == 0, cell,
        "all-null batch output entered the dense terminal");
    RequireCell(std::memcmp(output.bytes(), &no_output_before[0],
            no_output_before.size()) == 0,
        cell, "all-null batch encode modified output");

    SetPackedPointers(
        input.bytes(), output.bytes(), cell, original, recovery);
    ResetOutput(output.bytes(), cell);
    void* const omitted_output3 = recovery[3];
    recovery[3] = NULL;
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, cell.shard_bytes,
        &original[0], &recovery[0], scratch.data(), scratch.size()),
        LEO2_SUCCESS, cell, "execute sparse-output fallback");
    RequireCell(T8PackedCalls() == 0, cell,
        "sparse output entered the dense packed terminal");
    CheckParity(cell, generator, original, recovery);
    for (size_t i = 0; i < cell.shard_bytes; ++i)
    {
        RequireCell(static_cast<const uint8_t*>(omitted_output3)[i] == 0xa5,
            cell, "sparse fallback modified null parity shard three");
    }

    SetPackedPointers(
        input.bytes(), output.bytes(), cell, original, recovery);
    std::vector<uint8_t> detached_input(cell.shard_bytes);
    std::vector<uint8_t> detached_output(cell.shard_bytes, 0xa5);
    const unsigned detached_input_index = cell.original_count / 2U;
    const unsigned detached_output_index = cell.recovery_count / 2U;
    std::memcpy(&detached_input[0], original[detached_input_index],
        cell.shard_bytes);
    original[detached_input_index] = &detached_input[0];
    recovery[detached_output_index] = &detached_output[0];
    ResetOutput(output.bytes(), cell);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, cell.shard_bytes,
        &original[0], &recovery[0], scratch.data(), scratch.size()),
        LEO2_SUCCESS, cell, "execute detached-layout fallback");
    RequireCell(T8PackedCalls() == 0, cell,
        "detached layout entered the packed terminal");
    CheckParity(cell, generator, original, recovery);
    for (size_t i = 0; i < cell.shard_bytes; ++i)
        RequireCell(output.bytes()[
                static_cast<size_t>(detached_output_index) *
                    cell.shard_bytes + i] == 0xa5,
            cell, "detached-output fallback wrote the abandoned packed row");

    leo2_codec_destroy(codec);
}

void ExerciseValidationAtomicity(
    leo2_context* context,
    const Cell& cell)
{
    leo2_codec* codec = CreateCodec(context, cell);
    const size_t scratch_bytes = QueryScratch(codec, cell);
    AlignedBuffer scratch(scratch_bytes + leo2_scratch_alignment());
    AlignedBuffer input(
        static_cast<size_t>(cell.original_count) * cell.shard_bytes);
    AlignedBuffer output(
        static_cast<size_t>(cell.recovery_count) * cell.shard_bytes);
    std::vector<const void*> original(cell.original_count);
    std::vector<void*> recovery(cell.recovery_count);
    SetPackedPointers(
        input.bytes(), output.bytes(), cell, original, recovery);
    FillInput(input.bytes(), cell, UINT64_C(0x41544f4d49435434));

    ResetOutput(output.bytes(), cell);
    const std::vector<uint8_t> output_before(
        output.bytes(), output.bytes() +
            static_cast<size_t>(cell.recovery_count) * cell.shard_bytes);

    /* Scratch errors retain precedence even when pointer arrays are null. */
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, cell.shard_bytes,
        NULL, NULL, scratch.data(), scratch_bytes - 1U),
        LEO2_SCRATCH_TOO_SMALL, cell,
        "preserve scratch precedence over null pointer arrays");
    RequireCell(T8PackedCalls() == 0, cell,
        "invalid scratch/null arrays reached the packed terminal");
    RequireCell(std::memcmp(output.bytes(), &output_before[0],
            output_before.size()) == 0,
        cell, "scratch-precedence rejection modified output");

    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, cell.shard_bytes,
        NULL, &recovery[0], scratch.data(), scratch_bytes),
        LEO2_INVALID_ARGUMENT, cell, "reject null original pointer array");
    RequireCell(T8PackedCalls() == 0, cell,
        "null original array reached the packed terminal");
    RequireResult(leo2_encode(codec, cell.shard_bytes,
        &original[0], NULL, scratch.data(), scratch_bytes),
        LEO2_INVALID_ARGUMENT, cell, "reject null recovery pointer array");
    RequireCell(T8PackedCalls() == 0, cell,
        "null recovery array reached the packed terminal");
    RequireCell(std::memcmp(output.bytes(), &output_before[0],
            output_before.size()) == 0,
        cell, "null-array rejection modified output");

    original[4] = NULL;
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, cell.shard_bytes,
        &original[0], &recovery[0], scratch.data(), scratch_bytes),
        LEO2_INVALID_ARGUMENT, cell, "reject null source shard");
    RequireCell(T8PackedCalls() == 0, cell,
        "null source shard reached the packed terminal");
    RequireCell(std::memcmp(output.bytes(), &output_before[0],
            output_before.size()) == 0,
        cell, "null-source rejection modified output");
    SetPackedPointers(
        input.bytes(), output.bytes(), cell, original, recovery);

    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, cell.shard_bytes,
        &original[0], &recovery[0], scratch.data(), scratch_bytes - 1U),
        LEO2_SCRATCH_TOO_SMALL, cell, "reject undersized exact scratch");
    RequireCell(T8PackedCalls() == 0, cell,
        "undersized scratch reached the packed terminal");
    RequireCell(std::memcmp(output.bytes(), &output_before[0],
            output_before.size()) == 0,
        cell, "undersized scratch modified output");

    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, cell.shard_bytes,
        &original[0], &recovery[0], scratch.bytes() + 1, scratch_bytes),
        LEO2_BAD_ALIGNMENT, cell, "reject misaligned exact scratch");
    RequireCell(T8PackedCalls() == 0, cell,
        "misaligned scratch reached the packed terminal");
    RequireCell(std::memcmp(output.bytes(), &output_before[0],
            output_before.size()) == 0,
        cell, "misaligned scratch modified output");

    for (unsigned i = 0; i < cell.recovery_count; ++i)
        recovery[i] = input.bytes() +
            static_cast<size_t>(i) * cell.shard_bytes;
    const std::vector<uint8_t> input_before(
        input.bytes(), input.bytes() +
            static_cast<size_t>(cell.original_count) * cell.shard_bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, cell.shard_bytes,
        &original[0], &recovery[0], scratch.data(), scratch_bytes),
        LEO2_OVERLAP, cell, "reject packed input/output overlap");
    RequireCell(T8PackedCalls() == 0, cell,
        "input/output overlap reached the packed terminal");
    RequireCell(std::memcmp(input.bytes(), &input_before[0],
            input_before.size()) == 0,
        cell, "input/output-overlap rejection modified input");

    SetPackedPointers(
        input.bytes(), output.bytes(), cell, original, recovery);
    ResetOutput(output.bytes(), cell);
    recovery[3] = recovery[2];
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, cell.shard_bytes,
        &original[0], &recovery[0], scratch.data(), scratch_bytes),
        LEO2_OVERLAP, cell, "reject overlapping parity outputs");
    RequireCell(T8PackedCalls() == 0, cell,
        "output/output overlap reached the packed terminal");
    RequireCell(std::memcmp(output.bytes(), &output_before[0],
            output_before.size()) == 0,
        cell, "output/output-overlap rejection modified output");

    SetPackedPointers(
        input.bytes(), output.bytes(), cell, original, recovery);
    AlignedBuffer protected_storage(
        static_cast<size_t>(cell.recovery_count) * cell.shard_bytes);
    leo2_encode_batch_item* const protected_item =
        new (protected_storage.data()) leo2_encode_batch_item;
    protected_item->shard_bytes = cell.shard_bytes;
    protected_item->original = &original[0];
    protected_item->recovery = &recovery[0];
    protected_item->scratch = scratch.data();
    protected_item->scratch_bytes = scratch_bytes;
    for (unsigned i = 0; i < cell.recovery_count; ++i)
        recovery[i] = protected_storage.bytes() +
            static_cast<size_t>(i) * cell.shard_bytes;
    const std::vector<uint8_t> protected_before(
        protected_storage.bytes(),
        protected_storage.bytes() + protected_storage.size());
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch(codec, protected_item, 1),
        LEO2_OVERLAP, cell, "reject output/batch-metadata overlap");
    RequireCell(T8PackedCalls() == 0, cell,
        "batch-metadata overlap reached the packed terminal");
    RequireCell(std::memcmp(protected_storage.bytes(), &protected_before[0],
            protected_before.size()) == 0,
        cell, "batch-metadata-overlap rejection was not atomic");

    leo2_codec_destroy(codec);
}

void ExerciseScalarFallbacks()
{
    leo2_context_options options = {};
    options.struct_size = sizeof(options);
    options.backend = LEO2_BACKEND_SCALAR;
    options.thread_count = 1;
    leo2_context* context = NULL;
    Require(leo2_context_create(&options, &context) == LEO2_SUCCESS,
        "create scalar K=8/R=8/T=8 fallback context");
    static const Cell cells[] = {
        { 5, 5, 256 }, { 6, 6, 256 }, { 7, 7, 256 }, { 8, 8, 256 },
        { 5, 5, 1024 }, { 6, 6, 1024 }, { 7, 7, 1024 },
        { 8, 8, 1024 }, { 8, 8, 64 }
    };
    for (size_t i = 0; i < sizeof(cells) / sizeof(cells[0]); ++i)
        ExerciseCell(context, cells[i], false);
    leo2_context_destroy(context);
}

void ExerciseAdjacentShardSubobjects(leo2_context* context)
{
    /* Each row is a distinct array subobject even though the enclosing array
       gives the terminal the exact numeric 1024-byte stride it requires.
       This guards the direct leaf's integer-address implementation against a
       regression to cross-object C++ pointer arithmetic from row zero. */
    struct ShardObject
    {
        uint8_t bytes[1024];
    };
    static_assert(sizeof(ShardObject) == 1024,
        "direct-leaf subobject fixture gained padding");

    for (unsigned count = 7; count <= 8; ++count)
    {
        const Cell cell = { count, count, 1024 };
        leo2_codec* codec = CreateCodec(context, cell);
        AlignedBuffer scratch(QueryScratch(codec, cell));
        ShardObject input[8];
        ShardObject input_before[8];
        ShardObject output[8];
        std::vector<const void*> original(count);
        std::vector<void*> recovery(count);
        for (unsigned lane = 0; lane < 8; ++lane)
        {
            for (size_t i = 0; i < sizeof(input[lane].bytes); ++i)
            {
                input[lane].bytes[i] = static_cast<uint8_t>(
                    lane * 71U + i * 43U + (i >> 3) * 19U + 0x5dU);
                output[lane].bytes[i] = 0xa5;
            }
            if (lane < count)
            {
                original[lane] = input[lane].bytes;
                recovery[lane] = output[lane].bytes;
            }
        }
        std::memcpy(input_before, input, sizeof(input));

        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode(codec, cell.shard_bytes,
            &original[0], &recovery[0], scratch.data(), scratch.size()),
            LEO2_SUCCESS, cell, "encode adjacent shard subobjects");
        RequireCell(T8PackedCalls() == 1, cell,
            "adjacent shard subobjects missed the packed terminal");
        RequireCell(T8K7B1024DirectCalls() == (count == 7 ? 1U : 0U), cell,
            "adjacent shard subobjects selected the wrong K7 direct leaf");
        RequireCell(T8K8B1024DirectCalls() == (count == 8 ? 1U : 0U), cell,
            "adjacent shard subobjects selected the wrong K8 direct leaf");
        CheckParity(cell, MakeGenerator(cell), original, recovery);
        RequireCell(std::memcmp(input, input_before, sizeof(input)) == 0, cell,
            "adjacent shard subobject encode modified source data");
        leo2_codec_destroy(codec);
    }
}

void ExerciseDirectLinearBasis(leo2_context* context)
{
    /* The direct leaves are GF(2)-linear.  Exercise every source/bit basis
       vector against the independent generator oracle so the handwritten
       circuit is not accepted merely because a dense fixture happened to
       mask a transcribed coefficient. */
    struct ShardObject
    {
        uint8_t bytes[1024];
    };
    static_assert(sizeof(ShardObject) == 1024,
        "direct-leaf basis fixture gained padding");

    for (unsigned count = 7; count <= 8; ++count)
    {
        const Cell cell = { count, count, 1024 };
        leo2_codec* codec = CreateCodec(context, cell);
        AlignedBuffer scratch(QueryScratch(codec, cell));
        ShardObject input[8];
        ShardObject input_before[8];
        ShardObject output[8];
        std::vector<const void*> original(count);
        std::vector<void*> recovery(count);
        for (unsigned lane = 0; lane < count; ++lane)
        {
            original[lane] = input[lane].bytes;
            recovery[lane] = output[lane].bytes;
        }
        const leopard2_test::Matrix generator = MakeGenerator(cell);

        for (unsigned source = 0; source < count; ++source)
        {
            for (unsigned bit = 0; bit < 8; ++bit)
            {
                std::memset(input, 0, sizeof(input));
                std::memset(output, 0xa5, sizeof(output));
                const size_t offset =
                    (source * 131U + bit * 127U) % cell.shard_bytes;
                input[source].bytes[offset] =
                    static_cast<uint8_t>(1U << bit);
                std::memcpy(input_before, input, sizeof(input));

                leopard::ff8::TestOnlyResetHighEncodeCounts();
                RequireResult(leo2_encode(codec, cell.shard_bytes,
                    &original[0], &recovery[0],
                    scratch.data(), scratch.size()),
                    LEO2_SUCCESS, cell, "encode direct linear basis");
                RequireCell(T8PackedCalls() == 1, cell,
                    "linear basis missed the packed terminal");
                RequireCell(T8K7B1024DirectCalls() ==
                        (count == 7 ? 1U : 0U), cell,
                    "linear basis selected the wrong K7 direct leaf");
                RequireCell(T8K8B1024DirectCalls() ==
                        (count == 8 ? 1U : 0U), cell,
                    "linear basis selected the wrong K8 direct leaf");
                CheckParity(cell, generator, original, recovery);
                RequireCell(std::memcmp(
                        input, input_before, sizeof(input)) == 0, cell,
                    "linear-basis encode modified source data");
            }
        }
        leo2_codec_destroy(codec);
    }
}

void ExerciseAutoBackend(bool full_parity_terminal_available)
{
    leo2_context_options options = {};
    options.struct_size = sizeof(options);
    options.backend = LEO2_BACKEND_AUTO;
    options.thread_count = 1;
    leo2_context* context = NULL;
    Require(leo2_context_create(&options, &context) == LEO2_SUCCESS,
        "create AUTO K=8/R=8/T=8 context");
    const bool expect_terminal =
        leo2_context_backend(context) == LEO2_BACKEND_AVX2;
    static const Cell cells[] = {
        { 5, 5, 256 }, { 6, 6, 256 }, { 7, 7, 256 }, { 8, 8, 256 },
        { 5, 5, 1024 }, { 6, 6, 1024 }, { 7, 7, 1024 },
        { 8, 8, 1024 }, { 8, 8, 64 }
    };
    for (size_t i = 0; i < sizeof(cells) / sizeof(cells[0]); ++i)
    {
        const bool is_established_b64 =
            cells[i].original_count == 8 &&
            cells[i].recovery_count == 8 && cells[i].shard_bytes == 64;
        ExerciseCell(context, cells[i], expect_terminal &&
            (is_established_b64 || full_parity_terminal_available));
    }
    leo2_context_destroy(context);
}

} // namespace

int main()
{
    try
    {
        /* Scalar correctness remains useful on hosts without AVX2. */
        const bool full_parity_terminal_available =
            leopard2_internal::
                SetT8FullParityTerminalEnabledForDiagnostics(true);
        ExerciseScalarFallbacks();
        ExerciseAutoBackend(full_parity_terminal_available);

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
                "K=8/R=8/T=8 packed terminal AVX2 checks skipped; scalar checks passed\n");
            return 0;
        }
        Require(context_result == LEO2_SUCCESS,
            "create AVX2 K=8/R=8/T=8 terminal context");

        ExercisePromotedMatrix(context, full_parity_terminal_available);
        if (full_parity_terminal_available)
        {
            ExerciseAdjacentShardSubobjects(context);
            ExerciseDirectLinearBasis(context);
        }
        ExerciseNonPromotedCells(context);
        ExerciseForcedTransform(context, Cell{ 8, 8, 64 });
        ExerciseForcedTransform(context, Cell{ 7, 7, 1024 });
        ExerciseForcedTransform(context, Cell{ 8, 8, 1024 });
        static const Cell promoted_cells[] = {
            { 5, 5, 256 }, { 6, 6, 256 }, { 7, 7, 256 },
            { 8, 8, 256 }, { 5, 5, 1024 }, { 6, 6, 1024 },
            { 7, 7, 1024 }, { 8, 8, 1024 }, { 8, 8, 64 }
        };
        for (size_t i = 0;
             i < sizeof(promoted_cells) / sizeof(promoted_cells[0]); ++i)
        {
            ExerciseFallbackLayouts(context, promoted_cells[i]);
            ExerciseValidationAtomicity(context, promoted_cells[i]);
        }

        leo2_context_destroy(context);
        std::printf("K=8/R=8/T=8 packed terminal checks passed\n");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "K=8/R=8/T=8 packed terminal failure: %s\n",
            error.what());
        return 1;
    }
}
