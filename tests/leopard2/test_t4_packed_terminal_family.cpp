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

void CheckLegacyParity(
    const Cell& cell,
    const std::vector<const void*>& original,
    const std::vector<void*>& recovery)
{
    if ((cell.shard_bytes & 63U) != 0)
        return;

    const unsigned work_count = leo_encode_work_count(
        cell.original_count, cell.recovery_count);
    RequireCell(work_count >= cell.recovery_count, cell,
        "legacy work-count query failed");
    std::vector<std::vector<uint8_t> > work(
        work_count, std::vector<uint8_t>(cell.shard_bytes, 0));
    std::vector<void*> work_pointers(work_count);
    for (unsigned i = 0; i < work_count; ++i)
        work_pointers[i] = &work[i][0];
    RequireCell(leo_encode(cell.shard_bytes,
            cell.original_count, cell.recovery_count, work_count,
            &original[0], &work_pointers[0]) == Leopard_Success,
        cell, "legacy Leopard1 encode failed");
    for (unsigned parity = 0; parity < cell.recovery_count; ++parity)
    {
        if (!recovery[parity])
            continue;
        RequireCell(std::memcmp(recovery[parity], &work[parity][0],
                cell.shard_bytes) == 0,
            cell, "parity differs from legacy Leopard1");
    }
}

uint64_t T4PackedCalls()
{
    return leopard::ff8::TestOnlyGetHighEncodeCounts().t4_packed_calls;
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

void ExerciseCell(
    leo2_context* context,
    const Cell& cell,
    bool expect_terminal)
{
    leo2_codec* codec = CreateCodec(context, cell);
    if (cell.original_count == 8 && cell.shard_bytes == 512)
    {
        RequireCell(leopard2_internal::HighT4BatchMaximumBytes(
                cell.original_count, cell.recovery_count) == 0,
            cell, "K8 unexpectedly entered the reusable T=4 binding policy");
        leopard2_internal::CodecEncodePathInfo path = {};
        RequireCell(leopard2_internal::GetCodecEncodePathInfo(
                codec, cell.shard_bytes, cell.recovery_count, &path),
            cell, "query K8 encode-path metadata");
        RequireCell(!path.high_t4_batch_binding_selected, cell,
            "K8 ordinary terminal leaked into reusable T=4 binding");
    }
    const size_t scratch_bytes = QueryScratch(codec, cell);
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
    RequireCell(T4PackedCalls() == (expect_terminal ? 1U : 0U), cell,
        "packed public terminal selection mismatch");
    CheckParity(cell, generator, original, recovery);
    CheckLegacyParity(cell, original, recovery);
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
    RequireCell(T4PackedCalls() == (expect_terminal ? 1U : 0U), cell,
        "packed one-item batch terminal selection mismatch");
    CheckParity(cell, generator, original, recovery);
    CheckLegacyParity(cell, original, recovery);
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
    RequireCell(T4PackedCalls() == (expect_terminal ? 1U : 0U), cell,
        "unaligned packed terminal selection mismatch");
    CheckParity(cell, generator, original, recovery);
    CheckLegacyParity(cell, original, recovery);
    RequireSourceUnchanged(input, unaligned_input_before, cell);
    RequireOutputGuards(output, kGuardBytes + 3U, cell);

    FillGuards(output);
    item.shard_bytes = cell.shard_bytes;
    item.original = &original[0];
    item.recovery = &recovery[0];
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch(codec, &item, 1),
        LEO2_SUCCESS, cell, "execute unaligned packed one-item batch encode");
    RequireCell(T4PackedCalls() == (expect_terminal ? 1U : 0U), cell,
        "unaligned one-item batch terminal selection mismatch");
    CheckParity(cell, generator, original, recovery);
    CheckLegacyParity(cell, original, recovery);
    RequireSourceUnchanged(input, unaligned_input_before, cell);
    RequireOutputGuards(output, kGuardBytes + 3U, cell);

    leo2_codec_destroy(codec);
}

void ExercisePromotedMatrix(leo2_context* context)
{
    static const size_t bytes[] = { 64, 128, 256, 512 };
    for (unsigned k = 4; k <= 7; ++k)
    {
        for (unsigned r = 3; r <= 4; ++r)
        {
            for (size_t i = 0; i < sizeof(bytes) / sizeof(bytes[0]); ++i)
            {
                const Cell cell = { k, r, bytes[i] };
                ExerciseCell(context, cell, true);
            }
        }
    }
    static const size_t k8_bytes[] = { 64, 128, 256, 512 };
    for (unsigned r = 3; r <= 4; ++r)
    {
        for (size_t i = 0;
             i < sizeof(k8_bytes) / sizeof(k8_bytes[0]); ++i)
            ExerciseCell(context, Cell{ 8, r, k8_bytes[i] }, true);
    }
    ExerciseCell(context, Cell{ 4, 3, 1024 }, true);
    ExerciseCell(context, Cell{ 4, 4, 1024 }, true);
}

void ExerciseNonPromotedCells(leo2_context* context)
{
    static const Cell cells[] = {
        { 4, 3, 63 },
        { 5, 3, 1024 },
        { 7, 4, 1024 },
        { 8, 3, 63 },
        { 8, 4, 63 },
        { 8, 3, 65 },
        { 8, 4, 65 },
        { 8, 3, 127 },
        { 8, 4, 127 },
        { 8, 3, 129 },
        { 8, 4, 129 },
        { 8, 3, 255 },
        { 8, 4, 255 },
        { 8, 3, 257 },
        { 8, 4, 257 },
        { 8, 3, 1024 },
        { 8, 4, 1024 },
        { 8, 8, 64 },
        { 9, 4, 64 }
    };
    for (size_t i = 0; i < sizeof(cells) / sizeof(cells[0]); ++i)
        ExerciseCell(context, cells[i], false);
    for (unsigned k = 4; k <= 8; ++k)
    {
        for (unsigned r = 3; r <= 4; ++r)
        {
            ExerciseCell(context, Cell{ k, r, 511 }, false);
            ExerciseCell(context, Cell{ k, r, 513 }, false);
        }
    }
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
    RequireCell(T4PackedCalls() == 0, cell,
        "forced transform entered the packed T=4 terminal");
    CheckParity(cell, MakeGenerator(cell), original, recovery);
    CheckLegacyParity(cell, original, recovery);
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
    RequireCell(T4PackedCalls() == 0, cell,
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
    RequireCell(T4PackedCalls() == 0, cell,
        "all-null batch output entered the dense terminal");
    RequireCell(std::memcmp(output.bytes(), &no_output_before[0],
            no_output_before.size()) == 0,
        cell, "all-null batch encode modified output");

    SetPackedPointers(
        input.bytes(), output.bytes(), cell, original, recovery);
    ResetOutput(output.bytes(), cell);
    const unsigned omitted_index = cell.recovery_count / 2U;
    void* const omitted_output = recovery[omitted_index];
    recovery[omitted_index] = NULL;
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, cell.shard_bytes,
        &original[0], &recovery[0], scratch.data(), scratch.size()),
        LEO2_SUCCESS, cell, "execute sparse-output fallback");
    RequireCell(T4PackedCalls() == 0, cell,
        "sparse output entered the dense packed terminal");
    CheckParity(cell, generator, original, recovery);
    CheckLegacyParity(cell, original, recovery);
    for (size_t i = 0; i < cell.shard_bytes; ++i)
        RequireCell(static_cast<const uint8_t*>(omitted_output)[i] == 0xa5,
            cell, "sparse fallback modified a null parity shard");

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
    RequireCell(T4PackedCalls() == 0, cell,
        "detached layout entered the packed terminal");
    CheckParity(cell, generator, original, recovery);
    CheckLegacyParity(cell, original, recovery);
    for (size_t i = 0; i < cell.shard_bytes; ++i)
        RequireCell(output.bytes()[
                static_cast<size_t>(detached_output_index) *
                    cell.shard_bytes + i] == 0xa5,
            cell, "detached-output fallback wrote the abandoned packed row");

    /* Input aliases are part of the public contract, but the dense terminal
       requires distinct packed rows.  Exercise exact and partial aliases. */
    const unsigned aliased_input_index = cell.original_count / 2U;
    for (unsigned alias_kind = 0; alias_kind < 2; ++alias_kind)
    {
        SetPackedPointers(
            input.bytes(), output.bytes(), cell, original, recovery);
        original[aliased_input_index] = alias_kind == 0
            ? original[0] : input.bytes() + 1U;
        const std::vector<uint8_t> aliased_input_before = Snapshot(input);

        ResetOutput(output.bytes(), cell);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode(codec, cell.shard_bytes,
            &original[0], &recovery[0], scratch.data(), scratch.size()),
            LEO2_SUCCESS, cell, "execute input-alias public fallback");
        RequireCell(T4PackedCalls() == 0, cell,
            "input alias entered the dense packed terminal");
        CheckParity(cell, generator, original, recovery);
        CheckLegacyParity(cell, original, recovery);
        RequireSourceUnchanged(input, aliased_input_before, cell);

        ResetOutput(output.bytes(), cell);
        leo2_encode_batch_item aliased_item = {
            cell.shard_bytes, &original[0], &recovery[0],
            scratch.data(), scratch.size()
        };
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode_batch(codec, &aliased_item, 1),
            LEO2_SUCCESS, cell,
            "execute input-alias one-item batch fallback");
        RequireCell(T4PackedCalls() == 0, cell,
            "batched input alias entered the dense packed terminal");
        CheckParity(cell, generator, original, recovery);
        CheckLegacyParity(cell, original, recovery);
        RequireSourceUnchanged(input, aliased_input_before, cell);
    }

    leo2_codec_destroy(codec);
}

void ExerciseValidationAtomicity(leo2_context* context, const Cell& cell)
{
    RequireCell(cell.recovery_count == 4 &&
            (cell.shard_bytes == 64 || cell.shard_bytes == 512),
        cell, "invalid validation-atomicity fixture");
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
    RequireCell(T4PackedCalls() == 0, cell,
        "invalid scratch/null arrays reached the packed terminal");
    RequireCell(std::memcmp(output.bytes(), &output_before[0],
            output_before.size()) == 0,
        cell, "scratch-precedence rejection modified output");

    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, cell.shard_bytes,
        NULL, &recovery[0], scratch.data(), scratch_bytes),
        LEO2_INVALID_ARGUMENT, cell, "reject null original pointer array");
    RequireCell(T4PackedCalls() == 0, cell,
        "null original array reached the packed terminal");
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, cell.shard_bytes,
        &original[0], NULL, scratch.data(), scratch_bytes),
        LEO2_INVALID_ARGUMENT, cell, "reject null recovery pointer array");
    RequireCell(T4PackedCalls() == 0, cell,
        "null recovery array reached the packed terminal");
    RequireCell(std::memcmp(output.bytes(), &output_before[0],
            output_before.size()) == 0,
        cell, "null-array rejection modified output");

    original[cell.original_count / 2U] = NULL;
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, cell.shard_bytes,
        &original[0], &recovery[0], scratch.data(), scratch_bytes),
        LEO2_INVALID_ARGUMENT, cell, "reject null source shard");
    RequireCell(T4PackedCalls() == 0, cell,
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
    RequireCell(T4PackedCalls() == 0, cell,
        "undersized scratch reached the packed terminal");
    RequireCell(std::memcmp(output.bytes(), &output_before[0],
            output_before.size()) == 0,
        cell, "undersized scratch modified output");

    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, cell.shard_bytes,
        &original[0], &recovery[0], scratch.bytes() + 1, scratch_bytes),
        LEO2_BAD_ALIGNMENT, cell, "reject misaligned exact scratch");
    RequireCell(T4PackedCalls() == 0, cell,
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
    RequireCell(T4PackedCalls() == 0, cell,
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
    RequireCell(T4PackedCalls() == 0, cell,
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
    RequireCell(T4PackedCalls() == 0, cell,
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
        "create scalar T=4 fallback context");
    static const Cell cells[] = {
        { 4, 3, 64 },
        { 5, 4, 128 },
        { 7, 3, 256 },
        { 6, 4, 512 },
        { 4, 4, 1024 },
        { 8, 3, 128 },
        { 8, 4, 256 },
        { 8, 4, 512 }
    };
    for (size_t i = 0; i < sizeof(cells) / sizeof(cells[0]); ++i)
        ExerciseCell(context, cells[i], false);
    leo2_context_destroy(context);
}

} // namespace

int main()
{
    try
    {
        Require(leo_init() == Leopard_Success,
            "initialize legacy Leopard parity oracle");

        /* Scalar correctness remains useful on hosts without AVX2. */
        ExerciseScalarFallbacks();

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
                "T=4 packed terminal AVX2 checks skipped; scalar checks passed\n");
            return 0;
        }
        Require(context_result == LEO2_SUCCESS,
            "create AVX2 T=4 terminal context");

        ExercisePromotedMatrix(context);
        ExerciseNonPromotedCells(context);
        ExerciseForcedTransform(context, Cell{ 7, 4, 256 });
        ExerciseForcedTransform(context, Cell{ 7, 4, 512 });
        ExerciseForcedTransform(context, Cell{ 8, 3, 128 });
        ExerciseForcedTransform(context, Cell{ 8, 3, 512 });
        ExerciseFallbackLayouts(context, Cell{ 4, 3, 64 });
        ExerciseFallbackLayouts(context, Cell{ 7, 4, 256 });
        ExerciseFallbackLayouts(context, Cell{ 7, 4, 512 });
        ExerciseFallbackLayouts(context, Cell{ 8, 3, 64 });
        ExerciseFallbackLayouts(context, Cell{ 8, 4, 256 });
        ExerciseFallbackLayouts(context, Cell{ 8, 4, 512 });
        ExerciseValidationAtomicity(context, Cell{ 4, 4, 64 });
        ExerciseValidationAtomicity(context, Cell{ 7, 4, 512 });
        ExerciseValidationAtomicity(context, Cell{ 8, 4, 64 });
        ExerciseValidationAtomicity(context, Cell{ 8, 4, 512 });

        leo2_context_destroy(context);
        std::printf("T=4 packed terminal family checks passed\n");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "T=4 packed terminal family failure: %s\n",
            error.what());
        return 1;
    }
}
