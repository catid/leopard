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
#include "LeopardFF8.h"
#include "leopard.h"
#include "leopard2.h"

#include <algorithm>
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

static const unsigned kSide = 16;
static const size_t kShardBytes = 64;
static const size_t kPayloadBytes = kSide * kShardBytes;
static const size_t kGuardBytes = 64;
static const uint8_t kGuardValue = 0x6d;

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

void FillRandom(uint8_t* bytes, size_t count, uint64_t seed)
{
    uint64_t state = seed;
    for (size_t i = 0; i < count; ++i)
    {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        bytes[i] = static_cast<uint8_t>(state >> 24);
    }
}

void CheckParity(
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator,
    unsigned original_count,
    unsigned recovery_count,
    const std::vector<const void*>& original,
    const std::vector<void*>& recovery,
    size_t shard_bytes,
    const char* message)
{
    for (unsigned parity = 0; parity < recovery_count; ++parity)
    {
        if (!recovery[parity])
            continue;
        const std::vector<leopard2_test::Element>& row =
            generator[original_count + parity];
        const uint8_t* const output =
            static_cast<const uint8_t*>(recovery[parity]);
        for (size_t offset = 0; offset < shard_bytes; ++offset)
        {
            leopard2_test::Element expected = 0;
            for (unsigned source = 0; source < original_count; ++source)
            {
                const uint8_t value = static_cast<const uint8_t*>(
                    original[source])[offset];
                expected = field.add(expected,
                    field.multiply(row[source], value));
            }
            Require(output[offset] == expected, message);
        }
    }
}

leo2_codec* CreateCodec(
    leo2_context* context,
    unsigned original_count,
    unsigned recovery_count,
    leo2_profile profile)
{
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context,
        original_count, recovery_count, profile, LEO2_FIELD_GF8,
        NULL, &codec), LEO2_SUCCESS, "create T16 qualification codec");
    return codec;
}

void RequireRouteCounts(
    uint64_t three_pass_calls,
    uint64_t packed_calls,
    const char* message)
{
    const leopard::ff8::TestOnlyHighEncodeCounts counts =
        leopard::ff8::TestOnlyGetHighEncodeCounts();
    if (counts.t16_b64_three_pass_calls != three_pass_calls ||
        counts.balanced_b64_packed_calls != packed_calls)
        throw std::runtime_error(message);
}

void EncodePublic(
    leo2_codec* codec,
    const std::vector<const void*>& original,
    std::vector<void*>& recovery,
    AlignedBuffer& scratch,
    bool batch,
    uint64_t expected_three_pass_calls,
    uint64_t expected_packed_calls,
    const char* message)
{
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    leo2_result result;
    if (batch)
    {
        leo2_encode_batch_item item = {
            kShardBytes, &original[0], &recovery[0],
            scratch.data(), scratch.size()
        };
        result = leo2_encode_batch(codec, &item, 1);
    }
    else
    {
        result = leo2_encode(codec, kShardBytes,
            &original[0], &recovery[0], scratch.data(), scratch.size());
    }
    RequireResult(result, LEO2_SUCCESS, message);
    RequireRouteCounts(
        expected_three_pass_calls, expected_packed_calls, message);
}

void EncodeLegacyGeneric(
    const std::vector<const void*>& original,
    AlignedBuffer& legacy_storage,
    std::vector<void*>& legacy_work)
{
    const unsigned work_count = leo_encode_work_count(kSide, kSide);
    Require(work_count == legacy_work.size(),
        "legacy T16 work geometry changed during qualification");
    for (unsigned i = 0; i < work_count; ++i)
        legacy_work[i] = legacy_storage.bytes() + i * kShardBytes;
    std::memset(legacy_storage.bytes(), 0xa5, legacy_storage.size());
    Require(leopard::ff8::SetHighT16B64ThreePassEnabledForDiagnostics(false),
        "disable T16 candidate for legacy oracle");
    Require(leo_encode(kShardBytes, kSide, kSide, work_count,
        &original[0], &legacy_work[0]) == Leopard_Success,
        "legacy generic T16 encode failed");
}

void ExerciseBasisAndRandom(
    leo2_context* context,
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator)
{
    leo2_codec* codec = CreateCodec(
        context, kSide, kSide, LEO2_PROFILE_LEGACY_HIGH_V1);
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, kShardBytes, &scratch_bytes), LEO2_SUCCESS,
        "query T16 qualification scratch");
    AlignedBuffer scratch(scratch_bytes);
    AlignedBuffer input(kPayloadBytes);
    AlignedBuffer candidate(kPayloadBytes);
    AlignedBuffer generic(kPayloadBytes);
    const unsigned legacy_work_count = leo_encode_work_count(kSide, kSide);
    Require(legacy_work_count != 0, "query legacy T16 work count");
    AlignedBuffer legacy_storage(
        static_cast<size_t>(legacy_work_count) * kShardBytes);
    std::vector<void*> legacy_work(legacy_work_count);
    std::vector<const void*> original;
    std::vector<void*> candidate_recovery;
    std::vector<void*> generic_recovery;
    SetPackedPointers(input.bytes(), candidate.bytes(), kSide, kSide,
        kShardBytes, original, candidate_recovery);
    SetPackedPointers(input.bytes(), generic.bytes(), kSide, kSide,
        kShardBytes, original, generic_recovery);

    for (unsigned basis = 0; basis < kSide; ++basis)
    {
        std::memset(input.bytes(), 0, input.size());
        std::memset(input.bytes() + basis * kShardBytes, 1, kShardBytes);
        Require(leopard::ff8::SetHighT16B64ThreePassEnabledForDiagnostics(true),
            "enable T16 basis candidate");
        EncodePublic(codec, original, candidate_recovery, scratch,
            false, 1, 1, "encode T16 basis candidate");
        CheckParity(field, generator, kSide, kSide, original,
            candidate_recovery, kShardBytes,
            "T16 basis candidate differs from direct generator");

        Require(leopard::ff8::SetHighT16B64ThreePassEnabledForDiagnostics(false),
            "disable T16 basis candidate");
        EncodePublic(codec, original, generic_recovery, scratch,
            false, 0, 0, "encode T16 basis generic control");
        Require(std::memcmp(candidate.bytes(), generic.bytes(),
                kPayloadBytes) == 0,
            "T16 basis candidate differs from generic transform");

        EncodeLegacyGeneric(original, legacy_storage, legacy_work);
        Require(std::memcmp(candidate.bytes(), legacy_storage.bytes(),
                kPayloadBytes) == 0,
            "T16 basis candidate differs from legacy parity");
    }

    for (unsigned trial = 0; trial < 1024; ++trial)
    {
        FillRandom(input.bytes(), input.size(),
            UINT64_C(0x7431366236340000) ^ trial);
        Require(leopard::ff8::SetHighT16B64ThreePassEnabledForDiagnostics(true),
            "enable T16 random candidate");
        EncodePublic(codec, original, candidate_recovery, scratch,
            false, 1, 1, "encode T16 deterministic-random candidate");
        CheckParity(field, generator, kSide, kSide, original,
            candidate_recovery, kShardBytes,
            "T16 random candidate differs from direct generator");

        Require(leopard::ff8::SetHighT16B64ThreePassEnabledForDiagnostics(false),
            "disable T16 random candidate");
        EncodePublic(codec, original, generic_recovery, scratch,
            false, 0, 0, "encode T16 deterministic-random generic control");
        Require(std::memcmp(candidate.bytes(), generic.bytes(),
                kPayloadBytes) == 0,
            "T16 random candidate differs from generic transform");

        EncodeLegacyGeneric(original, legacy_storage, legacy_work);
        Require(std::memcmp(candidate.bytes(), legacy_storage.bytes(),
                kPayloadBytes) == 0,
            "T16 random candidate differs from legacy parity");
    }
    leo2_codec_destroy(codec);
}

void CheckGuards(
    const std::vector<uint8_t>& storage,
    size_t payload_offset,
    size_t payload_bytes,
    const char* message)
{
    for (size_t i = 0; i < payload_offset; ++i)
        Require(storage[i] == kGuardValue, message);
    for (size_t i = payload_offset + payload_bytes;
         i < storage.size(); ++i)
        Require(storage[i] == kGuardValue, message);
}

void ExerciseOffsetsAndBatch(
    leo2_context* context,
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator)
{
    leo2_codec* codec = CreateCodec(
        context, kSide, kSide, LEO2_PROFILE_LEGACY_HIGH_V1);
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, kShardBytes, &scratch_bytes), LEO2_SUCCESS,
        "query T16 offset scratch");
    AlignedBuffer scratch(scratch_bytes);

    for (size_t offset = 0; offset < 32; ++offset)
    {
        const size_t input_offset = kGuardBytes + offset;
        const size_t output_offset = kGuardBytes + (31U - offset);
        std::vector<uint8_t> input(
            input_offset + kPayloadBytes + kGuardBytes, kGuardValue);
        std::vector<uint8_t> output(
            output_offset + kPayloadBytes + kGuardBytes, kGuardValue);
        FillRandom(&input[input_offset], kPayloadBytes,
            UINT64_C(0x6f66667365740000) ^ offset);
        const std::vector<uint8_t> input_before = input;
        std::vector<const void*> original;
        std::vector<void*> recovery;
        SetPackedPointers(&input[input_offset], &output[output_offset],
            kSide, kSide, kShardBytes, original, recovery);

        Require(leopard::ff8::SetHighT16B64ThreePassEnabledForDiagnostics(true),
            "enable T16 offset candidate");
        EncodePublic(codec, original, recovery, scratch,
            false, 1, 1, "encode T16 offset candidate");
        Require(input == input_before,
            "T16 offset encode modified source data or guards");
        CheckGuards(output, output_offset, kPayloadBytes,
            "T16 offset encode modified output guards");
        CheckParity(field, generator, kSide, kSide, original, recovery,
            kShardBytes, "T16 offset parity differs from direct generator");

        std::fill(output.begin(), output.end(), kGuardValue);
        EncodePublic(codec, original, recovery, scratch,
            true, 1, 1, "encode T16 offset one-item batch");
        Require(input == input_before,
            "T16 offset batch modified source data or guards");
        CheckGuards(output, output_offset, kPayloadBytes,
            "T16 offset batch modified output guards");
        CheckParity(field, generator, kSide, kSide, original, recovery,
            kShardBytes,
            "T16 offset batch parity differs from direct generator");
    }
    leo2_codec_destroy(codec);
}

void ExerciseShape(
    leo2_context* context,
    unsigned original_count,
    unsigned recovery_count,
    leo2_profile profile,
    size_t shard_bytes,
    uint64_t expected_three_pass_calls,
    uint64_t expected_packed_calls,
    const char* message)
{
    leo2_codec* codec = CreateCodec(
        context, original_count, recovery_count, profile);
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, shard_bytes, &scratch_bytes), LEO2_SUCCESS,
        "query T16 neighbor scratch");
    AlignedBuffer scratch(scratch_bytes);
    AlignedBuffer input(
        static_cast<size_t>(original_count) * shard_bytes);
    AlignedBuffer output(
        static_cast<size_t>(recovery_count) * shard_bytes);
    FillRandom(input.bytes(), input.size(),
        UINT64_C(0x6e65696768626f72) ^ original_count ^
            (static_cast<uint64_t>(recovery_count) << 16) ^ shard_bytes);
    std::vector<const void*> original;
    std::vector<void*> recovery;
    SetPackedPointers(input.bytes(), output.bytes(), original_count,
        recovery_count, shard_bytes, original, recovery);
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            profile == LEO2_PROFILE_LOW_V1
                ? leopard2_test::kLow : leopard2_test::kLegacyHigh,
            original_count, recovery_count);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);

    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()), LEO2_SUCCESS, message);
    RequireRouteCounts(
        expected_three_pass_calls, expected_packed_calls, message);
    CheckParity(field, generator, original_count, recovery_count,
        original, recovery, shard_bytes, message);
    leo2_codec_destroy(codec);
}

void ExerciseNeighborRoutes(leo2_context* avx2_context)
{
    Require(leopard::ff8::SetHighT16B64ThreePassEnabledForDiagnostics(true),
        "enable T16 neighbor candidate");
    ExerciseShape(avx2_context, 16, 16, LEO2_PROFILE_LEGACY_HIGH_V1,
        63, 0, 0, "B63 entered a T16/B64 candidate route");
    ExerciseShape(avx2_context, 16, 16, LEO2_PROFILE_LEGACY_HIGH_V1,
        65, 0, 0, "B65 entered a T16/B64 candidate route");
    ExerciseShape(avx2_context, 16, 16, LEO2_PROFILE_LEGACY_HIGH_V1,
        128, 0, 0, "B128 entered a T16/B64 candidate route");
    ExerciseShape(avx2_context, 15, 16, LEO2_PROFILE_LEGACY_HIGH_V1,
        64, 0, 0, "K15/R16 entered a T16/B64 candidate route");
    ExerciseShape(avx2_context, 16, 15, LEO2_PROFILE_LEGACY_HIGH_V1,
        64, 0, 0, "K16/R15 entered a T16/B64 candidate route");
    ExerciseShape(avx2_context, 8, 8, LEO2_PROFILE_LEGACY_HIGH_V1,
        64, 0, 0, "T8 entered the T16/B64 candidate");
    ExerciseShape(avx2_context, 32, 32, LEO2_PROFILE_LEGACY_HIGH_V1,
        64, 0, 1, "T32 disagreed with the promoted balanced wrapper route");
    ExerciseShape(avx2_context, 16, 16, LEO2_PROFILE_LOW_V1,
        64, 0, 0, "low profile entered a T16/B64 candidate route");
}

void ExerciseBackendRoute(leo2_backend requested_backend)
{
    leo2_context_options options = {};
    options.struct_size = sizeof(options);
    options.backend = requested_backend;
    options.thread_count = 1;
    leo2_context* context = NULL;
    const leo2_result context_result = leo2_context_create(&options, &context);
    if (context_result == LEO2_UNSUPPORTED)
        return;
    RequireResult(context_result, LEO2_SUCCESS,
        "create T16 backend-neighbor context");
    const uint64_t expected =
        leo2_context_backend(context) == LEO2_BACKEND_AVX2 ? 1U : 0U;
    ExerciseShape(context, 16, 16, LEO2_PROFILE_LEGACY_HIGH_V1,
        64, expected, expected,
        "T16 candidate selection disagrees with backend");
    leo2_context_destroy(context);
}

void ExerciseIrregularLayouts(
    leo2_context* context,
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator)
{
    leo2_codec* codec = CreateCodec(
        context, kSide, kSide, LEO2_PROFILE_LEGACY_HIGH_V1);
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, kShardBytes, &scratch_bytes), LEO2_SUCCESS,
        "query irregular T16 scratch");
    AlignedBuffer scratch(scratch_bytes);
    AlignedBuffer input(kPayloadBytes);
    AlignedBuffer output(kPayloadBytes);
    FillRandom(input.bytes(), input.size(), UINT64_C(0x6972726567756c61));
    std::vector<const void*> original;
    std::vector<void*> recovery;
    SetPackedPointers(input.bytes(), output.bytes(), kSide, kSide,
        kShardBytes, original, recovery);

    std::vector<uint8_t> detached(kShardBytes);
    std::memcpy(&detached[0], original[7], kShardBytes);
    original[7] = &detached[0];
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "encode detached T16 input");
    RequireRouteCounts(0, 0,
        "detached T16 layout entered the packed generated terminal");
    CheckParity(field, generator, kSide, kSide, original, recovery,
        kShardBytes, "detached T16 parity differs from direct generator");

    original[7] = input.bytes() + 7U * kShardBytes;
    std::memset(output.bytes(), 0xa5, output.size());
    recovery[9] = NULL;
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "encode sparse-output T16 neighbor");
    RequireRouteCounts(0, 0,
        "sparse-output T16 neighbor entered a dense candidate route");
    CheckParity(field, generator, kSide, kSide, original, recovery,
        kShardBytes, "sparse-output T16 parity differs from direct generator");
    for (size_t i = 0; i < kShardBytes; ++i)
        Require(output.bytes()[9U * kShardBytes + i] == 0xa5,
            "sparse-output T16 encode modified a null destination");
    leo2_codec_destroy(codec);
}

} // namespace

int main()
{
    try
    {
        Require(leo_init() == 0,
            "initialize legacy Leopard T16 parity oracle");
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
                "T16/B64 three-pass AVX2 test skipped: AVX2 unavailable\n");
            return 0;
        }
        RequireResult(context_result, LEO2_SUCCESS,
            "create T16/B64 AVX2 context");

        const leopard2_test::BinaryField field =
            leopard2_test::make_legacy_gf8();
        const leopard2_test::ProfileLayout layout =
            leopard2_test::make_profile_layout(
                leopard2_test::kLegacyHigh, kSide, kSide);
        const leopard2_test::Matrix generator =
            leopard2_test::direct_systematic_generator(field, layout);

        ExerciseShape(context, 16, 16, LEO2_PROFILE_LEGACY_HIGH_V1,
            64, 0, 0, "default-off T16/B64 candidate route was selected");
        ExerciseBasisAndRandom(context, field, generator);
        ExerciseOffsetsAndBatch(context, field, generator);
        ExerciseNeighborRoutes(context);
        ExerciseIrregularLayouts(context, field, generator);
        leo2_context_destroy(context);

        Require(leopard::ff8::SetHighT16B64ThreePassEnabledForDiagnostics(true),
            "enable T16 backend-route candidate");
        ExerciseBackendRoute(LEO2_BACKEND_AUTO);
        ExerciseBackendRoute(LEO2_BACKEND_SCALAR);
        ExerciseBackendRoute(LEO2_BACKEND_SSSE3);
        ExerciseBackendRoute(LEO2_BACKEND_AVX512);
        ExerciseBackendRoute(LEO2_BACKEND_GFNI);
        Require(leopard::ff8::SetHighT16B64ThreePassEnabledForDiagnostics(false),
            "restore default-off T16 candidate");

        std::printf("T16/B64 three-pass AVX2 qualification checks passed\n");
        return 0;
    }
    catch (const std::exception& error)
    {
        leopard::ff8::SetHighT16B64ThreePassEnabledForDiagnostics(false);
        std::fprintf(stderr,
            "T16/B64 three-pass qualification failure: %s\n", error.what());
        return 1;
    }
}
