/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "direct_oracle.h"
#include "Leopard2Backend.h"
#include "Leopard2Direct.h"
#include "LeopardFF8.h"
#include "leopard.h"
#include "leopard2.h"

#include <chrono>
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

static const unsigned kOriginalCount = 32;
static const unsigned kRecoveryCount = 16;
static const size_t kShardBytes = 64;

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

void FillInput(uint8_t* input)
{
    uint64_t state = UINT64_C(0x5431365132423235);
    for (size_t i = 0; i < kOriginalCount * kShardBytes; ++i)
    {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        input[i] = static_cast<uint8_t>(state >> 24);
    }
}

void FillInput(uint8_t* input, size_t bytes, uint64_t state)
{
    for (size_t i = 0; i < bytes; ++i)
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
    uint32_t original_count,
    uint32_t recovery_count,
    size_t shard_bytes,
    const void** original,
    void** recovery)
{
    for (uint32_t i = 0; i < original_count; ++i)
        original[i] = input + static_cast<size_t>(i) * shard_bytes;
    for (uint32_t i = 0; i < recovery_count; ++i)
        recovery[i] = output + static_cast<size_t>(i) * shard_bytes;
}

uint64_t TwoBlockCalls()
{
    return leopard::ff8::TestOnlyGetHighEncodeCounts().two_block_calls;
}

void CheckDirectOracle(
    const void* const* original,
    void* const* recovery)
{
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh,
            kOriginalCount, kRecoveryCount);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);

    for (unsigned parity = 0; parity < kRecoveryCount; ++parity)
    {
        const std::vector<leopard2_test::Element>& row =
            generator[kOriginalCount + parity];
        const uint8_t* output =
            static_cast<const uint8_t*>(recovery[parity]);
        for (size_t offset = 0; offset < kShardBytes; ++offset)
        {
            leopard2_test::Element expected = 0;
            for (unsigned source = 0; source < kOriginalCount; ++source)
            {
                const uint8_t value =
                    static_cast<const uint8_t*>(original[source])[offset];
                expected = field.add(expected,
                    field.multiply(row[source], value));
            }
            Require(output[offset] == expected,
                "fused T16 q=2 parity differs from direct RS oracle");
        }
    }
}

void CheckDirectOracle(
    uint32_t original_count,
    uint32_t recovery_count,
    size_t shard_bytes,
    const void* const* original,
    void* const* recovery)
{
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh,
            original_count, recovery_count);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);

    for (uint32_t parity = 0; parity < recovery_count; ++parity)
    {
        if (!recovery[parity])
            continue;
        const std::vector<leopard2_test::Element>& row =
            generator[original_count + parity];
        const uint8_t* output =
            static_cast<const uint8_t*>(recovery[parity]);
        for (size_t offset = 0; offset < shard_bytes; ++offset)
        {
            leopard2_test::Element expected = 0;
            for (uint32_t source = 0; source < original_count; ++source)
            {
                const uint8_t value =
                    static_cast<const uint8_t*>(original[source])[offset];
                expected = field.add(expected,
                    field.multiply(row[source], value));
            }
            Require(output[offset] == expected,
                "public parity differs from direct RS oracle");
        }
    }
}

void CheckLegacyParity(
    const void* const* original,
    const uint8_t* expected)
{
    const unsigned work_count =
        leo_encode_work_count(kOriginalCount, kRecoveryCount);
    Require(work_count >= kRecoveryCount,
        "legacy work-count query rejected K32/R16");
    AlignedBuffer legacy_storage(work_count * kShardBytes);
    std::vector<void*> legacy_work(work_count);
    for (unsigned i = 0; i < work_count; ++i)
        legacy_work[i] = legacy_storage.bytes() + i * kShardBytes;
    Require(leo_encode(kShardBytes, kOriginalCount, kRecoveryCount,
            work_count, original, &legacy_work[0]) == Leopard_Success,
        "legacy K32/R16 encode failed");
    Require(std::memcmp(legacy_storage.bytes(), expected,
            kRecoveryCount * kShardBytes) == 0,
        "public terminal parity differs from legacy leo_encode");
}

double TimePublic(
    leo2_codec* codec,
    const void* const* original,
    void* const* recovery,
    void* scratch,
    size_t scratch_bytes,
    uint64_t iterations)
{
    const std::chrono::steady_clock::time_point start =
        std::chrono::steady_clock::now();
    for (uint64_t iteration = 0; iteration < iterations; ++iteration)
    {
        RequireResult(leo2_encode(codec, kShardBytes,
            original, recovery, scratch, scratch_bytes),
            LEO2_SUCCESS, "benchmark generic encode");
    }
    const std::chrono::steady_clock::time_point stop =
        std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::nano>(stop - start).count() /
        static_cast<double>(iterations);
}

void ExerciseAdjacentExclusions(leo2_context* context)
{
    struct Cell
    {
        uint32_t original_count;
        uint32_t recovery_count;
        size_t shard_bytes;
    };
    static const Cell kCells[] = {
        { 32, 16, 66 },
        { 32, 16, 65 },
        { 16, 16, 64 },
        { 33, 16, 64 },
        { 32, 8, 64 },
        { 32, 17, 64 }
    };

    for (size_t cell_index = 0;
         cell_index < sizeof(kCells) / sizeof(kCells[0]); ++cell_index)
    {
        const Cell& cell = kCells[cell_index];
        leo2_codec* codec = NULL;
        RequireResult(leo2_codec_create(context,
            cell.original_count, cell.recovery_count,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
            NULL, &codec), LEO2_SUCCESS, "create exclusion codec");
        size_t scratch_bytes = 0;
        RequireResult(leo2_encode_scratch_size(
            codec, cell.shard_bytes, &scratch_bytes), LEO2_SUCCESS,
            "query exclusion scratch");
        AlignedBuffer scratch(scratch_bytes);
        AlignedBuffer input(
            static_cast<size_t>(cell.original_count) * cell.shard_bytes);
        AlignedBuffer output(
            static_cast<size_t>(cell.recovery_count) * cell.shard_bytes);
        FillInput(input.bytes(), input.size(),
            UINT64_C(0x4558434c55444500) + cell_index);
        std::vector<const void*> original(cell.original_count);
        std::vector<void*> recovery(cell.recovery_count);
        SetPackedPointers(input.bytes(), output.bytes(),
            cell.original_count, cell.recovery_count, cell.shard_bytes,
            &original[0], &recovery[0]);

        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode(codec, cell.shard_bytes,
            &original[0], &recovery[0], scratch.data(), scratch.size()),
            LEO2_SUCCESS, "execute adjacent exclusion");
        Require(TwoBlockCalls() == 0,
            "adjacent K/R/B cell entered the T16 q=2 B64 terminal");
        CheckDirectOracle(cell.original_count, cell.recovery_count,
            cell.shard_bytes, &original[0], &recovery[0]);
        leo2_codec_destroy(codec);
    }
}

void ExerciseEligibleMatrix(leo2_context* context)
{
    static const size_t kEligibleByteCounts[] = {
        1, 2, 3, 7, 8, 15, 16, 17, 31, 32, 33, 63, 64
    };
    for (size_t byte_index = 0;
         byte_index < sizeof(kEligibleByteCounts) /
             sizeof(kEligibleByteCounts[0]); ++byte_index)
    {
        const size_t shard_bytes = kEligibleByteCounts[byte_index];
        for (uint32_t original_count = 17;
             original_count <= 32; ++original_count)
        {
            for (uint32_t recovery_count = 9;
                 recovery_count <= 16; ++recovery_count)
            {
                leo2_codec* codec = NULL;
                RequireResult(leo2_codec_create(context,
                    original_count, recovery_count,
                    LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
                    NULL, &codec), LEO2_SUCCESS,
                    "create eligible T16 q=2 codec");
                size_t scratch_bytes = 0;
                RequireResult(leo2_encode_scratch_size(
                    codec, shard_bytes, &scratch_bytes), LEO2_SUCCESS,
                    "query eligible T16 q=2 scratch");
                AlignedBuffer scratch(scratch_bytes);
                AlignedBuffer input(
                    static_cast<size_t>(original_count) * shard_bytes);
                AlignedBuffer output(
                    static_cast<size_t>(recovery_count) * shard_bytes);
                FillInput(input.bytes(), input.size(),
                    UINT64_C(0x5132543136423634) +
                        original_count * 257U + recovery_count);
                std::vector<const void*> original(original_count);
                std::vector<void*> recovery(recovery_count);
                SetPackedPointers(input.bytes(), output.bytes(),
                    original_count, recovery_count, shard_bytes,
                    &original[0], &recovery[0]);

                leopard::ff8::TestOnlyResetHighEncodeCounts();
                RequireResult(leo2_encode(codec, shard_bytes,
                    &original[0], &recovery[0],
                    scratch.data(), scratch.size()), LEO2_SUCCESS,
                    "execute eligible T16 q=2 terminal");
                Require(TwoBlockCalls() == 1,
                    "eligible T16 q=2 encode missed terminal");
                CheckDirectOracle(original_count, recovery_count,
                    shard_bytes, &original[0], &recovery[0]);
                const std::vector<uint8_t> candidate(
                    output.bytes(), output.bytes() + output.size());

                RequireResult(leo2_test_codec_set_encode_mode(
                    codec, LEO2_TEST_ENCODE_FORCE_TRANSFORM), LEO2_SUCCESS,
                    "force eligible T16 q=2 mature transform");
                std::memset(output.bytes(), 0xa5, output.size());
                leopard::ff8::TestOnlyResetHighEncodeCounts();
                RequireResult(leo2_encode(codec, shard_bytes,
                    &original[0], &recovery[0],
                    scratch.data(), scratch.size()), LEO2_SUCCESS,
                    "execute eligible T16 q=2 mature transform");
                Require(TwoBlockCalls() == 0,
                    "forced mature transform entered T16 q=2 terminal");
                Require(std::memcmp(output.bytes(), &candidate[0],
                        candidate.size()) == 0,
                    "eligible T16 q=2 terminal differs from mature transform");
                leo2_codec_destroy(codec);
            }
        }
    }
}

void ExerciseStagedTailInterfaces(leo2_context* context)
{
    static const uint32_t kTailOriginalCount = 17;
    static const uint32_t kTailRecoveryCount = 9;
    static const size_t kTailBytes = 33;
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context,
        kTailOriginalCount, kTailRecoveryCount,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        NULL, &codec), LEO2_SUCCESS, "create staged-tail codec");
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, kTailBytes, &scratch_bytes), LEO2_SUCCESS,
        "query staged-tail scratch");
    AlignedBuffer scratch(scratch_bytes);
    AlignedBuffer input(kTailOriginalCount * kTailBytes + 1);
    AlignedBuffer output(kTailRecoveryCount * kTailBytes + 3);
    FillInput(input.bytes() + 1, kTailOriginalCount * kTailBytes,
        UINT64_C(0x5354414745445441));
    std::vector<const void*> original(kTailOriginalCount);
    std::vector<void*> recovery(kTailRecoveryCount);
    SetPackedPointers(input.bytes() + 1, output.bytes() + 3,
        kTailOriginalCount, kTailRecoveryCount, kTailBytes,
        &original[0], &recovery[0]);

    leo2_encode_batch_item item = {
        kTailBytes, &original[0], &recovery[0],
        scratch.data(), scratch.size()
    };
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch(codec, &item, 1), LEO2_SUCCESS,
        "execute staged-tail batch terminal");
    Require(TwoBlockCalls() == 1,
        "staged-tail batch missed the fused terminal");
    CheckDirectOracle(kTailOriginalCount, kTailRecoveryCount,
        kTailBytes, &original[0], &recovery[0]);

    std::memset(output.bytes(), 0xa5, output.size());
    recovery[4] = NULL;
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kTailBytes,
        &original[0], &recovery[0], scratch.data(), scratch.size()),
        LEO2_SUCCESS, "execute staged-tail sparse fallback");
    Require(TwoBlockCalls() == 0,
        "staged-tail sparse output entered the packed terminal");
    CheckDirectOracle(kTailOriginalCount, kTailRecoveryCount,
        kTailBytes, &original[0], &recovery[0]);
    for (size_t i = 0; i < kTailBytes; ++i)
        Require(output.bytes()[3 + 4 * kTailBytes + i] == 0xa5,
            "staged-tail sparse fallback modified null output");
    recovery[4] = output.bytes() + 3 + 4 * kTailBytes;

    std::vector<uint8_t> detached(kTailBytes);
    std::memcpy(&detached[0], original.back(), kTailBytes);
    const void* const packed_last = original.back();
    original.back() = &detached[0];
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kTailBytes,
        &original[0], &recovery[0], scratch.data(), scratch.size()),
        LEO2_SUCCESS, "execute staged-tail detached fallback");
    Require(TwoBlockCalls() == 0,
        "staged-tail detached input entered the packed terminal");
    CheckDirectOracle(kTailOriginalCount, kTailRecoveryCount,
        kTailBytes, &original[0], &recovery[0]);
    original.back() = packed_last;

    const std::vector<uint8_t> output_before(
        output.bytes(), output.bytes() + output.size());
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kTailBytes,
        &original[0], &recovery[0], scratch.data(), scratch.size() - 1),
        LEO2_SCRATCH_TOO_SMALL, "reject undersized staged-tail scratch");
    Require(TwoBlockCalls() == 0,
        "undersized staged-tail scratch reached fused arithmetic");
    Require(std::memcmp(output.bytes(), &output_before[0],
            output_before.size()) == 0,
        "undersized staged-tail scratch modified output");

    for (uint32_t i = 0; i < kTailRecoveryCount; ++i)
        recovery[i] = input.bytes() + 1 + i * kTailBytes;
    const std::vector<uint8_t> input_before(
        input.bytes(), input.bytes() + input.size());
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kTailBytes,
        &original[0], &recovery[0], scratch.data(), scratch.size()),
        LEO2_OVERLAP, "reject overlapping staged-tail slabs");
    Require(TwoBlockCalls() == 0,
        "staged-tail overlap reached fused arithmetic");
    Require(std::memcmp(input.bytes(), &input_before[0],
            input_before.size()) == 0,
        "staged-tail overlap rejection modified input");
    leo2_codec_destroy(codec);
}

void ExerciseScalarExclusion()
{
    leo2_context_options options = {};
    options.struct_size = sizeof(options);
    options.backend = LEO2_BACKEND_SCALAR;
    options.thread_count = 1;
    leo2_context* context = NULL;
    RequireResult(leo2_context_create(&options, &context), LEO2_SUCCESS,
        "create scalar exclusion context");
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context,
        kOriginalCount, kRecoveryCount,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        NULL, &codec), LEO2_SUCCESS, "create scalar exclusion codec");
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, kShardBytes, &scratch_bytes), LEO2_SUCCESS,
        "query scalar exclusion scratch");
    AlignedBuffer scratch(scratch_bytes);
    AlignedBuffer input(kOriginalCount * kShardBytes);
    AlignedBuffer output(kRecoveryCount * kShardBytes);
    FillInput(input.bytes());
    const void* original[kOriginalCount];
    void* recovery[kRecoveryCount];
    SetPackedPointers(input.bytes(), output.bytes(),
        kOriginalCount, kRecoveryCount, kShardBytes, original, recovery);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes,
        original, recovery, scratch.data(), scratch.size()),
        LEO2_SUCCESS, "execute scalar exclusion");
    Require(TwoBlockCalls() == 0,
        "scalar context entered AVX2 K32/R16 terminal");
    CheckDirectOracle(original, recovery);
    leo2_codec_destroy(codec);
    leo2_context_destroy(context);
}

void Exercise(leo2_context* context, uint64_t benchmark_iterations)
{
    Require(leopard2_internal::
            SetHighT16Q2B64FusedEnabledForDiagnostics(true),
        "enable K32/R16 terminal");
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context,
        kOriginalCount, kRecoveryCount,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        NULL, &codec), LEO2_SUCCESS, "create K32/R16 codec");

    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, kShardBytes, &scratch_bytes), LEO2_SUCCESS,
        "query K32/R16 scratch");
    Require(scratch_bytes >= 3328,
        "K32/R16 scratch query lost the ordinary transform bound");
    AlignedBuffer scratch(scratch_bytes);
    AlignedBuffer input(kOriginalCount * kShardBytes);
    AlignedBuffer generic_output(kRecoveryCount * kShardBytes);
    AlignedBuffer fused_work(
        2 * kRecoveryCount * kShardBytes);
    FillInput(input.bytes());

    const void* original[kOriginalCount];
    void* generic_recovery[kRecoveryCount];
    for (unsigned row = 0; row < kOriginalCount; ++row)
        original[row] = input.bytes() + row * kShardBytes;
    for (unsigned row = 0; row < kRecoveryCount; ++row)
        generic_recovery[row] =
            generic_output.bytes() + row * kShardBytes;

    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes,
        original, generic_recovery, scratch.data(), scratch.size()),
        LEO2_SUCCESS, "execute public K32/R16 terminal");
    Require(TwoBlockCalls() == 1,
        "public K32/R16 encode missed the fused terminal");
    const std::vector<uint8_t> candidate_output(
        generic_output.bytes(),
        generic_output.bytes() + kRecoveryCount * kShardBytes);
    leopard::backend::AVX2FF8HighEncodeT16Q2B64Fused(
        input.bytes(), fused_work.bytes(),
        fused_work.bytes() + kRecoveryCount * kShardBytes,
        kOriginalCount, kRecoveryCount);

    Require(std::memcmp(&candidate_output[0], fused_work.bytes(),
        kRecoveryCount * kShardBytes) == 0,
        "public terminal parity differs from direct fused callback");
    void* fused_recovery[kRecoveryCount];
    for (unsigned row = 0; row < kRecoveryCount; ++row)
        fused_recovery[row] = fused_work.bytes() + row * kShardBytes;
    CheckDirectOracle(original, fused_recovery);
    CheckLegacyParity(original, generic_output.bytes());

    std::memset(generic_output.bytes(), 0xa5,
        kRecoveryCount * kShardBytes);
    leo2_encode_batch_item item = {
        kShardBytes, original, generic_recovery,
        scratch.data(), scratch.size()
    };
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch(codec, &item, 1), LEO2_SUCCESS,
        "execute one-item K32/R16 batch terminal");
    Require(TwoBlockCalls() == 1,
        "one-item K32/R16 batch missed the fused terminal");
    CheckDirectOracle(original, generic_recovery);

    AlignedBuffer unaligned_input(kOriginalCount * kShardBytes + 1);
    AlignedBuffer unaligned_output(kRecoveryCount * kShardBytes + 3);
    const void* unaligned_original[kOriginalCount];
    void* unaligned_recovery[kRecoveryCount];
    FillInput(unaligned_input.bytes() + 1,
        kOriginalCount * kShardBytes, UINT64_C(0x554e414c49474e45));
    SetPackedPointers(unaligned_input.bytes() + 1,
        unaligned_output.bytes() + 3, kOriginalCount, kRecoveryCount,
        kShardBytes, unaligned_original, unaligned_recovery);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes,
        unaligned_original, unaligned_recovery,
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute unaligned packed K32/R16 terminal");
    Require(TwoBlockCalls() == 1,
        "unaligned packed K32/R16 encode missed terminal");
    CheckDirectOracle(unaligned_original, unaligned_recovery);

    RequireResult(leo2_test_codec_set_encode_mode(
        codec, LEO2_TEST_ENCODE_FORCE_TRANSFORM), LEO2_SUCCESS,
        "force mature K32/R16 transform");
    std::memset(generic_output.bytes(), 0xa5,
        kRecoveryCount * kShardBytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes,
        original, generic_recovery, scratch.data(), scratch.size()),
        LEO2_SUCCESS, "execute forced mature K32/R16 transform");
    Require(TwoBlockCalls() == 0,
        "forced transform entered the K32/R16 terminal");
    Require(std::memcmp(generic_output.bytes(), &candidate_output[0],
            candidate_output.size()) == 0,
        "forced mature transform differs from fused terminal");
    CheckDirectOracle(original, generic_recovery);
    RequireResult(leo2_test_codec_set_encode_mode(
        codec, LEO2_TEST_ENCODE_AUTO), LEO2_SUCCESS,
        "restore K32/R16 AUTO encode mode");

    Require(leopard2_internal::
            SetHighT16Q2B64FusedEnabledForDiagnostics(false),
        "disable K32/R16 terminal control");
    std::memset(generic_output.bytes(), 0xa5,
        kRecoveryCount * kShardBytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes,
        original, generic_recovery, scratch.data(), scratch.size()),
        LEO2_SUCCESS, "execute same-binary K32/R16 control");
    Require(TwoBlockCalls() == 0,
        "same-binary control entered the K32/R16 terminal");
    Require(std::memcmp(generic_output.bytes(), &candidate_output[0],
            candidate_output.size()) == 0,
        "same-binary control differs from fused terminal");
    CheckDirectOracle(original, generic_recovery);
    Require(leopard2_internal::
            SetHighT16Q2B64FusedEnabledForDiagnostics(true),
        "restore K32/R16 terminal");

    std::memset(generic_output.bytes(), 0xa5,
        kRecoveryCount * kShardBytes);
    generic_recovery[3] = NULL;
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes,
        original, generic_recovery, scratch.data(), scratch.size()),
        LEO2_SUCCESS, "execute sparse K32/R16 fallback");
    Require(TwoBlockCalls() == 0,
        "sparse output entered the dense K32/R16 terminal");
    CheckDirectOracle(kOriginalCount, kRecoveryCount, kShardBytes,
        original, generic_recovery);
    for (size_t i = 0; i < kShardBytes; ++i)
    {
        Require(generic_output.bytes()[3 * kShardBytes + i] == 0xa5,
            "sparse fallback modified null K32/R16 parity");
    }
    generic_recovery[3] = generic_output.bytes() + 3 * kShardBytes;

    std::vector<uint8_t> detached(kShardBytes);
    std::memcpy(&detached[0], original[31], kShardBytes);
    const void* const packed_last = original[31];
    original[31] = &detached[0];
    std::memset(generic_output.bytes(), 0xa5,
        kRecoveryCount * kShardBytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes,
        original, generic_recovery, scratch.data(), scratch.size()),
        LEO2_SUCCESS, "execute non-packed K32/R16 fallback");
    Require(TwoBlockCalls() == 0,
        "non-packed input entered the packed K32/R16 terminal");
    CheckDirectOracle(original, generic_recovery);
    original[31] = packed_last;

    std::memset(generic_output.bytes(), 0xa5,
        kRecoveryCount * kShardBytes);
    const std::vector<uint8_t> output_before(
        generic_output.bytes(),
        generic_output.bytes() + kRecoveryCount * kShardBytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes,
        original, generic_recovery, scratch.data(), scratch.size() - 1),
        LEO2_SCRATCH_TOO_SMALL, "reject undersized K32/R16 scratch");
    Require(TwoBlockCalls() == 0,
        "undersized scratch reached K32/R16 terminal arithmetic");
    Require(std::memcmp(generic_output.bytes(), &output_before[0],
            output_before.size()) == 0,
        "undersized scratch modified K32/R16 output");

    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes,
        original, generic_recovery, scratch.bytes() + 1, scratch.size()),
        LEO2_BAD_ALIGNMENT, "reject misaligned K32/R16 scratch");
    Require(TwoBlockCalls() == 0,
        "misaligned scratch reached K32/R16 terminal arithmetic");

    for (unsigned i = 0; i < kRecoveryCount; ++i)
        generic_recovery[i] = input.bytes() + i * kShardBytes;
    const std::vector<uint8_t> input_before(
        input.bytes(), input.bytes() + input.size());
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, kShardBytes,
        original, generic_recovery, scratch.data(), scratch.size()),
        LEO2_OVERLAP, "reject overlapping K32/R16 slabs");
    Require(TwoBlockCalls() == 0,
        "overlap reached K32/R16 terminal arithmetic");
    Require(std::memcmp(input.bytes(), &input_before[0],
            input_before.size()) == 0,
        "overlap rejection modified K32/R16 input");
    for (unsigned i = 0; i < kRecoveryCount; ++i)
        generic_recovery[i] = generic_output.bytes() + i * kShardBytes;

    alignas(64) uint8_t protected_storage[
        kRecoveryCount * kShardBytes] = {};
    leo2_encode_batch_item* const protected_item =
        new (protected_storage) leo2_encode_batch_item;
    protected_item->shard_bytes = kShardBytes;
    protected_item->original = original;
    protected_item->recovery = generic_recovery;
    protected_item->scratch = scratch.data();
    protected_item->scratch_bytes = scratch.size();
    for (unsigned i = 0; i < kRecoveryCount; ++i)
        generic_recovery[i] = protected_storage + i * kShardBytes;
    const std::vector<uint8_t> protected_before(
        protected_storage, protected_storage + sizeof(protected_storage));
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch(codec, protected_item, 1),
        LEO2_OVERLAP, "reject K32/R16 output/batch-metadata overlap");
    Require(TwoBlockCalls() == 0,
        "batch metadata overlap reached K32/R16 terminal arithmetic");
    Require(std::memcmp(protected_storage, &protected_before[0],
            protected_before.size()) == 0,
        "metadata overlap rejection modified K32/R16 batch storage");
    for (unsigned i = 0; i < kRecoveryCount; ++i)
        generic_recovery[i] = generic_output.bytes() + i * kShardBytes;

    if (benchmark_iterations != 0)
    {
        static const uint64_t kWarmupIterations = 4096;
        static const unsigned kRounds = 9;
        leopard2_internal::
            SetHighT16Q2B64FusedEnabledForDiagnostics(false);
        TimePublic(codec, original, generic_recovery,
            scratch.data(), scratch.size(), kWarmupIterations);
        leopard2_internal::
            SetHighT16Q2B64FusedEnabledForDiagnostics(true);
        TimePublic(codec, original, generic_recovery,
            scratch.data(), scratch.size(), kWarmupIterations);
        for (unsigned round = 0; round < kRounds; ++round)
        {
            double generic_ns;
            double fused_ns;
            if ((round & 1U) == 0)
            {
                leopard2_internal::
                    SetHighT16Q2B64FusedEnabledForDiagnostics(false);
                generic_ns = TimePublic(codec, original, generic_recovery,
                    scratch.data(), scratch.size(), benchmark_iterations);
                leopard2_internal::
                    SetHighT16Q2B64FusedEnabledForDiagnostics(true);
                fused_ns = TimePublic(codec, original, generic_recovery,
                    scratch.data(), scratch.size(), benchmark_iterations);
            }
            else
            {
                leopard2_internal::
                    SetHighT16Q2B64FusedEnabledForDiagnostics(true);
                fused_ns = TimePublic(codec, original, generic_recovery,
                    scratch.data(), scratch.size(), benchmark_iterations);
                leopard2_internal::
                    SetHighT16Q2B64FusedEnabledForDiagnostics(false);
                generic_ns = TimePublic(codec, original, generic_recovery,
                    scratch.data(), scratch.size(), benchmark_iterations);
            }
            std::printf(
                "round=%u generic_ns=%.6f fused_ns=%.6f speedup=%.6f\n",
                round, generic_ns, fused_ns, generic_ns / fused_ns);
        }
        leopard2_internal::
            SetHighT16Q2B64FusedEnabledForDiagnostics(true);
    }
    leo2_codec_destroy(codec);
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        uint64_t benchmark_iterations = 0;
        if (argc == 3 && std::strcmp(argv[1], "--benchmark") == 0)
        {
            char* end = NULL;
            const unsigned long long parsed = std::strtoull(argv[2], &end, 10);
            Require(end != argv[2] && *end == '\0' && parsed != 0,
                "invalid benchmark iteration count");
            benchmark_iterations = static_cast<uint64_t>(parsed);
        }
        else
            Require(argc == 1, "usage: test [--benchmark iterations]");

        leo2_context_options options = {};
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AVX2;
        options.thread_count = 1;
        leo2_context* context = NULL;
        const leo2_result result = leo2_context_create(&options, &context);
        if (result == LEO2_UNSUPPORTED)
        {
            std::printf("T16 q=2 fused AVX2 test skipped: AVX2 unavailable\n");
            return 0;
        }
        RequireResult(result, LEO2_SUCCESS, "create AVX2 context");
        Exercise(context, benchmark_iterations);
        ExerciseEligibleMatrix(context);
        ExerciseStagedTailInterfaces(context);
        ExerciseAdjacentExclusions(context);
        leo2_context_destroy(context);
        ExerciseScalarExclusion();
        std::printf("T16 q=2 B64 fused AVX2 checks passed\n");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "T16 q=2 B64 fused failure: %s\n",
            error.what());
        return 1;
    }
}
