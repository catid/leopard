/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
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

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

namespace {

static const unsigned kOriginalCount = 9;
static const unsigned kMinimumRecoveryCount = 5;
static const unsigned kMaximumRecoveryCount = 8;

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

void FillInput(uint8_t* input, size_t shard_bytes)
{
    uint64_t state = UINT64_C(0x4b39523542323536);
    for (size_t i = 0; i < kOriginalCount * shard_bytes; ++i)
    {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        input[i] = static_cast<uint8_t>(state >> 24);
    }
}

void SetPackedPointers(
    uint8_t* input_base,
    uint8_t* output_base,
    const void** original,
    void** recovery,
    unsigned recovery_count,
    size_t shard_bytes)
{
    for (unsigned i = 0; i < kOriginalCount; ++i)
        original[i] = input_base + i * shard_bytes;
    for (unsigned i = 0; i < recovery_count; ++i)
        recovery[i] = output_base + i * shard_bytes;
}

void ResetOutput(
    uint8_t* output,
    unsigned recovery_count,
    size_t shard_bytes)
{
    std::memset(output, 0xa5, recovery_count * shard_bytes);
}

void CheckParity(
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator,
    const void* const* original,
    void* const* recovery,
    unsigned recovery_count,
    size_t shard_bytes)
{
    for (unsigned parity = 0; parity < recovery_count; ++parity)
    {
        if (!recovery[parity])
            continue;
        uint8_t* const output = static_cast<uint8_t*>(recovery[parity]);
        const std::vector<leopard2_test::Element>& row =
            generator[kOriginalCount + parity];
        for (size_t offset = 0; offset < shard_bytes; ++offset)
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
                "K9 parity differs from the independent oracle");
        }
    }
}

void RequireTerminalCalls(
    unsigned recovery_count,
    uint64_t expected,
    const char* message)
{
    const leopard::ff8::TestOnlyHighEncodeCounts counts =
        leopard::ff8::TestOnlyGetHighEncodeCounts();
    if (recovery_count == 5)
    {
        Require(counts.k9r5_tail_calls == expected, message);
        Require(counts.k9r6r8_tail_calls == 0,
            "R=5 encode entered the R=6..8 terminal counter");
    }
    else
    {
        Require(counts.k9r6r8_tail_calls == expected, message);
        Require(counts.k9r5_tail_calls == 0,
            "R=6..8 encode entered the R=5 terminal counter");
    }
}

bool SetTerminalEnabledForDiagnostics(
    unsigned recovery_count,
    size_t shard_bytes,
    bool enabled)
{
    if (shard_bytes == 1024)
    {
        Require(recovery_count == 5,
            "only K9/R5 owns the 1024-byte exact terminal");
        return leopard2_internal::
            SetK9R5B1024TerminalEnabledForDiagnostics(enabled);
    }
    return recovery_count == 5
        ? leopard2_internal::
            SetK9R5B256TerminalEnabledForDiagnostics(enabled)
        : leopard2_internal::
            SetK9R6R8B256TerminalEnabledForDiagnostics(enabled);
}

void ExerciseAVX2Terminal(
    leo2_context* context,
    unsigned recovery_count,
    size_t shard_bytes)
{
    Require(SetTerminalEnabledForDiagnostics(
        recovery_count, shard_bytes, true),
        "enable the K9 terminal for its focused test");
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context,
        kOriginalCount, recovery_count,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        NULL, &codec), LEO2_SUCCESS, "create K9 terminal codec");

    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, shard_bytes, &scratch_bytes), LEO2_SUCCESS,
        "query K9 terminal scratch");
    Require(scratch_bytes != 0, "K9 scratch query returned zero");
    AlignedBuffer scratch(scratch_bytes + leo2_scratch_alignment());
    AlignedBuffer input(kOriginalCount * shard_bytes + 8);
    AlignedBuffer output(recovery_count * shard_bytes + 8);
    FillInput(input.bytes(), shard_bytes);

    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh,
            kOriginalCount, recovery_count);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);

    alignas(64) const void* original[kOriginalCount];
    alignas(64) void* recovery[kMaximumRecoveryCount];
    SetPackedPointers(
        input.bytes(), output.bytes(), original, recovery, recovery_count,
        shard_bytes);

    ResetOutput(output.bytes(), recovery_count, shard_bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute packed public K9 terminal");
    RequireTerminalCalls(recovery_count, 1,
        "packed public encode did not select the K9 terminal");
    CheckParity(field, generator, original, recovery, recovery_count,
        shard_bytes);

    ResetOutput(output.bytes(), recovery_count, shard_bytes);
    leo2_encode_batch_item item = {
        shard_bytes, original, recovery, scratch.data(), scratch_bytes
    };
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch(codec, &item, 1), LEO2_SUCCESS,
        "execute packed one-item K9 batch terminal");
    RequireTerminalCalls(recovery_count, 1,
        "one-item batch did not select the K9 terminal");
    CheckParity(field, generator, original, recovery, recovery_count,
        shard_bytes);

    SetPackedPointers(
        input.bytes() + 1, output.bytes() + 3, original, recovery,
        recovery_count, shard_bytes);
    FillInput(input.bytes() + 1, shard_bytes);
    ResetOutput(output.bytes() + 3, recovery_count, shard_bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute unaligned packed K9 terminal");
    RequireTerminalCalls(recovery_count, 1,
        "unaligned packed encode missed the K9 terminal");
    CheckParity(field, generator, original, recovery, recovery_count,
        shard_bytes);

    SetPackedPointers(
        input.bytes(), output.bytes(), original, recovery, recovery_count,
        shard_bytes);
    FillInput(input.bytes(), shard_bytes);
    ResetOutput(output.bytes(), recovery_count, shard_bytes);
    RequireResult(leo2_test_codec_set_encode_mode(
        codec, LEO2_TEST_ENCODE_FORCE_TRANSFORM), LEO2_SUCCESS,
        "force mature K9 transform");
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute forced mature K9 transform");
    RequireTerminalCalls(recovery_count, 0,
        "forced transform unexpectedly entered the K9 terminal");
    CheckParity(field, generator, original, recovery, recovery_count,
        shard_bytes);
    RequireResult(leo2_test_codec_set_encode_mode(
        codec, LEO2_TEST_ENCODE_AUTO), LEO2_SUCCESS,
        "restore K9 AUTO encode mode");

    Require(SetTerminalEnabledForDiagnostics(
        recovery_count, shard_bytes, false),
        "disable the K9 benchmark terminal");
    ResetOutput(output.bytes(), recovery_count, shard_bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute same-binary K9 control route");
    RequireTerminalCalls(recovery_count, 0,
        "same-binary control unexpectedly entered the K9 terminal");
    CheckParity(field, generator, original, recovery, recovery_count,
        shard_bytes);

    leo2_encode_batch_binding* disabled_binding = NULL;
    if (shard_bytes == 1024)
    {
        ResetOutput(output.bytes(), recovery_count, shard_bytes);
        RequireResult(leo2_encode_batch_binding_create(
            codec, &item, 1, &disabled_binding), LEO2_SUCCESS,
            "create disabled K9/1024 reusable binding");
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode_batch_binding_execute(disabled_binding),
            LEO2_SUCCESS, "execute disabled K9/1024 reusable binding");
        RequireTerminalCalls(recovery_count, 0,
            "disabled reusable binding entered the K9/1024 terminal");
        CheckParity(field, generator, original, recovery, recovery_count,
            shard_bytes);
    }
    Require(SetTerminalEnabledForDiagnostics(
        recovery_count, shard_bytes, true),
        "restore the K9 benchmark terminal");

    if (shard_bytes == 1024)
    {
        /* A binding snapshots its route at construction.  Re-enabling the
           process-local diagnostic word must not mutate the old object. */
        ResetOutput(output.bytes(), recovery_count, shard_bytes);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode_batch_binding_execute(disabled_binding),
            LEO2_SUCCESS, "re-execute disabled K9/1024 reusable binding");
        RequireTerminalCalls(recovery_count, 0,
            "existing disabled binding observed the restored selector");
        CheckParity(field, generator, original, recovery, recovery_count,
            shard_bytes);
        leo2_encode_batch_binding_destroy(disabled_binding);

        ResetOutput(output.bytes(), recovery_count, shard_bytes);
        leo2_encode_batch_binding* enabled_binding = NULL;
        RequireResult(leo2_encode_batch_binding_create(
            codec, &item, 1, &enabled_binding), LEO2_SUCCESS,
            "create enabled K9/1024 reusable binding");
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode_batch_binding_execute(enabled_binding),
            LEO2_SUCCESS, "execute enabled K9/1024 reusable binding");
        RequireTerminalCalls(recovery_count, 1,
            "enabled reusable binding missed the K9/1024 terminal");
        CheckParity(field, generator, original, recovery, recovery_count,
            shard_bytes);

        /* The binding owns its route snapshot.  Disabling later one-shot and
           binding construction must not alter this existing object. */
        Require(SetTerminalEnabledForDiagnostics(
            recovery_count, shard_bytes, false),
            "disable selector after enabled K9/1024 binding creation");
        input.bytes()[0] ^= UINT8_C(0x5a);
        ResetOutput(output.bytes(), recovery_count, shard_bytes);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode_batch_binding_execute(enabled_binding),
            LEO2_SUCCESS,
            "execute enabled K9/1024 binding after selector disable");
        RequireTerminalCalls(recovery_count, 1,
            "enabled binding observed a later disabled selector");
        CheckParity(field, generator, original, recovery, recovery_count,
            shard_bytes);
        leo2_encode_batch_binding_destroy(enabled_binding);
        Require(SetTerminalEnabledForDiagnostics(
            recovery_count, shard_bytes, true),
            "restore selector after enabled K9/1024 binding test");
    }

    ResetOutput(output.bytes(), recovery_count, shard_bytes);
    recovery[3] = NULL;
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute sparse K9 fallback");
    RequireTerminalCalls(recovery_count, 0,
        "sparse output entered the dense K9 terminal");
    CheckParity(field, generator, original, recovery, recovery_count,
        shard_bytes);
    for (size_t i = 0; i < shard_bytes; ++i)
        Require(output.bytes()[3 * shard_bytes + i] == 0xa5,
            "sparse fallback modified a null parity shard");

    SetPackedPointers(
        input.bytes(), output.bytes(), original, recovery, recovery_count,
        shard_bytes);
    std::vector<uint8_t> detached(shard_bytes);
    std::memcpy(&detached[0], original[8], shard_bytes);
    original[8] = &detached[0];
    ResetOutput(output.bytes(), recovery_count, shard_bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute non-packed K9 fallback");
    RequireTerminalCalls(recovery_count, 0,
        "non-packed source entered the packed K9 terminal");
    CheckParity(field, generator, original, recovery, recovery_count,
        shard_bytes);

    SetPackedPointers(
        input.bytes(), output.bytes(), original, recovery, recovery_count,
        shard_bytes);
    ResetOutput(output.bytes(), recovery_count, shard_bytes);
    const std::vector<uint8_t> output_before(
        output.bytes(), output.bytes() + recovery_count * shard_bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, original, recovery,
        scratch.data(), scratch_bytes - 1), LEO2_SCRATCH_TOO_SMALL,
        "reject undersized K9 terminal scratch");
    RequireTerminalCalls(recovery_count, 0,
        "undersized scratch reached the K9 terminal");
    Require(std::memcmp(output.bytes(), &output_before[0],
            output_before.size()) == 0,
        "undersized scratch modified K9 output");

    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, original, recovery,
        scratch.bytes() + 1, scratch_bytes), LEO2_BAD_ALIGNMENT,
        "reject misaligned K9 terminal scratch");
    RequireTerminalCalls(recovery_count, 0,
        "misaligned scratch reached the K9 terminal");

    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, original, recovery,
        input.bytes(), scratch_bytes), LEO2_OVERLAP,
        "reject K9 terminal scratch overlapping packed input");
    RequireTerminalCalls(recovery_count, 0,
        "input-overlapping scratch reached the K9 terminal");
    Require(std::equal(output_before.begin(), output_before.end(),
            output.bytes()),
        "input-overlapping scratch modified K9 output");

    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, original, recovery,
        output.bytes(), scratch_bytes), LEO2_OVERLAP,
        "reject K9 terminal scratch overlapping packed output");
    RequireTerminalCalls(recovery_count, 0,
        "output-overlapping scratch reached the K9 terminal");
    Require(std::equal(output_before.begin(), output_before.end(),
            output.bytes()),
        "output-overlapping scratch modified K9 output");

    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, original, recovery,
        static_cast<void*>(original), scratch_bytes), LEO2_OVERLAP,
        "reject K9 terminal scratch overlapping input pointer metadata");
    RequireTerminalCalls(recovery_count, 0,
        "metadata-overlapping scratch reached the K9 terminal");
    Require(std::equal(output_before.begin(), output_before.end(),
            output.bytes()),
        "metadata-overlapping scratch modified K9 output");

    const uintptr_t overflow_address =
        std::numeric_limits<uintptr_t>::max() &
        ~static_cast<uintptr_t>(leo2_scratch_alignment() - 1U);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, original, recovery,
        reinterpret_cast<void*>(overflow_address), scratch_bytes),
        LEO2_INVALID_ARGUMENT,
        "reject overflowing K9 terminal scratch span");
    RequireTerminalCalls(recovery_count, 0,
        "overflowing scratch reached the K9 terminal");
    Require(std::equal(output_before.begin(), output_before.end(),
            output.bytes()),
        "overflowing scratch modified K9 output");

    for (unsigned i = 0; i < recovery_count; ++i)
        recovery[i] = input.bytes() + i * shard_bytes;
    const std::vector<uint8_t> input_before(
        input.bytes(), input.bytes() + kOriginalCount * shard_bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_OVERLAP,
        "reject overlapping packed K9 slabs");
    RequireTerminalCalls(recovery_count, 0,
        "overlapping slabs reached the K9 terminal");
    Require(std::memcmp(input.bytes(), &input_before[0],
            input_before.size()) == 0,
        "overlap rejection modified K9 input");

    SetPackedPointers(
        input.bytes(), output.bytes(), original, recovery, recovery_count,
        shard_bytes);
    AlignedBuffer protected_storage(
        kMaximumRecoveryCount * shard_bytes);
    leo2_encode_batch_item* const protected_item =
        new (protected_storage.bytes()) leo2_encode_batch_item;
    protected_item->shard_bytes = shard_bytes;
    protected_item->original = original;
    protected_item->recovery = recovery;
    protected_item->scratch = scratch.data();
    protected_item->scratch_bytes = scratch_bytes;
    for (unsigned i = 0; i < recovery_count; ++i)
        recovery[i] = protected_storage.bytes() + i * shard_bytes;
    const leo2_encode_batch_item protected_before = *protected_item;
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch(codec, protected_item, 1),
        LEO2_OVERLAP, "reject K9 output/batch-metadata overlap");
    RequireTerminalCalls(recovery_count, 0,
        "metadata overlap reached the K9 terminal");
    Require(std::memcmp(protected_item, &protected_before,
            sizeof(*protected_item)) == 0,
        "metadata rejection modified the batch descriptor");

    leo2_codec_destroy(codec);
}

void ExerciseScalarFallback(
    unsigned recovery_count,
    size_t shard_bytes)
{
    leo2_context_options options = {};
    options.struct_size = sizeof(options);
    options.backend = LEO2_BACKEND_SCALAR;
    options.thread_count = 1;
    leo2_context* context = NULL;
    RequireResult(leo2_context_create(&options, &context), LEO2_SUCCESS,
        "create scalar K9 fallback context");
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context,
        kOriginalCount, recovery_count,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        NULL, &codec), LEO2_SUCCESS, "create scalar K9 codec");

    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, shard_bytes, &scratch_bytes), LEO2_SUCCESS,
        "query scalar K9 scratch");
    AlignedBuffer scratch(scratch_bytes);
    AlignedBuffer input(kOriginalCount * shard_bytes);
    AlignedBuffer output(recovery_count * shard_bytes);
    FillInput(input.bytes(), shard_bytes);
    const void* original[kOriginalCount];
    void* recovery[kMaximumRecoveryCount];
    SetPackedPointers(
        input.bytes(), output.bytes(), original, recovery, recovery_count,
        shard_bytes);

    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, shard_bytes, original, recovery,
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute scalar K9 fallback");
    RequireTerminalCalls(recovery_count, 0,
        "scalar context entered the AVX2 K9 terminal");

    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh,
            kOriginalCount, recovery_count);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);
    CheckParity(field, generator, original, recovery, recovery_count,
        shard_bytes);
    leo2_codec_destroy(codec);
    leo2_context_destroy(context);
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
                "K9 exact AVX2 terminal test skipped: AVX2 unavailable\n");
            return 0;
        }
        RequireResult(context_result, LEO2_SUCCESS,
            "create AVX2 K9 terminal context");
        for (unsigned recovery_count = kMinimumRecoveryCount;
             recovery_count <= kMaximumRecoveryCount; ++recovery_count)
            ExerciseAVX2Terminal(context, recovery_count, 256);
        ExerciseAVX2Terminal(context, 5, 1024);
        leo2_context_destroy(context);
        for (unsigned recovery_count = kMinimumRecoveryCount;
             recovery_count <= kMaximumRecoveryCount; ++recovery_count)
            ExerciseScalarFallback(recovery_count, 256);
        ExerciseScalarFallback(5, 1024);
        std::printf(
            "K9/R=5..8/256 and K9/R5/1024 packed AVX2 terminal checks passed\n");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "K9 exact terminal failure: %s\n",
            error.what());
        return 1;
    }
}
