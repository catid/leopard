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

static const unsigned kRecoveryCount = 2;
static const size_t kK4TerminalMaximumBytes = 16U * 1024U * 1024U;

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

void FillInput(uint8_t* input, unsigned k, size_t bytes, uint64_t seed)
{
    uint64_t state = seed;
    for (size_t i = 0; i < static_cast<size_t>(k) * bytes; ++i)
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
    unsigned k,
    size_t bytes,
    const void** original,
    void** recovery)
{
    for (unsigned i = 0; i < k; ++i)
        original[i] = input + static_cast<size_t>(i) * bytes;
    for (unsigned i = 0; i < kRecoveryCount; ++i)
        recovery[i] = output + static_cast<size_t>(i) * bytes;
}

void CheckParity(
    unsigned k,
    size_t bytes,
    const void* const* original,
    void* const* recovery)
{
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            leopard2_test::kLegacyHigh, k, kRecoveryCount);
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(field, layout);
    for (unsigned parity = 0; parity < kRecoveryCount; ++parity)
    {
        if (!recovery[parity])
            continue;
        uint8_t* const output = static_cast<uint8_t*>(recovery[parity]);
        const std::vector<leopard2_test::Element>& row =
            generator[k + parity];
        for (size_t offset = 0; offset < bytes; ++offset)
        {
            leopard2_test::Element expected = 0;
            for (unsigned source = 0; source < k; ++source)
            {
                const uint8_t value =
                    static_cast<const uint8_t*>(original[source])[offset];
                expected = field.add(expected,
                    field.multiply(row[source], value));
            }
            Require(output[offset] == expected,
                "T=2 parity differs from the independent generator oracle");
        }
    }
}

uint64_t PackedCalls()
{
    return leopard::ff8::TestOnlyGetHighEncodeCounts().t2_packed_calls;
}

void ResetOutput(uint8_t* output, size_t bytes)
{
    std::memset(output, 0xa5, kRecoveryCount * bytes);
}

leo2_codec* CreateCodec(leo2_context* context, unsigned k)
{
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, k, kRecoveryCount,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        LEO2_SUCCESS, "create T=2 codec");
    return codec;
}

void ExercisePackedCell(
    leo2_context* context,
    unsigned k,
    size_t bytes,
    bool expect_terminal)
{
    leo2_codec* codec = CreateCodec(context, k);
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        LEO2_SUCCESS, "query T=2 scratch");
    AlignedBuffer scratch(scratch_bytes + leo2_scratch_alignment());
    AlignedBuffer input(static_cast<size_t>(k) * bytes + 8);
    AlignedBuffer output(kRecoveryCount * bytes + 8);
    FillInput(input.bytes() + 1, k, bytes,
        UINT64_C(0x54325041434b4544) ^ (static_cast<uint64_t>(k) << 32) ^ bytes);
    const void* original[16] = {};
    void* recovery[kRecoveryCount] = { NULL, NULL };
    SetPackedPointers(
        input.bytes() + 1, output.bytes() + 3, k, bytes,
        original, recovery);

    ResetOutput(output.bytes() + 3, bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, bytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute packed T=2 public encode");
    Require(PackedCalls() == (expect_terminal ? 1U : 0U),
        "packed public T=2 terminal selection mismatch");
    CheckParity(k, bytes, original, recovery);

    ResetOutput(output.bytes() + 3, bytes);
    leo2_encode_batch_item item = {
        bytes, original, recovery, scratch.data(), scratch_bytes
    };
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch(codec, &item, 1), LEO2_SUCCESS,
        "execute packed T=2 one-item batch");
    Require(PackedCalls() == (expect_terminal ? 1U : 0U),
        "packed batch T=2 terminal selection mismatch");
    CheckParity(k, bytes, original, recovery);

    if (k == 4 && (bytes == 64 || bytes == 65))
    {
        /* Exhaust every byte value in every source coordinate.  B=64 proves
           the AVX2 table maps, while B=65 places the value in the overlapping
           four-byte tail and proves the identical legacy representation. */
        const size_t offset = bytes == 64 ? 0U : 64U;
        for (unsigned source = 0; source < k; ++source)
        {
            for (unsigned value = 0; value < 256; ++value)
            {
                std::memset(input.bytes() + 1, 0,
                    static_cast<size_t>(k) * bytes);
                input.bytes()[1 + static_cast<size_t>(source) * bytes +
                    offset] = static_cast<uint8_t>(value);
                ResetOutput(output.bytes() + 3, bytes);
                leopard::ff8::TestOnlyResetHighEncodeCounts();
                RequireResult(leo2_encode(codec, bytes, original, recovery,
                    scratch.data(), scratch_bytes), LEO2_SUCCESS,
                    "execute exhaustive K4/R2 source value");
                Require(PackedCalls() == 1,
                    "K4/R2 source value missed the packed terminal");
                CheckParity(k, bytes, original, recovery);
            }
        }
    }

    if (expect_terminal)
    {
        AlignedBuffer second_output(kRecoveryCount * bytes + 16);
        AlignedBuffer second_scratch(
            scratch_bytes + leo2_scratch_alignment());
        std::vector<uint8_t> detached_source(bytes);
        std::memcpy(&detached_source[0], original[k - 1], bytes);
        for (size_t i = 0; i < bytes; ++i)
        {
            detached_source[i] ^= static_cast<uint8_t>(
                0x5bU + static_cast<unsigned>(i * 29U));
        }
        const void* second_original[16] = {};
        for (unsigned i = 0; i < k; ++i)
            second_original[i] = original[i];
        second_original[k - 1] = &detached_source[0];
        void* second_recovery[kRecoveryCount] = {
            second_output.bytes() + 2,
            second_output.bytes() + 2 + bytes + 5
        };
        ResetOutput(output.bytes() + 3, bytes);
        std::memset(second_recovery[0], 0xa5, bytes);
        std::memset(second_recovery[1], 0xa5, bytes);
        leo2_encode_batch_item items[2] = {
            { bytes, original, recovery, scratch.data(), scratch_bytes },
            { bytes, second_original, second_recovery,
                second_scratch.data(), scratch_bytes }
        };
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode_batch(codec, items, 2), LEO2_SUCCESS,
            "execute packed T=2 two-item batch");
        Require(PackedCalls() == 2,
            "prevalidated two-item T=2 routing mismatch");
        CheckParity(k, bytes, original, recovery);
        CheckParity(k, bytes, second_original, second_recovery);

        if (k == 4 && bytes == 64)
        {
            /* Multi-item preflight must reject cross-stripe writable aliasing
               before the prevalidated executor can write either stripe. */
            void* const independent_output = second_recovery[0];
            second_recovery[0] = recovery[0];
            ResetOutput(output.bytes() + 3, bytes);
            std::memset(independent_output, 0xa5, bytes);
            const std::vector<uint8_t> first_before(
                output.bytes() + 3,
                output.bytes() + 3 + kRecoveryCount * bytes);
            leopard::ff8::TestOnlyResetHighEncodeCounts();
            RequireResult(leo2_encode_batch(codec, items, 2),
                LEO2_OVERLAP,
                "reject cross-item K4/R2 writable alias");
            Require(PackedCalls() == 0,
                "cross-item alias reached the K4/R2 executor");
            Require(std::memcmp(output.bytes() + 3, &first_before[0],
                    first_before.size()) == 0,
                "cross-item alias rejection modified first output");
            second_recovery[0] = independent_output;
        }

        leo2_encode_batch_binding* two_item_binding = NULL;
        RequireResult(leo2_encode_batch_binding_create(
            codec, items, 2, &two_item_binding), LEO2_SUCCESS,
            "create two-item T=2 encode binding");
        ResetOutput(output.bytes() + 3, bytes);
        std::memset(second_recovery[0], 0xa5, bytes);
        std::memset(second_recovery[1], 0xa5, bytes);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode_batch_binding_execute(two_item_binding),
            LEO2_SUCCESS, "execute two-item T=2 encode binding");
        Require(PackedCalls() == 2,
            "prevalidated two-item binding T=2 routing mismatch");
        CheckParity(k, bytes, original, recovery);
        CheckParity(k, bytes, second_original, second_recovery);
        leo2_encode_batch_binding_destroy(two_item_binding);

        void* const second_parity_one = second_recovery[1];
        second_recovery[1] = NULL;
        ResetOutput(output.bytes() + 3, bytes);
        std::memset(second_recovery[0], 0xa5, bytes);
        std::memset(second_parity_one, 0xa5, bytes);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode_batch(codec, items, 2), LEO2_SUCCESS,
            "execute mixed full/sparse T=2 two-item batch");
        Require(PackedCalls() == 1,
            "mixed sparse T=2 batch routing mismatch");
        CheckParity(k, bytes, original, recovery);
        CheckParity(k, bytes, second_original, second_recovery);
        for (size_t i = 0; i < bytes; ++i)
            Require(static_cast<const uint8_t*>(second_parity_one)[i] == 0xa5,
                "mixed sparse T=2 batch modified null parity");

        leo2_encode_batch_binding* sparse_binding = NULL;
        RequireResult(leo2_encode_batch_binding_create(
            codec, items, 2, &sparse_binding), LEO2_SUCCESS,
            "create mixed sparse T=2 encode binding");
        ResetOutput(output.bytes() + 3, bytes);
        std::memset(second_recovery[0], 0xa5, bytes);
        std::memset(second_parity_one, 0xa5, bytes);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode_batch_binding_execute(sparse_binding),
            LEO2_SUCCESS, "execute mixed sparse T=2 encode binding");
        Require(PackedCalls() == 1,
            "mixed sparse binding T=2 routing mismatch");
        CheckParity(k, bytes, original, recovery);
        CheckParity(k, bytes, second_original, second_recovery);
        for (size_t i = 0; i < bytes; ++i)
            Require(static_cast<const uint8_t*>(second_parity_one)[i] == 0xa5,
                "mixed sparse T=2 binding modified null parity");
        leo2_encode_batch_binding_destroy(sparse_binding);
        second_recovery[1] = second_parity_one;
    }

    ResetOutput(output.bytes() + 3, bytes);
    leo2_encode_batch_binding* binding = NULL;
    RequireResult(leo2_encode_batch_binding_create(
        codec, &item, 1, &binding), LEO2_SUCCESS,
        "create T=2 encode binding");
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch_binding_execute(binding), LEO2_SUCCESS,
        "execute T=2 encode binding");
    Require(PackedCalls() == (expect_terminal ? 1U : 0U),
        "prevalidated T=2 terminal selection mismatch");
    CheckParity(k, bytes, original, recovery);
    leo2_encode_batch_binding_destroy(binding);

    ResetOutput(output.bytes() + 3, bytes);
    RequireResult(leo2_test_codec_set_encode_mode(
        codec, LEO2_TEST_ENCODE_FORCE_TRANSFORM), LEO2_SUCCESS,
        "force T=2 transform oracle");
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, bytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute forced T=2 transform oracle");
    Require(PackedCalls() == 0,
        "forced transform entered the packed T=2 terminal");
    CheckParity(k, bytes, original, recovery);

    leo2_codec_destroy(codec);
}

void ExerciseFallbackAndErrors(leo2_context* context, unsigned k)
{
    const size_t bytes = 64;
    leo2_codec* codec = CreateCodec(context, k);
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        LEO2_SUCCESS, "query T=2 contract scratch");
    AlignedBuffer scratch(scratch_bytes + leo2_scratch_alignment());
    AlignedBuffer input(k * bytes);
    AlignedBuffer output(kRecoveryCount * bytes);
    FillInput(input.bytes(), k, bytes, UINT64_C(0x5432434f4e545241));
    const void* original[16] = {};
    void* recovery[kRecoveryCount] = { NULL, NULL };
    SetPackedPointers(
        input.bytes(), output.bytes(), k, bytes, original, recovery);

    /* Sparse output and non-packed input remain valid mature fallbacks. */
    ResetOutput(output.bytes(), bytes);
    recovery[1] = NULL;
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, bytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute sparse T=2 fallback");
    Require(PackedCalls() == 0,
        "sparse output entered the packed T=2 terminal");
    CheckParity(k, bytes, original, recovery);
    for (size_t i = 0; i < bytes; ++i)
        Require(output.bytes()[bytes + i] == 0xa5,
            "sparse T=2 fallback modified null parity");

    SetPackedPointers(
        input.bytes(), output.bytes(), k, bytes, original, recovery);
    ResetOutput(output.bytes(), bytes);
    recovery[0] = NULL;
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, bytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute parity-one-only T=2 fallback");
    Require(PackedCalls() == 0,
        "parity-one-only output entered the packed T=2 terminal");
    CheckParity(k, bytes, original, recovery);
    for (size_t i = 0; i < bytes; ++i)
        Require(output.bytes()[i] == 0xa5,
            "parity-one-only T=2 fallback modified null parity");

    recovery[1] = NULL;
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, bytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute all-null T=2 no-op");
    Require(PackedCalls() == 0,
        "all-null T=2 no-op entered the packed circuit");

    SetPackedPointers(
        input.bytes(), output.bytes(), k, bytes, original, recovery);
    std::vector<uint8_t> detached(bytes);
    std::memcpy(&detached[0], original[1], bytes);
    original[1] = &detached[0];
    ResetOutput(output.bytes(), bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, bytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute non-packed T=2 fallback");
    Require(PackedCalls() == 0,
        "non-packed input entered the packed T=2 terminal");
    CheckParity(k, bytes, original, recovery);

    SetPackedPointers(
        input.bytes(), output.bytes(), k, bytes, original, recovery);
    ResetOutput(output.bytes(), bytes);
    const std::vector<uint8_t> output_before(
        output.bytes(), output.bytes() + kRecoveryCount * bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, bytes, original, recovery,
        scratch.data(), scratch_bytes - 1), LEO2_SCRATCH_TOO_SMALL,
        "reject undersized T=2 terminal scratch");
    Require(PackedCalls() == 0,
        "undersized scratch reached the packed T=2 circuit");
    Require(std::memcmp(output.bytes(), &output_before[0],
            output_before.size()) == 0,
        "undersized scratch modified T=2 output");

    RequireResult(leo2_encode(codec, bytes, original, recovery,
        scratch.bytes() + 1, scratch_bytes), LEO2_BAD_ALIGNMENT,
        "reject misaligned T=2 terminal scratch");
    Require(PackedCalls() == 0,
        "misaligned scratch reached the packed T=2 circuit");

    for (unsigned i = 0; i < kRecoveryCount; ++i)
        recovery[i] = input.bytes() + static_cast<size_t>(i) * bytes;
    const std::vector<uint8_t> input_before(
        input.bytes(), input.bytes() + k * bytes);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, bytes, original, recovery,
        scratch.data(), scratch_bytes), LEO2_OVERLAP,
        "reject overlapping packed T=2 slabs");
    Require(PackedCalls() == 0,
        "overlapping T=2 slabs reached the packed circuit");
    Require(std::memcmp(input.bytes(), &input_before[0],
            input_before.size()) == 0,
        "T=2 overlap rejection modified input");

    SetPackedPointers(
        input.bytes(), output.bytes(), k, bytes, original, recovery);
    alignas(64) uint8_t protected_storage[kRecoveryCount * bytes] = {};
    leo2_encode_batch_item* const protected_item =
        new (protected_storage) leo2_encode_batch_item;
    protected_item->shard_bytes = bytes;
    protected_item->original = original;
    protected_item->recovery = recovery;
    protected_item->scratch = scratch.data();
    protected_item->scratch_bytes = scratch_bytes;
    for (unsigned i = 0; i < kRecoveryCount; ++i)
        recovery[i] = protected_storage + static_cast<size_t>(i) * bytes;
    const leo2_encode_batch_item protected_before = *protected_item;
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode_batch(codec, protected_item, 1),
        LEO2_OVERLAP, "reject T=2 output/batch-metadata overlap");
    Require(PackedCalls() == 0,
        "metadata overlap reached the packed T=2 circuit");
    Require(std::memcmp(protected_item, &protected_before,
            sizeof(*protected_item)) == 0,
        "metadata-overlap rejection modified the batch descriptor");

    leo2_codec_destroy(codec);
}

void ExerciseNonAVX2Fallback(
    leo2_backend backend,
    const char* backend_name,
    unsigned k)
{
    leo2_context_options options = {};
    options.struct_size = sizeof(options);
    options.backend = backend;
    options.thread_count = 1;
    leo2_context* context = NULL;
    const leo2_result context_result = leo2_context_create(&options, &context);
    if (context_result == LEO2_UNSUPPORTED)
        return;
    RequireResult(context_result, LEO2_SUCCESS,
        "create non-AVX2 T=2 context");
    leo2_codec* codec = CreateCodec(context, k);
    const size_t bytes = 64;
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        LEO2_SUCCESS, "query scalar T=2 scratch");
    AlignedBuffer scratch(scratch_bytes);
    AlignedBuffer input(k * bytes);
    AlignedBuffer output(kRecoveryCount * bytes);
    FillInput(input.bytes(), k, bytes,
        UINT64_C(0x5343414c41525432) ^ k);
    const void* original[16] = {};
    void* recovery[kRecoveryCount];
    SetPackedPointers(
        input.bytes(), output.bytes(), k, bytes, original, recovery);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, bytes, original, recovery,
        scratch.data(), scratch.size()), LEO2_SUCCESS,
        "execute non-AVX2 T=2 fallback");
    Require(PackedCalls() == 0,
        "non-AVX2 context entered the AVX2 T=2 terminal");
    CheckParity(k, bytes, original, recovery);
    leo2_codec_destroy(codec);
    leo2_context_destroy(context);
    std::printf("%s K=%u T=2 fallback matched the generator oracle\n",
        backend_name, k);
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
            std::printf("T=2 packed AVX2 terminal test skipped: AVX2 unavailable\n");
            return 0;
        }
        RequireResult(result, LEO2_SUCCESS, "create AVX2 T=2 context");

        /* K=4 owns the generated circuit whose 16..31-byte remainder now
           uses an overlapping in-range vector.  Exhaust every exact length
           through the first two vector boundaries across public, batch, and
           reusable-binding entry points. */
        for (size_t bytes = 1; bytes <= 65; ++bytes)
            ExercisePackedCell(context, 4, bytes, true);

        static const size_t sizes[] = {
            1, 2, 3, 7, 15, 16, 17, 31, 32, 33, 63, 64, 65,
            127, 128, 129, 255, 256, 257, 1023, 1024, 1025,
            1984, 2016, 2047, 2048, 2049, 2111, 2112, 2113,
            4095, 4096, 4097, 8192, 16384, 16385,
            32768, 65536, 65537
        };
        for (unsigned k = 2; k <= 4; ++k)
        {
            for (size_t i = 0; i < sizeof(sizes) / sizeof(sizes[0]); ++i)
            {
                const size_t bytes = sizes[i];
                if (k == 4 && bytes <= 65)
                    continue;
                const bool selected = bytes <= 2048 ||
                    (k == 4 && bytes <= kK4TerminalMaximumBytes);
                ExercisePackedCell(context, k, bytes, selected);
            }
        }
        static const size_t multi_sizes[] = {
            32, 63, 64, 65, 128, 1024, 1984, 2047, 2048, 2049
        };
        for (unsigned k = 5; k <= 16; ++k)
        {
            for (size_t i = 0;
                 i < sizeof(multi_sizes) / sizeof(multi_sizes[0]); ++i)
            {
                const size_t bytes = multi_sizes[i];
                const bool selected = bytes <= 2048 && (bytes & 63U) == 0;
                ExercisePackedCell(context, k, bytes, selected);
            }
        }
        ExerciseFallbackAndErrors(context, 2);
        ExerciseFallbackAndErrors(context, 4);
        ExerciseFallbackAndErrors(context, 16);
        leo2_context_destroy(context);
        for (unsigned k = 2; k <= 4; k += 2)
        {
            ExerciseNonAVX2Fallback(LEO2_BACKEND_SCALAR, "scalar", k);
            ExerciseNonAVX2Fallback(LEO2_BACKEND_SSSE3, "SSSE3", k);
        }
        ExerciseNonAVX2Fallback(LEO2_BACKEND_AVX512, "AVX-512", 4);
        ExerciseNonAVX2Fallback(LEO2_BACKEND_GFNI, "GFNI", 4);
        std::printf("T=2 packed AVX2 terminal checks passed\n");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "T=2 packed terminal failure: %s\n", error.what());
        return 1;
    }
}
