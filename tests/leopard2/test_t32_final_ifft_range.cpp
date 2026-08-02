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
#include "Leopard2Backend.h"
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
        , bytes_(bytes == 0 ? 1 : bytes)
    {
#if defined(_MSC_VER)
        value_ = _aligned_malloc(bytes_, leo2_scratch_alignment());
#else
        if (posix_memalign(&value_, leo2_scratch_alignment(), bytes_) != 0)
            value_ = NULL;
#endif
        if (!value_)
            throw std::bad_alloc();
        std::memset(value_, 0, bytes_);
    }

    ~AlignedBuffer()
    {
#if defined(_MSC_VER)
        _aligned_free(value_);
#else
        std::free(value_);
#endif
    }

    void* data() const { return value_; }
    size_t size() const { return bytes_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);
    void* value_;
    size_t bytes_;
};

class Shards
{
public:
    Shards(unsigned count, size_t bytes, unsigned misalignment, uint8_t fill)
        : bytes_(bytes)
        , offset_(8U + misalignment)
        , rows_(count, std::vector<uint8_t>(
            bytes + offset_ + 8U, fill))
    {
    }

    unsigned count() const { return static_cast<unsigned>(rows_.size()); }
    size_t bytes() const { return bytes_; }

    uint8_t* data(unsigned row)
    {
        return &rows_[row][offset_];
    }

    const uint8_t* data(unsigned row) const
    {
        return &rows_[row][offset_];
    }

    std::vector<const void*> const_pointers() const
    {
        std::vector<const void*> result(count());
        for (unsigned row = 0; row < count(); ++row)
            result[row] = data(row);
        return result;
    }

    std::vector<void*> mutable_pointers()
    {
        std::vector<void*> result(count());
        for (unsigned row = 0; row < count(); ++row)
            result[row] = data(row);
        return result;
    }

    void fill_input(uint64_t seed)
    {
        uint64_t state = seed;
        for (unsigned row = 0; row < count(); ++row)
            for (size_t i = 0; i < bytes_; ++i)
            {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                data(row)[i] = static_cast<uint8_t>(state >> 24);
            }
    }

    bool guards_equal(uint8_t fill) const
    {
        for (unsigned row = 0; row < count(); ++row)
        {
            for (size_t i = 0; i < offset_; ++i)
                if (rows_[row][i] != fill)
                    return false;
            for (size_t i = offset_ + bytes_; i < rows_[row].size(); ++i)
                if (rows_[row][i] != fill)
                    return false;
        }
        return true;
    }

    bool equals(const Shards& other) const
    {
        return rows_ == other.rows_;
    }

private:
    size_t bytes_;
    size_t offset_;
    std::vector<std::vector<uint8_t> > rows_;
};

class DisableRangeScope
{
public:
    explicit DisableRangeScope(bool disabled)
    {
        leopard::ff8::TestOnlySetHighFinalIFFT2RangeDisabled(disabled);
    }

    ~DisableRangeScope()
    {
        leopard::ff8::TestOnlySetHighFinalIFFT2RangeDisabled(false);
    }

private:
    DisableRangeScope(const DisableRangeScope&);
    DisableRangeScope& operator=(const DisableRangeScope&);
};

leo2_codec* CreateCodec(
    leo2_context* context,
    unsigned k,
    unsigned r,
    leo2_profile profile)
{
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, k, r, profile, LEO2_FIELD_GF8,
        NULL, &codec), LEO2_SUCCESS, "create T32 test codec");
    return codec;
}

void CheckDirectParity(
    leopard2_test::ProfileKind kind,
    unsigned k,
    unsigned r,
    const Shards& original,
    const Shards& recovery,
    const std::vector<void*>& recovery_pointers)
{
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    struct CacheEntry
    {
        leopard2_test::ProfileKind kind;
        unsigned k;
        unsigned r;
        leopard2_test::Matrix generator;
    };
    static std::vector<CacheEntry> cache;
    const leopard2_test::Matrix* generator = NULL;
    for (size_t i = 0; i < cache.size(); ++i)
    {
        if (cache[i].kind == kind && cache[i].k == k && cache[i].r == r)
        {
            generator = &cache[i].generator;
            break;
        }
    }
    if (!generator)
    {
        CacheEntry entry;
        entry.kind = kind;
        entry.k = k;
        entry.r = r;
        const leopard2_test::ProfileLayout layout =
            leopard2_test::make_profile_layout(kind, k, r);
        entry.generator =
            leopard2_test::direct_systematic_generator(field, layout);
        cache.push_back(entry);
        generator = &cache.back().generator;
    }

    for (unsigned parity = 0; parity < r; ++parity)
    {
        if (!recovery_pointers[parity])
            continue;
        const std::vector<leopard2_test::Element>& row =
            (*generator)[k + parity];
        for (size_t offset = 0; offset < original.bytes(); ++offset)
        {
            leopard2_test::Element expected = 0;
            for (unsigned source = 0; source < k; ++source)
            {
                expected = field.add(expected, field.multiply(
                    row[source], original.data(source)[offset]));
            }
            Require(recovery.data(parity)[offset] == expected,
                "T32 parity differs from independent direct-RS oracle");
        }
    }
}

std::vector<std::vector<uint8_t> > LegacyControl(
    unsigned k,
    unsigned r,
    size_t bytes,
    const Shards& original)
{
    Require((bytes & 63U) == 0, "legacy control requires 64-byte alignment");
    const unsigned work_count = leo_encode_work_count(k, r);
    Require(work_count != 0, "legacy work-count query failed");
    std::vector<const void*> original_pointers = original.const_pointers();
    std::vector<std::vector<uint8_t> > work(
        work_count, std::vector<uint8_t>(bytes, 0));
    std::vector<void*> work_pointers(work_count);
    for (unsigned i = 0; i < work_count; ++i)
        work_pointers[i] = &work[i][0];

    DisableRangeScope disable(true);
    Require(leo_encode(bytes, k, r, work_count,
            &original_pointers[0], &work_pointers[0]) == Leopard_Success,
        "legacy mature-transform control failed");
    work.resize(r);
    return work;
}

void RunHighCase(
    leo2_context* context,
    unsigned k,
    unsigned r,
    size_t bytes,
    bool expect_range,
    unsigned input_misalignment = 0,
    unsigned output_misalignment = 0,
    bool batch = false,
    bool sparse = false,
    bool check_direct = true)
{
    leo2_codec* codec = CreateCodec(
        context, k, r, LEO2_PROFILE_LEGACY_HIGH_V1);
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        LEO2_SUCCESS, "query T32 encode scratch");
    AlignedBuffer scratch(scratch_bytes);
    Shards original(k, bytes, input_misalignment, 0x3c);
    Shards recovery(r, bytes, output_misalignment, 0xa5);
    original.fill_input(UINT64_C(0x74333266696e616c) ^
        static_cast<uint64_t>(k * 257U + r * 17U + bytes));
    const Shards original_before = original;
    std::vector<const void*> original_pointers = original.const_pointers();
    std::vector<void*> recovery_pointers = recovery.mutable_pointers();
    if (sparse)
        for (unsigned parity = 2; parity < r; parity += 5)
            recovery_pointers[parity] = NULL;

    leopard::ff8::TestOnlySetHighFinalIFFT2RangeDisabled(false);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    if (batch)
    {
        leo2_encode_batch_item item = {
            bytes, &original_pointers[0], &recovery_pointers[0],
            scratch.data(), scratch_bytes
        };
        RequireResult(leo2_encode_batch(codec, &item, 1), LEO2_SUCCESS,
            "execute one-item T32 batch");
    }
    else
    {
        RequireResult(leo2_encode(codec, bytes,
            &original_pointers[0], &recovery_pointers[0],
            scratch.data(), scratch_bytes), LEO2_SUCCESS,
            "execute T32 encode");
    }

    const uint64_t expected_calls = expect_range ? (k + 31U) / 32U : 0U;
    const leopard::ff8::TestOnlyHighEncodeCounts counts =
        leopard::ff8::TestOnlyGetHighEncodeCounts();
    Require(counts.final_ifft2_range_calls == expected_calls,
        "T32 final-IFFT range routing count mismatch");
    Require(original.equals(original_before),
        "T32 final-IFFT range modified a source shard or guard");
    Require(original.guards_equal(0x3c) && recovery.guards_equal(0xa5),
        "T32 final-IFFT range modified a guard byte");
    if (check_direct)
        CheckDirectParity(leopard2_test::kLegacyHigh,
            k, r, original, recovery, recovery_pointers);
    if (sparse)
    {
        for (unsigned parity = 2; parity < r; parity += 5)
            for (size_t i = 0; i < bytes; ++i)
                Require(recovery.data(parity)[i] == 0xa5,
                    "sparse T32 encode modified an unrequested parity row");
    }

    if (!sparse && (bytes & 63U) == 0)
    {
        const std::vector<std::vector<uint8_t> > legacy =
            LegacyControl(k, r, bytes, original);
        for (unsigned parity = 0; parity < r; ++parity)
            Require(std::memcmp(recovery.data(parity), &legacy[parity][0],
                    bytes) == 0,
                "T32 parity differs from disabled mature legacy transform");
    }

    leo2_codec_destroy(codec);
}

void TestTargetMatrix(leo2_context* context)
{
    static const unsigned original_counts[] = {
        63, 64, 65, 96, 100, 126, 127, 128
    };
    static const unsigned recovery_counts[] = { 28, 31, 32 };
    for (size_t k_i = 0;
         k_i < sizeof(original_counts) / sizeof(original_counts[0]); ++k_i)
        for (size_t r_i = 0;
             r_i < sizeof(recovery_counts) / sizeof(recovery_counts[0]); ++r_i)
        {
            // The mature legacy transform remains the byte-for-byte control for
            // the full matrix.  Keep the cubic direct-algebra oracle focused on
            // representative boundary and multi-block shapes so this regression
            // test does not become an exhaustive direct-RS campaign.
            const unsigned k = original_counts[k_i];
            const unsigned r = recovery_counts[r_i];
            const bool check_direct =
                (k == 63U && r == 28U) ||
                (k == 64U && r == 32U) ||
                (k == 65U && r == 31U) ||
                (k == 96U && r == 28U);
            RunHighCase(context, original_counts[k_i], recovery_counts[r_i],
                64, true, 0, 0, false, false, check_direct);
        }
}

void TestFallbacks(leo2_context* avx2_context)
{
    // Ragged entries execute one or more padded 64-byte internal passes, while
    // the larger aligned controls may execute 64-byte tiles.  Neither class may
    // inherit the exact-public-B=64 candidate policy.
    static const size_t byte_counts[] = {
        1, 31, 32, 63, 65, 127, 129, 256, 1024
    };
    for (size_t i = 0; i < sizeof(byte_counts) / sizeof(byte_counts[0]); ++i)
        RunHighCase(avx2_context, 64, 32, byte_counts[i], false);

    RunHighCase(avx2_context, 64, 31, 64, true, 1, 3);
    RunHighCase(avx2_context, 64, 31, 64, true, 0, 0, true, false);
    RunHighCase(avx2_context, 64, 31, 64, true, 0, 0, false, true);

    leo2_codec* low = CreateCodec(
        avx2_context, 31, 32, LEO2_PROFILE_LOW_V1);
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(low, 64, &scratch_bytes),
        LEO2_SUCCESS, "query low-profile guard scratch");
    AlignedBuffer scratch(scratch_bytes);
    Shards original(31, 64, 1, 0x3c);
    Shards recovery(32, 64, 3, 0xa5);
    original.fill_input(UINT64_C(0x6c6f772d67756172));
    std::vector<const void*> original_pointers = original.const_pointers();
    std::vector<void*> recovery_pointers = recovery.mutable_pointers();
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(low, 64, &original_pointers[0],
        &recovery_pointers[0], scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute low-profile T32 guard");
    Require(leopard::ff8::TestOnlyGetHighEncodeCounts()
            .final_ifft2_range_calls == 0,
        "low profile entered the high-only T32 range");
    CheckDirectParity(leopard2_test::kLow,
        31, 32, original, recovery, recovery_pointers);
    leo2_codec_destroy(low);

    leo2_context_options options = {};
    options.struct_size = sizeof(options);
    options.backend = LEO2_BACKEND_SCALAR;
    options.thread_count = 1;
    leo2_context* scalar = NULL;
    RequireResult(leo2_context_create(&options, &scalar), LEO2_SUCCESS,
        "create scalar T32 fallback context");
    RunHighCase(scalar, 64, 32, 64, false);
    leo2_context_destroy(scalar);
}

void TestSameBinaryControl(leo2_context* context)
{
    static const unsigned k = 64;
    static const unsigned r = 32;
    static const size_t bytes = 64;
    leo2_codec* codec = CreateCodec(
        context, k, r, LEO2_PROFILE_LEGACY_HIGH_V1);
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        LEO2_SUCCESS, "query same-binary control scratch");
    AlignedBuffer scratch(scratch_bytes);
    Shards original(k, bytes, 0, 0x3c);
    Shards candidate(r, bytes, 0, 0xa5);
    Shards control(r, bytes, 0, 0xa5);
    original.fill_input(UINT64_C(0x73616d652d62696e));
    std::vector<const void*> original_pointers = original.const_pointers();
    std::vector<void*> candidate_pointers = candidate.mutable_pointers();
    std::vector<void*> control_pointers = control.mutable_pointers();

    leopard::ff8::TestOnlySetHighFinalIFFT2RangeDisabled(false);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, bytes, &original_pointers[0],
        &candidate_pointers[0], scratch.data(), scratch_bytes), LEO2_SUCCESS,
        "execute enabled same-binary T32 path");
    Require(leopard::ff8::TestOnlyGetHighEncodeCounts()
            .final_ifft2_range_calls == 2,
        "enabled same-binary T32 path missed the range");

    {
        DisableRangeScope disable(true);
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        RequireResult(leo2_encode(codec, bytes, &original_pointers[0],
            &control_pointers[0], scratch.data(), scratch_bytes), LEO2_SUCCESS,
            "execute disabled same-binary T32 control");
        Require(leopard::ff8::TestOnlyGetHighEncodeCounts()
                .final_ifft2_range_calls == 0,
            "disabled same-binary control entered the T32 range");
    }
    for (unsigned parity = 0; parity < r; ++parity)
        Require(std::memcmp(candidate.data(parity), control.data(parity),
                bytes) == 0,
            "enabled T32 range differs from same-binary mature control");
    leo2_codec_destroy(codec);
}

void TestRejectedAliases(leo2_context* context)
{
    static const unsigned k = 64;
    static const unsigned r = 32;
    static const size_t bytes = 64;
    leo2_codec* codec = CreateCodec(
        context, k, r, LEO2_PROFILE_LEGACY_HIGH_V1);
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
        LEO2_SUCCESS, "query alias-test scratch");
    AlignedBuffer scratch(scratch_bytes);
    Shards original(k, bytes, 0, 0x3c);
    Shards recovery(r, bytes, 0, 0xa5);
    original.fill_input(UINT64_C(0x616c6961732d7433));
    const Shards original_before = original;
    const Shards recovery_before = recovery;
    std::vector<const void*> original_pointers = original.const_pointers();
    std::vector<void*> recovery_pointers = recovery.mutable_pointers();

    recovery_pointers[0] = const_cast<void*>(original_pointers[0]);
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, bytes, &original_pointers[0],
        &recovery_pointers[0], scratch.data(), scratch_bytes), LEO2_OVERLAP,
        "reject T32 recovery/source overlap");
    Require(leopard::ff8::TestOnlyGetHighEncodeCounts()
            .final_ifft2_range_calls == 0,
        "source overlap reached the T32 arithmetic range");
    Require(original.equals(original_before) && recovery.equals(recovery_before),
        "source-overlap rejection modified T32 storage");

    recovery_pointers = recovery.mutable_pointers();
    recovery_pointers[1] = recovery_pointers[0];
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, bytes, &original_pointers[0],
        &recovery_pointers[0], scratch.data(), scratch_bytes), LEO2_OVERLAP,
        "reject duplicate T32 recovery outputs");
    Require(leopard::ff8::TestOnlyGetHighEncodeCounts()
            .final_ifft2_range_calls == 0,
        "duplicate outputs reached the T32 arithmetic range");

    recovery_pointers = recovery.mutable_pointers();
    leopard::ff8::TestOnlyResetHighEncodeCounts();
    RequireResult(leo2_encode(codec, bytes, &original_pointers[0],
        &recovery_pointers[0], scratch.data(), scratch_bytes - 1),
        LEO2_SCRATCH_TOO_SMALL, "reject undersized T32 scratch");
    Require(leopard::ff8::TestOnlyGetHighEncodeCounts()
            .final_ifft2_range_calls == 0,
        "undersized scratch reached the T32 arithmetic range");
    Require(original.equals(original_before) && recovery.equals(recovery_before),
        "rejected T32 call modified storage");
    leo2_codec_destroy(codec);
}

void CheckBackendCapabilities()
{
    const leopard::backend::Ops* avx2 =
        leopard::backend::GetQualifiedOps(LEO2_BACKEND_AVX2);
    Require(avx2 && avx2->ff8_ifft_butterfly2_range,
        "qualified pure AVX2 backend omitted the T32 range callback");
    const leopard::backend::Ops* scalar =
        leopard::backend::GetQualifiedOps(LEO2_BACKEND_SCALAR);
    Require(scalar && !scalar->ff8_ifft_butterfly2_range,
        "scalar backend exposed the AVX2-only T32 range callback");
    const leopard::backend::Ops* ssse3 =
        leopard::backend::GetQualifiedOps(LEO2_BACKEND_SSSE3);
    if (ssse3)
        Require(!ssse3->ff8_ifft_butterfly2_range,
            "SSSE3 backend exposed the AVX2-only T32 range callback");
    const leopard::backend::Ops* avx512 =
        leopard::backend::GetQualifiedOps(LEO2_BACKEND_AVX512);
    if (avx512)
        Require(!avx512->ff8_ifft_butterfly2_range,
            "AVX-512 backend exposed the pure-AVX2 T32 range callback");
    const leopard::backend::Ops* gfni =
        leopard::backend::GetQualifiedOps(LEO2_BACKEND_GFNI);
    if (gfni)
        Require(!gfni->ff8_ifft_butterfly2_range,
            "GFNI backend exposed the nibble-table T32 range callback");
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
            std::printf("T32 final-IFFT range test skipped: AVX2 unavailable\n");
            return 0;
        }
        RequireResult(result, LEO2_SUCCESS, "create AVX2 T32 test context");
        CheckBackendCapabilities();
        TestTargetMatrix(context);
        TestFallbacks(context);
        TestSameBinaryControl(context);
        TestRejectedAliases(context);
        leopard::ff8::TestOnlySetHighFinalIFFT2RangeDisabled(false);
        leo2_context_destroy(context);
        std::printf("T32 final-IFFT range checks passed\n");
        return 0;
    }
    catch (const std::exception& error)
    {
        leopard::ff8::TestOnlySetHighFinalIFFT2RangeDisabled(false);
        std::fprintf(stderr, "T32 final-IFFT range failure: %s\n",
            error.what());
        return 1;
    }
}
