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
#include "leopard2.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

#if !defined(LEO2_EXPECT_HIGH_SPARSE_DIRECT)
#error "production high sparse-direct expectation must be explicit"
#endif
#if !defined(LEO2_EXPECT_HIGH_SPARSE_DIRECT_AUTO)
#error "production high sparse-direct AUTO expectation must be explicit"
#endif
#if !defined(LEO2_EXPECT_HIGH_DIRECT_PRODUCTION)
#error "production high direct-table expectation must be explicit"
#endif
#if defined(LEO2_ENABLE_TEST_HOOKS)
#error "the sparse-high production policy test must use the ordinary archive"
#endif

namespace {

typedef std::vector<uint8_t> Bytes;
typedef std::vector<Bytes> Shards;

struct Cell
{
    unsigned k;
    unsigned r;
    size_t bytes;
};

void Require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void RequireResult(leo2_result result, const char* message)
{
    if (result != LEO2_SUCCESS)
        throw std::runtime_error(
            std::string(message) + ": " + leo2_result_string(result));
}

class AlignedBuffer
{
public:
    explicit AlignedBuffer(size_t bytes)
        : value_(NULL)
        , bytes_(bytes)
    {
        Require(bytes != 0, "zero-size scratch allocation");
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

    void* data() const { return value_; }
    size_t size() const { return bytes_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);
    void* value_;
    size_t bytes_;
};

Shards MakeOriginal(unsigned k, size_t bytes, uint64_t seed)
{
    Shards original(k, Bytes(bytes, 0));
    uint64_t state = seed;
    for (unsigned shard = 0; shard < k; ++shard)
        for (size_t offset = 0; offset < bytes; ++offset)
        {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            original[shard][offset] = static_cast<uint8_t>(
                state + shard * 29U + offset * 131U);
        }
    return original;
}

Bytes OracleParity(
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator,
    const Shards& original,
    unsigned recovery_index)
{
    Bytes parity(original[0].size(), 0);
    const std::vector<leopard2_test::Element>& row =
        generator[original.size() + recovery_index];
    for (size_t offset = 0; offset < parity.size(); ++offset)
    {
        leopard2_test::Element value = 0;
        for (size_t source = 0; source < original.size(); ++source)
        {
            value = field.add(value, field.multiply(
                row[source], original[source][offset]));
        }
        parity[offset] = static_cast<uint8_t>(value);
    }
    return parity;
}

leo2_context* CreateContext(
    leo2_backend backend,
    uint32_t thread_count,
    bool prepare_tables,
    bool auto_select)
{
    leo2_context_options options = {};
    options.struct_size = sizeof(options);
    options.backend = backend;
    options.thread_count = thread_count;
    leo2_context* context = NULL;
    RequireResult(leo2_context_create(&options, &context),
        "sparse-high context create");
    Require(context != NULL, "context creation returned null");
    Require(
        leopard2_internal::SetContextHighSparseDirectEncodePolicyForDiagnostics(
            context, prepare_tables, auto_select),
        "sparse-high context policy rejected a valid state");
    return context;
}

leo2_codec* CreateCodec(
    leo2_context* context,
    unsigned k,
    unsigned r,
    uint32_t flags = 0)
{
    leo2_codec_options options = {};
    options.struct_size = sizeof(options);
    options.flags = flags;
    options.shard_layout = LEO2_SHARD_LAYOUT_NATIVE_V1;
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        &options, &codec), "sparse-high codec create");
    Require(codec != NULL, "codec creation returned null");
    return codec;
}

void CheckPath(
    const leo2_codec* codec,
    size_t bytes,
    unsigned requested_count,
    size_t expected_rows,
    bool expected_direct,
    const char* message)
{
    leopard2_internal::CodecEncodePathInfo path = {};
    Require(leopard2_internal::GetCodecEncodePathInfo(
        codec, bytes, requested_count, &path), message);
    Require(path.direct_generator_rows == expected_rows,
        "sparse-high table preparation differed from policy");
    Require(path.auto_direct_selected == expected_direct,
        "sparse-high AUTO selection differed from policy");
}

void EncodeOneAndCheck(
    const leo2_codec* codec,
    const Shards& original,
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator,
    unsigned r,
    unsigned parity)
{
    const size_t bytes = original[0].size();
    std::vector<const void*> input(original.size(), NULL);
    for (size_t i = 0; i < original.size(); ++i)
        input[i] = original[i].data();
    Shards recovery(r, Bytes(bytes, 0xa5));
    std::vector<void*> output(r, NULL);
    output[parity] = recovery[parity].data();
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, bytes, &scratch_bytes), "sparse-high scratch query");
    AlignedBuffer scratch(scratch_bytes);
    RequireResult(leo2_encode(codec, bytes, input.data(), output.data(),
        scratch.data(), scratch.size()), "sparse-high one-shot encode");
    Require(recovery[parity] == OracleParity(
        field, generator, original, parity),
        "sparse-high one-shot parity differs from oracle");
    for (unsigned i = 0; i < r; ++i)
    {
        if (i != parity)
            Require(recovery[i] == Bytes(bytes, 0xa5),
                "sparse-high one-shot modified a null output");
    }
}

void EncodeTwoAndCheck(
    const leo2_codec* codec,
    const Shards& original,
    const leopard2_test::BinaryField& field,
    const leopard2_test::Matrix& generator,
    unsigned r)
{
    const size_t bytes = original[0].size();
    std::vector<const void*> input(original.size(), NULL);
    for (size_t i = 0; i < original.size(); ++i)
        input[i] = original[i].data();
    Shards recovery(r, Bytes(bytes, 0xa5));
    std::vector<void*> output(r, NULL);
    output[0] = recovery[0].data();
    output[r - 1] = recovery[r - 1].data();
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, bytes, &scratch_bytes), "sparse-high Q2 scratch query");
    AlignedBuffer scratch(scratch_bytes);
    RequireResult(leo2_encode(codec, bytes, input.data(), output.data(),
        scratch.data(), scratch.size()), "sparse-high Q2 encode");
    Require(recovery[0] == OracleParity(field, generator, original, 0) &&
            recovery[r - 1] ==
                OracleParity(field, generator, original, r - 1),
        "sparse-high Q2 parity differs from oracle");
    for (unsigned i = 1; i + 1 < r; ++i)
        Require(recovery[i] == Bytes(bytes, 0xa5),
            "sparse-high Q2 modified a null output");
}

void CheckWitness(
    leo2_context* context,
    uint64_t expected_direct,
    uint64_t expected_transform)
{
    leopard2_internal::HighSparseEncodeRouteWitness witness = {};
    Require(
        leopard2_internal::ReadAndDisarmContextHighSparseEncodeRouteWitnessForDiagnostics(
            context, &witness),
        "sparse-high route witness could not be read");
    Require(witness.direct_calls == expected_direct &&
            witness.transform_calls == expected_transform,
        "sparse-high route witness observed the wrong executor");
}

void ExercisePolicyState(
    bool prepare_tables,
    bool auto_select,
    bool expected_direct)
{
    leo2_context* context = CreateContext(
        LEO2_BACKEND_AUTO, 1, prepare_tables, auto_select);
    Require(leo2_context_backend(context) == LEO2_BACKEND_AVX2,
        "AUTO did not resolve to the AVX2 production baseline");
    leo2_codec* codec = CreateCodec(context, 2, 16);
    CheckPath(codec, 4096, 1,
        (prepare_tables || LEO2_EXPECT_HIGH_DIRECT_PRODUCTION) ? 16 : 0,
        expected_direct, "sparse-high state path query");

    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(
            field, leopard2_test::make_profile_layout(
                leopard2_test::kLegacyHigh, 2, 16));
    const Shards original = MakeOriginal(
        2, 4096, UINT64_C(0x5354415445513141));
    Require(
        leopard2_internal::ArmContextHighSparseEncodeRouteWitnessForDiagnostics(context),
        "sparse-high route witness could not be armed");
    EncodeOneAndCheck(codec, original, field, generator, 16, 7);
    CheckWitness(context, expected_direct ? 1 : 0,
        expected_direct ? 0 : 1);

    leo2_codec_destroy(codec);
    leo2_context_destroy(context);
}

void ExerciseTransformWitness(
    leo2_context* context,
    leo2_codec* codec,
    unsigned k,
    unsigned r,
    size_t bytes,
    uint64_t seed)
{
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(
            field, leopard2_test::make_profile_layout(
                leopard2_test::kLegacyHigh, k, r));
    const Shards original = MakeOriginal(k, bytes, seed);
    Require(
        leopard2_internal::ArmContextHighSparseEncodeRouteWitnessForDiagnostics(context),
        "sparse-high exclusion witness could not be armed");
    EncodeOneAndCheck(codec, original, field, generator, r, r - 1);
    CheckWitness(context, 0, 1);
}

std::vector<Cell> CandidateCells()
{
    static const Cell cells[] = {
        { 2, 2, 4096 }, { 2, 4, 4096 },
        { 2, 8, 4096 }, { 2, 16, 4096 },
        { 3, 2, 4096 }, { 3, 4, 4096 },
        { 3, 8, 4096 }, { 3, 16, 4096 },
        { 4, 2, 4096 }, { 4, 4, 4096 },
        { 4, 8, 4096 }, { 4, 16, 4096 },
        { 8, 2, 4096 }, { 8, 4, 4096 },
        { 8, 8, 4096 }, { 8, 16, 4096 },
        { 12, 2, 4096 }, { 12, 4, 4096 },
        { 12, 8, 4096 }, { 12, 16, 4096 },
        { 16, 2, 4096 }, { 16, 4, 4096 },
        { 16, 8, 4096 }, { 16, 16, 4096 },
        { 2, 16, 1024 }, { 2, 16, 1088 },
        { 2, 16, 2048 }, { 2, 16, 4032 },
        { 2, 16, 4160 }, { 2, 16, 65536 },
        { 16, 2, 1024 }, { 16, 2, 1088 },
        { 16, 2, 2048 }, { 16, 2, 4032 },
        { 16, 2, 4160 }, { 16, 2, 65536 }
    };
    std::vector<Cell> result(cells, cells + sizeof(cells) / sizeof(cells[0]));
    for (size_t i = 0; i < result.size(); ++i)
        for (size_t j = i + 1; j < result.size(); ++j)
            Require(result[i].k != result[j].k ||
                    result[i].r != result[j].r ||
                    result[i].bytes != result[j].bytes,
                "sparse-high candidate table contains a duplicate");
    return result;
}

uint64_t ExerciseBatchEntry(
    leo2_context* context,
    unsigned k,
    unsigned r,
    size_t bytes,
    size_t item_count,
    bool use_binding);

uint64_t ExerciseCandidateCells(
    leo2_context* context,
    uint64_t& route_checks_out)
{
    const std::vector<Cell> cells = CandidateCells();
    Require(cells.size() == 36, "sparse-high candidate cell count changed");
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    uint64_t parity_checks = 0;
    for (size_t cell_index = 0; cell_index < cells.size(); ++cell_index)
    {
        const Cell& cell = cells[cell_index];
        leo2_codec* codec = CreateCodec(context, cell.k, cell.r);
        CheckPath(codec, cell.bytes, 1, cell.r, true,
            "sparse-high candidate path query");
        const leopard2_test::Matrix generator =
            leopard2_test::direct_systematic_generator(
                field, leopard2_test::make_profile_layout(
                    leopard2_test::kLegacyHigh, cell.k, cell.r));
        const Shards original = MakeOriginal(cell.k, cell.bytes,
            UINT64_C(0x43414e4449444154) + cell_index);
        Require(
            leopard2_internal::ArmContextHighSparseEncodeRouteWitnessForDiagnostics(
                context),
            "sparse-high cell witness could not be armed");
        for (unsigned parity = 0; parity < cell.r; ++parity)
        {
            EncodeOneAndCheck(
                codec, original, field, generator, cell.r, parity);
            ++parity_checks;
        }
        CheckWitness(context, cell.r, 0);
        leo2_codec_destroy(codec);
        route_checks_out += cell.r;
        route_checks_out += ExerciseBatchEntry(
            context, cell.k, cell.r, cell.bytes, cell.r, false);
        route_checks_out += ExerciseBatchEntry(
            context, cell.k, cell.r, cell.bytes, cell.r, true);
    }
    return parity_checks;
}

struct BatchStorage
{
    BatchStorage(
        unsigned k,
        unsigned r,
        size_t bytes,
        unsigned parity,
        size_t scratch_bytes,
        uint64_t seed)
        : original(MakeOriginal(k, bytes, seed))
        , recovery(r, Bytes(bytes, 0xa5))
        , input(k, NULL)
        , output(r, NULL)
        , scratch(new AlignedBuffer(scratch_bytes))
        , requested_parity(parity)
    {
        for (unsigned i = 0; i < k; ++i)
            input[i] = original[i].data();
        output[parity] = recovery[parity].data();
        descriptor.shard_bytes = bytes;
        descriptor.original = input.data();
        descriptor.recovery = output.data();
        descriptor.scratch = scratch->data();
        descriptor.scratch_bytes = scratch->size();
    }

    Shards original;
    Shards recovery;
    std::vector<const void*> input;
    std::vector<void*> output;
    std::unique_ptr<AlignedBuffer> scratch;
    unsigned requested_parity;
    leo2_encode_batch_item descriptor;
};

uint64_t ExerciseBatchEntry(
    leo2_context* context,
    unsigned k,
    unsigned r,
    size_t bytes,
    size_t item_count,
    bool use_binding)
{
    leo2_codec* codec = CreateCodec(context, k, r);
    size_t scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, bytes, &scratch_bytes), "sparse-high batch scratch query");
    std::vector<std::unique_ptr<BatchStorage> > storage;
    std::vector<leo2_encode_batch_item> items(item_count);
    for (size_t i = 0; i < item_count; ++i)
    {
        storage.push_back(std::unique_ptr<BatchStorage>(new BatchStorage(
            k, r, bytes, static_cast<unsigned>(i % r), scratch_bytes,
            UINT64_C(0x4241544348513141) + i)));
        items[i] = storage[i]->descriptor;
    }
    Require(
        leopard2_internal::ArmContextHighSparseEncodeRouteWitnessForDiagnostics(context),
        "sparse-high batch witness could not be armed");
    if (use_binding)
    {
        leo2_encode_batch_binding* binding = NULL;
        RequireResult(leo2_encode_batch_binding_create(
            codec, items.data(), items.size(), &binding),
            "sparse-high batch binding create");
        RequireResult(leo2_encode_batch_binding_execute(binding),
            "sparse-high batch binding execute");
        RequireResult(leo2_encode_batch_binding_execute(binding),
            "sparse-high reusable batch binding execute");
        leo2_encode_batch_binding_destroy(binding);
    }
    else
    {
        RequireResult(leo2_encode_batch(codec, items.data(), items.size()),
            "sparse-high batch encode");
    }
    const uint64_t expected_calls = use_binding
        ? 2U * item_count : item_count;
    CheckWitness(context, expected_calls, 0);

    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(
            field, leopard2_test::make_profile_layout(
                leopard2_test::kLegacyHigh, k, r));
    for (size_t i = 0; i < item_count; ++i)
    {
        const unsigned parity = storage[i]->requested_parity;
        Require(storage[i]->recovery[parity] == OracleParity(
            field, generator, storage[i]->original, parity),
            "sparse-high batch parity differs from oracle");
        for (unsigned other = 0; other < r; ++other)
        {
            if (other != parity)
                Require(storage[i]->recovery[other] == Bytes(bytes, 0xa5),
                    "sparse-high batch modified a null output");
        }
    }
    leo2_codec_destroy(codec);
    return expected_calls;
}

void ExerciseExclusions(leo2_context* context)
{
    leo2_codec* codec = CreateCodec(context, 2, 16);
    CheckPath(codec, 4096, 2, 16, false, "sparse-high Q2 exclusion");
    CheckPath(codec, 4096, 16, 16, false,
        "sparse-high full-output exclusion");
    CheckPath(codec, 4095, 1, 16, false,
        "sparse-high ragged-byte exclusion");
    CheckPath(codec, 1152, 1, 16, false,
        "sparse-high aligned-byte exclusion");
    const leopard2_test::BinaryField field =
        leopard2_test::make_legacy_gf8();
    const leopard2_test::Matrix generator =
        leopard2_test::direct_systematic_generator(
            field, leopard2_test::make_profile_layout(
                leopard2_test::kLegacyHigh, 2, 16));
    const Shards original = MakeOriginal(
        2, 4096, UINT64_C(0x4558434c55444551));
    Require(
        leopard2_internal::ArmContextHighSparseEncodeRouteWitnessForDiagnostics(context),
        "sparse-high Q2 exclusion witness could not be armed");
    EncodeTwoAndCheck(codec, original, field, generator, 16);
    CheckWitness(context, 0, 1);
    ExerciseTransformWitness(context, codec, 2, 16, 1152,
        UINT64_C(0x4558434c55444541));
    leo2_codec_destroy(codec);

    codec = CreateCodec(context, 3, 4);
    CheckPath(codec, 4032, 1, 4, false,
        "sparse-high non-boundary shape byte exclusion");
    CheckPath(codec, 4160, 1, 4, false,
        "sparse-high non-boundary shape upper exclusion");
    ExerciseTransformWitness(context, codec, 3, 4, 4032,
        UINT64_C(0x4558434c55444542));
    leo2_codec_destroy(codec);

    codec = CreateCodec(context, 5, 4);
    CheckPath(codec, 4096, 1,
        LEO2_EXPECT_HIGH_DIRECT_PRODUCTION ? 4 : 0, false,
        "sparse-high K-boundary exclusion");
    leo2_codec_destroy(codec);

    codec = CreateCodec(context, 4, 3);
    CheckPath(codec, 4096, 1,
        LEO2_EXPECT_HIGH_DIRECT_PRODUCTION ? 3 : 0, false,
        "sparse-high R-boundary exclusion");
    leo2_codec_destroy(codec);

    codec = CreateCodec(context, 2, 16, LEO2_CODEC_FORCE_GENERIC_DECODE);
    CheckPath(codec, 4096, 1, 16, false,
        "sparse-high codec-flags exclusion");
    ExerciseTransformWitness(context, codec, 2, 16, 4096,
        UINT64_C(0x4558434c55444546));
    leo2_codec_destroy(codec);
}

void ExerciseContextExclusions()
{
    leo2_context* explicit_context = CreateContext(
        LEO2_BACKEND_AVX2, 1, true, true);
    leo2_codec* codec = CreateCodec(explicit_context, 2, 16);
    CheckPath(codec, 4096, 1, 16, false,
        "sparse-high explicit-AVX2 request exclusion");
    ExerciseTransformWitness(explicit_context, codec, 2, 16, 4096,
        UINT64_C(0x4558504c49434954));
    leo2_codec_destroy(codec);
    leo2_context_destroy(explicit_context);

    leo2_context* threaded_context = CreateContext(
        LEO2_BACKEND_AUTO, 4, true, true);
    codec = CreateCodec(threaded_context, 2, 16);
    CheckPath(codec, 4096, 1, 16, false,
        "sparse-high thread-count exclusion");
    ExerciseTransformWitness(threaded_context, codec, 2, 16, 4096,
        UINT64_C(0x5448524541444544));
    leo2_codec_destroy(codec);
    leo2_context_destroy(threaded_context);
}

void CheckCompiledDefault()
{
    leo2_context_options options = {};
    options.struct_size = sizeof(options);
    options.backend = LEO2_BACKEND_AUTO;
    options.thread_count = 1;
    leo2_context* context = NULL;
    RequireResult(leo2_context_create(&options, &context),
        "sparse-high default context create");
    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(context, 2, 16,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
        "sparse-high default codec create");
    CheckPath(codec, 4096, 1,
        leo2_context_backend(context) == LEO2_BACKEND_AVX2 &&
            (LEO2_EXPECT_HIGH_DIRECT_PRODUCTION ||
             LEO2_EXPECT_HIGH_SPARSE_DIRECT) ? 16 : 0,
        leo2_context_backend(context) == LEO2_BACKEND_AVX2 &&
            LEO2_EXPECT_HIGH_SPARSE_DIRECT &&
            LEO2_EXPECT_HIGH_SPARSE_DIRECT_AUTO,
        "sparse-high compiled-default path query");
    leo2_codec_destroy(codec);
    leo2_context_destroy(context);
}

} // namespace

int main()
{
    try
    {
        Require(
            leopard2_internal::HighSparseDirectEncodeEnabled() ==
                (LEO2_EXPECT_HIGH_SPARSE_DIRECT != 0),
            "sparse-high support marker differs from build policy");
        Require(
            leopard2_internal::HighSparseDirectEncodeAutoEnabled() ==
                (LEO2_EXPECT_HIGH_SPARSE_DIRECT_AUTO != 0),
            "sparse-high AUTO marker differs from build policy");
        CheckCompiledDefault();

        if (!LEO2_EXPECT_HIGH_SPARSE_DIRECT)
        {
            Require(!leopard2_internal::
                    SetContextHighSparseDirectEncodePolicyForDiagnostics(
                        NULL, false, false),
                "uncompiled sparse-high policy accepted a null context");
            std::printf(
                "Leopard2 sparse-high production policy unavailable as expected\n");
            return 0;
        }

        Require(!leopard2_internal::
                SetContextHighSparseDirectEncodePolicyForDiagnostics(
                    NULL, false, false),
            "sparse-high policy accepted a null context");
        leo2_context_options invalid_options = {};
        invalid_options.struct_size = sizeof(invalid_options);
        invalid_options.backend = LEO2_BACKEND_AUTO;
        invalid_options.thread_count = 1;
        leo2_context* invalid_context = NULL;
        RequireResult(leo2_context_create(
            &invalid_options, &invalid_context),
            "sparse-high invalid-policy context create");
        Require(!leopard2_internal::
                SetContextHighSparseDirectEncodePolicyForDiagnostics(
                    invalid_context, false, true),
            "sparse-high policy accepted AUTO without tables");
        const bool auto_is_avx2 =
            leo2_context_backend(invalid_context) == LEO2_BACKEND_AVX2;
        leo2_context_destroy(invalid_context);
        if (!auto_is_avx2)
        {
            std::printf(
                "Leopard2 sparse-high production policy skipped: "
                "AUTO baseline is not AVX2\n");
            return 0;
        }

        ExercisePolicyState(false, false, false);
        ExercisePolicyState(true, false, false);
        ExercisePolicyState(true, true, true);

        leo2_context* context = CreateContext(
            LEO2_BACKEND_AUTO, 1, true, true);
        uint64_t route_checks = 0;
        const uint64_t parity_checks =
            ExerciseCandidateCells(context, route_checks);
        Require(parity_checks == 288,
            "sparse-high candidate parity-row count changed");
        Require(route_checks == 4U * parity_checks,
            "sparse-high candidate API route count changed");
        route_checks += ExerciseBatchEntry(
            context, 2, 16, 4096, 1, false);
        route_checks += ExerciseBatchEntry(
            context, 2, 16, 4096, 1, true);
        ExerciseExclusions(context);
        leo2_context_destroy(context);
        ExerciseContextExclusions();

        std::printf(
            "Leopard2 sparse-high production policy passed: "
            "states=3 tuples=36 parity_rows=%llu routes=%llu "
            "apis=one-shot,batch1,batch2,batch4,batch8,batch16,"
            "binding1,binding2,binding4,binding8,binding16 "
            "binding_reuse=2\n",
            static_cast<unsigned long long>(parity_checks),
            static_cast<unsigned long long>(route_checks));
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr,
            "Leopard2 sparse-high production policy failed: %s\n",
            error.what());
        return 1;
    }
}
