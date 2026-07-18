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

#include "leopard2.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

namespace {

typedef std::vector<uint8_t> Bytes;
typedef std::vector<Bytes> Shards;

void Require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void RequireResult(
    leo2_result actual,
    leo2_result expected,
    const std::string& operation)
{
    if (actual == expected)
        return;
    throw std::runtime_error(operation + ": got " +
        leo2_result_string(actual) + ", expected " +
        leo2_result_string(expected));
}

class AlignedBuffer
{
public:
    explicit AlignedBuffer(size_t bytes)
        : data_(NULL)
        , bytes_(bytes)
    {
        if (bytes_ == 0)
            return;
#if defined(_MSC_VER)
        data_ = _aligned_malloc(bytes_, leo2_scratch_alignment());
#else
        if (posix_memalign(&data_, leo2_scratch_alignment(), bytes_) != 0)
            data_ = NULL;
#endif
        if (!data_)
            throw std::bad_alloc();
        memset(data_, 0, bytes_);
    }

    ~AlignedBuffer()
    {
#if defined(_MSC_VER)
        _aligned_free(data_);
#else
        free(data_);
#endif
    }

    void* data() { return data_; }
    size_t size() const { return bytes_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);
    void* data_;
    size_t bytes_;
};

struct Fixture
{
    static const size_t kShortBytes = 17;
    static const size_t kLongBytes = 33;

    explicit Fixture(uint32_t thread_count)
        : context(NULL)
        , codec(NULL)
        , plan(NULL)
        , source_a(3, Bytes(kLongBytes, 0))
        , source_b(3, Bytes(kLongBytes, 0))
        , parity_a(2, Bytes(kLongBytes, 0))
        , parity_b(2, Bytes(kLongBytes, 0))
        , original_a(3, NULL)
        , original_b(3, NULL)
        , parity_a_mutable(2, NULL)
        , parity_b_mutable(2, NULL)
        , encode_short(0)
        , encode_long(0)
        , decode_short(0)
        , decode_long(0)
    {
        leo2_context_options options;
        memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_SCALAR;
        options.thread_count = thread_count;
        RequireResult(leo2_context_create(&options, &context),
            LEO2_SUCCESS, "context create");
        RequireResult(leo2_codec_create(context, 3, 2,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
            LEO2_SUCCESS, "codec create");

        Fill(source_a, 0x12345678u);
        Fill(source_b, 0x9abcdef0u);
        for (size_t i = 0; i < 3; ++i)
        {
            original_a[i] = &source_a[i][0];
            original_b[i] = &source_b[i][0];
        }
        for (size_t i = 0; i < 2; ++i)
        {
            parity_a_mutable[i] = &parity_a[i][0];
            parity_b_mutable[i] = &parity_b[i][0];
        }

        RequireResult(leo2_encode_scratch_size(codec, kShortBytes,
            &encode_short), LEO2_SUCCESS, "short encode scratch query");
        RequireResult(leo2_encode_scratch_size(codec, kLongBytes,
            &encode_long), LEO2_SUCCESS, "long encode scratch query");
        AlignedBuffer scratch_a(encode_long);
        AlignedBuffer scratch_b(encode_long);
        RequireResult(leo2_encode(codec, kLongBytes, &original_a[0],
            &parity_a_mutable[0], scratch_a.data(), scratch_a.size()),
            LEO2_SUCCESS, "encode source A");
        RequireResult(leo2_encode(codec, kLongBytes, &original_b[0],
            &parity_b_mutable[0], scratch_b.data(), scratch_b.size()),
            LEO2_SUCCESS, "encode source B");

        const uint8_t original_present[3] = { 0, 1, 1 };
        const uint8_t recovery_present[2] = { 1, 0 };
        RequireResult(leo2_decode_plan_create(codec, original_present,
            recovery_present, &plan), LEO2_SUCCESS, "plan create");
        RequireResult(leo2_decode_plan_scratch_size(plan, kShortBytes,
            &decode_short), LEO2_SUCCESS, "short decode scratch query");
        RequireResult(leo2_decode_plan_scratch_size(plan, kLongBytes,
            &decode_long), LEO2_SUCCESS, "long decode scratch query");
    }

    ~Fixture()
    {
        leo2_decode_plan_destroy(plan);
        leo2_codec_destroy(codec);
        leo2_context_destroy(context);
    }

    static void Fill(Shards& shards, uint32_t state)
    {
        for (size_t shard = 0; shard < shards.size(); ++shard)
            for (size_t i = 0; i < shards[shard].size(); ++i)
            {
                state = state * 1664525u + 1013904223u;
                shards[shard][i] = static_cast<uint8_t>(state >> 24);
            }
    }

    leo2_context* context;
    leo2_codec* codec;
    leo2_decode_plan* plan;
    Shards source_a;
    Shards source_b;
    Shards parity_a;
    Shards parity_b;
    std::vector<const void*> original_a;
    std::vector<const void*> original_b;
    std::vector<void*> parity_a_mutable;
    std::vector<void*> parity_b_mutable;
    size_t encode_short;
    size_t encode_long;
    size_t decode_short;
    size_t decode_long;

private:
    Fixture(const Fixture&);
    Fixture& operator=(const Fixture&);
};

void TestValidSharedInputs(Fixture& fixture)
{
    Shards output_a(2, Bytes(Fixture::kLongBytes, 0xa5));
    Shards output_b(2, Bytes(Fixture::kLongBytes, 0xa5));
    void* output_a_ptr[2] = { &output_a[0][0], &output_a[1][0] };
    void* output_b_ptr[2] = { &output_b[0][0], &output_b[1][0] };
    AlignedBuffer scratch_a(fixture.encode_long);
    AlignedBuffer scratch_b(fixture.encode_long);
    leo2_encode_batch_item encode_items[2] = {
        { Fixture::kLongBytes, &fixture.original_a[0], output_a_ptr,
          scratch_a.data(), scratch_a.size() },
        { Fixture::kLongBytes, &fixture.original_a[0], output_b_ptr,
          scratch_b.data(), scratch_b.size() }
    };
    RequireResult(leo2_encode_batch(fixture.codec, encode_items, 2),
        LEO2_SUCCESS, "shared-input encode batch");
    Require(output_a == fixture.parity_a && output_b == fixture.parity_a,
        "shared-input encode parity mismatch");

    const void* original[3] = {
        NULL, fixture.original_a[1], fixture.original_a[2]
    };
    const void* recovery[2] = { &fixture.parity_a[0][0], NULL };
    Bytes restored_a(Fixture::kLongBytes, 0xcc);
    Bytes restored_b(Fixture::kLongBytes, 0xcc);
    void* restored_a_ptr[3] = { &restored_a[0], NULL, NULL };
    void* restored_b_ptr[3] = { &restored_b[0], NULL, NULL };
    AlignedBuffer decode_scratch_a(fixture.decode_long);
    AlignedBuffer decode_scratch_b(fixture.decode_long);
    leo2_decode_batch_item decode_items[2] = {
        { Fixture::kLongBytes, original, recovery, restored_a_ptr,
          decode_scratch_a.data(), decode_scratch_a.size() },
        { Fixture::kLongBytes, original, recovery, restored_b_ptr,
          decode_scratch_b.data(), decode_scratch_b.size() }
    };
    RequireResult(leo2_decode_plan_execute_batch(
        fixture.plan, decode_items, 2), LEO2_SUCCESS,
        "shared-input decode batch");
    Require(restored_a == fixture.source_a[0] &&
            restored_b == fixture.source_a[0],
        "shared-input decode mismatch");
}

void TestLargeBatchCount(Fixture& fixture)
{
    /* Cross 8-bit loop/counter boundaries without imposing an arbitrary
       production batch cap.  All items deliberately share immutable inputs. */
    const size_t item_count = 257;
    std::vector<Bytes> outputs(
        item_count, Bytes(Fixture::kLongBytes, 0xa5));
    std::vector<std::vector<void*> > output_pointers(
        item_count, std::vector<void*>(2, NULL));
    std::vector<std::unique_ptr<AlignedBuffer> > scratches(item_count);
    std::vector<leo2_encode_batch_item> items(item_count);
    for (size_t i = 0; i < item_count; ++i)
    {
        output_pointers[i][0] = &outputs[i][0];
        scratches[i].reset(new AlignedBuffer(
            std::max(fixture.encode_long, fixture.decode_long)));
        items[i].shard_bytes = Fixture::kLongBytes;
        items[i].original = &fixture.original_a[0];
        items[i].recovery = &output_pointers[i][0];
        items[i].scratch = scratches[i]->data();
        items[i].scratch_bytes = scratches[i]->size();
    }
    RequireResult(leo2_encode_batch(
        fixture.codec, &items[0], items.size()), LEO2_SUCCESS,
        "257-item encode batch");
    for (size_t i = 0; i < item_count; ++i)
        Require(outputs[i] == fixture.parity_a[0],
            "257-item encode batch parity mismatch");

    const void* original[3] = {
        NULL, fixture.original_a[1], fixture.original_a[2]
    };
    std::vector<Bytes> restored(
        item_count, Bytes(Fixture::kLongBytes, 0xcc));
    std::vector<std::vector<const void*> > recovery_pointers(
        item_count, std::vector<const void*>(2, NULL));
    std::vector<std::vector<void*> > restored_pointers(
        item_count, std::vector<void*>(3, NULL));
    std::vector<leo2_decode_batch_item> decode_items(item_count);
    for (size_t i = 0; i < item_count; ++i)
    {
        recovery_pointers[i][0] = &outputs[i][0];
        restored_pointers[i][0] = &restored[i][0];
        decode_items[i].shard_bytes = Fixture::kLongBytes;
        decode_items[i].original = original;
        decode_items[i].recovery = &recovery_pointers[i][0];
        decode_items[i].restored_original = &restored_pointers[i][0];
        decode_items[i].scratch = scratches[i]->data();
        decode_items[i].scratch_bytes = scratches[i]->size();
    }
    RequireResult(leo2_decode_plan_execute_batch(
        fixture.plan, &decode_items[0], decode_items.size()),
        LEO2_SUCCESS, "257-item decode batch");
    for (size_t i = 0; i < item_count; ++i)
        Require(restored[i] == fixture.source_a[0],
            "257-item decode batch recovery mismatch");
}

void TestNoLossBatchIsTrueNoOp(Fixture& fixture)
{
    const uint8_t original_present[3] = { 1, 1, 1 };
    const uint8_t recovery_present[2] = { 1, 1 };
    leo2_decode_plan* no_loss = NULL;
    RequireResult(leo2_decode_plan_create(fixture.codec,
        original_present, recovery_present, &no_loss), LEO2_SUCCESS,
        "no-loss plan create");

    leo2_decode_batch_item items[2];
    memset(items, 0, sizeof(items));
    for (size_t i = 0; i < 2; ++i)
    {
        items[i].shard_bytes = 0;
        items[i].original = reinterpret_cast<const void* const*>(
            static_cast<uintptr_t>(1));
        items[i].recovery = reinterpret_cast<const void* const*>(
            static_cast<uintptr_t>(1));
        items[i].restored_original = reinterpret_cast<void* const*>(
            static_cast<uintptr_t>(1));
        items[i].scratch = reinterpret_cast<void*>(
            static_cast<uintptr_t>(1));
        items[i].scratch_bytes = 1;
    }
    const leo2_result result = leo2_decode_plan_execute_batch(
        no_loss, items, 2);
    leo2_decode_plan_destroy(no_loss);
    RequireResult(result, LEO2_SUCCESS,
        "no-loss batch inspected per-item execution state");
}

void TestSingleEncodeItemFastPath(Fixture& fixture)
{
    Shards output(2, Bytes(Fixture::kLongBytes, 0xa5));
    void* output_pointers[2] = { &output[0][0], &output[1][0] };
    AlignedBuffer scratch(fixture.encode_long);
    leo2_encode_batch_item item = {
        Fixture::kLongBytes, &fixture.original_a[0], output_pointers,
        scratch.data(), scratch.size()
    };
    RequireResult(leo2_encode_batch(fixture.codec, &item, 1),
        LEO2_SUCCESS, "single-item encode batch");
    Require(output == fixture.parity_a,
        "single-item encode batch parity mismatch");

    /* The fast path must protect the descriptor just like the general batch
       preflight.  Rejection must happen before either metadata or scratch is
       modified. */
    Bytes safe_output(Fixture::kLongBytes, 0x4d);
    alignas(64) leo2_encode_batch_item descriptor_item;
    void* descriptor_outputs[2] = { &descriptor_item, &safe_output[0] };
    descriptor_item.shard_bytes = Fixture::kLongBytes;
    descriptor_item.original = &fixture.original_a[0];
    descriptor_item.recovery = descriptor_outputs;
    descriptor_item.scratch = scratch.data();
    descriptor_item.scratch_bytes = scratch.size();
    const leo2_encode_batch_item descriptor_before = descriptor_item;
    RequireResult(leo2_encode_batch(
        fixture.codec, &descriptor_item, 1), LEO2_OVERLAP,
        "single-item output/descriptor overlap");
    Require(descriptor_item.shard_bytes == descriptor_before.shard_bytes &&
            descriptor_item.original == descriptor_before.original &&
            descriptor_item.recovery == descriptor_before.recovery &&
            descriptor_item.scratch == descriptor_before.scratch &&
            descriptor_item.scratch_bytes == descriptor_before.scratch_bytes,
        "single-item rejection modified descriptor metadata");

    alignas(64) leo2_encode_batch_item scratch_item;
    scratch_item.shard_bytes = Fixture::kLongBytes;
    scratch_item.original = &fixture.original_a[0];
    scratch_item.recovery = output_pointers;
    scratch_item.scratch = &scratch_item;
    scratch_item.scratch_bytes = fixture.encode_long;
    const leo2_encode_batch_item scratch_before = scratch_item;
    RequireResult(leo2_encode_batch(
        fixture.codec, &scratch_item, 1), LEO2_OVERLAP,
        "single-item scratch/descriptor overlap");
    Require(scratch_item.shard_bytes == scratch_before.shard_bytes &&
            scratch_item.original == scratch_before.original &&
            scratch_item.recovery == scratch_before.recovery &&
            scratch_item.scratch == scratch_before.scratch &&
            scratch_item.scratch_bytes == scratch_before.scratch_bytes,
        "single-item scratch rejection modified descriptor metadata");

    void* duplicate_outputs[2] = { &output[0][0], &output[0][0] };
    item.recovery = duplicate_outputs;
    RequireResult(leo2_encode_batch(fixture.codec, &item, 1),
        LEO2_OVERLAP, "single-item duplicate output overlap");

    void* input_output[2] = {
        const_cast<void*>(fixture.original_a[0]), &output[1][0]
    };
    item.recovery = input_output;
    RequireResult(leo2_encode_batch(fixture.codec, &item, 1),
        LEO2_OVERLAP, "single-item output/input overlap");
}

void TestEncodeConflicts(Fixture& fixture)
{
    void* no_output[2] = { NULL, NULL };

    /* A 17-byte output intersects only the second item's 33-byte input. */
    const Shards source_b_before = fixture.source_b;
    void* cross_input_output[2] = { &fixture.source_b[0][16], NULL };
    AlignedBuffer cross_scratch_a(fixture.encode_short);
    AlignedBuffer cross_scratch_b(fixture.encode_long);
    leo2_encode_batch_item cross_input_items[2] = {
        { Fixture::kShortBytes, &fixture.original_a[0], cross_input_output,
          cross_scratch_a.data(), cross_scratch_a.size() },
        { Fixture::kLongBytes, &fixture.original_b[0], no_output,
          cross_scratch_b.data(), cross_scratch_b.size() }
    };
    RequireResult(leo2_encode_batch(
        fixture.codec, cross_input_items, 2), LEO2_OVERLAP,
        "encode cross-item output/input overlap");
    Require(fixture.source_b == source_b_before,
        "rejected encode cross-input overlap changed input bytes");

    Bytes collision(64, 0x6c);
    const Bytes collision_before = collision;
    void* partial_a[2] = { &collision[0], NULL };
    void* partial_b[2] = { &collision[16], NULL };
    AlignedBuffer partial_scratch_a(fixture.encode_short);
    AlignedBuffer partial_scratch_b(fixture.encode_long);
    leo2_encode_batch_item partial_items[2] = {
        { Fixture::kShortBytes, &fixture.original_a[0], partial_a,
          partial_scratch_a.data(), partial_scratch_a.size() },
        { Fixture::kLongBytes, &fixture.original_b[0], partial_b,
          partial_scratch_b.data(), partial_scratch_b.size() }
    };
    RequireResult(leo2_encode_batch(fixture.codec, partial_items, 2),
        LEO2_OVERLAP, "encode partial output/output overlap");
    Require(collision == collision_before,
        "rejected encode partial output overlap changed bytes");
    partial_b[0] = partial_a[0];
    RequireResult(leo2_encode_batch(fixture.codec, partial_items, 2),
        LEO2_OVERLAP, "encode duplicate output/output overlap");

    Require(fixture.encode_short > leo2_scratch_alignment(),
        "encode scratch too small for partial-overlap fixture");
    AlignedBuffer shared_scratch(
        std::max(fixture.encode_short, fixture.encode_long) +
        leo2_scratch_alignment());
    Shards scratch_outputs(2, Bytes(Fixture::kLongBytes, 0xa5));
    void* scratch_output_a[2] = { &scratch_outputs[0][0], NULL };
    void* scratch_output_b[2] = { &scratch_outputs[1][0], NULL };
    leo2_encode_batch_item scratch_items[2] = {
        { Fixture::kShortBytes, &fixture.original_a[0], scratch_output_a,
          shared_scratch.data(), fixture.encode_short },
        { Fixture::kLongBytes, &fixture.original_b[0], scratch_output_b,
          static_cast<uint8_t*>(shared_scratch.data()) +
              leo2_scratch_alignment(), fixture.encode_long }
    };
    RequireResult(leo2_encode_batch(fixture.codec, scratch_items, 2),
        LEO2_OVERLAP, "encode partial scratch/scratch overlap");

    const void* original_metadata[3] = {
        fixture.original_b[0], fixture.original_b[1], fixture.original_b[2]
    };
    const void* original_metadata_before[3] = {
        original_metadata[0], original_metadata[1], original_metadata[2]
    };
    Bytes metadata_other_output(Fixture::kShortBytes, 0x4d);
    void* metadata_output[2] = {
        const_cast<void*>(static_cast<const void*>(original_metadata)),
        &metadata_other_output[0]
    };
    AlignedBuffer metadata_scratch_a(fixture.encode_short);
    AlignedBuffer metadata_scratch_b(fixture.encode_long);
    memset(metadata_scratch_a.data(), 0xa3, metadata_scratch_a.size());
    memset(metadata_scratch_b.data(), 0xb4, metadata_scratch_b.size());
    Bytes metadata_scratch_a_before(metadata_scratch_a.size());
    Bytes metadata_scratch_b_before(metadata_scratch_b.size());
    memcpy(&metadata_scratch_a_before[0], metadata_scratch_a.data(),
        metadata_scratch_a.size());
    memcpy(&metadata_scratch_b_before[0], metadata_scratch_b.data(),
        metadata_scratch_b.size());
    leo2_encode_batch_item metadata_items[2] = {
        { Fixture::kShortBytes, &fixture.original_a[0], metadata_output,
          metadata_scratch_a.data(), metadata_scratch_a.size() },
        { Fixture::kLongBytes, original_metadata, no_output,
          metadata_scratch_b.data(), metadata_scratch_b.size() }
    };
    RequireResult(leo2_encode_batch(fixture.codec, metadata_items, 2),
        LEO2_OVERLAP, "encode cross-item output/metadata overlap");
    Require(memcmp(original_metadata, original_metadata_before,
                sizeof(original_metadata)) == 0,
        "rejected encode metadata overlap changed pointer metadata");
    Require(memcmp(metadata_scratch_a.data(),
                &metadata_scratch_a_before[0], metadata_scratch_a.size()) == 0 &&
            memcmp(metadata_scratch_b.data(),
                &metadata_scratch_b_before[0], metadata_scratch_b.size()) == 0,
        "rejected encode metadata overlap changed scratch");

    Bytes untouched(Fixture::kLongBytes, 0xa5);
    const Bytes untouched_before = untouched;
    Shards later_output(2, Bytes(Fixture::kLongBytes, 0xa5));
    void* first_output[2] = { &untouched[0], NULL };
    void* later_output_ptr[2] = {
        &later_output[0][0], &later_output[1][0]
    };
    AlignedBuffer valid_scratch(fixture.encode_long);
    AlignedBuffer short_scratch(fixture.encode_long);
    leo2_encode_batch_item invalid_later[2] = {
        { Fixture::kLongBytes, &fixture.original_a[0], first_output,
          valid_scratch.data(), valid_scratch.size() },
        { Fixture::kLongBytes, &fixture.original_b[0], later_output_ptr,
          short_scratch.data(), fixture.encode_long - 1 }
    };
    RequireResult(leo2_encode_batch(fixture.codec, invalid_later, 2),
        LEO2_SCRATCH_TOO_SMALL, "encode later-item preflight failure");
    Require(untouched == untouched_before,
        "encode executed an earlier item before later validation failed");
}

void TestDecodeConflicts(Fixture& fixture)
{
    const void* original_a[3] = {
        NULL, fixture.original_a[1], fixture.original_a[2]
    };
    const void* recovery_a[2] = { &fixture.parity_a[0][0], NULL };
    const void* original_b[3] = {
        NULL, fixture.original_b[1], fixture.original_b[2]
    };
    const void* recovery_b[2] = { &fixture.parity_b[0][0], NULL };

    const Shards source_b_before = fixture.source_b;
    void* cross_output[3] = {
        &fixture.source_b[1][16], NULL, NULL
    };
    Bytes ordinary_output(Fixture::kLongBytes, 0xcc);
    void* ordinary_output_ptr[3] = { &ordinary_output[0], NULL, NULL };
    AlignedBuffer cross_scratch_a(fixture.decode_short);
    AlignedBuffer cross_scratch_b(fixture.decode_long);
    leo2_decode_batch_item cross_items[2] = {
        { Fixture::kShortBytes, original_a, recovery_a, cross_output,
          cross_scratch_a.data(), cross_scratch_a.size() },
        { Fixture::kLongBytes, original_b, recovery_b, ordinary_output_ptr,
          cross_scratch_b.data(), cross_scratch_b.size() }
    };
    RequireResult(leo2_decode_plan_execute_batch(
        fixture.plan, cross_items, 2), LEO2_OVERLAP,
        "decode cross-item output/input overlap");
    Require(fixture.source_b == source_b_before,
        "rejected decode cross-input overlap changed input bytes");

    Bytes collision(64, 0x6c);
    const Bytes collision_before = collision;
    void* partial_a[3] = { &collision[0], NULL, NULL };
    void* partial_b[3] = { &collision[16], NULL, NULL };
    AlignedBuffer partial_scratch_a(fixture.decode_short);
    AlignedBuffer partial_scratch_b(fixture.decode_long);
    leo2_decode_batch_item partial_items[2] = {
        { Fixture::kShortBytes, original_a, recovery_a, partial_a,
          partial_scratch_a.data(), partial_scratch_a.size() },
        { Fixture::kLongBytes, original_b, recovery_b, partial_b,
          partial_scratch_b.data(), partial_scratch_b.size() }
    };
    RequireResult(leo2_decode_plan_execute_batch(
        fixture.plan, partial_items, 2), LEO2_OVERLAP,
        "decode partial output/output overlap");
    Require(collision == collision_before,
        "rejected decode partial output overlap changed bytes");
    partial_b[0] = partial_a[0];
    RequireResult(leo2_decode_plan_execute_batch(
        fixture.plan, partial_items, 2), LEO2_OVERLAP,
        "decode duplicate output/output overlap");

    Require(fixture.decode_short > leo2_scratch_alignment(),
        "decode scratch too small for partial-overlap fixture");
    AlignedBuffer shared_scratch(
        std::max(fixture.decode_short, fixture.decode_long) +
        leo2_scratch_alignment());
    Bytes scratch_output_a(Fixture::kShortBytes, 0xcc);
    Bytes scratch_output_b(Fixture::kLongBytes, 0xcc);
    void* scratch_output_a_ptr[3] = {
        &scratch_output_a[0], NULL, NULL
    };
    void* scratch_output_b_ptr[3] = {
        &scratch_output_b[0], NULL, NULL
    };
    leo2_decode_batch_item scratch_items[2] = {
        { Fixture::kShortBytes, original_a, recovery_a,
          scratch_output_a_ptr, shared_scratch.data(), fixture.decode_short },
        { Fixture::kLongBytes, original_b, recovery_b,
          scratch_output_b_ptr,
          static_cast<uint8_t*>(shared_scratch.data()) +
              leo2_scratch_alignment(), fixture.decode_long }
    };
    RequireResult(leo2_decode_plan_execute_batch(
        fixture.plan, scratch_items, 2), LEO2_OVERLAP,
        "decode partial scratch/scratch overlap");

    const uint8_t metadata_original_present[3] = { 0, 0, 1 };
    const uint8_t metadata_recovery_present[2] = { 1, 1 };
    leo2_decode_plan* metadata_plan = NULL;
    RequireResult(leo2_decode_plan_create(fixture.codec,
        metadata_original_present, metadata_recovery_present,
        &metadata_plan), LEO2_SUCCESS, "metadata-overlap plan create");
    size_t metadata_decode_short = 0;
    size_t metadata_decode_long = 0;
    RequireResult(leo2_decode_plan_scratch_size(metadata_plan,
        Fixture::kShortBytes, &metadata_decode_short), LEO2_SUCCESS,
        "metadata-overlap short scratch query");
    RequireResult(leo2_decode_plan_scratch_size(metadata_plan,
        Fixture::kLongBytes, &metadata_decode_long), LEO2_SUCCESS,
        "metadata-overlap long scratch query");

    const void* metadata_original_a[3] = {
        NULL, NULL, fixture.original_a[2]
    };
    const void* metadata_recovery_a[2] = {
        &fixture.parity_a[0][0], &fixture.parity_a[1][0]
    };
    const void* original_metadata[3] = {
        NULL, NULL, fixture.original_b[2]
    };
    const void* metadata_recovery_b[2] = {
        &fixture.parity_b[0][0], &fixture.parity_b[1][0]
    };
    const void* original_metadata_before[3] = {
        original_metadata[0], original_metadata[1], original_metadata[2]
    };
    Bytes metadata_other_output(Fixture::kShortBytes, 0x5e);
    void* metadata_output[3] = {
        const_cast<void*>(static_cast<const void*>(original_metadata)),
        &metadata_other_output[0], NULL
    };
    Bytes metadata_output_b0(Fixture::kLongBytes, 0x6f);
    Bytes metadata_output_b1(Fixture::kLongBytes, 0x70);
    void* metadata_output_b[3] = {
        &metadata_output_b0[0], &metadata_output_b1[0], NULL
    };
    AlignedBuffer metadata_scratch_a(metadata_decode_short);
    AlignedBuffer metadata_scratch_b(metadata_decode_long);
    memset(metadata_scratch_a.data(), 0xc5, metadata_scratch_a.size());
    memset(metadata_scratch_b.data(), 0xd6, metadata_scratch_b.size());
    Bytes metadata_scratch_a_before(metadata_scratch_a.size());
    Bytes metadata_scratch_b_before(metadata_scratch_b.size());
    memcpy(&metadata_scratch_a_before[0], metadata_scratch_a.data(),
        metadata_scratch_a.size());
    memcpy(&metadata_scratch_b_before[0], metadata_scratch_b.data(),
        metadata_scratch_b.size());
    leo2_decode_batch_item metadata_items[2] = {
        { Fixture::kShortBytes, metadata_original_a, metadata_recovery_a,
          metadata_output,
          metadata_scratch_a.data(), metadata_scratch_a.size() },
        { Fixture::kLongBytes, original_metadata, metadata_recovery_b,
          metadata_output_b, metadata_scratch_b.data(),
          metadata_scratch_b.size() }
    };
    const leo2_result metadata_result = leo2_decode_plan_execute_batch(
        metadata_plan, metadata_items, 2);
    leo2_decode_plan_destroy(metadata_plan);
    RequireResult(metadata_result, LEO2_OVERLAP,
        "decode cross-item output/metadata overlap");
    Require(memcmp(original_metadata, original_metadata_before,
                sizeof(original_metadata)) == 0,
        "rejected decode metadata overlap changed pointer metadata");
    Require(memcmp(metadata_scratch_a.data(),
                &metadata_scratch_a_before[0], metadata_scratch_a.size()) == 0 &&
            memcmp(metadata_scratch_b.data(),
                &metadata_scratch_b_before[0], metadata_scratch_b.size()) == 0,
        "rejected decode metadata overlap changed scratch");

    Bytes untouched(Fixture::kLongBytes, 0xcc);
    const Bytes untouched_before = untouched;
    Bytes later_output(Fixture::kLongBytes, 0xcc);
    void* first_output[3] = { &untouched[0], NULL, NULL };
    void* later_output_ptr[3] = { &later_output[0], NULL, NULL };
    AlignedBuffer valid_scratch(fixture.decode_long);
    AlignedBuffer short_scratch(fixture.decode_long);
    leo2_decode_batch_item invalid_later[2] = {
        { Fixture::kLongBytes, original_a, recovery_a, first_output,
          valid_scratch.data(), valid_scratch.size() },
        { Fixture::kLongBytes, original_b, recovery_b, later_output_ptr,
          short_scratch.data(), fixture.decode_long - 1 }
    };
    RequireResult(leo2_decode_plan_execute_batch(
        fixture.plan, invalid_later, 2), LEO2_SCRATCH_TOO_SMALL,
        "decode later-item preflight failure");
    Require(untouched == untouched_before,
        "decode executed an earlier item before later validation failed");
}

void Run(uint32_t thread_count)
{
    Fixture fixture(thread_count);
    TestValidSharedInputs(fixture);
    TestSingleEncodeItemFastPath(fixture);
    TestEncodeConflicts(fixture);
    TestDecodeConflicts(fixture);
    TestValidSharedInputs(fixture);
    TestLargeBatchCount(fixture);
    TestNoLossBatchIsTrueNoOp(fixture);
}

} // namespace

int main()
{
    try
    {
        Run(1);
        Run(4);
        std::cout << "leopard2 batch aliasing passed: contexts=2 "
                  << "valid_shared=4 conflict_checks=20 "
                  << "large_encode_items=514 large_decode_items=514"
                  << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "leopard2 batch aliasing failed: "
                  << error.what() << std::endl;
        return 1;
    }
}
