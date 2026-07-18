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
#include "allocation_audit_config.h"

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
static std::atomic<bool> g_track_allocations(false);
static std::atomic<uint64_t> g_tracked_allocations(0);

#if defined(_MSC_VER)
#define LEO2_BATCH_TEST_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define LEO2_BATCH_TEST_NOINLINE __attribute__((noinline))
#else
#define LEO2_BATCH_TEST_NOINLINE
#endif

LEO2_BATCH_TEST_NOINLINE void* operator new(size_t bytes)
{
    if (g_track_allocations.load(std::memory_order_relaxed))
        g_tracked_allocations.fetch_add(1, std::memory_order_relaxed);
    void* result = malloc(bytes == 0 ? 1 : bytes);
    if (!result)
        throw std::bad_alloc();
    return result;
}

LEO2_BATCH_TEST_NOINLINE void* operator new[](size_t bytes)
{
    return ::operator new(bytes);
}

LEO2_BATCH_TEST_NOINLINE void* operator new(
    size_t bytes, const std::nothrow_t&) noexcept
{
    try { return ::operator new(bytes); }
    catch (...) { return NULL; }
}

LEO2_BATCH_TEST_NOINLINE void* operator new[](
    size_t bytes, const std::nothrow_t&) noexcept
{
    return ::operator new(bytes, std::nothrow);
}

LEO2_BATCH_TEST_NOINLINE void operator delete(void* pointer) noexcept
{
    free(pointer);
}
LEO2_BATCH_TEST_NOINLINE void operator delete[](void* pointer) noexcept
{
    free(pointer);
}
LEO2_BATCH_TEST_NOINLINE void operator delete(
    void* pointer, const std::nothrow_t&) noexcept
{
    ::operator delete(pointer);
}
LEO2_BATCH_TEST_NOINLINE void operator delete[](
    void* pointer, const std::nothrow_t&) noexcept
{
    ::operator delete[](pointer);
}

#undef LEO2_BATCH_TEST_NOINLINE
#endif

static void BeginAllocationAudit()
{
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
    g_tracked_allocations.store(0, std::memory_order_relaxed);
    g_track_allocations.store(true, std::memory_order_release);
#endif
}

static uint64_t EndAllocationAudit()
{
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
    g_track_allocations.store(false, std::memory_order_release);
    return g_tracked_allocations.load(std::memory_order_relaxed);
#else
    return 0;
#endif
}

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
    const void* data() const { return data_; }
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

static size_t BatchItemBytes(size_t item)
{
    return (item & 1u) == 0
        ? Fixture::kShortBytes : Fixture::kLongBytes;
}

struct ScalableEncodeBatch
{
    ScalableEncodeBatch(Fixture& fixture, size_t item_count)
        : fixture(fixture)
        , outputs(item_count, Shards(2, Bytes(Fixture::kLongBytes, 0xa5)))
        , original_pointers(item_count, fixture.original_a)
        , recovery_pointers(item_count, std::vector<void*>(2, NULL))
        , scratches(item_count)
        , items(item_count)
        , preflight_bytes(0)
    {
        for (size_t item = 0; item < item_count; ++item)
        {
            const size_t bytes = BatchItemBytes(item);
            recovery_pointers[item][0] = &outputs[item][0][0];
            recovery_pointers[item][1] = &outputs[item][1][0];
            const size_t scratch_bytes = bytes == Fixture::kShortBytes
                ? fixture.encode_short : fixture.encode_long;
            scratches[item].reset(new AlignedBuffer(scratch_bytes));
            items[item].shard_bytes = bytes;
            items[item].original = &original_pointers[item][0];
            items[item].recovery = &recovery_pointers[item][0];
            items[item].scratch = scratches[item]->data();
            items[item].scratch_bytes = scratches[item]->size();
        }
        RequireResult(leo2_encode_batch_preflight_scratch_size(
            fixture.codec, item_count, &preflight_bytes), LEO2_SUCCESS,
            "scalable encode preflight scratch query");
        preflight.reset(new AlignedBuffer(preflight_bytes));
    }

    leo2_result Execute()
    {
        return leo2_encode_batch_with_preflight_scratch(fixture.codec,
            items.empty() ? NULL : &items[0], items.size(),
            preflight->data(), preflight->size());
    }

    void Check() const
    {
        for (size_t item = 0; item < items.size(); ++item)
        {
            const size_t bytes = static_cast<size_t>(items[item].shard_bytes);
            for (size_t parity = 0; parity < outputs[item].size(); ++parity)
            {
                Require(std::equal(outputs[item][parity].begin(),
                            outputs[item][parity].begin() + bytes,
                            fixture.parity_a[parity].begin()),
                    "scalable encode parity mismatch");
            }
        }
    }

    Fixture& fixture;
    std::vector<Shards> outputs;
    std::vector<std::vector<const void*> > original_pointers;
    std::vector<std::vector<void*> > recovery_pointers;
    std::vector<std::unique_ptr<AlignedBuffer> > scratches;
    std::vector<leo2_encode_batch_item> items;
    size_t preflight_bytes;
    std::unique_ptr<AlignedBuffer> preflight;

private:
    ScalableEncodeBatch(const ScalableEncodeBatch&);
    ScalableEncodeBatch& operator=(const ScalableEncodeBatch&);
};

struct ScalableDecodeBatch
{
    ScalableDecodeBatch(
        Fixture& fixture,
        const ScalableEncodeBatch& encoded)
        : fixture(fixture)
        , restored(encoded.items.size(),
              Bytes(Fixture::kLongBytes, 0xcc))
        , original_pointers(encoded.items.size(),
              std::vector<const void*>(3, NULL))
        , recovery_pointers(encoded.items.size(),
              std::vector<const void*>(2, NULL))
        , restored_pointers(encoded.items.size(),
              std::vector<void*>(3, NULL))
        , scratches(encoded.items.size())
        , items(encoded.items.size())
        , preflight_bytes(0)
    {
        for (size_t item = 0; item < items.size(); ++item)
        {
            const size_t bytes =
                static_cast<size_t>(encoded.items[item].shard_bytes);
            original_pointers[item][0] = NULL;
            original_pointers[item][1] = fixture.original_a[1];
            original_pointers[item][2] = fixture.original_a[2];
            recovery_pointers[item][0] = &encoded.outputs[item][0][0];
            recovery_pointers[item][1] = NULL;
            restored_pointers[item][0] = &restored[item][0];
            size_t scratch_bytes = 0;
            RequireResult(leo2_decode_plan_scratch_size(
                fixture.plan, bytes, &scratch_bytes), LEO2_SUCCESS,
                "scalable decode item scratch query");
            scratches[item].reset(new AlignedBuffer(scratch_bytes));
            items[item].shard_bytes = bytes;
            items[item].original = &original_pointers[item][0];
            items[item].recovery = &recovery_pointers[item][0];
            items[item].restored_original = &restored_pointers[item][0];
            items[item].scratch = scratches[item]->data();
            items[item].scratch_bytes = scratches[item]->size();
        }
        RequireResult(leo2_decode_plan_batch_preflight_scratch_size(
            fixture.plan, items.size(), &preflight_bytes), LEO2_SUCCESS,
            "scalable decode preflight scratch query");
        preflight.reset(new AlignedBuffer(preflight_bytes));
    }

    leo2_result Execute()
    {
        return leo2_decode_plan_execute_batch_with_preflight_scratch(
            fixture.plan, items.empty() ? NULL : &items[0], items.size(),
            preflight->data(), preflight->size());
    }

    void Check() const
    {
        for (size_t item = 0; item < items.size(); ++item)
        {
            const size_t bytes = static_cast<size_t>(items[item].shard_bytes);
            Require(std::equal(restored[item].begin(),
                        restored[item].begin() + bytes,
                        fixture.source_a[0].begin()),
                "scalable decode result mismatch");
        }
    }

    Fixture& fixture;
    std::vector<Bytes> restored;
    std::vector<std::vector<const void*> > original_pointers;
    std::vector<std::vector<const void*> > recovery_pointers;
    std::vector<std::vector<void*> > restored_pointers;
    std::vector<std::unique_ptr<AlignedBuffer> > scratches;
    std::vector<leo2_decode_batch_item> items;
    size_t preflight_bytes;
    std::unique_ptr<AlignedBuffer> preflight;

private:
    ScalableDecodeBatch(const ScalableDecodeBatch&);
    ScalableDecodeBatch& operator=(const ScalableDecodeBatch&);
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

void TestScalableQueriesAndNoOp(Fixture& fixture)
{
    size_t bytes = 123;
    RequireResult(leo2_encode_batch_preflight_scratch_size(
        fixture.codec, 1, &bytes), LEO2_SUCCESS,
        "single-item encode preflight query");
    Require(bytes == 0, "single-item encode preflight was not zero");
    RequireResult(leo2_encode_batch_preflight_scratch_size(
        fixture.codec, 8, &bytes), LEO2_SUCCESS,
        "eight-item encode preflight query");
    Require(bytes == 0, "small encode batch did not retain compatibility path");
    RequireResult(leo2_encode_batch_preflight_scratch_size(
        fixture.codec, 9, &bytes), LEO2_SUCCESS,
        "nine-item encode preflight query");
    Require(bytes != 0 && bytes % leo2_scratch_alignment() == 0,
        "scalable encode preflight size is not aligned");
    const size_t encode_nine = bytes;
    RequireResult(leo2_encode_batch_preflight_scratch_size(
        fixture.codec, 64, &bytes), LEO2_SUCCESS,
        "64-item encode preflight query");
    Require(bytes > encode_nine,
        "encode preflight scratch did not scale with batch size");

    RequireResult(leo2_decode_plan_batch_preflight_scratch_size(
        fixture.plan, 8, &bytes), LEO2_SUCCESS,
        "eight-item decode preflight query");
    Require(bytes == 0, "small decode batch did not retain compatibility path");
    RequireResult(leo2_decode_plan_batch_preflight_scratch_size(
        fixture.plan, 9, &bytes), LEO2_SUCCESS,
        "nine-item decode preflight query");
    Require(bytes != 0 && bytes % leo2_scratch_alignment() == 0,
        "scalable decode preflight size is not aligned");

    bytes = 123;
    RequireResult(leo2_encode_batch_preflight_scratch_size(
        NULL, 9, &bytes), LEO2_INVALID_ARGUMENT,
        "null-codec encode preflight query");
    Require(bytes == 0, "failed encode preflight query did not clear output");
    bytes = 123;
    RequireResult(leo2_decode_plan_batch_preflight_scratch_size(
        NULL, 9, &bytes), LEO2_INVALID_ARGUMENT,
        "null-plan decode preflight query");
    Require(bytes == 0, "failed decode preflight query did not clear output");
    RequireResult(leo2_encode_batch_preflight_scratch_size(
        fixture.codec, 9, NULL), LEO2_INVALID_ARGUMENT,
        "null encode preflight output");
    RequireResult(leo2_decode_plan_batch_preflight_scratch_size(
        fixture.plan, 9, NULL), LEO2_INVALID_ARGUMENT,
        "null decode preflight output");

    if (sizeof(size_t) > sizeof(uint32_t))
    {
        const size_t too_many =
            static_cast<size_t>(0xffffffffu) + 1;
        RequireResult(leo2_encode_batch_preflight_scratch_size(
            fixture.codec, too_many, &bytes), LEO2_INVALID_ARGUMENT,
            "oversized encode batch preflight query");
        RequireResult(leo2_decode_plan_batch_preflight_scratch_size(
            fixture.plan, too_many, &bytes), LEO2_INVALID_ARGUMENT,
            "oversized decode batch preflight query");
    }

    const uint8_t all_original[3] = { 1, 1, 1 };
    const uint8_t all_recovery[2] = { 1, 1 };
    leo2_decode_plan* no_loss = NULL;
    RequireResult(leo2_decode_plan_create(fixture.codec,
        all_original, all_recovery, &no_loss), LEO2_SUCCESS,
        "scalable no-loss plan create");
    RequireResult(leo2_decode_plan_batch_preflight_scratch_size(
        no_loss, 1024, &bytes), LEO2_SUCCESS,
        "no-loss preflight scratch query");
    Require(bytes == 0, "no-loss batch requested preflight scratch");
    leo2_decode_batch_item junk[9];
    memset(junk, 0, sizeof(junk));
    for (size_t i = 0; i < 9; ++i)
    {
        junk[i].shard_bytes = 0;
        junk[i].original = reinterpret_cast<const void* const*>(
            static_cast<uintptr_t>(1));
        junk[i].recovery = reinterpret_cast<const void* const*>(
            static_cast<uintptr_t>(1));
        junk[i].restored_original = reinterpret_cast<void* const*>(
            static_cast<uintptr_t>(1));
        junk[i].scratch = reinterpret_cast<void*>(
            static_cast<uintptr_t>(1));
        junk[i].scratch_bytes = 1;
    }
    const leo2_result no_loss_result =
        leo2_decode_plan_execute_batch_with_preflight_scratch(
            no_loss, junk, 9, NULL, 0);
    leo2_decode_plan_destroy(no_loss);
    RequireResult(no_loss_result, LEO2_SUCCESS,
        "scalable no-loss batch inspected item state");
}

void TestScalableValidAndAllocationFree(Fixture& fixture)
{
    const size_t compatibility_counts[2] = { 1, 8 };
    for (size_t i = 0; i < 2; ++i)
    {
        ScalableEncodeBatch small(fixture, compatibility_counts[i]);
        RequireResult(small.Execute(), LEO2_SUCCESS,
            "small scalable encode compatibility path");
        small.Check();
        ScalableDecodeBatch small_decode(fixture, small);
        RequireResult(small_decode.Execute(), LEO2_SUCCESS,
            "small scalable decode compatibility path");
        small_decode.Check();
    }

    ScalableEncodeBatch encoded(fixture, 64);
    RequireResult(encoded.Execute(), LEO2_SUCCESS,
        "scalable encode batch");
    encoded.Check();

    BeginAllocationAudit();
    const leo2_result encode_result = encoded.Execute();
    const uint64_t encode_allocations = EndAllocationAudit();
    RequireResult(encode_result, LEO2_SUCCESS,
        "allocation-audited scalable encode");
    Require(encode_allocations == 0,
        "scalable encode preflight allocated memory");

    ScalableDecodeBatch decoded(fixture, encoded);
    RequireResult(decoded.Execute(), LEO2_SUCCESS,
        "scalable decode batch");
    decoded.Check();
    BeginAllocationAudit();
    const leo2_result decode_result = decoded.Execute();
    const uint64_t decode_allocations = EndAllocationAudit();
    RequireResult(decode_result, LEO2_SUCCESS,
        "allocation-audited scalable decode");
    Require(decode_allocations == 0,
        "scalable decode preflight allocated memory");
}

void TestScalableLargeBatch(Fixture& fixture)
{
    ScalableEncodeBatch encoded(fixture, 1024);
    RequireResult(encoded.Execute(), LEO2_SUCCESS,
        "1024-item scalable encode batch");
    encoded.Check();
    ScalableDecodeBatch decoded(fixture, encoded);
    RequireResult(decoded.Execute(), LEO2_SUCCESS,
        "1024-item scalable decode batch");
    decoded.Check();
}

void TestScalableEncodeConflicts(Fixture& fixture)
{
    {
        ScalableEncodeBatch batch(fixture, 9);
        batch.original_pointers[7] = fixture.original_b;
        batch.items[7].original = &batch.original_pointers[7][0];
        batch.recovery_pointers[0][0] = &fixture.source_b[0][16];
        const Shards before = fixture.source_b;
        RequireResult(batch.Execute(), LEO2_OVERLAP,
            "scalable encode output/input overlap");
        Require(fixture.source_b == before,
            "rejected scalable encode changed shared input");
    }
    {
        ScalableEncodeBatch batch(fixture, 9);
        Bytes collision(64, 0x5a);
        const Bytes before = collision;
        batch.recovery_pointers[0][0] = &collision[0];
        batch.recovery_pointers[1][0] = &collision[16];
        RequireResult(batch.Execute(), LEO2_OVERLAP,
            "scalable encode partial output/output overlap");
        Require(collision == before,
            "rejected scalable encode changed colliding outputs");
        batch.recovery_pointers[1][0] = batch.recovery_pointers[0][0];
        RequireResult(batch.Execute(), LEO2_OVERLAP,
            "scalable encode duplicate outputs");
    }
    {
        ScalableEncodeBatch batch(fixture, 9);
        AlignedBuffer shared(std::max(fixture.encode_short,
            fixture.encode_long) + leo2_scratch_alignment());
        batch.items[0].scratch = shared.data();
        batch.items[0].scratch_bytes = fixture.encode_short;
        batch.items[1].scratch = static_cast<uint8_t*>(shared.data()) +
            leo2_scratch_alignment();
        batch.items[1].scratch_bytes = fixture.encode_long;
        RequireResult(batch.Execute(), LEO2_OVERLAP,
            "scalable encode partial scratch overlap");
    }
    {
        ScalableEncodeBatch batch(fixture, 9);
        const std::vector<const void*> before = batch.original_pointers[8];
        batch.recovery_pointers[0][0] =
            const_cast<void*>(static_cast<const void*>(
                batch.original_pointers[8].data()));
        RequireResult(batch.Execute(), LEO2_OVERLAP,
            "scalable encode output/metadata overlap");
        Require(batch.original_pointers[8] == before,
            "rejected scalable encode changed metadata");
    }
    {
        /* A long immutable range contains a short immutable range, while a
           later writable range overlaps only the long one.  An adjacent-only
           sorted check misses this shape; the max-end sweep must reject it. */
        ScalableEncodeBatch batch(fixture, 9);
        Bytes backing(64, 0x3c);
        std::swap(batch.items[0].scratch, batch.items[1].scratch);
        std::swap(batch.items[0].scratch_bytes,
            batch.items[1].scratch_bytes);
        batch.items[0].shard_bytes = Fixture::kLongBytes;
        batch.items[1].shard_bytes = Fixture::kShortBytes;
        batch.original_pointers[0][0] = &backing[0];
        batch.original_pointers[1][0] = &backing[8];
        batch.recovery_pointers[2][0] = &backing[26];
        const Bytes before = backing;
        RequireResult(batch.Execute(), LEO2_OVERLAP,
            "scalable encode containing-input overlap");
        Require(backing == before,
            "containing-input rejection changed bytes");
    }
    {
        ScalableEncodeBatch batch(fixture, 9);
        Require(batch.preflight_bytes > 1,
            "scalable encode preflight unexpectedly tiny");
        RequireResult(leo2_encode_batch_with_preflight_scratch(
            fixture.codec, &batch.items[0], batch.items.size(), NULL,
            batch.preflight_bytes), LEO2_INVALID_ARGUMENT,
            "null scalable encode preflight scratch");
        RequireResult(leo2_encode_batch_with_preflight_scratch(
            fixture.codec, &batch.items[0], batch.items.size(),
            batch.preflight->data(), batch.preflight_bytes - 1),
            LEO2_SCRATCH_TOO_SMALL,
            "short scalable encode preflight scratch");
        AlignedBuffer unaligned(batch.preflight_bytes +
            leo2_scratch_alignment());
        RequireResult(leo2_encode_batch_with_preflight_scratch(
            fixture.codec, &batch.items[0], batch.items.size(),
            static_cast<uint8_t*>(unaligned.data()) + 1,
            batch.preflight_bytes), LEO2_BAD_ALIGNMENT,
            "unaligned scalable encode preflight scratch");
        const uintptr_t overflowing_workspace =
            (UINTPTR_MAX - batch.preflight_bytes / 2) &
            ~(static_cast<uintptr_t>(leo2_scratch_alignment()) - 1);
        RequireResult(leo2_encode_batch_with_preflight_scratch(
            fixture.codec, &batch.items[0], batch.items.size(),
            reinterpret_cast<void*>(overflowing_workspace),
            batch.preflight_bytes), LEO2_INVALID_ARGUMENT,
            "overflowing scalable encode preflight span");
        const uintptr_t overflowing_items = UINTPTR_MAX -
            sizeof(leo2_encode_batch_item) * 4;
        RequireResult(leo2_encode_batch_with_preflight_scratch(
            fixture.codec,
            reinterpret_cast<const leo2_encode_batch_item*>(
                overflowing_items),
            batch.items.size(), batch.preflight->data(),
            batch.preflight->size()), LEO2_INVALID_ARGUMENT,
            "overflowing scalable encode item span");
    }
    {
        ScalableEncodeBatch batch(fixture, 9);
        AlignedBuffer shared(std::max(
            batch.preflight_bytes, fixture.encode_short));
        memset(shared.data(), 0x91, shared.size());
        Bytes before(shared.size());
        memcpy(&before[0], shared.data(), shared.size());
        batch.items[0].scratch = shared.data();
        batch.items[0].scratch_bytes = fixture.encode_short;
        RequireResult(leo2_encode_batch_with_preflight_scratch(
            fixture.codec, &batch.items[0], batch.items.size(),
            shared.data(), batch.preflight_bytes), LEO2_OVERLAP,
            "scalable encode preflight/item-scratch overlap");
        Require(memcmp(shared.data(), &before[0], shared.size()) == 0,
            "workspace-overlap rejection modified workspace metadata");
    }
    {
        ScalableEncodeBatch batch(fixture, 9);
        AlignedBuffer shared(batch.preflight_bytes);
        const void** metadata = static_cast<const void**>(shared.data());
        for (size_t i = 0; i < fixture.original_a.size(); ++i)
            metadata[i] = fixture.original_a[i];
        const void* before[3] = {
            metadata[0], metadata[1], metadata[2]
        };
        batch.items[0].original = metadata;
        RequireResult(leo2_encode_batch_with_preflight_scratch(
            fixture.codec, &batch.items[0], batch.items.size(),
            shared.data(), shared.size()), LEO2_OVERLAP,
            "scalable encode preflight/metadata overlap");
        Require(metadata[0] == before[0] && metadata[1] == before[1] &&
                metadata[2] == before[2],
            "preflight/metadata rejection changed metadata");
    }
    {
        ScalableEncodeBatch batch(fixture, 9);
        const Bytes first_before = batch.outputs[0][0];
        Bytes preflight_before(batch.preflight->size());
        memcpy(&preflight_before[0], batch.preflight->data(),
            batch.preflight->size());
        --batch.items[8].scratch_bytes;
        RequireResult(batch.Execute(), LEO2_SCRATCH_TOO_SMALL,
            "scalable encode later-item validation failure");
        Require(batch.outputs[0][0] == first_before,
            "scalable encode ran an earlier item before rejection");
        Require(memcmp(batch.preflight->data(), &preflight_before[0],
                    preflight_before.size()) == 0,
            "read-only validation failure modified preflight scratch");
    }
}

void TestScalableDecodeConflicts(Fixture& fixture)
{
    ScalableEncodeBatch encoded(fixture, 9);
    RequireResult(encoded.Execute(), LEO2_SUCCESS,
        "decode-conflict source encode");
    {
        ScalableDecodeBatch batch(fixture, encoded);
        batch.original_pointers[7][1] = fixture.original_b[1];
        batch.original_pointers[7][2] = fixture.original_b[2];
        batch.restored_pointers[0][0] = &fixture.source_b[1][16];
        const Shards before = fixture.source_b;
        RequireResult(batch.Execute(), LEO2_OVERLAP,
            "scalable decode output/input overlap");
        Require(fixture.source_b == before,
            "rejected scalable decode changed input");
    }
    {
        ScalableDecodeBatch batch(fixture, encoded);
        Bytes collision(64, 0x6d);
        const Bytes before = collision;
        batch.restored_pointers[0][0] = &collision[0];
        batch.restored_pointers[1][0] = &collision[16];
        RequireResult(batch.Execute(), LEO2_OVERLAP,
            "scalable decode partial output/output overlap");
        Require(collision == before,
            "rejected scalable decode changed outputs");
    }
    {
        ScalableDecodeBatch batch(fixture, encoded);
        batch.restored_pointers[0][0] =
            const_cast<void*>(static_cast<const void*>(
                batch.original_pointers[8].data()));
        const std::vector<const void*> before = batch.original_pointers[8];
        RequireResult(batch.Execute(), LEO2_OVERLAP,
            "scalable decode output/metadata overlap");
        Require(batch.original_pointers[8] == before,
            "rejected scalable decode changed metadata");
    }
    {
        ScalableDecodeBatch batch(fixture, encoded);
        AlignedBuffer shared(batch.preflight_bytes);
        memset(shared.data(), 0x73, shared.size());
        Bytes before(shared.size());
        memcpy(&before[0], shared.data(), shared.size());
        batch.restored_pointers[0][0] = shared.data();
        RequireResult(
            leo2_decode_plan_execute_batch_with_preflight_scratch(
                fixture.plan, &batch.items[0], batch.items.size(),
                shared.data(), shared.size()), LEO2_OVERLAP,
            "scalable decode preflight/output overlap");
        Require(memcmp(shared.data(), &before[0], shared.size()) == 0,
            "decode workspace-overlap rejection modified bytes");
    }
    {
        ScalableDecodeBatch batch(fixture, encoded);
        const Bytes first_before = batch.restored[0];
        --batch.items[8].scratch_bytes;
        RequireResult(batch.Execute(), LEO2_SCRATCH_TOO_SMALL,
            "scalable decode later-item validation failure");
        Require(batch.restored[0] == first_before,
            "scalable decode ran an earlier item before rejection");
    }
}

void TestConcurrentScalableCalls(Fixture& fixture)
{
    static const size_t kCallers = 4;
    std::atomic<size_t> ready(0);
    std::atomic<bool> start(false);
    std::vector<std::string> errors(kCallers);
    std::vector<std::thread> threads;
    for (size_t caller = 0; caller < kCallers; ++caller)
    {
        threads.push_back(std::thread([&, caller]() {
            try
            {
                ScalableEncodeBatch encoded(fixture, 9);
                ScalableDecodeBatch decoded(fixture, encoded);
                ready.fetch_add(1, std::memory_order_release);
                while (!start.load(std::memory_order_acquire))
                    std::this_thread::yield();
                for (unsigned repetition = 0; repetition < 4; ++repetition)
                {
                    RequireResult(encoded.Execute(), LEO2_SUCCESS,
                        "concurrent scalable encode");
                    encoded.Check();
                    RequireResult(decoded.Execute(), LEO2_SUCCESS,
                        "concurrent scalable decode");
                    decoded.Check();
                }
            }
            catch (const std::exception& error)
            {
                errors[caller] = error.what();
            }
        }));
    }
    while (ready.load(std::memory_order_acquire) != kCallers)
        std::this_thread::yield();
    start.store(true, std::memory_order_release);
    for (size_t i = 0; i < threads.size(); ++i)
        threads[i].join();
    for (size_t i = 0; i < errors.size(); ++i)
        Require(errors[i].empty(), errors[i]);
}

void TestScalableDifferential(Fixture& fixture)
{
    uint32_t state = 0x31415926u;
    for (unsigned trial = 0; trial < 128; ++trial)
    {
        ScalableEncodeBatch batch(fixture, 9);
        Bytes arena(640, 0x4a);
        for (size_t item = 0; item < batch.items.size(); ++item)
        {
            state = state * 1664525u + 1013904223u;
            const size_t offset = (trial & 1u) == 0
                ? item * 64 + ((state >> 24) & 7u)
                : (state >> 24) % 96;
            batch.recovery_pointers[item][0] = &arena[offset];
        }
        const Bytes initial = arena;
        const leo2_result compatibility = leo2_encode_batch(
            fixture.codec, &batch.items[0], batch.items.size());
        const Bytes compatibility_bytes = arena;
        arena = initial;
        const leo2_result scalable = batch.Execute();
        RequireResult(scalable, compatibility,
            "randomized scalable encode/fallback differential");
        if (scalable == LEO2_SUCCESS)
            Require(arena == compatibility_bytes,
                "randomized scalable encode bytes differ from fallback");
    }

    ScalableEncodeBatch encoded(fixture, 9);
    RequireResult(encoded.Execute(), LEO2_SUCCESS,
        "randomized decode differential source encode");
    for (unsigned trial = 0; trial < 128; ++trial)
    {
        ScalableDecodeBatch batch(fixture, encoded);
        Bytes arena(640, 0x6b);
        for (size_t item = 0; item < batch.items.size(); ++item)
        {
            state = state * 1664525u + 1013904223u;
            const size_t offset = (trial & 1u) == 0
                ? item * 64 + ((state >> 24) & 7u)
                : (state >> 24) % 96;
            batch.restored_pointers[item][0] = &arena[offset];
        }
        const Bytes initial = arena;
        const leo2_result compatibility =
            leo2_decode_plan_execute_batch(
                fixture.plan, &batch.items[0], batch.items.size());
        const Bytes compatibility_bytes = arena;
        arena = initial;
        const leo2_result scalable = batch.Execute();
        RequireResult(scalable, compatibility,
            "randomized scalable decode/fallback differential");
        if (scalable == LEO2_SUCCESS)
            Require(arena == compatibility_bytes,
                "randomized scalable decode bytes differ from fallback");
    }
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
    TestScalableQueriesAndNoOp(fixture);
    TestScalableValidAndAllocationFree(fixture);
    TestScalableEncodeConflicts(fixture);
    TestScalableDecodeConflicts(fixture);
    TestScalableDifferential(fixture);
    TestScalableLargeBatch(fixture);
    if (thread_count > 1)
        TestConcurrentScalableCalls(fixture);
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
                  << " scalable_batch_sizes=9,64,1024"
                  << " allocation_audit=clean concurrent_scalable_calls=32"
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
