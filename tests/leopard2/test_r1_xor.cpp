/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"
#include "leopard.h"
#include "leopard2.h"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <thread>
#include <vector>

namespace {

typedef std::vector<uint8_t> Bytes;

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void require_result(leo2_result actual, leo2_result expected,
    const char* message)
{
    if (actual != expected)
    {
        char detail[256];
        std::snprintf(detail, sizeof(detail), "%s: got %s (%d), expected %s (%d)",
            message, leo2_result_string(actual), static_cast<int>(actual),
            leo2_result_string(expected), static_cast<int>(expected));
        throw std::runtime_error(detail);
    }
}

uint8_t pattern(uint64_t index, uint32_t salt)
{
    uint64_t value = index + UINT64_C(0x9e3779b97f4a7c15) * (salt + 1U);
    value ^= value >> 30;
    value *= UINT64_C(0xbf58476d1ce4e5b9);
    value ^= value >> 27;
    value *= UINT64_C(0x94d049bb133111eb);
    value ^= value >> 31;
    return static_cast<uint8_t>(value);
}

void fill(Bytes& bytes, uint32_t salt)
{
    for (size_t i = 0; i < bytes.size(); ++i)
        bytes[i] = pattern(i, salt);
}

std::vector<uint64_t> byte_counts()
{
    std::vector<uint64_t> result;
    for (uint64_t bytes = 0; bytes <= 257; ++bytes)
        result.push_back(bytes);
    const uint64_t boundaries[] = {
        511, 512, 513, 1023, 1024, 1025,
        4095, 4096, 4097, 65535, 65536, 65537,
        1048575, 1048576, 1048579
    };
    result.insert(result.end(), boundaries,
        boundaries + sizeof(boundaries) / sizeof(boundaries[0]));
    return result;
}

void test_primitive_case(
    const leopard::backend::Ops& ops,
    uint64_t byte_count,
    uint32_t salt)
{
    const size_t guard = 96;
    const size_t output_offset = guard + (salt * 7U + byte_count) % 47U;
    const size_t source0_offset = guard + (salt * 11U + byte_count) % 53U;
    const size_t source1_offset = guard + (salt * 13U + byte_count) % 59U;
    const size_t size = static_cast<size_t>(byte_count) + guard * 2U + 64U;
    Bytes output(size);
    Bytes source0(size);
    Bytes source1(size);
    fill(output, salt + 1U);
    fill(source0, salt + 2U);
    fill(source1, salt + 3U);
    const Bytes original_source0 = source0;
    const Bytes original_source1 = source1;
    Bytes expected = output;
    for (uint64_t i = 0; i < byte_count; ++i)
        expected[output_offset + static_cast<size_t>(i)] ^=
            source0[source0_offset + static_cast<size_t>(i)] ^
            source1[source1_offset + static_cast<size_t>(i)];
    ops.xor_memory_2to1(
        &output[output_offset],
        &source0[source0_offset],
        &source1[source1_offset], byte_count);
    require(output == expected, "fused XOR value, tail, or guard mismatch");
    require(source0 == original_source0 && source1 == original_source1,
        "fused XOR modified a read-only input");

    // Exact source aliasing is permitted by the public shard contract.  The
    // duplicated source cancels, so every destination byte must be unchanged.
    fill(output, salt + 4U);
    expected = output;
    ops.xor_memory_2to1(
        &output[output_offset],
        &source0[source0_offset],
        &source0[source0_offset], byte_count);
    require(output == expected, "identical aliased XOR sources did not cancel");

    // The contract permits overlapping read-only input shards as well as exact
    // aliases.  Exercise a one-byte displacement without overlapping output.
    Bytes shared(size + 1U);
    fill(shared, salt + 5U);
    fill(output, salt + 6U);
    expected = output;
    for (uint64_t i = 0; i < byte_count; ++i)
        expected[output_offset + static_cast<size_t>(i)] ^=
            shared[source0_offset + static_cast<size_t>(i)] ^
            shared[source0_offset + static_cast<size_t>(i) + 1U];
    const Bytes original_shared = shared;
    ops.xor_memory_2to1(
        &output[output_offset],
        &shared[source0_offset],
        &shared[source0_offset + 1U], byte_count);
    require(output == expected, "overlapping read-only XOR sources mismatch");
    require(shared == original_shared,
        "overlapping read-only XOR sources were modified");
}

void test_coarse_primitive_case(
    const leopard::backend::Ops& ops,
    uint64_t byte_count,
    uint32_t source_count,
    bool sparse,
    uint32_t salt)
{
    const size_t guard = 96;
    const size_t size = static_cast<size_t>(byte_count) + guard * 2U + 64U;
    const size_t output_offset = guard + (salt * 7U + source_count) % 31U;
    const size_t initial_offset = guard + (salt * 11U + source_count) % 29U;
    const size_t shared_offset = guard + (salt * 13U + source_count) % 27U;
    const size_t other_offset = guard + (salt * 17U + source_count) % 25U;

    Bytes output(size);
    Bytes initial(size);
    Bytes shared(size + 1U);
    Bytes other(size);
    fill(output, salt + 1U);
    fill(initial, salt + 2U);
    fill(shared, salt + 3U);
    fill(other, salt + 4U);
    const Bytes initial_before = initial;
    const Bytes shared_before = shared;
    const Bytes other_before = other;

    std::vector<const void*> sources(source_count, NULL);
    for (uint32_t i = 0; i < source_count; ++i)
    {
        if (sparse && ((i + salt) % 4U) == 0)
            continue;
        switch ((i + salt) % 7U)
        {
        case 0:
            // The initial source may also appear in the reduction list.
            sources[i] = &initial[initial_offset];
            break;
        case 1:
        case 2:
            // Exact read-only aliases exercise even/odd cancellation.
            sources[i] = &shared[shared_offset];
            break;
        case 3:
            // Partially overlapping read-only ranges are permitted.
            sources[i] = &shared[shared_offset + 1U];
            break;
        default:
            sources[i] = &other[other_offset];
            break;
        }
    }

    Bytes expected = output;
    for (uint64_t byte = 0; byte < byte_count; ++byte)
    {
        uint8_t value = initial[initial_offset + static_cast<size_t>(byte)];
        for (uint32_t source = 0; source < source_count; ++source)
            if (sources[source])
                value ^= static_cast<const uint8_t*>(
                    sources[source])[static_cast<size_t>(byte)];
        expected[output_offset + static_cast<size_t>(byte)] = value;
    }

    ops.xor_memory_sources(
        &output[output_offset], &initial[initial_offset],
        sources.empty() ? NULL : &sources[0], source_count, byte_count);
    require(output == expected,
        "coarse XOR value, source remainder, tail, or guard mismatch");
    require(initial == initial_before && shared == shared_before &&
            other == other_before,
        "coarse XOR modified a read-only input");
}

void test_primitive(const leopard::backend::Ops& ops)
{
    const std::vector<uint64_t> counts = byte_counts();
    for (size_t i = 0; i < counts.size(); ++i)
        test_primitive_case(ops, counts[i],
            static_cast<uint32_t>(i * 17U + ops.kind * 101U));

    static const uint64_t coarse_byte_counts[] = {
        0, 1, 2, 3, 7, 8, 15, 16, 17, 31, 32, 33,
        63, 64, 65, 127, 128, 129, 255, 256, 257,
        1023, 1024, 1025, 4095, 4096, 4097, 4113, 65537, 65553
    };
    static const uint32_t coarse_source_counts[] = {
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
        12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
        31, 32, 33
    };
    for (size_t byte_i = 0;
         byte_i < sizeof(coarse_byte_counts) / sizeof(coarse_byte_counts[0]);
         ++byte_i)
        for (size_t source_i = 0;
             source_i < sizeof(coarse_source_counts) /
                 sizeof(coarse_source_counts[0]);
             ++source_i)
            for (unsigned sparse = 0; sparse < 2; ++sparse)
                test_coarse_primitive_case(ops,
                    coarse_byte_counts[byte_i],
                    coarse_source_counts[source_i], sparse != 0,
                    static_cast<uint32_t>(
                        byte_i * 977U + source_i * 131U +
                        sparse * 17U + ops.kind * 1009U));

    // Exercise many complete backend groups without making every large-tail
    // case allocate a comparably large pointer list.
    test_coarse_primitive_case(
        ops, 257, 257, true, 0x243f6a88U + ops.kind);
}

void test_primitive_concurrency(const leopard::backend::Ops& ops)
{
    static const unsigned thread_count = 16;
    static const unsigned rounds = 48;
    std::atomic<bool> failed(false);
    std::vector<std::thread> threads;
    for (unsigned thread = 0; thread < thread_count; ++thread)
    {
        threads.push_back(std::thread([&ops, &failed, thread]() {
            try
            {
                for (unsigned round = 0; round < rounds; ++round)
                {
                    const uint64_t count = static_cast<uint64_t>(
                        (thread * 977U + round * 131U) % 8193U);
                    test_primitive_case(ops, count,
                        thread * 1009U + round * 37U);
                    test_coarse_primitive_case(ops, count,
                        (thread * 17U + round * 11U) % 41U,
                        (thread + round) % 3U == 0,
                        thread * 313U + round * 101U + ops.kind * 19U);
                }
            }
            catch (...)
            {
                failed.store(true, std::memory_order_relaxed);
            }
        }));
    }
    for (size_t i = 0; i < threads.size(); ++i)
        threads[i].join();
    require(!failed.load(std::memory_order_relaxed),
        "concurrent fused XOR execution failed");
}

class AlignedScratch
{
public:
    explicit AlignedScratch(size_t bytes)
        : storage_(bytes + leo2_scratch_alignment())
        , data_(NULL)
    {
        if (bytes == 0)
            return;
        const uintptr_t value = reinterpret_cast<uintptr_t>(&storage_[0]);
        const uintptr_t alignment = leo2_scratch_alignment();
        data_ = reinterpret_cast<void*>(
            (value + alignment - 1U) & ~(alignment - 1U));
    }

    void* data() { return data_; }

private:
    Bytes storage_;
    void* data_;
};

struct R1Fixture
{
    leo2_context* context;
    leo2_codec* codec;
    leo2_decode_plan* plan;
    uint32_t k;
    uint32_t missing;
    size_t bytes;
    std::vector<Bytes> storage;
    std::vector<const void*> original;
    Bytes recovery_storage;
    const void* recovery[1];

    R1Fixture(leo2_backend backend, uint32_t original_count,
        size_t shard_bytes, bool alias_inputs, leo2_field field,
        leo2_shard_layout layout, uint32_t missing_index = 0,
        uint32_t thread_count = 0)
        : context(NULL)
        , codec(NULL)
        , plan(NULL)
        , k(original_count)
        , missing(missing_index)
        , bytes(shard_bytes)
        , storage(original_count)
        , original(original_count, NULL)
        , recovery_storage(shard_bytes + 17U)
    {
        require(missing < k, "R=1 missing index is out of range");
        leo2_context_options options;
        std::memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = backend;
        options.thread_count = thread_count;
        require_result(leo2_context_create(&options, &context), LEO2_SUCCESS,
            "R=1 context create");
        leo2_codec_options codec_options;
        std::memset(&codec_options, 0, sizeof(codec_options));
        codec_options.struct_size = sizeof(codec_options);
        codec_options.shard_layout = layout;
        require_result(leo2_codec_create(context, k, 1,
            LEO2_PROFILE_LEGACY_HIGH_V1, field, &codec_options, &codec),
            LEO2_SUCCESS, "R=1 codec create");
        require(leo2_codec_field(codec) == field,
            "R=1 codec selected the wrong explicit field");
        require(leo2_codec_shard_layout(codec) == layout,
            "R=1 codec selected the wrong shard layout");
        for (uint32_t shard = 0; shard < k; ++shard)
        {
            const size_t offset = 1U + shard % 13U;
            storage[shard].resize(bytes + offset + 8U);
            fill(storage[shard], shard * 97U + static_cast<uint32_t>(bytes));
            original[shard] = &storage[shard][offset];
            if (layout == LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1)
                storage[shard][offset + bytes - 1U] = 0;
        }
        if (alias_inputs && k >= 3)
            original[1] = original[0];

        size_t scratch_bytes = 0;
        require_result(leo2_encode_scratch_size(codec, bytes, &scratch_bytes),
            LEO2_SUCCESS, "R=1 encode scratch query");
        AlignedScratch scratch(scratch_bytes);
        void* output[1] = { &recovery_storage[3] };
        require_result(leo2_encode(codec, bytes, &original[0], output,
            scratch.data(), scratch_bytes), LEO2_SUCCESS, "R=1 encode");
        recovery[0] = output[0];

        Bytes expected(bytes, 0);
        for (uint32_t shard = 0; shard < k; ++shard)
        {
            const uint8_t* input = static_cast<const uint8_t*>(original[shard]);
            for (size_t i = 0; i < bytes; ++i)
                expected[i] ^= input[i];
        }
        require(std::memcmp(recovery[0], &expected[0], bytes) == 0,
            "R=1 encode differs from independent XOR");

        if ((bytes & 63U) == 0)
        {
            const unsigned legacy_work_count =
                leo_encode_work_count(k, 1);
            require(legacy_work_count != 0,
                "legacy R=1 work-count query failed");
            std::vector<Bytes> legacy_storage(
                legacy_work_count, Bytes(bytes, 0xa5));
            std::vector<void*> legacy_work(legacy_work_count, NULL);
            for (unsigned i = 0; i < legacy_work_count; ++i)
                legacy_work[i] = &legacy_storage[i][0];
            require(leo_encode(bytes, k, 1, legacy_work_count,
                        &original[0], &legacy_work[0]) == Leopard_Success,
                "legacy R=1 encode failed");
            require(std::memcmp(
                        recovery[0], legacy_work[0], bytes) == 0,
                "R=1 encode differs from legacy Leopard wire bytes");
        }

        std::vector<uint8_t> original_present(k, 1);
        const uint8_t recovery_present[1] = { 1 };
        original_present[missing] = 0;
        require_result(leo2_decode_plan_create(codec, &original_present[0],
            recovery_present, &plan), LEO2_SUCCESS, "R=1 plan create");
    }

    ~R1Fixture()
    {
        leo2_decode_plan_destroy(plan);
        leo2_codec_destroy(codec);
        leo2_context_destroy(context);
    }

private:
    R1Fixture(const R1Fixture&);
    R1Fixture& operator=(const R1Fixture&);
};

void execute_and_check_decode(const R1Fixture& fixture)
{
    size_t scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(
        fixture.plan, fixture.bytes, &scratch_bytes), LEO2_SUCCESS,
        "R=1 decode scratch query");
    AlignedScratch scratch(scratch_bytes);
    std::vector<const void*> received = fixture.original;
    received[fixture.missing] = NULL;
    Bytes restored_storage(fixture.bytes + 11U, 0xa5);
    std::vector<void*> restored(fixture.k, NULL);
    restored[fixture.missing] = &restored_storage[5];
    require_result(leo2_decode_plan_execute(fixture.plan, fixture.bytes,
        &received[0], fixture.recovery, &restored[0], scratch.data(),
        scratch_bytes), LEO2_SUCCESS, "R=1 decode execute");
    require(std::memcmp(restored[fixture.missing],
                fixture.original[fixture.missing], fixture.bytes) == 0,
        "R=1 decode restored the wrong original");
    for (size_t i = 0; i < 5; ++i)
        require(restored_storage[i] == 0xa5, "R=1 decode changed a prefix guard");
    for (size_t i = 5 + fixture.bytes; i < restored_storage.size(); ++i)
        require(restored_storage[i] == 0xa5, "R=1 decode changed a suffix guard");

    std::fill(restored_storage.begin(), restored_storage.end(), 0x5a);
    leo2_decode_batch_item item = {
        fixture.bytes, &received[0], fixture.recovery, &restored[0],
        scratch.data(), scratch_bytes
    };
    require_result(leo2_decode_plan_execute_batch(fixture.plan, &item, 1),
        LEO2_SUCCESS, "R=1 single-item batch decode execute");
    require(std::memcmp(restored[fixture.missing],
                fixture.original[fixture.missing], fixture.bytes) == 0,
        "R=1 single-item batch restored the wrong original");
    for (size_t i = 0; i < 5; ++i)
        require(restored_storage[i] == 0x5a,
            "R=1 single-item batch changed a prefix guard");
    for (size_t i = 5 + fixture.bytes; i < restored_storage.size(); ++i)
        require(restored_storage[i] == 0x5a,
            "R=1 single-item batch changed a suffix guard");

    std::fill(restored_storage.begin(), restored_storage.end(), 0x3c);
    std::vector<uint8_t> original_present(fixture.k, 1);
    original_present[fixture.missing] = 0;
    const uint8_t recovery_present[1] = { 1 };
    require_result(leo2_decode(fixture.codec, fixture.bytes,
        &original_present[0], recovery_present, &received[0],
        fixture.recovery, &restored[0], scratch.data(), scratch_bytes),
        LEO2_SUCCESS, "R=1 one-shot decode execute");
    require(std::memcmp(restored[fixture.missing],
                fixture.original[fixture.missing], fixture.bytes) == 0,
        "R=1 one-shot decode restored the wrong original");
    for (size_t i = 0; i < 5; ++i)
        require(restored_storage[i] == 0x3c,
            "R=1 one-shot decode changed a prefix guard");
    for (size_t i = 5 + fixture.bytes; i < restored_storage.size(); ++i)
        require(restored_storage[i] == 0x3c,
            "R=1 one-shot decode changed a suffix guard");
}

void execute_public_r1_multi_item_batch(
    const R1Fixture& fixture, size_t batch_count, bool adversarial)
{
    static const size_t kOutputOffset = 5;
    size_t scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(
        fixture.plan, fixture.bytes, &scratch_bytes), LEO2_SUCCESS,
        "R=1 multi-item decode scratch query");

    std::vector<std::vector<const void*> > received(
        batch_count, fixture.original);
    std::vector<Bytes> restored_storage(
        batch_count, Bytes(fixture.bytes + 11U, 0x5a));
    std::vector<std::vector<void*> > restored(
        batch_count, std::vector<void*>(fixture.k, NULL));
    std::vector<std::unique_ptr<AlignedScratch> > scratch;
    std::vector<leo2_decode_batch_item> items(batch_count);
    scratch.reserve(batch_count);
    for (size_t item_i = 0; item_i < batch_count; ++item_i)
    {
        received[item_i][fixture.missing] = NULL;
        restored[item_i][fixture.missing] =
            &restored_storage[item_i][kOutputOffset];
        scratch.push_back(std::unique_ptr<AlignedScratch>(
            new AlignedScratch(scratch_bytes)));
        items[item_i].shard_bytes = fixture.bytes;
        items[item_i].original = &received[item_i][0];
        items[item_i].recovery = fixture.recovery;
        items[item_i].restored_original = &restored[item_i][0];
        items[item_i].scratch = scratch[item_i]->data();
        items[item_i].scratch_bytes = scratch_bytes;
    }

    require_result(leo2_decode_plan_execute_batch(
        fixture.plan, &items[0], items.size()), LEO2_SUCCESS,
        "R=1 multi-item batch decode execute");
    const auto check_outputs = [&]()
    {
        for (size_t item_i = 0; item_i < batch_count; ++item_i)
        {
            require(std::memcmp(restored[item_i][fixture.missing],
                        fixture.original[fixture.missing], fixture.bytes) == 0,
                "R=1 multi-item batch restored the wrong original");
            for (size_t i = 0; i < kOutputOffset; ++i)
                require(restored_storage[item_i][i] == 0x5a,
                    "R=1 multi-item batch changed an output prefix guard");
            for (size_t i = kOutputOffset + fixture.bytes;
                 i < restored_storage[item_i].size(); ++i)
                require(restored_storage[item_i][i] == 0x5a,
                    "R=1 multi-item batch changed an output suffix guard");
        }
    };
    check_outputs();

    /* Two or more items may use the allocation-free, caller-supplied interval
       preflight.  Legacy-high K=1,R=1 retains its measured nine-item cutoff.
       Compact R=1 plans omit all presence vectors, so exercise every compact
       policy boundary through that independent entry point instead of relying
       on the compatibility preflight above to cover reconstructed metadata. */
    const size_t scalable_minimum = fixture.k == 1 ? 9 : 2;
    if (batch_count >= scalable_minimum &&
        (fixture.k == 1 || fixture.k == 3 || fixture.k == 5 ||
         fixture.k == 6 || fixture.k == 9) && fixture.bytes == 4097)
    {
        size_t preflight_bytes = 0;
        require_result(leo2_decode_plan_batch_preflight_scratch_size(
            fixture.plan, batch_count, &preflight_bytes), LEO2_SUCCESS,
            "R=1 compact scalable preflight scratch query");
        require(preflight_bytes != 0,
            "R=1 compact scalable preflight scratch mismatch");
        AlignedScratch preflight(preflight_bytes);
        for (size_t item_i = 0; item_i < batch_count; ++item_i)
            std::fill(restored_storage[item_i].begin(),
                restored_storage[item_i].end(), 0x5a);
        require_result(leo2_decode_plan_execute_batch_with_preflight_scratch(
            fixture.plan, &items[0], items.size(), preflight.data(),
            preflight_bytes), LEO2_SUCCESS,
            "R=1 compact scalable batch decode execute");
        check_outputs();

        if (fixture.k == 1)
        {
            AlignedScratch optional_item_scratch(64);
            items[0].scratch = optional_item_scratch.data();
            items[0].scratch_bytes = 64;
            for (size_t item_i = 0; item_i < batch_count; ++item_i)
                std::fill(restored_storage[item_i].begin(),
                    restored_storage[item_i].end(), 0x5a);
            require_result(
                leo2_decode_plan_execute_batch_with_preflight_scratch(
                    fixture.plan, &items[0], items.size(), preflight.data(),
                    preflight_bytes), LEO2_SUCCESS,
                "K=1 scalable batch optional item scratch execute");
            check_outputs();
            items[0].scratch = NULL;
            items[0].scratch_bytes = 0;

            for (size_t item_i = 0; item_i < batch_count; ++item_i)
                std::fill(restored_storage[item_i].begin(),
                    restored_storage[item_i].end(), 0x6b);
            const std::vector<Bytes> scalable_before_failure =
                restored_storage;
            const void* const* const saved_recovery = items.back().recovery;
            items.back().recovery = NULL;
            require_result(
                leo2_decode_plan_execute_batch_with_preflight_scratch(
                    fixture.plan, &items[0], items.size(), preflight.data(),
                    preflight_bytes), LEO2_INVALID_ARGUMENT,
                "K=1 scalable batch late-item pointer rejection");
            require(restored_storage == scalable_before_failure,
                "K=1 scalable batch ran work before preflight failure");
            items.back().recovery = saved_recovery;
        }
    }

    for (size_t item_i = 0; item_i < batch_count; ++item_i)
        std::fill(restored_storage[item_i].begin(),
            restored_storage[item_i].end(), 0x5a);
    std::unique_ptr<AlignedScratch> binding_optional_scratch;
    if (batch_count == 1)
    {
        binding_optional_scratch.reset(new AlignedScratch(64));
        items[0].scratch = binding_optional_scratch->data();
        items[0].scratch_bytes = 64;
    }
    leo2_decode_batch_binding* binding = NULL;
    require_result(leo2_decode_batch_binding_create(
        fixture.plan, &items[0], items.size(), &binding),
        LEO2_SUCCESS, "R=1 decode binding create");
    require(binding != NULL &&
            leo2_decode_batch_binding_item_count(binding) == items.size(),
        "R=1 decode binding item-count mismatch");
    const void* const* const saved_binding_recovery = items[0].recovery;
    void* const* const saved_binding_restored = items[0].restored_original;
    items[0].recovery = NULL;
    items[0].restored_original = NULL;
    require_result(leo2_decode_batch_binding_execute(binding),
        LEO2_SUCCESS, "R=1 decode binding captured-metadata execute");
    items[0].recovery = saved_binding_recovery;
    items[0].restored_original = saved_binding_restored;
    check_outputs();

    uint8_t* const live_recovery =
        const_cast<uint8_t*>(static_cast<const uint8_t*>(fixture.recovery[0]));
    const uint8_t saved_recovery_byte = live_recovery[0];
    Bytes changed_expected(
        static_cast<const uint8_t*>(fixture.original[fixture.missing]),
        static_cast<const uint8_t*>(fixture.original[fixture.missing]) +
            fixture.bytes);
    changed_expected[0] ^= 0x5bu;
    live_recovery[0] ^= 0x5bu;
    for (size_t item_i = 0; item_i < batch_count; ++item_i)
        std::fill(restored_storage[item_i].begin(),
            restored_storage[item_i].end(), 0x5a);
    require_result(leo2_decode_batch_binding_execute(binding),
        LEO2_SUCCESS, "R=1 decode binding changed-input execute");
    for (size_t item_i = 0; item_i < batch_count; ++item_i)
        require(std::memcmp(restored[item_i][fixture.missing],
                    &changed_expected[0], fixture.bytes) == 0,
            "R=1 decode binding ignored changed received bytes");
    live_recovery[0] = saved_recovery_byte;
    for (size_t item_i = 0; item_i < batch_count; ++item_i)
        std::fill(restored_storage[item_i].begin(),
            restored_storage[item_i].end(), 0x5a);
    require_result(leo2_decode_batch_binding_execute(binding),
        LEO2_SUCCESS, "R=1 decode binding restored-input execute");
    check_outputs();
    leo2_decode_batch_binding_destroy(binding);
    if (batch_count == 1)
    {
        items[0].scratch = NULL;
        items[0].scratch_bytes = 0;
    }

    if (!adversarial)
        return;

    for (size_t item_i = 0; item_i < batch_count; ++item_i)
        std::fill(restored_storage[item_i].begin(),
            restored_storage[item_i].end(), 0x6b);
    const std::vector<Bytes> before_failure = restored_storage;
    if (items.back().scratch_bytes != 0)
    {
        --items.back().scratch_bytes;
        require_result(leo2_decode_plan_execute_batch(
            fixture.plan, &items[0], items.size()), LEO2_SCRATCH_TOO_SMALL,
            "R=1 compact batch late-item scratch rejection");
        require(restored_storage == before_failure,
            "R=1 compact batch ran an item before late validation failure");
        ++items.back().scratch_bytes;
    }

    const void* const* const saved_recovery = items.back().recovery;
    items.back().recovery = NULL;
    require_result(leo2_decode_plan_execute_batch(
        fixture.plan, &items[0], items.size()), LEO2_INVALID_ARGUMENT,
        "R=1 compact batch late-item pointer rejection");
    require(restored_storage == before_failure,
        "R=1 compact batch ran an item before late pointer rejection");
    items.back().recovery = saved_recovery;

    void* const saved_output = restored[0][fixture.missing];
    restored[0][fixture.missing] =
        const_cast<void*>(static_cast<const void*>(&received[1][0]));
    require_result(leo2_decode_plan_execute_batch(
        fixture.plan, &items[0], items.size()), LEO2_OVERLAP,
        "R=1 compact batch cross-item output/metadata overlap");
    require(restored_storage == before_failure,
        "R=1 compact batch metadata rejection modified an output");
    restored[0][fixture.missing] = saved_output;

    const Bytes aliased_input_before = fixture.k == 1
        ? fixture.recovery_storage
        : fixture.storage[fixture.missing + 1U];
    restored[0][fixture.missing] = fixture.k == 1
        ? const_cast<void*>(fixture.recovery[0])
        : const_cast<void*>(received[1][fixture.missing + 1U]);
    require_result(leo2_decode_plan_execute_batch(
        fixture.plan, &items[0], items.size()), LEO2_OVERLAP,
        "R=1 compact batch cross-item output/input overlap");
    require((fixture.k == 1
            ? fixture.recovery_storage
            : fixture.storage[fixture.missing + 1U]) == aliased_input_before,
        "R=1 compact batch alias rejection modified an input");
    require(restored_storage == before_failure,
        "R=1 compact batch alias rejection modified an output");
}

void execute_public_k1_encode_batch(
    const R1Fixture& fixture, size_t batch_count)
{
    require(fixture.k == 1, "K=1 encode batch received another code");
    static const size_t kOutputOffset = 5;
    size_t scratch_bytes = 1;
    size_t preflight_bytes = 1;
    require_result(leo2_encode_scratch_size(
        fixture.codec, fixture.bytes, &scratch_bytes), LEO2_SUCCESS,
        "K=1 encode batch scratch query");
    require_result(leo2_encode_batch_preflight_scratch_size(
        fixture.codec, batch_count, &preflight_bytes), LEO2_SUCCESS,
        "K=1 encode batch preflight scratch query");
    require(scratch_bytes == 0 &&
            (batch_count < 9 ? preflight_bytes == 0 : preflight_bytes != 0),
        "K=1 encode batch scratch query mismatch");
    AlignedScratch preflight(preflight_bytes);

    std::vector<Bytes> output_storage(
        batch_count, Bytes(fixture.bytes + 11U, 0x5a));
    std::vector<std::vector<void*> > outputs(
        batch_count, std::vector<void*>(1, NULL));
    std::vector<leo2_encode_batch_item> items(batch_count);
    for (size_t item_i = 0; item_i < batch_count; ++item_i)
    {
        outputs[item_i][0] = &output_storage[item_i][kOutputOffset];
        items[item_i].shard_bytes = fixture.bytes;
        items[item_i].original = &fixture.original[0];
        items[item_i].recovery = &outputs[item_i][0];
        items[item_i].scratch = NULL;
        items[item_i].scratch_bytes = 0;
    }
    const auto check_outputs = [&]()
    {
        for (size_t item_i = 0; item_i < batch_count; ++item_i)
        {
            require(std::memcmp(outputs[item_i][0], fixture.recovery[0],
                        fixture.bytes) == 0,
                "K=1 encode batch produced the wrong parity");
            for (size_t i = 0; i < kOutputOffset; ++i)
                require(output_storage[item_i][i] == 0x5a,
                    "K=1 encode batch changed an output prefix guard");
            for (size_t i = kOutputOffset + fixture.bytes;
                 i < output_storage[item_i].size(); ++i)
                require(output_storage[item_i][i] == 0x5a,
                    "K=1 encode batch changed an output suffix guard");
        }
    };
    require_result(leo2_encode_batch(
        fixture.codec, &items[0], items.size()), LEO2_SUCCESS,
        "K=1 encode batch execute");
    check_outputs();

    for (size_t item_i = 0; item_i < batch_count; ++item_i)
        std::fill(output_storage[item_i].begin(),
            output_storage[item_i].end(), 0x5a);
    require_result(leo2_encode_batch_with_preflight_scratch(
        fixture.codec, &items[0], items.size(), preflight.data(),
        preflight_bytes), LEO2_SUCCESS,
        "K=1 scalable encode batch execute");
    check_outputs();

    std::unique_ptr<AlignedScratch> binding_optional_scratch;
    if (batch_count == 1)
    {
        binding_optional_scratch.reset(new AlignedScratch(64));
        items[0].scratch = binding_optional_scratch->data();
        items[0].scratch_bytes = 64;
    }
    leo2_encode_batch_binding* binding = NULL;
    require_result(leo2_encode_batch_binding_create(
        fixture.codec, &items[0], items.size(), &binding),
        LEO2_SUCCESS, "K=1 encode binding create");
    require(binding != NULL &&
            leo2_encode_batch_binding_item_count(binding) == items.size(),
        "K=1 encode binding item-count mismatch");
    for (size_t item_i = 0; item_i < batch_count; ++item_i)
        std::fill(output_storage[item_i].begin(),
            output_storage[item_i].end(), 0x5a);
    require_result(leo2_encode_batch_binding_execute(binding),
        LEO2_SUCCESS, "K=1 encode binding execute");
    check_outputs();

    uint8_t* const live_source =
        const_cast<uint8_t*>(static_cast<const uint8_t*>(fixture.original[0]));
    const uint8_t saved_source_byte = live_source[0];
    live_source[0] ^= 0xa7u;
    for (size_t item_i = 0; item_i < batch_count; ++item_i)
        std::fill(output_storage[item_i].begin(),
            output_storage[item_i].end(), 0x5a);
    require_result(leo2_encode_batch_binding_execute(binding),
        LEO2_SUCCESS, "K=1 encode binding changed-source execute");
    for (size_t item_i = 0; item_i < batch_count; ++item_i)
        require(std::memcmp(outputs[item_i][0], live_source,
                    fixture.bytes) == 0,
            "K=1 encode binding ignored changed source bytes");
    live_source[0] = saved_source_byte;
    for (size_t item_i = 0; item_i < batch_count; ++item_i)
        std::fill(output_storage[item_i].begin(),
            output_storage[item_i].end(), 0x5a);
    require_result(leo2_encode_batch_binding_execute(binding),
        LEO2_SUCCESS, "K=1 encode binding restored-source execute");
    check_outputs();
    leo2_encode_batch_binding_destroy(binding);
    if (batch_count == 1)
    {
        items[0].scratch = NULL;
        items[0].scratch_bytes = 0;
    }

    if (batch_count >= 9 && fixture.bytes == 4097)
    {
        AlignedScratch optional_item_scratch(64);
        items[0].scratch = optional_item_scratch.data();
        items[0].scratch_bytes = 64;
        for (size_t item_i = 0; item_i < batch_count; ++item_i)
            std::fill(output_storage[item_i].begin(),
                output_storage[item_i].end(), 0x5a);
        require_result(leo2_encode_batch_with_preflight_scratch(
            fixture.codec, &items[0], items.size(), preflight.data(),
            preflight_bytes), LEO2_SUCCESS,
            "K=1 scalable encode optional item scratch execute");
        check_outputs();
        items[0].scratch = NULL;
        items[0].scratch_bytes = 0;

        for (size_t item_i = 0; item_i < batch_count; ++item_i)
            std::fill(output_storage[item_i].begin(),
                output_storage[item_i].end(), 0x6b);
        const std::vector<Bytes> scalable_before_failure = output_storage;
        const void* const* const saved_original = items.back().original;
        items.back().original = NULL;
        require_result(leo2_encode_batch_with_preflight_scratch(
            fixture.codec, &items[0], items.size(), preflight.data(),
            preflight_bytes), LEO2_INVALID_ARGUMENT,
            "K=1 scalable encode late-item pointer rejection");
        require(output_storage == scalable_before_failure,
            "K=1 scalable encode ran work before preflight failure");
        items.back().original = saved_original;
    }

    if (batch_count < 2)
        return;

    for (size_t item_i = 0; item_i < batch_count; ++item_i)
        std::fill(output_storage[item_i].begin(),
            output_storage[item_i].end(), 0x6b);
    const std::vector<Bytes> before_failure = output_storage;

    const void* const* const saved_original = items.back().original;
    items.back().original = NULL;
    require_result(leo2_encode_batch(
        fixture.codec, &items[0], items.size()), LEO2_INVALID_ARGUMENT,
        "K=1 encode batch late-item pointer rejection");
    require(output_storage == before_failure,
        "K=1 encode batch ran an item before late pointer rejection");
    items.back().original = saved_original;

    void* const saved_output = outputs[0][0];
    outputs[0][0] = &outputs[1][0];
    require_result(leo2_encode_batch(
        fixture.codec, &items[0], items.size()), LEO2_OVERLAP,
        "K=1 encode batch cross-item output/metadata overlap");
    require(output_storage == before_failure,
        "K=1 encode batch metadata rejection modified an output");
    outputs[0][0] = saved_output;

    outputs[0][0] = outputs[1][0];
    require_result(leo2_encode_batch(
        fixture.codec, &items[0], items.size()), LEO2_OVERLAP,
        "K=1 encode batch cross-item output overlap");
    require(output_storage == before_failure,
        "K=1 encode batch alias rejection modified an output");
}

void test_public_r1_multi_item_batch(leo2_backend backend)
{
    static const uint32_t counts[] = { 1, 3, 5, 6, 7, 9, 129 };
    static const size_t sizes[] = { 4096, 4097 };
    static const size_t batches[] = { 2, 8, 64 };
    static const uint32_t thread_counts[] = { 1, 4 };
    for (size_t thread_i = 0;
         thread_i < sizeof(thread_counts) / sizeof(thread_counts[0]);
         ++thread_i)
    {
        for (size_t count_i = 0;
             count_i < sizeof(counts) / sizeof(counts[0]); ++count_i)
        {
            for (size_t size_i = 0;
                 size_i < sizeof(sizes) / sizeof(sizes[0]); ++size_i)
            {
                const uint32_t k = counts[count_i];
                R1Fixture fixture(backend, k, sizes[size_i], false,
                    LEO2_FIELD_GF8, LEO2_SHARD_LAYOUT_NATIVE_V1, k / 2U,
                    thread_counts[thread_i]);
                for (size_t batch_i = 0;
                     batch_i < sizeof(batches) / sizeof(batches[0]);
                     ++batch_i)
                {
                    const bool adversarial =
                        thread_counts[thread_i] == 4 &&
                        (k == 1 || k == 3 || k == 9) &&
                        sizes[size_i] == 4096 && batches[batch_i] == 8;
                    execute_public_r1_multi_item_batch(
                        fixture, batches[batch_i], adversarial);
                    if (k == 1)
                        execute_public_k1_encode_batch(
                            fixture, batches[batch_i]);
                }
            }
        }
    }
}

void test_public_r1_overlap_rejection(leo2_backend backend)
{
    R1Fixture fixture(backend, 9, 65, false, LEO2_FIELD_GF8,
        LEO2_SHARD_LAYOUT_NATIVE_V1);

    size_t encode_scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(
        fixture.codec, fixture.bytes, &encode_scratch_bytes), LEO2_SUCCESS,
        "R=1 overlap encode scratch query");
    AlignedScratch encode_scratch(encode_scratch_bytes);
    void* overlapping_parity[1] = {
        const_cast<void*>(fixture.original[3])
    };
    require_result(leo2_encode(fixture.codec, fixture.bytes,
        &fixture.original[0], overlapping_parity, encode_scratch.data(),
        encode_scratch_bytes), LEO2_OVERLAP,
        "R=1 encode input/output overlap");

    size_t decode_scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(
        fixture.plan, fixture.bytes, &decode_scratch_bytes), LEO2_SUCCESS,
        "R=1 overlap decode scratch query");
    AlignedScratch decode_scratch(decode_scratch_bytes);
    std::vector<const void*> received = fixture.original;
    received[0] = NULL;
    std::vector<void*> restored(fixture.k, NULL);

    restored[0] = const_cast<void*>(received[3]);
    require_result(leo2_decode_plan_execute(fixture.plan, fixture.bytes,
        &received[0], fixture.recovery, &restored[0], decode_scratch.data(),
        decode_scratch_bytes), LEO2_OVERLAP,
        "R=1 decode input/output overlap");

    restored[0] = const_cast<void*>(fixture.recovery[0]);
    require_result(leo2_decode_plan_execute(fixture.plan, fixture.bytes,
        &received[0], fixture.recovery, &restored[0], decode_scratch.data(),
        decode_scratch_bytes), LEO2_OVERLAP,
        "R=1 decode recovery/output overlap");

    std::vector<uint8_t> original_present(fixture.k, 1);
    original_present[0] = 0;
    const uint8_t recovery_present[1] = { 1 };
    restored[0] = const_cast<void*>(received[3]);
    require_result(leo2_decode(fixture.codec, fixture.bytes,
        &original_present[0], recovery_present, &received[0],
        fixture.recovery, &restored[0], decode_scratch.data(),
        decode_scratch_bytes), LEO2_OVERLAP,
        "R=1 one-shot input/output overlap");

    restored[0] = const_cast<void*>(fixture.recovery[0]);
    require_result(leo2_decode(fixture.codec, fixture.bytes,
        &original_present[0], recovery_present, &received[0],
        fixture.recovery, &restored[0], decode_scratch.data(),
        decode_scratch_bytes), LEO2_OVERLAP,
        "R=1 one-shot recovery/output overlap");

    Bytes ordinary_output(fixture.bytes, 0x7d);
    restored[0] = &ordinary_output[0];
    uint8_t* const original_presence_in_scratch =
        static_cast<uint8_t*>(decode_scratch.data());
    std::memcpy(original_presence_in_scratch, &original_present[0],
        original_present.size());
    require_result(leo2_decode(fixture.codec, fixture.bytes,
        original_presence_in_scratch, recovery_present, &received[0],
        fixture.recovery, &restored[0], decode_scratch.data(),
        decode_scratch_bytes), LEO2_OVERLAP,
        "R=1 one-shot scratch/original-presence overlap");

    Bytes original_presence_output(fixture.bytes, 1);
    original_presence_output[0] = 0;
    restored[0] = &original_presence_output[0];
    require_result(leo2_decode(fixture.codec, fixture.bytes,
        &original_presence_output[0], recovery_present, &received[0],
        fixture.recovery, &restored[0], decode_scratch.data(),
        decode_scratch_bytes), LEO2_OVERLAP,
        "R=1 one-shot output/original-presence overlap");

    Bytes recovery_presence_output(fixture.bytes, 1);
    restored[0] = &recovery_presence_output[0];
    require_result(leo2_decode(fixture.codec, fixture.bytes,
        &original_present[0], &recovery_presence_output[0], &received[0],
        fixture.recovery, &restored[0], decode_scratch.data(),
        decode_scratch_bytes), LEO2_OVERLAP,
        "R=1 one-shot output/recovery-presence overlap");

    /* Repeat the alias-sensitive cases through the zero-scratch K=1
       specialization rather than relying on the general K=9 validator. */
    R1Fixture single(backend, 1, 65, false, LEO2_FIELD_GF8,
        LEO2_SHARD_LAYOUT_NATIVE_V1);
    size_t single_encode_scratch = 1;
    require_result(leo2_encode_scratch_size(single.codec, single.bytes,
        &single_encode_scratch), LEO2_SUCCESS,
        "K=1 R=1 overlap encode scratch query");
    require(single_encode_scratch == 0,
        "K=1 R=1 overlap encode unexpectedly requires scratch");
    AlignedScratch optional_encode_scratch(65);
    Bytes optional_parity_storage(single.bytes + 9U, 0);
    void* optional_parity[1] = { &optional_parity_storage[3] };
    require_result(leo2_encode(single.codec, single.bytes,
        &single.original[0], optional_parity, optional_encode_scratch.data(),
        64), LEO2_SUCCESS, "K=1 R=1 optional encode scratch");
    require_result(leo2_encode(single.codec, single.bytes,
        &single.original[0], optional_parity,
        static_cast<uint8_t*>(optional_encode_scratch.data()) + 1, 64),
        LEO2_BAD_ALIGNMENT, "K=1 R=1 misaligned optional encode scratch");
    void* scratch_parity[1] = { optional_encode_scratch.data() };
    require_result(leo2_encode(single.codec, single.bytes,
        &single.original[0], scratch_parity, optional_encode_scratch.data(),
        64), LEO2_OVERLAP, "K=1 R=1 encode output/scratch overlap");
    void* single_parity[1] = {
        const_cast<void*>(single.original[0])
    };
    require_result(leo2_encode(single.codec, single.bytes,
        &single.original[0], single_parity, NULL, 0), LEO2_OVERLAP,
        "K=1 R=1 encode input/output overlap");

    size_t single_decode_scratch = 1;
    require_result(leo2_decode_plan_scratch_size(single.plan, single.bytes,
        &single_decode_scratch), LEO2_SUCCESS,
        "K=1 R=1 overlap decode scratch query");
    require(single_decode_scratch == 0,
        "K=1 R=1 overlap decode unexpectedly requires scratch");
    const void* single_received[1] = { NULL };
    void* single_restored[1] = {
        const_cast<void*>(single.recovery[0])
    };
    require_result(leo2_decode_plan_execute(single.plan, single.bytes,
        single_received, single.recovery, single_restored, NULL, 0),
        LEO2_OVERLAP, "K=1 R=1 decode recovery/output overlap");

    AlignedScratch optional_decode_scratch(65);
    Bytes optional_restored_storage(single.bytes + 9U, 0);
    single_restored[0] = &optional_restored_storage[3];
    require_result(leo2_decode_plan_execute(single.plan, single.bytes,
        single_received, single.recovery, single_restored,
        optional_decode_scratch.data(), 64), LEO2_SUCCESS,
        "K=1 R=1 optional decode scratch");
    require_result(leo2_decode_plan_execute(single.plan, single.bytes,
        single_received, single.recovery, single_restored,
        static_cast<uint8_t*>(optional_decode_scratch.data()) + 1, 64),
        LEO2_BAD_ALIGNMENT, "K=1 R=1 misaligned optional decode scratch");
    single_restored[0] = optional_decode_scratch.data();
    require_result(leo2_decode_plan_execute(single.plan, single.bytes,
        single_received, single.recovery, single_restored,
        optional_decode_scratch.data(), 64), LEO2_OVERLAP,
        "K=1 R=1 decode output/scratch overlap");

    uint8_t single_original_present[65] = {};
    const uint8_t single_recovery_present[1] = { 1 };
    single_restored[0] = single_original_present;
    require_result(leo2_decode(single.codec, single.bytes,
        single_original_present, single_recovery_present, single_received,
        single.recovery, single_restored, NULL, 0), LEO2_OVERLAP,
        "K=1 R=1 one-shot output/original-presence overlap");

    uint8_t single_recovery_presence[65] = {};
    single_recovery_presence[0] = 1;
    const uint8_t missing_original[1] = { 0 };
    single_restored[0] = single_recovery_presence;
    require_result(leo2_decode(single.codec, single.bytes, missing_original,
        single_recovery_presence, single_received, single.recovery,
        single_restored, NULL, 0), LEO2_OVERLAP,
        "K=1 R=1 one-shot output/recovery-presence overlap");
}

void test_k2_avx2_terminal_contract(leo2_backend backend)
{
    R1Fixture fixture(backend, 2, 2048, false, LEO2_FIELD_GF8,
        LEO2_SHARD_LAYOUT_NATIVE_V1);
    const leo2_backend effective = leo2_context_backend(fixture.context);
    if (effective != LEO2_BACKEND_AVX2 && effective != LEO2_BACKEND_GFNI)
        return;

    size_t one_shot_scratch_bytes = 1;
    require_result(leo2_decode_scratch_size(
        fixture.codec, fixture.bytes, &one_shot_scratch_bytes),
        LEO2_SUCCESS, "K=2 terminal one-shot scratch query");
    require(one_shot_scratch_bytes == 0,
        "K=2 terminal one-shot unexpectedly requires data scratch");
    size_t reusable_scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(
        fixture.plan, fixture.bytes, &reusable_scratch_bytes),
        LEO2_SUCCESS, "K=2 terminal reusable scratch query");
    require(reusable_scratch_bytes != 0,
        "K=2 terminal unexpectedly changed reusable-plan scratch contract");

    const uint8_t recovery_present[1] = { 1 };
    for (uint32_t missing = 0; missing < 2; ++missing)
    {
        std::vector<const void*> received = fixture.original;
        received[missing] = NULL;
        uint8_t original_present[2] = { 1, 1 };
        original_present[missing] = 0;
        Bytes restored_storage(fixture.bytes + 9U, 0xa6);
        void* restored[2] = { NULL, NULL };
        restored[missing] = &restored_storage[3];
        require_result(leo2_decode(fixture.codec, fixture.bytes,
            original_present, recovery_present, &received[0],
            fixture.recovery, restored, NULL, 0), LEO2_SUCCESS,
            "K=2 zero-scratch one-shot decode");
        require(std::memcmp(restored[missing], fixture.original[missing],
                    fixture.bytes) == 0,
            "K=2 zero-scratch one-shot restored the wrong original");
        for (size_t i = 0; i < 3; ++i)
            require(restored_storage[i] == 0xa6,
                "K=2 zero-scratch one-shot changed a prefix guard");
        for (size_t i = 3 + fixture.bytes;
             i < restored_storage.size(); ++i)
            require(restored_storage[i] == 0xa6,
                "K=2 zero-scratch one-shot changed a suffix guard");
    }

    /* A no-loss decode remains a true no-op and must not inspect shard or
       output arrays.  Presence metadata itself is still validated. */
    const uint8_t all_originals_present[2] = { 1, 1 };
    const uint8_t parity_absent[1] = { 0 };
    require_result(leo2_decode(fixture.codec, fixture.bytes,
        all_originals_present, parity_absent, NULL, NULL, NULL, NULL, 0),
        LEO2_SUCCESS, "K=2 terminal no-loss no-op");
    uint8_t invalid_original_presence[2] = { 2, 1 };
    require_result(leo2_decode(fixture.codec, fixture.bytes,
        invalid_original_presence, parity_absent, NULL, NULL, NULL, NULL, 0),
        LEO2_INVALID_ARGUMENT, "K=2 invalid original presence");
    const uint8_t invalid_recovery_presence[1] = { 2 };
    require_result(leo2_decode(fixture.codec, fixture.bytes,
        all_originals_present, invalid_recovery_presence,
        NULL, NULL, NULL, NULL, 0), LEO2_INVALID_ARGUMENT,
        "K=2 invalid recovery presence");

    std::vector<const void*> received = fixture.original;
    received[0] = NULL;
    uint8_t original_present[2] = { 0, 1 };
    Bytes ordinary_output(fixture.bytes, 0x4d);
    void* restored[2] = { &ordinary_output[0], NULL };
    require_result(leo2_decode_plan_execute(fixture.plan, fixture.bytes,
        &received[0], fixture.recovery, restored, NULL, 0),
        LEO2_SCRATCH_TOO_SMALL,
        "K=2 terminal reusable plan accepted absent scratch");
    AlignedScratch optional_scratch(65);
    require_result(leo2_decode(fixture.codec, fixture.bytes,
        original_present, recovery_present, &received[0], fixture.recovery,
        restored, optional_scratch.data(), 64), LEO2_SUCCESS,
        "K=2 optional aligned scratch");
    require_result(leo2_decode(fixture.codec, fixture.bytes,
        original_present, recovery_present, &received[0], fixture.recovery,
        restored, static_cast<uint8_t*>(optional_scratch.data()) + 1, 64),
        LEO2_BAD_ALIGNMENT, "K=2 optional misaligned scratch");
    restored[0] = optional_scratch.data();
    require_result(leo2_decode(fixture.codec, fixture.bytes,
        original_present, recovery_present, &received[0], fixture.recovery,
        restored, optional_scratch.data(), 64), LEO2_OVERLAP,
        "K=2 output/optional-scratch overlap");
    restored[0] = &ordinary_output[0];
    uint8_t* presence_in_scratch =
        static_cast<uint8_t*>(optional_scratch.data());
    presence_in_scratch[0] = 0;
    presence_in_scratch[1] = 1;
    require_result(leo2_decode(fixture.codec, fixture.bytes,
        presence_in_scratch, recovery_present, &received[0], fixture.recovery,
        restored, optional_scratch.data(), 64), LEO2_OVERLAP,
        "K=2 presence/optional-scratch overlap");

    size_t encode_scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(
        fixture.codec, fixture.bytes, &encode_scratch_bytes),
        LEO2_SUCCESS, "K=2 terminal encode scratch query");
    AlignedScratch encode_scratch(encode_scratch_bytes);
    void* overlapping_parity[1] = {
        const_cast<void*>(fixture.original[1])
    };
    require_result(leo2_encode(fixture.codec, fixture.bytes,
        &fixture.original[0], overlapping_parity, encode_scratch.data(),
        encode_scratch_bytes), LEO2_OVERLAP,
        "K=2 terminal encode input/output overlap");

    leo2_encode_batch_item encode_item = {};
    void* batch_parity[1] = { &encode_item };
    encode_item.shard_bytes = fixture.bytes;
    encode_item.original = &fixture.original[0];
    encode_item.recovery = batch_parity;
    encode_item.scratch = encode_scratch.data();
    encode_item.scratch_bytes = encode_scratch_bytes;
    require_result(leo2_encode_batch(
        fixture.codec, &encode_item, 1), LEO2_OVERLAP,
        "K=2 terminal batch output/descriptor overlap");

    AlignedScratch reusable_scratch(reusable_scratch_bytes);
    leo2_decode_batch_item decode_item = {};
    void* batch_restored[2] = { &decode_item, NULL };
    decode_item.shard_bytes = fixture.bytes;
    decode_item.original = &received[0];
    decode_item.recovery = fixture.recovery;
    decode_item.restored_original = batch_restored;
    decode_item.scratch = reusable_scratch.data();
    decode_item.scratch_bytes = reusable_scratch_bytes;
    require_result(leo2_decode_plan_execute_batch(
        fixture.plan, &decode_item, 1), LEO2_OVERLAP,
        "K=2 terminal batch output/descriptor overlap");

    R1Fixture below_threshold(backend, 2, 2047, false, LEO2_FIELD_GF8,
        LEO2_SHARD_LAYOUT_NATIVE_V1);
    size_t below_scratch_bytes = 0;
    require_result(leo2_decode_scratch_size(below_threshold.codec,
        below_threshold.bytes, &below_scratch_bytes), LEO2_SUCCESS,
        "K=2 below-threshold one-shot scratch query");
    require(below_scratch_bytes != 0,
        "K=2 terminal escaped below its measured byte threshold");

    R1Fixture gf16_fixture(backend, 2, 2048, false, LEO2_FIELD_GF16,
        LEO2_SHARD_LAYOUT_NATIVE_V1);
    size_t gf16_scratch_bytes = 0;
    require_result(leo2_decode_scratch_size(gf16_fixture.codec,
        gf16_fixture.bytes, &gf16_scratch_bytes), LEO2_SUCCESS,
        "K=2 GF16 one-shot scratch query");
    require(gf16_scratch_bytes != 0,
        "GF8 K=2 terminal changed the GF16 scratch contract");
}

void test_public_r1(leo2_backend backend)
{
    /* K=1 is a copy-only code.  Its specialized validators do not stage any
       ranges, so both operations must accept the queried zero-byte scratch. */
    {
        R1Fixture fixture(backend, 1, 65, false, LEO2_FIELD_GF8,
            LEO2_SHARD_LAYOUT_NATIVE_V1);
        size_t encode_scratch = 1;
        size_t decode_scratch = 1;
        size_t one_shot_decode_scratch = 1;
        require_result(leo2_encode_scratch_size(
            fixture.codec, fixture.bytes, &encode_scratch), LEO2_SUCCESS,
            "K=1 R=1 encode scratch query");
        require_result(leo2_decode_plan_scratch_size(
            fixture.plan, fixture.bytes, &decode_scratch), LEO2_SUCCESS,
            "K=1 R=1 decode scratch query");
        require_result(leo2_decode_scratch_size(
            fixture.codec, fixture.bytes, &one_shot_decode_scratch),
            LEO2_SUCCESS, "K=1 R=1 one-shot decode scratch query");
        require(encode_scratch == 0 && decode_scratch == 0 &&
                one_shot_decode_scratch == 0,
            "K=1 R=1 copy-only path unexpectedly requires scratch");
        execute_and_check_decode(fixture);
    }

    // The production R=1 path needs only address-range validation scratch.
    // Its scratch footprint must therefore not grow with shard bytes as a
    // transform workspace would.
    {
        R1Fixture fixture(backend, 31, 64, false, LEO2_FIELD_GF8,
            LEO2_SHARD_LAYOUT_NATIVE_V1);
        size_t encode_small = 0;
        size_t encode_large = 0;
        size_t decode_small = 0;
        size_t decode_large = 0;
        require_result(leo2_encode_scratch_size(
            fixture.codec, 64, &encode_small), LEO2_SUCCESS,
            "R=1 small encode scratch query");
        require_result(leo2_encode_scratch_size(
            fixture.codec, 1024 * 1024, &encode_large), LEO2_SUCCESS,
            "R=1 large encode scratch query");
        require_result(leo2_decode_plan_scratch_size(
            fixture.plan, 64, &decode_small), LEO2_SUCCESS,
            "R=1 small decode scratch query");
        require_result(leo2_decode_plan_scratch_size(
            fixture.plan, 1024 * 1024, &decode_large), LEO2_SUCCESS,
            "R=1 large decode scratch query");
        require(encode_small == encode_large,
            "R=1 encode scratch still scales with shard bytes");
        require(decode_small == decode_large,
            "R=1 decode scratch still scales with shard bytes");
    }

    /* K=8 has seven live originals during one-loss repair and exercises the
       fused initial-plus-seven AVX2 reduction; K=9 exercises the full group
       of eight live sources. */
    static const uint32_t counts[] = { 1, 2, 3, 4, 5, 6, 8, 9, 10, 31 };
    static const size_t sizes[] = {
        1, 2, 3, 17, 31, 32, 33, 64, 65, 1024, 1025,
        4095, 4096, 4097, 65537, 1048579
    };
    for (size_t count_i = 0;
         count_i < sizeof(counts) / sizeof(counts[0]); ++count_i)
        for (size_t size_i = 0;
             size_i < sizeof(sizes) / sizeof(sizes[0]); ++size_i)
        {
            R1Fixture fixture(backend, counts[count_i], sizes[size_i],
                counts[count_i] == 3, LEO2_FIELD_GF8,
                LEO2_SHARD_LAYOUT_NATIVE_V1);
            execute_and_check_decode(fixture);
        }

    /* One-item reusable bindings are the public execution terminal for the
       K=1,R=1 copy code.  Cover each measured selection boundary, including
       both sides of the 64-KiB cutoff, with unaligned guarded buffers and an
       optional caller scratch span. */
    static const size_t k1_binding_sizes[] = {
        1, 64, 1024, 4096, 4097, 65536, 65537
    };
    for (size_t size_i = 0;
         size_i < sizeof(k1_binding_sizes) / sizeof(k1_binding_sizes[0]);
         ++size_i)
    {
        R1Fixture fixture(backend, 1, k1_binding_sizes[size_i], false,
            LEO2_FIELD_GF8, LEO2_SHARD_LAYOUT_NATIVE_V1);
        execute_public_r1_multi_item_batch(fixture, 1, false);
        execute_public_k1_encode_batch(fixture, 1);
    }

    /* K=2 enters the dense two-input reducer and K=7 the exact-arity
       fused-final reducer at 2 KiB on the GF8 AVX2 tier.  Recover every
       original on both sides of that boundary, including vector and page
       tails, so a null-source or unaligned-buffer error cannot hide behind
       the missing-zero route checks. */
    static const uint32_t tiny_threshold_counts[] = { 2, 7 };
    static const size_t tiny_threshold_sizes[] = {
        2047, 2048, 2049, 2079, 3072, 4095
    };
    for (size_t count_i = 0;
         count_i < sizeof(tiny_threshold_counts) /
             sizeof(tiny_threshold_counts[0]); ++count_i)
    {
        const uint32_t k = tiny_threshold_counts[count_i];
        for (uint32_t missing = 0; missing < k; ++missing)
            for (size_t size_i = 0;
                 size_i < sizeof(tiny_threshold_sizes) /
                     sizeof(tiny_threshold_sizes[0]); ++size_i)
            {
                R1Fixture fixture(backend, k,
                    tiny_threshold_sizes[size_i], false, LEO2_FIELD_GF8,
                    LEO2_SHARD_LAYOUT_NATIVE_V1, missing);
                execute_and_check_decode(fixture);
            }
    }

    /* K=3,5,6 cross to the exact-arity AVX2 reduction at 4 KiB.  Recover
       every original around that boundary, through an arbitrary vector tail,
       and beyond one MiB so null compaction and the unbounded large-shard
       policy cannot be validated only at missing index zero. */
    static const uint32_t small_fused_counts[] = { 3, 5, 6 };
    static const size_t small_fused_sizes[] = {
        4095, 4096, 4097, 4127, 65537, 1048593
    };
    for (size_t count_i = 0;
         count_i < sizeof(small_fused_counts) /
             sizeof(small_fused_counts[0]); ++count_i)
    {
        const uint32_t k = small_fused_counts[count_i];
        for (uint32_t missing = 0; missing < k; ++missing)
            for (size_t size_i = 0;
                 size_i < sizeof(small_fused_sizes) /
                     sizeof(small_fused_sizes[0]); ++size_i)
            {
                R1Fixture fixture(backend, k, small_fused_sizes[size_i],
                    false, LEO2_FIELD_GF8, LEO2_SHARD_LAYOUT_NATIVE_V1,
                    missing);
                execute_and_check_decode(fixture);
            }
    }

    /* These K values leave final live-source remainders three through six
       after one or two complete eight-source groups.  R1Fixture offsets each
       shard independently, so the loop covers threshold, tail, unaligned,
       encode, and one-loss decode behavior against the scalar XOR oracle. */
    static const uint32_t final_remainder_counts[] = {
        7, 12, 13, 14, 15, 20, 21, 22, 23
    };
    static const size_t final_remainder_sizes[] = {
        4096, 4097, 4113, 65553, 1048593
    };
    for (size_t count_i = 0;
         count_i < sizeof(final_remainder_counts) /
             sizeof(final_remainder_counts[0]); ++count_i)
        for (size_t size_i = 0;
             size_i < sizeof(final_remainder_sizes) /
                 sizeof(final_remainder_sizes[0]); ++size_i)
        {
            R1Fixture fixture(backend, final_remainder_counts[count_i],
                final_remainder_sizes[size_i], false, LEO2_FIELD_GF8,
                LEO2_SHARD_LAYOUT_NATIVE_V1,
                final_remainder_counts[count_i] / 2U);
            execute_and_check_decode(fixture);
        }

    /* K=7 is the smallest live remainder-six case and has no complete
       eight-source group, so the final group is the entire coarse reduction.
       Exercise the large shard sizes where one pass can matter most. */
    static const size_t k7_large_sizes[] = {
        1048576, 1048593,
        4194304, 4194321,
        8388608, 8388625,
        16777216, 16777233
    };
    for (size_t size_i = 0;
         size_i < sizeof(k7_large_sizes) / sizeof(k7_large_sizes[0]); ++size_i)
    {
        R1Fixture fixture(backend, 7, k7_large_sizes[size_i], false,
            LEO2_FIELD_GF8, LEO2_SHARD_LAYOUT_NATIVE_V1, 3);
        execute_and_check_decode(fixture);
    }

    /* Null positions before, at, and after a complete eight-source group must
       produce the same K-1 live count and final-remainder selection. */
    static const uint32_t null_boundary_counts[] = { 7, 12, 22, 23 };
    static const size_t null_boundary_sizes[] = { 4097, 65553 };
    for (size_t count_i = 0;
         count_i < sizeof(null_boundary_counts) /
             sizeof(null_boundary_counts[0]); ++count_i)
    {
        const uint32_t k = null_boundary_counts[count_i];
        static const uint32_t boundary_positions[] = { 0, 7, 8, 15, 16 };
        std::vector<uint32_t> positions;
        for (size_t boundary_i = 0;
             boundary_i < sizeof(boundary_positions) /
                 sizeof(boundary_positions[0]); ++boundary_i)
            if (boundary_positions[boundary_i] < k)
                positions.push_back(boundary_positions[boundary_i]);
        positions.push_back(k - 1U);
        std::sort(positions.begin(), positions.end());
        positions.erase(std::unique(positions.begin(), positions.end()),
            positions.end());
        for (size_t position_i = 0;
             position_i < positions.size();
             ++position_i)
            for (size_t size_i = 0;
                 size_i < sizeof(null_boundary_sizes) /
                     sizeof(null_boundary_sizes[0]); ++size_i)
            {
                R1Fixture fixture(backend, k, null_boundary_sizes[size_i],
                    false, LEO2_FIELD_GF8, LEO2_SHARD_LAYOUT_NATIVE_V1,
                    positions[position_i]);
                execute_and_check_decode(fixture);
            }
    }

    // GF16 uses the same bytewise XOR identity, but exercises a distinct codec
    // selection and shard-layout contract.  Even physical byte counts include
    // vector tails; padded-odd physical sizes represent an odd payload with a
    // required zero high byte in every field-symbol stream.
    static const size_t gf16_native_sizes[] = {
        2, 6, 34, 66, 1026, 65538
    };
    for (size_t size_i = 0;
         size_i < sizeof(gf16_native_sizes) / sizeof(gf16_native_sizes[0]);
         ++size_i)
    {
        R1Fixture fixture(backend, size_i % 2 == 0 ? 3U : 10U,
            gf16_native_sizes[size_i], size_i == 0, LEO2_FIELD_GF16,
            LEO2_SHARD_LAYOUT_NATIVE_V1);
        execute_and_check_decode(fixture);
    }
    R1Fixture gf16_single(backend, 1, 66, false, LEO2_FIELD_GF16,
        LEO2_SHARD_LAYOUT_NATIVE_V1);
    execute_and_check_decode(gf16_single);
    R1Fixture gf16_single_odd_payload(backend, 1, 34, false,
        LEO2_FIELD_GF16, LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1);
    execute_and_check_decode(gf16_single_odd_payload);
    static const uint32_t gf16_k1_threads[] = { 1, 4 };
    for (size_t thread_i = 0;
         thread_i < sizeof(gf16_k1_threads) / sizeof(gf16_k1_threads[0]);
         ++thread_i)
    {
        R1Fixture native_batch(backend, 1, 66, false, LEO2_FIELD_GF16,
            LEO2_SHARD_LAYOUT_NATIVE_V1, 0, gf16_k1_threads[thread_i]);
        execute_public_r1_multi_item_batch(native_batch, 9, false);
        execute_public_k1_encode_batch(native_batch, 9);
        R1Fixture padded_batch(backend, 1, 34, false, LEO2_FIELD_GF16,
            LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1, 0,
            gf16_k1_threads[thread_i]);
        execute_public_r1_multi_item_batch(padded_batch, 9, false);
        execute_public_k1_encode_batch(padded_batch, 9);
    }
    R1Fixture gf16_odd_payload(backend, 9, 34, true, LEO2_FIELD_GF16,
        LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1);
    execute_and_check_decode(gf16_odd_payload);
    R1Fixture gf16_boundary(backend, 256, 66, false, LEO2_FIELD_GF16,
        LEO2_SHARD_LAYOUT_NATIVE_V1);
    execute_and_check_decode(gf16_boundary);

    // Exact legacy comparisons at the GF8 parent boundary and the first GF16
    // R=1 parent supplement the independent XOR oracle above.
    R1Fixture gf8_legacy_boundary(backend, 129, 64, false, LEO2_FIELD_GF8,
        LEO2_SHARD_LAYOUT_NATIVE_V1);
    execute_and_check_decode(gf8_legacy_boundary);
    R1Fixture gf16_legacy_boundary(backend, 256, 64, false, LEO2_FIELD_GF16,
        LEO2_SHARD_LAYOUT_NATIVE_V1);
    execute_and_check_decode(gf16_legacy_boundary);

    R1Fixture interior_missing(backend, 9, 1025, false, LEO2_FIELD_GF8,
        LEO2_SHARD_LAYOUT_NATIVE_V1, 4);
    execute_and_check_decode(interior_missing);
    R1Fixture final_missing(backend, 9, 1025, false, LEO2_FIELD_GF8,
        LEO2_SHARD_LAYOUT_NATIVE_V1, 8);
    execute_and_check_decode(final_missing);

    test_public_r1_overlap_rejection(backend);
    test_k2_avx2_terminal_contract(backend);
    test_public_r1_multi_item_batch(backend);

    // Shared immutable mature and newly compact/fused plans must be race-free.
    static const uint32_t concurrent_counts[] = { 2, 5, 9 };
    for (size_t count_i = 0;
         count_i < sizeof(concurrent_counts) /
             sizeof(concurrent_counts[0]); ++count_i)
    {
        R1Fixture fixture(backend, concurrent_counts[count_i], 4097, true,
            LEO2_FIELD_GF8, LEO2_SHARD_LAYOUT_NATIVE_V1);
        std::atomic<bool> failed(false);
        std::vector<std::thread> threads;
        for (unsigned thread = 0; thread < 16; ++thread)
        {
            threads.push_back(std::thread([&fixture, &failed]() {
                try
                {
                    for (unsigned round = 0; round < 32; ++round)
                        execute_and_check_decode(fixture);
                }
                catch (...)
                {
                    failed.store(true, std::memory_order_relaxed);
                }
            }));
        }
        for (size_t i = 0; i < threads.size(); ++i)
            threads[i].join();
        require(!failed.load(std::memory_order_relaxed),
            "concurrent R=1 plan execution failed");
    }
}

bool public_backend_available(leo2_backend requested)
{
    leo2_context_options options;
    std::memset(&options, 0, sizeof(options));
    options.struct_size = sizeof(options);
    options.backend = requested;
    leo2_context* context = NULL;
    const leo2_result result = leo2_context_create(&options, &context);
    if (result == LEO2_UNSUPPORTED)
    {
        require(context == NULL,
            "unsupported explicit backend returned a context");
        return false;
    }
    require_result(result, LEO2_SUCCESS,
        "explicit backend availability probe");
    require(context != NULL, "available explicit backend omitted a context");
    require(leo2_context_backend(context) == requested,
        "explicit backend probe selected a different effective backend");
    leo2_context_destroy(context);
    return true;
}

} // namespace

int main()
{
    try
    {
        require(leo_init() == Leopard_Success, "Leopard initialization");
        const leo2_backend backends[] = {
            LEO2_BACKEND_SCALAR,
            LEO2_BACKEND_SSSE3,
            LEO2_BACKEND_AVX2,
            LEO2_BACKEND_AVX512,
            LEO2_BACKEND_GFNI
        };
        unsigned ops_tested = 0;
        for (size_t i = 0; i < sizeof(backends) / sizeof(backends[0]); ++i)
        {
            leopard::backend::QualificationStatus status =
                leopard::backend::QualificationAvailable;
            const leopard::backend::Ops* ops =
                leopard::backend::GetQualifiedOps(backends[i], &status);
            if (!ops)
            {
                require(backends[i] != LEO2_BACKEND_SCALAR,
                    "scalar fused XOR backend is unavailable");
                require(status == leopard::backend::QualificationUnavailable,
                    "fused XOR backend failed qualification");
                continue;
            }
            require(ops->kind == backends[i], "qualified the wrong XOR backend");
            require(ops->xor_memory_2to1 != NULL,
                "qualified backend omitted fused XOR");
            test_primitive(*ops);
            test_primitive_concurrency(*ops);
            ++ops_tested;
        }
        require(ops_tested != 0,
            "no qualified fused XOR backend was available for testing");

        // Runtime NEON intentionally rejects lower x86/scalar context requests
        // until its complete transform table is extracted.  Probe the public
        // contract rather than assuming every directly qualified fixed-ops
        // table is independently selectable on every architecture.
        std::vector<leo2_backend> public_requests(
            backends, backends + sizeof(backends) / sizeof(backends[0]));
        const leo2_backend effective = leopard::backend::ExecutionBackend();
        if (std::find(public_requests.begin(), public_requests.end(),
                effective) == public_requests.end())
            public_requests.push_back(effective);
        unsigned public_tested = 0;
        bool effective_tested = false;
        for (size_t i = 0; i < public_requests.size(); ++i)
        {
            if (!public_backend_available(public_requests[i]))
                continue;
            test_public_r1(public_requests[i]);
            ++public_tested;
            effective_tested = effective_tested ||
                public_requests[i] == effective;
        }
        require(public_tested != 0,
            "no public R=1 backend was available for testing");
        require(effective_tested,
            "effective runtime backend lacked public R=1 coverage");
        std::printf("Leopard2 R=1 fused XOR passed: backends=qualified "
            "tails=0..257 max_bytes=16777233 fields=gf8,gf16 "
            "exact_arity_k=3,5,6 "
            "final_remainder_k=7,12..15,20..23 "
            "null_boundaries=0,7,8,15,16,last "
            "batch=2,8,64 batch_k=3,5,6,7,9,129 "
            "batch_bytes=4096,4097 "
            "batch_threads=1,4 concurrency=pass\n");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "Leopard2 R=1 fused XOR failed: %s\n",
            error.what());
        return 1;
    }
}
