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
        throw std::runtime_error(message);
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

void test_primitive(const leopard::backend::Ops& ops)
{
    const std::vector<uint64_t> counts = byte_counts();
    for (size_t i = 0; i < counts.size(); ++i)
        test_primitive_case(ops, counts[i],
            static_cast<uint32_t>(i * 17U + ops.kind * 101U));
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
    size_t bytes;
    std::vector<Bytes> storage;
    std::vector<const void*> original;
    Bytes recovery_storage;
    const void* recovery[1];

    R1Fixture(leo2_backend backend, uint32_t original_count,
        size_t shard_bytes, bool alias_inputs, leo2_field field,
        leo2_shard_layout layout)
        : context(NULL)
        , codec(NULL)
        , plan(NULL)
        , k(original_count)
        , bytes(shard_bytes)
        , storage(original_count)
        , original(original_count, NULL)
        , recovery_storage(shard_bytes + 17U)
    {
        leo2_context_options options;
        std::memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = backend;
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

        std::vector<uint8_t> original_present(k, 1);
        const uint8_t recovery_present[1] = { 1 };
        original_present[0] = 0;
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
    received[0] = NULL;
    Bytes restored_storage(fixture.bytes + 11U, 0xa5);
    std::vector<void*> restored(fixture.k, NULL);
    restored[0] = &restored_storage[5];
    require_result(leo2_decode_plan_execute(fixture.plan, fixture.bytes,
        &received[0], fixture.recovery, &restored[0], scratch.data(),
        scratch_bytes), LEO2_SUCCESS, "R=1 decode execute");
    require(std::memcmp(restored[0], fixture.original[0], fixture.bytes) == 0,
        "R=1 decode restored the wrong original");
    for (size_t i = 0; i < 5; ++i)
        require(restored_storage[i] == 0xa5, "R=1 decode changed a prefix guard");
    for (size_t i = 5 + fixture.bytes; i < restored_storage.size(); ++i)
        require(restored_storage[i] == 0xa5, "R=1 decode changed a suffix guard");
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
}

void test_public_r1(leo2_backend backend)
{
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

    static const uint32_t counts[] = { 3, 9, 10, 31 };
    static const size_t sizes[] = {
        1, 2, 3, 17, 31, 32, 33, 64, 65, 1025, 65537, 1048579
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
    R1Fixture gf16_odd_payload(backend, 9, 34, true, LEO2_FIELD_GF16,
        LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1);
    execute_and_check_decode(gf16_odd_payload);
    R1Fixture gf16_boundary(backend, 256, 66, false, LEO2_FIELD_GF16,
        LEO2_SHARD_LAYOUT_NATIVE_V1);
    execute_and_check_decode(gf16_boundary);

    test_public_r1_overlap_rejection(backend);

    // Shared immutable codec/plan execution must remain race-free.
    R1Fixture fixture(backend, 9, 4097, true, LEO2_FIELD_GF8,
        LEO2_SHARD_LAYOUT_NATIVE_V1);
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
            LEO2_BACKEND_AVX2
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
            "tails=0..257 max_bytes=1048579 fields=gf8,gf16 "
            "concurrency=pass\n");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "Leopard2 R=1 fused XOR failed: %s\n",
            error.what());
        return 1;
    }
}
