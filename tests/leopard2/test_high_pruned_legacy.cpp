/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in the Leopard-RS
    BSD license are met.  See the repository LICENSE file.
*/

/*
    Exact Leopard1 diagnostic for the production Algorithm 5 schedule path.
    Old and new encoders must emit identical parity; old Leopard's full decoder
    and both new pruned workspaces must restore identical original bytes.
*/

#include "leopard.h"
#include "leopard2.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

namespace {

class AlignedBytes
{
public:
    explicit AlignedBytes(size_t bytes)
        : data_(NULL), bytes_(bytes)
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

    ~AlignedBytes()
    {
#if defined(_MSC_VER)
        _aligned_free(data_);
#else
        free(data_);
#endif
    }

    uint8_t* data() { return static_cast<uint8_t*>(data_); }
    const uint8_t* data() const { return static_cast<const uint8_t*>(data_); }
    size_t size() const { return bytes_; }

private:
    AlignedBytes(const AlignedBytes&);
    AlignedBytes& operator=(const AlignedBytes&);
    void* data_;
    size_t bytes_;
};

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void require_result(leo2_result result, const char* operation)
{
    if (result != LEO2_SUCCESS)
        throw std::runtime_error(std::string(operation) + ": " +
            leo2_result_string(result));
}

uint32_t mix(uint32_t value)
{
    value ^= value >> 16;
    value *= 0x7feb352du;
    value ^= value >> 15;
    value *= 0x846ca68bu;
    return value ^ (value >> 16);
}

struct Case
{
    uint32_t k;
    uint32_t r;
    leo2_field field;
    size_t bytes;
};

void run_case(leo2_context* context, const Case& test, size_t case_index,
    uint64_t& parity_bytes, uint64_t& restored_bytes)
{
    const size_t stride = (test.bytes + 63u) & ~size_t(63u);
    AlignedBytes originals(static_cast<size_t>(test.k) * stride);
    std::vector<const void*> original_input(test.k);
    for (uint32_t shard = 0; shard < test.k; ++shard)
    {
        uint8_t* data = originals.data() + static_cast<size_t>(shard) * stride;
        for (size_t byte = 0; byte < test.bytes; ++byte)
            data[byte] = static_cast<uint8_t>(
                mix(static_cast<uint32_t>(case_index * 65537u +
                    shard * 257u + byte)));
        original_input[shard] = data;
    }

    const unsigned old_encode_count = leo_encode_work_count(test.k, test.r);
    require(old_encode_count != 0, "legacy encode work count");
    AlignedBytes old_encode_storage(
        static_cast<size_t>(old_encode_count) * stride);
    std::vector<void*> old_encode_work(old_encode_count);
    for (unsigned i = 0; i < old_encode_count; ++i)
        old_encode_work[i] = old_encode_storage.data() +
            static_cast<size_t>(i) * stride;
    require(leo_encode(stride, test.k, test.r, old_encode_count,
        &original_input[0], &old_encode_work[0]) == Leopard_Success,
        "legacy encode");

    leo2_codec_options codec_options = {};
    codec_options.struct_size = sizeof(codec_options);
    codec_options.flags = LEO2_CODEC_FORCE_SPECIALIZED_DECODE |
        ((case_index & 1u) != 0 ? LEO2_CODEC_FORCE_TILED_DECODE
                                : LEO2_CODEC_FORCE_MATERIALIZED_DECODE);
    leo2_codec* codec = NULL;
    require_result(leo2_codec_create(context, test.k, test.r,
        LEO2_PROFILE_LEGACY_HIGH_V1, test.field, &codec_options, &codec),
        "codec create");

    AlignedBytes new_recovery_storage(static_cast<size_t>(test.r) * stride);
    std::vector<void*> new_recovery_output(test.r);
    for (uint32_t i = 0; i < test.r; ++i)
        new_recovery_output[i] = new_recovery_storage.data() +
            static_cast<size_t>(i) * stride;
    size_t encode_scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(
        codec, test.bytes, &encode_scratch_bytes), "encode scratch query");
    AlignedBytes encode_scratch(encode_scratch_bytes);
    require_result(leo2_encode(codec, test.bytes, &original_input[0],
        &new_recovery_output[0], encode_scratch.data(), encode_scratch.size()),
        "new encode");
    for (uint32_t i = 0; i < test.r; ++i)
    {
        require(memcmp(new_recovery_output[i], old_encode_work[i],
            test.bytes) == 0, "new parity differs from Leopard1");
        parity_bytes += test.bytes;
    }

    const uint32_t loss_count = std::min<uint32_t>(5,
        std::min(test.k, test.r));
    std::vector<uint32_t> missing(loss_count);
    for (uint32_t i = 0; i < loss_count; ++i)
        missing[i] = loss_count == 1 ? 0 :
            static_cast<uint32_t>((static_cast<uint64_t>(i) * (test.k - 1)) /
                (loss_count - 1));

    std::vector<const void*> old_original = original_input;
    std::vector<uint8_t> original_present(test.k, 1);
    for (uint32_t i = 0; i < loss_count; ++i)
    {
        old_original[missing[i]] = NULL;
        original_present[missing[i]] = 0;
    }
    std::vector<const void*> old_recovery(test.r);
    for (uint32_t i = 0; i < test.r; ++i)
        old_recovery[i] = old_encode_work[i];
    const unsigned old_decode_count = leo_decode_work_count(test.k, test.r);
    require(old_decode_count != 0, "legacy decode work count");
    AlignedBytes old_decode_storage(
        static_cast<size_t>(old_decode_count) * stride);
    std::vector<void*> old_decode_work(old_decode_count);
    for (unsigned i = 0; i < old_decode_count; ++i)
        old_decode_work[i] = old_decode_storage.data() +
            static_cast<size_t>(i) * stride;
    require(leo_decode(stride, test.k, test.r, old_decode_count,
        &old_original[0], &old_recovery[0], &old_decode_work[0]) ==
        Leopard_Success, "legacy decode");

    std::vector<uint8_t> recovery_present(test.r, 1);
    leo2_decode_plan* plan = NULL;
    require_result(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), "decode plan create");
    size_t decode_scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(
        plan, test.bytes, &decode_scratch_bytes), "decode scratch query");
    AlignedBytes decode_scratch(decode_scratch_bytes);
    AlignedBytes restored_storage(static_cast<size_t>(test.k) * stride);
    std::vector<void*> restored(test.k, NULL);
    for (uint32_t i = 0; i < loss_count; ++i)
        restored[missing[i]] = restored_storage.data() +
            static_cast<size_t>(missing[i]) * stride;
    std::vector<const void*> new_recovery(test.r);
    for (uint32_t i = 0; i < test.r; ++i)
        new_recovery[i] = new_recovery_output[i];
    require_result(leo2_decode_plan_execute(plan, test.bytes,
        &old_original[0], &new_recovery[0], &restored[0],
        decode_scratch.data(), decode_scratch.size()), "new decode");
    for (uint32_t i = 0; i < loss_count; ++i)
    {
        const uint32_t shard = missing[i];
        require(memcmp(restored[shard], old_decode_work[shard],
            test.bytes) == 0, "new decode differs from Leopard1");
        require(memcmp(restored[shard], original_input[shard],
            test.bytes) == 0, "restored original mismatch");
        restored_bytes += test.bytes;
    }

    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
}

} // namespace

int main()
{
    try
    {
        require(leo_init() == Leopard_Success, "Leopard initialization");
        leo2_context_options options = {};
        options.struct_size = sizeof(options);
        options.thread_count = 1;
        leo2_context* context = NULL;
        require_result(leo2_context_create(&options, &context),
            "context create");
        const Case cases[] = {
            { 3, 2, LEO2_FIELD_GF8, 65 },
            { 7, 3, LEO2_FIELD_GF8, 64 },
            { 17, 5, LEO2_FIELD_GF8, 129 },
            { 33, 9, LEO2_FIELD_GF8, 64 },
            { 65, 17, LEO2_FIELD_GF8, 257 },
            { 127, 31, LEO2_FIELD_GF8, 64 },
            { 191, 32, LEO2_FIELD_GF8, 129 },
            { 223, 16, LEO2_FIELD_GF8, 65 },
            { 240, 16, LEO2_FIELD_GF8, 64 },
            { 193, 65, LEO2_FIELD_GF16, 128 },
            { 257, 32, LEO2_FIELD_GF16, 128 },
            { 511, 64, LEO2_FIELD_GF16, 128 },
            { 1000, 200, LEO2_FIELD_GF16, 128 }
        };
        uint64_t parity_bytes = 0;
        uint64_t restored_bytes = 0;
        for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i)
            run_case(context, cases[i], i, parity_bytes, restored_bytes);
        leo2_context_destroy(context);
        std::cout << "PASS high_pruned_legacy cases="
                  << sizeof(cases) / sizeof(cases[0])
                  << " parity_bytes=" << parity_bytes
                  << " restored_bytes=" << restored_bytes
                  << " workspaces=materialized,tiled" << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "FAIL high_pruned_legacy: " << error.what() << std::endl;
        return 1;
    }
}
