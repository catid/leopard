/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in the Leopard-RS
    BSD license are met.  See the repository LICENSE file.
*/

#include <leopard2.h>

#include <stdint.h>

#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require_result(leo2_result result, const char* operation)
{
    if (result != LEO2_SUCCESS)
        throw std::runtime_error(
            std::string(operation) + ": " + leo2_result_string(result));
}

void* aligned_data(std::vector<unsigned char>& storage, size_t alignment)
{
    if (storage.empty())
        return NULL;
    const uintptr_t address = reinterpret_cast<uintptr_t>(&storage[0]);
    const uintptr_t aligned = (address + alignment - 1) & ~(alignment - 1);
    return reinterpret_cast<void*>(aligned);
}

uint64_t hash_bytes(uint64_t hash, const void* data, size_t bytes)
{
    const unsigned char* input = static_cast<const unsigned char*>(data);
    for (size_t i = 0; i < bytes; ++i)
    {
        hash ^= input[i];
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

struct Case
{
    uint32_t k;
    uint32_t r;
    leo2_profile profile;
    leo2_field field;
    uint64_t bytes;
};

void run_case(leo2_context* context, const Case& test_case, unsigned index)
{
    leo2_codec* codec = NULL;
    require_result(leo2_codec_create(context,
        test_case.k,
        test_case.r,
        test_case.profile,
        test_case.field,
        NULL,
        &codec), "leo2_codec_create");

    size_t encode_scratch_bytes = 0;
    size_t one_shot_decode_scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(
        codec, test_case.bytes, &encode_scratch_bytes),
        "leo2_encode_scratch_size");
    require_result(leo2_decode_scratch_size(
        codec, test_case.bytes, &one_shot_decode_scratch_bytes),
        "leo2_decode_scratch_size");

    std::vector<std::vector<unsigned char> > originals(
        test_case.k, std::vector<unsigned char>(test_case.bytes));
    std::vector<std::vector<unsigned char> > recovery(
        test_case.r, std::vector<unsigned char>(test_case.bytes));
    std::vector<const void*> original_ptrs(test_case.k);
    std::vector<void*> recovery_ptrs(test_case.r);
    for (uint32_t shard = 0; shard < test_case.k; ++shard)
    {
        for (uint64_t byte = 0; byte < test_case.bytes; ++byte)
        {
            originals[shard][static_cast<size_t>(byte)] =
                static_cast<unsigned char>(
                    17U * index + 29U * shard + 11U * byte +
                    (byte >> 2));
        }
        original_ptrs[shard] = &originals[shard][0];
    }
    for (uint32_t shard = 0; shard < test_case.r; ++shard)
        recovery_ptrs[shard] = &recovery[shard][0];

    const size_t alignment = leo2_scratch_alignment();
    std::vector<unsigned char> encode_scratch_storage(
        encode_scratch_bytes == 0 ? 0 : encode_scratch_bytes + alignment - 1);
    void* encode_scratch = aligned_data(encode_scratch_storage, alignment);
    require_result(leo2_encode(codec,
        test_case.bytes,
        &original_ptrs[0],
        &recovery_ptrs[0],
        encode_scratch,
        encode_scratch_bytes), "leo2_encode");

    uint64_t parity_hash = UINT64_C(1469598103934665603);
    for (uint32_t shard = 0; shard < test_case.r; ++shard)
        parity_hash = hash_bytes(
            parity_hash, &recovery[shard][0], static_cast<size_t>(test_case.bytes));

    std::vector<uint8_t> original_present(test_case.k, 1);
    std::vector<uint8_t> recovery_present(test_case.r, 1);
    original_present[0] = 0;
    original_ptrs[0] = NULL;

    leo2_decode_plan* plan = NULL;
    require_result(leo2_decode_plan_create(codec,
        &original_present[0],
        &recovery_present[0],
        &plan), "leo2_decode_plan_create");
    size_t plan_scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(
        plan, test_case.bytes, &plan_scratch_bytes),
        "leo2_decode_plan_scratch_size");

    std::vector<unsigned char> restored(test_case.bytes, 0);
    std::vector<void*> restored_ptrs(test_case.k, NULL);
    std::vector<const void*> recovery_inputs(test_case.r);
    restored_ptrs[0] = &restored[0];
    for (uint32_t shard = 0; shard < test_case.r; ++shard)
        recovery_inputs[shard] = &recovery[shard][0];
    std::vector<unsigned char> decode_scratch_storage(
        plan_scratch_bytes == 0 ? 0 : plan_scratch_bytes + alignment - 1);
    void* decode_scratch = aligned_data(decode_scratch_storage, alignment);
    require_result(leo2_decode_plan_execute(plan,
        test_case.bytes,
        &original_ptrs[0],
        &recovery_inputs[0],
        &restored_ptrs[0],
        decode_scratch,
        plan_scratch_bytes), "leo2_decode_plan_execute");
    if (std::memcmp(&restored[0], &originals[0][0],
            static_cast<size_t>(test_case.bytes)) != 0)
        throw std::runtime_error("restored original mismatch");

    const uint64_t restored_hash = hash_bytes(
        UINT64_C(1469598103934665603),
        &restored[0], static_cast<size_t>(test_case.bytes));
    std::cout << "case=" << index
              << ",profile=" << static_cast<unsigned>(leo2_codec_profile(codec))
              << ",field=" << static_cast<unsigned>(leo2_codec_field(codec))
              << ",parent=" << leo2_codec_parent_count(codec)
              << ",padded=" << leo2_codec_padded_side(codec)
              << ",encode_scratch=" << encode_scratch_bytes
              << ",decode_scratch=" << one_shot_decode_scratch_bytes
              << ",plan_scratch=" << plan_scratch_bytes
              << ",parity_hash=" << parity_hash
              << ",restored_hash=" << restored_hash << '\n';

    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
}

} // namespace

int main()
{
    try
    {
        leo2_context_options options;
        std::memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_SCALAR;
        options.thread_count = 1;
        leo2_context* context = NULL;
        require_result(leo2_context_create(&options, &context),
            "leo2_context_create");

        std::cout << "alignment=" << leo2_scratch_alignment()
                  << ",field_mask=" << leo2_context_field_mask(context)
                  << ",backend="
                  << static_cast<unsigned>(leo2_context_backend(context))
                  << '\n';
        const Case cases[] = {
            { 8, 3, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 31 },
            { 3, 5, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8, 17 },
            { 9, 3, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, 34 },
            { 3, 5, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16, 32 }
        };
        for (unsigned i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i)
            run_case(context, cases[i], i);
        leo2_context_destroy(context);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
