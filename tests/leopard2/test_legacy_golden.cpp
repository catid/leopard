/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in the Leopard-RS
    BSD license are met.  See the repository LICENSE file.
*/

#include "leopard2.h"
#include "legacy_golden_vectors.h"

#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

namespace {

class AlignedBuffer
{
public:
    explicit AlignedBuffer(size_t bytes)
        : pointer_(NULL)
        , bytes_(bytes)
    {
        if (bytes == 0)
            return;
#if defined(_MSC_VER)
        pointer_ = _aligned_malloc(bytes, leo2_scratch_alignment());
#else
        if (posix_memalign(&pointer_, leo2_scratch_alignment(), bytes) != 0)
            pointer_ = NULL;
#endif
    }

    ~AlignedBuffer()
    {
#if defined(_MSC_VER)
        _aligned_free(pointer_);
#else
        free(pointer_);
#endif
    }

    bool valid() const { return bytes_ == 0 || pointer_ != NULL; }
    void* data() { return pointer_; }
    size_t bytes() const { return bytes_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);

    void* pointer_;
    size_t bytes_;
};

static void Fail(const char* vector_name, const char* operation)
{
    fprintf(stderr, "legacy golden vector %s failed: %s\n", vector_name, operation);
    abort();
}

static void Require(bool condition, const char* vector_name, const char* operation)
{
    if (!condition)
        Fail(vector_name, operation);
}

static void RequireResult(
    leo2_result result,
    const char* vector_name,
    const char* operation)
{
    if (result != LEO2_SUCCESS)
    {
        fprintf(stderr, "legacy golden vector %s failed: %s: %s (%d)\n",
            vector_name, operation, leo2_result_string(result),
            static_cast<int>(result));
        abort();
    }
}

static uint8_t NextInputByte(uint64_t& state)
{
    state += UINT64_C(0x9e3779b97f4a7c15);
    uint64_t value = state;
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    value ^= value >> 31;
    return static_cast<uint8_t>(value >> 56);
}

static int HexValue(char value)
{
    if (value >= '0' && value <= '9')
        return value - '0';
    if (value >= 'a' && value <= 'f')
        return value - 'a' + 10;
    if (value >= 'A' && value <= 'F')
        return value - 'A' + 10;
    return -1;
}

static std::vector<uint8_t> DecodeParity(
    const leopard2_legacy_golden::GoldenVector& vector)
{
    const size_t parity_bytes = static_cast<size_t>(vector.recovery_count) *
        vector.shard_bytes;
    Require(strlen(vector.parity_hex) == parity_bytes * 2,
        vector.name, "artifact parity length");
    std::vector<uint8_t> parity(parity_bytes);
    for (size_t i = 0; i < parity_bytes; ++i)
    {
        const int high = HexValue(vector.parity_hex[i * 2]);
        const int low = HexValue(vector.parity_hex[i * 2 + 1]);
        Require(high >= 0 && low >= 0, vector.name, "artifact parity hex digit");
        parity[i] = static_cast<uint8_t>((high << 4) | low);
    }
    return parity;
}

static void RunVector(
    leo2_context* context,
    const leopard2_legacy_golden::GoldenVector& vector,
    uint64_t& parity_bytes_checked,
    uint64_t& recovered_bytes_checked)
{
    const leo2_field field = vector.field_bits == 8
        ? LEO2_FIELD_GF8 : LEO2_FIELD_GF16;
    Require(vector.field_bits == 8 || vector.field_bits == 16,
        vector.name, "artifact field");
    Require(vector.missing_count <= vector.recovery_count,
        vector.name, "artifact loss count");

    leo2_codec* codec = NULL;
    RequireResult(leo2_codec_create(
        context,
        vector.original_count,
        vector.recovery_count,
        LEO2_PROFILE_LEGACY_HIGH_V1,
        field,
        NULL,
        &codec), vector.name, "codec create");
    Require(codec != NULL, vector.name, "codec pointer");
    Require(leo2_codec_profile(codec) == LEO2_PROFILE_LEGACY_HIGH_V1,
        vector.name, "codec profile");
    Require(leo2_codec_field(codec) == field, vector.name, "codec field");
    if (field == LEO2_FIELD_GF16)
        Require(leo2_codec_parent_count(codec) > 256,
            vector.name, "GF16 parent boundary");

    std::vector<std::vector<uint8_t> > original(
        vector.original_count, std::vector<uint8_t>(vector.shard_bytes));
    std::vector<const void*> original_input(vector.original_count);
    uint64_t state = vector.seed;
    for (uint32_t shard = 0; shard < vector.original_count; ++shard)
    {
        for (uint32_t offset = 0; offset < vector.shard_bytes; ++offset)
            original[shard][offset] = NextInputByte(state);
        original_input[shard] = &original[shard][0];
    }

    const std::vector<uint8_t> golden_parity = DecodeParity(vector);
    std::vector<uint8_t> encoded_parity(golden_parity.size(), 0);
    std::vector<void*> encoded_output(vector.recovery_count);
    for (uint32_t shard = 0; shard < vector.recovery_count; ++shard)
        encoded_output[shard] = &encoded_parity[
            static_cast<size_t>(shard) * vector.shard_bytes];

    size_t encode_scratch_bytes = 0;
    RequireResult(leo2_encode_scratch_size(
        codec, vector.shard_bytes, &encode_scratch_bytes),
        vector.name, "encode scratch query");
    AlignedBuffer encode_scratch(encode_scratch_bytes);
    Require(encode_scratch.valid(), vector.name, "encode scratch allocation");
    RequireResult(leo2_encode(
        codec,
        vector.shard_bytes,
        &original_input[0],
        &encoded_output[0],
        encode_scratch.data(),
        encode_scratch.bytes()), vector.name, "encode");
    Require(encoded_parity == golden_parity,
        vector.name, "raw parity differs from frozen old data");
    parity_bytes_checked += static_cast<uint64_t>(golden_parity.size());

    std::vector<uint8_t> original_present(vector.original_count, 1);
    std::vector<uint8_t> recovery_present(vector.recovery_count, 1);
    for (size_t i = 0; i < vector.missing_count; ++i)
    {
        const uint32_t index = vector.missing_originals[i];
        Require(index < vector.original_count, vector.name, "artifact loss index");
        Require(original_present[index] != 0, vector.name, "duplicate loss index");
        original_present[index] = 0;
    }

    leo2_decode_plan* plan = NULL;
    RequireResult(leo2_decode_plan_create(codec,
        &original_present[0], &recovery_present[0], &plan),
        vector.name, "decode plan create");
    Require(plan != NULL, vector.name, "decode plan pointer");

    size_t decode_scratch_bytes = 0;
    RequireResult(leo2_decode_plan_scratch_size(
        plan, vector.shard_bytes, &decode_scratch_bytes),
        vector.name, "decode scratch query");
    AlignedBuffer decode_scratch(decode_scratch_bytes);
    Require(decode_scratch.valid(), vector.name, "decode scratch allocation");

    std::vector<const void*> received_original(vector.original_count, NULL);
    std::vector<const void*> received_recovery(vector.recovery_count, NULL);
    std::vector<std::vector<uint8_t> > restored(
        vector.original_count,
        std::vector<uint8_t>(vector.shard_bytes, static_cast<uint8_t>(0xa5)));
    std::vector<void*> restored_output(vector.original_count, NULL);
    for (uint32_t shard = 0; shard < vector.original_count; ++shard)
    {
        if (original_present[shard])
            received_original[shard] = &original[shard][0];
        else
            restored_output[shard] = &restored[shard][0];
    }
    for (uint32_t shard = 0; shard < vector.recovery_count; ++shard)
        received_recovery[shard] = &golden_parity[
            static_cast<size_t>(shard) * vector.shard_bytes];

    for (unsigned repeat = 0; repeat < 2; ++repeat)
    {
        for (size_t i = 0; i < vector.missing_count; ++i)
            std::fill(restored[vector.missing_originals[i]].begin(),
                restored[vector.missing_originals[i]].end(),
                static_cast<uint8_t>(0xa5));
        RequireResult(leo2_decode_plan_execute(
            plan,
            vector.shard_bytes,
            &received_original[0],
            &received_recovery[0],
            &restored_output[0],
            decode_scratch.data(),
            decode_scratch.bytes()), vector.name, "decode frozen old data");
        for (size_t i = 0; i < vector.missing_count; ++i)
        {
            const uint32_t index = vector.missing_originals[i];
            Require(restored[index] == original[index],
                vector.name, "recovered original bytes");
            recovered_bytes_checked += vector.shard_bytes;
        }
    }

    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
}

} // namespace

int main()
{
    leo2_context_options options;
    memset(&options, 0, sizeof(options));
    options.struct_size = sizeof(options);
    options.thread_count = 1;
    leo2_context* context = NULL;
    RequireResult(leo2_context_create(&options, &context),
        "context", "context create");

    uint64_t parity_bytes_checked = 0;
    uint64_t recovered_bytes_checked = 0;
    for (size_t i = 0;
        i < leopard2_legacy_golden::kGoldenVectorCount;
        ++i)
    {
        RunVector(context, leopard2_legacy_golden::kGoldenVectors[i],
            parity_bytes_checked, recovered_bytes_checked);
    }
    leo2_context_destroy(context);

    printf("leopard2 legacy raw golden vectors passed: vectors=%llu "
           "parity_bytes=%llu recovered_bytes=%llu\n",
        static_cast<unsigned long long>(
            leopard2_legacy_golden::kGoldenVectorCount),
        static_cast<unsigned long long>(parity_bytes_checked),
        static_cast<unsigned long long>(recovered_bytes_checked));
    return 0;
}
