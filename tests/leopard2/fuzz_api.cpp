/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in the Leopard-RS
    BSD license are met.  See the repository LICENSE file.
*/

#include "leopard2.h"

#include <algorithm>
#include <stdint.h>
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

class Random
{
public:
    explicit Random(uint64_t seed)
        : state_(seed ? seed : UINT64_C(0x9e3779b97f4a7c15))
    {}

    uint64_t next()
    {
        uint64_t value = state_;
        value ^= value >> 12;
        value ^= value << 25;
        value ^= value >> 27;
        state_ = value;
        return value * UINT64_C(2685821657736338717);
    }

private:
    uint64_t state_;
};

static void CheckImpl(bool condition, int line)
{
    if (!condition)
    {
        fprintf(stderr, "Leopard2 API fuzz invariant failed at line %d\n", line);
        abort();
    }
}

#define Check(condition) CheckImpl((condition), __LINE__)

static uint64_t HashInput(const uint8_t* data, size_t size)
{
    uint64_t hash = UINT64_C(1469598103934665603);
    for (size_t i = 0; i < size; ++i)
    {
        hash ^= data[i];
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

} // namespace

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size)
{
    if (!data || size < 6)
        return 0;

    try
    {
        const uint32_t k = 1u + data[0] % 32u;
        const uint32_t r = 1u + data[1] % 32u;
        size_t shard_bytes = 1u +
            (static_cast<size_t>(data[2]) |
             (static_cast<size_t>(data[3]) << 8)) % 1025u;
        const leo2_profile profile = (data[4] & 1u)
            ? LEO2_PROFILE_LEGACY_HIGH_V1 : LEO2_PROFILE_LOW_V1;

        leo2_context_options context_options;
        memset(&context_options, 0, sizeof(context_options));
        context_options.struct_size = sizeof(context_options);
        context_options.thread_count = 1;
        const leo2_backend requested_backends[] = {
            LEO2_BACKEND_AUTO,
            LEO2_BACKEND_SCALAR,
            LEO2_BACKEND_SSSE3,
            LEO2_BACKEND_AVX2,
            LEO2_BACKEND_AVX512,
            LEO2_BACKEND_GFNI
        };
        const leo2_backend requested_backend =
            requested_backends[(data[4] >> 2) % 6u];
        context_options.backend = requested_backend;
        leo2_context* context = NULL;
        leo2_result context_result = leo2_context_create(
            &context_options, &context);
        bool backend_fell_back = false;
        if (context_result == LEO2_UNSUPPORTED)
        {
            backend_fell_back = true;
            context_options.backend = LEO2_BACKEND_AUTO;
            context_result = leo2_context_create(&context_options, &context);
        }
        Check(context_result == LEO2_SUCCESS);
        Check(context != NULL);
        if (!backend_fell_back && requested_backend != LEO2_BACKEND_AUTO)
            Check(leo2_context_backend(context) == requested_backend);
        const unsigned field_mask = leo2_context_field_mask(context);
        Check(field_mask != 0);
        const leo2_field field = (field_mask & LEO2_FIELD_MASK_GF8)
            ? LEO2_FIELD_GF8 : LEO2_FIELD_GF16;
        if (field == LEO2_FIELD_GF16)
            shard_bytes = (shard_bytes + 1u) & ~static_cast<size_t>(1u);

        leo2_codec_options codec_options;
        memset(&codec_options, 0, sizeof(codec_options));
        codec_options.struct_size = sizeof(codec_options);
        if (data[4] & 2u)
            codec_options.flags = LEO2_CODEC_FORCE_GENERIC_DECODE;
        leo2_codec* codec = NULL;
        const leo2_result codec_result = leo2_codec_create(
            context, k, r, profile, field, &codec_options, &codec);
        Check(codec_result == LEO2_SUCCESS);
        Check(codec != NULL);

        Random random(HashInput(data, size));
        std::vector<std::vector<uint8_t> > original(
            k, std::vector<uint8_t>(shard_bytes));
        std::vector<std::vector<uint8_t> > recovery(
            r, std::vector<uint8_t>(shard_bytes));
        std::vector<const void*> original_input(k);
        std::vector<void*> recovery_output(r);
        for (uint32_t i = 0; i < k; ++i)
        {
            for (size_t j = 0; j < shard_bytes; ++j)
                original[i][j] = static_cast<uint8_t>(random.next());
            original_input[i] = &original[i][0];
        }
        for (uint32_t i = 0; i < r; ++i)
            recovery_output[i] = &recovery[i][0];

        size_t encode_scratch_bytes = 0;
        Check(leo2_encode_scratch_size(
            codec, shard_bytes, &encode_scratch_bytes) == LEO2_SUCCESS);
        AlignedBuffer encode_scratch(encode_scratch_bytes);
        Check(encode_scratch.valid());
        Check(leo2_encode(codec, shard_bytes, &original_input[0],
            &recovery_output[0], encode_scratch.data(),
            encode_scratch.bytes()) == LEO2_SUCCESS);

        std::vector<uint32_t> coordinates(k + r);
        for (uint32_t i = 0; i < k + r; ++i)
            coordinates[i] = i;
        for (size_t i = coordinates.size(); i > 1; --i)
            std::swap(coordinates[i - 1], coordinates[random.next() % i]);
        const uint32_t loss_count = data[5] % (r + 1u);
        std::vector<uint8_t> original_present(k, 1);
        std::vector<uint8_t> recovery_present(r, 1);
        for (uint32_t i = 0; i < loss_count; ++i)
        {
            if (coordinates[i] < k)
                original_present[coordinates[i]] = 0;
            else
                recovery_present[coordinates[i] - k] = 0;
        }

        leo2_decode_plan* plan = NULL;
        Check(leo2_decode_plan_create(codec, &original_present[0],
            &recovery_present[0], &plan) == LEO2_SUCCESS);
        size_t decode_scratch_bytes = 0;
        Check(leo2_decode_plan_scratch_size(
            plan, shard_bytes, &decode_scratch_bytes) == LEO2_SUCCESS);
        AlignedBuffer decode_scratch(decode_scratch_bytes);
        Check(decode_scratch.valid());

        std::vector<const void*> received_original(k);
        std::vector<const void*> received_recovery(r);
        std::vector<std::vector<uint8_t> > restored(
            k, std::vector<uint8_t>(shard_bytes, 0));
        std::vector<void*> restored_output(k, NULL);
        for (uint32_t i = 0; i < k; ++i)
        {
            if (original_present[i])
                received_original[i] = &original[i][0];
            else
                restored_output[i] = &restored[i][0];
        }
        for (uint32_t i = 0; i < r; ++i)
            if (recovery_present[i])
                received_recovery[i] = &recovery[i][0];

        for (unsigned repeat = 0; repeat < 2; ++repeat)
        {
            Check(leo2_decode_plan_execute(plan, shard_bytes,
                &received_original[0], &received_recovery[0],
                &restored_output[0], decode_scratch.data(),
                decode_scratch.bytes()) == LEO2_SUCCESS);
            for (uint32_t i = 0; i < k; ++i)
                if (!original_present[i])
                    Check(restored[i] == original[i]);
        }

        leo2_decode_plan_destroy(plan);
        leo2_codec_destroy(codec);
        leo2_context_destroy(context);
    }
    catch (...)
    {
        // The generated dimensions are deliberately bounded.  An exception
        // is therefore an unexpected harness or codec failure, not an input
        // rejection to hide from libFuzzer.
        abort();
    }
    return 0;
}
