/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in the Leopard-RS
    BSD license are met.  See the repository LICENSE file.
*/

#include "Leopard2Direct.h"
#include "direct_oracle.h"

#include <atomic>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

#ifndef LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE
#define LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE 0
#endif

namespace {

class AlignedBytes
{
public:
    explicit AlignedBytes(size_t bytes) : data_(NULL), bytes_(bytes)
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
    size_t size() const { return bytes_; }

private:
    AlignedBytes(const AlignedBytes&);
    AlignedBytes& operator=(const AlignedBytes&);
    void* data_;
    size_t bytes_;
};

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void require_result(leo2_result result, const std::string& operation)
{
    if (result != LEO2_SUCCESS)
        throw std::runtime_error(operation + ": " + leo2_result_string(result));
}

struct PolicyCase
{
    const char* name;
    leo2_backend backend;
    leo2_field field;
    leo2_profile profile;
    uint32_t k;
    uint32_t r;
    uint32_t losses;
    bool expect_skip;
};

std::vector<uint32_t> choose_missing(uint32_t k, uint32_t count)
{
    std::vector<uint32_t> result;
    std::vector<uint8_t> selected(k, 0);
    const uint32_t stride = k % 5U == 0 ? 7U : 5U;
    for (uint32_t step = 0; result.size() < count; ++step)
    {
        const uint32_t index = (step * stride + 3U) % k;
        if (!selected[index])
        {
            selected[index] = 1;
            result.push_back(index);
        }
    }
    return result;
}

leo2_context* create_context(leo2_backend backend)
{
    leo2_context_options options;
    memset(&options, 0, sizeof(options));
    options.struct_size = sizeof(options);
    options.backend = backend;
    options.thread_count = 1;
    leo2_context* context = NULL;
    const leo2_result result = leo2_context_create(&options, &context);
    if (result == LEO2_UNSUPPORTED)
        return NULL;
    require_result(result, "context create");
    return context;
}

void check_direct_parity(
    const PolicyCase& test,
    const std::vector<const void*>& original,
    const std::vector<void*>& recovery,
    size_t bytes)
{
    if (test.field != LEO2_FIELD_GF8)
        return;
    const leopard2_test::BinaryField field = leopard2_test::make_legacy_gf8();
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            test.profile == LEO2_PROFILE_LEGACY_HIGH_V1
                ? leopard2_test::kLegacyHigh : leopard2_test::kLow,
            test.k, test.r);
    std::vector<leopard2_test::Element> systematic_points(
        layout.parent_dimension);
    for (uint32_t i = 0; i < layout.parent_dimension; ++i)
        systematic_points[i] = static_cast<leopard2_test::Element>(
            layout.systematic_coordinates[i]);
    std::vector<leopard2_test::Element> values(
        layout.parent_dimension, 0);

    // Direct interpolation is intentionally slow and algebraically independent.
    // Sampling the beginning, middle, and vector-tail bytes covers byte handling
    // without constructing an O(K*R) generator matrix for this policy test.
    const size_t sample_bytes[] = { 0, bytes / 2, bytes - 1 };
    for (size_t sample = 0;
         sample < sizeof(sample_bytes) / sizeof(sample_bytes[0]); ++sample)
    {
        const size_t byte = sample_bytes[sample];
        if (sample != 0 && byte == sample_bytes[sample - 1])
            continue;
        for (uint32_t i = 0; i < test.k; ++i)
            values[i] = static_cast<const uint8_t*>(original[i])[byte];
        for (uint32_t i = 0; i < test.r; ++i)
        {
            const leopard2_test::Element expected =
                leopard2_test::lagrange_evaluate(field, systematic_points,
                    values, static_cast<leopard2_test::Element>(
                        layout.parity_coordinates[i]));
            require(static_cast<uint8_t*>(recovery[i])[byte] ==
                    expected,
                std::string(test.name) + " direct parity mismatch");
        }
    }
}

void execute_and_check(
    const PolicyCase& test,
    const leo2_decode_plan* plan,
    const std::vector<const void*>& decode_original,
    const std::vector<const void*>& expected_original,
    const std::vector<const void*>& recovery,
    const std::vector<uint32_t>& missing,
    size_t bytes,
    bool concurrent)
{
    const unsigned worker_count = concurrent ? 4U : 1U;
    size_t scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(
        plan, bytes, &scratch_bytes), "decode scratch query");

    std::vector<AlignedBytes*> scratch(worker_count, NULL);
    std::vector<AlignedBytes*> restored_storage(worker_count, NULL);
    std::vector<std::vector<void*> > restored(
        worker_count, std::vector<void*>(test.k, NULL));
    for (unsigned worker = 0; worker < worker_count; ++worker)
    {
        scratch[worker] = new AlignedBytes(scratch_bytes);
        restored_storage[worker] = new AlignedBytes(
            static_cast<size_t>(missing.size()) * bytes);
        for (size_t i = 0; i < missing.size(); ++i)
            restored[worker][missing[i]] =
                restored_storage[worker]->data() + i * bytes;
    }

    std::atomic<bool> success(true);
    const auto run = [&](unsigned worker) {
        if (leo2_decode_plan_execute(plan, bytes, decode_original.data(),
                recovery.data(), restored[worker].data(),
                scratch[worker]->data(), scratch[worker]->size()) !=
            LEO2_SUCCESS)
            success.store(false, std::memory_order_relaxed);
    };
    std::vector<std::thread> threads;
    for (unsigned worker = 0; worker < worker_count; ++worker)
        threads.push_back(std::thread(run, worker));
    for (size_t i = 0; i < threads.size(); ++i)
        threads[i].join();
    require(success.load(std::memory_order_relaxed),
        std::string(test.name) + " concurrent decode failed");

    for (unsigned worker = 0; worker < worker_count; ++worker)
    {
        for (size_t i = 0; i < missing.size(); ++i)
        {
            require(memcmp(restored[worker][missing[i]],
                    expected_original[missing[i]], bytes) == 0,
                std::string(test.name) + " restored original mismatch");
        }
        delete restored_storage[worker];
        delete scratch[worker];
    }
}

bool run_case(const PolicyCase& test)
{
    leo2_context* context = create_context(test.backend);
    if (!context)
        return false;
    const bool avx2_context =
        leo2_context_backend(context) == LEO2_BACKEND_AVX2;
    const bool production_skip = test.expect_skip &&
        (test.backend != LEO2_BACKEND_AUTO || avx2_context);
#if LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE != 0
    const bool experimental_skip = avx2_context &&
        test.profile == LEO2_PROFILE_LEGACY_HIGH_V1 &&
        test.field == LEO2_FIELD_GF8 &&
        test.k >= 5 && test.k <= 16 &&
        test.r >= 5 && test.r <= 8 &&
        test.losses >= 5 && test.losses <= 8;
#else
    const bool experimental_skip = false;
#endif
    const bool expect_skip = production_skip || experimental_skip;
    if ((test.field == LEO2_FIELD_GF8 &&
            (leo2_context_field_mask(context) & LEO2_FIELD_MASK_GF8) == 0) ||
        (test.field == LEO2_FIELD_GF16 &&
            (leo2_context_field_mask(context) & LEO2_FIELD_MASK_GF16) == 0))
    {
        leo2_context_destroy(context);
        return false;
    }

    leo2_codec_options codec_options;
    memset(&codec_options, 0, sizeof(codec_options));
    codec_options.struct_size = sizeof(codec_options);
    if (strcmp(test.name, "negative-t4") == 0)
        codec_options.flags = LEO2_CODEC_FORCE_SPECIALIZED_DECODE;
    leo2_codec* codec = NULL;
    require_result(leo2_codec_create(context, test.k, test.r,
        test.profile, test.field, &codec_options, &codec), "codec create");

    const size_t bytes = test.field == LEO2_FIELD_GF8 ? 65U : 66U;
    AlignedBytes original_storage(static_cast<size_t>(test.k) * bytes);
    AlignedBytes recovery_storage(static_cast<size_t>(test.r) * bytes);
    std::vector<const void*> original(test.k);
    std::vector<void*> recovery(test.r);
    for (uint32_t i = 0; i < test.k; ++i)
    {
        uint8_t* shard = original_storage.data() + static_cast<size_t>(i) * bytes;
        for (size_t byte = 0; byte < bytes; ++byte)
            shard[byte] = static_cast<uint8_t>(i * 29U + byte * 17U + 11U);
        original[i] = shard;
    }
    for (uint32_t i = 0; i < test.r; ++i)
        recovery[i] = recovery_storage.data() + static_cast<size_t>(i) * bytes;

    size_t encode_scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(
        codec, bytes, &encode_scratch_bytes), "encode scratch query");
    AlignedBytes encode_scratch(encode_scratch_bytes);
    require_result(leo2_encode(codec, bytes, original.data(), recovery.data(),
        encode_scratch.data(), encode_scratch.size()), "encode");
    // Keep an independent algebraic oracle on every explicitly selected shape
    // covered by the production policy.  Control cases exercise only policy
    // boundaries and would make this focused test unnecessarily expensive.
    if (expect_skip && avx2_context)
        check_direct_parity(test, original, recovery, bytes);

    const std::vector<uint32_t> missing =
        choose_missing(test.k, test.losses);
    std::vector<uint8_t> original_present(test.k, 1);
    std::vector<uint8_t> recovery_present(test.r, 1);
    std::vector<const void*> decode_original(original);
    std::vector<const void*> decode_recovery(test.r);
    for (uint32_t i = 0; i < test.r; ++i)
        decode_recovery[i] = recovery[i];
    for (size_t i = 0; i < missing.size(); ++i)
    {
        original_present[missing[i]] = 0;
        decode_original[missing[i]] = NULL;
    }

    leo2_decode_plan* plan = NULL;
    require_result(leo2_decode_plan_create(codec, original_present.data(),
        recovery_present.data(), &plan), "decode plan create");
    leopard2_internal::DecodePlanPrunedScheduleInfo info;
    require(leopard2_internal::GetDecodePlanPrunedScheduleInfo(plan, &info),
        std::string(test.name) + " schedule introspection failed");
    const size_t high_count =
        info.high_input_plan_count + info.high_output_plan_count;
    const size_t low_count =
        info.low_input_plan_count + info.low_output_plan_count;
    if (expect_skip)
    {
        require(high_count == 0 && low_count == 0,
            std::string(test.name) + " retained a pruned schedule");
    }
    else if (test.profile == LEO2_PROFILE_LEGACY_HIGH_V1)
        require(high_count != 0,
            std::string(test.name) + " did not compile its control schedule");
    else
        require(low_count != 0,
            std::string(test.name) + " did not compile its low-profile schedule");

    execute_and_check(test, plan, decode_original, original, decode_recovery,
        missing, bytes, expect_skip);
    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
    leo2_context_destroy(context);
    return true;
}

} // namespace

int main()
{
    const PolicyCase cases[] = {
        { "positive-t8-k2r", LEO2_BACKEND_AVX2, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 16, 8, 8, true },
        { "positive-t64-k2r", LEO2_BACKEND_AVX2, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 128, 64, 64, true },
        { "positive-t64-k3r", LEO2_BACKEND_AVX2, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 192, 64, 64, true },
        { "positive-t8-k2r-minus-one", LEO2_BACKEND_AVX2, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 15, 8, 8, true },
        { "positive-t8-r-minus-one", LEO2_BACKEND_AVX2, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 16, 7, 7, true },
        { "positive-t16-k2r", LEO2_BACKEND_AVX2, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 32, 16, 16, true },
        { "positive-t16-r-minus-one", LEO2_BACKEND_AVX2, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 32, 15, 15, true },
        { "positive-t64-k2r-minus-one", LEO2_BACKEND_AVX2,
          LEO2_FIELD_GF8, LEO2_PROFILE_LEGACY_HIGH_V1, 127, 64, 64, true },
        { "positive-t64-r-minus-one", LEO2_BACKEND_AVX2,
          LEO2_FIELD_GF8, LEO2_PROFILE_LEGACY_HIGH_V1, 128, 63, 63, true },
        { "negative-t4", LEO2_BACKEND_AVX2, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 8, 4, 4, false },
        { "negative-t32", LEO2_BACKEND_AVX2, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 64, 32, 32, false },
        { "negative-t32-r-minus-one", LEO2_BACKEND_AVX2, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 64, 31, 31, false },
        { "negative-k-below-boundary", LEO2_BACKEND_AVX2,
          LEO2_FIELD_GF8, LEO2_PROFILE_LEGACY_HIGH_V1, 14, 8, 8, false },
        { "negative-loss-below-r", LEO2_BACKEND_AVX2, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 16, 8, 7, false },
        { "negative-r-two-below-t", LEO2_BACKEND_AVX2, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 16, 6, 6, false },
        { "positive-auto-avx2", LEO2_BACKEND_AUTO, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 16, 8, 8, true },
        { "negative-scalar", LEO2_BACKEND_SCALAR, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 16, 8, 8, false },
        { "negative-gf16", LEO2_BACKEND_AVX2, LEO2_FIELD_GF16,
          LEO2_PROFILE_LEGACY_HIGH_V1, 16, 8, 8, false },
        { "negative-low-profile", LEO2_BACKEND_AVX2, LEO2_FIELD_GF8,
          LEO2_PROFILE_LOW_V1, 16, 8, 8, false }
    };

    unsigned executed = 0;
    unsigned skipped = 0;
    try
    {
        for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i)
        {
            if (run_case(cases[i]))
                ++executed;
            else
                ++skipped;
        }
        std::cout << "PASS dense_plan_policy cases=" << executed
                  << " unsupported=" << skipped
                  << " concurrent_plan_workers=12" << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "FAIL dense_plan_policy: " << error.what() << std::endl;
        return 1;
    }
}
