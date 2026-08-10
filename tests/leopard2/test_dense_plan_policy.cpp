/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in the Leopard-RS
    BSD license are met.  See the repository LICENSE file.
*/

#include "Leopard2Direct.h"
#include "Leopard2Dispatch.h"
#include "direct_oracle.h"

#include <algorithm>
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

#ifndef LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK
#define LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK 0
#endif

namespace {

const bool kExpectSmallDualRegularFallback =
    LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK != 0;

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

void execute_one_shot_and_check(
    const PolicyCase& test,
    const leo2_codec* codec,
    const std::vector<uint8_t>& original_present,
    const std::vector<uint8_t>& recovery_present,
    const std::vector<const void*>& decode_original,
    const std::vector<const void*>& expected_original,
    const std::vector<const void*>& recovery,
    const std::vector<uint32_t>& missing,
    size_t bytes,
    unsigned setup_mode)
{
    require(leopard2_internal::SetOneShotPlanSetupModeForDiagnostics(
            setup_mode),
        std::string(test.name) + " cannot select one-shot setup mode");
    size_t scratch_bytes = 0;
    require_result(leo2_decode_scratch_size(codec, bytes, &scratch_bytes),
        "one-shot scratch query");
    AlignedBytes scratch(scratch_bytes);
    AlignedBytes restored_storage(
        static_cast<size_t>(missing.size()) * bytes);
    std::vector<void*> restored(test.k, NULL);
    for (size_t i = 0; i < missing.size(); ++i)
        restored[missing[i]] = restored_storage.data() + i * bytes;

    require_result(leo2_decode(codec, bytes, original_present.data(),
        recovery_present.data(), decode_original.data(), recovery.data(),
        restored.data(), scratch.data(), scratch.size()),
        setup_mode == 0 ? "broad-metadata one-shot decode" :
                          "exact-byte one-shot decode");
    for (size_t i = 0; i < missing.size(); ++i)
    {
        require(memcmp(restored[missing[i]], expected_original[missing[i]],
                bytes) == 0,
            std::string(test.name) +
                " one-shot restored original mismatch in setup mode " +
                static_cast<char>('0' + setup_mode));
    }

    leo2_decode_plan* transient_plan = NULL;
    const leo2_result transient_result =
        leopard2_internal::CreateOneShotTransformPlanForDiagnostics(
            codec, bytes, original_present.data(), recovery_present.data(),
            &transient_plan);
    if (transient_result == LEO2_UNSUPPORTED)
        return;
    require_result(transient_result, "diagnostic transient plan create");
    size_t transient_scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(
        transient_plan, bytes, &transient_scratch_bytes),
        "diagnostic transient scratch query");
    require(transient_scratch_bytes <= scratch.size(),
        std::string(test.name) +
            " diagnostic transient plan exceeds one-shot scratch");
    memset(restored_storage.data(), 0, restored_storage.size());
    require_result(
        leopard2_internal::ExecuteOneShotTransformPlanForDiagnostics(
            transient_plan, bytes, original_present.data(),
            recovery_present.data(), decode_original.data(), recovery.data(),
            restored.data(), scratch.data(), scratch.size()),
        "diagnostic transient plan execute");
    for (size_t i = 0; i < missing.size(); ++i)
    {
        require(memcmp(restored[missing[i]], expected_original[missing[i]],
                bytes) == 0,
            std::string(test.name) +
                " diagnostic transient restored original mismatch");
    }
    leo2_decode_plan_destroy(transient_plan);
}

void check_one_shot_short_scratch_atomicity(
    const PolicyCase& test,
    const leo2_codec* codec,
    const std::vector<uint8_t>& original_present,
    const std::vector<uint8_t>& recovery_present,
    const std::vector<const void*>& decode_original,
    const std::vector<const void*>& recovery,
    const std::vector<uint32_t>& missing,
    size_t bytes)
{
    require(leopard2_internal::SetOneShotPlanSetupModeForDiagnostics(3),
        "cannot select production one-shot setup mode");
    size_t scratch_bytes = 0;
    require_result(leo2_decode_scratch_size(codec, bytes, &scratch_bytes),
        "atomicity scratch query");
    leo2_decode_plan* transient_plan = NULL;
    require_result(
        leopard2_internal::CreateOneShotTransformPlanForDiagnostics(
            codec, bytes, original_present.data(), recovery_present.data(),
            &transient_plan),
        "atomicity transient plan create");
    size_t transient_scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(
        transient_plan, bytes, &transient_scratch_bytes),
        "atomicity transient scratch query");
    leo2_decode_plan_destroy(transient_plan);
    require(transient_scratch_bytes != 0 &&
            transient_scratch_bytes <= scratch_bytes,
        "transform atomicity case unexpectedly needs no scratch");
    AlignedBytes scratch(scratch_bytes);
    memset(scratch.data(), 0x3c, scratch.size());
    AlignedBytes restored_storage(
        static_cast<size_t>(missing.size()) * bytes);
    memset(restored_storage.data(), 0xa5, restored_storage.size());
    std::vector<void*> restored(test.k, NULL);
    for (size_t i = 0; i < missing.size(); ++i)
        restored[missing[i]] = restored_storage.data() + i * bytes;
    const std::vector<uint8_t> original_present_before = original_present;
    const std::vector<uint8_t> recovery_present_before = recovery_present;

    const leo2_result result = leo2_decode(codec, bytes,
        original_present.data(), recovery_present.data(),
        decode_original.data(), recovery.data(), restored.data(),
        scratch.data(), transient_scratch_bytes - 1);
    require(result == LEO2_SCRATCH_TOO_SMALL,
        std::string(test.name) + " short scratch returned " +
            leo2_result_string(result));
    for (size_t i = 0; i < restored_storage.size(); ++i)
    {
        require(restored_storage.data()[i] == 0xa5,
            std::string(test.name) +
                " modified output before rejecting short scratch");
    }
    require(original_present == original_present_before &&
            recovery_present == recovery_present_before,
        std::string(test.name) +
            " modified a presence map before rejecting short scratch");
    for (size_t i = 0; i < scratch.size(); ++i)
    {
        require(scratch.data()[i] == 0x3c,
            std::string(test.name) +
                " modified scratch before rejecting short scratch");
    }
}

void check_one_shot_presence_aliases(
    const PolicyCase& test,
    const leo2_codec* codec,
    const std::vector<uint8_t>& original_present,
    const std::vector<uint8_t>& recovery_present,
    const std::vector<const void*>& decode_original,
    const std::vector<const void*>& recovery,
    const std::vector<uint32_t>& missing,
    size_t bytes)
{
    require(leopard2_internal::SetOneShotPlanSetupModeForDiagnostics(3),
        "cannot select production one-shot setup mode");
    size_t scratch_bytes = 0;
    require_result(leo2_decode_scratch_size(codec, bytes, &scratch_bytes),
        "presence-alias scratch query");
    require(scratch_bytes != 0 && !missing.empty(),
        "presence-alias case is not a transform decode");

    const size_t original_map_bytes = std::max<size_t>(
        scratch_bytes, original_present.size());
    const size_t recovery_map_bytes = std::max<size_t>(
        scratch_bytes, recovery_present.size());
    AlignedBytes original_map(original_map_bytes);
    AlignedBytes recovery_map(recovery_map_bytes);
    memcpy(original_map.data(), original_present.data(),
        original_present.size());
    memcpy(recovery_map.data(), recovery_present.data(),
        recovery_present.size());
    const std::vector<uint8_t> original_before(
        original_map.data(), original_map.data() + original_present.size());
    const std::vector<uint8_t> recovery_before(
        recovery_map.data(), recovery_map.data() + recovery_present.size());

    AlignedBytes scratch(scratch_bytes);
    AlignedBytes restored_storage(
        static_cast<size_t>(missing.size()) * bytes);
    std::vector<void*> restored(test.k, NULL);
    for (size_t i = 0; i < missing.size(); ++i)
        restored[missing[i]] = restored_storage.data() + i * bytes;

    const auto require_maps_unchanged = [&]() {
        require(memcmp(original_map.data(), original_before.data(),
                original_before.size()) == 0 &&
                memcmp(recovery_map.data(), recovery_before.data(),
                    recovery_before.size()) == 0,
            std::string(test.name) +
                " modified a protected presence map");
    };
    const auto expect_overlap = [&](const uint8_t* original_map_pointer,
                                    const uint8_t* recovery_map_pointer,
                                    void* scratch_pointer,
                                    std::vector<void*>& outputs,
                                    const char* operation) {
        const leo2_result result = leo2_decode(codec, bytes,
            original_map_pointer, recovery_map_pointer,
            decode_original.data(), recovery.data(), outputs.data(),
            scratch_pointer, scratch_bytes);
        require(result == LEO2_OVERLAP,
            std::string(test.name) + " " + operation + " returned " +
                leo2_result_string(result));
        require_maps_unchanged();
    };

    expect_overlap(original_map.data(), recovery_map.data(),
        original_map.data(), restored, "original-presence/scratch alias");
    expect_overlap(original_map.data(), recovery_map.data(),
        recovery_map.data(), restored, "recovery-presence/scratch alias");

    restored[missing[0]] = original_map.data();
    expect_overlap(original_map.data(), recovery_map.data(),
        scratch.data(), restored, "original-presence/output alias");
    restored[missing[0]] = recovery_map.data();
    expect_overlap(original_map.data(), recovery_map.data(),
        scratch.data(), restored, "recovery-presence/output alias");
}

void check_one_shot_concurrent_execution(
    const PolicyCase& test,
    const leo2_codec* codec,
    const std::vector<uint8_t>& original_present,
    const std::vector<uint8_t>& recovery_present,
    const std::vector<const void*>& decode_original,
    const std::vector<const void*>& expected_original,
    const std::vector<const void*>& recovery,
    const std::vector<uint32_t>& missing,
    size_t bytes)
{
    require(leopard2_internal::SetOneShotPlanSetupModeForDiagnostics(3),
        "cannot select production one-shot setup mode");
    size_t scratch_bytes = 0;
    require_result(leo2_decode_scratch_size(codec, bytes, &scratch_bytes),
        "concurrent one-shot scratch query");
    const unsigned worker_count = 4;
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
        if (leo2_decode(codec, bytes, original_present.data(),
                recovery_present.data(), decode_original.data(),
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
        std::string(test.name) + " concurrent one-shot decode failed");
    for (unsigned worker = 0; worker < worker_count; ++worker)
    {
        for (size_t i = 0; i < missing.size(); ++i)
        {
            require(memcmp(restored[worker][missing[i]],
                    expected_original[missing[i]], bytes) == 0,
                std::string(test.name) +
                    " concurrent one-shot output mismatch");
        }
        delete restored_storage[worker];
        delete scratch[worker];
    }
}

void check_transient_plan_rejections()
{
    leo2_context* context = create_context(LEO2_BACKEND_AVX2);
    if (!context)
        return;
    require(leopard2_internal::SetOneShotPlanSetupModeForDiagnostics(3),
        "cannot select transient rejection mode");

    const auto check = [&](const char* name, leo2_profile profile,
                           uint32_t k, uint32_t r,
                           uint32_t missing_count) {
        leo2_codec_options options;
        memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        leo2_codec* codec = NULL;
        require_result(leo2_codec_create(context, k, r, profile,
            LEO2_FIELD_GF8, &options, &codec),
            std::string(name) + " codec create");
        std::vector<uint8_t> original_present(k, 1);
        std::vector<uint8_t> recovery_present(r, 1);
        for (uint32_t i = 0; i < missing_count; ++i)
            original_present[i] = 0;
        leo2_decode_plan* plan = NULL;
        const leo2_result result =
            leopard2_internal::CreateOneShotTransformPlanForDiagnostics(
                codec, 64, original_present.data(),
                recovery_present.data(), &plan);
        const bool returned_plan = plan != NULL;
        leo2_decode_plan_destroy(plan);
        leo2_codec_destroy(codec);
        require(result == LEO2_UNSUPPORTED && !returned_plan,
            std::string(name) + " transient transform gate returned " +
                leo2_result_string(result));
    };

    check("no-loss", LEO2_PROFILE_LEGACY_HIGH_V1, 16, 8, 0);
    check("direct-xor", LEO2_PROFILE_LEGACY_HIGH_V1, 8, 1, 1);
    check("direct-copy", LEO2_PROFILE_LOW_V1, 1, 3, 1);
    check("direct-repair", LEO2_PROFILE_LEGACY_HIGH_V1, 4, 3, 2);
    leo2_context_destroy(context);
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
    const bool force_generic =
        strcmp(test.name, "forced-generic-256") == 0;
    if (strcmp(test.name, "negative-t4") == 0)
        codec_options.flags = LEO2_CODEC_FORCE_SPECIALIZED_DECODE;
    else if (force_generic)
        codec_options.flags = LEO2_CODEC_FORCE_GENERIC_DECODE;
    leo2_codec* codec = NULL;
    require_result(leo2_codec_create(context, test.k, test.r,
        test.profile, test.field, &codec_options, &codec), "codec create");

    size_t bytes = test.field == LEO2_FIELD_GF8 ? 65U : 66U;
    if (strcmp(test.name, "translated-low-256") == 0 || force_generic)
        bytes = 256;
    else if (strcmp(test.name, "target-k106r53-1024") == 0)
        bytes = 1024;
    else if (strcmp(test.name, "target-k106r53-1025") == 0)
        bytes = 1025;
    else if (strncmp(test.name, "tiny-direct-", 12) == 0)
        bytes = 63;
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
    if (strcmp(test.name, "mixed-original-parity-loss") == 0)
    {
        recovery_present[0] = 0;
        decode_recovery[0] = NULL;
    }

    leo2_decode_plan* plan = NULL;
    require_result(leo2_decode_plan_create(codec, original_present.data(),
        recovery_present.data(), &plan), "decode plan create");
    if (strncmp(test.name, "tiny-direct-", 12) == 0)
    {
        leopard2_internal::DecodePathInfo path;
        require_result(leopard2_internal::GetDecodePlanPathInfo(
            plan, bytes, false, &path),
            std::string(test.name) + " path introspection");
        const bool control =
            strcmp(test.name, "tiny-direct-k8r8l8-control") == 0;
        require((path.path == leopard2_internal::kDecodePathDirect) ==
                !control,
            std::string(test.name) + " selected the wrong tiny path");
        require(path.direct_executor == (control
                ? leopard2_internal::kDirectRepairExecutorNone
                : leopard2_internal::kDirectRepairExecutorSourceMajor),
            std::string(test.name) + " selected the wrong tiny executor");
    }
    leopard2_internal::DecodePlanPrunedScheduleInfo info;
    require(leopard2_internal::GetDecodePlanPrunedScheduleInfo(plan, &info),
        std::string(test.name) + " schedule introspection failed");
    const size_t high_count =
        info.high_input_plan_count + info.high_output_plan_count;
    const size_t low_count =
        info.low_input_plan_count + info.low_output_plan_count;
    const bool tiny_direct_case =
        strncmp(test.name, "tiny-direct-", 12) == 0;
    if (force_generic)
    {
        require(high_count == 0 && low_count == 0,
            std::string(test.name) +
                " forced generic plan compiled specialized schedules");
    }
    else if (tiny_direct_case)
    {
        // The assertions above fix the byte-aware execution route.  Depending
        // on shape, its retained transform fallback is either a full schedule
        // (reported as zero pruned entries) or a translated/pruned schedule.
    }
    else if (expect_skip)
    {
        require(high_count == 0 && low_count == 0,
            std::string(test.name) + " retained a pruned schedule");
    }
    else if (strcmp(test.name, "translated-low-256") == 0)
    {
        // The public wire profile remains legacy-high, while this balanced
        // decode executes Algorithm 4 in the translated low-profile view.
        require(high_count == 0 && low_count != 0,
            std::string(test.name) +
                " did not retain only its translated-low schedule");
    }
    else if (test.profile == LEO2_PROFILE_LEGACY_HIGH_V1)
        require(high_count != 0,
            std::string(test.name) + " did not compile its control schedule");
    else
        require(low_count != 0,
            std::string(test.name) + " did not compile its low-profile schedule");

    execute_and_check(test, plan, decode_original, original, decode_recovery,
        missing, bytes, expect_skip);
    execute_one_shot_and_check(test, codec, original_present,
        recovery_present, decode_original, original, decode_recovery,
        missing, bytes, 0);
    execute_one_shot_and_check(test, codec, original_present,
        recovery_present, decode_original, original, decode_recovery,
        missing, bytes, 1);
    execute_one_shot_and_check(test, codec, original_present,
        recovery_present, decode_original, original, decode_recovery,
        missing, bytes, 2);
    execute_one_shot_and_check(test, codec, original_present,
        recovery_present, decode_original, original, decode_recovery,
        missing, bytes, 3);
    require(leopard2_internal::SetOneShotPlanSetupModeForDiagnostics(3),
        "cannot restore one-shot production setup mode");
    if (strcmp(test.name, "target-k106r53-1024") == 0)
    {
        check_one_shot_short_scratch_atomicity(test, codec,
            original_present, recovery_present, decode_original,
            decode_recovery, missing, bytes);
        check_one_shot_presence_aliases(test, codec,
            original_present, recovery_present, decode_original,
            decode_recovery, missing, bytes);
        check_one_shot_concurrent_execution(test, codec,
            original_present, recovery_present, decode_original, original,
            decode_recovery, missing, bytes);
    }
    // Diagnostic-mode changes are process-local only.  Re-execute the public
    // immutable plan after the complete sequence and after overwriting the
    // caller's presence maps.  This proves that it retained owned metadata
    // and did not observe transient setup policy or caller-map lifetime.
    std::fill(original_present.begin(), original_present.end(), 0);
    std::fill(recovery_present.begin(), recovery_present.end(), 0);
    execute_and_check(test, plan, decode_original, original, decode_recovery,
        missing, bytes, false);
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
        { "tiny-direct-k5r5l5", LEO2_BACKEND_AVX2, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 5, 5, 5, false },
        { "tiny-direct-k8r8l5", LEO2_BACKEND_AVX2, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 8, 8, 5, false },
        { "tiny-direct-k8r8l8-control", LEO2_BACKEND_AVX2,
          LEO2_FIELD_GF8, LEO2_PROFILE_LEGACY_HIGH_V1, 8, 8, 8, false },
        { "tiny-direct-k16r8l8", LEO2_BACKEND_AVX2, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 16, 8, 8, false },
        { "negative-t32", LEO2_BACKEND_AVX2, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 64, 32, 32, false },
        { "negative-t32-r-minus-one", LEO2_BACKEND_AVX2, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 64, 31, 31, false },
        { "small-dual-k-below-boundary", LEO2_BACKEND_AVX2,
          LEO2_FIELD_GF8, LEO2_PROFILE_LEGACY_HIGH_V1, 14, 8, 8,
          kExpectSmallDualRegularFallback },
        { "mixed-original-parity-loss", LEO2_BACKEND_AVX2, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 16, 8, 7, false },
        { "small-dual-r-two-below-t", LEO2_BACKEND_AVX2, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 16, 6, 6,
          kExpectSmallDualRegularFallback },
        { "positive-auto-avx2", LEO2_BACKEND_AUTO, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 16, 8, 8, true },
        // GFNI executes the same AVX2-tier dense schedule and the setup
        // policy names it explicitly.  Unqualified hosts skip the row.
        { "positive-gfni-t8-k2r", LEO2_BACKEND_GFNI, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 16, 8, 8, true },
        { "negative-scalar", LEO2_BACKEND_SCALAR, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 16, 8, 8, false },
        { "negative-gf16", LEO2_BACKEND_AVX2, LEO2_FIELD_GF16,
          LEO2_PROFILE_LEGACY_HIGH_V1, 16, 8, 8, false },
        { "negative-low-profile", LEO2_BACKEND_AVX2, LEO2_FIELD_GF8,
          LEO2_PROFILE_LOW_V1, 16, 8, 8, false },
        { "translated-low-256", LEO2_BACKEND_AVX2, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 15, 15, 15, false },
        { "forced-generic-256", LEO2_BACKEND_AVX2, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 96, 32, 32, false },
        { "target-k106r53-1024", LEO2_BACKEND_AVX2, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 106, 53, 53, false },
        { "target-k106r53-1025", LEO2_BACKEND_AVX2, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 106, 53, 53, false },
        { "negative-ssse3", LEO2_BACKEND_SSSE3, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 16, 8, 8, false },
        { "negative-avx512", LEO2_BACKEND_AVX512, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 16, 8, 8, false },
        { "negative-neon", LEO2_BACKEND_NEON, LEO2_FIELD_GF8,
          LEO2_PROFILE_LEGACY_HIGH_V1, 16, 8, 8, false }
    };

    unsigned executed = 0;
    unsigned skipped = 0;
    try
    {
        require(!leopard2_internal::
                SetOneShotPlanSetupModeForDiagnostics(4),
            "one-shot setup accepted an unsupported diagnostic mode");
        require(leopard2_internal::
                SetOneShotPlanSetupModeForDiagnostics(3),
            "one-shot setup production mode unavailable");
        check_transient_plan_rejections();
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
