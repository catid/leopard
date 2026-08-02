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

#if defined(LEO2_TEST_ONE_SHOT_NO_LOSS_HOOK_BUILD)
#include "Leopard2Direct.h"
#endif

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <new>
#include <stdexcept>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

#ifndef LEO2_EXPERIMENT_ONE_SHOT_NO_LOSS_SHORT_CIRCUIT
#define LEO2_EXPERIMENT_ONE_SHOT_NO_LOSS_SHORT_CIRCUIT 1
#endif

#ifndef LEO2_TEST_EXPECT_GF8
#define LEO2_TEST_EXPECT_GF8 0
#endif

#ifndef LEO2_TEST_EXPECT_GF16
#define LEO2_TEST_EXPECT_GF16 0
#endif

#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
static std::atomic<bool> g_track_allocations(false);
static std::atomic<uint64_t> g_tracked_allocations(0);

#if defined(_MSC_VER)
#define LEO2_TEST_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define LEO2_TEST_NOINLINE __attribute__((noinline))
#else
#define LEO2_TEST_NOINLINE
#endif

LEO2_TEST_NOINLINE void* operator new(size_t bytes)
{
    if (g_track_allocations.load(std::memory_order_relaxed))
        g_tracked_allocations.fetch_add(1, std::memory_order_relaxed);
    void* const result = malloc(bytes == 0 ? 1 : bytes);
    if (!result)
        throw std::bad_alloc();
    return result;
}

LEO2_TEST_NOINLINE void* operator new[](size_t bytes)
{
    return ::operator new(bytes);
}

LEO2_TEST_NOINLINE void* operator new(
    size_t bytes, const std::nothrow_t&) noexcept
{
    try
    {
        return ::operator new(bytes);
    }
    catch (...)
    {
        return NULL;
    }
}

LEO2_TEST_NOINLINE void* operator new[](
    size_t bytes, const std::nothrow_t&) noexcept
{
    return ::operator new(bytes, std::nothrow);
}

LEO2_TEST_NOINLINE void operator delete(void* pointer) noexcept
{
    free(pointer);
}

LEO2_TEST_NOINLINE void operator delete[](void* pointer) noexcept
{
    free(pointer);
}

LEO2_TEST_NOINLINE void operator delete(
    void* pointer, const std::nothrow_t&) noexcept
{
    ::operator delete(pointer);
}

LEO2_TEST_NOINLINE void operator delete[](
    void* pointer, const std::nothrow_t&) noexcept
{
    ::operator delete[](pointer);
}

#undef LEO2_TEST_NOINLINE
#endif

namespace {

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

struct DecodeObservation
{
    leo2_result result;
    uint64_t allocations;
};

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void require_result(
    leo2_result actual,
    leo2_result expected,
    const char* operation)
{
    if (actual != expected)
    {
        std::cerr << operation << ": expected " << expected
                  << ", got " << actual << std::endl;
        throw std::runtime_error(operation);
    }
}

void begin_allocation_audit()
{
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
    g_tracked_allocations.store(0, std::memory_order_relaxed);
    g_track_allocations.store(true, std::memory_order_release);
#endif
}

uint64_t end_allocation_audit()
{
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
    g_track_allocations.store(false, std::memory_order_release);
    return g_tracked_allocations.load(std::memory_order_relaxed);
#else
    return 0;
#endif
}

DecodeObservation observe_decode(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    const uint8_t* original_present,
    const uint8_t* recovery_present,
    const void* const* original,
    const void* const* recovery,
    void* const* restored,
    void* scratch,
    size_t scratch_bytes)
{
    begin_allocation_audit();
    const leo2_result result = leo2_decode(codec, shard_bytes,
        original_present, recovery_present, original, recovery, restored,
        scratch, scratch_bytes);
    const uint64_t allocations = end_allocation_audit();
    const DecodeObservation observation = { result, allocations };
    return observation;
}

void require_no_loss_observation(
    const DecodeObservation& observation,
    const char* operation)
{
    require_result(observation.result, LEO2_SUCCESS, operation);
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
# if LEO2_EXPERIMENT_ONE_SHOT_NO_LOSS_SHORT_CIRCUIT
    require(observation.allocations == 0,
        "enabled one-shot no-loss path allocated");
# else
    require(observation.allocations > 0,
        "disabled one-shot no-loss control avoided plan allocation");
# endif
#else
    (void)observation;
#endif
}

leo2_codec* create_codec(
    leo2_context* context,
    uint32_t k,
    uint32_t r,
    leo2_profile profile,
    leo2_field field)
{
    leo2_codec* codec = NULL;
    require_result(leo2_codec_create(
        context, k, r, profile, field, NULL, &codec),
        LEO2_SUCCESS, "codec create");
    return codec;
}

void test_no_loss_shape(
    leo2_context* context,
    uint32_t k,
    uint32_t r,
    leo2_profile profile,
    leo2_field field,
    uint64_t odd_or_zero_bytes)
{
    leo2_codec* const codec = create_codec(
        context, k, r, profile, field);
    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 0);

    DecodeObservation observation = observe_decode(codec, 0,
        original_present.data(), recovery_present.data(), NULL, NULL, NULL,
        NULL, 0);
    require_no_loss_observation(
        observation, "zero-byte missing-parity no-loss decode");
    if (odd_or_zero_bytes != 0)
    {
        observation = observe_decode(codec, odd_or_zero_bytes,
            original_present.data(), recovery_present.data(), NULL, NULL,
            NULL, NULL, 0);
        require_no_loss_observation(
            observation, "odd-byte missing-parity no-loss decode");
    }

    std::fill(recovery_present.begin(), recovery_present.end(), 1);
    observation = observe_decode(codec,
        std::numeric_limits<uint64_t>::max(), original_present.data(),
        recovery_present.data(), NULL, NULL, NULL, NULL, 0);
    require_no_loss_observation(
        observation, "surplus-parity no-loss decode");

    for (uint32_t i = 0; i < r; ++i)
        recovery_present[i] = static_cast<uint8_t>(i & 1u);
    observation = observe_decode(codec, odd_or_zero_bytes,
        original_present.data(), recovery_present.data(), NULL, NULL, NULL,
        NULL, 0);
    require_no_loss_observation(
        observation, "mixed-parity no-loss decode");

    /*
        No-loss execution is a true no-op after presence validation.  Supply
        data that would fail ordinary alias and scratch checks and prove it is
        neither inspected nor modified.
    */
    uint8_t payload = 0x5a;
    std::vector<const void*> original(k, &payload);
    std::vector<const void*> recovery(r, &payload);
    std::vector<void*> restored(k, original_present.data());
    const std::vector<uint8_t> original_before = original_present;
    const std::vector<uint8_t> recovery_before = recovery_present;
    observation = observe_decode(codec,
        std::numeric_limits<uint64_t>::max(), original_present.data(),
        recovery_present.data(), original.data(), recovery.data(),
        restored.data(), original_present.data(),
        std::numeric_limits<size_t>::max());
    require_no_loss_observation(
        observation, "alias-heavy no-loss decode");
    require(original_present == original_before &&
            recovery_present == recovery_before && payload == 0x5a,
        "one-shot no-loss decode modified caller storage");

    leo2_codec_destroy(codec);
}

void require_preplan_failure(
    const DecodeObservation& observation,
    leo2_result expected,
    const char* operation)
{
    require_result(observation.result, expected, operation);
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
    require(observation.allocations == 0,
        "presence validation allocated before rejecting the call");
#else
    (void)observation;
#endif
}

void test_validation_and_failure_atomicity(
    leo2_context* context,
    leo2_field field)
{
    const uint32_t k = 9;
    const uint32_t r = 7;
    leo2_codec* const codec = create_codec(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, field);
    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 0);
    uint8_t output[17];
    memset(output, 0xa5, sizeof(output));
    void* restored[9] = { output, output, output, output, output,
                          output, output, output, output };

    original_present.back() = 2;
    require_preplan_failure(observe_decode(codec, sizeof(output),
        original_present.data(), recovery_present.data(), NULL, NULL,
        restored, output, sizeof(output)), LEO2_INVALID_ARGUMENT,
        "late invalid original presence");
    original_present.back() = 1;

    recovery_present.back() = 2;
    require_preplan_failure(observe_decode(codec, sizeof(output),
        original_present.data(), recovery_present.data(), NULL, NULL,
        restored, output, sizeof(output)), LEO2_INVALID_ARGUMENT,
        "late invalid recovery presence");
    recovery_present.back() = 0;

    /*
        The speculative scan stops at the first loss.  The canonical fallback
        must still validate every later presence byte.
    */
    original_present.front() = 0;
    original_present.back() = 2;
    require_preplan_failure(observe_decode(codec, sizeof(output),
        original_present.data(), recovery_present.data(), NULL, NULL,
        restored, output, sizeof(output)), LEO2_INVALID_ARGUMENT,
        "loss-prefix late invalid original presence");
    original_present.back() = 1;
    recovery_present.back() = 2;
    require_preplan_failure(observe_decode(codec, sizeof(output),
        original_present.data(), recovery_present.data(), NULL, NULL,
        restored, output, sizeof(output)), LEO2_INVALID_ARGUMENT,
        "loss-prefix late invalid recovery presence");
    recovery_present.back() = 0;

    require_preplan_failure(observe_decode(codec, sizeof(output),
        original_present.data(), recovery_present.data(), NULL, NULL,
        restored, output, sizeof(output)), LEO2_NEED_MORE_DATA,
        "loss-prefix insufficient data");
    original_present.front() = 1;

    const uintptr_t address_end =
        std::numeric_limits<uintptr_t>::max() - 1;
    const uint8_t* const overflowing =
        reinterpret_cast<const uint8_t*>(address_end);
    require_preplan_failure(observe_decode(codec, sizeof(output), overflowing,
        recovery_present.data(), NULL, NULL, restored, output,
        sizeof(output)), LEO2_INVALID_ARGUMENT,
        "overflowing original-presence span");
    require_preplan_failure(observe_decode(codec, sizeof(output),
        original_present.data(), overflowing, NULL, NULL, restored, output,
        sizeof(output)), LEO2_INVALID_ARGUMENT,
        "overflowing recovery-presence span");
    require_preplan_failure(observe_decode(NULL, sizeof(output),
        original_present.data(), recovery_present.data(), NULL, NULL,
        restored, output, sizeof(output)), LEO2_INVALID_ARGUMENT,
        "null codec");
    require_preplan_failure(observe_decode(codec, sizeof(output), NULL,
        recovery_present.data(), NULL, NULL, restored, output,
        sizeof(output)), LEO2_INVALID_ARGUMENT,
        "null original presence");
    require_preplan_failure(observe_decode(codec, sizeof(output),
        original_present.data(), NULL, NULL, NULL, restored, output,
        sizeof(output)), LEO2_INVALID_ARGUMENT,
        "null recovery presence");

    for (size_t i = 0; i < sizeof(output); ++i)
        require(output[i] == 0xa5,
            "rejected one-shot decode modified output or scratch");
    leo2_codec_destroy(codec);
}

void test_r1_presence_validation(
    leo2_context* context,
    leo2_field field)
{
    const uint32_t k = 9;
    leo2_codec* const codec = create_codec(context, k, 1,
        LEO2_PROFILE_LEGACY_HIGH_V1, field);
    std::vector<uint8_t> original_present(k, 1);
    uint8_t recovery_present[1] = { 0 };

    recovery_present[0] = 2;
    require_preplan_failure(observe_decode(codec, 0,
        original_present.data(), recovery_present, NULL, NULL, NULL, NULL, 0),
        LEO2_INVALID_ARGUMENT, "R=1 invalid recovery presence after no loss");

    recovery_present[0] = 0;
    original_present.back() = 2;
    require_preplan_failure(observe_decode(codec, 0,
        original_present.data(), recovery_present, NULL, NULL, NULL, NULL, 0),
        LEO2_INVALID_ARGUMENT, "R=1 late invalid original presence");

    recovery_present[0] = 1;
    original_present.front() = 0;
    require_preplan_failure(observe_decode(codec, 0,
        original_present.data(), recovery_present, NULL, NULL, NULL, NULL, 0),
        LEO2_INVALID_ARGUMENT,
        "R=1 invalid original presence after an earlier loss");

    original_present.back() = 1;
    recovery_present[0] = 2;
    require_preplan_failure(observe_decode(codec, 0,
        original_present.data(), recovery_present, NULL, NULL, NULL, NULL, 0),
        LEO2_INVALID_ARGUMENT,
        "R=1 invalid recovery presence after an earlier loss");

    leo2_codec_destroy(codec);
}

void test_loss_fallback(
    leo2_context* context,
    leo2_field field)
{
    const uint32_t k = 4;
    const uint32_t r = 3;
    const size_t bytes = 64;
    leo2_codec* const codec = create_codec(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, field);
    std::vector<std::vector<uint8_t> > originals(
        k, std::vector<uint8_t>(bytes));
    std::vector<std::vector<uint8_t> > parity(
        r, std::vector<uint8_t>(bytes));
    for (uint32_t shard = 0; shard < k; ++shard)
        for (size_t i = 0; i < bytes; ++i)
            originals[shard][i] = static_cast<uint8_t>(
                shard * 53u + i * 29u + 7u);
    const void* original_input[k];
    void* parity_output[r];
    for (uint32_t i = 0; i < k; ++i)
        original_input[i] = originals[i].data();
    for (uint32_t i = 0; i < r; ++i)
        parity_output[i] = parity[i].data();
    size_t encode_scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(
        codec, bytes, &encode_scratch_bytes),
        LEO2_SUCCESS, "encode scratch query");
    AlignedBuffer encode_scratch(encode_scratch_bytes);
    require_result(leo2_encode(codec, bytes, original_input, parity_output,
        encode_scratch.data(), encode_scratch.size()),
        LEO2_SUCCESS, "loss fallback fixture encode");

    size_t decode_scratch_bytes = 0;
    require_result(leo2_decode_scratch_size(
        codec, bytes, &decode_scratch_bytes),
        LEO2_SUCCESS, "decode scratch query");
    AlignedBuffer decode_scratch(decode_scratch_bytes);
    std::vector<uint8_t> restored_bytes(bytes, 0xa5);
    for (uint32_t missing_case = 0; missing_case < 2; ++missing_case)
    {
        const uint32_t missing = missing_case == 0 ? 0 : k - 1;
        uint8_t original_present[k];
        uint8_t recovery_present[r];
        const void* decode_original[k];
        const void* decode_recovery[r];
        void* restored[k];
        for (uint32_t i = 0; i < k; ++i)
        {
            original_present[i] = i == missing ? 0 : 1;
            decode_original[i] =
                i == missing ? NULL : originals[i].data();
            restored[i] =
                i == missing ? restored_bytes.data() : NULL;
        }
        for (uint32_t i = 0; i < r; ++i)
        {
            recovery_present[i] = 1;
            decode_recovery[i] = parity[i].data();
        }
        std::fill(restored_bytes.begin(), restored_bytes.end(), 0xa5);
        const DecodeObservation observation = observe_decode(
            codec, bytes, original_present, recovery_present,
            decode_original, decode_recovery, restored,
            decode_scratch.data(), decode_scratch.size());
        require_result(observation.result, LEO2_SUCCESS,
            "one-loss canonical fallback");
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
        require(observation.allocations > 0,
            "one-loss call bypassed canonical plan construction");
#endif
        require(restored_bytes == originals[missing],
            "one-loss canonical fallback restored wrong bytes");
    }
    leo2_codec_destroy(codec);
}

#if defined(LEO2_TEST_ONE_SHOT_NO_LOSS_HOOK_BUILD)
void test_forced_hook_modes(
    leo2_context* context,
    leo2_field field)
{
    leo2_codec* const codec = create_codec(context, 8, 8,
        LEO2_PROFILE_LEGACY_HIGH_V1, field);
    uint8_t original_present[8];
    uint8_t recovery_present[8];
    memset(original_present, 1, sizeof(original_present));
    memset(recovery_present, 0, sizeof(recovery_present));
    const leo2_test_decode_mode modes[] = {
        LEO2_TEST_DECODE_FORCE_TRANSLATED_LOW,
        LEO2_TEST_DECODE_FORCE_NATIVE_HIGH
    };
    for (size_t i = 0; i < sizeof(modes) / sizeof(modes[0]); ++i)
    {
        require_result(leo2_test_codec_set_decode_mode(codec, modes[i]),
            LEO2_SUCCESS, "select forced decode mode");
        const DecodeObservation observation = observe_decode(codec, 17,
            original_present, recovery_present, NULL, NULL, NULL, NULL, 0);
        require_result(observation.result, LEO2_SUCCESS,
            "forced-mode one-shot no-loss decode");
#if LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
        require(observation.allocations > 0,
            "forced hook mode bypassed canonical plan construction");
#endif
    }
    leo2_codec_destroy(codec);
}
#endif

} // namespace

int main()
{
    try
    {
        leo2_context_options options;
        memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AUTO;
        options.thread_count = 1;
        leo2_context* context = NULL;
        require_result(leo2_context_create(&options, &context),
            LEO2_SUCCESS, "context create");

#if LEO2_TEST_EXPECT_GF8
        test_no_loss_shape(context, 9, 7,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 0);
        test_no_loss_shape(context, 7, 9,
            LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8, 17);
        test_no_loss_shape(context, 65, 1,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 17);
        test_validation_and_failure_atomicity(
            context, LEO2_FIELD_GF8);
        test_r1_presence_validation(context, LEO2_FIELD_GF8);
        test_loss_fallback(context, LEO2_FIELD_GF8);
#endif
#if LEO2_TEST_EXPECT_GF16
        test_no_loss_shape(context, 9, 7,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, 65);
        test_no_loss_shape(context, 7, 9,
            LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16, 65);
        test_validation_and_failure_atomicity(
            context, LEO2_FIELD_GF16);
        test_r1_presence_validation(context, LEO2_FIELD_GF16);
        test_loss_fallback(context, LEO2_FIELD_GF16);
#endif

#if defined(LEO2_TEST_ONE_SHOT_NO_LOSS_HOOK_BUILD)
# if LEO2_TEST_EXPECT_GF8
        test_forced_hook_modes(context, LEO2_FIELD_GF8);
# elif LEO2_TEST_EXPECT_GF16
        test_forced_hook_modes(context, LEO2_FIELD_GF16);
# endif
#endif
        leo2_context_destroy(context);

        std::cout << "{\"one_shot_no_loss_short_circuit\":"
                  << LEO2_EXPERIMENT_ONE_SHOT_NO_LOSS_SHORT_CIRCUIT
                  << ",\"hook_build\":"
#if defined(LEO2_TEST_ONE_SHOT_NO_LOSS_HOOK_BUILD)
                  << 1
#else
                  << 0
#endif
                  << ",\"allocation_audit\":"
                  << LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE
                  << ",\"status\":\"ok\"}" << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << std::endl;
        return 1;
    }
}
