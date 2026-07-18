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

#include "leopard.h"
#include "leopard2.h"
#include "LeopardCommon.h"

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#if defined(__linux__)
#include <sched.h>
#endif

#if defined(_MSC_VER)
#include <malloc.h>
#endif

extern "C" {
uint32_t leo2_test_context_worker_count(const leo2_context* context);
void leo2_test_set_thread_start_fault(int enabled);
unsigned leo2_test_thread_start_fault_consumptions(void);
}

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
    const void* data() const { return data_; }
    size_t size() const { return bytes_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);

    void* data_;
    size_t bytes_;
};

class ContextOwner
{
public:
    ContextOwner() : context(NULL) {}
    ~ContextOwner() { leo2_context_destroy(context); }
    leo2_context* context;

private:
    ContextOwner(const ContextOwner&);
    ContextOwner& operator=(const ContextOwner&);
};

class CodecOwner
{
public:
    CodecOwner() : codec(NULL) {}
    ~CodecOwner() { leo2_codec_destroy(codec); }
    leo2_codec* codec;

private:
    CodecOwner(const CodecOwner&);
    CodecOwner& operator=(const CodecOwner&);
};

class PlanOwner
{
public:
    PlanOwner() : plan(NULL) {}
    ~PlanOwner() { leo2_decode_plan_destroy(plan); }
    leo2_decode_plan* plan;

private:
    PlanOwner(const PlanOwner&);
    PlanOwner& operator=(const PlanOwner&);
};

typedef std::vector<uint8_t> Bytes;
typedef std::vector<Bytes> Shards;

struct Counts
{
    Counts()
        : result_strings(0)
        , introspection_checks(0)
        , alias_checks(0)
        , scratch_checks(0)
        , scheduler_checks(0)
        , concurrent_batches(0)
        , batch_failure_checks(0)
        , legacy_checks(0)
    {}

    uint64_t result_strings;
    uint64_t introspection_checks;
    uint64_t alias_checks;
    uint64_t scratch_checks;
    uint64_t scheduler_checks;
    uint64_t concurrent_batches;
    uint64_t batch_failure_checks;
    uint64_t legacy_checks;
};

void require(bool condition, const std::string& operation)
{
    if (!condition)
        throw std::runtime_error(operation);
}

void require_result(
    leo2_result actual,
    leo2_result expected,
    const std::string& operation)
{
    if (actual == expected)
        return;
    throw std::runtime_error(operation + ": got " + leo2_result_string(actual) +
        " (" + std::to_string(static_cast<int>(actual)) + "), expected " +
        leo2_result_string(expected) + " (" +
        std::to_string(static_cast<int>(expected)) + ")");
}

std::vector<const void*> const_pointers(const Shards& shards)
{
    std::vector<const void*> pointers(shards.size(), NULL);
    for (size_t i = 0; i < shards.size(); ++i)
        pointers[i] = &shards[i][0];
    return pointers;
}

std::vector<void*> mutable_pointers(Shards& shards)
{
    std::vector<void*> pointers(shards.size(), NULL);
    for (size_t i = 0; i < shards.size(); ++i)
        pointers[i] = &shards[i][0];
    return pointers;
}

void fill_shards(Shards& shards, uint32_t seed)
{
    uint32_t state = seed;
    for (size_t shard = 0; shard < shards.size(); ++shard)
        for (size_t byte_i = 0; byte_i < shards[shard].size(); ++byte_i)
        {
            state = state * 1664525u + 1013904223u;
            shards[shard][byte_i] = static_cast<uint8_t>(
                (state >> 24) ^ (shard * 29u + byte_i * 7u));
        }
}

leo2_backend expected_backend_from_environment()
{
    const char* name = std::getenv("LEO2_EXPECT_BACKEND");
    if (!name || !name[0])
        return LEO2_BACKEND_AUTO;
    if (strcmp(name, "scalar") == 0)
        return LEO2_BACKEND_SCALAR;
    if (strcmp(name, "ssse3") == 0)
        return LEO2_BACKEND_SSSE3;
    if (strcmp(name, "avx2") == 0)
        return LEO2_BACKEND_AVX2;
    if (strcmp(name, "neon") == 0)
        return LEO2_BACKEND_NEON;
    if (strcmp(name, "avx512") == 0 || strcmp(name, "avx512vl") == 0)
        return LEO2_BACKEND_AVX512;
    throw std::runtime_error("unknown LEO2_EXPECT_BACKEND value");
}

void test_result_strings(Counts* counts)
{
    struct Expected
    {
        leo2_result result;
        const char* text;
    };
    const Expected expected[] = {
        { LEO2_SUCCESS, "Operation succeeded" },
        { LEO2_NEED_MORE_DATA, "Not enough received shards" },
        { LEO2_INVALID_ARGUMENT, "Invalid argument" },
        { LEO2_INVALID_COUNTS, "Invalid or unsupported shard counts" },
        { LEO2_UNSUPPORTED, "Requested profile, field, or backend is unsupported" },
        { LEO2_SCRATCH_TOO_SMALL, "Scratch buffer is too small" },
        { LEO2_BAD_ALIGNMENT, "Scratch buffer has insufficient alignment" },
        { LEO2_OVERLAP, "Unsupported input, output, or scratch overlap" },
        { LEO2_OUT_OF_MEMORY, "Allocation failed during setup" },
        { LEO2_INTERNAL_ERROR, "Internal initialization or execution error" }
    };
    for (size_t i = 0; i < sizeof(expected) / sizeof(expected[0]); ++i)
    {
        const char* actual = leo2_result_string(expected[i].result);
        require(actual != NULL && strcmp(actual, expected[i].text) == 0,
            "result string mismatch for value " +
            std::to_string(static_cast<int>(expected[i].result)));
        ++counts->result_strings;
    }
    require(strcmp(leo2_result_string(static_cast<leo2_result>(1)),
        "Unknown Leopard2 result") == 0, "positive unknown result string mismatch");
    require(strcmp(leo2_result_string(static_cast<leo2_result>(-100)),
        "Unknown Leopard2 result") == 0, "negative unknown result string mismatch");
    counts->result_strings += 2;
}

void verify_codec_identity(
    const leo2_codec* codec,
    uint32_t k,
    uint32_t r,
    uint32_t n,
    uint32_t padded,
    leo2_profile profile,
    leo2_field field,
    leo2_shard_layout layout,
    Counts* counts)
{
    require(leo2_codec_original_count(codec) == k, "original-count introspection");
    require(leo2_codec_recovery_count(codec) == r, "recovery-count introspection");
    require(leo2_codec_parent_count(codec) == n, "parent-count introspection");
    require(leo2_codec_padded_side(codec) == padded, "padded-side introspection");
    require(leo2_codec_profile(codec) == profile, "profile introspection");
    require(leo2_codec_field(codec) == field, "field introspection");
    require(leo2_codec_shard_layout(codec) == layout, "layout introspection");
    counts->introspection_checks += 7;
}

void test_introspection_and_null_contracts(
    leo2_context* context,
    Counts* counts)
{
    const size_t alignment = leo2_scratch_alignment();
    require(alignment != 0 && (alignment & (alignment - 1)) == 0,
        "scratch alignment is not a nonzero power of two");
    ++counts->introspection_checks;

    require(leo2_context_backend(NULL) == LEO2_BACKEND_AUTO,
        "null context backend sentinel");
    require(leo2_context_thread_count(NULL) == 0,
        "null context thread-count sentinel");
    require(leo2_codec_original_count(NULL) == 0, "null codec K sentinel");
    require(leo2_codec_recovery_count(NULL) == 0, "null codec R sentinel");
    require(leo2_codec_parent_count(NULL) == 0, "null codec N sentinel");
    require(leo2_codec_padded_side(NULL) == 0, "null codec padded sentinel");
    require(leo2_codec_profile(NULL) == LEO2_PROFILE_AUTO,
        "null codec profile sentinel");
    require(leo2_codec_field(NULL) == LEO2_FIELD_AUTO,
        "null codec field sentinel");
    require(leo2_codec_shard_layout(NULL) == LEO2_SHARD_LAYOUT_NATIVE_V1,
        "null codec layout sentinel");
    require(leo2_decode_plan_missing_original_count(NULL) == 0,
        "null plan missing-count sentinel");
    counts->introspection_checks += 10;

    const leo2_backend live_backend = leo2_context_backend(context);
    require(live_backend != LEO2_BACKEND_AUTO,
        "live context reports AUTO backend");
    require(leo2_context_thread_count(context) == 4,
        "explicit context thread count was not retained");
    require(live_backend >= LEO2_BACKEND_SCALAR &&
            live_backend <= LEO2_BACKEND_AVX512,
        "live context reports an out-of-range backend");
    const leo2_backend expected_backend = expected_backend_from_environment();
    if (expected_backend != LEO2_BACKEND_AUTO)
        require(live_backend == expected_backend,
            "forced build backend introspection mismatch");
    counts->introspection_checks += 3;

    leo2_context_options context_options;
    memset(&context_options, 0, sizeof(context_options));
    context_options.struct_size = sizeof(context_options);
    context_options.backend = leo2_context_backend(context);
    context_options.thread_count = 1;
    ContextOwner explicitly_selected;
    require_result(leo2_context_create(&context_options,
        &explicitly_selected.context), LEO2_SUCCESS,
        "explicit actual-backend context create");
    require(leo2_context_backend(explicitly_selected.context) ==
            leo2_context_backend(context) &&
        leo2_context_thread_count(explicitly_selected.context) == 1,
        "explicit backend/thread context introspection mismatch");
    counts->introspection_checks += 2;

    struct LargerContextOptions
    {
        leo2_context_options base;
        uint64_t ignored_tail;
    };
    LargerContextOptions larger_context;
    memset(&larger_context, 0, sizeof(larger_context));
    larger_context.base.struct_size = sizeof(larger_context);
    larger_context.base.backend = leo2_context_backend(context);
    larger_context.base.thread_count = 1;
    larger_context.ignored_tail = UINT64_C(0xfeedfacecafebeef);
    ContextOwner tail_context;
    require_result(leo2_context_create(&larger_context.base,
        &tail_context.context), LEO2_SUCCESS,
        "larger context options with ignored tail");
    require(leo2_context_thread_count(tail_context.context) == 1,
        "larger context options changed known fields");
    counts->introspection_checks += 2;

    leo2_context_options reserved_context = context_options;
    reserved_context.reserved = 1;
    leo2_context* rejected_context = reinterpret_cast<leo2_context*>(
        static_cast<uintptr_t>(1));
    const leo2_result reserved_context_result = leo2_context_create(
        &reserved_context, &rejected_context);
    if (reserved_context_result == LEO2_SUCCESS)
    {
        leo2_context_destroy(rejected_context);
        rejected_context = NULL;
    }
    require_result(reserved_context_result, LEO2_INVALID_ARGUMENT,
        "nonzero context reserved field");
    require(rejected_context == NULL,
        "reserved-context failure did not clear output");
    counts->introspection_checks += 2;

    require_result(leo2_context_create(&context_options, NULL),
        LEO2_INVALID_ARGUMENT, "null context output");
    leo2_codec* ignored_codec = reinterpret_cast<leo2_codec*>(
        static_cast<uintptr_t>(1));
    require_result(leo2_codec_create(NULL, 3, 2, LEO2_PROFILE_AUTO,
        LEO2_FIELD_AUTO, NULL, &ignored_codec), LEO2_INVALID_ARGUMENT,
        "null context codec create");
    require(ignored_codec == NULL,
        "null-context codec failure did not clear output");
    ignored_codec = reinterpret_cast<leo2_codec*>(static_cast<uintptr_t>(1));
    require_result(leo2_codec_create(context, 0, 2, LEO2_PROFILE_AUTO,
        LEO2_FIELD_AUTO, NULL, &ignored_codec), LEO2_INVALID_ARGUMENT,
        "zero-count codec create");
    require(ignored_codec == NULL,
        "zero-count codec failure did not clear output");
    require_result(leo2_codec_create(context, 3, 2, LEO2_PROFILE_AUTO,
        LEO2_FIELD_AUTO, NULL, NULL), LEO2_INVALID_ARGUMENT,
        "null codec output");
    leo2_context_destroy(NULL);
    leo2_codec_destroy(NULL);
    leo2_decode_plan_destroy(NULL);
    counts->introspection_checks += 10;

    CodecOwner high;
    require_result(leo2_codec_create(context, 9, 7, LEO2_PROFILE_AUTO,
        LEO2_FIELD_AUTO, NULL, &high.codec), LEO2_SUCCESS,
        "AUTO high codec create");
    verify_codec_identity(high.codec, 9, 7, 32, 8,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        LEO2_SHARD_LAYOUT_NATIVE_V1, counts);

    CodecOwner low;
    require_result(leo2_codec_create(context, 5, 11, LEO2_PROFILE_AUTO,
        LEO2_FIELD_AUTO, NULL, &low.codec), LEO2_SUCCESS,
        "AUTO low codec create");
    verify_codec_identity(low.codec, 5, 11, 32, 8,
        LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8,
        LEO2_SHARD_LAYOUT_NATIVE_V1, counts);

    leo2_codec_options options;
    memset(&options, 0, sizeof(options));
    options.struct_size = sizeof(options);
    options.shard_layout = LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1;
    CodecOwner padded;
    require_result(leo2_codec_create(context, 5, 3,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, &options,
        &padded.codec), LEO2_SUCCESS, "padded GF16 codec create");
    verify_codec_identity(padded.codec, 5, 3, 16, 4,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16,
        LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1, counts);

    struct LargerCodecOptions
    {
        leo2_codec_options base;
        uint64_t ignored_tail;
    };
    LargerCodecOptions larger_codec;
    memset(&larger_codec, 0, sizeof(larger_codec));
    larger_codec.base.struct_size = sizeof(larger_codec);
    larger_codec.ignored_tail = UINT64_C(0x0123456789abcdef);
    CodecOwner tail_codec;
    require_result(leo2_codec_create(context, 3, 2,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, &larger_codec.base,
        &tail_codec.codec), LEO2_SUCCESS,
        "larger codec options with ignored tail");
    verify_codec_identity(tail_codec.codec, 3, 2, 8, 2,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        LEO2_SHARD_LAYOUT_NATIVE_V1, counts);

    leo2_codec_options reserved_codec;
    memset(&reserved_codec, 0, sizeof(reserved_codec));
    reserved_codec.struct_size = sizeof(reserved_codec);
    reserved_codec.reserved = 1;
    CodecOwner rejected_codec;
    require_result(leo2_codec_create(context, 3, 2,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, &reserved_codec,
        &rejected_codec.codec), LEO2_INVALID_ARGUMENT,
        "nonzero codec reserved field");
    require(rejected_codec.codec == NULL,
        "reserved-codec failure did not clear output");
    counts->introspection_checks += 2;

    leo2_codec* experimental_codec = reinterpret_cast<leo2_codec*>(
        static_cast<uintptr_t>(1));
    require_result(leo2_codec_create(context, 3, 2,
        LEO2_PROFILE_EXACT_EXPERIMENTAL_V1, LEO2_FIELD_GF8, NULL,
        &experimental_codec), LEO2_UNSUPPORTED,
        "explicit unsupported exact profile");
    require(experimental_codec == NULL,
        "unsupported exact-profile failure did not clear output");
    counts->introspection_checks += 2;

    leo2_codec_options diagnostic_options;
    memset(&diagnostic_options, 0, sizeof(diagnostic_options));
    diagnostic_options.struct_size = sizeof(diagnostic_options);
    diagnostic_options.flags = LEO2_CODEC_FORCE_GENERIC_DECODE;
    CodecOwner diagnostic_codec;
    require_result(leo2_codec_create(context, 3, 2,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        &diagnostic_options, &diagnostic_codec.codec), LEO2_SUCCESS,
        "diagnostic decoder option");
    verify_codec_identity(diagnostic_codec.codec, 3, 2, 8, 2,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        LEO2_SHARD_LAYOUT_NATIVE_V1, counts);

    uint8_t valid_original_present[9];
    uint8_t valid_recovery_present[7];
    memset(valid_original_present, 1, sizeof(valid_original_present));
    memset(valid_recovery_present, 1, sizeof(valid_recovery_present));
    leo2_decode_plan* invalid_plan = reinterpret_cast<leo2_decode_plan*>(
        static_cast<uintptr_t>(1));
    require_result(leo2_decode_plan_create(NULL, valid_original_present,
        valid_recovery_present, &invalid_plan), LEO2_INVALID_ARGUMENT,
        "null codec plan create");
    require(invalid_plan == NULL,
        "null-codec plan creation did not clear output");
    invalid_plan = reinterpret_cast<leo2_decode_plan*>(
        static_cast<uintptr_t>(1));
    require_result(leo2_decode_plan_create(high.codec, NULL,
        valid_recovery_present, &invalid_plan), LEO2_INVALID_ARGUMENT,
        "null original-presence plan create");
    require(invalid_plan == NULL,
        "invalid plan creation did not clear output");
    invalid_plan = reinterpret_cast<leo2_decode_plan*>(
        static_cast<uintptr_t>(1));
    require_result(leo2_decode_plan_create(high.codec,
        valid_original_present, valid_recovery_present, NULL),
        LEO2_INVALID_ARGUMENT, "null plan output");
    counts->introspection_checks += 5;

    PlanOwner valid_plan;
    require_result(leo2_decode_plan_create(high.codec,
        valid_original_present, valid_recovery_present, &valid_plan.plan),
        LEO2_SUCCESS, "valid no-loss plan for scratch contracts");
    require_result(leo2_decode_plan_execute(valid_plan.plan, 0,
        NULL, NULL, NULL, NULL, 0), LEO2_SUCCESS,
        "zero-byte no-loss plan execute");
    require_result(leo2_decode(high.codec, 0,
        valid_original_present, valid_recovery_present,
        NULL, NULL, NULL, NULL, 0), LEO2_SUCCESS,
        "zero-byte one-shot no-loss decode");

    size_t scratch_bytes = 101;
    require_result(leo2_encode_scratch_size(NULL, 17, &scratch_bytes),
        LEO2_INVALID_ARGUMENT, "null codec encode scratch query");
    require(scratch_bytes == 0,
        "failed null-codec encode scratch query retained output");
    scratch_bytes = 102;
    require_result(leo2_encode_scratch_size(high.codec, 0, &scratch_bytes),
        LEO2_INVALID_ARGUMENT, "zero-byte encode scratch query");
    require(scratch_bytes == 0,
        "failed zero-byte encode scratch query retained output");
    require_result(leo2_encode_scratch_size(high.codec, 17, NULL),
        LEO2_INVALID_ARGUMENT, "null encode scratch output");
    scratch_bytes = 103;
    require_result(leo2_decode_scratch_size(NULL, 17, &scratch_bytes),
        LEO2_INVALID_ARGUMENT, "null codec decode scratch query");
    require(scratch_bytes == 0,
        "failed null-codec decode scratch query retained output");
    scratch_bytes = 104;
    require_result(leo2_decode_scratch_size(high.codec, 0, &scratch_bytes),
        LEO2_INVALID_ARGUMENT, "zero-byte decode scratch query");
    require(scratch_bytes == 0,
        "failed zero-byte decode scratch query retained output");
    require_result(leo2_decode_scratch_size(high.codec, 17, NULL),
        LEO2_INVALID_ARGUMENT, "null decode scratch output");
    scratch_bytes = 105;
    require_result(leo2_decode_plan_scratch_size(NULL, 17, &scratch_bytes),
        LEO2_INVALID_ARGUMENT, "null plan scratch query");
    require(scratch_bytes == 0,
        "failed null-plan scratch query retained output");
    scratch_bytes = 106;
    require_result(leo2_decode_plan_scratch_size(valid_plan.plan, 0,
        &scratch_bytes), LEO2_INVALID_ARGUMENT,
        "zero-byte plan scratch query");
    require(scratch_bytes == 0,
        "failed zero-byte plan scratch query retained output");
    require_result(leo2_decode_plan_scratch_size(valid_plan.plan, 17, NULL),
        LEO2_INVALID_ARGUMENT, "null plan scratch output");

    const uintptr_t invalid_size_address =
        std::numeric_limits<uintptr_t>::max() - (sizeof(size_t) - 1u);
    size_t* const invalid_size_output =
        reinterpret_cast<size_t*>(invalid_size_address);
    const uintptr_t invalid_u64_address =
        std::numeric_limits<uintptr_t>::max() - (sizeof(uint64_t) - 1u);
    uint64_t* const invalid_u64_output =
        reinterpret_cast<uint64_t*>(invalid_u64_address);
    require_result(leo2_codec_wire_shard_bytes(
        padded.codec, 1, invalid_u64_output), LEO2_INVALID_ARGUMENT,
        "unrepresentable wire-size output span");
    require_result(leo2_encode_scratch_size(
        high.codec, 17, invalid_size_output), LEO2_INVALID_ARGUMENT,
        "unrepresentable encode scratch output span");
    require_result(leo2_encode_batch_preflight_scratch_size(
        high.codec, 9, invalid_size_output), LEO2_INVALID_ARGUMENT,
        "unrepresentable encode-batch scratch output span");
    require_result(leo2_decode_plan_scratch_size(
        valid_plan.plan, 17, invalid_size_output), LEO2_INVALID_ARGUMENT,
        "unrepresentable plan scratch output span");
    require_result(leo2_decode_plan_batch_preflight_scratch_size(
        valid_plan.plan, 9, invalid_size_output), LEO2_INVALID_ARGUMENT,
        "unrepresentable decode-batch scratch output span");
    require_result(leo2_decode_scratch_size(
        high.codec, 17, invalid_size_output), LEO2_INVALID_ARGUMENT,
        "unrepresentable one-shot scratch output span");
    require_result(leo2_encode(NULL, 17, NULL, NULL, NULL, 0),
        LEO2_INVALID_ARGUMENT, "null codec encode");
    require_result(leo2_encode_batch(NULL, NULL, 0),
        LEO2_INVALID_ARGUMENT, "null codec empty batch");
    require_result(leo2_decode_plan_execute(NULL, 17, NULL, NULL, NULL,
        NULL, 0), LEO2_INVALID_ARGUMENT, "null plan execute");
    require_result(leo2_decode_plan_execute_batch(NULL, NULL, 0),
        LEO2_INVALID_ARGUMENT, "null plan empty batch");
    counts->introspection_checks += 29;
}

void test_default_affinity_thread_budget(Counts* counts)
{
#if defined(__linux__)
    cpu_set_t saved;
    CPU_ZERO(&saved);
    require(sched_getaffinity(0, sizeof(saved), &saved) == 0,
        "cannot read Linux affinity for default-thread test");

    cpu_set_t limited;
    CPU_ZERO(&limited);
    uint32_t selected = 0;
    for (int cpu = 0; cpu < CPU_SETSIZE && selected < 2; ++cpu)
    {
        if (CPU_ISSET(cpu, &saved))
        {
            CPU_SET(cpu, &limited);
            ++selected;
        }
    }
    require(selected != 0, "Linux affinity mask contains no CPUs");
    require(sched_setaffinity(0, sizeof(limited), &limited) == 0,
        "cannot constrain Linux affinity for default-thread test");

    leo2_context* raw_context = NULL;
    const leo2_result create_result = leo2_context_create(NULL, &raw_context);
    const int restore_result = sched_setaffinity(0, sizeof(saved), &saved);
    ContextOwner context;
    context.context = raw_context;
    require(restore_result == 0,
        "cannot restore Linux affinity after default-thread test");
    require_result(create_result, LEO2_SUCCESS,
        "affinity-constrained default context create");
    require(leo2_context_thread_count(context.context) == selected,
        "default context ignored the Linux allowed CPU set");
    require(leo2_test_context_worker_count(context.context) == 0,
        "unused default context eagerly started worker threads");
    counts->scheduler_checks += 3;
#else
    ContextOwner context;
    require_result(leo2_context_create(NULL, &context.context),
        LEO2_SUCCESS, "default context create");
    require(leo2_context_thread_count(context.context) >= 1 &&
            leo2_context_thread_count(context.context) <= 128,
        "default context thread budget is out of range");
    require(leo2_test_context_worker_count(context.context) == 0,
        "unused default context eagerly started worker threads");
    counts->scheduler_checks += 2;
#endif
}

void test_alias_and_scratch_contracts(leo2_context* context, Counts* counts)
{
    const size_t bytes = 17;
    CodecOwner xor_codec;
    require_result(leo2_codec_create(context, 3, 1,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL,
        &xor_codec.codec), LEO2_SUCCESS, "XOR codec create");

    Bytes a(bytes, 0);
    Bytes c(bytes, 0);
    for (size_t i = 0; i < bytes; ++i)
    {
        a[i] = static_cast<uint8_t>(i * 13u + 5u);
        c[i] = static_cast<uint8_t>(i * 29u + 11u);
    }
    const void* original[3] = { &a[0], &a[0], &c[0] };
    Bytes parity(bytes, 0xa5);
    void* recovery[1] = { &parity[0] };
    size_t encode_scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(xor_codec.codec, bytes,
        &encode_scratch_bytes), LEO2_SUCCESS, "encode scratch query");
    require(encode_scratch_bytes > 1, "encode scratch is unexpectedly tiny");
    AlignedBuffer encode_scratch(encode_scratch_bytes +
        leo2_scratch_alignment());
    require_result(leo2_encode(xor_codec.codec, bytes, original, recovery,
        encode_scratch.data(), encode_scratch_bytes), LEO2_SUCCESS,
        "allowed encode input-input alias");
    require(parity == c, "aliased XOR inputs produced the wrong parity");
    ++counts->alias_checks;

    require_result(leo2_encode(xor_codec.codec, bytes, original, recovery,
        encode_scratch.data(), encode_scratch_bytes - 1),
        LEO2_SCRATCH_TOO_SMALL, "short encode scratch");
    require_result(leo2_encode(xor_codec.codec, bytes, original, recovery,
        static_cast<uint8_t*>(encode_scratch.data()) + 1,
        encode_scratch_bytes), LEO2_BAD_ALIGNMENT,
        "unaligned encode scratch");
    counts->scratch_checks += 2;

    AlignedBuffer encode_collision(std::max(encode_scratch_bytes, bytes));
    memcpy(encode_collision.data(), &a[0], bytes);
    const void* overlapping_original[3] = {
        encode_collision.data(), &a[0], &c[0]
    };
    require_result(leo2_encode(xor_codec.codec, bytes, overlapping_original,
        recovery, encode_collision.data(), encode_scratch_bytes),
        LEO2_OVERLAP, "encode scratch/input overlap");
    void* overlapping_recovery[1] = { encode_collision.data() };
    require_result(leo2_encode(xor_codec.codec, bytes, original,
        overlapping_recovery, encode_collision.data(), encode_scratch_bytes),
        LEO2_OVERLAP, "encode scratch/output overlap");
    AlignedBuffer encode_tail_collision(encode_scratch_bytes + bytes);
    void* tail_recovery[1] = {
        static_cast<uint8_t*>(encode_tail_collision.data()) +
            encode_scratch_bytes
    };
    require_result(leo2_encode(xor_codec.codec, bytes, original,
        tail_recovery, encode_tail_collision.data(),
        encode_tail_collision.size()), LEO2_OVERLAP,
        "encode supplied scratch-tail/output overlap");
    counts->scratch_checks += 3;

    AlignedBuffer encode_metadata_scratch(encode_scratch_bytes);
    const void** scratch_original =
        static_cast<const void**>(encode_metadata_scratch.data());
    for (size_t i = 0; i < 3; ++i)
        scratch_original[i] = original[i];
    require_result(leo2_encode(xor_codec.codec, bytes, scratch_original,
        recovery, encode_metadata_scratch.data(),
        encode_metadata_scratch.size()), LEO2_OVERLAP,
        "encode scratch/original-descriptor overlap");
    void** scratch_recovery =
        static_cast<void**>(encode_metadata_scratch.data());
    scratch_recovery[0] = recovery[0];
    require_result(leo2_encode(xor_codec.codec, bytes, original,
        scratch_recovery, encode_metadata_scratch.data(),
        encode_metadata_scratch.size()), LEO2_OVERLAP,
        "encode scratch/recovery-descriptor overlap");

    void* output_overlaps_original_metadata[1] = {
        const_cast<void*>(static_cast<const void*>(original))
    };
    require_result(leo2_encode(xor_codec.codec, bytes, original,
        output_overlaps_original_metadata, encode_scratch.data(),
        encode_scratch_bytes), LEO2_OVERLAP,
        "encode output/original-descriptor overlap");
    void* output_overlaps_recovery_metadata[1];
    output_overlaps_recovery_metadata[0] =
        output_overlaps_recovery_metadata;
    require_result(leo2_encode(xor_codec.codec, bytes, original,
        output_overlaps_recovery_metadata, encode_scratch.data(),
        encode_scratch_bytes), LEO2_OVERLAP,
        "encode output/recovery-descriptor overlap");

    const uintptr_t near_address_end =
        std::numeric_limits<uintptr_t>::max() - sizeof(void*) + 1;
    const void* const* overflowing_original_metadata =
        reinterpret_cast<const void* const*>(near_address_end);
    require_result(leo2_encode(xor_codec.codec, bytes,
        overflowing_original_metadata, recovery, encode_scratch.data(),
        encode_scratch_bytes), LEO2_INVALID_ARGUMENT,
        "encode overflowing original-descriptor span");
    void* const* overflowing_recovery_metadata =
        reinterpret_cast<void* const*>(near_address_end);
    require_result(leo2_encode(xor_codec.codec, bytes, original,
        overflowing_recovery_metadata, encode_scratch.data(),
        encode_scratch_bytes), LEO2_INVALID_ARGUMENT,
        "encode overflowing recovery-descriptor span");
    require(parity == c,
        "rejected encode metadata overlap changed exact parity");
    counts->alias_checks += 7;

    CodecOwner codec;
    require_result(leo2_codec_create(context, 3, 2,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec.codec),
        LEO2_SUCCESS, "decode scratch codec create");
    Shards source(3, Bytes(bytes, 0));
    fill_shards(source, 0x12345678u);
    std::vector<const void*> source_pointers = const_pointers(source);
    Shards encoded(2, Bytes(bytes, 0));
    std::vector<void*> encoded_pointers = mutable_pointers(encoded);
    size_t source_scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(codec.codec, bytes,
        &source_scratch_bytes), LEO2_SUCCESS, "decode fixture encode query");
    AlignedBuffer source_scratch(source_scratch_bytes);
    require_result(leo2_encode(codec.codec, bytes, &source_pointers[0],
        &encoded_pointers[0], source_scratch.data(), source_scratch.size()),
        LEO2_SUCCESS, "decode fixture encode");

    const Shards exact_encoded = encoded;
    AlignedBuffer general_metadata_scratch(source_scratch_bytes);
    const void** general_scratch_original =
        static_cast<const void**>(general_metadata_scratch.data());
    void** general_scratch_recovery = reinterpret_cast<void**>(
        static_cast<uint8_t*>(general_metadata_scratch.data()) + 64);
    for (size_t i = 0; i < source_pointers.size(); ++i)
        general_scratch_original[i] = source_pointers[i];
    for (size_t i = 0; i < encoded_pointers.size(); ++i)
        general_scratch_recovery[i] = encoded_pointers[i];
    require_result(leo2_encode(codec.codec, bytes,
        general_scratch_original, general_scratch_recovery,
        general_metadata_scratch.data(), general_metadata_scratch.size()),
        LEO2_OVERLAP, "multi-output encode metadata/scratch overlap");
    require(encoded == exact_encoded,
        "rejected multi-output metadata overlap changed exact parity");
    counts->alias_checks += 2;

    uint8_t original_present[3] = { 0, 1, 1 };
    uint8_t recovery_present[2] = { 1, 0 };
    PlanOwner plan;
    require_result(leo2_decode_plan_create(codec.codec, original_present,
        recovery_present, &plan.plan), LEO2_SUCCESS, "scratch plan create");
    require(leo2_decode_plan_missing_original_count(plan.plan) == 1,
        "scratch plan missing count");
    size_t decode_scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(plan.plan, bytes,
        &decode_scratch_bytes), LEO2_SUCCESS, "decode scratch query");
    require(decode_scratch_bytes > 1, "decode scratch is unexpectedly tiny");
    AlignedBuffer decode_scratch(decode_scratch_bytes +
        leo2_scratch_alignment());
    const void* decode_original[3] = {
        NULL, source_pointers[1], source_pointers[2]
    };
    const void* decode_recovery[2] = { &encoded[0][0], NULL };
    Bytes restored(bytes, 0xcc);
    void* restored_original[3] = { &restored[0], NULL, NULL };

    require_result(leo2_decode_plan_execute(plan.plan, bytes, decode_original,
        decode_recovery, restored_original, decode_scratch.data(),
        decode_scratch_bytes - 1), LEO2_SCRATCH_TOO_SMALL,
        "short decode scratch");
    require_result(leo2_decode_plan_execute(plan.plan, bytes, decode_original,
        decode_recovery, restored_original,
        static_cast<uint8_t*>(decode_scratch.data()) + 1,
        decode_scratch_bytes), LEO2_BAD_ALIGNMENT,
        "unaligned decode scratch");
    counts->scratch_checks += 2;

    AlignedBuffer decode_collision(std::max(decode_scratch_bytes, bytes));
    memcpy(decode_collision.data(), source[1].data(), bytes);
    const void* overlapping_decode_original[3] = {
        NULL, decode_collision.data(), source_pointers[2]
    };
    require_result(leo2_decode_plan_execute(plan.plan, bytes,
        overlapping_decode_original, decode_recovery, restored_original,
        decode_collision.data(), decode_scratch_bytes), LEO2_OVERLAP,
        "decode scratch/input overlap");
    void* overlapping_restored[3] = {
        decode_collision.data(), NULL, NULL
    };
    require_result(leo2_decode_plan_execute(plan.plan, bytes, decode_original,
        decode_recovery, overlapping_restored, decode_collision.data(),
        decode_scratch_bytes), LEO2_OVERLAP,
        "decode scratch/output overlap");
    AlignedBuffer decode_tail_collision(decode_scratch_bytes + bytes);
    void* tail_restored[3] = {
        static_cast<uint8_t*>(decode_tail_collision.data()) +
            decode_scratch_bytes,
        NULL,
        NULL
    };
    require_result(leo2_decode_plan_execute(plan.plan, bytes, decode_original,
        decode_recovery, tail_restored, decode_tail_collision.data(),
        decode_tail_collision.size()), LEO2_OVERLAP,
        "decode supplied scratch-tail/output overlap");
    counts->scratch_checks += 3;

    AlignedBuffer decode_metadata_scratch(decode_scratch_bytes);
    const void** scratch_decode_original =
        static_cast<const void**>(decode_metadata_scratch.data());
    for (size_t i = 0; i < 3; ++i)
        scratch_decode_original[i] = decode_original[i];
    require_result(leo2_decode_plan_execute(plan.plan, bytes,
        scratch_decode_original, decode_recovery, restored_original,
        decode_metadata_scratch.data(), decode_metadata_scratch.size()),
        LEO2_OVERLAP, "decode scratch/original-descriptor overlap");
    const void** scratch_decode_recovery =
        static_cast<const void**>(decode_metadata_scratch.data());
    for (size_t i = 0; i < 2; ++i)
        scratch_decode_recovery[i] = decode_recovery[i];
    require_result(leo2_decode_plan_execute(plan.plan, bytes,
        decode_original, scratch_decode_recovery, restored_original,
        decode_metadata_scratch.data(), decode_metadata_scratch.size()),
        LEO2_OVERLAP, "decode scratch/recovery-descriptor overlap");
    void** scratch_decode_restored =
        static_cast<void**>(decode_metadata_scratch.data());
    for (size_t i = 0; i < 3; ++i)
        scratch_decode_restored[i] = restored_original[i];
    require_result(leo2_decode_plan_execute(plan.plan, bytes,
        decode_original, decode_recovery, scratch_decode_restored,
        decode_metadata_scratch.data(), decode_metadata_scratch.size()),
        LEO2_OVERLAP, "decode scratch/restored-descriptor overlap");

    void* restored_overlaps_original_metadata[3] = {
        const_cast<void*>(static_cast<const void*>(decode_original)),
        NULL, NULL
    };
    require_result(leo2_decode_plan_execute(plan.plan, bytes,
        decode_original, decode_recovery,
        restored_overlaps_original_metadata, decode_scratch.data(),
        decode_scratch_bytes), LEO2_OVERLAP,
        "decode output/original-descriptor overlap");
    void* restored_overlaps_recovery_metadata[3] = {
        const_cast<void*>(static_cast<const void*>(decode_recovery)),
        NULL, NULL
    };
    require_result(leo2_decode_plan_execute(plan.plan, bytes,
        decode_original, decode_recovery,
        restored_overlaps_recovery_metadata, decode_scratch.data(),
        decode_scratch_bytes), LEO2_OVERLAP,
        "decode output/recovery-descriptor overlap");
    void* restored_overlaps_own_metadata[3] = { NULL, NULL, NULL };
    restored_overlaps_own_metadata[0] = restored_overlaps_own_metadata;
    require_result(leo2_decode_plan_execute(plan.plan, bytes,
        decode_original, decode_recovery, restored_overlaps_own_metadata,
        decode_scratch.data(), decode_scratch_bytes), LEO2_OVERLAP,
        "decode output/restored-descriptor overlap");

    require_result(leo2_decode_plan_execute(plan.plan, bytes,
        reinterpret_cast<const void* const*>(near_address_end),
        decode_recovery, restored_original, decode_scratch.data(),
        decode_scratch_bytes), LEO2_INVALID_ARGUMENT,
        "decode overflowing original-descriptor span");
    require_result(leo2_decode_plan_execute(plan.plan, bytes,
        decode_original,
        reinterpret_cast<const void* const*>(near_address_end),
        restored_original, decode_scratch.data(), decode_scratch_bytes),
        LEO2_INVALID_ARGUMENT,
        "decode overflowing recovery-descriptor span");
    require_result(leo2_decode_plan_execute(plan.plan, bytes,
        decode_original, decode_recovery,
        reinterpret_cast<void* const*>(near_address_end),
        decode_scratch.data(), decode_scratch_bytes),
        LEO2_INVALID_ARGUMENT,
        "decode overflowing restored-descriptor span");
    counts->alias_checks += 9;

    require_result(leo2_decode_plan_execute(plan.plan, bytes, decode_original,
        decode_recovery, restored_original, decode_scratch.data(),
        decode_scratch_bytes), LEO2_SUCCESS, "valid decode after rejections");
    require(restored == source[0], "valid decode restored wrong bytes");
}

void test_one_shot_presence_alias_contracts(
    leo2_context* context,
    Counts* counts)
{
    const size_t bytes = 17;
    CodecOwner codec;
    require_result(leo2_codec_create(context, 4, 3,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec.codec),
        LEO2_SUCCESS, "one-shot presence codec create");

    Shards source(4, Bytes(bytes, 0));
    fill_shards(source, 0xc001d00du);
    std::vector<const void*> source_pointers = const_pointers(source);
    Shards encoded(3, Bytes(bytes, 0));
    std::vector<void*> encoded_pointers = mutable_pointers(encoded);
    size_t encode_scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(codec.codec, bytes,
        &encode_scratch_bytes), LEO2_SUCCESS,
        "one-shot presence encode scratch query");
    AlignedBuffer encode_scratch(encode_scratch_bytes);
    require_result(leo2_encode(codec.codec, bytes, &source_pointers[0],
        &encoded_pointers[0], encode_scratch.data(), encode_scratch.size()),
        LEO2_SUCCESS, "one-shot presence fixture encode");

    const uint8_t exact_original_present[4] = { 0, 0, 1, 1 };
    const uint8_t exact_recovery_present[3] = { 1, 1, 0 };
    const void* decode_original[4] = {
        NULL, NULL, source_pointers[2], source_pointers[3]
    };
    const void* decode_recovery[3] = {
        encoded_pointers[0], encoded_pointers[1], NULL
    };
    Bytes restored0(bytes, 0x91);
    Bytes restored1(bytes, 0x92);
    void* restored_original[4] = {
        &restored0[0], &restored1[0], NULL, NULL
    };
    size_t decode_scratch_bytes = 0;
    require_result(leo2_decode_scratch_size(codec.codec, bytes,
        &decode_scratch_bytes), LEO2_SUCCESS,
        "one-shot presence decode scratch query");
    AlignedBuffer decode_scratch(decode_scratch_bytes);

    /* A disjoint one-shot call remains legal and treats both presence arrays
       as immutable setup metadata. */
    uint8_t original_present[4];
    uint8_t recovery_present[3];
    memcpy(original_present, exact_original_present, sizeof(original_present));
    memcpy(recovery_present, exact_recovery_present, sizeof(recovery_present));
    require_result(leo2_decode(codec.codec, bytes, original_present,
        recovery_present, decode_original, decode_recovery,
        restored_original, decode_scratch.data(), decode_scratch.size()),
        LEO2_SUCCESS, "disjoint one-shot presence decode");
    require(restored0 == source[0] && restored1 == source[1],
        "disjoint one-shot presence decode restored wrong bytes");
    require(memcmp(original_present, exact_original_present,
                sizeof(original_present)) == 0 &&
            memcmp(recovery_present, exact_recovery_present,
                sizeof(recovery_present)) == 0,
        "disjoint one-shot decode changed presence metadata");
    ++counts->alias_checks;

    std::fill(restored0.begin(), restored0.end(), 0x81);
    std::fill(restored1.begin(), restored1.end(), 0x82);
    const Bytes exact_restored0 = restored0;
    const Bytes exact_restored1 = restored1;

    /* Multi-output validation uses scratch for sorted address ranges.  Reject
       before those writes when either presence array resides in scratch. */
    AlignedBuffer original_presence_scratch(decode_scratch_bytes + 1);
    memset(original_presence_scratch.data(), 0x6b,
        original_presence_scratch.size());
    uint8_t* const original_presence_inside_scratch =
        static_cast<uint8_t*>(original_presence_scratch.data()) + 1;
    memcpy(original_presence_inside_scratch, exact_original_present,
        sizeof(exact_original_present));
    Bytes original_presence_scratch_before(
        original_presence_scratch.size(), 0);
    memcpy(&original_presence_scratch_before[0],
        original_presence_scratch.data(), original_presence_scratch.size());
    require_result(leo2_decode(codec.codec, bytes,
        original_presence_inside_scratch,
        exact_recovery_present, decode_original, decode_recovery,
        restored_original, original_presence_scratch.data(),
        original_presence_scratch.size()), LEO2_OVERLAP,
        "one-shot scratch/original-presence overlap");
    require(memcmp(original_presence_scratch.data(),
                &original_presence_scratch_before[0],
                original_presence_scratch.size()) == 0 &&
            restored0 == exact_restored0 && restored1 == exact_restored1,
        "rejected scratch/original-presence overlap changed storage");

    AlignedBuffer recovery_presence_scratch(decode_scratch_bytes);
    memset(recovery_presence_scratch.data(), 0x7c,
        recovery_presence_scratch.size());
    memcpy(recovery_presence_scratch.data(), exact_recovery_present,
        sizeof(exact_recovery_present));
    Bytes recovery_presence_scratch_before(
        recovery_presence_scratch.size(), 0);
    memcpy(&recovery_presence_scratch_before[0],
        recovery_presence_scratch.data(), recovery_presence_scratch.size());
    require_result(leo2_decode(codec.codec, bytes, exact_original_present,
        static_cast<const uint8_t*>(recovery_presence_scratch.data()),
        decode_original, decode_recovery, restored_original,
        recovery_presence_scratch.data(), recovery_presence_scratch.size()),
        LEO2_OVERLAP, "one-shot scratch/recovery-presence overlap");
    require(memcmp(recovery_presence_scratch.data(),
                &recovery_presence_scratch_before[0],
                recovery_presence_scratch.size()) == 0 &&
            restored0 == exact_restored0 && restored1 == exact_restored1,
        "rejected scratch/recovery-presence overlap changed storage");

    /* A restored shard may not overwrite either presence bitmap.  Use a
       shard-sized backing allocation so the rejected range itself is valid. */
    Bytes original_presence_output(bytes, 0x4d);
    memcpy(&original_presence_output[0], exact_original_present,
        sizeof(exact_original_present));
    const Bytes original_presence_output_before = original_presence_output;
    void* output_on_original_presence[4] = {
        &original_presence_output[0], &restored1[0], NULL, NULL
    };
    require_result(leo2_decode(codec.codec, bytes,
        &original_presence_output[0], exact_recovery_present,
        decode_original, decode_recovery, output_on_original_presence,
        decode_scratch.data(), decode_scratch.size()), LEO2_OVERLAP,
        "one-shot output/original-presence overlap");
    require(original_presence_output == original_presence_output_before &&
            restored1 == exact_restored1,
        "rejected output/original-presence overlap changed storage");

    Bytes recovery_presence_output(bytes, 0x5e);
    memcpy(&recovery_presence_output[0], exact_recovery_present,
        sizeof(exact_recovery_present));
    const Bytes recovery_presence_output_before = recovery_presence_output;
    void* output_on_recovery_presence[4] = {
        &recovery_presence_output[0], &restored1[0], NULL, NULL
    };
    require_result(leo2_decode(codec.codec, bytes, exact_original_present,
        &recovery_presence_output[0], decode_original, decode_recovery,
        output_on_recovery_presence, decode_scratch.data(),
        decode_scratch.size()), LEO2_OVERLAP,
        "one-shot output/recovery-presence overlap");
    require(recovery_presence_output == recovery_presence_output_before &&
            restored1 == exact_restored1,
        "rejected output/recovery-presence overlap changed storage");

    /* Presence metadata joins, rather than replaces, the existing protected
       pointer-array metadata.  Exercise both scratch and output rejection in
       the one-shot wrapper and verify all caller-owned metadata is preserved. */
    AlignedBuffer descriptor_scratch(decode_scratch_bytes);
    memset(descriptor_scratch.data(), 0xa6, descriptor_scratch.size());
    const void** scratch_original =
        static_cast<const void**>(descriptor_scratch.data());
    for (size_t i = 0; i < 4; ++i)
        scratch_original[i] = decode_original[i];
    Bytes descriptor_scratch_before(descriptor_scratch.size(), 0);
    memcpy(&descriptor_scratch_before[0], descriptor_scratch.data(),
        descriptor_scratch.size());
    require_result(leo2_decode(codec.codec, bytes, exact_original_present,
        exact_recovery_present, scratch_original, decode_recovery,
        restored_original, descriptor_scratch.data(),
        descriptor_scratch.size()), LEO2_OVERLAP,
        "one-shot scratch/pointer-metadata overlap");
    require(memcmp(descriptor_scratch.data(), &descriptor_scratch_before[0],
                descriptor_scratch.size()) == 0 &&
            restored0 == exact_restored0 && restored1 == exact_restored1,
        "rejected scratch/pointer-metadata overlap changed storage");

    const void* decode_original_before[4];
    const void* decode_recovery_before[3];
    memcpy(decode_original_before, decode_original, sizeof(decode_original));
    memcpy(decode_recovery_before, decode_recovery, sizeof(decode_recovery));
    void* output_on_pointer_metadata[4] = {
        const_cast<void*>(static_cast<const void*>(decode_original)),
        &restored1[0], NULL, NULL
    };
    require_result(leo2_decode(codec.codec, bytes, exact_original_present,
        exact_recovery_present, decode_original, decode_recovery,
        output_on_pointer_metadata, decode_scratch.data(),
        decode_scratch.size()), LEO2_OVERLAP,
        "one-shot output/pointer-metadata overlap");
    require(memcmp(decode_original_before, decode_original,
                sizeof(decode_original)) == 0 &&
            memcmp(decode_recovery_before, decode_recovery,
                sizeof(decode_recovery)) == 0 &&
            restored1 == exact_restored1,
        "rejected output/pointer-metadata overlap changed storage");

    counts->alias_checks += 6;
}

void test_aligned_decode_input_staging_elision(
    leo2_context* context,
    Counts* counts)
{
    struct Case
    {
        uint32_t k;
        uint32_t r;
        leo2_profile profile;
        leo2_field field;
        uint32_t flags;
        size_t ragged_bytes;
    };
    const Case cases[] = {
        { 19, 5, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          LEO2_CODEC_FORCE_SPECIALIZED_DECODE, 63 },
        { 5, 19, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8,
          LEO2_CODEC_FORCE_SPECIALIZED_DECODE, 63 },
        { 19, 5, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16,
          LEO2_CODEC_FORCE_GENERIC_DECODE, 62 },
        { 5, 19, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16,
          LEO2_CODEC_FORCE_SPECIALIZED_DECODE, 62 }
    };

    for (size_t case_i = 0;
         case_i < sizeof(cases) / sizeof(cases[0]); ++case_i)
    {
        const Case& item = cases[case_i];
        leo2_codec_options options;
        memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.flags = item.flags;
        CodecOwner codec;
        require_result(leo2_codec_create(context, item.k, item.r,
            item.profile, item.field, &options, &codec.codec), LEO2_SUCCESS,
            "aligned-input scratch codec create");

        std::vector<uint8_t> original_present(item.k, 1);
        std::vector<uint8_t> recovery_present(item.r, 1);
        original_present[0] = 0;
        PlanOwner plan;
        require_result(leo2_decode_plan_create(codec.codec,
            &original_present[0], &recovery_present[0], &plan.plan),
            LEO2_SUCCESS, "aligned-input scratch plan create");

        size_t aligned_plan = 0;
        size_t ragged_plan = 0;
        size_t aligned_one_shot = 0;
        size_t ragged_one_shot = 0;
        require_result(leo2_decode_plan_scratch_size(
            plan.plan, 64, &aligned_plan), LEO2_SUCCESS,
            "aligned plan scratch query");
        require_result(leo2_decode_plan_scratch_size(
            plan.plan, item.ragged_bytes, &ragged_plan), LEO2_SUCCESS,
            "ragged plan scratch query");
        require_result(leo2_decode_scratch_size(
            codec.codec, 64, &aligned_one_shot), LEO2_SUCCESS,
            "aligned one-shot scratch query");
        require_result(leo2_decode_scratch_size(
            codec.codec, item.ragged_bytes, &ragged_one_shot), LEO2_SUCCESS,
            "ragged one-shot scratch query");

        const size_t staged_input_bytes =
            static_cast<size_t>(item.k + item.r) * 64;
        require(ragged_plan > aligned_plan &&
                ragged_plan - aligned_plan == staged_input_bytes,
            "aligned plan retained K+R input staging slots");
        require(ragged_one_shot > aligned_one_shot &&
                ragged_one_shot - aligned_one_shot == staged_input_bytes,
            "aligned one-shot retained K+R input staging slots");
        require(aligned_one_shot >= aligned_plan &&
                ragged_one_shot >= ragged_plan,
            "one-shot transform scratch does not cover the exact plan");
        const bool exact_high_output_retention =
            item.profile == LEO2_PROFILE_LEGACY_HIGH_V1 &&
            (item.flags & LEO2_CODEC_FORCE_GENERIC_DECODE) == 0;
        require(exact_high_output_retention
                ? aligned_one_shot > aligned_plan &&
                    ragged_one_shot > ragged_plan
                : aligned_one_shot == aligned_plan &&
                    ragged_one_shot == ragged_plan,
            "plan-specific transform scratch geometry is not deterministic");
        counts->scratch_checks += 4;

        const size_t byte_sizes[] = { 64, item.ragged_bytes };
        for (size_t size_i = 0;
             size_i < sizeof(byte_sizes) / sizeof(byte_sizes[0]); ++size_i)
        {
            const size_t bytes = byte_sizes[size_i];
            Shards source(item.k, Bytes(bytes, 0));
            fill_shards(source, static_cast<uint32_t>(
                0x51a9e000u + case_i * 17u + size_i));
            Shards recovery(item.r, Bytes(bytes, 0));
            std::vector<const void*> source_pointers = const_pointers(source);
            std::vector<void*> recovery_outputs = mutable_pointers(recovery);
            size_t encode_scratch_bytes = 0;
            require_result(leo2_encode_scratch_size(codec.codec, bytes,
                &encode_scratch_bytes), LEO2_SUCCESS,
                "input-staging fixture encode scratch query");
            AlignedBuffer encode_scratch(encode_scratch_bytes);
            require_result(leo2_encode(codec.codec, bytes,
                &source_pointers[0], &recovery_outputs[0],
                encode_scratch.data(), encode_scratch.size()), LEO2_SUCCESS,
                "input-staging fixture encode");

            const Shards source_before = source;
            const Shards recovery_before = recovery;
            std::vector<const void*> decode_original =
                const_pointers(source);
            std::vector<const void*> decode_recovery =
                const_pointers(recovery);
            decode_original[0] = NULL;
            Bytes restored(bytes, 0xcc);
            std::vector<void*> restored_original(item.k, NULL);
            restored_original[0] = &restored[0];
            const size_t execution_scratch_bytes =
                size_i == 0 ? aligned_plan : ragged_plan;
            AlignedBuffer execution_scratch(execution_scratch_bytes);
            require_result(leo2_decode_plan_execute(plan.plan, bytes,
                &decode_original[0], &decode_recovery[0],
                &restored_original[0], execution_scratch.data(),
                execution_scratch.size()), LEO2_SUCCESS,
                "input-staging transform decode");
            require(restored == source[0],
                "input-staging transform restored wrong bytes");
            require(source == source_before && recovery == recovery_before,
                "input-staging transform modified a caller input shard");
            ++counts->scratch_checks;
        }
    }
}

void test_tiled_decode_workspace_slopes(
    leo2_context* context,
    Counts* counts)
{
    struct Case
    {
        uint32_t k;
        uint32_t r;
        uint32_t missing;
        leo2_profile profile;
        leo2_field field;
        uint32_t flags;
    };
    const Case cases[] = {
        { 200, 9, 5, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          LEO2_CODEC_FORCE_SPECIALIZED_DECODE },
        { 9, 200, 5, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8,
          LEO2_CODEC_FORCE_SPECIALIZED_DECODE },
        { 1000, 17, 7, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16,
          LEO2_CODEC_FORCE_SPECIALIZED_DECODE },
        { 17, 1000, 7, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16,
          LEO2_CODEC_FORCE_SPECIALIZED_DECODE },
        { 128, 128, 128, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          LEO2_CODEC_FORCE_SPECIALIZED_DECODE },
        { 128, 128, 128, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8,
          LEO2_CODEC_FORCE_SPECIALIZED_DECODE },
        { 1000, 17, 7, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16,
          LEO2_CODEC_FORCE_GENERIC_DECODE }
    };

    for (size_t case_i = 0;
         case_i < sizeof(cases) / sizeof(cases[0]); ++case_i)
    {
        const Case& item = cases[case_i];
        leo2_codec_options options;
        memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.flags = item.flags;
        CodecOwner codec;
        require_result(leo2_codec_create(context, item.k, item.r,
            item.profile, item.field, &options, &codec.codec), LEO2_SUCCESS,
            "tiled-workspace codec create");

        std::vector<uint8_t> original_present(item.k, 1);
        std::vector<uint8_t> recovery_present(item.r, 1);
        for (uint32_t i = 0; i < item.missing; ++i)
            original_present[i] = 0;
        PlanOwner plan;
        require_result(leo2_decode_plan_create(codec.codec,
            &original_present[0], &recovery_present[0], &plan.plan),
            LEO2_SUCCESS, "tiled-workspace plan create");

        size_t plan_64 = 0;
        size_t plan_128 = 0;
        size_t one_shot_64 = 0;
        size_t one_shot_128 = 0;
        require_result(leo2_decode_plan_scratch_size(
            plan.plan, 64, &plan_64), LEO2_SUCCESS,
            "tiled-workspace 64-byte plan query");
        require_result(leo2_decode_plan_scratch_size(
            plan.plan, 128, &plan_128), LEO2_SUCCESS,
            "tiled-workspace 128-byte plan query");
        require_result(leo2_decode_scratch_size(
            codec.codec, 64, &one_shot_64), LEO2_SUCCESS,
            "tiled-workspace 64-byte one-shot query");
        require_result(leo2_decode_scratch_size(
            codec.codec, 128, &one_shot_128), LEO2_SUCCESS,
            "tiled-workspace 128-byte one-shot query");

        const size_t side = leo2_codec_padded_side(codec.codec);
        const size_t parent = leo2_codec_parent_count(codec.codec);
        const bool generic =
            (item.flags & LEO2_CODEC_FORCE_GENERIC_DECODE) != 0;
        const size_t tiled_plan_slots = item.profile == LEO2_PROFILE_LOW_V1
            ? side * 2
            : side * 2 + item.missing;
        const size_t tiled_one_shot_slots =
            item.profile == LEO2_PROFILE_LOW_V1
                ? side * 2
                : side * 2 + item.r;
        const size_t expected_plan_slots = generic
            ? parent
            : std::min(parent, tiled_plan_slots);
        const size_t expected_one_shot_slots = generic
            ? parent
            : std::min(parent, tiled_one_shot_slots);
        require(plan_128 > plan_64 &&
                plan_128 - plan_64 == expected_plan_slots * 64,
            "plan scratch does not have the declared side-sized slope");
        require(one_shot_128 > one_shot_64 &&
                one_shot_128 - one_shot_64 == expected_one_shot_slots * 64,
            "one-shot scratch does not have the declared worst-case slope");
        counts->scratch_checks += 2;
    }
}

void test_ragged_decode_tail_staging_slopes(
    leo2_context* context,
    Counts* counts)
{
    struct Case
    {
        uint32_t k;
        uint32_t r;
        uint32_t missing;
        leo2_profile profile;
        leo2_field field;
        uint32_t flags;
        size_t small_bytes;
        size_t large_bytes;
    };
    const Case cases[] = {
        { 200, 9, 5, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          LEO2_CODEC_FORCE_SPECIALIZED_DECODE, 65, 129 },
        { 9, 200, 5, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8,
          LEO2_CODEC_FORCE_SPECIALIZED_DECODE, 65, 129 },
        { 1000, 17, 7, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16,
          LEO2_CODEC_FORCE_SPECIALIZED_DECODE, 66, 130 },
        { 17, 1000, 7, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16,
          LEO2_CODEC_FORCE_SPECIALIZED_DECODE, 66, 130 },
        { 200, 9, 5, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
          LEO2_CODEC_FORCE_GENERIC_DECODE, 65, 129 }
    };

    for (size_t case_i = 0;
         case_i < sizeof(cases) / sizeof(cases[0]); ++case_i)
    {
        const Case& item = cases[case_i];
        leo2_codec_options options;
        memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.flags = item.flags;
        CodecOwner codec;
        require_result(leo2_codec_create(context, item.k, item.r,
            item.profile, item.field, &options, &codec.codec), LEO2_SUCCESS,
            "ragged-tail codec create");

        std::vector<uint8_t> original_present(item.k, 1);
        std::vector<uint8_t> recovery_present(item.r, 1);
        for (uint32_t i = 0; i < item.missing; ++i)
            original_present[i] = 0;
        PlanOwner plan;
        require_result(leo2_decode_plan_create(codec.codec,
            &original_present[0], &recovery_present[0], &plan.plan),
            LEO2_SUCCESS, "ragged-tail plan create");

        size_t small_plan = 0;
        size_t large_plan = 0;
        size_t small_one_shot = 0;
        size_t large_one_shot = 0;
        require_result(leo2_decode_plan_scratch_size(
            plan.plan, item.small_bytes, &small_plan), LEO2_SUCCESS,
            "ragged-tail small plan query");
        require_result(leo2_decode_plan_scratch_size(
            plan.plan, item.large_bytes, &large_plan), LEO2_SUCCESS,
            "ragged-tail large plan query");
        require_result(leo2_decode_scratch_size(
            codec.codec, item.small_bytes, &small_one_shot), LEO2_SUCCESS,
            "ragged-tail small one-shot query");
        require_result(leo2_decode_scratch_size(
            codec.codec, item.large_bytes, &large_one_shot), LEO2_SUCCESS,
            "ragged-tail large one-shot query");

        const size_t side = leo2_codec_padded_side(codec.codec);
        const size_t parent = leo2_codec_parent_count(codec.codec);
        const bool generic =
            (item.flags & LEO2_CODEC_FORCE_GENERIC_DECODE) != 0;
        const size_t tiled_plan_slots = item.profile == LEO2_PROFILE_LOW_V1
            ? side * 2
            : side * 2 + item.missing;
        const size_t tiled_one_shot_slots =
            item.profile == LEO2_PROFILE_LOW_V1
                ? side * 2
                : side * 2 + item.r;
        const size_t plan_slots = generic
            ? parent
            : std::min(parent, tiled_plan_slots);
        const size_t one_shot_slots = generic
            ? parent
            : std::min(parent, tiled_one_shot_slots);
        require(large_plan > small_plan &&
                large_plan - small_plan == plan_slots * 64,
            "ragged plan scratch still scales K+R staging with shard bytes");
        require(large_one_shot > small_one_shot &&
                large_one_shot - small_one_shot == one_shot_slots * 64,
            "ragged one-shot scratch still scales K+R staging with shard bytes");
        counts->scratch_checks += 2;

        const size_t byte_sizes[] = { item.small_bytes, item.large_bytes };
        for (size_t size_i = 0;
             size_i < sizeof(byte_sizes) / sizeof(byte_sizes[0]); ++size_i)
        {
            const size_t bytes = byte_sizes[size_i];
            Shards source(item.k, Bytes(bytes, 0));
            fill_shards(source, static_cast<uint32_t>(
                0x72b41000u + case_i * 37u + size_i));
            Shards recovery(item.r, Bytes(bytes, 0));
            std::vector<const void*> source_pointers = const_pointers(source);
            std::vector<void*> recovery_outputs = mutable_pointers(recovery);
            size_t encode_scratch_bytes = 0;
            require_result(leo2_encode_scratch_size(
                codec.codec, bytes, &encode_scratch_bytes), LEO2_SUCCESS,
                "ragged-tail encode scratch query");
            AlignedBuffer encode_scratch(encode_scratch_bytes);
            require_result(leo2_encode(codec.codec, bytes,
                &source_pointers[0], &recovery_outputs[0],
                encode_scratch.data(), encode_scratch.size()), LEO2_SUCCESS,
                "ragged-tail fixture encode");

            std::vector<const void*> decode_original =
                const_pointers(source);
            std::vector<const void*> decode_recovery =
                const_pointers(recovery);
            Shards restored(item.missing, Bytes(bytes, 0xcc));
            std::vector<void*> restored_original(item.k, NULL);
            for (uint32_t i = 0; i < item.missing; ++i)
            {
                decode_original[i] = NULL;
                restored_original[i] = &restored[i][0];
            }
            const size_t execution_scratch_bytes =
                size_i == 0 ? small_plan : large_plan;
            AlignedBuffer execution_scratch(execution_scratch_bytes);
            require_result(leo2_decode_plan_execute(plan.plan, bytes,
                &decode_original[0], &decode_recovery[0],
                &restored_original[0], execution_scratch.data(),
                execution_scratch.size()), LEO2_SUCCESS,
                "ragged-tail split decode");
            for (uint32_t i = 0; i < item.missing; ++i)
                require(restored[i] == source[i],
                    "ragged-tail split decode restored wrong bytes");
            ++counts->scratch_checks;
        }
    }
}

void test_ragged_encode_tail_staging_slopes(
    leo2_context* context,
    Counts* counts)
{
    struct Case
    {
        uint32_t k;
        uint32_t r;
        leo2_profile profile;
        leo2_field field;
        size_t small_bytes;
        size_t large_bytes;
    };
    const Case cases[] = {
        { 200, 9, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 65, 129 },
        { 9, 200, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8, 65, 129 },
        { 1000, 17, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, 66, 130 },
        { 17, 1000, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16, 66, 130 }
    };

    for (size_t case_i = 0;
         case_i < sizeof(cases) / sizeof(cases[0]); ++case_i)
    {
        const Case& item = cases[case_i];
        CodecOwner codec;
        require_result(leo2_codec_create(context, item.k, item.r,
            item.profile, item.field, NULL, &codec.codec), LEO2_SUCCESS,
            "ragged-encode codec create");

        size_t aligned = 0;
        size_t small = 0;
        size_t large = 0;
        require_result(leo2_encode_scratch_size(
            codec.codec, 64, &aligned), LEO2_SUCCESS,
            "ragged-encode aligned scratch query");
        require_result(leo2_encode_scratch_size(
            codec.codec, item.small_bytes, &small), LEO2_SUCCESS,
            "ragged-encode small scratch query");
        require_result(leo2_encode_scratch_size(
            codec.codec, item.large_bytes, &large), LEO2_SUCCESS,
            "ragged-encode large scratch query");

        const size_t work_slots =
            static_cast<size_t>(leo2_codec_padded_side(codec.codec)) * 2;
        const size_t staging_slots = item.k +
            (item.profile == LEO2_PROFILE_LOW_V1 ? item.r : 0);
        require(small > aligned &&
                small - aligned == staging_slots * 64,
            "ragged encode did not use fixed one-tile staging");
        require(large > small &&
                large - small == work_slots * 64,
            "ragged encode staging still scales with full shard bytes");
        counts->scratch_checks += 2;
    }
}

struct BatchFixture
{
    BatchFixture(leo2_context* context, size_t bytes)
        : bytes(bytes)
        , original(3, Bytes(bytes, 0))
        , recovery(2, Bytes(bytes, 0))
        , codec(NULL)
        , plan(NULL)
        , encode_scratch_bytes(0)
        , decode_scratch_bytes(0)
    {
        fill_shards(original, 0x9e3779b9u);
        original_pointers = const_pointers(original);
        recovery_pointers = mutable_pointers(recovery);
        require_result(leo2_codec_create(context, 3, 2,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec),
            LEO2_SUCCESS, "batch codec create");
        require_result(leo2_encode_scratch_size(codec, bytes,
            &encode_scratch_bytes), LEO2_SUCCESS, "batch encode scratch query");
        AlignedBuffer scratch(encode_scratch_bytes);
        require_result(leo2_encode(codec, bytes, &original_pointers[0],
            &recovery_pointers[0], scratch.data(), scratch.size()),
            LEO2_SUCCESS, "batch reference encode");

        const uint8_t original_present[3] = { 0, 1, 1 };
        const uint8_t recovery_present[2] = { 1, 0 };
        require_result(leo2_decode_plan_create(codec, original_present,
            recovery_present, &plan), LEO2_SUCCESS, "batch plan create");
        require_result(leo2_decode_plan_scratch_size(plan, bytes,
            &decode_scratch_bytes), LEO2_SUCCESS, "batch decode scratch query");
        require(encode_scratch_bytes > 1 && decode_scratch_bytes > 1,
            "batch scratch unexpectedly tiny");
    }

    ~BatchFixture()
    {
        leo2_decode_plan_destroy(plan);
        leo2_codec_destroy(codec);
    }

    const size_t bytes;
    Shards original;
    Shards recovery;
    std::vector<const void*> original_pointers;
    std::vector<void*> recovery_pointers;
    leo2_codec* codec;
    leo2_decode_plan* plan;
    size_t encode_scratch_bytes;
    size_t decode_scratch_bytes;

private:
    BatchFixture(const BatchFixture&);
    BatchFixture& operator=(const BatchFixture&);
};

void run_encode_batch_call(
    const BatchFixture& fixture,
    unsigned call_index,
    size_t item_count = 4)
{
    std::vector<Shards> output(item_count,
        Shards(2, Bytes(fixture.bytes, 0xa5)));
    std::vector<std::vector<void*> > output_pointers(item_count);
    std::vector<std::unique_ptr<AlignedBuffer> > scratch(item_count);
    std::vector<leo2_encode_batch_item> items(item_count);
    for (size_t item = 0; item < item_count; ++item)
    {
        output_pointers[item] = mutable_pointers(output[item]);
        scratch[item].reset(new AlignedBuffer(fixture.encode_scratch_bytes));
        items[item].shard_bytes = fixture.bytes;
        items[item].original = &fixture.original_pointers[0];
        items[item].recovery = &output_pointers[item][0];
        items[item].scratch = scratch[item]->data();
        items[item].scratch_bytes = scratch[item]->size();
    }
    require_result(leo2_encode_batch(fixture.codec, &items[0], items.size()),
        LEO2_SUCCESS, "concurrent encode batch " + std::to_string(call_index));
    for (size_t item = 0; item < item_count; ++item)
        require(output[item] == fixture.recovery,
            "concurrent encode batch parity mismatch");
}

void run_decode_batch_call(const BatchFixture& fixture, unsigned call_index)
{
    const size_t item_count = 4;
    const void* original[3] = {
        NULL, fixture.original_pointers[1], fixture.original_pointers[2]
    };
    const void* recovery[2] = { &fixture.recovery[0][0], NULL };
    std::vector<Bytes> restored(item_count, Bytes(fixture.bytes, 0xcc));
    std::vector<std::vector<void*> > output_pointers(item_count,
        std::vector<void*>(3, NULL));
    std::vector<std::unique_ptr<AlignedBuffer> > scratch(item_count);
    std::vector<leo2_decode_batch_item> items(item_count);
    for (size_t item = 0; item < item_count; ++item)
    {
        output_pointers[item][0] = &restored[item][0];
        scratch[item].reset(new AlignedBuffer(fixture.decode_scratch_bytes));
        items[item].shard_bytes = fixture.bytes;
        items[item].original = original;
        items[item].recovery = recovery;
        items[item].restored_original = &output_pointers[item][0];
        items[item].scratch = scratch[item]->data();
        items[item].scratch_bytes = scratch[item]->size();
    }
    require_result(leo2_decode_plan_execute_batch(fixture.plan, &items[0],
        items.size()), LEO2_SUCCESS,
        "concurrent decode batch " + std::to_string(call_index));
    for (size_t item = 0; item < item_count; ++item)
        require(restored[item] == fixture.original[0],
            "concurrent decode batch restored wrong bytes");
}

void test_concurrent_shared_context_batches(
    leo2_context* context,
    Counts* counts)
{
    const BatchFixture fixture(context, 33);
    require(leo2_test_context_worker_count(context) == 0,
        "context eagerly started workers before its first parallel batch");
    require_result(leo2_encode_batch(fixture.codec, NULL, 0),
        LEO2_SUCCESS, "empty batch before lazy start");
    run_encode_batch_call(fixture, 0, 1);
    require(leo2_test_context_worker_count(context) == 0,
        "single-item batch unnecessarily started worker threads");
    counts->scheduler_checks += 3;

    const unsigned caller_count = 8;
    const unsigned repetitions = 4;
    std::atomic<unsigned> ready(0);
    std::atomic<bool> start(false);
    std::vector<std::string> errors(caller_count);
    std::vector<std::thread> callers;
    for (unsigned caller = 0; caller < caller_count; ++caller)
    {
        callers.push_back(std::thread([&, caller]() {
            ready.fetch_add(1, std::memory_order_release);
            while (!start.load(std::memory_order_acquire))
                std::this_thread::yield();
            try
            {
                for (unsigned repetition = 0; repetition < repetitions;
                     ++repetition)
                {
                    const unsigned index = caller * repetitions + repetition;
                    if ((caller & 1u) == 0)
                        run_encode_batch_call(fixture, index);
                    else
                        run_decode_batch_call(fixture, index);
                }
            }
            catch (const std::exception& error)
            {
                errors[caller] = error.what();
            }
        }));
    }
    while (ready.load(std::memory_order_acquire) != caller_count)
        std::this_thread::yield();
    start.store(true, std::memory_order_release);
    for (size_t i = 0; i < callers.size(); ++i)
        callers[i].join();
    for (unsigned caller = 0; caller < caller_count; ++caller)
        require(errors[caller].empty(), errors[caller]);
    require(leo2_test_context_worker_count(context) == 3,
        "concurrent first use did not start exactly one persistent pool");
    ++counts->scheduler_checks;
    counts->concurrent_batches +=
        static_cast<uint64_t>(caller_count) * repetitions;
}

void test_lazy_thread_start_failure(Counts* counts)
{
    leo2_context_options options;
    memset(&options, 0, sizeof(options));
    options.struct_size = sizeof(options);
    options.backend = LEO2_BACKEND_AUTO;
    options.thread_count = 4;
    ContextOwner context;
    require_result(leo2_context_create(&options, &context.context),
        LEO2_SUCCESS, "lazy-failure context create");
    const BatchFixture fixture(context.context, 65);

    const size_t item_count = 2;
    std::vector<Shards> output(item_count,
        Shards(2, Bytes(fixture.bytes, 0xa5)));
    std::vector<std::vector<void*> > output_pointers(item_count);
    std::vector<std::unique_ptr<AlignedBuffer> > scratch(item_count);
    leo2_encode_batch_item items[item_count];
    for (size_t item = 0; item < item_count; ++item)
    {
        output_pointers[item] = mutable_pointers(output[item]);
        scratch[item].reset(new AlignedBuffer(fixture.encode_scratch_bytes));
        items[item].shard_bytes = fixture.bytes;
        items[item].original = &fixture.original_pointers[0];
        items[item].recovery = &output_pointers[item][0];
        items[item].scratch = scratch[item]->data();
        items[item].scratch_bytes = scratch[item]->size();
    }

    leo2_test_set_thread_start_fault(1);
    static const unsigned caller_count = 8;
    std::atomic<unsigned> ready(0);
    std::atomic<bool> start(false);
    std::vector<leo2_result> results(caller_count, LEO2_SUCCESS);
    std::vector<std::thread> callers;
    for (unsigned caller = 0; caller < caller_count; ++caller)
    {
        callers.push_back(std::thread([&, caller]() {
            ready.fetch_add(1, std::memory_order_release);
            while (!start.load(std::memory_order_acquire))
                std::this_thread::yield();
            results[caller] = leo2_encode_batch(
                fixture.codec, items, item_count);
        }));
    }
    while (ready.load(std::memory_order_acquire) != caller_count)
        std::this_thread::yield();
    start.store(true, std::memory_order_release);
    for (size_t caller = 0; caller < callers.size(); ++caller)
        callers[caller].join();
    for (unsigned caller = 0; caller < caller_count; ++caller)
        require_result(results[caller], LEO2_OUT_OF_MEMORY,
            "concurrent injected lazy thread-start failure");
    require(leo2_test_thread_start_fault_consumptions() == 1,
        "lazy thread-start fault was not consumed exactly once");
    require(leo2_test_context_worker_count(context.context) == 0,
        "failed lazy start retained partial worker threads");
    require_result(leo2_encode_batch(fixture.codec, items, item_count),
        LEO2_OUT_OF_MEMORY, "cached lazy thread-start failure");
    require(leo2_test_thread_start_fault_consumptions() == 1,
        "cached lazy failure retried thread creation");
    leo2_test_set_thread_start_fault(0);
    counts->scheduler_checks += 5;
}

void test_incremental_lazy_thread_growth(Counts* counts)
{
    leo2_context_options options;
    memset(&options, 0, sizeof(options));
    options.struct_size = sizeof(options);
    options.backend = LEO2_BACKEND_AUTO;
    options.thread_count = 8;
    ContextOwner context;
    require_result(leo2_context_create(&options, &context.context),
        LEO2_SUCCESS, "lazy-growth context create");
    const BatchFixture fixture(context.context, 129);

    run_encode_batch_call(fixture, 0, 2);
    require(leo2_test_context_worker_count(context.context) == 1,
        "two-item batch started unused workers");
    run_encode_batch_call(fixture, 1, 4);
    require(leo2_test_context_worker_count(context.context) == 3,
        "larger batch did not grow the persistent pool exactly");
    run_encode_batch_call(fixture, 2, 2);
    require(leo2_test_context_worker_count(context.context) == 3,
        "smaller batch discarded persistent workers");
    counts->scheduler_checks += 3;
}

void test_deterministic_batch_failures(leo2_context* context, Counts* counts)
{
    const BatchFixture fixture(context, 17);
    const size_t repetitions = 64;

#if SIZE_MAX > UINT32_MAX
    require_result(leo2_encode_batch(fixture.codec,
        reinterpret_cast<const leo2_encode_batch_item*>(
            static_cast<uintptr_t>(1)),
        static_cast<size_t>(UINT32_MAX) + 1u), LEO2_INVALID_ARGUMENT,
        "oversized encode batch");
    require_result(leo2_decode_plan_execute_batch(fixture.plan,
        reinterpret_cast<const leo2_decode_batch_item*>(
            static_cast<uintptr_t>(1)),
        static_cast<size_t>(UINT32_MAX) + 1u), LEO2_INVALID_ARGUMENT,
        "oversized decode batch");
    counts->batch_failure_checks += 2;
#endif

    const uintptr_t encode_item_address_end =
        std::numeric_limits<uintptr_t>::max() -
        sizeof(leo2_encode_batch_item) + 1;
    require_result(leo2_encode_batch(fixture.codec,
        reinterpret_cast<const leo2_encode_batch_item*>(
            encode_item_address_end), 1), LEO2_INVALID_ARGUMENT,
        "overflowing encode batch-item span");
    const uintptr_t decode_item_address_end =
        std::numeric_limits<uintptr_t>::max() -
        sizeof(leo2_decode_batch_item) + 1;
    require_result(leo2_decode_plan_execute_batch(fixture.plan,
        reinterpret_cast<const leo2_decode_batch_item*>(
            decode_item_address_end), 1), LEO2_INVALID_ARGUMENT,
        "overflowing decode batch-item span");
    counts->batch_failure_checks += 2;

    Shards encode_output(2, Bytes(fixture.bytes, 0));
    std::vector<void*> encode_output_pointers = mutable_pointers(encode_output);
    AlignedBuffer encode_short(fixture.encode_scratch_bytes);
    AlignedBuffer encode_unaligned(fixture.encode_scratch_bytes + 1);
    AlignedBuffer encode_valid(fixture.encode_scratch_bytes);
    leo2_encode_batch_item encode_items[3];
    memset(encode_items, 0, sizeof(encode_items));
    for (size_t i = 0; i < 3; ++i)
    {
        encode_items[i].shard_bytes = fixture.bytes;
        encode_items[i].original = &fixture.original_pointers[0];
        encode_items[i].recovery = &encode_output_pointers[0];
    }
    encode_items[0].scratch = encode_short.data();
    encode_items[0].scratch_bytes = fixture.encode_scratch_bytes - 1;
    encode_items[1].scratch =
        static_cast<uint8_t*>(encode_unaligned.data()) + 1;
    encode_items[1].scratch_bytes = fixture.encode_scratch_bytes;
    encode_items[2].scratch = encode_valid.data();
    encode_items[2].scratch_bytes = fixture.encode_scratch_bytes;
    for (size_t repeat = 0; repeat < repetitions; ++repeat)
    {
        require_result(leo2_encode_batch(fixture.codec, encode_items, 3),
            LEO2_SCRATCH_TOO_SMALL,
            "lowest-index encode batch failure");
        ++counts->batch_failure_checks;
    }

    leo2_encode_batch_item* encode_item_in_scratch =
        static_cast<leo2_encode_batch_item*>(encode_valid.data());
    encode_item_in_scratch->shard_bytes = fixture.bytes;
    encode_item_in_scratch->original = &fixture.original_pointers[0];
    encode_item_in_scratch->recovery = &encode_output_pointers[0];
    encode_item_in_scratch->scratch = encode_valid.data();
    encode_item_in_scratch->scratch_bytes = fixture.encode_scratch_bytes;
    require_result(leo2_encode_batch(
        fixture.codec, encode_item_in_scratch, 1), LEO2_OVERLAP,
        "encode batch-item/scratch overlap");
    leo2_encode_batch_item encode_item_output_overlap;
    memset(&encode_item_output_overlap, 0, sizeof(encode_item_output_overlap));
    void* encode_item_self_output[2] = {
        &encode_item_output_overlap, NULL
    };
    encode_item_output_overlap.shard_bytes = fixture.bytes;
    encode_item_output_overlap.original = &fixture.original_pointers[0];
    encode_item_output_overlap.recovery = encode_item_self_output;
    encode_item_output_overlap.scratch = encode_valid.data();
    encode_item_output_overlap.scratch_bytes =
        fixture.encode_scratch_bytes;
    require_result(leo2_encode_batch(
        fixture.codec, &encode_item_output_overlap, 1), LEO2_OVERLAP,
        "encode batch-item/output overlap");
    counts->batch_failure_checks += 2;

    const void* decode_original[3] = {
        NULL, fixture.original_pointers[1], fixture.original_pointers[2]
    };
    const void* decode_recovery[2] = { &fixture.recovery[0][0], NULL };
    Bytes restored(fixture.bytes, 0);
    void* restored_pointers[3] = { &restored[0], NULL, NULL };
    AlignedBuffer decode_unaligned(fixture.decode_scratch_bytes + 1);
    AlignedBuffer decode_short(fixture.decode_scratch_bytes);
    AlignedBuffer decode_valid(fixture.decode_scratch_bytes);
    leo2_decode_batch_item decode_items[3];
    memset(decode_items, 0, sizeof(decode_items));
    for (size_t i = 0; i < 3; ++i)
    {
        decode_items[i].shard_bytes = fixture.bytes;
        decode_items[i].original = decode_original;
        decode_items[i].recovery = decode_recovery;
        decode_items[i].restored_original = restored_pointers;
    }
    decode_items[0].scratch =
        static_cast<uint8_t*>(decode_unaligned.data()) + 1;
    decode_items[0].scratch_bytes = fixture.decode_scratch_bytes;
    decode_items[1].scratch = decode_short.data();
    decode_items[1].scratch_bytes = fixture.decode_scratch_bytes - 1;
    decode_items[2].scratch = decode_valid.data();
    decode_items[2].scratch_bytes = fixture.decode_scratch_bytes;
    for (size_t repeat = 0; repeat < repetitions; ++repeat)
    {
        require_result(leo2_decode_plan_execute_batch(
            fixture.plan, decode_items, 3), LEO2_BAD_ALIGNMENT,
            "lowest-index decode batch failure");
        ++counts->batch_failure_checks;
    }

    leo2_decode_batch_item* decode_item_in_scratch =
        static_cast<leo2_decode_batch_item*>(decode_valid.data());
    decode_item_in_scratch->shard_bytes = fixture.bytes;
    decode_item_in_scratch->original = decode_original;
    decode_item_in_scratch->recovery = decode_recovery;
    decode_item_in_scratch->restored_original = restored_pointers;
    decode_item_in_scratch->scratch = decode_valid.data();
    decode_item_in_scratch->scratch_bytes = fixture.decode_scratch_bytes;
    require_result(leo2_decode_plan_execute_batch(
        fixture.plan, decode_item_in_scratch, 1), LEO2_OVERLAP,
        "decode batch-item/scratch overlap");
    leo2_decode_batch_item decode_item_output_overlap;
    memset(&decode_item_output_overlap, 0, sizeof(decode_item_output_overlap));
    void* decode_item_self_output[3] = {
        &decode_item_output_overlap, NULL, NULL
    };
    decode_item_output_overlap.shard_bytes = fixture.bytes;
    decode_item_output_overlap.original = decode_original;
    decode_item_output_overlap.recovery = decode_recovery;
    decode_item_output_overlap.restored_original =
        decode_item_self_output;
    decode_item_output_overlap.scratch = decode_valid.data();
    decode_item_output_overlap.scratch_bytes =
        fixture.decode_scratch_bytes;
    require_result(leo2_decode_plan_execute_batch(
        fixture.plan, &decode_item_output_overlap, 1), LEO2_OVERLAP,
        "decode batch-item/output overlap");
    counts->batch_failure_checks += 2;

    Shards exact_batch_output(2, Bytes(fixture.bytes, 0xa5));
    std::vector<void*> exact_batch_output_pointers =
        mutable_pointers(exact_batch_output);
    leo2_encode_batch_item valid_encode_item;
    valid_encode_item.shard_bytes = fixture.bytes;
    valid_encode_item.original = &fixture.original_pointers[0];
    valid_encode_item.recovery = &exact_batch_output_pointers[0];
    valid_encode_item.scratch = encode_valid.data();
    valid_encode_item.scratch_bytes = fixture.encode_scratch_bytes;
    require_result(leo2_encode_batch(
        fixture.codec, &valid_encode_item, 1), LEO2_SUCCESS,
        "valid encode batch after metadata rejections");
    require(exact_batch_output == fixture.recovery,
        "valid encode batch parity changed after metadata rejections");

    Bytes exact_batch_restored(fixture.bytes, 0xcc);
    void* exact_batch_restored_pointers[3] = {
        &exact_batch_restored[0], NULL, NULL
    };
    leo2_decode_batch_item valid_decode_item;
    valid_decode_item.shard_bytes = fixture.bytes;
    valid_decode_item.original = decode_original;
    valid_decode_item.recovery = decode_recovery;
    valid_decode_item.restored_original = exact_batch_restored_pointers;
    valid_decode_item.scratch = decode_valid.data();
    valid_decode_item.scratch_bytes = fixture.decode_scratch_bytes;
    require_result(leo2_decode_plan_execute_batch(
        fixture.plan, &valid_decode_item, 1), LEO2_SUCCESS,
        "valid decode batch after metadata rejections");
    require(exact_batch_restored == fixture.original[0],
        "valid decode batch changed after metadata rejections");
}

void require_legacy_result(
    LeopardResult actual,
    LeopardResult expected,
    const std::string& operation)
{
    if (actual == expected)
        return;
    throw std::runtime_error(operation + ": got " + leo_result_string(actual) +
        ", expected " + leo_result_string(expected));
}

void test_legacy_negative_contract(Counts* counts)
{
    require(leopard::SIMDSafeAllocate(SIZE_MAX) == NULL,
        "legacy aligned allocator accepted maximum overflowing size");
    require(leopard::SIMDSafeAllocate(
                SIZE_MAX - leopard::kAlignmentBytes + 1) == NULL,
        "legacy aligned allocator accepted first overflowing size");

    const size_t aligned_probe_bytes = leopard::kAlignmentBytes + 1;
    uint8_t* aligned_probe = leopard::SIMDSafeAllocate(aligned_probe_bytes);
    require(aligned_probe != NULL,
        "legacy aligned allocator rejected a normal size");
    require(reinterpret_cast<uintptr_t>(aligned_probe) %
                leopard::kAlignmentBytes == 0,
        "legacy aligned allocator returned a misaligned pointer");
    for (size_t i = 0; i < aligned_probe_bytes; ++i)
        require(aligned_probe[i] == 0,
            "legacy aligned allocator did not zero-fill its payload");
    aligned_probe[aligned_probe_bytes - 1] = 0x5a;
    leopard::SIMDSafeFree(aligned_probe);

    const uint64_t simulated_32bit_limit = UINT32_MAX;
    const uint64_t largest_32bit_aligned =
        simulated_32bit_limit & ~UINT64_C(63);
    require(leopard::LegacyBufferBytesValidForAddressLimit(
                largest_32bit_aligned, simulated_32bit_limit),
        "legacy simulated 32-bit aligned size rejected");
    require(!leopard::LegacyBufferBytesValidForAddressLimit(
                simulated_32bit_limit + 1, simulated_32bit_limit),
        "legacy simulated 32-bit overflow size accepted");
    require(!leopard::LegacyBufferBytesValidForAddressLimit(
                0, simulated_32bit_limit) &&
            !leopard::LegacyBufferBytesValidForAddressLimit(
                65, simulated_32bit_limit),
        "legacy size predicate accepted invalid granularity");

#if SIZE_MAX < UINT64_MAX
    const uint64_t unaddressable_bytes =
        static_cast<uint64_t>(SIZE_MAX) + 1;
    require_legacy_result(leo_encode(
        unaddressable_bytes, 1, 1, 1, NULL, NULL),
        Leopard_InvalidSize, "legacy 32-bit encode size-width rejection");
    require_legacy_result(leo_decode(
        unaddressable_bytes, 1, 1, 1, NULL, NULL, NULL),
        Leopard_InvalidSize, "legacy 32-bit decode size-width rejection");
#endif

    require(LEO_VERSION == 2, "unexpected legacy API version");
    require(leo_init() == Leopard_Success, "legacy reinitialization failed");
    ++counts->legacy_checks;

    uint8_t shard[64] = { 0 };
    const void* original[3] = { shard, shard, shard };
    void* encode_work[4] = { shard, shard, shard, shard };
    require(leo_encode_work_count(3, 2) == 4,
        "legacy encode work-count changed");
    require(leo_encode_work_count(1, 1) == 1 &&
            leo_encode_work_count(65535, 1) == 1 &&
            leo_encode_work_count(32768, 32768) == 65536,
        "legacy encode work-count rejected a valid boundary");
    require(leo_encode_work_count(0, 1) == 0 &&
            leo_encode_work_count(2, 3) == 0 &&
            leo_encode_work_count(65536, 1) == 0 &&
            leo_encode_work_count(32769, 32767) == 0 &&
            leo_encode_work_count(std::numeric_limits<unsigned>::max(), 2) == 0,
        "legacy encode work-count accepted invalid counts");
    require_legacy_result(leo_encode(17, 3, 2, 4, original, encode_work),
        Leopard_InvalidSize, "legacy encode 64-byte restriction");
    require_legacy_result(leo_encode(64, 3, 4, 0, NULL, NULL),
        Leopard_InvalidCounts, "legacy encode R > K restriction");
    require_legacy_result(leo_encode(64, 3, 2, 4, NULL, encode_work),
        Leopard_InvalidInput, "legacy encode null input");
    require_legacy_result(leo_encode(64, 3, 2, 3, original, encode_work),
        Leopard_InvalidCounts, "legacy encode wrong work count");
    const uintptr_t overflowing_pointer_array_address =
        UINTPTR_MAX & ~(static_cast<uintptr_t>(alignof(void*)) - 1);
    const void* const* overflowing_const_pointer_array =
        reinterpret_cast<const void* const*>(
            overflowing_pointer_array_address);
    void** overflowing_mutable_pointer_array = reinterpret_cast<void**>(
        overflowing_pointer_array_address);
    require_legacy_result(leo_encode(64, 3, 2, 4,
        overflowing_const_pointer_array, encode_work),
        Leopard_InvalidInput,
        "legacy encode overflowing original pointer-array span");
    require_legacy_result(leo_encode(64, 3, 2, 4, original,
        overflowing_mutable_pointer_array), Leopard_InvalidInput,
        "legacy encode overflowing work pointer-array span");
    require_legacy_result(leo_encode(64, 65536, 1, 0, NULL, NULL),
        Leopard_TooMuchData, "legacy encode transmitted-count limit");
    require_legacy_result(leo_encode(64, 32769, 32767, 0, NULL, NULL),
        Leopard_TooMuchData, "legacy encode parent-count limit");
    require_legacy_result(leo_encode(64,
        std::numeric_limits<unsigned>::max(), 2, 0, NULL, NULL),
        Leopard_TooMuchData, "legacy encode overflow counts");
    counts->legacy_checks += 12;

    require(leo_decode_work_count(3, 2) == 8,
        "legacy decode work-count changed");
    require(leo_decode_work_count(1, 1) == 1 &&
            leo_decode_work_count(65535, 1) == 65535 &&
            leo_decode_work_count(32768, 32768) == 65536,
        "legacy decode work-count rejected a valid boundary");
    require(leo_decode_work_count(0, 1) == 0 &&
            leo_decode_work_count(2, 3) == 0 &&
            leo_decode_work_count(65536, 1) == 0 &&
            leo_decode_work_count(32769, 32767) == 0 &&
            leo_decode_work_count(std::numeric_limits<unsigned>::max(), 2) == 0,
        "legacy decode work-count accepted invalid counts");
    const void* recovery[2] = { shard, NULL };
    void* decode_work[8] = {
        shard, shard, shard, shard, shard, shard, shard, shard
    };
    require_legacy_result(leo_decode(17, 3, 2, 8, original, recovery,
        decode_work), Leopard_InvalidSize,
        "legacy decode 64-byte restriction");
    require_legacy_result(leo_decode(64, 3, 4, 0, NULL, NULL, NULL),
        Leopard_InvalidCounts, "legacy decode R > K restriction");
    require_legacy_result(leo_decode(64, 3, 2, 8, NULL, recovery,
        decode_work), Leopard_InvalidInput, "legacy decode null input");
    require_legacy_result(leo_decode(64, 3, 2, 8,
        overflowing_const_pointer_array, recovery, decode_work),
        Leopard_InvalidInput,
        "legacy decode overflowing original pointer-array span");
    require_legacy_result(leo_decode(64, 3, 2, 8, original,
        overflowing_const_pointer_array, decode_work),
        Leopard_InvalidInput,
        "legacy decode overflowing recovery pointer-array span");
    require_legacy_result(leo_decode(64, 3, 2, 8, original, recovery,
        overflowing_mutable_pointer_array), Leopard_InvalidInput,
        "legacy decode overflowing work pointer-array span");
    const void* incomplete_original[3] = { NULL, NULL, shard };
    require_legacy_result(leo_decode(64, 3, 2, 8, incomplete_original,
        recovery, decode_work), Leopard_NeedMoreData,
        "legacy decode insufficient recovery");
    require_legacy_result(leo_decode(64, 65536, 1, 0, NULL, NULL, NULL),
        Leopard_TooMuchData, "legacy decode transmitted-count limit");
    require_legacy_result(leo_decode(64, 32769, 32767, 0, NULL, NULL, NULL),
        Leopard_TooMuchData, "legacy decode parent-count limit");
    require_legacy_result(leo_decode(64,
        std::numeric_limits<unsigned>::max(), 2, 0, NULL, NULL, NULL),
        Leopard_TooMuchData, "legacy decode overflow counts");

    uint8_t single_source[64];
    uint8_t single_restored[64];
    for (size_t i = 0; i < sizeof(single_source); ++i)
        single_source[i] = static_cast<uint8_t>(i * 17u + 3u);
    memset(single_restored, 0xa5, sizeof(single_restored));
    const void* single_original[1] = { single_source };
    const void* single_recovery[1] = { NULL };
    void* single_work[1] = { single_restored };
    require_legacy_result(leo_decode(64, 1, 1, 1, single_original,
        single_recovery, single_work), Leopard_Success,
        "legacy one-original no-loss decode");
    require(memcmp(single_source, single_restored, sizeof(single_source)) == 0,
        "legacy one-original no-loss decode changed data");

    memset(single_restored, 0xa5, sizeof(single_restored));
    const void* missing_single_original[1] = { NULL };
    const void* received_single_recovery[1] = { single_source };
    require_legacy_result(leo_decode(64, 1, 1, 1,
        missing_single_original, received_single_recovery, single_work),
        Leopard_Success, "legacy one-original recovery decode");
    require(memcmp(single_source, single_restored, sizeof(single_source)) == 0,
        "legacy one-original recovery decode changed data");
    require_legacy_result(leo_decode(64, 1, 1, 1,
        missing_single_original, single_recovery, single_work),
        Leopard_NeedMoreData, "legacy one-original insufficient recovery");
    require_legacy_result(leo_decode(64, 1, 1, 0, single_original,
        single_recovery, single_work), Leopard_InvalidCounts,
        "legacy one-original wrong work count");
    require_legacy_result(leo_encode(64, 1, 1, 0, single_original,
        single_work), Leopard_InvalidCounts,
        "legacy one-original encode wrong work count");

    const void* one_original[1] = { shard };
    void* one_work[1] = { shard };
    require_legacy_result(leo_encode(64, 65535, 32768, 65536,
        one_original, one_work), Leopard_TooMuchData,
        "legacy encode field-order limit");
    counts->legacy_checks += 24;
#if SIZE_MAX < UINT64_MAX
    counts->legacy_checks += 2;
#endif
}

} // namespace

int main()
{
    try
    {
        Counts counts;
        test_result_strings(&counts);

        leo2_context_options options;
        memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AUTO;
        options.thread_count = 4;
        ContextOwner context;
        require_result(leo2_context_create(&options, &context.context),
            LEO2_SUCCESS, "contract context create");

        test_introspection_and_null_contracts(context.context, &counts);
        test_default_affinity_thread_budget(&counts);
        test_alias_and_scratch_contracts(context.context, &counts);
        test_one_shot_presence_alias_contracts(context.context, &counts);
        test_aligned_decode_input_staging_elision(context.context, &counts);
        test_tiled_decode_workspace_slopes(context.context, &counts);
        test_ragged_decode_tail_staging_slopes(context.context, &counts);
        test_ragged_encode_tail_staging_slopes(context.context, &counts);
        test_concurrent_shared_context_batches(context.context, &counts);
        test_incremental_lazy_thread_growth(&counts);
        test_lazy_thread_start_failure(&counts);
        test_deterministic_batch_failures(context.context, &counts);
        test_legacy_negative_contract(&counts);

        std::cout << "leopard2 public API contract passed:"
                  << " result_strings=" << counts.result_strings
                  << " introspection_checks=" << counts.introspection_checks
                  << " alias_checks=" << counts.alias_checks
                  << " scratch_checks=" << counts.scratch_checks
                  << " scheduler_checks=" << counts.scheduler_checks
                  << " concurrent_batches=" << counts.concurrent_batches
                  << " batch_failure_checks=" << counts.batch_failure_checks
                  << " legacy_checks=" << counts.legacy_checks
                  << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "leopard2 public API contract failed: "
                  << error.what() << std::endl;
        return 1;
    }
}
