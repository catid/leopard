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

#if defined(_MSC_VER)
#include <malloc.h>
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
        , concurrent_batches(0)
        , batch_failure_checks(0)
        , legacy_checks(0)
    {}

    uint64_t result_strings;
    uint64_t introspection_checks;
    uint64_t alias_checks;
    uint64_t scratch_checks;
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
            live_backend <= LEO2_BACKEND_NEON,
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
    require_result(leo2_encode(NULL, 17, NULL, NULL, NULL, 0),
        LEO2_INVALID_ARGUMENT, "null codec encode");
    require_result(leo2_encode_batch(NULL, NULL, 0),
        LEO2_INVALID_ARGUMENT, "null codec empty batch");
    require_result(leo2_decode_plan_execute(NULL, 17, NULL, NULL, NULL,
        NULL, 0), LEO2_INVALID_ARGUMENT, "null plan execute");
    require_result(leo2_decode_plan_execute_batch(NULL, NULL, 0),
        LEO2_INVALID_ARGUMENT, "null plan empty batch");
    counts->introspection_checks += 21;
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

    require_result(leo2_decode_plan_execute(plan.plan, bytes, decode_original,
        decode_recovery, restored_original, decode_scratch.data(),
        decode_scratch_bytes), LEO2_SUCCESS, "valid decode after rejections");
    require(restored == source[0], "valid decode restored wrong bytes");
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

void run_encode_batch_call(const BatchFixture& fixture, unsigned call_index)
{
    const size_t item_count = 4;
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
    counts->concurrent_batches +=
        static_cast<uint64_t>(caller_count) * repetitions;
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
    require(LEO_VERSION == 2, "unexpected legacy API version");
    require(leo_init() == Leopard_Success, "legacy reinitialization failed");
    ++counts->legacy_checks;

    uint8_t shard[64] = { 0 };
    const void* original[3] = { shard, shard, shard };
    void* encode_work[4] = { shard, shard, shard, shard };
    require(leo_encode_work_count(3, 2) == 4,
        "legacy encode work-count changed");
    require_legacy_result(leo_encode(17, 3, 2, 4, original, encode_work),
        Leopard_InvalidSize, "legacy encode 64-byte restriction");
    require_legacy_result(leo_encode(64, 3, 4, 0, NULL, NULL),
        Leopard_InvalidCounts, "legacy encode R > K restriction");
    require_legacy_result(leo_encode(64, 3, 2, 4, NULL, encode_work),
        Leopard_InvalidInput, "legacy encode null input");
    require_legacy_result(leo_encode(64, 3, 2, 3, original, encode_work),
        Leopard_InvalidCounts, "legacy encode wrong work count");
    counts->legacy_checks += 5;

    require(leo_decode_work_count(3, 2) == 8,
        "legacy decode work-count changed");
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
    const void* incomplete_original[3] = { NULL, NULL, shard };
    require_legacy_result(leo_decode(64, 3, 2, 8, incomplete_original,
        recovery, decode_work), Leopard_NeedMoreData,
        "legacy decode insufficient recovery");
    const void* one_original[1] = { shard };
    void* one_work[1] = { shard };
    require_legacy_result(leo_encode(64, 65535, 32768, 65536,
        one_original, one_work), Leopard_TooMuchData,
        "legacy encode field-order limit");
    counts->legacy_checks += 6;
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
        test_alias_and_scratch_contracts(context.context, &counts);
        test_aligned_decode_input_staging_elision(context.context, &counts);
        test_tiled_decode_workspace_slopes(context.context, &counts);
        test_concurrent_shared_context_batches(context.context, &counts);
        test_deterministic_batch_failures(context.context, &counts);
        test_legacy_negative_contract(&counts);

        std::cout << "leopard2 public API contract passed:"
                  << " result_strings=" << counts.result_strings
                  << " introspection_checks=" << counts.introspection_checks
                  << " alias_checks=" << counts.alias_checks
                  << " scratch_checks=" << counts.scratch_checks
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
