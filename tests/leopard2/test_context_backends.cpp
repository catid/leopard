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

#include "Leopard2Backend.h"
#include "Leopard2Direct.h"
#include "LeopardFF16.h"
#include "LeopardFF8.h"
#include "leopard.h"
#include "leopard2.h"

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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

typedef std::vector<uint8_t> Bytes;
typedef std::vector<Bytes> Shards;

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void require_result(
    leo2_result actual,
    leo2_result expected,
    const std::string& operation)
{
    if (actual != expected)
        throw std::runtime_error(operation + ": got " +
            leo2_result_string(actual) + ", expected " +
            leo2_result_string(expected));
}

class AlignedBuffer
{
public:
    explicit AlignedBuffer(size_t bytes)
        : data_(NULL)
    {
        if (bytes == 0)
            return;
#if defined(_MSC_VER)
        data_ = _aligned_malloc(bytes, leo2_scratch_alignment());
#else
        if (posix_memalign(&data_, leo2_scratch_alignment(), bytes) != 0)
            data_ = NULL;
#endif
        if (!data_)
            throw std::bad_alloc();
        std::memset(data_, 0, bytes);
    }

    ~AlignedBuffer()
    {
#if defined(_MSC_VER)
        _aligned_free(data_);
#else
        std::free(data_);
#endif
    }

    void* data() { return data_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);
    void* data_;
};

struct ContextEntry
{
    const leopard::backend::Ops* table;
    leo2_context* context;
};

struct TraceState
{
    std::atomic<const leopard::backend::Ops*> delegate;
    std::atomic<uint64_t> ff8_calls;
    std::atomic<uint64_t> ff16_calls;
    std::atomic<uint64_t> xor_calls;
    std::atomic<uint64_t> xor_two_to_one_calls;
    std::atomic<uint64_t> ff8_ifft_four_calls;
    std::atomic<uint64_t> ff8_fft_four_calls;
    std::atomic<uint64_t> ff16_ifft_four_calls;
    std::atomic<uint64_t> ff16_fft_four_calls;
    std::atomic<uint64_t> xor_four_calls;
};

TraceState g_trace;

const leopard::backend::Ops* trace_delegate()
{
    const leopard::backend::Ops* ops =
        g_trace.delegate.load(std::memory_order_acquire);
    if (!ops)
        std::abort();
    return ops;
}

void trace_ff8_multiply(void* destination, const void* source,
    uint16_t log, uint64_t bytes)
{
    g_trace.ff8_calls.fetch_add(1, std::memory_order_relaxed);
    trace_delegate()->ff8_multiply(destination, source, log, bytes);
}

void trace_ff8_multiply_add(void* destination, const void* source,
    uint16_t log, uint64_t bytes)
{
    g_trace.ff8_calls.fetch_add(1, std::memory_order_relaxed);
    trace_delegate()->ff8_multiply_add(destination, source, log, bytes);
}

void trace_ff16_multiply(void* destination, const void* source,
    uint16_t log, uint64_t bytes)
{
    g_trace.ff16_calls.fetch_add(1, std::memory_order_relaxed);
    trace_delegate()->ff16_multiply(destination, source, log, bytes);
}

void trace_ff16_multiply_add(void* destination, const void* source,
    uint16_t log, uint64_t bytes)
{
    g_trace.ff16_calls.fetch_add(1, std::memory_order_relaxed);
    trace_delegate()->ff16_multiply_add(destination, source, log, bytes);
}

void trace_xor(void* destination, const void* source, uint64_t bytes)
{
    g_trace.xor_calls.fetch_add(1, std::memory_order_relaxed);
    trace_delegate()->xor_memory(destination, source, bytes);
}

void trace_xor_2to1(
    void* destination,
    const void* source0,
    const void* source1,
    uint64_t bytes)
{
    g_trace.xor_two_to_one_calls.fetch_add(1, std::memory_order_relaxed);
    trace_delegate()->xor_memory_2to1(
        destination, source0, source1, bytes);
}

void trace_xor4(
    void* destination0, const void* source0,
    void* destination1, const void* source1,
    void* destination2, const void* source2,
    void* destination3, const void* source3,
    uint64_t bytes)
{
    g_trace.xor_four_calls.fetch_add(1, std::memory_order_relaxed);
    trace_delegate()->xor_memory4(
        destination0, source0, destination1, source1,
        destination2, source2, destination3, source3, bytes);
}

void trace_ff8_ifft(void* x, void* y, uint16_t log, uint64_t bytes)
{
    g_trace.ff8_calls.fetch_add(1, std::memory_order_relaxed);
    trace_delegate()->ff8_ifft_butterfly2(x, y, log, bytes);
}

void trace_ff8_fft(void* x, void* y, uint16_t log, uint64_t bytes)
{
    g_trace.ff8_calls.fetch_add(1, std::memory_order_relaxed);
    trace_delegate()->ff8_fft_butterfly2(x, y, log, bytes);
}

void trace_ff8_ifft_xor(const void* x, const void* y, void* x_output,
    void* y_output, uint16_t log, uint64_t bytes)
{
    g_trace.ff8_calls.fetch_add(1, std::memory_order_relaxed);
    trace_delegate()->ff8_ifft_butterfly2_xor(
        x, y, x_output, y_output, log, bytes);
}

void trace_ff16_ifft(void* x, void* y, uint16_t log, uint64_t bytes)
{
    g_trace.ff16_calls.fetch_add(1, std::memory_order_relaxed);
    trace_delegate()->ff16_ifft_butterfly2(x, y, log, bytes);
}

void trace_ff16_fft(void* x, void* y, uint16_t log, uint64_t bytes)
{
    g_trace.ff16_calls.fetch_add(1, std::memory_order_relaxed);
    trace_delegate()->ff16_fft_butterfly2(x, y, log, bytes);
}

void trace_ff8_ifft4(
    void* value0, void* value1, void* value2, void* value3,
    uint16_t log01, uint16_t log23, uint16_t log02, uint64_t bytes)
{
    g_trace.ff8_calls.fetch_add(1, std::memory_order_relaxed);
    g_trace.ff8_ifft_four_calls.fetch_add(1, std::memory_order_relaxed);
    trace_delegate()->ff8_ifft_butterfly4(
        value0, value1, value2, value3,
        log01, log23, log02, bytes);
}

void trace_ff8_fft4(
    void* value0, void* value1, void* value2, void* value3,
    uint16_t log01, uint16_t log23, uint16_t log02, uint64_t bytes)
{
    g_trace.ff8_calls.fetch_add(1, std::memory_order_relaxed);
    g_trace.ff8_fft_four_calls.fetch_add(1, std::memory_order_relaxed);
    trace_delegate()->ff8_fft_butterfly4(
        value0, value1, value2, value3,
        log01, log23, log02, bytes);
}

void trace_ff16_ifft4(
    void* value0, void* value1, void* value2, void* value3,
    uint16_t log01, uint16_t log23, uint16_t log02, uint64_t bytes)
{
    g_trace.ff16_calls.fetch_add(1, std::memory_order_relaxed);
    g_trace.ff16_ifft_four_calls.fetch_add(1, std::memory_order_relaxed);
    trace_delegate()->ff16_ifft_butterfly4(
        value0, value1, value2, value3,
        log01, log23, log02, bytes);
}

void trace_ff16_fft4(
    void* value0, void* value1, void* value2, void* value3,
    uint16_t log01, uint16_t log23, uint16_t log02, uint64_t bytes)
{
    g_trace.ff16_calls.fetch_add(1, std::memory_order_relaxed);
    g_trace.ff16_fft_four_calls.fetch_add(1, std::memory_order_relaxed);
    trace_delegate()->ff16_fft_butterfly4(
        value0, value1, value2, value3,
        log01, log23, log02, bytes);
}

class TraceOpsGuard
{
public:
    explicit TraceOpsGuard(const ContextEntry& entry)
        : entry_(entry), tracing_(*entry.table)
    {
        tracing_.ff8_multiply = trace_ff8_multiply;
        tracing_.ff8_multiply_add = trace_ff8_multiply_add;
        tracing_.ff16_multiply = trace_ff16_multiply;
        tracing_.ff16_multiply_add = trace_ff16_multiply_add;
        tracing_.xor_memory = trace_xor;
        tracing_.xor_memory_2to1 = trace_xor_2to1;
        tracing_.xor_memory4 = trace_xor4;
        tracing_.ff8_ifft_butterfly2 = trace_ff8_ifft;
        tracing_.ff8_fft_butterfly2 = trace_ff8_fft;
        tracing_.ff8_ifft_butterfly2_xor = trace_ff8_ifft_xor;
        tracing_.ff8_ifft_butterfly4 = trace_ff8_ifft4;
        tracing_.ff8_fft_butterfly4 = trace_ff8_fft4;
        tracing_.ff16_ifft_butterfly2 = trace_ff16_ifft;
        tracing_.ff16_fft_butterfly2 = trace_ff16_fft;
        tracing_.ff16_ifft_butterfly4 = trace_ff16_ifft4;
        tracing_.ff16_fft_butterfly4 = trace_ff16_fft4;
        g_trace.delegate.store(entry.table, std::memory_order_release);
        reset();
        leopard::backend::TestSetContextOps(entry.context, &tracing_);
    }

    ~TraceOpsGuard()
    {
        leopard::backend::TestSetContextOps(entry_.context, entry_.table);
        g_trace.delegate.store(NULL, std::memory_order_release);
    }

    void reset()
    {
        g_trace.ff8_calls.store(0, std::memory_order_relaxed);
        g_trace.ff16_calls.store(0, std::memory_order_relaxed);
        g_trace.xor_calls.store(0, std::memory_order_relaxed);
        g_trace.xor_two_to_one_calls.store(0, std::memory_order_relaxed);
        g_trace.ff8_ifft_four_calls.store(0, std::memory_order_relaxed);
        g_trace.ff8_fft_four_calls.store(0, std::memory_order_relaxed);
        g_trace.ff16_ifft_four_calls.store(0, std::memory_order_relaxed);
        g_trace.ff16_fft_four_calls.store(0, std::memory_order_relaxed);
        g_trace.xor_four_calls.store(0, std::memory_order_relaxed);
        leopard::ff8::TestOnlyResetTransformCallsiteCounts();
        leopard::ff16::TestOnlyResetTransformCallsiteCounts();
    }

    uint64_t ff8_calls() const
    {
        return g_trace.ff8_calls.load(std::memory_order_relaxed);
    }
    uint64_t ff16_calls() const
    {
        return g_trace.ff16_calls.load(std::memory_order_relaxed);
    }
    uint64_t xor_calls() const
    {
        return g_trace.xor_calls.load(std::memory_order_relaxed);
    }
    uint64_t xor_two_to_one_calls() const
    {
        return g_trace.xor_two_to_one_calls.load(
            std::memory_order_relaxed);
    }
    uint64_t ff8_ifft_four_calls() const
    {
        return g_trace.ff8_ifft_four_calls.load(std::memory_order_relaxed);
    }
    uint64_t ff8_fft_four_calls() const
    {
        return g_trace.ff8_fft_four_calls.load(std::memory_order_relaxed);
    }
    uint64_t ff16_ifft_four_calls() const
    {
        return g_trace.ff16_ifft_four_calls.load(std::memory_order_relaxed);
    }
    uint64_t ff16_fft_four_calls() const
    {
        return g_trace.ff16_fft_four_calls.load(std::memory_order_relaxed);
    }
    uint64_t xor_four_calls() const
    {
        return g_trace.xor_four_calls.load(std::memory_order_relaxed);
    }

private:
    TraceOpsGuard(const TraceOpsGuard&);
    TraceOpsGuard& operator=(const TraceOpsGuard&);
    const ContextEntry& entry_;
    leopard::backend::Ops tracing_;
};

void require_four_way_callsites(
    const TraceOpsGuard& trace,
    leo2_backend execution_backend,
    leo2_field field,
    size_t public_bytes,
    const std::string& operation,
    bool expect_ff8_accumulating_ifft,
    bool split_ragged_execution)
{
    if (field == LEO2_FIELD_GF8)
    {
        const leopard::ff8::TestOnlyTransformCallsiteCounts callsites =
            leopard::ff8::TestOnlyGetTransformCallsiteCounts();
        require(callsites.ifft_dit4 != 0,
            operation + " did not exercise GF8 IFFT_DIT4");
        if (expect_ff8_accumulating_ifft)
            require(callsites.ifft_dit4_xor != 0,
                operation + " did not exercise GF8 IFFT_DIT4_xor");
        else
            require(callsites.ifft_dit4_xor == 0,
                operation + " unexpectedly exercised GF8 IFFT_DIT4_xor");
        require(callsites.fft_dit4 != 0,
            operation + " did not exercise GF8 FFT_DIT4");
        require(trace.ff8_ifft_four_calls() ==
                callsites.ifft_dit4 + callsites.ifft_dit4_xor,
            operation + " has a GF8 inverse radix-four callsite bypass");
        require(trace.ff8_fft_four_calls() == callsites.fft_dit4,
            operation + " has a GF8 forward radix-four callsite bypass");
    }
    else
    {
        const leopard::ff16::TestOnlyTransformCallsiteCounts callsites =
            leopard::ff16::TestOnlyGetTransformCallsiteCounts();
        require(callsites.ifft_dit4 != 0,
            operation + " did not exercise GF16 IFFT_DIT4");
        require(callsites.fft_dit4 != 0,
            operation + " did not exercise GF16 FFT_DIT4");
        require(callsites.ifft_dit4 ==
                callsites.ifft_dit4_fused + callsites.ifft_dit4_split,
            operation + " has an unclassified GF16 inverse radix-four callsite");
        require(callsites.fft_dit4 ==
                callsites.fft_dit4_fused + callsites.fft_dit4_split,
            operation + " has an unclassified GF16 forward radix-four callsite");
        require(trace.ff16_ifft_four_calls() == callsites.ifft_dit4_fused,
            operation + " has a GF16 inverse fused-dispatch mismatch");
        require(trace.ff16_fft_four_calls() == callsites.fft_dit4_fused,
            operation + " has a GF16 forward fused-dispatch mismatch");

        const size_t transform_bytes = (public_bytes + 63U) & ~size_t(63U);
        const size_t prefix_bytes = public_bytes & ~size_t(63U);
        const bool transform_fused = transform_bytes == 64U ||
            (transform_bytes == 128U &&
             execution_backend == LEO2_BACKEND_AVX2);
        const bool prefix_fused = prefix_bytes == 64U ||
            (prefix_bytes == 128U &&
             execution_backend == LEO2_BACKEND_AVX2);
        const bool has_split_tail = split_ragged_execution &&
            (public_bytes & 63U) != 0;
        const bool mixed_split = has_split_tail && prefix_bytes != 0 &&
            !prefix_fused;
        const bool expect_all_fused = has_split_tail
            ? prefix_bytes == 0 || prefix_fused
            : transform_fused;
        if (mixed_split)
        {
            require(callsites.ifft_dit4_fused != 0 &&
                    callsites.ifft_dit4_split != 0,
                operation + " did not classify split-prefix/tiled-tail inverse calls");
            require(callsites.fft_dit4_fused != 0 &&
                    callsites.fft_dit4_split != 0,
                operation + " did not classify split-prefix/tiled-tail forward calls");
        }
        else if (expect_all_fused)
        {
            require(callsites.ifft_dit4_fused == callsites.ifft_dit4 &&
                    callsites.ifft_dit4_split == 0,
                operation + " did not fuse the qualified GF16 inverse size");
            require(callsites.fft_dit4_fused == callsites.fft_dit4 &&
                    callsites.fft_dit4_split == 0,
                operation + " did not fuse the qualified GF16 forward size");
        }
        else
        {
            require(callsites.ifft_dit4_split == callsites.ifft_dit4 &&
                    callsites.ifft_dit4_fused == 0,
                operation + " fused an unqualified GF16 inverse size");
            require(callsites.fft_dit4_split == callsites.fft_dit4 &&
                    callsites.fft_dit4_fused == 0,
                operation + " fused an unqualified GF16 forward size");
        }
    }
}

struct CodecCase
{
    uint32_t k;
    uint32_t r;
    leo2_profile profile;
    leo2_field field;
    size_t bytes;
};

std::vector<const void*> const_pointers(const Shards& shards)
{
    std::vector<const void*> pointers(shards.size(), NULL);
    for (size_t i = 0; i < shards.size(); ++i)
        pointers[i] = shards[i].empty() ? NULL : &shards[i][0];
    return pointers;
}

std::vector<void*> mutable_pointers(Shards& shards)
{
    std::vector<void*> pointers(shards.size(), NULL);
    for (size_t i = 0; i < shards.size(); ++i)
        pointers[i] = shards[i].empty() ? NULL : &shards[i][0];
    return pointers;
}

Shards make_originals(const CodecCase& test_case, uint32_t seed)
{
    Shards originals(test_case.k, Bytes(test_case.bytes));
    uint32_t state = seed;
    for (uint32_t shard = 0; shard < test_case.k; ++shard)
        for (size_t i = 0; i < test_case.bytes; ++i)
        {
            state = state * 1664525U + 1013904223U;
            originals[shard][i] = static_cast<uint8_t>(
                (state >> 24) ^ (shard * 31U + i * 17U));
        }
    return originals;
}

Shards execute_encode(
    const leo2_codec* codec,
    const CodecCase& test_case,
    const Shards& originals)
{
    size_t scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(codec, test_case.bytes,
        &scratch_bytes), LEO2_SUCCESS, "encode scratch query");
    AlignedBuffer scratch(scratch_bytes);
    Shards recovery(test_case.r, Bytes(test_case.bytes));
    const std::vector<const void*> original_pointers =
        const_pointers(originals);
    std::vector<void*> recovery_pointers = mutable_pointers(recovery);
    require_result(leo2_encode(codec, test_case.bytes,
        &original_pointers[0], &recovery_pointers[0], scratch.data(),
        scratch_bytes), LEO2_SUCCESS, "encode");
    return recovery;
}

Shards encode_case(
    leo2_context* context,
    const CodecCase& test_case,
    const Shards& originals,
    leo2_codec** codec_out)
{
    leo2_codec* codec = NULL;
    require_result(leo2_codec_create(context, test_case.k, test_case.r,
        test_case.profile, test_case.field, NULL, &codec), LEO2_SUCCESS,
        "codec create");
    require(codec != NULL, "codec create returned null");
    require(leo2_codec_profile(codec) == test_case.profile,
        "profile identity changed across backends");
    require(leo2_codec_field(codec) == test_case.field,
        "field identity changed across backends");
    *codec_out = codec;
    return execute_encode(codec, test_case, originals);
}

void execute_plan(
    const leo2_decode_plan* plan,
    const CodecCase& test_case,
    const Shards& originals,
    const Shards& recovery)
{
    const uint32_t first_missing = 0;
    const uint32_t second_missing = test_case.k > 2 ? test_case.k - 1 : 1;
    size_t scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(plan, test_case.bytes,
        &scratch_bytes), LEO2_SUCCESS, "decode scratch query");
    AlignedBuffer scratch(scratch_bytes);
    std::vector<const void*> original_pointers = const_pointers(originals);
    const std::vector<const void*> recovery_pointers =
        const_pointers(recovery);
    original_pointers[first_missing] = NULL;
    original_pointers[second_missing] = NULL;
    Shards restored(test_case.k, Bytes(test_case.bytes));
    std::vector<void*> restored_pointers(test_case.k, NULL);
    restored_pointers[first_missing] = &restored[first_missing][0];
    restored_pointers[second_missing] = &restored[second_missing][0];
    require_result(leo2_decode_plan_execute(plan, test_case.bytes,
        &original_pointers[0], &recovery_pointers[0], &restored_pointers[0],
        scratch.data(), scratch_bytes), LEO2_SUCCESS, "decode execute");
    require(restored[first_missing] == originals[first_missing],
        "first restored original differs across backend");
    require(restored[second_missing] == originals[second_missing],
        "second restored original differs across backend");
}

leo2_decode_plan* make_plan(
    const leo2_codec* codec,
    const CodecCase& test_case)
{
    std::vector<uint8_t> original_present(test_case.k, 1);
    std::vector<uint8_t> recovery_present(test_case.r, 1);
    original_present[0] = 0;
    original_present[test_case.k > 2 ? test_case.k - 1 : 1] = 0;
    leo2_decode_plan* plan = NULL;
    require_result(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &plan), LEO2_SUCCESS, "decode plan create");
    require(plan != NULL, "decode plan create returned null");
    return plan;
}

void decode_case(
    const leo2_codec* codec,
    const CodecCase& test_case,
    const Shards& originals,
    const Shards& recovery)
{
    leo2_decode_plan* plan = make_plan(codec, test_case);
    execute_plan(plan, test_case, originals, recovery);
    leo2_decode_plan_destroy(plan);
}

void test_concurrent_first_use()
{
    static const unsigned request_count = 3;
    static const unsigned lanes_per_request = 8;
    const leo2_backend requests[request_count] = {
        LEO2_BACKEND_SCALAR,
        LEO2_BACKEND_SSSE3,
        LEO2_BACKEND_AVX2
    };
    leo2_result results[request_count][lanes_per_request];
    for (unsigned request_i = 0; request_i < request_count; ++request_i)
        for (unsigned lane = 0; lane < lanes_per_request; ++lane)
            results[request_i][lane] = LEO2_INTERNAL_ERROR;

    const CodecCase test_case = {
        3, 2, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 17
    };
    std::atomic<unsigned> ready(0);
    std::atomic<bool> go(false);
    std::atomic<unsigned> failures(0);
    std::vector<std::thread> threads;
    for (unsigned request_i = 0; request_i < request_count; ++request_i)
        for (unsigned lane = 0; lane < lanes_per_request; ++lane)
        {
            threads.push_back(std::thread(
                [request_i, lane, &requests, &results, &test_case, &ready,
                 &go, &failures]() {
                leo2_context* context = NULL;
                leo2_codec* codec = NULL;
                try
                {
                    ready.fetch_add(1, std::memory_order_release);
                    while (!go.load(std::memory_order_acquire))
                        std::this_thread::yield();

                    leo2_context_options options;
                    std::memset(&options, 0, sizeof(options));
                    options.struct_size = sizeof(options);
                    options.backend = requests[request_i];
                    options.thread_count = 1;
                    const leo2_result result = leo2_context_create(
                        &options, &context);
                    results[request_i][lane] = result;
                    require(result == LEO2_SUCCESS ||
                            result == LEO2_UNSUPPORTED,
                        "first-use context returned an unexpected result");
                    if (result == LEO2_SUCCESS)
                    {
                        require(context != NULL,
                            "first-use context returned null");
                        require(leo2_context_backend(context) ==
                                requests[request_i],
                            "first-use context reported wrong backend");
                        const Shards originals = make_originals(test_case,
                            static_cast<uint32_t>(0x510e527fU +
                                request_i * 131U + lane));
                        const Shards recovery = encode_case(context,
                            test_case, originals, &codec);
                        decode_case(codec, test_case, originals, recovery);
                    }
                    else
                        require(context == NULL,
                            "unsupported first-use context was not cleared");
                }
                catch (...)
                {
                    failures.fetch_add(1, std::memory_order_relaxed);
                }
                leo2_codec_destroy(codec);
                leo2_context_destroy(context);
            }));
        }

    while (ready.load(std::memory_order_acquire) != threads.size())
        std::this_thread::yield();
    go.store(true, std::memory_order_release);
    for (size_t i = 0; i < threads.size(); ++i)
        threads[i].join();
    require(failures.load(std::memory_order_relaxed) == 0,
        "concurrent first-use context execution failed");
    for (unsigned request_i = 0; request_i < request_count; ++request_i)
        for (unsigned lane = 1; lane < lanes_per_request; ++lane)
            require(results[request_i][lane] == results[request_i][0],
                "first-use qualification result was nondeterministic");
}

std::vector<ContextEntry> create_contexts()
{
    std::vector<ContextEntry> contexts;
    const leo2_backend runtime = leopard::backend::ExecutionBackend();
    const leo2_backend requests[] = {
        LEO2_BACKEND_AUTO,
        LEO2_BACKEND_SCALAR,
        LEO2_BACKEND_SSSE3,
        LEO2_BACKEND_AVX2,
        LEO2_BACKEND_NEON
    };
    for (size_t i = 0; i < sizeof(requests) / sizeof(requests[0]); ++i)
    {
        leo2_context_options options;
        std::memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = requests[i];
        options.thread_count = 4;
        leo2_context* context = reinterpret_cast<leo2_context*>(
            static_cast<uintptr_t>(1));
        const leo2_result result = leo2_context_create(&options, &context);
        bool expected = requests[i] == LEO2_BACKEND_AUTO;
        if (requests[i] == LEO2_BACKEND_NEON)
            expected = runtime == LEO2_BACKEND_NEON;
        else if (requests[i] != LEO2_BACKEND_AUTO &&
                 runtime != LEO2_BACKEND_NEON)
            expected = leopard::backend::GetQualifiedOps(requests[i]) != NULL;
        require_result(result, expected ? LEO2_SUCCESS : LEO2_UNSUPPORTED,
            "explicit context result");
        if (!expected)
        {
            require(context == NULL, "failed context did not clear output");
            continue;
        }
        const leo2_backend effective = leo2_context_backend(context);
        require(effective == (requests[i] == LEO2_BACKEND_AUTO
                ? runtime : requests[i]),
            "context reports the wrong effective backend");
        const leopard::backend::Ops* table = effective == LEO2_BACKEND_NEON
            ? leopard::backend::GetQualifiedOps(LEO2_BACKEND_AUTO)
            : leopard::backend::GetQualifiedOps(effective);
        require(table != NULL, "context execution table is missing");
        const ContextEntry entry = { table, context };
        contexts.push_back(entry);
    }

    leo2_context_options invalid;
    std::memset(&invalid, 0, sizeof(invalid));
    invalid.struct_size = sizeof(invalid);
    invalid.backend = static_cast<uint32_t>(LEO2_BACKEND_NEON) + 1U;
    leo2_context* rejected = reinterpret_cast<leo2_context*>(
        static_cast<uintptr_t>(1));
    require_result(leo2_context_create(&invalid, &rejected),
        LEO2_INVALID_ARGUMENT, "out-of-range backend request");
    require(rejected == NULL, "invalid backend did not clear output");
    return contexts;
}

void test_process_default_immutable(
    const leopard::backend::Ops* process_default)
{
    require(process_default != NULL, "process default is missing");
    require(&leopard::backend::GetOps() == process_default,
        "legacy accessor differs from process default");
    for (unsigned raw = static_cast<unsigned>(LEO2_BACKEND_SCALAR);
         raw <= static_cast<unsigned>(LEO2_BACKEND_AVX2); ++raw)
    {
        (void)leopard::backend::GetQualifiedOps(
            static_cast<leo2_backend>(raw));
        require(&leopard::backend::GetDefaultOps() == process_default,
            "explicit qualification changed the process default");
    }
}

void test_traced_context_dispatch(const std::vector<ContextEntry>& contexts)
{
    const CodecCase transform_cases[] = {
        { 33, 16, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 257 },
        { 17, 33, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8, 129 },
        { 33, 17, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, 64 },
        { 17, 33, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16, 64 },
        { 33, 17, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, 66 },
        { 17, 33, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16, 66 },
        { 33, 17, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, 128 },
        { 17, 33, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16, 128 },
        { 33, 17, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, 130 },
        { 17, 33, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16, 130 },
        { 33, 17, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, 1026 },
        { 17, 33, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16, 1026 }
    };

    for (size_t context_i = 0; context_i < contexts.size(); ++context_i)
    {
        // Existing ARM NEON/SSE2NEON transforms intentionally precede the
        // scalar fallback table.  Their parity/recovery/concurrency coverage
        // remains below; exact table tracing becomes authoritative when the
        // native NEON ops extraction lands.
        if (leo2_context_backend(contexts[context_i].context) ==
                LEO2_BACKEND_NEON)
            continue;
        TraceOpsGuard trace(contexts[context_i]);
        for (size_t case_i = 0;
             case_i < sizeof(transform_cases) / sizeof(transform_cases[0]);
             ++case_i)
        {
            const CodecCase& test_case = transform_cases[case_i];
            const std::string profile_name =
                test_case.profile == LEO2_PROFILE_LEGACY_HIGH_V1 ?
                    "high-profile" : "low-profile";
            const Shards originals = make_originals(test_case,
                static_cast<uint32_t>(0x243f6a88U + case_i * 977U));
            leo2_codec* codec = NULL;
            trace.reset();
            const Shards recovery = encode_case(contexts[context_i].context,
                test_case, originals, &codec);
            if (test_case.field == LEO2_FIELD_GF8)
            {
                require(trace.ff8_calls() != 0,
                    "GF8 encode bypassed the context ops table");
            }
            else
            {
                require(trace.ff16_calls() != 0,
                    "GF16 encode bypassed the context ops table");
            }
            require_four_way_callsites(trace,
                leo2_context_backend(contexts[context_i].context),
                test_case.field,
                test_case.bytes,
                profile_name + " encode",
                test_case.field == LEO2_FIELD_GF8 &&
                    test_case.profile == LEO2_PROFILE_LEGACY_HIGH_V1,
                true);

            trace.reset();
            decode_case(codec, test_case, originals, recovery);
            if (test_case.field == LEO2_FIELD_GF8)
            {
                require(trace.ff8_calls() != 0,
                    "GF8 decode bypassed the context ops table");
            }
            else
            {
                require(trace.ff16_calls() != 0,
                    "GF16 decode bypassed the context ops table");
            }
            require_four_way_callsites(trace,
                leo2_context_backend(contexts[context_i].context),
                test_case.field,
                test_case.bytes,
                profile_name + " decode", false, true);
            require(trace.xor_four_calls() != 0,
                "decode bypassed the context grouped-XOR table");
            leo2_codec_destroy(codec);
        }

        const CodecCase generic_case = {
            33, 17, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 257
        };
        const Shards generic_originals = make_originals(
            generic_case, 0xa4093822U);
        leo2_codec_options generic_options;
        std::memset(&generic_options, 0, sizeof(generic_options));
        generic_options.struct_size = sizeof(generic_options);
        generic_options.flags = LEO2_CODEC_FORCE_GENERIC_DECODE;
        leo2_codec* generic_codec = NULL;
        require_result(leo2_codec_create(contexts[context_i].context,
            generic_case.k, generic_case.r, generic_case.profile,
            generic_case.field, &generic_options, &generic_codec),
            LEO2_SUCCESS, "traced generic codec create");
        const Shards generic_recovery = execute_encode(
            generic_codec, generic_case, generic_originals);
        trace.reset();
        decode_case(generic_codec, generic_case, generic_originals,
            generic_recovery);
        require(trace.ff8_calls() != 0,
            "generic decode bypassed the context ops table");
        require_four_way_callsites(trace,
            leo2_context_backend(contexts[context_i].context),
            generic_case.field,
            generic_case.bytes,
            "generic high-profile decode", false, true);
        require(trace.xor_four_calls() != 0,
            "generic decode bypassed the context grouped-XOR table");
        leo2_codec_destroy(generic_codec);

        const leo2_field direct_fields[] = {
            LEO2_FIELD_GF8, LEO2_FIELD_GF16
        };
        for (size_t field_i = 0; field_i < 2; ++field_i)
        {
            const CodecCase test_case = {
                3, 2, LEO2_PROFILE_LOW_V1, direct_fields[field_i],
                field_i == 0 ? 33U : 34U
            };
            const Shards originals = make_originals(test_case,
                static_cast<uint32_t>(0x85a308d3U + field_i));
            leo2_codec* codec = NULL;
            require_result(leo2_codec_create(contexts[context_i].context,
                test_case.k, test_case.r, test_case.profile, test_case.field,
                NULL, &codec), LEO2_SUCCESS, "traced direct codec create");
            require_result(leo2_test_codec_set_encode_mode(codec,
                LEO2_TEST_ENCODE_FORCE_DIRECT), LEO2_SUCCESS,
                "force traced direct encode");
            trace.reset();
            const Shards recovery = execute_encode(
                codec, test_case, originals);
            if (test_case.field == LEO2_FIELD_GF8)
                require(trace.ff8_calls() != 0,
                    "direct GF8 encode bypassed the context ops table");
            else
                require(trace.ff16_calls() != 0,
                    "direct GF16 encode bypassed the context ops table");

            trace.reset();
            decode_case(codec, test_case, originals, recovery);
            if (test_case.field == LEO2_FIELD_GF8)
                require(trace.ff8_calls() != 0,
                    "direct GF8 repair bypassed the context ops table");
            else
                require(trace.ff16_calls() != 0,
                    "direct GF16 repair bypassed the context ops table");
            leo2_codec_destroy(codec);
        }

        const CodecCase xor_case = {
            9, 1, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 33
        };
        const Shards xor_originals = make_originals(xor_case, 0x13198a2eU);
        leo2_codec* xor_codec = NULL;
        trace.reset();
        const Shards xor_recovery = encode_case(contexts[context_i].context,
            xor_case, xor_originals, &xor_codec);
        require(trace.xor_two_to_one_calls() == 4,
            "R=1 encode bypassed the fused context XOR table");
        require(trace.xor_calls() == 0,
            "even R=1 encode unexpectedly used a single-source XOR tail");

        std::vector<uint8_t> original_present(xor_case.k, 1);
        std::vector<uint8_t> recovery_present(xor_case.r, 1);
        original_present[0] = 0;
        leo2_decode_plan* xor_plan = NULL;
        require_result(leo2_decode_plan_create(xor_codec,
            &original_present[0], &recovery_present[0], &xor_plan),
            LEO2_SUCCESS, "R=1 plan create");
        size_t scratch_bytes = 0;
        require_result(leo2_decode_plan_scratch_size(xor_plan, xor_case.bytes,
            &scratch_bytes), LEO2_SUCCESS, "R=1 scratch query");
        AlignedBuffer scratch(scratch_bytes);
        std::vector<const void*> original_pointers =
            const_pointers(xor_originals);
        const std::vector<const void*> recovery_pointers =
            const_pointers(xor_recovery);
        original_pointers[0] = NULL;
        Bytes restored(xor_case.bytes);
        std::vector<void*> restored_pointers(xor_case.k, NULL);
        restored_pointers[0] = &restored[0];
        trace.reset();
        require_result(leo2_decode_plan_execute(xor_plan, xor_case.bytes,
            &original_pointers[0], &recovery_pointers[0],
            &restored_pointers[0], scratch.data(), scratch_bytes),
            LEO2_SUCCESS, "R=1 decode");
        require(restored == xor_originals[0], "R=1 restored data mismatch");
        require(trace.xor_two_to_one_calls() == 4,
            "R=1 decode bypassed the fused context XOR table");
        require(trace.xor_calls() == 0,
            "even R=1 decode unexpectedly used a single-source XOR tail");
        leo2_decode_plan_destroy(xor_plan);
        leo2_codec_destroy(xor_codec);
    }
}

void test_public_codecs(const std::vector<ContextEntry>& contexts)
{
    const CodecCase cases[] = {
        { 9, 7, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 257 },
        { 5, 11, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8, 129 },
        { 193, 65, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, 1024 },
        { 129, 129, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16, 1024 }
    };
    for (size_t case_i = 0; case_i < sizeof(cases) / sizeof(cases[0]); ++case_i)
    {
        const Shards originals = make_originals(cases[case_i],
            static_cast<uint32_t>(0x9e3779b9U + case_i * 101U));
        Shards reference;
        uint32_t reference_parent = 0;
        for (size_t context_i = 0; context_i < contexts.size(); ++context_i)
        {
            leo2_codec* codec = NULL;
            const Shards recovery = encode_case(contexts[context_i].context,
                cases[case_i], originals, &codec);
            if (context_i == 0)
            {
                reference = recovery;
                reference_parent = leo2_codec_parent_count(codec);
            }
            else
            {
                require(recovery == reference,
                    "parity differs between context backends");
                require(leo2_codec_parent_count(codec) == reference_parent,
                    "parent identity differs between context backends");
            }
            decode_case(codec, cases[case_i], originals, recovery);
            leo2_codec_destroy(codec);
        }
    }
}

void test_concurrent_public_codecs(
    const std::vector<ContextEntry>& contexts)
{
    const CodecCase test_cases[] = {
        { 33, 17, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 1025 },
        { 17, 33, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8, 1025 }
    };
    std::atomic<unsigned> ready(0);
    std::atomic<bool> go(false);
    std::atomic<unsigned> failures(0);
    std::vector<std::thread> threads;
    for (size_t context_i = 0; context_i < contexts.size(); ++context_i)
    {
        for (size_t case_i = 0;
             case_i < sizeof(test_cases) / sizeof(test_cases[0]); ++case_i)
        {
            for (unsigned lane = 0; lane < 2; ++lane)
            {
                leo2_context* const context = contexts[context_i].context;
                const CodecCase test_case = test_cases[case_i];
                const uint32_t seed = static_cast<uint32_t>(
                    0x243f6a88U + context_i * 257U + case_i * 67U +
                    lane * 17U);
                threads.push_back(std::thread(
                    [context, test_case, seed, &ready, &go, &failures]() {
                    try
                    {
                        const Shards originals = make_originals(
                            test_case, seed);
                        ready.fetch_add(1, std::memory_order_release);
                        while (!go.load(std::memory_order_acquire))
                            std::this_thread::yield();
                        for (unsigned iteration = 0; iteration < 8;
                             ++iteration)
                        {
                            leo2_codec* codec = NULL;
                            const Shards recovery = encode_case(
                                context, test_case, originals, &codec);
                            decode_case(
                                codec, test_case, originals, recovery);
                            leo2_codec_destroy(codec);
                        }
                    }
                    catch (...)
                    {
                        failures.fetch_add(1, std::memory_order_relaxed);
                    }
                }));
            }
        }
    }
    while (ready.load(std::memory_order_acquire) != threads.size())
        std::this_thread::yield();
    go.store(true, std::memory_order_release);
    for (size_t i = 0; i < threads.size(); ++i)
        threads[i].join();
    require(failures.load(std::memory_order_relaxed) == 0,
        "shared/mixed-context public execution failed");
}

void test_shared_codec_and_plan(
    const std::vector<ContextEntry>& contexts)
{
    const CodecCase test_cases[] = {
        { 33, 17, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, 1024 },
        { 17, 33, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16, 1026 }
    };
    const unsigned thread_count = 4;
    const unsigned repetitions = 4;
    for (size_t case_i = 0;
         case_i < sizeof(test_cases) / sizeof(test_cases[0]); ++case_i)
    {
        const CodecCase test_case = test_cases[case_i];
        const Shards originals = make_originals(test_case,
            static_cast<uint32_t>(0xb7e15162U + case_i * 131U));
        for (size_t context_i = 0; context_i < contexts.size(); ++context_i)
        {
            leo2_codec* codec = NULL;
            const Shards reference = encode_case(contexts[context_i].context,
                test_case, originals, &codec);
            leo2_decode_plan* plan = make_plan(codec, test_case);
            std::atomic<unsigned> ready(0);
            std::atomic<bool> go(false);
            std::atomic<unsigned> failures(0);
            std::vector<std::thread> threads;
            for (unsigned thread_i = 0; thread_i < thread_count; ++thread_i)
            {
                threads.push_back(std::thread(
                    [codec, plan, test_case, &originals, &reference, &ready,
                     &go, &failures]() {
                    try
                    {
                        ready.fetch_add(1, std::memory_order_release);
                        while (!go.load(std::memory_order_acquire))
                            std::this_thread::yield();
                        for (unsigned iteration = 0; iteration < repetitions;
                             ++iteration)
                        {
                            const Shards recovery = execute_encode(
                                codec, test_case, originals);
                            if (recovery != reference)
                                throw std::runtime_error(
                                    "shared codec parity mismatch");
                            execute_plan(
                                plan, test_case, originals, recovery);
                        }
                    }
                    catch (...)
                    {
                        failures.fetch_add(1, std::memory_order_relaxed);
                    }
                }));
            }
            while (ready.load(std::memory_order_acquire) != thread_count)
                std::this_thread::yield();
            go.store(true, std::memory_order_release);
            for (size_t thread_i = 0; thread_i < threads.size(); ++thread_i)
                threads[thread_i].join();
            require(failures.load(std::memory_order_relaxed) == 0,
                "shared codec/plan execution failed");
            leo2_decode_plan_destroy(plan);
            leo2_codec_destroy(codec);
        }
    }
}

void test_context_batches(const std::vector<ContextEntry>& contexts)
{
    const CodecCase test_case = {
        5, 11, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8, 129
    };
    const size_t batch_count = 8;
    for (size_t context_i = 0; context_i < contexts.size(); ++context_i)
    {
        leo2_codec* codec = NULL;
        require_result(leo2_codec_create(contexts[context_i].context,
            test_case.k, test_case.r, test_case.profile, test_case.field, NULL,
            &codec), LEO2_SUCCESS, "batch codec create");
        leo2_decode_plan* plan = make_plan(codec, test_case);
        size_t encode_scratch_bytes = 0;
        size_t decode_scratch_bytes = 0;
        require_result(leo2_encode_scratch_size(codec, test_case.bytes,
            &encode_scratch_bytes), LEO2_SUCCESS,
            "batch encode scratch query");
        require_result(leo2_decode_plan_scratch_size(plan, test_case.bytes,
            &decode_scratch_bytes), LEO2_SUCCESS,
            "batch decode scratch query");

        std::vector<Shards> originals(batch_count);
        std::vector<Shards> recovery(batch_count,
            Shards(test_case.r, Bytes(test_case.bytes)));
        std::vector<Shards> restored(batch_count,
            Shards(test_case.k, Bytes(test_case.bytes)));
        std::vector<std::vector<const void*> > encode_inputs(batch_count);
        std::vector<std::vector<void*> > encode_outputs(batch_count);
        std::vector<std::vector<const void*> > decode_originals(batch_count);
        std::vector<std::vector<const void*> > decode_recovery(batch_count);
        std::vector<std::vector<void*> > decode_outputs(batch_count,
            std::vector<void*>(test_case.k, NULL));
        std::vector<std::unique_ptr<AlignedBuffer> > encode_scratch;
        std::vector<std::unique_ptr<AlignedBuffer> > decode_scratch;
        std::vector<leo2_encode_batch_item> encode_items(batch_count);
        std::vector<leo2_decode_batch_item> decode_items(batch_count);
        encode_scratch.reserve(batch_count);
        decode_scratch.reserve(batch_count);
        const uint32_t second_missing = test_case.k - 1;

        for (size_t item = 0; item < batch_count; ++item)
        {
            originals[item] = make_originals(test_case,
                static_cast<uint32_t>(0x6a09e667U + item * 313U));
            encode_inputs[item] = const_pointers(originals[item]);
            encode_outputs[item] = mutable_pointers(recovery[item]);
            encode_scratch.push_back(std::unique_ptr<AlignedBuffer>(
                new AlignedBuffer(encode_scratch_bytes)));
            leo2_encode_batch_item encode_item = {
                test_case.bytes,
                &encode_inputs[item][0],
                &encode_outputs[item][0],
                encode_scratch[item]->data(),
                encode_scratch_bytes
            };
            encode_items[item] = encode_item;

            decode_originals[item] = const_pointers(originals[item]);
            decode_originals[item][0] = NULL;
            decode_originals[item][second_missing] = NULL;
            decode_recovery[item] = const_pointers(recovery[item]);
            decode_outputs[item][0] = &restored[item][0][0];
            decode_outputs[item][second_missing] =
                &restored[item][second_missing][0];
            decode_scratch.push_back(std::unique_ptr<AlignedBuffer>(
                new AlignedBuffer(decode_scratch_bytes)));
            leo2_decode_batch_item decode_item = {
                test_case.bytes,
                &decode_originals[item][0],
                &decode_recovery[item][0],
                &decode_outputs[item][0],
                decode_scratch[item]->data(),
                decode_scratch_bytes
            };
            decode_items[item] = decode_item;
        }

        require_result(leo2_encode_batch(codec, &encode_items[0], batch_count),
            LEO2_SUCCESS, "context encode batch");
        require_result(leo2_decode_plan_execute_batch(plan, &decode_items[0],
            batch_count), LEO2_SUCCESS, "context decode batch");
        for (size_t item = 0; item < batch_count; ++item)
        {
            require(restored[item][0] == originals[item][0],
                "batch first restored original mismatch");
            require(restored[item][second_missing] ==
                    originals[item][second_missing],
                "batch second restored original mismatch");
        }
        leo2_decode_plan_destroy(plan);
        leo2_codec_destroy(codec);
    }
}

void destroy_contexts(std::vector<ContextEntry>& contexts)
{
    for (size_t i = 0; i < contexts.size(); ++i)
        leo2_context_destroy(contexts[i].context);
    contexts.clear();
}

} // namespace

int main()
{
    std::vector<ContextEntry> contexts;
    try
    {
        require(leo_init() == Leopard_Success, "Leopard initialization");
        const leopard::backend::Ops* const process_default =
            &leopard::backend::GetDefaultOps();
        test_concurrent_first_use();
        contexts = create_contexts();
        require(!contexts.empty(), "no executable contexts");
        test_process_default_immutable(process_default);
        test_traced_context_dispatch(contexts);
        test_public_codecs(contexts);
        test_shared_codec_and_plan(contexts);
        test_context_batches(contexts);
        test_concurrent_public_codecs(contexts);
        destroy_contexts(contexts);
        std::printf("Leopard2 context backends passed: "
            "profiles=2 fields=2 isolation=pass\n");
        return 0;
    }
    catch (const std::exception& error)
    {
        destroy_contexts(contexts);
        std::fprintf(stderr, "Leopard2 context backends failed: %s\n",
            error.what());
        return 1;
    }
}
