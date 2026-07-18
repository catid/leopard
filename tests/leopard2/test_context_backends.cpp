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
    std::atomic<uint64_t> ff8_ifft_four_out_calls;
    std::atomic<uint64_t> ff8_weighted_ifft_four_calls;
    std::atomic<uint64_t> ff8_fft_four_calls;
    std::atomic<uint64_t> ff8_fft_two_out_calls;
    std::atomic<uint64_t> ff8_fft_four_out_calls;
    std::atomic<uint64_t> ff8_ifft_two_xor_calls;
    std::atomic<uint64_t> ff16_ifft_four_calls;
    std::atomic<uint64_t> ff16_ifft_four_out_calls;
    std::atomic<uint64_t> ff16_weighted_ifft_four_calls;
    std::atomic<uint64_t> ff16_fft_four_calls;
    std::atomic<uint64_t> ff16_fft_two_out_calls;
    std::atomic<uint64_t> ff16_ifft_two_xor_calls;
    std::atomic<uint64_t> ff16_fft_four_out_calls;
    std::atomic<uint64_t> xor_four_calls;
    std::atomic<uint64_t> ff8_ifft_four_range_calls;
    std::atomic<uint64_t> ff8_fft_four_range_calls;
    std::atomic<uint64_t> ff8_ifft_four_xor_range_calls;
    std::atomic<uint64_t> ff16_ifft_four_range_calls;
    std::atomic<uint64_t> ff16_fft_four_range_calls;
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

void trace_ff8_fft_out(const void* x, const void* y,
    void* x_output, void* y_output, uint16_t log, uint64_t bytes)
{
    g_trace.ff8_calls.fetch_add(1, std::memory_order_relaxed);
    g_trace.ff8_fft_two_out_calls.fetch_add(1, std::memory_order_relaxed);
    trace_delegate()->ff8_fft_butterfly2_out(
        x, y, x_output, y_output, log, bytes);
}

void trace_ff8_ifft_xor(const void* x, const void* y, void* x_output,
    void* y_output, uint16_t log, uint64_t bytes)
{
    g_trace.ff8_calls.fetch_add(1, std::memory_order_relaxed);
    g_trace.ff8_ifft_two_xor_calls.fetch_add(1, std::memory_order_relaxed);
    trace_delegate()->ff8_ifft_butterfly2_xor(
        x, y, x_output, y_output, log, bytes);
}

void trace_ff16_ifft(void* x, void* y, uint16_t log, uint64_t bytes)
{
    g_trace.ff16_calls.fetch_add(1, std::memory_order_relaxed);
    trace_delegate()->ff16_ifft_butterfly2(x, y, log, bytes);
}

void trace_ff16_ifft_xor(const void* x, const void* y,
    void* x_output, void* y_output, uint16_t log, uint64_t bytes)
{
    g_trace.ff16_calls.fetch_add(1, std::memory_order_relaxed);
    g_trace.ff16_ifft_two_xor_calls.fetch_add(1,
        std::memory_order_relaxed);
    trace_delegate()->ff16_ifft_butterfly2_xor(
        x, y, x_output, y_output, log, bytes);
}

void trace_ff16_fft(void* x, void* y, uint16_t log, uint64_t bytes)
{
    g_trace.ff16_calls.fetch_add(1, std::memory_order_relaxed);
    trace_delegate()->ff16_fft_butterfly2(x, y, log, bytes);
}

void trace_ff16_fft_out(const void* x, const void* y,
    void* x_output, void* y_output, uint16_t log, uint64_t bytes)
{
    g_trace.ff16_calls.fetch_add(1, std::memory_order_relaxed);
    g_trace.ff16_fft_two_out_calls.fetch_add(1, std::memory_order_relaxed);
    trace_delegate()->ff16_fft_butterfly2_out(
        x, y, x_output, y_output, log, bytes);
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

void trace_ff8_ifft4_out(
    const void* input0, const void* input1,
    const void* input2, const void* input3,
    void* output0, void* output1, void* output2, void* output3,
    uint16_t log01, uint16_t log23, uint16_t log02, uint64_t bytes)
{
    g_trace.ff8_calls.fetch_add(1, std::memory_order_relaxed);
    g_trace.ff8_ifft_four_out_calls.fetch_add(
        1, std::memory_order_relaxed);
    trace_delegate()->ff8_ifft_butterfly4_out(
        input0, input1, input2, input3,
        output0, output1, output2, output3,
        log01, log23, log02, bytes);
}

void trace_ff8_weighted_ifft4(
    const void* input0, const void* input1,
    const void* input2, const void* input3,
    void* output0, void* output1, void* output2, void* output3,
    uint16_t weight0, uint16_t weight1,
    uint16_t weight2, uint16_t weight3,
    uint8_t live_mask,
    uint16_t log01, uint16_t log23, uint16_t log02, uint64_t bytes)
{
    g_trace.ff8_calls.fetch_add(1, std::memory_order_relaxed);
    g_trace.ff8_weighted_ifft_four_calls.fetch_add(
        1, std::memory_order_relaxed);
    trace_delegate()->ff8_weighted_ifft_butterfly4(
        input0, input1, input2, input3,
        output0, output1, output2, output3,
        weight0, weight1, weight2, weight3, live_mask,
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

void trace_ff8_fft4_out(
    const void* input0, const void* input1,
    const void* input2, const void* input3,
    void* output0, void* output1, void* output2, void* output3,
    uint16_t log01, uint16_t log23, uint16_t log02, uint64_t bytes)
{
    g_trace.ff8_calls.fetch_add(1, std::memory_order_relaxed);
    g_trace.ff8_fft_four_out_calls.fetch_add(1, std::memory_order_relaxed);
    trace_delegate()->ff8_fft_butterfly4_out(
        input0, input1, input2, input3,
        output0, output1, output2, output3,
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

void trace_ff16_ifft4_out(
    const void* input0, const void* input1,
    const void* input2, const void* input3,
    void* output0, void* output1, void* output2, void* output3,
    uint16_t log01, uint16_t log23, uint16_t log02, uint64_t bytes)
{
    g_trace.ff16_calls.fetch_add(1, std::memory_order_relaxed);
    g_trace.ff16_ifft_four_out_calls.fetch_add(
        1, std::memory_order_relaxed);
    trace_delegate()->ff16_ifft_butterfly4_out(
        input0, input1, input2, input3,
        output0, output1, output2, output3,
        log01, log23, log02, bytes);
}

void trace_ff16_weighted_ifft4(
    const void* input0, const void* input1,
    const void* input2, const void* input3,
    void* output0, void* output1, void* output2, void* output3,
    uint16_t weight0, uint16_t weight1,
    uint16_t weight2, uint16_t weight3,
    uint8_t live_mask,
    uint16_t log01, uint16_t log23, uint16_t log02, uint64_t bytes)
{
    g_trace.ff16_calls.fetch_add(1, std::memory_order_relaxed);
    g_trace.ff16_weighted_ifft_four_calls.fetch_add(
        1, std::memory_order_relaxed);
    trace_delegate()->ff16_weighted_ifft_butterfly4(
        input0, input1, input2, input3,
        output0, output1, output2, output3,
        weight0, weight1, weight2, weight3, live_mask,
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

void trace_ff16_fft4_out(
    const void* input0, const void* input1,
    const void* input2, const void* input3,
    void* output0, void* output1, void* output2, void* output3,
    uint16_t log01, uint16_t log23, uint16_t log02, uint64_t bytes)
{
    g_trace.ff16_calls.fetch_add(1, std::memory_order_relaxed);
    g_trace.ff16_fft_four_out_calls.fetch_add(1, std::memory_order_relaxed);
    trace_delegate()->ff16_fft_butterfly4_out(
        input0, input1, input2, input3,
        output0, output1, output2, output3,
        log01, log23, log02, bytes);
}

void trace_ff8_ifft4_range(
    void* const* work, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t bytes, bool prefer_fused)
{
    g_trace.ff8_calls.fetch_add(distance, std::memory_order_relaxed);
    g_trace.ff8_ifft_four_calls.fetch_add(
        distance, std::memory_order_relaxed);
    g_trace.ff8_ifft_four_range_calls.fetch_add(
        1, std::memory_order_relaxed);
    trace_delegate()->ff8_ifft_butterfly4_range(
        work, distance, log01, log23, log02, bytes, prefer_fused);
}

void trace_ff8_fft4_range(
    void* const* work, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t bytes, bool prefer_fused)
{
    g_trace.ff8_calls.fetch_add(distance, std::memory_order_relaxed);
    g_trace.ff8_fft_four_calls.fetch_add(
        distance, std::memory_order_relaxed);
    g_trace.ff8_fft_four_range_calls.fetch_add(
        1, std::memory_order_relaxed);
    trace_delegate()->ff8_fft_butterfly4_range(
        work, distance, log01, log23, log02, bytes, prefer_fused);
}

void trace_ff8_ifft4_xor_range(
    void* const* work, void* const* xor_output, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t bytes)
{
    g_trace.ff8_calls.fetch_add(distance, std::memory_order_relaxed);
    g_trace.ff8_ifft_four_calls.fetch_add(
        distance, std::memory_order_relaxed);
    g_trace.xor_four_calls.fetch_add(distance, std::memory_order_relaxed);
    g_trace.ff8_ifft_four_xor_range_calls.fetch_add(
        1, std::memory_order_relaxed);
    trace_delegate()->ff8_ifft_butterfly4_xor_range(
        work, xor_output, distance,
        log01, log23, log02, bytes);
}

void trace_ff16_ifft4_range(
    void* const* work, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t bytes, bool prefer_fused)
{
    g_trace.ff16_calls.fetch_add(distance, std::memory_order_relaxed);
    if (prefer_fused)
        g_trace.ff16_ifft_four_calls.fetch_add(
            distance, std::memory_order_relaxed);
    g_trace.ff16_ifft_four_range_calls.fetch_add(
        1, std::memory_order_relaxed);
    trace_delegate()->ff16_ifft_butterfly4_range(
        work, distance, log01, log23, log02, bytes, prefer_fused);
}

void trace_ff16_fft4_range(
    void* const* work, unsigned distance,
    uint16_t log01, uint16_t log23, uint16_t log02,
    uint64_t bytes, bool prefer_fused)
{
    g_trace.ff16_calls.fetch_add(distance, std::memory_order_relaxed);
    if (prefer_fused)
        g_trace.ff16_fft_four_calls.fetch_add(
            distance, std::memory_order_relaxed);
    g_trace.ff16_fft_four_range_calls.fetch_add(
        1, std::memory_order_relaxed);
    trace_delegate()->ff16_fft_butterfly4_range(
        work, distance, log01, log23, log02, bytes, prefer_fused);
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
        tracing_.ff8_fft_butterfly2_out = trace_ff8_fft_out;
        tracing_.ff8_ifft_butterfly2_xor = trace_ff8_ifft_xor;
        tracing_.ff8_ifft_butterfly4 = trace_ff8_ifft4;
        tracing_.ff8_ifft_butterfly4_out = trace_ff8_ifft4_out;
        tracing_.ff8_weighted_ifft_butterfly4 = trace_ff8_weighted_ifft4;
        tracing_.ff8_fft_butterfly4 = trace_ff8_fft4;
        tracing_.ff8_fft_butterfly4_out = trace_ff8_fft4_out;
        tracing_.ff8_ifft_butterfly4_range = trace_ff8_ifft4_range;
        tracing_.ff8_fft_butterfly4_range = trace_ff8_fft4_range;
        tracing_.ff8_ifft_butterfly4_xor_range =
            trace_ff8_ifft4_xor_range;
        tracing_.ff16_ifft_butterfly2 = trace_ff16_ifft;
        tracing_.ff16_fft_butterfly2 = trace_ff16_fft;
        tracing_.ff16_fft_butterfly2_out = trace_ff16_fft_out;
        tracing_.ff16_ifft_butterfly2_xor = trace_ff16_ifft_xor;
        tracing_.ff16_ifft_butterfly4 = trace_ff16_ifft4;
        tracing_.ff16_ifft_butterfly4_out = trace_ff16_ifft4_out;
        tracing_.ff16_weighted_ifft_butterfly4 = trace_ff16_weighted_ifft4;
        tracing_.ff16_fft_butterfly4 = trace_ff16_fft4;
        tracing_.ff16_fft_butterfly4_out = trace_ff16_fft4_out;
        tracing_.ff16_ifft_butterfly4_range = trace_ff16_ifft4_range;
        tracing_.ff16_fft_butterfly4_range = trace_ff16_fft4_range;
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
        g_trace.ff8_ifft_four_out_calls.store(0, std::memory_order_relaxed);
        g_trace.ff8_weighted_ifft_four_calls.store(
            0, std::memory_order_relaxed);
        g_trace.ff8_fft_four_calls.store(0, std::memory_order_relaxed);
        g_trace.ff8_fft_two_out_calls.store(0, std::memory_order_relaxed);
        g_trace.ff8_fft_four_out_calls.store(0, std::memory_order_relaxed);
        g_trace.ff8_ifft_two_xor_calls.store(0, std::memory_order_relaxed);
        g_trace.ff16_ifft_four_calls.store(0, std::memory_order_relaxed);
        g_trace.ff16_ifft_four_out_calls.store(0, std::memory_order_relaxed);
        g_trace.ff16_weighted_ifft_four_calls.store(
            0, std::memory_order_relaxed);
        g_trace.ff16_fft_four_calls.store(0, std::memory_order_relaxed);
        g_trace.ff16_fft_two_out_calls.store(0, std::memory_order_relaxed);
        g_trace.ff16_ifft_two_xor_calls.store(
            0, std::memory_order_relaxed);
        g_trace.ff16_fft_four_out_calls.store(0, std::memory_order_relaxed);
        g_trace.xor_four_calls.store(0, std::memory_order_relaxed);
        g_trace.ff8_ifft_four_range_calls.store(
            0, std::memory_order_relaxed);
        g_trace.ff8_fft_four_range_calls.store(
            0, std::memory_order_relaxed);
        g_trace.ff8_ifft_four_xor_range_calls.store(
            0, std::memory_order_relaxed);
        g_trace.ff16_ifft_four_range_calls.store(
            0, std::memory_order_relaxed);
        g_trace.ff16_fft_four_range_calls.store(
            0, std::memory_order_relaxed);
        leopard::ff8::TestOnlyResetTransformCallsiteCounts();
        leopard::ff16::TestOnlyResetTransformCallsiteCounts();
        leopard::ff8::TestOnlyResetLowEncodeCounts();
        leopard::ff16::TestOnlyResetLowEncodeCounts();
        leopard::ff8::TestOnlyResetHighEncodeCounts();
        leopard::ff16::TestOnlyResetHighEncodeCounts();
        leopard::ff8::TestOnlyResetSparseEncodeCounts();
        leopard::ff16::TestOnlyResetSparseEncodeCounts();
        leopard::ff8::TestOnlyResetHighDecodeCounts();
        leopard::ff16::TestOnlyResetHighDecodeCounts();
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
    uint64_t ff8_ifft_four_out_calls() const
    {
        return g_trace.ff8_ifft_four_out_calls.load(
            std::memory_order_relaxed);
    }
    uint64_t ff8_weighted_ifft_four_calls() const
    {
        return g_trace.ff8_weighted_ifft_four_calls.load(
            std::memory_order_relaxed);
    }
    uint64_t ff8_fft_four_calls() const
    {
        return g_trace.ff8_fft_four_calls.load(std::memory_order_relaxed);
    }
    uint64_t ff8_fft_two_out_calls() const
    {
        return g_trace.ff8_fft_two_out_calls.load(std::memory_order_relaxed);
    }
    uint64_t ff8_fft_four_out_calls() const
    {
        return g_trace.ff8_fft_four_out_calls.load(std::memory_order_relaxed);
    }
    uint64_t ff8_ifft_two_xor_calls() const
    {
        return g_trace.ff8_ifft_two_xor_calls.load(
            std::memory_order_relaxed);
    }
    uint64_t ff16_ifft_four_calls() const
    {
        return g_trace.ff16_ifft_four_calls.load(std::memory_order_relaxed);
    }
    uint64_t ff16_ifft_four_out_calls() const
    {
        return g_trace.ff16_ifft_four_out_calls.load(
            std::memory_order_relaxed);
    }
    uint64_t ff16_weighted_ifft_four_calls() const
    {
        return g_trace.ff16_weighted_ifft_four_calls.load(
            std::memory_order_relaxed);
    }
    uint64_t ff16_fft_four_calls() const
    {
        return g_trace.ff16_fft_four_calls.load(std::memory_order_relaxed);
    }
    uint64_t ff16_fft_two_out_calls() const
    {
        return g_trace.ff16_fft_two_out_calls.load(std::memory_order_relaxed);
    }
    uint64_t ff16_ifft_two_xor_calls() const
    {
        return g_trace.ff16_ifft_two_xor_calls.load(
            std::memory_order_relaxed);
    }
    uint64_t ff16_fft_four_out_calls() const
    {
        return g_trace.ff16_fft_four_out_calls.load(std::memory_order_relaxed);
    }
    uint64_t xor_four_calls() const
    {
        return g_trace.xor_four_calls.load(std::memory_order_relaxed);
    }
    uint64_t ff8_ifft_four_range_calls() const
    {
        return g_trace.ff8_ifft_four_range_calls.load(
            std::memory_order_relaxed);
    }
    uint64_t ff8_fft_four_range_calls() const
    {
        return g_trace.ff8_fft_four_range_calls.load(
            std::memory_order_relaxed);
    }
    uint64_t ff8_ifft_four_xor_range_calls() const
    {
        return g_trace.ff8_ifft_four_xor_range_calls.load(
            std::memory_order_relaxed);
    }
    uint64_t ff16_ifft_four_range_calls() const
    {
        return g_trace.ff16_ifft_four_range_calls.load(
            std::memory_order_relaxed);
    }
    uint64_t ff16_fft_four_range_calls() const
    {
        return g_trace.ff16_fft_four_range_calls.load(
            std::memory_order_relaxed);
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
    bool split_ragged_execution,
    bool allow_pruned_schedule_calls = false)
{
    if (field == LEO2_FIELD_GF8)
    {
        const leopard::ff8::TestOnlyTransformCallsiteCounts callsites =
            leopard::ff8::TestOnlyGetTransformCallsiteCounts();
        require(allow_pruned_schedule_calls || callsites.ifft_dit4 != 0,
            operation + " did not exercise GF8 IFFT_DIT4");
        if (expect_ff8_accumulating_ifft)
            require(callsites.ifft_dit4_xor != 0,
                operation + " did not exercise GF8 IFFT_DIT4_xor");
        else
            require(callsites.ifft_dit4_xor == 0,
                operation + " unexpectedly exercised GF8 IFFT_DIT4_xor");
        require(allow_pruned_schedule_calls || callsites.fft_dit4 != 0,
            operation + " did not exercise GF8 FFT_DIT4");
        const uint64_t classified_ifft_calls =
            callsites.ifft_dit4 + callsites.ifft_dit4_xor;
        if (!allow_pruned_schedule_calls &&
            trace.ff8_ifft_four_range_calls() +
                trace.ff8_ifft_four_xor_range_calls() != 0)
        {
            const uint64_t inverse_ranges =
                trace.ff8_ifft_four_range_calls() +
                trace.ff8_ifft_four_xor_range_calls();
            require(inverse_ranges <= classified_ifft_calls,
                operation + " did not use the GF8 inverse stage boundary");
            require(trace.ff8_fft_four_range_calls() <=
                        callsites.fft_dit4,
                operation + " did not use the GF8 forward stage boundary");
        }
        require(allow_pruned_schedule_calls
                ? trace.ff8_ifft_four_calls() >= classified_ifft_calls
                : trace.ff8_ifft_four_calls() == classified_ifft_calls,
            operation + " has a GF8 inverse radix-four callsite bypass");
        require(allow_pruned_schedule_calls
                ? trace.ff8_fft_four_calls() >= callsites.fft_dit4
                : trace.ff8_fft_four_calls() == callsites.fft_dit4,
            operation + " has a GF8 forward radix-four callsite bypass");
    }
    else
    {
        const leopard::ff16::TestOnlyTransformCallsiteCounts callsites =
            leopard::ff16::TestOnlyGetTransformCallsiteCounts();
        require(allow_pruned_schedule_calls || callsites.ifft_dit4 != 0,
            operation + " did not exercise GF16 IFFT_DIT4");
        require(allow_pruned_schedule_calls || callsites.fft_dit4 != 0,
            operation + " did not exercise GF16 FFT_DIT4");
        require(callsites.ifft_dit4 ==
                callsites.ifft_dit4_fused + callsites.ifft_dit4_split,
            operation + " has an unclassified GF16 inverse radix-four callsite");
        require(callsites.fft_dit4 ==
                callsites.fft_dit4_fused + callsites.fft_dit4_split,
            operation + " has an unclassified GF16 forward radix-four callsite");
        if (!allow_pruned_schedule_calls &&
            trace.ff16_ifft_four_range_calls() != 0)
        {
            require(trace.ff16_ifft_four_range_calls() <=
                        callsites.ifft_dit4,
                operation + " did not use the GF16 inverse stage boundary");
            require(trace.ff16_fft_four_range_calls() <=
                        callsites.fft_dit4,
                operation + " did not use the GF16 forward stage boundary");
        }
        require(allow_pruned_schedule_calls
                ? trace.ff16_ifft_four_calls() >= callsites.ifft_dit4_fused
                : trace.ff16_ifft_four_calls() == callsites.ifft_dit4_fused,
            operation + " has a GF16 inverse fused-dispatch mismatch");
        require(allow_pruned_schedule_calls
                ? trace.ff16_fft_four_calls() >= callsites.fft_dit4_fused
                : trace.ff16_fft_four_calls() == callsites.fft_dit4_fused,
            operation + " has a GF16 forward fused-dispatch mismatch");

        // A C1 exact-mask schedule can replace every mature radix-four
        // callsite in one direction.  The tracing table above still proves
        // that the context-selected backend handled the pruned butterflies;
        // the fused/split counters below classify only mature calls.
        if (allow_pruned_schedule_calls &&
            (callsites.ifft_dit4 == 0 || callsites.fft_dit4 == 0))
            return;

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

void require_low_encode_no_copy(
    const TraceOpsGuard& trace,
    leo2_backend execution_backend,
    leo2_field field,
    size_t public_bytes,
    const std::string& operation)
{
    if (field == LEO2_FIELD_GF8)
    {
        const leopard::ff8::TestOnlyLowEncodeCounts counts =
            leopard::ff8::TestOnlyGetLowEncodeCounts();
        require(counts.fft_butterfly2_out_of_place +
                counts.fft_butterfly4_out_of_place != 0,
            operation + " did not enter an out-of-place first layer");
        require(trace.ff8_fft_two_out_calls() ==
                    counts.fft_butterfly2_out_of_place &&
                trace.ff8_fft_four_out_calls() ==
                    counts.fft_butterfly4_out_of_place,
            operation + " bypassed the context out-of-place ops table");
    }
    else
    {
        const leopard::ff16::TestOnlyLowEncodeCounts counts =
            leopard::ff16::TestOnlyGetLowEncodeCounts();
        require(counts.fft_butterfly2_out_of_place +
                counts.fft_butterfly4_out_of_place != 0,
            operation + " did not enter an out-of-place first layer");
        require(trace.ff16_fft_two_out_calls() ==
                    counts.fft_butterfly2_out_of_place &&
                trace.ff16_fft_four_out_calls() ==
                    counts.fft_butterfly4_out_of_place,
            operation + " bypassed the context out-of-place ops table");

        const size_t prefix_bytes = public_bytes & ~size_t(63U);
        const bool has_tail = (public_bytes & 63U) != 0;
        const bool prefix_fused = prefix_bytes == 64U ||
            (prefix_bytes == 128U &&
             execution_backend == LEO2_BACKEND_AVX2);
        // A ragged tail is a separate padded 64-byte transform and therefore
        // always uses the qualified fused first layer.  The aligned prefix
        // retains its own size/backend policy.
        const bool expect_all_fused = has_tail
            ? prefix_bytes == 0 || prefix_fused
            : prefix_fused;
        const bool expect_mixed = has_tail && prefix_bytes != 0 &&
            !prefix_fused;
        if (expect_mixed)
        {
            require(counts.fft_butterfly2_out_of_place != 0 &&
                    counts.fft_butterfly4_out_of_place != 0,
                operation + " did not preserve the split-prefix/tail policy");
        }
        else if (expect_all_fused)
        {
            require(counts.fft_butterfly2_out_of_place == 0 &&
                    counts.fft_butterfly4_out_of_place != 0,
                operation + " split a qualified GF16 first layer");
        }
        else
        {
            require(counts.fft_butterfly2_out_of_place != 0 &&
                    counts.fft_butterfly4_out_of_place == 0,
                operation + " fused an unqualified GF16 first layer");
        }
    }
}

void require_high_encode_source_staging(
    const TraceOpsGuard& trace,
    leo2_field field,
    const std::string& operation)
{
    uint64_t staged_groups = 0;
    uint64_t traced_groups = 0;
    if (field == LEO2_FIELD_GF8)
    {
        const leopard::ff8::TestOnlyHighEncodeCounts counts =
            leopard::ff8::TestOnlyGetHighEncodeCounts();
        staged_groups = counts.ifft_butterfly4_out_of_place;
        traced_groups = trace.ff8_ifft_four_out_calls();
    }
    else
    {
        const leopard::ff16::TestOnlyHighEncodeCounts counts =
            leopard::ff16::TestOnlyGetHighEncodeCounts();
        staged_groups = counts.ifft_butterfly4_out_of_place;
        traced_groups = trace.ff16_ifft_four_out_calls();
    }
    require(staged_groups != 0,
        operation + " did not stage caller sources through inverse radix four");
    require(traced_groups == staged_groups,
        operation + " bypassed the context inverse out-of-place ops table");
}

void require_coarse_stage_reduction(
    const TraceOpsGuard& trace,
    leo2_field field,
    const std::string& operation)
{
    if (field == LEO2_FIELD_GF8)
    {
        const leopard::ff8::TestOnlyTransformCallsiteCounts callsites =
            leopard::ff8::TestOnlyGetTransformCallsiteCounts();
        const uint64_t inverse_leaf_calls =
            callsites.ifft_dit4 + callsites.ifft_dit4_xor;
        const uint64_t inverse_range_calls =
            trace.ff8_ifft_four_range_calls() +
            trace.ff8_ifft_four_xor_range_calls();
        require(inverse_range_calls != 0 &&
                inverse_range_calls < inverse_leaf_calls,
            operation + " did not reduce GF8 inverse indirect dispatches");
        require(trace.ff8_fft_four_range_calls() != 0 &&
                trace.ff8_fft_four_range_calls() < callsites.fft_dit4,
            operation + " did not reduce GF8 forward indirect dispatches");
    }
    else
    {
        const leopard::ff16::TestOnlyTransformCallsiteCounts callsites =
            leopard::ff16::TestOnlyGetTransformCallsiteCounts();
        require(trace.ff16_ifft_four_range_calls() != 0 &&
                trace.ff16_ifft_four_range_calls() < callsites.ifft_dit4,
            operation + " did not reduce GF16 inverse indirect dispatches");
        require(trace.ff16_fft_four_range_calls() != 0 &&
                trace.ff16_fft_four_range_calls() < callsites.fft_dit4,
            operation + " did not reduce GF16 forward indirect dispatches");
    }
}

void require_high_decode_no_copy(
    const TraceOpsGuard& trace,
    leo2_backend execution_backend,
    leo2_field field,
    bool expect_pruned_accumulation,
    bool expect_receive_fusion,
    bool expect_weighted_boundary,
    bool expect_gf16_accumulation,
    bool expect_gf16_materialization,
    const std::string& operation)
{
    uint64_t output_blocks = 0;
    uint64_t butterfly2 = 0;
    uint64_t butterfly4 = 0;
    uint64_t copy_fallbacks = 0;
    uint64_t traced_butterfly2 = 0;
    uint64_t traced_butterfly4 = 0;
    uint64_t syndrome_accumulated_blocks = 0;
    uint64_t syndrome_materialized_blocks = 0;
    uint64_t syndrome_pruned_accumulated_blocks = 0;
    uint64_t syndrome_pruned_fallback_blocks = 0;
    uint64_t receive_butterfly4 = 0;
    uint64_t receive_copy_shards = 0;
    uint64_t receive_zero_shards = 0;
    uint64_t traced_receive_butterfly4 = 0;
    uint64_t locator_weighted_butterfly4 = 0;
    uint64_t locator_scale_rows_elided = 0;
    uint64_t locator_inactive_rows = 0;
    uint64_t traced_locator_weighted_butterfly4 = 0;
    if (field == LEO2_FIELD_GF8)
    {
        const leopard::ff8::TestOnlyHighDecodeCounts counts =
            leopard::ff8::TestOnlyGetHighDecodeCounts();
        output_blocks = counts.output_blocks;
        butterfly2 = counts.fft_butterfly2_out_of_place;
        butterfly4 = counts.fft_butterfly4_out_of_place;
        copy_fallbacks = counts.compatibility_copy_fallbacks;
        traced_butterfly2 = trace.ff8_fft_two_out_calls();
        traced_butterfly4 = trace.ff8_fft_four_out_calls();
        syndrome_accumulated_blocks = counts.syndrome_accumulated_blocks;
        syndrome_materialized_blocks = counts.syndrome_materialized_blocks;
        syndrome_pruned_accumulated_blocks =
            counts.syndrome_pruned_accumulated_blocks;
        syndrome_pruned_fallback_blocks =
            counts.syndrome_pruned_fallback_blocks;
        receive_butterfly4 =
            counts.receive_ifft_butterfly4_out_of_place;
        receive_copy_shards = counts.receive_copy_shards;
        receive_zero_shards = counts.receive_zero_shards;
        traced_receive_butterfly4 = trace.ff8_ifft_four_out_calls();
        locator_weighted_butterfly4 =
            counts.locator_weighted_ifft_butterfly4;
        locator_scale_rows_elided = counts.locator_scale_rows_elided;
        locator_inactive_rows = counts.locator_inactive_rows;
        traced_locator_weighted_butterfly4 =
            trace.ff8_weighted_ifft_four_calls();
    }
    else
    {
        const leopard::ff16::TestOnlyHighDecodeCounts counts =
            leopard::ff16::TestOnlyGetHighDecodeCounts();
        output_blocks = counts.output_blocks;
        butterfly2 = counts.fft_butterfly2_out_of_place;
        butterfly4 = counts.fft_butterfly4_out_of_place;
        copy_fallbacks = counts.compatibility_copy_fallbacks;
        traced_butterfly2 = trace.ff16_fft_two_out_calls();
        traced_butterfly4 = trace.ff16_fft_four_out_calls();
        syndrome_accumulated_blocks = counts.syndrome_accumulated_blocks;
        syndrome_materialized_blocks = counts.syndrome_materialized_blocks;
        syndrome_pruned_accumulated_blocks =
            counts.syndrome_pruned_accumulated_blocks;
        syndrome_pruned_fallback_blocks =
            counts.syndrome_pruned_fallback_blocks;
        receive_butterfly4 =
            counts.receive_ifft_butterfly4_out_of_place;
        receive_copy_shards = counts.receive_copy_shards;
        receive_zero_shards = counts.receive_zero_shards;
        traced_receive_butterfly4 = trace.ff16_ifft_four_out_calls();
        locator_weighted_butterfly4 =
            counts.locator_weighted_ifft_butterfly4;
        locator_scale_rows_elided = counts.locator_scale_rows_elided;
        locator_inactive_rows = counts.locator_inactive_rows;
        traced_locator_weighted_butterfly4 =
            trace.ff16_weighted_ifft_four_calls();
    }
    const uint64_t out_of_place_calls = butterfly2 + butterfly4;
    require(output_blocks != 0,
        operation + " did not evaluate a requested output block");
    require(out_of_place_calls != 0,
        operation + " did not enter an immutable-source first layer");
    require(copy_fallbacks == 0,
        operation + " unexpectedly entered the whole-block copy fallback");
    require(traced_butterfly2 == butterfly2 &&
            traced_butterfly4 == butterfly4,
        operation + " bypassed the selected context out-of-place table");
    require(traced_receive_butterfly4 == receive_butterfly4,
        operation + " bypassed the selected context inverse source table");
    require(traced_locator_weighted_butterfly4 ==
            locator_weighted_butterfly4,
        operation + " bypassed the selected context weighted boundary table");
    require(locator_scale_rows_elided ==
            locator_weighted_butterfly4 * 4U,
        operation + " did not account for every elided locator-scale row");
    require(locator_inactive_rows <= locator_scale_rows_elided,
        operation + " over-counted masked locator rows");
    if (expect_weighted_boundary)
        require(locator_weighted_butterfly4 != 0,
            operation + " did not remove the Algorithm 5 locator pass");
    else
        require(locator_weighted_butterfly4 == 0,
            operation + " widened the weighted boundary below four rows");
    if (expect_receive_fusion)
        require(receive_butterfly4 != 0,
            operation + " did not fuse a complete unpruned receive group");
    else
        require(receive_butterfly4 == 0,
            operation + " widened receive fusion beyond qualified backends");
    require(receive_butterfly4 != 0 ||
            receive_copy_shards + receive_zero_shards != 0,
        operation + " did not account for receive staging");
    require(syndrome_pruned_fallback_blocks == 0,
        operation + " rejected a compiled exact inverse sink");
    if (expect_pruned_accumulation)
        require(syndrome_pruned_accumulated_blocks != 0,
            operation + " did not select an eligible exact inverse sink");
    else
        require(syndrome_pruned_accumulated_blocks == 0,
            operation + " bypassed the measured materialized policy");
    if (field == LEO2_FIELD_GF8)
    {
        require(syndrome_accumulated_blocks != 0,
            operation + " did not fuse any later syndrome input block");
        require(trace.ff8_ifft_two_xor_calls() +
                    trace.ff8_ifft_four_xor_range_calls() != 0,
            operation + " bypassed the accumulating IFFT backend boundary");
    }
    else if (execution_backend == LEO2_BACKEND_SCALAR)
    {
        require(syndrome_accumulated_blocks == 0 &&
                syndrome_materialized_blocks != 0,
            operation + " did not retain the qualified scalar fallback");
        require(trace.ff16_ifft_two_xor_calls() == 0,
            operation + " unexpectedly used vector accumulation on scalar");
    }
    else
    {
        require(syndrome_accumulated_blocks != 0 ||
                syndrome_materialized_blocks != 0,
            operation + " did not classify any later GF16 syndrome block");
        if (expect_gf16_accumulation)
            require(syndrome_accumulated_blocks != 0,
                operation + " did not fuse any GF16 syndrome input block");
        if (expect_gf16_materialization)
            require(syndrome_materialized_blocks != 0,
                operation + " did not retain its GF16 materialized fallback");
        require((trace.ff16_ifft_two_xor_calls() != 0) ==
                    (syndrome_accumulated_blocks != 0),
            operation + " misclassified the GF16 accumulating IFFT route");
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

uint32_t transform_side(const CodecCase& test_case)
{
    const uint32_t active = test_case.profile == LEO2_PROFILE_LEGACY_HIGH_V1
        ? test_case.r : test_case.k;
    uint32_t result = 1;
    while (result < active)
        result <<= 1;
    return result;
}

bool final_inverse_layer_is_radix_four(uint32_t side)
{
    unsigned log2_side = 0;
    while ((uint32_t(1) << log2_side) < side)
        ++log2_side;
    return side >= 4 && (log2_side & 1U) == 0;
}

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
    const Shards& originals,
    bool sparse = false)
{
    size_t scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(codec, test_case.bytes,
        &scratch_bytes), LEO2_SUCCESS, "encode scratch query");
    AlignedBuffer scratch(scratch_bytes);
    Shards recovery(test_case.r, Bytes(test_case.bytes));
    const std::vector<const void*> original_pointers =
        const_pointers(originals);
    std::vector<void*> recovery_pointers = mutable_pointers(recovery);
    if (sparse)
        for (uint32_t i = 0; i < test_case.r; ++i)
            if (i != 0 && i + 1 != test_case.r && i % 7U != 3U)
                recovery_pointers[i] = NULL;
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
        // Multiple high-rate message blocks force the final GF8 IFFT layer to
        // accumulate into the first block.  Cover an exact vector tile and a
        // ragged public tail so scalar, SSSE3, and AVX2 must agree through the
        // fused IFFT-plus-XOR backend boundary.
        { 33, 16, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 64 },
        { 33, 16, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 129 },
        { 33, 16, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 257 },
        // Crossover boundaries: AVX2 selects the exact sink at 1 KiB, while
        // SSSE3 retains materialization until 64 KiB.
        { 33, 16, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 1024 },
        { 33, 16, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 65536 },
        // Algorithm 5 final-layer shapes: T=2 exercises the accumulating
        // radix-two boundary, T=4 the distance-one accumulating radix-four
        // boundary, and T=8 odd-log2 final radix two.
        { 17, 2, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 129 },
        { 17, 4, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 193 },
        { 17, 8, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 257 },
        { 17, 33, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8, 129 },
        // GF16 Algorithm 5 syndrome accumulation covers both possible final
        // layer shapes.  T=4 uses a split 192-byte prefix plus a compact
        // fused tail; T=16 similarly mixes a 1024-byte split prefix with a
        // padded tail.  Scalar remains the measured materialized fallback.
        { 17, 2, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, 130 },
        { 17, 4, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, 194 },
        { 17, 8, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, 258 },
        { 33, 16, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, 1026 },
        // The 64-KiB production target remains materialized on GF16 SSSE3
        // because its measured 4.86% win did not clear the 5% promotion rule;
        // AVX2 still selects the sink for the same case.
        { 33, 16, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, 65536 },
        { 33, 17, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, 64 },
        { 17, 33, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16, 64 },
        { 33, 17, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, 66 },
        { 17, 33, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16, 66 },
        { 33, 17, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, 128 },
        { 17, 33, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16, 128 },
        { 33, 17, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, 130 },
        { 17, 33, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16, 130 },
        { 33, 17, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, 1026 },
        // A 64-wide redundancy side ends on radix four rather than the odd
        // final layer above.  The unqualified 1024-byte prefix must retain
        // the split policy while fusing its second layer with accumulation;
        // the final two public bytes separately cover the padded tail.
        { 65, 63, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, 1026 },
        { 17, 33, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16, 1026 }
    };
    std::vector<Shards> reference_recovery(
        sizeof(transform_cases) / sizeof(transform_cases[0]));

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
            const uint32_t side = transform_side(test_case);
            const Shards originals = make_originals(test_case,
                static_cast<uint32_t>(0x243f6a88U + case_i * 977U));
            leo2_codec* codec = NULL;
            trace.reset();
            const Shards recovery = encode_case(contexts[context_i].context,
                test_case, originals, &codec);
            if (reference_recovery[case_i].empty())
                reference_recovery[case_i] = recovery;
            else
                require(recovery == reference_recovery[case_i],
                    "encode output changed across context backends");
            if (test_case.field == LEO2_FIELD_GF8)
            {
                require(trace.ff8_calls() != 0,
                    "GF8 encode bypassed the context ops table");
            }
            else
            {
                require(trace.ff16_calls() != 0,
                    "GF16 encode bypassed the context ops table");
                if (test_case.profile == LEO2_PROFILE_LEGACY_HIGH_V1 &&
                    leo2_context_backend(contexts[context_i].context) !=
                        LEO2_BACKEND_SCALAR)
                {
                    require(trace.ff16_ifft_two_xor_calls() != 0,
                        "GF16 high encode did not fuse its final inverse "
                        "layer with accumulation: backend=" +
                        std::to_string(static_cast<unsigned>(
                            leo2_context_backend(
                                contexts[context_i].context))) +
                        " K=" + std::to_string(test_case.k) +
                        " R=" + std::to_string(test_case.r) +
                        " bytes=" + std::to_string(test_case.bytes));
                }
            }
            if (side >= 4)
            {
                require_four_way_callsites(trace,
                    leo2_context_backend(contexts[context_i].context),
                    test_case.field,
                    test_case.bytes,
                    profile_name + " encode",
                    test_case.field == LEO2_FIELD_GF8 &&
                        test_case.profile == LEO2_PROFILE_LEGACY_HIGH_V1 &&
                        final_inverse_layer_is_radix_four(side),
                    true);
                if (side >= 16)
                    require_coarse_stage_reduction(
                        trace, test_case.field, profile_name + " encode");
            }
            else
            {
                require(test_case.profile == LEO2_PROFILE_LEGACY_HIGH_V1,
                    profile_name + " T=2 encode has an invalid profile");
                if (test_case.field == LEO2_FIELD_GF8)
                {
                    require(trace.ff8_ifft_two_xor_calls() != 0,
                        profile_name +
                            " T=2 encode missed accumulating GF8 IFFT2");
                }
                else if (leo2_context_backend(
                             contexts[context_i].context) ==
                         LEO2_BACKEND_SCALAR)
                {
                    require(trace.ff16_ifft_two_xor_calls() == 0,
                        profile_name +
                            " T=2 scalar encode used vector accumulation");
                }
                else
                {
                    require(trace.ff16_ifft_two_xor_calls() != 0,
                        profile_name +
                            " T=2 encode missed accumulating GF16 IFFT2");
                }
            }
            if (test_case.profile == LEO2_PROFILE_LOW_V1)
                require_low_encode_no_copy(
                    trace,
                    leo2_context_backend(contexts[context_i].context),
                    test_case.field, test_case.bytes,
                    profile_name + " encode");
            else if (side > 4)
                require_high_encode_source_staging(
                    trace, test_case.field, profile_name + " encode");

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
                profile_name + " decode",
                test_case.field == LEO2_FIELD_GF8 &&
                    test_case.profile == LEO2_PROFILE_LEGACY_HIGH_V1 &&
                    final_inverse_layer_is_radix_four(side),
                true,
                true);
            uint64_t exact_pruned_syndrome_blocks = 0;
            if (test_case.profile == LEO2_PROFILE_LEGACY_HIGH_V1)
            {
                if (test_case.field == LEO2_FIELD_GF8)
                    exact_pruned_syndrome_blocks = leopard::ff8::
                        TestOnlyGetHighDecodeCounts().
                            syndrome_pruned_accumulated_blocks;
                else
                    exact_pruned_syndrome_blocks = leopard::ff16::
                        TestOnlyGetHighDecodeCounts().
                            syndrome_pruned_accumulated_blocks;
            }
            if (side >= 4 && exact_pruned_syndrome_blocks == 0)
                require(trace.xor_four_calls() != 0,
                    "decode bypassed the context grouped-XOR table");
            if (test_case.profile == LEO2_PROFILE_LEGACY_HIGH_V1)
            {
                const bool expect_gf16_accumulation =
                    test_case.field == LEO2_FIELD_GF16 &&
                    ((test_case.k == 17 &&
                      (test_case.r == 2 || test_case.r == 4 ||
                       test_case.r == 8)) ||
                     (test_case.k == 33 && test_case.r == 16));
                const leo2_backend execution_backend =
                    leo2_context_backend(contexts[context_i].context);
                const bool expect_pruned_accumulation =
                    (execution_backend == LEO2_BACKEND_AVX2 &&
                     test_case.bytes >= 1024) ||
                    (execution_backend == LEO2_BACKEND_SSSE3 &&
                     test_case.field == LEO2_FIELD_GF8 &&
                     test_case.bytes >= 65536);
                const bool expect_gf16_materialization =
                    test_case.field == LEO2_FIELD_GF16 &&
                    ((test_case.k == 17 && test_case.r == 4) ||
                     (test_case.k == 33 && test_case.r == 16)) &&
                    (!expect_pruned_accumulation ||
                     (test_case.bytes & 63U) != 0);
                require_high_decode_no_copy(
                    trace,
                    execution_backend, test_case.field,
                    expect_pruned_accumulation,
                    test_case.field == LEO2_FIELD_GF8 &&
                        side > 4 && test_case.r == side &&
                        (execution_backend == LEO2_BACKEND_SSSE3 ||
                         execution_backend == LEO2_BACKEND_AVX2),
                    side >= 4,
                    expect_gf16_accumulation,
                    expect_gf16_materialization,
                    profile_name + " decode backend=" +
                        std::to_string(static_cast<unsigned>(
                            leo2_context_backend(
                                contexts[context_i].context))) +
                        " K=" + std::to_string(test_case.k) +
                        " R=" + std::to_string(test_case.r) +
                        " bytes=" + std::to_string(test_case.bytes));
            }
            leo2_codec_destroy(codec);
        }

        const CodecCase generic_cases[] = {
            { 33, 17, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 257 },
            { 33, 17, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, 258 }
        };
        for (size_t generic_i = 0;
             generic_i < sizeof(generic_cases) / sizeof(generic_cases[0]);
             ++generic_i)
        {
            const CodecCase& generic_case = generic_cases[generic_i];
            const Shards generic_originals = make_originals(
                generic_case, static_cast<uint32_t>(0xa4093822U + generic_i));
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
            if (generic_case.field == LEO2_FIELD_GF8)
                require(trace.ff8_calls() != 0,
                    "generic GF8 decode bypassed the context ops table");
            else
                require(trace.ff16_calls() != 0,
                    "generic GF16 decode bypassed the context ops table");
            require_four_way_callsites(trace,
                leo2_context_backend(contexts[context_i].context),
                generic_case.field,
                generic_case.bytes,
                "generic high-profile decode", false, true);
            require(trace.xor_four_calls() != 0,
                "generic decode bypassed the context grouped-XOR table");
            leo2_codec_destroy(generic_codec);
        }

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

void test_sparse_encode_contexts(
    const std::vector<ContextEntry>& contexts)
{
    const CodecCase cases[] = {
        { 65, 63, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, 257 },
        { 17, 65, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF8, 129 },
        { 257, 63, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, 66 },
        { 65, 129, LEO2_PROFILE_LOW_V1, LEO2_FIELD_GF16, 130 }
    };
    for (size_t case_i = 0; case_i < sizeof(cases) / sizeof(cases[0]); ++case_i)
    {
        const CodecCase& test_case = cases[case_i];
        const Shards originals = make_originals(test_case,
            static_cast<uint32_t>(0x13198a2eU + case_i * 193U));
        Shards reference;
        for (size_t context_i = 0; context_i < contexts.size(); ++context_i)
        {
            leo2_codec* codec = NULL;
            require_result(leo2_codec_create(contexts[context_i].context,
                test_case.k, test_case.r, test_case.profile, test_case.field,
                NULL, &codec), LEO2_SUCCESS, "sparse codec create");
            require_result(leo2_test_codec_set_encode_mode(
                codec, LEO2_TEST_ENCODE_FORCE_TRANSFORM),
                LEO2_SUCCESS, "force sparse transform");
            if (test_case.field == LEO2_FIELD_GF8)
                leopard::ff8::TestOnlyResetSparseEncodeCounts();
            else
                leopard::ff16::TestOnlyResetSparseEncodeCounts();
            const Shards actual = execute_encode(
                codec, test_case, originals, true);
            if (reference.empty())
                reference = actual;
            else
                require(actual == reference,
                    "sparse parity differs between context backends");

            uint64_t exact_blocks = 0;
            uint64_t prefix_butterflies = 0;
            uint64_t retained_butterflies = 0;
            if (test_case.field == LEO2_FIELD_GF8)
            {
                const leopard::ff8::TestOnlySparseEncodeCounts counts =
                    leopard::ff8::TestOnlyGetSparseEncodeCounts();
                exact_blocks = counts.exact_blocks;
                prefix_butterflies = counts.prefix_butterflies;
                retained_butterflies = counts.retained_butterflies;
            }
            else
            {
                const leopard::ff16::TestOnlySparseEncodeCounts counts =
                    leopard::ff16::TestOnlyGetSparseEncodeCounts();
                exact_blocks = counts.exact_blocks;
                prefix_butterflies = counts.prefix_butterflies;
                retained_butterflies = counts.retained_butterflies;
            }
            require(exact_blocks != 0 &&
                    retained_butterflies < prefix_butterflies,
                "forced sparse schedule did not reduce backend work");
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
            require_result(leo2_test_codec_set_encode_mode(
                codec, LEO2_TEST_ENCODE_FORCE_TRANSFORM),
                LEO2_SUCCESS, "force shared transform");
            const Shards sparse_reference = execute_encode(
                codec, test_case, originals, true);
            leo2_decode_plan* plan = make_plan(codec, test_case);
            std::atomic<unsigned> ready(0);
            std::atomic<bool> go(false);
            std::atomic<unsigned> failures(0);
            std::vector<std::thread> threads;
            for (unsigned thread_i = 0; thread_i < thread_count; ++thread_i)
            {
                threads.push_back(std::thread(
                    [codec, plan, test_case, &originals, &reference,
                     &sparse_reference, &ready, &go, &failures]() {
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
                            const Shards sparse = execute_encode(
                                codec, test_case, originals, true);
                            if (sparse != sparse_reference)
                                throw std::runtime_error(
                                    "shared codec sparse parity mismatch");
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
        test_sparse_encode_contexts(contexts);
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
