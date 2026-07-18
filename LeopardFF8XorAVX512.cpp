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

#include "LeopardCommon.h"
#include "LeopardFF8XorAVX512.h"

#if defined(LEO_HAS_FF8) && defined(LEO_FF8XOR_BUILD_AVX512)

#include <immintrin.h>

namespace leopard { namespace ff8xor {


// The generated source resolves XorValue at template instantiation time.  A
// distinct overload for each width keeps every circuit wire in a statically
// named vector value; there is no dynamically indexed gate interpreter.
static LEO_FORCE_INLINE __m256i XorValue(__m256i x, __m256i y)
{
    return _mm256_xor_si256(x, y);
}

static LEO_FORCE_INLINE __m512i XorValue(__m512i x, __m512i y)
{
    return _mm512_xor_si512(x, y);
}


}} // namespace leopard::ff8xor

#include "generated/LeopardFF8XorCircuits.inl"

namespace leopard { namespace ff8xor { namespace avx512 {


struct YmmTag {};
struct ZmmTag {};

template <typename Tag>
struct ValueIO;

template <>
struct ValueIO<YmmTag>
{
    typedef __m256i Value;
    static const unsigned kBytes = 32;

    static LEO_FORCE_INLINE Value Load(const uint8_t* source)
    {
        return _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(source));
    }

    static LEO_FORCE_INLINE void Store(uint8_t* destination, Value value)
    {
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(destination), value);
    }
};

template <>
struct ValueIO<ZmmTag>
{
    typedef __m512i Value;
    static const unsigned kBytes = 64;

    static LEO_FORCE_INLINE Value Load(const uint8_t* source)
    {
        return _mm512_loadu_si512(
            reinterpret_cast<const void*>(source));
    }

    static LEO_FORCE_INLINE void Store(uint8_t* destination, Value value)
    {
        _mm512_storeu_si512(reinterpret_cast<void*>(destination), value);
    }
};

template <unsigned Coefficient, typename Tag>
static LEO_FORCE_INLINE void MultiplyChunk(
    uint8_t* destination,
    const uint8_t* source,
    uint64_t plane_bytes,
    uint64_t offset)
{
    typedef typename ValueIO<Tag>::Value Value;
    Value x0 = ValueIO<Tag>::Load(source + offset);
    Value x1 = ValueIO<Tag>::Load(source + plane_bytes + offset);
    Value x2 = ValueIO<Tag>::Load(source + plane_bytes * 2 + offset);
    Value x3 = ValueIO<Tag>::Load(source + plane_bytes * 3 + offset);
    Value x4 = ValueIO<Tag>::Load(source + plane_bytes * 4 + offset);
    Value x5 = ValueIO<Tag>::Load(source + plane_bytes * 5 + offset);
    Value x6 = ValueIO<Tag>::Load(source + plane_bytes * 6 + offset);
    Value x7 = ValueIO<Tag>::Load(source + plane_bytes * 7 + offset);

    generated::MultiplyCircuit<Coefficient>::Apply(
        x0, x1, x2, x3, x4, x5, x6, x7);

    ValueIO<Tag>::Store(destination + offset, x0);
    ValueIO<Tag>::Store(destination + plane_bytes + offset, x1);
    ValueIO<Tag>::Store(destination + plane_bytes * 2 + offset, x2);
    ValueIO<Tag>::Store(destination + plane_bytes * 3 + offset, x3);
    ValueIO<Tag>::Store(destination + plane_bytes * 4 + offset, x4);
    ValueIO<Tag>::Store(destination + plane_bytes * 5 + offset, x5);
    ValueIO<Tag>::Store(destination + plane_bytes * 6 + offset, x6);
    ValueIO<Tag>::Store(destination + plane_bytes * 7 + offset, x7);
}

template <unsigned Skew, typename Tag, bool Inverse>
static LEO_FORCE_INLINE void ButterflyChunk(
    uint8_t* x_buffer,
    uint8_t* y_buffer,
    uint64_t plane_bytes,
    uint64_t offset)
{
    typedef typename ValueIO<Tag>::Value Value;
    Value x0 = ValueIO<Tag>::Load(x_buffer + offset);
    Value x1 = ValueIO<Tag>::Load(x_buffer + plane_bytes + offset);
    Value x2 = ValueIO<Tag>::Load(x_buffer + plane_bytes * 2 + offset);
    Value x3 = ValueIO<Tag>::Load(x_buffer + plane_bytes * 3 + offset);
    Value x4 = ValueIO<Tag>::Load(x_buffer + plane_bytes * 4 + offset);
    Value x5 = ValueIO<Tag>::Load(x_buffer + plane_bytes * 5 + offset);
    Value x6 = ValueIO<Tag>::Load(x_buffer + plane_bytes * 6 + offset);
    Value x7 = ValueIO<Tag>::Load(x_buffer + plane_bytes * 7 + offset);
    Value y0 = ValueIO<Tag>::Load(y_buffer + offset);
    Value y1 = ValueIO<Tag>::Load(y_buffer + plane_bytes + offset);
    Value y2 = ValueIO<Tag>::Load(y_buffer + plane_bytes * 2 + offset);
    Value y3 = ValueIO<Tag>::Load(y_buffer + plane_bytes * 3 + offset);
    Value y4 = ValueIO<Tag>::Load(y_buffer + plane_bytes * 4 + offset);
    Value y5 = ValueIO<Tag>::Load(y_buffer + plane_bytes * 5 + offset);
    Value y6 = ValueIO<Tag>::Load(y_buffer + plane_bytes * 6 + offset);
    Value y7 = ValueIO<Tag>::Load(y_buffer + plane_bytes * 7 + offset);

    if (Inverse)
    {
        generated::IFFTCircuit<Skew>::Apply(
            x0, x1, x2, x3, x4, x5, x6, x7,
            y0, y1, y2, y3, y4, y5, y6, y7);
    }
    else
    {
        generated::FFTCircuit<Skew>::Apply(
            x0, x1, x2, x3, x4, x5, x6, x7,
            y0, y1, y2, y3, y4, y5, y6, y7);
    }

    ValueIO<Tag>::Store(x_buffer + offset, x0);
    ValueIO<Tag>::Store(x_buffer + plane_bytes + offset, x1);
    ValueIO<Tag>::Store(x_buffer + plane_bytes * 2 + offset, x2);
    ValueIO<Tag>::Store(x_buffer + plane_bytes * 3 + offset, x3);
    ValueIO<Tag>::Store(x_buffer + plane_bytes * 4 + offset, x4);
    ValueIO<Tag>::Store(x_buffer + plane_bytes * 5 + offset, x5);
    ValueIO<Tag>::Store(x_buffer + plane_bytes * 6 + offset, x6);
    ValueIO<Tag>::Store(x_buffer + plane_bytes * 7 + offset, x7);
    ValueIO<Tag>::Store(y_buffer + offset, y0);
    ValueIO<Tag>::Store(y_buffer + plane_bytes + offset, y1);
    ValueIO<Tag>::Store(y_buffer + plane_bytes * 2 + offset, y2);
    ValueIO<Tag>::Store(y_buffer + plane_bytes * 3 + offset, y3);
    ValueIO<Tag>::Store(y_buffer + plane_bytes * 4 + offset, y4);
    ValueIO<Tag>::Store(y_buffer + plane_bytes * 5 + offset, y5);
    ValueIO<Tag>::Store(y_buffer + plane_bytes * 6 + offset, y6);
    ValueIO<Tag>::Store(y_buffer + plane_bytes * 7 + offset, y7);
}

template <unsigned Coefficient, typename Tag>
static uint64_t MultiplyWhole(
    void* destination_void,
    const void* source_void,
    uint64_t buffer_bytes)
{
    uint8_t* destination = reinterpret_cast<uint8_t*>(destination_void);
    const uint8_t* source = reinterpret_cast<const uint8_t*>(source_void);
    const uint64_t plane_bytes = buffer_bytes >> 3;
    uint64_t offset = 0;

    while (plane_bytes - offset >= ValueIO<Tag>::kBytes)
    {
        MultiplyChunk<Coefficient, Tag>(
            destination, source, plane_bytes, offset);
        offset += ValueIO<Tag>::kBytes;
    }
    return offset;
}

template <unsigned Skew, typename Tag, bool Inverse>
static uint64_t ButterflyWhole(
    void* x_void,
    void* y_void,
    uint64_t buffer_bytes)
{
    uint8_t* x_buffer = reinterpret_cast<uint8_t*>(x_void);
    uint8_t* y_buffer = reinterpret_cast<uint8_t*>(y_void);
    const uint64_t plane_bytes = buffer_bytes >> 3;
    uint64_t offset = 0;

    while (plane_bytes - offset >= ValueIO<Tag>::kBytes)
    {
        ButterflyChunk<Skew, Tag, Inverse>(
            x_buffer, y_buffer, plane_bytes, offset);
        offset += ValueIO<Tag>::kBytes;
    }
    return offset;
}

template <typename Tag>
static uint64_t XorWhole(
    void* destination_void,
    const void* source_void,
    uint64_t buffer_bytes)
{
    uint8_t* destination = reinterpret_cast<uint8_t*>(destination_void);
    const uint8_t* source = reinterpret_cast<const uint8_t*>(source_void);
    uint64_t offset = 0;

    while (buffer_bytes - offset >= ValueIO<Tag>::kBytes)
    {
        const typename ValueIO<Tag>::Value result = XorValue(
            ValueIO<Tag>::Load(destination + offset),
            ValueIO<Tag>::Load(source + offset));
        ValueIO<Tag>::Store(destination + offset, result);
        offset += ValueIO<Tag>::kBytes;
    }
    return offset;
}

template <typename Tag>
static uint64_t Xor2Whole(
    void* destination0_void,
    const void* source0_void,
    void* destination1_void,
    const void* source1_void,
    uint64_t buffer_bytes)
{
    uint8_t* destination0 = reinterpret_cast<uint8_t*>(destination0_void);
    const uint8_t* source0 = reinterpret_cast<const uint8_t*>(source0_void);
    uint8_t* destination1 = reinterpret_cast<uint8_t*>(destination1_void);
    const uint8_t* source1 = reinterpret_cast<const uint8_t*>(source1_void);
    uint64_t offset = 0;

    while (buffer_bytes - offset >= ValueIO<Tag>::kBytes)
    {
        const typename ValueIO<Tag>::Value result0 = XorValue(
            ValueIO<Tag>::Load(destination0 + offset),
            ValueIO<Tag>::Load(source0 + offset));
        const typename ValueIO<Tag>::Value result1 = XorValue(
            ValueIO<Tag>::Load(destination1 + offset),
            ValueIO<Tag>::Load(source1 + offset));
        ValueIO<Tag>::Store(destination0 + offset, result0);
        ValueIO<Tag>::Store(destination1 + offset, result1);
        offset += ValueIO<Tag>::kBytes;
    }
    return offset;
}

template <typename Tag>
static uint64_t Xor4Whole(
    void* destination0_void,
    const void* source0_void,
    void* destination1_void,
    const void* source1_void,
    void* destination2_void,
    const void* source2_void,
    void* destination3_void,
    const void* source3_void,
    uint64_t buffer_bytes)
{
    uint8_t* destination0 = reinterpret_cast<uint8_t*>(destination0_void);
    const uint8_t* source0 = reinterpret_cast<const uint8_t*>(source0_void);
    uint8_t* destination1 = reinterpret_cast<uint8_t*>(destination1_void);
    const uint8_t* source1 = reinterpret_cast<const uint8_t*>(source1_void);
    uint8_t* destination2 = reinterpret_cast<uint8_t*>(destination2_void);
    const uint8_t* source2 = reinterpret_cast<const uint8_t*>(source2_void);
    uint8_t* destination3 = reinterpret_cast<uint8_t*>(destination3_void);
    const uint8_t* source3 = reinterpret_cast<const uint8_t*>(source3_void);
    uint64_t offset = 0;

    while (buffer_bytes - offset >= ValueIO<Tag>::kBytes)
    {
        const typename ValueIO<Tag>::Value result0 = XorValue(
            ValueIO<Tag>::Load(destination0 + offset),
            ValueIO<Tag>::Load(source0 + offset));
        const typename ValueIO<Tag>::Value result1 = XorValue(
            ValueIO<Tag>::Load(destination1 + offset),
            ValueIO<Tag>::Load(source1 + offset));
        const typename ValueIO<Tag>::Value result2 = XorValue(
            ValueIO<Tag>::Load(destination2 + offset),
            ValueIO<Tag>::Load(source2 + offset));
        const typename ValueIO<Tag>::Value result3 = XorValue(
            ValueIO<Tag>::Load(destination3 + offset),
            ValueIO<Tag>::Load(source3 + offset));
        ValueIO<Tag>::Store(destination0 + offset, result0);
        ValueIO<Tag>::Store(destination1 + offset, result1);
        ValueIO<Tag>::Store(destination2 + offset, result2);
        ValueIO<Tag>::Store(destination3 + offset, result3);
        offset += ValueIO<Tag>::kBytes;
    }
    return offset;
}

bool KernelsBuilt()
{
    return true;
}

template <unsigned Coefficient>
uint64_t Multiply256(
    void* destination,
    const void* source,
    uint64_t buffer_bytes)
{
    return MultiplyWhole<Coefficient, YmmTag>(
        destination, source, buffer_bytes);
}

template <unsigned Coefficient>
uint64_t Multiply512(
    void* destination,
    const void* source,
    uint64_t buffer_bytes)
{
    return MultiplyWhole<Coefficient, ZmmTag>(
        destination, source, buffer_bytes);
}

template <unsigned Skew, bool Inverse>
uint64_t Butterfly256(
    void* x,
    void* y,
    uint64_t buffer_bytes)
{
    return ButterflyWhole<Skew, YmmTag, Inverse>(x, y, buffer_bytes);
}

template <unsigned Skew, bool Inverse>
uint64_t Butterfly512(
    void* x,
    void* y,
    uint64_t buffer_bytes)
{
    return ButterflyWhole<Skew, ZmmTag, Inverse>(x, y, buffer_bytes);
}

// The baseline object already specializes once per coefficient/skew and calls
// these symbols directly.  Explicit instantiation keeps implementation details
// isolated without adding a second runtime function-table dispatch.
#define LEO_FF8XOR_AVX512_INSTANTIATE(Value) \
    template uint64_t Multiply256<Value>( \
        void*, const void*, uint64_t); \
    template uint64_t Multiply512<Value>( \
        void*, const void*, uint64_t); \
    template uint64_t Butterfly256<Value, false>( \
        void*, void*, uint64_t); \
    template uint64_t Butterfly256<Value, true>( \
        void*, void*, uint64_t); \
    template uint64_t Butterfly512<Value, false>( \
        void*, void*, uint64_t); \
    template uint64_t Butterfly512<Value, true>( \
        void*, void*, uint64_t);
LEOPARD_FF8XOR_FOR_EACH_LOG(LEO_FF8XOR_AVX512_INSTANTIATE)
#undef LEO_FF8XOR_AVX512_INSTANTIATE

uint64_t Xor256(
    void* destination,
    const void* source,
    uint64_t buffer_bytes)
{
    return XorWhole<YmmTag>(destination, source, buffer_bytes);
}

uint64_t Xor512(
    void* destination,
    const void* source,
    uint64_t buffer_bytes)
{
    return XorWhole<ZmmTag>(destination, source, buffer_bytes);
}

uint64_t Xor2_256(
    void* destination0,
    const void* source0,
    void* destination1,
    const void* source1,
    uint64_t buffer_bytes)
{
    return Xor2Whole<YmmTag>(
        destination0, source0, destination1, source1, buffer_bytes);
}

uint64_t Xor2_512(
    void* destination0,
    const void* source0,
    void* destination1,
    const void* source1,
    uint64_t buffer_bytes)
{
    return Xor2Whole<ZmmTag>(
        destination0, source0, destination1, source1, buffer_bytes);
}

uint64_t Xor4_256(
    void* destination0,
    const void* source0,
    void* destination1,
    const void* source1,
    void* destination2,
    const void* source2,
    void* destination3,
    const void* source3,
    uint64_t buffer_bytes)
{
    return Xor4Whole<YmmTag>(
        destination0, source0, destination1, source1,
        destination2, source2, destination3, source3, buffer_bytes);
}

uint64_t Xor4_512(
    void* destination0,
    const void* source0,
    void* destination1,
    const void* source1,
    void* destination2,
    const void* source2,
    void* destination3,
    const void* source3,
    uint64_t buffer_bytes)
{
    return Xor4Whole<ZmmTag>(
        destination0, source0, destination1, source1,
        destination2, source2, destination3, source3, buffer_bytes);
}


}}} // namespace leopard::ff8xor::avx512

#else

namespace leopard { namespace ff8xor { namespace avx512 {


bool KernelsBuilt()
{
    return false;
}

uint64_t Xor256(void*, const void*, uint64_t)
{
    return 0;
}

uint64_t Xor512(void*, const void*, uint64_t)
{
    return 0;
}

uint64_t Xor2_256(void*, const void*, void*, const void*, uint64_t)
{
    return 0;
}

uint64_t Xor2_512(void*, const void*, void*, const void*, uint64_t)
{
    return 0;
}

uint64_t Xor4_256(
    void*, const void*, void*, const void*, void*, const void*,
    void*, const void*, uint64_t)
{
    return 0;
}

uint64_t Xor4_512(
    void*, const void*, void*, const void*, void*, const void*,
    void*, const void*, uint64_t)
{
    return 0;
}


}}} // namespace leopard::ff8xor::avx512

#endif
