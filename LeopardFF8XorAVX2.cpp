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
#include "LeopardFF8XorAVX2.h"

#if defined(LEO_HAS_FF8) && defined(LEO_FF8XOR_BUILD_AVX2)

#include <immintrin.h>

namespace leopard { namespace ff8xor {


// The generated source resolves this overload at template instantiation time.
// Every circuit wire remains a statically named YMM value; there is no indexed
// gate interpreter or register array in the payload loop.
static LEO_FORCE_INLINE __m256i XorValue(__m256i x, __m256i y)
{
    return _mm256_xor_si256(x, y);
}


}} // namespace leopard::ff8xor

#include "generated/LeopardFF8XorCircuits.inl"

namespace leopard { namespace ff8xor { namespace avx2 {


static const unsigned kVectorBytes = 32;

static LEO_FORCE_INLINE __m256i Load(const uint8_t* source)
{
    return _mm256_loadu_si256(reinterpret_cast<const __m256i*>(source));
}

static LEO_FORCE_INLINE void Store(uint8_t* destination, __m256i value)
{
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(destination), value);
}

template <unsigned Coefficient>
static LEO_FORCE_INLINE void MultiplyChunk(
    uint8_t* destination,
    const uint8_t* source,
    uint64_t plane_bytes,
    uint64_t offset)
{
    __m256i x0 = Load(source + offset);
    __m256i x1 = Load(source + plane_bytes + offset);
    __m256i x2 = Load(source + plane_bytes * 2 + offset);
    __m256i x3 = Load(source + plane_bytes * 3 + offset);
    __m256i x4 = Load(source + plane_bytes * 4 + offset);
    __m256i x5 = Load(source + plane_bytes * 5 + offset);
    __m256i x6 = Load(source + plane_bytes * 6 + offset);
    __m256i x7 = Load(source + plane_bytes * 7 + offset);

    generated::MultiplyCircuit<Coefficient>::Apply(
        x0, x1, x2, x3, x4, x5, x6, x7);

    Store(destination + offset, x0);
    Store(destination + plane_bytes + offset, x1);
    Store(destination + plane_bytes * 2 + offset, x2);
    Store(destination + plane_bytes * 3 + offset, x3);
    Store(destination + plane_bytes * 4 + offset, x4);
    Store(destination + plane_bytes * 5 + offset, x5);
    Store(destination + plane_bytes * 6 + offset, x6);
    Store(destination + plane_bytes * 7 + offset, x7);
}

template <unsigned Skew, bool Inverse>
static LEO_FORCE_INLINE void ButterflyChunk(
    uint8_t* x_buffer,
    uint8_t* y_buffer,
    uint64_t plane_bytes,
    uint64_t offset)
{
    __m256i x0 = Load(x_buffer + offset);
    __m256i x1 = Load(x_buffer + plane_bytes + offset);
    __m256i x2 = Load(x_buffer + plane_bytes * 2 + offset);
    __m256i x3 = Load(x_buffer + plane_bytes * 3 + offset);
    __m256i x4 = Load(x_buffer + plane_bytes * 4 + offset);
    __m256i x5 = Load(x_buffer + plane_bytes * 5 + offset);
    __m256i x6 = Load(x_buffer + plane_bytes * 6 + offset);
    __m256i x7 = Load(x_buffer + plane_bytes * 7 + offset);
    __m256i y0 = Load(y_buffer + offset);
    __m256i y1 = Load(y_buffer + plane_bytes + offset);
    __m256i y2 = Load(y_buffer + plane_bytes * 2 + offset);
    __m256i y3 = Load(y_buffer + plane_bytes * 3 + offset);
    __m256i y4 = Load(y_buffer + plane_bytes * 4 + offset);
    __m256i y5 = Load(y_buffer + plane_bytes * 5 + offset);
    __m256i y6 = Load(y_buffer + plane_bytes * 6 + offset);
    __m256i y7 = Load(y_buffer + plane_bytes * 7 + offset);

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

    Store(x_buffer + offset, x0);
    Store(x_buffer + plane_bytes + offset, x1);
    Store(x_buffer + plane_bytes * 2 + offset, x2);
    Store(x_buffer + plane_bytes * 3 + offset, x3);
    Store(x_buffer + plane_bytes * 4 + offset, x4);
    Store(x_buffer + plane_bytes * 5 + offset, x5);
    Store(x_buffer + plane_bytes * 6 + offset, x6);
    Store(x_buffer + plane_bytes * 7 + offset, x7);
    Store(y_buffer + offset, y0);
    Store(y_buffer + plane_bytes + offset, y1);
    Store(y_buffer + plane_bytes * 2 + offset, y2);
    Store(y_buffer + plane_bytes * 3 + offset, y3);
    Store(y_buffer + plane_bytes * 4 + offset, y4);
    Store(y_buffer + plane_bytes * 5 + offset, y5);
    Store(y_buffer + plane_bytes * 6 + offset, y6);
    Store(y_buffer + plane_bytes * 7 + offset, y7);
}

bool KernelsBuilt()
{
    return true;
}

template <unsigned Coefficient>
uint64_t Multiply(
    void* destination_void,
    const void* source_void,
    uint64_t buffer_bytes,
    uint64_t offset)
{
    uint8_t* destination = reinterpret_cast<uint8_t*>(destination_void);
    const uint8_t* source = reinterpret_cast<const uint8_t*>(source_void);
    const uint64_t plane_bytes = buffer_bytes >> 3;

    while (plane_bytes - offset >= kVectorBytes)
    {
        MultiplyChunk<Coefficient>(
            destination, source, plane_bytes, offset);
        offset += kVectorBytes;
    }
    return offset;
}

template <unsigned Skew, bool Inverse>
uint64_t Butterfly(
    void* x_void,
    void* y_void,
    uint64_t buffer_bytes,
    uint64_t offset)
{
    uint8_t* x_buffer = reinterpret_cast<uint8_t*>(x_void);
    uint8_t* y_buffer = reinterpret_cast<uint8_t*>(y_void);
    const uint64_t plane_bytes = buffer_bytes >> 3;

    while (plane_bytes - offset >= kVectorBytes)
    {
        ButterflyChunk<Skew, Inverse>(
            x_buffer, y_buffer, plane_bytes, offset);
        offset += kVectorBytes;
    }
    return offset;
}

// The baseline object specializes once per coefficient/skew and calls these
// symbols directly.  Explicit instantiation preserves whole-buffer dispatch
// while keeping all AVX2 instructions in this isolated translation unit.
#define LEO_FF8XOR_AVX2_INSTANTIATE(Value) \
    template uint64_t Multiply<Value>( \
        void*, const void*, uint64_t, uint64_t); \
    template uint64_t Butterfly<Value, false>( \
        void*, void*, uint64_t, uint64_t); \
    template uint64_t Butterfly<Value, true>( \
        void*, void*, uint64_t, uint64_t);
LEOPARD_FF8XOR_FOR_EACH_LOG(LEO_FF8XOR_AVX2_INSTANTIATE)
#undef LEO_FF8XOR_AVX2_INSTANTIATE

uint64_t Xor(
    void* destination_void,
    const void* source_void,
    uint64_t buffer_bytes,
    uint64_t offset)
{
    uint8_t* destination = reinterpret_cast<uint8_t*>(destination_void);
    const uint8_t* source = reinterpret_cast<const uint8_t*>(source_void);

    while (buffer_bytes - offset >= kVectorBytes)
    {
        Store(destination + offset,
            _mm256_xor_si256(
                Load(destination + offset),
                Load(source + offset)));
        offset += kVectorBytes;
    }
    return offset;
}


}}} // namespace leopard::ff8xor::avx2

#else

namespace leopard { namespace ff8xor { namespace avx2 {


bool KernelsBuilt()
{
    return false;
}

uint64_t Xor(void*, const void*, uint64_t, uint64_t start_offset)
{
    return start_offset;
}


}}} // namespace leopard::ff8xor::avx2

#endif
