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
#include "LeopardFF8XorAVX512Four.h"

#if defined(LEO_HAS_FF8) && defined(LEO_FF8XOR_BUILD_AVX512_FOUR)

#include <immintrin.h>

namespace leopard { namespace ff8xor {


// Tied destructive operands are intentional.  With all 32 architectural ZMM
// registers occupied by named planes, allowing an SSA result to choose a 33rd
// register turns a single XOR into a stack spill.  The +v constraint makes the
// destination both an input and output while retaining literal generated code.
static LEO_FORCE_INLINE __m512i XorValue(__m512i x, __m512i y)
{
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__(
        "vpxord %1, %0, %0"
        : "+v"(x)
        : "v"(y));
    return x;
#else
    return _mm512_xor_si512(x, y);
#endif
}

static LEO_FORCE_INLINE __m512i Xor3Value(
    __m512i x,
    __m512i y,
    __m512i z)
{
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__(
        "vpternlogd $0x96, %2, %1, %0"
        : "+v"(x)
        : "v"(y), "v"(z));
    return x;
#else
    return _mm512_ternarylogic_epi32(x, y, z, 0x96);
#endif
}


}} // namespace leopard::ff8xor

#include "generated/LeopardFF8XorFourBufferCircuits.inl"

namespace leopard { namespace ff8xor { namespace avx512four {


template <unsigned TupleIndex, bool Inverse, bool UseXor3>
static void ApplyWhole(
    void* a_void,
    void* b_void,
    void* c_void,
    void* d_void,
    uint64_t buffer_bytes)
{
    uint8_t* a_buffer = reinterpret_cast<uint8_t*>(a_void);
    uint8_t* b_buffer = reinterpret_cast<uint8_t*>(b_void);
    uint8_t* c_buffer = reinterpret_cast<uint8_t*>(c_void);
    uint8_t* d_buffer = reinterpret_cast<uint8_t*>(d_void);
    const uint64_t plane_bytes = buffer_bytes >> 3;

    for (uint64_t offset = 0; offset < plane_bytes; offset += 64)
    {
        __m512i a0 = _mm512_loadu_si512(a_buffer + offset);
        __m512i a1 = _mm512_loadu_si512(a_buffer + plane_bytes + offset);
        __m512i a2 = _mm512_loadu_si512(a_buffer + plane_bytes * 2 + offset);
        __m512i a3 = _mm512_loadu_si512(a_buffer + plane_bytes * 3 + offset);
        __m512i a4 = _mm512_loadu_si512(a_buffer + plane_bytes * 4 + offset);
        __m512i a5 = _mm512_loadu_si512(a_buffer + plane_bytes * 5 + offset);
        __m512i a6 = _mm512_loadu_si512(a_buffer + plane_bytes * 6 + offset);
        __m512i a7 = _mm512_loadu_si512(a_buffer + plane_bytes * 7 + offset);
        __m512i b0 = _mm512_loadu_si512(b_buffer + offset);
        __m512i b1 = _mm512_loadu_si512(b_buffer + plane_bytes + offset);
        __m512i b2 = _mm512_loadu_si512(b_buffer + plane_bytes * 2 + offset);
        __m512i b3 = _mm512_loadu_si512(b_buffer + plane_bytes * 3 + offset);
        __m512i b4 = _mm512_loadu_si512(b_buffer + plane_bytes * 4 + offset);
        __m512i b5 = _mm512_loadu_si512(b_buffer + plane_bytes * 5 + offset);
        __m512i b6 = _mm512_loadu_si512(b_buffer + plane_bytes * 6 + offset);
        __m512i b7 = _mm512_loadu_si512(b_buffer + plane_bytes * 7 + offset);
        __m512i c0 = _mm512_loadu_si512(c_buffer + offset);
        __m512i c1 = _mm512_loadu_si512(c_buffer + plane_bytes + offset);
        __m512i c2 = _mm512_loadu_si512(c_buffer + plane_bytes * 2 + offset);
        __m512i c3 = _mm512_loadu_si512(c_buffer + plane_bytes * 3 + offset);
        __m512i c4 = _mm512_loadu_si512(c_buffer + plane_bytes * 4 + offset);
        __m512i c5 = _mm512_loadu_si512(c_buffer + plane_bytes * 5 + offset);
        __m512i c6 = _mm512_loadu_si512(c_buffer + plane_bytes * 6 + offset);
        __m512i c7 = _mm512_loadu_si512(c_buffer + plane_bytes * 7 + offset);
        __m512i d0 = _mm512_loadu_si512(d_buffer + offset);
        __m512i d1 = _mm512_loadu_si512(d_buffer + plane_bytes + offset);
        __m512i d2 = _mm512_loadu_si512(d_buffer + plane_bytes * 2 + offset);
        __m512i d3 = _mm512_loadu_si512(d_buffer + plane_bytes * 3 + offset);
        __m512i d4 = _mm512_loadu_si512(d_buffer + plane_bytes * 4 + offset);
        __m512i d5 = _mm512_loadu_si512(d_buffer + plane_bytes * 5 + offset);
        __m512i d6 = _mm512_loadu_si512(d_buffer + plane_bytes * 6 + offset);
        __m512i d7 = _mm512_loadu_si512(d_buffer + plane_bytes * 7 + offset);

#if defined(__GNUC__) || defined(__clang__)
        // Break the compiler's equivalence between each loaded value and its
        // backing payload address.  Without these zero-instruction tied
        // barriers GCC may rematerialize a late circuit source by reading the
        // shard again.  Four groups stay below GCC's 30-operand asm limit
        // (each read/write operand counts twice); all 32 opaque values are then
        // live in the 32 architectural ZMM registers.
        __asm__ __volatile__(
            ""
            : "+v"(a0), "+v"(a1), "+v"(a2), "+v"(a3),
              "+v"(a4), "+v"(a5), "+v"(a6), "+v"(a7));
        __asm__ __volatile__(
            ""
            : "+v"(b0), "+v"(b1), "+v"(b2), "+v"(b3),
              "+v"(b4), "+v"(b5), "+v"(b6), "+v"(b7));
        __asm__ __volatile__(
            ""
            : "+v"(c0), "+v"(c1), "+v"(c2), "+v"(c3),
              "+v"(c4), "+v"(c5), "+v"(c6), "+v"(c7));
        __asm__ __volatile__(
            ""
            : "+v"(d0), "+v"(d1), "+v"(d2), "+v"(d3),
              "+v"(d4), "+v"(d5), "+v"(d6), "+v"(d7));
#endif

        if (Inverse)
        {
            if (UseXor3)
            {
                generated4::IFFT4CircuitXor3<TupleIndex>::Apply(
                    a0, a1, a2, a3, a4, a5, a6, a7,
                    b0, b1, b2, b3, b4, b5, b6, b7,
                    c0, c1, c2, c3, c4, c5, c6, c7,
                    d0, d1, d2, d3, d4, d5, d6, d7);
            }
            else
            {
                generated4::IFFT4CircuitXor2<TupleIndex>::Apply(
                    a0, a1, a2, a3, a4, a5, a6, a7,
                    b0, b1, b2, b3, b4, b5, b6, b7,
                    c0, c1, c2, c3, c4, c5, c6, c7,
                    d0, d1, d2, d3, d4, d5, d6, d7);
            }
        }
        else
        {
            if (UseXor3)
            {
                generated4::FFT4CircuitXor3<TupleIndex>::Apply(
                    a0, a1, a2, a3, a4, a5, a6, a7,
                    b0, b1, b2, b3, b4, b5, b6, b7,
                    c0, c1, c2, c3, c4, c5, c6, c7,
                    d0, d1, d2, d3, d4, d5, d6, d7);
            }
            else
            {
                generated4::FFT4CircuitXor2<TupleIndex>::Apply(
                    a0, a1, a2, a3, a4, a5, a6, a7,
                    b0, b1, b2, b3, b4, b5, b6, b7,
                    c0, c1, c2, c3, c4, c5, c6, c7,
                    d0, d1, d2, d3, d4, d5, d6, d7);
            }
        }

        _mm512_storeu_si512(a_buffer + offset, a0);
        _mm512_storeu_si512(a_buffer + plane_bytes + offset, a1);
        _mm512_storeu_si512(a_buffer + plane_bytes * 2 + offset, a2);
        _mm512_storeu_si512(a_buffer + plane_bytes * 3 + offset, a3);
        _mm512_storeu_si512(a_buffer + plane_bytes * 4 + offset, a4);
        _mm512_storeu_si512(a_buffer + plane_bytes * 5 + offset, a5);
        _mm512_storeu_si512(a_buffer + plane_bytes * 6 + offset, a6);
        _mm512_storeu_si512(a_buffer + plane_bytes * 7 + offset, a7);
        _mm512_storeu_si512(b_buffer + offset, b0);
        _mm512_storeu_si512(b_buffer + plane_bytes + offset, b1);
        _mm512_storeu_si512(b_buffer + plane_bytes * 2 + offset, b2);
        _mm512_storeu_si512(b_buffer + plane_bytes * 3 + offset, b3);
        _mm512_storeu_si512(b_buffer + plane_bytes * 4 + offset, b4);
        _mm512_storeu_si512(b_buffer + plane_bytes * 5 + offset, b5);
        _mm512_storeu_si512(b_buffer + plane_bytes * 6 + offset, b6);
        _mm512_storeu_si512(b_buffer + plane_bytes * 7 + offset, b7);
        _mm512_storeu_si512(c_buffer + offset, c0);
        _mm512_storeu_si512(c_buffer + plane_bytes + offset, c1);
        _mm512_storeu_si512(c_buffer + plane_bytes * 2 + offset, c2);
        _mm512_storeu_si512(c_buffer + plane_bytes * 3 + offset, c3);
        _mm512_storeu_si512(c_buffer + plane_bytes * 4 + offset, c4);
        _mm512_storeu_si512(c_buffer + plane_bytes * 5 + offset, c5);
        _mm512_storeu_si512(c_buffer + plane_bytes * 6 + offset, c6);
        _mm512_storeu_si512(c_buffer + plane_bytes * 7 + offset, c7);
        _mm512_storeu_si512(d_buffer + offset, d0);
        _mm512_storeu_si512(d_buffer + plane_bytes + offset, d1);
        _mm512_storeu_si512(d_buffer + plane_bytes * 2 + offset, d2);
        _mm512_storeu_si512(d_buffer + plane_bytes * 3 + offset, d3);
        _mm512_storeu_si512(d_buffer + plane_bytes * 4 + offset, d4);
        _mm512_storeu_si512(d_buffer + plane_bytes * 5 + offset, d5);
        _mm512_storeu_si512(d_buffer + plane_bytes * 6 + offset, d6);
        _mm512_storeu_si512(d_buffer + plane_bytes * 7 + offset, d7);
    }
}

typedef void (*FourBufferFunction)(
    void*, void*, void*, void*, uint64_t);

#define LEO_FF8XOR_FOUR_FUNCTION(Index, Skew01, Skew23, Skew02) \
    &ApplyWhole<Index, false, false>,
static const FourBufferFunction FFTXor2Functions[64] = {
    LEOPARD_FF8XOR_FOR_EACH_FOUR_TUPLE(LEO_FF8XOR_FOUR_FUNCTION)
};
#undef LEO_FF8XOR_FOUR_FUNCTION

#define LEO_FF8XOR_FOUR_FUNCTION(Index, Skew01, Skew23, Skew02) \
    &ApplyWhole<Index, false, true>,
static const FourBufferFunction FFTXor3Functions[64] = {
    LEOPARD_FF8XOR_FOR_EACH_FOUR_TUPLE(LEO_FF8XOR_FOUR_FUNCTION)
};
#undef LEO_FF8XOR_FOUR_FUNCTION

#define LEO_FF8XOR_FOUR_FUNCTION(Index, Skew01, Skew23, Skew02) \
    &ApplyWhole<Index, true, false>,
static const FourBufferFunction IFFTXor2Functions[64] = {
    LEOPARD_FF8XOR_FOR_EACH_FOUR_TUPLE(LEO_FF8XOR_FOUR_FUNCTION)
};
#undef LEO_FF8XOR_FOUR_FUNCTION

#define LEO_FF8XOR_FOUR_FUNCTION(Index, Skew01, Skew23, Skew02) \
    &ApplyWhole<Index, true, true>,
static const FourBufferFunction IFFTXor3Functions[64] = {
    LEOPARD_FF8XOR_FOR_EACH_FOUR_TUPLE(LEO_FF8XOR_FOUR_FUNCTION)
};
#undef LEO_FF8XOR_FOUR_FUNCTION

bool KernelsBuilt()
{
    return true;
}

unsigned GetTupleCount()
{
    return 64;
}

bool GetTuple(
    unsigned tuple_index,
    uint8_t& skew01,
    uint8_t& skew23,
    uint8_t& skew02)
{
    if (tuple_index >= 64)
        return false;
    const generated4::FourBufferSkewTuple& tuple =
        generated4::kFourBufferSkewTuples[tuple_index];
    skew01 = tuple.Skew01;
    skew23 = tuple.Skew23;
    skew02 = tuple.Skew02;
    return true;
}

int FindTupleIndex(uint8_t skew01, uint8_t skew23, uint8_t skew02)
{
    switch (skew01)
    {
#define LEO_FF8XOR_FOUR_CASE(Index, Skew01, Skew23, Skew02) \
    case Skew01: \
        return skew23 == Skew23 && skew02 == Skew02 ? Index : -1;
    LEOPARD_FF8XOR_FOR_EACH_FOUR_TUPLE(LEO_FF8XOR_FOUR_CASE)
#undef LEO_FF8XOR_FOUR_CASE
    default:
        return -1;
    }
}

bool Apply512(
    unsigned tuple_index,
    void* a,
    void* b,
    void* c,
    void* d,
    uint64_t buffer_bytes,
    bool inverse,
    bool use_xor3)
{
    const uint64_t plane_bytes = buffer_bytes >> 3;
    // A partial ZMM plane would require replaying the generated map on a
    // narrower tail.  Fall back to the established two-way kernels instead of
    // touching bytes outside a plane or paying an extra full memory pass.
    if (tuple_index >= 64 || plane_bytes < 64 || (plane_bytes & 63) != 0)
        return false;

    const FourBufferFunction* functions = inverse
        ? (use_xor3 ? IFFTXor3Functions : IFFTXor2Functions)
        : (use_xor3 ? FFTXor3Functions : FFTXor2Functions);
    functions[tuple_index](a, b, c, d, buffer_bytes);
    return true;
}


}}} // namespace leopard::ff8xor::avx512four

#else

namespace leopard { namespace ff8xor { namespace avx512four {


bool KernelsBuilt()
{
    return false;
}


unsigned GetTupleCount()
{
    return 0;
}

bool GetTuple(unsigned, uint8_t&, uint8_t&, uint8_t&)
{
    return false;
}

int FindTupleIndex(uint8_t, uint8_t, uint8_t)
{
    return -1;
}

bool Apply512(
    unsigned,
    void*,
    void*,
    void*,
    void*,
    uint64_t,
    bool,
    bool)
{
    return false;
}


}}} // namespace leopard::ff8xor::avx512four

#endif
