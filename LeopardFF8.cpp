/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the name of LHC-RS nor the names of its contributors may be
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

#include "LeopardFF8.h"

namespace leopard { namespace ff8 {


//------------------------------------------------------------------------------
// Datatypes and Constants

// LFSR Polynomial that generates the field elements
static const unsigned kPolynomial = 0x11D;

// Basis used for generating logarithm tables
static const ffe_t kBasis[kBits] = {
    1, 214, 152, 146, 86, 200, 88, 230 // Cantor basis
    // 1, 2, 4, 8, 16, 32, 64, 128 // Monomial basis
};


//------------------------------------------------------------------------------
// Field Operations

// Modulus for field operations
static const ffe_t kModulus = 255;

// z = x + y (mod kModulus)
static inline ffe_t AddMod(const ffe_t a, const ffe_t b)
{
    const unsigned sum = (unsigned)a + b;

    // Partial reduction step, allowing for kModulus to be returned
    return static_cast<ffe_t>(sum + (sum >> kBits));
}

// z = x - y (mod kModulus)
static inline ffe_t SubMod(const ffe_t a, const ffe_t b)
{
    const unsigned dif = (unsigned)a - b;

    // Partial reduction step, allowing for kModulus to be returned
    return static_cast<ffe_t>(dif + (dif >> kBits));
}


//------------------------------------------------------------------------------
// Logarithm Tables

static ffe_t LogLUT[kOrder];
static ffe_t ExpLUT[kOrder];


// Initialize LogLUT[], ExpLUT[]
static void InitializeLogarithmTables()
{
    // LFSR table generation:

    unsigned state = 1;
    for (unsigned i = 0; i < kModulus; ++i)
    {
        ExpLUT[state] = static_cast<ffe_t>(i);
        state <<= 1;
        if (state >= kOrder)
            state ^= kPolynomial;
    }
    ExpLUT[0] = kModulus;

    // Conversion to chosen basis:

    LogLUT[0] = 0;
    for (unsigned i = 0; i < kBits; ++i)
    {
        const ffe_t basis = kBasis[i];
        const unsigned width = static_cast<unsigned>(1UL << i);

        for (unsigned j = 0; j < width; ++j)
            LogLUT[j + width] = LogLUT[j] ^ basis;
    }

    for (unsigned i = 0; i < kOrder; ++i)
        LogLUT[i] = ExpLUT[LogLUT[i]];

    for (unsigned i = 0; i < kOrder; ++i)
        ExpLUT[LogLUT[i]] = i;

    ExpLUT[kModulus] = ExpLUT[0];
}


//------------------------------------------------------------------------------
// Fast Walsh-Hadamard Transform (FWHT) (mod kModulus)

#if defined(LEO_FF8_FWHT_OPTIMIZED)

// {a, b} = {a + b, a - b} (Mod Q)
static LEO_FORCE_INLINE void FWHT_2(ffe_t& LEO_RESTRICT a, ffe_t& LEO_RESTRICT b)
{
    const ffe_t sum = AddMod(a, b);
    const ffe_t dif = SubMod(a, b);
    a = sum;
    b = dif;
}

static LEO_FORCE_INLINE void FWHT_4(ffe_t* data)
{
    ffe_t t0 = data[0];
    ffe_t t1 = data[1];
    ffe_t t2 = data[2];
    ffe_t t3 = data[3];
    FWHT_2(t0, t1);
    FWHT_2(t2, t3);
    FWHT_2(t0, t2);
    FWHT_2(t1, t3);
    data[0] = t0;
    data[1] = t1;
    data[2] = t2;
    data[3] = t3;
}

static LEO_FORCE_INLINE void FWHT_4(ffe_t* data, unsigned s)
{
    unsigned x = 0;
    ffe_t t0 = data[x];  x += s;
    ffe_t t1 = data[x];  x += s;
    ffe_t t2 = data[x];  x += s;
    ffe_t t3 = data[x];
    FWHT_2(t0, t1);
    FWHT_2(t2, t3);
    FWHT_2(t0, t2);
    FWHT_2(t1, t3);
    unsigned y = 0;
    data[y] = t0;  y += s;
    data[y] = t1;  y += s;
    data[y] = t2;  y += s;
    data[y] = t3;
}

static inline void FWHT_8(ffe_t* data)
{
    ffe_t t0 = data[0];
    ffe_t t1 = data[1];
    ffe_t t2 = data[2];
    ffe_t t3 = data[3];
    ffe_t t4 = data[4];
    ffe_t t5 = data[5];
    ffe_t t6 = data[6];
    ffe_t t7 = data[7];
    FWHT_2(t0, t1);
    FWHT_2(t2, t3);
    FWHT_2(t4, t5);
    FWHT_2(t6, t7);
    FWHT_2(t0, t2);
    FWHT_2(t1, t3);
    FWHT_2(t4, t6);
    FWHT_2(t5, t7);
    FWHT_2(t0, t4);
    FWHT_2(t1, t5);
    FWHT_2(t2, t6);
    FWHT_2(t3, t7);
    data[0] = t0;
    data[1] = t1;
    data[2] = t2;
    data[3] = t3;
    data[4] = t4;
    data[5] = t5;
    data[6] = t6;
    data[7] = t7;
}

// Decimation in time (DIT) version
static void FWHT(ffe_t* data, const unsigned ldn)
{
    const unsigned n = (1UL << ldn);

    if (n <= 2)
    {
        if (n == 2)
            FWHT_2(data[0], data[1]);
        return;
    }

    for (unsigned ldm = ldn; ldm > 3; ldm -= 2)
    {
        unsigned m = (1UL << ldm);
        unsigned m4 = (m >> 2);
        for (unsigned r = 0; r < n; r += m)
            for (unsigned j = 0; j < m4; j++)
                FWHT_4(data + j + r, m4);
    }

    if (ldn & 1)
    {
        for (unsigned i0 = 0; i0 < n; i0 += 8)
            FWHT_8(data + i0);
    }
    else
    {
        for (unsigned i0 = 0; i0 < n; i0 += 4)
            FWHT_4(data + i0);
    }
}

#else // LEO_FF8_FWHT_OPTIMIZED

// Reference implementation
void FWHT(ffe_t* data, const unsigned bits)
{
    const unsigned size = (unsigned)(1UL << bits);
    for (unsigned width = 1; width < size; width <<= 1)
        for (unsigned i = 0; i < size; i += (width << 1))
            for (unsigned j = i; j < (width + i); ++j)
                FWHT_2(data[j], data[j + width]);
}

#endif // LEO_FF8_FWHT_OPTIMIZED

// Transform specialized for the finite field order
void FWHT(ffe_t data[kOrder])
{
    FWHT(data, kBits);
}


//------------------------------------------------------------------------------
// XOR Memory

void xor_mem(
    void * LEO_RESTRICT vx, const void * LEO_RESTRICT vy,
    unsigned bytes)
{
#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
    {
        LEO_M256 * LEO_RESTRICT x32 = reinterpret_cast<LEO_M256 *>(vx);
        const LEO_M256 * LEO_RESTRICT y32 = reinterpret_cast<const LEO_M256 *>(vy);
        do
        {
            const LEO_M256 x0 = _mm256_xor_si256(_mm256_loadu_si256(x32), _mm256_loadu_si256(y32));
            const LEO_M256 x1 = _mm256_xor_si256(_mm256_loadu_si256(x32 + 1), _mm256_loadu_si256(y32 + 1));
            const LEO_M256 x2 = _mm256_xor_si256(_mm256_loadu_si256(x32 + 2), _mm256_loadu_si256(y32 + 2));
            const LEO_M256 x3 = _mm256_xor_si256(_mm256_loadu_si256(x32 + 3), _mm256_loadu_si256(y32 + 3));
            _mm256_storeu_si256(x32, x0);
            _mm256_storeu_si256(x32 + 1, x1);
            _mm256_storeu_si256(x32 + 2, x2);
            _mm256_storeu_si256(x32 + 3, x3);
            bytes -= 128, x32 += 4, y32 += 4;
        } while (bytes >= 128);
        if (bytes > 0)
        {
            const LEO_M256 x0 = _mm256_xor_si256(_mm256_loadu_si256(x32), _mm256_loadu_si256(y32));
            const LEO_M256 x1 = _mm256_xor_si256(_mm256_loadu_si256(x32 + 1), _mm256_loadu_si256(y32 + 1));
            _mm256_storeu_si256(x32, x0);
            _mm256_storeu_si256(x32 + 1, x1);
        }
        return;
    }
#endif // LEO_TRY_AVX2
    LEO_M128 * LEO_RESTRICT x16 = reinterpret_cast<LEO_M128 *>(vx);
    const LEO_M128 * LEO_RESTRICT y16 = reinterpret_cast<const LEO_M128 *>(vy);
    do
    {
        const LEO_M128 x0 = _mm_xor_si128(_mm_loadu_si128(x16), _mm_loadu_si128(y16));
        const LEO_M128 x1 = _mm_xor_si128(_mm_loadu_si128(x16 + 1), _mm_loadu_si128(y16 + 1));
        const LEO_M128 x2 = _mm_xor_si128(_mm_loadu_si128(x16 + 2), _mm_loadu_si128(y16 + 2));
        const LEO_M128 x3 = _mm_xor_si128(_mm_loadu_si128(x16 + 3), _mm_loadu_si128(y16 + 3));
        _mm_storeu_si128(x16, x0);
        _mm_storeu_si128(x16 + 1, x1);
        _mm_storeu_si128(x16 + 2, x2);
        _mm_storeu_si128(x16 + 3, x3);
        bytes -= 64, x16 += 4, y16 += 4;
    } while (bytes > 0);
}

void xor_mem2(
    void * LEO_RESTRICT vx_0, const void * LEO_RESTRICT vy_0,
    void * LEO_RESTRICT vx_1, const void * LEO_RESTRICT vy_1,
    unsigned bytes)
{
#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
    {
        LEO_M256 * LEO_RESTRICT       x32_0 = reinterpret_cast<LEO_M256 *>      (vx_0);
        const LEO_M256 * LEO_RESTRICT y32_0 = reinterpret_cast<const LEO_M256 *>(vy_0);
        LEO_M256 * LEO_RESTRICT       x32_1 = reinterpret_cast<LEO_M256 *>      (vx_1);
        const LEO_M256 * LEO_RESTRICT y32_1 = reinterpret_cast<const LEO_M256 *>(vy_1);
        do
        {
            const LEO_M256 x0_0 = _mm256_xor_si256(_mm256_loadu_si256(x32_0),     _mm256_loadu_si256(y32_0));
            const LEO_M256 x1_0 = _mm256_xor_si256(_mm256_loadu_si256(x32_0 + 1), _mm256_loadu_si256(y32_0 + 1));
            const LEO_M256 x2_0 = _mm256_xor_si256(_mm256_loadu_si256(x32_0 + 2), _mm256_loadu_si256(y32_0 + 2));
            const LEO_M256 x3_0 = _mm256_xor_si256(_mm256_loadu_si256(x32_0 + 3), _mm256_loadu_si256(y32_0 + 3));
            const LEO_M256 x0_1 = _mm256_xor_si256(_mm256_loadu_si256(x32_1),     _mm256_loadu_si256(y32_1));
            const LEO_M256 x1_1 = _mm256_xor_si256(_mm256_loadu_si256(x32_1 + 1), _mm256_loadu_si256(y32_1 + 1));
            const LEO_M256 x2_1 = _mm256_xor_si256(_mm256_loadu_si256(x32_1 + 2), _mm256_loadu_si256(y32_1 + 2));
            const LEO_M256 x3_1 = _mm256_xor_si256(_mm256_loadu_si256(x32_1 + 3), _mm256_loadu_si256(y32_1 + 3));
            _mm256_storeu_si256(x32_0,     x0_0);
            _mm256_storeu_si256(x32_0 + 1, x1_0);
            _mm256_storeu_si256(x32_0 + 2, x2_0);
            _mm256_storeu_si256(x32_0 + 3, x3_0);
            _mm256_storeu_si256(x32_1,     x0_1);
            _mm256_storeu_si256(x32_1 + 1, x1_1);
            _mm256_storeu_si256(x32_1 + 2, x2_1);
            _mm256_storeu_si256(x32_1 + 3, x3_1);
            x32_0 += 4, y32_0 += 4;
            x32_1 += 4, y32_1 += 4;
            bytes -= 128;
        } while (bytes >= 128);
        if (bytes > 0)
        {
            const LEO_M256 x0_0 = _mm256_xor_si256(_mm256_loadu_si256(x32_0),     _mm256_loadu_si256(y32_0));
            const LEO_M256 x1_0 = _mm256_xor_si256(_mm256_loadu_si256(x32_0 + 1), _mm256_loadu_si256(y32_0 + 1));
            const LEO_M256 x0_1 = _mm256_xor_si256(_mm256_loadu_si256(x32_1),     _mm256_loadu_si256(y32_1));
            const LEO_M256 x1_1 = _mm256_xor_si256(_mm256_loadu_si256(x32_1 + 1), _mm256_loadu_si256(y32_1 + 1));
            _mm256_storeu_si256(x32_0,     x0_0);
            _mm256_storeu_si256(x32_0 + 1, x1_0);
            _mm256_storeu_si256(x32_1,     x0_1);
            _mm256_storeu_si256(x32_1 + 1, x1_1);
        }
        return;
    }
#endif // LEO_TRY_AVX2
    LEO_M128 * LEO_RESTRICT       x16_0 = reinterpret_cast<LEO_M128 *>      (vx_0);
    const LEO_M128 * LEO_RESTRICT y16_0 = reinterpret_cast<const LEO_M128 *>(vy_0);
    LEO_M128 * LEO_RESTRICT       x16_1 = reinterpret_cast<LEO_M128 *>      (vx_1);
    const LEO_M128 * LEO_RESTRICT y16_1 = reinterpret_cast<const LEO_M128 *>(vy_1);
    do
    {
        const LEO_M128 x0_0 = _mm_xor_si128(_mm_loadu_si128(x16_0),     _mm_loadu_si128(y16_0));
        const LEO_M128 x1_0 = _mm_xor_si128(_mm_loadu_si128(x16_0 + 1), _mm_loadu_si128(y16_0 + 1));
        const LEO_M128 x2_0 = _mm_xor_si128(_mm_loadu_si128(x16_0 + 2), _mm_loadu_si128(y16_0 + 2));
        const LEO_M128 x3_0 = _mm_xor_si128(_mm_loadu_si128(x16_0 + 3), _mm_loadu_si128(y16_0 + 3));
        const LEO_M128 x0_1 = _mm_xor_si128(_mm_loadu_si128(x16_1),     _mm_loadu_si128(y16_1));
        const LEO_M128 x1_1 = _mm_xor_si128(_mm_loadu_si128(x16_1 + 1), _mm_loadu_si128(y16_1 + 1));
        const LEO_M128 x2_1 = _mm_xor_si128(_mm_loadu_si128(x16_1 + 2), _mm_loadu_si128(y16_1 + 2));
        const LEO_M128 x3_1 = _mm_xor_si128(_mm_loadu_si128(x16_1 + 3), _mm_loadu_si128(y16_1 + 3));
        _mm_storeu_si128(x16_0,     x0_0);
        _mm_storeu_si128(x16_0 + 1, x1_0);
        _mm_storeu_si128(x16_0 + 2, x2_0);
        _mm_storeu_si128(x16_0 + 3, x3_0);
        _mm_storeu_si128(x16_1,     x0_1);
        _mm_storeu_si128(x16_1 + 1, x1_1);
        _mm_storeu_si128(x16_1 + 2, x2_1);
        _mm_storeu_si128(x16_1 + 3, x3_1);
        x16_0 += 4, y16_0 += 4;
        x16_1 += 4, y16_1 += 4;
        bytes -= 64;
    } while (bytes > 0);
}

void xor_mem3(
    void * LEO_RESTRICT vx_0, const void * LEO_RESTRICT vy_0,
    void * LEO_RESTRICT vx_1, const void * LEO_RESTRICT vy_1,
    void * LEO_RESTRICT vx_2, const void * LEO_RESTRICT vy_2,
    unsigned bytes)
{
#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
    {
        LEO_M256 * LEO_RESTRICT       x32_0 = reinterpret_cast<LEO_M256 *>      (vx_0);
        const LEO_M256 * LEO_RESTRICT y32_0 = reinterpret_cast<const LEO_M256 *>(vy_0);
        LEO_M256 * LEO_RESTRICT       x32_1 = reinterpret_cast<LEO_M256 *>      (vx_1);
        const LEO_M256 * LEO_RESTRICT y32_1 = reinterpret_cast<const LEO_M256 *>(vy_1);
        LEO_M256 * LEO_RESTRICT       x32_2 = reinterpret_cast<LEO_M256 *>      (vx_2);
        const LEO_M256 * LEO_RESTRICT y32_2 = reinterpret_cast<const LEO_M256 *>(vy_2);
        do
        {
            const LEO_M256 x0_0 = _mm256_xor_si256(_mm256_loadu_si256(x32_0),     _mm256_loadu_si256(y32_0));
            const LEO_M256 x1_0 = _mm256_xor_si256(_mm256_loadu_si256(x32_0 + 1), _mm256_loadu_si256(y32_0 + 1));
            const LEO_M256 x2_0 = _mm256_xor_si256(_mm256_loadu_si256(x32_0 + 2), _mm256_loadu_si256(y32_0 + 2));
            const LEO_M256 x3_0 = _mm256_xor_si256(_mm256_loadu_si256(x32_0 + 3), _mm256_loadu_si256(y32_0 + 3));
            const LEO_M256 x0_1 = _mm256_xor_si256(_mm256_loadu_si256(x32_1),     _mm256_loadu_si256(y32_1));
            const LEO_M256 x1_1 = _mm256_xor_si256(_mm256_loadu_si256(x32_1 + 1), _mm256_loadu_si256(y32_1 + 1));
            const LEO_M256 x2_1 = _mm256_xor_si256(_mm256_loadu_si256(x32_1 + 2), _mm256_loadu_si256(y32_1 + 2));
            const LEO_M256 x3_1 = _mm256_xor_si256(_mm256_loadu_si256(x32_1 + 3), _mm256_loadu_si256(y32_1 + 3));
            const LEO_M256 x0_2 = _mm256_xor_si256(_mm256_loadu_si256(x32_2),     _mm256_loadu_si256(y32_2));
            const LEO_M256 x1_2 = _mm256_xor_si256(_mm256_loadu_si256(x32_2 + 1), _mm256_loadu_si256(y32_2 + 1));
            const LEO_M256 x2_2 = _mm256_xor_si256(_mm256_loadu_si256(x32_2 + 2), _mm256_loadu_si256(y32_2 + 2));
            const LEO_M256 x3_2 = _mm256_xor_si256(_mm256_loadu_si256(x32_2 + 3), _mm256_loadu_si256(y32_2 + 3));
            _mm256_storeu_si256(x32_0,     x0_0);
            _mm256_storeu_si256(x32_0 + 1, x1_0);
            _mm256_storeu_si256(x32_0 + 2, x2_0);
            _mm256_storeu_si256(x32_0 + 3, x3_0);
            _mm256_storeu_si256(x32_1,     x0_1);
            _mm256_storeu_si256(x32_1 + 1, x1_1);
            _mm256_storeu_si256(x32_1 + 2, x2_1);
            _mm256_storeu_si256(x32_1 + 3, x3_1);
            _mm256_storeu_si256(x32_2,     x0_2);
            _mm256_storeu_si256(x32_2 + 1, x1_2);
            _mm256_storeu_si256(x32_2 + 2, x2_2);
            _mm256_storeu_si256(x32_2 + 3, x3_2);
            x32_0 += 4, y32_0 += 4;
            x32_1 += 4, y32_1 += 4;
            x32_2 += 4, y32_2 += 4;
            bytes -= 128;
        } while (bytes >= 128);
        if (bytes > 0)
        {
            const LEO_M256 x0_0 = _mm256_xor_si256(_mm256_loadu_si256(x32_0),     _mm256_loadu_si256(y32_0));
            const LEO_M256 x1_0 = _mm256_xor_si256(_mm256_loadu_si256(x32_0 + 1), _mm256_loadu_si256(y32_0 + 1));
            const LEO_M256 x0_1 = _mm256_xor_si256(_mm256_loadu_si256(x32_1),     _mm256_loadu_si256(y32_1));
            const LEO_M256 x1_1 = _mm256_xor_si256(_mm256_loadu_si256(x32_1 + 1), _mm256_loadu_si256(y32_1 + 1));
            const LEO_M256 x0_2 = _mm256_xor_si256(_mm256_loadu_si256(x32_2),     _mm256_loadu_si256(y32_2));
            const LEO_M256 x1_2 = _mm256_xor_si256(_mm256_loadu_si256(x32_2 + 1), _mm256_loadu_si256(y32_2 + 1));
            _mm256_storeu_si256(x32_0,     x0_0);
            _mm256_storeu_si256(x32_0 + 1, x1_0);
            _mm256_storeu_si256(x32_1,     x0_1);
            _mm256_storeu_si256(x32_1 + 1, x1_1);
            _mm256_storeu_si256(x32_2,     x0_2);
            _mm256_storeu_si256(x32_2 + 1, x1_2);
        }
        return;
    }
#endif // LEO_TRY_AVX2
    LEO_M128 * LEO_RESTRICT       x16_0 = reinterpret_cast<LEO_M128 *>      (vx_0);
    const LEO_M128 * LEO_RESTRICT y16_0 = reinterpret_cast<const LEO_M128 *>(vy_0);
    LEO_M128 * LEO_RESTRICT       x16_1 = reinterpret_cast<LEO_M128 *>      (vx_1);
    const LEO_M128 * LEO_RESTRICT y16_1 = reinterpret_cast<const LEO_M128 *>(vy_1);
    LEO_M128 * LEO_RESTRICT       x16_2 = reinterpret_cast<LEO_M128 *>      (vx_2);
    const LEO_M128 * LEO_RESTRICT y16_2 = reinterpret_cast<const LEO_M128 *>(vy_2);
    do
    {
        const LEO_M128 x0_0 = _mm_xor_si128(_mm_loadu_si128(x16_0),     _mm_loadu_si128(y16_0));
        const LEO_M128 x1_0 = _mm_xor_si128(_mm_loadu_si128(x16_0 + 1), _mm_loadu_si128(y16_0 + 1));
        const LEO_M128 x2_0 = _mm_xor_si128(_mm_loadu_si128(x16_0 + 2), _mm_loadu_si128(y16_0 + 2));
        const LEO_M128 x3_0 = _mm_xor_si128(_mm_loadu_si128(x16_0 + 3), _mm_loadu_si128(y16_0 + 3));
        const LEO_M128 x0_1 = _mm_xor_si128(_mm_loadu_si128(x16_1),     _mm_loadu_si128(y16_1));
        const LEO_M128 x1_1 = _mm_xor_si128(_mm_loadu_si128(x16_1 + 1), _mm_loadu_si128(y16_1 + 1));
        const LEO_M128 x2_1 = _mm_xor_si128(_mm_loadu_si128(x16_1 + 2), _mm_loadu_si128(y16_1 + 2));
        const LEO_M128 x3_1 = _mm_xor_si128(_mm_loadu_si128(x16_1 + 3), _mm_loadu_si128(y16_1 + 3));
        const LEO_M128 x0_2 = _mm_xor_si128(_mm_loadu_si128(x16_2),     _mm_loadu_si128(y16_2));
        const LEO_M128 x1_2 = _mm_xor_si128(_mm_loadu_si128(x16_2 + 1), _mm_loadu_si128(y16_2 + 1));
        const LEO_M128 x2_2 = _mm_xor_si128(_mm_loadu_si128(x16_2 + 2), _mm_loadu_si128(y16_2 + 2));
        const LEO_M128 x3_2 = _mm_xor_si128(_mm_loadu_si128(x16_2 + 3), _mm_loadu_si128(y16_2 + 3));
        _mm_storeu_si128(x16_0,     x0_0);
        _mm_storeu_si128(x16_0 + 1, x1_0);
        _mm_storeu_si128(x16_0 + 2, x2_0);
        _mm_storeu_si128(x16_0 + 3, x3_0);
        _mm_storeu_si128(x16_1,     x0_1);
        _mm_storeu_si128(x16_1 + 1, x1_1);
        _mm_storeu_si128(x16_1 + 2, x2_1);
        _mm_storeu_si128(x16_1 + 3, x3_1);
        _mm_storeu_si128(x16_2,     x0_2);
        _mm_storeu_si128(x16_2 + 1, x1_2);
        _mm_storeu_si128(x16_2 + 2, x2_2);
        _mm_storeu_si128(x16_2 + 3, x3_2);
        x16_0 += 4, y16_0 += 4;
        x16_1 += 4, y16_1 += 4;
        x16_2 += 4, y16_2 += 4;
        bytes -= 64;
    } while (bytes > 0);
}


//------------------------------------------------------------------------------
// Multiplies

// We require memory to be aligned since the SIMD instructions benefit from
// or require aligned accesses to the table data.
struct {
    LEO_ALIGNED LEO_M128 Lo[256];
    LEO_ALIGNED LEO_M128 Hi[256];
} Multiply128LUT;
#if defined(LEO_TRY_AVX2)
struct {
    LEO_ALIGNED LEO_M256 Lo[256];
    LEO_ALIGNED LEO_M256 Hi[256];
} Multiply256LUT;
#endif // LEO_TRY_AVX2

// Returns a * b
static ffe_t FFEMultiply(ffe_t a, ffe_t b)
{
    if (a == 0 || b == 0)
        return 0;
    return ExpLUT[AddMod(LogLUT[a], LogLUT[b])];
}

bool InitializeMultiplyTables()
{
    // Reuse aligned self test buffers to load table data
    uint8_t* lo = m_SelfTestBuffers.A;
    uint8_t* hi = m_SelfTestBuffers.B;

    for (int y = 0; y < 256; ++y)
    {
        for (unsigned char x = 0; x < 16; ++x)
        {
            lo[x] = FFEMultiply(x,      static_cast<uint8_t>(y));
            hi[x] = FFEMultiply(x << 4, static_cast<uint8_t>(y));
        }

        const LEO_M128 table_lo = _mm_loadu_si128((LEO_M128*)lo);
        const LEO_M128 table_hi = _mm_loadu_si128((LEO_M128*)hi);
        _mm_storeu_si128(Multiply128LUT.Lo + y, table_lo);
        _mm_storeu_si128(Multiply128LUT.Hi + y, table_hi);
#if defined(LEO_TRY_AVX2)
        if (CpuHasAVX2)
        {
            const LEO_M256 table_lo2 = _mm256_broadcastsi128_si256(table_lo);
            const LEO_M256 table_hi2 = _mm256_broadcastsi128_si256(table_hi);
            _mm256_storeu_si256(Multiply256LUT.Lo + y, table_lo2);
            _mm256_storeu_si256(Multiply256LUT.Hi + y, table_hi2);
        }
#endif // LEO_TRY_AVX2
    }

    return true;
}

// vx[] = vy[] * m
void mul_mem_set(
    void * LEO_RESTRICT vx, const void * LEO_RESTRICT vy,
    ffe_t m, unsigned bytes)
{
    if (m <= 1)
    {
        if (m == 1)
            memcpy(vx, vy, bytes);
        else
            memset(vx, 0, bytes);
        return;
    }

#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
    {
        const LEO_M256 table_lo_y = _mm256_loadu_si256(Multiply256LUT.Lo + m);
        const LEO_M256 table_hi_y = _mm256_loadu_si256(Multiply256LUT.Hi + m);

        const LEO_M256 clr_mask = _mm256_set1_epi8(0x0f);

        LEO_M256 * LEO_RESTRICT z32 = reinterpret_cast<LEO_M256 *>(vx);
        const LEO_M256 * LEO_RESTRICT x32 = reinterpret_cast<const LEO_M256 *>(vy);

        const unsigned count = bytes / 64;
        for (unsigned i = 0; i < count; ++i)
        {
            LEO_M256 x0 = _mm256_loadu_si256(x32 + i * 2);
            LEO_M256 l0 = _mm256_and_si256(x0, clr_mask);
            x0 = _mm256_srli_epi64(x0, 4);
            LEO_M256 h0 = _mm256_and_si256(x0, clr_mask);
            l0 = _mm256_shuffle_epi8(table_lo_y, l0);
            h0 = _mm256_shuffle_epi8(table_hi_y, h0);
            _mm256_storeu_si256(z32 + i * 2, _mm256_xor_si256(l0, h0));

            LEO_M256 x1 = _mm256_loadu_si256(x32 + i * 2 + 1);
            LEO_M256 l1 = _mm256_and_si256(x1, clr_mask);
            x1 = _mm256_srli_epi64(x1, 4);
            LEO_M256 h1 = _mm256_and_si256(x1, clr_mask);
            l1 = _mm256_shuffle_epi8(table_lo_y, l1);
            h1 = _mm256_shuffle_epi8(table_hi_y, h1);
            _mm256_storeu_si256(z32 + i * 2 + 1, _mm256_xor_si256(l1, h1));
        }
        return;
    }
#endif // LEO_TRY_AVX2

    const LEO_M128 table_lo_y = _mm_loadu_si128(Multiply128LUT.Lo + m);
    const LEO_M128 table_hi_y = _mm_loadu_si128(Multiply128LUT.Hi + m);

    const LEO_M128 clr_mask = _mm_set1_epi8(0x0f);

    LEO_M128 * LEO_RESTRICT       x16 = reinterpret_cast<LEO_M128 *>      (vx);
    const LEO_M128 * LEO_RESTRICT y16 = reinterpret_cast<const LEO_M128 *>(vy);

    do
    {
        LEO_M128 x3 = _mm_loadu_si128(y16 + 3);
        LEO_M128 l3 = _mm_and_si128(x3, clr_mask);
        x3 = _mm_srli_epi64(x3, 4);
        LEO_M128 h3 = _mm_and_si128(x3, clr_mask);
        l3 = _mm_shuffle_epi8(table_lo_y, l3);
        h3 = _mm_shuffle_epi8(table_hi_y, h3);

        LEO_M128 x2 = _mm_loadu_si128(y16 + 2);
        LEO_M128 l2 = _mm_and_si128(x2, clr_mask);
        x2 = _mm_srli_epi64(x2, 4);
        LEO_M128 h2 = _mm_and_si128(x2, clr_mask);
        l2 = _mm_shuffle_epi8(table_lo_y, l2);
        h2 = _mm_shuffle_epi8(table_hi_y, h2);

        LEO_M128 x1 = _mm_loadu_si128(y16 + 1);
        LEO_M128 l1 = _mm_and_si128(x1, clr_mask);
        x1 = _mm_srli_epi64(x1, 4);
        LEO_M128 h1 = _mm_and_si128(x1, clr_mask);
        l1 = _mm_shuffle_epi8(table_lo_y, l1);
        h1 = _mm_shuffle_epi8(table_hi_y, h1);

        LEO_M128 x0 = _mm_loadu_si128(y16);
        LEO_M128 l0 = _mm_and_si128(x0, clr_mask);
        x0 = _mm_srli_epi64(x0, 4);
        LEO_M128 h0 = _mm_and_si128(x0, clr_mask);
        l0 = _mm_shuffle_epi8(table_lo_y, l0);
        h0 = _mm_shuffle_epi8(table_hi_y, h0);

        _mm_storeu_si128(x16 + 3, _mm_xor_si128(l3, h3));
        _mm_storeu_si128(x16 + 2, _mm_xor_si128(l2, h2));
        _mm_storeu_si128(x16 + 1, _mm_xor_si128(l1, h1));
        _mm_storeu_si128(x16,     _mm_xor_si128(l0, h0));

        x16 += 4, y16 += 4;
        bytes -= 64;
    } while (bytes > 0);
}

// vx0[] *= m, vx1[] *= m
void mul_mem2_inplace(
    void * LEO_RESTRICT vx_0,
    void * LEO_RESTRICT vx_1,
    ffe_t m, unsigned bytes)
{
    if (m <= 1)
    {
        if (m == 0)
        {
            memset(vx_0, 0, bytes);
            memset(vx_1, 0, bytes);
        }
        return;
    }

#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
    {
        const LEO_M256 table_lo_y = _mm256_loadu_si256(Multiply256LUT.Lo + m);
        const LEO_M256 table_hi_y = _mm256_loadu_si256(Multiply256LUT.Hi + m);

        const LEO_M256 clr_mask = _mm256_set1_epi8(0x0f);

        LEO_M256 * LEO_RESTRICT x32_0 = reinterpret_cast<LEO_M256 *>(vx_0);
        LEO_M256 * LEO_RESTRICT x32_1 = reinterpret_cast<LEO_M256 *>(vx_1);

        do
        {
            LEO_M256 x0_0 = _mm256_loadu_si256(x32_0 + 1);
            LEO_M256 l0_0 = _mm256_and_si256(x0_0, clr_mask);
            x0_0 = _mm256_srli_epi64(x0_0, 4);
            LEO_M256 h0_0 = _mm256_and_si256(x0_0, clr_mask);
            l0_0 = _mm256_shuffle_epi8(table_lo_y, l0_0);
            h0_0 = _mm256_shuffle_epi8(table_hi_y, h0_0);
            l0_0 = _mm256_xor_si256(l0_0, h0_0);

            LEO_M256 x1_0 = _mm256_loadu_si256(x32_0);
            LEO_M256 l1_0 = _mm256_and_si256(x1_0, clr_mask);
            x1_0 = _mm256_srli_epi64(x1_0, 4);
            LEO_M256 h1_0 = _mm256_and_si256(x1_0, clr_mask);
            l1_0 = _mm256_shuffle_epi8(table_lo_y, l1_0);
            h1_0 = _mm256_shuffle_epi8(table_hi_y, h1_0);
            l1_0 = _mm256_xor_si256(l1_0, h1_0);

            LEO_M256 x0_1 = _mm256_loadu_si256(x32_1 + 1);
            LEO_M256 l0_1 = _mm256_and_si256(x0_1, clr_mask);
            x0_1 = _mm256_srli_epi64(x0_1, 4);
            LEO_M256 h0_1 = _mm256_and_si256(x0_1, clr_mask);
            l0_1 = _mm256_shuffle_epi8(table_lo_y, l0_1);
            h0_1 = _mm256_shuffle_epi8(table_hi_y, h0_1);
            l0_1 = _mm256_xor_si256(l0_1, h0_1);

            LEO_M256 x1_1 = _mm256_loadu_si256(x32_1);
            LEO_M256 l1_1 = _mm256_and_si256(x1_1, clr_mask);
            x1_1 = _mm256_srli_epi64(x1_1, 4);
            LEO_M256 h1_1 = _mm256_and_si256(x1_1, clr_mask);
            l1_1 = _mm256_shuffle_epi8(table_lo_y, l1_1);
            h1_1 = _mm256_shuffle_epi8(table_hi_y, h1_1);
            l1_1 = _mm256_xor_si256(l1_1, h1_1);

            _mm256_storeu_si256(x32_0 + 1, l0_0);
            _mm256_storeu_si256(x32_0, l1_0);
            _mm256_storeu_si256(x32_1 + 1, l0_1);
            _mm256_storeu_si256(x32_1, l1_1);

            x32_0 += 2;
            x32_1 += 2;
            bytes -= 64;
        } while (bytes > 0);
        return;
    }
#endif // LEO_TRY_AVX2

    const LEO_M128 table_lo_y = _mm_loadu_si128(Multiply128LUT.Lo + m);
    const LEO_M128 table_hi_y = _mm_loadu_si128(Multiply128LUT.Hi + m);

    const LEO_M128 clr_mask = _mm_set1_epi8(0x0f);

    LEO_M128 * LEO_RESTRICT x16_0 = reinterpret_cast<LEO_M128 *>(vx_0);
    LEO_M128 * LEO_RESTRICT x16_1 = reinterpret_cast<LEO_M128 *>(vx_1);

    do
    {
        LEO_M128 x3 = _mm_loadu_si128(x16_0 + 3);
        LEO_M128 l3 = _mm_and_si128(x3, clr_mask);
        x3 = _mm_srli_epi64(x3, 4);
        LEO_M128 h3 = _mm_and_si128(x3, clr_mask);
        l3 = _mm_shuffle_epi8(table_lo_y, l3);
        h3 = _mm_shuffle_epi8(table_hi_y, h3);

        LEO_M128 x2 = _mm_loadu_si128(x16_0 + 2);
        LEO_M128 l2 = _mm_and_si128(x2, clr_mask);
        x2 = _mm_srli_epi64(x2, 4);
        LEO_M128 h2 = _mm_and_si128(x2, clr_mask);
        l2 = _mm_shuffle_epi8(table_lo_y, l2);
        h2 = _mm_shuffle_epi8(table_hi_y, h2);

        LEO_M128 x1 = _mm_loadu_si128(x16_0 + 1);
        LEO_M128 l1 = _mm_and_si128(x1, clr_mask);
        x1 = _mm_srli_epi64(x1, 4);
        LEO_M128 h1 = _mm_and_si128(x1, clr_mask);
        l1 = _mm_shuffle_epi8(table_lo_y, l1);
        h1 = _mm_shuffle_epi8(table_hi_y, h1);

        LEO_M128 x0 = _mm_loadu_si128(x16_0);
        LEO_M128 l0 = _mm_and_si128(x0, clr_mask);
        x0 = _mm_srli_epi64(x0, 4);
        LEO_M128 h0 = _mm_and_si128(x0, clr_mask);
        l0 = _mm_shuffle_epi8(table_lo_y, l0);
        h0 = _mm_shuffle_epi8(table_hi_y, h0);

        _mm_storeu_si128(x16_0 + 3, _mm_xor_si128(l3, h3));
        _mm_storeu_si128(x16_0 + 2, _mm_xor_si128(l2, h2));
        _mm_storeu_si128(x16_0 + 1, _mm_xor_si128(l1, h1));
        _mm_storeu_si128(x16_0,     _mm_xor_si128(l0, h0));

        // FIXME: Add second one here

        x16_0 += 4;
        x16_1 += 4;
        bytes -= 64;
    } while (bytes > 0);
}


//------------------------------------------------------------------------------
// FFT Operations

// x[] ^= y[] * m, y[] ^= x[]
void mul_fft(
    void * LEO_RESTRICT x, void * LEO_RESTRICT y,
    ffe_t m, unsigned bytes)
{

}

// For i = {0, 1}: x_i[] ^= y_i[] * m, y_i[] ^= x_i[]
void mul_fft2(
    void * LEO_RESTRICT x_0, void * LEO_RESTRICT y_0,
    void * LEO_RESTRICT x_1, void * LEO_RESTRICT y_1,
    ffe_t m, unsigned bytes)
{

}

// For i = {0, 1, 2}: x_i[] ^= y_i[] * m, y_i[] ^= x_i[]
void mul_fft3(
    void * LEO_RESTRICT x_0, void * LEO_RESTRICT y_0,
    void * LEO_RESTRICT x_1, void * LEO_RESTRICT y_1,
    void * LEO_RESTRICT x_2, void * LEO_RESTRICT y_2,
    ffe_t m, unsigned bytes)
{

}


//------------------------------------------------------------------------------
// IFFT Operations

// y[] ^= x[], x[] ^= y[] * m
void mul_ifft(
    void * LEO_RESTRICT x, void * LEO_RESTRICT y,
    ffe_t m, unsigned bytes)
{

}

// For i = {0, 1}: y_i[] ^= x_i[], x_i[] ^= y_i[] * m
void mul_ifft2(
    void * LEO_RESTRICT x_0, void * LEO_RESTRICT y_0,
    void * LEO_RESTRICT x_1, void * LEO_RESTRICT y_1,
    ffe_t m, unsigned bytes)
{

}

// For i = {0, 1, 2}: y_i[] ^= x_i[], x_i[] ^= y_i[] * m
void mul_ifft3(
    void * LEO_RESTRICT x_0, void * LEO_RESTRICT y_0,
    void * LEO_RESTRICT x_1, void * LEO_RESTRICT y_1,
    void * LEO_RESTRICT x_2, void * LEO_RESTRICT y_2,
    ffe_t m, unsigned bytes)
{

}


//------------------------------------------------------------------------------
// API

static bool IsInitialized = false;

bool Initialize()
{
    if (IsInitialized)
        return true;

    if (!CpuHasSSSE3)
        return false;

    InitializeLogarithmTables();

    IsInitialized = true;
    return true;
}


}} // namespace leopard::ff8
