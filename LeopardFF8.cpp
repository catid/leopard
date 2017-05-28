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

#include "LeopardFF8.h"

#ifdef LEO_HAS_FF8

#include <string.h>

// Define this to enable the optimized version of FWHT()
#define LEO_FF8_FWHT_OPTIMIZED

namespace leopard { namespace ff8 {


//------------------------------------------------------------------------------
// Datatypes and Constants

// Modulus for field operations
static const ffe_t kModulus = 255;

// LFSR Polynomial that generates the field elements
static const unsigned kPolynomial = 0x11D;

// Basis used for generating logarithm tables
static const ffe_t kCantorBasis[kBits] = {
    1, 214, 152, 146, 86, 200, 88, 230
};

// Using the Cantor basis here enables us to avoid a lot of extra calculations
// when applying the formal derivative in decoding.


//------------------------------------------------------------------------------
// Field Operations

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

    // Conversion to Cantor basis:

    LogLUT[0] = 0;
    for (unsigned i = 0; i < kBits; ++i)
    {
        const ffe_t basis = kCantorBasis[i];
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
// Multiplies

// We require memory to be aligned since the SIMD instructions benefit from
// or require aligned accesses to the table data.
struct {
    LEO_ALIGNED LEO_M128 Lo[256];
    LEO_ALIGNED LEO_M128 Hi[256];
} static Multiply128LUT;
#if defined(LEO_TRY_AVX2)
struct {
    LEO_ALIGNED LEO_M256 Lo[256];
    LEO_ALIGNED LEO_M256 Hi[256];
} static Multiply256LUT;
#endif // LEO_TRY_AVX2

// Returns a * Log(b)
static ffe_t FFEMultiplyLog(ffe_t a, ffe_t log_b)
{
    if (a == 0)
        return 0;
    return ExpLUT[AddMod(LogLUT[a], log_b)];
}


bool InitializeMultiplyTables()
{
    for (int log_y = 0; log_y < 256; ++log_y)
    {
        uint8_t lo[16], hi[16];
        for (unsigned char x = 0; x < 16; ++x)
        {
            lo[x] = FFEMultiplyLog(x, static_cast<uint8_t>(log_y));
            hi[x] = FFEMultiplyLog(x << 4, static_cast<uint8_t>(log_y));
        }

        const LEO_M128 table_lo = _mm_loadu_si128((LEO_M128*)lo);
        const LEO_M128 table_hi = _mm_loadu_si128((LEO_M128*)hi);

        _mm_storeu_si128(Multiply128LUT.Lo + log_y, table_lo);
        _mm_storeu_si128(Multiply128LUT.Hi + log_y, table_hi);

#if defined(LEO_TRY_AVX2)
        if (CpuHasAVX2)
        {
            _mm256_storeu_si256(Multiply256LUT.Lo + log_y,
                _mm256_broadcastsi128_si256(table_lo));
            _mm256_storeu_si256(Multiply256LUT.Hi + log_y,
                _mm256_broadcastsi128_si256(table_hi));
        }
#endif // LEO_TRY_AVX2
    }

    return true;
}


void mul_mem(
    void * LEO_RESTRICT x, const void * LEO_RESTRICT y,
    ffe_t log_m, uint64_t bytes)
{
#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
    {
        const LEO_M256 table_lo_y = _mm256_loadu_si256(Multiply256LUT.Lo + log_m);
        const LEO_M256 table_hi_y = _mm256_loadu_si256(Multiply256LUT.Hi + log_m);

        const LEO_M256 clr_mask = _mm256_set1_epi8(0x0f);

        LEO_M256 * LEO_RESTRICT x32 = reinterpret_cast<LEO_M256 *>(x);
        const LEO_M256 * LEO_RESTRICT y32 = reinterpret_cast<const LEO_M256 *>(y);

        do
        {
#define LEO_MUL_256(x_ptr, y_ptr) { \
            LEO_M256 data = _mm256_loadu_si256(y_ptr); \
            LEO_M256 lo = _mm256_and_si256(data, clr_mask); \
            LEO_M256 hi = _mm256_srli_epi64(data, 4); \
            hi = _mm256_and_si256(hi, clr_mask); \
            lo = _mm256_shuffle_epi8(table_lo_y, lo); \
            hi = _mm256_shuffle_epi8(table_hi_y, hi); \
            _mm256_storeu_si256(x_ptr, _mm256_xor_si256(lo, hi)); }

            LEO_MUL_256(x32 + 1, y32 + 1);
            LEO_MUL_256(x32, y32);
            y32 += 2, x32 += 2;

            bytes -= 64;
        } while (bytes > 0);

        return;
    }
#endif // LEO_TRY_AVX2

    const LEO_M128 table_lo_y = _mm_loadu_si128(Multiply128LUT.Lo + log_m);
    const LEO_M128 table_hi_y = _mm_loadu_si128(Multiply128LUT.Hi + log_m);

    const LEO_M128 clr_mask = _mm_set1_epi8(0x0f);

    LEO_M128 * LEO_RESTRICT x16 = reinterpret_cast<LEO_M128 *>(x);
    const LEO_M128 * LEO_RESTRICT y16 = reinterpret_cast<const LEO_M128 *>(y);

    do
    {
#define LEO_MUL_128(x_ptr, y_ptr) { \
        LEO_M128 data = _mm_loadu_si128(y_ptr); \
        LEO_M128 lo = _mm_and_si128(data, clr_mask); \
        lo = _mm_shuffle_epi8(table_lo_y, lo); \
        LEO_M128 hi = _mm_srli_epi64(data, 4); \
        hi = _mm_and_si128(hi, clr_mask); \
        hi = _mm_shuffle_epi8(table_hi_y, hi); \
        _mm_storeu_si128(x_ptr, _mm_xor_si128(lo, hi)); }

        LEO_MUL_128(x16 + 3, y16 + 3);
        LEO_MUL_128(x16 + 2, y16 + 2);
        LEO_MUL_128(x16 + 1, y16 + 1);
        LEO_MUL_128(x16, y16);
        x16 += 4, y16 += 4;

        bytes -= 64;
    } while (bytes > 0);
}


//------------------------------------------------------------------------------
// FFT Operations

void fft_butterfly(
    void * LEO_RESTRICT x, void * LEO_RESTRICT y,
    ffe_t log_m, uint64_t bytes)
{
    if (log_m == kModulus)
    {
        xor_mem(y, x, bytes);
        return;
    }

#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
    {
        const LEO_M256 table_lo_y = _mm256_loadu_si256(Multiply256LUT.Lo + log_m);
        const LEO_M256 table_hi_y = _mm256_loadu_si256(Multiply256LUT.Hi + log_m);

        const LEO_M256 clr_mask = _mm256_set1_epi8(0x0f);

        LEO_M256 * LEO_RESTRICT x32 = reinterpret_cast<LEO_M256 *>(x);
        LEO_M256 * LEO_RESTRICT y32 = reinterpret_cast<LEO_M256 *>(y);

        do
        {
#define LEO_FFTB_256(x_ptr, y_ptr) { \
            LEO_M256 y_data = _mm256_loadu_si256(y_ptr); \
            LEO_M256 lo = _mm256_and_si256(y_data, clr_mask); \
            lo = _mm256_shuffle_epi8(table_lo_y, lo); \
            LEO_M256 hi = _mm256_srli_epi64(y_data, 4); \
            hi = _mm256_and_si256(hi, clr_mask); \
            hi = _mm256_shuffle_epi8(table_hi_y, hi); \
            LEO_M256 x_data = _mm256_loadu_si256(x_ptr); \
            x_data = _mm256_xor_si256(x_data, _mm256_xor_si256(lo, hi)); \
            y_data = _mm256_xor_si256(y_data, x_data); \
            _mm256_storeu_si256(x_ptr, x_data); \
            _mm256_storeu_si256(y_ptr, y_data); }

            LEO_FFTB_256(x32 + 1, y32 + 1);
            LEO_FFTB_256(x32, y32);
            y32 += 2, x32 += 2;

            bytes -= 64;
        } while (bytes > 0);

        return;
    }
#endif // LEO_TRY_AVX2

    const LEO_M128 table_lo_y = _mm_loadu_si128(Multiply128LUT.Lo + log_m);
    const LEO_M128 table_hi_y = _mm_loadu_si128(Multiply128LUT.Hi + log_m);

    const LEO_M128 clr_mask = _mm_set1_epi8(0x0f);

    LEO_M128 * LEO_RESTRICT x16 = reinterpret_cast<LEO_M128 *>(x);
    LEO_M128 * LEO_RESTRICT y16 = reinterpret_cast<LEO_M128 *>(y);

    do
    {
#define LEO_FFTB_128(x_ptr, y_ptr) { \
        LEO_M128 y_data = _mm_loadu_si128(y_ptr); \
        LEO_M128 lo = _mm_and_si128(y_data, clr_mask); \
        lo = _mm_shuffle_epi8(table_lo_y, lo); \
        LEO_M128 hi = _mm_srli_epi64(y_data, 4); \
        hi = _mm_and_si128(hi, clr_mask); \
        hi = _mm_shuffle_epi8(table_hi_y, hi); \
        LEO_M128 x_data = _mm_loadu_si128(x_ptr); \
        x_data = _mm_xor_si128(x_data, _mm_xor_si128(lo, hi)); \
        y_data = _mm_xor_si128(y_data, x_data); \
        _mm_storeu_si128(x_ptr, x_data); \
        _mm_storeu_si128(y_ptr, y_data); }

        LEO_FFTB_128(x16 + 3, y16 + 3);
        LEO_FFTB_128(x16 + 2, y16 + 2);
        LEO_FFTB_128(x16 + 1, y16 + 1);
        LEO_FFTB_128(x16, y16);
        x16 += 4, y16 += 4;

        bytes -= 64;
    } while (bytes > 0);
}


void fft_butterfly4(
    void * LEO_RESTRICT x_0, void * LEO_RESTRICT y_0,
    void * LEO_RESTRICT x_1, void * LEO_RESTRICT y_1,
    void * LEO_RESTRICT x_2, void * LEO_RESTRICT y_2,
    void * LEO_RESTRICT x_3, void * LEO_RESTRICT y_3,
    ffe_t log_m, uint64_t bytes)
{
    if (log_m == kModulus)
    {
        xor_mem4(
            y_0, x_0,
            y_1, x_1,
            y_2, x_2,
            y_3, x_3,
            bytes);
        return;
    }

#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
    {
        const LEO_M256 table_lo_y = _mm256_loadu_si256(Multiply256LUT.Lo + log_m);
        const LEO_M256 table_hi_y = _mm256_loadu_si256(Multiply256LUT.Hi + log_m);

        const LEO_M256 clr_mask = _mm256_set1_epi8(0x0f);

        LEO_M256 * LEO_RESTRICT x32_0 = reinterpret_cast<LEO_M256 *>(x_0);
        LEO_M256 * LEO_RESTRICT y32_0 = reinterpret_cast<LEO_M256 *>(y_0);
        LEO_M256 * LEO_RESTRICT x32_1 = reinterpret_cast<LEO_M256 *>(x_1);
        LEO_M256 * LEO_RESTRICT y32_1 = reinterpret_cast<LEO_M256 *>(y_1);
        LEO_M256 * LEO_RESTRICT x32_2 = reinterpret_cast<LEO_M256 *>(x_2);
        LEO_M256 * LEO_RESTRICT y32_2 = reinterpret_cast<LEO_M256 *>(y_2);
        LEO_M256 * LEO_RESTRICT x32_3 = reinterpret_cast<LEO_M256 *>(x_3);
        LEO_M256 * LEO_RESTRICT y32_3 = reinterpret_cast<LEO_M256 *>(y_3);

        do
        {
            LEO_FFTB_256(x32_0 + 1, y32_0 + 1);
            LEO_FFTB_256(x32_0, y32_0);
            y32_0 += 2, x32_0 += 2;

            LEO_FFTB_256(x32_1 + 1, y32_1 + 1);
            LEO_FFTB_256(x32_1, y32_1);
            y32_1 += 2, x32_1 += 2;

            LEO_FFTB_256(x32_2 + 1, y32_2 + 1);
            LEO_FFTB_256(x32_2, y32_2);
            y32_2 += 2, x32_2 += 2;

            LEO_FFTB_256(x32_3 + 1, y32_3 + 1);
            LEO_FFTB_256(x32_3, y32_3);
            y32_3 += 2, x32_3 += 2;

            bytes -= 64;
        } while (bytes > 0);

        return;
    }
#endif // LEO_TRY_AVX2

    const LEO_M128 table_lo_y = _mm_loadu_si128(Multiply128LUT.Lo + log_m);
    const LEO_M128 table_hi_y = _mm_loadu_si128(Multiply128LUT.Hi + log_m);

    const LEO_M128 clr_mask = _mm_set1_epi8(0x0f);

    LEO_M128 * LEO_RESTRICT x16_0 = reinterpret_cast<LEO_M128 *>(x_0);
    LEO_M128 * LEO_RESTRICT y16_0 = reinterpret_cast<LEO_M128 *>(y_0);
    LEO_M128 * LEO_RESTRICT x16_1 = reinterpret_cast<LEO_M128 *>(x_1);
    LEO_M128 * LEO_RESTRICT y16_1 = reinterpret_cast<LEO_M128 *>(y_1);
    LEO_M128 * LEO_RESTRICT x16_2 = reinterpret_cast<LEO_M128 *>(x_2);
    LEO_M128 * LEO_RESTRICT y16_2 = reinterpret_cast<LEO_M128 *>(y_2);
    LEO_M128 * LEO_RESTRICT x16_3 = reinterpret_cast<LEO_M128 *>(x_3);
    LEO_M128 * LEO_RESTRICT y16_3 = reinterpret_cast<LEO_M128 *>(y_3);

    do
    {
        LEO_FFTB_128(x16_0 + 3, y16_0 + 3);
        LEO_FFTB_128(x16_0 + 2, y16_0 + 2);
        LEO_FFTB_128(x16_0 + 1, y16_0 + 1);
        LEO_FFTB_128(x16_0, y16_0);
        x16_0 += 4, y16_0 += 4;

        LEO_FFTB_128(x16_1 + 3, y16_1 + 3);
        LEO_FFTB_128(x16_1 + 2, y16_1 + 2);
        LEO_FFTB_128(x16_1 + 1, y16_1 + 1);
        LEO_FFTB_128(x16_1, y16_1);
        x16_1 += 4, y16_1 += 4;

        LEO_FFTB_128(x16_2 + 3, y16_2 + 3);
        LEO_FFTB_128(x16_2 + 2, y16_2 + 2);
        LEO_FFTB_128(x16_2 + 1, y16_2 + 1);
        LEO_FFTB_128(x16_2, y16_2);
        x16_2 += 4, y16_2 += 4;

        LEO_FFTB_128(x16_3 + 3, y16_3 + 3);
        LEO_FFTB_128(x16_3 + 2, y16_3 + 2);
        LEO_FFTB_128(x16_3 + 1, y16_3 + 1);
        LEO_FFTB_128(x16_3, y16_3);
        x16_3 += 4, y16_3 += 4;

        bytes -= 64;
    } while (bytes > 0);
}


//------------------------------------------------------------------------------
// IFFT Operations

void ifft_butterfly(
    void * LEO_RESTRICT x, void * LEO_RESTRICT y,
    ffe_t log_m, uint64_t bytes)
{
    if (log_m == kModulus)
    {
        xor_mem(y, x, bytes);
        return;
    }

#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
    {
        const LEO_M256 table_lo_y = _mm256_loadu_si256(Multiply256LUT.Lo + log_m);
        const LEO_M256 table_hi_y = _mm256_loadu_si256(Multiply256LUT.Hi + log_m);

        const LEO_M256 clr_mask = _mm256_set1_epi8(0x0f);

        LEO_M256 * LEO_RESTRICT x32 = reinterpret_cast<LEO_M256 *>(x);
        LEO_M256 * LEO_RESTRICT y32 = reinterpret_cast<LEO_M256 *>(y);

        do
        {
#define LEO_IFFTB_256(x_ptr, y_ptr) { \
            LEO_M256 x_data = _mm256_loadu_si256(x_ptr); \
            LEO_M256 y_data = _mm256_loadu_si256(y_ptr); \
            y_data = _mm256_xor_si256(y_data, x_data); \
            _mm256_storeu_si256(y_ptr, y_data); \
            LEO_M256 lo = _mm256_and_si256(y_data, clr_mask); \
            lo = _mm256_shuffle_epi8(table_lo_y, lo); \
            LEO_M256 hi = _mm256_srli_epi64(y_data, 4); \
            hi = _mm256_and_si256(hi, clr_mask); \
            hi = _mm256_shuffle_epi8(table_hi_y, hi); \
            x_data = _mm256_xor_si256(x_data, _mm256_xor_si256(lo, hi)); \
            _mm256_storeu_si256(x_ptr, x_data); }

            LEO_IFFTB_256(x32 + 1, y32 + 1);
            LEO_IFFTB_256(x32, y32);
            y32 += 2, x32 += 2;

            bytes -= 64;
        } while (bytes > 0);

        return;
    }
#endif // LEO_TRY_AVX2

    const LEO_M128 table_lo_y = _mm_loadu_si128(Multiply128LUT.Lo + log_m);
    const LEO_M128 table_hi_y = _mm_loadu_si128(Multiply128LUT.Hi + log_m);

    const LEO_M128 clr_mask = _mm_set1_epi8(0x0f);

    LEO_M128 * LEO_RESTRICT x16 = reinterpret_cast<LEO_M128 *>(x);
    LEO_M128 * LEO_RESTRICT y16 = reinterpret_cast<LEO_M128 *>(y);

    do
    {
#define LEO_IFFTB_128(x_ptr, y_ptr) { \
        LEO_M128 x_data = _mm_loadu_si128(x_ptr); \
        LEO_M128 y_data = _mm_loadu_si128(y_ptr); \
        y_data = _mm_xor_si128(y_data, x_data); \
        _mm_storeu_si128(y_ptr, y_data); \
        LEO_M128 lo = _mm_and_si128(y_data, clr_mask); \
        lo = _mm_shuffle_epi8(table_lo_y, lo); \
        LEO_M128 hi = _mm_srli_epi64(y_data, 4); \
        hi = _mm_and_si128(hi, clr_mask); \
        hi = _mm_shuffle_epi8(table_hi_y, hi); \
        x_data = _mm_xor_si128(x_data, _mm_xor_si128(lo, hi)); \
        _mm_storeu_si128(x_ptr, x_data); }

        LEO_IFFTB_128(x16 + 3, y16 + 3);
        LEO_IFFTB_128(x16 + 2, y16 + 2);
        LEO_IFFTB_128(x16 + 1, y16 + 1);
        LEO_IFFTB_128(x16, y16);
        x16 += 4, y16 += 4;

        bytes -= 64;
    } while (bytes > 0);
}


void ifft_butterfly4(
    void * LEO_RESTRICT x_0, void * LEO_RESTRICT y_0,
    void * LEO_RESTRICT x_1, void * LEO_RESTRICT y_1,
    void * LEO_RESTRICT x_2, void * LEO_RESTRICT y_2,
    void * LEO_RESTRICT x_3, void * LEO_RESTRICT y_3,
    ffe_t log_m, uint64_t bytes)
{
    if (log_m == kModulus)
    {
        xor_mem4(
            y_0, x_0,
            y_1, x_1,
            y_2, x_2,
            y_3, x_3,
            bytes);
        return;
    }

#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
    {
        const LEO_M256 table_lo_y = _mm256_loadu_si256(Multiply256LUT.Lo + log_m);
        const LEO_M256 table_hi_y = _mm256_loadu_si256(Multiply256LUT.Hi + log_m);

        const LEO_M256 clr_mask = _mm256_set1_epi8(0x0f);

        LEO_M256 * LEO_RESTRICT x32_0 = reinterpret_cast<LEO_M256 *>(x_0);
        LEO_M256 * LEO_RESTRICT y32_0 = reinterpret_cast<LEO_M256 *>(y_0);
        LEO_M256 * LEO_RESTRICT x32_1 = reinterpret_cast<LEO_M256 *>(x_1);
        LEO_M256 * LEO_RESTRICT y32_1 = reinterpret_cast<LEO_M256 *>(y_1);
        LEO_M256 * LEO_RESTRICT x32_2 = reinterpret_cast<LEO_M256 *>(x_2);
        LEO_M256 * LEO_RESTRICT y32_2 = reinterpret_cast<LEO_M256 *>(y_2);
        LEO_M256 * LEO_RESTRICT x32_3 = reinterpret_cast<LEO_M256 *>(x_3);
        LEO_M256 * LEO_RESTRICT y32_3 = reinterpret_cast<LEO_M256 *>(y_3);

        do
        {
            LEO_IFFTB_256(x32_0 + 1, y32_0 + 1);
            LEO_IFFTB_256(x32_0, y32_0);
            y32_0 += 2, x32_0 += 2;

            LEO_IFFTB_256(x32_1 + 1, y32_1 + 1);
            LEO_IFFTB_256(x32_1, y32_1);
            y32_1 += 2, x32_1 += 2;

            LEO_IFFTB_256(x32_2 + 1, y32_2 + 1);
            LEO_IFFTB_256(x32_2, y32_2);
            y32_2 += 2, x32_2 += 2;

            LEO_IFFTB_256(x32_3 + 1, y32_3 + 1);
            LEO_IFFTB_256(x32_3, y32_3);
            y32_3 += 2, x32_3 += 2;

            bytes -= 64;
        } while (bytes > 0);

        return;
    }
#endif // LEO_TRY_AVX2

    const LEO_M128 table_lo_y = _mm_loadu_si128(Multiply128LUT.Lo + log_m);
    const LEO_M128 table_hi_y = _mm_loadu_si128(Multiply128LUT.Hi + log_m);

    const LEO_M128 clr_mask = _mm_set1_epi8(0x0f);

    LEO_M128 * LEO_RESTRICT x16_0 = reinterpret_cast<LEO_M128 *>(x_0);
    LEO_M128 * LEO_RESTRICT y16_0 = reinterpret_cast<LEO_M128 *>(y_0);
    LEO_M128 * LEO_RESTRICT x16_1 = reinterpret_cast<LEO_M128 *>(x_1);
    LEO_M128 * LEO_RESTRICT y16_1 = reinterpret_cast<LEO_M128 *>(y_1);
    LEO_M128 * LEO_RESTRICT x16_2 = reinterpret_cast<LEO_M128 *>(x_2);
    LEO_M128 * LEO_RESTRICT y16_2 = reinterpret_cast<LEO_M128 *>(y_2);
    LEO_M128 * LEO_RESTRICT x16_3 = reinterpret_cast<LEO_M128 *>(x_3);
    LEO_M128 * LEO_RESTRICT y16_3 = reinterpret_cast<LEO_M128 *>(y_3);

    do
    {
        LEO_IFFTB_128(x16_0 + 3, y16_0 + 3);
        LEO_IFFTB_128(x16_0 + 2, y16_0 + 2);
        LEO_IFFTB_128(x16_0 + 1, y16_0 + 1);
        LEO_IFFTB_128(x16_0, y16_0);
        x16_0 += 4, y16_0 += 4;

        LEO_IFFTB_128(x16_1 + 3, y16_1 + 3);
        LEO_IFFTB_128(x16_1 + 2, y16_1 + 2);
        LEO_IFFTB_128(x16_1 + 1, y16_1 + 1);
        LEO_IFFTB_128(x16_1, y16_1);
        x16_1 += 4, y16_1 += 4;

        LEO_IFFTB_128(x16_2 + 3, y16_2 + 3);
        LEO_IFFTB_128(x16_2 + 2, y16_2 + 2);
        LEO_IFFTB_128(x16_2 + 1, y16_2 + 1);
        LEO_IFFTB_128(x16_2, y16_2);
        x16_2 += 4, y16_2 += 4;

        LEO_IFFTB_128(x16_3 + 3, y16_3 + 3);
        LEO_IFFTB_128(x16_3 + 2, y16_3 + 2);
        LEO_IFFTB_128(x16_3 + 1, y16_3 + 1);
        LEO_IFFTB_128(x16_3, y16_3);
        x16_3 += 4, y16_3 += 4;

        bytes -= 64;
    } while (bytes > 0);
}


//------------------------------------------------------------------------------
// FFT

static ffe_t FFTSkew[kModulus]; // twisted factors used in FFT
static ffe_t LogWalsh[kOrder]; // factors used in the evaluation of the error locator polynomial

void FFTInitialize()
{
    ffe_t temp[kBits - 1];

    for (unsigned i = 1; i < kBits; ++i)
        temp[i - 1] = static_cast<ffe_t>(1UL << i);

    for (unsigned m = 0; m < (kBits - 1); ++m)
    {
        const unsigned step = 1UL << (m + 1);

        FFTSkew[(1UL << m) - 1] = 0;

        for (unsigned i = m; i < (kBits - 1); ++i)
        {
            const unsigned s = (1UL << (i + 1));

            for (unsigned j = (1UL << m) - 1; j < s; j += step)
                FFTSkew[j + s] = FFTSkew[j] ^ temp[i];
        }

        temp[m] = kModulus - LogLUT[FFEMultiplyLog(temp[m], LogLUT[temp[m] ^ 1])];

        for (unsigned i = m + 1; i < (kBits - 1); ++i)
        {
            const ffe_t sum = AddMod(LogLUT[temp[i] ^ 1], temp[m]);
            temp[i] = FFEMultiplyLog(temp[i], sum);
        }
    }

    for (unsigned i = 0; i < kOrder; ++i)
        FFTSkew[i] = LogLUT[FFTSkew[i]];

    // Precalculate FWHT(Log[i]):

    for (unsigned i = 0; i < kOrder; ++i)
        LogWalsh[i] = LogLUT[i];
    LogWalsh[0] = 0;
    FWHT(LogWalsh, kBits);
}


//------------------------------------------------------------------------------
// Encode

void Encode(
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned m,
    void* const * const data,
    void** work)
{
    // work <- data

    // FIXME: Unroll first loop to eliminate this
    for (unsigned i = 0; i < m; ++i)
        memcpy(work[i], data[i], buffer_bytes);

    // work <- IFFT(data, m, m)

    for (unsigned width = 1; width < m; width <<= 1)
    {
        const unsigned range = width << 1;
        const ffe_t* skewLUT = FFTSkew + m - 1;

        for (unsigned j = width; j < m; j += range)
        {
            const ffe_t skew = skewLUT[j];

            for (unsigned i = j - width; i < j; ++i)
                ifft_butterfly(work[i], work[i + width], skew, buffer_bytes);
        }
    }

    for (unsigned i = m; i + m <= original_count; i += m)
    {
        // temp <- data + i

        void** temp = work + m;

        // FIXME: Unroll first loop to eliminate this
        for (unsigned j = 0; j < m; ++j)
            memcpy(temp[j], data[j], buffer_bytes);

        // temp <- IFFT(temp, m, m + i)

        const ffe_t* skewLUT = FFTSkew + m + i - 1;

        for (unsigned width = 1; width < m; width <<= 1)
        {
            for (unsigned j = width; j < m; j += (width << 1))
            {
                const ffe_t skew = skewLUT[j];

                for (unsigned k = j - width; k < j; ++k)
                    ifft_butterfly(temp[k], temp[k + width], skew, buffer_bytes);
            }
        }

        // work <- work XOR temp

        // FIXME: Unroll last loop to eliminate this
        for (unsigned j = 0; j < m; ++j)
            xor_mem(work[j], temp[j], buffer_bytes);
    }

    const unsigned last_count = original_count % m;
    if (last_count != 0)
    {
        const unsigned i = original_count - last_count;

        // temp <- data + i

        void** temp = work + m;

        for (unsigned j = 0; j < last_count; ++j)
            memcpy(temp[j], data[j], buffer_bytes);
        for (unsigned j = last_count; j < m; ++j)
            memset(temp[j], 0, buffer_bytes);

        // temp <- IFFT(temp, m, m + i)

        for (unsigned width = 1, shift = 1; width < m; width <<= 1, ++shift)
        {
            // Calculate stop considering that the right is all zeroes
            const unsigned stop = ((last_count + width - 1) >> shift) << shift;
            const unsigned range = width << 1;
            const ffe_t* skewLUT = FFTSkew + m + i - 1;

            for (unsigned j = width; j < stop; j += range)
            {
                const ffe_t skew = skewLUT[j];

                for (unsigned k = j - width; k < j; ++k)
                    ifft_butterfly(temp[k], temp[k + width], skew, buffer_bytes);
            }
        }

        // work <- work XOR temp

        // FIXME: Unroll last loop to eliminate this
        for (unsigned j = 0; j < m; ++j)
            xor_mem(work[j], temp[j], buffer_bytes);
    }

    // work <- FFT(work, m, 0)

    for (unsigned width = (m >> 1); width > 0; width >>= 1)
    {
        const ffe_t* skewLUT = FFTSkew + width - 1;
        const unsigned range = width << 1;

        for (unsigned j = 0; j < m; j += range)
        {
            const ffe_t skew = skewLUT[j];

            for (unsigned k = j, count = j + width; k < count; ++k)
                fft_butterfly(data[k], data[k + width], skew, buffer_bytes);
        }
    }
}


//------------------------------------------------------------------------------
// Decode

void Decode(
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned m, // NextPow2(recovery_count)
    unsigned n, // NextPow2(m + original_count) = work_count
    void* const * const original, // original_count entries
    void* const * const recovery, // recovery_count entries
    void** work) // n entries
{
    // Fill in error locations

    ffe_t ErrorLocations[kOrder];
    for (unsigned i = 0; i < recovery_count; ++i)
        ErrorLocations[i] = recovery[i] ? 0 : 1;
    for (unsigned i = recovery_count; i < m; ++i)
        ErrorLocations[i] = 1;
    for (unsigned i = 0; i < original_count; ++i)
        ErrorLocations[i + m] = original[i] ? 0 : 1;
    memset(ErrorLocations + m + original_count, 0, (n - original_count - m) * sizeof(ffe_t));

    // Evaluate error locator polynomial

    FWHT(ErrorLocations, kBits);

    for (unsigned i = 0; i < kOrder; ++i)
        ErrorLocations[i] = ((unsigned)ErrorLocations[i] * (unsigned)LogWalsh[i]) % kModulus;

    FWHT(ErrorLocations, kBits);

    // work <- recovery data

    for (unsigned i = 0; i < recovery_count; ++i)
    {
        if (recovery[i])
            mul_mem(work[i], recovery[i], ErrorLocations[i], buffer_bytes);
        else
            memset(work[i], 0, buffer_bytes);
    }
    for (unsigned i = recovery_count; i < m; ++i)
        memset(work[i], 0, buffer_bytes);

    // work <- original data

    for (unsigned i = 0; i < original_count; ++i)
    {
        if (original[i])
            mul_mem(work[m + i], original[i], ErrorLocations[m + i], buffer_bytes);
        else
            memset(work[m + i], 0, buffer_bytes);
    }
    for (unsigned i = m + original_count; i < n; ++i)
        memset(work[i], 0, buffer_bytes);

    // work <- IFFT(work, n, 0)

    for (unsigned width = 1; width < n; width <<= 1)
    {
        const unsigned range = width << 1;

        for (unsigned j = width; j < n; j += range)
        {
            const ffe_t skew = FFTSkew[j - 1];

            for (unsigned i = j - width; i < j; ++i)
                ifft_butterfly(work[i], work[i + width], skew, buffer_bytes);
        }
    }

    // work <- FormalDerivative(work, n)

    for (unsigned i = 1; i < n; ++i)
    {
        const unsigned width = ((i ^ (i - 1)) + 1) >> 1;

        // If a large number of values are being XORed:
        for (unsigned j = i - width; j < i; ++j)
            xor_mem(work[j], work[j + width], buffer_bytes);
    }

    // work <- FFT(work, n, 0) truncated to m + original_count

    const unsigned output_count = m + original_count;
    for (unsigned width = (n >> 1); width > 0; width >>= 1)
    {
        const ffe_t* skewLUT = FFTSkew + width - 1;
        const unsigned range = width << 1;

        for (unsigned j = (m < range) ? 0 : m; j < output_count; j += range)
        {
            const ffe_t skew = skewLUT[j];

            for (unsigned i = j, count = j + width; i < count; ++i)
                fft_butterfly(work[i], work[i + width], skew, buffer_bytes);
        }
    }

    // Reveal erasures

    for (unsigned i = 0; i < original_count; ++i)
        if (!original[i])
            mul_mem(work[i], work[i + m], kModulus - ErrorLocations[i], buffer_bytes);
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
    FFTInitialize();

    IsInitialized = true;
    return true;
}


}} // namespace leopard::ff8

#endif // LEO_HAS_FF8
