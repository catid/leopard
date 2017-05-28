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

#include "LeopardFF16.h"

#ifdef LEO_HAS_FF16

#include <string.h>

namespace leopard { namespace ff16 {


//------------------------------------------------------------------------------
// Datatypes and Constants

// Basis used for generating logarithm tables
static const ffe_t kCantorBasis[kBits] = {
    0x0001, 0xACCA, 0x3C0E, 0x163E,
    0xC582, 0xED2E, 0x914C, 0x4012,
    0x6C98, 0x10D8, 0x6A72, 0xB900,
    0xFDB8, 0xFB34, 0xFF38, 0x991E
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

#if defined(LEO_FWHT_OPT)

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

static inline void FWHT_16(ffe_t* data)
{
    ffe_t t0 = data[0];
    ffe_t t1 = data[1];
    ffe_t t2 = data[2];
    ffe_t t3 = data[3];
    ffe_t t4 = data[4];
    ffe_t t5 = data[5];
    ffe_t t6 = data[6];
    ffe_t t7 = data[7];
    ffe_t t8 = data[8];
    ffe_t t9 = data[9];
    ffe_t t10 = data[10];
    ffe_t t11 = data[11];
    ffe_t t12 = data[12];
    ffe_t t13 = data[13];
    ffe_t t14 = data[14];
    ffe_t t15 = data[15];
    FWHT_2(t0, t1);
    FWHT_2(t2, t3);
    FWHT_2(t4, t5);
    FWHT_2(t6, t7);
    FWHT_2(t8, t9);
    FWHT_2(t10, t11);
    FWHT_2(t12, t13);
    FWHT_2(t14, t15);
    FWHT_2(t0, t2);
    FWHT_2(t1, t3);
    FWHT_2(t4, t6);
    FWHT_2(t5, t7);
    FWHT_2(t8, t10);
    FWHT_2(t9, t11);
    FWHT_2(t12, t14);
    FWHT_2(t13, t15);
    FWHT_2(t0, t4);
    FWHT_2(t1, t5);
    FWHT_2(t2, t6);
    FWHT_2(t3, t7);
    FWHT_2(t8, t12);
    FWHT_2(t9, t13);
    FWHT_2(t10, t14);
    FWHT_2(t11, t15);
    FWHT_2(t0, t8);
    FWHT_2(t1, t9);
    FWHT_2(t2, t10);
    FWHT_2(t3, t11);
    FWHT_2(t4, t12);
    FWHT_2(t5, t13);
    FWHT_2(t6, t14);
    FWHT_2(t7, t15);
    data[0] = t0;
    data[1] = t1;
    data[2] = t2;
    data[3] = t3;
    data[4] = t4;
    data[5] = t5;
    data[6] = t6;
    data[7] = t7;
    data[8] = t8;
    data[9] = t9;
    data[10] = t10;
    data[11] = t11;
    data[12] = t12;
    data[13] = t13;
    data[14] = t14;
    data[15] = t15;
}

static void FWHT_SmallData(ffe_t* data, unsigned ldn)
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

// Decimation in time (DIT) version
static void FWHT(ffe_t* data, const unsigned ldn)
{
    if (ldn <= 13)
    {
        FWHT_SmallData(data, ldn);
        return;
    }

    FWHT_2(data[2], data[3]);
    FWHT_4(data + 4);
    FWHT_8(data + 8);
    FWHT_16(data + 16);
    for (unsigned ldm = 5; ldm < ldn; ++ldm)
        FWHT(data + (unsigned)(1UL << ldm), ldm);

    for (unsigned ldm = 0; ldm < ldn; ++ldm)
    {
        const unsigned mh = (1UL << ldm);
        for (unsigned t1 = 0, t2 = mh; t1 < mh; ++t1, ++t2)
            FWHT_2(data[t1], data[t2]);
    }
}

#else // LEO_FWHT_OPT

// Reference implementation
void FWHT(ffe_t* data, const unsigned bits)
{
    const unsigned size = (unsigned)(1UL << bits);
    for (unsigned width = 1; width < size; width <<= 1)
        for (unsigned i = 0; i < size; i += (width << 1))
            for (unsigned j = i; j < (width + i); ++j)
                FWHT_2(data[j], data[j + width]);
}

#endif // LEO_FWHT_OPT

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

/*
    Muladd implementation notes:

    Specialize for 1-3 rows at a time since often times we're multiplying by
    the same (skew) value repeatedly, as the ISA-L library does here:

    https://github.com/01org/isa-l/blob/master/erasure_code/gf_3vect_mad_avx.asm#L258

    Except we should be doing it for 16-bit Galois Field.
    To implement that use the ALTMAP trick from Jerasure:

    http://lab.jerasure.org/jerasure/gf-complete/blob/master/src/gf_w16.c#L1140

    Except we should also support AVX2 since that is a 40% perf boost, so put
    the high and low bytes 32 bytes instead of 16 bytes apart.

    Also I think we should go ahead and precompute the multiply tables since
    it avoids a bunch of memory lookups for each muladd, and only costs 8 MB.
*/

// We require memory to be aligned since the SIMD instructions benefit from
// or require aligned accesses to the table data.
struct {
    LEO_ALIGNED LEO_M128 LUT[65536][4];
} static Multiply128LUT;
#if defined(LEO_TRY_AVX2)
struct {
    LEO_ALIGNED LEO_M256 LUT[65536][4];
} static Multiply256LUT;
#endif // LEO_TRY_AVX2

// Returns a * b
static ffe_t FFEMultiply(ffe_t a, ffe_t b)
{
    if (a == 0 || b == 0)
        return 0;
    return ExpLUT[AddMod(LogLUT[a], LogLUT[b])];
}

// Returns a * Log(b)
static ffe_t FFEMultiplyLog(ffe_t a, ffe_t log_b)
{
    if (a == 0)
        return 0;
    return ExpLUT[AddMod(LogLUT[a], b)];
}

bool InitializeMultiplyTables()
{
    for (int y = 0; y < 256; ++y)
    {
        uint8_t lo[16], hi[16];
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
            _mm256_storeu_si256(Multiply256LUT.Lo + y,
                _mm256_broadcastsi128_si256(table_lo));
            _mm256_storeu_si256(Multiply256LUT.Hi + y,
                _mm256_broadcastsi128_si256(table_hi));
        }
#endif // LEO_TRY_AVX2
    }

    return true;
}

// vx[] = vy[] * m
void mul_mem_set(
    void * LEO_RESTRICT vx, const void * LEO_RESTRICT vy,
    ffe_t m, uint64_t bytes)
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


//------------------------------------------------------------------------------
// FFT Operations

// x[] ^= y[] * m, y[] ^= x[]
void fft_butterfly(
    void * LEO_RESTRICT x, void * LEO_RESTRICT y,
    ffe_t m, uint64_t bytes)
{

}

// For i = {0, 1, 2, 3}: x_i[] ^= y_i[] * m, y_i[] ^= x_i[]
void fft_butterfly4(
    void * LEO_RESTRICT x_0, void * LEO_RESTRICT y_0,
    void * LEO_RESTRICT x_1, void * LEO_RESTRICT y_1,
    void * LEO_RESTRICT x_2, void * LEO_RESTRICT y_2,
    void * LEO_RESTRICT x_3, void * LEO_RESTRICT y_3,
    ffe_t m, uint64_t bytes)
{

}


//------------------------------------------------------------------------------
// IFFT Operations

// y[] ^= x[], x[] ^= y[] * m
void ifft_butterfly(
    void * LEO_RESTRICT x, void * LEO_RESTRICT y,
    ffe_t m, uint64_t bytes)
{

}

// For i = {0, 1, 2, 3}: y_i[] ^= x_i[], x_i[] ^= y_i[] * m
void ifft_butterfly4(
    void * LEO_RESTRICT x_0, void * LEO_RESTRICT y_0,
    void * LEO_RESTRICT x_1, void * LEO_RESTRICT y_1,
    void * LEO_RESTRICT x_2, void * LEO_RESTRICT y_2,
    void * LEO_RESTRICT x_3, void * LEO_RESTRICT y_3,
    ffe_t m, uint64_t bytes)
{

}


//------------------------------------------------------------------------------
// FFT

static ffe_t FFTSkew[kModulus]; // twisted factors used in FFT
static ffe_t LogWalsh[kOrder]; // factors used in the evaluation of the error locator polynomial

void FFTInitialize()
{
    ffe_t temp[kBits - 1];

    for (unsigned i = 1; i < kBits; ++i)
        temp[i - 1] = (ffe_t)((unsigned)1 << i);

    for (unsigned m = 0; m < (kBits - 1); ++m)
    {
        const unsigned step = (unsigned)1 << (m + 1);

        FFTSkew[((unsigned)1 << m) - 1] = 0;

        for (unsigned i = m; i < (kBits - 1); ++i)
        {
            const unsigned s = ((unsigned)1 << (i + 1));

            for (unsigned j = ((unsigned)1 << m) - 1; j < s; j += step)
                FFTSkew[j + s] = FFTSkew[j] ^ temp[i];
        }

        // TBD: This can be cleaned up
        temp[m] = kModulus - LogLUT[FFEMultiply(temp[m], temp[m] ^ 1)];

        for (unsigned i = m + 1; i < (kBits - 1); ++i)
            temp[i] = FFEMultiplyLog(temp[i], (LogLUT[temp[i] ^ 1] + temp[m]) % kModulus);
    }

    for (unsigned i = 0; i < kOrder; ++i)
        FFTSkew[i] = LogLUT[FFTSkew[i]];

    temp[0] = kModulus - temp[0];

    for (unsigned i = 1; i < (kBits - 1); ++i)
        temp[i] = (kModulus - temp[i] + temp[i - 1]) % kModulus;

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

    // TBD: Unroll first loop to eliminate this
    for (unsigned i = 0; i < m; ++i)
        memcpy(work[i], data[i], buffer_bytes);

    // work <- IFFT(data, m, m)

    for (unsigned width = 1; width < m; width <<= 1)
    {
        for (unsigned j = width; j < m; j += (width << 1))
        {
            const ffe_t skew = FFTSkew[j + m - 1];

            if (skew != kModulus)
            {
                for (unsigned i = j - width; i < j; ++i)
                    ifft_butterfly(work[i], work[i + width], skew, buffer_bytes);
            }
            else
            {
                for (unsigned i = j - width; i < j; ++i)
                    xor_mem(work[i + width], work[i], buffer_bytes);
            }
        }
    }

    for (unsigned i = m; i + m <= original_count; i += m)
    {
        // temp <- data + i

        void** temp = work + m;

        // TBD: Unroll first loop to eliminate this
        for (unsigned j = 0; j < m; ++j)
            memcpy(temp[j], data[j], buffer_bytes);

        // temp <- IFFT(temp, m, m + i)

        for (unsigned width = 1; width < m; width <<= 1)
        {
            for (unsigned j = width; j < m; j += (width << 1))
            {
                const ffe_t skew = FFTSkew[j + m + i - 1];

                if (skew != kModulus)
                {
                    for (unsigned k = j - width; k < j; ++k)
                        ifft_butterfly(temp[k], temp[k + width], skew, buffer_bytes);
                }
                else
                {
                    for (unsigned k = j - width; k < j; ++k)
                        xor_mem(temp[k + width], temp[k], buffer_bytes);
                }
            }
        }

        // work <- work XOR temp

        // TBD: Unroll last loop to eliminate this
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

            for (unsigned j = width; j < stop; j += (width << 1))
            {
                const ffe_t skew = FFTSkew[j + m + i - 1];

                if (skew != kModulus)
                {
                    for (unsigned k = j - width; k < j; ++k)
                        ifft_butterfly(temp[k], temp[k + width], skew, buffer_bytes);
                }
                else
                {
                    for (unsigned k = j - width; k < j; ++k)
                        xor_mem(temp[k + width], temp[k], buffer_bytes);
                }
            }
        }

        // work <- work XOR temp

        // TBD: Unroll last loop to eliminate this
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

            if (skew != kModulus)
            {
                for (unsigned k = j, count = j + width; k < count; ++k)
                    fft_butterfly(data[k], data[k + width], skew, buffer_bytes);
            }
            else
            {
                for (unsigned k = j, count = j + width; k < count; ++k)
                    xor_mem(work[k + width], work[k], buffer_bytes);
            }
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
            mul_mem_set(work[i], recovery[i], ErrorLocations[i], buffer_bytes);
        else
            memset(work[i], 0, buffer_bytes);
    }
    for (unsigned i = recovery_count; i < m; ++i)
        memset(work[i], 0, buffer_bytes);

    // work <- original data

    for (unsigned i = 0; i < original_count; ++i)
    {
        if (original[i])
            mul_mem_set(work[m + i], original[i], ErrorLocations[m + i], buffer_bytes);
        else
            memset(work[m + i], 0, buffer_bytes);
    }
    for (unsigned i = m + original_count; i < n; ++i)
        memset(work[i], 0, buffer_bytes);

    // work <- IFFT(work, n, 0)

    for (unsigned width = 1; width < n; width <<= 1)
    {
        for (unsigned j = width; j < n; j += (width << 1))
        {
            const ffe_t skew = FFTSkew[j - 1];

            if (skew != kModulus)
            {
                for (unsigned i = j - width; i < j; ++i)
                    ifft_butterfly(work[i], work[i + width], skew, buffer_bytes);
            }
            else
            {
                for (unsigned i = j - width; i < j; ++i)
                    xor_mem(work[i + width], work[i], buffer_bytes);
            }
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

            if (skew != kModulus)
            {
                for (unsigned i = j; i < j + width; ++i)
                    fft_butterfly(work[i], work[i + width], skew, buffer_bytes);
            }
            else
            {
                for (unsigned i = j; i < j + width; ++i)
                    xor_mem(work[i + width], work[i], buffer_bytes);
            }
        }
    }

    // Reveal erasures

    for (unsigned i = 0; i < original_count; ++i)
        if (!original[i])
            mul_mem_set(work[i], work[i + m], kModulus - ErrorLocations[i], buffer_bytes);
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


}} // namespace leopard::ff16

#endif // LEO_HAS_FF16
