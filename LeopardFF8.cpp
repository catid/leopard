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

#ifdef _MSC_VER
    #pragma warning(disable: 4752) // found Intel(R) Advanced Vector Extensions; consider using /arch:AVX
#endif

namespace leopard { namespace ff8 {


//------------------------------------------------------------------------------
// Datatypes and Constants

// Basis used for generating logarithm tables
static const ffe_t kCantorBasis[kBits] = {
    1, 214, 152, 146, 86, 200, 88, 230
};

// Using the Cantor basis {2} here enables us to avoid a lot of extra calculations
// when applying the formal derivative in decoding.


//------------------------------------------------------------------------------
// Field Operations

// z = x + y (mod kModulus)
static inline ffe_t AddMod(const ffe_t a, const ffe_t b)
{
    const unsigned sum = static_cast<unsigned>(a) + b;

    // Partial reduction step, allowing for kModulus to be returned
    return static_cast<ffe_t>(sum + (sum >> kBits));
}

// z = x - y (mod kModulus)
static inline ffe_t SubMod(const ffe_t a, const ffe_t b)
{
    const unsigned dif = static_cast<unsigned>(a) - b;

    // Partial reduction step, allowing for kModulus to be returned
    return static_cast<ffe_t>(dif + (dif >> kBits));
}


//------------------------------------------------------------------------------
// Fast Walsh-Hadamard Transform (FWHT) (mod kModulus)

// {a, b} = {a + b, a - b} (Mod Q)
static LEO_FORCE_INLINE void FWHT_2(ffe_t& LEO_RESTRICT a, ffe_t& LEO_RESTRICT b)
{
    const ffe_t sum = AddMod(a, b);
    const ffe_t dif = SubMod(a, b);
    a = sum;
    b = dif;
}

#if defined(LEO_FWHT_OPT)

static LEO_FORCE_INLINE void FWHT_4(ffe_t* data, unsigned s)
{
    const unsigned s2 = s << 1;

    ffe_t t0 = data[0];
    ffe_t t1 = data[s];
    ffe_t t2 = data[s2];
    ffe_t t3 = data[s2 + s];

    FWHT_2(t0, t1);
    FWHT_2(t2, t3);
    FWHT_2(t0, t2);
    FWHT_2(t1, t3);

    data[0] = t0;
    data[s] = t1;
    data[s2] = t2;
    data[s2 + s] = t3;
}

// Decimation in time (DIT) Fast Walsh-Hadamard Transform
// Unrolls pairs of layers to perform cross-layer operations in registers
// m_truncated: Number of elements that are non-zero at the front of data
static void FWHT(ffe_t* data, const unsigned m, const unsigned m_truncated)
{
    // Decimation in time: Unroll 2 layers at a time
    unsigned dist = 1, dist4 = 4;
    for (; dist4 <= m; dist = dist4, dist4 <<= 2)
    {
        // For each set of dist*4 elements:
        for (unsigned r = 0; r < m_truncated; r += dist4)
        {
            // For each set of dist elements:
            for (unsigned i = r; i < r + dist; ++i)
                FWHT_4(data + i, dist);
        }
    }

    // If there is one layer left:
    if (dist < m)
        for (unsigned i = 0; i < dist; ++i)
            FWHT_2(data[i], data[i + dist]);
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


//------------------------------------------------------------------------------
// Logarithm Tables

static ffe_t LogLUT[kOrder];
static ffe_t ExpLUT[kOrder];


// Returns a * Log(b)
static ffe_t MultiplyLog(ffe_t a, ffe_t log_b)
{
    /*
        Note that this operation is not a normal multiplication in a finite
        field because the right operand is already a logarithm.  This is done
        because it moves K table lookups from the Decode() method into the
        initialization step that is less performance critical.  The LogWalsh[]
        table below contains precalculated logarithms so it is easier to do
        all the other multiplies in that form as well.
    */
    if (a == 0)
        return 0;
    return ExpLUT[AddMod(LogLUT[a], log_b)];
}


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

    // Conversion to Cantor basis {2}:

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

    // Generate Exp table from Log table:

    for (unsigned i = 0; i < kOrder; ++i)
        ExpLUT[LogLUT[i]] = i;

    // Note: Handles modulus wrap around with LUT
    ExpLUT[kModulus] = ExpLUT[0];
}


//------------------------------------------------------------------------------
// Multiplies

/*
    The multiplication algorithm used follows the approach outlined in {4}.
    Specifically section 6 outlines the algorithm used here for 8-bit fields.
*/

struct {
    LEO_ALIGNED LEO_M128 Value[2];
} static Multiply128LUT[kOrder];

#if defined(LEO_TRY_AVX2)
struct {
    LEO_ALIGNED LEO_M256 Value[2];
} static Multiply256LUT[kOrder];
#endif // LEO_TRY_AVX2


void InitializeMultiplyTables()
{
    if (!CpuHasSSSE3)
        return;

    // For each value we could multiply by:
    for (unsigned log_m = 0; log_m < kOrder; ++log_m)
    {
        // For each 4 bits of the finite field width in bits:
        for (unsigned i = 0, shift = 0; i < 2; ++i, shift += 4)
        {
            // Construct 16 entry LUT for PSHUFB
            uint8_t lut[16];
            for (ffe_t x = 0; x < 16; ++x)
                lut[x] = MultiplyLog(x << shift, static_cast<ffe_t>(log_m));

            // Store in 128-bit wide table
            const LEO_M128 *v_ptr = reinterpret_cast<const LEO_M128 *>(&lut[0]);
            const LEO_M128 value = _mm_loadu_si128(v_ptr);
            _mm_storeu_si128(&Multiply128LUT[log_m].Value[i], value);

            // Store in 256-bit wide table
#if defined(LEO_TRY_AVX2)
            if (CpuHasAVX2)
            {
                _mm256_storeu_si256(&Multiply256LUT[log_m].Value[i],
                    _mm256_broadcastsi128_si256(value));
            }
#endif // LEO_TRY_AVX2
        }
    }
}


void mul_mem(
    void * LEO_RESTRICT x, const void * LEO_RESTRICT y,
    ffe_t log_m, uint64_t bytes)
{
#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
    {
        const LEO_M256 table_lo_y = _mm256_loadu_si256(&Multiply256LUT[log_m].Value[0]);
        const LEO_M256 table_hi_y = _mm256_loadu_si256(&Multiply256LUT[log_m].Value[1]);

        const LEO_M256 clr_mask = _mm256_set1_epi8(0x0f);

        LEO_M256 * LEO_RESTRICT x32 = reinterpret_cast<LEO_M256 *>(x);
        const LEO_M256 * LEO_RESTRICT y32 = reinterpret_cast<const LEO_M256 *>(y);

        do
        {
#define LEO_MUL_256(x_ptr, y_ptr) { \
            LEO_M256 data = _mm256_loadu_si256(y_ptr); \
            LEO_M256 lo = _mm256_and_si256(data, clr_mask); \
            lo = _mm256_shuffle_epi8(table_lo_y, lo); \
            LEO_M256 hi = _mm256_srli_epi64(data, 4); \
            hi = _mm256_and_si256(hi, clr_mask); \
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

    if (CpuHasSSSE3)
    {
        const LEO_M128 table_lo_y = _mm_loadu_si128(&Multiply128LUT[log_m].Value[0]);
        const LEO_M128 table_hi_y = _mm_loadu_si128(&Multiply128LUT[log_m].Value[1]);

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

        return;
    }

    // Reference version:
    ffe_t * LEO_RESTRICT x1 = reinterpret_cast<ffe_t *>(x);
    const ffe_t * LEO_RESTRICT y1 = reinterpret_cast<const ffe_t *>(y);

    do
    {
        for (unsigned j = 0; j < 64; ++j)
            x1[j] = MultiplyLog(y1[j], log_m);

        x1 += 64;
        y1 += 64;
        bytes -= 64;
    } while (bytes > 0);
}


//------------------------------------------------------------------------------
// FFT

// Twisted factors used in FFT
static ffe_t FFTSkew[kModulus];

// Factors used in the evaluation of the error locator polynomial
static ffe_t LogWalsh[kOrder];


static void FFTInitialize()
{
    ffe_t temp[kBits - 1];

    // Generate FFT skew vector {1}:

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

        temp[m] = kModulus - LogLUT[MultiplyLog(temp[m], LogLUT[temp[m] ^ 1])];

        for (unsigned i = m + 1; i < (kBits - 1); ++i)
        {
            const ffe_t sum = AddMod(LogLUT[temp[i] ^ 1], temp[m]);
            temp[i] = MultiplyLog(temp[i], sum);
        }
    }

    for (unsigned i = 0; i < kOrder; ++i)
        FFTSkew[i] = LogLUT[FFTSkew[i]];

    // Precalculate FWHT(Log[i]):

    for (unsigned i = 0; i < kOrder; ++i)
        LogWalsh[i] = LogLUT[i];
    LogWalsh[0] = 0;

    FWHT(LogWalsh, kOrder, kOrder);
}

void VectorFFTButterfly(
    const uint64_t bytes,
    unsigned count,
    void** x,
    void** y,
    const ffe_t log_m)
{
    if (log_m == kModulus)
    {
        VectorXOR(bytes, count, y, x);
        return;
    }

#ifdef LEO_USE_VECTOR4_OPT
    while (count >= 4)
    {
        fft_butterfly4(
            x[0], y[0],
            x[1], y[1],
            x[2], y[2],
            x[3], y[3],
            log_m, bytes);
        x += 4, y += 4;
        count -= 4;
    }
#endif // LEO_USE_VECTOR4_OPT

    for (unsigned i = 0; i < count; ++i)
        fft_butterfly(x[i], y[i], log_m, bytes);
}

/*
    Decimation in time IFFT:

    The decimation in time IFFT algorithm allows us to unroll 2 layers at a time,
    performing calculations on local registers and faster cache memory.

    Each ^___^ below indicates a butterfly between the associated indices.

    The ifft_butterfly(x, y) operation:

        if (log_m != kModulus)
            x[] ^= exp(log(y[]) + log_m)
        y[] ^= x[]

    Layer 0:
        0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7
        ^_^ ^_^ ^_^ ^_^ ^_^ ^_^ ^_^ ^_^

    Layer 1:
        0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7
        ^___^   ^___^   ^___^   ^___^
          ^___^   ^___^   ^___^   ^___^
  
    Layer 2:
        0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7 
        ^_______^       ^_______^
          ^_______^       ^_______^
            ^_______^       ^_______^
              ^_______^       ^_______^

    Layer 3:
        0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7
        ^_______________^
          ^_______________^
            ^_______________^
              ^_______________^
                ^_______________^
                  ^_______________^
                    ^_______________^
                      ^_______________^

    DIT layer 0-1 operations, grouped 4 at a time:
        {0-1, 2-3, 0-2, 1-3},
        {4-5, 6-7, 4-6, 5-7},

    DIT layer 1-2 operations, grouped 4 at a time:
        {0-2, 4-6, 0-4, 2-6},
        {1-3, 5-7, 1-5, 3-7},

    DIT layer 2-3 operations, grouped 4 at a time:
        {0-4, 0'-4', 0-0', 4-4'},
        {1-5, 1'-5', 1-1', 5-5'},
*/

void ifft_butterfly(
    void * LEO_RESTRICT x, void * LEO_RESTRICT y,
    ffe_t log_m, uint64_t bytes)
{
#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
    {
        const LEO_M256 table_lo_y = _mm256_loadu_si256(&Multiply256LUT[log_m].Value[0]);
        const LEO_M256 table_hi_y = _mm256_loadu_si256(&Multiply256LUT[log_m].Value[1]);

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

    if (CpuHasSSSE3)
    {
        const LEO_M128 table_lo_y = _mm_loadu_si128(&Multiply128LUT[log_m].Value[0]);
        const LEO_M128 table_hi_y = _mm_loadu_si128(&Multiply128LUT[log_m].Value[1]);

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

        return;
    }

    // Reference version:
    ffe_t * LEO_RESTRICT x1 = reinterpret_cast<ffe_t *>(x);
    ffe_t * LEO_RESTRICT y1 = reinterpret_cast<ffe_t *>(y);

    do
    {
        for (unsigned j = 0; j < 64; ++j)
        {
            ffe_t x_0 = x1[j];
            ffe_t y_0 = y1[j];
            y_0 ^= x_0;
            x_0 ^= MultiplyLog(y_0, log_m);
            x1[j] = x_0;
            y1[j] = y_0;
        }

        x1 += 64;
        y1 += 64;
        bytes -= 64;
    } while (bytes > 0);
}

// 4-way butterfly
static void IFFT_DIT4(
    uint64_t bytes,
    void** work,
    unsigned dist,
    const ffe_t log_m01,
    const ffe_t log_m23,
    const ffe_t log_m02)
{
#ifdef LEO_INTERLEAVE_BUTTERFLY4_OPT

#if defined(LEO_TRY_AVX2)

    if (CpuHasAVX2)
    {
        const LEO_M256 t01_lo = _mm256_loadu_si256(&Multiply256LUT[log_m01].Value[0]);
        const LEO_M256 t01_hi = _mm256_loadu_si256(&Multiply256LUT[log_m01].Value[1]);
        const LEO_M256 t23_lo = _mm256_loadu_si256(&Multiply256LUT[log_m23].Value[0]);
        const LEO_M256 t23_hi = _mm256_loadu_si256(&Multiply256LUT[log_m23].Value[1]);
        const LEO_M256 t02_lo = _mm256_loadu_si256(&Multiply256LUT[log_m02].Value[0]);
        const LEO_M256 t02_hi = _mm256_loadu_si256(&Multiply256LUT[log_m02].Value[1]);

        const LEO_M256 clr_mask = _mm256_set1_epi8(0x0f);

        LEO_M256 * LEO_RESTRICT work0 = reinterpret_cast<LEO_M256 *>(work[0]);
        LEO_M256 * LEO_RESTRICT work1 = reinterpret_cast<LEO_M256 *>(work[dist]);
        LEO_M256 * LEO_RESTRICT work2 = reinterpret_cast<LEO_M256 *>(work[dist * 2]);
        LEO_M256 * LEO_RESTRICT work3 = reinterpret_cast<LEO_M256 *>(work[dist * 3]);

        do
        {
#define LEO_IFFTB4_256(x_reg, y_reg, table_lo, table_hi) { \
                LEO_M256 lo = _mm256_and_si256(y_reg, clr_mask); \
                lo = _mm256_shuffle_epi8(table_lo, lo); \
                LEO_M256 hi = _mm256_srli_epi64(y_reg, 4); \
                hi = _mm256_and_si256(hi, clr_mask); \
                hi = _mm256_shuffle_epi8(table_hi, hi); \
                x_reg = _mm256_xor_si256(x_reg, _mm256_xor_si256(lo, hi)); }

            LEO_M256 work0_reg = _mm256_loadu_si256(work0);
            LEO_M256 work1_reg = _mm256_loadu_si256(work1);

            // First layer:
            work1_reg = _mm256_xor_si256(work0_reg, work1_reg);
            if (log_m01 != kModulus)
            {
                LEO_IFFTB4_256(work0_reg, work1_reg, t01_lo, t01_hi);
            }

            LEO_M256 work2_reg = _mm256_loadu_si256(work2);
            LEO_M256 work3_reg = _mm256_loadu_si256(work3);

            // First layer:
            work3_reg = _mm256_xor_si256(work2_reg, work3_reg);
            if (log_m23 != kModulus)
            {
                LEO_IFFTB4_256(work2_reg, work3_reg, t23_lo, t23_hi);
            }

            // Second layer:
            work2_reg = _mm256_xor_si256(work0_reg, work2_reg);
            work3_reg = _mm256_xor_si256(work1_reg, work3_reg);
            if (log_m02 != kModulus)
            {
                LEO_IFFTB4_256(work0_reg, work2_reg, t02_lo, t02_hi);
                LEO_IFFTB4_256(work1_reg, work3_reg, t02_lo, t02_hi);
            }

            _mm256_storeu_si256(work0, work0_reg);
            _mm256_storeu_si256(work1, work1_reg);
            _mm256_storeu_si256(work2, work2_reg);
            _mm256_storeu_si256(work3, work3_reg);

            work0++, work1++, work2++, work3++;

            bytes -= 32;
        } while (bytes > 0);

        return;
    }

#endif // LEO_TRY_AVX2

    if (CpuHasSSSE3)
    {
        const LEO_M128 t01_lo = _mm_loadu_si128(&Multiply128LUT[log_m01].Value[0]);
        const LEO_M128 t01_hi = _mm_loadu_si128(&Multiply128LUT[log_m01].Value[1]);
        const LEO_M128 t23_lo = _mm_loadu_si128(&Multiply128LUT[log_m23].Value[0]);
        const LEO_M128 t23_hi = _mm_loadu_si128(&Multiply128LUT[log_m23].Value[1]);
        const LEO_M128 t02_lo = _mm_loadu_si128(&Multiply128LUT[log_m02].Value[0]);
        const LEO_M128 t02_hi = _mm_loadu_si128(&Multiply128LUT[log_m02].Value[1]);

        const LEO_M128 clr_mask = _mm_set1_epi8(0x0f);

        LEO_M128 * LEO_RESTRICT work0 = reinterpret_cast<LEO_M128 *>(work[0]);
        LEO_M128 * LEO_RESTRICT work1 = reinterpret_cast<LEO_M128 *>(work[dist]);
        LEO_M128 * LEO_RESTRICT work2 = reinterpret_cast<LEO_M128 *>(work[dist * 2]);
        LEO_M128 * LEO_RESTRICT work3 = reinterpret_cast<LEO_M128 *>(work[dist * 3]);

        do
        {
#define LEO_IFFTB4_128(x_reg, y_reg, table_lo, table_hi) { \
                LEO_M128 lo = _mm_and_si128(y_reg, clr_mask); \
                lo = _mm_shuffle_epi8(table_lo, lo); \
                LEO_M128 hi = _mm_srli_epi64(y_reg, 4); \
                hi = _mm_and_si128(hi, clr_mask); \
                hi = _mm_shuffle_epi8(table_hi, hi); \
                x_reg = _mm_xor_si128(x_reg, _mm_xor_si128(lo, hi)); }

            LEO_M128 work0_reg = _mm_loadu_si128(work0);
            LEO_M128 work1_reg = _mm_loadu_si128(work1);

            // First layer:
            work1_reg = _mm_xor_si128(work0_reg, work1_reg);
            if (log_m01 != kModulus)
            {
                LEO_IFFTB4_128(work0_reg, work1_reg, t01_lo, t01_hi);
            }

            LEO_M128 work2_reg = _mm_loadu_si128(work2);
            LEO_M128 work3_reg = _mm_loadu_si128(work3);

            // First layer:
            work3_reg = _mm_xor_si128(work2_reg, work3_reg);
            if (log_m23 != kModulus)
            {
                LEO_IFFTB4_128(work2_reg, work3_reg, t23_lo, t23_hi);
            }

            // Second layer:
            work2_reg = _mm_xor_si128(work0_reg, work2_reg);
            work3_reg = _mm_xor_si128(work1_reg, work3_reg);
            if (log_m02 != kModulus)
            {
                LEO_IFFTB4_128(work0_reg, work2_reg, t02_lo, t02_hi);
                LEO_IFFTB4_128(work1_reg, work3_reg, t02_lo, t02_hi);
            }

            _mm_storeu_si128(work0, work0_reg);
            _mm_storeu_si128(work1, work1_reg);
            _mm_storeu_si128(work2, work2_reg);
            _mm_storeu_si128(work3, work3_reg);

            work0++, work1++, work2++, work3++;

            bytes -= 16;
        } while (bytes > 0);

        return;
    }

#endif // LEO_INTERLEAVE_BUTTERFLY4_OPT

    // First layer:
    if (log_m01 == kModulus)
        xor_mem(work[dist], work[0], bytes);
    else
        ifft_butterfly(work[0], work[dist], log_m01, bytes);

    if (log_m23 == kModulus)
        xor_mem(work[dist * 3], work[dist * 2], bytes);
    else
        ifft_butterfly(work[dist * 2], work[dist * 3], log_m23, bytes);

    // Second layer:
    if (log_m02 == kModulus)
    {
        xor_mem(work[dist * 2], work[0], bytes);
        xor_mem(work[dist * 3], work[dist], bytes);
    }
    else
    {
        ifft_butterfly(work[0], work[dist * 2], log_m02, bytes);
        ifft_butterfly(work[dist], work[dist * 3], log_m02, bytes);
    }
}

void IFFT_DIT(
    const uint64_t bytes,
    void* const* data,
    const unsigned m_truncated,
    void** work,
    void** xor_result,
    const unsigned m,
    const ffe_t* skewLUT)
{
    // FIXME: Roll into first layer
    if (data)
    {
        for (unsigned i = 0; i < m_truncated; ++i)
            memcpy(work[i], data[i], bytes);
        for (unsigned i = m_truncated; i < m; ++i)
            memset(work[i], 0, bytes);
    }

    // Decimation in time: Unroll 2 layers at a time
    unsigned dist = 1, dist4 = 4;
    for (; dist4 <= m; dist = dist4, dist4 <<= 2)
    {
        // For each set of dist*4 elements:
        for (unsigned r = 0; r < m_truncated; r += dist4)
        {
            const ffe_t log_m01 = skewLUT[r + dist];
            const ffe_t log_m23 = skewLUT[r + dist * 3];
            const ffe_t log_m02 = skewLUT[r + dist * 2];

            // For each set of dist elements:
            for (unsigned i = r; i < r + dist; ++i)
            {
                IFFT_DIT4(
                    bytes,
                    work + i,
                    dist,
                    log_m01,
                    log_m23,
                    log_m02);
            }
        }

        // I tried alternating sweeps left->right and right->left to reduce cache misses.
        // It provides about 1% performance boost when done for both FFT and IFFT, so it
        // does not seem to be worth the extra complexity.

        // Clear data after the first layer
        data = nullptr;
    }

    // If there is one layer left:
    if (dist < m)
    {
        const ffe_t log_m = skewLUT[dist];

        if (log_m == kModulus)
            VectorXOR(bytes, dist, work + dist, work);
        else
        {
            for (unsigned i = 0; i < dist; ++i)
            {
                ifft_butterfly(
                    work[i],
                    work[i + dist],
                    log_m,
                    bytes);
            }
        }
    }

    // FIXME: Roll into last layer
    if (xor_result)
        for (unsigned i = 0; i < m; ++i)
            xor_mem(xor_result[i], work[i], bytes);
}

/*
    Decimation in time FFT:

    The decimation in time FFT algorithm allows us to unroll 2 layers at a time,
    performing calculations on local registers and faster cache memory.

    Each ^___^ below indicates a butterfly between the associated indices.

    The fft_butterfly(x, y) operation:

        y[] ^= x[]
        if (log_m != kModulus)
            x[] ^= exp(log(y[]) + log_m)

    Layer 0:
        0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7
        ^_______________^
          ^_______________^
            ^_______________^
              ^_______________^
                ^_______________^
                  ^_______________^
                    ^_______________^
                      ^_______________^

    Layer 1:
        0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7 
        ^_______^       ^_______^
          ^_______^       ^_______^
            ^_______^       ^_______^
              ^_______^       ^_______^
  
    Layer 2:
        0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7
        ^___^   ^___^   ^___^   ^___^
          ^___^   ^___^   ^___^   ^___^

    Layer 3:
        0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7
        ^_^ ^_^ ^_^ ^_^ ^_^ ^_^ ^_^ ^_^

    DIT layer 0-1 operations, grouped 4 at a time:
        {0-0', 4-4', 0-4, 0'-4'},
        {1-1', 5-5', 1-5, 1'-5'},

    DIT layer 1-2 operations, grouped 4 at a time:
        {0-4, 2-6, 0-2, 4-6},
        {1-5, 3-7, 1-3, 5-7},

    DIT layer 2-3 operations, grouped 4 at a time:
        {0-2, 1-3, 0-1, 2-3},
        {4-6, 5-7, 4-5, 6-7},
*/

void fft_butterfly(
    void * LEO_RESTRICT x, void * LEO_RESTRICT y,
    ffe_t log_m, uint64_t bytes)
{
#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
    {
        const LEO_M256 table_lo_y = _mm256_loadu_si256(&Multiply256LUT[log_m].Value[0]);
        const LEO_M256 table_hi_y = _mm256_loadu_si256(&Multiply256LUT[log_m].Value[1]);

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

    if (CpuHasSSSE3)
    {
        const LEO_M128 table_lo_y = _mm_loadu_si128(&Multiply128LUT[log_m].Value[0]);
        const LEO_M128 table_hi_y = _mm_loadu_si128(&Multiply128LUT[log_m].Value[1]);

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

        return;
    }

    // Reference version:
    ffe_t * LEO_RESTRICT x1 = reinterpret_cast<ffe_t *>(x);
    ffe_t * LEO_RESTRICT y1 = reinterpret_cast<ffe_t *>(y);

    do
    {
        for (unsigned j = 0; j < 64; ++j)
        {
            ffe_t x_0 = x1[j];
            ffe_t y_0 = y1[j];
            x_0 ^= MultiplyLog(y_0, log_m);
            x1[j] = x_0;
            y1[j] = y_0 ^ x_0;
        }

        x1 += 64;
        y1 += 64;
        bytes -= 64;
    } while (bytes > 0);
}

#ifdef LEO_USE_VECTOR4_OPT

void fft_butterfly4(
    void * LEO_RESTRICT x_0, void * LEO_RESTRICT y_0,
    void * LEO_RESTRICT x_1, void * LEO_RESTRICT y_1,
    void * LEO_RESTRICT x_2, void * LEO_RESTRICT y_2,
    void * LEO_RESTRICT x_3, void * LEO_RESTRICT y_3,
    ffe_t log_m, uint64_t bytes)
{
#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
    {
        const LEO_M256 table_lo_y = _mm256_loadu_si256(&Multiply256LUT[log_m].Value[0]);
        const LEO_M256 table_hi_y = _mm256_loadu_si256(&Multiply256LUT[log_m].Value[1]);

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

    if (CpuHasSSSE3)
    {
        const LEO_M128 table_lo_y = _mm_loadu_si128(&Multiply128LUT[log_m].Value[0]);
        const LEO_M128 table_hi_y = _mm_loadu_si128(&Multiply128LUT[log_m].Value[1]);

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
}

#endif // LEO_USE_VECTOR4_OPT

static void FFT_DIT4(
    uint64_t bytes,
    void** work,
    unsigned dist,
    const ffe_t log_m01,
    const ffe_t log_m23,
    const ffe_t log_m02)
{
#ifdef LEO_INTERLEAVE_BUTTERFLY4_OPT

    if (CpuHasAVX2)
    {
        const LEO_M256 t01_lo = _mm256_loadu_si256(&Multiply256LUT[log_m01].Value[0]);
        const LEO_M256 t01_hi = _mm256_loadu_si256(&Multiply256LUT[log_m01].Value[1]);
        const LEO_M256 t23_lo = _mm256_loadu_si256(&Multiply256LUT[log_m23].Value[0]);
        const LEO_M256 t23_hi = _mm256_loadu_si256(&Multiply256LUT[log_m23].Value[1]);
        const LEO_M256 t02_lo = _mm256_loadu_si256(&Multiply256LUT[log_m02].Value[0]);
        const LEO_M256 t02_hi = _mm256_loadu_si256(&Multiply256LUT[log_m02].Value[1]);

        const LEO_M256 clr_mask = _mm256_set1_epi8(0x0f);

        LEO_M256 * LEO_RESTRICT work0 = reinterpret_cast<LEO_M256 *>(work[0]);
        LEO_M256 * LEO_RESTRICT work1 = reinterpret_cast<LEO_M256 *>(work[dist]);
        LEO_M256 * LEO_RESTRICT work2 = reinterpret_cast<LEO_M256 *>(work[dist * 2]);
        LEO_M256 * LEO_RESTRICT work3 = reinterpret_cast<LEO_M256 *>(work[dist * 3]);

        do
        {
#define LEO_FFTB4_256(x_reg, y_reg, table_lo, table_hi) { \
                LEO_M256 lo = _mm256_and_si256(y_reg, clr_mask); \
                lo = _mm256_shuffle_epi8(table_lo, lo); \
                LEO_M256 hi = _mm256_srli_epi64(y_reg, 4); \
                hi = _mm256_and_si256(hi, clr_mask); \
                hi = _mm256_shuffle_epi8(table_hi, hi); \
                x_reg = _mm256_xor_si256(x_reg, _mm256_xor_si256(lo, hi)); }

            LEO_M256 work0_reg = _mm256_loadu_si256(work0);
            LEO_M256 work2_reg = _mm256_loadu_si256(work2);
            LEO_M256 work1_reg = _mm256_loadu_si256(work1);
            LEO_M256 work3_reg = _mm256_loadu_si256(work3);

            // First layer:
            if (log_m02 != kModulus)
            {
                LEO_FFTB4_256(work0_reg, work2_reg, t02_lo, t02_hi);
                LEO_FFTB4_256(work1_reg, work3_reg, t02_lo, t02_hi);
            }
            work2_reg = _mm256_xor_si256(work0_reg, work2_reg);
            work3_reg = _mm256_xor_si256(work1_reg, work3_reg);

            // Second layer:
            if (log_m01 != kModulus)
            {
                LEO_FFTB4_256(work0_reg, work1_reg, t01_lo, t01_hi);
            }
            work1_reg = _mm256_xor_si256(work0_reg, work1_reg);

            _mm256_storeu_si256(work0, work0_reg);
            _mm256_storeu_si256(work1, work1_reg);

            // First layer:
            if (log_m23 != kModulus)
            {
                LEO_FFTB4_256(work2_reg, work3_reg, t23_lo, t23_hi);
            }
            work3_reg = _mm256_xor_si256(work2_reg, work3_reg);

            _mm256_storeu_si256(work2, work2_reg);
            _mm256_storeu_si256(work3, work3_reg);

            work0++, work1++, work2++, work3++;

            bytes -= 32;
        } while (bytes > 0);

        return;
    }

    if (CpuHasSSSE3)
    {
        const LEO_M128 t01_lo = _mm_loadu_si128(&Multiply128LUT[log_m01].Value[0]);
        const LEO_M128 t01_hi = _mm_loadu_si128(&Multiply128LUT[log_m01].Value[1]);
        const LEO_M128 t23_lo = _mm_loadu_si128(&Multiply128LUT[log_m23].Value[0]);
        const LEO_M128 t23_hi = _mm_loadu_si128(&Multiply128LUT[log_m23].Value[1]);
        const LEO_M128 t02_lo = _mm_loadu_si128(&Multiply128LUT[log_m02].Value[0]);
        const LEO_M128 t02_hi = _mm_loadu_si128(&Multiply128LUT[log_m02].Value[1]);

        const LEO_M128 clr_mask = _mm_set1_epi8(0x0f);

        LEO_M128 * LEO_RESTRICT work0 = reinterpret_cast<LEO_M128 *>(work[0]);
        LEO_M128 * LEO_RESTRICT work1 = reinterpret_cast<LEO_M128 *>(work[dist]);
        LEO_M128 * LEO_RESTRICT work2 = reinterpret_cast<LEO_M128 *>(work[dist * 2]);
        LEO_M128 * LEO_RESTRICT work3 = reinterpret_cast<LEO_M128 *>(work[dist * 3]);

        do
        {
#define LEO_FFTB4_128(x_reg, y_reg, table_lo, table_hi) { \
                LEO_M128 lo = _mm_and_si128(y_reg, clr_mask); \
                lo = _mm_shuffle_epi8(table_lo, lo); \
                LEO_M128 hi = _mm_srli_epi64(y_reg, 4); \
                hi = _mm_and_si128(hi, clr_mask); \
                hi = _mm_shuffle_epi8(table_hi, hi); \
                x_reg = _mm_xor_si128(x_reg, _mm_xor_si128(lo, hi)); }

            LEO_M128 work0_reg = _mm_loadu_si128(work0);
            LEO_M128 work2_reg = _mm_loadu_si128(work2);
            LEO_M128 work1_reg = _mm_loadu_si128(work1);
            LEO_M128 work3_reg = _mm_loadu_si128(work3);

            // First layer:
            if (log_m02 != kModulus)
            {
                LEO_FFTB4_128(work0_reg, work2_reg, t02_lo, t02_hi);
                LEO_FFTB4_128(work1_reg, work3_reg, t02_lo, t02_hi);
            }
            work2_reg = _mm_xor_si128(work0_reg, work2_reg);
            work3_reg = _mm_xor_si128(work1_reg, work3_reg);

            // Second layer:
            if (log_m01 != kModulus)
            {
                LEO_FFTB4_128(work0_reg, work1_reg, t01_lo, t01_hi);
            }
            work1_reg = _mm_xor_si128(work0_reg, work1_reg);

            _mm_storeu_si128(work0, work0_reg);
            _mm_storeu_si128(work1, work1_reg);

            // First layer:
            if (log_m23 != kModulus)
            {
                LEO_FFTB4_128(work2_reg, work3_reg, t23_lo, t23_hi);
            }
            work3_reg = _mm_xor_si128(work2_reg, work3_reg);

            _mm_storeu_si128(work2, work2_reg);
            _mm_storeu_si128(work3, work3_reg);

            work0++, work1++, work2++, work3++;

            bytes -= 16;
        } while (bytes > 0);

        return;
    }

#endif // LEO_INTERLEAVE_BUTTERFLY4_OPT

    // First layer:
    if (log_m02 == kModulus)
    {
        xor_mem(work[dist * 2], work[0], bytes);
        xor_mem(work[dist * 3], work[dist], bytes);
    }
    else
    {
        fft_butterfly(work[0], work[dist * 2], log_m02, bytes);
        fft_butterfly(work[dist], work[dist * 3], log_m02, bytes);
    }

    // Second layer:
    if (log_m01 == kModulus)
        xor_mem(work[dist], work[0], bytes);
    else
        fft_butterfly(work[0], work[dist], log_m01, bytes);

    if (log_m23 == kModulus)
        xor_mem(work[dist * 3], work[dist * 2], bytes);
    else
        fft_butterfly(work[dist * 2], work[dist * 3], log_m23, bytes);
}

void FFT_DIT(
    const uint64_t bytes,
    void** work,
    const unsigned m_truncated,
    const unsigned m,
    const ffe_t* skewLUT)
{
    // Decimation in time: Unroll 2 layers at a time
    unsigned dist4 = m, dist = m >> 2;
    for (; dist != 0; dist4 = dist, dist >>= 2)
    {
        // For each set of dist*4 elements:
        for (unsigned r = 0; r < m_truncated; r += dist4)
        {
            const ffe_t log_m01 = skewLUT[r + dist];
            const ffe_t log_m23 = skewLUT[r + dist * 3];
            const ffe_t log_m02 = skewLUT[r + dist * 2];

            // For each set of dist elements:
            for (unsigned i = r; i < r + dist; ++i)
            {
                FFT_DIT4(
                    bytes,
                    work + i,
                    dist,
                    log_m01,
                    log_m23,
                    log_m02);
            }
        }
    }

    // If there is one layer left:
    if (dist4 == 2)
    {
        for (unsigned r = 0; r < m_truncated; r += 2)
        {
            const ffe_t log_m = skewLUT[r + 1];

            if (log_m == kModulus)
                xor_mem(work[r + 1], work[r], bytes);
            else
            {
                fft_butterfly(
                    work[r],
                    work[r + 1],
                    log_m,
                    bytes);
            }
        }
    }
}


//------------------------------------------------------------------------------
// Reed-Solomon Encode

void ReedSolomonEncode(
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned m,
    void* const* data,
    void** work)
{
    // work <- IFFT(data, m, m)

    const ffe_t* skewLUT = FFTSkew + m - 1;

    IFFT_DIT(
        buffer_bytes,
        data,
        original_count < m ? original_count : m,
        work,
        nullptr, // No xor output
        m,
        skewLUT);

    if (m >= original_count)
        goto skip_body;

    // For sets of m data pieces:
    for (unsigned i = m; i + m <= original_count; i += m)
    {
        data += m;
        skewLUT += m;

        // work <- work xor IFFT(data + i, m, m + i)

        IFFT_DIT(
            buffer_bytes,
            data, // data source
            m,
            work + m, // temporary workspace
            work, // xor destination
            m,
            skewLUT);
    }

    // Handle final partial set of m pieces:
    const unsigned last_count = original_count % m;
    if (last_count != 0)
    {
        const unsigned i = original_count - last_count;

        data += m;
        skewLUT += m;

        // work <- work xor IFFT(data + i, m, m + i)

        IFFT_DIT(
            buffer_bytes,
            data, // data source
            last_count,
            work + m, // temporary workspace
            work, // xor destination
            m,
            skewLUT);
    }

skip_body:

    // work <- FFT(work, m, 0)
    FFT_DIT(
        buffer_bytes,
        work,
        recovery_count,
        m,
        FFTSkew - 1);
}


//------------------------------------------------------------------------------
// ErrorBitfield

#ifdef LEO_ERROR_BITFIELD_OPT

// Used in decoding to decide which final FFT operations to perform
class ErrorBitfield
{
    static const unsigned kWords = kOrder / 64;
    uint64_t Words[7][kWords] = {};

public:
    LEO_FORCE_INLINE void Set(unsigned i)
    {
        Words[0][i / 64] |= (uint64_t)1 << (i % 64);
    }

    void Prepare();

    LEO_FORCE_INLINE bool IsNeeded(unsigned mip_level, unsigned bit) const
    {
        if (mip_level >= 8)
            return true;
        return 0 != (Words[mip_level - 1][bit / 64] & ((uint64_t)1 << (bit % 64)));
    }
};

static const uint64_t kHiMasks[5] = {
    0xAAAAAAAAAAAAAAAAULL,
    0xCCCCCCCCCCCCCCCCULL,
    0xF0F0F0F0F0F0F0F0ULL,
    0xFF00FF00FF00FF00ULL,
    0xFFFF0000FFFF0000ULL,
};

void ErrorBitfield::Prepare()
{
    // First mip level is for final layer of FFT: pairs of data
    for (unsigned i = 0; i < kWords; ++i)
    {
        uint64_t w_i = Words[0][i];
        const uint64_t hi2lo0 = w_i | ((w_i & kHiMasks[0]) >> 1);
        const uint64_t lo2hi0 = ((w_i & (kHiMasks[0] >> 1)) << 1);
        Words[0][i] = w_i = hi2lo0 | lo2hi0;

        for (unsigned j = 1, bits = 2; j < 5; ++j, bits <<= 1)
        {
            const uint64_t hi2lo_j = w_i | ((w_i & kHiMasks[j]) >> bits);
            const uint64_t lo2hi_j = ((w_i & (kHiMasks[j] >> bits)) << bits);
            Words[j][i] = w_i = hi2lo_j | lo2hi_j;
        }
    }

    for (unsigned i = 0; i < kWords; ++i)
    {
        uint64_t w = Words[4][i];
        w |= w >> 32;
        w |= w << 32;
        Words[5][i] = w;
    }

    for (unsigned i = 0; i < kWords; i += 2)
        Words[6][i] = Words[6][i + 1] = Words[5][i] | Words[5][i + 1];
}

#endif // LEO_ERROR_BITFIELD_OPT


//------------------------------------------------------------------------------
// Reed-Solomon Decode

void ReedSolomonDecode(
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

#ifdef LEO_ERROR_BITFIELD_OPT
    ErrorBitfield ErrorBits;
#endif // LEO_ERROR_BITFIELD_OPT

    ffe_t ErrorLocations[kOrder] = {};
    for (unsigned i = 0; i < recovery_count; ++i)
        if (!recovery[i])
            ErrorLocations[i] = 1;
    for (unsigned i = recovery_count; i < m; ++i)
        ErrorLocations[i] = 1;
    for (unsigned i = 0; i < original_count; ++i)
    {
        if (!original[i])
        {
            ErrorLocations[i + m] = 1;
#ifdef LEO_ERROR_BITFIELD_OPT
            ErrorBits.Set(i + m);
#endif // LEO_ERROR_BITFIELD_OPT
        }
    }

#ifdef LEO_ERROR_BITFIELD_OPT
    ErrorBits.Prepare();
#endif // LEO_ERROR_BITFIELD_OPT

    // Evaluate error locator polynomial

    FWHT(ErrorLocations, kOrder, m + original_count);

    for (unsigned i = 0; i < kOrder; ++i)
        ErrorLocations[i] = ((unsigned)ErrorLocations[i] * (unsigned)LogWalsh[i]) % kModulus;

    FWHT(ErrorLocations, kOrder, kOrder);

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

    IFFT_DIT(
        buffer_bytes,
        nullptr,
        n,
        work,
        nullptr,
        n,
        FFTSkew - 1);

    // work <- FormalDerivative(work, n)

    for (unsigned i = 1; i < n; ++i)
    {
        const unsigned width = ((i ^ (i - 1)) + 1) >> 1;

        VectorXOR(
            buffer_bytes,
            width,
            work + i - width,
            work + i);
    }

    // work <- FFT(work, n, 0) truncated to m + original_count

    unsigned mip_level = LastNonzeroBit32(n);
    const unsigned output_count = m + original_count;
    for (unsigned width = (n >> 1); width > 0; width >>= 1, --mip_level)
    {
        const ffe_t* skewLUT = FFTSkew + width - 1;
        const unsigned range = width << 1;

#ifdef LEO_SCHEDULE_OPT
        for (unsigned j = (m < range) ? 0 : m; j < output_count; j += range)
#else
        for (unsigned j = 0; j < n; j += range)
#endif
        {
#ifdef LEO_ERROR_BITFIELD_OPT
            if (!ErrorBits.IsNeeded(mip_level, j))
                continue;
#endif // LEO_ERROR_BITFIELD_OPT

            VectorFFTButterfly(
                buffer_bytes,
                width,
                work + j,
                work + j + width,
                skewLUT[j]);
        }
    }

    // Reveal erasures

    for (unsigned i = 0; i < original_count; ++i)
        if (!original[i])
            mul_mem(work[i], work[i + m], kModulus - ErrorLocations[i + m], buffer_bytes);
}


//------------------------------------------------------------------------------
// API

static bool IsInitialized = false;

bool Initialize()
{
    if (IsInitialized)
        return true;

    InitializeLogarithmTables();
    InitializeMultiplyTables();
    FFTInitialize();

    IsInitialized = true;
    return true;
}


}} // namespace leopard::ff8

#endif // LEO_HAS_FF8
