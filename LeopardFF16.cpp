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

// {a, b} = {a + b, a - b} (Mod Q)
static LEO_FORCE_INLINE void FWHT_2(ffe_t& LEO_RESTRICT a, ffe_t& LEO_RESTRICT b)
{
    const ffe_t sum = AddMod(a, b);
    const ffe_t dif = SubMod(a, b);
    a = sum;
    b = dif;
}

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
    The multiplication algorithm used follows the approach outlined in {4}.
    Specifically section 7 outlines the algorithm used here for 16-bit fields.
    The ALTMAP memory layout is used since there is no need to convert in/out.
*/

struct {
    LEO_ALIGNED LEO_M128 Lo[4];
    LEO_ALIGNED LEO_M128 Hi[4];
} static Multiply128LUT[kOrder];
#if defined(LEO_TRY_AVX2)
struct {
    LEO_ALIGNED LEO_M256 Lo[4];
    LEO_ALIGNED LEO_M256 Hi[4];
} static Multiply256LUT[kOrder];
#endif // LEO_TRY_AVX2


void InitializeMultiplyTables()
{
    // For each value we could multiply by:
    for (unsigned log_m = 0; log_m < kOrder; ++log_m)
    {
        // For each 4 bits of the finite field width in bits:
        for (unsigned i = 0, shift = 0; i < 4; ++i, shift += 4)
        {
            // Construct 16 entry LUT for PSHUFB
            uint8_t prod_lo[16], prod_hi[16];
            for (ffe_t x = 0; x < 16; ++x)
            {
                const ffe_t prod = MultiplyLog(x << shift, static_cast<ffe_t>(log_m));
                prod_lo[x] = static_cast<uint8_t>(prod);
                prod_hi[x] = static_cast<uint8_t>(prod >> 8);
            }

            const LEO_M128 value_lo = _mm_loadu_si128((LEO_M128*)prod_lo);
            const LEO_M128 value_hi = _mm_loadu_si128((LEO_M128*)prod_hi);

            // Store in 128-bit wide table
            _mm_storeu_si128(&Multiply128LUT[log_m].Lo[i], value_lo);
            _mm_storeu_si128(&Multiply128LUT[log_m].Hi[i], value_hi);

            // Store in 256-bit wide table
#if defined(LEO_TRY_AVX2)
            if (CpuHasAVX2)
            {
                _mm256_storeu_si256(&Multiply256LUT[log_m].Lo[i],
                    _mm256_broadcastsi128_si256(value_lo));
                _mm256_storeu_si256(&Multiply256LUT[log_m].Hi[i],
                    _mm256_broadcastsi128_si256(value_hi));
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
#define LEO_MUL_TABLES_256() \
        const LEO_M256 T0_lo = _mm256_loadu_si256(&Multiply256LUT[log_m].Lo[0]); \
        const LEO_M256 T1_lo = _mm256_loadu_si256(&Multiply256LUT[log_m].Lo[1]); \
        const LEO_M256 T2_lo = _mm256_loadu_si256(&Multiply256LUT[log_m].Lo[2]); \
        const LEO_M256 T3_lo = _mm256_loadu_si256(&Multiply256LUT[log_m].Lo[3]); \
        const LEO_M256 T0_hi = _mm256_loadu_si256(&Multiply256LUT[log_m].Hi[0]); \
        const LEO_M256 T1_hi = _mm256_loadu_si256(&Multiply256LUT[log_m].Hi[1]); \
        const LEO_M256 T2_hi = _mm256_loadu_si256(&Multiply256LUT[log_m].Hi[2]); \
        const LEO_M256 T3_hi = _mm256_loadu_si256(&Multiply256LUT[log_m].Hi[3]);

        LEO_MUL_TABLES_256();

        const LEO_M256 clr_mask = _mm256_set1_epi8(0x0f);

        LEO_M256 * LEO_RESTRICT x32 = reinterpret_cast<LEO_M256 *>(x);
        const LEO_M256 * LEO_RESTRICT y32 = reinterpret_cast<const LEO_M256 *>(y);

        do
        {
#define LEO_MUL_256(x_ptr, y_ptr) { \
            const LEO_M256 A_lo = _mm256_loadu_si256(y_ptr); \
            const LEO_M256 A_hi = _mm256_loadu_si256(y_ptr + 1); \
            LEO_M256 data_0 = _mm256_and_si256(A_lo, clr_mask); \
            LEO_M256 data_1 = _mm256_srli_epi64(A_lo, 4); \
            data_1 = _mm256_and_si256(data_1, clr_mask); \
            LEO_M256 data_2 = _mm256_and_si256(A_hi, clr_mask); \
            LEO_M256 data_3 = _mm256_srli_epi64(A_hi, 4); \
            data_3 = _mm256_and_si256(data_3, clr_mask); \
            LEO_M256 output_lo = _mm256_shuffle_epi8(T0_lo, data_0); \
            output_lo = _mm256_xor_si256(output_lo, _mm256_shuffle_epi8(T1_lo, data_1)); \
            output_lo = _mm256_xor_si256(output_lo, _mm256_shuffle_epi8(T2_lo, data_2)); \
            output_lo = _mm256_xor_si256(output_lo, _mm256_shuffle_epi8(T3_lo, data_3)); \
            LEO_M256 output_hi = _mm256_shuffle_epi8(T0_hi, data_0); \
            output_hi = _mm256_xor_si256(output_hi, _mm256_shuffle_epi8(T1_hi, data_1)); \
            output_hi = _mm256_xor_si256(output_hi, _mm256_shuffle_epi8(T2_hi, data_2)); \
            output_hi = _mm256_xor_si256(output_hi, _mm256_shuffle_epi8(T3_hi, data_3)); \
            _mm256_storeu_si256(x_ptr, output_lo); \
            _mm256_storeu_si256(x_ptr + 1, output_hi); }

            LEO_MUL_256(x32, y32);
            y32 += 2, x32 += 2;

            bytes -= 64;
        } while (bytes > 0);

        return;
    }
#endif // LEO_TRY_AVX2

#define LEO_MUL_TABLES_128() \
    const LEO_M128 T0_lo = _mm_loadu_si128(&Multiply128LUT[log_m].Lo[0]); \
    const LEO_M128 T1_lo = _mm_loadu_si128(&Multiply128LUT[log_m].Lo[1]); \
    const LEO_M128 T2_lo = _mm_loadu_si128(&Multiply128LUT[log_m].Lo[2]); \
    const LEO_M128 T3_lo = _mm_loadu_si128(&Multiply128LUT[log_m].Lo[3]); \
    const LEO_M128 T0_hi = _mm_loadu_si128(&Multiply128LUT[log_m].Hi[0]); \
    const LEO_M128 T1_hi = _mm_loadu_si128(&Multiply128LUT[log_m].Hi[1]); \
    const LEO_M128 T2_hi = _mm_loadu_si128(&Multiply128LUT[log_m].Hi[2]); \
    const LEO_M128 T3_hi = _mm_loadu_si128(&Multiply128LUT[log_m].Hi[3]);

    LEO_MUL_TABLES_128();

    const LEO_M128 clr_mask = _mm_set1_epi8(0x0f);

    LEO_M128 * LEO_RESTRICT x16 = reinterpret_cast<LEO_M128 *>(x);
    const LEO_M128 * LEO_RESTRICT y16 = reinterpret_cast<const LEO_M128 *>(y);

    do
    {
#define LEO_MUL_128(x_ptr, y_ptr) { \
            const LEO_M128 A_lo = _mm_loadu_si128(y_ptr); \
            const LEO_M128 A_hi = _mm_loadu_si128(y_ptr + 2); \
            LEO_M128 data_0 = _mm_and_si128(A_lo, clr_mask); \
            LEO_M128 data_1 = _mm_srli_epi64(A_lo, 4); \
            data_1 = _mm_and_si128(data_1, clr_mask); \
            LEO_M128 data_2 = _mm_and_si128(A_hi, clr_mask); \
            LEO_M128 data_3 = _mm_srli_epi64(A_hi, 4); \
            data_3 = _mm_and_si128(data_3, clr_mask); \
            LEO_M128 output_lo = _mm_shuffle_epi8(T0_lo, data_0); \
            output_lo = _mm_xor_si128(output_lo, _mm_shuffle_epi8(T1_lo, data_1)); \
            output_lo = _mm_xor_si128(output_lo, _mm_shuffle_epi8(T2_lo, data_2)); \
            output_lo = _mm_xor_si128(output_lo, _mm_shuffle_epi8(T3_lo, data_3)); \
            LEO_M128 output_hi = _mm_shuffle_epi8(T0_hi, data_0); \
            output_hi = _mm_xor_si128(output_hi, _mm_shuffle_epi8(T1_hi, data_1)); \
            output_hi = _mm_xor_si128(output_hi, _mm_shuffle_epi8(T2_hi, data_2)); \
            output_hi = _mm_xor_si128(output_hi, _mm_shuffle_epi8(T3_hi, data_3)); \
            _mm_storeu_si128(x_ptr, output_lo); \
            _mm_storeu_si128(x_ptr + 2, output_hi); }

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
#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
    {
        LEO_MUL_TABLES_256();

        const LEO_M256 clr_mask = _mm256_set1_epi8(0x0f);

        LEO_M256 * LEO_RESTRICT x32 = reinterpret_cast<LEO_M256 *>(x);
        LEO_M256 * LEO_RESTRICT y32 = reinterpret_cast<LEO_M256 *>(y);

        do
        {
#define LEO_FFTB_256(x_ptr, y_ptr) { \
            LEO_M256 y_lo = _mm256_loadu_si256(y_ptr); \
            LEO_M256 data_0 = _mm256_and_si256(y_lo, clr_mask); \
            LEO_M256 data_1 = _mm256_srli_epi64(y_lo, 4); \
            data_1 = _mm256_and_si256(data_1, clr_mask); \
            LEO_M256 prod_lo = _mm256_shuffle_epi8(T0_lo, data_0); \
            prod_lo = _mm256_xor_si256(prod_lo, _mm256_shuffle_epi8(T1_lo, data_1)); \
            LEO_M256 prod_hi = _mm256_shuffle_epi8(T0_hi, data_0); \
            prod_hi = _mm256_xor_si256(prod_hi, _mm256_shuffle_epi8(T1_hi, data_1)); \
            LEO_M256 y_hi = _mm256_loadu_si256(y_ptr + 1); \
            data_0 = _mm256_and_si256(y_hi, clr_mask); \
            data_1 = _mm256_srli_epi64(y_hi, 4); \
            data_1 = _mm256_and_si256(data_1, clr_mask); \
            prod_lo = _mm256_xor_si256(prod_lo, _mm256_shuffle_epi8(T2_lo, data_0)); \
            prod_lo = _mm256_xor_si256(prod_lo, _mm256_shuffle_epi8(T3_lo, data_1)); \
            prod_hi = _mm256_xor_si256(prod_hi, _mm256_shuffle_epi8(T2_hi, data_0)); \
            prod_hi = _mm256_xor_si256(prod_hi, _mm256_shuffle_epi8(T3_hi, data_1)); \
            LEO_M256 x_lo = _mm256_loadu_si256(x_ptr); \
            LEO_M256 x_hi = _mm256_loadu_si256(x_ptr + 1); \
            x_lo = _mm256_xor_si256(prod_lo, x_lo); \
            _mm256_storeu_si256(x_ptr, x_lo); \
            x_hi = _mm256_xor_si256(prod_hi, x_hi); \
            _mm256_storeu_si256(x_ptr + 1, x_hi); \
            y_lo = _mm256_xor_si256(y_lo, x_lo); \
            _mm256_storeu_si256(y_ptr, y_lo); \
            y_hi = _mm256_xor_si256(y_hi, x_hi); \
            _mm256_storeu_si256(y_ptr + 1, y_hi); }

            LEO_FFTB_256(x32, y32);
            y32 += 2, x32 += 2;

            bytes -= 64;
        } while (bytes > 0);

        return;
    }
#endif // LEO_TRY_AVX2

    LEO_MUL_TABLES_128();

    const LEO_M128 clr_mask = _mm_set1_epi8(0x0f);

    LEO_M128 * LEO_RESTRICT x16 = reinterpret_cast<LEO_M128 *>(x);
    LEO_M128 * LEO_RESTRICT y16 = reinterpret_cast<LEO_M128 *>(y);

    do
    {
#define LEO_FFTB_128(x_ptr, y_ptr) { \
            LEO_M128 y_lo = _mm_loadu_si128(y_ptr); \
            LEO_M128 data_0 = _mm_and_si128(y_lo, clr_mask); \
            LEO_M128 data_1 = _mm_srli_epi64(y_lo, 4); \
            data_1 = _mm_and_si128(data_1, clr_mask); \
            LEO_M128 prod_lo = _mm_shuffle_epi8(T0_lo, data_0); \
            prod_lo = _mm_xor_si128(prod_lo, _mm_shuffle_epi8(T1_lo, data_1)); \
            LEO_M128 prod_hi = _mm_shuffle_epi8(T0_hi, data_0); \
            prod_hi = _mm_xor_si128(prod_hi, _mm_shuffle_epi8(T1_hi, data_1)); \
            LEO_M128 y_hi = _mm_loadu_si128(y_ptr + 2); \
            data_0 = _mm_and_si128(y_hi, clr_mask); \
            data_1 = _mm_srli_epi64(y_hi, 4); \
            data_1 = _mm_and_si128(data_1, clr_mask); \
            prod_lo = _mm_xor_si128(prod_lo, _mm_shuffle_epi8(T2_lo, data_0)); \
            prod_lo = _mm_xor_si128(prod_lo, _mm_shuffle_epi8(T3_lo, data_1)); \
            prod_hi = _mm_xor_si128(prod_hi, _mm_shuffle_epi8(T2_hi, data_0)); \
            prod_hi = _mm_xor_si128(prod_hi, _mm_shuffle_epi8(T3_hi, data_1)); \
            LEO_M128 x_lo = _mm_loadu_si128(x_ptr); \
            LEO_M128 x_hi = _mm_loadu_si128(x_ptr + 2); \
            x_lo = _mm_xor_si128(prod_lo, x_lo); \
            _mm_storeu_si128(x_ptr, x_lo); \
            x_hi = _mm_xor_si128(prod_hi, x_hi); \
            _mm_storeu_si128(x_ptr + 2, x_hi); \
            y_lo = _mm_xor_si128(y_lo, x_lo); \
            _mm_storeu_si128(y_ptr, y_lo); \
            y_hi = _mm_xor_si128(y_hi, x_hi); \
            _mm_storeu_si128(y_ptr + 2, y_hi); }

        LEO_FFTB_128(x16 + 1, y16 + 1);
        LEO_FFTB_128(x16, y16);
        x16 += 4, y16 += 4;

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
        LEO_MUL_TABLES_256();

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
            LEO_FFTB_256(x32_0, y32_0);
            y32_0 += 2, x32_0 += 2;

            LEO_FFTB_256(x32_1, y32_1);
            y32_1 += 2, x32_1 += 2;

            LEO_FFTB_256(x32_2, y32_2);
            y32_2 += 2, x32_2 += 2;

            LEO_FFTB_256(x32_3, y32_3);
            y32_3 += 2, x32_3 += 2;

            bytes -= 64;
        } while (bytes > 0);

        return;
    }
#endif // LEO_TRY_AVX2

    LEO_MUL_TABLES_128();

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
        LEO_FFTB_128(x16_0 + 1, y16_0 + 1);
        LEO_FFTB_128(x16_0, y16_0);
        x16_0 += 4, y16_0 += 4;

        LEO_FFTB_128(x16_1 + 1, y16_1 + 1);
        LEO_FFTB_128(x16_1, y16_1);
        x16_1 += 4, y16_1 += 4;

        LEO_FFTB_128(x16_2 + 1, y16_2 + 1);
        LEO_FFTB_128(x16_2, y16_2);
        x16_2 += 4, y16_2 += 4;

        LEO_FFTB_128(x16_3 + 1, y16_3 + 1);
        LEO_FFTB_128(x16_3, y16_3);
        x16_3 += 4, y16_3 += 4;

        bytes -= 64;
    } while (bytes > 0);
}

#endif // LEO_USE_VECTOR4_OPT


//------------------------------------------------------------------------------
// IFFT Operations

void ifft_butterfly(
    void * LEO_RESTRICT x, void * LEO_RESTRICT y,
    ffe_t log_m, uint64_t bytes)
{
#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
    {
        LEO_MUL_TABLES_256();

        const LEO_M256 clr_mask = _mm256_set1_epi8(0x0f);

        LEO_M256 * LEO_RESTRICT x32 = reinterpret_cast<LEO_M256 *>(x);
        LEO_M256 * LEO_RESTRICT y32 = reinterpret_cast<LEO_M256 *>(y);

        do
        {
#define LEO_IFFTB_256(x_ptr, y_ptr) { \
            LEO_M256 x_lo = _mm256_loadu_si256(x_ptr); \
            LEO_M256 y_lo = _mm256_loadu_si256(y_ptr); \
            y_lo = _mm256_xor_si256(y_lo, x_lo); \
            _mm256_storeu_si256(y_ptr, y_lo); \
            LEO_M256 data_0 = _mm256_and_si256(y_lo, clr_mask); \
            LEO_M256 data_1 = _mm256_srli_epi64(y_lo, 4); \
            data_1 = _mm256_and_si256(data_1, clr_mask); \
            LEO_M256 prod_lo = _mm256_shuffle_epi8(T0_lo, data_0); \
            prod_lo = _mm256_xor_si256(prod_lo, _mm256_shuffle_epi8(T1_lo, data_1)); \
            LEO_M256 prod_hi = _mm256_shuffle_epi8(T0_hi, data_0); \
            prod_hi = _mm256_xor_si256(prod_hi, _mm256_shuffle_epi8(T1_hi, data_1)); \
            LEO_M256 x_hi = _mm256_loadu_si256(x_ptr + 1); \
            LEO_M256 y_hi = _mm256_loadu_si256(y_ptr + 1); \
            y_hi = _mm256_xor_si256(y_hi, x_hi); \
            _mm256_storeu_si256(y_ptr + 1, y_hi); \
            data_0 = _mm256_and_si256(y_hi, clr_mask); \
            data_1 = _mm256_srli_epi64(y_hi, 4); \
            data_1 = _mm256_and_si256(data_1, clr_mask); \
            prod_lo = _mm256_xor_si256(prod_lo, _mm256_shuffle_epi8(T2_lo, data_0)); \
            prod_lo = _mm256_xor_si256(prod_lo, _mm256_shuffle_epi8(T3_lo, data_1)); \
            prod_hi = _mm256_xor_si256(prod_hi, _mm256_shuffle_epi8(T2_hi, data_0)); \
            prod_hi = _mm256_xor_si256(prod_hi, _mm256_shuffle_epi8(T3_hi, data_1)); \
            x_lo = _mm256_xor_si256(prod_lo, x_lo); \
            _mm256_storeu_si256(x_ptr, x_lo); \
            x_hi = _mm256_xor_si256(prod_hi, x_hi); \
            _mm256_storeu_si256(x_ptr + 1, x_hi); }

            LEO_IFFTB_256(x32, y32);
            y32 += 2, x32 += 2;

            bytes -= 64;
        } while (bytes > 0);

        return;
    }
#endif // LEO_TRY_AVX2

    LEO_MUL_TABLES_128();

    const LEO_M128 clr_mask = _mm_set1_epi8(0x0f);

    LEO_M128 * LEO_RESTRICT x16 = reinterpret_cast<LEO_M128 *>(x);
    LEO_M128 * LEO_RESTRICT y16 = reinterpret_cast<LEO_M128 *>(y);

    do
    {
#define LEO_IFFTB_128(x_ptr, y_ptr) { \
            LEO_M128 x_lo = _mm_loadu_si128(x_ptr); \
            LEO_M128 y_lo = _mm_loadu_si128(y_ptr); \
            y_lo = _mm_xor_si128(y_lo, x_lo); \
            _mm_storeu_si128(y_ptr, y_lo); \
            LEO_M128 data_0 = _mm_and_si128(y_lo, clr_mask); \
            LEO_M128 data_1 = _mm_srli_epi64(y_lo, 4); \
            data_1 = _mm_and_si128(data_1, clr_mask); \
            LEO_M128 prod_lo = _mm_shuffle_epi8(T0_lo, data_0); \
            prod_lo = _mm_xor_si128(prod_lo, _mm_shuffle_epi8(T1_lo, data_1)); \
            LEO_M128 prod_hi = _mm_shuffle_epi8(T0_hi, data_0); \
            prod_hi = _mm_xor_si128(prod_hi, _mm_shuffle_epi8(T1_hi, data_1)); \
            LEO_M128 x_hi = _mm_loadu_si128(x_ptr + 2); \
            LEO_M128 y_hi = _mm_loadu_si128(y_ptr + 2); \
            y_hi = _mm_xor_si128(y_hi, x_hi); \
            _mm_storeu_si128(y_ptr + 2, y_hi); \
            data_0 = _mm_and_si128(y_hi, clr_mask); \
            data_1 = _mm_srli_epi64(y_hi, 4); \
            data_1 = _mm_and_si128(data_1, clr_mask); \
            prod_lo = _mm_xor_si128(prod_lo, _mm_shuffle_epi8(T2_lo, data_0)); \
            prod_lo = _mm_xor_si128(prod_lo, _mm_shuffle_epi8(T3_lo, data_1)); \
            prod_hi = _mm_xor_si128(prod_hi, _mm_shuffle_epi8(T2_hi, data_0)); \
            prod_hi = _mm_xor_si128(prod_hi, _mm_shuffle_epi8(T3_hi, data_1)); \
            x_lo = _mm_xor_si128(prod_lo, x_lo); \
            _mm_storeu_si128(x_ptr, x_lo); \
            x_hi = _mm_xor_si128(prod_hi, x_hi); \
            _mm_storeu_si128(x_ptr + 2, x_hi); }

        LEO_IFFTB_128(x16 + 1, y16 + 1);
        LEO_IFFTB_128(x16, y16);
        x16 += 4, y16 += 4;

        bytes -= 64;
    } while (bytes > 0);
}

#ifdef LEO_USE_VECTOR4_OPT

void ifft_butterfly4(
    void * LEO_RESTRICT x_0, void * LEO_RESTRICT y_0,
    void * LEO_RESTRICT x_1, void * LEO_RESTRICT y_1,
    void * LEO_RESTRICT x_2, void * LEO_RESTRICT y_2,
    void * LEO_RESTRICT x_3, void * LEO_RESTRICT y_3,
    ffe_t log_m, uint64_t bytes)
{
#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
    {
        LEO_MUL_TABLES_256();

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
            LEO_IFFTB_256(x32_0, y32_0);
            y32_0 += 2, x32_0 += 2;

            LEO_IFFTB_256(x32_1, y32_1);
            y32_1 += 2, x32_1 += 2;

            LEO_IFFTB_256(x32_2, y32_2);
            y32_2 += 2, x32_2 += 2;

            LEO_IFFTB_256(x32_3, y32_3);
            y32_3 += 2, x32_3 += 2;

            bytes -= 64;
        } while (bytes > 0);

        return;
    }
#endif // LEO_TRY_AVX2

    LEO_MUL_TABLES_128();

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
        LEO_IFFTB_128(x16_0 + 1, y16_0 + 1);
        LEO_IFFTB_128(x16_0, y16_0);
        x16_0 += 4, y16_0 += 4;

        LEO_IFFTB_128(x16_1 + 1, y16_1 + 1);
        LEO_IFFTB_128(x16_1, y16_1);
        x16_1 += 4, y16_1 += 4;

        LEO_IFFTB_128(x16_2 + 1, y16_2 + 1);
        LEO_IFFTB_128(x16_2, y16_2);
        x16_2 += 4, y16_2 += 4;

        LEO_IFFTB_128(x16_3 + 1, y16_3 + 1);
        LEO_IFFTB_128(x16_3, y16_3);
        x16_3 += 4, y16_3 += 4;

        bytes -= 64;
    } while (bytes > 0);
}

#endif // LEO_USE_VECTOR4_OPT


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

    for (unsigned i = 0; i < kModulus; ++i)
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

void VectorIFFTButterfly(
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
        ifft_butterfly4(
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
        ifft_butterfly(x[i], y[i], log_m, bytes);
}


//------------------------------------------------------------------------------
// Reed-Solomon Encode

void ReedSolomonEncode(
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned m,
    const void* const * data,
    void** work)
{
    // work <- data

    // TBD: Unroll first loop to eliminate this
    unsigned first_end = m;
    if (original_count < m)
    {
        first_end = original_count;
        for (unsigned i = original_count; i < m; ++i)
            memset(work[i], 0, buffer_bytes);
    }
    for (unsigned i = 0; i < first_end; ++i)
        memcpy(work[i], data[i], buffer_bytes);

    // work <- IFFT(data, m, m)

    for (unsigned width = 1; width < m; width <<= 1)
    {
        const unsigned range = width << 1;
        const ffe_t* skewLUT = FFTSkew + width + m - 1;

#ifdef LEO_SCHEDULE_OPT
        for (unsigned j = 0; j < first_end; j += range)
#else
        for (unsigned j = 0; j < m; j += range)
#endif
        {
            VectorIFFTButterfly(
                buffer_bytes,
                width,
                work + j,
                work + j + width,
                skewLUT[j]);
        }
    }

    const unsigned last_count = original_count % m;
    if (m >= original_count)
        goto skip_body;

    for (unsigned i = m; i + m <= original_count; i += m)
    {
        // temp <- data + i

        data += m;
        void** temp = work + m;

        // TBD: Unroll first loop to eliminate this
        for (unsigned j = 0; j < m; ++j)
            memcpy(temp[j], data[j], buffer_bytes);

        // temp <- IFFT(temp, m, m + i)

        const ffe_t* skewLUT = FFTSkew + m + i - 1;

        for (unsigned width = 1; width < m; width <<= 1)
        {
            const unsigned range = width << 1;

            for (unsigned j = width; j < m; j += range)
            {
                VectorIFFTButterfly(
                    buffer_bytes,
                    width,
                    temp + j - width,
                    temp + j,
                    skewLUT[j]);
            }
        }

        // work <- work XOR temp

        // TBD: Unroll last loop to eliminate this
        VectorXOR(
            buffer_bytes,
            m,
            work,
            temp);
    }

    if (last_count != 0)
    {
        const unsigned i = original_count - last_count;

        // temp <- data + i

        data += m;
        void** temp = work + m;

        for (unsigned j = 0; j < last_count; ++j)
            memcpy(temp[j], data[j], buffer_bytes);
        for (unsigned j = last_count; j < m; ++j)
            memset(temp[j], 0, buffer_bytes);

        // temp <- IFFT(temp, m, m + i)

        for (unsigned width = 1, shift = 1; width < m; width <<= 1, ++shift)
        {
            const unsigned range = width << 1;
            const ffe_t* skewLUT = FFTSkew + width + m + i - 1;

#ifdef LEO_SCHEDULE_OPT
            // Calculate stop considering that the right is all zeroes
            const unsigned stop = ((last_count + range - 1) >> shift) << shift;
            for (unsigned j = 0; j < stop; j += range)
#else
            for (unsigned j = 0; j < m; j += range)
#endif
            {
                VectorIFFTButterfly(
                    buffer_bytes,
                    width,
                    temp + j,
                    temp + j + width,
                    skewLUT[j]);
            }
        }

        // work <- work XOR temp

        // TBD: Unroll last loop to eliminate this
        VectorXOR(
            buffer_bytes,
            m,
            work,
            temp);
    }

skip_body:

    // work <- FFT(work, m, 0)

    for (unsigned width = (m >> 1); width > 0; width >>= 1)
    {
        const ffe_t* skewLUT = FFTSkew + width - 1;
        const unsigned range = width << 1;

#ifdef LEO_SCHEDULE_OPT
        for (unsigned j = 0; j < recovery_count; j += range)
#else
        for (unsigned j = 0; j < m; j += range)
#endif
        {
            VectorFFTButterfly(
                buffer_bytes,
                width,
                work + j,
                work + j + width,
                skewLUT[j]);
        }
    }
}


//------------------------------------------------------------------------------
// ErrorBitfield

#ifdef LEO_ERROR_BITFIELD_OPT

// Used in decoding to decide which final FFT operations to perform
class ErrorBitfield
{
    static const unsigned kWordMips = 5;
    static const unsigned kWords = kOrder / 64;
    uint64_t Words[kWordMips][kWords] = {};

    static const unsigned kBigMips = 6;
    static const unsigned kBigWords = (kWords + 63) / 64;
    uint64_t BigWords[kBigMips][kBigWords] = {};

    static const unsigned kBiggestMips = 4;
    uint64_t BiggestWords[kBiggestMips] = {};

public:
    LEO_FORCE_INLINE void Set(unsigned i)
    {
        Words[0][i / 64] |= (uint64_t)1 << (i % 64);
    }

    void Prepare();

    LEO_FORCE_INLINE bool IsNeeded(unsigned mip_level, unsigned bit) const
    {
        if (mip_level >= 16)
            return true;
        if (mip_level >= 12)
        {
            bit /= 4096;
            return 0 != (BiggestWords[mip_level - 12] & ((uint64_t)1 << bit));
        }
        if (mip_level >= 6)
        {
            bit /= 64;
            return 0 != (BigWords[mip_level - 6][bit / 64] & ((uint64_t)1 << (bit % 64)));
        }
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

        for (unsigned j = 1, bits = 2; j < kWordMips; ++j, bits <<= 1)
        {
            const uint64_t hi2lo_j = w_i | ((w_i & kHiMasks[j]) >> bits);
            const uint64_t lo2hi_j = ((w_i & (kHiMasks[j] >> bits)) << bits);
            Words[j][i] = w_i = hi2lo_j | lo2hi_j;
        }
    }

    for (unsigned i = 0; i < kBigWords; ++i)
    {
        uint64_t w_i = 0;
        uint64_t bit = 1;
        const uint64_t* src = &Words[kWordMips - 1][i * 64];
        for (unsigned j = 0; j < 64; ++j, bit <<= 1)
        {
            const uint64_t w = src[j];
            w_i |= (w | (w >> 32) | (w << 32)) & bit;
        }
        BigWords[0][i] = w_i;

        for (unsigned j = 1, bits = 1; j < kBigMips; ++j, bits <<= 1)
        {
            const uint64_t hi2lo_j = w_i | ((w_i & kHiMasks[j - 1]) >> bits);
            const uint64_t lo2hi_j = ((w_i & (kHiMasks[j - 1] >> bits)) << bits);
            BigWords[j][i] = w_i = hi2lo_j | lo2hi_j;
        }
    }

    uint64_t w_i = 0;
    uint64_t bit = 1;
    const uint64_t* src = &BigWords[kBigMips - 1][0];
    for (unsigned j = 0; j < kBigWords; ++j, bit <<= 1)
    {
        const uint64_t w = src[j];
        w_i |= (w | (w >> 32) | (w << 32)) & bit;
    }
    BiggestWords[0] = w_i;

    for (unsigned j = 1, bits = 1; j < kBiggestMips; ++j, bits <<= 1)
    {
        const uint64_t hi2lo_j = w_i | ((w_i & kHiMasks[j - 1]) >> bits);
        const uint64_t lo2hi_j = ((w_i & (kHiMasks[j - 1] >> bits)) << bits);
        BiggestWords[j] = w_i = hi2lo_j | lo2hi_j;
    }
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
    const void* const * const original, // original_count entries
    const void* const * const recovery, // recovery_count entries
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

    const unsigned input_count = m + original_count;
    unsigned mip_level = 0;

    for (unsigned width = 1; width < n; width <<= 1, ++mip_level)
    {
        const unsigned range = width << 1;

        for (unsigned j = width; j < n; j += range)
        {
            VectorIFFTButterfly(
                buffer_bytes,
                width,
                work + j - width,
                work + j,
                FFTSkew[j - 1]);
        }
    }

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

    if (!CpuHasSSSE3)
        return false;

    InitializeLogarithmTables();
    InitializeMultiplyTables();
    FFTInitialize();

    IsInitialized = true;
    return true;
}


}} // namespace leopard::ff16

#endif // LEO_HAS_FF16
