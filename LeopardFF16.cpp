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

#ifdef _MSC_VER
    #pragma warning(disable: 4752) // found Intel(R) Advanced Vector Extensions; consider using /arch:AVX
#endif

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
#pragma omp parallel for
        for (int r = 0; r < (int)m_truncated; r += dist4)
        {
            // For each set of dist elements:
            const int i_end = r + dist;
            for (int i = r; i < i_end; ++i)
                FWHT_4(data + i, dist);
        }
    }

    // If there is one layer left:
    if (dist < m)
#pragma omp parallel for
        for (int i = 0; i < (int)dist; ++i)
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

struct Multiply128LUT_t
{
    LEO_M128 Lo[4];
    LEO_M128 Hi[4];
};

static const Multiply128LUT_t* Multiply128LUT = nullptr;

#define LEO_MUL_TABLES_128(table, log_m) \
    const LEO_M128 T0_lo_##table = _mm_loadu_si128(&Multiply128LUT[log_m].Lo[0]); \
    const LEO_M128 T1_lo_##table = _mm_loadu_si128(&Multiply128LUT[log_m].Lo[1]); \
    const LEO_M128 T2_lo_##table = _mm_loadu_si128(&Multiply128LUT[log_m].Lo[2]); \
    const LEO_M128 T3_lo_##table = _mm_loadu_si128(&Multiply128LUT[log_m].Lo[3]); \
    const LEO_M128 T0_hi_##table = _mm_loadu_si128(&Multiply128LUT[log_m].Hi[0]); \
    const LEO_M128 T1_hi_##table = _mm_loadu_si128(&Multiply128LUT[log_m].Hi[1]); \
    const LEO_M128 T2_hi_##table = _mm_loadu_si128(&Multiply128LUT[log_m].Hi[2]); \
    const LEO_M128 T3_hi_##table = _mm_loadu_si128(&Multiply128LUT[log_m].Hi[3]);

// 128-bit {prod_lo, prod_hi} = {value_lo, value_hi} * log_m
#define LEO_MUL_128(value_lo, value_hi, table) { \
            LEO_M128 data_1 = _mm_srli_epi64(value_lo, 4); \
            LEO_M128 data_0 = _mm_and_si128(value_lo, clr_mask); \
            data_1 = _mm_and_si128(data_1, clr_mask); \
            prod_lo = _mm_shuffle_epi8(T0_lo_##table, data_0); \
            prod_hi = _mm_shuffle_epi8(T0_hi_##table, data_0); \
            prod_lo = _mm_xor_si128(prod_lo, _mm_shuffle_epi8(T1_lo_##table, data_1)); \
            prod_hi = _mm_xor_si128(prod_hi, _mm_shuffle_epi8(T1_hi_##table, data_1)); \
            data_0 = _mm_and_si128(value_hi, clr_mask); \
            data_1 = _mm_srli_epi64(value_hi, 4); \
            data_1 = _mm_and_si128(data_1, clr_mask); \
            prod_lo = _mm_xor_si128(prod_lo, _mm_shuffle_epi8(T2_lo_##table, data_0)); \
            prod_hi = _mm_xor_si128(prod_hi, _mm_shuffle_epi8(T2_hi_##table, data_0)); \
            prod_lo = _mm_xor_si128(prod_lo, _mm_shuffle_epi8(T3_lo_##table, data_1)); \
            prod_hi = _mm_xor_si128(prod_hi, _mm_shuffle_epi8(T3_hi_##table, data_1)); }

// {x_lo, x_hi} ^= {y_lo, y_hi} * log_m
#define LEO_MULADD_128(x_lo, x_hi, y_lo, y_hi, table) { \
            LEO_M128 prod_lo, prod_hi; \
            LEO_MUL_128(y_lo, y_hi, table); \
            x_lo = _mm_xor_si128(x_lo, prod_lo); \
            x_hi = _mm_xor_si128(x_hi, prod_hi); }


#if defined(LEO_TRY_AVX2)

struct Multiply256LUT_t
{
    LEO_M256 Lo[4];
    LEO_M256 Hi[4];
};

static const Multiply256LUT_t* Multiply256LUT = nullptr;

#define LEO_MUL_TABLES_256(table, log_m) \
        const LEO_M256 T0_lo_##table = _mm256_loadu_si256(&Multiply256LUT[log_m].Lo[0]); \
        const LEO_M256 T1_lo_##table = _mm256_loadu_si256(&Multiply256LUT[log_m].Lo[1]); \
        const LEO_M256 T2_lo_##table = _mm256_loadu_si256(&Multiply256LUT[log_m].Lo[2]); \
        const LEO_M256 T3_lo_##table = _mm256_loadu_si256(&Multiply256LUT[log_m].Lo[3]); \
        const LEO_M256 T0_hi_##table = _mm256_loadu_si256(&Multiply256LUT[log_m].Hi[0]); \
        const LEO_M256 T1_hi_##table = _mm256_loadu_si256(&Multiply256LUT[log_m].Hi[1]); \
        const LEO_M256 T2_hi_##table = _mm256_loadu_si256(&Multiply256LUT[log_m].Hi[2]); \
        const LEO_M256 T3_hi_##table = _mm256_loadu_si256(&Multiply256LUT[log_m].Hi[3]);

// 256-bit {prod_lo, prod_hi} = {value_lo, value_hi} * log_m
#define LEO_MUL_256(value_lo, value_hi, table) { \
            LEO_M256 data_1 = _mm256_srli_epi64(value_lo, 4); \
            LEO_M256 data_0 = _mm256_and_si256(value_lo, clr_mask); \
            data_1 = _mm256_and_si256(data_1, clr_mask); \
            prod_lo = _mm256_shuffle_epi8(T0_lo_##table, data_0); \
            prod_hi = _mm256_shuffle_epi8(T0_hi_##table, data_0); \
            prod_lo = _mm256_xor_si256(prod_lo, _mm256_shuffle_epi8(T1_lo_##table, data_1)); \
            prod_hi = _mm256_xor_si256(prod_hi, _mm256_shuffle_epi8(T1_hi_##table, data_1)); \
            data_0 = _mm256_and_si256(value_hi, clr_mask); \
            data_1 = _mm256_srli_epi64(value_hi, 4); \
            data_1 = _mm256_and_si256(data_1, clr_mask); \
            prod_lo = _mm256_xor_si256(prod_lo, _mm256_shuffle_epi8(T2_lo_##table, data_0)); \
            prod_hi = _mm256_xor_si256(prod_hi, _mm256_shuffle_epi8(T2_hi_##table, data_0)); \
            prod_lo = _mm256_xor_si256(prod_lo, _mm256_shuffle_epi8(T3_lo_##table, data_1)); \
            prod_hi = _mm256_xor_si256(prod_hi, _mm256_shuffle_epi8(T3_hi_##table, data_1)); }

// {x_lo, x_hi} ^= {y_lo, y_hi} * log_m
#define LEO_MULADD_256(x_lo, x_hi, y_lo, y_hi, table) { \
            LEO_M256 prod_lo, prod_hi; \
            LEO_MUL_256(y_lo, y_hi, table); \
            x_lo = _mm256_xor_si256(x_lo, prod_lo); \
            x_hi = _mm256_xor_si256(x_hi, prod_hi); }

#endif // LEO_TRY_AVX2

// Stores the partial products of x * y at offset x + y * 65536
// Repeated accesses from the same y value are faster
struct Product16Table
{
    ffe_t LUT[4 * 16];
};
static const Product16Table* Multiply16LUT = nullptr;


// Reference version of muladd: x[] ^= y[] * log_m
static LEO_FORCE_INLINE void RefMulAdd(
    void* LEO_RESTRICT x,
    const void* LEO_RESTRICT y,
    ffe_t log_m,
    uint64_t bytes)
{
    const ffe_t* LEO_RESTRICT lut = Multiply16LUT[log_m].LUT;
    const uint8_t * LEO_RESTRICT y1 = reinterpret_cast<const uint8_t *>(y);
    uint8_t * LEO_RESTRICT x1 = reinterpret_cast<uint8_t *>(x);

    do
    {
        for (unsigned i = 0; i < 32; ++i)
        {
            const unsigned lo = y1[i];
            const unsigned hi = y1[i + 32];

            const ffe_t prod = \
                lut[(lo & 15)] ^ \
                lut[(lo >> 4) + 16] ^ \
                lut[(hi & 15) + 32] ^ \
                lut[(hi >> 4) + 48];

            x1[i] ^= (uint8_t)prod;
            x1[i + 32] ^= (uint8_t)(prod >> 8);
        }

        x1 += 64, y1 += 64;
        bytes -= 64;
    } while (bytes > 0);

}

// Reference version of mul: x[] = y[] * log_m
static LEO_FORCE_INLINE void RefMul(
    void* LEO_RESTRICT x,
    const void* LEO_RESTRICT y,
    ffe_t log_m,
    uint64_t bytes)
{
    const ffe_t* LEO_RESTRICT lut = Multiply16LUT[log_m].LUT;
    const uint8_t * LEO_RESTRICT y1 = reinterpret_cast<const uint8_t *>(y);
    uint8_t * LEO_RESTRICT x1 = reinterpret_cast<uint8_t *>(x);

    do
    {
        for (unsigned i = 0; i < 32; ++i)
        {
            const unsigned lo = y1[i];
            const unsigned hi = y1[i + 32];

            const ffe_t prod = \
                lut[(lo & 15)] ^ \
                lut[(lo >> 4) + 16] ^ \
                lut[(hi & 15) + 32] ^ \
                lut[(hi >> 4) + 48];

            x1[i] = (uint8_t)prod;
            x1[i + 32] = (uint8_t)(prod >> 8);
        }

        x1 += 64, y1 += 64;
        bytes -= 64;
    } while (bytes > 0);
}


static void InitializeMultiplyTables()
{
    // If we cannot use the PSHUFB instruction, generate Multiply8LUT:
    if (!CpuHasSSSE3)
    {
        Multiply16LUT = new Product16Table[65536];

        // For each log_m multiplicand:
#pragma omp parallel for
        for (int log_m = 0; log_m < (int)kOrder; ++log_m)
        {
            const Product16Table& lut = Multiply16LUT[log_m];

            for (unsigned nibble = 0, shift = 0; nibble < 4; ++nibble, shift += 4)
            {
                ffe_t* nibble_lut = (ffe_t*)&lut.LUT[nibble * 16];

                for (unsigned x_nibble = 0; x_nibble < 16; ++x_nibble)
                {
                    const ffe_t prod = MultiplyLog(x_nibble << shift, static_cast<ffe_t>(log_m));
                    nibble_lut[x_nibble] = prod;
                }
            }
        }

        return;
    }

#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
        Multiply256LUT = reinterpret_cast<const Multiply256LUT_t*>(SIMDSafeAllocate(sizeof(Multiply256LUT_t) * kOrder));
    else
#endif // LEO_TRY_AVX2
        Multiply128LUT = reinterpret_cast<const Multiply128LUT_t*>(SIMDSafeAllocate(sizeof(Multiply128LUT_t) * kOrder));

    // For each value we could multiply by:
#pragma omp parallel for
    for (int log_m = 0; log_m < (int)kOrder; ++log_m)
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
#if defined(LEO_TRY_AVX2)
            if (!CpuHasAVX2)
#endif // LEO_TRY_AVX2
            {
                _mm_storeu_si128((LEO_M128*)&Multiply128LUT[log_m].Lo[i], value_lo);
                _mm_storeu_si128((LEO_M128*)&Multiply128LUT[log_m].Hi[i], value_hi);
            }

            // Store in 256-bit wide table
#if defined(LEO_TRY_AVX2)
            if (CpuHasAVX2)
            {
                _mm256_storeu_si256((LEO_M256*)&Multiply256LUT[log_m].Lo[i],
                    _mm256_broadcastsi128_si256(value_lo));
                _mm256_storeu_si256((LEO_M256*)&Multiply256LUT[log_m].Hi[i],
                    _mm256_broadcastsi128_si256(value_hi));
            }
#endif // LEO_TRY_AVX2
        }
    }
}


static void mul_mem(
    void * LEO_RESTRICT x, const void * LEO_RESTRICT y,
    ffe_t log_m, uint64_t bytes)
{
#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
    {
        LEO_MUL_TABLES_256(0, log_m);

        const LEO_M256 clr_mask = _mm256_set1_epi8(0x0f);

        LEO_M256 * LEO_RESTRICT x32 = reinterpret_cast<LEO_M256 *>(x);
        const LEO_M256 * LEO_RESTRICT y32 = reinterpret_cast<const LEO_M256 *>(y);

        do
        {
#define LEO_MUL_256_LS(x_ptr, y_ptr) { \
            const LEO_M256 data_lo = _mm256_loadu_si256(y_ptr); \
            const LEO_M256 data_hi = _mm256_loadu_si256(y_ptr + 1); \
            LEO_M256 prod_lo, prod_hi; \
            LEO_MUL_256(data_lo, data_hi, 0); \
            _mm256_storeu_si256(x_ptr, prod_lo); \
            _mm256_storeu_si256(x_ptr + 1, prod_hi); }

            LEO_MUL_256_LS(x32, y32);
            y32 += 2, x32 += 2;

            bytes -= 64;
        } while (bytes > 0);

        return;
    }
#endif // LEO_TRY_AVX2

#if defined(LEO_TRY_SSSE3)
    if (CpuHasSSSE3)
    {
        LEO_MUL_TABLES_128(0, log_m);

        const LEO_M128 clr_mask = _mm_set1_epi8(0x0f);

        LEO_M128 * LEO_RESTRICT x16 = reinterpret_cast<LEO_M128 *>(x);
        const LEO_M128 * LEO_RESTRICT y16 = reinterpret_cast<const LEO_M128 *>(y);

        do
        {
#define LEO_MUL_128_LS(x_ptr, y_ptr) { \
                const LEO_M128 data_lo = _mm_loadu_si128(y_ptr); \
                const LEO_M128 data_hi = _mm_loadu_si128(y_ptr + 2); \
                LEO_M128 prod_lo, prod_hi; \
                LEO_MUL_128(data_lo, data_hi, 0); \
                _mm_storeu_si128(x_ptr, prod_lo); \
                _mm_storeu_si128(x_ptr + 2, prod_hi); }

            LEO_MUL_128_LS(x16 + 1, y16 + 1);
            LEO_MUL_128_LS(x16, y16);
            x16 += 4, y16 += 4;

            bytes -= 64;
        } while (bytes > 0);

        return;
    }
#endif // LEO_TRY_SSSE3

    RefMul(x, y, log_m, bytes);
}

// x[] ^= y[] * log_m for complete 64-byte ALTMAP tiles.
static void muladd_mem(
    void * LEO_RESTRICT x, const void * LEO_RESTRICT y,
    ffe_t log_m, uint64_t bytes)
{
#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
    {
        LEO_MUL_TABLES_256(0, log_m);
        const LEO_M256 clr_mask = _mm256_set1_epi8(0x0f);
        LEO_M256 * LEO_RESTRICT x32 = reinterpret_cast<LEO_M256 *>(x);
        const LEO_M256 * LEO_RESTRICT y32 = reinterpret_cast<const LEO_M256 *>(y);

        do
        {
            LEO_M256 x_lo = _mm256_loadu_si256(x32);
            LEO_M256 x_hi = _mm256_loadu_si256(x32 + 1);
            const LEO_M256 y_lo = _mm256_loadu_si256(y32);
            const LEO_M256 y_hi = _mm256_loadu_si256(y32 + 1);
            LEO_MULADD_256(x_lo, x_hi, y_lo, y_hi, 0);
            _mm256_storeu_si256(x32, x_lo);
            _mm256_storeu_si256(x32 + 1, x_hi);
            x32 += 2, y32 += 2;
            bytes -= 64;
        } while (bytes > 0);
        return;
    }
#endif // LEO_TRY_AVX2

#if defined(LEO_TRY_SSSE3)
    if (CpuHasSSSE3)
    {
        LEO_MUL_TABLES_128(0, log_m);
        const LEO_M128 clr_mask = _mm_set1_epi8(0x0f);
        LEO_M128 * LEO_RESTRICT x16 = reinterpret_cast<LEO_M128 *>(x);
        const LEO_M128 * LEO_RESTRICT y16 = reinterpret_cast<const LEO_M128 *>(y);

        do
        {
            for (unsigned i = 0; i < 2; ++i)
            {
                LEO_M128 x_lo = _mm_loadu_si128(x16 + i);
                LEO_M128 x_hi = _mm_loadu_si128(x16 + i + 2);
                const LEO_M128 y_lo = _mm_loadu_si128(y16 + i);
                const LEO_M128 y_hi = _mm_loadu_si128(y16 + i + 2);
                LEO_MULADD_128(x_lo, x_hi, y_lo, y_hi, 0);
                _mm_storeu_si128(x16 + i, x_lo);
                _mm_storeu_si128(x16 + i + 2, x_hi);
            }
            x16 += 4, y16 += 4;
            bytes -= 64;
        } while (bytes > 0);
        return;
    }
#endif // LEO_TRY_SSSE3

    RefMulAdd(x, y, log_m, bytes);
}


ffe_t MultiplyElements(ffe_t a, ffe_t b)
{
    if (a == 0 || b == 0)
        return 0;
    return MultiplyLog(a, LogLUT[b]);
}


ffe_t InverseElement(ffe_t value)
{
    LEO_DEBUG_ASSERT(value != 0);
    if (value == 0)
        return 0;
    return ExpLUT[kModulus - LogLUT[value]];
}


ffe_t ElementLog(ffe_t value)
{
    LEO_DEBUG_ASSERT(value != 0);
    return LogLUT[value];
}


void MultiplyBytes(
    void* destination,
    const void* source,
    ffe_t multiplier_log,
    uint64_t byte_count)
{
    const uint64_t complete = byte_count & ~static_cast<uint64_t>(63);
    if (complete != 0)
        mul_mem(destination, source, multiplier_log, complete);

    const uint64_t residual = byte_count - complete;
    LEO_DEBUG_ASSERT((residual & 1) == 0);
    const uint64_t symbols = residual / 2;
    uint8_t* output = reinterpret_cast<uint8_t*>(destination);
    const uint8_t* input = reinterpret_cast<const uint8_t*>(source);
    for (uint64_t i = 0; i < symbols; ++i)
    {
        const ffe_t value = static_cast<ffe_t>(input[complete + i] |
            (static_cast<unsigned>(input[complete + symbols + i]) << 8));
        const ffe_t product = MultiplyLog(value, multiplier_log);
        output[complete + i] = static_cast<uint8_t>(product);
        output[complete + symbols + i] = static_cast<uint8_t>(product >> 8);
    }
}


void MultiplyAddBytes(
    void* destination,
    const void* source,
    ffe_t multiplier_log,
    uint64_t byte_count)
{
    const uint64_t complete = byte_count & ~static_cast<uint64_t>(63);
    if (complete != 0)
        muladd_mem(destination, source, multiplier_log, complete);

    const uint64_t residual = byte_count - complete;
    LEO_DEBUG_ASSERT((residual & 1) == 0);
    const uint64_t symbols = residual / 2;
    uint8_t* output = reinterpret_cast<uint8_t*>(destination);
    const uint8_t* input = reinterpret_cast<const uint8_t*>(source);
    for (uint64_t i = 0; i < symbols; ++i)
    {
        const ffe_t value = static_cast<ffe_t>(input[complete + i] |
            (static_cast<unsigned>(input[complete + symbols + i]) << 8));
        const ffe_t product = MultiplyLog(value, multiplier_log);
        output[complete + i] ^= static_cast<uint8_t>(product);
        output[complete + symbols + i] ^= static_cast<uint8_t>(product >> 8);
    }
}


//------------------------------------------------------------------------------
// FFT

// Twisted factors used in FFT.  The transform kernels conventionally bias the
// pointer by one so that logical factor i is addressed as skewLUT[i + 1].
// Keep a real sentinel element before the logical table.  The transform base is
// passed as FFTSkewStorage directly; the historical one-before-begin pointer is
// never formed, and the sentinel itself is never consumed as a multiplier.
static ffe_t FFTSkewStorage[kOrder];
static ffe_t* const FFTSkew = FFTSkewStorage + 1;

// log(s_j(v_j)) for the normalized LCH basis.  General p_i values are the
// sum of the entries selected by the bits of i.
static ffe_t LchBasisNormalizerLog[kBits];

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

    for (unsigned bit = 0; bit < kBits; ++bit)
    {
        const unsigned width = 1UL << bit;
        ffe_t factor_log = 0;
        for (unsigned i = 0; i < width; ++i)
            factor_log = AddMod(factor_log, LogLUT[width ^ i]);
        LchBasisNormalizerLog[bit] = factor_log;
    }

    // Precalculate FWHT(Log[i]):

    for (unsigned i = 0; i < kOrder; ++i)
        LogWalsh[i] = LogLUT[i];
    LogWalsh[0] = 0;

    FWHT(LogWalsh, kOrder, kOrder);
}

/*
    Decimation in time IFFT:

    The decimation in time IFFT algorithm allows us to unroll 2 layers at a time,
    performing calculations on local registers and faster cache memory.

    Each ^___^ below indicates a butterfly between the associated indices.

    The ifft_butterfly(x, y) operation:

        y[] ^= x[]
        if (log_m != kModulus)
            x[] ^= exp(log(y[]) + log_m)

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

// 2-way butterfly
static void IFFT_DIT2(
    void * LEO_RESTRICT x, void * LEO_RESTRICT y,
    ffe_t log_m, uint64_t bytes)
{
#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
    {
        LEO_MUL_TABLES_256(0, log_m);

        const LEO_M256 clr_mask = _mm256_set1_epi8(0x0f);

        LEO_M256 * LEO_RESTRICT x32 = reinterpret_cast<LEO_M256 *>(x);
        LEO_M256 * LEO_RESTRICT y32 = reinterpret_cast<LEO_M256 *>(y);

        do
        {
#define LEO_IFFTB_256(x_ptr, y_ptr) { \
            LEO_M256 x_lo = _mm256_loadu_si256(x_ptr); \
            LEO_M256 x_hi = _mm256_loadu_si256(x_ptr + 1); \
            LEO_M256 y_lo = _mm256_loadu_si256(y_ptr); \
            LEO_M256 y_hi = _mm256_loadu_si256(y_ptr + 1); \
            y_lo = _mm256_xor_si256(y_lo, x_lo); \
            y_hi = _mm256_xor_si256(y_hi, x_hi); \
            _mm256_storeu_si256(y_ptr, y_lo); \
            _mm256_storeu_si256(y_ptr + 1, y_hi); \
            LEO_MULADD_256(x_lo, x_hi, y_lo, y_hi, 0); \
            _mm256_storeu_si256(x_ptr, x_lo); \
            _mm256_storeu_si256(x_ptr + 1, x_hi); }

            LEO_IFFTB_256(x32, y32);
            y32 += 2, x32 += 2;

            bytes -= 64;
        } while (bytes > 0);

        return;
    }
#endif // LEO_TRY_AVX2

#if defined(LEO_TRY_SSSE3)
    if (CpuHasSSSE3)
    {
        LEO_MUL_TABLES_128(0, log_m);

        const LEO_M128 clr_mask = _mm_set1_epi8(0x0f);

        LEO_M128 * LEO_RESTRICT x16 = reinterpret_cast<LEO_M128 *>(x);
        LEO_M128 * LEO_RESTRICT y16 = reinterpret_cast<LEO_M128 *>(y);

        do
        {
#define LEO_IFFTB_128(x_ptr, y_ptr) { \
                LEO_M128 x_lo = _mm_loadu_si128(x_ptr); \
                LEO_M128 x_hi = _mm_loadu_si128(x_ptr + 2); \
                LEO_M128 y_lo = _mm_loadu_si128(y_ptr); \
                LEO_M128 y_hi = _mm_loadu_si128(y_ptr + 2); \
                y_lo = _mm_xor_si128(y_lo, x_lo); \
                y_hi = _mm_xor_si128(y_hi, x_hi); \
                _mm_storeu_si128(y_ptr, y_lo); \
                _mm_storeu_si128(y_ptr + 2, y_hi); \
                LEO_MULADD_128(x_lo, x_hi, y_lo, y_hi, 0); \
                _mm_storeu_si128(x_ptr, x_lo); \
                _mm_storeu_si128(x_ptr + 2, x_hi); }

            LEO_IFFTB_128(x16 + 1, y16 + 1);
            LEO_IFFTB_128(x16, y16);
            x16 += 4, y16 += 4;

            bytes -= 64;
        } while (bytes > 0);

        return;
    }
#endif // LEO_TRY_SSSE3

    // Reference version:
    xor_mem(y, x, bytes);
    RefMulAdd(x, y, log_m, bytes);
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
        LEO_MUL_TABLES_256(01, log_m01);
        LEO_MUL_TABLES_256(23, log_m23);
        LEO_MUL_TABLES_256(02, log_m02);

        const LEO_M256 clr_mask = _mm256_set1_epi8(0x0f);

        LEO_M256 * LEO_RESTRICT work0 = reinterpret_cast<LEO_M256 *>(work[0]);
        LEO_M256 * LEO_RESTRICT work1 = reinterpret_cast<LEO_M256 *>(work[dist]);
        LEO_M256 * LEO_RESTRICT work2 = reinterpret_cast<LEO_M256 *>(work[dist * 2]);
        LEO_M256 * LEO_RESTRICT work3 = reinterpret_cast<LEO_M256 *>(work[dist * 3]);

        do
        {
            LEO_M256 work_reg_lo_0 = _mm256_loadu_si256(work0);
            LEO_M256 work_reg_hi_0 = _mm256_loadu_si256(work0 + 1);
            LEO_M256 work_reg_lo_1 = _mm256_loadu_si256(work1);
            LEO_M256 work_reg_hi_1 = _mm256_loadu_si256(work1 + 1);

            // First layer:
            work_reg_lo_1 = _mm256_xor_si256(work_reg_lo_0, work_reg_lo_1);
            work_reg_hi_1 = _mm256_xor_si256(work_reg_hi_0, work_reg_hi_1);
            if (log_m01 != kModulus)
                LEO_MULADD_256(work_reg_lo_0, work_reg_hi_0, work_reg_lo_1, work_reg_hi_1, 01);

            LEO_M256 work_reg_lo_2 = _mm256_loadu_si256(work2);
            LEO_M256 work_reg_hi_2 = _mm256_loadu_si256(work2 + 1);
            LEO_M256 work_reg_lo_3 = _mm256_loadu_si256(work3);
            LEO_M256 work_reg_hi_3 = _mm256_loadu_si256(work3 + 1);

            work_reg_lo_3 = _mm256_xor_si256(work_reg_lo_2, work_reg_lo_3);
            work_reg_hi_3 = _mm256_xor_si256(work_reg_hi_2, work_reg_hi_3);
            if (log_m23 != kModulus)
                LEO_MULADD_256(work_reg_lo_2, work_reg_hi_2, work_reg_lo_3, work_reg_hi_3, 23);

            // Second layer:
            work_reg_lo_2 = _mm256_xor_si256(work_reg_lo_0, work_reg_lo_2);
            work_reg_hi_2 = _mm256_xor_si256(work_reg_hi_0, work_reg_hi_2);
            work_reg_lo_3 = _mm256_xor_si256(work_reg_lo_1, work_reg_lo_3);
            work_reg_hi_3 = _mm256_xor_si256(work_reg_hi_1, work_reg_hi_3);
            if (log_m02 != kModulus)
            {
                LEO_MULADD_256(work_reg_lo_0, work_reg_hi_0, work_reg_lo_2, work_reg_hi_2, 02);
                LEO_MULADD_256(work_reg_lo_1, work_reg_hi_1, work_reg_lo_3, work_reg_hi_3, 02);
            }

            _mm256_storeu_si256(work0, work_reg_lo_0);
            _mm256_storeu_si256(work0 + 1, work_reg_hi_0);
            _mm256_storeu_si256(work1, work_reg_lo_1);
            _mm256_storeu_si256(work1 + 1, work_reg_hi_1);
            _mm256_storeu_si256(work2, work_reg_lo_2);
            _mm256_storeu_si256(work2 + 1, work_reg_hi_2);
            _mm256_storeu_si256(work3, work_reg_lo_3);
            _mm256_storeu_si256(work3 + 1, work_reg_hi_3);

            work0 += 2, work1 += 2, work2 += 2, work3 += 2;

            bytes -= 64;
        } while (bytes > 0);

        return;
    }

#endif // LEO_TRY_AVX2

#if defined(LEO_TRY_SSSE3)
    if (CpuHasSSSE3)
    {
        LEO_MUL_TABLES_128(01, log_m01);
        LEO_MUL_TABLES_128(23, log_m23);
        LEO_MUL_TABLES_128(02, log_m02);

        const LEO_M128 clr_mask = _mm_set1_epi8(0x0f);

        LEO_M128 * LEO_RESTRICT work0 = reinterpret_cast<LEO_M128 *>(work[0]);
        LEO_M128 * LEO_RESTRICT work1 = reinterpret_cast<LEO_M128 *>(work[dist]);
        LEO_M128 * LEO_RESTRICT work2 = reinterpret_cast<LEO_M128 *>(work[dist * 2]);
        LEO_M128 * LEO_RESTRICT work3 = reinterpret_cast<LEO_M128 *>(work[dist * 3]);

        do
        {
            for (unsigned i = 0; i < 2; ++i)
            {
                LEO_M128 work_reg_lo_0 = _mm_loadu_si128(work0);
                LEO_M128 work_reg_hi_0 = _mm_loadu_si128(work0 + 2);
                LEO_M128 work_reg_lo_1 = _mm_loadu_si128(work1);
                LEO_M128 work_reg_hi_1 = _mm_loadu_si128(work1 + 2);

                // First layer:
                work_reg_lo_1 = _mm_xor_si128(work_reg_lo_0, work_reg_lo_1);
                work_reg_hi_1 = _mm_xor_si128(work_reg_hi_0, work_reg_hi_1);
                if (log_m01 != kModulus)
                    LEO_MULADD_128(work_reg_lo_0, work_reg_hi_0, work_reg_lo_1, work_reg_hi_1, 01);

                LEO_M128 work_reg_lo_2 = _mm_loadu_si128(work2);
                LEO_M128 work_reg_hi_2 = _mm_loadu_si128(work2 + 2);
                LEO_M128 work_reg_lo_3 = _mm_loadu_si128(work3);
                LEO_M128 work_reg_hi_3 = _mm_loadu_si128(work3 + 2);

                work_reg_lo_3 = _mm_xor_si128(work_reg_lo_2, work_reg_lo_3);
                work_reg_hi_3 = _mm_xor_si128(work_reg_hi_2, work_reg_hi_3);
                if (log_m23 != kModulus)
                    LEO_MULADD_128(work_reg_lo_2, work_reg_hi_2, work_reg_lo_3, work_reg_hi_3, 23);

                // Second layer:
                work_reg_lo_2 = _mm_xor_si128(work_reg_lo_0, work_reg_lo_2);
                work_reg_hi_2 = _mm_xor_si128(work_reg_hi_0, work_reg_hi_2);
                work_reg_lo_3 = _mm_xor_si128(work_reg_lo_1, work_reg_lo_3);
                work_reg_hi_3 = _mm_xor_si128(work_reg_hi_1, work_reg_hi_3);
                if (log_m02 != kModulus)
                {
                    LEO_MULADD_128(work_reg_lo_0, work_reg_hi_0, work_reg_lo_2, work_reg_hi_2, 02);
                    LEO_MULADD_128(work_reg_lo_1, work_reg_hi_1, work_reg_lo_3, work_reg_hi_3, 02);
                }

                _mm_storeu_si128(work0, work_reg_lo_0);
                _mm_storeu_si128(work0 + 2, work_reg_hi_0);
                _mm_storeu_si128(work1, work_reg_lo_1);
                _mm_storeu_si128(work1 + 2, work_reg_hi_1);
                _mm_storeu_si128(work2, work_reg_lo_2);
                _mm_storeu_si128(work2 + 2, work_reg_hi_2);
                _mm_storeu_si128(work3, work_reg_lo_3);
                _mm_storeu_si128(work3 + 2, work_reg_hi_3);

                work0++, work1++, work2++, work3++;
            }

            work0 += 2, work1 += 2, work2 += 2, work3 += 2;
            bytes -= 64;
        } while (bytes > 0);

        return;
    }
#endif // LEO_TRY_SSSE3

#endif // LEO_INTERLEAVE_BUTTERFLY4_OPT

    // First layer:
    if (log_m01 == kModulus)
        xor_mem(work[dist], work[0], bytes);
    else
        IFFT_DIT2(work[0], work[dist], log_m01, bytes);

    if (log_m23 == kModulus)
        xor_mem(work[dist * 3], work[dist * 2], bytes);
    else
        IFFT_DIT2(work[dist * 2], work[dist * 3], log_m23, bytes);

    // Second layer:
    if (log_m02 == kModulus)
    {
        xor_mem(work[dist * 2], work[0], bytes);
        xor_mem(work[dist * 3], work[dist], bytes);
    }
    else
    {
        IFFT_DIT2(work[0], work[dist * 2], log_m02, bytes);
        IFFT_DIT2(work[dist], work[dist * 3], log_m02, bytes);
    }
}


// Unrolled IFFT for encoder
static void IFFT_DIT_Encoder(
    const uint64_t bytes,
    const void* const* data,
    const unsigned m_truncated,
    void** work,
    void** xor_result,
    const unsigned m,
    const ffe_t* skewLUT)
{
    // I tried rolling the memcpy/memset into the first layer of the FFT and
    // found that it only yields a 4% performance improvement, which is not
    // worth the extra complexity.
#pragma omp parallel for
    for (int i = 0; i < (int)m_truncated; ++i)
        memcpy(work[i], data[i], bytes);
#pragma omp parallel for
    for (int i = m_truncated; i < (int)m; ++i)
        memset(work[i], 0, bytes);

    // I tried splitting up the first few layers into L3-cache sized blocks but
    // found that it only provides about 5% performance boost, which is not
    // worth the extra complexity.

    // Decimation in time: Unroll 2 layers at a time
    unsigned dist = 1, dist4 = 4;
    for (; dist4 <= m; dist = dist4, dist4 <<= 2)
    {
        // For each set of dist*4 elements:
#pragma omp parallel for
        for (int r = 0; r < (int)m_truncated; r += dist4)
        {
            const unsigned i_end = r + dist;
            const ffe_t log_m01 = skewLUT[i_end];
            const ffe_t log_m02 = skewLUT[i_end + dist];
            const ffe_t log_m23 = skewLUT[i_end + dist * 2];

            // For each set of dist elements:
            for (int i = r; i < (int)i_end; ++i)
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
    }

    // If there is one layer left:
    if (dist < m)
    {
        // Assuming that dist = m / 2
        LEO_DEBUG_ASSERT(dist * 2 == m);

        const ffe_t log_m = skewLUT[dist];

        if (log_m == kModulus)
            VectorXOR_Threads(bytes, dist, work + dist, work);
        else
        {
#pragma omp parallel for
            for (int i = 0; i < (int)dist; ++i)
            {
                IFFT_DIT2(
                    work[i],
                    work[i + dist],
                    log_m,
                    bytes);
            }
        }
    }

    // I tried unrolling this but it does not provide more than 5% performance
    // improvement for 16-bit finite fields, so it's not worth the complexity.
    if (xor_result)
        VectorXOR_Threads(bytes, m, xor_result, work);
}


// Basic no-frills version for decoder
static void IFFT_DIT_Decoder(
    const uint64_t bytes,
    const unsigned m_truncated,
    void** work,
    const unsigned m,
    const ffe_t* skewLUT)
{
    // Decimation in time: Unroll 2 layers at a time
    unsigned dist = 1, dist4 = 4;
    for (; dist4 <= m; dist = dist4, dist4 <<= 2)
    {
        // For each set of dist*4 elements:
#pragma omp parallel for
        for (int r = 0; r < (int)m_truncated; r += dist4)
        {
            const unsigned i_end = r + dist;
            const ffe_t log_m01 = skewLUT[i_end];
            const ffe_t log_m02 = skewLUT[i_end + dist];
            const ffe_t log_m23 = skewLUT[i_end + dist * 2];

            // For each set of dist elements:
            for (int i = r; i < (int)i_end; ++i)
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
    }

    // If there is one layer left:
    if (dist < m)
    {
        // Assuming that dist = m / 2
        LEO_DEBUG_ASSERT(dist * 2 == m);

        const ffe_t log_m = skewLUT[dist];

        if (log_m == kModulus)
            VectorXOR_Threads(bytes, dist, work + dist, work);
        else
        {
#pragma omp parallel for
            for (int i = 0; i < (int)dist; ++i)
            {
                IFFT_DIT2(
                    work[i],
                    work[i + dist],
                    log_m,
                    bytes);
            }
        }
    }
}

/*
    Decimation in time FFT:

    The decimation in time FFT algorithm allows us to unroll 2 layers at a time,
    performing calculations on local registers and faster cache memory.

    Each ^___^ below indicates a butterfly between the associated indices.

    The fft_butterfly(x, y) operation:

        if (log_m != kModulus)
            x[] ^= exp(log(y[]) + log_m)
        y[] ^= x[]

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

// 2-way butterfly
static void FFT_DIT2(
    void * LEO_RESTRICT x, void * LEO_RESTRICT y,
    ffe_t log_m, uint64_t bytes)
{
#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
    {
        LEO_MUL_TABLES_256(0, log_m);

        const LEO_M256 clr_mask = _mm256_set1_epi8(0x0f);

        LEO_M256 * LEO_RESTRICT x32 = reinterpret_cast<LEO_M256 *>(x);
        LEO_M256 * LEO_RESTRICT y32 = reinterpret_cast<LEO_M256 *>(y);

        do
        {
#define LEO_FFTB_256(x_ptr, y_ptr) { \
            LEO_M256 x_lo = _mm256_loadu_si256(x_ptr); \
            LEO_M256 x_hi = _mm256_loadu_si256(x_ptr + 1); \
            LEO_M256 y_lo = _mm256_loadu_si256(y_ptr); \
            LEO_M256 y_hi = _mm256_loadu_si256(y_ptr + 1); \
            LEO_MULADD_256(x_lo, x_hi, y_lo, y_hi, 0); \
            _mm256_storeu_si256(x_ptr, x_lo); \
            _mm256_storeu_si256(x_ptr + 1, x_hi); \
            y_lo = _mm256_xor_si256(y_lo, x_lo); \
            y_hi = _mm256_xor_si256(y_hi, x_hi); \
            _mm256_storeu_si256(y_ptr, y_lo); \
            _mm256_storeu_si256(y_ptr + 1, y_hi); }

            LEO_FFTB_256(x32, y32);
            y32 += 2, x32 += 2;

            bytes -= 64;
        } while (bytes > 0);

        return;
    }
#endif // LEO_TRY_AVX2

#if defined(LEO_TRY_SSSE3)
    if (CpuHasSSSE3)
    {
        LEO_MUL_TABLES_128(0, log_m);

        const LEO_M128 clr_mask = _mm_set1_epi8(0x0f);

        LEO_M128 * LEO_RESTRICT x16 = reinterpret_cast<LEO_M128 *>(x);
        LEO_M128 * LEO_RESTRICT y16 = reinterpret_cast<LEO_M128 *>(y);

        do
        {
#define LEO_FFTB_128(x_ptr, y_ptr) { \
                LEO_M128 x_lo = _mm_loadu_si128(x_ptr); \
                LEO_M128 x_hi = _mm_loadu_si128(x_ptr + 2); \
                LEO_M128 y_lo = _mm_loadu_si128(y_ptr); \
                LEO_M128 y_hi = _mm_loadu_si128(y_ptr + 2); \
                LEO_MULADD_128(x_lo, x_hi, y_lo, y_hi, 0); \
                _mm_storeu_si128(x_ptr, x_lo); \
                _mm_storeu_si128(x_ptr + 2, x_hi); \
                y_lo = _mm_xor_si128(y_lo, x_lo); \
                y_hi = _mm_xor_si128(y_hi, x_hi); \
                _mm_storeu_si128(y_ptr, y_lo); \
                _mm_storeu_si128(y_ptr + 2, y_hi); }

            LEO_FFTB_128(x16 + 1, y16 + 1);
            LEO_FFTB_128(x16, y16);
            x16 += 4, y16 += 4;

            bytes -= 64;
        } while (bytes > 0);

        return;
    }
#endif // LEO_TRY_SSSE3

    // Reference version:
    RefMulAdd(x, y, log_m, bytes);
    xor_mem(y, x, bytes);
}


// 4-way butterfly
static void FFT_DIT4(
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
        LEO_MUL_TABLES_256(01, log_m01);
        LEO_MUL_TABLES_256(23, log_m23);
        LEO_MUL_TABLES_256(02, log_m02);

        const LEO_M256 clr_mask = _mm256_set1_epi8(0x0f);

        LEO_M256 * LEO_RESTRICT work0 = reinterpret_cast<LEO_M256 *>(work[0]);
        LEO_M256 * LEO_RESTRICT work1 = reinterpret_cast<LEO_M256 *>(work[dist]);
        LEO_M256 * LEO_RESTRICT work2 = reinterpret_cast<LEO_M256 *>(work[dist * 2]);
        LEO_M256 * LEO_RESTRICT work3 = reinterpret_cast<LEO_M256 *>(work[dist * 3]);

        do
        {
            LEO_M256 work_reg_lo_0 = _mm256_loadu_si256(work0);
            LEO_M256 work_reg_hi_0 = _mm256_loadu_si256(work0 + 1);
            LEO_M256 work_reg_lo_1 = _mm256_loadu_si256(work1);
            LEO_M256 work_reg_hi_1 = _mm256_loadu_si256(work1 + 1);
            LEO_M256 work_reg_lo_2 = _mm256_loadu_si256(work2);
            LEO_M256 work_reg_hi_2 = _mm256_loadu_si256(work2 + 1);
            LEO_M256 work_reg_lo_3 = _mm256_loadu_si256(work3);
            LEO_M256 work_reg_hi_3 = _mm256_loadu_si256(work3 + 1);

            // First layer:
            if (log_m02 != kModulus)
            {
                LEO_MULADD_256(work_reg_lo_0, work_reg_hi_0, work_reg_lo_2, work_reg_hi_2, 02);
                LEO_MULADD_256(work_reg_lo_1, work_reg_hi_1, work_reg_lo_3, work_reg_hi_3, 02);
            }
            work_reg_lo_2 = _mm256_xor_si256(work_reg_lo_0, work_reg_lo_2);
            work_reg_hi_2 = _mm256_xor_si256(work_reg_hi_0, work_reg_hi_2);
            work_reg_lo_3 = _mm256_xor_si256(work_reg_lo_1, work_reg_lo_3);
            work_reg_hi_3 = _mm256_xor_si256(work_reg_hi_1, work_reg_hi_3);

            // Second layer:
            if (log_m01 != kModulus)
                LEO_MULADD_256(work_reg_lo_0, work_reg_hi_0, work_reg_lo_1, work_reg_hi_1, 01);
            work_reg_lo_1 = _mm256_xor_si256(work_reg_lo_0, work_reg_lo_1);
            work_reg_hi_1 = _mm256_xor_si256(work_reg_hi_0, work_reg_hi_1);

            _mm256_storeu_si256(work0, work_reg_lo_0);
            _mm256_storeu_si256(work0 + 1, work_reg_hi_0);
            _mm256_storeu_si256(work1, work_reg_lo_1);
            _mm256_storeu_si256(work1 + 1, work_reg_hi_1);

            if (log_m23 != kModulus)
                LEO_MULADD_256(work_reg_lo_2, work_reg_hi_2, work_reg_lo_3, work_reg_hi_3, 23);
            work_reg_lo_3 = _mm256_xor_si256(work_reg_lo_2, work_reg_lo_3);
            work_reg_hi_3 = _mm256_xor_si256(work_reg_hi_2, work_reg_hi_3);

            _mm256_storeu_si256(work2, work_reg_lo_2);
            _mm256_storeu_si256(work2 + 1, work_reg_hi_2);
            _mm256_storeu_si256(work3, work_reg_lo_3);
            _mm256_storeu_si256(work3 + 1, work_reg_hi_3);

            work0 += 2, work1 += 2, work2 += 2, work3 += 2;

            bytes -= 64;
        } while (bytes > 0);

        return;
    }

#endif // LEO_TRY_AVX2

#if defined(LEO_TRY_SSSE3)
    if (CpuHasSSSE3)
    {
        LEO_MUL_TABLES_128(01, log_m01);
        LEO_MUL_TABLES_128(23, log_m23);
        LEO_MUL_TABLES_128(02, log_m02);

        const LEO_M128 clr_mask = _mm_set1_epi8(0x0f);

        LEO_M128 * LEO_RESTRICT work0 = reinterpret_cast<LEO_M128 *>(work[0]);
        LEO_M128 * LEO_RESTRICT work1 = reinterpret_cast<LEO_M128 *>(work[dist]);
        LEO_M128 * LEO_RESTRICT work2 = reinterpret_cast<LEO_M128 *>(work[dist * 2]);
        LEO_M128 * LEO_RESTRICT work3 = reinterpret_cast<LEO_M128 *>(work[dist * 3]);

        do
        {
            for (unsigned i = 0; i < 2; ++i)
            {
                LEO_M128 work_reg_lo_0 = _mm_loadu_si128(work0);
                LEO_M128 work_reg_hi_0 = _mm_loadu_si128(work0 + 2);
                LEO_M128 work_reg_lo_1 = _mm_loadu_si128(work1);
                LEO_M128 work_reg_hi_1 = _mm_loadu_si128(work1 + 2);
                LEO_M128 work_reg_lo_2 = _mm_loadu_si128(work2);
                LEO_M128 work_reg_hi_2 = _mm_loadu_si128(work2 + 2);
                LEO_M128 work_reg_lo_3 = _mm_loadu_si128(work3);
                LEO_M128 work_reg_hi_3 = _mm_loadu_si128(work3 + 2);

                // First layer:
                if (log_m02 != kModulus)
                {
                    LEO_MULADD_128(work_reg_lo_0, work_reg_hi_0, work_reg_lo_2, work_reg_hi_2, 02);
                    LEO_MULADD_128(work_reg_lo_1, work_reg_hi_1, work_reg_lo_3, work_reg_hi_3, 02);
                }
                work_reg_lo_2 = _mm_xor_si128(work_reg_lo_0, work_reg_lo_2);
                work_reg_hi_2 = _mm_xor_si128(work_reg_hi_0, work_reg_hi_2);
                work_reg_lo_3 = _mm_xor_si128(work_reg_lo_1, work_reg_lo_3);
                work_reg_hi_3 = _mm_xor_si128(work_reg_hi_1, work_reg_hi_3);

                // Second layer:
                if (log_m01 != kModulus)
                    LEO_MULADD_128(work_reg_lo_0, work_reg_hi_0, work_reg_lo_1, work_reg_hi_1, 01);
                work_reg_lo_1 = _mm_xor_si128(work_reg_lo_0, work_reg_lo_1);
                work_reg_hi_1 = _mm_xor_si128(work_reg_hi_0, work_reg_hi_1);

                _mm_storeu_si128(work0, work_reg_lo_0);
                _mm_storeu_si128(work0 + 2, work_reg_hi_0);
                _mm_storeu_si128(work1, work_reg_lo_1);
                _mm_storeu_si128(work1 + 2, work_reg_hi_1);

                if (log_m23 != kModulus)
                    LEO_MULADD_128(work_reg_lo_2, work_reg_hi_2, work_reg_lo_3, work_reg_hi_3, 23);
                work_reg_lo_3 = _mm_xor_si128(work_reg_lo_2, work_reg_lo_3);
                work_reg_hi_3 = _mm_xor_si128(work_reg_hi_2, work_reg_hi_3);

                _mm_storeu_si128(work2, work_reg_lo_2);
                _mm_storeu_si128(work2 + 2, work_reg_hi_2);
                _mm_storeu_si128(work3, work_reg_lo_3);
                _mm_storeu_si128(work3 + 2, work_reg_hi_3);

                work0++, work1++, work2++, work3++;
            }

            work0 += 2, work1 += 2, work2 += 2, work3 += 2;
            bytes -= 64;
        } while (bytes > 0);

        return;
    }
#endif // LEO_TRY_SSSE3

#endif // LEO_INTERLEAVE_BUTTERFLY4_OPT

    // First layer:
    if (log_m02 == kModulus)
    {
        xor_mem(work[dist * 2], work[0], bytes);
        xor_mem(work[dist * 3], work[dist], bytes);
    }
    else
    {
        FFT_DIT2(work[0], work[dist * 2], log_m02, bytes);
        FFT_DIT2(work[dist], work[dist * 3], log_m02, bytes);
    }

    // Second layer:
    if (log_m01 == kModulus)
        xor_mem(work[dist], work[0], bytes);
    else
        FFT_DIT2(work[0], work[dist], log_m01, bytes);

    if (log_m23 == kModulus)
        xor_mem(work[dist * 3], work[dist * 2], bytes);
    else
        FFT_DIT2(work[dist * 2], work[dist * 3], log_m23, bytes);
}


// In-place FFT for encoder and decoder
static void FFT_DIT(
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
#pragma omp parallel for
        for (int r = 0; r < (int)m_truncated; r += dist4)
        {
            const unsigned i_end = r + dist;
            const ffe_t log_m01 = skewLUT[i_end];
            const ffe_t log_m02 = skewLUT[i_end + dist];
            const ffe_t log_m23 = skewLUT[i_end + dist * 2];

            // For each set of dist elements:
            for (int i = r; i < (int)i_end; ++i)
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
#pragma omp parallel for
        for (int r = 0; r < (int)m_truncated; r += 2)
        {
            const ffe_t log_m = skewLUT[r + 1];

            if (log_m == kModulus)
                xor_mem(work[r + 1], work[r], bytes);
            else
            {
                FFT_DIT2(
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
    const void* const * data,
    void** work)
{
    // work <- IFFT(data, m, m)

    const ffe_t* skewLUT = FFTSkewStorage + m;

    IFFT_DIT_Encoder(
        buffer_bytes,
        data,
        original_count < m ? original_count : m,
        work,
        nullptr, // No xor output
        m,
        skewLUT);

    const unsigned last_count = original_count % m;
    if (m >= original_count)
        goto skip_body;

    // For sets of m data pieces:
    for (unsigned i = m; i + m <= original_count; i += m)
    {
        data += m;
        skewLUT += m;

        // work <- work xor IFFT(data + i, m, m + i)

        IFFT_DIT_Encoder(
            buffer_bytes,
            data, // data source
            m,
            work + m, // temporary workspace
            work, // xor destination
            m,
            skewLUT);
    }

    // Handle final partial set of m pieces:
    if (last_count != 0)
    {
        data += m;
        skewLUT += m;

        // work <- work xor IFFT(data + i, m, m + i)

        IFFT_DIT_Encoder(
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
        FFTSkewStorage);
}


void ReedSolomonEncodeLow(
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned p,
    const void* const* data,
    void* const* recovery,
    void** work)
{
    LEO_DEBUG_ASSERT(p >= 2 && p <= kOrder / 2);
    LEO_DEBUG_ASSERT(original_count <= p);
    LEO_DEBUG_ASSERT(recovery_count <= kOrder - p);

    // Interpolate the systematic prefix.  Coordinates original_count..p-1
    // are shortened and therefore supplied to the transform as known zeroes.
    IFFT_DIT_Encoder(
        buffer_bytes,
        data,
        original_count,
        work,
        nullptr,
        p,
        FFTSkewStorage);

    // Evaluate the coefficient block on each requested parity coset.  The
    // coefficient copy keeps work[0..p-1] reusable across all parity blocks.
    for (unsigned recovery_offset = 0;
        recovery_offset < recovery_count;
        recovery_offset += p)
    {
        unsigned block_count = recovery_count - recovery_offset;
        if (block_count > p)
            block_count = p;

        // The transform supports prefix pruning.  A nullable output is skipped
        // safely; trim a completely unused suffix and avoid empty blocks.
        unsigned requested_count = block_count;
        while (requested_count > 0 && !recovery[recovery_offset + requested_count - 1])
            --requested_count;
        if (requested_count == 0)
            continue;

#pragma omp parallel for
        for (int i = 0; i < (int)p; ++i)
            memcpy(work[p + i], work[i], buffer_bytes);

        FFT_DIT(
            buffer_bytes,
            work + p,
            requested_count,
            p,
            FFTSkewStorage + p + recovery_offset);

#pragma omp parallel for
        for (int i = 0; i < (int)requested_count; ++i)
        {
            void* const output = recovery[recovery_offset + i];
            if (output)
                memcpy(output, work[p + i], buffer_bytes);
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


static void FFT_DIT_ErrorBits(
    const uint64_t bytes,
    void** work,
    const unsigned n_truncated,
    const unsigned n,
    const ffe_t* skewLUT,
    const ErrorBitfield& error_bits)
{
    unsigned mip_level = LastNonzeroBit32(n);

    // Decimation in time: Unroll 2 layers at a time
    unsigned dist4 = n, dist = n >> 2;
    for (; dist != 0; dist4 = dist, dist >>= 2, mip_level -=2)
    {
        // For each set of dist*4 elements:
#pragma omp parallel for
        for (int r = 0; r < (int)n_truncated; r += dist4)
        {
            if (!error_bits.IsNeeded(mip_level, r))
                continue;

            const unsigned i_end = r + dist;
            const ffe_t log_m01 = skewLUT[i_end];
            const ffe_t log_m02 = skewLUT[i_end + dist];
            const ffe_t log_m23 = skewLUT[i_end + dist * 2];

            // For each set of dist elements:
#pragma omp parallel for
            for (int i = r; i < (int)i_end; ++i)
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
#pragma omp parallel for
        for (int r = 0; r < (int)n_truncated; r += 2)
        {
            if (!error_bits.IsNeeded(mip_level, r))
                continue;

            const ffe_t log_m = skewLUT[r + 1];

            if (log_m == kModulus)
                xor_mem(work[r + 1], work[r], bytes);
            else
            {
                FFT_DIT2(
                    work[r],
                    work[r + 1],
                    log_m,
                    bytes);
            }
        }
    }
}

#endif // LEO_ERROR_BITFIELD_OPT


//------------------------------------------------------------------------------
// Reed-Solomon Decode

void PrepareDecodeWalshReference(
    unsigned n,
    const uint8_t* erasures,
    ffe_t* locator_logs)
{
    LEO_DEBUG_ASSERT(n >= 2 && n <= kOrder);

    ffe_t error_locations[kOrder] = {};
    for (unsigned i = 0; i < n; ++i)
        error_locations[i] = erasures[i] ? 1 : 0;

    FWHT(error_locations, kOrder, n);

#pragma omp parallel for
    for (int i = 0; i < (int)kOrder; ++i)
    {
        error_locations[i] = static_cast<ffe_t>(
            ((unsigned)error_locations[i] * (unsigned)LogWalsh[i]) % kModulus);
    }

    FWHT(error_locations, kOrder, kOrder);
    memcpy(locator_logs, error_locations, n * sizeof(ffe_t));
}


static unsigned CountErasures(
    unsigned n,
    const uint8_t* erasures,
    const uint8_t* excluded)
{
    unsigned count = 0;
    for (unsigned i = 0; i < n; ++i)
        if (erasures[i] && (!excluded || !excluded[i]))
            ++count;
    return count;
}


bool IsDirectLocatorPreferred(unsigned n, unsigned erasure_count)
{
    // Whole-locator measurements show that the cache-friendly Walsh kernels
    // cross the scalar direct loop well before their nominal operation counts
    // meet.  Four direct contributions per field coordinate is a conservative
    // portable cutoff; the codec-level calibration table may tighten it later.
    const uint64_t direct_work = (uint64_t)n * erasure_count;
    const uint64_t walsh_work = (uint64_t)kOrder * 4;
    return direct_work <= walsh_work;
}


static void AddDirectLocatorContributions(
    unsigned n,
    const uint8_t* erasures,
    const uint8_t* excluded,
    ffe_t* locator_logs)
{
    // omega_i + omega_j = omega_(i xor j) in the Cantor coordinate order.
    // At an erased coordinate the self factor is omitted, yielding the
    // locator derivative needed by the recovery formula.  This is identical
    // to the LogWalsh[0] = 0 convention in the reference convolution.
    for (unsigned erased = 0; erased < n; ++erased)
    {
        if (!erasures[erased] || (excluded && excluded[erased]))
            continue;
        for (unsigned i = 0; i < n; ++i)
            if (i != erased)
                locator_logs[i] = AddMod(locator_logs[i], LogLUT[i ^ erased]);
    }
}


void PrepareDecode(
    unsigned n,
    const uint8_t* erasures,
    ffe_t* locator_logs)
{
    LEO_DEBUG_ASSERT(n >= 2 && n <= kOrder);

    const unsigned erasure_count = CountErasures(n, erasures, nullptr);
    if (!IsDirectLocatorPreferred(n, erasure_count))
    {
        PrepareDecodeWalshReference(n, erasures, locator_logs);
        return;
    }

    memset(locator_logs, 0, n * sizeof(ffe_t));
    AddDirectLocatorContributions(n, erasures, nullptr, locator_logs);
}


void PrepareDecodeWithPermanent(
    unsigned n,
    const uint8_t* erasures,
    const uint8_t* permanent_erasures,
    const ffe_t* permanent_locator_logs,
    ffe_t* locator_logs)
{
    LEO_DEBUG_ASSERT(n >= 2 && n <= kOrder);

    const unsigned dynamic_count = CountErasures(
        n, erasures, permanent_erasures);
    if (!IsDirectLocatorPreferred(n, dynamic_count))
    {
        PrepareDecode(n, erasures, locator_logs);
        return;
    }

    memcpy(locator_logs, permanent_locator_logs, n * sizeof(ffe_t));
    AddDirectLocatorContributions(
        n, erasures, permanent_erasures, locator_logs);
}


// mul_mem() uses restrict-qualified input and output pointers.  Keep its
// contract intact while revealing a result in its existing coordinate slot.
static void mul_mem_inplace(void* data, ffe_t log_m, uint64_t bytes)
{
    alignas(64) uint8_t source[64];
    uint8_t* output = reinterpret_cast<uint8_t*>(data);

    for (uint64_t offset = 0; offset < bytes; offset += sizeof(source))
    {
        memcpy(source, output + offset, sizeof(source));
        mul_mem(output + offset, source, log_m, sizeof(source));
    }
}


static void AddFormalDerivative(
    uint64_t buffer_bytes,
    unsigned n,
    void** work)
{
    // In the normalized LCH basis this triangular XOR adds f' to f.  The
    // decoder evaluates the result only at locator roots, where f is zero, so
    // f + f' and the formal derivative have the same requested values.
    for (unsigned i = 1; i < n; ++i)
    {
        const unsigned width = ((i ^ (i - 1)) + 1) >> 1;

        if (width < 8)
            VectorXOR(buffer_bytes, width, work + i - width, work + i);
        else
            VectorXOR_Threads(buffer_bytes, width, work + i - width, work + i);
    }
}


void ReedSolomonDecodePrepared(
    uint64_t buffer_bytes,
    unsigned n,
    const void* const* coordinate_data,
    const uint8_t* requested_outputs,
    const ffe_t* locator_logs,
    void** work)
{
    LEO_DEBUG_ASSERT(n >= 2 && n <= kOrder);

    unsigned input_count = n;
    while (input_count > 0 && !coordinate_data[input_count - 1])
        --input_count;

#ifdef LEO_ERROR_BITFIELD_OPT
    ErrorBitfield error_bits;
    for (unsigned i = 0; i < n; ++i)
        if (requested_outputs[i])
            error_bits.Set(i);
    error_bits.Prepare();
#endif // LEO_ERROR_BITFIELD_OPT

#pragma omp parallel for
    for (int i = 0; i < (int)n; ++i)
    {
        if (coordinate_data[i])
            mul_mem(work[i], coordinate_data[i], locator_logs[i], buffer_bytes);
        else
            memset(work[i], 0, buffer_bytes);
    }

    IFFT_DIT_Decoder(buffer_bytes, input_count, work, n, FFTSkewStorage);

    AddFormalDerivative(buffer_bytes, n, work);

#ifdef LEO_ERROR_BITFIELD_OPT
    FFT_DIT_ErrorBits(buffer_bytes, work, n, n, FFTSkewStorage, error_bits);
#else
    FFT_DIT(buffer_bytes, work, n, n, FFTSkewStorage);
#endif // LEO_ERROR_BITFIELD_OPT

#pragma omp parallel for
    for (int i = 0; i < (int)n; ++i)
    {
        if (requested_outputs[i])
            mul_mem_inplace(work[i], kModulus - locator_logs[i], buffer_bytes);
    }
}


static ffe_t SubspaceDerivativeLog(unsigned size)
{
    ffe_t result = 0;
    for (unsigned i = 1; i < size; ++i)
        result = AddMod(result, LogLUT[i]);
    return result;
}


static ffe_t SubspaceAtLog(unsigned size, unsigned shift)
{
    ffe_t result = 0;
    for (unsigned i = 0; i < size; ++i)
        result = AddMod(result, LogLUT[shift ^ i]);
    return result;
}


static ffe_t LchNormalizerLog(unsigned index)
{
    ffe_t result = 0;
    for (unsigned bit = 0; bit < kBits; ++bit)
    {
        if ((index & (1UL << bit)) == 0)
            continue;
        result = AddMod(result, LchBasisNormalizerLog[bit]);
    }
    return result;
}


void PrepareLowDecode(
    unsigned n,
    unsigned p,
    ffe_t* block_factors)
{
    LEO_DEBUG_ASSERT(p >= 2 && p <= n && n <= kOrder);

    const ffe_t numerator_log = SubspaceDerivativeLog(p);

    const unsigned block_count = n / p;
    for (unsigned block = 1; block < block_count; ++block)
    {
        const unsigned shift = block * p;
        const ffe_t denominator_log = SubspaceAtLog(p, shift);
        block_factors[block - 1] = SubMod(numerator_log, denominator_log);
    }
}


void PrepareHighDecode(
    unsigned n,
    unsigned t,
    ffe_t* output_factors)
{
    LEO_DEBUG_ASSERT(t >= 2 && t < n && n <= kOrder);

    // c_n = product(V_n \ {0}).  R10's full-field form silently has c_m=1;
    // retaining this numerator is required for an active parent in GF16.
    const ffe_t active_derivative_log = SubspaceDerivativeLog(n);

    // p_(N-T) from the normalized LCH basis.
    const ffe_t normalization_log = LchNormalizerLog(n - t);

    const unsigned block_count = n / t;
#pragma omp parallel for
    for (int block = 1; block < (int)block_count; ++block)
    {
        const unsigned coordinate = (unsigned)block * t;
        const ffe_t subspace_log = SubspaceAtLog(t, coordinate);
        const ffe_t factor = SubMod(
            active_derivative_log,
            AddMod(normalization_log, subspace_log));
        for (unsigned i = 0; i < t; ++i)
            output_factors[coordinate + i] = factor;
    }
}


#if defined(LEO2_ENABLE_TEST_HOOKS)

static ffe_t NonzeroElementFromLog(ffe_t value_log)
{
    return ExpLUT[value_log];
}


ffe_t TestOnlyFFTMultiplier(unsigned logical_index)
{
    LEO_DEBUG_ASSERT(logical_index < kModulus);
    const ffe_t multiplier_log = FFTSkew[logical_index];
    return multiplier_log == kModulus ? 0 : ExpLUT[multiplier_log];
}


ffe_t TestOnlySubspaceDerivative(unsigned size)
{
    LEO_DEBUG_ASSERT(size >= 1 && size <= kOrder);
    return NonzeroElementFromLog(SubspaceDerivativeLog(size));
}


ffe_t TestOnlySubspaceAt(unsigned size, unsigned shift)
{
    LEO_DEBUG_ASSERT(size >= 1 && size <= kOrder);
    LEO_DEBUG_ASSERT((shift & (size - 1)) == 0 && shift + size <= kOrder);
    if (shift < size)
        return 0;
    return NonzeroElementFromLog(SubspaceAtLog(size, shift));
}


ffe_t TestOnlyLchNormalizer(unsigned index)
{
    LEO_DEBUG_ASSERT(index < kOrder);
    return NonzeroElementFromLog(LchNormalizerLog(index));
}


void TestOnlyLchForward(
    uint64_t buffer_bytes,
    unsigned size,
    unsigned shift,
    unsigned requested_output_count,
    void** work)
{
    LEO_DEBUG_ASSERT(size >= 1 && size <= kOrder);
    LEO_DEBUG_ASSERT((size & (size - 1)) == 0);
    LEO_DEBUG_ASSERT((shift & (size - 1)) == 0 && shift + size <= kOrder);
    LEO_DEBUG_ASSERT(requested_output_count <= size);
    FFT_DIT(buffer_bytes, work, requested_output_count, size,
        FFTSkewStorage + shift);
}


void TestOnlyLchInverse(
    uint64_t buffer_bytes,
    unsigned size,
    unsigned shift,
    unsigned known_input_count,
    void** work)
{
    LEO_DEBUG_ASSERT(size >= 1 && size <= kOrder);
    LEO_DEBUG_ASSERT((size & (size - 1)) == 0);
    LEO_DEBUG_ASSERT((shift & (size - 1)) == 0 && shift + size <= kOrder);
    LEO_DEBUG_ASSERT(known_input_count <= size);
    IFFT_DIT_Decoder(buffer_bytes, known_input_count, work, size,
        FFTSkewStorage + shift);
}


void TestOnlyAddFormalDerivative(
    uint64_t buffer_bytes,
    unsigned size,
    void** work)
{
    AddFormalDerivative(buffer_bytes, size, work);
}

#endif // LEO2_ENABLE_TEST_HOOKS


void ReedSolomonDecodeLowPrepared(
    uint64_t buffer_bytes,
    unsigned n,
    unsigned p,
    const void* const* coordinate_data,
    const uint8_t* requested_outputs,
    const ffe_t* locator_logs,
    const ffe_t* block_factors,
    void** work)
{
    LEO_DEBUG_ASSERT(p >= 2 && p <= n && n <= kOrder);

#ifdef LEO_ERROR_BITFIELD_OPT
    ErrorBitfield error_bits;
    for (unsigned i = 0; i < p; ++i)
        if (requested_outputs[i])
            error_bits.Set(i);
    error_bits.Prepare();
#endif

#pragma omp parallel for
    for (int i = 0; i < (int)n; ++i)
    {
        if (coordinate_data[i])
            mul_mem(work[i], coordinate_data[i], locator_logs[i], buffer_bytes);
        else
            memset(work[i], 0, buffer_bytes);
    }

    const unsigned block_count = n / p;
    for (unsigned block = 0; block < block_count; ++block)
    {
        const unsigned offset = block * p;
        unsigned input_count = p;
        while (input_count > 0 && !coordinate_data[offset + input_count - 1])
            --input_count;
        IFFT_DIT_Decoder(
            buffer_bytes, input_count, work + offset, p,
            FFTSkewStorage + offset);
    }

    // R10 Corollary 1 derivative term.  The omitted scalar-sum contribution
    // is f-hat^(0), which vanishes at every requested erasure because
    // f(e) * Lambda(e) = 0; leaving slot zero unchanged is therefore exact.
    AddFormalDerivative(buffer_bytes, p, work);

    // Weighted block reduction by c_k / s_k(omega_(block*P)).
    for (unsigned block = 1; block < block_count; ++block)
    {
        const unsigned offset = block * p;
#pragma omp parallel for
        for (int i = 0; i < (int)p; ++i)
        {
            mul_mem_inplace(work[offset + i], block_factors[block - 1], buffer_bytes);
            xor_mem(work[i], work[offset + i], buffer_bytes);
        }
    }

#ifdef LEO_ERROR_BITFIELD_OPT
    FFT_DIT_ErrorBits(buffer_bytes, work, p, p, FFTSkewStorage, error_bits);
#else
    FFT_DIT(buffer_bytes, work, p, p, FFTSkewStorage);
#endif

#pragma omp parallel for
    for (int i = 0; i < (int)p; ++i)
        if (requested_outputs[i])
            mul_mem_inplace(work[i], kModulus - locator_logs[i], buffer_bytes);
}


void ReedSolomonDecodeHighPrepared(
    uint64_t buffer_bytes,
    unsigned n,
    unsigned t,
    const void* const* coordinate_data,
    const uint8_t* requested_outputs,
    const ffe_t* locator_logs,
    const ffe_t* output_factors,
    void** work)
{
    LEO_DEBUG_ASSERT(t >= 2 && t < n && n <= kOrder);

#pragma omp parallel for
    for (int i = 0; i < (int)n; ++i)
    {
        if (coordinate_data[i])
            memcpy(work[i], coordinate_data[i], buffer_bytes);
        else
            memset(work[i], 0, buffer_bytes);
    }

    const unsigned block_count = n / t;
    for (unsigned block = 0; block < block_count; ++block)
    {
        const unsigned offset = block * t;
        unsigned input_count = t;
        while (input_count > 0 && !coordinate_data[offset + input_count - 1])
            --input_count;
        IFFT_DIT_Decoder(
            buffer_bytes, input_count, work + offset, t,
            FFTSkewStorage + offset);
        if (block != 0)
        {
            if (t < 8)
                VectorXOR(buffer_bytes, t, work, work + offset);
            else
                VectorXOR_Threads(buffer_bytes, t, work, work + offset);
        }
    }

    // h on V_t, then z = h * Lambda on V_t.
    FFT_DIT(buffer_bytes, work, t, t, FFTSkewStorage);
#pragma omp parallel for
    for (int i = 0; i < (int)t; ++i)
    {
        if (coordinate_data[i])
            mul_mem_inplace(work[i], locator_logs[i], buffer_bytes);
        else
            memset(work[i], 0, buffer_bytes);
    }
    IFFT_DIT_Decoder(buffer_bytes, t, work, t, FFTSkewStorage);

    // Evaluate z only on message blocks that contain a requested original.
    for (unsigned block = 1; block < block_count; ++block)
    {
        const unsigned offset = block * t;
        unsigned requested_count = t;
        while (requested_count > 0 && !requested_outputs[offset + requested_count - 1])
            --requested_count;
        if (requested_count == 0)
            continue;
#pragma omp parallel for
        for (int i = 0; i < (int)t; ++i)
            memcpy(work[offset + i], work[i], buffer_bytes);
        FFT_DIT(
            buffer_bytes, work + offset, requested_count, t,
            FFTSkewStorage + offset);
#pragma omp parallel for
        for (int i = 0; i < (int)requested_count; ++i)
        {
            const unsigned coordinate = offset + i;
            if (requested_outputs[coordinate])
            {
                const ffe_t reveal_log = SubMod(
                    output_factors[coordinate], locator_logs[coordinate]);
                mul_mem_inplace(work[coordinate], reveal_log, buffer_bytes);
            }
        }
    }
}

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
    ErrorBitfield error_bits;
#endif // LEO_ERROR_BITFIELD_OPT

    ffe_t error_locations[kOrder] = {};
    for (unsigned i = 0; i < recovery_count; ++i)
        if (!recovery[i])
            error_locations[i] = 1;
    for (unsigned i = recovery_count; i < m; ++i)
        error_locations[i] = 1;
    for (unsigned i = 0; i < original_count; ++i)
    {
        if (!original[i])
        {
            error_locations[i + m] = 1;
#ifdef LEO_ERROR_BITFIELD_OPT
            error_bits.Set(i + m);
#endif // LEO_ERROR_BITFIELD_OPT
        }
    }

#ifdef LEO_ERROR_BITFIELD_OPT
    error_bits.Prepare();
#endif // LEO_ERROR_BITFIELD_OPT

    // Evaluate error locator polynomial

    FWHT(error_locations, kOrder, m + original_count);

#pragma omp parallel for
    for (int i = 0; i < (int)kOrder; ++i)
        error_locations[i] = ((unsigned)error_locations[i] * (unsigned)LogWalsh[i]) % kModulus;

    FWHT(error_locations, kOrder, kOrder);

    // work <- recovery data

#pragma omp parallel for
    for (int i = 0; i < (int)recovery_count; ++i)
    {
        if (recovery[i])
            mul_mem(work[i], recovery[i], error_locations[i], buffer_bytes);
        else
            memset(work[i], 0, buffer_bytes);
    }
#pragma omp parallel for
    for (int i = recovery_count; i < (int)m; ++i)
        memset(work[i], 0, buffer_bytes);

    // work <- original data

#pragma omp parallel for
    for (int i = 0; i < (int)original_count; ++i)
    {
        if (original[i])
            mul_mem(work[m + i], original[i], error_locations[m + i], buffer_bytes);
        else
            memset(work[m + i], 0, buffer_bytes);
    }
#pragma omp parallel for
    for (int i = m + original_count; i < (int)n; ++i)
        memset(work[i], 0, buffer_bytes);

    // work <- IFFT(work, n, 0)

    IFFT_DIT_Decoder(
        buffer_bytes,
        m + original_count,
        work,
        n,
        FFTSkewStorage);

    // work <- FormalDerivative(work, n)

    for (unsigned i = 1; i < n; ++i)
    {
        const unsigned width = ((i ^ (i - 1)) + 1) >> 1;

        if (width < 8)
        {
            VectorXOR(
                buffer_bytes,
                width,
                work + i - width,
                work + i);
        }
        else
        {
            VectorXOR_Threads(
                buffer_bytes,
                width,
                work + i - width,
                work + i);
        }
    }

    // work <- FFT(work, n, 0) truncated to m + original_count

    const unsigned output_count = m + original_count;

#ifdef LEO_ERROR_BITFIELD_OPT
    FFT_DIT_ErrorBits(
        buffer_bytes, work, output_count, n, FFTSkewStorage, error_bits);
#else
    FFT_DIT(buffer_bytes, work, output_count, n, FFTSkewStorage);
#endif

    // Reveal erasures

    for (unsigned i = 0; i < original_count; ++i)
        if (!original[i])
            mul_mem(work[i], work[i + m], kModulus - error_locations[i + m], buffer_bytes);
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


}} // namespace leopard::ff16

#endif // LEO_HAS_FF16
