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
#include "Leopard2Backend.h"

#ifdef LEO_HAS_FF16

#include <string.h>

#if defined(LEO2_ENABLE_TEST_HOOKS)
#include <atomic>
#endif

#ifdef _MSC_VER
    #pragma warning(disable: 4752) // found Intel(R) Advanced Vector Extensions; consider using /arch:AVX
#endif

namespace leopard { namespace ff16 {

#if defined(LEO2_ENABLE_TEST_HOOKS)
static std::atomic<uint64_t> TestIFFTDIT4Calls(0);
static std::atomic<uint64_t> TestFFTDIT4Calls(0);
static std::atomic<uint64_t> TestIFFTDIT4FusedCalls(0);
static std::atomic<uint64_t> TestIFFTDIT4SplitCalls(0);
static std::atomic<uint64_t> TestFFTDIT4FusedCalls(0);
static std::atomic<uint64_t> TestFFTDIT4SplitCalls(0);
static std::atomic<uint64_t> TestLowFFTButterfly2OutCalls(0);
static std::atomic<uint64_t> TestLowFFTButterfly4OutCalls(0);
static std::atomic<uint64_t> TestLowDirectOutputBlocks(0);
static std::atomic<uint64_t> TestLowDirectFFTButterfly2OutCalls(0);
static std::atomic<uint64_t> TestLowDirectFFTButterfly4OutCalls(0);
static std::atomic<uint64_t> TestHighIFFTButterfly4OutCalls(0);
static std::atomic<uint64_t> TestHighInputCopyShards(0);
static std::atomic<uint64_t> TestHighOutputBlocks(0);
static std::atomic<uint64_t> TestHighFFTButterfly2OutCalls(0);
static std::atomic<uint64_t> TestHighFFTButterfly4OutCalls(0);
static std::atomic<uint64_t> TestHighCompatibilityCopyFallbacks(0);
static std::atomic<uint64_t> TestHighPrunedOutputBlocks(0);
static std::atomic<uint64_t> TestHighMatureOutputBlocks(0);
static std::atomic<bool> TestForceHighDecodeCopyFallback(false);
static std::atomic<uint64_t> TestHighSyndromeAccumulatedBlocks(0);
static std::atomic<uint64_t> TestHighSyndromeMaterializedBlocks(0);
static std::atomic<uint64_t> TestHighSyndromePrunedAccumulatedBlocks(0);
static std::atomic<uint64_t> TestHighSyndromePrunedFallbackBlocks(0);
static std::atomic<uint64_t> TestHighSyndromeBlockZeroIFFTElisions(0);
static std::atomic<uint64_t> TestHighSyndromeForwardTransforms(0);
static std::atomic<uint64_t> TestHighSyndromeForwardTransformElisions(0);
static std::atomic<uint64_t> TestHighSyndromeBlockZeroXorShards(0);
static std::atomic<uint64_t> TestHighReceiveCopyShards(0);
static std::atomic<uint64_t> TestHighReceiveZeroShards(0);
static std::atomic<uint64_t> TestSparseExactBlocks(0);
static std::atomic<uint64_t> TestSparsePrefixButterflies(0);
static std::atomic<uint64_t> TestSparseRetainedButterflies(0);
static std::atomic<uint64_t> TestSparseRequestedOutputCopies(0);
#endif

static uint16_t PrunedMultiplierLog(
    const void* context,
    uint32_t storage_index);

#if defined(LEO2_ENABLE_TEST_HOOKS)
static void TestRecordSparseEncodeBlock(
    unsigned transform_size,
    unsigned requested_prefix,
    unsigned requested_outputs,
    const leopard2_internal::SparseForwardPlanBatchView* plans,
    unsigned block)
{
    if (!plans || block >= plans->block_count)
        return;
    TestSparsePrefixButterflies.fetch_add(
        leopard2_internal::PrefixForwardButterflyCount(
            transform_size, requested_prefix),
        std::memory_order_relaxed);
    TestSparseRequestedOutputCopies.fetch_add(
        requested_outputs, std::memory_order_relaxed);
    TestSparseExactBlocks.fetch_add(1, std::memory_order_relaxed);
    TestSparseRetainedButterflies.fetch_add(
        leopard2_internal::CountSparseForwardRetainedButterflies(
            transform_size,
            plans->operation_masks +
                static_cast<size_t>(block) * plans->operation_stride,
            plans->operation_stride),
        std::memory_order_relaxed);
}
#endif

static inline bool UseFusedButterfly4(
    const backend::Ops& ops,
    uint64_t bytes)
{
    // Production retains the existing one-tile operation across backend
    // tables.  Final-source x86 evidence also qualifies two tiles for AVX2,
    // while the SSSE3 two-tile kernel regresses compact public tails that stage
    // to 128 bytes.  A positive AVX2 check keeps the unmeasured native-ARM
    // facade on the conservative schedule even though that facade currently
    // carries a scalar Ops table.
    // Larger working sets retain the two-layer schedule, which has lower
    // table/register pressure despite performing more backend calls.
    return bytes == 64 ||
        (bytes == 128 &&
         (ops.kind == LEO2_BACKEND_AVX2 ||
          ops.kind == LEO2_BACKEND_AVX512));
}


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
template<bool EnableParallel>
static void FWHT(ffe_t* data, const unsigned m, const unsigned m_truncated)
{
    // Decimation in time: Unroll 2 layers at a time
    unsigned dist = 1, dist4 = 4;
    for (; dist4 <= m; dist = dist4, dist4 <<= 2)
    {
        // For each set of dist*4 elements:
        LEO_OPENMP_PARALLEL_FOR_IF(EnableParallel)
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
        LEO_OPENMP_PARALLEL_FOR_IF(EnableParallel)
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

#if defined(LEO_TRY_SSSE3)
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
#endif // LEO_TRY_SSSE3

static bool InitializeMultiplyTables()
{
#if !defined(LEO_TRY_SSSE3)
    // Portable CMake builds use the selected backend's process-lifetime FF16
    // table.  The former 8-MiB legacy scalar table was write-only.
    return true;
#else
    if (!CpuHasSSSE3)
        return true;

    void* table = nullptr;
#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
        table = SIMDSafeAllocate(sizeof(Multiply256LUT_t) * kOrder);
    else
#endif // LEO_TRY_AVX2
        table = SIMDSafeAllocate(sizeof(Multiply128LUT_t) * kOrder);
    if (!table)
        return false;
#if defined(LEO_TRY_AVX2)
    if (CpuHasAVX2)
        Multiply256LUT = reinterpret_cast<const Multiply256LUT_t*>(table);
    else
#endif // LEO_TRY_AVX2
        Multiply128LUT = reinterpret_cast<const Multiply128LUT_t*>(table);

    // For each value we could multiply by:
    // Table construction is one-time setup.  Keep it serial so leo_init() and
    // leo2_context_create() do not instantiate a persistent OpenMP worker team
    // before the application asks Leopard to execute byte-heavy work.
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
    return true;
#endif // LEO_TRY_SSSE3
}


static void mul_mem(
    const backend::Ops& ops,
    void * LEO_RESTRICT x, const void * LEO_RESTRICT y,
    ffe_t log_m, uint64_t bytes)
{
#if defined(LEO_TRY_AVX2)
    if (&ops == &backend::GetDefaultOps() && CpuHasAVX2)
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
    if (&ops == &backend::GetDefaultOps() && CpuHasSSSE3)
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

    ops.ff16_multiply(x, y, log_m, bytes);
}

// x[] ^= y[] * log_m for complete 64-byte ALTMAP tiles.
static void muladd_mem(
    const backend::Ops& ops,
    void * LEO_RESTRICT x, const void * LEO_RESTRICT y,
    ffe_t log_m, uint64_t bytes)
{
#if defined(LEO_TRY_AVX2)
    if (&ops == &backend::GetDefaultOps() && CpuHasAVX2)
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
    if (&ops == &backend::GetDefaultOps() && CpuHasSSSE3)
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

    ops.ff16_multiply_add(x, y, log_m, bytes);
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


ffe_t MultiplyLogElement(ffe_t value, ffe_t multiplier_log)
{
    return MultiplyLog(value, multiplier_log);
}


void MultiplyBytes(
    const backend::Ops& ops,
    void* destination,
    const void* source,
    ffe_t multiplier_log,
    uint64_t byte_count)
{
    const uint64_t complete = byte_count & ~static_cast<uint64_t>(63);
    if (complete != 0)
        mul_mem(ops, destination, source, multiplier_log, complete);

    const uint64_t residual = byte_count - complete;
    LEO_DEBUG_ASSERT((residual & 1) == 0);
    if (residual != 0)
        ops.ff16_multiply(
            static_cast<uint8_t*>(destination) + complete,
            static_cast<const uint8_t*>(source) + complete,
            multiplier_log, residual);
}


void MultiplyBytes(
    void* destination,
    const void* source,
    ffe_t multiplier_log,
    uint64_t byte_count)
{
    MultiplyBytes(backend::GetDefaultOps(), destination, source,
        multiplier_log, byte_count);
}


void MultiplyAddBytes(
    const backend::Ops& ops,
    void* destination,
    const void* source,
    ffe_t multiplier_log,
    uint64_t byte_count)
{
    const uint64_t complete = byte_count & ~static_cast<uint64_t>(63);
    if (complete != 0)
        muladd_mem(ops, destination, source, multiplier_log, complete);

    const uint64_t residual = byte_count - complete;
    LEO_DEBUG_ASSERT((residual & 1) == 0);
    if (residual != 0)
        ops.ff16_multiply_add(
            static_cast<uint8_t*>(destination) + complete,
            static_cast<const uint8_t*>(source) + complete,
            multiplier_log, residual);
}


void MultiplyAddBytes(
    void* destination,
    const void* source,
    ffe_t multiplier_log,
    uint64_t byte_count)
{
    MultiplyAddBytes(backend::GetDefaultOps(), destination, source,
        multiplier_log, byte_count);
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

// Walsh transforms of LogLUT restricted to each proper active additive
// subspace.  A table of length n starts at n - 2, because the preceding
// power-of-two lengths sum to n - 2.  The full-field table is LogWalsh.
static ffe_t ActiveLogWalsh[kOrder - 2];


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

    FWHT<false>(LogWalsh, kOrder, kOrder);

    // On V_n, XOR convolution is diagonalized by an n-point Walsh transform.
    // Unlike the full-field case, applying that transform twice multiplies by
    // n modulo kModulus rather than by kOrder == 1 (mod kModulus).  Scale the
    // fixed transformed kernel by n^-1 = kOrder / n (mod kModulus) once here.
    for (unsigned n = 2; n < kOrder; n <<= 1)
    {
        ffe_t* active = ActiveLogWalsh + n - 2;
        for (unsigned i = 0; i < n; ++i)
            active[i] = LogLUT[i];
        active[0] = 0;
        FWHT<false>(active, n, n);

        const unsigned inverse_n = kOrder / n;
        for (unsigned i = 0; i < n; ++i)
        {
            active[i] = static_cast<ffe_t>(
                (static_cast<unsigned>(active[i]) * inverse_n) % kModulus);
        }
    }
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
    const backend::Ops& ops,
    void * LEO_RESTRICT x, void * LEO_RESTRICT y,
    ffe_t log_m, uint64_t bytes)
{
#if defined(LEO_TRY_AVX2)
    if (&ops == &backend::GetDefaultOps() && CpuHasAVX2)
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
    if (&ops == &backend::GetDefaultOps() && CpuHasSSSE3)
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

    ops.ff16_ifft_butterfly2(x, y, log_m, bytes);
}


// 4-way butterfly
static void IFFT_DIT4(
    const backend::Ops& ops,
    uint64_t bytes,
    void** work,
    unsigned dist,
    const ffe_t log_m01,
    const ffe_t log_m23,
    const ffe_t log_m02)
{
    const bool use_fused_four_way = UseFusedButterfly4(ops, bytes);
#if defined(LEO2_ENABLE_TEST_HOOKS)
    TestIFFTDIT4Calls.fetch_add(1, std::memory_order_relaxed);
    (use_fused_four_way ? TestIFFTDIT4FusedCalls : TestIFFTDIT4SplitCalls)
        .fetch_add(1, std::memory_order_relaxed);
#endif
#ifdef LEO_INTERLEAVE_BUTTERFLY4_OPT

#if defined(LEO_TRY_AVX2)

    if (use_fused_four_way &&
        &ops == &backend::GetDefaultOps() && CpuHasAVX2)
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
    if (use_fused_four_way &&
        &ops == &backend::GetDefaultOps() && CpuHasSSSE3)
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

    if (use_fused_four_way)
    {
        ops.ff16_ifft_butterfly4(
            work[0], work[dist], work[dist * 2], work[dist * 3],
            log_m01, log_m23, log_m02, bytes);
        return;
    }

    // First layer:
    if (log_m01 == kModulus)
        xor_mem(ops, work[dist], work[0], bytes);
    else
        IFFT_DIT2(ops, work[0], work[dist], log_m01, bytes);

    if (log_m23 == kModulus)
        xor_mem(ops, work[dist * 3], work[dist * 2], bytes);
    else
        IFFT_DIT2(ops, work[dist * 2], work[dist * 3], log_m23, bytes);

    // Second layer:
    if (log_m02 == kModulus)
    {
        xor_mem(ops, work[dist * 2], work[0], bytes);
        xor_mem(ops, work[dist * 3], work[dist], bytes);
    }
    else
    {
        IFFT_DIT2(ops, work[0], work[dist * 2], log_m02, bytes);
        IFFT_DIT2(ops, work[dist], work[dist * 3], log_m02, bytes);
    }
}

static void IFFT_DIT4_Range(
    const backend::Ops& ops,
    uint64_t bytes,
    void** work,
    unsigned dist,
    const ffe_t log_m01,
    const ffe_t log_m23,
    const ffe_t log_m02)
{
    if (dist == 1)
    {
        IFFT_DIT4(ops, bytes, work, dist,
            log_m01, log_m23, log_m02);
        return;
    }
#if defined(LEO_TRY_AVX2) || defined(LEO_TRY_SSSE3)
    // Preserve the legacy default-context in-field SIMD route on targets
    // where the field translation unit intentionally owns it.  Portable x86
    // production disables these macros and dispatches into an isolated ISA
    // backend below.
    if (&ops == &backend::GetDefaultOps())
    {
        for (unsigned i = 0; i < dist; ++i)
            IFFT_DIT4(ops, bytes, work + i, dist,
                log_m01, log_m23, log_m02);
        return;
    }
#endif
    const bool use_fused_four_way = UseFusedButterfly4(ops, bytes);
#if defined(LEO2_ENABLE_TEST_HOOKS)
    TestIFFTDIT4Calls.fetch_add(dist, std::memory_order_relaxed);
    (use_fused_four_way ? TestIFFTDIT4FusedCalls : TestIFFTDIT4SplitCalls)
        .fetch_add(dist, std::memory_order_relaxed);
#endif
    ops.ff16_ifft_butterfly4_range(
        work, dist, log_m01, log_m23, log_m02,
        bytes, use_fused_four_way);
}


// {x_out, y_out} ^= IFFT_DIT2({x_in, y_in}).  The transform workspace is
// dead after the encoder's final stage, so retaining the transformed values
// there only creates a store/reload pass before accumulation.
static void IFFT_DIT2_xor(
    const backend::Ops& ops,
    const void* x_in,
    const void* y_in,
    void* x_out,
    void* y_out,
    const ffe_t log_m,
    const uint64_t bytes)
{
    if (log_m == kModulus)
    {
        // Zero skew suppresses multiplication but retains y ^= x.
        ops.xor_memory(x_out, x_in, bytes);
        ops.xor_memory_2to1(y_out, x_in, y_in, bytes);
        return;
    }
    ops.ff16_ifft_butterfly2_xor(
        x_in, y_in, x_out, y_out, log_m, bytes);
}


// Accumulating form of the established two-layer GF16 schedule.  This is
// selected only where UseFusedButterfly4() already rejects the four-way
// kernel.  The first layer remains in place; the second layer writes directly
// into the accumulator and avoids materializing its dead workspace result.
static void IFFT_DIT4_xor_Split_Range(
    const backend::Ops& ops,
    const uint64_t bytes,
    void** work,
    void** xor_output,
    const unsigned dist,
    const ffe_t log_m01,
    const ffe_t log_m23,
    const ffe_t log_m02)
{
#if defined(LEO2_ENABLE_TEST_HOOKS)
    TestIFFTDIT4Calls.fetch_add(dist, std::memory_order_relaxed);
    TestIFFTDIT4SplitCalls.fetch_add(dist, std::memory_order_relaxed);
#endif
    for (unsigned i = 0; i < dist; ++i)
    {
        if (log_m01 == kModulus)
            xor_mem(ops, work[i + dist], work[i], bytes);
        else
            IFFT_DIT2(ops, work[i], work[i + dist], log_m01, bytes);

        if (log_m23 == kModulus)
            xor_mem(ops, work[i + dist * 3U], work[i + dist * 2U], bytes);
        else
            IFFT_DIT2(ops,
                work[i + dist * 2U], work[i + dist * 3U],
                log_m23, bytes);

        IFFT_DIT2_xor(ops,
            work[i], work[i + dist * 2U],
            xor_output[i], xor_output[i + dist * 2U],
            log_m02, bytes);
        IFFT_DIT2_xor(ops,
            work[i + dist], work[i + dist * 3U],
            xor_output[i + dist], xor_output[i + dist * 3U],
            log_m02, bytes);
    }
}


// Unrolled IFFT for encoder.  The compile-time accumulation choice keeps the
// scalar oracle's stage loop identical after optimization while sharing one
// source definition with the qualified vector schedule.
template<bool FuseAccumulation>
static void IFFT_DIT_Encoder_Impl(
    const backend::Ops& ops,
    const uint64_t bytes,
    const uint64_t source_policy_bytes,
    const void* const* data,
    const unsigned m_truncated,
    void** work,
    void** xor_result,
    const bool high_profile_transform,
    const unsigned m,
    const ffe_t* skewLUT)
{
    // Decimation in time: Unroll 2 layers at a time
    unsigned dist = 1, dist4 = 4;
    // Large AVX2 legacy-high blocks cross back below the mature copy-first
    // schedule once each shard no longer fits the useful cache region.  Apply
    // that policy to the first block as well as later accumulated blocks: the
    // first block accounts for most of the observed large-T regression.  The
    // low-profile interpolation remains source-staged.  Keep this deterministic
    // wire-independent crossover out of the byte loop.
    const bool stage_sources = m > 4 && !(
        high_profile_transform &&
        (ops.kind == LEO2_BACKEND_AVX2 ||
         ops.kind == LEO2_BACKEND_AVX512) &&
        m >= 256 && source_policy_bytes > 16U * 1024U);
    if (stage_sources)
    {
        // Consume complete encoder-source groups directly through the first
        // two inverse layers.  This replaces the former whole-active-prefix
        // copy plus first transform read/write with one read and one store.
        // Later layers still require the inactive suffix to contain zeros.
LEO_OPENMP_PARALLEL_FOR
        for (int i = m_truncated; i < (int)m; ++i)
            memset(work[i], 0, bytes);

        const unsigned full_group_count = m_truncated / 4;
LEO_OPENMP_PARALLEL_FOR
        for (int group = 0; group < (int)full_group_count; ++group)
        {
            const unsigned r = static_cast<unsigned>(group) * 4U;
#if defined(LEO2_ENABLE_TEST_HOOKS)
            TestHighIFFTButterfly4OutCalls.fetch_add(
                1, std::memory_order_relaxed);
#endif
            ops.ff16_ifft_butterfly4_out(
                data[r], data[r + 1], data[r + 2], data[r + 3],
                work[r], work[r + 1], work[r + 2], work[r + 3],
                skewLUT[r + 1], skewLUT[r + 3], skewLUT[r + 2],
                bytes);
        }

        const unsigned r = full_group_count * 4U;
        if (r < m_truncated)
        {
#if defined(LEO2_ENABLE_TEST_HOOKS)
            TestHighInputCopyShards.fetch_add(
                m_truncated - r, std::memory_order_relaxed);
#endif
            for (unsigned i = r; i < m_truncated; ++i)
                memcpy(work[i], data[i], bytes);
            IFFT_DIT4_Range(
                ops, bytes, work + r, 1,
                skewLUT[r + 1], skewLUT[r + 3], skewLUT[r + 2]);
        }
        dist = 4;
        dist4 = 16;
    }
    else
    {
        // Tiny transforms can fuse their only inverse stage with the final
        // accumulator.  The qualified large legacy-high AVX2 region above
        // also retains this established copy-first path.
#if defined(LEO2_ENABLE_TEST_HOOKS)
        TestHighInputCopyShards.fetch_add(
            m_truncated, std::memory_order_relaxed);
#endif
LEO_OPENMP_PARALLEL_FOR
        for (int i = 0; i < (int)m_truncated; ++i)
            memcpy(work[i], data[i], bytes);
LEO_OPENMP_PARALLEL_FOR
        for (int i = m_truncated; i < (int)m; ++i)
            memset(work[i], 0, bytes);
    }

    bool accumulated = false;
    for (; dist4 <= m; dist = dist4, dist4 <<= 2)
    {
        const bool accumulate_split_stage = xor_result &&
            FuseAccumulation && dist4 == m &&
            !UseFusedButterfly4(ops, bytes);
        // For each set of dist*4 elements:
LEO_OPENMP_PARALLEL_FOR
        for (int r = 0; r < (int)m_truncated; r += dist4)
        {
            const unsigned i_end = r + dist;
            const ffe_t log_m01 = skewLUT[i_end];
            const ffe_t log_m02 = skewLUT[i_end + dist];
            const ffe_t log_m23 = skewLUT[i_end + dist * 2];

            if (accumulate_split_stage)
            {
                IFFT_DIT4_xor_Split_Range(
                    ops, bytes, work + r, xor_result + r, dist,
                    log_m01, log_m23, log_m02);
            }
            else
            {
                IFFT_DIT4_Range(
                    ops, bytes, work + r, dist,
                    log_m01, log_m23, log_m02);
            }
        }
        if (accumulate_split_stage)
            accumulated = true;

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

        if (xor_result && FuseAccumulation)
        {
LEO_OPENMP_PARALLEL_FOR
            for (int i = 0; i < (int)dist; ++i)
            {
                IFFT_DIT2_xor(ops,
                    work[i], work[i + dist],
                    xor_result[i], xor_result[i + dist],
                    log_m, bytes);
            }
            accumulated = true;
        }
        else if (log_m == kModulus)
            VectorXOR_Threads(ops, bytes, dist, work + dist, work);
        else
        {
LEO_OPENMP_PARALLEL_FOR
            for (int i = 0; i < (int)dist; ++i)
                IFFT_DIT2(ops, work[i], work[i + dist], log_m, bytes);
        }
    }

    // I tried unrolling this but it does not provide more than 5% performance
    // improvement for 16-bit finite fields, so it's not worth the complexity.
    if (xor_result && !accumulated)
        VectorXOR_Threads(ops, bytes, m, xor_result, work);
}


static void IFFT_DIT_Encoder(
    const backend::Ops& ops,
    const uint64_t bytes,
    const uint64_t source_policy_bytes,
    const void* const* data,
    const unsigned m_truncated,
    void** work,
    void** xor_result,
    const bool high_profile_transform,
    const unsigned m,
    const ffe_t* skewLUT)
{
    // Matched-source crossover evidence qualifies the accumulating schedule
    // for the isolated x86 vector tables.  The scalar implementation remains
    // a correctness oracle and retains the materialize-then-XOR schedule,
    // which is faster for its table arithmetic.  Future native vector tables
    // require their own end-to-end qualification before opting in.
    if (ops.kind == LEO2_BACKEND_SSSE3 ||
        ops.kind == LEO2_BACKEND_AVX2 ||
        ops.kind == LEO2_BACKEND_AVX512)
    {
        IFFT_DIT_Encoder_Impl<true>(ops, bytes, source_policy_bytes,
            data, m_truncated,
            work, xor_result, high_profile_transform, m, skewLUT);
    }
    else
    {
        IFFT_DIT_Encoder_Impl<false>(ops, bytes, source_policy_bytes,
            data, m_truncated,
            work, xor_result, high_profile_transform, m, skewLUT);
    }
}


// Basic no-frills version for decoder.  When xor_result is non-null, qualified
// vector backends accumulate the final inverse layer directly into those
// disjoint buffers.  The source workspace is dead after that layer.  Scalar
// and compact fused-radix-four cases retain the established materialize-then-
// XOR schedule because that is their measured crossover winner.
template<bool FuseAccumulation>
static bool IFFT_DIT_DecoderImpl(
    const backend::Ops& ops,
    const uint64_t bytes,
    const unsigned m_truncated,
    void** work,
    const unsigned m,
    const ffe_t* skewLUT,
    void** xor_result)
{
    // Decimation in time: Unroll 2 layers at a time
    unsigned dist = 1, dist4 = 4;
    bool accumulated = false;
    for (; dist4 <= m; dist = dist4, dist4 <<= 2)
    {
        const bool accumulate_split_stage = xor_result &&
            FuseAccumulation && dist4 == m &&
            !UseFusedButterfly4(ops, bytes);
        // For each set of dist*4 elements:
LEO_OPENMP_PARALLEL_FOR
        for (int r = 0; r < (int)m_truncated; r += dist4)
        {
            const unsigned i_end = r + dist;
            const ffe_t log_m01 = skewLUT[i_end];
            const ffe_t log_m02 = skewLUT[i_end + dist];
            const ffe_t log_m23 = skewLUT[i_end + dist * 2];

            if (accumulate_split_stage)
            {
                IFFT_DIT4_xor_Split_Range(
                    ops, bytes, work + r, xor_result + r, dist,
                    log_m01, log_m23, log_m02);
            }
            else
            {
                IFFT_DIT4_Range(
                    ops, bytes, work + r, dist,
                    log_m01, log_m23, log_m02);
            }
        }
        if (accumulate_split_stage && m_truncated != 0)
            accumulated = true;
    }

    // If there is one layer left:
    if (dist < m)
    {
        // Assuming that dist = m / 2
        LEO_DEBUG_ASSERT(dist * 2 == m);

        const ffe_t log_m = skewLUT[dist];

        if (xor_result && FuseAccumulation)
        {
LEO_OPENMP_PARALLEL_FOR
            for (int i = 0; i < (int)dist; ++i)
            {
                IFFT_DIT2_xor(ops,
                    work[i], work[i + dist],
                    xor_result[i], xor_result[i + dist],
                    log_m, bytes);
            }
            accumulated = true;
        }
        else if (log_m == kModulus)
            VectorXOR_Threads(ops, bytes, dist, work + dist, work);
        else
        {
LEO_OPENMP_PARALLEL_FOR
            for (int i = 0; i < (int)dist; ++i)
            {
                IFFT_DIT2(
                    ops,
                    work[i],
                    work[i + dist],
                    log_m,
                    bytes);
            }
        }
    }

    if (xor_result && !accumulated)
    {
        if (m < 8)
            VectorXOR(ops, bytes, m, xor_result, work);
        else
            VectorXOR_Threads(ops, bytes, m, xor_result, work);
    }
    return accumulated;
}


static void IFFT_DIT_Decoder(
    const backend::Ops& ops,
    const uint64_t bytes,
    const unsigned m_truncated,
    void** work,
    const unsigned m,
    const ffe_t* skewLUT)
{
    (void)IFFT_DIT_DecoderImpl<false>(
        ops, bytes, m_truncated, work, m, skewLUT, NULL);
}


static bool IFFT_DIT_DecoderAccumulate(
    const backend::Ops& ops,
    const uint64_t bytes,
    const unsigned m_truncated,
    void** work,
    void** xor_result,
    const unsigned m,
    const ffe_t* skewLUT)
{
    LEO_DEBUG_ASSERT(xor_result != NULL);
    if (ops.kind == LEO2_BACKEND_SSSE3 ||
        ops.kind == LEO2_BACKEND_AVX2 ||
        ops.kind == LEO2_BACKEND_AVX512)
    {
        return IFFT_DIT_DecoderImpl<true>(
            ops, bytes, m_truncated, work, m, skewLUT, xor_result);
    }
    return IFFT_DIT_DecoderImpl<false>(
        ops, bytes, m_truncated, work, m, skewLUT, xor_result);
}


// GF16 source fusion was below the production threshold.  Keep the measured
// copy-first schedule while centralizing test accounting.  Algorithm 5 calls
// this for contributing later T-row blocks, and for raw block zero only when
// no later block contributes.
static void StageHighDecodeSources(
    const uint64_t bytes,
    const void* const* sources,
    void** work,
    const unsigned count)
{
#if defined(LEO2_ENABLE_TEST_HOOKS)
    uint64_t copy_count = 0;
    for (unsigned i = 0; i < count; ++i)
        copy_count += sources[i] != NULL;
    TestHighReceiveCopyShards.fetch_add(
        copy_count, std::memory_order_relaxed);
    TestHighReceiveZeroShards.fetch_add(
        count - copy_count, std::memory_order_relaxed);
#endif
LEO_OPENMP_PARALLEL_FOR
    for (int i = 0; i < (int)count; ++i)
    {
        if (sources[i])
            memcpy(work[i], sources[i], bytes);
        else
            memset(work[i], 0, bytes);
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
    const backend::Ops& ops,
    void * LEO_RESTRICT x, void * LEO_RESTRICT y,
    ffe_t log_m, uint64_t bytes)
{
#if defined(LEO_TRY_AVX2)
    if (&ops == &backend::GetDefaultOps() && CpuHasAVX2)
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
    if (&ops == &backend::GetDefaultOps() && CpuHasSSSE3)
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

    ops.ff16_fft_butterfly2(x, y, log_m, bytes);
}


// 4-way butterfly
static void FFT_DIT4(
    const backend::Ops& ops,
    uint64_t bytes,
    void** work,
    unsigned dist,
    const ffe_t log_m01,
    const ffe_t log_m23,
    const ffe_t log_m02)
{
    const bool use_fused_four_way = UseFusedButterfly4(ops, bytes);
#if defined(LEO2_ENABLE_TEST_HOOKS)
    TestFFTDIT4Calls.fetch_add(1, std::memory_order_relaxed);
    (use_fused_four_way ? TestFFTDIT4FusedCalls : TestFFTDIT4SplitCalls)
        .fetch_add(1, std::memory_order_relaxed);
#endif
#ifdef LEO_INTERLEAVE_BUTTERFLY4_OPT

#if defined(LEO_TRY_AVX2)

    if (use_fused_four_way &&
        &ops == &backend::GetDefaultOps() && CpuHasAVX2)
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
    if (use_fused_four_way &&
        &ops == &backend::GetDefaultOps() && CpuHasSSSE3)
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

    if (use_fused_four_way)
    {
        ops.ff16_fft_butterfly4(
            work[0], work[dist], work[dist * 2], work[dist * 3],
            log_m01, log_m23, log_m02, bytes);
        return;
    }

    // First layer:
    if (log_m02 == kModulus)
    {
        xor_mem(ops, work[dist * 2], work[0], bytes);
        xor_mem(ops, work[dist * 3], work[dist], bytes);
    }
    else
    {
        FFT_DIT2(ops, work[0], work[dist * 2], log_m02, bytes);
        FFT_DIT2(ops, work[dist], work[dist * 3], log_m02, bytes);
    }

    // Second layer:
    if (log_m01 == kModulus)
        xor_mem(ops, work[dist], work[0], bytes);
    else
        FFT_DIT2(ops, work[0], work[dist], log_m01, bytes);

    if (log_m23 == kModulus)
        xor_mem(ops, work[dist * 3], work[dist * 2], bytes);
    else
        FFT_DIT2(ops, work[dist * 2], work[dist * 3], log_m23, bytes);
}


static void FFT_DIT4_Range(
    const backend::Ops& ops,
    uint64_t bytes,
    void** work,
    unsigned dist,
    const ffe_t log_m01,
    const ffe_t log_m23,
    const ffe_t log_m02);


// In-place FFT for encoder and decoder
static void FFT_DIT(
    const backend::Ops& ops,
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
LEO_OPENMP_PARALLEL_FOR
        for (int r = 0; r < (int)m_truncated; r += dist4)
        {
            const unsigned i_end = r + dist;
            const ffe_t log_m01 = skewLUT[i_end];
            const ffe_t log_m02 = skewLUT[i_end + dist];
            const ffe_t log_m23 = skewLUT[i_end + dist * 2];

            FFT_DIT4_Range(
                ops, bytes, work + r, dist,
                log_m01, log_m23, log_m02);
        }
    }

    // If there is one layer left:
    if (dist4 == 2)
    {
LEO_OPENMP_PARALLEL_FOR
        for (int r = 0; r < (int)m_truncated; r += 2)
        {
            const ffe_t log_m = skewLUT[r + 1];

            if (log_m == kModulus)
                xor_mem(ops, work[r + 1], work[r], bytes);
            else
            {
                FFT_DIT2(
                    ops,
                    work[r],
                    work[r + 1],
                    log_m,
                    bytes);
            }
        }
    }
}

static void FFT_DIT4_Range(
    const backend::Ops& ops,
    uint64_t bytes,
    void** work,
    unsigned dist,
    const ffe_t log_m01,
    const ffe_t log_m23,
    const ffe_t log_m02)
{
    if (dist == 1)
    {
        FFT_DIT4(ops, bytes, work, dist,
            log_m01, log_m23, log_m02);
        return;
    }
#if defined(LEO_TRY_AVX2) || defined(LEO_TRY_SSSE3)
    if (&ops == &backend::GetDefaultOps())
    {
        for (unsigned i = 0; i < dist; ++i)
            FFT_DIT4(ops, bytes, work + i, dist,
                log_m01, log_m23, log_m02);
        return;
    }
#endif
    const bool use_fused_four_way = UseFusedButterfly4(ops, bytes);
#if defined(LEO2_ENABLE_TEST_HOOKS)
    TestFFTDIT4Calls.fetch_add(dist, std::memory_order_relaxed);
    (use_fused_four_way ? TestFFTDIT4FusedCalls : TestFFTDIT4SplitCalls)
        .fetch_add(dist, std::memory_order_relaxed);
#endif
    ops.ff16_fft_butterfly4_range(
        work, dist, log_m01, log_m23, log_m02,
        bytes, use_fused_four_way);
}

enum SourceEvaluationCallsite
{
    SourceEvaluationLowEncode,
    SourceEvaluationHighDecode
};

// Out-of-place first FFT layer(s) for reusable coefficient evaluation.  This
// keeps the coefficient block immutable while the backend combines the old
// copy with the first transform pass.  A complete dense low-encode block may
// also write its final layer directly to caller output.
static void FFT_DIT_FromCoefficients(
    const backend::Ops& ops,
    const uint64_t bytes,
    void* const* coefficients,
    void** evaluation_work,
    const unsigned m_truncated,
    const unsigned m,
    const ffe_t* skewLUT,
    void* const* final_outputs,
    SourceEvaluationCallsite callsite)
{
    (void)callsite;
    LEO_DEBUG_ASSERT(m >= 2);
    LEO_DEBUG_ASSERT(!final_outputs || m_truncated == m);
#if defined(LEO2_ENABLE_TEST_HOOKS)
    if (callsite == SourceEvaluationHighDecode)
        TestHighMatureOutputBlocks.fetch_add(
            1, std::memory_order_relaxed);
    if (callsite == SourceEvaluationHighDecode &&
        TestForceHighDecodeCopyFallback.load(std::memory_order_relaxed))
    {
LEO_OPENMP_PARALLEL_FOR
        for (int i = 0; i < (int)m; ++i)
            memcpy(evaluation_work[i], coefficients[i], bytes);
        TestHighCompatibilityCopyFallbacks.fetch_add(
            1, std::memory_order_relaxed);
        FFT_DIT(ops, bytes, evaluation_work, m_truncated, m, skewLUT);
        return;
    }
#endif
    if (m == 2)
    {
#if defined(LEO2_ENABLE_TEST_HOOKS)
        (callsite == SourceEvaluationLowEncode
            ? TestLowFFTButterfly2OutCalls
            : TestHighFFTButterfly2OutCalls).fetch_add(
                1, std::memory_order_relaxed);
#endif
        ops.ff16_fft_butterfly2_out(
            coefficients[0], coefficients[1],
            final_outputs ? final_outputs[0] : evaluation_work[0],
            final_outputs ? final_outputs[1] : evaluation_work[1],
            skewLUT[1], bytes);
        return;
    }

    const unsigned first_dist = m >> 2;
    const ffe_t log_m01 = skewLUT[first_dist];
    const ffe_t log_m02 = skewLUT[first_dist * 2];
    const ffe_t log_m23 = skewLUT[first_dist * 3];
    const bool use_fused_four_way = UseFusedButterfly4(ops, bytes);
    if (final_outputs && m == 4)
    {
        if (use_fused_four_way)
        {
#if defined(LEO2_ENABLE_TEST_HOOKS)
            (callsite == SourceEvaluationLowEncode
                ? TestLowFFTButterfly4OutCalls
                : TestHighFFTButterfly4OutCalls).fetch_add(
                    1, std::memory_order_relaxed);
#endif
            ops.ff16_fft_butterfly4_out(
                coefficients[0], coefficients[1], coefficients[2],
                coefficients[3], final_outputs[0], final_outputs[1],
                final_outputs[2], final_outputs[3], log_m01, log_m23,
                log_m02, bytes);
        }
        else
        {
#if defined(LEO2_ENABLE_TEST_HOOKS)
            (callsite == SourceEvaluationLowEncode
                ? TestLowFFTButterfly2OutCalls
                : TestHighFFTButterfly2OutCalls).fetch_add(
                    2, std::memory_order_relaxed);
            if (callsite == SourceEvaluationLowEncode)
                TestLowDirectFFTButterfly2OutCalls.fetch_add(
                    2, std::memory_order_relaxed);
#endif
            ops.ff16_fft_butterfly2_out(
                coefficients[0], coefficients[2], evaluation_work[0],
                evaluation_work[2], log_m02, bytes);
            ops.ff16_fft_butterfly2_out(
                coefficients[1], coefficients[3], evaluation_work[1],
                evaluation_work[3], log_m02, bytes);
            ops.ff16_fft_butterfly2_out(
                evaluation_work[0], evaluation_work[1], final_outputs[0],
                final_outputs[1], log_m01, bytes);
            ops.ff16_fft_butterfly2_out(
                evaluation_work[2], evaluation_work[3], final_outputs[2],
                final_outputs[3], log_m23, bytes);
        }
        return;
    }
LEO_OPENMP_PARALLEL_FOR
    for (int i = 0; i < (int)first_dist; ++i)
    {
        if (use_fused_four_way)
        {
#if defined(LEO2_ENABLE_TEST_HOOKS)
            (callsite == SourceEvaluationLowEncode
                ? TestLowFFTButterfly4OutCalls
                : TestHighFFTButterfly4OutCalls).fetch_add(
                1, std::memory_order_relaxed);
#endif
            ops.ff16_fft_butterfly4_out(
                coefficients[i], coefficients[i + first_dist],
                coefficients[i + first_dist * 2],
                coefficients[i + first_dist * 3],
                evaluation_work[i], evaluation_work[i + first_dist],
                evaluation_work[i + first_dist * 2],
                evaluation_work[i + first_dist * 3],
                log_m01, log_m23, log_m02, bytes);
            continue;
        }

        // Keep the established GF16 size/backend policy for the first two
        // layers.  The out-of-place two-way operations still combine the old
        // coefficient copy with the first layer; the second layer resumes
        // in-place on the evaluation workspace.
#if defined(LEO2_ENABLE_TEST_HOOKS)
        (callsite == SourceEvaluationLowEncode
            ? TestLowFFTButterfly2OutCalls
            : TestHighFFTButterfly2OutCalls).fetch_add(
            2, std::memory_order_relaxed);
#endif
        ops.ff16_fft_butterfly2_out(
            coefficients[i], coefficients[i + first_dist * 2],
            evaluation_work[i], evaluation_work[i + first_dist * 2],
            log_m02, bytes);
        ops.ff16_fft_butterfly2_out(
            coefficients[i + first_dist],
            coefficients[i + first_dist * 3],
            evaluation_work[i + first_dist],
            evaluation_work[i + first_dist * 3],
            log_m02, bytes);

        if (log_m01 == kModulus)
            xor_mem(ops, evaluation_work[i + first_dist],
                evaluation_work[i], bytes);
        else
            FFT_DIT2(ops, evaluation_work[i],
                evaluation_work[i + first_dist], log_m01, bytes);
        if (log_m23 == kModulus)
            xor_mem(ops, evaluation_work[i + first_dist * 3],
                evaluation_work[i + first_dist * 2], bytes);
        else
            FFT_DIT2(ops, evaluation_work[i + first_dist * 2],
                evaluation_work[i + first_dist * 3], log_m23, bytes);
    }

    unsigned dist4 = first_dist;
    unsigned dist = first_dist >> 2;
    for (; dist != 0; dist4 = dist, dist >>= 2)
    {
#if defined(LEO2_ENABLE_TEST_HOOKS)
        if (final_outputs && dist == 1 &&
            callsite == SourceEvaluationLowEncode)
        {
            if (use_fused_four_way)
                TestLowDirectFFTButterfly4OutCalls.fetch_add(
                    m / 4, std::memory_order_relaxed);
            else
                TestLowDirectFFTButterfly2OutCalls.fetch_add(
                    m / 2, std::memory_order_relaxed);
        }
#endif
LEO_OPENMP_PARALLEL_FOR
        for (int r = 0; r < (int)m_truncated; r += dist4)
        {
            const unsigned i_end = r + dist;
            const ffe_t remaining_log_m01 = skewLUT[i_end];
            const ffe_t remaining_log_m02 = skewLUT[i_end + dist];
            const ffe_t remaining_log_m23 = skewLUT[i_end + dist * 2];
            if (final_outputs && dist == 1)
            {
                if (use_fused_four_way)
                    ops.ff16_fft_butterfly4_out(
                        evaluation_work[r], evaluation_work[r + 1],
                        evaluation_work[r + 2], evaluation_work[r + 3],
                        final_outputs[r], final_outputs[r + 1],
                        final_outputs[r + 2], final_outputs[r + 3],
                        remaining_log_m01, remaining_log_m23,
                        remaining_log_m02, bytes);
                else
                {
                    if (remaining_log_m02 == kModulus)
                    {
                        xor_mem(ops, evaluation_work[r + 2],
                            evaluation_work[r], bytes);
                        xor_mem(ops, evaluation_work[r + 3],
                            evaluation_work[r + 1], bytes);
                    }
                    else
                    {
                        FFT_DIT2(ops, evaluation_work[r],
                            evaluation_work[r + 2], remaining_log_m02,
                            bytes);
                        FFT_DIT2(ops, evaluation_work[r + 1],
                            evaluation_work[r + 3], remaining_log_m02,
                            bytes);
                    }
                    ops.ff16_fft_butterfly2_out(
                        evaluation_work[r], evaluation_work[r + 1],
                        final_outputs[r], final_outputs[r + 1],
                        remaining_log_m01, bytes);
                    ops.ff16_fft_butterfly2_out(
                        evaluation_work[r + 2], evaluation_work[r + 3],
                        final_outputs[r + 2], final_outputs[r + 3],
                        remaining_log_m23, bytes);
                }
            }
            else
                FFT_DIT4_Range(
                    ops, bytes, evaluation_work + r, dist,
                    remaining_log_m01, remaining_log_m23,
                    remaining_log_m02);
        }
        if (final_outputs && dist == 1)
            return;
    }
    if (dist4 == 2)
    {
#if defined(LEO2_ENABLE_TEST_HOOKS)
        if (final_outputs && callsite == SourceEvaluationLowEncode)
            TestLowDirectFFTButterfly2OutCalls.fetch_add(
                m / 2, std::memory_order_relaxed);
#endif
LEO_OPENMP_PARALLEL_FOR
        for (int r = 0; r < (int)m_truncated; r += 2)
        {
            const ffe_t log_m = skewLUT[r + 1];
            if (final_outputs)
                ops.ff16_fft_butterfly2_out(
                    evaluation_work[r], evaluation_work[r + 1],
                    final_outputs[r], final_outputs[r + 1], log_m, bytes);
            else if (log_m == kModulus)
                xor_mem(ops, evaluation_work[r + 1], evaluation_work[r], bytes);
            else
                FFT_DIT2(ops, evaluation_work[r], evaluation_work[r + 1],
                    log_m, bytes);
        }
    }
}


//------------------------------------------------------------------------------
// Reed-Solomon Encode

void ReedSolomonEncodeWithSourcePolicy(
    const backend::Ops& ops,
    uint64_t buffer_bytes,
    uint64_t source_policy_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned requested_output_count,
    unsigned m,
    const void* const * data,
    void** work,
    const leopard2_internal::SparseForwardPlanBatchView* sparse_plans)
{
#if !defined(LEO2_ENABLE_TEST_HOOKS)
    (void)requested_output_count;
#endif
    // work <- IFFT(data, m, m)

    const ffe_t* skewLUT = FFTSkewStorage + m;

    IFFT_DIT_Encoder(
        ops,
        buffer_bytes,
        source_policy_bytes,
        data,
        original_count < m ? original_count : m,
        work,
        nullptr, // No xor output
        true, // Legacy-high transform
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
            ops,
            buffer_bytes,
            source_policy_bytes,
            data, // data source
            m,
            work + m, // temporary workspace
            work, // xor destination
            true, // Legacy-high transform
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
            ops,
            buffer_bytes,
            source_policy_bytes,
            data, // data source
            last_count,
            work + m, // temporary workspace
            work, // xor destination
            true, // Legacy-high transform
            m,
            skewLUT);
    }

skip_body:

    // work <- FFT(work, m, 0)
#if defined(LEO2_ENABLE_TEST_HOOKS)
    TestRecordSparseEncodeBlock(
        m, recovery_count, requested_output_count, sparse_plans, 0);
#endif
    if (sparse_plans && sparse_plans->block_count == 1)
    {
        const bool executed =
            leopard2_internal::ExecuteSparseForwardPlan(
                ops, buffer_bytes, m, 0, kModulus,
                sparse_plans->operation_masks,
                sparse_plans->operation_stride,
                PrunedMultiplierLog, FFTSkewStorage, work);
        LEO_DEBUG_ASSERT(executed);
        (void)executed;
    }
    else
    {
        FFT_DIT(
            ops,
            buffer_bytes,
            work,
            recovery_count,
            m,
            FFTSkewStorage);
    }
}


void ReedSolomonEncode(
    const backend::Ops& ops,
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned requested_output_count,
    unsigned m,
    const void* const * data,
    void** work,
    const leopard2_internal::SparseForwardPlanBatchView* sparse_plans)
{
    ReedSolomonEncodeWithSourcePolicy(
        ops, buffer_bytes, buffer_bytes, original_count,
        recovery_count, requested_output_count, m, data, work, sparse_plans);
}


void ReedSolomonEncode(
    const backend::Ops& ops,
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned m,
    const void* const* data,
    void** work)
{
    ReedSolomonEncode(ops, buffer_bytes, original_count, recovery_count,
        recovery_count, m, data, work, NULL);
}


void ReedSolomonEncode(
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned m,
    const void* const* data,
    void** work)
{
    ReedSolomonEncode(backend::GetDefaultOps(), buffer_bytes, original_count,
        recovery_count, m, data, work);
}


void ReedSolomonEncodeLow(
    const backend::Ops& ops,
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned p,
    const void* const* data,
    void* const* recovery,
    void** work,
    const leopard2_internal::SparseForwardPlanBatchView* sparse_plans)
{
    LEO_DEBUG_ASSERT(p >= 2 && p <= kOrder / 2);
    LEO_DEBUG_ASSERT(original_count <= p);
    LEO_DEBUG_ASSERT(recovery_count <= kOrder - p);

    // Interpolate the systematic prefix.  Coordinates original_count..p-1
    // are shortened and therefore supplied to the transform as known zeroes.
    IFFT_DIT_Encoder(
        ops,
        buffer_bytes,
        buffer_bytes,
        data,
        original_count,
        work,
        nullptr,
        false, // Low-profile interpolation
        p,
        FFTSkewStorage);

    // Evaluate the immutable coefficient block on each requested parity coset.
    // The out-of-place first layer writes work[p..2p) directly, avoiding the
    // former whole-P coefficient-copy pass for every nonempty parity block.
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

        const unsigned block = recovery_offset / p;
#if defined(LEO2_ENABLE_TEST_HOOKS)
        unsigned requested_outputs = 0;
        for (unsigned i = 0; i < requested_count; ++i)
            requested_outputs += recovery[recovery_offset + i] != NULL;
        TestRecordSparseEncodeBlock(
            p, requested_count, requested_outputs, sparse_plans, block);
#endif
        const bool use_sparse =
            sparse_plans && block < sparse_plans->block_count;
        bool direct_output = !use_sparse && requested_count == p;
        for (unsigned i = 0; direct_output && i < p; ++i)
            direct_output = recovery[recovery_offset + i] != NULL;
        if (use_sparse)
        {
            const bool executed =
                leopard2_internal::ExecuteSparseForwardPlanFromSources(
                    ops, buffer_bytes, p, p + recovery_offset, kModulus,
                    sparse_plans->operation_masks +
                        static_cast<size_t>(block) *
                            sparse_plans->operation_stride,
                    sparse_plans->operation_stride,
                    PrunedMultiplierLog, FFTSkewStorage, work, work + p);
            LEO_DEBUG_ASSERT(executed);
            (void)executed;
        }
        else
        {
            FFT_DIT_FromCoefficients(
                ops,
                buffer_bytes,
                work,
                work + p,
                requested_count,
                p,
                FFTSkewStorage + p + recovery_offset,
                direct_output ? recovery + recovery_offset : NULL,
                SourceEvaluationLowEncode);
        }

        if (direct_output)
        {
#if defined(LEO2_ENABLE_TEST_HOOKS)
            TestLowDirectOutputBlocks.fetch_add(
                1, std::memory_order_relaxed);
#endif
            continue;
        }
LEO_OPENMP_PARALLEL_FOR
        for (int i = 0; i < (int)requested_count; ++i)
        {
            void* const output = recovery[recovery_offset + i];
            if (output)
                memcpy(output, work[p + i], buffer_bytes);
        }
    }
}


void ReedSolomonEncodeLow(
    const backend::Ops& ops,
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned p,
    const void* const* data,
    void* const* recovery,
    void** work)
{
    ReedSolomonEncodeLow(ops, buffer_bytes, original_count, recovery_count,
        p, data, recovery, work, NULL);
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
    ReedSolomonEncodeLow(backend::GetDefaultOps(), buffer_bytes,
        original_count, recovery_count, p, data, recovery, work);
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


template<class OutputDependencies>
static void FFT_DIT_ErrorBits(
    const backend::Ops& ops,
    const uint64_t bytes,
    void** work,
    const unsigned n_truncated,
    const unsigned n,
    const ffe_t* skewLUT,
    const OutputDependencies& error_bits)
{
    unsigned mip_level = LastNonzeroBit32(n);

    // Decimation in time: Unroll 2 layers at a time
    unsigned dist4 = n, dist = n >> 2;
    for (; dist != 0; dist4 = dist, dist >>= 2, mip_level -=2)
    {
        // For each set of dist*4 elements:
LEO_OPENMP_PARALLEL_FOR
        for (int r = 0; r < (int)n_truncated; r += dist4)
        {
            if (!error_bits.IsNeeded(mip_level, r))
                continue;

            const unsigned i_end = r + dist;
            const ffe_t log_m01 = skewLUT[i_end];
            const ffe_t log_m02 = skewLUT[i_end + dist];
            const ffe_t log_m23 = skewLUT[i_end + dist * 2];

            // For each set of dist elements:
LEO_OPENMP_PARALLEL_FOR
            for (int i = r; i < (int)i_end; ++i)
            {
                FFT_DIT4(
                    ops,
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
LEO_OPENMP_PARALLEL_FOR
        for (int r = 0; r < (int)n_truncated; r += 2)
        {
            if (!error_bits.IsNeeded(mip_level, r))
                continue;

            const ffe_t log_m = skewLUT[r + 1];

            if (log_m == kModulus)
                xor_mem(ops, work[r + 1], work[r], bytes);
            else
            {
                FFT_DIT2(
                    ops,
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

    // Locator construction is setup work.  Keep it serial so creating a codec
    // or decode plan cannot instantiate a persistent OpenMP worker team.
    FWHT<false>(error_locations, kOrder, n);

    for (int i = 0; i < (int)kOrder; ++i)
    {
        error_locations[i] = static_cast<ffe_t>(
            ((unsigned)error_locations[i] * (unsigned)LogWalsh[i]) % kModulus);
    }

    FWHT<false>(error_locations, kOrder, kOrder);
    memcpy(locator_logs, error_locations, n * sizeof(ffe_t));
}


static void PrepareDecodeWalshActiveCombined(
    unsigned n,
    const uint8_t* erasures,
    const uint8_t* additional_erasures,
    ffe_t* locator_logs)
{
    LEO_DEBUG_ASSERT(n >= 2 && n <= kOrder);

    for (unsigned i = 0; i < n; ++i)
    {
        locator_logs[i] =
            (erasures[i] || (additional_erasures && additional_erasures[i])) ?
                1 : 0;
    }

    // This routine is reached only while constructing immutable setup state.
    // Execution-time transforms retain their parallel FWHT and byte kernels.
    FWHT<false>(locator_logs, n, n);

    const ffe_t* transformed_kernel =
        n == kOrder ? LogWalsh : ActiveLogWalsh + n - 2;
    for (unsigned i = 0; i < n; ++i)
    {
        locator_logs[i] = static_cast<ffe_t>(
            (static_cast<unsigned>(locator_logs[i]) *
                static_cast<unsigned>(transformed_kernel[i])) % kModulus);
    }

    FWHT<false>(locator_logs, n, n);
}


void PrepareDecodeWalshActive(
    unsigned n,
    const uint8_t* erasures,
    ffe_t* locator_logs)
{
    PrepareDecodeWalshActiveCombined(n, erasures, nullptr, locator_logs);
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
    // Pinned whole-locator sweeps show that the modulo-65535 Walsh butterflies
    // have a parent-size-dependent fixed cost.  Ambiguous cells remain direct:
    // the first active cell had at least a 10% observed margin across the frozen
    // layouts on the calibration host.  Parents of at most 32 points stay direct.
    unsigned direct_cutoff;
    if (n <= 32)
        direct_cutoff = n;
    else if (n == 64)
        direct_cutoff = 34;
    else if (n == 128)
        direct_cutoff = 24;
    else if (n == 256 || n == 512)
        direct_cutoff = 16;
    else
        direct_cutoff = 14;
    return erasure_count <= direct_cutoff;
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


void PrepareDecodeDirect(
    unsigned n,
    const uint8_t* erasures,
    ffe_t* locator_logs)
{
    LEO_DEBUG_ASSERT(n >= 2 && n <= kOrder);

    memset(locator_logs, 0, n * sizeof(ffe_t));
    AddDirectLocatorContributions(n, erasures, nullptr, locator_logs);
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
        PrepareDecodeWalshActive(n, erasures, locator_logs);
        return;
    }

    PrepareDecodeDirect(n, erasures, locator_logs);
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
        PrepareDecodeWalshActiveCombined(
            n, erasures, permanent_erasures, locator_logs);
        return;
    }

    memcpy(locator_logs, permanent_locator_logs, n * sizeof(ffe_t));
    AddDirectLocatorContributions(
        n, erasures, permanent_erasures, locator_logs);
}


// mul_mem() uses restrict-qualified input and output pointers.  Keep its
// contract intact while revealing a result in its existing coordinate slot.
static void mul_mem_inplace(
    const backend::Ops& ops,
    void* data,
    ffe_t log_m,
    uint64_t bytes)
{
    alignas(64) uint8_t source[64];
    uint8_t* output = reinterpret_cast<uint8_t*>(data);

    for (uint64_t offset = 0; offset < bytes; offset += sizeof(source))
    {
        memcpy(source, output + offset, sizeof(source));
        mul_mem(ops, output + offset, source, log_m, sizeof(source));
    }
}


static void AddFormalDerivative(
    const backend::Ops& ops,
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
            VectorXOR(ops, buffer_bytes, width, work + i - width, work + i);
        else
            VectorXOR_Threads(
                ops, buffer_bytes, width, work + i - width, work + i);
    }
}


void ReedSolomonDecodePrepared(
    const backend::Ops& ops,
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

LEO_OPENMP_PARALLEL_FOR
    for (int i = 0; i < (int)n; ++i)
    {
        if (coordinate_data[i])
            mul_mem(ops, work[i], coordinate_data[i], locator_logs[i],
                buffer_bytes);
        else
            memset(work[i], 0, buffer_bytes);
    }

    IFFT_DIT_Decoder(ops, buffer_bytes, input_count, work, n,
        FFTSkewStorage);

    AddFormalDerivative(ops, buffer_bytes, n, work);

#ifdef LEO_ERROR_BITFIELD_OPT
    FFT_DIT_ErrorBits(
        ops, buffer_bytes, work, n, n, FFTSkewStorage, error_bits);
#else
    FFT_DIT(ops, buffer_bytes, work, n, n, FFTSkewStorage);
#endif // LEO_ERROR_BITFIELD_OPT

LEO_OPENMP_PARALLEL_FOR
    for (int i = 0; i < (int)n; ++i)
    {
        if (requested_outputs[i])
            mul_mem_inplace(
                ops, work[i], kModulus - locator_logs[i], buffer_bytes);
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
    ReedSolomonDecodePrepared(backend::GetDefaultOps(), buffer_bytes, n,
        coordinate_data, requested_outputs, locator_logs, work);
}


void ReedSolomonDecodePlanned(
    const backend::Ops& ops,
    uint64_t buffer_bytes,
    unsigned n,
    const void* const* coordinate_data,
    unsigned input_count,
    const uint32_t* requested_coordinates,
    unsigned requested_count,
    const leopard2_internal::OutputDependencyView& output_dependencies,
    const ffe_t* locator_logs,
    bool reveal_outputs,
    void** work)
{
    LEO_DEBUG_ASSERT(n >= 2 && n <= kOrder);
    LEO_DEBUG_ASSERT(input_count <= n);

LEO_OPENMP_PARALLEL_FOR
    for (int i = 0; i < (int)n; ++i)
    {
        if (coordinate_data[i])
            mul_mem(ops, work[i], coordinate_data[i], locator_logs[i],
                buffer_bytes);
        else
            memset(work[i], 0, buffer_bytes);
    }

    if (input_count != 0)
        IFFT_DIT_Decoder(ops, buffer_bytes, input_count, work, n,
            FFTSkewStorage);

    AddFormalDerivative(ops, buffer_bytes, n, work);

#ifdef LEO_ERROR_BITFIELD_OPT
    FFT_DIT_ErrorBits(ops, buffer_bytes, work, n, n, FFTSkewStorage,
        output_dependencies);
#else
    (void)output_dependencies;
    FFT_DIT(ops, buffer_bytes, work, n, n, FFTSkewStorage);
#endif

    if (reveal_outputs)
    {
LEO_OPENMP_PARALLEL_FOR
        for (int i = 0; i < (int)requested_count; ++i)
        {
            const uint32_t coordinate = requested_coordinates[i];
            LEO_DEBUG_ASSERT(coordinate < n);
            mul_mem_inplace(ops, work[coordinate],
                kModulus - locator_logs[coordinate], buffer_bytes);
        }
    }
}


void ReedSolomonDecodePlanned(
    uint64_t buffer_bytes,
    unsigned n,
    const void* const* coordinate_data,
    unsigned input_count,
    const uint32_t* requested_coordinates,
    unsigned requested_count,
    const leopard2_internal::OutputDependencyView& output_dependencies,
    const ffe_t* locator_logs,
    void** work)
{
    ReedSolomonDecodePlanned(backend::GetDefaultOps(), buffer_bytes, n,
        coordinate_data, input_count, requested_coordinates, requested_count,
        output_dependencies, locator_logs, true, work);
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

static uint16_t PrunedMultiplierLog(
    const void* context,
    uint32_t storage_index)
{
    LEO_DEBUG_ASSERT(context != NULL && storage_index < kOrder);
    return static_cast<const ffe_t*>(context)[storage_index];
}


bool PrepareSparseForwardPlan(
    unsigned size,
    unsigned shift,
    uint8_t* dependency_workspace,
    size_t dependency_bytes,
    uint8_t* operation_masks,
    size_t retained_bytes,
    leopard2_internal::SparseForwardPlanStats& stats)
{
    return leopard2_internal::CompileSparseForwardPlan(
        kOrder, kModulus, size, shift,
        dependency_workspace, dependency_bytes,
        operation_masks, retained_bytes,
        PrunedMultiplierLog, FFTSkewStorage, stats);
}


bool ExecuteSparseForwardPlan(
    const backend::Ops& ops,
    uint64_t buffer_bytes,
    unsigned size,
    unsigned shift,
    const uint8_t* operation_masks,
    size_t retained_bytes,
    void** work)
{
    return leopard2_internal::ExecuteSparseForwardPlan(
        ops, buffer_bytes, size, shift, kModulus,
        operation_masks, retained_bytes,
        PrunedMultiplierLog, FFTSkewStorage, work);
}


bool ExecuteSparseForwardPlanFromSources(
    const backend::Ops& ops,
    uint64_t buffer_bytes,
    unsigned size,
    unsigned shift,
    const uint8_t* operation_masks,
    size_t retained_bytes,
    void* const* source,
    void** work)
{
    return leopard2_internal::ExecuteSparseForwardPlanFromSources(
        ops, buffer_bytes, size, shift, kModulus,
        operation_masks, retained_bytes,
        PrunedMultiplierLog, FFTSkewStorage, source, work);
}


bool PreparePrunedTransformPlan(
    unsigned size,
    unsigned shift,
    bool inverse,
    const uint8_t* input_mask,
    const uint8_t* output_mask,
    leopard2_internal::PrunedTransformPlan& plan)
{
    return leopard2_internal::CompilePrunedTransformPlan(
        kOrder, kModulus, size, shift, inverse,
        input_mask, output_mask, PrunedMultiplierLog,
        FFTSkewStorage, plan);
}


#if defined(LEO2_ENABLE_TEST_HOOKS)

void TestOnlyResetTransformCallsiteCounts()
{
    TestIFFTDIT4Calls.store(0, std::memory_order_relaxed);
    TestFFTDIT4Calls.store(0, std::memory_order_relaxed);
    TestIFFTDIT4FusedCalls.store(0, std::memory_order_relaxed);
    TestIFFTDIT4SplitCalls.store(0, std::memory_order_relaxed);
    TestFFTDIT4FusedCalls.store(0, std::memory_order_relaxed);
    TestFFTDIT4SplitCalls.store(0, std::memory_order_relaxed);
}


TestOnlyTransformCallsiteCounts TestOnlyGetTransformCallsiteCounts()
{
    TestOnlyTransformCallsiteCounts result;
    result.ifft_dit4 = TestIFFTDIT4Calls.load(std::memory_order_relaxed);
    result.fft_dit4 = TestFFTDIT4Calls.load(std::memory_order_relaxed);
    result.ifft_dit4_fused =
        TestIFFTDIT4FusedCalls.load(std::memory_order_relaxed);
    result.ifft_dit4_split =
        TestIFFTDIT4SplitCalls.load(std::memory_order_relaxed);
    result.fft_dit4_fused =
        TestFFTDIT4FusedCalls.load(std::memory_order_relaxed);
    result.fft_dit4_split =
        TestFFTDIT4SplitCalls.load(std::memory_order_relaxed);
    return result;
}


void TestOnlyResetLowEncodeCounts()
{
    TestLowFFTButterfly2OutCalls.store(0, std::memory_order_relaxed);
    TestLowFFTButterfly4OutCalls.store(0, std::memory_order_relaxed);
    TestLowDirectOutputBlocks.store(0, std::memory_order_relaxed);
    TestLowDirectFFTButterfly2OutCalls.store(0, std::memory_order_relaxed);
    TestLowDirectFFTButterfly4OutCalls.store(0, std::memory_order_relaxed);
}


TestOnlyLowEncodeCounts TestOnlyGetLowEncodeCounts()
{
    TestOnlyLowEncodeCounts result;
    result.fft_butterfly2_out_of_place =
        TestLowFFTButterfly2OutCalls.load(std::memory_order_relaxed);
    result.fft_butterfly4_out_of_place =
        TestLowFFTButterfly4OutCalls.load(std::memory_order_relaxed);
    result.direct_output_blocks =
        TestLowDirectOutputBlocks.load(std::memory_order_relaxed);
    result.direct_output_butterfly2_out_of_place =
        TestLowDirectFFTButterfly2OutCalls.load(std::memory_order_relaxed);
    result.direct_output_butterfly4_out_of_place =
        TestLowDirectFFTButterfly4OutCalls.load(std::memory_order_relaxed);
    return result;
}


void TestOnlyResetHighEncodeCounts()
{
    TestHighIFFTButterfly4OutCalls.store(0, std::memory_order_relaxed);
    TestHighInputCopyShards.store(0, std::memory_order_relaxed);
}


TestOnlyHighEncodeCounts TestOnlyGetHighEncodeCounts()
{
    TestOnlyHighEncodeCounts result;
    result.ifft_butterfly4_out_of_place =
        TestHighIFFTButterfly4OutCalls.load(std::memory_order_relaxed);
    result.input_copy_shards =
        TestHighInputCopyShards.load(std::memory_order_relaxed);
    return result;
}


void TestOnlyResetSparseEncodeCounts()
{
    TestSparseExactBlocks.store(0, std::memory_order_relaxed);
    TestSparsePrefixButterflies.store(0, std::memory_order_relaxed);
    TestSparseRetainedButterflies.store(0, std::memory_order_relaxed);
    TestSparseRequestedOutputCopies.store(0, std::memory_order_relaxed);
}


TestOnlySparseEncodeCounts TestOnlyGetSparseEncodeCounts()
{
    TestOnlySparseEncodeCounts result;
    result.exact_blocks =
        TestSparseExactBlocks.load(std::memory_order_relaxed);
    result.prefix_butterflies =
        TestSparsePrefixButterflies.load(std::memory_order_relaxed);
    result.retained_butterflies =
        TestSparseRetainedButterflies.load(std::memory_order_relaxed);
    result.requested_output_copies =
        TestSparseRequestedOutputCopies.load(std::memory_order_relaxed);
    return result;
}


void TestOnlyResetHighDecodeCounts()
{
    TestHighOutputBlocks.store(0, std::memory_order_relaxed);
    TestHighFFTButterfly2OutCalls.store(0, std::memory_order_relaxed);
    TestHighFFTButterfly4OutCalls.store(0, std::memory_order_relaxed);
    TestHighCompatibilityCopyFallbacks.store(0, std::memory_order_relaxed);
    TestHighPrunedOutputBlocks.store(0, std::memory_order_relaxed);
    TestHighMatureOutputBlocks.store(0, std::memory_order_relaxed);
    TestHighSyndromeAccumulatedBlocks.store(0, std::memory_order_relaxed);
    TestHighSyndromeMaterializedBlocks.store(0, std::memory_order_relaxed);
    TestHighSyndromePrunedAccumulatedBlocks.store(
        0, std::memory_order_relaxed);
    TestHighSyndromePrunedFallbackBlocks.store(
        0, std::memory_order_relaxed);
    TestHighSyndromeBlockZeroIFFTElisions.store(
        0, std::memory_order_relaxed);
    TestHighSyndromeForwardTransforms.store(0, std::memory_order_relaxed);
    TestHighSyndromeForwardTransformElisions.store(
        0, std::memory_order_relaxed);
    TestHighSyndromeBlockZeroXorShards.store(
        0, std::memory_order_relaxed);
    TestHighReceiveCopyShards.store(0, std::memory_order_relaxed);
    TestHighReceiveZeroShards.store(0, std::memory_order_relaxed);
}


void TestOnlySetHighDecodeCopyFallback(bool enabled)
{
    TestForceHighDecodeCopyFallback.store(enabled, std::memory_order_relaxed);
}


bool TestOnlyHighDecodeCopyFallbackEnabled()
{
    return TestForceHighDecodeCopyFallback.load(std::memory_order_relaxed);
}


TestOnlyHighDecodeCounts TestOnlyGetHighDecodeCounts()
{
    TestOnlyHighDecodeCounts result;
    result.output_blocks = TestHighOutputBlocks.load(std::memory_order_relaxed);
    result.fft_butterfly2_out_of_place =
        TestHighFFTButterfly2OutCalls.load(std::memory_order_relaxed);
    result.fft_butterfly4_out_of_place =
        TestHighFFTButterfly4OutCalls.load(std::memory_order_relaxed);
    result.compatibility_copy_fallbacks =
        TestHighCompatibilityCopyFallbacks.load(std::memory_order_relaxed);
    result.pruned_output_blocks =
        TestHighPrunedOutputBlocks.load(std::memory_order_relaxed);
    result.mature_output_blocks =
        TestHighMatureOutputBlocks.load(std::memory_order_relaxed);
    result.syndrome_accumulated_blocks =
        TestHighSyndromeAccumulatedBlocks.load(std::memory_order_relaxed);
    result.syndrome_materialized_blocks =
        TestHighSyndromeMaterializedBlocks.load(std::memory_order_relaxed);
    result.syndrome_pruned_accumulated_blocks =
        TestHighSyndromePrunedAccumulatedBlocks.load(
            std::memory_order_relaxed);
    result.syndrome_pruned_fallback_blocks =
        TestHighSyndromePrunedFallbackBlocks.load(
            std::memory_order_relaxed);
    result.syndrome_block_zero_ifft_elisions =
        TestHighSyndromeBlockZeroIFFTElisions.load(
            std::memory_order_relaxed);
    result.syndrome_forward_transforms =
        TestHighSyndromeForwardTransforms.load(std::memory_order_relaxed);
    result.syndrome_forward_transform_elisions =
        TestHighSyndromeForwardTransformElisions.load(
            std::memory_order_relaxed);
    result.syndrome_block_zero_xor_shards =
        TestHighSyndromeBlockZeroXorShards.load(
            std::memory_order_relaxed);
    result.receive_ifft_butterfly4_out_of_place = 0;
    result.receive_copy_shards =
        TestHighReceiveCopyShards.load(std::memory_order_relaxed);
    result.receive_zero_shards =
        TestHighReceiveZeroShards.load(std::memory_order_relaxed);
    return result;
}

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
    const backend::Ops& ops,
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
    FFT_DIT(ops, buffer_bytes, work, requested_output_count, size,
        FFTSkewStorage + shift);
}


void TestOnlyLchForward(
    uint64_t buffer_bytes,
    unsigned size,
    unsigned shift,
    unsigned requested_output_count,
    void** work)
{
    TestOnlyLchForward(
        backend::GetDefaultOps(), buffer_bytes, size, shift,
        requested_output_count, work);
}


void TestOnlyLchInverse(
    const backend::Ops& ops,
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
    IFFT_DIT_Decoder(ops, buffer_bytes, known_input_count, work, size,
        FFTSkewStorage + shift);
}


void TestOnlyLchInverse(
    uint64_t buffer_bytes,
    unsigned size,
    unsigned shift,
    unsigned known_input_count,
    void** work)
{
    TestOnlyLchInverse(
        backend::GetDefaultOps(), buffer_bytes, size, shift,
        known_input_count, work);
}


void TestOnlyAddFormalDerivative(
    uint64_t buffer_bytes,
    unsigned size,
    void** work)
{
    AddFormalDerivative(backend::GetDefaultOps(), buffer_bytes, size, work);
}

#endif // LEO2_ENABLE_TEST_HOOKS


void ReedSolomonDecodeLowPrepared(
    const backend::Ops& ops,
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

LEO_OPENMP_PARALLEL_FOR
    for (int i = 0; i < (int)n; ++i)
    {
        if (coordinate_data[i])
            mul_mem(ops, work[i], coordinate_data[i], locator_logs[i],
                buffer_bytes);
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
            ops,
            buffer_bytes, input_count, work + offset, p,
            FFTSkewStorage + offset);
    }

    // R10 Corollary 1 derivative term.  The omitted scalar-sum contribution
    // is f-hat^(0), which vanishes at every requested erasure because
    // f(e) * Lambda(e) = 0; leaving slot zero unchanged is therefore exact.
    AddFormalDerivative(ops, buffer_bytes, p, work);

    // Weighted block reduction by c_k / s_k(omega_(block*P)).
    // Each source block is dead after this reduction, so accumulate the
    // scaled value directly instead of scaling and XORing in two memory passes.
    for (unsigned block = 1; block < block_count; ++block)
    {
        const unsigned offset = block * p;
LEO_OPENMP_PARALLEL_FOR
        for (int i = 0; i < (int)p; ++i)
            muladd_mem(ops, work[i], work[offset + i],
                block_factors[block - 1], buffer_bytes);
    }

#ifdef LEO_ERROR_BITFIELD_OPT
    FFT_DIT_ErrorBits(
        ops, buffer_bytes, work, p, p, FFTSkewStorage, error_bits);
#else
    FFT_DIT(ops, buffer_bytes, work, p, p, FFTSkewStorage);
#endif

LEO_OPENMP_PARALLEL_FOR
    for (int i = 0; i < (int)p; ++i)
        if (requested_outputs[i])
            mul_mem_inplace(
                ops, work[i], kModulus - locator_logs[i], buffer_bytes);
}


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
    ReedSolomonDecodeLowPrepared(backend::GetDefaultOps(), buffer_bytes, n,
        p, coordinate_data, requested_outputs, locator_logs, block_factors,
        work);
}


static void ReedSolomonDecodeLowPrunedPlannedImpl(
    const backend::Ops& ops,
    uint64_t buffer_bytes,
    unsigned n,
    unsigned p,
    const void* const* coordinate_data,
    const uint16_t* block_input_counts,
    const uint32_t* requested_coordinates,
    unsigned requested_count,
    const leopard2_internal::OutputDependencyView& output_dependencies,
    const ffe_t* locator_logs,
    const ffe_t* block_factors,
    const leopard2_internal::PrunedTransformBlock* input_plans,
    unsigned input_plan_count,
    const leopard2_internal::PrunedTransformPlan* output_plan,
    bool reveal_outputs,
    void** work)
{
    LEO_DEBUG_ASSERT(p >= 2 && p <= n && n <= kOrder);
    LEO_DEBUG_ASSERT(input_plan_count == 0 || input_plans != NULL);

LEO_OPENMP_PARALLEL_FOR
    for (int i = 0; i < (int)n; ++i)
    {
        if (coordinate_data[i])
            mul_mem(ops, work[i], coordinate_data[i], locator_logs[i],
                buffer_bytes);
        else
            memset(work[i], 0, buffer_bytes);
    }

    const unsigned block_count = n / p;
    unsigned input_plan_index = 0;
    for (unsigned block = 0; block < block_count; ++block)
    {
        const unsigned input_count = block_input_counts[block];
        LEO_DEBUG_ASSERT(input_count <= p);
        if (input_count != 0)
        {
            const unsigned offset = block * p;
            const leopard2_internal::PrunedTransformPlan* pruned = NULL;
            if (input_plan_index < input_plan_count && input_plans &&
                input_plans[input_plan_index].block == block)
            {
                pruned = &input_plans[input_plan_index].plan;
                ++input_plan_index;
            }
            if (pruned)
            {
                LEO_DEBUG_ASSERT(pruned->size == p &&
                    pruned->shift == offset && pruned->inverse);
                const bool executed =
                    leopard2_internal::ExecutePrunedTransformPlan(
                        ops, buffer_bytes, *pruned, work + offset);
                LEO_DEBUG_ASSERT(executed);
                (void)executed;
            }
            else
                IFFT_DIT_Decoder(
                    ops, buffer_bytes, input_count, work + offset, p,
                    FFTSkewStorage + offset);
        }
    }
    LEO_DEBUG_ASSERT(input_plan_index == input_plan_count);

    AddFormalDerivative(ops, buffer_bytes, p, work);

    // The transformed nonzero block is dead after its weighted contribution.
    for (unsigned block = 1; block < block_count; ++block)
    {
        if (block_input_counts[block] == 0)
            continue;
        const unsigned offset = block * p;
LEO_OPENMP_PARALLEL_FOR
        for (int i = 0; i < (int)p; ++i)
            muladd_mem(ops, work[i], work[offset + i],
                block_factors[block - 1], buffer_bytes);
    }

    if (output_plan && output_plan->size != 0)
    {
        LEO_DEBUG_ASSERT(output_plan->size == p &&
            output_plan->shift == 0 && !output_plan->inverse);
        const bool executed = leopard2_internal::ExecutePrunedTransformPlan(
            ops, buffer_bytes, *output_plan, work);
        LEO_DEBUG_ASSERT(executed);
        (void)executed;
    }
    else
    {
#ifdef LEO_ERROR_BITFIELD_OPT
        FFT_DIT_ErrorBits(ops, buffer_bytes, work, p, p, FFTSkewStorage,
            output_dependencies);
#else
        (void)output_dependencies;
        FFT_DIT(ops, buffer_bytes, work, p, p, FFTSkewStorage);
#endif
    }

    if (reveal_outputs)
    {
LEO_OPENMP_PARALLEL_FOR
        for (int i = 0; i < (int)requested_count; ++i)
        {
            const uint32_t coordinate = requested_coordinates[i];
            LEO_DEBUG_ASSERT(coordinate < p);
            mul_mem_inplace(ops, work[coordinate],
                kModulus - locator_logs[coordinate], buffer_bytes);
        }
    }
}


void ReedSolomonDecodeLowPrunedPlanned(
    const backend::Ops& ops,
    uint64_t buffer_bytes,
    unsigned n,
    unsigned p,
    const void* const* coordinate_data,
    const uint16_t* block_input_counts,
    const uint32_t* requested_coordinates,
    unsigned requested_count,
    const leopard2_internal::OutputDependencyView& output_dependencies,
    const ffe_t* locator_logs,
    const ffe_t* block_factors,
    const leopard2_internal::PrunedTransformBlock* input_plans,
    unsigned input_plan_count,
    const leopard2_internal::PrunedTransformPlan* output_plan,
    void** work)
{
    ReedSolomonDecodeLowPrunedPlannedImpl(
        ops, buffer_bytes, n, p, coordinate_data, block_input_counts,
        requested_coordinates, requested_count, output_dependencies,
        locator_logs, block_factors, input_plans, input_plan_count,
        output_plan, true, work);
}


void ReedSolomonDecodeLowPrunedPlannedUnrevealed(
    const backend::Ops& ops,
    uint64_t buffer_bytes,
    unsigned n,
    unsigned p,
    const void* const* coordinate_data,
    const uint16_t* block_input_counts,
    const uint32_t* requested_coordinates,
    unsigned requested_count,
    const leopard2_internal::OutputDependencyView& output_dependencies,
    const ffe_t* locator_logs,
    const ffe_t* block_factors,
    const leopard2_internal::PrunedTransformBlock* input_plans,
    unsigned input_plan_count,
    const leopard2_internal::PrunedTransformPlan* output_plan,
    void** work)
{
    ReedSolomonDecodeLowPrunedPlannedImpl(
        ops, buffer_bytes, n, p, coordinate_data, block_input_counts,
        requested_coordinates, requested_count, output_dependencies,
        locator_logs, block_factors, input_plans, input_plan_count,
        output_plan, false, work);
}


void ReedSolomonDecodeLowPlanned(
    const backend::Ops& ops,
    uint64_t buffer_bytes,
    unsigned n,
    unsigned p,
    const void* const* coordinate_data,
    const uint16_t* block_input_counts,
    const uint32_t* requested_coordinates,
    unsigned requested_count,
    const leopard2_internal::OutputDependencyView& output_dependencies,
    const ffe_t* locator_logs,
    const ffe_t* block_factors,
    void** work)
{
    ReedSolomonDecodeLowPrunedPlanned(
        ops, buffer_bytes, n, p, coordinate_data, block_input_counts,
        requested_coordinates, requested_count, output_dependencies,
        locator_logs, block_factors, NULL, 0, NULL, work);
}


void ReedSolomonDecodeLowPlanned(
    uint64_t buffer_bytes,
    unsigned n,
    unsigned p,
    const void* const* coordinate_data,
    const uint16_t* block_input_counts,
    const uint32_t* requested_coordinates,
    unsigned requested_count,
    const leopard2_internal::OutputDependencyView& output_dependencies,
    const ffe_t* locator_logs,
    const ffe_t* block_factors,
    void** work)
{
    ReedSolomonDecodeLowPlanned(backend::GetDefaultOps(), buffer_bytes, n, p,
        coordinate_data, block_input_counts, requested_coordinates,
        requested_count, output_dependencies, locator_logs, block_factors,
        work);
}


static void ReedSolomonDecodeLowTiledPrunedPlannedImpl(
    const backend::Ops& ops,
    uint64_t buffer_bytes,
    unsigned n,
    unsigned p,
    const void* const* coordinate_data,
    const uint16_t* block_input_counts,
    const uint32_t* requested_coordinates,
    unsigned requested_count,
    const leopard2_internal::OutputDependencyView& output_dependencies,
    const ffe_t* locator_logs,
    const ffe_t* block_factors,
    const leopard2_internal::PrunedTransformBlock* input_plans,
    unsigned input_plan_count,
    const leopard2_internal::PrunedTransformPlan* output_plan,
    bool reveal_outputs,
    void** work)
{
    LEO_DEBUG_ASSERT(p >= 2 && p <= n && n <= kOrder);
    LEO_DEBUG_ASSERT(n % p == 0);
    LEO_DEBUG_ASSERT(input_plan_count == 0 || input_plans != NULL);

    void** const accumulator = work;
    void** const tile = work + p;
    const unsigned block_count = n / p;

    const unsigned first_input_count = block_input_counts[0];
    LEO_DEBUG_ASSERT(first_input_count <= p);
LEO_OPENMP_PARALLEL_FOR
    for (int i = 0; i < (int)p; ++i)
    {
        if (coordinate_data[i])
            mul_mem(
                ops, accumulator[i], coordinate_data[i], locator_logs[i],
                buffer_bytes);
        else
            memset(accumulator[i], 0, buffer_bytes);
    }
    unsigned input_plan_index = 0;
    if (first_input_count != 0)
    {
        const leopard2_internal::PrunedTransformPlan* pruned = NULL;
        if (input_plan_index < input_plan_count && input_plans &&
            input_plans[input_plan_index].block == 0)
        {
            pruned = &input_plans[input_plan_index].plan;
            ++input_plan_index;
        }
        if (pruned)
        {
            LEO_DEBUG_ASSERT(pruned->size == p &&
                pruned->shift == 0 && pruned->inverse);
            const bool executed =
                leopard2_internal::ExecutePrunedTransformPlan(
                    ops, buffer_bytes, *pruned, accumulator);
            LEO_DEBUG_ASSERT(executed);
            (void)executed;
        }
        else
            IFFT_DIT_Decoder(
                ops, buffer_bytes, first_input_count, accumulator, p,
                FFTSkewStorage);
    }

    AddFormalDerivative(ops, buffer_bytes, p, accumulator);

    // Reuse the tile as an immutable multiply-add source for one memory pass.
    for (unsigned block = 1; block < block_count; ++block)
    {
        const unsigned input_count = block_input_counts[block];
        LEO_DEBUG_ASSERT(input_count <= p);
        if (input_count == 0)
            continue;
        const unsigned offset = block * p;
LEO_OPENMP_PARALLEL_FOR
        for (int i = 0; i < (int)p; ++i)
        {
            const unsigned coordinate = offset + i;
            if (coordinate_data[coordinate])
                mul_mem(
                    ops, tile[i], coordinate_data[coordinate],
                    locator_logs[coordinate], buffer_bytes);
            else
                memset(tile[i], 0, buffer_bytes);
        }
        const leopard2_internal::PrunedTransformPlan* pruned = NULL;
        if (input_plan_index < input_plan_count && input_plans &&
            input_plans[input_plan_index].block == block)
        {
            pruned = &input_plans[input_plan_index].plan;
            ++input_plan_index;
        }
        if (pruned)
        {
            LEO_DEBUG_ASSERT(pruned->size == p &&
                pruned->shift == offset && pruned->inverse);
            const bool executed =
                leopard2_internal::ExecutePrunedTransformPlan(
                    ops, buffer_bytes, *pruned, tile);
            LEO_DEBUG_ASSERT(executed);
            (void)executed;
        }
        else
            IFFT_DIT_Decoder(
                ops, buffer_bytes, input_count, tile, p,
                FFTSkewStorage + offset);
LEO_OPENMP_PARALLEL_FOR
        for (int i = 0; i < (int)p; ++i)
            muladd_mem(ops, accumulator[i], tile[i],
                block_factors[block - 1], buffer_bytes);
    }
    LEO_DEBUG_ASSERT(input_plan_index == input_plan_count);

    if (output_plan && output_plan->size != 0)
    {
        LEO_DEBUG_ASSERT(output_plan->size == p &&
            output_plan->shift == 0 && !output_plan->inverse);
        const bool executed = leopard2_internal::ExecutePrunedTransformPlan(
            ops, buffer_bytes, *output_plan, accumulator);
        LEO_DEBUG_ASSERT(executed);
        (void)executed;
    }
    else
    {
#ifdef LEO_ERROR_BITFIELD_OPT
        FFT_DIT_ErrorBits(
            ops, buffer_bytes, accumulator, p, p, FFTSkewStorage,
            output_dependencies);
#else
        (void)output_dependencies;
        FFT_DIT(ops, buffer_bytes, accumulator, p, p, FFTSkewStorage);
#endif
    }

    if (reveal_outputs)
    {
LEO_OPENMP_PARALLEL_FOR
        for (int i = 0; i < (int)requested_count; ++i)
        {
            const uint32_t coordinate = requested_coordinates[i];
            LEO_DEBUG_ASSERT(coordinate < p);
            mul_mem_inplace(
                ops, accumulator[coordinate],
                kModulus - locator_logs[coordinate], buffer_bytes);
        }
    }
}


void ReedSolomonDecodeLowTiledPrunedPlanned(
    const backend::Ops& ops,
    uint64_t buffer_bytes,
    unsigned n,
    unsigned p,
    const void* const* coordinate_data,
    const uint16_t* block_input_counts,
    const uint32_t* requested_coordinates,
    unsigned requested_count,
    const leopard2_internal::OutputDependencyView& output_dependencies,
    const ffe_t* locator_logs,
    const ffe_t* block_factors,
    const leopard2_internal::PrunedTransformBlock* input_plans,
    unsigned input_plan_count,
    const leopard2_internal::PrunedTransformPlan* output_plan,
    void** work)
{
    ReedSolomonDecodeLowTiledPrunedPlannedImpl(
        ops, buffer_bytes, n, p, coordinate_data, block_input_counts,
        requested_coordinates, requested_count, output_dependencies,
        locator_logs, block_factors, input_plans, input_plan_count,
        output_plan, true, work);
}


void ReedSolomonDecodeLowTiledPrunedPlannedUnrevealed(
    const backend::Ops& ops,
    uint64_t buffer_bytes,
    unsigned n,
    unsigned p,
    const void* const* coordinate_data,
    const uint16_t* block_input_counts,
    const uint32_t* requested_coordinates,
    unsigned requested_count,
    const leopard2_internal::OutputDependencyView& output_dependencies,
    const ffe_t* locator_logs,
    const ffe_t* block_factors,
    const leopard2_internal::PrunedTransformBlock* input_plans,
    unsigned input_plan_count,
    const leopard2_internal::PrunedTransformPlan* output_plan,
    void** work)
{
    ReedSolomonDecodeLowTiledPrunedPlannedImpl(
        ops, buffer_bytes, n, p, coordinate_data, block_input_counts,
        requested_coordinates, requested_count, output_dependencies,
        locator_logs, block_factors, input_plans, input_plan_count,
        output_plan, false, work);
}


void ReedSolomonDecodeLowTiledPlanned(
    const backend::Ops& ops,
    uint64_t buffer_bytes,
    unsigned n,
    unsigned p,
    const void* const* coordinate_data,
    const uint16_t* block_input_counts,
    const uint32_t* requested_coordinates,
    unsigned requested_count,
    const leopard2_internal::OutputDependencyView& output_dependencies,
    const ffe_t* locator_logs,
    const ffe_t* block_factors,
    void** work)
{
    ReedSolomonDecodeLowTiledPrunedPlanned(
        ops, buffer_bytes, n, p, coordinate_data, block_input_counts,
        requested_coordinates, requested_count, output_dependencies,
        locator_logs, block_factors, NULL, 0, NULL, work);
}


// Algorithm 5 line 3 contributes IFFT_0(F_0) for the zero-shift receive
// block, and line 4 immediately applies FFT_0 to the complete sum.  Linearity
// gives FFT_0(IFFT_0(F_0) + H_later) = F_0 + FFT_0(H_later).  Preserve the
// immutable received values F_0, transform only later-block coefficients, and
// stage F_0 directly when no later block contributes.
static void FinishHighDecodeSyndrome(
    const backend::Ops& ops,
    uint64_t buffer_bytes,
    const void* const* block_zero_sources,
    unsigned block_zero_input_count,
    bool have_later_contribution,
    unsigned t,
    void** work)
{
    (void)block_zero_input_count;
#if defined(LEO2_ENABLE_TEST_HOOKS)
    if (block_zero_input_count != 0)
        TestHighSyndromeBlockZeroIFFTElisions.fetch_add(
            1, std::memory_order_relaxed);
#endif
    if (!have_later_contribution)
    {
        StageHighDecodeSources(
            buffer_bytes, block_zero_sources, work, t);
#if defined(LEO2_ENABLE_TEST_HOOKS)
        TestHighSyndromeForwardTransformElisions.fetch_add(
            1, std::memory_order_relaxed);
#endif
        return;
    }

    FFT_DIT(ops, buffer_bytes, work, t, t, FFTSkewStorage);
#if defined(LEO2_ENABLE_TEST_HOOKS)
    TestHighSyndromeForwardTransforms.fetch_add(
        1, std::memory_order_relaxed);
    uint64_t xor_shards = 0;
    for (unsigned i = 0; i < t; ++i)
        xor_shards += block_zero_sources[i] != NULL;
    TestHighSyndromeBlockZeroXorShards.fetch_add(
        xor_shards, std::memory_order_relaxed);
#endif
LEO_OPENMP_PARALLEL_FOR
    for (int i = 0; i < (int)t; ++i)
        if (block_zero_sources[i])
            xor_mem(ops, work[i], block_zero_sources[i], buffer_bytes);
}


void ReedSolomonDecodeHighPrepared(
    const backend::Ops& ops,
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

    const unsigned block_count = n / t;
    unsigned block_zero_input_count = t;
    while (block_zero_input_count > 0 &&
           !coordinate_data[block_zero_input_count - 1])
        --block_zero_input_count;
    bool have_later_contribution = false;
    for (unsigned block = 1; block < block_count; ++block)
    {
        const unsigned offset = block * t;
        unsigned input_count = t;
        while (input_count > 0 && !coordinate_data[offset + input_count - 1])
            --input_count;
        if (input_count == 0)
            continue;
        void** const block_work = have_later_contribution
            ? work + offset : work;
        StageHighDecodeSources(
            buffer_bytes, coordinate_data + offset, block_work, t);
        if (!have_later_contribution)
        {
            IFFT_DIT_Decoder(
                ops, buffer_bytes, input_count, block_work, t,
                FFTSkewStorage + offset);
#if defined(LEO2_ENABLE_TEST_HOOKS)
            TestHighSyndromeMaterializedBlocks.fetch_add(
                1, std::memory_order_relaxed);
#endif
        }
        else
        {
            const bool accumulated = IFFT_DIT_DecoderAccumulate(
                ops, buffer_bytes, input_count, block_work, work, t,
                FFTSkewStorage + offset);
#if defined(LEO2_ENABLE_TEST_HOOKS)
            (accumulated ? TestHighSyndromeAccumulatedBlocks
                         : TestHighSyndromeMaterializedBlocks)
                .fetch_add(1, std::memory_order_relaxed);
#else
            (void)accumulated;
#endif
        }
        have_later_contribution = true;
    }

    // h on V_t, then z = h * Lambda on V_t.
    FinishHighDecodeSyndrome(
        ops, buffer_bytes, coordinate_data, block_zero_input_count,
        have_later_contribution, t, work);
LEO_OPENMP_PARALLEL_FOR
    for (int i = 0; i < (int)t; ++i)
    {
        if (coordinate_data[i])
            mul_mem_inplace(ops, work[i], locator_logs[i], buffer_bytes);
        else
            memset(work[i], 0, buffer_bytes);
    }
    IFFT_DIT_Decoder(ops, buffer_bytes, t, work, t, FFTSkewStorage);

    // Evaluate z only on message blocks that contain a requested original.
    for (unsigned block = 1; block < block_count; ++block)
    {
        const unsigned offset = block * t;
        unsigned requested_count = t;
        while (requested_count > 0 && !requested_outputs[offset + requested_count - 1])
            --requested_count;
        if (requested_count == 0)
            continue;
#if defined(LEO2_ENABLE_TEST_HOOKS)
        TestHighOutputBlocks.fetch_add(1, std::memory_order_relaxed);
#endif
        FFT_DIT_FromCoefficients(
            ops, buffer_bytes, work, work + offset, requested_count, t,
            FFTSkewStorage + offset, NULL, SourceEvaluationHighDecode);
LEO_OPENMP_PARALLEL_FOR
        for (int i = 0; i < (int)requested_count; ++i)
        {
            const unsigned coordinate = offset + i;
            if (requested_outputs[coordinate])
            {
                const ffe_t reveal_log = SubMod(
                    output_factors[coordinate], locator_logs[coordinate]);
                mul_mem_inplace(
                    ops, work[coordinate], reveal_log, buffer_bytes);
            }
        }
    }
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
    ReedSolomonDecodeHighPrepared(backend::GetDefaultOps(), buffer_bytes, n,
        t, coordinate_data, requested_outputs, locator_logs, output_factors,
        work);
}


static bool ShouldUsePrunedHighSyndromeSink(
    const backend::Ops& ops,
    uint64_t buffer_bytes)
{
    // Pinned production ABBA measurements promote AVX2 from the safe 1-KiB
    // crossover.  SSSE3 was faster, but its 64-KiB target cell improved only
    // 4.86%, below the project's default 5% complexity threshold, so it keeps
    // the mature materialized schedule.
    return (ops.kind == LEO2_BACKEND_AVX2 ||
            ops.kind == LEO2_BACKEND_AVX512) &&
        buffer_bytes >= 1024;
}


void ReedSolomonDecodeHighPrunedPlanned(
    const backend::Ops& ops,
    uint64_t buffer_bytes,
    unsigned n,
    unsigned t,
    const void* const* coordinate_data,
    const uint16_t* block_input_counts,
    const uint32_t* requested_coordinates,
    const leopard2_internal::DecodeOutputBlock* output_blocks,
    unsigned output_block_count,
    const ffe_t* locator_logs,
    const ffe_t* output_factors,
    const leopard2_internal::PrunedTransformBlock* input_plans,
    unsigned input_plan_count,
    const leopard2_internal::PrunedTransformPlan* output_plans,
    unsigned output_plan_count,
    void* const* systematic_output,
    void** work)
{
    LEO_DEBUG_ASSERT(t >= 2 && t < n && n <= kOrder);
    LEO_DEBUG_ASSERT(input_plan_count == 0 || input_plans != NULL);
    LEO_DEBUG_ASSERT(output_plan_count == 0 || output_plans != NULL);
    LEO_DEBUG_ASSERT(output_plan_count == 0 ||
        output_plan_count == output_block_count);

    const unsigned block_count = n / t;
    unsigned input_plan_index = 0;
    if (input_plan_index < input_plan_count && input_plans &&
        input_plans[input_plan_index].block == 0)
    {
        LEO_DEBUG_ASSERT(input_plans[input_plan_index].plan.size == t &&
            input_plans[input_plan_index].plan.shift == 0 &&
            input_plans[input_plan_index].plan.inverse);
        ++input_plan_index;
    }
    bool have_later_contribution = false;
    for (unsigned block = 1; block < block_count; ++block)
    {
        const unsigned offset = block * t;
        const unsigned input_count = block_input_counts[block];
        LEO_DEBUG_ASSERT(input_count <= t);
        if (input_count == 0)
            continue;
        const leopard2_internal::PrunedTransformPlan* pruned = NULL;
        if (input_plan_index < input_plan_count && input_plans &&
            input_plans[input_plan_index].block == block)
        {
            pruned = &input_plans[input_plan_index].plan;
            ++input_plan_index;
        }
        void** const block_work = have_later_contribution
            ? work + offset : work;
        StageHighDecodeSources(
            buffer_bytes, coordinate_data + offset, block_work, t);
        if (pruned)
        {
            LEO_DEBUG_ASSERT(pruned->size == t &&
                pruned->shift == offset && pruned->inverse);
            const bool use_accumulating_sink =
                ShouldUsePrunedHighSyndromeSink(ops, buffer_bytes);
            if (have_later_contribution && use_accumulating_sink)
            {
                const bool accumulated = leopard2_internal::
                    ExecutePrunedInverseTransformPlanAccumulate(
                        ops, buffer_bytes, *pruned,
                        block_work, work);
                if (accumulated)
                {
#if defined(LEO2_ENABLE_TEST_HOOKS)
                    TestHighSyndromeAccumulatedBlocks.fetch_add(
                        1, std::memory_order_relaxed);
                    TestHighSyndromePrunedAccumulatedBlocks.fetch_add(
                        1, std::memory_order_relaxed);
#endif
                    have_later_contribution = true;
                    continue;
                }
            }
#if defined(LEO2_ENABLE_TEST_HOOKS)
            if (have_later_contribution && use_accumulating_sink)
                TestHighSyndromePrunedFallbackBlocks.fetch_add(
                    1, std::memory_order_relaxed);
#endif
            const bool executed =
                leopard2_internal::ExecutePrunedTransformPlan(
                    ops, buffer_bytes, *pruned, block_work);
            LEO_DEBUG_ASSERT(executed);
            (void)executed;
            if (have_later_contribution)
            {
                if (t < 8)
                    VectorXOR(ops, buffer_bytes, t, work, block_work);
                else
                    VectorXOR_Threads(
                        ops, buffer_bytes, t, work, block_work);
            }
#if defined(LEO2_ENABLE_TEST_HOOKS)
            TestHighSyndromeMaterializedBlocks.fetch_add(
                1, std::memory_order_relaxed);
#endif
        }
        else
        {
            bool accumulated = false;
            if (have_later_contribution)
                accumulated = IFFT_DIT_DecoderAccumulate(
                    ops, buffer_bytes, input_count, block_work, work, t,
                    FFTSkewStorage + offset);
            else
                IFFT_DIT_Decoder(
                    ops, buffer_bytes, input_count, block_work, t,
                    FFTSkewStorage + offset);
#if defined(LEO2_ENABLE_TEST_HOOKS)
            (have_later_contribution && accumulated
                ? TestHighSyndromeAccumulatedBlocks
                : TestHighSyndromeMaterializedBlocks)
                .fetch_add(1, std::memory_order_relaxed);
#else
            (void)accumulated;
#endif
        }
        have_later_contribution = true;
    }
    LEO_DEBUG_ASSERT(input_plan_index == input_plan_count);

    FinishHighDecodeSyndrome(
        ops, buffer_bytes, coordinate_data, block_input_counts[0],
        have_later_contribution, t, work);
LEO_OPENMP_PARALLEL_FOR
    for (int i = 0; i < (int)t; ++i)
    {
        if (coordinate_data[i])
            mul_mem_inplace(ops, work[i], locator_logs[i], buffer_bytes);
        else
            memset(work[i], 0, buffer_bytes);
    }
    IFFT_DIT_Decoder(ops, buffer_bytes, t, work, t, FFTSkewStorage);

    for (unsigned output_block = 0;
         output_block < output_block_count;
         ++output_block)
    {
        const leopard2_internal::DecodeOutputBlock& descriptor =
            output_blocks[output_block];
        LEO_DEBUG_ASSERT(descriptor.block > 0 && descriptor.block < block_count);
        LEO_DEBUG_ASSERT(descriptor.requested_prefix > 0 &&
            descriptor.requested_prefix <= t);
        LEO_DEBUG_ASSERT(descriptor.requested_begin < descriptor.requested_end);
        const unsigned offset = descriptor.block * t;
#if defined(LEO2_ENABLE_TEST_HOOKS)
        TestHighOutputBlocks.fetch_add(1, std::memory_order_relaxed);
#endif
        const leopard2_internal::PrunedTransformPlan* pruned =
            output_plans && output_plan_count == output_block_count &&
            output_plans[output_block].size != 0
                ? &output_plans[output_block]
                : NULL;
        if (pruned)
        {
#if defined(LEO2_ENABLE_TEST_HOOKS)
            TestHighPrunedOutputBlocks.fetch_add(
                1, std::memory_order_relaxed);
#endif
            LEO_DEBUG_ASSERT(pruned->size == t &&
                pruned->shift == offset && !pruned->inverse);
#if defined(LEO2_ENABLE_TEST_HOOKS)
            const bool force_copy_fallback =
                TestForceHighDecodeCopyFallback.load(
                std::memory_order_relaxed);
            if (!force_copy_fallback)
                TestHighFFTButterfly2OutCalls.fetch_add(
                    t / 2, std::memory_order_relaxed);
            const bool executed = !force_copy_fallback &&
                leopard2_internal::ExecutePrunedForwardTransformPlanFromSources(
                    ops, buffer_bytes, *pruned, work, work + offset);
            LEO_DEBUG_ASSERT(executed || force_copy_fallback);
#else
            const bool executed =
                leopard2_internal::ExecutePrunedForwardTransformPlanFromSources(
                    ops, buffer_bytes, *pruned, work, work + offset);
            LEO_DEBUG_ASSERT(executed);
#endif
            if (!executed)
            {
#if defined(LEO2_ENABLE_TEST_HOOKS)
                TestHighCompatibilityCopyFallbacks.fetch_add(
                    1, std::memory_order_relaxed);
#endif
LEO_OPENMP_PARALLEL_FOR
                for (int i = 0; i < (int)t; ++i)
                    memcpy(work[offset + i], work[i], buffer_bytes);
                const bool fallback =
                    leopard2_internal::ExecutePrunedTransformPlan(
                        ops, buffer_bytes, *pruned, work + offset);
                LEO_DEBUG_ASSERT(fallback);
                (void)fallback;
            }
        }
        else
            FFT_DIT_FromCoefficients(
                ops, buffer_bytes, work, work + offset,
                descriptor.requested_prefix, t,
                FFTSkewStorage + offset, NULL, SourceEvaluationHighDecode);
LEO_OPENMP_PARALLEL_FOR
        for (int i = (int)descriptor.requested_begin;
             i < (int)descriptor.requested_end;
             ++i)
        {
            const uint32_t coordinate = requested_coordinates[i];
            LEO_DEBUG_ASSERT(coordinate >= offset && coordinate < offset + t);
            const ffe_t reveal_log = SubMod(
                output_factors[coordinate], locator_logs[coordinate]);
            if (systematic_output)
            {
                LEO_DEBUG_ASSERT(coordinate >= t);
                void* const destination = systematic_output[coordinate - t];
                LEO_DEBUG_ASSERT(destination != NULL);
                mul_mem(ops, destination, work[coordinate], reveal_log,
                    buffer_bytes);
            }
            else
                mul_mem_inplace(
                    ops, work[coordinate], reveal_log, buffer_bytes);
        }
    }
}


void ReedSolomonDecodeHighPlanned(
    const backend::Ops& ops,
    uint64_t buffer_bytes,
    unsigned n,
    unsigned t,
    const void* const* coordinate_data,
    const uint16_t* block_input_counts,
    const uint32_t* requested_coordinates,
    const leopard2_internal::DecodeOutputBlock* output_blocks,
    unsigned output_block_count,
    const ffe_t* locator_logs,
    const ffe_t* output_factors,
    void** work)
{
    ReedSolomonDecodeHighPrunedPlanned(
        ops, buffer_bytes, n, t, coordinate_data, block_input_counts,
        requested_coordinates, output_blocks, output_block_count,
        locator_logs, output_factors, NULL, 0, NULL, 0, NULL, work);
}


void ReedSolomonDecodeHighPlanned(
    uint64_t buffer_bytes,
    unsigned n,
    unsigned t,
    const void* const* coordinate_data,
    const uint16_t* block_input_counts,
    const uint32_t* requested_coordinates,
    const leopard2_internal::DecodeOutputBlock* output_blocks,
    unsigned output_block_count,
    const ffe_t* locator_logs,
    const ffe_t* output_factors,
    void** work)
{
    ReedSolomonDecodeHighPlanned(backend::GetDefaultOps(), buffer_bytes, n, t,
        coordinate_data, block_input_counts, requested_coordinates,
        output_blocks, output_block_count, locator_logs, output_factors,
        work);
}


void ReedSolomonDecodeHighTiledPrunedPlanned(
    const backend::Ops& ops,
    uint64_t buffer_bytes,
    unsigned n,
    unsigned t,
    const void* const* coordinate_data,
    const uint16_t* block_input_counts,
    const uint32_t* requested_coordinates,
    const leopard2_internal::DecodeOutputBlock* output_blocks,
    unsigned output_block_count,
    const ffe_t* locator_logs,
    const ffe_t* output_factors,
    void* const* requested_output,
    const leopard2_internal::PrunedTransformBlock* input_plans,
    unsigned input_plan_count,
    const leopard2_internal::PrunedTransformPlan* output_plans,
    unsigned output_plan_count,
    void** work)
{
    LEO_DEBUG_ASSERT(t >= 2 && t < n && n <= kOrder);
    LEO_DEBUG_ASSERT(input_plan_count == 0 || input_plans != NULL);
    LEO_DEBUG_ASSERT(output_plan_count == 0 || output_plans != NULL);
    LEO_DEBUG_ASSERT(output_plan_count == 0 ||
        output_plan_count == output_block_count);
    LEO_DEBUG_ASSERT(n % t == 0);

    void** const accumulator = work;
    void** const tile = work + t;
    const unsigned block_count = n / t;

    unsigned input_plan_index = 0;
    LEO_DEBUG_ASSERT(block_input_counts[0] <= t);
    if (input_plan_index < input_plan_count && input_plans &&
        input_plans[input_plan_index].block == 0)
    {
        LEO_DEBUG_ASSERT(input_plans[input_plan_index].plan.size == t &&
            input_plans[input_plan_index].plan.shift == 0 &&
            input_plans[input_plan_index].plan.inverse);
        ++input_plan_index;
    }

    bool have_later_contribution = false;
    for (unsigned block = 1; block < block_count; ++block)
    {
        const unsigned input_count = block_input_counts[block];
        LEO_DEBUG_ASSERT(input_count <= t);
        if (input_count == 0)
            continue;
        const unsigned offset = block * t;
        const leopard2_internal::PrunedTransformPlan* pruned = NULL;
        if (input_plan_index < input_plan_count && input_plans &&
            input_plans[input_plan_index].block == block)
        {
            pruned = &input_plans[input_plan_index].plan;
            ++input_plan_index;
        }
        void** const block_work = have_later_contribution
            ? tile : accumulator;
        StageHighDecodeSources(
            buffer_bytes, coordinate_data + offset, block_work, t);
        if (pruned)
        {
            LEO_DEBUG_ASSERT(pruned->size == t &&
                pruned->shift == offset && pruned->inverse);
            const bool use_accumulating_sink =
                ShouldUsePrunedHighSyndromeSink(ops, buffer_bytes);
            const bool accumulated = have_later_contribution &&
                use_accumulating_sink &&
                leopard2_internal::
                    ExecutePrunedInverseTransformPlanAccumulate(
                        ops, buffer_bytes, *pruned,
                        block_work, accumulator);
            if (accumulated)
            {
#if defined(LEO2_ENABLE_TEST_HOOKS)
                TestHighSyndromeAccumulatedBlocks.fetch_add(
                    1, std::memory_order_relaxed);
                TestHighSyndromePrunedAccumulatedBlocks.fetch_add(
                    1, std::memory_order_relaxed);
#endif
                have_later_contribution = true;
                continue;
            }
#if defined(LEO2_ENABLE_TEST_HOOKS)
            if (have_later_contribution && use_accumulating_sink)
                TestHighSyndromePrunedFallbackBlocks.fetch_add(
                    1, std::memory_order_relaxed);
#endif
            const bool executed =
                leopard2_internal::ExecutePrunedTransformPlan(
                    ops, buffer_bytes, *pruned, block_work);
            LEO_DEBUG_ASSERT(executed);
            (void)executed;
#if defined(LEO2_ENABLE_TEST_HOOKS)
            TestHighSyndromeMaterializedBlocks.fetch_add(
                1, std::memory_order_relaxed);
#endif
        }
        else
        {
            bool accumulated = false;
            if (have_later_contribution)
                accumulated = IFFT_DIT_DecoderAccumulate(
                    ops, buffer_bytes, input_count,
                    block_work, accumulator, t,
                    FFTSkewStorage + offset);
            else
                IFFT_DIT_Decoder(
                    ops, buffer_bytes, input_count, block_work, t,
                    FFTSkewStorage + offset);
#if defined(LEO2_ENABLE_TEST_HOOKS)
            (have_later_contribution && accumulated
                ? TestHighSyndromeAccumulatedBlocks
                : TestHighSyndromeMaterializedBlocks)
                .fetch_add(1, std::memory_order_relaxed);
#else
            (void)accumulated;
#endif
        }
        if (have_later_contribution && pruned)
        {
            if (t < 8)
                VectorXOR(
                    ops, buffer_bytes, t, accumulator, block_work);
            else
                VectorXOR_Threads(
                    ops, buffer_bytes, t, accumulator, block_work);
        }
        have_later_contribution = true;
    }
    LEO_DEBUG_ASSERT(input_plan_index == input_plan_count);

    FinishHighDecodeSyndrome(
        ops, buffer_bytes, coordinate_data, block_input_counts[0],
        have_later_contribution, t, accumulator);
LEO_OPENMP_PARALLEL_FOR
    for (int i = 0; i < (int)t; ++i)
    {
        if (coordinate_data[i])
            mul_mem_inplace(
                ops, accumulator[i], locator_logs[i], buffer_bytes);
        else
            memset(accumulator[i], 0, buffer_bytes);
    }
    IFFT_DIT_Decoder(
        ops, buffer_bytes, t, accumulator, t, FFTSkewStorage);

    for (unsigned output_block = 0;
         output_block < output_block_count;
         ++output_block)
    {
        const leopard2_internal::DecodeOutputBlock& descriptor =
            output_blocks[output_block];
        LEO_DEBUG_ASSERT(descriptor.block > 0 && descriptor.block < block_count);
        LEO_DEBUG_ASSERT(descriptor.requested_prefix > 0 &&
            descriptor.requested_prefix <= t);
        LEO_DEBUG_ASSERT(descriptor.requested_begin < descriptor.requested_end);
        const unsigned offset = descriptor.block * t;
#if defined(LEO2_ENABLE_TEST_HOOKS)
        TestHighOutputBlocks.fetch_add(1, std::memory_order_relaxed);
#endif
        const leopard2_internal::PrunedTransformPlan* pruned =
            output_plans && output_plan_count == output_block_count &&
            output_plans[output_block].size != 0
                ? &output_plans[output_block]
                : NULL;
        if (pruned)
        {
#if defined(LEO2_ENABLE_TEST_HOOKS)
            TestHighPrunedOutputBlocks.fetch_add(
                1, std::memory_order_relaxed);
#endif
            LEO_DEBUG_ASSERT(pruned->size == t &&
                pruned->shift == offset && !pruned->inverse);
#if defined(LEO2_ENABLE_TEST_HOOKS)
            const bool force_copy_fallback =
                TestForceHighDecodeCopyFallback.load(
                std::memory_order_relaxed);
            if (!force_copy_fallback)
                TestHighFFTButterfly2OutCalls.fetch_add(
                    t / 2, std::memory_order_relaxed);
            const bool executed = !force_copy_fallback &&
                leopard2_internal::ExecutePrunedForwardTransformPlanFromSources(
                    ops, buffer_bytes, *pruned, accumulator, tile);
            LEO_DEBUG_ASSERT(executed || force_copy_fallback);
#else
            const bool executed =
                leopard2_internal::ExecutePrunedForwardTransformPlanFromSources(
                    ops, buffer_bytes, *pruned, accumulator, tile);
            LEO_DEBUG_ASSERT(executed);
#endif
            if (!executed)
            {
#if defined(LEO2_ENABLE_TEST_HOOKS)
                TestHighCompatibilityCopyFallbacks.fetch_add(
                    1, std::memory_order_relaxed);
#endif
LEO_OPENMP_PARALLEL_FOR
                for (int i = 0; i < (int)t; ++i)
                    memcpy(tile[i], accumulator[i], buffer_bytes);
                const bool fallback =
                    leopard2_internal::ExecutePrunedTransformPlan(
                        ops, buffer_bytes, *pruned, tile);
                LEO_DEBUG_ASSERT(fallback);
                (void)fallback;
            }
        }
        else
            FFT_DIT_FromCoefficients(
                ops, buffer_bytes, accumulator, tile,
                descriptor.requested_prefix, t,
                FFTSkewStorage + offset, NULL, SourceEvaluationHighDecode);
LEO_OPENMP_PARALLEL_FOR
        for (int i = (int)descriptor.requested_begin;
             i < (int)descriptor.requested_end;
             ++i)
        {
            const uint32_t coordinate = requested_coordinates[i];
            LEO_DEBUG_ASSERT(coordinate >= offset && coordinate < offset + t);
            LEO_DEBUG_ASSERT(requested_output[i] != NULL);
            const ffe_t reveal_log = SubMod(
                output_factors[coordinate], locator_logs[coordinate]);
            mul_mem(
                ops, requested_output[i], tile[coordinate - offset],
                reveal_log, buffer_bytes);
        }
    }
}


void ReedSolomonDecodeHighTiledPlanned(
    const backend::Ops& ops,
    uint64_t buffer_bytes,
    unsigned n,
    unsigned t,
    const void* const* coordinate_data,
    const uint16_t* block_input_counts,
    const uint32_t* requested_coordinates,
    const leopard2_internal::DecodeOutputBlock* output_blocks,
    unsigned output_block_count,
    const ffe_t* locator_logs,
    const ffe_t* output_factors,
    void* const* requested_output,
    void** work)
{
    ReedSolomonDecodeHighTiledPrunedPlanned(
        ops, buffer_bytes, n, t, coordinate_data, block_input_counts,
        requested_coordinates, output_blocks, output_block_count,
        locator_logs, output_factors, requested_output,
        NULL, 0, NULL, 0, work);
}

void ReedSolomonDecode(
    const backend::Ops& ops,
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

    FWHT<true>(error_locations, kOrder, m + original_count);

LEO_OPENMP_PARALLEL_FOR
    for (int i = 0; i < (int)kOrder; ++i)
        error_locations[i] = ((unsigned)error_locations[i] * (unsigned)LogWalsh[i]) % kModulus;

    FWHT<true>(error_locations, kOrder, kOrder);

    // work <- recovery data

LEO_OPENMP_PARALLEL_FOR
    for (int i = 0; i < (int)recovery_count; ++i)
    {
        if (recovery[i])
            mul_mem(ops, work[i], recovery[i], error_locations[i],
                buffer_bytes);
        else
            memset(work[i], 0, buffer_bytes);
    }
LEO_OPENMP_PARALLEL_FOR
    for (int i = recovery_count; i < (int)m; ++i)
        memset(work[i], 0, buffer_bytes);

    // work <- original data

LEO_OPENMP_PARALLEL_FOR
    for (int i = 0; i < (int)original_count; ++i)
    {
        if (original[i])
            mul_mem(ops, work[m + i], original[i],
                error_locations[m + i], buffer_bytes);
        else
            memset(work[m + i], 0, buffer_bytes);
    }
LEO_OPENMP_PARALLEL_FOR
    for (int i = m + original_count; i < (int)n; ++i)
        memset(work[i], 0, buffer_bytes);

    // work <- IFFT(work, n, 0)

    IFFT_DIT_Decoder(
        ops,
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
                ops,
                buffer_bytes,
                width,
                work + i - width,
                work + i);
        }
        else
        {
            VectorXOR_Threads(
                ops,
                buffer_bytes,
                width,
                work + i - width,
                work + i);
        }
    }

    // work <- FFT(work, n, 0) truncated to m + original_count

    const unsigned output_count = m + original_count;

#ifdef LEO_ERROR_BITFIELD_OPT
    FFT_DIT_ErrorBits(ops, buffer_bytes, work, output_count, n,
        FFTSkewStorage, error_bits);
#else
    FFT_DIT(ops, buffer_bytes, work, output_count, n, FFTSkewStorage);
#endif

    // Reveal erasures

    for (unsigned i = 0; i < original_count; ++i)
        if (!original[i])
            mul_mem(ops, work[i], work[i + m],
                kModulus - error_locations[i + m], buffer_bytes);
}


void ReedSolomonDecode(
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned m,
    unsigned n,
    const void* const* original,
    const void* const* recovery,
    void** work)
{
    ReedSolomonDecode(backend::GetDefaultOps(), buffer_bytes, original_count,
        recovery_count, m, n, original, recovery, work);
}


//------------------------------------------------------------------------------
// API

static bool IsInitialized = false;

bool Initialize()
{
    if (IsInitialized)
        return true;

    InitializeLogarithmTables();
    if (!InitializeMultiplyTables())
        return false;
    FFTInitialize();

    IsInitialized = true;
    return true;
}


}} // namespace leopard::ff16

#endif // LEO_HAS_FF16
