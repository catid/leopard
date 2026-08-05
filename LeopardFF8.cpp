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
#include "Leopard2Backend.h"

#ifdef LEO_HAS_FF8

#include <atomic>
#include <string.h>

#ifndef LEO2_EXPERIMENT_HIGH_HALF_TAIL_COLUMN
#define LEO2_EXPERIMENT_HIGH_HALF_TAIL_COLUMN 1
#endif

/* Same-text attribution control for the exact-64-byte T=32 final-IFFT range.
   The macro changes only a nonzero initialized data word. */
#ifndef LEO2_DIAGNOSTIC_DISABLE_T32_FINAL_IFFT2_RANGE
#define LEO2_DIAGNOSTIC_DISABLE_T32_FINAL_IFFT2_RANGE 0
#endif
#if LEO2_DIAGNOSTIC_DISABLE_T32_FINAL_IFFT2_RANGE < 0 || \
    LEO2_DIAGNOSTIC_DISABLE_T32_FINAL_IFFT2_RANGE > 1
#error "LEO2_DIAGNOSTIC_DISABLE_T32_FINAL_IFFT2_RANGE must be 0 or 1"
#endif

/*
    Same-text control for direct final-layer output from a dense partial P=16
    LOW_V1 parity block.  Both values are nonzero so candidate/control builds
    retain the same branch and transform code; only initialized data changes.
*/
#ifndef LEO2_EXPERIMENT_LOW_P16_PARTIAL_DIRECT_OUTPUT
#define LEO2_EXPERIMENT_LOW_P16_PARTIAL_DIRECT_OUTPUT 1
#endif
#if LEO2_EXPERIMENT_LOW_P16_PARTIAL_DIRECT_OUTPUT < 0 || \
    LEO2_EXPERIMENT_LOW_P16_PARTIAL_DIRECT_OUTPUT > 1
#error "LEO2_EXPERIMENT_LOW_P16_PARTIAL_DIRECT_OUTPUT must be 0 or 1"
#endif

#if defined(LEO2_WEIGHTED_LOCATOR_BOUNDARY_ABBA)
#include <stdlib.h>
#endif

#ifdef _MSC_VER
    #pragma warning(disable: 4752) // found Intel(R) Advanced Vector Extensions; consider using /arch:AVX
#endif

namespace leopard { namespace ff8 {

static volatile uint32_t g_high_final_ifft2_range_mode =
    1U + LEO2_DIAGNOSTIC_DISABLE_T32_FINAL_IFFT2_RANGE;
static std::atomic<uint32_t> g_low_p16_partial_direct_output_mode(
    2U - LEO2_EXPERIMENT_LOW_P16_PARTIAL_DIRECT_OUTPUT);
static thread_local bool g_low_p16_partial_direct_output_route_selected = false;

static inline unsigned LowP16PartialDirectOutputMode()
{
    return g_low_p16_partial_direct_output_mode.load(
        std::memory_order_relaxed);
}

#if defined(LEO2_ENABLE_TEST_HOOKS)
static std::atomic<uint64_t> TestIFFTDIT4Calls(0);
static std::atomic<uint64_t> TestIFFTDIT4XorCalls(0);
static std::atomic<uint64_t> TestFFTDIT4Calls(0);
static std::atomic<uint64_t> TestLowFFTButterfly2OutCalls(0);
static std::atomic<uint64_t> TestLowFFTButterfly4OutCalls(0);
static std::atomic<uint64_t> TestLowFFTButterfly8OutCalls(0);
static std::atomic<uint64_t> TestLowDirectOutputBlocks(0);
static std::atomic<uint64_t> TestLowDirectPartialOutputBlocks(0);
static std::atomic<uint64_t> TestLowDirectOutputShards(0);
static std::atomic<uint64_t> TestLowAvoidedScatterBytes(0);
static std::atomic<uint64_t> TestLowScatterBytes(0);
static std::atomic<uint64_t> TestLowDirectFFTButterfly2OutCalls(0);
static std::atomic<uint64_t> TestLowDirectFFTButterfly4OutCalls(0);
static std::atomic<uint64_t> TestHighIFFTButterfly4OutCalls(0);
static std::atomic<uint64_t> TestHighInputCopyShards(0);
static std::atomic<uint64_t> TestHighForwardFusedCalls(0);
static std::atomic<uint64_t> TestHighWholeTransformCalls(0);
static std::atomic<uint64_t> TestHighK5R5PartialCalls(0);
static std::atomic<uint64_t> TestHighTwoBlockCalls(0);
static std::atomic<uint64_t> TestHighSmallTransformCalls(0);
static std::atomic<uint64_t> TestHighT2PackedCalls(0);
static std::atomic<uint64_t> TestHighT4PackedCalls(0);
static std::atomic<uint64_t> TestHighT8PackedCalls(0);
static std::atomic<uint64_t> TestHighT8K7B1024DirectCalls(0);
static std::atomic<uint64_t> TestHighT8K8B1024DirectCalls(0);
static std::atomic<uint64_t> TestHighBalancedB64PackedCalls(0);
static std::atomic<uint64_t> TestHighFinalIFFT2RangeCalls(0);
static std::atomic<uint64_t> TestHighTailColumnCalls(0);
static std::atomic<uint64_t> TestHighHalfTailColumnCalls(0);
static std::atomic<uint64_t> TestHighK9R5TailCalls(0);
static std::atomic<uint64_t> TestHighK9R6R8TailCalls(0);
static std::atomic<uint64_t> TestHighOutputBlocks(0);
static std::atomic<uint64_t> TestHighFFTButterfly2OutCalls(0);
static std::atomic<uint64_t> TestHighFFTButterfly4OutCalls(0);
static std::atomic<uint64_t> TestHighCompatibilityCopyFallbacks(0);
static std::atomic<uint64_t> TestHighPrunedOutputBlocks(0);
static std::atomic<uint64_t> TestHighMatureOutputBlocks(0);
static std::atomic<bool> TestDisableHighFinalIFFT2Range(false);
static std::atomic<bool> TestForceHighDecodeCopyFallback(false);
static std::atomic<uint64_t> TestHighSyndromeAccumulatedBlocks(0);
static std::atomic<uint64_t> TestHighSyndromeMaterializedBlocks(0);
static std::atomic<uint64_t> TestHighSyndromePrunedAccumulatedBlocks(0);
static std::atomic<uint64_t> TestHighSyndromePrunedFallbackBlocks(0);
static std::atomic<uint64_t> TestHighSyndromeBlockZeroIFFTElisions(0);
static std::atomic<uint64_t> TestHighSyndromeForwardTransforms(0);
static std::atomic<uint64_t> TestHighSyndromeForwardTransformElisions(0);
static std::atomic<uint64_t> TestHighSyndromeBlockZeroXorShards(0);
static std::atomic<uint64_t> TestHighReceiveIFFTButterfly4OutCalls(0);
static std::atomic<uint64_t> TestHighReceiveCopyShards(0);
static std::atomic<uint64_t> TestHighReceiveZeroShards(0);
static std::atomic<uint64_t> TestHighLocatorWeightedIFFTButterfly4Calls(0);
static std::atomic<uint64_t> TestHighLocatorScaleRowsElided(0);
static std::atomic<uint64_t> TestHighLocatorInactiveRows(0);
static std::atomic<uint64_t> TestHighLocatorFallbackRows(0);
static std::atomic<uint64_t> TestSparseExactBlocks(0);
static std::atomic<uint64_t> TestSparsePrefixButterflies(0);
static std::atomic<uint64_t> TestSparseRetainedButterflies(0);
static std::atomic<uint64_t> TestSparseRequestedOutputCopies(0);
#endif

static uint16_t PrunedMultiplierLog(
    const void* context,
    uint32_t storage_index);

static void mul_mem_inplace(
    const backend::Ops& ops,
    void* data,
    ffe_t log_m,
    uint64_t bytes);

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
// For legacy-high K=T+d,R=T codes with d in {1,2}, the first T message shards
// form one complete inverse-transform block.  The remaining source(s) are
// systematic generator columns in parity-coordinate space.  Entries for each
// supported T occupy [T,2T), so the dyadic ranges are disjoint.
static ffe_t HighTailGeneratorLogs[2][kOrder / 2];
#if LEO2_EXPERIMENT_HIGH_HALF_TAIL_COLUMN
static uint16_t HighHalfTailGeneratorLogs[kOrder / 2];
#endif


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

#if defined(LEO_TRY_SSSE3)
struct Multiply128LUT_t
{
    LEO_M128 Value[2];
};

static const Multiply128LUT_t* Multiply128LUT = nullptr;

// 128-bit x_reg ^= y_reg * log_m
#define LEO_MULADD_128(x_reg, y_reg, table_lo, table_hi) { \
                LEO_M128 lo = _mm_and_si128(y_reg, clr_mask); \
                lo = _mm_shuffle_epi8(table_lo, lo); \
                LEO_M128 hi = _mm_srli_epi64(y_reg, 4); \
                hi = _mm_and_si128(hi, clr_mask); \
                hi = _mm_shuffle_epi8(table_hi, hi); \
                x_reg = _mm_xor_si128(x_reg, _mm_xor_si128(lo, hi)); }

#if defined(LEO_TRY_AVX2)

struct Multiply256LUT_t
{
    LEO_M256 Value[2];
};

static const Multiply256LUT_t* Multiply256LUT = nullptr;

// 256-bit x_reg ^= y_reg * log_m
#define LEO_MULADD_256(x_reg, y_reg, table_lo, table_hi) { \
                LEO_M256 lo = _mm256_and_si256(y_reg, clr_mask); \
                lo = _mm256_shuffle_epi8(table_lo, lo); \
                LEO_M256 hi = _mm256_srli_epi64(y_reg, 4); \
                hi = _mm256_and_si256(hi, clr_mask); \
                hi = _mm256_shuffle_epi8(table_hi, hi); \
                x_reg = _mm256_xor_si256(x_reg, _mm256_xor_si256(lo, hi)); }

#endif // LEO_TRY_AVX2
#endif // LEO_TRY_SSSE3

static bool InitializeMultiplyTables()
{
#if !defined(LEO_TRY_SSSE3)
    // Portable CMake builds route scalar multiplication through the immutable
    // backend ops table.  The former 64-KiB legacy scalar table had no reader.
    return true;
#else
    // Whole-translation-unit diagnostic/legacy SIMD builds retain their
    // existing nibble tables.  A CPU without PSHUFB falls through to ops.
    if (!CpuHasSSSE3)
        return true;

    void* table = nullptr;
#ifdef LEO_TRY_AVX2
    if (CpuHasAVX2)
        table = SIMDSafeAllocate(sizeof(Multiply256LUT_t) * kOrder);
    else
#endif // LEO_TRY_AVX2
        table = SIMDSafeAllocate(sizeof(Multiply128LUT_t) * kOrder);
    if (!table)
        return false;
#ifdef LEO_TRY_AVX2
    if (CpuHasAVX2)
        Multiply256LUT = reinterpret_cast<const Multiply256LUT_t*>(table);
    else
#endif // LEO_TRY_AVX2
        Multiply128LUT = reinterpret_cast<const Multiply128LUT_t*>(table);

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

            const LEO_M128 *v_ptr = reinterpret_cast<const LEO_M128 *>(&lut[0]);
            const LEO_M128 value = _mm_loadu_si128(v_ptr);

            // Store in 128-bit wide table
#if defined(LEO_TRY_AVX2)
            if (!CpuHasAVX2)
#endif // LEO_TRY_AVX2
                _mm_storeu_si128((LEO_M128*)&Multiply128LUT[log_m].Value[i], value);

            // Store in 256-bit wide table
#if defined(LEO_TRY_AVX2)
            if (CpuHasAVX2)
            {
                _mm256_storeu_si256((LEO_M256*)&Multiply256LUT[log_m].Value[i],
                    _mm256_broadcastsi128_si256(value));
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

#if defined(LEO_TRY_SSSE3)
    if (&ops == &backend::GetDefaultOps() && CpuHasSSSE3)
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
#endif // LEO_TRY_SSSE3

    ops.ff8_multiply(x, y, log_m, bytes);
}

// x[] ^= y[] * log_m for complete SIMD tiles.
static void muladd_mem(
    const backend::Ops& ops,
    void * LEO_RESTRICT x, const void * LEO_RESTRICT y,
    ffe_t log_m, uint64_t bytes)
{
#if defined(LEO_TRY_AVX2)
    if (&ops == &backend::GetDefaultOps() && CpuHasAVX2)
    {
        const LEO_M256 table_lo_y = _mm256_loadu_si256(&Multiply256LUT[log_m].Value[0]);
        const LEO_M256 table_hi_y = _mm256_loadu_si256(&Multiply256LUT[log_m].Value[1]);
        const LEO_M256 clr_mask = _mm256_set1_epi8(0x0f);
        LEO_M256 * LEO_RESTRICT x32 = reinterpret_cast<LEO_M256 *>(x);
        const LEO_M256 * LEO_RESTRICT y32 = reinterpret_cast<const LEO_M256 *>(y);

        do
        {
            LEO_M256 x0 = _mm256_loadu_si256(x32);
            LEO_M256 y0 = _mm256_loadu_si256(y32);
            LEO_MULADD_256(x0, y0, table_lo_y, table_hi_y);
            _mm256_storeu_si256(x32, x0);
            LEO_M256 x1 = _mm256_loadu_si256(x32 + 1);
            LEO_M256 y1 = _mm256_loadu_si256(y32 + 1);
            LEO_MULADD_256(x1, y1, table_lo_y, table_hi_y);
            _mm256_storeu_si256(x32 + 1, x1);
            x32 += 2, y32 += 2;
            bytes -= 64;
        } while (bytes > 0);
        return;
    }
#endif // LEO_TRY_AVX2

#if defined(LEO_TRY_SSSE3)
    if (&ops == &backend::GetDefaultOps() && CpuHasSSSE3)
    {
        const LEO_M128 table_lo_y = _mm_loadu_si128(&Multiply128LUT[log_m].Value[0]);
        const LEO_M128 table_hi_y = _mm_loadu_si128(&Multiply128LUT[log_m].Value[1]);
        const LEO_M128 clr_mask = _mm_set1_epi8(0x0f);
        LEO_M128 * LEO_RESTRICT x16 = reinterpret_cast<LEO_M128 *>(x);
        const LEO_M128 * LEO_RESTRICT y16 = reinterpret_cast<const LEO_M128 *>(y);

        do
        {
            for (unsigned i = 0; i < 4; ++i)
            {
                LEO_M128 x_data = _mm_loadu_si128(x16 + i);
                const LEO_M128 y_data = _mm_loadu_si128(y16 + i);
                LEO_MULADD_128(x_data, y_data, table_lo_y, table_hi_y);
                _mm_storeu_si128(x16 + i, x_data);
            }
            x16 += 4, y16 += 4;
            bytes -= 64;
        } while (bytes > 0);
        return;
    }
#endif // LEO_TRY_SSSE3

    ops.ff8_multiply_add(x, y, log_m, bytes);
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


const ffe_t* ElementLogTable()
{
    return LogLUT;
}

static bool InitializeHighTailGeneratorLogs()
{
    // Legacy-high K=T+d,R=T for d in {1,2} uses parent [4T,3T], with
    // systematic coordinates T..4T-1 and tail public sources at x=2T+d-1.
    // Compute those sources' systematic Lagrange columns at parity
    // coordinates 0..T-1:
    //
    //   L_x(y) = Z(y) / ((y + x) Z'(x)),
    //   Z(t) = product_{s=T}^{4T-1} (t + s).
    //
    // Setup is byte-independent and runs once after the legacy GF8 logarithm
    // tables exist.  Keeping the derivation here avoids embedding a generated
    // field-specific constant list in the execution path.
    // T=128 would require a 512-coordinate parent and cannot exist in GF8.
    static const unsigned kTransformSizes[] = { 8, 16, 32, 64 };
    for (unsigned size_i = 0;
         size_i < sizeof(kTransformSizes) / sizeof(kTransformSizes[0]);
         ++size_i)
    {
        const unsigned transform_size = kTransformSizes[size_i];
        for (unsigned tail = 0; tail < 2; ++tail)
        {
            const unsigned source_coordinate =
                transform_size * 2U + tail;
            ffe_t derivative = 1;
            for (unsigned coordinate = transform_size;
                 coordinate < transform_size * 4U; ++coordinate)
            {
                if (coordinate != source_coordinate)
                {
                    derivative = MultiplyElements(derivative,
                        static_cast<ffe_t>(
                            coordinate ^ source_coordinate));
                }
            }
            if (derivative == 0)
                return false;

            for (unsigned parity = 0; parity < transform_size; ++parity)
            {
                ffe_t vanishing = 1;
                for (unsigned coordinate = transform_size;
                     coordinate < transform_size * 4U; ++coordinate)
                {
                    vanishing = MultiplyElements(vanishing,
                        static_cast<ffe_t>(parity ^ coordinate));
                }
                const ffe_t denominator = MultiplyElements(derivative,
                    static_cast<ffe_t>(parity ^ source_coordinate));
                if (vanishing == 0 || denominator == 0)
                    return false;
                const ffe_t coefficient = MultiplyElements(
                    vanishing, InverseElement(denominator));
                if (coefficient == 0)
                    return false;
                HighTailGeneratorLogs[tail][transform_size + parity] =
                    ElementLog(coefficient);
            }
        }
    }

#if LEO2_EXPERIMENT_HIGH_HALF_TAIL_COLUMN
    // Legacy-high K=65,R=65 has T=128 and parent [256,128].  The first
    // 64 public sources occupy systematic coordinates 128..191, while the
    // one source that activates the upper half is coordinate 192.  Retain
    // the regular half-block transform and precompute that source's
    // systematic Lagrange column at parity coordinates 0..127:
    //
    //   L_x(y) = Z(y) / ((y + x) Z'(x)),
    //   Z(t) = product_{s=128}^{255} (t + s).
    //
    // This is an execution choice only; it evaluates the exact same legacy
    // parent code and therefore must remain byte-identical.
    static const unsigned kHalfSide = 64;
    static const unsigned kTransformSize = 128;
    static const unsigned kSourceCoordinate =
        kTransformSize + kHalfSide;
    ffe_t derivative = 1;
    for (unsigned coordinate = kTransformSize;
         coordinate < kOrder; ++coordinate)
    {
        if (coordinate != kSourceCoordinate)
        {
            derivative = MultiplyElements(derivative,
                static_cast<ffe_t>(coordinate ^ kSourceCoordinate));
        }
    }
    if (derivative == 0)
        return false;
    for (unsigned parity = 0; parity < kTransformSize; ++parity)
    {
        ffe_t vanishing = 1;
        for (unsigned coordinate = kTransformSize;
             coordinate < kOrder; ++coordinate)
        {
            vanishing = MultiplyElements(vanishing,
                static_cast<ffe_t>(parity ^ coordinate));
        }
        const ffe_t denominator = MultiplyElements(derivative,
            static_cast<ffe_t>(parity ^ kSourceCoordinate));
        if (vanishing == 0 || denominator == 0)
            return false;
        const ffe_t coefficient = MultiplyElements(
            vanishing, InverseElement(denominator));
        if (coefficient == 0)
            return false;
        HighHalfTailGeneratorLogs[parity] = ElementLog(coefficient);
    }
#endif

    return true;
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
    {
#ifdef LEO_TARGET_MOBILE
        // The legacy mobile RefMul loop is XOR-based because its existing
        // callers provide zeroed work.  Preserve the assignment contract of
        // this direct-repair helper for arbitrary caller output buffers.
        memset(destination, 0, static_cast<size_t>(complete));
#endif
        mul_mem(ops, destination, source, multiplier_log, complete);
    }
    const uint64_t residual = byte_count - complete;
    if (residual != 0)
        ops.ff8_multiply(
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
    MultiplyBytes(
        backend::GetDefaultOps(), destination, source, multiplier_log, byte_count);
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
    if (residual != 0)
        ops.ff8_multiply_add(
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
    MultiplyAddBytes(
        backend::GetDefaultOps(), destination, source, multiplier_log, byte_count);
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

    FWHT(LogWalsh, kOrder, kOrder);

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
        FWHT(active, n, n);

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
            LEO_MULADD_256(x_data, y_data, table_lo_y, table_hi_y); \
            _mm256_storeu_si256(x_ptr, x_data); }

            LEO_IFFTB_256(x32 + 1, y32 + 1);
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
            LEO_MULADD_128(x_data, y_data, table_lo_y, table_hi_y); \
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
#endif // LEO_TRY_SSSE3

    ops.ff8_ifft_butterfly2(x, y, log_m, bytes);
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
#if defined(LEO2_ENABLE_TEST_HOOKS)
    TestIFFTDIT4Calls.fetch_add(1, std::memory_order_relaxed);
#endif
#ifdef LEO_INTERLEAVE_BUTTERFLY4_OPT

#if defined(LEO_TRY_AVX2)

    if (&ops == &backend::GetDefaultOps() && CpuHasAVX2)
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
            // First layer:
            LEO_M256 work0_reg = _mm256_loadu_si256(work0);
            LEO_M256 work1_reg = _mm256_loadu_si256(work1);

            work1_reg = _mm256_xor_si256(work0_reg, work1_reg);
            if (log_m01 != kModulus)
                LEO_MULADD_256(work0_reg, work1_reg, t01_lo, t01_hi);

            LEO_M256 work2_reg = _mm256_loadu_si256(work2);
            LEO_M256 work3_reg = _mm256_loadu_si256(work3);

            work3_reg = _mm256_xor_si256(work2_reg, work3_reg);
            if (log_m23 != kModulus)
                LEO_MULADD_256(work2_reg, work3_reg, t23_lo, t23_hi);

            // Second layer:
            work2_reg = _mm256_xor_si256(work0_reg, work2_reg);
            work3_reg = _mm256_xor_si256(work1_reg, work3_reg);
            if (log_m02 != kModulus)
            {
                LEO_MULADD_256(work0_reg, work2_reg, t02_lo, t02_hi);
                LEO_MULADD_256(work1_reg, work3_reg, t02_lo, t02_hi);
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

#if defined(LEO_TRY_SSSE3)
    if (&ops == &backend::GetDefaultOps() && CpuHasSSSE3)
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
            // First layer:
            LEO_M128 work0_reg = _mm_loadu_si128(work0);
            LEO_M128 work1_reg = _mm_loadu_si128(work1);

            work1_reg = _mm_xor_si128(work0_reg, work1_reg);
            if (log_m01 != kModulus)
                LEO_MULADD_128(work0_reg, work1_reg, t01_lo, t01_hi);

            LEO_M128 work2_reg = _mm_loadu_si128(work2);
            LEO_M128 work3_reg = _mm_loadu_si128(work3);

            work3_reg = _mm_xor_si128(work2_reg, work3_reg);
            if (log_m23 != kModulus)
                LEO_MULADD_128(work2_reg, work3_reg, t23_lo, t23_hi);

            // Second layer:
            work2_reg = _mm_xor_si128(work0_reg, work2_reg);
            work3_reg = _mm_xor_si128(work1_reg, work3_reg);
            if (log_m02 != kModulus)
            {
                LEO_MULADD_128(work0_reg, work2_reg, t02_lo, t02_hi);
                LEO_MULADD_128(work1_reg, work3_reg, t02_lo, t02_hi);
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
#endif // LEO_TRY_SSSE3

#endif // LEO_INTERLEAVE_BUTTERFLY4_OPT

    ops.ff8_ifft_butterfly4(
        work[0], work[dist], work[dist * 2], work[dist * 3],
        log_m01, log_m23, log_m02, bytes);
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
    // Non-x86 production builds can retain the legacy in-field SIMD path
    // around the scalar Ops table.  Keep that default-context route native;
    // x86 production field TUs disable these macros and use the isolated
    // backend range below.
    if (&ops == &backend::GetDefaultOps())
    {
        for (unsigned i = 0; i < dist; ++i)
            IFFT_DIT4(ops, bytes, work + i, dist,
                log_m01, log_m23, log_m02);
        return;
    }
#endif
#if defined(LEO2_ENABLE_TEST_HOOKS)
    TestIFFTDIT4Calls.fetch_add(dist, std::memory_order_relaxed);
#endif
    ops.ff8_ifft_butterfly4_range(
        work, dist, log_m01, log_m23, log_m02, bytes, true);
}


// {x_out, y_out} ^= IFFT_DIT2( {x_in, y_in} )
static void IFFT_DIT2_xor(
    const backend::Ops& ops,
    void * LEO_RESTRICT x_in, void * LEO_RESTRICT y_in,
    void * LEO_RESTRICT x_out, void * LEO_RESTRICT y_out,
    const ffe_t log_m, uint64_t bytes)
{
#if defined(LEO_TRY_AVX2)
    if (&ops == &backend::GetDefaultOps() && CpuHasAVX2)
    {
        const LEO_M256 table_lo_y = _mm256_loadu_si256(&Multiply256LUT[log_m].Value[0]);
        const LEO_M256 table_hi_y = _mm256_loadu_si256(&Multiply256LUT[log_m].Value[1]);

        const LEO_M256 clr_mask = _mm256_set1_epi8(0x0f);

        const LEO_M256 * LEO_RESTRICT x32_in = reinterpret_cast<const LEO_M256 *>(x_in);
        const LEO_M256 * LEO_RESTRICT y32_in = reinterpret_cast<const LEO_M256 *>(y_in);
        LEO_M256 * LEO_RESTRICT x32_out = reinterpret_cast<LEO_M256 *>(x_out);
        LEO_M256 * LEO_RESTRICT y32_out = reinterpret_cast<LEO_M256 *>(y_out);

        do
        {
#define LEO_IFFTB_256_XOR(x_ptr_in, y_ptr_in, x_ptr_out, y_ptr_out) { \
            LEO_M256 x_data_out = _mm256_loadu_si256(x_ptr_out); \
            LEO_M256 y_data_out = _mm256_loadu_si256(y_ptr_out); \
            LEO_M256 x_data_in = _mm256_loadu_si256(x_ptr_in); \
            LEO_M256 y_data_in = _mm256_loadu_si256(y_ptr_in); \
            y_data_in = _mm256_xor_si256(y_data_in, x_data_in); \
            y_data_out = _mm256_xor_si256(y_data_out, y_data_in); \
            _mm256_storeu_si256(y_ptr_out, y_data_out); \
            LEO_MULADD_256(x_data_in, y_data_in, table_lo_y, table_hi_y); \
            x_data_out = _mm256_xor_si256(x_data_out, x_data_in); \
            _mm256_storeu_si256(x_ptr_out, x_data_out); }

            LEO_IFFTB_256_XOR(x32_in + 1, y32_in + 1, x32_out + 1, y32_out + 1);
            LEO_IFFTB_256_XOR(x32_in, y32_in, x32_out, y32_out);
            y32_in += 2, x32_in += 2, y32_out += 2, x32_out += 2;

            bytes -= 64;
        } while (bytes > 0);

        return;
    }
#endif // LEO_TRY_AVX2

#if defined(LEO_TRY_SSSE3)
    if (&ops == &backend::GetDefaultOps() && CpuHasSSSE3)
    {
        const LEO_M128 table_lo_y = _mm_loadu_si128(&Multiply128LUT[log_m].Value[0]);
        const LEO_M128 table_hi_y = _mm_loadu_si128(&Multiply128LUT[log_m].Value[1]);

        const LEO_M128 clr_mask = _mm_set1_epi8(0x0f);

        const LEO_M128 * LEO_RESTRICT x16_in = reinterpret_cast<const LEO_M128 *>(x_in);
        const LEO_M128 * LEO_RESTRICT y16_in = reinterpret_cast<const LEO_M128 *>(y_in);
        LEO_M128 * LEO_RESTRICT x16_out = reinterpret_cast<LEO_M128 *>(x_out);
        LEO_M128 * LEO_RESTRICT y16_out = reinterpret_cast<LEO_M128 *>(y_out);

        do
        {
#define LEO_IFFTB_128_XOR(x_ptr_in, y_ptr_in, x_ptr_out, y_ptr_out) { \
            LEO_M128 x_data_out = _mm_loadu_si128(x_ptr_out); \
            LEO_M128 y_data_out = _mm_loadu_si128(y_ptr_out); \
            LEO_M128 x_data_in = _mm_loadu_si128(x_ptr_in); \
            LEO_M128 y_data_in = _mm_loadu_si128(y_ptr_in); \
            y_data_in = _mm_xor_si128(y_data_in, x_data_in); \
            y_data_out = _mm_xor_si128(y_data_out, y_data_in); \
            _mm_storeu_si128(y_ptr_out, y_data_out); \
            LEO_MULADD_128(x_data_in, y_data_in, table_lo_y, table_hi_y); \
            x_data_out = _mm_xor_si128(x_data_out, x_data_in); \
            _mm_storeu_si128(x_ptr_out, x_data_out); }

            LEO_IFFTB_128_XOR(x16_in + 3, y16_in + 3, x16_out + 3, y16_out + 3);
            LEO_IFFTB_128_XOR(x16_in + 2, y16_in + 2, x16_out + 2, y16_out + 2);
            LEO_IFFTB_128_XOR(x16_in + 1, y16_in + 1, x16_out + 1, y16_out + 1);
            LEO_IFFTB_128_XOR(x16_in, y16_in, x16_out, y16_out);
            y16_in += 4, x16_in += 4, y16_out += 4, x16_out += 4;

            bytes -= 64;
        } while (bytes > 0);

        return;
    }
#endif // LEO_TRY_SSSE3

    ops.ff8_ifft_butterfly2_xor(
        x_in, y_in, x_out, y_out, log_m, bytes);
}


#if defined(LEO_TRY_AVX2) || defined(LEO_TRY_SSSE3)
// Legacy in-field SIMD route for xor_result ^= IFFT_DIT4(work).  Portable
// builds use the context-selected accumulating range backend, including the
// distance-one T=4 case.
static void IFFT_DIT4_xor(
    const backend::Ops& ops,
    uint64_t bytes,
    void** work_in,
    void** xor_out,
    unsigned dist,
    const ffe_t log_m01,
    const ffe_t log_m23,
    const ffe_t log_m02)
{
#if defined(LEO2_ENABLE_TEST_HOOKS)
    TestIFFTDIT4XorCalls.fetch_add(1, std::memory_order_relaxed);
#endif
#ifdef LEO_INTERLEAVE_BUTTERFLY4_OPT

#if defined(LEO_TRY_AVX2)

    if (&ops == &backend::GetDefaultOps() && CpuHasAVX2)
    {
        const LEO_M256 t01_lo = _mm256_loadu_si256(&Multiply256LUT[log_m01].Value[0]);
        const LEO_M256 t01_hi = _mm256_loadu_si256(&Multiply256LUT[log_m01].Value[1]);
        const LEO_M256 t23_lo = _mm256_loadu_si256(&Multiply256LUT[log_m23].Value[0]);
        const LEO_M256 t23_hi = _mm256_loadu_si256(&Multiply256LUT[log_m23].Value[1]);
        const LEO_M256 t02_lo = _mm256_loadu_si256(&Multiply256LUT[log_m02].Value[0]);
        const LEO_M256 t02_hi = _mm256_loadu_si256(&Multiply256LUT[log_m02].Value[1]);

        const LEO_M256 clr_mask = _mm256_set1_epi8(0x0f);

        const LEO_M256 * LEO_RESTRICT work0 = reinterpret_cast<const LEO_M256 *>(work_in[0]);
        const LEO_M256 * LEO_RESTRICT work1 = reinterpret_cast<const LEO_M256 *>(work_in[dist]);
        const LEO_M256 * LEO_RESTRICT work2 = reinterpret_cast<const LEO_M256 *>(work_in[dist * 2]);
        const LEO_M256 * LEO_RESTRICT work3 = reinterpret_cast<const LEO_M256 *>(work_in[dist * 3]);
        LEO_M256 * LEO_RESTRICT xor0 = reinterpret_cast<LEO_M256 *>(xor_out[0]);
        LEO_M256 * LEO_RESTRICT xor1 = reinterpret_cast<LEO_M256 *>(xor_out[dist]);
        LEO_M256 * LEO_RESTRICT xor2 = reinterpret_cast<LEO_M256 *>(xor_out[dist * 2]);
        LEO_M256 * LEO_RESTRICT xor3 = reinterpret_cast<LEO_M256 *>(xor_out[dist * 3]);

        do
        {
            // First layer:
            LEO_M256 work0_reg = _mm256_loadu_si256(work0);
            LEO_M256 work1_reg = _mm256_loadu_si256(work1);
            work0++, work1++;

            work1_reg = _mm256_xor_si256(work0_reg, work1_reg);
            if (log_m01 != kModulus)
                LEO_MULADD_256(work0_reg, work1_reg, t01_lo, t01_hi);

            LEO_M256 work2_reg = _mm256_loadu_si256(work2);
            LEO_M256 work3_reg = _mm256_loadu_si256(work3);
            work2++, work3++;

            work3_reg = _mm256_xor_si256(work2_reg, work3_reg);
            if (log_m23 != kModulus)
                LEO_MULADD_256(work2_reg, work3_reg, t23_lo, t23_hi);

            // Second layer:
            work2_reg = _mm256_xor_si256(work0_reg, work2_reg);
            work3_reg = _mm256_xor_si256(work1_reg, work3_reg);
            if (log_m02 != kModulus)
            {
                LEO_MULADD_256(work0_reg, work2_reg, t02_lo, t02_hi);
                LEO_MULADD_256(work1_reg, work3_reg, t02_lo, t02_hi);
            }

            work0_reg = _mm256_xor_si256(work0_reg, _mm256_loadu_si256(xor0));
            work1_reg = _mm256_xor_si256(work1_reg, _mm256_loadu_si256(xor1));
            work2_reg = _mm256_xor_si256(work2_reg, _mm256_loadu_si256(xor2));
            work3_reg = _mm256_xor_si256(work3_reg, _mm256_loadu_si256(xor3));

            _mm256_storeu_si256(xor0, work0_reg);
            _mm256_storeu_si256(xor1, work1_reg);
            _mm256_storeu_si256(xor2, work2_reg);
            _mm256_storeu_si256(xor3, work3_reg);
            xor0++, xor1++, xor2++, xor3++;

            bytes -= 32;
        } while (bytes > 0);

        return;
    }

#endif // LEO_TRY_AVX2

#if defined(LEO_TRY_SSSE3)
    if (&ops == &backend::GetDefaultOps() && CpuHasSSSE3)
    {
        const LEO_M128 t01_lo = _mm_loadu_si128(&Multiply128LUT[log_m01].Value[0]);
        const LEO_M128 t01_hi = _mm_loadu_si128(&Multiply128LUT[log_m01].Value[1]);
        const LEO_M128 t23_lo = _mm_loadu_si128(&Multiply128LUT[log_m23].Value[0]);
        const LEO_M128 t23_hi = _mm_loadu_si128(&Multiply128LUT[log_m23].Value[1]);
        const LEO_M128 t02_lo = _mm_loadu_si128(&Multiply128LUT[log_m02].Value[0]);
        const LEO_M128 t02_hi = _mm_loadu_si128(&Multiply128LUT[log_m02].Value[1]);

        const LEO_M128 clr_mask = _mm_set1_epi8(0x0f);

        const LEO_M128 * LEO_RESTRICT work0 = reinterpret_cast<const LEO_M128 *>(work_in[0]);
        const LEO_M128 * LEO_RESTRICT work1 = reinterpret_cast<const LEO_M128 *>(work_in[dist]);
        const LEO_M128 * LEO_RESTRICT work2 = reinterpret_cast<const LEO_M128 *>(work_in[dist * 2]);
        const LEO_M128 * LEO_RESTRICT work3 = reinterpret_cast<const LEO_M128 *>(work_in[dist * 3]);
        LEO_M128 * LEO_RESTRICT xor0 = reinterpret_cast<LEO_M128 *>(xor_out[0]);
        LEO_M128 * LEO_RESTRICT xor1 = reinterpret_cast<LEO_M128 *>(xor_out[dist]);
        LEO_M128 * LEO_RESTRICT xor2 = reinterpret_cast<LEO_M128 *>(xor_out[dist * 2]);
        LEO_M128 * LEO_RESTRICT xor3 = reinterpret_cast<LEO_M128 *>(xor_out[dist * 3]);

        do
        {
            // First layer:
            LEO_M128 work0_reg = _mm_loadu_si128(work0);
            LEO_M128 work1_reg = _mm_loadu_si128(work1);
            work0++, work1++;

            work1_reg = _mm_xor_si128(work0_reg, work1_reg);
            if (log_m01 != kModulus)
                LEO_MULADD_128(work0_reg, work1_reg, t01_lo, t01_hi);

            LEO_M128 work2_reg = _mm_loadu_si128(work2);
            LEO_M128 work3_reg = _mm_loadu_si128(work3);
            work2++, work3++;

            work3_reg = _mm_xor_si128(work2_reg, work3_reg);
            if (log_m23 != kModulus)
                LEO_MULADD_128(work2_reg, work3_reg, t23_lo, t23_hi);

            // Second layer:
            work2_reg = _mm_xor_si128(work0_reg, work2_reg);
            work3_reg = _mm_xor_si128(work1_reg, work3_reg);
            if (log_m02 != kModulus)
            {
                LEO_MULADD_128(work0_reg, work2_reg, t02_lo, t02_hi);
                LEO_MULADD_128(work1_reg, work3_reg, t02_lo, t02_hi);
            }

            work0_reg = _mm_xor_si128(work0_reg, _mm_loadu_si128(xor0));
            work1_reg = _mm_xor_si128(work1_reg, _mm_loadu_si128(xor1));
            work2_reg = _mm_xor_si128(work2_reg, _mm_loadu_si128(xor2));
            work3_reg = _mm_xor_si128(work3_reg, _mm_loadu_si128(xor3));

            _mm_storeu_si128(xor0, work0_reg);
            _mm_storeu_si128(xor1, work1_reg);
            _mm_storeu_si128(xor2, work2_reg);
            _mm_storeu_si128(xor3, work3_reg);
            xor0++, xor1++, xor2++, xor3++;

            bytes -= 16;
        } while (bytes > 0);

        return;
    }
#endif // LEO_TRY_SSSE3

#endif // LEO_INTERLEAVE_BUTTERFLY4_OPT

    ops.ff8_ifft_butterfly4(
        work_in[0], work_in[dist], work_in[dist * 2], work_in[dist * 3],
        log_m01, log_m23, log_m02, bytes);
    ops.xor_memory4(
        xor_out[0], work_in[0],
        xor_out[dist], work_in[dist],
        xor_out[dist * 2], work_in[dist * 2],
        xor_out[dist * 3], work_in[dist * 3], bytes);
}
#endif

static void IFFT_DIT4_xor_Range(
    const backend::Ops& ops,
    uint64_t bytes,
    void** work,
    void** xor_output,
    unsigned dist,
    const ffe_t log_m01,
    const ffe_t log_m23,
    const ffe_t log_m02)
{
#if defined(LEO_TRY_AVX2) || defined(LEO_TRY_SSSE3)
    if (&ops == &backend::GetDefaultOps())
    {
        for (unsigned i = 0; i < dist; ++i)
            IFFT_DIT4_xor(ops, bytes, work + i, xor_output + i, dist,
                log_m01, log_m23, log_m02);
        return;
    }
#endif
#if defined(LEO2_ENABLE_TEST_HOOKS)
    TestIFFTDIT4XorCalls.fetch_add(dist, std::memory_order_relaxed);
#endif
    ops.ff8_ifft_butterfly4_xor_range(
        work, xor_output, dist,
        log_m01, log_m23, log_m02, bytes);
}


enum EncoderIFFTCallsite
{
    EncoderIFFTHigh,
    EncoderIFFTLow
};

// Unrolled IFFT for encoder
static void IFFT_DIT_Encoder(
    const backend::Ops& ops,
    const uint64_t bytes,
    const void* const* data,
    const unsigned m_truncated,
    void** work,
    void** xor_result,
    const unsigned m,
    const ffe_t* skewLUT,
    bool contiguous_zero_suffix,
    EncoderIFFTCallsite callsite,
    bool allow_t32_final_ifft_range)
{
    // Decimation in time: Unroll 2 layers at a time
    unsigned dist = 1, dist4 = 4;
    if (m > 4)
    {
        // The encoder block transform consumes caller data exactly once. Complete
        // four-shard groups feed the first two inverse layers directly and
        // are stored to work only after both butterflies.  Qualified AVX2
        // ragged groups similarly consume their live prefix directly; other
        // backends retain the compatibility copy into zero-padded work.  Any
        // complete suffix groups still must be initialized because later
        // stages combine them with the active prefix.
        // Whole-codec ABBA measurements promote the mask-specialized path
        // only where eliminating the four-row staging pass clears the 5%
        // gate.  At wider sides those four rows are too small a fraction of
        // the transform; keep the mature schedule there.
        const bool measured_ragged_shape =
            (m == 8 && (bytes >= 4U * 1024U || m_truncated == 1)) ||
            (m == 16 && bytes >= 64U * 1024U);
        const bool direct_ragged_group = measured_ragged_shape &&
            (ops.kind == LEO2_BACKEND_AVX2 ||
             ops.kind == LEO2_BACKEND_GFNI) &&
            ops.ff8_weighted_ifft_butterfly4 != NULL &&
            (m_truncated & 3U) != 0;
        const unsigned initialized_prefix = direct_ragged_group
            ? (m_truncated + 3U) & ~3U
            : m_truncated;
        bool return_after_radix8 = false;
        if (contiguous_zero_suffix && initialized_prefix < m)
        {
            const size_t row_bytes = static_cast<size_t>(bytes);
            const size_t suffix_rows = m - initialized_prefix;
            const bool aggregate_zero_safe =
                row_bytes <= SIZE_MAX / suffix_rows;
            LEO_DEBUG_ASSERT(aggregate_zero_safe);
            if (!aggregate_zero_safe)
            {
                for (unsigned i = initialized_prefix; i < m; ++i)
                    memset(work[i], 0, bytes);
            }
            else
            {
#if !defined(NDEBUG)
                /* The aggregate form requires rows from one backing
                   allocation, not merely numerically adjacent unrelated
                   objects.  Integer stepping keeps this diagnostic free of
                   cross-object pointer arithmetic. */
                uintptr_t expected =
                    reinterpret_cast<uintptr_t>(work[initialized_prefix]);
                for (unsigned i = initialized_prefix + 1; i < m; ++i)
                {
                    LEO_DEBUG_ASSERT(expected <= UINTPTR_MAX - row_bytes);
                    expected += row_bytes;
                    LEO_DEBUG_ASSERT(
                        reinterpret_cast<uintptr_t>(work[i]) == expected);
                }
                static_cast<void>(expected);
#endif
                memset(work[initialized_prefix], 0,
                    suffix_rows * row_bytes);
            }
        }
        else
        {
            for (unsigned i = initialized_prefix; i < m; ++i)
                memset(work[i], 0, bytes);
        }

        unsigned r = 0;
        /*
            When the backend publishes an out-of-place radix-eight group, stage
            three inverse layers per load/store round instead of two.  That
            kernel is bound by load/store throughput rather than arithmetic, so
            a third less traffic per transform layer is a direct win; see the
            contract in Leopard2Backend.h.

            Restricted to a whole number of eight-coordinate groups.  A ragged
            tail cannot use it, because the distance-four layer of a short final
            group would pair coordinates that belong to the zero-padded suffix
            rather than to the group, and those pairings are the following
            stage's work.
        */
        /*
            The window is exactly m == 32, which is where every measured cell
            gained: 1.084x at 1 KiB, 1.109x at 4 KiB, 1.094x at 16 KiB and
            1.279x at 64 KiB per shard for K=224, plus 1.133x and 1.146x for
            K=200 and K=160.

            It is not wider for two separate reasons.  At m == 8 a radix-eight
            group would consume every layer, so the accumulating final layer
            that folds a message block into xor_result would never run and every
            block after the first would be dropped; a randomized cross-build
            differential caught exactly that.  At m == 16 three layers leave
            only one, which displaces the fused dist4 == m accumulating
            radix-four with a plain radix-two sweep plus a separate xor, and
            that measured 0.88x to 0.94x.  At m == 128 the staging is a smaller
            share of the call and the result is marginal, 1.05x to 1.07x with
            one neutral cell, so it stays out until it has its own evidence.
        */
        if (ops.ff8_ifft_butterfly8_out != NULL && m == 32 &&
            (m_truncated & 7U) == 0 && m_truncated >= 8)
        {
            for (; r + 8 <= m_truncated; r += 8)
            {
                const void* inputs[8];
                void* outputs[8];
                for (unsigned lane = 0; lane < 8; ++lane)
                {
                    inputs[lane] = data[r + lane];
                    outputs[lane] = work[r + lane];
                }
                const uint8_t skews[7] = {
                    static_cast<uint8_t>(skewLUT[r + 1]),
                    static_cast<uint8_t>(skewLUT[r + 3]),
                    static_cast<uint8_t>(skewLUT[r + 5]),
                    static_cast<uint8_t>(skewLUT[r + 7]),
                    static_cast<uint8_t>(skewLUT[r + 2]),
                    static_cast<uint8_t>(skewLUT[r + 6]),
                    static_cast<uint8_t>(skewLUT[r + 4])
                };
                ops.ff8_ifft_butterfly8_out(inputs, outputs, skews, bytes);
            }
            dist = 8;
            dist4 = 32;
            return_after_radix8 = true;
        }
        for (; !return_after_radix8 && r + 4 <= m_truncated; r += 4)
        {
#if defined(LEO2_ENABLE_TEST_HOOKS)
            TestHighIFFTButterfly4OutCalls.fetch_add(
                1, std::memory_order_relaxed);
#endif
            ops.ff8_ifft_butterfly4_out(
                data[r], data[r + 1], data[r + 2], data[r + 3],
                work[r], work[r + 1], work[r + 2], work[r + 3],
                skewLUT[r + 1], skewLUT[r + 3], skewLUT[r + 2],
                bytes);
        }
        if (!return_after_radix8 && r < m_truncated)
        {
            const unsigned remaining = m_truncated - r;
            if (direct_ragged_group)
            {
                const void* inputs[4] = { NULL, NULL, NULL, NULL };
                for (unsigned lane = 0; lane < remaining; ++lane)
                    inputs[lane] = data[r + lane];
                ops.ff8_weighted_ifft_butterfly4(
                    inputs[0], inputs[1], inputs[2], inputs[3],
                    work[r], work[r + 1], work[r + 2], work[r + 3],
                    0, 0, 0, 0,
                    static_cast<uint8_t>((1U << remaining) - 1U),
                    skewLUT[r + 1], skewLUT[r + 3], skewLUT[r + 2],
                    bytes);
            }
            else
            {
#if defined(LEO2_ENABLE_TEST_HOOKS)
                TestHighInputCopyShards.fetch_add(
                    remaining, std::memory_order_relaxed);
#endif
                for (unsigned i = r; i < m_truncated; ++i)
                    memcpy(work[i], data[i], bytes);
                IFFT_DIT4_Range(
                    ops, bytes, work + r, 1,
                    skewLUT[r + 1], skewLUT[r + 3], skewLUT[r + 2]);
            }
        }
        if (!return_after_radix8)
        {
            dist = 4;
            dist4 = 16;
        }
    }
    else
    {
        // Preserve the final-stage accumulating path for m <= 4: staging
        // directly into work and then fusing IFFT with xor_result avoids an
        // extra write/read pair that an out-of-place source transform would
        // reintroduce for these tiny transforms.
#if defined(LEO2_ENABLE_TEST_HOOKS)
        TestHighInputCopyShards.fetch_add(
            m_truncated, std::memory_order_relaxed);
#endif
        for (unsigned i = 0; i < m_truncated; ++i)
            memcpy(work[i], data[i], bytes);
        for (unsigned i = m_truncated; i < m; ++i)
            memset(work[i], 0, bytes);
    }

    for (; dist4 <= m; dist = dist4, dist4 <<= 2)
    {
        // For each set of dist*4 elements:
        for (unsigned r = 0; r < m_truncated; r += dist4)
        {
            const unsigned i_end = r + dist;
            const ffe_t log_m01 = skewLUT[i_end];
            const ffe_t log_m02 = skewLUT[i_end + dist];
            const ffe_t log_m23 = skewLUT[i_end + dist * 2];

            if (dist4 == m && xor_result)
            {
                IFFT_DIT4_xor_Range(
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
        bool use_prepared_final_range =
            callsite == EncoderIFFTHigh && m == 32U && bytes == 64U &&
            allow_t32_final_ifft_range &&
            g_high_final_ifft2_range_mode == 1U &&
            log_m != kModulus &&
            ops.ff8_ifft_butterfly2_range != NULL;
#if defined(LEO2_ENABLE_TEST_HOOKS)
        use_prepared_final_range = use_prepared_final_range &&
            !TestDisableHighFinalIFFT2Range.load(std::memory_order_relaxed);
#endif

        if (xor_result)
        {
            if (log_m == kModulus)
            {
                for (unsigned i = 0; i < dist; ++i)
                    xor_mem_2to1(
                        ops, xor_result[i], work[i], work[i + dist], bytes);
            }
            else if (use_prepared_final_range)
            {
#if defined(LEO2_ENABLE_TEST_HOOKS)
                TestHighFinalIFFT2RangeCalls.fetch_add(
                    1, std::memory_order_relaxed);
#endif
                ops.ff8_ifft_butterfly2_range(
                    work, xor_result, dist, log_m, bytes);
            }
            else
            {
                for (unsigned i = 0; i < dist; ++i)
                {
                    IFFT_DIT2_xor(
                        ops,
                        work[i],
                        work[i + dist],
                        xor_result[i],
                        xor_result[i + dist],
                        log_m,
                        bytes);
                }
            }
        }
        else
        {
            if (log_m == kModulus)
                VectorXOR(ops, bytes, dist, work + dist, work);
            else if (use_prepared_final_range)
            {
#if defined(LEO2_ENABLE_TEST_HOOKS)
                TestHighFinalIFFT2RangeCalls.fetch_add(
                    1, std::memory_order_relaxed);
#endif
                ops.ff8_ifft_butterfly2_range(
                    work, NULL, dist, log_m, bytes);
            }
            else
            {
                for (unsigned i = 0; i < dist; ++i)
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
    }
}


// Basic no-frills version for decoder.  When xor_result is non-null, the
// final inverse layer accumulates directly into those disjoint buffers.  The
// caller must not observe work afterwards: the accumulating backend need not
// materialize the final transformed values back into work.
static void IFFT_DIT_DecoderImpl(
    const backend::Ops& ops,
    const uint64_t bytes,
    const unsigned m_truncated,
    void** work,
    const unsigned m,
    const ffe_t* skewLUT,
    void** xor_result,
    const unsigned initial_dist)
{
    // Decimation in time: Unroll 2 layers at a time
    unsigned dist = initial_dist, dist4 = initial_dist << 2;
    for (; dist4 <= m; dist = dist4, dist4 <<= 2)
    {
        // For each set of dist*4 elements:
        for (unsigned r = 0; r < m_truncated; r += dist4)
        {
            const unsigned i_end = r + dist;
            const ffe_t log_m01 = skewLUT[i_end];
            const ffe_t log_m02 = skewLUT[i_end + dist];
            const ffe_t log_m23 = skewLUT[i_end + dist * 2];

            if (dist4 == m && xor_result)
            {
                IFFT_DIT4_xor_Range(
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
    }

    // If there is one layer left:
    if (dist < m)
    {
        // Assuming that dist = m / 2
        LEO_DEBUG_ASSERT(dist * 2 == m);

        const ffe_t log_m = skewLUT[dist];

        if (xor_result)
        {
            if (log_m == kModulus)
            {
                for (unsigned i = 0; i < dist; ++i)
                {
                    xor_mem(ops, xor_result[i], work[i], bytes);
                    xor_mem_2to1(
                        ops, xor_result[i + dist],
                        work[i], work[i + dist], bytes);
                }
            }
            else
            {
                for (unsigned i = 0; i < dist; ++i)
                {
                    IFFT_DIT2_xor(
                        ops,
                        work[i],
                        work[i + dist],
                        xor_result[i],
                        xor_result[i + dist],
                        log_m,
                        bytes);
                }
            }
        }
        else if (log_m == kModulus)
            VectorXOR(ops, bytes, dist, work + dist, work);
        else
        {
            for (unsigned i = 0; i < dist; ++i)
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
    else if (xor_result && m == 1 && m_truncated != 0)
        xor_mem(ops, xor_result[0], work[0], bytes);
}


static void IFFT_DIT_Decoder(
    const backend::Ops& ops,
    const uint64_t bytes,
    const unsigned m_truncated,
    void** work,
    const unsigned m,
    const ffe_t* skewLUT)
{
    IFFT_DIT_DecoderImpl(
        ops, bytes, m_truncated, work, m, skewLUT, NULL, 1);
}


// Algorithm 5's h*Lambda boundary is immediately followed by an inverse LCH
// transform.  Fold the four independent locator weights (or masked zeros) into
// the first two inverse layers so the workspace is read and written once
// instead of once for scaling and again for the boundary butterfly.
static void IFFT_DIT_DecoderWeightedLocator(
    const backend::Ops& ops,
    const uint64_t bytes,
    const void* const* coordinate_data,
    const ffe_t* locator_logs,
    void** work,
    const unsigned m,
    const ffe_t* skewLUT)
{
#if defined(LEO2_WEIGHTED_LOCATOR_BOUNDARY_ABBA)
    const char* const diagnostic_selection =
        getenv("LEO2_WEIGHTED_LOCATOR_BOUNDARY");
    const bool select_weighted = diagnostic_selection != NULL &&
        diagnostic_selection[0] == '1' && diagnostic_selection[1] == '\0' &&
        ops.ff8_weighted_ifft_butterfly4 != NULL;
#else
    // Same-binary, isolated ABBA qualification found a durable win only for
    // the AVX2 GF8 boundary in this deliberately conservative region.  Keep
    // scalar/SSSE3, small redundancy sides, sparse receive sets, and very
    // small/large shard passes on the established scale-then-IFFT path.
    const bool statically_qualified =
        (ops.kind == LEO2_BACKEND_AVX2 ||
         ops.kind == LEO2_BACKEND_AVX512 ||
         ops.kind == LEO2_BACKEND_GFNI) &&
        ops.ff8_weighted_ifft_butterfly4 != NULL &&
        m >= 64 &&
        bytes >= 16U * 1024U && bytes <= 256U * 1024U;
    unsigned live_count = 0;
    if (statically_qualified)
    {
        const unsigned required_live = (m + 1U) / 2U;
        for (unsigned i = 0; i < m && live_count < required_live; ++i)
            live_count += coordinate_data[i] != NULL;
    }
    const bool select_weighted =
        statically_qualified && live_count >= (m + 1U) / 2U;
#endif
    if (m < 4 || !select_weighted)
    {
        for (unsigned i = 0; i < m; ++i)
        {
            if (coordinate_data[i])
                mul_mem_inplace(ops, work[i], locator_logs[i], bytes);
            else
                memset(work[i], 0, bytes);
        }
#if defined(LEO2_ENABLE_TEST_HOOKS)
        TestHighLocatorFallbackRows.fetch_add(
            m, std::memory_order_relaxed);
#endif
        IFFT_DIT_Decoder(ops, bytes, m, work, m, skewLUT);
        return;
    }

    LEO_DEBUG_ASSERT((m & 3U) == 0);
#if defined(LEO2_ENABLE_TEST_HOOKS)
    uint64_t inactive_rows = 0;
#endif
    for (unsigned r = 0; r < m; r += 4)
    {
        uint8_t live_mask = 0;
        for (unsigned lane = 0; lane < 4; ++lane)
        {
            if (coordinate_data[r + lane])
                live_mask |= static_cast<uint8_t>(1U << lane);
#if defined(LEO2_ENABLE_TEST_HOOKS)
            else
                ++inactive_rows;
#endif
        }
        ops.ff8_weighted_ifft_butterfly4(
            work[r], work[r + 1], work[r + 2], work[r + 3],
            work[r], work[r + 1], work[r + 2], work[r + 3],
            locator_logs[r], locator_logs[r + 1],
            locator_logs[r + 2], locator_logs[r + 3], live_mask,
            skewLUT[r + 1], skewLUT[r + 3], skewLUT[r + 2], bytes);
    }
#if defined(LEO2_ENABLE_TEST_HOOKS)
    TestHighLocatorWeightedIFFTButterfly4Calls.fetch_add(
        m / 4, std::memory_order_relaxed);
    TestHighLocatorScaleRowsElided.fetch_add(m, std::memory_order_relaxed);
    TestHighLocatorInactiveRows.fetch_add(
        inactive_rows, std::memory_order_relaxed);
#endif
    IFFT_DIT_DecoderImpl(
        ops, bytes, m, work, m, skewLUT, NULL, 4);
}


static void StageHighDecodeSources(
    const uint64_t bytes,
    const void* const* sources,
    void** work,
    const unsigned count)
{
#if defined(LEO2_ENABLE_TEST_HOOKS)
    uint64_t copy_count = 0;
#endif
    for (unsigned i = 0; i < count; ++i)
    {
        if (sources[i])
        {
            memcpy(work[i], sources[i], bytes);
#if defined(LEO2_ENABLE_TEST_HOOKS)
            ++copy_count;
#endif
        }
        else
            memset(work[i], 0, bytes);
    }
#if defined(LEO2_ENABLE_TEST_HOOKS)
    TestHighReceiveCopyShards.fetch_add(
        copy_count, std::memory_order_relaxed);
    TestHighReceiveZeroShards.fetch_add(
        count - copy_count, std::memory_order_relaxed);
#endif
}


// Complete live GF8 receive groups enter the first two inverse layers without
// a staging pass only on the isolated and qualified SSSE3/AVX2-family
// backends.
// Pruned blocks are staged separately so their exact sink schedule remains
// authoritative; scalar, NEON, and other unmeasured backends copy first.
static void IFFT_DIT_DecoderFromSources(
    const backend::Ops& ops,
    const uint64_t bytes,
    const unsigned m_truncated,
    const void* const* sources,
    void** work,
    void** xor_result,
    const unsigned m,
    const ffe_t* skewLUT)
{
    LEO_DEBUG_ASSERT(m_truncated <= m);
    LEO_DEBUG_ASSERT(m <= 4 || (m & 3U) == 0);
    const bool qualified_source_backend =
        ops.kind == LEO2_BACKEND_SSSE3 ||
        ops.kind == LEO2_BACKEND_AVX2 ||
        ops.kind == LEO2_BACKEND_AVX512 ||
        ops.kind == LEO2_BACKEND_GFNI;
    if (m <= 4 || !qualified_source_backend)
    {
        StageHighDecodeSources(bytes, sources, work, m);
        IFFT_DIT_DecoderImpl(
            ops, bytes, m_truncated, work, m, skewLUT, xor_result, 1);
        return;
    }

    for (unsigned r = 0; r < m; r += 4)
    {
        const bool live0 = sources[r] != NULL;
        const bool live1 = sources[r + 1] != NULL;
        const bool live2 = sources[r + 2] != NULL;
        const bool live3 = sources[r + 3] != NULL;
        if (live0 && live1 && live2 && live3)
        {
            ops.ff8_ifft_butterfly4_out(
                sources[r], sources[r + 1], sources[r + 2], sources[r + 3],
                work[r], work[r + 1], work[r + 2], work[r + 3],
                skewLUT[r + 1], skewLUT[r + 3], skewLUT[r + 2], bytes);
#if defined(LEO2_ENABLE_TEST_HOOKS)
            TestHighReceiveIFFTButterfly4OutCalls.fetch_add(
                1, std::memory_order_relaxed);
#endif
            continue;
        }

        const bool live[4] = { live0, live1, live2, live3 };
        unsigned live_count = 0;
        for (unsigned lane = 0; lane < 4; ++lane)
        {
            if (live[lane])
            {
                memcpy(work[r + lane], sources[r + lane], bytes);
                ++live_count;
            }
            else
                memset(work[r + lane], 0, bytes);
        }
#if defined(LEO2_ENABLE_TEST_HOOKS)
        TestHighReceiveCopyShards.fetch_add(
            live_count, std::memory_order_relaxed);
        TestHighReceiveZeroShards.fetch_add(
            4U - live_count, std::memory_order_relaxed);
#endif
        if (live_count != 0)
            IFFT_DIT4_Range(
                ops, bytes, work + r, 1,
                skewLUT[r + 1], skewLUT[r + 3], skewLUT[r + 2]);
    }

    IFFT_DIT_DecoderImpl(
        ops, bytes, m_truncated, work, m, skewLUT, xor_result, 4);
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
        const LEO_M256 table_lo_y = _mm256_loadu_si256(&Multiply256LUT[log_m].Value[0]);
        const LEO_M256 table_hi_y = _mm256_loadu_si256(&Multiply256LUT[log_m].Value[1]);

        const LEO_M256 clr_mask = _mm256_set1_epi8(0x0f);

        LEO_M256 * LEO_RESTRICT x32 = reinterpret_cast<LEO_M256 *>(x);
        LEO_M256 * LEO_RESTRICT y32 = reinterpret_cast<LEO_M256 *>(y);

        do
        {
#define LEO_FFTB_256(x_ptr, y_ptr) { \
            LEO_M256 y_data = _mm256_loadu_si256(y_ptr); \
            LEO_M256 x_data = _mm256_loadu_si256(x_ptr); \
            LEO_MULADD_256(x_data, y_data, table_lo_y, table_hi_y); \
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

#if defined(LEO_TRY_SSSE3)
    if (&ops == &backend::GetDefaultOps() && CpuHasSSSE3)
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
            LEO_M128 x_data = _mm_loadu_si128(x_ptr); \
            LEO_MULADD_128(x_data, y_data, table_lo_y, table_hi_y); \
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
#endif // LEO_TRY_SSSE3

    ops.ff8_fft_butterfly2(x, y, log_m, bytes);
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
#if defined(LEO2_ENABLE_TEST_HOOKS)
    TestFFTDIT4Calls.fetch_add(1, std::memory_order_relaxed);
#endif
#ifdef LEO_INTERLEAVE_BUTTERFLY4_OPT

#if defined(LEO_TRY_AVX2)
    if (&ops == &backend::GetDefaultOps() && CpuHasAVX2)
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
            LEO_M256 work0_reg = _mm256_loadu_si256(work0);
            LEO_M256 work2_reg = _mm256_loadu_si256(work2);
            LEO_M256 work1_reg = _mm256_loadu_si256(work1);
            LEO_M256 work3_reg = _mm256_loadu_si256(work3);

            // First layer:
            if (log_m02 != kModulus)
            {
                LEO_MULADD_256(work0_reg, work2_reg, t02_lo, t02_hi);
                LEO_MULADD_256(work1_reg, work3_reg, t02_lo, t02_hi);
            }
            work2_reg = _mm256_xor_si256(work0_reg, work2_reg);
            work3_reg = _mm256_xor_si256(work1_reg, work3_reg);

            // Second layer:
            if (log_m01 != kModulus)
                LEO_MULADD_256(work0_reg, work1_reg, t01_lo, t01_hi);
            work1_reg = _mm256_xor_si256(work0_reg, work1_reg);

            _mm256_storeu_si256(work0, work0_reg);
            _mm256_storeu_si256(work1, work1_reg);
            work0++, work1++;

            if (log_m23 != kModulus)
                LEO_MULADD_256(work2_reg, work3_reg, t23_lo, t23_hi);
            work3_reg = _mm256_xor_si256(work2_reg, work3_reg);

            _mm256_storeu_si256(work2, work2_reg);
            _mm256_storeu_si256(work3, work3_reg);
            work2++, work3++;

            bytes -= 32;
        } while (bytes > 0);

        return;
    }
#endif // LEO_TRY_AVX2

#if defined(LEO_TRY_SSSE3)
    if (&ops == &backend::GetDefaultOps() && CpuHasSSSE3)
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
            LEO_M128 work0_reg = _mm_loadu_si128(work0);
            LEO_M128 work2_reg = _mm_loadu_si128(work2);
            LEO_M128 work1_reg = _mm_loadu_si128(work1);
            LEO_M128 work3_reg = _mm_loadu_si128(work3);

            // First layer:
            if (log_m02 != kModulus)
            {
                LEO_MULADD_128(work0_reg, work2_reg, t02_lo, t02_hi);
                LEO_MULADD_128(work1_reg, work3_reg, t02_lo, t02_hi);
            }
            work2_reg = _mm_xor_si128(work0_reg, work2_reg);
            work3_reg = _mm_xor_si128(work1_reg, work3_reg);

            // Second layer:
            if (log_m01 != kModulus)
                LEO_MULADD_128(work0_reg, work1_reg, t01_lo, t01_hi);
            work1_reg = _mm_xor_si128(work0_reg, work1_reg);

            _mm_storeu_si128(work0, work0_reg);
            _mm_storeu_si128(work1, work1_reg);
            work0++, work1++;

            if (log_m23 != kModulus)
                LEO_MULADD_128(work2_reg, work3_reg, t23_lo, t23_hi);
            work3_reg = _mm_xor_si128(work2_reg, work3_reg);

            _mm_storeu_si128(work2, work2_reg);
            _mm_storeu_si128(work3, work3_reg);
            work2++, work3++;

            bytes -= 16;
        } while (bytes > 0);

        return;
    }
#endif // LEO_TRY_SSSE3

#endif // LEO_INTERLEAVE_BUTTERFLY4_OPT

    ops.ff8_fft_butterfly4(
        work[0], work[dist], work[dist * 2], work[dist * 3],
        log_m01, log_m23, log_m02, bytes);
}

static void FFT_DIT4_Range(
    const backend::Ops& ops,
    uint64_t bytes,
    void** work,
    unsigned dist,
    const ffe_t log_m01,
    const ffe_t log_m23,
    const ffe_t log_m02,
    bool prefer_fused = false)
{
    if (dist == 1 && !prefer_fused)
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
#if defined(LEO2_ENABLE_TEST_HOOKS)
    TestFFTDIT4Calls.fetch_add(dist, std::memory_order_relaxed);
#endif
    ops.ff8_fft_butterfly4_range(
        work, dist, log_m01, log_m23, log_m02, bytes, prefer_fused);
}


// In-place FFT for encoder and decoder
static void FFT_DIT(
    const backend::Ops& ops,
    const uint64_t bytes,
    void** work,
    const unsigned m_truncated,
    const unsigned m,
    const ffe_t* skewLUT,
    bool prefer_fused = false)
{
    // Decimation in time: Unroll 2 layers at a time
    unsigned dist4 = m, dist = m >> 2;
    for (; dist != 0; dist4 = dist, dist >>= 2)
    {
        // For each set of dist*4 elements:
        for (unsigned r = 0; r < m_truncated; r += dist4)
        {
            const unsigned i_end = r + dist;
            const ffe_t log_m01 = skewLUT[i_end];
            const ffe_t log_m02 = skewLUT[i_end + dist];
            const ffe_t log_m23 = skewLUT[i_end + dist * 2];

            FFT_DIT4_Range(
                ops, bytes, work + r, dist,
                log_m01, log_m23, log_m02, prefer_fused);
        }
    }

    // If there is one layer left:
    if (dist4 == 2)
    {
        for (unsigned r = 0; r < m_truncated; r += 2)
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

enum SourceEvaluationCallsite
{
    SourceEvaluationLowEncode,
    SourceEvaluationHighDecode
};

// Out-of-place first FFT layer(s) for evaluation from reusable coefficients.
// The coefficient pointers remain read-only across every parity or message
// coset.  The backend combines the former copy pass with the first butterfly
// pass.  Remaining layers operate in-place on evaluation_work, except that a
// dense full block, or a qualified aligned partial low-encode block, may
// write its final layer to caller output.
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
    LEO_DEBUG_ASSERT(!final_outputs || m_truncated == m ||
        (callsite == SourceEvaluationLowEncode &&
         (LowP16PartialDirectOutputMode() == 1U ||
          LowP16PartialDirectOutputMode() == 3U) &&
         ops.kind == LEO2_BACKEND_AVX2 && m == 16 &&
         (m_truncated == 4 || m_truncated == 8 || m_truncated == 12)));
#if defined(LEO2_ENABLE_TEST_HOOKS)
    if (callsite == SourceEvaluationHighDecode)
        TestHighMatureOutputBlocks.fetch_add(
            1, std::memory_order_relaxed);
    if (callsite == SourceEvaluationHighDecode &&
        TestForceHighDecodeCopyFallback.load(std::memory_order_relaxed))
    {
        for (unsigned i = 0; i < m; ++i)
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
        ops.ff8_fft_butterfly2_out(
            coefficients[0], coefficients[1],
            final_outputs ? final_outputs[0] : evaluation_work[0],
            final_outputs ? final_outputs[1] : evaluation_work[1],
            skewLUT[1], bytes);
        return;
    }

    /*
        A radix-eight first stage carries three forward layers out of place
        instead of two.  At m == 32 that is what removes the trailing radix-two
        sweep: five layers become one radix-eight round plus one radix-four
        round, where the radix-four schedule needs two radix-four rounds plus a
        radix-two.  A radix-two round costs 2.0 shard touches per shard-layer
        against radix-eight's 0.67, and a low-rate profile at K=32 R=192 64 KiB
        attributes 21% of encode to the in-place radix-two butterfly.

        Skew for a layer at distance D inside a group is
        skew[subgroup base + D], the same rule the radix-four stage uses, so
        with d = m >> 3 the seven skews are skew[4d]; skew[2d], skew[6d]; and
        skew[d], skew[3d], skew[5d], skew[7d].

        m == 16 is excluded because ceil(4/2) == ceil(4/3): it would save no
        touches while still displacing a fused radix-four.
    */
    /*
        A radix-eight first stage carries three forward layers out of place
        instead of two.  At m == 32 that removes the trailing radix-two sweep:
        five layers become one radix-eight round plus one radix-four round,
        where the radix-four schedule needs two radix-four rounds plus a
        radix-two.  A radix-two round costs 2.0 shard touches per shard-layer
        against radix-eight's 0.67, and a low-rate profile at K=32 R=192 64 KiB
        attributes 21% of encode to the in-place radix-two butterfly.

        Skew for a layer at distance D inside a group is skew[subgroup base + D],
        the rule the radix-four stage below uses, so with d = m >> 3 the seven
        skews are skew[4d]; skew[2d], skew[6d]; then skew[d], skew[3d],
        skew[5d], skew[7d].

        m == 16 is excluded: ceil(4/2) == ceil(4/3), so it would save no touches
        while still displacing a fused radix-four.
    */
    unsigned first_dist = m >> 2;
    if (ops.ff8_fft_butterfly8_out != NULL && m == 32)
    {
        const unsigned d = m >> 3;
        const uint8_t skews[7] = {
            static_cast<uint8_t>(skewLUT[d * 4]),
            static_cast<uint8_t>(skewLUT[d * 2]),
            static_cast<uint8_t>(skewLUT[d * 6]),
            static_cast<uint8_t>(skewLUT[d]),
            static_cast<uint8_t>(skewLUT[d * 3]),
            static_cast<uint8_t>(skewLUT[d * 5]),
            static_cast<uint8_t>(skewLUT[d * 7])
        };
        for (unsigned i = 0; i < d; ++i)
        {
#if defined(LEO2_ENABLE_TEST_HOOKS)
            if (callsite == SourceEvaluationLowEncode)
                TestLowFFTButterfly8OutCalls.fetch_add(
                    1, std::memory_order_relaxed);
#endif
            const void* inputs[8];
            void* outputs[8];
            for (unsigned lane = 0; lane < 8; ++lane)
            {
                inputs[lane] = coefficients[i + lane * d];
                outputs[lane] = evaluation_work[i + lane * d];
            }
            ops.ff8_fft_butterfly8_out(inputs, outputs, skews, bytes);
        }
        first_dist = d;
    }
    else
    {
    const ffe_t log_m01 = skewLUT[first_dist];
    const ffe_t log_m02 = skewLUT[first_dist * 2];
    const ffe_t log_m23 = skewLUT[first_dist * 3];
    for (unsigned i = 0; i < first_dist; ++i)
    {
#if defined(LEO2_ENABLE_TEST_HOOKS)
        (callsite == SourceEvaluationLowEncode
            ? TestLowFFTButterfly4OutCalls
            : TestHighFFTButterfly4OutCalls).fetch_add(
                1, std::memory_order_relaxed);
#endif
        void* const* const output = final_outputs && m == 4
            ? final_outputs + i : evaluation_work + i;
        ops.ff8_fft_butterfly4_out(
            coefficients[i], coefficients[i + first_dist],
            coefficients[i + first_dist * 2],
            coefficients[i + first_dist * 3],
            output[0], output[first_dist],
            output[first_dist * 2], output[first_dist * 3],
            log_m01, log_m23, log_m02, bytes);
    }
    }
    if (final_outputs && m == 4)
        return;

    unsigned dist4 = first_dist;
    unsigned dist = first_dist >> 2;
    for (; dist != 0; dist4 = dist, dist >>= 2)
    {
#if defined(LEO2_ENABLE_TEST_HOOKS)
        if (final_outputs && dist == 1 &&
            callsite == SourceEvaluationLowEncode)
            TestLowDirectFFTButterfly4OutCalls.fetch_add(
                m_truncated / 4, std::memory_order_relaxed);
#endif
        for (unsigned r = 0; r < m_truncated; r += dist4)
        {
            const unsigned i_end = r + dist;
            const ffe_t remaining_log_m01 = skewLUT[i_end];
            const ffe_t remaining_log_m02 = skewLUT[i_end + dist];
            const ffe_t remaining_log_m23 = skewLUT[i_end + dist * 2];
            if (final_outputs && dist == 1)
                ops.ff8_fft_butterfly4_out(
                    evaluation_work[r], evaluation_work[r + 1],
                    evaluation_work[r + 2], evaluation_work[r + 3],
                    final_outputs[r], final_outputs[r + 1],
                    final_outputs[r + 2], final_outputs[r + 3],
                    remaining_log_m01, remaining_log_m23,
                    remaining_log_m02, bytes);
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
        for (unsigned r = 0; r < m_truncated; r += 2)
        {
            const ffe_t log_m = skewLUT[r + 1];
            if (final_outputs)
                ops.ff8_fft_butterfly2_out(
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

#if LEO2_EXPERIMENT_HIGH_HALF_TAIL_COLUMN
#if defined(_MSC_VER)
#define LEO2_HIGH_HALF_TAIL_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define LEO2_HIGH_HALF_TAIL_NOINLINE __attribute__((noinline))
#else
#define LEO2_HIGH_HALF_TAIL_NOINLINE
#endif

static LEO2_HIGH_HALF_TAIL_NOINLINE void AccumulateHighHalfTailColumn(
    const backend::Ops& ops,
    void** work,
    unsigned output_count,
    const void* source,
    uint64_t buffer_bytes)
{
    unsigned output = 0;
    if (ops.ff8_multiply_add_outputs)
    {
        for (; output + 1 < output_count;)
        {
            const unsigned remaining = output_count - output;
            const unsigned group_count = remaining < 8 ? remaining : 8;
            void* destinations[8];
            uint16_t logs[8];
            for (unsigned lane = 0; lane < group_count; ++lane)
            {
                destinations[lane] = work[output + lane];
                logs[lane] = HighHalfTailGeneratorLogs[output + lane];
            }
            ops.ff8_multiply_add_outputs(
                destinations, source, logs, group_count, buffer_bytes);
            output += group_count;
        }
    }
    for (; output < output_count; ++output)
    {
        MultiplyAddBytes(ops, work[output], source,
            static_cast<ffe_t>(HighHalfTailGeneratorLogs[output]),
            buffer_bytes);
    }
}

#undef LEO2_HIGH_HALF_TAIL_NOINLINE
#endif

void ReedSolomonEncode(
    const backend::Ops& ops,
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned requested_output_count,
    unsigned m,
    const void* const* data,
    void** work,
    const leopard2_internal::SparseForwardPlanBatchView* sparse_plans,
    bool exact_public_64_byte_shard,
    bool allow_sub_2k_register_kernels,
    bool contiguous_temporary_work)
{
    LEO_DEBUG_ASSERT(!exact_public_64_byte_shard || buffer_bytes == 64U);
#if !defined(LEO2_ENABLE_TEST_HOOKS)
    (void)requested_output_count;
#endif
    const bool whole_transform_size =
        m == 8 || m == 16 || m == 32 || m == 64;
    const bool dense_schedule = !sparse_plans ||
        sparse_plans->block_count == 0;
    // Exact-main AVX2 crossover measurements qualify every valid T=2 shape
    // from 2 KiB.  Below that boundary the generated packed terminal handles
    // dense public and prevalidated layouts without perturbing this generic
    // transform.  For T=4, complete sub-2-KiB AVX2 vectors can use the
    // existing register-resident K=3..7/K=9..11 kernels without staging the
    // four-row accumulator.  Larger calls retain the established T=4
    // thresholds and K=8/K>=12 direct-input callback.  Keep the boundaries
    // explicit: buffer_bytes is the current execution pass, so a padded tail
    // cannot inherit its aligned prefix's decision.
    const bool sub_2k_register_t4 =
        allow_sub_2k_register_kernels &&
        ops.kind == LEO2_BACKEND_AVX2 &&
        m == 4 && buffer_bytes >= 32U && buffer_bytes < 2U * 1024U &&
        (buffer_bytes & 31U) == 0 &&
        ((original_count >= 3 && original_count <= 7) ||
         (original_count >= 9 && original_count <= 11));
    const bool small_transform_shape =
        m == 2 ? original_count >= 2 && buffer_bytes >= 2U * 1024U :
        m == 4 && original_count >= 3 &&
            (buffer_bytes >= 2U * 1024U || sub_2k_register_t4);
    if ((ops.kind == LEO2_BACKEND_AVX2 ||
         ops.kind == LEO2_BACKEND_GFNI) &&
        ops.ff8_high_encode_small && small_transform_shape &&
        requested_output_count == recovery_count &&
        dense_schedule)
    {
#if defined(LEO2_ENABLE_TEST_HOOKS)
        TestHighSmallTransformCalls.fetch_add(1, std::memory_order_relaxed);
        const bool fused_t4 = sub_2k_register_t4 ||
            (m == 4 &&
             (((original_count == 5 || original_count == 6 ||
                original_count >= 9) && buffer_bytes >= 2U * 1024U) ||
             ((original_count == 3 || original_count == 7) &&
                buffer_bytes >= 4U * 1024U) ||
             (original_count == 4 && buffer_bytes >= 8U * 1024U)) &&
            ((original_count >= 3 && original_count <= 7) ||
             (original_count >= 9 && original_count <= 11)) &&
            (buffer_bytes & 31U) == 0);
        TestHighInputCopyShards.fetch_add(
            fused_t4 || original_count <= m ||
                    (m == 4 && original_count % m == 3)
                ? 0U : original_count % m,
            std::memory_order_relaxed);
#endif
        ops.ff8_high_encode_small(
            data, original_count, work, m, FFTSkewStorage + m,
            FFTSkewStorage, buffer_bytes);
        return;
    }
#if LEO2_EXPERIMENT_HIGH_HALF_TAIL_COLUMN
    // Frozen same-source AVX2 ABBA screens selected the 8-KiB floor:
    // 1.073x at 8 KiB and 1.064x at 64 KiB.  The 4-KiB cell was only
    // 1.003x, so it deliberately retains the mature transform.  Every
    // adjacent K/R shape remained within 0.4% of its control.
    if (ops.kind == LEO2_BACKEND_AVX2 &&
        ops.ff8_multiply_add_outputs &&
        m == 128 && original_count == 65 &&
        recovery_count == 65 &&
        buffer_bytes >= 8192 &&
        requested_output_count == recovery_count)
    {
        // The 65th source is the only dynamic input in the upper half of this
        // T=128 systematic block.  Transform the regular 64-source component
        // and add its independently derived generator column in source-major
        // groups, loading the source once per group instead of activating the
        // complete upper-half inverse-transform dependency tree.
        ReedSolomonEncode(
            ops, buffer_bytes, 64, recovery_count, requested_output_count,
            m, data, work, NULL, exact_public_64_byte_shard);
#if defined(LEO2_ENABLE_TEST_HOOKS)
        TestHighHalfTailColumnCalls.fetch_add(
            1, std::memory_order_relaxed);
#endif
        AccumulateHighHalfTailColumn(
            ops, work, recovery_count, data[64], buffer_bytes);
        return;
    }
#endif
    const bool direct_tail_column_size =
        m == 8 || m == 16 || m == 32 || m == 64;
    const bool direct_one_tail =
        direct_tail_column_size && original_count == m + 1U;
    const bool direct_two_tails =
        m == 32U && original_count == m + 2U && buffer_bytes >= 16384U;
    if ((ops.kind == LEO2_BACKEND_AVX2 ||
         ops.kind == LEO2_BACKEND_GFNI) &&
        (direct_one_tail || direct_two_tails) &&
        recovery_count == m &&
        requested_output_count == m)
    {
        // Encode the complete first message block through the mature
        // transform, then add the precomputed systematic generator column for
        // message coordinates 2T and (when selected) 2T+1 directly in
        // parity-coordinate space.
        ReedSolomonEncode(
            ops, buffer_bytes, m, m, m, m,
            data, work, NULL, exact_public_64_byte_shard);
        const unsigned tail_count = original_count - m;
        if (tail_count == 2 &&
            ops.ff8_multiply_add_2_sources_2_outputs)
        {
            for (unsigned output = 0; output < m; output += 2)
            {
                void* destinations[2] = {
                    work[output], work[output + 1]
                };
                const uint16_t logs0[2] = {
                    HighTailGeneratorLogs[0][m + output],
                    HighTailGeneratorLogs[0][m + output + 1]
                };
                const uint16_t logs1[2] = {
                    HighTailGeneratorLogs[1][m + output],
                    HighTailGeneratorLogs[1][m + output + 1]
                };
                ops.ff8_multiply_add_2_sources_2_outputs(
                    destinations, data[m], data[m + 1], logs0, logs1,
                    buffer_bytes);
            }
        }
        else
        {
            for (unsigned output = 0; output < m; ++output)
            {
                for (unsigned tail = 0; tail < tail_count; ++tail)
                {
                    MultiplyAddBytes(ops, work[output], data[m + tail],
                        HighTailGeneratorLogs[tail][m + output],
                        buffer_bytes);
                }
            }
        }
#if defined(LEO2_ENABLE_TEST_HOOKS)
        TestHighTailColumnCalls.fetch_add(1, std::memory_order_relaxed);
#endif
        return;
    }
    bool use_one_block_callback =
        (ops.kind == LEO2_BACKEND_AVX512 && buffer_bytes > 1024) ||
        (ops.kind == LEO2_BACKEND_AVX2 &&
         m == 8 && original_count == 8 && recovery_count == 8 &&
         buffer_bytes == 64);
    if (use_one_block_callback &&
        ops.ff8_high_encode_one_block &&
        (ops.ff8_high_encode_one_block_sides & m) != 0 &&
        whole_transform_size &&
        (original_count == m || original_count + 1U == m) &&
        (recovery_count == m || recovery_count + 1U == m) &&
        requested_output_count == recovery_count && dense_schedule)
    {
#if defined(LEO2_ENABLE_TEST_HOOKS)
        TestHighIFFTButterfly4OutCalls.fetch_add(
            m / 4, std::memory_order_relaxed);
        TestHighForwardFusedCalls.fetch_add(1, std::memory_order_relaxed);
        TestHighWholeTransformCalls.fetch_add(1, std::memory_order_relaxed);
#endif
        ops.ff8_high_encode_one_block(
            data, work, m |
                (original_count == m
                    ? 0U : backend::kFF8HighEncodeShortenedInput),
            FFTSkewStorage + m,
            FFTSkewStorage, buffer_bytes);
        return;
    }
    // work <- IFFT(data, m, m)

    const ffe_t* skewLUT = FFTSkewStorage + m;

    IFFT_DIT_Encoder(
        ops,
        buffer_bytes,
        data,
        original_count < m ? original_count : m,
        work,
        nullptr, // No xor output
        m,
        skewLUT,
        false,
        EncoderIFFTHigh,
        exact_public_64_byte_shard);

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
            data, // data source
            m,
            work + m, // temporary workspace
            work, // xor destination
            m,
            skewLUT,
            false,
            EncoderIFFTHigh,
            exact_public_64_byte_shard);
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
            data, // data source
            last_count,
            work + m, // temporary workspace
            work, // xor destination
            m,
            skewLUT,
            contiguous_temporary_work,
            EncoderIFFTHigh,
            exact_public_64_byte_shard);
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
        // The AVX2 leaf already uses the fused radix-four kernel through
        // 1 KiB.  Above that boundary, whole-forward measurements favor
        // retaining the same one-pass schedule for these explicit transform
        // sizes.  T=8/32/128 retain a final radix-two layer; T=16/64 are
        // complete radix-four stacks.  Keep this an allowlist so an unmeasured
        // parent size or another backend cannot silently inherit the policy.
        const bool fused_forward_size =
            m == 8 || m == 16 || m == 32 || m == 64 || m == 128;
        const bool prefer_fused_forward =
            (ops.kind == LEO2_BACKEND_AVX2 ||
             ops.kind == LEO2_BACKEND_GFNI) && fused_forward_size &&
            buffer_bytes > 1024;
#if defined(LEO2_ENABLE_TEST_HOOKS)
        if (prefer_fused_forward)
            TestHighForwardFusedCalls.fetch_add(1, std::memory_order_relaxed);
#endif
        FFT_DIT(
            ops,
            buffer_bytes,
            work,
            recovery_count,
            m,
            FFTSkewStorage,
            prefer_fused_forward);
    }
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
        recovery_count, m, data, work, NULL, buffer_bytes == 64U);
}


void ReedSolomonEncode(
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned m,
    const void* const* data,
    void** work)
{
    ReedSolomonEncode(
        backend::GetDefaultOps(), buffer_bytes, original_count, recovery_count,
        m, data, work);
}

void ReedSolomonEncodeOneBlockT8(
    const backend::Ops& ops,
    const void* const* data,
    void* const* work,
    uint64_t byte_count)
{
    LEO_DEBUG_ASSERT(ops.kind == LEO2_BACKEND_AVX2);
    LEO_DEBUG_ASSERT(ops.ff8_high_encode_one_block != NULL);
    LEO_DEBUG_ASSERT((ops.ff8_high_encode_one_block_sides & 8U) != 0);
#if defined(LEO2_ENABLE_TEST_HOOKS)
    TestHighIFFTButterfly4OutCalls.fetch_add(2, std::memory_order_relaxed);
    TestHighForwardFusedCalls.fetch_add(1, std::memory_order_relaxed);
    TestHighWholeTransformCalls.fetch_add(1, std::memory_order_relaxed);
#endif
    ops.ff8_high_encode_one_block(
        data, work, 8U, FFTSkewStorage + 8, FFTSkewStorage, byte_count);
}

void ReedSolomonEncodeK5R5T8(
    const backend::Ops& ops,
    const void* const* data,
    void* const* work,
    uint64_t byte_count)
{
    LEO_DEBUG_ASSERT(ops.kind == LEO2_BACKEND_AVX2);
    LEO_DEBUG_ASSERT(ops.ff8_high_encode_one_block != NULL);
    LEO_DEBUG_ASSERT((ops.ff8_high_encode_one_block_sides & 8U) != 0);
#if defined(LEO2_ENABLE_TEST_HOOKS)
    TestHighIFFTButterfly4OutCalls.fetch_add(2, std::memory_order_relaxed);
    TestHighForwardFusedCalls.fetch_add(1, std::memory_order_relaxed);
    TestHighWholeTransformCalls.fetch_add(1, std::memory_order_relaxed);
    TestHighK5R5PartialCalls.fetch_add(1, std::memory_order_relaxed);
#endif
    ops.ff8_high_encode_one_block(
        data, work, 8U | backend::kFF8HighEncodeK5R5Partial,
        FFTSkewStorage + 8, FFTSkewStorage, byte_count);
}

void ReedSolomonEncodeK6R6T8(
    const backend::Ops& ops,
    const void* const* data,
    void* const* work,
    uint64_t byte_count)
{
    LEO_DEBUG_ASSERT(ops.kind == LEO2_BACKEND_AVX2);
    LEO_DEBUG_ASSERT(ops.ff8_high_encode_one_block != NULL);
    LEO_DEBUG_ASSERT((ops.ff8_high_encode_one_block_sides & 8U) != 0);
#if defined(LEO2_ENABLE_TEST_HOOKS)
    TestHighIFFTButterfly4OutCalls.fetch_add(2, std::memory_order_relaxed);
    TestHighForwardFusedCalls.fetch_add(1, std::memory_order_relaxed);
    TestHighWholeTransformCalls.fetch_add(1, std::memory_order_relaxed);
#endif
    ops.ff8_high_encode_one_block(
        data, work, 8U | backend::kFF8HighEncodeK6R6Partial,
        FFTSkewStorage + 8, FFTSkewStorage, byte_count);
}

void ReedSolomonEncodeK7R7T8(
    const backend::Ops& ops,
    const void* const* data,
    void* const* work,
    uint64_t byte_count)
{
    LEO_DEBUG_ASSERT(ops.kind == LEO2_BACKEND_AVX2);
    LEO_DEBUG_ASSERT(ops.ff8_high_encode_one_block != NULL);
    LEO_DEBUG_ASSERT((ops.ff8_high_encode_one_block_sides & 8U) != 0);
#if defined(LEO2_ENABLE_TEST_HOOKS)
    TestHighIFFTButterfly4OutCalls.fetch_add(2, std::memory_order_relaxed);
    TestHighForwardFusedCalls.fetch_add(1, std::memory_order_relaxed);
    TestHighWholeTransformCalls.fetch_add(1, std::memory_order_relaxed);
#endif
#if defined(LEO2_HAVE_AVX2_T8_K8_B1024_DIRECT) && \
    !defined(LEO2_GFNI_VARIANT)
    if (byte_count == 1024)
    {
#if defined(LEO2_ENABLE_TEST_HOOKS)
        TestHighT8K7B1024DirectCalls.fetch_add(
            1, std::memory_order_relaxed);
#endif
        backend::AVX2FF8HighEncodeK7R7T8B1024(data, work);
        return;
    }
#endif
    ops.ff8_high_encode_one_block(
        data, work, 8U | backend::kFF8HighEncodeK7R7Partial,
        FFTSkewStorage + 8, FFTSkewStorage, byte_count);
}

void ReedSolomonEncodeK8R8T8B1024(
    const backend::Ops& ops,
    const void* const* data,
    void* const* work)
{
    LEO_DEBUG_ASSERT(ops.kind == LEO2_BACKEND_AVX2);
    (void)ops;
#if defined(LEO2_ENABLE_TEST_HOOKS)
    TestHighIFFTButterfly4OutCalls.fetch_add(2, std::memory_order_relaxed);
    TestHighForwardFusedCalls.fetch_add(1, std::memory_order_relaxed);
    TestHighWholeTransformCalls.fetch_add(1, std::memory_order_relaxed);
#endif
#if defined(LEO2_HAVE_AVX2_T8_K8_B1024_DIRECT) && \
    !defined(LEO2_GFNI_VARIANT)
#if defined(LEO2_ENABLE_TEST_HOOKS)
    TestHighT8K8B1024DirectCalls.fetch_add(
        1, std::memory_order_relaxed);
#endif
    backend::AVX2FF8HighEncodeK8R8T8B1024(data, work);
#else
    LEO_DEBUG_ASSERT(ops.ff8_high_encode_one_block != NULL);
    LEO_DEBUG_ASSERT((ops.ff8_high_encode_one_block_sides & 8U) != 0);
    ops.ff8_high_encode_one_block(
        data, work, 8U, FFTSkewStorage + 8, FFTSkewStorage, 1024);
#endif
}

void ReedSolomonEncodeT4Batch(
    const backend::Ops& ops,
    const void* const* data,
    void* const* recovery,
    uint32_t item_count,
    uint32_t original_count,
    uint32_t recovery_count,
    uint64_t byte_count)
{
    LEO_DEBUG_ASSERT(ops.kind == LEO2_BACKEND_AVX2);
    LEO_DEBUG_ASSERT(ops.ff8_high_encode_t4_batch != NULL);
#if defined(LEO2_ENABLE_TEST_HOOKS)
    TestHighSmallTransformCalls.fetch_add(
        item_count, std::memory_order_relaxed);
#endif
    ops.ff8_high_encode_t4_batch(
        data, recovery, item_count, original_count, recovery_count,
        FFTSkewStorage + 4, FFTSkewStorage, byte_count);
}

#if LEO2_EXPERIMENT_HIGH_T8_TWO_BLOCK_BINDING
void ReedSolomonEncodeK9T8(
    const backend::Ops& ops,
    const void* const* data,
    void* const* work,
    uint32_t recovery_count,
    uint64_t byte_count)
{
    LEO_DEBUG_ASSERT(ops.kind == LEO2_BACKEND_AVX2);
    LEO_DEBUG_ASSERT(ops.ff8_high_encode_one_block != NULL);
    LEO_DEBUG_ASSERT((ops.ff8_high_encode_one_block_sides & 8U) != 0);
    LEO_DEBUG_ASSERT(recovery_count >= 5 && recovery_count <= 8);
    LEO_DEBUG_ASSERT(byte_count == 256);
#if defined(LEO2_ENABLE_TEST_HOOKS)
    TestHighIFFTButterfly4OutCalls.fetch_add(2, std::memory_order_relaxed);
    TestHighForwardFusedCalls.fetch_add(1, std::memory_order_relaxed);
    TestHighWholeTransformCalls.fetch_add(1, std::memory_order_relaxed);
    TestHighTailColumnCalls.fetch_add(1, std::memory_order_relaxed);
    if (recovery_count == 5)
        TestHighK9R5TailCalls.fetch_add(1, std::memory_order_relaxed);
    else
        TestHighK9R6R8TailCalls.fetch_add(1, std::memory_order_relaxed);
#endif
    ops.ff8_high_encode_one_block(
        data, work,
        8U | backend::kFF8HighEncodeK9Tail |
            (recovery_count <<
                backend::kFF8HighEncodeK9OutputCountShift),
        FFTSkewStorage + 8, FFTSkewStorage, byte_count);
}

void ReedSolomonEncodeTwoBlocksT8(
    const backend::Ops& ops,
    const void* const* data,
    void* const* work,
    uint64_t byte_count)
{
    LEO_DEBUG_ASSERT(ops.kind == LEO2_BACKEND_AVX2);
    LEO_DEBUG_ASSERT(ops.ff8_high_encode_two_blocks_t8 != NULL);
#if defined(LEO2_ENABLE_TEST_HOOKS)
    /*
        Each fused callback replaces two eight-point inverse transforms,
        containing two radix-four groups apiece, followed by one forward
        transform.  Preserve the generic path's operation-count semantics so
        route tests can distinguish this callback from a silent fallback.
    */
    TestHighIFFTButterfly4OutCalls.fetch_add(4, std::memory_order_relaxed);
    TestHighForwardFusedCalls.fetch_add(1, std::memory_order_relaxed);
    TestHighWholeTransformCalls.fetch_add(1, std::memory_order_relaxed);
    TestHighTwoBlockCalls.fetch_add(1, std::memory_order_relaxed);
#endif
    ops.ff8_high_encode_two_blocks_t8(
        data, work, FFTSkewStorage + 8, FFTSkewStorage + 16,
        FFTSkewStorage, byte_count);
}

void ReedSolomonEncodeTwoBlocksT8Tiny(
    const backend::Ops& ops,
    const void* const* data,
    void* const* work,
    uint64_t byte_count)
{
    LEO_DEBUG_ASSERT(ops.kind == LEO2_BACKEND_AVX2);
    LEO_DEBUG_ASSERT(ops.ff8_high_encode_two_blocks_t8_tiny != NULL);
    LEO_DEBUG_ASSERT(byte_count >= 1 && byte_count < 64);
    if (!ops.ff8_high_encode_two_blocks_t8_tiny ||
        byte_count == 0 || byte_count >= 64)
        return;

#if defined(LEO2_ENABLE_TEST_HOOKS)
    TestHighIFFTButterfly4OutCalls.fetch_add(4, std::memory_order_relaxed);
    TestHighForwardFusedCalls.fetch_add(1, std::memory_order_relaxed);
    TestHighWholeTransformCalls.fetch_add(1, std::memory_order_relaxed);
#endif
    ops.ff8_high_encode_two_blocks_t8_tiny(
        data, work, FFTSkewStorage + 8, FFTSkewStorage + 16,
        FFTSkewStorage, byte_count);
}
#endif


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
        data,
        original_count,
        work,
        nullptr,
        p,
        FFTSkewStorage,
        false,
        EncoderIFFTLow,
        false);

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
        const bool dense_partial_candidate =
            ops.kind == LEO2_BACKEND_AVX2 && p == 16 &&
            requested_count == block_count && requested_count < p &&
            (requested_count == 4 || requested_count == 8 ||
             requested_count == 12);
        const unsigned partial_mode = dense_partial_candidate
            ? LowP16PartialDirectOutputMode()
            : 2U;
        const bool dense_partial_direct_output =
            dense_partial_candidate &&
            (partial_mode == 1U || partial_mode == 3U);
        bool direct_output = !use_sparse &&
            (requested_count == p || dense_partial_direct_output);
        for (unsigned i = 0; direct_output && i < requested_count; ++i)
            direct_output = recovery[recovery_offset + i] != NULL;
        if (dense_partial_candidate &&
            (partial_mode == 3U || partial_mode == 4U))
        {
            g_low_p16_partial_direct_output_route_selected =
                direct_output && dense_partial_direct_output;
        }
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
            TestLowDirectOutputShards.fetch_add(
                requested_count, std::memory_order_relaxed);
            TestLowAvoidedScatterBytes.fetch_add(
                static_cast<uint64_t>(requested_count) * buffer_bytes * 2U,
                std::memory_order_relaxed);
            if (requested_count < p)
                TestLowDirectPartialOutputBlocks.fetch_add(
                    1, std::memory_order_relaxed);
#endif
            continue;
        }
        for (unsigned i = 0; i < requested_count; ++i)
        {
            void* const output = recovery[recovery_offset + i];
            if (output)
            {
                memcpy(output, work[p + i], buffer_bytes);
#if defined(LEO2_ENABLE_TEST_HOOKS)
                TestLowScatterBytes.fetch_add(
                    buffer_bytes * 2U, std::memory_order_relaxed);
#endif
            }
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
    ReedSolomonEncodeLow(
        backend::GetDefaultOps(), buffer_bytes, original_count, recovery_count,
        p, data, recovery, work);
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
        for (unsigned r = 0; r < n_truncated; r += dist4)
        {
            if (!error_bits.IsNeeded(mip_level, r))
                continue;

            const ffe_t log_m01 = skewLUT[r + dist];
            const ffe_t log_m23 = skewLUT[r + dist * 3];
            const ffe_t log_m02 = skewLUT[r + dist * 2];

            // For each set of dist elements:
            for (unsigned i = r; i < r + dist; ++i)
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
        for (unsigned r = 0; r < n_truncated; r += 2)
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

    FWHT(error_locations, kOrder, n);

    for (unsigned i = 0; i < kOrder; ++i)
    {
        error_locations[i] = static_cast<ffe_t>(
            ((unsigned)error_locations[i] * (unsigned)LogWalsh[i]) % kModulus);
    }

    FWHT(error_locations, kOrder, kOrder);
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

    FWHT(locator_logs, n, n);

    const ffe_t* transformed_kernel =
        n == kOrder ? LogWalsh : ActiveLogWalsh + n - 2;
    for (unsigned i = 0; i < n; ++i)
    {
        locator_logs[i] = static_cast<ffe_t>(
            (static_cast<unsigned>(locator_logs[i]) *
                static_cast<unsigned>(transformed_kernel[i])) % kModulus);
    }

    FWHT(locator_logs, n, n);
}


void PrepareDecodeWalshActive(
    unsigned n,
    const uint8_t* erasures,
    ffe_t* locator_logs)
{
    PrepareDecodeWalshActiveCombined(n, erasures, nullptr, locator_logs);
}


void PrepareDecodeWalshActiveWithBackend(
    const backend::Ops& ops,
    unsigned n,
    const uint8_t* erasures,
    ffe_t* locator_logs)
{
    LEO_DEBUG_ASSERT(n >= 2 && n <= kOrder);

    if (n >= 32 && ops.ff8_walsh_locator)
    {
        const ffe_t* transformed_kernel =
            n == kOrder ? LogWalsh : ActiveLogWalsh + n - 2;
        ops.ff8_walsh_locator(
            erasures, transformed_kernel, locator_logs, n);
        return;
    }
    PrepareDecodeWalshActive(n, erasures, locator_logs);
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
    // Pinned multi-pattern whole-locator sweeps show that branch and transform
    // overhead move the crossover for small parents.  Ambiguous cells remain
    // direct: the first active cell had at least a 10% observed margin across
    // the frozen layouts on the calibration host.
    unsigned direct_cutoff;
    if (n <= 8)
        direct_cutoff = n;
    else if (n == 16)
        direct_cutoff = 8;
    else if (n == 32 || n == 128)
        direct_cutoff = 9;
    else if (n == 64)
        direct_cutoff = 8;
    else
        direct_cutoff = 7;
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


void PrepareDecodeKnownCount(
    unsigned n,
    const uint8_t* erasures,
    unsigned erasure_count,
    ffe_t* locator_logs)
{
    LEO_DEBUG_ASSERT(n >= 2 && n <= kOrder);
    LEO_DEBUG_ASSERT(erasure_count <= n);

    if (!IsDirectLocatorPreferred(n, erasure_count))
    {
        PrepareDecodeWalshActive(n, erasures, locator_logs);
        return;
    }

    PrepareDecodeDirect(n, erasures, locator_logs);
}


void PrepareDecode(
    unsigned n,
    const uint8_t* erasures,
    ffe_t* locator_logs)
{
    PrepareDecodeKnownCount(
        n, erasures, CountErasures(n, erasures, nullptr), locator_logs);
}


void PrepareDecodeWithPermanentKnownCount(
    unsigned n,
    const uint8_t* erasures,
    const uint8_t* permanent_erasures,
    const ffe_t* permanent_locator_logs,
    unsigned dynamic_erasure_count,
    ffe_t* locator_logs)
{
    LEO_DEBUG_ASSERT(n >= 2 && n <= kOrder);
    LEO_DEBUG_ASSERT(dynamic_erasure_count <= n);

    if (!IsDirectLocatorPreferred(n, dynamic_erasure_count))
    {
        PrepareDecodeWalshActiveCombined(
            n, erasures, permanent_erasures, locator_logs);
        return;
    }

    memcpy(locator_logs, permanent_locator_logs, n * sizeof(ffe_t));
    AddDirectLocatorContributions(
        n, erasures, permanent_erasures, locator_logs);
}


void PrepareDecodeWithPermanent(
    unsigned n,
    const uint8_t* erasures,
    const uint8_t* permanent_erasures,
    const ffe_t* permanent_locator_logs,
    ffe_t* locator_logs)
{
    PrepareDecodeWithPermanentKnownCount(
        n, erasures, permanent_erasures, permanent_locator_logs,
        CountErasures(n, erasures, permanent_erasures), locator_logs);
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
        VectorXOR(ops, buffer_bytes, width, work + i - width, work + i);
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

    for (unsigned i = 0; i < n; ++i)
    {
        if (coordinate_data[i])
            mul_mem(
                ops, work[i], coordinate_data[i], locator_logs[i], buffer_bytes);
        else
            memset(work[i], 0, buffer_bytes);
    }

    IFFT_DIT_Decoder(ops, buffer_bytes, input_count, work, n, FFTSkewStorage);

    AddFormalDerivative(ops, buffer_bytes, n, work);

#ifdef LEO_ERROR_BITFIELD_OPT
    FFT_DIT_ErrorBits(
        ops, buffer_bytes, work, n, n, FFTSkewStorage, error_bits);
#else
    FFT_DIT(ops, buffer_bytes, work, n, n, FFTSkewStorage);
#endif // LEO_ERROR_BITFIELD_OPT

    for (unsigned i = 0; i < n; ++i)
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
    ReedSolomonDecodePrepared(
        backend::GetDefaultOps(), buffer_bytes, n, coordinate_data,
        requested_outputs, locator_logs, work);
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

    for (unsigned i = 0; i < n; ++i)
    {
        if (coordinate_data[i])
            mul_mem(
                ops, work[i], coordinate_data[i], locator_logs[i], buffer_bytes);
        else
            memset(work[i], 0, buffer_bytes);
    }

    if (input_count != 0)
        IFFT_DIT_Decoder(
            ops, buffer_bytes, input_count, work, n, FFTSkewStorage);

    AddFormalDerivative(ops, buffer_bytes, n, work);

#ifdef LEO_ERROR_BITFIELD_OPT
    FFT_DIT_ErrorBits(
        ops, buffer_bytes, work, n, n, FFTSkewStorage, output_dependencies);
#else
    (void)output_dependencies;
    FFT_DIT(ops, buffer_bytes, work, n, n, FFTSkewStorage);
#endif

    if (reveal_outputs)
        for (unsigned i = 0; i < requested_count; ++i)
        {
            const uint32_t coordinate = requested_coordinates[i];
            LEO_DEBUG_ASSERT(coordinate < n);
            mul_mem_inplace(
                ops, work[coordinate], kModulus - locator_logs[coordinate],
                buffer_bytes);
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
    ReedSolomonDecodePlanned(
        backend::GetDefaultOps(), buffer_bytes, n, coordinate_data, input_count,
        requested_coordinates, requested_count, output_dependencies,
        locator_logs, true, work);
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

    // c_n = product(V_n \ {0}).  R10's full-field form silently has c_m=1
    // because s_m(x) = x^(2^m) - x gives s_m' = 1 (eq. (63)).  This numerator
    // is retained for representation independence, NOT because production
    // needs it: for the declared Cantor bases c_j = 1 and p_i = 1 for every j
    // in both GF8 and GF16, so the factor reduces to 1 here.  The general
    // active-parent form is what the exhaustive GF(2^4) oracle proves, using
    // a plain polynomial basis where c_n != 1
    // (tests/leopard2/direct_oracle.cpp, test_direct_oracle.cpp).
    const ffe_t active_derivative_log = SubspaceDerivativeLog(n);

    // p_(N-T) from the normalized LCH basis.
    const ffe_t normalization_log = LchNormalizerLog(n - t);

    const unsigned block_count = n / t;
    for (unsigned block = 1; block < block_count; ++block)
    {
        const unsigned coordinate = block * t;
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


bool SetLowP16PartialDirectOutputEnabledForDiagnostics(bool enabled)
{
    g_low_p16_partial_direct_output_route_selected = false;
    g_low_p16_partial_direct_output_mode.store(
        enabled ? 3U : 4U, std::memory_order_release);
    return true;
}


unsigned LowP16PartialDirectOutputModeForDiagnostics()
{
    const unsigned mode = g_low_p16_partial_direct_output_mode.load(
        std::memory_order_acquire);
    return mode == 3U ? 1U : mode == 4U ? 2U : mode;
}


bool LowP16PartialDirectOutputRouteSelectedForDiagnostics()
{
    return g_low_p16_partial_direct_output_route_selected;
}


bool FinishLowP16PartialDirectOutputRouteProbeForDiagnostics()
{
    const unsigned mode = g_low_p16_partial_direct_output_mode.load(
        std::memory_order_acquire);
    if (mode != 3U && mode != 4U)
        return false;
    g_low_p16_partial_direct_output_mode.store(
        mode == 3U ? 1U : 2U, std::memory_order_release);
    return true;
}


#if defined(LEO2_ENABLE_TEST_HOOKS)

void TestOnlyResetTransformCallsiteCounts()
{
    TestIFFTDIT4Calls.store(0, std::memory_order_relaxed);
    TestIFFTDIT4XorCalls.store(0, std::memory_order_relaxed);
    TestFFTDIT4Calls.store(0, std::memory_order_relaxed);
}


TestOnlyTransformCallsiteCounts TestOnlyGetTransformCallsiteCounts()
{
    TestOnlyTransformCallsiteCounts result;
    result.ifft_dit4 = TestIFFTDIT4Calls.load(std::memory_order_relaxed);
    result.ifft_dit4_xor =
        TestIFFTDIT4XorCalls.load(std::memory_order_relaxed);
    result.fft_dit4 = TestFFTDIT4Calls.load(std::memory_order_relaxed);
    return result;
}


void TestOnlyResetLowEncodeCounts()
{
    TestLowFFTButterfly2OutCalls.store(0, std::memory_order_relaxed);
    TestLowFFTButterfly4OutCalls.store(0, std::memory_order_relaxed);
    TestLowFFTButterfly8OutCalls.store(0, std::memory_order_relaxed);
    TestLowDirectOutputBlocks.store(0, std::memory_order_relaxed);
    TestLowDirectPartialOutputBlocks.store(0, std::memory_order_relaxed);
    TestLowDirectOutputShards.store(0, std::memory_order_relaxed);
    TestLowAvoidedScatterBytes.store(0, std::memory_order_relaxed);
    TestLowScatterBytes.store(0, std::memory_order_relaxed);
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
    result.fft_butterfly8_out_of_place =
        TestLowFFTButterfly8OutCalls.load(std::memory_order_relaxed);
    result.direct_output_blocks =
        TestLowDirectOutputBlocks.load(std::memory_order_relaxed);
    result.direct_partial_output_blocks =
        TestLowDirectPartialOutputBlocks.load(std::memory_order_relaxed);
    result.direct_output_shards =
        TestLowDirectOutputShards.load(std::memory_order_relaxed);
    result.avoided_scatter_bytes =
        TestLowAvoidedScatterBytes.load(std::memory_order_relaxed);
    result.scatter_bytes =
        TestLowScatterBytes.load(std::memory_order_relaxed);
    result.direct_output_butterfly2_out_of_place =
        TestLowDirectFFTButterfly2OutCalls.load(std::memory_order_relaxed);
    result.direct_output_butterfly4_out_of_place =
        TestLowDirectFFTButterfly4OutCalls.load(std::memory_order_relaxed);
    return result;
}


void TestOnlySetLowP16PartialDirectOutputEnabled(bool enabled)
{
    g_low_p16_partial_direct_output_mode.store(
        enabled ? 1U : 2U, std::memory_order_relaxed);
}


bool TestOnlyLowP16PartialDirectOutputEnabled()
{
    const unsigned mode = LowP16PartialDirectOutputMode();
    return mode == 1U || mode == 3U;
}


void TestOnlyResetHighEncodeCounts()
{
    TestHighIFFTButterfly4OutCalls.store(0, std::memory_order_relaxed);
    TestHighInputCopyShards.store(0, std::memory_order_relaxed);
    TestHighForwardFusedCalls.store(0, std::memory_order_relaxed);
    TestHighWholeTransformCalls.store(0, std::memory_order_relaxed);
    TestHighK5R5PartialCalls.store(0, std::memory_order_relaxed);
    TestHighTwoBlockCalls.store(0, std::memory_order_relaxed);
    TestHighSmallTransformCalls.store(0, std::memory_order_relaxed);
    TestHighT2PackedCalls.store(0, std::memory_order_relaxed);
    TestHighT4PackedCalls.store(0, std::memory_order_relaxed);
    TestHighT8PackedCalls.store(0, std::memory_order_relaxed);
    TestHighT8K7B1024DirectCalls.store(0, std::memory_order_relaxed);
    TestHighT8K8B1024DirectCalls.store(0, std::memory_order_relaxed);
    TestHighBalancedB64PackedCalls.store(0, std::memory_order_relaxed);
    TestHighFinalIFFT2RangeCalls.store(0, std::memory_order_relaxed);
    TestHighTailColumnCalls.store(0, std::memory_order_relaxed);
    TestHighHalfTailColumnCalls.store(0, std::memory_order_relaxed);
    TestHighK9R5TailCalls.store(0, std::memory_order_relaxed);
    TestHighK9R6R8TailCalls.store(0, std::memory_order_relaxed);
}


TestOnlyHighEncodeCounts TestOnlyGetHighEncodeCounts()
{
    TestOnlyHighEncodeCounts result;
    result.ifft_butterfly4_out_of_place =
        TestHighIFFTButterfly4OutCalls.load(std::memory_order_relaxed);
    result.input_copy_shards =
        TestHighInputCopyShards.load(std::memory_order_relaxed);
    result.forward_fused_calls =
        TestHighForwardFusedCalls.load(std::memory_order_relaxed);
    result.whole_transform_calls =
        TestHighWholeTransformCalls.load(std::memory_order_relaxed);
    result.k5r5_partial_calls =
        TestHighK5R5PartialCalls.load(std::memory_order_relaxed);
    result.two_block_calls =
        TestHighTwoBlockCalls.load(std::memory_order_relaxed);
    result.small_transform_calls =
        TestHighSmallTransformCalls.load(std::memory_order_relaxed);
    result.t2_packed_calls =
        TestHighT2PackedCalls.load(std::memory_order_relaxed);
    result.t4_packed_calls =
        TestHighT4PackedCalls.load(std::memory_order_relaxed);
    result.t8_packed_calls =
        TestHighT8PackedCalls.load(std::memory_order_relaxed);
    result.t8_k7_b1024_direct_calls =
        TestHighT8K7B1024DirectCalls.load(std::memory_order_relaxed);
    result.t8_k8_b1024_direct_calls =
        TestHighT8K8B1024DirectCalls.load(std::memory_order_relaxed);
    result.balanced_b64_packed_calls =
        TestHighBalancedB64PackedCalls.load(std::memory_order_relaxed);
    result.final_ifft2_range_calls =
        TestHighFinalIFFT2RangeCalls.load(std::memory_order_relaxed);
    result.tail_column_calls =
        TestHighTailColumnCalls.load(std::memory_order_relaxed);
    result.half_tail_column_calls =
        TestHighHalfTailColumnCalls.load(std::memory_order_relaxed);
    result.k9r5_tail_calls =
        TestHighK9R5TailCalls.load(std::memory_order_relaxed);
    result.k9r6r8_tail_calls =
        TestHighK9R6R8TailCalls.load(std::memory_order_relaxed);
    return result;
}

void TestOnlyRecordT2PackedCall()
{
    TestHighT2PackedCalls.fetch_add(1, std::memory_order_relaxed);
}

void TestOnlyRecordT4PackedCall()
{
    TestHighT4PackedCalls.fetch_add(1, std::memory_order_relaxed);
}

void TestOnlyRecordT8PackedCall()
{
    TestHighT8PackedCalls.fetch_add(1, std::memory_order_relaxed);
}

void TestOnlyRecordBalancedB64PackedCall()
{
    TestHighBalancedB64PackedCalls.fetch_add(1, std::memory_order_relaxed);
}


void TestOnlySetHighFinalIFFT2RangeDisabled(bool disabled)
{
    TestDisableHighFinalIFFT2Range.store(disabled, std::memory_order_relaxed);
}


bool TestOnlyHighFinalIFFT2RangeDisabled()
{
    return TestDisableHighFinalIFFT2Range.load(std::memory_order_relaxed);
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
    TestHighReceiveIFFTButterfly4OutCalls.store(
        0, std::memory_order_relaxed);
    TestHighReceiveCopyShards.store(0, std::memory_order_relaxed);
    TestHighReceiveZeroShards.store(0, std::memory_order_relaxed);
    TestHighLocatorWeightedIFFTButterfly4Calls.store(
        0, std::memory_order_relaxed);
    TestHighLocatorScaleRowsElided.store(0, std::memory_order_relaxed);
    TestHighLocatorInactiveRows.store(0, std::memory_order_relaxed);
    TestHighLocatorFallbackRows.store(0, std::memory_order_relaxed);
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
    result.receive_ifft_butterfly4_out_of_place =
        TestHighReceiveIFFTButterfly4OutCalls.load(
            std::memory_order_relaxed);
    result.receive_copy_shards =
        TestHighReceiveCopyShards.load(std::memory_order_relaxed);
    result.receive_zero_shards =
        TestHighReceiveZeroShards.load(std::memory_order_relaxed);
    result.locator_weighted_ifft_butterfly4 =
        TestHighLocatorWeightedIFFTButterfly4Calls.load(
            std::memory_order_relaxed);
    result.locator_scale_rows_elided =
        TestHighLocatorScaleRowsElided.load(std::memory_order_relaxed);
    result.locator_inactive_rows =
        TestHighLocatorInactiveRows.load(std::memory_order_relaxed);
    result.locator_fallback_rows =
        TestHighLocatorFallbackRows.load(std::memory_order_relaxed);
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
    const backend::Ops& ops,
    uint64_t buffer_bytes,
    unsigned size,
    void** work)
{
    AddFormalDerivative(ops, buffer_bytes, size, work);
}


void TestOnlyAddFormalDerivative(
    uint64_t buffer_bytes,
    unsigned size,
    void** work)
{
    TestOnlyAddFormalDerivative(
        backend::GetDefaultOps(), buffer_bytes, size, work);
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

    for (unsigned i = 0; i < n; ++i)
    {
        if (coordinate_data[i])
            mul_mem(
                ops, work[i], coordinate_data[i], locator_logs[i], buffer_bytes);
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
            ops, buffer_bytes, input_count, work + offset, p,
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
        for (unsigned i = 0; i < p; ++i)
            muladd_mem(ops, work[i], work[offset + i],
                block_factors[block - 1], buffer_bytes);
    }

#ifdef LEO_ERROR_BITFIELD_OPT
    FFT_DIT_ErrorBits(
        ops, buffer_bytes, work, p, p, FFTSkewStorage, error_bits);
#else
    FFT_DIT(ops, buffer_bytes, work, p, p, FFTSkewStorage);
#endif

    for (unsigned i = 0; i < p; ++i)
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
    ReedSolomonDecodeLowPrepared(
        backend::GetDefaultOps(), buffer_bytes, n, p, coordinate_data,
        requested_outputs, locator_logs, block_factors, work);
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

    const unsigned block_count = n / p;
    unsigned seed_block = 0;
    if (block_input_counts[0] == 0)
    {
        seed_block = 1;
        while (seed_block < block_count &&
               block_input_counts[seed_block] == 0)
            ++seed_block;
    }
    const bool seed_from_later = seed_block > 0 && seed_block < block_count;
    unsigned input_plan_index = 0;
    if (!seed_from_later)
    {
        // Preserve the ordinary block-zero-live path without seed-specific
        // branches in its byte and transform loops.
        for (unsigned i = 0; i < n; ++i)
        {
            if (coordinate_data[i])
                mul_mem(ops, work[i], coordinate_data[i], locator_logs[i],
                    buffer_bytes);
            else
                memset(work[i], 0, buffer_bytes);
        }

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

        AddFormalDerivative(ops, buffer_bytes, p, work);

        // The transformed nonzero block is dead after its weighted
        // contribution.
        for (unsigned block = 1; block < block_count; ++block)
        {
            if (block_input_counts[block] == 0)
                continue;
            const unsigned offset = block * p;
            for (unsigned i = 0; i < p; ++i)
                muladd_mem(ops, work[i], work[offset + i],
                    block_factors[block - 1], buffer_bytes);
        }
    }
    else
    {
        const unsigned seed_offset = seed_block * p;

        // Block zero is empty.  Seed its accumulator directly from the first
        // live later block.  Folding that block's scalar into the input
        // multiplication is exact because scalar multiplication commutes with
        // the inverse LCH transform.  Its ordinary materialized slot is dead.
        for (unsigned i = 0; i < p; ++i)
        {
            const unsigned coordinate = seed_offset + i;
            if (coordinate_data[coordinate])
                mul_mem(ops, work[i], coordinate_data[coordinate],
                    AddMod(locator_logs[coordinate],
                        block_factors[seed_block - 1]),
                    buffer_bytes);
            else
                memset(work[i], 0, buffer_bytes);
        }
        for (unsigned block = 1; block < block_count; ++block)
        {
            if (block == seed_block || block_input_counts[block] == 0)
                continue;
            const unsigned offset = block * p;
            for (unsigned i = 0; i < p; ++i)
            {
                const unsigned coordinate = offset + i;
                if (coordinate_data[coordinate])
                    mul_mem(ops, work[coordinate],
                        coordinate_data[coordinate], locator_logs[coordinate],
                        buffer_bytes);
                else
                    memset(work[coordinate], 0, buffer_bytes);
            }
        }

        for (unsigned block = 1; block < block_count; ++block)
        {
            const unsigned input_count = block_input_counts[block];
            LEO_DEBUG_ASSERT(input_count <= p);
            if (input_count == 0)
                continue;
            const unsigned offset = block * p;
            void** const block_work =
                block == seed_block ? work : work + offset;
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
                        ops, buffer_bytes, *pruned, block_work);
                LEO_DEBUG_ASSERT(executed);
                (void)executed;
            }
            else
                IFFT_DIT_Decoder(ops, buffer_bytes, input_count, block_work, p,
                    FFTSkewStorage + offset);
        }

        // Algorithm 4 applies the derivative only to block zero.  That term is
        // zero here, so the weighted later-block seed receives no derivative.
        for (unsigned block = 1; block < block_count; ++block)
        {
            if (block == seed_block || block_input_counts[block] == 0)
                continue;
            const unsigned offset = block * p;
            for (unsigned i = 0; i < p; ++i)
                muladd_mem(ops, work[i], work[offset + i],
                    block_factors[block - 1], buffer_bytes);
        }
    }
    LEO_DEBUG_ASSERT(input_plan_index == input_plan_count);

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
        FFT_DIT_ErrorBits(
            ops, buffer_bytes, work, p, p, FFTSkewStorage,
            output_dependencies);
#else
        (void)output_dependencies;
        FFT_DIT(ops, buffer_bytes, work, p, p, FFTSkewStorage);
#endif
    }

    if (reveal_outputs)
        for (unsigned i = 0; i < requested_count; ++i)
        {
            const uint32_t coordinate = requested_coordinates[i];
            LEO_DEBUG_ASSERT(coordinate < p);
            mul_mem_inplace(
                ops, work[coordinate], kModulus - locator_logs[coordinate],
                buffer_bytes);
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

#if LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL
bool ReedSolomonDecodeLowP32B64TerminalExperimental(
    const backend::Ops& ops,
    const void* const* coordinate_data,
    const uint32_t* requested_coordinates,
    unsigned requested_count,
    const ffe_t* locator_logs,
    ffe_t block_factor,
    void* const* restored,
    void** work)
{
    // Keep coordinate and normalization ownership in the field layer.  The
    // backend receives the same shift-zero and shift-P tables used by the
    // mature Algorithm 4 implementation, so the experiment changes only the
    // execution schedule and never the code definition.
    return backend::AVX2FF8LowP32B64Terminal(
        ops, coordinate_data, requested_coordinates, requested_count,
        locator_logs, block_factor,
        FFTSkewStorage, FFTSkewStorage + 32, FFTSkewStorage,
        restored, work);
}

bool ReedSolomonDecodeLowP128B64TerminalExperimental(
    const backend::Ops& ops,
    const void* const* coordinate_data,
    const uint32_t* requested_coordinates,
    unsigned requested_count,
    const ffe_t* locator_logs,
    ffe_t block_factor,
    void* const* restored,
    void** work)
{
    return backend::AVX2FF8LowP128B64Terminal(
        ops, coordinate_data, requested_coordinates, requested_count,
        locator_logs, block_factor,
        FFTSkewStorage, FFTSkewStorage + 128, FFTSkewStorage,
        restored, work);
}
#endif


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
    ReedSolomonDecodeLowPlanned(
        backend::GetDefaultOps(), buffer_bytes, n, p, coordinate_data,
        block_input_counts, requested_coordinates, requested_count,
        output_dependencies, locator_logs, block_factors, work);
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
    unsigned seed_block = 0;
    if (first_input_count == 0)
    {
        seed_block = 1;
        while (seed_block < block_count &&
               block_input_counts[seed_block] == 0)
            ++seed_block;
    }
    const bool seed_from_later = seed_block > 0 && seed_block < block_count;

    unsigned input_plan_index = 0;
    if (!seed_from_later)
    {
        for (unsigned i = 0; i < p; ++i)
        {
            if (coordinate_data[i])
                mul_mem(
                    ops, accumulator[i], coordinate_data[i], locator_logs[i],
                    buffer_bytes);
            else
                memset(accumulator[i], 0, buffer_bytes);
        }
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

        // Preserve the ordinary block-zero-live reduction without
        // seed-specific branches in its byte loop.
        for (unsigned block = 1; block < block_count; ++block)
        {
            const unsigned input_count = block_input_counts[block];
            LEO_DEBUG_ASSERT(input_count <= p);
            if (input_count == 0)
                continue;
            const unsigned offset = block * p;
            for (unsigned i = 0; i < p; ++i)
            {
                const unsigned coordinate = offset + i;
                if (coordinate_data[coordinate])
                    mul_mem(ops, tile[i], coordinate_data[coordinate],
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
                IFFT_DIT_Decoder(ops, buffer_bytes, input_count, tile, p,
                    FFTSkewStorage + offset);
            for (unsigned i = 0; i < p; ++i)
                muladd_mem(ops, accumulator[i], tile[i],
                    block_factors[block - 1], buffer_bytes);
        }
    }
    else
    {
        // Block zero's derivative term is zero.  Transform the first later
        // live block directly in the accumulator with its scalar folded into
        // the input multiplication; use the ordinary tile for later blocks.
        for (unsigned block = 1; block < block_count; ++block)
        {
            const unsigned input_count = block_input_counts[block];
            LEO_DEBUG_ASSERT(input_count <= p);
            if (input_count == 0)
                continue;
            const unsigned offset = block * p;
            const bool seed_accumulator = block == seed_block;
            void** const block_work = seed_accumulator ? accumulator : tile;
            for (unsigned i = 0; i < p; ++i)
            {
                const unsigned coordinate = offset + i;
                if (coordinate_data[coordinate])
                {
                    ffe_t multiplier_log = locator_logs[coordinate];
                    if (seed_accumulator)
                        multiplier_log = AddMod(
                            multiplier_log, block_factors[block - 1]);
                    mul_mem(ops, block_work[i], coordinate_data[coordinate],
                        multiplier_log, buffer_bytes);
                }
                else
                    memset(block_work[i], 0, buffer_bytes);
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
                        ops, buffer_bytes, *pruned, block_work);
                LEO_DEBUG_ASSERT(executed);
                (void)executed;
            }
            else
                IFFT_DIT_Decoder(ops, buffer_bytes, input_count, block_work, p,
                    FFTSkewStorage + offset);
            if (!seed_accumulator)
                for (unsigned i = 0; i < p; ++i)
                    muladd_mem(ops, accumulator[i], tile[i],
                        block_factors[block - 1], buffer_bytes);
        }
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
        for (unsigned i = 0; i < requested_count; ++i)
        {
            const uint32_t coordinate = requested_coordinates[i];
            LEO_DEBUG_ASSERT(coordinate < p);
            mul_mem_inplace(
                ops, accumulator[coordinate],
                kModulus - locator_logs[coordinate], buffer_bytes);
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


// Algorithm 5 constructs h by adding one inverse transform per T-point
// receive block and then evaluating that sum on the zero-shift block.  The
// block-zero contribution is an exact inverse/forward pair:
//
//   FFT_0(IFFT_0(F_0) + H_later) = F_0 + FFT_0(H_later).
//
// Keep F_0 in its received evaluation form.  When later blocks contributed,
// work contains H_later: transform only that sum and XOR the immutable F_0
// sources into it.  With no later contribution, staging F_0 directly removes
// both transforms.  Null sources are the selected-codeword zero values.
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
        IFFT_DIT_DecoderFromSources(
            ops, buffer_bytes, input_count, coordinate_data + offset,
            block_work, have_later_contribution ? work : NULL, t,
            FFTSkewStorage + offset);
#if defined(LEO2_ENABLE_TEST_HOOKS)
        (have_later_contribution ? TestHighSyndromeAccumulatedBlocks
                                 : TestHighSyndromeMaterializedBlocks)
            .fetch_add(1, std::memory_order_relaxed);
#endif
        have_later_contribution = true;
    }

    // h on V_t, then z = h * Lambda on V_t.
    FinishHighDecodeSyndrome(
        ops, buffer_bytes, coordinate_data, block_zero_input_count,
        have_later_contribution, t, work);
    IFFT_DIT_DecoderWeightedLocator(
        ops, buffer_bytes, coordinate_data, locator_logs,
        work, t, FFTSkewStorage);

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
        for (unsigned i = 0; i < requested_count; ++i)
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
    ReedSolomonDecodeHighPrepared(
        backend::GetDefaultOps(), buffer_bytes, n, t, coordinate_data,
        requested_outputs, locator_logs, output_factors, work);
}


static bool ShouldUsePrunedHighSyndromeSink(
    const backend::Ops& ops,
    uint64_t buffer_bytes)
{
    // The immutable-schedule preflight is fixed per execution.  Pinned ABBA
    // crossover sweeps found AVX2 neutral-to-faster from 1 KiB, while GF8
    // SSSE3 did not become a credible win until 64 KiB.  Keep smaller shards
    // on the materialize-then-XOR oracle instead of paying that fixed cost.
    if (ops.kind == LEO2_BACKEND_AVX2 ||
        ops.kind == LEO2_BACKEND_AVX512 ||
        ops.kind == LEO2_BACKEND_GFNI)
        return buffer_bytes >= 1024;
    if (ops.kind == LEO2_BACKEND_SSSE3)
        return buffer_bytes >= 65536;
    return false;
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
        if (pruned)
        {
            StageHighDecodeSources(
                buffer_bytes, coordinate_data + offset, block_work, t);
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
#if defined(LEO2_ENABLE_TEST_HOOKS)
            TestHighSyndromeMaterializedBlocks.fetch_add(
                1, std::memory_order_relaxed);
#endif
            if (have_later_contribution)
                VectorXOR(ops, buffer_bytes, t, work, block_work);
        }
        else
        {
            // Source fusion is intentionally confined to unpruned blocks;
            // exact plans above retain their independently qualified sink.
            IFFT_DIT_DecoderFromSources(
                ops, buffer_bytes, input_count, coordinate_data + offset,
                block_work, have_later_contribution ? work : NULL, t,
                FFTSkewStorage + offset);
#if defined(LEO2_ENABLE_TEST_HOOKS)
            (have_later_contribution ? TestHighSyndromeAccumulatedBlocks
                                     : TestHighSyndromeMaterializedBlocks)
                .fetch_add(1, std::memory_order_relaxed);
#endif
        }
        have_later_contribution = true;
    }
    LEO_DEBUG_ASSERT(input_plan_index == input_plan_count);

    FinishHighDecodeSyndrome(
        ops, buffer_bytes, coordinate_data, block_input_counts[0],
        have_later_contribution, t, work);
    IFFT_DIT_DecoderWeightedLocator(
        ops, buffer_bytes, coordinate_data, locator_logs,
        work, t, FFTSkewStorage);

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
                for (unsigned i = 0; i < t; ++i)
                    memcpy(work[offset + i], work[i], buffer_bytes);
                const bool fallback =
                    leopard2_internal::ExecutePrunedTransformPlan(
                        ops, buffer_bytes, *pruned, work + offset);
                LEO_DEBUG_ASSERT(fallback);
                (void)fallback;
            }
        }
        else
        {
            FFT_DIT_FromCoefficients(
                ops, buffer_bytes, work, work + offset,
                descriptor.requested_prefix, t,
                FFTSkewStorage + offset, NULL, SourceEvaluationHighDecode);
        }
        for (uint32_t i = descriptor.requested_begin;
             i < descriptor.requested_end;
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
    ReedSolomonDecodeHighPlanned(
        backend::GetDefaultOps(), buffer_bytes, n, t, coordinate_data,
        block_input_counts, requested_coordinates, output_blocks,
        output_block_count, locator_logs, output_factors, work);
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
        if (pruned)
        {
            StageHighDecodeSources(
                buffer_bytes, coordinate_data + offset, block_work, t);
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
            if (have_later_contribution)
                VectorXOR(
                    ops, buffer_bytes, t, accumulator, block_work);
        }
        else
        {
            // The first later block becomes H_later directly.  Subsequent
            // dead tiles accumulate their final inverse layer into it.
            IFFT_DIT_DecoderFromSources(
                ops, buffer_bytes, input_count, coordinate_data + offset,
                block_work,
                have_later_contribution ? accumulator : NULL,
                t, FFTSkewStorage + offset);
#if defined(LEO2_ENABLE_TEST_HOOKS)
            (have_later_contribution ? TestHighSyndromeAccumulatedBlocks
                                     : TestHighSyndromeMaterializedBlocks)
                .fetch_add(1, std::memory_order_relaxed);
#endif
        }
        have_later_contribution = true;
    }
    LEO_DEBUG_ASSERT(input_plan_index == input_plan_count);

    FinishHighDecodeSyndrome(
        ops, buffer_bytes, coordinate_data, block_input_counts[0],
        have_later_contribution, t, accumulator);
    IFFT_DIT_DecoderWeightedLocator(
        ops, buffer_bytes, coordinate_data, locator_logs,
        accumulator, t, FFTSkewStorage);

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
                for (unsigned i = 0; i < t; ++i)
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
        for (uint32_t i = descriptor.requested_begin;
             i < descriptor.requested_end;
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

    FWHT(error_locations, kOrder, m + original_count);

    for (unsigned i = 0; i < kOrder; ++i)
        error_locations[i] = ((unsigned)error_locations[i] * (unsigned)LogWalsh[i]) % kModulus;

    FWHT(error_locations, kOrder, kOrder);

    // work <- recovery data

    for (unsigned i = 0; i < recovery_count; ++i)
    {
        if (recovery[i])
            mul_mem(
                ops, work[i], recovery[i], error_locations[i], buffer_bytes);
        else
            memset(work[i], 0, buffer_bytes);
    }
    for (unsigned i = recovery_count; i < m; ++i)
        memset(work[i], 0, buffer_bytes);

    // work <- original data

    for (unsigned i = 0; i < original_count; ++i)
    {
        if (original[i])
            mul_mem(
                ops, work[m + i], original[i], error_locations[m + i],
                buffer_bytes);
        else
            memset(work[m + i], 0, buffer_bytes);
    }
    for (unsigned i = m + original_count; i < n; ++i)
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

        VectorXOR(
            ops,
            buffer_bytes,
            width,
            work + i - width,
            work + i);
    }

    // work <- FFT(work, n, 0) truncated to m + original_count

    const unsigned output_count = m + original_count;

#ifdef LEO_ERROR_BITFIELD_OPT
    FFT_DIT_ErrorBits(
        ops, buffer_bytes, work, output_count, n, FFTSkewStorage, error_bits);
#else
    FFT_DIT(ops, buffer_bytes, work, output_count, n, FFTSkewStorage);
#endif

    // Reveal erasures

    for (unsigned i = 0; i < original_count; ++i)
        if (!original[i])
            mul_mem(
                ops, work[i], work[i + m],
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
    ReedSolomonDecode(
        backend::GetDefaultOps(), buffer_bytes, original_count, recovery_count,
        m, n, original, recovery, work);
}


//------------------------------------------------------------------------------
// API

static bool IsInitialized = false;

bool Initialize()
{
    if (IsInitialized)
        return true;

    InitializeLogarithmTables();
    if (!InitializeHighTailGeneratorLogs())
        return false;
    if (!InitializeMultiplyTables())
        return false;
    FFTInitialize();

    IsInitialized = true;
    return true;
}


}} // namespace leopard::ff8

#endif // LEO_HAS_FF8
